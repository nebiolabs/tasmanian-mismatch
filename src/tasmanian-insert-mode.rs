use clap::Parser;
use rayon::prelude::*;
use rust_htslib::bam::{record::Aux, FetchDefinition, IndexedReader, Read, Reader, Record};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use rustmanian_mismatch::{
	base_to_char, compute_read_len_max_from_sample_bam, load_reference_genome, ReferenceGenome, complement
};

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct InsertKey {
	base_change: String,
	read_num: u8,
	base_position: usize,
	/// 1 = first read in reference coordinates, 2 = second
	reference_order: u8,
}

#[derive(Debug, Clone)]
struct SoftclipComparison {
	read_pos: usize,
	ref_pos: usize,
	read_base: char,
	ref_base: char,
}

#[derive(Debug, Clone, clap::ValueEnum, PartialEq)]
enum OverlapMode {
	/// Split the overlap at its midpoint; read 1 keeps [ov_start, mid), read 2 keeps [mid, ov_end)
	Cut,
	/// Scale each read's positions to fill its half of [1, 2*max_read_len] via a cubic spline; no bases are dropped
	Stretch,
}

#[derive(Debug, Clone, Copy, clap::ValueEnum, PartialEq, Eq)]
enum PositionMode {
	Read,
	Insert,
}

#[derive(Parser, Debug)]
#[command(name = "tasmanian-insert-mode")]
#[command(version, about = "Parallel indexed BAM mismatch caller by read and position")]
struct Args {
	/// Coordinate-sorted BAM file path (requires .bai)
	bam_path: String,

	/// Reference FASTA path
	reference_path: String,

	/// Number of threads to use (0 keeps rayon default)
	#[arg(short = 't', long, default_value_t = 0)]
	threads: usize,

	/// Region size in bp for indexed chunking
	#[arg(short = 'r', long, default_value_t = 1_000_000)]
	region_size: usize,

	/// Minimum base quality
	#[arg(short = 'q', long, default_value_t = 20)]
	min_base_quality: u8,

	/// Minimum mapping quality
	#[arg(short = 'm', long, default_value_t = 10)]
	min_map_quality: u8,

	/// SAM flags that must be present
	#[arg(short = 'f', default_value_t = 0)]
	required_flags: u16,

	/// SAM flags that, if present, skip read
	#[arg(short = 'F', default_value_t = 0)]
	filter_flags: u16,

	/// SAM flags that, if all present, skip read
	#[arg(short = 'G', default_value_t = 0)]
	excl_flags: u16,

	/// How to handle read-pair overlapping positions: cut (split at midpoint) or stretch (cubic spline scale)
	#[arg(long, value_enum, default_value = "cut")]
	overlap_mode: OverlapMode,

	/// Write output to file instead of stdout
	#[arg(short = 'o', long)]
	output_file: Option<String>,

	/// Methylation mode for C>T and G>A changes: if true, treat them as C>C and G>G to avoid confounding by bisulfite conversion
	#[arg(long, default_value_t = false)]
	methylation_mode: bool,

	/// minimum fragment length to be accepted for insert mode counting;
	#[arg(long, default_value_t = 0)]
	min_fragment_length: usize,

	/// maximum fragment length for insert mode counting
	#[arg(long, default_value_t = 1500)]
	max_fragment_length: usize,

	/// read position or fragment (insert) position mode
	#[arg(long, value_enum, default_value = "insert")]
	position_mode: PositionMode,
}

fn configure_thread_pool(threads: usize) {
	if threads > 0 {
		rayon::ThreadPoolBuilder::new()
			.num_threads(threads)
			.build_global()
			.expect("Failed to configure rayon thread pool");
	}
}

fn build_tid_map_and_regions(
	bam_path: &str,
	region_size: usize,
) -> (Arc<HashMap<i32, String>>, Vec<(i32, i64, i64)>) {
	let bam = Reader::from_path(bam_path).expect("Failed to open BAM file");

	let tid_to_name: Arc<HashMap<i32, String>> = Arc::new(
		(0..bam.header().target_count())
			.map(|i| {
				(
					i as i32,
					String::from_utf8_lossy(bam.header().tid2name(i)).to_string(),
				)
			})
			.collect(),
	);

	let mut regions = Vec::new();
	for &tid in tid_to_name.keys() {
		let chr_len = bam.header().target_len(tid as u32).unwrap_or(0) as i64;
		let mut start = 0i64;
		while start < chr_len {
			let end = std::cmp::min(start + region_size as i64, chr_len);
			regions.push((tid, start, end));
			start = end;
		}
	}

	(tid_to_name, regions)
}

/// Parse the MC tag and return the mate's exclusive end position, or `None`
/// when the MC tag is absent, the record is unpaired / mate-unmapped, or the
/// reads are on different contigs.
fn mc_mate_end(record: &Record) -> Option<i64> {
	if !record.is_paired() || record.is_mate_unmapped() {
		return None;
	}
	if record.tid() < 0 || record.mtid() < 0 || record.tid() != record.mtid() {
		return None;
	}
	let mate_start = record.mpos();
	let mate_span = match record.aux(b"MC".as_ref()) {
		Ok(Aux::String(mc)) => i64::try_from(cigar_reference_span(mc)?).ok()?,
		_ => return None,
	};
	Some(mate_start + mate_span)
}

fn estimated_fragment_length(record: &Record, mate_end: Option<i64>) -> Option<usize> {
	let mate_end = mate_end?;

	let read_start = record.pos();
	let read_end = record.cigar().end_pos();
	let mate_start = record.mpos();

	let fragment_start = std::cmp::min(read_start, mate_start);
	let fragment_end = std::cmp::max(read_end, mate_end);
	if fragment_end <= fragment_start {
		return None;
	}

	usize::try_from(fragment_end - fragment_start).ok()
}

fn should_skip_record(record: &Record, args: &Args, mate_end: Option<i64>) -> bool {
	if record.is_unmapped() {
		return true;
	}

	if record.mapq() < args.min_map_quality {
		return true;
	}

	let flags = record.flags();
	if args.required_flags != 0 && (flags & args.required_flags) != args.required_flags {
		return true;
	}
	if args.filter_flags != 0 && (flags & args.filter_flags) != 0 {
		return true;
	}
	if args.excl_flags != 0 && (flags & args.excl_flags) == args.excl_flags {
		return true;
	}

	let fragment_length = estimated_fragment_length(record, mate_end)
		.unwrap_or_else(|| record.insert_size().abs() as usize);

//	if fragment_length <= args.min_fragment_length || fragment_length > args.max_fragment_length {
//		eprintln!("skipping read: {}\tframent_length: {}\tmax and min set are {} and {}", String::from_utf8_lossy(record.qname()), fragment_length, args.min_fragment_length, args.max_fragment_length);
//		return true;
//	}

	false
}

fn record_read_num(record: &Record) -> u8 {
	if record.is_first_in_template() {
		1
	} else if record.is_last_in_template() {
		2
	} else {
		1
	}
}

fn read_is_first_in_reference(record: &Record) -> bool {
	if !record.is_paired() || record.is_mate_unmapped() {
		return true;
	}

	let read_coord = (record.tid(), record.pos());
	let mate_coord = (record.mtid(), record.mpos());

	if read_coord < mate_coord {
		true
	} else if read_coord > mate_coord {
		false
	} else {
		record_read_num(record) == 1
	}
}

fn insert_mode_read_position(
	is_first: bool,
	read_pos: usize,
	seq_len: usize,
	max_read_len: usize,
	stretch: bool,
	fragment_len: Option<usize>,
) -> usize {
	let short_fragment = fragment_len.map_or(false, |fl| fl <= max_read_len);
	let trailing_bases = seq_len - (read_pos + 1);

	if stretch && seq_len > 1 {
		// Cubic smooth-step: s = t²(3 − 2t), maps [0,1]→[0,1] with zero
		// derivatives at both ends, stretching each read to fill its half
		// of the virtual insert space regardless of actual read length.
		if is_first {
			let t = read_pos as f64 / (seq_len - 1) as f64;
			let s = t * t * (3.0 - 2.0 * t);
			1 + (s * (max_read_len - 1) as f64).round() as usize
		} else {
			let trailing = seq_len - 1 - read_pos;
			let t = trailing as f64 / (seq_len - 1) as f64;
			let s = t * t * (3.0 - 2.0 * t);
			2 * max_read_len - (s * (max_read_len - 1) as f64).round() as usize
		}
	} else if is_first {
		if read_pos >= seq_len / 2 && short_fragment {
			2 * max_read_len - (seq_len - read_pos) + 1
		} else {
			read_pos + 1
		}
	} else if read_pos < seq_len / 2 && short_fragment {
		read_pos + 1
	} else {
		2 * max_read_len - trailing_bases
	}
}

fn oriented_read_position(pos: usize, seq_len: usize, is_reverse: bool, max_read_len: usize) -> usize {
	let half = (seq_len + 1) /2;

	let returnit = match (is_reverse, pos <= half) {
		(true, true) => max_read_len - pos,
		(true, false) => seq_len - pos,
		(false, false) => pos + (max_read_len - seq_len) + 1,
		(false, true) => pos + 1,
	};

	returnit
}

fn base_position_for_mode(
	position_mode: PositionMode,
	read_pos: usize,
	seq_len: usize,
	is_reverse: bool,
	is_first_in_reference: bool,
	max_read_len: usize,
	stretch: bool,
	fragment_len: Option<usize>,
) -> usize {

	match position_mode {
		PositionMode::Read => oriented_read_position(read_pos, seq_len, is_reverse, max_read_len),
		PositionMode::Insert => {
			insert_mode_read_position(is_first_in_reference, read_pos, seq_len, max_read_len, stretch, fragment_len)
		}
	}
}

fn cigar_reference_span(cigar: &str) -> Option<usize> {
	let mut ref_bases = 0usize;
	let mut run_len = 0usize;

	for b in cigar.bytes() {
		if b.is_ascii_digit() {
			run_len = run_len * 10 + (b - b'0') as usize;
			continue;
		}

		match b {
			b'M' | b'=' | b'X' | b'D' | b'N' => {
				ref_bases += run_len;
			}
			b'I' | b'S' | b'H' | b'P' => {}
			_ => return None,
		}
		run_len = 0;
	}

	if run_len != 0 {
		return None;
	}

	Some(ref_bases)
}

fn softclip_side_comparisons(
	record: &Record,
	ref_seq: &[u8],
	read_start: usize,
	clip_len: usize,
	ref_start: usize,
) -> Vec<SoftclipComparison> {
	if clip_len == 0 || read_start + clip_len > record.seq().len() || ref_start + clip_len > ref_seq.len() {
		return Vec::new();
	}

	let read_seq = record.seq();
	let mut comparisons = Vec::with_capacity(clip_len);

	for offset in 0..clip_len {
		let read_index = read_start + offset;
		let ref_index = ref_start + offset;

		let Some(read_base) = base_to_char(read_seq[read_index]) else {
			continue;
		};
		let Some(ref_base) = base_to_char(ref_seq[ref_index].to_ascii_uppercase()) else {
			continue;
		};

		comparisons.push(SoftclipComparison {
			read_pos: read_index,
			ref_pos: ref_index,
			read_base,
			ref_base,
		});
	}

	comparisons
}

fn softclip_identity(comparisons: &[SoftclipComparison]) -> Option<f64> {
	if comparisons.is_empty() {
		None
	} else {
		let matches = comparisons
			.iter()
			.filter(|comparison| comparison.read_base == comparison.ref_base)
			.count();
		Some(matches as f64 / comparisons.len() as f64)
	}
}

fn qualifying_softclip_comparisons(
	record: &Record,
	ref_seq: &[u8],
	min_identity: f64,
) -> Vec<SoftclipComparison> {
	let cigar = record.cigar();
	let aligned_start = record.pos();
	let aligned_end = cigar.end_pos();
	let mut comparisons = Vec::new();

	if let Some(rust_htslib::bam::record::Cigar::SoftClip(len)) = cigar.iter().next() {
		if aligned_start >= *len as i64 {
			let side = softclip_side_comparisons(
				record,
				ref_seq,
				0,
				*len as usize,
				(aligned_start - *len as i64) as usize,
			);
			if softclip_identity(&side).is_some_and(|identity| identity >= min_identity) {
				comparisons.extend(side);
			}
		}
	}

	if let Some(rust_htslib::bam::record::Cigar::SoftClip(len)) = cigar.iter().last() {
		let side = softclip_side_comparisons(
			record,
			ref_seq,
			record.seq().len().saturating_sub(*len as usize),
			*len as usize,
			aligned_end as usize,
		);
		if softclip_identity(&side).is_some_and(|identity| identity >= min_identity) {
			comparisons.extend(side);
		}
	}

	comparisons
}

fn softclip_has_min_reference_identity(
	record: &Record,
	reference: &ReferenceGenome,
	tid_to_name: &HashMap<i32, String>,
	min_identity: f64,
) -> bool {
	let tid = record.tid();
	if tid < 0 {
		return false;
	}

	let Some(chr_name) = tid_to_name.get(&tid) else {
		return false;
	};
	let Some(ref_seq) = reference.get(chr_name.as_str()) else {
		return false;
	};

	!qualifying_softclip_comparisons(record, ref_seq, min_identity).is_empty()
}

fn softclip_has_reference_identity_at_least_66(
	record: &Record,
	reference: &ReferenceGenome,
	tid_to_name: &HashMap<i32, String>,
) -> bool {
	softclip_has_min_reference_identity(record, reference, tid_to_name, 0.66)
}

fn overlap_interval(record: &Record, mate_end: Option<i64>) -> Option<(usize, usize)> {
	let mate_end = mate_end?;

	let read_start = record.pos();
	let read_end = record.cigar().end_pos();
	let mate_start = record.mpos();

	let ov_start = std::cmp::max(read_start, mate_start);
	let ov_end = std::cmp::min(read_end, mate_end);
	if ov_start < ov_end {
		Some((ov_start as usize, ov_end as usize))
	} else {
		None
	}
}

fn build_base_change(read_num: u8, ref_base: char, read_base: char, is_reverse: bool, methylation_mode: bool) -> String {
	let mut base_change = match methylation_mode {
		false => format!("{}>{}", ref_base, read_base),
		true => match (read_num, ref_base, read_base, is_reverse) {
			(1, 'C', 'T', false) 
			| (2, 'G', 'A', false)
			| (1, 'G', 'A', true)
			| (2, 'C', 'T', true) => format!("{}>{}", ref_base, ref_base),
			_ => match is_reverse {
				true => format!("{}>{}", complement(ref_base), complement(read_base)),
				false => format!("{}>{}", ref_base, read_base),
			},
		},
	};

	if is_reverse && !methylation_mode {
		base_change = format!("{}>{}", complement(ref_base), complement(read_base));
	}

	base_change
}

fn compare_record_to_reference(
	record: &Record,
	reference: &ReferenceGenome,
	tid_to_name: &HashMap<i32, String>,
	min_base_quality: u8,
	max_read_len: usize,
	position_mode: PositionMode,
	overlap_mode: &OverlapMode,
	methylation_mode: bool,
	mate_end: Option<i64>,
	local_counts: &mut HashMap<InsertKey, usize>,
) {
	let tid = record.tid();
	if tid < 0 {
		return;
	}

	let Some(chr_name) = tid_to_name.get(&tid) else {
		return;
	};
	let Some(ref_seq) = reference.get(chr_name.as_str()) else {
		return;
	};

	let seq = record.seq();
	let qual = record.qual();
	let seq_len = seq.len();
	let read_num = record_read_num(record);
	let is_first_in_reference = read_is_first_in_reference(record);
	let overlap = overlap_interval(record, mate_end);
	let softclip_comparisons = qualifying_softclip_comparisons(record, ref_seq, 0.66);
	let stretch = *overlap_mode == OverlapMode::Stretch;

	let mut read_pos = 0usize;
	let mut ref_pos = record.pos() as usize;



	let seq_str: String = (0..seq_len)		
		.map(|i| base_to_char(seq[i]).unwrap_or('N'))
		.collect();




	use rust_htslib::bam::record::Cigar::*;
	for op in record.cigar().iter() {
		match op {
			Match(len) | Equal(len) | Diff(len) => {
				for i in 0..*len as usize {
					let rp = read_pos + i;
					let gp = ref_pos + i;

					if *overlap_mode == OverlapMode::Cut {
						if let Some((ov_start, ov_end)) = overlap {
							if gp >= ov_start && gp < ov_end && read_num == 2 {
								continue; // In overlap, count only read1 and skip read2.
							}
						}
					}

					if rp >= seq_len || rp >= qual.len() || gp >= ref_seq.len() {
						continue;
					}
					if qual[rp] < min_base_quality {
						continue;
					}

					let Some(read_base) = base_to_char(seq[rp]) else {
						continue;
					};

					let ref_byte = ref_seq[gp].to_ascii_uppercase();
					let Some(ref_base) = base_to_char(ref_byte) else {
						continue;
					};

					let reference_order = if is_first_in_reference { 1u8 } else { 2u8 };
					let ref_base_change = build_base_change(
						read_num, ref_base, read_base, record.is_reverse(), methylation_mode
					);

					let base_position = base_position_for_mode(
						position_mode,
						rp,
						seq_len,
						record.is_reverse(),
						is_first_in_reference,
						max_read_len,
						stretch,
						estimated_fragment_length(record, mate_end),
					);

					let key = InsertKey {
						base_change: ref_base_change,
						read_num,
						base_position,
						reference_order,
					};
					*local_counts.entry(key).or_insert(0) += 1;
				}

				read_pos += *len as usize;
				ref_pos += *len as usize;
			}
			Ins(len) | SoftClip(len) => {
				read_pos += *len as usize;
			}
			Del(len) | RefSkip(len) => {
				ref_pos += *len as usize;
			}
			HardClip(_) | Pad(_) => {}
		}
	}

	for comparison in softclip_comparisons {
		if comparison.read_pos >= qual.len() || qual[comparison.read_pos] < min_base_quality {
			continue;
		}

		let reference_order = if is_first_in_reference { 1u8 } else { 2u8 };
		let ref_base_change = build_base_change(
			read_num,
			comparison.ref_base,
			comparison.read_base,
			record.is_reverse(),
			methylation_mode,
		);
		let base_position = base_position_for_mode(
			position_mode,
			comparison.read_pos,
			seq_len,
			record.is_reverse(),
			is_first_in_reference,
			max_read_len,
			stretch,
			estimated_fragment_length(record, mate_end),
		);
		let key = InsertKey {
			base_change: ref_base_change,
			read_num,
			base_position,
			reference_order,
		};
		*local_counts.entry(key).or_insert(0) += 1;

		let _ = comparison.ref_pos;
	}
}

fn process_region(
	bam_path: &str,
	tid: i32,
	start: i64,
	end: i64,
	reference: &ReferenceGenome,
	tid_to_name: &HashMap<i32, String>,
	args: &Args,
	max_read_len: usize,
) -> (HashMap<InsertKey, usize>, usize) {
	let mut bam = IndexedReader::from_path(bam_path).expect("Failed to open indexed BAM");
	if let Err(error) = bam.fetch(FetchDefinition::Region(tid, start, end)) {
		eprintln!(
			"Warning: failed to fetch region tid {}:{}-{}: {}",
			tid, start, end, error
		);
		return (HashMap::new(), 0);
	}

	let mut local_counts: HashMap<InsertKey, usize> = HashMap::new();
	let mut local_read_count = 0usize;

	for result in bam.records() {
		let record = match result {
			Ok(r) => r,
			Err(_) => continue,
		};

		// Ensure each read is counted once across chunks.
		let read_start = record.pos();
		if read_start < start || read_start >= end {
			continue;
		}

		let mate_end = mc_mate_end(&record);

		if should_skip_record(&record, args, mate_end) {
			continue;
		}

		compare_record_to_reference(
			&record,
			reference,
			tid_to_name,
			args.min_base_quality,
			max_read_len,
			args.position_mode,
			&args.overlap_mode,
			args.methylation_mode,
			mate_end,
			&mut local_counts,
		);
		local_read_count += 1;
	}

	(local_counts, local_read_count)
}

fn write_output(
	counts: &HashMap<InsertKey, usize>,
	output_file: Option<&str>,
	position_mode: PositionMode,
) -> io::Result<()> {
	let position_label = match position_mode {
		PositionMode::Read => "read_position",
		PositionMode::Insert => "fragment_position",
	};

	let mut rows: Vec<(&InsertKey, &usize)> = counts.iter().collect();
	rows.sort_by(|(a, _), (b, _)| {
		a.reference_order
			.cmp(&b.reference_order)
			.then(a.read_num.cmp(&b.read_num))
			.then(a.base_position.cmp(&b.base_position))
			.then(a.base_change.cmp(&b.base_change))
	});

	if let Some(path) = output_file {
		let file = File::create(path)?;
		let mut writer = BufWriter::new(file);
		writeln!(
			writer,
			"base_change\tread_num\treference_order\t{}\tcount",
			position_label
		)?;
		for (key, count) in rows {
			writeln!(
				writer,
				"{}\t{}\t{}\t{}\t{}",
				key.base_change, key.read_num, key.reference_order, key.base_position, count
			)?;
		}
		return Ok(());
	}

	println!(
		"base_change\tread_num\treference_order\t{}\tcount",
		position_label
	);
	for (key, count) in rows {
		println!(
			"{}\t{}\t{}\t{}\t{}",
			key.base_change, key.read_num, key.reference_order, key.base_position, count
		);
	}
	Ok(())
}

fn main() {
	let args = Args::parse();
	eprintln!("using position mode {:?} and overlap mode {:?}\n", args.position_mode, args.overlap_mode);

	configure_thread_pool(args.threads);

	let max_read_len = compute_read_len_max_from_sample_bam(&args.bam_path, 10_000);
	eprintln!("Sampled max read length is {}", max_read_len);

	let reference = Arc::new(load_reference_genome(&args.reference_path));
	let (tid_to_name, regions) = build_tid_map_and_regions(&args.bam_path, args.region_size);
	eprintln!("Processing {} indexed regions", regions.len());

	let total_reads = AtomicUsize::new(0);
	let global_counts: Arc<Mutex<HashMap<InsertKey, usize>>> = Arc::new(Mutex::new(HashMap::new()));

	regions.par_iter().for_each(|(tid, start, end)| {
		let (region_counts, region_reads) = process_region(
			&args.bam_path,
			*tid,
			*start,
			*end,
			&reference,
			&tid_to_name,
			&args,
			max_read_len,
		);

		if !region_counts.is_empty() {
			let mut counts = global_counts.lock().expect("Lock poisoned");
			for (key, count) in region_counts {
				*counts.entry(key).or_insert(0) += count;
			}
		}

		let prev = total_reads.fetch_add(region_reads, Ordering::Relaxed);
		if (prev + region_reads) / 100_000 > prev / 100_000 {
			eprintln!("Processed {} reads", prev + region_reads);
		}
	});

	eprintln!("Total reads processed: {}", total_reads.load(Ordering::Relaxed));

	let counts = global_counts.lock().expect("Lock poisoned");
	write_output(&counts, args.output_file.as_deref(), args.position_mode)
		.expect("Failed to write output");
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn cigar_reference_span_counts_only_reference_consuming_ops() {
		assert_eq!(cigar_reference_span("10M2I3D4S5="), Some(18));
		assert_eq!(cigar_reference_span("7M1N2X3H"), Some(10));
	}

	#[test]
	fn cigar_reference_span_rejects_invalid_cigar_strings() {
		assert_eq!(cigar_reference_span("10M2"), None);
		assert_eq!(cigar_reference_span("10M2Z"), None);
	}

	#[test]
	fn softclip_identity_returns_none_for_empty_input() {
		assert_eq!(softclip_identity(&[]), None);
	}

	#[test]
	fn softclip_identity_computes_match_fraction() {
		let comparisons = vec![
			SoftclipComparison {
				read_pos: 0,
				ref_pos: 10,
				read_base: 'A',
				ref_base: 'A',
			},
			SoftclipComparison {
				read_pos: 1,
				ref_pos: 11,
				read_base: 'C',
				ref_base: 'T',
			},
			SoftclipComparison {
				read_pos: 2,
				ref_pos: 12,
				read_base: 'G',
				ref_base: 'G',
			},
		];

		assert_eq!(softclip_identity(&comparisons), Some(2.0 / 3.0));
	}

	#[test]
	fn insert_mode_read_position_handles_first_and_second_reads_without_stretch() {
		assert_eq!(insert_mode_read_position(true, 0, 5, 8, false, None), 1);
		assert_eq!(insert_mode_read_position(true, 4, 5, 8, false, None), 5);
		assert_eq!(insert_mode_read_position(false, 0, 5, 8, false, None), 12);
		assert_eq!(insert_mode_read_position(false, 4, 5, 8, false, None), 16);
	}

	#[test]
	fn insert_mode_read_position_stretch_hits_expected_endpoints() {
		assert_eq!(insert_mode_read_position(true, 0, 5, 8, true, None), 1);
		assert_eq!(insert_mode_read_position(true, 4, 5, 8, true, None), 8);
		assert_eq!(insert_mode_read_position(false, 0, 5, 8, true, None), 9);
		assert_eq!(insert_mode_read_position(false, 4, 5, 8, true, None), 16);
	}

	#[test]
	fn base_position_for_mode_dispatches_to_read_and_insert_modes() {
		assert_eq!(
			base_position_for_mode(PositionMode::Read, 1, 5, false, true, 8, false, None),
			2
		);
		assert_eq!(
			base_position_for_mode(PositionMode::Insert, 1, 5, false, true, 8, false, None),
			2
		);
		assert_eq!(
			base_position_for_mode(PositionMode::Insert, 1, 5, false, false, 8, false, None),
			13
		);
	}

	#[test]
	fn build_base_change_preserves_raw_mismatches_when_methylation_mode_is_off() {
		assert_eq!(build_base_change(1, 'C', 'T', false, false), "C>T");
		assert_eq!(build_base_change(2, 'G', 'A', true, false), "G>A");
	}

	#[test]
	fn build_base_change_collapses_bisulfite_signatures_in_methylation_mode() {
		assert_eq!(build_base_change(1, 'C', 'T', false, true), "C>C");
		assert_eq!(build_base_change(2, 'G', 'A', false, true), "G>G");
		assert_eq!(build_base_change(1, 'G', 'A', true, true), "G>G");
		assert_eq!(build_base_change(2, 'C', 'T', true, true), "C>C");
		assert_eq!(build_base_change(1, 'C', 'A', false, true), "C>A");
	}
}

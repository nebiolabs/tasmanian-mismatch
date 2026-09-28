use crate::io::with_output_writer;
use crate::types::InsertKey;
use crate::utils::{BASES, ratio, split_base_change};
use std::collections::{BTreeMap, HashMap};

/// One summary row for a genomic window: a mismatch class for one read, totaled over every
/// read/fragment position and reference order in the window.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowRow {
    pub read_num: u8,
    /// The mismatch class, `ref_base>read_base`.
    pub ref_base: char,
    pub read_base: char,
    /// Observations of this mismatch class.
    pub count: u64,
    /// Observations of any read base against this class's reference base (the rate's
    /// denominator, matching `normalize_mismatch_counts` but pooled over positions).
    pub ref_total: u64,
    /// Observations of any A/C/G/T base on this read in the window (`rate_all`'s denominator),
    /// so the classes' `rate_all` values add up to the window's overall mismatch rate.
    pub all_bases: u64,
}

impl WindowRow {
    pub fn rate(&self) -> f64 {
        ratio(self.count, self.ref_total)
    }

    pub fn rate_all(&self) -> f64 {
        ratio(self.count, self.all_bases)
    }
}

/// A window's coordinates (contig `tid`, 0-based, half-open) and its summary rows.
pub struct WindowSummary {
    pub tid: i32,
    pub start: i64,
    pub end: i64,
    pub rows: Vec<WindowRow>,
}

/// Collapse one window's position-level counts into per-(read_num, base_change) mismatch rows.
///
/// Every A/C/G/T mismatch class whose reference base was observed gets a row, with a zero
/// count if it never occurred, so windows can be compared without filling gaps. Other
/// observed mismatch classes (e.g. involving `N`) are kept as they appear. Matches
/// (`C>C`) only contribute to `ref_total`.
pub fn summarize_window_counts(counts: &HashMap<InsertKey, usize>) -> Vec<WindowRow> {
    let mut ref_totals: HashMap<(u8, char), u64> = HashMap::new();
    let mut mismatches: BTreeMap<(u8, char, char), u64> = BTreeMap::new();

    for (key, &count) in counts {
        let Some((ref_base, read_base)) = split_base_change(&key.base_change) else {
            continue;
        };
        *ref_totals.entry((key.read_num, ref_base)).or_insert(0) += count as u64;
        if ref_base != read_base {
            *mismatches
                .entry((key.read_num, ref_base, read_base))
                .or_insert(0) += count as u64;
        }
    }

    let mut all_bases: HashMap<u8, u64> = HashMap::new();
    for (&(read_num, ref_base), &total) in &ref_totals {
        if !BASES.contains(&ref_base) {
            continue;
        }
        *all_bases.entry(read_num).or_insert(0) += total;
        for read_base in BASES.into_iter().filter(|&b| b != ref_base) {
            mismatches
                .entry((read_num, ref_base, read_base))
                .or_insert(0);
        }
    }

    mismatches
        .into_iter()
        .map(|((read_num, ref_base, read_base), count)| WindowRow {
            read_num,
            ref_base,
            read_base,
            count,
            ref_total: ref_totals[&(read_num, ref_base)],
            all_bases: all_bases.get(&read_num).copied().unwrap_or(0),
        })
        .collect()
}

/// Write window summaries as TSV, in the order given, naming contigs via `tid_to_name`.
pub fn write_window_output(
    windows: &[WindowSummary],
    tid_to_name: &HashMap<i32, String>,
    output_file: Option<&str>,
) -> std::io::Result<()> {
    with_output_writer(output_file, |w| {
        writeln!(
            w,
            "chrom\tstart\tend\tread_num\tbase_change\tcount\tref_total\trate\tall_bases\trate_all"
        )?;
        for window in windows {
            let chrom = &tid_to_name[&window.tid];
            for row in &window.rows {
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}>{}\t{}\t{}\t{:.6}\t{}\t{:.6}",
                    chrom,
                    window.start,
                    window.end,
                    row.read_num,
                    row.ref_base,
                    row.read_base,
                    row.count,
                    row.ref_total,
                    row.rate(),
                    row.all_bases,
                    row.rate_all()
                )?;
            }
        }
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::ReferenceOrder;

    fn key(base_change: &str, read_num: u8, base_position: usize) -> InsertKey {
        InsertKey {
            base_change: base_change.to_string(),
            read_num,
            base_position,
            reference_order: ReferenceOrder::First,
        }
    }

    fn row<'a>(rows: &'a [WindowRow], read_num: u8, base_change: &str) -> &'a WindowRow {
        let (ref_base, read_base) = split_base_change(base_change).expect("valid base change");
        rows.iter()
            .find(|r| r.read_num == read_num && (r.ref_base, r.read_base) == (ref_base, read_base))
            .unwrap_or_else(|| panic!("no {base_change} row for read {read_num}"))
    }

    #[test]
    fn pools_positions_and_fills_unobserved_mismatch_classes() {
        let counts = HashMap::from([
            (key("C>T", 1, 1), 3),
            (key("C>T", 1, 2), 1),
            (key("C>C", 1, 1), 90),
            (key("C>C", 1, 2), 6),
            (key("G>A", 2, 1), 2),
            (key("G>G", 2, 1), 48),
        ]);
        let rows = summarize_window_counts(&counts);

        let ct = row(&rows, 1, "C>T");
        assert_eq!((ct.count, ct.ref_total), (4, 100));
        assert!((ct.rate() - 0.04).abs() < 1e-12);
        // Read 1 only observed C bases, so all bases = C bases; read 2 only G.
        assert_eq!(ct.all_bases, 100);
        assert_eq!(row(&rows, 2, "G>A").all_bases, 50);
        // Unobserved C mismatches on read 1 get explicit zero rows over the same total.
        assert_eq!(row(&rows, 1, "C>A").count, 0);
        assert_eq!(row(&rows, 1, "C>G").ref_total, 100);
        assert_eq!(row(&rows, 2, "G>A").ref_total, 50);

        // Matches only feed the denominator, and unobserved reference bases get no rows.
        assert!(rows.iter().all(|r| r.ref_base != r.read_base));
        assert!(rows.iter().all(|r| r.ref_base != 'A'));
        assert_eq!(rows.len(), 6);
    }

    #[test]
    fn empty_window_has_no_rows() {
        assert!(summarize_window_counts(&HashMap::new()).is_empty());
    }

    #[test]
    fn rate_all_divides_by_every_base_on_the_read() {
        let counts = HashMap::from([
            (key("C>T", 1, 1), 2),
            (key("C>C", 1, 1), 98),
            (key("A>G", 1, 1), 3),
            (key("A>A", 1, 1), 297),
        ]);
        let rows = summarize_window_counts(&counts);
        let ct = row(&rows, 1, "C>T");
        let ag = row(&rows, 1, "A>G");
        assert_eq!((ct.all_bases, ag.all_bases), (400, 400));
        assert!((ct.rate_all() - 2.0 / 400.0).abs() < 1e-12);
        // The read's classes' rate_all values sum to its overall mismatch rate, 5 / 400.
        let total: f64 = rows
            .iter()
            .filter(|r| r.read_num == 1)
            .map(WindowRow::rate_all)
            .sum();
        assert!((total - 5.0 / 400.0).abs() < 1e-12);
    }
}

//! Tasmanian-specific SAM record wrapper

use super::flags::{OntFlags, ProperPairFlags, ReadNumber, Strand};
use anyhow::Result;
use noodles::sam::{
    alignment::{
        record::cigar::op::Kind as CigarOp,
        Record,
    },
    Header,
};

/// Extended SAM record with Tasmanian-specific metadata
#[derive(Debug, Clone)]
pub struct TasmanianRecord {
    pub name: String,
    pub flag: u16,
    pub chrom: String,
    pub start: u64, // 0-based
    pub mapq: u8,
    pub cigar_str: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub tlen: i64,

    // Intersection metadata (set by intersection processor)
    pub tm_tag: Option<(u64, u64)>,    // Intersection bounds from tm:Z tag
    pub confidence_bases: Option<u32>, // Complement bases from tc:i tag

    // Computed fields
    seq_len: usize,
}

impl TasmanianRecord {
    /// Create a TasmanianRecord from a noodles Record
    pub fn from_noodles<R: Record>(record: &R, header: &Header) -> Result<Self> {
        let name = record
            .name()
            .map(|n| {
                let bytes: Vec<u8> = n.iter().copied().collect();
                String::from_utf8_lossy(&bytes).to_string()
            })
            .unwrap_or_default();

        let flag = record.flags()?.bits();

        let chrom = record
            .reference_sequence_id(header)
            .transpose()?
            .and_then(|id| header.reference_sequences().get_index(id))
            .map(|(name, _)| name.to_string())
            .unwrap_or_default();

        let start = record
            .alignment_start()
            .transpose()?
            .map(|pos| pos.get() as u64 - 1) // Convert to 0-based
            .unwrap_or(0);

        let mapq = record.mapping_quality().transpose()?.map_or(255, |q| q.get());

        // Build CIGAR string
        let cigar = record.cigar();
        let cigar_str = cigar
            .iter()
            .filter_map(|op_result| {
                op_result.ok().map(|op| format!("{}{}", op.len(), cigar_op_char(op.kind())))
            })
            .collect::<Vec<_>>()
            .join("");

        let seq: Vec<u8> = record
            .sequence()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let seq_len = seq.len();

        let qual: Vec<u8> = record
            .quality_scores()
            .iter()
            .filter_map(|q| q.ok().map(|v| v + 33))
            .collect();

        let tlen = record.template_length()?;

        // Parse Tasmanian tags
        let (tm_tag, confidence_bases) = parse_tasmanian_tags(record)?;

        Ok(TasmanianRecord {
            name,
            flag,
            chrom,
            start,
            mapq,
            cigar_str,
            seq,
            qual,
            tlen: tlen as i64,
            tm_tag,
            confidence_bases,
            seq_len,
        })
    }

    /// Get sequence length
    pub fn seq_len(&self) -> usize {
        self.seq_len
    }

    /// Get end position (0-based, exclusive)
    pub fn end(&self) -> u64 {
        self.start + self.seq_len as u64
    }

    /// Check if this is a proper pair for standard paired-end data
    pub fn is_proper_pair(&self) -> bool {
        ProperPairFlags::is_proper_pair(self.flag)
    }

    /// Check if this is a valid ONT read
    pub fn is_valid_ont(&self) -> bool {
        OntFlags::is_valid_ont(self.flag)
    }

    /// Get the read number (1 or 2)
    pub fn read_number(&self, ont_mode: bool) -> Option<ReadNumber> {
        if ont_mode {
            // ONT reads are all "read 1"
            if OntFlags::is_valid_ont(self.flag) {
                Some(ReadNumber::First)
            } else {
                None
            }
        } else {
            ProperPairFlags::read_number(self.flag)
        }
    }

    /// Get the strand
    pub fn strand(&self, ont_mode: bool) -> Option<Strand> {
        if ont_mode {
            OntFlags::strand(self.flag)
        } else {
            ProperPairFlags::strand(self.flag)
        }
    }

    /// Check if CIGAR is valid (not '*')
    pub fn has_valid_cigar(&self) -> bool {
        !self.cigar_str.is_empty() && self.cigar_str != "*"
    }

    /// Check if this read has indels in its CIGAR
    pub fn has_indels(&self) -> bool {
        self.cigar_str.contains('I') || self.cigar_str.contains('D')
    }

    /// Get absolute template length
    pub fn abs_tlen(&self) -> u64 {
        self.tlen.unsigned_abs()
    }

    /// Check if the read has the tasmanian intersection tag
    pub fn has_tasmanian_tag(&self) -> bool {
        self.tm_tag.is_some()
    }

    /// Check if a position is in the intersection region
    pub fn is_in_intersection(&self, pos: usize) -> bool {
        if let Some((start, end)) = self.tm_tag {
            pos >= start as usize && pos < end as usize
        } else {
            false
        }
    }
}

/// Parse Tasmanian-specific tags from a SAM record
fn parse_tasmanian_tags<R: Record>(record: &R) -> Result<(Option<(u64, u64)>, Option<u32>)> {
    use noodles::sam::alignment::record::data::field::Value;

    let mut tm_tag = None;
    let mut tc_tag = None;

    let data = record.data();
    for result in data.iter() {
        let (tag, value) = result?;
        // Tag is [u8; 2], use as_ref() to get slice
        let tag_bytes: &[u8] = tag.as_ref();

        if tag_bytes == b"tm" {
            // tm:Z:a.b;start.end
            if let Value::String(s) = value {
                let s_bytes: Vec<u8> = s.iter().copied().collect();
                let s = String::from_utf8_lossy(&s_bytes);
                // Parse "a.b;..." format - we only need a.b
                if let Some(ab_part) = s.split(';').next() {
                    let parts: Vec<&str> = ab_part.split('.').collect();
                    if parts.len() >= 2 {
                        if let (Ok(a), Ok(b)) = (parts[0].parse::<u64>(), parts[1].parse::<u64>()) {
                            tm_tag = Some((a, b));
                        }
                    }
                }
            }
        } else if tag_bytes == b"tc" {
            // tc:i:N
            match value {
                Value::Int32(n) => {
                    tc_tag = Some(n as u32);
                }
                Value::UInt32(n) => {
                    tc_tag = Some(n);
                }
                Value::Int16(n) => {
                    tc_tag = Some(n as u32);
                }
                Value::UInt16(n) => {
                    tc_tag = Some(n as u32);
                }
                Value::Int8(n) => {
                    tc_tag = Some(n as u32);
                }
                Value::UInt8(n) => {
                    tc_tag = Some(n as u32);
                }
                _ => {}
            }
        }
    }

    Ok((tm_tag, tc_tag))
}

/// Convert CIGAR operation kind to character
fn cigar_op_char(op: CigarOp) -> char {
    match op {
        CigarOp::Match => 'M',
        CigarOp::Insertion => 'I',
        CigarOp::Deletion => 'D',
        CigarOp::Skip => 'N',
        CigarOp::SoftClip => 'S',
        CigarOp::HardClip => 'H',
        CigarOp::Pad => 'P',
        CigarOp::SequenceMatch => '=',
        CigarOp::SequenceMismatch => 'X',
    }
}

/// Parse CIGAR string into operations
pub fn parse_cigar(cigar: &str) -> Vec<(CigarOp, u32)> {
    let mut ops = Vec::new();
    let mut num = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            let len: u32 = num.parse().unwrap_or(0);
            num.clear();

            let op = match c {
                'M' => CigarOp::Match,
                'I' => CigarOp::Insertion,
                'D' => CigarOp::Deletion,
                'N' => CigarOp::Skip,
                'S' => CigarOp::SoftClip,
                'H' => CigarOp::HardClip,
                'P' => CigarOp::Pad,
                '=' => CigarOp::SequenceMatch,
                'X' => CigarOp::SequenceMismatch,
                _ => continue,
            };

            ops.push((op, len));
        }
    }

    ops
}

/// Compute reference and sequence index mappings from CIGAR
///
/// Returns (seq_mask, ref_mask) where each mask indicates which positions to use
pub fn compute_cigar_mappings(cigar: &str, seq_len: usize) -> (Vec<bool>, Vec<bool>) {
    let ops = parse_cigar(cigar);
    let mut seq_mask = vec![true; seq_len];
    let mut ref_mask = vec![true; seq_len];
    let mut pos = 0usize;

    for (op, len) in ops {
        let len = len as usize;
        match op {
            CigarOp::Match | CigarOp::SequenceMatch | CigarOp::SequenceMismatch => {
                pos += len;
            }
            CigarOp::Insertion => {
                // Insertion: bases in seq but not in ref
                for i in pos..pos.saturating_add(len).min(seq_mask.len()) {
                    seq_mask[i] = false;
                }
                pos += len;
            }
            CigarOp::Deletion | CigarOp::Skip => {
                // Deletion: bases in ref but not in seq
                for i in pos..pos.saturating_add(len).min(ref_mask.len()) {
                    ref_mask[i] = false;
                }
                pos += len;
            }
            CigarOp::SoftClip => {
                pos += len;
            }
            CigarOp::HardClip | CigarOp::Pad => {
                // Hard clips don't consume sequence
            }
        }
    }

    (seq_mask, ref_mask)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cigar() {
        let ops = parse_cigar("10M5I10M");
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], (CigarOp::Match, 10));
        assert_eq!(ops[1], (CigarOp::Insertion, 5));
        assert_eq!(ops[2], (CigarOp::Match, 10));
    }

    #[test]
    fn test_parse_cigar_softclip() {
        let ops = parse_cigar("5S90M5S");
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], (CigarOp::SoftClip, 5));
        assert_eq!(ops[1], (CigarOp::Match, 90));
        assert_eq!(ops[2], (CigarOp::SoftClip, 5));
    }
}

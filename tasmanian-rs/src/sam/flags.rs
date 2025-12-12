//! SAM flag interpretation

/// Strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// Read number in pair
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadNumber {
    First,
    Second,
}

impl ReadNumber {
    /// Get the 1-based read number for output
    pub fn as_number(&self) -> u8 {
        match self {
            ReadNumber::First => 1,
            ReadNumber::Second => 2,
        }
    }
}

/// SAM flags for properly paired reads
/// These are the flags we care about for artifact analysis
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProperPairFlags;

impl ProperPairFlags {
    /// Flag 99: read 1, forward strand, mate reverse
    pub const READ1_FWD: u16 = 99;
    /// Flag 163: read 2, forward strand, mate reverse
    pub const READ2_FWD: u16 = 163;
    /// Flag 83: read 1, reverse strand, mate forward
    pub const READ1_REV: u16 = 83;
    /// Flag 147: read 2, reverse strand, mate forward
    pub const READ2_REV: u16 = 147;

    /// Check if a flag represents a proper pair
    pub fn is_proper_pair(flag: u16) -> bool {
        matches!(
            flag,
            Self::READ1_FWD | Self::READ2_FWD | Self::READ1_REV | Self::READ2_REV
        )
    }

    /// Get the read number from a flag
    pub fn read_number(flag: u16) -> Option<ReadNumber> {
        match flag {
            Self::READ1_FWD | Self::READ1_REV => Some(ReadNumber::First),
            Self::READ2_FWD | Self::READ2_REV => Some(ReadNumber::Second),
            _ => None,
        }
    }

    /// Get the strand from a flag
    pub fn strand(flag: u16) -> Option<Strand> {
        match flag {
            Self::READ1_FWD | Self::READ2_FWD => Some(Strand::Forward),
            Self::READ1_REV | Self::READ2_REV => Some(Strand::Reverse),
            _ => None,
        }
    }
}

/// ONT (Oxford Nanopore) specific flags
#[derive(Debug, Clone, Copy)]
pub struct OntFlags;

impl OntFlags {
    /// Primary alignment, forward
    pub const PRIMARY_FWD: u16 = 0;
    /// Supplementary alignment, forward
    pub const SUPP_FWD: u16 = 2048;
    /// Primary alignment, reverse
    pub const PRIMARY_REV: u16 = 16;
    /// Supplementary alignment, reverse
    pub const SUPP_REV: u16 = 2064;

    /// Check if this is a valid ONT alignment
    pub fn is_valid_ont(flag: u16) -> bool {
        matches!(
            flag,
            Self::PRIMARY_FWD | Self::SUPP_FWD | Self::PRIMARY_REV | Self::SUPP_REV
        )
    }

    /// Get the strand from an ONT flag
    pub fn strand(flag: u16) -> Option<Strand> {
        match flag {
            Self::PRIMARY_FWD | Self::SUPP_FWD => Some(Strand::Forward),
            Self::PRIMARY_REV | Self::SUPP_REV => Some(Strand::Reverse),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_proper_pair_detection() {
        assert!(ProperPairFlags::is_proper_pair(99));
        assert!(ProperPairFlags::is_proper_pair(163));
        assert!(ProperPairFlags::is_proper_pair(83));
        assert!(ProperPairFlags::is_proper_pair(147));
        assert!(!ProperPairFlags::is_proper_pair(0));
        assert!(!ProperPairFlags::is_proper_pair(4)); // unmapped
    }

    #[test]
    fn test_read_number() {
        assert_eq!(
            ProperPairFlags::read_number(99),
            Some(ReadNumber::First)
        );
        assert_eq!(
            ProperPairFlags::read_number(163),
            Some(ReadNumber::Second)
        );
        assert_eq!(
            ProperPairFlags::read_number(83),
            Some(ReadNumber::First)
        );
        assert_eq!(
            ProperPairFlags::read_number(147),
            Some(ReadNumber::Second)
        );
    }

    #[test]
    fn test_strand() {
        assert_eq!(ProperPairFlags::strand(99), Some(Strand::Forward));
        assert_eq!(ProperPairFlags::strand(83), Some(Strand::Reverse));
    }
}

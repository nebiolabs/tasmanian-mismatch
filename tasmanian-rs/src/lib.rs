//! Tasmanian - A tool for analysis of reference mismatches in high throughput sequencing data
//!
//! Unlike other tools, Tasmanian can evaluate portions of reads that overlap with specified
//! genomic regions (e.g., repeats) separately from non-overlapping portions.

pub mod bed;
pub mod cli;
pub mod error_table;
pub mod reference;
pub mod sam;

mod analysis;
mod intersection;
mod output;

pub use analysis::{analyze, AnalysisConfig};
pub use bed::BedIndex;
pub use error_table::ErrorTables;
pub use intersection::IntersectionProcessor;
pub use output::write_csv;
pub use reference::Reference;

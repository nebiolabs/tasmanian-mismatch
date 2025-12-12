//! CSV output generation
//!
//! Outputs mismatch tables in the format expected by downstream analysis.

use crate::error_table::{mismatch_column_name, ErrorTables, MismatchTypeIter};
use anyhow::Result;
use std::io::Write;

/// Write error tables to CSV format
///
/// Output format:
/// read,position,Ia_a,Ia_t,...,Ca_a,...,Na_a,...,cIa_a,...,cCa_a,...
///
/// Where:
/// - I = intersection
/// - C = complement
/// - N = non-intersecting (unrelated)
/// - cI = confident intersection
/// - cC = confident complement
pub fn write_csv<W: Write>(writer: &mut W, tables: &ErrorTables, read_length: usize) -> Result<()> {
    // Write header
    let mut header = vec!["read".to_string(), "position".to_string()];

    // Add columns for each table
    for prefix in ["I", "C", "N", "cI", "cC"] {
        for (ref_base, alt_base) in MismatchTypeIter::new() {
            header.push(mismatch_column_name(prefix, ref_base, alt_base));
        }
    }

    writeln!(writer, "{}", header.join(","))?;

    // Write data rows
    for read_num in [1u8, 2u8] {
        for pos in 0..read_length {
            let mut row = vec![read_num.to_string(), (pos + 1).to_string()];

            // Intersection table
            add_table_columns(&mut row, &tables.intersection, read_num, pos);

            // Complement table
            add_table_columns(&mut row, &tables.complement, read_num, pos);

            // Unrelated table
            add_table_columns(&mut row, &tables.unrelated, read_num, pos);

            // Confident intersection table
            add_table_columns(&mut row, &tables.intersection_confident, read_num, pos);

            // Confident complement table
            add_table_columns(&mut row, &tables.complement_confident, read_num, pos);

            writeln!(writer, "{}", row.join(","))?;
        }
    }

    Ok(())
}

/// Add columns from an error table to a row
fn add_table_columns(
    row: &mut Vec<String>,
    table: &crate::error_table::ErrorTable,
    read_num: u8,
    pos: usize,
) {
    for (ref_base, alt_base) in MismatchTypeIter::new() {
        let count = table.get(read_num, pos, ref_base, alt_base);
        row.push(count.to_string());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_csv_header() {
        let tables = ErrorTables::new(false);
        let mut output = Vec::new();

        write_csv(&mut output, &tables, 2).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.lines().collect();

        // Check header
        assert!(lines[0].starts_with("read,position,"));
        assert!(lines[0].contains("Ia_a"));
        assert!(lines[0].contains("Ca_a"));
        assert!(lines[0].contains("Na_a"));
        assert!(lines[0].contains("cIa_a"));
        assert!(lines[0].contains("cCa_a"));

        // Check we have correct number of data rows (2 reads * 2 positions)
        assert_eq!(lines.len(), 5); // 1 header + 4 data rows
    }

    #[test]
    fn test_csv_data() {
        let mut tables = ErrorTables::new(false);

        // Add some test data
        tables.unrelated.increment(1, 0, b'A', b'T');
        tables.unrelated.increment(1, 0, b'A', b'T');
        tables.intersection.increment(1, 0, b'C', b'G');

        let mut output = Vec::new();
        write_csv(&mut output, &tables, 2).unwrap();

        let output_str = String::from_utf8(output).unwrap();

        // The unrelated A->T count should be 2
        // The intersection C->G count should be 1
        assert!(output_str.contains("1,1,")); // read 1, position 1
    }
}

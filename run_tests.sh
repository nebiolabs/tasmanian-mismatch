#!/bin/bash
# Quick test script for tasmanian-mismatch

set -e  # Exit on error

echo "=================================="
echo "Running Unit tests (rust functions)"
echo "=================================="
echo -e "[0/5] testing methods..."
cargo test --test unit_tests -- --nocapture

echo "=================================="
echo "Running functional tests (main and arguments)"
echo "=================================="

echo -e "\n[1/5] Creating minimal test data..."
./create_minimal_test.sh

# Check if test files were created
if [ ! -f "test_reads.bam" ]; then
    echo "Error: test_reads.bam not found"
    exit 1
fi

if [ ! -f "test_ref.fa" ]; then
    echo "Error: test_ref.fa not found"
    exit 1
fi

echo -e "\n[2/5] Building project..."
cargo build --release

echo -e "\n[3/5] Running basic test (no options)..."
cargo run --release test_reads.bam test_ref.fa > test_output_basic.csv 2> test_output_basic.log

# Check output contents. Compare to expected values.
LINES=$(wc -l < test_output_basic.csv)
echo "  Output lines: $LINES"

FIRST_ROW=$(sed -n '2p' test_output_basic.csv)
FIRST_VALUE=$(echo $FIRST_ROW | cut -d',' -f3)  # A>A column
echo "  First A>A value: $FIRST_VALUE (expected: 1)"

SECOND_ROW=$(sed -n '3p' test_output_basic.csv)
SECOND_VALUE=$(echo $SECOND_ROW | cut -d',' -f7)  # C>C column
echo "  Second C>C value: $SECOND_VALUE (expected: 1)"

if [ "$FIRST_VALUE" != "1" ]; then
    echo "ERROR: First A>A value mismatch. Expected 1, got $FIRST_VALUE"
    exit 1
fi

if [ "$SECOND_VALUE" != "1" ]; then
    echo "ERROR: Second C>C value mismatch. Expected 1, got $SECOND_VALUE"
    exit 1
fi

# Run methylation test
echo -e "\n[4/5] Running methylation mode test..."
cargo run --release test_reads.bam test_ref.fa --methylation --cpg-only > test_output_meth.csv 2> test_output_meth.log

LINES_METH=$(wc -l < test_output_meth.csv)
echo "  Output lines: $LINES_METH"

# Compare values from methylation output
FIRST_ROW_METH=$(sed -n '2p' test_output_meth.csv)
FIRST_VALUE_METH=$(echo $FIRST_ROW_METH | cut -d',' -f3)  # A>A column
echo "  First A>A value: $FIRST_VALUE_METH (expected: 1)"

SECOND_ROW_METH=$(sed -n '3p' test_output_meth.csv)
SECOND_VALUE_METH=$(echo $SECOND_ROW_METH | cut -d',' -f6)  # C>C column
echo "  Second C>C value: $SECOND_VALUE_METH (expected: 1)"

if [ "$FIRST_VALUE_METH" != "1" ]; then
    echo "ERROR: Methylation first A>A value mismatch. Expected 1, got $FIRST_VALUE_METH"
    exit 1
fi

if [ "$SECOND_VALUE_METH" != "0" ]; then
    echo "ERROR: Methylation second C>C value mismatch. Expected 1, got $SECOND_VALUE_METH"
    exit 1
fi

# Run with threads
echo -e "\n[5/5] Running with threading..."
cargo run --release test_reads.bam test_ref.fa --threads 4 > test_output_threads.csv 2> test_output_threads.log

LINES_THREADS=$(wc -l < test_output_threads.csv)
echo "  Output lines: $LINES_THREADS"

# Compare output values
FIRST_ROW_THREADS=$(sed -n '2p' test_output_threads.csv)
FIRST_VALUE_THREADS=$(echo $FIRST_ROW_THREADS | cut -d',' -f3)  # A>A column
echo "  First A>A value: $FIRST_VALUE_THREADS (expected: 1)"

SECOND_ROW_THREADS=$(sed -n '3p' test_output_threads.csv)
SECOND_VALUE_THREADS=$(echo $SECOND_ROW_THREADS | cut -d',' -f6)  # C>C column
echo "  Second C>C value: $SECOND_VALUE_THREADS (expected: 1)"

if [ "$FIRST_VALUE_THREADS" != "1" ]; then
    echo "ERROR: Threads first A>A value mismatch. Expected 1, got $FIRST_VALUE_THREADS"
    exit 1
fi

if [ "$SECOND_VALUE_THREADS" != "0" ]; then
    echo "ERROR: Threads second C>C value mismatch. Expected 1, got $SECOND_VALUE_THREADS"
    exit 1
fi

# Check outputs exist and have data
if [ $LINES -lt 2 ]; then
    echo "ERROR: Basic output has too few lines"
    exit 1
fi

if [ $LINES_METH -lt 2 ]; then
    echo "ERROR: Methylation output has too few lines"
    exit 1
fi

echo -e "\n=================================="
echo "All tests passed!"
echo "=================================="
echo "Test output files created:"
echo "  - test_output_basic.csv (and .log)"
echo "  - test_output_meth.csv (and .log)"
echo "  - test_output_threads.csv (and .log)"

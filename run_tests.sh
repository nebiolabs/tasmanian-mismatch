#!/bin/bash
# Quick test script for tasmanian-mismatch

set -e  # Exit on error

echo "=================================="
echo "Running tasmanian-mismatch tests"
echo "=================================="

# unit tests
echo -e "[0/5] testing methods..."
cargo test --test unit_tests -- --nocapture


echo "=================================="
echo "Running functional tests"
echo "=================================="


# Create minimal test data
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

# Build the project
echo -e "\n[2/5] Building project..."
cargo build --release

# Run basic test
echo -e "\n[3/5] Running basic test (no options)..."
cargo run --release test_reads.bam test_ref.fa > test_output_basic.csv 2> test_output_basic.log

# Count lines in output
LINES=$(wc -l < test_output_basic.csv)
echo "  Output lines: $LINES"

# Run methylation test
echo -e "\n[4/5] Running methylation mode test..."
cargo run --release test_reads.bam test_ref.fa --methylation --cpg-only > test_output_meth.csv 2> test_output_meth.log

LINES_METH=$(wc -l < test_output_meth.csv)
echo "  Output lines: $LINES_METH"

# Run with threads
echo -e "\n[5/5] Running with threading..."
cargo run --release test_reads.bam test_ref.fa --threads 4 > test_output_threads.csv 2> test_output_threads.log

LINES_THREADS=$(wc -l < test_output_threads.csv)
echo "  Output lines: $LINES_THREADS"

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

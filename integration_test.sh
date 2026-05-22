#!/bin/bash
# Integration test - runs the tool and validates output

set -e

echo "Running integration test..."

# Create minimal test data
./create_minimal_test.sh

# Run the tool on test data
echo ""
echo "Testing basic run..."
cargo run --release -- test_reads.bam test_ref.fa > test_output.csv 2> test_output.log

# Validate output has header and data
LINES=$(wc -l < test_output.csv)
echo "Output has $LINES lines"

if [ $LINES -lt 2 ]; then
    echo "FAIL: Expected at least 2 lines (header + data)"
    exit 1
fi

# Check for expected columns
if ! head -1 test_output.csv | grep -q "Read,Position"; then
    echo "FAIL: Expected 'Read,Position' in header"
    exit 1
fi

echo "Integration test PASSED"

# Show summary from stderr
echo ""
echo "Tool output:"
tail -10 test_output.log

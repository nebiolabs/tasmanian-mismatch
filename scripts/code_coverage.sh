#!/bin/bash
# Code coverage script using grcov
set -e

# Check if grcov is installed
if ! command -v grcov &> /dev/null; then
    echo "Error: grcov not installed. Run: cargo install grcov"
    exit 1
fi

echo "Cleaning previous coverage data..."
rm -rf target/coverage/
mkdir -p target/coverage

# Set environment variables for coverage
export CARGO_INCREMENTAL=0
export RUSTFLAGS="-Cinstrument-coverage"
# Use an absolute path so subprocess invocations (integration tests that run
# the tasmanian-* binaries as child processes) write profraw files to the same
# directory regardless of their working directory.  %p = PID, %m = binary module
# ID — together they guarantee unique filenames across parallel test runs.
export LLVM_PROFILE_FILE="$(pwd)/target/coverage/tasmanian-%p-%m.profraw"
# Run at Debug log level so that debug!() bodies and log_enabled!(Debug) branches
# are exercised.  Without this, those branches are always-false at the default
# Info level and show up as uncovered lines.  Scoped to this crate's library and
# binary targets to avoid noise from dependencies.
export RUST_LOG=rustmanian_mismatch=debug,tasmanian_mismatch=debug,tasmanian_diagnostics=debug,tasmanian_rescale_quality=debug

# Let `cargo test` build the instrumented binaries AND run the tests in one step.
# The integration tests locate the binaries via Cargo's CARGO_BIN_EXE_* env vars,
# so no manual binary path needs to be exported here.  Building separately would
# change the module IDs and prevent grcov from matching subprocess profraw files
# against the binaries it finds on disk.

# Run all tests with instrumented binaries
echo "Building instrumented binary and running tests..."
cargo test

# Generate coverage report
echo "Generating coverage report with grcov..."
grcov target/coverage \
    --binary-path target/debug/ \
    --source-dir . \
    --output-types html,lcov,markdown \
    --branch \
    --ignore-not-existing \
    --ignore "/*" \
    --ignore "target/*" \
    --ignore "tests/*" \
    --excl-line "COV_EXCL_LINE" \
    --excl-start "COV_EXCL_START" \
    --excl-stop "COV_EXCL_STOP" \
    --output-path target/coverage/

echo ""
echo "Coverage Summary:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "$(grep 'Total coverage:' target/coverage/markdown.md)"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "HTML report: target/coverage/html/index.html"
echo "LCOV report: target/coverage/lcov"
echo "Markdown report: target/coverage/markdown.md"

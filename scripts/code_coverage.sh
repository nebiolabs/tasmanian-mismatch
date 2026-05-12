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
# picodup as a child process) write profraw files to the same directory
# regardless of their working directory.  %p = PID, %m = binary module ID —
# together they guarantee unique filenames across parallel test runs.
export LLVM_PROFILE_FILE="$(pwd)/target/coverage/picodup-%p-%m.profraw"
# Run at Debug log level so that debug!() bodies and log_enabled!(Debug) branches
# are exercised.  Without this, those branches are always-false at the default
# Info level and show up as uncovered lines.
export RUST_LOG=picodup=debug

# Let `cargo test` build the instrumented binary AND run the tests in one step.
# Doing a separate `cargo build --bin picodup` beforehand is tempting but wrong:
# `cargo test` will rebuild the binary (different module ID), making grcov unable
# to match the subprocess profraw files against the binary it finds on disk.
export PICODUP_BIN="./target/debug/picodup"

# Run all tests with instrumented binary
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
    --ignore "src/bin/*" \
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

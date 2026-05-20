#!/usr/bin/env python3
"""
Generate a synthetic tasmanian-mismatch TSV for visualization testing.

Produces realistic-looking mismatch count data with all 12 mismatch types,
both read numbers, both reference orders, across 150 fragment positions.

Usage:
    python scripts/generate_test_data.py > test_viz_data.tsv
    python scripts/generate_test_data.py -o test_viz_data.tsv
"""

import argparse
import math
import random
import sys

MISMATCHES = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G",
]

# Bases that appear as the reference base in each mismatch
REF_BASE = {m: m[0] for m in MISMATCHES}

# Give C>T and G>A elevated counts (common deamination artifacts)
MISMATCH_WEIGHTS = {
    "A>C": 0.3, "A>G": 0.5, "A>T": 0.2,
    "C>A": 0.4, "C>G": 0.2, "C>T": 2.5,
    "G>A": 2.0, "G>C": 0.2, "G>T": 0.3,
    "T>A": 0.2, "T>C": 0.4, "T>G": 0.3,
}

POSITIONS = 100


def gaussian(pos, center, sigma, amplitude):
    return amplitude * math.exp(-((pos - center) ** 2) / (2 * sigma ** 2))


def count_for(mismatch, read_num, ref_order, pos):
    """Synthetic count: end-bias bumps + flat background + noise."""
    weight = MISMATCH_WEIGHTS[mismatch]

    # Fragment-position model:
    # ref_order=1 → mismatch enriched near fragment start (low positions)
    # ref_order=2 → mismatch enriched near fragment end (high positions)
    if ref_order == 1:
        end_signal = gaussian(pos, 5, 8, 80 * weight)
    else:
        end_signal = gaussian(pos, POSITIONS - 5, 8, 80 * weight)

    background = weight * 15
    noise = random.gauss(0, weight * 3)

    raw = end_signal + background + noise
    return max(0, int(round(raw)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    random.seed(args.seed)

    out = open(args.output, "w") if args.output else sys.stdout

    header = "base_change\tread_num\treference_order\tfragment_position\tcount"
    print(header, file=out)

    for mismatch in MISMATCHES:
        for read_num in [1, 2]:
            for ref_order in [1, 2]:
                for pos in range(1, POSITIONS + 1):
                    count = count_for(mismatch, read_num, ref_order, pos)
                    if count > 0:
                        print(f"{mismatch}\t{read_num}\t{ref_order}\t{pos}\t{count}", file=out)

    if args.output:
        out.close()
        print(f"Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

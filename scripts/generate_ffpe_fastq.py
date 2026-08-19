#!/usr/bin/env python3
"""
Generate synthetic paired-end FASTQ files with FFPE artifacts from a reference genome.

FFPE artifacts manifest as deamination damage:
  - G>A substitutions accumulating at the END of read1
  - C>T substitutions accumulating at the BEGINNING of read2

This simulates the characteristic pattern of formalin-induced cytosine deamination
at the ends of DNA fragments in FFPE-processed samples.

Usage:
    python generate_ffpe_fastq.py \
        --reference /bioinfo/ref/t2t_chm13_v2/fasta/t2t_chm13_v2.fa \
        --read-length 100 \
        --min-insert 300 \
        --max-insert 700 \
        --num-pairs 1000000 \
        --output-prefix ffpe_simulated \
        --damage-rate 0.3 \
        --damage-decay 0.12
"""

import argparse
import random
import sys
import math
from pathlib import Path
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate paired-end FASTQ files with FFPE artifacts"
    )
    parser.add_argument(
        "--reference", "-r", required=True,
        help="Path to reference FASTA file"
    )
    parser.add_argument(
        "--read-length", "-l", type=int, default=100,
        help="Read length in bp (default: 100)"
    )
    parser.add_argument(
        "--min-insert", type=int, default=300,
        help="Minimum insert size in bp (default: 300)"
    )
    parser.add_argument(
        "--max-insert", type=int, default=700,
        help="Maximum insert size in bp (default: 700)"
    )
    parser.add_argument(
        "--num-pairs", "-n", type=int, default=1_000_000,
        help="Number of read pairs to generate (default: 1000000)"
    )
    parser.add_argument(
        "--output-prefix", "-o", default="ffpe_simulated",
        help="Output file prefix (default: ffpe_simulated)"
    )
    parser.add_argument(
        "--damage-rate", "-d", type=float, default=0.3,
        help="Maximum damage probability at fragment ends (default: 0.3)"
    )
    parser.add_argument(
        "--damage-decay", type=float, default=0.12,
        help="Exponential decay rate for damage probability moving away from "
             "fragment ends (default: 0.12). Higher = faster decay."
    )
    parser.add_argument(
        "--base-error-rate", type=float, default=0.001,
        help="Background sequencing error rate (default: 0.001)"
    )
    parser.add_argument(
        "--seed", "-s", type=int, default=None,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--chromosomes", "-c", nargs="*", default=None,
        help="Restrict to specific chromosomes (default: all)"
    )
    return parser.parse_args()


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def load_reference(fasta_path, chromosomes=None):
    """Load reference sequences into memory, optionally filtering by chromosome."""
    print(f"Loading reference from {fasta_path}...", file=sys.stderr)
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        if chromosomes and record.id not in chromosomes:
            continue
        seq = str(record.seq).upper()
        # Only keep sequences longer than max possible insert
        if len(seq) > 1000:
            sequences[record.id] = seq
            print(f"  Loaded {record.id}: {len(seq):,} bp", file=sys.stderr)
    if not sequences:
        print("ERROR: No sequences loaded from reference!", file=sys.stderr)
        sys.exit(1)
    print(f"  Total: {len(sequences)} sequences loaded", file=sys.stderr)
    return sequences


def get_damage_probability(position_from_end, damage_rate, decay):
    """
    Calculate damage probability as exponential decay from fragment end.

    position_from_end: 0-based distance from the fragment end
    damage_rate: maximum damage probability (at position 0)
    decay: exponential decay constant
    """
    return damage_rate * math.exp(-decay * position_from_end)


def apply_ffpe_damage_read1(seq, damage_rate, decay):
    """
    Apply FFPE damage to read1: G>A at the END of the read.

    In FFPE, deamination occurs at fragment ends. For read1 (forward strand),
    this manifests as G>A substitutions accumulating toward the 3' end (end of read).
    """
    seq_list = list(seq)
    read_len = len(seq_list)

    for i in range(read_len):
        # Distance from the end of the read
        dist_from_end = (read_len - 1) - i
        prob = get_damage_probability(dist_from_end, damage_rate, decay)

        if seq_list[i] == 'G' and random.random() < prob:
            seq_list[i] = 'A'

    return ''.join(seq_list)


def apply_ffpe_damage_read2(seq, damage_rate, decay):
    """
    Apply FFPE damage to read2: C>T at the BEGINNING of the read.

    In FFPE, deamination occurs at fragment ends. For read2 (reverse strand),
    this manifests as C>T substitutions accumulating toward the 5' end (beginning of read).
    """
    seq_list = list(seq)

    for i in range(len(seq_list)):
        # Distance from the beginning of the read
        dist_from_start = i
        prob = get_damage_probability(dist_from_start, damage_rate, decay)

        if seq_list[i] == 'C' and random.random() < prob:
            seq_list[i] = 'T'

    return ''.join(seq_list)


def apply_sequencing_errors(seq, error_rate):
    """Apply random sequencing errors at a low background rate."""
    if error_rate <= 0:
        return seq
    bases = list("ACGT")
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if seq_list[i] in bases and random.random() < error_rate:
            alternatives = [b for b in bases if b != seq_list[i]]
            seq_list[i] = random.choice(alternatives)
    return ''.join(seq_list)


def generate_quality_string(read_length, is_damaged_end=False):
    """
    Generate a realistic quality string.

    Base quality is high (~30-37) with slight degradation at read ends.
    Damaged positions get slightly lower quality.
    """
    quals = []
    for i in range(read_length):
        # Base quality around Q30-Q37
        base_qual = random.gauss(35, 2)
        # Slight quality drop at read ends
        dist_from_center = abs(i - read_length / 2) / (read_length / 2)
        base_qual -= dist_from_center * 3
        # Clamp to valid Phred range
        qual = max(2, min(40, int(base_qual)))
        quals.append(chr(qual + 33))
    return ''.join(quals)


def sample_insert_size(min_insert, max_insert):
    """
    Sample insert size from a normal distribution truncated to [min_insert, max_insert].

    The mean is the midpoint and sd is chosen so that ~95% falls within range.
    """
    mean = (min_insert + max_insert) / 2
    sd = (max_insert - min_insert) / 4  # ~95% within range
    while True:
        size = int(random.gauss(mean, sd))
        if min_insert <= size <= max_insert:
            return size


def generate_reads(args):
    """Main read generation loop."""
    sequences = load_reference(args.reference, args.chromosomes)

    # Build a weighted list of chromosomes by length for random sampling
    chrom_names = list(sequences.keys())
    chrom_lengths = [len(sequences[c]) for c in chrom_names]
    total_length = sum(chrom_lengths)
    chrom_weights = [l / total_length for l in chrom_lengths]

    r1_path = Path(f"{args.output_prefix}_R1.fastq")
    r2_path = Path(f"{args.output_prefix}_R2.fastq")

    print(f"Generating {args.num_pairs:,} read pairs...", file=sys.stderr)
    print(f"  Read length: {args.read_length} bp", file=sys.stderr)
    print(f"  Insert size: {args.min_insert}-{args.max_insert} bp", file=sys.stderr)
    print(f"  Damage rate: {args.damage_rate} (decay: {args.damage_decay})", file=sys.stderr)
    print(f"  Output: {r1_path}, {r2_path}", file=sys.stderr)

    generated = 0
    skipped = 0

    with open(r1_path, 'w') as f1, open(r2_path, 'w') as f2:
        while generated < args.num_pairs:
            # Select a random chromosome weighted by length
            chrom = random.choices(chrom_names, weights=chrom_weights, k=1)[0]
            chrom_seq = sequences[chrom]
            chrom_len = len(chrom_seq)

            # Sample insert size
            insert_size = sample_insert_size(args.min_insert, args.max_insert)

            # Sample a random position (ensuring fragment fits)
            max_start = chrom_len - insert_size
            if max_start <= 0:
                continue
            start = random.randint(0, max_start)

            # Extract fragment
            fragment = chrom_seq[start:start + insert_size]

            # Skip fragments with Ns
            if 'N' in fragment:
                skipped += 1
                continue

            # Read1: first read_length bases of fragment (forward strand)
            read1_seq = fragment[:args.read_length]
            # Read2: last read_length bases of fragment, reverse complemented
            read2_seq = reverse_complement(fragment[-args.read_length:])

            if len(read1_seq) < args.read_length or len(read2_seq) < args.read_length:
                skipped += 1
                continue

            # Apply FFPE damage
            read1_seq = apply_ffpe_damage_read1(read1_seq, args.damage_rate, args.damage_decay)
            read2_seq = apply_ffpe_damage_read2(read2_seq, args.damage_rate, args.damage_decay)

            # Apply background sequencing errors
            read1_seq = apply_sequencing_errors(read1_seq, args.base_error_rate)
            read2_seq = apply_sequencing_errors(read2_seq, args.base_error_rate)

            # Generate quality strings
            qual1 = generate_quality_string(args.read_length)
            qual2 = generate_quality_string(args.read_length)

            # Write FASTQ entries
            generated += 1
            read_name = f"simulated_{generated}_{chrom}_{start}_{insert_size}"

            f1.write(f"@{read_name}/1\n{read1_seq}\n+\n{qual1}\n")
            f2.write(f"@{read_name}/2\n{read2_seq}\n+\n{qual2}\n")

            if generated % 100_000 == 0:
                print(f"  Generated {generated:,} / {args.num_pairs:,} pairs...",
                      file=sys.stderr)

    print(f"\nDone! Generated {generated:,} read pairs "
          f"(skipped {skipped:,} fragments with Ns)", file=sys.stderr)
    print(f"Output files:\n  {r1_path}\n  {r2_path}", file=sys.stderr)


def main():
    args = parse_args()
    if args.seed is not None:
        random.seed(args.seed)
    generate_reads(args)


if __name__ == "__main__":
    main()

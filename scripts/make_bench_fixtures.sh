#!/usr/bin/env bash
# Rebuild the benchmark fixture: 25% of the read pairs in chr20:57,000,001-58,000,000 of the
# EM-seq paper's 10 ng FFPE library, moved onto a 1 Mb contig so the committed reference is
# 1 MB instead of all of chr20. Requires samtools 1.20+.
#
# The input is the BAM submitted to SRA for run SRR13076821 (10ng_EM-seq_FFPE_1, study
# SRP233253): bwa-meth against the reference below, duplicate-marked and coordinate-sorted.
# Download the original file, 10ng_EM_seq_FFPE_1.md.bam, from the run's Data access tab:
#
#   https://www.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR13076821&display=alignment
#
# Then index it and run this script from the repository root:
#
#   samtools index 10ng_EM_seq_FFPE_1.md.bam
#   scripts/make_bench_fixtures.sh 10ng_EM_seq_FFPE_1.md.bam
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: $0 <SRR13076821 BAM>; see the header of this script for where to get it" >&2
    exit 1
fi
BAM=$1
if [ ! -e "$BAM.bai" ] && [ ! -e "${BAM%.bam}.bai" ]; then
    echo "$BAM has no index; run: samtools index $BAM" >&2
    exit 1
fi
REF_URL=https://neb-em-seq-sra.s3.us-east-1.amazonaws.com/grch38_core%2Bmeth_controls.fa.gz
OUT=${OUT:-tests/fixtures}
START=57000001
END=58000000

# The reference is plain gzip with no index, so stream it and keep only chr20. awk exits at
# the next contig, which ends the download there. curl then reports a write failure (23),
# so success is judged by chr20's length instead of the pipeline's status.
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
curl -fsS "$REF_URL" | gzip -dc |
    awk '/^>/ { if (keep) exit; keep = ($1 == ">chr20") } keep' > "$tmp/chr20.fa" || true
if ! samtools faidx "$tmp/chr20.fa" || [ "$(cut -f2 "$tmp/chr20.fa.fai")" != 64444167 ]; then
    echo "failed to fetch all of chr20 from $REF_URL" >&2
    exit 1
fi
samtools faidx "$tmp/chr20.fa" "chr20:$START-$END" | sed '1s/.*/>chr20_57M/' > "$OUT/chr20_57M.fa"

# Primary, QC-passed pairs with both mates starting in the window, in window coordinates.
{
    printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr20_57M\tLN:%s\n' $((END - START + 1))
    samtools view -H "$BAM" | grep '^@RG'
    samtools view -F 0xB00 "$BAM" "chr20:$START-$END" |
        awk -F'\t' -v OFS='\t' -v s="$START" -v e=$((END - 500)) '
            $7 == "=" && $4 >= s && $8 >= s && $4 <= e && $8 <= e {
                $3 = "chr20_57M"; $4 -= s - 1; $8 -= s - 1; print
            }'
} | samtools view -b --no-PG --subsample 0.25 --subsample-seed 42 -o "$OUT/emseq_ffpe_chr20_57M.bam" -
samtools index "$OUT/emseq_ffpe_chr20_57M.bam"

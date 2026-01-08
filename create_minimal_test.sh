#!/bin/bash
# Create a minimal test BAM file with known patterns

echo "Creating minimal test files for tasmanian-mismatch..."

# Create a tiny reference (100bp)
cat > test_ref.fa << 'EOF'
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF

# Create a SAM file with a few reads (some with mismatches)
cat > test_reads.sam << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:100
read1	99	chr1	1	60	20M	=	30	49	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	NM:i:0
read1	147	chr1	30	60	20M	=	1	-49	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	NM:i:0
read2	99	chr1	10	60	20M	=	40	50	ACGTACGTACGTACGTAGGT	IIIIIIIIIIIIIIIIIIII	NM:i:1
read2	147	chr1	40	60	20M	=	10	-50	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	NM:i:0
read3	0	chr1	50	60	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	NM:i:0
EOF

# Convert SAM to BAM
if command -v samtools &> /dev/null; then
    echo "Converting SAM to BAM..."
    samtools view -Sb test_reads.sam | samtools sort -o test_reads.bam
    samtools index test_reads.bam
    echo "Created test_reads.bam and test_reads.bam.bai"
else
    echo "WARNING: samtools not found. Please install samtools and run:"
    echo "  samtools view -Sb test_reads.sam > test_reads.bam"
    echo "  samtools index test_reads.bam"
fi

echo "Created test_ref.fa and test_reads.sam"
echo ""

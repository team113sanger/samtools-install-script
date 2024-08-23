#!/bin/bash
set -uo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Check if version argument is provided
if [ $# -eq 0 ]; then
    echo -e "${RED}Error: Expected samtools version must be provided as an argument.${NC}"
    echo "Usage: $0 <version>"
    exit 1
fi

# Expected samtools version
EXPECTED_VERSION="$1"

# Function to print test results
print_result() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}PASS${NC}: $1"
    else
        echo -e "${RED}FAIL${NC}: $1"
        exit 1
    fi
}

# Function to compare version numbers
# e.g. version_ge "1.10" "1.9" -- exit status 0
# e.g. version_ge "1.10" "1.11" -- exit status 1
version_ge() {
    test "$(echo "$@" | tr " " "\n" | sort -rV | head -n 1)" == "$1"
}


echo "Starting samtools test suite..."
echo "Expected samtools version: $EXPECTED_VERSION"

# Test 1: Check if samtools is in PATH
echo "Test 1: Checking if samtools is in PATH..."
which samtools > /dev/null
print_result "samtools is in PATH"

# Test 2: Check samtools version
echo "Test 2: Checking samtools version..."
version_output=$(samtools --version | head -n 1)
if [[ "$version_output" == *"$EXPECTED_VERSION"* ]]; then
    print_result "samtools version matches expected version ($EXPECTED_VERSION)"
else
    echo -e "${RED}FAIL${NC}: Version mismatch. Expected $EXPECTED_VERSION, got $version_output"
    exit 1
fi

# Test 3: Check for required programs in PATH
echo "Test 3: Checking for required programs in PATH..."
# For 1.10 only a different list of programs is required. 1.10 is the oldest
# version we support, everything newer should be fine.
case "$EXPECTED_VERSION" in
    "1.10")
        required_programs=(
            "ace2sam"
            "maq2sam-long"
            "maq2sam-short"
            "md5fa"
            "md5sum-lite"
            "plot-bamstats"
            "samtools"
            "wgsim"
        )
        ;;
    *)
        required_programs=(
            "ace2sam"
            "maq2sam-long"
            "maq2sam-short"
            "md5fa"
            "md5sum-lite"
            "plot-ampliconstats"
            "plot-bamstats"
            "samtools"
            "wgsim"
        )
        ;;
esac
for program in "${required_programs[@]}"; do
    which "$program" > /dev/null
    print_result "$program is in PATH"
done

# Test 4: Index the reference FASTA
echo "Test 4: Indexing the reference FASTA..."
samtools faidx ex1.fa
[ -f ex1.fa.fai ] && print_result "Indexed reference FASTA"

# Test 5: Convert SAM to BAM
echo "Test 5: Converting SAM to BAM..."
samtools view -S -b -t ex1.fa.fai -o ex1.bam ex1.sam.gz
[ -f ex1.bam ] && print_result "Converted SAM to BAM"

# Test 6: Index BAM file
echo "Test 6: Indexing BAM file..."
samtools index ex1.bam
[ -f ex1.bam.bai ] && print_result "Indexed BAM file"

# Test 7: View a portion of the BAM file
echo "Test 7: Viewing a portion of the BAM file..."
samtools view ex1.bam seq2:450-550 > view_output.txt
[ -s view_output.txt ] && print_result "Viewed portion of BAM file"

# Test 8: Generate pileup
echo "Test 8: Generating pileup..."
samtools mpileup -f ex1.fa ex1.bam > ex1.pileup
[ -s ex1.pileup ] && print_result "Generated pileup"

# Test 9: Generate uncompressed VCF
if version_ge "$EXPECTED_VERSION" "1.15"; then
    echo "Test 9: Skipped for version $EXPECTED_VERSION (>= 1.15)"
else
    echo "Test 9: Generating uncompressed VCF..."
    samtools mpileup -vu -f ex1.fa ex1.bam > ex1.vcf
    [ -s ex1.vcf ] && print_result "Generated uncompressed VCF"
fi

# Test 10: Generate compressed BCF
if version_ge "$EXPECTED_VERSION" "1.15"; then
    echo "Test 10: Skipped for version $EXPECTED_VERSION (>= 1.15)"
else
    echo "Test 10: Generating compressed BCF..."
    samtools mpileup -g -f ex1.fa ex1.bam > ex1.bcf
    [ -s ex1.bcf ] && print_result "Generated compressed BCF"
fi

# Test 11: Convert BAM to SAM
echo "Test 11: Converting BAM back to SAM..."
samtools view -h -o ex1_new.sam ex1.bam
[ -f ex1_new.sam ] && print_result "Converted BAM back to SAM"

# Test 12: Compare original and new SAM files (ignoring headers)
echo "Test 12: Comparing original and new SAM files (ignoring headers)..."
gunzip -c ex1.sam.gz | sort > ex1_sorted.sam  # Compressed SAM file has no headers
samtools view ex1_new.sam | sort > ex1_new_sorted.sam
diff ex1_sorted.sam ex1_new_sorted.sam
print_result "Original and new SAM files are identical (ignoring headers)"


# Test 13: Test tview (this just checks if it runs without error)
echo "Test 13: Testing tview..."
echo "q" | samtools tview -p seq2:450 ex1.bam ex1.fa > /dev/null 2>&1
print_result "tview ran successfully"

echo -e "\n${GREEN}All tests passed successfully!${NC}"

# Cleanup
rm -f ex1.fa.fai ex1.bam ex1.bam.bai view_output.txt ex1.pileup ex1.vcf ex1.bcf ex1_new.sam ex1_original.sam ex1_sorted.sam ex1_new_sorted.sam

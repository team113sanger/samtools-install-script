#!/bin/bash
set -uo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# Check if version argument is provided
if [ $# -eq 0 ]; then
    echo -e "${RED}Error: Expected htslib version must be provided as an argument.${NC}"
    echo "Usage: $0 <version>"
    exit 1
fi

# Check if `EXPECTS_HTSLIB` is set (allow any value) - if absent, exit early
if [ -z ${EXPECTS_HTSLIB+x} ]; then
    echo -e "${YELLOW}Warning: Skipping htslib tests as EXPECTS_HTSLIB is not set. Set to any value to run the tests.${NC}"
    exit 0
fi

# Expected htslib version
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

echo "Starting htslib test suite..."
echo "Expected htslib version: $EXPECTED_VERSION"

# Test 1: Check if pkg-config can find htslib
echo "Test 1: Checking if pkg-config can find htslib..."
pkg-config --exists htslib
print_result "pkg-config can find htslib"

# Test 2: Compile test program
echo "Test 2: Compiling test program..."
gcc -o test_htslib test_htslib.c $(pkg-config --cflags --libs htslib)
print_result "Test program compiled successfully"

# Test 3: Run test program and check version
echo "Test 3: Running test program and checking version..."
output=$(./test_htslib)
echo "$output"
if [[ "$output" == "htslib version: $EXPECTED_VERSION" ]]; then
    print_result "htslib version matches expected version ($EXPECTED_VERSION)"
else
    echo -e "${RED}FAIL${NC}: Version mismatch. Expected $EXPECTED_VERSION, got ${output#htslib version: }"
    exit 1
fi

# Test 4: Check for htslib executables
echo "Test 4: Checking for htslib executables..."
required_programs=(
    "bgzip"
    "tabix"
    "htsfile"
)
for program in "${required_programs[@]}"; do
    which "$program" > /dev/null
    print_result "$program is in PATH"
done

# Test 5: Basic functionality test of bgzip
echo "Test 5: Testing basic functionality of bgzip..."
echo "Hello, htslib!" > test.txt
bgzip test.txt
[ -f test.txt.gz ] && print_result "bgzip compressed the file"
bgzip -d test.txt.gz
[ -f test.txt ] && [ "$(cat test.txt)" == "Hello, htslib!" ] && print_result "bgzip decompressed the file correctly"

# Test 6: Basic functionality test of tabix
echo "Test 6: Testing basic functionality of tabix..."
echo -e "chr1\t100\t200\tfeature1\nchr1\t150\t250\tfeature2" > test.bed
bgzip test.bed
tabix -p bed test.bed.gz
[ -f test.bed.gz.tbi ] && print_result "tabix indexed the bgzipped BED file"

# Test 7: Basic functionality test of htsfile
echo "Test 7: Testing basic functionality of htsfile..."
htsfile test.bed.gz > htsfile_output.txt
grep -qP "test.bed.gz:\tBED BGZF-compressed genomic region data" htsfile_output.txt
print_result "htsfile correctly identified the bgzipped BED file"

echo -e "\n${GREEN}All tests passed successfully!${NC}"

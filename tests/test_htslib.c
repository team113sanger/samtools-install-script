#include <stdio.h>
#include <htslib/hts.h>

int main() {
    printf("htslib version: %s\n", hts_version());
    return 0;
}

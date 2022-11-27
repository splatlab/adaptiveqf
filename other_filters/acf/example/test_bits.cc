#include "cuckoofilter.h"
#include "bitsutil.h"

#include <assert.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <vector>

using cuckoofilter::CuckooFilter;

void ulltobin(char buffer[65], uint64_t n) {
	for (int i = 63; i >= 0; i--) {
		buffer[i] = n & 1 ? '1' : '0';
		n >>= 1;
	}
	buffer[64] = '\0';
}

double rand_zipfian(double s, double max) {
        double p = (double)rand() / RAND_MAX;

        double pD = p * (12 * (pow(max, -s + 1) - 1) / (1 - s) + 6 + 6 * pow(max, -s) + s - s * pow(max, -s + 1));
        double x = max / 2;
        while (true) {
                double m = pow(x, -s - 2);
                double mx = m * x;
                double mxx = mx * x;
                double mxxx = mxx * x;

                double b = 12 * (mxxx - 1) / (1 - s) + 6 + 6 * mxx + s - (s * mx) - pD;
                double c = 12 * mxx - (6 * s * mx) + (m * s * (s + 1));
                double newx = x - b / c > 1 ? x - b / c : 1;
                if (abs(newx - x) <= 0.01) { // this is the tolerance for approximation
                        return newx;
                }
                x = newx;
        }
}

int main(int argc, char **argv) {
	assert(argc >= 3);
	char buffer[65];
	printf("%llx\t%llx\n", strtoull(argv[1], NULL, 2), strtoull(argv[2], NULL, 2));
	ulltobin(buffer, hasvalue16(strtoull(argv[1], NULL, 2), strtoull(argv[2], NULL, 2)));
	printf("%s\n", buffer);
	ulltobin(buffer, hasvalue4(0x0000000000000000ULL, 1));
	printf("%s\n", buffer);

	return 0;
}

#include <stdint.h>
#include <time.h>
#include <x86intrin.h>

#include "Color.h"
#include "Util.h"

unsigned long long time_diff(double (*func)(Color, Color), struct Color col1,
			     struct Color col2, const long int RUNS) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(col1, col2);
	}

	unsigned long long total = 0;
	for (long int i = 0; i < RUNS; ++i) {

		unsigned long long tic = __rdtsc();
		func(col1, col2);
		unsigned long long toc = __rdtsc();
		total += toc - tic;
	}

	return total / RUNS;
}

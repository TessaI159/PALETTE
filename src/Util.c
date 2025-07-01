#include <x86intrin.h>

#include "Color.h"
#include "Util.h"

uint64_t time_diff(float (*func)(Color*, Color*), struct Color col1,
			     struct Color col2, const uint64_t RUNS) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(&col1, &col2);
	}

	unsigned long long total = 0;
	for (uint64_t i = 0; i < RUNS; ++i) {

		unsigned long long tic = __rdtsc();
		func(&col1, &col2);
		unsigned long long toc = __rdtsc();
		total += toc - tic;
	}

	return total / RUNS;
}

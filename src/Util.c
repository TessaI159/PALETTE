#if defined(__x86_64__)
#include <x86intrin.h>
#endif
#include "Color.h"
#include "Util.h"

uint64_t time_diff(float (*func)(Color *, Color *), struct Color col1,
		   struct Color col2, const uint64_t RUNS) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(&col1, &col2);
	}

	unsigned long long total = 0;
	unsigned long long tic	 = __rdtsc();
#pragma omp parallel for reduction(+ : total)

	for (uint64_t i = 0; i < RUNS; ++i) {

		func(&col1, &col2);
	}
	unsigned long long toc = __rdtsc();

	return (toc - tic) / RUNS;
}

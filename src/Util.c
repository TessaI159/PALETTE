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

	uint64_t tic = __rdtsc();

	for (uint64_t i = 0; i < RUNS; ++i) {

		func(&col1, &col2);
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

uint64_t time_conv(void (*func)(Color *, enum ColorSpace), Color col,
		   enum ColorSpace s, const uint64_t RUNS) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(&col, s);
	}
	uint64_t tic = __rdtsc();

	for (uint64_t i = 0; i < RUNS; ++i) {

		func(&col, s);
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

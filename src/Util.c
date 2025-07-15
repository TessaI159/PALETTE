#include <stdio.h>
#include <string.h>
#if defined(__x86_64__)
#include <x86intrin.h>
#endif
#include "Align.h"
#include "Color.h"
#include "Util.h"

uint64_t time_diff(float (*func)(Color *, Color *), struct Color *col1,
		   struct Color *col2, const uint64_t RUNS) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(col1, col2);
	}

	uint64_t tic = __rdtsc();
	for (uint64_t i = 0; i < RUNS; ++i) {

		func(col1, col2);
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

uint64_t time_conv(void (*func)(Color *, enum ColorSpace), Color *col,
		   enum ColorSpace f, enum ColorSpace t, const uint64_t RUNS) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(col, t);
	}

	uint64_t tic = __rdtsc();
	for (uint64_t i = 0; i < RUNS; ++i) {
		func(col, t);
		col->current_space = f;
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

uint64_t time_avg_cie(Color (*func)(const cielab_SoA *, uint16_t),
		  cielab_SoA *colors, uint16_t num_col, uint64_t RUNS) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(colors, num_col);
	}

	uint64_t tic = __rdtsc();
	for (uint64_t i = 0; i < RUNS; ++i) {
		func(colors, num_col);
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

uint64_t time_avg_ok(Color (*func)(const oklab_SoA *, uint16_t),
		  oklab_SoA *colors, uint16_t num_col, uint64_t RUNS) {
	/* Warm the function up */
	for (int i = 0; i < 1000; ++i) {
		func(colors, num_col);
	}

	uint64_t tic = __rdtsc();
	for (uint64_t i = 0; i < RUNS; ++i) {
		func(colors, num_col);
	}
	uint64_t toc = __rdtsc();

	return (toc - tic) / RUNS;
}

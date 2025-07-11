#include <stdint.h>
#include <stdio.h>
#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "Color.h"

struct Color lab_avg_avx2(const struct Color **colors, uint16_t num_col) {
	puts("AVX2 average\n");
	return Color_create(0, 0, 0);
}

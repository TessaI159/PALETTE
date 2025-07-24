#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "Average.h"
#include "Average_Internal.h"
#include "Color.h"
#include "FeatureDetection.h"
#include "Parameters.h"

#define INST_SETS 3
#define SPACES 2
#define BOOLS 2

static void (*color_avg_intern)(const struct Color *restrict, struct Color *restrict, uint8_t);

static void (*avg_table[INST_SETS][SPACES][BOOLS])(const struct Color *restrict,
						   struct Color *restrict, uint8_t) = {
    [AVX] = {[OK]  = {[true] = oklab_avg_avx_cw, [false] = oklab_avg_avx},
	     [CIE] = {[true] = cielab_avg_avx_cw, [false] = cielab_avg_avx}},
    [SSE] = {[OK]  = {[true] = oklab_avg_sse_cw, [false] = oklab_avg_sse},
	     [CIE] = {[true] = cielab_avg_sse_cw, [false] = cielab_avg_sse}},
    [FB]  = {[OK]  = {[true] = oklab_avg_fb_cw, [false] = oklab_avg_fb},
	     [CIE] = {[true] = cielab_avg_fb_cw, [false] = cielab_avg_fb}}};

static inline void average_init() {
	enum InstSet inst;

	if (!features.initialized) {
		query_features(&features);
	}

	if (features.avx) {
		inst = AVX;
	} else if (features.sse) {
		inst = SSE;
	} else {
		inst = FB;
	}
	color_avg_intern = avg_table[inst][g_params.space][g_params.cw];
}

static inline void check_initialized() {
	static bool average_initialized = false;
	if (!average_initialized) {
		average_init();
		average_initialized = true;
	}
}

void color_avg(const struct Color *restrict colors, struct Color *restrict cents,
	       uint8_t which_cent) {
	check_initialized();
	return color_avg_intern(colors, cents, which_cent);
}

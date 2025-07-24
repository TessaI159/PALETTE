#include <stdlib.h>

#include "Color.h"
#include "Difference.h"
#include "Difference_Internal.h"
#include "FeatureDetection.h"
#include "Parameters.h"

#define INST_SETS 3
#define SPACES 2
#define DIFF_FUNCS 3

static void (*delta_diff_intern)(const struct Color *restrict,
				 const struct Color *restrict, float *restrict);

static void (*diff_table[INST_SETS][SPACES][DIFF_FUNCS])(
    const struct Color *restrict, const struct Color *restrict,
    float *restrict) = {
    [AVX] = {[OK]  = {[EUCLIDEAN] = delta_lab_euc_diff_avx,
		      [D94]	  = NULL,
		      [D2000]	  = NULL},
	     [CIE] = {[EUCLIDEAN] = delta_lab_euc_diff_avx,
		      [D94]	  = delta_cie94_diff_avx,
		      [D2000]	  = delta_ciede2000_diff_avx}},
    [SSE] = {[OK]  = {[EUCLIDEAN] = delta_lab_euc_diff_sse,
		      [D94]	  = NULL,
		      [D2000]	  = NULL},
	     [CIE] = {[EUCLIDEAN] = delta_lab_euc_diff_sse,
		      [D94]	  = delta_cie94_diff_sse,
		      [D2000]	  = delta_ciede2000_diff_sse}},

    [FB] = {[OK]  = {[EUCLIDEAN] = delta_lab_euc_diff_fallback,
		     [D94]	 = NULL,
		     [D2000]	 = NULL},
	    [CIE] = {[EUCLIDEAN] = delta_lab_euc_diff_fallback,
		     [D94]	 = delta_cie94_diff_fallback,
		     [D2000]	 = delta_ciede2000_diff_fallback}}};

static inline void difference_init() {
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
	delta_diff_intern = diff_table[inst][g_params.space][g_params.diff_fn];
}

static inline void check_initialized() {
	static bool diff_initd = false;
	if (!diff_initd) {
		difference_init();
		diff_initd = true;
	}
}

void delta_diff(const struct Color *restrict colors,
		const struct Color *restrict cents, float *restrict diffs) {
	check_initialized();
	return delta_diff_intern(colors, cents, diffs);
}

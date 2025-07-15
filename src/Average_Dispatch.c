#include <stdio.h>
#include <stdlib.h>

#include "Average.h"
#include "Average_Internal.h"
#include "Color.h"
#include "FeatureDetection.h"
#include "Parameters.h"

static struct Color (*cielab_avg_intern)(
    const struct cielab_SoA *__restrict colors,
    uint16_t num_col) = cielab_avg_fallback;

static struct Color (*oklab_avg_intern)(
    const struct oklab_SoA *__restrict colors,
    uint16_t num_col) = oklab_avg_fallback;

static inline void average_init() {
	if (!features.initialized) {
		query_features(&features);
	}
	if (features.avx) {
		if (global_parameters.chroma_weight) {
			cielab_avg_intern = cielab_avg_avx_cw;
			oklab_avg_intern  = oklab_avg_avx_cw;
		} else {
			cielab_avg_intern = cielab_avg_avx;
			oklab_avg_intern  = oklab_avg_avx;
		}

	} else if (features.sse) {
		if (global_parameters.chroma_weight) {
			cielab_avg_intern = cielab_avg_sse_cw;
			oklab_avg_intern  = oklab_avg_sse_cw;
		} else {
			cielab_avg_intern = cielab_avg_sse;
			oklab_avg_intern  = oklab_avg_sse;
		}
	} else {
		if (global_parameters.chroma_weight) {
			cielab_avg_intern = cielab_avg_fallback_cw;
			oklab_avg_intern  = oklab_avg_fallback_cw;
		} else {
			cielab_avg_intern = cielab_avg_fallback;
			oklab_avg_intern  = oklab_avg_fallback;
		}
	}
}

static inline void check_initialized() {
	static bool average_initialized = false;
	if (!average_initialized) {
		average_init();
		average_initialized = true;
	}
}

struct Color color_avg(const void *__restrict colors, uint16_t num_col) {
	check_initialized();
	switch (global_parameters.working_colorspace) {
	case WORKING_CIELAB:
		return cielab_avg_intern((const struct cielab_SoA *)colors,
					 num_col);
		break;
	case WORKING_OKLAB:
		return oklab_avg_intern((const struct oklab_SoA *)colors,
					num_col);
		break;
	default:
		fprintf(stderr, "Unknown working color space");
		exit(1);
	}
}

#include "Average.h"
#include "Color.h"
#include "FeatureDetection.h"
#include "Flags.h"

static struct Color (*cielab_avg_intern)(
    const struct cielab_SoA *__restrict colors,
    uint16_t num_col) = cielab_avg_fallback;

static struct Color (*oklab_avg_intern)(
    const struct oklab_SoA *__restrict colors,
    uint16_t num_col) = oklab_avg_fallback;

static bool check_initialized() {
	static bool average_initialized = false;
	if (!average_initialized) {
		average_init();
		average_initialized = true;
	}
	return average_initialized;
}

struct Color cielab_avg(const struct cielab_SoA *__restrict colors,
			uint16_t num_col) {
	check_initialized();
	return cielab_avg_intern(colors, num_col);
}

struct Color oklab_avg(const struct oklab_SoA *__restrict colors,
		       uint16_t num_col) {
	check_initialized();
	return oklab_avg_intern(colors, num_col);
}

void average_init() {
	if (!features.initialized) {
		query_features(&features);
	}
	if (features.avx) {
		if (global_flags.chroma_average) {
			cielab_avg_intern = cielab_avg_avx_cw;
			oklab_avg_intern  = oklab_avg_avx_cw;
		} else {
			cielab_avg_intern = cielab_avg_avx;
			oklab_avg_intern  = oklab_avg_avx;
		}

	} else if (features.sse) {
		if (global_flags.chroma_average) {
			cielab_avg_intern = cielab_avg_sse_cw;
			oklab_avg_intern  = oklab_avg_sse_cw;
		} else {
			cielab_avg_intern = cielab_avg_sse;
			oklab_avg_intern  = oklab_avg_sse;
		}
	} else {
		if (global_flags.chroma_average) {
			cielab_avg_intern = cielab_avg_fallback_cw;
			oklab_avg_intern  = oklab_avg_fallback_cw;
		} else {
			cielab_avg_intern = cielab_avg_fallback;
			oklab_avg_intern  = oklab_avg_fallback;
		}
	}
}

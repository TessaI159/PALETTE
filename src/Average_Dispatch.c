#include "Average.h"
#include "Color.h"
#include "FeatureDetection.h"

static struct Color (*cielab_avg_intern)(const struct Color *colors,
				      uint16_t num_col) = cielab_avg_fallback;
void average_init() {
	if (!features.initialized) {
		query_features(&features);
	}
	if (features.avx2) {
		cielab_avg_intern = cielab_avg_avx2;
	} else if (features.avx) {
		cielab_avg_intern = cielab_avg_avx;
	} else {
		cielab_avg_intern = cielab_avg_fallback;
	}
}

struct Color cielab_avg(const struct Color *colors, uint16_t num_col) {
  static bool average_initialized = false;
	if (!average_initialized) {
		average_init();
		average_initialized = true;
	}
	return cielab_avg_intern(colors, num_col);
}

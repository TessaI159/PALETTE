#include "Average.h"
#include "Color.h"
#include "FeatureDetection.h"

static struct Color (*lab_avg_intern)(const struct Color *colors,
				      uint16_t num_col) = lab_avg_fallback;
void average_init() {
	if (!features.initialized) {
		query_features(&features);
	}
	if (features.avx2) {
		lab_avg_intern = lab_avg_avx2;
	} else if (features.avx) {
		lab_avg_intern = lab_avg_avx;
	} else {
		lab_avg_intern = lab_avg_fallback;
	}
}

struct Color lab_avg(const struct Color *colors, uint16_t num_col) {
  static bool average_initialized = false;
	if (!average_initialized) {
		average_init();
		average_initialized = true;
	}
	return lab_avg_intern(colors, num_col);
}

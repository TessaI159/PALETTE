#include <stdio.h>
#include <stdlib.h>

#include "Average.h"
#include "Average_Internal.h"
#include "Color.h"
#include "FeatureDetection.h"
#include "Parameters.h"

static struct Color (*color_avg_intern)(const struct Color);

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

struct Color color_avg(const struct Color colors) {
	check_initialized();
	return color_avg_intern(colors);
}

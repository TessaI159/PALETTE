#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Align.h"
#include "Average.h"
#include "Color.h"
#include "FeatureDetection.h"
#include "Util.h"
#include "Video.h"

/* TODO (Tess): Consider throwing GPU compute capabilities in here somewhere */

#define RUNS pow(2, 22)
#define RUNS_AVG pow(2, 15)
#define NUM_TEST_COL 10000

extern struct system_features_t features;
extern struct Color cielab_avg_fallback(const struct cielab_SoA *colors, uint16_t num_col);
extern struct Color cielab_avg_sse3(const struct cielab_SoA *colors,
				    uint16_t		     num_col);

static inline int width_from_pixels(int pixels) {
	return round(sqrt(pixels) * 4.0f / 3.0f);
}

static inline int height_from_width(int width) {
	return round(9.0f / 16.0f * width);
}

int main(int argc, char **argv) {
	srand(0);
	query_features(&features);

	printf("Color is %zu bytes\n", sizeof(struct Color));

	Color col1 = Color_create(rand() % 255, rand() % 255, rand() % 255);
	Color col2 = Color_create(rand() % 255, rand() % 255, rand() % 255);

	printf("Delta ok took %lu cycles on average.\n",
	       time_diff(delta_ok_diff, &col1, &col2, RUNS));
	printf("Delta cie76 took %lu cycles on average.\n",
	       time_diff(delta_cie76_diff, &col1, &col2, RUNS));
	printf("Delta cie94 took %lu cycles on average.\n",
	       time_diff(delta_cie94_diff, &col1, &col2, RUNS));
	printf("Delta ciede2000 took %lu cycles on average.\n",
	       time_diff(delta_ciede2000_diff, &col1, &col2, RUNS));

	printf("%" PRIu64 " bytes of l1 cache\n", features.l[1]);
	printf("%" PRIu64 " bytes of l2 cache\n", features.l[2]);
	printf("%" PRIu64 " bytes of l3 cache\n", features.l[3]);

	uint64_t cache	= features.l[1];
	uint64_t pixels = cache / sizeof(Color);
	int	 width	= width_from_pixels(pixels);
	int	 height = height_from_width(width);
	printf("%" PRIu64 " bytes of l1 cache can hold ~%" PRIu64
	       " colors. This means that if we want each frame to be held "
	       "entirely in l1 cache at a 16:9 ratio, the frame must be shrunk "
	       "to ~%dx%d\n",
	       cache, pixels, width, height);

	cache += features.l2_shared ? features.l[2]
				    : features.l[2] / features.physical_cores;
	pixels = cache / sizeof(Color);
	width  = width_from_pixels(pixels);
	height = height_from_width(width);

	printf("%" PRIu64 " bytes of l1 and l2 cache can hold ~%" PRIu64
	       " colors. This means that if we want each frame to be held "
	       "entirely in l1 and l2 cache at a 16:9 ratio, the frame must be "
	       "shrunk to ~%dx%d\n",
	       cache, pixels, width, height);

	cache += features.l[3] / features.physical_cores;
	pixels = cache / sizeof(Color);
	width  = width_from_pixels(pixels);
	height = height_from_width(width);

	printf("%" PRIu64 " bytes of l1, l2, and l3 cache can hold ~%" PRIu64
	       " colors. This means that if we want each frame to be held "
	       "entirely in l1, l2, and l3 cache at a 16:9 ratio, the frame "
	       "must be shrunk to ~%dx%d\n",
	       cache, pixels, width, height);

	/* struct Video in_video; */
	/* open_video_file(&in_video, "vid.webm"); */
	struct Color test_colors[NUM_TEST_COL] = {};
	for (size_t i = 0; i < NUM_TEST_COL; ++i) {
		test_colors[i] =
		    Color_create(rand() % 255, rand() % 255, rand() % 255);
		convert_to(&test_colors[i], COLOR_CIELAB);
	}
	cielab_SoA soa32;
	soa32.l = aligned_malloc(32, NUM_TEST_COL * sizeof(float));
	soa32.a = aligned_malloc(32, NUM_TEST_COL * sizeof(float));
	soa32.b = aligned_malloc(32, NUM_TEST_COL * sizeof(float));
	
	uint64_t avg =
	    time_avg(cielab_avg_fallback, &soa32, NUM_TEST_COL, RUNS_AVG);
	printf("Fallback averaging %" PRIu16
	       " colors took %lu cycles on average, or about %.3f "
	       "cycles per color on average.\n",
	       NUM_TEST_COL, avg, (double)avg / (double)NUM_TEST_COL);
	avg = time_avg(cielab_avg_sse3, &soa32, NUM_TEST_COL, RUNS_AVG);
	printf("SSE3 averaging %" PRIu16
	       " colors took %lu cycles on average, or about %.3f "
	       "cycles per color on average.\n",
	       NUM_TEST_COL, avg, (double)avg / (double)NUM_TEST_COL);
	avg = time_avg(cielab_avg_avx, &soa32, NUM_TEST_COL, RUNS_AVG);
	printf("AVX averaging %" PRIu16
	       " colors took %lu cycles on average, or about %.3f "
	       "cycles per color on average.\n",
	       NUM_TEST_COL, avg, (double)avg / (double)NUM_TEST_COL);
	
	aligned_free(soa32.l);
	aligned_free(soa32.a);
	aligned_free(soa32.b);

	return 0;
}

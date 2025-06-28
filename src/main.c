#include <SDL3/SDL_thread.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "Color.h"
#include "FeatureDetection.h"
#include "Util.h"

#define RUNS 32768

int main(int argc, char** argv) {
	srand(0);

	struct system_features_t features;
	memset(&features, 0, sizeof(struct system_features_t));

	query_cache_sizes(&features);
	query_cpu_features(&features);
	query_logical_cores(&features);
	query_physical_cores(&features);

	printf("sRGB Euclidean diff takes an average of %.5f cycles\n",
	       (double)test_diff_speed(euclidean_diff, RUNS) / (double)RUNS);

	printf("RedMean diff takes an average of %.5f cycles\n",
	       (double)test_diff_speed(redmean_diff, RUNS) / (double)RUNS);

	printf("deltaok diff takes an average of %.5f cycles\n",
	       (double)test_diff_speed(delta_ok_diff, RUNS) / (double)RUNS);

	printf("cie76 diff takes an average of %.5f cycles\n",
	       (double)test_diff_speed(delta_cie76_diff, RUNS) / (double)RUNS);

	printf("cie94 diff takes an average of %.5f cycles\n",
	       (double)test_diff_speed(delta_cie94_diff, RUNS) / (double)RUNS);

	printf("ciede2000 diff takes an average of %.5f cycles\n",
	       (double)test_diff_speed(delta_ciede2000_diff, RUNS) /
		   (double)RUNS);

	printf("Each color is %" PRIu64 " bytes.\n", Color_size());

	printf("Your system %s AVX.\n",
	       features.avx ? "supports" : "does not support");
	printf("Your system %s AVX2.\n",
	       features.avx2 ? "supports" : "does not support");
	printf("Your system %s fma3.\n",
	       features.fma3 ? "supports" : "does not support");
	printf("Your system %s neon\n",
	       features.neon ? "supports" : "does not support");
	printf("Your system has %" PRIu64 " KiB l1 per core, %" PRIu64
	       " KiB l2 per core, %" PRIu64 " KiB l3 total\n.",
	       features.l1 / 1024, features.l2 / 1024, features.l3 / 1024);
	printf("Your system has %" PRIu8 " cores and can run %" PRIu8
	       " threads.\n",
	       features.physical_cores, features.logical_cores);

	return 0;
}

/* Begin color averaging functions */

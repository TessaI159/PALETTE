#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Color.h"
#include "FeatureDetection.h"
#include "Util.h"

#define RUNS 131072

int main(int argc, char **argv) {
	srand(0);

	struct system_features_t features;
	memset(&features, 0, sizeof(struct system_features_t));

	query_cache_sizes(&features);
	query_cpu_features(&features);
	query_logical_cores(&features);
	query_physical_cores(&features);

	printf("Color is %zu bytes\n", sizeof(struct Color));

	Color col1 = Color_create(rand() * 255, rand() * 255, rand() * 255);
	Color col2 = Color_create(rand() * 255, rand() * 255, rand() * 255);
	Color_calc_spaces(&col1);
	Color_calc_spaces(&col2);
	printf("Euclidean took %lu cycles on average.\n",
	       time_diff(euclidean_diff, col1, col2, RUNS));
	printf("Redmean took %lu cycles on average.\n",
	       time_diff(redmean_diff, col1, col2, RUNS));
	printf("Delta ok took %lu cycles on average.\n",
	       time_diff(delta_ok_diff, col1, col2, RUNS));
	printf("Delta cie76 took %lu cycles on average.\n",
	       time_diff(delta_cie76_diff, col1, col2, RUNS));
	printf("Delta cie94 took %lu cycles on average.\n",
	       time_diff(delta_cie94_diff, col1, col2, RUNS));
	printf("Delta ciede2000 took %lu cycles on average.\n",
	       time_diff(delta_ciede2000_diff, col1, col2, RUNS));

	return 0;
}

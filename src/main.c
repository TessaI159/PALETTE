#include <SDL3/SDL_thread.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>

#include "FeatureDetection.h"



int main(int argc, char** argv) {
	srand(0);

	struct system_features_t features;
	memset(&features, 0, sizeof(struct system_features_t));

	query_cache_sizes(&features);
	query_cpu_features(&features);
	query_logical_cores(&features);
	query_physical_cores(&features);


	return 0;
}


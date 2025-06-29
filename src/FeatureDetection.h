#ifndef FEATURE_DETECTION_H
#define FEATURE_DETECTION_H

#include <stdint.h>

struct system_features_t {
	int	 avx;
	int	 avx2;
	int	 fma3;
	uint64_t l1;
	uint64_t l2;
	uint64_t l3;
	uint8_t	 logical_cores;
	uint8_t	 physical_cores;
	uint16_t line_size;
};

void query_cache_sizes(struct system_features_t *features);
void query_cpu_features(struct system_features_t *features);
void query_logical_cores(struct system_features_t *features);
void query_physical_cores(struct system_features_t *features);

#endif

#ifndef FEATURE_DETECTION_H
#define FEATURE_DETECTION_H

#include <stdint.h>

struct system_features_t {
	int	 avx;
	int	 avx2;
	int	 fma3;
	uint64_t l[3];
	uint8_t	 logical_cores;
	uint8_t	 physical_cores;
	uint16_t line_size;
	int	 l2_shared;
};

void query_features(struct system_features_t *features);

#endif

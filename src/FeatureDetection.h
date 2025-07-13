#ifndef FEATURE_DETECTION_H
#define FEATURE_DETECTION_H

#include <stdbool.h>
#include <stdint.h>

struct system_features_t {
	bool	 sse2;
	bool	 avx;
	bool	 avx2;
	bool	 fma3;
	uint64_t l[4];
	uint8_t	 logical_cores;
	uint8_t	 physical_cores;
	uint16_t line_size;
	bool	 l2_shared;
	bool	 initialized;
};

extern struct system_features_t features;

void query_features(struct system_features_t *features);

#endif

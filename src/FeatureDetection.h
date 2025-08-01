#ifndef FEATURE_DETECTION_H
#define FEATURE_DETECTION_H

#include <stdbool.h>
#include <stdint.h>

enum InstSet { AVX, SSE, FB };

struct system_features_t {
	bool	     sse;
	bool	     avx;
	bool	     fma3;
	uint64_t     l[4];
	uint8_t	     logical_cores;
	uint8_t	     physical_cores;
	uint16_t     line_size;
	bool	     l2_shared;
	bool	     initialized;
	enum InstSet instructions;
};

extern struct system_features_t features;
extern enum InstSet		inst_set;

void query_features(struct system_features_t *features);

#endif

#ifndef FEATURE_DETECTION_H
#define FEATURE_DETECTION_H

#include <stdbool.h>
#include <stdint.h>

enum InstSet { AVX, SSE, FB };

struct system_features_t {
  bool         sse;            /*<< Whether or not the system supports SSE */
  bool         avx;            /*<< Whether or not the system supports AVX */
  bool         fma3;           /*<< Whether or not the system supports FMA */
  uint64_t     l[4];           /*<< Cache sizes in bytes */
  uint8_t      logical_cores;  /*<< Logical core count */
  uint8_t      physical_cores; /*<< Physical core count */
  uint16_t     line_size;      /*<< Size of a cache line */
  bool         l2_shared;      /*<< Indicates if cores share an L2 pool */
  bool         initialized;
  enum InstSet instructions;
};

extern struct system_features_t features;
extern enum InstSet             inst_set;

/**
 * @brief queries the system to determine what features in system_features_t are
 * available
 * @param features A point to a system_features_t object
 */
void query_features(struct system_features_t *features);

#endif

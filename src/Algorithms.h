#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <stdbool.h>
#include <stdint.h>

#include "Color.h"

struct kmeans_params {
  uint8_t k;          /**< Number of centroids */
  uint32_t max_iters; /**< Max clustering iterations */
  float tolerance;    /**< Centroid movement threshold (e.g. 1e-4f) */
  uint32_t seed;      /**< PRNG seed for initialization  */
};

struct kmeans_result {
  struct Color centroids; /**< Extracted k centroids in SoA */
  uint8_t *assignments;   /**< Cluster id [0, k-1] for each input pixel */
  uint32_t iters_run;     /**< Total iterations required for convergence */
  float inertia;          /**< Sum of squared distances to centroids */
  uint32_t num_pixels;
};

/**
 * @brief Runs KMeans clustering on SoA color data.
 * @param colors Input color array in SoA format.
 * @param params Clustering parameters.
 * @param result Output struct containing centroids and assignments.
 * @return true on success, false on error.
 */
bool kmeans_run(const struct Color *colors, const struct kmeans_params *params,
                struct kmeans_result *result);

/**
 *@brief Frees all memory allocated inside a kmeans_result
 */
void kmeans_result_free(struct kmeans_result *result);

#endif

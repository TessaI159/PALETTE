#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Algorithms.h"
#include "Average.h"
#include "Difference.h"

#define SQUARE(x) ((x) * (x))

/**
 * @brief Psuedo-random integer generator
 * @note See (Marsaglia, 2003)
 */

/* TODO (Tess): If you're ever bored one day, generalize the KMeans algorithm
   to run on data of arbitrary dimension with an arbitrary(ish) distance
   function */

static inline uint32_t xorshift32(uint32_t *state) {
  uint32_t x = *state;
  if (x == 0)
    x = 2463534242U;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *state = x;
  return x;
}

static inline void assign_pixels(struct kmeans_result *result, float *diffs,
                                 uint32_t n, uint8_t k) {
  float total_inertia = 0.0f;
  for (uint32_t i = 0; i < n; ++i) {
    uint8_t best_j = 0;
    float min_d = diffs[0 * n + i];
    for (uint8_t j = 1; j < k; ++j) {
      float d = diffs[j * n + i];
      if (d < min_d) {
        min_d = d;
        best_j = j;
      }
    }
    result->assignments[i] = best_j;
    total_inertia += SQUARE(min_d);
  }
  result->inertia = total_inertia;
}

bool kmeans_run(const struct Color *colors, const struct kmeans_params *params,
                struct kmeans_result *result) {
  if (!colors || colors->num == 0 || !params || params->k == 0 || !result)
    return false;

  uint32_t n = colors->num;
  uint8_t k = params->k;
  if (k > n)
    k = (uint8_t)n;

  uint32_t prng_state = params->seed ? params->seed : 123456789;

  /* Memory initialization */
  memset(result, 0, sizeof(*result));
  if (!(Color_alloc(&result->centroids, k)))
    return false;

  result->assignments = (uint8_t *)malloc(n * sizeof(uint8_t));
  if (!result->assignments) {
    kmeans_result_free(result);
    return false;
  }
  result->num_pixels = n;

  float *diffs = (float *)malloc((size_t)k * n * sizeof(float));
  if (!diffs) {
    kmeans_result_free(result);
    return false;
  }

  struct Color new_cents;
  if (!Color_alloc(&new_cents, k)) {
    free(diffs);
    kmeans_result_free(result);
    return false;
  }

  struct Color cluster_scratch;
  if (!Color_alloc(&cluster_scratch, n)) {
    Color_free(&new_cents);
    free(diffs);
    kmeans_result_free(result);
    return false;
  }

  /* KMeans++ initialization  */
  uint32_t first_idx = xorshift32(&prng_state) % n;
  result->centroids.alpha[0] = colors->alpha[first_idx];
  result->centroids.beta[0] = colors->beta[first_idx];
  result->centroids.gamma[0] = colors->gamma[first_idx];

  for (uint8_t c = 1; c < k; ++c) {
    struct Color active_cents = {.alpha = result->centroids.alpha,
                                 .beta = result->centroids.beta,
                                 .gamma = result->centroids.gamma,
                                 .num = c};
    delta_diff(colors, &active_cents, diffs);

    double sum_dist_sq = 0.0;
    for (uint32_t i = 0; i < n; ++i) {
      float min_d = diffs[0 * n + i];
      for (uint8_t j = 1; j < c; ++j) {
        float d = diffs[j * n + i];
        if (d < min_d)
          min_d = d;
      }
      sum_dist_sq += (double)(SQUARE(min_d));
    }
    double threshold =
        ((double)(xorshift32(&prng_state) % 1000000) / 1000000.0) * sum_dist_sq;
    double running = 0.0;
    uint32_t chosen = n - 1;
    for (uint32_t i = 0; i < n; ++i) {
      float min_d = diffs[0 * n + i];
      for (uint8_t j = 1; j < c; ++j) {
        float d = diffs[j * n + i];
        if (d < min_d)
          min_d = d;
      }
      running += (double)(SQUARE(min_d));
      if (running > threshold) {
        chosen = i;
        break;
      }
    }
    result->centroids.alpha[c] = colors->alpha[chosen];
    result->centroids.beta[c] = colors->beta[chosen];
    result->centroids.gamma[c] = colors->gamma[chosen];
  }

  /* The actual KMeans algorithm */
  uint32_t max_iters = params->max_iters > 0 ? params->max_iters : 50;
  float tol = params->tolerance > 0.0f ? params->tolerance : 1e-4f;
  float tol_sq = SQUARE(tol);

  for (uint32_t iter = 0; iter < max_iters; ++iter) {
    result->iters_run = iter + 1;

    /* Step 1: Compute all pairwise distances */
    delta_diff(colors, &result->centroids, diffs);

    /* Step 2: Assign each pixel to its nearest centroid  */
    assign_pixels(result, diffs, n, k);

    /* Step 3: Update centroids  */
    float max_movement_sq = 0.0f;
    for (uint8_t j = 0; j < k; ++j) {
      uint32_t count = 0;
      for (uint32_t i = 0; i < n; ++i) {
        if (result->assignments[i] == j) {
          cluster_scratch.alpha[count] = colors->alpha[i];
          cluster_scratch.beta[count] = colors->beta[i];
          cluster_scratch.gamma[count] = colors->gamma[i];
          count++;
        }
      }
      if (count > 0) {
        cluster_scratch.num = count;
        color_avg(&cluster_scratch, &new_cents, j);
      } else {
        /*  Cluster is empty. Re-seed to a new random pixel. */
        /* TODO (Tess): Make sure the new random color doesn't accidentaly
           match another centroid exactly. */
        uint32_t rand_idx = xorshift32(&prng_state) % n;
        new_cents.alpha[j] = colors->alpha[rand_idx];
        new_cents.beta[j] = colors->beta[rand_idx];
        new_cents.gamma[j] = colors->gamma[rand_idx];
      }

      float movement_sq =
          SQUARE(new_cents.alpha[j] - result->centroids.alpha[j]) +
          SQUARE(new_cents.beta[j] - result->centroids.beta[j]) +
          SQUARE(new_cents.gamma[j] - result->centroids.gamma[j]);
      if (movement_sq > max_movement_sq) {
        max_movement_sq = movement_sq;
      }
    }
    /* Copy newly calculated centroids to result  */
    memcpy(result->centroids.alpha, new_cents.alpha, k * sizeof(float));
    memcpy(result->centroids.beta, new_cents.beta, k * sizeof(float));
    memcpy(result->centroids.gamma, new_cents.gamma, k * sizeof(float));

    /* Check for convergence  */
    if (max_movement_sq < tol_sq) {
      break;
    }
  }
  delta_diff(colors, &result->centroids, diffs);
  assign_pixels(result, diffs, n, k);
  Color_free(&cluster_scratch);
  Color_free(&new_cents);
  free(diffs);
  return true;
}

void kmeans_result_free(struct kmeans_result *result) {
  if (!result)
    return;
  Color_free(&result->centroids);
  if (result->assignments) {
    free(result->assignments);
    result->assignments = NULL;
  }
  result->num_pixels = 0;
  result->iters_run = 0;
  result->inertia = 0.0f;
}

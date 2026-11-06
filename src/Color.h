#ifndef COLOR_H
#define COLOR_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

struct Color {
  float *alpha; /* e.g. R in sRGB, L in Lab */
  float *beta;  /* e.g. G in sRGB, a in Lab */
  float *gamma; /* e.g. B in sRGB, b in Lab */
  uint32_t num;
};

/**
 * @brief Allocates 32-byte aligned memory for SoA color buffers.
 * @return true on success, false on failure.
 */
bool Color_alloc(struct Color *c, uint32_t num);

/**
 * @brief Frees aligned SoA buffers and zeroes out the struct.
 */
void Color_free(struct Color *c);

/**
 * @brief Deep clones an SoA Color struct into an allocated destination.
 */
bool Color_clone(struct Color *dest, const struct Color *src);

/**
 * @brief Sets values at index idx in the SoA buffers.
 */
static inline void Color_set(struct Color *c, uint32_t idx, float a, float b,
                             float g) {
  c->alpha[idx] = a;
  c->beta[idx] = b;
  c->gamma[idx] = g;
}

/**
 * @brief Gets values at index idx in the Soa buffers.
 */
static inline void Color_get(struct Color *c, uint32_t idx, float *a, float *b,
                             float *g) {
  if (a)
    *a = c->alpha[idx];
  if (b)
    *b = c->beta[idx];
  if (g)
    *g = c->gamma[idx];
}

#endif

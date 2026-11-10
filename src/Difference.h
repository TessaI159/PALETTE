#ifndef DIFFERENCE_H
#define DIFFERENCE_H

#include <stdint.h>

/* TODO (Tess): I regret calling these differences and not distances. Change? */

typedef struct Color Color;

/**
 * @brief Uses a dispatch table to find all pairwise differences of color and
 * cents
 * @param colors An SoA of color objects
 * @param cents An Soa of color objects
 * @param diffs An array capable of holding colors->num*cents->num floats
 * @note The float at [a * colors->num + b] is the difference between centroid a
 * and color b
 */
void delta_diff(const struct Color *restrict colors,
                const struct Color *restrict cents, float *restrict diffs);

#endif

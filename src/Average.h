#ifndef AVERAGE_H
#define AVERAGE_H

#include <stdint.h>

typedef struct Color Color;

/**
 * @brief Finds the average color of all Colors in colors and updates
 * cents[which_cent] to that color
 * @param colors An SoA of Color objects
 * @param cents An SoA of Color objects representing centroids. Modifies the
 * which_centh centroid to the average of colors
 * @param which_cent An array index indicating which centroid you would like to
 * update
 */
void color_avg(const struct Color *restrict colors,
               struct Color *restrict cents, uint8_t which_cent);

#endif

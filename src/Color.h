#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

typedef struct Color Color;
enum ColorSpace;

uint64_t Color_size();

/* Caller must call Color_destroy when finished with this color.
 * Returns NULL on failure */
Color *Color_create(uint8_t a, uint8_t b, uint8_t c);

/* Caller must set color = NULL after calling */
void Color_destroy(Color *color);

/* Calculates all spaces. Used for testing */
void Color_calc_spaces(Color *color);

/* Euclidean srgb diff. ~30 cycles*/
double euclidean_diff(Color *sam, Color *ref);
/* Euclidean diff with color weights. ~30 cycles */
double redmean_diff(Color *sam, Color *ref);
/* Euclidean oklab diff. ~60 cycles */
double delta_ok_diff(Color *sam, Color *ref);
/* Euclidean cielab diff. ~500 cycles */
double delta_cie76_diff(Color *sam, Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIE94 ~550 cycles */
double delta_cie94_diff(Color *sam, Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIEDE2000 ~1600 cycles */
double delta_ciede2000_diff(Color *sam, Color *ref);
#endif

#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

enum ColorSpace { COLOR_OKLAB = 0, COLOR_CIELAB, COLOR_SRGB, COLOR_GRAY };

/* Supported color spaces */
struct okLAB {
	float l, a, b;
};

/* non-linear */
struct sRGB {
	float r, g, b;
};


struct cieLAB {
	float l, a, b;
};

/* non-linear */
struct Grayscale {
	float l;
};

struct Color {
	struct okLAB	 oklab;
	struct cieLAB	 cielab;
	struct sRGB	 srgb;
	struct Grayscale grayscale;
	struct sRGB	 lsrgb;
	uint8_t		 valid_spaces;
};

/* Must be constructed with rgb for now */
struct Color Color_create(const uint8_t r, const uint8_t g, const uint8_t b);

/* Calculates all spaces. Used for testing */
void Color_calc_spaces(struct Color *color);

/* Euclidean srgb diff. ~30 cycles*/
float euclidean_diff(struct Color *sam, struct Color *ref);
/* Euclidean diff with color weights. ~35 cycles */
float redmean_diff(struct Color *sam, struct Color *ref);
/* Euclidean oklab diff. ~35 cycles */
float delta_ok_diff(struct Color *sam, struct Color *ref);
/* Euclidean cielab diff. ~45 cycles */
float delta_cie76_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIE94 ~90 cycles */
float delta_cie94_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIEDE2000 ~900 cycles */
float delta_ciede2000_diff(struct Color *sam, struct Color *ref);
#endif

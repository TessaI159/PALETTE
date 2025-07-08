#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

enum ColorSpace {
	COLOR_OKLAB  = 1U,
	COLOR_CIELAB = 2U,
	COLOR_SRGB   = 4U,
	COLOR_GRAY   = 8U
};

/* Supported color spaces */
struct okLAB {
	double l, a, b;
};

/* Stored linear and normalized */
struct sRGB {
	double r, g, b;
};

struct cieLAB {
	double l, a, b;
};

struct Grayscale {
	double l;
};

struct Color {
	struct okLAB	 oklab;
	struct cieLAB	 cielab;
	struct sRGB	 srgb;
	struct Grayscale grayscale;
	uint8_t		 valid_spaces;
};

/* Must be constructed with rgb for now */
struct Color Color_create(const double r, const double g, const double b);

struct Color Color_create_norm(const double r, const double g, const double b);

/* Calculates all spaces. Used for testing */
void Color_calc_spaces(struct Color *color);

/* Euclidean oklab diff. ~35 cycles */
double delta_ok_diff(struct Color *sam, struct Color *ref);
/* Euclidean cielab diff. ~45 cycles */
double delta_cie76_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIE94 ~100 cycles */
double delta_cie94_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIEDE2000 ~450 cycles */
double delta_ciede2000_diff(struct Color *sam, struct Color *ref);
#endif

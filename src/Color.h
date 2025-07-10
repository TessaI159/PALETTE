#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

enum ColorSpace { COLOR_SRGB, COLOR_OKLAB, COLOR_CIELAB, COLOR_GRAY };

/* Supported color spaces */
struct okLAB {
	float l, a, b;
};

/* Stored linear and normalized */
struct sRGB {
	float r, g, b;
};

struct cieLAB {
	float l, a, b;
};

struct Grayscale {
	float l;
};


struct Color {
	struct okLAB	 oklab;
	struct cieLAB	 cielab;
	struct sRGB	 srgb;
	struct Grayscale grayscale;
	uint8_t		 valid_spaces;
};

/* Must be constructed with rgb for now */
struct Color Color_create(const float r, const float g, const float b);

struct Color Color_create_norm(const float r, const float g, const float b);

/* Calculates all spaces. Used for testing */
void Color_calc_spaces(struct Color *color);


/* Euclidean oklab diff. ~35 cycles */
float delta_ok_diff(struct Color *sam, struct Color *ref);
/* Euclidean cielab diff. ~45 cycles */
float delta_cie76_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIE94 ~100 cycles */
float delta_cie94_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIEDE2000 ~450 cycles */
float delta_ciede2000_diff(struct Color *sam, struct Color *ref);
#endif

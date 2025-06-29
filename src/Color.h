#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

enum ColorSpace {
	COLOR_OKLAB = 1,
	COLOR_CIELAB,
	COLOR_SRGB,
	COLOR_XYZ,
};

/* Supported color spaces */
struct okLAB {
	double l, a, b;
};

struct sRGB {
	double r, g, b;
};

struct XYZ {
	double x, y, z;
};

struct cieLAB {
	double l, a, b;
};

struct Color {
	struct okLAB  oklab;
	struct cieLAB cielab;
	struct sRGB   srgb;
	struct XYZ    xyz;
	uint32_t      valid_spaces;
};

/* Must be constructed with rgb for now */
struct Color Color_create(const char r, const char g, const char b);

/* Calculates all spaces. Used for testing */
void Color_calc_spaces(struct Color *color);

/* Euclidean srgb diff. ~30 cycles*/
double euclidean_diff(const struct Color sam, const struct Color ref);
/* Euclidean diff with color weights. ~30 cycles */
double redmean_diff(const struct Color sam, const struct Color ref);
/* Euclidean oklab diff. ~60 cycles */
double delta_ok_diff(struct Color sam, struct Color ref);
/* Euclidean cielab diff. ~500 cycles */
double delta_cie76_diff(struct Color sam, struct Color ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIE94 ~550 cycles */
double delta_cie94_diff(struct Color sam, struct Color ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIEDE2000 ~1600 cycles */
double delta_ciede2000_diff(struct Color sam, struct Color ref);
#endif

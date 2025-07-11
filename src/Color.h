#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

enum ColorSpace { COLOR_SRGB, COLOR_OKLAB, COLOR_CIELAB };

typedef struct E2000_diff E200_diff;

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

struct Color {
	enum ColorSpace current_space;
	union {
		struct okLAB  oklab;
		struct cieLAB cielab;
		struct sRGB   srgb;
	} data;
};

#ifdef PALETTE_DEBUG
void  Color_print(struct Color *color);
float delta_ciede2000_diff_fast(struct E2000_diff *diff);
void  convert_to(struct Color *color, enum ColorSpace t);
#endif

/* Must be constructed with rgb for now */
struct Color Color_create(const float r, const float g, const float b);
struct Color Color_create_norm(const float r, const float g, const float b);
struct Color Color_create_lab(const float l, const float a, const float b);

/* Euclidean oklab diff. ~35 cycles */
float delta_ok_diff(struct Color *sam, struct Color *ref);
/* Euclidean cielab diff. ~45 cycles */
float delta_cie76_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIE94 ~100 cycles */
float delta_cie94_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIEDE2000 ~450 cycles */
float delta_ciede2000_diff(struct Color *sam, struct Color *ref);
#endif

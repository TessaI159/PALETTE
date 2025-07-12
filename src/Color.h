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

/* For debugging purposes only */
#ifdef PALETTE_DEBUG
void  Color_print(struct Color *color);
float delta_ciede2000_diff_fast(struct E2000_diff *diff);
void  convert_to(struct Color *color, enum ColorSpace t);
#endif

/**
 *@brief Creates a new instance of a Color struct. No dynamic memory.
 *
 * Creates and returns a new Color struct instantiated with sRGB values r, g,
 * b. sRGB values are always stored linearized and normalized.
 *
 *@param r A floating point representation of the r channel in [0,255]
 *
 *@param g A floating point representation of the g channel in [0,255]
 *
 *@param b A floating point representation of the b channel in [0,255]
 *
 *@return A Color struct instantiated with sRGB values r, g, b
 **/
struct Color Color_create(const float r, const float g, const float b);

/**
 *@brief Creates a new instance of a Color struct. No dynamic memory.
 *
 * Creates and returns a new Color struct instantiated with sRGB values r, g,
 * b. sRGB values are always stored linearized and normalized.
 *
 *@param r A floating point representation of the r channel in [0,1]
 *
 *@param g A floating point representation of the g channel in [0,1]
 *
 *@param b A floating point representation of the b channel in [0,1]
 *
 *@return A Color struct instantiated with sRGB values r, g, b
 **/
struct Color Color_create_norm(const float r, const float g, const float b);
struct Color Color_create_cielab(const float l, const float a, const float b);
struct Color Color_create_oklab(const float l, const float a, const float b);

/* Euclidean oklab diff. ~35 cycles */
float delta_ok_diff(struct Color *sam, struct Color *ref);
/* Euclidean cielab diff. ~45 cycles */
float delta_cie76_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIE94 ~100 cycles */
float delta_cie94_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIEDE2000 ~450 cycles */
float delta_ciede2000_diff(struct Color *sam, struct Color *ref);
#endif

#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>
#include "Align.h"

enum ColorSpace { COLOR_SRGB, COLOR_OKLAB, COLOR_CIELAB };

typedef struct E2000_diff E200_diff;

/* Supported color spaces */
struct okLAB {
	float l, a, b;
	float _pad;
};

/* Stored linear and normalized */
struct sRGB {
	float r, g, b;
	float _pad;
};

struct cieLAB {
	float l, a, b;
	float _pad;
};

struct Color {
	union {
		struct okLAB  oklab;
		struct cieLAB cielab;
		struct sRGB   srgb;
	} data;
	enum ColorSpace current_space;
};

struct cielab_SoA {
  float *l;
  float *a;
  float *b;
};

struct oklab_SoA {
  float *l;
  float *a;
  float *b;
};

/* For debugging purposes only */
#ifdef PALETTE_DEBUG
void  Color_print(struct Color *color);
float delta_ciede2000_diff_fast(struct E2000_diff *diff);
#endif

void convert_to(struct Color *color, enum ColorSpace t);

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

/**
 *@brief Creates a new instance of a Color struct. No dynamic memory.
 *
 * Creates and returns a new Color struct instantiated with cieLAB values l, a,
 * b.
 *
 *@param l A floating point representation of the l channel in [0,100]
 *
 *@param a A floating point representation of the a channel: unbounded
 *
 *@param b A floating point representation of the b channel: unbounded
 *
 *@return A Color struct instantiated with cieLAB values l, a, b
 **/
struct Color Color_create_cielab(const float l, const float a, const float b);

/**
 *@brief Creates a new instance of a Color struct. No dynamic memory.
 *
 * Creates and returns a new Color struct instantiated with okLAB values l, a,
 * b.
 *
 *@param l A floating point representation of the l channel in [0,100]
 *
 *@param a A floating point representation of the a channel: unbounded
 *
 *@param b A floating point representation of the b channel: unbounded
 *
 *@return A Color struct instantiated with okLAB values l, a, b
 **/
struct Color Color_create_oklab(const float l, const float a, const float b);

/* Euclidean oklab diff. ~80 cycles debug, ~5 cycles release */
float delta_ok_diff(struct Color *sam, struct Color *ref);
/* Euclidean cielab diff. ~80 cycles debug, ~5 cycles release */
float delta_cie76_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIE94 ~160 cycles debug, ~40
   cycles releasex */
float delta_cie94_diff(struct Color *sam, struct Color *ref);
/* https://en.wikipedia.org/wiki/Color_difference#CIEDE2000 ~685 cycles debug,
   ~500 cycles release */
float delta_ciede2000_diff(struct Color *sam, struct Color *ref);
#endif

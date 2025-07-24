#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "Color.h"

#define SQUARE(x) ((x) * (x))

void cielab_avg_fb(const struct Color *restrict colors,
		   struct Color *restrict cents, uint8_t which_cent) {
	uint16_t num_col = colors->num;
	float	 l1, l2, l3, l4;
	float	 a1, a2, a3, a4;
	float	 b1, b2, b3, b4;
	l1 = l2 = l3 = l4 = 0.0f;
	a1 = a2 = a3 = a4 = 0.0f;
	b1 = b2 = b3 = b4 = 0.0f;

	for (size_t i = 0; i + 3 < num_col; i += 4) {
		l1 += colors->alpha[i];
		l2 += colors->alpha[i + 1];
		l3 += colors->alpha[i + 2];
		l4 += colors->alpha[i + 3];

		a1 += colors->beta[i];
		a2 += colors->beta[i + 1];
		a3 += colors->beta[i + 2];
		a4 += colors->beta[i + 3];

		b1 += colors->gamma[i];
		b2 += colors->gamma[i + 1];
		b3 += colors->gamma[i + 2];
		b4 += colors->gamma[i + 3];
	}
	float l = l1 + l2 + l3 + l4;
	float a = a1 + a2 + a3 + a4;
	float b = b1 + b2 + b3 + b4;

	const float scale	= 1.0f / num_col;
	cents->alpha[which_cent] = l * scale;
	cents->beta[which_cent]	= a * scale;
	cents->gamma[which_cent] = b * scale;
}

void cielab_avg_fb_cw(const struct Color *restrict colors,
		      struct Color *restrict cents, uint8_t which_cent) {
	uint16_t num_col = colors->num;
	float	 l	 = 0.0f;
	float	 a	 = 0.0f;
	float	 b	 = 0.0f;
	float	 w	 = 0.0f;
	float	 chroma	 = 0.0f;
	float	 tmpa;
	float	 tmpb;

	size_t i = 0;
	for (; i < num_col; ++i) {
		tmpa   = colors->beta[i];
		tmpb   = colors->gamma[i];
		chroma = sqrtf(SQUARE(tmpa) + SQUARE(tmpb));
		l += colors->alpha[i];
		a += tmpa * chroma;
		b += tmpb * chroma;
		w += chroma;
	}

	l /= num_col;
	if (w == 0.0f) {
		a = 0.0f;
		b = 0.0f;
	} else {
		a /= w;
		b /= w;
	}

	cents->alpha[which_cent] = l;
	cents->beta[which_cent]	= a;
	cents->gamma[which_cent] = b;
}

void oklab_avg_fb(const struct Color *restrict colors,
		  struct Color *restrict cents, uint8_t which_cent) {
	uint16_t num_col = colors->num;
	float	 l1, l2, l3, l4;
	float	 a1, a2, a3, a4;
	float	 b1, b2, b3, b4;
	l1 = l2 = l3 = l4 = 0.0f;
	a1 = a2 = a3 = a4 = 0.0f;
	b1 = b2 = b3 = b4 = 0.0f;

	for (size_t i = 0; i + 3 < num_col; i += 4) {
		l1 += colors->alpha[i];
		l2 += colors->alpha[i + 1];
		l3 += colors->alpha[i + 2];
		l4 += colors->alpha[i + 3];

		a1 += colors->beta[i];
		a2 += colors->beta[i + 1];
		a3 += colors->beta[i + 2];
		a4 += colors->beta[i + 3];

		b1 += colors->gamma[i];
		b2 += colors->gamma[i + 1];
		b3 += colors->gamma[i + 2];
		b4 += colors->gamma[i + 3];
	}
	float l	    = l1 + l2 + l3 + l4;
	float a	    = a1 + a2 + a3 + a4;
	float b	    = b1 + b2 + b3 + b4;
	float scale = 1.0f / num_col;

	cents->alpha[which_cent] = l * scale;
	cents->beta[which_cent]	= a * scale;
	cents->gamma[which_cent] = b * scale;
}

void oklab_avg_fb_cw(const struct Color *restrict colors,
		     struct Color *restrict cents, uint8_t which_cent) {
	uint16_t num_col = colors->num;
	float	 l	 = 0.0f;
	float	 a	 = 0.0f;
	float	 b	 = 0.0f;
	float	 w	 = 0.0f;
	float	 chroma	 = 0.0f;
	float	 tmpa;
	float	 tmpb;

	size_t i = 0;
	for (; i < num_col; ++i) {
		tmpa   = colors->beta[i];
		tmpb   = colors->gamma[i];
		chroma = sqrtf(SQUARE(tmpa) + SQUARE(tmpb));
		l += colors->alpha[i + 0];
		a += tmpa * chroma;
		b += tmpb * chroma;
		w += chroma;
	}

	l /= num_col;
	if (w == 0.0f) {
		a = 0.0f;
		b = 0.0f;
	} else {
		a /= w;
		b /= w;
	}
	cents->alpha[which_cent] = l;
	cents->beta[which_cent]	= a;
	cents->gamma[which_cent] = b;
}

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "Color.h"

#define SQUARE(x) ((x) * (x))

struct Color cielab_avg_fallback(const struct cielab_SoA *__restrict colors,
				 uint16_t num_col) {
	float l1, l2, l3, l4;
	float a1, a2, a3, a4;
	float b1, b2, b3, b4;
	l1 = l2 = l3 = l4 = 0.0f;
	a1 = a2 = a3 = a4 = 0.0f;
	b1 = b2 = b3 = b4 = 0.0f;

	for (size_t i = 0; i + 3 < num_col; i += 4) {
		l1 += colors->l[i];
		l2 += colors->l[i + 1];
		l3 += colors->l[i + 2];
		l4 += colors->l[i + 3];

		a1 += colors->a[i];
		a2 += colors->a[i + 1];
		a3 += colors->a[i + 2];
		a4 += colors->a[i + 3];

		b1 += colors->b[i];
		b2 += colors->b[i + 1];
		b3 += colors->b[i + 2];
		b4 += colors->b[i + 3];
	}
	float l = l1 + l2 + l3 + l4;
	float a = a1 + a2 + a3 + a4;
	float b = b1 + b2 + b3 + b4;
	l /= num_col;
	a /= num_col;
	b /= num_col;
	return Color_create_cielab(l, a, b);
}

struct Color cielab_avg_fallback_cw(const struct cielab_SoA *__restrict colors,
				    uint16_t num_col) {
	float l	     = 0.0f;
	float a	     = 0.0f;
	float b	     = 0.0f;
	float w	     = 0.0f;
	float chroma = 0.0f;
	float tmpa;
	float tmpb;

	size_t i = 0;
	for (; i < num_col; ++i) {
		tmpa   = colors->a[i];
		tmpb   = colors->b[i];
		chroma = sqrtf(SQUARE(tmpa) + SQUARE(tmpb));
		l += colors->l[i];
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

	return Color_create_cielab(l, a, b);
}

struct Color oklab_avg_fallback(const struct oklab_SoA *__restrict colors,
				 uint16_t num_col) {
	float l1, l2, l3, l4;
	float a1, a2, a3, a4;
	float b1, b2, b3, b4;
	l1 = l2 = l3 = l4 = 0.0f;
	a1 = a2 = a3 = a4 = 0.0f;
	b1 = b2 = b3 = b4 = 0.0f;

	for (size_t i = 0; i + 3 < num_col; i += 4) {
		l1 += colors->l[i];
		l2 += colors->l[i + 1];
		l3 += colors->l[i + 2];
		l4 += colors->l[i + 3];

		a1 += colors->a[i];
		a2 += colors->a[i + 1];
		a3 += colors->a[i + 2];
		a4 += colors->a[i + 3];

		b1 += colors->b[i];
		b2 += colors->b[i + 1];
		b3 += colors->b[i + 2];
		b4 += colors->b[i + 3];
	}
	float l = l1 + l2 + l3 + l4;
	float a = a1 + a2 + a3 + a4;
	float b = b1 + b2 + b3 + b4;
	l /= num_col;
	a /= num_col;
	b /= num_col;
	return Color_create_oklab(l, a, b);
}

struct Color oklab_avg_fallback_cw(const struct oklab_SoA *__restrict colors,
				   uint16_t num_col) {
	float l	     = 0.0f;
	float a	     = 0.0f;
	float b	     = 0.0f;
	float w	     = 0.0f;
	float chroma = 0.0f;
	float tmpa;
	float tmpb;

	size_t i = 0;
	for (; i < num_col; ++i) {
		tmpa   = colors->a[i];
		tmpb   = colors->b[i];
		chroma = sqrtf(SQUARE(tmpa) + SQUARE(tmpb));
		l += colors->l[i + 0];
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

	return Color_create_oklab(l, a, b);
}

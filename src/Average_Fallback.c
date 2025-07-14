#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "Color.h"

struct Color cielab_avg_fallback(const struct cielab_SoA *__restrict colors,
				 uint16_t num_col) {
	float l, a, b;
	l = a = b = 0.0f;
	for (size_t i = 0; i + 3 < num_col; i += 4) {
		l += colors->l[i];
		l += colors->l[i + 1];
		l += colors->l[i + 2];
		l += colors->l[i + 3];

		a += colors->a[i];
		a += colors->a[i + 1];
		a += colors->a[i + 2];
		a += colors->a[i + 3];
		
		b += colors->b[i];
		b += colors->b[i + 1];
		b += colors->b[i + 2];
		b += colors->b[i + 3];
	}

	l /= num_col;
	a /= num_col;
	b /= num_col;
	return Color_create_cielab(l, a, b);
}

struct Color oklab_avg_fallback(const struct oklab_SoA *__restrict colors,
				 uint16_t num_col) {
	float l, a, b;
	l = a = b = 0.0f;
	for (size_t i = 0; i + 3 < num_col; i += 4) {
		l += colors->l[i];
		l += colors->l[i + 1];
		l += colors->l[i + 2];
		l += colors->l[i + 3];

		a += colors->a[i];
		a += colors->a[i + 1];
		a += colors->a[i + 2];
		a += colors->a[i + 3];
		
		b += colors->b[i];
		b += colors->b[i + 1];
		b += colors->b[i + 2];
		b += colors->b[i + 3];
	}

	l /= num_col;
	a /= num_col;
	b /= num_col;
	return Color_create_oklab(l, a, b);
}



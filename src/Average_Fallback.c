#include <stdint.h>
#include <stdio.h>

#include "Color.h"

struct Color lab_avg_fallback(struct Color *colors, uint16_t num_col) {
	float l = 0.0f;
	float a = 0.0f;
	float b = 0.0f;
	for (size_t i = 0; i < num_col; ++i) {
	  l += colors[i].cielab.l;
	  a += colors[i].cielab.a;
	  b += colors[i].cielab.b;
	}
	l /= (float)num_col;
	a /= (float)num_col;
	b /= (float)num_col;
	
	return Color_create_lab(l, a, b);
}

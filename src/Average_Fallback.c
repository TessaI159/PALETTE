#include <stdint.h>
#include <stdio.h>

#include "Color.h"

struct Color lab_avg_fallback(struct Color *colors, uint16_t num_col) {
	float total_l = 0.0f;
	float total_a = 0.0f;
	float total_b = 0.0f;
	for (size_t i = 0; i < num_col; ++i) {
	  total_l += colors[i].cielab.l;
	  total_a += colors[i].cielab.a;
	  total_b += colors[i].cielab.b;
	}
	total_l /= (float)num_col;
	total_a /= (float)num_col;
	total_b /= (float)num_col;
	/* a;lsk djf;alskdjf ;oaksdj fp;oiaeuhdpfiuh  */
	return Color_create(0,0,0);
}

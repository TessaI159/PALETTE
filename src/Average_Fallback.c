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
	}

	l /= num_col;
	a /= num_col;
	b /= num_col;
	return Color_create_cielab(l, a, b);
}

/* NOTE (Tess): This is too slow */
/* TODO (Tess): Potentially change this to use less trig, or faster trig. */
/* struct Color cielab_avg_fallback(struct Color *colors, uint16_t num_col) { */
/* 	float sum_cos = 0.0f, sum_sin = 0.0f, sum_l = 0.0f, sum_c = 0.0f; */
/* 	for (size_t i = 0; i < num_col; ++i) { */
/* 		if (colors[i].current_space != COLOR_CIELAB) { */
/* 			convert_to(&colors[i], COLOR_CIELAB); */
/* 		} */
/* 		float a = colors[i].data.cielab.a; */
/* 		float b = colors[i].data.cielab.b; */
/* 		float h = atan2f(b, a); */
/* 		float c = sqrt(SQUARE(a) + SQUARE(b)); */
/* 		sum_cos += cosf(h); */
/* 		sum_sin += sinf(h); */
/* 		sum_c += c; */
/* 		sum_l += colors[i].data.cielab.l; */
/* 	} */

/* 	float avg_h = atan2f(sum_sin, sum_cos); */
/* 	float avg_c = sum_c / num_col; */
/* 	float l	    = sum_l / num_col; */
/* 	float a	    = avg_c * cosf(avg_h); */
/* 	float b	    = avg_c * sinf(avg_h); */

/* 	return Color_create_cielab(l, a, b); */
/* } */

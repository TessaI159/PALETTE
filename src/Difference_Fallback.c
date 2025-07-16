#include <math.h>
#include <stdlib.h>

#include "Color.h"
#include "Difference_Internal.h"

#define SQUARE(x) ((x) * (x))

/* row * num_colu + colu */
/* Colors are columns, Centroids are rows */
float *delta_ok_diff_fallback(const struct Color colors,
			      const struct Color cents) {

	uint32_t num_col  = colors.num;
	uint32_t num_cent = cents.num;
	float	*diffs	  = malloc(sizeof(float) * num_col * num_cent);

	for (size_t i = 0; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			diffs[j * num_col + i] =
			    sqrtf(SQUARE(colors.alpha[i] - cents.alpha[j]) +
				  SQUARE(colors.beta[i] - cents.beta[j]) +
				  SQUARE(colors.gamma[i] - cents.gamma[j]));
		}
	}
	return diffs;
}

float *delta_cie76_diff_fallback(const struct Color colors,
				 const struct Color cents) {
	uint32_t num_col  = colors.num;
	uint32_t num_cent = cents.num;
	float	*diffs	  = malloc(sizeof(float) * num_col * num_cent);

	for (size_t i = 0; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			diffs[j * num_col + i] =
			    sqrtf(SQUARE(colors.alpha[i] - cents.alpha[j]) +
				  SQUARE(colors.beta[i] - cents.beta[j]) +
				  SQUARE(colors.gamma[i] - cents.gamma[j]));
		}
	}
	return diffs;
}

float *delta_cie94_diff_fallback(const struct Color colors,
				 const struct Color cents) {
	uint32_t	    num_col  = colors.num;
	uint32_t	    num_cent = cents.num;
	float		   *diffs = malloc(sizeof(float) * num_col * num_cent);
	static const double K1	  = 0.045;
	static const double K2	  = 0.015;

	for (size_t i = 0; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			const double dl = colors.alpha[i] - cents.alpha[j];
			const double da = colors.beta[i] - cents.beta[j];
			const double db = colors.gamma[i] - cents.gamma[j];
		}
	}
}

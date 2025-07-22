#include <math.h>
#include <stdlib.h>

#include "Color.h"
#include "Difference_Internal.h"

#define SQUARE(x) ((x) * (x))

/* row * num_col + col */
/* Colors are columns, Centroids are rows */
void delta_ok_diff_fallback(const struct Color colors,
			      const struct Color cents, float *diffs) {
	uint32_t num_col  = colors.num;
	uint32_t num_cent = cents.num;

	for (size_t i = 0; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			diffs[j * num_col + i] =
			    sqrtf(SQUARE(colors.alpha[i] - cents.alpha[j]) +
				  SQUARE(colors.beta[i] - cents.beta[j]) +
				  SQUARE(colors.gamma[i] - cents.gamma[j]));
		}
	}
}

void delta_cie76_diff_fallback(const struct Color colors,
				 const struct Color cents, float *diffs) {
	uint32_t num_col  = colors.num;
	uint32_t num_cent = cents.num;

	for (size_t i = 0; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			diffs[j * num_col + i] =
			    sqrtf(SQUARE(colors.alpha[i] - cents.alpha[j]) +
				  SQUARE(colors.beta[i] - cents.beta[j]) +
				  SQUARE(colors.gamma[i] - cents.gamma[j]));
		}
	}
}

void delta_cie94_diff_fallback(const struct Color colors,
				 const struct Color cents, float *diffs) {
	uint32_t num_col  = colors.num;
	uint32_t num_cent = cents.num;

	static const double K1 = 0.045;
	static const double K2 = 0.015;

	for (size_t i = 0; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			const double dl = colors.alpha[i] - cents.alpha[j];
			const double da = colors.beta[i] - cents.beta[j];
			const double db = colors.gamma[i] - cents.gamma[j];
			const double C1 = sqrt(SQUARE(colors.beta[i]) +
					       SQUARE(colors.gamma[i]));
			const double C2 = sqrt(SQUARE(cents.beta[i]) +
					       SQUARE(cents.gamma[i]));
			const double dC = C1 - C2;

			const double dH2 = SQUARE(da) + SQUARE(db) + SQUARE(dC);

			const double sl = 1.0;
			const double sc = fma(K1, C1, 1.0);
			const double sh = fma(K2, C1, 1.0);

			const double term_L = dl / sl;
			const double term_C = dC / sc;
			const double term_H = sqrt(dH2 < 0.0 ? 0.0 : dH2) / sh;

			diffs[j * num_col + i] = sqrt(
			    SQUARE(term_L) + SQUARE(term_C) + SQUARE(term_H));
		}
	}
}

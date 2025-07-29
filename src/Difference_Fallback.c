#include <math.h>
#include <stdlib.h>

#include "Color.h"
#include "Difference_Internal.h"

#define SQUARE(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
#define POW7(x) ((x) * (x) * (x) * (x) * (x) * (x) * (x))
#define DEG2RAD(x) (((x) / (180.0)) * (M_PI))
#define RAD2DEG(x) (((x) / M_PI) * (180.0))

/* row * num_col + col */
/* Colors are columns, Centroids are rows */
void delta_lab_euc_diff_fallback(const struct Color *restrict colors,
				 const struct Color *restrict cents,
				 float *restrict diffs) {
	const uint32_t num_col	= colors->num;
	const uint32_t num_cent = cents->num;

	for (size_t i = 0; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			diffs[j * num_col + i] =
			    sqrtf(SQUARE(colors->alpha[i] - cents->alpha[j]) +
				  SQUARE(colors->beta[i] - cents->beta[j]) +
				  SQUARE(colors->gamma[i] - cents->gamma[j]));
		}
	}
}

void delta_cie94_diff_fallback(const struct Color *restrict colors,
			       const struct Color *restrict cents,
			       float *restrict diffs) {
	const uint32_t num_col	= colors->num;
	const uint32_t num_cent = cents->num;

	const double K1 = 0.045;
	const double K2 = 0.015;

	for (size_t i = 0; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			const double dl = colors->alpha[i] - cents->alpha[j];
			const double da = colors->beta[i] - cents->beta[j];
			const double db = colors->gamma[i] - cents->gamma[j];
			const double C1 = sqrt(SQUARE(colors->beta[i]) +
					       SQUARE(colors->gamma[i]));
			const double C2 = sqrt(SQUARE(cents->beta[i]) +
					       SQUARE(cents->gamma[i]));
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

void delta_ciede2000_diff_fallback(const struct Color *restrict colors,
				   const struct Color *restrict cents,
				   float *restrict diffs) {
	const uint32_t num_col	= colors->num;
	const uint32_t num_cent = cents->num;

	for (size_t i = 0; i < num_col; ++i) {
		for (size_t j = 0; j < num_cent; ++j) {
			double L1 = colors->alpha[i];
			double a1 = colors->beta[i];
			double b1 = colors->gamma[i];
			double L2 = cents->alpha[i];
			double a2 = cents->beta[i];
			double b2 = cents->gamma[i];

			/* Step 1 */
			const double C1 = sqrt(SQUARE(a1) + SQUARE(b1));
			const double C2 = sqrt(SQUARE(a2) + SQUARE(b2));

			const double avgC = (C1 + C2) / 2.0;
			/* Saves calculating this twice */
			const double G =
			    0.5 * (1.0 - sqrt(POW7(avgC) /
					      (POW7(avgC) + 6103515625.0)));


			const double a1p = (1.0 + G) * a1;
			const double a2p = (1.0 + G) * a2;
			const double C1p = sqrt(SQUARE(a1p) + SQUARE(b1));
			const double C2p = sqrt(SQUARE(a2p) + SQUARE(b2));

			const double at1 = RAD2DEG(atan2(b1, a1p));
			const double h1p = (b1 == 0.0 && a1p == 0.0) ? 0
					   : at1 < 0		     ? at1 + 360
								     : at1;
			const double at2 = RAD2DEG(atan2(b2, a2p));
			const double h2p = (b2 == 0.0 && a2p == 0.0) ? 0
					   : at2 < 0		     ? at2 + 360
								     : at2;

			/* Step 2 */
			const double dLp = L2 - L1;
			const double dCp = C2p - C1p;

			double dhp = 0.0;
			if (C1p == 0.0 || C2p == 0.0) {
				dhp = 0.0;
			} else {
				if (fabs(h2p - h1p) <= 180.0) {
					dhp = h2p - h1p;
				} else if (h2p - h1p > 180.0) {
					dhp = h2p - h1p - 360.0;
				} else if (h2p - h1p < -180.0) {
					dhp = h2p - h1p + 360.0;
				} else {
					diffs[j * num_col + i] = 0.0;
				}
			}

			const double dHp =
			    2.0 * sqrt(C1p * C2p) * sin(DEG2RAD(dhp / 2.0));

			/* Step 3 */
			const double avgLp = (L1 + L2) / 2.0;
			const double avgCp = (C1p + C2p) / 2.0;

			double avghp;
			if (C1p == 0.0 || C2p == 0.0) {
				avghp = h1p + h2p;
			} else {
				if (fabs(h1p - h2p) <= 180.0) {
					avghp = (h1p + h2p) / 2.0;
				} else {
					if (h1p + h2p < 360.0) {
						avghp =
						    (h1p + h2p + 360.0) / 2.0;
					} else {
						avghp =
						    (h1p + h2p - 360.0) / 2.0;
					}
				}
			}

			const double T =
			    1.0 - 0.17 * cos(DEG2RAD(avghp - 30.0)) +
			    0.24 * cos(DEG2RAD(2.0 * avghp)) +
			    0.32 * cos(DEG2RAD(3.0 * avghp + 6.0)) -
			    0.20 * cos(DEG2RAD(4.0 * avghp - 63.0));

			const double dtheta =
			    30.0 * exp(-(SQUARE((avghp - 275.0) / 25.0)));
			const double RC =
			    2.0 *
			    sqrt(POW7(avgCp) / (POW7(avgCp) + 6103515625.0));
			const double SL =
			    1.0 + ((0.015 * (SQUARE(avgLp - 50.0))) /
				   sqrt(20.0 + (SQUARE(avgLp - 50.0))));
			const double SC = 1.0 + 0.045 * avgCp;
			const double SH = 1.0 + 0.015 * avgCp * T;
			const double RT = -sin(DEG2RAD(2.0 * dtheta)) * RC;
			const double terml = dLp / SL;
			const double termc = dCp / SC;
			const double termh = dHp / SH;

			diffs[j * num_col + i] =
			    sqrt(SQUARE(terml) + SQUARE(termc) + SQUARE(termh) +
				 RT * termc * termh);
		}
	}
}

#ifndef DIFFERENCE_INTERNAL_H
#define DIFFERENCE_INTERNAL_H

#include <stdint.h>
typedef struct Color Color;


void delta_lab_euc_diff_fallback(const struct Color *restrict colors,
			    const struct Color *restrict cents,
				  float *restrict diffs);
void delta_lab_euc_diff_avx(const struct Color *restrict colors,
		       const struct Color *restrict cents,
		       float *restrict diffs);
void delta_lab_euc_diff_sse(const struct Color *restrict colors,
			    const struct Color *restrict cents,
			    float *restrict diffs);


void delta_cie94_diff_fallback(const struct Color *restrict colors,
			       const struct Color *restrict cents,
			       float *restrict diffs);
void delta_cie94_diff_sse(const struct Color *restrict colors,
			  const struct Color *restrict cents,
			  float *restrict diffs);
void delta_cie94_diff_avx(const struct Color *restrict colors,
			  const struct Color *restrict cents,
			  float *restrict diffs);

void delta_ciede2000_diff_fallback(const struct Color *restrict colors,
				   const struct Color *restrict cents,
				   float *restrict diffs);
void delta_ciede2000_diff_sse(const struct Color *restrict colors,
			      const struct Color *restrict cents,
			      float *restrict diffs);
void delta_ciede2000_diff_avx(const struct Color *restrict colors,
			      const struct Color *restrict cents,
			      float *restrict diffs);

#endif

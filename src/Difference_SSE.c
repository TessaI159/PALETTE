#if defined(__x86_64__) || defined(_M_X64)
#include <pmmintrin.h>
#endif

#include "Color.h"

/* row * num_col + col */
/* Colors are columns, Centroids are rows */
void delta_ok_diff_sse(const struct Color colors, const struct Color cents,
		       float *diffs) {
	uint16_t num_col  = colors.num;
	uint16_t num_cent = cents.num;
	for (size_t col = 0; col < colors.num; ++col) {
	  
	}
}

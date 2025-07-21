#ifndef DIFFERENCE_INTERNAL_H
#define DIFFERENCE_INTERNAL_H

#include <stdint.h>
typedef struct Color Color;

float *delta_ok_diff_fallback(const struct Color colors,
			      const struct Color centroids, float *diffs);
float *delta_ok_diff_sse(const struct Color colors,
			 const struct Color centroids);
float *delta_ok_diff_avx(const struct Color colors,
			 const struct Color centroids);

float *delta_cie76_diff_fallback(const struct Color colors,
				 const struct Color centroids, float *diffs);
float *delta_cie76_diff_sse(const struct Color colors,
			    const struct Color centroids);
float *delta_cie76_diff_avx(const struct Color colors,
			    const struct Color centroids);

float *delta_cie94_diff_fallback(const struct Color colors,
				 const struct Color centroids, float *diffs);
float *delta_cie94_diff_sse(const struct Color colors,
			    const struct Color centroids);
float *delta_cie94_diff_avx(const struct Color colors,
			    const struct Color centroids);

float *delta_ciede2000_diff_fallback(const struct Color colors,
				     const struct Color centroids, float *diffs);
float *delta_ciede2000_diff_sse(const struct Color colors,
				const struct Color centroids);
float *delta_ciede2000_diff_avx(const struct Color colors,
				const struct Color centroids);

#endif

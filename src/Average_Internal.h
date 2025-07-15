#ifndef AVERAGE_INTERNAL_H
#define AVERAGE_INTERNAL_H

#include <stdint.h>

typedef struct cielab_SoA cielab_SoA;
typedef struct oklab_SoA  oklab_SoA;

/* None of these are static because the functions utilizing SSE/AVX have to be
   in their own translation units to avoid blowing up a computer that doesn't
   support them */

struct Color cielab_avg_avx(const struct cielab_SoA *__restrict colors,
			    uint16_t num_col);
struct Color cielab_avg_sse(const struct cielab_SoA *__restrict colors,
			    uint16_t num_col);
struct Color cielab_avg_fallback(const struct cielab_SoA *__restrict colors,
				 uint16_t num_col);
struct Color cielab_avg_avx_cw(const struct cielab_SoA *__restrict colors,
			       uint16_t num_col);
struct Color cielab_avg_sse_cw(const struct cielab_SoA *__restrict colors,
			       uint16_t num_col);
struct Color cielab_avg_fallback_cw(const struct cielab_SoA *__restrict colors,
				    uint16_t num_col);

struct Color oklab_avg_avx(const struct oklab_SoA *__restrict colors,
			   uint16_t num_col);
struct Color oklab_avg_sse(const struct oklab_SoA *__restrict colors,
			   uint16_t num_col);
struct Color oklab_avg_fallback(const struct oklab_SoA *__restrict colors,
				uint16_t num_col);
struct Color oklab_avg_avx_cw(const struct oklab_SoA *__restrict colors,
			      uint16_t num_col);
struct Color oklab_avg_sse_cw(const struct oklab_SoA *__restrict colors,
			      uint16_t num_col);
struct Color oklab_avg_fallback_cw(const struct oklab_SoA *__restrict colors,
				   uint16_t num_col);

#endif

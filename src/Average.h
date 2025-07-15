#ifndef AVERAGE_H
#define AVERAGE_H

#include <stdint.h>

void			  average_init();
typedef struct cielab_SoA cielab_SoA;
typedef struct oklab_SoA  oklab_SoA;

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
struct Color cielab_avg(const struct cielab_SoA *colors, uint16_t num_col);

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
struct Color oklab_avg(const struct oklab_SoA *colors, uint16_t num_col);

#endif

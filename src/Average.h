#ifndef AVERAGE_H
#define AVERAGE_H

#include <stdint.h>

void average_init();

struct Color cielab_avg_avx2(struct Color *colors, uint16_t num_col);
struct Color cielab_avg_avx(struct Color *colors, uint16_t num_col);
struct Color cielab_avg_sse2(struct Color *colors, uint16_t num_col);
struct Color cielab_avg_fallback(struct Color *__restrict colors, uint16_t num_col);

struct Color cielab_avg(struct Color *colors, uint16_t num_col);

#endif

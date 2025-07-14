#ifndef AVERAGE_H
#define AVERAGE_H

#include <stdint.h>

void average_init();
typedef struct cielab_SoA cielabSoA;

struct Color cielab_avg_avx2(const struct cielab_SoA *colors, uint16_t num_col);
struct Color cielab_avg_avx(const struct cielab_SoA *colors, uint16_t num_col);
struct Color cielab_avg_sse3(const struct cielab_SoA *colors, uint16_t num_col);
struct Color cielab_avg_fallback(const struct cielab_SoA *__restrict colors, uint16_t num_col);

struct Color cielab_avg(const struct cielab_SoA *colors, uint16_t num_col);

#endif

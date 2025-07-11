#ifndef AVERAGE_H
#define AVERAGE_H

#include <stdint.h>

void average_init();
struct Color lab_avg_fallback(const struct Color **colors, uint16_t num_col);
struct Color lab_avg_avx(const struct Color **colors, uint16_t num_col);
struct Color lab_avg_avx2(const struct Color **colors, uint16_t num_col);
struct Color lab_avg(const struct Color **colors, uint16_t num_col);

#endif

#ifndef UTIL_H
#define UTIL_H

#include <stdint.h>
typedef struct Color	Color;
typedef enum ColorSpace ColorSpace;

uint64_t time_diff(float (*func)(Color *, Color *), Color col1, Color col2,
		   const uint64_t RUNS);
uint64_t time_conv(void (*func)(Color *, enum ColorSpace), Color col,
		   enum ColorSpace f, enum ColorSpace t, const uint64_t RUNS);
uint64_t time_avg(Color (*func)(Color *, uint16_t), Color *colors,
		  uint16_t num_col, uint64_t RUNS);

#endif

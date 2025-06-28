#ifndef UTIL_H
#define UTIL_H
#include <stddef.h>
#include <stdint.h>
#include <x86intrin.h>

typedef struct Color Color;


uint_fast64_t measure_cycles(double (*func)(Color *, Color *), Color *,
			     Color *);
uint_fast64_t test_diff_speed(double (*func)(Color *, Color *),
			      const uint64_t RUNS);

void *aligned_malloc(uint64_t alignment, uint64_t size);
void aligned_free(void* ptr);

#endif

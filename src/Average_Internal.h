#ifndef AVERAGE_INTERNAL_H
#define AVERAGE_INTERNAL_H

#include <stdint.h>

typedef struct Color Color;

/* None of these are static because the functions utilizing SSE/AVX have to be
   in their own translation units to avoid blowing up a computer that doesn't
   support them */

void cielab_avg_avx(const struct Color *restrict colors,
		    struct Color *restrict cents, uint8_t which_cent);
void cielab_avg_sse(const struct Color *restrict colors,
		    struct Color *restrict cents, uint8_t which_cent);
void cielab_avg_fb(const struct Color *restrict colors,
		   struct Color *restrict cents, uint8_t which_cent);
void cielab_avg_avx_cw(const struct Color *restrict colors,
		       struct Color *restrict cents, uint8_t which_cent);
void cielab_avg_sse_cw(const struct Color *restrict colors,
		       struct Color *restrict cents, uint8_t which_cent);
void cielab_avg_fb_cw(const struct Color *restrict colors,
		      struct Color *restrict cents, uint8_t which_cent);

void oklab_avg_avx(const struct Color *restrict colors,
		   struct Color *restrict cents, uint8_t which_cent);
void oklab_avg_sse(const struct Color *restrict colors,
		   struct Color *restrict cents, uint8_t which_cent);
void oklab_avg_fb(const struct Color *restrict colors,
		  struct Color *restrict cents, uint8_t which_cent);
void oklab_avg_avx_cw(const struct Color *restrict colors,
		      struct Color *restrict cents, uint8_t which_cent);
void oklab_avg_sse_cw(const struct Color *restrict colors,
		      struct Color *restrict cents, uint8_t which_cent);
void oklab_avg_fb_cw(const struct Color *restrict colors,
		     struct Color *restrict cents, uint8_t which_cent);

#endif

#ifndef CONVERT_H
#define CONVERT_H

typedef struct Color Color;

void convert_cielab_to_srgb(struct Color colors);
void convert_oklab_to_srgb(struct Color colors);
void convert_srgb_to_cielab(struct Color colors);
void convert_srgb_to_oklab(struct Color colors);

#endif

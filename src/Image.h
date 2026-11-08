#ifndef IMAGE_H
#define IMAGE_H

#include "Color.h"

#include <stdbool.h>
#include <stdint.h>

struct Image {
  uint32_t w;
  uint32_t h;
  uint8_t *rgb; /* Interleaved RGB24 buffer: w * h * 3 bytes */
};

/**
 * @brief Frees the internal buffer of an Image.
 */
void Image_free(struct Image *img);

/**
 * @brief Loads an image from a PPM (P6 or P3) or 24-bit BMP file.
 */
bool Image_load(struct Image *img, const char *filename);

/**
 * @brief Saves an Image to a binary PPM (P6) file
 */
bool Image_save_ppm(const struct Image *img, const char *filename);

/**
 * @brief Ingests an Image into an SoA struct Color with optional downsampling.
 * @param img Source Image
 * @param max_samples Maximum number of pixels to ingest (e.g. 10000). If 0 or
 * >= w*h, ingests all.
 * @param out_colors Allocated destination Color SoA with linear sRGB values [0,
 * 1].
 * @return true on success, false on failure.
 */
bool Image_to_color_soa(const struct Image *img, uint32_t max_samples,
                        struct Color *out_colors);

/**
 * @brief Saves an extracted palette to a PPM swatch visualizer
 * @param filename Output .ppm file path.
 * @param palette Extracted centroids in sRGB space.
 * @param swatch_width Width of each color swatch.
 * @param swatch_height Height of each color swatch.
 */
bool Image_save_palette_swatch_ppm(const char *filename,
                                   const struct Color *restrict palette,
                                   uint32_t swatch_width,
                                   uint32_t swatch_height);
#endif

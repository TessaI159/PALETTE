#include "Image.h"
#include "Conversion.h"

#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool Image_alloc(struct Image *img, uint32_t w, uint32_t h) {
  if (!img)
    return false;
  img->w = 0;
  img->h = 0;
  img->rgb = NULL;

  if (w == 0 || h == 0)
    return true;

  /* Guard against size_t overflow in w * h * 3 */
  if ((size_t)w > SIZE_MAX / 3 / h)
    return false;

  size_t size = (size_t)w * h * 3;
  img->rgb = (uint8_t *)malloc(size);
  if (!img->rgb)
    return false;

  img->w = w;
  img->h = h;
  return true;
}

void Image_free(struct Image *img) {
  if (!img)
    return;
  if (img->rgb) {
    free(img->rgb);
    img->rgb = NULL;
  }
  img->w = 0;
  img->h = 0;
}

/* PPM loading routine */
static inline void skip_ppm_comments(FILE *fp) {
  int ch;
  while ((ch = fgetc(fp)) != EOF) {
    if (isspace(ch))
      continue;
    if (ch == '#') {
      /* Skip till the end of the line */
      while ((ch = fgetc(fp)) != EOF && ch != '\n')
        ;
    } else {
      ungetc(ch, fp);
      break;
    }
  }
}

/* Rescale a sample in [0, maxval] to [0, 255], with rounding. */
static inline uint8_t scale_sample(uint32_t v, uint32_t maxval) {
  if (v > maxval)
    v = maxval;
  return (uint8_t)((v * 255u + maxval / 2) / maxval);
}

static bool load_ppm(struct Image *img, FILE *fp, int format) {
  if (format != 3 && format != 6)
    return false;

  uint32_t w = 0, h = 0, maxval = 0;
  skip_ppm_comments(fp);
  if (fscanf(fp, "%" SCNu32, &w) != 1)
    return false;
  skip_ppm_comments(fp);
  if (fscanf(fp, "%" SCNu32, &h) != 1)
    return false;
  skip_ppm_comments(fp);
  if (fscanf(fp, "%" SCNu32, &maxval) != 1)
    return false;

  if (fgetc(fp) == EOF) /* Consume single trailing whitespace */
    return false;
  if (w == 0 || h == 0 || maxval == 0 || maxval > 65535)
    return false;
  if (!Image_alloc(img, w, h))
    return false;

  size_t total_pixels = (size_t)w * h;
  size_t samples = total_pixels * 3;

  if (format == 6) {
    /* P6: binary. maxval > 255 means 2 bytes per sample, big-endian. */
    size_t bps = maxval > 255 ? 2 : 1;

    if (bps == 1 && maxval == 255) {
      /* Fast path: bytes are already 8-bit RGB */
      if (fread(img->rgb, 1, samples, fp) != samples) {
        Image_free(img);
        return false;
      }
      return true;
    }

    if (samples > SIZE_MAX / bps) {
      Image_free(img);
      return false;
    }
    size_t   raw_size = samples * bps;
    uint8_t *raw = (uint8_t *)malloc(raw_size);
    if (!raw) {
      Image_free(img);
      return false;
    }
    if (fread(raw, 1, raw_size, fp) != raw_size) {
      free(raw);
      Image_free(img);
      return false;
    }
    for (size_t i = 0; i < samples; ++i) {
      uint32_t v =
          (bps == 1) ? raw[i] : (((uint32_t)raw[2 * i] << 8) | raw[2 * i + 1]);
      img->rgb[i] = scale_sample(v, maxval);
    }
    free(raw);
  } else {
    /* P3: ASCII */
    for (size_t i = 0; i < samples; ++i) {
      uint32_t val;
      if (fscanf(fp, "%" SCNu32, &val) != 1) {
        Image_free(img);
        return false;
      }
      img->rgb[i] = scale_sample(val, maxval);
    }
  }
  return true;
}

/* BMP loading routine */
/* https://en.wikipedia.org/wiki/BMP_file_format#File_header */
/* Note: BMP fields are little-endian; this reads them directly, so it assumes
   a little-endian host. */
#pragma pack(push, 1)
struct bmp_file_header {
  uint16_t bfType;
  uint32_t bfSize;
  uint16_t bfReserved1;
  uint16_t bfReserved2;
  uint32_t bfOffBits;
};

struct bmp_info_header {
  uint32_t biSize;
  int32_t  biWidth;
  int32_t  biHeight;
  uint16_t biPlanes;
  uint16_t biBitCount;
  uint32_t biCompression;
  uint32_t biSizeImage;
  int32_t  biXPelsPerMeter;
  int32_t  biYPelsPerMeter;
  uint32_t biClrUsed;
  uint32_t biClrImportant;
};
#pragma pack(pop)

static bool load_bmp(struct Image *img, FILE *fp) {
  struct bmp_file_header bfh;
  struct bmp_info_header bih;

  if (fread(&bfh, sizeof(bfh), 1, fp) != 1)
    return false;
  if (bfh.bfType != 0x4D42) /* 4D42 is the magic number indicating bitmap */
    /* https://jacobfilipp.com/DrDobbs/articles/DDJ/1994/9409/9409a/9409a.htm#00c2_0049
     */
    return false;

  if (fread(&bih, sizeof(bih), 1, fp) != 1)
    return false;
  if (bih.biSize < sizeof(bih))
    return false;
  if (bih.biBitCount != 24 && bih.biBitCount != 32)
    return false;
  if (bih.biCompression != 0) /* No compressed images */
    return false;

  /* Width must be positive; height may be negative (top-down), but
     INT32_MIN can't be negated safely. */
  if (bih.biWidth <= 0 || bih.biHeight == 0 || bih.biHeight == INT32_MIN)
    return false;

  uint32_t w = (uint32_t)bih.biWidth;
  uint32_t h = (uint32_t)(bih.biHeight > 0 ? bih.biHeight : -bih.biHeight);
  bool     top_down = bih.biHeight < 0;

  uint32_t bytes_per_pixel = bih.biBitCount / 8;
  if ((size_t)w > (SIZE_MAX - 3) / bytes_per_pixel)
    return false;
  size_t row_stride = (((size_t)w * bytes_per_pixel + 3) / 4) * 4;

  if (!Image_alloc(img, w, h))
    return false;

  if (fseek(fp, (long)bfh.bfOffBits, SEEK_SET) != 0) {
    Image_free(img);
    return false;
  }

  uint8_t *row_buf = (uint8_t *)malloc(row_stride);
  if (!row_buf) {
    Image_free(img);
    return false;
  }

  for (uint32_t y = 0; y < h; ++y) {
    uint32_t target_y = top_down ? y : (h - 1 - y);
    if (fread(row_buf, 1, row_stride, fp) != row_stride) {
      free(row_buf);
      Image_free(img);
      return false;
    }

    /* row_buf and img->rgb are separate allocations, so the accesses
       through src and dest never overlap. */
    const uint8_t *restrict src = row_buf;
    uint8_t *restrict dest = img->rgb + (size_t)target_y * w * 3;
    for (uint32_t x = 0; x < w; ++x) {
      /* BMP stores B, G, R (and optionally A); we want R, G, B */
      dest[x * 3 + 0] = src[x * bytes_per_pixel + 2];
      dest[x * 3 + 1] = src[x * bytes_per_pixel + 1];
      dest[x * 3 + 2] = src[x * bytes_per_pixel + 0];
    }
  }
  free(row_buf);
  return true;
}

bool Image_load(struct Image *img, const char *filename) {
  if (!img || !filename)
    return false;

  /* Leave img in a valid empty state on any failure. */
  img->w = 0;
  img->h = 0;
  img->rgb = NULL;

  FILE *fp = fopen(filename, "rb");
  if (!fp)
    return false;

  char magic[2];
  if (fread(magic, 1, 2, fp) != 2) {
    fclose(fp);
    return false;
  }

  bool ok = false;
  if (magic[0] == 'P' && (magic[1] == '6' || magic[1] == '3')) {
    /* Already positioned just past the 2-byte magic */
    ok = load_ppm(img, fp, magic[1] - '0');
  } else if (magic[0] == 'B' && magic[1] == 'M') {
    rewind(fp); /* BMP header parsing starts at byte 0 */
    ok = load_bmp(img, fp);
  }

  fclose(fp);
  return ok;
}

bool Image_save_ppm(const struct Image *img, const char *filename) {
  if (!img || !img->rgb || img->w == 0 || img->h == 0 || !filename) {
    return false;
  }
  FILE *fp = fopen(filename, "wb");
  if (!fp)
    return false;

  size_t total_pixels = (size_t)img->w * img->h;
  bool   ok =
      fprintf(fp, "P6\n%" PRIu32 " %" PRIu32 "\n255\n", img->w, img->h) > 0;
  if (ok)
    ok = fwrite(img->rgb, 3, total_pixels, fp) == total_pixels;
  if (fclose(fp) != 0)
    ok = false;
  return ok;
}

bool Image_to_color_soa(const struct Image *img, uint32_t max_samples,
                        struct Color *out_colors) {
  if (!img || !img->rgb || img->w == 0 || img->h == 0 || !out_colors)
    return false;

  size_t   total_pixels = (size_t)img->w * img->h;
  uint32_t n;
  if (max_samples > 0 && max_samples < total_pixels)
    n = max_samples;
  else
    n = total_pixels > UINT32_MAX ? UINT32_MAX : (uint32_t)total_pixels;

  if (!Color_alloc(out_colors, n))
    return false;

  /* Local restrict pointer: lets the compiler know nothing else writes to
     the pixel data while we read it. */
  const uint8_t *restrict rgb = img->rgb;

  double step = (double)total_pixels / (double)n;
  for (uint32_t i = 0; i < n; ++i) {
    size_t pixel_idx = (size_t)(i * step);
    if (pixel_idx >= total_pixels)
      pixel_idx = total_pixels - 1;

    float r = (float)rgb[pixel_idx * 3 + 0] / 255.0f;
    float g = (float)rgb[pixel_idx * 3 + 1] / 255.0f;
    float b = (float)rgb[pixel_idx * 3 + 2] / 255.0f;

    Color_set(out_colors, i, r, g, b);
  }

  /* Linearize */
  srgb_linearize(out_colors);
  return true;
}
bool Image_save_palette_swatch_ppm(const char         *filename,
                                   const struct Color *palette,
                                   uint32_t            swatch_width,
                                   uint32_t            swatch_height) {
  if (!filename || !palette || palette->num == 0)
    return false;
  if (swatch_width == 0)
    swatch_width = 80;
  if (swatch_height == 0)
    swatch_height = 80;

  uint32_t k = palette->num;
  uint32_t total_w = k * swatch_width;
  uint32_t total_h = swatch_height;

  struct Image img;

  if (!Image_alloc(&img, total_w, total_h))
    return false;

  /* Clone palette to delinearize without mutating input */
  struct Color srgb_pal;
  if (!Color_clone(&srgb_pal, palette)) {
    Image_free(&img);
    return false;
  }
  srgb_gamma_correct(&srgb_pal);

  for (uint32_t j = 0; j < k; ++j) {
    int r = (int)roundf(fminf(fmaxf(srgb_pal.alpha[j], 0.0f), 1.0f) * 255.0f);
    int g = (int)roundf(fminf(fmaxf(srgb_pal.beta[j], 0.0f), 1.0f) * 255.0f);
    int b = (int)roundf(fminf(fmaxf(srgb_pal.gamma[j], 0.0f), 1.0f) * 255.0f);

    for (uint32_t y = 0; y < total_h; ++y) {
      for (uint32_t x = j * swatch_width; x < (j + 1) * swatch_width; ++x) {
        size_t idx = (y * total_w + x) * 3;
        img.rgb[idx + 0] = (uint8_t)r;
        img.rgb[idx + 1] = (uint8_t)g;
        img.rgb[idx + 2] = (uint8_t)b;
      }
    }
  }

  Color_free(&srgb_pal);
  bool ok = Image_save_ppm(&img, filename);
  Image_free(&img);
  return ok;
}

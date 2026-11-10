#ifndef CONVERSION_H
#define CONVERSION_H

typedef struct Color Color;

/* Conversions between linear normalized sRGB [0,1] and perceptual spaces */

/**
 * @brief Converts linear, normalized sRGB colors to oklab.
 * @param colors Input color array in SoA format. Modified in place.
 */
void srgb_to_oklab(struct Color *restrict colors);
/**
 * @brief Converts linear, normalized sRGB colors to CIELAB.
 * @param colors Input color array in SoA format. Modified in place.
 */
void srgb_to_cielab(struct Color *restrict colors);
/**
 * @brief Converts OKLab colors to linear, normalized sRGB.
 * @param colors Input color array in SoA format. Modified in place.
 */
void oklab_to_srgb(struct Color *restrict colors);
/**
 * @brief Converts CIELAB colors to linear, normalized sRGB.
 * @param colors Input color array in SoA format. Modified in place.
 */
void cielab_to_srgb(struct Color *restrict colors);

//* Gamma encode and decode. [0, 1] */
/**
 * @brief linearizes (gamma decodes) gamma-corrected, normalized sRGB colors.
 * @param colors Input color array in SoA format. Modified in place.
 */
void srgb_linearize(struct Color *restrict colors);
/**
 * @brief Gamma corrects (delinearizes) lineariazed, normalized sRGB colors.
 * @param colors Input color array in SoA format. Modified in place.
 */
void srgb_gamma_correct(struct Color *restrict colors);

#endif

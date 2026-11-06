#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Algorithms.h"
#include "Average.h"
#include "Color.h"
#include "Conversion.h"
#include "Difference.h"
#include "Parameters.h"
#include "unity.h"

#define TOL 1e-3f
#define COLOR_TOL 1e-6f

void test_color_alloc_and_lifecycle(void) {
  struct Color c;

  /* Allocation */
  TEST_ASSERT_TRUE(Color_alloc(&c, 100));
  TEST_ASSERT_EQUAL_UINT(100, c.num);
  TEST_ASSERT_NOT_NULL(c.alpha);
  TEST_ASSERT_NOT_NULL(c.beta);
  TEST_ASSERT_NOT_NULL(c.gamma);

  /* 32-byte alignment */
  TEST_ASSERT_EQUAL_PTR(0, (uintptr_t)c.alpha % 32);
  TEST_ASSERT_EQUAL_PTR(0, (uintptr_t)c.beta % 32);
  TEST_ASSERT_EQUAL_PTR(0, (uintptr_t)c.gamma % 32);

  /* Set and get */
  Color_set(&c, 0, 0.1f, 0.2f, 0.3f);
  Color_set(&c, 99, 0.7f, 0.8f, 0.9f);

  float a, b, g;
  Color_get(&c, 0, &a, &b, &g);
  TEST_ASSERT_FLOAT_WITHIN(COLOR_TOL, 0.1f, a);
  TEST_ASSERT_FLOAT_WITHIN(COLOR_TOL, 0.2f, b);
  TEST_ASSERT_FLOAT_WITHIN(COLOR_TOL, 0.3f, g);

  Color_get(&c, 99, &a, &b, &g);
  TEST_ASSERT_FLOAT_WITHIN(COLOR_TOL, 0.7f, a);
  TEST_ASSERT_FLOAT_WITHIN(COLOR_TOL, 0.8f, b);
  TEST_ASSERT_FLOAT_WITHIN(COLOR_TOL, 0.9f, g);

  /* Clone */
  struct Color clone;
  TEST_ASSERT_TRUE(Color_clone(&clone, &c));
  TEST_ASSERT_EQUAL_UINT(100, clone.num);
  Color_get(&clone, 99, &a, &b, &b);
  TEST_ASSERT_FLOAT_WITHIN(COLOR_TOL, 0.7f, a);

  /* Change original, not the clone */
  Color_set(&c, 99, 0.0f, 0.0f, 0.0f);
  Color_get(&clone, 99, &a, &b, &b);
  TEST_ASSERT_FLOAT_WITHIN(COLOR_TOL, 0.7f, a);

  /* Free */
  Color_free(&c);
  Color_free(&clone);
  TEST_ASSERT_NULL(c.alpha);
  TEST_ASSERT_EQUAL_UINT(0, c.num);
}

void test_srgb_oklab_roundtrip(void) {
  /* Reference test colors in linear sRGB [0, 1] */
  const uint32_t num = 8;
  float test_rgb[8][3] = {
      {1.0f, 1.0f, 1.0f}, /* White */
      {0.0f, 0.0f, 0.0f}, /* Black */
      {1.0f, 0.0f, 0.0f}, /* Red */
      {0.0f, 1.0f, 0.0f}, /* Green */
      {0.0f, 0.0f, 1.0f}, /* Blue */
      {0.5f, 0.5f, 0.5f}, /* Mid Gray */
      {0.8f, 0.4f, 0.2f}, /* Orange/Brown */
      {0.2f, 0.7f, 0.9f}  /* Cyan/Sky */
  };
  struct Color c;
  Color_alloc(&c, num);
  for (uint32_t i = 0; i < num; ++i) {
    Color_set(&c, i, test_rgb[i][0], test_rgb[i][1], test_rgb[i][2]);
  }
  /* sRGB -> Oklab */
  srgb_to_oklab(&c);

  /* Verify white: L ~ 1.0, a ~ 0.0, b ~ 0.0 */
  TEST_ASSERT_FLOAT_WITHIN(0.01f, 1.0f, c.alpha[0]);
  TEST_ASSERT_FLOAT_WITHIN(0.01f, 0.0f, c.beta[0]);
  TEST_ASSERT_FLOAT_WITHIN(0.01f, 0.0f, c.gamma[0]);

  /* Verify black: L ~ 0.0, a ~ 0.0, b ~ 0.0 */
  TEST_ASSERT_FLOAT_WITHIN(0.01f, 0.0f, c.alpha[1]);
  TEST_ASSERT_FLOAT_WITHIN(0.01f, 0.0f, c.beta[1]);
  TEST_ASSERT_FLOAT_WITHIN(0.01f, 0.0f, c.gamma[1]);

  /* Oklab -> sRGB */
  oklab_to_srgb(&c);
  for (uint32_t i = 0; i < num; ++i) {
    TEST_ASSERT_FLOAT_WITHIN(TOL, test_rgb[i][0], c.alpha[i]);
    TEST_ASSERT_FLOAT_WITHIN(TOL, test_rgb[i][1], c.beta[i]);
    TEST_ASSERT_FLOAT_WITHIN(TOL, test_rgb[i][2], c.gamma[i]);
  }
  Color_free(&c);
}

void test_srgb_cielab_roundtrip(void) {

  /* Reference test colors in linear sRGB [0, 1] */
  const uint32_t num = 8;
  float test_rgb[8][3] = {
      {1.0f, 1.0f, 1.0f}, /* White */
      {0.0f, 0.0f, 0.0f}, /* Black */
      {1.0f, 0.0f, 0.0f}, /* Red */
      {0.0f, 1.0f, 0.0f}, /* Green */
      {0.0f, 0.0f, 1.0f}, /* Blue */
      {0.5f, 0.5f, 0.5f}, /* Mid Gray */
      {0.8f, 0.4f, 0.2f}, /* Orange/Brown */
      {0.2f, 0.7f, 0.9f}  /* Cyan/Sky */
  };

  struct Color c;
  Color_alloc(&c, num);
  for (uint32_t i = 0; i < num; ++i) {
    Color_set(&c, i, test_rgb[i][0], test_rgb[i][1], test_rgb[i][2]);
  }

  /* sRGB -> CIElab */

  srgb_to_cielab(&c);

  /* Verify white */
  TEST_ASSERT_FLOAT_WITHIN(0.5f, 100.0f, c.alpha[0]);
  TEST_ASSERT_FLOAT_WITHIN(0.5f, 0.0f, c.beta[0]);
  TEST_ASSERT_FLOAT_WITHIN(0.5f, 0.0f, c.gamma[0]);

  /* verify black */
  TEST_ASSERT_FLOAT_WITHIN(0.5f, 0.0f, c.alpha[0]);
  TEST_ASSERT_FLOAT_WITHIN(0.5f, 0.0f, c.beta[0]);
  TEST_ASSERT_FLOAT_WITHIN(0.5f, 0.0f, c.gamma[0]);

  /* CIElab -> sRGB */
  cielab_to_srgb(&c);

  for (uint32_t i = 0; i < num; ++i) {
    TEST_ASSERT_FLOAT_WITHIN(TOL, test_rgb[i][0], c.alpha[i]);
    TEST_ASSERT_FLOAT_WITHIN(TOL, test_rgb[i][1], c.beta[i]);
    TEST_ASSERT_FLOAT_WITHIN(TOL, test_rgb[i][2], c.gamma[i]);
  }
  Color_free(&c);
}

/* TODO (Tess): Use Sharma here for CIEDE2000*/

void test_difference_metrics(void) {
  struct Color colors;
  struct Color cents;
  Color_alloc(&colors, 2);
  Color_alloc(&cents, 2);

  /* Point 1: (50, 10, 20), Point 2: (60, -10, 30) */
  Color_set(&colors, 0, 50.0f, 10.0f, 20.0f);
  Color_set(&colors, 2, 60.0f, -10.0f, 30.0f);

  Color_set(&cents, 0, 50.0f, 10.0f, 20.0f);
  Color_set(&cents, 2, 60.0f, -10.0f, 30.0f);

  float diffs[4];

  g_params.space = OK;
  g_params.diff_fn = EUCLIDEAN;
  delta_diff(&colors, &cents, diffs);

  /* Check metric identity d(0,0)=d(1,1)=0 */
  TEST_ASSERT_FLOAT_WITHIN(1e-5f, 0.0f, diffs[0 * 2 + 0]);
  TEST_ASSERT_FLOAT_WITHIN(1e-5f, 0.0f, diffs[1 * 2 + 1]);

  /* Check metric symmetry d(0,1)=d(1,0) */
  TEST_ASSERT_FLOAT_WITHIN(1e-5f, diffs[0 * 2 + 1], diffs[1 * 2 + 0]);

  /* Expected Euclidean ! 24.4949 */
  TEST_ASSERT_FLOAT_WITHIN(1e-3f, 24.4949f, diffs[0 * 2 + 1]);

  /* Test CIE94 distance */
  g_params.space = CIE;
  g_params.diff_fn = D94;
  delta_diff(&colors, &cents, diffs);

  TEST_ASSERT_FLOAT_WITHIN(1e-5f, 0.0f, diffs[0 * 2 + 0]);
  TEST_ASSERT_FLOAT_WITHIN(1e-5f, 0.0f, diffs[1 * 2 + 1]);
  TEST_ASSERT_TRUE(diffs[0 * 2 + 1 > 0.0f]);

  /* Test CIEDE2000 distance */
  g_params.space = CIE;
  g_params.diff_fn = D2000;
  delta_diff(&colors, &cents, diffs);

  TEST_ASSERT_FLOAT_WITHIN(1e-4f, 0.0f, diffs[0 * 2 + 0]);
  TEST_ASSERT_FLOAT_WITHIN(1e-4f, 0.0f, diffs[1 * 2 + 1]);
  TEST_ASSERT_TRUE(diffs[0 * 2 + 1] > 0.0f);

  Color_free(&colors);
  Color_free(&cents);
}

void test_color_average(void) {
  struct Color colors;
  struct Color cents;
  Color_alloc(&colors, 5);
  Color_alloc(&cents, 1);

  /* 5 Samples with a known mean */
  Color_set(&colors, 0, 10.0f, 20.0f, 30.0f);
  Color_set(&colors, 1, 20.0f, 30.0f, 40.0f);
  Color_set(&colors, 2, 30.0f, 40.0f, 50.0f);
  Color_set(&colors, 3, 40.0f, 50.0f, 60.0f);
  Color_set(&colors, 4, 50.0f, 60.0f, 70.0f);

  g_params.space = CIE;
  g_params.cw = false;
  color_avg(&colors, &cents, 0);

  /* Mean should be 30, 40, 50 under Euclidean average */
  TEST_ASSERT_FLOAT_WITHIN(1e-4f, 30.0f, cents.alpha[0]);
  TEST_ASSERT_FLOAT_WITHIN(1e-4f, 40.0f, cents.beta[0]);
  TEST_ASSERT_FLOAT_WITHIN(1e-4f, 50.0f, cents.gamma[0]);

  Color_free(&colors);
  Color_free(&cents);
}

void test_kmeans_clustering(void) {
  /* Generate 300 points in 3 tight distinct clusters */

  const uint32_t num = 300;
  struct Color data;
  Color_alloc(&data, num);

  float centers[3][3] = {
      {10.0f, 50.0f, -20.0f}, {70.0f, -30.0f, 40.0f}, {40.0f, 10.0f, 80.0f}};

  srand(1234);
  for (uint32_t i = 0; i < num; ++i) {
    int c = i % 3;
    float noise_a = ((float)(rand() % 20) - 10.0f) / 10.0f;
    float noise_b = ((float)(rand() % 20) - 10.0f) / 10.0f;
    float noise_g = ((float)(rand() % 20) - 10.0f) / 10.0f;
    Color_set(&data, i, centers[c][0] + noise_a, centers[c][1] + noise_b,
              centers[c][2] + noise_g);
  }

  g_params.space = OK;
  g_params.diff_fn = EUCLIDEAN;
  g_params.cw = false;

  struct kmeans_params kp = {
      .k = 3, .max_iters = 30, .tolerance = 1e-4f, .seed = 999};

  struct kmeans_result res;
  TEST_ASSERT_TRUE(kmeans_run(&data, &kp, &res));
  TEST_ASSERT_EQUAL_UINT(3, res.centroids.num);
  TEST_ASSERT_TRUE(res.iters_run > 0);
  TEST_ASSERT_TRUE(res.inertia < 500.0f);

  kmeans_result_free(&res);
  Color_free(&data);
}

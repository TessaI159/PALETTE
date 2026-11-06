#include "unity.h"

extern void test_queue_setup(void);
extern void test_queue_push(void);
extern void test_queue_pop(void);
extern void test_queue_counters(void);

/* Color & Algorithm Tests */
extern void test_color_alloc_and_lifecycle(void);
extern void test_srgb_oklab_roundtrip(void);
extern void test_srgb_cielab_roundtrip(void);
extern void test_color_average(void);
extern void test_kmeans_clustering(void);
extern void test_difference_metrics(void);

void setUp(void) {};
void tearDown(void) {};

int main(void) {
  UNITY_BEGIN();

  /* Queue unit tests */
  RUN_TEST(test_queue_setup);
  RUN_TEST(test_queue_push);
  RUN_TEST(test_queue_pop);
  RUN_TEST(test_queue_counters);

  /* Color & SoA Tests */
  RUN_TEST(test_color_alloc_and_lifecycle);
  RUN_TEST(test_srgb_oklab_roundtrip);
  RUN_TEST(test_srgb_cielab_roundtrip);
  RUN_TEST(test_difference_metrics);
  RUN_TEST(test_color_average);
  RUN_TEST(test_kmeans_clustering);

  return UNITY_END();
}

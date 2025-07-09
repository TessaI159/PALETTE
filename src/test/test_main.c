#include "unity.h"

extern void test_queue_setup(void);
extern void test_queue_push(void);
extern void test_queue_pop(void);
extern void test_queue_counters(void);
extern void test_color_create(void);
extern void test_cie94_diff(void);
extern void test_cie76_diff(void);
extern void test_oklab_diff(void);
extern void test_ciede2000_diff(void);
extern void test_color_check_flags(void);

int main(void) {
	UNITY_BEGIN();
	RUN_TEST(test_queue_setup);
	RUN_TEST(test_queue_push);
	RUN_TEST(test_queue_pop);
	RUN_TEST(test_queue_counters);
	RUN_TEST(test_color_check_flags);
	RUN_TEST(test_color_create);
	RUN_TEST(test_cie76_diff);
	RUN_TEST(test_cie94_diff);
	RUN_TEST(test_ciede2000_diff);
	RUN_TEST(test_oklab_diff);
	return UNITY_END();
}

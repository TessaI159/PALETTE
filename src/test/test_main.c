#include "unity.h"

void setUp(void) {
}
void tearDown(void) {
}

extern void test_queue_setup(void);
extern void test_queue_push(void);
extern void test_queue_pop(void);
extern void test_queue_counters(void);

int main(void) {
	UNITY_BEGIN();
	RUN_TEST(test_queue_setup);
	RUN_TEST(test_queue_push);
	RUN_TEST(test_queue_pop);
	RUN_TEST(test_queue_counters);
	return UNITY_END();
}

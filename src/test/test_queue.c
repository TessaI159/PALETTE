#include <stdlib.h>

#include "Queue.h"
#include "unity.h"

void test_queue_setup(void) {
	struct Queue q;
	TEST_ASSERT_TRUE(Queue_create(&q, 8192));
	TEST_ASSERT_EQUAL_UINT_MESSAGE(8192, q.cap, "Wrong initial capacity");
	TEST_ASSERT_EQUAL_UINT_MESSAGE(0, q.head, "Wrong head value");
	TEST_ASSERT_EQUAL_UINT_MESSAGE(0, q.tail, "Wrong tail value");
	Queue_destroy(&q);
	Queue_create(&q, 8100);
	TEST_ASSERT_EQUAL_UINT_MESSAGE(8192, q.cap,
				       "Wrong rounded initial capacity");
	Queue_destroy(&q);
	Queue_create(&q, -1);
	TEST_ASSERT_EQUAL_UINT_MESSAGE(2, q.cap,
				       "Wrong negative initial capacity");
}

void test_queue_push(void) {
	struct Queue q;
	int	     l = 0;
	Queue_create(&q, 8192);
	for (size_t i = 0; i < q.mask; ++i) {
		int *val = malloc(sizeof(int));
		*val	 = i;
		TEST_ASSERT_TRUE_MESSAGE(Queue_push(&q, (void *)val),
					 "Unable to push");
	}
	TEST_ASSERT_FALSE_MESSAGE(Queue_push(&q, (void *)&l),
				  "Push should have failed");
	free(Queue_pop(&q));
	TEST_ASSERT_FALSE_MESSAGE(Queue_push(&q, NULL),
				  "Queue should not allow pushing NULL");
	Queue_destroy(&q);
}

void test_queue_pop(void) {
	struct Queue q;
	Queue_create(&q, 8192);
	for (size_t i = 0; i < q.mask; ++i) {
		int *val = malloc(sizeof(int));
		*val	 = i;
		Queue_push(&q, (void *)val);
	}

	for (size_t i = 0; i < q.cap - 2; ++i) {
		free(Queue_pop(&q));
	}
	int *test = (int *)Queue_pop(&q);
	TEST_ASSERT_EQUAL_INT_MESSAGE(q.cap - 2, *test,
				      "Popped the wrong value");
	TEST_ASSERT_NULL_MESSAGE(Queue_pop(&q),
				 "Pop should have returned NULL");
	free(test);
	Queue_destroy(&q);
}

void test_queue_counters(void) {
	struct Queue q;
	Queue_create(&q, 8192);
	TEST_ASSERT_TRUE_MESSAGE(Queue_is_empty(&q),
				 "Queue should report empty before loops");
	TEST_ASSERT_FALSE_MESSAGE(Queue_is_full(&q),
				  "Queue should not report full before loops");
	for (size_t i = 0; i < q.mask; ++i) {
		TEST_ASSERT_FALSE_MESSAGE(
		    Queue_is_full(&q),
		    "Queue should not report full in loop in increment loop");
		TEST_ASSERT_EQUAL_UINT_MESSAGE(
		    i, Queue_size(&q),
		    "Queue is reporting the wrong size in increment loop");
		TEST_ASSERT_EQUAL_UINT_MESSAGE(
		    q.cap - i - 1, Queue_empty_slots(&q),
		    "Queue is reporting the wrong number of empty slots in "
		    "increment loop");
		int *val = malloc(sizeof(int));
		*val	 = i;
		Queue_push(&q, val);
		TEST_ASSERT_FALSE_MESSAGE(
		    Queue_is_empty(&q),
		    "Queue should not report empty in increment loop");
		TEST_ASSERT_EQUAL_UINT_MESSAGE(
		    i + 1, Queue_size(&q),
		    "Queue reporting wrong size in increment loop");
	}
	TEST_ASSERT_TRUE_MESSAGE(Queue_is_full(&q),
				 "Queue should report full between loops");
	TEST_ASSERT_FALSE_MESSAGE(
	    Queue_is_empty(&q), "Queue should not report empty between loops");

	size_t q_size = Queue_size(&q);
	for (size_t i = 0; i < q_size; ++i) {
		TEST_ASSERT_FALSE_MESSAGE(
		    Queue_is_empty(&q),
		    "Queue should not report empty in decrement loop");
		free(Queue_pop(&q));
		TEST_ASSERT_EQUAL_UINT_MESSAGE(
		    q_size - i - 1, Queue_size(&q),
		    "Queue is reporting the wrong size in increment loop");
		TEST_ASSERT_EQUAL_UINT_MESSAGE(
		    i + 1, Queue_empty_slots(&q),
		    "Queue is reporting the wrong number of empty slots in "
		    "increment loop");
		TEST_ASSERT_FALSE_MESSAGE(
		    Queue_is_full(&q),
		    "Queue should not report full in loop in decrement loop");
	}

	TEST_ASSERT_TRUE_MESSAGE(
	    Queue_is_empty(&q),
	    "Queue should be empty after removing Queue_size() elements");
	Queue_destroy(&q);
}

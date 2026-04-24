#include <vector>
#include <stdexcept>

#include <math.h>
#include <stdlib.h>

#include <gtest/gtest.h>

#include "ewma_search/interval.h"
#include "ewma_search/ewma.h"
#include "ewma_search/ewma_search.h"

// Do some tests

const double test_tolerance = 1e-10;

TEST(EWMASearch, SameSizeKnapsack) {

	// Check for off-by-one when matching needles to haystacks
	// of the same size.

	std::vector<int> haystack = {4, 3};
	std::vector<int> needle = haystack;
	double alpha = 0;

	EXPECT_EQ(ewma_search(haystack, needle, alpha, test_tolerance).first,
		0);
}

TEST(EWMASearch, DivisionByZeroTest) {

	// If needle is zero at any point, taking it at face value will lead
	// to a divide by zero when calculating the ideal clock.

	std::vector<int> haystack = {1, 2, 3, 4};
	std::vector<int> needle = {0};

	EXPECT_THROW(
		ewma_search(haystack, needle, 0.5, test_tolerance),
		std::invalid_argument);
}

TEST(EWMASearch, LongerNeedleThanHaystack) {
	// This caused a memory access error due to a bug reading
	// outside the haystack. Fixed.

	std::vector<int> needle = std::vector<int>({1, 2, 3, 4});
	std::vector<int> haystack = std::vector<int>({0});

	// Expect no match.

	EXPECT_EQ(
		ewma_search(haystack, needle, 0.5, test_tolerance),
		NO_EWMA_MATCH);
}

TEST(EWMASearch, VerySmallAlphaHang) {
	// Extreme values of alpha could lead get_boundary to hang
	// when determining the lower or upper boundaries or its
	// midpoint.

	double alpha_bugged = 2.64611e-313;
	std::vector<int> haystack = std::vector<int>({16580592, 0});
	std::vector<int> needle = std::vector<int>({19});

	// Expect the first match to be on index zero, because it's
	// always possible to adjust the clocks so that a length one
	// needle matches.

	EXPECT_EQ( ewma_search(haystack, needle, alpha_bugged,
		test_tolerance).first, 0);
}

TEST(EWMASearch, NegativeNeedlesNotHandled) {

	// Negative needles can cause misses due to the way bounds
	// assume that needles are nonnegative. Thus we skip them
	// entirely.

	std::vector<int> haystack = std::vector<int>({-1, 16777215});
	std::vector<int> needle = std::vector<int>({-1});

	EXPECT_THROW(
		ewma_search(haystack, needle, 0.1, test_tolerance),
		std::invalid_argument);
}

TEST(EWMASearch, EmptyNeedleNotMatched) {

	// By convention, I don't want empty needles to match.
	// (Maybe change this later; it shouldn't matter much in
	// any case. What's important is that it doesn't hang or
	// crash anything.)

	std::vector<int> haystack = std::vector<int>({1});
	std::vector<int> needle = std::vector<int>({});

	EXPECT_EQ(ewma_search(haystack, needle, 0.1, test_tolerance),
		NO_EWMA_MATCH);
}


TEST(EWMASearch, HandleSingleInterval) {

	// ewma_search should return the single value of alpha
	// that will match when only a single value will work.
	// This tests that get_boundary will find an initial
	// valid point down to machine epsilon, not just to the
	// given tolerance.

	std::vector<int> haystack = {21, 21, 0};
	std::vector<int> needle = {4, 3};
	double alpha = -2.00089e-11;

	EXPECT_EQ(ewma_search(haystack, needle, alpha, test_tolerance).first,
		0);
}

TEST(EWMASearch, HandleSingleIntervalSecond) {

	// Same, but with a different pattern.

	std::vector<int> haystack = {21, 21};
	std::vector<int> needle = {4, 3};
	double alpha = -8.3259e-16;

	EXPECT_EQ(ewma_search(haystack, needle, alpha, test_tolerance).first,
		0);
}

TEST(EWMASearch, HandleSingleIntervalMarginZero) {

	// Handle a single interval case where the EWMA margin function
	// returns zero on a valid solution even with long double.
	// The bug is fixed by explicitly checking if a clock is
	// valid when one of the bounds returned by the EWMA margin
	// function is zero.

	// I suspect that in an idealized mathematical world, any epsilon
	// below zero would work in this example. But in the real world,
	// we're limited by numerical precision. This is not a problem as
	// long as the brute-force approach is similarly limited, though.

	std::vector<int> haystack = {21, 21};
	std::vector<int> needle = {4, 3};
	double alpha = -1.2e-16;

	EXPECT_EQ(ewma_search(haystack, needle, alpha, test_tolerance).first,
		0);
}
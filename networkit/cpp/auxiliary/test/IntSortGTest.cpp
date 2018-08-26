/*
 * IntSortGTest.cpp
 *
 *  Created on: 08. Aug. 2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
#include <random>

#include "../Random.h"
#include "../IntSort.h"

namespace NetworKit {

class IntSortGTest : public ::testing::TestWithParam<size_t> {};

static void check_intsort(const std::vector<unsigned>& input, unsigned max_key, bool vector_interface) {
	auto sorted = input; // copy input

	if (vector_interface) {
		Aux::intsort(sorted, [] (unsigned x) {return x;}, max_key);
	} else {
		Aux::intsort(sorted.begin(), sorted.end(), [] (unsigned x) {return x;}, max_key);
	}

	auto ref = input;
	std::sort(ref.begin(), ref.end());

	ASSERT_EQ(sorted, ref);
}

static std::vector<unsigned> generate_sequence(size_t n) {
	std::vector<unsigned> vector(n);
	int tmp = 0;
	for(auto& x : vector) {x = tmp++;}
	return vector;
}

TEST_P(IntSortGTest, AlreadySorted) {
	auto n = GetParam();
	auto vec = generate_sequence(n);
	check_intsort(vec, n-1, false);
	check_intsort(vec, n-1, true);
}

TEST_P(IntSortGTest, RandomOrderNoDuplicates) {
	auto n = GetParam();
	auto vec = generate_sequence(n);
	std::shuffle(vec.begin(), vec.end(), Aux::Random::getURNG());
	check_intsort(vec, n-1, false);
	check_intsort(vec, n-1, true);
}

TEST_P(IntSortGTest, RandomSamplesFromNQuarter) {
	auto n = GetParam();
	std::vector<unsigned> vec(n);

	auto& prng = Aux::Random::getURNG();
	std::uniform_int_distribution<unsigned> distr(0, tlx::div_ceil(n-1, 4));
	for(auto& x : vec) {x = distr(prng);}

	check_intsort(vec, n-1, false);
	check_intsort(vec, n-1, true);
}

TEST_P(IntSortGTest, RandomSamplesFromN) {
	auto n = GetParam();
	std::vector<unsigned> vec(n);

	auto& prng = Aux::Random::getURNG();
	std::uniform_int_distribution<unsigned> distr(0, n-1);
	for(auto& x : vec) {x = distr(prng);}

	check_intsort(vec, n-1, false);
	check_intsort(vec, n-1, true);
}

TEST_P(IntSortGTest, RandomSamplesFromMaxInt) {
	auto n = GetParam();
	std::vector<unsigned> vec(n);

	auto& prng = Aux::Random::getURNG();
	std::uniform_int_distribution<unsigned> distr;
	for(auto& x : vec) {x = distr(prng);}

	check_intsort(vec, std::numeric_limits<unsigned>::max(), false);
	check_intsort(vec, std::numeric_limits<unsigned>::max(), true);
}

INSTANTIATE_TEST_CASE_P(IntSortGTest, IntSortGTest, ::testing::Values(0, 1, 8, 1000, 9232, 123451, 1233456),);
}

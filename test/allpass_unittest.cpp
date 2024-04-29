// This sample shows how to write a simple unit test for a function,
// using Google C++ testing framework.
// It is based on https://github.com/google/googletest/blob/main/docs/index.md

#include <limits.h>
#include "dy4.h"
#include "stereo.h"
#include "genfunc.h"
#include "gtest/gtest.h"

namespace {

class AllPass_Fixture: public ::testing::Test {

	public:

	AllPass_Fixture( ) {
	}

	void SetUp( ) {
		
	}

	void TearDown( ) {
	}

	~AllPass_Fixture( )  {
	}

	unsigned short int max_value = 10;
	unsigned char precision = 2;
	float epsilon = 10-2;
};

TEST_F(AllPass_Fixture, AllPassTest){

	std::vector<float> s = {0.0, 1.0};
	std::vector<float> b = {2.0, 3.0, 4.0, 5.0};

	std::vector<float> es = {4.0, 5.0};
	std::vector<float> eo = {0.0, 1.0, 2.0, 3.0};

	std::vector<float> o;

	fast_allpass(o, b, s);

	ASSERT_EQ(b.size(), o.size()) << "The two block vectors are of unequal length";

	for (unsigned int i = 0; i < o.size(); ++i) {
		EXPECT_NEAR(std::abs(eo[i]), std::abs(o[i]), epsilon) << "Unexpected value in block at index " << i;
	}

	for (unsigned int i = 0; i < s.size(); ++i) {
		EXPECT_NEAR(std::abs(es[i]), std::abs(s[i]), epsilon) << "unexpected value in state at index " << i;
	}

}


} // end of namespace

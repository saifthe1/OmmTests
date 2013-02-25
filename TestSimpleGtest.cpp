#include<iostream>
#include<stdio.h>
#include "gtest/gtest.h"

int simpleMultiply(int num)
{
	int n = num*num;
	return n;
}

TEST(SimpleNumberTest,One)
{
	EXPECT_EQ(16,simpleMultiply(4));
}

int main(int argc, char **argv){
	testing::InitGoogleTest(&argc,argv);
	return RUN_ALL_TESTS();
}

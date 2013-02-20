#include<iostream>
#include<stdio.h>
#include "gtest/gtest.h"

int simpleMultiply()
{
	int n = 4*4;
	printf("Multiply %d",n);
	return n;
}

TEST(SimpleNumberTest,One)
{
	EXPECT_EQ(16,simpleMultiply());
}

int main(int argc, char **argv){
	testing::InitGoogleTest(&argc,argv);
	return RUN_ALL_TESTS();
}

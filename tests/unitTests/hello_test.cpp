/**
 * This function is adapted from Google Test's Quickstart Guide.
 * Source: Google Test Quickstart for Google Test and CMake
 * URL: https://google.github.io/googletest/quickstart-cmake.html
 */

#include <gtest/gtest.h>

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

// Main function running the tests
int main(int argc, char **argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

}

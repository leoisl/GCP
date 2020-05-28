#include "gtest/gtest.h"

#include "GenotypeConfidencePercentiler.h"

TEST(BuildPercentiler, NotEnoughData_Fails){
    std::vector<GenotypeConfidence> v{1.2};
    EXPECT_THROW(GenotypeConfidencePercentiler p(v), NotEnoughData);
}



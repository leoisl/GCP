#include "gtest/gtest.h"

#include "GenotypeConfidencePercentiler.h"

TEST(BuildPercentiler, NotEnoughData_Fails){
    std::vector<GenotypeConfidence> v{1.2};
    EXPECT_THROW(GenotypeConfidencePercentiler p(v), NotEnoughData);
}

std::vector<GenotypeConfidence > ordered_v{1.1, 1.5, 5., 5., 10., 12.2, 14., 16., 20., 100.};
GenotypeConfidencePercentiler p(ordered_v);

TEST(PercentileQueries, QueryBelowMin_Get0){
   auto result = p.get_confidence_percentile(1.0);
   EXPECT_FLOAT_EQ(result, 0.);
}

TEST(PercentileQueries, QueryAboveMax_Get100){
    auto result = p.get_confidence_percentile(1000);
    EXPECT_FLOAT_EQ(result, 100.);
}

TEST(PercentileQueries, QueryKnownConfidence_CorrectPercentile){
    auto result = p.get_confidence_percentile(10.);
    EXPECT_FLOAT_EQ(result, 50.);
}

TEST(PercentileQueries, QueryKnownDuplicatedConfidence_CorrectAveragedPercentile){
    auto result = p.get_confidence_percentile(5.);
    EXPECT_FLOAT_EQ(result, 35.);
}


TEST(PercentileQueries, QueryUnKnownConfidence_CorrectInterpolatedPercentile){
    auto result = p.get_confidence_percentile(15.);
    // 15 exactly in between 14 and 16, so expect percentile exactly in between their
    // percentiles of 70 and 80 (resp.)
    EXPECT_FLOAT_EQ(result, 75.);
}

TEST(PercentileQueries, QueryMinElement_GetFirstQuantile){
    auto result = p.get_confidence_percentile(1.1);
    EXPECT_FLOAT_EQ(result, 10.);
}


std::vector<GenotypeConfidence> unordered_v{10., 20., 4., 15.};
GenotypeConfidencePercentiler p2(unordered_v);

TEST(PercentileQueries, QueryKnownConfidenceInUnorderedEntries_CorrectPercentile){
   auto result = p2.get_confidence_percentile(10.);
   EXPECT_FLOAT_EQ(result, 50.);
}

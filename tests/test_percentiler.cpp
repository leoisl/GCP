#include "gtest/gtest.h"

#include "GCP.h"

using namespace GCP;

void check_if_vector_of_doubles_are_equal(const std::vector<double> &v1, const std::vector<double> &v2) {
  EXPECT_EQ(v1.size(), v2.size());
  for (uint32_t i=0; i < v1.size(); ++i) {
    EXPECT_DOUBLE_EQ(v1[i], v2[i]);
  }
}


TEST(BuildPercentiler, NotEnoughData_Fails){
    std::vector<GenotypeConfidence> v{1.2};
    EXPECT_THROW(Percentiler p(v), NotEnoughData);
}

std::vector<GenotypeConfidence > ordered_confidences{1.1, 1.5, 5., 5., 10., 12.2, 14., 16., 20., 100.};
Percentiler ordered_percentiler(ordered_confidences);

TEST(PercentileQueries, QueryBelowMin_Get0){
   auto result = ordered_percentiler.get_confidence_percentile(1.0);
   EXPECT_DOUBLE_EQ(result, 0.);
}

TEST(PercentileQueries, QueryAboveMax_Get100){
    auto result = ordered_percentiler.get_confidence_percentile(1000);
    EXPECT_DOUBLE_EQ(result, 100.);
}

TEST(PercentileQueries, QueryKnownConfidence_CorrectPercentile){
    auto result = ordered_percentiler.get_confidence_percentile(10.);
    EXPECT_DOUBLE_EQ(result, 50.);
}

TEST(PercentileQueries, QueryKnownDuplicatedConfidence_CorrectAveragedPercentile){
    auto result = ordered_percentiler.get_confidence_percentile(5.);
    EXPECT_DOUBLE_EQ(result, 35.);
}


TEST(PercentileQueries, QueryUnKnownConfidence_CorrectInterpolatedPercentile){
    auto result = ordered_percentiler.get_confidence_percentile(15.);
    // 15 exactly in between 14 and 16, so expect percentile exactly in between their
    // percentiles of 70 and 80 (resp.)
    EXPECT_DOUBLE_EQ(result, 75.);
}

TEST(PercentileQueries, QueryUnKnownConfidence_CorrectInterpolatedPercentile2){
  auto result = ordered_percentiler.get_confidence_percentile(14.5);
  // 15 exactly in between 14 and 16, so expect percentile exactly in between their
  // percentiles of 70 and 80 (resp.)
  EXPECT_DOUBLE_EQ(result, 72.5);
}


TEST(PercentileQueries, QueryMinElement_GetFirstQuantile){
    auto result = ordered_percentiler.get_confidence_percentile(1.1);
    EXPECT_DOUBLE_EQ(result, 10.);
}


std::vector<GenotypeConfidence> many_duplicates_confidences{1., 1., 1., 1., 1., 10., 10., 10., 10., 15.};
Percentiler percentile_with_duplicates(many_duplicates_confidences);

TEST(PercentileQueries, QueryPercentileWithDuplicates){
    std::vector<GenotypeConfidence> genotype_confidences({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});

    std::vector<GenotypeConfidence> actual_percentiles =
        percentile_with_duplicates.get_confidence_percentiles(genotype_confidences.begin(), genotype_confidences.end());
    std::vector<GenotypeConfidence> expected_percentiles({0.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 100.0, 100.0, 100.0, 100.0});

    check_if_vector_of_doubles_are_equal(actual_percentiles, expected_percentiles);
}


    auto result = p2.get_confidence_percentile(10.);
    EXPECT_FLOAT_EQ(result, 75.);
}

std::vector<GenotypeConfidence> unordered_v{10., 20., 15., 12., 4.};
Percentiler p3(unordered_v);

TEST(PercentileQueries, QueryKnownConfidenceInUnorderedEntries_CorrectPercentile){
   auto result = p3.get_confidence_percentile(10.);
   EXPECT_FLOAT_EQ(result, 40.);
}

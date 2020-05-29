#include "gtest/gtest.h"

#include "GCP.h"
#include "utils_test.h"

using namespace GCP;

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



TEST(PercentileQueries, QueryPercentileWithDuplicates){
  test_several_confidences(
      {1., 1., 1., 1., 1., 10., 10., 10., 10., 15.},
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
      {0.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 100.0, 100.0, 100.0, 100.0}
      );
}

TEST(PercentileQueries, QueryPercentileWithDuplicates2){
  test_several_confidences(
      {5, 20, 25, 25, 50, 60, 60, 60, 80},
      {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100},
      {0.0, 11.11111111111111, 14.814814814814813, 18.51851851851852, 22.22222222222222, 38.888888888888886, 42.22222222222222, 45.55555555555556, 48.88888888888889, 52.22222222222223, 55.55555555555556, 66.66666666666666, 77.77777777777777, 83.33333333333333, 88.88888888888889, 94.44444444444444, 100.0, 100.0, 100.0, 100.0, 100.0}
  );
}

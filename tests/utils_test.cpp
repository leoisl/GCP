#include "utils_test.h"

void check_if_vector_of_doubles_are_equal(const std::vector<double> &v1, const std::vector<double> &v2) {
  EXPECT_EQ(v1.size(), v2.size());
  for (uint32_t i=0; i < v1.size(); ++i) {
    EXPECT_DOUBLE_EQ(v1[i], v2[i]);
  }
}

void test_several_confidences(
    const std::vector<GenotypeConfidence> &simulated_genotype_confidences,
    const std::vector<GenotypeConfidence> &queried_genotype_confidences,
    const std::vector<GenotypePercentile> &expected_percentiles
) {
  Percentiler percentiler(simulated_genotype_confidences);
  std::vector<GenotypeConfidence> actual_percentiles =
      percentiler.get_confidence_percentiles(queried_genotype_confidences.begin(),
                                             queried_genotype_confidences.end());
  check_if_vector_of_doubles_are_equal(actual_percentiles, expected_percentiles);
}
#ifndef GCP_UTILS_TEST_H
#define GCP_UTILS_TEST_H

#include "gtest/gtest.h"
#include "GCP.h"

using namespace GCP;

void check_if_vector_of_doubles_are_equal(const std::vector<double> &v1, const std::vector<double> &v2);

void test_several_confidences(
    const std::vector<GenotypeConfidence> &simulated_genotype_confidences,
    const std::vector<GenotypeConfidence> &queried_genotype_confidences,
    const std::vector<GenotypePercentile> &expected_percentiles
);

#endif // GCP_UTILS_TEST_H

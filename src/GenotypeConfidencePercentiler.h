#ifndef GCP_GENOTYPECONFIDENCEPERCENTILER_H
#define GCP_GENOTYPECONFIDENCEPERCENTILER_H

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "GenotypeConfidenceSimulator.h"
#include "ranker.h"
#include "custom_types.h"


class NotEnoughData : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

/**
 * Class responsible for assigning confidence percentiles to raw genotype confidences.
 * This class is loosely coupled from the other classes (to satisfy https://github.com/leoisl/GCP/issues/6)
 */
class GenotypeConfidencePercentiler {
public:
  /**
   * Builds a GenotypeConfidencePercentiler from a vector of simulated genotype confidences
   */
  explicit GenotypeConfidencePercentiler(const std::vector<GenotypeConfidence> &simulated_genotype_confidences) :
      simulated_genotype_confidences(simulated_genotype_confidences),
      confidence_to_percentile(get_genotype_confidence_to_percentile_map()),
      min_simulated_genotype_confidence(*std::min_element(simulated_genotype_confidences.begin(), simulated_genotype_confidences.end())),
      max_simulated_genotype_confidence(*std::max_element(simulated_genotype_confidences.begin(), simulated_genotype_confidences.end())) {
    bool we_have_at_least_two_simulated_genotype_confidences = simulated_genotype_confidences.size() >= 2;
    if (not we_have_at_least_two_simulated_genotype_confidences) {
      throw NotEnoughData("Please provide at least two simulated genotype confidences.");
    }
  }


  /**
   * Get the confidence percentile given a genotype confidence.
   */
  GenotypePercentile get_confidence_percentile(GenotypeConfidence genotype_confidence) {
    bool genotype_confidence_is_in_confidence_to_percentile = confidence_to_percentile.find(genotype_confidence) != confidence_to_percentile.end();
    if (not genotype_confidence_is_in_confidence_to_percentile) {
      GenotypePercentile percentile = get_confidence_percentile_through_linear_interpolation(genotype_confidence);
      confidence_to_percentile[genotype_confidence] = percentile;
    }

    GenotypePercentile percentile = confidence_to_percentile[genotype_confidence];
    return percentile;
  }


  // destructor
  virtual ~GenotypeConfidencePercentiler() = default;

  // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
  GenotypeConfidencePercentiler(const GenotypeConfidencePercentiler& other) = delete;
  GenotypeConfidencePercentiler(GenotypeConfidencePercentiler&& other) = delete;
  GenotypeConfidencePercentiler& operator=(const GenotypeConfidencePercentiler& other) = delete;
  GenotypeConfidencePercentiler& operator=(GenotypeConfidencePercentiler&& other) = delete;

private:
  const std::vector<GenotypeConfidence> &simulated_genotype_confidences;
  std::map<GenotypeConfidence, GenotypePercentile> confidence_to_percentile;
  GenotypeConfidence min_simulated_genotype_confidence, max_simulated_genotype_confidence;

  std::vector<GenotypePercentile> get_genotype_confidence_percentiles() const {
    std::vector<double> ranks;
    ranks.reserve(simulated_genotype_confidences.size());
    rank(simulated_genotype_confidences, ranks);

    std::vector<GenotypePercentile> percentiles;
    percentiles.reserve(simulated_genotype_confidences.size());
    for (double a_rank : ranks) {
      GenotypePercentile percentile = 100 * a_rank / simulated_genotype_confidences.size();
      percentiles.push_back(percentile);
    }

    return percentiles;
  }

  std::map<GenotypeConfidence, GenotypePercentile> get_genotype_confidence_to_percentile_map() const {
    std::vector<GenotypePercentile> percentiles = get_genotype_confidence_percentiles();

    std::map<GenotypeConfidence, GenotypePercentile> confidence_to_percentile;
    for (uint32_t i=0; i< simulated_genotype_confidences.size(); ++i) {
      double genotype_confidence = simulated_genotype_confidences[i];
      double percentile = percentiles[i];
      confidence_to_percentile[genotype_confidence] = percentile;
    }

    return confidence_to_percentile;
  }

  GenotypePercentile get_confidence_percentile_through_linear_interpolation(GenotypeConfidence genotype_confidence) const {
    if (genotype_confidence <= min_simulated_genotype_confidence) {
      return 0.0;
    }
    if (genotype_confidence >= max_simulated_genotype_confidence) {
      return 100.0;
    }

    // here, genotype_confidence is between the min and max simulated_genotype_confidence
    // so we are guaranteed a just greater and just smaller elements
    auto just_greater_iterator = confidence_to_percentile.upper_bound(genotype_confidence);
    GenotypeConfidence just_greater_confidence = just_greater_iterator->first;
    GenotypePercentile just_greater_percentile = just_greater_iterator->second;
    auto just_smaller_iterator = just_greater_iterator;
    --just_smaller_iterator;
    GenotypeConfidence just_smaller_confidence = just_smaller_iterator->first;
    GenotypePercentile just_smaller_percentile = just_smaller_iterator->second;

    GenotypePercentile percentile = linear_interpolation(just_greater_confidence, just_greater_percentile, just_smaller_confidence,
                                                         just_smaller_percentile, genotype_confidence);

    return percentile;
  }

  static double linear_interpolation(double big_x, double big_y, double small_x, double small_y, double x) {
    // see https://en.wikipedia.org/wiki/Linear_interpolation for this formula
    return (small_y*(big_x-x) + big_y*(x-small_x)) / (big_x-small_x);
  }
};

#endif // GCP_GENOTYPECONFIDENCEPERCENTILER_H

#ifndef GCP_GENOTYPECONFIDENCEPERCENTILER_H
#define GCP_GENOTYPECONFIDENCEPERCENTILER_H

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "custom_types.h"


class NotEnoughData : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

/**
 * Class responsible for assigning confidence percentiles to raw genotype confidences.
 * This class is decoupled from the others (to satisfy https://github.com/leoisl/GCP/issues/6)
 */
class GenotypeConfidencePercentiler {
public:
  /**
   * Builds a GenotypeConfidencePercentiler from a vector of simulated genotype confidences
   */
  explicit GenotypeConfidencePercentiler(const std::vector<GenotypeConfidence> &simulated_genotype_confidences) :
          entries(simulated_genotype_confidences),
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
      auto lo = std::lower_bound(entries.begin(), entries.end(), genotype_confidence);
      if (lo == entries.end()) return 100.0;
      else if (*lo == genotype_confidence){
          auto hi = std::upper_bound(entries.begin(), entries.end(), genotype_confidence);
          if (lo == hi - 1) return iterator_to_percentile(lo);
          // Case: multiple identical entries, take average
          auto distance = std::distance(lo, hi);
          auto lo_percentile = iterator_to_percentile(lo);
          auto hi_percentile = iterator_to_percentile(--hi);
          return (hi_percentile + lo_percentile) / distance;
      }
      else {
          if (lo == entries.begin()) return 0.0;
          // Case: need to interpolate
          auto hi = lo;
          --lo;
          auto lo_percentile = iterator_to_percentile(lo);
          auto hi_percentile = iterator_to_percentile(hi);
          return linear_interpolation(*hi, hi_percentile, *lo, lo_percentile, genotype_confidence);
      }
  }

  // destructor
  virtual ~GenotypeConfidencePercentiler() = default;

  // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
  GenotypeConfidencePercentiler(const GenotypeConfidencePercentiler& other) = delete;
  GenotypeConfidencePercentiler(GenotypeConfidencePercentiler&& other) = delete;
  GenotypeConfidencePercentiler& operator=(const GenotypeConfidencePercentiler& other) = delete;
  GenotypeConfidencePercentiler& operator=(GenotypeConfidencePercentiler&& other) = delete;

private:
    const std::vector<GenotypeConfidence> &entries;
    using GC_it = std::vector<GenotypeConfidence>::const_iterator;
    GenotypeConfidence min_simulated_genotype_confidence, max_simulated_genotype_confidence;

    std::size_t get_rank(GC_it const& it){
        return std::distance(entries.begin(), it) + 1;
    }

    GenotypePercentile rank_to_percentile(std::size_t const rank) {
        return 100 * rank / entries.size();
    }

    GenotypePercentile iterator_to_percentile(GC_it const& it){
        return rank_to_percentile(get_rank(it));
    };



  static double linear_interpolation(double big_x, double big_y, double small_x, double small_y, double x) {
    // see https://en.wikipedia.org/wiki/Linear_interpolation for this formula
    return (small_y*(big_x-x) + big_y*(x-small_x)) / (big_x-small_x);
  }
};

#endif // GCP_GENOTYPECONFIDENCEPERCENTILER_H

#ifndef GCP_GENOTYPECONFIDENCEPERCENTILER_H
#define GCP_GENOTYPECONFIDENCEPERCENTILER_H

#include "Genotyper.h"
#include "Simulator.h"
#include <vector>
#include <map>
#include "ranker.h"
#include <iostream>

class GenotypeConfidencePercentiler {
public:
  /**
   * Builds a GenotypeConfidencePercentiler given your custom genotyper and simulator.
   * @param genotyper
   * @param simulator
   */
  GenotypeConfidencePercentiler(const Genotyper *genotyper, const Simulator *simulator) :
  genotyper(genotyper), simulator(simulator) {}

  virtual ~GenotypeConfidencePercentiler() = default;

  // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
  GenotypeConfidencePercentiler(const GenotypeConfidencePercentiler& other) = delete;
  GenotypeConfidencePercentiler(GenotypeConfidencePercentiler&& other) = delete;
  GenotypeConfidencePercentiler& operator=(const GenotypeConfidencePercentiler& other) = delete;
  GenotypeConfidencePercentiler& operator=(GenotypeConfidencePercentiler&& other) = delete;

  /**
   * Get the confidence percentile given a genotype confidence.
   * @param genotype_confidence
   * @return confidence percentile
   */
  double get_confidence_percentile(double genotype_confidence) const {
    // TODO
  }


protected:
  const Genotyper *genotyper;
  const Simulator *simulator;

  /**
   * Adapted from https://www.geeksforgeeks.org/rounding-floating-point-number-two-decimal-places-c-c/
   * TODO: refactor to maybe remove this function? It is not guaranteed that the double will always
   *    be rounded to two decimal points, due to FP representation.
   */
  static double round_to_two_decimal_points(double number)
  {
    char buffer[1024];
    std::sprintf(buffer, "%.2f", number);

    double rounded_number;
    std::sscanf(buffer, "%f", &rounded_number);

    return rounded_number;
  }

  std::vector<double> get_genotype_confidence_percentiles(const std::vector<double> &genotype_confidences) const {
    std::vector<double> ranks;
    ranks.reserve(genotype_confidences.size());
    rank(genotype_confidences, ranks);

    std::vector<double> percentiles;
    percentiles.reserve(genotype_confidences.size());
    for (double a_rank : ranks) {
      double percentile = 100 * a_rank / genotype_confidences.size();
      percentiles.push_back(percentile);
    }

    return percentiles;
  }

  std::map<double, double> get_genotype_confidence_to_confidence_percentile_map(
      const std::vector<double> &genotype_confidences) const {
    std::vector<double> percentiles = get_genotype_confidence_percentiles(genotype_confidences);
    std::map<double, double> confidence_to_percentile;

    for (uint32_t i=0; i<genotype_confidences.size(); ++i) {
      double genotype_confidence = genotype_confidences[i];
      double percentile = percentiles[i];
      confidence_to_percentile[genotype_confidence] = round_to_two_decimal_points(percentile);
    }

    return confidence_to_percentile;
  }


};

#endif // GCP_GENOTYPECONFIDENCEPERCENTILER_H

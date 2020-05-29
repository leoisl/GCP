#include <vector>
#include <iostream>

#include "GCP.h"

using namespace GCP;

/**
 * This exemplify using GenotypeConfidencePercentiler without using any of the other classes/concepts;
 * This addresses https://github.com/leoisl/GCP/issues/6
 */

int main() {
  std::vector<GenotypeConfidence> simulated_genotype_confidences =
      {0.5, 2.0, 2.5, 2.5, 5.0, 6.0, 6.0, 6.0, 8.0};
  Percentiler genotype_confidence_percentiler(simulated_genotype_confidences);

  for (double confidence = 0.0; confidence <= 10.0; confidence += 0.5) {
    std::cout << "Genotype confidence: " << confidence
              << "; Percentile: " << genotype_confidence_percentiler.get_confidence_percentile(confidence)
              << std::endl;
  }

  return 0;
}

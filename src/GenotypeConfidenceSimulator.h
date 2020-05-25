#ifndef GCP_GENOTYPECONFIDENCESIMULATOR_H
#define GCP_GENOTYPECONFIDENCESIMULATOR_H

#include "Genotyper.h"
#include "GenotyperInputProducer.h"
#include "Model.h"
#include "custom_types.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

/**
 * Class responsible for simulating genotype confidences given a Model, a GenotyperInputProducer and a Genotyper.
 */
template <typename ModelData, typename GenotyperInput>
class GenotypeConfidenceSimulator {
public:
  // ctor/dtor
  GenotypeConfidenceSimulator(
            Model<ModelData> *model,
            GenotyperInputProducer<ModelData, GenotyperInput> *genotyper_input_producer,
            Genotyper<GenotyperInput> * genotyper) :
      model(model), genotyper_input_producer(genotyper_input_producer), genotyper(genotyper) {
  }
  virtual ~GenotypeConfidenceSimulator() = default;

  // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
  GenotypeConfidenceSimulator(const GenotypeConfidenceSimulator & other) = delete;
  GenotypeConfidenceSimulator(GenotypeConfidenceSimulator && other) = delete;
  GenotypeConfidenceSimulator & operator=(const GenotypeConfidenceSimulator & other) = delete;
  GenotypeConfidenceSimulator & operator=(GenotypeConfidenceSimulator && other) = delete;


  // main method
  std::vector<GenotypeConfidence> simulate(uint32_t iterations=10000) {
    std::vector<GenotypeConfidence> genotype_confidences;
    genotype_confidences.reserve(iterations);

    for (uint32_t iterations_done = 0; iterations_done < iterations; ++iterations_done) {
      ModelData model_data = model->produce_data();
      GenotyperInput genotyper_input = genotyper_input_producer->produce_input(model_data);
      GenotypeConfidence genotype_confidence = genotyper->get_genotype_confidence(genotyper_input);
      genotype_confidences.push_back(genotype_confidence);
    }

    std::sort(genotype_confidences.begin(), genotype_confidences.end());
    return genotype_confidences;
  }

private:
  Model<ModelData> *model;
  GenotyperInputProducer<ModelData, GenotyperInput> *genotyper_input_producer;
  Genotyper<GenotyperInput> * genotyper;
};

#endif // GCP_GENOTYPECONFIDENCESIMULATOR_H

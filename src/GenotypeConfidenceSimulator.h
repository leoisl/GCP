#ifndef GCP_GENOTYPECONFIDENCESIMULATOR_H
#define GCP_GENOTYPECONFIDENCESIMULATOR_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "Model.h"
#include "custom_types.h"

/**
 * Class responsible for simulating genotype confidences given a Model, a GenotyperInputProducer and a Genotyper.
 */
template <typename ModelData, typename Genotyper>
class GenotypeConfidenceSimulator {
public:
    // ctor/dtor
    explicit GenotypeConfidenceSimulator(Model<ModelData> *model) : model(model) {}
    virtual ~GenotypeConfidenceSimulator() = default;

    // main method
    std::vector<GenotypeConfidence> simulate(uint32_t iterations) {
        std::vector<GenotypeConfidence> genotype_confidences;
        genotype_confidences.reserve(iterations);

        for (uint32_t iterations_done = 0; iterations_done < iterations; ++iterations_done) {
            ModelData model_data = model->produce_data();
            Genotyper gtyped(model_data);
            genotype_confidences.push_back(gtyped.get_genotype_confidence());
        }

        std::sort(genotype_confidences.begin(), genotype_confidences.end());
        return genotype_confidences;
    }


    // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
    GenotypeConfidenceSimulator(const GenotypeConfidenceSimulator & other) = delete;
    GenotypeConfidenceSimulator(GenotypeConfidenceSimulator && other) = delete;
    GenotypeConfidenceSimulator & operator=(const GenotypeConfidenceSimulator & other) = delete;
    GenotypeConfidenceSimulator & operator=(GenotypeConfidenceSimulator && other) = delete;

private:
    Model<ModelData> *model;
};

#endif // GCP_GENOTYPECONFIDENCESIMULATOR_H

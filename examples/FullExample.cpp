#include "Model.h"
#include "GenotypeConfidencePercentiler.h"
#include "GenotypeConfidenceSimulator.h"
#include "Genotyper.h"
#include "custom_types.h"
#include <cmath>
#include <iostream>
#include <string>


/**
 * Defines the data produced by our custom Model (you can define whatever class/struct here).
 * Pandora will produce a set of data in its model, gramtools another set of data, etc...
 * In this example, we model the kmer coverage.
 */
struct MyModelData {
  uint32_t correct_coverage, incorrect_coverage;
  MyModelData(uint32_t correct_coverage, uint32_t incorrect_coverage) :
      correct_coverage(correct_coverage), incorrect_coverage(incorrect_coverage) {}
};


/**
 * This inheritance and implementation will teach the library how to produce data from your model.
 * You can define whatever Model you want here.
 * This is a simple random Model.
 */
class MyModel : public Model<MyModelData> {
  using Model<MyModelData>::Model;
  virtual MyModelData produce_data() override {
    uint32_t correct_coverage = random_number_generator()%1000 + 1;
    uint32_t incorrect_coverage = random_number_generator()%100;
    return MyModelData(correct_coverage, incorrect_coverage);
  }
};


/**
 * This data structure represents whatever Data Structure the C++ implementation
 * of the variant caller uses. This is the data that is given to the genotyper
 * of the variant caller to be genotyped. It is built from the model data and any other data.
 * Each variant caller represent this data in different way (Pandora's
 * data structure is different from gramtools', which is different from whatever
 * other caller.
 * In this simple example, we just have the likelihoods of two alleles.
 */
struct MyGenotyperInput {
  double likelihood_1, likelihood_2;
  MyGenotyperInput(double likelihood_1, double likelihood_2) :
     likelihood_1(likelihood_1), likelihood_2(likelihood_2) {}
};


/**
 * This inheritance and implementation will teach the library how to build
 * an input to the genotyper using the model data and any other data.
 * This bridges the Model to the Genotyper.
 * This is a very dummy example of a GenotyperInputProducer.
 */
class MyGenotyperInputProducer : public GenotyperInputProducer<MyModelData, MyGenotyperInput> {
public:
  virtual MyGenotyperInput produce_input(const MyModelData &model_data) override {
    return MyGenotyperInput(model_data.correct_coverage, model_data.incorrect_coverage);
  }
};


/**
 * This inheritance and implementation will teach the library how to use the
 * variant caller's genotyper to get the genotype confidence.
 *
 * THIS IS REQUIRED FOR THE CLIENT TO IMPLEMENT.
 */
class MyGenotyper : public Genotyper<MyGenotyperInput> {
public:
  virtual double get_genotype_confidence (const MyGenotyperInput& data) override {
    return std::abs(data.likelihood_1 - data.likelihood_2);
  }
};



/**
 * Usage of the library
 */
int main() {
  // create the Model, the GenotyperInputProducer, and the Genotyper
  // we use polymorphism  for our library to execute the previously user defined code (thus the pointers)
  Model<MyModelData>* my_model = new MyModel();
  GenotyperInputProducer<MyModelData, MyGenotyperInput>* my_genotyper_input_producer = new MyGenotyperInputProducer();
  Genotyper<MyGenotyperInput>* my_genotyper = new MyGenotyper();

  // create the GenotypeConfidenceSimulator
  GenotypeConfidenceSimulator<MyModelData, MyGenotyperInput> genotype_confidence_simulator(my_model, my_genotyper_input_producer, my_genotyper);

  // simulate 10000 confidences according to our Model and our Genotyper
  std::vector<GenotypeConfidence> simulated_confidences = genotype_confidence_simulator.simulate();

  // create the GenotypeConfidencePercentiler
  GenotypeConfidencePercentiler genotype_confidence_percentiler(simulated_confidences);

  // query some genotype confidence percentiles
  for (double confidence = 0.0; confidence <= 1000.0; confidence += 100.0) {
    std::cout << "Genotype confidence: " << confidence
              << "; Percentile: " << genotype_confidence_percentiler.get_confidence_percentile(confidence)
              << std::endl;
  }

  // cleanup
  delete my_model;
  delete my_genotyper_input_producer;
  delete my_genotyper;

  return 0;
}
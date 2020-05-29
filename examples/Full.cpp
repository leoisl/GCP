#include <iostream>
#include <string>

#include "GCP.h"

/**
 * Define data struct/class required for genotyping to run, and
 * that will get produced by simulation.
 * Assumption: pandora/gramtools have defined this internally.
 */
struct MyModelData {
  uint32_t correct_coverage, incorrect_coverage;
  MyModelData(uint32_t correct_coverage, uint32_t incorrect_coverage) :
      correct_coverage(correct_coverage), incorrect_coverage(incorrect_coverage) {}
};


/**
 * Assumption: pandora/gramtools have defined this internally
 */
class MyGenotyper{
private:
    double gt_conf = 0.;
public:
    /**
     * This must be defined for the Simulator to use
     */
    MyGenotyper(MyModelData input){
        double likelihood1 = input.correct_coverage;
        double likelihood2 = input.incorrect_coverage;
        likelihood1 >= likelihood2 ? gt_conf = likelihood1 - likelihood2
                : gt_conf = likelihood2 - likelihood1;
    }
    /**
     * This must be defined for the Simulator to use
     */
    double get_genotype_confidence() {return gt_conf;}
};


/**
 * REQUIRED: Teach the library how to produce simulated data.
 * This example samples coverage uniformly at random.
 */
class MyModel : public Model<MyModelData> {
  using Model<MyModelData>::Model;
  MyModelData produce_data() override {
    uint32_t correct_coverage = random_number_generator()%1000 + 1;
    uint32_t incorrect_coverage = random_number_generator()%100;
    return MyModelData(correct_coverage, incorrect_coverage);
  }
};


/**
 * Usage of the library
 */
int main() {
  // create the Model, the GenotyperInputProducer, and the Genotyper
  // we use polymorphism  for our library to execute the previously user defined code (thus the pointers)
  Model<MyModelData>* my_model = new MyModel();

  // simulate 10000 confidences according to our Model and our Genotyper
  Simulator<MyModelData, MyGenotyper> genotype_confidence_simulator(my_model);
  std::vector<GenotypeConfidence> simulated_confidences = genotype_confidence_simulator.simulate(10000);

  // create the Percentiler
  Percentiler genotype_confidence_percentiler(simulated_confidences);

  // query some genotype confidence percentiles
  for (double confidence = 0.0; confidence <= 1000.0; confidence += 100.0) {
    std::cout << "Genotype confidence: " << confidence
              << "; Percentile: " << genotype_confidence_percentiler.get_confidence_percentile(confidence)
              << std::endl;
  }

  // cleanup
  delete my_model;

  return 0;
}
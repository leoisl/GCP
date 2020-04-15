#include "Genotyper.h"
#include "Simulator.h"
#include "GenotypeConfidencePercentiler.h"
#include <string>
#include <iostream>
#include <cmath>


/**
 * This data structure represents whatever Data Structure the C++ implementation
 * of the variant caller uses. This is the data that is given to the genotyper
 * of the variant caller to be genotyped.
 * This is just for illustration purposes, this Data Structure does not need to
 * be defined by the client of the library (it is probably already defined
 * in the variant caller code).
 * Each variant caller represent this data in different way (Pandora's
 * data structure is different from gramtools', which is different from whatever
 * other caller.
 * In this simple example, we just have REF and ALT coverage.
 *
 * THIS IS NOT REQUIRED FOR THE CLIENT TO IMPLEMENT (USUALLY ALREADY IMPLEMENTED).
 */
class MyDataToBeGenotyped {
public:
  MyDataToBeGenotyped(uint32_t ref_coverage, uint32_t alt_coverage):
  ref_coverage(ref_coverage), alt_coverage(alt_coverage) {}
  const uint32_t ref_coverage;
  const uint32_t alt_coverage;
};


/**
 * This inheritance and implementation will teach the library how to build
 * and destroy simulated data to be given to the genotyper of this variant caller.
 * Simulated data should be the same data structure that the genotyper of this
 * variant caller receives.
 * In this example, it is an instance of MyDataToBeGenotyped.
 *
 * THIS IS REQUIRED FOR THE CLIENT TO IMPLEMENT.
 *
 * TODO: CHECK IF THIS INHERITANCE GETS SIMPLIFIED TO THE CLIENT THROUGH TEMPLATES.
 */
class MySimulator : public Simulator {
public:
  using Simulator::Simulator;
  virtual void* build_simulation_data_for_genotyper(
      uint32_t correct_coverage, uint32_t incorrect_coverage) const {
    MyDataToBeGenotyped* data = new MyDataToBeGenotyped(correct_coverage, incorrect_coverage);
    return static_cast<void*>(data);
  }

  virtual void destroy_simulation_data(void *data) const  {
    MyDataToBeGenotyped* casted_data = static_cast<MyDataToBeGenotyped*>(data);
    delete casted_data;
  }
};


/**
 * This inheritance and implementation will teach the library how to use the
 * variant caller's genotyper to get the genotype confidence.
 * The genotyper is expected to receive a data structure of the same type
 * of the simulated data. In this example, we receive an instance of
 * MyDataToBeGenotyped, and return the genotype confidence.
 *
 * THIS IS REQUIRED FOR THE CLIENT TO IMPLEMENT.
 */
class MyGenotyper : public Genotyper {
public:
  virtual double get_genotype_confidence (const void* data) const {
    MyDataToBeGenotyped* data_to_be_genotyped = (MyDataToBeGenotyped*)(data);
    return std::abs(double(data_to_be_genotyped->ref_coverage - data_to_be_genotyped->alt_coverage));
  }
};





/**
 * Usage of the library
 */
int main() {
  // create Genotyper, Simulator and  GenotypeConfidencePercentiler
  // we use polymorphism  for our library to execute the previously user defined code
  Genotyper* genotyper = new MyGenotyper();
  double mean_depth=50, variance_depth=100, error_rate=0.11;
  Simulator* simulator = new MySimulator(mean_depth, variance_depth, error_rate);
  GenotypeConfidencePercentiler gcp(genotyper, simulator);

  // query some genotype confidence percentiles
  std::cout << gcp.get_confidence_percentile(10) << std::endl;
  std::cout << gcp.get_confidence_percentile(50) << std::endl;
  std::cout << gcp.get_confidence_percentile(100) << std::endl;
  std::cout << gcp.get_confidence_percentile(150) << std::endl;
  std::cout << gcp.get_confidence_percentile(200) << std::endl;

  // cleanup
  delete genotyper;
  delete simulator;

  return 0;
}
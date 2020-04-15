#ifndef GCP_SIMULATOR_H
#define GCP_SIMULATOR_H

#include <iostream>
#include <algorithm>
#include <random>
#include <cmath>
#include <vector>
#include "Genotyper.h"

class Simulator {
public:
  Simulator(double mean_depth, double variance_depth, double error_rate, uint32_t seed=42, uint32_t iterations=10000) :
      mean_depth(mean_depth), variance_depth(variance_depth), error_rate(error_rate), random_number_generator(seed), iterations(iterations){
    fix_mean_and_variance_depth();
  }
  virtual ~Simulator() = default;

  // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
  Simulator(const Simulator& other) = delete;
  Simulator(Simulator&& other) = delete;
  Simulator& operator=(const Simulator& other) = delete;
  Simulator& operator=(Simulator&& other) = delete;


  std::vector<double> simulate_genotype_confidences(const Genotyper* genotyper) {
    std::vector<double> genotype_confidences;
    genotype_confidences.reserve(iterations);

    double number_of_successes = get_number_of_successes();
    double probability_of_success = get_probability_of_success();
    std::negative_binomial_distribution<uint32_t> negative_binomial_distribution(number_of_successes,
                                                                                 probability_of_success);
    std::binomial_distribution<uint32_t> binomial_distribution(
        number_of_successes, probability_of_success);

    for (uint32_t iterations_done = 0; iterations_done < iterations; ) {
      uint32_t correct_coverage = negative_binomial_distribution(random_number_generator);
      uint32_t incorrect_coverage = binomial_distribution(random_number_generator);

      bool simulated_coverage_is_zero = correct_coverage==0 && incorrect_coverage==0;
      if (simulated_coverage_is_zero) {
        continue;
      }

      void* simulation_data = build_simulation_data_for_genotyper(correct_coverage, incorrect_coverage);
      double genotype_confidence = genotyper->get_genotype_confidence(simulation_data);
      destroy_simulation_data(simulation_data);

      genotype_confidences.push_back(genotype_confidence);
      iterations_done += 1;
    }

    std::sort(genotype_confidences.begin(), genotype_confidences.end());
    return genotype_confidences;
  }

protected:
  double mean_depth;
  double variance_depth;
  double error_rate;
  std::default_random_engine random_number_generator;
  uint32_t iterations;

  /**
   * To be inherited and implemented by the client.
   * Build the simulation data for the genotyper, which is inherited and also implementd by the client.
   * Remember that you also have access to instance variables such as mean_depth,
   *    variance_depth, error_rate.
   * @param correct_coverage: correct coverage from negative binomial distribution
   * @param incorrect_coverage: incorrect coverage from negative binomial distribution
   * @return data to be given to the genotyper.
   */
  virtual void* build_simulation_data_for_genotyper(
      uint32_t correct_coverage, uint32_t incorrect_coverage) const = 0;

  /**
   * To be inherited and implemented by the client.
   * Destroy the simulation data to cleanup memory (e.g. delete, or free).
   * TODO: is there a way to not need the client to do this? Maybe use templates instead of void pointers?
   * @param data
   */
  virtual void destroy_simulation_data(void *data) const = 0;


  void fix_mean_and_variance_depth() {
    if (variance_depth < mean_depth) {
      variance_depth = 2.0 * mean_depth;
      std::cerr << "Variance in read depth is smaller than mean read depth. "
                  "Setting variance = 2 * mean, so that variant simulations can run. "
                  "GT_CONF_PERCENTILE may not be very useful as a result of this.";
    }
  }


  double get_number_of_successes() const {
    return std::pow(mean_depth, 2) / (variance_depth - mean_depth);
  }


  double get_probability_of_success() const {
    return 1.0 - (variance_depth - mean_depth) / variance_depth;
  }
};

#endif // GCP_SIMULATOR_H

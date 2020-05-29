#ifndef GCP_MODEL_H
#define GCP_MODEL_H

#include <random>

/**
 * Class responsible for producing data.
 * This data will be used to produce inputs to a Genotyper.
 * The client is responsible for implementing `produce_data`;
 * eg in Pandora, we model k-mer coverage with a negative binomial distribution.
 */
template <typename ModelData>
class Model {
public:
  explicit Model(uint32_t seed=42) : random_number_generator(seed) {}

  /**
   * To be inherited and implemented by the client.
   * Produce ModelData to be fed to the GenotyperInputProducer to produce input to the Genotyper.
   * Note that you have access to a random number generator.
   * As the ModelData returned is a template, this class can model anything.
   */
  virtual ModelData produce_data() = 0;

  // destructor
  virtual ~Model() = default;

  // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
  Model(const Model& other) = delete;
  Model(Model&& other) = delete;
  Model& operator=(const Model& other) = delete;
  Model& operator=(Model&& other) = delete;

protected:
  std::default_random_engine random_number_generator;
};

#endif // GCP_MODEL_H

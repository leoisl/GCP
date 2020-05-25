#ifndef GCP_GENOTYPERINPUTPRODUCER_H
#define GCP_GENOTYPERINPUTPRODUCER_H

/**
 * Class responsible for producing input data to the genotyper, given a ModelData.
 */
template <typename ModelData, typename GenotyperInput>
class GenotyperInputProducer {
public:
  /**
   * To be inherited and implemented by the client.
   * Abstract method to produce a GenotyperInput given a ModelData.
   */
  virtual GenotyperInput produce_input(const ModelData &model_data) = 0;

  // ctor/dtor
  GenotyperInputProducer() = default;
  virtual ~GenotyperInputProducer() = default;

  // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
  GenotyperInputProducer(const GenotyperInputProducer & other) = delete;
  GenotyperInputProducer(GenotyperInputProducer && other) = delete;
  GenotyperInputProducer & operator=(const GenotyperInputProducer & other) = delete;
  GenotyperInputProducer & operator=(GenotyperInputProducer && other) = delete;
};

#endif // GCP_GENOTYPERINPUTPRODUCER_H

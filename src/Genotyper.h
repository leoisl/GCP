#ifndef GCP_GENOTYPER_H
#define GCP_GENOTYPER_H

#include "GenotyperInputProducer.h"

/**
 * Generic class that represents a genotyper that genotypes data.
 * GenotyperInput: type of the data the genotyper accepts.
 */
template <typename GenotyperInput>
class Genotyper {
public:
  /**
   * To be inherited and implemented by the client.
   * Compute the genotype confidence given the data.
   */
  virtual double get_genotype_confidence (const GenotyperInput& data) = 0;

  // ctor/dtor
  Genotyper() = default;
  virtual ~Genotyper() = default;

  // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
  Genotyper(const Genotyper& other) = delete;
  Genotyper(Genotyper&& other) = delete;
  Genotyper& operator=(const Genotyper& other) = delete;
  Genotyper& operator=(Genotyper&& other) = delete;
};

#endif // GCP_GENOTYPER_H

#ifndef GCP_GENOTYPER_H
#define GCP_GENOTYPER_H


/**
 * Generic class that represents a genotyper that genotypes VCF records or generic data.
 */
class Genotyper {
public:
  /**
   * Concrete genotyper configuration must be initialised here.
   */
  Genotyper() = default;

  /**
   * Concrete genotyper cleanup.
   */
  virtual ~Genotyper() = default;

  // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
  Genotyper(const Genotyper& other) = delete;
  Genotyper(Genotyper&& other) = delete;
  Genotyper& operator=(const Genotyper& other) = delete;
  Genotyper& operator=(Genotyper&& other) = delete;

  /**
   * Compute the genotype confidence given the data.
   * TODO: there are safer ways to do this instead of using void*, which can be error prone,
   *        but all options involve adding dependencies or updating C++ version, e.g.:
   *        1. https://en.cppreference.com/w/cpp/utility/any (requires C++17);
   *        2. https://www.boost.org/doc/libs/1_72_0/doc/html/any.html#id-1.3.5.3 (requires boost);
   */
  virtual double get_genotype_confidence (const void* data) const = 0;
};

#endif // GCP_GENOTYPER_H

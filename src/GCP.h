#ifndef GCP_H
#define GCP_H

#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>

namespace GCP {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* ________INTERFACE________ */

/*
 * This is the Interface of the Genotype Confidence Percentiler header-only library.
 * You just need to look at the interface, no need to look at the implementation.
 */

////////////////////////////////////////////////////////////////////////////////
// typedefs
using GenotypeConfidence = double;
using GenotypePercentile = double;
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
/* ________Model________ */

/**
* Class responsible for producing data used to be genotyped by a Genotyper.
* The client is responsible for implementing `produce_data`;
* eg in Pandora, we model k-mer coverage with a negative binomial distribution.
*/
template<typename ModelData>
class Model {
public:
    explicit inline Model(uint32_t seed = 42);

    /**
     * To be inherited and implemented by the client.
     * Note that you have access to a random number generator.
     * As the ModelData returned is a template, this class can model anything.
     */
    virtual ModelData produce_data() = 0;

    // destructor
    virtual ~Model() = default;

    // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
    Model(const Model &other) = delete;
    Model(Model &&other) = delete;
    Model &operator=(const Model &other) = delete;
    Model &operator=(Model &&other) = delete;

protected:
    std::default_random_engine random_number_generator;
};
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
/* ________Simulator________ */

/**
* Class responsible for simulating genotype confidences given a Model and a Genotyper.
*/
template<typename ModelData, typename Genotyper>
class Simulator {
public:
    // ctor/dtor
    explicit inline Simulator(Model<ModelData> * model);
    virtual ~Simulator() = default;

    // main method
    inline std::vector<GenotypeConfidence> simulate(uint32_t iterations);

    // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
    Simulator(const Simulator &other) = delete;
    Simulator(Simulator &&other) = delete;
    Simulator &operator=(const Simulator &other) = delete;
    Simulator &operator=(Simulator &&other) = delete;

private:
    Model<ModelData> *model;
};
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
/* ________Percentiler________ */

class NotEnoughData : public std::runtime_error {
    using std::runtime_error::runtime_error;
};


/**
* Class responsible for assigning confidence percentiles to raw genotype confidences.
*/
class Percentiler {
public:
    /**
     * Builds a GenotypeConfidencePercentiler from a vector of simulated genotype confidences
     */
    explicit inline Percentiler(const std::vector<GenotypeConfidence> &unsorted_genotype_confidences);


    /**
     * Get the confidence percentile given a genotype confidence.
     */
    inline GenotypePercentile get_confidence_percentile(GenotypeConfidence queried_confidence) const;

    /**
     * Utility method: get several confidence percentiles.
     */
    template<typename Iterator_Type>
    inline std::vector<GenotypePercentile> get_confidence_percentiles(Iterator_Type begin, Iterator_Type end) const;


    // destructor
    virtual ~Percentiler() = default;

    // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
    Percentiler(const Percentiler &other) = delete;
    Percentiler(Percentiler &&other) = delete;
    Percentiler &operator=(const Percentiler &other) = delete;
    Percentiler &operator=(Percentiler &&other) = delete;

private:
    std::map<GenotypeConfidence, GenotypePercentile> confidence_to_percentile;

    // helper methods
    static inline std::map<GenotypeConfidence, GenotypePercentile> create_confidence_to_percentile_map (const std::vector<GenotypeConfidence> &unsorted_genotype_confidences);
    static inline GenotypePercentile iterator_to_percentile(const std::vector<GenotypeConfidence>::const_iterator &it,
            const std::vector<GenotypeConfidence> &genotype_confidences);
    static inline double linear_interpolation(double x1, double x2, double y1, double y2, double x);
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




/*
 * YOU HAVE REACHED THE END OF THE INTERFACE.
 * BELOW ARE IMPLEMENTATION DETAILS USERS OF THIS LIBRARY DOES NOT NEED TO LOOK.
 */

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////






























////////////////////////////////////////////////////////////////////////////////
// Implementation details
////////////////////////////////////////////////////////////////////////////////


/* ________Model________ */
template<typename ModelData>
Model<ModelData>::Model(uint32_t seed) : random_number_generator(seed) {}


/* ________Simulator________ */
template<typename ModelData, typename Genotyper>
Simulator<ModelData, Genotyper>::Simulator(Model<ModelData> * model) : model(model) {}

template<typename ModelData, typename Genotyper>
std::vector<GenotypeConfidence> Simulator<ModelData, Genotyper>::simulate(uint32_t iterations) {
  std::vector<GenotypeConfidence> genotype_confidences;
  genotype_confidences.reserve(iterations);

  for (uint32_t iterations_done = 0; iterations_done < iterations; ++iterations_done) {
    ModelData model_data = model->produce_data();
    Genotyper genotyper(model_data);
    GenotypeConfidence genotype_confidence = genotyper.get_genotype_confidence();
    genotype_confidences.push_back(genotype_confidence);
  }

  std::sort(genotype_confidences.begin(), genotype_confidences.end());
  return genotype_confidences;
}


/* ________Percentiler________ */
Percentiler::Percentiler(const std::vector<GenotypeConfidence> &unsorted_genotype_confidences) {
  bool enough_data = unsorted_genotype_confidences.size() >= 2;
  if (not enough_data) {
    throw NotEnoughData("Please provide at least two simulated genotype confidences.");
  }
  this->confidence_to_percentile = create_confidence_to_percentile_map(unsorted_genotype_confidences);

}

GenotypePercentile Percentiler::get_confidence_percentile(GenotypeConfidence queried_confidence) const {
  auto upper_bound = confidence_to_percentile.upper_bound(queried_confidence);

  bool queried_confidence_is_lower_or_equal_to_min_confidence = upper_bound == confidence_to_percentile.begin();
  if (queried_confidence_is_lower_or_equal_to_min_confidence) {
    return 0.0;
  }

  bool queried_confidence_is_larger_than_max_confidence = upper_bound == confidence_to_percentile.end();
  if (queried_confidence_is_larger_than_max_confidence) {
    return 100.0;
  }

  double confidence = upper_bound->first;
  double percentile = upper_bound->second;
  bool queried_confidence_is_found = confidence == queried_confidence;
  if (queried_confidence_is_found) {
    return percentile;
  }

  // interpolation has to be done now
  auto lower_bound = upper_bound;
  --lower_bound;
  return linear_interpolation(upper_bound->first, lower_bound->first,
                              upper_bound->second, lower_bound->second,
                              queried_confidence);
}

template<typename Iterator_Type>
std::vector<GenotypePercentile> Percentiler::get_confidence_percentiles(Iterator_Type begin, Iterator_Type end) const {
  std::vector<GenotypePercentile> confidence_percentiles;
  for(Iterator_Type it = begin; it != end; ++it) {
    GenotypePercentile percentile = this->get_confidence_percentile(*it);
    confidence_percentiles.push_back(percentile);
  }
  return confidence_percentiles;
}

std::map<GenotypeConfidence, GenotypePercentile> Percentiler::create_confidence_to_percentile_map (const std::vector<GenotypeConfidence> &unsorted_genotype_confidences) {
  std::map<GenotypeConfidence, GenotypePercentile> confidence_to_percentile;

  // sort the confidences
  std::vector<GenotypeConfidence> genotype_confidences(unsorted_genotype_confidences);
  std::sort(genotype_confidences.begin(), genotype_confidences.end());

  auto it_current_confidence = genotype_confidences.begin();
  while (it_current_confidence != genotype_confidences.end()) {
    // it_confidence_just_greater_than_current points to the first confidence greater than it_current_confidence
    auto it_confidence_just_greater_than_current = std::upper_bound(genotype_confidences.begin(),
                                                                    genotype_confidences.end(), *it_current_confidence);
    auto it_current_confidence_highest_percentile = it_confidence_just_greater_than_current - 1;

    GenotypePercentile lowest_percentile_for_this_confidence = iterator_to_percentile(it_current_confidence, genotype_confidences);
    GenotypePercentile highest_percentile_for_this_confidence = iterator_to_percentile(it_current_confidence_highest_percentile, genotype_confidences);
    GenotypePercentile average_percentile_for_this_Confidence =
        lowest_percentile_for_this_confidence + (highest_percentile_for_this_confidence - lowest_percentile_for_this_confidence) / 2;
    confidence_to_percentile[*it_current_confidence] = average_percentile_for_this_Confidence;


    it_current_confidence = it_confidence_just_greater_than_current;
  }

  return confidence_to_percentile;
}

GenotypePercentile Percentiler::iterator_to_percentile(const std::vector<GenotypeConfidence>::const_iterator &it,
                                                 const std::vector<GenotypeConfidence> &genotype_confidences) {
  auto rank = std::distance(genotype_confidences.begin(), it) + 1;
  return 100. * rank / genotype_confidences.size();
};

double Percentiler::linear_interpolation(double x1, double x2, double y1, double y2, double x) {
  double slope = (y2 - y1) / (x2 - x1);
  return y1 + slope * (x - x1);
}

} // end namespace GCP
#endif // GCP_H

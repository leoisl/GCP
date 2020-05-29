#ifndef GCP_H
#define GCP_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>

namespace GCP {

using GenotypeConfidence = double;
using GenotypePercentile = double;

/* ________Model________ */

/**
* Class responsible for producing data used by a Genotyper.
* The client is responsible for implementing `produce_data`;
* eg in Pandora, we model k-mer coverage with a negative binomial distribution.
*/
template<typename ModelData>
class Model {
public:
    explicit Model(uint32_t seed = 42) : random_number_generator(seed) {}

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


/* ________Simulator________ */

/**
* Class responsible for simulating genotype confidences given a Model and a Genotyper.
*/
template<typename ModelData, typename Genotyper>
class Simulator {
public:
    // ctor/dtor
    explicit Simulator(Model<ModelData> *model) : model(model) {}

    virtual ~Simulator() = default;

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
    Simulator(const Simulator &other) = delete;

    Simulator(Simulator &&other) = delete;

    Simulator &operator=(const Simulator &other) = delete;

    Simulator &operator=(Simulator &&other) = delete;

private:
    Model<ModelData> *model;
};


/* ________Percentiler________ */

class NotEnoughData : public std::runtime_error {
    using std::runtime_error::runtime_error;
};


/**
* Class responsible for assigning confidence percentiles to raw genotype confidences.
* This class is decoupled from the others (to satisfy https://github.com/leoisl/GCP/issues/6)
*/
class Percentiler {
public:
    /**
     * Builds a GenotypeConfidencePercentiler from a vector of simulated genotype confidences
     */
    explicit Percentiler(const std::vector<GenotypeConfidence> &input_entries) :
            entries(input_entries),
            min_entry(*std::min_element(input_entries.begin(), input_entries.end())),
            max_entry(*std::max_element(input_entries.begin(), input_entries.end())) {

        bool enough_data = input_entries.size() >= 2;
        if (not enough_data) {
            throw NotEnoughData("Please provide at least two simulated genotype confidences.");
        }
    }


    /**
     * Get the confidence percentile given a genotype confidence.
     */
    GenotypePercentile get_confidence_percentile(GenotypeConfidence query) {
        auto lo = std::lower_bound(entries.begin(), entries.end(), query);
        if (lo == entries.end()) return 100.0;
        else if (*lo == query) {
            auto hi = std::upper_bound(entries.begin(), entries.end(), query);
            if (lo == hi - 1) return iterator_to_percentile(lo);
            // Case: multiple identical entries, take average
            auto lo_percentile = iterator_to_percentile(lo);
            auto hi_percentile = iterator_to_percentile(--hi);
            return lo_percentile + (hi_percentile - lo_percentile) / 2;
        } else {
            if (lo == entries.begin()) return 0.0;
            // Case: need to interpolate
            auto hi = lo;
            --lo;
            auto lo_percentile = iterator_to_percentile(lo);
            auto hi_percentile = iterator_to_percentile(hi);
            return linear_interpolation(*hi, hi_percentile, *lo, lo_percentile, query);
        }
    }

    // destructor
    virtual ~Percentiler() = default;

    // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
    Percentiler(const Percentiler &other) = delete;

    Percentiler(Percentiler &&other) = delete;

    Percentiler &operator=(const Percentiler &other) = delete;

    Percentiler &operator=(Percentiler &&other) = delete;

private:
    const std::vector<GenotypeConfidence> &entries;
    using GC_it = std::vector<GenotypeConfidence>::const_iterator;
    GenotypeConfidence min_entry, max_entry;

    std::size_t get_rank(GC_it const &it) {
        return std::distance(entries.begin(), it) + 1;
    }

    GenotypePercentile rank_to_percentile(std::size_t const rank) {
        return 100 * rank / entries.size();
    }

    GenotypePercentile iterator_to_percentile(GC_it const &it) {
        return rank_to_percentile(get_rank(it));
    };


    static double linear_interpolation(double big_x, double big_y, double small_x, double small_y, double x) {
        // see https://en.wikipedia.org/wiki/Linear_interpolation for this formula
        return (small_y * (big_x - x) + big_y * (x - small_x)) / (big_x - small_x);
    }
};

} // end namespace GCP
#endif // GCP_H

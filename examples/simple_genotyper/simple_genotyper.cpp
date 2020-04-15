#include "Genotyper.h"
#include <map>
#include <string>
#include <iostream>

/**
 * Simple genotyper, the genotype information is given as a map,
 * and all it does is to compute max likelihood - 2nd max likelihood.
 */
using GenotypingDataMap = const std::map<std::string, double>;
class SimpleGenotyper : public Genotyper {
public:
  SimpleGenotyper() = default;
  virtual ~SimpleGenotyper() = default;
  
  virtual double get_genotype_confidence (const void* data) const {
    GenotypingDataMap* genotyping_data_map = (GenotypingDataMap*)(data);

    double max_likelihood = genotyping_data_map->at("max_likelihood");
    double second_max_likelihood = genotyping_data_map->at("second_max_likelihood");
    double genotype_confidence = max_likelihood - second_max_likelihood;
    return genotype_confidence;
  }
};


int main() {
  GenotypingDataMap genotyping_data_map_1 = {
      {"max_likelihood", 15},
      {"second_max_likelihood", 8.5}
  };
  GenotypingDataMap genotyping_data_map_2 = {
      {"max_likelihood", 20},
      {"second_max_likelihood", 19.5}
  };

  SimpleGenotyper simple_genotyper;
  std::cout << "1st gt conf = " << simple_genotyper.get_genotype_confidence(&genotyping_data_map_1) << std::endl;
  std::cout << "2nd gt conf = " << simple_genotyper.get_genotype_confidence(&genotyping_data_map_2) << std::endl;

  return 0;
}
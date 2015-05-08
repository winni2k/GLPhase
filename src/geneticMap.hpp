/* @(#)geneticMap.hpp
 */

#ifndef _GENETICMAP_HPP
#define _GENETICMAP_HPP 1

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

#include <vector>
#include <algorithm>
#include <string>
#include <exception>
#include <utility>
#include "utils.hpp"
#include <boost/algorithm/string.hpp>

class GeneticMap {

private:
  std::vector<std::pair<unsigned, double>> m_sortedMap;

public:
  GeneticMap(){};
  GeneticMap(std::string &fileName);
  double GeneticLocation(unsigned genomLoc);

  // genetic distance in Morgans
  double GeneticDistance(unsigned first, unsigned second);
  bool empty();
};

#endif /* _GENETICMAP_HPP */

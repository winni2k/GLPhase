#include "geneticMap.hpp"

using namespace std;
GeneticMap::GeneticMap(string &fileName) {

  // populate genetic map vector
  ifile gmFD(fileName);
  if (!gmFD.isGood())
    throw std::runtime_error("[GeneticMap] Could not open file: " + fileName);

  // read header
  string buffer;
  getline(gmFD, buffer, '\n');
  if (buffer != "position COMBINED_rate(cM/Mb) Genetic_Map(cM)")
    throw std::runtime_error("header of genetic map does not look right. "
                             "Expected 'position COMBINED_rate(cM/Mb) "
                             "Genetic_Map(cM)', but got '" +
                             buffer + "'");

  // read in each line
  vector<string> cols;
  unsigned rowNum = 2;
  assert(m_sortedMap.empty());
  while (getline(gmFD, buffer, '\n')) {
    cols.clear();
    boost::split(cols, buffer, boost::is_any_of(" "));
    if (cols.size() != 3)
      throw std::runtime_error("Could not find three columns in line " +
                               to_string(rowNum));

    // check the next position is greater than the previous one
    unsigned position = stoul(cols[0]);
    if (!m_sortedMap.empty() && position <= m_sortedMap.back().first)
      throw std::runtime_error(
          "Genetic map does not appear to be sorted by position at line " +
          to_string(rowNum));

    m_sortedMap.push_back(make_pair(position, stod(cols[2])));

    ++rowNum;
  }
  if (m_sortedMap.empty())
    throw std::runtime_error("Genetic map is empty.");
}

double GeneticMap::GeneticLocation(unsigned genomLoc) {

  // test assumptions
  if (m_sortedMap.empty())
    return 0;
  if (genomLoc < m_sortedMap.front().first ||
      genomLoc > m_sortedMap.back().first)
    throw std::runtime_error("Genomic location " + to_string(genomLoc) +
                             " is not on genetic map");

  // find position in map that is greater or equal to genomLoc
  auto high = std::lower_bound(
      m_sortedMap.begin(), m_sortedMap.end(), genomLoc,
      [](const pair<unsigned, double> &a, unsigned b) { return a.first < b; });

  // genomLoc is equal to high
  if (high->first == genomLoc)
    return high->second;

  // interpolate if genomLoc is less than high
  auto low = high - 1;
  double diffGenD = high->second - low->second;
  double diffGenomD = high->first - low->first;

  return low->second + (genomLoc - low->first) / diffGenomD * diffGenD;
}

double GeneticMap::GeneticDistance(unsigned first, unsigned second) {

  if (first > second)
    throw std::runtime_error("first is larger than second");

  return (GeneticLocation(second) - GeneticLocation(first)) / 100;
}

/* @(#)bcfReader.hpp
 */

#ifndef _BCFREADER_HPP
#define _BCFREADER_HPP 1

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

// necessary for use of htslib
#include <string>
#include <vector>
#include <exception>
#include <stdexcept>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <cfloat>
#include <cmath>
#include "bio.hpp"
#include <stdint.h>

namespace BCFReaderHelper {

enum class extract_t { GL, Haps };

double phred2Prob(double phred);
}

class BCFReader {
private:
  std::string m_extractString;
  std::vector<std::string> m_sampNames;
  std::vector<std::vector<double> > m_GLs;
  std::vector<std::vector<char> > m_haps;
  std::vector<Bio::snp> m_sites;
  const BCFReaderHelper::extract_t m_extractType =
      BCFReaderHelper::extract_t::GL;

  std::vector<char> ExtractRecAlleles(bcf1_t *rec, bcf_hdr_t *hdr);
  std::vector<double> ExtractRecGLs(bcf1_t *rec, bcf_hdr_t *hdr,
                                    const std::string &extractString);

public:
  explicit BCFReader(std::string fileName,
                     BCFReaderHelper::extract_t extractType) {
    BCFReader(fileName, extractType, "");
  }
  explicit BCFReader(std::string fileName,
                     BCFReaderHelper::extract_t extractType,
                     std::string region);

  std::vector<std::string> GetSampNames() const { return m_sampNames; }
  std::vector<Bio::snp> GetSites() const { return m_sites; }
  double GetSiteGL(size_t siteNum, size_t glNum) {
    return m_GLs.at(siteNum).at(glNum);
  }
};

#endif /* _BCFREADER_HPP */

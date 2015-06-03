/* @(#)glReader.hpp
 */

#ifndef _GLREADER_HPP
#define _GLREADER_HPP 1

#include <algorithm>
#include <unordered_set>
#include <string>
#include <utility>
#include <vector>
#include "utils.hpp"
#include "bio.hpp"
#include "htslibpp/vcf.hpp"
#include "htslibpp/synced_bcf_reader.hpp"

namespace Bio {

namespace GLHelper {

enum class gl_t { BCF, STBIN };
enum class gl_ret_t { STANDARD, ST_DROP_FIRST };

struct init {
  std::string glFile;
  std::string nameFile;
  std::string sampleSubsetFile;
  unsigned size_limit = 0; // 0 = no limit on number of sites to read
  gl_t glType = gl_t::STBIN;
  gl_ret_t glRetType = gl_ret_t::ST_DROP_FIRST;
  Bio::Region targetRegion;
  std::vector<std::string> glTags = {"GL", "PL"};
};
}

class GLReader {
private:
  GLHelper::init m_init;
  Bio::snp_storage_ordered m_sites;
  std::vector<std::string> m_names;
  std::vector<size_t> m_filteredNameIDXs;
  std::unordered_set<std::string> m_keepNames;
  std::vector<float> m_gls;
  size_t m_numNotRead = 0;
  const std::vector<std::string> knownGLTags = {"GL", "PL"};

  void LoadNames();
  void LoadFilteredNameIDXs();
  void LoadKeepNames();
  void LoadGLs();
  void LoadSTBinNames();
  void LoadBCFNames();
  void LoadSTBinGLs();
  void LoadBCFGLs();
  std::pair<std::vector<float>, Bio::snp_storage_ordered> GetSTBinGLs();

  std::string chooseTag(bcf_hdr_t &hdr);

  // utility functions
  template <typename INT_T>
  std::vector<float> convert_int_to_float(bcf_hdr_t &hdr,
                                          const std::string &tag,
                                          bcf1_extended<false> &rec) {

    const int numVals = 3;
    auto gls = rec.get_format_int<INT_T>(hdr, tag);
    if (gls.second != m_names.size() * numVals)
      throw std::runtime_error("Returned number of values is not correct: " +
                               std::to_string(gls.second));
    INT_T *glFirst = gls.first.get();
    std::vector<float> out_gls;
    for (size_t idx : m_filteredNameIDXs) {
      INT_T *p = glFirst + 3 * idx;
      float homR = phred2prob<float, INT_T>(*p);
      float het = phred2prob<float, INT_T>(*(p + 1));
      float homA = phred2prob<float, INT_T>(*(p + 2));

      if (m_init.glRetType != GLHelper::gl_ret_t::ST_DROP_FIRST) {
        out_gls.push_back(homR);
        out_gls.push_back(het);
        out_gls.push_back(homA);
      } else {
        float sum = homR + het + homA;
        out_gls.push_back(het / sum);
        out_gls.push_back(homA / sum);
      }
    }
    return out_gls;
  }

public:
  GLReader(){};
  GLReader(GLHelper::init init) : m_init(std::move(init)){};
  void clear() {
    m_sites.clear();
    m_names.clear();
    m_filteredNameIDXs.clear();
    m_keepNames.clear();
    m_gls.clear();
    m_numNotRead = 0;
  };

  // Setters
  void SetArgs(GLHelper::init init) {
    m_init = std::move(init);
    clear();
  }
  void SetRetGLType(GLHelper::gl_ret_t type) {
    m_init.glRetType = type;
    clear();
  }
  void SetSiteReadLimit(size_t limit) {
    m_init.size_limit = limit;
    clear();
  }
  void SetRegion(Bio::Region targetRegion) {
    m_init.targetRegion = std::move(targetRegion);
    clear();
  }
  void SetSamplesFile(std::string samplesFile) {
    m_init.sampleSubsetFile = std::move(samplesFile);
    clear();
  }

  // Getters
  std::pair<std::vector<float>, Bio::snp_storage_ordered> GetGLs();
  std::vector<std::string> GetNames();
  std::string GetGLFile() const { return m_init.glFile; };
  size_t GetNumNotRead() const { return m_numNotRead; };
};
}
#endif /* _GLREADER_HPP */

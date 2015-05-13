/* @(#)glReader.hpp
 */

#ifndef _GLREADER_HPP
#define _GLREADER_HPP 1

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
  unsigned size_limit = 0;  // 0 = no limit on number of sites to read
  gl_t glType = gl_t::STBIN;
  gl_ret_t glRetType = gl_ret_t::ST_DROP_FIRST;
  Bio::Region targetRegion;
};

}

class GLReader {
private:
  GLHelper::init m_init;
  Bio::snp_storage_ordered m_sites;
  std::vector<std::string> m_names;
  std::vector<float> m_gls;
  size_t m_numNotRead = 0;

  void LoadNames();
  void LoadGLs();
  void LoadSTBinNames();
  void LoadBCFNames();
  void LoadSTBinGLs();
  void LoadBCFGLs();
  std::pair<std::vector<float>, Bio::snp_storage_ordered> GetSTBinGLs();

public:
  GLReader(){};
  GLReader(GLHelper::init init) : m_init(std::move(init)){};
  void clear() {
    m_sites.clear();
    m_names.clear();
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

  // Getters
  std::pair<std::vector<float>, Bio::snp_storage_ordered> GetGLs();
  std::vector<std::string> GetNames();
  std::string GetGLFile() const { return m_init.glFile; };
  size_t GetNumNotRead() const { return m_numNotRead; };
};
}
#endif /* _GLREADER_HPP */

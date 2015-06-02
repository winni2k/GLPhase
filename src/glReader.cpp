#include "glReader.hpp"

using namespace std;
using namespace Bio;
using namespace htspp;

std::pair<std::vector<float>, Bio::snp_storage_ordered> GLReader::GetGLs() {
  LoadGLs();
  return make_pair(move(m_gls), move(m_sites));
}
void GLReader::LoadGLs() {
  if (m_sites.empty() || m_gls.empty()) {
    switch (m_init.glType) {
    case GLHelper::gl_t::STBIN:
      LoadSTBinGLs();
      break;
    case GLHelper::gl_t::BCF:
      LoadBCFGLs();
      break;
    }
  }
}

void GLReader::LoadNames() {
  if (!m_names.empty())
    return;
  switch (m_init.glType) {
  case GLHelper::gl_t::STBIN:
    LoadSTBinNames();
    break;
  case GLHelper::gl_t::BCF:
    LoadBCFNames();
    break;
  }
}

void GLReader::LoadSTBinNames() {

  m_names.clear();
  const string &binFile = m_init.nameFile;
  ifile binFD(binFile, false, "gz");

  // parse header
  string buffer;
  vector<string> tokens;
  if (!binFD.isGood())
    throw std::runtime_error("Error opening 'nameFile' file [" + binFile + "]");
  getline(binFD, buffer);
  boost::split(tokens, buffer, boost::is_any_of("\t"));
  if (tokens.size() < 4)
    throw std::runtime_error("Input bin file [" + binFile +
                             "] does not contain any sample information");
  m_names.reserve(tokens.size() - 3);
  for (auto it = tokens.begin() + 3; it != tokens.end(); ++it)
    m_names.push_back(std::move(*it));
}

void GLReader::LoadBCFNames() {
  m_names.clear();
  bcfFile_cpp bcf(m_init.nameFile, "r");
  bcf_hdr hdr(*bcf_hdr_read(bcf.data()));
  m_names = get_sample_names(*(hdr.data()));
}

vector<string> GLReader::GetNames() {
  LoadNames();
  return std::move(m_names);
}

std::string GLReader::chooseTag(bcf_hdr_t &hdr) {

  for (auto const &t : m_init.glTags)
    if (std::find(knownGLTags.begin(), knownGLTags.end(), t) !=
            knownGLTags.end() &&
        bcf_hdr_get_hrec(&hdr, BCF_HL_FMT, "ID", t.c_str(), NULL))
      return t;

  throw std::runtime_error("[GLReader] Could not find a known tag");
}

void GLReader::LoadBCFGLs() {
  if (m_names.empty())
    LoadBCFNames();
  m_sites.clear();
  m_gls.clear();
  m_numNotRead = 0; // reset num not read counter

  bcf_srs sr;
  sr.add_reader(m_init.glFile);
  if (!m_init.targetRegion.empty()) {
    if (sr.all_indexed())
      sr.set_region(m_init.targetRegion.AsString());
    else
      sr.set_target(m_init.targetRegion.AsString());
  }
  //  bcfFile_cpp bcf(m_init.glFile, "r");

  bcf1_extended<false> rec;
  // for now, extract GLs from GL tag
  const int numVals = 3;
  bcf_hdr_t *hdr = sr.get_header(0);

  // decide
  string tag = chooseTag(*hdr);
  while (sr.next_line() > 0) {

    // stop reading sites, but count how many not read
    if (m_init.size_limit != 0 && m_sites.size() >= m_init.size_limit) {
      ++m_numNotRead;
      continue;
    }

    rec.acquire_wrap(*(sr.get_line(0)));

    //    while (rec.bcf_read(bcf, hdr) >= 0) {

    // read in site
    vector<string> alleles = rec.alleles();
    if (alleles.size() < 2)
      throw std::runtime_error("Too few alleles in BCF record");
    if (alleles.size() > 2)
      throw std::runtime_error("More than two alleles per record are not "
                               "supported. Please break BCF into biallelics "
                               "using bcftools norm -m -");
    m_sites.push_back(snp(move(rec.chromName(*hdr)), rec.pos1(),
                          move(alleles[0]), move(alleles[1])));

    // read in format
    if (tag == "GL") {
      auto gls = rec.get_format_float(*hdr, tag);
      if (gls.second != m_names.size() * numVals)
        throw std::runtime_error("Returned number of values is not correct: " +
                                 to_string(gls.second));
      float *p = gls.first.get();
      for (size_t sampNum = 0; sampNum < m_names.size(); ++sampNum, p += 3) {
        float homR = gl2prob(*p);
        float het = gl2prob(*(p + 1));
        float homA = gl2prob(*(p + 2));

        if (m_init.glRetType != GLHelper::gl_ret_t::ST_DROP_FIRST) {
          m_gls.push_back(homR);
          m_gls.push_back(het);
          m_gls.push_back(homA);
        } else {
          float sum = homR + het + homA;
          m_gls.push_back(het / sum);
          m_gls.push_back(homA / sum);
        }
      }
    } else if (tag == "PL") {
      auto gls = rec.get_format_int32(*hdr, tag);
      if (gls.second != m_names.size() * numVals)
        throw std::runtime_error("Returned number of values is not correct: " +
                                 to_string(gls.second));
      int32_t *p = gls.first.get();
      for (size_t sampNum = 0; sampNum < m_names.size(); ++sampNum, p += 3) {
        float homR = phred2prob<float, int32_t>(*p);
        float het = phred2prob<float, int32_t>(*(p + 1));
        float homA = phred2prob<float, int32_t>(*(p + 2));

        if (m_init.glRetType != GLHelper::gl_ret_t::ST_DROP_FIRST) {
          m_gls.push_back(homR);
          m_gls.push_back(het);
          m_gls.push_back(homA);
        } else {
          float sum = homR + het + homA;
          m_gls.push_back(het / sum);
          m_gls.push_back(homA / sum);
        }
      }
    }
  }
}

void GLReader::LoadSTBinGLs() {
  if (m_names.empty())
    LoadSTBinNames();
  m_sites.clear();
  m_gls.clear();
  m_numNotRead = 0; // reset not read counter

  const string &binFile = m_init.glFile;
  ifile inputFD(binFile, false, "gz");

  // parse body
  if (!inputFD.isGood())
    throw std::runtime_error("Error reading from file [" + binFile + "]");
  size_t lineNum = 1;
  string buffer;
  vector<string> tokens;
  // discard first line
  getline(inputFD, buffer);

  auto &tr = m_init.targetRegion;
  // read input gls
  while (getline(inputFD, buffer)) {
    ++lineNum;
    boost::split(tokens, buffer, boost::is_any_of("\t"));
    if (tokens.size() != 3 + m_names.size())
      throw std::runtime_error(
          "Input line " + to_string(lineNum) +
          " does not have the correct number of columns [" +
          to_string(3 + m_names.size()) + "]");

    // skip GLs not in target region
    unsigned pos = stoul(tokens[1]);
    if (!tr.empty())
      if (!tr.chrom_eq(tokens[0]) || tr.startBP() > pos || tr.endBP() < pos)
        continue;

    // stop reading sites, but count how many not read
    if (m_init.size_limit != 0 && m_sites.size() >= m_init.size_limit) {
      ++m_numNotRead;
      continue;
    }

    // allow split on space for non snps
    string ref, alt;
    if (tokens[2].size() > 2) {
      size_t space = tokens[2].find_first_of(" ");
      if (space == string::npos)
        throw std::runtime_error("Could not parse alleles [" + tokens[2] + "]");

      ref = tokens[2].substr(0, space);
      alt = tokens[2].substr(space + 1, string::npos);
    } else {
      ref = tokens[2][0];
      alt = tokens[2][1];
    }

    // save site
    m_sites.push_back(snp(move(tokens[0]), pos, move(ref), move(alt)));

    // parse two likelihoods in each column
    size_t idx = 0;
    for (size_t sampNum = 0; sampNum != m_names.size(); ++sampNum) {
      float hetProb = stof(tokens[3 + sampNum], &idx);
      float homAltProb = stof(tokens[3 + sampNum].substr(idx));
      assert(hetProb + homAltProb <= 1);
      if (m_init.glRetType != GLHelper::gl_ret_t::ST_DROP_FIRST)
        m_gls.push_back(max(0.0f, 1 - hetProb - homAltProb));
      m_gls.push_back(hetProb);
      m_gls.push_back(homAltProb);
    }
  }
}

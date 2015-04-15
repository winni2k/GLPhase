#include "glReader.hpp"

using namespace std;
using namespace Bio;

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
      throw logic_error(
          "Tried to read BCF GL file, but this method is not implemented yet");
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
    throw logic_error(
        "Tried to read BCF GL file, but this method is not implemented yet");
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
    throw std::runtime_error("Error opening file [" + binFile + "]");
  getline(binFD, buffer);
  boost::split(tokens, buffer, boost::is_any_of("\t"));
  if (tokens.size() < 4)
    throw std::runtime_error("Input bin file [" + binFile +
                             "] does not contain any sample information");
  m_names.reserve(tokens.size() - 3);
  for (auto it = tokens.begin() + 3; it != tokens.end(); ++it)
    m_names.push_back(std::move(*it));
}

vector<string> GLReader::GetNames() {
  LoadNames();
  return std::move(m_names);
}

void GLReader::LoadSTBinGLs() {
  if (m_names.empty())
    LoadSTBinNames();
  m_sites.clear();
  m_gls.clear();

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
  while (getline(inputFD, buffer)) {
    ++lineNum;
    boost::split(tokens, buffer, boost::is_any_of("\t"));
    if (tokens.size() != 3 + m_names.size())
      throw std::runtime_error(
          "Input line " + to_string(lineNum) +
          " does not have the correct number of columns [" +
          to_string(3 + m_names.size()) + "]");

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
    m_sites.push_back(
        snp(move(tokens[0]), stoi(tokens[1]), move(ref), move(alt)));

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

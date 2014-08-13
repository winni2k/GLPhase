/*
  C++ wrapper to tabix indexer

  This code was copied from Erik Garison's tabixpp library
  commit: c2d6c12eb827805fb13db4bab20f74b212b8b6d0

*/

#include "tabix.hpp"

Tabix::Tabix(const string &file) : m_fname(file) {

  m_fp = hts_open(m_fname.c_str(), "r");
  if (!m_fp)
    throw std::runtime_error("Could not read " + m_fname);
  m_tbx = tbx_index_load(m_fname.c_str());
  if (!m_tbx)
    throw std::runtime_error("Could not load .tbi index of " + m_fname);
  m_itr = tbx_itr_queryi(m_tbx, HTS_IDX_START, 0, 0);
}

Tabix::~Tabix(void) {
  hts_close(m_fp);
  hts_itr_destroy(m_itr);
  tbx_destroy(m_tbx);
  free(m_str.s);
}

string Tabix::getHeader() {
  string header;
  while (hts_getline(m_fp, KS_SEP_LINE, &m_str) >= 0) {
    if (!m_str.l || m_str.s[0] != m_tbx->conf.meta_char)
      break;
    header += string(m_str.s) + "\n";
  }
  return header;
}

void Tabix::setRegion(const string &region) {
  tbx_itr_destroy(m_itr);
  m_itr = tbx_itr_querys(m_tbx, region.c_str());
  if (!m_itr)
    throw std::runtime_error("Could not seq to region: " + region);
}

bool Tabix::getNextLine(string &line) {

  if (tbx_itr_next(m_fp, m_tbx, m_itr, &m_str) >= 0) {
    line = m_str.s;
    return true;
  } else
    return false;
}

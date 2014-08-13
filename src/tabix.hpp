/*
  C++ wrapper to tabix indexer
  
  This code was copied from Erik Garison's tabixpp library
  commit: c2d6c12eb827805fb13db4bab20f74b212b8b6d0

*/


#ifndef __TABIX_HPP
#define __TABIX_HPP

#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/kseq.h>
#include <iostream>
#include <exception>
#include <stdexcept>


using namespace std;

class Tabix {

    htsFile *m_fp;
    hts_itr_t *m_itr;
    tbx_t * m_tbx;
    kstring_t m_str = {0,0,0};

public:

    string m_fname;

    Tabix(void) = delete;
    explicit Tabix(const string& file);
    ~Tabix(void);

    std::string getHeader();
    void setRegion(const string& region);
    bool getNextLine(string& line);

};

#endif


/*

# Tue Apr 14 10:53:44 BST 2015

This code was written by Erik Garison <erik.garrison@gmail.com>

This code comes from this commit:
https://github.com/ekg/tabixpp/tree/5a4e4a25c7c07198ea0ab1ab875c601da435587a

- Warren Kretzschmar <wkretzsch@gmail.com>
*/

#ifndef _TABIX_HPP
#define _TABIX_HPP 1

#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/kseq.h>
#include <iostream>
#include <cstring>
#include <vector>

class Tabix {

  htsFile *fn;
  tbx_t *tbx;
  hts_itr_t *iter;
  const tbx_conf_t *idxconf;
  std::string firstline;
  bool has_jumped;
  std::vector<std::string>::iterator current_chrom;

public:
  std::string filename;
  std::vector<std::string> chroms;

  Tabix(void);
  Tabix(std::string &file);
  ~Tabix(void);

  void getHeader(std::string &header);
  bool setRegion(std::string &region);
  bool getNextLine(std::string &line);
};

#endif /* _TABIX_HPP */

/* @(#)impute.hpp
 */

#ifndef _IMPUTE_H
#define _IMPUTE_H 1

#include "globals.h"
#include <gsl/gsl_rng.h>
#include <sys/time.h>
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cfloat>
#include <cstdio>
#include <cmath>
#include <vector>
#include <ctime>
//#include	<omp.h>
#include <set>
#include <zlib.h>
#include "bio.hpp"

#define MaskHA 0x80
#define MaskHB 0x40
#define MaskG 0x3f
typedef float fast;
typedef double real;

struct Site {
  std::string chr, all;
  unsigned pos;
};

class Impute {

protected:
  gsl_rng *rng;
  unsigned hn; // 2 * number of individuals = number of haplotypes
  unsigned pn; // 3 * number of sites = number of gls
  unsigned en; // 4 * number of sites = number of possible emissions
  unsigned wn; // this is the number of blocks of size 64 to save haps in
  std::vector<uint64_t> haps, hnew;
  std::vector<unsigned> hsum;
  std::vector<fast> tran, emit;

  fast pc[4][4]; // mutation matrix

  // *P is pointer to haplotype (in bits, i.e. each uint64_t contains 64 sites)
  // set1() sets the Ith bit in P to 1
  void set1(uint64_t *P, unsigned I) {
    // I >> Uint64_TShift is moving along the array according to which uint64_t
    // I is in
    // e.g. I <= 63 is first uint64_t, and so forth in blocks of 64 bits
    P[I >> WORDSHIFT] |= static_cast<uint64_t>(1) << (I & WORDMOD);
  }

  void set0(uint64_t *P, unsigned I) {
    P[I >> WORDSHIFT] &= ~(static_cast<uint64_t>(1) << (I & WORDMOD));
  }

  // test if bit I is 1
  uint64_t test(uint64_t *P, unsigned I) {
    return (P[I >> WORDSHIFT] >> (I & WORDMOD)) & static_cast<uint64_t>(1);
  }

  fast hmm(unsigned I, unsigned *P, fast S);

  void hmm_work(unsigned I, unsigned *P, fast S);

  virtual fast solve(unsigned I, unsigned &N, fast S);

  void replace(unsigned I);

  void result(void);

  // for the pointer to an array shift the site number bitwise by 6 (ie divide
  // by 64) then shift again by number of bits overlapping with 000111111, ie.
  // it keeps the unique address of I.
  // this is particularly designed for 1024 max sites... also converts it back
  // to type uint64_t from int
  virtual fast hmm_like(unsigned I, unsigned *P);
  bool load_bin(const char *F);

public:
  unsigned in; // number of samples
  unsigned mn; // number of sites

  // number of burnin, sampling iterations and folds
  static unsigned bn, sn, nn;
  static real density, conf;
  static std::vector<std::string> vcf_file;
  static std::set<std::string> male;
  static bool is_x, is_y;

  static bool gender(char *F);

  static void document(void);

  std::vector<Bio::snp> m_glSites;
  std::vector<std::string> name;
  std::vector<real> posi;
  std::vector<fast> prob;
  //    std::vector <std::uint16_t> pare;

  Impute() {
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(NULL));
  } // set default seed?
  ~Impute() { gsl_rng_free(rng); }

  // needs to be updated to work with htslib
  // unsigned load_vcf(const char *F);

  void initialize(void);

  void estimate(void);

  void save_vcf(const char *F);

  //    void save_pare(const char *F);
};

#endif /* _IMPUTE_H */

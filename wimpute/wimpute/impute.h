/* @(#)impute.h
 */

#ifndef _IMPUTE_H
#define _IMPUTE_H 1

#include	<gsl/gsl_rng.h>
#include	<sys/time.h>
#include	<algorithm>
#include	<iostream>
#include	<stdint.h>
#include	<unistd.h>
#include	<cstdlib>
#include	<cstring>
#include	<fstream>
#include	<sstream>
#include	<tabix.h>
#include	<cfloat>
#include	<cstdio>
#include	<vector>
#include	<zlib.h>
#include	<cmath>
#include	<ctime>
#include	<omp.h>
#include	<set>

#define    WordShift    6
#define    WordMod    63
#define    MaskHA    0x80
#define    MaskHB    0x40
#define    MaskG    0x3f
typedef float fast;
typedef double real;
typedef unsigned uint;
typedef uint64_t word;

using    namespace    std;

struct Site {
    string chr, all;
    uint pos;
};

class Impute {
private:
    gsl_rng *rng;
  uint hn; // 2 * number of individuals = number of haplotypes
  uint pn; // 3 * number of sites (number of transitions)
  uint en;
  uint wn; // this is the number of blocks of size 64 to save haps in
    vector<word> haps, hnew;
    vector<uint> hsum;
    vector<fast> tran, emit;

    fast pc[4][4];      // mutation matrix
    void set1(word *P, uint I) {
        P[I >> WordShift] |= static_cast<word>(1) << (I & WordMod);
    }

    void set0(word *P, uint I) {
        P[I >> WordShift] &= ~(static_cast<word>(1) << (I & WordMod));
    }

    word test(word *P, uint I) {
        return (P[I >> WordShift] >> (I & WordMod)) & static_cast<word>(1);
    }

    // for the pointer to an array shift the site number bitwise by 6 (ie divide by 64) then shift again by number of bits overlapping with 000111111, ie. it keeps the unique address of I.
    // this is particularly designed for 1024 max sites... also converts it back to type word from int
    fast hmm_like(uint I, uint *P);

    fast hmm(uint I, uint *P, fast S);

    void hmm_work(uint I, uint *P, fast S);

    fast solve(uint I, uint    &N, fast S, bool P);

    void replace(uint I);

    void result(void);

public:
  uint in; // number of samples
  uint mn; // number of sites
    static uint bn, sn, nn;
    static real density, conf;
    static vector <string> vcf_file;
    static set <string> male;
    static bool is_x, is_y;

    static bool gender(char *F);

    static void document(void);

    vector<Site> site;
    vector <string> name;
    vector<real> posi;
    vector<fast> prob;
    vector <uint16_t> pare;

    Impute() {
        rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(rng, time(NULL));
    }  // set default seed?
    ~Impute() {
        gsl_rng_free(rng);
    }

    bool load_bin(const char *F);

    uint load_vcf(const char *F);

    void initialize(void);

    void estimate(void);

    void save_vcf(const char *F);

    void save_pare(const char *F);
};

#endif /* _IMPUTE_H */


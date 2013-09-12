#include "haplotype.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

void Haplotype::Set(unsigned uSite, bool bAllele){

    assert(uSite < m_uNumAlleles);
    
    // set uSite to 1
    if(bAllele)
        m_utHap[uSite >> h_WordShift] |=  static_cast< uint64_t >(1) << (uSite & h_WordMod);
    // set uSite to 0
    else
        m_utHap[uSite >> h_WordShift] &= ~(static_cast< uint64_t >(1) <<( uSite & h_WordMod));
}

bool Haplotype::TestSite(unsigned uSite) const {

    assert(uSite < m_uNumAlleles);

//    uint64_t uWordNum = uSite >> h_WordShift;
//    return (m_utHap[uWordNum] >> (uSite & h_WordMod)) & static_cast< uint64_t >(1);
    return (m_utHap[uSite >> h_WordShift] >> (uSite & h_WordMod)) & static_cast< uint64_t >(1);
}

unsigned Haplotype::HammingDist(const Haplotype &oCompareHap) const {

    assert(oCompareHap.m_uNumAlleles == m_uNumAlleles);

    unsigned uHammingDist = 0;
    for( unsigned uWord = 0; uWord < m_utHap.size(); uWord ++){
        uint64_t uDiff = m_utHap[uWord] xor oCompareHap.GetWord(uWord);
        uHammingDist += __builtin_popcountll(uDiff);
    }

    return uHammingDist;
}

unsigned Haplotype::HammingDist(const uint64_t *upHap) const {

    unsigned uHammingDist = 0;
    for( unsigned uWord = 0; uWord < m_utHap.size(); uWord ++, upHap++){
        uint64_t uDiff = m_utHap[uWord] xor *upHap;
        uHammingDist += __builtin_popcountll(uDiff);
    }

    return uHammingDist;
}


unsigned Haplotype::HammingDist(const  uint64_t  *upHap1, const uint64_t  *upHap2) const {

    unsigned uHammingDist = 0;
    for( unsigned uWord = 0; uWord < m_utHap.size(); uWord ++, upHap1++, upHap2++){
        uint64_t uDiff = *upHap2 xor *upHap1;
        uHammingDist += __builtin_popcountll(uDiff);
    }

    return uHammingDist;
}



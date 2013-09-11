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

bool Haplotype::TestSite(unsigned uSite){

    assert(uSite < m_uNumAlleles);
        
    return (m_utHap[uSite >> h_WordShift] >> (uSite & h_WordMod)) & static_cast< uint64_t >(1);
}

unsigned Haplotype::HammingDist(Haplotype &oCompareHap){

    assert(oCompareHap.m_uNumAlleles == m_uNumAlleles);

    unsigned uHammingDist = 0;
    for( unsigned uWord; uWord < m_utHap.size(); uWord ++){
        uint64_t uDiff = m_utHap[uWord] xor oCompareHap.GetWord(uWord);
        uHammingDist += __builtin_popcountll(uDiff);
    }

    return uHammingDist;
}

unsigned Haplotype::HammingDist(uint64_t *upHap){

    unsigned uHammingDist = 0;
    for( unsigned uWord; uWord < m_utHap.size(); uWord ++, upHap++){
        uint64_t uDiff = m_utHap[uWord] xor *upHap;
        uHammingDist += __builtin_popcountll(uDiff);
    }

    return uHammingDist;
    


}
    

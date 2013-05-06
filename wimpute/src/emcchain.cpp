#include "emcchain.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

//default constructor
EMCChain::EMCChain(const fast fTemp, const fast fSelectTemp, const uint uI, const uint individualNum, const uint uID ){

    m_uChainID = uID;
    // set the random number generator
    m_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(m_rng, time(NULL));

    // set temperature and selection temperature
    setTemp(fTemp);
    m_fSelectTemp = fSelectTemp;
    
    // set variables associated with number of individuals
    assert(individualNum > 0);
    m_uLastIndividualIndex = individualNum - 1;
    m_uHapNum = 2 * individualNum;

    // set current individual
    assert(uI < individualNum);
    m_uI = uI;

    // initialize parent haplotypes
    RandomizeParents();
}    

void EMCChain::RandomizeParents () {

    for (uint j = 0; j < 4; j++) {
        do m_auParents[j] = gsl_rng_get(m_rng) % m_uHapNum;
        while (m_auParents[j] / 2 == m_uI);
    }

}

void EMCChain::setLike ( const fast fLike ){
    m_fCurr = fLike;
    m_fCrossoverProb = expf( -m_fCurr / m_fSelectTemp );
}

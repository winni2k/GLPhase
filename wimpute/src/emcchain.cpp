#include "emcchain.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

void EMCChain::EMCChain(const fast fTemp, const uint uI, const uint individualNum ) : m_uChainID(++s_uChainIDCounter){

    // set the random number generator
    m_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(m_rng, time(NULL));

    // set temperature
    setTemp(fTemp);

    // set variables associated with number of individuals
    assert(individualNum > 0);
    m_uLastIndividualIndex = individualNum - 1;
    m_uHapNum = 2 * individualNum;

    // set current individual
    assert(uI < individualNum);
    m_uI = uI;

    // initialize parent haplotypes
    randomizeParents();
}    

void EMCChain::RandomizeParents () {

    for (uint j = 0; j < 4; j++) {
        do m_uParents[j] = gsl_rng_get(m_rng) % m_uHapNum;
        while (m_uParents[j] / 2 == m_uI);
    }

}

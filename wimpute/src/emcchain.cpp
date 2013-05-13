#include "emcchain.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

//default constructor
EMCChain::EMCChain(const fast fTemp, const fast fSelectTemp, const uint uI, const uint individualNum, const uint uID ){

    m_uChainID = uID;
//    std::cerr << "EMCCchain " << m_uChainID << " allocated!" << std::endl;

    // set the random number generator
    //   m_rng = gsl_rng_alloc(gsl_rng_default);
    //    gsl_rng_set(m_rng, time(NULL));

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
}    

void EMCChain::setLike ( const fast fLike ){
    m_fCurr = fLike;
    m_fSelection = expf( m_fCurr / m_fSelectTemp );
    assert(m_fSelection > 0);
    assert(m_fSelection <= 1);
}

void EMCChain::setParent( const uint uParentIndex, const uint uParentNum ){
    assert(uParentNum < m_uHapNum);
    assert(uParentIndex < 5);
    m_auParents[uParentIndex] = uParentIndex;
}

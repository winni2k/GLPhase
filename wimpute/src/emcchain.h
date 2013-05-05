/* @(#)emcchain.h
 */

#ifndef _EMCCHAIN_H
#define _EMCCHAIN_H 1

#include	<gsl/gsl_rng.h>
#include	<sys/time.h>

// used in estimate_EMC() for now
// this might have to become a class as I go along
// Evolutionary Monte Carlo Chain
class EMCChain {

private:
    gsl_rng *m_rng; //random number generator object
    const m_uChainID; // chain ID, set at cunstruction
    static s_uChainIDCounter = 0; // counts number of chains in existence
    fast m_fTemp; // temperature of chain
    uint m_auParents[4]; // integer of Parents (not current individual)
    uint m_uI;  // current individual; 0-based
    uint m_uLastIndividualIndex; // 0-based number of individuals
    uint m_uHapNum; // no. of haps;
    
    // temp set function
    void setTemp (const fast fTemp){ m_fTemp = fTemp };
//    fast getTemp (void){ return m_fTemp };        

public:
//    static uint s_uI; // individual number of all EMCChains


    void EMCChain(const fast fTemp, const uint uI, const numIndividuals )
    : m_uChainID( ++s_uChainIDCounter );

    void randomizeParents();

}


#endif /* _EMCCHAIN_H */


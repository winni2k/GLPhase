/* @(#)emcchain.h
 */

#ifndef _EMCCHAIN_H
#define _EMCCHAIN_H 1

#include	<gsl/gsl_rng.h>
#include	<sys/time.h> // time()
#include	<ctime> // time() ?
#include        <cassert> // assert()
#include	<stdint.h>
#include	<cmath>
#include	<cstdio>
#include	<iostream>

typedef float fast;
typedef unsigned uint;

// used in estimate_EMC() for now
// this might have to become a class as I go along
// Evolutionary Monte Carlo Chain
class EMCChain {
private:
    fast m_fCurr; // current likelihood
    fast m_fCrossoverProb; //current crossover probability
    fast m_fSelectTemp; // selection temperature
    
public:
    gsl_rng *m_rng; //random number generator object
    const uint m_uChainID; // chain ID, set at cunstruction
    static uint s_uChainCounter; // counts number of chains in existence
    fast m_fTemp; // temperature of chain
    uint m_auParents[4]; // integer of Parents (not current individual)
    uint m_uI;  // current individual; 0-based
    uint m_uLastIndividualIndex; // 0-based number of individuals
    uint m_uHapNum; // no. of haps;
    fast m_fProp; // proposed likelihood
    
    // temp set function
    void setTemp (const fast fTemp){ m_fTemp = fTemp; };
//    fast getTemp (void){ return m_fTemp };        

    void setLike (const fast fLike);
    fast getLike() const { return m_fCurr; };
    fast getCrossProb() const { return m_fCrossoverProb; };
    
    // default constructor
    EMCChain(const fast fTemp, const fast fSelectTemp, const uint uI, const uint numIndividuals );
//    ~EMCChain() { --s_uChainCounter; };

    void RandomizeParents();

};


#endif /* _EMCCHAIN_H */


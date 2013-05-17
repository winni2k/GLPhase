/* @(#)wimpute.h
 */

#ifndef _WIMPUTE_H
#define _WIMPUTE_H 1
/*
  #include <string>
  #include <fstream>
  #include <iostream>
*/
#include <memory>
#include "impute.h"
#include "emcchain.h"
#include <limits>

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

class Wimpute : public Impute
{

private:
    ofstream m_ofsLogFileStream;
    gzFile m_gzLogFileStream;
    bool m_bLogIsGz;
    string m_sLogFile;
    uint m_nIteration;
    uint m_uCycles;

    // numerator and denominator of relationship matrix
    // numerator = number of accepted proposals
    // denominator = number of total proposals
    // first index is individual for which proposal was made
    // second index is individual from which was copied
    vector<vector< unsigned >> m_2dRelationshipMatNum;
    vector<vector< unsigned >> m_2dRelationshipMatDen;

    
    // Wimpute redefinition of hmm_like
    // so far it only adds logging
    virtual  fast hmm_like(uint I, uint *P) override;

    virtual fast solve(uint I, uint    &N, fast S, bool P) override;

    fast solve_EMC(uint I, uint    N, fast S, bool P);

    fast solve_AMH(uint I, uint    N, fast S, bool P);

public:

    // print out usage
    static void document(void);
    static int s_iEstimator; // see main.cpp and document for documentation
    static uint s_uParallelChains; // see main.cpp and document for documentation
    static uint s_uCycles; // see main.cpp and document for documentation
    static bool s_bIsLogging; // true if logging

    // are we logging?
    bool LogOn(void) { return s_bIsLogging; }

    // set and open log file
    void SetLog( const string &sLogFile);

    void WriteToLog( const string & tInput );
    void WriteToLog( const EMCChain & rcChain, const bool bMutate );

    bool load_bin(const char *F);

    void estimate(void);

    // EMC version of estimate()
    void estimate_EMC(void);

    // AMH version of estimate()
    void estimate_AMH(void);

    // Roulette Wheel Selection, returns index of chain selected
    int RWSelection( const vector <EMCChain> & rvcChains);

    unsigned RJSelection( const vector<unsigned> & vuRetMatNum, const vector<unsigned> & vuRetMatDen, unsigned I, unsigned hn, gsl_rng * rng);

};

#endif /* _WIMPUTE_H */


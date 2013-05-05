/* @(#)wimpute.h
 */

#ifndef _WIMPUTE_H
#define _WIMPUTE_H 1
/*
#include <string>
#include <fstream>
#include <iostream>
*/
#include "impute.h"
#include "emcchain.h"

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
    
    // Wimpute redefinition of hmm_like
    // so far it only adds logging
    virtual  fast hmm_like(uint I, uint *P) override;

    virtual fast solve(uint I, uint    &N, fast S, bool P) override;

    fast solve_EMC(const uint I, const uint    &N, const fast S, const bool P);    
    
public:

    // print out usage
    static void document(void);
    static int s_iEstimator; // see main.cpp and document for documentation
    static uint s_uParallelChains; // see main.cpp and document for documentation
    static uint s_uCycles; // see main.cpp and document for documentation
    
    // are we logging?
    bool LogOn() { return !m_sLogFile.empty(); }

    // set and open log file
    void SetLog( const string &sLogFile);

    void WriteToLog( const stringstream & tInput );

    bool load_bin(const char *F);
    
    void estimate(void);

    // EMC version of estimate()
    void estimate_EMC(void);
    
};

#endif /* _WIMPUTE_H */


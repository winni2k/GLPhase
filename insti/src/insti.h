/* @(#)insti.h
 */

#ifndef _INSTI_H
#define _INSTI_H 1
/*
  #include <string>
  #include <fstream>
  #include <iostream>
*/
#include <memory>
#include <limits>
#include <cassert>
#include <stdint.h>
#include <utility>
#include "impute.h"
#include "emcchain.h"
#include "utils.h"
#include "relationship.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

class Insti : public Impute
{

private:
    ofstream m_ofsLogFileStream;
    gzFile m_gzLogFileStream;
    bool m_bLogIsGz;
    string m_sLogFile;
    unsigned m_nIteration;
    unsigned m_uCycles;
    bool m_bUsingRefHaps = false;

    // keep track of relationship graph
    Relationship m_oRelationship;

    // reference haplotypes
    vector<uint64_t> m_vRefHaps;
    unsigned m_uNumRefHaps = 0;
    
    // Insti redefinition of hmm_like
    // so far it only adds logging
    virtual  fast hmm_like(unsigned I, unsigned *P) override;

    fast solve(unsigned I, unsigned N, fast pen, Relationship  &oRelationship);
    virtual fast solve(unsigned I, unsigned    &N, fast pen) override { cerr << I << N << pen; exit(1); }

    fast solve_EMC(unsigned I, unsigned    N, fast S);

    // returns the number of a hap that is not owned by individual I
    unsigned SampleHap( unsigned I, bool bUseRefPanel, gsl_rng * rng, unsigned hn, unsigned m_uNumRefHaps);
    
public:

    Insti()
        : m_oRelationship(Insti::s_uNumClusters, Insti::s_uClusterType) {};
    
    bool load_refPanel( string legendFile, string hapsFile );
    
    // print out usage
    static void document(void);
    static int s_iEstimator; // see main.cpp and document for documentation
    static unsigned s_uParallelChains; // see main.cpp and document for documentation
    static unsigned s_uCycles; // see main.cpp and document for documentation
    static bool s_bIsLogging; // true if logging
    static unsigned s_uNumClusters; // number of clusters to use
    static unsigned s_uClusterType; // what type of clustering
    static unsigned s_uSABurninGen; // number of simulated annealing burnin generations
    static unsigned s_uNonSABurninGen; // number of non-SA burning generations
    static unsigned s_uStartClusterGen; // 0-based generation number at which to start clustering
        
    // bool flag to keep track if we want to phase samples from ref haps only in first round
    static bool s_bKickStartFromRef;

    static string s_sLegendFile; //location of legend file
    static string s_sRefHapsFile; // location of reference haplotypes file
    
    unsigned GetNumWords() { return wn; }

    // returns allele of hap number hap at sites number site
    bool TestRefHap(uint hap, uint site){ return test(&m_vRefHaps[ hap * wn], site) == 1; }
    
    // are we logging?
    bool LogOn(void) { return s_bIsLogging; }

    // set and open log file
    void SetLog( const string &sLogFile);

    void WriteToLog( const string & tInput );
    void WriteToLog( const EMCChain & rcChain, const bool bMutate );

    bool load_bin(const char *F);

    void initialize();
    
    void estimate();

    // EMC version of estimate()
    void estimate_EMC();

    // AMH version of estimate()
    void estimate_AMH(unsigned uRelMatType);

    // Roulette Wheel Selection, returns index of chain selected
    int RWSelection( const vector <EMCChain> & rvcChains);

    unsigned RJSelection( const vector<unsigned> & vuRetMatNum, const vector<unsigned> & vuRetMatDen, unsigned I, unsigned hn, gsl_rng * rng);

    void save_relationship_graph ( string sOutputFile );

};

#endif /* _INSTI_H */


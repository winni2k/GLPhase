/* @(#)insti.hpp
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
#include <map>
#include <unordered_map>
#include "impute.hpp"
#include "emcchain.hpp"
#include "utils.hpp"
#include "relationship.hpp"
#include "version.hpp"
#include "enums.hpp"
#include "snp.hpp"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

class Insti : public Impute {

private:
    ofstream m_ofsLogFileStream;
    gzFile m_gzLogFileStream;
    bool m_bLogIsGz;
    string m_sLogFile;
    unsigned m_nIteration;
    unsigned m_uCycles;
    bool m_bUsingRefHaps = false;
    bool m_bUsingScaffold = false;

    // keep track of relationship graph
    Relationship m_oRelationship;

    // keep track of GL sites, names as unordered map of snps
    unordered_map<unsigned, snp> m_sitesUnordered;

    // name -> name index in global variable "name"
    unordered_map<string, unsigned> m_namesUnordered;
    
    // reference haplotypes
    vector<uint64_t> m_vRefHaps;
    unsigned m_uNumRefHaps = 0;

    // scaffold haplotypes
    vector<uint64_t> m_ScaffoldHaps;
    vector<unsigned> m_vuScaffoldPositions;
    vector<std::string> m_ScaffoldIDs;
    unsigned m_uNumScaffoldHaps = 0;

    // Insti redefinition of hmm_like
    // so far it only adds logging
    virtual  fast hmm_like(unsigned I, unsigned *P) override;

    fast solve(unsigned I, unsigned N, fast pen, Relationship  &oRelationship);
    virtual fast solve(unsigned I, unsigned    &N, fast pen) override
    {
        cerr << I << N << pen;
        exit(1);
    }

    fast solve_EMC(unsigned I, unsigned    N, fast S);

    // returns the number of a hap that is not owned by individual I
    unsigned SampleHap(unsigned I, bool bUseRefPanel, gsl_rng * rng, unsigned hn,
                       unsigned m_uNumRefHaps);

    // data loading
    vector<vector<char> > OpenHap(string hapFile);
    vector<snp> OpenLegend(string legendFile);
    void OpenHaps(string hapFile, vector<vector<char> > & loadHaps,
                  vector <snp> & sites);
    void OpenSample(string sampleFile, vector<string> & IDs);
    void MatchSamples(const vector<std::string> & IDs, unsigned numHaps);
    void SubsetSamples(vector<string> & loadIDs, vector<vector< char> > & loadHaps);
    void OrderSamples(vector<string> & loadIDs,
                      vector<vector<char> > & loadHaps);
    bool SwapMatch(const snp & loadSite,
                      const Site & existSite, vector<char>  & loadHap,
                       vector<char> & existHap);

public:

    unsigned GetScaffoldNumWords(){
        assert(m_ScaffoldHaps.size() > 0);
        assert(m_ScaffoldHaps.size() % m_uNumScaffoldHaps == 0);
        return (m_ScaffoldHaps.size() / m_uNumScaffoldHaps);
    }
    string GetScaffoldID(unsigned idx){
        assert(idx < m_ScaffoldIDs.size());
        return m_ScaffoldIDs[idx];
    }
    bool TestScaffoldSite(unsigned hapNum, unsigned siteNum){
        assert(hapNum * GetScaffoldNumWords() < m_ScaffoldHaps.size());
        assert(siteNum < GetScaffoldNumWords() * (WordMod + 1));
        return test(&m_ScaffoldHaps[hapNum * GetScaffoldNumWords()], siteNum);
    }

    
    Insti()
        : m_oRelationship(Insti::s_uNumClusters, Insti::s_uClusterType) {};

    // data loading
    // haps/sample
    bool LoadHapsSamp(string hapsFile, string sampleFile, PanelType panelType);

    // hap/leg/samp
    bool LoadHapLegSamp(string legendFile, string hapFile, string sampleFile,
                        PanelType panelType);

    // filter out sites that aren't in main gl data
    void FilterSites(vector<vector<char> > & loadHaps,
                     vector<snp> & loadSites, vector<vector<char> > & filtHaps,
                     vector<snp> & filtSites, PanelType panelType);

    // copy haplotypes to space in program where they actually are supposed to go
    void LoadHaps(vector<vector<char> > & inHaps,  PanelType panelType);



    void CheckPanelPrereqs(PanelType panelType);

    void Char2BitVec(const vector<vector<char> > & inHaps, double numWords, vector<uint64_t> & storeHaps){
        assert(numWords >= 0);
        Char2BitVec(inHaps, static_cast<unsigned>(numWords), storeHaps);
    }
    
    void Char2BitVec(const vector<vector<char> > & inHaps, unsigned numWords, vector<uint64_t> & storeHaps);

    void CalculateVarAfs();



    // print out usage
    static void document(void);
    static int s_iEstimator; // see main.cpp and document for documentation

    // see main.cpp and document for documentation
    static unsigned s_uParallelChains;
    static unsigned s_uCycles; // see main.cpp and document for documentation
    static bool s_bIsLogging; // true if logging
    static unsigned s_uNumClusters; // number of clusters to use
    static unsigned s_uClusterType; // what type of clustering

    // number of simulated annealing burnin generations
    static unsigned s_uSABurninGen;
    static unsigned s_uNonSABurninGen; // number of non-SA burning generations

    // 0-based generation number at which to start clustering
    static unsigned    s_uStartClusterGen;

    // bool flag to keep track if we want to phase samples from ref haps only in first round
    static bool s_bKickStartFromRef;

    static string s_sRefLegendFile; //location of sample file
    static string s_sRefHapFile; // location of reference haplotypes file

    static string s_scaffoldHapsFile; // location of scaffold haps file
    static string s_scaffoldSampleFile; // location of scaffold sample file
    static double s_scaffoldFreqCutoff; // cutoff MAF for what to fix in scaffold

    unsigned GetNumWords()
    {
        return wn;
    }

    // returns allele of hap number hap at sites number site
    bool TestRefHap(uint hap, uint site)
    {
        return test(&m_vRefHaps[ hap * wn], site) == 1;
    }

    // are we logging?
    bool LogOn(void)
    {
        return s_bIsLogging;
    }

    // set and open log file
    void SetLog(const string &sLogFile);

    void WriteToLog(const string & tInput);
    void WriteToLog(const EMCChain & rcChain, const bool bMutate);

    bool load_bin(const char *F);

    void initialize();

    void estimate();

    // EMC version of estimate()
    void estimate_EMC();

    // AMH version of estimate()
    void estimate_AMH(unsigned uRelMatType);

    // Roulette Wheel Selection, returns index of chain selected
    int RWSelection(const vector <EMCChain> & rvcChains);

    unsigned RJSelection(const vector<unsigned> & vuRetMatNum,
                         const vector<unsigned> & vuRetMatDen, unsigned I, unsigned hn, gsl_rng * rng);

    void save_relationship_graph(string sOutputFile);
    void save_vcf(const char *F);

    bool UsingScaffold()
    {
        return (m_uNumScaffoldHaps > 0);
    };
};

#endif /* _INSTI_H */




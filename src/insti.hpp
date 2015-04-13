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
#include <unordered_map>
#include "impute.hpp"
#include "emcchain.hpp"
#include "utils.hpp"
#include "sampler.hpp"
#include "version.hpp"
#include "hapPanel.hpp"
#include "MHSampler.hpp"
#include "kMedoids.hpp"
#include "kNN.hpp"
#include "glPack.hpp"
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include <boost/algorithm/string.hpp>
#include <cfloat>
#include "tabix.hpp"
#include "vcf_parser.hpp"
#include "bio.hpp"
#include "globals.h"
#include <omp.h>
#include "geneticMap.hpp"

#ifndef NCUDA
#include "hmmLike.hpp"
#endif

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
//#include <boost/spirit/include/phoenix_stl.hpp>

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

namespace InstiHelper {
struct Init {
  unsigned estimator = 0;
  std::string scaffoldHapsFile;
  std::string scaffoldSampleFile; // location of scaffold sample file
  double scaffoldFreqLB =
      0.05; // cutoff variant AF for what to cluster on in scaffold
  double scaffoldFreqUB =
      0.95; // cutoff variant AF for what to cluster in scaffold
  bool scaffoldUsingMAF = false;
  bool initPhaseFromScaffold = false;
  size_t reclusterEveryNGen = 0; // 0 means don't recluster
  size_t numThreads = 1;

  // genetic map file name
  std::string geneticMap;
};
}
// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

enum class InstiPanelType { REFERENCE, SCAFFOLD };

class Insti : public Impute {

private:
  std::ofstream m_ofsLogFileStream;
  gzFile m_gzLogFileStream;
  bool m_bLogIsGz;
  std::string m_sLogFile;
  unsigned m_nIteration;
  unsigned m_uCycles;
  bool m_bUsingRefHaps = false;

  //    bool m_bUsingScaffold = false;

  // keep track of relationship graph
  std::shared_ptr<Sampler> m_sampler;

  // keep track of GL sites, names as unordered map of snps
  std::unordered_map<std::string, Bio::snp> m_sitesUnordered;

  // name -> name index in global variable "name"
  std::unordered_map<std::string, unsigned> m_namesUnordered;

  // reference haplotypes
  std::vector<uint64_t> m_vRefHaps;
  unsigned m_uNumRefHaps = 0;

  // scaffold haplotypes
  HapPanel m_scaffold;
  size_t m_reclusterEveryNGen = 0;

  // see InstiHelper::init for default values
  const std::string m_scaffoldHapsFile;   // location of scaffold haps file
  const std::string m_scaffoldSampleFile; // location of scaffold sample file
  double m_scaffoldFreqCutoff; // cutoff MAF for what to fix in scaffold
  const bool m_initPhaseFromScaffold;
  GeneticMap m_geneticMap;
  const InstiHelper::Init m_init;

  // holds chrom, start and end position, etc.
  Bio::Region m_runRegion;

  // Insti redefinition of hmm_like
  // so far it only adds logging
  virtual fast hmm_like(unsigned I, unsigned *P) override;

  fast solve(unsigned I, unsigned N, fast pen,
             std::shared_ptr<Sampler> &sampler);
#ifndef NCUDA
  /*
    update haps sampleStride haplotypes at a time
  */
  fast cudaSolve(HMMLike &hapSampler, unsigned sampleStride, fast pen);
#endif

  virtual fast solve(unsigned I, unsigned &N, fast pen) override {
    std::cerr << I << N << pen;
    exit(1);
  }

  fast solve_EMC(unsigned I, unsigned N, fast S);

  // returns the number of a hap that is not owned by individual I
  unsigned SampleHap(unsigned I, bool bUseRefPanel, gsl_rng *rng, unsigned hn,
                     unsigned m_uNumRefHaps);

  // data loading
  std::vector<std::vector<char>> OpenHap(std::string hapFile);
  std::vector<Bio::snp> OpenLegend(std::string legendFile);
  void OpenVCFGZ(const std::string &vcf, const std::string &region,
                 std::vector<std::vector<char>> &haps,
                 std::vector<Bio::snp> &sites, std::vector<std::string> &ids);
  void OpenHaps(const std::string &hapFile,
                std::vector<std::vector<char>> &loadHaps,
                std::vector<Bio::snp> &sites);
  void OpenTabHaps(const std::string &hapFile,
                   std::vector<std::vector<char>> &loadHaps,
                   std::vector<Bio::snp> &sites);
  void OpenSample(const std::string &sampleFile,
                  std::vector<std::string> &fillSampleIDs);
  void MatchSamples(const std::vector<std::string> &IDs, unsigned numHaps);
  void SubsetSamples(std::vector<std::string> &loadIDs,
                     std::vector<std::vector<char>> &loadHaps);
  void OrderSamples(std::vector<std::string> &loadIDs,
                    std::vector<std::vector<char>> &loadHaps);
  bool SwapMatch(const Bio::snp &loadSite, const Site &existSite,
                 std::vector<char> &loadHap, std::vector<char> &existHap);

  // expects scaffold to have been initialized
  void SetHapsAccordingToScaffold();

public:
  // uuid
  const boost::uuids::uuid m_tag;

  Insti() = delete;

  Insti(InstiHelper::Init &init);

  unsigned GetScaffoldNumWordsPerHap() { return m_scaffold.NumWordsPerHap(); }
  std::string GetScaffoldID(unsigned idx) { return m_scaffold.GetID(idx); }
  bool TestScaffoldSite(unsigned hapNum, unsigned siteNum) {
    assert(hapNum < m_scaffold.NumHaps());
    assert(siteNum < m_scaffold.MaxSites());
    std::vector<uint64_t> *scaffoldHaps = m_scaffold.Haplotypes();
    uint64_t *scaffoldHapsPointer = scaffoldHaps->data();
    return test(&scaffoldHapsPointer[hapNum * m_scaffold.NumWordsPerHap()],
                siteNum);
  }
  unsigned GetScaffoldNumHaps() { return m_scaffold.NumHaps(); }
  unsigned GetScaffoldNumSites() { return m_scaffold.NumSites(); }

  // data loading
  // haps/sample
  void LoadHapsSamp(const std::string &hapsFile, const std::string &sampleFile,
                    InstiPanelType panelType);

  // hap/leg/samp
  bool LoadHapLegSamp(const std::string &legendFile, const std::string &hapFile,
                      const std::string &sampleFile, InstiPanelType panelType);

  // vcf.gz file
  void LoadVCFGZ(const std::string &vcf, InstiPanelType panel_t,
                 const std::string &region);

  // filter out sites that aren't in main gl data
  void FilterSites(std::vector<std::vector<char>> &loadHaps,
                   std::vector<Bio::snp> &loadSites,
                   std::vector<std::vector<char>> &filtHaps,
                   std::vector<Bio::snp> &filtSites, InstiPanelType panelType);

  // copy haplotypes to space in program where they actually are supposed to
  // go
  // along with list of sites if applicable
  void LoadHaps(std::vector<std::vector<char>> &inHaps,
                std::vector<Bio::snp> &inSites,
                std::vector<std::string> &inSampleIDs,
                InstiPanelType panelType);

  void CheckPanelPrereqs(InstiPanelType panelType);

  std::vector<uint64_t>
  Char2BitVec(const std::vector<std::vector<char>> &inHaps, double numWords) {
    assert(numWords >= 0);
    return Char2BitVec(inHaps, static_cast<unsigned>(numWords));
  }

  std::vector<uint64_t>
  Char2BitVec(const std::vector<std::vector<char>> &inHaps, unsigned numWords);

  void CalculateVarAfs();

  // Metropolis Hastings with Annealing is default
  unsigned m_estimator = 0;

  // see main.cpp and document for documentation
  static unsigned s_uParallelChains;
  static unsigned s_uCycles;      // see main.cpp and document for documentation
  static bool s_bIsLogging;       // true if logging
  static unsigned s_uNumClusters; // number of clusters to use
  static unsigned s_uClusterType; // what type of clustering
  static kNNDistT s_clusterDistanceMetric; // what type of clustering

  // number of simulated annealing burnin generations
  static unsigned s_uSABurninGen;
  static unsigned s_uNonSABurninGen; // number of non-SA burning generations

  // 0-based generation number at which to start clustering
  static unsigned s_uStartClusterGen;

  // bool flag to keep track if we want to phase samples from ref haps only in
  // first round
  static bool s_bKickStartFromRef;

  static std::string s_sRefLegendFile; // location of sample file
  static std::string s_sRefHapFile;    // location of reference haplotypes file

  static MHType s_MHSamplerType;

  // print out usage
  static void document();

  unsigned GetNumWords() { return wn; }

  // returns allele of hap number hap at sites number site
  bool TestRefHap(uint hap, uint site) {
    return test(&m_vRefHaps[hap * wn], site) == 1;
  }

  // test if bit I is 1 in hap numebr hapNum
  // provided for unit testing only
  uint64_t TestMainHap_(unsigned hapNum, unsigned I) {
    uint64_t *ha = &haps[hapNum * wn];
    return test(ha, I);
  }

  // are we logging?
  bool LogOn(void) { return s_bIsLogging; }

  // set and open log file
  void SetLog(const std::string &sLogFile);

  void WriteToLog(const std::string &tInput);
  void WriteToLog(const EMCChain &rcChain, const bool bMutate);

  void load_bin(const std::string &binFile);

  void initialize();

  void estimate();

  // EMC version of estimate()
  void estimate_EMC();

  // AMH version of estimate()
  void estimate_AMH();

  // Roulette Wheel Selection, returns index of chain selected
  int RWSelection(const std::vector<EMCChain> &rvcChains);

  unsigned RJSelection(const std::vector<unsigned> &vuRetMatNum,
                       const std::vector<unsigned> &vuRetMatDen, unsigned I,
                       unsigned hn, gsl_rng *rng);

  void save_relationship_graph(std::string sOutputFile);
  void save_vcf(const char *F, std::string commandLine);

  bool UsingScaffold() { return (m_scaffold.Initialized()); };
};

#endif /* _INSTI_H */

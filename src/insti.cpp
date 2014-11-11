#include "insti.hpp"

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;
#ifndef NCUDA
using namespace HMMLikeCUDA;
#endif
using namespace Bio;
namespace HPH = HapPanelHelper;

/*
  #define DEBUG
  #define DEBUG2
  #define DEBUG3
*/

#ifdef DEBUG
#define DEBUG_MSG(str)                                                         \
  do {                                                                         \
    std::cerr << str;                                                          \
  } while (false)
#else
#define DEBUG_MSG(str)                                                         \
  do {                                                                         \
  } while (false)
#endif

#ifdef DEBUG2
#define DEBUG_MSG2(str)                                                        \
  do {                                                                         \
    std::cerr << str;                                                          \
  } while (false)
#else
#define DEBUG_MSG2(str)                                                        \
  do {                                                                         \
  } while (false)
#endif

#ifdef DEBUG3
#define DEBUG_MSG3(str)                                                        \
  do {                                                                         \
    std::cerr << str;                                                          \
  } while (false)
#else
#define DEBUG_MSG3(str)                                                        \
  do {                                                                         \
  } while (false)
#endif

// initializyng static member variables

// number of parallel chains to use for parallel estimators
unsigned Insti::s_uParallelChains = 5;
unsigned Insti::s_uCycles = 0;
bool Insti::s_bIsLogging = false;
bool Insti::s_bKickStartFromRef = false;
unsigned Insti::s_uNumClusters = 0;
unsigned Insti::s_uClusterType = 0;
kNNDistT Insti::s_clusterDistanceMetric = kNNDistT::hamming;
unsigned Insti::s_uSABurninGen = 28;
unsigned Insti::s_uNonSABurninGen = 28;
MHType Insti::s_MHSamplerType = MHType::MH;

// start clustering after simulated annealing burnin
// need to reset this in main.cpp -- should really move to a better option
// handling approach...
unsigned Insti::s_uStartClusterGen = Insti::s_uSABurninGen;

namespace InstiHelper {
unordered_set<string> LoadSamplesUOSet(const string &sampleFile) {

  std::unordered_set<std::string> inputSamples;

  if (sampleFile.empty())
    return inputSamples;

  ifile samplesFD(sampleFile);

  if (!samplesFD.isGood())
    throw std::runtime_error("[LoadSamplesUOSet] Could not open file: [" +
                             sampleFile + "]");

  string buffer;
  while (getline(samplesFD, buffer, '\n')) {
    auto ret = inputSamples.insert(buffer);
    if (ret.second == false)
      throw runtime_error(
          "[LoadSamplesUOSet] Error while loading samples from file [" +
          sampleFile + "]: The sample " + *ret.first +
          " was encountered twice");
  }
  return inputSamples;
}
}

Insti::Insti(InstiHelper::Init init)
    : m_geneticMap(init.geneticMap), m_init{ std::move(init) },
      m_tag(boost::uuids::random_generator()()) {

  omp_set_num_threads(init.numThreads);
  // bounds check estimator input
  if (init.estimator < 4)
    m_estimator = init.estimator;
  else
    throw std::runtime_error("[insti] Estimator needs to be less than 4");

  // load keep samples
  auto keepSamples = InstiHelper::LoadSamplesUOSet(m_init.inputSamplesKeep);

  // load input GLs according to file type
  LoadGLs(keepSamples);

  // setting range object based on load bin
  if (m_glSites.empty())
    throw std::runtime_error(
        "[insti] Loaded data seems to contain no sites in file: " +
        string(init.inputGLFile));
  m_runRegion =
      Region(m_glSites[0].chr, m_glSites[0].pos, m_glSites.back().pos);

  // setting number of cycles to use
  // here is best place to do it because in is defined in load_bin()
  if (s_uCycles > 0)
    m_uCycles = s_uCycles;
  else {
    m_uCycles = nn * in; // this was how snptools does it
  }
}

// return the probability of the model given the input haplotypes P and
// emission and transition matrices of individual I
// call Impute::hmm_like and print out result
fast Insti::hmm_like(unsigned I, uint *P) { return Impute::hmm_like(I, P); }

// take an individual number I, a set of four haplotype indices P, and
// penalty S to update haplotypes of individual I
// All our tests are with the penalty set to 1
void Insti::hmm_work(unsigned I, unsigned *P, fast S) {

  // set up pointers to the beginning of four hapltoypes
  // m0 is "mother haplotype 0", m1 is "mother haplotype 1"
  // f0 is "father haplotype 0", f1 is "father haplotype 1"

  // wn stands for "word number".  This is the number of unsigned
  // integers that a haplotype is stored in

  uint64_t *f0 = &haps[P[0] * wn], *f1 = &haps[P[1] * wn],
           *m0 = &haps[P[2] * wn], *m1 = &haps[P[3] * wn];

  //	backward recursion
  // mn is the number of sites in the chunk being phased
  // fast is a typedef for float (not my choice)
  vector<fast> beta(mn * 4);

  // create pointers that point to last set of elements of emit, tran and beta
  // en is the number of emission matrix elements per sample
  fast *e = &emit[(I + 1) * en - 4], *t = &tran[(mn - 1) * 3], sum,
       *b = &beta[(mn - 1) * 4];
  fast l00 = 0, l01 = 0, l10 = 0, l11 = 0;
  fast b00, b01, b10, b11;

  // initial state of backward sampler
  b[0] = b[1] = b[2] = b[3] = 1;

  // fill b with the forward probabilites
  // test(f0, m) returns 1 if haplotype f0 is 1 at site m and 0 otherwise
  for (unsigned m = mn - 1; m; m--, e -= 4, t -= 3) {

    b00 = b[0] * e[(test(f0, m) << 1) | test(m0, m)];
    b01 = b[1] * e[(test(f0, m) << 1) | test(m1, m)];
    b10 = b[2] * e[(test(f1, m) << 1) | test(m0, m)];
    b11 = b[3] * e[(test(f1, m) << 1) | test(m1, m)];
    b -= 4;
    b[0] = b00 * t[0] + (b01 + b10) * t[1] + b11 * t[2];
    b[1] = b01 * t[0] + (b00 + b11) * t[1] + b10 * t[2];
    b[2] = b10 * t[0] + (b00 + b11) * t[1] + b01 * t[2];
    b[3] = b11 * t[0] + (b01 + b10) * t[1] + b00 * t[2];
    sum = (b[0] + b[1] + b[2] + b[3]);
    if (sum == 0)
      throw runtime_error(
          "Chosen haplotypes are incompatible with emission matrix");
    sum = 1.0f / sum;
    b[0] *= sum;
    b[1] *= sum;
    b[2] *= sum;
    b[3] *= sum;
  }

  //	forward sampling
  // walk through b
  uint64_t *ha = nullptr;
  if (m_init.serializeHapUpdate)
    ha = &haps[I * 2 * wn];
  else
    ha = &hnew[I * 2 * wn];
  uint64_t *hb = ha + wn;

  fast *p = &prob[I * pn];
  e = &emit[I * en];
  t = &tran[0];
  b = &beta[0];
  unsigned s00, s01, s10, s11;
  for (unsigned m = 0; m < mn; m++, e += 4, t += 3, p += 3, b += 4) {
    s00 = (test(f0, m) << 1) | test(m0, m);
    s01 = (test(f0, m) << 1) | test(m1, m);
    s10 = (test(f1, m) << 1) | test(m0, m);
    s11 = (test(f1, m) << 1) | test(m1, m);
    if (m) {
      b00 = l00 * t[0] + l01 * t[1] + l10 * t[1] + l11 * t[2],
      b01 = l00 * t[1] + l01 * t[0] + l10 * t[2] + l11 * t[1];
      b10 = l00 * t[1] + l01 * t[2] + l10 * t[0] + l11 * t[1],
      b11 = l00 * t[2] + l01 * t[1] + l10 * t[1] + l11 * t[0];
      l00 = b00 * e[s00];
      l01 = b01 * e[s01];
      l10 = b10 * e[s10];
      l11 = b11 * e[s11];
    } else {
      l00 = 0.25f * e[s00];
      l01 = 0.25f * e[s01];
      l10 = 0.25f * e[s10];
      l11 = 0.25f * e[s11];
    }
    sum = 1.0f / (l00 + l01 + l10 + l11);
    l00 *= sum;
    l01 *= sum;
    l10 *= sum;
    l11 *= sum;

    if (m_sitesAtWhichPhaseIsFixed.empty() ||
        test(m_sitesAtWhichPhaseIsFixed.data(), m) == 0) {

      fast p00 = l00 * b[0], p01 = l01 * b[1], p10 = l10 * b[2],
           p11 = l11 * b[3];

      // powf effectively inflates the importance of small numbers
      // while S is < 1 (first bn/2 iterations)
      fast c00 = powf(p[0] * (p00 * pc[s00][0] + p01 * pc[s01][0] +
                              p10 * pc[s10][0] + p11 * pc[s11][0]),
                      S);
      fast c01 = powf(p[1] * (p00 * pc[s00][1] + p01 * pc[s01][1] +
                              p10 * pc[s10][1] + p11 * pc[s11][1]),
                      S);
      fast c10 = powf(p[1] * (p00 * pc[s00][2] + p01 * pc[s01][2] +
                              p10 * pc[s10][2] + p11 * pc[s11][2]),
                      S);
      fast c11 = powf(p[2] * (p00 * pc[s00][3] + p01 * pc[s01][3] +
                              p10 * pc[s10][3] + p11 * pc[s11][3]),
                      S);

      // randomly choose new haplotypes at this site weighted by c
      sum = gsl_rng_uniform(rng) * (c00 + c01 + c10 + c11);
      if (sum < c00) {
        set0(ha, m);
        set0(hb, m);
      } else if (sum < c00 + c01) {
        set0(ha, m);
        set1(hb, m);
      } else if (sum < c00 + c01 + c10) {
        set1(ha, m);
        set0(hb, m);
      } else {
        set1(ha, m);
        set1(hb, m);
      }
    }
    // just copy over phase from current haplotype estimates for this individual
    else {
      if (!m_init.serializeHapUpdate) {
        if (test(&haps[I * 2 * wn], m))
          set1(ha, m);
        else
          set0(ha, m);
        if (test(&haps[I * 2 * wn + wn], m))
          set1(hb, m);
        else
          set0(hb, m);
      }
    }
  }
}

void Insti::SetLog(const string &sLogFile) {

  s_bIsLogging = true;
  m_sLogFile = sLogFile;
  m_bLogIsGz = false;
  size_t found = sLogFile.find(".gz", sLogFile.length() - 3);

  if (found != std::string::npos) {
    m_bLogIsGz = 1;
    m_gzLogFileStream = gzopen(m_sLogFile.c_str(), "wt");
  } else {

    // open logFileStream for append if it's not yet open
    if (!m_ofsLogFileStream)
      m_ofsLogFileStream.open(m_sLogFile);

    // otherwise close and open again
    else {
      m_ofsLogFileStream.close();
      m_ofsLogFileStream.open(m_sLogFile);
    }

    // exit if the file cannot be opened
    if (!m_ofsLogFileStream.is_open()) {
      cerr << "could not open log file " << m_sLogFile << " for writing"
           << endl;
      exit(1);
    }
  }

  cout << "Logging to:\t" << m_sLogFile << endl;
};

void Insti::WriteToLog(const string &tInput) {
  if (!s_bIsLogging)
    return;

  if (m_bLogIsGz)
    gzprintf(m_gzLogFileStream, tInput.c_str());
  else {

    // open logFileStream for append if it's not yet open
    if (!m_ofsLogFileStream)
      m_ofsLogFileStream.open(m_sLogFile);

    // exit if the file cannot be opened
    if (!m_ofsLogFileStream.is_open()) {
      cerr << "could not open log file " << m_sLogFile << " for writing"
           << endl;
      exit(1);
    }

    // write input to log file
    m_ofsLogFileStream << tInput;
    m_ofsLogFileStream.flush();
  }

  DEBUG_MSG("wrote something" << endl);
};

// what to log when given an EMCChain
void Insti::WriteToLog(const EMCChain &rcChain, const bool bMutate) {
  stringstream message;
  message << m_nIteration << "\t" << rcChain.m_uI << "\t" << rcChain.getLike()
          << "\t" << rcChain.m_uChainID << "\t" << rcChain.m_fTemp << "\t"
          << bMutate << endl;
  WriteToLog(message.str());
}

// Roulette Wheel Selection, returns index of chain selected
int Insti::RWSelection(const vector<EMCChain> &rvcChains) {
  double dTotalProb = 0; // always positive, could  be larger than 1

  for (const auto &icChain : rvcChains) {
    dTotalProb += icChain.getSelection();
    DEBUG_MSG2("\ttotal prob:\t" << dTotalProb << endl);
  }

  assert(dTotalProb > 0);
  assert(dTotalProb < std::numeric_limits<double>::max());
  double dStopPoint = gsl_rng_uniform(rng) * dTotalProb;
  assert(dStopPoint > 0);
  DEBUG_MSG2("\t..." << endl);
  DEBUG_MSG2("\ttotalProb:\t" << dTotalProb << endl);
  int iChainIndex = 0;

  do
    dStopPoint -= rvcChains[iChainIndex].getSelection();

  while (dStopPoint > 0 && ++iChainIndex);

  assert(iChainIndex >= 0);
  assert(iChainIndex < static_cast<int>(rvcChains.size()));
  return iChainIndex;
}

void Insti::LoadGLs(const unordered_set<string> &keepSamples) {

  if (m_init.inputGLFileType == "bin")
    LoadGLBin();
  else if (m_init.inputGLFileType == "bcf")
    LoadGLBCF(keepSamples);
  else
    throw std::runtime_error("Unknown GL file type specified: " +
                             m_init.inputGLFileType);
}

void Insti::LoadGLBin() {
  if (m_init.inputGLFile.substr(m_init.inputGLFile.size() - 4, 4) != ".bin")
    throw runtime_error("[insti] Expecting input GL file to be in .bin format "
                        "and have file name ending .bin [" +
                        m_init.inputGLFile + "]");
  bool bRetVal = Impute::load_bin(m_init.inputGLFile.c_str());
  if (bRetVal == false)
    throw std::runtime_error("[insti] Unable to load input .bin file: " +
                             m_init.inputGLFile);
}

void Insti::LoadGLBCF(const unordered_set<string> &keepSamples) {

  string bcfFile = m_init.inputGLFile;

  clog << m_tag << ": Reading " + bcfFile << "\n";

  // first clear all the data that will be filled
  name.clear();
  m_glSites.clear();
  prob.clear();
  posi.clear();

  BCFReader bcf(bcfFile, BCFReaderHelper::extract_t::GL, m_init.inputGLRegion);
  vector<string> glSampNames = bcf.GetSampNames();

  // subset the samples according to subset list
  vector<char> keepSamp;
  keepSamp.reserve(glSampNames.size());
  name.reserve(glSampNames.size());
  unsigned numSampsKept = 0;
  for (auto samp : glSampNames) {
    auto ret = keepSamples.find(samp);
    if (!keepSamples.empty() && ret == keepSamples.end())
      keepSamp.push_back(0);
    else {
      name.push_back(samp);
      keepSamp.push_back(1);
      ++numSampsKept;
    }
  }
  if (numSampsKept != keepSamples.size() && !keepSamples.empty())
    throw runtime_error("[Insti::LoadGLBCF] Number of samples loaded is not "
                        "equal to the number of samples in keep file");
  assert(name.size() == numSampsKept);

  m_glSites = bcf.GetSites();
  mn = m_glSites.size(); // mn = number of sites

  // cutting away up to 5% of markers if region contains too many sites
  if (mn > NUMSITES) {
    if (mn > ceil(NUMSITES * 1.05))
      throw std::runtime_error(
          "Region " + m_init.inputGLRegion + " of the GL file (" +
          m_init.inputGLFile +
          ") contains more sites than the maximum chunk size " +
          to_string(NUMSITES) + ". If the extra sites were dropped, then more "
                                "than 5% of the chunk size would be dropped.");
    else {
      clog << m_tag
           << ": [Insti] Chunk contains too many sites.  Dropping the last " +
                  to_string(mn - NUMSITES) + " sites.\n";
      m_glSites.erase(m_glSites.begin() + NUMSITES, m_glSites.end());
      mn = NUMSITES;
    }
  }

  posi.reserve(mn);

  // extract GLs of interest
  prob.reserve(mn * numSampsKept * 2);
  for (size_t siteNum = 0; siteNum < mn; ++siteNum) {
    posi.push_back(m_glSites[siteNum].pos);
    for (size_t sampNum = 0; sampNum < keepSamp.size(); ++sampNum) {
      if (keepSamp[sampNum] == 1) {
        // save the het GL
        prob.push_back(bcf.GetSiteGL(siteNum, 3 * sampNum + 1));
        // save the hom alt GL
        prob.push_back(bcf.GetSiteGL(siteNum, 3 * sampNum + 2));
      }
    }
  }

  // clean up
  in = name.size();
  clog << m_tag << ": [Insti] sites\t" << mn << "\n";
  clog << m_tag << ": [Insti] sample\t" << in << "\n"; // which sample
}

void Insti::LoadVCFGZ(string vcf, InstiPanelType panel_t,
                      const Bio::Region &region) {

  cout << "[insti] Loading haplotypes from VCF: " << vcf << endl;
  if (vcf.size() < 7 || vcf.compare(vcf.size() - 7, 7, ".vcf.gz") != 0)
    throw std::runtime_error(
        "[insti] input VCF file does not end in '.vcf.gz': " + vcf);

  HapPanelHelper::Init init;
  init.hapFileType = HapPanelHelper::HapFileType::VCF;
  init.keepRegion = region;
  vector<string> inFiles{ vcf };
  LoadHaps(move(inFiles), move(init), panel_t);
};

void Insti::LoadHapsSamp(string hapsFile, string sampleFile,
                         InstiPanelType panelType, Bio::Region loadRegion) {

  // make sure both files are defined
  if (sampleFile.empty())
    throw std::runtime_error("[insti] Error while loading scaffold haps/sample "
                             "files: Need to define a sample file");

  if (hapsFile.empty())
    throw std::runtime_error("[insti] Need to define a haps file");

  // load haps file
  CheckPanelPrereqs(panelType);

  HapPanelHelper::Init init;
  init.keepRegion = std::move(loadRegion);
  // read the haps and sites from a haps file
  if (hapsFile.size() > 11 &&
      hapsFile.compare(hapsFile.size() - 11, 11, ".tabhaps.gz") == 0)
    init.hapFileType = HPH::HapFileType::WTCCC_TBX;
  else
    init.hapFileType = HPH::HapFileType::WTCCC;

  vector<string> inFiles{ hapsFile, sampleFile };
  LoadHaps(move(inFiles), move(init), panelType);
}

// swap two haps if they match and return true on swap, false otherwise
bool Insti::SwapMatch(const snp &loadSite, const Site &existSite,
                      vector<char> &loadHap, vector<char> &existHap) {
  if (loadSite.pos == existSite.pos && loadSite.chr == existSite.chr &&
      loadSite.ref + loadSite.alt == existSite.all) {
    std::swap(loadHap, existHap);

    return 1;
  }

  return 0;
}

// put the haplotypes in the right place in the program structure
void Insti::LoadHaps(vector<string> inFiles, HapPanelHelper::Init init,
                     InstiPanelType panelType) {

  CheckPanelPrereqs(panelType);

  try {

    // store the haplotypes in the correct place based on what type of panel we
    // are loading
    if (panelType == InstiPanelType::REFERENCE) {
      init.matchAllSites = true;
      init.keepSites = m_glSites;
      HapPanel input(inFiles, init);
      vector<uint64_t> storeHaps = input.Haplotypes_uint64_t(GetNumWords());
      storeHaps.swap(m_vRefHaps);
      m_uNumRefHaps = input.NumHaps();
      cout << "Reference panel haplotypes\t" << m_uNumRefHaps << endl;
      return;
    } else if (panelType == InstiPanelType::SCAFFOLD) {
      static_assert(WORDMOD >= 0,
                    "WORDMOD global is not greater or equal to 0");
      init.keepAllSitesInRegion = true;
      init.keepSampleIDs = name;
      init.allowReorderOfSitesAtPos = true;
      m_scaffold = HapPanel(inFiles, std::move(init));

      cout << "Scaffold haplotypes\t" << m_scaffold.NumHaps() << endl;

      if (m_scaffold.NumHaps() != hn)
        throw runtime_error(
            "Error while reading scaffold: Scaffold needs to have two "
            "haplotypes for every input sample");

      return;
    } else
      assert(false); // this should not happen
  }
  catch (exception &e) {

    cerr << "Error loading files [";
    for (auto f : inFiles)
      cerr << " " << f;
    cerr << " ]" << endl << e.what() << endl;
    exit(2);
  }
}

void Insti::CheckPanelPrereqs(InstiPanelType panelType) {
  switch (panelType) {
  case InstiPanelType::REFERENCE:
    m_bUsingRefHaps = true;
    assert(m_vRefHaps.size() == 0);
    break;

  case InstiPanelType::SCAFFOLD:
    assert(m_scaffold.empty());
    break;

  default:
    assert(false); // this should not happen
  }
}

void Insti::LoadRefPanel(string commaSeparatedRefPanelFiles, string fileType,
                         Bio::Region loadRegion) {

  // make sure required files are defined
  vector<string> inFiles;
  boost::split(inFiles, commaSeparatedRefPanelFiles, boost::is_any_of(","));

  HapPanelHelper::Init init;

  if (inFiles.size() == 0)
    throw runtime_error("Could not read input file [" +
                        commaSeparatedRefPanelFiles + "]");

  if (fileType == "VCF") {
    if (inFiles.size() != 1)
      throw runtime_error("Unexpected number of files given in input: [" +
                          commaSeparatedRefPanelFiles + "]");
    init.hapFileType = HPH::HapFileType::VCF;
  } else if (fileType == "WTCCC") {
    if (inFiles.size() != 2)
      throw runtime_error("Number of input files given is not 2: [" +
                          commaSeparatedRefPanelFiles + "]");

    if (inFiles[0].size() > 11 &&
        inFiles[0].compare(inFiles[0].size() - 11, 11, ".tabhaps.gz") == 0)
      init.hapFileType = HPH::HapFileType::WTCCC_TBX;
    else
      init.hapFileType = HPH::HapFileType::WTCCC;
  } else if (fileType == "IMPUTE2") {
    if (inFiles.size() != 3)
      throw runtime_error("Number of input files given is not 3: [" +
                          commaSeparatedRefPanelFiles + "]");

    init.hapFileType = HPH::HapFileType::IMPUTE2;
  } else
    throw runtime_error("Unknown ref panel input file type: " + fileType);

  init.keepRegion = std::move(loadRegion);
  LoadHaps(move(inFiles), move(init), InstiPanelType::REFERENCE);
}

/* CHANGES from impute.cpp
   support for reference haplotypes
*/

void Insti::initialize() {

  // copied from SNPTools Impute
  // all haplotypes are saved in 64 bit unsigned ints (a word), where each bit
  // represents a position
  // first, figure out how many words we'll need to store a hap and save in wn
  // if total sites overlaps with 00111111,  then wn = mn number of sites
  // shifter to right....
  // we define a minimum block size of 64.
  wn = (mn & WORDMOD) ? (mn >> WORDSHIFT) + 1 : (mn >> WORDSHIFT);
  hn = in * 2;          // number of haps
  haps.resize(hn * wn); // space to store all haplotypes
  if (!m_init.serializeHapUpdate)
    hnew.resize(hn * wn);  // number of haplotypes = 2 * number of samples  ...
                           // haps mn is # of sites,
  hsum.assign(hn * mn, 0); // one unsigned for every hap's site - what for?  To
                           // estimate allele probs

  //    pare.assign(in * in, 0);  // in x in matrix, one uint16 for every pair
  // of individuals
  pn = 3 *
       mn; // set the number of transitions.  three transitions for every site
  tran.resize(pn); // tran looks like the transition matrix,

  // i.e. recombination rate
  // transitions between 3 types of genotypes P(RR), P(RA) P(AA)
  vector<fast> temp(in * pn);

  // initialize emission matrix
  // 4 emissions, for each site 4 emissions * number of samples.  (0|0, 1|1 0|1
  // 1|0)
  en = 4 * mn;
  emit.resize(in * en);

  // is_par defines whether a site is in the paralogous region
  // if the site is on chromosome X
  // these magic numbers should probably get their own #define statement...
  // move away from vector of bools to vector of chars
  vector<bool> is_par(mn);

  if (posi.size() == mn) {
    for (unsigned m = 0; m < mn; m++)
      is_par[m] = (posi[m] >= 60001 && posi[m] <= 2699520) ||
                  (posi[m] >= 154931044 && posi[m] <= 155270560);
  }

  if (posi.size() != mn) {
    posi.resize(mn);

    for (unsigned m = 0; m < mn; m++)
      posi[m] = m;
  } // if sites not stored

  // initialize the mutation rate mu:
  // If S is a harmonic series of length hn (number of haplotypes),
  // then mu = 1/S ( hn + 1/S)
  // initialize recombination rate rho based on SNP density
  fast mu = 0, rho;

  for (unsigned i = 1; i < hn; i++)
    mu += 1.0 / i;

  // mu * (mn-1) looks like the Watterson estimator of the population mutation
  // rate theta
  // 0.5 / (posi[mn - 1] - posi[0]) / density looks like a correction for finite
  // sites
  mu = 1 / mu;
  rho = 0.5 * mu * (mn - 1) / (posi[mn - 1] - posi[0]) / density;
  mu = mu / (hn + mu); // rho is recombination rate?  mu is mutation rate

  // initialzie the site transition matrix tran
  // posi is recombination between its position and previous
  // r < 1; the larger the number of haplotypes, the smaller r gets
  // tran is a site's recombination probability matrix
  // r therefore must be a recombination rate estimate
  for (unsigned m = mn - 1; m; m--) {

    // genetic distance is in centiMorgans
    float r = 0;
    try {
      r = m_geneticMap.GeneticDistance(posi[m - 1], posi[m]) / 100;
    }
    // use old estimate of recombination rate if map does not exist
    catch (GeneticMapHelper::GenomPosOutsideMap &e) {
      clog << m_tag << ": [GeneticMap] " << e.what() << "\n";
      // replaced ad-hoc genetic distance estimate by genetic map
      unsigned genomDist = (posi[m] - posi[m - 1]) * rho;
      r = genomDist / (genomDist + hn);
    }
    catch (...) {
      throw;
    }

    r = r <= 1 ? r : 1;

    //  4 state HMM with three transitions at each position
    // for each position, transition.  r= recombination,
    // 1-r= no recombination
    tran[m * 3] = (1 - r) * (1 - r);
    tran[m * 3 + 1] = r * (1 - r);
    tran[m * 3 + 2] = r * r;
  }

  // initialize site mutation probability matrix
  // diagonal is chance of no mutation
  // the diagonal "rotated by 90 degrees" is the chance of both positions
  // mutating
  // all other entries are chance of just one mutation
  pc[0][0] = pc[1][1] = pc[2][2] = pc[3][3] =
      (1 - mu) * (1 - mu); //	  probability of mutating no positions for each
  // parents haplotype
  pc[0][1] = pc[0][2] = pc[1][0] = pc[1][3] = pc[2][0] = pc[2][3] = pc[3][1] =
      pc[3][2] = mu * (1 - mu); //	  probability of mutating one position
  // for each parental haplotype
  pc[0][3] = pc[1][2] = pc[2][1] = pc[3][0] =
      mu * mu; //	  probability of mutating both positions for each
  // parental haplotype

  // initialize individual haplotypes
  for (unsigned i = 0; i < in; i++) {

    // define pointers to an individual's two haplotypes: ha and hb
    uint64_t *ha = &haps[i * 2 * wn], *hb = ha + wn;

    // here an individual's transition and probability matrices are
    // pulled out of the set of all individuals' matrices
    // t = genotype transition matrix? or genotype probability matrix?
    // e = phase emission matrix
    // p = a site's genotype probability (initially GL)
    // now iterate through each site
    fast *t = &temp[i * pn], *e = &emit[i * en], *p = &prob[i * 2];

    for (unsigned m = 0; m < mn; m++, t += 3, e += 4, p += hn) {

      // initialize genotype probabilities as genotype likelihoods
      if (is_x && male.find(name[i]) != male.end() &&
          !is_par[m]) { /// treat it differently for genders
        t[0] = max(1 - p[0] - p[1], 0.0f);
        t[1] = 0;
        t[2] = max(p[1], 0.0f);

        if (t[0] + t[2]) {
          t[0] /= t[0] + t[2];
          t[2] = 1 - t[0];
        } else
          t[0] = t[2] = 0.5;
      } else {

        // initial prob is the GL.
        t[0] = max(1 - p[0] - p[1], 0.0f);
        t[1] = max(p[0], 0.0f);
        t[2] = max(p[1], 0.0f);
      }

      // set each hap's bit randomly to 0 or 1 according to GLs, assuming flat
      // prior for GLs
      // leave both 0 if we have a ref hom (most common)
      if (gsl_rng_uniform(rng) < t[0])
        ;

      // choose one of the haps to be 1 and the other 0 if het
      // random phase
      else if (gsl_rng_uniform(rng) < t[1] / (t[1] + t[2])) {
        if (gsl_rng_uniform(rng) < 0.5)
          set1(ha, m);
        else
          set1(hb, m);
      }

      // hom
      else {
        set1(ha, m);
        set1(hb, m);
      }

      // initial emit assumes all states are by random mutation only.  basic
      // state is ref/ref
      for (unsigned j = 0; j < 4; j++)
        e[j] = pc[j][0] * t[0] + pc[j][1] * t[1] + pc[j][2] * t[1] +
               pc[j][3] * t[2];
    }
  }

  swap(temp, prob); // swap the assignments to each vector

  // create an unordered map version of site
  assert(m_sitesUnordered.size() == 0);

  for (auto oneSite : m_glSites)
    m_sitesUnordered.insert(std::make_pair(oneSite.pos, oneSite));

  assert(m_sitesUnordered.size() == m_glSites.size());

  // create unordered map version of names
  assert(m_namesUnordered.size() == 0);

  for (unsigned nameIdx = 0; nameIdx < name.size(); nameIdx++)
    m_namesUnordered.insert(std::make_pair(name[nameIdx], nameIdx));

  assert(m_namesUnordered.size() == name.size());

  // end of copy from SNPTools Impute
  // load ref haps
  if (!m_init.refPanelFiles.empty())
    LoadRefPanel(m_init.refPanelFiles, m_init.refPanelFilesType, Region());

  if (m_bUsingRefHaps) {

    // add ref haplotypes to sample haps
    haps.insert(haps.end(), m_vRefHaps.begin(), m_vRefHaps.end());

    // enlarge hnew so haps and hnew can be swapped
    // the ref haps are never updated, so they'll stick around forever
    if (!m_init.serializeHapUpdate)
      hnew.insert(hnew.end(), m_vRefHaps.begin(), m_vRefHaps.end());
  }

  // load the scaffold
  if (!m_init.scaffoldFiles.at("h").empty()) {
    LoadScaffold();

    // change phase and genotype of main haps to that from scaffold
    if (m_init.initPhaseFromScaffold)
      SetHapsAccordingToScaffold();

    // fix haplotypes from scaffold
    if (m_init.fixPhaseFromScaffold)
      FixEmitAccordingToScaffold();
  }
}

#ifndef NCUDA
// update the sampler to be used
// and replace the haplotypes when done

fast Insti::cudaSolve(HMMLike &hapSampler, unsigned sampleStride, fast pen) {

  assert(sampleStride == in);

  // sample four haps for N samples
  unsigned firstSampIdx{ 0 }, lastSampIdx{ 0 };
  vector<unsigned> propHaps =
      hapSampler.RunHMMOnSamples(firstSampIdx, lastSampIdx);

  // make sure the correct number of samples were updated
  assert(firstSampIdx == 0);
  assert(lastSampIdx == sampleStride - 1);

  // do a final round of haplotype estimation
  // this could be parallelized using threads perhaps? or an openmp pragma
  assert(m_init.serializeHapUpdate == false);
#pragma omp parallel for
  for (unsigned sampNum = firstSampIdx; sampNum <= lastSampIdx; ++sampNum)
    hmm_work(sampNum, &propHaps[sampNum * 4], pen);

  // replace old by new haplotypes
  swap(hnew, haps);

  // haplotypes in the hmmLike are updated automagically next time
  // RunHMMOnSamples is called

  // this should be the sum of likelihoods eventually when I feel like
  // implementing it...
  return 0;
}
#endif

// this part of the code seems to be responsible for:
// A - finding a set of four haps that are close to the current individual
// B - running the HMM and udating the individual I's haplotypes
// A takes much longer than B

/* CHANGES from impute.cpp:
   moved logging to solve from hmm_like()
   cycles is now stored in a private member variable and defined after
   load_bin() time
*/

// solve(individual, number of cycles, penalty, sampler)
fast Insti::solve(unsigned I, unsigned N, fast pen,
                  shared_ptr<Sampler> &sampler) {

  // write log header
  if (s_bIsLogging) {
    stringstream message;
    message << "##"
               "iteration\tindividual\tproposal\taccepted\tind1\tind2\tind3\t"
               "ind4" << endl;
    WriteToLog(message.str());
  }

  // pick 4 haplotype indices at random not from individual
  vector<unsigned> propHaps(4);

  for (unsigned j = 0; j < propHaps.size(); j++)
    propHaps[j] = sampler->SampleHap(I);

  // get a probability of the model for individual I given p
  fast curr = hmm_like(I, propHaps.data());

  // pick a random haplotype to replace with another one from all
  // haplotypes.  calculate the new probability of the model given
  // those haplotypes.
  // accept new set if probability has increased.
  // otherwise, accept with penalized probability
  MHSampler<unsigned> mhSampler(rng, curr, s_MHSamplerType, pen);
  for (unsigned n = 0; n < N; n++) { // sample haps for each sample N cycles
    unsigned rp = gsl_rng_get(rng) & 3, oh = propHaps[rp];

    // kickstart phasing and imputation by only sampling haps
    // from ref panel in first round
    if (n == 0 && s_bKickStartFromRef)
      propHaps[rp] = sampler->SampleHap(I, true);
    else
      propHaps[rp] = sampler->SampleHap(I, false);

    // get likelihood of proposed hap set
    fast prop = hmm_like(I, propHaps.data());

    // log all proposals and individuals proposed
    if (s_bIsLogging) {
      stringstream message;
      message << m_nIteration << "\t" << I << "\t" << prop;
      WriteToLog(message.str());
    }

    // do a MH acceptance of proposal
    mhSampler.SampleHap(propHaps[rp], oh, prop);

    if (s_bIsLogging) {
      stringstream message;
      message << "\t" << mhSampler.accepted();
      for (unsigned i = 0; i < 4; ++i)
        message << "\t" << propHaps[i];
      message << "\n";
      WriteToLog(message.str());
    }

    // Update relationship graph with proportion pen
    sampler->UpdatePropDistProp(propHaps, mhSampler.accepted(), I, pen);
  }

  // if we have passed the burnin cycles (n >= bn)
  // start sampling the haplotypes for output
  /*
  if (P) {
      uint16_t *pa = &pare[I * in];
      for (unsigned i = 0; i < 4; i++) pa[p[i] / 2]++;
  }
  */
  hmm_work(I, propHaps.data(), pen);
  return curr;
}

/* CHANGES from impute.cpp
   added a member variable for n so
*/

void Insti::estimate() {

  // choose a sampling scheme
  switch (m_estimator) {

  // simulated annealing MH
  case 0:

    // cluster!
    if (Insti::s_uNumClusters > 0) {
      switch (Insti::s_uClusterType) {

      // use kmedoids
      case 0:
      case 1:
        m_sampler = make_shared<KMedoids>(s_uClusterType, s_uNumClusters, haps,
                                          wn, mn, rng);
        break;

      // use kNN
      case 2:
        if (UsingScaffold())
          m_sampler = make_shared<KNN>(
              s_uNumClusters, m_scaffold.Haplotypes_uint64_t(),
              m_scaffold.NumWordsPerHap(), m_scaffold.NumSites(),
              m_init.scaffoldFreqLB, m_init.scaffoldFreqUB,
              m_init.scaffoldUsingMAF, rng, Insti::s_clusterDistanceMetric);
        else
          throw myException(
              "kNN sampler with no scaffold is not implemented yet");

        break;

      default:
        cout << "unexpected cluster type: " << s_uClusterType << endl;
        document();
      }
    }

    // just sample uniformly
    else
      m_sampler = make_shared<UnifSampler>(rng, in, hn + m_uNumRefHaps);

    break;

  case 1:
    estimate_EMC();
    return;

  case 2:

    // create a relationshipGraph object
    // initialize relationship matrix
    // create an in x uSamplingInds matrix
    m_sampler = make_shared<GraphSampler>(rng, in, hn + m_uNumRefHaps,
                                          GraphSampT::sampSamp);
    estimate_AMH();
    return;

  case 3:
    m_sampler = make_shared<GraphSampler>(rng, in, hn + m_uNumRefHaps,
                                          GraphSampT::sampHap);
    estimate_AMH();
    return;

  default:
    document();
  }

  // the default sampler shall be sampling randomly
  shared_ptr<Sampler> iterationSampler =
      make_shared<UnifSampler>(rng, in, hn + m_uNumRefHaps);

#ifndef NCUDA

  // create a cuda hapsampler object
  // repack gls for device
  const unsigned sampleStride = in;
  GLPackHelper::Init init(prob, *rng);
  init.numSamps = in;
  init.sampleStride = sampleStride;
  init.useVQ = true;
  init.numBitsPerGL = BITSPERCODE;
  GLPack glPack(init);

  // create hap sampler object
  HMMLike cudaHapSampler(haps, hn, glPack, m_uCycles, tran, &pc,
                         iterationSampler, *rng);
#endif

  timeval startTime, currentTime;
  gettimeofday(&startTime, NULL);

  cout.setf(ios::fixed);
  cout.precision(3);
  cout << m_tag << ":\titer\tpress\tlike\tfold\trunTime\texpectedRunTime"
       << endl;

  // log sampler state
  if (s_bIsLogging) {
    assert(!m_sLogFile.empty());
    try {
      m_sampler->Save(m_sLogFile, name);
    }
    catch (exception &e) {
      clog << m_tag << ": " << e.what() << "\n";
    }
  }

  // n is number of cycles = burnin + sampling cycles
  // increase penalty from 2/bn to 1 as we go through burnin
  // iterations.
  for (m_nIteration = 0; m_nIteration < bn + sn; ++m_nIteration) {
    fast sum = 0, iter = 0;
    fast pen = min<fast>((m_nIteration + 1.0f) / Insti::s_uSABurninGen, 1);
    pen *= pen; // pen = 1 after bn/2 iterations

    // recluster if it's the right generation
    // yes, don't recluster on first iteration
    bool reclustThisGen = false;
    if (m_init.reclusterEveryNGen > 0) {
      switch (m_init.reclusterStage[0]) {
      // recluster in sampling generation
      case 's':
        if (m_nIteration > bn &&
            (m_nIteration - bn) % m_init.reclusterEveryNGen == 0)
          reclustThisGen = true;
        break;
      // recluster in burnin generation
      case 'b':
        if (m_nIteration < bn && m_nIteration > 0 &&
            m_nIteration % m_init.reclusterEveryNGen == 0)
          reclustThisGen = true;
        break;
      // recluster in all gens
      case 'a':
        if (m_nIteration > 0 && m_nIteration % m_init.reclusterEveryNGen == 0)
          reclustThisGen = true;
        break;
      default:
        throw runtime_error("unexpected recluster stage option given: " +
                            string(m_init.reclusterStage));
      }
    }

    if (reclustThisGen) {
      if (std::dynamic_pointer_cast<KNN>(iterationSampler)) {
        cout << m_tag << ":\tReclustering haplotypes\n";
        m_sampler = make_shared<KNN>(
            s_uNumClusters, haps, WN, NUMSITES, m_init.scaffoldFreqLB,
            m_init.scaffoldFreqUB, m_init.scaffoldUsingMAF, rng,
            Insti::s_clusterDistanceMetric);
        cout << m_tag << ":\tReclustering complete.\n";
      } else
        throw runtime_error("Tried to recluster when not using KNN sampler");
    }

    // the uniform relationship "graph" will be used until
    // -M option says not to.
    if (m_nIteration < Insti::s_uStartClusterGen)
      iterationSampler = make_shared<UnifSampler>(rng, in, hn + m_uNumRefHaps);
    else
      iterationSampler = m_sampler;

// Phase each individual based on the rest of the individuals
// but don't update haps between updates??? !!! !!! !!!
#ifndef NCUDA
    // update the sampler, because it might have changed
    cudaHapSampler.UpdateSampler(iterationSampler);
    sum = cudaSolve(cudaHapSampler, sampleStride, pen);
#else
    if (m_init.serializeHapUpdate) {

      // need to update samples in random order
      std::vector<size_t> sampIndices;
      sampIndices.reserve(in);
      for (size_t i = 0; i < in; ++i)
        sampIndices.push_back(i);
      std::random_shuffle(sampIndices.begin(), sampIndices.end());

      // update samples one at a time
      for (unsigned i = 0; i < in; i++) {

        // re-update graph based on current haplotypes
        sum += solve(sampIndices.at(i), m_uCycles, pen, iterationSampler);
        iter += m_uCycles;
      }
    } else {
#pragma omp parallel for
      for (unsigned i = 0; i < in; i++) {

        // re-update graph based on current haplotypes
        sum += solve(i, m_uCycles, pen, iterationSampler);
        iter += m_uCycles;
      }
      swap(hnew, haps);
    }
#endif

    if (m_nIteration >= bn)
      for (unsigned i = 0; i < in; i++)
        replace(i); // call replace

    m_sampler->UpdatePropDistHaps(haps);

    // give an update
    gettimeofday(&currentTime, NULL);
    double runTimeInMin =
        static_cast<double>(currentTime.tv_sec - startTime.tv_sec) / 60;
    cout << m_tag << ":\t" << m_nIteration << '\t' << pen << '\t'
         << sum / in / mn << '\t' << iter / in / in << "\t" << runTimeInMin
         << "m\t" << runTimeInMin / (m_nIteration + 1) * (bn + sn) << "m\t"
         << endl;
  }

  cout << endl;
  result(); // call result
}

/* estimate_EMC -- Evolutionary Monte Carlo
   Here we try to increase the speed of convergence by
   running a parallel chain evolutionary monte carlo scheme.

   The reference for this implementation is
   "Advanced Markov Chain Monte Carlo Methods" by Liang, Liu and Carroll
   first edition?, 2010, pp. 128-132
*/

// solve(individual, number of cycles, penalty, burnin?)
fast Insti::solve_EMC(unsigned I, unsigned N, fast S) {
  DEBUG_MSG("Entering solve_EMC..." << endl);

  // for lack of a better place, define free parameters here
  fast fMutationRate = 0.3f; // see p.134 of Liang et al.
  fast fSelectTemp = 10000;
  fast fMaxTemp = Insti::s_uParallelChains;

  // write log header
  stringstream message;
  message << "##iteration\tindividual\tproposal\tchainID\tchainTemp"
          << "\tmutation" << endl;
  WriteToLog(message.str());

  // initialize emc chains with increasing temperatures
  vector<EMCChain> vcChains;
  vector<uint>
  vuChainTempHierarchy; // index of Chains sorted by temperature, ascending

  for (unsigned i = 0; i < Insti::s_uParallelChains; i++) {
    vcChains.push_back(EMCChain((i + 1) * fMaxTemp / Insti::s_uParallelChains,
                                fSelectTemp, I, in, i));
    DEBUG_MSG2("\tlength of vcChains\t" << vcChains.size() << endl);

    // initialize current likelihood
    // randomize parent haps
    for (unsigned j = 0; j < 4; j++) {
      do {
        unsigned uChosenHap = gsl_rng_get(rng) % vcChains[i].m_uHapNum;
        vcChains[i].setParent(j, uChosenHap);
        DEBUG_MSG3("\t\tchosen hap:\t" << uChosenHap << endl
                                       << "\t\tchosen parent:\t"
                                       << vcChains[i].getParent(j) << endl);
      } while (vcChains[i].getParent(j) / 2 == vcChains[i].m_uI);
    }

    DEBUG_MSG2("\tsetting likelihood" << endl;);

    // set likelihood
    vcChains[i].setLike(hmm_like(vcChains[i].m_uI, vcChains[i].getParents()));
    vuChainTempHierarchy.push_back(i);
  }

  // pick a random haplotype to replace with another one from all
  // haplotypes.  calculate the new probability of the model given
  // those haplotypes.
  // accept new set if probability has increased.
  // otherwise, accept with penalized probability
  for (unsigned n = 0; n < N; n++) { // fixed number of iterations
    DEBUG_MSG2("\tIteration " << n << endl);

    //  now choose whether to mutate or crossover
    bool bMutate = gsl_rng_uniform(rng) > fMutationRate;
    fast prop;

    if (bMutate) {

      // mutate
      DEBUG_MSG2("\tMutating...");

      // choose chain randomly (uniform)
      unsigned j = gsl_rng_get(rng) % Insti::s_uParallelChains;
      DEBUG_MSG2("\t"
                 << "Temp: " << vcChains[j].m_fTemp);
      fast curr = vcChains[j].getLike();

      // choose parent hap (rp) to mutate
      // replaced hap is stored in oh (Original Hap)
      unsigned rp = gsl_rng_get(rng) & 3, oh = vcChains[j].getParent(rp);

      // mutate parent hap
      do
        vcChains[j].setParent(rp, gsl_rng_get(rng) % hn);

      while (vcChains[j].getParent(rp) / 2 == I);

      // calculate acceptance probability
      prop = hmm_like(vcChains[j].m_uI, vcChains[j].getParents());

      if (prop > curr ||
          gsl_rng_uniform(rng) < exp((curr - prop) / vcChains[j].m_fTemp)) {
        vcChains[j].setLike(prop);
        WriteToLog(vcChains[j], bMutate);
      } else
        vcChains[j].setParent(rp, oh);
    } else {

      // crossover
      DEBUG_MSG2("\tCrossing Over...");

      // 1. choose random chain to work on by roulette wheel selection
      int iFirstChainIndex = RWSelection(vcChains);
      DEBUG_MSG2("\t\tFirst Chain:\t" << vcChains[iFirstChainIndex].m_uChainID
                                      << endl);
      DEBUG_MSG2("\t\tSelecting second chain:");

      // 2. select second chain at random (uniform) from remaining chains
      int iSecondChain;

      do {
        iSecondChain = gsl_rng_get(rng) % Insti::s_uParallelChains;
      } while (iSecondChain == iFirstChainIndex);

      //---- only crossover with certain probability depending
      //---- on crossover probabilities of parents
      // A. save parents and likelihoods of original chains
      EMCChain cFirstOrigChain = vcChains[iFirstChainIndex];
      EMCChain cSecondOrigChain = vcChains[iSecondChain];
      bool bFirstOrigChainHigherLike =
          cFirstOrigChain.getLike() > cSecondOrigChain.getLike();

      // B. cross over
      // uniform crossover: find haps to keep
      uint64_t cSelection = static_cast<uint64_t>(gsl_rng_get(rng) & 15);

      for (unsigned i = 0; i < 4; i++) {

        // if bit at location i is 1, exchange haps
        if (test(&cSelection, i)) {
          unsigned oh = vcChains[iFirstChainIndex].getParent(i);
          vcChains[iFirstChainIndex].setParent(
              i, vcChains[iSecondChain].getParent(i));
          vcChains[iSecondChain].setParent(i, oh);
        }
      }

      // update likelihoods of crossed over chains
      auto &rcFirstChain = vcChains[iFirstChainIndex];
      rcFirstChain.setLike(
          hmm_like(rcFirstChain.m_uI, rcFirstChain.getParents()));
      auto &rcSecondChain = vcChains[iSecondChain];
      rcSecondChain.setLike(
          hmm_like(rcSecondChain.m_uI, rcSecondChain.getParents()));
      bool const bFirstChainHigherLike =
          cFirstOrigChain.getLike() > cSecondOrigChain.getLike();

      // C. deterimen if cross-over was successful
      bool bCrossAccepted = false;

      // the order of likelihoods is not the same
      if (bFirstOrigChainHigherLike ^ bFirstChainHigherLike) {
        bCrossAccepted =
            gsl_rng_uniform(rng) <=
            exp(((cSecondOrigChain.getLike() - rcFirstChain.getLike()) /
                 cSecondOrigChain.m_fTemp) +
                (cFirstOrigChain.getLike() - rcSecondChain.getLike()) /
                    cFirstOrigChain.m_fTemp);
      }

      // the order of the likelihoods matches
      else {
        bCrossAccepted =
            gsl_rng_uniform(rng) <=
            exp(((cSecondOrigChain.getLike() - rcSecondChain.getLike()) /
                 cSecondOrigChain.m_fTemp) +
                (cFirstOrigChain.getLike() - rcFirstChain.getLike()) /
                    cFirstOrigChain.m_fTemp);
      }

      // replace old chains if cross-over was not accepted
      if (!bCrossAccepted) {
        vcChains[iFirstChainIndex] = cFirstOrigChain;
        vcChains[iSecondChain] = cSecondOrigChain;
        stringstream message;
        message << "# Unsuccessful Crossover\tChainIDs:\t"
                << rcFirstChain.m_uChainID << "\t" << rcSecondChain.m_uChainID
                << endl;
        WriteToLog(message.str());
      } else {

        // otherwise log changes to likelihood
        WriteToLog(vcChains[iFirstChainIndex], bMutate);
        WriteToLog(vcChains[iSecondChain], bMutate);
      }
    }

    // now try Insti::s_uParallelChains exchanges
    DEBUG_MSG2("\tExchanging..." << endl);
    unsigned uNumExchanges = 0;

    for (unsigned i = 0; i < Insti::s_uParallelChains; i++) {
      unsigned uFirstChainIndex = gsl_rng_get(rng) % Insti::s_uParallelChains;
      DEBUG_MSG3("\t\tfirstChainIndex " << uFirstChainIndex);
      unsigned uFirstChainHierarchyIndex =
          vuChainTempHierarchy[uFirstChainIndex];
      DEBUG_MSG3("\tfirst chain: " << vuChainTempHierarchy[uFirstChainIndex]);

      // selecting second chain
      unsigned uSecondChainIndex;

      if (uFirstChainIndex == 0)
        uSecondChainIndex = uFirstChainIndex + 1;
      else if (uFirstChainIndex == Insti::s_uParallelChains - 1)
        uSecondChainIndex = Insti::s_uParallelChains - 2;
      else if (gsl_rng_get(rng) & 1)
        uSecondChainIndex = uFirstChainIndex - 1;
      else
        uSecondChainIndex = uFirstChainIndex + 1;

      unsigned uSecondCHI = vuChainTempHierarchy[uSecondChainIndex];
      DEBUG_MSG3("\tsecond chain: " << vuChainTempHierarchy[uSecondChainIndex]);

      // MH step for exchange
      fast fAcceptProb =
          min<fast>(exp((vcChains[uFirstChainHierarchyIndex].getLike() -
                         vcChains[uSecondCHI].getLike()) *
                        ((1 / vcChains[uFirstChainHierarchyIndex].m_fTemp) -
                         (1 / vcChains[uSecondCHI].m_fTemp))),
                    1);
      DEBUG_MSG3("\taccept prob: " << fAcceptProb);

      // exchange with acceptance probability
      if (gsl_rng_uniform(rng) < fAcceptProb) {

        // exchange temperatures
        fast fTemp = vcChains[uFirstChainHierarchyIndex].m_fTemp;
        vcChains[uFirstChainHierarchyIndex].setTemp(
            vcChains[uSecondCHI].m_fTemp);
        vcChains[uSecondCHI].setTemp(fTemp);

        // exchange location in vcChains
        std::swap(vuChainTempHierarchy[uFirstChainIndex],
                  vuChainTempHierarchy[uSecondChainIndex]);
        ++uNumExchanges;
      }

      DEBUG_MSG3("\tnumExchanges: " << uNumExchanges << endl);
    }

    // keep track of number of exchanges
    stringstream message;
    message << "# Number of Exchanges out of total:\t" << uNumExchanges << "\t"
            << Insti::s_uParallelChains << endl;
    WriteToLog(message.str());
  }

  // now select a chain for sampling according to roulette wheel selection
  int iFirstChainIndex = RWSelection(vcChains);
  auto &rcFirstChain = vcChains[iFirstChainIndex];

  // if we have passed the burnin cycles (n >= bn)
  // start sampling the haplotypes for output
  /*
  if (P) {
      uint16_t *pa = &pare[I * in];
      for (unsigned i = 0; i < 4; i++) pa[rcFirstChain.getParent(i) / 2]++;
  }
  */
  DEBUG_MSG("Updating individual " << I << "\n");

  // update haplotypes of I
  //    cerr << "parents (2 lines):" << endl;
  //    cerr << rcFirstChain.m_auParents[1] << endl;
  //    cerr << vcChains[ iFirstChainIndex ].m_auParents[1] << endl;
  hmm_work(I, rcFirstChain.getParents(), S);
  return rcFirstChain.getLike();
}

/* estimate_EMC -- Evolutionary Monte Carlo
   Here we try to increase the speed of convergence by
   running a parallel chain evolutionary monte carlo scheme.

   The reference for this implementation is
   "Advanced Markov Choin Monte Carlo Methods" by Liang, Liu and Carroll
   first edition?, 2010, pp. 128-132
*/

void Insti::estimate_EMC() {
  cout.setf(ios::fixed);
  cout.precision(3);
  cout << "Running Evolutionary Monte Carlo\n";
  cout << "iter\tpress\tlike\tfold\n";

  // n is number of cycles = burnin + sampling cycles
  // increase penalty from 2/bn to 1 as we go through burnin
  // iterations.
  for (unsigned n = 0; n < bn + sn; n++) {
    m_nIteration = n;
    fast sum = 0, pen = min<fast>(2 * (n + 1.0f) / bn, 1), iter = 0;
    pen *= pen; // pen = 1 after bn/2 iterations

    for (unsigned i = 0; i < in; i++) {
      sum += solve_EMC(i, m_uCycles,
                       pen); // call solve=> inputs the sample number,
      iter += m_uCycles;
    }

    if (!m_init.serializeHapUpdate)
      swap(hnew, haps);

    if (n >= bn)
      for (unsigned i = 0; i < in; i++)
        replace(i); // call replace

    cout << n << '\t' << pen << '\t' << sum / in / mn << '\t' << iter / in / in
         << '\n';
  }

  cout << endl;
  result(); // call result
}

/* estimate_AMH -- Adaptive Metropolis Hastings
   Here we try to increase the speed of convergence by
   keeping track of which individuals tend to copy from each other

   The reference for parts of this implementation is
   "Advanced Markov Choin Monte Carlo Methods" by Liang, Liu and Carroll
   first edition?, 2010, pp. 309
*/
void Insti::estimate_AMH() {
  cout.setf(ios::fixed);
  cout.precision(3);
  cout << "Running Adaptive Metropolis Hastings\n";
  cout << "iter\tpress\tlike\tfold" << endl;

  // n is number of cycles = burnin + sampling cycles
  // increase penalty from 2/bn to 1 as we go through burnin
  // iterations.
  for (unsigned n = 0; n < bn + sn; n++) {

    //        cout << "iter\t" << n << endl;
    m_nIteration = n;
    fast sum = 0, pen = min<fast>(2 * (n + 1.0f) / bn, 1), iter = 0;
    pen *= pen; // pen = 1 after bn/2 iterations

    // update all individuals once
    for (unsigned i = 0; i < in; i++) {

      //            if( i % 1000 == 0)
      //                cout << "cycle\t" << i << endl;
      sum += solve(i, m_uCycles, pen, m_sampler);
      iter += m_uCycles;
    }

    if (!m_init.serializeHapUpdate)
      swap(hnew, haps);

    if (n >= bn)
      for (unsigned i = 0; i < in; i++)
        replace(i); // call replace

    cout << n << '\t' << pen << '\t' << sum / in / mn << '\t' << iter / in / in
         << endl;
  }

  cout << endl;
  result(); // call result
}

void Insti::save_relationship_graph(string sOutputFile) {
  vector<string> vsSampNames;
  vsSampNames.insert(vsSampNames.end(), name.begin(), name.end());

  for (unsigned i = 0; i < ceil(m_uNumRefHaps / 2); i++)
    vsSampNames.push_back(string("refSamp") + sutils::uint2str(i));

  m_sampler->Save(sOutputFile, vsSampNames);
}

void Insti::save_vcf(const char *F, string commandLine) {
  string temp = F;
  temp += ".vcf.gz";
  ofile vcfFD(temp);
  vcfFD << std::setprecision(3);
  vcfFD << "##fileformat=VCFv4.0\n";
  vcfFD << "##source=WTCHG:INSTIv" << VERSION_MAJOR << "." << VERSION_MINOR
        << "." << VERSION_XSTR(VERSION_REVISION) << "\n";

  vcfFD << "##phaseAndImputeCall=" << commandLine << "\n";
  vcfFD << "##iteration=%u\n";
  vcfFD << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
  vcfFD << "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Phred-scaled "
        << "genotype posterior probability\">\n";
  vcfFD << "##FORMAT=<ID=APP,Number=2,Type=Float,Description=\"Phred-scaled "
        << "allelic posterior probability, P(Allele=1|Haplotype)\">\n";
  vcfFD << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

  for (unsigned i = 0; i < in; i++)
    vcfFD << "\t" << name[i];

  for (unsigned m = 0; m < mn; m++) {
    vcfFD << "\n" << m_glSites[m].chr << "\t" << m_glSites[m].pos << "\t.\t"
          << m_glSites[m].ref << "\t" << m_glSites[m].alt
          << "\t100\tPASS\t.\tGT:GP:APP";

    fast *p = &prob[m * hn];

    for (unsigned i = 0; i < in; i++, p += 2) {
      fast prr = (1 - p[0]) * (1 - p[1]),
           pra = (1 - p[0]) * p[1] + p[0] * (1 - p[1]), paa = p[0] * p[1];

      if (prr >= pra && prr >= paa)
        vcfFD << "\t0|0";
      else if (pra >= prr && pra >= paa) {
        if (p[0] > p[1])
          vcfFD << "\t1|0";
        else
          vcfFD << "\t0|1";
      } else
        vcfFD << "\t1|1";

      // test for any p being zero
      vector<fast> vfProb = { prr, pra, paa, p[0], p[1] };

      for (auto &phred : vfProb) {
        if (phred == 0)
          phred = FLT_MAX_10_EXP * 10;
        else if (phred == 1)
          phred = 0;
        else if (phred > 1 || phred < 0) {
          cerr << "error: probability " << phred << " is not between 0 and 1";
          exit(1);
        } else
          phred = -10 * log10(phred);
      }

      // print a sample's GP and APP fields
      for (unsigned i = 0; i < vfProb.size(); ++i) {
        if (i == 0 || i == 3)
          vcfFD << ":";
        vcfFD << vfProb[i];
        if (i != 2 && i != 4)
          vcfFD << ",";
      }
    }
  }

  vcfFD << "\n";
  vcfFD.close();
}

void Insti::SetHapsAccordingToScaffold() {

  if (m_scaffold.NumHaps() == 0)
    return;

  SetHapsAccordingToHapPanel(m_scaffold);
}

void Insti::SetHapsAccordingToHapPanel(const HapPanel &hapPanel) {

  assert(!hapPanel.empty());
  assert(hapPanel.NumHaps() == hn);

  // pull out each hap and change it to match information contained in the
  // hapPanel
  vector<uint64_t> hapPanelHaps = hapPanel.Haplotypes_uint64_t();
  size_t nwph = hapPanel.NumWordsPerHap();
  assert(2 * in == hapPanel.NumHaps());
  for (unsigned hapNum = 0; hapNum < 2 * in; hapNum++) {

    // define pointers to an individual's two haplotype hap
    uint64_t *hap = &haps[hapNum * wn];
    uint64_t *hapPanelHap = &hapPanelHaps[hapNum * nwph];

    // check each site if it needs to be replaced
    unsigned hapPanelSiteIdx = 0;

    unsigned siteIdx = 0;
    for (unsigned hapPanelPositionIdx = 0;
         hapPanelPositionIdx < hapPanel.NumSites(); hapPanelPositionIdx++) {

      // find the site in site that matches the hapPanel position
      for (; siteIdx < m_glSites.size(); siteIdx++) {
        if (hapPanel.snp_is_equal(m_glSites[siteIdx], hapPanelPositionIdx))
          break;
      }

      // if this is not true, then the hapPanel position could not be found in
      // site
      assert(siteIdx < m_glSites.size());

      // set the site to what is in the hapPanel
      if (test(hapPanelHap, hapPanelSiteIdx))
        set1(hap, siteIdx);
      else
        set0(hap, siteIdx);

      ++hapPanelSiteIdx;
    }
  }
}

void Insti::LoadScaffold() {

  assert(!m_init.scaffoldFiles.at("h").empty());

  // redefine region if necessary
  Bio::Region scaffoldRegion = m_runRegion;
  if (m_init.scaffoldExtraRegionSize > 0) {
    unsigned overhang = m_init.scaffoldExtraRegionSize;
    unsigned start = m_glSites[0].pos;
    start = start < overhang ? 0 : start - overhang;
    unsigned end = m_glSites.back().pos + overhang;
    scaffoldRegion = Bio::Region(m_glSites[0].chr, start, end);
  }

  if (!m_init.scaffoldFiles.at("s").empty())
    LoadHapsSamp(m_init.scaffoldFiles.at("h"), m_init.scaffoldFiles.at("s"),
                 InstiPanelType::SCAFFOLD, scaffoldRegion);

  // then this might be a VCF instead!
  else if (m_init.scaffoldFiles.at("h").size() >= 7 &&
           m_init.scaffoldFiles.at("h").compare(
               m_init.scaffoldFiles.at("h").size() - 7, 7, ".vcf.gz") == 0)
    LoadVCFGZ(m_init.scaffoldFiles.at("h"), InstiPanelType::SCAFFOLD,
              scaffoldRegion);
  else {
    throw std::runtime_error("[insti] Error while loading scaffold haps/sample "
                             "files: Need to define a sample file");
  }
}

void Insti::FixEmitAccordingToScaffold() {

  clog << m_tag << ": [insti] Fixing emission matrix according to scaffold"
       << "\n";
  const unordered_map<string, string> &files = m_init.scaffoldFiles;

  // figure out how to what type of scaffold we are loading
  if (files.at("h").empty())
    throw runtime_error("Fixing of phase according to scaffold was requested, "
                        "but no haplotypes were supplied (-h)");

  // prepare loading of haplotypes
  HapPanelHelper::Init init;
  init.keepSites = m_glSites;
  init.keepSampleIDs = name;
  init.allowReorderOfSitesAtPos = true;
  vector<string> inFiles;

  // assume vcfgz
  if (files.at("s").empty()) {
    clog << m_tag << ": [insti] Loading haplotypes in GZVCF for fixing of phase"
         << "\n";
    inFiles.push_back(files.at("h"));
    init.hapFileType = HapPanelHelper::HapFileType::VCF;
  }
  // assume WTCCC
  else if (files.at("l").empty()) {
    clog << m_tag
         << ": [insti] Loading haplotypes in WTCCC format fixing of phase"
         << "\n";
    inFiles.push_back(files.at("h"));
    inFiles.push_back(files.at("s"));
    if (inFiles[0].size() > 11 &&
        inFiles[0].compare(inFiles[0].size() - 11, 11, ".tabhaps.gz") == 0)
      init.hapFileType = HPH::HapFileType::WTCCC_TBX;
    else
      init.hapFileType = HPH::HapFileType::WTCCC;

  }

  // assume IMPUTE2
  else {
    clog << m_tag
         << ": [insti] Loading haplotypes in IMPUTE2 format fixing of phase"
         << "\n";
    inFiles.push_back(files.at("h"));
    inFiles.push_back(files.at("s"));
    inFiles.push_back(files.at("l"));
    init.hapFileType = HPH::HapFileType::IMPUTE2;
  }

  // set include variants
  if (!m_init.fixPhaseAlwaysKeepStrandFile.empty())
    init.alwaysKeepStrandFile = m_init.fixPhaseAlwaysKeepStrandFile;

  // load scaffold and filter sites
  HapPanel haplotypes(move(inFiles), move(init));
  haplotypes.FilterSitesOnAlleleFrequency(
      m_init.fixPhaseFrequencyLowerBound, m_init.fixPhaseFrequencyUpperBound,
      m_init.fixPhaseReferenceAlleleFrequency);

  // fix emit according to scaffold
  // and build indicator vector of which sites are fixed
  assert(m_sitesAtWhichPhaseIsFixed.empty());
  m_sitesAtWhichPhaseIsFixed.resize(wn);
  for (unsigned i = 0; i < in; i++) {
    fast *e = &emit[i * en];

    for (unsigned m = 0; m < mn; ++m, e += 4) {
      int scaffoldSiteIdx = haplotypes.FindSiteIndex(m_glSites[m]);
      if (scaffoldSiteIdx > -1) {
        const size_t phase =
            (haplotypes.GetAllele(static_cast<size_t>(scaffoldSiteIdx), i * 2)
             << 1) |
            (haplotypes.GetAllele(static_cast<size_t>(scaffoldSiteIdx),
                                  i * 2 + 1));
        assert(phase < 4);
        for (unsigned j = 0; j < 4; j++)
          e[j] = pc[j][phase];

        // mark this site as fixed
        set1(m_sitesAtWhichPhaseIsFixed.data(), m);
      }
    }
  }

  // and now initialize haps at sites that are fixed to pre-existing haplotype
  // phase
  SetHapsAccordingToHapPanel(haplotypes);
}

void Insti::document() {
  cerr << "Author\tWarren W. Kretzschmar @ Marchini Group @ Universiy of "
          "Oxford - Statistics";
  cerr << "\n\nThis code is based on SNPTools impute:";
  cerr << "\nhaplotype imputation by cFDSL distribution";
  cerr << "\nAuthor\tYi Wang @ Fuli Yu' Group @ BCM-HGSC";
  cerr << "\n\nusage\tinsti [options] GLFile";
  cerr << "\n\t-n <fold>       sample size*fold of nested MH sampler "
          "iteration "
          "(2)";
  cerr << "\n\t-o <name>\tPrefix to use for output files. \".phased\" will be "
          "added to to the prefix.";
  cerr << "\n\t-P <thread>     number of threads (0=MAX,default=1)";
  cerr << "\n\t-x <gender>     impute x chromosome data";
  cerr << "\n\t-e <file>       write log to file";
  cerr << "\n\t-g <file>       genetic map (required)";
  cerr << "\n\t-F <string>     Input GL file type (bin). Valid values are "
          "\"bin\" or \"bcf\"";
  cerr << "\n\t-R <string>     Region to run imputation on.  Only works with "
          "VCF/BCF input GL file. Region should be of the form 20:100-1000 "
          "for "
          "chromosome 20 starting at genomic position 100 and ending at "
          "genomic position 1000 (inclusive).";
  cerr << "\n\t-S <file>       list of samples to keep from GLs";

  cerr << "\n\n    GENERATION OPTIONS";
  cerr << "\n\t-m <mcmc>       sampling generations (200)";
  cerr << "\n\t-C <integer>    number of cycles to estimate an individual's "
          "parents before updating";
  cerr << "\n\t-B <integer>    number of simulated annealing generations (28)";
  cerr << "\n\t-i <integer>    number of non-simulated annealing burnin "
          "generations (28)";
  cerr << "\n\t-M <integer>    generation number at which to start "
          "clustering, "
          "0-based (28)";

  cerr << "\n\n    HAPLOTYPE ESTIMATION OPTIONS";
  cerr << "\n\t-E <integer>    choice of estimation algorithm (0)";
  cerr << "\n\t                0 - Metropolis Hastings with simulated "
          "annealing";
  cerr << "\n\t                1 - Evolutionary Monte Carlo with -p parallel "
          "chains";
  cerr << "\n\t                2 - Adaptive Metropolis Hastings - "
          "sample/sample matrix";
  cerr << "\n\t                3 - Adaptive Metropolis Hastings - "
          "sample/haplotype matrix";
  cerr << "\n\t-p <integer>    number of parallel chains to use in parallel "
          "estimation algorithms";
  cerr << "\n\t                (at least 2, default 5)";
  cerr << "\n\t-K <integer>    number of clusters to use for haplotypes "
          "clustering (0 = option is off).";
  cerr << "\n\t                Does not currently work with -k option";
  cerr << "\n\t-t <integer>    Cluster type (0)";
  cerr << "\n\t                0 - k-Medoids -- PAM";
  cerr << "\n\t                1 - k-Medoids -- Park and Jun 2008";
  cerr << "\n\t                2 - k-Nearest Neighbors -- IMPUTE2 (-K is the "
          "number of haplotypes to keep)";
  cerr << "\n\t-T              Use shared tract length as distance metric for "
          "clustering";
  cerr << "\n\t-r <integer>    Recluster every -r generations. Only works "
          "when "
          "-t=2.  Will start at generation -M";
  cerr << "\n\t-u <char>       Only recluster in burnin ('b'), sampling ('s') "
          "or all ('a') generations.  Respects -M and -r options. Default is "
          "'s'.";

  cerr << "\n\n    REFERENCE PANEL OPTIONS";
  cerr << "\n\t-H <file>       Comma separated haplotype files "
          "(VCF,HAP/SAMP,HAP/LEG/SAMP)";
  cerr << "\n\t-L <file>       haplotype format of -H (one of "
          "VCF,WTCCC,IMPUTE2)";
  cerr << "\n\t-k              Kickstart phasing by using only ref panel in "
          "first iteration";
  cerr << "\n\n    SCAFFOLD OPTIONS";
  cerr << "\n\t-h <file>       WTCCC style HAPS file";
  cerr << "\n\t-s <file>       WTCCC style SAMPLE file";
  cerr << "\n\t-q <float>      Lower bound of variant allele frequency "
          "([0-1], "
          "default 0.05"
       << ") "
          "above which sites are used for clustering from scaffold.";
  cerr << "\n\t-Q <float>      Upper bound of variant allele frequency "
          "([0-1], "
          "default 0.95"
       << ") "
          "below which sites are used for clustering from scaffold.";
  cerr << "\n\t-a              Use minor allele frequency instead of variant "
          "allele frequency for clustering and applying -q and -Q.";
  cerr << "\n\t-f              Initialize phase according to scaffold (default "
          "off).";
  cerr << "\n\t-O <integer>    Size in genomic coordinates of the regions "
          "past "
          "the regions specified by the GLs to include in the scaffold "
          "(default 0).";

  cerr << "\n\n    FIXING PHASE FROM SCAFFOLD OPTIONS";
  cerr << "\n\tSpecifying any of the options below turns on fixing of phase "
          "from scaffold.";
  cerr << "\n\tThe scaffold is specified with the -h, -s and -l options. ";
  cerr << "\n\t-w <float>      Set lower bound of allele frequency of sites to "
          "fix ([0-1], default 0).";
  cerr << "\n\t-W <float>      Set upper bound of allele frequency of sites to "
          "fix ([0-1], default 1).";
  cerr << "\n\t-A              Use minor allele frequency instead of variant "
          "allele frequency for applying -w and -W";
  cerr << "\n\t-I <file>       Strand file of sites to always keep";
  //  cerr << "\n\t-I <file>       Site include list (format: "
  //          "'CHROM\tPOS\tREF\tALT')";
  cerr << "\n\n";
  exit(1);
}

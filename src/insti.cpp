#include "insti.hpp"

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;
#ifndef NCUDA
using namespace HMMLikeCUDA;
#endif
using namespace Bio;

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
string Insti::s_sRefLegendFile = "";
string Insti::s_sRefHapFile = "";
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

Insti::Insti(InstiHelper::Init &init)
    : m_reclusterEveryNGen(init.reclusterEveryNGen),
      m_scaffoldHapsFile(init.scaffoldHapsFile),
      m_scaffoldSampleFile(init.scaffoldSampleFile),
      m_initPhaseFromScaffold(init.initPhaseFromScaffold),
      m_geneticMap(init.geneticMap), m_init(init),
      m_tag(boost::uuids::random_generator()()) {

  omp_set_num_threads(init.numThreads);
  // bounds check estimator input
  if (init.estimator < 4)
    m_estimator = init.estimator;
  else
    throw std::runtime_error("[insti] Estimator needs to be less than 4");
}

// return the probability of the model given the input haplotypes P and
// emission and transition matrices of individual I
// call Impute::hmm_like and print out result
fast Insti::hmm_like(unsigned I, uint *P) { return Impute::hmm_like(I, P); }

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

void Insti::load_bin(const string &binFile) {
  bool bRetVal = Impute::load_bin(binFile.c_str());

  if (bRetVal == false)
    throw std::runtime_error("[insti] Unable to load input .bin file: " +
                             string(binFile));

  // setting range object based on load bin
  if (site.empty())
    throw std::runtime_error(
        "[insti] Loaded data seems to contain no sites in file: " +
        string(binFile));
  m_runRegion = Region(site[0].chr, site[0].pos, site.back().pos);

  // setting number of cycles to use
  // here is best place to do it because in is defined in load_bin()
  if (s_uCycles > 0)
    m_uCycles = s_uCycles;
  else {
    m_uCycles = nn * in; // this was how snptools does it
  }
}

// HAPS/SAMPLE sample file
void Insti::OpenSample(const string &sampleFile,
                       vector<string> &fillSampleIDs) {

  cout << m_tag << ": [insti] Loading samples file: " << sampleFile << endl;

  // read in sample file
  ifile sampleFD(sampleFile);

  if (!sampleFD.isGood())
    throw myException("Could not open file: " + sampleFile);

  string buffer;

  // discard header
  unsigned lineNum = 0;
  unsigned numCols = 0;

  while (getline(sampleFD, buffer, '\n')) {
    if (buffer.size() == 0)
      throw myException("Error in sample file " + sampleFile +
                        ": empty lines detected.");

    lineNum++;
    vector<string> tokens;
    sutils::tokenize(buffer, tokens);

    // make sure header start is correct
    if (lineNum == 1) {
      vector<string> headerTokenized;
      string header = "ID_1 ID_2 missing";
      sutils::tokenize(header, headerTokenized);

      for (unsigned i = 0; i != headerTokenized.size(); i++) {
        if (tokens[i] != headerTokenized[i])
          throw myException(
              "Error in sample file (" + sampleFile +
              "): header start does not match:\n\t" + header +
              "\n.  Instead the first line of the header is:\n\t" + buffer +
              "\n");
      }

      numCols = tokens.size();
      continue;
    }

    if (tokens.size() != numCols)
      throw myException(
          "Error in sample file (" + sampleFile + ") at line number " +
          sutils::uint2str(lineNum) + ": number of columns (" +
          sutils::uint2str(tokens.size()) + ") does not match previous line");

    // ignore second line of samples file
    if (lineNum == 2)
      continue;

    // now add sample id to vector
    fillSampleIDs.push_back(tokens[0]);
  }
}

void Insti::OpenTabHaps(const string &hapsFile, vector<vector<char> > &loadHaps,
                        vector<snp> &loadSites) {
  cout << m_tag << ": [insti] Loading tabixed haps file: " << hapsFile << endl;

  // clear all the containers that are going to be filled up
  vector<snp> fillSites;
  vector<vector<char> > fillHaps;
  fillSites.reserve(site.size());
  fillHaps.reserve(site.size());

  assert(m_sitesUnordered.size() == site.size());

  string file(hapsFile);
  Tabix tabix(file);
  string region(m_runRegion.AsString());
  tabix.setRegion(region);

  // start parsing haplotype lines
  string buffer;
  unsigned numHaps = 0;
  unsigned previousPos = 0;
  unsigned matchSites = 0;
  unsigned glSiteIdx =
      0; // used for checking sites against already loaded glSites
  while (tabix.getNextLine(buffer) && glSiteIdx != site.size()) {

    // split on tab
    vector<string> tokens;
    boost::split(tokens, buffer, boost::is_any_of("\t"));
    if (tokens[0] != site[glSiteIdx].chr)
      throw std::runtime_error("[insti] Found site on chromosome '" +
                               tokens[0] + "' when chromosome '" + site[0].chr +
                               "' was expected");
    if (tokens.size() < 6)
      throw std::runtime_error("[insti] Site has no haplotypes at site: " +
                               buffer.substr(0, 50));

    // parse position
    size_t lastChar;
    unsigned pos = stoul(tokens[2], &lastChar);
    if (lastChar != tokens[2].size())
      throw std::runtime_error(
          "[insti] Malformed site: could not parse position: " +
          buffer.substr(0, 50));

    // check to see if this site exists in our gl panel
    auto canSiteR = m_sitesUnordered.equal_range(pos);

    // the load site is not in the GL sites
    if (distance(canSiteR.first, canSiteR.second) == 0)
      continue;

    Bio::snp newSite(tokens[0], pos, tokens[3], tokens[4]);
    // test each of the sites to see if it matches
    bool match = false;
    for (auto candidateSite = canSiteR.first; candidateSite != canSiteR.second;
         ++candidateSite) {
      const snp &cSite = candidateSite->second;
      if ((newSite.pos == cSite.pos) && (newSite.chr == cSite.chr) &&
          (newSite.ref == cSite.ref) && newSite.alt == cSite.alt) {
        match = true;
        break;
      }
    }
    if (!match)
      continue;

    if (pos == previousPos)
      throw std::runtime_error("[insti] Multiple variants at the same "
                               "position are not implemented yet");
    previousPos = pos;
    ++matchSites;

    // check to make sure the correct number of haplotypes exist in line
    if (numHaps == 0)
      numHaps = tokens.size() - 5;
    else if (numHaps != tokens.size() - 5)
      throw std::runtime_error("[insti] Wrong number of alleles at site: " +
                               buffer.substr(0, 50));

    // stuff haplotypes into vector
    vector<char> alleles;
    alleles.reserve(numHaps);
    for (unsigned allNum = 0; allNum < numHaps; ++allNum)
      if (tokens[allNum + 5] == "0")
        alleles.push_back(0);
      else if (tokens[allNum + 5] == "1")
        alleles.push_back(1);
      else
        throw std::runtime_error(
            "[insti] Malformed site: allele is neither 0 or 1 at site: " +
            buffer.substr(0, 50));

    // stuff that vector into vector of char vectors...
    fillHaps.push_back(move(alleles));

    // keep the site because it matches
    fillSites.push_back(move(newSite));
  }

  // we are done do a swap to fill the incoming vectors
  swap(loadSites, fillSites);
  swap(loadHaps, fillHaps);
}

// read in the haps file
// store haps and sites
void Insti::OpenHaps(const string &hapsFile, vector<vector<char> > &loadHaps,
                     vector<snp> &sites) {

  cout << m_tag << ": [insti] Loading haps file: " << hapsFile << endl;
  ifile hapsFD(hapsFile);

  if (!hapsFD.isGood())
    throw myException("Could not open file: " + hapsFile);

  string buffer;
  unsigned lineNum = 0;
  unsigned keptSites = 0;
  unsigned numHaps = 0;
  sites.clear();
  loadHaps.clear();
  assert(m_sitesUnordered.size() == site.size());

  // for making sure input file is sorted by position
  unsigned previousPos = 0;
  // create a map of site positions
  while (getline(hapsFD, buffer, '\n')) {
    if (keptSites == m_sitesUnordered.size())
      break;

    //    if (lineNum % 1000 == 0)
    //      cout << "Sites kept:\t" << keptSites << " / " << lineNum << "\n";

    lineNum++;

    //// only read position
    // figure out where the first three spaces are located
    vector<size_t> spaceIdxs;
    size_t nextStartIdx = 0;
    for (unsigned i = 0; i < 3; ++i) {
      if (i == 0)
        spaceIdxs.push_back(buffer.find_first_of(" "));
      else
        spaceIdxs.push_back(buffer.find_first_of(" ", nextStartIdx));
      if (spaceIdxs[i] == string::npos)
        throw myException("Space in Haps file not found where expected");
      nextStartIdx = spaceIdxs[i] + 1;
    }

    // load chrom and check input haps is correct chrom
    string chr = buffer.substr(0, spaceIdxs[0]);
    if (chr != site[0].chr)
      throw myException("Found site on chromosome '" + chr +
                        "' when chromosome '" + site[0].chr + "' was expected");

    // move ahead and extract position from third field
    size_t endReadIdx = 0;
    unsigned pos = stoul(
        buffer.substr(spaceIdxs[1] + 1, spaceIdxs[2] - (spaceIdxs[1] + 1)),
        &endReadIdx, 0);

    // make sure the whole field was parsed!
    if (endReadIdx != spaceIdxs[2] - (spaceIdxs[1] + 1))
      throw myException("Input haplotypes line field three was read as " +
                        sutils::uint2str(pos) +
                        " and does not seem to be an unsigned integer.\n" +
                        "Read was " + sutils::uint2str(endReadIdx) +
                        " char(s) long" + "\nThe field was actually " +
                        sutils::uint2str(spaceIdxs[2]) + " char(s) long");

    // make sure input sites are sorted by position
    if (pos < previousPos)
      throw myException("Input haplotypes file " + hapsFile +
                        " needs to be sorted by position");
    previousPos = pos;

    // start loading only once we hit the first site
    if (pos < site[0].pos)
      continue;

    // stop loading sites if the current site is past the last site position in
    // the GLs
    if (pos > site.back().pos)
      break;

    // only keep sites that we know of
    auto foundSite = m_sitesUnordered.find(pos);

    if (foundSite == m_sitesUnordered.end())
      continue;

    // split entire line for processing
    vector<string> tokens;
    boost::split(tokens, buffer, boost::is_any_of(" "));

    // make sure header start is correct
    if (numHaps == 0) {
      if (tokens.size() <= 5)
        throw myException("haps file " + hapsFile +
                          " contains too few columns (" +
                          sutils::uint2str(tokens.size()) + ")");

      // count number of haps
      numHaps = tokens.size() - 5;

    } else {
      if (tokens.size() - 5 != numHaps)
        throw myException(
            "Every row of haplotypes file must have the same number of "
            "columns");
    }

    // store site of haplotype
    sites.push_back(snp(chr, pos, tokens[3], tokens[4]));
    vector<char> loadSite;
    loadSite.reserve(numHaps);

    for (unsigned i = 5; i != tokens.size(); i++) {
      int val = atoi(tokens[i].c_str());

      if (val == 0)
        loadSite.push_back(0);
      else if (val == 1)
        loadSite.push_back(1);
      else {
        throw myException("All alleles are not 0 or 1 In haplotypes file " +
                          hapsFile + " at line number " +
                          sutils::uint2str(lineNum) + " (0-based)");
      }
    }

    loadHaps.push_back(loadSite);

    // keeping this site
    keptSites++;
  }

  cout << "Sites kept:\t" << keptSites << " / " << lineNum << "\n";

  if (numHaps == 0) {
    cerr << "Number of haplotypes in haps file is 0.  Haps file empty?\n";
    exit(1);
  }

  assert(loadHaps[0].size() == numHaps);
}

void Insti::OpenVCFGZ(const string &vcf, const string &region,
                      vector<vector<char> > &loadHaps, vector<snp> &loadSites,
                      vector<string> &ids) {
  loadHaps.clear();
  loadSites.clear();
  ids.clear();

  // open vcf using tabixpp
  string inVCF = vcf;
  string inRegion = region;
  Tabix tbx(inVCF);
  if (tbx.setRegion(inRegion) != true)
    throw std::runtime_error("Malformed region string: " + region);

  // get the header
  string header;
  tbx.getHeader(header);

  using qi::phrase_parse;
  using qi::char_;
  using qi::omit;
  using qi::lit;
  //  using qi::int_;
  //  using qi::_1;
  using qi::as_string;
  using phoenix::push_back;

  // parse header for sample names
  auto first = header.begin();
  auto last = header.end();
  qi::parse(first, last,

            //  Begin grammar
            omit[*("##" > +(qi::graph | ' ' | '\t') > lit('\n'))] >
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" >
                +(qi::graph | ' ') % omit[lit('\t')] > omit[lit("\n")],
            //  End grammar
            ids);

  if (first != last)
    throw std::runtime_error(
        "Could not find #CHROM line in VCF header of file: " + vcf +
        "\nHeader found:\n" + header);

  // start churning through the region of interest
  unsigned siteNum = 0;
  string line;
  GTParser::vcf_grammar<string::const_iterator> grammar(siteNum);
  while (tbx.getNextLine(line)) {
    auto cfirst = line.cbegin();
    vector<char> loadSite;
    loadSite.reserve(ids.size() * 2);

    // load values from string into these objects
    //    string chr = "", ref = "", alt = "";
    unsigned genomic_pos;
    vector<string> firstCols(8);
    vector<t_genotype> genotypes;
    qi::parse(cfirst, line.cend(), grammar, firstCols[0], genomic_pos,
              firstCols[1], firstCols[2], firstCols[3], firstCols[4],
              firstCols[5], firstCols[6], firstCols[7], genotypes);

    /*
        phrase_parse(
            first, line.end(),

            //  Begin grammar
            (as_string[+(char_ - lit("\t"))][ref(chr) = _1] > int_[ref(pos) =
       _1] >
             omit[+(char_ - lit("\t"))] >
             as_string[+(char_ - lit("\t"))][ref(ref) = _1] >
             as_string[+(char_ - lit("\t"))][ref(alt) = _1] >
             omit[repeat(4)[+(char_ - lit("\t"))]] >
             *(char_('0', '1')[push_back(phoenix::ref(site), atoi(_1))] >
       lit('|') >
               char_('0', '1')[push_back(phoenix::ref(site), atoi(_1))]) >
             lit("\n")),
            //  End grammar
            lit("\t"));
    */

    // make sure parser parsed everything
    if (cfirst != line.cend())
      throw std::runtime_error(
          "Malformed VCF line at line " + to_string(siteNum + 1) +
          " of VCF body in file: " + vcf + "\nLine: " + line);
    /*    if (first != line.end())
          throw std::runtime_error("VCF line was not completely parsed at line "
       +
                                   to_string(siteNum + 1) +
                                   " of VCF body in file: " + vcf);
    */
    // make sure the correct number of haps were found
    if (genotypes.size() != ids.size())
      throw std::runtime_error(
          "Malformed VCF line: incorrect number of columns at line " +
          to_string(siteNum + 1) + " of VCF body in file: " + vcf);

    // copy alleles to site
    for (auto geno : genotypes) {
      loadSite.push_back(geno.allele1);
      loadSite.push_back(geno.allele2);
    }

    loadSites.push_back(
        snp(firstCols[0], genomic_pos, firstCols[2], firstCols[3]));
    loadHaps.push_back(move(loadSite));
    ++siteNum;
  }
}

void Insti::LoadVCFGZ(const string &vcf, InstiPanelType panel_t,
                      const string &region) {

  cout << "[insti] Loading haplotypes from VCF: " << vcf << endl;
  if (vcf.size() < 7 || vcf.compare(vcf.size() - 7, 7, ".vcf.gz") != 0)
    throw std::runtime_error(
        "[insti] input VCF file does not end in '.vcf.gz': " + vcf);

  // load haps file
  CheckPanelPrereqs(panel_t);

  vector<vector<char> > loadHaps;
  vector<snp> loadSites;
  vector<string> scaffoldSampleIDs;

  // read the haps and sites from a haps file
  OpenVCFGZ(vcf, region, loadHaps, loadSites, scaffoldSampleIDs);
  if (loadHaps.empty())
    throw myException("VCF contains no haplotypes: " + vcf);

  // get sample information if we are using a scaffold
  // make sure samples match
  if (panel_t == InstiPanelType::SCAFFOLD) {
    SubsetSamples(scaffoldSampleIDs, loadHaps);
    OrderSamples(scaffoldSampleIDs, loadHaps);

    assert(!loadHaps.empty());
    MatchSamples(scaffoldSampleIDs, loadHaps[0].size());
  }

  assert(!loadHaps.empty());
  vector<vector<char> > filtHaps;
  vector<snp> filtSites;

  FilterSites(loadHaps, loadSites, filtHaps, filtSites, panel_t);

  // loading haplotypes into place
  LoadHaps(filtHaps, filtSites, scaffoldSampleIDs, panel_t);
};

void Insti::LoadHapsSamp(const string &hapsFile, const string &sampleFile,
                         InstiPanelType panelType) {

  // make sure both files are defined
  if (sampleFile.empty() && panelType == InstiPanelType::SCAFFOLD)
    throw std::runtime_error("[insti] Error while loading scaffold haps/sample "
                             "files: Need to define a sample file");

  if (hapsFile.size() == 0)
    throw std::runtime_error("[insti] Need to define a haps file");

  // load haps file
  CheckPanelPrereqs(panelType);

  vector<vector<char> > loadHaps;
  vector<snp> loadSites;

  // read the haps and sites from a haps file
  if (hapsFile.size() > 11 &&
      hapsFile.compare(hapsFile.size() - 11, 11, ".tabhaps.gz") == 0)
    OpenTabHaps(hapsFile, loadHaps, loadSites);
  else
    OpenHaps(hapsFile, loadHaps, loadSites);
  if (loadHaps.empty())
    throw myException("Haplotypes file is empty: " + hapsFile);

  // get sample information if we are using a scaffold
  // make sure samples match
  vector<string> scaffoldSampleIDs;
  if (panelType == InstiPanelType::SCAFFOLD) {
    OpenSample(sampleFile, scaffoldSampleIDs);
    SubsetSamples(scaffoldSampleIDs, loadHaps);
    OrderSamples(scaffoldSampleIDs, loadHaps);

    assert(!loadHaps.empty());
    MatchSamples(scaffoldSampleIDs, loadHaps[0].size());
  }

  assert(!loadHaps.empty());
  vector<vector<char> > filtHaps;
  vector<snp> filtSites;

  FilterSites(loadHaps, loadSites, filtHaps, filtSites, panelType);

  // loading haplotypes into place
  LoadHaps(filtHaps, filtSites, scaffoldSampleIDs, panelType);
}

void Insti::OrderSamples(vector<string> &loadIDs,
                         vector<vector<char> > &loadHaps) {

  assert(loadIDs.size() == m_namesUnordered.size());

  if (!loadHaps.empty())
    assert(loadIDs.size() * 2 == loadHaps[0].size());

  vector<unsigned> orderedLoadNameIdxs;
  orderedLoadNameIdxs.reserve(loadIDs.size());

  for (unsigned idIdx = 0; idIdx < loadIDs.size(); idIdx++) {
    auto foundID = m_namesUnordered.find(loadIDs[idIdx]);

    if (foundID != m_namesUnordered.end()) {
      assert(foundID->second < m_namesUnordered.size());
      orderedLoadNameIdxs.push_back(foundID->second);
    } else
      assert(false); // programming error
  }

  // sort names according to orderedLoadNameIdxs
  vector<string> tempIDs(orderedLoadNameIdxs.size());

  for (unsigned nameOrderIdx = 0; nameOrderIdx < orderedLoadNameIdxs.size();
       nameOrderIdx++)
    tempIDs[orderedLoadNameIdxs[nameOrderIdx]] = loadIDs[nameOrderIdx];

  std::swap(tempIDs, loadIDs);

  // sort haplotypes according to orderedLoadNameIdxs
  for (unsigned siteIdx = 0; siteIdx < loadHaps.size(); siteIdx++) {

    // subset haps
    vector<char> tempSite(orderedLoadNameIdxs.size() * 2);

    for (unsigned nameOrderIdx = 0; nameOrderIdx < orderedLoadNameIdxs.size();
         nameOrderIdx++) {
      tempSite[orderedLoadNameIdxs[nameOrderIdx] * 2] =
          loadHaps[siteIdx][nameOrderIdx * 2];
      tempSite[orderedLoadNameIdxs[nameOrderIdx] * 2 + 1] =
          loadHaps[siteIdx][nameOrderIdx * 2 + 1];
    }

    assert(tempSite.size() == orderedLoadNameIdxs.size() * 2);
    std::swap(tempSite, loadHaps[siteIdx]);
    assert(loadHaps[siteIdx].size() == orderedLoadNameIdxs.size() * 2);
  }

  if (!loadHaps.empty())
    assert(loadIDs.size() * 2 == loadHaps[0].size());
}
void Insti::SubsetSamples(vector<string> &loadIDs,
                          vector<vector<char> > &loadHaps) {

  assert(loadIDs.size() >= m_namesUnordered.size());

  if (!loadHaps.empty())
    assert(loadIDs.size() * 2 == loadHaps[0].size());

  vector<unsigned> idIdxsToKeep;

  for (unsigned idIdx = 0; idIdx < loadIDs.size(); idIdx++) {
    auto foundID = m_namesUnordered.find(loadIDs[idIdx]);

    if (foundID != m_namesUnordered.end())
      idIdxsToKeep.push_back(idIdx);
  }

  if (idIdxsToKeep.size() != name.size())
    throw myException(
        "Samples file does not contain right number of matching IDs to GL IDs "
        "(" +
        sutils::uint2str(idIdxsToKeep.size()) + "/" +
        sutils::uint2str(name.size()) + ")");

  // subset haplotypes
  vector<string> tempIDs;
  tempIDs.reserve(idIdxsToKeep.size());

  for (auto keepIdx : idIdxsToKeep)
    tempIDs.push_back(loadIDs[keepIdx]);

  std::swap(tempIDs, loadIDs);

  for (unsigned siteIdx = 0; siteIdx < loadHaps.size(); siteIdx++) {

    // subset haps
    vector<char> tempSite;
    tempSite.reserve(idIdxsToKeep.size() * 2);

    for (auto keepIdx : idIdxsToKeep) {
      tempSite.push_back(loadHaps[siteIdx][keepIdx * 2]);
      tempSite.push_back(loadHaps[siteIdx][keepIdx * 2 + 1]);
    }

    assert(tempSite.size() == idIdxsToKeep.size() * 2);
    std::swap(tempSite, loadHaps[siteIdx]);
    assert(loadHaps[siteIdx].size() == idIdxsToKeep.size() * 2);
  }

  if (!loadHaps.empty())
    assert(loadIDs.size() * 2 == loadHaps[0].size());
}

// only keep sites in main gl set
void Insti::FilterSites(vector<vector<char> > &loadHaps, vector<snp> &loadSites,
                        vector<vector<char> > &filtHaps, vector<snp> &filtSites,
                        InstiPanelType panelType) {

  assert(loadSites.size() > 0);
  assert(loadSites.size() == loadHaps.size());
  assert(loadHaps[0].size() > 0);
  filtHaps.clear();
  filtSites.clear();

  // initializing filtHaps
  //    for(auto & iter  : filtHaps)
  //      iter.resize(loadSites.size());

  const unsigned numHaplotypesLoaded = loadHaps[0].size();

  if (panelType == InstiPanelType::REFERENCE) {

    // the reference panel may not have less sites than the GLs
    assert(site.size() <= loadSites.size());
    filtHaps.reserve(site.size());
  } else if (panelType == InstiPanelType::SCAFFOLD) {

    // the scaffold may have less or more sites than the GLs
  } else
    assert(false); // programming error

  // now keep all sites that are in the GL set of sites
  // look up each load site to see if its in the GL sites
  unsigned previousPos = 0;
  unsigned matchSites = 0;
  for (unsigned loadSiteIdx = 0; loadSiteIdx < loadSites.size();
       ++loadSiteIdx) {
    const snp &loadSite = loadSites[loadSiteIdx];
    auto canSiteR = m_sitesUnordered.equal_range(loadSite.pos);

    // the load site is not in the GL sites
    if (distance(canSiteR.first, canSiteR.second) == 0)
      continue;

    // test each of the sites to see if it matches
    bool match = false;
    for (auto candidateSite = canSiteR.first; candidateSite != canSiteR.second;
         ++candidateSite) {
      const snp &cSite = candidateSite->second;
      if ((loadSite.pos == cSite.pos) &&
          (loadSite.chr.empty() ? true : loadSite.chr == cSite.chr) &&
          (loadSite.ref == cSite.ref) && loadSite.alt == cSite.alt) {
        match = true;
        break;
      }
    }
    if (match) {
      if (loadSite.pos < previousPos)
        throw std::runtime_error(
            "[insti] Haplotypes file is not ordered by position: ");
      if (loadSite.pos == previousPos)
        throw std::runtime_error("[insti] Multiple variants at the same "
                                 "position are not implemented yet");
      previousPos = loadSite.pos;
      ++matchSites;

      // keep the site because it matches
      filtSites.push_back(loadSite);
      filtHaps.push_back(move(loadHaps[loadSiteIdx]));
    }
  }

  cout << "[insti] Number of loaded haplotypes: " << numHaplotypesLoaded
       << endl;
  cout << "[insti] Number of haplotypes left after filtering: "
       << filtHaps[0].size() << endl;
  cout << "[insti] Number of candidate sites: " << site.size() << endl;
  cout << "[insti] Number of sites kept: " << filtSites.size() << endl;

  if (filtHaps[0].size() != numHaplotypesLoaded)
    throw std::runtime_error(
        "Error: No sites in loaded panel matched already loaded sites. Is your "
        "input file sorted?");

  // we always want to find at most as many matches as are in GL sites
  if (matchSites > site.size())
    throw std::runtime_error("[insti] too many matching sites found");

  assert(matchSites > 0);

  if (panelType == InstiPanelType::REFERENCE)
    if (site.size() != matchSites)
      throw std::runtime_error(
          "[insti] Could not find all GL sites in loaded haplotype panel");
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

// make sure samples loaded and samples in main GLs match
// would only be used for scaffold at the moment
void Insti::MatchSamples(const vector<std::string> &IDs, unsigned numHaps) {

  // make sure ids match names pulled from bin files
  if (name.size() != IDs.size())
    throw myException("number of IDs loaded (" + sutils::uint2str(IDs.size()) +
                      ") does not match number of IDs in GLs (" +
                      sutils::uint2str(name.size()) + ").");

  for (unsigned i = 0; i < name.size(); i++) {
    if (name[i] != IDs[i])
      throw myException(
          "ID number " + sutils::uint2str(i) +
          " (0-based): bin and IDs from sample file do not match (" + name[i] +
          " and " + IDs[i] + " respectively) do not agree.");
  }

  // make sure number of haplotypes in haplotypes matches number of samples *2
  if (IDs.size() * 2 != numHaps)
    throw myException("Number of haplotypes according to haps file (" +
                      sutils::uint2str(numHaps) + ") and sample file (" +
                      sutils::uint2str(IDs.size() * 2) + ") do not match.");
}

// put the haplotypes in the right place in the program structure
void Insti::LoadHaps(vector<vector<char> > &inHaps, vector<snp> &inSites,
                     vector<string> &inSampleIDs, InstiPanelType panelType) {

  assert(inHaps.size() == inSites.size());
  // convert char based haplotypes to bitvector
  unsigned numHaps = inHaps[0].size();

  // store the haplotypes in the correct place based on what type of panel we
  // are loading
  switch (panelType) {
  case InstiPanelType::REFERENCE: {
    HapPanel temp;
    vector<uint64_t> storeHaps =
        temp.Char2BitVec(inHaps, GetNumWords(), WORDMOD + 1);
    storeHaps.swap(m_vRefHaps);
    m_uNumRefHaps = numHaps;
    cout << "Reference panel haplotypes\t" << m_uNumRefHaps << endl;
    return;
  }

  case InstiPanelType::SCAFFOLD: {
    static_assert(WORDMOD >= 0, "WORDMOD global is not greater or equal to 0");
    m_scaffold.Init(inHaps, inSites, inSampleIDs);
    cout << "Scaffold haplotypes\t" << m_scaffold.NumHaps() << endl;

    try {
      if (m_scaffold.NumHaps() != hn)
        throw myException(
            "Error while reading scaffold: Scaffold needs to have two "
            "haplotypes for every input sample");
    }
    catch (exception &e) {
      cout << e.what() << endl;
      exit(1);
    };

    return;
  }

  default:
    assert(false); // this should not happen
  }
}

void Insti::CheckPanelPrereqs(InstiPanelType panelType) {
  switch (panelType) {
  case InstiPanelType::REFERENCE:
    m_bUsingRefHaps = true;
    assert(m_vRefHaps.size() == 0);
    break;

  case InstiPanelType::SCAFFOLD:
    assert(!m_scaffold.Initialized());
    break;

  default:
    assert(false); // this should not happen
  }
}

vector<snp> Insti::OpenLegend(string legendFile) {

  // read in legend file
  ifile legendFD(legendFile);

  if (!legendFD.isGood())
    throw myException("Could not open file: " + legendFile);

  string buffer;
  vector<snp> loadLeg;

  // discard header
  unsigned lineNum = 0;

  while (getline(legendFD, buffer, '\n')) {
    if (buffer.size() == 0)
      throw myException("Error in legend file " + legendFile +
                        ": empty lines detected.");

    lineNum++;
    vector<string> tokens;
    sutils::tokenize(buffer, tokens);

    // make sure header start is correct
    if (lineNum == 1) {
      vector<string> headerTokenized;
      string header = "id position a0 a1";
      sutils::tokenize(header, headerTokenized);

      for (unsigned i = 0; i != headerTokenized.size(); i++) {
        if (tokens[i] != headerTokenized[i])
          throw myException(
              "Error in legend file " + legendFile +
              ": header start does not match:\n\t" + header +
              "\n.  Instead the first line of the header is:\n\t" + buffer +
              "\n");
      }

      continue;
    }

    // add each site to loadLeg
    loadLeg.push_back(
        snp("", strtoul(tokens[1].c_str(), NULL, 0), tokens[2], tokens[3]));
  }

  return loadLeg;
}

vector<vector<char> > Insti::OpenHap(string hapFile) {

  // read in the hap file
  ifile hapFD(hapFile);

  if (!hapFD.isGood())
    throw myException("Could not open file: " + hapFile);

  string buffer;
  int lineNum = -1;
  unsigned uNumHaps = 0;
  vector<vector<char> > loadHaps;

  while (getline(hapFD, buffer, '\n')) {
    lineNum++;

    // read line
    vector<string> tokens;
    sutils::tokenize(buffer, tokens);

    // make sure header start is correct
    if (lineNum == 0) {

      // count number of haps
      uNumHaps = tokens.size();
    }

    // make sure the number of haps does not change
    if (tokens.size() != uNumHaps)
      throw myException(
          "Every row of hap file must have the same number of columns");

    // store haplotypes
    vector<char> inSite;
    inSite.reserve(uNumHaps);

    for (unsigned i = 0; i != tokens.size(); i++) {
      int val = atoi(tokens[i].c_str());

      if (val == 0)
        inSite.push_back(0);
      else if (val == 1)
        inSite.push_back(1);
      else {
        throw myException("All alleles are not 0 or 1 In haplotypes file " +
                          hapFile + " at line number " +
                          sutils::uint2str(lineNum) + " (0-based)");
      }
    }

    loadHaps.push_back(inSite);
  }

  if (uNumHaps == 0)
    throw myException("num haps is 0.  Haps file empty?\n");

  return loadHaps;
}

bool Insti::LoadHapLegSamp(const string &legendFile, const string &hapFile,
                           const string &sampleFile, InstiPanelType panelType) {

  // make sure required files are defined
  if (legendFile.size() == 0) {
    cerr << "Need to define a legend file if defining a hap file\n";
    document();
  }

  if (hapFile.size() == 0) {
    cerr << "Need to define a hap file if defining a legend file\n";
    document();
  }

  if (sampleFile.size() == 0 && panelType == InstiPanelType::SCAFFOLD) {
    cerr << "Need to define a sample file if using scaffold\n";
    document();
  }

  // for now sample loading is not implemented
  if (sampleFile.size() > 0) {
    cerr << "for now sample loading is not implemented in the sampleghap "
            "paradigm";
    document();
  }
  CheckPanelPrereqs(panelType);

  // Load Samples (not implemented)
  vector<string> sampleIDs;

  // Load the site list in the legend
  cout << "Loading legend file: " << legendFile << endl;

  try {
    auto legend = OpenLegend(legendFile);

    /*catch (exception& e) {
        cout << e.what() << " in legend file " << legendFile << endl;
        exit(1);
        }*/

    cout << "Loading hap file: " << hapFile << endl;
    auto loadHaps = OpenHap(hapFile);

    vector<vector<char> > filtHaps;
    vector<snp> filtSites;
    FilterSites(loadHaps, legend, filtHaps, filtSites, panelType);

    LoadHaps(filtHaps, filtSites, sampleIDs, panelType);
  }
  catch (exception &e) {
    cerr << "Error loading haplotypes file " << hapFile << ": " << e.what()
         << endl;
    exit(2);
  }

  return true;
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
  hn = in * 2;             // number of haps
  haps.resize(hn * wn);    // space to store all haplotypes
  hnew.resize(hn * wn);    // number of haplotypes = 2 * number of samples  ...
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
  fast mu = 0; //, rho;

  for (unsigned i = 1; i < hn; i++)
    mu += 1.0 / i;

  // mu * (mn-1) looks like the Watterson estimator of the population mutation
  // rate theta
  // 0.5 / (posi[mn - 1] - posi[0]) / density looks like a correction for finite
  // sites
  mu = 1 / mu;
  //  rho = 0.5 * mu * (mn - 1) / (posi[mn - 1] - posi[0]) / density;
  mu = mu / (hn + mu); // rho is recombination rate?  mu is mutation rate

  // initialzie the site transition matrix tran
  // posi is recombination between its position and previous
  // r < 1; the larger the number of haplotypes, the smaller r gets
  // tran is a site's recombination probability matrix
  // r therefore must be a recombination rate estimate
  for (unsigned m = mn - 1; m; m--) {

    // replaced ad-hoc genetic distance estimate by genetic map
    //    posi[m] = (posi[m] - posi[m - 1]) * rho;
    //    fast r = posi[m] / (posi[m] + hn);
    float r = m_geneticMap.GeneticDistance(posi[m - 1], posi[m]);

    // threshold genetic map to max value of 1
    r = r > 1 ? 1 : r;

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

  for (auto oneSite : site)
    m_sitesUnordered.insert(std::make_pair(
        oneSite.pos, snp(oneSite.chr, oneSite.pos, oneSite.all.substr(0, 1),
                         oneSite.all.substr(1, 1))));

  assert(m_sitesUnordered.size() == site.size());

  // create unordered map version of names
  assert(m_namesUnordered.size() == 0);

  for (unsigned nameIdx = 0; nameIdx < name.size(); nameIdx++)
    m_namesUnordered.insert(std::make_pair(name[nameIdx], nameIdx));

  assert(m_namesUnordered.size() == name.size());

  // end of copy from SNPTools Impute
  // load ref haps
  if (s_sRefLegendFile.size() > 0 || s_sRefHapFile.size() > 0)
    LoadHapLegSamp(s_sRefLegendFile, s_sRefHapFile, "",
                   InstiPanelType::REFERENCE);

  if (m_bUsingRefHaps) {

    // add ref haplotypes to sample haps
    haps.insert(haps.end(), m_vRefHaps.begin(), m_vRefHaps.end());

    // enlarge hnew so haps and hnew can be swapped
    // the ref haps are never updated, so they'll stick around forever
    hnew.insert(hnew.end(), m_vRefHaps.begin(), m_vRefHaps.end());
  }

  // load the scaffold
  if (!m_scaffoldHapsFile.empty()) {
    if (!m_scaffoldSampleFile.empty())
      LoadHapsSamp(m_scaffoldHapsFile, m_scaffoldSampleFile,
                   InstiPanelType::SCAFFOLD);

    // then this might be a VCF instead!
    else if (m_scaffoldHapsFile.size() >= 7 &&
             m_scaffoldHapsFile.compare(m_scaffoldHapsFile.size() - 7, 7,
                                        ".vcf.gz") == 0)
      LoadVCFGZ(m_scaffoldHapsFile, InstiPanelType::SCAFFOLD,
                m_runRegion.AsString());
    else {
      throw std::runtime_error(
          "[insti] Error while loading scaffold haps/sample "
          "files: Need to define a sample file");
    }
  }
  // change phase and genotype of main haps to that from scaffold
  if (m_initPhaseFromScaffold) {
    SetHapsAccordingToScaffold();
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
              s_uNumClusters, (*m_scaffold.Haplotypes()),
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
      cerr << e.what() << endl;
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
    // yes, don't recluster on iteration bn
    if (m_reclusterEveryNGen > 0 && m_nIteration > bn)
      if ((m_nIteration - bn) % m_reclusterEveryNGen == 0)
        if (std::dynamic_pointer_cast<KNN>(iterationSampler)) {
          cout << m_tag << ":\tReclustering haplotypes\n";
          m_sampler = make_shared<KNN>(
              s_uNumClusters, haps, WN, NUMSITES, m_init.scaffoldFreqLB,
              m_init.scaffoldFreqUB, m_init.scaffoldUsingMAF, rng,
              Insti::s_clusterDistanceMetric);
          cout << m_tag << ":\tReclustering complete.\n";
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
    for (unsigned i = 0; i < in; i++) {

      // re-update graph based on current haplotypes
      sum += solve(i, m_uCycles, pen, iterationSampler);
      iter += m_uCycles;
    }

    swap(hnew, haps);
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
    vcfFD << "\n" << site[m].chr << "\t" << site[m].pos << "\t.\t"
          << site[m].all[0] << "\t" << site[m].all[1]
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

  assert(m_scaffold.Initialized());

  if (m_scaffold.NumHaps() == 0)
    return;

  assert(m_scaffold.NumHaps() == hn);

  // pull out each hap and change it to match information contained in the
  // scaffold
  for (unsigned hapNum = 0; hapNum < 2 * in; hapNum++) {

    // define pointers to an individual's two haplotype hap
    uint64_t *hap = &haps[hapNum * wn];
    uint64_t *scaffHap = m_scaffold.Hap(hapNum);

    // check each site if it needs to be replaced
    unsigned scaffoldSiteIdx = 0;

    unsigned siteIdx = 0;
    for (unsigned scaffoldPositionIdx = 0;
         scaffoldPositionIdx < m_scaffold.NumSites(); scaffoldPositionIdx++) {

      // find the site in site that matches the scaffold position
      for (; siteIdx < site.size(); siteIdx++) {
        if (site[siteIdx].pos == m_scaffold.Position(scaffoldPositionIdx))
          break;
      }

      // if this is not true, then the scaffold position could not be found in
      // site
      assert(siteIdx < site.size());

      // set the site to what is in the scaffold
      if (test(scaffHap, scaffoldSiteIdx))
        set1(hap, siteIdx);
      else
        set0(hap, siteIdx);

      ++scaffoldSiteIdx;
    }
  }
}

void Insti::document() {
  cerr << "Author\tWarren W. Kretzschmar @ Marchini Group @ Universiy of "
          "Oxford - Statistics";
  cerr << "\n\nThis code is based on SNPTools impute:";
  cerr << "\nhaplotype imputation by cFDSL distribution";
  cerr << "\nAuthor\tYi Wang @ Fuli Yu' Group @ BCM-HGSC";
  cerr << "\n\nusage\timpute [options] 1.bin 2.bin ...";
  //  cerr << "\n\t-d <density>    relative SNP density to Sanger sequencing
  // (1)";

  //    cerr << "\n\t-b <burn>       burn-in generations (56)";
  cerr << "\n\t-l <file>       list of input files";
  cerr << "\n\t-n <fold>       sample size*fold of nested MH sampler "
          "iteration "
          "(2)";
  cerr << "\n\t-o <name>\tPrefix to use for output files";
  cerr << "\n\t-P <thread>     number of threads (0=MAX,default=1)";
  cerr << "\n\t-v <vcf>        integrate known genotype in VCF format";
  cerr << "\n\t-c <conf>       confidence of known genotype (0.9998)";
  cerr << "\n\t-x <gender>     impute x chromosome data";
  cerr << "\n\t-e <file>       write log to file";
  cerr << "\n\t-g <file>       genetic map (required)";

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
  cerr << "\n\t-r <integer>    Recluster every -r generations. Only works when "
          "-t=2";

  cerr << "\n\n    REFERENCE PANEL OPTIONS";
  cerr << "\n\t-H <file>       IMPUTE2 style HAP file";
  cerr << "\n\t-L <file>       IMPUTE2 style LEGEND file";
  cerr << "\n\t-k              Kickstart phasing by using only ref panel in "
          "first iteration";
  cerr << "\n\n    SCAFFOLD OPTIONS";
  cerr << "\n\t-h <file>       WTCCC style HAPS file";
  cerr << "\n\t-s <file>       WTCCC style SAMPLE file";
  cerr << "\n\t-q <float>      Lower bound of variant allele frequency ([0-1], "
          "default 0.05) "
          "above which sites are used for clustering from scaffold.";
  cerr << "\n\t-Q <float>      Upper bound of variant allele frequency ([0-1], "
          "default 0.05) "
          "below which sites are used for clustering from scaffold.";
  cerr << "\n\t-a              Use minor allele frequency instead of variant "
          "allele frequency for clustering and applying -q and -Q."
          "below which sites are used for clustering from scaffold.";
  cerr << "\n\t-f              Fix phase according to scaffold (default off).";
  cerr << "\n\n";
  exit(1);
}

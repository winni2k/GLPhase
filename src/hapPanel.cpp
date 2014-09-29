#include "hapPanel.hpp"

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;
using namespace Bio;

namespace HapPanelHelper {

namespace qi = boost::spirit::qi;

// converts a vector of vectors of chars (0,1) to a bitvector
vector<uint64_t> Char2BitVec(const vector<vector<char> > &inHaps,
                             unsigned numWords, unsigned wordSize) {

  assert(ceil(static_cast<double>(inHaps.size()) / wordSize) == numWords);

  unsigned numSites = inHaps.size();
  unsigned numHaps = inHaps[0].size();
  vector<uint64_t> storeHaps;
  storeHaps.resize(numWords * numHaps);

  for (unsigned siteNum = 0; siteNum < numSites; siteNum++) {
    for (unsigned hapNum = 0; hapNum < numHaps; hapNum++) {
      if (inHaps[siteNum][hapNum] == 0)
        HapPanelHelper::set0(&storeHaps[hapNum * numWords], siteNum);
      else if (inHaps[siteNum][hapNum] == 1)
        HapPanelHelper::set1(&storeHaps[hapNum * numWords], siteNum);
      else {
        cout << "programming error";
        throw 1;
      }
    }
  }
  return storeHaps;
}

void set1(uint64_t *P, unsigned I) {

  // I >> Uint64_TShift is moving along the array according to which uint64_t
  // I is in
  // e.g. I <= 63 is first uint64_t, and so forth in blocks of 64 bits
  P[I >> WORDSHIFT] |= static_cast<uint64_t>(1) << (I & WORDMOD);
}

void set0(uint64_t *P, unsigned I) {
  P[I >> WORDSHIFT] &= ~(static_cast<uint64_t>(1) << (I & WORDMOD));
}

// make sure samples loaded and samples in main GLs match
// would only be used for scaffold at the moment
void MatchSamples(const std::vector<std::string> &canonicalIDs,
                  const std::vector<std::string> &IDs, unsigned numHaps) {

  // make sure ids match names pulled from bin files
  if (canonicalIDs.size() != IDs.size())
    throw runtime_error("number of IDs loaded (" +
                        sutils::uint2str(IDs.size()) +
                        ") does not match number of IDs expected (" +
                        sutils::uint2str(canonicalIDs.size()) + ").");

  for (unsigned i = 0; i < canonicalIDs.size(); i++) {
    if (canonicalIDs[i] != IDs[i])
      throw runtime_error(
          "ID number " + sutils::uint2str(i) +
          " (0-based): expected IDs and IDs from sample file do not match (" +
          canonicalIDs[i] + " and " + IDs[i] + " respectively) do not agree.");
  }

  // make sure number of haplotypes in haplotypes matches number of samples *2
  if (IDs.size() * 2 != numHaps)
    throw runtime_error("Number of haplotypes according to haps file (" +
                        sutils::uint2str(numHaps) + ") and sample file (" +
                        sutils::uint2str(IDs.size() * 2) + ") do not match.");
}

// HAPS/SAMPLE sample file
vector<string> OpenSample(const string &sampleFile) {

  vector<string> fillSampleIDs;
  cout << "[HapPanel] Loading samples file: " << sampleFile << endl;

  // read in sample file
  ifile sampleFD(sampleFile);

  if (!sampleFD.isGood())
    throw runtime_error("Could not open sample file: [" + sampleFile + "]");

  string buffer;

  // discard header
  unsigned lineNum = 0;
  unsigned numCols = 0;

  while (getline(sampleFD, buffer, '\n')) {
    if (buffer.size() == 0)
      throw runtime_error("Error in sample file " + sampleFile +
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
          throw runtime_error(
              "Error in sample file (" + sampleFile +
              "): header start does not match:\n\t" + header +
              "\n.  Instead the first line of the header is:\n\t" + buffer +
              "\n");
      }

      numCols = tokens.size();
      continue;
    }

    if (tokens.size() != numCols)
      throw runtime_error(
          "Error in sample file (" + sampleFile + ") at line number " +
          sutils::uint2str(lineNum) + ": number of columns (" +
          sutils::uint2str(tokens.size()) + ") does not match previous line");

    // ignore second line of samples file
    if (lineNum == 2)
      continue;

    // now add sample id to vector
    fillSampleIDs.push_back(tokens[0]);
  }
  clog << "[HapPanel] Loaded " << fillSampleIDs.size() << " sample IDs" << endl;

  return fillSampleIDs;
}

// read in the haps file
// store haps and sites
void OpenHaps(const string &hapsFile, const Region &region,
              vector<vector<char> > &loadHaps, vector<snp> &loadSites) {

  cout << "[HapPanel] Loading haps file: " << hapsFile << endl;
  ifile hapsFD(hapsFile);

  if (!hapsFD.isGood())
    throw runtime_error("Could not open haps file [" + hapsFile + "]");

  string buffer;
  unsigned lineNum = 0;
  unsigned keptSites = 0;
  unsigned numHaps = 0;
  loadSites.clear();
  loadHaps.clear();

  // for making sure input file is sorted by position
  unsigned previousPos = 0;
  // create a map of site positions
  while (getline(hapsFD, buffer, '\n')) {

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
        throw runtime_error("Space in Haps file not found where expected");
      nextStartIdx = spaceIdxs[i] + 1;
    }

    // load chrom and check input haps is correct chrom
    string chr = buffer.substr(0, spaceIdxs[0]);
    if (chr != region.GetChrom())
      throw runtime_error("Found site on chromosome '" + chr +
                          "' when chromosome '" + region.GetChrom() +
                          "' was expected");

    // move ahead and extract position from third field
    size_t endReadIdx = 0;
    unsigned pos = stoul(
        buffer.substr(spaceIdxs[1] + 1, spaceIdxs[2] - (spaceIdxs[1] + 1)),
        &endReadIdx, 0);

    // make sure the whole field was parsed!
    if (endReadIdx != spaceIdxs[2] - (spaceIdxs[1] + 1))
      throw runtime_error("Input haplotypes line field three was read as " +
                          sutils::uint2str(pos) +
                          " and does not seem to be an unsigned integer.\n" +
                          "Read was " + sutils::uint2str(endReadIdx) +
                          " char(s) long" + "\nThe field was actually " +
                          sutils::uint2str(spaceIdxs[2]) + " char(s) long");

    // make sure input sites are sorted by position
    if (pos < previousPos)
      throw runtime_error("Input haplotypes file " + hapsFile +
                          " needs to be sorted by position");
    previousPos = pos;

    // start loading only once we hit the first site
    if (pos < region.GetStartBP())
      continue;

    // stop loading sites if the current site is past the last site position in
    // the GLs
    if (pos > region.GetEndBP())
      break;

    // split entire line for processing
    vector<string> tokens;
    boost::split(tokens, buffer, boost::is_any_of(" "));

    // make sure header start is correct
    if (numHaps == 0) {
      if (tokens.size() <= 5)
        throw runtime_error("haps file " + hapsFile +
                            " contains too few columns (" +
                            sutils::uint2str(tokens.size()) + ")");

      // count number of haps
      numHaps = tokens.size() - 5;

    } else {
      if (tokens.size() - 5 != numHaps)
        throw runtime_error(
            "Every row of haplotypes file must have the same number of "
            "columns");
    }

    // store site of haplotype
    loadSites.push_back(snp(chr, pos, tokens[3], tokens[4]));
    vector<char> loadSite;
    loadSite.reserve(numHaps);

    for (unsigned i = 5; i != tokens.size(); i++) {
      int val = atoi(tokens[i].c_str());

      if (val == 0)
        loadSite.push_back(0);
      else if (val == 1)
        loadSite.push_back(1);
      else {
        throw runtime_error("All alleles are not 0 or 1 In haplotypes file " +
                            hapsFile + " at line number " +
                            sutils::uint2str(lineNum) + " (0-based)");
      }
    }

    loadHaps.push_back(std::move(loadSite));

    // keeping this site
    keptSites++;
  }

  cout << "[HapPanel] Sites kept:\t" << keptSites << " / " << lineNum << "\n";

  if (numHaps == 0) {
    cerr << "[HapPanel] Number of haplotypes in haps file is 0.  Haps file "
            "empty?\n";
    exit(1);
  }
}

void OpenTabHaps(const string &hapsFile, const Region &region,
                 vector<vector<char> > &loadHaps, vector<snp> &loadSites) {
  cout << "[HapPanel] Loading tabixed haps file: " << hapsFile << endl;

  // clear all the containers that are going to be filled up
  vector<snp> fillSites;
  vector<vector<char> > fillHaps;

  string file(hapsFile);
  Tabix tabix(file);
  tabix.setRegion(region.AsString());

  // start parsing haplotype lines
  string buffer;
  unsigned numHaps = 0;
  while (tabix.getNextLine(buffer)) {

    // split on tab
    vector<string> tokens;
    boost::split(tokens, buffer, boost::is_any_of("\t"));
    if (tokens[0] != region.GetChrom())
      throw std::runtime_error("[insti] Found site on chromosome '" +
                               tokens[0] + "' when chromosome '" +
                               region.GetChrom() + "' was expected");
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

    Bio::snp newSite(tokens[0], pos, tokens[3], tokens[4]);

    // check to make sure the correct number of haplotypes exist in line
    if (tokens.size() < 6)
      throw std::runtime_error("[HapPanel] no haplotypes in haps file: " +
                               hapsFile);
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

vector<snp> OpenLegend(const string &legendFile, string chrom) {

  // read in legend file
  ifile legendFD(legendFile);

  if (!legendFD.isGood())
    throw runtime_error("Could not open legend file [" + legendFile + "]");

  string buffer;
  vector<snp> loadLeg;

  // discard header
  unsigned lineNum = 0;

  while (getline(legendFD, buffer, '\n')) {
    if (buffer.size() == 0)
      throw runtime_error("Error in legend file " + legendFile +
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
          throw runtime_error(
              "Error in legend file " + legendFile +
              ": header start does not match:\n\t" + header +
              "\n.  Instead the first line of the header is:\n\t" + buffer +
              "\n");
      }

      continue;
    }

    // add each site to loadLeg
    loadLeg.push_back(
        snp(chrom, strtoul(tokens[1].c_str(), NULL, 0), tokens[2], tokens[3]));
  }

  return loadLeg;
}

vector<vector<char> > OpenHap(const string &hapFile) {

  // read in the hap file
  ifile hapFD(hapFile);

  if (!hapFD.isGood())
    throw runtime_error("Could not open hap file [" + hapFile + "]");

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
      throw runtime_error(
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
        throw runtime_error("All alleles are not 0 or 1 In haplotypes file " +
                            hapFile + " at line number " +
                            sutils::uint2str(lineNum) + " (0-based)");
      }
    }

    loadHaps.push_back(inSite);
  }

  if (uNumHaps == 0)
    throw runtime_error("num haps is 0.  Haps file empty?\n");

  return loadHaps;
}

vector<string> OpenSamp(const string &sampleFile) {
  ifile sampFD{ sampleFile };
  if (!sampFD.isGood())
    throw runtime_error("Could not open sample file [" + sampleFile + "]");

  string buffer;
  vector<string> tokens;
  // read header
  getline(sampFD, buffer, '\n');
  boost::split(tokens, buffer, boost::is_any_of(" "));
  if (tokens[0] != "sample")
    throw runtime_error("Input sample file does not have a header with first "
                        "column equal to 'sample': " +
                        sampleFile);

  // read body
  vector<string> sampIDs;
  while (getline(sampFD, buffer, '\n')) {
    boost::split(tokens, buffer, boost::is_any_of(" "));
    sampIDs.push_back(move(tokens[0]));
  }
  clog << "[HapPanel] Loaded " << sampIDs.size() << " sample IDs" << endl;
  return sampIDs;
}

void OpenVCFGZ(string vcf, const Region &region,
               vector<vector<char> > &loadHaps, vector<snp> &loadSites,
               vector<string> &ids) {
  loadHaps.clear();
  loadSites.clear();
  ids.clear();

  // open vcf using tabixpp
  Tabix tbx(vcf);
  tbx.setRegion(region.AsString());

  // get the header
  string header{ tbx.getHeader() };

  using qi::omit;
  using qi::lit;

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

    // make sure parser parsed everything
    if (cfirst != line.cend())
      throw std::runtime_error(
          "Malformed VCF line at line " + to_string(siteNum + 1) +
          " of VCF body in file: " + vcf + "\nLine: " + line);

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
}

using namespace HapPanelHelper;

HapPanel::HapPanel(const vector<string> &inFiles, HapPanelHelper::Init init)
    : m_init{ std::move(init), } {

  // push sampleIDs into unordered map
  m_keepSampIDsUnorderedMap.reserve(m_init.keepSampleIDs.size());
  for (size_t idx = 0; idx < m_init.keepSampleIDs.size(); ++idx)
    m_keepSampIDsUnorderedMap.insert(make_pair(m_init.keepSampleIDs[idx], idx));

  if (!(m_init.keepSites.empty() xor m_init.keepRegion.empty()))
    throw runtime_error(
        "Please only specify site list or region for filtering");

  if (!m_init.alwaysKeepVariantsFile.empty())
    LoadAlwaysKeepVariantsFile(m_init.alwaysKeepVariantsFile);

  // handle region/sites to keep
  if (!m_init.keepSites.empty())
    LoadKeepSites(m_init.keepSites);
  else {
    assert(!m_init.matchAllSites);
    assert(m_init.keepAllSitesInRegion);
    assert(m_init.alwaysKeepVariantsFile.empty());
    m_keepRegion = m_init.keepRegion;
  }

  // load haplotypes and samples
  vector<vector<char> > loadHaps;
  vector<snp> loadSites;
  vector<string> loadIDs;
  if (m_init.hapFileType == HapPanelHelper::HapFileType::IMPUTE2) {
    assert(inFiles.size() == 3);
    loadSites =
        std::move(OpenLegend(std::move(inFiles[1]), m_keepRegion.GetChrom()));
    loadHaps = std::move(OpenHap(std::move(inFiles[0])));
    loadIDs = std::move(OpenSamp(std::move(inFiles[2])));
  } else if (m_init.hapFileType == HapPanelHelper::HapFileType::WTCCC) {
    assert(inFiles.size() == 2);
    OpenHaps(move(inFiles[0]), m_keepRegion, loadHaps, loadSites);
    loadIDs = move(OpenSample(inFiles[1]));
  } else if (m_init.hapFileType == HapPanelHelper::HapFileType::WTCCC_TBX) {
    assert(inFiles.size() == 2);
    OpenTabHaps(inFiles[0], m_keepRegion, loadHaps, loadSites);
    loadIDs = move(OpenSample(inFiles[1]));
  } else if (m_init.hapFileType == HapPanelHelper::HapFileType::VCF) {
    assert(inFiles.size() == 1);
    OpenVCFGZ(inFiles[0], m_keepRegion, loadHaps, loadSites, loadIDs);
  } else
    assert(false);

  if (loadHaps.empty())
    throw std::runtime_error("Input files contain no haplotypes: " +
                             inFiles[0]);

  // filter samples
  if (!m_init.keepSampleIDs.empty()) {
    SubsetSamples(loadIDs, loadHaps);
    OrderSamples(loadIDs, loadHaps);
  }

  // filter on sites
  assert(!loadHaps.empty());
  if (!m_init.keepAllSitesInRegion && !m_keepSites.empty())
    FilterSites(loadHaps, loadSites);

  std::unordered_map<Bio::snp, size_t, Bio::snpKeyHasher> sitesUnorderedMap;
  sitesUnorderedMap.reserve(loadSites.size());
  for (size_t i = 0; i < loadSites.size(); ++i)
    sitesUnorderedMap.insert(make_pair(loadSites[i], i));

  std::swap(m_haps, loadHaps);
  std::swap(m_sites, loadSites);
  std::swap(m_sampleIDs, loadIDs);
  std::swap(m_sitesUnorderedMap, sitesUnorderedMap);

  CheckPanel();
};

void HapPanel::SubsetSamples(vector<string> &loadIDs,
                             vector<vector<char> > &loadHaps) {

  if (m_init.keepSampleIDs.empty())
    return;

  assert(loadIDs.size() >= m_init.keepSampleIDs.size());

  if (!loadHaps.empty())
    assert(loadIDs.size() * 2 == loadHaps[0].size());

  vector<unsigned> idIdxsToKeep;

  for (unsigned idIdx = 0; idIdx < loadIDs.size(); idIdx++) {
    auto foundID = m_keepSampIDsUnorderedMap.find(loadIDs[idIdx]);

    if (foundID != m_keepSampIDsUnorderedMap.end())
      idIdxsToKeep.push_back(idIdx);
  }

  // subset haplotypes
  vector<string> tempIDs;
  tempIDs.reserve(idIdxsToKeep.size());

  for (auto keepIdx : idIdxsToKeep)
    tempIDs.push_back(loadIDs[keepIdx]);

  std::swap(tempIDs, loadIDs);

  size_t numSitesLoaded = 0;
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
    ++numSitesLoaded;
  }

  if (!loadHaps.empty())
    assert(loadIDs.size() * 2 == loadHaps[0].size());

  cout << "[HapPanel] Number of loaded sites: " << numSitesLoaded << endl;
  cout << "[HapPanel] Number of haplotypes left after filtering: "
       << loadHaps[0].size() << endl;
}

// Filter sites based on keep sites
void HapPanel::FilterSites(vector<vector<char> > &loadHaps,
                           vector<snp> &loadSites) {

  assert(loadSites.size() > 0);
  assert(loadSites.size() == loadHaps.size());
  assert(loadHaps[0].size() > 0);

  // initializing filtHaps
  //    for(auto & iter  : filtHaps)
  //      iter.resize(loadSites.size());

  // now keep all sites that are in the keepSites
  unsigned previousIdx = 0;
  typedef pair<std::unordered_map<Bio::snp, HapPanelHelper::siteMeta,
                                  Bio::snpKeyHasher>::const_iterator,
               size_t> snpIterator_index_pair;
  vector<snpIterator_index_pair> keepSitesIterators;
  for (unsigned loadSiteIdx = 0; loadSiteIdx < loadSites.size();
       ++loadSiteIdx) {

    auto got = m_keepSites.find(loadSites[loadSiteIdx]);
    if (got == m_keepSites.end())
      continue;
    if (!keepSitesIterators.empty()) {
      if (got->second.index <= previousIdx)
        if (got->first.pos != loadSites[loadSiteIdx].pos ||
            !m_init.allowReorderOfSitesAtPos) {
          Bio::snp &site = loadSites[loadSiteIdx];
          throw std::runtime_error("Haplotypes file is not ordered in same "
                                   "order as site list at position " +
                                   site.chr + ":" + to_string(site.pos) + " " +
                                   site.ref + " " + site.alt);
        }
      previousIdx = got->second.index;
    }
    keepSitesIterators.push_back(make_pair(got, loadSiteIdx));
  }

  // sort the load site indices by the order that the sites appear in
  // m_keepSites
  sort(keepSitesIterators.begin(), keepSitesIterators.end(),
       [](const snpIterator_index_pair &a, const snpIterator_index_pair &b) {
    return a.first->second.index < b.first->second.index;
  });

  // move sites to keep from load sites into filtered sites
  vector<snp> filtSites;
  vector<vector<char> > filtHaps;
  for (size_t keepSiteIteratorIdx = 0;
       keepSiteIteratorIdx < keepSitesIterators.size(); ++keepSiteIteratorIdx) {

    assert(loadSites[keepSitesIterators[keepSiteIteratorIdx].second] ==
           keepSitesIterators[keepSiteIteratorIdx].first->first);

    // keep the sites that match the sites in m_keepSites
    filtSites.push_back(
        std::move(loadSites[keepSitesIterators[keepSiteIteratorIdx].second]));
    filtHaps.push_back(
        std::move(loadHaps[keepSitesIterators[keepSiteIteratorIdx].second]));
  }

  cout << "[HapPanel] Number of candidate sites: " << loadSites.size() << endl;

  std::swap(filtSites, loadSites);
  std::swap(filtHaps, loadHaps);

  cout << "[HapPanel] Number of sites kept: " << loadSites.size() << endl;

  if (keepSitesIterators.empty())
    throw std::runtime_error("Error: No sites in loaded panel matched "
                             "already loaded sites. Is your "
                             "input file sorted?");

  // we always want to find at most as many matches as are in GL sites
  if (keepSitesIterators.size() > m_keepSites.size())
    throw std::runtime_error("[HapPanel] too many matching sites found. There "
                             "are duplicate sites in input.");
}

void HapPanel::OrderSamples(vector<string> &loadIDs,
                            vector<vector<char> > &loadHaps) {

  assert(loadIDs.size() == m_keepSampIDsUnorderedMap.size());

  if (!loadHaps.empty())
    assert(loadIDs.size() * 2 == loadHaps[0].size());

  vector<size_t> orderedLoadNameIdxs;
  orderedLoadNameIdxs.reserve(loadIDs.size());

  for (unsigned idIdx = 0; idIdx < loadIDs.size(); idIdx++) {
    auto foundID = m_keepSampIDsUnorderedMap.find(loadIDs[idIdx]);

    if (foundID == m_keepSampIDsUnorderedMap.end())
      assert(false);

    assert(foundID->second < m_init.keepSampleIDs.size());
    orderedLoadNameIdxs.push_back(foundID->second);
  }

  // sort names according to orderedLoadNameIdxs
  vector<string> tempIDs(orderedLoadNameIdxs.size());

  for (unsigned nameOrderIdx = 0; nameOrderIdx < orderedLoadNameIdxs.size();
       nameOrderIdx++)
    tempIDs[orderedLoadNameIdxs[nameOrderIdx]] = loadIDs[nameOrderIdx];

  std::swap(tempIDs, loadIDs);
  for (size_t i = 0; i < loadIDs.size(); ++i)
    assert(loadIDs[i] == m_init.keepSampleIDs.at(i));

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

void HapPanel::CheckPanel() const {

  if (m_init.matchAllSites)
    if (m_sites.size() != m_keepSites.size())
      throw runtime_error("Input file should match all sites in keep list (" +
                          to_string(m_keepSites.size()) +
                          "), but it only matches a subset of those sites (" +
                          to_string(m_sites.size()) + ")");

  assert(!m_haps.empty());
  if (!m_init.keepSampleIDs.empty())
    MatchSamples(m_init.keepSampleIDs, m_sampleIDs, m_haps[0].size());
}

void HapPanel::FilterSitesOnAlleleFrequency(double lowerBound,
                                            double upperBound,
                                            bool useAlternateAlleleFrequency) {

  clog << "[HapPanel] Filtering on allele frequency [" + to_string(lowerBound) +
              "," + to_string(upperBound) + "]" << endl
       << "\tuse alternate allele frequency = " << useAlternateAlleleFrequency
       << endl;
  assert(lowerBound >= 0);
  assert(upperBound <= 1);
  assert(lowerBound <= upperBound);

  // calculate allele frequencies
  vector<double> alleleFrequencies;
  alleleFrequencies.reserve(m_sites.size());
  for (auto siteAlleles : m_haps)
    alleleFrequencies.push_back(
        std::count(siteAlleles.begin(), siteAlleles.end(), 1) /
        static_cast<double>(siteAlleles.size()));

  // convert allele frequencies to the correct scale
  if (!useAlternateAlleleFrequency)
    for (auto &af : alleleFrequencies)
      af = abs(af - 0.5) * -1 + 0.5;

  // move sites that have correct allele frequency or are marked as
  // "alwaysKeep==true" to new filtered list of sites and haplotypes
  vector<vector<char> > filteredHaps;
  vector<snp> filteredSites;
  {
    auto got = m_alwaysKeepSites.end();
    for (size_t i = 0; i < alleleFrequencies.size(); ++i) {

      // logic for figuring out if site is in the always keep list
      bool alwaysKeep = false;
      if (!m_alwaysKeepSites.empty()) {
        got = m_alwaysKeepSites.find(m_sites[i]);
        if (got != m_alwaysKeepSites.end())
          alwaysKeep = true;
      }

      // move site across if it is to be kept
      if (alwaysKeep == true || (alleleFrequencies[i] >= lowerBound &&
                                 alleleFrequencies[i] <= upperBound)) {
        filteredHaps.push_back(m_haps[i]);
        filteredSites.push_back(m_sites[i]);
      }
    }
  }

  clog << "[HapPanel] Number of sites kept: " << filteredSites.size() << "/"
       << m_sites.size() << endl;

  if (filteredSites.empty())
    throw runtime_error(
        "[HapPanel] No sites left after filtering on allele frequency");

  // put the filtered haplotypes back in to position
  std::unordered_map<Bio::snp, size_t, Bio::snpKeyHasher> sitesUnorderedMap;
  sitesUnorderedMap.reserve(filteredSites.size());
  for (size_t i = 0; i < filteredSites.size(); ++i)
    sitesUnorderedMap[filteredSites[i]] = i;

  std::swap(m_haps, filteredHaps);
  std::swap(m_sites, filteredSites);
  std::swap(m_sitesUnorderedMap, sitesUnorderedMap);
}

void HapPanel::LoadAlwaysKeepVariantsFile(std::string alwaysKeepVariantsFile) {

  string varFile = "[" + alwaysKeepVariantsFile + "]";
  clog << "[HapPanel] Loading variants to always keep file " << varFile << endl;
  assert(!alwaysKeepVariantsFile.empty());

  ifile variantsFD(alwaysKeepVariantsFile);

  if (!variantsFD.isGood())
    throw runtime_error("[HapPanel] Could not open variants file " + varFile);

  string buffer;
  // read header
  getline(variantsFD, buffer, '\n');
  if (buffer != "CHROM\tPOS\tREF\tALT")
    throw runtime_error(
        "[HapPanel] First line in variants file " + varFile +
        " does not match expected header 'CHROM\tPOS\tREF\tALT'");

  // read body
  vector<string> tokens;
  int lineNum = 0;
  while (getline(variantsFD, buffer, '\n')) {
    ++lineNum;
    boost::split(tokens, buffer, boost::is_any_of("\t"));
    if (tokens.size() != 4)
      throw runtime_error("[HapPanel] Line " + to_string(lineNum) +
                          " of variants file " + varFile +
                          " does not contain 4 columns");

    m_alwaysKeepSites.insert(
        Bio::snp(tokens[0], stoul(tokens[1]), tokens[2], tokens[3]));
  }

  // now set the always keep flags at the appropriate positions if possible
  if (!m_keepSites.empty())
    for (auto aks : m_alwaysKeepSites) {
      auto got = m_keepSites.find(aks);
      if (got != m_keepSites.end())
        got->second.alwaysKeep = true;
    }
}

void HapPanel::LoadKeepSites(const std::vector<Bio::snp> &keepSites) {

  assert(!keepSites.empty());
  // make sure keep sites are sorted on position
  if (!is_sorted(keepSites.begin(), keepSites.end(), Bio::snpPosComp()))
    throw runtime_error("[HapPanel] input sites are not sorted by position");

  // set region to extract
  m_keepRegion = Region(string(keepSites.front().chr), keepSites.front().pos,
                        keepSites.back().pos);

  // create unordered map from keepsites
  // store some meta information as well
  m_keepSites.reserve(keepSites.size());
  for (size_t siteIdx = 0; siteIdx < keepSites.size(); ++siteIdx) {
    HapPanelHelper::siteMeta sMeta;
    sMeta.index = siteIdx;
    if (!m_alwaysKeepSites.empty() &&
        m_alwaysKeepSites.find(keepSites[siteIdx]) != m_alwaysKeepSites.end())
      sMeta.alwaysKeep = true;
    m_keepSites.insert(std::make_pair(keepSites[siteIdx], std::move(sMeta)));
  }
}

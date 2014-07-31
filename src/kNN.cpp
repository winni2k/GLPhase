#include "kNN.hpp"
#include <utility>

using namespace std;
using namespace KNNhelper;

// initialization for kNN clustering
KNN::KNN(unsigned numClust, const std::vector<uint64_t> &haplotypes,
         unsigned numWordsPerHap, unsigned numSites, double freqLB,
         double freqUB, bool usingMAF, gsl_rng *rng, kNNDistT distMetric)
    : Sampler(rng, haplotypes.size() / numWordsPerHap / 2,
              haplotypes.size() / numWordsPerHap),
      m_numWordsPerHap(numWordsPerHap), m_numSites(numSites),
      m_numClusters(numClust), m_distMetric(distMetric), m_freqLB(freqLB),
      m_freqUB(freqUB), m_usingMAF(usingMAF), m_usingScaffold(true) {

    // check bounds 
  if (m_freqLB < 0)
    throw std::range_error("[KNN] freq lower bound is less than 0: " +
                           to_string(m_freqLB));
  if (m_usingMAF)
    if (m_freqUB > 0.5)
      throw std::range_error("[KNN] freq upper bound is greater than 0.5: " +
                             to_string(m_freqUB));
    else if (freqUB > 1)
      throw std::range_error("[KNN] freq upper bound is greater than 0.5: " +
                             to_string(m_freqUB));
  if (m_freqLB > m_freqUB)
    throw std::range_error("[KNN] freq upper bound(" + to_string(m_freqUB) +
                           ") less than lower bound(" + to_string(m_freqLB) +
                           ")");

  assert(haplotypes.size() == m_numHaps * m_numWordsPerHap);

  if (m_numClusters > m_numHaps - 2) {
    m_numClusters = m_numHaps - 2;
    cout << "[KNN]: Too few haplotypes. Reducing number of clusters to "
         << to_string(m_numClusters) << "\n";
  }

  cout << "[KNN]: Initialization start\n";
  cout << "[KNN]: Distance metric = ";
  switch (distMetric) {
  case kNNDistT::hamming:
    cout << "hamming\n";
    break;
  case kNNDistT::tracLen:
    cout << "tract length\n";
    break;
  default:
    throw std::runtime_error("Unexpected distance metric suplied");
  }

  assert(m_numWordsPerHap > 0);

  // make sure uNumSites is a sensible number
  assert(m_numSites > 0);
  assert(m_numWordsPerHap * (WORDMOD + 1) >= m_numSites);

  // dereferencing the vector pointer
  assert(haplotypes.size() % m_numWordsPerHap == 0);
  assert(m_numHaps % 2 == 0);

  if (UsingScaffold())
    CalculateVarAfs(haplotypes);

  CutDownHaps();
  ClusterHaps(haplotypes);

  cout << "[KNN]: Initialization complete\n";
}

unsigned KNN::SampleHap(unsigned uInd) {

  // this should be out of bounds if left unchanged
  unsigned uPropHap = m_numHaps;

  if (UsingScaffold())
    uPropHap =
        m_neighbors[uInd][gsl_rng_uniform_int(m_rng, m_neighbors[uInd].size())];
  else {
    cout << "kNN without a scaffold is not implemented yet" << endl;
    exit(4);
  }

  // make sure the return value is sensible
  assert(floor(uPropHap / 2) != uInd);
  return uPropHap;
}

vector<unsigned> KNN::Neighbors(unsigned uHapNum) {
  std::vector<unsigned> returnHapNums;
  for (auto iter : m_neighbors[floor(uHapNum / 2)])
    returnHapNums.push_back(iter);
  return returnHapNums;
};

void KNN::ClusterHaps(const vector<uint64_t> &haplotypes) {

  if (!UsingScaffold()) {
    cout << "kNN without a scaffold is not implemented yet" << endl;
    exit(4);
  }

  m_neighbors.clear();
  m_neighbors.reserve(m_numHaps / 2);
  unsigned numCommonSites = m_vuCommonSiteNums.size();

  // create vector of comparison haplotypes
  vector<Haplotype> commonSiteHaps;

  for (unsigned hapNum = 0; hapNum < m_numHaps; hapNum++) {
    Haplotype compHap(numCommonSites);
    AssignHap(compHap, haplotypes.data() + hapNum * m_numWordsPerHap);
    commonSiteHaps.push_back(compHap);
  }

  // let user know where we stand
  cout << "[KNN::ClusterHaps]: Finding " << m_numClusters
       << " nearest neighbors for " << m_numHaps << " haplotypes" << endl;

  for (unsigned uSampNum = 0; uSampNum < m_numHaps / 2; ++uSampNum) {

    // assign only the common sites of each hap to new haps
    Haplotype hHapA(numCommonSites);
    Haplotype hHapB(numCommonSites);
    AssignHap(hHapA, haplotypes.data() + uSampNum * 2 * m_numWordsPerHap);
    AssignHap(hHapB, haplotypes.data() + (uSampNum * 2 + 1) * m_numWordsPerHap);

    // calculate neighbors of this sample
    // hap number then distance
    vector<dist> distsHapA(m_numHaps - 2);
    vector<dist> distsHapB(m_numHaps - 2);

#pragma omp parallel for
    for (unsigned hapNum = 0; hapNum < m_numHaps; hapNum++) {
      if (hapNum >> 1 == uSampNum)
        continue;

      dist distHapA;
      dist distHapB;
      distHapA.idx = hapNum;
      distHapB.idx = hapNum;

      // measure distance between the pair of haps using hamming or trac length
      // as distance metric
      if (m_distMetric == kNNDistT::hamming) {
        distHapA.distance = hHapA.HammingDist(commonSiteHaps[hapNum]);
        distHapB.distance = hHapB.HammingDist(commonSiteHaps[hapNum]);
      }

      // the distance is actually the total number of sites minus the max tract
      // length
      else if (m_distMetric == kNNDistT::tracLen) {
        assert(hHapA.MaxTractLen(commonSiteHaps[hapNum]) <= numCommonSites);
        assert(hHapB.MaxTractLen(commonSiteHaps[hapNum]) <= numCommonSites);
        distHapA.distance =
            numCommonSites - hHapA.MaxTractLen(commonSiteHaps[hapNum]);
        distHapB.distance =
            numCommonSites - hHapB.MaxTractLen(commonSiteHaps[hapNum]);
      }

      size_t outHapNum = (hapNum >> 1) > uSampNum ? hapNum - 2 : hapNum;
      distsHapA[outHapNum] = distHapA;
      distsHapB[outHapNum] = distHapB;
    }

    // sort the haplotypes by distance
    sort(distsHapA.begin(), distsHapA.end(), by_dist());
    sort(distsHapB.begin(), distsHapB.end(), by_dist());

    // keep only the k samples with smallest distance
    vector<unsigned> vuDists(2 * ceil(m_numClusters / 2.0));

    unordered_map<unsigned, bool> sampleList;
    for (unsigned idx = 0; idx != ceil(m_numClusters / 2.0); idx++) {
      // insert hapA
      vuDists[2 * idx] = distsHapA[idx].idx;
      sampleList.emplace(distsHapA[idx].idx, true);
    }

    // insert only unique elements
    unsigned hapIdx = 0;
    for (unsigned idx = 0; idx != ceil(m_numClusters / 2.0); idx++) {

      bool inserted = false;
      while (!inserted) {
        auto retVal = sampleList.emplace(distsHapB[hapIdx].idx, true);
        ++hapIdx;
        if (retVal.second)
          inserted = true;
      }

      // insert hapB
      vuDists[2 * idx + 1] = distsHapB[hapIdx - 1].idx;
    }

    // put the closest k in m_neighbors
    m_neighbors.push_back(vuDists);
    //    m_distsHapA.push_back(std::move(distsHapA));
    //    m_distsHapB.push_back(std::move(distsHapB));
  }
}

void KNN::AssignHap(Haplotype &hHap, const uint64_t *oa) {
  unsigned uSiteNum = 0;

  for (auto iHapNum : m_vuCommonSiteNums) {
    hHap.Set(uSiteNum, hHap.TestSite(iHapNum, oa));
    uSiteNum++;
  }
}

void KNN::CutDownHaps() {

  for (unsigned uPosNum = 0; uPosNum < m_numSites; uPosNum++) {
    double AF = m_varAFs[uPosNum];
    if (m_usingMAF)
      AF = abs(abs(AF - 0.5) - 0.5);
    if (AF >= m_freqLB && AF <= m_freqUB)
      m_vuCommonSiteNums.push_back(uPosNum);
  }
}

/*
  calculate the variant allele frequency from scaffold for each position and
  store it in m_varAFs
*/
void KNN::CalculateVarAfs(const vector<uint64_t> &vuScaffoldHaps) {

  m_varAFs.resize(m_numSites);
  Haplotype hTestHap(m_numSites);

  //    cout << "num words per hap: " << m_numWordsPerHap << endl;
  for (unsigned uPosNum = 0; uPosNum < m_numSites; uPosNum++) {

    unsigned uAltAllNum = 0;

    for (unsigned hapNum = 0; hapNum < m_numHaps; hapNum++) {
      if (hTestHap.TestSite(uPosNum,
                            vuScaffoldHaps.data() + hapNum * m_numWordsPerHap))
        uAltAllNum++;
    }

    //        cout << "alternate allele count: " << uAltAllNum << endl;

    m_varAFs[uPosNum] = static_cast<double>(uAltAllNum) / m_numHaps;
  }

  /*    for(unsigned i = 0; i < m_numHaps; i++){
          cout << "hap num: " << i << ": ";
          for (unsigned j = 0; j < m_numSites; j ++){
              cout << hTestHap.TestSite(j, vuScaffoldHaps.data() + i *
     m_numWordsPerHap);
          }
          cout << endl;
      }
  */
}

void KNN::Save(const string &fileName, const vector<string> &sampNames) {

  cerr << "[KNN]: Saving neighbors to prefix " << fileName << ".kNN.neighbors"
       << endl;

  if (sampNames.size() != m_neighbors.size())
    throw std::runtime_error("number of sample names (" +
                             to_string(sampNames.size()) +
                             ") and number of samples with neighbors (" +
                             to_string(m_neighbors.size()) + ") do not match");

  ofile out(fileName + ".kNN.neighbors");
  out << "#Sample\tneighbors\n";
  for (size_t sampNum = 0; sampNum != sampNames.size(); ++sampNum) {
    out << sampNames[sampNum];
    for (size_t neighborNum = 0; neighborNum != m_neighbors[sampNum].size();
         ++neighborNum) {
      out << "\t" << sampNames[m_neighbors[sampNum][neighborNum] / 2];
      if (m_neighbors[sampNum][neighborNum] & 1)
        out << ".hapB";
      else
        out << ".hapA";
    }
    out << "\n";
  }
  out.close();
  /*
    this was for debugging
    // also save distance between every pair of haps
    ofile dists(fileName + ".kNN.dists.gz");
    dists << "#Hap1\tHap2\tdist\n";
    for (size_t sampNum = 0; sampNum != sampNames.size(); ++sampNum) {
      for (size_t distNum = 0; distNum != m_distsHapA[sampNum].size();
           ++distNum) {
        {
          const unsigned idx = m_distsHapA[sampNum][distNum].idx;
          const unsigned d = m_distsHapA[sampNum][distNum].distance;
          dists << sampNames[sampNum] << ".hapA\t" << sampNames[idx / 2];
          if ((idx & 1) == 0)
            dists << ".hapA\t";
          else
            dists << ".hapB\t";
          dists << to_string(d) << "\n";
        }
        {
          const unsigned idx = m_distsHapB[sampNum][distNum].idx;
          const unsigned d = m_distsHapB[sampNum][distNum].distance;
          dists << sampNames[sampNum] << ".hapB\t" << sampNames[idx / 2];
          if ((idx & 1) == 0)
            dists << ".hapA\t";
          else
            dists << ".hapB\t";
          dists << to_string(d) << "\n";
        }
      }
    }
  */
}

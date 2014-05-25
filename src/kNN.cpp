#include "kNN.hpp"
#include <utility>

using namespace std;

// initialization for kNN clustering
KNN::KNN(unsigned numClust, const std::vector<uint64_t> &haplotypes,
         unsigned numWordsPerHap, unsigned numSites, double freqCutoff,
         gsl_rng *rng, kNNDistT distMetric)
    : Sampler(rng, haplotypes.size() / numWordsPerHap / 2,
              haplotypes.size() / numWordsPerHap),
      m_numWordsPerHap(numWordsPerHap), m_numSites(numSites),
      m_numClusters(numClust), m_distMetric(distMetric),
      m_freqCutoff(freqCutoff) {

  assert(haplotypes.size() == m_numHaps * m_numWordsPerHap);

  cout << "[KNN]: Initialization start\n";
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

// returns the neighbors of haplotype hapNum
struct dist {
  unsigned idx;
  unsigned distance;
};

struct by_dist {
  bool operator()(dist const &a, dist const &b) {
    return a.distance < b.distance;
  }
};

void KNN::ClusterHaps(const vector<uint64_t> &haplotypes) {

  if (!UsingScaffold()) {
    cout << "kNN without a scaffold is not implemented yet" << endl;
    exit(4);
  }

  m_neighbors.clear();
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

  for (unsigned uSampNum = 0; uSampNum < m_numHaps / 2; uSampNum++) {

    // assign only the common sites of each hap to new haps
    Haplotype hHapA(numCommonSites);
    Haplotype hHapB(numCommonSites);
    AssignHap(hHapA, haplotypes.data() + uSampNum * 2 * m_numWordsPerHap);
    AssignHap(hHapB, haplotypes.data() + (uSampNum * 2 + 1) * m_numWordsPerHap);

    // calculate neighbors of this sample
    // hap number then distance
    vector<dist> distsHapA;
    vector<dist> distsHapB;

    for (unsigned hapNum = 0; hapNum < m_numHaps; hapNum++) {
      if (floor(hapNum / 2) == uSampNum)
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
      } else if (m_distMetric == kNNDistT::tracLen) {
        distHapA.distance = hHapA.MaxTractLen(commonSiteHaps[hapNum]);
        distHapB.distance = hHapB.MaxTractLen(commonSiteHaps[hapNum]);
      }

      distsHapA.push_back(distHapA);
      distsHapB.push_back(distHapB);
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
        hapIdx++;
        if (retVal.second)
          inserted = true;
      }

      // insert hapB
      vuDists[2 * idx + 1] = distsHapB[hapIdx - 1].idx;
    }

    // put the closest k in m_neighbors
    m_neighbors.push_back(vuDists);
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
    if (m_vdVarAfs[uPosNum] >= m_freqCutoff)
      m_vuCommonSiteNums.push_back(uPosNum);
  }
}

/*
  calculate the variant allele frequency from scaffold for each position and
  store it in m_vdVarAfs
*/
void KNN::CalculateVarAfs(const vector<uint64_t> &vuScaffoldHaps) {

  m_vdVarAfs.resize(m_numSites);
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

    m_vdVarAfs[uPosNum] = static_cast<double>(uAltAllNum) / m_numHaps;
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

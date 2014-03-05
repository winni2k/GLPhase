#include "kNN.hpp"
#include <utility>

using namespace std;

// initialization for kNN clustering
void KNN::init(const vector<uint64_t> &vuHaplotypes, unsigned uNumWordsPerHap,
               unsigned uNumSites, double dFreqCutoff) {

  cout << "[KNN]: Initialization start\n";
  assert(m_bInitialized == false);
  assert(uNumWordsPerHap > 0);

  // make sure uNumSites is a sensible number
  assert(uNumSites > 0);
  assert(uNumWordsPerHap * (h_WordMod + 1) >= uNumSites);
  m_uNumSites = uNumSites;

  m_dFreqCutoff = dFreqCutoff;

  // dereferencing the vector pointer
  assert(vuHaplotypes.size() % uNumWordsPerHap == 0);
  m_uNumHaps = vuHaplotypes.size() / uNumWordsPerHap;
  assert(m_uNumHaps % 2 == 0);
  m_uNumWordsPerHap = uNumWordsPerHap;

  if (UsingScaffold())
    CalculateVarAfs(vuHaplotypes);

  CutDownHaps();
  ClusterHaps(vuHaplotypes);

  m_bInitialized = true;
  cout << "[KNN]: Initialization complete\n";
}

unsigned KNN::SampleHap(unsigned uInd, gsl_rng *rng) {

  // this should be out of bounds if left unchanged
  unsigned uPropHap = m_uNumHaps;

  if (UsingScaffold())
    uPropHap =
        m_neighbors[uInd][gsl_rng_uniform_int(rng, m_neighbors[uInd].size())];
  else {
    cout << "kNN without a scaffold is not implemented yet" << endl;
    exit(4);
  }

  // make sure the return value is sensible
  assert(floor(uPropHap / 2) != uInd);
  return uPropHap;
}

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

void KNN::ClusterHaps(const vector<uint64_t> &vuHaplotypes) {

  if (!UsingScaffold()) {
    cout << "kNN without a scaffold is not implemented yet" << endl;
    exit(4);
  }

  m_neighbors.clear();
  unsigned numCommonSites = m_vuCommonSiteNums.size();

  // create vector of comparison haplotypes
  vector<Haplotype> commonSiteHaps;

  for (unsigned hapNum = 0; hapNum < m_uNumHaps; hapNum++) {
    Haplotype compHap(numCommonSites);
    AssignHap(compHap, vuHaplotypes.data() + hapNum * m_uNumWordsPerHap);
    commonSiteHaps.push_back(compHap);
  }

  // let user know where we stand
  cout << "[KNN::ClusterHaps]: Finding " << m_uNumClusters
       << " nearest neighbors for " << m_uNumHaps / 2 - 1 << " samples" << endl;

  for (unsigned uSampNum = 0; uSampNum < m_uNumHaps / 2; uSampNum++) {

    // assign only the common sites of each hap to new haps
    Haplotype hHapA(numCommonSites);
    Haplotype hHapB(numCommonSites);
    AssignHap(hHapA, vuHaplotypes.data() + uSampNum * 2 * m_uNumWordsPerHap);
    AssignHap(hHapB,
              vuHaplotypes.data() + (uSampNum * 2 + 1) * m_uNumWordsPerHap);

    // calculate neighbors of this sample
    // hap number then distance
    vector<dist> distsHapA;
    vector<dist> distsHapB;

    for (unsigned hapNum = 0; hapNum < m_uNumHaps; hapNum++) {
      if (floor(hapNum / 2) == uSampNum)
        continue;

      dist distHapA;
      dist distHapB;
      distHapA.idx = hapNum;
      distHapB.idx = hapNum;
      distHapA.distance = hHapA.HammingDist(commonSiteHaps[hapNum]);
      distHapB.distance = hHapB.HammingDist(commonSiteHaps[hapNum]);

      distsHapA.push_back(distHapA);
      distsHapB.push_back(distHapB);
    }

    // sort the haplotypes by distance
    sort(distsHapA.begin(), distsHapA.end(), by_dist());
    sort(distsHapB.begin(), distsHapB.end(), by_dist());

    // keep only the k samples with smallest distance
    vector<unsigned> vuDists(2 * ceil(m_uNumClusters / 2.0));

    unordered_map<unsigned, bool> sampleList;
    for (unsigned idx = 0; idx != ceil(m_uNumClusters / 2.0); idx++) {
      // insert hapA
      vuDists[2 * idx] = distsHapA[idx].idx;
      sampleList.emplace(distsHapA[idx].idx, true);
    }

    // insert only unique elements
    unsigned hapIdx = 0;
    for (unsigned idx = 0; idx != ceil(m_uNumClusters / 2.0); idx++) {

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

  for (unsigned uPosNum = 0; uPosNum < m_uNumSites; uPosNum++) {
    if (m_vdVarAfs[uPosNum] >= m_dFreqCutoff)
      m_vuCommonSiteNums.push_back(uPosNum);
  }
}

/*
  calculate the variant allele frequency from scaffold for each position and
  store it in m_vdVarAfs
*/
void KNN::CalculateVarAfs(const vector<uint64_t> &vuScaffoldHaps) {

  m_vdVarAfs.resize(m_uNumSites);
  Haplotype hTestHap(m_uNumSites);

  //    cout << "num words per hap: " << m_uNumWordsPerHap << endl;
  for (unsigned uPosNum = 0; uPosNum < m_uNumSites; uPosNum++) {

    unsigned uAltAllNum = 0;

    for (unsigned hapNum = 0; hapNum < m_uNumHaps; hapNum++) {
      if (hTestHap.TestSite(uPosNum,
                            vuScaffoldHaps.data() + hapNum * m_uNumWordsPerHap))
        uAltAllNum++;
    }

    //        cout << "alternate allele count: " << uAltAllNum << endl;

    m_vdVarAfs[uPosNum] = static_cast<double>(uAltAllNum) / m_uNumHaps;
  }

  /*    for(unsigned i = 0; i < m_uNumHaps; i++){
          cout << "hap num: " << i << ": ";
          for (unsigned j = 0; j < m_uNumSites; j ++){
              cout << hTestHap.TestSite(j, vuScaffoldHaps.data() + i *
     m_uNumWordsPerHap);
          }
          cout << endl;
      }
  */
}

#include "kMedoids.hpp"

using namespace std;

// initialization for medoids clustering
KMedoids::KMedoids(unsigned uClusterType, unsigned uNumClust,
                   const std::vector<uint64_t> &pvuHaplotypes,
                   unsigned uNumWordsPerHap, unsigned uNumSites, gsl_rng *rng)
    : Sampler(rng, 0, pvuHaplotypes.size() / uNumWordsPerHap),
      m_uNumWordsPerHap(uNumWordsPerHap), m_uNumSites(uNumSites),
      m_uClusterType(uClusterType){
      
  m_vuMedoidHapNum.resize(uNumClust);

  assert(m_vuMedoidHapNum.size() > 0);
  assert(m_uNumWordsPerHap > 0);
  assert(m_uNumSites > 0);

  // dereferencing the vector pointer
  assert(pvuHaplotypes.size() % m_uNumWordsPerHap == 0);

  // actually using the k-medoids algorithm for clustering
  // initialize with random medoids -- they may be the same haplotype
  m_vuHapMedNum.resize(m_numHaps, 0);
  m_vuHapHammingDist.resize(m_numHaps, 0);

  for (unsigned uClustNum = 0; uClustNum < m_vuMedoidHapNum.size(); uClustNum++)
    m_vuMedoidHapNum[uClustNum] = gsl_rng_uniform_int(rng, m_numHaps);

  // find closest medoid for each haplotype
  AssignHapsToBestMedoids(pvuHaplotypes);

  // Now update medoids to best place;
  UpdateMedoids(pvuHaplotypes);
}

void KMedoids::InputTesting(const vector<uint64_t> &pvuHaplotypes) {
  assert(pvuHaplotypes.size() % m_uNumWordsPerHap == 0);
  assert(pvuHaplotypes.size() % m_vuHapMedNum.size() == 0);
}

void KMedoids::AssignHapsToBestMedoids(const vector<uint64_t> &pvuHaplotypes) {
  InputTesting(pvuHaplotypes);
  cerr << "\t    Assigning haplotypes to medoids..\n";

  Haplotype hTester(m_uNumSites);

  for (unsigned uHapNum = 0; uHapNum < m_vuHapMedNum.size(); uHapNum++) {

    // initialize to first medoid
    unsigned uBestMedoidNum = 0;
    unsigned uHamming = hTester.HammingDist(
        &pvuHaplotypes[uHapNum * m_uNumWordsPerHap],
        &pvuHaplotypes[m_vuMedoidHapNum[uBestMedoidNum] * m_uNumWordsPerHap]);

    // also skip first medoid, cause we already did that
    for (unsigned uMedNum = 1; uMedNum < m_vuMedoidHapNum.size(); uMedNum++) {

      // check to see if this medoid is closer to hap. If so, keep it as the
      // best guess
      unsigned uNewHamming = hTester.HammingDist(
          &pvuHaplotypes[uHapNum * m_uNumWordsPerHap],
          &pvuHaplotypes[m_vuMedoidHapNum[uMedNum] * m_uNumWordsPerHap]);

      if (uNewHamming < uHamming) {
        uHamming = uNewHamming;
        uBestMedoidNum = uMedNum;
      }
    }

    m_vuHapMedNum[uHapNum] = uBestMedoidNum;
    m_vuHapHammingDist[uHapNum] = uHamming;
  }
}

void KMedoids::UpdateMedoids(const vector<uint64_t> &pvuHaplotypes) {

  cerr << "\tUpdating medoids..\n";
  InputTesting(pvuHaplotypes);

  double dBestLoss;

  if (m_uClusterType == 1)
    dBestLoss = MedoidLoss(pvuHaplotypes, 1.0);
  else
    dBestLoss = MedoidLoss(pvuHaplotypes, 2.0);

  double dLastLoss = dBestLoss + 1;

  //// permute medoid locations until loss change is smaller than m_dDelta
  while (abs(dLastLoss - dBestLoss) > m_dDelta) {

    dLastLoss = dBestLoss;
    vector<unsigned> vuPrevMedoidHapNum(m_vuMedoidHapNum);

    // update all medoids once
    for (unsigned uMedNum = 0; uMedNum < m_vuMedoidHapNum.size(); uMedNum++) {

      //            cerr << "\t    Permuting medoid locations for medoid:\t" <<
      // uMedNum << "\r";

      if (m_uClusterType == 0)
        dBestLoss = UpdateMedoidPAM(pvuHaplotypes, dBestLoss, uMedNum);
      else if (m_uClusterType == 1)
        UpdateMedoidParkJun(pvuHaplotypes, 1.0);
    }

    // reassign haps to medoids
    AssignHapsToBestMedoids(pvuHaplotypes);

    // calculate total loss
    if (m_uClusterType == 1)
      dBestLoss = MedoidLoss(pvuHaplotypes, 1.0);
    else
      dBestLoss = MedoidLoss(pvuHaplotypes, 2.0);

    // revert back if best loss is larger than last loss
    if (dLastLoss <= dBestLoss) {
      dBestLoss = dLastLoss;
      m_vuMedoidHapNum = vuPrevMedoidHapNum;
    }

    // let's see difference between cluster and actual loss
    cerr << "\t    Best loss: " << dBestLoss << "; Last loss: " << dLastLoss
         << endl;
  }

  cerr << "\n\tUpdating complete.\n";
}

double KMedoids::UpdateMedoidPAM(const vector<uint64_t> &pvuHaplotypes,
                                 double dBestLoss, unsigned uMedNum) {

  // pick a haplotype
  unsigned uOrigMedHapNum = m_vuMedoidHapNum[uMedNum];

  for (unsigned uHapNum = 0; uHapNum < m_vuHapMedNum.size(); uHapNum++) {

    // only look at moves where the original medoid hap
    // and target hap are not the same
    if (uOrigMedHapNum == uHapNum)
      continue;

    // try moving medoid to new hap
    unsigned uPrevMedHapNum = m_vuMedoidHapNum[uMedNum];
    m_vuMedoidHapNum[uMedNum] = uHapNum;
    double dLoss = MedoidLoss(pvuHaplotypes, 2.0);

    // update loss if successful
    // revert to previous configuration if not
    if (dLoss <= dBestLoss)
      m_vuMedoidHapNum[uMedNum] = uPrevMedHapNum;
    else
      dBestLoss = dLoss;
  }

  return dBestLoss;
}

double KMedoids::UpdateMedoidParkJun(const vector<uint64_t> &pvuHaplotypes,
                                     unsigned uMedNum) {

  // pick a haplotype
  unsigned uOrigMedHapNum = m_vuMedoidHapNum[uMedNum];
  vector<unsigned> vuMedHaps;
  unsigned uHapNum = 0;

  // pull out the haps in around this medoid
  for (auto iHapMedNum : m_vuHapMedNum) {
    if (iHapMedNum == uMedNum)
      vuMedHaps.push_back(uHapNum);

    uHapNum++;
  }

  // figure out what current loss is
  double dBestLoss = MedoidLoss(pvuHaplotypes, vuMedHaps, 1.0);

  // only look at haps around medoid
  for (auto iMedHap : vuMedHaps) {

    // only look at moves where the original medoid hap
    // and target hap are not the same
    if (uOrigMedHapNum == iMedHap)
      continue;

    // try moving medoid to new hap
    unsigned uPrevMedHapNum = m_vuMedoidHapNum[uMedNum];
    m_vuMedoidHapNum[uMedNum] = iMedHap;
    double dLoss = MedoidLoss(pvuHaplotypes, vuMedHaps, 1.0);

    // update loss if successful
    // revert to previous configuration if not
    if (dLoss <= dBestLoss)
      m_vuMedoidHapNum[uMedNum] = uPrevMedHapNum;
    else
      dBestLoss = dLoss;
  }

  return dBestLoss;
}

double KMedoids::MedoidLoss(const vector<uint64_t> &pvuHaplotypes,
                            const std::vector<unsigned> &vuMedHaps,
                            double dPower) {

  // compute squared hamming distance between every haplotype and its medoid
  Haplotype hTester(m_uNumSites);
  double dLoss = 0;

  // perform loss across all medoids
  if (vuMedHaps.size() == 0)
    for (unsigned uHapNum = 0; uHapNum < m_vuHapMedNum.size(); uHapNum++)
      dLoss += pow(hTester.HammingDist(
                       &pvuHaplotypes[uHapNum * m_uNumWordsPerHap],
                       &pvuHaplotypes[m_vuMedoidHapNum[m_vuHapMedNum[uHapNum]] *
                                      m_uNumWordsPerHap]),
                   dPower);
  else {
    unsigned uMedoidHapNum = m_vuMedoidHapNum[m_vuHapMedNum[vuMedHaps[0]]];

    for (auto uHapNum : vuMedHaps)
      dLoss += pow(hTester.HammingDist(
                       &pvuHaplotypes[uHapNum * m_uNumWordsPerHap],
                       &pvuHaplotypes[uMedoidHapNum * m_uNumWordsPerHap]),
                   dPower);
  }

  return dLoss;
}

unsigned KMedoids::SampleHap(unsigned uInd) {

  // constrain cols by using clustering
  vector<unsigned> vuHapsOfInterest;

  unsigned uIndHap1Medoid = m_vuHapMedNum[uInd * 2];
  unsigned uIndHap2Medoid = m_vuHapMedNum[uInd * 2 + 1];
  unsigned uHapCounter = 0;

  for (auto iHapMedNum : m_vuHapMedNum) {
    if (iHapMedNum == uIndHap1Medoid || iHapMedNum == uIndHap2Medoid)
      vuHapsOfInterest.push_back(uHapCounter);

    uHapCounter++;
  }

  assert(vuHapsOfInterest.size() > 2); // this avoids nasty business of only
                                       // sampling haps from the individual of
                                       // interest

  unsigned uPropHap =
      m_numHaps; // this should be out of bounds if left unchanged

  while (1) {
    uPropHap = gsl_rng_uniform_int(m_rng, vuHapsOfInterest.size());

    if (vuHapsOfInterest[uPropHap] / 2 != uInd)
      break;
  }

  // make sure the return value is sensible
  assert(uPropHap < m_numHaps);
  return uPropHap;
}

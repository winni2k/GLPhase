#include "kMedoids.h"

using namespace std;
// initialization for medoids clustering
void KMedoids::init(const vector< uint64_t > * pvuHaplotypes, unsigned uNumWordsPerHap, unsigned uNumSites, gsl_rng *rng) {

    assert(m_bInitialized == false);
    assert(m_vuMedoidHapNum.size() > 0);
    assert(uNumWordsPerHap > 0);

    assert(uNumSites > 0);
    m_uNumSites = uNumSites;
    
    // dereferencing the vector pointer
    assert((*pvuHaplotypes).size() % uNumWordsPerHap == 0);
    m_uNumHaps = (*pvuHaplotypes).size() / uNumWordsPerHap;
    m_uNumWordsPerHap = uNumWordsPerHap;
    
    // actually using the k-medoids algorithm for clustering
    // initialize with random medoids -- they may be the same haplotype
    m_vuHapMedNum.resize(m_uNumHaps, 0);
    m_vuHapHammingDist.resize(m_uNumHaps, 0);
    for( unsigned uClustNum = 0; uClustNum < m_vuMedoidHapNum.size(); uClustNum ++){
        m_vuMedoidHapNum[uClustNum] = gsl_rng_uniform_int(rng, m_uNumHaps);
    }

    // initialize all haplotypes as being closest to the first medoid
    for( unsigned uHapNum = 0; uHapNum < m_vuHapMedNum.size(); uHapNum++)
        m_vuHapMedNum[uHapNum] = m_vuMedoidHapNum[0];

    m_bInitialized = true;

    // assign each haplotype to the closest medoid
    UpdateMedoids(pvuHaplotypes);

}

void KMedoids::InputTesting(const vector< uint64_t > * pvuHaplotypes ){
    assert(m_bInitialized == true);    
    assert( (*pvuHaplotypes).size() % m_uNumWordsPerHap == 0 );
    assert( (*pvuHaplotypes).size() % m_vuHapMedNum.size() == 0 );
}

void KMedoids::AssignHapsToBestMedoids(const vector< uint64_t > * pvuHaplotypes ){
    InputTesting(pvuHaplotypes);
    cerr << "\t    Assigning haplotypes to medoids..\n";

    Haplotype hTester(m_uNumSites);
    for( unsigned uHapNum = 0; uHapNum < m_vuHapMedNum.size(); uHapNum++){

        // initialize to first medoid
        unsigned uBestMedoidNum = 0;
        unsigned uHamming = hTester.HammingDist(&( *pvuHaplotypes)[uHapNum], &( *pvuHaplotypes)[m_vuMedoidHapNum[uBestMedoidNum] * m_uNumWordsPerHap]);

        // also skip first medoid, cause we already did that
        for( unsigned uMedNum = 1; uMedNum < m_vuMedoidHapNum.size(); uMedNum ++){

            // check to see if this medoid is closer to hap. If so, keep it as the best guess
            unsigned uNewHamming = hTester.HammingDist(&( *pvuHaplotypes)[uHapNum], &( *pvuHaplotypes)[m_vuMedoidHapNum[uMedNum] * m_uNumWordsPerHap]);
            if(uNewHamming < uHamming){
                uHamming = uNewHamming;
                uBestMedoidNum = uMedNum;
            }
        }

        m_vuHapMedNum[uHapNum] = uBestMedoidNum;
        m_vuHapHammingDist[uHapNum] = uHamming;
    }
}

void KMedoids::UpdateMedoids(const vector< uint64_t > * pvuHaplotypes ){

    cerr << "\tUpdating medoids..\n";
    InputTesting(pvuHaplotypes);        

    // first figure out which medoid is closest to each haplotype
    AssignHapsToBestMedoids(pvuHaplotypes);

    double dBestLoss = MedoidLoss(pvuHaplotypes);
    double dLastLoss = dBestLoss+1;
    
    //// permute medoid locations until loss change is smaller than m_dDelta    
    while(dLastLoss - dBestLoss > m_dDelta){

        dLastLoss = dBestLoss;
        
        // pick a medoid
        for( unsigned uMedNum = 0; uMedNum < m_vuMedoidHapNum.size(); uMedNum ++){
            cerr << "\t    Permuting medoid locations for medoid:\t" << uMedNum << "\r";

            if(m_uClusterType == 0)
                dBestLoss = UpdateMedoidPAM(pvuHaplotypes, dBestLoss, uMedNum);
            else if(m_uClusterType == 1)
                dBestLoss = UpdateMedoidParkJun(pvuHaplotypes, dBestLoss, uMedNum);
        }
    }
    cerr << "\n\tUpdating complete.\n";
}

double KMedoids::UpdateMedoidPAM(const vector< uint64_t > * pvuHaplotypes, double dBestLoss, unsigned uMedNum ){

    // pick a haplotype
    unsigned uOrigMedHapNum = m_vuMedoidHapNum[uMedNum];
    for ( unsigned uHapNum = 0; uHapNum < m_vuHapMedNum.size(); uHapNum ++){
            
        // only look at moves where the original medoid hap
        // and target hap are not the same
        if ( uOrigMedHapNum == uHapNum ) continue;

        // try moving medoid to new hap
        unsigned uPrevMedHapNum = m_vuMedoidHapNum[uMedNum];            
        m_vuMedoidHapNum[uMedNum] = uHapNum;
        double dLoss = MedoidLoss(pvuHaplotypes);

        // update loss if successful
        // revert to previous configuration if not
        if(dLoss <= dBestLoss)
            m_vuMedoidHapNum[uMedNum] = uPrevMedHapNum;
        else
            dBestLoss = dLoss;
    }
    return dBestLoss;
}

double KMedoids::UpdateMedoidParkJun(const vector< uint64_t > * pvuHaplotypes, double dBestLoss, unsigned uMedNum ){
    // pick a haplotype
    unsigned uOrigMedHapNum = m_vuMedoidHapNum[uMedNum];
    vector< unsigned > vuMedHaps;
    unsigned uHapNum = 0;

    // pull out the haps in around this medoid
    for (auto iHapMedNum : m_vuHapMedNum){
        if(iHapMedNum == uMedNum)
            vuMedHaps.push_back(uHapNum);
        uHapNum ++;
    }

    // only look at haps around medoid
    for ( auto iMedHap : vuMedHaps){
            
        // only look at moves where the original medoid hap
        // and target hap are not the same
        if ( uOrigMedHapNum == iMedHap ) continue;

        // try moving medoid to new hap
        unsigned uPrevMedHapNum = m_vuMedoidHapNum[uMedNum];            
        m_vuMedoidHapNum[uMedNum] = iMedHap;
        double dLoss = MedoidLoss(pvuHaplotypes);

        // update loss if successful
        // revert to previous configuration if not
        if(dLoss <= dBestLoss)
            m_vuMedoidHapNum[uMedNum] = uPrevMedHapNum;
        else
            dBestLoss = dLoss;
    }
    return dBestLoss;

}

double KMedoids::MedoidLoss(const vector< uint64_t > * pvuHaplotypes ){

    // compute squared hamming distance between every haplotype and its medoid
    Haplotype hTester(m_uNumSites);
    double dLoss = 0;
    for( unsigned uHapNum = 0; uHapNum < m_vuHapMedNum.size(); uHapNum ++)
        dLoss += pow(
            hTester.HammingDist(
                &( *pvuHaplotypes)[uHapNum * m_uNumWordsPerHap],
                &( *pvuHaplotypes)[m_vuMedoidHapNum[m_vuHapMedNum[uHapNum]] * m_uNumWordsPerHap]
                ),
            2.0);

    return dLoss;
}

unsigned KMedoids::SampleHap(unsigned uInd, gsl_rng *rng){
        
    // constrain cols by using clustering
    vector < unsigned > vuHapsOfInterest;

    unsigned uIndHap1Medoid = m_vuHapMedNum[uInd * 2];
    unsigned uIndHap2Medoid = m_vuHapMedNum[uInd * 2 + 1];
    unsigned uHapCounter = 0;
    for(auto iHapMedNum: m_vuHapMedNum){
        if( iHapMedNum == uIndHap1Medoid || iHapMedNum == uIndHap2Medoid)
            vuHapsOfInterest.push_back(uHapCounter);
        uHapCounter ++;
    }
    assert(vuHapsOfInterest.size() > 2); // this avoids nasty business of only sampling haps from the individual of interest

    unsigned uPropHap = m_uNumHaps; // this should be out of bounds if left unchanged
    while(1){
        uPropHap = gsl_rng_uniform_int(rng, vuHapsOfInterest.size());
        if ( vuHapsOfInterest[uPropHap] / 2 != uInd ) break;
    }
    // make sure the return value is sensible
    assert(uPropHap < m_uNumHaps);
    return uPropHap;

}


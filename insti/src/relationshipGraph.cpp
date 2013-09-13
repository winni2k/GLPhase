#include "relationshipGraph.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

void RelationshipGraph::init(int iGraphType, unsigned uSamples, unsigned uHaplotypes){

    // only allow initialization once. Construct new object otherwise.
    assert(m_bInitialized == false);
    
    m_bUsingRelMat = (iGraphType < 2);

    m_iGraphType = iGraphType;

    // make sure we have at least as many haplotypes as samples
    assert(uSamples * 2  <= uHaplotypes);

    // figure out how many rows and columns our matrix is going to have
    m_uRows = uSamples;
    m_uCols = 0;

    switch (m_iGraphType){
    case 0:
        m_uCols = ceil(uHaplotypes/2);
        m_bUsingHaps = false;
        break;
    case 1:
        m_uCols = uHaplotypes;
        m_bUsingHaps = true;
        break;
    case 2:
        m_uCols = ceil(uHaplotypes/2);
        m_bUsingHaps = false;
        break;
    default:
        std::cerr << "Unknown graph type selected: " << m_iGraphType << std::endl;
        assert(false);
    }

    if(m_bUsingRelMat){
        // resize the relationship matrixes to the appropriate size
        // assign all numerators and denominators a count of 1
        m_2dRelationshipMatNum.resize(m_uRows);
        m_2dRelationshipMatDen.resize(m_uRows);
        for(unsigned i = 0; i < m_uRows; i++){
            m_2dRelationshipMatNum[i].resize(m_uCols,1);
            m_2dRelationshipMatDen[i].resize(m_uCols,1);
        }
    }

    m_bInitialized = true;
}

// initialization for medoids clustering
void RelationshipGraph::init(int iGraphType, unsigned uNumClust, const vector< uint64_t > * pvuHaplotypes, unsigned uNumWordsPerHap, unsigned uNumSites, gsl_rng *rng) {

    assert(m_bInitialized == false);
    assert(iGraphType == 3);
    assert(uNumClust > 0);
    assert(uNumWordsPerHap > 0);

    assert(uNumSites > 0);
    m_uNumSites = uNumSites;
    
    // dereferencing the vector pointer
    assert((*pvuHaplotypes).size() % (h_WordMod + 1) == 0);
    unsigned uNumHaps = (*pvuHaplotypes).size() / uNumWordsPerHap;
    m_uNumWordsPerHap = uNumWordsPerHap;
    
    m_bUsingCluster = true;
    m_iGraphType = iGraphType;

    // actually using the k-medoids algorithm for clustering
    // initialize with random medoids -- they may be the same haplotype
    m_vuMedoidHapNum.resize(uNumClust, 0);
    m_vuHapMedNum.resize(uNumHaps, 0);
    m_vuHapHammingDist.resize(uNumHaps, 0);
    for( unsigned uClustNum = 0; uClustNum < uNumClust; uClustNum ++){
        m_vuMedoidHapNum[uClustNum] = gsl_rng_uniform_int(rng, uNumHaps);
    }

    // initialize all haplotypes as being closest to the first medoid
    for( unsigned uHapNum = 0; uHapNum < m_vuHapMedNum.size(); uHapNum++)
        m_vuHapMedNum[uHapNum] = m_vuMedoidHapNum[0];

    m_bInitialized = true;

    // assign each haplotype to the closest medoid
    UpdateMedoids(pvuHaplotypes);

}

void RelationshipGraph::UpdateMedoids(const vector< uint64_t > * pvuHaplotypes ){

    assert(m_bInitialized == true);
    assert(m_bUsingCluster == true);
    assert( (*pvuHaplotypes).size() % m_uNumWordsPerHap == 0 );
    assert( (*pvuHaplotypes).size() % m_vuHapMedNum.size() == 0 );

    // first figure out which medoid is closest to each haplotype
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

    auto dBestLoss = MedoidLoss(pvuHaplotypes);

    //// permute medoid locations
    // pick a medoid
    for( unsigned uMedNum = 0; uMedNum < m_vuMedoidHapNum.size(); uMedNum ++){

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
    }
}


double RelationshipGraph::MedoidLoss(const vector< uint64_t > * pvuHaplotypes ){

    assert(m_bInitialized == true);
    assert(m_bUsingCluster == true);

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

    

// takes sample num and haplotype num as well as graph type
// samples can be updated, while haplotypes can only be copied from
// every sample has two haplotypes
// therefore: any haplotypes that don't match a sample are reference haps

// numerator and denominator of relationship matrix
// numerator = number of accepted proposals
// denominator = number of total proposals
// first index is individual for which proposal was made
// second index is individual from which was copied

// iGraphType:
// -1 = undefined
// 0 = sample/sample graph
// 1 = sample/haplotype graph
// 2 = no graph - all samples are equally related

unsigned RelationshipGraph::Hap2Col(unsigned uHap){

    assert(m_bInitialized == true);
    
    switch (m_iGraphType){
    case 0:
        return(uHap/2);
    case 1:
        return(uHap);
    case 2:
        return(uHap/2);
    default:
        cerr << "Unknown graph type selected: " << m_iGraphType << endl;
        assert(false);
    }
}

unsigned RelationshipGraph::Col2Hap(unsigned uCol){

    assert(m_bInitialized == true);
    
    switch (m_iGraphType){
    case 0:
        return(uCol * 2);
    case 1:
        return(uCol);
    case 2:
        return(uCol * 2);
    default:
        cerr << "Unknown graph type selected: " << m_iGraphType << endl;
        assert(false);
    }
}

// update graph medoids
void RelationshipGraph::UpdateGraph(const vector< uint64_t > * pvuHaplotypes ){

    UpdateMedoids(pvuHaplotypes);
    
}

// update graph probabilistically
void RelationshipGraph::UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float dUpdateProb, gsl_rng *rng){

    assert(m_bInitialized == true);
        
    // don't update the graph if we aren't using it
    if(m_iGraphType == 2 || m_iGraphType == 3)
        return;

    // update relationship graph with probability S
    if(gsl_rng_uniform(rng) < dUpdateProb){
        UpdateGraph( p, bAccepted, uInd);
    }

    
}

// update graph with ratio 1
void RelationshipGraph::UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd){

    assert(m_bInitialized == true);
    
    RelationshipGraph::UpdateGraph( p, bAccepted, uInd, 1.0f);
}


// update graph with weight of fRatio
void RelationshipGraph::UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float fRatio){

    assert(m_bInitialized == true);
    
    // don't update the graph if we aren't using it
    if(m_iGraphType == 2 || m_iGraphType == 3)
        return;

    // update relationship matrix
    for( unsigned i = 0; i<4; i++){
        unsigned uPropHap = p[i];

        // convert to correct index based on cols representing haps vs samples
        unsigned uProp = Hap2Col(uPropHap);

        // update proposal hap or sample for current sample
        m_2dRelationshipMatDen[uInd][uProp] += fRatio;
        if(bAccepted)
            m_2dRelationshipMatNum[uInd][uProp] += fRatio;

        // update current sample for proposal hap or sample if proposal is not
        // from reference panel
        if(  m_2dRelationshipMatDen.size() > uProp ){
            if(m_bUsingHaps){
                m_2dRelationshipMatDen[uProp][Col2Hap(uInd)] += fRatio;
                m_2dRelationshipMatDen[uProp][Col2Hap(uInd)+1] += fRatio;
                if(bAccepted){
                    m_2dRelationshipMatNum[uProp][Col2Hap(uInd)] += fRatio;
                    m_2dRelationshipMatNum[uProp][Col2Hap(uInd)+1] += fRatio;
                }
            }
            else{
                m_2dRelationshipMatDen[uProp][uInd] += fRatio;
                if(bAccepted)
                    m_2dRelationshipMatNum[uProp][uInd] += fRatio;
            }
        }
    }
}

// sample a haplotype based on the relationship graph that does not come from individual uInd
unsigned RelationshipGraph::SampleHap(unsigned uInd, gsl_rng *rng, bool bOnlyFromRef){

    assert(m_bInitialized == true);
    
    if(bOnlyFromRef){
        unsigned uRetHap;

        // no guarantees that cluster will work with kickstart,
        // so I'm backing away from that for now
        if( ! m_bUsingCluster){
            while(1){
                uRetHap = RelationshipGraph::SampleHap( uInd, rng);
                if(uRetHap >= m_uRows * 2)
                    break;
            }
        }
        else {
            uRetHap = RelationshipGraph::SampleHap( uInd, rng);
        }
        return(uRetHap);
    }
    else{
        cerr << "programming error";
        exit(1);
    }
}


// sample a haplotype based on the relationship graph that does not come from individual uInd
unsigned RelationshipGraph::SampleHap(unsigned uInd, gsl_rng *rng){

    assert(m_bInitialized == true);
    
    assert(uInd <= m_uRows - 1);

    // if we are clustering, then constrain the number of haplotypes we can look at
    unsigned uPropHap = Col2Hap(m_uCols); // this hap is out of bounds

    // sample uniformly if no graph has been created
    if(m_iGraphType == 2){        
        while(1){
            // m_uCols is 1 based, but gsl_rng makes 0 based choice
            uPropHap = gsl_rng_uniform_int(rng, m_uCols);
            if ( Hap2Col(uPropHap) != uInd) break;
        }
    }

    // sample uniformly from subset out of clustering
    else if(m_iGraphType == 3){
        
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

        while(1){
            uPropHap = gsl_rng_uniform_int(rng, vuHapsOfInterest.size());
            if ( Hap2Col(vuHapsOfInterest[uPropHap]) != uInd ) break;
        }
    }

    
    // otherwise use the relationship graph
    else{

        vector<float> & vuRelRowNum = m_2dRelationshipMatNum[uInd];
        vector<float> & vuRelRowDen = m_2dRelationshipMatDen[uInd];
        unsigned uPropInd = 0;
        unsigned uProp = 0;
//        unsigned uTryNum = 0;
        while(1){
//            uTryNum ++;
            
            // m_uCols is 1 based, but % makes choice 0 based
            uPropHap = gsl_rng_get(rng) % Col2Hap( m_uCols );
            uPropInd = uPropHap / 2;
            uProp = Hap2Col(uPropHap);

            // resample if individual chosen is the same as
            // current individual uInd
            if(uPropInd == uInd)
                continue;

            assert(vuRelRowNum[uProp] / vuRelRowDen[uProp] > 0 );
            assert(vuRelRowNum[uProp] / vuRelRowDen[uProp] <= 1 );
            // resample if individual does not pass rejection sample
            if( gsl_rng_uniform(rng) <= vuRelRowNum[uProp] / vuRelRowDen[uProp] )
                break;
        }
    }

    // make sure the return value is sensible
    assert(uPropHap < Col2Hap(m_uCols));
    return uPropHap;
}


// based on code from SNPTools::Impute
// name should be a list of sample names that includes reference panel sample names
// fileName should be the basename for output
void    RelationshipGraph::Save(string fileName, const vector<string> & name) {

    assert(m_bInitialized == true);
    
    // don't print graph if we don't have one...
    if(m_iGraphType == 2) return;

    cerr << "Saving Relationship graph to prefix " << fileName << " ..." << endl;
    
    unsigned uExpectedNumNames = m_2dRelationshipMatDen[0].size();
    if(m_bUsingHaps) uExpectedNumNames = ceil(uExpectedNumNames/2);
//    cerr << "name.size() = " << name.size() << endl;
//    cerr << "uExpectedNumNames = " << uExpectedNumNames << endl;
    assert(name.size() == uExpectedNumNames);
    
    string numeratorFile = fileName + ".relGraph.num.gz";
    string denominatorFile = fileName + ".relGraph.den.gz";
    string ratioFile = fileName + ".relGraph.gz";

    typedef std::shared_ptr<ofile> ofile_ptr;
    ofile_ptr numFile(new ofile(numeratorFile));
    ofile_ptr denFile(new ofile(denominatorFile));
    ofile_ptr ratFile(new ofile(ratioFile));


    vector<ofile_ptr> ofiles;
    ofiles.push_back( numFile );
    ofiles.push_back( denFile );
    ofiles.push_back( ratFile );

    //start printing header
    *numFile << "Numerator";
    *denFile << "Denominator";
    *ratFile << "Numerator/Denominator";

    for ( unsigned uFileNum = 0; uFileNum < ofiles.size(); uFileNum++){

        // finish printing header
        for (auto sName: name){
            *ofiles[uFileNum] << "\t" << sName;
            if(m_bUsingHaps) *ofiles[uFileNum] << ".hapA" << "\t" << sName << ".hapB";
        }
        *ofiles[uFileNum] << endl;
        
        // print data rows
        // cycle through samples
        for (unsigned uRowNum = 0; uRowNum < m_uRows; uRowNum++) {

            // print sample name
            *ofiles[uFileNum] << name[uRowNum];

            // print rest of row
            for ( unsigned uColNum = 0; uColNum < m_uCols; uColNum++){
                float fPrintVal = 0;
                switch ( uFileNum ){
                case 0:
                    fPrintVal = m_2dRelationshipMatNum[uRowNum][uColNum];
                    assert(fPrintVal != 0);
                    break;
                case 1:
                    fPrintVal = m_2dRelationshipMatDen[uRowNum][uColNum];
                    assert(fPrintVal != 0);
                    break;
                case 2:
                    fPrintVal =  m_2dRelationshipMatNum[uRowNum][uColNum] /  m_2dRelationshipMatDen[uRowNum][uColNum];
                    assert(fPrintVal != 0);
                    break;
                default:
                    cerr << "Programming error: RelationshipGraph::Save: Unknown file number: " << uFileNum << endl;
                    assert(false);
                }
                *ofiles[uFileNum] << "\t" << fPrintVal;
            }
            *ofiles[uFileNum] << endl;
        }
    }
    cerr << "\t...saving of relationship graph complete." << endl;
}



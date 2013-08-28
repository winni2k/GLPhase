#include "relationshipGraph.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

void RelationshipGraph::init(int iGraphType, unsigned uSamples, unsigned uHaplotypes){

    // only allow initialization once. Construct new object otherwise.
    assert(m_bInitialized == false);
    
    m_bUsingRelMat = (iGraphType != 2);
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

void RelationshipGraph::UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float dUpdateProb, gsl_rng *rng){

    // don't update the graph if we aren't using it
    if(m_iGraphType == 2)
        return;

    // update relationship graph with probability S
    if(gsl_rng_uniform(rng) < dUpdateProb){
        UpdateGraph( p, bAccepted, uInd);
    }

    
}

void RelationshipGraph::UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd){

    RelationshipGraph::UpdateGraph( p, bAccepted, uInd, 1.0f);
}

void RelationshipGraph::UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float fRatio){

    // don't update the graph if we aren't using it
    if(m_iGraphType == 2)
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

    if(bOnlyFromRef){
        unsigned uRetHap;
        while(1){
            uRetHap = RelationshipGraph::SampleHap( uInd, rng);
            if(uRetHap >= m_uRows * 2)
                break;
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

    // don't print graph if we don't have one...
    if(m_iGraphType == 2) return(0);

    assert(uInd <= m_uRows - 1);
    unsigned uPropHap = Col2Hap(m_uCols); // this hap is out of bounds

    // sample uniformly if no graph has been created
    if(m_iGraphType == 2){
        while(1){
            // m_uCols is 1 based, but % makes choice 0 based
            uPropHap = gsl_rng_get(rng) % Col2Hap( m_uCols );
            if ( Hap2Col(uPropHap) != uInd) break;
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
void    RelationshipGraph::Save(string fileName, const vector<string> & name, unsigned uNumSamples) {

    unsigned uExpectedNumNames = m_2dRelationshipMatDen[0].size();
    if(m_bUsingHaps) uExpectedNumNames = ceil(uExpectedNumNames/2);
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
        for (unsigned i = 0; i < uNumSamples; i++){            
            *ofiles[uFileNum] << "\t" << name[i];
            if(m_bUsingHaps) *ofiles[uFileNum] << "\t" << name[i];
        }
        
        // print data rows
        // cycle through samples
        for (unsigned uRowNum = 0; uRowNum < m_uRows; uRowNum++) {

            // print sample name
            *ofiles[uFileNum] << "\t" << name[uRowNum];
            for ( unsigned uColNum = 0; uColNum < m_uCols; uColNum++){
                float fPrintVal = 0;
                switch (uFileNum){
                case 0:
                    fPrintVal = m_2dRelationshipMatNum[uRowNum][uColNum];
                    assert(fPrintVal != 0);
                case 1:
                    fPrintVal = m_2dRelationshipMatDen[uRowNum][uColNum];
                    assert(fPrintVal != 0);
                case 2:
                    fPrintVal =  m_2dRelationshipMatNum[uRowNum][uColNum] /  m_2dRelationshipMatDen[uRowNum][uColNum];
                    assert(fPrintVal != 0);
                default:
                    cerr << "Programming error: RelationshipGraph::Save: Unknown file number:" << uFileNum << endl;
                    assert(false);
                }
                *ofiles[uFileNum] << "\t" << fPrintVal;
            }
            *ofiles[uFileNum] << endl;
        }
    }
}



#include "relationshipGraph.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

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

void RelationshipGraph::UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, double dUpdateProb, gsl_rng *rng){

    // don't update the graph if we aren't using it
    if(m_iGraphType == 2)
        return;

    // update relationship graph with probability S
    if(gsl_rng_uniform(rng) < dUpdateProb){
        UpdateGraph( p, bAccepted, uInd);
    }

    
}

void RelationshipGraph::UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd){

    // don't update the graph if we aren't using it
    if(m_iGraphType == 2)
        return;

    // update relationship matrix
    for( unsigned i = 0; i<4; i++){
        unsigned uPropHap = p[i];

        // convert to correct index based on cols representing haps vs samples
        unsigned uProp = Hap2Col(uPropHap);

        // update proposal hap or sample for current sample
        m_2dRelationshipMatDen[uInd][uProp] ++;
        if(bAccepted)
            m_2dRelationshipMatNum[uInd][uProp] ++;

        // update current sample for proposal hap or sample if proposal is not
        // from reference panel
        if(  m_2dRelationshipMatDen.size() > uProp ){
            if(m_bUsingHaps){
                m_2dRelationshipMatDen[uProp][Col2Hap(uInd)]++;
                m_2dRelationshipMatDen[uProp][Col2Hap(uInd)+1]++;
                if(bAccepted){
                    m_2dRelationshipMatNum[uProp][Col2Hap(uInd)] ++;
                    m_2dRelationshipMatNum[uProp][Col2Hap(uInd)+1] ++;
                }
            }
            else{
                m_2dRelationshipMatDen[uProp][uInd]++;
                if(bAccepted)
                    m_2dRelationshipMatNum[uProp][uInd] ++;
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

        vector<unsigned> & vuRelRowNum = m_2dRelationshipMatNum[uInd];
        vector<unsigned> & vuRelRowDen = m_2dRelationshipMatDen[uInd];
        unsigned uPropInd = 0;
        unsigned uProp = 0;
        unsigned uTryNum = 0;
        while(1){

            uTryNum ++;
            
            // m_uCols is 1 based, but % makes choice 0 based
            uPropHap = gsl_rng_get(rng) % Col2Hap( m_uCols );
            uPropInd = uPropHap / 2;
            uProp = Hap2Col(uPropHap);

            // resample if individual chosen is the same as
            // current individual uInd
            if(uPropInd == uInd)
                continue;

            // resample if individual does not pass rejection sample
            if (vuRelRowNum[uProp] / vuRelRowDen[uProp] < 1){
                cerr << "Try "<< uTryNum <<"\tNumerator: "<< vuRelRowNum[uProp] << "\tDenominator: "<< vuRelRowDen[uProp] <<"\tProposal: " << uProp << "\r";
            }
            if( gsl_rng_uniform(rng) <= vuRelRowNum[uProp] / vuRelRowDen[uProp] )
                if (vuRelRowNum[uProp] / vuRelRowDen[uProp] < 1)
                    cerr << endl;
                break;
        }
    }

    // make sure the return value is sensible
    assert(uPropHap < Col2Hap(m_uCols));
    return uPropHap;
}


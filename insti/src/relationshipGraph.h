/* @(#)relationshipGraph.h
 */

#ifndef _RELATIONSHIPGRAPH_H
#define _RELATIONSHIPGRAPH_H 1

#include        <gsl/gsl_rng.h>
#include        <vector>
#include        <cassert>
#include        <cmath>
#include        <iostream>
#include        "utils.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

class RelationshipGraph{

private:

    // 1 based number of rows and columns
    unsigned m_uRows;
    unsigned m_uCols;

//    unsigned m_bUsingRefPanel;
    bool m_bUsingRelMat = true;

    std::vector< std::vector< float > > m_2dRelationshipMatNum;
    std::vector< std::vector< float > > m_2dRelationshipMatDen;

    // -1 = undefined
    // 0 = sample/sample graph
    // 1 = sample/haplotype graph
    // 2 = no graph - all samples are equally related
    int m_iGraphType = -1;
    bool m_bUsingHaps = false;

    // hap to column index converter
    unsigned Hap2Col(unsigned uHap);

    // column to hap index converter
    unsigned Col2Hap(unsigned uCol);

public:

    // empty constructor makes my life easier sometimes
    RelationshipGraph();
    
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

    // m_bUsingRefPanel(uSamples * 2 < uHaplotypes)
    RelationshipGraph(int iGraphType, unsigned uSamples, unsigned uHaplotypes)
        : m_bUsingRelMat(iGraphType != 2), m_iGraphType(iGraphType) {

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

    }

    // returns a haplotype sampled using the relationship graph
    unsigned SampleHap(unsigned uInd, gsl_rng *rng);

    // returns a haplotype sampled using the relationship graph, but only from the reference haplotypes
    unsigned SampleHap(unsigned uInd, gsl_rng *rng, bool bOnlyFromRef);

    // update graph with proposal
    void UpdateGraph( unsigned * p, bool bAccepted, unsigned uInd);

    // update graph with probability dUpdateProb
    void UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float dUpdateProb, gsl_rng *rng);

    // update graph with number fRatio instead of 1
    void UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float fRatio);

    void Save(std::string fileName);

};

#endif /* _RELATIONSHIPGRAPH_H */


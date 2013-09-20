/* @(#)relationship.h
 */

#ifndef _RELATIONSHIP_H
#define _RELATIONSHIP_H 1

#include        <vector>
#include        "relationshipGraph.hpp"
#include        "kMedoids.hpp"


//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");


class Relationship {

private:

    unsigned m_iGraphType;
    KMedoids m_oKMedoids;
    RelationshipGraph m_oRelGraph;
    
public:

    Relationship(unsigned uNumClusters, unsigned uClusterType)
        : m_oKMedoids(uClusterType, uNumClusters), m_oRelGraph() {};

    Relationship(int iGraphType, unsigned uSamples, unsigned uHaplotypes)
        : m_oKMedoids(0, 0), m_oRelGraph() {
        m_oRelGraph.init(iGraphType, uSamples, uHaplotypes);
    }
    
    Relationship()
        : m_oKMedoids(0, 0), m_oRelGraph() {};

    void init(int iGraphType, unsigned uSamples, unsigned uHaplotypes)
        { m_iGraphType = iGraphType; m_oRelGraph.init(iGraphType, uSamples, uHaplotypes); };

    void init(int iGraphType, const vector< uint64_t > * pvuHaplotypes, unsigned uNumWordsPerHap, unsigned uNumSites, gsl_rng *rng)
        { assert(iGraphType == 3); m_iGraphType = iGraphType;  m_oKMedoids.init(pvuHaplotypes, uNumWordsPerHap, uNumSites, rng); }
    
    // sample from cluster graph
    unsigned SampleHap(unsigned uInd, gsl_rng *rng){
        return SampleHap(uInd, rng, false);
    }

    // returns a haplotype sampled using the relationship graph, but only from the reference haplotypes
    unsigned SampleHap(unsigned uInd, gsl_rng *rng, bool bOnlyFromRef){
        if(bOnlyFromRef){
            assert(m_iGraphType != 3);
            return m_oRelGraph.SampleHap(uInd, rng, bOnlyFromRef);
        }
        else if(m_iGraphType == 3)
            return m_oKMedoids.SampleHap(uInd, rng);
        else
            return m_oRelGraph.SampleHap(uInd, rng);
    }

    // update graph with proposal
    void UpdateGraph( unsigned * p, bool bAccepted, unsigned uInd){ UpdateGraph(p, bAccepted, uInd, 1.0f); };

    // update graph with probability dUpdateProb
//    void UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float dUpdateProb, gsl_rng *rng){ m_oRelGraph.UpdateGraph(p, bAccepted, uInd, dUpdateProb, rng); };

    // update graph with number fRatio instead of 1
    void UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float fRatio){
        if(m_iGraphType == 3)
            return;
        else
            m_oRelGraph.UpdateGraph(p, bAccepted, uInd, fRatio);
    };
        
    // update medoids
    void UpdateGraph( const vector< uint64_t > * pvuHaplotypes ){
        if(m_iGraphType == 3)
            m_oKMedoids.UpdateMedoids(pvuHaplotypes);
        else
            return;
    };

    void Save(std::string fileName, const vector<std::string> & name){
        if(m_iGraphType == 3)
            return;
        else
            m_oRelGraph.Save(fileName, name);
    };
    
};

#endif /* _RELATIONSHIP_H */


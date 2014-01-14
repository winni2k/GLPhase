/* @(#)relationship.hpp
 */

#ifndef _RELATIONSHIP_H
#define _RELATIONSHIP_H 1

#include        <vector>
#include        "relationshipGraph.hpp"
#include        "kMedoids.hpp"
<<<<<<< HEAD
#include        "kNN.hpp"
=======
#include        "kMeans.hpp"
>>>>>>> refs/heads/synced/master


//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");


class Relationship {

  private:

    /*
     -1 = undefined
     0 = sample/sample graph
     1 = sample/haplotype graph
     2 = no graph - all samples are equally related
     3 = kmedoids
<<<<<<< HEAD
     4 = kNN
=======
     4 = kmeans
>>>>>>> refs/heads/synced/master
    */
    int m_graphType = -1;
    KMedoids m_oKMedoids;
    RelationshipGraph m_oRelGraph;
<<<<<<< HEAD
    KNN m_oKNN;
=======
    KMeans m_oKMeans;
>>>>>>> refs/heads/synced/master

  public:

    // call for clustering
    Relationship(unsigned uNumClusters, unsigned uClusterType)
        : m_oKMedoids(uClusterType, uNumClusters), m_oRelGraph(),
<<<<<<< HEAD
          m_oKNN(uNumClusters) {};

    // call for relationship graph - ignore clustering
    Relationship(int iGraphType, unsigned uSamples, unsigned uHaplotypes)
        : m_oKMedoids(0, 0), m_oRelGraph(), m_oKNN(0) {
=======
          m_oKMeans(uNumClusters) {};

    // call for relationship graph - ignore clustering
    Relationship(int iGraphType, unsigned uSamples, unsigned uHaplotypes)
        : m_oKMedoids(0, 0), m_oRelGraph(), m_oKMeans(0) {
>>>>>>> refs/heads/synced/master
        SetIGraphType(iGraphType);
        m_oRelGraph.init(iGraphType, uSamples, uHaplotypes);
    }

    void SetIGraphType(int iGraphType){
        assert(iGraphType < 5);
        assert(iGraphType >= 0);
        assert(m_graphType == -1);
        m_graphType = iGraphType;
    }
    
    // initialize relationship graph
    void init(int iGraphType, unsigned uSamples, unsigned uHaplotypes) {
        SetIGraphType(iGraphType);
        m_oRelGraph.init(iGraphType, uSamples, uHaplotypes);
    };

<<<<<<< HEAD
    // initialize kmedoids or kNN
=======
    // initialize kmedoids or kmeans
>>>>>>> refs/heads/synced/master
    void init(int iGraphType, const vector< uint64_t > & vuHaplotypes,
              unsigned uNumWordsPerHap, unsigned uNumSites, gsl_rng *rng) {

        SetIGraphType(iGraphType);
        switch(iGraphType){
        case 3:
            m_oKMedoids.init(&vuHaplotypes, uNumWordsPerHap, uNumSites, rng);
            break;
        case 4:
<<<<<<< HEAD
            m_oKNN.init(vuHaplotypes, uNumWordsPerHap, uNumSites, -1);
=======
            m_oKMeans.init(vuHaplotypes, uNumWordsPerHap, uNumSites, -1);
>>>>>>> refs/heads/synced/master
            break;
        default:
            cout << "unexpected graph type: " << iGraphType << endl;
            exit(3);
        }
    }

<<<<<<< HEAD
    // initialize kNN with scaffold
=======
    // initialize kmeans with scaffold
>>>>>>> refs/heads/synced/master
    void init(int iGraphType,  const vector< uint64_t > & vuScaffoldHaplotypes, unsigned uNumWordsPerHap, unsigned uNumSites, double dScaffoldFreqCutoff) {

        SetIGraphType(iGraphType);
        assert(iGraphType == 4);
<<<<<<< HEAD
        m_oKNN.init(vuScaffoldHaplotypes, uNumWordsPerHap, uNumSites,
=======
        m_oKMeans.init(vuScaffoldHaplotypes, uNumWordsPerHap, uNumSites,
>>>>>>> refs/heads/synced/master
                       dScaffoldFreqCutoff);
    }

    // sample from cluster graph
    unsigned SampleHap(unsigned uInd, gsl_rng *rng) {
        return SampleHap(uInd, rng, false);
    }

    // returns a haplotype sampled using the relationship graph, but only from the reference haplotypes
    unsigned SampleHap(unsigned uInd, gsl_rng *rng, bool bOnlyFromRef) {
        switch(m_graphType){
        case 3:
            return m_oKMedoids.SampleHap(uInd, rng);
        case 4:
<<<<<<< HEAD
            return m_oKNN.SampleHap(uInd, rng);
=======
            return m_oKMeans.SampleHap(uInd, rng);
>>>>>>> refs/heads/synced/master
        default:
            return m_oRelGraph.SampleHap(uInd, rng, bOnlyFromRef);
        }
    }

    // update graph with proposal
    void UpdateGraph(unsigned * p, bool bAccepted, unsigned uInd) {
        UpdateGraph(p, bAccepted, uInd, 1.0f);
    };

    // update graph with probability dUpdateProb
    //    void UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float dUpdateProb, gsl_rng *rng){ m_oRelGraph.UpdateGraph(p, bAccepted, uInd, dUpdateProb, rng); };

    // update graph with number fRatio instead of 1
    void UpdateGraph(unsigned *p, bool bAccepted, unsigned uInd, float fRatio) {
        if (m_graphType < 3)
            m_oRelGraph.UpdateGraph(p, bAccepted, uInd, fRatio);
    };

    // update medoids
    void UpdateGraph(const vector< uint64_t > * pvuHaplotypes) {
        if (m_graphType == 3)
            m_oKMedoids.UpdateMedoids(pvuHaplotypes);
    };

    void Save(std::string fileName, const vector<std::string> & name) {
        if (m_graphType < 3)
            m_oRelGraph.Save(fileName, name);
    };

};

#endif /* _RELATIONSHIP_H */



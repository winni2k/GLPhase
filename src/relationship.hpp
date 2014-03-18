/* @(#)relationship.hpp
 */

#ifndef _RELATIONSHIP_H
#define _RELATIONSHIP_H 1

#include <vector>
#include "relationshipGraph.hpp"
#include "kMedoids.hpp"
#include "kNN.hpp"
#include "hapPanel.hpp"
#include "enums.hpp"

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");


// this class manages 
class Relationship {

private:
  /*
   -1 = undefined
   0 = sample/sample graph
   1 = sample/haplotype graph
   2 = no graph - all samples are equally related
   3 = kmedoids
   4 = kNN
  */
  RelT m_graphType = RelT::undefined;
  KMedoids m_oKMedoids;
  RelationshipGraph m_oRelGraph;
  KNN m_oKNN;

public:
  // call for clustering
  Relationship(unsigned uNumClusters, unsigned uClusterType)
      : m_oKMedoids(uClusterType, uNumClusters), m_oRelGraph(),
        m_oKNN(uNumClusters) {};

  // call for relationship graph - ignore clustering
  Relationship(RelT graphType, unsigned uSamples, unsigned uHaplotypes)
      : m_oKMedoids(0, 0), m_oRelGraph(), m_oKNN(0) {
    SetGraphType(graphType);
    m_oRelGraph.init(graphType, uSamples, uHaplotypes);
  }

  void SetGraphType(RelT graphType) {
    assert(graphType != RelT::undefined);
    assert(m_graphType == RelT::undefined);
    m_graphType = graphType;
  }

  // initialize relationship graph
  void init(RelT graphType, unsigned uSamples, unsigned uHaplotypes) {
    SetGraphType(graphType);
    m_oRelGraph.init(graphType, uSamples, uHaplotypes);
  };

  // initialize kmedoids or kNN
  void init(RelT graphType, const std::vector<uint64_t> &vuHaplotypes,
            unsigned uNumWordsPerHap, unsigned uNumSites, gsl_rng *rng) {

    SetGraphType(graphType);

    switch (graphType) {
    case RelT::kMedoids:
      m_oKMedoids.init(&vuHaplotypes, uNumWordsPerHap, uNumSites, rng);
      break;

    case RelT::kNN:
      m_oKNN.init(vuHaplotypes, uNumWordsPerHap, uNumSites, -1);
      break;

    default:
      std::cerr << "unexpected graph type " << std::endl;
      exit(3);
    }
  }

  // initialize kNN with scaffold
  void init(RelT graphType, HapPanel &scaffold, double dScaffoldFreqCutoff) {

    assert(graphType == RelT::kNN);
    SetGraphType(graphType);
    m_oKNN.init(*scaffold.Haplotypes(), scaffold.NumWordsPerHap(),
                scaffold.NumSites(), dScaffoldFreqCutoff);
  }

  // sample from cluster graph
  unsigned SampleHap(unsigned uInd, gsl_rng *rng) {
    return SampleHap(uInd, rng, false);
  }

  // returns a haplotype sampled using the relationship graph, but only from the
  // reference haplotypes
  unsigned SampleHap(unsigned uInd, gsl_rng *rng, bool bOnlyFromRef) {
    switch (m_graphType) {
    case RelT::kMedoids:
      return m_oKMedoids.SampleHap(uInd, rng);

    case RelT::kNN:
      return m_oKNN.SampleHap(uInd, rng);

    default:
      return m_oRelGraph.SampleHap(uInd, rng, bOnlyFromRef);
    }
  }

  // update graph with proposal
  void UpdateGraph(unsigned *p, bool bAccepted, unsigned uInd) {
    UpdateGraph(p, bAccepted, uInd, 1.0f);
  };

  // update graph with probability dUpdateProb
  //    void UpdateGraph( unsigned *p, bool bAccepted, unsigned uInd, float
  // dUpdateProb, gsl_rng *rng){ m_oRelGraph.UpdateGraph(p, bAccepted, uInd,
  // dUpdateProb, rng); };

  // update graph with number fRatio instead of 1
  void UpdateGraph(unsigned *p, bool bAccepted, unsigned uInd, float fRatio) {
    switch (m_graphType) {
    case RelT::sampSampGraph:
      ;
    case RelT::sampHapGraph:
      ;
    case RelT::noGraph:
      m_oRelGraph.UpdateGraph(p, bAccepted, uInd, fRatio);
      break;
    default:
      throw "unexpected graphType";
    };
  };

  // update medoids
  void UpdateGraph(const std::vector<uint64_t> *pvuHaplotypes) {
    if (m_graphType == RelT::kMedoids)
      m_oKMedoids.UpdateMedoids(pvuHaplotypes);
  };

  void Save(std::string fileName, const std::vector<std::string> &name) {
    switch (m_graphType) {
    case RelT::sampSampGraph:
      ;
    case RelT::sampHapGraph:
      ;
    case RelT::noGraph:
      m_oRelGraph.Save(fileName, name);
      break;
    default:
      throw "unexpected graphType";
    }
  };
};

#endif /* _RELATIONSHIP_H */

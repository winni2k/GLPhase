#include <globals.h>
#include <objects/chaplotype.h>
#include <containers/genhap_set.h>

// HCHelper = haplotypeClusterHelper
namespace HCHelper {
set<int> sampleIndex(int k, int n);
}
class haplotypeCluster {
public:
  int start, stop, min_cluster_size, max_cluster_size, nhap, ncond;
  genhap_set *g;
  haplotypeCluster(genhap_set *);
  ~haplotypeCluster();
  void build(int a, int b, float c1, float c2, int ncond);
  int mask(int curr_ind, vector<bool> &output);

private:
  int cstart, cstop, nsnp, K;
  float c1, c2;
  int niteration;
  map<unsigned int, vector<int> > clusters;
  vector<unsigned int> clustering;
  vector<bool> cluster_mask;
  vector<vector<vector<float> > > dlook;
  vector<vector<float> > mu;
  int summariseKmeans();
  map<unsigned int, int> cluster_size;

  int hap2float(int idx, vector<float> &output);
  float euc(int idx, vector<float> &mu);
  float euc(int i, int j, vector<vector<float> > &mu,
            vector<vector<vector<float> > > &dlook);
  pair<int, int> kmeans(unsigned int cluster, unsigned int depth);
  int bifurc_kmeans();
  int bifurc_kmeans(unsigned int cluster, int n, int depth);
  int getIdx(int idx);
  int updateMeans(pair<unsigned int, unsigned int> &cluster_id,
                  pair<int, int> &cluster_size);
  int kplusplus(int K, vector<int> &output);
  int kplusplus(int K, vector<int> &input, vector<int> &output);
  int topupClusters();
  int closestN(int N, int clusteridx, vector<pair<int, float> > &output);
};

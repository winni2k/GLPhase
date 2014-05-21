#include <objects/haplotypeCluster.h>

haplotypeCluster::haplotypeCluster(genhap_set *g) { this->g = g; }

haplotypeCluster::~haplotypeCluster() {
  // indexing.clear();
  clustering.clear();
}

void haplotypeCluster::build(int a, int b, float c1, float c2, int ncond) {

  this->c1 = c1;
  this->c2 = c2;
  this->ncond = ncond;
  K = 2;
  min_cluster_size = max((int)c1 * ncond, (int)ncond);
  max_cluster_size = max((int)c2 * ncond, (int)ncond);
  start = a;
  stop = b;
  cstart = 8 * (a / 8);
  cstop = 8 * (b / 8);
  nsnp = (cstop - cstart);
  nhap = g->vecH.size();
  clustering.assign(nhap, 0);
  cluster_mask.resize(nhap);
  niteration = 20;

  if (_DEBUG > 1) {
    cout << "c1 = " << c1 << endl;
    cout << "c2 = " << c2 << endl;
    cout << "min_cluster_size = " << min_cluster_size << endl;
    cout << "max_cluster_size = " << max_cluster_size << endl;
    cout << "nhap = " << nhap << endl;
    cout << "[cstart, cstop] = [" << a << "," << b << "]" << endl;
  }

  assert(a < b);
  assert(c1 <= c2);

  bifurc_kmeans(); // divisive K-means clustering routine - my current
                   // preference
  topupClusters(); // ensures cluster size >= min_cluster_size

  for (map<unsigned int, vector<int> >::iterator it = clusters.begin();
       it != clusters.end(); it++)
    assert(it->second.size() >= min_cluster_size);

  if (_DEBUG > 0)
    summariseKmeans();
}

set<int> sampleIndex(int k, int n) {
  set<int> ret;
  int idx;

  for (int i = 0; i < k; i++) {
    idx = putils::getRandom(n);
    while (ret.count(idx))
      idx = putils::getRandom(n);
    ret.insert(idx);
  }
  return (ret);
}

vector<int> set2vec(set<int> tmp) {
  vector<int> ret(tmp.size());
  int i = 0;
  for (set<int>::iterator idx1 = tmp.begin(); idx1 != tmp.end(); idx1++) {
    ret[i] = *idx1;
    i++;
  }
  return (ret);
}

int haplotypeCluster::mask(int curr_ind, vector<bool> &output) {
  assert(curr_ind < nhap / 2);
  if (output.size() != nhap)
    output.resize(nhap, false);
  int h1 = 2 * curr_ind;
  int h2 = 2 * curr_ind + 1;
  int count = 0;

  vector<int>::iterator idx2 = clusters[clustering[h1]].begin();
  while (idx2 != clusters[clustering[h1]].end()) {
    if (!output[*idx2]) {
      output[*idx2] = true;
      count++;
    }
    idx2++;
  }

  idx2 = clusters[clustering[h2]].begin();
  while (idx2 != clusters[clustering[h2]].end()) {
    if (!output[*idx2]) {
      output[*idx2] = true;
      count++;
    }
    idx2++;
  }

  output[2 * curr_ind] = false;
  output[2 * curr_ind + 1] = false;
  count -= 2;
  assert(count >= ncond);
  return (0);
}

// RANDOM MASKING - EXPERIMENTS ONLY
/*
  int haplotypeCluster::mask(int curr_ind,vector<bool> & output) {
  assert(curr_ind<nhap/2);
  if(output.size() != nhap) output.resize(nhap,false);
  set<int> idx1 = sampleIndex(c1*ncond,nhap);
  for(set<int>::iterator idx2=idx1.begin();idx2!=idx1.end();idx2++) {
  output[*idx2] = true;
  }
  //assert(count>0);
  output[2*curr_ind] = false;
  output[2*curr_ind+1] = false;

  return(0);
  }
 */

int haplotypeCluster::topupClusters() {
  if (_DEBUG > 1)
    cout << "Topping up clusters." << endl;
  for (map<unsigned int, vector<int> >::iterator it = clusters.begin();
       it != clusters.end(); it++) {
    if (_DEBUG > 5)
      cout << it->first << " " << it->second.size() << endl;
    if (it->second.size() < min_cluster_size) {
      if (_DEBUG > 5)
        cout << "Topping up " << it->first << endl;
      int k_needed = min_cluster_size - it->second.size();
      while (k_needed > 0) {
        int addme = putils::getRandom(nhap);
        if (clustering[addme] != it->first) {
          it->second.push_back(addme);
          k_needed--;
          if (_DEBUG > 5)
            cout << addme << endl;
        }
      }
    }
  }

  return (0);
}

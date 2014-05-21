#include <objects/haplotypeCluster.h>

pair<int, int> argmin2(vector<float> &input) {

  pair<int, int> ret;
  if (input[0] < input[1]) {
    ret.first = 0;
    ret.second = 1;
  } else {
    ret.first = 1;
    ret.second = 0;
  }

  for (int i = 2; i < input.size(); i++) {
    if (input[i] < input[ret.second]) {
      if (input[i] < input[ret.first]) {
        ret.second = ret.first;
        ret.first = i;
      } else {
        ret.second = i;
      }
    }
  }
  return (ret);
}

// entry call to bifurc_kmeans
int haplotypeCluster::bifurc_kmeans() {
  clustering.assign(nsnp, 0);
  pair<int, int> n = kmeans(0, 0);
  assert(n.first > 0 && n.second > 0);
  if (_DEBUG > 0)
    cout << "Initial bifurc " << n.first << " " << n.second << endl;
  bifurc_kmeans(0, n.first, 1);
  bifurc_kmeans(1, n.second, 1);

  int ncluster = 0;
  for (int i = 0; i < nhap; i++) {
    if (cluster_size.count(clustering[i]))
      cluster_size[clustering[i]]++;
    else
      cluster_size[clustering[i]] = 0;
  }
  for (map<unsigned int, int>::iterator it = cluster_size.begin();
       it != cluster_size.end(); it++) {
    clusters[it->first] = vector<int>();
    clusters[it->first].reserve(it->second);
  }

  for (int i = 0; i < nhap; i++) {
    clusters[clustering[i]].push_back(i);
  }

  return (ncluster);
}

int haplotypeCluster::bifurc_kmeans(unsigned int cluster, int n, int depth) {
  // cout << "bifurc received " << input.size() << endl;
  if (n <= max_cluster_size) {
    // cout <<"bifurc terminating"<<endl;
    return (0);
  } else {
    depth++;
    pair<int, int> n2 = kmeans(cluster, depth);
    assert(n2.first > 0);
    assert(n2.second > 0);
    if (_DEBUG > 0)
      cout << "\tbifurc " << n2.first << " " << n2.second << endl;
    bifurc_kmeans(cluster, n2.first, depth);
    bifurc_kmeans(cluster + pow(2, depth), n2.second, depth);
    return (0);
  }
}

int haplotypeCluster::hap2float(int idx, vector<float> &output) {
  for (int j = 0; j < (cstop - cstart); j++)
    output[j] = (float)(g->vecH[idx][cstart + j]);
  return (0);
}

void printmu(vector<vector<float> > &mu) {
  int nprint = 30;
  for (int i = 0; i < mu.size(); i++) {
    cout << i << ": ";
    for (int j = 0; j < nprint; j++)
      cout << mu[i][j] << " ";
    cout << endl;
  }
}

int argmin(vector<float> &input) {
  int minind = 0;
  for (int i = 1; i < input.size(); i++)
    if (input[i] < input[minind])
      minind = i;
  return (minind);
}

int argmin(vector<int> &input) {
  int minind = 0;
  for (int i = 1; i < input.size(); i++)
    if (input[i] < input[minind])
      minind = i;
  return (minind);
}

int haplotypeCluster::kplusplus(int K, vector<int> &input,
                                vector<int> &output) {
  int ncand = 100;
  vector<float> d;
  float maxdist, mindist;
  set<int>::iterator idx1;
  int cand;

  for (int i = 0; i < K; i++) {
    set<int> candidates = sampleIndex(ncand, input.size());
    if (i == 0) {
      cand = input[*candidates.begin()];
      output[i] = cand;
    } else {
      output[i] = -1;
      maxdist = -1.0;
      d.assign(i, 0.0);
      for (idx1 = candidates.begin(); idx1 != candidates.end(); idx1++) {
        cand = input[*idx1];
        for (int j = 0; j < i; j++)
          d[j] = g->vecH[output[j]].hamming(g->vecH[cand], cstart, cstop);

        mindist = *min_element(d.begin(), d.begin() + i);
        if (mindist > maxdist) {
          output[i] = cand;
          maxdist = mindist;
        }
      }
    }
  }

  return (0);
}

int haplotypeCluster::kplusplus(int K, vector<int> &output) {
  int ncand = 100;
  vector<float> d;
  float maxdist, mindist;
  set<int>::iterator idx1;

  for (int i = 0; i < K; i++) {
    set<int> candidates = sampleIndex(ncand, nhap);
    if (i == 0) {
      output[i] = *candidates.begin();
    } else {
      output[i] = -1;
      maxdist = -1.0;
      d.assign(i, 0.0);
      for (idx1 = candidates.begin(); idx1 != candidates.end(); idx1++) {
        for (int j = 0; j < i; j++)
          d[j] = g->vecH[output[j]].hamming(g->vecH[*idx1], cstart, cstop);

        mindist = *min_element(d.begin(), d.begin() + i);
        if (mindist > maxdist) {
          output[i] = *idx1;
          maxdist = mindist;
        }
      }
    }
  }

  return (0);
}

int haplotypeCluster::updateMeans(pair<unsigned int, unsigned int> &cluster_id,
                                  pair<int, int> &cluster_size) {

  if (_DEBUG > 3) {
    printmu(mu);
    cout << endl;
  }

  mu.assign(K, vector<float>(nsnp, 0.0));
  cluster_size.first = 0;
  cluster_size.second = 0;

  for (int i = 0; i < nhap; i++) {
    if (cluster_mask[i]) {
      int k = clustering[i] != cluster_id.first;
      for (int j = 0; j < nsnp; j++)
        mu[k][j] += (float)(g->vecH[i][cstart + j]);
      if (clustering[i] == cluster_id.first)
        cluster_size.first++;
      else
        cluster_size.second++;
    }
  }

  if (_DEBUG > 1)
    cout << cluster_size.first << "\t" << cluster_size.second << endl;

  for (int j = 0; j < nsnp; j++) {
    mu[0][j] /= cluster_size.first;
    mu[1][j] /= cluster_size.second;
  }

  if (_DEBUG > 3) {
    printmu(mu);
    cout << endl;
  }

  return (0);
}

int haplotypeCluster::summariseKmeans() {
  int minc = nhap;
  int maxc = 0;
  int totc = 0;
  unsigned long ncalc = 0;
  int l;
  for (map<unsigned int, vector<int> >::iterator canopy = clusters.begin();
       canopy != clusters.end(); canopy++) {
    l = canopy->second.size();
    if (l < minc)
      minc = l;
    if (l > maxc)
      maxc = l;
    totc += l;
    ncalc += (l * (l - 1)) / 2;
  }

  cout << clusters.size() << " clusters." << endl;
  ncalc = 0;
  cout << "Max: " << maxc << "\tMin: " << minc
       << "\tMean: " << totc / clusters.size() << endl;
  return (0);
}

float haplotypeCluster::euc(int i, int j, vector<vector<float> > &mu,
                            vector<vector<vector<float> > > &dlook) {
  float ret = 0.0;
  for (int k = 0; k < nsnp / 8; k++) {
    unsigned char rawval = g->vecH[i].data[cstart / 8 + k];
    if (dlook[j][k][rawval] < 0.0) {
      bitset<8> bs(rawval);
      dlook[j][k][rawval] = 0.0;
      for (int l = 0; l < 8; l++)
        dlook[j][k][rawval] +=
            pow((double)(mu[j][k * 8 + l] - (float)bs[l]), 2.0);
    }
    ret += dlook[j][k][rawval];
  }

  return (ret);
}

float haplotypeCluster::euc(int idx, vector<float> &mu) {
  float ret = 0.0;
  for (int j = 0; j < (cstop - cstart); j++) {
    ret += pow((double)(mu[j] - (float)(g->vecH[idx][cstart + j])), 2.0);
  }
  //	cout << ret  <<" ";
  return (ret);
}

int haplotypeCluster::getIdx(int idx) {
  int count = 0;
  int i = 0;
  while (count <= idx) {
    count += cluster_mask[i];
    i++;
  }
  i--;
  assert(i < nhap);
  assert(cluster_mask[i]);
  return (i);
}

pair<int, int> haplotypeCluster::kmeans(unsigned int cluster,
                                        unsigned int depth) {
  int nchanged = 0;
  int debug = 0;

  if (_DEBUG > 1)
    cout << "Clustering (K-means) K = " << K << ". NSNP = " << nsnp << endl;

  pair<unsigned int, unsigned int> cluster_id(cluster, cluster + pow(2, depth));

  int N = 0;
  for (int i = 0; i < nhap; i++) {
    if (clustering[i] == cluster) {
      cluster_mask[i] = true;
      N++;
    } else
      cluster_mask[i] = false;
  }

  mu.assign(K, vector<float>(nsnp, -1));
  dlook.resize(K);
  pair<int, int> cluster_size;
  vector<float> d(K);
  vector<float> tmp(nsnp);
  float ss;
  //	kplusplus(K,tmp2,initial_centroids);
  int idx1 = getIdx(putils::getRandom(N));

  int idx2 = getIdx(putils::getRandom(N));
  float maxdist = g->vecH[idx1].hamming(g->vecH[idx2], cstart, cstop);

  for (int i = 0; i < 100; i++) {
    int tmp1 = getIdx(putils::getRandom(N));
    float tmp2 = g->vecH[idx1].hamming(g->vecH[tmp1], cstart, cstop);
    if (tmp2 > maxdist) {
      idx2 = tmp1;
      maxdist = tmp2;
    }
  }

  for (int j = 0; j < nsnp; j++) {
    mu[0][j] = (float)(g->vecH[idx1][cstart + j]);
    mu[1][j] = (float)(g->vecH[idx2][cstart + j]);
  }

  if (_DEBUG > 1) {
    printmu(mu);
    cout << endl;
  }
  int closest_cluster;

  // LLOYDS ALGORITHM - standard K-means routine
  for (int iteration = 0; iteration < niteration; iteration++) {
    nchanged = 0;
    for (int i = 0; i < K; i++)
      dlook[i].assign(nsnp / 8, vector<float>(256, -1.0));

    ss = 0.0;
    for (int i = 0; i < nhap; i++) {
      if (cluster_mask[i]) {
        for (int j = 0; j < K; j++)
          d[j] = euc(i, j, mu, dlook);
        if (d[0] < d[1]) {
          ss += d[0];
          closest_cluster = cluster_id.first;
        } else {
          ss += d[1];
          closest_cluster = cluster_id.second;
        }
        if (clustering[i] != closest_cluster)
          nchanged++;

        clustering[i] = closest_cluster;
      }
    }

    updateMeans(cluster_id, cluster_size);
    if (_DEBUG > 1)
      cout << "LLOYDS K-MEANS ITERATION " << iteration
           << " Mean SS = " << ss / (float)N << "\t" << nchanged
           << " haps changed clusters." << endl;
    if (nchanged == 0)
      break;
  }

  //	if(_DEBUG>0) cout << "Lloyds Total SS = " << SS() << "\t"<<nchanged << "
  //haps changed clusters on last iteration" << endl;

  int minind = cluster_size.first > cluster_size.second;
  int minsize = min(cluster_size.first, cluster_size.second);

  if (minsize < ncond && K == 2) { // not partioning well. do a random split.
    vector<pair<int, float> > closest;
    if (_DEBUG > 0)
      cout << "WARNING: " << N << " did not partition well. "
           << cluster_size.first << " " << cluster_size.second
           << ". Regrouping." << endl;
    closestN(N / 2, minind, closest);
    for (int j = 0; j < closest.size();
         j++) { // change cluster of closet N2 guys.
      if (minind == 0) {
        clustering[closest[j].first] = cluster_id.first;
        cluster_size.first++;
        cluster_size.second--;
      } else {
        clustering[closest[j].first] = cluster_id.second;
        cluster_size.first--;
        cluster_size.second++;
      }
    }
  }

  return (cluster_size);
};

int haplotypeCluster::closestN(int N, int clusteridx,
                               vector<pair<int, float> > &output) {

  output.resize(N);
  int count = 0;
  float d, maxdist = -1.0;
  int maxind = -1;

  for (int i = 0; i < nhap; i++) {
    if (cluster_mask[i]) {
      d = euc(i, clusteridx, mu, dlook);
      pair<int, float> obs(i, d);
      if (count < N) {
        output[count] = obs;
        if (obs.second > maxdist) {
          maxdist = obs.second;
          maxind = count;
        }
      } else if (obs.second > maxdist) {
        output[maxind] = obs;
        maxind = 0;
        for (int j = 1; j < N; j++)
          if (output[j].second > output[maxind].second)
            maxind = j;
        maxdist = output[maxind].second;
      }
      count++;
    }
  }
  return 0;
}

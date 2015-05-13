#include "glVQ.hpp"

//#define DEBUG 1

#ifdef DEBUG
#include <iostream>
#endif

using namespace std;

// returns pair of normalized to 1 GLs for Ref/Ref and Het
pair<float, float> GLVQHelper::GLTrio2RRHet(vector<float>::const_iterator gl) {

  float sum = *gl + *(gl + 1) + *(gl + 2);
  assert(sum > 0);
  return make_pair(max(*gl / sum, 0.0f), max(*(gl + 1) / sum, 0.0f));
}

// create code book a la wikipedia
// https://en.wikipedia.org/wiki/Vector_quantization
GLVQ::GLVQ(const vector<float> &inGLs, gsl_rng &rng, unsigned codeBookSize)
    : m_rng(rng) {
  BuildCodeBook(codeBookSize, inGLs);
}

void GLVQ::BuildCodeBook(unsigned codeBookSize, const vector<float> &inGLs) {

  assert(inGLs.size() % 3 == 0);
  const unsigned numGLs = inGLs.size() / 3;
  assert(numGLs > 0);

  assert(codeBookSize >= 4);
  if (codeBookSize - 1 > numeric_limits<unsigned char>::max())
    throw std::runtime_error(
        "codebook size is too large: " + to_string(codeBookSize) + " > " +
        to_string(std::numeric_limits<unsigned char>::max() + 1));
  m_codeBook.reserve(codeBookSize);

  // initialize code book centroids randomly, but add in four positions that I
  // think are important
  m_codeBook.push_back(make_pair(0, 1));
  m_codeBook.push_back(make_pair(0, 0));
  m_codeBook.push_back(make_pair(1, 0));
  m_codeBook.push_back(make_pair(1.0 / 3, 1.0 / 3));

  // initialize randomly in 0,1 space
  for (unsigned entryNum = 4; entryNum < codeBookSize; ++entryNum) {
    const float RR = gsl_rng_uniform(&m_rng);
    const float Het = gsl_rng_uniform(&m_rng) * (1 - RR);
    m_codeBook.push_back(make_pair(RR, Het));
  }

  // create random set of points to train on
  vector<tuple<float, float, unsigned char>> trainingSet;

  // training set size seems to make a difference
  // let's try all gls
  const size_t trainingSetSize = numGLs;
  trainingSet.reserve(trainingSetSize);
  for (unsigned iter = 0; iter < trainingSetSize; ++iter) {
    unsigned glIdx = 0;
    if (trainingSetSize < numGLs)
      glIdx = gsl_rng_uniform_int(&m_rng, numGLs);
    else
      glIdx = iter;

    float RR = inGLs[glIdx * 3];
    float Het = inGLs[glIdx * 3 + 1];
    float sum = inGLs[glIdx * 3 + 2] + RR + Het;
    RR /= sum;
    Het /= sum;
    trainingSet.push_back(make_tuple(RR, Het, 0));
  }

  // train codebook
  double previousVar = numeric_limits<double>::max();
  double currentVar = AssignPoints(trainingSet);
  size_t iter = 0;
  const size_t maxIter = 100;
  while (abs(currentVar - previousVar) > currentVar / 100) {
    if (++iter > maxIter)
      throw std::runtime_error("[GLVQ] Clustering failed to converge within " +
                               to_string(maxIter) + " iterations");
    previousVar = currentVar;
    m_codeBook = UpdateCodeBook(trainingSet);
    currentVar = AssignPoints(trainingSet);
  }
}

// assign the points in the training set to the nearest centroid of the codeBook
double
GLVQ::AssignPoints(vector<tuple<float, float, unsigned char>> &points) const {

  // find nearest centroid
  double var = 0;
  for (auto &p : points) {
    GLVQHelper::EuclidianDist comp(get<0>(p), get<1>(p));
    get<2>(p) =
        distance(m_codeBook.cbegin(),
                 min_element(m_codeBook.cbegin(), m_codeBook.cend(), comp));
    var += pow(get<0>(p) - m_codeBook[get<2>(p)].first, 2) +
           pow(get<1>(p) - m_codeBook[get<2>(p)].second, 2);
  }
  return var;
}

vector<pair<float, float>>
GLVQ::UpdateCodeBook(const vector<tuple<float, float, unsigned char>> &points) {

  // new codebook
  vector<pair<float, float>> newCodeBook;
  newCodeBook.reserve(m_codeBook.size());

  // keeps track of how many points each cluster has
  vector<size_t> numPoints;
  numPoints.reserve(m_codeBook.size());

  // initialize clusters to close to 0
  for (size_t i = m_codeBook.size(); i != 0; --i) {
    newCodeBook.push_back(make_pair(0, 0));
    numPoints.push_back(0);
  }

  // sum up points
  for (auto p : points) {
    const size_t entryNum = get<2>(p);
    newCodeBook[entryNum].first += get<0>(p);
    newCodeBook[entryNum].second += get<1>(p);
    ++numPoints[entryNum];
  }

  // average results
  for (size_t i = 0; i < numPoints.size(); ++i)
    if (numPoints[i] > 0) {
      newCodeBook[i].first /= numPoints[i];
      newCodeBook[i].second /= numPoints[i];
    }
    // if no point is closest to this centroid, then assign it to a random point
    else {
      const auto &point = points[gsl_rng_uniform_int(&m_rng, points.size())];
      newCodeBook[i].first = get<0>(point);
      newCodeBook[i].second = get<1>(point);
    }

  // make sure the first four points are kept at important points
  assert(m_codeBook.size() > 3);
  /*  newCodeBook[0] = make_pair(0, 1);
    newCodeBook[1] = make_pair(0, 0);
    newCodeBook[2] = make_pair(1, 0);
    newCodeBook[3] = make_pair(1.0 / 3, 1.0 / 3);
  */
  return newCodeBook;
}

// use code book to find code corresponding to a pair of floats
unsigned char GLVQ::FindGLCode(float RR, float Het) const {

  // the gl code is the distance between the codebook element with the least
  // distance between it and the input gls, and the first element
  // yay for function objects!
  assert(m_codeBook.size() > 0);
  GLVQHelper::EuclidianDist comp(RR, Het);
  return distance(m_codeBook.cbegin(),
                  min_element(m_codeBook.cbegin(), m_codeBook.cend(), comp));
}

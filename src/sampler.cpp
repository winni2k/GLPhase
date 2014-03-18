#include "sampler.hpp"

unsigned Sampler::Sampler(gsl_rng *rng, unsigned numSamples,
                          unsigned numHaplotypes) {

  // make sure we have at least as many haplotypes as samples
  assert(numSamples * 2 <= numHaplotypes);
}

unsigned Sampler::SampleHap(unsigned excludeInd, bool onlyFromRef) {

  if (onlyFromRef) {
    unsigned retHap{ 0 };

    // so I'm backing away from that for now
    while (1) {
      retHap = SampleHap(excludeInd);
      if (retHap >= m_numSamples * 2)
        break;
    }
    return (retHap);
  } else {
    return SampleHap(excludeInd);
  }
}

unsigned UnifSampler::SampleHap(unsigned excludeInd) {

  assert(excludeInd < m_numSamples);
  while (1) {

    // m_uCols is 1 based, but gsl_rng makes 0 based choice
    unsigned propHap = gsl_rng_uniform_int(m_rng, m_numHaps);
    if (propHap / 2 != excludeInd)
      return propHap;
  }
}

void GraphSampler::GraphSampler(gsl_rng *rng, unsigned numSamples,
                                unsigned numHaps, GraphSampT graphType) {

  // figure out how many rows and columns our matrix is going to have
  switch (graphType) {
  case GraphSampT::sampSamp:
    m_usingHaps = false;
    break;
  case GraphSampT::sampHap:
    m_usingHaps = true;
    break;
  default:
    throw "unexpected graph type encountered";
  }

  if (m_usingHaps)
    m_colHapFactor = 2;
  else
    m_colHapFactor = 1;

  m_cols =
      ceil(static_cast<float>(m_numHaps) / static_cast<float>(m_colHapFactor));

  // resize the relationship matrixes to the appropriate size
  // assign all numerators and denominators a count of 1
  m_2dRelationshipMatNum.resize(m_rows);
  m_2dRelationshipMatDen.resize(m_rows);
  for (unsigned i = 0; i < m_uRows; i++) {
    m_2dRelationshipMatNum[i].resize(m_cols, 1);
    m_2dRelationshipMatDen[i].resize(m_cols, 1);
  }
}

unsigned GraphSampler::SampleHap(unsigned excludeInd) {

  // make sure input is sensible
  assert(excludeInd < m_numSamples);

  vector<float> &vuRelRowNum = m_2dRelationshipMatNum[uInd];
  vector<float> &vuRelRowDen = m_2dRelationshipMatDen[uInd];

  // initialize propHap to out of range value
  unsigned propHap = m_numHaps;
  while (1) {

    // sample a haplotype
    propHap = gsl_rng_get(m_rng) % m_numHaps;

    // resample if individual chosen is the same as
    // the individual not to sample from
    if (propHap / 2 == excludeInd)
      continue;

    // get hap corresponding column index
    unsigned propCol = Hap2Col(propHap);

    assert(vuRelRowNum[propCol] / vuRelRowDen[propCol] > 0);
    assert(vuRelRowNum[propCol] / vuRelRowDen[propCol] <= 1);

    // resample if individual does not pass rejection sample
    if (gsl_rng_uniform(rng) <= vuRelRowNum[propCol] / vuRelRowDen[propCol])
      break;
  }

  // make sure the return value is sensible
  assert(propHap < m_numHaps);

  return propHap;
}

void GraphSampler::UpdatePropDistProp(const vector<unsigned> &propHaps,
                                      unsigned updateIndNum, bool accepted,
                                      float penRatio = 1) {

  // update relationship matrix
  for (auto &propHap : propHaps) {

    // convert to correct index based on cols representing haps vs samples
    unsigned propCol = Hap2Col(propHap);

    // update proposal hap or sample for current sample
    m_2dRelationshipMatDen[updateIndNum][propCol] += penRatio;
    if (bAccepted)
      m_2dRelationshipMatNum[updateIndNum][propCol] += penRatio;

    // update current sample for proposal hap or sample if proposal is not
    // from reference panel
    if (m_2dRelationshipMatDen.size() > propCol) {
      if (m_bUsingHaps) {
        m_2dRelationshipMatDen[propCol][Col2Hap(updateIndNum)] += penRatio;
        m_2dRelationshipMatDen[propCol][Col2Hap(updateIndNum) + 1] += penRatio;
        if (bAccepted) {
          m_2dRelationshipMatNum[propCol][Col2Hap(updateIndNum)] += penRatio;
          m_2dRelationshipMatNum[propCol][Col2Hap(updateIndNum) + 1] +=
              penRatio;
        }
      } else {
        m_2dRelationshipMatDen[propCol][updateIndNum] += penRatio;
        if (bAccepted)
          m_2dRelationshipMatNum[propCol][updateIndNum] += penRatio;
      }
    }
  }
}

// based on code from SNPTools::Impute
// name should be a list of sample names that includes reference panel sample
// names
// fileName should be the basename for output
void GraphSampler::Save(string fileName, const vector<string> &name) {

  cerr << "Saving Relationship graph to prefix " << fileName << " ..." << endl;

  {
    unsigned expectedNumNames = m_2dRelationshipMatDen[0].size();
    if (m_usingHaps)
      expectedNumNames = expectedNumNames / 2 + expectedNumNames % 2;

    assert(name.size() == expectedNumNames);
  }

  string numeratorFile = fileName + ".relGraph.num.gz";
  string denominatorFile = fileName + ".relGraph.den.gz";
  string ratioFile = fileName + ".relGraph.gz";

  typedef std::shared_ptr<ofile> ofile_ptr;
  ofile_ptr numFile(new ofile(numeratorFile));
  ofile_ptr denFile(new ofile(denominatorFile));
  ofile_ptr ratFile(new ofile(ratioFile));

  vector<ofile_ptr> ofiles;
  ofiles.push_back(numFile);
  ofiles.push_back(denFile);
  ofiles.push_back(ratFile);

  // start printing header
  *numFile << "Numerator";
  *denFile << "Denominator";
  *ratFile << "Numerator/Denominator";

  for (unsigned uFileNum = 0; uFileNum < ofiles.size(); uFileNum++) {

    // finish printing header
    for (auto sName : name) {
      *ofiles[uFileNum] << "\t" << sName;
      if (m_bUsingHaps)
        *ofiles[uFileNum] << ".hapA"
                          << "\t" << sName << ".hapB";
    }
    *ofiles[uFileNum] << endl;

    // print data rows
    // cycle through samples
    for (unsigned uRowNum = 0; uRowNum < m_uRows; uRowNum++) {

      // print sample name
      *ofiles[uFileNum] << name[uRowNum];

      // print rest of row
      for (unsigned uColNum = 0; uColNum < m_uCols; uColNum++) {
        float fPrintVal = 0;
        switch (uFileNum) {
        case 0:
          fPrintVal = m_2dRelationshipMatNum[uRowNum][uColNum];
          assert(fPrintVal != 0);
          break;
        case 1:
          fPrintVal = m_2dRelationshipMatDen[uRowNum][uColNum];
          assert(fPrintVal != 0);
          break;
        case 2:
          fPrintVal = m_2dRelationshipMatNum[uRowNum][uColNum] /
                      m_2dRelationshipMatDen[uRowNum][uColNum];
          assert(fPrintVal != 0);
          break;
        default:
          cerr << "Programming error: RelationshipGraph::Save: Unknown file "
                  "number: " << uFileNum << endl;
          assert(false);
        }
        *ofiles[uFileNum] << "\t" << fPrintVal;
      }
      *ofiles[uFileNum] << endl;
    }
  }
  cerr << "\t...saving of relationship graph complete." << endl;
}

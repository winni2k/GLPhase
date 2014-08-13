/* @(#)bcfReader.cpp
 */

#include "bcfReader.hpp"

using namespace std;

namespace BCFReaderHelper {
double phred2Prob(double phred) {
  assert(phred >= 0);

  if (phred >= DBL_MAX_10_EXP * 10)
    return DBL_MIN;
  else
    return std::pow(10, -phred / 10);
}
}

BCFReader::BCFReader(string fileName, BCFReaderHelper::extract_t extractType,
                     string region)
    : m_extractType(extractType) {

  assert(fileName.size() > 0);
  bcf_srs_t *readers = bcf_sr_init();

  // Set the region in the bcf to iterate over
  if (!region.empty())
    if (bcf_sr_set_regions(readers, region.c_str(), 0) < 0)
      throw std::runtime_error("[BCFReader] Failed "
                               "to read the "
                               "regions: " +
                               region);

  // read header and check it
  if (!bcf_sr_add_reader(readers, fileName.c_str()))
    throw std::runtime_error("[BCFReader] Could not open file: " + fileName);

  bcf_hdr_t *hdr = readers->readers[0].header;

  // parse #CHROM header line for samples
  if (bcf_hdr_nsamples(hdr) < 1)
    throw std::runtime_error("[BCFReader] No samples in input file: " +
                             fileName);
  assert(m_sampNames.empty());
  m_sampNames.reserve(bcf_hdr_nsamples(hdr));
  for (int32_t sampNum = 0; sampNum < bcf_hdr_nsamples(hdr); ++sampNum)
    m_sampNames.push_back(hdr->samples[sampNum]);

  // extract each line of data
  unsigned lineNum = 0;
  string extractString;
  while (bcf_sr_next_line(readers)) {

    bcf1_t *rec = readers->readers[0].buffer[0];

    // store site information
    string chr = bcf_hdr_id2name(hdr, rec->rid);
    int pos = rec->pos + 1;
    string a1(rec->d.allele[0]);
    string a2(rec->d.allele[1]);

    // define GL type to search for
    if (lineNum == 0) {
      if (m_extractType == BCFReaderHelper::extract_t::Haps)
        extractString = "GT";
      else if (m_extractType == BCFReaderHelper::extract_t::GL) {
        if (bcf_get_fmt(hdr, rec, "GL"))
          extractString = "GL";
        else if (bcf_get_fmt(hdr, rec, "PL"))
          extractString = "PL";
        else
          throw std::runtime_error("Could not find GL or PL field in VCF/BCF");
      } else
        throw std::logic_error("unexpected BCFReaderHelper::extract_t type");
    }

    // now check if expected field exists
    if (!bcf_get_fmt(hdr, rec, extractString.c_str()))
      throw std::runtime_error("expected " + extractString +
                               " field in VCF/BCF");

    // read haps
    try {
      if (m_extractType == BCFReaderHelper::extract_t::Haps)
        m_haps.push_back(ExtractRecAlleles(rec, hdr));

      // read GLs
      else if (m_extractType == BCFReaderHelper::extract_t::GL)
        m_GLs.push_back(ExtractRecGLs(rec, hdr, extractString));

      // could not figure out what to read
      else
        throw std::logic_error("unexpected BCFReaderHelper::extract_t type");
    }
    catch (std::runtime_error &e) {
      throw std::runtime_error(string(e.what()) + " at line " +
                               to_string(lineNum) + " and site " + chr + ":" +
                               to_string(pos));
    }

    m_sites.push_back(std::move(Bio::snp(chr, pos, a1, a2)));
    ++lineNum;
  }

  // clean up
  bcf_sr_destroy(readers);
}

vector<char> BCFReader::ExtractRecAlleles(bcf1_t *rec, bcf_hdr_t *hdr) {

  assert(m_extractType == BCFReaderHelper::extract_t::Haps);
  int *gt_arr = nullptr, ngt_arr = 0;
  const int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
  if (ngt != 2 * static_cast<int>(m_sampNames.size())) {
    free(gt_arr);
    throw std::runtime_error("Malformed VCF. Too few or too many "
                             "GT fields");
  }

  vector<char> siteAlleles;
  siteAlleles.reserve(ngt);
  for (int gtNum = 0; gtNum < ngt; ++gtNum) {
    assert(gt_arr[gtNum] != bcf_gt_missing);
    assert(gt_arr[gtNum] != bcf_int32_vector_end);

    if ((gtNum & 1) == 1 && !bcf_gt_is_phased(gt_arr[gtNum])) {
      free(gt_arr);
      throw std::runtime_error("Error in GT data, genotype is not phased.");
    }

    assert(bcf_gt_allele(gt_arr[gtNum]) < 2);
    siteAlleles.push_back(bcf_gt_allele(gt_arr[gtNum]) == 1);
  }

  // store site alleles
  free(gt_arr);
  return siteAlleles;
}

vector<double> BCFReader::ExtractRecGLs(bcf1_t *rec, bcf_hdr_t *hdr,
                                        const string &extractString) {

  assert(m_extractType == BCFReaderHelper::extract_t::GL);
  int m_arr = 0;
  float *arr = NULL;
  int stride = 3;
  int n_arr =
      bcf_get_format_float(hdr, rec, extractString.c_str(), &arr, &m_arr);
  if (n_arr / stride != bcf_hdr_nsamples(hdr)) {
    free(arr);
    throw std::runtime_error("Malformed VCF. Too few or too many "
                             "GT fields");
  }

  // convert GL to double
  vector<double> siteGLs;
  siteGLs.reserve(n_arr);
  if (extractString == "GL")
    for (int glNum = 0; glNum != n_arr; ++glNum)
      siteGLs.push_back(pow(10.0f, arr[glNum]));
  else if (extractString == "PL")
    for (int glNum = 0; glNum != n_arr; ++glNum)
      siteGLs.push_back(BCFReaderHelper::phred2Prob(arr[glNum]));

  // store gls
  free(arr);
  return siteGLs;
}

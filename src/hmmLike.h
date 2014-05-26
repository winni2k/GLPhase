/* @(#)hmmLike.h
 */

// this file shall only contain gcc compatible namespace definitions

#ifndef _HMMLIKE_H
#define _HMMLIKE_H 1

#include "globals.h"
#include <stdint.h>
#include <utility>
#include <thrust/host_vector.h>
#include <exception>
#include <stdexcept>
#include <sstream>
#include <string>

// this is the number of haplotypes we use in the HMM
// That is, the HMM has NUMSITES * NUMHAPS states
#define NUMHAPS 4

// these functions are implemented in HMMLike.cu
namespace HMMLikeCUDA {
extern void CheckDevice();
extern void CopyTranToDevice(const std::vector<float> &tran);
extern void CopyMutationMatToDevice(const float (*mutMat)[4][4]);
extern void CopyPackedGLsToDevice(const std::vector<uint32_t> &packedGLs);
extern void
CopyCodeBookToDevice(const std::vector<std::pair<float, float> > &codeBook);
extern void CopyHapPanelToDevice(const std::vector<uint64_t> &hapPanel);
extern void SetUpRNGs(size_t numSamples, unsigned long seed);
extern void RunHMMOnDevice(const std::vector<uint64_t> &hapPanel,
                           const std::vector<unsigned> &extraPropHaps,
                           unsigned numSites, unsigned numSamples,
                           unsigned numCycles, std::vector<unsigned> &hapIdxs);

extern void Cleanup();
}

#endif /* _HMMLIKE_H */

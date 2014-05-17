/* @(#)hmmLike.h
 */

// this file shall only contain gcc compatible namespace definitions

#ifndef _HMMLIKE_H
#define _HMMLIKE_H 1

// this is the number of sites that the HMM will run on
#define NUMSITES 1024
#define WN 16 // word num = NUMSITES / 64

#define WordShift 6
// input wordsize is 64 because we are working with uint64_t
#define WORDSIZE 64
#define WordMod 63

// this is the number of haplotypes we use in the HMM
// That is, the HMM has NUMSITES * NUMHAPS states
#define NUMHAPS 4


// these functions are implemented in HMMLike.cu
namespace HMMLikeCUDA {
extern void CheckDevice();
extern void CopyTranToDevice(const std::vector<float> &tran);
extern void CopyMutationMatToDevice(const float (*mutMat)[4][4]);
extern void CopyPackedGLsToDevice(const std::vector<char> &packedGLs);
extern void RunHMMOnDevice(const std::vector<uint64_t> &hapPanel,
                           const std::vector<unsigned> &extraPropHaps,
                           unsigned numSites, unsigned numSamples,
                           unsigned numCycles, std::vector<unsigned> &hapIdxs,
                           unsigned long seed);

extern void Cleanup();
}



#endif /* _HMMLIKE_H */

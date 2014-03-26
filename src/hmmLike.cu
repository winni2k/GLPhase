#include "hmmLike.hcu"

#include <iostream>

using namespace std;

namespace HMMLikeCUDA {

__constant__ float transitionMat[NUMSITES * 3];
__constant__ float mutationMat[4 * 4];

__device__ void UnpackGLs(char GLset, float *GLs) {

  GLs[0] = (((GLset >> 4) & 15) + 0.5f) / 16;
  GLs[1] = ((GLset & 15) + 0.5f) / 16;
  GLs[2] = max(1 - GLs[0] - GLs[1], 0.0f);
}

// definition of HMM kernel
__global__ void hmmLike(const char *__restrict__ d_packedGLs,
                        const uint64_t *__restrict__ d_hapPanel,
                        const unsigned *__restrict__ d_hapIdxs,
                        const unsigned *__restrict__ d_extraPropHaps,
                        unsigned *d_chosenHapIdxs, unsigned numSamples,
                        unsigned numCycles, unsigned idxOffset) {

  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  int sampIdx = idx + idxOffset;

  if (idx < numSamples) {

    // define emission matrix

    for (int cycle = 0; cycle < numCycles; ++cycle) {
      for (int site = 0; site < NUMSITES; ++site) {

        /*
          Convert packed GLs back to floats
        */
        float GLs[3];
        UnpackGLs(d_packedGLs[idx + site * NUMSITES], GLs);

        /*
          return test result based on GLs
         */
        d_chosenHapIdxs[idx * 4 + site] = idx * 4 + site;
        if (site == 4)
          return;

      } // end of site
    }   // end of cycle
  }
}

extern "C" void CheckDevice() {

  cudaError_t err;
  int deviceCount = 0;

  // note your project will need to link with cuda.lib files on windows
  cout << "Querrying CUDA Device(s)...\n\n";

  err = cudaGetDeviceCount(&deviceCount);

  if (err != cudaSuccess) {
    cout << "cuDeviceGetCount returned " << err << "\n";
    cout << "Result = FAIL\n";
    exit(EXIT_FAILURE); //
  }

  cout << "Found " << deviceCount << " device(s)" << endl;

  for (int i = 0; i < deviceCount; ++i) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    cout << "Device number: " << i << endl;
    cout << "Name: " << prop.name << endl;
    cout << "Compute capability: " << prop.major << '.' << prop.minor << endl;

    if (prop.major < 3 || prop.major == 3 && prop.minor < 5) {
      cerr << "Compute capability >= 3.5 required" << endl;
      exit(EXIT_FAILURE);
    }

    // set device to use more L1 cache than shared mem
    cout << "\nSetting device to prefer L1 cache over shared mem\n";
    if (cudaDeviceSetCacheConfig(cudaFuncCachePreferL1) != cudaSuccess) {
      cerr << "Could not set device caching preference to L1\n";
    }

    cudaFuncCache pCacheConfig;
    if (cudaDeviceGetCacheConfig(&pCacheConfig) != cudaSuccess) {
      cerr << "Could not get device caching preference\n";
    }
    if (pCacheConfig == cudaFuncCachePreferL1)
      cout << "Device caching preference is set to prefer L1" << endl;
    else
      cout << "Device caching preference is not set to prefer L1" << endl;
  }
}

extern "C" cudaError_t CopyTranToDevice(const vector<float> &tran) {

  assert(tran.size() == NUMSITES * 3);
  return cudaMemcpyToSymbol(transitionMat, tran.data(),
                            sizeof(float) * NUMSITES * 3, 0,
                            cudaMemcpyHostToDevice);
}

extern "C" cudaError_t CopyMutationMatToDevice(const float (*mutMat)[4][4]) {

  vector<float> h_mutMat(4 * 4);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      h_mutMat[i + 4 * j] = ((*mutMat))[i][j];

  return cudaMemcpyToSymbol(mutationMat, h_mutMat.data(), sizeof(float) * 4 * 4,
                            0, cudaMemcpyHostToDevice);
}

extern "C" void RunHMMOnDevice(const vector<char> &packedGLs,
                               const vector<uint64_t> &hapPanel,
                               const vector<unsigned> &extraPropHaps,
                               unsigned numSites, unsigned numSamples,
                               unsigned numCycles, unsigned sampIdxOffset,
                               vector<unsigned> &hapIdxs) {
  assert(numSites == NUMSITES);
  assert(numSites * numSamples == packedGLs.size());
  assert(hapIdxs.size() == numSamples * 4);
  assert(extraPropHaps.size() == numSamples * numCycles);

  /*
    copy packedGLs to device memory
  */
  cout << "Copying packed GLs to device\n";
  size_t glSize = packedGLs.size() * sizeof(char);

  // allocate memory on device
  char *d_packedGLs = NULL;
  cudaError_t err = cudaMalloc((void **)&d_packedGLs, glSize);
  if (err != cudaSuccess) {
    cerr << "Failed to allocate packed GLs on device\n";
    exit(EXIT_FAILURE);
  }

  // copy data across
  err =
      cudaMemcpy(d_packedGLs, packedGLs.data(), glSize, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    cerr << "Failed to copy packed GLs to device\n";
    exit(EXIT_FAILURE);
  }

  /*
    copy haplotypes to device memory
  */
  cout << "Copying HapPanel to device\n";
  size_t hapPanelSize = hapPanel.size() * sizeof(uint64_t);

  // allocate memory on device
  uint64_t *d_hapPanel = NULL;
  err = cudaMalloc((void **)&d_hapPanel, hapPanelSize);
  if (err != cudaSuccess) {
    cerr << "Failed to allocate hapPanel on device\n";
    exit(EXIT_FAILURE);
  }

  // copy data across
  err = cudaMemcpy(d_hapPanel, hapPanel.data(), hapPanelSize,
                   cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    cerr << "Failed to copy hapPanel to device\n";
    exit(EXIT_FAILURE);
  }

  /*
    copy initial hap indices to device memory
  */
  cout << "Copying hap indices to device\n";
  size_t hapIdxsSize = hapIdxs.size() * sizeof(unsigned);

  // allocate memory on device
  unsigned *d_hapIdxs = NULL;
  err = cudaMalloc((void **)&d_hapIdxs, hapIdxsSize);
  if (err != cudaSuccess) {
    cerr << "Failed to allocate hap indices on device\n";
    exit(EXIT_FAILURE);
  }

  // copy data across
  err = cudaMemcpy(d_hapIdxs, hapIdxs.data(), hapIdxsSize,
                   cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    cerr << "Failed to copy hap indices to device\n";
    exit(EXIT_FAILURE);
  }

  /*
    copy extra proposal haps to device memory
  */
  cout << "Copying extra proposal haps to device\n";
  size_t extraPropHapsSize = extraPropHaps.size() * sizeof(unsigned);

  // allocate memory on device
  unsigned *d_extraPropHaps = NULL;
  err = cudaMalloc((void **)&d_extraPropHaps, extraPropHapsSize);
  if (err != cudaSuccess) {
    cerr << "Failed to allocate extra prop haps on device\n";
    exit(EXIT_FAILURE);
  }

  // copy data across
  err = cudaMemcpy(d_extraPropHaps, extraPropHaps.data(), extraPropHapsSize,
                   cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    cerr << "Failed to copy extra prop haps to device\n";
    exit(EXIT_FAILURE);
  }

  /*
    allocate device memory for results
  */
  unsigned *d_chosenHapIdxs = NULL;
  err = cudaMalloc((void **)&d_chosenHapIdxs, hapIdxsSize);
  if (err != cudaSuccess) {
    cerr << "Failed to allocate memory for result hap idxs on device\n";
    exit(EXIT_FAILURE);
  }

  // determine thread and block size
  int threadsPerBlock = 32;
  int blocksPerRun = (numSamples + threadsPerBlock - 1) / threadsPerBlock;
  cout << "Running with " << threadsPerBlock << " threads per block in "
       << blocksPerRun << " thread blocks\n";

  /*
    run kernel
  */
  hmmLike << <blocksPerRun, threadsPerBlock>>>
      (d_packedGLs, d_hapPanel, d_hapIdxs, d_extraPropHaps, d_chosenHapIdxs,
       numSamples, numCycles, sampIdxOffset);

  /*
    copy result hap indices back to host into hapIdxs
  */
  err = cudaMemcpy(hapIdxs.data(), d_chosenHapIdxs, hapIdxsSize,
                   cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    cerr << "Failed to copy chosen indices to host\n";
    exit(EXIT_FAILURE);
  }

  // return nothing as the return data is stored in hapIdxs
  return;
};
}

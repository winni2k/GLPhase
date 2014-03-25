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
__global__ void hmmLike(const char *GLs, const unsigned *outSet, int idxOffset,
                        int numElements, int numCycles, int numSites, float mu,
                        float rho) {

  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  int sampIdx = idx + idxOffset;

  if (idx < numElements) {

    // define emission matrix

    for (int cycle = 0; cycle < numCycles; ++cycle) {
      for (int site = 0; site < numSites; ++site) {

        /*
          Convert packed GLs back to floats
        */
        float GLs[3];
        for (int i = 0; i < 3; ++i)
          GLs[i] = 0;
        UnpackGLs(GLs[idx], GLs);

      } // end of site
    }   // end of cycle
  }
}

void CheckDevice() {

  cudaError_t err;
  int deviceCount = 0;

  // note your project will need to link with cuda.lib files on windows
  cout << "Querrying CUDA Device(s)...\n\n";

  err = cudaGetDeviceCount(&deviceCount);

  if (err != cudaSuccess) {
    cout << "cuDeviceGetCount returned " << err << "\n";
    cout << "Result = FAIL\n";
    exit(EXIT_FAILURE);
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
  }
}

cudaError_t CopyTranToDevice(const vector<float> &tran) {

  assert(tran.size() == NUMSITES * 3);
  return cudaMemcpyToSymbol(transitionMat, tran.data(),
                            sizeof(float) * NUMSITES * 3, 0,
                            cudaMemcpyHostToDevice);
}
cudaError_t CopyMutationMatToDevice(const float (*mutMat)[4][4]) {

  vector<float> h_mutMat(4 * 4);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
        h_mutMat[i + 4 * j] = ((*mutMat))[i][j];

  return cudaMemcpyToSymbol(mutationMat, h_mutMat.data(), sizeof(float) * 4 * 4,
                            0, cudaMemcpyHostToDevice);
}
}

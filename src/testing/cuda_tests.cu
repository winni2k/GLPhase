//#include "cuda.h"
//#include "cuda_runtime.h"
#include "hmmLike.hcu"
#include <iostream>

using namespace std;
using namespace HMMLikeCUDA;

namespace HMMLikeCUDATest {
__global__ void HostUnpackGLs(char GLset, float *GLs) { UnpackGLs(GLset, GLs); }

bool UnpackGLs(char GLset, float *GLs) {

  cudaError_t err = cudaSuccess;

  // figure out how big output will be
  size_t size = 3 * sizeof(float);

  // Allocate the device input GLset
  float *d_GLs = NULL;
  err = cudaMalloc((void **)&d_GLs, size);

  if (err != cudaSuccess) {
    cerr << "Failed to allocate device GLset vector (error code "
         << cudaGetErrorString(err) << ")!\n";

    exit(EXIT_FAILURE);
  }

  HostUnpackGLs << <1, 1>>> (GLset, d_GLs);

  const cudaError_t retErr = cudaGetLastError();

  err = cudaMemcpy(GLs, d_GLs, size, cudaMemcpyDeviceToHost);

  if (err != cudaSuccess) {
    cerr << "Failed to copy vector GLs from device to host (error code "
         << cudaGetErrorString(err) << ")!\n";

    exit(EXIT_FAILURE);
  }

  // deallocate memory
  cudaFree(d_GLs);

  if (retErr == cudaSuccess)
    return true;
  else
    return false;
}

cudaError_t CopyTranToHost(vector<float> &tran) {

  assert(tran.size() == NUMSITES * 3);
  return cudaMemcpyFromSymbol(tran.data(), transitionMat,
                              sizeof(float) * NUMSITES * 3, 0,
                              cudaMemcpyDeviceToHost);
}
cudaError_t CopyMutMatToHost(vector<float> &mutMat) {

  assert(mutMat.size() == 4 * 4);
  return cudaMemcpyFromSymbol(mutMat.data(), mutationMat, sizeof(float) * 4 * 4,
                              0, cudaMemcpyDeviceToHost);
}
}

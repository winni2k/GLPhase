//#include "cuda.h"
//#include "cuda_runtime.h"
#include "hmmLike.hcu"
#include <iostream>

using namespace std;
using namespace HMMLikeCUDA;

namespace HMMLikeCUDATest {

__global__ void GlobalFillEmit(const float *GLs, float *emit) {
  float nEmit[4];
  float nGLs[3];
  for (int i = 0; i < 3; ++i)
    nGLs[i] = GLs[i];
  FillEmit(nGLs, nEmit);
  for (int i = 0; i < 4; ++i)
    emit[i] = nEmit[i];
}

extern "C" void FillEmit(const vector<float> &GLs, vector<float> &emit) {

  assert(GLs.size() == 3);
  assert(emit.size() == 4);
  cudaError_t err = cudaSuccess;
  float *d_GLs;
  err = cudaMalloc(&d_GLs, 3 * sizeof(float));
  assert(err == cudaSuccess);

  err = cudaMemcpy(d_GLs, GLs.data(), 3 * sizeof(float), cudaMemcpyHostToDevice);
  assert(err == cudaSuccess);

  float *d_emit;
  err = cudaMalloc(&d_emit, 4 * sizeof(float));
  assert(err == cudaSuccess);

  GlobalFillEmit << <1, 1>>> (d_GLs, d_emit);

  err = cudaMemcpy(emit.data(), d_emit, 4 * sizeof(float), cudaMemcpyDeviceToHost);
  assert(err == cudaSuccess);
};

__global__ void GlobalUnpackGLs(char GLset, float *GLs) {
  float nGLs[3];
  UnpackGLs(GLset, nGLs);
  for (int i = 0; i < 3; ++i)
    GLs[i] = nGLs[i];
}

extern "C" bool UnpackGLs(char GLset, float *GLs) {

  cudaError_t err = cudaSuccess;

  // figure out how big output will be
  size_t size = 3 * sizeof(float);

  // Allocate the device input GLset
  float *d_GLs;
  err = cudaMalloc(&d_GLs, size);

  if (err != cudaSuccess) {
    cerr << "Failed to allocate device GLset vector (error code "
         << cudaGetErrorString(err) << ")!\n";

    exit(EXIT_FAILURE);
  }

  GlobalUnpackGLs << <1, 1>>> (GLset, d_GLs);

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
  else {
    cerr << "Failed to unpack gls. Error: " << cudaGetErrorString(retErr)
         << "\nError code: " << retErr << "\n";
    exit(EXIT_FAILURE);
  }
}

extern "C" cudaError_t CopyTranToHost(vector<float> &tran) {

  assert(tran.size() == NUMSITES * 3);
  return cudaMemcpyFromSymbol(tran.data(), transitionMat,
                              sizeof(float) * NUMSITES * 3, 0,
                              cudaMemcpyDeviceToHost);
}

extern "C" cudaError_t CopyMutMatToHost(vector<float> &mutMat) {

  assert(mutMat.size() == 4 * 4);
  return cudaMemcpyFromSymbol(mutMat.data(), mutationMat, sizeof(float) * 4 * 4,
                              0, cudaMemcpyDeviceToHost);
}
}

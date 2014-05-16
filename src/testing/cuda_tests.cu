//#include "cuda.h"
//#include "cuda_runtime.h"
#include "hmmLike.hcu"

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

void FillEmit(const vector<float> &GLs, vector<float> &emit) {

  assert(GLs.size() == 3);
  assert(emit.size() == 4);
  cudaError_t err = cudaSuccess;
  float *d_GLs;
  err = cudaMalloc(&d_GLs, 3 * sizeof(float));
  assert(err == cudaSuccess);

  err =
      cudaMemcpy(d_GLs, GLs.data(), 3 * sizeof(float), cudaMemcpyHostToDevice);
  assert(err == cudaSuccess);

  float *d_emit;
  err = cudaMalloc(&d_emit, 4 * sizeof(float));
  assert(err == cudaSuccess);

  GlobalFillEmit << <1, 1>>> (d_GLs, d_emit);

  err = cudaMemcpy(emit.data(), d_emit, 4 * sizeof(float),
                   cudaMemcpyDeviceToHost);
  assert(err == cudaSuccess);
};

__global__ void GlobalHmmLike(unsigned idx, const unsigned (*hapIdxs)[4],
                              const char *__restrict__ d_packedGLs,
                              unsigned packedGLStride,
                              const uint64_t *__restrict__ d_hapPanel,
                              float *retLike) {
  *retLike = hmmLike(idx, *hapIdxs, d_packedGLs, packedGLStride, d_hapPanel);

  return;
}

float CallHMMLike(unsigned idx, const unsigned (*hapIdxs)[4],
                  const vector<char> &packedGLs, unsigned packedGLStride,
                  const vector<uint64_t> &h_hapPanel) {

  cudaError_t err = cudaSuccess;

  /*
  copy initial hap indices to device memory
*/
  cout << "Copying hap indices to device\n";
  // allocate memory on device
  unsigned(*d_hapIdxs)[4];
  err = cudaMalloc(&d_hapIdxs, 4 * sizeof(unsigned));
  if (err != cudaSuccess) {
    cerr << "Failed to allocate\n";
    exit(EXIT_FAILURE);
  }

  // copy data across
  err = cudaMemcpy(d_hapIdxs, hapIdxs, 4 * sizeof(unsigned),
                   cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    cerr << "Failed to copy\n";
    exit(EXIT_FAILURE);
  }

  /*
    copy packedGLs to device memory
  */
  cout << "Copying packed GLs to device\n";
  size_t glSize = packedGLs.size() * sizeof(char);

  // allocate memory on device
  char *d_packedGLs;
  err = cudaMalloc(&d_packedGLs, glSize);
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
  size_t hapPanelSize = h_hapPanel.size() * sizeof(uint64_t);

  // allocate memory on device
  uint64_t *d_hapPanel;
  err = cudaMalloc(&d_hapPanel, hapPanelSize);
  if (err != cudaSuccess) {
    cerr << "Failed to allocate hapPanel on device\n";
    exit(EXIT_FAILURE);
  }

  // copy data across
  err = cudaMemcpy(d_hapPanel, h_hapPanel.data(), hapPanelSize,
                   cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    cerr << "Failed to copy hapPanel to device\n";
    exit(EXIT_FAILURE);
  }

  // create result on device
  thrust::device_vector<float> d_like(1, 1);
  float *d_likePtr = thrust::raw_pointer_cast(&d_like[0]);
  /*
    run kernel
  */
  GlobalHmmLike << <1, 1>>>
      (idx, d_hapIdxs, d_packedGLs, packedGLStride, d_hapPanel, d_likePtr);
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    cerr << "Failed to run HMMLike kernel: " << cudaGetErrorString(err) << "\n";
    exit(EXIT_FAILURE);
  }

  thrust::host_vector<float> h_like = d_like;
  return d_like[0];
};

__global__ void GlobalUnpackGLs(char GLset, float *GLs) {
  float nGLs[3];
  UnpackGLs(GLset, nGLs);
  for (int i = 0; i < 3; ++i)
    GLs[i] = nGLs[i];
}

bool UnpackGLs(char GLset, float *GLs) {

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

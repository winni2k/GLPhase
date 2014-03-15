#ifndef _ENUMS_H
#define _ENUMS_H 1

enum class PanelType { REFERENCE, SCAFFOLD };

enum class RelT {
  sampSampGraph,
  sampHapGraph,
  noGraph,
  kMedoids,
  kNN,
  undefined
};

#endif

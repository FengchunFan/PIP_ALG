//reference: https://github.com/ucrparlay/PIP-algorithms/blob/main/sequence.h#L606
#include <math.h>

#include "parallel.h"
using namespace parlay;

using Type = long long;

template <class T>
T sequ_filter(T *A, T *output, size_t n);
template <class T>
T seq_filter(T *A, size_t n);

template <class T>
T pal_filter(T *A, size_t n) {
  T total = 0;
  Type *output = (Type *)malloc(n * sizeof(Type));
  total = sequ_filter(A, output, n);
  return total;
}

template <class T>
T seq_filter(T *A, size_t n) {
  T total = 0;
  Type *output = (Type *)malloc(n * sizeof(Type));
  total = sequ_filter(A, output, n);
  return total;
}

template <class T>
T sequ_filter(T *A, T *output, size_t n) {
  int k = 0;
	for (size_t i = 0; i < n; i++){
    if (A[i] % 10 == 0) output[k++] = A[i];
  }
  parallel_for(0, n, [&](size_t i) {
    A[i] = output[i];
  });
  return k;
}
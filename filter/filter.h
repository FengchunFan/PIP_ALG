//reference: https://github.com/ucrparlay/PIP-algorithms/blob/main/sequence.h#L606
#include <math.h>

#include "parallel.h"
using namespace parlay;

using Type = long long;

#define nblocks(_n,_bsize) (1 + ((_n)-1)/(_bsize))

template <class T>
T sequ_filter(T *A, T *output, size_t n);
template <class T>
T seq_filter(T *A, size_t n);

template <class T>
T pal_filter(T *A, size_t n) {
  T total = 0;
  if(n < 1e2){
    total = sequ_filter(A, A, n);
    return total;
  }
  
  size_t b = sqrt(n); //block size
  size_t l = nblocks(n, b); //number of blocks

  Type *Sums = (Type *)malloc(l * sizeof(Type));

  parallel_for(0, l, [&](size_t i) {
    size_t s = i * b; //starting point
    size_t e = std::min(s + b, n); //end point, do not exceed n
    size_t k = s; //keep track of number of selected element

    for (size_t j = s; j < e; j++) //loop through the block
			if (A[j] % 10 == 0) A[k++] = A[j]; //if selected, push it to k index of the current block
		Sums[i] = k - s;
  });

  for(size_t i = 0; i < l; ++i){
    total += Sums[i];
  }

  free(Sums);

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
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
template <typename T>
T scan_up_in_place(T *A, size_t s, size_t t);
template <typename T>
void scan_down_in_place(T *A, size_t s, size_t t, T offset);
template <typename T>
T scan_in_place(T *A, size_t n);

template <class T>
T reduce(T *A, size_t n) {
  if (n == 0) {
    return 0;
  } else if (n == 1) {
    return A[0];
  } else {
    T v1, v2;
    auto f1 = [&]() { v1 = reduce(A, n / 2); };
    auto f2 = [&]() { v2 = reduce(A + n / 2, n - n / 2); };
    par_do(f1, f2);
    return v1 + v2;
  }
}

template <typename T>
T pal_filter_in_place(T *A, size_t n) {
  T total = 0;
  if(n < 1e5){
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

  total = reduce(Sums,l);

  size_t start = 0;
  for(size_t i = 0; i < l; ++i){
    for(size_t j = 0; j < (size_t)Sums[i]; ++j){
      A[start] = A[i * b + j];
      start++;
    }
  }

  free(Sums);

  return total;
}

template <typename T>
T pal_filter_non_in_place(T *A, size_t n) {
  T total = 0;
  if(n < 1e5){
    total = sequ_filter(A, A, n);
    return total;
  }
  
  Type *flag = (Type *)malloc(n * sizeof(Type));
  Type *B = (Type *)malloc(n * sizeof(Type));

  parallel_for(0, n, [&](size_t i) {
    flag[i] = (A[i] % 10 == 0);
  });

  scan_in_place((T*)flag, n);

  parallel_for(0, n, [&](size_t i) {
    if(A[i] % 10 == 0){
      B[flag[i]] = A[i];
    }
  });

  parallel_for(0, flag[n-1], [&](size_t i) {
    A[i] = B[i];
  });

  total = flag[n-1];

  free(flag);
  free(B);
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

template <typename T>
T scan_in_place(T *A, size_t n) {
  T total = 0;
  total = scan_up_in_place(A, 0, n-1);
  scan_down_in_place(A, 0, n-1, (T)0);
  return total;
}

template <typename T>
T scan_up_in_place(T *A, size_t s, size_t t) {
  T temp = A[t];
  if(t - s + 1 < 1e5){
    for(size_t i = s; i < t; ++i){
      temp += A[i];
    }
    A[t] = temp; //right here sum is inclusive
    return temp;
  }
  //if(s==t)return A[s];
  T v1, v2;
  auto f1 = [&]() { v1 = scan_up_in_place(A, s, (s+t)/2); };
  auto f2 = [&]() { v2 = scan_up_in_place(A, (s+t)/2 + 1, t); };
  par_do(f1, f2);
  A[t] = v1 + v2;
  return v1 + v2;
}

template <typename T>
void scan_down_in_place(T *A, size_t s, size_t t, T offset) {
  if(t - s + 1 < 1e5){
    T temp = 0;
    T p = offset;
    for (size_t i = s; i <= t; i++) {
      p += temp;
      temp = A[i];
      A[i] = p;
    }
    return;
  }
  /*if(s == t){
    A[s] = offset;
    return;
  }*/
  T leftsum = A[(s+t)/2];
  auto f1 = [&]() {scan_down_in_place(A, s, (s+t)/2, offset);};
  auto f2 = [&]() {scan_down_in_place(A, (s+t)/2 + 1, t, offset + leftsum);};
  par_do(f1, f2);
}
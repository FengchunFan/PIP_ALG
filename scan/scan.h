//reference: https://github.com/ucrparlay/PIP-algorithms/blob/main/sequence.h#L606
#include <math.h>

#include "parallel.h"
using namespace parlay;

using Type = long long;

template <typename T>
T scan_up(T *A, size_t s, size_t t);
template <typename T>
void scan_down(T *A, size_t s, size_t t, T offset);

/**
 * Inplace block scan DIRECTLY without filter
 */
template <typename T>
T scan(T *A, size_t n) {
  /**
   * REPLACE THE BODY OF THE FOLLOWING
   * FUNCTION WITH YOUR PARALLEL IMPLEMENTATION
  */
  
  // note: 1) Extra functions are allowed 
  //       2) DO NOT change the name of this `T scan(T *A, size_t n)`
  
  //Try the divide and conquer approach

  //coarsening
  T total = 0;

  total = scan_up(A, 0, n-1);
  scan_down(A, 0, n-1, (T)0);
  return total;
}

template <typename T>
T scan_up(T *A, size_t s, size_t t) {
  T temp = A[t];
  if(t - s + 1 < 1e4){
    for(size_t i = s; i < t; ++i){
      temp += A[i];
    }
    A[t] = temp;
    return temp;
  }
  if(s==t)return A[s];
  T v1, v2;
  auto f1 = [&]() { v1 = scan_up(A, s, (s+t)/2); };
  auto f2 = [&]() { v2 = scan_up(A, (s+t)/2 + 1, t); };
  par_do(f1, f2);
  A[t] = v1 + v2;
  return v1 + v2;
}

template <typename T>
void scan_down(T *A, size_t s, size_t t, T offset) {
  if(t - s + 1 < 1e4){
    T temp = 0;
    T p = offset;
    for (size_t i = s; i <= t; i++) {
      p += temp;
      temp = A[i];
      A[i] = p;
    }
    return;
  }
  if(s == t){
    A[s] = offset;
    return;
  }
  T leftsum = A[(s+t)/2];
  auto f1 = [&]() {scan_down(A, s, (s+t)/2, offset);};
  auto f2 = [&]() {scan_down(A, (s+t)/2 + 1, t, offset + leftsum);};
  par_do(f1, f2);
}
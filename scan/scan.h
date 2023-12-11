//reference: https://github.com/ucrparlay/PIP-algorithms/blob/main/sequence.h#L606
#include <math.h>

#include "parallel.h"
using namespace parlay;

using Type = long long;

template <typename T>
T scan_up_in_place(T *A, size_t s, size_t t);
template <typename T>
void scan_down_in_place(T *A, size_t s, size_t t, T offset);

template <typename T>
T scan_up(T *A, T *LS, size_t n);
template <typename T>
void scan_down(T *A, T *LS, size_t n, T offset);

/**
 * Inplace block scan DIRECTLY without filter
 */
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

template <typename T>
T scan(T *A, size_t n) {
  T total = 0;
  Type *LS = (Type *)malloc((n-1) * sizeof(Type));

  total = scan_up(A, LS, n);
  scan_down(A, LS, n, (T)0);

  free(LS);
  return total;
}

//builds	a	tree	of	sums, only left branch involved
//similar to a reduce function but store information
template <typename T>
T scan_up(T *A, T *LS, size_t n) {
  if (n < 1e5){
    T total = 0;
    for (size_t i = 0; i < n; i++) {
      total += A[i];
    }
    return total;
  }else{
    size_t m = n/2;
    T v1, v2;
    auto f1 = [&]() { v1 = scan_up(A, LS, m); };
    auto f2 = [&]() { v2 = scan_up(A + m, LS + m, n - m); };
    par_do(f1, f2);
    LS[m-1] = v1;
    return v1 + v2;
  }
}

//traverses	the	tree to	compute	prefixes
template <typename T>
void scan_down(T *A, T *LS, size_t n, T offset) {
  if(n < 1e5){
    T total = 0;
    for (size_t i = 0; i < n; i++) {
      T tmp = A[i];
      A[i] = total + offset;
      total += tmp;
    }
  }else{
    size_t m = n/2;
    auto f1 = [&]() {scan_down(A, LS, m, offset);};
    auto f2 = [&]() {scan_down(A + m, LS + m, n - m, offset + LS[m-1]);};
    par_do(f1, f2);
  }
}
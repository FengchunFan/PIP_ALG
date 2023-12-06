//reference: https://github.com/ucrparlay/PIP-algorithms/blob/main/sequence.h#L606
#include <math.h>
#include "utils.h"
#include "parallel.h"
using namespace parlay;

using Type = long long;

template <class T, class size_t>
struct getA {
		T* A;
		getA(T* AA) : A(AA) {}
		T operator() (size_t i) { return A[i]; }
	};

#define nblocks(_n,_bsize) (1 + ((_n)-1)/(_bsize))

#define blocked_for(_i, _s, _e, _bsize, _body)  {	\
    size_t _ss = _s;					\
    size_t _ee = _e;					\
    size_t _n = _ee-_ss;					\
    size_t _l = nblocks(_n,_bsize);			\
    parallel_for(0, _l, [&](size_t _i){     \
      size_t _s = _ss + _i * (_bsize);			\
      size_t _e = min(_s + (_bsize), _ee);			\
      _body						\
	});						\
  }

template <class T, class size_t, class F, class G>
	T reduceSerial(size_t s, size_t e, F f, G g) {
		T r = g(s);
		for (size_t j = s + 1; j < e; j++) r = f(r, g(j));
		return r;
	}

template <class T>
T sequ_filter(T *A, T *output, size_t n);
template <class T>
T seq_filter(T *A, size_t n);
template <class T>
int binary_search(T *A, size_t &target, size_t &start, size_t &end);

template <class T, class size_t, class F, class G>
	T scan(T* Out, T s, T e, F f, G g, T zero, bool inclusive, bool back) {
		size_t n = e - s;
		size_t l = nblocks(n, sqrt(n));
		if (l <= 2) return scanSerial(Out, s, e, f, g, zero, inclusive, back);
    Type *Sums = (Type *)malloc(nblocks(n, sqrt(n)) * sizeof(Type));
		blocked_for(i, s, e, 1e4,
			Sums[i] = reduceSerial<T>(s, e, f, g););
		T total = scan(Sums, (size_t)0, l, f, getA<T, size_t>(Sums), zero, false, back);
		blocked_for(i, s, e, 1e4,
			scanSerial(Out, s, e, f, g, Sums[i], inclusive, back););
		free(Sums);
		return total;
	}

template <class T, class size_t>
T plusScan(T *In, T* Out, size_t n) {
  return scan(Out, (T)0, n, utils::addF<T>(), getA<T, size_t>(In), (T)0, false, false);
}

template <class T>
T pal_filter(T *A, size_t n) {
  T total = 0;
  if(n < 1e4){
    total = sequ_filter(A, A, n);
  }
  
  size_t b = sqrt(n); //block size
  size_t l = nblocks(n, b);

  Type *Sums = (Type *)malloc(l+1 * sizeof(Type));

  parallel_for(0, l, [&](size_t i) {
    size_t s = i * b; //starting point
    size_t e = std::min((size_t)s + (size_t)b, (size_t)n); //end point, do not exceed n
    size_t k = s;

    for (size_t j = s; j < e; j++)
			if (A[j] % 10 == 0) A[k++] = A[j];
		Sums[i] = k - s;
  });

  size_t m = plusScan(Sums, Sums, l);
	Sums[l] = m;

	size_t startblock = 1;
	while (startblock < l) {
		size_t startplace = b * startblock;
		size_t endblock = binary_search(Sums, startplace, startblock, l);
		parallel_for(startblock, endblock+1, [&](size_t i) {
			size_t S = i*b;
			size_t D = Sums[i];
			size_t L = Sums[i + 1] - D;
			for (size_t j = 0; j < L; j++) {
				A[D + j] = A[S + j];
			}
		});
		startblock = endblock + 1;
			//reportTime();
	}
	free(Sums);
	return m;
}

template <class T>
int binary_search(T *A, size_t &target, size_t &start, size_t &end){
  if (start + 1 >= end)	return start;
	T temp = A[(start + end) / 2 + 1];
	if (temp > target) {
		return binary_search(A, target, start, (start + end) / 2);
	}	else {
		return binary_search(A, target, (start + end) / 2, end);
	}
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
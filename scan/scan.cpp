#include <cstdlib>
#include <iostream>

#include "scan.h"
#include "get_time.h"
#include "scan_seq.h"

using Type = long long;

int main(int argc, char *argv[]) {
  size_t n = 1e9;
  //size_t n = 1e4;
  int num_rounds = 3;
  if (argc >= 2) {
    n = atoll(argv[1]);
  }
  if (argc >= 3) {
    num_rounds = atoi(argv[2]);
  }
  Type *seq_scan_A = (Type *)malloc(n * sizeof(Type));
  Type *pal_scan_A = (Type *)malloc(n * sizeof(Type));
  Type *pal_scan_B = (Type *)malloc(n * sizeof(Type));
  parallel_for(0, n, [&](size_t i) {
    seq_scan_A[i] = i;
    pal_scan_A[i] = i;
    pal_scan_B[i] = i;
  });

  double total_time = 0;
  std::cout << "******* Parallel In Place Results: ******* " << std::endl;
  for (int i = 0; i <= num_rounds; i++) {
    parlay::timer t;
    long long exclusive_sum = scan_in_place(pal_scan_A, n);
    t.stop();
    if (i == 0) {
      std::cout << "Exclusive sum (Parallel In Place): " << exclusive_sum << std::endl;
      std::cout << "Warmup round running time: " << t.total_time() << std::endl;
    } else {
      std::cout << "Round " << i << " running time: " << t.total_time()
                << std::endl;
      total_time += t.total_time();
    }
  }

  double total_time_2 = 0;
  std::cout << "******* Parallel Non In Place Results: ******* " << std::endl;
  for (int i = 0; i <= num_rounds; i++) {
    parlay::timer t;
    long long exclusive_sum = scan(pal_scan_B, n);
    t.stop();
    if (i == 0) {
      std::cout << "Exclusive sum (Parallel non in place): " << exclusive_sum << std::endl;
      std::cout << "Warmup round running time: " << t.total_time() << std::endl;
    } else {
      std::cout << "Round " << i << " running time: " << t.total_time()
                << std::endl;
      total_time_2 += t.total_time();
    }
  }

  std::cout << "\n******* Sequential Results: ******* " << std::endl;
  double seq_total_time = 0;
  for (int i = 0; i <= num_rounds; i++) {
    parlay::timer seq_t;
    long long seq_ans = seq_scan(seq_scan_A, n);
    seq_t.stop();

    if (i == 0) {
      std::cout << "Exclusive sum (Sequential): " << seq_ans << std::endl;
      std::cout << "Warmup round running time (Sequential): "
                << seq_t.total_time() << std::endl;
    } else {
      std::cout << "Round " << i
                << " running time (sequential): " << seq_t.total_time()
                << std::endl;
      seq_total_time += seq_t.total_time();
    }
  }
  std::cout << "\n******* Summary: ******* " << std::endl;
  std::cout << "Average sequential running time: "
            << seq_total_time / num_rounds << std::endl;
  std::cout << "Average parallel in place running time: " << total_time / num_rounds
            << std::endl;
  std::cout << "Average parallel non in place running time: " << total_time_2 / num_rounds
            << std::endl;
  
  /*std::cout << "seq results: ";
  for (size_t i = 0; i < n; i++) {
    std::cout << seq_scan_A[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "pal results: ";
  for (size_t i = 0; i < n; i++) {
    std::cout << pal_scan_A[i] << " ";
  }
  std::cout << std::endl;*/
  for (size_t i = 0; i < n; i++) {
    if (seq_scan_A[i] != pal_scan_A[i]) {
      std::cout << "***********************" << std::endl;
      std::cout << "**** Wrong answer 1! ****" << std::endl;
      std::cout << "***** index: "<< i << " ****" << std::endl;
      std::cout << "***********************" << std::endl;
      break;
    }
  }

  for (size_t i = 0; i < n; i++) {
    if (pal_scan_A[i] != pal_scan_B[i]) {
      std::cout << "***********************" << std::endl;
      std::cout << "**** Wrong answer 2! ****" << std::endl;
      std::cout << "***** index: "<< i << " ****" << std::endl;
      std::cout << "***********************" << std::endl;
      break;
    }
  }

  free(seq_scan_A);
  free(pal_scan_A);
  free(pal_scan_B);
  return 0;
}

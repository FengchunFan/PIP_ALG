#include <cstdlib>
#include <iostream>

//#include "filter_real.h"
#include "filter.h"
#include "get_time.h"

using Type = long long;

int main(int argc, char *argv[]) {
  size_t n = 1e9;
  //size_t n = 1e4;
  int num_rounds = 0;
  long long pal_ans = 0;
  long long pal2_ans = 0;
  long long seq_ans = 0;
  if (argc >= 2) {
    n = atoll(argv[1]);
  }
  if (argc >= 3) {
    num_rounds = atoi(argv[2]);
  }
  Type *seq_filter_A = (Type *)malloc(n * sizeof(Type));
  Type *pal_filter_A = (Type *)malloc(n * sizeof(Type));
  Type *pal_filter_B = (Type *)malloc(n * sizeof(Type));

  parallel_for(0, n, [&](size_t i) {
    seq_filter_A[i] = i;
    pal_filter_A[i] = i;
    pal_filter_B[i] = i;
  });

  double total_time = 0;
  std::cout << "******* Parallel In Place Results: ******* " << std::endl;
  for (int i = 0; i <= num_rounds; i++) {
    parlay::timer t;
    pal_ans = pal_filter_in_place(pal_filter_A, n);
    t.stop();
    if (i == 0) {
      std::cout << "Number of selected elements (in place Parallel): " << pal_ans << std::endl;
      std::cout << "Warmup round running time: " << t.total_time() << std::endl;
    } else {
      std::cout << "Round " << i << " running time: " << t.total_time()
                << std::endl;
      total_time += t.total_time();
    }
  }

  double total_time_2 = 0;
  std::cout << "\n******* Parallel Non In Place Results: ******* " << std::endl;
  for (int i = 0; i <= num_rounds; i++) {
    parlay::timer t;
    pal2_ans = pal_filter_non_in_place(pal_filter_B, n);
    t.stop();
    if (i == 0) {
      std::cout << "Number of selected elements (non in place Parallel): " << pal2_ans << std::endl;
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
    seq_ans = seq_filter(seq_filter_A, n);
    seq_t.stop();

    if (i == 0) {
      std::cout << "Number of selected elements (Sequential): " << seq_ans << std::endl;
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
            << seq_total_time / 1 << std::endl;
  std::cout << "Average parallel in place running time: " << total_time / 1
            << std::endl;
  std::cout << "Average parallel non in place running time: " << total_time_2 / 1
            << std::endl;
  
  /*std::cout << "seq results: ";
  for (size_t i = 0; i < 20; i++) {
    std::cout << seq_filter_A[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "pal results: ";
  for (size_t i = 0; i < 20; i++) {
    std::cout << pal_filter_A[i] << " ";
  }*/
  /*std::cout << std::endl;
  std::cout << "seq ans: " << seq_ans << ". pal ans: " << pal_ans;
  std::cout << std::endl;*/

  if(seq_ans == pal_ans && pal_ans == pal2_ans){
    for (size_t i = 0; i < (size_t)seq_ans - 1; i++) {
      if (seq_filter_A[i] != pal_filter_A[i] || pal_filter_B[i] != pal_filter_A[i]) {
        std::cout << "***********************" << std::endl;
        std::cout << "**** Wrong answer! ****" << std::endl;
        std::cout << "***** index: "<< i << " ****" << std::endl;
        std::cout << "***********************" << std::endl;
        break;
      }
    }
  }else{
    std::cout << "**** Wrong answer! ****" << std::endl;
  }

  free(seq_filter_A);
  free(pal_filter_A);
  free(pal_filter_B);
  return 0;
}

#include <bits/stdc++.h>

#include "matrix/column_major_matrix.h"
#include "matrix/row_major_matrix.h"

using namespace std;

#define N 500
#define M 4000
#define L 2000
int main(int argc, char **argv) {

  Row_Major_Matrix<int> rr(N, M);
  Column_Major_Matrix<int> cc(M, L);

  auto par_bg = chrono::high_resolution_clock::now();
  Row_Major_Matrix<int> rr2_par = rr % cc;
  auto par_ed = chrono::high_resolution_clock::now();
  cout << chrono::duration_cast<chrono::duration<double>>(par_ed - par_bg).count() << '\n';

  auto seq_bg = chrono::high_resolution_clock::now();
  Row_Major_Matrix<int> rr2_seq = rr * cc;
  auto seq_ed = chrono::high_resolution_clock::now();
  cout << chrono::duration_cast<chrono::duration<double>>(seq_ed - seq_bg).count() << '\n';

  for (int i = 0; i < N; i++)
    for (int j = 0; j < L; j++)
      assert(rr2_seq[i][j] == rr2_par[i][j]);
}

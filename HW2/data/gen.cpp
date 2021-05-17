#include <bits/stdc++.h>
using namespace std;
int main() {
  int n = 1000;
  int m = 1000;
  cout << ">Seq1\n";
  while (n--)
    cout << "ATCG"[rand() % 4];
  cout << '\n';
  cout << ">Seq2\n";
  while (m--)
    cout << "ATCG"[rand() % 4];
  cout << '\n';
}

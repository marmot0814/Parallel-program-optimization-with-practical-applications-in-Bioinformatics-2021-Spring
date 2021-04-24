#include <bits/stdc++.h>
using namespace std;

template <typename T>
vector<int> discretization(const vector<T> &vec) {
  vector<T> pool = vec;
  sort(pool.begin(), pool.end());
  pool.erase(unique(pool.begin(), pool.end()), pool.end());

  vector<int> ret;
  for (auto &v : vec)
    ret.push_back(
      lower_bound(pool.begin(), pool.end(), v) - pool.begin() + 1
    );
  return ret;
}

#define L_TYPE 0
#define S_TYPE 1

bool is_LMS(int i, const vector<bool> &type) {
  return i > 0 and type[i - 1] == L_TYPE and type[i] == S_TYPE;
}

void induceLMS(
  vector<int> &sa,
  const vector<int> &s,
  const vector<bool> &type,
  const vector<int> &cnt
) {
  int len = s.size(), sig = cnt.size();
  vector<int> ptr(sig, 0);
  for (int i = 1; i < sig; i++)
    ptr[i] = cnt[i] - 1;

  for (int i = 0; i < len; i++)
    if (is_LMS(i, type))
      sa[ptr[s[i]]--] = i;
}
void induceL(
  vector<int> &sa,
  const vector<int> &s,
  const vector<bool> &type,
  const vector<int> &cnt
) {
  int len = s.size(), sig = cnt.size();
  vector<int> ptr(sig, 0);
  for (int i = 1; i < sig; i++)
    ptr[i] = cnt[i - 1];

  for (int i = 0; i < len; i++)
    if (sa[i] > 0 && type[sa[i] - 1] == L_TYPE)
      sa[ptr[s[sa[i] - 1]]++] = sa[i] - 1;
}
void induceS(
  vector<int> &sa,
  const vector<int> &s,
  const vector<bool> &type,
  const vector<int> &cnt
) {
  int len = s.size(), sig = cnt.size();
  vector<int> ptr(sig, 0);
  for (int i = 1; i < sig; i++)
    ptr[i] = cnt[i] - 1;

  for (int i = len - 1; i >= 0; i--)
    if (sa[i] > 0 && type[sa[i] - 1] == S_TYPE)
      sa[ptr[s[sa[i] - 1]]--] = sa[i] - 1;
}

template <typename T>
vector<int> sais(const vector<T> &vec) {
  auto s = discretization(vec);
  s.push_back(0);

  int len = s.size();
  
  vector<bool> type(len, S_TYPE);
  for (int i = len - 2; i >= 0; i--) {
    if (s[i] < s[i + 1])
      type[i] = S_TYPE;
    else if (s[i] > s[i + 1])
      type[i] = L_TYPE;
    else
      type[i] = type[i + 1];
  }

  int sig = max_element(s.begin(), s.end()) - s.begin() + 1;
  vector<int> cnt(sig, 0);
  for (int i = 0; i < len; i++)
    cnt[s[i]]++;
  for (int i = 1; i < sig; i++)
    cnt[i] += cnt[i - 1];

  vector<int> sa(len, -1);
  induceLMS(sa, s, type, cnt);
  induceL(sa, s, type, cnt);
  induceS(sa, s, type, cnt);
  return sa;
}

vector<int> sais(const string &s) {
  return sais(vector<char>(s.begin(), s.end()));
}

int main() {
  auto sa = sais("immissiissippi");
  for (auto &v : sa)
    cout << v + 1 << ' ';
  cout << '\n';

}

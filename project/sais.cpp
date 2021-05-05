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

vector<int> get_s_ptr(const vector<int> &cnt) {
  vector<int> ptr(cnt.size(), 0);
  for (int i = 1; i < (int)ptr.size(); i++)
    ptr[i] = cnt[i] - 1;
  return ptr;
}

vector<int> get_l_ptr(const vector<int> &cnt) {
  vector<int> ptr(cnt.size(), 0);
  for (int i = 1; i < (int)ptr.size(); i++)
    ptr[i] = cnt[i - 1];
  return ptr;
}

bool is_lms(int i, const vector<bool> &type) {
  return i > 0 and type[i - 1] == L_TYPE and type[i] == S_TYPE;
}

vector<bool> get_suffix_type(const vector<int> &s) {
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
  return type;
}


/* === print === */
void print_sa_header(auto &sa, auto &cnt) {
  cerr << "Bucket:\t";
  auto l_ptr = get_l_ptr(cnt);
  int ptr = 0;
  for (int i = 0; i < (int)sa.size(); i++) {
    if (l_ptr[ptr] == i)
      cerr << ptr++;
    else
      cerr << ' ';
    cerr << "   ";
  }
  cerr << '\n';
}

void print_sa_content(auto &sa, auto &cnt) {
  cerr << "SA:\t";
  auto l_ptr = get_l_ptr(cnt);
  int ptr = 0;
  auto s_ptr = get_s_ptr(cnt);
  for (int i = 0; i < (int)sa.size(); i++) {
    if (l_ptr[ptr] == i) {
      cerr << "{";
    } else {
      cerr << " ";
    }
    cerr << setw(2) << setfill('0') << sa[i];
    if (s_ptr[ptr] == i) {
      cerr << "}", ptr++;
    } else {
      cerr << " ";
    }
  }
  cerr << '\n';
}

void print_sa_ptr(auto &sa, vector<pair<int, char>> idx) {
  cerr << '\t';
  for (int i = 0; i < (int)sa.size(); i++) {
    cerr << " ";
    bool printed = false;
    for (auto &p : idx) {
      if (i == p.first) {
        cerr << p.second;
        printed = true;
      }
    }
    if (not printed)
      cerr << " ";
    cerr << "  ";
  }
  cerr << '\n';
}

void print_s_data(const vector<int> &s) {
  cerr << "Index:\t";
  for (int i = 0; i < (int)s.size(); i++) {
    cerr << " ";
    cerr << setw(2) << setfill('0') << i;
    cerr << " ";
  }
  cerr << '\n';
  cerr << "S:\t";
  for (int i = 0; i < (int)s.size(); i++) {
    cerr << "  ";
    cerr << s[i];
    cerr << " ";
  }
  cerr << '\n';

  auto type = get_suffix_type(s);
  cerr << "T:\t";
  for (int i = 0; i < (int)s.size(); i++) {
    cerr << "  ";
    cerr << (type[i] == L_TYPE ? "L" : "S");
    cerr << " ";
  }
  cerr << '\n';

  cerr << "LMS:\t";
  for (int i = 0; i < (int)s.size(); i++) {
    cerr << "  ";
    cerr << (is_lms(i, type) ? '*' : ' ');
    cerr << " ";
  }
  cerr << '\n';
}
/* === print === */

vector<int> get_suffix_cnt(const vector<int> &s) {
  int sigma = *max_element(s.begin(), s.end()) + 1;
  vector<int> cnt(sigma, 0);
  for (int i = 0; i < (int)s.size(); i++)
    cnt[s[i]]++;
  for (int i = 1; i < sigma; i++)
    cnt[i] += cnt[i - 1];
  return cnt;
}

vector<int> get_lms_suffix(
  const vector<int> &s,
  const vector<bool> &type
) {
  vector<int> lms;
  for (int i = 0; i < (int)s.size(); i++)
    if (is_lms(i, type))
      lms.push_back(i);
  return lms;
}

void induce_L(
  vector<int> &sa,
  const vector<int> &s,
  const vector<bool> &type,
  const vector<int> &cnt
) {
  cerr << "\ninduce L:\n";
  print_s_data(s);
  print_sa_header(sa, cnt);
  auto l_ptr = get_l_ptr(cnt);
  for (int i = 0; i < (int)s.size(); i++) {
    if (sa[i] > 0 && type[sa[i] - 1] == L_TYPE) {
      sa[l_ptr[s[sa[i] - 1]]++] = sa[i] - 1;
      print_sa_content(sa, cnt);
      print_sa_ptr(sa, {{l_ptr[s[sa[i] - 1]] - 1, '^'}, {i, '@'}});
    }
  }
}

void induce_S(
  vector<int> &sa,
  const vector<int> &s,
  const vector<bool> &type,
  const vector<int> &cnt
) {
  cerr << "\ninduce S:\n";
  print_s_data(s);
  print_sa_header(sa, cnt);
  auto s_ptr = get_s_ptr(cnt);
  for (int i = (int)s.size() - 1; i >= 0; i--) {
    if (sa[i] > 0 && type[sa[i] - 1] == S_TYPE) {
      sa[s_ptr[s[sa[i] - 1]]--] = sa[i] - 1;
      print_sa_content(sa, cnt);
      print_sa_ptr(sa, {{s_ptr[s[sa[i] - 1]] + 1, '^'}, {i, '@'}});
    }
  }
}

bool is_equal_lms(
  const vector<int> &s,
  int x, int y,
  const vector<bool> &type
) {
  do {
    if (s[x] != s[y])
      return false;
    x++, y++;
  } while (not is_lms(x, type) and not is_lms(y, type));
  return s[x] == s[y];
}

vector<int> get_s1(
  const vector<int> &sa,
  const vector<int> &s,
  const vector<bool> &type
) {
  vector<int> label(s.size(), -1);
  int prev_x = -1, label_cnt = 1;
  for (int i = 1; i < (int)s.size(); i++) {
    int x = sa[i];
    if (not is_lms(x, type))
      continue;

    if (prev_x >= 0 && not is_equal_lms(s, x, prev_x, type))
      label_cnt++;

    label[x] = label_cnt;
    prev_x = x;
  }
  label.back() = 0;

  vector<int> s1;
  for (int i = 0; i < (int)s.size(); i++) {
    if (label[i] == -1)
      continue;
    s1.push_back(label[i]);
  }
  return s1;
}

vector<int> sais_core(const vector<int> &s) {

  cerr << "\nsais_core\n";
  auto type = get_suffix_type(s);
  auto cnt = get_suffix_cnt(s);

  auto lmss = get_lms_suffix(s, type);

  vector<int> sa(s.size(), -1);
  auto s_ptr = get_s_ptr(cnt);

  cerr << "induce lms\n";
  print_s_data(s);
  print_sa_header(sa, cnt);
  for (auto &lms : lmss) {
    sa[s_ptr[s[lms]]--] = lms;
    print_sa_content(sa, cnt);
    print_sa_ptr(sa, {{s_ptr[s[lms]] + 1, '^'}});
  }

  induce_L(sa, s, type, cnt);
  induce_S(sa, s, type, cnt);

  auto s1 = get_s1(sa, s, type);

  vector<int> sa1(s1.size(), -1);
  for (int i = 0; i < (int)s1.size(); i++) {
    if (sa1[s1[i]] != -1) {
      sa1 = sais_core(s1);
      break;
    }
    sa1[s1[i]] = i;
  }

  fill(sa.begin(), sa.end(), -1);
  s_ptr = get_s_ptr(cnt);

  cerr << "induce lms\n";
  print_s_data(s);
  print_sa_header(sa, cnt);
  for (int i = (int)lmss.size() - 1; i >= 0; i--) {
    sa[s_ptr[s[lmss[sa1[i]]]]--] = lmss[sa1[i]];
    print_sa_content(sa, cnt);
    print_sa_ptr(sa, {{s_ptr[s[lmss[sa1[i]]]] + 1, '^'}});
  }

  induce_L(sa, s, type, cnt);
  induce_S(sa, s, type, cnt);

  return sa;
}

template <typename T>
vector<int> sais(const vector<T> &vec) {
  auto s = discretization(vec);
  s.push_back(0);

  return sais_core(s);
}

vector<int> sais(const string &s) {
  return sais(vector<char>(s.begin(), s.end()));
}

int main() {
  string s; cin >> s;
  auto sa = sais(s);
}
// GATCCGCATCGCGATCG

// GATCCGATCGATCG


// CATGGTATGGT

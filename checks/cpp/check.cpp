// check.cpp -- independent check of a certificate: every relation is an exact identity that
// places x_p in the interior of a nondegenerate simplex, and every k-subset of the points
// contains the support of some relation (Proposition II.5 of the paper).
//
//   c++ -O2 -std=c++17 -o check check.cpp
//   ./check ../../highdim/witness_r4_11.json ../../highdim/cert_r4_11.json 4 8
//
// Both files are read as flat lists of integers: the points as n*d coordinates, the
// certificate as records (p, s_0..s_d, m_0..m_d, orient), which is the order of the keys
// in the JSON files of the repository.  Indices are 0-based.  Exit status 0 iff all
// relations and the covering condition hold.
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

typedef __int128 i128;

static std::vector<long long> ints(const char* path) {
  std::ifstream in(path);
  if (!in) { std::cerr << "cannot read " << path << "\n"; std::exit(2); }
  std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  std::vector<long long> v;
  for (size_t i = 0; i < s.size();) {
    if (std::isdigit((unsigned char)s[i])) {
      bool neg = i > 0 && s[i - 1] == '-';
      long long x = 0;
      while (i < s.size() && std::isdigit((unsigned char)s[i])) x = x * 10 + (s[i++] - '0');
      v.push_back(neg ? -x : x);
    } else {
      ++i;
    }
  }
  return v;
}

// determinant by cofactor expansion along the first row (matrices of size at most 6)
static i128 det(const std::vector<std::vector<i128>>& a) {
  size_t m = a.size();
  if (m == 1) return a[0][0];
  i128 s = 0;
  for (size_t c = 0; c < m; ++c) {
    std::vector<std::vector<i128>> b;
    for (size_t r = 1; r < m; ++r) {
      std::vector<i128> row;
      for (size_t j = 0; j < m; ++j) if (j != c) row.push_back(a[r][j]);
      b.push_back(row);
    }
    i128 t = a[0][c] * det(b);
    s += (c % 2) ? -t : t;
  }
  return s;
}

int main(int argc, char** argv) {
  if (argc != 5) { std::cerr << "usage: check points.json cert.json d k\n"; return 2; }
  int d = std::atoi(argv[3]), k = std::atoi(argv[4]);
  std::vector<long long> pv = ints(argv[1]), cv = ints(argv[2]);
  if (d < 1 || pv.size() % d) { std::cerr << "bad points file\n"; return 2; }
  int n = (int)(pv.size() / d);
  size_t rec = 1 + 2 * (d + 1) + 1;
  if (n > 30 || cv.size() % rec) { std::cerr << "bad certificate file\n"; return 2; }
  auto X = [&](int i, int j) { return (i128)pv[(size_t)i * d + j]; };
  auto rowdet = [&](const std::vector<int>& idx) {
    std::vector<std::vector<i128>> a;
    for (int i : idx) {
      std::vector<i128> row{1};
      for (int j = 0; j < d; ++j) row.push_back(X(i, j));
      a.push_back(row);
    }
    return det(a);
  };

  int bad = 0;
  std::vector<unsigned> supports;
  for (size_t r = 0; r < cv.size(); r += rec) {
    int p = (int)cv[r];
    std::vector<int> S;
    std::vector<i128> m;
    for (int j = 0; j <= d; ++j) S.push_back((int)cv[r + 1 + j]);
    for (int j = 0; j <= d; ++j) m.push_back(cv[r + 2 + d + j]);
    i128 stored = cv[r + rec - 1], M = 0;
    unsigned mask = 1u << p;
    bool good = p >= 0 && p < n;
    for (int j = 0; j <= d; ++j) {
      good = good && S[j] >= 0 && S[j] < n && S[j] != p && !(mask >> S[j] & 1) && m[j] > 0;
      if (good) mask |= 1u << S[j];
      M += m[j];
    }
    for (int c = 0; good && c < d; ++c) {
      i128 lhs = 0;
      for (int j = 0; j <= d; ++j) lhs += m[j] * X(S[j], c);
      good = lhs == M * X(p, c);
    }
    i128 D = good ? rowdet(S) : 0;
    good = good && D != 0 && D == stored && D % M == 0;
    // Cramer: the numerators N_j (x_p in place of x_{s_j}) satisfy m_j D = M N_j, all of the sign of D
    for (int j = 0; good && j <= d; ++j) {
      std::vector<int> T = S;
      T[j] = p;
      i128 N = rowdet(T);
      good = m[j] * D == M * N && (N > 0) == (D > 0);
    }
    if (!good) { ++bad; std::cerr << "bad relation at record " << r / rec << "\n"; }
    supports.push_back(mask);
  }

  long long subsets = 0, uncovered = 0;
  for (unsigned T = 0; T < (1u << n); ++T) {
    if (__builtin_popcount(T) != k) continue;
    ++subsets;
    bool hit = false;
    for (unsigned s : supports) if ((s & T) == s) { hit = true; break; }
    if (!hit) ++uncovered;
  }
  std::printf("%d points in R^%d, %zu relations: %d bad; %lld %d-subsets, %lld without a support\n",
              n, d, supports.size(), bad, subsets, k, uncovered);
  bool ok = bad == 0 && uncovered == 0 && !supports.empty();
  std::printf("%s\n", ok ? "OK: certificate valid" : "FAILED");
  return ok ? 0 : 1;
}

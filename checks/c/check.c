/* check.c -- independent brute-force check, in exact integer arithmetic, that a set of
 * integer points in R^d is in general position and has no k points in convex position.
 * It does not read the certificate.
 *
 *   cc -O2 -o check check.c && ./check ../../highdim/witness_r4_11.json 4 8
 *
 * The points file is read as a flat list of integers (the JSON list of lists of the
 * repository); its length must be a multiple of d.  A subset Y is in convex position iff
 * no point of Y lies in the interior of a simplex of d+1 other points of Y (Caratheodory;
 * Lemma II.1 of the paper).  A point p lies in the interior of the simplex S iff replacing
 * any vertex of S by p keeps the nonzero orientation sign of S.  As a control the program
 * also counts the (k-1)-subsets in convex position, which must be positive.
 * Exit status 0 iff general position holds, there is no k-gon and the control is positive. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAXN 20
#define MAXD 6
typedef __int128 i128;

static int n, d, k;
static long long P[MAXN][MAXD];
static signed char *sgn; /* orientation sign of each sorted (d+1)-subset, indexed by bitmask */

static i128 det(i128 a[MAXD + 1][MAXD + 1], int m) { /* Bareiss, exact */
  i128 prev = 1; int s = 1;
  for (int i = 0; i < m - 1; i++) {
    if (a[i][i] == 0) {
      int r;
      for (r = i + 1; r < m && a[r][i] == 0; r++) ;
      if (r == m) return 0;
      for (int c = 0; c < m; c++) { i128 t = a[i][c]; a[i][c] = a[r][c]; a[r][c] = t; }
      s = -s;
    }
    for (int r = i + 1; r < m; r++)
      for (int c = i + 1; c < m; c++)
        a[r][c] = (a[r][c] * a[i][i] - a[r][i] * a[i][c]) / prev;
    prev = a[i][i];
  }
  return s * a[m - 1][m - 1];
}

static int orient(const int *idx) {
  i128 a[MAXD + 1][MAXD + 1];
  for (int r = 0; r <= d; r++) {
    a[r][0] = 1;
    for (int c = 0; c < d; c++) a[r][c + 1] = P[idx[r]][c];
  }
  i128 v = det(a, d + 1);
  return v > 0 ? 1 : (v < 0 ? -1 : 0);
}

/* orientation of an arbitrary ordered (d+1)-tuple from the table of sorted ones */
static int orient_tuple(const int *t) {
  int mask = 0, inv = 0;
  for (int i = 0; i <= d; i++) {
    mask |= 1 << t[i];
    for (int j = i + 1; j <= d; j++) inv += t[i] > t[j];
  }
  return (inv & 1) ? -sgn[mask] : sgn[mask];
}

static int inside(int p, const int *S) {
  int s0 = orient_tuple(S), T[MAXD + 1];
  for (int j = 0; j <= d; j++) {
    memcpy(T, S, sizeof(int) * (d + 1));
    T[j] = p;
    if (orient_tuple(T) != s0) return 0;
  }
  return 1;
}

/* does the set Y (size m) contain a point inside a simplex of d+1 others? */
static int not_convex(const int *Y, int m) {
  int S[MAXD + 1], rest[MAXN];
  for (int a = 0; a < m; a++) {
    int r = 0;
    for (int b = 0; b < m; b++) if (b != a) rest[r++] = Y[b];
    int c[MAXD + 1];
    for (int i = 0; i <= d; i++) c[i] = i;
    for (;;) {
      for (int i = 0; i <= d; i++) S[i] = rest[c[i]];
      if (inside(Y[a], S)) return 1;
      int i = d;
      while (i >= 0 && c[i] == r - d - 1 + i) i--;
      if (i < 0) break;
      c[i]++;
      for (int j = i + 1; j <= d; j++) c[j] = c[j - 1] + 1;
    }
  }
  return 0;
}

static long long count_convex(int m) {
  int c[MAXN], Y[MAXN];
  long long cnt = 0;
  for (int i = 0; i < m; i++) c[i] = i;
  for (;;) {
    for (int i = 0; i < m; i++) Y[i] = c[i];
    if (!not_convex(Y, m)) cnt++;
    int i = m - 1;
    while (i >= 0 && c[i] == n - m + i) i--;
    if (i < 0) break;
    c[i]++;
    for (int j = i + 1; j < m; j++) c[j] = c[j - 1] + 1;
  }
  return cnt;
}

int main(int argc, char **argv) {
  if (argc != 4) { fprintf(stderr, "usage: %s points.json d k\n", argv[0]); return 2; }
  d = atoi(argv[2]); k = atoi(argv[3]);
  if (d < 1 || d > MAXD) { fprintf(stderr, "d out of range\n"); return 2; }
  FILE *f = fopen(argv[1], "r");
  if (!f) { perror(argv[1]); return 2; }
  long long v[MAXN * MAXD]; int nv = 0, ch, neg = 0, in = 0; long long cur = 0;
  while ((ch = fgetc(f)) != EOF) {
    if (isdigit(ch)) { cur = cur * 10 + (ch - '0'); in = 1; continue; }
    if (in) { if (nv == MAXN * MAXD) { fprintf(stderr, "too many numbers\n"); return 2; } v[nv++] = neg ? -cur : cur; }
    in = 0; cur = 0; neg = (ch == '-');
  }
  if (in) v[nv++] = neg ? -cur : cur;
  fclose(f);
  if (nv % d) { fprintf(stderr, "%d numbers is not a multiple of d=%d\n", nv, d); return 2; }
  n = nv / d;
  if (n > MAXN || k > n) { fprintf(stderr, "bad n=%d\n", n); return 2; }
  for (int i = 0; i < n; i++) for (int j = 0; j < d; j++) P[i][j] = v[i * d + j];

  sgn = calloc((size_t)1 << n, 1);
  int c[MAXD + 1], zero = 0, total = 0;
  for (int i = 0; i <= d; i++) c[i] = i;
  for (;;) {
    int s = orient(c), mask = 0;
    for (int i = 0; i <= d; i++) mask |= 1 << c[i];
    sgn[mask] = (signed char)s; total++;
    if (s == 0) zero++;
    int i = d;
    while (i >= 0 && c[i] == n - d - 1 + i) i--;
    if (i < 0) break;
    c[i]++;
    for (int j = i + 1; j <= d; j++) c[j] = c[j - 1] + 1;
  }
  printf("%d points in R^%d: %d of %d simplices degenerate\n", n, d, zero, total);
  long long gons = count_convex(k), control = count_convex(k - 1);
  printf("convex %d-subsets: %lld (control: convex %d-subsets: %lld)\n", k, gons, k - 1, control);
  int good = zero == 0 && gons == 0 && control > 0;
  printf("%s\n", good ? "OK: general position, no k-gon" : "FAILED");
  free(sgn);
  return good ? 0 : 1;
}

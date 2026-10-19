/* extend.c -- one-point extensions of a fixed point set in R^d with no k points in convex
   position.  If X has no k-gon, then X u {y} has a k-gon only among the sets T u {y} with T a
   convex (k-1)-subset of X, and whether T u {y} is in convex position depends only on the
   orientations det(1, y, x_S) for the d-subsets S of X, that is, on the cell of the
   hyperplane arrangement of X that contains y.  The program samples y from several
   distributions, computes those orientations once per sample, and counts the bad sets
   T u {y} by table lookup.  It prints every new best count and writes the best y found.

   Usage: extend points.txt d k seed seconds [out.txt]
   points.txt holds one point per line (integers or decimals). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#define MAXN 20
#define MAXD 6
static int n, d, k;
static double X[MAXN][MAXD];
static uint64_t rs = 88172645463325252ULL;
static uint64_t rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double unif(void){ return (rnd() >> 11) * (1.0/9007199254740992.0); }
static double gauss(void){ double u = unif() + 1e-300, v = unif(); return sqrt(-2*log(u))*cos(6.283185307179586*v); }

static double det(double A[MAXD+1][MAXD+1], int m){
  double M[MAXD+1][MAXD+1]; memcpy(M, A, sizeof(M)); double r = 1;
  for (int i = 0; i < m; i++){
    int piv = i; for (int j = i+1; j < m; j++) if (fabs(M[j][i]) > fabs(M[piv][i])) piv = j;
    if (M[piv][i] == 0) return 0;
    if (piv != i){ for (int c = 0; c < m; c++){ double t = M[i][c]; M[i][c] = M[piv][c]; M[piv][c] = t; } r = -r; }
    r *= M[i][i];
    for (int j = i+1; j < m; j++){ double f = M[j][i]/M[i][i]; for (int c = i; c < m; c++) M[j][c] -= f*M[i][c]; }
  }
  return r;
}
/* orientation of the points pts[0..d] (rows (1, x)) */
static double orient(double pts[MAXD+1][MAXD]){
  double A[MAXD+1][MAXD+1];
  for (int r = 0; r <= d; r++){ A[r][0] = 1; for (int c = 0; c < d; c++) A[r][c+1] = pts[r][c]; }
  return det(A, d+1);
}

/* sign table for (d+1)-subsets of the n+1 points X u {y}; y has index n.
   Indexed by the bitmask of the subset through a hash of sorted indices. */
static long C[MAXN+2][MAXN+2];
static long rank_of(const int *s, int m){ long r = 0; for (int i = 0; i < m; i++) r += C[s[i]][i+1]; return r; }
static signed char *sgn;     /* sgn[rank of sorted (d+1)-subset of 0..n] */
static double *val;

static void set_orient(const int *s){ /* s sorted, d+1 entries */
  double pts[MAXD+1][MAXD];
  for (int i = 0; i <= d; i++) for (int c = 0; c < d; c++) pts[i][c] = X[s[i]][c];
  double o = orient(pts); long r = rank_of(s, d+1);
  val[r] = o; sgn[r] = o > 0 ? 1 : (o < 0 ? -1 : 0);
}
/* sign of the orientation of the ordered tuple t[0..d] */
static int osign(const int *t){
  int s[MAXD+1]; memcpy(s, t, sizeof(int)*(d+1)); int par = 0;
  for (int i = 0; i <= d; i++) for (int j = 0; j < d - i; j++) if (s[j] > s[j+1]){ int u = s[j]; s[j] = s[j+1]; s[j+1] = u; par ^= 1; }
  int g = sgn[rank_of(s, d+1)];
  return par ? -g : g;
}
/* is p in the interior of the simplex s[0..d]? (all points distinct, general position) */
static int inside(int p, const int *s){
  int o = osign(s); if (!o) return 0;
  int t[MAXD+1];
  for (int j = 0; j <= d; j++){ memcpy(t, s, sizeof(int)*(d+1)); t[j] = p; if (osign(t) != o) return 0; }
  return 1;
}
/* convex position of the set S (m points)? */
static int convex(const int *S, int m){
  int rest[MAXN], sub[MAXD+1], c[MAXD+1];
  for (int a = 0; a < m; a++){
    int r = 0; for (int b = 0; b < m; b++) if (b != a) rest[r++] = S[b];
    for (int i = 0; i <= d; i++) c[i] = i;
    while (1){
      for (int i = 0; i <= d; i++) sub[i] = rest[c[i]];
      if (inside(S[a], sub)) return 0;
      int i = d; while (i >= 0 && c[i] == r - d - 1 + i) i--;
      if (i < 0) break;
      c[i]++; for (int j = i+1; j <= d; j++) c[j] = c[j-1] + 1;
    }
  }
  return 1;
}

static int ncv; static int (*cv)[MAXN];   /* convex (k-1)-subsets of X */
static double y[MAXD], best[MAXD];

static int count_bad(void){
  for (int c = 0; c < d; c++) X[n][c] = y[c];
  /* orientations of the (d+1)-subsets containing n */
  int s[MAXD+1], c2[MAXD];
  for (int i = 0; i < d; i++) c2[i] = i;
  double minabs = 1e300;
  while (1){
    for (int i = 0; i < d; i++) s[i] = c2[i];
    s[d] = n;
    set_orient(s);
    double a = fabs(val[rank_of(s, d+1)]); if (a < minabs) minabs = a;
    int i = d-1; while (i >= 0 && c2[i] == n - d + i) i--;
    if (i < 0) break;
    c2[i]++; for (int j = i+1; j < d; j++) c2[j] = c2[j-1] + 1;
  }
  if (minabs < 1e-9) return 1 << 30;   /* not in general position */
  int bad = 0, S[MAXN];
  for (int t = 0; t < ncv; t++){
    memcpy(S, cv[t], sizeof(int)*(k-1)); S[k-1] = n;
    if (convex(S, k)) bad++;
  }
  return bad;
}

int main(int argc, char **argv){
  if (argc < 6){ fprintf(stderr, "usage: extend points.txt d k seed seconds [out.txt]\n"); return 2; }
  d = atoi(argv[2]); k = atoi(argv[3]); rs ^= (uint64_t)atoll(argv[4]) * 0x9E3779B97F4A7C15ULL;
  for (int i = 0; i < 20; i++) rnd();
  double secs = atof(argv[5]);
  FILE *f = fopen(argv[1], "r"); if (!f){ perror(argv[1]); return 2; }
  n = 0; double v;
  for (int i = 0; fscanf(f, "%lf", &v) == 1; i++){ X[i/d][i%d] = v; n = i/d + 1; }
  fclose(f);
  if (n + 1 > MAXN || d > MAXD){ fprintf(stderr, "too large\n"); return 2; }
  for (int i = 0; i <= MAXN+1; i++){ C[i][0] = 1; for (int j = 1; j <= i && j <= MAXN+1; j++) C[i][j] = C[i-1][j-1] + (j < i ? C[i-1][j] : 0); }
  long ns = C[n+1][d+1];
  sgn = calloc(ns, 1); val = calloc(ns, sizeof(double));
  { int s[MAXD+1]; for (int i = 0; i <= d; i++) s[i] = i;
    while (1){ set_orient(s);
      int i = d; while (i >= 0 && s[i] == n - 1 - d + i) i--;
      if (i < 0) break;
      s[i]++; for (int j = i+1; j <= d; j++) s[j] = s[j-1] + 1; } }
  /* convex (k-1)-subsets of X (and check: no convex k-subset) */
  cv = malloc(sizeof(*cv) * C[n][k-1]); ncv = 0;
  { int s[MAXN]; for (int i = 0; i < k-1; i++) s[i] = i;
    while (1){ if (convex(s, k-1)){ memcpy(cv[ncv], s, sizeof(int)*(k-1)); ncv++; }
      int i = k-2; while (i >= 0 && s[i] == n - k + 1 + i) i--;
      if (i < 0) break;
      s[i]++; for (int j = i+1; j < k-1; j++) s[j] = s[j-1] + 1; } }
  double lo[MAXD], hi[MAXD];
  for (int c = 0; c < d; c++){ lo[c] = hi[c] = X[0][c]; for (int i = 1; i < n; i++){ if (X[i][c] < lo[c]) lo[c] = X[i][c]; if (X[i][c] > hi[c]) hi[c] = X[i][c]; } }
  double span = 0; for (int c = 0; c < d; c++) if (hi[c] - lo[c] > span) span = hi[c] - lo[c];
  fprintf(stderr, "%d points in R^%d, %d convex %d-subsets\n", n, d, ncv, k-1);
  int bestbad = 1 << 30; long samples = 0; time_t t0 = time(NULL);
  int hist[64] = {0};
  while (difftime(time(NULL), t0) < secs){
    int mode = rnd() % 4;
    if (mode == 0){ /* near a point */
      int i = rnd() % n; double sc = span * pow(10, -3 * unif());
      for (int c = 0; c < d; c++) y[c] = X[i][c] + sc * gauss();
    } else if (mode == 1){ /* random convex combination */
      double w[MAXN], sw = 0; for (int i = 0; i < n; i++){ w[i] = -log(unif() + 1e-300); w[i] = pow(w[i], 1 + 4*unif()); sw += w[i]; }
      for (int c = 0; c < d; c++){ y[c] = 0; for (int i = 0; i < n; i++) y[c] += w[i]*X[i][c]/sw; }
    } else if (mode == 2){ /* box, possibly enlarged */
      double e = pow(10, 2*unif()) - 1;
      for (int c = 0; c < d; c++) y[c] = lo[c] - e*span/2 + (hi[c] - lo[c] + e*span) * unif();
    } else { /* near the best so far */
      if (bestbad == 1 << 30){ samples++; continue; }
      double sc = span * pow(10, -4 * unif());
      for (int c = 0; c < d; c++) y[c] = best[c] + sc * gauss();
    }
    int b = count_bad(); samples++;
    if (b < 64) hist[b]++;
    if (b < bestbad){
      bestbad = b; memcpy(best, y, sizeof(best));
      printf("samples %ld best %d  y =", samples, b); for (int c = 0; c < d; c++) printf(" %.10g", y[c]); printf("\n"); fflush(stdout);
      if (argc > 6){ FILE *g = fopen(argv[6], "w"); for (int c = 0; c < d; c++) fprintf(g, "%.17g ", y[c]); fprintf(g, "\n"); fclose(g); }
      if (b == 0) break;
    }
  }
  printf("samples %ld, best %d; histogram of bad counts:", samples, bestbad);
  for (int b = 0; b < 64; b++) if (hist[b]) printf(" %d:%d", b, hist[b]);
  printf("\n");
  return bestbad == 0 ? 0 : 1;
}

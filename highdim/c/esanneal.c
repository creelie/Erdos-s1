/* Simulated annealing for point sets in R^d (general position) with no k
   points in convex position.  Cost = number of k-subsets in convex position.
   A subset is in convex position iff no member lies in the simplex spanned by
   d+1 other members (Caratheodory).  Convex position is hereditary, so convex
   k-sets are enumerated by DFS over convex sets.
   Usage: esanneal n d k seed seconds [init.txt]   (writes best_<n>_<d>_<k>_<seed>.txt) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#define MAXN 40
#define MAXD 8
typedef uint64_t u64;
static int n, d, k;
static double P[MAXN][MAXD];
static u64 *W[MAXN]; static int Wn[MAXN], Wcap[MAXN];
static u64 rs = 88172645463325252ULL;
static u64 rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double unif(void){ return (rnd() >> 11) * (1.0/9007199254740992.0); }
static double gauss(void){ double u=unif()+1e-300, v=unif(); return sqrt(-2*log(u))*cos(2*M_PI*v); }

static double det(double A[MAXD+1][MAXD+1], int m){
  double M[MAXD+1][MAXD+1]; memcpy(M, A, sizeof(M)); double r = 1;
  for (int i=0;i<m;i++){
    int piv=i; for(int j=i+1;j<m;j++) if (fabs(M[j][i])>fabs(M[piv][i])) piv=j;
    if (M[piv][i]==0) return 0;
    if (piv!=i){ for(int c=0;c<m;c++){double t=M[i][c];M[i][c]=M[piv][c];M[piv][c]=t;} r=-r; }
    r *= M[i][i];
    for (int j=i+1;j<m;j++){ double f=M[j][i]/M[i][i]; for(int c=i;c<m;c++) M[j][c]-=f*M[i][c]; }
  }
  return r;
}
static double orient(const int *idx){
  double A[MAXD+1][MAXD+1];
  for (int c=0;c<=d;c++){ A[0][c]=1; for(int r=0;r<d;r++) A[r+1][c]=P[idx[c]][r]; }
  return det(A, d+1);
}
static void addW(int p, u64 m){
  if (Wn[p]==Wcap[p]){ Wcap[p]=Wcap[p]?2*Wcap[p]:256; W[p]=realloc(W[p], Wcap[p]*sizeof(u64)); }
  W[p][Wn[p]++]=m;
}
static double minabs;   /* smallest |orientation|, to keep general position */
static void build(void){
  for (int p=0;p<n;p++) Wn[p]=0;
  minabs = 1e300;
  int idx[MAXD+1]; for (int i=0;i<=d;i++) idx[i]=i;
  while (1){
    double o = orient(idx); if (fabs(o)<minabs) minabs=fabs(o);
    u64 m=0; for(int i=0;i<=d;i++) m|=1ULL<<idx[i];
    for (int p=0;p<n;p++){
      if (m>>p&1) continue;
      int ok=1, t[MAXD+1];
      for (int j=0;j<=d && ok;j++){ memcpy(t,idx,sizeof(int)*(d+1)); t[j]=p; double o2=orient(t); if ((o2>0)!=(o>0)) ok=0; }
      if (ok) addW(p, m);
    }
    int i=d; while (i>=0 && idx[i]==n-d-1+i) i--; if (i<0) break;
    idx[i]++; for (int j=i+1;j<=d;j++) idx[j]=idx[j-1]+1;
  }
}
static long cnt;
static int convex_with(u64 T, int q){ /* T convex, q not in T: is T+q convex? */
  u64 U = T | (1ULL<<q);
  for (int i=0;i<Wn[q];i++) if ((W[q][i] & T)==W[q][i]) return 0;
  for (u64 s=T; s; s&=s-1){ int p=__builtin_ctzll(s);
    for (int i=0;i<Wn[p];i++){ u64 m=W[p][i]; if ((m>>q&1) && (m & U)==m) return 0; } }
  return 1;
}
static void dfs(u64 T, int size, int start){
  if (size==k){ cnt++; return; }
  for (int q=start;q<n;q++) if (n-q >= k-size && convex_with(T,q)) dfs(T|(1ULL<<q), size+1, q+1);
}
static long cost(void){ build(); cnt=0; dfs(0,0,0); return cnt; }
static void save(const char *fn){ FILE *f=fopen(fn,"w"); for(int i=0;i<n;i++){ for(int j=0;j<d;j++) fprintf(f,"%.17g ",P[i][j]); fprintf(f,"\n"); } fclose(f); }

int main(int argc, char **argv){
  n=atoi(argv[1]); d=atoi(argv[2]); k=atoi(argv[3]); int seed=atoi(argv[4]); double secs=atof(argv[5]);
  rs ^= (u64)seed*0x9E3779B97F4A7C15ULL; for(int i=0;i<10;i++) rnd();
  int m0=0;
  if (argc>6){ FILE *f=fopen(argv[6],"r"); while (m0<n && fscanf(f,"%lf",&P[m0][0])==1){ for(int j=1;j<d;j++) fscanf(f,"%lf",&P[m0][j]); m0++; } fclose(f); }
  for (int i=m0;i<n;i++) for(int j=0;j<d;j++) P[i][j] = m0? P[rnd()%m0][j]+0.05*gauss() : gauss();
  char fn[128]; sprintf(fn,"best_%d_%d_%d_%d.txt",n,d,k,seed);
  long c=cost(), best=c; double Best[MAXN][MAXD]; memcpy(Best,P,sizeof(P));
  double T=2.0; long it=0, lastimp=0; time_t t0=time(NULL);
  printf("start cost %ld\n", c); fflush(stdout);
  static const double sc[4]={0.01,0.03,0.1,0.3};
  while (best>0 && difftime(time(NULL),t0)<secs){
    it++;
    double old[MAXD]; int i=rnd()%n; memcpy(old,P[i],sizeof(old));
    double r=unif();
    if (r<0.8){ double s=sc[rnd()%4]; for(int j=0;j<d;j++) P[i][j]+=s*gauss(); }
    else if (r<0.95){ int o=rnd()%n; for(int j=0;j<d;j++) P[i][j]=P[o][j]+0.2*gauss(); }
    else { double s=(double[]){0.5,1,3}[rnd()%3]; for(int j=0;j<d;j++) P[i][j]=s*gauss(); }
    long c2=cost();
    if (minabs<1e-9){ memcpy(P[i],old,sizeof(old)); continue; }
    if (c2<=c || unif()<exp((c-c2)/T)) c=c2; else memcpy(P[i],old,sizeof(old));
    if (c<best){ best=c; memcpy(Best,P,sizeof(P)); lastimp=it; save(fn);
      printf("%.0fs it %ld best %ld\n", difftime(time(NULL),t0), it, best); fflush(stdout); }
    T = T*0.999; if (T<0.1) T=0.1;
    if (it-lastimp>4000){ memcpy(P,Best,sizeof(P)); c=best; T=1.5; lastimp=it; }
  }
  memcpy(P,Best,sizeof(P)); save(fn);
  printf(best==0? "FOUND\n" : "# end best %ld\n", best); printf("# it %ld\n", it);
  return 0;
}

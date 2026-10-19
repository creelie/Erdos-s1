/* Make a point set with no k points in convex position more robust.
   For a k-subset T, r(T) is the largest distance from a point p of T to the
   boundary of a simplex spanned by d+1 other points of T that contains p
   (0 if there is none).  The margin is min over T of r(T), divided by the
   diameter of the set.  A positive margin certifies that no k points are in
   convex position, and the larger it is, the coarser the integer grid the set
   can be rounded to.  Hill climbing on the margin.
   Usage: polish n d k seconds in.txt out.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#define MAXN 24
#define MAXD 6
typedef uint64_t u64;
static int n, d, k;
static double P[MAXN][MAXD];
static double *DM[MAXN];
static long Cb[MAXN+1][MAXN+1];
static u64 rs = 0x2545F4914F6CDD1DULL;
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
static double facetnorm(const int *idx, int skip){  /* sqrt of the Gram determinant of the facet's edge vectors */
  int v[MAXD+1], t=0; for (int i=0;i<=d;i++) if (i!=skip) v[t++]=idx[i];
  double E[MAXD][MAXD], G[MAXD+1][MAXD+1];
  for (int a=1;a<d;a++) for (int c=0;c<d;c++) E[a-1][c]=P[v[a]][c]-P[v[0]][c];
  for (int a=0;a<d-1;a++) for (int b=0;b<d-1;b++){ double s=0; for (int c=0;c<d;c++) s+=E[a][c]*E[b][c]; G[a][b]=s; }
  double g = d>1 ? det(G, d-1) : 1; return sqrt(fabs(g));
}
static long rank(u64 m){ long r=0; int i=0; for (u64 s=m; s; s&=s-1){ int p=__builtin_ctzll(s); r+=Cb[p][i+1]; i++; } return r; }
static double margin(void){
  double D=0; for (int a=0;a<n;a++) for (int b=a+1;b<n;b++){ double s=0; for (int c=0;c<d;c++) s+=(P[a][c]-P[b][c])*(P[a][c]-P[b][c]); if (s>D) D=s; }
  D=sqrt(D);
  int idx[MAXD+1]; for (int i=0;i<=d;i++) idx[i]=i;
  while (1){
    u64 m=0; for (int i=0;i<=d;i++) m|=1ULL<<idx[i]; long r=rank(m);
    double o=orient(idx), F[MAXD+1];
    for (int j=0;j<=d;j++) F[j]=facetnorm(idx,j);
    for (int p=0;p<n;p++){ if (m>>p&1) continue; DM[p][r]=0;
      if (fabs(o)<1e-300) continue;
      double dist=1e300; int ok=1, t[MAXD+1];
      for (int j=0;j<=d && ok;j++){ memcpy(t,idx,sizeof(int)*(d+1)); t[j]=p; double oj=orient(t);
        if ((oj>0)!=(o>0) || oj==0) ok=0; else { double h=fabs(oj)/F[j]; if (h<dist) dist=h; } }
      if (ok) DM[p][r]=dist/D;
    }
    int i=d; while (i>=0 && idx[i]==n-d-1+i) i--; if (i<0) break;
    idx[i]++; for (int j=i+1;j<=d;j++) idx[j]=idx[j-1]+1;
  }
  double worst=1e300; int T[MAXN]; for (int i=0;i<k;i++) T[i]=i;
  while (1){
    double best=0;
    for (int a=0;a<k;a++){ int p=T[a], b[MAXN], r=0; for (int c=0;c<k;c++) if (c!=a) b[r++]=T[c];
      int c[MAXD+1]; for (int i=0;i<=d;i++) c[i]=i;
      while (1){ u64 m=0; for (int i=0;i<=d;i++) m|=1ULL<<b[c[i]]; double v=DM[p][rank(m)]; if (v>best) best=v;
        int i=d; while (i>=0 && c[i]==r-d-1+i) i--; if (i<0) break; c[i]++; for (int j=i+1;j<=d;j++) c[j]=c[j-1]+1; }
    }
    if (best<worst) worst=best;
    if (worst<=0) return 0;
    int i=k-1; while (i>=0 && T[i]==n-k+i) i--; if (i<0) break;
    T[i]++; for (int j=i+1;j<k;j++) T[j]=T[j-1]+1;
  }
  return worst;
}
static void save(const char *fn){ FILE *f=fopen(fn,"w"); for(int i=0;i<n;i++){ for(int j=0;j<d;j++) fprintf(f,"%.17g ",P[i][j]); fprintf(f,"\n"); } fclose(f); }
int main(int argc, char **argv){
  if (argc<7){ fprintf(stderr,"usage: polish n d k seconds in.txt out.txt\n"); return 2; }
  n=atoi(argv[1]); d=atoi(argv[2]); k=atoi(argv[3]); double secs=atof(argv[4]);
  for (int i=0;i<=MAXN;i++){ Cb[i][0]=1; for (int j=1;j<=MAXN;j++) Cb[i][j]= i? Cb[i-1][j-1]+Cb[i-1][j] : 0; }
  for (int p=0;p<n;p++) DM[p]=calloc(Cb[n][d+1],sizeof(double));
  FILE *f=fopen(argv[5],"r"); if (!f){ perror(argv[5]); return 2; }
  for (int i=0;i<n;i++) for (int j=0;j<d;j++) if (fscanf(f,"%lf",&P[i][j])!=1){ fprintf(stderr,"short input\n"); return 2; }
  fclose(f);
  double D=0; for (int a=0;a<n;a++) for (int b=a+1;b<n;b++){ double s=0; for (int c=0;c<d;c++) s+=(P[a][c]-P[b][c])*(P[a][c]-P[b][c]); if (s>D) D=s; } D=sqrt(D);
  double cur=margin(); printf("start margin %.3e\n", cur); fflush(stdout);
  if (cur<=0){ printf("input has a convex %d-subset\n",k); return 1; }
  double step=0.02; long it=0, since=0; time_t t0=time(NULL);
  while (difftime(time(NULL),t0)<secs){
    it++; int i=rnd()%n; double old[MAXD]; memcpy(old,P[i],sizeof(old));
    double s=step*D*(unif()<0.5?1:0.2);
    for (int j=0;j<d;j++) P[i][j]+=s*gauss();
    double m=margin();
    if (m>=cur){ if (m>cur*1.0001) since=0; cur=m; } else { memcpy(P[i],old,sizeof(old)); since++; }
    if (since>3000){ step*=0.5; since=0; if (step<1e-5) step=0.02; }
    if (it%2000==0){ printf("%.0fs it %ld margin %.3e step %.1e\n", difftime(time(NULL),t0), it, cur, step); fflush(stdout); save(argv[6]); }
  }
  save(argv[6]); printf("# end margin %.3e\n", cur); return 0;
}

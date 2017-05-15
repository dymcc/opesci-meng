/* Generated from acoustic_skew_tile_parameterized.cloog by CLooG 0.18.4-b2b4f77 gmp bits in 0.02s. */
/* DON'T FORGET TO USE -lm OPTION TO COMPILE. */

/* Useful headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Parameter value. */
#define PARVAL1 1
#define PARVAL2 1
#define PARVAL3 1
#define PARVAL4 1
/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#ifdef TIME 
#define IF_TIME(foo) foo; 
#else
#define IF_TIME(foo)
#endif

/* Statement macros (please set). */
#define S1(i,j,k,l,m,n,o) {total++; printf("S1 %d %d %d %d %d %d %d\n",i,j,k,l,m,n,o);}

int main() {
  /* Scattering iterators. */
  int t, xx, yy, zz, x, y, z;
  /* Original iterators. */
  int i, j, k, l, m, n, o;
  /* Parameters. */
  int M=PARVAL1, N=PARVAL2, O=PARVAL3, P=PARVAL4;
  int total=0;

  int t, xx, yy, zz, x, y, z;
  for (t=0;t<=82;t++) {
    for (xx=0;xx<=5;xx++) {
      for (yy=0;yy<=5;yy++) {
        for (zz=0;zz<=5;zz++) {
          for (x=max(4*t+4,4*t+8*xx);x<=4*t+8*xx+7;x++) {
            for (y=max(4*t+4,4*t+8*yy);y<=4*t+8*yy+7;y++) {
              for (z=max(4*t+4,4*t+8*zz);z<=4*t+8*zz+7;z++) {
                S1(t,xx,yy,zz,(-4*t+x),(-4*t+y),(-4*t+z));
              }
            }
          }
        }
      }
    }
  }

  printf("Number of integral points: %d.\n",total);
  return 0;
}

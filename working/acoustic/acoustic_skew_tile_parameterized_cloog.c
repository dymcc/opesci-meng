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
  int i4, xx, yy, zz, i1, i2, i3;
  /* Original iterators. */
  int i, j, k, l, m, n, o;
  /* Parameters. */
  int M=PARVAL1, N=PARVAL2, O=PARVAL3, P=PARVAL4;
  int total=0;

  for (i4=0;i4<=82;i4++) {
    for (xx=1;xx<=13;xx++) {
      for (yy=1;yy<=13;yy++) {
        for (zz=1;zz<=13;zz++) {
          for (i1=i4+2*xx;i1<=i4+2*xx+1;i1++) {
            for (i2=i4+2*yy;i2<=i4+2*yy+1;i2++) {
              for (i3=i4+2*zz;i3<=i4+2*zz+1;i3++) {
                S1(i4,xx,yy,zz,(-i4+i1),(-i4+i2),(-i4+i3));
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

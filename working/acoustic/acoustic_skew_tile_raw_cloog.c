/* Generated from acoustic_skew_tile.cloog by CLooG 0.18.4-b2b4f77 gmp bits in 0.02s. */
/* DON'T FORGET TO USE -lm OPTION TO COMPILE. */

/* Useful headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Parameter value. */
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
  int ii1, ii2, ii3, i4, i1, i2, i3;
  /* Original iterators. */
  int i, j, k, l, m, n, o;
  int total=0;

  for (ii1=1;ii1<=13;ii1++) {
    for (ii2=1;ii2<=13;ii2++) {
      for (ii3=1;ii3<=13;ii3++) {
        for (i4=0;i4<=82;i4++) {
          for (i1=2*ii1+i4;i1<=2*ii1+i4+1;i1++) {
            for (i2=2*ii2+i4;i2<=2*ii2+i4+1;i2++) {
              for (i3=2*ii3+i4;i3<=2*ii3+i4+1;i3++) {
                S1(ii1,ii2,ii3,i4,(-i4+i1),(-i4+i2),(-i4+i3));
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

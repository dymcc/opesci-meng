/* Generated from acoustic_skew_tile.cloog by CLooG 0.18.4-b2b4f77 gmp bits in 0.01s. */
/* DON'T FORGET TO USE -lm OPTION TO COMPILE. */

/* Useful headers. */
#include <string.h>
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
  int a[200][80][80];
  memset(a, 0, sizeof a);
  int total=0;

  for (ii1=0;ii1<=9;ii1++) {
      for (ii2=0;ii2<=9;ii2++) {
          for (i4=0;i4<=192;i4++) {
              for (i1=1*ii1+i4;i1<=min(i4+75,8*ii1+i4+7);i1++) {
                  for (i2=1*ii2+i4;i2<=min(i4+75,8*ii2+i4+7);i2++) {
                      /*S1(ii1,ii2,i4,(-i4+i1),(-i4+i2));*/
                      a[i4+2][i1-i4][i2-i4] = a[i4+1][i1-i4][i2-i4] + a[i4+1][i1-i4][i2-i4] + a[i4+1][i1-i4][i2-i4 - 1] + a[i4+1][i1-i4][i2-i4 + 1] + a[i4+1][i1-i4 - 1][i2-i4] + a[i4+1][i1-i4 + 1][i2-i4] - a[i4+1][i1-i4][i2-i4] + a[i4][i1-i4][i2-i4] + a[i4+1][i1-i4][i2-i4];
                  }
              }
          }
      }
  }
  return 0;
}

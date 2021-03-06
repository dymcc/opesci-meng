/* Generated from tile_simplestencil.cloog by CLooG 0.18.4-b2b4f77 gmp bits in 0.00s. */
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
#define S1(i,j,k,l) {total++; printf("S1 %d %d %d %d \n",i,j,k,l);}

int main() {
  /* Original iterators. */
  int i, j, k, l;
  int total=0;

  for (i=0;i<=49;i++) {
    for (j=0;j<=49;j++) {
      for (k=0;k<=49;k++) {
        for (l=0;l<=49;l++) {
                S1(i,j,k,l);
        }
      }
    }
  }

  printf("Number of integral points: %d.\n",total);
  return 0;
}

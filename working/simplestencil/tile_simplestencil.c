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
#define S1(i,j,k,l,m,n,o) {total++; printf("S1 %d %d %d %d %d %d %d\n",i,j,k,l,m,n,o);}

int main() {
  /* Original iterators. */
  int i, j, k, l, m, n, o;
  int total=0;

  for (i=0;i<=12;i++) {
    for (j=0;j<=12;j++) {
      for (k=0;k<=12;k++) {
        for (l=0;l<=49;l++) {
          for (m=4*i;m<=min(49,4*i+3);m++) {
            for (n=4*j;n<=min(49,4*j+3);n++) {
              for (o=4*k;o<=min(49,4*k+3);o++) {
                S1(i,j,k,l,m,n,o);
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

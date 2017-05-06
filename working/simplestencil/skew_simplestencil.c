/* Generated from skew_simplestencil.cloog by CLooG 0.18.4-b2b4f77 gmp bits in 0.00s. */
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
#define S1(i,j,k,l) {total++; printf("S1 %d %d %d %d\n",i,j,k,l);}

int main() {
  /* Scattering iterators. */
  int t1, t2, t3, t4;
  /* Original iterators. */
  int i, j, k, l;
  int total=0;

  for (t1=0;t1<=49;t1++) {
    for (t2=2*t1;t2<=2*t1+49;t2++) {
      for (t3=2*t1;t3<=2*t1+49;t3++) {
        for (t4=2*t1;t4<=2*t1+49;t4++) {
          S1(t1,(-2*t1+t2),(-2*t1+t3),(-2*t1+t4));
        }
      }
    }
  }

  printf("Number of integral points: %d.\n",total);
  return 0;
}

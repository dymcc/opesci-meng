/* Generated from skew_tile_simplestencil.cloog by CLooG 0.18.4-b2b4f77 gmp bits in 0.02s. */
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
  int xx, yy, zz, t, x, y, z;
  /* Original iterators. */
  int i, j, k, l, m, n, o;
  int total=0;

  for (xx=0;xx<=12;xx++) {
    for (yy=0;yy<=12;yy++) {
      for (zz=0;zz<=12;zz++) {
        for (t=0;t<=49;t++) {
          for (x=4*xx+2*t;x<=min(2*t+49,4*xx+2*t+3);x++) {
            for (y=4*yy+2*t;y<=min(2*t+49,4*yy+2*t+3);y++) {
              for (z=4*zz+2*t;z<=min(2*t+49,4*zz+2*t+3);z++) {
                S1(xx,yy,zz,t,(-2*t+x),(-2*t+y),(-2*t+z));
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

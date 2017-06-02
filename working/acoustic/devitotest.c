/* Generated from cloog_inputs/grid_512/acoustic_skew4_tile16.cloog by CLooG 0.18.5-56daab3 gmp bits in 0.02s. */
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
#define S1(i,j,k,l,m,n,o) {total++;}//printf("S1 %d %d %d %d %d %d %d\n",i,j,k,l,m,n,o);}

int main() {
  /* Scattering iterators. */
  int time, xx, yy, zz, x, y, z;
  /* Original iterators. */
  int i, j, k, l, m, n, o;
  /* Parameters. */
  int M=PARVAL1, N=PARVAL2, O=PARVAL3, P=PARVAL4;
  int total=0;


  for (time=1;time<=31;time++) {
      for (int x = 4; x < 532 - 4; x += 1)
      {
          for (int y = 4; y < 532 - 4; y += 1)
          {
              for (int z = 4; z < 532 - 4; z += 1)
                  S1(time,xx,yy,zz,(-4*time+x),(-4*time+y),(-4*time+z));
          }
      }
  }

  printf("Number of integral points: %d.\n",total);
  return 0;
}

/* Generated from ACTUALLYSKEWEDTILING.cloog by CLooG 0.18.5-56daab3 gmp bits in 0.04s. */
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
  int xx, yy, zz, time, x, y, z;
  /* Original iterators. */
  int i, j, k, l, m, n, o;
  /* Parameters. */
  int time_size=PARVAL1, ub=PARVAL2, lb=PARVAL3, ts=PARVAL4;
  int total=0;

  for (xx=0;xx<=floord(time_size+65,4);xx++) {
    for (yy=max(0,xx-17);yy<=min(floord(time_size+65,4),xx+17);yy++) {
      for (zz=max(max(0,xx-17),yy-17);zz<=min(min(floord(time_size+65,4),xx+17),yy+17);zz++) {
        for (time=max(max(max(1,4*xx-67),4*yy-67),4*zz-67);time<=min(min(min(time_size-2,4*xx+2),4*yy+2),4*zz+2);time++) {
          for (x=max(16*xx,4*time+4);x<=min(4*time+271,16*xx+15);x++) {
            for (y=max(16*yy,4*time+4);y<=min(4*time+271,16*yy+15);y++) {
              for (z=max(16*zz,4*time+4);z<=min(4*time+271,16*zz+15);z++) {
                S1(time,xx,yy,zz,x,y,z);
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

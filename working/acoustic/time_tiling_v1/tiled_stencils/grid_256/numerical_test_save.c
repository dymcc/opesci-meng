#define _POSIX_C_SOURCE 200809L
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include "xmmintrin.h"
#include "pmmintrin.h"

#define DIM 276
#define TIME 10

#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

static float u_grid[TIME][DIM][DIM][DIM] __attribute__((aligned(64)));
static float v_grid[TIME][DIM][DIM][DIM] __attribute__((aligned(64)));

void stencil(float u[TIME][DIM][DIM][DIM], int time, int x, int y, int z){
    int t0 = time;
    int t1 = time + 1;
    int t2 = time - 1;
    u[t1][x][y][z] = 1e-1F +
                     (3.04e-1F*u[t2][x][y][z]
                   - 8.25e-5F*(u[t0][x][y][z - 4] +
                               u[t0][x][y][z + 4] +
                               u[t0][x][y - 4][z] +
                               u[t0][x][y + 4][z] +
                               u[t0][x - 4][y][z] +
                               u[t0][x + 4][y][z])
                   + 1.17e-3F*(u[t0][x][y][z - 3] +
                               u[t0][x][y][z + 3] +
                               u[t0][x][y - 3][z] +
                               u[t0][x][y + 3][z] +
                               u[t0][x - 3][y][z] +
                               u[t0][x + 3][y][z])
                   - 9.24e-3F*(u[t0][x][y][z - 2] +
                               u[t0][x][y][z + 2] +
                               u[t0][x][y - 2][z] +
                               u[t0][x][y + 2][z] +
                               u[t0][x - 2][y][z] +
                               u[t0][x + 2][y][z])
                + 7.39e-2F*(u[t0][x][y][z - 1] +
                               u[t0][x][y][z + 1] +
                               u[t0][x][y - 1][z] +
                               u[t0][x][y + 1][z] +
                               u[t0][x - 1][y][z] +
                               u[t0][x + 1][y][z])
                    + 4.11e-3F*u[t0][x][y][z]
                    - 3.94e-1F*u[t0][x][y][z])/(3.04e-1F);
}

void forward_skewed_tiled_time(float u[TIME][DIM][DIM][DIM], const int time_size, const int x_size, const int y_size, const int z_size){
    for (int xx=0;xx<=floord(time_size+65,4);xx++) {
        for (int yy=max(0,xx-17);yy<=min(floord(time_size+65,4),xx+17);yy++) {
            for (int zz=max(max(0,xx-17),yy-17);zz<=min(min(floord(time_size+65,4),xx+17),yy+17);zz++) {
                #pragma omp parallel
                {
                    #pragma omp for schedule(static)
                    for (int time=max(max(max(1,4*xx-67),4*yy-67),4*zz-67);time<=min(min(min(time_size-2,4*xx+2),4*yy+2),4*zz+2);time++) {
                        int skew = 4*time;
                        for (int x=max(16*xx,4*time+4);x<=min(4*time+271,16*xx+15);x++) {
                            for (int y=max(16*yy,4*time+4);y<=min(4*time+271,16*yy+15);y++) {
                            #pragma ivdep
                            #pragma omp simd
                                for (int z=max(16*zz,4*time+4);z<=min(4*time+271,16*zz+15);z++) {
                                stencil(u, time, x-skew, y-skew, z-skew);
                            }
                        }
                    }
                }
            }
        }
    }
    }
}

void forward_skewed_tiled(float u[TIME][DIM][DIM][DIM], const int time_size, const int x_size, const int y_size, const int z_size){
    for (int time = 1; time < time_size - 1; time += 1){
    int skew = 4*time;
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for (int xx=ceild(time-2,4);xx<=floord(time+67,4);xx++) {
            for (int yy=ceild(time-2,4);yy<=floord(time+67,4);yy++) {
                for (int zz=ceild(time-2,4);zz<=floord(time+67,4);zz++) {
                    for (int x=max(16*xx,4*time+4);x<=min(4*time+271,16*xx+15);x++) {
                        for (int y=max(16*yy,4*time+4);y<=min(4*time+271,16*yy+15);y++) {
                            #pragma ivdep
                            #pragma omp simd
                            for (int z=max(16*zz,4*time+4);z<=min(4*time+271,16*zz+15);z++) {
                                stencil(u, time, x-skew, y-skew, z-skew);
                            }
                        }
                    }
                }
            }
        }
    }
    }
}

void forward_control(float u[TIME][DIM][DIM][DIM], const int time_size, const int x_size, const int y_size, const int z_size)
{
  for (int time = 1; time < time_size - 1; time += 1)
  {
    #pragma omp parallel
    {
      #pragma omp for schedule(static)
      for (int x = 4; x < x_size - 4; x += 1)
      {
        for (int y = 4; y < y_size - 4; y += 1)
        {
          #pragma ivdep
          #pragma omp simd
          for (int z = 4; z < z_size - 4; z += 1)
          {
              stencil(u, time, x, y, z);
          }
        }
      }
    }
    }
}

int compare_output(const int time_size, const int x_size, const int y_size, const int z_size){
    int equal = 0;
    for (int time = 1; time < TIME; time += 1)
    {
        {
            #pragma omp for schedule(static)
            for (int x = 0; x < x_size - 0; x += 1)
            {
                for (int y = 0; y < y_size - 0; y += 1)
                {
                    for (int z = 0; z < z_size - 0; z += 1)
                    {
                        equal += fabs(u_grid[time][x][y][z] - v_grid[time][x][y][z]) >= 10e-6;
                    }
                }
            }
        }
    }
    return equal;
}

int main(){
    /*v_grid[0][4][4][4] = 1.013e-3F;*/
    /*v_grid[0][4][4][5] = 1.013e-3F;*/
    /*v_grid[0][4][7][5] = 1.013e-3F;*/
    forward_control(u_grid, TIME, DIM, DIM, DIM);
    forward_skewed_tiled_time(v_grid, TIME, DIM, DIM, DIM);
    printf("%d", compare_output(TIME, DIM, DIM, DIM));
    /*forward_tiled(v, TIME, DIM, DIM, DIM);*/
}


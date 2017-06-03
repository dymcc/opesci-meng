#define _POSIX_C_SOURCE 200809L
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include "xmmintrin.h"
#include "pmmintrin.h"

#define DIM 276
#define TIME 10
#define TIME_SIZE 83

#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

static float u_grid[TIME][DIM][DIM][DIM] __attribute__((aligned(64)));
static float v_grid[TIME][DIM][DIM][DIM] __attribute__((aligned(64)));

void forward_skewed_tiled(float u[TIME][DIM][DIM][DIM], const int time_size, const int x_size, const int y_size, const int z_size){
    for (int time = 1; time < time_size - 1; time += 1){
    int skew = 4*time;
    int t0 = (time) % 8;
    int t1 = (time + 1) % 8;
    int t2 = (time - 1) % 8;
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for(int xx=0;xx<=11;xx++) {
            for(int yy=0;yy<=11;yy++) {
                for(int zz=0;zz<=11;zz++) {
                        for(int x=max(4*time+4,24*xx+4*time);x<=min(4*time+271,24*xx+4*time+23);x++) {
                            for(int y=max(4*time+4,24*yy+4*time);y<=min(4*time+271,24*yy+4*time+23);y++) {
                                #pragma ivdep
                                #pragma omp simd
                                for(int z=max(4*time+4,24*zz+4*time);z<=min(4*time+271,24*zz+4*time+23);z++) {
                                    u[t1][x-skew][y-skew][z-skew] = 1e-1F +
                                                     (3.04e-1F*u[t2][x-skew][y-skew][z-skew]
                                                   - 8.25e-5F*(u[t0][x-skew][y-skew][z-skew - 4] +
                                                               u[t0][x-skew][y-skew][z-skew + 4] +
                                                               u[t0][x-skew][y-skew - 4][z-skew] +
                                                               u[t0][x-skew][y-skew + 4][z-skew] +
                                                               u[t0][x-skew - 4][y-skew][z-skew] +
                                                               u[t0][x-skew + 4][y-skew][z-skew])
                                                   + 1.17e-3F*(u[t0][x-skew][y-skew][z-skew - 3] +
                                                               u[t0][x-skew][y-skew][z-skew + 3] +
                                                               u[t0][x-skew][y-skew - 3][z-skew] +
                                                               u[t0][x-skew][y-skew + 3][z-skew] +
                                                               u[t0][x-skew - 3][y-skew][z-skew] +
                                                               u[t0][x-skew + 3][y-skew][z-skew])
                                                   - 9.24e-3F*(u[t0][x-skew][y-skew][z-skew - 2] +
                                                               u[t0][x-skew][y-skew][z-skew + 2] +
                                                               u[t0][x-skew][y-skew - 2][z-skew] +
                                                               u[t0][x-skew][y-skew + 2][z-skew] +
                                                               u[t0][x-skew - 2][y-skew][z-skew] +
                                                               u[t0][x-skew + 2][y-skew][z-skew])
                                                + 7.39e-2F*(u[t0][x-skew][y-skew][z-skew - 1] +
                                                               u[t0][x-skew][y-skew][z-skew + 1] +
                                                               u[t0][x-skew][y-skew - 1][z-skew] +
                                                               u[t0][x-skew][y-skew + 1][z-skew] +
                                                               u[t0][x-skew - 1][y-skew][z-skew] +
                                                               u[t0][x-skew + 1][y-skew][z-skew])
                                                    + 4.11e-3F*u[t0][x-skew][y-skew][z-skew]
                                                    - 3.94e-1F*u[t0][x-skew][y-skew][z-skew])/(3.04e-1F);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
void forward_tiled_mod(float u[TIME][DIM][DIM][DIM], const int time_size, const int x_size, const int y_size, const int z_size){
    for (int tt=0;tt<=floord(time_size-2,4);tt++) {
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for(int xx=0;xx<=11;xx++) {
            for(int yy=0;yy<=11;yy++) {
                for(int zz=0;zz<=11;zz++) {
                    for (int time=max(1,4*tt);time<=min(time_size-2,4*tt+3);time++) {
                        int skew = 4*time;
                        int t0 = (time) % 8;
                        int t1 = (time + 1) % 8;
                        int t2 = (time - 1) % 8;
                        for(int x=max(4*time+4,24*xx+4*time);x<=min(4*time+271,24*xx+4*time+23);x++) {
                            for(int y=max(4*time+4,24*yy+4*time);y<=min(4*time+271,24*yy+4*time+23);y++) {
                                #pragma ivdep
                                #pragma omp simd
                                for(int z=max(4*time+4,24*zz+4*time);z<=min(4*time+271,24*zz+4*time+23);z++) {
                                    u[t1][x-skew][y-skew][z-skew] = 1e-1F +
                                                     (3.04e-1F*u[t2][x-skew][y-skew][z-skew]
                                                   - 8.25e-5F*(u[t0][x-skew][y-skew][z-skew - 4] +
                                                               u[t0][x-skew][y-skew][z-skew + 4] +
                                                               u[t0][x-skew][y-skew - 4][z-skew] +
                                                               u[t0][x-skew][y-skew + 4][z-skew] +
                                                               u[t0][x-skew - 4][y-skew][z-skew] +
                                                               u[t0][x-skew + 4][y-skew][z-skew])
                                                   + 1.17e-3F*(u[t0][x-skew][y-skew][z-skew - 3] +
                                                               u[t0][x-skew][y-skew][z-skew + 3] +
                                                               u[t0][x-skew][y-skew - 3][z-skew] +
                                                               u[t0][x-skew][y-skew + 3][z-skew] +
                                                               u[t0][x-skew - 3][y-skew][z-skew] +
                                                               u[t0][x-skew + 3][y-skew][z-skew])
                                                   - 9.24e-3F*(u[t0][x-skew][y-skew][z-skew - 2] +
                                                               u[t0][x-skew][y-skew][z-skew + 2] +
                                                               u[t0][x-skew][y-skew - 2][z-skew] +
                                                               u[t0][x-skew][y-skew + 2][z-skew] +
                                                               u[t0][x-skew - 2][y-skew][z-skew] +
                                                               u[t0][x-skew + 2][y-skew][z-skew])
                                                + 7.39e-2F*(u[t0][x-skew][y-skew][z-skew - 1] +
                                                               u[t0][x-skew][y-skew][z-skew + 1] +
                                                               u[t0][x-skew][y-skew - 1][z-skew] +
                                                               u[t0][x-skew][y-skew + 1][z-skew] +
                                                               u[t0][x-skew - 1][y-skew][z-skew] +
                                                               u[t0][x-skew + 1][y-skew][z-skew])
                                                    + 4.11e-3F*u[t0][x-skew][y-skew][z-skew]
                                                    - 3.94e-1F*u[t0][x-skew][y-skew][z-skew])/(3.04e-1F);
                                }
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
    int t0 = (time) % 8;
    int t1 = (time + 1) % 8;
    int t2 = (time - 1) % 8;
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
        }
      }
    }
    }
}

int compare_output(const int time_size, const int x_size, const int y_size, const int z_size){
    int equal = 0;
    for (int time = 1; time < TIME -1; time += 1)
    {
        {
            #pragma omp for schedule(static)
            for (int x = 4; x < x_size - 4; x += 1)
            {
                for (int y = 4; y < y_size - 4; y += 1)
                {
                    for (int z = 4; z < z_size - 4; z += 1)
                    {
                        equal += fabs(u_grid[time][x][y][z] - v_grid[time][x][y][z]) >= 10e-5;
                    }
                }
            }
        }
    }
    return equal;
}

int main(){
    forward_control(u_grid, TIME, DIM, DIM, DIM);
    forward_skewed_tiled(v_grid, TIME, DIM, DIM, DIM);
    printf("%d", compare_output(TIME, DIM, DIM, DIM));
    /*forward_tiled(v, TIME, DIM, DIM, DIM);*/
}


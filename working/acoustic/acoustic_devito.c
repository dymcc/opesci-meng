#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include "xmmintrin.h"
#include "pmmintrin.h"

#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define X_TILE_SIZE 16
#define Y_TILE_SIZE 16
#define Z_TILE_SIZE 16
#define SKEW 4*time

struct profiler
{
  double loop_x_0;
  double loop_p_src_1;
  double loop_p_rec_2;
} ;

void f_1_0(float *restrict damp_vec, const int x_size, const int x, const int y_size, const int y, const int z_size, const int zz, float *restrict m_vec, float *restrict u_vec, const int time_size, const int time)
{
  float (*restrict damp)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) damp_vec;
  float (*restrict m)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) m_vec;
  float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
  #pragma omp simd
  for (int z=max(4*time+4,4*time+27*zz);z<=min(4*time+211,4*time+27*zz+26);z++) {
    u[time + 1][x][y][z-4*time] = ((3.04F*damp[x][y][z-4*time] - 2*m[x][y][z-4*time])*u[time - 1][x][y][z-4*time] - 8.25142857142857e-5F*(u[time][x][y][z-4*time - 4] + u[time][x][y][z-4*time + 4] + u[time][x][y - 4][z-4*time] + u[time][x][y + 4][z-4*time] + u[time][x - 4][y][z-4*time] + u[time][x + 4][y][z-4*time]) + 1.17353650793651e-3F*(u[time][x][y][z-4*time - 3] + u[time][x][y][z-4*time + 3] + u[time][x][y - 3][z-4*time] + u[time][x][y + 3][z-4*time] + u[time][x - 3][y][z-4*time] + u[time][x + 3][y][z-4*time]) - 9.2416e-3F*(u[time][x][y][z-4*time - 2] + u[time][x][y][z-4*time + 2] + u[time][x][y - 2][z-4*time] + u[time][x][y + 2][z-4*time] + u[time][x - 2][y][z-4*time] + u[time][x + 2][y][z-4*time]) + 7.39328e-2F*(u[time][x][y][z-4*time - 1] + u[time][x][y][z-4*time + 1] + u[time][x][y - 1][z-4*time] + u[time][x][y + 1][z-4*time] + u[time][x - 1][y][z-4*time] + u[time][x + 1][y][z-4*time]) + 4*m[x][y][z-4*time]*u[time][x][y][z-4*time] - 3.94693333333333e-1F*u[time][x][y][z-4*time])/(3.04F*damp[x][y][z-4*time] + 2*m[x][y][z-4*time]);
  }
}
void f_1_1(float *restrict m_vec, const int x_size, const int y_size, const int z_size, float *restrict src_vec, const int time_size, const int time, float *restrict src_coords_vec, const int d_size, float *restrict u_vec)
{
  float (*restrict m)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) m_vec;
  float (*restrict src)[1] __attribute__((aligned(64))) = (float (*)[1]) src_vec;
  float (*restrict src_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) src_coords_vec;
  float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
  for (int p_src = 0; p_src < 1; p_src += 1)
  {
    int temp11 = (int)(floor(5.0e-2F*src_coords[p_src][1]));
    int temp14 = (int)(floor(5.0e-2F*src_coords[p_src][2]));
    int temp8 = (int)(floor(5.0e-2F*src_coords[p_src][0]));
    float temp18 = (float)(-2.0e+1F*temp8 + src_coords[p_src][0]);
    float temp19 = (float)(-2.0e+1F*temp11 + src_coords[p_src][1]);
    float temp20 = (float)(-2.0e+1F*temp14 + src_coords[p_src][2]);
    u[time + 1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(-1.25e-4F*temp18*temp19*temp20 - 5.0e-2F*(temp18 + temp19 + temp20) + 2.5e-3F*(temp18*temp19 + temp18*temp20 + temp19*temp20) + 1)*src[time][p_src]/m[temp8 + 10][temp11 + 10][temp14 + 10] + u[time + 1][temp8 + 10][temp11 + 10][temp14 + 10];
    u[time + 1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(1.25e-4F*temp18*temp19*temp20 + 5.0e-2F*temp19 - 2.5e-3F*(temp18*temp19 + temp19*temp20))*src[time][p_src]/m[temp8 + 10][temp11 + 11][temp14 + 10] + u[time + 1][temp8 + 10][temp11 + 11][temp14 + 10];
    u[time + 1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(1.25e-4F*temp18*temp19*temp20 + 5.0e-2F*temp18 - 2.5e-3F*(temp18*temp19 + temp18*temp20))*src[time][p_src]/m[temp8 + 11][temp11 + 10][temp14 + 10] + u[time + 1][temp8 + 11][temp11 + 10][temp14 + 10];
    u[time + 1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 9.2416F*(1.25e-4F*temp18*temp19*temp20 + 5.0e-2F*temp20 - 2.5e-3F*(temp18*temp20 + temp19*temp20))*src[time][p_src]/m[temp8 + 10][temp11 + 10][temp14 + 11] + u[time + 1][temp8 + 10][temp11 + 10][temp14 + 11];
    u[time + 1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(-1.25e-4F*temp18*temp19*temp20 + 2.5e-3F*temp18*temp19)*src[time][p_src]/m[temp8 + 11][temp11 + 11][temp14 + 10] + u[time + 1][temp8 + 11][temp11 + 11][temp14 + 10];
    u[time + 1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 9.2416F*(-1.25e-4F*temp18*temp19*temp20 + 2.5e-3F*temp19*temp20)*src[time][p_src]/m[temp8 + 10][temp11 + 11][temp14 + 11] + u[time + 1][temp8 + 10][temp11 + 11][temp14 + 11];
    u[time + 1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 9.2416F*(-1.25e-4F*temp18*temp19*temp20 + 2.5e-3F*temp18*temp20)*src[time][p_src]/m[temp8 + 11][temp11 + 10][temp14 + 11] + u[time + 1][temp8 + 11][temp11 + 10][temp14 + 11];
    u[time + 1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 1.1552e-3F*temp18*temp19*temp20*src[time][p_src]/m[temp8 + 11][temp11 + 11][temp14 + 11] + u[time + 1][temp8 + 11][temp11 + 11][temp14 + 11];
  }
}
void f_1_2(float *restrict rec_vec, const int time_size, const int time, float *restrict rec_coords_vec, const int d_size, float *restrict u_vec, const int x_size, const int y_size, const int z_size)
{
  float (*restrict rec)[101] __attribute__((aligned(64))) = (float (*)[101]) rec_vec;
  float (*restrict rec_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) rec_coords_vec;
  float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
  for (int p_rec = 0; p_rec < 101; p_rec += 1)
  {
    int temp41 = (int)(floor(5.0e-2F*rec_coords[p_rec][0]));
    int temp44 = (int)(floor(5.0e-2F*rec_coords[p_rec][1]));
    int temp47 = (int)(floor(5.0e-2F*rec_coords[p_rec][2]));
    float temp42 = (float)(-2.0e+1F*temp41 + rec_coords[p_rec][0]);
    float temp45 = (float)(-2.0e+1F*temp44 + rec_coords[p_rec][1]);
    float temp48 = (float)(-2.0e+1F*temp47 + rec_coords[p_rec][2]);
    rec[time][p_rec] = 1.25e-4F*temp42*temp45*temp48*u[time][temp41 + 11][temp44 + 11][temp47 + 11] + (-1.25e-4F*temp42*temp45*temp48 + 2.5e-3F*temp42*temp45)*u[time][temp41 + 11][temp44 + 11][temp47 + 10] + (-1.25e-4F*temp42*temp45*temp48 + 2.5e-3F*temp42*temp48)*u[time][temp41 + 11][temp44 + 10][temp47 + 11] + (-1.25e-4F*temp42*temp45*temp48 + 2.5e-3F*temp45*temp48)*u[time][temp41 + 10][temp44 + 11][temp47 + 11] + (1.25e-4F*temp42*temp45*temp48 + 5.0e-2F*temp42 - 2.5e-3F*(temp42*temp45 + temp42*temp48))*u[time][temp41 + 11][temp44 + 10][temp47 + 10] + (1.25e-4F*temp42*temp45*temp48 + 5.0e-2F*temp45 - 2.5e-3F*(temp42*temp45 + temp45*temp48))*u[time][temp41 + 10][temp44 + 11][temp47 + 10] + (1.25e-4F*temp42*temp45*temp48 + 5.0e-2F*temp48 - 2.5e-3F*(temp42*temp48 + temp45*temp48))*u[time][temp41 + 10][temp44 + 10][temp47 + 11] + (-1.25e-4F*temp42*temp45*temp48 - 5.0e-2F*(temp42 + temp45 + temp48) + 2.5e-3F*(temp42*temp45 + temp42*temp48 + temp45*temp48) + 1)*u[time][temp41 + 10][temp44 + 10][temp47 + 10];
  }
}

int Forward(float *restrict damp_vec, float *restrict m_vec, float *restrict rec_vec, float *restrict rec_coords_vec, float *restrict src_vec, float *restrict src_coords_vec, float *restrict u_vec, const int d_size, const int time_size, const int x_size, const int y_size, const int z_size, struct profiler *timings)
{
  float (*restrict damp)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) damp_vec;
  float (*restrict m)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) m_vec;
  float (*restrict rec)[101] __attribute__((aligned(64))) = (float (*)[101]) rec_vec;
  float (*restrict rec_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) rec_coords_vec;
  float (*restrict src)[1] __attribute__((aligned(64))) = (float (*)[1]) src_vec;
  float (*restrict src_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) src_coords_vec;
  float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
  /* Flush denormal numbers to zero in hardware */
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  int time, xx, yy, zz, x, y, z;
  for (time=1;time<=81;time++) {
    struct timeval start_loop_x_0, end_loop_x_0;
    gettimeofday(&start_loop_x_0, NULL);
    for (xx=0;xx<=7;xx++) {
      for (yy=0;yy<=7;yy++) {
        for (zz=0;zz<=7;zz++) {
          for (x=max(4*time+4,4*time+27*xx);x<=min(4*time+211,4*time+27*xx+26);x++) {
            for (y=max(4*time+4,4*time+27*yy);y<=min(4*time+211,4*time+27*yy+26);y++) {
              f_1_0(damp_vec,x_size,x-4*time,y_size,y-4*time,z_size,zz,m_vec,u_vec,time_size,time);
            }
          }
        }
      }
    }
      gettimeofday(&end_loop_x_0, NULL);
      timings->loop_x_0 += (double)(end_loop_x_0.tv_sec-start_loop_x_0.tv_sec)+(double)(end_loop_x_0.tv_usec-start_loop_x_0.tv_usec)/1000000;
      struct timeval start_loop_p_src_1, end_loop_p_src_1;
      gettimeofday(&start_loop_p_src_1, NULL);
      /* noinline? */
      f_1_1(m_vec,x_size,y_size,z_size,src_vec,time_size,time,src_coords_vec,d_size,u_vec);
      gettimeofday(&end_loop_p_src_1, NULL);
      timings->loop_p_src_1 += (double)(end_loop_p_src_1.tv_sec-start_loop_p_src_1.tv_sec)+(double)(end_loop_p_src_1.tv_usec-start_loop_p_src_1.tv_usec)/1000000;
      struct timeval start_loop_p_rec_2, end_loop_p_rec_2;
      gettimeofday(&start_loop_p_rec_2, NULL);
      /* noinline? */
      f_1_2(rec_vec,time_size,time,rec_coords_vec,d_size,u_vec,x_size,y_size,z_size);
      gettimeofday(&end_loop_p_rec_2, NULL);
      timings->loop_p_rec_2 += (double)(end_loop_p_rec_2.tv_sec-start_loop_p_rec_2.tv_sec)+(double)(end_loop_p_rec_2.tv_usec-start_loop_p_rec_2.tv_usec)/1000000;
    }
  return 0;
}

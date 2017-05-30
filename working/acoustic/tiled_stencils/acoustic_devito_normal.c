#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "sys/time.h"
#include "xmmintrin.h"
#include "pmmintrin.h"

struct profiler
{
  double loop_x_0;
  double loop_p_src_1;
  double loop_p_rec_2;
} ;


int Forward(float *restrict damp_vec, float *restrict m_vec, float *restrict rec_vec, float *restrict rec_coords_vec, float *restrict src_vec, float *restrict src_coords_vec, float *restrict u_vec, const int d_size, const int p_rec_size, const int p_src_size, const int t_size, const int time_size, const int x_size, const int y_size, const int z_size, struct profiler *timings)
{
  float (*restrict damp)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) damp_vec;
  float (*restrict m)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) m_vec;
  float (*restrict rec)[p_rec_size] __attribute__((aligned(64))) = (float (*)[p_rec_size]) rec_vec;
  float (*restrict rec_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) rec_coords_vec;
  float (*restrict src)[p_src_size] __attribute__((aligned(64))) = (float (*)[p_src_size]) src_vec;
  float (*restrict src_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) src_coords_vec;
  float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
  /* Flush denormal numbers to zero in hardware */
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  printf("%d %d %d %d", time_size, x_size, y_size, z_size);
  for (int time = 1; time < time_size - 1; time += 1)
  {
    int t0 = (time) % 3;
    int t1 = (time + 1) % 3;
    int t2 = (time - 1) % 3;
    struct timeval start_loop_x_0, end_loop_x_0;
    gettimeofday(&start_loop_x_0, NULL);
    for (int x = 4; x < x_size - 4; x += 1)
    {
      for (int y = 4; y < y_size - 4; y += 1)
      {
        for (int z = 4; z < z_size - 4; z += 1)
        {
          float tcse0 = 3.04F*damp[x][y][z];
          u[t1][x][y][z] = ((tcse0 - 2*m[x][y][z])*u[t2][x][y][z] - 8.25142857142857e-5F*(u[t0][x][y][z - 4] + u[t0][x][y][z + 4] + u[t0][x][y - 4][z] + u[t0][x][y + 4][z] + u[t0][x - 4][y][z] + u[t0][x + 4][y][z]) + 1.17353650793651e-3F*(u[t0][x][y][z - 3] + u[t0][x][y][z + 3] + u[t0][x][y - 3][z] + u[t0][x][y + 3][z] + u[t0][x - 3][y][z] + u[t0][x + 3][y][z]) - 9.2416e-3F*(u[t0][x][y][z - 2] + u[t0][x][y][z + 2] + u[t0][x][y - 2][z] + u[t0][x][y + 2][z] + u[t0][x - 2][y][z] + u[t0][x + 2][y][z]) + 7.39328e-2F*(u[t0][x][y][z - 1] + u[t0][x][y][z + 1] + u[t0][x][y - 1][z] + u[t0][x][y + 1][z] + u[t0][x - 1][y][z] + u[t0][x + 1][y][z]) + 4*m[x][y][z]*u[t0][x][y][z] - 3.94693333333333e-1F*u[t0][x][y][z])/(tcse0 + 2*m[x][y][z]);
        }
      }
    }
    gettimeofday(&end_loop_x_0, NULL);
    timings->loop_x_0 += (double)(end_loop_x_0.tv_sec-start_loop_x_0.tv_sec)+(double)(end_loop_x_0.tv_usec-start_loop_x_0.tv_usec)/1000000;
    struct timeval start_loop_p_src_1, end_loop_p_src_1;
    gettimeofday(&start_loop_p_src_1, NULL);
    for (int p_src = 0; p_src < p_src_size; p_src += 1)
    {
      float tcse0 = (float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p_src][1])) + src_coords[p_src][1]);
      float tcse1 = (float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p_src][2])) + src_coords[p_src][2]);
      float tcse2 = (float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p_src][0])) + src_coords[p_src][0]);
      float tcse3 = -2.5e-3F*tcse1*tcse2;
      float tcse4 = -2.5e-3F*tcse0*tcse1;
      float tcse5 = -2.5e-3F*tcse0*tcse2;
      float tcse6 = 2.5e-3F*tcse0*tcse1;
      float tcse7 = 2.5e-3F*tcse1*tcse2;
      float tcse8 = 2.5e-3F*tcse0*tcse2;
      float tcse9 = 1.25e-4F*tcse0*tcse1*tcse2;
      float tcse10 = -1.25e-4F*tcse0*tcse1*tcse2;
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(-5.0e-2F*tcse0 - 5.0e-2F*tcse1 + tcse10 - 5.0e-2F*tcse2 + tcse6 + tcse7 + tcse8 + 1)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(5.0e-2F*tcse0 + tcse4 + tcse5 + tcse9)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(5.0e-2F*tcse2 + tcse3 + tcse5 + tcse9)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 9.2416F*(5.0e-2F*tcse1 + tcse3 + tcse4 + tcse9)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(tcse10 + tcse8)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 9.2416F*(tcse10 + tcse6)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 9.2416F*(tcse10 + tcse7)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 1.1552e-3F*tcse0*tcse1*tcse2*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
    }
    gettimeofday(&end_loop_p_src_1, NULL);
    timings->loop_p_src_1 += (double)(end_loop_p_src_1.tv_sec-start_loop_p_src_1.tv_sec)+(double)(end_loop_p_src_1.tv_usec-start_loop_p_src_1.tv_usec)/1000000;
    struct timeval start_loop_p_rec_2, end_loop_p_rec_2;
    gettimeofday(&start_loop_p_rec_2, NULL);
    for (int p_rec = 0; p_rec < p_rec_size; p_rec += 1)
    {
      float tcse0 = (float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + rec_coords[p_rec][0]);
      float tcse1 = (float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + rec_coords[p_rec][2]);
      float tcse2 = (float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + rec_coords[p_rec][1]);
      float tcse3 = -2.5e-3F*tcse1*tcse2;
      float tcse4 = -2.5e-3F*tcse0*tcse1;
      float tcse5 = -2.5e-3F*tcse0*tcse2;
      float tcse6 = 2.5e-3F*tcse1*tcse2;
      float tcse7 = 2.5e-3F*tcse0*tcse1;
      float tcse8 = 2.5e-3F*tcse0*tcse2;
      float tcse9 = 1.25e-4F*tcse0*tcse1*tcse2;
      float tcse10 = -1.25e-4F*tcse0*tcse1*tcse2;
      rec[time][p_rec] = 1.25e-4F*tcse0*tcse1*tcse2*u[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 11] + (tcse10 + tcse6)*u[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 11] + (tcse10 + tcse7)*u[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 11] + (tcse10 + tcse8)*u[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 10] + (5.0e-2F*tcse0 + tcse4 + tcse5 + tcse9)*u[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 10] + (5.0e-2F*tcse1 + tcse3 + tcse4 + tcse9)*u[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 11] + (5.0e-2F*tcse2 + tcse3 + tcse5 + tcse9)*u[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 10] + (-5.0e-2F*tcse0 - 5.0e-2F*tcse1 + tcse10 - 5.0e-2F*tcse2 + tcse6 + tcse7 + tcse8 + 1)*u[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 10];
    }
    gettimeofday(&end_loop_p_rec_2, NULL);
    timings->loop_p_rec_2 += (double)(end_loop_p_rec_2.tv_sec-start_loop_p_rec_2.tv_sec)+(double)(end_loop_p_rec_2.tv_usec-start_loop_p_rec_2.tv_usec)/1000000;
  }
  return 0;
}

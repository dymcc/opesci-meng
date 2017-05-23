#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
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

void f_1_0(float *restrict damp_vec, const int x_size, const int x, const int y_size, const int y, const int z_size, float *restrict m_vec, float *restrict u_vec, const int t_size, const int t0, const int t1, const int t2)
{
  float (*restrict damp)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) damp_vec;
  float (*restrict m)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) m_vec;
  float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
  #pragma ivdep
  #pragma omp simd
  for (int z = 4; z < z_size - 4; z += 1)
  {
    u[t1][x][y][z] = ((3.04F*damp[x][y][z] - 2*m[x][y][z])*u[t2][x][y][z] - 8.25142857142857e-5F*(u[t0][x][y][z - 4] + u[t0][x][y][z + 4] + u[t0][x][y - 4][z] + u[t0][x][y + 4][z] + u[t0][x - 4][y][z] + u[t0][x + 4][y][z]) + 1.17353650793651e-3F*(u[t0][x][y][z - 3] + u[t0][x][y][z + 3] + u[t0][x][y - 3][z] + u[t0][x][y + 3][z] + u[t0][x - 3][y][z] + u[t0][x + 3][y][z]) - 9.2416e-3F*(u[t0][x][y][z - 2] + u[t0][x][y][z + 2] + u[t0][x][y - 2][z] + u[t0][x][y + 2][z] + u[t0][x - 2][y][z] + u[t0][x + 2][y][z]) + 7.39328e-2F*(u[t0][x][y][z - 1] + u[t0][x][y][z + 1] + u[t0][x][y - 1][z] + u[t0][x][y + 1][z] + u[t0][x - 1][y][z] + u[t0][x + 1][y][z]) + 4*m[x][y][z]*u[t0][x][y][z] - 3.94693333333333e-1F*u[t0][x][y][z])/(3.04F*damp[x][y][z] + 2*m[x][y][z]);
  }
}
void f_1_1(float *restrict m_vec, const int x_size, const int y_size, const int z_size, float *restrict src_vec, const int time_size, float *restrict src_coords_vec, const int d_size, float *restrict u_vec, const int t_size, const int t1, const int time)
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
    u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(-1.25e-4F*temp18*temp19*temp20 - 5.0e-2F*(temp18 + temp19 + temp20) + 2.5e-3F*(temp18*temp19 + temp18*temp20 + temp19*temp20) + 1)*src[time][p_src]/m[temp8 + 10][temp11 + 10][temp14 + 10] + u[t1][temp8 + 10][temp11 + 10][temp14 + 10];
    u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(1.25e-4F*temp18*temp19*temp20 + 5.0e-2F*temp19 - 2.5e-3F*(temp18*temp19 + temp19*temp20))*src[time][p_src]/m[temp8 + 10][temp11 + 11][temp14 + 10] + u[t1][temp8 + 10][temp11 + 11][temp14 + 10];
    u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(1.25e-4F*temp18*temp19*temp20 + 5.0e-2F*temp18 - 2.5e-3F*(temp18*temp19 + temp18*temp20))*src[time][p_src]/m[temp8 + 11][temp11 + 10][temp14 + 10] + u[t1][temp8 + 11][temp11 + 10][temp14 + 10];
    u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 9.2416F*(1.25e-4F*temp18*temp19*temp20 + 5.0e-2F*temp20 - 2.5e-3F*(temp18*temp20 + temp19*temp20))*src[time][p_src]/m[temp8 + 10][temp11 + 10][temp14 + 11] + u[t1][temp8 + 10][temp11 + 10][temp14 + 11];
    u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = 9.2416F*(-1.25e-4F*temp18*temp19*temp20 + 2.5e-3F*temp18*temp19)*src[time][p_src]/m[temp8 + 11][temp11 + 11][temp14 + 10] + u[t1][temp8 + 11][temp11 + 11][temp14 + 10];
    u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 9.2416F*(-1.25e-4F*temp18*temp19*temp20 + 2.5e-3F*temp19*temp20)*src[time][p_src]/m[temp8 + 10][temp11 + 11][temp14 + 11] + u[t1][temp8 + 10][temp11 + 11][temp14 + 11];
    u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 9.2416F*(-1.25e-4F*temp18*temp19*temp20 + 2.5e-3F*temp18*temp20)*src[time][p_src]/m[temp8 + 11][temp11 + 10][temp14 + 11] + u[t1][temp8 + 11][temp11 + 10][temp14 + 11];
    u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = 1.1552e-3F*temp18*temp19*temp20*src[time][p_src]/m[temp8 + 11][temp11 + 11][temp14 + 11] + u[t1][temp8 + 11][temp11 + 11][temp14 + 11];
  }
}
void f_1_2(float *restrict rec_vec, const int time_size, float *restrict rec_coords_vec, const int d_size, float *restrict u_vec, const int t_size, const int x_size, const int y_size, const int z_size, const int t0, const int time)
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
    rec[time][p_rec] = 1.25e-4F*temp42*temp45*temp48*u[t0][temp41 + 11][temp44 + 11][temp47 + 11] + (-1.25e-4F*temp42*temp45*temp48 + 2.5e-3F*temp42*temp45)*u[t0][temp41 + 11][temp44 + 11][temp47 + 10] + (-1.25e-4F*temp42*temp45*temp48 + 2.5e-3F*temp42*temp48)*u[t0][temp41 + 11][temp44 + 10][temp47 + 11] + (-1.25e-4F*temp42*temp45*temp48 + 2.5e-3F*temp45*temp48)*u[t0][temp41 + 10][temp44 + 11][temp47 + 11] + (1.25e-4F*temp42*temp45*temp48 + 5.0e-2F*temp42 - 2.5e-3F*(temp42*temp45 + temp42*temp48))*u[t0][temp41 + 11][temp44 + 10][temp47 + 10] + (1.25e-4F*temp42*temp45*temp48 + 5.0e-2F*temp45 - 2.5e-3F*(temp42*temp45 + temp45*temp48))*u[t0][temp41 + 10][temp44 + 11][temp47 + 10] + (1.25e-4F*temp42*temp45*temp48 + 5.0e-2F*temp48 - 2.5e-3F*(temp42*temp48 + temp45*temp48))*u[t0][temp41 + 10][temp44 + 10][temp47 + 11] + (-1.25e-4F*temp42*temp45*temp48 - 5.0e-2F*(temp42 + temp45 + temp48) + 2.5e-3F*(temp42*temp45 + temp42*temp48 + temp45*temp48) + 1)*u[t0][temp41 + 10][temp44 + 10][temp47 + 10];
  }
}

int Forward(float *restrict damp_vec, float *restrict m_vec, float *restrict rec_vec, float *restrict rec_coords_vec, float *restrict src_vec, float *restrict src_coords_vec, float *restrict u_vec, const int d_size, const int t_size, const int time_size, const int x_size, const int y_size, const int z_size, struct profiler *timings)
{
  float (*restrict damp)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) damp_vec;
  float (*restrict m)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) m_vec;
  float (*restrict rec)[101] __attribute__((aligned(64))) = (float (*)[101]) rec_vec;
  float (*restrict rec_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) rec_coords_vec;
  float (*restrict src)[1] __attribute__((aligned(64))) = (float (*)[1]) src_vec;
  float (*restrict src_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) src_coords_vec;
  float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
  /* DLE: moved denormals flag */
  for (int time = 1; time < time_size - 1; time += 1)
  {
    int t0 = (time) % 3;
    int t1 = (time + 1) % 3;
    int t2 = (time - 1) % 3;
    struct timeval start_loop_x_0, end_loop_x_0;
    gettimeofday(&start_loop_x_0, NULL);
    #pragma omp parallel
    {
      /* Flush denormal numbers to zero in hardware */
      _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
      _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
      #pragma omp for schedule(static)
      for (int x = 4; x < x_size - 4; x += 1)
      {
        for (int y = 4; y < y_size - 4; y += 1)
        {
          #pragma noinline
          f_1_0(damp_vec,x_size,x,y_size,y,z_size,m_vec,u_vec,t_size,t0,t1,t2);
        }
      }
    }
    gettimeofday(&end_loop_x_0, NULL);
    timings->loop_x_0 += (double)(end_loop_x_0.tv_sec-start_loop_x_0.tv_sec)+(double)(end_loop_x_0.tv_usec-start_loop_x_0.tv_usec)/1000000;
    struct timeval start_loop_p_src_1, end_loop_p_src_1;
    gettimeofday(&start_loop_p_src_1, NULL);
    #pragma noinline
    f_1_1(m_vec,x_size,y_size,z_size,src_vec,time_size,src_coords_vec,d_size,u_vec,t_size,t1,time);
    gettimeofday(&end_loop_p_src_1, NULL);
    timings->loop_p_src_1 += (double)(end_loop_p_src_1.tv_sec-start_loop_p_src_1.tv_sec)+(double)(end_loop_p_src_1.tv_usec-start_loop_p_src_1.tv_usec)/1000000;
    struct timeval start_loop_p_rec_2, end_loop_p_rec_2;
    gettimeofday(&start_loop_p_rec_2, NULL);
    #pragma noinline
    f_1_2(rec_vec,time_size,rec_coords_vec,d_size,u_vec,t_size,x_size,y_size,z_size,t0,time);
    gettimeofday(&end_loop_p_rec_2, NULL);
    timings->loop_p_rec_2 += (double)(end_loop_p_rec_2.tv_sec-start_loop_p_rec_2.tv_sec)+(double)(end_loop_p_rec_2.tv_usec-start_loop_p_rec_2.tv_usec)/1000000;
  }
  return 0;
}

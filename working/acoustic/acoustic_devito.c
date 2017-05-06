#include "assert.h"
#include "stdlib.h"
#include "math.h"
#include "stdio.h"
#include "string.h"
#include "sys/time.h"
#include "sys/param.h"
struct profiler
{
  double post_stencils0;
  double post_stencils1;
  double loop_body;
} ;
int ForwardOperator(float *u_vec, float *damp_vec, float *m_vec, float *src_vec, float *src_coords_vec, float *rec_vec, float *rec_coords_vec, struct profiler *timings)
{
  float (*u)[320][320][320] = (float (*)[320][320][320]) u_vec;
  float (*damp)[320][320] = (float (*)[320][320]) damp_vec;
  float (*m)[320][320] = (float (*)[320][320]) m_vec;
  float (*src)[1] = (float (*)[1]) src_vec;
  float (*src_coords)[3] = (float (*)[3]) src_coords_vec;
  float (*rec)[101] = (float (*)[101]) rec_vec;
  float (*rec_coords)[3] = (float (*)[3]) rec_coords_vec;
  {
    int t0;
    int t1;
    int t2;
      for (int i4 = 0; i4<83; i4++)
      {
        {
          t0 = i4;
          t1 = t0 + 1;
          t2 = t1 + 1;
        }
        struct timeval start_loop_body, end_loop_body;
        gettimeofday(&start_loop_body, NULL);
        {
          for (int i1 = 3; i1<317; i1++)
            for (int i2 = 3; i2<317; i2++)
            {
              #pragma GCC ivdep
              for (int i3 = 3; i3<317; i3++)
              {
                float temp1 = 3.04F*damp[i1][i2][i3];
                float temp3 = 2*m[i1][i2][i3];
                float temp6 = 2.50000000000000e-3F;
                u[t2][i1][i2][i3] = (1.84832e+1F*temp6*(u[t1][i1][i2][i3 - 1] + u[t1][i1][i2][i3 + 1] + u[t1][i1][i2 - 1][i3] + u[t1][i1][i2 + 1][i3] + u[t1][i1 - 1][i2][i3] + u[t1][i1 + 1][i2][i3]) - 1.108992e+2F*temp6*u[t1][i1][i2][i3] + (temp1 - temp3)*u[t0][i1][i2][i3] + 4*m[i1][i2][i3]*u[t1][i1][i2][i3])/(temp1 + temp3);
              }
            }
        }
        {
          gettimeofday(&end_loop_body, NULL);
          timings->loop_body += (double)(end_loop_body.tv_sec-start_loop_body.tv_sec)+(double)(end_loop_body.tv_usec-start_loop_body.tv_usec)/1000000;
        }
        {
          struct timeval start_post_stencils0, end_post_stencils0;
          gettimeofday(&start_post_stencils0, NULL);
          for (int p1 = 0; p1 < 1; p1++)
          {
            u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10] = 9.2416F*(-1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) - 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) - 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1]) - 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) + 1)*src[i4][p1]/m[(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10] + u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10];
            u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10] = 9.2416F*(1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) + 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1]))*src[i4][p1]/m[(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10] + u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10];
            u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10] = 9.2416F*(1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) + 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0]))*src[i4][p1]/m[(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10] + u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10];
            u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11] = 9.2416F*(1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) + 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]))*src[i4][p1]/m[(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11] + u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11];
            u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10] = 9.2416F*(-1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1]))*src[i4][p1]/m[(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10] + u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 10];
            u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11] = 9.2416F*(-1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]))*src[i4][p1]/m[(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11] + u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 10][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11];
            u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11] = 9.2416F*(-1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2]))*src[i4][p1]/m[(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11] + u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 10][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11];
            u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11] = 1.1552e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][0])) + src_coords[p1][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][1])) + src_coords[p1][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p1][2])) + src_coords[p1][2])*src[i4][p1]/m[(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11] + u[t2][(int)(floor(5.0e-2F*src_coords[p1][0])) + 11][(int)(floor(5.0e-2F*src_coords[p1][1])) + 11][(int)(floor(5.0e-2F*src_coords[p1][2])) + 11];
          }
          {
            gettimeofday(&end_post_stencils0, NULL);
            timings->post_stencils0 += (double)(end_post_stencils0.tv_sec-start_post_stencils0.tv_sec)+(double)(end_post_stencils0.tv_usec-start_post_stencils0.tv_usec)/1000000;
          }
          struct timeval start_post_stencils1, end_post_stencils1;
          gettimeofday(&start_post_stencils1, NULL);
          for (int p2 = 0; p2 < 101; p2++)
          {
            rec[i4][p2] = (-1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1]))*u[t2][(int)(floor(5.0e-2F*rec_coords[p2][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p2][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p2][2])) + 10] + (-1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]))*u[t2][(int)(floor(5.0e-2F*rec_coords[p2][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p2][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p2][2])) + 11] + (-1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]))*u[t2][(int)(floor(5.0e-2F*rec_coords[p2][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p2][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p2][2])) + 11] + (1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) + 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0]))*u[t2][(int)(floor(5.0e-2F*rec_coords[p2][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p2][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p2][2])) + 10] + (1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) + 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1]))*u[t2][(int)(floor(5.0e-2F*rec_coords[p2][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p2][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p2][2])) + 10] + (1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) - 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) + 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]))*u[t2][(int)(floor(5.0e-2F*rec_coords[p2][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p2][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p2][2])) + 11] + (-1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) - 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0]) + 2.5e-3F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) - 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1]) - 5.0e-2F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2]) + 1)*u[t2][(int)(floor(5.0e-2F*rec_coords[p2][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p2][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p2][2])) + 10] + 1.25e-4F*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][0])) + rec_coords[p2][0])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][1])) + rec_coords[p2][1])*(float)(-2.0e+1F*(int)(floor(5.0e-2F*rec_coords[p2][2])) + rec_coords[p2][2])*u[t2][(int)(floor(5.0e-2F*rec_coords[p2][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p2][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p2][2])) + 11];
          }
          {
            gettimeofday(&end_post_stencils1, NULL);
            timings->post_stencils1 += (double)(end_post_stencils1.tv_sec-start_post_stencils1.tv_sec)+(double)(end_post_stencils1.tv_usec-start_post_stencils1.tv_usec)/1000000;
          }
        }
      }
  }
  return 0;
}

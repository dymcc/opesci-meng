#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "stdio.h"
#include "signal.h"
#include "math.h"
#include "sys/time.h"
#include "xmmintrin.h"
#include "pmmintrin.h"

#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

struct profiler
{
  double loop_x_0;
  double loop_p_rec_1;
} ;


int Forward(float *restrict damp_vec, float *restrict m_vec, float *restrict rec_vec, float *restrict rec_coords_vec, float *restrict u_vec, const int d_size, const int p_rec_size, const int t_size, const int time_size, const int x_size, const int y_size, const int z_size, struct profiler *timings)
{
    float (*restrict damp)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) damp_vec;
    float (*restrict m)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) m_vec;
    float (*restrict rec)[p_rec_size] __attribute__((aligned(64))) = (float (*)[p_rec_size]) rec_vec;
    float (*restrict rec_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) rec_coords_vec;
    float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
    /* DLE: moved denormals flag */
    struct timeval start_loop_x_0, end_loop_x_0;
    gettimeofday(&start_loop_x_0, NULL);
    for (int tt=0;tt<=floord(time_size-2,8);tt++) {
        for (int xx=ceild(tt-3,4);xx<=min(floord(tt+9,4),floord(time_size+65,32));xx++) {
            for (int yy=ceild(tt-3,4);yy<=min(floord(tt+9,4),floord(time_size+65,32));yy++) {
                for (int time=max(max(max(1,8*tt),32*xx-67),32*yy-67);time<=min(min(min(time_size-2,8*tt+7),32*xx+30),32*yy+30);time++) {
                        int skew = 4*time;
                        int t0 = (time) % 8;
                        int t1 = (time + 1) % 8;
                        int t2 = (time - 1) % 8;
                        #pragma omp parallel
                        {
                            /* Flush denormal numbers to zero in hardware */
                            _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
                            _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
                            #pragma omp for schedule(static)
                            for (int x=max(128*xx,4*time+4);x<=min(4*time+271,128*xx+127);x++) {
                                for (int y=max(128*yy,4*time+4);y<=min(4*time+271,128*yy+127);y++) {
                                    #pragma ivdep
                                    #pragma omp simd
                                    for (int z=4*time+4;z<=4*time+271;z++) {
                                        float tcse0 = 3.04F*damp[x-skew][y-skew][z-skew];
                                        u[t1][x-skew][y-skew][z-skew] = ((tcse0 - 2*m[x-skew][y-skew][z-skew])*u[t2][x-skew][y-skew][z-skew] - 8.25142857142857e-5F*(u[t0][x-skew][y-skew][z-skew - 4] + u[t0][x-skew][y-skew][z-skew + 4] + u[t0][x-skew][y-skew - 4][z-skew] + u[t0][x-skew][y-skew + 4][z-skew] + u[t0][x-skew - 4][y-skew][z-skew] + u[t0][x-skew + 4][y-skew][z-skew]) + 1.17353650793651e-3F*(u[t0][x-skew][y-skew][z-skew - 3] + u[t0][x-skew][y-skew][z-skew + 3] + u[t0][x-skew][y-skew - 3][z-skew] + u[t0][x-skew][y-skew + 3][z-skew] + u[t0][x-skew - 3][y-skew][z-skew] + u[t0][x-skew + 3][y-skew][z-skew]) - 9.2416e-3F*(u[t0][x-skew][y-skew][z-skew - 2] + u[t0][x-skew][y-skew][z-skew + 2] + u[t0][x-skew][y-skew - 2][z-skew] + u[t0][x-skew][y-skew + 2][z-skew] + u[t0][x-skew - 2][y-skew][z-skew] + u[t0][x-skew + 2][y-skew][z-skew]) + 7.39328e-2F*(u[t0][x-skew][y-skew][z-skew - 1] + u[t0][x-skew][y-skew][z-skew + 1] + u[t0][x-skew][y-skew - 1][z-skew] + u[t0][x-skew][y-skew + 1][z-skew] + u[t0][x-skew - 1][y-skew][z-skew] + u[t0][x-skew + 1][y-skew][z-skew]) + 4*m[x-skew][y-skew][z-skew]*u[t0][x-skew][y-skew][z-skew] - 3.94693333333333e-1F*u[t0][x-skew][y-skew][z-skew])/(tcse0 + 2*m[x-skew][y-skew][z-skew]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
    }
    gettimeofday(&end_loop_x_0, NULL);
    timings->loop_x_0 += (double)(end_loop_x_0.tv_sec-start_loop_x_0.tv_sec)+(double)(end_loop_x_0.tv_usec-start_loop_x_0.tv_usec)/1000000;
    struct timeval start_loop_p_rec_1, end_loop_p_rec_1;
    gettimeofday(&start_loop_p_rec_1, NULL);
    for (int time = 1; time < time_size - 1; time += 1)
    {
        int t0 = (time) % 8;
        int t1 = (time + 1) % 8;
        int t2 = (time - 1) % 8;
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
    }
    gettimeofday(&end_loop_p_rec_1, NULL);
    timings->loop_p_rec_1 += (double)(end_loop_p_rec_1.tv_sec-start_loop_p_rec_1.tv_sec)+(double)(end_loop_p_rec_1.tv_usec-start_loop_p_rec_1.tv_usec)/1000000;
    return 0;
}

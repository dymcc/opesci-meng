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

struct profiler
{
    double loop_time_0;
} ;


int Forward(float *restrict damp_vec, float *restrict m_vec, float *restrict u_vec, const int time_size, const int x_size, const int y_size, const int z_size, struct profiler *timings)
{
    float (*restrict damp)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) damp_vec;
    float (*restrict m)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) m_vec;
    float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
    /* DLE: moved denormals flag */
    struct timeval start_loop_time_0, end_loop_time_0;
    gettimeofday(&start_loop_time_0, NULL);
    #pragma omp parallel
    {
        /* Flush denormal numbers to zero in hardware */
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
        #pragma omp for schedule(static)
        for (int xx=0;xx<=16;xx++) {
            for (int yy=0;yy<=16;yy++) {
                for (int zz=0;zz<=16;zz++) {
                    for (int time = 1; time < time_size - 1; time += 1)
                    {
                        int skew = 4*time;
                        for (int x=max(skew+4,skew+16*xx);x<=skew+16*xx+15;x++) {
                            for (int y=max(skew+4,skew+16*yy);y<=skew+16*yy+15;y++) {
                                #pragma ivdep
                                #pragma omp simd
                                for (int z=max(skew+4,skew+16*zz);z<=skew+16*zz+15;z++) {
                                    float tcse0 = 3.04F*damp[x-skew][y-skew][z-skew];
                                    u[time + 1][x-skew][y-skew][z-skew] = ((tcse0 - 2*m[x-skew][y-skew][z-skew])*u[time - 1][x-skew][y-skew][z-skew] - 8.25142857142857e-5F*(u[time][x-skew][y-skew][z-skew - 4] + u[time][x-skew][y-skew][z-skew + 4] + u[time][x-skew][y-skew - 4][z-skew] + u[time][x-skew][y-skew + 4][z-skew] + u[time][x-skew - 4][y-skew][z-skew] + u[time][x-skew + 4][y-skew][z-skew]) + 1.17353650793651e-3F*(u[time][x-skew][y-skew][z-skew - 3] + u[time][x-skew][y-skew][z-skew + 3] + u[time][x-skew][y-skew - 3][z-skew] + u[time][x-skew][y-skew + 3][z-skew] + u[time][x-skew - 3][y-skew][z-skew] + u[time][x-skew + 3][y-skew][z-skew]) - 9.2416e-3F*(u[time][x-skew][y-skew][z-skew - 2] + u[time][x-skew][y-skew][z-skew + 2] + u[time][x-skew][y-skew - 2][z-skew] + u[time][x-skew][y-skew + 2][z-skew] + u[time][x-skew - 2][y-skew][z-skew] + u[time][x-skew + 2][y-skew][z-skew]) + 7.39328e-2F*(u[time][x-skew][y-skew][z-skew - 1] + u[time][x-skew][y-skew][z-skew + 1] + u[time][x-skew][y-skew - 1][z-skew] + u[time][x-skew][y-skew + 1][z-skew] + u[time][x-skew - 1][y-skew][z-skew] + u[time][x-skew + 1][y-skew][z-skew]) + 4*m[x-skew][y-skew][z-skew]*u[time][x-skew][y-skew][z-skew] - 3.94693333333333e-1F*u[time][x-skew][y-skew][z-skew])/(tcse0 + 2*m[x-skew][y-skew][z-skew]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    gettimeofday(&end_loop_time_0, NULL);
    timings->loop_time_0 += (double)(end_loop_time_0.tv_sec-start_loop_time_0.tv_sec)+(double)(end_loop_time_0.tv_usec-start_loop_time_0.tv_usec)/1000000;
    return 0;
}

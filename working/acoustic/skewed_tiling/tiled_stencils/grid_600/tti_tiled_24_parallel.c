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
  double loop_x_0;
  double loop_x_1;
  double loop_p_src_2;
  double loop_p_rec_3;
} ;


int Kernel(float *restrict damp_vec, float *restrict delta_vec, float *restrict epsilon_vec, float *restrict m_vec, float *restrict phi_vec, float *restrict rec_vec, float *restrict rec_coords_vec, float *restrict src_vec, float *restrict src_coords_vec, float *restrict theta_vec, float *restrict u_vec, float *restrict v_vec, const int d_size, const int p_rec_size, const int p_src_size, const int t_size, const int time_size, const int x_size, const int y_size, const int z_size, struct profiler *timings)
{
  float (*restrict damp)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) damp_vec;
  float (*restrict delta)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) delta_vec;
  float (*restrict epsilon)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) epsilon_vec;
  float (*restrict m)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) m_vec;
  float (*restrict phi)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) phi_vec;
  float (*restrict rec)[p_rec_size] __attribute__((aligned(64))) = (float (*)[p_rec_size]) rec_vec;
  float (*restrict rec_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) rec_coords_vec;
  float (*restrict src)[p_src_size] __attribute__((aligned(64))) = (float (*)[p_src_size]) src_vec;
  float (*restrict src_coords)[d_size] __attribute__((aligned(64))) = (float (*)[d_size]) src_coords_vec;
  float (*restrict theta)[y_size][z_size] __attribute__((aligned(64))) = (float (*)[y_size][z_size]) theta_vec;
  float (*restrict u)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) u_vec;
  float (*restrict v)[x_size][y_size][z_size] __attribute__((aligned(64))) = (float (*)[x_size][y_size][z_size]) v_vec;
  float (*r0)[60][60];
  float (*r1)[60][60];
  float (*r2)[60][60];
  float (*r3)[60][60];
  posix_memalign((void**)&r0, 64, sizeof(float[60][60][60]));
  posix_memalign((void**)&r1, 64, sizeof(float[60][60][60]));
  posix_memalign((void**)&r2, 64, sizeof(float[60][60][60]));
  posix_memalign((void**)&r3, 64, sizeof(float[60][60][60]));
  /* DLE: moved denormals flag */
  struct timeval start_loop_x_0, end_loop_x_0;
  gettimeofday(&start_loop_x_0, NULL);
int  xx, yy, zz, x, y, z;
  for (xx=0;xx<=25;xx++) {
      for (yy=0;yy<=25;yy++) {
          for (zz=0;zz<=25;zz++) {
              for (x=max(1,24*xx);x<=min(617,24*xx+23);x++) {
                  for (y=max(1,24*yy);y<=min(617,24*yy+23);y++) {
                      for (z=max(1,24*zz);z<=min(617,24*zz+23);z++) {
                          r0[x][y][z] = sin(phi[x][y][z]);
                          r1[x][y][z] = sin(theta[x][y][z]);
                          r2[x][y][z] = cos(phi[x][y][z]);
                          r3[x][y][z] = cos(theta[x][y][z]);
                      }
                  }
              }
          }
      }
  }
  gettimeofday(&end_loop_x_0, NULL);
  timings->loop_x_0 += (double)(end_loop_x_0.tv_sec-start_loop_x_0.tv_sec)+(double)(end_loop_x_0.tv_usec-start_loop_x_0.tv_usec)/1000000;
  for (int time = 1; time < time_size - 1; time += 1)
  {
    int t0 = (time) % 3;
    int t1 = (time + 1) % 3;
    int t2 = (time - 1) % 3;
    struct timeval start_loop_x_1, end_loop_x_1;
    gettimeofday(&start_loop_x_1, NULL);
    int  xx, yy, zz, x, y, z;
    #pragma omp parallel
    {
      /* Flush denormal numbers to zero in hardware */
      _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
      _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
      #pragma omp for schedule(static)
      for (xx=0;xx<=25;xx++) {
          for (yy=0;yy<=25;yy++) {
              for (zz=0;zz<=25;zz++) {
                  for (x=max(3*time+4,3*time+24*xx);x<=min(3*time+615,3*time+24*xx+23);x++) {
                      for (y=max(3*time+4,3*time+24*yy);y<=min(3*time+615,3*time+24*yy+23);y++) {
                          #pragma ivdep
                          #pragma omp simd aligned(damp,m,u:32)
                          /*for (int z = 4; z < z_size - 4; z += 1)*/
                          for (z=max(3*time+4,3*time+24*zz);z<=min(3*time+615,3*time+24*zz+23);z++) {
                              {
                                  float ti1 = r1[x][y][z]*r2[x][y][z];
                                  float ti2 = r0[x][y][z]*r1[x][y][z];
                                  float ti4 = r1[x][y][z + 1]*r2[x][y][z + 1];
                                  float ti5 = r0[x][y][z + 1]*r1[x][y][z + 1];
                                  float ti7 = r1[x][y][z - 2]*r2[x][y][z - 2];
                                  float ti8 = r0[x][y][z - 2]*r1[x][y][z - 2];
                                  float ti10 = r1[x][y][z - 1]*r2[x][y][z - 1];
                                  float ti11 = r0[x][y][z - 1]*r1[x][y][z - 1];
                                  float ti13 = r0[x][y][z + 2]*r1[x][y][z + 2];
                                  float ti14 = r1[x][y][z + 2]*r2[x][y][z + 2];
                                  float ti16 = r1[x][y][z - 3]*r2[x][y][z - 3];
                                  float ti17 = r0[x][y][z - 3]*r1[x][y][z - 3];
                                  float ti19 = r0[x + 1][y][z]*r1[x + 1][y][z];
                                  float ti20 = r1[x + 1][y][z]*r2[x + 1][y][z];
                                  float ti22 = r1[x][y + 1][z]*r2[x][y + 1][z];
                                  float ti23 = r0[x][y + 1][z]*r1[x][y + 1][z];
                                  float ti25 = r0[x - 2][y][z]*r1[x - 2][y][z];
                                  float ti26 = r1[x - 2][y][z]*r2[x - 2][y][z];
                                  float ti28 = r1[x][y - 2][z]*r2[x][y - 2][z];
                                  float ti29 = r0[x][y - 2][z]*r1[x][y - 2][z];
                                  float ti31 = r0[x - 1][y][z]*r1[x - 1][y][z];
                                  float ti32 = r1[x - 1][y][z]*r2[x - 1][y][z];
                                  float ti34 = r1[x][y - 1][z]*r2[x][y - 1][z];
                                  float ti35 = r0[x][y - 1][z]*r1[x][y - 1][z];
                                  float ti37 = r0[x - 3][y][z]*r1[x - 3][y][z];
                                  float ti38 = r1[x - 3][y][z]*r2[x - 3][y][z];
                                  float ti40 = r0[x][y + 2][z]*r1[x][y + 2][z];
                                  float ti41 = r1[x][y + 2][z]*r2[x][y + 2][z];
                                  float ti43 = r1[x + 2][y][z]*r2[x + 2][y][z];
                                  float ti44 = r0[x + 2][y][z]*r1[x + 2][y][z];
                                  float ti46 = r1[x][y - 3][z]*r2[x][y - 3][z];
                                  float ti47 = r0[x][y - 3][z]*r1[x][y - 3][z];
                                  float ti51 = r2[x][y][z - 1]*r3[x][y][z - 1];
                                  float ti52 = r0[x][y][z - 1]*r3[x][y][z - 1];
                                  float ti64 = r0[x][y][z + 2]*r3[x][y][z + 2];
                                  float ti65 = r2[x][y][z + 2]*r3[x][y][z + 2];
                                  float ti67 = r2[x][y][z - 3]*r3[x][y][z - 3];
                                  float ti68 = r0[x][y][z - 3]*r3[x][y][z - 3];
                                  float ti72 = r2[x][y][z]*r3[x][y][z];
                                  float ti73 = r0[x][y][z]*r3[x][y][z];
                                  float ti77 = r2[x][y][z + 1]*r3[x][y][z + 1];
                                  float ti78 = r0[x][y][z + 1]*r3[x][y][z + 1];
                                  float ti82 = r2[x][y][z - 2]*r3[x][y][z - 2];
                                  float ti83 = r0[x][y][z - 2]*r3[x][y][z - 2];
                                  float ti89 = r0[x + 1][y][z]*r3[x + 1][y][z];
                                  float ti90 = r2[x + 1][y][z]*r3[x + 1][y][z];
                                  float ti92 = r2[x][y + 1][z]*r3[x][y + 1][z];
                                  float ti93 = r0[x][y + 1][z]*r3[x][y + 1][z];
                                  float ti95 = r0[x - 2][y][z]*r3[x - 2][y][z];
                                  float ti96 = r2[x - 2][y][z]*r3[x - 2][y][z];
                                  float ti98 = r2[x][y - 2][z]*r3[x][y - 2][z];
                                  float ti99 = r0[x][y - 2][z]*r3[x][y - 2][z];
                                  float ti101 = r0[x - 1][y][z]*r3[x - 1][y][z];
                                  float ti102 = r2[x - 1][y][z]*r3[x - 1][y][z];
                                  float ti104 = r2[x][y - 1][z]*r3[x][y - 1][z];
                                  float ti105 = r0[x][y - 1][z]*r3[x][y - 1][z];
                                  float ti107 = r0[x - 3][y][z]*r3[x - 3][y][z];
                                  float ti108 = r2[x - 3][y][z]*r3[x - 3][y][z];
                                  float ti110 = r0[x][y + 2][z]*r3[x][y + 2][z];
                                  float ti111 = r2[x][y + 2][z]*r3[x][y + 2][z];
                                  float ti113 = r2[x + 2][y][z]*r3[x + 2][y][z];
                                  float ti114 = r0[x + 2][y][z]*r3[x + 2][y][z];
                                  float ti116 = r2[x][y - 3][z]*r3[x][y - 3][z];
                                  float ti117 = r0[x][y - 3][z]*r3[x][y - 3][z];
                                  float tcse0 = -4.16666666666667e-2F*u[t0][x][y][z - 2];
                                  float tcse1 = -4.16666666666667e-2F*u[t0][x][y][z + 1];
                                  float tcse2 = 3.33333333333333e-2F*u[t0][x][y + 1][z + 1];
                                  float tcse3 = 7.5e-2F*u[t0][x + 1][y][z + 1];
                                  float tcse4 = -1.25e-2F*u[t0][x - 1][y][z - 1];
                                  float tcse5 = -3.33333333333333e-2F*u[t0][x][y - 1][z - 1];
                                  float tcse6 = 4.16666666666667e-3F*u[t0][x][y - 2][z - 2];
                                  float tcse7 = -2.5e-2F*u[t0][x + 2][y][z + 2];
                                  float tcse8 = -4.16666666666667e-3F*u[t0][x][y + 2][z + 2];
                                  float tcse9 = -4.16666666666667e-2F*u[t0][x][y][z - 3];
                                  float tcse10 = 4.16666666666667e-3F*u[t0][x - 2][y][z - 2];
                                  float tcse11 = 7.5e-2F*u[t0][x][y + 1][z + 1];
                                  float tcse12 = 3.33333333333333e-2F*u[t0][x + 1][y][z + 1];
                                  float tcse13 = -1.25e-2F*u[t0][x][y - 1][z - 1];
                                  float tcse14 = -4.16666666666667e-2F*u[t0][x][y][z - 1];
                                  float tcse15 = -3.33333333333333e-2F*u[t0][x - 1][y][z - 1];
                                  float tcse16 = 7.5e-2F*v[t0][x][y + 1][z + 1];
                                  float tcse17 = -1.25e-2F*v[t0][x][y - 1][z - 1];
                                  float tcse18 = -4.16666666666667e-2F*v[t0][x][y - 3][z];
                                  float tcse19 = -4.16666666666667e-2F*v[t0][x][y][z - 3];
                                  float tcse20 = -4.16666666666667e-2F*v[t0][x + 2][y][z];
                                  float tcse21 = -4.16666666666667e-2F*v[t0][x][y][z];
                                  float tcse22 = -1.25e-2F*v[t0][x][y][z];
                                  float tcse23 = -4.16666666666667e-2F*v[t0][x + 1][y][z];
                                  float tcse24 = 3.33333333333333e-2F*v[t0][x + 1][y][z + 1];
                                  float tcse25 = 3.33333333333333e-2F*v[t0][x + 1][y + 1][z];
                                  float tcse26 = -2.5e-2F*v[t0][x][y][z];
                                  float tcse27 = -4.16666666666667e-2F*v[t0][x - 2][y][z];
                                  float tcse28 = 4.16666666666667e-3F*v[t0][x - 2][y][z - 2];
                                  float tcse29 = 4.16666666666667e-3F*v[t0][x - 2][y - 2][z];
                                  float tcse30 = -4.16666666666667e-2F*v[t0][x][y - 2][z];
                                  float tcse31 = -4.16666666666667e-2F*v[t0][x][y][z - 2];
                                  float tcse32 = -4.16666666666667e-3F*v[t0][x][y][z];
                                  float tcse33 = 4.16666666666667e-3F*v[t0][x][y - 2][z - 2];
                                  float tcse34 = -4.16666666666667e-2F*v[t0][x][y + 1][z];
                                  float tcse35 = 7.5e-2F*v[t0][x + 1][y + 1][z];
                                  float tcse36 = -4.16666666666667e-2F*v[t0][x][y][z + 1];
                                  float tcse37 = 7.5e-2F*v[t0][x + 1][y][z + 1];
                                  float tcse38 = -3.33333333333333e-2F*v[t0][x][y][z];
                                  float tcse39 = 3.33333333333333e-2F*v[t0][x][y + 1][z + 1];
                                  float tcse40 = -4.16666666666667e-2F*v[t0][x - 1][y][z];
                                  float tcse41 = 7.5e-2F*v[t0][x][y][z];
                                  float tcse42 = -3.33333333333333e-2F*v[t0][x - 1][y][z - 1];
                                  float tcse43 = -3.33333333333333e-2F*v[t0][x - 1][y - 1][z];
                                  float tcse44 = -2.5e-2F*v[t0][x + 2][y][z + 2];
                                  float tcse45 = -2.5e-2F*v[t0][x + 2][y + 2][z];
                                  float tcse46 = 4.16666666666667e-3F*v[t0][x][y][z];
                                  float tcse47 = -4.16666666666667e-3F*v[t0][x][y + 2][z + 2];
                                  float tcse48 = -1.25e-2F*v[t0][x - 1][y - 1][z];
                                  float tcse49 = -4.16666666666667e-2F*v[t0][x][y - 1][z];
                                  float tcse50 = -1.25e-2F*v[t0][x - 1][y][z - 1];
                                  float tcse51 = -4.16666666666667e-2F*v[t0][x][y][z - 1];
                                  float tcse52 = 3.33333333333333e-2F*v[t0][x][y][z];
                                  float tcse53 = -3.33333333333333e-2F*v[t0][x][y - 1][z - 1];
                                  float tcse54 = -4.16666666666667e-2F*u[t0][x][y][z];
                                  float tcse55 = -1.25e-2F*u[t0][x][y][z];
                                  float tcse56 = -2.5e-2F*u[t0][x][y][z];
                                  float tcse57 = -4.16666666666667e-2F*u[t0][x - 2][y][z];
                                  float tcse58 = -4.16666666666667e-2F*u[t0][x][y - 3][z];
                                  float tcse59 = -4.16666666666667e-2F*u[t0][x][y - 2][z];
                                  float tcse60 = -4.16666666666667e-2F*u[t0][x + 1][y][z];
                                  float tcse61 = -4.16666666666667e-2F*u[t0][x][y - 1][z];
                                  float tcse62 = -4.16666666666667e-2F*u[t0][x][y + 1][z];
                                  float tcse63 = 7.5e-2F*u[t0][x + 1][y + 1][z];
                                  float tcse64 = 7.5e-2F*u[t0][x][y][z];
                                  float tcse65 = -2.5e-2F*u[t0][x + 2][y + 2][z];
                                  float tcse66 = -4.16666666666667e-2F*u[t0][x + 2][y][z];
                                  float tcse67 = -1.25e-2F*u[t0][x - 1][y - 1][z];
                                  float tcse68 = -4.16666666666667e-2F*u[t0][x - 1][y][z];
                                  float tcse69 = 3.33333333333333e-2F*u[t0][x + 1][y + 1][z];
                                  float tcse70 = 4.16666666666667e-3F*u[t0][x - 2][y - 2][z];
                                  float tcse71 = -4.16666666666667e-3F*u[t0][x][y][z];
                                  float tcse72 = -3.33333333333333e-2F*u[t0][x][y][z];
                                  float tcse73 = -3.33333333333333e-2F*u[t0][x - 1][y - 1][z];
                                  float tcse74 = 4.16666666666667e-3F*u[t0][x][y][z];
                                  float tcse75 = 3.33333333333333e-2F*u[t0][x][y][z];
                                  float tcse76 = 1.4161764656229F*damp[x][y][z];
                                  float tcse77 = tcse76 - 2.0F*m[x][y][z];
                                  float tcse78 = tcse76 + 2.0F*m[x][y][z];
                                  float tcse79 = tcse75 + 4.16666666666667e-3F*(u[t0][x - 3][y][z] - u[t0][x + 1][y][z]) - 3.33333333333333e-2F*u[t0][x - 2][y][z];
                                  float tcse80 = 4.16666666666667e-3F*(u[t0][x - 3][y - 2][z] - u[t0][x - 3][y + 2][z]) + 3.33333333333333e-2F*(-u[t0][x - 3][y - 1][z] + u[t0][x - 3][y + 1][z]);
                                  float tcse81 = tcse74 + 3.33333333333333e-2F*(-u[t0][x + 1][y][z] + u[t0][x + 3][y][z]) - 4.16666666666667e-3F*u[t0][x + 4][y][z];
                                  float tcse82 = tcse73 + 4.16666666666667e-3F*(u[t0][x - 1][y - 2][z] - u[t0][x - 1][y + 2][z]) + 3.33333333333333e-2F*u[t0][x - 1][y + 1][z];
                                  float tcse83 = tcse72 + 4.16666666666667e-3F*(u[t0][x][y - 1][z] - u[t0][x][y + 3][z]) + 3.33333333333333e-2F*u[t0][x][y + 2][z];
                                  float tcse84 = tcse75 + 4.16666666666667e-3F*(u[t0][x][y - 3][z] - u[t0][x][y + 1][z]) - 3.33333333333333e-2F*u[t0][x][y - 2][z];
                                  float tcse85 = tcse72 + 4.16666666666667e-3F*(u[t0][x - 1][y][z] - u[t0][x + 3][y][z]) + 3.33333333333333e-2F*u[t0][x + 2][y][z];
                                  float tcse86 = tcse71 + 3.33333333333333e-2F*(-u[t0][x][y - 3][z] + u[t0][x][y - 1][z]) + 4.16666666666667e-3F*u[t0][x][y - 4][z];
                                  float tcse87 = tcse74 + 3.33333333333333e-2F*(-u[t0][x][y + 1][z] + u[t0][x][y + 3][z]) - 4.16666666666667e-3F*u[t0][x][y + 4][z];
                                  float tcse88 = 4.16666666666667e-3F*(u[t0][x - 2][y - 3][z] - u[t0][x + 2][y - 3][z]) + 3.33333333333333e-2F*(-u[t0][x - 1][y - 3][z] + u[t0][x + 1][y - 3][z]);
                                  float tcse89 = tcse71 + 3.33333333333333e-2F*(-u[t0][x - 3][y][z] + u[t0][x - 1][y][z]) + 4.16666666666667e-3F*u[t0][x - 4][y][z];
                                  float tcse90 = tcse70 + 3.33333333333333e-2F*(-u[t0][x - 1][y - 2][z] + u[t0][x + 1][y - 2][z]) - 4.16666666666667e-3F*u[t0][x + 2][y - 2][z];
                                  float tcse91 = tcse70 + 3.33333333333333e-2F*(-u[t0][x - 2][y - 1][z] + u[t0][x - 2][y + 1][z]) - 4.16666666666667e-3F*u[t0][x - 2][y + 2][z];
                                  float tcse92 = tcse69 + 4.16666666666667e-3F*(u[t0][x - 2][y + 1][z] - u[t0][x + 2][y + 1][z]) - 3.33333333333333e-2F*u[t0][x - 1][y + 1][z];
                                  float tcse93 = tcse69 + 4.16666666666667e-3F*(u[t0][x + 1][y - 2][z] - u[t0][x + 1][y + 2][z]) - 3.33333333333333e-2F*u[t0][x + 1][y - 1][z];
                                  float tcse94 = 4.16666666666667e-3F*(u[t0][x][y - 2][z] - u[t0][x][y + 2][z]) + 3.33333333333333e-2F*(-u[t0][x][y - 1][z] + u[t0][x][y + 1][z]);
                                  float tcse95 = tcse73 + 4.16666666666667e-3F*(u[t0][x - 2][y - 1][z] - u[t0][x + 2][y - 1][z]) + 3.33333333333333e-2F*u[t0][x + 1][y - 1][z];
                                  float tcse96 = 4.16666666666667e-3F*(u[t0][x - 2][y][z] - u[t0][x + 2][y][z]) + 3.33333333333333e-2F*(-u[t0][x - 1][y][z] + u[t0][x + 1][y][z]);
                                  float tcse97 = tcse67 + tcse68 + 7.5e-2F*u[t0][x - 1][y + 1][z] - 2.5e-2F*u[t0][x - 1][y + 2][z] + 4.16666666666667e-3F*u[t0][x - 1][y + 3][z];
                                  float tcse98 = tcse74 - 1.25e-2F*u[t0][x - 4][y][z] - 4.16666666666667e-2F*u[t0][x - 3][y][z] + 7.5e-2F*u[t0][x - 2][y][z] - 2.5e-2F*u[t0][x - 1][y][z];
                                  float tcse99 = tcse65 + tcse66 - 1.25e-2F*u[t0][x + 2][y - 1][z] + 7.5e-2F*u[t0][x + 2][y + 1][z] + 4.16666666666667e-3F*u[t0][x + 2][y + 3][z];
                                  float tcse100 = tcse64 + tcse68 - 1.25e-2F*u[t0][x - 2][y][z] - 2.5e-2F*u[t0][x + 1][y][z] + 4.16666666666667e-3F*u[t0][x + 2][y][z];
                                  float tcse101 = tcse62 + tcse63 - 1.25e-2F*u[t0][x - 1][y + 1][z] - 2.5e-2F*u[t0][x + 2][y + 1][z] + 4.16666666666667e-3F*u[t0][x + 3][y + 1][z];
                                  float tcse102 = tcse61 + tcse67 + 7.5e-2F*u[t0][x + 1][y - 1][z] - 2.5e-2F*u[t0][x + 2][y - 1][z] + 4.16666666666667e-3F*u[t0][x + 3][y - 1][z];
                                  float tcse103 = tcse60 + tcse63 - 1.25e-2F*u[t0][x + 1][y - 1][z] - 2.5e-2F*u[t0][x + 1][y + 2][z] + 4.16666666666667e-3F*u[t0][x + 1][y + 3][z];
                                  float tcse104 = tcse59 - 1.25e-2F*u[t0][x - 1][y - 2][z] + 7.5e-2F*u[t0][x + 1][y - 2][z] - 2.5e-2F*u[t0][x + 2][y - 2][z] + 4.16666666666667e-3F*u[t0][x + 3][y - 2][z];
                                  float tcse105 = tcse65 - 4.16666666666667e-2F*u[t0][x][y + 2][z] - 1.25e-2F*u[t0][x - 1][y + 2][z] + 7.5e-2F*u[t0][x + 1][y + 2][z] + 4.16666666666667e-3F*u[t0][x + 3][y + 2][z];
                                  float tcse106 = tcse58 + tcse74 - 1.25e-2F*u[t0][x][y - 4][z] + 7.5e-2F*u[t0][x][y - 2][z] - 2.5e-2F*u[t0][x][y - 1][z];
                                  float tcse107 = tcse57 - 1.25e-2F*u[t0][x - 2][y - 1][z] + 7.5e-2F*u[t0][x - 2][y + 1][z] - 2.5e-2F*u[t0][x - 2][y + 2][z] + 4.16666666666667e-3F*u[t0][x - 2][y + 3][z];
                                  float tcse108 = tcse56 + tcse59 - 1.25e-2F*u[t0][x][y - 3][z] + 7.5e-2F*u[t0][x][y - 1][z] + 4.16666666666667e-3F*u[t0][x][y + 1][z];
                                  float tcse109 = tcse56 + tcse57 - 1.25e-2F*u[t0][x - 3][y][z] + 7.5e-2F*u[t0][x - 1][y][z] + 4.16666666666667e-3F*u[t0][x + 1][y][z];
                                  float tcse110 = tcse55 + tcse62 + 7.5e-2F*u[t0][x][y + 2][z] - 2.5e-2F*u[t0][x][y + 3][z] + 4.16666666666667e-3F*u[t0][x][y + 4][z];
                                  float tcse111 = tcse55 + tcse60 + 7.5e-2F*u[t0][x + 2][y][z] - 2.5e-2F*u[t0][x + 3][y][z] + 4.16666666666667e-3F*u[t0][x + 4][y][z];
                                  float tcse112 = tcse54 - 1.25e-2F*u[t0][x - 1][y][z] + 7.5e-2F*u[t0][x + 1][y][z] - 2.5e-2F*u[t0][x + 2][y][z] + 4.16666666666667e-3F*u[t0][x + 3][y][z];
                                  float tcse113 = tcse61 + tcse64 - 1.25e-2F*u[t0][x][y - 2][z] - 2.5e-2F*u[t0][x][y + 1][z] + 4.16666666666667e-3F*u[t0][x][y + 2][z];
                                  float tcse114 = tcse54 - 1.25e-2F*u[t0][x][y - 1][z] + 7.5e-2F*u[t0][x][y + 1][z] - 2.5e-2F*u[t0][x][y + 2][z] + 4.16666666666667e-3F*u[t0][x][y + 3][z];
                                  float tcse115 = ti10*(tcse50 + tcse51 + 7.5e-2F*v[t0][x + 1][y][z - 1] - 2.5e-2F*v[t0][x + 2][y][z - 1] + 4.16666666666667e-3F*v[t0][x + 3][y][z - 1]) + ti11*(tcse53 + 4.16666666666667e-3F*(v[t0][x][y - 2][z - 1] - v[t0][x][y + 2][z - 1]) + 3.33333333333333e-2F*v[t0][x][y + 1][z - 1]) + (tcse52 + 4.16666666666667e-3F*(v[t0][x][y][z - 3] - v[t0][x][y][z + 1]) - 3.33333333333333e-2F*v[t0][x][y][z - 2])*r3[x][y][z - 1];
                                  float tcse116 = ti34*(tcse48 + tcse49 + 7.5e-2F*v[t0][x + 1][y - 1][z] - 2.5e-2F*v[t0][x + 2][y - 1][z] + 4.16666666666667e-3F*v[t0][x + 3][y - 1][z]) + ti35*(tcse52 + 4.16666666666667e-3F*(v[t0][x][y - 3][z] - v[t0][x][y + 1][z]) - 3.33333333333333e-2F*v[t0][x][y - 2][z]) + (tcse53 + 4.16666666666667e-3F*(v[t0][x][y - 1][z - 2] - v[t0][x][y - 1][z + 2]) + 3.33333333333333e-2F*v[t0][x][y - 1][z + 1])*r3[x][y - 1][z];
                                  float tcse117 = ti40*(tcse46 + 3.33333333333333e-2F*(-v[t0][x][y + 1][z] + v[t0][x][y + 3][z]) - 4.16666666666667e-3F*v[t0][x][y + 4][z]) + ti41*(tcse45 - 4.16666666666667e-2F*v[t0][x][y + 2][z] - 1.25e-2F*v[t0][x - 1][y + 2][z] + 7.5e-2F*v[t0][x + 1][y + 2][z] + 4.16666666666667e-3F*v[t0][x + 3][y + 2][z]) + (tcse47 + 3.33333333333333e-2F*(-v[t0][x][y + 2][z - 1] + v[t0][x][y + 2][z + 1]) + 4.16666666666667e-3F*v[t0][x][y + 2][z - 2])*r3[x][y + 2][z];
                                  float tcse118 = ti13*(tcse47 + 3.33333333333333e-2F*(-v[t0][x][y - 1][z + 2] + v[t0][x][y + 1][z + 2]) + 4.16666666666667e-3F*v[t0][x][y - 2][z + 2]) + ti14*(tcse44 - 4.16666666666667e-2F*v[t0][x][y][z + 2] - 1.25e-2F*v[t0][x - 1][y][z + 2] + 7.5e-2F*v[t0][x + 1][y][z + 2] + 4.16666666666667e-3F*v[t0][x + 3][y][z + 2]) + (tcse46 + 3.33333333333333e-2F*(-v[t0][x][y][z + 1] + v[t0][x][y][z + 3]) - 4.16666666666667e-3F*v[t0][x][y][z + 4])*r3[x][y][z + 2];
                                  float tcse119 = ti37*(4.16666666666667e-3F*(v[t0][x - 3][y - 2][z] - v[t0][x - 3][y + 2][z]) + 3.33333333333333e-2F*(-v[t0][x - 3][y - 1][z] + v[t0][x - 3][y + 1][z])) + ti38*(tcse46 - 1.25e-2F*v[t0][x - 4][y][z] - 4.16666666666667e-2F*v[t0][x - 3][y][z] + 7.5e-2F*v[t0][x - 2][y][z] - 2.5e-2F*v[t0][x - 1][y][z]) + (4.16666666666667e-3F*(v[t0][x - 3][y][z - 2] - v[t0][x - 3][y][z + 2]) + 3.33333333333333e-2F*(-v[t0][x - 3][y][z - 1] + v[t0][x - 3][y][z + 1]))*r3[x - 3][y][z];
                                  float tcse120 = ti31*(tcse43 + 4.16666666666667e-3F*(v[t0][x - 1][y - 2][z] - v[t0][x - 1][y + 2][z]) + 3.33333333333333e-2F*v[t0][x - 1][y + 1][z]) + ti32*(tcse40 + tcse41 - 1.25e-2F*v[t0][x - 2][y][z] - 2.5e-2F*v[t0][x + 1][y][z] + 4.16666666666667e-3F*v[t0][x + 2][y][z]) + (tcse42 + 4.16666666666667e-3F*(v[t0][x - 1][y][z - 2] - v[t0][x - 1][y][z + 2]) + 3.33333333333333e-2F*v[t0][x - 1][y][z + 1])*r3[x - 1][y][z];
                                  float tcse121 = ti4*(tcse36 + tcse37 - 1.25e-2F*v[t0][x - 1][y][z + 1] - 2.5e-2F*v[t0][x + 2][y][z + 1] + 4.16666666666667e-3F*v[t0][x + 3][y][z + 1]) + ti5*(tcse39 + 4.16666666666667e-3F*(v[t0][x][y - 2][z + 1] - v[t0][x][y + 2][z + 1]) - 3.33333333333333e-2F*v[t0][x][y - 1][z + 1]) + (tcse38 + 4.16666666666667e-3F*(v[t0][x][y][z - 1] - v[t0][x][y][z + 3]) + 3.33333333333333e-2F*v[t0][x][y][z + 2])*r3[x][y][z + 1];
                                  float tcse122 = ti22*(tcse34 + tcse35 - 1.25e-2F*v[t0][x - 1][y + 1][z] - 2.5e-2F*v[t0][x + 2][y + 1][z] + 4.16666666666667e-3F*v[t0][x + 3][y + 1][z]) + ti23*(tcse38 + 4.16666666666667e-3F*(v[t0][x][y - 1][z] - v[t0][x][y + 3][z]) + 3.33333333333333e-2F*v[t0][x][y + 2][z]) + (tcse39 + 4.16666666666667e-3F*(v[t0][x][y + 1][z - 2] - v[t0][x][y + 1][z + 2]) - 3.33333333333333e-2F*v[t0][x][y + 1][z - 1])*r3[x][y + 1][z];
                                  float tcse123 = ti7*(tcse31 - 1.25e-2F*v[t0][x - 1][y][z - 2] + 7.5e-2F*v[t0][x + 1][y][z - 2] - 2.5e-2F*v[t0][x + 2][y][z - 2] + 4.16666666666667e-3F*v[t0][x + 3][y][z - 2]) + ti8*(tcse33 + 3.33333333333333e-2F*(-v[t0][x][y - 1][z - 2] + v[t0][x][y + 1][z - 2]) - 4.16666666666667e-3F*v[t0][x][y + 2][z - 2]) + (tcse32 + 3.33333333333333e-2F*(-v[t0][x][y][z - 3] + v[t0][x][y][z - 1]) + 4.16666666666667e-3F*v[t0][x][y][z - 4])*r3[x][y][z - 2];
                                  float tcse124 = ti28*(tcse30 - 1.25e-2F*v[t0][x - 1][y - 2][z] + 7.5e-2F*v[t0][x + 1][y - 2][z] - 2.5e-2F*v[t0][x + 2][y - 2][z] + 4.16666666666667e-3F*v[t0][x + 3][y - 2][z]) + ti29*(tcse32 + 3.33333333333333e-2F*(-v[t0][x][y - 3][z] + v[t0][x][y - 1][z]) + 4.16666666666667e-3F*v[t0][x][y - 4][z]) + (tcse33 + 3.33333333333333e-2F*(-v[t0][x][y - 2][z - 1] + v[t0][x][y - 2][z + 1]) - 4.16666666666667e-3F*v[t0][x][y - 2][z + 2])*r3[x][y - 2][z];
                                  float tcse125 = ti25*(tcse29 + 3.33333333333333e-2F*(-v[t0][x - 2][y - 1][z] + v[t0][x - 2][y + 1][z]) - 4.16666666666667e-3F*v[t0][x - 2][y + 2][z]) + ti26*(tcse26 + tcse27 - 1.25e-2F*v[t0][x - 3][y][z] + 7.5e-2F*v[t0][x - 1][y][z] + 4.16666666666667e-3F*v[t0][x + 1][y][z]) + (tcse28 + 3.33333333333333e-2F*(-v[t0][x - 2][y][z - 1] + v[t0][x - 2][y][z + 1]) - 4.16666666666667e-3F*v[t0][x - 2][y][z + 2])*r3[x - 2][y][z];
                                  float tcse126 = ti19*(tcse25 + 4.16666666666667e-3F*(v[t0][x + 1][y - 2][z] - v[t0][x + 1][y + 2][z]) - 3.33333333333333e-2F*v[t0][x + 1][y - 1][z]) + ti20*(tcse22 + tcse23 + 7.5e-2F*v[t0][x + 2][y][z] - 2.5e-2F*v[t0][x + 3][y][z] + 4.16666666666667e-3F*v[t0][x + 4][y][z]) + (tcse24 + 4.16666666666667e-3F*(v[t0][x + 1][y][z - 2] - v[t0][x + 1][y][z + 2]) - 3.33333333333333e-2F*v[t0][x + 1][y][z - 1])*r3[x + 1][y][z];
                                  float tcse127 = ti1*(tcse21 - 1.25e-2F*v[t0][x - 1][y][z] + 7.5e-2F*v[t0][x + 1][y][z] - 2.5e-2F*v[t0][x + 2][y][z] + 4.16666666666667e-3F*v[t0][x + 3][y][z]) + ti2*(4.16666666666667e-3F*(v[t0][x][y - 2][z] - v[t0][x][y + 2][z]) + 3.33333333333333e-2F*(-v[t0][x][y - 1][z] + v[t0][x][y + 1][z])) + (4.16666666666667e-3F*(v[t0][x][y][z - 2] - v[t0][x][y][z + 2]) + 3.33333333333333e-2F*(-v[t0][x][y][z - 1] + v[t0][x][y][z + 1]))*r3[x][y][z];
                                  float tcse128 = ti31*(tcse40 + tcse48 + 7.5e-2F*v[t0][x - 1][y + 1][z] - 2.5e-2F*v[t0][x - 1][y + 2][z] + 4.16666666666667e-3F*v[t0][x - 1][y + 3][z]) + ti32*(tcse52 + 4.16666666666667e-3F*(v[t0][x - 3][y][z] - v[t0][x + 1][y][z]) - 3.33333333333333e-2F*v[t0][x - 2][y][z]) + (tcse40 + tcse50 + 7.5e-2F*v[t0][x - 1][y][z + 1] - 2.5e-2F*v[t0][x - 1][y][z + 2] + 4.16666666666667e-3F*v[t0][x - 1][y][z + 3])*r3[x - 1][y][z];
                                  float tcse129 = ti43*(tcse46 + 3.33333333333333e-2F*(-v[t0][x + 1][y][z] + v[t0][x + 3][y][z]) - 4.16666666666667e-3F*v[t0][x + 4][y][z]) + ti44*(tcse20 + tcse45 - 1.25e-2F*v[t0][x + 2][y - 1][z] + 7.5e-2F*v[t0][x + 2][y + 1][z] + 4.16666666666667e-3F*v[t0][x + 2][y + 3][z]) + (tcse20 + tcse44 - 1.25e-2F*v[t0][x + 2][y][z - 1] + 7.5e-2F*v[t0][x + 2][y][z + 1] + 4.16666666666667e-3F*v[t0][x + 2][y][z + 3])*r3[x + 2][y][z];
                                  float tcse130 = ti16*(4.16666666666667e-3F*(v[t0][x - 2][y][z - 3] - v[t0][x + 2][y][z - 3]) + 3.33333333333333e-2F*(-v[t0][x - 1][y][z - 3] + v[t0][x + 1][y][z - 3])) + ti17*(tcse19 - 1.25e-2F*v[t0][x][y - 1][z - 3] + 7.5e-2F*v[t0][x][y + 1][z - 3] - 2.5e-2F*v[t0][x][y + 2][z - 3] + 4.16666666666667e-3F*v[t0][x][y + 3][z - 3]) + (tcse19 + tcse46 - 1.25e-2F*v[t0][x][y][z - 4] + 7.5e-2F*v[t0][x][y][z - 2] - 2.5e-2F*v[t0][x][y][z - 1])*r3[x][y][z - 3];
                                  float tcse131 = ti46*(4.16666666666667e-3F*(v[t0][x - 2][y - 3][z] - v[t0][x + 2][y - 3][z]) + 3.33333333333333e-2F*(-v[t0][x - 1][y - 3][z] + v[t0][x + 1][y - 3][z])) + ti47*(tcse18 + tcse46 - 1.25e-2F*v[t0][x][y - 4][z] + 7.5e-2F*v[t0][x][y - 2][z] - 2.5e-2F*v[t0][x][y - 1][z]) + (tcse18 - 1.25e-2F*v[t0][x][y - 3][z - 1] + 7.5e-2F*v[t0][x][y - 3][z + 1] - 2.5e-2F*v[t0][x][y - 3][z + 2] + 4.16666666666667e-3F*v[t0][x][y - 3][z + 3])*r3[x][y - 3][z];
                                  float tcse132 = ti10*(tcse42 + 4.16666666666667e-3F*(v[t0][x - 2][y][z - 1] - v[t0][x + 2][y][z - 1]) + 3.33333333333333e-2F*v[t0][x + 1][y][z - 1]) + ti11*(tcse17 + tcse51 + 7.5e-2F*v[t0][x][y + 1][z - 1] - 2.5e-2F*v[t0][x][y + 2][z - 1] + 4.16666666666667e-3F*v[t0][x][y + 3][z - 1]) + (tcse41 + tcse51 - 1.25e-2F*v[t0][x][y][z - 2] - 2.5e-2F*v[t0][x][y][z + 1] + 4.16666666666667e-3F*v[t0][x][y][z + 2])*r3[x][y][z - 1];
                                  float tcse133 = ti34*(tcse43 + 4.16666666666667e-3F*(v[t0][x - 2][y - 1][z] - v[t0][x + 2][y - 1][z]) + 3.33333333333333e-2F*v[t0][x + 1][y - 1][z]) + ti35*(tcse41 + tcse49 - 1.25e-2F*v[t0][x][y - 2][z] - 2.5e-2F*v[t0][x][y + 1][z] + 4.16666666666667e-3F*v[t0][x][y + 2][z]) + (tcse17 + tcse49 + 7.5e-2F*v[t0][x][y - 1][z + 1] - 2.5e-2F*v[t0][x][y - 1][z + 2] + 4.16666666666667e-3F*v[t0][x][y - 1][z + 3])*r3[x][y - 1][z];
                                  float tcse134 = ti19*(tcse23 + tcse35 - 1.25e-2F*v[t0][x + 1][y - 1][z] - 2.5e-2F*v[t0][x + 1][y + 2][z] + 4.16666666666667e-3F*v[t0][x + 1][y + 3][z]) + ti20*(tcse38 + 4.16666666666667e-3F*(v[t0][x - 1][y][z] - v[t0][x + 3][y][z]) + 3.33333333333333e-2F*v[t0][x + 2][y][z]) + (tcse23 + tcse37 - 1.25e-2F*v[t0][x + 1][y][z - 1] - 2.5e-2F*v[t0][x + 1][y][z + 2] + 4.16666666666667e-3F*v[t0][x + 1][y][z + 3])*r3[x + 1][y][z];
                                  float tcse135 = ti25*(tcse27 - 1.25e-2F*v[t0][x - 2][y - 1][z] + 7.5e-2F*v[t0][x - 2][y + 1][z] - 2.5e-2F*v[t0][x - 2][y + 2][z] + 4.16666666666667e-3F*v[t0][x - 2][y + 3][z]) + ti26*(tcse32 + 3.33333333333333e-2F*(-v[t0][x - 3][y][z] + v[t0][x - 1][y][z]) + 4.16666666666667e-3F*v[t0][x - 4][y][z]) + (tcse27 - 1.25e-2F*v[t0][x - 2][y][z - 1] + 7.5e-2F*v[t0][x - 2][y][z + 1] - 2.5e-2F*v[t0][x - 2][y][z + 2] + 4.16666666666667e-3F*v[t0][x - 2][y][z + 3])*r3[x - 2][y][z];
                                  float tcse136 = ti7*(tcse28 + 3.33333333333333e-2F*(-v[t0][x - 1][y][z - 2] + v[t0][x + 1][y][z - 2]) - 4.16666666666667e-3F*v[t0][x + 2][y][z - 2]) + ti8*(tcse31 - 1.25e-2F*v[t0][x][y - 1][z - 2] + 7.5e-2F*v[t0][x][y + 1][z - 2] - 2.5e-2F*v[t0][x][y + 2][z - 2] + 4.16666666666667e-3F*v[t0][x][y + 3][z - 2]) + (tcse26 + tcse31 - 1.25e-2F*v[t0][x][y][z - 3] + 7.5e-2F*v[t0][x][y][z - 1] + 4.16666666666667e-3F*v[t0][x][y][z + 1])*r3[x][y][z - 2];
                                  float tcse137 = ti28*(tcse29 + 3.33333333333333e-2F*(-v[t0][x - 1][y - 2][z] + v[t0][x + 1][y - 2][z]) - 4.16666666666667e-3F*v[t0][x + 2][y - 2][z]) + ti29*(tcse26 + tcse30 - 1.25e-2F*v[t0][x][y - 3][z] + 7.5e-2F*v[t0][x][y - 1][z] + 4.16666666666667e-3F*v[t0][x][y + 1][z]) + (tcse30 - 1.25e-2F*v[t0][x][y - 2][z - 1] + 7.5e-2F*v[t0][x][y - 2][z + 1] - 2.5e-2F*v[t0][x][y - 2][z + 2] + 4.16666666666667e-3F*v[t0][x][y - 2][z + 3])*r3[x][y - 2][z];
                                  float tcse138 = ti4*(tcse24 + 4.16666666666667e-3F*(v[t0][x - 2][y][z + 1] - v[t0][x + 2][y][z + 1]) - 3.33333333333333e-2F*v[t0][x - 1][y][z + 1]) + ti5*(tcse16 + tcse36 - 1.25e-2F*v[t0][x][y - 1][z + 1] - 2.5e-2F*v[t0][x][y + 2][z + 1] + 4.16666666666667e-3F*v[t0][x][y + 3][z + 1]) + (tcse22 + tcse36 + 7.5e-2F*v[t0][x][y][z + 2] - 2.5e-2F*v[t0][x][y][z + 3] + 4.16666666666667e-3F*v[t0][x][y][z + 4])*r3[x][y][z + 1];
                                  float tcse139 = ti22*(tcse25 + 4.16666666666667e-3F*(v[t0][x - 2][y + 1][z] - v[t0][x + 2][y + 1][z]) - 3.33333333333333e-2F*v[t0][x - 1][y + 1][z]) + ti23*(tcse22 + tcse34 + 7.5e-2F*v[t0][x][y + 2][z] - 2.5e-2F*v[t0][x][y + 3][z] + 4.16666666666667e-3F*v[t0][x][y + 4][z]) + (tcse16 + tcse34 - 1.25e-2F*v[t0][x][y + 1][z - 1] - 2.5e-2F*v[t0][x][y + 1][z + 2] + 4.16666666666667e-3F*v[t0][x][y + 1][z + 3])*r3[x][y + 1][z];
                                  float tcse140 = ti1*(4.16666666666667e-3F*(v[t0][x - 2][y][z] - v[t0][x + 2][y][z]) + 3.33333333333333e-2F*(-v[t0][x - 1][y][z] + v[t0][x + 1][y][z])) + ti2*(tcse21 - 1.25e-2F*v[t0][x][y - 1][z] + 7.5e-2F*v[t0][x][y + 1][z] - 2.5e-2F*v[t0][x][y + 2][z] + 4.16666666666667e-3F*v[t0][x][y + 3][z]) + (tcse21 - 1.25e-2F*v[t0][x][y][z - 1] + 7.5e-2F*v[t0][x][y][z + 1] - 2.5e-2F*v[t0][x][y][z + 2] + 4.16666666666667e-3F*v[t0][x][y][z + 3])*r3[x][y][z];
                                  float tcse141 = tcse114*ti73 + tcse96*ti72 - (tcse54 - 1.25e-2F*u[t0][x][y][z - 1] + 7.5e-2F*u[t0][x][y][z + 1] - 2.5e-2F*u[t0][x][y][z + 2] + 4.16666666666667e-3F*u[t0][x][y][z + 3])*r1[x][y][z];
                                  float tcse142 = 2.08333333333333e-2F*(tcse141*ti73 - tcse141*r1[x][y][z] + ti72*(tcse112*ti72 + tcse94*ti73 - (4.16666666666667e-3F*(u[t0][x][y][z - 2] - u[t0][x][y][z + 2]) + 3.33333333333333e-2F*(-u[t0][x][y][z - 1] + u[t0][x][y][z + 1]))*r1[x][y][z]) + (tcse112*r0[x][y][z] - tcse94*r2[x][y][z])*r0[x][y][z] - (-tcse114*r2[x][y][z] + tcse96*r0[x][y][z])*r2[x][y][z]) + 3.75e-2F*(-ti102*(tcse100*ti102 + tcse82*ti101 - (tcse15 + 4.16666666666667e-3F*(u[t0][x - 1][y][z - 2] - u[t0][x - 1][y][z + 2]) + 3.33333333333333e-2F*u[t0][x - 1][y][z + 1])*r1[x - 1][y][z]) - ti105*(tcse113*ti105 + tcse95*ti104 - (tcse13 + tcse61 + 7.5e-2F*u[t0][x][y - 1][z + 1] - 2.5e-2F*u[t0][x][y - 1][z + 2] + 4.16666666666667e-3F*u[t0][x][y - 1][z + 3])*r1[x][y - 1][z]) - (tcse100*r0[x - 1][y][z] - tcse82*r2[x - 1][y][z])*r0[x - 1][y][z] + (-tcse113*r2[x][y - 1][z] + tcse95*r0[x][y - 1][z])*r2[x][y - 1][z] + (ti51*(tcse15 + 4.16666666666667e-3F*(u[t0][x - 2][y][z - 1] - u[t0][x + 2][y][z - 1]) + 3.33333333333333e-2F*u[t0][x + 1][y][z - 1]) + ti52*(tcse13 + tcse14 + 7.5e-2F*u[t0][x][y + 1][z - 1] - 2.5e-2F*u[t0][x][y + 2][z - 1] + 4.16666666666667e-3F*u[t0][x][y + 3][z - 1]) - (tcse14 + tcse64 - 1.25e-2F*u[t0][x][y][z - 2] - 2.5e-2F*u[t0][x][y][z + 1] + 4.16666666666667e-3F*u[t0][x][y][z + 2])*r1[x][y][z - 1])*r1[x][y][z - 1]) + 6.25e-3F*(ti90*(tcse111*ti90 + tcse93*ti89 - (tcse12 + 4.16666666666667e-3F*(u[t0][x + 1][y][z - 2] - u[t0][x + 1][y][z + 2]) - 3.33333333333333e-2F*u[t0][x + 1][y][z - 1])*r1[x + 1][y][z]) + ti93*(tcse110*ti93 + tcse92*ti92 - (tcse11 + tcse62 - 1.25e-2F*u[t0][x][y + 1][z - 1] - 2.5e-2F*u[t0][x][y + 1][z + 2] + 4.16666666666667e-3F*u[t0][x][y + 1][z + 3])*r1[x][y + 1][z]) - (-tcse110*r2[x][y + 1][z] + tcse92*r0[x][y + 1][z])*r2[x][y + 1][z] + (tcse111*r0[x + 1][y][z] - tcse93*r2[x + 1][y][z])*r0[x + 1][y][z] - (ti77*(tcse12 + 4.16666666666667e-3F*(u[t0][x - 2][y][z + 1] - u[t0][x + 2][y][z + 1]) - 3.33333333333333e-2F*u[t0][x - 1][y][z + 1]) + ti78*(tcse1 + tcse11 - 1.25e-2F*u[t0][x][y - 1][z + 1] - 2.5e-2F*u[t0][x][y + 2][z + 1] + 4.16666666666667e-3F*u[t0][x][y + 3][z + 1]) - (tcse1 + tcse55 + 7.5e-2F*u[t0][x][y][z + 2] - 2.5e-2F*u[t0][x][y][z + 3] + 4.16666666666667e-3F*u[t0][x][y][z + 4])*r1[x][y][z + 1])*r1[x][y][z + 1]) + 1.25e-2F*(ti96*(tcse109*ti96 + tcse91*ti95 - (tcse10 + 3.33333333333333e-2F*(-u[t0][x - 2][y][z - 1] + u[t0][x - 2][y][z + 1]) - 4.16666666666667e-3F*u[t0][x - 2][y][z + 2])*r1[x - 2][y][z]) + ti99*(tcse108*ti99 + tcse90*ti98 - (tcse59 - 1.25e-2F*u[t0][x][y - 2][z - 1] + 7.5e-2F*u[t0][x][y - 2][z + 1] - 2.5e-2F*u[t0][x][y - 2][z + 2] + 4.16666666666667e-3F*u[t0][x][y - 2][z + 3])*r1[x][y - 2][z]) - (-tcse108*r2[x][y - 2][z] + tcse90*r0[x][y - 2][z])*r2[x][y - 2][z] + (tcse109*r0[x - 2][y][z] - tcse91*r2[x - 2][y][z])*r0[x - 2][y][z] - (ti82*(tcse10 + 3.33333333333333e-2F*(-u[t0][x - 1][y][z - 2] + u[t0][x + 1][y][z - 2]) - 4.16666666666667e-3F*u[t0][x + 2][y][z - 2]) + ti83*(tcse0 - 1.25e-2F*u[t0][x][y - 1][z - 2] + 7.5e-2F*u[t0][x][y + 1][z - 2] - 2.5e-2F*u[t0][x][y + 2][z - 2] + 4.16666666666667e-3F*u[t0][x][y + 3][z - 2]) - (tcse0 + tcse56 - 1.25e-2F*u[t0][x][y][z - 3] + 7.5e-2F*u[t0][x][y][z - 1] + 4.16666666666667e-3F*u[t0][x][y][z + 1])*r1[x][y][z - 2])*r1[x][y][z - 2]) + 1.66666666666667e-2F*(-ti102*(tcse79*ti102 + tcse97*ti101 - (tcse4 + tcse68 + 7.5e-2F*u[t0][x - 1][y][z + 1] - 2.5e-2F*u[t0][x - 1][y][z + 2] + 4.16666666666667e-3F*u[t0][x - 1][y][z + 3])*r1[x - 1][y][z]) - ti105*(tcse102*ti104 + tcse84*ti105 - (tcse5 + 4.16666666666667e-3F*(u[t0][x][y - 1][z - 2] - u[t0][x][y - 1][z + 2]) + 3.33333333333333e-2F*u[t0][x][y - 1][z + 1])*r1[x][y - 1][z]) + ti90*(tcse103*ti89 + tcse85*ti90 - (tcse3 + tcse60 - 1.25e-2F*u[t0][x + 1][y][z - 1] - 2.5e-2F*u[t0][x + 1][y][z + 2] + 4.16666666666667e-3F*u[t0][x + 1][y][z + 3])*r1[x + 1][y][z]) + ti93*(tcse101*ti92 + tcse83*ti93 - (tcse2 + 4.16666666666667e-3F*(u[t0][x][y + 1][z - 2] - u[t0][x][y + 1][z + 2]) - 3.33333333333333e-2F*u[t0][x][y + 1][z - 1])*r1[x][y + 1][z]) - (tcse101*r0[x][y + 1][z] - tcse83*r2[x][y + 1][z])*r2[x][y + 1][z] + (tcse102*r0[x][y - 1][z] - tcse84*r2[x][y - 1][z])*r2[x][y - 1][z] + (-tcse103*r2[x + 1][y][z] + tcse85*r0[x + 1][y][z])*r0[x + 1][y][z] - (tcse79*r0[x - 1][y][z] - tcse97*r2[x - 1][y][z])*r0[x - 1][y][z] + (ti51*(tcse14 + tcse4 + 7.5e-2F*u[t0][x + 1][y][z - 1] - 2.5e-2F*u[t0][x + 2][y][z - 1] + 4.16666666666667e-3F*u[t0][x + 3][y][z - 1]) + ti52*(tcse5 + 4.16666666666667e-3F*(u[t0][x][y - 2][z - 1] - u[t0][x][y + 2][z - 1]) + 3.33333333333333e-2F*u[t0][x][y + 1][z - 1]) - (tcse75 + 4.16666666666667e-3F*(u[t0][x][y][z - 3] - u[t0][x][y][z + 1]) - 3.33333333333333e-2F*u[t0][x][y][z - 2])*r1[x][y][z - 1])*r1[x][y][z - 1] - (ti77*(tcse1 + tcse3 - 1.25e-2F*u[t0][x - 1][y][z + 1] - 2.5e-2F*u[t0][x + 2][y][z + 1] + 4.16666666666667e-3F*u[t0][x + 3][y][z + 1]) + ti78*(tcse2 + 4.16666666666667e-3F*(u[t0][x][y - 2][z + 1] - u[t0][x][y + 2][z + 1]) - 3.33333333333333e-2F*u[t0][x][y - 1][z + 1]) - (tcse72 + 4.16666666666667e-3F*(u[t0][x][y][z - 1] - u[t0][x][y][z + 3]) + 3.33333333333333e-2F*u[t0][x][y][z + 2])*r1[x][y][z + 1])*r1[x][y][z + 1]) + 2.08333333333333e-3F*(-ti108*(tcse80*ti107 + tcse98*ti108 - (4.16666666666667e-3F*(u[t0][x - 3][y][z - 2] - u[t0][x - 3][y][z + 2]) + 3.33333333333333e-2F*(-u[t0][x - 3][y][z - 1] + u[t0][x - 3][y][z + 1]))*r1[x - 3][y][z]) - ti110*(tcse105*ti111 + tcse87*ti110 - (tcse8 + 3.33333333333333e-2F*(-u[t0][x][y + 2][z - 1] + u[t0][x][y + 2][z + 1]) + 4.16666666666667e-3F*u[t0][x][y + 2][z - 2])*r1[x][y + 2][z]) - ti113*(tcse81*ti113 + tcse99*ti114 - (tcse66 + tcse7 - 1.25e-2F*u[t0][x + 2][y][z - 1] + 7.5e-2F*u[t0][x + 2][y][z + 1] + 4.16666666666667e-3F*u[t0][x + 2][y][z + 3])*r1[x + 2][y][z]) - ti117*(tcse106*ti117 + tcse88*ti116 - (tcse58 - 1.25e-2F*u[t0][x][y - 3][z - 1] + 7.5e-2F*u[t0][x][y - 3][z + 1] - 2.5e-2F*u[t0][x][y - 3][z + 2] + 4.16666666666667e-3F*u[t0][x][y - 3][z + 3])*r1[x][y - 3][z]) + ti96*(tcse107*ti95 + tcse89*ti96 - (tcse57 - 1.25e-2F*u[t0][x - 2][y][z - 1] + 7.5e-2F*u[t0][x - 2][y][z + 1] - 2.5e-2F*u[t0][x - 2][y][z + 2] + 4.16666666666667e-3F*u[t0][x - 2][y][z + 3])*r1[x - 2][y][z]) + ti99*(tcse104*ti98 + tcse86*ti99 - (tcse6 + 3.33333333333333e-2F*(-u[t0][x][y - 2][z - 1] + u[t0][x][y - 2][z + 1]) - 4.16666666666667e-3F*u[t0][x][y - 2][z + 2])*r1[x][y - 2][z]) - (tcse104*r0[x][y - 2][z] - tcse86*r2[x][y - 2][z])*r2[x][y - 2][z] + (tcse105*r0[x][y + 2][z] - tcse87*r2[x][y + 2][z])*r2[x][y + 2][z] + (-tcse106*r2[x][y - 3][z] + tcse88*r0[x][y - 3][z])*r2[x][y - 3][z] + (-tcse107*r2[x - 2][y][z] + tcse89*r0[x - 2][y][z])*r0[x - 2][y][z] - (-tcse80*r2[x - 3][y][z] + tcse98*r0[x - 3][y][z])*r0[x - 3][y][z] - (tcse81*r0[x + 2][y][z] - tcse99*r2[x + 2][y][z])*r0[x + 2][y][z] + (ti64*(tcse8 + 3.33333333333333e-2F*(-u[t0][x][y - 1][z + 2] + u[t0][x][y + 1][z + 2]) + 4.16666666666667e-3F*u[t0][x][y - 2][z + 2]) + ti65*(tcse7 - 4.16666666666667e-2F*u[t0][x][y][z + 2] - 1.25e-2F*u[t0][x - 1][y][z + 2] + 7.5e-2F*u[t0][x + 1][y][z + 2] + 4.16666666666667e-3F*u[t0][x + 3][y][z + 2]) - (tcse74 + 3.33333333333333e-2F*(-u[t0][x][y][z + 1] + u[t0][x][y][z + 3]) - 4.16666666666667e-3F*u[t0][x][y][z + 4])*r1[x][y][z + 2])*r1[x][y][z + 2] + (ti67*(4.16666666666667e-3F*(u[t0][x - 2][y][z - 3] - u[t0][x + 2][y][z - 3]) + 3.33333333333333e-2F*(-u[t0][x - 1][y][z - 3] + u[t0][x + 1][y][z - 3])) + ti68*(tcse9 - 1.25e-2F*u[t0][x][y - 1][z - 3] + 7.5e-2F*u[t0][x][y + 1][z - 3] - 2.5e-2F*u[t0][x][y + 2][z - 3] + 4.16666666666667e-3F*u[t0][x][y + 3][z - 3]) - (tcse74 + tcse9 - 1.25e-2F*u[t0][x][y][z - 4] + 7.5e-2F*u[t0][x][y][z - 2] - 2.5e-2F*u[t0][x][y][z - 1])*r1[x][y][z - 3])*r1[x][y][z - 3] - (ti82*(tcse0 - 1.25e-2F*u[t0][x - 1][y][z - 2] + 7.5e-2F*u[t0][x + 1][y][z - 2] - 2.5e-2F*u[t0][x + 2][y][z - 2] + 4.16666666666667e-3F*u[t0][x + 3][y][z - 2]) + ti83*(tcse6 + 3.33333333333333e-2F*(-u[t0][x][y - 1][z - 2] + u[t0][x][y + 1][z - 2]) - 4.16666666666667e-3F*u[t0][x][y + 2][z - 2]) - (tcse71 + 3.33333333333333e-2F*(-u[t0][x][y][z - 3] + u[t0][x][y][z - 1]) + 4.16666666666667e-3F*u[t0][x][y][z - 4])*r1[x][y][z - 2])*r1[x][y][z - 2]);
                                  u[t1][x][y][z] = 1.0F*(tcse77*u[t2][x][y][z] + 4.01111156356834F*(tcse142*epsilon[x][y][z] + (-3.75e-2F*(tcse120*ti32 + tcse132*r3[x][y][z - 1] + tcse133*ti35) + 1.25e-2F*(tcse125*ti26 + tcse136*r3[x][y][z - 2] + tcse137*ti29) + 6.25e-3F*(tcse126*ti20 + tcse138*r3[x][y][z + 1] + tcse139*ti23) + 2.08333333333333e-2F*(tcse127*ti1 + tcse140*ti2 + tcse140*r3[x][y][z]) + 1.66666666666667e-2F*(-tcse115*r3[x][y][z - 1] - tcse116*ti35 + tcse121*r3[x][y][z + 1] + tcse122*ti23 - tcse128*ti32 + tcse134*ti20) + 2.08333333333333e-3F*(-tcse117*ti40 - tcse118*r3[x][y][z + 2] - tcse119*ti38 + tcse123*r3[x][y][z - 2] + tcse124*ti29 - tcse129*ti43 - tcse130*r3[x][y][z - 3] - tcse131*ti47 + tcse135*ti26))*delta[x][y][z]) + 4.0F*m[x][y][z]*u[t0][x][y][z])/tcse78;
                                  v[t1][x][y][z] = 1.0F*(4.01111156356834F*tcse142*delta[x][y][z] + tcse77*v[t2][x][y][z] - 1.50416683633813e-1F*(tcse120*ti32 + tcse132*r3[x][y][z - 1] + tcse133*ti35) + 5.01388945446042e-2F*(tcse125*ti26 + tcse136*r3[x][y][z - 2] + tcse137*ti29) + 2.50694472723021e-2F*(tcse126*ti20 + tcse138*r3[x][y][z + 1] + tcse139*ti23) + 8.35648242410071e-2F*(tcse127*ti1 + tcse140*ti2 + tcse140*r3[x][y][z]) + 6.68518593928056e-2F*(-tcse115*r3[x][y][z - 1] - tcse116*ti35 + tcse121*r3[x][y][z + 1] + tcse122*ti23 - tcse128*ti32 + tcse134*ti20) + 8.3564824241007e-3F*(-tcse117*ti40 - tcse118*r3[x][y][z + 2] - tcse119*ti38 + tcse123*r3[x][y][z - 2] + tcse124*ti29 - tcse129*ti43 - tcse130*r3[x][y][z - 3] - tcse131*ti47 + tcse135*ti26) + 4.0F*m[x][y][z]*v[t0][x][y][z])/tcse78;
                              }
                          }
                      }
                  }
              }
          }
      }
    }
    gettimeofday(&end_loop_x_1, NULL);
    timings->loop_x_1 += (double)(end_loop_x_1.tv_sec-start_loop_x_1.tv_sec)+(double)(end_loop_x_1.tv_usec-start_loop_x_1.tv_usec)/1000000;
    struct timeval start_loop_p_src_2, end_loop_p_src_2;
    gettimeofday(&start_loop_p_src_2, NULL);
    for (int p_src = 0; p_src < p_src_size; p_src += 1)
    {
      float tcse0 = (float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p_src][2])) + src_coords[p_src][2]);
      float tcse1 = (float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p_src][1])) + src_coords[p_src][1]);
      float tcse2 = (float)(-2.0e+1F*(int)(floor(5.0e-2F*src_coords[p_src][0])) + src_coords[p_src][0]);
      float tcse3 = -2.5e-3F*tcse1*tcse2;
      float tcse4 = -2.5e-3F*tcse0*tcse1;
      float tcse5 = -2.5e-3F*tcse0*tcse2;
      float tcse6 = 2.5e-3F*tcse1*tcse2;
      float tcse7 = 2.5e-3F*tcse0*tcse1;
      float tcse8 = 2.5e-3F*tcse0*tcse2;
      float tcse9 = 1.25e-4F*tcse0*tcse1*tcse2;
      float tcse10 = -1.25e-4F*tcse0*tcse1*tcse2;
      float tcse11 = 2.50694472723021e-4F*tcse0*tcse1*tcse2*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      float tcse12 = 2.00555578178417F*(tcse10 + tcse8)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      float tcse13 = 2.00555578178417F*(tcse10 + tcse7)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      float tcse14 = 2.00555578178417F*(tcse10 + tcse6)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      float tcse15 = 2.00555578178417F*(5.0e-2F*tcse0 + tcse4 + tcse5 + tcse9)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      float tcse16 = 2.00555578178417F*(5.0e-2F*tcse2 + tcse3 + tcse5 + tcse9)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      float tcse17 = 2.00555578178417F*(5.0e-2F*tcse1 + tcse3 + tcse4 + tcse9)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      float tcse18 = 2.00555578178417F*(-5.0e-2F*tcse0 - 5.0e-2F*tcse1 + tcse10 - 5.0e-2F*tcse2 + tcse6 + tcse7 + tcse8 + 1)*src[time][p_src]/m[(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = tcse18 + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = tcse17 + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = tcse16 + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = tcse15 + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = tcse14 + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = tcse13 + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = tcse12 + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = tcse11 + u[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = tcse18 + v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = tcse17 + v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = tcse16 + v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = tcse15 + v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10] = tcse14 + v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 10];
      v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = tcse13 + v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = tcse12 + v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 10][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
      v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11] = tcse11 + v[t1][(int)(floor(5.0e-2F*src_coords[p_src][0])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][1])) + 11][(int)(floor(5.0e-2F*src_coords[p_src][2])) + 11];
    }
    gettimeofday(&end_loop_p_src_2, NULL);
    timings->loop_p_src_2 += (double)(end_loop_p_src_2.tv_sec-start_loop_p_src_2.tv_sec)+(double)(end_loop_p_src_2.tv_usec-start_loop_p_src_2.tv_usec)/1000000;
    struct timeval start_loop_p_rec_3, end_loop_p_rec_3;
    gettimeofday(&start_loop_p_rec_3, NULL);
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
      rec[time][p_rec] = 1.25e-4F*tcse0*tcse1*tcse2*v[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 11] + (tcse10 + tcse6)*v[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 11] + (tcse10 + tcse7)*v[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 11] + (tcse10 + tcse8)*v[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 10] + (5.0e-2F*tcse0 + tcse4 + tcse5 + tcse9)*v[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 10] + (5.0e-2F*tcse1 + tcse3 + tcse4 + tcse9)*v[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 11] + (5.0e-2F*tcse2 + tcse3 + tcse5 + tcse9)*v[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 11][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 10] + (-5.0e-2F*tcse0 - 5.0e-2F*tcse1 + tcse10 - 5.0e-2F*tcse2 + tcse6 + tcse7 + tcse8 + 1)*v[t0][(int)(floor(5.0e-2F*rec_coords[p_rec][0])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][1])) + 10][(int)(floor(5.0e-2F*rec_coords[p_rec][2])) + 10];
    }
    gettimeofday(&end_loop_p_rec_3, NULL);
    timings->loop_p_rec_3 += (double)(end_loop_p_rec_3.tv_sec-start_loop_p_rec_3.tv_sec)+(double)(end_loop_p_rec_3.tv_usec-start_loop_p_rec_3.tv_usec)/1000000;
  }
  free(r0);
  free(r1);
  free(r2);
  free(r3);
  return 0;
}

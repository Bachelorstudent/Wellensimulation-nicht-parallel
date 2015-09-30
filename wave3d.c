
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
//#include <omp.h>

#include <emmintrin.h>
#include <immintrin.h>
#include <mm_malloc.h>
#include "wavestate3d.h"

/*
zum ausführen:
./wave3d nx ny nz iterations
*/
int
main(int argc, char **argv)
{
  float time;
  clock_t s, e;

  pwavestate3d wv;
  pwavestate3dAllign wvVec;
  real delta;
  int cnt;
  int nx, ny, nz, steps, iterations;

  nx = atoi(argv[argc-4]);
  ny = atoi(argv[argc-3]);
  nz = atoi(argv[argc-2]);
  iterations = atoi(argv[argc-1]);
  steps = 1;
  delta = 0.0001;

  wv = new_wavestate3d(nx, ny, nz, 1.0/(nx+1), 1.0, steps);
  wvVec = new_wavestate3d_allign(nx, ny, nz, 1.0/(nx+1), 1.0, steps);

  zero_wavestate3d(wv, steps);
  zero_wavestate3d_allign(wvVec, steps);
  
  boundary_wavestate3d(wv, 0.1, steps);
  boundary_wavestate3d_allign(wvVec, 0.1, steps);

  printf("start leapfrog\n");
  s = clock();
  for(cnt=0;cnt<iterations;++cnt){
    leapfrog_wavestate3d(wv, delta, steps);
  }
  e = clock();
  time = (double) (e - s) / (double) CLOCKS_PER_SEC;
  printf("time leapfrog:\t\t %.3f s\n", time);
  
  s = clock();
  for(cnt=0;cnt<iterations;++cnt){
    leapfrog_wavestate3d_v2(wv, delta, steps);
  }
  e = clock();
  time = (double) (e - s) / (double) CLOCKS_PER_SEC;
  printf("time leapfrog v2:\t %.3f s\n", time);
  
  s = clock();
  for(cnt=0;cnt<iterations;++cnt){
    leapfrog_wavestate3d_sse2(wvVec, delta, steps);
  }
  e = clock();
  time = (double) (e - s) / (double) CLOCKS_PER_SEC;
  printf("time leapfrog sse2:\t %.3f s\n", time);
/*
  s = clock();
  for(cnt=0;cnt<iterations;++cnt){
    leapfrog_wavestate3d_avx(wvVec, delta, steps);
  }
  e = clock();
  time = (double) (e - s) / (double) CLOCKS_PER_SEC;
  printf("time leapfrog avx:\t %.3f s\n", time);
*/

  del_wavestate3d(wv);
  del_wavestate3d_allign(wvVec);
  
  return 0;
}

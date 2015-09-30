#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include <omp.h>

#include <emmintrin.h>
#include <immintrin.h>
#include <mm_malloc.h>
#include "wavestate3d.h"

/* ------------------------------------------------------------
   Create a wavestate3d object
   ------------------------------------------------------------ */

pwavestate3d
new_wavestate3d(int nx, int ny, int nz, real h, real crho, int steps)
{
  pwavestate3d wv;
  int dimensionLengthX = nx + steps * 2;
  int dimensionLengthY = ny + steps * 2;
  int dimensionLengthZ = nz + steps * 2;  
  //counter
  int i,j;

  wv = (pwavestate3d) malloc(sizeof(wavestate3d));
  
  wv->x = (real ***) malloc((size_t) sizeof(real) * dimensionLengthZ);
  for(i = 0; i < dimensionLengthZ; i++){
    wv->x[i] = (real **) malloc((size_t) sizeof(real) * dimensionLengthY);
    for(j = 0; j < dimensionLengthY; j++){
      wv->x[i][j] = (real *) malloc((size_t) sizeof(real) * dimensionLengthX);
    }
  }

  wv->v = (real ***) malloc((size_t) sizeof(real) * dimensionLengthZ);
  for(i = 0; i < dimensionLengthZ; i++){
    wv->v[i] = (real **) malloc((size_t) sizeof(real) * dimensionLengthY);
    for(j = 0; j < dimensionLengthY; j++){
      wv->v[i][j] = (real *) malloc((size_t) sizeof(real) * dimensionLengthX);
    }
  }

  wv->nx = nx;
  wv->ny = ny;
  wv->nz = nz;
  wv->h = h;
  wv->crho = crho;
  
  return wv;
}

pwavestate3dAllign
new_wavestate3d_allign(int nx, int ny, int nz, real h, real crho, int steps)
{
  pwavestate3dAllign wvVec;
  int dimensionLength = (nz + steps * 2) * (ny + steps * 2) * (nx + steps * 2);

  wvVec = (pwavestate3dAllign) malloc(sizeof(wavestate3dAllign));

  wvVec->x = (real*) _mm_malloc(dimensionLength * sizeof(real), 32);
  wvVec->v = (real*) _mm_malloc(dimensionLength * sizeof(real), 32);

  wvVec->nx = nx;
  wvVec->ny = ny;
  wvVec->nz = nz;
  wvVec->h = h;
  wvVec->crho = crho;
  
  return wvVec;
}
/* ------------------------------------------------------------
   Delete a wavestate3d object
   ------------------------------------------------------------ */

void
del_wavestate3d(pwavestate3d wv)
{
  free(wv->v);
  free(wv->x);

  wv->v = 0;			/* Safety measure */
  wv->x = 0;

  free(wv);
}

void
del_wavestate3d_allign(pwavestate3dAllign wv)
{
  free(wv->v);
  free(wv->x);

  wv->v = 0;			/* Safety measure */
  wv->x = 0;

  free(wv);
}
/* ------------------------------------------------------------
   Set displacements and velocities to zero
   ------------------------------------------------------------ */

void
zero_wavestate3d(pwavestate3d wv, int steps)
{
  real ***x = wv->x;
  real ***v = wv->v;
  int nx = wv->nx + steps*2;
  int ny = wv->ny + steps*2;
  int nz = wv->nz + steps*2;
  int i,j,k;

  for(i=0; i<nz; i++){
    for(j=0; j<ny; j++){
      for(k=0; k<nx; k++){
	x[i][j][k] = 0.0;
      }
    }
  }

  for(i=0; i<nz; i++){
    for(j=0; j<ny; j++){
      for(k=0; k<nx; k++){
	v[i][j][k] = 0.0;
      }
    }
  }
  
}
void
zero_wavestate3d_allign(pwavestate3dAllign wv, int steps)
{
  real *x = wv->x;
  real *v = wv->v;
  int nx = wv->nx + steps*2;
  int ny = wv->ny + steps*2;
  int nz = wv->nz + steps*2;
  int i,j,k;

  for(i=0; i<nz; i++){
    for(j=0; j<ny; j++){
      for(k=0; k<nx; k++){
	x[i*ny*nx + j*nx + k] = 0.0;
      }
    }
  }

  for(i=0; i<nz; i++){
    for(j=0; j<ny; j++){
      for(k=0; k<nx; k++){
	v[i*ny*nx + j*nx + k] = 0.0;
      }
    }
  }
  
}

/* ------------------------------------------------------------
   Set interesting boundary values
   ------------------------------------------------------------ */

void
boundary_wavestate3d(pwavestate3d wv, real t, int steps)
{
  real ***x = wv->x;

  if(0.0 < t && t < 0.25)
  {
    x[0][0][steps] = sin(M_PI * t / 0.125);
    x[0][steps][0] = sin(M_PI * t / 0.125);
    x[steps][0][0] = sin(M_PI * t / 0.125);
  }
}

void
boundary_wavestate3d_allign(pwavestate3dAllign wv, real t, int steps)
{
  real *x = wv->x;
  
  int arrayCntX = wvVec->nx + steps * 2;
  int arrayCntXY = wvVec->ny + steps * 2;
  arrayCntXY *= arrayCntX;
  
  if(0.0 < t && t < 0.25)
    x[arrayCntXY + arrayCntX + 1] = sin(M_PI * t / 0.125);
}

/* ------------------------------------------------------------
   Perform one step of the leapfrog method
   ------------------------------------------------------------ */

void
leapfrog_wavestate3d(pwavestate3d wv, real delta, int steps)
{
  real ***x = wv->x;
  real ***v = wv->v;
  int lastImportantField[3];
  lastImportantField[0] = wv->nz + steps*2 - 1;
  lastImportantField[1] = wv->ny + steps*2 - 1;
  lastImportantField[2] = wv->nx + steps*2 - 1;
  real h = wv->h;
  real crho = wv->crho;
  int stepsLocal, i,j,k;

  for(stepsLocal=0; stepsLocal<steps && lastImportantField[0]>=wv->nz+steps && lastImportantField[1]>=wv->ny+steps && lastImportantField[2]>=wv->nx+steps; stepsLocal++, lastImportantField[0]--, lastImportantField[1]--, lastImportantField[2]--){
    for(i=stepsLocal+1; i<lastImportantField[0]; i++){
      for(j=stepsLocal+1; j<lastImportantField[1]; j++){
        for(k=stepsLocal+1; k<lastImportantField[2]; k++){
          //Update velocities
	  v[i][j][k] += delta * crho * (x[i+1][j][k] + x[i-1][j][k] + x[i][j+1][k] + x[i][j-1][k] + x[i][j][k+1] + x[i][j][k-1] - 6.0 * x[i][j][k]) / h / h;
          //Update displacements	
          x[i][j][k] += delta * v[i][j][k];
        }
      }
    }
  }
}

void
leapfrog_wavestate3d_v2(pwavestate3d wv, real delta, int steps)
{
  real ***x = wv->x;
  real ***v = wv->v;
  int lastImportantField[3];
  lastImportantField[0] = wv->nz + steps*2 - 1;
  lastImportantField[1] = wv->ny + steps*2 - 1;
  lastImportantField[2] = wv->nx + steps*2 - 1;
  real h = wv->h;
  real crho = wv->crho;
  int stepsLocal, i,j,k;
  real constVal = (crho*delta)/(h*h);
  
  for(stepsLocal=0; stepsLocal<steps && lastImportantField[0]>=wv->nz+steps && lastImportantField[1]>=wv->ny+steps && lastImportantField[2]>=wv->nx+steps; stepsLocal++, lastImportantField[0]--, lastImportantField[1]--, lastImportantField[2]--){
    for(i=stepsLocal+1; i<lastImportantField[0]; i++){
      for(j=stepsLocal+1; j<lastImportantField[1]; j++){
        for(k=stepsLocal+1; k<lastImportantField[2]; k++){
          //Update velocities
	  v[i][j][k] += (x[i][j][k+1] + x[i][j][k-1] + x[i][j+1][k] + x[i][j-1][k] + x[i+1][j][k] + x[i-1][j][k] - 6.0 * x[i][j][k]) * constVal;
          //Update displacements
	  x[i][j][k] += delta * v[i][j][k];
        }
      }
    }
  }
}

void
leapfrog_wavestate3d_sse2(pwavestate3dAllign wvVec, real delta, int steps)
{
  real *x = wvVec->x;
  real *v = wvVec->v;
  int lastImportantField[3];
  lastImportantField[0] = wvVec->nz + steps*2 - 1;
  lastImportantField[1] = wvVec->ny + steps*2 - 1;
  lastImportantField[2] = wvVec->nx + steps*2 - 1;
  real h = wvVec->h;
  real crho = wvVec->crho;
  int stepsLocal, i,j,k;
  const real constVal = (crho*delta)/(h*h);
  const real const6 = 6.0;

  int arrayCntX = wvVec->nx + steps * 2;
  int arrayCntXY = wvVec->ny + steps * 2;
  arrayCntXY *= arrayCntX;

#ifdef __SSE2__
__m128d rCenter, rLeft, rRight, rConstVal1, rConstVal2, rResult, rConst6;
#endif

  rConstVal1 = _mm_load_pd1(&constVal);
  rConstVal2 = _mm_load_pd1(&delta);
  rConst6 = _mm_load_pd1(&const6);

  for(stepsLocal=0; stepsLocal<steps && lastImportantField[0]>=wvVec->nz+steps && lastImportantField[1]>=wvVec->ny+steps && lastImportantField[2]>=wvVec->nx+steps; stepsLocal++, lastImportantField[0]--, lastImportantField[1]--, lastImportantField[2]--){
    for(i=stepsLocal+1; i<lastImportantField[0]; i++){
      for(j=stepsLocal+1; j<lastImportantField[1]; j++){

        k=stepsLocal+1;
        if((i*arrayCntXY + j*arrayCntX + k) % 2){
          v[i*arrayCntXY + j*arrayCntX + k] += (x[i*arrayCntXY + j*arrayCntX + k + 1] + x[i*arrayCntXY + j*arrayCntX + k - 1] + x[i*arrayCntXY + (j+1)*arrayCntX + k] + x[i*arrayCntXY + (j-1)*arrayCntX + k] + x[(i+1)*arrayCntXY + j*arrayCntX + k] + x[(i-1)*arrayCntXY + j*arrayCntX + k] - 6.0 * x[i*arrayCntXY + j*arrayCntX + k]) * constVal;
          x[i*arrayCntXY + j*arrayCntX + k] += delta * v[i*arrayCntXY + j*arrayCntX + k];
          k++;
        }

        for(; k+1<lastImportantField[2]; k+=2){
          rCenter = _mm_load_pd(&x[i*arrayCntXY + j*arrayCntX + k]);
          
          rLeft = _mm_load_pd(&x[i*arrayCntXY + j*arrayCntX + k - 2]);
          rRight = _mm_load_pd(&x[i*arrayCntXY + j*arrayCntX + k + 2]);
          
          rLeft = _mm_shuffle_pd(rLeft, rCenter, 0x00010000);
          rRight = _mm_shuffle_pd(rCenter, rRight, 0x00010000);
          rCenter = _mm_mul_pd(rCenter, rConst6);

          rResult = _mm_add_pd(rRight, rLeft);//printf("werte: %d %d %d %d   %d\n",i,j,k,arrayCnt,i*arrayCnt*arrayCnt + (j-1)*arrayCnt + k);
          rResult = _mm_add_pd(rResult, _mm_load_pd(&x[i*arrayCntXY + (j-1)*arrayCntX + k]));
          rResult = _mm_add_pd(rResult, _mm_load_pd(&x[i*arrayCntXY + (j+1)*arrayCntX + k]));
          rResult = _mm_add_pd(rResult, _mm_load_pd(&x[(i+1)*arrayCntXY + j*arrayCntX + k]));
          rResult = _mm_add_pd(rResult, _mm_load_pd(&x[(i-1)*arrayCntXY + j*arrayCntX + k]));

          rResult = _mm_sub_pd(rResult, rCenter);
          rResult = _mm_mul_pd(rResult, rConstVal1);
          rResult = _mm_add_pd(rResult, _mm_load_pd(&v[i*arrayCntXY + j*arrayCntX + k]));          
          _mm_store_pd(&v[i*arrayCntXY + j*arrayCntX + k], rResult);


          rResult = _mm_load_pd(&v[i*arrayCntXY + j*arrayCntX + k]);
          rResult = _mm_mul_pd(rResult, rConstVal2);
          rCenter = _mm_load_pd(&x[i*arrayCntXY + j*arrayCntX + k]);
          rResult = _mm_add_pd(rResult, rCenter);
          _mm_store_pd(&x[i*arrayCntXY + j*arrayCntX + k], rResult);
        }
 
        if(k+1 == lastImportantField[2]){
          v[i*arrayCntXY + j*arrayCntX + k] += (x[i*arrayCntXY + j*arrayCntX + k + 1] + x[i*arrayCntXY + j*arrayCntX + k - 1] + x[i*arrayCntXY + (j+1)*arrayCntX + k] + x[i*arrayCntXY + (j-1)*arrayCntX + k] + x[(i+1)*arrayCntXY + j*arrayCntX + k] + x[(i-1)*arrayCntXY + j*arrayCntX + k] - 6.0 * x[i*arrayCntXY + j*arrayCntX + k]) * constVal;
          x[i*arrayCntXY + j*arrayCntX + k] += delta * v[i*arrayCntXY + j*arrayCntX + k];
        }
      }
    }
  }
}

/*
void
leapfrog_wavestate3d_avx(pwavestate3dAllign wvVec, real delta, int steps)
{
  real *x = wvVec->x;
  real *v = wvVec->v;
  int lastImportantField[3];
  lastImportantField[0] = wvVec->nz + steps*2 - 1;
  lastImportantField[1] = wvVec->ny + steps*2 - 1;
  lastImportantField[2] = wvVec->nx + steps*2 - 1;
  real h = wvVec->h;
  real crho = wvVec->crho;
  int stepsLocal, i,j,k;
  const real constVal = (crho*delta)/(h*h);
  const real const6 = 6.0;

  int arrayCntX = wvVec->nx + steps * 2;
  int arrayCntXY = wvVec->ny + steps * 2;
  arrayCntXY *= arrayCntX;

#ifdef __AVX__
  __m256d rCenter, rLeft, rRight, rConstVal1, rConstVal2, rResult, rConst6;
#endif

  rConstVal1 = _mm256_broadcast_sd(&constVal);
  rConstVal2 = _mm256_broadcast_sd(&delta);
  rConst6 = _mm256_broadcast_sd(&const6);

  for(stepsLocal=0; stepsLocal<steps && lastImportantField[0]>=wvVec->nz+steps && lastImportantField[1]>=wvVec->ny+steps && lastImportantField[2]>=wvVec->nx+steps; stepsLocal++, lastImportantField[0]--, lastImportantField[1]--, lastImportantField[2]--){
    for(i=stepsLocal+1; i<lastImportantField[0]; i++){
      for(j=stepsLocal+1; j<lastImportantField[1]; j++){

        k=stepsLocal+1;
        while((i*arrayCntXY + j*arrayCntX + k) % 4){
          v[i*arrayCntXY + j*arrayCntX + k] += (x[i*arrayCntXY + j*arrayCntX + k + 1] + x[i*arrayCntXY + j*arrayCntX + k - 1] + x[i*arrayCntXY + (j+1)*arrayCntX + k] + x[i*arrayCntXY + (j-1)*arrayCntX + k] + x[(i+1)*arrayCntXY + j*arrayCntX + k] + x[(i-1)*arrayCntXY + j*arrayCntX + k] - 6.0 * x[i*arrayCntXY + j*arrayCntX + k]) * constVal;
          x[i*arrayCntXY + j*arrayCntX + k] += delta * v[i*arrayCntXY + j*arrayCntX + k];
          k++;
        }
		
        //avx zusatz
        if(k+3>=lastImportantField[2]){
          for(; k<lastImportantField[2]; ++k){
            v[i*arrayCntXY + j*arrayCntX + k] += (x[i*arrayCntXY + j*arrayCntX + k + 1] + x[i*arrayCntXY + j*arrayCntX + k - 1] + x[i*arrayCntXY + (j+1)*arrayCntX + k] + x[i*arrayCntXY + (j-1)*arrayCntX + k] + x[(i+1)*arrayCntXY + j*arrayCntX + k] + x[(i-1)*arrayCntXY + j*arrayCntX + k] - 6.0 * x[i*arrayCntXY + j*arrayCntX + k]) * constVal;
            x[i*arrayCntXY + j*arrayCntX + k] += delta * v[i*arrayCntXY + j*arrayCntX + k];
          }
        }

        for(; k+3<lastImportantField[2]; k+=4){
          rCenter = _mm256_load_pd(&x[i*arrayCntXY + j*arrayCntX + k]);
          
          rLeft = _mm256_load_pd(&x[i*arrayCntXY + j*arrayCntX + k - 1]);
          rRight = _mm256_load_pd(&x[i*arrayCntXY + j*arrayCntX + k + 1]);
          
          rCenter = _mm256_mul_pd(rCenter, rConst6);

          rResult = _mm256_add_pd(rRight, rLeft);//printf("werte: %d %d %d %d   %d\n",i,j,k,arrayCnt,i*arrayCnt*arrayCnt + (j-1)*arrayCnt + k);
          rResult = _mm256_add_pd(rResult, _mm256_load_pd(&x[i*arrayCntXY + (j-1)*arrayCntX + k]));
          rResult = _mm256_add_pd(rResult, _mm256_load_pd(&x[i*arrayCntXY + (j+1)*arrayCntX + k]));
          rResult = _mm256_add_pd(rResult, _mm256_load_pd(&x[(i+1)*arrayCntXY + j*arrayCntX + k]));
          rResult = _mm256_add_pd(rResult, _mm256_load_pd(&x[(i-1)*arrayCntXY + j*arrayCntX + k]));

          rResult = _mm256_sub_pd(rResult, rCenter);
          rResult = _mm256_mul_pd(rResult, rConstVal1);
          rResult = _mm256_add_pd(rResult, _mm256_load_pd(&v[i*arrayCntXY + j*arrayCntX + k]));          
          _mm256_store_pd(&v[i*arrayCntXY + j*arrayCntX + k], rResult);


          rResult = _mm256_load_pd(&v[i*arrayCntXY + j*arrayCntX + k]);
          rResult = _mm256_mul_pd(rResult, rConstVal2);
          rCenter = _mm256_load_pd(&x[i*arrayCntXY + j*arrayCntX + k]);
          rResult = _mm256_add_pd(rResult, rCenter);
          _mm256_store_pd(&x[i*arrayCntXY + j*arrayCntX + k], rResult);
        }
		
       for(; k<lastImportantField[2]; ++k){
          v[i*arrayCntXY + j*arrayCntX + k] += (x[i*arrayCntXY + j*arrayCntX + k + 1] + x[i*arrayCntXY + j*arrayCntX + k - 1] + x[i*arrayCntXY + (j+1)*arrayCntX + k] + x[i*arrayCntXY + (j-1)*arrayCntX + k] + x[(i+1)*arrayCntXY + j*arrayCntX + k] + x[(i-1)*arrayCntXY + j*arrayCntX + k] - 6.0 * x[i*arrayCntXY + j*arrayCntX + k]) * constVal;
          x[i*arrayCntXY + j*arrayCntX + k] += delta * v[i*arrayCntXY + j*arrayCntX + k];
        }		
      }
    }
  }
}
*/

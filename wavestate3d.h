
#ifndef WAVESTATE3D_H
#define WAVESTATE3D_H

typedef struct _wavestate3d wavestate3d;
typedef wavestate3d *pwavestate3d;
typedef const wavestate3d *pcwavestate3d;

typedef struct _wavestate3dAllign wavestate3dAllign;
typedef wavestate3dAllign *pwavestate3dAllign;
typedef const wavestate3dAllign *pcwavestate3dAllign;

#include "settings.h"

struct _wavestate3d {
  int nx;			/* Number of point masses in x direction */
  int ny;			/* Number of point masses in y direction */
  int nz;			/* Number of point masses in z direction */

  real h;			/* Step width */
  real crho;			/* Elasticity / density constant */

  real ***x;			/* Displacements */
  real ***v;			/* Velocities */
};

struct _wavestate3dAllign {
  int nx;			/* Number of point masses in x direction */
  int ny;			/* Number of point masses in y direction */
  int nz;			/* Number of point masses in z direction */

  real h;			/* Step width */
  real crho;			/* Elasticity / density constant */

  real *x;			/* Displacements */
  real *v;			/* Velocities */
};

/* Create a wavestate3d object */
pwavestate3d
new_wavestate3d(int nx, int ny, int nz, real h, real crho, int steps);

pwavestate3dAllign
new_wavestate3d_allign(int nx, int ny, int nz, real h, real crho, int steps);

/* Delete a wavestate3d object */
void
del_wavestate3d(pwavestate3d wv);

void
del_wavestate3d_allign(pwavestate3dAllign wv);


/* Set displacements and velocities to zero */
void
zero_wavestate3d(pwavestate3d wv, int steps);

void
zero_wavestate3d_allign(pwavestate3dAllign wv, int steps);

/* Set interesting boundary values */
void
boundary_wavestate3d(pwavestate3d wv, real t);

void
boundary_wavestate3d_allign(pwavestate3dAllign wv, real t);

/* Perform one step of the leapfrog method */
void
leapfrog_wavestate3d(pwavestate3d wv, real delta, int steps);
void
leapfrog_wavestate3d_v2(pwavestate3d wv, real delta, int steps);
void
leapfrog_wavestate3d_sse2(pwavestate3dAllign wvVec, real delta, int steps);
//void
//leapfrog_wavestate3d_avx(pwavestate3dAllign wvVec, real delta, int steps);
#endif

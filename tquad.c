#include <stdio.h>
#include <math.h>

#include <glib.h>

#include "intsincos.h"

#define M411  0
#define M412  1
#define M413  2
#define M414  3
#define M421  4
#define M422  5
#define M423  6
#define M424  7
#define M431  8
#define M432  9
#define M433 10
#define M434 11
#define M441 12
#define M442 13
#define M443 14
#define M444 15

static void invert4x4(gdouble *Ai, gdouble *A)

{
  gdouble det ;

  det = 
    A[M411]*A[M422]*A[M433]*A[M444] + A[M411]*A[M423]*A[M434]*A[M442] + 
    A[M411]*A[M424]*A[M432]*A[M443] +
    A[M412]*A[M421]*A[M434]*A[M443] + A[M412]*A[M423]*A[M431]*A[M444] + 
    A[M412]*A[M424]*A[M433]*A[M441] +
    A[M413]*A[M421]*A[M432]*A[M444] + A[M413]*A[M422]*A[M434]*A[M441] + 
    A[M413]*A[M424]*A[M431]*A[M442] +
    A[M414]*A[M421]*A[M433]*A[M442] + A[M414]*A[M422]*A[M431]*A[M443] + 
    A[M414]*A[M423]*A[M432]*A[M441] -
    A[M411]*A[M422]*A[M434]*A[M443] - A[M411]*A[M423]*A[M432]*A[M444] -
    A[M411]*A[M424]*A[M433]*A[M442] -
    A[M412]*A[M421]*A[M433]*A[M444] - A[M412]*A[M423]*A[M434]*A[M441] -
    A[M412]*A[M424]*A[M431]*A[M443] -
    A[M413]*A[M421]*A[M434]*A[M442] - A[M413]*A[M422]*A[M431]*A[M444] -
    A[M413]*A[M424]*A[M432]*A[M441] -
    A[M414]*A[M421]*A[M432]*A[M443] - A[M414]*A[M422]*A[M433]*A[M441] -
    A[M414]*A[M423]*A[M431]*A[M442] ;

  det = 1.0/det ;

  Ai[M411] = 
    A[M422]*A[M433]*A[M444] + A[M423]*A[M434]*A[M442] + 
    A[M424]*A[M432]*A[M443] -
    A[M422]*A[M434]*A[M443] - A[M423]*A[M432]*A[M444] - 
    A[M424]*A[M433]*A[M442] ;
  Ai[M412] = 
    A[M412]*A[M434]*A[M443] + A[M413]*A[M432]*A[M444] + 
    A[M414]*A[M433]*A[M442] -
    A[M412]*A[M433]*A[M444] - A[M413]*A[M434]*A[M442] - 
    A[M414]*A[M432]*A[M443] ;
  Ai[M413] = 
    A[M412]*A[M423]*A[M444] + A[M413]*A[M424]*A[M442] + 
    A[M414]*A[M422]*A[M443] -
    A[M412]*A[M424]*A[M443] - A[M413]*A[M422]*A[M444] - 
    A[M414]*A[M423]*A[M442] ;
  Ai[M414] = 
    A[M412]*A[M424]*A[M433] + A[M413]*A[M422]*A[M434] + 
    A[M414]*A[M423]*A[M432] -
    A[M412]*A[M423]*A[M434] - A[M413]*A[M424]*A[M432] - 
    A[M414]*A[M422]*A[M433] ;
  Ai[M421] = 
    A[M421]*A[M434]*A[M443] + A[M423]*A[M431]*A[M444] + 
    A[M424]*A[M433]*A[M441] -
    A[M421]*A[M433]*A[M444] - A[M423]*A[M434]*A[M441] - 
    A[M424]*A[M431]*A[M443] ;
  Ai[M422] = 
    A[M411]*A[M433]*A[M444] + A[M413]*A[M434]*A[M441] + 
    A[M414]*A[M431]*A[M443] -
    A[M411]*A[M434]*A[M443] - A[M413]*A[M431]*A[M444] - 
    A[M414]*A[M433]*A[M441] ;
  Ai[M423] = 
    A[M411]*A[M424]*A[M443] + A[M413]*A[M421]*A[M444] +
    A[M414]*A[M423]*A[M441] -
    A[M411]*A[M423]*A[M444] - A[M413]*A[M424]*A[M441] - 
    A[M414]*A[M421]*A[M443] ;
  Ai[M424] = 
    A[M411]*A[M423]*A[M434] + A[M413]*A[M424]*A[M431] + 
    A[M414]*A[M421]*A[M433] -
    A[M411]*A[M424]*A[M433] - A[M413]*A[M421]*A[M434] - 
    A[M414]*A[M423]*A[M431] ;
  Ai[M431] = 
    A[M421]*A[M432]*A[M444] + A[M422]*A[M434]*A[M441] + 
    A[M424]*A[M431]*A[M442] -
    A[M421]*A[M434]*A[M442] - A[M422]*A[M431]*A[M444] - 
    A[M424]*A[M432]*A[M441] ;
  Ai[M432] = 
    A[M411]*A[M434]*A[M442] + A[M412]*A[M431]*A[M444] +
    A[M414]*A[M432]*A[M441] -
    A[M411]*A[M432]*A[M444] - A[M412]*A[M434]*A[M441] - 
    A[M414]*A[M431]*A[M442] ;
  Ai[M433] = 
    A[M411]*A[M422]*A[M444] + A[M412]*A[M424]*A[M441] +
    A[M414]*A[M421]*A[M442] -
    A[M411]*A[M424]*A[M442] - A[M412]*A[M421]*A[M444] - 
    A[M414]*A[M422]*A[M441] ;
  Ai[M434] = 
    A[M411]*A[M424]*A[M432] + A[M412]*A[M421]*A[M434] +
    A[M414]*A[M422]*A[M431] -
    A[M411]*A[M422]*A[M434] - A[M412]*A[M424]*A[M431] - 
    A[M414]*A[M421]*A[M432] ;
  Ai[M441] = 
    A[M421]*A[M433]*A[M442] + A[M422]*A[M431]*A[M443] + 
    A[M423]*A[M432]*A[M441] -
    A[M421]*A[M432]*A[M443] - A[M422]*A[M433]*A[M441] - 
    A[M423]*A[M431]*A[M442] ;
  Ai[M442] = 
    A[M411]*A[M432]*A[M443] + A[M412]*A[M433]*A[M441] + 
    A[M413]*A[M431]*A[M442] -
    A[M411]*A[M433]*A[M442] - A[M412]*A[M431]*A[M443] - 
    A[M413]*A[M432]*A[M441] ;
  Ai[M443] = 
    A[M411]*A[M423]*A[M442] + A[M412]*A[M421]*A[M443] + 
    A[M413]*A[M422]*A[M441] -
    A[M411]*A[M422]*A[M443] - A[M412]*A[M423]*A[M441] - 
    A[M413]*A[M421]*A[M442] ;
  Ai[M444] = 
    A[M411]*A[M422]*A[M433] + A[M412]*A[M423]*A[M431] + 
    A[M413]*A[M421]*A[M432] -
    A[M411]*A[M423]*A[M432] - A[M412]*A[M421]*A[M433] - 
    A[M413]*A[M422]*A[M431] ;

  Ai[M411] *= det ; Ai[M412] *= det ; Ai[M413] *= det ; Ai[M414] *= det ;
  Ai[M421] *= det ; Ai[M422] *= det ; Ai[M423] *= det ; Ai[M424] *= det ;
  Ai[M431] *= det ; Ai[M432] *= det ; Ai[M433] *= det ; Ai[M434] *= det ;
  Ai[M441] *= det ; Ai[M442] *= det ; Ai[M443] *= det ; Ai[M444] *= det ;

  return ;
}

gdouble tet_quad_analytical(gdouble *xtet, gdouble *y,
			    gint gm,
			    gint mx, gint my, gint mz, 
			    int_sin_cos_workspace_t *w0,
			    int_sin_cos_workspace_t *w1,
			    gdouble *work, gdouble *I)

{
  gdouble J, x[12], z, co[3], cs[3], ct[3], cn[3], dJ, sc = 1.0 ;
  gint idx[] = {0, 1, 2}, N, nx = 3, inc[3], i, grad = 0 ;
  gboolean scale = TRUE ;

  N = MAX(mx + my + mz + 1, 4) ;

  /* N = 8 ; */

  J = 0 ;
  tet_projection(y, &(xtet[3*1]), &(xtet[3*2]), &(xtet[3*3]),
		 co, cs, ct, cn, 
		 &(x[3*0]), &(x[3*1]), &(x[3*2]), &z) ;
  if ( fabs(z) > 1e-12 ) {
    pyramid_quad_mn(x, idx, 3, nx, z, gm, N, I, work, w0, w1) ;
    dJ = pyramid_quad_general(gm, N, I, co, cs, ct, cn, mx, my, mz, z, grad) ;
    J += dJ ;
    /* fprintf(stderr, "%lg ", dJ) ; */
  }

  tet_projection(y, &(xtet[3*0]), &(xtet[3*3]), &(xtet[3*2]),
		 co, cs, ct, cn, 
		 &(x[3*0]), &(x[3*1]), &(x[3*2]), &z) ;
  if ( fabs(z) > 1e-12 ) {
    pyramid_quad_mn(x, idx, 3, nx, z, gm, N, I, work, w0, w1) ;
    dJ = pyramid_quad_general(gm, N, I, co, cs, ct, cn, mx, my, mz, z, grad) ;
    J += dJ ;
    /* fprintf(stderr, "%lg ", dJ) ; */
  }

  tet_projection(y, &(xtet[3*3]), &(xtet[3*0]), &(xtet[3*1]),
		 co, cs, ct, cn, 
		 &(x[3*0]), &(x[3*1]), &(x[3*2]), &z) ;
  if ( fabs(z) > 1e-12 ) {
    pyramid_quad_mn(x, idx, 3, nx, z, gm, N, I, work, w0, w1) ;
    dJ = pyramid_quad_general(gm, N, I, co, cs, ct, cn, mx, my, mz, z, grad) ;
    J += dJ ;
    /* fprintf(stderr, "%lg ", dJ) ; */
  }

  tet_projection(y, &(xtet[3*2]), &(xtet[3*1]), &(xtet[3*0]),
		 co, cs, ct, cn, 
		 &(x[3*0]), &(x[3*1]), &(x[3*2]), &z) ;
  if ( fabs(z) > 1e-12 ) {
    pyramid_quad_mn(x, idx, 3, nx, z, gm, N, I, work, w0, w1) ;
    dJ = pyramid_quad_general(gm, N, I, co, cs, ct, cn, mx, my, mz, z, grad) ;
    J += dJ ;
    /* fprintf(stderr, "%lg ", dJ) ; */
  }
  /* fprintf(stderr, "\n") ; */

  return J ;
}

gdouble tri_quad_analytical(gdouble *xtri, gdouble *y,
			    gint gm,
			    gint mx, gint my, gint mz, 
			    int_sin_cos_workspace_t *w0,
			    int_sin_cos_workspace_t *w1,
			    gdouble *work, gdouble *I)

{
  gdouble J, x[12], z, co[3], cs[3], ct[3], cn[3] ;
  gint idx[] = {0, 1, 2}, N, nx = 3 ;

  N = MAX(mx + my + mz, 4) ;

  J = 0 ;
  tet_projection(y, &(xtri[3*0]), &(xtri[3*1]), &(xtri[3*2]),
		 co, cs, ct, cn, 
		 &(x[3*0]), &(x[3*1]), &(x[3*2]), &z) ;
  poly_quad_mn(x, idx, 3, nx, z, gm, N, I, work, w0, w1) ;
  J += polygon_quad_general(gm, N, I, co, cs, ct, cn, mx, my, mz) ;

  return J ;
}

gint tet_quad_linear_analytical(gdouble *x0, gdouble *x1, 
				gdouble *x2, gdouble *x3,
				gdouble *f0, gdouble *f1, 
				gdouble *f2, gdouble *f3,
				gint gm,
				gint nc, gdouble *y,
				int_sin_cos_workspace_t *w0,
				int_sin_cos_workspace_t *w1,
				gdouble *work, 
				gdouble *I)

/*
  integral \int f R^gm dV for tet x 0, 1, 2, 3 with f interpolated
  linearly between vertices, nc number of components (e.g. three for
  vector problems), y field point, on exit I contains integral for
  each component of f 
*/

{
  gdouble J[4], A[16], Ai[16], co[3], cs[3], ct[3], cn[3], z, x[12], f[4] ;
  gdouble Ic[2048] ;
  gint nx = 3, idx[] = {0, 1, 2}, N = 4, grad = 0, i ;

  A[4*0+0] = 1.0 ; A[4*0+1] = x0[0] ; A[4*0+2] = x0[1] ; A[4*0+3] = x0[2] ;
  A[4*1+0] = 1.0 ; A[4*1+1] = x1[0] ; A[4*1+2] = x1[1] ; A[4*1+3] = x1[2] ;
  A[4*2+0] = 1.0 ; A[4*2+1] = x2[0] ; A[4*2+2] = x2[1] ; A[4*2+3] = x2[2] ;
  A[4*3+0] = 1.0 ; A[4*3+1] = x3[0] ; A[4*3+2] = x3[1] ; A[4*3+3] = x3[2] ;

  invert4x4(Ai, A) ;

  J[0] = J[1] = J[2] = J[3] = 0.0 ;

  tet_projection(y, x1, x2, x3, co, cs, ct, cn, 
		 &(x[3*0]), &(x[3*1]), &(x[3*2]), &z) ;
  if ( fabs(z) > 1e-12 ) {
    pyramid_quad_mn(x, idx, 3, nx, z, gm, N, Ic, work, w0, w1) ;
    J[0] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 0, 0, z, grad) ;
    J[1] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 1, 0, 0, z, grad) ;
    J[2] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 1, 0, z, grad) ;
    J[3] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 0, 1, z, grad) ;
  }

  tet_projection(y, x0, x3, x2, co, cs, ct, cn, 
		 &(x[3*0]), &(x[3*1]), &(x[3*2]), &z) ;
  if ( fabs(z) > 1e-12 ) {
    pyramid_quad_mn(x, idx, 3, nx, z, gm, N, Ic, work, w0, w1) ;
    J[0] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 0, 0, z, grad) ;
    J[1] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 1, 0, 0, z, grad) ;
    J[2] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 1, 0, z, grad) ;
    J[3] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 0, 1, z, grad) ;
  }

  tet_projection(y, x3, x0, x1, co, cs, ct, cn, 
		 &(x[3*0]), &(x[3*1]), &(x[3*2]), &z) ;
  if ( fabs(z) > 1e-12 ) {
    pyramid_quad_mn(x, idx, 3, nx, z, gm, N, Ic, work, w0, w1) ;
    J[0] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 0, 0, z, grad) ;
    J[1] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 1, 0, 0, z, grad) ;
    J[2] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 1, 0, z, grad) ;
    J[3] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 0, 1, z, grad) ;
  }

  tet_projection(y, x2, x1, x0, co, cs, ct, cn, 
		 &(x[3*0]), &(x[3*1]), &(x[3*2]), &z) ;
  if ( fabs(z) > 1e-12 ) {
    pyramid_quad_mn(x, idx, 3, nx, z, gm, N, Ic, work, w0, w1) ;
    J[0] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 0, 0, z, grad) ;
    J[1] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 1, 0, 0, z, grad) ;
    J[2] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 1, 0, z, grad) ;
    J[3] += pyramid_quad_general(gm, N, Ic, co, cs, ct, cn, 0, 0, 1, z, grad) ;
  }

  for ( i = 0 ; i < nc ; i ++ ) {
    f[0] = Ai[4*0+0]*f0[i]+Ai[4*0+1]*f1[i]+Ai[4*0+2]*f2[i]+Ai[4*0+3]*f3[i] ;
    f[1] = Ai[4*1+0]*f0[i]+Ai[4*1+1]*f1[i]+Ai[4*1+2]*f2[i]+Ai[4*1+3]*f3[i] ;
    f[2] = Ai[4*2+0]*f0[i]+Ai[4*2+1]*f1[i]+Ai[4*2+2]*f2[i]+Ai[4*2+3]*f3[i] ;
    f[3] = Ai[4*3+0]*f0[i]+Ai[4*3+1]*f1[i]+Ai[4*3+2]*f2[i]+Ai[4*3+3]*f3[i] ;

    I[i] = f[0]*J[0] + f[1]*J[1] + f[2]*J[2] + f[3]*J[3] ;
  }

  return 0 ;
}

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <glib.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_gamma.h>

#include "intsincos.h"

/*not needed here because it is picked up from GTS via GTV*/
/* #include "predicates.h" */

/*binomials up to 20*/
static gint _BINOMIALS[] = {
  1,
  1, 1,
  1, 2,  1,
  1, 3,  3,  1,
  1, 4,  6,  4,   1,
  1, 5, 10, 10,   5,   1,
  1, 6, 15, 20,  15,   6,   1,
  1, 7, 21, 35,  35,  21,   7,    1,
  1, 8, 28, 56,  70,  56,  28,    8,    1,
  1, 9, 36, 84, 126, 126,  84,   36,    9,    1,
  1,10, 45,120, 210, 252, 210,  120,   45,   10,   1,
  1,11, 55,165, 330, 462, 462,  330,  165,   55,  11,   1,
  1,12, 66,220, 495, 792, 924,  792,  495,  220,  66,  12,   1,
  1,13, 78,286, 715,1287,1716, 1716, 1287,  715, 286,  78,  13,  1,
  1,14, 91,364,1001,2002,3003, 3432, 3003, 2002,1001, 364,  91, 14,  1,
  1,15,105,455,1365,3003,5005, 6435, 6435, 5005,3003,1365, 455,105, 15, 1,
  1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368,1820,560,120,16,1,
  1,17,136,680,2380,6188,12376,19448,24310,24310,19448,12376,6188,2380,680,136,17,1,
  1,18,153,816,3060,8568,18564,31824,43758,48620,43758,31824,18564,8568,3060,816,153,18,1,
  1,19,171,969,3876,11628,27132,50388,75582,92378,92378,75582,50388,27132,11628,3876,969,171,19,1,
  1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1
} ;

#define _binomial(_m,_k) (_BINOMIALS[(_m)*((_m)+1)/2+(_k)])
#define IS_EVEN(_n) ((_n)==2*((_n)/2))

#define _ellip_inc_E(_t, _k) gsl_sf_ellint_E((_t), (_k), GSL_PREC_DOUBLE)
#define _ellip_inc_F(_t, _k) gsl_sf_ellint_F((_t), (_k), GSL_PREC_DOUBLE)

static gint in_polygon(gdouble *x, gint *idx, gint str, gint nx, gdouble *y)
  
/*
  test to check if y lies inside or outside polygon defined by x (see
  poly_quad_mn for definition of polygon)

  return -1 for outside polygon, -2 for (strictly) inside, i for point
  on edge i
 */

{
  gint i ;
  gdouble sgn, ont ;

  sgn = orient2d(&(x[idx[0]*str]), &(x[idx[1]*str]), y) ;
  if ( sgn == 0 ) return 0 ;

  for ( i = 1 ; i < nx-1 ; i ++ ) {
    ont = sgn*orient2d(&(x[idx[1]*str]), &(x[idx[i+1]*str]), y) ;
    if ( ont == 0 ) return 1 ;
    if ( ont < 0 ) return -1 ;
  }
  ont = sgn*orient2d(&(x[idx[nx-1]*str]), &(x[idx[0]*str]), y) ;
  if ( ont == 0 ) return nx-1 ;
  if ( ont < 0 ) return -1 ;

  return -2 ;
}

/*
  recursions and initializations all come from Gradshteyn & Ryzhik

  table ordering: (m,n) = (0,-1) (1,-1), (2,-1) followed by 
  (m,n) (m>=0,n>=0) at _table_offset(m+n)+n 
*/

static gdouble recursion_n_down(gint m, gint n, gint p,
				gdouble S, gdouble C, gdouble d,
				gdouble k2, gdouble kd2, 
				gint n4, gdouble In4,
				gint n2, gdouble In2)

{
  gdouble I = 0.0 ;
  gint i ;

  g_assert(n4 == n2 + 2) ;

  if ( n == n2 ) return In2 ;
  if ( n == n4 ) return In4 ;

  for ( i = n2 - 2 ; i >= n ; i -= 2 ) {
    I = (m+i+4+p)*k2*In4 - ((i+p+3)*k2 - (m+i+2)*kd2)*In2 -
      pow(S,m+1)*pow(C,i+1)*pow(d,p+2) ;
    I /= (i+1)*kd2 ;
    In4 = In2 ; In2 = I ;
  }

  return I ;
}

static gdouble recursion_n_up(gint m, gint n, gint p,
			      gdouble S, gdouble C, gdouble d,
			      gdouble k2, gdouble kd2, 
			      gint n4, gdouble In4,
			      gint n2, gdouble In2)

{
  gdouble I = 0.0 ;
  gint i ;

  g_assert(n4 == n2-2) ;

  if ( n == n2 ) return In2 ; 
  if ( n == n4 ) return In4 ;

  for ( i = n2 + 2 ; i <= n ; i += 2 ) {
    I = pow(S,m+1)*pow(C,i-3)*pow(d,p+2) + 
      ((i + p - 1)*k2 - (m + i - 2)*kd2)*In2 + (i - 3)*kd2*In4 ;
    I /= (m + i + p)*k2 ;
    In4 = In2 ; In2 = I ;
  }

  return I ;
}

static gdouble recursion_m_down(gint m, gint n, gint p,
				gdouble S, gdouble C, gdouble d,
				gdouble k2, gdouble kd2, 
				gint m4, gdouble Im4,
				gint m2, gdouble Im2)

{
  gdouble I = 0.0 ;
  gint i ;

  g_assert(m4 == m2+2) ;

  if ( m == m2 ) return Im2 ; 
  if ( m == m4 ) return Im4 ;

  for ( i = m2 - 2 ; i >= m ; i -= 2 ) {
    I = pow(S,i+1)*pow(C,n+1)*pow(d,p+2) + 
      ((i + n + 2) + (i + p + 3)*k2)*Im2 - (i+n+p+4)*k2*Im4 ;
    I /= (i + 1) ;
    Im4 = Im2 ; Im2 = I ;
  }

  return I ;
}

static gdouble recursion_m_up(gint m, gint n, gint p,
			      gdouble S, gdouble C, gdouble d,
			      gdouble k2, gdouble kd2, 
			      gint m4, gdouble Im4,
			      gint m2, gdouble Im2)

{
  gdouble I = 0.0 ;
  gint i ;

  g_assert(m4 == m2-2) ;

  if ( m == m2 ) return Im2 ; 
  if ( m == m4 ) return Im4 ;

  for ( i = m2 + 2 ; i <= m ; i += 2 ) {
    I = pow(S,i-3)*pow(C,n+1)*pow(d,p+2) + 
      ((i + n - 2) + (i + p - 1)*k2)*Im2 - (i - 3)*Im4 ;
    I /= (i + n + p)*k2 ;
    Im4 = Im2 ; Im2 = I ;
  }

  return I ;
}

static gint int_sin_cos_delta_init_p1(gint N, gdouble t,
				      gdouble S, gdouble C, gdouble d,
				      gdouble k, gdouble k2, 
				      gdouble kd, gdouble kd2,
				      gdouble *Imn)

{
  gint off ;
  gdouble F, E, aS, logCp, S2, C2, k4, sc, logSk, logkd ;
  
  aS = asin(k*S) ; S2 = S*S ;
  logSk = log((d+kd*S)/(d-kd*S)) ; logkd = log((d+kd)/(d-kd)) ;
  Imn[_table_index(-1,-1)] = 0.5*log((1.0-d)/(1.0+d)) + 
    0.5*kd*log((d+kd)/(d-kd)) ;
  Imn[_table_index(0,-1)] = 0.5*kd*logSk + k*aS ;
  Imn[_table_index(1,-1)] = -d + kd/2*logkd ;
  Imn[_table_index(2,-1)] = -0.5*d*S + (2.0*k2-1)/2/k*aS + 0.5*kd*logSk ;
  Imn[_table_index(3,-1)] = -(k2*S2+3*k2-1)/3/k2*d + 0.5*kd*logkd ;

  Imn[_table_index(-1,0)] = -0.5*log((d+C)/(d-C)) + k*log(k*(k*C+d)) ;
  Imn[_table_index(-1,1)] = d + 0.5*log((1.0-d)/(1.0+d)) ;
  Imn[_table_index(-1,2)] = 0.5*d*C + 0.5*(k2+1)/k*log(k*C+d) -
    0.5*log((d+C)/(d-C)) ;
  Imn[_table_index(-1,3)] = -(k2*S2 - 3*k2 - 1)/3/k2*d + 
    0.5*log((1.0-d)/(1.0+d)) ;

  off = _table_offset(0) ;
  E = _ellip_inc_E(t,k) ;
  Imn[off+0] = E ;
  if ( N == 0 ) return 0 ;

  off = _table_offset(1) ;
  logCp = log(d+k*C) ;
  Imn[off+0] = -d*C/2 - kd2/2/k*logCp ;
  Imn[off+1] = S*d/2 + 1.0/2/k*aS ;
  if ( N == 1 ) return 0 ;

  off = _table_offset(2) ;
  F = _ellip_inc_F(t,k) ;
  Imn[off+0] = (-d*S*C + (kd2*F + (2*k2-1)*E)/k2)/3 ;
  Imn[off+1] = -d*d*d/3/k2 ;
  Imn[off+2] = (d*S*C - (kd2*F - (k2+1)*E)/k2)/3 ;
  if ( N == 2 ) return 0 ;
  
  off = _table_offset(3) ;
  C2 = C*C ; sc = 1.0/8/k2 ;
  Imn[off+0] = (-(2*k2*S2+3*k2-1)*d*C + 
		(-1+k2*(-2+3*k2))/k*logCp)*sc ;
  Imn[off+1] = ((2*k2*S2-1)*d*S + 1.0/k*aS)*sc ;
  Imn[off+2] = (-(2*k2*C2+kd2)*d*C + kd2*kd2/k*log(k*C+d))*sc ;
  Imn[off+3] = ((2*k2*C2+2*k2+1)*d*S + (4*k2-1)/k*aS)*sc ;
  if ( N == 3 ) return 0 ;

  off = _table_offset(4) ;
  k4 = k2*k2 ; sc = 1.0/15/k2 ;
  Imn[off+0] = (-(3*k2*S2 + 4*k2 - 1)*d*S*C - 
		(2*(-1 + k2*(-1 +2*k2))*F - (-2 + k2*(-3 + 8*k2))*E)/k2)*sc ;
  Imn[off+1] = (-2 + k2*S2*(-1 + 3*k2*S2))/k2*d*sc ;
  Imn[off+2] = (-(3*k2*C2 - 2*k2 + 1)*d*S*C -
		(kd2*(1+kd2)*F - 2*(1 - k2 + k4)*E)/k2)*sc ;
  Imn[off+3] = -(3*k4*S2*S2 - k2*(5*k2+1)*S2 + 5*k2 - 2)/k2*d*sc ;
  Imn[off+4] = ((3*k2*C2 + 3*k2 + 1)*d*S*C +
		(2*kd2*(kd2 - 2*k2)*F + (-2 + k2*(7 + 3*k2))*E)/k2)*sc ;

  return 0 ;
}

static gint int_sin_cos_delta_init_0(gint N, gdouble t,
				      gdouble S, gdouble C, gdouble d,
				      gdouble k, gdouble k2, 
				      gdouble kd, gdouble kd2,
				      gdouble *Imn)

{
  gint off ;
  gdouble S2, C2 ;
  
  Imn[_table_index(-1,-1)] = log(S/C) ;
  Imn[_table_index(0,-1)] = log(sqrt((1.0+S)/(1.0-S))) ;
  Imn[_table_index(1,-1)] = -log(C) ;
  Imn[_table_index(2,-1)] = -S + log(tan(0.25*M_PI+0.5*t)) ;
  Imn[_table_index(3,-1)] = 0.5*C*C - log(C) ;

  Imn[_table_index(-1,0)] = log(tan(0.5*t)) ;
  Imn[_table_index(-1,1)] = log(S) ;
  Imn[_table_index(-1,2)] = C + log(tan(0.5*t)) ;
  Imn[_table_index(-1,3)] = 0.5*C*C + log(S) ;

  off = _table_offset(0) ;
  Imn[off+0] = t ;
  if ( N == 0 ) return 0 ;

  off = _table_offset(1) ;
  Imn[off+0] = -C ;
  Imn[off+1] = S ;
  if ( N == 1 ) return 0 ;

  off = _table_offset(2) ;
  Imn[off+0] = -0.5*S*C + 0.5*t ;
  Imn[off+1] = 0.5*S*S ;
  Imn[off+2] = 0.5*S*C + 0.5*t ;
  if ( N == 2 ) return 0 ;
  
  off = _table_offset(3) ;
  C2 = C*C ; S2 = S*S ;
  Imn[off+0] = C*C2/3.0 - C ;
  Imn[off+1] = S*S2/3.0 ;
  Imn[off+2] = -C*C2/3.0 ;
  Imn[off+3] = S - S*S2/3.0 ;
  if ( N == 3 ) return 0 ;

  off = _table_offset(4) ;
  Imn[off+0] = -0.375*S*C - 0.25*S*S2*C + 0.375*t ;
  Imn[off+1] = 0.25*S2*S2 ;
  Imn[off+2] = -0.125*C*S + 0.25*C*S*S2 + 0.125*t ;
  Imn[off+3] = -0.25*C2*C2 ;
  Imn[off+4] = 0.375*t + 0.375*S*C + 0.25*S*C*C2 ;

  return 0 ;
}


static gint int_sin_cos_delta_init_m1(gint N, gdouble t,
				      gdouble S, gdouble C, gdouble d,
				      gdouble k, gdouble k2, 
				      gdouble kd, gdouble kd2,
				      gdouble *Imn)

{
  gint off ;
  gdouble F, E, aS, logCp, logCm, sc, logSk ;
  
  aS = asin(k*S) ; logSk = log((d-kd*S)/(d+kd*S)) ;
  Imn[_table_index(-1,-1)] = 0.5*log((1.0-d)/(1.0+d)) + 
    0.5/kd*log((d+kd)/(d-kd)) ;
  Imn[_table_index(0,-1)] = -0.5/kd*logSk ;
  Imn[_table_index(1,-1)] = 0.5/kd*log((d+kd)/(d-kd)) ;
  Imn[_table_index(2,-1)] = -0.5/kd*logSk - aS/k ;
  Imn[_table_index(3,-1)] = d/k2 + 0.5/kd*log((d+kd)/(d-kd)) ;

  Imn[_table_index(-1,0)] = -0.5*log((d+C)/(d-C)) ;
  Imn[_table_index(-1,1)] = 0.5*log((1.0-d)/(1.0+d)) ;
  Imn[_table_index(-1,2)] = -0.5*log((d+C)/(d-C)) + log(k*C+d)/k ;
  Imn[_table_index(-1,3)] = d/k2 - 0.5*log((1.0+d)/(1.0-d));

  off = _table_offset(0) ;
  Imn[off+0] = F = _ellip_inc_F(t,k) ;
  if ( N == 0 ) return 0 ;

  off = _table_offset(1) ;
  sc = 1.0/k ;
  logCp = log(d+k*C) ; logCm = log(d-k*C) ;
  Imn[off+0] = 0.5*(logCm - logCp)*sc ;
  Imn[off+1] = aS*sc ;
  if ( N == 1 ) return 0 ;

  off = _table_offset(2) ;
  E = _ellip_inc_E(t,k) ; sc *= sc ;
  Imn[off+0] = (F - E)*sc ;
  Imn[off+1] = -d*sc ;
  Imn[off+2] = (E - kd2*F)*sc ;
  if ( N == 2 ) return 0 ;
  
  off = _table_offset(3) ;
  sc *= 0.5 ;
  Imn[off+0] = (C*d - (1+k2)/k*logCp)*sc ;
  Imn[off+1] = (-S*d + aS/k)*sc ;
  Imn[off+2] = (-C*d + kd2/k*logCp)*sc ;
  Imn[off+3] = (S*d + (2*k2 - 1)/k*aS)*sc ;
  if ( N == 3 ) return 0 ;

  off = _table_offset(4) ;
  sc *= 2.0/3 ;
  Imn[off+0] = (S*C*d + ((2+k2)*F - 2*(1+k2)*E)/k2)*sc ;
  Imn[off+1] = -3.0*(2+k2*S*S)*d*sc*sc ;
  Imn[off+2] = (-S*C*d + ((2-k2)*E + (2*k2-2)*F)/k2)*sc ;
  Imn[off+3] = -3.0*(k2*C*C - 2*kd2)*d*sc*sc ;
  Imn[off+4] = (S*C*d + ((4*k2-2)*E + (2 + k2*(-5 + 3*k2))*F)/k2)*sc ;

  return 0 ;
}

static gint int_sin_cos_delta_init(gint N, gint r, gdouble t,
				   gdouble S, gdouble C, gdouble d,
				   gdouble k, gdouble k2, 
				   gdouble kd, gdouble kd2,
				   gdouble *Imn)

{
  switch ( r ) {
  default: g_error("%s: case r == %d not implemented", __FUNCTION__, r) ;
  case -1: 
    return int_sin_cos_delta_init_m1(N, t, S, C, d, k, k2, kd, kd2, Imn) ;
    break ;
  case  0: 
    return int_sin_cos_delta_init_0(N, t, S, C, d, k, k2, kd, kd2, Imn) ;
    break ;
  case  1: 
    return int_sin_cos_delta_init_p1(N, t, S, C, d, k, k2, kd, kd2, Imn) ;
    break ;
  }

  return 0 ;
}

static gint int_sin_cos_delta_fill_mn(gint N, gint r, gdouble t,
				      gdouble S, gdouble C, gdouble d,
				      gdouble k, gdouble k2, 
				      gdouble kd, gdouble kd2,
				      gdouble *Imn)
/*
  Indefinite integral:

  \int^{t}\sin^{m} t \cos^{n} t (1-k^{2}\sin^{2}t)^{r/2} d t.

  r odd or 0, r >= -3; 0 <= m+n <= N

  Gradshteyn & Ryzhik 2.581 et seq.
*/

{
  int_sin_cos_delta_init(N, r, t, S, C, d, k, k2, kd, kd2, Imn) ;

  if ( N > 4 ) 
    int_sin_cos_delta_expand_mn(4, N, r, t, S, C, d, k, k2, kd, kd2, Imn) ;

  return 0 ;
}

gint int_sin_cos_delta_expand_mn(gint M, gint N, gint r, gdouble t,
				 gdouble S, gdouble C, gdouble d,
				 gdouble k, gdouble k2, 
				 gdouble kd, gdouble kd2,
				 gdouble *Imn)
/*
  Expand table of integrals using recursion relations (see
  int_sin_cos_delta_fill_mn)

  Indefinite integral:

  \int^{t}\sin^{m} t \cos^{n} t (1-k^{2}\sin^{2}t)^{r/2} d t.

  r odd or 0, r >= -3; 0 <= m+n <= N

  Gradshteyn & Ryzhik 2.581 et seq.
*/

{
  gint off, off2, off4, m, n, i ;
  gdouble cft, Cn ;

  /*recursion based on pre-existing table, on m for the
    first half of each row, then on n*/
  for ( i = M+1 ; i <= N ; i ++ ) {
    off = _table_offset(i) ; off2 = _table_offset(i-2) ; 
    off4 = _table_offset(i-4) ;

    Cn = C*pow(d,r+2) ;
    for ( n = 0 ; n <= i/2 ; n ++ ) {
      m = i - n ;
      cft = pow(S, m-3)*Cn ;
      Imn[off+n] = cft + 
	((m + n - 2) + (m + r - 1)*k2)*Imn[off2+n] - 
	(m - 3)*Imn[off4+n] ;
      Imn[off+n] /= (m + n + r)*k2 ;
      Cn *= C ;
    }
    Cn = pow(C,i/2+1-3)*pow(d,r+2) ;
    for ( n = i/2+1 ; n <= i ; n ++ ) {
      m = i - n ;
      cft = pow(S, m+1)*Cn ;
      Imn[off+n] = cft + 
	((n + r - 1)*k2 - (m + n - 2)*kd2)*Imn[off2+n-2] +
	(n - 3)*kd2*Imn[off4+n-4] ;
      Imn[off+n] /= (m + n + r)*k2 ;

      Cn *= C ;
    }
  }
  
  return 0 ;
}

gint int_sin_cos_workspace_init(int_sin_cos_workspace_t *w,
				gdouble t, gdouble k, gint N)

{
  g_assert(N <= w->Nmax) ;

  w->N = N ;
  w->t = t ; w->S = sin(t) ; w->C = cos(t) ; 
  w->k = k ; w->k2 = k*k ; w->kd2 = 1.0-k*k ; w->kd = sqrt(w->kd2) ;
  w->d = sqrt(1.0-k*k*w->S*w->S) ;

  g_assert(k >= 0 && k < 1) ;

  int_sin_cos_delta_fill_mn(N, -1, w->t, w->S, w->C, w->d, 
			    w->k, w->k2, w->kd, w->kd2, w->Im1) ;
  int_sin_cos_delta_fill_mn(N,  0, w->t, w->S, w->C, w->d, 
			    w->k, w->k2, w->kd, w->kd2, w->I0) ;
  int_sin_cos_delta_fill_mn(N,  1, w->t, w->S, w->C, w->d, 
			    w->k, w->k2, w->kd, w->kd2, w->Ip1) ;
  return 0 ;
}

int_sin_cos_workspace_t *int_sin_cos_workspace_alloc(gint Nmax)

{
  int_sin_cos_workspace_t *w ;

  w = (int_sin_cos_workspace_t *)g_malloc(sizeof(int_sin_cos_workspace_t)) ;

  w->Nmax = Nmax ;

  w->Im1 = (gdouble *)g_malloc(_table_offset(Nmax+1)*sizeof(gdouble)) ;
  w->I0 = (gdouble *)g_malloc(_table_offset(Nmax+1)*sizeof(gdouble)) ;
  w->Ip1 = (gdouble *)g_malloc(_table_offset(Nmax+1)*sizeof(gdouble)) ;

  return w ;
}

gdouble int_sin_cos_delta_mn(int_sin_cos_workspace_t *w,
			     gint m, gint n, gint p)

{
  gdouble I4, I2 ;
  gint m2, m4, n2, n4 ;

  /* fprintf(stderr, "%d %d %d\n", m, n, p) ; */

  /* g_assert(m >= 0) ; */
  g_assert(p == 0 || p != 2*(p/2)) ;
  g_assert(m+n <= w->Nmax) ;

  /* if ( n < -4 && w->k < 0.05 && p != 0 )  */
  /*   return int_sin_cos_delta_mn_series(w, m, n, p) ; */

  if ( p > 1 ) 
    return int_sin_cos_delta_mn(w, m, n, p-2) -
      w->k2*int_sin_cos_delta_mn(w, m+2, n, p-2) ;

  if ( p < -1 ) {
    g_assert( !IS_EVEN(p) ) ;
    I4 = int_sin_cos_delta_mn(w, m-2, n,   p+2) ;
    I2 = int_sin_cos_delta_mn(w, m,   n-2, p+2) ;
    return (-pow(w->S,m-1)*pow(w->C,n-1)*pow(w->d,p+2) +
	    I4*(m-1) - I2*(n-1))/(p+2)/w->k2 ;
  }

  if ( m+n > w->N ) {
    int_sin_cos_delta_expand_mn(w->N, m+n, -1, w->t, w->S, w->C, w->d, 
				w->k, w->k2, w->kd, w->kd2, w->Im1) ;
    int_sin_cos_delta_expand_mn(w->N, m+n,  0, w->t, w->S, w->C, w->d,
    				w->k, w->k2, w->kd, w->kd2, w->I0) ;
    int_sin_cos_delta_expand_mn(w->N, m+n,  1, w->t, w->S, w->C, w->d, 
				w->k, w->k2, w->kd, w->kd2, w->Ip1) ;
    w->N = m+n ;
  }

  if ( ((n >= 0 && m >= 0 && (m+n) <= w->N)) || 
       ((m <=  3 && m >= -1) && (n == -1)) ||
       ((m == -1) && (n <=  3 && n >= -1)) ) {
    if ( p == -1 ) return w->Im1[_table_index(m,n)] ;
    if ( p ==  0 ) return w->I0[_table_index(m,n)] ;
    if ( p ==  1 ) return w->Ip1[_table_index(m,n)] ;
  }

  if ( m == -1 ) {
    if ( n >= 0 ) {
      n4 = ( IS_EVEN(n) ? 0 : 1 ) ; n2 = n4 + 2 ;
      return recursion_n_up(m, n, p, w->S, w->C, w->d, w->k2, w->kd2, 
			    n4, int_sin_cos_delta_mn(w, m, n4, p),
			    n2, int_sin_cos_delta_mn(w, m, n2, p)) ;
    }
    n4 = ( IS_EVEN(n) ? 2 : 1 ) ; n2 = n4 - 2 ;

    return recursion_n_down(m, n, p, w->S, w->C, w->d, w->k2, w->kd2, 
			    n4, int_sin_cos_delta_mn(w, m, n4, p),
			    n2, int_sin_cos_delta_mn(w, m, n2, p)) ;    
  }

  if ( n == -1 ) {
    if ( m >= 0 ) {
      m4 = ( IS_EVEN(m) ? 0 : 1 ) ; m2 = m4 + 2 ;
      return recursion_m_up(m, n, p, w->S, w->C, w->d, w->k2, w->kd2,
			    m4, int_sin_cos_delta_mn(w, m4, n, p),
			    m2, int_sin_cos_delta_mn(w, m2, n, p)) ;
    }
    m4 = ( IS_EVEN(m) ? 2 : 1 ) ; m2 = m4 - 2 ;

    return recursion_m_down(m, n, p, w->S, w->C, w->d, w->k2, w->kd2, 
			    m4, int_sin_cos_delta_mn(w, m4, n, p),
			    m2, int_sin_cos_delta_mn(w, m2, n, p)) ;    
  }

  if ( n < 0 && m >= 0 ) {
    n4 = ( IS_EVEN(n) ? 2 : 1 ) ; n2 = n4 - 2 ;
    m4 = ( IS_EVEN(m) ? 0 : 1 ) ; m2 = m4 + 2 ;

    I2 = recursion_m_up(m, n2, p, w->S, w->C, w->d, w->k2, w->kd2, 
			m4, int_sin_cos_delta_mn(w, m4, n2, p),
			m2, int_sin_cos_delta_mn(w, m2, n2, p)) ;
    I4 = recursion_m_up(m, n4, p, w->S, w->C, w->d, w->k2, w->kd2, 
			m4, int_sin_cos_delta_mn(w, m4, n4, p),
			m2, int_sin_cos_delta_mn(w, m2, n4, p)) ;
    return recursion_n_down(m, n, p, w->S, w->C, w->d, w->k2, w->kd2, 
			    n4, I4, n2, I2) ;
  }

  if ( m < 0 && n >= 0 ) {
    m4 = ( IS_EVEN(m) ? 2 : 1 ) ; m2 = m4 - 2 ;

    return recursion_m_down(m, n, p, w->S, w->C, w->d, w->k2, w->kd2, 
			    m4, int_sin_cos_delta_mn(w, m4, n, p),
			    m2, int_sin_cos_delta_mn(w, m2, n, p)) ;
  }

  g_assert(m < -1 && n < -1) ;

  m4 = ( IS_EVEN(m) ? 2 : 1 ) ; m2 = m4 - 2 ;
  n4 = ( IS_EVEN(n) ? 2 : 1 ) ; n2 = n4 - 2 ;

  I4 = recursion_m_down(m, n4, p, w->S, w->C, w->d, w->k2, w->kd2, 
			m4, int_sin_cos_delta_mn(w, m4, n4, p),
			m2, int_sin_cos_delta_mn(w, m2, n4, p)) ;
  I2 = recursion_m_down(m, n2, p, w->S, w->C, w->d, w->k2, w->kd2, 
			m4, int_sin_cos_delta_mn(w, m4, n2, p),
			m2, int_sin_cos_delta_mn(w, m2, n2, p)) ;

  return recursion_n_down(m, n, p, w->S, w->C, w->d, w->k2, w->kd2, 
			  n4, I4, n2, I2) ;

  return 0.0 ;
}

static gdouble quad_mn_log(gint m, gint n, int_sin_cos_workspace_t *w)

/*indefinite integral 
  \int^{t}\sin^{m} t \cos^{n} t\log((\Delta+\alpha')/(\Delta-\alpha')) dt*/

{
  gdouble I, L, d, kd, S, C ;
  gint q, md, nd, sgn ;

  if ( m+n == 2*((m+n)/2) ) 
    g_error("%s: m+n must be even (%d+%d=%d)", __FUNCTION__, m, n, m+n) ;

  S = w->S ; C = w->C ; d = w->d ; kd = w->kd ;
  /* k = w->k ; */
  L = log((d+kd)/(d-kd)) ;

  I = 0.0 ; sgn = 1 ;
  if ( m != 2*(m/2) ) {
    md = m/2 ;
    for ( q = 0 ; q <= md ; (q ++), (sgn = -sgn) ) {
      I += (gdouble)_binomial(md,q)*sgn/(n+2*q+1)*
	(-L*pow(C,n+2*q+1) + 
	 2*kd*int_sin_cos_delta_mn(w, 1, n+2*q, -1)) ;
    }
    return I ;
  }

  nd = n/2 ;
  for ( q = 0 ; q <= nd ; (q ++), (sgn = -sgn) ) {
    I += (gdouble)_binomial(nd,q)*sgn/(m+2*q+1)*
      (L*pow(S,m+2*q+1) -
       2*kd*int_sin_cos_delta_mn(w, m+2*q+2, -1, -1)) ;
  }

  return I ;
}

static gdouble quad_mn_Rp(gint m, gint n, gint p, gdouble z,
			  int_sin_cos_workspace_t *w)

/*
  indefinite integral 
  \int^{t} \sin^{m} t \cos^{n} t R^{p} dt

  m+n odd, p odd, p>= -1
*/

{
  gint u, q ;
  gdouble I, L, aq, aL, bt ;

  g_assert(m+n != 2*((m+n)/2)) ;
  g_assert(p != 2*(p/2)) ;
  /* g_assert(p >= -1) ; */

  bt = z/w->k ;

  if ( p >= -1 ) {
    u = (p-1)/2 ;
    L = quad_mn_log(m, n, w)*0.5 ;
    if ( u == -1 ) return L ;
    
    aq = 1.0 ; aL = 1.0 ;
    I = int_sin_cos_delta_mn(w, m, n-2*u-2, 2*u+1)*pow(bt,2*u+2) ;
    for ( q = 0 ; q <= u-1 ; q ++ ) {
      aq *= (gdouble)(2*u-2*q+1)/(u-q)*z*z/2 ;
      I += aq*pow(bt,2*u-2*q)*int_sin_cos_delta_mn(w, m, n-2*u+2*q,2*u-2*q-1) ;
    }
    
    I *= w->kd/(u+1)*0.5 ;
    
    for ( q = 0 ; q <= u ; q ++ ) aL *= (gdouble)(2*u+1-2*q)/(u+1-q) ;
    
    I += aL*L*pow(z*z/2,u+1) ;

    return I ;
  }

  /*Gradshteyn & Ryzhik, 2.263.4*/
  u = -(p-1)/2 - 1 ; aq = 1.0 ;

  I = int_sin_cos_delta_mn(w, m, n+2*u-2, -(2*u-1))/pow(bt, 2*u-1) ;
  for ( q = 1 ; q <= u-1 ; q ++ ) {
    aq *= (u-q)/(2*u-2*q-1)*(2/z/z) ;
    I += 
      aq*int_sin_cos_delta_mn(w, m, n+2*u-2*q-2, -(2*u-2*q-1))/
      pow(bt, 2*u-2*q-1) ;
  }

  I *= bt*w->kd/z/z/(2*u-1) ;

  return I ;
}

gint subtriangle_quad_mn(gdouble r1, gdouble r2, gdouble th,
			 gdouble psi, gdouble z,
			 gint N, gint gm, gdouble *I,
			 int_sin_cos_workspace_t *w0,
			 int_sin_cos_workspace_t *w1)

/*
  integral of x_{1}^{n}y_{1}^{m}R^{\gamma} over triangle 
  (0,0), r_{1}(\cos(\psi),\sin(\psi)), 
  r_{2}(\cos(\theta+\psi),\sin(\theta+\psi)),

  field point at (0,0,z)

  orientation of triangle vertices is respected through sign of \theta
*/

/*need to fix problems with high order terms at small z*/

{
  gdouble al, bt, phi, a, tmp, sgn, buf[64], Cp, Sp, z2k ;
  gint i, m, n, k, q, s, off ;

  a = (r2*cos(th) - r1)/r2/sin(th) ;
  phi = atan(a) ;

  bt = (r1*r1 + z*z*(1.0+a*a))/(1.0+a*a) ;
  al = sqrt(z*z/bt) ;
  bt = sqrt(bt) ;

  if ( al > 1.0 - 1e-8 ) {
    for ( i = 0 ; i <= N ; i ++ ) {
      off = i*(i+1)/2 ;
      for ( n = 0 ; n <= i ; n ++ ) 
	I[off+n] = 0.0 ;
    }
    return 0 ;
  }

  int_sin_cos_workspace_init(w1, th+phi, al, N) ;
  int_sin_cos_workspace_init(w0, phi,    al, N) ; 

  for ( q = 0 ; q <= N/2 ; q ++ ) {
    /*m+n even to start*/
    off = (2*q)*(2*q+1)/2 ;
    for ( n = 0 ; n <= 2*q ; n ++ ) {
      m = 2*q - n ;
      I[off+n] = 0.0 ; z2k = 1.0 ; sgn = 1 ;
      for ( k = 0 ; k <= q ; (k ++), (sgn = -sgn) ) {
	tmp = (int_sin_cos_delta_mn(w1, m, n+2*k-2*q-gm-2, 2*q-2*k+gm+2) -
	       int_sin_cos_delta_mn(w0, m, n+2*k-2*q-gm-2, 2*q-2*k+gm+2))*
	  pow(bt,2*q-2*k+gm+2) ;
	tmp -= (int_sin_cos_delta_mn(w1, m, n, 0) -
		int_sin_cos_delta_mn(w0, m, n, 0))*pow(fabs(z),2*q-2*k+gm+2) ;
	I[off+n] += 1.0/(2*q-2*k+gm+2)*_binomial(q,k)*tmp*z2k ;
	z2k *= -z*z ;
      }
    }
    /*m+n odd*/
    off = (2*q+1)*(2*q+2)/2 ;
    for ( n = 0 ; n <= 2*q+1 ; n ++ ) {
      m = 2*q+1 - n ;
      I[off+n] = 0.0 ; z2k = 1.0 ;
      for ( k = q+1 ; k >= 0 ; k -- ) {
	I[off+n] += 
	  _binomial(q+1,k)*(quad_mn_Rp(m, n, gm+2*k, z, w1) -
			    quad_mn_Rp(m, n, gm+2*k, z, w0))*z2k ;
	z2k *= -z*z ;
      }
    }
  }
  
  Cp = cos(phi-psi) ; Sp = sin(phi-psi) ;
  /*compute the monomial integrals over the source triangle*/
  for ( i = N ; i >= 0 ; i -- ) {
    off = i*(i+1)/2 ;
    memcpy(buf, &(I[off]), (i+1)*sizeof(gdouble)) ;
    memset(&(I[off]), 0, (i+1)*sizeof(gdouble)) ;
    for ( n = 0 ; n <= i ; n ++ ) {
      m = i - n ; sgn = 1 ;
      for ( q = 0 ; q <= m ; (q ++), (sgn = -sgn) ) {
	for ( s = 0 ; s <= n ; s ++ ) {
	  I[off+n] += _binomial(m,q)*_binomial(n,s)*sgn*
	    pow(Cp,i-q-s)*pow(Sp,q+s)*buf[n+q-s] ;
	}
      }
    }
  }

  return 0 ;
}

gint poly_quad_mn(gdouble *x, gint *idx, gint str, 
		  gint nx, gdouble z, gint gm, gint N,
		  gdouble *I, gdouble *work,
		  int_sin_cos_workspace_t *w0,
		  int_sin_cos_workspace_t *w1)

/*
  Integrals of x^{n}y^{m}R^{\gamma} (n+m <= N) over polygon x[idx[0]],
  x[idx[1]], ..., x[idx[nx-1]], x[idx[0]]

  vertex coordinates are x[idx[i]*str+0,1]

  integrals are packed in I as I[(n+m)*(n+m+1)/2+n]
*/

{
  gdouble *buf, *rt, xo[2] = {0,0} ;
  gint i, j, nc, inp ;

  rt = work ; buf = &(work[2*(nx+2)]) ;

  for ( i = 0 ; i < nx ; i ++ ) {
    rt[2*i+0] = sqrt(x[idx[i]*str+0]*x[idx[i]*str+0] +
		     x[idx[i]*str+1]*x[idx[i]*str+1]) ;
    rt[2*i+1] = atan2(x[idx[i]*str+1], x[idx[i]*str+0]) ;
    if ( rt[2*i+1] < 0 ) rt[2*i+1] += 2*M_PI ;
  }

  rt[2*nx+0] = rt[2*0+0] ; rt[2*nx+1] = rt[2*0+1] ;

  /*this needs a test to set up the angles properly, e.g. for
    (0,0)--(0,-1)--(1,0)*/
  if ( (inp = in_polygon(x, idx, str, nx, xo)) == -2 ) 
    rt[2*nx+1] += 2*M_PI ;
  else {
    /*check for being on an edge*/
    /* g_assert(inp < 0) ; */
  }
  

  nc = (N+1)*(N+2)/2 ; memset(I, 0, nc*sizeof(gdouble)) ;
  for ( i = 0 ; i < nx ; i ++ ) {
    /* if ( i != inp ) { */
    if ( (rt[2*(i+1)+1] != rt[2*i+1]) && 
	   (rt[2*i+0] != 0) && ( rt[2*(i+1)+0] != 0) ) {
	subtriangle_quad_mn(rt[2*i+0], rt[2*(i+1)+0], 
			    rt[2*(i+1)+1]-rt[2*i+1], rt[2*i+1],
			    z, N, gm, buf, w0, w1) ;
	for ( j = 0 ; j < nc ; j ++ ) I[j] += buf[j] ;
    }
    /*   if ( inp >= 0 ) fprintf(stderr, "X %lg\n", I[0]) ; */
    /* } */
  }

  return 0 ;
}

/* static gdouble quad_z(gdouble z, gint n) */

/* { */
/*   fprintf(stderr, "Hello\n") ; */
/*   if ( n == -1 ) { */
/*     return z*z*(log(fabs(z))-1) ; */
/*     /\* return -z*z - z*log(fabs(z)) ; *\/ */
/* /\* log(fabs(z)) ; *\/ */
/*   } */

/*   g_assert(n > -1) ; */

/*   return pow(z,n+1)/(n+1) ; */
/* } */

gint pyramid_quad_mn(gdouble *x, gint *idx, gint str, 
		     gint nx, gdouble z, gint gm, gint N,
		     gdouble *I, gdouble *work,
		     int_sin_cos_workspace_t *w0,
		     int_sin_cos_workspace_t *w1)

/*
  Integrals of x^{n}y^{m}z^{p}R^{\gamma} (n+m <= N) over pyramid with
  base formed by polygon x[idx[0]], x[idx[1]], ..., x[idx[nx-1]],
  x[idx[0]], in the plane z=0, and apex at (0,0,z)

  base vertex coordinates are (x[idx[i]*str+0,1],0.0)

  integrals are found from I as 
  I[(n+m)*(n+m+1)/2+n]*I[offset + p*(N+gm+2+1) + (n+m)+gm+2]

  offset = (N+1)*(N+2)/2 + 128
  p = 0, ..., N
*/

{
  gint off, p, q ;
  gdouble cft, zp ;

  /*not worked out yet for higher order singularities*/
  /* g_assert(gm >= -1) ; */

  /*integrals on the polygon base*/
  poly_quad_mn(x, idx, str, nx, z, gm, N, I, work, w0, w1) ;
  
  /*scaling factors for powers of z*/
  off = (N+1)*(N+2)/2 + 128 ;

  zp = 1.0 ;
  for ( p = 0 ; p <= N+1 ; p ++ ) {
    cft = 1.0/(p+1) ; zp *= z ;
    for ( q = 0 ; q <= N + gm + 2 ; q ++ ) {
      I[off + p*(N + gm + 2 + 1) + q] = zp*cft ;
      cft *= (gdouble)(q+1)/(gdouble)(p+q+2) ;
    }
  }

  p = 0 ; q = -1 ;
  I[off + p*(N + gm + 2 +1) + q] = z*log(fabs(z)) ;
/* z*z*(log(fabs(z))-1) ; */
/* -z*fabs(z)*(1-log(fabs(z))) ; */

  return 0 ;
}

gdouble pyramid_quad_general(gint gm, gint N, gdouble *I,
			     gdouble o[],
			     gdouble s[], gdouble t[], gdouble n[],
			     gint mx, gint my, gint mz, 
			     gdouble z, gint grad)

/*
  convert integrals on pyramid in reference position to integrals in
  global coordinates using coordinate system from tet_projection()
 */

{
  gdouble Im ;
  gint ix, jx, kx, iy, jy, ky, iz, jz, kz, off ;

  Im = 0.0 ;

  off = (N+1)*(N+2)/2 + 128 ;
  for ( ix = 0 ; ix <= mx ; ix ++ ) {
    for ( jx = 0 ; jx <= ix ; jx ++ ) {
      for ( kx = 0 ; kx <= mx-ix ; kx ++ ) {
	for ( iy = 0 ; iy <= my ; iy ++ ) {
	  for ( jy = 0 ; jy <= iy ; jy ++ ) {
	    for ( ky = 0 ; ky <= my-iy ; ky ++ ) {
	      for ( iz = 0 ; iz <= mz ; iz ++ ) {
		for ( jz = 0 ; jz <= iz ; jz ++ ) {
		  for ( kz = 0 ; kz <= mz-iz ; kz ++ ) {
		    Im += 
		      _binomial(mx,ix)*_binomial(ix,jx)*_binomial(mx-ix,kx)*
		      _binomial(my,iy)*_binomial(iy,jy)*_binomial(my-iy,ky)*
		      _binomial(mz,iz)*_binomial(iz,jz)*_binomial(mz-iz,kz)*
		      pow(s[0],jx)*pow(t[0],ix-jx)*pow(n[0],kx)*
		      pow(o[0],mx-ix-kx)*
		      pow(s[1],jy)*pow(t[1],iy-jy)*pow(n[1],ky)*
		      pow(o[1],my-iy-ky)*
		      pow(s[2],jz)*pow(t[2],iz-jz)*pow(n[2],kz)*
		      pow(o[2],mz-iz-kz)*
		      I[((ix+iy+iz)*(ix+iy+iz+1)/2)+jx+jy+jz]*
		      I[off+(kx+ky+kz)*(N+gm+2+1) + (ix+iy+iz)+gm+2] ;
		  } /*kz*/
		} /*jz*/
	      } /*iz*/
	    } /*ky*/
	  } /*jy*/
	} /*iy*/
      } /*kx*/
    } /*jx*/
  } /*ix*/

  return Im ;
}

gdouble polygon_quad_general(gint gm, gint N, gdouble *I,
			     gdouble o[],
			     gdouble s[], gdouble t[], gdouble n[],
			     gint mx, gint my, gint mz)

/*
  convert integrals on polygon in reference position to integrals in
  global coordinates using coordinate system from tet_projection()
 */

{
  gdouble Im ;
  gint ix, jx ;

  Im = 0.0 ;

  /* off = (N+1)*(N+2)/2 + 128 ; */
  for ( ix = 0 ; ix <= mx ; ix ++ ) {
    for ( jx = 0 ; jx <= ix ; jx ++ ) {
      Im += _binomial(mx,ix)*_binomial(ix,jx)*
	pow(s[0],jx)*pow(t[0],ix-jx)*pow(o[0],mx-ix)*
	I[(ix*(ix+1)/2)+jx] ;
    } /*jx*/
  } /*ix*/

  return Im ;
}

gint tet_projection(gdouble *x0, gdouble *x1, gdouble *x2, gdouble *x3,
		    gdouble *o, gdouble *s, gdouble *t, gdouble *n,
		    gdouble *y1, gdouble *y2, gdouble *y3,
		    gdouble *z)

/*
  Transformation of general tetrahedron to reference coordinate
  system, (s,t,n) with
  x = o + y[0]*s + y[1]*t + y[2]*n ;
  x1,2,3 map to y1,2,3, x0 maps to o + z*n
 */

{
  gdouble len ;

  s[0] = x3[0] - x1[0] ; s[1] = x3[1] - x1[1] ; s[2] = x3[2] - x1[2] ;
  t[0] = x2[0] - x1[0] ; t[1] = x2[1] - x1[1] ; t[2] = x2[2] - x1[2] ;

  n[0] = t[1]*s[2] - t[2]*s[1] ; n[1] = t[2]*s[0] - t[0]*s[2] ;
  n[2] = t[0]*s[1] - t[1]*s[0] ;

  t[0] = n[1]*s[2] - n[2]*s[1] ; t[1] = n[2]*s[0] - n[0]*s[2] ;
  t[2] = n[0]*s[1] - n[1]*s[0] ;
  
  len = sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]) ;
  s[0] /= len ; s[1] /= len ; s[2] /= len ;

  len = sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]) ;
  t[0] /= len ; t[1] /= len ; t[2] /= len ;

  len = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]) ;
  n[0] /= len ; n[1] /= len ; n[2] /= len ;

  *z = (x0[0]-x1[0])*n[0] + (x0[1]-x1[1])*n[1] + (x0[2]-x1[2])*n[2] ;
  o[0] = x0[0] - *z*n[0] ; o[1] = x0[1] - *z*n[1] ; o[2] = x0[2] - *z*n[2] ;

  y1[0] = (x1[0]-o[0])*s[0] + (x1[1]-o[1])*s[1] + (x1[2]-o[2])*s[2] ;
  y1[1] = (x1[0]-o[0])*t[0] + (x1[1]-o[1])*t[1] + (x1[2]-o[2])*t[2] ;
  y1[2] = (x1[0]-o[0])*n[0] + (x1[1]-o[1])*n[1] + (x1[2]-o[2])*n[2] ;

  y2[0] = (x2[0]-o[0])*s[0] + (x2[1]-o[1])*s[1] + (x2[2]-o[2])*s[2] ;
  y2[1] = (x2[0]-o[0])*t[0] + (x2[1]-o[1])*t[1] + (x2[2]-o[2])*t[2] ;
  y2[2] = (x2[0]-o[0])*n[0] + (x2[1]-o[1])*n[1] + (x2[2]-o[2])*n[2] ;

  y3[0] = (x3[0]-o[0])*s[0] + (x3[1]-o[1])*s[1] + (x3[2]-o[2])*s[2] ;
  y3[1] = (x3[0]-o[0])*t[0] + (x3[1]-o[1])*t[1] + (x3[2]-o[2])*t[2] ;
  y3[2] = (x3[0]-o[0])*n[0] + (x3[1]-o[1])*n[1] + (x3[2]-o[2])*n[2] ;

  return 0 ;
}

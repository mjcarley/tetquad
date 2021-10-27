#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <gqr.h>

#include <blaswrap.h>

#include "duffy.h"

/*
 * Duffy transformation methods for integration on tetrahedra from 
 * 
 * Jia-He Lv, Yu-Yong Jiao, Xia-Ting Feng, Peter Wriggers, Xiao-Ying
 * Zhuang, Timon Rabczuk, A series of Duffy-distance transformation
 * for integrating 2D and 3D vertex singularities.
 * 
 * DOI: 10.1002/nme.6016
 */

gdouble orient3d(gdouble *pa, gdouble *pb, gdouble *pc, gdouble *pd);

static gint plane_project(gdouble *x1, gdouble *x2, gdouble *x3, gdouble *x,
			  gdouble *p)

/*
 * project point x onto the plane defined by x123
 */
  
{
  gdouble n[3], r1[3], r2[3], t ;

  duffy_vector3(r1,x1,x2) ;
  duffy_vector3(r2,x2,x3) ;
  duffy_cross3(n,r1,r2) ;

  t = sqrt(duffy_length3sqr(n)) ;
  n[0] /= t ; n[1] /= t ; n[2] /= t ; 

  t = duffy_diff_scalar3(n,x,x1) ;

  p[0] = x[0] + t*n[0] ; 
  p[1] = x[1] + t*n[1] ; 
  p[2] = x[2] + t*n[2] ; 
  
  return 0 ;
}

static gint interp(gdouble *x1, gdouble *x2, gdouble *x3, gdouble *x4,
		   gdouble xi, gdouble eta, gdouble zeta, gdouble *x)

{
  gdouble a[4] ;
  
  a[0] = 1.0 - xi ; a[1] = xi - eta ; a[2] = eta - zeta ; a[3] = zeta ;
  x[0] = a[0]*x1[0] + a[1]*x2[0] + a[2]*x3[0] + a[3]*x4[0] ;
  x[1] = a[0]*x1[1] + a[1]*x2[1] + a[2]*x3[1] + a[3]*x4[1] ;
  x[2] = a[0]*x1[2] + a[1]*x2[2] + a[2]*x3[2] + a[3]*x4[2] ;
  
  return 0 ;
}

/* static gboolean point_in_tet(gdouble *x1, gdouble *x2, */
/* 			     gdouble *x3, gdouble *x4, */
/* 			     gdouble *p) */

/* { */
/*   gint i, s ; */
/*   gdouble min, max ; */

/*   /\*bounding box test here to eliminate most predicate calls*\/ */
/*   for ( i = 0 ; i < 3 ; i ++ ) { */
/*     min = MIN(x1[i], MIN(x2[i], MIN(x3[i], x4[i]))) ; */
/*     max = MAX(x1[i], MAX(x2[i], MAX(x3[i], x4[i]))) ; */
/*     if ( p[i] < min ) return FALSE ; */
/*     if ( p[i] > max ) return FALSE ; */
/*   } */
  
/*   /\*add a circumsphere test here to speed things up a bit?*\/ */
/*   if ( (s = SIGN(orient3d(x1, x2, x3, p))) == 0 ) */
/*     return TRUE ;   */
/*   if ( s != SIGN(orient3d(x1, x2, x3, x4)) ) */
/*     return FALSE ; */

/*   if ( (s = SIGN(orient3d(x2, x3, x4, p))) == 0 ) */
/*     return TRUE ;   */
/*   if ( s != SIGN(orient3d(x2, x3, x4, x1)) ) */
/*     return FALSE ; */

/*   if ( (s = SIGN(orient3d(x3, x4, x1, p))) == 0 ) */
/*     return TRUE ;   */
/*   if ( s != SIGN(orient3d(x3, x4, x1, x2)) ) */
/*     return FALSE ; */

/*   if ( (s = SIGN(orient3d(x4, x1, x2, p))) == 0 ) */
/*     return TRUE ;   */
/*   if ( s != SIGN(orient3d(x4, x1, x2, x3)) ) */
/*     return FALSE ; */
  
/*   return TRUE ; */
/* } */

gint duffy_ref_tet_quad_plain(gdouble *x1, gdouble *x2,
			      gdouble *x3, gdouble *x4,
			      gdouble V,
			      gdouble al, gdouble bt,
			      gdouble *qu, gint ustr,
			      gdouble *wu, gint wustr, gint nu,
			      gdouble *qv, gint vstr,
			      gdouble *wv, gint wvstr, gint nv,
			      gdouble *qw, gint wstr,
			      gdouble *ww, gint wwstr, gint nw,
			      duffy_func_t func, gpointer data,
			      gdouble *q, gint nq)

{
  gint i, j, k ;
  gdouble ub, du, vb, dv, wb, dw, u, v, w, y[3], yT[3], r1T[3] ;
  gdouble xi, eta, zeta, l1T, J, wt ;
  
  ub = 0.5 ; du = 0.5 ;
  vb = 0.5 ; dv = 0.5 ;
  wb = 0.5 ; dw = 0.5 ;
  
  for ( j = 0 ; j < nv ; j ++ ) {
    v = vb + dv*qv[j*vstr] ;
	
    for ( k = 0 ; k < nw ; k ++ ) {
      w = wb + dw*qw[k*wstr] ;

      yT[0] = (1.0-v)*x2[0] + v*(1.0 - w)*x3[0] + v*w*x4[0] ; 
      yT[1] = (1.0-v)*x2[1] + v*(1.0 - w)*x3[1] + v*w*x4[1] ; 
      yT[2] = (1.0-v)*x2[2] + v*(1.0 - w)*x3[2] + v*w*x4[2] ; 

      duffy_vector3(r1T, x1, yT) ;
      l1T = sqrt(duffy_length3sqr(r1T)) ;

      for ( i = 0 ; i < nu ; i ++ ) {
      	u = ub + du*qu[i*ustr] ;
	J = V*bt*pow(u, 3*bt-1-al*bt)*v*du*dv*dw ;
	xi = pow(u, bt) ;
	eta = xi*v ;
	zeta = eta*w ;

	interp(x1, x2, x3, x4, xi, eta, zeta, y) ;

	wt = J*wu[i*wustr]*wv[j*wvstr]*ww[k*wwstr]/pow(l1T,al) ;
	func(y, x1, wt, q, nq, data) ;
      }
    }
  }
  
  return 0 ;
}

gint duffy_ref_tet_quad_vector(gdouble *x1, gdouble *x2,
			       gdouble *x3, gdouble *x4,
			       gdouble V,
			       gdouble al, gdouble bt,
			       gdouble *qu, gint ustr,
			       gdouble *wu, gint wustr, gint nu,
			       gdouble *qv, gint vstr,
			       gdouble *wv, gint wvstr, gint nv,
			       gdouble *qw, gint wstr,
			       gdouble *ww, gint wwstr, gint nw,
			       duffy_func_t func, gpointer data,
			       gdouble *q, gint nq)

{
  gint i, j, k ;
  gdouble ub, du, vb, dv, wb, dw, u, v, y[3], r1T[3*4] ;
  gdouble xi, eta, zeta, l1T, J, wt ;
  gdouble bA[9], vw[3*4] ;
  gdouble d1 = 1.0, d0 = 0.0 ;
  gint i1 = 1, i3 = 3, i4 = 4, m ;

  ub = 0.5 ; du = 0.5 ;
  vb = 0.5 ; dv = 0.5 ;
  wb = 0.5 ; dw = 0.5 ;

  /*base interpolation matrix*/
  /* bA[0] = x2[0] ; bA[1] = x3[0] - x2[0] ; bA[2] = x4[0] - x3[0] ;  */
  /* bA[3] = x2[1] ; bA[4] = x3[1] - x2[1] ; bA[5] = x4[1] - x3[1] ;  */
  /* bA[6] = x2[2] ; bA[7] = x3[2] - x2[2] ; bA[8] = x4[2] - x3[2] ;  */
  bA[0] = x1[0] - x2[0] ; bA[3] = x2[0] - x3[0] ; bA[6] = x3[0] - x4[0] ; 
  bA[1] = x1[1] - x2[1] ; bA[4] = x2[1] - x3[1] ; bA[7] = x3[1] - x4[1] ; 
  bA[2] = x1[2] - x2[2] ; bA[5] = x2[2] - x3[2] ; bA[8] = x3[2] - x4[2] ; 

  /* vw[0] = vw[1] = vw[2] = vw[3] = 1.0 ; */
  vw[0] = vw[3] = vw[6] = vw[9] = 1.0 ;
  for ( j = 0 ; j < nv ; j ++ ) {
    v = vb + dv*qv[j*vstr] ;

    vw[1] = vw[4] = vw[7] = vw[10] = v ;
    
    for ( k = 0 ; k < nw ; k += 4 ) {
      /* w = wb + dw*qw[k*wstr] ; */

      blaswrap_daxpy(i4, dw, &(qw[k]), i1, &(vw[2]), i3) ;
      /* vw[8] += wb ; vw[9] += wb ; vw[10] += wb ; vw[11] += wb ;  */
      vw[2] += wb ; vw[5] += wb ; vw[8] += wb ; vw[11] += wb ; 
      blaswrap_dscal(i4, v, &(vw[2]), i3) ;
      
      blaswrap_dgemm(FALSE, FALSE, i4, i3, i3, d1, vw, i3, bA, i3, d0,
		     r1T, i3) ;

      for ( m = 0 ; m < 4 ; m ++ ) {
	l1T = sqrt(duffy_length3sqr(&(r1T[3*m]))) ;

	for ( i = 0 ; i < nu ; i ++ ) {
	  u = ub + du*qu[i*ustr] ;
	  J = V*bt*pow(u, 3*bt-1-al*bt)*v*du*dv*dw ;
	  xi = pow(u, bt) ;
	  /* eta = xi*vw[4] ; */
	  /* zeta = eta*vw[8+m] ; */
	  eta = xi*vw[1] ;
	  /* zeta = eta*vw[8+m] ; */
	  zeta = xi*vw[3*m+2] ;
	  J = V*bt*pow(u, 3*bt-1-al*bt)*vw[1]*du*dv*dw ;
	  
	  interp(x1, x2, x3, x4, xi, eta, zeta, y) ;
	  
	  /* wt = J*wu[i*wustr]*wv[j*wvstr]*ww[k*wwstr]/pow(l1T,al) ; */
	  wt = J*wu[i*wustr]*wv[j*wvstr]*ww[(k+m)*wwstr]/pow(l1T,al) ;
	  func(y, x1, wt, q, nq, data) ;
	}
      }
    }
  }
  
  return 0 ;
}


gint duffy_ref_tet_quad_ch(gdouble *x1, gdouble *x2,
			   gdouble *x3, gdouble *x4,
			   gdouble V,
			   gdouble al, gdouble bt,
			   gdouble *qu, gint ustr,
			   gdouble *wu, gint wustr, gint nu,
			   gdouble *qv, gint vstr, gdouble *wv,
			   gint wvstr, gint nv,
			   gdouble *qw, gint wstr, gdouble *ww,
			   gint wwstr, gint nw,
			   duffy_func_t func, gpointer data,
			   gdouble *q, gint nq)
/*
 * Section 4.2: Near singularity hidden in cell height 
 */
  
{
  gint i, j, k ;
  gdouble ub, du, rb, dr, tb, dt, u, v, w, y[3] ;
  gdouble xi, eta, zeta, r12, d2, u1[3], u2[3] ;
  gdouble th, rh, rhbar, A, C, S, rbmin, rbmax ;
  gdouble J, Jw, lu1, lu2, u12, wt ;

  ub = 0.5 ; du = 0.5 ;
  tb = M_PI/8.0 ; dt = M_PI/8.0 ;

  duffy_vector3(u1, x2, x3) ;
  duffy_vector3(u2, x3, x4) ;

  lu1 = duffy_length3sqr(u1) ;
  lu2 = duffy_length3sqr(u2) ;
  u12 = duffy_scalar3(u1,u2) ;
  r12 = duffy_distance3sqr(x1,x2) ;

  for ( j = 0 ; j < nv ; j ++ ) {
    th = tb + dt*qv[j*vstr] ;
    C = cos(th) ; S = sin(th) ;

    A = lu1*C*C + lu2*S*S + 2.0*u12*C*S ;
    d2 = r12/A ;
    A = sqrt(A) ;
    
    rbmin = 0.5*log(d2) ;
    /* rbmax = 0.5*log(1.0/C/C + d2) ; */
    rbmax = 0.5*log(1.0/C/C + d2) ;
    rb = 0.5*(rbmax + rbmin) ; dr = 0.5*(rbmax - rbmin) ;
    /* g_assert(rbmax > rbmin) ; */
    
    for ( k = 0 ; k < nw ; k ++ ) {
      rhbar = rb + dr*qw[k*wstr] ;
      rh = sqrt(exp(2.0*rhbar) - d2) ;
      v = rh*C ; w = rh*S ;

      Jw = exp((2.0-al)*rhbar)*dr*V*bt*du*dt ;

      for ( i = 0 ; i < nu ; i ++ ) {
      	u = ub + du*qu[i*ustr] ;
	xi = pow(u, bt) ;	
	eta = xi*v ;
	zeta = xi*w ;

	J = Jw*pow(u,3*bt-1-al*bt) ;

	interp(x1, x2, x3, x4, xi, eta, zeta, y) ;
	wt = J*wu[i*wustr]*wv[j*wvstr]*ww[k*wwstr]/pow(A,al) ;
	func(y, x1, wt, q, nq, data) ;
      }
    }
  }
  
  return 0 ;
}

gint duffy_tet_quad(gdouble *x1, gdouble *x2,
		    gdouble *x3, gdouble *x4,
		    gdouble al,  gdouble bt,
		    gdouble *qu, gint ustr, gdouble *wu, gint wustr, gint nu,
		    gdouble *qv, gint vstr, gdouble *wv, gint wvstr, gint nv,
		    gdouble *qw, gint wstr, gdouble *ww, gint wwstr, gint nw,
		    duffy_func_t func, gpointer data,
		    gdouble *q, gint nq)

/*
 * integration on general tetrahedron with singular point at x1
 */
  
{
  gdouble V, *t1, *t2, *t3, *t4, xp[3] ;
  gdouble tol ;

  tol = 1e-12 ;

  V = orient3d(x1, x2, x3, x4) ;
  if ( fabs(V) < tol ) return 0 ;

  if ( V > 0.0 ) {
    t1 = x1 ; t2 = x2 ; t3 = x3 ; t4 = x4 ;
  } else {
    /*permute vertices to make V positive on tetrahedron*/
    t1 = x1 ; t2 = x3 ; t3 = x2 ; t4 = x4 ;
    V = -V ;
  }

  /*projection onto plane of tetrahedron base*/
  plane_project(t2, t3, t4, t1, xp) ;

#if 0
  V = orient3d(t1, xp, t2, t3) ;
  duffy_ref_tet_quad_plain(t1, xp, t2, t3, V, al, bt,
			   qu, ustr, wu, wustr, nu,
			   qv, vstr, wv, wvstr, nv,
			   qw, wstr, ww, wwstr, nw,
			   func, data, q, nq) ;
  V = orient3d(t1, xp, t3, t4) ;
  duffy_ref_tet_quad_plain(t1, xp, t3, t4, V, al, bt,
			   qu, ustr, wu, wustr, nu,
			   qv, vstr, wv, wvstr, nv,
			   qw, wstr, ww, wwstr, nw,
			   func, data, q, nq) ;
  V = orient3d(t1, xp, t4, t2) ;
  duffy_ref_tet_quad_plain(t1, xp, t4, t2, V, al, bt,
			   qu, ustr, wu, wustr, nu,
			   qv, vstr, wv, wvstr, nv,
			   qw, wstr, ww, wwstr, nw,
			   func, data, q, nq) ;
#else  
  V = orient3d(t1, xp, t2, t3) ;
  duffy_ref_tet_quad_ch(t1, xp, t2, t3, V, al, bt,
			qu, ustr, wu, wustr, nu,
			qv, vstr, wv, wvstr, nv,
			qw, wstr, ww, wwstr, nw,
			func, data, q, nq) ;
  V = orient3d(t1, xp, t3, t4) ;
  duffy_ref_tet_quad_ch(t1, xp, t3, t4, V, al, bt,
			qu, ustr, wu, wustr, nu,
			qv, vstr, wv, wvstr, nv,
			qw, wstr, ww, wwstr, nw,
			func, data, q, nq) ;
  V = orient3d(t1, xp, t4, t2) ;
  duffy_ref_tet_quad_ch(t1, xp, t4, t2, V, al, bt,
			qu, ustr, wu, wustr, nu,
			qv, vstr, wv, wvstr, nv,
			qw, wstr, ww, wwstr, nw,
			func, data, q, nq) ;
#endif
  
  return 0 ;
}

gint duffy_tet_quad_general(gdouble *x1, gdouble *x2,
			    gdouble *x3, gdouble *x4,
			    gdouble *x,
			    gdouble al,  gdouble bt,
			    gdouble *qu, gint ustr,
			    gdouble *wu, gint wustr, gint nu,
			    gdouble *qv, gint vstr,
			    gdouble *wv, gint wvstr, gint nv,
			    gdouble *qw, gint wstr,
			    gdouble *ww, gint wwstr, gint nw,
			    duffy_func_t func, gpointer data,
			    gdouble *q, gint nq)

{
  
  return 0 ;
}

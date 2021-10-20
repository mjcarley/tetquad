/* This file is part of Tetquad, a set of functions for evaluating
 * integrals on tetrahedra.
 *
 * Copyright (C) 2021 Michael Carley
 *
 * Tetquad is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. Tetquad is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tetquad.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <glib.h>

#include <gqr.h>

#include "tetquad.h"
#include "intsincos.h"
#include "duffy.h"

gchar *test_names[] = {"singular",
		       "moments",
		       ""} ;

gint tq_quadrature_js_select(gint p, gdouble **q, gint *n) ;

static gint parse_test_name(gchar *str)

{
  gint i ;

  i = 0 ; 

  while ( strlen(test_names[i]) != 0 ) {
    if ( strcmp(test_names[i], str) == 0 ) return i ;
    i ++ ;
  } ;
  
  return -1 ;
}

static gint volume_test_func(gdouble *y, gdouble *s, gdouble R, gdouble wt,
			     gdouble *q, gint nq, gpointer data[])

{
  gint N, nn, n, m, p, i ;

  N = *((gint *)data[0]) ;
  
  i = 0 ;
  for ( nn = 0 ; nn <= N ; nn ++ ) {
    for ( n = 0 ; n <= nn ; n ++ ) {
      for ( m = 0 ; m <= nn - n ; m ++ ) {
	p = nn - m - n ;

	q[i] += wt*pow(y[0], n)*pow(y[1], m)*pow(y[2], p) ;

	i ++ ;
	if ( i == nq ) return 0 ;
      }
    }
  }

  return 0 ;
}

static gint volume_test_func_d(gdouble *y, gdouble *x, gdouble wt, gdouble *q,
			       gint nq, gpointer data[])

{
  gint N, nn, n, m, p, i ;

  N = *((gint *)data[0]) ;
  i = 0 ;
  for ( nn = 0 ; nn <= N ; nn ++ ) {
    for ( n = 0 ; n <= nn ; n ++ ) {
      for ( m = 0 ; m <= nn - n ; m ++ ) {
	p = nn - m - n ;

	q[i] += wt*pow(y[0], n)*pow(y[1], m)*pow(y[2], p) ;

	i ++ ;
	if ( i == nq ) return 0 ;
      }
    }
  }

  return 0 ;
}

/* static gdouble tetvol4(gdouble *x) */

/* { */
/*   gdouble V, xq[3], x1[3], x2[3], x3[3] ; */
  
/*   tq_vector_init(x1, &(x[3*0]), &(x[3*1])) ; */
/*   tq_vector_init(x2, &(x[3*0]), &(x[3*2])) ; */
/*   tq_vector_init(x3, &(x[3*0]), &(x[3*3])) ; */
/*   tq_vector_cross(xq, x1, x2) ; */

/*   V = fabs(tq_vector_scalar(xq,x3))/6.0 ; */

/*   return V ; */
/* } */

static gint ref_quad4(gdouble *x, gdouble *qr, gint nqr,
		      tq_tetquad_func_t qfunc, gpointer qdata,
		      gdouble *q, gint nq)

{
  gint i ;
  gdouble xq[3], V, x1[3], x2[3], x3[3], R ;
  
  for ( i = 0 ; i < nq ; i ++ ) q[i] = 0.0 ;

  tq_vector_init(x1, &(x[3*0]), &(x[3*1])) ;
  tq_vector_init(x2, &(x[3*0]), &(x[3*2])) ;
  tq_vector_init(x3, &(x[3*0]), &(x[3*3])) ;

  tq_vector_cross(xq, x1, x2) ;
  V = fabs(tq_vector_scalar(xq,x3))/6.0 ;

  for ( i = 0 ; i < nqr ; i ++ ) {
    xq[0] =
      x[3*0+0]*qr[4*i+0] + x[3*1+0]*qr[4*i+1] + x[3*2+0]*qr[4*i+2] +
      x[3*3+0]*(1.0 - qr[4*i+0] - qr[4*i+1] - qr[4*i+2]) ;
    xq[1] =
      x[3*0+1]*qr[4*i+0] + x[3*1+1]*qr[4*i+1] + x[3*2+1]*qr[4*i+2] +
      x[3*3+1]*(1.0 - qr[4*i+0] - qr[4*i+1] - qr[4*i+2]) ;
    xq[2] =
      x[3*0+2]*qr[4*i+0] + x[3*1+2]*qr[4*i+1] + x[3*2+2]*qr[4*i+2] +
      x[3*3+2]*(1.0 - qr[4*i+0] - qr[4*i+1] - qr[4*i+2]) ;

    R = sqrt(tq_vector_distance2(x, xq)) ;
    
    qfunc(xq, NULL, R, V*qr[4*i+3], q, nq, qdata) ;
  }
  
  return 0 ;
}

static gint poly_test_func(gdouble *y, gdouble *s, gdouble R, gdouble wt,
			   gdouble *q, gint nq, gpointer *data)

{
  gint N = *((gint *)data[0]) ;
  gint gm = *((gint *)data[1]) ;
  gint i, n, m, p, nn ;
  
  wt *= pow(R,gm) ;
  i = 0 ;
  for ( nn = 0 ; nn <= N ; nn ++ ) {
    for ( n = 0 ; n <= nn ; n ++ ) {
      for ( m = 0 ; m <= nn - n ; m ++ ) {
	p = nn - m - n ;
	
	q[i] += pow(y[0], n)*pow(y[1], m)*pow(y[2], p)*wt ;
	i ++ ;
      }
    }
  }
  
  return 0 ;
}

static gint analytical_test(gdouble *xt, gint N, gdouble z, gint gm,
			    gint ngp, gint ngr, gint ngt, gint dmax,
			    gdouble tol, gboolean use_duffy)

{
  gqr_rule_t *qp, *qt, *qr ;
  gdouble q[1024], qa[1024], qd[1024], qq[1024], I[1024], J[1024], err[4] ;
  gdouble *work, *qrule, *qref ;
  gpointer data[4] ;
  gint nq, idx[3], n, m, p, off, i, nn, sign, d, order, nqref, Nc ;
  int_sin_cos_workspace_t *w0, *w1 ;
  tq_tetquad_func_t func = (tq_tetquad_func_t)poly_test_func ;
  duffy_func_t duffy_func = (duffy_func_t)volume_test_func_d ;
    
  fprintf(stderr, "check against analytical integration\n") ;
  fprintf(stderr, "====================================\n") ;
  fprintf(stderr, "N     = %d\n", N) ;
  fprintf(stderr, "gamma = %d\n", gm) ;
  fprintf(stderr, "z     = %lg\n", z) ;

  Nc = MAX(N, 3) ;
  
  xt[3*0+2] = z ; 

  qp = gqr_rule_alloc(ngp) ;
  qr = gqr_rule_alloc(ngr) ;
  qt = gqr_rule_alloc(ngt) ;

  order = 20 ;
  
  gqr_rule_select(qp, GQR_GAUSS_LEGENDRE, ngp, NULL) ;
  gqr_rule_select(qr, GQR_GAUSS_LEGENDRE, ngr, NULL) ;
  gqr_rule_select(qt, GQR_GAUSS_LEGENDRE, ngt, NULL) ;

  data[0] = &Nc ; data[1] = &gm ; nq = (Nc+1)*(Nc+2)*(Nc+3)/6 ;
  work = (gdouble *)g_malloc(4*nq*(dmax+1)*sizeof(gdouble)) ;

  memset(q, 0, nq*sizeof(gdouble)) ;
  tq_tet_quad(&(xt[3*0]), &(xt[3*1]), &(xt[3*2]), &(xt[3*3]),
	      &(gqr_rule_abscissa(qp,0)), 1, &(gqr_rule_weight(qp,0)), 1,
	      ngp,
	      &(gqr_rule_abscissa(qt,0)), 1, &(gqr_rule_weight(qt,0)), 1,
	      ngt,
	      &(gqr_rule_abscissa(qr,0)), 1, &(gqr_rule_weight(qr,0)), 1,
	      ngr,
	      func, data, q, nq) ;

  memset(qd, 0, nq*sizeof(gdouble)) ;
  duffy_tet_quad(&(xt[3*0]), &(xt[3*1]), &(xt[3*2]), &(xt[3*3]),
		 1.0, 1.0, 
		 &(gqr_rule_abscissa(qp,0)), 1, &(gqr_rule_weight(qp,0)), 1,
		 ngp,
		 &(gqr_rule_abscissa(qt,0)), 1, &(gqr_rule_weight(qt,0)), 1,
		 ngt,
		 &(gqr_rule_abscissa(qr,0)), 1, &(gqr_rule_weight(qr,0)), 1,
		 ngr,
		 duffy_func, data, qd, nq) ;


  memset(qa, 0, nq*sizeof(gdouble)) ;
  d = tq_tet_quad_adaptive(&(xt[3*0]), &(xt[3*1]), &(xt[3*2]), &(xt[3*3]),
			   &(gqr_rule_abscissa(qp,0)), 1,
			   &(gqr_rule_weight(qp,0)), 1,
			   ngp,
			   &(gqr_rule_abscissa(qt,0)), 1,
			   &(gqr_rule_weight(qt,0)), 1,
			   ngt,
			   &(gqr_rule_abscissa(qr,0)), 1,
			   &(gqr_rule_weight(qr,0)), 1,
			   ngr,
			   func, data,
			   tol, dmax, qa, nq, work) ;

  memset(qq, 0, nq*sizeof(gdouble)) ;
  tq_quadrature_js_select(order, &qrule, &nqref) ;
  ref_quad4(xt, qrule, nqref, func, data, qq, nq) ;

  w0 = int_sin_cos_workspace_alloc(4*Nc) ;
  w1 = int_sin_cos_workspace_alloc(4*Nc) ;
  idx[0] = 1 ; idx[1] = 2 ; idx[2] = 3 ;

  sign = 1 ;
  pyramid_quad_mn(xt, idx, 3, 3, z, gm, Nc, I, work, w0, w1) ;
  /*check and correct for negative volume case*/
  if ( I[0] < 0.0 ) sign = -1 ;

  if ( use_duffy ) qref = qd ; else qref = J ;
  
  err[0] = err[1] = err[2] = err[3] = 0.0 ;
  off = (Nc+1)*(Nc+2)/2 + 128 ; i = 0 ;
  for ( nn = 0 ; nn <= N ; nn ++ ) {
    for ( n = 0 ; n <= nn ; n ++ ) {
      for ( m = 0 ; m <= nn - n ; m ++ ) {
	p = nn - m - n ;
	
	J[i] = sign*I[(n+m)*(n+m+1)/2+n]*I[off + p*(Nc+gm+2+1) + (n+m)+gm+2] ;
	fprintf(stdout, "%d %d %d %d %lg (%lg, %lg, %lg, %lg)\n",
		i, n, m, p, qref[i],
		fabs(qref[i] - q[i]), fabs(qref[i] - qa[i]),
		fabs(qref[i] - qd[i]), fabs(qref[i] - qq[i])) ;
	err[0] = MAX(err[0], fabs(qref[i]-q [i])) ;
	err[1] = MAX(err[1], fabs(qref[i]-qa[i])) ;
	err[2] = MAX(err[2], fabs(qref[i]-qd[i])) ;
	err[3] = MAX(err[3], fabs(qref[i]-qq[i])) ;

	i ++ ;
      }
    }
  }

  /* for ( i = 0 ; i < 4 ; i ++ ) err[i] /= J[0] ; */
  
  fprintf(stderr, "errors: %lg %lg %lg %lg\n",
	  log10(err[0]), log10(err[1]), log10(err[2]), log10(err[3])) ;
  fprintf(stderr, "recursion depth: %d\n", d) ;
  
  return 0 ;
}

static gint moment_test(gdouble *xt, gint N, gint order,
			gint ngp, gint ngr, gint ngt)

{
  gqr_rule_t *qp, *qt, *qr ;
  gdouble q[1024], J[1024], *qrule, err ;
  gint nq, n, m, p, i, j, nn, nqref ;
  gint rot[] = {0, 1, 2, 3, 0, 1, 2, 3} ;
  gpointer data[4] ;
  
  fprintf(stderr, "monomial integration with rotation of tet\n") ;
  fprintf(stderr, "=========================================\n") ;
  fprintf(stderr, "N     = %d\n", N) ;

  qp = gqr_rule_alloc(ngp) ;
  qr = gqr_rule_alloc(ngr) ;
  qt = gqr_rule_alloc(ngt) ;

  gqr_rule_select(qp, GQR_GAUSS_LEGENDRE, ngp, NULL) ;
  gqr_rule_select(qr, GQR_GAUSS_LEGENDRE, ngr, NULL) ;
  gqr_rule_select(qt, GQR_GAUSS_LEGENDRE, ngt, NULL) ;

  nq = (N+1)*(N+2)*(N+3)/6 ;

  data[0] = &N ;
  
  tq_quadrature_js_select(order, &qrule, &nqref) ;
  ref_quad4(xt, qrule, nqref, volume_test_func, data, J, nq) ;
  
  for ( j = 0 ; j < 4 ; j ++ ) {
    memset(q, 0, nq*sizeof(gdouble)) ;
    tq_tet_quad(&(xt[3*rot[j+0]]), &(xt[3*rot[j+1]]),
		&(xt[3*rot[j+2]]), &(xt[3*rot[j+3]]),
		&(gqr_rule_abscissa(qp,0)), 1, &(gqr_rule_weight(qp,0)), 1,
		ngp,
		&(gqr_rule_abscissa(qt,0)), 1, &(gqr_rule_weight(qt,0)), 1,
		ngt,
		&(gqr_rule_abscissa(qr,0)), 1, &(gqr_rule_weight(qr,0)), 1,
		ngr,
		(tq_tetquad_func_t)volume_test_func, data, q, nq) ;

    i = 0 ; err = 0.0 ;
    for ( nn = 0 ; nn <= N ; nn ++ ) {
      for ( n = 0 ; n <= nn ; n ++ ) {
	for ( m = 0 ; m <= nn - n ; m ++ ) {
	  p = nn - m - n ;
	  
	  fprintf(stdout, "%d %d %d %d %d %lg %lg (%lg)\n",
		  j, i, n, m, p, J[i], q[i], fabs(J[i] - q[i])) ;
	  err = MAX(err, fabs(J[i]-q[i])) ;
	  i ++ ;
	}
      }
    }
    fprintf(stderr, "rotation %d, error: %lg\n", j, err) ;
  }
  
  return 0 ;
}

static void print_help_message(FILE *f, gchar *name)

{
  fprintf(f,
	  "Usage:\n"
	  "  %s [options] -i [input]\n\n",
	  name) ;

  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -D use Duffy transformation\n"
	  "  -d # maximum recursion depth in adaptive quadrature\n"
	  "  -e # error tolerance for adaptive quadrature\n"
	  "  -i # input file name\n"
	  "  -N # maximum order of monomials in test integrands\n"
	  "  -p # number of nodes in quadrature rule for phi\n"
	  "  -r # number of nodes in quadrature rule for rho\n"
	  "  -t # number of nodes in quadrature rule for theta\n"
	  "  -T # test case to run\n"
	  "  -z # tetrahedron height\n") ;
    
    return ;
}

gint main(gint argc, gchar **argv)

{
  gdouble x[12], z, tol ;
  gint ngp, ngt, ngr, order, N, test, i, dmax ;
  FILE *input ;
  /* , *output ; */
  gchar ch, *progname, *ipfile ;
  gboolean use_duffy ;
  
  input = stdin ;
  /* output = stdout ; */
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  use_duffy = FALSE ;
  z = 1.0 ; N = 4 ;

  dmax = 4 ; tol = 1e-6 ;
  
  ngp = 4 ; ngr = 4 ; ngt = 4 ;
  order = 10 ;
  ipfile = NULL ;
  
  test = -1 ;
  while ( (ch = getopt(argc, argv, "hDd:e:i:N:p:r:T:t:z:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h': print_help_message(stderr, progname) ; return 0 ; break ;
    case 'D': use_duffy = TRUE ; break ;
    case 'd': dmax = atoi(optarg) ; break ;
    case 'e': tol = atof(optarg) ; break ;
    case 'i': ipfile = g_strdup(optarg) ; break ;
    case 'N': N = atoi(optarg) ; break ;
    case 'p': ngp = atoi(optarg) ; break ;
    case 'r': ngr = atoi(optarg) ; break ;
    case 't': ngt = atoi(optarg) ; break ;
    case 'T': test = parse_test_name(optarg) ; break ;
    case 'z': z = atof(optarg) ; break ;
    }
  }

  if ( test == -1 ) {
    fprintf(stderr, "%s: unrecognized or undefined test\n", progname) ;

    return 1 ;
  }

  x[3*0+0] = 0.0 ; x[3*0+1] = 0.0 ; x[3*0+2] = z ; 
  x[3*1+0] = 0.0 ; x[3*1+1] = 0.0 ; x[3*1+2] = 0.0 ;
  x[3*2+0] = 1.0 ; x[3*2+1] = 0.0 ; x[3*2+2] = 0.0 ;
  x[3*3+0] = 0.0 ; x[3*3+1] = 1.0 ; x[3*3+2] = 0.0 ;

  if ( ipfile != NULL ) {
    input = fopen(ipfile, "r") ;
    if ( input == NULL ) {
      fprintf(stderr, "%s: cannot open input file %s\n",
	      progname, ipfile) ;
      return 1 ;
    }

    for ( i = 0 ; i < 12 ; i ++ ) fscanf(input, "%lg", &(x[i])) ;

    fclose(input) ;
  }

  if ( test == 0 ) {
    analytical_test(x, N, z, -1, ngp, ngr, ngt, dmax, tol, use_duffy) ;

    return 0 ;
  }

  if ( test == 1 ) {
    moment_test(x, N, order, ngp, ngr, ngt) ;

    return 0 ;
  }

  return 0 ;
}

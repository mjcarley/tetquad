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
#include <string.h>
#include <math.h>

#include <glib.h>

#include <blaswrap.h>

#include "tetquad.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#endif /*HAVE_AVX_INSTRUCTIONS*/

#define TQ_ADATA_SIZE  32
#define TQ_ADATA_QTH    0
#define TQ_ADATA_QTSTR  1
#define TQ_ADATA_WTH    2
#define TQ_ADATA_NTH    3
#define TQ_ADATA_WTSTR  4
#define TQ_ADATA_QPH    5
#define TQ_ADATA_QPSTR  6
#define TQ_ADATA_WPH    7
#define TQ_ADATA_WPSTR  8
#define TQ_ADATA_NPH    9
#define TQ_ADATA_QR    10
#define TQ_ADATA_QRSTR 11
#define TQ_ADATA_WR    12
#define TQ_ADATA_WRSTR 13
#define TQ_ADATA_NR    14
#define TQ_ADATA_QFUNC 15
#define TQ_ADATA_QDATA 16

static void adaptive_quad_eval(gdouble *x1, gdouble *x2,
			       gdouble *x3, gdouble *x4,
			       gpointer adata[], gdouble *work, gint nq)

{
  gdouble *qth = adata[TQ_ADATA_QTH  ] ;
  gint qtstr = *((gint *)adata[TQ_ADATA_QTSTR]) ;
  gdouble *wth = adata[TQ_ADATA_WTH  ] ; 
  gint wtstr = *((gint *)adata[TQ_ADATA_WTSTR]) ; 
  gint nth = *((gint *)adata[TQ_ADATA_NTH]) ; 
  gdouble *qph = adata[TQ_ADATA_QPH  ] ; 
  gint qpstr = *((gint *)adata[TQ_ADATA_QPSTR]) ; 
  gdouble *wph = adata[TQ_ADATA_WPH  ] ; 
  gint wpstr = *((gint *)adata[TQ_ADATA_WPSTR]) ; 
  gint nph = *((gint *)adata[TQ_ADATA_NPH]) ; 
  gdouble *qr =	adata[TQ_ADATA_QR   ] ; 
  gint qrstr = *((gint *)adata[TQ_ADATA_QRSTR]) ; 
  gdouble *wr =	adata[TQ_ADATA_WR   ] ; 
  gint wrstr = *((gint *)adata[TQ_ADATA_WRSTR]) ; 
  gint nr = *((gint *)adata[TQ_ADATA_NR]) ; 
  tq_tetquad_func_t qfunc = adata[TQ_ADATA_QFUNC] ; 
  gpointer qdata = adata[TQ_ADATA_QDATA] ; 

  tq_tet_quad(x1, x2, x3, x4,
	      qph, qpstr, wph, wpstr, nph,
	      qth, qtstr, wth, wtstr, nth,
	      qr , qrstr, wr , wrstr, nr ,
	      qfunc, qdata, work, nq) ;
  
  return ;
}

static void adaptive_quad_recursion(gdouble *x1, gdouble *x2,
				    gdouble *x3, gdouble *x4,
				    gpointer adata[], gdouble *work, gint nq,
				    gint tri,
				    gint d, gint dmax, gdouble tol)

{
  gdouble x23[3], x34[3], x42[3] ;
  gint off0, off1, i ;
  gboolean recurse ;
  
  if ( d == dmax ) return ;
  /*split the face and calculate the integral on the resulting tetrahedra*/
  tq_vector_interp(x23, x2, x3, 0.5) ;
  tq_vector_interp(x34, x3, x4, 0.5) ;
  tq_vector_interp(x42, x4, x2, 0.5) ;

  /*workspace offsets for this depth of recursion*/
  off0 = 4*nq*(d-1)+tri*nq ; off1 = 4*nq*d ;

  memset(&(work[off1]), 0, 4*nq*sizeof(gdouble)) ;
  adaptive_quad_eval(x1, x2 , x23, x42, adata, &(work[off1+0*nq]), nq) ;
  adaptive_quad_eval(x1, x23,  x3, x34, adata, &(work[off1+1*nq]), nq) ;
  adaptive_quad_eval(x1, x34,  x4, x42, adata, &(work[off1+2*nq]), nq) ;
  adaptive_quad_eval(x1, x23, x34, x42, adata, &(work[off1+3*nq]), nq) ;

  recurse = FALSE ;
  for ( i = 0 ; i < nq ; i ++ ) {
    if ( fabs(work[off0+i] -
	      work[off1+0*nq+i] - work[off1+1*nq+i] -
	      work[off1+2*nq+i] - work[off1+3*nq+i]) > tol ) {
      recurse = TRUE ; break ;
    }	      
  }

  if ( !recurse ) return ;
  /* tol *= 0.25 ; */

  adaptive_quad_recursion(x1,  x2, x23, x42, adata, work, nq, 0, d+1, dmax,
			  tol) ;
  adaptive_quad_recursion(x1, x23,  x3, x34, adata, work, nq, 1, d+1, dmax,
			  tol) ;
  adaptive_quad_recursion(x1, x34,  x4, x42, adata, work, nq, 2, d+1, dmax,
			  tol) ;
  adaptive_quad_recursion(x1, x23, x34, x42, adata, work, nq, 3, d+1, dmax,
			  tol) ;

  for ( i = 0 ; i < nq ; i ++ ) {
    work[off0+i] = work[off1+0*nq+i] + work[off1+1*nq+i] +
      work[off1+2*nq+i] + work[off1+3*nq+i] ;
  }
  
  return ;
}

gint tq_tet_quad_adaptive(gdouble *x1, gdouble *x2, gdouble *x3, gdouble *x4,
			  gdouble *qph, gint qpstr, gdouble *wph, gint wpstr,
			  gint nph,
			  gdouble *qth, gint qtstr, gdouble *wth, gint wtstr,
			  gint nth,
			  gdouble *qr, gint qrstr, gdouble *wr, gint wrstr,
			  gint nr,
			  tq_tetquad_func_t qfunc, gpointer qdata,
			  gdouble tol, gint dmax,
			  gdouble *q, gint nq, gdouble *work)

{
  gpointer adata[TQ_ADATA_SIZE] ;
  gint i ;
  
  adata[TQ_ADATA_QTH  ] = qth ;
  adata[TQ_ADATA_QTSTR] = &qtstr ;
  adata[TQ_ADATA_WTH  ] = wth ;
  adata[TQ_ADATA_WTSTR] = &wtstr ;
  adata[TQ_ADATA_NTH  ] = &nth ;
  adata[TQ_ADATA_QPH  ] = qph ;
  adata[TQ_ADATA_QPSTR] = &qpstr ;
  adata[TQ_ADATA_WPH  ] = wph ;
  adata[TQ_ADATA_WPSTR] = &wpstr ;
  adata[TQ_ADATA_NPH  ] = &nph ;
  adata[TQ_ADATA_QR   ] = qr ;
  adata[TQ_ADATA_QRSTR] = &qrstr ;
  adata[TQ_ADATA_WR   ] = wr ;
  adata[TQ_ADATA_WRSTR] = &wrstr ;
  adata[TQ_ADATA_NR  ] = &nr ;
  adata[TQ_ADATA_QFUNC] = qfunc ;
  adata[TQ_ADATA_QDATA] = qdata ;
  
  memset(work, 0, 4*nq*dmax*sizeof(gdouble)) ;

  adaptive_quad_eval(x1, x2, x3, x4, adata, work, nq) ;

  adaptive_quad_recursion(x1, x2, x3, x4, adata, work, nq, 0, 1, dmax, tol) ;
  
  for ( i = 0 ; i < nq ; i ++ ) q[i] += work[i] ;
  
  return 0 ;
}

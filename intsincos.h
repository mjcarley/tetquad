#ifndef _INTSINCOS_H_INCLUDED_
#define _INTSINCOS_H_INCLUDED_

#define INTSINCOS_ORDER 12

double orient2d(double *pa, double *pb, double *pc) ;

#define _table_offset(_n) ( ((_n) == -1 ? 3 : (9 + (_n)*((_n)+1)/2)))
#define _table_index(_m,_n) (((_n)==-1) && ((_m)==-1) ? 0 :		\
			     ((_n) == -1 ? ((_m)+1) :			\
			      (((_m) == -1 ? ((_n)+5) :			\
				_table_offset(((_m)+(_n)))+(_n)))))

typedef struct {
  gint N, Nmax ;
  gdouble t, S, C, d, k, k2, kd, kd2 ;
  gdouble *Im1, *I0, *Ip1 ;
} int_sin_cos_workspace_t ;

int_sin_cos_workspace_t *int_sin_cos_workspace_alloc(gint Nmax) ;
gint int_sin_cos_workspace_init(int_sin_cos_workspace_t *w,
				gdouble t, gdouble k, gint N) ;
gdouble int_sin_cos_delta_mn(int_sin_cos_workspace_t *w,
			     gint m, gint n, gint p) ;
gint int_sin_cos_delta_expand_mn(gint M, gint N, gint r, gdouble t,
				 gdouble S, gdouble C, gdouble d,
				 gdouble k, gdouble k2, 
				 gdouble kd, gdouble kd2,
				 gdouble *Imn) ;
gint subtriangle_quad_mn(gdouble r1, gdouble r2, gdouble th,
			 gdouble psi, gdouble z,
			 gint N, gint gm, gdouble *I,
			 int_sin_cos_workspace_t *w0,
			 int_sin_cos_workspace_t *w1) ;
gint poly_quad_mn(gdouble *x, gint *idx, gint str, 
		  gint nx, gdouble z, gint gm, gint N,
		  gdouble *I, gdouble *work,
		  int_sin_cos_workspace_t *w0,
		  int_sin_cos_workspace_t *w1) ;
gint pyramid_quad_mn(gdouble *x, gint *idx, gint str, 
		     gint nx, gdouble z, gint gm, gint N,
		     gdouble *I, gdouble *work,
		     int_sin_cos_workspace_t *w0,
		     int_sin_cos_workspace_t *w1) ;
gdouble pyramid_quad_general(gint gm, gint N, gdouble *I,
			     gdouble o[],
			     gdouble s[], gdouble t[], gdouble n[],
			     gint mx, gint my, gint mz, 
			     gdouble z, gint grad) ;
gdouble polygon_quad_general(gint gm, gint N, gdouble *I,
			     gdouble o[],
			     gdouble s[], gdouble t[], gdouble n[],
			     gint mx, gint my, gint mz) ;
gint tet_projection(gdouble *x0, gdouble *x1, gdouble *x2, gdouble *x3,
		    gdouble *o, gdouble *s, gdouble *t, gdouble *n,
		    gdouble *y1, gdouble *y2, gdouble *y3,
		    gdouble *z) ;

#endif /*_INTSINCOS_H_INCLUDED_*/

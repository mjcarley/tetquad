Tetquad is basic (and currently quite crude) code for integrating
singular functions on tetrahedra. The implementation here is
developmental code which contains a lot of material which will be
removed or rewritten. To compile the test program, you will also need
the gqr library for Gaussian Quadrature Rules:

https://github.com/mjcarley/gqr

The test code will evaluate integrals on tetrahedra read from file or
on an internally defined geometry, and give an estimate of the error
compared to an analytical calculation. 

The algorithm is described in

https://arxiv.org/abs/2110.07270

The method is implemented in the files tetquad.c and tetquad.h, which
can stand alone in another code. The other source files are included
for the calculation of reference integrals using an analytical method,
or a Duffy transformation. In order to use the quadrature function in
your own code, you will need to supply precomputed Gauss-Legendre
quadrature rules and a function for the integrand as a
tq_tetquad_func_t, defined in tetquad.h

typedef gint (*tq_tetquad_func_t)(gdouble *y, gdouble *s,
				  gdouble R, gdouble wt,
				  gdouble *q, gint nq, gpointer data) ;

An integrand function of this type should interpret the arguments as follows:

y:    evaluation point of integrand;
s:    unit magnitude vector from singular point to y;
R:    distance from singular point to y;
wt:   quadrature weight;
q:    integral;
nq:   number of terms to be integrated;
data: user data to pass to function;

and should increment the nq entries of q as follows

q[i] += wt*f_i(y), 0<= i < nq

In the function poly_test_func in test.c, this is demonstrated for the
integration of monomials divided by R.

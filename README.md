Tetquad is basic (and currently quite crude) code for integrating
singular functions on tetrahedra. The implementation here is
developmental code which contains a lot of material which will be
removed or rewritten. To compile the test program, you will also need
the gqr library for Gaussian Quadrature Rules:

https://github.com/mjcarley/gqr

The test code will evaluate integrals on tetrahedra read from file or
on an internally defined geometry, and give an estimate of the error
compared to an analytical calculation. 

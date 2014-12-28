
 ===== readme for fmm =====

 Routines from the book Computer Methods for Mathematical
 Computations, by Forsythe, Malcolm, and Moler.
 Developed by George Forsythe, Mike Malcolm, and Cleve Moler.




Downloaded from http://www.netlib.no/netlib/fmm/
• file: fmm/readme
	o for: overview of fmm 
• file: fmm/decomp.f
	o for: decomposes a matrix by Gaussian elimination and estimates the , condition of the matrix 
	o prec: double 
• file: fmm/solve.f
	o for: solution of linear system, A*x = b, do not use if (fmm/decomp) has , detected a singularity 
	o prec: double 
• file: fmm/quanc8.f
	o for: estimate the integral of f(x) in a finite interval, user provided , tolerance using an automatic adaptive routine based on the 8-panel , Newton-Cotes rule 
	o prec: double 
• file: fmm/rkf45.f
	o for: Fehlberg fourth-fifth order Runge-Kutta method 
	o prec: double 
• file: fmm/spline.f
	o for: compute the coefficients for a cubic interpolating spline 
	o prec: double 
• file: fmm/seval.f
	o for: evaluate a cubic interpolating spline 
	o prec: double 
• file: fmm/svd.f
	o for: determines the singular value decomposition, SVD, of a real , rectangular matrix, using Householder bidiagonalization and a variant , of the QR algorithm 
	o prec: double 
• file: fmm/fmin.f
	o for: an approximation to the point where a user function attains a minimum , on an interval is determined 
	o prec: double 
• file: fmm/urand.f
	o for: is a uniform random number generator based on theory and suggestions , given in D.E. Knuth (1969), Vol 2 
	o prec: double 
• file: fmm/zeroin.f
	o for: find a zero of a user function in an interval 
	o prec: double 

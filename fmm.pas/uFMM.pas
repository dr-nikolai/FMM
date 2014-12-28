{@abstract(Computer Methods for Mathematical Computation.
           Pascal translation of the famous Forsythe-Malcolm-Moler collection
           of Fortran procedures)
@author(Nikolai Shokhirev <nikolai@u.arizona.edu>
                          http://www.shokhirev.com/nikolai.html)
@created(2004.12.12)
@lastmod(2004.12.24)
©Nikolai V. Shokhirev, 2001-2004

  REFERENCES
  1. "Computer Methods for Mathematical Computations", by George E. Forsythe,
      Michael A. Malcolm, and Cleve B. Moler. Prentice-Hall, 1977.
  2. Netlib Repository http://www.netlib.org/fmm/
  3. Public Domain Aeronautical Software (PDAS) http://www.pdas.com/free.htm
  4. ForToPas - FORTRAN to Pascal converter, Nikolai Shokhirev,
     http://www.shokhirev.com/nikolai/programs/tools.html

From PDAS comments:

  PURPOSE - A collection of Fortran procedures for mathematical computation
    based on the procedures from the book [1].

  COLLECTION AUTHORS - George E. Forsythe, Michael A. Malcolm, &
    Cleve B. Moler, Stanford University (1977)
    Ralph L. Carmichael, Public Domain Aeronautical Software

  ORIGINAL AUTHORS -
    Decomp,Solve - Forsythe,Malcolm,Moler
    Spline,Seval - Forsythe,Malcolm,Moler
    Quanc8 - Forsythe,Malcolm,Moler
    RKF45 - H.A. Watts AND L.F. Shampine  (Sandia)
    ZEROIN - Richard Brent
    FMIN - Van Winjngaarden,Decker,Brent, Et Al
    SVD - Golub and Reinsch                       }
unit uFMM;

interface

uses
  uMatTypes,
  uDynArrays;

type
  TDiffEqProc = procedure(const t: TFloat; const y, yp: IFArr1D);

{
  Pascal version by Nikolai V. Shokhirev
PURPOSE:
  Solution of linear system, a*x = b - translation from FORTRAN            

  *** Do not use if FMMDecomp has detected singularity. ***
input:                                                                     
  a[1..n,1..n] - triangularized matrix obtained from FMMDecomp.            
  b[1..n] - right hand side vector.                                        
  ipvt[1..n] - pivot vector obtained from FMMDecomp .                      
output:                                                                    
  b[1..n] - solution vector, x .                                         }
procedure FMMSolve(const a: IFArr2D; const b: IFArr1D; const ipvt: IIArr1D);


{
  Pascal version by Nikolai V. Shokhirev
PURPOSE:
  Decomposes a matrix by Gaussian elimination and estimates the condition  
  of the matrix. - translation from FORTRAN                                

  *** Use FMMSolve to compute solutions to linear systems. ***
input:                                                                     
  a[1..n,1..n] - matrix to be triangularized.                              
output:                                                                    
  a - contains an upper triangular matrix u and a permuted version of a lower
      triangular matrix i-l so that (permutation matrix)*a = l*u .         
  cond - an estimate of the condition of a .                               
      For the linear system a*x = b, changes in a and b may cause changes cond
      times as large in x. If cond+1.0 = cond , a is singular to working   
      precision. cond is set to 1.0e+32 if exact singularity is detected.  
  ipvt - the pivot vector.                                                 
      ipvt[k] = the index of the k-th pivot row                            
      ipvt[n] = (-1)**(number of interchanges)                             
  work space..                                                             
     the vector work must be declared and included in the call. Its input  
     contents are ignored. Its output contents are usually unimportant.    
   The determinant of a can be obtained on output by                       
     det(a) = ipvt[n] * a[1,1] * a[2,2] * ... * a[n,n].                  }
procedure FMMDecomp(const a: IFArr2D; var cond: TFloat; const ipvt: IIArr1D);

{
  Pascal version by Nikolai V. Shokhirev
PURPOSE:
  Compute the coefficients b,c,d for a cubic interpolating spline so that the
  interpolated value is given by
    s[x] = y[k] + b[k]*(x-x[k]) + c[k]*(x-x[k])**2 + d[k]*(x-x[k])**3
      when x[k] <= x <= x[k+1]
  The end conditions match the third derivatives of the interpolated curve to
  the third derivatives of the unique polynomials through the first four and
  last four points.

  Use Seval or Seval3 to evaluate the spline.
INPUT:
   x  -  abscissas of knots
   y  -  ordinates of knots
OUTPUT:
   b  -  linear coeff
   c  -  quadratic coeff.
   d  -  cubic coeff.
  All arrays have the same limits       }
procedure FMMSpline(const x, y, b, c, d: IFArr1D);

{
  Pascal version by Nikolai V. Shokhirev                                   
PURPOSE:
  Construct the natural spline through a set of points
NOTES:
  A natural spline has zero second derivative at both endpoints.
INPUT:
  x, y - coordinates of knots
OUTPUT:
   b, c, d  -  linear, quadratic and cubic coefficients          }
procedure NaturalSpline(const x, y, b, c, d: IFArr1D);

{
  Pascal version by Nikolai V. Shokhirev                                   
PURPOSE:
  Evaluate the cubic spline function
     Seval=y[i]+b[i]!(u-x[i])+c[i]*(u-x[i])**2+d[i]*(u-x[i])**3
           where  x[i] <= u < x[i+1]
NOTES:
  if u<x[Lo], i = Lo is used;if u>x[Hi], i = Hi is used
INPUT:
  u - abscissa at which the spline is to be evaluated
  x - abscissas of knots
  y - ordinates of knots
  b, c, d - linear, quadratic, cubic coeff
IN-OUT
  i - last used index; should be set to x.Lo1 before the 1st call  }
function Seval(const u: TFloat; var i: TInt; const x, y, b, c, d: IFArr1D): TFloat;

{
  Pascal version by Nikolai V. Shokhirev                                   
PURPOSE:
  Evaluate the cubic spline function

     Seval=y[i]+b[i]!(u-x[i])+c[i]*(u-x[i])**2+d[i]*(u-x[i])**3
           where  x[i] <= u < x[i+1]
NOTES:
  if u<x[Lo], i=1 is used;if u>x[Hi], i=Dim is used

  u - abscissa at which the spline is to be evaluated
  x - abscissas of knots
  y - ordinates of knots
  b, c, d - linear,quadratic,cubic coeff
  f, fp, fpp, fppp - function, 1st, 2nd, 3rd derivatives
IN-OUT
  i - last used index; should be set to x.Lo1 before the 1st call }
procedure Seval3(const u: TFloat; var i: TInt; const x, y, b, c, d: IFArr1D;
                       var f, fp, fpp, fppp: TFloat);

{
  Pascal version by Nikolai V. Shokhirev                                   
PURPOSE:
  Estimate the integral of FUN(X) from A to B to a user provided tolerance.
  This is an automatic adaptive routine based on the 8-panel Newton-Cotes rule.
INPUT:
  a - lower limit of integration
  b - upper limit of integration
  abserr - absolute error tolerance (>0)
  relerr - relative error tolerance (>0)
OUTPUT:
  result - approx. value of integral
  errest - estimate of actual error
  nofun  - no. of function evaluations
  flag   - reliability indicator. ( = 0 is O.K.)   }
procedure Quanc8(Fun: TFloatFunction1D; const a, b, abserr, relerr: TFloat;
                            var result, errest, flag: TFloat; var nofun: TInt);

{
  Pascal version by Nikolai V. Shokhirev                                   
PURPOSE:
  Integrate a system of neqn first order differential equations by the Fehlberg
  fourth-fifth order Runge-Kutta method.
NOTES:
  RKF45 is primarily designed to solve non-stiff and mildly stiff differential
  equations when derivative evaluations are inexpensive.
  Rkf45 should generally not be used when the user is demanding high accuracy.

  The arrays yp, f1, f2, f3, f4, and f5 (of dimension at least neqn) and the
  variables h, savre, savae, nfe, kop, init,jflag,and kflag are used internally
  by the code and appear in the call list to eliminate local retention of
  variables between calls. Accordingly, they should not be altered. Items of
  possible interest are
    yp - derivative of solution vector at t
    h - an appropriate stepsize to be used for the next step
    nfe- counter on the number of derivative function evaluations

ABSTRACT:
  Subroutine Rkf45 integrates a system of neqn first order ordinary differential
  equations of the form

 	  dy[i]/dt = f(t,y[1],y[2],...,y[neqn])

 	where the y[i] are given at t .
  Typically the subroutine is used to integrate from t to tout but it can be
  used as a one-step integrator to advance the solution a Double step in the
  direction of tout. On return, the parameters in the call list are set for
  continuing the integration. The user has only to call Rkf45 again (and perhaps
  define a new value for tout). Rkfs45 calls subroutine Fehl which computes an
  approximate solution over one step.

  Rkf45 uses the Runge-Kutta-Fehlberg(4,5) method described in the reference
  E. Fehlberg:  Low-order Classical Runge-Kutta formulas with Stepsize 
  Control, NASA TR R-315.

  The performance of Rkf45 is illustrated in the reference
  L.F.Shampine,H.A.Watts,S.Davenport: "Solving Non-stiff Ordinary
  Differential Equations - The State of the Art.
   Sandia Laboratories Report SAND75-0182, to appear in SIAM Review. 

  First call to Rkf45:

  The user must provide storage in his calling program for the arrays in
  the call list y, yp, f1, f2, f3, f4, f5 - 1D float arrays with the same
  limits and dimension = neqn, supply procedure F(t,y,yp) and initialize the
  following parameters:
    y - vector of initial conditions
    t - starting point of integration , must be a variable
    tout - output point at which solution is desired.
    t = tout is allowed on the first call only, in which case Rkf45 returns with
      iflag=2 if continuation is possible.
    relerr, abserr - relative and absolute local error tolerances which must be
      non-negative. relerr must be a variable while abserr may be a constant.
      The code should normally not be used with relative error control smaller
      than about i.e-8 to avoid limiting precision difficulties the code
      requires relerr to be larger than an internally computed relative error
      parameter which is machine dependent. In particular, pure absolute error
      is not permitted. If a smaller than allowable value of relerr is
      attempted, rkf45 increases relerr appropriately and returns control to the
      user before  continuing the integration.
    iflag - +1,-1 indicator to initialize the code for each new problem.
      Normal input is +1. the user should set iflag = -1 only when one-step
      integrator control is essential. In this case, Rkf45 attempts to advance
      the solution a Double step in the direction of tout each time it is
      called. Since this mode of operation results in extra computing overhead,
      it should be avoided unless needed.

  Output from Rkf45

  y - solution at t
  t - last point reached in integration.
  iflag =  2 - integration reached tout. Indicates successful return and is the
               normal mode for continuing integration.
       	= -2 - a Double successful step in the direction of tout has been taken.
               Normal mode for continuing integration one step at a time.
 	      =  3 - integration was not completed because relative error tolerance
               was too small. relerr has been increased	appropriately for
               continuing.
 	      =  4 - integration was not completed because more than 3000 derivative
               evaluations were needed. This is approximately 500 steps.
 	      =  5 - integration was not completed because solution vanished making
               a pure relative error test impossible. Must use non-zero abserr
               to continue. Using the one-step integration mode for one step
               is a good way to proceed.
 	      =  6 - integration was not completed because requested accuracy could
               not be achieved using smallest allowable stepsize. User must
               increase the error tolerance before continued integration can
               be attempted.
 	      =  7 - it is likely that Rkf45 is inefficient for solving this problem.
               Too much output is restricting the natural stepsize choice.
               Use the one-step integrator mode.
 	      = 8  - invalid input parameters. This indicator occurs if any of the
               following is satisfied:
                 t = tout and iflag <> +1 or -1
 	               relerr or abserr < 0.
                 iflag = 0 or < -2 or > 8
  h, nfe, yp, f1, f2, f3, f4, f5 - information which is usually of no interest
  to the user but necessary for subsequent calls.
  - yp contain the first derivatives of the solution vector y at t.
  - h contains the stepsize to be attempted on the next step.
  - nfe contains the derivative evaluation counter.

  Subsequent calls to Rkf45

  Subroutine Rkf45 returns with all information needed to continue the
  integration. If the integration reached tout, the user need only define
  a new tout and call Rkf45 again. In the one-step integrator mode (iflag=-2)
  the user must keep in mind that each step taken is in the direction of the
  current tout. Upon reaching tout (indicated by changing iflag to 2),the 
  user must then define a new tout and reset iflag to -2 to continue in the 
  one-step integrator mode. 

  If the integration was not completed but the user still wants to continue 
  (iflag = 3,4 cases), he just calls Rkf45 again. With iflag = 3, the relerr
  parameter has been adjusted appropriately for continuing the integration.
  In the case of iflag=4 the function counter will be reset to 0 and 
  another 3000 function evaluations are allowed. 

  However,in the case iflag = 5, the user must first alter the error criterion
  to use a positive value of abserr before integration can proceed. If he 
  does not, execution is terminated. 

  Also, in the case iflag = 6, it is necessary for the user to reset iflag
  to 2 (or -2 when the one-step integration mode is being used) as well as 
  increasing either abserr,relerr or both before the integration can be 
  continued. If this is not done, execution will be terminated. The 
  occurrence of iflag=6 indicates a trouble spot (solution is changing 
  rapidly,singularity may be present) and it often is inadvisable to continue.

  If (flag = 7 is encountered, the user should use the one-step integration mode
  with the stepsize determined by the code or consider switching to the 
  Adams codes DE/STEP, INTRP. If the user insists upon continuing the
  integration with Rkf45, he must reset iflag to 2 before calling Rkf45 again.
  Otherwise, execution will be terminated.

  If iflag = 8 is obtained, integration can not be continued unless the invalid
  input parameters are corrected. 

  It should be noted that the arrays work,iwork contain information required
  for subsequent integration. Accordingly, work and iwork should not be 
  altered.

  This interfacing routine merely relieves the user of a long calling list 
  via the splitting apart of two working storage arrays.

INPUT
  neqn - number of equations
  y - solution vector at t
  t - independent variable
  tout - output point at which solution is desired
  relerr - relative error tolerance
  abserr -  absolute error tolerance
OUTPUT
  iflag - indicator for status of work
  h - step size
  savre;
  savae;
  nfe - number of function evaluations
  kop;
IN, OUT
  yp - derivatives
  f1, f2, f3, f4, f5
  init;
  jflag;
  kflag;                               }
procedure Rkfs45(F: TDiffEqProc; var t: TFloat; const tout, abserr: TFloat;
                 var relerr: TFloat; var iflag: TInt;
                 const y, yp, f1, f2, f3, f4, f5: IFArr1D;
                 var h, savre, savae: TFloat;
                 var nfe, kop, init, jflag, kflag: TInt);

{
  Pascal version by Nikolai V. Shokhirev                                   
PURPOSE:
  A zero of the function f(x) is computed in the interval ax, bx .

INPUT:
  ax   left endpoint of initial interval
  bx   right endpoint of initial interval
  f   function subprogram which evaluates f(x) for any x in the interval ax, bx
  tol  desired length of the interval of uncertainty of the final result ( >= 0)

OUTPUT:
  zeroin abcissa approximating a zero of f in the interval ax, bx

    It is assumed that  f(ax)  and  f(bx)  have opposite signs without a check.
  zeroin returns a zero x in the given interval [ax, bx] to within a tolerance
  4*macheps*abs(x) + tol, where macheps is the relative machine precision.
   This function subprogram is a slightly modified translation of the algol 60
  procedure zero given in richard brent, algorithms for minimization without
  derivatives, prentice - hall, inc. (1973).                   }
function zeroin(ax, bx: TFloat; f: TFloatFunction1D; tol: TFloat): TFloat;


{
  Pascal version by Nikolai V. Shokhirev                                   
PUTPOSE:
  An approximation x to the point where f attains a minimum on
  the interval (ax,bx) is determined.

INPUT:
  ax  left endpoint of initial interval
  bx  right endpoint of initial interval
  f   function subprogram which evaluates f(x) for any x in the interval (ax,bx)
  tol  desired length of the interval of uncertainty of the final result ( >= 0)

OUTPUT:
  fmin abcissa approximating the point where f attains a minimum

    The method used is a combination of golden section search and successive
  parabolic interpolation. Convergence is never much slower than that for
  a fibonacci search. If f has a continuous second derivative which is positive
  at the minimum (which is not at ax or bx), then convergence is superlinear,
  and usually of the order of about 1.324....
    The function f is never evaluated at two points closer together than
  eps*abs(fmin) + (tol/3), where eps is approximately the square root of
  the relative machine precision.  if  f  is a unimodal function and the
  computed values of  f  are always unimodal when separated by at least
  eps*abs(x) + (tol/3), then fmin approximates the abcissa of the global minimum
  of f on the interval ax,bx with  an error less than 3*eps*abs(fmin) + tol.
    If  f  is not unimodal, then fmin may approximate a local, but perhaps
  non-global, minimum to the same accuracy.
    This function subprogram is a slightly modified version of the algol 60
  procedure localmin given in richard brent, algorithms for minimization without
  derivatives, prentice - hall, inc. (1973). }
function FMMfmin(ax, bx: TFloat; f: TFloatFunction1D; tol: TFloat): TFloat;


{
  pascal version by Eugene B. Krissinel, 1991                              
  Singular Value Decomposition - translation from FORTRAN                  
  Ref: G.E.Forsythe, M.A.Malcolm, C.B.Moler.                               
       Computer methods for mathematical computations. Prentice-Hall, 1977. 
PURPOSE:
  This subroutine determines the singular value decomposition              
  A  =  U * W * VT of a real M by N rectangular matrix.                    
  Householder bidiagonalization and a variant of the QR algorithm are used.
  SVD of the matrix A:  A  =    U  *  W  * VT                              
                      N1xN2   N1xN2 N2xN2 N2xN2                            
  Input:                                                                   
    Case N1 > N2                                                           
      NA = M = N1, N = N2                                                  
      Dimensions: A[M,N],W[N],U[M,N],V[M,N],RV1[N]                         
    Case N1 < N2                                                           
      NA = N = N1, M = N2                                                  
      Dimensions: A[N,M],W[N],U[M,N],V[M,N],RV1[N]                         
    *** Always M >= N                                                      
    A contains the rectangular input matrix to be decomposed.              
    MatU should be set to true if the U matrix in the                      
         decomposition is desired, and to false otherwise.                 
    MatV should be set to true if the V matrix in the                      
         decomposition is desired, and to false otherwise.                 
  Output:                                                                  
    A is unaltered (unless overwritten by U or V).                         
    W contains the N (non-negative) singular values of A (the diagonal     
      elements of s).  They are unordered.  If an error exit is made,      
      the singular values should be correct for indices ierr+1,ierr+2,...,N.
    U contains N orthogonal column U-vectors of the decomposition,         
      if MatU has been set to true. Otherwise U is used as a temporary array.
      U may coincide with A.                                               
      If an error exit is made, the columns of U corresponding             
       to indices of correct singular values should be correct.            
      *** In the case N1 < N2 the last M-N rows contain zero               
    V contains N orthogonal column V-vectors of the decomposition,         
      if MatV has been set to true. Otherwise V is not referenced.         
      V may also coincide with A if U is not needed.                       
      If an error exit is made, the columns of V corresponding             
      to indices of correct singular values should be correct.             
      *** In the case N1 > N2 the last M-N rows are zero                   
    RetCode                                                                
           = 0  for normal exit,                                           
           = k  if the k-th singular value has not been                    
                                          determined after 30 iterations.  
    RV1[N] is a temporary storage array.                                   
  This is a modified version of a routine from the eispack collection
  by the NATS project modified to eliminate machep                       }
procedure FMMSVD(NA, M, N: integer; const A,U,V: IFArr2D; const W,RV1: IFArr1D;
                                      MatU,MatV: boolean; var RetCode: integer);

implementation

uses
  math;

procedure FMMSolve(const a: IFArr2D; const b: IFArr1D; const ipvt: IIArr1D);
var
  n, kb, km1, nm1, kp1, i, k, m: TInt;
  t: TFloat;
begin
  // no Lim checks
  n := a.Dim1;
  if (n > 1) then
  begin
    // forward elimination
    nm1 := n-1;
    for k := 1 to nm1 do
    begin
      kp1 := k+1;
      m := ipvt[k];
      t := b[m];
      b[m] := b[k];
      b[k] := t;
      for i := kp1 to n do
        b[i] := b[i] + a[i,k]*t;
    end;
    // back substitution
    for kb := 1 to nm1 do
    begin
       km1 := n-kb;
       k := km1+1;
       b[k] := b[k]/a[k,k];
       t := -b[k];
       for i := 1 to km1 do
         b[i] := b[i] + a[i,k]*t;
    end;
  end;
  b[1] := b[1]/a[1,1];
end;

procedure FMMDecomp(const a: IFArr2D; var cond: TFloat; const ipvt: IIArr1D);
var
  work: IFArr1D;
  ek, t, anorm, ynorm, znorm: TFloat;
  n, nm1, i, j, k, kp1, kb, km1, m: TInt;
begin
  // no Lim checks
  n := a.Dim1;
  ipvt[n] := 1;
  if (n = 1) then
  begin// 1-by-1
    if (a[1,1] <> 0.0) then cond := 1.0 else cond := 1.0e+32;
    exit;
  end;

  work := TFArr1D.Create(n); // work space..
  //     compute 1-norm of a
  anorm := 0.0;
  for j := 1 to n do
  begin
    t := 0.0;
    for i := 1 to n do
      t := t + abs(a[i,j]);
    if (t  >  anorm) then anorm := t;
  end;

  //     gaussian elimination with partial pivoting
  nm1 := n - 1;
  for k := 1 to nm1 do
  begin
    kp1:= k+1;

    // find pivot
    m := k;
    for i := kp1 to n do
    begin
      if (abs(a[i,k]) > abs(a[m,k])) then m := i;
    end;
    ipvt[k] := m;
    if (m <> k) then ipvt[n] := -ipvt[n];
    t := a[m,k];
    a[m,k] := a[k,k];
    a[k,k] := t;
//     if (t  =  0.0) goto 35;
    // skip step if pivot is zero
    if (t <> 0.0) then
    begin
      // compute multipliers
      for i := kp1 to n do
        a[i,k] := -a[i,k]/t;
        // interchange and eliminate by columns
       for j := kp1 to n do
       begin
         t := a[m,j];
         a[m,j] := a[k,j];
         a[k,j] := t;
         if (t <> 0.0) then
           for i := kp1 to n do
              a[i,j] := a[i,j] + a[i,k]*t;
      end;
    end;
  end;
{ cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)
  Estimate obtained by one step of inverse iteration for the small singular
  vector. This involves solving two systems of equations,
    (a-transpose)*y = e  and  a*z = y
  where e is a vector of +1 or -1 chosen to cause growth in y.
  estimate = (1-norm of z)/(1-norm of y) }

//     solve (a-transpose)*y = e
  for k := 1 to  n do
  begin
    t := 0.0;
    if (k > 1) then
    begin
      km1 := k-1;
      for i := 1 to  km1 do
        t := t + a[i,k]*work[i];
    end;
    ek := 1.0;
    if (t < 0.0) then ek := -1.0;
    if (a[k,k] = 0.0) then
    begin
      cond := 1.0e+32;
      exit;
    end;
    work[k] := -(ek + t)/a[k,k];
  end;

  for kb := 1 to nm1 do
  begin
    k := n - kb;
    t := 0.0;
    kp1 := k+1;
    for i := kp1 to  n do
      t := t + a[i,k]*work[i];
    work[k] := t + work[k];
    m := ipvt[k];
    if (m <> k) then
    begin
       t := work[m];
       work[m] := work[k];
       work[k] := t;
    end;
  end;

  ynorm := 0.0;
  for i := 1 to  n do
    ynorm := ynorm + abs(work[i]);

//     solve a*z = y
  FMMSolve(a, work, ipvt);

  znorm := 0.0;
  for i := 1 to  n do
    znorm := znorm + abs(work[i]);

  // estimate condition
  cond := anorm*znorm/ynorm;
  if (cond < 1.0) then cond := 1.0;

end;


procedure FMMSpline(const x, y, b, c, d: IFArr1D);
const
  ZERO  = 0.0;
  TWO   = 2.0;
  THREE = 3.0;
var
  k, Dim, Lo, Hi: TInt;
  t: TFloat;
begin
  Dim := x.Dim1;
  Lo := x.Lo1;
  Hi := x.Hi1;

  if (Dim < 3) then
  begin // Straight line - special case for n < 3
    b[Lo] := ZERO;
    if (Dim = 2) then b[Lo] := (y[Lo+1]-y[Lo])/(x[Lo+1]-x[Lo]);
    c[Lo] := ZERO;
    d[Lo] := ZERO;
    if (Dim < 2) then exit;
    b[Lo+1] := b[Lo];
    c[Lo+1] := ZERO;
    d[Lo+1] := ZERO;
    exit;
  end;

//.....Set up tridiagonal system.........................................
//.     b = diagonal, d = offdiagonal, c = right-hand side
  d[Lo]   := x[Lo+1]-x[Lo];
  c[Lo+1] := (y[Lo+1]-y[Lo])/d[Lo];
//  for k  :=  2 to n-1 do
  for k := Lo+1 to Hi-1 do
  begin
    d[k]   := x[k+1]-x[k];
    b[k]   := TWO*(d[k-1]+d[k]);
    c[k+1] := (y[k+1]-y[k])/d[k];
    c[k]   := c[k+1]-c[k];
  end;

//.....End conditions.  third derivatives at x[Lo] and x[Hi] obtained
//.       from divided differences.......................................
  b[Lo] := -d[Lo];
  b[Hi] := -d[Hi-1];
  c[Lo] := ZERO;
  c[Hi] := ZERO;
  if (Dim > 3) then
  begin
    c[Lo] := c[Lo+2]/(x[Lo+3]-x[Lo+1]) - c[Lo+1]/(x[Lo+2]-x[Lo]);
    c[Hi] := c[Hi-1]/(x[Hi]-x[Hi-2]) - c[Hi-2]/(x[Hi-1]-x[Hi-3]);
    c[Lo] := c[Lo]*d[Lo]*d[Lo]/(x[Lo+3]-x[Lo]);
    c[Hi] := -c[Hi]*d[Hi-1]*d[Hi-1]/(x[Hi]-x[Hi-3]);
  end;

  for  k  := Lo+1 to Hi do
  begin // forward elimination
    t := d[k-1]/b[k-1];
    b[k] := b[k]-t*d[k-1];
    c[k] := c[k]-t*c[k-1];
  end;

  c[Hi] := c[Hi]/b[Hi];
// back substitution ( makes c the sigma of text)
  for  k := Hi-1 downto Lo do
  begin
    c[k] := (c[k]-d[k]*c[k+1])/b[k];
  end;

//.....Compute polynomial coefficients...................................
  b[Hi] := (y[Hi]-y[Hi-1])/d[Hi-1]+d[Hi-1]*(c[Hi-1]+c[Hi]+c[Hi]);
  for  k := Lo to Hi-1  do
  begin
    b[k] := (y[k+1]-y[k])/d[k]-d[k]*(c[k+1]+c[k]+c[k]);
    d[k] := (c[k+1]-c[k])/d[k];
    c[k] := THREE*c[k];
  end;
  c[Hi] := THREE*c[Hi];
  d[Hi] := d[Hi-1];

  exit;
end;// of FMMspline ---------------------------------------------------------


procedure NaturalSpline(const x, y, b, c, d: IFArr1D);
const
  ZERO  = 0.0;
  TWO   = 2.0;
  THREE = 3.0;
var
  k, Dim, Lo, Hi: TInt;
//  t: TFloat;
begin
  Lo := x.Lo1;
  Hi := x.Hi1;
  Dim := x.Dim1;

  if (Dim < 3) then
  begin // Straight line - special case for Dim < 3
    b[Lo] := ZERO;
    if (Dim = 2) then b[Lo]:=(y[Lo+1]-y[Lo])/(x[Lo+1]-x[Lo]);
    c[Lo] := ZERO;
    d[Lo] := ZERO;
    b[Lo+1] := b[Lo];
    c[Lo+1] := ZERO;
    d[Lo+1] := ZERO;
    exit;
  end;

//  d[1:n-1] := x[2:n]-x[1:n-1]; - Put the h-array of the text into array d
  for k := Lo to Hi-1 do
    d[k] := x[k+1]-x[k];
//.....Set up the upper triangular system in locations 2 through Dim-1 of
//        arrays b and c. B holds the diagonal and c the right hand side.
  b[Lo+1] := TWO*(d[Lo]+d[Lo+1]);
  c[Lo+1] := (y[Lo+2]-y[Lo+1])/d[Lo+1]-(y[Lo+1]-y[Lo])/d[Lo];
  for k := Lo+2 to Hi-1  do
  begin
    b[k] := TWO*(d[k-1]+d[k])-d[k-1]*d[k-1]/b[k-1];
    c[k] := (y[k+1]-y[k])/d[k]-(y[k]-y[k-1])/d[k-1]-d[k-1]*c[k-1]/b[k-1];
  end;

  c[Hi-1] := c[Hi-1]/b[Hi-1];
// Back substitute to get c-array
  for  k := Hi-2 downto Lo+1  do begin
    c[k] := (c[k]-d[k]*c[k+1])/b[k];
  end;
  c[Lo] := ZERO;
  c[Hi] := ZERO   ;
// c now holds the sigma array of the text


//.....Compute polynomial coefficients ..................................
  b[Hi] := (y[Hi]-y[Hi-1])/d[Hi-1]+d[Hi-1]*(c[Hi-1]+c[Hi]+c[Hi]);
  for  k := Lo to Hi-1 do
  begin
    b[k] := (y[k+1]-y[k])/d[k]-d[k]*(c[k+1]+c[k]+c[k]);
    d[k] := (c[k+1]-c[k])/d[k];
    c[k] := THREE*c[k];
  end;
  c[Hi] := THREE*c[Hi];
  d[Hi] := d[Hi-1];
end;//  --------------------------------------------

function Seval(const u: TFloat; var i: TInt; const x, y, b, c, d: IFArr1D): TFloat;
var
  k, Lo, Hi, j: TInt;
  dx: TFloat;
begin
//  Dim := SIZE(x);
  Lo := x.Lo1;
  Hi := x.Hi1;

//.....First check if u is in the same interval found on the last call to Seval
  if (i < Lo) then  i := Lo;
  if (i >= Hi) then  i := Hi;
  if ( (u < x[i]) or (u >= x[i+1]) ) then
  begin
    i := Lo ;
    j := Hi+1;
    // binary search
    repeat
      k := (i+j) div 2;
      if (u < x[k]) then j := k else i := k;
    until (j <= i+1);
  end;
  dx := u-x[i];
// evaluate the spline
  Result := y[i]+dx*(b[i]+dx*(c[i]+dx*d[i]));
end;//  -------------------------------------------------------

procedure Seval3(const u: TFloat; var i: TInt; const x, y, b, c, d: IFArr1D;
                       var f, fp, fpp, fppp: TFloat);
const
  TWO = 2.0;
  THREE = 3.0;
  SIX = 6.0;
var
  k, Lo, Hi,j: TInt;
  dx: TFloat;
begin
//  Dim := SIZE(x);
  Lo := x.Lo1;
  Hi := x.Hi1;

//.....First check if u is in the same interval found on the
//        last call to Seval.............................................
  if (i < Lo) then  i := Lo;
  if (i >= Hi) then  i := Hi;
  if ( (u < x[i]) or (u >= x[i+1]) ) then
  begin
    i := Lo ;
    j := Hi+1;
    // binary search
    repeat
      k := (i+j) div 2;
      if (u < x[k]) then j := k else i := k;
    until (j <= i+1);
  end;
  dx := u-x[i]   ;
// evaluate the spline
  f:=y[i]+dx*(b[i]+dx*(c[i]+dx*d[i]));
  fp:=b[i]+dx*(TWO*c[i] + dx*THREE*d[i]);
  fpp:=TWO*c[i] + dx*SIX*d[i];
  fppp:=SIX*d[i];
end;//// -------------------------------------------------------

{
PURPOSE
  Estimate the integral of FUN(X) from A to B to a user provided tolerance.
  This is an automatic adaptive routine based on the 8-panel Newton-Cotes rule.
INPUT
  a - lower limit of integration
  b - upper limit of integration
  abserr - absolute error tolerance (>0)
  relerr - relative error tolerance (>0)
OUTPUT
  result - approx. value of integral
  errest - estimate of actual error
  nofun  - no. of function evaluations
  flag   - reliability indicator. (=0 is O.K.)   }
procedure Quanc8(Fun: TFloatFunction1D; const a, b, abserr, relerr: TFloat;
                            var result, errest, flag: TFloat; var nofun: TInt);
label
  30, 50, 70, 80;

const
  ZERO: TFloat = 0.0;
  HALF = 0.5;
  ONE = 1.0;
  SIXTEEN =16.0;

  LEVMIN = 1;
  LEVOUT = 6;
  NOMAX = 5000;
{
  W0 =  0.27908289241622574955908289241623; // =   3956.0/14175.0;
  W1 =  1.6615167548500881834215167548501;  // =  23552.0/14175.0;
  W2 = -0.26186948853615520282186948853616; // =  -3712.0/14175.0;
  W3 =  2.9618342151675485008818342151675;  // =  41984.0/14175.0;
  W4 = -1.2811287477954144620811287477954;  // = -18160.0/14175.0;
}
var
  area: TFloat;
  cor11: TFloat;
  esterr: TFloat;
  f, x: array[0..16] of TFloat;
  fsave, xsave: array[1..8,1..30] of TFloat;
  i, j, jj: TInt;
  lev, nim: TInt;
  NOFIN: TInt;  //...Trouble if NOFUN reaches NOFIN
  qprev, qnow, qdiff, qleft: TFloat;
  stone, step: TFloat;
  tolerr: TFloat;
  qright: array[1..31] of TFloat;
  LEVMAX: TInt;
  W0, W1, W2, W3, W4: TFloat;


begin
//..... S t a g e   1  (General initialization)..........................
  LEVMAX  :=  30;
  W0 :=   3956.0/14175.0;
  W1 :=  23552.0/14175.0;
  W2 :=  -3712.0/14175.0;
  W3 :=  41984.0/14175.0;
  W4 := -18160.0/14175.0;
//   NOFIN := NOMAX-8*(LEVMAX-LEVOUT+2**(LEVOUT+1));
  NOFIN := NOMAX - 8*(LEVMAX - LEVOUT + Ipower(2, LEVOUT+1));

//...initialize running sums to zero...
  area  := ZERO;
  cor11 := ZERO;
  flag  := ZERO;
  result := ZERO;
  errest := ZERO;
  nofun  := 0;
  if (a = b) then exit;

//..... S T A G E   2   (initialization for first interval)..............
  lev := 0;
  nim := 1;
  x[0] := a;
  x[16] := b;
  qprev := ZERO;
  stone := (b-a)/SIXTEEN;
  x[8] :=  HALF*(x[0]+x[16]);
  x[4] :=  HALF*(x[0]+x[8]);
  x[12] := HALF*(x[8]+x[16]);
  x[2] :=  HALF*(x[0]+x[4]);
  x[6] :=  HALF*(x[4]+x[8]);
  x[10] := HALF*(x[8]+x[12]);
  x[14] := HALF*(x[12]+x[16]);

//  for  j := 0 to 16 step 2 do
//  f[j] := Fun(x[j]);
//  for  j := 0 to 8 do
  for  j := 0 to 8 do
    f[j*2] := Fun(x[j*2]);
  nofun := 9;

//..... S T A G E   3   (central calculation)............................
//...requires qprev, x0, x2, x4,....x[16], f0, f2,f 4,....f16.
//...calculates x1, x3,....x15, f1, f3,....f15, qleft, qright, qnow, qdiff, area

// 30for  j := 1 to 15 step 2 do begin
30: j := 1;
    while j <= 15 do
    begin // keeps coming back here until successful
      x[j] := HALF*(x[j-1] + x[j+1]);
      f[j] := Fun(x[j]);
      inc(j,2);
    end;
    nofun := nofun+8;
    step := (x[16]-x[0])/SIXTEEN;
    qleft := (w0*(f[0]+f[8]) +w1*(f[1]+f[7])+ w2*(f[2]+f[6])+
              w3*(f[3]+f[5]) +w4*f[4])*step;
    qright[lev+1] := (w0*(f[8] +f[16]) +w1*(f[9]+f[15]) +w2*(f[10]+f[14])+
                      w3*(f[11]+f[13]) +w4*f[12])*step;
    qnow := qleft + qright[lev+1];
    qdiff := qnow - qprev;
    area := area + qdiff;

//..... S T A G E   4   (interval convergence test) .....................
    esterr := ABS(qdiff)/1023.0;
    tolerr := MAX(abserr,relerr*ABS(area))*(step/stone);
    if (lev < LEVMIN) then goto 50;

    if (lev >= LEVMAX) then
    begin
      flag := flag+ONE;
      goto 70;
    end;

    if (nofun > NOFIN) then
    begin //..... S T A G E   6  (trouble section)
      nofin := nofin + nofin;
// Number of function values is about to exceed limit
      levmax := levout;
      flag := flag+(b-x[0])/(b-a);
      goto 70;
    end;

    if (esterr <= TOLERR) then goto 70;

//..... S T A G E   5   (no convergence).................................
//...locate next interval...
50: nim := nim + nim;
    lev := lev + 1;

// store right hand elements for future use
//    fsave[1:8,lev] := f[9:16];
//    xsave[1:8,lev] := x[9:16];
    for i := 1 to 8 do
    begin
      fsave[i,lev] := f[i+8];
      xsave[i,lev] := x[i+8];
    end;

    qprev := qleft;
    for  I  :=  1 to 8  do
    begin // assemble left hand elements for immediate use
      f[18-2*i] := f[9-i];
      x[18-2*i] := x[9-i];
    end;
    goto 30;

//..... S T A G E   7   (integral converged) ............................
//...add contributions into running sums...
70: result := result+qnow;
   errest := errest+esterr;
   cor11 := cor11+qdiff/1023.0;

  while odd(nim) do
  begin // 72 loop in text
//    if (nim = 2*(nim/2) then) EXIT   ;
    nim := nim div 2;
    lev := lev-1;
  end;
  nim := nim + 1;

  if (lev <= 0) then goto 80;
// looks like success!
  qprev := qright[lev]   ;
// assemble elements required for next interval
  x[0] := x[16];
  f[0] := f[16];
  for  i  :=  1 to 8  do
  begin
    f[2*i] := fsave[i,lev];
    x[2*i] := xsave[i,lev];
  end;

  goto 30;

//..... S T A G E   8  (Finalize and return).............................
80: result := result+cor11;
  if (errest <> ZERO) then
  begin  // make sure ERREST is not < roundoff level ....
    while (abs(result)+errest) = abs(RESULT) do
    begin
      errest := errest + errest;
    end;
  end;

end;// ---------------------------------------------------

{
  Pascal version by Nikolai V. Shokhirev
Ref.
  George E. Forsythe, Michael A. Malcolm, and Cleve B. Moler.
  "Computer Methods for Mathematical Computations". Prentice-Hall, 1977.
PURPOSE:
  Integrate a system of neqn first-order ordinary differential equations
  of the form dy/dt=F(t,y) [y is a vector], where the initial values of y
  and the initial values of yp, the derivatives are specified at the starting
  point t.  Fehl advances the solution over the fixed step h and returns the
  fifth order (sixth-order locally) solution approximation at t+h in array s.
  neqn: TInt; // number of equations
INPUT
  t - starting point
  h - step size
  y - array of length neqn; function at t
  yp - array of length neqn; derivatives at t
OUTPUT
  f1 - array of length neqn for internal use
  f2 - array of length neqn for internal use
  f3 - array of length neqn for internal use
  f4 - array of length neqn for internal use
  f5 - array of length neqn for internal use
  s  - array of length neqn; the results         }
procedure Fehl(F: TDiffEqProc; const t, h: TFloat;
                               const y, yp, f1, f2, f3, f4, f5, s : IFarr1D);
const
  C1  = 0.25;                // = 0.25;
  C2  = 0.09375;             // = 3.0/32.0;
  C3  = 3.0;                 // = 3.0;
  C4  = 0.375;               // = 3.0/8.0;
  C5  = 4.55166135639508e-4; // = 1.0/2197.0;
  C6  = 1932.0;              // = 1932.0;
  C7  = 7296.0;              // = 7296.0;
  C8  = -7200.0;             // = -7200.0;
  C9  = 0.923076923076923;   // = 12.0/13.0;
  C10 = 2.43664717348928e-4; // = 1.0/4104.0;
  C11 = 8341.0;              // = 8341.0;
  C12 = -845.0;              // = -845.0;
  C13 = 29440.0;             // = 29440.0;
  C14 = -32832.0;            // = -32832.0;
  C15 = 4.87329434697856e-5; // = 1.0/20520.0;
  C16 = -6080.0;             // = -6080.0;
  C17 = 9295.0;              // = 9295.0;
  C18 = -5643.0;             // = -5643.0;
  C19 = 41040.0;             // = 41040.0;
  C20 = -28352.0;            // = -28352.0;
  C21 = 0.5;                 // = 0.5;
  C22 = 1.31267187797402e-7; // = 1.0/7618050.0;
  C23 = 902880.0;            // = 902880.0;
  C24 = 3855735.0;           // = 3855735.0;
  C25 = -1371249.0;          // = -1371249.0;
  C26 = 3953664.0;           // = 3953664.0;
  C27 = 277020.0;            // = 277020.0;

var
//  neqn: TInt; // number of equations
//  t: TFloat;  // starting point
//  h: TFloat;  // step size

  ch: TFloat;
  i, Lo, Hi: TInt;
{
  y: IFarr1D; // array of length neqn; function at t
  yp: IFarr1D; // array of length neqn; derivatives at t
  f1: IFarr1D; // array of length neqn for internal use
  f2: IFarr1D; // array of length neqn for internal use
  f3: IFarr1D; // array of length neqn for internal use
  f4: IFarr1D; // array of length neqn for internal use
  f5: IFarr1D; // array of length neqn for internal use
  s : IFarr1D; // array of length neqn; the results
  C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14,
  C15,  C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26, C27: TFloat;
  C1  := 0.25;
  C2  := 3.0/32.0;
  C3  := 3.0;
  C4  := 3.0/8.0;
  C5  := 1.0/2197.0;
  C6  := 1932.0;
  C7  := 7296.0;
  C8  := -7200.0;
  C9  := 12.0/13.0;
  C10 := 1.0/4104.0;
  C11 := 8341.0;
  C12 := -845.0;
  C13 := 29440.0;
  C14 := -32832.0;
  C15 := 1.0/20520.0;
  C16 := -6080.0;
  C17 := 9295.0;
  C18 := -5643.0;
  C19 := 41040.0;
  C20 := -28352.0;
  C21 := 0.5;
  C22 := 1.0/7618050.0;
  C23 := 902880.0;
  C24 := 3855735.0;
  C25 := -1371249.0;
  C26 := 3953664.0;
  C27 := 277020.0;   }
begin  //-----------------------------------------------------------------------
  ch := C1*h;

  Lo := y.Lo1;
  Hi := y.Hi1;
  for i := Lo to Hi do
    f5[i] := y[i]+ch*yp[i];
  F(t+ch, f5, f1);

  ch := C2*h;
  for i := Lo to Hi do
    f5[i] := y[i] + ch*(yp[i]+C3*f1[i]);
  F(t+C4*h, f5, f2);

  ch := C5*h;
  for i := Lo to Hi do
    f5[i] := y[i] + ch*(C6*yp[i] + (C7*f2[i]+C8*f1[i]) );
  F(t+C9*h, f5, f3);

  ch := C10*h;
  for i := Lo to Hi do
    f5[i] := y[i] + ch*( (C11*yp[i] + C12*f3[i]) + (C13*f2[i]+C14*f1[i]) );
  F(t+h, f5, f4);

  ch := C15*h;
  for i := Lo to Hi do
    f1[i] := y[i] + ch*( C16*yp[i] + (C17*f3[i]+C18*f4[i]) + (C19*f1[i]+C20*f2[i]) );
  F(t+C21*h, f1, f5);

//.....compute approximate solution at t+h
  ch := C22*h;
  for i := Lo to Hi do
    s[i] := y[i] + ch*( C23*yp[i] + (C24*f3[i]+C25*f4[i]) + (C26*f2[i]+C27*f5[i]) );

end;// of Fehl -----------------------------------------------------


{
PURPOSE:
  Integrate a system of neqn first order differential equations by the Fehlberg
  fourth-fifth order Runge-Kutta method.
NOTES:
  RKF45 is primarily designed to solve non-stiff and mildly stiff differential
  equations when derivative evaluations are inexpensive.
  Rkf45 should generally not be used when the user is demanding high accuracy.

  Rkfs integrates a system of first order ordinary differential equations as
  described in the comments for Rkf45. The arrays yp, f1, f2, f3, f4, and f5
  (of dimension neqn) and the variables h, savre, savae, nfe, kop,
  init, jflag, and kflag are used internally by the code and appear in the call
  list to eliminate local retention of variables between calls. Accordingly,
  they should not be altered. Items of possible interest are
    yp - derivative of solution vector at t
    h - an appropriate stepsize to be used for the next step
    nfe- counter on the number of derivative function evaluations

ABSTRACT:
  Subroutine Rkf45 integrates a system of neqn first order ordinary differential
  equations of the form
 	  dy[i]/dt = f(t,y[1],y[2],...,y[neqn])
 	where the y[i] are given at t .
  Typically the subroutine is used to integrate from t to tout but it can be
  used as a one-step integrator to advance the solution a Double step in the
  direction of tout. On return, the parameters in the call list are set for
  continuing the integration. The user has only to call Rkf45 again (and perhaps
  define a new value for tout). Rkfs45 calls subroutine Fehl which computes an
  approximate solution over one step.

  Rkf45 uses the Runge-Kutta-Fehlberg(4,5) method described in the reference
  E. Fehlberg:  Low-order Classical Runge-Kutta formulas with Stepsize 
  Control, NASA TR R-315.

  The performance of Rkf45 is illustrated in the reference
  L.F.Shampine,H.A.Watts,S.Davenport: "Solving Non-stiff Ordinary
  Differential Equations - The State of the Art.
   Sandia Laboratories Report SAND75-0182, to appear in SIAM Review. 

  First call to Rkf45:

  The user must provide storage in his calling program for the arrays in
  the call list y, yp, f1, f2, f3, f4, f5 - 1D float arrays with the same
  limits and dimension = neqn, supply procedure F(t,y,yp) and initialize the
  following parameters:
    y - vector of initial conditions
    t - starting point of integration , must be a variable
    tout - output point at which solution is desired.
    t = tout is allowed on the first call only, in which case Rkf45 returns with
      iflag=2 if continuation is possible.
    relerr, abserr - relative and absolute local error tolerances which must be
      non-negative. relerr must be a variable while abserr may be a constant.
      The code should normally not be used with relative error control smaller
      than about i.e-8 to avoid limiting precision difficulties the code
      requires relerr to be larger than an internally computed relative error
      parameter which is machine dependent. In particular, pure absolute error
      is not permitted. If a smaller than allowable value of relerr is
      attempted, rkf45 increases relerr appropriately and returns control to the
      user before  continuing the integration.
    iflag - +1,-1 indicator to initialize the code for each new problem.
      Normal input is +1. the user should set iflag = -1 only when one-step
      integrator control is essential. In this case, Rkf45 attempts to advance
      the solution a Double step in the direction of tout each time it is
      called. Since this mode of operation results in extra computing overhead,
      it should be avoided unless needed.

  Output from Rkf45

  y - solution at t
  t - last point reached in integration.
  iflag =  2 - integration reached tout. Indicates successful return and is the
               normal mode for continuing integration.
       	= -2 - a Double successful step in the direction of tout has been taken.
               Normal mode for continuing integration one step at a time.
 	      =  3 - integration was not completed because relative error tolerance
               was too small. relerr has been increased	appropriately for
               continuing.
 	      =  4 - integration was not completed because more than 3000 derivative
               evaluations were needed. This is approximately 500 steps.
 	      =  5 - integration was not completed because solution vanished making
               a pure relative error test impossible. Must use non-zero abserr
               to continue. Using the one-step integration mode for one step
               is a good way to proceed.
 	      =  6 - integration was not completed because requested accuracy could
               not be achieved using smallest allowable stepsize. User must
               increase the error tolerance before continued integration can
               be attempted.
 	      =  7 - it is likely that Rkf45 is inefficient for solving this problem.
               Too much output is restricting the natural stepsize choice.
               Use the one-step integrator mode.
 	      =  8 - invalid input parameters. This indicator occurs if any of the
               following is satisfied:
                 t = tout and iflag <> +1 or -1
 	               relerr or abserr < 0.
                 iflag = 0 or < -2 or > 8
  h, nfe, yp, f1, f2, f3, f4, f5 - information which is usually of no interest
  to the user but necessary for subsequent calls.
  - yp contain the first derivatives of the solution vector y at t.
  - h contains the stepsize to be attempted on the next step.
  - nfe contains the derivative evaluation counter.

  Subsequent calls to Rkf45

  Subroutine Rkf45 returns with all information needed to continue the
  integration. If the integration reached tout, the user need only define
  a new tout and call Rkf45 again. In the one-step integrator mode (iflag=-2)
  the user must keep in mind that each step taken is in the direction of the
  current tout. Upon reaching tout (indicated by changing iflag to 2),the 
  user must then define a new tout and reset iflag to -2 to continue in the 
  one-step integrator mode. 

  If the integration was not completed but the user still wants to continue 
  (iflag = 3,4 cases), he just calls Rkf45 again. With iflag = 3, the relerr
  parameter has been adjusted appropriately for continuing the integration.
  In the case of iflag=4 the function counter will be reset to 0 and 
  another 3000 function evaluations are allowed. 

  However,in the case iflag = 5, the user must first alter the error criterion
  to use a positive value of abserr before integration can proceed. If he 
  does not, execution is terminated. 

  Also, in the case iflag = 6, it is necessary for the user to reset iflag
  to 2 (or -2 when the one-step integration mode is being used) as well as 
  increasing either abserr,relerr or both before the integration can be 
  continued. If this is not done, execution will be terminated. The 
  occurrence of iflag=6 indicates a trouble spot (solution is changing 
  rapidly,singularity may be present) and it often is inadvisable to continue.

  If (flag = 7 is encountered, the user should use the one-step integration mode
  with the stepsize determined by the code or consider switching to the 
  Adams codes DE/STEP, INTRP. If the user insists upon continuing the
  integration with Rkf45, he must reset iflag to 2 before calling Rkf45 again.
  Otherwise, execution will be terminated.

  If iflag = 8 is obtained, integration can not be continued unless the invalid
  input parameters are corrected. 

  It should be noted that the arrays work,iwork contain information required
  for subsequent integration. Accordingly, work and iwork should not be 
  altered.

  This interfacing routine merely relieves the user of a long calling list 
  via the splitting apart of two working storage arrays.
INPUT
  y - solution vector at t
  t - independent variable
  tout - output point at which solution is desired
  relerr - relative error tolerance
  abserr -  absolute error tolerance
OUTPUT
  iflag - indicator for status of work
  h - step size
  savre, savae
  nfe - number of function evaluations
  kop;
IN, OUT
  yp - derivatives
  f1, f2, f3, f4, f5
  init;
  jflag;
  kflag;                               }
procedure Rkfs45(F: TDiffEqProc; var t: TFloat; const tout, abserr: TFloat;
                 var relerr: TFloat; var iflag: TInt;
                 const y, yp, f1, f2, f3, f4, f5: IFArr1D;
                 var h, savre, savae: TFloat;
                 var nfe, kop, init, jflag, kflag: TInt);
label
  25, 30, 40, 45, 50, 55, 60, 65, 80, 85, 95, 100, 150, 200, 220, 240, 260, 300;
const
  ZERO = 0.0;
  ONE = 1.0;
  REMIN = 1E-12; { - minimum acceptable value of relerr.
                   Attempts to obtain higher accuracy with this routine are
                   usually very expensive and unsuccessful.}
  MAXNFE = 3000; { - the expense is controlled by restricting the number of
                   function evaluations to be approximately MAXNFE.
                   As set, this corresponds to about 500 steps. }
var
  rer, s          : TFloat;
  scale, tol, toln: TFloat;
  ypk, eps, u26: TFloat;
  a, ae, dt, ee, eeoet, esttol, et: TFloat;
  hmin: TFloat;
  k, mflag: TInt;
  Lo, Hi, neqn: TInt;
  hfaild, output: Boolean;
begin
//.....check input parameters............................................
  Lo := y.Lo1;
  Hi := y.Hi1;
  neqn := y.Dim1;
  eps := cMachEps;
  u26 := 26.0*eps;
//  if (neqn < 1) or (relerr < ZERO) or (abserr < ZERO) then
  if (relerr < ZERO) or (abserr < ZERO) then
  begin
    iflag := 8;
    exit;
  end;
  mflag := ABS(iflag);
  if (mflag = 0) or (mflag > 8) then
  begin
    iflag := 8;
    exit;
  end;

  if (mflag = 1) then goto 50;

//.....check continuation posibilities...................................
  if ((t = tout) and (kflag <> 3)) then
  begin
    iflag := 8; // invalid flag handler
    exit;
  end;
  if (mflag <> 2) then goto 25;

//...you get here if iflag=+2 or -2
  if ((kflag = 3) or  (init = 0)) then goto 45;
  if (kflag = 4) then goto 40;
  if ((kflag = 5) and  (abserr = ZERO)) then  goto 30;
  if ((kflag = 6) and  (relerr <= savre)  and  (abserr <= savae)) then goto 30;
  goto 50;

//...you get here if iflag = 3, 4, 5, 6, 7, or 8.....
25: if (iflag = 3) then goto 45;
    if (iflag = 4) then goto 40;
    if ((iflag = 5) and (abserr > ZERO)) then goto 45;

//.....integration cannot be continued since user did not respond to the
//..... instructions pertaining to iflag = 5, 6, 7, or 8
30: exit;
//!! maybe this is too strong

//.....reset function evaluation counter.................................
40: nfe := 0;
    if (mflag = 2) then goto 50;

//.....reset flag value from previous call...............................
45: iflag := jflag;
    if (kflag = 3) then mflag := ABS(iflag);

50: jflag := iflag;  // save for subsequent calls
    kflag := 0;      // set the continuation flag
    savre := relerr; // save for subsequent calls
    savae := abserr; // save for subsequent calls

//.....restrict relative error tolerance to be at least as large as
//       2*eps+remin to avoid limiting precision difficulties arising
//       from impossible accuracy requests.
    rer := eps + eps + remin;
    if (relerr >= rer) then goto 55;

//.....relative error tolerance too small................................
    relerr := rer;
    iflag := 3;
    kflag := 3;
    exit;

55: dt := tout-t;
    if (mflag = 1) then goto 60;
    if (init = 0) then goto 65;
    goto 80;

//.....initialization....................................................
60: init := 0;
    kop := 0;
    a := t;
    F(a,y,yp);
    nfe := 1;
    if (t <> tout) then goto 65;
    iflag := 2;
    exit;

65: init := 1;
    h := ABS(dt);
    toln := ZERO;
    for  k := Lo to Hi do
    begin
      tol := relerr*ABS(y[k])+abserr;
      if (tol <= ZERO) then continue;
      toln := tol;
      ypk := ABS(yp[k]);
      if (ypk*IntPower(h, 5) > tol) then h := Power((tol/ypk), 0.2);
   end;

   if (toln <= ZERO) then h := ZERO;
   h := max(h, u26*MAX(ABS(t), ABS(dt)));
   jflag := ISign2(2,iflag);

//.....set stepsize for integration in the direction from t to tout
80: h := SIGN2(h,dt);

//.....test to see if rkf45 is being severely impacted by too many outputs
   if (ABS(h) >= ABS(dt+dt)) then kop := kop + 1;
   if (kop <> 100) then goto 85;
   kop := 0;
   iflag := 7;
   exit;

85: if (ABS(dt) > u26*ABS(t)) then goto 95;

//.....if too close to output point, extrapolate and return
    for  k := Lo to Hi  do
    begin
      y[k] := y[k]+dt*yp[k];
    end;
    a := tout;
    F(a, y, yp);
    nfe := nfe +1 ;
    goto 300;

//.....initialize output point indicator.................................
95: output := false ;

//.....scale the error tolerances to avoid premature underflow in the
//        error tolerance
    scale := 2.0/relerr;
    ae := scale*abserr;


//.....step by step integration..........................................
100:hfaild :=  false ;
//...set smallest allowable step size...
    hmin := u26*ABS(t);

//.....adjust stepsize if necessary to hit the output point.
//.....look ahead two steps to avoid drastic changes in the stepsize and
//.....thus lessen the impact of output points on the code.
    dt := tout-t;
    if (ABS(dt) >= ABS(h+h)) then goto 200;
    if (ABS(dt) > ABS(h)) then goto 150;


//.....the next successful step will complete the integration to the
//     output point
    output :=  true ;
    h := dt;
    goto 200;

150:h := 0.5*dt;

//.....core integrator for taking a Double step..........................

200:if (nfe <= maxnfe) then goto 220;

//.....too much work...
    iflag := 4;
    kflag := 4;
    exit;

//.....advance an approximate solution over one step of length h
220:Fehl(F,t,h,y,yp,f1,f2,f3,f4,f5,f1);
    nfe := nfe+5;

//.....compute and test allowable tolerances versus local error estimates
//.....  and remove scaling of tolerances. note that relative error is
//.....  is measured with respect to the average of the magnitudes of the
//.....   solution at the beginning and end of the step.
    eeoet := ZERO;
    for  k := Lo to Hi  do
    begin
      et := ABS(y[k])+ABS(f1[k])+ae;
      if (et > ZERO) then goto 240;
//.....Inappropriate error tolerance.....................................
      iflag := 5;
      exit;
240:  ee := ABS( (-2090.0*yp[k]+(21970.0*f3[k]-15048.0*f4[k]))+
                 (22528.0*f2[k]-27360.0*f5[k]) );
      eeoet := MAX(eeoet,ee/et);
    end;
    esttol := ABS(h)*eeoet*scale/752400.0;
    if (esttol <= ONE) then goto 260;

//.....unsuccessful step. reduce the stepsize and start again............
//     (decrease is limited to a factor of 1/10)
    hfaild :=  true ;
    output :=  false ;
    s := 0.1;
    if (esttol < 59049.0) then s := 0.9/Power(esttol,0.2);
    h := s*h;
    if (ABS(h) > hmin) then goto 200;

//.....requested error unobtainable at smallest allowable stepsize.......
    iflag := 6;
    kflag := 6;
    exit;

//.....successful step.   store solution at t+h and evaluate derivatives
//     there..........................
260:t := t+h;
    for  k := Lo to Hi  do
    begin
      y[k] := f1[k];
    end;
    a := t;
    F(a,y,yp);
    nfe := nfe + 1;

// Choose next stepsize.  The increase is limited to a factor of 5.
// If step failure has just occured, next stepsize is not allowed
//  to increase.
    s := 5.0;
    if (esttol > 1.889568E-4) then s := 0.9/Power(esttol,0.2);
    if (hfaild) then s := MIN(s,ONE);
    h := SIGN2(MAX(s*ABS(h),hmin),h);

//.....E N D   O F   C O R E   I N T E G R A T O R   ....................

//.....Should we take another step? .....................................
    if (output) then goto 300;
    if (iflag > 0) then goto 100;


//.....integration successfully completed................................
//...one step mode...
    iflag := -2;
    exit;

//...interval mode...
300:t := tout;
    iflag := 2;
    exit;
end;// of Rkfs45 ---------------------------------------------------------------


function zeroin(ax,bx: TFloat; f: TFloatFunction1D; tol: TFloat): TFloat;
const
  IterMax = 9999;
var
  a, b, c, d, e, eps, fa, fb, fc: TFloat;
  tol1, xm, p, q, r, s: TFloat;
  iter: TInt;
begin
  eps := MachEps;

// initialization
  a := ax;
  b := bx;
  fa := f(a);
  fb := f(b);

  c := a;
  fc := fa;
  d := b - a;
  e := d;

  for iter := 1 to IterMax do
  begin
    if (abs(fc) < abs(fb)) then
    begin
      a := b;
      b := c;
      c := a;
      fa := fb;
      fb := fc;
      fc := fa;
    end;

    // convergence test
    tol1 := 2.0*eps*abs(b) + 0.5*tol;
    xm := 0.5*(c - b);
    if (abs(xm) <= tol1) or (fb = 0.0) then break;

    // is bisection necessary
    if (abs(e) >= tol1) or (abs(fa) > abs(fb)) then
    begin // bisection
      d := xm;
      e := d;
    end else
    begin
      if (a = c) then
      begin // inverse quadratic interpolation
        q := fa/fc;
        r := fb/fc;
        s := fb/fa;
        p := s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
        q := (q - 1.0)*(r - 1.0)*(s - 1.0);
      end else
      begin // linear interpolation
        s := fb/fa;
        p := 2.0*xm*s;
        q := 1.0 - s;
      end;

      // adjust signs
      if (p  >  0.0) then q := -q;
      p := abs(p);

      // is interpolation acceptable
      if ((2.0*p) < (3.0*xm*q - abs(tol1*q))) or (p < abs(0.5*e*q)) then
      begin // bisection
        d := xm;
        e := d;
      end else
      begin
        e := d;
        d := p/q;
      end;
    end;

    // complete step
    a := b;
    fa := fb;
    if (abs(d) > tol1) then b := b + d;
    if (abs(d) <= tol1) then b := b + sign2(tol1, xm);
    fb := f(b);

    if ((fb*(fc/abs(fc))) > 0.0) then
    begin
      c := a;
      fc := fa;
      d := b - a;
      e := d;
    end;
  end;
  result := b;
end;


function FMMfmin(ax, bx: TFloat; f: TFloatFunction1D; tol: TFloat): TFloat;
const
  IterMax = 9999;
var
  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w: TFloat;
  fu,fv,fw,fx,x: TFloat;
  iter: TInt;
begin

//  c is the squared inverse of the golden ratio

      c := 0.5*(3.0 - sqrt(5.0));

//  eps is approximately the square root of the relative machine
//  precision.

//  eps := 1.00;
//  repeat
//    eps := eps/2.00;
//    tol1 := 1.0 + eps;
//  until tol1  =  1.00;
//  eps := sqrt(eps);
  eps := SqrtMachEps;

//  initialization

  a  := ax;
  b  := bx;
  v  := a + c*(b - a);
  w  := v;
  x  := v;
  e  := 0.0;
  fx := f(x);
  fv := fx;
  fw := fx;

//  main loop starts here
  for iter := 1 to IterMax do
  begin
    xm := 0.5*(a + b);
    tol1 := eps*abs(x) + tol/3.0;
    tol2 := 2.0*tol1;

//  check stopping criterion
    if (abs(x - xm) <= (tol2 - 0.5*(b - a))) then break;

// is golden-section necessary
    if (abs(e) > tol1) then
    begin //  fit parabola
      r := (x - w)*(fx - fv);
      q := (x - v)*(fx - fw);
      p := (x - v)*q - (x - w)*r;
      q := 2.00*(q - r);
      if (q  >  0.0) then p := -p;
      q :=  abs(q);
      r := e;
      e := d;
    end;

    if (abs(e) <= tol1) or (abs(p) >= abs(0.5*q*r)) or
       (p <= q*(a - x)) or (p >= q*(b - x)) then
    begin
//  a golden-section step
      if (x >= xm) then e := a - x;
      if (x < xm) then e := b - x;
      d := c*e;
    end else
    begin
//  a parabolic interpolation step
      d := p/q;
      u := x + d;

//  f must not be evaluated too close to ax or bx
      if ((u - a) < tol2) then d := sign2(tol1, xm - x);
      if ((b - u) < tol2) then d := sign2(tol1, xm - x);
    end;

//  f must not be evaluated too close to x
    if (abs(d) >= tol1) then u := x + d;
    if (abs(d) < tol1) then u := x + sign2(tol1, d);
    fu := f(u);

//  update  a, b, v, w, and x
    if (fu  <=  fx) then
    begin
      if (u >= x) then a := x;
      if (u < x) then b := x;
      v := w;
      fv := fw;
      w := x;
      fw := fx;
      x := u;
      fx := fu;
    end else
    begin
      if (u < x) then a := u;
      if (u >= x) then b := u;
      if (fu <= fw) or (w  =  x) then
      begin
        v := w;
        fv := fw;
        w := u;
        fw := fu;
      end else
      if (fu <= fv) or (v = x) or (v = w) then
      begin
        v := u;
        fv := fu;
      end;
    end;
  end;
  result := x;
end;


procedure FMMSVD(NA,M,N: integer; const A,U,V: IFArr2D;
                const W,RV1: IFArr1D; MatU,MatV: boolean; var RetCode: integer);
const
  ItnLimit = 30;

var
  i,j,k,l,i1,k1,l1,its,mn,ExitKey: integer;
  C,G,F,X,S,H,Y,Z,Scale,ANorm,GG: TFloat;

begin
  RetCode := 0;

  if NA=M then
    for i := 1 to M  do
      for j := 1 to N  do
        U[i,j] := A[i,j]
  else
    for i := 1 to M  do
      for j := 1 to N  do
        U[i,j] := A[j,i];

  G :=  0.0;
  Scale :=  0.0;
  ANorm :=  0.0;
  for i:=1 to N  do
  begin
    l:= i+1;
    RV1[i] := Scale*G;
    G:= 0.0;
    S:= 0.0;
    Scale:= 0.0;
    if i<=M then
    begin
      for k:=1 to M  do
        Scale := Scale+abs(U[k,i]);
      if Scale<>0.0 then
      begin
        for k:=i to M  do
        begin
          U[k,i] := U[k,i]/Scale;
          S := S+sqr(U[k,i]);
        end;
        F := U[i,i];
        G := -Sign2(sqrt(S),F);
        H := F*G-S;
        U[i,i] := F-G;
        if i<>N then
          for j := l to N  do
          begin
            S := 0.0;
            for k := i to M  do
              S := S+U[k,i]*U[k,j];
            F := S/H;
            for k := i to M  do
              U[k,j] := U[k,j] + F*U[k,i];
          end;// j
          for k:=i to M  do
            U[k,i] := Scale*U[k,i];
        end;// if i<>N
      end;// if i<=M

      W[i] := Scale*G;
      G := 0.0;
      S := 0.0;
      Scale := 0.0;

      if (i<=M) and (i<>N) then
      begin
        for k:=l to N  do
          Scale := Scale+abs(U[i,k]);
        if Scale <> 0.0 then
        begin
          for k:=l to N  do
          begin
            U[i,k] := U[i,k]/Scale;
            S:= S+sqr(U[i,k]);
          end;
          F := U[i,l];
          G := -Sign2(sqrt(S),F);
          H := F*G-S;
          U[i,l] := F-G;
          for k:=l to N  do
            RV1[k] := U[i,k]/H;
          if i <> M then
            for j:=l to M  do
            begin
              S := 0.0;
              for k:=l to N  do
                S := S+U[j,k]*U[i,k];
              for k:=l to N  do
                U[j,k] := U[j,k]+S*RV1[k];
            end;
          for k:=l to N  do
            U[i,k] := Scale*U[i,k];
        end;
    end;
    ANorm := Max(ANorm,abs(W[i])+abs(RV1[i]));
  end;    {  of the loop  on  i   }

  {   Accumulation of the right-hand transformations   }
  if MatV then
    for i:= N downto 1  do
    begin
      if i<>N then
      begin
        if G<>0.0 then
        begin
          for j:=l to N  do
            V[j,i] := (U[i,j]/U[i,l]) / G;
          for j:=l to N  do
          begin
            S := 0.0;
            for k:=l to N  do
              S := S+U[i,k]*V[k,j];
            for k:=l to N  do
              V[k,j] := V[k,j]+S*V[k,i];
          end;
        end;
        for j:=l to N  do
        begin
          V[i,j] := 0.0;
          V[j,i] := 0.0;
        end;
      end;

      V[i,i] := 1.0;
      G:= RV1[i];
      l:= i;

    end; // i
  {   Accumulation of the left-hand transformations   }
  if  MatU then
  begin
    mn :=  min(M,N);
    for i:=mn downto 1  do
    begin
      l := i+1;
      G := W[i];
      if i<>N then
        for j:=l to N  do
          U[i,j] := 0.0;
      if G<>0.0 then
      begin
        if i<>mn then
          for j:=l to N  do
          begin
            S := 0.0;
            for k:=l to M  do
              S := S+U[k,i]*U[k,j];
            F := (S/U[i,i]) / G;
            for k:=i to M  do
              U[k,j] := U[k,j]+F*U[k,i];
          end;
        for j:=i to M  do
          U[j,i] := U[j,i]/G;
      end else
        for j:=i to M  do
          U[j,i] := 0.0;

      U[i,i] := U[i,i]+1.0;
    end;
  end;
  { Diagonalization of the two-diagonal form. }
  for k := N downto 1  do
  begin
    k1:= k-1;
    its := 0;

    repeat
      ExitKey:= 0;
      l:= k+1;
      while (ExitKey=0) and (l>1)  do
      begin
        dec(l);
        l1 := l-1;
        if  abs(RV1[l])+ANorm = ANorm then ExitKey:=1
        else if  l1>0 then
          if abs(W[l1])+ANorm = ANorm then ExitKey:=2;
      end;

      if ExitKey<>1 then
      begin
        C := 0.0;
        S := 1.0;
        ExitKey := 0;
        i := l;
        while (ExitKey=0) and (i<=k)  do
        begin
          F :=  S*RV1[i];
          RV1[i] := C*RV1[i];
          if abs(F)+ANorm = ANorm then ExitKey := 1 else
          begin
            G := W[i];
            H := sqroot2(F,G);
            W[i] := H;
            C := G/H;
            S := -F/H;
            if  MatU then
              for j:=1 to M  do
              begin
                Y :=  U[j,l1];
                Z :=  U[j,i];
                U[j,l1]:=  Y*C+Z*S;
                U[j,i]:=  -Y*S+Z*C
              end;
            inc(i);
          end;
        end;
      end;

      { Convergence Checking }
      Z := W[k];
      if l<>k then
        begin
          if its=Itnlimit then
          begin
            RetCode := k;
            Exit;
          end;
          inc(its);
          X:=  W[l];
          Y:=  W[k1];
          G:=  RV1[k1];
          H:=  RV1[k];
          F:=  ((Y-Z)*(Y+Z) + (G-H)*(G+H)) / (2.0*H*Y);
          if  abs(F)<=1.0 then   GG := Sign2(sqrt(F*F+1.0),F)
                           else   GG := F*sqrt(1.0+1.0/F/F);
          F:=  ((X-Z)*(X+Z) + H*(Y/(F+GG)-H)) / X;

     { Next  QR - Transformation }
          C:=  1.0;
          S:=  1.0;
          for i1:=l to k1  do
          begin
            i := i1+1;
            G := RV1[i];
            Y := W[i];
            H := S*G;
            G := C*G;
            Z := sqroot2(F,H);
            RV1[i1] := Z;
            C := F/Z;
            S := H/Z;
            F := X*C+G*S;
            G := -X*S+G*C;
            H := Y*S;
            Y := Y*C;
            if  MatV then
              for j:=1 to N  do
              begin
                X :=  V[j,i1];
                Z :=  V[j,i];
                V[j,i1]:=  X*C+Z*S;
                V[j,i]:=  -X*S+Z*C;
              end;

            Z := sqroot2(F,H);
            W[i1] := Z;
            if Z<>0.0 then
            begin
              C := F/Z;
              S := H/Z;
            end;
            F := C*G+S*Y;
            X := -S*G+C*Y;
            if  MatU then
              for J:=1 to M  do
              begin
                Y :=  U[j,i1];
                Z :=  U[j,i];
                U[j,i1]:=  Y*C+Z*S;
                U[j,i]:=  -Y*S+Z*C;
              end;
          end;

          RV1[l] := 0.0;
          RV1[k] := F;
          W[k]:= X;
        end

        else if Z<0.0 then
        begin
          W[k] := -Z;
          if  MatV then
            for j:=1 to N  do
              V[j,k] := -V[j,k];
        end;
    until  l=k;
  end; { of kk - loop }

  if NA = N then
    for i:=1 to M  do
      for j:=1 to N  do
      begin
        G := V[i,j];
        V[i,j] := U[i,j];
        U[i,j] := G;
      end;
end; { of FMMSVD --------------------------------------------------------------}


initialization

end.

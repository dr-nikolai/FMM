{@abstract(Forsythe-Malcolm-Moler collection tests)
@author(Nikolai Shokhirev <nikolai@shokhirev.com> <nikolai@u.arizona.edu>
                          http://www.shokhirev.com/nikolai.html)
@created(2004.12.12)
@lastmod(2004.12.12)
©Nikolai V. Shokhirev, 2003-2004}
unit Test_uFMM;

interface

uses
  TestFramework,
  SysUtils,
  uMatTypes,
//  uOptimization,
//  uOptimUtils, uZeroin,
//  uLevenbergMarquardt,
  uDynArrays;

type

Check_DecompSolve = class(TTestCase)
private
  V1, V2: IFArr1D;
  M1: IFArr2D;
  d: integer;
  msg: string;
public
  procedure setUp;  override;
  procedure tearDown; override;
published
  procedure VerifyDecompSolve;
end;

Check_RKF45 = class(TTestCase)
private
public
   procedure setUp;  override;
   procedure tearDown; override;
published
   procedure VerifyOrbit;
end;

Check_Spline = class(TTestCase)
private
  x, y, b, c, d: IFArr1D;
  n: TInt;
  msg: string;
public
   procedure setUp;  override;
   procedure tearDown; override;
published
   procedure VerifyFMMSpline;
   procedure VerifyNaturalSpline;
end;

Check_Quanc8 = class(TTestCase)
private
  x1, x2, abserr, relerr, eps: TFloat;
  r0, result, errest, flag: TFloat;
  nofun: TInt;
public
   procedure setUp;  override;
   procedure tearDown; override;
published
   procedure VerifyIntSin;
   procedure VerifyPoly08;
   procedure VerifySqrt;
end;

Check_ZeroMin = class(TTestCase)
private
  x0, xmin, xmax, tol: TFloat;
  x1, x2: TFloat;
  msg: string;
public
   procedure setUp;  override;
   procedure tearDown; override;
published
   procedure VerifyFmin;
   procedure VerifyZeroin;
end;

Check_SVD = class(TTestCase)
private
  W, T: IFArr1D;
  A, U, V, A1: IFArr2D;
  err, N1, N2: integer;
  del: TFloat;
  msg: string;
public
   procedure setUp;  override;
   procedure tearDown; override;
published
  procedure VerifyFMMSVD1;
  procedure VerifyFMMSVD2;
end;

function Suite : ITestSuite;

implementation

uses
  uDynArrUtils,
  uDynArrUtilsX,
  uDynObjAlg,
  uDynArrIO,
  uFMM,
  math;

var
  nexp: TInt;
  alfasqD: TFloat;

function Suite : ITestSuite;
begin
  result := TTestSuite.Create('uFMM Tests');

  result.addTest(Check_DecompSolve.Suite);
  result.addTest(Check_Spline.Suite);
  result.addTest(Check_Quanc8.Suite);
  result.addTest(Check_RKF45.Suite);
  result.addTest(Check_ZeroMin.Suite);
  result.addTest(Check_SVD.Suite);
end;

{ Check_Quanc8 }

procedure Check_Quanc8.setUp;
begin

end;

procedure Check_Quanc8.tearDown;
begin

end;

function F5(const x: TFloat): TFloat;
begin
  if x = 0.0 then
    result := 1.0
  else
    result := sin(x)/x;
end;

procedure Check_Quanc8.VerifyIntSin;
begin
  x1 := 0.0;
  x2 := 2.0*Pi;
  r0 := 1.4181515761326282;  // from Mathematica
  abserr := 2.0e-12;
  relerr := 1.0e-12;
  eps := max(abserr, relerr*r0);
  Quanc8(F5, x1, x2, abserr, relerr, result, errest, flag, nofun);
  Check(SameValue(r0,result,eps),'VerifyIntSin faild; flag = '+FloatToStr(flag));
end;


function FPn(const x: TFloat): TFloat;
begin
    result := (nexp+1.0)*IntPower(x,nexp);;
end;

procedure Check_Quanc8.VerifyPoly08;
var
  n: TInt;
begin
  x1 := 0.0;
  x2 := 1.0;
  r0 := 1.0;
  abserr := 1.0e-12;
  relerr := 1.0e-12;
  eps := max(abserr, relerr*r0);
  for n := 0 to 8 do
  begin
    nexp := n;
    Quanc8(FPn, x1, x2, abserr, relerr, result, errest, flag, nofun);
    Check(SameValue(r0, result, eps), 'VerifyPoly08 faild n = '+IntToStr(n));
  end;
end;

function FSqrt(const x: TFloat): TFloat;
begin
  result := sqrt(x);
end;

procedure Check_Quanc8.VerifySqrt;
begin
  x1 := 0.0;
  x2 := 4.0;
  r0 := 16.0/3.0;
  abserr := 1.0e-12;
  relerr := 1.0e-12;
  eps := max(abserr, relerr*r0);
  Quanc8(FSqrt, x1, x2, abserr, relerr, result, errest, flag, nofun);
  Check(SameValue(r0,result,eps),'VerifySqrt faild; flag = '+FloatToStr(flag));
end;

{ Check_SVD }

procedure Check_SVD.setUp;
begin
  inherited;

end;

procedure Check_SVD.tearDown;
begin
  inherited;

end;

procedure Check_SVD.VerifyFMMSVD1;
var
  NA: TInt;
  x : TFloat;
  MatU, MatV: boolean;
begin
  N1 := 6;  N2 := 4; // M = N1   N = N2
  x := N1;  x := x*x*x;
  A  := TFArr2D.Create(N1,N2); // M > N
  A1 := TFArr2D.Create(N1,N2);
  U  := TFArr2D.Create(N1,N2);
  V  := TFArr2D.Create(N2,N2);
  W  := TFArr1D.Create(N2);
  T  := TFArr1D.Create(N2);
  RandArr2D(A,-2.0,2.0,0);
  A1.Assign(A);
  NA := N1;
  MatU := true;
  MatV := true;
  FMMSVD(NA, N1, N2, A, U, V, W, T, MatU, MatV, err);
  Check((err=0),'Fail FMMSVD1 - FMMSVD:  err = '+IntToStr(err));
  if err <> 0 then exit;

  A := Mt1xDxMt2T(U,W,V);
  Check((err=0),'Fail FMMSVD1 - Mt1xDxMt2T:  err = '+IntToStr(err));
  if err <> 0 then exit;

  Check(SameArr(A, A1, msg, x*MachEps) and (err=0),'Fail FMMSVD1:'+msg);
end;

procedure Check_SVD.VerifyFMMSVD2;
var
  i, j, NA: TInt;
  x : TFloat;
  MatU, MatV: boolean;
begin
  N1 := 4;  N2 := 6;
  x := N2;  x := x*x*x;
  A  := TFArr2D.Create(N1,N2); // M < N
  A1 := TFArr2D.Create(N1,N2);
  U  := TFArr2D.Create(N2,N1);
  V  := TFArr2D.Create(N2,N1);
  W  := TFArr1D.Create(N1);
  T  := TFArr1D.Create(N1);
  RandArr2D(A,-2.0,2.0,0);
  A1.Assign(A);
  NA := N1;
  MatU := true;
  MatV := true;
  FMMSVD(NA, N2, N1, A, U, V, W, T, MatU, MatV, err);
  Check((err=0),'Fail FMMSVD1 - FMMSVD:  err = '+IntToStr(err));
  if err <> 0 then exit;

  A := Mt1xDxMt2T(U,W,V);
  Check((err=0),'Fail FMMSVD1 - Mt1xDxMt2T:  err = '+IntToStr(err));
  if err <> 0 then exit;
  del := 0.0;
  for i := 1 to N1 do
    for j := 1 to N2 do
      del := max(del, abs(A[i,j]-A1[i,j]));
  Check((del < x*MachEps) and (err=0),'Fail FMMSVD1: del = '+FloatToStr(del)+'  err = '+IntToStr(err));
end;

{ Check_RKF45 }

procedure Check_RKF45.setUp;
begin
  inherited;

end;

procedure Check_RKF45.tearDown;
begin
  inherited;

end;

function f3(const x: TFloat):TFloat;
begin
  result := 1.0/(x+1.0)+1/(1.0-x);
end;


// Sample6
procedure Orbit(const t: TFloat; const y, yp: IFArr1D);                            
var
  r: TFloat;
begin
//----------------------------------------------------------------------------
  r:=y[1]*y[1]+y[2]*y[2];
  r:=r*SQRT(r)/alfasqD; // alfasqD is a module variable, used here and Sample6
  yp[1]:=y[3];
  yp[2]:=y[4];
  yp[3]:=-y[1]/r;
  yp[4]:=-y[2]/r;
  exit ;
end;//

// Sample6  from PDAS: http://www.pdas.com/free.htm
procedure Check_RKF45.VerifyOrbit;
var
  i, NEQN, iflag, nfe, kop, init, jflag, kflag : TInt;
  y, yp, f1, f2, f3, f4, f5: IFArr1D;
  y1, y2: IFArr1D;
  t,  tout,  tprint,  tfinal,  alfa,  ecc,  relerr,  abserr: TFloat;
  err, h, savre, savae: TFloat;
  ss: string;
begin
  NEQN := 4;
  y  := TFArr1D.Create(NEQN);
  yp := TFArr1D.Create(NEQN);
  f1 := TFArr1D.Create(NEQN);
  f2 := TFArr1D.Create(NEQN);
  f3 := TFArr1D.Create(NEQN);
  f4 := TFArr1D.Create(NEQN);
  f5 := TFArr1D.Create(NEQN);

  y1 := TFArr1D.Create(25);
  ss := '(0.750000000,0.619766772,0.294413418,-0.105183668,-0.490308911,'+
        '-0.813951671,-1.054037690,-1.200736640,-1.249995470,-1.200723410,'+
        '-1.054012660,-0.813917398,-0.490270019,-0.105147965,0.294436693,'+
        '0.619762182,0.749918401,0.619556844,0.294083089,-0.105574392,'+
        '-0.490691930,-0.814271748,-1.054254890,-1.200823430,-1.249934080)';
  StrToVector(y1, ss);
  y2 := TFArr1D.Create(25);
  ss := '(0.000000000,0.477790743,0.812175512,0.958031714,0.939863801,'+
        '0.799574375,0.575683534,0.300133795,-0.000028848,-0.300188512,'+
        '-0.575729311,-0.799605250,-0.939875066,-0.958015144,-0.812122822,'+
        '-0.477689207,0.000150567,0.477936178,0.812232316,0.957946718,'+
        '0.939628422,0.799207032,0.575215459,0.299602211,-0.000583312)';
  StrToVector(y2, ss);
  iflag := 1;
  tprint := 0.5;
  tfinal := 12.0;
  alfa := 3.141592653589/4.0;
  relerr := 1.0E-5;
  abserr := 1.0E-5;

  ecc:=0.25;
  alfasqD:=alfa*alfa;
  y[1]:=1.0-ecc;
  y[2]:=0.0;
  y[3]:=0.0;
  y[4]:=alfa*SQRT((1.0+ecc)/(1.0-ecc));
  t:=0.0;
  tout:=t;
  ss := '';
  err := 0.0;
  i := 0;
  repeat
    inc(i);
    Rkfs45(Orbit,t,tout, abserr, relerr, iflag, y, yp, f1, f2, f3, f4, f5,
           h, savre, savae, nfe, kop, init, jflag, kflag);
    err := max(err,abs(y[1]-y1[i])+abs(y[2]-y2[i]));
    CASE (iflag) of
    1: begin
         ss := ss + ' IMPROPER CALL';
         break;
       end;
    3: ss := ss + ' Tolerances reset due to iflag = 3, relerr = '+
                   FloatToStr(relerr)+' abserr = '+ FloatToStr(abserr);
    4: ss := ss + 'Warning.  Many steps...';
    5: begin
         abserr := 1.0E-9;
         ss := ss + ' Tolerances reset due to iflag = 5, relerr = '+
                     FloatToStr(relerr)+' abserr = '+ FloatToStr(abserr);
       end;
    6: begin
         relerr := 10.0*relerr;
         ss := ss + ' Tolerances reset due to iflag = 6, relerr = '+
                     FloatToStr(relerr)+' abserr = '+ FloatToStr(abserr);
       end;
    7: begin
         ss := ss + ' WARNING. Much output';
         iflag := 2;
       end;
    else
        tout := t+ tprint;
    end;
  until (t >= tfinal);
  ss := ss + ' err = '+FloatToStr(err);
  ss := ss + ' Successful completion of Sample6';
  Check((iflag=2) and (err < 2.0E-5), 'faild VerifyOrbit ' + ss);
end;


{ Check_Spline }

procedure Check_Spline.setUp;
begin
  inherited;
  n := 10;
  x := TFArr1D.Create(2,n+1);
  y := TFArr1D.Create(2,n+1);
  b := TFArr1D.Create(2,n+1);
  c := TFArr1D.Create(2,n+1);
  d := TFArr1D.Create(2,n+1);
end;

procedure Check_Spline.tearDown;
begin
  inherited;

end;

procedure Check_Spline.VerifyNaturalSpline;
var
  i, k : TInt;
  t, u, f, fp, fpp, fppp: TFloat;
begin
  t := 1;
  for i := x.Lo1 to x.Hi1 do
  begin
    x[i] := t;
    y[i] := t*t*t;
    t := t+1.0;
  end;
  NaturalSpline(x, y, b, c, d);
  i := x.Lo1; // 1st time
  t := 0.0;
  for k := x.Lo1 to x.Hi1 do
  begin
    u := x[k];
    f := Seval(u, i, x, y, b, c, d);
    t := max(t,abs(f-y[k]));
  end;
  msg :=  'faild VerifyNaturalSpline Seval at knots err='+FloatToStr(t);
  Check((t < 1.0e-11), msg);
  //Chek only at knots because NaturalSpline does not interpolate!
end;

procedure Check_Spline.VerifyFMMSpline;
var
  i : TInt;
  t, u, f, fp, fpp, fppp: TFloat;
begin
  t := 1;
  for i := x.Lo1 to x.Hi1 do
  begin
    x[i] := t;
    y[i] := t*t*t;
    t := t+1.0;
  end;
  FMMSpline(x, y, b, c, d);
  u := 2.5;
  i := x.Lo1; // 1st time
  f := Seval(u, i, x, y, b, c, d);
  t :=abs(f-15.625);
  msg :=  'faild VerifyFMMSpline Seval err='+FloatToStr(t);
  Check((t < 1.0e-11), msg);
  i := x.Lo1; // 1st time
  Seval3(u, i, x, y, b, c, d, f, fp, fpp, fppp);
  t :=abs(f-15.625);
  msg :=  'faild VerifyFMMSpline Seval3 f err='+FloatToStr(t);
  Check((t < 1.0e-11), msg);
  t :=abs(fp-18.75);
  msg :=  'faild VerifyFMMSpline Seval3 fp err='+FloatToStr(t);
  Check((t < 1.0e-11), msg);
  t :=abs(fpp-15.0);
  msg :=  'faild VerifyFMMSpline Seval3 fpp err='+FloatToStr(t);
  Check((t < 1.0e-11), msg);
  t :=abs(fppp-6.0);
  msg :=  'faild VerifyFMMSpline Seval3 fppp err='+FloatToStr(t);
  Check((t < 1.0e-11), msg);
end;


{ Check_DecompSolve }

procedure Check_DecompSolve.setUp;
begin
  inherited;

end;

procedure Check_DecompSolve.tearDown;
begin
  inherited;

end;

procedure Check_DecompSolve.VerifyDecompSolve;
var
  x : TFloat;
  cond: TFloat;
  ipvt: IIArr1D;
begin
  d := 6;
  x := d; x := x*x*x;
  M1 := TFArr2D.Create(d,d);
  ipvt := TIArr1D.Create(d);
  RandArr2D(M1,-1.0,1.0,0);
  MtShiftDiag(M1, 3.0);
  V1 := TFArr1D.Create(d);
  RandArr1D(V1,-1.0,1.0);
  V2 := MtxVt(M1,V1);
  FMMDecomp(M1, cond, ipvt);
  FMMSolve(M1, V2, ipvt);
   Check(SameArr(V2, V1, msg, x*MachEps),'Fail VerifyDecompSolve:'+msg);
end;

{ Check_ZeroMin }

procedure Check_ZeroMin.setUp;
begin
  inherited;

end;

procedure Check_ZeroMin.tearDown;
begin
  inherited;

end;

function f(const x: TFloat): TFloat;
begin
  result := exp(2.0-x-x)-2.0*exp(1.0-x);
end;

procedure Check_ZeroMin.VerifyFmin;
var
  i : TInt;
begin
  xmin := 0.0;
  xmax := 2.0;
  tol := 0.000001;
  for i := 1 to 99 do
  begin
    x1 := RndFloat(xmin,0.8);
    x2 := RndFloat(1.2,xmax);
    x0 := FMMfmin(x1, x2, f, tol);
    msg := 'faild VerifyFmin:'+', x0='+FloatToStr(x0)+', i='+IntToStr(i)+
           ', x1='+FloatToStr(x1)+', x2='+FloatToStr(x2);
    Check(abs(x0-1.0)<tol, msg);
  end;
end;

function dfdx(const x: TFloat): TFloat;
begin
  result := -2.0*(exp(2.0-x-x)-exp(1.0-x));
end;

procedure Check_ZeroMin.VerifyZeroin;
var
  i : TInt;
begin
  xmin := 0.0;
  xmax := 2.0;
  tol := 0.000001;
  for i := 1 to 99 do
  begin
    x1 := RndFloat(xmin,0.8);
    x2 := RndFloat(1.2,xmax);
    x0 := zeroin(x1, x2, dfdx, tol);
    msg := 'faild VerifyZeroin:'+', x0='+FloatToStr(x0)+', i='+IntToStr(i)+
           ', x1='+FloatToStr(x1)+', x2='+FloatToStr(x2);
    Check(abs(x0-1.0)<tol, msg);
  end;
end;

initialization

  RegisterTest('uFMM Tests',Suite);

end.

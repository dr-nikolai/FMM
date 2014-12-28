{@abstract(Math types and constants.
 Basic Math Types, Interfaces, Constants and Math utilities)
@author(Nikolai Shokhirev <nikolai@shokhirev.com> <nikolai@u.arizona.edu>
        http://www.shokhirev.com/nikolai.html
        http://www.chem.arizona.edu/~shokhirn/nikolai.html)
@created(January 01, 2003)
@lastmod(December 12, 2004)
©Nikolai V. Shokhirev, 2001-2003 }
{
Januaty 2004:
Incorporated some definitions, names and functions from Earl F. Glynn's 
ComplexMathlibrary unit (www.efg2.com/Lab/Mathematics/ComplexMath.htm)
They alre labeled "From EFG Unit" below.
NOTE: the type Complex differs from TComplex in EFG Unit}
unit  uMatTypes;

interface

uses
  SysUtils;

type
  TArrayType =(atGeneric,    // [Lo..Hi] - arbitrary limits
               atNatural,    // [1..Dim] - one-based
               atZeroBased,  // [0..Dim-1] - zero-based
               atCentered);  // [-((Dim-1) div 2)..Dim div 2]

const
  { for double = 2.22044604925031e-16}
  cMachEps = 2.22044604925031e-16;
  { for double = 1.49011611938476e-8}
  cSqrtMachEps = 2.22044604925031e-16;
  { for double = 1.7e-308}
  MinFloat = 1.7e-308;
  { for double = 1.7e+308 }
  MaxFloat = 1.7e+308;
  { for double = -1.7e+308}
  NegMaxFloat = -1.7e+308;
  { for double = 1.3e-154}
  SqrtMinFloat = 1.3e-154;
  { for double = 1.3e+154}
  SqrtMaxFloat = 1.3e+154;
  { for double = 308.0}
  MaxExp =  308.0;
  { for double = 709.0 ~ 308*ln(10) }
  ExpArg = 709.0;
  { for double = 1.1e-77}
  SqrtSqrtMinFloat = 1.1e-77;
(*
  { for Extended = 3.4e-4932 ; actual 3.6 x 10^–4951 for Extended}
  MinFloat = 3.4e-4932;
  { for Extended = 1.1e+4932 ; actual 1.1 x 10^4932 for Extended}
  MaxFloat = 1.1e+4932;
  NegMaxFloat = -1.1e+4932;
  { for Extended = 1.9e-2466}
  SqrtMinFloat = 1.9e-2466;
  { for Extended = 1.0e+2466}
  SqrtMaxFloat = 1.0e+2466;
  { for Extended = 4932.3}
  MaxExp =  4932.3;
  { for Extended = 11356.3 = 4932*ln(10) }
  ExpArg = 11356.3;   // 4932*ln(10)
*)
  { HighInt = 2147483647}
  HighInt = 2147483647;
  { LowInt = -2147483647}
  LowInt = -2147483647;
  { Safe Factor = 1024.0}
  SafeFactor = 1024.0;

  // TwoPI = 2.0*PI;  [From EFG Unit]
  TwoPI = 6.283185307179586; // 6.283185307179586476925286766559;
  // HalfLn2PI = 0.5*(Ln(TwoPI));  [From EFG Unit]
  HalfLn2PI = 0.9189385332046727;  // 0.91893853320467274178032973640562
type

  { Interface for the Restore objects }
  IRestore = interface
  ['{DF81E0EC-94AC-4DBA-933F-2471667C2796}']
  end;

  { Interface for Comment: string }
  IComment = interface
  ['{C9CB72E4-C14C-466A-B485-36B8C78553E7}']
    function GetComment: string;
    procedure SetComment(const Value: string);
    property Comment: string read GetComment write SetComment;
  end;

  { Object Comment: string }
  TComment = class(TInterfacedObject, IComment)
  protected
    fComment: string;
    function GetComment: string;
    procedure SetComment(const Value: string);
  public
    property Comment: string read GetComment write SetComment;
  end;

  { currently = longint }
  TInt   = longint;
  { currently = double }
  TFloat  =  double;

  { complex number }
  Complex = record
              Re,Im: TFloat;
            end;

  { Float function of x }
  TFloatFunction1D = function(const x: TFloat): TFloat;
  { Float function of x, y }
  TFloatFunction2D = function(x, y: TFloat): TFloat;
  { Float function of x, y, z}
  TFloatFunction3D = function(x, y, z: TFloat): TFloat;
  { Float function of n, x}
  TBasisFunction = function(n: TInt; x: TFloat): TFloat;

  TComplexfunction = function (const z: Complex): Complex; // From EFG Unit

  EComplexLnZero     = class(exception);  // From EFG Unit
  EComplexZerodivide = class(exception);  // From EFG Unit
  EComplexinvalidop  = class(exception);  // From EFG Unit
  EComplexZerotoZero = class(exception);  // From EFG Unit

  { Conversion functions for Complex }

  { Basic Conversion functions for Complex }
  function ComplexToStr(const z: Complex): string;
  { Conversion: zero Im is not diaplayed }
  function CmplxToStr(const z: Complex; Width: TInt; Decimals: TInt): string;
  { Conversion: zero Im is diaplayed }
  function CmplxToStr0(const z: Complex; Width: TInt; Decimals: TInt): string;
  { Conversion }
  function StrToCmplx(const s: string): Complex;

  { Polar To Complex Conversion: (r, angle) -> (Re, Im) }
  function PolarToComplex(r, angle: TFloat): Complex;
  { Complex To Polar Conversion: (Re, Im) -> (r, angle) }
  procedure ComplexToPolar(const z: Complex; var r, angle: TFloat);

  { Basic functions for Complex Calculations }

  { result = max(|z1.re-z2.re|,|z1.im-z2.im|) }
  function Diff0(const z1, z2: Complex): TFloat;
  { abs(c1.Re-c2.Re) < eps And abs(c1.Im-c2.Im) < eps }
  function SameComplex(const c1, c2: Complex; eps: TFloat): boolean;
  { result = x + i*y}
  function cmplx(x: TFloat; y: TFloat = 0.0): Complex;
  { result = x + i*y - same as cmplx [From EFG Unit] }
  function Cset (x: TFloat; y: TFloat = 0.0): Complex;
  { result = sqrt(sqr(z.Re)+sqr(x.Im)) = |z|  }
  function Cabs (const z: Complex): TFloat;
  { max(abs(z.re),abs(z.im)) }
  function Cabs0(const z: Complex): TFloat;
  { result = abs(z.Re) + abs(z.Im)    }
  function Cabs1(const z: Complex): TFloat;
  { result = sqr(z.Re)+ sqr(z.Im) = |z|^2  }
  function Cabs2(const z: Complex): TFloat;
  { result = sqr(z.Re)+ sqr(z.Im) = |z|^2 same as Cabs2 [From EFG Unit] }
  function CAbsSqr (const z: Complex): TFloat;
  { result = -z }
  function Cneg(const z: Complex): Complex;
  { result = 1/z }
  function Cinv(const z: Complex): Complex;
  { result = z.Re -i*z.Im }
  function conjug(const z: Complex): Complex;
  { result = z.Re -i*z.Im same as conjug [From EFG Unit]}
  function Cconjugate(const z: Complex): Complex;
  { result = z.Re }
  function Re(const z: Complex): TFloat;
  { result = z.Im }
  function Im(const z: Complex): TFloat;
  { result = z/r r := Cabs (z) }
  function Cunit(const z: Complex; var r: TFloat): Complex;
  { result = z + w }
  function Cadd(const z: Complex;  w: Complex): Complex; overload;
  { result = z + x }
  function Cadd(const z: Complex;  x: TFloat): Complex; overload;
  { result = x + w }
  function Cadd(x: TFloat; const w: Complex): Complex; overload;
  { result = z - w }
  function Csub(const z, w: Complex): Complex; overload;
  { result = z - x }
  function Csub(const z: Complex;  x: TFloat): Complex; overload;
  { result = x - w }
  function Csub(x: TFloat; const w: Complex): Complex; overload;
  { result = z * w }
  function Cmul(const z, w: Complex): Complex; overload;
  { result = x * w }
  function Cmul(x: TFloat; const w: Complex): Complex; overload;
  { result = z * x }
  function Cmul(const z: Complex;  x: TFloat): Complex; overload;
  { result = z / w }
  function Cdiv(const z, w: Complex): Complex; overload;
  { result = x / w }
  function Cdiv(x: TFloat; const w: Complex): Complex; overload;
  { result = z / x }
  function Cdiv(const z: Complex;  x: TFloat): Complex; overload;
  { result = exp(v)  }
  function Cexp(const z: Complex): Complex;
  { result = exp(i*f)  }
  function CexpIm(f: TFloat): Complex;
  { result = sqrt(z), result.Re > 0 }
  function Csqrt(const z: Complex): Complex; overload;
  { result = sqrt(z), result.Re > 0 }
  function Csqrt(const x: TFloat): Complex; overload;
  // result = SQR(a)  [From EFG Unit]
  function CSqr  (const z: Complex): Complex;

  { result = sqrt(x) }
  //function Croot(x: TFloat): Complex;

  // -PI < theta <= PI     [From EFG Unit]
  function FixAngle (const theta: TFloat): TFloat;
  { complex natural log: result = ln(a)
    NOTE: principal value only [From EFG Unit] }
  function CLn (const z: Complex): Complex;

  { s := s + z }
  procedure Csum(var s: Complex; const z: Complex); overload;
  { s := s + x }
  procedure Csum(var s: Complex; x: TFloat); overload;

  { r1 <= re < r2, i1 <= im < i2 }
  function RndComplex(r1, r2, i1, i2: TFloat): Complex;

  { Special functions for Complex Calculations }

  // complex Bessel functions of order zero [From EFG Unit]
//  function CI0 (const z: Complex): Complex;        // z = I0(a)
  // complex Bessel functions of order zero [From EFG Unit]
//   function CJ0 (const z: Complex): Complex;        // z = J0(a)

//  function CLnGamma (const a: Complex): Complex;
//  function CGamma   (const a: Complex): Complex;
  { result = er(z) }
  function cer(const z: Complex): Complex;

const
  // cmplx(0.0, 0.0)
  Cmplx0: Complex = (Re: 0.0; Im: 0.0);
  // cmplx(1.0, 0.0)
  Cmplx1: Complex = (Re: 1.0; Im: 0.0);
  // cmplx(0.0, 1.0)
  CmplxIm1: Complex = (Re: 0.0; Im: 1.0);

  // cmplx(0.0, 0.0) [From EFG Unit]
  ComplexZero: Complex = (Re: 0.0; Im: 0.0);
  // cmplx(1.0, 0.0)  [From EFG Unit]
  ComplexOne: Complex = (Re: 1.0; Im: 0.0);

var
  MachEps, SqrtMachEps: TFloat;

  { boolean to string conversion }
  function BoolToStr(Value: boolean): string;
  { string to boolean conversion }
  function StrToBool(Value: string): boolean;

  { Calculation of the machine epsilon  }
  procedure MachinEps ( var MachEps  : TFloat );
  { Swap of float numbers }
  procedure Swap(var f1, f2: TFloat);
  { Swap of integer numbers }
  procedure iSwap(var i1, i2: integer);
  { random true/false }
  function RndBool: boolean;
  { random number r1 <= r < r2 }
  function RndFloat(r1, r2: TFloat): TFloat;
  { random number N1 <= result <= N2 }
  function RndNum(N1, N2: integer): Integer;
  { random number N1 <= r <= N2, r <> r1 }
  function RndNumEx(N1, N2, r1: integer): Integer;
  { if x > 0 then result := 1 else if x = 0 then result := 0 else result := -1;}
  function iSign(x: TInt)  : TInt;
  { Kroneker delta : if x1 = x2 then result := 1 else result := 0 }
  function delta(x1: TInt; x2: TInt = 0): TInt;
  { if x > 0.0 then result := 1.0 <br>
    else if x < 0.0 then result := -1.0 else result := 0.0; }
  function Sign0(x: TFloat): TFloat;
  {   if  X >= 0.0  then  result := 1.0 else  result := -1.0;}
  function Sign1(x: TFloat) : TFloat;
  { FORTRAN signum: |a|*signum(b) }
  function Sign2( a, b : TFloat) : TFloat;
  // FORTRAN signum
  function ISign2( a, b : TInt) : TInt;
  { n!/(m!(n-m)! }
  function Combinations(m, n: TInt): TInt;
  { n! }
  function Fac(m: TInt): TInt;
  { m**n }
  function IPower(m, n: TInt): TInt;
  { obsolete exp1 = exp}
  function exp1(X: TFloat): TFloat;
  { 10**x or 10^x = Power(10.0, x) [Math]}
  function power10(x: TFloat)   : TFloat;
  { lg(x) = log10(x) [Math] }
  function lg(x: TFloat)        : TFloat;
  { FORTRAN amin2(x1, x2) = min(x1, x2) [Math] }
  function amin2(x1, x2: TFloat): TFloat;
  { FORTRAN amax2(x1, x2) = max(x1, x2) [Math] }
  function amax2(x1, x2: TFloat): TFloat;
  { safe sqrt(x + y); it is assumed that (x+y) > 0 }
  function sqroot1(x, y: TFloat): TFloat;
  { safe sqrt(x*x+y*y); functionally the same as pythag}
  function sqroot2(x, y: TFloat): TFloat;
  { finds dsqrt(a**2+b**2) without overflow or destructive underflow;
    functionally the same as sqroot2}
  function pythag(a, b: TFloat): TFloat;

implementation { ============================================================= }

uses
  Math;

{ TComment }

function TComment.GetComment: string;
begin
  result := fComment;
end;

procedure TComment.SetComment(const Value: string);
begin
  fComment := Value;
end;

{ Unit Functions }

function BoolToStr(Value: boolean): string;
begin
  if Value then result := 'True' else result := 'False';
end;

function StrToBool(Value: string): boolean;
begin
  result := UpperCase(Value) = 'TRUE' ;
end;

procedure MachinEps(var MachEps: TFloat);
var
  x: TFloat;
begin
  MachEps := 1.0;
  repeat
    MachEps := MachEps/2.0;
    x := 1.0 + MachEps;
  until not (x > 1.0);
  MachEps := MachEps*2.0;
// MachEps = 2.2204460492e-16 for double
// MachEps = 1.0842021725e-19 for extended
end;

procedure MachEpsMinus(var MachEps  : TFloat );
begin
  MachEps := 1.0;
  while (1.0-MachEps/2.0) < 1.0  do
    MachEps := MachEps/2.0;
// MachEpsMinus = 1.11022302462516E-16 for double
end;

{ Swap of float numbers }
procedure Swap(var f1, f2: TFloat);
var
  f: TFloat;
begin
  f := f1;
  f1 := f2;
  f2 := f;
end;

{ Swap of integer numbers }
procedure iSwap(var i1, i2: integer);
var
  i: integer;
begin
  i := i1;
  i1 := i2;
  i2 := i;
end;

function RndBool: boolean;
begin
  result := (random(2) = 0);
end;

{ r1 <= r < r2 }
function RndFloat(r1, r2: TFloat): TFloat;
begin
  result := r1 + (r2-r1)*random;
end;

{ N1 <= r <= N2 }
function RndNum(N1, N2: integer): Integer;
begin
  result := N1 + random(N2-N1+1);
end;

{ N1 <= r <= N2, r <> r1 }
function RndNumEx(N1, N2, r1: integer): Integer;
begin
  repeat
    result := N1 + random(N2-N1+1);
  until result <> r1;
end;

function  iSign (x: TInt): TInt;
begin
  if x > 0 then result := 1 else
    if x = 0 then result := 0 else
      result := -1;
end;

function delta(x1: TInt; x2: TInt): TInt;
begin
  if x1 = x2 then result := 1 else result := 0;
end;

function sign0(x: TFloat): TFloat;
begin
  if x > 0.0 then result := 1.0
  else if x < 0.0 then result := -1.0
                  else result := 0.0;
end;

function  Sign1 ( X : TFloat ) : TFloat;
begin
  if  X >= 0.0  then  result := 1.0
                else  result := -1.0;
end;

// FORTRAN signum
function Sign2( a, b : TFloat) : TFloat;
begin
  if b < 0.0 then result := - abs(a)
             else result :=   abs(a) ;
end;

// FORTRAN signum
function ISign2( a, b : TInt) : TInt;
begin
  if b < 0 then result := - abs(a)
           else result :=   abs(a) ;
end;

function Combinations(m, n: TInt): TInt;
var
  i, C   : TInt;
begin
  C := 1;
  if m > 0  then
    for i := 1 to m  do
      C := (C*(n-i+1)) div i;
  Combinations := C;
end;

function Fac(m: TInt): TInt;
var
  i, z: TInt;
begin
  z := 1;
  if m > 1  then
    for i := 2 to m  do
      z := z*i;
  Fac := z;
end;

{ m**n }
function IPower(m, n: TInt): TInt;
var
  n1, n2: TInt;
begin
  if n = 0 then
  begin
    result := 1;
    exit;
  end;
  if n = 1 then
  begin
    result := m;
    exit;
  end;
  if n = 2 then
  begin
    result := m*m;
    exit;
  end;
  n1 := n div 2;
  n2 := n - n1;
  result := IPower(m, n1)*IPower(m, n2);
end;

// obsolete
function exp1(X: TFloat): TFloat;{ exp with underflow protection}
begin
  if X = 0.0 then exp1 := 1.0
             else if X > -ExpArg then  exp1 := exp(X)
                                 else  exp1 := 0.0;
end;

function power10(x: TFloat): TFloat;
begin
if x > -MaxExp then power10 := exp1(x*ln(10.0))
               else power10 := 1.0;
end;

function lg(x: TFloat): TFloat;
begin
  lg := ln(x)/ln(10.0);
end;

function amin2(x1,x2: TFloat): TFloat;
begin
if x1 < x2 then amin2 := x1
           else amin2 := x2;
end;

function amax2(x1,x2: TFloat): TFloat;
begin
if x1 > x2 then amax2 := x1
           else amax2 := x2;
end;

{ safe sqrt: sqroot = sqrt(x + y) }
function sqroot1(x, y: TFloat) : TFloat;
begin
  if (x + y) <= 0 then result := 0.0
  else if x > y then
    result := sqrt(x)*sqrt(1.0 + y/x)
  else
    result := sqrt(y)*sqrt(1.0 + x/y);
end;

{ safe sqrt: sqroot2 = sqrt(x*x+y*y); functionally the same as pythag}
function sqroot2(x, y: TFloat) : TFloat;
begin
  if abs(x) > abs(y) then
    result := abs(x)*sqrt(1.0+sqr(y/x))
  else if abs(y) = 0.0 then
    result := 0.0
  else
    result := abs(y)*sqrt(1.0+sqr(x/y));
end;

{ finds dsqrt(a**2+b**2) without overflow or destructive underflow;
  functionally the same as sqroot2}
function pythag(a, b: TFloat): TFloat;
var
  p, r: TFloat;
begin
  result := 0.0;
  p := max(abs(a),abs(b));
  if p = 0.0 then exit;
  r := min(abs(a),abs(b));
  if r < SqrtMinFloat*p then result := p
  else
    result := p*sqrt(1.0 + sqr(r/p));
{ original EISPACK algorithm:
var
  p,r,s,t,u: TFloat;
  1, 2: label;
begin
    p := max(abs(a),abs(b));
    if (p  =  0.0) then goto 2
    r := sqr(min(abs(a),abs(b))/p);
1:t := 4.0 + r;
  if (t  =  4.0) goto 2;
  s := r/t;
  u := 1.0 + 2.0*s;
  p := u*p;
  r := (s/u)**2 * r;
  goto 1;
2:result := p;
end;}
end;

{ Basic Conversion functions for Complex }
function ComplexToStr(const z: Complex): string;
begin
  result := FloatToStr(z.Re);
  if z.Im > 0.0 then
    result := result + ' +i*'+FloatToStr(z.Im)
  else
    result := result + ' -i*'+FloatToStr(z.Im);
end;

{ Conversion functions for Complex }
function CmplxToStr(const z: Complex; Width: TInt; Decimals: TInt): string;
var
  s: string;
begin
  Str(z.Re:Width:Decimals, s);
  result := s;
  if z.Im <> 0.0 then
  begin
    Str(abs(z.Im):Width:Decimals, s);
    s := Trim(s);
    if z.Im > 0.0 then
      result := result + ' +i*'+s
    else
      result := result + ' -i*'+s;
  end;
end;

{ Conversion functions for Complex }
function CmplxToStr0(const z: Complex; Width: TInt; Decimals: TInt): string;
var
  s: string;
begin
  Str(z.Re:Width:Decimals, s);
  result := s;
//  if z.Im <> 0.0 then
  begin
    Str(abs(z.Im):Width:Decimals, s);
    s := Trim(s);
    if z.Im > 0.0 then
      result := result + ' +i*'+s
    else
      result := result + ' -i*'+s;
  end;
end;

function StrToCmplx(const s: string): Complex;
//  2.2 + i*1.1
// +2.2 - I*1.1
// -2.2 + j 1.1
//  2.2 + J*1.1
var
  r, i, ss: string;
  n: TInt;
  sig: TFloat;
begin
  ss := StringReplace(s, ' ', '', [rfReplaceAll, rfIgnoreCase]);
  n := Pos('i',ss);
  if n = 0 then n := Pos('I',ss);
  if n = 0 then n := Pos('j',ss);
  if n = 0 then n := Pos('J',ss);
  if n = 0 then // real
  begin
    r := ss;
    i := '';
  end else
  begin
    case n of
    1: begin // imaginary
         sig := 1.0;
         r := '';
         if ss[n+1] = '*' then i := copy(ss, n+2,512)
         else i := copy(ss, n+1,512);
       end;
    2: begin // imaginary
         r := '';
         if ss[1] = '-' then sig := -1.0 else sig := 1.0;
         if ss[n+1] = '*' then i := copy(ss, n+2,512)
         else i := copy(ss, n+1,512);
       end;
    else
      sig := 1.0;
      if ss[n-1] = '-' then
      begin
        sig := -1.0;
        r := copy(ss, 1, n-2);
      end else
        if ss[n-1] = '+' then
          r := copy(ss, 1, n-2)
        else // no sign = +
          r := copy(ss, 1, n-1);

      if ss[n+1] = '*' then i := copy(ss, n+2,512)
      else i := copy(ss, n+1,512);
    end;
  end;

  if r ='' then
    result.Re := 0.0
  else
    result.Re := StrToFloat(r);

  if i ='' then
    result.Im := 0.0
  else
    result.Im := sig*StrToFloat(i);
end;

{ Polar To Complex Conversion: (r, angle) -> (Re, Im) }
function PolarToComplex(r, angle: TFloat): Complex;
begin
  result.Re := r*cos(angle);
  result.Im := r*sin(angle);
end;

{ Complex To Polar Conversion: (Re, Im) -> (r, angle) }
procedure ComplexToPolar(const z: Complex; var r, angle: TFloat);
begin
  r := Cabs(z);
  angle := Arctan2(z.Im, z.Re);
end;

{ Basic functions for Complex Calculations }

{ result = max(abs(z1.re-z2.re),abs(z1.im-z2.im)) }
function Diff0(const z1,z2: Complex): TFloat;
begin
  result := max(abs(z1.re-z2.re),abs(z1.im-z2.im));
end;

{ abs(c1.Re-c2.Re) < eps And abs(c1.Im-c2.Im) < eps }
function SameComplex(const c1, c2: Complex; eps: TFloat): boolean;
begin
  result := (abs(c1.Re-c2.Re) < eps) and (abs(c1.Im-c2.Im) < eps);
end;

{ result = x + i*y }
function cmplx(x: TFloat; y: TFloat = 0.0): Complex;
begin
  result.Re := x;
  result.Im := y;
end;

{ result = x + i*y }
function Cset(x: TFloat; y: TFloat = 0.0): Complex;
begin
  result.Re := x;
  result.Im := y;
end;

function Cabs (const z: Complex ): TFloat;
{ sqrt(sqr(Re)+sqr(Im)) }
begin
  result := sqroot2(z.Re, z.Im);
end;

function Cabs0(const z: Complex): TFloat;
{ max(abs(z.re),abs(z.im)) }
begin
  result := max(abs(z.re),abs(z.im));
end;

function Cabs1(const z: Complex ): TFloat;
{ abs(z.Re) + abs(z.Im)    }
begin
  result := abs(z.re) + abs(z.im);
end;

function Cabs2(const z: Complex ): TFloat;
{ sqr(z.Re)+ sqr(z.Im)     }
begin
  result := sqr(z.re) + sqr(z.im);
end;

function CAbsSqr(const z: Complex ): TFloat;
{ sqr(z.Re)+ sqr(z.Im)     }
begin
  result := sqr(z.re) + sqr(z.im);
end;

function Re(const z: Complex):TFloat;
begin
  result := z.re;
end;

function Im(const z: Complex):TFloat;
begin
  result := z.im;
end;

function Cneg(const z: Complex): Complex;
{ result = -z }
begin
  result.Re := -z.Re;
  result.Im := -z.Im;
end;

function Cinv(const z: Complex): Complex;
{ result = 1/z }
var
  r: TFloat;
begin
  if Cabs0(z) > MinFloat then
  begin
    r := Cabs2(z);
    result.Re :=  z.Re/r;
    result.Im := -z.Im/r;
  end else
  begin
    result.Re :=  sign(z.Re)*MaxFloat;
    result.Im := -sign(z.Im)*MaxFloat;
  end;
end;

function Conjug(const z: Complex): Complex;
begin
  result.Re :=  z.Re;
  result.Im := -z.Im;
end;

function Cconjugate(const z: Complex): Complex;
begin
  result.Re :=  z.Re;
  result.Im := -z.Im;
end;

{ result = z/r r := Cabs (z) }
function Cunit(const z: Complex; var r: TFloat): Complex;
begin
  r := Cabs(z);
  if r > MinFloat then
  begin
    result.Re := z.Re/r;
    result.Im := z.Im/r;
  end else
  begin
    r := 0.0;
    result.Re := 0.0;
    result.Im := 0.0;
  end;
end;

function Cadd(const z: Complex;  w: Complex): Complex;
begin
  result.Re :=  z.Re + w.Re;
  result.Im :=  z.Im + w.Im;
end;

function Cadd(const z: Complex;  x: TFloat): Complex;
begin
  result.Re :=  z.Re + x;
  result.Im :=  z.Im;
end;

function Cadd(x: TFloat; const w: Complex): Complex;
begin
  result.Re :=  w.Re + x;
  result.Im :=  w.Im;
end;

function Csub(const z, w: Complex): Complex;
begin
  result.Re :=  z.Re - w.Re;
  result.Im :=  z.Im - w.Im;
end;

function Csub(const z: Complex;  x: TFloat): Complex;
begin
  result.Re :=  z.Re - x;
  result.Im :=  z.Im;
end;

function Csub(x: TFloat; const w: Complex): Complex;
begin
  result.Re :=  x - w.Re;
  result.Im :=    - w.Im;
end;

function Cmul(const z, w: Complex): Complex;
begin
  result.Re :=  z.Re*w.Re - z.Im*w.Im;
  result.Im :=  z.Im*w.Re + z.Re*w.Im;
end;

function Cmul(x: TFloat; const w: Complex): Complex;
begin
  result.Re :=  x*w.Re;
  result.Im :=  x*w.Im;
end;

function Cmul(const z: Complex; x: TFloat): Complex;
begin
  result.Re :=  z.Re*x;
  result.Im :=  z.Im*x;
end;

function Cdiv(const z, w: Complex): Complex;
{ result = z / w }
var
  rz, rw: TFloat;
  zu, wu: Complex;
begin
  wu := Cunit(w, rw);
  if rw > MinFloat then
  begin
    zu := Cunit(z, rz);
    try
      result := Cmul(Cmul(zu, conjug(wu)),rz/rw);
    except
      on E:Exception do
        raise EComplexZerodivide.Create('Complex zero division');
    end;
  end;
end;

function Cdiv(x: TFloat; const w: Complex): Complex;
var
  rw: TFloat;
  wu: Complex;
begin
  wu := Cunit(w, rw);
  if rw > MinFloat then
  begin
    try
      result := Cmul(x/rw, conjug(wu));
    except
      on E:Exception do
        raise EComplexZerodivide.Create('Complex zero division');
    end;
  end;
end;

function Cdiv(const z: Complex;  x: TFloat): Complex;
begin
  try
    result.Re :=  z.Re/x;
    result.Im :=  z.Im/x;
  except
    on E:Exception do
        raise EComplexZerodivide.Create('Complex zero division');
  end;
end;

{ z := exp(v)  }
function Cexp(const z: Complex): Complex;
var
  r: TFloat;
begin
  r := exp1(z.re);
  result.Re := r*cos(z.im);
  result.Im := r*sin(z.im);
end;

{ z := exp(i*f)  }
function CexpIm(f: TFloat): Complex;
begin
  result.Re := cos(f);
  result.Im := sin(f);
end;

{ result = sqrt(z), result.Re > 0 }
function Csqrt(const z: Complex): Complex;
begin
  result.Re :=  sqrt(abs(0.5*(Cabs(z) + z.re)));
  result.Im :=  sqrt(abs(0.5*(Cabs(z) - z.re)));
end;

{ result = sqrt(z), result.Re > 0 }
function Csqrt(const x: TFloat): Complex;
begin
  if x < 0.0 then result := cmplx(0.0, sqrt(-x))
             else result := cmplx(sqrt(x));
end;

// z = SQR(a)
function CSqr(const z: Complex): Complex;
begin
  result.Re := SQR(z.Re) - SQR(z.Im);
  result.Im := 2.0*z.Re*z.Im
end {CSqr};

{ result = sqrt(x) }
(*
function Croot(x: TFloat): Complex;
begin
  if x < 0.0 then result := cmplx(0.0, sqrt(-x))
             else result := cmplx(sqrt(x));
end;
*)

// // -PI < theta <= PI     [From EFG Unit]
function FixAngle (const theta: TFloat): TFloat;
begin
  result := theta;

  while result > PI do
    result := result - TwoPI;

  while result <= -PI do
    result := result + TwoPI;
end {FixAngle};

// complex natural log, exponential: z = ln(a)  NOTE: principal value only
function CLn (const z: Complex): Complex;
var
  r, a: TFloat;
begin  // Abramowitz formula 4.1.2 on p. 67
  ComplexToPolar(z, r, a);
  try
    result.Re := Ln(r);
    result.Im := FixAngle(a)
  except
    On EZerodivide do   // Ln(0) raises EZerodivide
      raise EComplexLnZero.Create('CLn(0)')
  end
end {CLn};

{ s := s + z }
procedure Csum(var s: Complex; const z: Complex);
begin
  s := Cadd(s, z);
end;

{ s := s + x }
procedure Csum(var s: Complex; x: TFloat);
begin
  s := Cadd(s, x);
end;

{ r1 <= re < r2, i1 <= im < i2 }
function RndComplex(r1, r2, i1, i2: TFloat): Complex;
begin
  result.Re := r1 + (r2-r1)*random;
  result.Im := i1 + (i2-i1)*random;
end;

{ Special functions for Complex Calculations }

function cer(const z: Complex): Complex;
{ result = er(z) }

  function cerser(z: Complex): Complex;
  { series expansion of er(z) }
  const
    epsilon = 1.0E-14;
         sp = 1.7724538509055160;{ sqrt(pi) }
  var
    s, t, m, zz, zz2: Complex;
  begin
    m := cmplx1;
    t := cmplx1;
    s := cmplx0;
    zz := Cmul(z, z);
    zz2 := Cmul(zz, 2.0);
    repeat
      Csum(s, t);
      Csum(m, 2.0);
      t := Cmul(t, zz2);
      t := Cdiv(t, m);
    until (Cabs0(t) < epsilon);
    s := Cmul(s, z);
    s := Cmul(s, 2.0/sp);
    t := Cexp(zz);
    result := Csub(t, s);
  end;

  function cerabr(const z: Complex): Complex;
  { approximation of er(z): asimptotics from Abramowitz }
  const
    a1 = 0.4613135;    b1 = 0.1901635;
    a2 = 0.09999216;   b2 = 1.7844927;
    a3 = 0.002883894;  b3 = 5.5253437;
  var
    u, v, zz                : Complex;
  begin
    zz := Cmul(z, z);
    u := Cadd(zz, b1);    v := Cdiv(a1, u);
    u := Cadd(zz, b2);    u := Cdiv(a2, u);
    Csum(v, u);
    u := Cadd(zz, b3);    u := Cdiv(a3, u);
    Csum(v, u);
    result := Cmul(v, z);
  end;

begin
  if (abs(z.re) > 3.0) or (abs(z.im) > 3.9)
     then result := cerabr(z)
     else result := cerser(z);
end;{ of cer }

(*
{ TRestore1D }

constructor TRestore1D.Create(Lim: array of ILimits1D);
var
  n, i: TInt;
begin
  inherited Create;
  n := Length(Lim);
  SetLength(fLim, n);
  SetLength(fBase, n);
  for i := 0 to n-1 do
  begin
    fLim[i]:= Lim[i];
    fBase[i] := Lim[i].Base;
    Lim[i].Base := 1;
  end;
end;

destructor TRestore1D.Destroy;
var
  i: TInt;
begin
  for i := 0 to Length(fLim)-1 do
  begin
    fLim[i].Base := fBase[i];
  end;
  fLim := nil;
  fBase := nil;
  inherited;
end;

{ TRestore2D }

constructor TRestore2D.Create(Lim: array of ILimits2D);
var
  n, i: TInt;
begin
  inherited Create;
  n := Length(Lim);
  SetLength(fLim, n);
  SetLength(fBase1, n);
  SetLength(fBase2, n);
  for i := 0 to n-1 do
  begin
    fLim[i]:= Lim[i];
    fBase1[i] := Lim[i].Base1;
    fBase2[i] := Lim[i].Base2;
    Lim[i].Base1 := 1;
    Lim[i].Base2 := 1;
  end;
end;

destructor TRestore2D.Destroy;
var
  i: TInt;
begin
  for i := 0 to Length(fLim)-1 do
  begin
    fLim[i].Base1 := fBase1[i];
    fLim[i].Base2 := fBase2[i];
  end;
  fLim := nil;
  fBase1 := nil;
  fBase2 := nil;
  inherited;
end;
*)
resourcestring
  RS_RangeError = 'Demension %d is illegal';
  RS_IndexError = 'Subscript %d out of range [%d, %d]';

{unit variable setting}

initialization
  Randomize;

  MachinEps(MachEps);
  SqrtMachEps := sqrt(MachEps);
end.

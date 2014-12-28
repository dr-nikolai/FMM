{@abstract(Dynamic arrays version December, 2004)
@author(Nikolai Shokhirev <nikolai@shokhirev.com> <nikolai@u.arizona.edu>
        http://www.shokhirev.com/nikolai.html
        http://www.chem.arizona.edu/~shokhirn/nikolai.html)
@created(2003.01.01)
@lastmod(2004.12.12)
©Nikolai V. Shokhirev, 2003-2004 }
unit uDynArrays;

interface

uses
  SysUtils, uMatTypes;

const
  _= MaxInt;  // open index

type
  TFArr = array of TFloat;

  TFArrFarr = array of array of TFloat;

  TArrayType =(atZeroBased,  // [0..Dim-1] - zero-based
               atNatural,    // [1..Dim] - one-based
               atCentered,   // [-((Dim-1) div 2)..(Dim div 2)]
               atGeneric);   // [Lo..Hi] - arbitrary limits

  TSliceType = (_Col,_Row);

  EInvalidType = class(Exception);

  ELimMismatch = class(Exception);

  ENonSquareMatrix = class(Exception);

  ESingularMatrix = class(Exception);

{ Interface for 1D Limits }
  ILim1D = interface(IComment)
  ['{CBB895D7-1C56-41F7-BC7A-746B43D5771B}']
    function GetHi1: TInt;
    function GetLo1: TInt;
    function GetDim1: TInt;
    procedure SetLo1(const Value: TInt);
    property Lo1: TInt read GetLo1 write SetLo1;
    property Hi1: TInt read GetHi1;
    property Dim1: TInt read GetDim1;
  end;

{ Interface for 2D Limits }
  ILim2D = interface(ILim1D)
  ['{12E901FD-8B67-45F5-A34E-A3FC97E4DFA3}']
    function GetHi2: TInt;
    function GetLo2: TInt;
    function GetDim2: TInt;
    procedure SetLo2(const Value: TInt);
    property Lo2: TInt read GetLo2 write SetLo2;
    property Hi2: TInt read GetHi2;
    property Dim2: TInt read GetDim2;
  end;

{ Interface for 3D Limits }
  ILim3D = interface(ILim2D)
  ['{19AB1137-82DB-490F-8F99-A102BDA7EBAF}']
    function GetHi3: TInt;
    function GetLo3: TInt;
    function GetDim3: TInt;
    procedure SetLo3(const Value: TInt);
    property Lo3: TInt read GetLo3 write SetLo3;
    property Hi3: TInt read GetHi3;
    property Dim3: TInt read GetDim3;
  end;

{ Interface for 1D string Array }
  ISArr1D = interface(ILim1D)
  ['{0A65650A-AD8B-432B-B599-17A1716A037E}']
    procedure Assign(const A: ISArr1D);
    procedure Swap(i, j: TInt);
    function GetValue(i1: TInt): string;
    procedure SetValue(i1: TInt; const Value: string);
    property Value[i1: TInt]: string read GetValue write SetValue; default;
  end;

{ Interface for 1D Boolean Array }
  IBArr1D = interface(ILim1D)
  ['{ECB13ED0-8352-4431-B0AB-7C8FDE84B33E}']
    procedure Fill(C: boolean);
    procedure Assign(const A: IBArr1D);
    procedure Swap(i, j: TInt);
    function GetValue(i1: TInt): boolean;
    procedure SetValue(i1: TInt; const Value: boolean);
    property Value[i1: TInt]: boolean read GetValue write SetValue; default;
  end;

{ Interface for 2D Boolean Array }
  IBArr2D = interface(ILim2D)
  ['{FEA13DBE-68EC-4179-9D35-7EEB79A27CD6}']
    function GetSlice(i1, i2: TInt): IBArr1D;
    procedure SetSlice(i1, i2: TInt; const aValue: IBArr1D);
    procedure Fill(C: boolean);
    procedure Assign(const A: IBArr2D);
    procedure Swap(i, j: TInt; st: TSliceType);
    function GetValue(i1, i2: TInt): boolean;
    procedure SetValue(i1, i2: TInt; const Value: boolean);
    property Value[i1, i2: TInt]: boolean read GetValue write SetValue; default;
  end;

{ Interface for 1D Boolean Array }
  IBArr3D = interface(ILim3D)
  ['{C50A22C5-4CBC-415D-B4C6-732F1FCA470B}']
    procedure Fill(C: boolean);
    procedure Assign(const A: IBArr3D);
    procedure Swap(i, j: TInt; st: TSliceType);
    function GetValue(i1, i2, i3: TInt): boolean;
    procedure SetValue(i1, i2, i3: TInt; const Value: boolean);
    property Value[i1, i2, i3: TInt]: boolean read GetValue write SetValue; default;
  end;

{ Interface for 1D integer Array }
  IIArr1D = interface(ILim1D)
  ['{C245F081-177C-47B4-AF97-11C91A7A5C09}']
    procedure Fill(C: TInt);
    procedure Assign(const A: IIArr1D);
    procedure Swap(i, j: TInt);
    function GetValue(i1: TInt): TInt;
    procedure SetValue(i1: TInt; const Value: TInt);
    property Value[i1: TInt]: TInt read GetValue write SetValue; default;
  end;

{ Interface for 2D integer Array }
  IIArr2D = interface(ILim2D)
  ['{530C289C-7B12-426D-9018-58D4D1DF3FAC}']
    function GetSlice(i1, i2: TInt): IIArr1D;
    procedure SetSlice(i1, i2: TInt; const aValue: IIArr1D);
    procedure Fill(C: TInt);
    procedure Assign(const A: IIArr2D);
    procedure Swap(i, j: TInt; st: TSliceType);
    function GetValue(i1, i2: TInt): TInt;
    procedure SetValue(i1, i2: TInt; const Value: TInt);
    property Value[i1, i2: TInt]: TInt read GetValue write SetValue; default;
  end;

{ Interface for 3D integer Array }
  IIArr3D = interface(ILim3D)
  ['{83EBC67A-D1EF-4F8D-B56B-82AF1584FF3F}']
    procedure Fill(C: TInt);
    procedure Assign(const A: IIArr3D);
    function GetValue(i1, i2, i3: TInt): TInt;
    procedure SetValue(i1, i2, i3: TInt; const Value: TInt);
    property Value[i1, i2, i3: TInt]: TInt read GetValue write SetValue; default;
  end;

{ Interface for 1D float Array }
  IFArr1D = interface(ILim1D)
  ['{C245F081-177C-47B4-AF97-11C91A7A5C09}']
    procedure Fill(C: TFloat);
    procedure Times(C: TFloat);
    procedure Assign(const A: IFArr1D);
    function Norm(Normalize: boolean = false): TFloat;
    function MaxAbs: TFloat;
    function Dot(const A: IFArr1D): TFloat;
    procedure Swap(i, j: TInt);
    function GetValue(i1: TInt): TFloat;
    procedure SetValue(i1: TInt; const Value: TFloat);
    property Value[i1: TInt]: TFloat read GetValue write SetValue; default;
  end;

{ Interface for 2D float Array }
  IFArr2D = interface(ILim2D)
  ['{530C289C-7B12-426D-9018-58D4D1DF3FAC}']
    function GetSlice(i1, i2: TInt): IFArr1D;
    procedure SetSlice(i1, i2: TInt; const aValue: IFArr1D);
    procedure Fill(C: TFloat);
    procedure Times(C: TFloat);
    procedure Assign(const A: IFArr2D);
    function Norm2(i1, i2: TInt): TFloat;
    function Norm(i1, i2: TInt): TFloat;
    procedure Swap(i, j: TInt; st: TSliceType);
    function GetValue(i1, i2: TInt): TFloat;
    procedure SetValue(i1, i2: TInt; const Value: TFloat);
    property Value[i1, i2: TInt]: TFloat read GetValue write SetValue; default;
  end;

{ Interface for 3D float Array }
  IFArr3D = interface(ILim3D)
  ['{83EBC67A-D1EF-4F8D-B56B-82AF1584FF3F}']
    procedure Fill(C: TFloat);
    procedure Times(C: TFloat);
    procedure Assign(const A: IFArr3D);
    function GetValue(i1, i2, i3: TInt): TFloat;
    procedure SetValue(i1, i2, i3: TInt; const Value: TFloat);
    property Value[i1, i2, i3: TInt]: TFloat read GetValue write SetValue; default;
  end;

{<pre> Interface for Dynamic Complex 1D Array:<br>
  CVector = ReVector + i*ImVector </pre>}
  ICArr1D = interface(ILim1D)
  ['{0ADE7BFF-EFC1-4614-BE71-716509288B64}']
    procedure Conjugate;
    procedure Fill(C: Complex);
    procedure Times(C: Complex);
    procedure Assign(const A: ICArr1D);
    procedure Swap(i, j: TInt);
    function MaxAbs: TFloat;
    procedure SetRe(const Value: IFArr1D);
    procedure SetIm(const Value: IFArr1D);
    function  GetRe: IFArr1D;
    function  GetIm: IFArr1D;
    function GetValue(i: TInt): Complex;
    procedure SetValue(i: TInt; const Value: Complex);
    property Re: IFArr1D read GetRe write SetRe;
    property Im: IFArr1D read GetIm write SetIm;
    property Value[i1: TInt]: Complex read GetValue write SetValue; default;
  end;

{<pre> Interface for Dynamic Complex 2D Array:<br>
  CMatrix = ReMatrix + i*ImMatrix </pre>}
  ICArr2D = interface(ILim2D)
  ['{4A1D1265-E5F9-49CE-815F-35762A81E895}']
    procedure Conjugate;
    procedure Fill(C: Complex);
    procedure Times(C: Complex);
    procedure Assign(const A: ICArr2D);
    procedure Swap(i, j: TInt; st: TSliceType);
    function GetValue(i1, i2: TInt): Complex;
    procedure SetValue(i1, i2: TInt; const Value: Complex);
    procedure SetRe(const Value: IFArr2D);
    procedure SetIm(const Value: IFArr2D);
    function GetRe: IFArr2D;
    function GetIm: IFArr2D;
    property Re: IFArr2D read GetRe write SetRe;
    property Im: IFArr2D read GetIm write SetIm;
    property Value[i1, i2: TInt]: Complex read GetValue write SetValue; default;
  end;

  IEigenSys = interface(IComment)
  ['{99B2FEC7-F026-46C2-A6D1-D6D32EE9086D}']
    procedure SetNames(const Value: ISArr1D);
    procedure SetValues(const Value: IFArr1D);
    procedure SetVectors(const Value: IFArr2D);
    function  GetNames: ISArr1D;
    function  GetValues: IFArr1D;
    function  GetVectors: IFArr2D;
    property Names: ISArr1D read GetNames write SetNames;
    property Values: IFArr1D read GetValues write SetValues;
    property Vectors: IFArr2D read GetVectors write SetVectors;
  end;
  
  IHEigenSys = interface(IComment)
  ['{A5276476-AFD8-42AB-B428-E65A821F0234}']
    procedure SetNames(const Value: ISArr1D);
    procedure SetValues(const Value: IFArr1D);
    procedure SetVectors(const Value: ICArr2D);
    function  GetNames: ISArr1D;
    function  GetValues: IFArr1D;
    function  GetVectors: ICArr2D;
    property Names: ISArr1D read GetNames write SetNames;
    property Values: IFArr1D read GetValues write SetValues;
    property Vectors: ICArr2D read GetVectors write SetVectors;
  end;

  ISVDSys = interface(IComment)
  ['{891A2307-1A3D-4D96-9027-84ABA06FEA83}']
    procedure SetValues(const Value: IFArr1D);
    procedure SetVvectors(const Value: IFArr2D);
    procedure SetUvectors(const Value: IFArr2D);
    function  GetValues: IFArr1D;
    function  GetVvectors: IFArr2D;
    function  GetUvectors: IFArr2D;
    property Values: IFArr1D read GetValues write SetValues;
    property Uvectors: IFArr2D read GetUvectors write SetUvectors;
    property Vvectors: IFArr2D read GetVvectors write SetVvectors;
  end;

  { Object for 1D Limits }
  TLim1D = class(TComment, ILim1D, IRestore)
  protected
    fHi1: TInt;
    fLo1: TInt;
    fDim1: TInt;
    function GetHi1: TInt;
    function GetLo1: TInt;
    function GetDim1: TInt;
    procedure SetLo1(const Value: TInt);
  public
    constructor Create(aHi1: TInt); overload;
    constructor Create(aLo1, aHi1: TInt); overload;
    property Hi1: TInt read GetHi1;
    property Lo1: TInt read GetLo1 write SetLo1;
    property Dim1: TInt read GetDim1;
  end;

{ Object for 2D Limits }
  TLim2D = class(TLim1D, ILim2D, IRestore)
  protected
    fHi2: TInt;
    fLo2: TInt;
    fDim2: TInt;
    function GetHi2: TInt;
    function GetLo2: TInt;
    function GetDim2: TInt;
    procedure SetLo2(const Value: TInt);
  public
    constructor Create(aHi1, aHi2: TInt); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2: TInt); overload;
    property Hi2: TInt read GetHi2;
    property Lo2: TInt read GetLo2 write SetLo2;
    property Dim2: TInt read GetDim2;
  end;

{ Object for 2D Limits }
  TLim3D = class(TLim2D, ILim3D)
  protected
    fHi3: TInt;
    fLo3: TInt;
    fDim3: TInt;
    function GetHi3: TInt;
    function GetLo3: TInt;
    function GetDim3: TInt;
    procedure SetLo3(const Value: TInt);
  public
    constructor Create(aHi1, aHi2, aHi3: TInt); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2, aLo3, aHi3: TInt); overload;
    property Hi3: TInt read GetHi3;
    property Lo3: TInt read GetLo3 write SetLo3;
    property Dim3: TInt read GetDim3;
  end;

  { Limits1D Restore object }
  TRestore1D = class(TinterfacedObject, IRestore)
  private
    fLim: array of ILim1D;
    fBase1: array of integer;
  public
    constructor Create(Lim: array of ILim1D);
    destructor Destroy; override;
  end;

  { Limits2D Restore object }
  TRestore2D = class(TinterfacedObject, Irestore)
  private
    fLim: array of ILim2D;
    fBase1, fBase2: array of integer;
  public
    constructor Create(Lim: array of ILim2D);
    destructor Destroy; override;
  end;

{ Object for 1D string Array }
  TSArr1D = class(TLim1D, ISArr1D)
  private
    fValue: array of string;
  protected
    function GetValue(i1: TInt): string;
    procedure SetValue(i1: TInt; const Value: string);
  public
    constructor Create(A: ISArr1D; CopyData: boolean = false); overload;
    constructor Create(aHi1: TInt); overload;
    constructor Create(aLo1, aHi1: TInt); overload;
    destructor Destroy; override;
    procedure Assign(const A: ISArr1D);
    procedure Swap(i, j: TInt);
    property Value[i1: TInt]: string read GetValue write SetValue; default;
  end;

{ Object for 1D boolean Array }
  TBArr1D = class(TLim1D, IBArr1D)
  private
    fValue: array of boolean;
  protected
    function GetValue(i1: TInt): boolean;
    procedure SetValue(i1: TInt; const Value: boolean);
  public
    constructor Create(A: IBArr1D; CopyData: boolean = false); overload;
    constructor Create(aHi1: TInt); overload;
    constructor Create(aLo1, aHi1: TInt); overload;
    destructor Destroy; override;
    procedure Fill(C: boolean);
    procedure Assign(const A: IBArr1D);
    procedure Swap(i, j: TInt);
    property Value[i1: TInt]: boolean read GetValue write SetValue; default;
  end;

{ Object for 2D boolean Array }
  TBArr2D = class(TLim2D, IBArr2D)
  private
    fValue: array of array of boolean;
  protected
    function GetSlice(i1, i2: TInt): IBArr1D;
    procedure SetSlice(i1, i2: TInt; const aValue: IBArr1D);
    function GetValue(i1, i2: TInt): boolean;
    procedure SetValue(i1, i2: TInt; const Value: boolean);
  public
    constructor Create(A: IBArr2D; CopyData: boolean = false); overload;
    constructor Create(aHi1, aHi2: TInt); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2: TInt); overload;
    destructor Destroy; override;
    procedure Fill(C: boolean);
    procedure Assign(const A: IBArr2D);
    procedure Swap(i, j: TInt; st: TSliceType);
    property Value[i1, i2: TInt]: boolean read GetValue write SetValue; default;
  end;

{ Object for 1D integer Array }
  TIArr1D = class(TLim1D, IIArr1D)
  private
    fValue: array of TInt;
  protected
    function GetValue(i1: TInt): TInt;
    procedure SetValue(i1: TInt; const Value: TInt);
  public
    constructor Create(A: IIArr1D; CopyData: boolean = false); overload;
    constructor Create(aHi1: TInt); overload;
    constructor Create(aLo1, aHi1: TInt); overload;
    destructor Destroy; override;
    procedure Fill(C: TInt);
    procedure Assign(const A: IIArr1D);
    procedure Swap(i, j: TInt);
    property Value[i1: TInt]: TInt read GetValue write SetValue; default;
  end;

{ Object for 2D integer Array }
  TIArr2D = class(TLim2D, IIArr2D)
  private
    fValue: array of array of TInt;
  protected
    function GetSlice(i1, i2: TInt): IIArr1D;
    procedure SetSlice(i1, i2: TInt; const aValue: IIArr1D);
    function GetValue(i1, i2: TInt): TInt;
    procedure SetValue(i1, i2: TInt; const Value: TInt);
  public
    constructor Create(A: IIArr2D; CopyData: boolean = false); overload;
    constructor Create(aHi1, aHi2: TInt); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2: TInt); overload;
    destructor Destroy; override;
    procedure Fill(C: TInt);
    procedure Assign(const A: IIArr2D);
    procedure Swap(i, j: TInt; st: TSliceType);
    property Value[i1, i2: TInt]: TInt read GetValue write SetValue; default;
  end;

{ Object for 3D integer Array }
  TIArr3D = class(TLim3D, IIArr3D)
  private
    fValue: array of array of array of TInt;
  protected
    function GetValue(i1, i2, i3: TInt): TInt;
    procedure SetValue(i1, i2, i3: TInt; const Value: TInt);
  public
    constructor Create(A: IIArr3D; CopyData: boolean = false); overload;
    constructor Create(aHi1, aHi2, aHi3: TInt); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2, aLo3, aHi3: TInt); overload;
    destructor Destroy; override;
    procedure Fill(C: TInt);
    procedure Assign(const A: IIArr3D);
    property Value[i1, i2, i3: TInt]: TInt read GetValue write SetValue; default;
  end;

{ Object for 1D float Array }
  TFArr1D = class(TLim1D, IFArr1D)
  protected
    fInternalValue: boolean;
    fValue: TFArr;
    function GetValue(i1: TInt): TFloat; virtual;
    procedure SetValue(i1: TInt; const Value: TFloat); virtual;
  public
    constructor Create(A: IFArr1D; CopyData: boolean = false); overload;
    constructor Create(aHi1: TInt); overload;
    constructor Create(aLo1, aHi1: TInt); overload;
    constructor Create(aLo1, aHi1: TInt; var A: TFArr); overload;
    destructor Destroy; override;
    procedure Fill(C: TFloat);
    procedure Times(C: TFloat);
    procedure Assign(const A: IFArr1D);
    function Norm(Normalize: boolean = false): TFloat;
    function MaxAbs: TFloat;
    function Dot(const A: IFArr1D): TFloat;
    procedure Swap(i, j: TInt);
    property Value[i1: TInt]: TFloat read GetValue write SetValue; default;
  end;

{ Object for 2D float Array }
  TFArr2D = class(TLim2D, IFArr2D)
  private
    fValue: array of array of TFloat;
  protected
    function GetSlice(i1, i2: TInt): IFArr1D;
    procedure SetSlice(i1, i2: TInt; const aValue: IFArr1D);
    function GetValue(i1, i2: TInt): TFloat;
    procedure SetValue(i1, i2: TInt; const Value: TFloat);
  public
    constructor Create(A: IFArr2D; CopyData: boolean = false); overload;
    constructor Create(aHi1, aHi2: TInt); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2: TInt); overload;
    destructor Destroy; override;
    procedure Fill(C: TFloat);
    procedure Times(C: TFloat);
    procedure Assign(const A: IFArr2D);
    function Norm2(i1, i2: TInt): TFloat;
    function Norm(i1, i2: TInt): TFloat;
    procedure Swap(i, j: TInt; st: TSliceType);
    property Value[i1, i2: TInt]: TFloat read GetValue write SetValue; default;
  end;

{ Object for 2D float Array }
  TFArr3D = class(TLim3D, IFArr3D)
  private
    fValue: array of array of array of TFloat;
  protected
    function GetValue(i1, i2, i3: TInt): TFloat;
    procedure SetValue(i1, i2, i3: TInt; const Value: TFloat);
  public
    constructor Create(A: IFArr3D; CopyData: boolean = false); overload;
    constructor Create(aHi1, aHi2, aHi3: TInt); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2, aLo3, aHi3: TInt); overload;
    destructor Destroy; override;
    procedure Fill(C: TFloat);
    procedure Times(C: TFloat);
    procedure Assign(const A: IFArr3D);
    property Value[i1, i2, i3: TInt]: TFloat read GetValue write SetValue; default;
  end;

    { Object for 1D complex Array }
  TCArr1D = class(TComment, ICArr1D)
  private
    fRe: IFArr1D;
    fIm: IFArr1D;
  protected
    function GetHi1: TInt;
    function GetLo1: TInt;
    function GetDim1: TInt;
    procedure SetLo1(const Value: TInt);
    procedure SetRe(const Value: IFArr1D);
    procedure SetIm(const Value: IFArr1D);
    function GetRe: IFArr1D;
    function GetIm: IFArr1D;
    function GetValue(i1: TInt): Complex;
    procedure SetValue(i1: TInt; const Value: Complex);
  public
    constructor Create(A: ICArr1D; CopyData: boolean = false); overload;
    constructor Create(aHi1: TInt); overload;
    constructor Create(aLo1, aHi1: TInt); overload;
    destructor Destroy; override;
    procedure Conjugate;
    procedure Fill(C: Complex);
    procedure Times(C: Complex);
    procedure Assign(const A: ICArr1D);
    procedure Swap(i, j: TInt);
    function MaxAbs: TFloat;
    property Re: IFArr1D read GetRe write SetRe;
    property Im: IFArr1D read GetIm write SetIm;
    property Hi1: TInt read GetHi1;
    property Lo1: TInt read GetLo1 write SetLo1;
    property Dim1: TInt read GetDim1;
    property Value[i1: TInt]: Complex read GetValue write SetValue; default;
  end;

{ Object for 2D complex Array }
  TCArr2D = class(TComment, ICArr2D)
  private
    fRe: IFArr2D;
    fIm: IFArr2D;
  protected
    function GetHi1: TInt;
    function GetLo1: TInt;
    function GetDim1: TInt;
    procedure SetLo1(const Value: TInt);
    function GetHi2: TInt;
    function GetLo2: TInt;
    function GetDim2: TInt;
    procedure SetLo2(const Value: TInt);
    procedure SetRe(const Value: IFArr2D);
    procedure SetIm(const Value: IFArr2D);
    function  GetRe: IFArr2D;
    function  GetIm: IFArr2D;
    function GetValue(i1, i2: TInt): Complex;
    procedure SetValue(i1, i2: TInt; const Value: Complex);
  public
    constructor Create(A: ICArr2D; CopyData: boolean = false); overload;
    constructor Create(aHi1, aHi2: TInt); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2: TInt); overload;
    destructor Destroy; override;
    procedure Conjugate;
    procedure Fill(C: Complex);
    procedure Times(C: Complex);
    procedure Assign(const A: ICArr2D);
    procedure Swap(i, j: TInt; st: TSliceType);
    property Re: IFArr2D read GetRe write SetRe;
    property Im: IFArr2D read GetIm write SetIm;
    property Hi1: TInt read GetHi1;
    property Lo1: TInt read GetLo1 write SetLo1;
    property Dim1: TInt read GetDim1;
    property Hi2: TInt read GetHi2;
    property Lo2: TInt read GetLo2 write SetLo2;
    property Dim2: TInt read GetDim2;
    property Value[i1, i2: TInt]: Complex read GetValue write SetValue; default;
  end;

{ Object for Real eigensystem }
  TEigenSys = class(TComment, IEigenSys)
  private
    fValues: IFArr1D;
    fVectors: IFArr2D;
    fNames: ISArr1D;
  protected
    procedure SetNames(const Value: ISArr1D);
    procedure SetValues(const Value: IFArr1D);
    procedure SetVectors(const Value: IFArr2D);
    function  GetNames: ISArr1D;
    function  GetValues: IFArr1D;
    function  GetVectors: IFArr2D;
  public
    constructor Create(const a: ILim2D); overload;
    constructor Create(aHi1: TInt); overload;
    constructor Create(aLo1, aHi1: TInt); overload;
    destructor Destroy; override;
    property Names: ISArr1D read GetNames write SetNames;
    property Values: IFArr1D read GetValues write SetValues;
    property Vectors: IFArr2D read GetVectors write SetVectors;
  end;

{ Object for Hermitian eigensystem }
  THEigenSys = class(TComment, IHEigenSys)
  private
    fValues: IFArr1D;
    fVectors: ICArr2D;
    fNames: ISArr1D;
  protected
    procedure SetNames(const Value: ISArr1D);
    procedure SetValues(const Value: IFArr1D);
    procedure SetVectors(const Value: ICArr2D);
    function  GetNames: ISArr1D;
    function  GetValues: IFArr1D;
    function  GetVectors: ICArr2D;
  public
    constructor Create(const a: ILim2D); overload;
    constructor Create(aHi1: TInt); overload;
    constructor Create(aLo1, aHi1: TInt); overload;
    destructor Destroy; override;
    property Names: ISArr1D read GetNames write SetNames;
    property Values: IFArr1D read GetValues write SetValues;
    property Vectors: ICArr2D read GetVectors write SetVectors;
  end;

{ Object for SVD system }
  TSVDSys = class(TComment, ISVDSys)
  private
    fValues: IFArr1D;
    fVvectors: IFArr2D;
    fUvectors: IFArr2D;
  protected
    procedure SetValues(const Value: IFArr1D);
    procedure SetVvectors(const Value: IFArr2D);
    procedure SetUvectors(const Value: IFArr2D);
    function  GetValues: IFArr1D;
    function  GetVvectors: IFArr2D;
    function  GetUvectors: IFArr2D;
  public
    constructor Create(const A: IFArr2D; fullmat: boolean = false); overload;
    constructor Create(aHi1, aHi2: TInt; fullmat: boolean = false); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2: TInt; fullmat: boolean = false); overload;
    property Values: IFArr1D read GetValues write SetValues;
    property Uvectors: IFArr2D read GetUvectors write SetUvectors;
    property Vvectors: IFArr2D read GetVvectors write SetVvectors;
  end;

(*
{ Object for 2D float Array }
  TSparse2D = class(TLim2D, IFArr2D)
  private
    fValue: array of array of TFloat;
  protected
    function GetSlice(i1, i2: TInt): IFArr1D;
    procedure SetSlice(i1, i2: TInt; const aValue: IFArr1D);
    function GetValue(i1, i2: TInt): TFloat;
    procedure SetValue(i1, i2: TInt; const Value: TFloat);
  public
    constructor Create(A: IFArr2D; CopyData: boolean = false); overload;
    constructor Create(aHi1, aHi2: TInt); overload;
    constructor Create(aLo1, aHi1, aLo2, aHi2: TInt); overload;
    destructor Destroy; override;
    procedure Fill(C: TFloat);
    procedure Times(C: TFloat);
    procedure Assign(const A: IFArr2D; MatchLim: boolean = true);
    function Norm2(i1, i2: TInt): TFloat;
    function Norm(i1, i2: TInt): TFloat;
    procedure Swap(i, j: TInt; st: TSliceType);
    property Value[i1, i2: TInt]: TFloat read GetValue write SetValue; default;
//    property Vector[i1, i2: TInt]: IFArr1D read GetSlice write SetSlice;
  end;
*)

{ A.Lo1 = A.Lo2 and A.Hi1 = A.Hi2 }
function IsSquare(const A: ILim2D): boolean;

{ A1.Lo1 = A2.Lo1 and A1.Hi1 = A2.Hi1 }
function SameLim(const A1, A2: ILim1D): boolean; overload;
{ A1.Lo1 = A2.Lo1 and A1.Hi1 = A2.Hi1 and A1.Lo2 = A2.Lo2 and A1.Hi2 = A2.Hi2 }
function SameLim(const A1, A2: ILim2D): boolean; overload;
{ A1.Lo1 = A2.Lo1 and A1.Hi1 = A2.Hi1 and A1.Lo2 = A2.Lo2 and A1.Hi2 = A2.Hi2
  and A1.Lo3 = A2.Lo3 and A1.Hi3 = A2.Hi3  }
function SameLim(const A1, A2: ILim3D): boolean; overload;
{ Same limits with transposed:    <br>
  A1.Lo1 = A2.Lo2 and A1.Hi1 = A2.Hi2 and A1.Lo2 = A2.Lo1 and A1.Hi2 = A2.Hi1 }
function SameLimT(const A1, A2: ILim2D): boolean;

{<pre>
  TArrayType =(atZeroBased,  // [0..Dim-1] - zero-based         <br>
               atNatural,    // [1..Dim] - one-based            <br>
               atCentered,   // [-((Dim-1) div 2)..(Dim div 2)] <br>
               atGeneric);   // [Lo..Hi] - arbitrary limits     </pre>}
function ArrType(const A: ILim1D): TArrayType;

ResourceString
  RS_LimMismatch = 'Lim Mismatch';

implementation

uses
  Math;

ResourceString
  RS_IndexError = '%s index = %d : out of range [%d, %d]';
  RS_IndexOutOfRange = 'Index out of range';
//  RS_LimMismatch = 'Lim Mismatch';
  RS_ZeroVector = 'Cannot Normalize: Zero Vector';
  RS_TwoBlnkIndices = 'Two Blnk Indices: [_,_] are Illegal';
  RS_NoBlnkIndices = 'One index must be blank';

{ A.Lo1 = A.Lo2 and A.Hi1 = A.Hi2 }
function IsSquare(const A: ILim2D): boolean;
begin
  result := (A.Lo1 = A.Lo2) and (A.Hi1 = A.Hi2);
end;
{ A1.Lo1 = A2.Lo1 and A1.Hi1 = A2.Hi1 }
function SameLim(const A1, A2: ILim1D): boolean;
begin
  result := (A1.Lo1 = A2.Lo1) and (A1.Hi1 = A2.Hi1);
end;

{ A1.Lo1 = A2.Lo1 and A1.Hi1 = A2.Hi1 and A1.Lo2 = A2.Lo2 and A1.Hi2 = A2.Hi2 }
function SameLim(const A1, A2: ILim2D): boolean;
begin
  result := (A1.Lo1 = A2.Lo1) and (A1.Hi1 = A2.Hi1) and
            (A1.Lo2 = A2.Lo2) and (A1.Hi2 = A2.Hi2);
end;

{ A1.Lo1 = A2.Lo1 and A1.Hi1 = A2.Hi1 and A1.Lo2 = A2.Lo2 and A1.Hi2 = A2.Hi2
  and A1.Lo3 = A2.Lo3 and A1.Hi3 = A2.Hi3  }
function SameLim(const A1, A2: ILim3D): boolean;
begin
  result := (A1.Lo1 = A2.Lo1) and (A1.Hi1 = A2.Hi1) and
            (A1.Lo2 = A2.Lo2) and (A1.Hi2 = A2.Hi2) and
            (A1.Lo3 = A2.Lo3) and (A1.Hi3 = A2.Hi3);
end;

{ A1.Lo1 = A2.Lo2 and A1.Hi1 = A2.Hi2 and A1.Lo2 = A2.Lo1 and A1.Hi2 = A2.Hi1 }
function SameLimT(const A1, A2: ILim2D): boolean;
begin
  result := (A1.Lo1 = A2.Lo2) and (A1.Hi1 = A2.Hi2) and
            (A1.Lo2 = A2.Lo1) and (A1.Hi2 = A2.Hi1);
end;

{<pre>
  TArrayType =(atZeroBased,  // [0..Dim-1] - zero-based         <br>
               atNatural,    // [1..Dim] - one-based            <br>
               atCentered,   // [-((Dim-1) div 2)..(Dim div 2)] <br>
               atGeneric);   // [Lo..Hi] - arbitrary limits     </pre>}
function ArrType(const A: ILim1D): TArrayType;
begin
  case A.Lo1 of
  0: result := atZeroBased;
  1: result := atNatural
  else
    if (A.Lo1 = -((A.Dim1-1) div 2) ) and (A.Hi1=(A.Dim1 div 2)) then
      result := atCentered
    else
      result := atGeneric;
  end;
end;

{ TLim1D }

constructor TLim1D.Create(aLo1, aHi1: TInt);
begin
  inherited Create;
  fHi1 := aHi1;
  fLo1 := aLo1;
  fDim1 := fHi1 - fLo1 + 1;
end;

constructor TLim1D.Create(aHi1: TInt);
begin
  Create(1, aHi1);
end;

function TLim1D.GetDim1: TInt;
begin
  result := fHi1 - fLo1 + 1;
end;

function TLim1D.GetHi1: TInt;
begin
  result := fHi1;
end;

function TLim1D.GetLo1: TInt;
begin
  result := fLo1;
end;

procedure TLim1D.SetLo1(const Value: TInt);
begin
  fLo1 := Value;
  fHi1 := fLo1 + fDim1 - 1;
end;

{ TLim2D }

constructor TLim2D.Create(aLo1, aHi1, aLo2, aHi2: TInt);
begin
  inherited Create(aLo1, aHi1);
  fHi2 := aHi2;
  fLo2 := aLo2;
  fDim2 := fHi2 - fLo2 + 1;
end;

constructor TLim2D.Create(aHi1, aHi2: TInt);
begin
  Create(1, aHi1,1, aHi2);
end;

function TLim2D.GetDim2: TInt;
begin
  result := fHi2 - fLo2 + 1;
end;

function TLim2D.GetHi2: TInt;
begin
  result := fHi2;
end;

function TLim2D.GetLo2: TInt;
begin
  result := fLo2;
end;

procedure TLim2D.SetLo2(const Value: TInt);
begin
  fLo2 := Value;
  fHi2 := fLo2 + fDim2 - 1;
end;

{ TLim3D }

constructor TLim3D.Create(aLo1, aHi1, aLo2, aHi2, aLo3, aHi3: TInt);
begin
  inherited Create(aLo1, aHi1, aLo2, aHi2);
  fHi3 := aHi3;
  fLo3 := aLo3;
  fDim3 := fHi3 - fLo3 + 1;
end;

constructor TLim3D.Create(aHi1, aHi2, aHi3: TInt);
begin
  Create(1, aHi1, 1, aHi2, 1, aHi3);
end;

function TLim3D.GetDim3: TInt;
begin
  result := fHi3 - fLo3 + 1;
end;

function TLim3D.GetHi3: TInt;
begin
  result := fHi3;
end;

function TLim3D.GetLo3: TInt;
begin
  result := fLo3;
end;

procedure TLim3D.SetLo3(const Value: TInt);
begin
  fLo3 := Value;
  fHi3 := fLo3 + fDim3 - 1;
end;

{ TRestore1D }

constructor TRestore1D.Create(Lim: array of ILim1D);
var
  n, i: TInt;
begin
  inherited Create;
  n := Length(Lim);
  SetLength(fLim, n);
  SetLength(fBase1, n);
  for i := 0 to n-1 do
  begin
    fLim[i]:= Lim[i];
    fBase1[i] := Lim[i].Lo1;
    Lim[i].Lo1 := 1;
  end;
end;

destructor TRestore1D.Destroy;
var
  i: TInt;
begin
  for i := 0 to Length(fLim)-1 do
  begin
    fLim[i].Lo1 := fBase1[i];
  end;
  fLim := nil;
  fBase1 := nil;
  inherited;
end;

{ TRestore2D }

constructor TRestore2D.Create(Lim: array of ILim2D);
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
    fBase1[i] := Lim[i].Lo1;
    fBase2[i] := Lim[i].Lo2;
    Lim[i].Lo1 := 1;
    Lim[i].Lo2 := 1;
  end;
end;

destructor TRestore2D.Destroy;
var
  i: TInt;
begin
  for i := 0 to Length(fLim)-1 do
  begin
    fLim[i].Lo1 := fBase1[i];
    fLim[i].Lo2 := fBase2[i];
  end;
  fLim := nil;
  fBase1 := nil;
  fBase2 := nil;
  inherited;
end;

{ TSArr1D }

constructor TSArr1D.Create(aHi1: TInt);
begin
  Create(1, aHi1);
end;

procedure TSArr1D.Assign(const A: ISArr1D);
var
  i1: TInt;
begin
  Assert( (Hi1=A.Hi1) and (Lo1=A.Lo1), RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
      Value[i1] := A[i1];
end;

constructor TSArr1D.Create(aLo1, aHi1: TInt);
begin
  inherited;
  SetLength(fValue, fHi1-fLo1+1);
end;

constructor TSArr1D.Create(A: ISArr1D; CopyData: boolean = false);
begin
  Create(A.Lo1, A.Hi1);
  if CopyData then Assign(A);
end;

destructor TSArr1D.Destroy;
begin
  SetLength(fValue, 0);
  fValue := nil;
  inherited;
end;

function TSArr1D.GetValue(i1: TInt): string;
begin
  Assert( (i1 >= fLo1) and (i1 <= fHi1), RS_IndexOutOfRange);
  result := fValue[i1-fLo1];
end;

procedure TSArr1D.SetValue(i1: TInt; const Value: string);
begin
  Assert( (i1 >= fLo1) and (i1 <= fHi1), RS_IndexOutOfRange);
  fValue[i1-fLo1] := Value;
end;

procedure TSArr1D.Swap(i, j: TInt);
var
  v: string;
begin
  v := Value[i];
  Value[i] := Value[j];
  Value[j] := v
end;

{ TBArr1D }

constructor TBArr1D.Create(aHi1: TInt);
begin
  Create(1, aHi1);
end;

procedure TBArr1D.Assign(const A: IBArr1D);
var
  i1: TInt;
begin
  Assert( (Hi1=A.Hi1) and (Lo1=A.Lo1), RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
      Value[i1] := A[i1];
end;

constructor TBArr1D.Create(aLo1, aHi1: TInt);
begin
  inherited;
  SetLength(fValue, fHi1-fLo1+1);
end;

constructor TBArr1D.Create(A: IBArr1D; CopyData: boolean = false);
begin
  Create(A.Lo1, A.Hi1);
  if CopyData then Assign(A);
end;

destructor TBArr1D.Destroy;
begin
  SetLength(fValue, 0);
  fValue := nil;
  inherited;
end;

function TBArr1D.GetValue(i1: TInt): boolean;
begin
  Assert( (i1 >= fLo1) and (i1 <= fHi1), RS_IndexOutOfRange);
  result := fValue[i1-fLo1];
end;

procedure TBArr1D.SetValue(i1: TInt; const Value: boolean);
begin
  Assert( (i1 >= fLo1) and (i1 <= fHi1), RS_IndexOutOfRange);
  fValue[i1-fLo1] := Value;
end;

procedure TBArr1D.Swap(i, j: TInt);
var
  v: boolean;
begin
  v := Value[i];
  Value[i] := Value[j];
  Value[j] := v
end;

procedure TBArr1D.Fill(C: boolean);
var
  i1: TInt;
begin
  for i1 := Lo1 to Hi1 do
    Value[i1] := C;
end;

{ TBArr2D }

constructor TBArr2D.Create(aHi1, aHi2: TInt);
begin
  Create(1,aHi1,1,aHi2);
end;

procedure TBArr2D.Assign(const A: IBArr2D);
var
  i1, i2: TInt;
begin
  Assert((Hi1=A.Hi1) and (Lo1=A.Lo1) and (Hi2=A.Hi2) and (Lo2=A.Lo2),RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
    for i2 := A.Lo2 to A.Hi2 do
      Value[i1,i2] := A[i1,i2];
end;

constructor TBArr2D.Create(aLo1, aHi1, aLo2, aHi2: TInt);
begin
  inherited;
  SetLength(fValue, fHi2-fLo2+1, fHi1-fLo1+1);
end;

constructor TBArr2D.Create(A: IBArr2D; CopyData: boolean = false);
begin
  Create(A.Lo1, A.Hi1, A.Lo2, A.Hi2);
  if CopyData then Assign(A);
end;

destructor TBArr2D.Destroy;
begin
  SetLength(fValue, 0, 0);
  fValue := nil;
  inherited;
end;

function TBArr2D.GetValue(i1, i2: TInt): boolean;
begin
  Assert((i1>=fLo1) and (i1<=fHi1) and (i2>=fLo2) and (i2 <= fHi2),RS_IndexOutOfRange);
  result := fValue[i2-fLo2,i1-fLo1];
end;

procedure TBArr2D.SetValue(i1, i2: TInt; const Value: boolean);
begin
  Assert((i1>=fLo1)and(i1<=fHi1)and(i2>=fLo2)and(i2 <= fHi2),RS_IndexOutOfRange);
  fValue[i2-fLo2,i1-fLo1] := Value;
end;

procedure TBArr2D.Swap(i, j: TInt; st: TSliceType);
var
  k: TInt;
  v: boolean;
begin
  case st of
  _Row:for k := Lo2 to Hi2 do
        begin
          v := Value[i,k];
          Value[i,k] := Value[j,k];
          Value[j,k] := v;
        end;
  _Col:for k := Lo1 to Hi1 do
        begin
          v := Value[k,i];
          Value[k,i] := Value[k,j];
          Value[k,j] := v;
        end;
  else
//    Raise ERangeError.CreateFmt(RS_DimError,[idx, 2]);
  end;
end;

procedure TBArr2D.Fill(C: boolean);
var
  i1, i2: TInt;
begin
  for i1 := Lo1 to Hi1 do
    for i2 := Lo2 to Hi2 do
        Value[i1,i2] := C;
end;

function TBArr2D.GetSlice(i1, i2: TInt): IBArr1D;
var
  t: IBArr1D;
  k: TInt;
begin
  if i1 = _ then
  begin
    if i2 = _ then
      Raise ERangeError.Create(RS_TwoBlnkIndices);
    t := TBArr1D.Create(Lo1,Hi1);
    for k := t.Lo1 to t.Hi1 do t[k] := Value[k,i2];
  end else
  if i2 = _ then
  begin
    t := TBArr1D.Create(Lo2,Hi2);
    for k := t.Lo1 to t.Hi1 do t[k] := Value[i1,k];
  end else
    Raise ERangeError.Create(RS_NoBlnkIndices);
  result := t;
end;

procedure TBArr2D.SetSlice(i1, i2: TInt; const aValue: IBArr1D);
var
  k: TInt;
begin
  if i1 = _ then
  begin
    if i2 = _ then
      Raise ERangeError.Create(RS_TwoBlnkIndices);
    if (aValue.Lo1<>Lo1) or (aValue.Hi1 <> Hi1) then
      Raise ERangeError.Create(RS_LimMismatch);
    for k := aValue.Lo1 to aValue.Hi1 do Value[k,i2] := aValue[k];
  end else
  if i2 = _ then
  begin
    if (aValue.Lo1<>Lo2) or (aValue.Hi1 <> Hi2) then
      Raise ERangeError.Create(RS_LimMismatch);
    for k := aValue.Lo1 to aValue.Hi1 do Value[i1,k] := aValue[k];
  end else
      Raise ERangeError.Create(RS_NoBlnkIndices);
end;

{ TIArr1D }

procedure TIArr1D.Assign(const A: IIArr1D);
var
  i1: TInt;
begin
  Assert( (Hi1=A.Hi1) and (Lo1=A.Lo1), RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
      Value[i1] := A[i1];
end;

constructor TIArr1D.Create(aLo1, aHi1: TInt);
begin
  inherited;
  SetLength(fValue, fHi1-fLo1+1);
end;

constructor TIArr1D.Create(aHi1: TInt);
begin
  Create(1, aHi1);
end;

constructor TIArr1D.Create(A: IIArr1D; CopyData: boolean = false);
begin
  Create(A.Lo1, A.Hi1);
  if CopyData then Assign(A);
end;

destructor TIArr1D.Destroy;
begin
  SetLength(fValue, 0);
  fValue := nil;
  inherited;
end;

procedure TIArr1D.Fill(C: TInt);
var
  i1: TInt;
begin
  for i1 := Lo1 to Hi1 do
    Value[i1] := C;
end;

function TIArr1D.GetValue(i1: TInt): TInt;
begin
  Assert((i1 >= fLo1) and (i1 <= fHi1), RS_IndexOutOfRange);
  result := fValue[i1-fLo1];
end;

procedure TIArr1D.SetValue(i1: TInt; const Value: TInt);
begin
  Assert((i1 >= fLo1) and (i1 <= fHi1), RS_IndexOutOfRange);
  fValue[i1-fLo1] := Value;
end;

procedure TIArr1D.Swap(i, j: TInt);
var
  v: TInt;
begin
  v := Value[i];
  Value[i] := Value[j];
  Value[j] := v
end;

{ TIArr2D }

procedure TIArr2D.Assign(const A: IIArr2D);
var
  i1, i2: TInt;
begin
  Assert((Hi1=A.Hi1) and (Lo1=A.Lo1) and (Hi2=A.Hi2) and (Lo2=A.Lo2),RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
    for i2 := A.Lo2 to A.Hi2 do
      Value[i1,i2] := A[i1,i2];
end;

constructor TIArr2D.Create(aLo1, aHi1, aLo2, aHi2: TInt);
begin
  inherited;
  SetLength(fValue, fHi2-fLo2+1, fHi1-fLo1+1);
end;

constructor TIArr2D.Create(aHi1, aHi2: TInt);
begin
  Create(1, aHi1, 1, aHi2);
end;

constructor TIArr2D.Create(A: IIArr2D; CopyData: boolean = false);
begin
  Create(A.Lo1, A.Hi1, A.Lo2, A.Hi2);
  if CopyData then Assign(A);
end;

destructor TIArr2D.Destroy;
begin
  SetLength(fValue, 0, 0);
  fValue := nil;
  inherited;
end;

procedure TIArr2D.Fill(C: TInt);
var
  i1, i2: TInt;
begin
  for i1 := Lo1 to Hi1 do
    for i2 := Lo2 to Hi2 do
        Value[i1,i2] := C;
end;

function TIArr2D.GetValue(i1, i2: TInt): TInt;
begin
  Assert((i1>=fLo1) and (i1<=fHi1) and (i2>=fLo2) and (i2 <= fHi2),RS_IndexOutOfRange);
  result := fValue[i2-fLo2,i1-fLo1];
end;

procedure TIArr2D.SetValue(i1, i2: TInt; const Value: TInt);
begin
  Assert((i1>=fLo1) and (i1<=fHi1) and (i2>=fLo2) and (i2 <= fHi2),RS_IndexOutOfRange);
  fValue[i2-fLo2,i1-fLo1] := Value;
end;

function TIArr2D.GetSlice(i1, i2: TInt): IIArr1D;
var
  t: IIArr1D;
  k: TInt;
begin
  if i1 = _ then
  begin
    if i2 = _ then
      Raise ERangeError.Create(RS_TwoBlnkIndices);
    t := TIArr1D.Create(Lo1,Hi1);
    for k := t.Lo1 to t.Hi1 do t[k] := Value[k,i2];
  end else
  if i2 = _ then
  begin
    t := TIArr1D.Create(Lo2,Hi2);
    for k := t.Lo1 to t.Hi1 do t[k] := Value[i1,k];
  end else
    Raise ERangeError.Create(RS_NoBlnkIndices);
  result := t;
end;

procedure TIArr2D.SetSlice(i1, i2: TInt; const aValue: IIArr1D);
var
  k: TInt;
begin
  if i1 = _ then
  begin
    if i2 = _ then
      Raise ERangeError.Create(RS_TwoBlnkIndices);
    if (aValue.Lo1<>Lo1) or (aValue.Hi1 <> Hi1) then
      Raise ERangeError.Create(RS_LimMismatch);
    for k := aValue.Lo1 to aValue.Hi1 do Value[k,i2] := aValue[k];
  end else
  if i2 = _ then
  begin
    if (aValue.Lo1<>Lo2) or (aValue.Hi1 <> Hi2) then
      Raise ERangeError.Create(RS_LimMismatch);
    for k := aValue.Lo1 to aValue.Hi1 do Value[i1,k] := aValue[k];
  end else
      Raise ERangeError.Create(RS_NoBlnkIndices);
end;

procedure TIArr2D.Swap(i, j: TInt; st: TSliceType);
var
  k, v: TInt;
begin
  case st of
  _Row:for k := Lo2 to Hi2 do
        begin
          v := Value[i,k];
          Value[i,k] := Value[j,k];
          Value[j,k] := v;
        end;
  _Col:for k := Lo1 to Hi1 do
        begin
          v := Value[k,i];
          Value[k,i] := Value[k,j];
          Value[k,j] := v;
        end;
  else
//    Raise ERangeError.CreateFmt(RS_DimError,[idx, 2]);
  end;
end;

{ TIArr3D }

procedure TIArr3D.Assign(const A: IIArr3D);
var
  i1, i2, i3: TInt;
begin
  Assert((Hi1=A.Hi1) and (Lo1=A.Lo1) and (Hi2=A.Hi2) and (Lo2=A.Lo2) and
         (Hi3=A.Hi3) and (Lo3=A.Lo3),RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
    for i2 := A.Lo2 to A.Hi2 do
      for i3 := A.Lo3 to A.Hi3 do
        Value[i1,i2,i3] := A[i1,i2,i3];
end;

constructor TIArr3D.Create(aLo1, aHi1, aLo2, aHi2, aLo3, aHi3: TInt);
begin
  inherited;
  SetLength(fValue, fHi3-fLo3+1, fHi2-fLo2+1, fHi1-fLo1+1);
end;

constructor TIArr3D.Create(aHi1, aHi2, aHi3: TInt);
begin
  Create(1, aHi1, 1, aHi2, 1, aHi3);
end;

constructor TIArr3D.Create(A: IIArr3D; CopyData: boolean = false);
begin
  Create(A.Lo1, A.Hi1, A.Lo2, A.Hi2, A.Lo3, A.Hi3);
  if CopyData then Assign(A);
end;

destructor TIArr3D.Destroy;
begin
  SetLength(fValue, 0, 0, 0);
  fValue := nil;
  inherited;
end;

procedure TIArr3D.Fill(C: TInt);
var
  i1, i2, i3: TInt;
begin
  for i1 := Lo1 to Hi1 do
    for i2 := Lo2 to Hi2 do
      for i3 := Lo3 to Hi3 do
        Value[i1,i2,i3] := C;
end;

function TIArr3D.GetValue(i1, i2, i3: TInt): TInt;
begin
  Assert((i1 >= fLo1) and (i1 <= fHi1) and (i2 >= fLo2) and (i2 <= fHi2) and
         (i3 >= fLo3) and (i3 <= fHi3), RS_IndexOutOfRange);
  result := fValue[i3-fLo3,i2-fLo2,i1-fLo1];
end;

procedure TIArr3D.SetValue(i1, i2, i3: TInt; const Value: TInt);
begin
  Assert((i1 >= fLo1) and (i1 <= fHi1) and (i2 >= fLo2) and (i2 <= fHi2) and
         (i3 >= fLo3) and (i3 <= fHi3), RS_IndexOutOfRange);
  fValue[i3-fLo3,i2-fLo2,i1-fLo1] := Value;
end;

{ TFArr1D }

procedure TFArr1D.Assign(const A: IFArr1D);
var
  i1: TInt;
begin
  Assert( (Hi1=A.Hi1) and (Lo1=A.Lo1), RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
      Value[i1] := A[i1];
end;

constructor TFArr1D.Create(aLo1, aHi1: TInt);
begin
  inherited;
  fInternalValue := true;
  SetLength(fValue, fHi1-fLo1+1);
end;

constructor TFArr1D.Create(aHi1: TInt);
begin
  Create(1, aHi1);
end;

constructor TFArr1D.Create(A: IFArr1D; CopyData: boolean = false);
begin
  Create(A.Lo1, A.Hi1);
  if CopyData then Assign(A);
end;

constructor TFArr1D.Create(aLo1, aHi1: TInt; var  A: TFArr);
var
  d: TInt;
begin
  fInternalValue := false;
  if aHi1 = _ then
  begin
    d := Length(A);
    inherited Create(aLo1, aLo1+d-1);
    fValue := A;
  end else
  begin
    inherited Create(aLo1, aHi1);
    SetLength(A, fHi1-fLo1+1);
    fValue := A;
  end;
end;

destructor TFArr1D.Destroy;
begin
  if fInternalValue then
  begin
    SetLength(fValue, 0);
    fValue := nil;
  end;
  inherited;
end;

function TFArr1D.Dot(const A: IFArr1D): TFloat;
var
  i1: TInt;
  sum: TFloat;
begin
  if (Hi1<>A.Hi1) then
    Raise ERangeError.Create(RS_LimMismatch);
  sum := 0.0;
  for i1 := A.Lo1 to A.Hi1 do
    sum := sum + Value[i1]*A[i1];
  result := sum;
end;

procedure TFArr1D.Fill(C: TFloat);
var
  i1: TInt;
begin
  for i1 := Lo1 to Hi1 do
    Value[i1] := C;
end;

function TFArr1D.GetValue(i1: TInt): TFloat;
begin
  Assert((i1>=fLo1) and (i1<=fHi1), RS_IndexOutOfRange);
  result := fValue[i1-fLo1];
end;

function TFArr1D.Norm(Normalize: boolean): TFloat;
var
  i1: TInt;
  sum: TFloat;
begin
  sum := 0.0;
  for i1 := Lo1 to Hi1 do
    sum := sum + sqr(Value[i1]);
  result := sqrt(sum);
  if Normalize then
  begin
    if result > SqrtMinFloat then
      for i1 := Lo1 to Hi1 do
        Value[i1] := Value[i1]/result
    else
      Raise EMathError.Create(RS_ZeroVector);
  end;
end;

function TFArr1D.MaxAbs: TFloat;
var
  i1: TInt;
  t: TFloat;
begin
  result := 0.0;
  for i1 := Lo1 to Hi1 do
  begin
    t := abs(Value[i1]);
    if t > result then
     result := t;
  end;
end;

procedure TFArr1D.SetValue(i1: TInt; const Value: TFloat);
begin
  If (i1 < fLo1) or (i1 > fHi1) then
     Raise ERangeError.CreateFmt(RS_IndexError,['',i1, fLo1, fHi1]);
  fValue[i1-fLo1] := Value;
end;

procedure TFArr1D.Swap(i, j: TInt);
var
  v: TFloat;
begin
  v := Value[i];
  Value[i] := Value[j];
  Value[j] := v
end;

procedure TFArr1D.Times(C: TFloat);
var
  i1: TInt;
begin
  for i1 := Lo1 to Hi1 do
    Value[i1] := Value[i1]*C;
end;

{ TFArr2D }

procedure TFArr2D.Assign(const A: IFArr2D);
var
  i1, i2: TInt;
begin
  Assert((Hi1=A.Hi1) and (Lo1=A.Lo1) and (Hi2=A.Hi2) and (Lo2=A.Lo2),RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
    for i2 := A.Lo2 to A.Hi2 do
      Value[i1,i2] := A[i1,i2];
end;

constructor TFArr2D.Create(aLo1, aHi1, aLo2, aHi2: TInt);
begin
  inherited;
  SetLength(fValue, fHi2-fLo2+1, fHi1-fLo1+1);
end;

constructor TFArr2D.Create(aHi1, aHi2: TInt);
begin
  Create(1, aHi1, 1, aHi2);
end;

constructor TFArr2D.Create(A: IFArr2D; CopyData: boolean = false);
begin
  Create(A.Lo1, A.Hi1, A.Lo2, A.Hi2);
  if CopyData then Assign(A);
end;

destructor TFArr2D.Destroy;
begin
  SetLength(fValue, 0, 0);
  fValue := nil;
  inherited;
end;

procedure TFArr2D.Fill(C: TFloat);
var
  i1, i2: TInt;
begin
  for i1 := Lo1 to Hi1 do
    for i2 := Lo2 to Hi2 do
        Value[i1,i2] := C;
end;

function TFArr2D.GetValue(i1, i2: TInt): TFloat;
begin
  Assert((i1>=fLo1) and (i1<=fHi1) and (i2>=fLo2) and (i2<=fHi2),RS_IndexOutOfRange);
  result := fValue[i2-fLo2,i1-fLo1];
end;

procedure TFArr2D.SetValue(i1, i2: TInt; const Value: TFloat);
begin
  Assert((i1>=fLo1) and (i1<=fHi1) and (i2>=fLo2) and (i2<=fHi2),RS_IndexOutOfRange);
  fValue[i2-fLo2,i1-fLo1] := Value;
end;

{ M.Norm(i1,_) = Norm2Row(M, i1)
  M.Norm(_,i2) = NormC2ol(M, i2)
  M.Norm(_,_) = Norm2(M)         }
function TFArr2D.Norm2(i1, i2: TInt): TFloat;
var
  k, l: TInt;
  s: TFloat;
begin
   s := 0.0;
  if i1 = _ then
  begin
    if i2 = _ then
    begin  // Full norm
      for k := fLo1 to fHi1 do
        for l := fLo2 to fHi2 do
          s := s + sqr(Value[k,l]);
    end else
    begin // i2-th col
      for k := fLo1 to fHi1 do
        s := s + sqr(Value[k,i2]);
    end;
  end else
  if i2 = _ then
  begin  // i1-th row
    for l := fLo2 to fHi2 do
      s := s + sqr(Value[i1,l]);
  end else
    Raise ERangeError.Create(RS_NoBlnkIndices);
  result := s;
end;

{ M.Norm(i1,_) = NormRow(M, i1)
  M.Norm(_,i2) = NormCol(M, i2)
  M.Norm(_,_) = Norm(M)     }
function TFArr2D.Norm(i1, i2: TInt): TFloat;
begin
  result := sqrt(Norm2(i1, i2))
end;

{ M.GetSlice(i1,_) = GetRow(M, i1)
  M.GetSlice(_,i2) = GetCol(M, i2) }
function TFArr2D.GetSlice(i1, i2: TInt): IFArr1D;
var
  t: IFArr1D;
  k: TInt;
begin
  if i1 = _ then
  begin
    if i2 = _ then
      Raise ERangeError.Create(RS_TwoBlnkIndices);
    t := TFArr1D.Create(Lo1,Hi1);  // col
    for k := t.Lo1 to t.Hi1 do t[k] := Value[k,i2];
  end else
  if i2 = _ then
  begin
    t := TFArr1D.Create(Lo2,Hi2);  // row
    for k := t.Lo1 to t.Hi1 do t[k] := Value[i1,k];
  end else
    Raise ERangeError.Create(RS_NoBlnkIndices);
  result := t;
end;

procedure TFArr2D.SetSlice(i1, i2: TInt; const aValue: IFArr1D);
var
  k: TInt;
begin
  if i1 = _ then
  begin
    if i2 = _ then
      Raise ERangeError.Create(RS_TwoBlnkIndices);
    if (aValue.Lo1<>Lo1) or (aValue.Hi1 <> Hi1) then
      Raise ERangeError.Create(RS_LimMismatch);
    for k := aValue.Lo1 to aValue.Hi1 do Value[k,i2] := aValue[k];
  end else
  if i2 = _ then
  begin
    if (aValue.Lo1<>Lo2) or (aValue.Hi1 <> Hi2) then
      Raise ERangeError.Create(RS_LimMismatch);
    for k := aValue.Lo1 to aValue.Hi1 do Value[i1,k] := aValue[k];
  end else
      Raise ERangeError.Create(RS_NoBlnkIndices);
end;

procedure TFArr2D.Swap(i, j: TInt; st: TSliceType);
var
  k: TInt;
  v: TFloat;
begin
  case st of
  _Row:for k := Lo2 to Hi2 do
        begin
          v := Value[i,k];
          Value[i,k] := Value[j,k];
          Value[j,k] := v;
        end;
  _Col:for k := Lo1 to Hi1 do
        begin
          v := Value[k,i];
          Value[k,i] := Value[k,j];
          Value[k,j] := v;
        end;
  else
//    Raise ERangeError.CreateFmt(RS_DimError,[idx, 2]);
  end;
end;

procedure TFArr2D.Times(C: TFloat);
var
  i1, i2: TInt;
begin
  for i1 := Lo1 to Hi1 do
    for i2 := Lo2 to Hi2 do
      Value[i1,i2] := Value[i1,i2]*C;
end;

{ TFArr3D }

procedure TFArr3D.Assign(const A: IFArr3D);
var
  i1, i2, i3: TInt;
begin
  Assert( (Hi1=A.Hi1) and (Lo1=A.Lo1) and (Hi2=A.Hi2) and (Lo2=A.Lo2) and
          (Hi3=A.Hi3) and (Lo3=A.Lo3), RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
    for i2 := A.Lo2 to A.Hi2 do
      for i3 := A.Lo3 to A.Hi3 do
        Value[i1,i2,i3] := A[i1,i2,i3];
end;

constructor TFArr3D.Create(aLo1, aHi1, aLo2, aHi2, aLo3, aHi3: TInt);
begin
  inherited;
  SetLength(fValue, fHi3-fLo3+1, fHi2-fLo2+1, fHi1-fLo1+1);
end;

constructor TFArr3D.Create(aHi1, aHi2, aHi3: TInt);
begin
  Create(1, aHi1, 1, aHi2, 1, aHi3);
end;

constructor TFArr3D.Create(A: IFArr3D; CopyData: boolean = false);
begin
  Create(A.Lo1, A.Hi1, A.Lo2, A.Hi2, A.Lo3, A.Hi3);
  if CopyData then Assign(A);
end;

destructor TFArr3D.Destroy;
begin
  SetLength(fValue, 0, 0, 0);
  fValue := nil;
  inherited;
end;

procedure TFArr3D.Fill(C: TFloat);
var
  i1, i2, i3: TInt;
begin
  for i1 := Lo1 to Hi1 do
    for i2 := Lo2 to Hi2 do
      for i3 := Lo3 to Hi3 do
        Value[i1,i2,i3] := C;
end;

function TFArr3D.GetValue(i1, i2, i3: TInt): TFloat;
begin
  Assert((i1 >= fLo1) and (i1 <= fHi1) and (i2 >= fLo2) and (i2 <= fHi2) and
         (i3 >= fLo3) and (i3 <= fHi3), RS_IndexOutOfRange);
  result := fValue[i3-fLo3,i2-fLo2,i1-fLo1];
end;

procedure TFArr3D.SetValue(i1, i2, i3: TInt; const Value: TFloat);
begin
  Assert((i1 >= fLo1) and (i1 <= fHi1) and (i2 >= fLo2) and (i2 <= fHi2) and
         (i3 >= fLo3) and (i3 <= fHi3), RS_IndexOutOfRange);
  fValue[i3-fLo3,i2-fLo2,i1-fLo1] := Value;
end;

procedure TFArr3D.Times(C: TFloat);
var
  i1, i2, i3: TInt;
begin
  for i1 := Lo1 to Hi1 do
    for i2 := Lo2 to Hi2 do
      for i3 := Lo3 to Hi3 do
        Value[i1,i2,i3] := Value[i1,i2,i3]*C;
end;

{ TCArr1D }

constructor TCArr1D.Create(aLo1, aHi1: TInt);
begin
  inherited Create;
  fRe := TFArr1D.Create(aLo1, aHi1);
  fIm := TFArr1D.Create(aLo1, aHi1);
end;

constructor TCArr1D.Create(aHi1: TInt);
begin
  Create(1, aHi1);
end;

constructor TCArr1D.Create(A: ICArr1D; CopyData: boolean);
begin
  Create(A.Lo1, A.Hi1);
  if CopyData then Assign(A);
end;

destructor TCArr1D.Destroy;
begin
  fRe := nil;
  fIm := nil;
  inherited;
end;

procedure TCArr1D.Assign(const A: ICArr1D);
begin
  fRe.Assign(A.Re);
  fIm.Assign(A.Im);
end;

procedure TCArr1D.Fill(C: Complex);
begin
  fRe.Fill(C.Re);
  fIm.Fill(C.Im);
end;

function TCArr1D.GetIm: IFArr1D;
begin
  result := fIm;
end;

function TCArr1D.GetRe: IFArr1D;
begin
  result := fRe;
end;

function TCArr1D.GetValue(i1: TInt): Complex;
begin
  result.Re := fRe[i1];
  result.Im := fIm[i1];
end;

procedure TCArr1D.SetIm(const Value: IFArr1D);
begin
  fIm := Value;
end;

procedure TCArr1D.SetRe(const Value: IFArr1D);
begin
  fRe := Value;
end;

procedure TCArr1D.SetValue(i1: TInt; const Value: Complex);
begin
  fRe[i1] := Value.Re;
  fIm[i1] := Value.Im;
end;

procedure TCArr1D.Swap(i, j: TInt);
begin
  fRe.Swap(i, j);
  fIm.Swap(i, j);
end;

procedure TCArr1D.Times(C: Complex);
var
  i1: TInt;
  x: TFloat;
begin
  for i1 := Lo1 to Hi1 do
    begin
        x := fRe[i1]*C.Re - fIm[i1]*C.Im;
        fIm[i1] := fRe[i1]*C.Im + fIm[i1]*C.Re;
        fRe[i1] := x;
    end;
end;

procedure TCArr1D.Conjugate;
begin
  Im.Times(-1.0);
end;

function TCArr1D.GetDim1: TInt;
begin
  result := fRe.Dim1;
end;

function TCArr1D.GetHi1: TInt;
begin
  result := fRe.Hi1;
end;

function TCArr1D.GetLo1: TInt;
begin
  result := fRe.Lo1;
end;

procedure TCArr1D.SetLo1(const Value: TInt);
begin
  fRe.Lo1 := Value;
  fIm.Lo1 := Value;
end;

function TCArr1D.MaxAbs: TFloat;
var
  i1: TInt;
  t: TFloat;
begin
  result := 0.0;
  for i1 := Lo1 to Hi1 do
  begin
    t := abs(fRe[i1])+abs(fIm[i1]);
    if t > result then
     result := t;
  end;
end;

{ TCArr2D }

procedure TCArr2D.Assign(const A: ICArr2D);
var
  i1, i2: TInt;
begin
  Assert((Hi1=A.Hi1) and (Lo1=A.Lo1) and (Hi2=A.Hi2) and (Lo2=A.Lo2),RS_LimMismatch);
  for i1 := A.Lo1 to A.Hi1 do
    for i2 := A.Lo2 to A.Hi2 do
      Value[i1,i2] := A[i1,i2];
end;

constructor TCArr2D.Create(aLo1, aHi1, aLo2, aHi2: TInt);
begin
  inherited Create;
  fRe := TFArr2D.Create(aLo1, aHi1, aLo2, aHi2);
  fIm := TFArr2D.Create(aLo1, aHi1, aLo2, aHi2);
end;

constructor TCArr2D.Create(aHi1, aHi2: TInt);
begin
  Create(1, aHi1, 1, aHi2);
end;

constructor TCArr2D.Create(A: ICArr2D; CopyData: boolean);
begin
  Create(A.Lo1, A.Hi1, A.Lo2, A.Hi2);
  if CopyData then Assign(A);
end;

destructor TCArr2D.Destroy;
begin
  fRe := nil;
  fIm := nil;
  inherited;
end;

procedure TCArr2D.Fill(C: Complex);
var
  i1, i2: TInt;
begin
  for i1 := Lo1 to Hi1 do
    for i2 := Lo2 to Hi2 do
    begin
      fRe[i1,i2] := C.Re;
      fIm[i1,i2] := C.Im;
    end;
end;

function TCArr2D.GetIm: IFArr2D;
begin
  result := fIm;
end;

function TCArr2D.GetRe: IFArr2D;
begin
  result := fRe;
end;

function TCArr2D.GetValue(i1, i2: TInt): Complex;
begin
  result.Re := fRe[i1, i2];
  result.Im := fIm[i1, i2];
end;

procedure TCArr2D.SetIm(const Value: IFArr2D);
begin
  fIm := Value;
end;

procedure TCArr2D.SetRe(const Value: IFArr2D);
begin
  fRe := Value;
end;

procedure TCArr2D.SetValue(i1, i2: TInt; const Value: Complex);
begin
  fRe[i1, i2] := Value.Re;
  fIm[i1, i2] := Value.Im;
end;

procedure TCArr2D.Swap(i, j: TInt; st: TSliceType);
begin
  fRe.Swap(i, j, st);
  fIm.Swap(i, j, st);
end;

procedure TCArr2D.Times(C: Complex);
var
  i1, i2: TInt;
  x: TFloat;
begin
  for i1 := Lo1 to Hi1 do
    for i2 := Lo2 to Hi2 do
    begin
        x := fRe[i1, i2]*C.Re - fIm[i1, i2]*C.Im;
        fIm[i1, i2] := fRe[i1, i2]*C.Im + fIm[i1, i2]*C.Re;
        fRe[i1, i2] := x;
    end;
end;

procedure TCArr2D.Conjugate;
begin
  Im.Times(-1.0);
end;

function TCArr2D.GetDim1: TInt;
begin
  result := fRe.Dim1;
end;

function TCArr2D.GetDim2: TInt;
begin
  result := fRe.Dim2;
end;

function TCArr2D.GetHi1: TInt;
begin
  result := fRe.Hi1;
end;

function TCArr2D.GetHi2: TInt;
begin
  result := fRe.Hi2;
end;

function TCArr2D.GetLo1: TInt;
begin
  result := fRe.Lo1;
end;

function TCArr2D.GetLo2: TInt;
begin
  result := fRe.Lo2;
end;

procedure TCArr2D.SetLo1(const Value: TInt);
begin
  fRe.Lo1 := Value;
  fIm.Lo1 := Value;
end;

procedure TCArr2D.SetLo2(const Value: TInt);
begin
  fRe.Lo2 := Value;
  fIm.Lo2 := Value;
end;

{ TEigenSys }

constructor TEigenSys.Create(aLo1, aHi1: TInt);
begin
  inherited Create;
  fVectors := TFArr2D.Create(aLo1, aHi1, aLo1, aHi1);
  fValues := TFArr1D.Create(aLo1, aHi1);
  fNames := TSArr1D.Create(aLo1, aHi1);
end;

constructor TEigenSys.Create(aHi1: TInt);
begin
  Create(1, aHi1);
end;

constructor TEigenSys.Create(const a: ILim2D);
begin
  if not IsSquare(a) then
    raise ENonSquareMatrix.Create('Matrix is not square');
  Create(a.Lo1, a.Hi1);
end;

destructor TEigenSys.Destroy;
begin
  fVectors := nil;
  fValues := nil;
  fNames := nil;
  inherited;
end;

function TEigenSys.GetNames: ISArr1D;
begin
  result := fNames;
end;

function TEigenSys.GetValues: IFArr1D;
begin
  result := fValues;
end;

function TEigenSys.GetVectors: IFArr2D;
begin
  result := fVectors;
end;

procedure TEigenSys.SetNames(const Value: ISArr1D);
begin
  fNames := Value;
end;

procedure TEigenSys.SetValues(const Value: IFArr1D);
begin
  fValues := Value;
end;

procedure TEigenSys.SetVectors(const Value: IFArr2D);
begin
  fVectors := Value;
end;

{ THEigenSys }

constructor THEigenSys.Create(aLo1, aHi1: TInt);
begin
  inherited Create;
  fVectors := TCArr2D.Create(aLo1, aHi1, aLo1, aHi1);
  fValues := TFArr1D.Create(aLo1, aHi1);
  fNames := TSArr1D.Create(aLo1, aHi1);
end;

constructor THEigenSys.Create(aHi1: TInt);
begin
  Create(1, aHi1);
end;

constructor THEigenSys.Create(const a: ILim2D);
begin
  if not IsSquare(a) then
    raise ENonSquareMatrix.Create('Matrix is not square');
  Create(a.Lo1, a.Hi1);
end;

destructor THEigenSys.Destroy;
begin
  fVectors := nil;
  fValues := nil;
  fNames := nil;
  inherited;
end;

function THEigenSys.GetNames: ISArr1D;
begin
  result := fNames;
end;

function THEigenSys.GetValues: IFArr1D;
begin
  result := fValues;
end;

function THEigenSys.GetVectors: ICArr2D;
begin
  result := fVectors;
end;

procedure THEigenSys.SetNames(const Value: ISArr1D);
begin
  fNames := Value;
end;

procedure THEigenSys.SetValues(const Value: IFArr1D);
begin
  fValues := Value;
end;

procedure THEigenSys.SetVectors(const Value: ICArr2D);
begin
  fVectors := Value;
end;

{ TSVDSys }

constructor TSVDSys.Create(aLo1, aHi1, aLo2, aHi2: TInt; fullmat: boolean);
var
  n, m, nm, nm1, nm2: TInt;
begin
  inherited Create;
  m := aHi1 - aLo1 + 1;
  n := aHi2 - aLo2 + 1;
  nm := min(n,m);
  if fullmat then
  begin
    nm1 := m;
    nm2 := n;
  end else
  begin
    nm1 := nm;
    nm2 := nm;
  end;
  fUvectors := TFArr2D.Create(aLo1, aHi1, 1, nm1);
  fVvectors := TFArr2D.Create(aLo2, aHi2, 1, nm2);
  fValues := TFArr1D.Create(1, nm);
end;

constructor TSVDSys.Create(const A: IFArr2D; fullmat: boolean);
begin
  Create(A.Lo1, A.Hi1, A.Lo2, A.Hi2, fullmat);
end;

constructor TSVDSys.Create(aHi1, aHi2: TInt; fullmat: boolean);
begin
  Create(1, aHi1, 1, aHi2, fullmat);
end;

function TSVDSys.GetUvectors: IFArr2D;
begin
  result := fUvectors;
end;

function TSVDSys.GetValues: IFArr1D;
begin
  result := fValues;
end;

function TSVDSys.GetVvectors: IFArr2D;
begin
  result := fVvectors;
end;

procedure TSVDSys.SetUvectors(const Value: IFArr2D);
begin
  fUvectors := Value;
end;

procedure TSVDSys.SetValues(const Value: IFArr1D);
begin
  fValues := Value;
end;

procedure TSVDSys.SetVvectors(const Value: IFArr2D);
begin
  fVvectors := Value;
end;

end.

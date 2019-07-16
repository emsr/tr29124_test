unit mp_ccalc;

{Parse and evaluate mp_complex expressions}

interface

{$i STD.INC}

{$X+}
{$ifdef BIT16}
  {$N+}
{$endif}


uses
  BTypes, mp_types;

{$i mp_conf.inc}


(*************************************************************************

 DESCRIPTION   :  Parse and evaluate mp_complex expressions

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA :  (mp_types)

 MEMORY USAGE  :  heap

 DISPLAY MODE  :  ---

 REFERENCES    :  [19] T. Norvell: Parsing Expressions by Recursive Descent,
                       http://www.engr.mun.ca/~theo/Misc/exp_parsing.htm
                  [20] G. Toal's tutorial pages OperatorPrecedence.html,CompilersOneOhOne.html,
                       GrahamToalsCompilerDemo.html at http://www.gtoal.com/software/
                  [21] T.R. Nicely's parser.c in factor1.zip from http://www.trnicely.net
                       derived from GMP demo pexpr.c (see [15])

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.0.10   12.11.18  W.Ehrhardt  First BP version derived from mp_rcalc
 0.0.11   12.11.18  we          constant I
 0.0.12   12.11.18  we          _ARG, _CONJ, _RE, _IM
 0.0.13   12.11.18  we          _LOG10
 0.0.14   13.11.18  we          check some functions for invalid argument
 0.0.15   13.11.18  we          _ARCCOTHC
**************************************************************************)



(*-------------------------------------------------------------------------
 (C) Copyright 2008-2018 Wolfgang Ehrhardt

 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from
 the use of this software.

 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it
 freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software in
    a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

 2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

 3. This notice may not be removed or altered from any source distribution.
----------------------------------------------------------------------------*)


{#Z+} {Turn off reference for help}

{Parse errors >0, eval errors < 0}
const
  Err_Missing_LeftBracket  =  1;   {Missing "("}
  Err_Missing_Comma        =  2;   {Missing argument separator}
  Err_Missing_RightBracket =  3;   {Missing ")"}
  Err_Unknown_Function     =  4;   {Unknown function}
  Err_Unknown_Element      =  5;   {Unknown element}
  Err_Trailing_Garbage     =  6;   {Trailing garbage}
  Err_Invalid_Number       =  7;   {Invalid number}
  Err_Unknown_Operation    =  8;   {Unknown operation}
  Err_Invalid_HexNumber    =  9;   {Invalid hex number}

  Err_MPERR_Eval           = -1;   {MP_Err <> MP_OK}
  Err_Division_by_zero     = -2;   {Division by zero}
  Err_Invalid_argument     = -3;   {Invalid arg, e.g. ln(0)}
  Err_X_not_init           = -4;   {Variable X not initialized}
  Err_Y_not_init           = -5;   {Variable Y not initialized}
  Err_Z_not_init           = -6;   {Variable Z not initialized}
  Err_Overflow             = -8;   {Overflow}
{#Z-}


type
  TFOperation = (_CONST, _CHS, _ABS, _ADD, _SUB, _MUL, _DIV, _EXPT, _SQRT,
                 _SQR, _X, _Y, _Z, _I, _ARG, _CONJ, _RE, _IM, _ARCCOSH,
                 _AGM, _ARCCOS, _ARCSIN, _ARCTAN, _ARCTAN2, _ARCSINH, _ARCTANH,
                 _COS, _COSH, _EXP, _EXPM1, _LN, _LN1P, _LOG10,
                 _SIN, _SINH, _TAN, _TANH, _COT, _CSC, _SEC, _COTH, _CSCH, _SECH,
                 _ARCCOT, _ARCCOTC, _ARCCSC, _ARCSEC, _ARCCOTH, _ARCCOTHC,
                 _ARCCSCH, _ARCSECH, _NROOT,
                 _TEST);
                {implemented operators, functions, and variables}

type
  PFExpr  = ^TFExpr;                       {Expression node pointer}
  TFExpr  = record                         {binary tree node}
              op:  TFOperation;            {operation/function/variable}
              nn:  byte;                   {number of nodes}
              case integer of
                0: (Value: mp_complex);    {value if op=_CONST}
                1: (SNode: PFExpr);        {expr = (SNode op) or op(SNode)}
                2: (LNode, RNode: PFExpr;) {expr = LNode op RNode or op(LNode,RNode)}
            end;

  TFEval  = record                         {Evaluation record}
              X  : mp_complex;             {Variable X}
              Y  : mp_complex;             {Variable Y}
              Z  : mp_complex;             {Variable Z}
              Res: mp_complex;             {Evaluation result}
              Err: integer;                {Eval error code}
            end;


function  mpc_parse(psz: pchar8; var e: PFExpr; var Err: integer): pchar8;
  {-Parse string psz into expression tree e, if OK Err=0 and result^=#0,}
  { else Err=Err_xx and result points to error position}

procedure mpc_eval(e: PFExpr; var evr: TFEval);
  {-Evaluate expression tree e, result in evr}

procedure mpc_clear_expr(var e: PFExpr);
  {-Release memory used by e}

procedure mpc_calculate(psz: pchar8; var evr: TFEval; var EPos: integer);
  {-Parse and evaluate string psz}

function  mpc_calc_errorstr(Err: integer): mp_string;
  {-Translate known error codes}

procedure mpc_init_eval(var evr: TFEval);
  {-Initialize the mp_complex records of evr}

procedure mpc_clear_eval(var evr: TFEval);
  {-Clear the mp_complex records of evr}


implementation


uses
  mp_base, mp_numth, mp_real, mp_cmplx;


(* Approx. grammar for expression parser:

  <Expr>    ::=   <Term> '+' <Term>
                | <Term> '-' <Term>
                | '+'<Term>
                | '-'<Term>;
  <Term>    ::=   <Factor> '*' <Factor>
                | <Factor> '/' <Factor>
                | <Factor>;
  <Factor>  ::=   <Element> '^' <Factor>
                | <Element>;
  <Element> ::= <Func> | <Var> | <number> | 'i' | 'pi' | '(' <Expr> ')';
  <Func>    ::= <Ident> '(' <Arglist> ')';
  <Var>     ::= 'X'..'Z' | 'x'..'z';
  <Arglist  ::= <Expr> | <Expr> ',' <Expr>;
  <Ident>   ::= <alpha> {<alpha> | <digit>};
  <intnum>  ::= <digit> { <digit> };
  <expo>    ::=   'e' ['+' | '-'] <intnum>
                | 'E' ['+' | '-'] <intnum>
  <number>  ::= ['+' | '-'] [<intnum>] ['.' [<intnum>] <expo>];
  <digit>   ::= '0'..'9';
  <alpha>   ::= 'A'..'Z'| 'a'..'z';
*)



type
  str10  = string[10];
  TFunc  = record
               op: TFOperation;
             arg2: boolean;
             name: str10;
           end;
type
  TOpVar = _X.._Z;

const
  MaxFun = 41;

const
  FuncTab : array[1..MaxFun] of TFunc = (
             (op: _SQRT      ; arg2: false; name: 'SQRT'),
             (op: _SQR       ; arg2: false; name: 'SQR'),
             (op: _ABS       ; arg2: false; name: 'ABS'),
             (op: _ARG       ; arg2: false; name: 'ARG'),
             (op: _CONJ      ; arg2: false; name: 'CONJ'),
             (op: _RE        ; arg2: false; name: 'RE'),
             (op: _IM        ; arg2: false; name: 'IM'),
             (op: _ARCCOSH   ; arg2: false; name: 'ARCCOSH'),
             (op: _AGM       ; arg2: true ; name: 'AGM'),
             (op: _ARCCOS    ; arg2: false; name: 'ARCCOS'),
             (op: _ARCSIN    ; arg2: false; name: 'ARCSIN'),
             (op: _ARCTAN    ; arg2: false; name: 'ARCTAN'),
             (op: _ARCSINH   ; arg2: false; name: 'ARCSINH'),
             (op: _ARCTANH   ; arg2: false; name: 'ARCTANH'),
             (op: _COS       ; arg2: false; name: 'COS'),
             (op: _COSH      ; arg2: false; name: 'COSH'),
             (op: _EXP       ; arg2: false; name: 'EXP'),
             (op: _EXPM1     ; arg2: false; name: 'EXPM1'),
             (op: _LN        ; arg2: false; name: 'LN'),
             (op: _LN1P      ; arg2: false; name: 'LN1P'),
             (op: _LOG10     ; arg2: false; name: 'LOG10'),
             (op: _SIN       ; arg2: false; name: 'SIN'),
             (op: _SINH      ; arg2: false; name: 'SINH'),
             (op: _TAN       ; arg2: false; name: 'TAN'),
             (op: _TANH      ; arg2: false; name: 'TANH'),
             (op: _COT       ; arg2: false; name: 'COT'),
             (op: _CSC       ; arg2: false; name: 'CSC'),
             (op: _SEC       ; arg2: false; name: 'SEC'),
             (op: _COTH      ; arg2: false; name: 'COTH'),
             (op: _CSCH      ; arg2: false; name: 'CSCH'),
             (op: _SECH      ; arg2: false; name: 'SECH'),
             (op: _ARCCOT    ; arg2: false; name: 'ARCCOT'),
             (op: _ARCCOTC   ; arg2: false; name: 'ARCCOTC'),
             (op: _ARCCSC    ; arg2: false; name: 'ARCCSC'),
             (op: _ARCSEC    ; arg2: false; name: 'ARCSEC'),
             (op: _ARCCOTH   ; arg2: false; name: 'ARCCOTH'),
             (op: _ARCCOTHC  ; arg2: false; name: 'ARCCOTHC'),
             (op: _ARCCSCH   ; arg2: false; name: 'ARCCSCH'),
             (op: _ARCSECH   ; arg2: false; name: 'ARCSECH'),
             (op: _NROOT     ; arg2: true;  name: 'NROOT'),
             (op: _TEST      ; arg2: false; name: 'TEST')  {used for tests, development etc}
            );


{---------------------------------------------------------------------------}
function SkipWhite(psz: pchar8): pchar8;
  {-Skip white space}
begin
  while psz^ in [' ',#13,#10,#9] do inc(psz);
  SkipWhite := psz;
end;


{---------------------------------------------------------------------------}
procedure mkNode(var r: PFExpr; op: TFOperation; nn: byte; e1, e2: PFExpr);
  {-Make a new expression node for (e1 op e2) or (e1 op)}
begin
  r := mp_alloc(sizeof(TFExpr));
  r^.nn := nn;
  r^.op := op;
  r^.LNode := e1;
  r^.RNode := e2;
end;


{---------------------------------------------------------------------------}
procedure mkNode0c(var r: PFExpr);
  {-Make a new node for a constant}
begin
  {alloc Value node initialize mp_complex}
  r := mp_alloc(sizeof(TFExpr));
  r^.op := _CONST;
  r^.nn := 0;
  mpc_init(r^.Value);
end;


{---------------------------------------------------------------------------}
procedure mkNode0v(var r: PFExpr; const v: TOpVar);
  {-Make a new node for a variable v = _X, _Y, or _Z}
begin
  {alloc Value node initialize mp_complex}
  r := mp_alloc(sizeof(TFExpr));
  r^.op := v;
  r^.nn := 0;
end;


{---------------------------------------------------------------------------}
procedure mkNode1(var r: PFExpr; op: TFOperation; e: PFExpr);
  {-Make a new expression node for r := e op}
begin
  mkNode(r,op,1,e,nil);
end;


{---------------------------------------------------------------------------}
procedure mkNode2(var r: PFExpr; op: TFOperation; e1, e2: PFExpr);
  {-Make a new expression node for r := e1 op e2}
begin
  mkNode(r,op,2,e1,e2);
end;


{---------------------------------------------------------------------------}
function GetIdent(psz: pchar8): mp_string;
  {-Gather next identifier, break if not in [A-Z, a-z, 0-9], result is uppercase}
  { first character in must be in [A-Z, a-z]}
var
  s: mp_string;
begin
  s := '';
  if psz^ in ['A'..'Z', 'a'..'z'] then begin
    while psz^ in ['A'..'Z', 'a'..'z', '0'..'9'] do begin
      s := s+upcase(psz^);
      inc(psz);
    end;
  end;
  GetIdent := s;
end;


{---------------------------------------------------------------------------}
function GetFuncIndex(const s: mp_string): integer;
  {-Test if s is a known function name, if yes return index else 0}
var
  i: integer;
begin
  for i:=1 to MaxFun do begin
    if FuncTab[i].name=s then begin
      GetFuncIndex := i;
      exit;
    end;
  end;
  GetFuncIndex := 0;
end;


{---------------------------------------------------------------------------}
function Expr(psz: pchar8; var e: PFExpr; var Err: integer): pchar8; forward;
  {-Parse string psz into expression tree}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function Func(psz: pchar8; idx: integer; var e: PFExpr; var Err: integer): pchar8;
  {-Build function expression, psz points between function name and "("}
var
  e1,e2: PFExpr;
const
  na: array[boolean] of byte = (1,2);

  procedure clear;
    {-Clear local PFExpr if error}
  begin
    if e1<>nil then mpc_clear_expr(e1);
    if e2<>nil then mpc_clear_expr(e2);
  end;

begin
  e1 := nil;
  e2 := nil;
  e  := nil;
  Func := psz;
  if Err<>0 then exit;
  psz := SkipWhite(psz);
  if psz^ <> '(' then begin
    Func := psz;
    Err := Err_Missing_LeftBracket;
    exit;
  end;
  {get first argument}
  psz := Expr(psz+1,e1, Err);
  Func := psz;
  if Err<>0 then begin
    clear;
    exit;
  end;
  psz := SkipWhite(psz);
  if FuncTab[idx].arg2 then begin
    {if two arguments search for ","}
    if psz^ <> mp_arg_sep{','} then begin
      Func := psz;
      Err := Err_Missing_Comma;
      clear;
      exit;
    end;
    {evaluate second argument}
    psz := Expr(psz+1,e2,Err);
    Func := psz;
    if Err<>0 then begin
      clear;
      exit;
    end;
    psz := SkipWhite(psz);
  end
  else e2:=nil;
  {search for closing ")"}
  if psz^ <> ')' then begin
    Func := psz;
    Err := Err_Missing_RightBracket;
    clear;
    exit;
  end;
  inc(psz);
  with FuncTab[idx] do mkNode(e, op, na[arg2], e1, e2);
  Func := psz;
end;


{---------------------------------------------------------------------------}
function Element(psz: pchar8; var e: PFExpr; var Err: integer): pchar8;
  {-Parse an Element}
var
  res: PFExpr;
  s,sc: pchar8;
  lsc: word;
  i: integer;
  s0: char8;
  neg: boolean;
  id: str255;
begin
  Element := psz;
  e := nil;
  if Err<>0 then exit;
  psz := SkipWhite(psz);
  s0 := psz^;
  if s0 in ['A'..'Z', 'a'..'z'] then begin
    id := GetIdent(psz);
    i := GetFuncIndex(id);
    if i=0 then begin
      if id='X' then begin mkNode0v(e, _X); inc(psz); end
      else if id='Y' then begin mkNode0v(e, _Y); inc(psz); end
      else if id='Z' then begin mkNode0v(e, _Z); inc(psz); end
      else if id='I' then begin
        mkNode0c(res);
        mpc_seti(res^.Value);
        inc(psz,1);
        e := res;
      end
      else if id='PI' then begin
        mkNode0c(res);
        with res^.Value do begin
          mpf_set_pi(re);
          mpf_set0(im);
        end;
        inc(psz,2);
        e := res;
      end
      else if id='LN2' then begin
        mkNode0c(res);
        with res^.Value do begin
          mpf_set_ln2(re);
          mpf_set0(im);
        end;
        inc(psz,3);
        e := res;
      end
    {$ifdef MPC_E1Ln10Tab}
      else if id='LN10' then begin
        mkNode0c(res);
        with res^.Value do begin
          mpf_set_ln10(re);
          mpf_set0(im);
        end;
        inc(psz,4);
        e := res;
      end
      else if id='E' then begin
        mkNode0c(res);
        with res^.Value do begin
          mpf_set_exp1(re);
          mpf_set0(im);
        end;
        inc(psz,1);
        e := res;
      end
    {$endif}
      else begin
        Err := Err_Unknown_Function;
        Element := psz;
        exit;
      end;
    end
    else begin
      inc(psz, length(id));
      psz := Func(psz,i,e, Err);
    end;
  end
  else if s0='(' then begin
    psz := Expr(psz+1,e,Err);
    if Err<>0 then exit;
    psz := SkipWhite(psz);
    if psz^<>')' then begin
      Err := Err_Missing_RightBracket;
      Element := psz;
      exit;
    end;
    inc(psz);
  end
  else if s0 in ['0'..'9','+','-',mp_fract_sep{'.'}] then begin
    {get sign}
    neg := false;
    if psz^ in ['+','-'] then begin
      if psz^='-' then neg := true;
      inc(psz);
    end;

    {count decimal characters}
    s := psz;
    {integer part}
    while psz^ in ['0'..'9'] do inc(psz);

    if psz^=mp_fract_sep{'.'} then begin
      inc(psz);
      {fractional part}
      while psz^ in ['0'..'9'] do inc(psz);
    end;

    if (upcase(psz^)='E') and (s<>psz) then begin
       {exponent part}
       inc(psz);
       if psz^ in ['+','-'] then inc(psz);
       while psz^ in ['0'..'9'] do inc(psz);
    end;

    {alloc and move digit string to temp storage}
    if s=psz then begin
      {empty digit string}
      Err := Err_Invalid_Number;
      Element := psz;
      exit;
    end;

    lsc := psz-s+1;
    sc := mp_getmem(lsc);
    move(s^,sc^,lsc-1);
    sc[psz-s] := #0;
    {Make a constant node}
    mkNode0c(res);
    {convert digit string to mp_float}
    mpf_read_decimal(res^.Value.re,sc);
    {apply sign}
    if neg then mpf_chs(res^.Value.re, res^.Value.re);
    mpf_set0(res^.Value.im);
    e := res;
    mp_freemem(pointer(sc),lsc);
  end
  else begin
    Err := Err_Unknown_Element;
    Element := psz;
    exit;
  end;
  Element := psz;
end;


{---------------------------------------------------------------------------}
function Factor(psz: pchar8; var e: PFExpr; var Err: integer): pchar8;
  {-Parse a Factor}
var
  t: PFExpr;
begin
  e := nil;
  if Err<>0 then begin
    Factor := psz;
    exit;
  end;
  psz := Element(psz,e,Err);
  if Err<>0 then begin
    Factor := psz;
    exit;
  end;
  {Look for optional power part of Factor}
  psz := SkipWhite(psz);
  if psz^='^' then begin
    inc(psz);
    t := nil;
    psz := Factor(psz, t, Err);
    if Err=0 then mkNode2(e, _EXPT, e, t)
    else if t<>nil then mpc_clear_expr(t);
  end;
  Factor := psz;
end;


{---------------------------------------------------------------------------}
function Term(psz: pchar8; var e: PFExpr; var Err: integer): pchar8;
  {-Parse a Term}
var
  t: PFExpr;
begin
  e := nil;
  if Err<>0 then begin
    Term := psz;
    exit;
  end;
  t := nil;
  psz := Factor(psz,e, Err);
  while Err=0 do begin
    psz := SkipWhite(psz);
    case upcase(psz[0]) of
      '*' : begin
              psz := Factor(psz+1, t, Err);
              if Err=0 then mkNode2(e, _MUL, e, t);
            end;
      '/':  begin
              psz := Factor(psz+1, t, Err);
              if Err=0 then mkNode2(e, _DIV, e, t);
            end;
     else   break;
    end;
  end;
  if (Err<>0) and (t<>nil) then mpc_clear_expr(t);
  Term := psz;
end;


{---------------------------------------------------------------------------}
function Expr(psz: pchar8; var e: PFExpr; var Err: integer): pchar8;
  {-Parse an Expr}
var
  t: PFExpr;
  c: char8;
begin
  e := nil;
  if Err<>0 then begin
    Expr := psz;
    exit;
  end;
  psz := SkipWhite(psz);
  c  := psz^;
  {Unary prefix operators}
  if c='+' then begin
    {just skip the +}
    psz := Term(psz+1, e, Err)
  end
  else if c='-' then begin
    psz := Term(psz+1, e, Err);
    if Err<>0 then begin
      Expr := psz;
      exit;
    end;
    mkNode1(e, _CHS, e);
  end
  else psz := Term(psz, e, Err);

  t := nil;
  while Err=0 do begin
    psz := SkipWhite(psz);
    case psz^ of
      '+': begin
             psz := Term(psz+1, t, Err);
             if Err=0 then mkNode2(e, _ADD, e, t);
           end;
      '-': begin
             psz := Term(psz+1, t, Err);
             if Err=0 then mkNode2(e, _SUB, e, t);
           end;
      else break;
    end; {case}
  end;
  if (Err<>0) and (t<>nil) then mpc_clear_expr(t);
  Expr := psz;
end;


{---------------------------------------------------------------------------}
function mpc_parse(psz: pchar8; var e: PFExpr; var Err: integer): pchar8;
  {-Parse string psz into expression tree e, Err=0 and psz^=#0 if OK,}
  { else Err=Err_xx and result points to error postions}
begin
  Err := 0;
  psz := Expr(psz, e, Err);
  if (Err=0) and (psz^<>#0) then Err := Err_Trailing_Garbage;
  mpc_parse := psz;
end;



{---------------------------------------------------------------------------}
{------  The next functions are essentially copies from AMath  -------------}
{---------------------------------------------------------------------------}

const
  x80000000 = longint($80000000);

type
  THexDblA = packed array[0..7] of byte;  {Double   as array of bytes}
  THexSglA = packed array[0..3] of byte;  {Single   as array of bytes}
  THexExtW = packed array[0..4] of word;  {Extended as array of word}
  THexDblW = packed array[0..3] of word;  {Double   as array of word}
  THexSglW = packed array[0..1] of word;  {Single   as array of word}

type
  TDblRec  = packed record     {Double as sign, exponent, significand}
               lm: longint;    {low  32 bit of significand}
               hm: longint;    {high bits of significand, biased exponent and sign}
             end;

const MaxDblHex   : THexDblW = ($ffff,$ffff,$ffff,$7fef); {1.797693134862315E+308}
const MaxSglHex   : THexSglA = ($ff,$ff,$7f,$7f);  {3.4028234E+38}


{---------------------------------------------------------------------------}
function predd(d: double): double;
  {-Return next representable double after d in the direction -Inf}
begin
  with TDblRec(d) do begin
    if THexDblW(d)[3] and $7FF0=$7FF0 then begin
      {Inf or Nan}
      if (hm and $7FFFFFFF=$7FF00000) and (lm=0) then begin
        {d is +- Inf}
        if d>0.0 then d := double(MaxDblHex);
      end;
    end
    else begin
      {finite number}
      if d=0.0 then begin
        hm := x80000000;
        lm := 1;
      end
      else if d<0.0 then begin
        {d<0: increment significand}
        inc(lm);
        if lm=0 then inc(hm);
      end
      else begin
        {d>0: decrement significand}
        if lm=0 then dec(hm);
        dec(lm);
      end;
    end;
    predd := d;
  end;
end;


{---------------------------------------------------------------------------}
function preds(s: single): single;
  {-Return next representable single after s in the direction -Inf}
var
  L: longint absolute s;
begin
  if L and $7F800000 = $7F800000 then begin
    {Inf or Nan, don't change Nan or -Inf}
    if L and $7FFFFF = 0 then begin
      {s is +- Inf}
      if s>0.0 then s := single(MaxSglHex);
    end;
  end
  else begin
    {finite number}
    if s=0.0 then L := x80000000
    else if s<0.0 then inc(L)
    else dec(L);
  end;
  preds := s;
end;


{---------------------------------------------------------------------------}
function succd(d: double): double;
  {-Return next representable double after d in the direction +Inf}
begin
  with TDblRec(d) do begin
    if THexDblW(d)[3] and $7FF0=$7FF0 then begin
      {Inf or Nan}
      if (hm and $7FFFFFFF=$7FF00000) and (lm=0) then begin
        {d is +- Inf}
        if d<0.0 then d := -double(MaxDblHex);
      end;
    end
    else begin
      {finite number}
      if d=0.0 then begin
        hm := 0;
        lm := 1;
      end
      else if d>0.0 then begin
        {d>0: increment significand}
        inc(lm);
        {hm < $7FF00000, so inc(hm) cannot overflow and will give}
        {the correct result succd(predd(PosInf_d)) = PosInf_d}
        if lm=0 then inc(hm);
      end
      else begin
        {d<0: decrement significand}
        if lm=0 then dec(hm);
        dec(lm);
      end;
    end;
    succd := d;
  end;
end;


{---------------------------------------------------------------------------}
function succs(s: single): single;
  {-Return next representable single after s in the direction +Inf}
var
  L: longint absolute s;
begin
  if L and $7F800000 = $7F800000 then begin
    {Inf or Nan, don't change Nan or +Inf}
    if L and $7FFFFF = 0 then begin
      {s is +- Inf}
      if s<0.0 then s := -single(MaxSglHex);
    end;
  end
  else begin
    {finite number}
    if s=0.0 then L := 1
    else if s>0.0 then inc(L)
    else dec(L);
  end;
  succs := s;
end;


{$ifndef EXT64}

{Only implemented for hardware 80-bit extended}
type
  TExtRec  = packed record     {Extended as sign, exponent, significand}
               lm: longint;    {low  32 bit of significand}
               hm: longint;    {high 32 bit of significand}
               xp: word;       {biased exponent and sign  }
             end;

const MinExtHex   : THexExtW = ($0000,$0000,$0000,$8000,$0001); {MinExtended as Hex}
const MaxExtHex   : THexExtW = ($ffff,$ffff,$ffff,$ffff,$7ffe); {MaxExtended as Hex}


{---------------------------------------------------------------------------}
function predx(x: extended): extended;
  {-Return next representable extended after x in the direction -Inf}
begin
  with TExtRec(x) do begin
    if xp and $7FFF=$7FFF then begin
      {Inf or Nan}
      if (hm=x80000000) and (lm=0) then begin
        {x is +- Inf}
        if x>0.0 then x := extended(MaxExtHex);
      end;
      predx := x;
      exit;
    end;

    if xp and $7FFF = 0 then begin
      {Handle pseudo-denormal: Set exponent to +/- 1, significand is unchanged}
      if hm<0 then xp := xp or 1;
    end
    else if hm>=0 then begin
      {don't touch unnormals}
      predx := x;
      exit;
    end;

    {finite number}
    if x=0 then begin
      xp := $8000;
      lm := 1;
    end
    else if xp and $8000 <> 0 then begin
      {x<0: increment significand}
      inc(lm);
      if lm=0 then begin
        inc(hm);
        if (hm=0) or ((xp=$8000) and (hm=x80000000))  then begin
          inc(xp);
          hm := hm or x80000000;
          if xp=$FFFF then x := -DblNegInf;
        end;
      end;
    end
    else begin
      {x>0: decrement significand}
      if hm<0 then begin
        if lm=0 then begin
          if hm=x80000000 then begin
            dec(xp);
            dec(hm);
            if xp>0 then hm := hm or x80000000;
          end
          else dec(hm);
        end;
        dec(lm);
      end
      else begin
        {denormal}
        if lm=0 then dec(hm);
        dec(lm);
      end;
    end;
  end;
  predx := x;
end;


{---------------------------------------------------------------------------}
function succx(x: extended): extended;
  {-Return next representable extended after x in the direction +Inf}
begin
  with TExtRec(x) do begin

    if xp and $7FFF=$7FFF then begin
      {Inf or Nan}
      if (hm=x80000000) and (lm=0) then begin
        {x is +- Inf}
        if x<0.0 then x := -extended(MaxExtHex);
      end;
      succx := x;
      exit;
    end;
    if xp and $7FFF = 0 then begin
      {Handle pseudo-denormal: Set exponent to +/- 1, significand is unchanged}
      if hm<0 then xp := xp or 1;
    end
    else if hm>=0 then begin
      {don't touch unnormals}
      succx := x;
      exit;
    end;

    {finite number}
    if x=0.0 then begin
      xp := 0;
      lm := 1;
    end
    else if xp and $8000 = 0 then begin
      {x>0: increment significand}
      inc(lm);
      if lm=0 then begin
        inc(hm);
        if (hm=0) or ((xp=0) and (hm=x80000000)) then begin
          inc(xp);
          hm := hm or x80000000;
          if xp=$7FFF then x :=  DblPosInf;;
        end;
      end;
    end
    else begin
      {x<0: decrement significand}
      if lm=0 then begin
        if (hm>=0) or (hm=x80000000)  then begin
          dec(hm);
          dec(xp);
          if xp and $7FFF > 0 then hm := hm or x80000000;
        end
        else dec(hm);
      end;
      dec(lm);
    end;
  end;
  succx := x;
end;

{$endif}


{---------------------------------------------------------------------------}
procedure eval(e: PFExpr; var r: mp_complex; var evr: TFEval);
  {-(internal) evaluate expression tree e, result in r}
var
  v1,v2: mp_complex;
  lr,li,lx: longint;
  td: double;
begin
  if evr.Err<>0 then exit;

  if e^.nn=0 then begin
    case e^.op of
      _CONST: begin
                mpc_copy(e^.Value,r);
                if MP_Error<>0 then evr.Err := Err_MPERR_Eval;
              end;
          _X: begin
                if mpc_not_init(evr.X) then evr.Err := Err_X_not_init
                else mpc_copy(evr.X,r);
              end;
          _Y: begin
                if mpc_not_init(evr.Y) then evr.Err := Err_Y_not_init
                else mpc_copy(evr.Y,r);
              end;
          _Z: begin
                if mpc_not_init(evr.Z) then evr.Err := Err_Z_not_init
                else mpc_copy(evr.Z,r);
              end;
         else evr.Err := Err_Unknown_Operation;
    end;
    exit;
  end;

  {always initialize two mp_complex although some ops need only 0 or 1}
  mpc_initp2(v1,v2,mpf_get_default_prec);
  if MP_Error<>MP_OKAY then begin
    evr.Err := Err_MPERR_Eval;
    exit;
  end;

  if e^.nn=1 then begin
    eval(e^.LNode, v1, evr);
    if evr.Err=0 then begin
      lr := s_mpf_ldx(v1.re);
      li := s_mpf_ldx(v1.im);
      case e^.op of
            _CHS:  mpc_chs(v1,r);
            _ABS:  begin
                     mpc_abs(v1,r.re);
                     mpf_set0(r.im);
                   end;
            _ARG:  begin
                     mpc_arg(v1,r.re);
                     mpf_set0(r.im);
                   end;
             _RE:  begin
                     mpf_copyp(v1.re,r.re);
                     mpf_set0(r.im);
                   end;
             _IM:  begin
                     mpf_copyp(v1.im,r.re);
                     mpf_set0(r.im);
                   end;
           _CONJ:  begin
                     mpc_conj(v1,r);
                   end;
            _SQR:  begin
                     if (lr<MaxLongint) and (li<MaxLongint) then mpc_sqr(v1,r)
                     else evr.Err := Err_Overflow;
                   end;
           _SQRT:  begin
                     mpc_sqrt(v1,r)
                   end;

        _ARCCOSH:  begin
                     mpc_arccosh(v1,r);
                   end;

         _ARCCOS:  begin
                     mpc_arccos(v1,r);
                   end;

         _ARCSIN:  begin
                     mpc_arcsin(v1,r);
                   end;

         _ARCTAN:  begin
                     if mpc_is_ia(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arctan(v1,r);
                   end;

        _ARCSINH:  begin
                     mpc_arcsinh(v1,r);
                   end;

        _ARCTANH:  begin
                     if mpc_is1a(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arctanh(v1,r)
                   end;

            _COS:  begin
                     if li>30 then evr.Err := Err_Overflow
                     else mpc_cos(v1,r);
                   end;

           _COSH:  begin
                     if lr>30 then evr.Err := Err_Overflow
                     else mpc_cosh(v1,r);
                   end;

            _EXP:  begin
                     if lr>30 then evr.Err := Err_Overflow
                     else mpc_exp(v1,r);
                   end;

          _EXPM1:  begin
                     if lr>30 then evr.Err := Err_Overflow
                     else mpc_expm1(v1,r);
                   end;

             _LN:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_ln(v1,r)
                   end;

           _LN1P:  begin
                     if (mpf_is1a(v1.re) and s_mpf_is_neg(v1.re)) then evr.Err := Err_Invalid_argument
                     else mpc_ln1p(v1,r)
                   end;

          _LOG10:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_log10(v1,r)
                   end;

            _SIN:  begin
                     if li>30 then evr.Err := Err_Overflow
                     else mpc_sin(v1,r);
                   end;

           _SINH:  begin
                     if lr>30 then evr.Err := Err_Overflow
                     else mpc_sinh(v1,r);
                   end;

            _TAN:  begin
                     mpc_tan(v1,r);
                   end;

           _TANH:  begin
                     mpc_tanh(v1,r);
                   end;

            _COT:  begin
                     {indirect calc to avoid div by 0}
                     mpc_tan(v1,v2);
                     if mpc_is0(v2) then evr.Err := Err_Invalid_argument
                     else mpc_inv(v2,r);
                   end;

            _CSC:  begin
                     {indirect calc to avoid div by 0}
                     mpc_sin(v1,v2);
                     if mpc_is0(v2) then evr.Err := Err_Invalid_argument
                     else mpc_inv(v2,r);
                   end;

            _SEC:  begin
                     {indirect calc to avoid div by 0}
                     mpc_cos(v1,v2);
                     if mpc_is0(v2) then evr.Err := Err_Invalid_argument
                     else mpc_inv(v2,r);
                   end;

           _COTH:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_coth(v1,r);
                   end;

           _CSCH:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_csch(v1,r);
                   end;

           _SECH:  begin
                     mpc_sech(v1,r);
                   end;

         _ARCCOT:  begin
                     if mpc_is_ia(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arccot(v1,r);
                   end;

        _ARCCOTC:  begin
                     if mpc_is_ia(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arccotc(v1,r);
                   end;

         _ARCCSC:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arccsc(v1,r);
                   end;

         _ARCSEC:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arcsec(v1,r)
                   end;

        _ARCCOTH:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arccoth(v1,r);
                   end;

       _ARCCOTHC:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arccothc(v1,r);
                   end;

        _ARCCSCH:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arccsch(v1,r);
                   end;

        _ARCSECH:  begin
                     if mpc_is0(v1) then evr.Err := Err_Invalid_argument
                     else mpc_arcsech(v1,r);
                   end;

           _TEST:  begin
                   end;

            else   evr.Err := Err_Unknown_Operation
      end;
    end;
  end
  else begin
    eval(e^.LNode, v1, evr);
    if evr.Err=0 then begin
      eval(e^.RNode, v2, evr);
    end;
    if evr.Err=0 then begin
      case e^.op of
          _ADD:  mpc_add(v1,v2,r);
          _SUB:  mpc_sub(v1,v2,r);
          _MUL:  begin
                   mpc_mul(v1,v2,r);
                 end;
          _DIV:  begin
                   if mpc_is0(v2) then evr.Err := Err_Division_by_zero
                   else begin
                     mpc_div(v1,v2,r);
                   end;
                 end;

         _EXPT:  begin
                   mpc_pow(v1,v2,r);
                 end;

          _AGM:  begin
                   mpc_agm(v1,v2,r);
                 end;

        _NROOT:  begin
                   evr.Err := Err_Invalid_argument;
                   if mpf_is0(v2.im) then begin
                     td := mpf_todouble(v2.re);
                     if (abs(td)<>DblPosInf) and (frac(td)=0) then begin
                       lx := round(td);
                       evr.Err := 0;
                       mpc_nroot(v1,lx,r);
                     end
                   end;
                 end;

          else   evr.Err := Err_Unknown_Operation
      end;
    end;
  end;
  if (evr.ERR=0) and (MP_Error<>MP_OKAY) then evr.Err := Err_MPERR_Eval;
  mpc_clear2(v1,v2);
end;


{---------------------------------------------------------------------------}
procedure mpc_clear_expr(var e: PFExpr);
  {-Release memory used by e and clear mp_complex values}
begin
  if e<>nil then with e^ do begin
    case nn of
         0:  if op=_CONST then mpc_clear(Value);
         1:  if SNode<>nil then mpc_clear_expr(SNode);
       else  begin
               if e^.LNode<>nil then mpc_clear_expr(e^.LNode);
               if e^.RNode<>nil then mpc_clear_expr(e^.RNode);
             end;
    end;
    mp_freemem(pointer(e),sizeof(TFExpr));
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_eval(e: PFExpr; var evr: TFEval);
  {-Evaluate expression tree e, result in evr}
begin
  with evr do begin
    Err := 0;
    eval(e,Res,evr);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_calculate(psz: pchar8; var evr: TFEval; var EPos: integer);
  {-Parse and evaluate string psz}
var
  e: PFExpr;
  pc: pchar8;
begin
  e := nil;
  pc := mpc_parse(psz,e,evr.Err);
  EPos := pc-psz;
  if evr.Err=0 then mpc_eval(e, evr);
  mpc_clear_expr(e);
end;


{---------------------------------------------------------------------------}
function mpc_calc_errorstr(Err: integer): mp_string;
  {-Translate known error codes}
var
  s: string[20];
begin
  case Err of
    Err_Missing_LeftBracket  : mpc_calc_errorstr := 'Missing "("';
  {$ifdef Unicode}
    Err_Missing_Comma        : mpc_calc_errorstr := mp_string('Missing argument separator ("'+mp_arg_sep+'")');
  {$else}
    Err_Missing_Comma        : mpc_calc_errorstr := 'Missing argument separator ("'+mp_arg_sep+'")';
  {$endif}
    Err_Missing_RightBracket : mpc_calc_errorstr := 'Missing ")"';
    Err_Unknown_Function     : mpc_calc_errorstr := 'Unknown function';
    Err_Unknown_Element      : mpc_calc_errorstr := 'Unknown element';
    Err_Trailing_Garbage     : mpc_calc_errorstr := 'Trailing garbage';
    Err_Invalid_Number       : mpc_calc_errorstr := 'Invalid number';
    Err_Unknown_Operation    : mpc_calc_errorstr := 'Unknown operation';
    Err_MPERR_Eval           : mpc_calc_errorstr := 'MP_Err <> MP_OK';
    Err_Division_by_zero     : mpc_calc_errorstr := 'Division by zero';
    Err_Invalid_argument     : mpc_calc_errorstr := 'Invalid argument(s)';
    Err_X_not_init           : mpc_calc_errorstr := 'Variable X not initialized';
    Err_Y_not_init           : mpc_calc_errorstr := 'Variable Y not initialized';
    Err_Z_not_init           : mpc_calc_errorstr := 'Variable Z not initialized';
    Err_Overflow             : mpc_calc_errorstr := 'Overflow';
    else begin
      str(Err,s);
      mpc_calc_errorstr := 'mpf_calc error '+s;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_init_eval(var evr: TFEval);
  {-Initialize the mp_complex records of evr}
begin
  with evr do mpc_initp4(X,Y,Z,Res,mpf_get_default_prec);
  evr.Err := 0;
end;


{---------------------------------------------------------------------------}
procedure mpc_clear_eval(var evr: TFEval);
  {-Clear the mp_complex records of evr}
begin
  with evr do mpc_clear4(X,Y,Z,Res);
end;


end.


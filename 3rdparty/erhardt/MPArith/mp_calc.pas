unit mp_calc;

{Parse and evaluate mp_int expressions}

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

 DESCRIPTION   :  Parse and evaluate mp_int expressions

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25S, FPC, VP, WDOSX

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
 0.0.10   22.09.06  W.Ehrhardt  First working version
 0.0.11   22.09.06  we          uses mp_getmem, mp_freemem
 0.0.12   23.09.06  we          nn field, different make node routines
 0.0.13   23.09.06  we          only one sign character
 0.0.14   23.09.06  we          GetTL with LMax parameter
 0.0.15   23.09.06  we          Eval routines with Err variable, unit mp_calc
 0.0.16   24.09.06  we          mp_calc_errorstr, eval errors are negative
 0.0.17   25.09.06  we          Functions: abs, kronecker
 0.0.18   30.09.06  we          TEval record, mp_init_eval, mp_clear_eval
 0.0.19   30.09.06  we          Variables X,Y,Z; function sqr
 0.0.20   14.03.07  we          'parse error' changed to 'mp_calc error'
 0.0.21   15.03.07  we          Return Err_MPERR_Eval if MP_Error in eval/eval_mod
 0.0.22   15.03.07  we          Err_Invalid_argument for 0^-n

 1.0.00   07.05.07  we          Handle case v1^(-v2) mod m in eval_mod
 1.0.01   08.05.07  we          Error check for root() and jacobi()
 1.0.02   09.05.07  we          invmod via mp_xgcd avoids error or pre-check

 1.2.00   07.09.07  we          binomial function
 1.2.01   08.09.07  we          Error check: estimate bitsize of binomial(n,k)
 1.2.02   08.09.07  we          <expr> mod 1 is set to zero without evaluation
 1.2.03   08.09.07  we          estimate bitsize(mp_binomial) using lnfac

 1.3.00   14.10.07  we          corrected imprecise lnfac from Num. Recipes (1.ed)
 1.3.01   24.10.07  we          hex numbers in expressions
 1.3.02   24.10.07  we          sqrtmod

 1.5.00   16.01.08  we          clear local PExpr if error
 1.5.01   30.01.08  we          boolean functions and, or, xor

 1.6.00   17.05.08  we          simple invalid check for _EXPT
 1.6.01   05.06.08  we          _ISPPRIME = is probable prime

 1.7.00   13.09.08  we          eval: _INVMOD with mp_invmodf
 1.7.01   17.09.08  we          eval: GetTL with mp_is_longint
 1.7.02   24.09.08  we          string replaced by mp_string

 1.8.00   05.11.08  we          second argument in root must be > 0

 1.9.00   02.12.08  we          Uses BTypes: char8, pchar8

 1.10.00  09.01.09  we          eval/_EXPT with mp_is_longint and mp_expt_int
 1.10.01  21.01.09  we          changes related to (s)mp_divrem
 1.10.02  02.02.09  we          cbrtmod
 1.10.03  03.02.09  we          cbrtmod for prime powers

 1.11.00  23.03.09  we          initialize e:=nil in Element,Expr,Factor,Func,Term
 1.11.01  30.03.09  we          removed redefinition of str255

 1.15.00  13.05.10  we          mp_arg_sep

 1.17.00  06.10.10  we          _DFACT: double factorial
 1.17.01  27.12.10  we          changed parsing of !!
 1.17.02  30.12.10  we          binok for k < 0
 1.17.03  02.01.11  we          avoid D12 warning in mp_calc_errorstr

 1.19.00  19.11.11  we          Check if 2nd argument of sqrt/cbrtmod is prime power

 1.21.00  15.07.12  we          _PRIMEPI: return primepi32
 1.21.01  19.07.12  we          _CATALAN
 1.21.02  20.07.12  we          _POPCOUNT, _RANDPRIME, _MAURER, _PERM, _POCH
 1.21.03  27.07.12  we          _EULERPHI, _CARMICHAEL, _MOEBIUS, _PRIMROOT, _ORDER
 1.21.04  30.07.12  we          Special case 0^b, b>=0
 1.21.05  30.07.12  we          _SIGMA

 1.22.00  07.08.12  we          _PRIME
 1.22.01  09.09.12  we          _QNR

 1.24.00  30.10.12  we          _VAL
 1.24.01  03.01.13  we          Fix typo: Err_Trailing_Gargabe -> Err_Trailing_Garbage

 1.25.00  02.02.13  we          Fix _SQRTMOD and _CBRTMOD if v2=prime^1

 1.29.00  16.07.14  we          _SAFEPRIME

 1.30.00  27.09.14  we          _RANDOM uses mp_random
 1.30.01  02.10.14  we          Keep some compilers happy and initialize tl := 0

 1.31.00  05.11.14  we          dfact for n=-1, -3

 1.33.00  11.08.16  we          fix error for JACOBI
 1.33.01  08.09.16  we          _ISCARMICHAEL, function names now string[12]

 1.34.00  09.12.16  we          _ISCARMICHAEL uses mp_is_pcarmichael

 1.35.00  26.07.17  we          _SWING

 1.36.00  11.09.17  we          _PSP
 1.36.01  12.09.17  we          _NEXTPOWER
 1.36.02  12.09.17  we          _DIGITSUM

 1.37.00  14.10.17  we          input in radix from (2..MaxRadix)
 1.37.01  17.10.17  we          Fix value for Err_Invalid_RadixNumber
 1.37.02  16.05.18  we          Avoid implicit string conversion warning

**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2006-2018 Wolfgang Ehrhardt

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
  Err_Invalid_RadixNumber  = 10;   {Invalid radix number}

  Err_MPERR_Eval           = -1;   {MP_Err <> MP_OK}
  Err_Division_by_zero     = -2;   {Division by zero}
  Err_Invalid_argument     = -3;   {Invalid arg, e.g. sqrt(-2)}
  Err_X_not_init           = -4;   {Variable X not initialized}
  Err_Y_not_init           = -5;   {Variable Y not initialized}
  Err_Z_not_init           = -6;   {Variable Z not initialized}
  Err_no_solution          = -7;   {(sqrtmod:) no solution}
{#Z-}


type
  TOperation = (_CONST, _CHS, _ABS, _ADD, _SUB, _MUL, _DIV, _MOD, _EXPT,
                _SQRT, _SQR, _ROOT, _GCD, _LCM, _INVMOD, _RANDOM, _MIN, _MAX,
                _FACT, _DFACT, _PRIMORIAL, _NEXTPRIME, _PREVPRIME, _JACOBI, _KRONECKER,
                _FIB, _LUC, _FERMAT, _MERSENNE, _X, _Y, _Z, _BINOMIAL,
                _SQRTMOD, _AND, _OR, _XOR, _ISPPRIME, _CBRTMOD, _PRIMEPI,
                _CATALAN, _POPCOUNT, _RANDPRIME, _MAURER, _PERM, _POCH,
                _EULERPHI, _CARMICHAEL, _MOEBIUS, _PRIMROOT, _ORDER, _SPSP,
                _SIGMA, _PRIME, _QNR, _VAL, _SAFEPRIME, _ISCARMIKE,
                _SWING, _PSP, _NEXTPOWER, _DIGITSUM, _TEST);
                {implemented operators, functions, and variables}

type
  PExpr  = ^TExpr;                       {Expression node pointer}
  TExpr  = record                        {binary tree node}
             op:  TOperation;            {operation/function/variable}
             nn:  byte;                  {number of nodes}
             case integer of
               0: (Value: mp_int);       {value if op=_CONST}
               1: (SNode: PExpr);        {expr = (SNode op) or op(SNode)}
               2: (LNode, RNode: PExpr;) {expr = LNode op RNode or op(LNode,RNode)}
           end;

  TEval  = record                        {Evaluation record}
             X  : mp_int;                {Variable X}
             Y  : mp_int;                {Variable Y}
             Z  : mp_int;                {Variable Z}
             Res: mp_int;                {Evaluation result}
             Err: integer;               {Eval error code}
           end;


function  mp_parse(psz: pchar8; var e: PExpr; var Err: integer): pchar8;
  {-Parse string psz into expression tree e, if OK Err=0 and result^=#0,}
  { else Err=Err_xx and result points to error position}

procedure mp_eval(e: PExpr; var evr: TEval);
  {-Evaluate expression tree e, result in evr}

procedure mp_eval_mod(e: PExpr; const m: mp_int; var evr: TEval);
  {-Evaluate expression tree e mod m, result in evr}

procedure mp_clear_expr(var e: PExpr);
  {-Release memory used by e}

procedure mp_calculate(psz: pchar8; var evr: TEval; var EPos: integer);
  {-Parse and evaluate string psz}

function  mp_calc_errorstr(Err: integer): mp_string;
  {-Translate known error codes}

procedure mp_init_eval(var evr: TEval);
  {-Initialize the mp_ints of evr}

procedure mp_clear_eval(var evr: TEval);
  {-Clear the mp_ints of evr}


implementation


uses
  mp_base, mp_prime, mp_modul, mp_numth, mp_pfu;


(* Approx. grammar for expression parser:

  <Expr>    ::=   <Term> '+' <Term>
                | <Term> '-' <Term>
                | '+'<Term>
                | '-'<Term>;
  <Term>    ::=   <Factor> '*' <Factor>
                | <Factor> '/' <Factor>
                | <Factor> '%' <Factor>
                | <Factor> 'div' <Factor>
                | <Factor> 'mod' <Factor>
                | <Factor>;
  <Factor>  ::=   <Element> '^' <Factor>
                | <Element> {'!' | '!!' | '#'}
                | <Element>;
  <Element> ::= <Func> | <Var> | <number> | <hexnum> | radixnum | '(' <Expr> ')';
  <Func>    ::= <Ident> '(' <Arglist> ')';
  <Var>     ::= 'X'..'Z' | 'x'..'z';
  <Arglist  ::= <Expr> | <Expr> ',' <Expr>;
  <Ident>   ::= <alpha> {<alpha>};
  <number>  ::= ['+' | '-'] <digit> { <digit> };
  <digit>   ::= '0'..'9';
  <hexnum>  ::= '$' <hexdigit> { <hexdigit> };
  <radixnum>::= '[' <number> ']' <radixdigit> { <radixdigit> };
  <hexdigit>::= '0'..'9'| 'A'..'F'| 'a'..'f';
  <alpha>   ::= 'A'..'Z'| 'a'..'z';
*)



type
  str12  = string[12];
  TFunc  = record
               op: TOperation;
             arg2: boolean;
             name: str12;
           end;
const
  MaxFun = 50;

const
  FuncTab : array[1..MaxFun] of TFunc = (
             (op: _SQRT      ; arg2:false ;name: 'SQRT'),
             (op: _SQR       ; arg2:false ;name: 'SQR'),
             (op: _ABS       ; arg2:false ;name: 'ABS'),
             (op: _ROOT      ; arg2:true  ;name: 'ROOT'),           {root(a,b)}
             (op: _GCD       ; arg2:true  ;name: 'GCD'),
             (op: _LCM       ; arg2:true  ;name: 'LCM'),
             (op: _INVMOD    ; arg2:true  ;name: 'INVMOD'),
             (op: _NEXTPRIME ; arg2:false ;name: 'NEXTPRIME'),
             (op: _PREVPRIME ; arg2:false ;name: 'PREVPRIME'),
             (op: _FIB       ; arg2:false ;name: 'FIB'),
             (op: _LUC       ; arg2:false ;name: 'LUC'),
             (op: _FERMAT    ; arg2:false ;name: 'FERMAT'),
             (op: _MERSENNE  ; arg2:false ;name: 'MERSENNE'),
             (op: _RANDOM    ; arg2:false ;name: 'RANDOM'),
             (op: _MIN       ; arg2:true  ;name: 'MIN'),
             (op: _MAX       ; arg2:true  ;name: 'MAX'),
             (op: _JACOBI    ; arg2:true  ;name: 'JACOBI'),
             (op: _KRONECKER ; arg2:true  ;name: 'KRONECKER'),
             (op: _BINOMIAL  ; arg2:true  ;name: 'BINOMIAL'),
             (op: _SQRTMOD   ; arg2:true  ;name: 'SQRTMOD'),
             (op: _AND       ; arg2:true  ;name: 'AND'),
             (op: _OR        ; arg2:true  ;name: 'OR'),
             (op: _XOR       ; arg2:true  ;name: 'XOR'),
             (op: _CBRTMOD   ; arg2:true  ;name: 'CBRTMOD'),
             (op: _ISPPRIME  ; arg2:false ;name: 'ISPPRIME'),
             (op: _DFACT     ; arg2:false ;name: 'DFACT'),
             (op: _PRIMEPI   ; arg2:false ;name: 'PRIMEPI'),
             (op: _POPCOUNT  ; arg2:false ;name: 'POPCOUNT'),
             (op: _MAURER    ; arg2:false ;name: 'MAURER'),
             (op: _RANDPRIME ; arg2:false ;name: 'RANDPRIME'),
             (op: _CATALAN   ; arg2:false ;name: 'CATALAN'),
             (op: _PERM      ; arg2:true  ;name: 'PERM'),
             (op: _POCH      ; arg2:true  ;name: 'POCH'),
             (op: _EULERPHI  ; arg2:false ;name: 'EULERPHI'),
             (op: _CARMICHAEL; arg2:false ;name: 'CARMICHAEL'),
             (op: _MOEBIUS   ; arg2:false ;name: 'MOEBIUS'),
             (op: _PRIMROOT  ; arg2:false ;name: 'PRIMROOT'),
             (op: _ORDER     ; arg2:true  ;name: 'ORDER'),
             (op: _SPSP      ; arg2:true  ;name: 'SPSP'),
             (op: _SIGMA     ; arg2:true  ;name: 'SIGMA'),
             (op: _PRIME     ; arg2:false ;name: 'PRIME'),
             (op: _QNR       ; arg2:false ;name: 'QNR'),
             (op: _VAL       ; arg2:true  ;name: 'VAL'),
             (op: _SAFEPRIME ; arg2:false ;name: 'SAFEPRIME'),
             (op: _ISCARMIKE ; arg2:false ;name: 'ISCARMICHAEL'),
             (op: _SWING     ; arg2:false ;name: 'SWING'),
             (op: _PSP       ; arg2:true  ;name: 'PSP'),
             (op: _NEXTPOWER ; arg2:false ;name: 'NEXTPOWER'),
             (op: _DIGITSUM  ; arg2:true  ;name: 'DIGITSUM'),
             (op: _TEST      ; arg2:true  ;name: 'TEST'));


type
  TOpVar = _X.._Z;

{---------------------------------------------------------------------------}
function SkipWhite(psz: pchar8): pchar8;
  {-Skip white space}
begin
  while psz^ in [' ',#13,#10,#9] do inc(psz);
  SkipWhite := psz;
end;


{---------------------------------------------------------------------------}
procedure mkNode(var r: PExpr; op: TOperation; nn: byte; e1, e2: PExpr);
  {-Make a new expression node for (e1 op e2) or (e1 op)}
begin
  r := mp_alloc(sizeof(TExpr));
  r^.nn   := nn;
  r^.op   := op;
  r^.LNode := e1;
  r^.RNode := e2;
end;


{---------------------------------------------------------------------------}
procedure mkNode0c(var r: PExpr);
  {-Make a new node for a constant}
begin
  {alloc Value node initialize mp_int}
  r := mp_alloc(sizeof(TExpr));
  r^.op := _CONST;
  r^.nn := 0;
  mp_init(r^.Value);
end;


{---------------------------------------------------------------------------}
procedure mkNode0v(var r: PExpr; const v: TOpVar);
  {-Make a new node for a variable v = _X, _Y, or _Z}
begin
  {alloc Value node initialize mp_int}
  r := mp_alloc(sizeof(TExpr));
  r^.op := v;
  r^.nn := 0;
end;


{---------------------------------------------------------------------------}
procedure mkNode1(var r: PExpr; op: TOperation; e: PExpr);
  {-Make a new expression node for r := e op}
begin
  mkNode(r,op,1,e,nil);
end;


{---------------------------------------------------------------------------}
procedure mkNode2(var r: PExpr; op: TOperation; e1, e2: PExpr);
  {-Make a new expression node for r := e1 op e2}
begin
  mkNode(r,op,2,e1,e2);
end;


{---------------------------------------------------------------------------}
function GetIdent(psz: pchar8): mp_string;
  {-Gather next identifier, break if not in [A-Z, a-z], result is uppercase}
var
  s: mp_string;
begin
  s := '';
  while psz^ in ['A'..'Z', 'a'..'z'] do begin
    s := s+upcase(psz^);
    inc(psz);
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
function  Expr(psz: pchar8; var e: PExpr; var Err: integer): pchar8; forward;
  {-Parse string psz into expression tree}
procedure eval_mod(e: PExpr; const m: mp_int; var r: mp_int; var evr: TEval); forward;
  {-Evaluate r = e mod m}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function Func(psz: pchar8; idx: integer; var e: PExpr; var Err: integer): pchar8;
  {-Build function expression, psz points between function name and "("}
var
  e1,e2: PExpr;
const
  na: array[boolean] of byte = (1,2);

  procedure clear;
    {-Clear local PExpr if error}
  begin
    if e1<>nil then mp_clear_expr(e1);
    if e2<>nil then mp_clear_expr(e2);
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
    clear;    {V1.5.00}
    exit;
  end;
  psz := SkipWhite(psz);
  if FuncTab[idx].arg2 then begin
    {if two arguments search for ","}
    if psz^ <> mp_arg_sep{','} then begin
      Func := psz;
      Err := Err_Missing_Comma;
      clear; {V1.5.00}
      exit;
    end;
    {evaluate second argument}
    psz := Expr(psz+1,e2,Err);
    Func := psz;
    if Err<>0 then begin
      clear;    {V1.5.00}
      exit;
    end;
    psz := SkipWhite(psz);
  end
  else e2:=nil;
  {search for closing ")"}
  if psz^ <> ')' then begin
    Func := psz;
    Err := Err_Missing_RightBracket;
    clear; {V1.5.00}
    exit;
  end;
  inc(psz);
  with FuncTab[idx] do mkNode(e, op, na[arg2], e1, e2);
  Func := psz;
end;


{---------------------------------------------------------------------------}
function Element(psz: pchar8; var e: PExpr; var Err: integer): pchar8;
  {-Parse an Element}
var
  res: PExpr;
  s,sc: pchar8;
  lsc: word;
  i: integer;
  r: longint;
  s0: char8;
  neg: boolean;
  id: str255;
begin
  e := nil;
  Element := psz;
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
  else if s0 in ['0'..'9','+','-'] then begin
    {get sign}
    neg := false;
    if psz^ in ['+','-'] then begin
      if psz^='-' then neg := true;
      inc(psz);
    end;
    {count decimal digits}
    s := psz;
    while psz^ in ['0'..'9'] do inc(psz);
    {alloc and move digit string to temp storage}
    if s=psz then begin
      {empty digit string}
      Err := Err_Invalid_Number;
      Element := psz;
      exit;
    end;
    lsc := psz-s+1;
    sc := mp_getmem(lsc);
    move(s^,sc^,pred(lsc));
    sc[psz-s] := #0;
    {Make a constant node}
    mkNode0c(res);
    {convert digit string to mp_int}
    mp_read_decimal(res^.Value,sc);
    {apply sign}
    if neg then mp_chs(res^.Value, res^.Value);
    e := res;
    mp_freemem(pointer(sc),lsc);
  end
  else if s0='$' then begin
    sc := psz;
    inc(psz);
    {count hex digits}
    s := psz;
    while upcase(psz^) in ['0'..'9','A'..'F'] do inc(psz);
    {alloc and move hex string to temp storage}
    if s=psz then begin
      {empty hex string}
      Err := Err_Invalid_HexNumber;
      Element := sc;
      exit;
    end;
    lsc := psz-s+1;
    sc := mp_getmem(lsc);
    move(s^,sc^,pred(lsc));
    sc[psz-s] := #0;
    {Make a constant node}
    mkNode0c(res);
    {convert hex string to mp_int}
    mp_read_radix(res^.Value,sc,16);
    e := res;
    mp_freemem(pointer(sc),lsc);
  end
  else if s0='[' then begin
    {format [rr]xxxxxxxxxxxx, rr=00..36}
    sc := psz;
    inc(psz);
    i := 0;
    id := '';
    repeat
      s0 := psz^;
      inc(psz);
      if (s0=']') or (s0=#0) then break;
      inc(i);
      id := id + s0;
    until i=255;
    {$ifdef D12Plus}
      val(string(id),r,i);
    {$else}
      val(id,r,i);
    {$endif}
    if (s0<>']') or (i<>0) or (r<2) or (r>MAXRadix) then begin
      {empty radix number string}
      Err :=  Err_Invalid_RadixNumber;
      Element := sc;
      exit;
    end;
    {count radix digits}
    s := psz;
    while pos(upcase(psz^), mp_ucrmap)>0 do inc(psz);
    {alloc and move radix string to temp storage}
    if s=psz then begin
      {empty hex string}
      Err := Err_Invalid_RadixNumber;
      Element := sc;
      exit;
    end;
    lsc := psz-s+1;
    sc := mp_getmem(lsc);
    move(s^,sc^,pred(lsc));
    sc[psz-s] := #0;
    {Make a constant node}
    mkNode0c(res);
    {convert hex string to mp_int}
    mp_read_radix(res^.Value,sc,r);
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
function Factor(psz: pchar8; var e: PExpr; var Err: integer): pchar8;
  {-Parse a Factor}
var
  t: PExpr;
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
  {Look for optional factorial, primorial, power part of Factor. }
  {Note: Starting with V1.17.00 n!! is now parsed as (n!!) not as}
  {      (n!)!, n!!! as (n!!)!, n!!!! as (n!!)!! etc             }
  while psz^ in ['!', '#'] do begin
    if psz^='#' then begin
      mkNode1(e, _PRIMORIAL, e);
      inc(psz);
    end
    else begin
      inc(psz);
      if psz^='!' then begin
        mkNode1(e, _DFACT, e);
        inc(psz);
      end
      else mkNode1(e, _FACT, e);
    end;
  end;
  psz := SkipWhite(psz);
  if psz^='^' then begin
    inc(psz);
    t := nil;
    psz := Factor(psz, t, Err);
    if Err=0 then mkNode2(e, _EXPT, e, t)
    else if t<>nil then mp_clear_expr(t); {V1.5.00}
  end;
  Factor := psz;
end;


{---------------------------------------------------------------------------}
function Term(psz: pchar8; var e: PExpr; var Err: integer): pchar8;
  {-Parse a Term}
var
  t: PExpr;
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
      '/': begin
              psz := Factor(psz+1, t, Err);
              if Err=0 then mkNode2(e, _DIV, e, t);
            end;
      'D' : if GetIdent(psz)='DIV' then begin
              psz := Factor(psz+3, t, Err);
              if Err=0 then mkNode2(e, _DIV, e, t);
            end
            else break;
      '%' : begin
              psz := Factor(psz+1, t, Err);
              if Err=0 then mkNode2(e, _MOD, e, t);
            end;
      'M' : if GetIdent(psz)='MOD' then begin
              psz := Factor(psz+3, t, Err);
              if Err=0 then mkNode2(e, _MOD, e, t);
            end
            else break;

     else   break;
    end;
  end;
  if (Err<>0) and (t<>nil) then mp_clear_expr(t); {V1.5.00}
  Term := psz;
end;


{---------------------------------------------------------------------------}
function Expr(psz: pchar8; var e: PExpr; var Err: integer): pchar8;
  {-Parse an Expr}
var
  t: PExpr;
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
  if (Err<>0) and (t<>nil) then mp_clear_expr(t); {V1.5.00}
  Expr := psz;
end;


{---------------------------------------------------------------------------}
function mp_parse(psz: pchar8; var e: PExpr; var Err: integer): pchar8;
  {-Parse string psz into expression tree e, Err=0 and psz^=#0 if OK,}
  { else Err=Err_xx and result points to error position}
begin
  Err := 0;
  psz := Expr(psz, e, Err);
  if (Err=0) and (psz^<>#0) then Err := Err_Trailing_Garbage;
  mp_parse := psz;
end;


{---------------------------------------------------------------------------}
function lnfac(n: longint): double;
  {-Approx. value for ln(n!) using Lanczos complex Gamma approximation}
const
  {Paul Godfrey, http://my.fit.edu/~gabdo/gammacoeff.txt, g=5}
  coef: array[0..6] of double = ( 1.0000000001900148, 76.180091729471463,
                                 -86.505320329416768, 24.014098240830910,
                                 -1.2317395724501554, 0.1208650973866179E-2,
                                 -0.5395239384953128E-5);
var
  t,s: double;
  i: integer;
begin
  t := n+5.5;
  t := (n+0.5)*ln(t)-t+0.918938533204672742; {ln(sqrt(2*pi))}
  s := coef[0];
  for i:=1 to 6 do s := s+coef[i]/(n+i);
  lnfac := t+ln(s);
end;


{---------------------------------------------------------------------------}
function binok(n,k: longint): boolean;
  {-Check if bitsize(binomial(n,k)) is OK}
var
  t: longint;
begin
  binok := true;
  if k<0 then begin
    if (n>=0) or (n<k) then exit;
    t := n-k;
    n := -succ(k);
    k := t;
  end;
  if n<0 then n := k-n-1;
  if (k>=n) or (k<=0) or (n<2) then exit;
  binok := MaxMersenne > (lnfac(n) - lnfac(k) - lnfac(n-k))/ln(2.0);
end;


{---------------------------------------------------------------------------}
function permok(n,k: longint): boolean;
  {-Check if bitsize(perm(n,k)) is OK}
begin
  if (n<0) or (k<0) then permok := false
  else if k>n then permok := true
  else permok := MaxMersenne > (lnfac(n) - lnfac(n-k))/ln(2.0);
end;


{---------------------------------------------------------------------------}
function sigmaok(k,n: longint): boolean;
var
  a: double;
{$ifdef MPC_NOHALT}
const
  am = MaxMersenne*0.693;
{$else}
const
  am = MaxMersenne*0.69;  {greater safety margin if error will halt program}
{$endif}
begin
  sigmaok := true;
  n := abs(n);
  if (n<=1) or (k<=0) then exit;
  a := k*ln(n);
  sigmaok := a < am;
end;


{---------------------------------------------------------------------------}
function pochok(n,k: longint): boolean;
  {-Check if bitsize(poch(n,k)) is OK}
begin
  pochok := true;
  if k<2 then exit;
  if (n <= 0) and (n+k >= 1) then exit;
  if n<0 then n := 1-n-k; {use absolute value}
  {poch(n,k) = gamma(n+k)/gamma(n) = (n+k-1)!/(n-1)! }
  pochok := MaxMersenne > (lnfac(n+k-1) - lnfac(n-1))/ln(2.0);
end;


{---------------------------------------------------------------------------}
procedure eval(e: PExpr; var r: mp_int; var evr: TEval);
  {-(internal) evaluate expression tree e, result in r}
var
  v1,v2: mp_int;
  tl,tl2: longint;

  {---------------------------------------------------}
  function GetTL(const a: mp_int; LMin,LMax: longint): boolean;
    {-Return true and value in tl if a has at most 31 bit and is in [LMin, LMax]}
  var
    OK: boolean;
  begin
    OK := mp_is_longint(a,tl) and (tl>=LMin) and (tl<=LMax);
    if (not OK) and (evr.Err=0) then evr.Err := Err_Invalid_argument;
    GetTL := OK;
  end;

begin
  if evr.Err<>0 then exit;

  if e^.nn=0 then begin
    case e^.op of
      _CONST: begin
                mp_copy(e^.Value,r);
                if MP_Error<>0 then evr.Err := Err_MPERR_Eval;
              end;
          _X: begin
                if mp_not_init(evr.X) then evr.Err := Err_X_not_init
                else mp_copy(evr.X,r);
              end;
          _Y: begin
                if mp_not_init(evr.Y) then evr.Err := Err_Y_not_init
                else mp_copy(evr.Y,r);
              end;
          _Z: begin
                if mp_not_init(evr.Z) then evr.Err := Err_Z_not_init
                else mp_copy(evr.Z,r);
              end;
         else evr.Err := Err_Unknown_Operation;
    end;
    exit;
  end;

  {always initialize 2 mp_int although some ops need only 0 or 1}
  mp_init2(v1,v2);
  if MP_Error<>MP_OKAY then begin
    evr.Err := Err_MPERR_Eval;
    exit;
  end;

  tl := 0; {keep some compilers happy}
  if e^.nn=1 then begin
    eval(e^.LNode, v1, evr);
    if evr.Err=0 then begin
      case e^.op of
            _CHS:  mp_chs(v1,r);
            _ABS:  mp_abs(v1,r);
           _FACT:  if GetTL(v1, 0, MaxFact)      then mp_fact(tl,r);
          _DFACT:  if GetTL(v1, -3, 15*MaxFact div 8) and (tl<>-2) then mp_dfact(tl,r);
      _PRIMORIAL:  if GetTL(v1, 0, MaxPrimorial) then mp_primorial(tl, r);
          _SWING:  if GetTL(v1, 0, MaxLongint)   then mp_swing(tl,r);
        _PRIMEPI:  if GetTL(v1, 0, MaxLongint)   then mp_set_int(r, primepi32(tl));
       _EULERPHI:  if GetTL(v1, 0, MaxLongint)   then mp_set_int(r, EulerPhi32(tl));
     _CARMICHAEL:  if GetTL(v1, 0, MaxLongint)   then mp_set_int(r, Carmichael32(tl));
   {  _ISCARMIKE:  if GetTL(v1, 0, MaxLongint)   then mp_set(r,ord(is_Carmichael32(tl)) and 1);}
      _ISCARMIKE:  mp_set(r,ord(mp_is_pcarmichael(v1)) and 1);
        _MOEBIUS:  if GetTL(v1, 0, MaxLongint)   then mp_set_int(r, Moebius32(tl));
       _PRIMROOT:  if GetTL(v1, 0, MaxLongint)   then mp_set_int(r, primroot32(tl));
          _PRIME:  if GetTL(v1, 1, 105097565)    then mp_set_int(r, prime32(tl));

         _FERMAT:  if GetTL(v1, 0, MaxFermat)    then mp_fermat(tl,r);
        _CATALAN:  if GetTL(v1, 0, MaxCatalan)   then mp_catalan(tl,r);
       _MERSENNE:  if GetTL(v1, 0, MaxMersenne)  then mp_mersenne(tl,r);
            _FIB:  if GetTL(v1, -MaxFibonacci, MaxFibonacci) then mp_fib(tl,r);
            _LUC:  if GetTL(v1, -MaxLucas, MaxLucas) then mp_lucas(tl,r);
            _SQR:  begin
                     if 2*mp_bitsize(v1) < MP_MAXBIT then mp_sqr(v1,r)
                     else evr.Err := Err_Invalid_argument;
                   end;
           _SQRT:  begin
                     if v1.sign=MP_ZPOS then mp_sqrt(v1,r)
                     else evr.Err := Err_Invalid_argument;
                   end;
      _NEXTPRIME:  begin
                     mp_exch(v1,r);
                     mp_nextprime(r);
                   end;
                   {
           _TEST:  begin
                     mp_set(r,ord(mp_is_slpsp(v1)) and 1);
                   end;
                   }

      _NEXTPOWER:  begin
                     mp_exch(v1,r);
                     mp_nextpower(r);
                   end;
      _PREVPRIME:  begin
                     mp_exch(v1,r);
                     mp_prevprime(r);
                   end;
      _SAFEPRIME:  begin
                     mp_exch(v1,r);
                     mp_safeprime(r);
                   end;
       _ISPPRIME:  begin
                     mp_set(r,ord(mp_is_pprime(v1)) and 1);
                   end;
       _POPCOUNT:  begin
                     mp_set_int(r,mp_popcount(v1));
                   end;
            _QNR:  begin
                     if v1.sign=MP_ZPOS then mp_set_int(r,mp_qnr(v1))
                     else evr.Err := Err_Invalid_argument;
                   end;
         _MAURER:  if GetTL(v1, 2, MP_MAXBIT-1) then begin
                     mp_provable_prime(tl,r);
                   end;
      _RANDPRIME:  if GetTL(v1, 2, MP_MAXBIT-1) then begin
                     mp_rand_prime(tl,pt_normal,r);
                   end;
         _RANDOM:  begin
                     {
                     mp_abs(v1,v1);
                     tl := mp_bitsize(v1);
                     mp_rand_bits(r, tl);
                     if (tl>0) and mp_is_ge(r,v1) then mp_mod(r,v1,r);
                     }
                     mp_random(v1,r);
                   end;
            else   evr.Err := Err_Unknown_Operation
      end;
    end;
  end
  else begin
    if e^.op = _MOD then begin
      {special optimization for mod: First evaluate right expression}
      {then evaluate left expression using modular operations}
      eval(e^.RNode, v2, evr);
      if evr.Err=0 then begin
        {if right expression = 1, always result is always 0}
        if mp_is1(v2) then mp_zero(r)
        else eval_mod(e^.LNode, v2, r, evr);
      end;
    end
    else begin
      eval(e^.LNode, v1, evr);
      if evr.Err=0 then begin
        eval(e^.RNode, v2, evr);
        if (evr.Err=0) and mp_is0(v2) and (e^.op in [_DIV, _MOD, _INVMOD]) then evr.Err := Err_Division_by_zero;
      end;
      if evr.Err=0 then begin
        case e^.op of
            _ADD:  mp_add(v1,v2,r);
            _SUB:  mp_sub(v1,v2,r);
            _MUL:  begin
                     tl  := mp_bitsize(v1);
                     tl2 := mp_bitsize(v1);
                     if tl+tl2 < MP_MAXBIT then mp_mul(v1,v2,r)
                     else evr.Err := Err_Invalid_argument;
                   end;
            _DIV:  mp_div(v1,v2,r);
           _ROOT:  if GetTL(v2,1,MaxLongint) then begin
                     if (not odd(tl)) and (v1.sign=MP_NEG) then evr.Err := Err_Invalid_argument
                     else mp_n_root(v1,tl,r);
                   end;
            _GCD:  mp_gcd(v1,v2,r);
            _LCM:  mp_lcm(v1,v2,r);
         _INVMOD:  if not mp_invmodf(v1,v2,r) then evr.Err := Err_no_solution;
        _SQRTMOD:  begin
                     evr.Err := Err_Invalid_argument;
                     mp_is_power_max(v2,r,tl);
                     if (tl>=1) and mp_is_pprime(r) then begin
                       mp_sqrtmodpk(v1,r,tl,r,evr.Err);
                       if evr.Err=-1 then evr.Err := Err_no_solution
                       else if evr.Err<>0 then evr.Err := Err_Invalid_argument;
                     end;
                   end;
        _CBRTMOD:  begin
                     evr.Err := Err_Invalid_argument;
                     mp_is_power_max(v2,r,tl);
                     if (tl>=1) and mp_is_pprime(r) then begin
                       mp_cbrtmodpk(v1,r,tl,r,evr.Err);
                       if evr.Err=1 then evr.Err := Err_Invalid_argument
                       else if evr.Err<>0 then evr.Err := Err_no_solution;
                     end;
                   end;
            _MIN:  if mp_is_lt(v1,v2) then mp_exch(v1,r) else mp_exch(v2,r);
            _MAX:  if mp_is_gt(v1,v2) then mp_exch(v1,r) else mp_exch(v2,r);
         _JACOBI:  if (mp_cmp_d(v2, 3)=MP_LT) or mp_iseven(v2) then evr.Err := Err_Invalid_argument
                   else mp_set_int(r,mp_jacobi(v1,v2));
      _KRONECKER:  mp_set_int(r,mp_kronecker(v1,v2));
       _BINOMIAL:  if GetTL(v2,-MaxLongint,MaxLongint) then begin
                     tl2 := tl;
                     if GetTL(v1,-MaxLongint,MaxLongint) then begin
                       if binok(tl,tl2) then mp_binomial(tl,tl2,r)
                       else evr.Err := Err_Invalid_argument;
                     end;
                   end;
            _VAL:  if GetTL(v2, 2, MaxLongint) then begin
                     mp_set_int(r,mp_val(v1,tl));
                   end
                   else evr.Err := Err_Invalid_argument;

           _PERM:  if GetTL(v2,0,MaxLongint) then begin
                     tl2 := tl;
                     if GetTL(v1,0,MaxLongint) then begin
                       if permok(tl,tl2) then mp_perm(tl,tl2,r)
                       else evr.Err := Err_Invalid_argument;
                     end;
                   end;
           _POCH:  if GetTL(v2,-MaxLongint,MaxLongint) then begin
                     tl2 := tl;
                     if GetTL(v1,-MaxLongint,MaxLongint) then begin
                       if pochok(tl,tl2) then mp_poch(tl,tl2,r)
                       else evr.Err := Err_Invalid_argument;
                     end;
                   end;
          _SIGMA:  if GetTL(v2,-MaxLongint,MaxLongint) then begin
                     tl2 := tl;
                     if GetTL(v1,0,MaxLongint) then begin
                       if sigmaok(tl,tl2) then mp_sigmak(tl,tl2,r)
                       else evr.Err := Err_Invalid_argument;
                     end;
                   end;
           _ORDER:  if GetTL(v2,2,MaxLongint) then begin
                     tl2 := tl;
                     if GetTL(v1,-MaxLongint,MaxLongint) then begin
                       mp_set_int(r, order32(tl,tl2));
                     end;
                   end;
           _EXPT:  if v2.sign=MP_NEG then begin
                     case mp_cmp_mag_d(v1,1) of
                       MP_GT : mp_zero(r);
                       MP_EQ : if mp_isodd(v2) then mp_copy(v1,r) else mp_abs(v1,r);
                       MP_LT : evr.Err := Err_Invalid_argument;
                     end;
                   end
                   else begin
                     if mp_is1(v1) then mp_set1(r)
                     else if mp_is0(v1) then begin
                       if mp_is0(v2) then mp_set1(r) else mp_zero(r);
                     end
                     else if mp_is_longint(v2,tl) then begin
                       {v1>=2, therefore v2 must be < 2^31. The next check}
                       {does not catch all cases, e.g. 2^(MaxMersenne-1)}
                       if tl<1+trunc(MaxMersenne/(mp_bitsize(v1)-1.0)) then mp_expt_int(v1,tl,r)
                       else evr.Err := Err_Invalid_argument;
                     end
                     else evr.Err := Err_Invalid_argument;
                   end;
            _AND:  mp_and(v1,v2,r);
             _OR:  mp_or(v1,v2,r);
            _XOR:  mp_xor(v1,v2,r);

           _SPSP:  begin
                     if mp_cmp_d(v2, 1)=MP_GT then begin
                       mp_set(r,ord(mp_is_spsp(v1,v2)) and 1);
                     end
                     else evr.Err := Err_Invalid_argument;
                   end;

            _PSP:  begin
                     if mp_cmp_d(v2, 1)=MP_GT then begin
                       mp_set(r,ord(mp_is_psp(v1,v2)) and 1);
                     end
                     else evr.Err := Err_Invalid_argument;
                   end;
       _DIGITSUM:  if GetTL(v2,2,MaxRadix) then begin
                     mp_set_int(r, mp_digitsum(v1,tl));
                   end;
           _TEST:  begin
                     if GetTL(v2, 2, MaxRadix) then mp_reverse(v1, tl, r);
                   end;
            else   evr.Err := Err_Unknown_Operation
        end;
      end;
    end;
  end;
  if (evr.ERR=0) and (MP_Error<>MP_OKAY) then evr.Err := Err_MPERR_Eval;
  mp_clear2(v1,v2);
end;


{---------------------------------------------------------------------------}
procedure eval_mod(e: PExpr; const m: mp_int; var r: mp_int; var evr: TEval);
  {-(internal) evaluate expression tree e mod m, result in r}
var
  v1, v2: mp_int;
begin
  if evr.Err<>0 then exit;
  if mp_is0(m) then begin
    evr.Err := Err_Division_by_zero;
    exit;
  end;

  mp_init2(v1,v2);
  if MP_Error<>MP_OKAY then begin
    evr.Err := Err_MPERR_Eval;
    exit;
  end;

  case e^.op of
      _EXPT: begin
              {Note: a^-b mod 1 gives error, evaluated as (a mod 1)^-b mod 1 = 0^-b mod 1}
              eval_mod(e^.LNode, m, v1, evr);
              if evr.Err=0 then eval(e^.RNode, v2, evr);
              if evr.Err=0 then begin
                if v2.sign=MP_NEG then begin
                  {v2<0, we need v1^-1 mod m, error if gcd(v1,m)<>1}
                  if mp_gcd1(v1,m,r) and (not mp_is0(v1)) then mp_exptmod(v1,v2,m,r)
                  else evr.Err := Err_Invalid_argument
                end
                else mp_exptmod(v1,v2,m,r);
              end;
            end;
      _ADD: begin
              eval_mod(e^.LNode, m, v1, evr);
              if evr.Err=0 then eval_mod(e^.RNode, m, v2, evr);
              if evr.Err=0 then mp_addmod(v1,v2,m,r);
            end;
      _SUB: begin
              eval_mod(e^.LNode, m, v1, evr);
              if evr.Err=0 then eval_mod(e^.RNode, m, v2, evr);
              if evr.Err=0 then mp_submod(v1,v2,m,r);
            end;
      _MUL: begin
              eval_mod(e^.LNode, m, v1, evr);
              if evr.Err=0 then eval_mod(e^.RNode, m, v2, evr);
              if evr.Err=0 then mp_mulmod(v1,v2,m,r);
            end;
      _SQR: begin
              eval_mod(e^.LNode, m, v1, evr);
              if evr.Err=0 then mp_sqrmod(v1,m,r);
            end;
      else  begin
              eval(e, v1, evr);
              if evr.Err=0 then mp_mod(v1,m,r);
            end;
  end;
  if (evr.ERR=0) and (MP_Error<>MP_OKAY) then evr.Err := Err_MPERR_Eval;
  mp_clear2(v1,v2);
end;


{---------------------------------------------------------------------------}
procedure mp_clear_expr(var e: PExpr);
  {-Release memory used by e and clear mp_int values}
begin
  if e<>nil then with e^ do begin
    case nn of
         0:  if op=_CONST then mp_clear(Value);
         1:  if SNode<>nil then mp_clear_expr(SNode);
       else  begin
               if e^.LNode<>nil then mp_clear_expr(e^.LNode);
               if e^.RNode<>nil then mp_clear_expr(e^.RNode);
             end;
    end;
    mp_freemem(pointer(e),sizeof(TExpr));
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_eval(e: PExpr; var evr: TEval);
  {-Evaluate expression tree e, result in evr}
begin
  with evr do begin
    Err := 0;
    eval(e,Res,evr);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_eval_mod(e: PExpr; const m: mp_int; var evr: TEval);
  {-Evaluate expression tree e mod m, result in evr}
begin
  with evr do begin
    Err := 0;
    eval_mod(e,m,Res,evr);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_calculate(psz: pchar8; var evr: TEval; var EPos: integer);
  {-Parse and evaluate string psz}
var
  e: PExpr;
  pc: pchar8;
begin
  e := nil;
  pc := mp_parse(psz,e,evr.Err);
  EPos := pc-psz;
  if evr.Err=0 then mp_eval(e, evr);
  mp_clear_expr(e);
end;


{---------------------------------------------------------------------------}
function mp_calc_errorstr(Err: integer): mp_string;
  {-Translate known error codes}
var
  s: string[22];
begin
  case Err of
    Err_Missing_LeftBracket  : mp_calc_errorstr := 'Missing "("';
  {$ifdef Unicode}
    Err_Missing_Comma        : mp_calc_errorstr := mp_string('Missing argument separator ("'+mp_arg_sep+'")');
  {$else}
    Err_Missing_Comma        : mp_calc_errorstr := 'Missing argument separator ("'+mp_arg_sep+'")';
  {$endif}
    Err_Missing_RightBracket : mp_calc_errorstr := 'Missing ")"';
    Err_Unknown_Function     : mp_calc_errorstr := 'Unknown function';
    Err_Unknown_Element      : mp_calc_errorstr := 'Unknown element';
    Err_Trailing_Garbage     : mp_calc_errorstr := 'Trailing garbage';
    Err_Invalid_Number       : mp_calc_errorstr := 'Invalid number';
    Err_Invalid_HexNumber    : mp_calc_errorstr := 'Invalid hex number';
    Err_Invalid_RadixNumber  : mp_calc_errorstr := 'Invalid radix number';
    Err_Unknown_Operation    : mp_calc_errorstr := 'Unknown operation';
    Err_MPERR_Eval           : mp_calc_errorstr := 'MP_Err <> MP_OK';
    Err_Division_by_zero     : mp_calc_errorstr := 'Division by zero';
    Err_Invalid_argument     : mp_calc_errorstr := 'Invalid argument(s)';
    Err_X_not_init           : mp_calc_errorstr := 'Variable X not initialized';
    Err_Y_not_init           : mp_calc_errorstr := 'Variable Y not initialized';
    Err_Z_not_init           : mp_calc_errorstr := 'Variable Z not initialized';
    Err_no_solution          : mp_calc_errorstr := 'No solution';
    else begin
      str(Err,s);
      mp_calc_errorstr := 'mp_calc error '+s;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_init_eval(var evr: TEval);
  {-Initialize the mp_ints of evr}
begin
  with evr do mp_init4(X,Y,Z,Res);
  evr.Err := 0;
end;


{---------------------------------------------------------------------------}
procedure mp_clear_eval(var evr: TEval);
  {-Clear the mp_ints of evr}
begin
  with evr do mp_clear4(X,Y,Z,Res);
end;


end.

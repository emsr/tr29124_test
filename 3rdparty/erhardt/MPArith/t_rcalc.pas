{Test program for mp_rcalc unit, (c) W.Ehrhardt 2008-2018}

program t_rcalc;

{$i STD.INC}
{$i mp_conf.inc}

{$x+}  {pchar I/O}
{$i+}  {RTE on I/O error}

{$ifdef BIT16}
{$N+}
{$endif}

{$ifdef Delphi}
{$J+}
{$endif}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

(*
2008-01-17: prec [n], display prec
2008-04-19: help: added pi, removed !
2008-09-16: bin, adjust nfd if radix changes, help: log(a,b)
2008-09-25: HexLong, removed mem_util
2008-11-04: Initial Starttimer
2008-12-02: Uses BTypes
2009-01-09: display sign with ; line terminator
2009-07-07: D12 fixes, display new prec
2009-11-10: additional elementary transcendental functions
2010-06-03: display MPF_MAX_PREC related info
2010-01-09: XH command
2010-09-05: DDH, SCI, ALT commands
2011-06-21: ASX/ASD commands
2011-06-25: NFD command, updated help screen
2012-01-17: Fix XH for EXT64 (extended=double)
2012-08-02: updated help display
2013-01-12: DH command
2017-09-29: Radix command
2018-06-19: Soft ASX for EXT64, some ??_ext -> ??_dbl
2018-11-15: SH command
*)

{
TODO:
 - edit last line
}

uses
  BTypes, HRTimer,
  {$ifdef MPC_Diagnostic}
     mp_supp,
  {$endif}
  mp_types, mp_base, mp_real, mp_rcalc;


{---------------------------------------------------------------------------}
function HexLong(L: longint): mp_string;
  {-longint as hex string, LSB first}
var
  i: integer;
  s: string[8];
begin
  s := '';
  for i:=0 to 7 do begin
    s := mp_ucrmap[L and $F] + s;
    L := L shr 4;
  end;
  HexLong := s;
end;


{---------------------------------------------------------------------------}
function HexW(w: word): mp_string;
  {-longint as hex string, LSB first}
var
  i: integer;
  s: string[8];
begin
  s := '';
  for i:=0 to 3 do begin
    s := mp_ucrmap[w and $F] + s;
    w := w shr 4;
  end;
  HexW := s;
end;


{$ifndef EXT64}
{---------------------------------------------------------------------------}
function HexExtended(x: extended): mp_string;
  {-Extended as Hex array of word}
var
  xa: TMPHexExtW;
begin
  xa := TMPHexExtW(x);
  HexExtended := '($'+HexW(xa[0])+',$'+HexW(xa[1])+',$'+HexW(xa[2])+',$'+HexW(xa[3])+',$'+HexW(xa[4])+')';
end;
{$endif}


{---------------------------------------------------------------------------}
function HexDoubleDouble(ddh,ddl: double): mp_string;
var
  h: TMPHexDblW;
  l: TMPHexDblW;
begin
  h := TMPHexDblW(ddh);
  l := TMPHexDblW(ddl);
  HexDoubleDouble := '($'+HexW(l[0])+',$'+HexW(l[1])+',$'+HexW(l[2])+',$'+HexW(l[3])
                    +',$'+HexW(h[0])+',$'+HexW(h[1])+',$'+HexW(h[2])+',$'+HexW(h[3])+')';
end;


{---------------------------------------------------------------------------}
function HexDouble(d: double): mp_string;
var
  da: TMPHexDblW;
begin
  da := TMPHexDblW(d);
  HexDouble := '($'+HexW(da[0])+',$'+HexW(da[1])+',$'+HexW(da[2])+',$'+HexW(da[3])+')';
end;


{---------------------------------------------------------------------------}
function HexSingle(s: single): mp_string;
var
  si: longint absolute s;
begin
  HexSingle := HexLong(si);
end;


{---------------------------------------------------------------------------}
procedure ShowInfo;
begin
  writeln;
{$ifdef MPC_E1Ln10Tab}
  writeln('Operators:  +  -  *  /  ^             Constants: pi, e, ln2, ln10');
{$else}
  writeln('Operators:  +  -  *  /  ^             Constants: pi, ln2');
{$endif}
  writeln;
  writeln('Functions:  abs(a)       agm(a,b)     arccos(a)     arccosh(a)   arccosh1p(a)');
  writeln('            arccot(a)    arccotc(a)   arccoth(a)    arccsc(a)    arccsch(a)');
  writeln('            arcgd(a)     archav(a)    arcsec(a)     arcsech(a)   arcsin(a)');
  writeln('            arcsinh(a)   arctan(a)    arctan2(y,x)  arctanh(a)   asd(a)');
  writeln('            ass(a)       asx(a)       CE(a)         CK(a)        cos(a)');
  writeln('            cosh(a)      coshm1(a)    cot(a)        coth(a)      csc(a)');
  writeln('            csch(a)      EE(a)        EK(a)         exp(a)       expm1(a)');
  writeln('            frac(a)      gd(a)        hav(a)        hypot(a,b)   int(a)');
  writeln('            lambertw(a)  ln(a)        ln1p(a)       log10(a)     log2(a)');
  writeln('            log(a,b)     max(a,b)     min(a,b)      nroot(a,n)   numbpart(a)');
  writeln('            predd(a)     preds(a)     predx(a)      random(a)    round(a)');
  writeln('            sec(a)       sech(a)      sin(a)        sinc(a)      sinh(a)');
  writeln('            sqr(a)       sqrt(a)      succd(a)      succs(a)     succx(a)');
  writeln('            tan(a)       tanh(a)      vers(a)');
  writeln;
  writeln('Variables:  x  y  z,  Syntax: Var=expression[;] or Var=_[;]');
  writeln;
  writeln('Other    :  line terminator ";" shows only sign/ldx/chksum/time of result');
  writeln('            bin, dec, hex: set radix for result display');
  writeln('            xh,sh,dh, ddh: last result as hex extended/single/double(double)');
  writeln('            sci,  alt:  display results using scientific/alternative format');
  writeln('            prec  [n]:  display/set bit precision');
  writeln('            nfd   [n]:  display/set number of digits in fractional part');
  writeln('            radix [n]:  set radix for result display');
  writeln('            as[x|d|s]:  set result to extended/double/single value');
  writeln('            "_<enter>"  re-displays last result');
{$ifdef CPUARM}
  writeln('            ".<enter>"  displays time for last calculation; the ARM version');
  writeln('                        does not use TSC and is not high-precision.');
{$else}
  writeln('            ".<enter>"  displays time for last calculation');
{$endif}
  writeln;
end;


{---------------------------------------------------------------------------}
procedure HelpLine;
begin
  writeln('Type "?<enter>" to get some info about commands, "\q" or "quit" to end.');
end;


const
  nfd0: word = 60;
  maxnfd = 200;
  MaxSglHex: packed array[0..3] of byte = ($ff,$ff,$7f,$7f);  {3.4028234E+38}

var
  evr: TFEval;
{$ifdef VirtualPascal}
  ir: longint;
{$else}
  ir: integer;
{$endif}
  EP,i: integer;
  HR: THRTimer;
  dp,ctime,ddh,ddl: double;
  ts: single;
  ac: char8;
  use_sci,print,script: boolean;
  s: string[255];
  ss: string[2];
  cmd: string[20];
  radix,nfd: word;
  rc: string[2];
  prec,n: longint;
  tmp: mp_float;
{$ifndef EXT64}
  tx: extended;
{$endif}
begin

  {mp_verbose := 3;}

  StartTimer(HR);
  mpf_init_eval(evr);
  mpf_init(tmp);

  mp_uppercase := true;
  use_sci := false;
  script := false;
  ctime  := 0.0;
  rc := 'D';
  radix := 10;
  nfd := nfd0;

  for i:=1 to paramcount do begin
    if (paramstr(i)='/s') or (paramstr(i)='/S') then script := true;
  end;

  writeln('Test of MPArith V', MP_VERSION, ' [mp_rcalc]   (c) W.Ehrhardt 2008-2018');
  writeln('Karatsuba  cutoffs: mul/sqr = ',mp_mul_cutoff,'/',mp_sqr_cutoff);
  writeln('Toom-3, BZ cutoffs: mul/sqr = ',mp_t3m_cutoff,'/',mp_t3s_cutoff,  ',  div = ',mp_bz_cutoff);
  if not script then HelpLine;
  prec := mpf_get_default_prec;
  dp :=prec*ln(2)/ln(10);
  writeln('Current bit precision = ', prec, ' (max:',MPF_MAX_PREC, '),  decimal precision = ', dp:1:1);
  writeln;

  repeat
    ac := #0;
    if not script then write('[',rc,']:=> ');
    readln(s);
    while (s<>'') and (s[1]=' ') do delete(s,1,1);
    while (s<>'') and (s[length(s)]=' ') do delete(s,length(s),1);
    if s='' then continue;

    if s='?' then begin
      ShowInfo;
      continue;
    end;
    if s='.' then begin
      writeln('Time = ', ctime:1:3, ' ms');
      continue;
    end;
    if s[1]='"' then begin
      delete(s,1,1);
      writeln(s);
      continue;
    end;

    cmd := s;
    for i:=1 to length(cmd) do cmd[i] := upcase(cmd[i]);

    if (cmd='\Q') or (cmd='Q') or (cmd='QUIT') or (cmd='\@') then break;
    if cmd='HEX' then begin
      radix := 16;
      rc := 'H';
      dp := prec/4.0;
      if dp<nfd0 then nfd := trunc(dp) else nfd := nfd0;
      if script then writeln(cmd);
      continue;
    end;
    if cmd='DEC' then begin
      radix := 10;
      rc := 'D';
      dp := prec*ln(2)/ln(10);
      if dp<nfd0 then nfd := trunc(dp) else nfd := nfd0;
      if script then writeln(cmd);
      continue;
    end;
    if cmd='BIN' then begin
      radix := 2;
      rc := 'B';
      if prec<nfd0 then nfd := prec else nfd := nfd0;
      if script then writeln(cmd);
      continue;
    end;
    if copy(cmd,1,5)='RADIX' then begin
      if length(cmd)=5 then begin
        writeln('Radix = ', radix);
        continue;
      end;
      delete(cmd,1,5);
      while (cmd<>'') and (cmd[1]=' ') do delete(cmd,1,1);
      {$ifdef D12Plus}
        val(string(cmd),n,i);
      {$else}
        val(cmd,n,i);
      {$endif}
      if (i=0) and (n>1) and (n<=MaxRadix) then begin
        radix := n;
        str(n,rc);
      end
      else writeln('Invalid radix: 2 .. ', Maxradix);
      continue;
    end;
    if (cmd='ALT') or (cmd='SCI') then begin
      use_sci := cmd='SCI';
      continue;
    end;
    if copy(cmd,1,4)='PREC' then begin
      if length(cmd)=4 then begin
        prec := mpf_get_default_prec;
        dp := prec*ln(2)/ln(10);
        writeln('Current bit precision = ', prec, ',  decimal precision = ', dp:1:1);
        continue;
      end;
      delete(cmd,1,4);
      while (cmd<>'') and (cmd[1]=' ') do delete(cmd,1,1);
      {$ifdef D12Plus}
        val(string(cmd),prec,ir);
      {$else}
        val(cmd,prec,ir);
      {$endif}
      if (ir=0) and (prec>=MPF_MIN_PREC) then begin
        if prec>MPF_MAX_PREC then begin
          prec := MPF_MAX_PREC;
          write('** prec > ',MPF_MAX_PREC, ' - ');
        end;
        mpf_set_default_prec(prec);
        prec := mpf_get_default_prec;
        dp :=prec*ln(2)/ln(10);
        if dp<nfd0 then nfd := trunc(dp) else nfd := nfd0;
        mpf_chg_prec(evr.Res, prec);
        writeln('New bit precision = ', prec, ',  decimal precision = ', dp:1:1);
        continue;
      end;
    end;
    if copy(cmd,1,3)='NFD' then begin
      if length(cmd)=3 then begin
        writeln('Number of digits in fractional part = ', nfd0);
        continue;
      end;
      delete(cmd,1,3);
      while (cmd<>'') and (cmd[1]=' ') do delete(cmd,1,1);
      {$ifdef D12Plus}
        val(string(cmd),n,ir);
      {$else}
        val(cmd,n,ir);
      {$endif}
      if (ir=0) and (n>=10) and (n<=maxnfd) then begin
        nfd := n;
        nfd0 := n;
        writeln('Number of digits in fractional part = ', nfd0);
      end
      else writeln('Invalid number [10 <= nfd <= ',maxnfd);
      continue;
    end;
    if cmd='ASX' then begin
      {$ifdef EXT64}
        {Software simulation for 64 mantissa}
        {Todo: denormal with bitprec < 64}
        n := Evr.Res.bitprec;
        s_mpf_normalizep(Evr.Res, 64);
        if Evr.Res.exponent > 16320 then Evr.Err:= Err_Overflow
        else if Evr.Res.exponent < -16508 then mpf_set0(Evr.Res);
        s_mpf_normalizep(Evr.Res, n);
      {$else}
        if (Evr.Err=0) and (Evr.Res.bitprec>64) then begin
          tx := mpf_toextended(Evr.Res);
          if abs(tx)<>DblPosInf then begin
            mpf_set_ext(Evr.Res,tx);
          end
          else Evr.Err:= Err_Overflow;
        end;
      {$endif}
      s := '_';
    end;
    if cmd='ASD' then begin
      if (Evr.Err=0) and (Evr.Res.bitprec>64) then begin
        ddh := mpf_todouble(Evr.Res);
        if abs(ddh)<>DblPosInf then begin
          mpf_set_dbl(Evr.Res,ddh);
        end
        else Evr.Err:= Err_Overflow;
      end;
      s := '_';
    end;
    if cmd='ASS' then begin
      if (Evr.Err=0) and (Evr.Res.bitprec>32) then begin
        ddh := mpf_todouble(Evr.Res);
        if (abs(ddh)<>DblPosInf) and (abs(ddh)<Single(MaxSglHex)) then begin
          ts := ddh;
          mpf_set_dbl(Evr.Res,ts);
        end
        else Evr.Err:= Err_Overflow;
      end;
      s := '_';
    end;
{$ifndef EXT64}
    if cmd='XH' then begin
      if (Evr.Err=0) and (Evr.Res.bitprec>=64) then begin
        tx := mpf_toextended(Evr.Res);
        if abs(tx)<> DblPosInf then begin
          mpf_set_ext(tmp,tx);
          write('Result = ', HexExtended(tx));
          writeln('  {',mpf_decimal(tmp,20),'}');
        end
        else writeln('** INF');
      end;
      continue;
    end;
{$endif}
    if cmd='DH' then begin
      if (Evr.Err=0) and (Evr.Res.bitprec>=53) then begin
        ddh := mpf_todouble(Evr.Res);
        if abs(ddh)<>DblPosInf then begin
          write('Result = ', HexDouble(ddh));
          writeln('  {',ddh,'}');
        end
        else writeln('** INF');
      end;
      continue;
    end;
    if cmd='SH' then begin
      if (Evr.Err=0) and (Evr.Res.bitprec>=24) then begin
        ddh := mpf_todouble(Evr.Res);
        if (abs(ddh)<>DblPosInf) and (abs(ddh)<Single(MaxSglHex)) then begin
          ts := ddh;
          write('Result = ', HexSingle(ts));
          writeln('  {',ts,'}');
        end
        else writeln('** INF');
      end;
      continue;
    end;
    if cmd='DDH' then begin
      if (Evr.Err=0) and (Evr.Res.bitprec>=106) then begin
        ddh := mpf_todouble(Evr.Res);
        if abs(ddh)<>DblPosInf then begin
          mpf_set_dbl(tmp,ddh);
          mpf_sub(Evr.Res,tmp,tmp);
          ddl := mpf_todouble(tmp);
          write('Result = ', HexDoubleDouble(ddh,ddl));
          writeln('  {',ddh,',', ddl,'}');
        end
        else writeln('** INF');
      end;
      continue;
    end;

    {Echo input in script mode}
    if script then writeln('[Expr] = ',s);

    {Check for trailing ";", i.e. dont print the result}
    if (s<>'') and (s[length(s)]=';') then begin
      delete(s,length(s),1);
      while (s<>'') and (s[length(s)]=' ') do delete(s,length(s),1);
      print := false;
    end
    else print := true;
    if s='' then continue;

    if s<>'_' then begin
      s := s + #0;
      ac := #0;
      if (length(s)>1) and (upcase(s[1]) in ['X','Y','Z']) then begin
        {Check if we have an assignment to x,y,z. Analyse first non-blank}
        {char; next while/if statements are OK because s has trailing #0}
        while (s[2]=' ') and (s[3]=' ') do delete(s,2,1);
        if (s[2]=' ') and (s[3]='=') then delete(s,2,1);
        if s[2]='=' then begin
          ac := upcase(s[1]);
          delete(s,1,2);
          while (s<>'') and (s[1]=' ') do delete(s,1,1);
        end;
      end;
      if s<>'_'#0 then begin
        RestartTimer(HR);
        mpf_calculate(@s[1],evr,EP);
        ctime := 1000.0*ReadSeconds(HR);
        if evr.Err=Err_MPERR_Eval then set_mp_error(MP_OKAY);
      end;
    end
    else if evr.Err>0 then continue;

    if evr.Err>0 then begin
      writeln('Error ', evr.Err,', ',mpf_calc_errorstr(evr.Err),':  <',copy(s,EP+1,length(s)-EP-1), '>');
      if evr.Err=Err_Unknown_Function then HelpLine;
    end
    else if evr.Err<0 then writeln(#7'Error: ', evr.Err, ' [',mpf_calc_errorstr(evr.Err),']')
    else begin
      if ac in ['X', 'Y', 'Z'] then begin
        case ac of
          'X':  mpf_copy(evr.Res,evr.X);
          'Y':  mpf_copy(evr.Res,evr.Y);
          'Z':  mpf_copy(evr.Res,evr.Z);
        end;
        write(ac,' = ');
      end
      else write('Result = ');
      if print then begin
        if use_sci then mpf_output_radix(evr.Res, radix, nfd)
        else mpf_output_radix_alt(evr.Res, radix, nfd);
      end
      else begin
        case mp_sign(evr.Res.Mantissa) of
          +1 : ss := '>0';
           0 : ss := '=0';
          else ss := '<0'
        end;
        write('[ ',ss,',  ldx=',s_mpf_ldx(evr.Res), ',  chksum=$',
              HexLong(mpf_checksum(evr.Res)),',  time=', ctime:1:3, ' ms]');
      end;
    end;
    evr.Err := 0;
    writeln;
  until false;
  mpf_clear(tmp);
  mpf_clear_eval(evr);
  {$ifdef MPC_Diagnostic}
     mp_dump_meminfo;
     mp_dump_diagctr;
  {$endif}
end.


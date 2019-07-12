{Test program for mp_ccalc unit, (c) W.Ehrhardt 2018}

program t_ccalc;

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

{$ifndef MPC_NOHALT}
  {$message error 't_ccalc should be compiled with -dMPC_NOHALT'}
{$endif}


uses
  BTypes, HRTimer,
  {$ifdef MPC_Diagnostic}
     mp_supp,
  {$endif}
  mp_types, mp_base, mp_real, mp_cmplx, mp_ccalc;


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
procedure ShowInfo;
begin
  writeln;
{$ifdef MPC_E1Ln10Tab}
  writeln('Operators:  +  -  *  /  ^             Constants: i, pi, e, ln2, ln10');
{$else}
  writeln('Operators:  +  -  *  /  ^             Constants: i, pi, ln2');
{$endif}
  writeln;
  writeln('Functions:  abs(a)       agm(a,b)     arccos(a)    arccosh(a)   arccot(a)');
  writeln('            arccotc(a)   arccoth(a)   arccsc(a)    arccsch(a)   arcsec(a)');
  writeln('            arcsech(a)   arcsin(a)    arcsinh(a)   arctan(a)    arctanh(a)');
  writeln('            arg(a)       conj(a)      cos(a)       cosh(a)      cot(a)');
  writeln('            coth(a)      csc(a)       csch(a)      exp(a)       expm1(a)');
  writeln('            im(a)        ln(a)        ln1p(a)      log10(a)     nroot(a,n)');
  writeln('            re(a)        sec(a)       sech(a)      sin(a)       sinh(a)');
  writeln('            sqr(a)       sqrt(a)      tan(a)       tanh(a)');
  writeln;
  writeln('Variables:  x  y  z,  Syntax: Var=expression[;] or Var=_[;]');
  writeln;
  writeln('Other    :  line terminator ";" shows only chksum/time of result');
  writeln('            sci,  alt:  display results using scientific/alternative format');
  writeln('            prec  [n]:  display/set bit precision');
  writeln('            nfd   [n]:  display/set number of digits in fractional part');
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
  nfd0: word = 50;
  maxnfd = 200;
var
  evr: TFEval;
{$ifdef VirtualPascal}
  ir: longint;
{$else}
  ir: integer;
{$endif}
  EP,i: integer;
  HR: THRTimer;
  dp,ctime: double;
  ac: char8;
  use_sci,print: boolean;
  s: string[255];
  resprefix: string[8];
  cmd: string[20];
  nfd: word;
  prec,n: longint;
begin

  {mp_verbose := 3;}

  StartTimer(HR);
  mpc_init_eval(evr);

  mp_uppercase := true;
  use_sci := false;
  ctime := 0.0;
  nfd := nfd0;

  writeln('Test of MPArith V', MP_VERSION, ' [mp_ccalc]   (c) W.Ehrhardt 2008-2018');
  writeln('Karatsuba  cutoffs: mul/sqr = ',mp_mul_cutoff,'/',mp_sqr_cutoff);
  writeln('Toom-3, BZ cutoffs: mul/sqr = ',mp_t3m_cutoff,'/',mp_t3s_cutoff,  ',  div = ',mp_bz_cutoff);
  HelpLine;
  prec := mpf_get_default_prec;
  dp :=prec*ln(2)/ln(10);
  writeln('Current bit precision = ', prec, ' (max:',MPF_MAX_PREC, '),  decimal precision = ', dp:1:1);
  writeln;

  repeat
    ac := #0;
    write('=> ');
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
        mpf_chg_prec(evr.Res.re, prec);
        mpf_chg_prec(evr.Res.im, prec);
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
        mpc_calculate(@s[1],evr,EP);
        ctime := 1000.0*ReadSeconds(HR);
        if evr.Err=Err_MPERR_Eval then set_mp_error(MP_OKAY);
      end;
    end
    else if evr.Err>0 then continue;

    if evr.Err>0 then begin
      writeln('Error ', evr.Err,', ',mpc_calc_errorstr(evr.Err),':  <',copy(s,EP+1,length(s)-EP-1), '>');
      if evr.Err=Err_Unknown_Function then HelpLine;
    end
    else if evr.Err<0 then writeln(#7'Error: ', evr.Err, ' [',mpc_calc_errorstr(evr.Err),']')
    else begin
      resprefix := 'Result';
      if ac in ['X', 'Y', 'Z'] then begin
        case ac of
          'X':  mpc_copy(evr.Res,evr.X);
          'Y':  mpc_copy(evr.Res,evr.Y);
          'Z':  mpc_copy(evr.Res,evr.Z);
        end;
        resprefix := ac;
      end;
      if print then begin
        if use_sci then begin
          write(resprefix,'.re = '); mpf_output_decimal(evr.Res.re, nfd); writeln;
          write(resprefix,'.re = '); mpf_output_decimal(evr.Res.im, nfd); writeln;
        end
        else begin
          write(resprefix,'.re = '); mpf_output_decimal_alt(evr.Res.re, nfd); writeln;
          write(resprefix,'.im = '); mpf_output_decimal_alt(evr.Res.im, nfd); writeln;
        end;
      end
      else begin
        write('[ chksum=$', HexLong(mpc_checksum(evr.Res)),',  time=', ctime:1:3, ' ms]');
      end;
    end;
    evr.Err := 0;
    writeln;
  until false;
  mpc_clear_eval(evr);
  {$ifdef MPC_Diagnostic}
     mp_dump_meminfo;
     mp_dump_diagctr;
  {$endif}
end.


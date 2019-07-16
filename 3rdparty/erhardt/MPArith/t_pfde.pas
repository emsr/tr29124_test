{-Test prog for factorization, we 22.08.04}

{
 Vers.  Date      Author      Modification
 -----  --------  -------     ------------------------------------------
  0.10  09.08.04  W.Ehrhardt  Initial version, incl Pollard rho and (p-1)
  0.20  15.08.05  we          William's (p+1)
  0.30  17.09.05  we          Brent's ECM
  0.31  17.09.05  we          uses unit pfdu
  0.32  22.09.06  we          with mp_calc
  0.33  01.01.12  we          for D12Plus
}

program t_pfde;

{$i STD.INC}

{$x+,i+}

{$ifndef FPC}
{$N+}
{$endif}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  crt,
  mp_types, mp_base, mp_supp, mp_calc, mp_numth,
  pfdu, pfdu_crt, HRTimer;
var
  evr: TEval;
  s: mp_string;
  EP: integer;
  ctx: pfd_ctx;
  HR: THRTimer;
begin
{
  PP1Bnd := 0;
  PM1Bnd := 0;
  RhoCnt := 0;
  PBrmax := 0;
}
  StartTimer(HR);
  pfd_initialize(ctx);
  mp_init_eval(evr);
  textmode(CO80+Font8x8);
  textattr := green;
  writeln('Test of MPArith V', MP_VERSION, ' [mp_calc/pfdu]  (c) W.Ehrhardt 2006-2014');
  writeln('To exit the test program entering an empty expressions string.');
  writeln;
  {s := '1606938044258990275541962091078582604694884852567413570404019';}{Fermat/RSA}
  {s := '18548676741817250104151622545580576823736636896432849057'+
        '10984160646722888555430591384041316374473729421512365598'+
        '29709849969346650897776687202384767704706338162219624578'+
        '777915220190863619885201763980069247978050169295918863';}  {HOLF}
  s := {$ifdef D12Plus}mp_string{$endif}(paramstr(1));
  if s='/+240' then begin
    lim240 := false;
    s := '';
  end;
  while mp_error=MP_OKAY do begin
    if s='' then begin
      textattr := cyan;
      write('Expression to factor: ');
      textattr := lightgray;
      readln(s);
      if (s='') or (s='q') then break;
    end;
    textattr := lightgray;
    s := s + #0;
    mp_calculate(@s[1],evr,EP);
    if evr.Err>0 then writeln('Error ', evr.Err,' at pos ', EP+1,'  <',copy(s,EP+1,length(s)-EP-1), '>')
    else if evr.Err<0 then writeln('Error: ', evr.Err)
    else begin
      write('N=');
      mp_output_decimal(evr.Res);
      writeln;
      pfd_reset(ctx);
      RestartTimer(HR);
      pfd_factor(ctx, evr.Res);
      writeln(' Time: ',ReadSeconds(HR):1:3, 's');
    end;
    {$ifndef VirtualPascal}
      Sound(800);
      delay(100);
      Nosound;
    {$endif}
    s := '';
  end;
  mp_clear_eval(evr);
  pfd_finalize(ctx);
  {$ifdef MPC_Diagnostic}
    mp_dump_meminfo;
    mp_dump_diagctr;
  {$endif}
end.

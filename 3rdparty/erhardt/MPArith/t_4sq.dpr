{Demo program for MPArith library, (c) W.Ehrhardt 2012}

program t_4sq;

{$i STD.INC}
{$i mp_conf.inc}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

{$ifndef FPC}
  {$N+}
{$endif}

uses
  hrtimer,
  mp_types, mp_base, mp_numth;

var
  a,b,c,d,e: mp_int;
  HR: THRTimer;

begin
  mp_init5(a,b,c,d,e);
  writeln('Test of MP library version ', MP_VERSION, '   (c) W.Ehrhardt 2012');
  StartTimer(HR);
  if paramcount=0 then mp_lucas(1000,a)
  else mp_read_decimal_astr(a,{$ifdef D12Plus}mp_string{$endif}(paramstr(1)));
  write(mp_adecimal(a),'=');
  mp_4sq_sd(a,b,c,d,e);
  write(mp_adecimal(b),'^2 + ');
  write(mp_adecimal(c),'^2 + ');
  write(mp_adecimal(d),'^2 + ');
  writeln(mp_adecimal(e),'^2');
  writeln('Time [s]: ',ReadSeconds(HR):1:3);
  mp_clear5(a,b,c,d,e);
end.

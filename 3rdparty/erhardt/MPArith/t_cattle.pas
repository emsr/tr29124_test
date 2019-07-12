{Demo program for MPArith library, (c) W.Ehrhardt 2012}
{Solution of the Cattle Problem of Archimedes         }
{Needs MP_MAXBIT > 700000, i.e. a 32/64 bit compiler  }

program t_cattle;

{$i STD.INC}
{$i mp_conf.inc}

{$x+}  {pchar I/O}
{$i+}  {RTE on I/O error}

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
  i: longint;
  HR: THRTimer;
  st: ansistring;
  tf: text;
begin

  {Ref: H.C.Williams, R.A.German, C.R.Zarnke, Solution of the Cattle Problem of}
  {Archimedes, Math.Comp. 19 (1965), 671-674, http://dx.doi.org/10.2307/2003954}

  {See also: http://en.wikipedia.org/wiki/Archimedes%27_cattle_problem}

  writeln('Demo of MP library version ', MP_VERSION, '   (c) W.Ehrhardt 2012');
  writeln('** Solution of the Cattle Problem of Archimedes *');

  if longint(MAXDigits)*lv_digit_bit < 700000 then begin
    writeln('MP_MAXBIT < 700000');
    halt;
  end;
  mp_init5(a,b,c,d,e);

  StartTimer(HR);
  writeln('Calculation ...');
  mp_set_int(c,2*3*7*11*29*353);
  mp_pell(c,a,b);
  mp_powerd(a,b,c,2329,d,e);
  mp_div_int(e,2*4657,@e,i);
  mp_read_decimal(a,'224571490814418'); {'...'=2*3*11*29*41*107*4657*5743}
  mp_sqr(e,d);
  mp_mul(d,a,c);
  writeln('Calculation time [s]: ',ReadSeconds(HR):1:3);
  RestartTimer(HR);

  writeln('Formatting and I/O ...');
  st := mp_adecimal(c);
  assign(tf,'t_cattle.txt');
  rewrite(tf);
  writeln(tf,'Archimedes'' cattle number T = ');
  writeln(tf, st);
  close(tf);

  writeln('Archimedes'' cattle number T (206000+ digits) = ');
  writeln(copy(st,1,36),'...',copy(st,length(st)-35,36));
  writeln('Complete T written to t_cattle.txt');
  writeln('Formatting and I/O time [s]: ',ReadSeconds(HR):1:3);

  mp_clear5(a,b,c,d,e);
end.

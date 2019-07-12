{Simple test for pasrand unit, we 12.2017}

program t_rnd_06;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
{$endif}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  hrtimer,
  pasrand;

{$ifdef BIT16}
const
  NMAX = 5000000;
{$else}
const
  NMAX = 100000000;
{$endif}

var
  n: longint;
  rng: pasrand_ctx;
  HR: THRTimer;
  cps,sec: double;
begin
  writeln('Test for pasrand unit: ', NMAX, ' pasrand_next calls     (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', pasrand_selftest);
  pasrand_init(rng,1);
  StartTimer(HR);
  for n:=1 to NMAX do begin
    pasrand_next(rng);
  end;
  sec := readseconds(HR);
  cps := NMAX/sec;
  writeln('Done. Time ', sec:1:3, 's,  calls/sec ', cps:1:0);
end.

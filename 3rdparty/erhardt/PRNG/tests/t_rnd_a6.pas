{Simple test for well1024a unit, we 12.2017}

program t_rnd_96;

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
  well1024;

{$ifdef BIT16}
const
  NMAX = 5000000;
{$else}
const
  NMAX = 100000000;
{$endif}

var
  n: longint;
  rng: well1024a_ctx;
  HR: THRTimer;
  cps,sec: double;
begin
  writeln('Test for well1024a unit: ', NMAX, ' well1024a_next calls     (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', well1024a_selftest);
  well1024a_init(rng,1);
  StartTimer(HR);
  for n:=1 to NMAX do begin
    well1024a_next(rng);
  end;
  sec := readseconds(HR);
  cps := NMAX/sec;
  writeln('Done. Time ', sec:1:3, 's,  calls/sec ', cps:1:0);
end.

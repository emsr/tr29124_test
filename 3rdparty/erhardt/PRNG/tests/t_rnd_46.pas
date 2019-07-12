{Simple test for mt19937 unit, we Aug.2017}

program t_rnd_46;

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
  mt19937;

{$ifdef BIT16}
const
  NMAX = 5000000;
{$else}
const
  NMAX = 100000000;
{$endif}

var
  n: longint;
  rng: mt19937_ctx;
  HR: THRTimer;
  cps,sec: double;
begin
  writeln('Test for mt19937 unit: ', NMAX, ' mt19937_next calls     (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', mt19937_selftest);
  mt19937_init(rng,1);
  StartTimer(HR);
  for n:=1 to NMAX do begin
    mt19937_next(rng);
  end;
  sec := readseconds(HR);
  cps := NMAX/sec;
  writeln('Done. Time ', sec:1:3, 's,  calls/sec ', cps:1:0);
end.

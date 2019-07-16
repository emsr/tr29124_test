{Simple test for tt800 unit, we Aug.2017}

program t_rnd_26;

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
  tt800;

{$ifdef BIT16}
const
  NMAX = 5000000;
{$else}
const
  NMAX = 100000000;
{$endif}

var
  n: longint;
  rng: tt800_ctx;
  HR: THRTimer;
  cps,sec: double;
begin
  writeln('Test for tt800 unit: ', NMAX, ' tt800_next calls     (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', tt800_selftest);
  tt800_init(rng,1);
  StartTimer(HR);
  for n:=1 to NMAX do begin
    tt800_next(rng);
  end;
  sec := readseconds(HR);
  cps := NMAX/sec;
  writeln('Done. Time ', sec:1:3, 's,  calls/sec ', cps:1:0);
end.

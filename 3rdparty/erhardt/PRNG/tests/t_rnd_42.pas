{Simple test for mt19937 unit, we July 2005}

program t_rnd_42;

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
  MT19937;

{$ifdef BIT32}
const
  NMAX = 100000000;
{$else}
const
  NMAX = 100000000 div 50;
{$endif}
var
  n: longint;
  rng: mt19937_ctx;
begin
  writeln('Test for mt19937 unit: ', NMAX, ' mt19937_next calls     (c) 2005 W.Ehrhardt');
  mt19937_init0(rng);
  for n:=1 to NMAX do begin
    mt19937_next(rng);
  end;
  writeln('Done');
end.

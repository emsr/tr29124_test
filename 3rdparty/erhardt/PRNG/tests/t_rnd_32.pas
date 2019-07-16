{Simple test for ISAAC unit, we July 2005}

program t_rnd_32;

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
  isaac;

{$ifdef BIT32}
const
  NMAX = 100000000;
{$else}
const
  NMAX = 100000000 div 20;
{$endif}
var
  n: longint;
  rng: isaac_ctx;
begin
  writeln('Test for ISAAC unit: ', NMAX, ' isaac_next calls     (c) 2005 W.Ehrhardt');
  isaac_init0(rng);
  for n:=1 to NMAX do begin
    isaac_next(rng);
  end;
  writeln('Done');
end.

{Simple test for salsar unit, we Apr.2006}

program t_rnd_82;

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
  salsar;

const
  NMAX = 100000000;
var
  n: longint;
  rng: salsar_ctx;
begin
  writeln('Test for salsar unit: ', NMAX, ' salsar_next calls     (c) 2006 W.Ehrhardt');
  salsar_init(rng,1);
  for n:=1 to NMAX do begin
    salsar_next(rng);
  end;
  writeln('Done');
end.

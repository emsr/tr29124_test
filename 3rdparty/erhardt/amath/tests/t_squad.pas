{Test program for AMTools/quadsolve,  (c) W.Ehrhardt 2010}
program t_squad;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}


uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_squadm;

begin
  test_quadsolve;
  test_cubsolve;
end.

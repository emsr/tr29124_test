{Test program for AMTools/quadsolve,  (c) W.Ehrhardt 2010}
program t_root;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}


uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_rootm;

begin
  test_quadsolve;
  test_cubsolve;
  test_PolyRoots;
end.

{Test program for DAMTools/quadsolve,  (c) W.Ehrhardt 2013}
program t_quadt;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}


uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_dquadm;

begin
  test_quadsolve;
  test_cubsolve;
end.

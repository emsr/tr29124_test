{Main test unit for DAMTools quadsolve/Polyroots  (c) W.Ehrhardt 2013-2017}
program t_root;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}


uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_drootm;

begin
  test_quadsolve;
  test_cubsolve;
  test_PolyRoots;
end.

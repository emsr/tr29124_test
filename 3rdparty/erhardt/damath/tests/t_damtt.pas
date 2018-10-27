{Test program DAMTools: Root finding, function minimization, integration, ...  (c) W.Ehrhardt 2010}
program t_damtt;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_damttm;
begin
  test_amtools;
  test_amint;
  test_amconvacc;
end.

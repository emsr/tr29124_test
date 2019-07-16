{T_DAMATH - regression test program for DAMATH unit  (c) 2013  W.Ehrhardt}

program t_damath;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_damatm;

begin
  test_damath_main;
end.

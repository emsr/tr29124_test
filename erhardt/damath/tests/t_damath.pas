{T_DAMATH - regression test program for DAMATH unit  (c) 2013  W.Ehrhardt}

{Many test cases were calculated with Maple V R4, Digits:=20}
{others with Pari/GP 2.3.4 and \p 20, some with Cephes qcalc}

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

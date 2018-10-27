{T_AMATH - regression test program for AMATH unit  (c) 2009-2011  W.Ehrhardt}

{Many test cases were calculated with Maple V R4, Digits:=20}
{others with Pari/GP 2.3.4 and \p 20, some with Cephes qcalc}

program t_amath;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_amathm;

begin
  test_amath_main;
end.

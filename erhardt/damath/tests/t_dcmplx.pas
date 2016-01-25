{Test/dev program for DAMCmplx  (c) W.Ehrhardt 2013-2014}
program t_dcmplx;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

{$ifdef BIT16}
{$N+,F+}
{$endif}

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  DAMath, DAMCmplx, t_dcmplu;

begin
  test_all_complex;
end.

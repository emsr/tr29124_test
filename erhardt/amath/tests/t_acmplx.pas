{Test/dev program for AMCmplx  (c) W.Ehrhardt 2013-2014}
program t_acmplx;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  AMath, AMCmplx, t_acmplu;
begin
  test_all_complex;
end.

{Test/dev program for AMath/FPU  (c) W.Ehrhardt 2017}
program t_fpu;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

{$ifdef BIT16}
{$N+,X+}
{$endif}

uses
  t_fpum;

begin
  test_fpu;
end.


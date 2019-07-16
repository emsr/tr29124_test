{Test/dev program for DFPU  (c) W.Ehrhardt 2017}
program t_dfpu;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

{$ifdef BIT16}
{$N+,X+}
{$endif}

uses
  t_dfpum;

begin
  test_fpu;
end.


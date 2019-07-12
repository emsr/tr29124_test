{Test program for AMQuat,  (c) W.Ehrhardt 2017}
program t_quat;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}


uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_quatm;

begin
  test_quat;
end.

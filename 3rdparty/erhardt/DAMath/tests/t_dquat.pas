{Test program for DAMQuat,  (c) W.Ehrhardt 2017}
program t_dquat;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}


uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_dquatm;

begin
  test_quat;
end.

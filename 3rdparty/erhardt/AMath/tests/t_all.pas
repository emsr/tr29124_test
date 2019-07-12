{T_ALL - all-in-one test program for AMATH package  (c) 2011-2017  W.Ehrhardt}

{For 32 bit, code segment too large for 16 bit}

program t_all;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

uses
  t_amathm, t_sfx_m, t_sfd_m, t_rootm, t_amtoom, t_acmplu, t_quatm, t_fpum;


procedure write_sep;
begin
  writeln;
  writeln('********************************************************');
  writeln('********************************************************');
end;


begin
  writeln('--------------------------------------------------');
  writeln('All-in-one test AMath package  (c) W.Ehrhardt 2018');
  writeln('--------------------------------------------------');

  test_amath_main;
  write_sep;

  test_sfx_main;
  write_sep;

  test_sfd_main;
  write_sep;

  test_quadsolve;
  test_cubsolve;
  test_PolyRoots;
  write_sep;

  test_amtools;
  test_amint;
  test_amconvacc;
  write_sep;

  test_all_complex;
  write_sep;

  test_quat;
  write_sep;

  test_fpu;

end.



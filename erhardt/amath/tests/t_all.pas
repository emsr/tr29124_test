{T_ALL - all-in-one test program for AMATH package  (c) 2011-2013  W.Ehrhardt}

{For 32 bit, code segment too large for 16 bit}

program t_all;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}


{$ifdef BIT16}
  {$ifndef DPMI}
   {$ifndef Windows}
     {$define USEOVR}
   {$endif}
  {$endif}
{$endif}


uses
  {$ifdef USEOVR}
    Overlay,
  {$endif}
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  t_amathm, t_sfx_m, t_sfd_m, t_squadm, t_amtoom, t_acmplu;

{$ifdef USEOVR}

  {$o t_sfx_m}
  {$o t_sfx1}
  {$o t_sfx2}
  {$o t_sfx1a}
  {$o t_sfx3}
  {$o t_sfx3a}
  {$o t_sfx3b}
  {$o t_sfx3c}
  {$o t_sfx4}
  {$o t_sfx4a}
  {$o t_sfx5}
  {$o t_sfx6}
  {$o t_sfx6a}
  {$o t_sfx7}

  {$o t_sfd1}
  {$o t_sfd2}
  {$o t_sfd1a}
  {$o t_sfd3}
  {$o t_sfd3a}
  {$o t_sfd3b}
  {$o t_sfd3c}
  {$o t_sfd4}
  {$o t_sfd4a}
  {$o t_sfd5}
  {$o t_sfd6}
  {$o t_sfd6a}
  {$o t_sfd7}
{$endif}

procedure write_sep;
begin
  writeln;
  writeln('********************************************************');
  writeln('********************************************************');
end;


begin
  writeln('--------------------------------------------------');
  writeln('All-in-one test AMath package  (c) W.Ehrhardt 2012');
  writeln('--------------------------------------------------');
  {$ifdef USEOVR}
    OvrInit('T_ALL.OVR');
    if OvrResult <> ovrOk then begin
      case OvrResult of
        ovrError: Writeln('Program has no overlays.');
        ovrNotFound: Writeln('Overlay file not found.');
      end;
      Halt(1);
    end;
  {$endif}
  test_amath_main;
  write_sep;

  test_sfx_main;
  write_sep;

  test_sfd_main;
  write_sep;

  test_quadsolve;
  test_cubsolve;
  write_sep;

  test_amtools;
  test_amint;
  test_amconvacc;
  write_sep;

  test_all_complex; 

end.

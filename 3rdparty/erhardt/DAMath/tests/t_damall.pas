{T_DAMALL - all-in-one test program for DAMATH package  (c) 2013  W.Ehrhardt}

{For BP7D or better}

program t_damall;

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
  t_damatm, t_drootm, t_damttm, t_dcmplu, t_sfd_m, t_dquatm, t_dfpum;

  {$ifdef USEOVR}
    {$o t_sfd_m}
    {$o t_sfd1}
    {$o t_sfd1a}
    {$o t_sfd1b}
    {$o t_sfd2}
    {$o t_sfd3}
    {$o t_sfd3a}
    {$o t_sfd3b}
    {$o t_sfd3c}
    {$o t_sfd4}
    {$o t_sfd4a}
    {$o t_sfd5}
    {$o t_sfd5a}
    {$o t_sfd6}
    {$o t_sfd6a}
    {$o t_sfd7}
    {$o t_sfd7a}
    {$o t_sfd8}
    {$o t_sfd8a}

    {$o t_damatm}
    {$o t_damat0}
    {$o t_damat1}

    {$o t_drootm}
    {$o t_dquatm}
    {$o t_dfpum}

  {$endif}

{---------------------------------------------------------------------------}
procedure write_sep;
begin
  writeln;
  writeln('********************************************************');
  writeln('********************************************************');
end;


begin

  writeln('---------------------------------------------------');
  writeln('All-in-one test DAMath package  (c) W.Ehrhardt 2018');
  writeln('---------------------------------------------------');

  {$ifdef USEOVR}
    OvrInit('T_DAMALL.OVR');
    if OvrResult <> ovrOk then begin
      case OvrResult of
        ovrError: Writeln('Program has no overlays.');
        ovrNotFound: Writeln('Overlay file not found.');
      end;
      Halt(1);
    end;
  {$endif}

  test_damath_main;
  write_sep;

  test_quadsolve;
  test_cubsolve;
  test_PolyRoots;
  write_sep;

  test_amtools;
  test_amint;
  test_amconvacc;
  write_sep;

  test_sfd_main;
  write_sep;

  test_all_complex;
  write_sep;

  test_quat;
  write_sep;

  test_fpu;

end.

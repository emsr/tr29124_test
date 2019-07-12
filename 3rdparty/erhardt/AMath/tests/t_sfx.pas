{T_SFX - Regression test for SPECFUNX unit  (c) 2010-2018  W.Ehrhardt}

{Many test cases were calculated with Maple V R4, Digits>=20}
{others with Pari/GP 2.3.4 and \p 20, some with Cephes qcalc}

{Note that the 16 bit compilers are not very accurate when compiling}
{20+ digit literals, therefore in most case only 19 digits are used.}

program t_sfx;
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
  {$ifdef debug}
  sfbasic,
  {$endif}
  t_sfx_m;

  {$ifdef USEOVR}
    {$o t_sfx1}
    {$o t_sfx1a}
    {$o t_sfx1b}
    {$o t_sfx2}
    {$o t_sfx2a}
    {$o t_sfx2b}
    {$o t_sfx2c}
    {$o t_sfx3}
    {$o t_sfx3a}
    {$o t_sfx3b}
    {$o t_sfx3c}
    {$o t_sfx3d}
    {$o t_sfx4}
    {$o t_sfx4a}
    {$o t_sfx5}
    {$o t_sfx5a}
    {$o t_sfx6}
    {$o t_sfx6a}
    {$o t_sfx6b}
    {$o t_sfx7}
    {$o t_sfx7a}
    {$o t_sfx8}
    {$o t_sfx8a}
    {$o t_sfx8b}
    {$o t_sfx9}
    {$o t_sfx9a}
  {$endif}

begin
  {$ifdef USEOVR}
    OvrInit('T_SFX.OVR');
    if OvrResult <> ovrOk then begin
      case OvrResult of
        ovrError: Writeln('Program has no overlays.');
        ovrNotFound: Writeln('Overlay file not found.');
      end;
      Halt(1);
    end;
  {$endif}
  test_sfx_main;
  {$ifdef debug}
    sfc_dump_diagctr;
  {$endif}
end.


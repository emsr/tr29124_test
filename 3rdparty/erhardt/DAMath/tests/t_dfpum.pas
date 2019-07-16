{Test unit for DFPU functions  (c) W.Ehrhardt 2017-2018}

unit t_dfpum;

interface

{$i STD.INC}

{$ifdef BIT16}
  {$N+,X+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}


procedure test_fpu;
  {-Test of DFPU functions}

implementation


uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  dfpu,
  damath;


{---------------------------------------------------------------------------}
procedure test_rp_double;
var
  rm: TFPURoundingMode;
  eps,x,y: double;
begin
  SetPrecisionMode(pmDouble);
  writeln;
  writeln('Double precision');
  writeln('----------------');
  for rm := rmNearest to rmTruncate do begin
    SetRoundMode(rm);
    writeln(' Rounding mode ', RoundingNames[GetRoundMode]);
    eps := 1.0;
    repeat
      x := 0.5*eps;
      if rm=rmUp then y := 1-x else y := 1+x;
      if y=1.0 then break;
      eps := x;
    until false;
    writeln('  eps = ', eps);
    eps := 1e-20;
    x := +1 - eps;
    writeln('  +1 - 1e-20 = ', x:25:16, '    ', dbl2hex(x));
    x := +1 + eps;
    writeln('  +1 + 1e-20 = ', x:25:16, '    ', dbl2hex(x));
    x := -1 - eps;
    writeln('  -1 - 1e-20 = ', x:25:16, '    ', dbl2hex(x));
    x := -1 + eps;
    writeln('  -1 + 1e-20 = ', x:25:16, '    ', dbl2hex(x));
    x := ldexpd(1,-1200);
    writeln('  ldexpd(+1,-1200) = ', x:26, '    ', dbl2hex(x));
    x := ldexpd(-1,-1200);
    writeln('  ldexpd(-1,-1200) = ', x:26, '    ', dbl2hex(x));
  end;
end;

{---------------------------------------------------------------------------}
procedure test_rp_single;
var
  rm: TFPURoundingMode;
  eps,x,y: single;
begin
  SetPrecisionMode(pmDouble);
  writeln;
  writeln('Single precision');
  writeln('----------------');
  for rm := rmNearest to rmTruncate do begin
    SetRoundMode(rm);
    writeln(' Rounding mode ', RoundingNames[GetRoundMode]);
    eps := 1.0;
    repeat
      x := 0.5*eps;
      if rm=rmUp then y := 1-x else y := 1+x;
      if y=1.0 then break;
      eps := x;
    until false;
    writeln('  eps = ', eps);
    eps := 1e-20;
    x := +1 - eps;
    writeln('  +1 - 1e-20 = ', x:15:7, '    ', sgl2hex(x));
    x := +1 + eps;
    writeln('  +1 + 1e-20 = ', x:15:7, '    ', sgl2hex(x));
    x := -1 - eps;
    writeln('  -1 - 1e-20 = ', x:15:7, '    ', sgl2hex(x));
    x := -1 + eps;
    writeln('  -1 + 1e-20 = ', x:15:7, '    ', sgl2hex(x));
    x := ldexpd(1,-200);
    writeln('  ldexps(+1,-200) = ', x:20, '    ', sgl2hex(x));
    x := ldexpd(-1,-200);
    writeln('  ldexps(-1,-200) = ', x:20, '    ', sgl2hex(x));
  end;
end;


{---------------------------------------------------------------------------}
procedure test_snan_exeptions;
var
  emask: byte;
  x,y: double;
begin
  writeln;
  writeln('Test signaling NaNs and exceptions');
  writeln('----------------------------------');
  writeln('     sNan_d = ', Dbl2Hex(sNaN_d));
  writeln('     sNan_s = ', Sgl2Hex(sNaN_s));

  GetExceptionMask(emask);
  SetExceptionMask(emask or $3F);

  writeln(' 1 + sNan_d = ', Dbl2Hex(1 + sNaN_d));
  writeln(' 1 + sNan_s = ', Sgl2Hex(1 + sNaN_s));

  x := -Pi;
  y := ln(x);
  writeln('    ln(-Pi) = ', Dbl2Hex(y), '  ',y:10);

  y := sqrt(x);
  writeln('  sqrt(-Pi) = ', Dbl2Hex(y), '  ',y:10);

  x := PosInf_d;
  y := x-x;
  writeln('  Inf - Inf = ', Dbl2Hex(y), '  ',y:10);

  y := NegInf_d/PosInf_d;
  writeln(' -Inf / Inf = ', Dbl2Hex(y), '  ',y:10);

  SetExceptionMask(emask);

end;



{---------------------------------------------------------------------------}
procedure test_fpu;
  {-Test of DFPU functions}
var
  init_rm: TFPURoundingMode;
  init_pm: TFPUPrecisionMode;
begin

  writeln('T_FPU - test program for DFPU functions   (c) 2017-2018 W.Ehrhardt');
  writeln;

  init_rm := GetRoundMode;
  init_pm := GetPrecisionMode;

  writeln('Check rounding to nearest');
  writeln('         Default: ', isRMNearest);

  SetRoundMode(rmTruncate);
  writeln('After rmTruncate: ', isRMNearest);

  SetRoundMode(rmUp);
  writeln('      After rmUp: ', isRMNearest);

  SetRoundMode(rmDown);
  writeln('    After rmDown: ', isRMNearest);

  SetRoundMode(rmNearest);
  writeln(' After rmNearest: ', isRMNearest);
  writeln;

  write('Test pmSingle  : ');
  SetPrecisionMode(pmSingle);
  writeln(PrecisionNames[GetPrecisionMode]);

  write('Test pmReserved: ');
  SetPrecisionMode(pmReserved);
  writeln(PrecisionNames[GetPrecisionMode]);

  write('Test pmDouble  : ');
  SetPrecisionMode(pmDouble);
  writeln(PrecisionNames[GetPrecisionMode]);

  write('Test pmExtended: ');
  SetPrecisionMode(pmExtended);
  writeln(PrecisionNames[GetPrecisionMode]);
  writeln;

  SetPrecisionMode(init_pm);

  test_rp_double;
  test_rp_single;

  SetRoundMode(init_rm);
  SetPrecisionMode(init_pm);

  test_snan_exeptions;

end;

end.


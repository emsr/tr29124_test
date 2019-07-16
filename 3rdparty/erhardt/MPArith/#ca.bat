if not exist m#add.dat goto ende

@echo off

rem generic compile batch file for typical console compilers
rem pause after each compiler if %1 not empty
rem call self on low env space condition
rem (c) 2003-2014 W.Ehrhardt

rem test file
rem =========

set SRC=t_check
set SRC=t_chkrat
::set SRC=t_chkfp
::set SRC=t_chkmpc
::set SRC=t_intro
::set SRC=t_cattle

::copy t_4sq.dpr t_4sq.pas
::set SRC=t_4sq


rem log file (may be con or nul)
rem ============================
::set LOG=nul
set LOG=%SRC%.LOG


rem parameters for test file
rem ========================
set PARA=

rem build options
rem ========================
::set OPT=-dFullData
::set OPT=-ddebug
::set OPT=


rem test whether enough space in environment
rem ========================================
set PCB=A_rather_long_environment_string_for_testing
if  (%PCB%)==(A_rather_long_environment_string_for_testing) goto OK


rem call self with 4096 byte env
rem ============================
set PCB=
%COMSPEC% /E:4096 /C %0 %1 %2 %3 %4 %5 %6 %7 %8 %9
goto ende


:OK

echo Test %SRC% for many compilers  >%LOG%
ver >>%LOG%

set PCB=call p7 -b -dMPC_MAXRadix64
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

set PCB=call p7d -CP -b
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

set PCB=call fpc1 -B -TGO32V2
rem WE: FPC 1.0.10
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

:: set PCB=call fpc2 -B -dMPC_MAXRadix64
:: rem WE: FPC 2.0.4
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

:: set PCB=call fpc22 -B
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

:: set PCB=call fpc222 -B
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

set PCB=call fpc224 -B -dMPC_FPrec30K
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

:: set PCB=call fpc240 -B
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

:: set PCB=call fpc242 -B
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

set PCB=call fpc244 -B
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

:: set PCB=call fpc260 -B
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

:: set PCB=call fpc262 -B
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

set PCB=call fpc264d -B -dMPC_MAXRadix64
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

set PCB=call fpc264 -B
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

:: set PCB=call fpc300 -B
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

:: set PCB=call fpc302 -B -dMPC_MAXRadix64
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

set PCB=call fpc304 -B -dMPC_MAXRadix64
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

set PCB=call fpc311 -B
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

set PCB=call vpc -b -dMPC_MAXRadix64
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

:: set PCB=D:\DMX\M2\DCC32.EXE -b
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

set PCB=D:\DMX\M3\DCC32.EXE -b -dMPC_MAXRadix64
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

:: set PCB=D:\DMX\M4\DCC32.EXE -b -dMPC_FPrec30K
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

:: set PCB=D:\DMX\M5\DCC32.EXE -b  -dMPC_MAXRadix64 -dMPC_FPrec30K
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

set PCB=D:\DMX\M6\DCC32.EXE -b
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

set PCB=D:\DMX\M7\DCC32.EXE -b  -dMPC_MAXRadix64 -dMPC_FPrec30K
del %SRC%.exe >nul
%PCB% %OPT% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

:: set PCB=D:\DMX\M9\DCC32.EXE -b
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

:: set PCB=D:\DMX\M10\DCC32.EXE -b -dMPC_MAXRadix64
:: del %SRC%.exe >nul
:: %PCB% %OPT% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

set PCB=D:\DMX\M12\DCC32.EXE -b
del %SRC%.exe >nul
%PCB% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%

:: set PCB=D:\DMX\M17\DCC32.EXE -b
:: del %SRC%.exe >nul
:: %PCB% %SRC%.pas
:: if not (%1%)==() pause
:: echo. >>%LOG%
:: echo Results for %PCB% >>%LOG%
:: %SRC%.exe %PARA% >>%LOG%

set PCB=D:\DMX\M18\DCC32.EXE -b   -dMPC_MAXRadix64
del %SRC%.exe >nul
%PCB% %SRC%.pas
if not (%1%)==() pause
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
%SRC%.exe %PARA% >>%LOG%


echo.
echo **** Log file: %LOG%

:ende

set PCB=
set SRC=
set LOG=
set PARA=
set OPT=


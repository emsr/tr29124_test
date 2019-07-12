@echo off

set SRC=%1

if not (%SRC%)==() goto foundsrc
echo *************** NO SOURCE FILE ********************
goto end

:foundsrc
set PARA=
set LOG=%SRC%.log
set ERROR=

echo Test %SRC% for (most) console compilers >%LOG%
ver >>%LOG%

set PCB=bpc -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=bpc -CP -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=call p5 -b -l
call $ca1
if not (%ERROR%)==() goto end

set PCB=call p55 -b -l
call $ca1
if not (%ERROR%)==() goto end

set PCB=call p6 -l
call $ca1
if not (%ERROR%)==() goto end

set PCB=call vpc -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc1 -B -TGO32V2
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc2 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc22 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc222 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc224 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc240 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc242 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc244 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc260 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc262 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc264 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc264d -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc300 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc302 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc304 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc311 -B
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M2\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M3\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M4\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M5\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M6\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M7\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M9\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

call wdosx %SRC%.exe
echo. >>%LOG%
echo Results for WDOSX >>%LOG%
%SRC%.exe %PARA% >>%LOG%

set PCB=D:\DMX\M10\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M12\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M17\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

set PCB=D:\DMX\M18\DCC32.EXE -b
call $ca1
if not (%ERROR%)==() goto end

:end

echo.
echo **** Log file: %LOG%

set PCB=
set SRC=
set LOG=
set PARA=
set ERROR=


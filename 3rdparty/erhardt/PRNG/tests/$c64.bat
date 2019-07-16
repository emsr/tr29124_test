@echo off

set SRC=%1

if not (%SRC%)==() goto foundsrc
echo *************** NO SOURCE FILE ********************
goto end

:foundsrc
set PARA=
set LOG=%SRC%.l64
set ERROR=

echo Test %SRC% for (most) console compilers >%LOG%
ver >>%LOG%

set PCB=call fpc64264 -B -O3
call $ca1
if not (%ERROR%)==() goto end

set PCB=call fpc64304 -B -O3
call $ca1
if not (%ERROR%)==() goto end

set PCB=call b1764
call $ca1
if not (%ERROR%)==() goto end

set PCB=call b1864
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


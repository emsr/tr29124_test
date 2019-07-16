@echo off


if (%SRC%)==() goto end
if (%LOG%)==() goto end
::if (%PCB%)==() goto end

echo. >>%LOG%
echo. >>%LOG%
echo Results for %PCB% >>%LOG%
echo --------------------------------- >>%LOG%
del %SRC%.exe >nul
if exist %SRC%.exe goto delerr
%PCB% %SRC%.pas
if errorlevel 1 goto comperr
%SRC%.exe %PARA% >>%LOG%
goto end

:delerr
echo *************  Error deleting %SRC%.exe ************** >>%LOG%
goto end

:delerr
echo *************  Compile error %SRC%.pas ************** >>%LOG%
goto end

:end


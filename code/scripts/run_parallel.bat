@echo off
setlocal EnableExtensions

for %%I in ("%~dp0..\..") do set "REPOSITORY_ROOT=%%~fI"
set "OUTPUT_ROOT=%REPOSITORY_ROOT%\results\parallel-experiments"
set "CHECK_ONLY=0"
set "CHECK_OPTION="

:parse_arguments
if "%~1"=="" goto arguments_done
if /I "%~1"=="--output-dir" goto parse_output
if /I "%~1"=="--check-only" goto parse_check
if /I "%~1"=="--help" goto help
if /I "%~1"=="-h" goto help
echo ERROR: Unknown argument: "%~1"
goto usage_error

:parse_output
if "%~2"=="" goto missing_output
set "OUTPUT_ROOT=%~2"
shift
shift
goto parse_arguments

:parse_check
set "CHECK_ONLY=1"
set "CHECK_OPTION=--check-only"
shift
goto parse_arguments

:arguments_done
for %%I in ("%OUTPUT_ROOT%") do set "OUTPUT_ROOT=%%~fI"
set "PRECPACK_REQUIRE_GUROBI_RUNTIME=1"
if "%CHECK_ONLY%"=="1" goto gurobi_check_done
call "%~dp0run_batch.bat" --help | findstr /C:"Optional Gurobi root strengthening: enabled" >nul
if errorlevel 1 goto missing_gurobi

:gurobi_check_done

set "ITEM_DIRECTORY=%REPOSITORY_ROOT%\data\instances\otto\n_0100"
set "GRAPH_01=%REPOSITORY_ROOT%\data\bpp-gp-graphs\separation-01\n_0100"
set "GRAPH_03=%REPOSITORY_ROOT%\data\bpp-gp-graphs\separation-03\n_0100"
set /a ITEM_COUNT=0
set /a GRAPH_01_COUNT=0
set /a GRAPH_03_COUNT=0
for %%F in ("%ITEM_DIRECTORY%\*.txt") do if exist "%%~fF" set /a ITEM_COUNT+=1 >nul
for %%F in ("%GRAPH_01%\*.graph") do if exist "%%~fF" set /a GRAPH_01_COUNT+=1 >nul
for %%F in ("%GRAPH_03%\*.graph") do if exist "%%~fF" set /a GRAPH_03_COUNT+=1 >nul
if not "%ITEM_COUNT%"=="525" goto invalid_matrix
if not "%GRAPH_01_COUNT%"=="525" goto invalid_matrix
if not "%GRAPH_03_COUNT%"=="525" goto invalid_matrix

echo PrecPack controlled parallel experiments
echo   datasets       : 4
echo   threads        : 2, 4, 8
echo   configurations : 12
echo   total runs     : 6300
echo   output         : "%OUTPUT_ROOT%"

:run_matrix
for %%T in (2 4 8) do (
    call :run_configuration salbp-i salbp-i 350 %%T ""
    if errorlevel 1 goto run_failed
    call :run_configuration bpp-p bpp-p 1000 %%T ""
    if errorlevel 1 goto run_failed
    call :run_configuration bpp-gp-separation-01 bpp-gp 75 %%T "%GRAPH_01%"
    if errorlevel 1 goto run_failed
    call :run_configuration bpp-gp-separation-03 bpp-gp 75 %%T "%GRAPH_03%"
    if errorlevel 1 goto run_failed
)
if "%CHECK_ONLY%"=="1" echo Parallel experiment validation completed.
if "%CHECK_ONLY%"=="0" echo Parallel experiment runs completed.
goto success

:run_configuration
echo.
echo [%~1, threads=%~4]
if "%~5"=="" goto run_without_graph
call "%~dp0run_batch.bat" --problem %~2 --input "%ITEM_DIRECTORY%" --graph-dir "%~5" --time-limit %~3 --memory-limit-mb 24576 --threads %~4 --output-dir "%OUTPUT_ROOT%\%~1\threads-%~4" %CHECK_OPTION%
exit /b %ERRORLEVEL%

:run_without_graph
call "%~dp0run_batch.bat" --problem %~2 --input "%ITEM_DIRECTORY%" --time-limit %~3 --memory-limit-mb 24576 --threads %~4 --output-dir "%OUTPUT_ROOT%\%~1\threads-%~4" %CHECK_OPTION%
exit /b %ERRORLEVEL%

:missing_output
echo ERROR: --output-dir requires a value.
goto usage_error

:missing_gurobi
echo ERROR: The controlled experiment requires a Gurobi-enabled build.
echo Rebuild with code\scripts\build.bat --gurobi on.
exit /b 1

:invalid_matrix
echo ERROR: Each n=100 experiment data set must contain 525 files.
exit /b 1

:run_failed
echo ERROR: A parallel experiment configuration failed.
exit /b 1

:success
exit /b 0

:usage_error
echo Usage: %~nx0 [--output-dir DIR] [--check-only]
exit /b 2

:help
echo Usage: %~nx0 [--output-dir DIR] [--check-only]
exit /b 0

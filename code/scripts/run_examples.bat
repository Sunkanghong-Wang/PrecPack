@echo off
setlocal EnableExtensions

for %%I in ("%~dp0..\..") do set "REPOSITORY_ROOT=%%~fI"
set "THREADS=1"
set "OUTPUT_ROOT=%REPOSITORY_ROOT%\results\examples"

:parse_arguments
if "%~1"=="" goto arguments_done
if /I "%~1"=="--threads" goto parse_threads
if /I "%~1"=="--output-dir" goto parse_output
if /I "%~1"=="--help" goto help
if /I "%~1"=="-h" goto help
echo ERROR: Unknown argument: "%~1"
goto usage_error

:parse_threads
if "%~2"=="" goto missing_threads
set "THREADS=%~2"
shift
shift
goto parse_arguments

:parse_output
if "%~2"=="" goto missing_output
set "OUTPUT_ROOT=%~2"
shift
shift
goto parse_arguments

:arguments_done
for %%I in ("%OUTPUT_ROOT%") do set "OUTPUT_ROOT=%%~fI"
set "BASE_INSTANCE=%REPOSITORY_ROOT%\data\instances\otto\n_0020\instance_n=20_1.txt"

call "%~dp0run_batch.bat" --problem bpp-p --input "%BASE_INSTANCE%" --time-limit 10 --threads "%THREADS%" --output-dir "%OUTPUT_ROOT%\bpp-p"
if errorlevel 1 exit /b %ERRORLEVEL%
call "%~dp0run_batch.bat" --problem salbp-i --input "%REPOSITORY_ROOT%\data\instances\scholl\Bowman\Bowman_c20.txt" --time-limit 10 --threads "%THREADS%" --output-dir "%OUTPUT_ROOT%\salbp-i"
if errorlevel 1 exit /b %ERRORLEVEL%
call "%~dp0run_batch.bat" --problem bpp-gp --input "%BASE_INSTANCE%" --graph-dir "%REPOSITORY_ROOT%\data\bpp-gp-graphs\separation-03\n_0020\instance_n=20_1.graph" --time-limit 10 --threads "%THREADS%" --output-dir "%OUTPUT_ROOT%\bpp-gp"
exit /b %ERRORLEVEL%

:missing_threads
echo ERROR: --threads requires a value.
goto usage_error

:missing_output
echo ERROR: --output-dir requires a value.
goto usage_error

:usage_error
echo Usage: %~nx0 [--threads N] [--output-dir DIR]
exit /b 2

:help
echo Usage: %~nx0 [--threads N] [--output-dir DIR]
exit /b 0

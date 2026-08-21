@echo off
setlocal EnableExtensions

for %%I in ("%~dp0..\..") do set "REPOSITORY_ROOT=%%~fI"
set "BINARY=%REPOSITORY_ROOT%\build\Release\precpack.exe"
if exist "%BINARY%" goto binary_ready
set "BINARY=%REPOSITORY_ROOT%\build\precpack.exe"
if exist "%BINARY%" goto binary_ready
echo ERROR: PrecPack is not built. Run code\scripts\build.bat first.
exit /b 1

:binary_ready
set "OMP_NUM_THREADS=1"
set "MKL_NUM_THREADS=1"
set "OPENBLAS_NUM_THREADS=1"
set "BLIS_NUM_THREADS=1"
set "VECLIB_MAXIMUM_THREADS=1"
if defined GUROBI_HOME set "PATH=%GUROBI_HOME%\bin;%PATH%"
set "PRECPACK_REPOSITORY_ROOT=%REPOSITORY_ROOT%"

"%BINARY%" --batch %*
exit /b %ERRORLEVEL%

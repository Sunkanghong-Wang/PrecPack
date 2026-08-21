@echo off
setlocal EnableExtensions

for %%I in ("%~dp0..\..") do set "REPOSITORY_ROOT=%%~fI"
set "CODE_DIRECTORY=%REPOSITORY_ROOT%\code"
set "BUILD_DIRECTORY=%REPOSITORY_ROOT%\build"
set "GUROBI_MODE=AUTO"

:parse_arguments
if "%~1"=="" goto arguments_done
if /I "%~1"=="--help" goto help
if /I "%~1"=="-h" goto help
if /I "%~1"=="--gurobi" goto parse_gurobi
echo ERROR: Unknown argument: "%~1"
goto usage_error

:parse_gurobi
if "%~2"=="" goto missing_gurobi_value
set "GUROBI_MODE=%~2"
shift
shift
goto parse_arguments

:arguments_done
if /I "%GUROBI_MODE%"=="auto" set "GUROBI_MODE=AUTO"
if /I "%GUROBI_MODE%"=="on" set "GUROBI_MODE=ON"
if /I "%GUROBI_MODE%"=="off" set "GUROBI_MODE=OFF"
if not "%GUROBI_MODE%"=="AUTO" if not "%GUROBI_MODE%"=="ON" if not "%GUROBI_MODE%"=="OFF" goto invalid_gurobi_value

echo PrecPack Release build
echo   Source       : "%CODE_DIRECTORY%"
echo   Build        : "%BUILD_DIRECTORY%"
echo   Gurobi mode  : %GUROBI_MODE%
echo.

if defined CMAKE_BIN goto configured_cmake
where cmake >nul 2>nul
if errorlevel 1 goto cmake_missing
set "CMAKE_COMMAND=cmake"
goto cmake_ready

:configured_cmake
set "CMAKE_COMMAND=%CMAKE_BIN%"

:cmake_ready
echo   CMake        : "%CMAKE_COMMAND%"

set "OMP_NUM_THREADS=1"
set "MKL_NUM_THREADS=1"
set "OPENBLAS_NUM_THREADS=1"
set "BLIS_NUM_THREADS=1"
set "VECLIB_MAXIMUM_THREADS=1"
if defined GUROBI_HOME set "PATH=%GUROBI_HOME%\bin;%PATH%"

echo.
echo [1/3] Configuring PrecPack...
"%CMAKE_COMMAND%" -S "%CODE_DIRECTORY%" -B "%BUILD_DIRECTORY%" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DPRECPACK_NATIVE_OPTIMIZATION=ON -DPRECPACK_ENABLE_IPO=ON -DPRECPACK_ENABLE_SANITIZERS=OFF -DPRECPACK_GUROBI=%GUROBI_MODE% "-DGUROBI_ROOT=%GUROBI_HOME%"
if errorlevel 1 goto configure_failed

echo.
echo [2/3] Building PrecPack...
"%CMAKE_COMMAND%" --build "%BUILD_DIRECTORY%" --config Release --target precpack --parallel
if errorlevel 1 goto build_failed

echo.
echo [3/3] Running the test suite...
"%CMAKE_COMMAND%" --build "%BUILD_DIRECTORY%" --config Release --target check
if errorlevel 1 goto test_failed

echo.
echo PrecPack was built and tested successfully.
if exist "%BUILD_DIRECTORY%\Release\precpack.exe" echo   Executable: "%BUILD_DIRECTORY%\Release\precpack.exe"
if exist "%BUILD_DIRECTORY%\precpack.exe" echo   Executable: "%BUILD_DIRECTORY%\precpack.exe"
exit /b 0

:cmake_missing
echo ERROR: CMake was not found.
echo Install CMake 3.20 or newer and make cmake available on PATH,
echo or set CMAKE_BIN to the full path of cmake.exe.
exit /b 1

:missing_gurobi_value
echo ERROR: --gurobi requires one of: auto, on, off.
goto usage_error

:invalid_gurobi_value
echo ERROR: Invalid --gurobi value: "%GUROBI_MODE%"
goto usage_error

:configure_failed
echo.
echo ERROR: PrecPack configuration failed. Review the CMake messages above.
exit /b 1

:build_failed
echo.
echo ERROR: PrecPack compilation failed. Review the compiler messages above.
exit /b 1

:test_failed
echo.
echo ERROR: PrecPack tests failed. Review the test output above.
exit /b 1

:usage_error
echo Usage: %~nx0 [--gurobi auto^|on^|off]
exit /b 2

:help
echo Usage: %~nx0 [--gurobi auto^|on^|off]
echo.
echo Configure, build, and test PrecPack in Release mode.
echo This Windows launcher calls CMake directly and does not require Python.
exit /b 0

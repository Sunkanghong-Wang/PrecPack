@echo off
setlocal

where py >nul 2>nul
if errorlevel 1 goto python
py -3 "%~dpn0.py" %*
exit /b %errorlevel%

:python
where python >nul 2>nul
if errorlevel 1 goto missing
python "%~dpn0.py" %*
exit /b %errorlevel%

:missing
echo error: Python 3 was not found. 1>&2
exit /b 1

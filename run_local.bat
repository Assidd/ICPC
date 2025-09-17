@echo off
setlocal EnableExtensions EnableDelayedExpansion

set "IN=%~1"
if "%IN%"=="" set "IN=in.txt"
if not exist "%IN%" (
  echo Input file "%IN%" not found. Put a test into "in.txt" or pass a path as arg.
  exit /b 1
)

echo [1/5] Build solver
g++ -O3 -pipe -s -std=c++17 -static main.cpp -o solver.exe 2>nul
if errorlevel 1 g++ -O3 -pipe -s -std=c++17 main.cpp -o solver.exe || exit /b 1

echo [2/5] Run solver
type "%IN%" | solver.exe > out.txt || exit /b 1

echo [3/5] Build checker
g++ -O2 -pipe -s -static helps\checker.cpp -o checker.exe 2>nul
if errorlevel 1 g++ -O2 -pipe -s helps\checker.cpp -o checker.exe || exit /b 1

echo [4/5] Build scorer
g++ -O2 -pipe -s -static helps\scorer.cpp -o scorer.exe 2>nul
if errorlevel 1 g++ -O2 -pipe -s helps\scorer.cpp -o scorer.exe || exit /b 1

echo [5/5] Validate ^& score
checker.exe "%IN%" out.txt || exit /b 1
scorer.exe  "%IN%" out.txt

echo Done.
exit /b 0

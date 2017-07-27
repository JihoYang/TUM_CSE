rem
rem This script creates the mysolver.def...
rem
rem
cl /c /Ox /MD /Fo"tyr369.obj" "tyr369.c"
rem
set OPATH=%PATH%
PATH=C:\PROGRA~1\MSC~1.SOF\MD_Adams\R3\utils\win32;%OPATH%
echo EXPORTS > mysolver.def
%PYTHON_EXE% %SHORT_TOPDIR%utils\bld_gendef.py -e -o tmp.def "tyr369.obj"
type tmp.def >> mysolver.def
del tmp.def
set PATH=%OPATH%
set OPATH=

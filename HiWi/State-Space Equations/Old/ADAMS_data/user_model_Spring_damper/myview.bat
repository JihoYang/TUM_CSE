rem
rem This script creates the myview.def...
rem
rem
cl /c /Ox /MD /Fo"tyr369.obj" "tyr369.c"
cl /c /Ox /MD /Fo"rpf369.obj" "rpf369.c"
rem
set OPATH=%PATH%
PATH=C:\PROGRA~1\MSC~1.SOF\MD_Adams\R3\utils\win32;%OPATH%
echo EXPORTS > myview.def
%PYTHON_EXE% %SHORT_TOPDIR%utils\bld_gendef.py -e -o tmp.def "tyr369.obj"
type tmp.def >> myview.def
del tmp.def
%PYTHON_EXE% %SHORT_TOPDIR%utils\bld_gendef.py -e -o tmp.def "rpf369.obj"
type tmp.def >> myview.def
del tmp.def
set PATH=%OPATH%
set OPATH=

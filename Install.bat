@echo off
title Installing JRI-NativeLibrary
echo Please Type in the path to your R\bin\x64 directory. It should be something like:
echo "C:\Users\Me\Documents\R\R-3.4.2\bin\x64"
set /p rPath=
set rInstallPath=%rPath%\R.exe
%rInstallPath% CMD BATCH %~dp0inst.R
pause
set rPath=%rPath%\*.dll
echo Please Type in the path to your jri\bin\x64 directory. It should be something like:
echo "C:\Users\Me\Documents\R\R-3.4.2\library\rJava\jri\x64"
set /p jriPath=
set jriPath=%jriPath%\*.dll
echo Please Type in the path to your jre\bin directory. It should be something like:
echo "C:\Program Files\Java\jre1.8.0_121\bin"
set /p jrePath=
set jrePath=%jrePath%\.
echo Transferring dlls to JRE home directory . . .
xcopy %rPath% %jrePath% /w
xcopy %jriPath% %jrePath% 

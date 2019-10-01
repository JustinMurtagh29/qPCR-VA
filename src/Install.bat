@echo off
echo Installing R-packages . . .
R CMD BATCH Inst.R

title Installing JRI-NativeLibrary
for /r "C:\" %%i in (R.*) do (
	if "%%~xi" == ".exe" (
		set "var=%%i"
		set "vari=%%~pi"
	)
)
set rPath=C:%vari%
set rInstallPath=%var%

%rInstallPath% CMD BATCH %~dp0inst.R
set rPath="%rPath%*.dll"

for /r "C:\" %%i in (jri.*) do (
	if "%%~xi" == ".dll" (
		set "jvari=%%~pi"
	)
)
set jriPath=C:%jvari%
set jriPath="%jriPath%*.dll"

for /r "C:\Program Files" %%i in (java.*) do (
	if "%%~xi" == ".exe" (
		set "javari=%%~pi"
	)
)
set jrePath=C:%javari%
set jrePath="%jrePath%."
echo Transferring dlls to JRE home directory . . .

xcopy %rPath% %jrePath% /w
xcopy %jriPath% %jrePath% 
pause

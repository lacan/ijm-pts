@ECHO OFF

ECHO Copying BIOP_Pizza.ijm to current folder
copy C:\Fiji\plugins\ActionBar\BIOP_Pizza.ijm D:\Projects\Fiji\biop-pizza-target

ECHO Creating JAR FILE
jar cf BIOP_PTS.jar plugins.config *.ijm

ECHO Copying BIOP_PTS.jar to Fiji folder
copy BIOP_PTS.jar C:\Fiji\plugins\BIOP

PAUSE
# Microsoft Developer Studio Project File - Name="JMLib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=JMLib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "JMLib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "JMLib.mak" CFG="JMLib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "JMLib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "JMLib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "JMLib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x412 /d "NDEBUG"
# ADD RSC /l 0x412 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "JMLib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x412 /d "_DEBUG"
# ADD RSC /l 0x412 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "JMLib - Win32 Release"
# Name "JMLib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\annexb.c
# End Source File
# Begin Source File

SOURCE=.\biaridecod.c
# End Source File
# Begin Source File

SOURCE=.\block.c
# End Source File
# Begin Source File

SOURCE=.\cabac.c
# End Source File
# Begin Source File

SOURCE=.\context_ini.c
# End Source File
# Begin Source File

SOURCE=.\erc_api.c
# End Source File
# Begin Source File

SOURCE=.\erc_do_i.c
# End Source File
# Begin Source File

SOURCE=.\erc_do_p.c
# End Source File
# Begin Source File

SOURCE=.\errorconcealment.c
# End Source File
# Begin Source File

SOURCE=.\filehandle.c
# End Source File
# Begin Source File

SOURCE=.\fmo.c
# End Source File
# Begin Source File

SOURCE=.\header.c
# End Source File
# Begin Source File

SOURCE=.\image.c
# End Source File
# Begin Source File

SOURCE=.\ldecod.c
# End Source File
# Begin Source File

SOURCE=.\leaky_bucket.c
# End Source File
# Begin Source File

SOURCE=.\loopFilter.c
# End Source File
# Begin Source File

SOURCE=.\macroblock.c
# End Source File
# Begin Source File

SOURCE=.\mb_access.c
# End Source File
# Begin Source File

SOURCE=.\mbuffer.c
# End Source File
# Begin Source File

SOURCE=.\memalloc.c
# End Source File
# Begin Source File

SOURCE=.\nal.c
# End Source File
# Begin Source File

SOURCE=.\nal_part.c
# End Source File
# Begin Source File

SOURCE=.\nalu.c
# End Source File
# Begin Source File

SOURCE=.\nalucommon.c
# End Source File
# Begin Source File

SOURCE=.\output.c
# End Source File
# Begin Source File

SOURCE=.\parset.c
# End Source File
# Begin Source File

SOURCE=.\parsetcommon.c
# End Source File
# Begin Source File

SOURCE=.\sei.c
# End Source File
# Begin Source File

SOURCE=.\transform8x8.c
# End Source File
# Begin Source File

SOURCE=.\vlc.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\annexb.h
# End Source File
# Begin Source File

SOURCE=.\biaridecod.h
# End Source File
# Begin Source File

SOURCE=.\block.h
# End Source File
# Begin Source File

SOURCE=.\cabac.h
# End Source File
# Begin Source File

SOURCE=.\context_ini.h
# End Source File
# Begin Source File

SOURCE=.\contributors.h
# End Source File
# Begin Source File

SOURCE=.\ctx_tables.h
# End Source File
# Begin Source File

SOURCE=.\defines.h
# End Source File
# Begin Source File

SOURCE=.\elements.h
# End Source File
# Begin Source File

SOURCE=.\erc_api.h
# End Source File
# Begin Source File

SOURCE=.\erc_do.h
# End Source File
# Begin Source File

SOURCE=.\erc_globals.h
# End Source File
# Begin Source File

SOURCE=.\errorconcealment.h
# End Source File
# Begin Source File

SOURCE=.\fmo.h
# End Source File
# Begin Source File

SOURCE=.\global.h
# End Source File
# Begin Source File

SOURCE=.\header.h
# End Source File
# Begin Source File

SOURCE=.\image.h
# End Source File
# Begin Source File

SOURCE=.\ldecod.h
# End Source File
# Begin Source File

SOURCE=.\leaky_bucket.h
# End Source File
# Begin Source File

SOURCE=.\loopfilter.h
# End Source File
# Begin Source File

SOURCE=.\macroblock.h
# End Source File
# Begin Source File

SOURCE=.\mb_access.h
# End Source File
# Begin Source File

SOURCE=.\mbuffer.h
# End Source File
# Begin Source File

SOURCE=.\memalloc.h
# End Source File
# Begin Source File

SOURCE=.\nalu.h
# End Source File
# Begin Source File

SOURCE=.\nalucommon.h
# End Source File
# Begin Source File

SOURCE=.\output.h
# End Source File
# Begin Source File

SOURCE=.\parset.h
# End Source File
# Begin Source File

SOURCE=.\parsetcommon.h
# End Source File
# Begin Source File

SOURCE=.\sei.h
# End Source File
# Begin Source File

SOURCE=.\transform8x8.h
# End Source File
# Begin Source File

SOURCE=.\vlc.h
# End Source File
# End Group
# Begin Source File

SOURCE=.\decoder.cfg
# End Source File
# End Target
# End Project

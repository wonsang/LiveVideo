// HMCPT.h : main header file for the HMCPT application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"       // main symbols

// view mode
#define ANALYZER_MODE		4010

// CHMCPTApp:
// See HMCPT.cpp for the implementation of this class
//

class CHMCPTApp : public CWinApp
{
public:
	CHMCPTApp();


// Overrides
public:
	virtual BOOL InitInstance();

// Implementation

public:
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CHMCPTApp theApp;
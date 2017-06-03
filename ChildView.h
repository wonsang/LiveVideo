// ChildView.h : interface of the CChildView class
//


#pragma once
#include "maindlg.h"
#include "bkdextdlg.h"

// CChildView window

class CChildView : public CWnd
{
// Construction
public:
	CChildView();

// Attributes
public:
	CMainDlg mDlg;

// Operations
public:

// Overrides
	protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

// Implementation
public:
	virtual ~CChildView();

	// Generated message map functions
protected:
	afx_msg void OnPaint();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnToolsAnalyzer();
	
//	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnToolsBkdExtractor();
	CBkdextDlg bDlg;
};


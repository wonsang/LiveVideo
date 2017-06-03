// ChildView.cpp : implementation of the CChildView class
//

#include "stdafx.h"
#include "HMCPT.h"
#include "ChildView.h"

//custom
#include "MainFrm.h"
#include "MainDlg.h"
#include "BkdextDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CChildView

CChildView::CChildView()
{
}

CChildView::~CChildView()
{
}


BEGIN_MESSAGE_MAP(CChildView, CWnd)
	ON_WM_PAINT()
	ON_COMMAND(ID_TOOLS_ANALYZER, &CChildView::OnToolsAnalyzer)
//	ON_WM_CREATE()
ON_COMMAND(ID_TOOLS_BKDEXTRACTOR, &CChildView::OnToolsBkdExtractor)
END_MESSAGE_MAP()



// CChildView message handlers

BOOL CChildView::PreCreateWindow(CREATESTRUCT& cs) 
{
	if (!CWnd::PreCreateWindow(cs))
		return FALSE;

	cs.dwExStyle |= WS_EX_CLIENTEDGE;
	cs.style &= ~WS_BORDER;
	cs.lpszClass = AfxRegisterWndClass(CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, 
		::LoadCursor(NULL, IDC_ARROW), reinterpret_cast<HBRUSH>(COLOR_WINDOW+1), NULL);

	return TRUE;
}

void CChildView::OnPaint() 
{
	CPaintDC dc(this); // device context for painting
	
	// TODO: Add your message handler code here
	
	// Do not call CWnd::OnPaint() for painting messages
}

void CChildView::OnToolsAnalyzer()
{
	// TODO: 여기에 명령 처리기 코드를 추가합니다.
	CMainFrame *pFrame = (CMainFrame *)AfxGetMainWnd();
	mDlg.DoModal();
}
void CChildView::OnToolsBkdExtractor()
{
	// TODO: 여기에 명령 처리기 코드를 추가합니다.
	CMainFrame *pFrame = (CMainFrame *)AfxGetMainWnd();
	bDlg.DoModal();
}

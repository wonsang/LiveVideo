#pragma once
#include "afxwin.h"

// Custom
#include "Convert.h"

#define clip(a)		((a)<0?0:((a)>255?255:(a)))

// CBkdextDlg 대화 상자입니다.

class CBkdextDlg : public CDialog
{
	DECLARE_DYNAMIC(CBkdextDlg)

public:
	CBkdextDlg(CWnd* pParent = NULL);   // 표준 생성자입니다.
	virtual ~CBkdextDlg();

// 대화 상자 데이터입니다.
	enum { IDD = IDD_DIALOG_BKDEXT };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 지원입니다.

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonOpen();
	afx_msg void OnBnClickedButtonBkd();
	CString PathName;
	int width;
	int height;
	// Picture Control
	CStatic Pic;
	afx_msg void OnPaint();
	LPBITMAPINFO BmpInfo;
	HANDLE hloc;
	LPBYTE RGBbuf;
	int ShowFlag;
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	ColorSpaceConversions pCon;
};

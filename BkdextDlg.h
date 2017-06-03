#pragma once
#include "afxwin.h"

// Custom
#include "Convert.h"

#define clip(a)		((a)<0?0:((a)>255?255:(a)))

// CBkdextDlg ��ȭ �����Դϴ�.

class CBkdextDlg : public CDialog
{
	DECLARE_DYNAMIC(CBkdextDlg)

public:
	CBkdextDlg(CWnd* pParent = NULL);   // ǥ�� �������Դϴ�.
	virtual ~CBkdextDlg();

// ��ȭ ���� �������Դϴ�.
	enum { IDD = IDD_DIALOG_BKDEXT };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV �����Դϴ�.

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

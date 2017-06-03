#pragma once
#include "afxwin.h"
#include "afxcmn.h"

// Custom
#include "Convert.h"

#define clip(a)					((a)<0?0:((a)>255?255:(a)))
#define IsNearlyZero(a,tol)		(((a)>=0)&&((a)<=tol))||(((a)<0)&&((a)>=-tol))

// CMainDlg 대화 상자입니다.

class CMainDlg : public CDialog
{
	DECLARE_DYNAMIC(CMainDlg)

public:
	CMainDlg(CWnd* pParent = NULL);   // 표준 생성자입니다.
	virtual ~CMainDlg();

// 대화 상자 데이터입니다.
	enum { IDD = IDD_DIALOG_MAIN };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 지원입니다.

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonAnal();
	// The main function of decoding H.264/AVC videos
	void MainDec(CString filename);
	afx_msg void OnBnClickedButtonOpen();
	// Input path Name
	CString PathName;
	// image width
	int width;
	// image height
	int height;
	// Max frame number
	int MaxFrmNum;
	// Allocate buffers
	void AllocBuf(void);
	// Free buffers
	void FreeBuf(void);
	// Write an output video file
	void WriteDec(void);
	void write_yuv_file(unsigned char *srcY, unsigned char *srcU, unsigned char *srcV, FILE* Write_FILE);
	void arrcpy(unsigned char* DupArr, unsigned char* OrgArr, int size);
	CButton m_OpenBtn;
	CButton m_AnalBtn;
	afx_msg void OnBnClickedButtonBkd();
	afx_msg void OnBnClickedButtonGo();
	void WriteRes(void);
	void WriteBks(void);
	CString BkdPathName;
	afx_msg void OnBnClickedButtonNext();
	afx_msg void OnBnClickedButtonPrev();
	afx_msg void OnPaint();
	void DisplayPic(int DisplayPicNum);
	LPBYTE* RGBbuf;
	int ShowFlag;
	LPBITMAPINFO BmpInfo;
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	HANDLE hloc;
	unsigned char** srcY;
	unsigned char** srcU;
	unsigned char** srcV;
	ColorSpaceConversions pCon;
	// Original picture
	CStatic OriPic;
	// Background Subtracted Picture
	CStatic BksPic;
	// Residual Picture
	CStatic ResPic;
	// Original Motion Vectors
	CStatic OmvPic;
	// Filtered Motion Vectors
	CStatic FmvPic;
	int DisplayFrmNum;
	FILE* DECFILE;
	FILE* BKSFILE;
	FILE* RESFILE;
	FILE* OMVFILE;
	FILE* MBTFILE;
	afx_msg void OnBnClickedButtonShow();
	afx_msg void OnStnClickedStaticOmv();
};

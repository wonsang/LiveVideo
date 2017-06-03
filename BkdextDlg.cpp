// BkdextDlg.cpp : 구현 파일입니다.
//

#include "stdafx.h"
#include "HMCPT.h"
#include "BkdextDlg.h"

// CBkdextDlg 대화 상자입니다.

IMPLEMENT_DYNAMIC(CBkdextDlg, CDialog)

CBkdextDlg::CBkdextDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CBkdextDlg::IDD, pParent)
	, PathName(_T(""))
	, width(0)
	, height(0)
	, RGBbuf(0)
	, ShowFlag(0)
{
	width = 352;
	height = 240;
	ShowFlag = 0;

	RGBbuf = (unsigned char*)malloc(width*height*3);
}

CBkdextDlg::~CBkdextDlg()
{
	free(RGBbuf);
}

void CBkdextDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_STATIC_PIC, Pic);
}


BEGIN_MESSAGE_MAP(CBkdextDlg, CDialog)
	ON_BN_CLICKED(IDC_BUTTON_OPEN, &CBkdextDlg::OnBnClickedButtonOpen)
	ON_BN_CLICKED(IDC_BUTTON_BKD, &CBkdextDlg::OnBnClickedButtonBkd)
	ON_WM_PAINT()
	ON_WM_CREATE()
END_MESSAGE_MAP()


// CBkdextDlg 메시지 처리기입니다.

void CBkdextDlg::OnBnClickedButtonOpen()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	char szFilter[] = "YUV video format (*.yuv)|*.yuv|H.264/AVC video format (*.264)|*.264|all files (*.*)|*.*|";
	CFileDialog dlg(TRUE, NULL, NULL, OFN_FILEMUSTEXIST|OFN_PATHMUSTEXIST, szFilter);
		
	if(dlg.DoModal() == IDOK)
	{
		PathName = dlg.GetPathName();
	}
}

void CBkdextDlg::OnBnClickedButtonBkd()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	char szFilter[] = "YUV format (*.cif)|*.cif|YUV format (*.yuv)|*.yuv|All Files (*.*)|*.*||";
	CFileDialog dlg(TRUE, NULL, NULL, OFN_HIDEREADONLY, szFilter);

	UpdateData(true);

	if(dlg.DoModal() == IDOK)
	{
		
		// Save the image
		int FrameSize = width * height * 1.5;

		unsigned char* FrameBuf;
		FrameBuf = new unsigned char [FrameSize];		

		FILE *Video = fopen(PathName, "rb");
		fread(FrameBuf, sizeof(unsigned char), FrameSize, Video);

 		CString FileName = dlg.GetFileName();
		FILE *FileWrite = fopen(FileName, "wb");
		fwrite(FrameBuf, sizeof(unsigned char), FrameSize, FileWrite);
		fclose(FileWrite);		

		// Copy the image to buffer for showing
		unsigned char* srcY;
		unsigned char* srcU;
		unsigned char* srcV;
		srcY = (unsigned char*)malloc(FrameSize);
		srcU = (unsigned char*)malloc(FrameSize/4);
		srcV = (unsigned char*)malloc(FrameSize/4);

		int CbOffSet = height * width;
		int CrOffSet = CbOffSet + (CbOffSet / 4);
		for(int y=0; y < height; y++)
		{
			for(int x=0; x<width; x++)
			{
				srcY[y*width + x] = FrameBuf[y*width + x];
				srcU[(y>>1)*(width>>1) + (x>>1)] = FrameBuf[CbOffSet + (y>>1)*(width>>1) + (x>>1)];
				srcV[(y>>1)*(width>>1) + (x>>1)] = FrameBuf[CrOffSet + (y>>1)*(width>>1) + (x>>1)];
			}
		}

		pCon.YV12_to_RGB24(srcY,srcU,srcV,RGBbuf,width,height);

		free(srcY);
		free(srcU);
		free(srcV);
		free(FrameBuf);		

		ShowFlag = 1;
	}

	UpdateData(false);
}

void CBkdextDlg::OnPaint()
{
	CPaintDC dc(this); // device context for painting
	// TODO: 여기에 메시지 처리기 코드를 추가합니다.
	// 그리기 메시지에 대해서는 CDialog::OnPaint()을(를) 호출하지 마십시오.

	CDC *pDC = Pic.GetDC();
	CRect pRect;
	Pic.GetWindowRect(pRect);

	if(ShowFlag==1)
	{
		BmpInfo->bmiHeader.biBitCount = 24;
		pDC->SetStretchBltMode(STRETCH_DELETESCANS);
		StretchDIBits(	pDC->m_hDC,
						0,0,pRect.Width(),pRect.Height(),				// destination;
						0,0,width,height,				// source;
						RGBbuf, BmpInfo,
						DIB_RGB_COLORS, SRCCOPY);
	}
}

int CBkdextDlg::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CDialog::OnCreate(lpCreateStruct) == -1)
		return -1;

	// TODO:  여기에 특수화된 작성 코드를 추가합니다.
	// allocates memory block into device-independent bitmap structure
	hloc = GlobalAlloc(GMEM_ZEROINIT | GMEM_MOVEABLE,
   		sizeof(BITMAPINFOHEADER) + (sizeof(RGBQUAD) * 256));
	BmpInfo = (LPBITMAPINFO) GlobalLock(hloc);

	// initializes information about the dimensions and color format of a DIB
	BmpInfo->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	BmpInfo->bmiHeader.biPlanes = 1;
	BmpInfo->bmiHeader.biBitCount = 24;
	BmpInfo->bmiHeader.biCompression = BI_RGB;
	BmpInfo->bmiHeader.biWidth = width;
	BmpInfo->bmiHeader.biHeight = height;

	return 0;
}


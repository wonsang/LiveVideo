// MainDlg.cpp : 구현 파일입니다.
//

#include "stdafx.h"
#include "HMCPT.h"
#include "MainDlg.h"

#define NO				0
#define YES				1

// Additional include files
extern "C" {
#include "lib/JM/ldecod.h"
//#include "lib/JM/image.h"
}


// CMainDlg 대화 상자입니다.

IMPLEMENT_DYNAMIC(CMainDlg, CDialog)

CMainDlg::CMainDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CMainDlg::IDD, pParent)
	, PathName(_T(""))
	, width(0)
	, height(0)
	, MaxFrmNum(0)
	, BkdPathName(_T(""))
	, RGBbuf(NULL)
	, ShowFlag(0)
	, DisplayFrmNum(0)
	, DECFILE(NULL)
{
	width = 320;
	height = 240;
	MaxFrmNum = 200;
	ShowFlag = 0;

	int YSize = width*height;
	int UVSize = YSize/4;
	RGBbuf = (unsigned char**)malloc(5);
	srcY = (unsigned char**)malloc(5);
	srcU = (unsigned char**)malloc(5);
	srcV = (unsigned char**)malloc(5);
	for(int i=0; i<5; i++)
	{
		RGBbuf[i] = (unsigned char*)malloc(width*height*3);
		srcY[i] = (unsigned char*)malloc(YSize);
		srcU[i] = (unsigned char*)malloc(UVSize);
		srcV[i] = (unsigned char*)malloc(UVSize);
	}
}

CMainDlg::~CMainDlg()
{
	free_mem2D(RGBbuf);
	free_mem2D(srcY);
	free_mem2D(srcU);
	free_mem2D(srcV);

	fclose(DECFILE);
	fclose(BKSFILE);
	fclose(RESFILE);
	fclose(OMVFILE);
	//fclose(FMVFILE);
}

void CMainDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_FRMNUM, MaxFrmNum);
	DDX_Control(pDX, IDC_BUTTON_OPEN, m_OpenBtn);
	DDX_Control(pDX, IDC_BUTTON_ANAL, m_AnalBtn);
	DDX_Control(pDX, IDC_STATIC_ORIPIC, OriPic);
	DDX_Control(pDX, IDC_STATIC_BKSPIC, BksPic);
	DDX_Control(pDX, IDC_STATIC_RESPIC, ResPic);
	DDX_Control(pDX, IDC_STATIC_OMVPIC, OmvPic);
	DDX_Control(pDX, IDC_STATIC_FMVPIC, FmvPic);
	DDX_Text(pDX, IDC_EDIT_FRM, DisplayFrmNum);
}


BEGIN_MESSAGE_MAP(CMainDlg, CDialog)
	ON_BN_CLICKED(IDC_BUTTON_ANAL, &CMainDlg::OnBnClickedButtonAnal)
	ON_BN_CLICKED(IDC_BUTTON_OPEN, &CMainDlg::OnBnClickedButtonOpen)
	ON_BN_CLICKED(IDC_BUTTON_BKD, &CMainDlg::OnBnClickedButtonBkd)
	ON_BN_CLICKED(IDC_BUTTON_GO, &CMainDlg::OnBnClickedButtonGo)
	ON_BN_CLICKED(IDC_BUTTON_NEXT, &CMainDlg::OnBnClickedButtonNext)
	ON_BN_CLICKED(IDC_BUTTON_PREV, &CMainDlg::OnBnClickedButtonPrev)
	ON_WM_PAINT()
	ON_WM_CREATE()
	ON_BN_CLICKED(IDC_BUTTON_SHOW, &CMainDlg::OnBnClickedButtonShow)
	ON_STN_CLICKED(IDC_STATIC_OMV, &CMainDlg::OnStnClickedStaticOmv)
END_MESSAGE_MAP()


// CMainDlg 메시지 처리기입니다.

void CMainDlg::OnBnClickedButtonAnal()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	UpdateData(true);

	int i,j;

	//--------------------------------------------
	//설정
	//--------------------------------------------

	//일반 사항
	iheight							= height;
	iwidth							= width;
	max_motion_picture_num			= MaxFrmNum;	
	start_frame_num					= 0;

	//기존 알고리즘
	full_decode_flag				= 0;			//for the basic algorithm
	occlusion_handling				= 0;
	modify_tracking_flag			= 0;
	tracking_flag					= 0;			//0: disable / 1: tracking

	//새로운 알고리즘
	segmentation_flag				= 1;			//perform motion segmentation
	object_recognition_flag			= 1;			//allow the object recognition in the motion segmentation process
	trajectory_optimization_flag	= 1;			//allow trajectory optimization in every I frame by background subtraction
	decision_window_size			= 8;
	occur_th						= 1;			//the threshold for occurrence probability (값이 작을수록, 확률 임계값은 높아진다, 즉, 0이면 확률값은 1)

	//디코더 설정
	decode_residue_flag				= 0;			//decode residual data
	bkd_ref_flag					= 1;			//allow to use the background image for Intra prediction
	partial_decode_flag				= 1;			//for the segmentation algorithm

	//결과 분석 설정
	analysis_mode					= 0;
	video_save_flag					= 1;
	extract_flag					= 1;			//extract motion vectors and residual data
	res_anal_flag					= 0;			//extract residual data as a text file (just Y)
	read_mb_type_flag				= 1;			//extract and draw macroblock types
	extract_trajectory_flag			= 1;			//extract and draw the motion trajectory	
	output_trajectory_flag			= 1;			//record the trajectory information of objects (position, width, and height)
	output_group_buffer_flag		= 1;			//output the status of group buffers throughout all frames
	output_spatial_noise_flag		= 1;			//output the spatial noise rate throughout all frames
	time_measure_flag				= 0;			//for the segmentation algorithm

	//간편 설정 모음
	partial_decode_flag				= NO;			//부분 디코딩?
	time_measure_flag				= NO;			//처리시간 측정?
	
	if(time_measure_flag==1)	//시간 측정을 위한 최적화
	{
		extract_flag				= 0;
		res_anal_flag				= 0;
		extract_trajectory_flag		= 0;
		output_trajectory_flag		= 0;
		analysis_mode				= 1;
		output_group_buffer_flag	= 0;
		output_spatial_noise_flag	= 0;
	}

	//디스플레이 설정
	block_point_flag				= 0;			//show block points
	draw_textured_blocks			= 1;			//draw textured blocks in macroblock type images
	draw_motion_blocks				= 0;
	draw_textured_motion_blocks		= 0;
	multiple_block_drawing_flag		= 1;			//allow to draw blocks with multiple colors in macroblock type images
	residual_reverse_contrast_flag	= 1;			//make reverse contrast images of residual data

	//MB-Type 그리기 설정
	DrawMbTypeEnb[0]				= 0;
	DrawMbTypeEnb[1]				= 1;
	DrawMbTypeEnb[2]				= 1;
	DrawMbTypeEnb[3]				= 1;
	DrawMbTypeEnb[4]				= 1;
	DrawMbTypeEnb[5]				= 1;
	DrawMbTypeEnb[6]				= 1;
	DrawMbTypeEnb[7]				= 1;
	DrawMbTypeEnb[8]				= 1;
	DrawMbTypeEnb[9]				= 0;
	DrawMbTypeEnb[10]				= 0;
	DrawMbTypeEnb[11]				= 0;
	DrawMbTypeEnb[12]				= 0;
	DrawMbTypeEnb[13]				= 0;
	DrawMbTypeEnb[14]				= 0;
	DrawMbTypeEnb[15]				= 0;

	//기타 사항
	//th_textured						= 0.2;			//block texture property measure (1: no texture, 0: full texture)	

	//초기화
	max_buf_group_num				= 0;

	// allocate variables
	AllocBuf();

	//group buffer initialization
	for(i=0; i<MAX_BUF_NUM; i++)
	{
		gbuf[i].DetectFrame = -1;
		gbuf[i].Status = ERROR_GROUP;
		gbuf[i].CheckCnt = -1;
		gbuf[i].FrameCnt = -1;
		gbuf[i].ChkOn = 0;
		gbuf[i].DominantGroupIdx = -1;
		gbuf[i].OccurCnt = 0;
		for(j=0; j<MaxFrmNum; j++)
			gbuf[i].MbNum[j] = 0;
	}

	if(res_anal_flag==1)
		RESANAL_FILE = fopen("residual.txt","wb");
	if(extract_flag==1)
	{
		OMV_FILE = fopen("omv.yuv", "wb");
		RES_FILE = fopen("res.yuv", "wb");
	}
	if(video_save_flag==1)
		BKS_FILE = fopen("bks.yuv", "wb");
	if(read_mb_type_flag==1)
		MBT_FILE = fopen("mbt.yuv", "wb");
	if(extract_trajectory_flag==1)
		MTR_FILE = fopen("mtr.yuv", "wb");
	if(output_trajectory_flag==1)
		TRA_FILE = fopen("tra.dat", "wb");
	if(output_group_buffer_flag==1)
		BUF_FILE = fopen("buf.dat", "wb");
	if(output_spatial_noise_flag==1)
		SPA_FILE = fopen("spa.dat", "wb");

	//open analysis files
	if(analysis_mode==1)
	{
		//ET_FILE = fopen("energy.txt","wb");
		//WT_FILE = fopen("weight.txt","wb");
		TT_FILE = fopen("time.txt","wb");
	}

	if(video_save_flag==1)
	{
		int FrameSize = width*height;
		FILE *BKD_FILE = fopen(BkdPathName, "rb");
		fread(BkdBuf[0], sizeof(unsigned char), FrameSize, BKD_FILE);
		fread(BkdBuf[1], sizeof(unsigned char), FrameSize/4, BKD_FILE);
		fread(BkdBuf[2], sizeof(unsigned char), FrameSize/4, BKD_FILE);
		fclose(BKD_FILE);
	}
	
	// decode and extract information
	MainDec(PathName);

	// write result data
	//WriteDec();
	//WriteBks();
	//WriteRes();
	//WriteFmv();
	//WriteMvf();

	if(res_anal_flag==1)
		fclose(RESANAL_FILE);
	if(extract_flag==1)
	{
		fclose(OMV_FILE);
		fclose(RES_FILE);
	}
	if(video_save_flag==1)
		fclose(BKS_FILE);
	if(read_mb_type_flag==1)
		fclose(MBT_FILE);
	if(extract_trajectory_flag==1)
		fclose(MTR_FILE);
	if(output_trajectory_flag==1)
		fclose(TRA_FILE);
	if(output_group_buffer_flag==1)
		fclose(BUF_FILE);
	if(output_spatial_noise_flag==1)
		fclose(SPA_FILE);

	//close analysis files
	if(analysis_mode==1)
	{
		//fclose(ET_FILE);
		//fclose(WT_FILE);
		fclose(TT_FILE);
	}


	// free variables
	FreeBuf();

	AfxMessageBox("The analysis of an input video is finished.");

	UpdateData(false);
}

// The main function of decoding H.264/AVC videos
void CMainDlg::MainDec(CString filename)
{
	UpdateData(true);

	// ===========================================================================================================
	// Initialization of Decoding
	// ===========================================================================================================

	input = new struct inp_par [1];
	snr = new struct snr_par [1];
	img = new struct img_par [1];

	//Allocate memory for the structures
	char *config_filename=NULL;
	
	char fn[100];
	//strcpy(fn,((CFile *)ar.GetFile())->GetFileName());
	strcpy(fn, filename);
	strcpy(input->infile,fn);      //! set default bitstream name
	strcpy(input->outfile,"outdec.yuv"); //! set default output file name
	strcpy(input->reffile,"outrec.yuv"); //! set default reference file name
	input->FileFormat = PAR_OF_ANNEXB;
	input->ref_offset=0;
	input->poc_scale=1;
	
#ifdef _LEAKYBUCKET_
	input->R_decoder=500000;          //! Decoder rate
	input->B_decoder=104000;          //! Decoder buffer size
	input->F_decoder=73000;           //! Decoder initial delay
	strcpy(input->LeakyBucketParamFile,"leakybucketparam.cfg");    // file where Leaky Bucket params (computed by encoder) are stored
#endif

	config_filename = "decoder.cfg";
	init_conf(input, config_filename);		//Read input from configuration file

#if TRACEA
	if ((p_trace=fopen(TRACEFILE,"w"))==0)             // append new statistic at the end
	{
		snprintf(errortext, ET_SIZE, "Error open file %s!",TRACEFILE);
		error(errortext,500);
	}
#endif

	if ((p_out=open(input->outfile, OPENFLAGS_WRITE, OPEN_PERMISSIONS))==-1)
	{
		snprintf(errortext, ET_SIZE, "Error open file %s ",input->outfile);
		error(errortext,500);
	}

	p_ref=open(input->reffile,OPENFLAGS_READ);

#ifdef _LEAKYBUCKET_
	calc_buffer(input);
#endif

	init_old_slice();

	// Initialization of decoding
	bitss=fopen(input->infile, "rb");
	fseek(bitss, 0, SEEK_SET);

	// ===========================================================================================================
	// Main Decoding
	// ===========================================================================================================

	first_object_tracking_main();
	last_object_tracking_main();

	// ===========================================================================================================
	// Close Decoding
	// ===========================================================================================================

	active_sps = NULL;
	Co_located = NULL;
	erc_errorVar = NULL;

	UpdateData(false);
	Invalidate(false);
}

void CMainDlg::OnBnClickedButtonOpen()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	char szFilter[] = "h.264 (*.264)|*.264|all files (*.*)|*.*|";
	CFileDialog dlg(TRUE, NULL, NULL, OFN_FILEMUSTEXIST|OFN_PATHMUSTEXIST, szFilter);
		
	if(dlg.DoModal() == IDOK)
	{
		PathName = dlg.GetPathName();
	}
}

// Allocate buffers
void CMainDlg::AllocBuf(void)
{
	int i,j;
	int max_mb_nr;
	int b4_height, b4_width;
	int intra_period = 4;	//IPPP

	max_mb_nr = height*width/16/16;
	blk_dec_enb = (int*)calloc(max_mb_nr,sizeof(int));
	blk_list = (int**)calloc(max_mb_nr,sizeof(int*));
	for(i=0; i<max_mb_nr; i++)
	{
		blk_list[i] = (int*)calloc(16,sizeof(int));
	}

	if(tracking_flag==1)
	{
		upd_con = (unsigned char***)calloc(MaxFrmNum,sizeof(unsigned char**));
		for(i=0; i<MaxFrmNum; i++)
		{
			upd_con[i] = (unsigned char**)calloc(init_con_nr,sizeof(unsigned char*));
			for(j=0; j<init_con_nr; j++)
				upd_con[i][j] = (unsigned char*)calloc(2,sizeof(unsigned char));
		}

		object_center = (double**)calloc(MaxFrmNum,sizeof(double*));
		for(i=0; i<MaxFrmNum; i++)
			object_center[i] = (double*)calloc(2,sizeof(double));

		pred_con = (double**)calloc(init_con_nr,sizeof(double*));
		for(i=0; i<init_con_nr; i++)
			pred_con[i] = (double*)calloc(2,sizeof(double));
	}

	//output pictures
	if(video_save_flag==1)
	{
/*
		dec_pic = (unsigned char***)calloc(MaxFrmNum,sizeof(unsigned char**));
		for(i=0; i<MaxFrmNum; i++)
		{
			dec_pic[i] = (unsigned char**)calloc(3,sizeof(unsigned char*));
			for(j=0; j<3; j++)
				dec_pic[i][j] = (unsigned char*)calloc(height*width,sizeof(unsigned char));
		}
*/
		/*
		Y = new unsigned char[height*width];
		U = new unsigned char[height*width/4];
		V = new unsigned char[height*width/4];
		*/

		curr_frm = (unsigned char***)calloc(height,sizeof(unsigned char**));
		for(i=0; i<height; i++)
		{
			curr_frm[i] = (unsigned char**)calloc(width,sizeof(unsigned char*));
			for(j=0; j<width; j++)
				curr_frm[i][j] = (unsigned char*)calloc(3,sizeof(unsigned char));
		}

		prev_frm = (unsigned char***)calloc(height,sizeof(unsigned char**));
		for(i=0; i<height; i++)
		{
			prev_frm[i] = (unsigned char**)calloc(width,sizeof(unsigned char*));
			for(j=0; j<width; j++)
				prev_frm[i][j] = (unsigned char*)calloc(3,sizeof(unsigned char));
		}

		curr_res = (unsigned char***)calloc(height,sizeof(unsigned char**));
		for(i=0; i<height; i++)
		{
			curr_res[i] = (unsigned char**)calloc(width,sizeof(unsigned char*));
			for(j=0; j<width; j++)
				curr_res[i][j] = (unsigned char*)calloc(3,sizeof(unsigned char));
		}

		BkdBuf = (unsigned char**)calloc(3,sizeof(unsigned char*));
		BkdBuf[0] = (unsigned char*)calloc(height*width,sizeof(unsigned char));
		BkdBuf[1] = (unsigned char*)calloc(height*width/4,sizeof(unsigned char));
		BkdBuf[2] = (unsigned char*)calloc(height*width/4,sizeof(unsigned char));

		bkd_str = (unsigned char**)calloc(3,sizeof(unsigned char*));
		bkd_str[0] = (unsigned char*)calloc(height*width,sizeof(unsigned char));
		bkd_str[1] = (unsigned char*)calloc(height*width/4,sizeof(unsigned char));
		bkd_str[2] = (unsigned char*)calloc(height*width/4,sizeof(unsigned char));
	}

	if(extract_flag==1)
	{
/*
		dec_res = (unsigned char***)calloc(MaxFrmNum,sizeof(unsigned char**));
		for(i=0; i<MaxFrmNum; i++)
		{
			dec_res[i] = (unsigned char**)calloc(3,sizeof(unsigned char*));
			for(j=0; j<3; j++)
				dec_res[i][j] = (unsigned char*)calloc(height*width,sizeof(unsigned char));
		}
*/
		res_str = (unsigned char**)calloc(3,sizeof(unsigned char*));
		res_str[0] = (unsigned char*)calloc(height*width,sizeof(unsigned char));
		res_str[1] = (unsigned char*)calloc(height*width/4,sizeof(unsigned char));
		res_str[2] = (unsigned char*)calloc(height*width/4,sizeof(unsigned char));

		frm_omv = (unsigned char**)calloc(3,sizeof(unsigned char*));
		frm_omv[0] = (unsigned char*)calloc(width*height,sizeof(unsigned char));
		frm_omv[1] = (unsigned char*)calloc(width*height/4,sizeof(unsigned char));
		frm_omv[2] = (unsigned char*)calloc(width*height/4,sizeof(unsigned char));
	}

	if((extract_flag==1)||(tracking_flag==1))
	{
		b4_height = height/4;
		b4_width = width/4;

		upd_mv = (double***)calloc(b4_height,sizeof(double**));
		for(i=0; i<b4_height; i++)
		{
			upd_mv[i] = (double**)calloc(b4_width,sizeof(double*));
			for(j=0; j<b4_width; j++)
				upd_mv[i][j] = (double*)calloc(2,sizeof(double));
		}

		dec_mv = (double***)calloc(b4_height,sizeof(double**));
		for(i=0; i<b4_height; i++)
		{
			dec_mv[i] = (double**)calloc(b4_width,sizeof(double*));
			for(j=0; j<b4_width; j++)
				dec_mv[i][j] = (double*)calloc(2,sizeof(double));
		}
	}

	if(read_mb_type_flag==1)
	{
		MbType = (int*)calloc(max_mb_nr,sizeof(int));
		B8Mode = (int**)calloc(max_mb_nr,sizeof(int*));
		for(i=0; i<max_mb_nr; i++)
			B8Mode[i] = (int*)calloc(4,sizeof(int));

		frm_mbt = (unsigned char**)calloc(3,sizeof(unsigned char*));
		frm_mbt[0] = (unsigned char*)calloc(width*height,sizeof(unsigned char));
		frm_mbt[1] = (unsigned char*)calloc(width*height/4,sizeof(unsigned char));
		frm_mbt[2] = (unsigned char*)calloc(width*height/4,sizeof(unsigned char));
	}

	if((multiple_block_drawing_flag==1)&&(time_measure_flag==0))
	{
		DrawBlockFlag = (int**)calloc(10,sizeof(int*));
		for(i=0; i<10; i++)
			DrawBlockFlag[i] = (int*)calloc(b4_height*b4_width,sizeof(int));
	}

	if(segmentation_flag==1)
	{
		mblabel = (int**)calloc(MaxFrmNum,sizeof(int*));
		for(i=0; i<MaxFrmNum; i++)
			mblabel[i] = (int*)calloc(max_mb_nr,sizeof(int));

		for(i=0; i<MAX_BUF_NUM; i++)
		{
			gbuf[i].OccurProb = (double*)calloc(decision_window_size+1,sizeof(double));
			gbuf[i].FrmGroupIdx = (int*)calloc(MaxFrmNum,sizeof(int));
			gbuf[i].PosX = (int*)calloc(MaxFrmNum,sizeof(int));
			gbuf[i].PosY = (int*)calloc(MaxFrmNum,sizeof(int));
			gbuf[i].Width = (int*)calloc(MaxFrmNum,sizeof(int));
			gbuf[i].Height = (int*)calloc(MaxFrmNum,sizeof(int));
			gbuf[i].MbNum = (int*)calloc(MaxFrmNum,sizeof(int));
		}

		frm_buf = (unsigned char***)calloc(intra_period,sizeof(unsigned char**));
		for(i=0; i<intra_period; i++)
		{
			frm_buf[i] = (unsigned char**)calloc(3,sizeof(unsigned char*));
			frm_buf[i][0] = (unsigned char*)calloc(width*height,sizeof(unsigned char));
			frm_buf[i][1] = (unsigned char*)calloc(width*height/4,sizeof(unsigned char));
			frm_buf[i][2] = (unsigned char*)calloc(width*height/4,sizeof(unsigned char));
		}
	}
}

// Free buffers
void CMainDlg::FreeBuf(void)
{
	int max_mb_nr, b4_height, b4_width;
	int intra_period = 4;	//IPPP

	max_mb_nr = iheight*iwidth/16/16;
	b4_height = iheight/4;
	b4_width = iwidth/4;

	if(tracking_flag==1)
	{
		free(blk_dec_enb);
		free_mem2Dint(blk_list);
		free_mem2Ddouble(pred_con);
		free_mem3D(upd_con,max_motion_picture_num);
		free_mem2Ddouble(object_center);
	}

	if(extract_flag==1)
	{
		free_mem3D(curr_res,iheight);
		//free_mem3D(dec_res,max_motion_picture_num);
		free_mem2D(res_str);
		free_mem2D(frm_omv);
	}

	if((extract_flag==1)||(tracking_flag==1))
	{
		free_mem3Ddouble(upd_mv,b4_height);
		free_mem3Ddouble(dec_mv,b4_height);
	}

	if(video_save_flag==1)
	{
		free_mem3D(curr_frm,iheight);
		free_mem3D(prev_frm,iheight);
		//free_mem3D(dec_pic,max_motion_picture_num);
		free_mem2D(BkdBuf);
		free_mem2D(bkd_str);
		
		/*
		Y= NULL;
		U= NULL;
		V= NULL;
		delete Y;
		delete U;
		delete V;
		*/
	}

	if(read_mb_type_flag==1)
	{
		free(MbType);
		free_mem2Dint(B8Mode);
		free_mem2D(frm_mbt);
	}

	if((multiple_block_drawing_flag==1)&&(time_measure_flag==0))
	{
		free_mem2Dint(DrawBlockFlag);
	}

	if(segmentation_flag==1)
	{
		free_mem2Dint(mblabel);
		free_mem3D(frm_buf,intra_period);
	}
}

// Write an output video file
void CMainDlg::WriteDec(void)
{
	int i;

	if(video_save_flag==1)
	{
		/*
		if(fixed_box_mode==1)
		{
			ColHalfInterval = 25;
			RowHalfInterval = 25;
		}
		*/

		FILE* Write_FILE = fopen("dec.yuv","wb"); //write file (yuv sequence)
		//for(i=0; i<MaxFrmNum; i++)
		//	write_yuv_file(dec_pic[i][0],dec_pic[i][1],dec_pic[i][2],Write_FILE);
		fclose(Write_FILE);
	}
}

void CMainDlg::write_yuv_file(unsigned char *srcY, unsigned char *srcU, unsigned char *srcV, FILE* Write_FILE)
{
	int i;
/*
	int p;
	int cx=0, cy=0;
	int px[4], py[4];
	int thickness, color;
*/

	//arrcpy(Y,srcY,height*width);
	//arrcpy(U,srcU,(height>>1)*(width>>1));
	//arrcpy(V,srcV,(height>>1)*(width>>1));
	
/*
	cx = (int)(object_center[frame_number][0]+0.5);
	cy = (int)(object_center[frame_number][1]+0.5);

	if(track_box_flag==1)
	{
		if(fixed_box_mode==1)
		{
			px[0] = clipx(cx-ColHalfInterval); py[0] = clipy(cy-RowHalfInterval);
			px[1] = clipx(cx+ColHalfInterval); py[1] = clipy(cy-RowHalfInterval);
			px[2] = clipx(cx+ColHalfInterval); py[2] = clipy(cy+RowHalfInterval);
			px[3] = clipx(cx-ColHalfInterval); py[3] = clipy(cy+RowHalfInterval);
		}
		else
		{
			px[0] = clipx(cx-left_dis) ; py[0] = clipy(cy-top_dis);
			px[1] = clipx(cx-right_dis); py[1] = clipy(cy-top_dis);
			px[2] = clipx(cx-right_dis); py[2] = clipy(cy-bottom_dis);
			px[3] = clipx(cx-left_dis) ; py[3] = clipy(cy-bottom_dis);
		}

		//Draw the Network of Feature points
		thickness = 10; color = 100;
		for(p=0; p<init_con_nr; p++)
		{
			WritePoint(7,30, Y,U,V,upd_con[frame_number][p][0], upd_con[frame_number][p][1]);
			if(p<init_con_nr-1)
			{
				WriteLine(thickness, color, Y,U,V,frame_number,
					upd_con[frame_number][p][0], upd_con[frame_number][p][1],
					upd_con[frame_number][p+1][0], upd_con[frame_number][p+1][1]);
			}
		}
		
		//Draw Object Indicator Box
		thickness = 10; color = 200;
		WriteLine(thickness, color, Y,U,V,frame_number, px[0], py[0], px[1], py[1]);
		WriteLine(thickness, color, Y,U,V,frame_number, px[1], py[1], px[2], py[2]);
		WriteLine(thickness, color, Y,U,V,frame_number, px[2], py[2], px[3], py[3]);
		WriteLine(thickness, color, Y,U,V,frame_number, px[3], py[3], px[0], py[0]);
	}
*/

	//Write a file
	for(i = 0; i < height ; i++)
		fwrite(&(srcY[i*width]),1,width,Write_FILE );
		
	for(i = 0; i < (height/2); i++)
		fwrite(&(srcU[i*width/2]),1,width/2,Write_FILE );

	for(i = 0; i < (height/2); i++)
		fwrite(&(srcV[i*width/2]),1,width/2,Write_FILE );
}

void CMainDlg::arrcpy(unsigned char* DupArr, unsigned char* OrgArr, int size)
{
	for(int i=0; i<size; i++)
		DupArr[i] = (unsigned char) OrgArr[i];
}

void CMainDlg::OnBnClickedButtonBkd()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	char szFilter[] = "CIF Picture Format (*.cif)|*.cif|QCIF Picture Format (*.qcif)|*.qcif|YUV Picture Format (*.yuv)|*.yuv|all files (*.*)|*.*|";
	CFileDialog dlg(TRUE, NULL, NULL, OFN_FILEMUSTEXIST|OFN_PATHMUSTEXIST, szFilter);
		
	if(dlg.DoModal() == IDOK)
	{
		BkdPathName = dlg.GetPathName();
	}
}

void CMainDlg::OnBnClickedButtonGo()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	UpdateData(true);
	DisplayPic(DisplayFrmNum);
	OnPaint();
	UpdateData(false);
}

void CMainDlg::WriteRes(void)
{
	int i;

	if(video_save_flag==1)
	{
		if(extract_flag==1)
		{
			FILE* Write_FILE = fopen("res.yuv","wb"); //write file (yuv sequence)
			//for(i=0; i<MaxFrmNum; i++)
			//	write_yuv_file(dec_res[i][0],dec_res[i][1],dec_res[i][2],Write_FILE);
			fclose(Write_FILE);
		}
	}
}

void CMainDlg::WriteBks(void)
{
	if(video_save_flag==1)
	{
		int i, j;
		int YSize = width * height;
		int FrameSize = YSize * 1.5;
		int UVSize = YSize /4;
		int Voffset = YSize + UVSize;
		int tol = 25;

		// Get a background picture
		unsigned char* BkdBuf;
		BkdBuf = new unsigned char [FrameSize];	
		FILE *Bkd = fopen(BkdPathName, "rb");
		fread(BkdBuf, sizeof(unsigned char), FrameSize, Bkd);

		unsigned char* DecBuf;
		DecBuf = new unsigned char [FrameSize];	
		FILE *Video = fopen("dec.yuv", "rb");

		unsigned char** bkd_str;
		bkd_str = (unsigned char**)calloc(3,sizeof(unsigned char*));
		bkd_str[0] = (unsigned char*)calloc(height*width,sizeof(unsigned char));
		bkd_str[1] = (unsigned char*)calloc(height*width/4,sizeof(unsigned char));
		bkd_str[2] = (unsigned char*)calloc(height*width/4,sizeof(unsigned char));

		FILE* Write_FILE = fopen("bks.yuv","wb"); //write file (yuv sequence)

		for(i=0; i<MaxFrmNum; i++)
		{
			fread(DecBuf, sizeof(unsigned char), FrameSize, Video);

			for(j=0; j<YSize; j++)
			{
				if(IsNearlyZero(DecBuf[j] - BkdBuf[j],tol))
					bkd_str[0][j] = 16; //clip(DecBuf[j] - BkdBuf[j] +16);
				else
					bkd_str[0][j] = DecBuf[j];
			}
			for(j=0; j<UVSize; j++)
			{
				if(IsNearlyZero(DecBuf[j+YSize] - BkdBuf[j+YSize],tol))
					bkd_str[1][j] = 128; //clip(DecBuf[j+YSize] - BkdBuf[j+YSize] +128);
				else
					bkd_str[1][j] = DecBuf[j+YSize];
			}
			for(j=0; j<UVSize; j++)
			{
				if(IsNearlyZero(DecBuf[j+Voffset] - BkdBuf[j+Voffset],tol))
					bkd_str[2][j] = 128; //clip(DecBuf[j+Voffset] - BkdBuf[j+Voffset] +128);
				else
					bkd_str[2][j] = DecBuf[j+Voffset];
			}

			write_yuv_file(bkd_str[0],bkd_str[1],bkd_str[2],Write_FILE);
		}

		free(BkdBuf);
		free(DecBuf);
		free_mem2D(bkd_str);

		fclose(Write_FILE);
	}
}

void CMainDlg::OnBnClickedButtonNext()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	if(DisplayFrmNum<MaxFrmNum-1)
	{
		UpdateData(true);
		DisplayFrmNum++;
		DisplayPic(DisplayFrmNum);
		OnPaint();
		UpdateData(false);
	}
}

void CMainDlg::OnBnClickedButtonPrev()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	if(DisplayFrmNum>0)
	{
		UpdateData(true);
		DisplayFrmNum--;
		DisplayPic(DisplayFrmNum);
		OnPaint();
		UpdateData(false);
	}
}

void CMainDlg::OnPaint()
{
	CPaintDC dc(this); // device context for painting
	// TODO: 여기에 메시지 처리기 코드를 추가합니다.
	// 그리기 메시지에 대해서는 CDialog::OnPaint()을(를) 호출하지 마십시오.

	if(ShowFlag==1)
	{
		// Show an original picture
		CDC *pDC = OriPic.GetDC();
		CRect pRect;
		OriPic.GetWindowRect(pRect);

		BmpInfo->bmiHeader.biBitCount = 24;
		pDC->SetStretchBltMode(STRETCH_DELETESCANS);
		StretchDIBits(	pDC->m_hDC,
						0,0,pRect.Width(),pRect.Height(),				// destination;
						0,0,width,height,								// source;
						RGBbuf[0], BmpInfo,
						DIB_RGB_COLORS, SRCCOPY);

		// Show a background subtracted picture
		pDC = BksPic.GetDC();
		BksPic.GetWindowRect(pRect);

		BmpInfo->bmiHeader.biBitCount = 24;
		pDC->SetStretchBltMode(STRETCH_DELETESCANS);
		StretchDIBits(	pDC->m_hDC,
						0,0,pRect.Width(),pRect.Height(),				// destination;
						0,0,width,height,								// source;
						RGBbuf[1], BmpInfo,
						DIB_RGB_COLORS, SRCCOPY);

		// Show a residual picture
		pDC = ResPic.GetDC();
		ResPic.GetWindowRect(pRect);

		BmpInfo->bmiHeader.biBitCount = 24;
		pDC->SetStretchBltMode(STRETCH_DELETESCANS);
		StretchDIBits(	pDC->m_hDC,
						0,0,pRect.Width(),pRect.Height(),				// destination;
						0,0,width,height,								// source;
						RGBbuf[2], BmpInfo,
						DIB_RGB_COLORS, SRCCOPY);

		// Show original motion vectors
		pDC = OmvPic.GetDC();
		OmvPic.GetWindowRect(pRect);

		BmpInfo->bmiHeader.biBitCount = 24;
		pDC->SetStretchBltMode(STRETCH_DELETESCANS);
		StretchDIBits(	pDC->m_hDC,
						0,0,pRect.Width(),pRect.Height(),				// destination;
						0,0,width,height,								// source;
						RGBbuf[3], BmpInfo,
						DIB_RGB_COLORS, SRCCOPY);

		// Show filtered motion vectors
		pDC = FmvPic.GetDC();
		FmvPic.GetWindowRect(pRect);

		BmpInfo->bmiHeader.biBitCount = 24;
		pDC->SetStretchBltMode(STRETCH_DELETESCANS);
		StretchDIBits(	pDC->m_hDC,
						0,0,pRect.Width(),pRect.Height(),				// destination;
						0,0,width,height,								// source;
						RGBbuf[4], BmpInfo,
						DIB_RGB_COLORS, SRCCOPY);

	}
}

void CMainDlg::DisplayPic(int DisplayPicNum)
{
	int YSize = height*width;
	int UVSize = YSize/4;
	int FrameSize = YSize*1.5;

	// Decoded picture
	fseek(DECFILE, DisplayPicNum*FrameSize, SEEK_SET);
	fread(srcY[0], sizeof(unsigned char), YSize, DECFILE);
	fread(srcU[0], sizeof(unsigned char), UVSize, DECFILE);
	fread(srcV[0], sizeof(unsigned char), UVSize, DECFILE);
	pCon.YV12_to_RGB24(srcY[0],srcU[0],srcV[0],RGBbuf[0],width,height);

	// Background subtracted picture
	fseek(BKSFILE, DisplayPicNum*FrameSize, SEEK_SET);
	fread(srcY[1], sizeof(unsigned char), YSize, BKSFILE);
	fread(srcU[1], sizeof(unsigned char), UVSize, BKSFILE);
	fread(srcV[1], sizeof(unsigned char), UVSize, BKSFILE);
	pCon.YV12_to_RGB24(srcY[1],srcU[1],srcV[1],RGBbuf[1],width,height);

	// Residual picture
	fseek(RESFILE, DisplayPicNum*FrameSize, SEEK_SET);
	fread(srcY[2], sizeof(unsigned char), YSize, RESFILE);
	fread(srcU[2], sizeof(unsigned char), UVSize, RESFILE);
	fread(srcV[2], sizeof(unsigned char), UVSize, RESFILE);
	pCon.YV12_to_RGB24(srcY[2],srcU[2],srcV[2],RGBbuf[2],width,height);

	// Original motion vectors
	fseek(OMVFILE, DisplayPicNum*FrameSize, SEEK_SET);
	fread(srcY[3], sizeof(unsigned char), YSize, OMVFILE);
	fread(srcU[3], sizeof(unsigned char), UVSize, OMVFILE);
	fread(srcV[3], sizeof(unsigned char), UVSize, OMVFILE);
	pCon.YV12_to_RGB24(srcY[3],srcU[3],srcV[3],RGBbuf[3],width,height);

	// Macroblock types
	fseek(MBTFILE, DisplayPicNum*FrameSize, SEEK_SET);
	fread(srcY[4], sizeof(unsigned char), YSize, MBTFILE);
	fread(srcU[4], sizeof(unsigned char), UVSize, MBTFILE);
	fread(srcV[4], sizeof(unsigned char), UVSize, MBTFILE);
	pCon.YV12_to_RGB24(srcY[4],srcU[4],srcV[4],RGBbuf[4],width,height);

	ShowFlag = 1;
}

int CMainDlg::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CDialog::OnCreate(lpCreateStruct) == -1)
		return -1;

	// TODO:  여기에 특수화된 작성 코드를 추가합니다.
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

void CMainDlg::OnBnClickedButtonShow()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	DECFILE = fopen("dec.yuv", "rb");
	BKSFILE = fopen("bks.yuv", "rb");
	RESFILE = fopen("res.yuv", "rb");
	OMVFILE = fopen("omv.yuv", "rb");
	MBTFILE = fopen("mbt.yuv", "rb");
	//FMVFILE = fopen("fmv.yuv", "rb");

	UpdateData(true);
	DisplayPic(0);
	OnPaint();
	UpdateData(false);
}

void CMainDlg::OnStnClickedStaticOmv()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
}


/*!
 ***********************************************************************
 * \file image.c
 *
 * \brief
 *    Decode a Slice
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Byeong-Moon Jeon                <jeonbm@lge.com>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Gabi Blaettermann               <blaetter@hhi.de>
 *    - Ye-Kui Wang                     <wyk@ieee.org>
 *    - Antti Hallapuro                 <antti.hallapuro@nokia.com>
 *    - Alexis Tourapis                 <alexismt@ieee.org>
 *    - Jill Boyce                      <jill.boyce@thomson.net>
 *    - Saurav K Bandyopadhyay          <saurav@ieee.org>
 *    - Zhenyu Wu                       <Zhenyu.Wu@thomson.net
 *    - Purvin Pandit                   <Purvin.Pandit@thomson.net>
 *
 ***********************************************************************
 */

#include "contributors.h"

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "global.h"
#include "errorconcealment.h"
#include "image.h"
#include "mbuffer.h"
#include "fmo.h"
#include "nalu.h"
#include "parsetcommon.h"
#include "parset.h"
#include "header.h"
#include "rtpp.h"
#include "sei.h"
#include "output.h"
#include "biaridecod.h"
#include "mb_access.h"
#include "memalloc.h"
#include "annexb.h"

#include "context_ini.h"
#include "cabac.h"
#include "loopfilter.h"

#include "vlc.h"

//Added by Wonsang You at Nov-06-2006
#include <math.h>
#include <stdio.h>
#include <windows.h>
//Added-End

#include "erc_api.h"
extern objectBuffer_t *erc_object_list;
extern ercVariables_t *erc_errorVar;
extern frame erc_recfr;
extern int erc_mvperMB;
extern struct img_par *erc_img;

//extern FILE *p_out2;

extern StorablePicture **listX[6];
extern ColocatedParams *Co_located;

StorablePicture *dec_picture;

OldSliceParams old_slice;

//Added by Wonsang You at FEB 10, 2007
struct FeaturePoint
{
	int x;
	int y;
	double te;
	double me;
	double fe;
};
//Added-End

void MbAffPostProc()
{
  imgpel temp[16][32];

  imgpel ** imgY  = dec_picture->imgY;
  imgpel ***imgUV = dec_picture->imgUV;

  int i, x, y, x0, y0, uv;
  for (i=0; i<(int)dec_picture->PicSizeInMbs; i+=2)
  {
    if (dec_picture->mb_field[i])
    {
      get_mb_pos(i, &x0, &y0);
      for (y=0; y<(2*MB_BLOCK_SIZE);y++)
        for (x=0; x<MB_BLOCK_SIZE; x++)
          temp[x][y] = imgY[y0+y][x0+x];

      for (y=0; y<MB_BLOCK_SIZE;y++)
        for (x=0; x<MB_BLOCK_SIZE; x++)
        {
          imgY[y0+(2*y)][x0+x]   = temp[x][y];
          imgY[y0+(2*y+1)][x0+x] = temp[x][y+MB_BLOCK_SIZE];
		  //Added by Wonsang You at FEB 6, 2007
		  curr_frm[y0+(2*y)][x0+x][0] = (unsigned char)temp[x][y];
		  curr_frm[y0+(2*y+1)][x0+x][0] = (unsigned char)temp[x][y+MB_BLOCK_SIZE];
		  //Added-End
        }

      if (dec_picture->chroma_format_idc != YUV400)
      {
        x0 = x0 / (16/img->mb_cr_size_x);
        y0 = y0 / (16/img->mb_cr_size_y);

        for (uv=0; uv<2; uv++)
        {
          for (y=0; y<(2*img->mb_cr_size_y);y++)
            for (x=0; x<img->mb_cr_size_x; x++)
              temp[x][y] = imgUV[uv][y0+y][x0+x];
          
          for (y=0; y<img->mb_cr_size_y;y++)
            for (x=0; x<img->mb_cr_size_x; x++)
            {
              imgUV[uv][y0+(2*y)][x0+x]   = temp[x][y];
              imgUV[uv][y0+(2*y+1)][x0+x] = temp[x][y+img->mb_cr_size_y];
			  //Added by Wonsang You at FEB 6, 2007
			  curr_frm[y0+(2*y)][x0+x][uv+1] = (unsigned char)temp[x][y];
			  curr_frm[y0+(2*y+1)][x0+x][uv+1] = (unsigned char)temp[x][y+img->mb_cr_size_y];
			  //Added-End
            }
        }
      }
    }
  }
}

/*!
 ***********************************************************************
 * \brief
 *    decodes one I- or P-frame
 *
 ***********************************************************************
 */

int decode_one_frame(struct img_par *img,struct inp_par *inp, struct snr_par *snr)
{
  int current_header;
  Slice *currSlice = img->currentSlice;

  img->current_slice_nr = 0;
  img->current_mb_nr = -4711;     // initialized to an impossible value for debugging -- correct value is taken from slice header
  currSlice->next_header = -8888; // initialized to an impossible value for debugging -- correct value is taken from slice header
  img->num_dec_mb = 0;
  img->newframe = 1;

  while ((currSlice->next_header != EOS && currSlice->next_header != SOP))
  {
    current_header = read_new_slice();

    if (current_header == EOS)
    {
      exit_picture();
      return EOS;
    }

    decode_slice(img, inp, current_header);

    img->newframe = 0;
    img->current_slice_nr++;

  }

  exit_picture();
    
  return (SOP);
}


/*!
 ***********************************************************************
 * \brief
 *    partial decoding of one I- or P-frame
 *	  Added by Wonsang You at JAN 31, 2007
 *
 ***********************************************************************
 */

int partial_decode_one_frame(struct img_par *img,struct inp_par *inp, struct snr_par *snr)
{
  int current_header;
  Slice *currSlice = img->currentSlice;

  img->current_slice_nr = 0;
  img->current_mb_nr = -4711;     // initialized to an impossible value for debugging -- correct value is taken from slice header
  currSlice->next_header = -8888; // initialized to an impossible value for debugging -- correct value is taken from slice header
  img->num_dec_mb = 0;
  img->newframe = 1;

  //Modified by Wonsang You at JAN 31, 2007
  while ((currSlice->next_header != EOS && currSlice->next_header != SOP))
  {
	current_header = partial_read_new_slice();

	if (current_header == EOS)
	{
	  exit_picture();
	  return EOS;
	}

	gen_decblock_list();
	partial_decode_slice(img, inp, current_header);		//partial decoding
	pred_pos();
	decide_pos();
	set_buffer_PD();
	
	img->newframe = 0;
	img->current_slice_nr++;

  }
  //Modified-End

  exit_picture();
    
  return (SOP);

}

int last_partial_decode_one_frame(struct img_par *img,struct inp_par *inp, struct snr_par *snr)
{
  int i,j;
  int current_header;
  Slice *currSlice = img->currentSlice;
  LARGE_INTEGER freq, g1, g2;
  double diff;

  int intra_period = INTRA_PERIOD;	//IPPP
  int curr_frame_type = motion_picture_num%intra_period;
  int max_mb_nr = (int)(iwidth*iheight/16/16);
 
  if(tracking_flag==1)
  {
	  pred_pos();
	  decide_pos();
  }
  else
  {
	if(extract_flag==1)
	{
		if(motion_picture_num>0)
			partial_store_motion_into_buffer();
	}
  }

  if(extract_flag==1)
	extract_fmv();
  if(segmentation_flag==1)
  {
	//initialize decoded block flag
	for(i=0; i<max_mb_nr; i++)
	{
		blk_dec_enb[i] = 0;
		for(j=0; j<16; j++)
			blk_list[i][j] = 0;
	}

    motion_segment();

	if(time_measure_flag==0)
	{
		//Write a file
		for(i=0; i<iheight; i++)
		{
			for(j=0; j<iwidth; j++)
			{
				frm_buf[0][0][i*iwidth+j] = curr_frm[i][j][0];
				frm_buf[0][1][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][1];
				frm_buf[0][2][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][2];
			}
		}

		for(i = 0; i < iheight ; i++)
			fwrite(&(frm_buf[0][0][i*iwidth]),1,iwidth,MTR_FILE );

		for(i = 0; i < (iheight/2); i++)
			fwrite(&(frm_buf[0][1][i*iwidth/2]),1,iwidth/2,MTR_FILE );

		for(i = 0; i < (iheight/2); i++)
			fwrite(&(frm_buf[0][2][i*iwidth/2]),1,iwidth/2,MTR_FILE );
	}

  }
  else
  {
	  if((read_mb_type_flag==1)&&(time_measure_flag==0))
		draw_mb_type();
  }

  set_buffer_PD();

  img->newframe = 0;
  img->current_slice_nr++;

  //Added by Wonsang You at FEB 25, 2007
  first_frame = 0;
  //Added-End

  //Modified by Wonsang You at FEB 26, 2007
  while ((currSlice->next_header != EOS && currSlice->next_header != SOP)&&(motion_picture_num<max_motion_picture_num-1))
  {
	current_header = partial_read_new_slice();

	if (current_header == EOS)
	{
	  exit_picture();
	  return EOS;
	}

	//START: measure the computational time
	if(analysis_mode==1)
	{
		if(modify_tracking_flag==0)
		{
			QueryPerformanceFrequency(&freq);
			QueryPerformanceCounter(&g1);
		}
	}

	if((tracking_flag==1)&&(segmentation_flag==0))
		gen_decblock_list();
	else if((tracking_flag==0)&&(segmentation_flag==1))
	{
		curr_frame_type = motion_picture_num%intra_period;
		if(curr_frame_type==0)
			pred_segment();		//객체 위치 및 영역 예측
		gen_segment_decblock_list();
	}

	partial_decode_slice(img, inp, current_header);		//partial decoding

	//---the another way for random access
	//if(motion_picture_num==3)
	//	partial_decode_slice(img, inp, current_header);		//partial decoding
	//else	
	//	exit_slice();
	
	if(tracking_flag==1)
	{
		pred_pos();
		decide_pos();
	}
	else
	{
		if(extract_flag==1)
			partial_store_motion_into_buffer();
	}

	if(extract_flag==1)
		extract_fmv();
	if(segmentation_flag==1)
		motion_segment();
	else
	{
		if((read_mb_type_flag==1)&&(time_measure_flag==0))
			draw_mb_type();
	}

	//END: measure the computational time
	if(analysis_mode==1)
	{
		if(modify_tracking_flag==0)
		{
			QueryPerformanceCounter(&g2);
			diff = (double)(g2.QuadPart - g1.QuadPart)/ freq.QuadPart;	
			fprintf(TT_FILE,"%d\t%f\n", motion_picture_num, diff);
		}	
	}
	
	set_buffer_PD();

	img->newframe = 0;
	img->current_slice_nr++;

	if(motion_picture_num==max_motion_picture_num-1)
	{
	  return EOS;
	}
  }
  //Modified-End

  exit_picture();
    
  return (SOP);

}

int first_partial_decode_one_frame(struct img_par *img,struct inp_par *inp, struct snr_par *snr)
{
  int i,j;
  int current_header;
  Slice *currSlice = img->currentSlice;
 
  int tol = 20;
  int YSize = iheight*iwidth;
  int Voffset = YSize + YSize/4;

  img->current_slice_nr = 0;
  img->current_mb_nr = -4711;     // initialized to an impossible value for debugging -- correct value is taken from slice header
  currSlice->next_header = -8888; // initialized to an impossible value for debugging -- correct value is taken from slice header
  img->num_dec_mb = 0;
  img->newframe = 1;

  //Added by Wonsang You at FEB 25, 2007
  first_frame = 1;
  //Added-End

  current_header = partial_read_new_slice();

  if (current_header == EOS)
  {
    exit_picture();
    return EOS;
  }

  if((tracking_flag==1)&&(segmentation_flag==0))
		gen_decblock_list();

  partial_decode_slice(img, input, current_header);		//partial decoding

  if(video_save_flag==1)
  {
	  /*
	  // Save the background subtracted picture
		for(i=0; i<iheight; i++)
		{
			for(j=0; j<iwidth; j++)
			{
				if(IsNearlyZero(curr_frm[i][j][0] - BkdBuf[0][i*iwidth+j],tol))
					bkd_str[0][i*iwidth+j] = 16; //clip(DecBuf[j] - BkdBuf[j] +16);
				else
					bkd_str[0][i*iwidth+j] = curr_frm[i][j][0];

				if(IsNearlyZero(curr_frm[i>>1][j>>1][1] - BkdBuf[0][(i>>1)*(iwidth>>1)+(j>>1)+YSize],tol))
					bkd_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = 128; //clip(DecBuf[j+YSize] - BkdBuf[j+YSize] +128);
				else
					bkd_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][1];

				if(IsNearlyZero(curr_frm[i>>1][j>>1][2] - BkdBuf[0][(i>>1)*(iwidth>>1)+(j>>1)+Voffset],tol))
					bkd_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = 128; //clip(DecBuf[j+Voffset] - BkdBuf[j+Voffset] +128);
				else
					bkd_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][2];
			}
		}

		for(i = 0; i < iheight ; i++)
			fwrite(&(bkd_str[0][i*iwidth]),1,iwidth,BKS_FILE );				
		for(i = 0; i < (iheight/2); i++)
			fwrite(&(bkd_str[1][i*iwidth/2]),1,iwidth/2,BKS_FILE );
		for(i = 0; i < (iheight/2); i++)
			fwrite(&(bkd_str[2][i*iwidth/2]),1,iwidth/2,BKS_FILE );
		*/

	  if(extract_flag==1)
	  {
		  for(i=0; i<iheight; i++)
		  {
			for(j=0; j<iwidth; j++)
			{
				if(residual_reverse_contrast_flag==1)
				{
					res_str[0][i*iwidth+j] = max(0,min(255,ContrastStretch(10,curr_res[i][j][0]) +16));
					res_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = 128;
					res_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = 128;
				}
				else
				{
					res_str[0][i*iwidth+j] = max(0,min(255,curr_res[i][j][0] +16));
					res_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = max(0,min(255,curr_res[i>>1][j>>1][1] +128));
					res_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = max(0,min(255,curr_res[i>>1][j>>1][2] +128));
				}
			}
		  }

		  for(i = 0; i < iheight ; i++)
			fwrite(&(res_str[0][i*iwidth]),1,iwidth,RES_FILE );				
		  for(i = 0; i < (iheight/2); i++)
			fwrite(&(res_str[1][i*iwidth/2]),1,iwidth/2,RES_FILE );
		  for(i = 0; i < (iheight/2); i++)
			fwrite(&(res_str[2][i*iwidth/2]),1,iwidth/2,RES_FILE );
	  }
  }
  
  return 1;

}


/*!
 ***********************************************************************
 * \brief
 *    allocate buffers for partial decoding
 *	  Added by Wonsang You at FEB 2, 2007
 *
 ***********************************************************************
 */

void alloc_buffer_PD()
{
	int max_mb_nr;
	int b4_height, b4_width;
	int i,j;

	max_mb_nr = iheight*iwidth/16/16;
	b4_height = iheight/4;
	b4_width = iwidth/4;

	blk_dec_enb = (int*)calloc(max_mb_nr,sizeof(int));
	blk_list = (int**)calloc(max_mb_nr,sizeof(int*));
	for(i=0; i<max_mb_nr; i++)
	{
		blk_list[i] = (int*)calloc(16,sizeof(int));
	}

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

	curr_frm = (unsigned char***)calloc(iheight,sizeof(unsigned char**));
	for(i=0; i<iheight; i++)
	{
		curr_frm[i] = (unsigned char**)calloc(iwidth,sizeof(unsigned char*));
		for(j=0; j<iwidth; j++)
			curr_frm[i][j] = (unsigned char*)calloc(3,sizeof(unsigned char));
	}

	prev_frm = (unsigned char***)calloc(iheight,sizeof(unsigned char**));
	for(i=0; i<iheight; i++)
	{
		prev_frm[i] = (unsigned char**)calloc(iwidth,sizeof(unsigned char*));
		for(j=0; j<iwidth; j++)
			prev_frm[i][j] = (unsigned char*)calloc(3,sizeof(unsigned char));
	}

	pred_con = (double**)calloc(init_con_nr,sizeof(double*));
    for(i=0; i<init_con_nr; i++)
		pred_con[i] = (double*)calloc(2,sizeof(double));

}


/*!
 ***********************************************************************
 * \brief
 *    delete buffers for partial decoding
 *	  Added by Wonsang You at FEB 2, 2007
 *
 ***********************************************************************
 */

void free_buffer_PD()
{
	int max_mb_nr, b4_height, b4_width;

	max_mb_nr = iheight*iwidth/16/16;
	b4_height = iheight/4;
	b4_width = iwidth/4;

	free(blk_dec_enb);
	free_mem2Dint(blk_list);
	free_mem3Ddouble(upd_mv,b4_height);	//double
	free_mem3Ddouble(dec_mv,b4_height);	//short
	free_mem3D(curr_frm,iheight);
	free_mem3D(prev_frm,iheight);
	free_mem2Ddouble(pred_con);		//double
	

}


/*!
 ***********************************************************************
 * \brief
 *    get forward motion vectors and predict positions of feature points
 *	  Added by Wonsang You at FEB 26, 2007
 *
 ***********************************************************************
 */

int pred_pos()
{
	int point;
	int block_x, block_y;
	double tmp_x, tmp_y;

	//predict positions of feature points
	if(motion_picture_num>0)
	{
		partial_store_motion_into_buffer();	

		point = 0;
		for(point=0; point<init_con_nr; point++)
		{
			block_x = (int)(upd_con[motion_picture_num-1][point][0]/4);
			block_y = (int)(upd_con[motion_picture_num-1][point][1]/4);

			tmp_x = upd_con[motion_picture_num-1][point][0] + upd_mv[block_y][block_x][0];
			tmp_y = upd_con[motion_picture_num-1][point][1] + upd_mv[block_y][block_x][1];

			pred_con[point][0] = clipx(tmp_x);
			pred_con[point][1] = clipy(tmp_y);
		}
	}

	return 1;
}

double clipx(double x)
{
	double out;

	if(x<0)
		out=0;
	else if(x>=iwidth)
		out=iwidth-1;
	else
		out=x;

	return out;
}

double clipy(double y)
{
	double out;

	if(y<0)
		out=0;
	else if(y>=iheight)
		out=iheight-1;
	else
		out=y;

	return out;
}

double clipbx(double bx)
{
	double b_width = iwidth/4;		//320%4
	double out;

	if(bx<0)
		out=0;
	else if(bx>=b_width)
		out=b_width-1;
	else
		out=bx;

	return out;
}

double clipby(double by)
{
	double b_height = iheight/4;		//240%4
	double out;

	if(by<0)
		out=0;
	else if(by>=b_height)
		out=b_height-1;
	else
		out=by;

	return out;
}


/*!
 ***********************************************************************
 * \brief
 *    generate the decoded block list
 *	  Added by Wonsang You at JAN 31, 2007
 *
 ***********************************************************************
 */

int gen_decblock_list()
{
	int f,p,i,j,l,r;
	int intra_period, curr_frame_type, max_mb_nr;
	int mb_nr, x4_nr;	
	int multi, dif_mv_multi, unit_half_interval;
	int start_block_x, start_block_y, end_block_x, end_block_y;
	double dif_mv[2];
	double center_px, center_py;
	int orig_bx, orig_by;
	int SearchErr;
	double pred_con_x, pred_con_y;

	//set parameters
	intra_period = INTRA_PERIOD;
	curr_frame_type = motion_picture_num%intra_period;
	max_mb_nr = (int)(iwidth*iheight/16/16);

	///////////////////////////////////////////////////////////////////////
	//main algorithm
	///////////////////////////////////////////////////////////////////////
	if((curr_frame_type==0)||(full_decode_flag==1)||(tracking_flag==0))		//full decoded blocks in I-frames
	{
		for(l=0; l<max_mb_nr; l++)
		{
			blk_dec_enb[l] = 1;
			for(r=0; r<16; r++)
				blk_list[l][r] = 1;
		}
	}
	else		//partial decoded blocks in P-frames
	{
		//Modified Simple Algorithm
		if(motion_picture_num==1)		//The first P-frame has no forward motion vectors.
		{
			for(l=0; l<max_mb_nr; l++)
			{
				blk_dec_enb[l] = 1;
				for(r=0; r<16; r++)
					blk_list[l][r] = 1;
			}
		}
		else
		{
			for(f=curr_frame_type; f<intra_period; f++)
			{
				multi = intra_period - f;
				dif_mv_multi = (intra_period-1) - f;
				if(f==intra_period-1)
					SearchErr = 0;
				else
					SearchErr = dec_search_err;
				unit_half_interval = multi*(search_half_interval+SearchErr)+sad_half_interval;

				for(p=0; p<init_con_nr; p++)
				{
					orig_bx = (int)(upd_con[motion_picture_num-2][p][0]/4);
					orig_by = (int)(upd_con[motion_picture_num-2][p][1]/4);
					pred_con_x = (double)(upd_con[motion_picture_num-1][p][0] + upd_mv[orig_by][orig_bx][0]);
					pred_con_y = (double)(upd_con[motion_picture_num-1][p][1] + upd_mv[orig_by][orig_bx][1]);
					dif_mv[0] = (upd_mv[orig_by][orig_bx][0] + dec_mv[orig_by][orig_bx][0]/4);
					dif_mv[1] = (upd_mv[orig_by][orig_bx][1] + dec_mv[orig_by][orig_bx][1]/4);
					center_px = pred_con_x + dif_mv_multi*dif_mv[0];
					center_py = pred_con_y + dif_mv_multi*dif_mv[1];

					start_block_x = (int)clipbx((center_px - unit_half_interval)/4);
					start_block_y = (int)clipby((center_py - unit_half_interval)/4);
					end_block_x = (int)clipbx((center_px + unit_half_interval)/4);
					end_block_y = (int)clipby((center_py + unit_half_interval)/4);

					for(i=start_block_x; i<=end_block_x; i++)
					{
						for(j=start_block_y; j<=end_block_y; j++)
						{
							trans_block_idx(i,j,&mb_nr,&x4_nr);
							blk_dec_enb[mb_nr] = 1;
							blk_list[mb_nr][x4_nr] = 1;
						}
					}
				}
			}
		}
	}

	return 1;

}


/*!
 ***********************************************************************
 * \brief
 *    generate the decoded block list for motion segmentation
 *	  Added by Wonsang You at OCT 15, 2007
 *
 ***********************************************************************
 */

int gen_segment_decblock_list()
{
	int k,r,l,i,j;
	int min_x, max_x, min_y, max_y;
	int min_mbx, max_mbx, min_mby, max_mby;
	int mbx, mby, mb_nr;

	//set parameters
	int hor_decoding_margin = HOR_DECODING_MARGIN;
	int ver_decoding_margin = VER_DECODING_MARGIN;	
	int intra_period = INTRA_PERIOD;
	int curr_frame_type = motion_picture_num%intra_period;
	int mbwidth = (int)(iwidth/16);
	int max_mb_nr = (int)(iwidth*iheight/16/16);

	if(curr_frame_type==0)		//for I-frames
	{
		if(partial_decode_flag==1)
		{
			//initialize decoded block flag
			for(i=0; i<max_mb_nr; i++)
			{
				blk_dec_enb[i] = 0;
				for(j=0; j<16; j++)
					blk_list[i][j] = 0;
			}

			for(k=0; k<max_buf_group_num; k++)
			{
				if(gbuf[k].Status==OBJECT_GROUP)		//각 object에 대하여
				{
					//왼쪽 위 좌표의 매크로블록 좌표, 오른쪽 아래 좌표의 매크로블록 좌표 얻기
					min_x = ClipX(gbuf[k].PosX[gbuf[k].FrameCnt+1] - gbuf[k].Width[gbuf[k].FrameCnt+1]/2 - hor_decoding_margin);
					max_x = ClipX(gbuf[k].PosX[gbuf[k].FrameCnt+1] + gbuf[k].Width[gbuf[k].FrameCnt+1]/2 + hor_decoding_margin);
					min_y = ClipY(gbuf[k].PosY[gbuf[k].FrameCnt+1] - gbuf[k].Height[gbuf[k].FrameCnt+1]/2 - ver_decoding_margin);
					max_y = ClipY(gbuf[k].PosY[gbuf[k].FrameCnt+1] + gbuf[k].Height[gbuf[k].FrameCnt+1]/2 + ver_decoding_margin);

					min_mbx = (int)(min_x/16);
					max_mbx = (int)(max_x/16);
					min_mby = (int)(min_y/16);
					max_mby = (int)(max_y/16);

					for(mby=min_mby; mby<=max_mby; mby++)
					{
						for(mbx=min_mbx; mbx<=max_mbx; mbx++)
						{
							//객체 내의 각 매크로블록에 대하여 매크로블록 번호 추출
							mb_nr = mby*mbwidth + mbx;
							blk_dec_enb[mb_nr] = 1;
							for(r=0; r<16; r++)
								blk_list[mb_nr][r] = 1;
						}
					}				
				}
			}
		}
		else		//full decoding
		{
			for(l=0; l<max_mb_nr; l++)
			{
				blk_dec_enb[l] = 1;
				for(r=0; r<16; r++)
					blk_list[l][r] = 1;
			}
		}
	}
	else		//for P-frames
	{
		for(l=0; l<max_mb_nr; l++)
		{
			blk_dec_enb[l] = 1;
			for(r=0; r<16; r++)
				blk_list[l][r] = 1;
		}
	}

	return 1;

}


/*!
 ***********************************************************************
 * \brief
 *    translate block coordiate into block index
 *	  Added by Wonsang You at FEB 01, 2007
 *
 ***********************************************************************
 */

void trans_block_idx(int block_x, int block_y, int *mb_nr, int *x4_nr)
{
	int mbx, mby, mb_rem_x, mb_rem_y, x8_blk_x, x8_blk_y, x8_rem_x, x8_rem_y, x8_nr;

	mbx = (int)(block_x / 4);
	mby = (int)(block_y / 4);
	mb_rem_x = (int)(block_x % 4);
	mb_rem_y = (int)(block_y % 4);
	x8_blk_x = (int)(mb_rem_x / 2);
	x8_blk_y = (int)(mb_rem_y / 2);
	x8_rem_x = (int)(mb_rem_x % 2);
	x8_rem_y = (int)(mb_rem_y % 2);
	x8_nr = (int)(2*x8_blk_y + x8_blk_x);

	*x4_nr = (int)(x8_nr*4 + 2*x8_rem_y + x8_rem_x);
	*mb_nr = (int)(iwidth/16)*mby + mbx;
}


/*!
 ***********************************************************************
 * \brief
 *    translate block index into block coordiate
 *	  Added by Wonsang You at FEB 01, 2007
 *
 ***********************************************************************
 */

void trans_block_coordinate(int mb_nr, int x4_nr, int *block_x, int *block_y)
{
	int mbx, mby, x8_nr, x8_rem_nr, x8_blk_x, x8_blk_y, x8_rem_x, x8_rem_y;

	mbx = (int)(mb_nr%(iwidth/16));
	mby = (int)(mb_nr/(iwidth/16));

	x8_nr = (int)(x4_nr/4);
	x8_rem_nr = (int)(x4_nr%4);

	x8_blk_x = (int)(x8_nr%2);
	x8_blk_y = (int)(x8_nr/2);
	x8_rem_x = (int)(x8_rem_nr%2);
	x8_rem_y = (int)(x8_rem_nr/2);

	*block_x = (int)(mbx*4 + x8_blk_x*2 + x8_rem_x);
	*block_y = (int)(mbx*4 + x8_blk_y*2 + x8_rem_y);

}


/*!
 ***********************************************************************
 * \brief
 *    get reference blocks for one decoded block
 *	  Added by Wonsang You at FEB 01, 2007
 *
 ***********************************************************************
 */

void get_ref_blocks(int block_x, int block_y, int **neighbor)
{
	double ori[4][2], c[4][2];
	int mid_line_bx, mid_line_by, mid_line_x, mid_line_y, check_4multi_x, check_4multi_y;

	//Get four corners of an original block
	ori[0][0] = (double)block_x*4;
	ori[0][1] = (double)block_y*4;
	ori[1][0] = (double)(block_x+1)*4;
	ori[1][1] = (double)block_y*4;
	ori[2][0] = (double)block_x*4;
	ori[2][1] = (double)(block_y+1)*4;
	ori[3][0] = (double)(block_x+1)*4;
	ori[3][1] = (double)(block_y+1)*4;

	//Get four corners of a predicted block
	c[0][0] = ori[0][0] + (double)dec_mv[block_y][block_x][0]/4;
	c[0][1] = ori[0][1] + (double)dec_mv[block_y][block_x][1]/4;
	c[1][0] = ori[1][0] + (double)dec_mv[block_y][block_x][0]/4;
	c[1][1] = ori[1][1] + (double)dec_mv[block_y][block_x][1]/4;
	c[2][0] = ori[2][0] + (double)dec_mv[block_y][block_x][0]/4;
	c[2][1] = ori[2][1] + (double)dec_mv[block_y][block_x][1]/4;
	c[3][0] = ori[3][0] + (double)dec_mv[block_y][block_x][0]/4;
	c[3][1] = ori[3][1] + (double)dec_mv[block_y][block_x][1]/4;

	//Check which blocks are overlapped with a predicted block
	mid_line_bx = (int)((int)c[1][0]/4);
	mid_line_by = (int)((int)c[2][1]/4);
	mid_line_x = mid_line_bx*4;
	mid_line_y = mid_line_by*4;

	check_4multi_x = (int)((int)c[0][0]) % 4;
	check_4multi_y = (int)((int)c[0][1]) % 4;

	if(check_4multi_x == 0)
	{
		if(check_4multi_y != 0)
		{
			neighbor[0][0] = (mid_line_bx - 1);
			neighbor[0][1] = (mid_line_by - 1);
			neighbor[1][0] = (mid_line_bx - 1);
			neighbor[1][1] = (mid_line_by);
			neighbor[2][0] = -1;
			neighbor[2][1] = -1;
			neighbor[3][0] = -1;
			neighbor[3][1] = -1;
		}
		else
		{
			neighbor[0][0] = (mid_line_bx - 1);
			neighbor[0][1] = (mid_line_by - 1);
			neighbor[1][0] = -1;
			neighbor[1][1] = -1;
			neighbor[2][0] = -1;
			neighbor[2][1] = -1;
			neighbor[3][0] = -1;
			neighbor[3][1] = -1;
		}
	}	
	else if(check_4multi_y == 0)
	{
		neighbor[0][0] = (mid_line_bx - 1);
		neighbor[0][1] = (mid_line_by - 1);
		neighbor[1][0] = (mid_line_bx);
		neighbor[1][1] = (mid_line_by - 1);
		neighbor[2][0] = -1;
		neighbor[2][1] = -1;
		neighbor[3][0] = -1;
		neighbor[3][1] = -1;
	}
	else if((mid_line_x > c[0][0])&&(mid_line_y > c[0][1]))
	{
		neighbor[0][0] = (mid_line_bx - 1);
		neighbor[0][1] = (mid_line_by - 1);
		neighbor[1][0] = (mid_line_bx);
		neighbor[1][1] = (mid_line_by - 1);
		neighbor[2][0] = (mid_line_bx - 1);
		neighbor[2][1] = (mid_line_by);
		neighbor[3][0] = (mid_line_bx);
		neighbor[3][1] = (mid_line_by);
	}
	else
	{
		//AfxMessageBox("Neighbor blocks are not detected for a motion vector");
	}

	//Availability
	neighbor[0][2] = block_availability(neighbor[0][0],neighbor[0][1]);
	neighbor[1][2] = block_availability(neighbor[1][0],neighbor[1][1]);
	neighbor[2][2] = block_availability(neighbor[2][0],neighbor[2][1]);
	neighbor[3][2] = block_availability(neighbor[3][0],neighbor[3][1]);
}


/*!
 ***********************************************************************
 * \brief
 *    decide the best position of each feature point
 *	  Added by Wonsang You at FEB 1, 2007
 *
 ***********************************************************************
 */

int decide_pos()
{
	int x,y,xn,yn,ii,jj,k;
	int point;
	int int_x, int_y, orig_x, orig_y, min;
	double sum;
	int pred_xx, pred_yy, orig_xx, orig_yy;
	int sum_x, sum_y;
	int search_hi;
	int intra_period, curr_frame_type;
	int max_search_point_nr;
	double center_distance, motion_parameter, square_vec_diff, reliability, min_energy;
	double local_form_energy,local_energy,global_min_energy;
	int curr_dist[2]={0,0};
	int prev_dist[2]={0,0};
	int selected_idx,h,min_search_idx,min_path_idx;
	int int_bx,int_by,orig_bx,orig_by;
	double global_min_texture_energy, global_min_motion_energy, global_min_form_energy;
	double weight_sum, neuron_value, multi, learning_const;
	double te_control_factor,me_control_factor,fe_control_factor;
	double squared_err;
	double avg_reliability;

	struct FeaturePoint **feature;
	double **local_min_energy;
	int **local_min_path;

	//for each feature point
	if(occlusion_handling == 1)
	{
		if(motion_picture_num>0)
		{
			//Calculation of AvgDiffPF
			
			//whether occlusion is finished
			if((occlusion_flag == FALSE)||((occlusion_flag == TRUE)&&(AvgDiffPF < th_mvd)))
			{
				intra_period = INTRA_PERIOD;
				curr_frame_type = motion_picture_num%intra_period;

				if(curr_frame_type==0)
					search_hi = search_half_interval + intra_search_err;
				else
					search_hi = search_half_interval;
				max_search_point_nr = (2*search_hi+1)*(2*search_hi+1);

				//Feature Buffers
				feature = (struct FeaturePoint**)calloc(init_con_nr,sizeof(struct FeaturePoint*));
				for(k=0; k<init_con_nr; k++)
					feature[k] = (struct FeaturePoint*)calloc(max_search_point_nr,sizeof(struct FeaturePoint));
				local_min_energy = (double**)calloc(init_con_nr,sizeof(double*));
				for(k=0; k<init_con_nr; k++)
					local_min_energy[k] = (double*)calloc(max_search_point_nr,sizeof(double));
				local_min_path = (int**)calloc(init_con_nr,sizeof(int*));
				for(k=0; k<init_con_nr; k++)
					local_min_path[k] = (int*)calloc(max_search_point_nr,sizeof(int));

				//Control Parameter
				te_control_factor = 1.0;
				fe_control_factor = 0.01;
				me_control_factor = 0.005;
				motion_parameter = 1.0;
				learning_const = 5;
				avg_reliability = 0;

				for(point=0; point<init_con_nr; point++)
				{
					//Fitting
					int_x = (int)(pred_con[point][0] +0.5);
					int_y = (int)(pred_con[point][1] +0.5);
					orig_x = (int)(upd_con[motion_picture_num-1][point][0] +0.5);
					orig_y = (int)(upd_con[motion_picture_num-1][point][1] +0.5);

					int_bx = (int)(pred_con[point][0]/4);
					int_by = (int)(pred_con[point][1]/4);
					orig_bx = (int)(upd_con[motion_picture_num-1][point][0]/4);
					orig_by = (int)(upd_con[motion_picture_num-1][point][1]/4);

					min = 100000000;

					//Check each search point
					square_vec_diff = (upd_mv[orig_by][orig_bx][0] + dec_mv[int_by][int_bx][0]/4)*(upd_mv[orig_by][orig_bx][0] + dec_mv[int_by][int_bx][0]/4)
									+ (upd_mv[orig_by][orig_bx][1] + dec_mv[int_by][int_bx][1]/4)*(upd_mv[orig_by][orig_bx][1] + dec_mv[int_by][int_bx][1]/4);
					reliability = exp(-(square_vec_diff/2/motion_parameter));
					avg_reliability += reliability;

					k=0;
					for(x=int_x-search_hi; x<=int_x+search_hi; x++)
					{
						for(y=int_y-search_hi; y<=int_y+search_hi; y++)
						{
							sum = 0;

							xn = (int)clipx((double)x);
							yn = (int)clipy((double)y);

							feature[point][k].x = xn;
							feature[point][k].y = yn;
													
							//Texture Energy
							for(ii=-sad_half_interval; ii<=sad_half_interval; ii++)
								for(jj=-sad_half_interval; jj<=sad_half_interval; jj++)
								{
									pred_xx = (int)clipx((double)(xn + ii));
									pred_yy = (int)clipy((double)(yn + jj));
									orig_xx = (int)clipx((double)(orig_x + ii));
									orig_yy = (int)clipy((double)(orig_y + jj));
									sum += ((curr_frm[pred_yy][pred_xx][0] - prev_frm[orig_yy][orig_xx][0])*(curr_frm[pred_yy][pred_xx][0] - prev_frm[orig_yy][orig_xx][0]));
									//sum += ((curr_frm[pred_yy][pred_xx][1] - prev_frm[orig_yy][orig_xx][1])*(curr_frm[pred_yy][pred_xx][1] - prev_frm[orig_yy][orig_xx][1]));
									//sum += ((curr_frm[pred_yy][pred_xx][2] - prev_frm[orig_yy][orig_xx][2])*(curr_frm[pred_yy][pred_xx][2] - prev_frm[orig_yy][orig_xx][2]));							
									//sum += fabs(curr_frm[pred_yy][pred_xx][0] - prev_frm[orig_yy][orig_xx][0]);
								}
							feature[point][k].te = te_control_factor*sum /255/255/(2*sad_half_interval+1)/(2*sad_half_interval+1);

							//Motion Energy
							//dist_parameter = 100;
							//bxn = (int)(xn/4);
							//byn = (int)(yn/4);
							//square_vec_diff = (upd_mv[orig_by][orig_bx][0] + dec_mv[byn][bxn][0])*(upd_mv[orig_by][orig_bx][0] + dec_mv[byn][bxn][0])
							//				+ (upd_mv[orig_by][orig_bx][1] + dec_mv[byn][bxn][1])*(upd_mv[orig_by][orig_bx][1] + dec_mv[byn][bxn][1]);
							//reliability = exp(-(square_vec_diff/2/motion_parameter));
							center_distance = sqrt((xn-int_x)*(xn-int_x) + (yn-int_y)*(yn-int_y));
							feature[point][k].me = me_control_factor*reliability*center_distance;
							//feature[point][k].me = me_control_factor*(1-reliability);
							//feature[point][k].me = (1 - reliability)*exp(center_distance/dist_parameter);

							//Form Energy and Dynamic Programming
							if(point==0)
							{
								feature[point][k].fe = 0;
								local_min_energy[point][k] = wt*feature[point][k].te + wm*feature[point][k].me;
							}
							else
							{
								min_energy = 100000000;
								for(h=0; h<max_search_point_nr; h++)
								{
									curr_dist[0] = xn - feature[point-1][h].x;
									curr_dist[1] = yn - feature[point-1][h].y;

									prev_dist[0] = upd_con[motion_picture_num-1][point][0] - upd_con[motion_picture_num-1][point-1][0];
									prev_dist[1] = upd_con[motion_picture_num-1][point][1] - upd_con[motion_picture_num-1][point-1][1];

									local_form_energy = sqrt((curr_dist[0]-prev_dist[0])*(curr_dist[0]-prev_dist[0])
													+ (curr_dist[1]-prev_dist[1])*(curr_dist[1]-prev_dist[1]));
									local_form_energy = fe_control_factor*sqrt(local_form_energy);
									//local_form_energy = local_form_energy/4/(prev_dist[0]*prev_dist[0]+prev_dist[1]*prev_dist[1]);
									local_energy = wf*local_form_energy + wt*feature[point][k].te + wm*feature[point][k].me
													+ local_min_energy[point-1][h];

									if(local_energy < min_energy)
									{
										min_energy = local_energy;
										feature[point][k].fe = local_form_energy;
										selected_idx = h;
									}
								}
								local_min_energy[point][k] = min_energy;
								local_min_path[point][k] = selected_idx;
							}

							k++;
						}
					}
				}

				//Dynamic Programming
				min_energy = 100000000;
				for(k=0; k<max_search_point_nr; k++)
				{
					if(local_min_energy[init_con_nr-1][k]<min_energy)
					{
						min_energy = local_min_energy[init_con_nr-1][k];
						min_search_idx = k;
					}
				}
				global_min_energy = min_energy;		//enable to delete this command
				global_min_texture_energy = feature[init_con_nr-1][min_search_idx].te;
				global_min_motion_energy = feature[init_con_nr-1][min_search_idx].me;
				global_min_form_energy = feature[init_con_nr-1][min_search_idx].fe;
				upd_con[motion_picture_num][init_con_nr-1][0] = feature[init_con_nr-1][min_search_idx].x;
				upd_con[motion_picture_num][init_con_nr-1][1] = feature[init_con_nr-1][min_search_idx].y;

				for(point=init_con_nr-1; point>0; point--)
				{
					min_path_idx = local_min_path[point][min_search_idx];
					upd_con[motion_picture_num][point-1][0] = feature[point-1][min_path_idx].x;
					upd_con[motion_picture_num][point-1][1] = feature[point-1][min_path_idx].y;
					global_min_texture_energy += feature[point-1][min_path_idx].te;
					global_min_motion_energy += feature[point-1][min_path_idx].me;
					global_min_form_energy += feature[point-1][min_path_idx].fe;

					min_search_idx = min_path_idx;			
				}


				//Update Weight Factors (Neural Network)
				weight_sum = wt*global_min_texture_energy + wm*global_min_motion_energy + wf*global_min_form_energy - w0;
				//neuron_value = (2/(1+exp(-weight_sum)))-1;
				//multi = -learning_const*neuron_value*(1-neuron_value)*(1-neuron_value)/2;
				weight_sum=weight_sum;
				neuron_value = 1/(1+exp(-weight_sum));
				multi = learning_const*(0.5-neuron_value)*neuron_value*(1-neuron_value);
				w0 += -multi;
				wt += global_min_texture_energy*multi;
				wf += global_min_form_energy*multi;
				wm += global_min_motion_energy*multi;
				squared_err = (0.5-neuron_value)*(0.5-neuron_value);

				//Calculation of AvgDiffPL

				//Occlusion Detection
				if(((global_min_form_energy>th_form)&&(AvgDiffPE > th_mvd))
					||((global_min_texture_energy>th_texture)&&(avg_reliability<th_reliability)))
				{
						occlusion_flag = TRUE;

						//3D parameter estimation

						//calculate pred_mv

						//prediction of last_mv from 3D parameter

						//calculate indices of feature points
				}


				//Report Results
				//avg_reliability = avg_reliability / init_con_nr;
				//fprintf(ET_FILE,"%f	%f	%f	%f	%f\n",global_min_texture_energy,global_min_form_energy,global_min_motion_energy,squared_err,avg_reliability);
				//fprintf(WT_FILE,"%f	%f	%f	%f	%f	%f\n",wt,wf,wm,global_min_texture_energy*multi,global_min_form_energy*multi,global_min_motion_energy*multi);

				free(feature);
				free_mem2Ddouble(local_min_energy);
				free_mem2Dint(local_min_path);
			}
			else
			{
				//3D parameter estimation

				//calculate pred_mv

				//prediction of last_mv from 3D parameter

				//calculate indices of feature points
			}
		}

		//Object Center
		sum_x = 0;
		sum_y = 0;
		for(point=0; point<init_con_nr; point++)
		{
			sum_x += (int)upd_con[motion_picture_num][point][0];
			sum_y += (int)upd_con[motion_picture_num][point][1];
		}
		object_center[motion_picture_num][0] = sum_x / init_con_nr;
		object_center[motion_picture_num][1] = sum_y / init_con_nr;
	}
	else
	{
		//for each feature point
		if(motion_picture_num>0)
		{
			intra_period = INTRA_PERIOD;
			curr_frame_type = motion_picture_num%intra_period;

			if(curr_frame_type==0)
				search_hi = search_half_interval + intra_search_err;
			else
				search_hi = search_half_interval;
			max_search_point_nr = (2*search_hi+1)*(2*search_hi+1);

			//Feature Buffers
			feature = (struct FeaturePoint**)calloc(init_con_nr,sizeof(struct FeaturePoint*));
			for(k=0; k<init_con_nr; k++)
				feature[k] = (struct FeaturePoint*)calloc(max_search_point_nr,sizeof(struct FeaturePoint));
			local_min_energy = (double**)calloc(init_con_nr,sizeof(double*));
			for(k=0; k<init_con_nr; k++)
				local_min_energy[k] = (double*)calloc(max_search_point_nr,sizeof(double));
			local_min_path = (int**)calloc(init_con_nr,sizeof(int*));
			for(k=0; k<init_con_nr; k++)
				local_min_path[k] = (int*)calloc(max_search_point_nr,sizeof(int));

			//Control Parameter
			te_control_factor = 1.0;
			fe_control_factor = 0.01;
			me_control_factor = 0.005;
			motion_parameter = 1.0;
			learning_const = 5;
			avg_reliability = 0;

			for(point=0; point<init_con_nr; point++)
			{
				//Fitting
				int_x = (int)(pred_con[point][0] +0.5);
				int_y = (int)(pred_con[point][1] +0.5);
				orig_x = (int)(upd_con[motion_picture_num-1][point][0] +0.5);
				orig_y = (int)(upd_con[motion_picture_num-1][point][1] +0.5);

				int_bx = (int)(pred_con[point][0]/4);
				int_by = (int)(pred_con[point][1]/4);
				orig_bx = (int)(upd_con[motion_picture_num-1][point][0]/4);
				orig_by = (int)(upd_con[motion_picture_num-1][point][1]/4);

				min = 100000000;

				//Check each search point
				square_vec_diff = (upd_mv[orig_by][orig_bx][0] + dec_mv[int_by][int_bx][0]/4)*(upd_mv[orig_by][orig_bx][0] + dec_mv[int_by][int_bx][0]/4)
								+ (upd_mv[orig_by][orig_bx][1] + dec_mv[int_by][int_bx][1]/4)*(upd_mv[orig_by][orig_bx][1] + dec_mv[int_by][int_bx][1]/4);
				reliability = exp(-(square_vec_diff/2/motion_parameter));
				avg_reliability += reliability;

				k=0;
				for(x=int_x-search_hi; x<=int_x+search_hi; x++)
				{
					for(y=int_y-search_hi; y<=int_y+search_hi; y++)
					{
						sum = 0;

						xn = (int)clipx((double)x);
						yn = (int)clipy((double)y);

						feature[point][k].x = xn;
						feature[point][k].y = yn;
												
						//Texture Energy
						for(ii=-sad_half_interval; ii<=sad_half_interval; ii++)
							for(jj=-sad_half_interval; jj<=sad_half_interval; jj++)
							{
								pred_xx = (int)clipx((double)(xn + ii));
								pred_yy = (int)clipy((double)(yn + jj));
								orig_xx = (int)clipx((double)(orig_x + ii));
								orig_yy = (int)clipy((double)(orig_y + jj));
								sum += ((curr_frm[pred_yy][pred_xx][0] - prev_frm[orig_yy][orig_xx][0])*(curr_frm[pred_yy][pred_xx][0] - prev_frm[orig_yy][orig_xx][0]));
								//sum += ((curr_frm[pred_yy][pred_xx][1] - prev_frm[orig_yy][orig_xx][1])*(curr_frm[pred_yy][pred_xx][1] - prev_frm[orig_yy][orig_xx][1]));
								//sum += ((curr_frm[pred_yy][pred_xx][2] - prev_frm[orig_yy][orig_xx][2])*(curr_frm[pred_yy][pred_xx][2] - prev_frm[orig_yy][orig_xx][2]));							
								//sum += fabs(curr_frm[pred_yy][pred_xx][0] - prev_frm[orig_yy][orig_xx][0]);
							}
						feature[point][k].te = te_control_factor*sum /255/255/(2*sad_half_interval+1)/(2*sad_half_interval+1);

						//Motion Energy
						//dist_parameter = 100;
						//bxn = (int)(xn/4);
						//byn = (int)(yn/4);
						//square_vec_diff = (upd_mv[orig_by][orig_bx][0] + dec_mv[byn][bxn][0])*(upd_mv[orig_by][orig_bx][0] + dec_mv[byn][bxn][0])
						//				+ (upd_mv[orig_by][orig_bx][1] + dec_mv[byn][bxn][1])*(upd_mv[orig_by][orig_bx][1] + dec_mv[byn][bxn][1]);
						//reliability = exp(-(square_vec_diff/2/motion_parameter));
						center_distance = sqrt((xn-int_x)*(xn-int_x) + (yn-int_y)*(yn-int_y));
						feature[point][k].me = me_control_factor*reliability*center_distance;
						//feature[point][k].me = me_control_factor*(1-reliability);
						//feature[point][k].me = (1 - reliability)*exp(center_distance/dist_parameter);

						//Form Energy and Dynamic Programming
						if(point==0)
						{
							feature[point][k].fe = 0;
							local_min_energy[point][k] = wt*feature[point][k].te + wm*feature[point][k].me;
						}
						else
						{
							min_energy = 100000000;
							for(h=0; h<max_search_point_nr; h++)
							{
								curr_dist[0] = xn - feature[point-1][h].x;
								curr_dist[1] = yn - feature[point-1][h].y;

								prev_dist[0] = upd_con[motion_picture_num-1][point][0] - upd_con[motion_picture_num-1][point-1][0];
								prev_dist[1] = upd_con[motion_picture_num-1][point][1] - upd_con[motion_picture_num-1][point-1][1];

								local_form_energy = sqrt((curr_dist[0]-prev_dist[0])*(curr_dist[0]-prev_dist[0])
												+ (curr_dist[1]-prev_dist[1])*(curr_dist[1]-prev_dist[1]));
								local_form_energy = fe_control_factor*sqrt(local_form_energy);
								//local_form_energy = local_form_energy/4/(prev_dist[0]*prev_dist[0]+prev_dist[1]*prev_dist[1]);
								local_energy = wf*local_form_energy + wt*feature[point][k].te + wm*feature[point][k].me
												+ local_min_energy[point-1][h];

								if(local_energy < min_energy)
								{
									min_energy = local_energy;
									feature[point][k].fe = local_form_energy;
									selected_idx = h;
								}
							}
							local_min_energy[point][k] = min_energy;
							local_min_path[point][k] = selected_idx;
						}

						k++;
					}
				}
			}

			//Dynamic Programming
			min_energy = 100000000;
			for(k=0; k<max_search_point_nr; k++)
			{
				if(local_min_energy[init_con_nr-1][k]<min_energy)
				{
					min_energy = local_min_energy[init_con_nr-1][k];
					min_search_idx = k;
				}
			}
			global_min_energy = min_energy;		//enable to delete this command
			global_min_texture_energy = feature[init_con_nr-1][min_search_idx].te;
			global_min_motion_energy = feature[init_con_nr-1][min_search_idx].me;
			global_min_form_energy = feature[init_con_nr-1][min_search_idx].fe;
			upd_con[motion_picture_num][init_con_nr-1][0] = feature[init_con_nr-1][min_search_idx].x;
			upd_con[motion_picture_num][init_con_nr-1][1] = feature[init_con_nr-1][min_search_idx].y;

			for(point=init_con_nr-1; point>0; point--)
			{
				min_path_idx = local_min_path[point][min_search_idx];
				upd_con[motion_picture_num][point-1][0] = feature[point-1][min_path_idx].x;
				upd_con[motion_picture_num][point-1][1] = feature[point-1][min_path_idx].y;
				global_min_texture_energy += feature[point-1][min_path_idx].te;
				global_min_motion_energy += feature[point-1][min_path_idx].me;
				global_min_form_energy += feature[point-1][min_path_idx].fe;

				min_search_idx = min_path_idx;			
			}


			//Update Weight Factors (Neural Network)
			weight_sum = wt*global_min_texture_energy + wm*global_min_motion_energy + wf*global_min_form_energy - w0;
			//neuron_value = (2/(1+exp(-weight_sum)))-1;
			//multi = -learning_const*neuron_value*(1-neuron_value)*(1-neuron_value)/2;
			weight_sum=weight_sum;
			neuron_value = 1/(1+exp(-weight_sum));
			multi = learning_const*(0.5-neuron_value)*neuron_value*(1-neuron_value);
			w0 += -multi;
			wt += global_min_texture_energy*multi;
			wf += global_min_form_energy*multi;
			wm += global_min_motion_energy*multi;
			squared_err = (0.5-neuron_value)*(0.5-neuron_value);


			//Report Results
			if(analysis_mode==1)
			{
				if(modify_tracking_flag==0)
				{
					avg_reliability = avg_reliability / init_con_nr;
					fprintf(ET_FILE,"%f	%f	%f	%f	%f\n",global_min_texture_energy,global_min_form_energy,global_min_motion_energy,squared_err,avg_reliability);
					fprintf(WT_FILE,"%f	%f	%f	%f	%f	%f\n",wt,wf,wm,global_min_texture_energy*multi,global_min_form_energy*multi,global_min_motion_energy*multi);
				}
			}

			free(feature);
			free_mem2Ddouble(local_min_energy);
			free_mem2Dint(local_min_path);
		}

		//Object Center
		sum_x = 0;
		sum_y = 0;
		for(point=0; point<init_con_nr; point++)
		{
			sum_x += (int)upd_con[motion_picture_num][point][0];
			sum_y += (int)upd_con[motion_picture_num][point][1];
		}
		object_center[motion_picture_num][0] = sum_x / init_con_nr;
		object_center[motion_picture_num][1] = sum_y / init_con_nr;
	}
	

	return 1;
}


/*!
 ***********************************************************************
 * \brief
 *    set buffers
 *	  Added by Wonsang You at FEB 1, 2007
 *
 ***********************************************************************
 */

void set_buffer_PD()
{
	int i,j,bx,by;
	int b_height = (int)iheight / 4;
	int b_width = (int)iwidth / 4;
	int max_mb_nr = b_height*b_width;
	int tol = 20;
	int YSize = iheight*iwidth;
	int Voffset = YSize + YSize/4;
	int FrameSize = YSize*1.5;

	//save the current frame as the previous frame
	for(i=0; i<iheight; i++)
	{
		for(j=0; j<iwidth; j++)
		{
			prev_frm[i][j][0] = curr_frm[i][j][0];
			prev_frm[i][j][1] = curr_frm[i][j][1];
			prev_frm[i][j][2] = curr_frm[i][j][2];
		}
	}

	//save partial decoded pictures into display buffer
	if(video_save_flag==1)
	{
		if(motion_picture_num>0)
		{			
			/*
			// Save the background subtracted picture
			for(i=0; i<iheight; i++)
			{
				for(j=0; j<iwidth; j++)
				{
					if(IsNearlyZero(curr_frm[i][j][0] - BkdBuf[0][i*iwidth+j],tol))
						bkd_str[0][i*iwidth+j] = 16; //clip(DecBuf[j] - BkdBuf[j] +16);
					else
						bkd_str[0][i*iwidth+j] = curr_frm[i][j][0];

					if(IsNearlyZero(curr_frm[i>>1][j>>1][1] - BkdBuf[0][(i>>1)*(iwidth>>1)+(j>>1)+YSize],tol))
						bkd_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = 128; //clip(DecBuf[j+YSize] - BkdBuf[j+YSize] +128);
					else
						bkd_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][1];

					if(IsNearlyZero(curr_frm[i>>1][j>>1][2] - BkdBuf[0][(i>>1)*(iwidth>>1)+(j>>1)+Voffset],tol))
						bkd_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = 128; //clip(DecBuf[j+Voffset] - BkdBuf[j+Voffset] +128);
					else
						bkd_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][2];
				}
			}

			for(i = 0; i < iheight ; i++)
				fwrite(&(bkd_str[0][i*iwidth]),1,iwidth,BKS_FILE );				
			for(i = 0; i < (iheight/2); i++)
				fwrite(&(bkd_str[1][i*iwidth/2]),1,iwidth/2,BKS_FILE );
			for(i = 0; i < (iheight/2); i++)
				fwrite(&(bkd_str[2][i*iwidth/2]),1,iwidth/2,BKS_FILE );
			*/

			if(extract_flag==1)
			{
				for(i=0; i<iheight; i++)
				{
					for(j=0; j<iwidth; j++)
					{
						if(residual_reverse_contrast_flag==1)
						{
							res_str[0][i*iwidth+j] = max(0,min(255,ContrastStretch(10,curr_res[i][j][0]) +16));
							res_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = 128;
							res_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = 128;
						}
						else
						{
							res_str[0][i*iwidth+j] = max(0,min(255,curr_res[i][j][0] +16));
							res_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = max(0,min(255,curr_res[i>>1][j>>1][1] +128));
							res_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = max(0,min(255,curr_res[i>>1][j>>1][2] +128));
						}
					}
				}

				for(i = 0; i < iheight ; i++)
					fwrite(&(res_str[0][i*iwidth]),1,iwidth,RES_FILE );				
				for(i = 0; i < (iheight/2); i++)
					fwrite(&(res_str[1][i*iwidth/2]),1,iwidth/2,RES_FILE );
				for(i = 0; i < (iheight/2); i++)
					fwrite(&(res_str[2][i*iwidth/2]),1,iwidth/2,RES_FILE );
			}
		}
	}

	//initialize the current frame buffer
	for(i=0; i<iheight; i++)
	{
		for(j=0; j<iwidth; j++)
		{
			curr_frm[i][j][0] = 0;
			curr_frm[i][j][1] = 0;
			curr_frm[i][j][2] = 0;
		}
	}

	//initialize the current residual buffer
	if(video_save_flag==1)
	{
		for(i=0; i<iheight; i++)
		{
			for(j=0; j<iwidth; j++)
			{
				curr_res[i][j][0] = 0;
				curr_res[i][j][1] = 0;
				curr_res[i][j][2] = 0;
			}
		}
	}

	//initialize decoded block flag
	if(tracking_flag==1)
	{
		for(i=0; i<max_mb_nr; i++)
		{
			blk_dec_enb[i] = 0;
			for(j=0; j<16; j++)
				blk_list[i][j] = 0;
		}
	}

	//Initialize the updated forward motion vectors
	if((extract_flag==1)||(tracking_flag==1))
	{
		for(by=0; by<b_height; by++)
		{
			for(bx=0; bx<b_width; bx++)
			{
				dec_mv[by][bx][0] = 0;
				dec_mv[by][bx][1] = 0;
			}
		}
	}

	//initialize the current frame buffer
	if((multiple_block_drawing_flag==1)&&(time_measure_flag==0))
	{
		for(i=0; i<10; i++)
			for(j=0; j<b_height*b_width; j++)
			{
				DrawBlockFlag[i][j] = 0;
			}
	}
}


/*!
 ************************************************************************
 * \brief
 *    Convert file read buffer to source picture structure
 * \param imgX
 *    Pointer to image plane
 * \param buf
 *    Buffer for file output
 * \param size_x
 *    horizontal image size in pixel
 * \param size_y
 *    vertical image size in pixel
 * \param symbol_size_in_bytes
 *    number of bytes used per pel
 ************************************************************************
 */
void buf2img (imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes)
{
  int i,j;

  unsigned short tmp16, ui16;
  unsigned long  tmp32, ui32;

  if (symbol_size_in_bytes> sizeof(imgpel))
  {
    error ("Source picture has higher bit depth than imgpel data type. Please recompile with larger data type for imgpel.", 500);
  }

  if (( sizeof(char) == sizeof (imgpel)) && ( sizeof(char) == symbol_size_in_bytes))
  {
    // imgpel == pixel_in_file == 1 byte -> simple copy
    for(j=0;j<size_y;j++)
      memcpy(imgX[j], buf+j*size_x, size_x);
  }
  else
  {
    // sizeof (imgpel) > sizeof(char)
    if (testEndian())
    {
      // big endian
      switch (symbol_size_in_bytes)
      {
      case 1:
        {
          for(j=0;j<size_y;j++)
            for(i=0;i<size_x;i++)
            {
              imgX[j][i]= buf[i+j*size_x];
            }
          break;
        }
      case 2:
        {
          for(j=0;j<size_y;j++)
            for(i=0;i<size_x;i++)
            {
              memcpy(&tmp16, buf+((i+j*size_x)*2), 2);
              ui16  = (tmp16 >> 8) | ((tmp16&0xFF)<<8);
              imgX[j][i] = (imgpel) ui16;
            }
          break;
        }
      case 4:
        {
          for(j=0;j<size_y;j++)
            for(i=0;i<size_x;i++)
            {
              memcpy(&tmp32, buf+((i+j*size_x)*4), 4);
              ui32  = ((tmp32&0xFF00)<<8) | ((tmp32&0xFF)<<24) | ((tmp32&0xFF0000)>>8) | ((tmp32&0xFF000000)>>24);
              imgX[j][i] = (imgpel) ui32;
            }
        }
      default:
        {
           error ("reading only from formats of 8, 16 or 32 bit allowed on big endian architecture", 500);
           break;
        }
      }

    }
    else
    {
      // little endian
      for (j=0; j < size_y; j++)
        for (i=0; i < size_x; i++)
        {
          imgX[j][i]=0;
          memcpy(&(imgX[j][i]), buf +((i+j*size_x)*symbol_size_in_bytes), symbol_size_in_bytes);
        }

    }
  }
}


/*!
************************************************************************
* \brief
*    Find PSNR for all three components.Compare decoded frame with
*    the original sequence. Read inp->jumpd frames to reflect frame skipping.
************************************************************************
*/
void find_snr(
              struct snr_par  *snr,   //!< pointer to snr parameters
              StorablePicture *p,     //!< picture to be compared
              int p_ref)              //!< open reference YUV file
{
  int SubWidthC  [4]= { 1, 2, 2, 1};
  int SubHeightC [4]= { 1, 2, 1, 1};
  int crop_left, crop_right, crop_top, crop_bottom;

  int i,j;
  int64 diff_y,diff_u,diff_v;
  int uv;
  int64  status;
  int symbol_size_in_bytes = img->pic_unit_bitsize_on_disk/8;
  int size_x, size_y;
  int size_x_cr, size_y_cr;
  int64 framesize_in_bytes;
  unsigned int max_pix_value_sqd = img->max_imgpel_value * img->max_imgpel_value;
  unsigned int max_pix_value_sqd_uv = img->max_imgpel_value_uv * img->max_imgpel_value_uv;
  Boolean rgb_output = (active_sps->vui_seq_parameters.matrix_coefficients==0);
  unsigned char *buf;

  // picture error concealment
  char yuv_types[4][6]= {"4:0:0","4:2:0","4:2:2","4:4:4"};

  // calculate frame number
  int  psnrPOC = active_sps->mb_adaptive_frame_field_flag ? p->poc /(input->poc_scale) : p->poc/(input->poc_scale);

  // cropping for luma
  if (p->frame_cropping_flag)
  {
    crop_left   = SubWidthC[p->chroma_format_idc] * p->frame_cropping_rect_left_offset;
    crop_right  = SubWidthC[p->chroma_format_idc] * p->frame_cropping_rect_right_offset;
    crop_top    = SubHeightC[p->chroma_format_idc]*( 2 - p->frame_mbs_only_flag ) *  p->frame_cropping_rect_top_offset;
    crop_bottom = SubHeightC[p->chroma_format_idc]*( 2 - p->frame_mbs_only_flag ) *   p->frame_cropping_rect_bottom_offset;
  }
  else
  {
    crop_left = crop_right = crop_top = crop_bottom = 0;
  }

  size_x = p->size_x - crop_left - crop_right;
  size_y = p->size_y - crop_top - crop_bottom;

  // cropping for chroma
  if (p->frame_cropping_flag)
  {
    crop_left   = p->frame_cropping_rect_left_offset;
    crop_right  = p->frame_cropping_rect_right_offset;
    crop_top    = ( 2 - p->frame_mbs_only_flag ) *  p->frame_cropping_rect_top_offset;
    crop_bottom = ( 2 - p->frame_mbs_only_flag ) *   p->frame_cropping_rect_bottom_offset;
  }
  else
  {
    crop_left = crop_right = crop_top = crop_bottom = 0;
  }

  if ((p->chroma_format_idc==YUV400) && input->write_uv)
  {
    size_x_cr = p->size_x/2;
    size_y_cr = p->size_y/2;
  }
  else
  {
    size_x_cr = p->size_x_cr - crop_left - crop_right;
    size_y_cr = p->size_y_cr - crop_top  - crop_bottom;
  }

  framesize_in_bytes = (((int64)size_y*size_x) + ((int64)size_y_cr*size_x_cr)*2) * symbol_size_in_bytes;

  if (psnrPOC==0 && img->psnr_number)
    img->idr_psnr_number = img->number*img->ref_poc_gap/(input->poc_scale);

  img->psnr_number=max(img->psnr_number,img->idr_psnr_number+psnrPOC);

  frame_no = img->idr_psnr_number+psnrPOC;

  // KS: this buffer should actually be allocated only once, but this is still much faster than the previous version
  buf = malloc ( size_y * size_x * symbol_size_in_bytes );

  if (NULL == buf)
  {
    no_mem_exit("find_snr: buf");
  }

  status = lseek (p_ref, framesize_in_bytes * frame_no, SEEK_SET);
  if (status == -1)
  {
    fprintf(stderr, "Error in seeking frame number: %d\n", frame_no);
    free (buf);
    return;
  }

  if(rgb_output)
    lseek (p_ref, framesize_in_bytes/3, SEEK_CUR);

  read(p_ref, buf, size_y * size_x * symbol_size_in_bytes);
  buf2img(imgY_ref, buf, size_x, size_y, symbol_size_in_bytes);

  if (p->chroma_format_idc != YUV400)
  {
    for (uv=0; uv < 2; uv++)
    {
      if(rgb_output && uv==1)
        lseek (p_ref, -framesize_in_bytes, SEEK_CUR);
    
      read(p_ref, buf, size_y_cr * size_x_cr*symbol_size_in_bytes);
      buf2img(imgUV_ref[uv], buf, size_x_cr, size_y_cr, symbol_size_in_bytes);
    }
  }

   if(rgb_output) 
     lseek (p_ref, framesize_in_bytes*2/3, SEEK_CUR);
  
  free (buf);

  img->quad[0]=0;
  diff_y=0;
  for (j=0; j < size_y; ++j)
  {
    for (i=0; i < size_x; ++i)
    {
      diff_y += img->quad[abs(p->imgY[j][i]-imgY_ref[j][i])];
    }
  }

  // Chroma
  diff_u=0;
  diff_v=0;
  
  if (p->chroma_format_idc != YUV400)
  {
    for (j=0; j < size_y_cr; ++j)
    {
      for (i=0; i < size_x_cr; ++i)
      {
        diff_u += img->quad[abs(imgUV_ref[0][j][i]-p->imgUV[0][j][i])];
        diff_v += img->quad[abs(imgUV_ref[1][j][i]-p->imgUV[1][j][i])];
      }
    }
  }

#if ZEROSNR
  if (diff_y == 0)
    diff_y = 1;
  if (diff_u == 0)
    diff_u = 1;
  if (diff_v == 0)
    diff_v = 1; 
#endif

  // Collecting SNR statistics
  if (diff_y != 0)
    snr->snr_y=(float)(10*log10(max_pix_value_sqd*(double)((double)(size_x)*(size_y) / diff_y)));        // luma snr for current frame
  else
    snr->snr_y=0.0;
  if (diff_u != 0)
    snr->snr_u=(float)(10*log10(max_pix_value_sqd_uv*(double)((double)(size_x_cr)*(size_y_cr) / (diff_u))));    //  chroma snr for current frame
  else
    snr->snr_u=0.0;
  if (diff_v != 0)
    snr->snr_v=(float)(10*log10(max_pix_value_sqd_uv*(double)((double)(size_x_cr)*(size_y_cr) / (diff_v))));    //  chroma snr for current frame
  else
    snr->snr_v=0;

  if (img->number == 0) // first
  {
    snr->snr_ya=snr->snr_y1=snr->snr_y;                                                        // keep luma snr for first frame
    snr->snr_ua=snr->snr_u1=snr->snr_u;                                                        // keep chroma snr for first frame
    snr->snr_va=snr->snr_v1=snr->snr_v;                                                        // keep chroma snr for first frame
  
  }
  else
  {
    snr->snr_ya=(float)(snr->snr_ya*(snr->frame_ctr)+snr->snr_y)/(snr->frame_ctr+1); // average snr chroma for all frames
    snr->snr_ua=(float)(snr->snr_ua*(snr->frame_ctr)+snr->snr_u)/(snr->frame_ctr+1); // average snr luma for all frames
    snr->snr_va=(float)(snr->snr_va*(snr->frame_ctr)+snr->snr_v)/(snr->frame_ctr+1); // average snr luma for all frames
  } 

  // picture error concealment
  if(p->concealed_pic)
  {
      fprintf(stdout,"%04d(P)  %8d %5d %5d %7.4f %7.4f %7.4f  %s %5d\n", 
          frame_no, p->frame_poc, p->pic_num, p->qp, 
          snr->snr_y, snr->snr_u, snr->snr_v, yuv_types[p->chroma_format_idc], 0);      

  }
}


/*!
 ************************************************************************
 * \brief
 *    Interpolation of 1/4 subpixel
 ************************************************************************
 */
void get_block(int ref_frame, StorablePicture **list, int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE])
{

  int dx, dy;
  int x, y;
  int i, j;
  int maxold_x,maxold_y;
  int result;
  int pres_x;
  int pres_y; 
  int tmp_res[4][9];
  static const int COEF[6] = {    1, -5, 20, 20, -5, 1  };

  dx = x_pos&3;
  dy = y_pos&3;
  x_pos = (x_pos-dx)/4;
  y_pos = (y_pos-dy)/4;

  maxold_x = dec_picture->size_x-1;
  maxold_y = dec_picture->size_y-1;

  if (dec_picture->mb_field[img->current_mb_nr])
    maxold_y = dec_picture->size_y/2 - 1;

  if (dx == 0 && dy == 0) 
  {  /* fullpel position */
    for (j = 0; j < BLOCK_SIZE; j++)
      for (i = 0; i < BLOCK_SIZE; i++)
        block[i][j] = list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i))];
  }
  else 
  { /* other positions */

    if (dy == 0) 
    { /* No vertical interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) 
      {
        for (i = 0; i < BLOCK_SIZE; i++) 
        {
          for (result = 0, x = -2; x < 4; x++)
            result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
          block[i][j] = max(0, min(img->max_imgpel_value, (result+16)/32));
        }
      }

      if ((dx&1) == 1) 
      {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block[i][j] = (block[i][j] + list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+dx/2))] +1 )/2;
      }
    }
    else if (dx == 0) 
    {  /* No horizontal interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) 
      {
        for (i = 0; i < BLOCK_SIZE; i++) 
        {
          for (result = 0, y = -2; y < 4; y++)
            result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
          block[i][j] = max(0, min(img->max_imgpel_value, (result+16)/32));
        }
      }

      if ((dy&1) == 1) 
      {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
           block[i][j] = (block[i][j] + list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+dy/2))][max(0,min(maxold_x,x_pos+i))] +1 )/2;
      }
    }
    else if (dx == 2) 
    {  /* Vertical & horizontal interpolation */

      for (j = -2; j < BLOCK_SIZE+3; j++) 
      {
        for (i = 0; i < BLOCK_SIZE; i++)
          for (tmp_res[i][j+2] = 0, x = -2; x < 4; x++)
            tmp_res[i][j+2] += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
      }

      for (j = 0; j < BLOCK_SIZE; j++) 
      {
        for (i = 0; i < BLOCK_SIZE; i++) 
        {
          for (result = 0, y = -2; y < 4; y++)
            result += tmp_res[i][j+y+2]*COEF[y+2];
          block[i][j] = max(0, min(img->max_imgpel_value, (result+512)/1024));
        } 
      }

      if ((dy&1) == 1)
      {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block[i][j] = (block[i][j] + max(0, min(img->max_imgpel_value, (tmp_res[i][j+2+dy/2]+16)/32)) +1 )/2;
      }
    }
    else if (dy == 2)
    {  /* Horizontal & vertical interpolation */

      for (j = 0; j < BLOCK_SIZE; j++)
      {
        for (i = -2; i < BLOCK_SIZE+3; i++)
          for (tmp_res[j][i+2] = 0, y = -2; y < 4; y++)
            tmp_res[j][i+2] += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
      }

      for (j = 0; j < BLOCK_SIZE; j++)
      {
        for (i = 0; i < BLOCK_SIZE; i++)
        {
          for (result = 0, x = -2; x < 4; x++)
            result += tmp_res[j][i+x+2]*COEF[x+2];
          block[i][j] = max(0, min(img->max_imgpel_value, (result+512)/1024));
        }
      }

      if ((dx&1) == 1)
      {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block[i][j] = (block[i][j] + max(0, min(img->max_imgpel_value, (tmp_res[j][i+2+dx/2]+16)/32))+1)/2;
      }
    }
    else
    {  /* Diagonal interpolation */

      for (j = 0; j < BLOCK_SIZE; j++)
      {
        for (i = 0; i < BLOCK_SIZE; i++)
        {
          pres_y = dy == 1 ? y_pos+j : y_pos+j+1;
          pres_y = max(0,min(maxold_y,pres_y));
          for (result = 0, x = -2; x < 4; x++)
            result += list[ref_frame]->imgY[pres_y][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
          block[i][j] = max(0, min(img->max_imgpel_value, (result+16)/32));
        }
      }

      for (j = 0; j < BLOCK_SIZE; j++)
      {
        for (i = 0; i < BLOCK_SIZE; i++)
        {
          pres_x = dx == 1 ? x_pos+i : x_pos+i+1;
          pres_x = max(0,min(maxold_x,pres_x));
          for (result = 0, y = -2; y < 4; y++)
            result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][pres_x]*COEF[y+2];
          block[i][j] = (block[i][j] + max(0, min(img->max_imgpel_value, (result+16)/32)) +1 ) / 2;
        }
      }

    }
  }
}


void reorder_lists(int currSliceType, Slice * currSlice)
{

  if ((currSliceType != I_SLICE)&&(currSliceType != SI_SLICE))
  {
    if (currSlice->ref_pic_list_reordering_flag_l0)
    {
      reorder_ref_pic_list(listX[0], &listXsize[0], 
                           img->num_ref_idx_l0_active - 1, 
                           currSlice->reordering_of_pic_nums_idc_l0, 
                           currSlice->abs_diff_pic_num_minus1_l0, 
                           currSlice->long_term_pic_idx_l0);
    }
    if (NULL == listX[0][img->num_ref_idx_l0_active-1])
    {
      error("RefPicList0[ num_ref_idx_l0_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
    }
    // that's a definition
    listXsize[0] = img->num_ref_idx_l0_active;
  }
  if (currSliceType == B_SLICE)
  {
    if (currSlice->ref_pic_list_reordering_flag_l1)
    {
      reorder_ref_pic_list(listX[1], &listXsize[1], 
                           img->num_ref_idx_l1_active - 1, 
                           currSlice->reordering_of_pic_nums_idc_l1, 
                           currSlice->abs_diff_pic_num_minus1_l1, 
                           currSlice->long_term_pic_idx_l1);
    }
    if (NULL == listX[1][img->num_ref_idx_l1_active-1])
    {
      error("RefPicList1[ num_ref_idx_l1_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
    }
    // that's a definition
    listXsize[1] = img->num_ref_idx_l1_active;
  }

  free_ref_pic_list_reordering_buffer(currSlice);
}


/*!
 ************************************************************************
 * \brief
 *    initialize ref_pic_num array
 ************************************************************************
 */
void set_ref_pic_num()
{
  int i,j;

  int slice_id=img->current_slice_nr;

  for (i=0;i<listXsize[LIST_0];i++)
  {
    dec_picture->ref_pic_num        [slice_id][LIST_0][i]=listX[LIST_0][i]->poc * 2 + ((listX[LIST_0][i]->structure==BOTTOM_FIELD)?1:0) ; 
    dec_picture->frm_ref_pic_num    [slice_id][LIST_0][i]=listX[LIST_0][i]->frame_poc * 2; 
    dec_picture->top_ref_pic_num    [slice_id][LIST_0][i]=listX[LIST_0][i]->top_poc * 2; 
    dec_picture->bottom_ref_pic_num [slice_id][LIST_0][i]=listX[LIST_0][i]->bottom_poc * 2 + 1; 
    //printf("POCS %d %d %d %d ",listX[LIST_0][i]->frame_poc,listX[LIST_0][i]->bottom_poc,listX[LIST_0][i]->top_poc,listX[LIST_0][i]->poc);
    //printf("refid %d %d %d %d\n",(int) dec_picture->frm_ref_pic_num[LIST_0][i],(int) dec_picture->top_ref_pic_num[LIST_0][i],(int) dec_picture->bottom_ref_pic_num[LIST_0][i],(int) dec_picture->ref_pic_num[LIST_0][i]);
  }

  for (i=0;i<listXsize[LIST_1];i++)
  {
    dec_picture->ref_pic_num        [slice_id][LIST_1][i]=listX[LIST_1][i]->poc  *2 + ((listX[LIST_1][i]->structure==BOTTOM_FIELD)?1:0);
    dec_picture->frm_ref_pic_num    [slice_id][LIST_1][i]=listX[LIST_1][i]->frame_poc * 2; 
    dec_picture->top_ref_pic_num    [slice_id][LIST_1][i]=listX[LIST_1][i]->top_poc * 2; 
    dec_picture->bottom_ref_pic_num [slice_id][LIST_1][i]=listX[LIST_1][i]->bottom_poc * 2 + 1; 
  }

  if (!active_sps->frame_mbs_only_flag)
  {
    if (img->structure==FRAME)
      for (j=2;j<6;j++)
        for (i=0;i<listXsize[j];i++)
        {
          dec_picture->ref_pic_num        [slice_id][j][i] = listX[j][i]->poc * 2 + ((listX[j][i]->structure==BOTTOM_FIELD)?1:0);
          dec_picture->frm_ref_pic_num    [slice_id][j][i] = listX[j][i]->frame_poc * 2 ;
          dec_picture->top_ref_pic_num    [slice_id][j][i] = listX[j][i]->top_poc * 2 ;
          dec_picture->bottom_ref_pic_num [slice_id][j][i] = listX[j][i]->bottom_poc * 2 + 1;
        }
  }

}


/*!
 ************************************************************************
 * \brief
 *    Reads new slice from bit_stream
 ************************************************************************
 */
int read_new_slice()
{
  NALU_t *nalu = AllocNALU(MAX_CODED_FRAME_SIZE);
  int current_header = 0;
  int ret;
  int BitsUsedByHeader;
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream;

  int slice_id_a, slice_id_b, slice_id_c;
  int redundant_pic_cnt_b, redundant_pic_cnt_c;
  long ftell_position, expected_slice_type;
  
//  int i;
  expected_slice_type = NALU_TYPE_DPA;

  while (1)
  {
	ftell_position = ftell(bits);

    if (input->FileFormat == PAR_OF_ANNEXB)
      ret=GetAnnexbNALU (nalu);
//    else
//      ret=GetRTPNALU (nalu);

    //In some cases, zero_byte shall be present. If current NALU is a VCL NALU, we can't tell
    //whether it is the first VCL NALU at this point, so only non-VCL NAL unit is checked here.
    CheckZeroByteNonVCL(nalu, &ret);

    NALUtoRBSP(nalu);
//    printf ("nalu->len %d\n", nalu->len);
    
    if (ret < 0)
      printf ("Error while getting the NALU in file format %s, exit\n", input->FileFormat==PAR_OF_ANNEXB?"Annex B":"RTP");
    if (ret == 0)
    {
//      printf ("read_new_slice: returning %s\n", "EOS");
      if(expected_slice_type != NALU_TYPE_DPA)
      {
        /* oops... we found the next slice, go back! */
        fseek(bits, ftell_position, SEEK_SET);
        FreeNALU(nalu);
        return current_header;
      }
      else
        return EOS;
    }

    // Got a NALU
    if (nalu->forbidden_bit)
    {
      printf ("Found NALU w/ forbidden_bit set, bit error?  Let's try...\n");
    }

    switch (nalu->nal_unit_type)
    {
      case NALU_TYPE_SLICE:
      case NALU_TYPE_IDR:
        img->idr_flag = (nalu->nal_unit_type == NALU_TYPE_IDR);
        img->nal_reference_idc = nalu->nal_reference_idc;
        img->disposable_flag = (nalu->nal_reference_idc == NALU_PRIORITY_DISPOSABLE);
        currSlice->dp_mode = PAR_DP_1;
        currSlice->max_part_nr = 1;
        currSlice->ei_flag = 0;
        currStream = currSlice->partArr[0].bitstream;
        currStream->ei_flag = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        // Some syntax of the Slice Header depends on the parameter set, which depends on
        // the parameter set ID of the SLice header.  Hence, read the pic_parameter_set_id
        // of the slice header first, then setup the active parameter sets, and then read
        // the rest of the slice header
        BitsUsedByHeader = FirstPartOfSliceHeader();
        UseParameterSet (currSlice->pic_parameter_set_id);
        BitsUsedByHeader+= RestOfSliceHeader ();

        FmoInit (active_pps, active_sps);

        AssignQuantParam (active_pps, active_sps);

        if(is_new_picture())
        {
		  //Added by Wonsang You at Oct-25-2006
		  //if(motion_picture_num>0)
			//store_motion_into_buffer();
		  motion_picture_num++;
		  //Added-End

          init_picture(img, input);
          
          current_header = SOP;
          //check zero_byte if it is also the first NAL unit in the access unit
          CheckZeroByteVCL(nalu, &ret);
        }
        else
          current_header = SOS;
  
        init_lists(img->type, img->currentSlice->structure);
        reorder_lists (img->type, img->currentSlice);

        if (img->structure==FRAME)
        {
          init_mbaff_lists();
        }

/*        if (img->frame_num==1) // write a reference list
        {
          count ++;
          if (count==1)
            for (i=0; i<listXsize[0]; i++)
              write_picture(listX[0][i], p_out2);
        }
*/

        // From here on, active_sps, active_pps and the slice header are valid
        if (img->MbaffFrameFlag)
          img->current_mb_nr = currSlice->start_mb_nr << 1;
        else
          img->current_mb_nr = currSlice->start_mb_nr;

        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset%8 != 0) 
          {
            ByteStartPosition++;
          }
          arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
        }
// printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
        FreeNALU(nalu);
        return current_header;
        break;
      case NALU_TYPE_DPA:
        //! The state machine here should follow the same ideas as the old readSliceRTP()
        //! basically:
        //! work on DPA (as above)
        //! read and process all following SEI/SPS/PPS/PD/Filler NALUs
        //! if next video NALU is dpB, 
        //!   then read and check whether it belongs to DPA, if yes, use it
        //! else
        //!   ;   // nothing
        //! read and process all following SEI/SPS/PPS/PD/Filler NALUs
        //! if next video NALU is dpC
        //!   then read and check whether it belongs to DPA (and DPB, if present), if yes, use it, done
        //! else
        //!   use the DPA (and the DPB if present)

        /* 
            LC: inserting the code related to DP processing, mainly copying some of the parts
            related to NALU_TYPE_SLICE, NALU_TYPE_IDR.
        */

        if(expected_slice_type != NALU_TYPE_DPA)
        {
          /* oops... we found the next slice, go back! */
          fseek(bits, ftell_position, SEEK_SET);
          FreeNALU(nalu);
          return current_header;
        }

        img->idr_flag          = (nalu->nal_unit_type == NALU_TYPE_IDR);
        img->nal_reference_idc = nalu->nal_reference_idc;
        img->disposable_flag   = (nalu->nal_reference_idc == NALU_PRIORITY_DISPOSABLE);
        currSlice->dp_mode     = PAR_DP_3;
        currSlice->max_part_nr = 3;
        currSlice->ei_flag     = 0;
        currStream             = currSlice->partArr[0].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
        
        BitsUsedByHeader     = FirstPartOfSliceHeader();
        UseParameterSet (currSlice->pic_parameter_set_id);
        BitsUsedByHeader    += RestOfSliceHeader ();
        
        FmoInit (active_pps, active_sps);
        
        if(is_new_picture())
        {
          init_picture(img, input);
          current_header = SOP;
          CheckZeroByteVCL(nalu, &ret);		  
        }
        else
          current_header = SOS;

        
        init_lists(img->type, img->currentSlice->structure);
        reorder_lists (img->type, img->currentSlice);
        
        if (img->structure==FRAME)
        {
          init_mbaff_lists();
        }

        // From here on, active_sps, active_pps and the slice header are valid
        if (img->MbaffFrameFlag)
          img->current_mb_nr = currSlice->start_mb_nr << 1;
        else
          img->current_mb_nr = currSlice->start_mb_nr;


        /* 
           LC:
              Now I need to read the slice ID, which depends on the value of 
              redundant_pic_cnt_present_flag (pag.49). 
        */
        
        slice_id_a  = ue_v("NALU:SLICE_A slice_idr", currStream);
        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset%8 != 0) 
          {
            ByteStartPosition++;
          }
          arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
        }
// printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
        break;
      case NALU_TYPE_DPB:
        /* LC: inserting the code related to DP processing */

        currStream             = currSlice->partArr[1].bitstream;
        currStream->ei_flag    = 0;
        currSlice->dp_mode     = PAR_DP_3;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        slice_id_b  = ue_v("NALU:SLICE_B slice_idr", currStream);
        if (active_pps->redundant_pic_cnt_present_flag)
          redundant_pic_cnt_b = ue_v("NALU:SLICE_B redudand_pic_cnt", currStream);
        else
          redundant_pic_cnt_b = 0;
        
        /*  LC: Initializing CABAC for the current data stream. */

        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset % 8 != 0) 
            ByteStartPosition++;
          
          arideco_start_decoding (&currSlice->partArr[1].de_cabac, currStream->streamBuffer, 
            ByteStartPosition, &currStream->read_len, img->type);
          
        }

        /* LC: resilience code to be inserted */
        /*         FreeNALU(nalu); */
        /*         return current_header; */

        break;
      case NALU_TYPE_DPC:
        /* LC: inserting the code related to DP processing */
        currSlice->dp_mode     = PAR_DP_3;
        currStream             = currSlice->partArr[2].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
        
        slice_id_c  = ue_v("NALU:SLICE_C slice_idr", currStream);
        if (active_pps->redundant_pic_cnt_present_flag)
          redundant_pic_cnt_c = ue_v("NALU:SLICE_C redudand_pic_cnt", currStream);
        else
          redundant_pic_cnt_c = 0;
        
        /* LC: Initializing CABAC for the current data stream. */

        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset % 8 != 0) 
            ByteStartPosition++;
          
          arideco_start_decoding (&currSlice->partArr[2].de_cabac, currStream->streamBuffer, 
            ByteStartPosition, &currStream->read_len, img->type);
        }

        /* LC: resilience code to be inserted */

        FreeNALU(nalu);
        return current_header;

        break;
      case NALU_TYPE_SEI:
        printf ("read_new_slice: Found NALU_TYPE_SEI, len %d\n", nalu->len);
        InterpretSEIMessage(nalu->buf,nalu->len,img);
        break;
      case NALU_TYPE_PPS:
        ProcessPPS(nalu);
        break;

      case NALU_TYPE_SPS:
        ProcessSPS(nalu);
        break;
      case NALU_TYPE_AUD:
//        printf ("read_new_slice: Found 'Access Unit Delimiter' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_EOSEQ:
//        printf ("read_new_slice: Found 'End of Sequence' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_EOSTREAM:
//        printf ("read_new_slice: Found 'End of Stream' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_FILL:
        printf ("read_new_slice: Found NALU_TYPE_FILL, len %d\n", nalu->len);
        printf ("Skipping these filling bits, proceeding w/ next NALU\n");
        break;
      default:
        printf ("Found NALU type %d, len %d undefined, ignore NALU, moving on\n", nalu->nal_unit_type, nalu->len);
    }
  }
  FreeNALU(nalu);

  return  current_header;
}


/*!
 ************************************************************************
 * \brief
 *    Partially reads new slice from bit_stream
 ************************************************************************
 */
int partial_read_new_slice()
{
  NALU_t *nalu = AllocNALU(MAX_CODED_FRAME_SIZE);
  int current_header = 0;
  int ret;
  int BitsUsedByHeader;
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream;

  int slice_id_a, slice_id_b, slice_id_c;
  int redundant_pic_cnt_b, redundant_pic_cnt_c;
  long ftell_position, expected_slice_type;

  //Modified by Wonsang You at FEB 25, 2007
  int i=0;
  //Modified-End

  //Modified by Wonsang You at MAY 7, 2007
  unsigned char *Buf;
  LARGE_INTEGER freq, g1, g2;
  double diff;
  //Modified-End
  
  expected_slice_type = NALU_TYPE_DPA;

  while (1)
  {
	//Modified by Wonsang You at FEB 8, 2007
	ftell_position = ftell(bitss);
	//Modified-End

	//Modified by Wonsang You at FEB 25, 2007
	if (input->FileFormat == PAR_OF_ANNEXB)
	{
	  if(first_frame)
	  {
		  if ((Buf = (unsigned char*)calloc (nalu->max_size , sizeof(char))) == NULL) no_mem_exit("GetAnnexbNALU: Buf");

		  /*
		  //START: measure the computational time		  
		  if(analysis_mode==1)
		  {
			  if(modify_tracking_flag==0)
			  {
					QueryPerformanceFrequency(&freq);
					QueryPerformanceCounter(&g1);
			  }
		  }
		  */
		  
		  i=0;
		  while(i<=start_frame_num)
		  {
			ret=partial_GetAnnexbNALU (nalu,Buf);
			ftell_position += nalu->len + 4;
			if((nalu->nal_unit_type==NALU_TYPE_SLICE)||(nalu->nal_unit_type==NALU_TYPE_IDR))
				i++;
			else
			{
				CheckZeroByteNonVCL(nalu, &ret);

				NALUtoRBSP(nalu);
				//    printf ("nalu->len %d\n", nalu->len);
    
				if (ret < 0)
				  printf ("Error while getting the NALU in file format %s, exit\n", input->FileFormat==PAR_OF_ANNEXB?"Annex B":"RTP");
				if (ret == 0)
				{
				  //      printf ("read_new_slice: returning %s\n", "EOS");
				  if(expected_slice_type != NALU_TYPE_DPA)
				  {
					/* oops... we found the next slice, go back! */
					//Modified by Wonsang You at FEB 8, 2007
					fseek(bitss, ftell_position, SEEK_SET);
					//Modified-End
					FreeNALU(nalu);
					return current_header;
				  }
				  else
					return EOS;
				}

				// Got a NALU
				if (nalu->forbidden_bit)
				{
				  printf ("Found NALU w/ forbidden_bit set, bit error?  Let's try...\n");
				}

				if(nalu->nal_unit_type==NALU_TYPE_SPS)
					ProcessSPS(nalu);
				else if(nalu->nal_unit_type==NALU_TYPE_PPS)
					ProcessPPS(nalu);
			}
		  }

		  /*
		  //END: measure the computational time
		  if(analysis_mode==1)
		  {
			  if(modify_tracking_flag==0)
			  {
					QueryPerformanceCounter(&g2);
					diff = (double)(g2.QuadPart - g1.QuadPart)/ freq.QuadPart;	
					fprintf(TT_FILE,"%d\t%f\n", motion_picture_num, diff); 
			  }
		  }
		  */
		  
		  motion_picture_num = -1;
		  free(Buf);
	  }
	  else
	  {
		  if ((Buf = (unsigned char*)calloc (nalu->max_size , sizeof(char))) == NULL) no_mem_exit("GetAnnexbNALU: Buf");
		  ret=partial_GetAnnexbNALU (nalu,Buf);
		  free(Buf);
	  }
	}
//    else
//      ret=GetRTPNALU (nalu);
	//Modified-End

    //In some cases, zero_byte shall be present. If current NALU is a VCL NALU, we can't tell
    //whether it is the first VCL NALU at this point, so only non-VCL NAL unit is checked here.
    CheckZeroByteNonVCL(nalu, &ret);

    NALUtoRBSP(nalu);
//    printf ("nalu->len %d\n", nalu->len);
    
    if (ret < 0)
      printf ("Error while getting the NALU in file format %s, exit\n", input->FileFormat==PAR_OF_ANNEXB?"Annex B":"RTP");
    if (ret == 0)
    {
//      printf ("read_new_slice: returning %s\n", "EOS");
      if(expected_slice_type != NALU_TYPE_DPA)
      {
        /* oops... we found the next slice, go back! */
		//Modified by Wonsang You at FEB 8, 2007
        fseek(bitss, ftell_position, SEEK_SET);
		//Modified-End
        FreeNALU(nalu);
        return current_header;
      }
      else
        return EOS;
    }

    // Got a NALU
    if (nalu->forbidden_bit)
    {
      printf ("Found NALU w/ forbidden_bit set, bit error?  Let's try...\n");
    }

    switch (nalu->nal_unit_type)
    {
      case NALU_TYPE_SLICE:
      case NALU_TYPE_IDR:
        img->idr_flag = (nalu->nal_unit_type == NALU_TYPE_IDR);
        img->nal_reference_idc = nalu->nal_reference_idc;
        img->disposable_flag = (nalu->nal_reference_idc == NALU_PRIORITY_DISPOSABLE);
        currSlice->dp_mode = PAR_DP_1;
        currSlice->max_part_nr = 1;
        currSlice->ei_flag = 0;
        currStream = currSlice->partArr[0].bitstream;
        currStream->ei_flag = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        // Some syntax of the Slice Header depends on the parameter set, which depends on
        // the parameter set ID of the SLice header.  Hence, read the pic_parameter_set_id
        // of the slice header first, then setup the active parameter sets, and then read
        // the rest of the slice header
        BitsUsedByHeader = FirstPartOfSliceHeader();
        UseParameterSet (currSlice->pic_parameter_set_id);
        BitsUsedByHeader+= RestOfSliceHeader ();

        FmoInit (active_pps, active_sps);

        AssignQuantParam (active_pps, active_sps);

        if(is_new_picture())
        {
		  //Added by Wonsang You at FEB 6, 2007
		  motion_picture_num++;
		  //if(motion_picture_num>0)		
			//partial_store_motion_into_buffer();			  
		  //Added-End

          init_picture(img, input);
          
          current_header = SOP;
          //check zero_byte if it is also the first NAL unit in the access unit
          CheckZeroByteVCL(nalu, &ret);
        }
        else
          current_header = SOS;
  
        init_lists(img->type, img->currentSlice->structure);
        reorder_lists (img->type, img->currentSlice);

        if (img->structure==FRAME)
        {
          init_mbaff_lists();
        }

/*        if (img->frame_num==1) // write a reference list
        {
          count ++;
          if (count==1)
            for (i=0; i<listXsize[0]; i++)
              write_picture(listX[0][i], p_out2);
        }
*/

        // From here on, active_sps, active_pps and the slice header are valid
        if (img->MbaffFrameFlag)
          img->current_mb_nr = currSlice->start_mb_nr << 1;
        else
          img->current_mb_nr = currSlice->start_mb_nr;

        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset%8 != 0) 
          {
            ByteStartPosition++;
          }
          arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
        }
// printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
        FreeNALU(nalu);
        return current_header;
        break;
      case NALU_TYPE_DPA:
        //! The state machine here should follow the same ideas as the old readSliceRTP()
        //! basically:
        //! work on DPA (as above)
        //! read and process all following SEI/SPS/PPS/PD/Filler NALUs
        //! if next video NALU is dpB, 
        //!   then read and check whether it belongs to DPA, if yes, use it
        //! else
        //!   ;   // nothing
        //! read and process all following SEI/SPS/PPS/PD/Filler NALUs
        //! if next video NALU is dpC
        //!   then read and check whether it belongs to DPA (and DPB, if present), if yes, use it, done
        //! else
        //!   use the DPA (and the DPB if present)

        /* 
            LC: inserting the code related to DP processing, mainly copying some of the parts
            related to NALU_TYPE_SLICE, NALU_TYPE_IDR.
        */

        if(expected_slice_type != NALU_TYPE_DPA)
        {
          /* oops... we found the next slice, go back! */
		  //Modified by Wonsang You at FEB 8, 2007
          fseek(bitss, ftell_position, SEEK_SET);
		  //Modified-End
          FreeNALU(nalu);
          return current_header;
        }

        img->idr_flag          = (nalu->nal_unit_type == NALU_TYPE_IDR);
        img->nal_reference_idc = nalu->nal_reference_idc;
        img->disposable_flag   = (nalu->nal_reference_idc == NALU_PRIORITY_DISPOSABLE);
        currSlice->dp_mode     = PAR_DP_3;
        currSlice->max_part_nr = 3;
        currSlice->ei_flag     = 0;
        currStream             = currSlice->partArr[0].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
        
        BitsUsedByHeader     = FirstPartOfSliceHeader();
        UseParameterSet (currSlice->pic_parameter_set_id);
        BitsUsedByHeader    += RestOfSliceHeader ();
        
        FmoInit (active_pps, active_sps);
        
        if(is_new_picture())
        {
          init_picture(img, input);
          current_header = SOP;
          CheckZeroByteVCL(nalu, &ret);		  
        }
        else
          current_header = SOS;

        
        init_lists(img->type, img->currentSlice->structure);
        reorder_lists (img->type, img->currentSlice);
        
        if (img->structure==FRAME)
        {
          init_mbaff_lists();
        }

        // From here on, active_sps, active_pps and the slice header are valid
        if (img->MbaffFrameFlag)
          img->current_mb_nr = currSlice->start_mb_nr << 1;
        else
          img->current_mb_nr = currSlice->start_mb_nr;


        /* 
           LC:
              Now I need to read the slice ID, which depends on the value of 
              redundant_pic_cnt_present_flag (pag.49). 
        */
        
        slice_id_a  = ue_v("NALU:SLICE_A slice_idr", currStream);
        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset%8 != 0) 
          {
            ByteStartPosition++;
          }
          arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
        }
// printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
        break;
      case NALU_TYPE_DPB:
        /* LC: inserting the code related to DP processing */

        currStream             = currSlice->partArr[1].bitstream;
        currStream->ei_flag    = 0;
        currSlice->dp_mode     = PAR_DP_3;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        slice_id_b  = ue_v("NALU:SLICE_B slice_idr", currStream);
        if (active_pps->redundant_pic_cnt_present_flag)
          redundant_pic_cnt_b = ue_v("NALU:SLICE_B redudand_pic_cnt", currStream);
        else
          redundant_pic_cnt_b = 0;
        
        /*  LC: Initializing CABAC for the current data stream. */

        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset % 8 != 0) 
            ByteStartPosition++;
          
          arideco_start_decoding (&currSlice->partArr[1].de_cabac, currStream->streamBuffer, 
            ByteStartPosition, &currStream->read_len, img->type);
          
        }

        /* LC: resilience code to be inserted */
        /*         FreeNALU(nalu); */
        /*         return current_header; */

        break;
      case NALU_TYPE_DPC:
        /* LC: inserting the code related to DP processing */
        currSlice->dp_mode     = PAR_DP_3;
        currStream             = currSlice->partArr[2].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
        
        slice_id_c  = ue_v("NALU:SLICE_C slice_idr", currStream);
        if (active_pps->redundant_pic_cnt_present_flag)
          redundant_pic_cnt_c = ue_v("NALU:SLICE_C redudand_pic_cnt", currStream);
        else
          redundant_pic_cnt_c = 0;
        
        /* LC: Initializing CABAC for the current data stream. */

        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset % 8 != 0) 
            ByteStartPosition++;
          
          arideco_start_decoding (&currSlice->partArr[2].de_cabac, currStream->streamBuffer, 
            ByteStartPosition, &currStream->read_len, img->type);
        }

        /* LC: resilience code to be inserted */

        FreeNALU(nalu);
        return current_header;

        break;
      case NALU_TYPE_SEI:
        printf ("read_new_slice: Found NALU_TYPE_SEI, len %d\n", nalu->len);
        InterpretSEIMessage(nalu->buf,nalu->len,img);
        break;
      case NALU_TYPE_PPS:
        ProcessPPS(nalu);
        break;

      case NALU_TYPE_SPS:
        ProcessSPS(nalu);
        break;
      case NALU_TYPE_AUD:
//        printf ("read_new_slice: Found 'Access Unit Delimiter' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_EOSEQ:
//        printf ("read_new_slice: Found 'End of Sequence' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_EOSTREAM:
//        printf ("read_new_slice: Found 'End of Stream' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_FILL:
        printf ("read_new_slice: Found NALU_TYPE_FILL, len %d\n", nalu->len);
        printf ("Skipping these filling bits, proceeding w/ next NALU\n");
        break;
      default:
        printf ("Found NALU type %d, len %d undefined, ignore NALU, moving on\n", nalu->nal_unit_type, nalu->len);
    }
  }
  FreeNALU(nalu);

  return  current_header;
}


/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new picture
 ************************************************************************
 */
void init_picture(struct img_par *img, struct inp_par *inp)
{
  int i,k,l;
  Slice *currSlice = img->currentSlice;

  if (dec_picture)
  {
    // this may only happen on slice loss
    exit_picture();
  }

  if (img->frame_num != img->pre_frame_num && img->frame_num != (img->pre_frame_num + 1) % img->MaxFrameNum) 
  {
    if (active_sps->gaps_in_frame_num_value_allowed_flag == 0)
    {
      // picture error concealment
      if(inp->conceal_mode !=0)
      {
        if((img->frame_num) < ((img->pre_frame_num + 1) % img->MaxFrameNum))
        {
          // Conceal lost IDR frames and any frames immediately 
          //   following the IDR. Use frame copy for these since 
          //   lists cannot be formed correctly for motion copy
          img->conceal_mode = 1; 
          img->IDR_concealment_flag = 1;
          conceal_lost_frames(img);
          //reset to original concealment mode for future drops
          img->conceal_mode = inp->conceal_mode;
        }
        else
        {
          //reset to original concealment mode for future drops
          img->conceal_mode = inp->conceal_mode;

          img->IDR_concealment_flag = 0;
          conceal_lost_frames(img);
        }
      }
      else
      {   // Advanced Error Concealment would be called here to combat unintentional loss of pictures.
        error("An unintentional loss of pictures occurs! Exit\n", 100);
      }
    }
    if(img->conceal_mode == 0)
      fill_frame_num_gap(img);
  }

  if(img->nal_reference_idc)
  {
    img->pre_frame_num = img->frame_num;
  }

  //img->num_dec_mb = 0;
  
  //calculate POC
  decode_poc(img);

  if(img->nal_reference_idc)
    img->last_ref_pic_poc = img->framepoc;

  //  dumppoc (img);

  if (img->structure==FRAME ||img->structure==TOP_FIELD)
  {
#ifdef WIN32
    _ftime (&(img->tstruct_start));             // start time ms
#else
    ftime (&(img->tstruct_start));              // start time ms
#endif
    time( &(img->ltime_start));                // start time s
  }

  dec_picture = alloc_storable_picture (img->structure, img->width, img->height, img->width_cr, img->height_cr);
  dec_picture->top_poc=img->toppoc;
  dec_picture->bottom_poc=img->bottompoc;
  dec_picture->frame_poc=img->framepoc;
  dec_picture->qp=img->qp;
  dec_picture->slice_qp_delta=currSlice->slice_qp_delta;
  dec_picture->chroma_qp_offset[0] = active_pps->chroma_qp_index_offset;
  dec_picture->chroma_qp_offset[1] = active_pps->second_chroma_qp_index_offset;

  // reset all variables of the error concealment instance before decoding of every frame.
  // here the third parameter should, if perfectly, be equal to the number of slices per frame.
  // using little value is ok, the code will allocate more memory if the slice number is larger
  ercReset(erc_errorVar, img->PicSizeInMbs, img->PicSizeInMbs, dec_picture->size_x);
  erc_mvperMB = 0;

  switch (img->structure )
  {
  case TOP_FIELD:
    {
      dec_picture->poc=img->toppoc;
      img->number *= 2;
      break;
    }
  case BOTTOM_FIELD:
    {
      dec_picture->poc=img->bottompoc;
      img->number++;
      break;
    }
  case FRAME:
    {
      dec_picture->poc=img->framepoc;
      break;
    }
  default:
    error("img->structure not initialized", 235);
  }
    
  img->current_slice_nr=0;

  if (img->type > SI_SLICE)
  {
    set_ec_flag(SE_PTYPE);
    img->type = P_SLICE;  // concealed element
  }

  // CAVLC init
  for (i=0;i < (int)img->PicSizeInMbs; i++)
    for (k=0;k<4;k++)
      for (l=0;l<(4 + img->num_blk8x8_uv);l++)
        img->nz_coeff[i][k][l]=-1;  // CAVLC

  if(active_pps->constrained_intra_pred_flag)
  {
    for (i=0; i<(int)img->PicSizeInMbs; i++)
    {
      img->intra_block[i] = 1;
    }
  }

  // Set the slice_nr member of each MB to -1, to ensure correct when packet loss occurs
  // TO set Macroblock Map (mark all MBs as 'have to be concealed')
  for(i=0; i<(int)img->PicSizeInMbs; i++)
  {
    img->mb_data[i].slice_nr = -1; 
    img->mb_data[i].ei_flag = 1;
  }

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0; // define vertical positions
  img->block_x = img->pix_x = img->pix_c_x = 0; // define horizontal positions

  dec_picture->slice_type = img->type;
  dec_picture->used_for_reference = (img->nal_reference_idc != 0);
  dec_picture->idr_flag = img->idr_flag;
  dec_picture->no_output_of_prior_pics_flag = img->no_output_of_prior_pics_flag;
  dec_picture->long_term_reference_flag = img->long_term_reference_flag;
  dec_picture->adaptive_ref_pic_buffering_flag = img->adaptive_ref_pic_buffering_flag;

  dec_picture->dec_ref_pic_marking_buffer = img->dec_ref_pic_marking_buffer;
  img->dec_ref_pic_marking_buffer = NULL;

  dec_picture->MbaffFrameFlag = img->MbaffFrameFlag;
  dec_picture->PicWidthInMbs = img->PicWidthInMbs;
  dec_picture->pic_num = img->frame_num;
  dec_picture->frame_num = img->frame_num;
  dec_picture->coded_frame = (img->structure==FRAME);

  dec_picture->chroma_format_idc = active_sps->chroma_format_idc;

  dec_picture->frame_mbs_only_flag = active_sps->frame_mbs_only_flag;
  dec_picture->frame_cropping_flag = active_sps->frame_cropping_flag;

  if (dec_picture->frame_cropping_flag)
  {
    dec_picture->frame_cropping_rect_left_offset   = active_sps->frame_cropping_rect_left_offset;
    dec_picture->frame_cropping_rect_right_offset  = active_sps->frame_cropping_rect_right_offset;
    dec_picture->frame_cropping_rect_top_offset    = active_sps->frame_cropping_rect_top_offset;
    dec_picture->frame_cropping_rect_bottom_offset = active_sps->frame_cropping_rect_bottom_offset;
  }
}


/*!
 ************************************************************************
 * \brief
 *    finish decoding of a picture, conceal errors and store it 
 *    into the DPB
 ************************************************************************
 */
void exit_picture()
{
  char yuv_types[4][6]= {"4:0:0","4:2:0","4:2:2","4:4:4"};
  int ercStartMB;
  int ercSegment;
  frame recfr;
  unsigned int i;
  int structure, frame_poc, slice_type, refpic, qp, pic_num, chroma_format_idc;

  int tmp_time;                   // time used by decoding the last frame
  char yuvFormat[10];

  // return if the last picture has already been finished
  if (dec_picture==NULL)
  {
    return;
  }

  //deblocking for frame or field
  DeblockPicture( img, dec_picture );

  if (dec_picture->MbaffFrameFlag)
    MbAffPostProc();

  recfr.yptr = &dec_picture->imgY[0][0];
  if (dec_picture->chroma_format_idc != YUV400)
  {
    recfr.uptr = &dec_picture->imgUV[0][0][0];
    recfr.vptr = &dec_picture->imgUV[1][0][0];
  }

  //! this is always true at the beginning of a picture
  ercStartMB = 0;
  ercSegment = 0;

  //! mark the start of the first segment
  if (!dec_picture->MbaffFrameFlag)
  {
    ercStartSegment(0, ercSegment, 0 , erc_errorVar);
    //! generate the segments according to the macroblock map
    for(i = 1; i<dec_picture->PicSizeInMbs; i++)
    {
      if(img->mb_data[i].ei_flag != img->mb_data[i-1].ei_flag)
      {
        ercStopSegment(i-1, ercSegment, 0, erc_errorVar); //! stop current segment
        
        //! mark current segment as lost or OK
        if(img->mb_data[i-1].ei_flag)
          ercMarkCurrSegmentLost(dec_picture->size_x, erc_errorVar);
        else
          ercMarkCurrSegmentOK(dec_picture->size_x, erc_errorVar);
        
        ercSegment++;  //! next segment
        ercStartSegment(i, ercSegment, 0 , erc_errorVar); //! start new segment
        ercStartMB = i;//! save start MB for this segment 
      }
    }
    //! mark end of the last segment
    ercStopSegment(dec_picture->PicSizeInMbs-1, ercSegment, 0, erc_errorVar);
    if(img->mb_data[i-1].ei_flag)
      ercMarkCurrSegmentLost(dec_picture->size_x, erc_errorVar);
    else
      ercMarkCurrSegmentOK(dec_picture->size_x, erc_errorVar);
    
    //! call the right error concealment function depending on the frame type.
    erc_mvperMB /= dec_picture->PicSizeInMbs;
    
    erc_img = img;
    if(dec_picture->slice_type == I_SLICE || dec_picture->slice_type == SI_SLICE) // I-frame
      ercConcealIntraFrame(&recfr, dec_picture->size_x, dec_picture->size_y, erc_errorVar);
    else
      ercConcealInterFrame(&recfr, erc_object_list, dec_picture->size_x, dec_picture->size_y, erc_errorVar, dec_picture->chroma_format_idc);
  }

  if (img->structure == FRAME)         // buffer mgt. for frame mode
    frame_postprocessing(img, input);
  else
    field_postprocessing(img, input);   // reset all interlaced variables

  structure  = dec_picture->structure;
  slice_type = dec_picture->slice_type;
  frame_poc  = dec_picture->frame_poc;
  refpic     = dec_picture->used_for_reference;
  qp         = dec_picture->qp;
  pic_num    = dec_picture->pic_num;

  chroma_format_idc= dec_picture->chroma_format_idc;

  store_picture_in_dpb(dec_picture);
  dec_picture=NULL;

  if (img->last_has_mmco_5)
  {
    img->pre_frame_num = 0;
  }

  if ((structure==FRAME)||structure==BOTTOM_FIELD)
  {
    
#ifdef WIN32
    _ftime (&(img->tstruct_end));             // start time ms
#else
    ftime (&(img->tstruct_end));              // start time ms
#endif
    
    time( &(img->ltime_end));                // start time s

    tmp_time=(img->ltime_end*1000+img->tstruct_end.millitm) - (img->ltime_start*1000+img->tstruct_start.millitm);
    tot_time=tot_time + tmp_time;

    sprintf(yuvFormat,"%s", yuv_types[chroma_format_idc]);

	/*    
    if(slice_type == I_SLICE) // I picture
      fprintf(stdout,"%04d(I)  %8d %5d %5d %7.4f %7.4f %7.4f  %s %5d\n",
      frame_no, frame_poc, pic_num, qp, snr->snr_y, snr->snr_u, snr->snr_v, yuvFormat, tmp_time);
    else if(slice_type == P_SLICE) // P pictures
      fprintf(stdout,"%04d(P)  %8d %5d %5d %7.4f %7.4f %7.4f  %s %5d\n",
      frame_no, frame_poc, pic_num, qp, snr->snr_y, snr->snr_u, snr->snr_v, yuvFormat, tmp_time);
    else if(slice_type == SP_SLICE) // SP pictures
      fprintf(stdout,"%04d(SP) %8d %5d %5d %7.4f %7.4f %7.4f  %s %5d\n",
      frame_no, frame_poc, pic_num, qp, snr->snr_y, snr->snr_u, snr->snr_v, yuvFormat, tmp_time);
    else if (slice_type == SI_SLICE)
      fprintf(stdout,"%04d(SI) %8d %5d %5d %7.4f %7.4f %7.4f  %s %5d\n",
      frame_no, frame_poc, pic_num, qp, snr->snr_y, snr->snr_u, snr->snr_v, yuvFormat, tmp_time);
    else if(refpic) // stored B pictures
      fprintf(stdout,"%04d(RB) %8d %5d %5d %7.4f %7.4f %7.4f  %s %5d\n",
      frame_no, frame_poc, pic_num, qp, snr->snr_y, snr->snr_u, snr->snr_v, yuvFormat, tmp_time);
    else // B pictures
      fprintf(stdout,"%04d(B)  %8d %5d %5d %7.4f %7.4f %7.4f  %s %5d\n",
      frame_no, frame_poc, pic_num, qp, snr->snr_y, snr->snr_u, snr->snr_v, yuvFormat, tmp_time);

    fflush(stdout);
*/
    if(slice_type == I_SLICE || slice_type == SI_SLICE || slice_type == P_SLICE || refpic)   // I or P pictures
      img->number++;
    else
      Bframe_ctr++;    // B pictures
    snr->frame_ctr++;

    g_nFrame++;
  }

  img->current_mb_nr = -4712;   // impossible value for debugging, StW
  img->current_slice_nr = 0;

}

/*!
 ************************************************************************
 * \brief
 *    write the encoding mode and motion vectors of current 
 *    MB to the buffer of the error concealment module.
 ************************************************************************
 */

void ercWriteMBMODEandMV(struct img_par *img,struct inp_par *inp)
{
  extern objectBuffer_t *erc_object_list;
  int i, ii, jj, currMBNum = img->current_mb_nr;
  int mbx = xPosMB(currMBNum,dec_picture->size_x), mby = yPosMB(currMBNum,dec_picture->size_x);
  objectBuffer_t *currRegion, *pRegion;
  Macroblock *currMB = &img->mb_data[currMBNum];
  short***  mv;

  currRegion = erc_object_list + (currMBNum<<2);

  if(img->type != B_SLICE) //non-B frame
  {
    for (i=0; i<4; i++)
    {
      pRegion             = currRegion + i;
      pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
                             currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  :
                             currMB->b8mode[i]==0      ? REGMODE_INTER_COPY :
                             currMB->b8mode[i]==1      ? REGMODE_INTER_PRED : REGMODE_INTER_PRED_8x8);
      if (currMB->b8mode[i]==0 || currMB->b8mode[i]==IBLOCK)  // INTRA OR COPY
      {
        pRegion->mv[0]    = 0;
        pRegion->mv[1]    = 0;
        pRegion->mv[2]    = 0;
      }
      else
      {
        ii              = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
        jj              = 4*mby + (i/2)*2;
        if (currMB->b8mode[i]>=5 && currMB->b8mode[i]<=7)  // SMALL BLOCKS
        {
          pRegion->mv[0]  = (dec_picture->mv[LIST_0][jj][ii][0] + dec_picture->mv[LIST_0][jj][ii+1][0] + dec_picture->mv[LIST_0][jj+1][ii][0] + dec_picture->mv[LIST_0][jj+1][ii+1][0] + 2)/4;
          pRegion->mv[1]  = (dec_picture->mv[LIST_0][jj][ii][1] + dec_picture->mv[LIST_0][jj][ii+1][1] + dec_picture->mv[LIST_0][jj+1][ii][1] + dec_picture->mv[LIST_0][jj+1][ii+1][1] + 2)/4;
        }
        else // 16x16, 16x8, 8x16, 8x8
        {
          pRegion->mv[0]  = dec_picture->mv[LIST_0][jj][ii][0];
          pRegion->mv[1]  = dec_picture->mv[LIST_0][jj][ii][1];
//          pRegion->mv[0]  = dec_picture->mv[LIST_0][4*mby+(i/2)*2][4*mbx+(i%2)*2+BLOCK_SIZE][0];
//          pRegion->mv[1]  = dec_picture->mv[LIST_0][4*mby+(i/2)*2][4*mbx+(i%2)*2+BLOCK_SIZE][1];
        }
        erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);
        pRegion->mv[2]    = dec_picture->ref_idx[LIST_0][jj][ii];
      }
    }
  }
  else  //B-frame
  {
    for (i=0; i<4; i++)
    {
      ii                  = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
      jj                  = 4*mby + (i/2)*2;
      pRegion             = currRegion + i;
      pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
                             currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  : REGMODE_INTER_PRED_8x8);
      if (currMB->mb_type==I16MB || currMB->b8mode[i]==IBLOCK)  // INTRA
      {
        pRegion->mv[0]    = 0;
        pRegion->mv[1]    = 0;
        pRegion->mv[2]    = 0;
      }
      else
      {
        int idx = (dec_picture->ref_idx[0][jj][ii]<0)?1:0;
//        int idx = (currMB->b8mode[i]==0 && currMB->b8pdir[i]==2 ? LIST_0 : currMB->b8pdir[i]==1 ? LIST_1 : LIST_0);
//        int idx = currMB->b8pdir[i]==0 ? LIST_0 : LIST_1;
        mv                = dec_picture->mv[idx];
        pRegion->mv[0]    = (mv[jj][ii][0] + mv[jj][ii+1][0] + mv[jj+1][ii][0] + mv[jj+1][ii+1][0] + 2)/4;
        pRegion->mv[1]    = (mv[jj][ii][1] + mv[jj][ii+1][1] + mv[jj+1][ii][1] + mv[jj+1][ii+1][1] + 2)/4;
        erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);

        pRegion->mv[2]  = (dec_picture->ref_idx[idx][jj][ii]);
/*        
        if (currMB->b8pdir[i]==0 || (currMB->b8pdir[i]==2 && currMB->b8mode[i]!=0)) // forward or bidirect
        {
          pRegion->mv[2]  = (dec_picture->ref_idx[LIST_0][jj][ii]);
          ///???? is it right, not only "img->fw_refFrArr[jj][ii-4]"
        }
        else
        {
          pRegion->mv[2]  = (dec_picture->ref_idx[LIST_1][jj][ii]);
//          pRegion->mv[2]  = 0;
        }
        */
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    set defaults for old_slice
 *    NAL unit of a picture"
 ************************************************************************
 */
void init_old_slice()
{
  old_slice.field_pic_flag = 0;

  old_slice.pps_id = INT_MAX;

  old_slice.frame_num = INT_MAX;

  old_slice.nal_ref_idc = INT_MAX;
  
  old_slice.idr_flag = 0;

  old_slice.pic_oder_cnt_lsb          = UINT_MAX;
  old_slice.delta_pic_oder_cnt_bottom = INT_MAX;

  old_slice.delta_pic_order_cnt[0] = INT_MAX;
  old_slice.delta_pic_order_cnt[1] = INT_MAX;

}

/*!
 ************************************************************************
 * \brief
 *    save slice parameters that are needed for checking of "first VCL
 *    NAL unit of a picture"
 ************************************************************************
 */
void exit_slice()
{

  old_slice.pps_id = img->currentSlice->pic_parameter_set_id;

  old_slice.frame_num = img->frame_num;

  old_slice.field_pic_flag = img->field_pic_flag;

  if(img->field_pic_flag)
  {
    old_slice.bottom_field_flag = img->bottom_field_flag;
  }

  old_slice.nal_ref_idc   = img->nal_reference_idc;
  
  old_slice.idr_flag = img->idr_flag;
  if (img->idr_flag)
  {
    old_slice.idr_pic_id = img->idr_pic_id;
  }

  if (active_sps->pic_order_cnt_type == 0)
  {
    old_slice.pic_oder_cnt_lsb          = img->pic_order_cnt_lsb;
    old_slice.delta_pic_oder_cnt_bottom = img->delta_pic_order_cnt_bottom;
  }

  if (active_sps->pic_order_cnt_type == 1)
  {
    old_slice.delta_pic_order_cnt[0] = img->delta_pic_order_cnt[0];
    old_slice.delta_pic_order_cnt[1] = img->delta_pic_order_cnt[1];
  }
}

/*!
 ************************************************************************
 * \brief
 *    detect if current slice is "first VCL NAL unit of a picture"
 ************************************************************************
 */
int is_new_picture()
{
  int result=0;

  result |= (old_slice.pps_id != img->currentSlice->pic_parameter_set_id);

  result |= (old_slice.frame_num != img->frame_num);

  result |= (old_slice.field_pic_flag != img->field_pic_flag);

  if(img->field_pic_flag && old_slice.field_pic_flag)
  {
    result |= (old_slice.bottom_field_flag != img->bottom_field_flag);
  }

  result |= (old_slice.nal_ref_idc != img->nal_reference_idc) && ((old_slice.nal_ref_idc == 0) || (img->nal_reference_idc == 0));
  
  result |= ( old_slice.idr_flag != img->idr_flag);

  if (img->idr_flag && old_slice.idr_flag)
  {
    result |= (old_slice.idr_pic_id != img->idr_pic_id);
  }

  if (active_sps->pic_order_cnt_type == 0)
  {
    result |=  (old_slice.pic_oder_cnt_lsb          != img->pic_order_cnt_lsb);
    result |=  (old_slice.delta_pic_oder_cnt_bottom != img->delta_pic_order_cnt_bottom);
  }

  if (active_sps->pic_order_cnt_type == 1)
  {
    result |= (old_slice.delta_pic_order_cnt[0] != img->delta_pic_order_cnt[0]);
    result |= (old_slice.delta_pic_order_cnt[1] != img->delta_pic_order_cnt[1]);
  }

  return result;
}


/*!
 ************************************************************************
 * \brief
 *    decodes one slice
 ************************************************************************
 */
void decode_one_slice(struct img_par *img,struct inp_par *inp)
{

  Boolean end_of_slice = FALSE;
  int read_flag;
  img->cod_counter=-1;

  set_ref_pic_num();

  if (img->type == B_SLICE)
      compute_colocated(Co_located, listX);

  //reset_ec_flags();

  while (end_of_slice == FALSE) // loop over macroblocks
  {

#if TRACEA
  fprintf(p_trace,"\n*********** POC: %i (I/P) MB: %i Slice: %i Type %d **********\n", img->ThisPOC, img->current_mb_nr, img->current_slice_nr, img->type);
#endif

	// Initializes the current macroblock
    start_macroblock(img,inp, img->current_mb_nr);
    // Get the syntax elements from the NAL
    read_flag = read_one_macroblock(img,inp);
    decode_one_macroblock(img,inp);

	//Added by Wonsang You at Oct-25-2006
	//store_mbinfo_into_buffer(img);
	//Added-End

    if(img->MbaffFrameFlag && dec_picture->mb_field[img->current_mb_nr])
    {
      img->num_ref_idx_l0_active >>= 1;
      img->num_ref_idx_l1_active >>= 1;
    }

    ercWriteMBMODEandMV(img,inp);

    end_of_slice=exit_macroblock(img,inp,(!img->MbaffFrameFlag||img->current_mb_nr%2));
  }

  exit_slice();
  //reset_ec_flags();

}

void partial_decode_one_slice(struct img_par *img,struct inp_par *inp)
{
  Boolean end_of_slice = FALSE;
  int read_flag;
  img->cod_counter=-1;

  set_ref_pic_num();

  if (img->type == B_SLICE)
      compute_colocated(Co_located, listX);

  //reset_ec_flags();

  while (end_of_slice == FALSE) // loop over macroblocks
  {

#if TRACEA
  fprintf(p_trace,"\n*********** POC: %i (I/P) MB: %i Slice: %i Type %d **********\n", img->ThisPOC, img->current_mb_nr, img->current_slice_nr, img->type);
#endif

	// Initializes the current macroblock
	start_macroblock(img,inp, img->current_mb_nr);
	// Get the syntax elements from the NAL
	read_flag = read_one_macroblock(img,inp);
	
	//Modified by Wonsang You at AUG-26-2007
	if(blk_dec_enb[img->current_mb_nr] == 1)
	{	
		partial_decode_one_macroblock(img,inp,blk_list[img->current_mb_nr]);
	}

	//Modified-End
	
	if(img->MbaffFrameFlag && dec_picture->mb_field[img->current_mb_nr])
	{
	  img->num_ref_idx_l0_active >>= 1;
	  img->num_ref_idx_l1_active >>= 1;
	}

	ercWriteMBMODEandMV(img,inp);

	end_of_slice=exit_macroblock(img,inp,(!img->MbaffFrameFlag||img->current_mb_nr%2));
	
  }

  exit_slice();
  //reset_ec_flags();

}


void decode_slice(struct img_par *img,struct inp_par *inp, int current_header)
{
  Slice *currSlice = img->currentSlice;

  if (active_pps->entropy_coding_mode_flag)
  {
    init_contexts (img);
    cabac_new_slice();
  }

  if ( (active_pps->weighted_bipred_idc > 0  && (img->type == B_SLICE)) || (active_pps->weighted_pred_flag && img->type !=I_SLICE))
    fill_wp_params(img);

  //printf("frame picture %d %d %d\n",img->structure,img->ThisPOC,img->direct_spatial_mv_pred_flag);
  

  // decode main slice information
  if ((current_header == SOP || current_header == SOS) && currSlice->ei_flag == 0)
    decode_one_slice(img,inp);
    
  // setMB-Nr in case this slice was lost
//  if(currSlice->ei_flag)  
//    img->current_mb_nr = currSlice->last_mb_nr + 1;

}


void partial_decode_slice(struct img_par *img,struct inp_par *inp, int current_header)
{
  Slice *currSlice = img->currentSlice;

  if (active_pps->entropy_coding_mode_flag)
  {
    init_contexts (img);
    cabac_new_slice();
  }

  if ( (active_pps->weighted_bipred_idc > 0  && (img->type == B_SLICE)) || (active_pps->weighted_pred_flag && img->type !=I_SLICE))
    fill_wp_params(img);

  //printf("frame picture %d %d %d\n",img->structure,img->ThisPOC,img->direct_spatial_mv_pred_flag);
  

  // decode main slice information
  if ((current_header == SOP || current_header == SOS) && currSlice->ei_flag == 0)
    partial_decode_one_slice(img,inp);
    
  // setMB-Nr in case this slice was lost
//  if(currSlice->ei_flag)  
//    img->current_mb_nr = currSlice->last_mb_nr + 1;

}

/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after frame decoding
 ************************************************************************
 */
void frame_postprocessing(struct img_par *img, struct inp_par *inp)
{
}

/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after field decoding
 ************************************************************************
 */
void field_postprocessing(struct img_par *img, struct inp_par *inp)
{
  img->number /= 2;
}



void reset_wp_params(struct img_par *img)
{
  int i,comp;
  int log_weight_denom;

  for (i=0; i<MAX_REFERENCE_PICTURES; i++)
  {
    for (comp=0; comp<3; comp++)
    {
      log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
      img->wp_weight[0][i][comp] = 1<<log_weight_denom;
      img->wp_weight[1][i][comp] = 1<<log_weight_denom;
    }
  }
}


void fill_wp_params(struct img_par *img)
{
  int i, j, k;
  int comp;
  int log_weight_denom;
  int tb, td;
  int bframe = (img->type==B_SLICE);
  int max_bwd_ref, max_fwd_ref;
  int tx,DistScaleFactor;

  max_fwd_ref = img->num_ref_idx_l0_active;
  max_bwd_ref = img->num_ref_idx_l1_active;

  if (active_pps->weighted_bipred_idc == 2 && bframe)
  {
    img->luma_log2_weight_denom = 5;
    img->chroma_log2_weight_denom = 5;
    img->wp_round_luma = 16;
    img->wp_round_chroma = 16;

    for (i=0; i<MAX_REFERENCE_PICTURES; i++)
    {
      for (comp=0; comp<3; comp++)
      {
        log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
        img->wp_weight[0][i][comp] = 1<<log_weight_denom;
        img->wp_weight[1][i][comp] = 1<<log_weight_denom;
        img->wp_offset[0][i][comp] = 0;
        img->wp_offset[1][i][comp] = 0;
      }
    }
  }

  if (bframe)
  {
    for (i=0; i<max_fwd_ref; i++)
    {
      for (j=0; j<max_bwd_ref; j++)
      {
        for (comp = 0; comp<3; comp++)
        {
          log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
          if (active_pps->weighted_bipred_idc == 1)
          {
            img->wbp_weight[0][i][j][comp] =  img->wp_weight[0][i][comp];
            img->wbp_weight[1][i][j][comp] =  img->wp_weight[1][j][comp];
          }
          else if (active_pps->weighted_bipred_idc == 2)
          {
            td = Clip3(-128,127,listX[LIST_1][j]->poc - listX[LIST_0][i]->poc);
            if (td == 0 || listX[LIST_1][j]->is_long_term || listX[LIST_0][i]->is_long_term)
            {
              img->wbp_weight[0][i][j][comp] =   32;
              img->wbp_weight[1][i][j][comp] =   32;
            }
            else
            {
              tb = Clip3(-128,127,img->ThisPOC - listX[LIST_0][i]->poc);

              tx = (16384 + abs(td/2))/td;
              DistScaleFactor = Clip3(-1024, 1023, (tx*tb + 32 )>>6);
              
              img->wbp_weight[1][i][j][comp] = DistScaleFactor >> 2;
              img->wbp_weight[0][i][j][comp] = 64 - img->wbp_weight[1][i][j][comp];
              if (img->wbp_weight[1][i][j][comp] < -64 || img->wbp_weight[1][i][j][comp] > 128)
              {
                img->wbp_weight[0][i][j][comp] = 32;
                img->wbp_weight[1][i][j][comp] = 32;
                img->wp_offset[0][i][comp] = 0;
                img->wp_offset[1][i][comp] = 0;
              }
            }
          }
        }
     }
   }
 }

  if (bframe && img->MbaffFrameFlag)
  {
    for (i=0; i<2*max_fwd_ref; i++)
    {
      for (j=0; j<2*max_bwd_ref; j++)
      {
        for (comp = 0; comp<3; comp++)
        {
          for (k=2; k<6; k+=2)
          {
            img->wp_offset[k+0][i][comp] = img->wp_offset[0][i/2][comp];
            img->wp_offset[k+1][i][comp] = img->wp_offset[1][i/2][comp];

            log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
            if (active_pps->weighted_bipred_idc == 1)
            {
              img->wbp_weight[k+0][i][j][comp] =  img->wp_weight[0][i/2][comp];
              img->wbp_weight[k+1][i][j][comp] =  img->wp_weight[1][j/2][comp];
            }
            else if (active_pps->weighted_bipred_idc == 2)
            {
              td = Clip3(-128,127,listX[k+LIST_1][j]->poc - listX[k+LIST_0][i]->poc);
              if (td == 0 || listX[k+LIST_1][j]->is_long_term || listX[k+LIST_0][i]->is_long_term)
              {
                img->wbp_weight[k+0][i][j][comp] =   32;
                img->wbp_weight[k+1][i][j][comp] =   32;
              }
              else
              {
                tb = Clip3(-128,127,((k==2)?img->toppoc:img->bottompoc) - listX[k+LIST_0][i]->poc);
                
                tx = (16384 + abs(td/2))/td;
                DistScaleFactor = Clip3(-1024, 1023, (tx*tb + 32 )>>6);

                img->wbp_weight[k+1][i][j][comp] = DistScaleFactor >> 2;
                img->wbp_weight[k+0][i][j][comp] = 64 - img->wbp_weight[k+1][i][j][comp];
                if (img->wbp_weight[k+1][i][j][comp] < -64 || img->wbp_weight[k+1][i][j][comp] > 128)
                {
                  img->wbp_weight[k+1][i][j][comp] = 32;
                  img->wbp_weight[k+0][i][j][comp] = 32;
                  img->wp_offset[k+0][i][comp] = 0;
                  img->wp_offset[k+1][i][comp] = 0;
                }
              }
            }
          }
        }
      }
    }
  }
}

void partial_store_motion_into_buffer()
{
	//extern int img_height;
	//extern int img_width;

	int bx, by;
	
	int b_height = (int)iheight / 4;
	int b_width = (int)iwidth / 4;

	//Initialize the updated forward motion vectors
	for(by=0; by<b_height; by++)
	{
		for(bx=0; bx<b_width; bx++)
		{
			upd_mv[by][bx][0] = 0;
			upd_mv[by][bx][1] = 0;
		}
	}
	
	//Update forward motion vectors
	for(by=0; by<b_height; by++)
	{
		for(bx=0; bx<b_width; bx++)
			partial_update_motion_flow(bx,by);
	}

}

void partial_update_motion_flow(int bx, int by)
{
	extern double ***upd_mv;
	extern double ***dec_mv;

	int neighbor[4][3];
	double block_cont[4];
	double ori[4][2], c[4][2];
	int block_x=0;
	int block_y=0;

	double block_area = 16;
	int sgn_x;
	int sgn_y;
	int mid_line_x, mid_line_y, mid_line_bx, mid_line_by;
	int check_4multi_x, check_4multi_y;
	double cont_ratio;
	int intra_period, curr_frame_type;
	//int iframe_num;

	intra_period = INTRA_PERIOD;
	curr_frame_type = motion_picture_num%intra_period;
	//iframe_num = motion_picture_num - curr_frame_type;

	if(curr_frame_type==0)		//the last P-frame in one GOP
	{
		if(motion_picture_num!=0)
		{
			upd_mv[by][bx][0] = (double)(-dec_mv[by][bx][0]/4);		//from previous motion vector
			upd_mv[by][bx][1] = (double)(-dec_mv[by][bx][1]/4);
		}	
	}
	else					//otherwise
	{
		if(dec_mv[by][bx][0]>0)
			sgn_x = 1;
		else if(dec_mv[by][bx][0] == 0)
			sgn_x = 0;
		else
			sgn_x = -1;
		if(dec_mv[by][bx][1]>0)
			sgn_y = 1;
		else if(dec_mv[by][bx][1] == 0)
			sgn_y = 0;
		else
			sgn_y = -1;

		//////////////////////////////////////////////////////
		// Get Neighbor blocks and its Contributions
		//////////////////////////////////////////////////////

		//Get four corners of an original block
		ori[0][0] = bx*4;
		ori[0][1] = by*4;
		ori[1][0] = (bx+1)*4;
		ori[1][1] = by*4;
		ori[2][0] = bx*4;
		ori[2][1] = (by+1)*4;
		ori[3][0] = (bx+1)*4;
		ori[3][1] = (by+1)*4;

		//Get four corners of a predicted block
		c[0][0] = ori[0][0] + (double)dec_mv[by][bx][0]/4;
		c[0][1] = ori[0][1] + (double)dec_mv[by][bx][1]/4;
		c[1][0] = ori[1][0] + (double)dec_mv[by][bx][0]/4;
		c[1][1] = ori[1][1] + (double)dec_mv[by][bx][1]/4;
		c[2][0] = ori[2][0] + (double)dec_mv[by][bx][0]/4;
		c[2][1] = ori[2][1] + (double)dec_mv[by][bx][1]/4;
		c[3][0] = ori[3][0] + (double)dec_mv[by][bx][0]/4;
		c[3][1] = ori[3][1] + (double)dec_mv[by][bx][1]/4;
		
		//Check which blocks are overlapped with a predicted block
		mid_line_bx = (int)((int)c[1][0]/4);
		mid_line_by = (int)((int)c[2][1]/4);
		mid_line_x = mid_line_bx*4;
		mid_line_y = mid_line_by*4;

		check_4multi_x = (int)((int)c[0][0]) % 4;
		check_4multi_y = (int)((int)c[0][1]) % 4;

		if(check_4multi_x == 0)
		{
			if(check_4multi_y != 0)
			{
				neighbor[0][0] = (mid_line_bx - 1);
				neighbor[0][1] = (mid_line_by - 1);
				neighbor[1][0] = (mid_line_bx - 1);
				neighbor[1][1] = (mid_line_by);
				neighbor[2][0] = -1;
				neighbor[2][1] = -1;
				neighbor[3][0] = -1;
				neighbor[3][1] = -1;

				block_cont[0] = (double)(4.0 * fabs(mid_line_y - c[0][1]));
				block_cont[1] = (double)(4.0 * fabs(mid_line_y - c[2][1]));
				block_cont[2] = 0;
				block_cont[3] = 0;
			}
			else
			{
				neighbor[0][0] = (mid_line_bx - 1);
				neighbor[0][1] = (mid_line_by - 1);
				neighbor[1][0] = -1;
				neighbor[1][1] = -1;
				neighbor[2][0] = -1;
				neighbor[2][1] = -1;
				neighbor[3][0] = -1;
				neighbor[3][1] = -1;

				block_cont[0] = 16;
				block_cont[1] = 0;
				block_cont[2] = 0;
				block_cont[3] = 0;
			}
		}	
		else if(check_4multi_y == 0)
		{
			neighbor[0][0] = (mid_line_bx - 1);
			neighbor[0][1] = (mid_line_by - 1);
			neighbor[1][0] = (mid_line_bx);
			neighbor[1][1] = (mid_line_by - 1);
			neighbor[2][0] = -1;
			neighbor[2][1] = -1;
			neighbor[3][0] = -1;
			neighbor[3][1] = -1;

			block_cont[0] = (double)(fabs(mid_line_x - c[0][0]) * 4);
			block_cont[1] = (double)(fabs(mid_line_x - c[1][0]) * 4);
			block_cont[2] = 0;
			block_cont[3] = 0;
		}
		else if((mid_line_x > c[0][0])&&(mid_line_y > c[0][1]))
		{
			neighbor[0][0] = (mid_line_bx - 1);
			neighbor[0][1] = (mid_line_by - 1);
			neighbor[1][0] = (mid_line_bx);
			neighbor[1][1] = (mid_line_by - 1);
			neighbor[2][0] = (mid_line_bx - 1);
			neighbor[2][1] = (mid_line_by);
			neighbor[3][0] = (mid_line_bx);
			neighbor[3][1] = (mid_line_by);

			block_cont[0] = (double)fabs(mid_line_x - c[0][0]) * fabs(mid_line_y - c[0][1]);
			block_cont[1] = (double)fabs(mid_line_x - c[1][0]) * fabs(mid_line_y - c[1][1]);
			block_cont[2] = (double)fabs(mid_line_x - c[2][0]) * fabs(mid_line_y - c[2][1]);
			block_cont[3] = (double)fabs(mid_line_x - c[3][0]) * fabs(mid_line_y - c[3][1]);
		}
		else
		{
			//AfxMessageBox("Neighbor blocks are not detected for a motion vector");
		}

		//Availability
		neighbor[0][2] = block_availability(neighbor[0][0],neighbor[0][1]);
		neighbor[1][2] = block_availability(neighbor[1][0],neighbor[1][1]);
		neighbor[2][2] = block_availability(neighbor[2][0],neighbor[2][1]);
		neighbor[3][2] = block_availability(neighbor[3][0],neighbor[3][1]);


		/////////////////////////////////////////////////////////
		// Update Neighbor Motion Vectors
		/////////////////////////////////////////////////////////

		if(neighbor[0][2])
		{
			block_x = neighbor[0][0];
			block_y = neighbor[0][1];
			cont_ratio = (double)block_cont[0] / block_area;
			upd_mv[block_y][block_x][0] += (double)(cont_ratio) * (double)(-sgn_x) * (fabs((double)dec_mv[by][bx][0])/4);
			upd_mv[block_y][block_x][1] += (double)(cont_ratio) * (double)(-sgn_y) * (fabs((double)dec_mv[by][bx][1])/4);
		}
		if(neighbor[1][2])
		{
			block_x = neighbor[1][0];
			block_y = neighbor[1][1];
			cont_ratio = (double)block_cont[1] / block_area;
			upd_mv[block_y][block_x][0] += (double)(cont_ratio) * (double)(-sgn_x) * (fabs((double)dec_mv[by][bx][0])/4);
			upd_mv[block_y][block_x][1] += (double)(cont_ratio) * (double)(-sgn_y) * (fabs((double)dec_mv[by][bx][1])/4);
		}
		if(neighbor[2][2])
		{
			block_x = neighbor[2][0];
			block_y = neighbor[2][1];
			cont_ratio = (double)block_cont[2] / block_area;
			upd_mv[block_y][block_x][0] += (double)(cont_ratio) * (double)(-sgn_x) * (fabs((double)dec_mv[by][bx][0])/4);
			upd_mv[block_y][block_x][1] += (double)(cont_ratio) * (double)(-sgn_y) * (fabs((double)dec_mv[by][bx][1])/4);
		}
		if(neighbor[3][2])
		{
			block_x = neighbor[3][0];
			block_y = neighbor[3][1];
			cont_ratio = (double)block_cont[3] / block_area;
			upd_mv[block_y][block_x][0] += (double)(cont_ratio) * (double)(-sgn_x) * (fabs((double)dec_mv[by][bx][0])/4);
			upd_mv[block_y][block_x][1] += (double)(cont_ratio) * (double)(-sgn_y) * (fabs((double)dec_mv[by][bx][1])/4);
		}
	}
			
}

int block_availability(int bx, int by)
{
	//extern int img_height;
	//extern int img_width;

	int b4_height = (int)iheight / 4;
	int b4_width = (int)iwidth / 4;

	if((bx>=0)&&(bx<b4_width)&&(by>=0)&&(by<b4_height))
	{
		return 1;
	}
	else
	{
		return 0;
	}
	
}

void write_point(int thickness, int R, int G, int B, unsigned char *Y, unsigned char *U, unsigned char *V, int px, int py)
{
	int xi,yi;

	int inpY = R*0.299 + G*0.587 + B*0.114;
	int inpU = B*0.5 + 128 - R*0.169 - G*0.332;
	int inpV = R*0.5 + 128 - G*0.419 - B*0.0813;

	int half_thickness = (int)(thickness/2);

	for(yi=py-half_thickness; yi<py+half_thickness; yi++)
	{		
		for(xi=px-half_thickness; xi<px+half_thickness; xi++)
		{
			if((xi < 0)||(yi < 0))
				continue;
			Y[yi*iwidth+xi] = inpY;
			U[yi*iwidth/4+xi/2] = inpU;
			V[yi*iwidth/4+xi/2] = inpV;
		}
	}
}

void write_line(int thickness, int R, int G, int B, unsigned char *Y, unsigned char *U, unsigned char *V, int px1, int py1, int px2, int py2)
{
	int xi, yi;
	int half_thickness = (int)(thickness/2);
	int x1 = px1-half_thickness;
	int y1 = py1-half_thickness;
	int i;
	int dx,dy;
	int stepx, stepy;
	int fraction;

	int inpY = R*0.299 + G*0.587 + B*0.114;
	int inpU = B*0.5 + 128 - R*0.169 - G*0.332;
	int inpV = R*0.5 + 128 - G*0.419 - B*0.0813;

	if(px1==px2)
	{
		if(py1>py2)
		{
			for(i=-half_thickness; i<=half_thickness; i++)
			{
				for(yi=py2; yi<py1; yi++)
				{
					if((px1+i >= 0)&&(px1+i < iwidth)&&(yi >= 0)&&(yi <iheight))
					{
						Y[yi*iwidth+px1+i] = inpY;
						U[(yi>>1)*(iwidth>>1)+((px1+i)>>1)] = inpU;
						V[(yi>>1)*(iwidth>>1)+((px1+i)>>1)] = inpV;
					}
				}
			}
		}
		else if(py1<py2)
		{
			for(i=-half_thickness; i<=half_thickness; i++)
			{
				for(yi=py1; yi<py2; yi++)
				{
					if((px1+i >= 0)&&(px1+i < iwidth)&&(yi >= 0)&&(yi <iheight))
					{
						Y[yi*iwidth+px1+i] = inpY;
						U[(yi>>1)*(iwidth>>1)+((px1+i)>>1)] = inpU;
						V[(yi>>1)*(iwidth>>1)+((px1+i)>>1)] = inpV;
					}
				}
			}
		}
	}
	else if(py1==py2)
	{
		if(px1>px2)
		{
			for(i=-half_thickness; i<=half_thickness; i++)
			{
				for(xi=px2; xi<px1; xi++)
				{
					if((py1+i >= 0)&&(py1+i < iheight)&&(xi >= 0)&&(xi <iwidth))
					{
						Y[(py1+i)*iwidth+xi] = inpY;
						U[((py1+i)>>1)*(iwidth>>1)+(xi>>1)] = inpU;
						V[((py1+i)>>1)*(iwidth>>1)+(xi>>1)] = inpV;
					}
				}
			}
		}
		else if(px1<px2)
		{
			for(i=-half_thickness; i<=half_thickness; i++)
			{
				for(xi=px1; xi<px2; xi++)
				{
					if((py1+i >= 0)&&(py1+i < iheight)&&(xi >= 0)&&(xi <iwidth))
					{
						Y[(py1+i)*iwidth+xi] = inpY;
						U[((py1+i)>>1)*(iwidth>>1)+(xi>>1)] = inpU;
						V[((py1+i)>>1)*(iwidth>>1)+(xi>>1)] = inpV;
					}
				}
			}
		}
	}
	else		//Fast Bresenham Algorithm
	{
		dy = py2 - py1; 
		dx = px2 - px1; 
		
		if (dy < 0) {dy = -dy; stepy = -1;} else {stepy = 1;} 
		if (dx < 0) {dx = -dx; stepx = -1;} else {stepx = 1;} 
		dy <<= 1;                                      
		dx <<= 1;                                      
		Y[py1*iwidth+px1] = inpY;
		U[(py1>>1)*(iwidth>>1)+(px1>>1)] = inpU;  
		V[(py1>>1)*(iwidth>>1)+(px1>>1)] = inpV;  
		if (dx > dy)
		{ 
			fraction = dy - (dx >> 1);     
			while (px1 != px2)
			{ 
				if (fraction >= 0)
				{ 
					py1 += stepy; 
					fraction -= dx;                 
				} 
				px1 += stepx; 
				fraction += dy;                     
				Y[py1*iwidth+px1] = inpY;
				U[(py1>>1)*(iwidth>>1)+(px1>>1)] = inpU;  
				V[(py1>>1)*(iwidth>>1)+(px1>>1)] = inpV;
			} 
		}
		else
		{ 
			fraction = dx - (dy >> 1); 
			while (py1 != py2) 
			{ 
				if (fraction >= 0)
				{ 
					px1 += stepx; 
					fraction -= dy; 
				} 
				py1 += stepy; 
				fraction += dx;      
				Y[py1*iwidth+px1] = inpY;
				U[(py1>>1)*(iwidth>>1)+(px1>>1)] = inpU;  
				V[(py1>>1)*(iwidth>>1)+(px1>>1)] = inpV;
			} 
		} 
	}

}

void extract_fmv()
{
	int i,j;
	int bx,by,cx,cy;
	int b_height = iheight / 4;
	int b_width = iwidth / 4;
	int cx_dest, cy_dest;

	// Get a decoded frame
	for(i=0; i<iheight; i++)
	{
		for(j=0; j<iwidth; j++)
		{
			frm_omv[0][i*iwidth+j] = curr_frm[i][j][0];
			frm_omv[1][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][1];
			frm_omv[2][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][2];
		}
	}

	for(by=0; by<b_height; by++)
	{
		for(bx=0; bx<b_width; bx++)
		{
			// Get a center point
			cx = bx*4+2;
			cy = by*4+2;				

			// Draw points in the decoded frame
			if(block_point_flag==1)
				write_point(3, 0, 0, 0, frm_omv[0], frm_omv[1], frm_omv[2], cx, cy);

			// Draw motion vectors in the decoded frame
			if(motion_picture_num>0)
			{
				//if(BTex[by*b_width+bx]<=th_textured)		//high textured area
				//	write_line(1, 0, 0, 0, frm_omv[0], frm_omv[1], frm_omv[2], cx, cy, cx_dest, cy_dest);

				if(((dec_mv[by][bx][0]!=0)||(dec_mv[by][bx][1]!=0))) //&&((cx_dest==cx)&&(cy_dest==cy))
				{
					cx_dest = (int)ClipX(cx + dec_mv[by][bx][0]);
					cy_dest = (int)ClipY(cy + dec_mv[by][bx][1]);

				//	if(BTex[by*b_width+bx]<=th_textured)		//high textured area
					write_line(1, 0, 0, 0, frm_omv[0], frm_omv[1], frm_omv[2], cx, cy, cx_dest, cy_dest);
					write_point(2, 0, 0, 0, frm_omv[0], frm_omv[1], frm_omv[2], cx, cy);
				//	else
				//		write_point(3, 180, 50, 180, frm_omv[0], frm_omv[1], frm_omv[2], cx, cy);
				}
			}
		}
	}

	//Write a file
	for(i = 0; i < iheight ; i++)
		fwrite(&(frm_omv[0][i*iwidth]),1,iwidth,OMV_FILE );

	for(i = 0; i < (iheight/2); i++)
		fwrite(&(frm_omv[1][i*iwidth/2]),1,iwidth/2,OMV_FILE );

	for(i = 0; i < (iheight/2); i++)
		fwrite(&(frm_omv[2][i*iwidth/2]),1,iwidth/2,OMV_FILE );
}

void draw_mb_type()
{
	int i, j, r, c, n, m;
	int bcol, brow, b4col, b4row, px, py, bx, by, b4x, b4y;
	int max_mb_nr = (int)(iheight*iwidth/16/16);
	int intra_period = INTRA_PERIOD;	//IPPP
	int curr_frame_type = motion_picture_num%intra_period;
	int bwidth = iwidth/16;
	int width, height, TopLeft_x, TopLeft_y, BottomRight_x, BottomRight_y;
	int px_dest, py_dest;
	int thickness = 1;
	int R = 0;
	int G = 0;
	int B = 0;
	int b4width = iwidth/4;
	int b4height = iheight/4;
	int ResE, ResEnergyTh, MaxCoeffIdx;
	double mag = 0;
	double MvMagTh;
	int group_index;
	BYTE zigzag_scan[16] = {0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15};
	int group_color[10][3] = {{255,128,128},{255,255,0},{255,128,0},{0,255,0},{0,128,0},{0,255,255},
								{0,0,255},{255,0,255},{128,128,0},{128,0,255}};

	int ShowMbTypeBlock = 1;
	int ShowTexturedBlock = 1;
	int ShowMotionBlock = 0;
	int ShowTexturedMotionBlock = 0;
	
	// Get a decoded frame
	for(i=0; i<iheight; i++)
	{
		for(j=0; j<iwidth; j++)
		{
			frm_mbt[0][i*iwidth+j] = curr_frm[i][j][0];
			frm_mbt[1][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][1];
			frm_mbt[2][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][2];
		}
	}

	for(i=0; i<max_mb_nr; i++)
	{
		bcol = i%bwidth;
		brow = (int)(i/bwidth);
		px = bcol*16;
		py = brow*16;

		if(curr_frame_type==0)		//I Frames
		{
			if((MbType[i]==I4MB)&&(DrawMbTypeEnb[9]==1))	//I4x4
			{
				for(r=0; r<16; r+=4)
					for(c=0; c<16; c+=4)
					{
						select_blocks(0, px+c, py+r, ClipX(px+c+4), ClipY(py+r+4));
						draw_rectangle(ShowMbTypeBlock, 1, 30, 130, 130, frm_mbt[0], frm_mbt[1], frm_mbt[2], px+c, py+r, ClipX(px+c+4), ClipY(py+r+4));
					}
			}
			else if((MbType[i]==IPCM)&&(DrawMbTypeEnb[14]==1))	//IPCM
			{
				select_blocks(0, px, py, ClipX(px+16), ClipY(py+16));
				draw_rectangle(ShowMbTypeBlock, 1, 0, 150, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
			}
			else if((MbType[i]==I16MB)&&(DrawMbTypeEnb[10]==1))	//I16x16
			{
				select_blocks(0,px, py, ClipX(px+16), ClipY(py+16));
				draw_rectangle(ShowMbTypeBlock, 1, 30, 130, 130, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
			}
			else	//I16x16
			{
				//select_blocks(0, px, py, ClipX(px+16), ClipY(py+16));
				//draw_rectangle(ShowMbTypeBlock, 1, 30, 130, 130, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
			}
		}
		else						//P Frames
		{
			if((MbType[i]==0)&&(DrawMbTypeEnb[0]==1))	//SKIP
			{
				select_blocks(0, px, py, ClipX(px+16), ClipY(py+16));
				draw_rectangle(ShowMbTypeBlock, 1, 180, 50, 180, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
			}
			else if((MbType[i]==1)&&(DrawMbTypeEnb[1]==1))	//P16x16
			{
				select_blocks(0, px, py, ClipX(px+16), ClipY(py+16));
				draw_rectangle(ShowMbTypeBlock, 1, 0, 0, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));

				if((segmentation_flag==1)&&(IS_ASSIGNED(mblabel[motion_picture_num][i])==1)&&(label_size[mblabel[motion_picture_num][i]]>1))
				{
					group_index = display_group_idx[mblabel[motion_picture_num][i]];
					draw_rectangle(ShowMbTypeBlock, 1, group_color[group_index][0], group_color[group_index][1], group_color[group_index][2], frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
				}
			}
			else if((MbType[i]==2)&&(DrawMbTypeEnb[2]==1))	//P16x8
			{
				for(r=0; r<16; r+=8)
					for(c=0; c<16; c+=16)
					{
						select_blocks(0, px+c, py+r, ClipX(px+c+16), ClipY(py+r+8));
						draw_rectangle(ShowMbTypeBlock, 1, 0, 0, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], px+c, py+r, ClipX(px+c+16), ClipY(py+r+8));
					}

				if((segmentation_flag==1)&&(IS_ASSIGNED(mblabel[motion_picture_num][i])==1)&&(label_size[mblabel[motion_picture_num][i]]>1))
				{
					group_index = display_group_idx[mblabel[motion_picture_num][i]];
					draw_rectangle(ShowMbTypeBlock, 1, group_color[group_index][0], group_color[group_index][1], group_color[group_index][2], frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
				}
			}
			else if((MbType[i]==3)&&(DrawMbTypeEnb[3]==1))	//P8x16
			{
				for(r=0; r<16; r+=16)
					for(c=0; c<16; c+=8)
					{
						select_blocks(0, px+c, py+r, ClipX(px+c+8), ClipY(py+r+16));
						draw_rectangle(ShowMbTypeBlock, 1, 0, 0, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], px+c, py+r, ClipX(px+c+8), ClipY(py+r+16));
					}

				if((segmentation_flag==1)&&(IS_ASSIGNED(mblabel[motion_picture_num][i])==1)&&(label_size[mblabel[motion_picture_num][i]]>1))
				{
					group_index = display_group_idx[mblabel[motion_picture_num][i]];
					draw_rectangle(ShowMbTypeBlock, 1, group_color[group_index][0], group_color[group_index][1], group_color[group_index][2], frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
				}
			}
			else if(MbType[i]==P8x8)	
			{
				if(DrawMbTypeEnb[8]==1)
				{
					for(r=0; r<16; r+=8)
						for(c=0; c<16; c+=8)
						{
							select_blocks(0, px+c, py+r, ClipX(px+c+8), ClipY(py+r+8));
							draw_rectangle(ShowMbTypeBlock, 1, 0, 0, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], px+c, py+r, ClipX(px+c+8), ClipY(py+r+8));
						}
				}

				for(j=0; j<4; j++)		//for each 8x8 block
				{
					b4col = j%2;
					b4row = (int)(j/2);
					bx = px + b4col*8;
					by = py + b4row*8;

					if((B8Mode[i][j]==4)&&(DrawMbTypeEnb[4]==1))	//P8x8
					{
						select_blocks(0, bx, by, ClipX(bx+8), ClipY(by+8));
						draw_rectangle(ShowMbTypeBlock, 1, 0, 0, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], bx, by, ClipX(bx+8), ClipY(by+8));
					}
					else if((B8Mode[i][j]==5)&&(DrawMbTypeEnb[5]==1))	//P8x4
					{
						for(r=0; r<8; r+=4)
							for(c=0; c<8; c+=8)
							{
								select_blocks(0, bx+c, by+r, ClipX(bx+c+8), ClipY(by+r+4));
								draw_rectangle(ShowMbTypeBlock, 1, 0, 0, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], bx+c, by+r, ClipX(bx+c+8), ClipY(by+r+4));
							}
					}
					else if((B8Mode[i][j]==6)&&(DrawMbTypeEnb[6]==1))	//P4x8
					{
						for(r=0; r<8; r+=8)
							for(c=0; c<8; c+=4)
							{
								select_blocks(0, bx+c, by+r, ClipX(bx+c+4), ClipY(by+r+8));
								draw_rectangle(ShowMbTypeBlock, 1, 0, 0, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], bx+c, by+r, ClipX(bx+c+4), ClipY(by+r+8));
							}
					}
					else if((B8Mode[i][j]==7)&&(DrawMbTypeEnb[7]==1))	//P4x4
					{
						for(r=0; r<8; r+=4)
							for(c=0; c<8; c+=4)
							{
								select_blocks(0, bx+c, by+r, ClipX(bx+c+4), ClipY(by+r+4));
								draw_rectangle(ShowMbTypeBlock, 1, 0, 0, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], bx+c, by+r, ClipX(bx+c+4), ClipY(by+r+4));
							}
					}
				}

				if((segmentation_flag==1)&&(IS_ASSIGNED(mblabel[motion_picture_num][i])==1)&&(label_size[mblabel[motion_picture_num][i]]>1))
				{
					group_index = display_group_idx[mblabel[motion_picture_num][i]];
					draw_rectangle(ShowMbTypeBlock, 1, group_color[group_index][0], group_color[group_index][1], group_color[group_index][2], frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
				}
			}
			else if((MbType[i]==I4MB)&&(DrawMbTypeEnb[9]==1))
			{
				for(r=0; r<16; r+=4)
					for(c=0; c<16; c+=4)
					{
						select_blocks(0, px+c, py+r, ClipX(px+c+4), ClipY(py+r+4));
						draw_rectangle(ShowMbTypeBlock, 1, 30, 130, 130, frm_mbt[0], frm_mbt[1], frm_mbt[2], px+c, py+r, ClipX(px+c+4), ClipY(py+r+4));
					}
			}
			else if((MbType[i]==IPCM)&&(DrawMbTypeEnb[14]==1))
			{
				select_blocks(0, px, py, ClipX(px+16), ClipY(py+16));
				draw_rectangle(ShowMbTypeBlock, 1, 0, 150, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
			}
			else if((MbType[i]==I16MB)&&(DrawMbTypeEnb[10]==1))
			{
				select_blocks(0, px, py, ClipX(px+16), ClipY(py+16));
				draw_rectangle(ShowMbTypeBlock, 1, 30, 130, 130, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
			}
			else	//I16x16
			{
				//select_blocks(0, px, py, ClipX(px+16), ClipY(py+16));
				//draw_rectangle(ShowMbTypeBlock, 1, 30, 130, 130, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+16), ClipY(py+16));
			}
		}
	}

	//Draw textured blocks
	if(draw_textured_blocks==1)
	{
		//Setting and Initialization
		ResEnergyTh = 1;	//Residual Energy Threshold
		MaxCoeffIdx = 10;	//Maximum Coefficient Index in the zigzag scan

		for(by=0; by<b4height; by++)
			for(bx=0; bx<b4width; bx++)
			{
				ResE = 0;

				px = bx*4;
				py = by*4;

				//calculate the residual energy
				for(i=0; i<=MaxCoeffIdx; i++)
				{
					c = (zigzag_scan[i] & 3);
					r = ((zigzag_scan[i] >> 2) & 3);
					ResE +=(int)(curr_res[py+r][px+c][0]*curr_res[py+r][px+c][0]);
				}

				if(ResE>=ResEnergyTh)
				{
					if((multiple_block_drawing_flag==1)&&(DrawBlockFlag[0][by*b4width+bx]==1))
					{
						//calculate the mean of motion vectors in the neighborhood
						//calculate the variance of motion vectors in the neighborhood

						//select blocks which have motion vectors similar with its neighbor blocks
						select_blocks(1, px, py, ClipX(px+4), ClipY(py+4));
						draw_rectangle(ShowTexturedBlock, 1, 150, 150, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+4), ClipY(py+4));		//금색
					}
				}
			}
	}

	//Draw motion blocks
	if(draw_motion_blocks==1)
	{
		//Setting and Initialization
		MvMagTh = 4;	//Motion Vector Magnitude Threshold

		for(by=0; by<b4height; by++)
			for(bx=0; bx<b4width; bx++)
			{
				px = bx*4;
				py = by*4;

				//calculate the squared magnitude of the corresponding motion vector
				mag = dec_mv[by][bx][0]*dec_mv[by][bx][0] + dec_mv[by][bx][1]*dec_mv[by][bx][1];

				if(mag>MvMagTh)
				{
					if((multiple_block_drawing_flag==1)&&(DrawBlockFlag[0][by*b4width+bx]==1))
					{
						select_blocks(2, px, py, ClipX(px+4), ClipY(py+4));
						draw_rectangle(ShowMotionBlock, 1, 250, 100, 50, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+4), ClipY(py+4));		//주황
					}
				}
			}
	}

	//Draw textured motion blocks
	if(draw_textured_motion_blocks==1)
	{
		for(by=0; by<b4height; by++)
			for(bx=0; bx<b4width; bx++)
			{
				px = bx*4;
				py = by*4;

				if((DrawBlockFlag[0][by*b4width+bx]==1)&&(DrawBlockFlag[1][by*b4width+bx]==1)&&(DrawBlockFlag[2][by*b4width+bx]==1))
				{
					select_blocks(3, px, py, ClipX(px+4), ClipY(py+4));
					draw_rectangle(ShowTexturedMotionBlock, 1, 200, 0, 0, frm_mbt[0], frm_mbt[1], frm_mbt[2], px, py, ClipX(px+4), ClipY(py+4));		//빨강
				}
			}
	}

	//Write a file
	for(i = 0; i < iheight ; i++)
		fwrite(&(frm_mbt[0][i*iwidth]),1,iwidth,MBT_FILE );

	for(i = 0; i < (iheight/2); i++)
		fwrite(&(frm_mbt[1][i*iwidth/2]),1,iwidth/2,MBT_FILE );

	for(i = 0; i < (iheight/2); i++)
		fwrite(&(frm_mbt[2][i*iwidth/2]),1,iwidth/2,MBT_FILE );
}

void draw_rectangle(int ShowFlag, int thickness, int R, int G, int B, unsigned char *Y, unsigned char *U, unsigned char *V, int TopLeft_x, int TopLeft_y, int BottomRight_x, int BottomRight_y)
{
	int width = BottomRight_x - TopLeft_x;
	int height = BottomRight_y - TopLeft_y;
	
	if(ShowFlag==1)
	{
		write_line(thickness, R, G, B, Y, U, V, TopLeft_x, TopLeft_y, TopLeft_x+width, TopLeft_y);
		write_line(thickness, R, G, B, Y, U, V, TopLeft_x+width, TopLeft_y, TopLeft_x+width, TopLeft_y+height);
		write_line(thickness, R, G, B, Y, U, V, TopLeft_x+width, TopLeft_y+height, TopLeft_x, TopLeft_y+height);
		write_line(thickness, R, G, B, Y, U, V, TopLeft_x, TopLeft_y+height, TopLeft_x, TopLeft_y);
	}
}

void select_blocks(int DrawBlockLevel, int TopLeft_x, int TopLeft_y, int BottomRight_x, int BottomRight_y)
{
	int r,c;
	int b4x, b4y;
	int b4width = iwidth/4;

	if(multiple_block_drawing_flag==1)
	{
		for(r=TopLeft_y; r<BottomRight_y; r+=4)
			for(c=TopLeft_x; c<BottomRight_x; c+=4)
			{
				b4x = (int)(c/4);
				b4y = (int)(r/4);
				DrawBlockFlag[DrawBlockLevel][b4y*b4width+b4x] = 1;
			}
	}
}

void motion_segment()
{
	int k,h,i,j,n,m;
	int mbx, mby;
	int intra_period = INTRA_PERIOD;	//IPPP
	int curr_frame_type = motion_picture_num%intra_period;
	int mbwidth = iwidth/16;
	int mbheight = iheight/16;
	int max_mb_nr = (int)(iheight*iwidth/16/16);
	int max_label_idx = -1;
	int flip_flag, label_change_flag, curr_mb_label, max_neighbor_label_size, max_group_num;
	int neighbor_block_idx[8], neighbor_block_label[8], best_neighbor_label;
	int recurrent_cnt;
	int px,py,start_px,start_py;
	int ResGroupON[UNKNOWN_LABEL+1];
	int overlap_check_flag, curr_buf_address, pros_buf_address, new_buf_address, overlap_mb_cnt;
	//int inv_display_group_idx[UNKNOWN_LABEL+1];
	int TmpX, TmpY, mb_cnt;
	int max_mbx, min_mbx, max_mby, min_mby;
	double total_occur_prob;
	int prev_dist;
	int picture_num, frame_type, cnt;
	int tmp_gbuf_Status[MAX_BUF_NUM], tmp_gbuf_FrmGroupIdx[MAX_BUF_NUM];
	int prev_max_buf_group_num;
	int overlap_buf_group_idx, assigned_buf_idx[UNKNOWN_LABEL+1];
	int mb_label;
	int min_x, max_x, min_y, max_y, min_obj_x, max_obj_x, min_obj_y, max_obj_y;
	int tol = BKD_STR_TOL;	//배경제거를 위한 파라메터
	int YSize = iheight*iwidth;
	int Voffset = YSize + YSize/4;
	double delta_x, delta_y, delta_width, delta_height;
	double mid_pos_x, mid_pos_y, mid_width, mid_height;
	int hor_decoding_margin = HOR_DECODING_MARGIN;
	int ver_decoding_margin = VER_DECODING_MARGIN;
	int th_delta_width = DELTA_WIDTH_MARGIN;
	int th_delta_height = DELTA_HEIGHT_MARGIN;
	int diff_width, diff_height;
	int th_diff_width = DIFF_WIDTH_MARGIN;
	int th_diff_height = DIFF_HEIGHT_MARGIN;
	int left_top_mbx, left_top_mby, right_bottom_mbx, right_bottom_mby;
	int rejected_group_num, total_block_groups;
	double spatial_noise_rate;
	int overlap_flag[UNKNOWN_LABEL];

	int *bkd_obj;
	bkd_obj = (int*)calloc(YSize,sizeof(int));

	if(curr_frame_type==1)
		prev_dist = 2;
	else if(curr_frame_type>1)
		prev_dist = 1;

	if(motion_picture_num==659)
			i=0;

	if(curr_frame_type==0)		//I-frames
	{
		if(trajectory_optimization_flag==1)
		{
			if(motion_picture_num>=intra_period)
			{
				//----------
				//background subtraction (배경 제거)
				//----------

				/*
				for(i=0; i<iheight; i++)
				{
					for(j=0; j<iwidth; j++)
					{
						if(IsNearlyZero(curr_frm[i][j][0] - BkdBuf[0][i*iwidth+j],tol))
						{
							bkd_str[0][i*iwidth+j] = 16; //clip(DecBuf[j] - BkdBuf[j] +16);
							bkd_obj[i*iwidth+j] = 0;
						}
						else
						{
							bkd_str[0][i*iwidth+j] = curr_frm[i][j][0];
							bkd_obj[i*iwidth+j] = 1;
						}

						if(IsNearlyZero(curr_frm[i>>1][j>>1][1] - BkdBuf[0][(i>>1)*(iwidth>>1)+(j>>1)+YSize],tol))
							bkd_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = 128; //clip(DecBuf[j+YSize] - BkdBuf[j+YSize] +128);
						else
							bkd_str[1][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][1];

						if(IsNearlyZero(curr_frm[i>>1][j>>1][2] - BkdBuf[0][(i>>1)*(iwidth>>1)+(j>>1)+Voffset],tol))
							bkd_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = 128; //clip(DecBuf[j+Voffset] - BkdBuf[j+Voffset] +128);
						else
							bkd_str[2][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][2];
					}
				}
				*/

				//----------
				// Calculation of Object Position and Region in I frames (객체 위치 및 영역 계산)
				//----------

				//initialization
				for(i=0; i<YSize; i++)		bkd_obj[i] = 0;
				for(i=0; i<YSize; i++)		bkd_str[0][i] = 0;
				for(i=0; i<YSize/4; i++)	bkd_str[1][i] = 0;
				for(i=0; i<YSize/4; i++)	bkd_str[2][i] = 0;

				for(i=0; i<max_buf_group_num; i++)		//for each group of the current frame (현재 프레임의 각 그룹에 대하여)
				{
					//increase FrameCnt for the object group (객체 그룹에 대하여 FrameCnt를 증가시킨다.)
					gbuf_Status[i][intra_period-1] = gbuf[i].Status;

					if((gbuf[i].Status==OBJECT_GROUP)&&(gbuf[i].ObjectFrame+gbuf[i].FrameCnt>=intra_period))
					{
						(gbuf[i].FrameCnt)++;

						//extract the color information of objects (객체 컬러정보 추출)
						//1. background subtraction (예측된 객체 영역에서 배경 제거)
						min_x = ClipX(gbuf[i].PosX[gbuf[i].FrameCnt] - gbuf[i].Width[gbuf[i].FrameCnt]/2 -hor_decoding_margin);
						max_x = ClipX(gbuf[i].PosX[gbuf[i].FrameCnt] + gbuf[i].Width[gbuf[i].FrameCnt]/2 +hor_decoding_margin);
						min_y = ClipY(gbuf[i].PosY[gbuf[i].FrameCnt] - gbuf[i].Height[gbuf[i].FrameCnt]/2 -ver_decoding_margin);
						max_y = ClipY(gbuf[i].PosY[gbuf[i].FrameCnt] + gbuf[i].Height[gbuf[i].FrameCnt]/2 +ver_decoding_margin);

						for(m=min_y; m<=max_y; m++)			
						{
							for(n=min_x; n<=max_x; n++)
							{
								if(IsNearlyZero(curr_frm[m][n][0] - BkdBuf[0][m*iwidth+n],tol))
								{
									bkd_str[0][m*iwidth+n] = 16; //clip(DecBuf[j] - BkdBuf[j] +16);
									bkd_obj[m*iwidth+n] = 0;
								}
								else
								{
									bkd_str[0][m*iwidth+n] = curr_frm[m][n][0];
									bkd_obj[m*iwidth+n] = 1;
								}

								if(IsNearlyZero(curr_frm[m>>1][n>>1][1] - BkdBuf[0][(m>>1)*(iwidth>>1)+(n>>1)+YSize],tol))
									bkd_str[1][(m>>1)*(iwidth>>1)+(n>>1)] = 128; //clip(DecBuf[j+YSize] - BkdBuf[j+YSize] +128);
								else
									bkd_str[1][(m>>1)*(iwidth>>1)+(n>>1)] = curr_frm[m>>1][n>>1][1];

								if(IsNearlyZero(curr_frm[m>>1][n>>1][2] - BkdBuf[0][(m>>1)*(iwidth>>1)+(n>>1)+Voffset],tol))
									bkd_str[2][(m>>1)*(iwidth>>1)+(n>>1)] = 128; //clip(DecBuf[j+Voffset] - BkdBuf[j+Voffset] +128);
								else
									bkd_str[2][(m>>1)*(iwidth>>1)+(n>>1)] = curr_frm[m>>1][n>>1][2];
							}
						}					
						
						//2. region optimization (객체 픽셀 위치의 최대, 최소값 추출 후 객체 영역 최적화)
						max_obj_x = min_x-1;
						min_obj_x = max_x+1;
						max_obj_y = min_y-1;
						min_obj_y = max_y+1;
						for(m=min_y; m<=max_y; m++)			
						{
							for(n=min_x; n<=max_x; n++)
							{
								if(bkd_obj[m*iwidth+n]==1)
								{
									if(max_obj_x<n)
										max_obj_x = n;
									if(min_obj_x>n)
										min_obj_x = n;
									if(max_obj_y<m)
										max_obj_y = m;
									if(min_obj_y>m)
										min_obj_y = m;
								}
							}
						}

						// calculation of object location and region in the current I-frame (현재 I-frame에서의 객체 위치 및 영역 계산)
						if((min_obj_x<max_obj_x)&&(min_obj_y<max_obj_y))
						{
							gbuf[i].PosX[gbuf[i].FrameCnt] = (int)(min_obj_x+(max_obj_x - min_obj_x)/2);
							gbuf[i].PosY[gbuf[i].FrameCnt] = (int)(min_obj_y+(max_obj_y - min_obj_y)/2);
							gbuf[i].Width[gbuf[i].FrameCnt] = (int)(max_obj_x - min_obj_x);
							gbuf[i].Height[gbuf[i].FrameCnt] = (int)(max_obj_y - min_obj_y);
						}

						// detecting the disappearance of objects (객체의 사라짐 관찰)
						left_top_mbx = (int)(ClipX(gbuf[i].PosX[gbuf[i].FrameCnt] - gbuf[i].Width[gbuf[i].FrameCnt]/2)/16);								
						left_top_mby = (int)(ClipY(gbuf[i].PosY[gbuf[i].FrameCnt] - gbuf[i].Height[gbuf[i].FrameCnt]/2)/16);
						right_bottom_mbx = (int)(ClipX(gbuf[i].PosX[gbuf[i].FrameCnt] + gbuf[i].Width[gbuf[i].FrameCnt]/2)/16);
						right_bottom_mby = (int)(ClipY(gbuf[i].PosY[gbuf[i].FrameCnt] + gbuf[i].Height[gbuf[i].FrameCnt]/2)/16);

						if((left_top_mbx==right_bottom_mbx)&&(left_top_mby==right_bottom_mby))
							gbuf[i].Status = HIDDEN_GROUP;

						//motion interpolation in P frmaes of the previous GOP (이전 GOP에서의 객체 위치 및 영역 계산)
						if(gbuf[i].FrameCnt>=intra_period)
						{
							delta_x = ((double)gbuf[i].PosX[gbuf[i].FrameCnt-intra_period] - (double)gbuf[i].PosX[gbuf[i].FrameCnt])/(double)intra_period;
							delta_y = ((double)gbuf[i].PosY[gbuf[i].FrameCnt-intra_period] - (double)gbuf[i].PosY[gbuf[i].FrameCnt])/(double)intra_period;
							delta_width = ((double)gbuf[i].Width[gbuf[i].FrameCnt-intra_period] - (double)gbuf[i].Width[gbuf[i].FrameCnt])/(double)intra_period;
							delta_height = ((double)gbuf[i].Height[gbuf[i].FrameCnt-intra_period] - (double)gbuf[i].Height[gbuf[i].FrameCnt])/(double)intra_period;

							/*
							if((delta_width>th_delta_width)||(delta_width<-th_delta_width))
							{
								//gbuf[i].Width[gbuf[i].FrameCnt] = gbuf[i].Width[gbuf[i].FrameCnt-intra_period] - intra_period*th_delta_width;
								//delta_width = th_delta_width;
								gbuf[i].Width[gbuf[i].FrameCnt] = gbuf[i].Width[gbuf[i].FrameCnt-intra_period];
								delta_width = 0;
							}
							//else if(delta_width<-th_delta_width)
							//{
							//	gbuf[i].Width[gbuf[i].FrameCnt] = gbuf[i].Width[gbuf[i].FrameCnt-intra_period] + intra_period*th_delta_width;
							//	delta_width = -th_delta_width;
							//}
							if((delta_height>th_delta_height)||(delta_height<-th_delta_height))
							{
								//gbuf[i].Height[gbuf[i].FrameCnt] = gbuf[i].Height[gbuf[i].FrameCnt-intra_period] - intra_period*th_delta_height;
								//delta_height = th_delta_height;
								gbuf[i].Height[gbuf[i].FrameCnt] = gbuf[i].Height[gbuf[i].FrameCnt-intra_period];
								delta_height = 0;
							}
							//else if(delta_height<-th_delta_height)
							//{
							//	gbuf[i].Height[gbuf[i].FrameCnt] = gbuf[i].Height[gbuf[i].FrameCnt-intra_period] + intra_period*th_delta_height;
							//	delta_height = -th_delta_height;
							//}						
							*/

							mid_pos_x = (double)gbuf[i].PosX[gbuf[i].FrameCnt];
							mid_pos_y = (double)gbuf[i].PosY[gbuf[i].FrameCnt];
							mid_width = (double)gbuf[i].Width[gbuf[i].FrameCnt];
							mid_height = (double)gbuf[i].Height[gbuf[i].FrameCnt];
							for(h=1; h<intra_period; h++)
							{
								mid_pos_x += delta_x;
								mid_pos_y += delta_y;
								mid_width += delta_width;
								mid_height += delta_height;

								gbuf[i].PosX[gbuf[i].FrameCnt-h] = (int)(mid_pos_x+0.5);
								gbuf[i].PosY[gbuf[i].FrameCnt-h] = (int)(mid_pos_y+0.5);
								gbuf[i].Width[gbuf[i].FrameCnt-h] = (int)(mid_width+0.5);
								gbuf[i].Height[gbuf[i].FrameCnt-h] = (int)(mid_height+0.5);							
							}
						}
					}
				}

				//drawing the motion trajectory (이전 GOP에 대한 결과 동영상 쓰기)
				if((extract_trajectory_flag==1)&&(object_recognition_flag==1))
					draw_trajectory();

				// Output trajectory data (Trajectory 데이터 출력하기)
				if((output_trajectory_flag==1)&&(time_measure_flag!=1))
				{
					for(h=intra_period-1; h>=0; h--)
					{
						for(i=0; i<max_buf_group_num; i++)		//for each group of the current frame (현재 프레임의 각 그룹에 대하여)
						{
							if(gbuf[i].Status==OBJECT_GROUP)
							{
								left_top_mbx = (int)(ClipX(gbuf[i].PosX[gbuf[i].FrameCnt-h] - gbuf[i].Width[gbuf[i].FrameCnt-h]/2)/16);								
								left_top_mby = (int)(ClipY(gbuf[i].PosY[gbuf[i].FrameCnt-h] - gbuf[i].Height[gbuf[i].FrameCnt-h]/2)/16);
								right_bottom_mbx = (int)(ClipX(gbuf[i].PosX[gbuf[i].FrameCnt-h] + gbuf[i].Width[gbuf[i].FrameCnt-h]/2)/16);
								right_bottom_mby = (int)(ClipY(gbuf[i].PosY[gbuf[i].FrameCnt-h] + gbuf[i].Height[gbuf[i].FrameCnt-h]/2)/16);

								if(left_top_mbx<right_bottom_mbx)
								{
									left_top_mbx = left_top_mbx +1;									
									right_bottom_mbx = right_bottom_mbx -1;									
								}
								if(left_top_mby<right_bottom_mby)
								{
									left_top_mby = left_top_mby +1;
									right_bottom_mby = right_bottom_mby -1;
								}								
							}
							else
							{
								left_top_mbx = 0;
								left_top_mby = 0;
								right_bottom_mbx = 0;
								right_bottom_mby = 0;							
							}
							fprintf(TRA_FILE,"%d	%d	%d	%d	",left_top_mbx,left_top_mby,right_bottom_mbx,right_bottom_mby);
						}
						fprintf(TRA_FILE,"\n");
					}
				}

				if(time_measure_flag==0)
				{
					//Writing the background subtracted images
					for(i = 0; i < iheight ; i++)
						fwrite(&(bkd_str[0][i*iwidth]),1,iwidth,BKS_FILE );
					for(i = 0; i < (iheight/2); i++)
						fwrite(&(bkd_str[1][i*iwidth/2]),1,iwidth/2,BKS_FILE );
					for(i = 0; i < (iheight/2); i++)
						fwrite(&(bkd_str[2][i*iwidth/2]),1,iwidth/2,BKS_FILE );
				}
			}
		}
	}
	else		//P-frames
	{
		//initialization
		for(i=0; i<UNKNOWN_LABEL+1; i++)
			label_size[i] = 0;

		//----------
		//select inter blocks: delete skip blocks (인터 매크로블록의 추출) 
		//----------
		for(mby=0; mby<mbheight; mby++)
			for(mbx=0; mbx<mbwidth; mbx++)
			{
				k= mby*mbwidth + mbx;

				//if((MbType[k]==1)||(MbType[k]==2)||(MbType[k]==3)||(MbType[k]==P8x8))
				if(MbType[k]>0)
					mblabel[motion_picture_num][k] = UNKNOWN_LABEL;
				else
					mblabel[motion_picture_num][k] = -1;
			}

		//----------
		//macroblock grouping: block group labeling (매크로블록의 그룹핑)
		//----------

		rejected_group_num = 0;		
		recurrent_cnt = 0;
		do
		{
			flip_flag = 0;

			for(k=0; k<max_mb_nr; k++)
			{
				if(mblabel[motion_picture_num][k]>=0)	//no SKIP block
				{
					mbx = k%mbwidth;
					mby = (int)(k/mbwidth);

					curr_mb_label = mblabel[motion_picture_num][k];

					if(curr_mb_label==UNKNOWN_LABEL)
						max_neighbor_label_size = 0;
					else if(curr_mb_label<UNKNOWN_LABEL)
						max_neighbor_label_size = label_size[curr_mb_label];

					if(mbx+1>=mbwidth)
						neighbor_block_label[0] = -1;
					else
					{
						neighbor_block_idx[0] = mby*mbwidth+(mbx+1);
						neighbor_block_label[0] = mblabel[motion_picture_num][neighbor_block_idx[0]];
					}
					if((mbx+1>=mbwidth)||(mby+1>=mbheight))
						neighbor_block_label[1] = -1;
					else
					{
						neighbor_block_idx[1] = (mby+1)*mbwidth+(mbx+1);
						neighbor_block_label[1] = mblabel[motion_picture_num][neighbor_block_idx[1]];
					}
					if(mby+1>=mbheight)
						neighbor_block_label[2] = -1;
					else
					{
						neighbor_block_idx[2] = (mby+1)*mbwidth+mbx;
						neighbor_block_label[2] = mblabel[motion_picture_num][neighbor_block_idx[2]];
					}
					if((mby+1>=mbheight)||(mbx-1<0))
						neighbor_block_label[3] = -1;
					else
					{
						neighbor_block_idx[3] = (mby+1)*mbwidth+(mbx-1);
						neighbor_block_label[3] = mblabel[motion_picture_num][neighbor_block_idx[3]];
					}
					if(mbx-1<0)
						neighbor_block_label[4] = -1;
					else
					{
						neighbor_block_idx[4] = mby*mbwidth+(mbx-1);
						neighbor_block_label[4] = mblabel[motion_picture_num][neighbor_block_idx[4]];
					}
					if((mbx-1<0)||(mby-1<0))
						neighbor_block_label[5] = -1;
					else
					{
						neighbor_block_idx[5] = (mby-1)*mbwidth+(mbx-1);
						neighbor_block_label[5] = mblabel[motion_picture_num][neighbor_block_idx[5]];
					}
					if(mby-1<0)
						neighbor_block_label[6] = -1;
					else
					{
						neighbor_block_idx[6] = (mby-1)*mbwidth+mbx;
						neighbor_block_label[6] = mblabel[motion_picture_num][neighbor_block_idx[6]];
					}
					if((mby-1<0)||(mbx+1>=mbwidth))
						neighbor_block_label[7] = -1;
					else
					{
						neighbor_block_idx[7] = (mby-1)*mbwidth+(mbx+1);
						neighbor_block_label[7] = mblabel[motion_picture_num][neighbor_block_idx[7]];
					}

					label_change_flag = 0;
					if((neighbor_block_label[0]<0)&&(neighbor_block_label[1]<0)&&(neighbor_block_label[2]<0)&&(neighbor_block_label[3]<0)
						&&(neighbor_block_label[4]<0)&&(neighbor_block_label[5]<0)&&(neighbor_block_label[6]<0)&&(neighbor_block_label[7]<0))
					{
						mblabel[motion_picture_num][k] = -1;		//changing the label of isolated macroblocks into SKIP mode (고립된 매크로블록의 레이블을 SKIP 모드로 전환한다.)
						rejected_group_num++;
					}
					else		//group consisting of several non-skip blocks
					{
						//no neighbor blocks are assigned
						if((IS_ASSIGNED(neighbor_block_label[0])==0)&&(IS_ASSIGNED(neighbor_block_label[1])==0)&&(IS_ASSIGNED(neighbor_block_label[2])==0)&&(IS_ASSIGNED(neighbor_block_label[3])==0)
							&&(IS_ASSIGNED(neighbor_block_label[4])==0)&&(IS_ASSIGNED(neighbor_block_label[5])==0)&&(IS_ASSIGNED(neighbor_block_label[6])==0)&&(IS_ASSIGNED(neighbor_block_label[7])==0))
						{
							if(curr_mb_label==UNKNOWN_LABEL)
							{
								//new label assigned
								max_label_idx++;
								mblabel[motion_picture_num][k] = max_label_idx;
								label_size[max_label_idx] = 1;							
								flip_flag = 1;
							}
						}
						else	//there are more than one assigned neighbor block.
						{
							best_neighbor_label = curr_mb_label;

							for(n=0; n<8; n++)
							{
								if(IS_ASSIGNED(neighbor_block_label[n])==1)
								{
									if((best_neighbor_label!=neighbor_block_label[n])&&(max_neighbor_label_size<label_size[neighbor_block_label[n]]))
									{
										if(!((curr_mb_label<UNKNOWN_LABEL)&&(max_neighbor_label_size==label_size[curr_mb_label])&&(curr_mb_label==neighbor_block_label[n])))
										{
											best_neighbor_label = neighbor_block_label[n];
											max_neighbor_label_size = label_size[neighbor_block_label[n]];
										}
									}
								}
							}
							
							if(curr_mb_label!=best_neighbor_label)
							{
								if((curr_mb_label<UNKNOWN_LABEL)&&(label_size[curr_mb_label]>0))
									label_size[curr_mb_label] -= 1;
								mblabel[motion_picture_num][k] = best_neighbor_label;
								label_size[best_neighbor_label] += 1;		
								
								flip_flag = 1;
							}
						}
					}
				}
			}

			recurrent_cnt++;

		}while(flip_flag!=0);

		total_block_groups = rejected_group_num + max_label_idx +1;

		//check the residual existence of each group
		for(i=0; i<=max_label_idx; i++)
			ResGroupON[i] = 0;

		for(mby=0; mby<mbheight; mby++)
			for(mbx=0; mbx<mbwidth; mbx++)
			{
				k= mby*mbwidth + mbx;
				if(mblabel[motion_picture_num][k]>=0)	//no SKIP block
				{
					//Check the Residual Existence
					start_px = mbx*16;
					start_py = mby*16;
					for(py=start_py; py<start_py+16; py++)
						for(px=start_px; px<start_px+16; px++)
						{
							if(curr_res[py][px][0]!=0)
								ResGroupON[mblabel[motion_picture_num][k]] = 1;
						}
				}
			}

		for(i=0; i<=max_label_idx; i++)
		{
			if(ResGroupON[i]==0)
			{
				label_size[i] = 0;
				rejected_group_num++;
			}
		}

		//reindexing each group
		h=0;
		for(i=0; i<=max_label_idx; i++)
		{
			if(label_size[i]>1)
			{
				//inv_display_group_idx[h] = i;
				display_group_idx[i] = h;
				h++;
			}
			else
			{
				display_group_idx[i] = -1;
			}
		}
		max_group_num = h;

		//-----------
		//Spatial and Temporal Filtering: detection of object groups and removal of error groups (객체 그룹의 인식과 에러 그룹의 제거)
		//-----------

		if(object_recognition_flag==1)
		{
			//initialization
			for(i=0; i<MAX_BUF_NUM; i++)
			{
				gbuf[i].ChkOn = 0;
				gbuf[i].DominantGroupIdx = -1;
				tmp_gbuf_Status[i] = gbuf[i].Status;		//buffer status in the previous frame (이전 프레임에서 각 버퍼의 상태)
				if((gbuf[i].CheckCnt>=0)&&((gbuf[i].Status==CANDIDATE_GROUP)||(gbuf[i].Status==OBJECT_GROUP)))
					tmp_gbuf_FrmGroupIdx[i] = gbuf[i].FrmGroupIdx[gbuf[i].CheckCnt];		//group index in the previous frame (이전 프레임에서 각 프레임 그룹 인덱스)
				else
					tmp_gbuf_FrmGroupIdx[i] = -1;
			}
			prev_max_buf_group_num = max_buf_group_num;
			
			for(i=0; i<UNKNOWN_LABEL; i++)
			{
				assigned_buf_idx[i] = -1;
				overlap_flag[i] = 0;
			}

			for(i=0; i<UNKNOWN_LABEL; i++)        //for each group of the current frame: i => block group label (현재 프레임의 각 그룹에 대하여)
			{
				if(label_size[i]>0)
				{
					overlap_check_flag = 0;		// check whether there exists the buffer address corresponding to the current frame (해당 그룹에 일치하는 버퍼 주소가 있는지 여부를 체크한다.)
					curr_buf_address = 0;
					pros_buf_address = -1;
					while((overlap_check_flag!=1)&&(curr_buf_address<prev_max_buf_group_num))        // if there is no checked group and all buffer addresses are not being searched, (체크된 그룹이 없고 버퍼의 모든 주소가 아직 검색되지 않는 동안)
					{
						   if(((tmp_gbuf_Status[curr_buf_address]==CANDIDATE_GROUP)||(tmp_gbuf_Status[curr_buf_address]==OBJECT_GROUP))&&(gbuf[curr_buf_address].CheckCnt>=0))        //if block group is already Candidate or Object group, (버퍼의 주소가 추천 그룹 또는 객체 그룹이라면)
						   {
								  //check if there exists overlap with any group in the current frame (현재 프레임의 그룹과 겹치는 부분이 있는지 체크한다.)
								  overlap_mb_cnt = 0;
								  for(m=0; m<max_mb_nr; m++)        //for all macroblocks in the group of the current frame (현재 프레임의 그룹 내의 모든 매크로블록에 대하여)
								  {
										 if((mblabel[motion_picture_num][m]==i)&&(mblabel[motion_picture_num-prev_dist][m]==tmp_gbuf_FrmGroupIdx[curr_buf_address]))      // if the label of the current macroblock of the current frame is identical to the group label of the current buffer address, (현재 프레임에서 해당 매크로블록의 레이블이 현재 버퍼 주소의 그룹 레이블과 일치하다면)
											   overlap_mb_cnt++;
								  }

								  if(overlap_mb_cnt>0)        //if overlapping (겹치는 부분이 있다면)
								  {
										 overlap_check_flag = 1;

										 if((gbuf[curr_buf_address].ChkOn==1)&&(gbuf[curr_buf_address].DominantGroupIdx>=0)&&(gbuf[curr_buf_address].DominantGroupIdx!=i)
											 &&(gbuf[curr_buf_address].CheckCnt>0))        //if the buffer address was checked, (버퍼의 주소가 체크된 주소라면)
										 {
											   for(m=0; m<max_mb_nr; m++)        //for all macroblocks in the group of the current frame, (현재 프레임에서 그룹 내의 모든 매크로블록에 대하여)
											   {
													 if(mblabel[motion_picture_num][m]==i)
														   mblabel[motion_picture_num][m] = gbuf[curr_buf_address].DominantGroupIdx;      //changing the label of all macroblocks included in this group of the current frame into priority label index of the overlapped existing group. (현재 프레임의 이 그룹에 포함된 모든 매크로블록의 레이블을 겹쳐진 기존 그룹의 최우선 레이블 인덱스로 전환한다.)
											   }

											   assigned_buf_idx[i] = curr_buf_address;		//address of buffer group assigned to the current group. (현재 그룹에 할당된 버퍼 그룹의 주소)
											   gbuf[curr_buf_address].MbNum[gbuf[curr_buf_address].CheckCnt] = label_size[i] + gbuf[curr_buf_address].MbNum[gbuf[curr_buf_address].CheckCnt];
											   label_size[gbuf[curr_buf_address].DominantGroupIdx] = gbuf[curr_buf_address].MbNum[gbuf[curr_buf_address].CheckCnt];
											   label_size[i] = 0;											   
										 }
										 else      //if the buffer address was not checked, (버퍼의 주소가 체크된 주소가 아니라면)
										 {
												gbuf[curr_buf_address].DominantGroupIdx = i;      //setting the group of the current frame up as the priority label index corresponding to the current buffer address. (현재 프레임의 그룹을 현재 버퍼 주소에 대응하는 최우선 레이블 인덱스로 선정한다.)
												assigned_buf_idx[i] = curr_buf_address;		//the address of buffer group assigned the current group. (현재 그룹에 할당된 버퍼 그룹의 주소)
												(gbuf[curr_buf_address].CheckCnt)++;      //update check counter. (check Counter를 업데이트한다.)
												if(tmp_gbuf_Status[curr_buf_address]==OBJECT_GROUP)
													(gbuf[curr_buf_address].FrameCnt)++;
												gbuf[curr_buf_address].MbNum[gbuf[curr_buf_address].CheckCnt] = label_size[i];
										 }
	                                    
										 //judging that the group of the current group has the identical label to the buffer group, update the parameters of the buffer group. (현재 프레임의 그룹이 버퍼의 그룹과 동일한 레이블을 지니는 것으로 판단하고, 버퍼의 그룹에 대한 파라메터를 업데이트한다.)
										 // Calculating a partial Occurrence Probability

										 if(gbuf[curr_buf_address].CheckCnt<=decision_window_size)
										 {
											 if(gbuf[curr_buf_address].MbNum[gbuf[curr_buf_address].CheckCnt-1]>0)
													gbuf[curr_buf_address].OccurProb[gbuf[curr_buf_address].CheckCnt] = (double)overlap_mb_cnt / (double)gbuf[curr_buf_address].MbNum[gbuf[curr_buf_address].CheckCnt-1];      //update occurrence probability (occurrence probability를 업데이트한다.)
										 }
										 gbuf[curr_buf_address].FrmGroupIdx[gbuf[curr_buf_address].CheckCnt] = gbuf[curr_buf_address].DominantGroupIdx;      //update frame group pointer. (frame group pointer를 업데이트한다.)
										 gbuf[curr_buf_address].ChkOn = 1;                                   //update group check status. (group check status를 업데이트한다.)
										 gbuf[curr_buf_address].OccurCnt++;

										 if((tmp_gbuf_Status[curr_buf_address] == CANDIDATE_GROUP)&&(gbuf[curr_buf_address].CheckCnt==decision_window_size))	//if all frames between window interval were checked, (window 간격 사이의 모든 프레임이 체크된 경우,)
										 {
											 //calculate total occurrence probability. (total occurrence probability를 계산한다.)
											  total_occur_prob = 0;
											  for(n=1; n<=decision_window_size; n++)
											  {
												  if((gbuf[curr_buf_address].OccurProb[n]>0)&&(gbuf[curr_buf_address].OccurProb[n]<=1))
														total_occur_prob += -(double)log(gbuf[curr_buf_address].OccurProb[n]);
											  }											  
											  if((gbuf[curr_buf_address].OccurCnt>0)&&(gbuf[curr_buf_address].OccurCnt<=decision_window_size))
											  {
												  for(n=0; n<decision_window_size-gbuf[curr_buf_address].OccurCnt; n++)
														total_occur_prob += -(double)log((double)gbuf[curr_buf_address].OccurCnt/(double)decision_window_size);
											  }
											  else
													total_occur_prob = occur_th +10;

											  if(total_occur_prob<=occur_th)
											  {
													 gbuf[curr_buf_address].Status = OBJECT_GROUP;
													 gbuf[curr_buf_address].FrameCnt = 0;
													 gbuf[curr_buf_address].ObjectFrame = motion_picture_num;
											  }
											  else
											  {
													 gbuf[curr_buf_address].DetectFrame = -1;
													 gbuf[curr_buf_address].Status = ERROR_GROUP;
													 gbuf[curr_buf_address].CheckCnt = -1;
													 gbuf[curr_buf_address].FrameCnt = -1;
													 gbuf[curr_buf_address].ChkOn = 0;
													 gbuf[curr_buf_address].DominantGroupIdx = -1;
													 gbuf[curr_buf_address].OccurCnt = 0;
													 
													 //initialize the label of all macroblocks included in the error group over all window frames. (윈도우 전체 프레임에 대하여, 이 에러 그룹에 포함된 모든 매크로블록의 레이블을 초기화한다.)
													 //picture_num = motion_picture_num;
													 //n = decision_window_size;
													 //while(n>=0)
													 //{
													 //	 frame_type = picture_num%intra_period;
													 //	 if(frame_type==0)	picture_num--;

													 //	 for(m=0; m<max_mb_nr; m++)
													 //	 {
													 //		 if(mblabel[picture_num][m]==gbuf[curr_buf_address].FrmGroupIdx[n])
													 //			 mblabel[picture_num][m] = -1;
													 //	 }
													 //	 picture_num--;
													 //	 n--;
													 //}
											  }
										 }

										 //if((gbuf[curr_buf_address].Status==OBJECT_GROUP)&&(gbuf[curr_buf_address].FrameCnt<intra_period-(gbuf[curr_buf_address].ObjectFrame%intra_period)))
										 if(gbuf[curr_buf_address].Status==OBJECT_GROUP)
										 {
											 mb_cnt = 0;
											 max_mbx = -1;
											 min_mbx = mbwidth+1;
											 max_mby = -1;
											 min_mby = mbheight+1;
											 TmpX = 0;
											 TmpY = 0;
											 for(m=0; m<max_mb_nr; m++)        //for all macroblocks in the current frame groups. (현재 프레임의 그룹 내의 모든 매크로블록에 대하여)
											 {
													 if(mblabel[motion_picture_num][m]==gbuf[curr_buf_address].DominantGroupIdx)
													 {
															 mbx = m%mbwidth;
															 mby = (int)(m/mbwidth);
															 TmpX += mbx;
															 TmpY += mby;

															 if(mbx>max_mbx)	max_mbx = mbx;
															 if(mbx<min_mbx)	min_mbx = mbx;
															 if(mby>max_mby)   max_mby = mby;
															 if(mby<min_mby)   min_mby = mby;

															 mb_cnt++;   
													 }                                            
											 }

											 //Calculate the position and size of object groups
											 if(mb_cnt>0)
											 {
													 gbuf[curr_buf_address].Width[gbuf[curr_buf_address].FrameCnt] = (max_mbx - min_mbx +1)*16;           //update group size의 width
													 gbuf[curr_buf_address].Height[gbuf[curr_buf_address].FrameCnt] = (max_mby - min_mby +1)*16;          //update group size의 height

													 if(gbuf[curr_buf_address].FrameCnt==0)
													 {
														 gbuf[curr_buf_address].PosX[gbuf[curr_buf_address].FrameCnt] = TmpX*16 / mb_cnt +8;         //update group position x coordinate.
														 gbuf[curr_buf_address].PosY[gbuf[curr_buf_address].FrameCnt] = TmpY*16 / mb_cnt +8;         //update group position x coordinate.
													 }
													 else if(gbuf[curr_buf_address].FrameCnt>0)
													 {
														 diff_width = abs(gbuf[curr_buf_address].Width[gbuf[curr_buf_address].FrameCnt] 
																		- gbuf[curr_buf_address].Width[gbuf[curr_buf_address].FrameCnt-1]);
														 diff_height = abs(gbuf[curr_buf_address].Height[gbuf[curr_buf_address].FrameCnt] 
																		- gbuf[curr_buf_address].Height[gbuf[curr_buf_address].FrameCnt-1]);

														 if(diff_width<th_diff_width)
														 {														 
															 gbuf[curr_buf_address].PosX[gbuf[curr_buf_address].FrameCnt] = TmpX*16 / mb_cnt +8;         //update group position x coordinate.
														 }
														 else
														 {
															 gbuf[curr_buf_address].Width[gbuf[curr_buf_address].FrameCnt] = gbuf[curr_buf_address].Width[gbuf[curr_buf_address].FrameCnt-1];
															 gbuf[curr_buf_address].PosX[gbuf[curr_buf_address].FrameCnt] = gbuf[curr_buf_address].PosX[gbuf[curr_buf_address].FrameCnt-1];
														 }

														 if(diff_height<th_diff_height)
														 {
															gbuf[curr_buf_address].PosY[gbuf[curr_buf_address].FrameCnt] = TmpY*16 / mb_cnt +8;         //group position x coordinate.
														 }
														 else
														 {
															 gbuf[curr_buf_address].Height[gbuf[curr_buf_address].FrameCnt] = gbuf[curr_buf_address].Height[gbuf[curr_buf_address].FrameCnt-1];
															 gbuf[curr_buf_address].PosY[gbuf[curr_buf_address].FrameCnt] = gbuf[curr_buf_address].PosY[gbuf[curr_buf_address].FrameCnt-1];
														 }
													 }
											 }											 
										 }
								  }
						   }
						   else        //if the buffer address is in the unused state, (버퍼의 주소가 미사용 상태라면)
						   {
								  if(pros_buf_address<0)        //if the priority buffer address was not set up, (최우선 버퍼 주소가 설정되어 있지 않다면)
										 pros_buf_address = curr_buf_address;        //when a new group is registered, set up the priority buffer address (pros_buf_address). ( 새로운 그룹이 등록될 때 위치할 최우선 버퍼 주소(pros_buf_address)를 설정한다.)
						   }

						   curr_buf_address++;
					}

					if(overlap_check_flag==0)      //if there is no checked address, (in other words, if not overlapped) (체크된 주소가 없다면, 즉 오버랩된 주소가 없다면)
					{
						  //register teh current group as a recommended group of the new address of buffer, and set up related parameters. (프레임의 현재 그룹을 버퍼의 새로운 주소에 추천 그룹으로 등록하고 관련 파라메터를 설정한다.)
						  if(pros_buf_address<0)
							  new_buf_address = max_buf_group_num;
						  else
							  new_buf_address = pros_buf_address;

						  assigned_buf_idx[i] = new_buf_address;
						  gbuf[new_buf_address].CheckCnt = 0;
						  gbuf[new_buf_address].FrameCnt = -1;
						  gbuf[new_buf_address].MbNum[0] = label_size[i];
						  gbuf[new_buf_address].DetectFrame = motion_picture_num;
						  gbuf[new_buf_address].FrmGroupIdx[0] = i;
						  gbuf[new_buf_address].DominantGroupIdx = i;      //set up the group of the current frame as the priority label index corresponding to the current buffer address. (현재 프레임의 그룹을 현재 버퍼 주소에 대응하는 최우선 레이블 인덱스로 선정한다.)
						  
						  /*
						  mb_cnt = 0;
						  max_mbx = -1;
						  min_mbx = mbwidth+1;
						  max_mby = -1;
						  min_mby = mbheight+1;
						  TmpX = 0;
						  TmpY = 0;
						  for(m=0; m<max_mb_nr; m++)        //for all macroblocks in the current frame group. (현재 프레임의 그룹 내의 모든 매크로블록에 대하여)
						  {
								  if(mblabel[motion_picture_num][m]==i)
								  {
										  mbx = m%mbwidth;
										  mby = (int)(m/mbwidth);
										  TmpX += mbx;        
										  TmpY += mby;    

										  if(mbx>max_mbx)	max_mbx = mbx;
										  if(mbx<min_mbx)	min_mbx = mbx;
										  if(mby>max_mby)   max_mby = mby;
										  if(mby<min_mby)   min_mby = mby;

										  mb_cnt++;   
								  }                                            
						  }
						  if(mb_cnt>0)
						  {
								  gbuf[new_buf_address].PosX[0] = TmpX*16 / mb_cnt +8;         //update group position x coordinate.
								  gbuf[new_buf_address].PosY[0] = TmpY*16 / mb_cnt +8;         //update group position x coordinate.							  
						  }
						  gbuf[new_buf_address].Width[0] = (max_mbx - min_mbx +1)*16;           //update group size의 width.
						  gbuf[new_buf_address].Height[0] = (max_mby - min_mby +1)*16;          //update group size의 height.
						  */

						  gbuf[new_buf_address].ChkOn = 1;                                   //update group check status.
						  gbuf[new_buf_address].Status = CANDIDATE_GROUP;

						  if(pros_buf_address<0)
								max_buf_group_num++;
					}
				}
			}

			//체크되지 않은 버퍼 주소에 대하여 파라메터를 업데이트한다.
			//for(m=0; m<max_mb_nr; m++)
			//	buf_mblabel[m] = mblabel[motion_picture_num][m];	//macroblock label in the previous P frame (이전 P 프레임에서의 매크로블록 레이블)

			//Check Block Group Buffer (BGB)

			for(i=0; i<max_buf_group_num; i++)        //for each group of the current frame : i => block group buffer label (현재 프레임의 각 그룹에 대하여)
			{
				  if(gbuf[i].ChkOn!=1)
				  {
					  if((tmp_gbuf_Status[i] == CANDIDATE_GROUP)||(tmp_gbuf_Status[i] == OBJECT_GROUP))
					  {
						  //check whether there is overlap with the current frame group (현재 프레임의 그룹과 겹치는 부분이 있는지 체크한다.)
						  overlap_buf_group_idx = -1;
						  for(m=0; m<max_mb_nr; m++)        //for all macroblock in the current frame group (현재 프레임의 그룹 내의 모든 매크로블록에 대하여)
						  {
								 mb_label = mblabel[motion_picture_num][m];
								 if((mb_label>=0)&&(mb_label<UNKNOWN_LABEL)&&(label_size[mb_label]>0))
								 {
									 if((mblabel[motion_picture_num-prev_dist][m]==tmp_gbuf_FrmGroupIdx[i])&&(assigned_buf_idx[mb_label]!=i)&&(assigned_buf_idx[mb_label]>=0))      //if the buffer group is overlapped with the current group, and it is different with the buffer group assigned previously, (버퍼 그룹이 현재 그룹과 overlap되며, 기존의 할당된 버퍼 그룹과 다르다면,)
									 {
										   overlap_buf_group_idx = assigned_buf_idx[mb_label];
										   if(tmp_gbuf_Status[i] == OBJECT_GROUP)
												overlap_flag[mb_label] = 1;
										   break;
									 }
								 }
						  }
						  
						  if((overlap_buf_group_idx>=0)&&(!((gbuf[overlap_buf_group_idx].Status==OBJECT_GROUP)&&(tmp_gbuf_Status[i]==OBJECT_GROUP)))
							  &&(gbuf[overlap_buf_group_idx].Status!=ERROR_GROUP)&&(tmp_gbuf_Status[i]!=ERROR_GROUP))        //if there is overlap, but the previous group and overlapped group are not simultaneously object groups, (겹치는 부분이 있되, 이전 그룹과 overlap 그룹이 한꺼번에 객체 그룹이 되지 않는 경우,)
						  {
							  if(tmp_gbuf_Status[i]==CANDIDATE_GROUP)		//if previous group is recommended group, (이전 그룹이 추천 그룹이라면)
							  {
								  //recognize the previous group as error group (이전 그룹을 에러 그룹으로 처리한다.)
								  gbuf[i].DetectFrame = -1;
								  gbuf[i].Status = ERROR_GROUP;
								  gbuf[i].CheckCnt = -1;
								  gbuf[i].FrameCnt = -1;
								  gbuf[i].ChkOn = 0;
								  gbuf[i].DominantGroupIdx = -1;
								  gbuf[i].OccurCnt = 0;
							  }
							  else if(tmp_gbuf_Status[i]==OBJECT_GROUP)
							  {
								  if((gbuf[overlap_buf_group_idx].Status==CANDIDATE_GROUP)&&(gbuf[i].CheckCnt>=0)&&(gbuf[overlap_buf_group_idx].CheckCnt>=0))		//if the previous group is object group, and overlap group is recommended group, (이전 그룹이 객체 그룹이고, overlap 그룹이 추천 그룹이라면)
								  {
									  //copy all parameters of the existing overlap group into the previous group (기존 overlap 그룹이 갖고 있는 모든 파라메터를 그대로 이전 그룹으로 복사한다.)
									  (gbuf[i].CheckCnt)++;
									  (gbuf[i].FrameCnt)++;
									  gbuf[i].MbNum[gbuf[i].CheckCnt] = gbuf[overlap_buf_group_idx].MbNum[gbuf[overlap_buf_group_idx].CheckCnt];
									  if(gbuf[i].CheckCnt<=decision_window_size)
											gbuf[i].OccurProb[gbuf[i].CheckCnt] = gbuf[overlap_buf_group_idx].OccurProb[gbuf[overlap_buf_group_idx].CheckCnt];
									  gbuf[i].FrmGroupIdx[gbuf[i].CheckCnt] = gbuf[overlap_buf_group_idx].FrmGroupIdx[gbuf[overlap_buf_group_idx].CheckCnt];
									  //if(gbuf[i].FrameCnt<intra_period-(gbuf[i].ObjectFrame%intra_period))
									  //{
											 mb_cnt = 0;
											 max_mbx = -1;
											 min_mbx = mbwidth+1;
											 max_mby = -1;
											 min_mby = mbheight+1;
											 TmpX = 0;
											 TmpY = 0;
											 for(m=0; m<max_mb_nr; m++)        //for all macroblocks in the current frame group (현재 프레임의 그룹 내의 모든 매크로블록에 대하여)
											 {
													 if(mblabel[motion_picture_num][m]==gbuf[overlap_buf_group_idx].FrmGroupIdx[gbuf[overlap_buf_group_idx].CheckCnt])
													 {
															 mbx = m%mbwidth;
															 mby = (int)(m/mbwidth);
															 TmpX += mbx;
															 TmpY += mby;

															 if(mbx>max_mbx)	max_mbx = mbx;
															 if(mbx<min_mbx)	min_mbx = mbx;
															 if(mby>max_mby)   max_mby = mby;
															 if(mby<min_mby)   min_mby = mby;

															 mb_cnt++;   
													 }                                            
											 }

											 if(mb_cnt>0)
											 {
													 gbuf[i].Width[gbuf[i].FrameCnt] = (max_mbx - min_mbx +1)*16;           //update group size의 width.
													 gbuf[i].Height[gbuf[i].FrameCnt] = (max_mby - min_mby +1)*16;          //update group size의 height.

													 if(gbuf[i].FrameCnt==0)
													 {
														 gbuf[i].PosX[gbuf[i].FrameCnt] = TmpX*16 / mb_cnt +8;         //update group position x coordinate.
														 gbuf[i].PosY[gbuf[i].FrameCnt] = TmpY*16 / mb_cnt +8;         //update group position x coordinate.
													 }
													 else if(gbuf[i].FrameCnt>0)
													 {
														 diff_width = abs(gbuf[i].Width[gbuf[i].FrameCnt] - gbuf[i].Width[gbuf[i].FrameCnt-1]);
														 diff_height = abs(gbuf[i].Height[gbuf[i].FrameCnt] - gbuf[i].Height[gbuf[i].FrameCnt-1]);

														 if(diff_width<th_diff_width)
														 {														 
															 gbuf[i].PosX[gbuf[i].FrameCnt] = TmpX*16 / mb_cnt +8;         //update group position x coordinate.
														 }
														 else
														 {
															 gbuf[i].Width[gbuf[i].FrameCnt] = gbuf[i].Width[gbuf[i].FrameCnt-1];
															 gbuf[i].PosX[gbuf[i].FrameCnt] = gbuf[i].PosX[gbuf[i].FrameCnt-1];
														 }

														 if(diff_height<th_diff_height)
														 {
															 gbuf[i].PosY[gbuf[i].FrameCnt] = TmpY*16 / mb_cnt +8;         //update group position x coordinate.
														 }
														 else
														 {
															 gbuf[i].Height[gbuf[i].FrameCnt] = gbuf[i].Height[gbuf[i].FrameCnt-1];
															 gbuf[i].PosY[gbuf[i].FrameCnt] = gbuf[i].PosY[gbuf[i].FrameCnt-1];
														 }
													 }
											 }	
									  //}
									  gbuf[i].ChkOn = 1;
									  gbuf[i].OccurCnt++;
									  gbuf[i].DominantGroupIdx = gbuf[overlap_buf_group_idx].DominantGroupIdx;
									  
									  //recognize this overlap group as error group (이 overlap 그룹은 에러 그룹으로 처리한다.)
									  gbuf[overlap_buf_group_idx].DetectFrame = -1;
									  gbuf[overlap_buf_group_idx].Status = ERROR_GROUP;
									  gbuf[overlap_buf_group_idx].CheckCnt = -1;
									  gbuf[overlap_buf_group_idx].FrameCnt = -1;
									  gbuf[overlap_buf_group_idx].ChkOn = 0;
									  gbuf[overlap_buf_group_idx].DominantGroupIdx = -1;
									  gbuf[overlap_buf_group_idx].OccurCnt = 0;
								  }
							  }
						  }
						  else		// if there is no overlap and both the previous group and overlap group are object groups, (겹치는 부분이 없거나, 이전 그룹과 overlap 그룹이 모두 객체 그룹인 경우,)
						  {
							  if(gbuf[i].CheckCnt>=0)
							  {
								  (gbuf[i].CheckCnt)++;      //update check Counter.
								  if(tmp_gbuf_Status[i] == OBJECT_GROUP)
									  (gbuf[i].FrameCnt)++;

								  //include the virtual MB group to the temporary macoblock label buffer. ( 임시 매크로블록 레이블 버퍼에 가상의 매크로블록 그룹을 포함시킨다.)
								  max_label_idx++;
								  cnt=0;
								  for(m=0; m<max_mb_nr; m++)        //for all MBs in the current frame group. (현재 프레임의 그룹 내의 모든 매크로블록에 대하여)
								  {
									  if(mblabel[motion_picture_num-prev_dist][m]==tmp_gbuf_FrmGroupIdx[i])
									  {
										  mblabel[motion_picture_num][m] = max_label_idx;
										  cnt++;
									  }
								  }
								  label_size[max_label_idx] = cnt;

								  gbuf[i].FrmGroupIdx[gbuf[i].CheckCnt] = max_label_idx;      //updte frame group pointer.
								  gbuf[i].MbNum[gbuf[i].CheckCnt] = gbuf[i].MbNum[gbuf[i].CheckCnt-1];
								  if(gbuf[i].CheckCnt<=decision_window_size)
										gbuf[i].OccurProb[gbuf[i].CheckCnt] = 1;      //updte occurrence probability.								  
								  gbuf[i].ChkOn = 1;                                   //updte group check status.

								  if((tmp_gbuf_Status[i] == CANDIDATE_GROUP)&&(gbuf[i].CheckCnt==decision_window_size))	//if all frames in the window interval were checked. (window 간격 사이의 모든 프레임이 체크된 경우,)
								  {
									   //calculate total occurrence probability. (total occurrence probability를 계산한다.)
									   total_occur_prob = 0;
									   for(n=1; n<=decision_window_size; n++)
									   {
										   if((gbuf[i].OccurProb[n]>0)&&(gbuf[i].OccurProb[n]<=1))
												total_occur_prob += -(double)log(gbuf[i].OccurProb[n]);
									   }
									   if((gbuf[i].OccurCnt>0)&&(gbuf[i].OccurCnt<=decision_window_size))
											total_occur_prob += -(double)log((double)gbuf[i].OccurCnt/(double)decision_window_size);
									   else
										   total_occur_prob = occur_th +10;

									   if(total_occur_prob<=occur_th)
									   {
											 gbuf[i].Status = OBJECT_GROUP;
											 gbuf[i].FrameCnt = 0;
											 gbuf[i].ObjectFrame = motion_picture_num;
									   }
									   else
									   {
											 gbuf[i].DetectFrame = -1;
											 gbuf[i].Status = ERROR_GROUP;
											 gbuf[i].CheckCnt = -1;
											 gbuf[i].FrameCnt = -1;
											 gbuf[i].ChkOn = 0;
											 gbuf[i].DominantGroupIdx = -1;
											 gbuf[i].OccurCnt = 0;

											 //initialize all MB labels in the error group over all window frames (윈도우 전체 프레임에 대하여, 이 에러 그룹에 포함된 모든 매크로블록의 레이블을 초기화한다.)
											 //picture_num = motion_picture_num;
											 //n = decision_window_size;
											 //while(n>=0)
											 //{
											 //	 frame_type = picture_num%intra_period;
											 //	 if(frame_type==0)	picture_num--;

											 //	 for(m=0; m<max_mb_nr; m++)
											 //	 {
											 //		 if(mblabel[picture_num][m]==gbuf[i].FrmGroupIdx[n])
											 //			 mblabel[picture_num][m] = -1;
											 //	 }
											 //	 picture_num--;
											 //	 n--;
											 //}
									   }
								  }

								  //if((gbuf[i].Status==OBJECT_GROUP)&&(gbuf[i].FrameCnt<intra_period-(gbuf[i].ObjectFrame%intra_period)))
								  if(tmp_gbuf_Status[i]==OBJECT_GROUP)
								  {
									  gbuf[i].PosX[gbuf[i].FrameCnt] = gbuf[i].PosX[gbuf[i].FrameCnt-1];         //update group position x coordinate.
									  gbuf[i].PosY[gbuf[i].FrameCnt] = gbuf[i].PosY[gbuf[i].FrameCnt-1];         //update group position x coordinate.
									  gbuf[i].Width[gbuf[i].FrameCnt] = gbuf[i].Width[gbuf[i].FrameCnt-1];           //update group size의 width.
									  gbuf[i].Height[gbuf[i].FrameCnt] = gbuf[i].Height[gbuf[i].FrameCnt-1];          //update group size의 height.
								  }
								  /*
								  else
								  {
									  mb_cnt = 0;
									  max_mbx = -1;
									  min_mbx = mbwidth+1;
									  max_mby = -1;
									  min_mby = mbheight+1;
									  TmpX = 0;
									  TmpY = 0;
									  for(m=0; m<max_mb_nr; m++)        //for all Mbs in the current frame group. (현재 프레임의 그룹 내의 모든 매크로블록에 대하여)
									  {
											 if(mblabel[motion_picture_num][m]==max_label_idx)
											 {
													 mbx = m%mbwidth;
													 mby = (int)(m/mbwidth);
													 TmpX += mbx;
													 TmpY += mby;

													 if(mbx>max_mbx)	max_mbx = mbx;
													 if(mbx<min_mbx)	min_mbx = mbx;
													 if(mby>max_mby)   max_mby = mby;
													 if(mby<min_mby)   min_mby = mby;

													 mb_cnt++;   
											 }                                            
									  }

									  if(mb_cnt>0)
									  {
											 gbuf[i].Width[gbuf[i].FrameCnt] = (max_mbx - min_mbx +1)*16;           //update group size의 width.
											 gbuf[i].Height[gbuf[i].FrameCnt] = (max_mby - min_mby +1)*16;          //update group size의 height.

											 diff_width = abs(gbuf[i].Width[gbuf[i].FrameCnt] - gbuf[i].Width[gbuf[i].FrameCnt-1]);
											 diff_height = abs(gbuf[i].Height[gbuf[i].FrameCnt] - gbuf[i].Height[gbuf[i].FrameCnt-1]);

											 if(diff_width<th_delta_width)
											 {														 
												 gbuf[i].PosX[gbuf[i].FrameCnt] = TmpX*16 / mb_cnt +8;         //update group position x coordinate.
											 }
											 else
											 {
												 gbuf[i].Width[gbuf[i].FrameCnt] = gbuf[i].Width[gbuf[i].FrameCnt-1];
												 gbuf[i].PosX[gbuf[i].FrameCnt] = gbuf[i].PosX[gbuf[i].FrameCnt-1];
											 }

											 if(diff_height<th_delta_height)
											 {
												gbuf[i].PosY[gbuf[i].FrameCnt] = TmpY*16 / mb_cnt +8;         //update group position x coordinate.
											 }
											 else
											 {
												 gbuf[i].Height[gbuf[i].FrameCnt] = gbuf[i].Height[gbuf[i].FrameCnt-1];
												 gbuf[i].PosY[gbuf[i].FrameCnt] = gbuf[i].PosY[gbuf[i].FrameCnt-1];
											 }
									  }
								  }*/
							  }
						  }
					  }
				  }
			}

			//overlapped active groups의 처리
			for(i=0; i<UNKNOWN_LABEL; i++)
			{
				if(overlap_flag[i]==1)
				{
					for(m=0; m<max_mb_nr; m++)        //for all MBs in the current frame group (현재 프레임의 그룹 내의 모든 매크로블록에 대하여)
					{
						  if((mblabel[motion_picture_num][m]==i)&&(motion_picture_num>=prev_dist)
							  &&((mblabel[motion_picture_num-prev_dist][m]!=tmp_gbuf_FrmGroupIdx[assigned_buf_idx[i]])))
							  mblabel[motion_picture_num][m] = -1;
					}
				}
			}

			for(i=0; i<max_buf_group_num; i++)
				gbuf_Status[i][curr_frame_type-1] = gbuf[i].Status;

			//for(i=0; i<max_buf_group_num; i++)        //for each group of the current frame (현재 프레임의 각 그룹에 대하여)
			//{
			//	if(gbuf[i].Status==OBJECT_GROUP)
			//	{
			//		if(motion_picture_num!=gbuf[i].ObjectFrame+gbuf[i].FrameCnt)
			//			n=0;
			//	}
			//}
		}

		//buffer for P-frame (P-frame을 위한 버퍼)
		if(time_measure_flag==0)
		{
			for(i=0; i<iheight; i++)
			{
				for(j=0; j<iwidth; j++)
				{
					frm_buf[curr_frame_type-1][0][i*iwidth+j] = curr_frm[i][j][0];
					frm_buf[curr_frame_type-1][1][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][1];
					frm_buf[curr_frame_type-1][2][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][2];
				}
			}
		}
		
	}

	//drawing each group in terms of color
	if((read_mb_type_flag==1)&&(time_measure_flag==0))
		draw_mb_type();

	//drawing the motion trajectory
	if(trajectory_optimization_flag!=1)
	{
		if((extract_trajectory_flag==1)&&(object_recognition_flag==1)&&(time_measure_flag==0))
			draw_trajectory();
	}

	//output the status of group buffers
	if((output_group_buffer_flag==1)&&(time_measure_flag==0))
	{
		if(max_buf_group_num==0)
		{
			fprintf(BUF_FILE,"%d\n",motion_picture_num);
		}
		else if(max_buf_group_num>0)
		{
			fprintf(BUF_FILE,"%d\t",motion_picture_num);
			for(i=0; i<max_buf_group_num; i++)
			{
				if(i<max_buf_group_num-1)
					fprintf(BUF_FILE,"%d\t",gbuf[i].Status);
				else if(i==max_buf_group_num-1)
					fprintf(BUF_FILE,"%d\n",gbuf[i].Status);
			}
		}
	}

	//output the spatial noise rate
	if((output_spatial_noise_flag==1)&&(time_measure_flag==0))
	{
		if(curr_frame_type==0)		//I-frames
		{
			fprintf(SPA_FILE,"%d\n",motion_picture_num);
		}
		else
		{
			if(total_block_groups>0)
				spatial_noise_rate = (double)rejected_group_num/(double)total_block_groups;
			else
				spatial_noise_rate = 0;

			fprintf(SPA_FILE,"%d\t%d\t%d\t%.5f\n",motion_picture_num,rejected_group_num,total_block_groups,spatial_noise_rate);
		}
	}

	free(bkd_obj);
}

void draw_trajectory()
{
	int i,j,k,f;
	int R,G,B;
	int group_x, group_y, group_width, group_height;
	int curr_idx;
	int intra_period = INTRA_PERIOD;	//IPPP
	int curr_frame_type = motion_picture_num%intra_period;

	R=30; G=130; B=130;

	// Get a decoded frame
	if(curr_frame_type==0)
	{
		//해당 I-frame의 저장
		for(i=0; i<iheight; i++)
		{
			for(j=0; j<iwidth; j++)
			{
				frm_buf[intra_period-1][0][i*iwidth+j] = curr_frm[i][j][0];
				frm_buf[intra_period-1][1][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][1];
				frm_buf[intra_period-1][2][(i>>1)*(iwidth>>1)+(j>>1)] = curr_frm[i>>1][j>>1][2];
			}
		}
		
		if(motion_picture_num>0)
		{
			for(f=0; f<intra_period; f++)		//GOP의 각 프레임에 대하여
			{
				//그룹 버퍼로부터 객체 인덱스, 위치, 크기 정보를 추출한다.
				for(k=0; k<max_buf_group_num; k++)        //현재 프레임의 각 그룹에 대하여
				{
					curr_idx = gbuf[k].FrameCnt -(intra_period-1) +f;
					if(((gbuf_Status[k][f]==OBJECT_GROUP))&&(curr_idx>=0))
					{
						if(curr_idx>=0)
						{
							group_x = gbuf[k].PosX[curr_idx];
							group_y = gbuf[k].PosY[curr_idx];
							group_width = gbuf[k].Width[curr_idx];
							group_height = gbuf[k].Height[curr_idx];
						
							draw_rectangle(1, 1, R, G, B, frm_buf[f][0], frm_buf[f][1], frm_buf[f][2], group_x-group_width/2, group_y-group_height/2, group_x+group_width/2, group_y+group_height/2);		//black
						}
					}
				}

				//Write a file
				for(i = 0; i < iheight ; i++)
					fwrite(&(frm_buf[f][0][i*iwidth]),1,iwidth,MTR_FILE );

				for(i = 0; i < (iheight/2); i++)
					fwrite(&(frm_buf[f][1][i*iwidth/2]),1,iwidth/2,MTR_FILE );

				for(i = 0; i < (iheight/2); i++)
					fwrite(&(frm_buf[f][2][i*iwidth/2]),1,iwidth/2,MTR_FILE );
			}
		}
		else if(motion_picture_num==0)
		{
			//Write a file
			for(i = 0; i < iheight ; i++)
				fwrite(&(frm_buf[intra_period-1][0][i*iwidth]),1,iwidth,MTR_FILE );

			for(i = 0; i < (iheight/2); i++)
				fwrite(&(frm_buf[intra_period-1][1][i*iwidth/2]),1,iwidth/2,MTR_FILE );

			for(i = 0; i < (iheight/2); i++)
				fwrite(&(frm_buf[intra_period-1][2][i*iwidth/2]),1,iwidth/2,MTR_FILE );
		}
	}
}

void pred_segment()
{
	//This function predicts the location and region of objects in the I-frame based on some previous frames. 
	//Also, partial decoding is performed for these regions. The location of object is computed via linear extrapolation, 
	//and the size of region is decided as the maximum width and height among the previous three frames plus margin.
	//일정한 개수의 이전 프레임으로부터 I 프레임에서 객체의 위치와 영역을 예측한다.
	//이 영역에 대하여 부분적인 디코딩을 수행한다.
	//객체의 위치는 linear extrapolation에 의하여 구하고, 영역 크기는 이전 3개 프레임의 객체 영역에 대한 최대 width 및 height의 값에 margin을 더한 것으로 한다.

	int k,i;
	double slope, pram;
	int max_width, max_height;
	int intra_period = INTRA_PERIOD;	//IPPP
	int curr_frame_type = motion_picture_num%intra_period;
	int window_size = 3;
	int prev_x[2], prev_y[2];
	int mbx, mby;
	int min_x, max_x, min_y, max_y, min_mbx, max_mbx, min_mby, max_mby;
	int mb_height = (int)(iheight/16);
	int mb_width = (int)(iwidth/16);

	if(curr_frame_type==0)		//I-frame에 대해서
	{
		if(motion_picture_num>0)
		{
			for(k=0; k<max_buf_group_num; k++)        //현재 프레임의 각 그룹에 대하여
			{
				if(gbuf[k].Status==OBJECT_GROUP)
				{
					/*
					if(gbuf[k].FrameCnt>=window_size-1)
					{
						//위치를 예측한다
						prev_x[0] = gbuf[k].PosX[gbuf[k].FrameCnt-1];
						prev_y[0] = gbuf[k].PosY[gbuf[k].FrameCnt-1];
						prev_x[1] = gbuf[k].PosX[gbuf[k].FrameCnt];
						prev_y[1] = gbuf[k].PosY[gbuf[k].FrameCnt];

						gbuf[k].PosX[gbuf[k].FrameCnt+1] = 2*prev_x[1] - prev_x[0];
						gbuf[k].PosY[gbuf[k].FrameCnt+1] = 2*prev_y[1] - prev_y[0];

						//영역 크기를 예측한다.
						max_width = -1;
						max_height = -1;
						for(i=0; i<window_size; i++)
						{
							if(max_width<gbuf[k].Width[gbuf[k].CheckCnt-i])
								max_width = gbuf[k].Width[gbuf[k].CheckCnt-i];
							if(max_height<gbuf[k].Height[gbuf[k].CheckCnt-i])
								max_height = gbuf[k].Height[gbuf[k].CheckCnt-i];
						}

						gbuf[k].Width[gbuf[k].FrameCnt+1] = max_width;
						gbuf[k].Height[gbuf[k].FrameCnt+1] = max_height;
					}
					else
					{
					*/

						gbuf[k].PosX[gbuf[k].FrameCnt+1] = gbuf[k].PosX[gbuf[k].FrameCnt];
						gbuf[k].PosY[gbuf[k].FrameCnt+1] = gbuf[k].PosY[gbuf[k].FrameCnt];
						gbuf[k].Width[gbuf[k].FrameCnt+1] = gbuf[k].Width[gbuf[k].FrameCnt];
						gbuf[k].Height[gbuf[k].FrameCnt+1] = gbuf[k].Height[gbuf[k].FrameCnt];

						/*
						min_x = ClipX(gbuf[k].PosX[gbuf[k].FrameCnt+1] - gbuf[k].Width[gbuf[k].FrameCnt+1]/2 - HOR_DECODING_MARGIN);
						max_x = ClipX(gbuf[k].PosX[gbuf[k].FrameCnt+1] + gbuf[k].Width[gbuf[k].FrameCnt+1]/2 + HOR_DECODING_MARGIN);
						min_y = ClipY(gbuf[k].PosY[gbuf[k].FrameCnt+1] - gbuf[k].Height[gbuf[k].FrameCnt+1]/2 - VER_DECODING_MARGIN);
						max_y = ClipY(gbuf[k].PosY[gbuf[k].FrameCnt+1] + gbuf[k].Height[gbuf[k].FrameCnt+1]/2 + VER_DECODING_MARGIN);

						min_mbx = (int)(min_x/16);
						max_mbx = (int)(max_x/16);
						min_mby = (int)(min_y/16);
						max_mby = (int)(max_y/16);

						for(mby=0; mby<mb_height; mby++)			
						{
							for(mbx=0; mbx<mb_width; mbx++)			
							{
								if((mby>=min_mby)&&(mby<=max_mby)&&(mbx>=min_mbx)&&(mbx<=max_mbx))
									blk_dec_enb[mby*mb_width+mbx] = 1;
								else
									blk_dec_enb[mby*mb_width+mbx] = 0;
							}
						}
						*/
					//}
				}
			}
		}
	}

}
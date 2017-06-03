
/*!
 ************************************************************************
 * \file image.h
 *
 * \brief
 *    prototypes for image.c
 *
 ************************************************************************
 */

#ifdef __cplusplus
extern "C" { //visible only to a C++ compiler
#endif

#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "mbuffer.h"

// ws 070829
#define IsNearlyZero(a,tol)			(((a)>=0)&&((a)<=(tol)))||(((a)<0)&&((a)>=-(tol)))
#define ClipX(x)					(x<0?0:(x>=iwidth?iwidth-1:x))
#define ClipY(y)					(y<0?0:(y>=iheight?iheight-1:y))
#define ClipMBX(mbx)				(mbx<0?0:(mbx>=iwidth/16?iwidth/16-1:mbx))
#define ClipMBY(mby)				(mby<0?0:(mby>=iheight/16?iheight/16-1:mby))
//#define ContrastStretch(x)		(x==0?255:0)
//#define ContrastStretch(x)		((255-x)*(255-x)/255)
#define ContrastStretch(a,x)		(x<=(255/(1+a))?(255-a*x):((255-x)/a))
#define IS_ASSIGNED(x)				(x<0?0:(x>=UNKNOWN_LABEL?0:1))

extern StorablePicture *dec_picture;

void find_snr(struct snr_par *snr, StorablePicture *p, int p_ref);
void get_block(int ref_frame, StorablePicture **list, int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE]);
int  picture_order(struct img_par *img);

//Added by Wonsang You at JAN 26 2007
void partial_store_motion_into_buffer();
void partial_update_motion_flow(int bx, int by);
int block_availability(int bx, int by);
void exit_slice();
void write_point(int thickness, int R, int G, int B, unsigned char *Y, unsigned char *U, unsigned char *V, int px, int py);
void write_line(int thickness, int R, int G, int B, unsigned char *Y, unsigned char *U, unsigned char *V, int px1, int py1, int px2, int py2);
void extract_fmv();
void draw_mb_type();
void draw_rectangle(int ShowFlag, int thickness, int R, int G, int B, unsigned char *Y, unsigned char *U, unsigned char *V, int TopLeft_x, int TopLeft_y, int BottomRight_x, int BottomRight_y);
void motion_segment();
void select_blocks(int DrawBlockLevel, int TopLeft_x, int TopLeft_y, int BottomRight_x, int BottomRight_y);
void draw_trajectory();
void pred_segment();
int gen_segment_decblock_list();
//Added-End

#endif

#ifdef __cplusplus
} //visible only to a C++ compiler
#endif
/*!
 ************************************************************************
 *  \file
 *     loopfilter.h
 *  \brief
 *     external loop filter interface
 ************************************************************************
 */

#ifdef __cplusplus
extern "C" { //visible only to a C++ compiler
#endif

#ifndef _LOOPFILTER_H_
#define _LOOPFILTER_H_

#include "global.h"
#include "mbuffer.h"

void DeblockPicture(struct img_par *img, StorablePicture *p) ;

#endif //_LOOPFILTER_H_

#ifdef __cplusplus
} //visible only to a C++ compiler
#endif
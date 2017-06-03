
/*!
 **************************************************************************************
 * \file
 *    parset.h
 * \brief
 *    Picture and Sequence Parameter Sets, encoder operations
 *    This code reflects JVT version xxx
 *  \date 25 November 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details) 
 *      - Stephan Wenger        <stewe@cs.tu-berlin.de>
 ***************************************************************************************
 */

#ifdef __cplusplus
extern "C" { //visible only to a C++ compiler
#endif


#ifndef _NALU_H_
#define _NALU_H_

#include <stdio.h>
#include "nalucommon.h"

extern FILE *bits;
//Added by Wonsang You at FEB 7, 2007
extern FILE *bitss;
//Added-End

int GetAnnexbNALU (NALU_t *nalu);
int NALUtoRBSP (NALU_t *nalu);

#endif

#ifdef __cplusplus
} //visible only to a C++ compiler
#endif

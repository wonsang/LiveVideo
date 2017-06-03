
/*!
 *************************************************************************************
 * \file annexb.h
 *
 * \brief
 *    Annex B byte stream buffer handling.
 *
 *************************************************************************************
 */

#ifdef __cplusplus
extern "C" { //visible only to a C++ compiler
#endif

#ifndef _ANNEXB_H_
#define _ANNEXB_H_

#include "nalucommon.h"

extern int IsFirstByteStreamNALU;
extern int LastAccessUnitExists;
extern int NALUCount;

int  GetAnnexbNALU (NALU_t *nalu);
void OpenBitstreamFile (char *fn);
void CloseBitstreamFile();
void CheckZeroByteNonVCL(NALU_t *nalu, int * ret);
void CheckZeroByteVCL(NALU_t *nalu, int * ret);

//Added by Wonsang You at FEB 8, 2007
int  partial_GetAnnexbNALU (NALU_t *nalu, unsigned char *Buf);
void partial_CloseBitstreamFile();
//Added-End

#endif

#ifdef __cplusplus
} //visible only to a C++ compiler
#endif
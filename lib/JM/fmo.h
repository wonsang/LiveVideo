
/*!
 ***************************************************************************
 *
 * \file fmo.h
 *
 * \brief
 *    Support for Flexilble Macroblock Ordering (FMO)
 *
 * \date
 *    19 June, 2002
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 **************************************************************************/

#ifdef __cplusplus
extern "C" { //visible only to a C++ compiler
#endif

#ifndef _FMO_H_
#define _FMO_H_


int FmoInit (pic_parameter_set_rbsp_t* pps, seq_parameter_set_rbsp_t* sps);
int FmoFinit ();

int FmoGetNumberOfSliceGroup();
int FmoGetLastMBOfPicture();
int FmoGetLastMBInSliceGroup(int SliceGroup);
int FmoGetSliceGroupId (int mb);
int FmoGetNextMBNr (int CurrentMbNr);

#endif

#ifdef __cplusplus
} //visible only to a C++ compiler
#endif
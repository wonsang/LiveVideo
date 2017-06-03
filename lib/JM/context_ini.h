
/*!
 *************************************************************************************
 * \file context_ini.h
 *
 * \brief
 *    CABAC context initializations
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Heiko Schwarz                   <hschwarz@hhi.de>
 **************************************************************************************
 */

#ifdef __cplusplus
extern "C" { //visible only to a C++ compiler
#endif


#ifndef _CONTEXT_INI_
#define _CONTEXT_INI_

void  init_contexts  (struct img_par* img);

#endif

#ifdef __cplusplus
} //visible only to a C++ compiler
#endif


/*!
 ****************************************************************************
 * \file errorconcealment.h
 *
 * \brief
 *    Header file for errorconcealment.c
 *
 ****************************************************************************
 */

#ifdef __cplusplus
extern "C" { //visible only to a C++ compiler
#endif

#ifndef _ERRORCONCEALMENT_H_
#define _ERRORCONCEALMENT_H_

int set_ec_flag(int se);
void reset_ec_flags();
int get_concealed_element(SyntaxElement *sym);

#endif

#ifdef __cplusplus
} //visible only to a C++ compiler
#endif
/* divide_data.h */
/*********************************************/
/*	data -> explained & explanatory				*/
/*********************************************/
#ifndef DIVIDE_DATA_H
#define DIVIDE_DATA_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*	data -> explained & explanatory */

extern
void DataToExpl(	/* input */
						int		n,
						int		p,
						int		i_ex,
						double	*data,		/* [n*(p+1)] */
						/* output */
						double	*explained,		/* [n] */ 
						double	*explanatory	/* [n*p] */ 
						);

#ifdef __cpluscplus
}
#endif



/* template */
/*
#ifndef F_
#define F_*/
/* */
/*
#endif
*/

#endif

/* read_data.h */
/*********************************************/
/*	read data ( dim of matrix )					*/
/*	read data ( data )								*/
/*********************************************/
#ifndef READ_DATA_H
#define READ_DATA_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

extern
void ReadFail(	int status
					);

extern
void ReadDim(	const char	*filename,
					int			*n,			/* number of data */ 
					int			*p,			/* number of explanatory variable */ 
					int			*i_ex			/* index of explained varaible */
					);

extern
void ReadData(	const char	*filename,
					int			n,
					int			p,
					double		*data
					);

#ifdef __cplusplus
}
#endif

#endif

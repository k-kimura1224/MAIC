/* divide_data.h */
/*********************************************/
/*	data -> explained & explanatory				*/
/*********************************************/

#include <stdio.h>

#include "divide_data.h"

/*	data -> explained & explanatory */
void DataToExpl(	/* input */
						int		n,
						int		p,
						int		i_ex,
						double	*data,		/* [n*(p+1)] */
						/* output */
						double	*explained,		/* [n] */
						double	*explanatory	/* [n*p] */
						)
{
	int i,j;
	int ct=0;

	for(i=0; i<n; ++i){
		*(explained+i) = *(data+(i*(p+1))+(i_ex-1));
	}

	for(i=0; i<n; ++i){
		for(j=0; j<(p+1); ++j){
			if( j!=(i_ex-1) ){
				*(explanatory+ct) = *(data+(i*(p+1))+j);
				ct++;
			}
		}
	}
}



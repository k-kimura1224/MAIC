/* matrix.h */
/*********************************************/
/*		generate Identity Matrix					*/
/* 	generate zero Matrix 						*/
/* 	generate zero Matrix (integer)			*/
/*		compute matrix AB								*/
/*		compute vector Ax								*/
/*		compute transposed matrix					*/
/*		change vector a <-> b						*/
/*		compute matrix aA+bB (n,m)					*/
/*		compute matrix aX, a:scalar				*/
/*		compute 	Ax = b 								*/
/*********************************************/

#include <stdio.h>

#include "matrix.h"

void GenerateZeroMatInt(	int		n,
									int		m,
	/* integer */				int		*O
									)
{
	int i,j;
	int ct=0;
	for(i=0; i<n; ++i){
		for(j=0; j<m; ++j){
			*(O+ct)	=	0;
			ct++;
		}
	}
}

/* transpose A */
void TraMat(	int 		n,
					int		m,
					double	*A,	/* [n*m] or [n][m] */
					double	*T		/* (m,n) */
){
	int i,j;

	for(i=0; i<m; ++i){
		for(j=0; j<n; ++j){
			*(T+(n*i)+j) = *(A+(j*m)+i);
		}
	}
}


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

#ifndef __SCIP_MATRIX__
#define __SCIP_MATRIX__

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

extern
void GenerateZeroMatInt(	int		n,
									int		m,
	/* integer */				int		*O
									);

/* transpose A */
extern
void TraMat(	int 		n,
					int		m,
					double	*A,	/* [n*m] or [n][m] */
					double	*T		/* (m,n) */
);

#ifdef __cplusplus
}
#endif

#endif

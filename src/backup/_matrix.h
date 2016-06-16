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
//		compute 	Ax = b 
/*********************************************/
#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
//#if defined(__APPLE__)
//#include "clapack.h"
//#endif
#include "vector.h"

#ifndef F_GENERATEIDENMAT
#define F_GENERATEIDENMAT
/* generate Identity Matrix */
void GenerateIdenMat(	int		n,
								double	I[n][n]
){
	int i,j;
	for(i=0; i<n; ++i){
		for(j=0; j<n; ++j){
			I[i][j] = 0;
		}
		I[i][i] = 1;
	}
}
#endif

#ifndef F_GENERATEZEROMAT
#define F_GENERATEZEROMAT
/* generate zero Matrix */
void GenerateZeroMat(	int		n,
								double	O[n][n]
){
	int	i,j;
	for(i=0; i<n; ++i){
		for(j=0; j<n; ++j){
			O[i][j] = 0;
		}
	}
}
#endif

#ifndef F_GENERATEZEROMATINT
#define F_GENERATEZEROMATINT
void GenerateZeroMatInt(	int		n,
									int		m,
	/* integer */				int	O[n][m]
									)
{
	int i,j;
	for(i=0; i<n; ++i){
		for(j=0; j<m; ++j){
			O[i][j]	=	0;
		}
	}
}
#endif

#ifndef F_MATXMAT
#define F_MATXMAT
/* matrix AB */
void MatXMat(	int		n,
					int		m,
					int		l,
					double	A[n][m],
					double	B[m][l],
					double	*C	/* (n,l) */
){
	int 		i,j,k;

	for(i=0; i<n; ++i){
		for(j=0; j<l; ++j){
			*(C+(l*i)+j) = 0;
			for(k=0; k<m; ++k){
				*(C+(l*i)+j) = *(C+(l*i)+j) + (A[i][k]*B[k][j]);
			}
		}
	}
}
#endif

#ifndef F_MATXVEC
#define F_MATXVEC
/* vector Ax */
void MatXVec(	int		n,
					int		m,
					double	A[n][m],
					double	x[m],
					double	*y	/* (n,1) */
){
	int	i;

	for(i=0; i<n; ++i){
		*(y+i) = innerproduct(A[i], x, m);
	}
}
#endif
#ifndef F_TRAMAT
#define F_TRAMAT
/* transpose A */
void TraMat(	int 		n,
					int		m,
					double	A[n][m],
					double	*T	/* (m,n) */
){
	int i,j;

	for(i=0; i<m; ++i){
		for(j=0; j<n; ++j){
			*(T+(n*i)+j) = A[j][i];
		}
	}
}
#endif

#ifndef F_CHANGEVEC
#define F_CHANGEVEC
/* change vector a <-> b */
void ChangeVec(	int		n,
						double	a[n],
						double	b[n]
){
	int i;	
	double buf[n];

	for(i=0; i<n; ++i){
		buf[i] = a[i];
		*(a+i) = b[i];
		*(b+i) = buf[i];
	}

}
#endif

#ifndef F_MATLINECOMB
#define F_MATLINECOMB
/* Compute matrix aA+bB */
void MatLineComb(	int	 	n,
						int		m,
						double	A[n][m],
						double	B[n][m],
						double	a,
						double	b,
						double	*C
){
	int	i,j;

	for(i=0; i<n; ++i){
		for(j=0; j<m; ++j){
			*(C+(m*i)+j) = a*A[i][j] + b*B[i][j];
		}
	}
}
#endif

#ifndef F_SCALARMAT
#define F_SCALARMAT
// aX, a:scalar
void ScalarMat(
	//	input
	int		n,
	int		m,
	double	X[n][m],
	double	a,
	//	output
	double	*Y
	)
{	
	int	i,j,ct=0;
	for(i=0; i<n; ++i){
		for(j=0; j<m; ++j){
		*(Y+ct)	=	a*X[i][j];
		ct++;
		}
	}
}
#endif

/* template */
/*
#ifndef F_
#define F_
*/
/* */
/*
#endif
*/

#endif

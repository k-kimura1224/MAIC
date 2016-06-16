/* cblapack.h */
/*********************************************/
/* supposition--------------------				*/
/* |  matirx A is 					|				*/
/* |		double, n * m, 			|				*/
/*	|		one dimensional array	|				*/
/* -------------------------------				*/
/* mat:				get A(i,j)						*/
/*	printM:			print Matrix					*/
/* printv:			print vector					*/
/*	dcopy:			double copy vecto				*/
/*	pickupcolumn:	generate new mat B			*/
/*	_dgelsy:			least suquares					*/
/* _dposv:			solve linear equation		*/
/* dgemv:			Ax + y							*/
/* dgemv_t:			(A^t)x							*/
/* dgemm_t:			(A^t)(A)							*/			
/* ddot:				dot product						*/
/* dscal:         (alpha)x							*/
/* daxpy:			(alpha)x + (beta)y			*/
/* dnrm2:			|| x ||							*/
/* dger:				(alpha)xy^t + A				*/
/*********************************************/
#ifndef CBLAPACK_H
#define CBLAPACK_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(__APPLE__)
#include "cblas.h"
#include "clapack.h"
#else
#define dposv dposv_
#define dgelsy dgelsy_
#define dcopy cblas_dcopy
#define dgemv cblas_dgemv
#define ddot cblas_ddot
#define dgemm cblas_dgemm
#define dscal cblas_dscal
#define daxpy cblas_daxpy
#define dnrm2 cblas_dnrm2
#define dger cblas_dger
#endif

#ifndef F_MAT_
#define F_MAT_
/* return A(i,j)	*/
/* A[n*m], Real	*/
double	mat_(
	const double	*A,
	const int		n,
	const int		i,
	const int		j
	)
{
	return *(A+i+(j*n));
}
#endif

#ifndef F_PRINTM_
#define F_PRINTM_
/* print Matrix	*/
/* A[n*m], Real 	*/
void	printM_(
	double	*A,
	int		n,
	int		m
	)
{
	int	i,j;

	printf("\n");

	for(i=0; i<n; ++i){
		for(j=0; j<m; ++j){
			printf(" %f,", mat_(A,n,i,j));
		}
		printf("\n");
	}
}
#endif

#ifndef F_PRINTV_
#define F_PRINTV_
/* print vector x[n] */
void printv_(
	double	*x,
	int		n
	)
{
	int i;

	printf("\n");

	for(i=0; i<n; ++i){
		printf(" %f,", *(x+i));
	}
	printf("\n");


}
#endif

#ifndef F_DCOPY_
#define F_DCOPY_
/* copy the vector */
/* y := x */
void dcopy_(
	const double	*x,
	double			*y,
	const int		n		/* size */
	)
{
	int inc=1;	/* increment */

	cblas_dcopy( n, x, inc, y, inc); 
}
#endif

#ifndef F_PICKUPCOLUMN
#define F_PICKUPCOLUMN
/* pick up column vectors of mat A */
/* generate new matrix B[n*l] */
void pickupcolumn_(
	const double	*A,	/* matrix A is n * m */
	const int		n,
	const int		m,
	const int		*x,	/* x[m], 1(if necessary) or 0(o.w.) */
	double			*B		/* B[n*l] */
	)
{
	int	i;
	int	ct=0;

	for(i=0; i<m; ++i){
		if( *(x+i)==1 ){
			dcopy_( A+(i*n), B+(ct*n), n);
			ct++;
		}
	}

}
#endif

#ifndef F__DGELSY_
#define F__DGELSY_
void _dgelsy_(
	int		n,	
	int		m,
	const	double	*y,	/* y[n] */
	const	double	*X,	/* X[n*m] */
	double			*a,	/* a[m] */
	int				*info	/* info = 0 if successful */
	)
{
	int one = 1;
	int i;
	
	if(n<m){
		printf("It is unexpected error n<m");
		exit(1);
	}

	double	XX[n*m];
	double	yy[n];
	int		JPVT[m];
	double	RCOND	=	1.0E-8;
	int		RANK;
	int		nb		=	50;
	int		LWORK	=	3*m + nb*(m+1);
	double	WORK[LWORK];

	dcopy_( X, XX, n*m);
	dcopy_( y, yy, n);
	for(i=0; i<m; ++i) JPVT[i]=0;
	
	dgelsy_( &n, &m, &one, XX, &n, yy, &n,
				JPVT, &RCOND, &RANK, WORK, &LWORK, info);

	if( LWORK<=(int)WORK[i] ){
		printf("error of dgelsy \n");
		printf("LWORK = %d", LWORK);
		printf("%d are necessary for optimal performance\n", (int)WORK[0]);
		stop();
	}
	
	dcopy_( yy, a, m);
}
#endif

#ifndef F__DPOSV_
#define F__DPOSV_
/* compute Ax=b with Cholesky decomposition */
/* A: positive definite symmetric matrix */
/* if successful, return 0 */
int _dposv_(
	const double	*A,		/* A[n*n] */
	const double	*b,		/* b[n]	*/
	int		n,			/* size */
	double			*x			/* solution */
	)
{
	char	uplo[1]	=	"L";
	int	one		=	1;
	int	info;
	double	AA[n*n];
	
	dcopy_( A, AA, n*n);
	dcopy_( b, x, n);

	dposv_( uplo, &n, &one, AA, &n, x, &n, &info);

	return info;

}
#endif

#ifndef F_DGEMV_
#define F_DGEMV_
/* z:= (alpha)Ax + (beta)y */
void dgemv_(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const double	*x,	/* x[m] */
	const double	*y,	/* y[n] */
	const double	al,	/* alpha */
	const double	be,	/* beta */
	double			*z		/* return z[n] */
	)
{
	int one=1;

	dcopy_( y, z, n);
	
	cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, al, A, n, x, one, be, z, one);
}
#endif

#ifndef F_DGEMV_T
#define F_DGEMV_T
/* z := (A^t)(x) */
void dgemv_t(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const	double	*x,	/* x[n] */
	double			*z		/* z[m] */
	)
{
	int one=1;

	cblas_dgemv(CblasRowMajor, CblasNoTrans,
					m, n, 1.0, A, n, x, one, 0.0, z, one);
}
#endif

#ifndef F_DGEMM_T
#define F_DGEMM_T
/* B := (A^t)(A) */
void dgemm_t(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	double			*B	/* B[m*m] */
	)
{
	double	zero=0.0;
	double	one=1.0;
	

	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
					m, m, n, one, A, n, A, n, zero, B, m);

}
#endif

#ifndef F_DDOT_
#define F_DDOT_
double ddot_(
	const double	*x,	/* vector */
	const double	*y,	/* vector */
	const int		n		/* size */
	)
{
	int 		one=1;
	double 	xy;

	xy = cblas_ddot( n, x, one, y, one);

	return xy;
}
#endif

#ifndef F_DSCAL_
#define F_DSCAL_
/* y := (alpha)x */
void dscal_(
	const double	*x,		/* vector */
	const int		n,			/* size */
	const double	alpha,	/* scalar */
	double			*y
	)
{
	int one = 1;
	dcopy_( x, y, n);

	cblas_dscal( n, alpha, y, one);

}
#endif

#ifndef F_DAXPY_
#define F_DAXPY_
/* z := (alpha)x + (beta)y */
void daxpy_(
	const double	*x,		/* vector */
	const double	*y,		/* vector */
	const int		n,			/* size */
	const double	alpha,	/* scalar */
	const double 	beta,		/* scalar */
	double			*z
	)
{
	int one=1;

	dscal_( y, n, beta, z);

	cblas_daxpy( n, alpha, x, one, z, one);

}
#endif

#ifndef F_DNRM2_
#define F_DNRM2_
/* return ||x|| */
double dnrm2_(
	const double	*x,		/* vector */
	const int		n			/* size */
	)
{
	int one = 1;
	return cblas_dnrm2( n, x, one);
}
#endif

#ifndef F_DGER_
#define F_DGER_
/* B := (alpha)xy^t + A */
void dger_(
	const double	*A,	/* A[n*m] */
	const double	*x,	/* x[n] */
	const double	*y,	/* y[m] */
	const int		n,		/* size */
	const int 		m,		/* size */
	const	double	alpha,/* scalar */
	double			*B		/* B[n*m} */
	)
{
	int one = 1;

	dcopy_( A, B, n*m);

	cblas_dger( CblasColMajor, n, m, alpha, x, one, y, one, B, n);
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

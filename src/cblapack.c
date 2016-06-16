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
/* dgemv_1:			Ax + y							*/
/* dgemv_2:			Ax 								*/
/* dgemv_t:			(A^t)x							*/
/* dgemm_t:			(A^t)(A)							*/			
/* ddot:				dot product						*/
/* dscal:         (alpha)x							*/
/* daxpy:			(alpha)x + (beta)y			*/
/* dnrm2:			|| x ||							*/
/* dger:				(alpha)xy^t + A				*/
/*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(__APPLE__)
#include "cblas.h"
#include "clapack.h"
#else
#include <cblas.h>
#ifdef __cpluscplus
extern "C" {
#endif
int dposv_( char* uplo, int* n, int* nrhs, 
	double* A, int* lda, double* x,
	int* ldb, int* info);
#ifdef __cpluscplus
}
#endif
#endif

#include "cblapack.h"

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

/* copy the vector */
/* y := x */
void mydcopy_(
	const double	*x,
	double			*y,
	const int		n		/* size */
	)
{
	int inc=1;	/* increment */

	cblas_dcopy( n, x, inc, y, inc); 
}

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
			mydcopy_( A+(i*n), B+(ct*n), n);
			ct++;
		}
	}

}

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
	char	uplo[2]	=	"L";
	int	one		=	1;
	int	info;
	double	*AA;

	/* alloc */
	AA	=	(double *)malloc( sizeof(double) * (n*n) );
	if( AA == NULL ){
		printf("error in cblack.c\n");
		exit(1);
	}
	
	mydcopy_( A, AA, n*n);
	mydcopy_( b, x, n);

	dposv_( uplo, &n, &one, AA, &n, x, &n, &info);

	/* free */
	free(AA);
	
	return info;

}

/* z:= (alpha)Ax + (beta)y */
void dgemv_1(
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

	mydcopy_( y, z, n);
	
	cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, al, A, n, x, one, be, z, one);
}

/* z:= Ax */
void dgemv_2(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const double	*x,	/* x[m] */
	double			*z		/* return z[n] */
	)
{
	int one=1;

	int i;
	for(i=0; i<m; i++){
		z[i] = 0;
	}
	
	cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, 1.0, A, n, x, one, 0.0, z, one);
}

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

	int i;
	for(i=0; i<m; i++){
		z[i] = 0;
	}

	cblas_dgemv(CblasRowMajor, CblasNoTrans,
					m, n, 1.0, A, n, x, one, 0.0, z, one);
}

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

double myddot_(
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

/* y := (alpha)x */
void mydscal_(
	const double	*x,		/* vector */
	const int		n,			/* size */
	const double	alpha,	/* scalar */
	double			*y
	)
{
	int one = 1;
	mydcopy_( x, y, n);

	cblas_dscal( n, alpha, y, one);

}

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

	mydscal_( y, n, beta, z);

	cblas_daxpy( n, alpha, x, one, z, one);

}

/* return ||x|| */
double dnrm2_(
	const double	*x,		/* vector */
	const int		n			/* size */
	)
{
	int one = 1;
	return cblas_dnrm2( n, x, one);
}

/* B := (alpha)xy^t + A */
void dger_1(
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

	mydcopy_( A, B, n*m);

	cblas_dger( CblasColMajor, n, m, alpha, x, one, y, one, B, n);
}

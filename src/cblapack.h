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
#ifndef CBLAPACK_H
#define CBLAPACK_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*
#if defined(__APPLE__)
#include "cblas.h"
#include "clapack.h"
#else
#include <cblas.h>
#endif
*/

/* return A(i,j)	*/
/* A[n*m], Real	*/
extern
double	mat_(
	const double	*A,
	const int		n,
	const int		i,
	const int		j
	);

/* print Matrix	*/
/* A[n*m], Real 	*/
extern
void	printM_(
	double	*A,
	int		n,
	int		m
	);

/* print vector x[n] */
extern
void printv_(
	double	*x,
	int		n
	);

/* copy the vector */
/* y := x */
extern
void mydcopy_(
	const double	*x,
	double			*y,
	const int		n		/* size */
	);

/* pick up column vectors of mat A */
/* generate new matrix B[n*l] */
extern
void pickupcolumn_(
	const double	*A,	/* matrix A is n * m */
	const int		n,
	const int		m,
	const int		*x,	/* x[m], 1(if necessary) or 0(o.w.) */
	double			*B		/* B[n*l] */
	);

/* compute Ax=b with Cholesky decomposition */
/* A: positive definite symmetric matrix */
/* if successful, return 0 */
extern
int _dposv_(
	const double	*A,		/* A[n*n] */
	const double	*b,		/* b[n]	*/
	int		n,			/* size */
	double			*x			/* solution */
	);

/* z:= (alpha)Ax + (beta)y */
extern
void dgemv_1(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const double	*x,	/* x[m] */
	const double	*y,	/* y[n] */
	const double	al,	/* alpha */
	const double	be,	/* beta */
	double			*z		/* return z[n] */
	);

/* z:= Ax */
void dgemv_2(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const double	*x,	/* x[m] */
	double			*z		/* return z[n] */
	);

/* z := (A^t)(x) */
extern
void dgemv_t(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const double	*x,	/* x[n] */
	double			*z		/* z[m] */
	);

/* B := (A^t)(A) */
extern
void dgemm_t(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	double			*B	/* B[m*m] */
	);

extern
double myddot_(
	const double	*x,	/* vector */
	const double	*y,	/* vector */
	const int		n		/* size */
	);

/* y := (alpha)x */
extern
void mydscal_(
	const double	*x,		/* vector */
	const int		n,			/* size */
	const double	alpha,	/* scalar */
	double			*y
	);

/* z := (alpha)x + (beta)y */
extern
void daxpy_(
	const double	*x,		/* vector */
	const double	*y,		/* vector */
	const int		n,			/* size */
	const double	alpha,	/* scalar */
	const double 	beta,		/* scalar */
	double			*z
	);

/* return ||x|| */
extern
double dnrm2_(
	const double	*x,		/* vector */
	const int		n			/* size */
	);

/* B := (alpha)xy^t + A */
extern
void dger_1(
	const double	*A,	/* A[n*m] */
	const double	*x,	/* x[n] */
	const double	*y,	/* y[m] */
	const int		n,		/* size */
	const int 		m,		/* size */
	const	double	alpha,/* scalar */
	double			*B		/* B[n*m} */
	);

#endif

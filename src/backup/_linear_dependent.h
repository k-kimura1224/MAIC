/* linear_dependent.h */
/*********************************************/
/*	linear_dependent for col.vec			 */
/*	generate linear_dependent group          */
/*********************************************/
#ifndef LINEAR_DEPENDENT_H
#define LINEAR_DEPENDENT_H

#include <stdio.h>
#include <stdlib.h>

#include "vector.h"
#include "cblapack.h"

//extern void		stop(void);
//extern void		SortingBubble(int n, double x[n], double *y, int *z);
///* function of cblapack.h */
//extern double	mat_(const double *A, const int n, const int i, const int j);
//extern void		printM_(double *A, int n, int m);
//extern void		printv_(double *x, int n);
//extern double	ddot_(const double *x, const double *y, const int n);
//extern void 	dscal_(const double *x, const int n,
//					const double alpha, double *y);
//extern void		daxpy_(const double *x, const double *y, const int n,
//					const double alpha, const double beta, double *z);
//extern double	dnrm2_(const double *x, const int n);
//extern void		dcopy_(const double *x, double *y, const int n);
//extern void 	dger_(const double *A, const double *x, const double *y,
//					const int n, const int m, const double alpha, double *B);
//extern void		dgemv_t(const double *A, const int n, const int m,
//					const double *x, double *z);
//extern void		dgemv_(const double *A, const int n, const int m,
//					const double *x, const double *y,
//					const double al, const double be, double *z);
//extern int		_dposv_(const double *A, const double *b,
//					int	n, double *x);


#ifndef F_LINEAR_DEPENDENT
#define F_LINEAR_DEPENDENT
/*	check linear_dependent for col.vec */	
void LinearDependent(	/* input */
								const int		n,
								const int		m,
								const double	*A,
								/* output  dep->1 , indep->0 */
								int		*result  /* result[m] */
/**
 * ex. A={a_1, a_2, a_3}
 *  a_1, a_2 is linear independent
 *  a_1, a_2, a_3, is linear dependent
 * -----> result = { 0, 0, 1 } 
**/
){
	int		i,j;
	int		memo = 0; 
	int		add  = 0;	
	int		rank = 1;				/* rank of sub */
	double	norm;						/* norm of b */
	double 	e	= 	10E-8;

	double	*sub		=	NULL;
	double	*inv		=	NULL;
	double	*x			=	NULL;	
	double	*b			=	NULL;		/* b := a_(j+1) - A_jx */
	double	*c			=	NULL;

	double	*sub_new		=	NULL;
	double	*inv_new		=	NULL;
	double	*x_new		=	NULL;
	double	*b_new		=	NULL;
	
	double	*P	=	NULL;
	double	*q	=	NULL;
	double	r;

	*result = 0;

	/**
	 * initialize 
	 *    sub := (a_1), inv = (A_1'A_1)-1,
	 *		x = inv sub a_i, b = a_i - sub x
	**/

	x	=	(double *)malloc(sizeof(double)*rank);
	b	=	(double *)malloc(sizeof(double)*n);
	if( x==NULL || b==NULL ){
		printf("malloc error\n");
		stop();
	}

	double a_0_norm = ddot_( A, A, n);

	for(i=1; i<m; ++i){
		memo = i;
		*x		=	ddot_( A, A+(i*n), n)/a_0_norm;
		daxpy_( A+(i*n), A, n, 1.0, -(*x), b);
		norm = dnrm2_(b, n);
		
		if( norm < e ){
			*(result+i) = 1;
		}else{
			*(result+i) = 0;
			break;
		}
	}
	
	sub	=	(double *)malloc(sizeof(double)*(n*rank));		/* rank=1 */
	inv	=	(double *)malloc(sizeof(double)*(rank*rank));
	b_new	=	(double *)malloc(sizeof(double)*n);
	if( sub==NULL || inv==NULL || b_new==NULL ){
		printf("malloc error\n");
		stop();
	}
	
	dcopy_( A, sub, n);
	*inv = 1/a_0_norm;

	/*********************************************/
	/**
	 *   (subA_new'subA_new)^-1
	 *		=	( P		q )
	 			( q^t		r )
	**/

	add = 1;

	for(i=(memo+1); i<m; ++i){
		
		if( add ){
			/* update 1 */
			P	=	(double *)malloc(sizeof(double)*(rank*rank));
			q	=	(double *)malloc(sizeof(double)*rank);
			rank++;
			sub_new	=	(double *)malloc(sizeof(double)*(n*rank));
			inv_new	=	(double *)malloc(sizeof(double)*(rank*rank));
			x_new		=	(double *)malloc(sizeof(double)*rank);
			c			=	(double *)malloc(sizeof(double)*rank);
			if( sub_new==NULL || inv_new==NULL || x_new==NULL || c==NULL ){
				printf("malloc error\n");
				stop();
			}
			
			dcopy_( sub, sub_new, n*(rank-1));
			dcopy_( A+(memo*n), sub_new+(n*(rank-1)), n);
			add = 0;
		}
		
		/* r, q, P*/
		r = 1/ddot_( A+(memo*n), b, n);
		dscal_( x, rank-1, -r, q);
		dger_( inv, x, q, rank-1, rank-1, -1.0, P);
		
		/* generate inv_new from r,q,P */
		for(j=0; j<rank-1; ++j){
			dcopy_( P+(j*(rank-1)), inv_new+(j*rank), rank-1);
			inv_new[(rank-1)+(j*rank)] = q[j];
		}
		dcopy_( q, inv_new+(rank*(rank-1)), rank-1);
		inv_new[(rank*rank)-1] = r;
		
		/**
		 *	compute 
		 *	c		=	sub_new^t a[i] 
		 *	x_new	=	inv_new c
		 * b_new	=	a[i] - sub_new x_new
		**/
		dgemv_t( sub_new, n, rank, A+(i*n), c);
		dgemv_( inv_new, rank, rank, c, c, 1.0, 0.0, x_new);
		dgemv_( sub_new, n, rank, x_new, A+(i*n), -1.0, 1.0, b_new);

		/* compute norm of b_new */
		norm = dnrm2_( b_new, n);

		if( norm < e ){	/* if linear dep */
			*(result+i) = 1;
		}else{
			*(result+i) = 0;
			memo	= 	i;
			add	=	1;

			/* update 2 */
			if( i<(m-1) ){
				free(inv);
				free(x);
				free(sub);
			
				inv	=	(double *)malloc(sizeof(double)*(rank*rank));
				x		=	(double *)malloc(sizeof(double)*rank);
				sub	=	(double *)malloc(sizeof(double)*(n*rank));

				if( inv==NULL || x==NULL || sub==NULL ){
					printf("malloc error\n");
					stop();
				}

				dcopy_( inv_new, inv, rank*rank);	/* inv	<-	inv_new */
				dcopy_( x_new, x, rank);				/* x		<-	x_new */	
				dcopy_( sub_new, sub, n*rank);		/* sub	<- sub_new */
				dcopy_( b_new, b, n);					/* b		<-	b_new */

				free(inv_new);
				free(x_new);
				free(sub_new);
				free(P);
				free(q);
				free(c);
			}
		}
	}

	

	free(sub);
	free(inv);
	free(x);
	free(b);
	free(sub_new);
	free(inv_new);
	free(x_new);
	free(b_new);
	free(P);
	free(q);
	free(c);
}
#endif

#ifndef F_LINEAR_DEPENDENTGROUP
#define F_LINEAR_DEPENDENTGROUP
/*	generate linear_dependent group */
void LinearDependentGroup(	/* input */
								const int		m,
								const double	*Q,	/* A^t A [m*m] */
								const int 		ndep,
								const int		Mdep[ndep],
								const int		D[m],
								int				groupX[ndep][m]
/**
 * ex. A={a_1, a_2, a_3}, ndep=1
 *  a_1, a_2 is linear independent
 *  a_1, a_2, a_3, is linear dependent
 * -----> groupX[0] = { 1, 1, 1 } 
**/
){
	int	i,j,k,t;
	int	ct1	=	0;
	int	ct2	=	0;
	int	rank	=	1;

	double 	e	= 	10E-8;

	double	*subQ	=	NULL;
	double	*d		=	NULL;
	double	*x		=	NULL;
	double	*sort	=	NULL;
	int		*num	=	NULL;
	

	for(i=1; i<m; ++i){
		if( i==Mdep[ct1] && rank>1 ){
		/**
		 * solve:
		 * 	subQ x = d
		 *		
		 *		subQ 	:= subA'subA
		 *		d		:= subA'a_i
		**/
			
			/* malloc */
			subQ	=	(double *)malloc(sizeof(double)*(rank*rank));
			d		=	(double *)malloc(sizeof(double)*rank);
			x		=	(double *)malloc(sizeof(double)*rank);
			sort	=	(double *)malloc(sizeof(double)*rank);
			num	=	(int *)malloc(sizeof(int)*rank);
			if( subQ==NULL || d==NULL || x==NULL ){
				printf("malloc error\n");
				stop();
			}
			
			/* generate subQ */
			ct2 = 0;
			for(j=0; j<m; ++j){
				if( D[j]==0 ){
					for(k=0; k<m; ++k){
						if( D[k]==0 )	subQ[ct2++] = mat_( Q, m, j, k);
						else if( k==Mdep[ct1] ) break;
					}
				}else if( j==Mdep[ct1] ) break;
			}

			/* generate d */
			ct2 = 0;
			for(j=0; j<m; ++j){
				if( D[j]==0 )	d[ct2++] = mat_( Q, m, j, i);
				else if( j==Mdep[ct1] ) break;
			}

			/* solve subQ x = d with dposv */
			_dposv_( subQ, d, rank, x);
			
			/* square */
			for(j=0; j<rank; ++j) x[j] *= x[j];
			
			/* sort */
			SortingBubble( rank, x, sort, num);

			/* check whether 0 */
			for(j=0; j<rank; ++j){
				if( sort[j+1]/sort[j] < e ) break;
			}

			/* groupX <- 1 */
			for(k=0; k<(j+1); ++k){
				ct2=0;
				for(t=0; t<m; ++t){
					if( D[t]==0 ){
						if( num[k]==ct2 ){
							groupX[ct1][t] = 1;
							break;
						}
						ct2++;
					}
				}
			}

			ct1++;
			/* free */
			free(subQ);
			free(d);
			free(x);
			free(sort);
			free(num);
		}else if( i==Mdep[ct1] && rank==1 ){
			groupX[ct1][0] = 1;
			ct1++;
		}else{
			rank++;
		}
	}
}
#endif

#ifndef F_PRINTLINEDEPGROUP
#define F_PRINTLINEDEPGROUP
void printLineDepGroup(
	int	ndep,
	int	p,
	int	groupX[ndep][p]
	)
{
	int i,j;

	if( ndep==0 )
		printf("linear independent\n");
	else
		printf("linear dependent sets\n");

	for(i=0; i<ndep; ++i){
		printf(" X_%d = { ", i+1);
		for(j=0; j<p; ++j){
			if( groupX[i][j]==1 )
				printf("%d ", j+1);
		}
		printf("}\n");
	}
	printf("\n");
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

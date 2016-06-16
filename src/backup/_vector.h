/* vector.h */
/*********************************************/
/*		print matrix(x,y)								*/
// 	print integer matrix
/*		starat a new line 							*/
/* 	line												*/
/*		print scalar a 								*/
//		print integer number
/* 	print vector(x)								*/
/* 	print int vector(x)							*/
/*		vector sum										*/
/*		vector sum (integer)							*/
/*		generate zero vector 						*/
/*		generate zero vector (integer)			*/
/*		inner product									*/
/*		vector ax+yb									*/
// 	vector ax = z scalar:a
//		check 0 or 1
/* 	check vector whether zero or nonzero	*/
/*		compute norm of vector						*/
//		compute unit vector
/*		compute ^2 for vector`s elements			*/
/*		bubble sort										*/
/*		print from 1 to n								*/
/*		stop												*/
/*		error												*/
/*********************************************/
#ifndef BASE_H
#define BASE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#if defined(__APPLE__)
//#include "cblas.h"
//#endif

#ifndef F_PRINTM
#define F_PRINTM
/* print matrix(x,y) */
void printM(	int		x,
					int		y,
					double	Mat[x][y]
){
	int i;
	int j;
	printf("\n");
	for(i = 0; i < x; i++){
		for(j = 0; j < y; j++){
			printf(" %f,",Mat[i][j]);
		}
		printf("\n");
	}
}
#endif

#ifndef F_PRINTNM
#define F_PRINTNM
/* print matrix(x,y) , integer */
void printNM(	int		x,
					int		y,
					int		Mat[x][y]
){
	int i;
	int j;
	printf("\n");
	for(i = 0; i < x; i++){
		for(j = 0; j < y; j++){
			printf(" %d ,",Mat[i][j]);
		}
		printf("\n");
	}
}
#endif


#ifndef F_NEWLINE
#define F_NEWLINE
/* start a new line */
void newline(void){
	printf("\n");
}
#endif

#ifndef F_LONGLINE
#define F_LOMGLINE
/* line */
void Longline(void){
	printf("*------------------------------------------------------------*\n");
}
#endif

#ifndef F_PRINTS
#define F_PRINTS
/* print scalar a */
void prints(	double	a
){
	newline();
	printf("%f", a);
	newline();
}
#endif

#ifndef F_PRINTN
#define F_PRINTN
// print integer number
void printN(	int	a
){
	newline();
	printf("%d", a);
	newline();
}

#endif

#ifndef F_PRINTV
#define F_PRINTV
/* print vector(x) */
void printv(	int		x,
					double	vec[x]
){
	int i;
	newline();
	for(i = 0; i < x; i++){
		printf("%f ,", vec[i]);
	}
	newline();
}
#endif

#ifndef F_PRINTINTV
#define F_PRINTINTV
/* print int vector(x) */
void printintv(	int		x,
					int		vec[x]
){
	int i;
	newline();
	printf(" ");
	for(i = 0; i < x; i++){
		printf("%d , ", vec[i]);
	}
	newline();
}
#endif



#ifndef F_SUM
#define F_SUM
/* vector sum */
double sum(	double 	x[],
				int		i
){
	int j;
	double y=0.0;
	for(j = 0; j < i; j++){
      y = y + x[j];
   }
   return y;
}
#endif

#ifndef F_SUMINT
#define F_SUMINT
/* vector sum (integer) */
int sumint(	int x[],
				int i
){
	int j;
	int y=0;
	for( j = 0; j < i; j++){
		y = y + x[j];
	}
	return y;
}
#endif

#ifndef F_GENERATEZEROVEC
#define F_GENERATEZEROVEC
/* generate zero vector */
void GenerateZeroVec(	int		n,
								double	x[n]
){
	int i;
	for(i=0; i<n; ++i){
		x[i] = 0;
	}
}
#endif

#ifndef F_GENERATEZEROVECINT
#define F_GENERATEZEROVECINT
/* generate zero vector */
void GenerateZeroVecInt(	int		n,
									int		x[n]
){
	int i;
	for(i=0; i<n; ++i){
		x[i] = 0;
	}
}
#endif

#ifndef F_NAISEKI
#define F_NAISEKI
/* inner product */
double innerproduct(	double	a[],
						double	b[],
						int		x
){
	int i;
	double c[x];
	for(i = 0; i < x; i++){
		c[i] = a[i] * b[i];
	}
	return sum(c, x);
}
#endif

#ifndef F_VECLINECOMB
#define F_VECLINECOMB
/* vector ax + by = z */
void VecLineComb(	int		n,
						double	x[n],
						double	y[n],
						double	a,
						double	b,
						double	*z
){
	int i;
	for(i=0; i<n; ++i){
		*(z+i) = a*x[i] + b*y[i];
	}
}

#endif

#ifndef F_SCALAR_V
#define F_SCALAR_V
// vector ax = z scalar:a
void Scalar_v(
	int		n,
	double	x[n],
	double	a,
	double	*z
	)
{
	int i;
	for(i=0; i<n; ++i){
		*(z+i) = a*x[i];
	}
}
#endif

#ifndef F_STOP
#define F_STOP
/* stop */
void stop(void){
	newline();
	printf("stop!");
	newline();
	exit(1);
}
#endif

#ifndef F_CHECK01
#define F_CHECK01
// 0 --> 0, o.w. --> 1
int Check01(
	double x
	)
{
	if( x*x < 0.000000000001 )
		return 0;
	else
		return 1;
}
#endif

#ifndef F_NORMOFVEC
#define F_NORMOFVEC
/* norm of vector x */
double NormOfVec(	int		n,
						double	x[n]
){
	int		i;
	double	norm;
	
	norm = 0;
	
	for(i=0; i<n; i++){
		norm += x[i]*x[i];
	}

	return sqrt(norm);
}
#endif

#ifndef F_CHECKVEC01
#define F_CHECKVEC01
int CheckVec01(	int		n,
						double	x[n]
){
	double	norm = NormOfVec( n, x); 
	
	if( norm < 0.00000001 )
		return 0;
	else
		return 1;
}
#endif

#ifndef F_UNITVECTOR
#define F_UNITVECTOR
//		compute unit vector
void UnitVector(
// input
	int		n,
	double	x[n],
// output
	double	*e
	)
{
	double norm = NormOfVec( n, x);
	
	Scalar_v( n, x, 1/norm, e);
}
#endif

#ifndef F_FINDMINCOMPVEC
#define F_FINDMINCOMPVEC
/* find minimal component of vector */
double FindMinCompVec(	
	// input
	int		n,
	double	x[n],
	// output
	int		*num
	)
{
	int		i;
	int		y_num;
	double	y;

	y		= x[0];
	y_num	= 0;

	for(i=1; i<n; ++i){
		if( y >= x[i] ){
			y 		= x[i];
			y_num	= i;
		}
	}
	
	*num 	= y_num;
	return y;
}
#endif

#ifndef F_FINDMININDEXVEC
#define F_FINDMININDEXVEC
/* find index of minimal component of vector */
int FindMinIndexVec(	int		n,
							double	x[n]
){
	int		i;
	double	y;
	int		z;

	y = x[0];
	z = 0;

	for(i=1; i<n; ++i){
			if( y >= x[i] ){
				y = x[i];
				z = i;
			}
	}

	return z;
}
#endif

#ifndef F_SQUAREELE
#define F_SQUAREELE
/*	compute ^2 for vector`s elements	*/
void SquareEle(
	// input 
	int 		n,
	double	a[n],
	// return
	double	*b
	)
{
	int i;

	for(i=0; i<n; ++i)
		*(b+i)	=	a[i]*a[i];
}
#endif

#ifndef F_SORTINGBUBBLE
#define F_SORTINGBUBBLE
/* bubble sort */
void SortingBubble(
	// input 
	int 		n,
	double	x[n],
	// output
	double	*y,
	int		*z
	)
{
	int i,j;

	double	a[n];
	int		b[n];

	double 	buf1;
	int		buf2;

	for(i=0; i<n; ++i){
		a[i]	=	x[i];
		b[i]	=	i;
	}

	// bubble sort
	for(i=0; i<(n-1); ++i){
		for(j=(n-1); j>i; --j){
			if( a[j]>a[j-1] ){
				buf1	=	a[j];
				a[j]	=	a[j-1];
				a[j-1]=	buf1;

				buf2	=	b[j];
				b[j]	=	b[j-1];
				b[j-1]=	buf2;
			}
		}
	}

	for(i=0; i<n; ++i){
		*(y+i) = a[i];
		*(z+i) = b[i];
	}
}
#endif

#ifndef F_ONETON
#define F_ONETON
// print from 1 to n
void oneton(
	int n
	)
{
	int i;

	printf(" ");
	for(i=0; i<n; ++i){
		if( i<9 )
			printf("%d , ", i+1);
		else
			printf("%d, ", i+1);
	}
	printf("\n");
}
#endif

//flag
#ifndef F_ERROR
#define F_ERROR
/* error */
int error(	int n){
	if( n==0 ){
		printf("\nerror\n");
		return -1;
	}else{
		printf("\nerror-%d\n",n);
		n += 1;
		return -1;
	}
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

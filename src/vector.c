/* vector.h */
/*********************************************/
/*		print matrix(x,y)								*/
/* 	print integer matrix 						*/
/*		starat a new line 							*/
/* 	line												*/
/*		print scalar a 								*/
/*		print integer number 						*/
/* 	print vector(x)								*/
/* 	print int vector(x)							*/
/*		vector sum										*/
/*		vector sum (integer)							*/
/*		generate zero vector 						*/
/*		generate zero vector (integer)			*/
/*		inner product									*/
/*		vector ax+yb									*/
/* 	vector ax = z scalar:a						*/
/*		check 0 or 1									*/
/* 	check vector whether zero or nonzero	*/
/*		compute norm of vector						*/
/*		compute unit vector							*/
/*		compute ^2 for vector`s elements			*/
/*		bubble sort										*/
/*		print from 1 to n								*/
/*		stop												*/
/*		error												*/
/*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "vector.h"

/* print matrix(x,y) */
void printM(	int		x,
					int		y,
					double	*Mat
){
	int i,j;
	for(i = 0; i < x; i++){
		for(j=0; j<y; j++){
			printf(" %f,", *(Mat+(i*y)+j));
		}
		printf("\n");
	}
}

/* print matrix(x,y) , integer */
void printNM(	int		x,
					int		y,
					int		*Mat
){
	int i;
	int j;
	for(i = 0; i < x; i++){
		for(j = 0; j < y; j++){
			printf(" %d,", *(Mat+(i*y)+j));
		}
		printf("\n");
	}
}

/* start a new line */
void newline(void){
	printf("\n");
}

/* line */
void Longline(void){
	printf("*------------------------------------------------------------*\n");
}

/* print scalar a */
void prints(	double	a
){
	newline();
	printf("%f", a);
	newline();
}

/* print integer number */
void printN(	int	a
){
	printf("%d", a);
	newline();
}

/* print vector(x) */
void printv(	int		x,
					double	*vec
){
	int i;
	for(i = 0; i < x; i++){
		printf("%f,", vec[i]);
	}
	newline();
}

/* print int vector(x) */
void printintv(	int		x,
						int		*vec
){
	int i;
	newline();
	printf(" ");
	for(i = 0; i < x; i++){
		printf("%d , ", vec[i]);
	}
	newline();
}

/* vector sum */
double sum(	double 	*x,
				int		i
){
	int j;
	double y=0.0;
	for(j = 0; j < i; j++){
      y = y + x[j];
   }
   return y;
}

/* vector sum (integer) */
int sumint(	int *x,
				int i
){
	int j;
	int y=0;
	for( j = 0; j < i; j++){
		y = y + x[j];
	}
	return y;
}

/* generate zero vector */
void GenerateZeroVec(	int		n,
								double	*x
){
	int i;
	for(i=0; i<n; ++i){
		x[i] = 0;
	}
}

/* generate zero vector */
void GenerateZeroVecInt(	int		n,
									int		*x
){
	int i;
	for(i=0; i<n; ++i){
		x[i] = 0;
	}
}

/* vector ax + by = z */
void VecLineComb(	int		n,
						double	*x,
						double	*y,
						double	a,
						double	b,
						double	*z
){
	int i;
	for(i=0; i<n; ++i){
		*(z+i) = a*x[i] + b*y[i];
	}
}

/* vector ax = z scalar:a */
void Scalar_v(
	int		n,
	double	*x,
	double	a,
	double	*z
	)
{
	int i;
	for(i=0; i<n; ++i){
		*(z+i) = a*x[i];
	}
}

/* stop */
void stop(void){
	newline();
	printf("stop!");
	newline();
	exit(1);
}

/* 0 --> 0, o.w. --> 1 */
int Check01(
	double x
	)
{
	if( x*x < 0.000000000001 )
		return 0;
	else
		return 1;
}

/* norm of vector x */
double NormOfVec(	int		n,
						double	*x
){
	int		i;
	double	norm;
	
	norm = 0;
	
	for(i=0; i<n; i++){
		norm += x[i]*x[i];
	}

	return sqrt(norm);
}

int CheckVec01(	int		n,
						double	*x
){
	double	norm = NormOfVec( n, x); 
	
	if( norm < 0.00000001 )
		return 0;
	else
		return 1;
}

/*	compute unit vector */
void UnitVector(
/* input */
	int		n,
	double	*x,
/* output */
	double	*e
	)
{
	double norm = NormOfVec( n, x);
	
	Scalar_v( n, x, 1/norm, e);
}

/* find minimal component of vector */
double FindMinCompVec(	
	/* input */
	int		n,
	double	*x,
	/* output */
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

/* find index of minimal component of vector */
int FindMinIndexVec(	int		n,
							double	*x
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

/*	compute ^2 for vector`s elements	*/
void SquareEle(
	/* input */
	int 		n,
	double	*a,
	/* output */
	double	*b
	)
{
	int i;

	for(i=0; i<n; ++i)
		*(b+i)	=	a[i]*a[i];
}

/* bubble sort */
void SortingBubble(
	/* input */
	int 		n,
	double	*x,
	/* output */
	double	*y,
	int		*z
	)
{
	int i,j;

	double	*a;
	int		*b;

	double 	buf1;
	int		buf2;

	/* alloc */
	a	=	(double *)malloc(sizeof(double) * n );
	b	=	(int *)malloc(sizeof(int) * n );

	if( ( a == NULL ) || ( b == NULL ) ){
		printf("error in vector.c\n");
		exit(1);
	}

	for(i=0; i<n; ++i){
		a[i]	=	x[i];
		b[i]	=	i;
	}

	/* bubble sort */
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

	/* free */
	free(a);
	free(b);
}

/* print from 1 to n */
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

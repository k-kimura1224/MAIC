/* vector.h */
/*********************************************/
/*		print matrix(x,y)								*/
/* 	print integer matrix							*/
/*		starat a new line 							*/
/* 	line												*/
/*		print scalar a 								*/
/*		print integer number							*/
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
#ifndef BASE_H
#define BASE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* print matrix(x,y) */
extern
void printM(	int		x,
					int		y,
					double	*Mat
);

/* print matrix(x,y) , integer */
extern
void printNM(	int		x,
					int		y,
					int		*Mat
);

/* start a new line */
extern
void newline(void);

/* line */
extern
void Longline(void);

/* print scalar a */
extern
void prints(double a);

/* print integer number */
extern
void printN(int a);

/* print vector(x) */
extern
void printv(	int		x,
					double	*vec
);

/* print int vector(x) */
extern
void printintv(	int		x,
						int		*vec
);

/* vector sum */
extern
double sum(	double 	*x,
				int		i
);

/* vector sum (integer) */
extern
int sumint(	int *x,
				int i
);

/* generate zero vector */
extern
void GenerateZeroVec(	int		n,
								double	*x
);

/* generate zero vector */
extern
void GenerateZeroVecInt(	int		n,
									int		*x
);

/* vector ax + by = z */
extern
void VecLineComb(	int		n,
						double	*x,
						double	*y,
						double	a,
						double	b,
						double	*z
);

/* vector ax = z scalar:a */
extern
void Scalar_v(
	int		n,
	double	*x,
	double	a,
	double	*z
	);

/* stop */
extern
void stop(void);

/* 0 --> 0, o.w. --> 1*/ 
extern
int Check01(
	double x
	);

/* norm of vector x */
extern
double NormOfVec(	int		n,
						double	*x
);

extern
int CheckVec01(	int		n,
						double	*x
);

/*	compute unit vector */
extern
void UnitVector(
/* input */
	int		n,
	double	*x,
/* output */
	double	*e
	);

/* find minimal component of vector */
extern
double FindMinCompVec(	
/* input */
	int		n,
	double	*x,
/* output */
	int		*num
	);

/* find index of minimal component of vector */
extern
int FindMinIndexVec(	int		n,
							double	*x
);

/*	compute ^2 for vector`s elements	*/
extern
void SquareEle(
	/* input */
	int 		n,
	double	*a,
	/* return */
	double	*b
	);

/* bubble sort */
extern
void SortingBubble(
	/* input */
	int 		n,
	double	*x,
	/* output */
	double	*y,
	int		*z
	);

/* print from 1 to n */
extern
void oneton(
	int n
	);

/* error */
extern
int error(	int n);

#endif

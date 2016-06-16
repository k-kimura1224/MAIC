/* linear_dependent.h */
/*********************************************/
/*	linear_dependent for col.vec			 */
/*	generate linear_dependent group          */
/*********************************************/
#ifndef LINEAR_DEPENDENT_H
#define LINEAR_DEPENDENT_H

#include <stdio.h>
#include <stdlib.h>

/*	check linear_dependent for col.vec */	
extern
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
);

/*	generate linear_dependent group */
extern
void LinearDependentGroup(	
								/* input */
								const int		m,
								const double	*Q,		/* A^t A [m*m] */
								const int 		ndep,
								const int		*Mdep,	/* [ndep] */
								const int		*D,		/* [m] */
								/* output */
								int				*groupX	/* [ndep*m] */
/**
 * ex. A={a_1, a_2, a_3}, ndep=1
 *  a_1, a_2 is linear independent
 *  a_1, a_2, a_3, is linear dependent
 * -----> groupX[0] = { 1, 1, 1 } 
**/
);

extern
void printLineDepGroup(
	int	ndep,
	int	p,
	int	*groupX	/* [ndep*p] */
	);

#endif

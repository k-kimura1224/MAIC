/* normalization.h */
/*********************************************/
/*	generate mean										*/
/* generate variance									*/
/* normalization										*/
/*********************************************/
#ifndef NORMALIZATION_H
#define NORMALIZATION_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* generate mean */
extern
void GeneMean(	int		n,
					int		p,
					double	*data, /* [n*(p+1)] */
					double	*m		 /* [p+1] */
					);

/* generate variance */
extern
void GeneVar(	int		n,
					int		p,
					double	*data,	/* [n*(p+1)] */
					double	*m,		/* [p+1] */
					double	*v			/* [p+1] */
					);

/* normalization */
extern
void Normalization(	int		n,
							int		p,
							double	*data	/* [n*(p+1)] */
	);

#ifdef __cplusplus
}
#endif

#endif

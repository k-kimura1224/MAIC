/* normalization.h */
/*********************************************/
/*	generate mean										*/
/* generate variance									*/
/* normalization										*/
/*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "normalization.h"

/* generate mean */
void GeneMean(	int		n,
					int		p,
					double	*data, /* [n*(p+1)] */
					double	*m		 /* [p+1] */
					)
{
	int i,j;

	for(i=0; i<p+1; ++i){
		*(m+i) = 0.0;
		for(j=0; j<n; ++j){
			*(m+i) += *(data+(j*(p+1))+i);
		}
		*(m+i) *= 1/(double)n;
	}
}

/* generate variance */
void GeneVar(	int		n,
					int		p,
					double	*data,	/* [n*(p+1)] */
					double	*m,		/* [p+1] */
					double	*v			/* [p+1] */
					)
{
	int i,j;
	
	for(i=0; i<p+1; ++i){
		*(v+i) = 0.0;
		for(j=0; j<n; ++j){
			*(v+i) += pow((*(data+(j*(p+1))+i)-m[i]),2.0);
		}
		*(v+i) *= 1/(((double)n)-1.0);
	}
}

/* normalization */
void Normalization(	int		n,
							int		p,
							double	*data	/* [n*(p+1)] */
){
	int 		i,j;
	double	*m;
	double	*v;

	m	=	(double *)malloc( sizeof(double) * (p+1));
	v	=	(double *)malloc( sizeof(double) * (p+1));

	if( ( m == NULL ) || ( v == NULL ) ){
		printf("error in normalization.c\n");
		exit(1);
	}

	GeneMean( n, p, data, m);
	GeneVar(  n, p, data, m, v);

	for(i=0; i<n; i++){
		for(j=0; j<p+1; j++){
			*(data+(i*(p+1))+j) = (*(data+(i*(p+1))+j)-m[j])/sqrt(v[j]);
		}
	}

	free(m);
	free(v);
}


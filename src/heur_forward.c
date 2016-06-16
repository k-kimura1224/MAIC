
/**@file   heur_forward.h
 * @brief  Forward selection
 *
 * Forward selection involves starting with no variables in the model.
 * this algorithm calculates AIC of the model that increased each variable,
 * selects adding the variable that improves the model the most,
 * and repeat this process until none improves the model.
 * 
 *
 */
#include <assert.h>
#include <string.h>
#include <math.h>

#include "heur_forward.h"
#include "probdata_linereg.h"
#include "matrix.h"
#include "vector.h"
#include "set_myparameter.h"
#include "cblapack.h"

#define HEUR_NAME             "forward"
#define HEUR_DESC             "primal heuristic using forward selection"
#define HEUR_DISPCHAR         'f'
#define HEUR_PRIORITY         1000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH        	10
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE

#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define MYPARA_LOG				0
/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
	int	a;
};


/*
 * Local methods
 */
static
SCIP_Real RSSvalue(
	int				dim,
	SCIP_Real		*a,	/* [dim] */
	SCIP_Real		*Xy,	/* [dim] */
	SCIP_Real		r
	)
{
	return - myddot_( a, Xy, dim) + r;
}

static SCIP_Real AICvalue(
	int				n,
	int				dim,
	SCIP_Real		RSS
	)
{
	return (double)n * log(RSS) + 2.0 * (double)dim;
}

static SCIP_RETCODE	Calculate_a(
   SCIP*          scip, /**< SCIP data structure */
	int				n,
	int				dim,
	SCIP_Real*		X_,	/* [n*(dim-1)] */
	SCIP_Real*		x,		/* x_i [n]*/
	SCIP_Real*		Y_,	/* [(dim-1)*(dim-1)] */
	SCIP_Real*		Xy_new,/* [dim] */
	/* return */
	SCIP_Real*		a,		/* [dim] */
	SCIP_Real*		W_		/* [dim*dim] */
	)
{
	SCIP_Real*	b;
	SCIP_Real*	c;
	SCIP_Real*	d;
	SCIP_Real	u;
	SCIP_Real*	v;
	SCIP_Real*	V_;
	SCIP_Real*	Z_;

	int	j,t;
	/* alloc */
	SCIP_CALL( SCIPallocBufferArray( scip, &b, dim-1));
	SCIP_CALL( SCIPallocBufferArray( scip, &c, dim-1));
	SCIP_CALL( SCIPallocBufferArray( scip, &d, n));
	SCIP_CALL( SCIPallocBufferArray( scip, &v, dim-1));
	SCIP_CALL( SCIPallocBufferArray( scip, &V_, (dim-1)*(dim-1)));
	SCIP_CALL( SCIPallocBufferArray( scip, &Z_, (dim)*(dim)));

	/* 1. b <- X^t x_i */
	dgemv_t( X_, n, dim-1, x, b);

	/* 2. c <- Y b */
	dgemv_2( Y_, dim-1, dim-1, b, c);

	/* 3. d <- - X c + x_i */
	dgemv_1( X_, n, dim-1, c, x, -1.0, 1.0, d);

	/* 4. u <- 1/<x_i, d> */
	u = 1.0 / myddot_( x, d, n);

	/* 5. v <- - u c */
	mydscal_( c, dim-1, -u, v);

	/* 6. V <- Y + u c c^t */
	dger_1( Y_, c, c, dim-1, dim-1, u, V_);

	/* 7. Z */
	/* V */
	for(j=0; j<(dim-1); j++){
		for(t=0; t<(dim-1); t++){
			*(Z_ + j + (t*dim) ) = mat_( V_, dim-1, j, t);
		}
	}
	/* v */
	for(j=0; j<(dim-1); j++){
		*(Z_ + dim-1 + (j*dim) )  = v[j];
		*(Z_ + j + ((dim-1)*dim)) = v[j];
	}

	/* u */
	*(Z_ + dim-1 + ((dim-1)*dim)) = u;

	/* 8. a_old <- Z (Xy) */
	dgemv_2( Z_, dim, dim, Xy_new, a);

	for(j=0; j<(dim*dim); j++){
		W_[j] = Z_[j];
	}

	SCIPfreeBufferArray( scip, &b);
	SCIPfreeBufferArray( scip, &c);
	SCIPfreeBufferArray( scip, &d);
	SCIPfreeBufferArray( scip, &v);
	SCIPfreeBufferArray( scip, &V_);
	SCIPfreeBufferArray( scip, &Z_);

	return SCIP_OKAY;
}

	
/*
 * Callback methods of primal heuristic
 */


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyForward)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurForward(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeForward)
{   /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecForward)
{  /*lint --e{715}*/
	
   SCIP_PROBDATA* probdata;
	int	n;
	int	p;
	int	ndep;

	/* "_" means the matrix for blas */
	SCIP_Real*	y;				/* [n] */
	SCIP_Real*	orig_X_;		/* [n*p] */
	SCIP_Real*	orig_Q_;		/* [p*p] <- (X^t) X */
	SCIP_Real*	orig_q;		/* [p]   <- (X^t) y */
	SCIP_Real	r;

	int*	Mdep;					/* [ndep] */
	int*	groupX;				/* [ndep*p] */

	int	dim;
	int*	list;					/* [p] */

	SCIP_Real	RSS;
	SCIP_Real	RSS_new;

	SCIP_Real	AIC = 1e+06;
	SCIP_Real	AIC_new;

	SCIP_Real*	a_old;		/* value of sol */
	SCIP_Real*	a;		/* value of sol */
	SCIP_Real*	a_new;		/* value of sol */
	SCIP_Real*	Xy;	/* sub vector of orig_q */ 
	SCIP_Real*	Xyn;	/* Xy_new */
	SCIP_Real*	X_;	
	SCIP_Real*	Y_;	
	SCIP_Real*	W_;	
	SCIP_Real*	Wnew_;	

	/* variables */
	SCIP_VAR**	var_a;		/* [p] continuous variables */
	SCIP_VAR**	var_z;		/* [p] 01 variables */
	SCIP_VAR**	var_ep;		/* [n] continuous variables */
	SCIP_VAR*	var_rss;		/* continuous variable, residual sum of squares */
	SCIP_VAR*	var_log;		/* continuous variable, log(rss) */

	/* branching info */
	int	ublb;
	int	*Branchz;		/* [3*p] */

	/* set solution */ 
	SCIP_Real *ep;
	
	int	nsols;
	int	store;
	SCIP_SOL**	sols;
	SCIP_Real	objval;

	SCIP_SOL*	sol;
	SCIP_Real*	solvals;
	SCIP_Bool	success;
	int			nvars	=	SCIPgetNVars(scip);
	SCIP_VAR**	vars;

	int	i,j,k,ct;
	int	memo;
	int	checkld;

#if MYPARA_LOG
	printf("forward selection");
	Longline();
#endif

   /* get heuristic data */
	/*
   SCIP_HEURDATA* heurdata;
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(lastsolindices != NULL);
	*/

	/* get values from probdata */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

	n	=	SCIPprobdataGetNdatas(probdata);
	p	=	SCIPprobdataGetNexvars(probdata);
	ndep	=	SCIPprobdataGetNdep(probdata);

	y	=	SCIPprobdataGety(probdata);
	orig_X_	=	SCIPprobdataGetX(probdata);
	orig_Q_	=	SCIPprobdataGetQ(probdata);
	orig_q	=	SCIPprobdataGetq(probdata);
	r	=	SCIPprobdataGetr(probdata);

	if( ndep ){
		Mdep		=	SCIPprobdataGetMdep(probdata);
		groupX	=	SCIPprobdataGetgroupX(probdata);
	}else{
		Mdep		=	NULL;
		groupX	=	NULL;
	}

	/* variables */
	var_a		=	SCIPprobdataGetVars_a(probdata);
	var_z		=	SCIPprobdataGetVars_z(probdata);
	var_ep	=	SCIPprobdataGetVars_ep(probdata);
	var_rss	=	SCIPprobdataGetVar_rss(probdata);
	var_log	=	SCIPprobdataGetVar_log(probdata);

	/* get branching info */
	/* alloc */
	SCIP_CALL( SCIPallocBufferArray(scip, &Branchz, 3*p));
	
	GenerateZeroVecInt( 3*p, Branchz);

	for(i=0; i<p; ++i){
		ublb					=	SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i]) 
								+	SCIPcomputeVarLbLocal(scip, var_z[i]));
		*(Branchz+(ublb*p)+i) 	= 	1;
	}

#if MYPARA_LOG
	for(i=0; i<3; i++){
		for(j=0; j<p; j++){
			printf("%d, ", *(Branchz+(i*p)+j));
		}
		newline();
	}
#endif


#if MYPARA_LOG
	printf("step1:\n");
#endif

	/* alloc */
	SCIP_CALL( SCIPallocBufferArray(scip, &list, p));
	SCIP_CALL( SCIPallocBufferArray(scip, &X_, n*1));
	SCIP_CALL( SCIPallocBufferArray(scip, &Y_, 1));
	SCIP_CALL( SCIPallocBufferArray(scip, &Xy, 1));
	SCIP_CALL( SCIPallocBufferArray(scip, &a_old, 1));

	GenerateZeroVecInt( p, list);
	dim = 0;
	for(i=0; i<p; i++){
		if( Branchz[i]==1 ){ /* if z_i is fixed to 0 */
			list[i] = -1;
		}else if( Branchz[p+i]==1 ){/* if z_i is unfixed */
			list[i] = 0;
		}else if ( Branchz[(2*p)+i]==1 ){
			dim++;
			list[i] = dim;

			if( dim==1 ){
				/* if dim == 1 */
				a_old[0] = orig_q[i] / mat_( orig_Q_, p, i, i);
				RSS = RSSvalue( 1, a_old, &orig_q[i], r);
				AIC = AICvalue( n, dim, RSS);

				/* update X_ and Xy */
				mydcopy_( &orig_X_[n * i], &X_[0], n);
				Xy[0] = orig_q[i];

				/* generate Y */ 
				Y_[0] = 1 / mat_( orig_Q_, p, i, i);

				/* if dim == 1 */
			}else{
				/* realloc */
				SCIP_CALL( SCIPreallocBufferArray( scip, &a_old, dim));
				SCIP_CALL( SCIPreallocBufferArray( scip, &Y_, dim*dim));	

				/* define Xy_new */
				SCIP_CALL( SCIPallocBufferArray( scip, &Xyn, dim));	
				mydcopy_( Xy, Xyn, dim-1);
				Xyn[dim-1] = orig_q[i];

				/* calculate a_old and update Y_ */
				SCIP_CALL( Calculate_a( scip, n, dim, X_, &orig_X_[n*i], Y_, Xyn, a_old, Y_));
				
				RSS = RSSvalue( dim, a_old, Xyn, r);
				AIC = AICvalue( n, dim, RSS);
				
				SCIPfreeBufferArray( scip, &Xyn);

				/* update X_ and Xy */
				SCIP_CALL( SCIPreallocBufferArray( scip, &X_, n*dim));	
				SCIP_CALL( SCIPreallocBufferArray( scip, &Xy, dim));	
				mydcopy_( &orig_X_[n*i], &X_[n*(dim-1)], n);
				Xy[dim-1] = orig_q[i];

			}

#if MYPARA_LOG
			printf("---> %dth variable, AIC:%f\n", i, AIC);
#endif
		}else{
			printf("error:heur_forward.c\n");
			stop();
		}
	}

#if MYPARA_LOG
	printf("list:");
	for(i=0; i<p; i++){
		if( i%10 == 0 ){
			newline();
		}
		printf("%d, ", list[i]);
	}
	newline();
	printf("step2:\n");
#endif

	if( dim==0 ){
		memo = -1;
		RSS = 1e+06;
		dim++;
		for(i=0; i<p; i++){
			if( list[i]==0 ){
				a_old[0] = orig_q[i] / mat_( orig_Q_, p, i, i);
				RSS_new = RSSvalue( 1, a_old, &orig_q[i], r);
				if( RSS_new < RSS ){
					RSS = RSS_new;
					memo = i;
				}
			}
#if MYPARA_LOG
			printf("%d: RSS = %f\n", i, RSS_new);
#endif
		}

		if( memo < 0 || memo >= p){
			printf("error in heur_forward.c\n");
			stop();
		}

		AIC = AICvalue( n, dim, RSS);
		list[memo] = dim;

		/* update X_ and Xy */
		mydcopy_( &orig_X_[n * memo], &X_[0], n);
		Xy[0] = orig_q[memo];

		/* generate Y */ 
		Y_[0] = 1 / mat_( orig_Q_, p, memo, memo);

#if MYPARA_LOG
	printf("---> %dth variable, AIC:%f\n", memo, AIC);
#endif
	}

#if MYPARA_LOG
	printf("list:");
	for(i=0; i<p; i++){
		if( i%10 == 0 ){
			newline();
		}
		printf("%d, ", list[i]);
	}
	newline();
#endif
	
	if( AIC >= 1e+06 ){
		printf("error:heur_forward.c\n");
		stop();

		/* free */
		SCIPfreeBufferArray(scip, &Branchz);
		SCIPfreeBufferArray(scip, &list);
		SCIPfreeBufferArray(scip, &X_);
		SCIPfreeBufferArray(scip, &Y_);
		SCIPfreeBufferArray(scip, &Xy);
		SCIPfreeBufferArray(scip, &a_old);
	
		*result = SCIP_DIDNOTFIND;
		return SCIP_OKAY;
	}

	while(1){
		dim++;
		memo = -1;
		RSS = 1e+06;

#if MYPARA_LOG
		printf("(dim=%d) ", dim);
		Longline();
#endif

		/* alloc */
		SCIP_CALL( SCIPallocBufferArray( scip, &a_new, dim));
		SCIP_CALL( SCIPallocBufferArray( scip, &a, dim));
		SCIP_CALL( SCIPallocBufferArray( scip, &W_, (dim)*(dim)));
		SCIP_CALL( SCIPallocBufferArray( scip, &Wnew_, (dim)*(dim)));
		SCIP_CALL( SCIPallocBufferArray( scip, &Xyn, dim));	
		mydcopy_( Xy, Xyn, dim-1);
		
		for(i=0; i<p; i++){

			checkld = 1;
			if( ndep ){
				for(j=0; j<ndep; j++){
					for(k=0; k<p; k++){
						if( groupX[(j*p)+k]==1 ){
							if( (i!=k) && (list[k]<=0) ) break;
							if( k==Mdep[j] ) checkld = 0;
						}
					}
					if( checkld==0 ) break;
				}
			}

			if( list[i]==0 && checkld ){
				/* define Xy_new */
				Xyn[dim-1] = orig_q[i];

				/* calculate a_old and store W_ for updating Y_ */
				SCIP_CALL( Calculate_a( scip, n, dim, X_, &orig_X_[n*i], Y_, Xyn, a_new, Wnew_));

				/* test */
				RSS_new = RSSvalue( dim, a_new, Xyn, r);
				if( RSS_new < RSS ){
					RSS = RSS_new;
					memo = i;
					mydcopy_( Wnew_, W_, dim*dim);
					mydcopy_( a_new, a, dim);
				}
#if MYPARA_LOG
				printf("%d: RSS = %f\n", i, RSS_new);
#endif
			}
		}

#if 0
		if( memo < 0 || memo >= p ){
			if( memo == -1 ){
				for(i=0; i<p; i++){
					if( list[i] == 0 ){
						memo = i;
						break;
					}
				}
				if( memo != -1 ){
					printf("error in heur_forward.c[memo==%d]\n", memo);
					stop();
				}
			}else{
				printf("error in heur_forward.c\n");
				stop();
			}
		}
#endif

		AIC_new = AICvalue( n, dim, RSS);
		if( AIC_new < AIC ){
			AIC = AIC_new;
			list[memo] = dim;

#if MYPARA_LOG
			printf("---> %dth variable, AIC:%f\n", memo, AIC);
#endif
			/* update */
			SCIP_CALL( SCIPreallocBufferArray(scip, &Y_, dim*dim));
			mydcopy_( W_, Y_, dim*dim);

			SCIP_CALL( SCIPreallocBufferArray(scip, &a_old, dim));
			mydcopy_( a, a_old, dim);

			/* update X_ and Xy */
			SCIP_CALL( SCIPreallocBufferArray( scip, &X_, n*dim));	
			SCIP_CALL( SCIPreallocBufferArray( scip, &Xy, dim));	
			mydcopy_( &orig_X_[n*memo], &X_[n*(dim-1)], n);
			Xy[dim-1] = orig_q[memo];

		}else{
			memo = -1;
#if MYPARA_LOG
			printf("--> no selection, (AIC:%f)\n", AIC_new);
#endif
		}

		/* free */
		SCIPfreeBufferArray(scip, &a_new);
		SCIPfreeBufferArray(scip, &a);
		SCIPfreeBufferArray(scip, &W_);
		SCIPfreeBufferArray(scip, &Wnew_);
		SCIPfreeBufferArray( scip, &Xyn);

		if( memo == -1 ){
			dim--;
			break;
		}

	}

	/* check object value of solution */
	nsols = SCIPgetNSols(scip);
	
	if( nsols < MP_NUM_SOL ){
		store = 1;
	}else{
		sols = SCIPgetSols(scip);
		objval = AIC;
		nsols = MP_NUM_SOL;

		if( objval < SCIPgetSolOrigObj(scip,sols[nsols-1]) ){
			store = 1;
		}else{
			store = 0;
		}
	}

	if( store ){
		/*  generate solution  */
		/* alloc */
		SCIP_CALL( SCIPallocBufferArray(scip, &ep, n));
		dgemv_1( X_, n, dim, a_old, y, -1.0, 1.0, ep);
		
	
		/* set solution */
		/* alloc */
		SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars));
		SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars));
	
		ct=0;
		
		/* a */
		for(i=0; i<p; ++i){
			vars[ct] = var_a[i];
			if( list[i] > 0 ){
				solvals[ct] = a_old[list[i]-1];
			}else{
				solvals[ct] = 0.0;
			}
			ct++;
		}
	
		/* z */
		for(i=0; i<p; i++){
			vars[ct] = var_z[i];
			if( list[i] > 0 ){
				solvals[ct] = 1.0;
			}else{
				solvals[ct] = 0.0;
			}
			ct++;
		}
	
		/* ep */
		for(i=0; i<n; ++i){
			vars[ct]		=	var_ep[i];
			solvals[ct]	=	ep[i];
			ct++;
		}
	
		vars[ct]		=	var_rss;
		solvals[ct] =	myddot_( ep, ep, n);
		ct++;
	
		vars[ct]		=	var_log;
		solvals[ct]	=	log(myddot_( ep, ep, n));
		ct++;
	
		if( ct!=nvars ){
			SCIPerrorMessage("It is unexpected error in set sol,");
			printf("( ct, nvars) = ( %d, %d)", ct, nvars);
			stop();
		}
	
		SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
		SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));
		SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, FALSE, TRUE, TRUE, &success));

		/* free */
		SCIPfreeBufferArray(scip, &ep);
		SCIPfreeBufferArray(scip, &solvals);
		SCIPfreeBufferArray(scip, &vars);
	}

	/* free */
	SCIPfreeBufferArray(scip, &Branchz);
	SCIPfreeBufferArray(scip, &list);
	SCIPfreeBufferArray(scip, &X_);
	SCIPfreeBufferArray(scip, &Y_);
	SCIPfreeBufferArray(scip, &Xy);
	SCIPfreeBufferArray(scip, &a_old);
	
	*result = SCIP_FOUNDSOL;
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */


/** creates the local primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurForward(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Local primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

	/*
   heurdata->
	*/

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecForward, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyForward) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeForward) );

   return SCIP_OKAY;
}

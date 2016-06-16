
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
#include "cblapack.h"

#define HEUR_NAME             "forward"
#define HEUR_DESC             "primal heuristic using forward selection"
#define HEUR_DISPCHAR         'f'
#define HEUR_PRIORITY         1000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH        	0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE

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
	SCIP_Real*	orig_X_;			/* [n*p] */
	SCIP_Real*	orig_Q_;		/* [p*p] <- (X^t) X */
	SCIP_Real*	orig_q;		/* [p]   <- (X^t) y */
	SCIP_Real	r;

	int*	Mdep;					/* [ndep] */
	//int*	groupX;				/* [ndep*p] */

	/* for forward selection */
	int	dim;
	int*	list;					/* [p] */
	SCIP_Real*	a;				/* [dim] */
	SCIP_Real*	a_old;		/* [dim-1] */
	SCIP_Real*	a_new;		/* [dim] */
	SCIP_Real	RSS;			/* residual sum of square */
	SCIP_Real	RSS_new;
	SCIP_Real	AIC;
	SCIP_Real	AIC_new;

	/*
	 *	X: sub matrix of orig_X_ 
	 *	Y:	(X^t X)^-1 
	 * X_new = (X, x_i);
	 * Z: (X_new ^t X_new)^-1
	 *		= ( V   v
	 			 v^t u )
	 */

	SCIP_Real*	Xy;	/* sub vector of orig_q */ 
	SCIP_Real*	X_;	
	SCIP_Real*	Y_;	/* [(dim-1)*(dim-1)] */
	SCIP_Real*	Z_;	/* [dim*dim] */
	SCIP_Real*	W_;	/* [dim*dim] */
	SCIP_Real*	V_;	/* [(dim-1)*(dim-1)] */
	SCIP_Real*	v;		/* [dim-1] */
	SCIP_Real	u;
	
	SCIP_Real*	b;		/* [dim-1] */
	SCIP_Real*	c;		/* [dim-1] */
	SCIP_Real*	d;		/* [n] */

	/* variables */
	SCIP_VAR**	var_a;		/* [p] continuous variables */
	SCIP_VAR**	var_z;		/* [p] 01 variables */
	SCIP_VAR**	var_ep;		/* [n] continuous variables */
	SCIP_VAR*	var_rss;		/* continuous variable, residual sum of squares */
	SCIP_VAR*	var_log;		/* continuous variable, log(rss) */

	/* set solution */ 
	SCIP_Real *ep;
	
	SCIP_SOL*	sol;
	SCIP_Real*	solvals;
	SCIP_Bool	success;
	int			nvars	=	SCIPgetNVars(scip);
	SCIP_VAR**	vars;

	int 	i,j,t,ct;
	int	memo;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);

#if MYPARA_LOG
	printf("forward selection!\n");
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
		//groupX	=	SCIPprobdataGetgroupX(probdata);
	}else{
		Mdep		=	NULL;
		//groupX	=	NULL;
	}

	/* dim = 1 */
	dim = 1;
	
	/* alloc */
	SCIP_CALL( SCIPallocBufferArray(scip, &a_old, dim));
	SCIP_CALL( SCIPallocBufferArray(scip, &X_, n*p));
	SCIP_CALL( SCIPallocBufferArray(scip, &Xy, p));
	SCIP_CALL( SCIPallocBufferArray(scip, &d, n));
	SCIP_CALL( SCIPallocBufferArray(scip, &list, p));
	GenerateZeroVecInt( p, list);

	if( ndep ){
		for(i=0; i<ndep; i++){
			list[Mdep[i]] = -1;
		}
	}

#if MYPARA_LOG
	printf("(dim=%d) ", dim);
	Longline();
#endif

	RSS = 1e+06;
	memo = -1;

	for(i=0; i<p; i++){
		if( list[i] == 0 ){
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
	
	if( memo < 0 || memo >= p ){
		printf("error in heur_forward.c\n");
		stop();
	}

	AIC = AICvalue( n, dim, RSS);
	list[memo] = dim;
	
	/* update X_ and Xy */
	mydcopy_( &orig_X_[n * memo], &X_[n * (dim-1)], n);
	Xy[dim-1] = orig_q[memo];

	/* generate Y ( dim = 1 ) */
	SCIP_CALL( SCIPallocBufferArray( scip, &Y_, dim*dim));
	Y_[0] = 1 / mat_( orig_Q_, p, memo, memo);

#if MYPARA_LOG
	printf("---> %dth variable, AIC:%f\n", memo, AIC);
#endif


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
		SCIP_CALL( SCIPallocBufferArray( scip, &b, dim-1));
		SCIP_CALL( SCIPallocBufferArray( scip, &c, dim-1));
		SCIP_CALL( SCIPallocBufferArray( scip, &v, dim-1));
		SCIP_CALL( SCIPallocBufferArray( scip, &V_, (dim-1)*(dim-1)));
		SCIP_CALL( SCIPallocBufferArray( scip, &Z_, (dim)*(dim)));
		SCIP_CALL( SCIPallocBufferArray( scip, &W_, (dim)*(dim)));
		
		for(i=0; i<p; i++){
			/*
			 * 1. b <- X^t x_i
			 * 2.	c <- Y b
			 * 3. d <- - X c + x_i
			 * 4. u <- 1 / <x_i, d> 
			 * 5. v <- - u c
			 * 6. V <- Y + u c c^t
			 * 7. Z <- ( V    v
			 	          v^t  u )
			 * 8. a_new <- Z (Xy)
			 */

			if( list[i]==0 ){

				/* 1. b <- X^t x_i */
				dgemv_t( X_, n, dim-1, &orig_X_[n * i], b);
				//printv( dim-1, b);

				/* 2. c <- Y b */
				dgemv_2( Y_, dim-1, dim-1, b, c);
				//printv( dim-1, c);

				/* 3. d <- - X c + x_i */
				dgemv_1( X_, n, dim-1, c, &orig_X_[n * i], -1.0, 1.0, d);
				//printv( n, d);

				/* 4. u <- 1/<x_i, d> */
				u = 1.0 / myddot_( &orig_X_[n * i], d, n);
				//prints(u);
				
				/* 5. v <- - u c */
				mydscal_( c, dim-1, -u, v);
				//printv( dim-1, v);

				/* 6. V <- Y + u c c^t */
				dger_1( Y_, c, c, dim-1, dim-1, u, V_);
				//printM_( V_, dim-1, dim-1);

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
				//printM_( Z_, dim, dim);

				/* 8. a_new <- Z (Xy) */
				Xy[dim-1] = orig_q[i];
				dgemv_2( Z_, dim, dim, Xy, a_new);
				//printv( dim, a_new);

				/* test */
				RSS_new = RSSvalue( dim, a_new, Xy, r);
				if( RSS_new < RSS ){
					RSS = RSS_new;
					memo = i;
					mydcopy_( Z_, W_, dim*dim);
					mydcopy_( a_new, a, dim);
				}

#if MYPARA_LOG
				printf("%d: RSS = %f\n", i, RSS_new);
#endif

			}
		}

		if( memo < 0 || memo >= p ){
			printf("error in heur_forward.c\n");
			stop();
		}

		AIC_new = AICvalue( n, dim, RSS);
		if( AIC_new < AIC ){
			AIC = AIC_new;
			list[memo] = dim;

#if MYPARA_LOG
			printf("---> %dth variable, AIC:%f\n", memo, AIC);
#endif

			/* copy and free */
			SCIPfreeBufferArray(scip, &Y_);
			SCIP_CALL( SCIPallocBufferArray(scip, &Y_, dim*dim));
			mydcopy_( W_, Y_, dim*dim);

			SCIPfreeBufferArray(scip, &a_old);
			SCIP_CALL( SCIPallocBufferArray(scip, &a_old, dim));
			mydcopy_( a, a_old, dim);

			/* update X_ and Xy */
			mydcopy_( &orig_X_[n * memo], &X_[n * (dim-1)], n);
			Xy[dim-1] = orig_q[memo];

		}else{
			memo = -1;
			SCIPfreeBufferArray(scip, Y_);
#if MYPARA_LOG
			printf("--> no selection, (AIC:%f)\n", AIC_new);
#endif
		}

		/* free */
		SCIPfreeBufferArray(scip, &a_new);
		SCIPfreeBufferArray(scip, &a);
		SCIPfreeBufferArray(scip, &b);
		SCIPfreeBufferArray(scip, &c);
		SCIPfreeBufferArray(scip, &v);
		SCIPfreeBufferArray(scip, &V_);
		SCIPfreeBufferArray(scip, &Z_);
		SCIPfreeBufferArray(scip, &W_);

		if( memo == -1 ){
			dim--;
			break;
		}
	}

	/* variables */
	var_a		=	SCIPprobdataGetVars_a(probdata);
	var_z		=	SCIPprobdataGetVars_z(probdata);
	var_ep	=	SCIPprobdataGetVars_ep(probdata);
	var_rss	=	SCIPprobdataGetVar_rss(probdata);
	var_log	=	SCIPprobdataGetVar_log(probdata);

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
	SCIPfreeBufferArray(scip, &a_old);
	SCIPfreeBufferArray(scip, &list);
	SCIPfreeBufferArray(scip, &d);
	SCIPfreeBufferArray(scip, &X_);
	SCIPfreeBufferArray(scip, &Xy);
	SCIPfreeBufferArray(scip, &ep);
	SCIPfreeBufferArray(scip, &solvals);
	SCIPfreeBufferArray(scip, &vars);

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

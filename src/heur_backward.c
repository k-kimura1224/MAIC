
/**@file   heur_backward.c
 * @brief  Backward selection
 *
 * 
 *
 */
#include <assert.h>
#include <string.h>
#include <math.h>

#include "heur_backward.h"
#include "probdata_linereg.h"
#include "matrix.h"
#include "vector.h"
#include "set_myparameter.h"
#include "cblapack.h"

#define HEUR_NAME             "backward"
#define HEUR_DESC             "primal heuristic using backward selection"
#define HEUR_DISPCHAR         'b'
#define HEUR_PRIORITY         1000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         10
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
/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyBackward)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurBackward(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeBackward)
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
SCIP_DECL_HEUREXEC(heurExecBackward)
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
	int*	groupX;				/* [ndep*p] */

	/* for backward selection */
	int	dim;
	int*	list;					/* [p] */
	SCIP_Real	RSS;			/* residual sum of square */
	SCIP_Real	RSS_new;
	SCIP_Real	AIC;
	SCIP_Real	AIC_new;
	SCIP_Real*	a_new;		/* [dim] */
	SCIP_Real*	a;				/* [dim] */
	SCIP_Real*	a_old;		/* [dim] */

	int	ublb;
	int	*Branchz;		/* [3*p] */

	SCIP_Real*	Q_;	/* sub matrix of orig_Q_ */
	SCIP_Real*	Xy;	/* sub vector of orig_q */ 

	/* variables */
	SCIP_VAR**	var_a;		/* [p] continuous variables */
	SCIP_VAR**	var_z;		/* [p] 01 variables */
	SCIP_VAR**	var_ep;		/* [n] continuous variables */
	SCIP_VAR*	var_rss;		/* continuous variable, residual sum of squares */
	SCIP_VAR*	var_log;		/* continuous variable, log(rss) */

	/* set solution */ 
	SCIP_Real *ep;
	SCIP_Real*	X_;	
	
	SCIP_SOL*	sol;
	SCIP_Real*	solvals;
	SCIP_Bool	success;
	int			nvars	=	SCIPgetNVars(scip);
	SCIP_VAR**	vars;

	int	nsols;
	SCIP_SOL**	sols;
	int	store;
	SCIP_Real	objval;

	int 	i,j,t,ct,ct_a;
	int	memo;
	int	dpv;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);

#if MYPARA_LOG
	printf("backward selection!\n");
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

	if( ndep ){
		for(i=0; i<ndep; i++){
			memo = -1; 
			for(j=0; j<p; j++){
				if( *(groupX+(i*p)+j)==1 ){
					if( *(Branchz+j)==1 ) break;
					if( *(Branchz+p+j)==1 ) memo=j;
					if( j==Mdep[i] ){
						if( memo==-1 ){
							printf("error in heur_backward.c\n");
							stop();
						}
						*(Branchz+p+memo) = 0;
						*(Branchz+memo) = 1;
						break;
					}
				}
			}
		}
	}
	
#if MYPARA_LOG
	printf("linear dependent\n");
	if( ndep ){
		for(i=0; i<3; i++){
			for(j=0; j<p; j++){
				printf("%d, ", *(Branchz+(i*p)+j));
			}
			newline();
		}
	}
#endif


	/* alloc */
	SCIP_CALL( SCIPallocBufferArray(scip, &list, p));

	/* list */
	for(i=0; i<p; i++){
		list[i] = 1 - Branchz[i];
	}

	dim = p - sumint( &Branchz[0], p);
	AIC = 1e+06;
	SCIP_CALL( SCIPallocBufferArray(scip, &a_old, dim));

	while(1){
		dim--;
		memo = -1;
		RSS = 1e+06;
#if MYPARA_LOG
		printf("(dim=%d) ", dim);
		Longline();
#endif

		/* alloc */
		SCIP_CALL( SCIPallocBufferArray( scip, &a_new, dim));
		SCIP_CALL( SCIPallocBufferArray( scip, &Q_, dim*dim));
		SCIP_CALL( SCIPallocBufferArray( scip, &Xy, dim));
		SCIP_CALL( SCIPallocBufferArray( scip, &a, dim));
		
		for(i=0; i<p; i++){
			/*
			 * solve
			 *		Q a_new = Xy
			 *
			 */

			if( (*(Branchz+p+i)==1) && (list[i]==1) ){

				/* generate Q and Xy */
				/* Q */
				ct=0;
				for(j=0; j<p; j++){
					if( (list[j]==1) && (j != i) ){
						for(t=0; t<p; t++){
							if( (list[t]==1) && (t != i)){
								Q_[ct++] = mat_( orig_Q_, p, j, t);
							}
						}
					}
				}

				if( ct != (dim*dim) ){
					printf("error in heur_backward.c\n");
					stop();
				}
				
				/* Xy */
				ct = 0;
				for(j=0; j<p; j++){
					if( (list[j]==1) && (j != i) ){
						Xy[ct++] = orig_q[j];
					}
				}

				if( ct != dim ){
					printf("error in heur_backward.c\n");
					stop();
				}

				dpv = _dposv_( Q_, Xy, dim, a_new );

				if( dpv == 0 ){

					/* test */
					RSS_new = RSSvalue( dim, a_new, Xy, r);
					if( RSS_new < RSS ){
						RSS = RSS_new;
						memo = i;
						mydcopy_( a_new, a, dim);
					}
#if MYPARA_LOG
					printf("%d: RSS = %f\n", i, RSS_new);
#endif
				}
			}
		}

		if( memo < 0 || memo >= p ){
			/* free */
			SCIPfreeBufferArray(scip, &a_new);
			SCIPfreeBufferArray(scip, &Q_);
			SCIPfreeBufferArray(scip, &Xy);
			SCIPfreeBufferArray(scip, &a);
			break;
		}

		AIC_new = AICvalue( n, dim, RSS);
		if( AIC_new < AIC ){
			AIC = AIC_new;
			list[memo] = 0;

#if MYPARA_LOG
			printf("---> %dth variable, AIC:%f\n", memo, AIC);
#endif

			/* copy */
			SCIPfreeBufferArray(scip, &a_old);
			SCIP_CALL(SCIPallocBufferArray(scip, &a_old, dim));
			mydcopy_( a, a_old, dim);

		}else{
			memo = -1;
#if MYPARA_LOG
			printf("--> no selection, (AIC:%f)\n", AIC_new);
#endif
		}

		/* free */
		SCIPfreeBufferArray(scip, &a_new);
		SCIPfreeBufferArray(scip, &Q_);
		SCIPfreeBufferArray(scip, &Xy);
		SCIPfreeBufferArray(scip, &a);

		if( memo == -1 ){
			dim++;
			break;
		}else if( sumint(list,p) == sumint(&Branchz[2*p],p) ){
			break;
		}
	}

	if( AIC >= 1e+06 ){
		/* free */
		SCIPfreeBufferArray(scip, &Branchz);
		SCIPfreeBufferArray(scip, &list);
		SCIPfreeBufferArray(scip, &a_old);
	   return SCIP_OKAY;
	}

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
		SCIP_CALL( SCIPallocBufferArray(scip, &X_, n*dim));
	
		/* ep */
		pickupcolumn_( orig_X_, n, p, list, X_);
		dgemv_1( X_, n, dim, a_old, y, -1.0, 1.0, ep);
	
		/* set solution */
		/* alloc */
		SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars));
		SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars));
	
		ct=0;
		ct_a=0;	
		/* a */
		for(i=0; i<p; ++i){
			vars[ct] = var_a[i];
			if( list[i] == 1 ){
				solvals[ct] = a_old[ct_a++];
			}else{
				solvals[ct] = 0.0;
			}
			ct++;
		}
	
		if( ct_a != dim ){
			printf("error in heur_backward.c\n");
			stop();
		}
	
		/* z */
		for(i=0; i<p; i++){
			vars[ct] = var_z[i];
			solvals[ct] = (double)list[i];
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
		SCIPfreeBufferArray(scip, &X_);
		SCIPfreeBufferArray(scip, &solvals);
		SCIPfreeBufferArray(scip, &vars);
	}

	/* free */
	SCIPfreeBufferArray(scip, &Branchz);
	SCIPfreeBufferArray(scip, &list);
	SCIPfreeBufferArray(scip, &a_old);
	
	*result = SCIP_FOUNDSOL;
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */


/** creates the local primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurBackward(
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
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecBackward, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyBackward) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeBackward) );

   return SCIP_OKAY;
}

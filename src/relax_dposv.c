
/**@file   relax_dposv.c
 * @brief  
 */

#include <math.h>
#include <assert.h>

#include "relax_dposv.h"
#include "probdata_linereg.h"
#include "matrix.h"
#include "vector.h"
#include "get_mybranchvarz.h"
#include "cblapack.h"
#include "set_myparameter.h"
#include "linear_dependent.h"

#define RELAX_NAME             "myrelaxator_dposv"
#define RELAX_DESC             "relaxator dposv"
#define RELAX_PRIORITY         1000
#define RELAX_FREQ             1

#define MYPARA_LOG				0
/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
	int test;
};

/*
 * Local methods
 */


/*
 * Callback methods of relaxator
 */

/** copy method for relaxator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_RELAXCOPY(relaxCopyDposv)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(relax != NULL);
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);

   /* call inclusion method of relaxator */
	SCIP_CALL( SCIPincludeRelaxDposv(scip));

   return SCIP_OKAY;
}

/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeDposv)
{   /*lint --e{715}*/
	SCIP_RELAXDATA* relaxdata;

   assert(relax != NULL);
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);
   assert(scip != NULL);

   /* free relaxator data */
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   SCIPfreeMemory(scip, &relaxdata);
   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}

/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecDposv)
{

   SCIP_PROBDATA* probdata;
	int	n;
	int	p;
	int	ndep;

	SCIP_Real*	y;				/* [n] */
	SCIP_Real*	X;				/* [n*p] */
	SCIP_Real*	orig_Q;		/* [p*p] */
	SCIP_Real*	orig_q;		/* [p] */
	SCIP_Real	r;
	
	int*	Mdep;					/* [ndep] */
	int*	groupX;				/* [ndep*p] */
	
	/* variables */
	SCIP_VAR**	var_a;		/* [p] continuous variables */
	SCIP_VAR**	var_z;		/* [p] 01 variables */
	SCIP_VAR**	var_ep;		/* [n] continuous variables */
	SCIP_VAR*	var_rss;		/* continuous variable, residual sum of squares */
	SCIP_VAR*	var_log;		/* continuous variable, log(rss) */
	
	int	ublb;
	int	*Branchz;		/* [3*p] */
	
	SCIP_NODE* node;
	int	parentbr;			/* branching variable number */
	SCIP_Real	ParentDual;

	int	sum_Branchz[3];
	int	dim;					/* the number of z_j that we can use */
	int 	*z01;					/* size:p z01[j]=1 if z_j that we can use */
	int	info;					/* if successful, info = 0 */

	/* for while */
	SCIP_Real	*Q;			/* Q[dim*dim] is sub matrix of Q[p*p] */
	SCIP_Real	*q;			/* q[dim]     is sub vector of q[p]   */
	SCIP_Real	*a_z01;		/* a[dim] */
	SCIP_Real	*sub_x;
	int	*sub_D;
	int	sub_ndep;
	int	*sub_Mdep;
	int	*group;
	int	check;
	
	/* update local dual bound */
	SCIP_Real	rss;  		/* residual sum of square */
	SCIP_Real	log_rss;		/* log ( rss ) */
	SCIP_Real	Dual;			/* local dual bound */

	/* set the primal solution */ 
	SCIP_Real *a;
	SCIP_Real *ep;
	SCIP_Real *x;
	
	SCIP_SOL*	sol;
	SCIP_Real*	solvals;
	SCIP_HEUR*	heur;
	SCIP_Bool	success;
	int			nvars	=	SCIPgetNVars(scip);
	SCIP_VAR**	vars;
	
	int	nsols;
	SCIP_SOL**	sols;
	int	store;
	SCIP_Real	primal;

#if 0
	/* set the solution of relaxation */
	SCIP_Real*	relaxsolvals;
	SCIP_Real	bigM = 100.0;
#endif

	int	i,j,ct;

   assert(relax != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);
   assert(result != NULL);


   /* get relaxator data */
	/*
	SCIP_RELAXDATA* relaxdata;
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
	*/

	/* get from probdata */	
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

	n	=	SCIPprobdataGetNdatas(probdata);
	p	=	SCIPprobdataGetNexvars(probdata);
	ndep	=	SCIPprobdataGetNdep(probdata);

	y	=	SCIPprobdataGety(probdata);
	X	=	SCIPprobdataGetX(probdata);
	orig_Q	=	SCIPprobdataGetQ(probdata);
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
	SCIP_CALL( SCIPallocBufferArray(scip, &Branchz, 3*p));
	
	GenerateZeroVecInt( 3*p, Branchz);

	for(i=0; i<p; ++i){
		ublb					=	SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i]) 
								+	SCIPcomputeVarLbLocal(scip, var_z[i]));
		*(Branchz+(ublb*p)+i) 	= 	1;
	}
	
	/* use parent node info */
	node = SCIPgetCurrentNode(scip);

#if MYPARA_LOG
	printf("(%lld)", node->number);
	Longline();
	for(i=0; i<3; i++){
		for(j=0; j<p; j++){
			printf("%d, ", *(Branchz+(i*p)+j));
		}
		newline();
	}
#endif

	parentbr = -1;

	if( sumint(&Branchz[p],p)<p){ 
		parentbr = SCIPparentGetBranchingZ(scip, node, p, var_z);
		if( parentbr == -1 ){
			
#if MYPARA_LOG || 0
				printf("branching var:?????\n");
#endif

#if 0
			if( sumint(&Branchz[p],p)==0 ){
				
				/* free */
				SCIPfreeBufferArray(scip, &Branchz);

				*result=SCIP_CUTOFF;
				return SCIP_OKAY;
			}else{
#if MYPARA_LOG
				printf("branching var:?????\n");
#endif
				
				/* free */
				SCIPfreeBufferArray(scip, &Branchz);

				*result=SCIP_SUCCESS;
				return SCIP_OKAY;
				/*stop();*/
			 }
#endif
		}
	}

#if MYPARA_LOG 
	if( parentbr>=0 ){
		printf("branching var:z_");
		printN(parentbr+1);
	}else{
		printf("\n");
	}
#endif

	if( parentbr>=0 ){
		if( *(Branchz+(2*p)+parentbr)==1 ){

			SCIP_NODE*	parent;
			parent	=	SCIPnodeGetParent(node);
			ParentDual = SCIPgetNodeLowerbound(scip,parent);
			SCIP_CALL( SCIPupdateLocalLowerbound(scip,ParentDual+2.0));
			
			/* free */
			SCIPfreeBufferArray(scip, &Branchz);

			*result=SCIP_SUCCESS;
			return SCIP_OKAY;
		}
	}

	/* adjust fixed z from linear dependent */
	if( ndep ){
		int	dummy;
		for(i=0; i<ndep; ++i){
			dummy=-1;
			for(j=0; j<p; ++j){
				if( *(groupX+(i*p)+j)==1 ){
					if( *(Branchz+j)==1 ) break;
					if( *(Branchz+p+j)==1 ) dummy=j;
					if( j==*(Mdep+i) ){
						if( dummy==-1 ) stop();
						*(Branchz+p+dummy) = 0;
						*(Branchz+dummy) = 1;
						break;
					}
				}
			}
		}
	}
 
#if MYPARA_LOG 
	if( ndep ){
		for(i=0; i<3; i++){
			for(j=0; j<p; j++){
				printf(" %d,", *(Branchz+(i*p)+j));
			}
			newline();
		}
	}
#endif
		
	SCIP_CALL( SCIPallocBufferArray(scip, &z01, p));

	while(1){
		for(i=0; i<3; ++i) sum_Branchz[i] = sumint(&Branchz[i*p], p);
		dim = sum_Branchz[1] + sum_Branchz[2];

		if( dim==0 ){
			SCIPerrorMessage("It is an unexpected error1");
			for(i=0; i<3; i++){
				for(j=0; j<p; j++){
					printf(" %d,", *(Branchz+(i*p)+j));
				}
				newline();
			}
			stop();
		}

		/* alloc */
		SCIP_CALL( SCIPallocBufferArray(scip, &Q, dim*dim));
		SCIP_CALL( SCIPallocBufferArray(scip, &q, dim));
		SCIP_CALL( SCIPallocBufferArray(scip, &a_z01, dim));

		/**
		 * define
		 *		Q[dim*dim], q[dim]
		 *	with the branching info
		**/

		for(i=0; i<p; ++i) z01[i] = 1 - *(Branchz+i);
			
		ct=0;
		for(i=0; i<p; ++i){
			for(j=0; j<p; ++j){
				if( (z01[i]==1) && (z01[j]==1) ){
					Q[ct++] = orig_Q[(i*p)+j];
				}
			}
		}

		ct=0;
		for(i=0; i<p; ++i){
			if( z01[i]==1 ) q[ct++] = orig_q[i];
		}
		
		/* solve with dposv of clapack */
		info = _dposv_( Q, q, dim, a_z01);

		if( info==0 )	break;	/* solving is successful */ 

		SCIP_CALL( SCIPallocBufferArray(scip, &sub_x, n*dim));
		SCIP_CALL( SCIPallocBufferArray(scip, &sub_D, dim));

		pickupcolumn_(X, n, p, z01, sub_x);
		LinearDependent( n, dim, sub_x, sub_D);	
		sub_ndep = sumint( sub_D, dim);
			
		if( sub_ndep==0 ){
			SCIPerrorMessage("It is an unexpected error2");
			stop();
		}
	
		SCIP_CALL( SCIPallocBufferArray(scip, &sub_Mdep, sub_ndep));
		SCIP_CALL( SCIPallocBufferArray(scip, &group, sub_ndep*dim));
		
		ct=0;
		for(i=0; i<dim; ++i){
			if( sub_D[i]==1 ){
				sub_Mdep[ct++] = i;
			}
		}

		GenerateZeroVecInt( sub_ndep*dim, group);
		
		for(i=0; i<sub_ndep; ++i){
			*(group+(i*dim)+sub_Mdep[i]) = 1;
		}

		LinearDependentGroup( dim, Q, sub_ndep, sub_Mdep, sub_D, group);

		check = 0;

		for(i=0; i<sub_ndep; ++i){
			ct=0;
			for(j=0; j<p; ++j){
				if( z01[j]==1 ){
					if( *(group+(i*dim)+ct)==1 ){
						if( *(Branchz+j)==1 ) break;
						if( *(Branchz+p+j)==1 ){
							*(Branchz+p+j) = 0;
							*(Branchz+j) = 1;
							check = 1;
							break;
						}
					}
					ct++;
				}
			}
		}

		/* free */
		SCIPfreeBufferArray(scip, &sub_x);
		SCIPfreeBufferArray(scip, &sub_D);
		SCIPfreeBufferArray(scip, &sub_Mdep);
		SCIPfreeBufferArray(scip, &group);

		if(!check){
			/* free */
			SCIPfreeBufferArray(scip, &Branchz);

			SCIPfreeBufferArray(scip, &z01);

			SCIPfreeBufferArray(scip, &Mdep);
			SCIPfreeBufferArray(scip, &groupX);
			SCIPfreeBufferArray(scip, &Q);
			SCIPfreeBufferArray(scip, &q);
			SCIPfreeBufferArray(scip, &a_z01);
			*result=SCIP_CUTOFF;
			return SCIP_OKAY;
		}

		/* free */
		SCIPfreeBufferArray(scip, &Q);
		SCIPfreeBufferArray(scip, &q);
		SCIPfreeBufferArray(scip, &a_z01);

	}
			
	/* update local lower bound */
	rss		=	r - myddot_( q, a_z01, dim);
	log_rss	=	log(rss);
	Dual = (double)n * log_rss + 2.0*(double)sum_Branchz[2];

	SCIP_CALL( SCIPupdateLocalLowerbound(scip, Dual));


	if( rss<=0 ){
		prints(rss);
		SCIPerrorMessage("rss <= 0");
		stop();
	}

	nsols = SCIPgetNSols(scip);
	
	if( nsols < MP_NUM_SOL ){
		store = 1;
	}else{
		sols = SCIPgetSols(scip);
		primal = Dual + 2.0*(double)sum_Branchz[1];
		nsols = MP_NUM_SOL;
		
		if( primal < SCIPgetSolOrigObj(scip,sols[nsols-1]) ){
			store = 1;
		}else{
			store = 0;
		}
	}
	
	if( store ){
		/*  generate solution from branching info */
		/* alloc */
		SCIP_CALL( SCIPallocBufferArray(scip, &a, p));
		SCIP_CALL( SCIPallocBufferArray(scip, &ep, n));
		SCIP_CALL( SCIPallocBufferArray(scip, &x, n*dim));
	
		ct=0;
		for(i=0; i<p; ++i){
			if( *(Branchz+i)==0 )	a[i]	=	a_z01[ct++];
			else							a[i]	=	0.0;
		}
	
		if( ct!=dim ){
			printN(dim);
			SCIPerrorMessage("It is unexpected error2");
			stop();
		}
	
		pickupcolumn_(X, n, p, z01, x);
		dgemv_1( x, n, dim, a_z01, y, -1.0, 1.0, ep);
	
		/* set solution */
		SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars));
		SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars));
	
		ct=0;
	
		for(i=0; i<p; ++i){
			vars[ct]		=	var_a[i];
			solvals[ct]	=	a[i];
			ct++;
		}
	
		for(i=0; i<p; ++i){
			vars[ct]		=	var_z[i];
			solvals[ct]	=	(double)z01[i];
			ct++;
		}
		
		for(i=0; i<n; ++i){
			vars[ct]		=	var_ep[i];
			solvals[ct]	=	(double)ep[i];
			ct++;
		}
		
		vars[ct]		=	var_rss;
		solvals[ct]	=	rss;
		ct++;
	
		vars[ct]		=	var_log;
		solvals[ct]	=	log_rss;
		ct++;
	
		if( ct!=nvars ){
			SCIPerrorMessage("It is unexpected error in set sol,");
			printf("( ct, nvars) = ( %d, %d)", ct, nvars);
			stop();
		}
	
		heur	=	SCIPfindHeur(scip, "trysol");
		SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
		SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));
		SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, FALSE, TRUE, TRUE, &success));

		/* free */
		SCIPfreeBufferArray(scip, &a);
		SCIPfreeBufferArray(scip, &ep);
		SCIPfreeBufferArray(scip, &x);
		SCIPfreeBufferArray(scip, &solvals);
		SCIPfreeBufferArray(scip, &vars);
	}

#if 0
	/* set the solution of relaxation */
	SCIP_CALL( SCIPallocBufferArray(scip, &relaxsolvals, nvars));
	mydcopy_( solvals, relaxsolvals, nvars);

	for(i=0; i<p; ++i){
		ublb	=	SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i]) 
					+	SCIPcomputeVarLbLocal(scip, var_z[i]));

		if( ublb == 1 ){
			relaxsolvals[p+i] = fabs(relaxsolvals[i]) / bigM;
		}
	}
	SCIP_CALL( SCIPsetRelaxSolVals( scip, nvars, vars, relaxsolvals));
#endif

#if 0
	for(i=0; i<p; ++i){
		ublb	=	SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i]) 
					+	SCIPcomputeVarLbLocal(scip, var_z[i]));

		if( ublb == 1 ){
			SCIP_CALL( SCIPaddExternBranchCand( scip, var_z[i], 1000, 
			fabs(solvals[i]) / bigM));
		}
	}
#endif

	
	/* free */
	SCIPfreeBufferArray(scip, &Branchz);
	SCIPfreeBufferArray(scip, &z01);
	SCIPfreeBufferArray(scip, &Q);
	SCIPfreeBufferArray(scip, &q);
	SCIPfreeBufferArray(scip, &a_z01);
#if 0
	SCIPfreeBufferArray(scip, &relaxsolvals);
#endif
	*result=SCIP_SUCCESS;
   return SCIP_OKAY;
}

/*
 * relaxator  specific interface methods
 */


/** creates the myrelaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxDposv(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   /* create myrelaxator data */
   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );


   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecDposv, relaxdata) );

   assert(relax != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetRelaxCopy(scip, relax, relaxCopyDposv) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeDposv) );


   return SCIP_OKAY;
}

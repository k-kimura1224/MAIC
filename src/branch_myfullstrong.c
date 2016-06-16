/**@file   branch_myfullstrong.c
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "branch_myfullstrong.h"
#include "probdata_linereg.h"
#include "set_myparameter.h"
#include "vector.h"
#include "matrix.h"
#include "cblapack.h"

#define BRANCHRULE_NAME            "mfs"
#define BRANCHRULE_DESC            "mybranching rule "
#define BRANCHRULE_PRIORITY        200000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define MYPARA_LOG	0

/*
 * Data structures
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

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
	int a;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of branching rule
 */

/* TODO: Implement all necessary branching rule methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyMyfullstrong)
{  /*lint --e{715}*/

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleMyfullstrong(scip) ) ;

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeMyfullstrong)
{  /*lint --e{715}*/

   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitMyfullstrong)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 1
static
SCIP_DECL_BRANCHEXIT(branchExitMyfullstrong)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}
#else
#define branchExitMyfullstrong NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolMyfullstrong)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myfullstrong branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolMyfullstrong NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolMyfullstrong)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myfullstrong branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolMyfullstrong NULL
#endif


/** branching execution method for fractional LP solutions */
#if 0
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMyfullstrong)
{  /*lint --e{715}*/

	printf("mybranching rule!\n");
   return SCIP_OKAY;
}
#else
#define branchExeclpMyfullstrong NULL
#endif

/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextMyfullstrong)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myfullstrong branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextMyfullstrong NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 1
static
SCIP_DECL_BRANCHEXECPS(branchExecpsMyfullstrong)
{  /*lint --e{715}*/

	SCIP_PROBDATA*	probdata;
   SCIP_VAR**	cands;
   int	ncands;
	SCIP_NODE*	childnode_0;		/* z_j = 0 */
	SCIP_NODE*	childnode_1;		/* z_j = 1 */

	/* probdata */
	int	p;
	int	ndep;
	SCIP_VAR**	var_z;		/* [p] 01 variables */

	/* "_" means the matrix for blas */
	SCIP_Real*	orig_Q_;		/* [p*p] <- (X^t) X */
	SCIP_Real*	orig_q;		/* [p]   <- (X^t) y */
	SCIP_Real	r;

	int*	Mdep;					/* [ndep] */
	int*	groupX;				/* [ndep*p] */

	int	dim;
	SCIP_Real	RSS;			/* residual sum of square */
	SCIP_Real	RSS_new;
	SCIP_Real*	a;				/* [dim] */

	int	ublb;
	int	*Branchz;			/* [3*p] */
	int	*Branchz_new;		/* [3*p] */

	SCIP_Real*	Q_;	/* sub matrix of orig_Q_ */
	SCIP_Real*	Xy;	/* sub vector of orig_q */ 

	int*	list;			/* list of candidate variables */

	int	i,j,t,memo,ct;
	int	ind;
	int	dpv;

#if MYPARA_LOG
	printf("[myfullstrong brnaching]");
	Longline();
#endif

   /* get branching rule data */
	/*
   SCIP_BRANCHRULEDATA* branchruledata;
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
	*/
	
	/* get problem data*/
	probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

	p	=	SCIPprobdataGetNexvars(probdata);
	ndep	=	SCIPprobdataGetNdep(probdata);

	orig_Q_	=	SCIPprobdataGetQ(probdata);
	orig_q	=	SCIPprobdataGetq(probdata);
	r	=	SCIPprobdataGetr(probdata);
	var_z		=	SCIPprobdataGetVars_z(probdata);

	if( ndep ){
		Mdep		=	SCIPprobdataGetMdep(probdata);
		groupX	=	SCIPprobdataGetgroupX(probdata);
	}else{
		Mdep		=	NULL;
		groupX	=	NULL;
	}

	/* alloc */
	SCIP_CALL( SCIPallocBufferArray(scip, &list, p));
	SCIP_CALL( SCIPallocBufferArray(scip, &Branchz, 3*p));
	SCIP_CALL( SCIPallocBufferArray(scip, &Branchz_new, 3*p));

	
	GenerateZeroVecInt( p, list);
	GenerateZeroVecInt( 3*p, Branchz);

   /* get pseudo candidates (non-fixed integer variables) */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands) );
	
	for(i=0; i<ncands; i++){
		for(j=0; j<p; j++){
			if( cands[i]==var_z[j] ){
				list[j] = 1;
				break;	
			}
		}
	}

#if MYPARA_LOG
	printf("list:");
	printintv( p, list);
#endif

	/* get branching info */
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
	
	RSS = -1.0;
	ind = -1;

	for(i=0; i<p; i++){

		/* copy */
		for(j=0; j<(3*p); j++){
			Branchz_new[j] = Branchz[j];
		}

		/* 
		 * solve 
		 *   Q a = Xy
		 */

		if( list[i] == 1 ){

			Branchz_new[i] = 1;
			Branchz_new[p+i] = 0;
			
			if( ndep ){
				for(t=0; t<ndep; t++){
					memo = -1; 
					for(j=0; j<p; j++){
						if( *(groupX+(t*p)+j)==1 ){
							if( Branchz_new[j]==1 ) break;
							if( Branchz_new[p+j]==1 ) memo=j;
							if( j==Mdep[t] ){
								if( memo==-1 ){
									printf("error in branch_myfullstrong.c\n");
									stop();
								}
								*(Branchz_new+p+memo) = 0;
								*(Branchz_new+memo) = 1;
								break;
							}
						}
					}
				}
			}

			dim = p - sumint( &Branchz_new[0], p);
	
			/* alloc */
			SCIP_CALL( SCIPallocBufferArray( scip, &a, dim));
			SCIP_CALL( SCIPallocBufferArray( scip, &Q_, dim*dim));
			SCIP_CALL( SCIPallocBufferArray( scip, &Xy, dim));

			/* generate Q and Xy */
			/* Q */
			ct = 0;
			for(j=0; j<p; j++){
				if( (Branchz_new[j]==0) && (j != i) ){
					for(t=0; t<p; t++){
						if( (Branchz_new[t]==0) && (t != i ) ){
							Q_[ct++] = mat_( orig_Q_, p, j, t);
						}
					}
				}
			}

			if( ct != (dim*dim) ){
				printf("error in branch_myfullstrong.c\n");
				stop();
			}

			/* Xy */
			ct = 0;
			for(j=0; j<p; j++){
				if( (Branchz_new[j]==0) && (j != i) ){
					Xy[ct++] = orig_q[j];
				}
			}

			if( ct != dim ){
				printf("error in branch_myfullstrong.c\n");
				stop();
			}

			dpv = _dposv_( Q_, Xy, dim, a);

			if( dpv == 0 ){
				/* test */
				RSS_new = RSSvalue( dim, a, Xy, r);
				if( RSS_new > RSS ){
					RSS = RSS_new;
					ind = i;
				}
#if MYPARA_LOG
				printf("%d: RSS = %f\n", i, RSS_new);
#endif
			}

			/* free */
			SCIPfreeBufferArray(scip, &Q_);
			SCIPfreeBufferArray(scip, &Xy);
			SCIPfreeBufferArray(scip, &a);

		 }
	}

#if MYPARA_LOG
	printf("max->%dth var. \n", ind);
#endif

	if( ind == -1 ){
		/* free */
		SCIPfreeBufferArray(scip, &list);
		SCIPfreeBufferArray(scip, &Branchz);
		SCIPfreeBufferArray(scip, &Branchz_new);

   	*result = SCIP_DIDNOTRUN;
		return SCIP_OKAY;
	}

	SCIP_CALL( SCIPbranchVar( scip, var_z[ind], &childnode_0, NULL, &childnode_1));

	/* free */
	SCIPfreeBufferArray(scip, &list);
	SCIPfreeBufferArray(scip, &Branchz);
	SCIPfreeBufferArray(scip, &Branchz_new);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}
#else
#define branchExecpsMyfullstrong NULL
#endif


/*
 * branching rule specific interface methods
 */

/** creates the myfullstrong branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMyfullstrong(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create myfullstrong branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
	
   /* TODO: (optional) create branching rule specific data here */


   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsMyfullstrong) );
#if 0
   SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpMyfullstrong) );
#endif

   return SCIP_OKAY;
}

/**@file   branch_myrule.c
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "branch_myrule.h"
#include "probdata_linereg.h"
#include "set_myparameter.h"
#include "vector.h"

#define BRANCHRULE_NAME            "myrule"
#define BRANCHRULE_DESC            "mybranching rule "
#define BRANCHRULE_PRIORITY        100000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define MYPARA_LOG	0

/*
 * Data structures
 */

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
SCIP_DECL_BRANCHCOPY(branchCopyMyrule)
{  /*lint --e{715}*/

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleMyrule(scip) ) ;

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeMyrule)
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
SCIP_DECL_BRANCHINIT(branchInitMyrule)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 1
static
SCIP_DECL_BRANCHEXIT(branchExitMyrule)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}
#else
#define branchExitMyrule NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolMyrule)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myrule branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolMyrule NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolMyrule)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myrule branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolMyrule NULL
#endif


/** branching execution method for fractional LP solutions */
#if 0
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMyrule)
{  /*lint --e{715}*/

	printf("mybranching rule!\n");
   return SCIP_OKAY;
}
#else
#define branchExeclpMyrule NULL
#endif

/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextMyrule)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myrule branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextMyrule NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 1
static
SCIP_DECL_BRANCHEXECPS(branchExecpsMyrule)
{  /*lint --e{715}*/

	SCIP_PROBDATA*	probdata;
   SCIP_VAR**	cands;
   int	ncands;
	SCIP_NODE*	childnode_0;		/* z_j = 0 */
	SCIP_NODE*	childnode_1;		/* z_j = 1 */

	/* probdata */
	int	p;
	SCIP_VAR**	var_z;		/* [p] 01 variables */

	int*	list;			/* list of candidate variables */

	int	nsols;
	SCIP_SOL** sols;
	int	z_val;
	int*	n_z;			/* the number of times that z is 1 */
	int	max;
	int	max_z,ind;

	int	i,j;

#if MYPARA_LOG
	printf("[myrule brnaching]");
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
	var_z		=	SCIPprobdataGetVars_z(probdata);

	/* alloc */
	SCIP_CALL( SCIPallocBufferArray(scip, &list, p));
	SCIP_CALL( SCIPallocBufferArray(scip, &n_z, p));
	
	GenerateZeroVecInt( p, list);
	GenerateZeroVecInt( p, n_z);

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

	nsols = SCIPgetNSols(scip);
	sols = SCIPgetSols(scip);

	if( nsols == 0 ){
		printf("error in branch_myrule.c\n");
		stop();
	}

	if( nsols >= MP_NUM_SOL ) max = MP_NUM_SOL;
	else 					max = nsols;

	for(i=0; i<max; i++){
		for(j=0; j<p; j++){
			if( list[j] == 1 ){
				z_val = SCIPgetSolVal( scip, sols[i], var_z[j]);
				n_z[j] += z_val;
			}
		}
	}

	max_z = 0;
	ind = -1;
	for(i=0; i<p; i++){
		if( list[i] == 1 ){
			if( max_z < n_z[i] ){
				max_z = n_z[i];
				ind = i;
			}
			if( max_z == max ) break;
		}
	}


#if MYPARA_LOG
	printintv( p, n_z);
	printN(ind);
#endif

	if( ind == -1 ){
		for(i=0; i<p; i++){
			if( list[i] == 1 ){
				ind = i;
				break;
			}
		}
	}

	if( ind == -1 ){
		printf("error in branch_myrule.c[ind:%d]\n", ind);
		stop();
	}

	SCIP_CALL( SCIPbranchVar( scip, var_z[ind], &childnode_0, NULL, &childnode_1));

	/* free */
	SCIPfreeBufferArray(scip, &list);
	SCIPfreeBufferArray(scip, &n_z);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}
#else
#define branchExecpsMyrule NULL
#endif


/*
 * branching rule specific interface methods
 */

/** creates the myrule branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMyrule(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create myrule branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
	
   /* TODO: (optional) create branching rule specific data here */


   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyMyrule) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeMyrule) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitMyrule) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitMyrule) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsMyrule) );
#if 0
   SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolMyrule) );
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolMyrule) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextMyrule) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpMyrule) );
#endif

   return SCIP_OKAY;
}


/**@file   cmain.c
 */
#include <stdio.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "reader_linereg.h"
#include "relax_dposv.h"
#include "heur_forward.h"
#include "heur_backward.h"
#include "branch_myrule.h"
#include "branch_myfullstrong.h"
#include "probdata_linereg.h"

#include "matrix.h"
#include "vector.h"
#include "cblapack.h"
#include "linear_dependent.h"
#include "set_myparameter.h"

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */

static
SCIP_RETCODE OutputSolution(
	SCIP*				scip		
)
{
		
	SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
	int	n	=	SCIPprobdataGetNdatas(probdata);
	int	p	=	SCIPprobdataGetNexvars(probdata);
	SCIP_SOL*	sol;
	SCIP_VAR*	var;
	char			varname[SCIP_MAXSTRLEN];
	SCIP_Real 	solval;
	int			solvalint;
	SCIP_Real 	AIC;
	int			ct=0;

	int i;

   assert(probdata != NULL);

   printf("\nSolution:\n");
	
	/* print a part of solution */

	sol = SCIPgetBestSol(scip);
	AIC = (double)n*log(2.0*M_PI) 
		 + (double)n + 2.0;
	

	printf("\tz\t\t\ta\n\n");
	for(i=0; i<p; i++){
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "a%d", i+1);
		var			= SCIPfindVar(scip, varname);
		solval		= SCIPgetSolVal(scip, sol, var);
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z%d", i+1);
		var			= SCIPfindVar(scip, varname);
		solvalint 	= SCIPgetSolVal(scip, sol, var);
		printf("%2d:\t%d\t%20.15g\n", i+1, solvalint, solval);
		if( solvalint==1 ) ct++;

		AIC += 2*(double)solvalint;
	}
	

	(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "rss", NULL);
	var		=	SCIPfindVar(scip, varname);
	solval	=	SCIPgetSolVal(scip, sol, var);
	AIC 		+=	(double)n*log(solval/(double)n);
	
	printf("\nrss : %20.15g\n", solval);
	printf("AIC : %20.15g\n", AIC);
	printf("k   :\t  %d\n\n", ct);
	
	return SCIP_OKAY;
}
	

static
SCIP_RETCODE runShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array containing shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include steiner tree reader */
   SCIP_CALL( SCIPincludeReaderLinereg(scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include relaxator */
	SCIP_CALL( SCIPincludeRelaxDposv(scip) );

	/* include primal heuristics */
	SCIP_CALL( SCIPincludeHeurForward(scip));
	
	/* include primal heuristics */
	SCIP_CALL( SCIPincludeHeurBackward(scip));

	/* include branching rule */
	SCIP_CALL( SCIPincludeBranchruleMyrule(scip));

	/* include branching rule */
	SCIP_CALL( SCIPincludeBranchruleMyfullstrong(scip));

   /* set LINEREG-specific default parameters */

	SCIP_CALL( SCIPsetMyParameter( scip ));

   /**********************************
    * Process command line arguments *
    **********************************/
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   if( SCIPgetNSols(scip) > 0 )
   {
		OutputSolution(scip);
	}

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                   argc,               /**< number of shell parameters */
   char**                argv                /**< array containing shell parameters */
   )
{
   SCIP_RETCODE retcode;

   retcode = runShell(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}

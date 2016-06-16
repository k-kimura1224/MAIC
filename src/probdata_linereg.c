
/**@file   probdata_linereg.c
 * @brief  Problem data for +***** problems
 * @author 
 *
 *
 */


#include <string.h>
#include <stdio.h>

#include "probdata_linereg.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/misc.h"
#include "scip/struct_misc.h"

#include "read_data.h"
#include "normalization.h"
#include "divide_data.h"
#include "matrix.h"
#include "vector.h"
#include "cblapack.h"
#include "linear_dependent.h"

#define CENTER_OK    0           /**< do nothing */
#define CENTER_DEG   1           /**< find maximum degree */
#define CENTER_SUM   2           /**< find the minimum distance sum */
#define CENTER_MIN   3           /**< find the minimum largest distance */
#define CENTER_ALL   4           /**< find the minimum distance sum to all knots */

#define MODE_CUT    0           /**< branch and cut */
#define MODE_FLOW   1           /**< use flow model */
#define MODE_PRICE  2           /**< branch and price */


struct SCIP_ProbData
{
	int		mode;			/**< solving mode selected by the user (Cut, Price, Flow) */

   /** for FiberSCIP **/
   SCIP_Bool             ug;                 /**< indicates if this ug dual bound is set or not */
   int                   nSolvers;           /**< the number of solvers */
   SCIP_Real             ugDual;             /**< dual bound set by ug */
	/*******************/

	int				n;				/* the number of datas */
	int				p;				/*	the number of explain vars */
	SCIP_Real*		y;				/* [n] */
	SCIP_Real*		x;				/* [n*p] */
	SCIP_Real*		Q;				/* [p*p] */
	SCIP_Real*		q;				/* [p] */
	SCIP_Real		r;
	int				ndep;			/*	number of groups */
	int*				Mdep;			/* [ndep] max number in each groups */
	int*				groupX;		/* [ndep][p] -> member of group[i]  */

	SCIP_VAR**		a;				/* [p] continuous variables */
	SCIP_VAR**		z;				/* [p] 01 variables */
	SCIP_VAR**		ep;			/* [n] continuous variables */
	SCIP_VAR*		rss;			/* continuous variable, residual sum of squares */
	SCIP_VAR*		log_rss;		/* continuous variable, log(rss) */
	int				nvars;		/* number of variables */

};

/**@name Local methods
 *
 * @{
 */


/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata           /**< pointer to problem data */
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, probdata) );

   (*probdata)->ug = FALSE;
   (*probdata)->nSolvers =0;
   (*probdata)->ugDual = 0.0;

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
	int	i;
	int	n = (*probdata)->n;
	int	p = (*probdata)->p;

	SCIPdebugMessage("probdataFree \n");
	assert(scip != NULL);
	assert(probdata != NULL);

	SCIPfreeMemoryArrayNull(scip, &(*probdata)->y);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->x);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->Q);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->q);
	if( (*probdata)->ndep ){
		SCIPfreeMemoryArrayNull(scip, &(*probdata)->Mdep);
		SCIPfreeMemoryArrayNull(scip, &(*probdata)->groupX);
	}

	/* release variables */
	for(i=0; i<p; i++){
		SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->a[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->z[i]) );
	}
	for(i=0; i<n; i++){
		SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->ep[i]) );
	}
	SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->rss) );
	SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->log_rss) );

	SCIPfreeMemoryArrayNull(scip, &(*probdata)->a);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->z);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->ep);

	/* constraints */
	/*
	SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->defrss));
	SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->deflog));
	for(i=0; i<n; ++i) SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->defep[i]));
	for(i=0; i<p; ++i){
		SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->aM1[i]));
		SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->aM2[i]));
	}

	SCIPfreeMemoryArrayNull(scip, &(*probdata)->defep);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->aM1);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->aM2);

	if( (*probdata)->ndep ){
		for(i=0; i<(*probdata)->ndep; i++){
			SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->ld[i]));
		}
		SCIPfreeMemoryArrayNull(scip, &(*probdata)->ld);
	}
	*/
	
	/* free probdata */
	SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


/** create constraints (in Flow or Price Mode) */
static
SCIP_RETCODE createConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,            /**< problem data */
	int*						ng
   )
{
	int	i,j,ct;
	int	n = probdata->n;
	int	p = probdata->p;

	/* constraint */
	SCIP_CONS*	defrss;
	SCIP_CONS*	deflog;
	SCIP_CONS**	defep;
	SCIP_CONS**	aM1;
	SCIP_CONS**	aM2;

	char consname[SCIP_MAXSTRLEN];
	SCIP_Real	one = 1.0;
	SCIP_Real	minusone = -1.0;
	SCIP_Real	bigM = 100.0;


	assert(scip != NULL);
	assert(probdata != NULL);

	SCIPdebugMessage("createConstraints \n");

	/* alloc */
	SCIP_CALL( SCIPallocMemoryArray(scip, &defep, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &aM1, p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &aM2, p));
	/**
	 * quadratic constraint defrss:
	 *		ep1^2 + .. + epn^2 - rss = 0.0
	**/
	{
		SCIP_Real*	onen;
		SCIP_CALL( SCIPallocMemoryArray(scip, &onen, n));
		for(i=0; i<n; ++i) onen[i] = 1.0;

		/* create quadratic constraint */
		SCIP_CALL( SCIPcreateConsBasicQuadratic( scip, &defrss, "def_rss",
						1, &probdata->rss, &minusone, n, probdata->ep, probdata->ep, onen, 0.0, 0.0));
		
		SCIPfreeMemoryArrayNull(scip, &onen);
	}

	/**
	 * nonlinear constraint deflog:
	 * 	log( rss ) - log_rss = 0.0
	**/
	{
		SCIP_EXPR*		rssexpr;
		SCIP_EXPR*		expr;
		SCIP_VAR*		var1[1];
		SCIP_VAR*		var2[1];
		SCIP_Real		coef[1];
		SCIP_EXPRTREE*	exprtree;

		var1[0] = probdata->rss;
		var2[0] = probdata->log_rss;
		coef[0] = -1;

		/* setup expression */
		SCIP_CALL( SCIPexprCreate( SCIPblkmem(scip), &rssexpr, SCIP_EXPR_VARIDX, 0));

		/* expression for expr : log( rss ) */
		SCIP_CALL( SCIPexprCreate( SCIPblkmem(scip), &expr, SCIP_EXPR_LOG, rssexpr));

		/* expression tree from expr */
		SCIP_CALL( SCIPexprtreeCreate( SCIPblkmem(scip), &exprtree, expr, 1, 0, NULL));
		SCIP_CALL( SCIPexprtreeSetVars(exprtree, 1, var1));

		/* create nonlinear constraint for exprtree - log_rss = 0.0 */
		SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &deflog, "def_log",
						1, var2, coef, 1, &exprtree, &one, 0.0, 0.0));

		SCIP_CALL( SCIPexprtreeFree(&exprtree));

	}

	/**
	 * linear constraint defep
	 *		y - ax -ep = 0
	**/
	{
		SCIP_Real*	coef;
		SCIP_CALL( SCIPallocMemoryArray( scip, &coef, p));

		for(i=0; i<n; ++i){
			(void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "def_ep%d", i+1);
			for(j=0; j<p; ++j) coef[j] = probdata->x[i+(n*j)];
			SCIP_CALL( SCIPcreateConsBasicLinear( scip, &defep[i], consname,
							p, probdata->a, coef, probdata->y[i], probdata->y[i]));
			SCIP_CALL( SCIPaddCoefLinear( scip, defep[i], probdata->ep[i], one));
		}

		SCIPfreeMemoryArrayNull(scip, &coef);
	}

	/**
	 * linear constraint aM1:
	 * 0 <= a + Mz
	**/
	{
		for(i=0; i<p; ++i){
			(void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "aM1_%d", i+1);
			SCIP_CALL( SCIPcreateConsBasicLinear( scip, &aM1[i], consname,
							0, NULL, NULL, 0.0, SCIPinfinity(scip)));
			SCIP_CALL( SCIPaddCoefLinear( scip, aM1[i], probdata->a[i], one));
			SCIP_CALL( SCIPaddCoefLinear( scip, aM1[i], probdata->z[i], bigM));
		}
	}

	/**
	 * linear constraint aM2:
	 * a - Mz <= 0
	**/
	{
		for(i=0; i<p; ++i){
			(void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "aM2_%d", i+1);
			SCIP_CALL( SCIPcreateConsBasicLinear( scip, &aM2[i], consname,
							0, NULL, NULL, -SCIPinfinity(scip), 0.0));
			SCIP_CALL( SCIPaddCoefLinear( scip, aM2[i], probdata->a[i], one));
			SCIP_CALL( SCIPaddCoefLinear( scip, aM2[i], probdata->z[i], -bigM));
		}
	}

	/* add constraints to problem */
	SCIP_CALL( SCIPaddCons(scip, defrss));
	SCIP_CALL( SCIPaddCons(scip, deflog));
	for(i=0; i<n; ++i) SCIP_CALL( SCIPaddCons(scip, defep[i]));
	for(i=0; i<p; ++i){
		SCIP_CALL( SCIPaddCons(scip, aM1[i]));
		SCIP_CALL( SCIPaddCons(scip, aM2[i]));
	}

	/* constraints */
	SCIP_CALL( SCIPreleaseCons(scip, &defrss));
	SCIP_CALL( SCIPreleaseCons(scip, &deflog));
	for(i=0; i<n; ++i) SCIP_CALL( SCIPreleaseCons(scip, &defep[i]));
	for(i=0; i<p; ++i){
		SCIP_CALL( SCIPreleaseCons(scip, &aM1[i]));
		SCIP_CALL( SCIPreleaseCons(scip, &aM2[i]));
	}

	/**
	 * linearly dependent constraint 
	 * ex.
	 *  dep = { 1, 2, 3, 4}
	 *  add constraint:
	 *  "ld"
	 *		0 <= z1 + z2 + z3 + z4 <= 3
	**/
	if( probdata->ndep ){
		/* constraint ld */
		
		SCIP_CONS**	ld;
		SCIP_CALL( SCIPallocMemoryArray( scip, &ld, probdata->ndep));

		ct=0;
		
		for(i=0; i<(probdata->ndep); ++i){
			(void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "ld_%d", i+1);
			SCIP_CALL( SCIPcreateConsBasicLinear( scip, &ld[i], consname,
							0, NULL, NULL, 0.0, ng[i]-1));
			for(j=0; j<p; ++j){
				if( probdata->groupX[ct]==1 ){
					SCIP_CALL( SCIPaddCoefLinear( scip, ld[i], probdata->z[j], one));
				}
				ct++;
			}
			SCIP_CALL( SCIPaddCons(scip, ld[i]));
			SCIP_CALL( SCIPreleaseCons(scip, &ld[i]));
		}
		
		SCIPfreeMemoryArrayNull(scip, &ld);
	}

	SCIPfreeMemoryArrayNull(scip, &defep);
	SCIPfreeMemoryArrayNull(scip, &aM1);
	SCIPfreeMemoryArrayNull(scip, &aM2);

   return SCIP_OKAY;
}

/** create initial columns */
static
SCIP_RETCODE createVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata           /**< problem data */
   )
{
	int	i;
	int	n = probdata->n;
	int	p = probdata->p;
	char	varname[SCIP_MAXSTRLEN];

	assert(scip != NULL);
	assert(probdata != NULL);

	/* create variables */
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->a, p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->z, p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->ep, n));

	SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->rss, "rss",
					0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS));
	SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->log_rss, "log_rss",
					-SCIPinfinity(scip), SCIPinfinity(scip),
					(double)n, SCIP_VARTYPE_CONTINUOUS));

	for(i=0; i<p; ++i){
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "a%d", i+1);
		SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->a[i], varname,
						-SCIPinfinity(scip), SCIPinfinity(scip),
						0.0, SCIP_VARTYPE_CONTINUOUS));
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z%d", i+1);
		SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->z[i], varname,
						0.0, 1.0, 2.0, SCIP_VARTYPE_BINARY));
	}

	for(i=0; i<n; ++i){
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "ep%d", i+1);
		SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->ep[i], varname,
						-SCIPinfinity(scip), SCIPinfinity(scip),
						0.0, SCIP_VARTYPE_CONTINUOUS));
	}

	/* add variables to problem */
	SCIP_CALL( SCIPaddVar(scip, probdata->rss));
	SCIP_CALL( SCIPaddVar(scip, probdata->log_rss));
	for(i=0; i<p; ++i){
		SCIP_CALL( SCIPaddVar(scip, probdata->a[i]));
		SCIP_CALL( SCIPaddVar(scip, probdata->z[i]));
	}
	for(i=0; i<n; ++i) SCIP_CALL( SCIPaddVar(scip, probdata->ep[i]));

	probdata->nvars = 2 + p + p + n;

	/* branching variables */
	for(i=0; i<p; ++i){
		SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->z[i], 1000));
	}

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** copies user data of source SCIP for the target SCIP */
static
SCIP_DECL_PROBCOPY(probcopyLinereg)
{

	int	n		=	sourcedata->n;
	int	p		=	sourcedata->p;
	int	ndep	=	sourcedata->ndep;
	int	i;
	SCIP_Bool success;

	SCIPdebugMessage("########################## probcopy ###########################\n");

	SCIP_CALL( probdataCreate(scip, targetdata) );

	(*targetdata)->n		=	sourcedata->n;
	(*targetdata)->p		=	sourcedata->p;
	(*targetdata)->ndep	=	sourcedata->ndep;


	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->y, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->x, n*p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->Q, p*p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->q, p));
	(*targetdata)->r = sourcedata->r;

	if( ndep ){
		SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->Mdep, ndep));
		SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->groupX, ndep*p));
	}

	for(i=0; i<n; i++)
		(*targetdata)->y[i] = sourcedata->y[i];

	for(i=0; i<n*p; i++)
		(*targetdata)->x[i] = sourcedata->x[i];

	for(i=0; i<p*p; i++)
		(*targetdata)->Q[i] = sourcedata->Q[i];

	for(i=0; i<p; i++)
		(*targetdata)->q[i] = sourcedata->q[i];

	if( ndep ){
		for(i=0; i<ndep; i++)
			(*targetdata)->Mdep[i] = sourcedata->Mdep[i];

		for(i=0; i<(ndep*p); i++)
			(*targetdata)->groupX[i] = sourcedata->groupX[i];
	}else{
		(*targetdata)->Mdep = NULL;
		(*targetdata)->groupX = NULL;
	}

	/* variables */

	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->a, p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->z, p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->ep, n));

	for(i=0; i<p; i++){
		SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->a[i],
			&((*targetdata)->a[i]), varmap, consmap, global, &success) );
		assert(success);
		SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->a[i]) );
	
		SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->z[i],
			&((*targetdata)->z[i]), varmap, consmap, global, &success) );
		assert(success);
		SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->z[i]) );
	}
	
	for(i=0; i<n; i++){
		SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->ep[i],
			&((*targetdata)->ep[i]), varmap, consmap, global, &success) );
		assert(success);
		SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->ep[i]) );
	}

	SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->rss,
		&((*targetdata)->rss), varmap, consmap, global, &success) );
	assert(success);
	SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->rss) );

	SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->log_rss,
		&((*targetdata)->log_rss), varmap, consmap, global, &success) );
	assert(success);
	SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->log_rss) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigLinereg)
{
   SCIPdebugMessage("probdelorigLinereg \n");

   SCIPdebugMessage("free original problem data\n");

   /* free the (original) probdata */
   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransLinereg)
{
	SCIP_Real timelimit;
	SCIP_Bool update;

	int	n		=	sourcedata->n;
	int	p		=	sourcedata->p;
	int	ndep	=	sourcedata->ndep;
	int	i;

	SCIPdebugMessage("probtransLinereg \n");

	SCIP_CALL( SCIPgetBoolParam(scip, "linereg/countpresoltime", &update) );

	/* adjust time limit to take into account reading time */
	if( update )
	{
	   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
	   timelimit -= SCIPgetReadingTime(scip);
	   timelimit = MAX(0.0,timelimit);
	   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit) );
	}

	/* create transform probdata */
	SCIP_CALL( probdataCreate(scip, targetdata) );

	(*targetdata)->n		=	sourcedata->n;
	(*targetdata)->p		=	sourcedata->p;
	(*targetdata)->ndep	=	sourcedata->ndep;


	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->y, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->x, n*p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->Q, p*p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->q, p));
	(*targetdata)->r = sourcedata->r;
	if( ndep ){
		SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->Mdep, ndep));
		SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->groupX, ndep*p));
	}

	for(i=0; i<n; i++)
		(*targetdata)->y[i] = sourcedata->y[i];

	for(i=0; i<n*p; i++)
		(*targetdata)->x[i] = sourcedata->x[i];

	for(i=0; i<p*p; i++)
		(*targetdata)->Q[i] = sourcedata->Q[i];

	for(i=0; i<p; i++)
		(*targetdata)->q[i] = sourcedata->q[i];

	if( ndep ){
		for(i=0; i<ndep; i++)
			(*targetdata)->Mdep[i] = sourcedata->Mdep[i];

		for(i=0; i<(ndep*p); i++)
			(*targetdata)->groupX[i] = sourcedata->groupX[i];
	}else{
		(*targetdata)->Mdep = NULL;
		(*targetdata)->groupX = NULL;
	}

	/* variables */
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->a, p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->z, p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->ep, n));

	SCIP_CALL( SCIPtransformVars(scip, p, sourcedata->a, (*targetdata)->a) );
	SCIP_CALL( SCIPtransformVars(scip, p, sourcedata->z, (*targetdata)->z) );
	SCIP_CALL( SCIPtransformVars(scip, n, sourcedata->ep, (*targetdata)->ep) );
	SCIP_CALL( SCIPtransformVars(scip, 1, &sourcedata->rss, &(*targetdata)->rss) );
	SCIP_CALL( SCIPtransformVars(scip, 1, &sourcedata->log_rss, &(*targetdata)->log_rss) );

   return SCIP_OKAY;
}


static
SCIP_DECL_PROBEXITSOL(probexitsolLinreg)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(probdata != NULL);

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransLinereg)
{
   SCIPdebugMessage("free transformed problem data\n");

	SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** get problem name 
 *
 *  Return NULL on error
 */
static
SCIP_RETCODE getProblemName(
	const char*	filename,	/*	input filename			  */
	char*			probname,	/*	output problemname	  */
	int			maxSize		/* maximum size of p.name */
	)
{
	int	i=0;
	int	j=0;
	int	l;

	/*	first find end of string */
	while( filename[i]!=0 )
		++i;
	l = i;

	/* go back until '.' or '/' or '\' appears */
	while( (i>0) && (filename[i]!='.') && (filename[i]!='/') && (filename[i]!='\\'))
		--i;

	/* if we found '.', search for '/' or '\\' */
	if( filename[i]=='.' ){
		l = i;
		while( (i>0) && (filename[i]!='/') && (filename[i]!='\\') )
			--i;
	}

	/* crrect counter */
	if( (filename[i]=='/') || (filename[i]=='\\') )
		++i;

	/* copy name */
	while( (i<l) && (filename[i]!=0) ){
		probname[j++] = filename[i++];
		if( j>maxSize-1)
			return SCIP_ERROR;
	}
	probname[j] = 0;

	return SCIP_OKAY;

}

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
	SCIP_PROBDATA* probdata;
	int	n, p, i_ex, ndep;
	char	mode;
	SCIP_Real*	data;	/* [n*(p+1)] */
	SCIP_Real*	X;		/* [n*p] */
	int*	D;				/* [p] */
	int*	ng;			/* ng: number of groupX[i] */


	int	i;
	int	ct=0;
	char	probname[SCIP_MAXSTRLEN];

	/* create problem data */
	 SCIP_CALL( probdataCreate(scip, &probdata) );
	
	/* read data */
	ReadDim( filename, &n, &p, &i_ex);
	SCIP_CALL( SCIPallocMemoryArray( scip, &data, n*(p+1)));
	ReadData( filename, n, p, data);
	probdata->n = n;
	probdata->p = p;
	
	/* adjust data	*/
	Normalization( n, p, data);

	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->y, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &X, n*p));
	DataToExpl( n, p, i_ex, data, probdata->y, X);

	SCIPfreeMemoryArrayNull( scip, &data);

	/* output information of problem */
	SCIP_CALL( getProblemName( filename, probname, SCIP_MAXSTRLEN));

	SCIPinfoMessage( scip, NULL, "File name\t:\t%s\n", filename);
	SCIPinfoMessage( scip, NULL, "Problem name \t:\t%s\n", probname);
	SCIPinfoMessage( scip, NULL, "Number of data\t:\t%d\n", n);
	SCIPinfoMessage( scip, NULL, "Number of var\t:\t%d\n", p);

	/**
	 * for using lapack and blas,
	 * define
	 *		x[n*p]	:=	X[n][p]
	 *		Q[p*p]	:=	(X^t)(X)
	 *		q[p]		:=	(X^t)(y)
	 *		r			:=	(y^t)(y)
	**/


	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->x, n*p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->Q, p*p));
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->q, p));
	GenerateZeroVec( p*p, probdata->Q);
	GenerateZeroVec( p, probdata->q);
	
	TraMat( n, p, X, probdata->x);		
	dgemm_t( probdata->x, n, p, probdata->Q); 
	dgemv_t( probdata->x, n, p, probdata->y, probdata->q);
	probdata->r = myddot_( probdata->y, probdata->y, n);

	SCIPfreeMemoryArrayNull( scip, &X);

	/*	linear dependent */
	SCIP_CALL( SCIPallocMemoryArray(scip, &D, p));
	
	LinearDependent( n, p, probdata->x, D);
	ndep	=	sumint( D, p);
	probdata->ndep	=	ndep;

	if( ndep ){
		SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->Mdep, ndep));
		ct=0;

		for(i=0; i<p; ++i){
			if( D[i]==1 ){
				probdata->Mdep[ct++] = i;
			}
		}
		
		SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->groupX, ndep*p));
		GenerateZeroMatInt( ndep, p, probdata->groupX);
		
		
		for(i=0; i<ndep; ++i){
			*(probdata->groupX+(i*p)+probdata->Mdep[i]) = 1;
		}

		LinearDependentGroup( p, probdata->Q, ndep, probdata->Mdep, D, probdata->groupX);
		printLineDepGroup( ndep, p, probdata->groupX);

		SCIP_CALL( SCIPallocBufferArray(scip, &ng, ndep));
		for(i=0; i<ndep; ++i){
			ng[i] = sumint( &(probdata->groupX[i*p]), p);
		}

	}else{
		printf("linear independent\n\n");

		probdata->Mdep		=	NULL;
		probdata->groupX	=	NULL;
		ng						=	NULL;
	}

	SCIPfreeMemoryArrayNull(scip, &D);

   /* get parameters */
	SCIP_CALL( SCIPgetCharParam(scip, "linereg/mode", &mode) );

	/* create a problem in SCIP and add non-NULL callbacks via setter functions */
	SCIP_CALL( SCIPcreateProbBasic(scip, probname) );
	SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigLinereg) );
	SCIP_CALL( SCIPsetProbTrans(scip, probtransLinereg) );
	SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransLinereg) );
	SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolLinreg) );
	SCIP_CALL( SCIPsetProbCopy(scip, probcopyLinereg) );

	/* set objective sense */
	SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

	/* set user problem data */
	SCIP_CALL( SCIPsetProbData(scip, probdata) );

	/* disable sub-SCIP heuristics */
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rens/freq", -1) );
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rins/freq", -1) );
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/dins/freq", -1) );
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/crossover/freq", -1) );
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/mutation/freq", -1) );
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/clique/freq", -1) );
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/vbounds/freq", -1) );
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/bound/freq", -1) );
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/zeroobj/freq", -1) );

	/* create and add variables */
	SCIP_CALL( createVariables(scip, probdata) );

	/* create and add constraints */
	SCIP_CALL( createConstraints(scip, probdata, ng) );

	/* free */
	if( ndep ){
		SCIPfreeBufferArray(scip, &ng);
	}

   return SCIP_OKAY;
}

/** returns the number of datas */
int SCIPprobdataGetNdatas(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->n;
}

/** returns the number of explain vars */
int SCIPprobdataGetNexvars(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->p;
}

/** returns the number of linedep groups */
int SCIPprobdataGetNdep(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->ndep;
}

/** returns y */
extern
SCIP_Real*	SCIPprobdataGety(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->y;
}

/** returns X */
extern
SCIP_Real*	SCIPprobdataGetX(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->x;
}

/** returns Q */
extern
SCIP_Real*	SCIPprobdataGetQ(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->Q;
}

/** returns q */
extern
SCIP_Real*	SCIPprobdataGetq(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->q;
}

/** returns r */
extern
SCIP_Real	SCIPprobdataGetr(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->r;
}

/** returns Mdep */
extern
int*	SCIPprobdataGetMdep(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->Mdep;
}

/** returns groupX */
extern
int*	SCIPprobdataGetgroupX(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->groupX;
}

/** returns var_a */
extern
SCIP_VAR**	SCIPprobdataGetVars_a(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->a;
}

/** returns var_z */
extern
SCIP_VAR**	SCIPprobdataGetVars_z(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->z;
}

/** returns var_ep */
extern
SCIP_VAR**	SCIPprobdataGetVars_ep(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->ep;
}

/** returns var_rss */
extern
SCIP_VAR*	SCIPprobdataGetVar_rss(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->rss;
}

/** returns var_log */
extern
SCIP_VAR*	SCIPprobdataGetVar_log(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->log_rss;
}

/** writes end of log file */
/*
SCIP_RETCODE SCIPprobdataWriteLogfileEnd(
   SCIP*                 scip     
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   if( probdata->logfile != NULL )
   {
      int success;
      SCIP_Real factor = 1.0;

      if( probdata->stp_type ==  STP_MAX_NODE_WEIGHT )
         factor = -1.0;

      SCIPprobdataWriteLogLine(scip, "End\n");
      SCIPprobdataWriteLogLine(scip, "\n");
      SCIPprobdataWriteLogLine(scip, "SECTION Run\n");
      if( probdata->ug )
      {
         SCIPprobdataWriteLogLine(scip, "Threads %d\n", probdata->nSolvers);
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
         SCIPprobdataWriteLogLine(scip, "Dual %16.9f\n", factor * probdata->ugDual);
      }
      else
      {
         SCIPprobdataWriteLogLine(scip, "Threads 1\n");
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
         SCIPprobdataWriteLogLine(scip, "Dual %16.9f\n", factor * SCIPgetDualbound(scip));
      }
      SCIPprobdataWriteLogLine(scip, "Primal %16.9f\n", factor * SCIPgetPrimalbound(scip));
      SCIPprobdataWriteLogLine(scip, "End\n");

      if( SCIPgetNSols(scip) > 0 )
      {
         SCIPprobdataWriteLogLine(scip, "\n");
         SCIPprobdataWriteLogLine(scip, "SECTION Finalsolution\n");

         SCIP_CALL( SCIPprobdataWriteSolution(scip, probdata->logfile) );
         SCIPprobdataWriteLogLine(scip, "End\n");
      }

      success = fclose(probdata->logfile);
      if( success != 0 )
      {
         SCIPerrorMessage("An error occurred while closing file <%s>\n", probdata->logfile);
         return SCIP_FILECREATEERROR;
      }

      probdata->logfile = NULL;
   }


   return SCIP_OKAY;
}
*/

/** writes end of log file */
void SCIPprobdataSetDualBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dual                /**< dual bound */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   probdata->ug = TRUE;
   probdata->ugDual = dual;
}

/** writes end of log file */
void SCIPprobdataSetNSolvers(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nSolvers            /**< the number of solvers */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   probdata->nSolvers = nSolvers;
}

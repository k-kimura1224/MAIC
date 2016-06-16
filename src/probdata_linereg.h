
/**@file   probdata_linereg.h
 * @brief  Problem data for stp problem
 * @author 
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_LINEREG__
#define __SCIP_PROBDATA_LINEREG__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** sets up the problem data */
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   );

/** returns the number of datas */
extern
int SCIPprobdataGetNdatas(
   SCIP_PROBDATA*			probdata
	);

/** returns the number of explain vars */
extern
int SCIPprobdataGetNexvars(
   SCIP_PROBDATA*			probdata
	);

/** returns the number of linedep groups */
extern
int SCIPprobdataGetNdep(
   SCIP_PROBDATA*			probdata
	);

/** returns y */
extern
SCIP_Real*	SCIPprobdataGety(
   SCIP_PROBDATA*			probdata
	);

/** returns X */
extern
SCIP_Real*	SCIPprobdataGetX(
   SCIP_PROBDATA*			probdata
	);

/** returns Q */
extern
SCIP_Real*	SCIPprobdataGetQ(
   SCIP_PROBDATA*			probdata
	);

/** returns q */
extern
SCIP_Real*	SCIPprobdataGetq(
   SCIP_PROBDATA*			probdata
	);

/** returns r */
extern
SCIP_Real	SCIPprobdataGetr(
   SCIP_PROBDATA*			probdata
	);

/** returns Mdep */
extern
int*	SCIPprobdataGetMdep(
   SCIP_PROBDATA*			probdata
	);

/** returns groupX */
extern
int*	SCIPprobdataGetgroupX(
   SCIP_PROBDATA*			probdata
	);

/** returns var_a */
extern
SCIP_VAR**	SCIPprobdataGetVars_a(
   SCIP_PROBDATA*			probdata
	);

/** returns var_z */
extern
SCIP_VAR**	SCIPprobdataGetVars_z(
   SCIP_PROBDATA*			probdata
	);

/** returns var_ep */
extern
SCIP_VAR**	SCIPprobdataGetVars_ep(
   SCIP_PROBDATA*			probdata
	);

/** returns var_rss */
extern
SCIP_VAR*	SCIPprobdataGetVar_rss(
   SCIP_PROBDATA*			probdata
	);

/** returns var_log */
extern
SCIP_VAR*	SCIPprobdataGetVar_log(
   SCIP_PROBDATA*			probdata
	);

/** set dual bound by ug */
extern
void SCIPprobdataSetDualBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dual
   );

/** set the number of solvers */
extern
void SCIPprobdataSetNSolvers(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nSolvers            /**< the number of solvers */
   );

#ifdef __cplusplus
}
#endif

#endif

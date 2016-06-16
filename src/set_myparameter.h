
/**@file   set_myparameter.h
 * @ingroup 
 * @brief  
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SET_MYPARAMETER_H__
#define __SCIP_SET_MYPARAMETER_H__

#include "scip/scip.h"

/**
 * define the value of my parameter
**/
#define MP_RESTARTDFS	300000	/* default:10000 */
#define MP_TIMELIMIT		5000.0
#define MP_NUM_SOL		100		/* number of stored solution (df:10) */
#define MP_LPRELAX		-1			/* if -1, disable the LP relaxation */
#define MP_DFPRIMALHEUR	-1			/* if -1, disable the all default heuristic */
#define MP_PRINTORIGP	0			/* if nonzero, print original problem.  */
#define MP_PRINTREFOP	0			/* if nonzero, print reformulated problem */

#ifdef __cplusplus
extern "C" {
#endif
/** creates  and includes in SCIP */
EXTERN
SCIP_RETCODE SCIPsetMyParameter(
   SCIP*                scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif

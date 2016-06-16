
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

#ifndef __SCIP_HEUR_FORWARD_H__
#define __SCIP_HEUR_FORWARD_H__


#include "scip/scip.h"
#ifdef __cplusplus
extern "C" {
#endif

/** creates the local primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurForward(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif

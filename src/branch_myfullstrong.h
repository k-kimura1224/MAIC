
/**@file   branch_myfullstrong.h
 * @ingroup BRANCHINGRULES
 * @brief mybranchingrule  - full strong -
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_MYFULLSTRONG_H__
#define __SCIP_BRANCH_MYFULLSTRONG_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the xyz branching rule and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleMyfullstrong(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif

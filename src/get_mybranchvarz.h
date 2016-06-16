
/**@file     get_mybranchvarz.h
 * @ingroup 
 * @brief  
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NAME_H__
#define __SCIP_NAME_H__


#include "scip/scip.h"
/**
 * in parent node,
 * branching varizble z is returned
**/

#ifdef __cplusplus
extern "C" {
#endif

/** creates  and includes in SCIP */
extern
int SCIPparentGetBranchingZ(
	SCIP*			scip,
	SCIP_NODE*	node,
	int			p,	
	SCIP_VAR**	var_z
	);

#ifdef __cplusplus
}
#endif

#endif

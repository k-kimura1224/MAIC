
/**
 * my prameter(MP) is defined 
 * in set_myparameter.h
**/

#include "set_myparameter.h"

SCIP_RETCODE SCIPsetMyParameter(
	SCIP*		scip
	)
{
	
	/* if -1, disable the LP relaxation and only use my custom relaxation */
	SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", MP_LPRELAX));

	/* depth first search with periodical selection of the best node */
	SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/restartdfs/stdpriority", MP_RESTARTDFS));

	/* timelimit */
	SCIP_CALL( SCIPsetRealParam(scip, "limits/time", MP_TIMELIMIT));

	/**
	 * maximal number of solutions candidates to store 
	 * in the solution storage of the original problem [default:10]
	**/
	SCIP_CALL( SCIPsetIntParam(scip, "limits/maxorigsol", MP_NUM_SOL));
	SCIP_CALL( SCIPsetIntParam(scip, "limits/maxsol", MP_NUM_SOL));

	/* if -1, disable the all primal heuristic */
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/actconsdiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/clique/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/coefdiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/crossover/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/dins/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/dualval/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/feaspump/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/fixandinfer/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/fracdiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/guideddiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/intdiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/intshifting/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/linesearchdiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/localbranching/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/mutation/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/nlpdiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/objpscostdiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/octane/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/oneopt/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/proximity/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/pscostdiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/randrounding/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rens/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rins/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rootsoldiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rounding/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/shiftandpropagate/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/shifting/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/simplerounding/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/trivial/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/twoopt/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/undercover/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/vbounds/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/veclendiving/freq", MP_DFPRIMALHEUR));
	SCIP_CALL( SCIPsetIntParam(scip, "heuristics/zeroobj/freq", MP_DFPRIMALHEUR));
	

	return SCIP_OKAY;
}

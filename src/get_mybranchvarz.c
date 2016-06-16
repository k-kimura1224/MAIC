 
/**
 * in parent node,
 * branching varizble z is returned
 */

#include "get_mybranchvarz.h"
#include "vector.h"


extern
int SCIPparentGetBranchingZ(
	SCIP*			scip,
	SCIP_NODE*	node,
	int			p,	
	SCIP_VAR**	var_z
	)
{
	char	name[SCIP_MAXSTRLEN];
	
	int i;
	int j;
	SCIP_VAR**            branchvars;         /* array of variables on which the branchings has been performed in all ancestors */
	SCIP_VAR*	lastbranchvar = NULL;
   SCIP_Real*            branchbounds;       /* array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes;         /* array of boundtypes which the branchings in all ancestors set */
   int*                  nodeswitches;       /* marks, where in the arrays the branching decisions of the next node on the path start
                                              * branchings performed at the parent of node always start at position 0. For single variable branching,
                                              * nodeswitches[i] = i holds */
   int                   nbranchvars;        /* number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int                   branchvarssize;     /* available slots in arrays */
   int                   nnodes;             /* number of nodes in the nodeswitch array */
   int                   nodeswitchsize;     /* available slots in node switch array */
	int						result;

   branchvarssize = SCIPnodeGetDepth(node);
   nodeswitchsize = branchvarssize;
   
   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &branchbounds, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodeswitches, nodeswitchsize) );

   SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize );
   
   /* if the arrays were to small, we have to reallocate them and recall SCIPnodeGetAncestorBranchingPath */
   if( nbranchvars > branchvarssize || nnodes > nodeswitchsize )
   {
      branchvarssize = nbranchvars;
      nodeswitchsize = nnodes;

      /* memory reallocation */
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchvars, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchbounds, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &boundtypes, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &nodeswitches, nodeswitchsize) );
      
      SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize);
      assert(nbranchvars == branchvarssize);
   }
  
   /* we only want to create output, if branchings were performed */
   if( nbranchvars >= 1 )
   {

      /* print all nodes, starting from the root, which is last in the arrays */
      for( j = nnodes-1; j >= 0; --j)
      {
         int end;
         if(j == nnodes-1){
            end =  nbranchvars;
         }else{ 
            end =  nodeswitches[j+1];
         }

         for( i = nodeswitches[j]; i < end; ++i )
         {
            if( i > nodeswitches[j] ){
					/* error,if two branchings at one node */
					SCIPerrorMessage("It is unexpected error1 in get_mybranchvarz ");
					stop();
            }
			}

			if( j==0 ){
				(void) SCIPsnprintf(name, SCIP_MAXSTRLEN,SCIPvarGetName(branchvars[nodeswitches[j]]),NULL);
				lastbranchvar = branchvars[nodeswitches[j]];
			}

		}
   }
   
	result = -1;
	for(i=0; i<p; i++){
		if( lastbranchvar == var_z[i] ){
			result = i;
			break;
		}
	}

#if 0
	if( result == -1 ){
		printf("%s\n", name);
	}
		printf("%s\n", name);
#endif

   /* free all local memory */
   SCIPfreeBufferArray(scip, &nodeswitches);
   SCIPfreeBufferArray(scip, &boundtypes);
   SCIPfreeBufferArray(scip, &branchbounds);
   SCIPfreeBufferArray(scip, &branchvars);

	return result;

}

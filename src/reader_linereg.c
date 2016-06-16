
/**@file   reader_linereg.c
 * @brief  **** problem reader file reader
 * @author 
 *
 * This file implements the reader used to read and write ******* problems.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_linereg.h"
#include "reader_linereg.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "lineregreader"
#define READER_DESC             "file reader for linear regression data format"
#define READER_EXTENSION        "linereg"

#define   DEFAULT_COMPCENTRAL  1             /**< selection type for the root (for undirected STPs) */
#define   DEFAULT_EMITGRAPH    FALSE         /**< emit graph? */
#define   DEFAULT_COUNTPRESOLTIME  TRUE      /**< count presolving time as part of overall solution time? */
#define   DEFAULT_REDUCTION    1             /**< reduction mode to apply */
#define   DEFAULT_MINELIMS     5             /**< minimal number of eliminations to be achieved for reiteration of reduction methods */
#define   DEFAULT_PRETIMELIMIT -1.0          /**< presolving time limit */

#define LINEREG_MODES "cfp" /**< valid values for user parameter 'linereg/mode' */

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyLinereg)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderLinereg(scip) );

   return SCIP_OKAY;
}

/** problem reading method of the reader */
static
SCIP_DECL_READERREAD(readerReadLinereg)
{  /*lint --e{715}*/
   SCIP_RETCODE          retcode;
   SCIP_PROBDATA*        probdata;
   char                  mode;

   *result = SCIP_DIDNOTRUN;

   /* get solving mode parameter */
   SCIP_CALL( SCIPgetCharParam(scip, "linereg/mode", &mode) );

   retcode = SCIPprobdataCreate(scip, filename);

   if( retcode == SCIP_READERROR )
      return SCIP_READERROR;

   SCIP_CALL( retcode );

   probdata = SCIPgetProbData(scip);
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT ||  probdata == NULL )
      return SCIP_READERROR;

#if 0
		SCIPinfoMessage(scip, NULL, "Original problem:\n");
		SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );
   	SCIPinfoMessage(scip, NULL, "\n");
#endif

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** problem writing method of the reader */
static
SCIP_DECL_READERWRITE(readerWriteLinereg)
{  /*lint --e{715}*/

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** includes the stp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderLinereg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyLinereg) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadLinereg) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteLinereg) );

   /* include user parameters */

	SCIP_CALL( SCIPaddBoolParam(scip,
         "linereg/countpresoltime",
        "count presolving time to solving time?",
         NULL, FALSE, DEFAULT_COUNTPRESOLTIME, NULL, NULL) );

	SCIP_CALL( SCIPaddCharParam(scip,
	 "linereg/mode",
         "Solving mode: 'c'ut, 'f'low ,'p'rice",
         NULL, FALSE, 'c', LINEREG_MODES, NULL, NULL) );

   return SCIP_OKAY;
}

/**@} */

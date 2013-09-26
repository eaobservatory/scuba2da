#ifndef HEADGEN____sc2dalib_h
#define HEADGEN____sc2dalib_h 
 
 
/*+ my_fclose
*/
void my_fclose(FILE **fp);


/**
 * my_closeFiles
 * Closes myInfo-> (fpData, fpMcecmd, fpBatch, Strchart, fpOtheruse)
 *  Flushes log file.
 *  Stat: used 2 times
 */
/*+ my_closeFiles
*/
void my_closeFiles(dasInfoStruct_t *myInfo);

// =======sc2dalib_a*******
//====================//

/**
 * \fn void sc2dalib_actionfileEnd(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *     int flag, StatusType *status)
 *
 * \brief function
 *  end Drama action and close file 
 *
 * \param con     SDSU_CONTXT struct pointer
 * \param myInfo dasInfoStruct_t struct pointer
 * \param flag   int  1: check IN_SEQ, 0, don't
 * \param status  StatusType.  given and return
 */
/*+ sc2dalib_actionfileEnd 
*/
void sc2dalib_actionfileEnd
(
SDSU_CONTEXT    *con,           
dasInfoStruct_t *myInfo,      
int             flag,
StatusType      *status
);

/*+ sc2dalib_allctsharedMem   
 */
void sc2dalib_allctsharedMem
(
SDSU_CONTEXT       *con,         
dasInfoStruct_t    *myInfo,
int                startSeq,
int                endSeq,
StatusType         *status
);

/*+ sc2dalib_batchInit
*/
void sc2dalib_batchInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_batchParser
*/
void sc2dalib_batchParser
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,   
dasCmdInfo_t          *mycmdInfo, 
char                  *dateTime,
struct mcexml_struct  *mceInxpt,
StatusType            *status
);

/*+ sc2dalib_callbkdrvMsg 
*/
void sc2dalib_callbkdrvMsg
(
USED4DITS      *argptr,    
StatusType     *status
);

/*+ sc2dalib_callbkMsg 
*/
void sc2dalib_callbkMsg
(
USED4DITS    *argptr,    
StatusType   *status
);

/*+ sc2dalib_closesharedMem   
 */
void sc2dalib_closesharedMem
(
dasInfoStruct_t    *myInfo,
int                whichShared,
StatusType         *status
);

/*+ sc2dalib_changeSQVal
 */
void sc2dalib_changeSQVal
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                  *dateTime,
ARRAYSET              *setup,
StatusType            *status
);

/*+ sc2dalib_Cmd 
 */
void sc2dalib_Cmd
(
SDSU_CONTEXT      *con,
dasInfoStruct_t   *myInfo,
dasCmdInfo_t      *mycmdInfo,
char              *dateTime,
StatusType        *status
);

/*+ sc2dalib_chkChecksum 
*/
void sc2dalib_chkChecksum
(
dasInfoStruct_t  *myInfo,   
dasCmdInfo_t     *mycmdInfo,
StatusType       *status
);

/*+ sc2dalib_chksdsID 
*/
void sc2dalib_chksdsID
(
dasInfoStruct_t  *myInfo,   
StatusType       *status
);

/*+ sc2dalib_coadddatawriteDisp   
 */
void sc2dalib_coadddatawriteDisp
(
SDSU_CONTEXT      *con,
dasInfoStruct_t   *myInfo,
char              *scanfile,
char              *coaddPtr,
double            *imagedataPtr,    
StatusType        *status
);

/*+ sc2dalib_configInit
*/
void sc2dalib_configInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
struct mcexml_struct  *mceInxpt, 
char                  *dateTime,
char                  *mode,
char                  *config,
StatusType            *status
);

/*+ sc2dalib_configInitENG
*/
void sc2dalib_configInitENG
(
dasInfoStruct_t       *myInfo,    
char                  *mode,
char                  *config,
StatusType            *status
);

/*+ sc2dalib_configInitRTSC
*/
void sc2dalib_configInitRTSC
(
dasInfoStruct_t       *myInfo,    
char                  *config,
StatusType            *status
);

/*+ sc2dalib_configInitSet
*/
void sc2dalib_configInitSet
(
dasInfoStruct_t       *myInfo,  
char                  *dateTime,  
char                  *mode,
char                  *config,
StatusType            *status
);

/*+ sc2dalib_createSDP
 */
void sc2dalib_createSDP
(
dasInfoStruct_t  *myInfo,
StatusType       *status
);

/*+ sc2dalib_createQLSDS
 */
void sc2dalib_createQLSDS
(
dasInfoStruct_t  *myInfo,
StatusType       *status
);

/*+ sc2dalib_updateQLSDS
 */
void sc2dalib_updateQLSDS
(
SDSU_CONTEXT     *con,
dasInfoStruct_t  *myInfo,
int              *coaddnumPtr,
double           *timePtr,
char             *fitsPtr,
uint32           *fitsDim,    // fits dimension [0] size_n, [1] num_n
double           *coadddataPtr,
uint32           *dataDim,      // data dimension [0] col_n, [1] row_n
char             *scanfileName,
StatusType       *status
);

/*+ sc2dalib_cnvtbatch2Cmd
*/
void sc2dalib_cnvtbatch2Cmd
(
char            *setFile,
char            *cmd,
StatusType      *status
);

/*+ sc2dalib_dispInfo
*/
void sc2dalib_dispInfo
(
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
int             inSeq,
StatusType      *status
);

/*+ sc2dalib_dispResults - 
*/
void sc2dalib_dispResults
(
dasInfoStruct_t       *myInfo,   
dasCmdInfo_t          *mycmdInfo, 
StatusType            *status
);

/*+ sc2dalib_dispReply - 
*/
void sc2dalib_dispReply
(
dasCmdInfo_t  *mycmdInfo,
int           howMany, 
StatusType    *status
);

/*+ sc2dalib_downld2pciInit
*/
void sc2dalib_downld2pciInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *file,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_download2PCI
*/
void sc2dalib_download2PCI
(
SDSU_CONTEXT *con, 
char        *lodfile,
StatusType *status
);

/*+ sc2dalib_endAction - 
*/
void sc2dalib_endAction
(
SDSU_CONTEXT    *con,           
dasInfoStruct_t *myInfo,       
StatusType      *status
);

/*+ int sc2dalib_finddateTime 
 */
int sc2dalib_finddateTime
(
short         which,
char          *dateArray
);

/*+ sc2dalib_frametakeInit   
 */
void sc2dalib_frametakeInit
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
int                   startSeq,
int                   endSeq,
struct mcexml_struct  *mceInxpt,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_fluxJump
 */
void sc2dalib_fluxJump
(
dasInfoStruct_t  *myInfo,
ARRAYSET         *setup,
int              *array,
int              ch,
StatusType       *status
);

/*+ sc2dalib_getcmdBuf 
 */
void sc2dalib_getcmdBuf
(
dasInfoStruct_t        *myInfo,
dasCmdInfo_t           *mycmdInfo,
struct mcexml_struct   *mceInxpt,
StatusType             *status
);

/*+  sc2dalib_healthStatus
 */
void sc2dalib_healthStatus
(
int        which,
int        health,
StatusType *status
);

/*+ sc2dalib_initHeatBiasHandling
*/
void sc2dalib_initHeatBiasHandling
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,   
int                   *flag,
StatusType            *status
);

/*+ sc2dalib_initInit
*/
void sc2dalib_initInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *xmlfilename,
StatusType            *status
);

/*+ sc2dalib_initInitENG
*/
void sc2dalib_initInitENG
(
SDSU_CONTEXT          *con, 
dasInfoStruct_t       *myInfo, 
char                  *xmlfilename,
StatusType            *status
);

/*+ sc2dalib_initInitRTSC
*/
void sc2dalib_initInitRTSC
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *xmlfilename,
StatusType            *status
);

/*+ sc2dalib_initInitReadXML
*/
void sc2dalib_initInitReadXML
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *xmlfilename,
StatusType            *status
);

/*+ sc2dalib_initgetArrayID
*/
void sc2dalib_initgetArrayID
(
SDSU_CONTEXT          *con,
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceinxPtr,
char                  *name,
char                  *mcexml,
int                   *mcePort,
StatusType            *status
);

/*+ sc2dalib_isFile 
*/
void sc2dalib_isFile
(
const char *fname,
StatusType *status
);

/*+ sc2dalib_heaterslopeInit
*/
void sc2dalib_heaterslopeInit
(
SDSU_CONTEXT      *con,
dasInfoStruct_t   *myInfo, 
char              *dateTime,
StatusType         *status
);

/*+ sc2dalib_heaterSlope
*/
void sc2dalib_heaterSlope
(
dasInfoStruct_t  *myInfo,
int              *heaterMask,
double           *heaterSlope,
StatusType        *status
);

/*+ sc2dalib_heaterslopeRead
*/
void sc2dalib_heaterslopeRead
(
dasInfoStruct_t   *myInfo,    
double            *pixelSlope,
StatusType        *status
);

/*+ sc2dalib_heaterslopeSave
*/
void sc2dalib_heaterslopeSave
(
dasInfoStruct_t *myInfo,
double          *heaterSlope,
StatusType      *status

);

/*+ sc2dalib_mcecmdInit
*/
void sc2dalib_mcecmdInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *mymceCmd,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_mceonflycmdInit
*/
void sc2dalib_mceonflycmdInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *mymceCmd,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_mceerrRep - 
*/
void sc2dalib_mceerrRep
(
dasInfoStruct_t *myInfo,
dasCmdInfo_t    *mycmdInfo,
StatusType      *status
);

/*+ sc2dalib_mceerrnewRep - 
*/
void sc2dalib_mceerrnewRep
(
dasInfoStruct_t *myInfo,
dasCmdInfo_t    *mycmdInfo,
StatusType      *status
);

/*+ sc2dalib_mceerreachCard 
*/
void sc2dalib_mceerreachCard
(
dasInfoStruct_t *myInfo,
char            *mceCmd,
char            *cardName,
int             errNo,
StatusType      *status
);

/*+ sc2dalib_mcestatusInit
*/
void sc2dalib_mcestatusInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_msgprintSave
*/
void sc2dalib_msgprintSave
(
dasInfoStruct_t *myInfo, 
char            *string, 
char            *string2,
int             flag,
StatusType      *status
);

/*+ sc2dalib_openFiles
*/
void sc2dalib_openFiles
(
dasInfoStruct_t *myInfo,       
int             getData,
int             getBatch,
StatusType      *status
);

/*+ sc2dalib_pcicmdInit
*/
void sc2dalib_pcicmdInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *pciCmd,
char                  *memType,
long                  *memAddr,
long                  *memValue,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_pciblkInit
*/
void sc2dalib_pciblkInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *pciCmd,
char                  *memType,
long                  *startAddr,
long                  *blkSize,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_pcierrRep -
*/
void sc2dalib_pcierrRep
(
PCI_CMD *pciptr,
char * cmd,
StatusType *status
);

/*+ sc2dalib_pcistatusRep - 
*/
void sc2dalib_pcistatusRep
(
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
StatusType      *status
);

/*+   sc2dalib_pcisendCmd
*/
void  sc2dalib_pcisendCmd
(
SDSU_CONTEXT    *con,
char            *cmd,
dasInfoStruct_t *myInfo,
long            memAddr,
long            memValue,
char            *memType,
StatusType      *status
);

/*+   sc2dalib_pcisendblkCmd
*/
void  sc2dalib_pcisendblkCmd
(
SDSU_CONTEXT    *con,
char            *cmd,
dasInfoStruct_t *myInfo,
long            startAddr,
long            blkSize,
char            *memType,
StatusType      *status
);

/*+ sc2dalib_pixelmonInit
*/
void sc2dalib_pixelmonInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
char                  *dateTime,
ARRAYSET              *pixelset,
StatusType            *status
);

/*+ sc2dalib_readframeRate
*/
void sc2dalib_readframeRate
(
SDSU_CONTEXT          *con,
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceinxPtr,
long                   *rate,
StatusType            *status
);

/*+ sc2dalib_readmceVal
*/
void sc2dalib_readmceVal
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceinxPtr,
char                  *cmd,
int                   *val,
int                   howMany,
StatusType            *status
);

/*+ sc2dalib_readpixelMask
*/
void sc2dalib_readpixelMask
(
dasInfoStruct_t   *myInfo,    
int               *pixelMask,
SdsIdType         maskId,
StatusType        *status
);

/*+ sc2dalib_updateDebug 
*/
void sc2dalib_updateDebug
(
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
StatusType      *status    
);

/*+ sc2dalib_saveframeData 
*/
void sc2dalib_saveframeData
(
SDSU_CONTEXT    *con,       
dasInfoStruct_t *myInfo,
FILE            *fp,    
char            *byte,
int             datasize,
StatusType      *status
);

/*+ sc2dalib_savepixelData 
*/
void sc2dalib_savepixelData
(
SDSU_CONTEXT    *con,       
dasInfoStruct_t *myInfo,    
char            *byte,
int             datasize,
ARRAYSET        *pixelset,
StatusType      *status
);

/*+ sc2dalib_savecmdBuf - save command buffer into a file
*/
void sc2dalib_savecmdBuf
(
FILE *fpout,          
int   mxbuf,
char  buffer[mxbuf]  
);

/*+ sc2dalib_savemceReply -save MCE reply into a file
*/
void sc2dalib_savemceReply
(
FILE *fpout,    //store the cmdBuf sent to PCI, in byte order  
int   mxbuf,
char  *byte,    // pointer to reply Buffer 
uint32 count    // the count_Th reply 
);

/*+ sc2dalib_savemcerepData 
*/
void sc2dalib_savemcerepData
(
dasInfoStruct_t  *myInfo,
FILE             *fpout,   //store the replied data  
int              mxbuf,
char             *byte    // pointer to reply Buffer 
);

/*+ sc2dalib_savemcerepData 
*/
void sc2dalib_saveMceRepToHash
(
dasInfoStruct_t  *myInfo,
FILE             *fpout,      // file pointer where to store the replied data  
int              mxbuf,       // How many words came back from the MCE
char             *byte,       // Pointer to reply Buffer
char            *mycmd       // Pointer to command which caused the data
);

/*+ sc2dalib_sendCmd   
 */
void sc2dalib_sendCmd
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
dasCmdInfo_t          *mycmdInfo,
struct mcexml_struct  *mceInxpt,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_setmceVal
*/
void sc2dalib_setmceVal
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceinxPtr,
char                  *cmd,
StatusType            *status
);

/*+ sc2dalib_setseqInit   
 */
void sc2dalib_setseqInit
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_setseqInitENG
 */
void sc2dalib_setseqInitENG
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
char                  *whichCard,
StatusType            *status
);

/*+ sc2dalib_setseqInitRTSC
 */
void sc2dalib_setseqInitRTSC
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,
char                  *dateTime,
char                  *whichCard,
StatusType            *status
);

/*+ sc2dalib_fillInBeam   
 */
void sc2dalib_fillInBeam
(
long              bb_inbeam,         
long              fts_inbeam,         
long              pol_inbeam,
double            shutter,         
char             *inBeamArray,
StatusType       *status
);

/*+ sc2dalib_seqInit   
 */
void sc2dalib_seqInit
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
char                  *dateTime,
int                   **lookupTable,
StatusType            *status
);

/*+ sc2dalib_seqcallsc2headmanseqStart   
 */
void sc2dalib_seqcallsc2headmanseqStart
(
dasInfoStruct_t       *myInfo,
long                  startSeq,
long                  endSeq,
StatusType            *status
);

/*+ sc2dalib_seqcallsc2headmanmainHead   
 */
void sc2dalib_seqcallsc2headmanmainHead
(
dasInfoStruct_t       *myInfo,
long                  startSeq,
long                  endSeq,
StatusType            *status
);

/*+ sc2dalib_seqchkQL   
 */
void sc2dalib_seqchkQL
(
SDSU_CONTEXT          *con,
DRAMA_INNERMSG        *dhmsg,         
dasInfoStruct_t       *myInfo,
DitsReasonType        entReason,
int                   *localTable,
StatusType            *status
);

/*+ sc2dalib_seqsendQL   
 */
void sc2dalib_seqsendQL
(
SDSU_CONTEXT          *con,
DRAMA_INNERMSG        *dhmsg,         
dasInfoStruct_t       *myInfo,
StatusType            *status
);

/*+ sc2dalib_seqchkEnd   
 */
void sc2dalib_seqchkEnd
(
SDSU_CONTEXT          *con,
DRAMA_INNERMSG        *dramamsg,         
dasInfoStruct_t       *myInfo,
StatusType            *status
);

/*+ sc2dalib_seqEnd - 
*/
void sc2dalib_seqEnd
(
SDSU_CONTEXT    *con,           
dasInfoStruct_t *myInfo,       
int             *lookupTable,
StatusType      *status
);

/*+ sc2dalib_setsharedMem
*/
void sc2dalib_setsharedMem
(                  
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
int             whichShared,
int             sharedmemSize,  
StatusType      *status
);

/*+ sc2dalib_stepHeaterCurrent
*/
void sc2dalib_stepHeaterCurrent
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,   
StatusType            *status
);

/*+ sc2dalib_stepTESBias
*/
void sc2dalib_stepTESBias
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,   
StatusType            *status
);

/*+ sc2dalib_stopFrame
*/
void sc2dalib_stopFrame
(                  
SDSU_CONTEXT         *con,
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceinxPtr,
char                 *dateTime,
StatusType           *status
);

/*+ sc2dalib_trkheatInit
*/
void sc2dalib_trkheatInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_trkheatUpdate
 */
void sc2dalib_trkheatUpdate
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
int                  *data,
double               *slope,
StatusType           *status
);

/*+ sc2dalib_getStats
 */
void sc2dalib_getStats
(
 double theArray[],
 int index[],
 double *theMedian,
 int *theIndex,
 int n,
 double *avg,
 double *sigma
);

/*+ sc2dalib_trkheatpixelUpdate
 */
void sc2dalib_trkheatpixelUpdate
(
dasInfoStruct_t      *myInfo,
int                  heat,
int                  *data,
int                  *fdbkRef,
int                  numActivePixels,
int                  *activePixels,
double               *slope,
double               *heatVal,
int                  first,
int                  *pixelData,
StatusType           *status
);

/*+ sc2dalib_getsq2fbparaInit
*/
void sc2dalib_getsq2fbparaInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
char                  *dateTime,
ARRAYSET              *arrayset,
StatusType            *status
);

/*+ sc2dalib_trksq2fbInit
*/
void sc2dalib_trksq2fbInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
char                  *dateTime,
ARRAYSET              *arrayset,
long                  *option,
int                   *sq2fboptInit,
int                   *sq1biasOpt,
StatusType            *status
);

/*+ sc2dalib_trksq2FB
 */
void sc2dalib_trksq2FB
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                  *dateTime,
char                  *data,
ARRAYSET              *setup,
int                   *sq2fbOpt,
int                   option,
int                   skipDo,
StatusType            *status
);

/*+ sc2dalib_updatesq2fbVal
 */
void sc2dalib_updatesq2fbVal
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                  *dateTime,
char                  *data,
ARRAYSET              *setup,
int                   *sq2fbOpt,
StatusType            *status
);

/*+ sc2dalib_changesq2fbVal
 */
void sc2dalib_changesq2fbVal
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                 *dateTime,
int                  *valArray,
StatusType           *status
);

/*+ sc2dalib_sendarrayCmd
 */
void sc2dalib_sendarrayCmd
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                  *dateTime,
char                  *mceCmd,
int                   *array,
int                   howMany,
StatusType            *status
);

/*+ sc2dalib_trksq2fbupdateGain
 */
void sc2dalib_trksq2fbupdateGain
(
dasInfoStruct_t      *myInfo,
ARRAYSET              *setup,
int                   ch,
int                   delta,
double                *gain,
int                   *flag,
StatusType            *status
);

/*+ sc2dalib_readOPT
*/
void sc2dalib_readOPT
(
dasInfoStruct_t *myInfo,
int             *optArray,
int             dim,
char            *optFile,
StatusType      *status
);

/*+ sc2dalib_trksq2fbsavefbZG
*/
void sc2dalib_trksq2fbsavefbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dalib_savelckPoints1
*/
void sc2dalib_savelckPoints1
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
void       *paramarray,
int        flag,
StatusType *status
);

/*+ sc2dalib_trksq2fbsavesq2fbOPT
*/
void sc2dalib_trksq2fbsavesq2fbOPT
(
dasInfoStruct_t *myInfo,
int             *sq2fbOpt,
StatusType      *status
);

/*+ sc2dalib_variablesInit
*/
void sc2dalib_variablesInit
(
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
long            cols,
long            rows,
int             *pixelMask,
StatusType      *status
);

/*+ sc2dalib_versionInit
*/
void sc2dalib_versionInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *dateTime,
StatusType            *status
);

/*+ sc2dalib_versionInfo 
*/
void sc2dalib_versionInfo
(
SDSU_CONTEXT    *localcon,
dasInfoStruct_t *myInfo,
StatusType      *status
);

/*+ sc2dalib_callorphanHandler
*/
void sc2dalib_callorphanHandler
(
StatusType      *status
);

/*+ sc2dalib_sq1optptsInit
*/
void sc2dalib_sq1optptsInit
(
dasInfoStruct_t       *myInfo,
int                    *option,
StatusType            *status
);

/*+ sc2dalib_sq1optptsreadInput
*/
void sc2dalib_sq1optptsreadInput
(
FILE       *fp,
int        *bias,
int        *fdbk,
int        *refbias,
int        *reffdbk,
StatusType *status
);

/*+ sc2dalib_sq1biasfdbkPoints
*/
void sc2dalib_sq1biasfdbkPoints
(
char       *delimiters,
char        *name,
int        *paramarray,
int         row,
int         setVal,
StatusType *status
);

/*+ sc2dalib_sq1biasfdbk4notUsed
*/
void sc2dalib_sq1biasfdbk4notUsed
(
char        *name,
int        *paramarray,
int         row,
int         Val,
StatusType *status
);

/*+ sc2dalib_sq1optwritesq1bsq2fOut
*/
void sc2dalib_sq1optwritesq1bsq2fOut
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
int              *sq1bias,       //  bias for each SQ1 (given) 
int              *sq2feedback,    //SQ2 feedback at lock point for sq1bias(given) 
int              *sq1refbias,     // reference bias for each SQ1 (given) 
int              *sq2reffeedback, // SQ2 feedback at lock point for sq1refbias (given)
int               flag,            // = 1 allow sq2fb< 0;  =0 set sq2fb<0 pixel =BAD
StatusType *status
);

/*+ sc2dalib_chksinglePixel
*/
void sc2dalib_chksinglePixel
(
ARRAYSET *setup,
int      *sq1comp,       
int      *sq2comp,    
int      *sq1,     
int      *sq2, 
int      *flag,            
StatusType *status
);

/*+ sc2dalib_sq2fbeachcolOutlier
*/
void sc2dalib_sq2fbeachcolOutlier
(
ARRAYSET   *setup,
OPT_STRUCT *sq2fbHist,
int        *data,
int        col,
StatusType *status
);

/*+ sc2dalib_sq2fbhalfFlux
*/
void sc2dalib_sq2fbhalfFlux
(
ARRAYSET   *setup,
int        *data,        
int        which,
int        peakVal,
int        halfFlux,     
int        *changeInx,   
StatusType *status
);

/*+ sc2dalib_sq2fbOutlier
*/
void sc2dalib_sq2fbOutlier
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
char             *filename,
int              *sq2feedback,        
int              *sq2reffeedback,        
StatusType       *status
);

/*+ sc2dalib_sq1biassq2fdbkClipmean
*/
void sc2dalib_sq1biassq2fdbkClipmean
(
int  *sq1bias,
int  *sq1boptimal,
int  *sq2fdbk,
int  *sq2foptimal,
int  colNo,
int  rowNo,
StatusType *status
);

/*+ sc2dalib_sq1optAlgrm
*/
void sc2dalib_sq1optAlgrm
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
int              *sq1bcomp,
int              *sq2fcomp,
int              *sq1bias,         
int              *sq2feedback,    
int              *sq1refbias,     
int              *sq2reffeedback, 
int              *bqual,         
int              *fqual,         
int              *squal,         
double           *sq1scale,      
int              flag,           
StatusType *status
);

/*+ sc2dalib_writepixelInx
*/
void sc2dalib_writepixelInx
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
char             *filename,
void             *array,
int              flag,
StatusType       *status
);

/*+ writesq1biasfdbkOut
*/
void sc2dalib_writequalOut
(
FILE       *fd,
double     *sq1scale, 
int        *squal,
int        *fqual,
int        *bqual,
StatusType *status
);

/*+ writesq1biasfdbkOut
*/
void sc2dalib_writecompOut
(
FILE       *fd,
int        *sq1bcomp,
int        *sq2fcomp,
StatusType *status
);

/*+ sc2dalib_sc2dafindalter4Opt
*/
void sc2dalib_sc2dafindalter4Opt
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup,
int             *optArray, 
int             *orgArray, 
int             flag,  
int             which,           
StatusType      *status
);

/*+ sc2dalib_sq1optsave2File
*/
void sc2dalib_sq1optsave2File
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
int              *sq1bcomp,
int              *sq2fcomp,
StatusType *status
);

/*+ sc2dalib_sq1optsavesq1bComp
*/
void sc2dalib_sq1optsavesq1bComp
(
FILE       *fp,
char       *name,
int        *paramarray,
int        flag,
StatusType *status
);

/*+ sc2dalib_sq1optsavesq2fComp
*/
void sc2dalib_sq1optsavesq2fComp
(
FILE       *fp,
char       *name,
int        *paramarray,
int        flag,
StatusType *status
);

/*+ sc2dalib_writearray2setupFile
*/
void sc2dalib_writearray2setupFile
(
FILE       *fp,
char       *item,
void       *array,
int        flag,
StatusType *status
);

 
 
#endif

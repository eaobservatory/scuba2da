#ifndef HEADGEN____sc2dalibservo_h
#define HEADGEN____sc2dalibservo_h 
 
 
/*+ sc2daservoAlgrm
*/
void sc2daservoAlgrm
(
SDSU_CONTEXT     *con,
char            *byte,
dasInfoStruct_t *myInfo,
ARRAYSET       *setup,
int             *flagPtr,
int             fdbkInx,
char            *databuf,
StatusType      *status
);

/*+ sc2daservodataheadInit
*/
void sc2daservodataheadInit
(
ARRAYSET              *setup,
char                  *servoData,
StatusType            *status
);

/*+ sc2daservodataheadarrayVal
*/
void sc2daservodataheadarrayVal
(
char       *infoPtr,
void        *paramarray,
char        *name,
int         howMany,
int         flag,
StatusType *status
);

/*+ sc2daservodataheadarrayfloatVal
*/
void sc2daservodataheadarrayfloatVal
(
char       *infoPtr,
double     *paramarray,
char        *name,
int         howMany,
StatusType *status
);

/*+ sc2daservodatamemInit
*/
void sc2daservodatamemInit
(
SDSU_CONTEXT     *con,      
dasInfoStruct_t  *myInfo,    
ARRAYSET         *arrayset,
char             **dataBuf,
StatusType       *status
);

/*+ sc2daservobindataoutInit
*/
void sc2daservobindataoutInit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,  
int             colrow,   
StatusType      *status
);

/*+ sc2daservoFstGO   
 */
void sc2daservoFstGO
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
int                startSeq,
int                endSeq,
struct mcexml_struct  *mceInxpt,
char                  *dateTime,
StatusType            *status
);

/*+ sc2daservoEnd   
 */
void sc2daservoEnd
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
ARRAYSET             *setup,
char                  *dataBuf,
char                  *ssalckdataBuf,
int                   flag,
StatusType            *status
);

/*+ sc2daservoInit
*/
void sc2daservoInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
ARRAYSET              *arrayset,
char                  *dateTime,
char                  **dataBuf,
ARRAYSET              *ssalckset,
char                  **ssalckBuf,
StatusType            *status
);

/*+ sc2daservosetupInit
*/
void sc2daservosetupInit
(
dasInfoStruct_t       *myInfo,    
ARRAYSET              *setup,
StatusType            *status
);

/*+ sc2daservossalckInit
*/
void sc2daservossalckInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
ARRAYSET              *arrayset,
ARRAYSET              *ssalckset,
char                  **ssalckBuf,
StatusType            *status
);

/*+ sc2daservoinitValue
*/
void sc2daservoinitValue
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
ARRAYSET              *setup,
struct mcexml_struct  *mceinxPtr,
char                  *servoData,
char                  *dateTime,
StatusType            *status
);

/*+ sc2dareadmceparVal
*/
void sc2dareadmceparVal
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
ARRAYSET              *setup,
struct mcexml_struct  *mceinxPtr,
StatusType            *status
);

/*+ sc2daservoinitSSA
*/
void sc2daservoinitSSA
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
ARRAYSET              *setup,
struct mcexml_struct  *mceInxpt,
char                  *servoData,
char                  *dateTime,
StatusType            *status
);

/*+ sc2daservoinitSQ2
*/
void sc2daservoinitSQ2
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
ARRAYSET              *setup,
struct mcexml_struct  *mceInxpt,
char                  *servoData,
char                  *dateTime,
StatusType            *status
);

/*+ sc2daservoinitSQ1
*/
void sc2daservoinitSQ1
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
ARRAYSET              *setup,
struct mcexml_struct  *mceInxpt,
char                  *servoData,
char                  *dateTime,
StatusType            *status
);

/*+ sc2daservoprintMsg
*/
void sc2daservoprintMsg
(
ARRAYSET   *setup, 
FILE       *fp,
int        flag,
void       *paramarray,
char       *name,
StatusType *status
);

/*+ sc2daservoFunc
*/
void sc2daservoFunc
(
SDSU_CONTEXT         *con,
char                 *byte,
dasInfoStruct_t      *myInfo,
ARRAYSET             *setup,
struct mcexml_struct *mceInxpt,
int                  *isEnd,
char                 *dataBuf,
char                 *dateTime,
ARRAYSET             *ssalckset,
char                 *ssalckBuf,
StatusType           *status
);

/*+ sc2daservoFunc2
*/
void sc2daservoFunc2
(
SDSU_CONTEXT         *con,
char                 *byte,
dasInfoStruct_t      *myInfo,
ARRAYSET             *setup,
struct mcexml_struct *mceInxpt,
int                  *isEnd,
char                 *dataBuf,
char                 *dateTime,
ARRAYSET             *ssalckset,
char                 *ssalckBuf,
StatusType           *status
);

/*+ sc2daservossalckFunc
*/
void sc2daservossalckFunc
(
SDSU_CONTEXT         *con,
char                 *byte,
dasInfoStruct_t      *myInfo,
ARRAYSET             *setup,
ARRAYSET             *ssalckset,
struct mcexml_struct *mceInxpt,
char                 *dataBuf,
char                 *dateTime,
int                  *initFlag,
BIAS_PEAK            *peakbiasPtr,
int                   biasInx,
char                 *servoData,
StatusType           *status
);

/*+ sc2dasetBiasFdbk
*/
void sc2dasetBiasFdbk
(
SDSU_CONTEXT           *con,
dasInfoStruct_t        *myInfo,
ARRAYSET               *setup,
struct  mcexml_struct  *mceInxpt,
char                   *dateTime,
int                    fdbkInx,
StatusType             *status
);

/*+ sc2daservoresetFdbk
*/
void sc2daservoresetFdbk
(
SDSU_CONTEXT           *con,
dasInfoStruct_t        *myInfo,
ARRAYSET               *setup,
struct  mcexml_struct  *mceInxpt,
char                   *dateTime,
int                    *biasFlag,
StatusType             *status
);

/*+ sc2dasetbiasfdbkSSA
*/
void sc2dasetbiasfdbkSSA
(
SDSU_CONTEXT           *con,
dasInfoStruct_t        *myInfo,
ARRAYSET               *setup,
struct  mcexml_struct  *mceInxpt,
char                   *dateTime,
int                    fdbkInx,
StatusType             *status
);

/*+ sc2dacableCorrection
*/
void sc2dacableCorrection
(
SDSU_CONTEXT           *con,
dasInfoStruct_t        *myInfo,
ARRAYSET               *setup,
struct  mcexml_struct  *mceInxpt,
char                   *dateTime,
StatusType             *status
);

/*+ sc2dasetbiasfdbkSQ2
*/
void sc2dasetbiasfdbkSQ2
(
SDSU_CONTEXT           *con,
dasInfoStruct_t        *myInfo,
ARRAYSET               *setup,
struct  mcexml_struct  *mceInxpt,
char                   *dateTime,
int                    fdbkInx,
StatusType             *status
);

/*+ sc2dasetbiasfdbkSQ1
*/
void sc2dasetbiasfdbkSQ1
(
SDSU_CONTEXT           *con,
dasInfoStruct_t        *myInfo,
ARRAYSET               *setup,
struct  mcexml_struct  *mceInxpt,
char                   *dateTime,
int                    fdbkInx,
StatusType             *status
);

/*+ sc2dasetFDBK
*/
void sc2dasetFDBK
(
SDSU_CONTEXT           *con,
dasInfoStruct_t        *myInfo,
ARRAYSET               *setup,
struct  mcexml_struct  *mceInxpt,
char                   *dateTime,
int                    step,
StatusType             *status
);

/*+ sc2dasetstepdownFB
*/
void sc2dasetstepdownFB
(
SDSU_CONTEXT           *con,
dasInfoStruct_t        *myInfo,
ARRAYSET               *setup,
struct  mcexml_struct  *mceInxpt,
char                   *dateTime,
int                    manySteps,
StatusType             *status
);

/*+ sc2daservosaveData
*/
void sc2daservosaveData
(
SDSU_CONTEXT     *con,
dasInfoStruct_t  *myInfo,
ARRAYSET         *setup,
char             *byte,
int              *data,
int              fdbkInx,
StatusType       *status
);

/*+ sc2daservosavebinData
*/
void sc2daservosavebinData
(
dasInfoStruct_t  *myInfo,
ARRAYSET         *setup,
int              *data,
int              fdbkInx,
StatusType       *status
);

/*+ sc2daservosavebindata2File
*/
void sc2daservosavebindata2File
(
dasInfoStruct_t  *myInfo,
ARRAYSET         *setup,
StatusType       *status
);

/*+ sc2datransitFunc
*/
void sc2datransitFunc
(
SDSU_CONTEXT         *con,
char                 *byte,
dasInfoStruct_t      *myInfo,
ARRAYSET             *setup,
struct mcexml_struct *mceInxpt,
int                  *isEnd,
char                 *dataBuf,
char                 *dateTime,
StatusType           *status
);

/*+ sc2dawavegetData
*/
void sc2dawavegetData
(
char            *byte,
dasInfoStruct_t *myInfo,
ARRAYSET       *setup,
int             *flagPtr,
char            *databuf,
StatusType      *status
);

/*+ sc2datransitTES
*/
void sc2datransitTES
(
SDSU_CONTEXT           *con,
dasInfoStruct_t        *myInfo,
ARRAYSET               *setup,
struct  mcexml_struct  *mceInxpt,
char                   *dateTime,
StatusType             *status
);

/*+ sc2dafindTransit
*/
void sc2dafindTransit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *tesData,
StatusType      *status
);

 
 
#endif

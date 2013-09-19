#ifndef HEADGEN____sc2dalibsetup_h
#define HEADGEN____sc2dalibsetup_h 
 
 
/*+ sc2dalibsetup_applychFilter
*/
void sc2dalibsetup_applychFilter
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
double          *chData,
int              whichCh,
int              row,
int              whichRow,
int             dataLen,
int             *filtedchData,
StatusType      *status
);

/*+ sc2dalibsetup_linearConv
*/
void sc2dalibsetup_chlinearConv
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup,
double     xData,
double     *filteData,
int        row,
int        whichRow,
int        whichCh,
int        ith,
int        filtedBeg,
int        filtedEnd,
StatusType *status
);

/*+ sc2dalibsetup_sq2fbeachcolOutlier
*/
void sc2dalibsetup_sq2fbeachcolOutlier
(
ARRAYSET   *setup,
OPT_STRUCT *sq2fbHist,
int        *data,
int        col,
StatusType *status
);

/*+ sc2dalibsetup_sq2fbhalfFlux
*/
void sc2dalibsetup_sq2fbhalfFlux
(
ARRAYSET   *setup,
int        *data,        
int        which,
int        peakVal,
int        halfFlux,     
int        *changeInx,   
StatusType *status
);

/*+ sc2dalibsetup_sq2fbOutlier
*/
void sc2dalibsetup_sq2fbOutlier
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
char             *filename,
int              *sq2feedback,        
int              *sq2reffeedback,        
StatusType       *status
);

/*+ sc2dalibsetup_writepixelInx
*/
void sc2dalibsetup_writepixelInx
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
char             *filename,
void             *array,
int              flag,
StatusType       *status
);

/*+ sc2dalibsetup_savesq1biasoptPts
*/
void sc2dalibsetup_savesq1biasoptPts
(
FILE       *fp,
char       *name,
int        *paramarray,
int        flag,
StatusType *status
);

/*+ sc2dalibsetup_savesq2fboptPts
*/
void sc2dalibsetup_savesq2fboptPts
(
FILE       *fp,
char       *name,
int        *paramarray,
int        flag,
StatusType *status
);

/*+ sc2dalibsetup_readheaterSetup
*/
void sc2dalibsetup_readheaterSetup
(
dasInfoStruct_t *myInfo,
StatusType      *status
);

/*+ sc2dalibsetup_readheaterData
*/
void  sc2dalibsetup_readheaterData
(
dasInfoStruct_t  *myInfo,
int              *frameNo,
double           *chData,
int              *pixelData,
int              *heaterMask,
StatusType       *status
);

/*+ sc2dalibsetup_servoreadsetupWrap
*/
void sc2dalibsetup_servoreadsetupWrap
(
dasInfoStruct_t       *myInfo,
StatusType            *status
);

/*+ sc2dalibsetup_servoreadSetup
*/
void sc2dalibsetup_servoreadSetup
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dalibsetup_servochkSetup
*/
void sc2dalibsetup_servochkSetup
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
int             paramNo,
StatusType      *status
);

/*+ sc2dalibsetup_servosetCard
*/
void sc2dalibsetup_servosetCard
(
char       *delimiters,
char       *name,
int        *count,
StatusType *status
);

/*+ sc2dalibsetup_servosetparamVal
*/
void sc2dalibsetup_servosetparamVal
(
char       *delimiters,
int        *param,
int        *count,
StatusType *status
);

/*+ sc2dalibsetup_servosetparamarrayVal
*/
void sc2dalibsetup_servosetparamarrayVal
(
char       *delimiters,
int        *paramarray,
char        *name,
int         *index,
int         *count,
int         howMany,
StatusType *status
);

/*+ sc2dalibsetup_servosetparfloatVal
*/
void sc2dalibsetup_servosetparfloatVal
(
char       *delimiters,
double     *paramarray,
char        *name,
int         *index,
int         *count,
StatusType *status
);

/*+ sc2dalibsetup_servosetwhichRC
*/
void sc2dalibsetup_servosetwhichRC
(
dasInfoStruct_t *myInfo,
char       *delimiters,
ARRAYSET   *setup,
int        *count,
StatusType *status
);

/*+ sc2dalibsetup_servosetGain
*/
void sc2dalibsetup_servosetGain
(
char       *delimiters,
ARRAYSET   *setup,
int         *index,
int        *count,
StatusType *status
);

/*+ sc2dalibsetup_servosq1fdbkreadSetup
*/
void sc2dalibsetup_servosq1fdbkreadSetup
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dalibsetup_sq1biasfdbk4notUsed
*/
void sc2dalibsetup_sq1biasfdbk4notUsed
(
char        *name,
int        *paramarray,
int         row,
int         Val,
StatusType *status
);

/*+ sc2dalibsetup_sq1biasfdbkPoints
*/
void sc2dalibsetup_sq1biasfdbkPoints
(
char       *delimiters,
char        *name,
int        *paramarray,
int         row,
int         setVal,
StatusType *status
);

/*+ sc2dalibsetup_chksinglePixel
*/
void sc2dalibsetup_chksinglePixel
(
ARRAYSET *setup,
int      *sq1comp,       
int      *sq2comp,    
int      *sq1,     
int      *sq2, 
int      *flag,            
StatusType *status
);

/*+ sc2dalibsetup_sq1biassq2fdbkClipmean
*/
void sc2dalibsetup_sq1biassq2fdbkClipmean
(
int  *sq1bias,
int  *sq1boptimal,
int  *sq2fdbk,
int  *sq2foptimal,
int  colNo,
int  rowNo,
StatusType *status
);

/*+ sc2dalibsetup_writesq1biasfdbkOut
*/
void sc2dalibsetup_writesq1biasfdbkOut
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
char *filename,
int  *sq1bias,       // best bias for each SQ1 (given) 
int *sq2feedback,    //SQ2 feedback at lock point for sq1bias(given) 
int *sq1refbias,     // reference bias for each SQ1 (given) 
int *sq2reffeedback, // SQ2 feedback at lock point for sq1refbias (given)
int flag,            // = 1 allow sq2fb< 0;  flg=0 set  sq2fb<0 pixel =BAD
StatusType *status
);

/*+ sc2dalibsetup_sc2writequalOut
*/
void sc2dalibsetup_sc2writequalOut
(
FILE       *fd,
double     *sq1scale, 
int        *squal,
int        *fqual,
int        *bqual,
StatusType *status
);

/*+ sc2dalibsetup_findalter4Opt
*/
void sc2dalibsetup_findalter4Opt
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup,
int             *optArray, 
int             *orgArray, 
int             flag,  
int             which,           
StatusType      *status
);

/*+ sc2dalibsetup_writesq1biasfdbkOut
*/
void sc2dalibsetup_writecompOut
(
FILE       *fd,
int        *sq1bcomp,
int        *sq2fcomp,
StatusType *status
);

/*+ sc2dalibsetup_readsq1Data
*/
void sc2dalibsetup_readsq1Data
(
FILE       *fp,
int        *bias,
int        *fdbk,
int        *refbias,
int        *reffdbk,
StatusType      *status
);

/*+ sc2dalibsetup_whichRC
*/
void sc2dalibsetup_whichRC
(
dasInfoStruct_t  *myInfo,    
char             *cardName,  
StatusType       *status
);

/*+ sc2dalibsetup_readFile*/
void sc2dalibsetup_readFile
(
dasInfoStruct_t *myInfo,
char            **filememPtr, 
int             *fileLen,
StatusType      *status
);

/*+ sc2dalibsetup_servodataoutInit
*/
void sc2dalibsetup_servodataoutInit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,  
int             colrow,   
StatusType      *status
);

/*+ sc2dalibsetup_heaterGen
*/
void sc2dalibsetup_heaterGen
(
ARRAYSET        *setup,
int             funcFlag,
StatusType      *status
);

/*+ sc2dalibsetup_pixelTransit
*/
void sc2dalibsetup_pixelTransit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *data,
int             whichRow,
int             heatLvl,
int            flag,
FILE           *fp,
StatusType      *status
);

/*+ sc2dalibsetup_heater2Power
*/
void sc2dalibsetup_heater2Power
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
double          *pwer,
StatusType      *status
);

/*+ sc2dalibsetup_heater2Power
*/
void sc2dalibsetup_sq1fb2Current
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
double          *current,
StatusType      *status
);

/*+ sc2dalibsetup_singExt
*/
void sc2dalibsetup_signExt
(
dasInfoStruct_t *myInfo,
int             *input,
int             signBit,
int             dataMask,
int             signextMask,
StatusType      *status
);

/*+ sc2dalibsetup_correlateCoeff
*/
void sc2dalibsetup_correlateCoeff
(
int np,               /* number of points (given) */
double *x,           /* X data (given) */
double *y,           /* Y data (given) */
double *coeff,       /* correlation coefficient (returned) */
StatusType *status   /* global status (given and returned) */
);

/*+ sc2dalibsetup_transitResult
*/
void sc2dalibsetup_transitResult
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
FILE            *fp,
StatusType      *status
);

/*+ sc2dalibsetup_saveTransit
*/
void sc2dalibsetup_saveTransit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
int             row,
int             whichRow,
StatusType      *status
);

 
 
#endif

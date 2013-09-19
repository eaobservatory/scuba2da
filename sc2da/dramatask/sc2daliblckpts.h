#ifndef HEADGEN____sc2daliblckpts_h
#define HEADGEN____sc2daliblckpts_h 
 
 
/*+ sc2daapplyFilter
*/
void sc2daapplyFilter
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int              whichCh,
int             *biasData,
int             *filtedbiasData,
StatusType      *status
);

/*+ sc2dachkDeriv
*/
void sc2dachkDeriv
(
ARRAYSET   *setup, 
void       *chdata, 
double     *derivt, 
int         ch,
int         flag,
StatusType *status
);

/*+ sc2dafilterFunc
*/
void sc2dafilterFunc
(
dasInfoStruct_t    *myInfo,    
double             fc,
double             sf,
double             *impulseRep, 
int                filterOrder, 
int                useotherFlag,
StatusType         *status
);

/*+ sc2dafindcableadjVal
*/
void sc2dafindcableadjVal
(
ARRAYSET   *setup, 
int        biasInx,
StatusType *status
);

/*+ sc2dafindlockPts
*/
void sc2dafindlockPts
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fplck,
char            *servoFile,
StatusType      *status
);

/*+ sc2dafindinitlckPts
*/
void sc2dafindinitlckPts
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup, 
char       *databuf,
StatusType *status
);

/*+ sc2dafindlckptsGain
*/
void sc2dafindlckptsGain
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fplck,
char            *servoFile,
int             *pixelMask,
int              flag,
StatusType      *status
);

/*+ sc2dafindflatBottom
*/
void sc2dafindflatBottom
(
SERVO_INFO    *myservoInfo,
ARRAYSET       *setup,
char           *databuf,
StatusType     *status
);

/*+ sc2dafindflatDeriv
*/
void sc2dafindflatDeriv
(
ARRAYSET   *setup, 
int        *chdata, 
int        fbStart,
int        *derivt, 
double     *deriVal,
double     *thdVal,
int        ch,
StatusType *status
);

/*+ sc2dafindframelckptsGain
*/
void sc2dafindframelckptsGain
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int             *pixelMask,
StatusType      *status
);

/*+ sc2dafindpixellckptsGain2
*/
void sc2dafindpixellckptsGain2
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int        row,
int        col, 
BIAS_PEAK *peakInfo,
int       *pixelMask,
StatusType *status
);

/*+ sc2dafindsq2fbGain
*/
void sc2dafindsq2fbGain
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int             *pixelMask,
StatusType      *status
);

/*+ sc2dafindpixelfluxPeriod
*/
void sc2dafindpixelfluxPeriod
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup,
int             pixel,
StatusType      *status
);

/*+ sc2dafindmaxMin
*/
void sc2dafindmaxMin
(
double     *chdata, 
BIAS_PEAK  *peakInfo,
ARRAYSET   *setup,
int        ch,
StatusType *status
);

/*+ sc2dafindmeanVal
*/
void sc2dafindmeanVal
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fplck,
char            *servoFile,
StatusType      *status
);

/*+ sc2dafindminMax
*/
void sc2dafindminMax
(
double     *chdata, 
BIAS_PEAK  *peakInfo,
ARRAYSET   *setup,
int        ch,
StatusType *status
);

/*+ sc2dafindmaxModul
*/
void sc2dafindmaxModul
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup, 
char       *databuf,
char       *fileName,
StatusType *status
);

/*+ sc2dafindrefmaxModul
*/
void sc2dafindrefmaxModul
(
dasInfoStruct_t *myInfo,
ARRAYSET       *setup,
char           *databuf,
StatusType     *status
);

/*+ sc2dafindrefPoint
*/
void sc2dafindrefPoint
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup,
int        ch, 
int        *chdatarefbiasPtr,
StatusType *status
);

/*+ sc2daleastsquareFit
*/
void sc2daleastsquareFit
(
ARRAYSET   *setup, 
char       *databuf,
StatusType *status
);

/*+ sc2dalinearConv
*/
void sc2dalinearConv
(
ARRAYSET   *setup,
double     xData,
double     *filteData,
StatusType *status
);

/*+ sc2dalookeachCh
*/
void sc2dalookeachCh
(
ARRAYSET   *setup, 
int     *biasData,
BIAS_PEAK  *peakInfo,
int         whichCh,
StatusType *status
);

/*+ sc2dalookabsMax
*/
void sc2dalookabsMax
(
void     *chdata, 
ARRAYSET   *setup,
int        ch,
int        *thdFlag,
int         flag,
BIAS_PEAK  *peakInfo,
StatusType *status
);

/*+ sc2dalookMax
*/
void sc2dalookMax
(
double     *chdata, 
BIAS_PEAK  *peakInfo,
ARRAYSET   *setup,
int        ch,
int        *dataInx,
int         order,
StatusType *status
);

/*+ sc2dalookmeaneachCh
*/
void sc2dalookmeaneachCh
(
ARRAYSET   *setup, 
int     *biasData,
BIAS_PEAK  *peakInfo,
int         whichCh,
StatusType *status
);

/*+ sc2dalookMin
*/
void sc2dalookMin
(
double     *chdata, 
BIAS_PEAK  *peakInfo,
ARRAYSET   *setup,
int        ch,
int        *dataInx,
int        order,
StatusType *status
);

/*+ sc2dareadbinaryfile2Mem*/
void sc2dareadbinaryfile2Mem
(
dasInfoStruct_t *myInfo,      
FILE            *fp, 
char            **memPtr,       
int             *fileLen,             
StatusType      *status
);

/*+ sc2dasaveservoData
*/
void sc2dasaveservoData
(
ARRAYSET   *setup, 
char       *databuf,
char       *fileName,
StatusType *status
);

/*+ sc2dasavelockPts
*/
void sc2dasavelockPts
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
);

/*+ sc2dasaveframelockPts
*/
void sc2dasaveframelockPts
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2daservocopydataBuf
*/
void sc2daservocopydataBuf
(
dasInfoStruct_t       *myInfo,    
ARRAYSET              *setup,
char                  *dataBuf,
StatusType            *status
);

/*+ sc2daservogetfiltedData
*/
void sc2daservogetfiltedData
(
dasInfoStruct_t       *myInfo,    
ARRAYSET              *setup,
FILE                  *fp,
StatusType            *status
);

/*+ sc2daservolckPoints
*/
void sc2daservolckPoints
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
int        *paramarray,
int         flag,
StatusType *status
);

/*+ sc2dasavessalckDivbyn
*/
void sc2dasavessalckDivbyn
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
int        *paramarray,
int         flag,
StatusType *status
);

/*+ sc2daservofindfluxPeriod
*/
void sc2daservofindfluxPeriod
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int        ch, 
int        *chdatamaxbiasPtr,
StatusType *status
);

/*+ sc2daservosq1lckPoints
*/
void sc2daservosq1lckPoints
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
int        *paramarray,
int         flag,
StatusType *status
);

/*+ sc2daservolckgainPoints
*/
void sc2daservolckgainPoints
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
double     *paramarray,
int         flag,
StatusType *status
);

/*+ sc2daservolckPoints1
*/
void sc2daservolckPoints1
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
int        *paramarray,
StatusType *status
);

/*+ sc2daservolckgainPoints1
*/
void sc2daservolckgainPoints1
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
double     *paramarray,
StatusType *status
);

/*+ sc2daservopixellckpts4MCE
*/
void sc2daservopixellckpts4MCE
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
FILE       *fp,
int        flag,
char       *name,
void       *paramarray,
StatusType *status
);

/*+ sc2daservopixellckPoints
*/
void sc2daservopixellckPoints
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
FILE       *fp,
int        flag,
void       *paramarray,
StatusType *status
);

/*+ sc2dasq1maxrefPoints
*/
void sc2dasq1maxrefPoints
(
ARRAYSET   *setup,
FILE       *fp,
char        *name,
int        *paramarray,
StatusType *status
);

/*+ sc2daservosaveFilted
*/
void sc2daservosaveFilted
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int             *biasData,
int             *filtedbiasData,
int              bias,
FILE             *fp,
StatusType      *status
);

/*+ sc2dasave4nextStep
*/
void sc2dasave4nextStep
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
);

/*+ sc2dasave4cableStep
*/
void sc2dasave4cableStep
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
);

/*+ sc2dasave4ssaStep
*/
void sc2dasave4ssaStep
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
);

/*+ sc2dasave4sq2Step
*/
void sc2dasave4sq2Step
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
);

/*+ sc2dasave4sq1Step
*/
void sc2dasave4sq1Step
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
);

/*+ sc2dassabiasoffset4batchfile
*/
void sc2dassabiasoffset4batchfile
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasavecableadjInit
*/
void sc2dasavecableadjInit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dassabias4mcelockSSA
*/
void sc2dassabias4mcelockSSA
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasavessaoutmidVal
*/
void sc2dasavessaoutmidVal
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasavessabiasOffset
*/
void sc2dasavessabiasOffset
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasavessabiasfbZG
*/
void sc2dasavessabiasfbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+  sc2dassafbZGAppd4SQ2lock
*/
void  sc2dassafbZGAppd4SQ2lock
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasq2bias4batchfile
*/
void sc2dasq2bias4batchfile
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasavesq2biasfbZG
*/
void sc2dasavesq2biasfbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasq2bias4ssaLock
*/
void sc2dasq2bias4ssaLock
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dassafb4batchfile
*/
void sc2dassafb4batchfile
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dassafbsq2bias4mcelockSQ2
*/
void sc2dassafbsq2bias4mcelockSQ2
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasavesq2biasFB
*/
void sc2dasavesq2biasFB
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasavesq2Flux
*/
void sc2dasavesq2Flux
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasavesq2open4pfbZG
*/
void sc2dasavesq2open4pfbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasavesq2openbiasfbZG
*/
void sc2dasavesq2openbiasfbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2dasq2openZ4SQ1open
*/
void sc2dasq2openZ4SQ1open
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
);

/*+ sc2daservofindslopGain
*/
void sc2daservofindslopGain
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup,
int        ch, 
int        *chdatamaxbiasPtr,
StatusType *status
);

/*+ sc2damath_linfit - straight line fit */

void sc2damath_linfit
(
int np,               /* number of points (given) */
double x[],           /* X data (given) */
double y[],           /* Y data (given) */
double wt[],          /* weights (given) */
double *grad,         /* slope (returned) */
double *cons,         /* offset (returned) */
int *status           /* global status (given and returned) */
);

/*+ sc2dafindpixeltwoPeaks
*/
void sc2dafindpixeltwoPeaks
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int        pixel,
int         flag,
StatusType *status
);

/*+ sc2dafindpixelmaxPeak
*/
void sc2dafindpixelmaxPeak
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int         pixel, 
int         inx,
int         startFB,
int         searchPerd,
StatusType *status
);

/*+ sc2dafindpixelminPeak
*/
void sc2dafindpixelminPeak
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int         pixel, 
int         inx,
int         startFB,
int         searchPerd,
StatusType *status
);

/*+ sc2daliblckpts_applycolFilter
*/
void sc2daliblckpts_applycolFilter
(
double   *impulse, 
int      filterOrder,
int      ydim,
int      *chData,
int      *filtedData,
double   *convolPtr, 
StatusType *status
);

/*+ sc2daliblckpts_applyFilter
*/
void sc2daliblckpts_applyFilter
(
double   *impulse, 
int      filterOrder,
int      ydim,
int      whichCh,
int      totalCh,
int      *orgData,
int      *filtedData,
double   *convolPtr, 
StatusType *status
);

/*+ sc2daliblckpts_filterFunc
*/
void sc2daliblckpts_filterFunc
(
double   *impulse, 
int      filterOrder, 
int      useFlag,
StatusType *status
);

/*+ sc2daliblckpts_linearConv
*/
void sc2daliblckpts_linearConv
(
double    *impulsePtr,
double    *convolPtr,
int       filterOrder,
double    xData,
double    *filtedData,
StatusType *status
);

/*+ sc2daliblckpts_finddataSlope
*/
void sc2daliblckpts_finddataSlope
(
int    *data,
int    start,
double *slope,
StatusType *status
);

/*+ sc2daliblckpts_findFiltered
*/
void sc2daliblckpts_findFiltered
(
ARRAYSET   *setup,
int        *chData,
int        *chfiltedData,
StatusType *status
);

/*+ sc2daliblckpts_chkPeaks
*/                                                      
void sc2daliblckpts_chkPeaks
(
int            *data,
int             ydim,
int             *peakNo,
int             servoFlag,
int             ch,
int             *status
);

/*+ sc2daliblckpts_findsearchPt
*/                                                      
void sc2daliblckpts_findsearchPt
(
int            *chData,
int             ydim,
int             *fbStart,
int             *value,
int             *p2pVal,
int             servoFlag,
int             ch,
StatusType      *status
);

/*+ sc2daliblckpts_findfluxPeriod
*/
void sc2daliblckpts_findfluxPeriod
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup,
int             ch, 
int             *data,
int             searchVal,
int             start,
StatusType      *status
);

/*+ sc2daliblckpts_findframelckptsGain
*/
void sc2daliblckpts_findframelckptsGain
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int             *pixelMask,
StatusType      *status
);

/*+ sc2daliblckpts_findinitlckPts
*/
void sc2daliblckpts_findinitlckPts
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup, 
char       *databuf,
StatusType *status
);

/*+ sc2daliblckpts_findlckptsGain
*/
void sc2daliblckpts_findlckptsGain
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fplck,
char            *servoFile,
int             *pixelMask,
int             flag,
StatusType      *status
);

/*+ sc2daliblckpts_findlockPts
*/
void sc2daliblckpts_findlockPts
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fplck,
char            *servoFile,
StatusType      *status
);

/*+ sc2daliblckpts_findmidVal
*/
void sc2daliblckpts_findmidVal
(
int    *data,
int    *midVal,
int    *midInx,
int    start,
int    stop,
StatusType *status
);

/*+ sc2daliblckpts_findmidvalslpFlux
*/
void sc2daliblckpts_findmidvalslpFlux
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int        ch, 
int        *chData,
StatusType *status
);

/*+ sc2daliblckpts_findslpGain
*/
void sc2daliblckpts_findslpGain
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup,
int        ch, 
int        *chdatamaxbiasPtr,
StatusType *status
);

/*+ sc2daliblckpts_findsq2fbGain
*/
void sc2daliblckpts_findsq2fbGain
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int             *pixelMask,
StatusType      *status
);

/*+ sc2daliblckpts_findpixellckptsGain
*/
void sc2daliblckpts_findpixellckptsGain
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int        row,
int        col, 
BIAS_PEAK  *peakInfo,
int        *pixelMask,
StatusType *status
);

/*+ sc2daliblckpts_readBinary
*/
void sc2daliblckpts_readBinary
(
dasInfoStruct_t *myInfo,
char            **frameData,
StatusType      *status
);

/*+ sc2daliblckpts_savefilteredData
*/
void sc2daliblckpts_savefilteredData
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
StatusType      *status
);

/*+ sc2daliblckpts_findMax
*/
void sc2daliblckpts_findMax
(
int   *data,
int   ydim, 
int   *maxVal, 
int   *maxInx,
int   *startOffset,
int   searchPerd,
int   ith,
StatusType *status
);

/*+ sc2daliblckpts_findMin
*/
void sc2daliblckpts_findMin
(
int   *data,
int   ydim, 
int   *minVal,
int   *minInx,
int   *startOffset,
int   searchPerd,
int   ith,
StatusType *status
);

/*+ sc2daliblckpts_getThreshold
*/
void sc2daliblckpts_getThreshold
(
double *thredVal,
int    *peakVal, 
int    servoFlag,
int    *status
);

/*+ sc2daliblckpts_findMaxMin
*/
void sc2daliblckpts_findMaxMin
(
dasInfoStruct_t *myInfo,    
int   *data,
int   ydim, 
int   *maxVal, 
int   *maxInx,
int   *minVal,
int   *minInx,
int   *fbstartStp,
int   startOffset,
int   servoFlag,
StatusType *status
);

/*+ sc2daliblckpts_findmaxWanted
*/                                                      
void sc2daliblckpts_findmaxWanted
(
FILE            *fp,
int             *p2pData,
int             whichCh,
int             ydim,
int             *bias,
int             *maxp2pVal,
int             *maxp2pInx,
double          thresHold,
int             servoFlag,
StatusType      *status
);

/*+ sc2daliblckpts_findmaxModul
*/
void sc2daliblckpts_findmaxModul
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup, 
char       *databuf,
char       *fileName,
StatusType *status
);

 
 
#endif

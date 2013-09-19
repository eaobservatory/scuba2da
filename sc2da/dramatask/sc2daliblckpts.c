/**
 * \file sc2daliblckpts.c
 *
 * \brief 
 *  collection of sub-functions for sc2da.c
 *  find lock points  
 *        
 * \author Xiaofeng Gao, UKATC (xg@roe.ac.uk)
 *
 * \version 
 *
 * \date 
 *
 * Copyright (c) 2005.  UK Astronomy Technology Centre (UK ATC), 
 * An establishment of the Particle Physics and Astronomy Research Council
 * PPARC).
 *
 * Web: www.roe.ac.uk
 *
 *
 *  $Log: sc2daliblckpts.c,v $
 *  Revision 1.39  2012/06/05 23:23:40  cwalther
 *  Cleared the heatTackFailed flag somewhere where I knew would execute during setup
 *
 *  Revision 1.38  2011/09/02 19:16:20  bgorges
 *  Standardized entry checks on all fuctions/actions that should quietly return if entry status is bad.
 *
 *  Revision 1.37  2011/06/02 21:31:22  bgorges
 *  Checking all of fopen and fopen64 to make sure a file is actually opened.
 *
 *  Revision 1.36  2011/04/22 00:25:47  cwalther
 *  I put a considerable number of comments in this code and it has the changes made to normalize thresholds when the number of samples changes
 *
 *  Revision 1.35  2011/03/01 23:47:28  cwalther
 *  Added test for maximum gaini, set to zero if greater than maximum
 *
 *  Revision 1.34  2010/12/04 00:29:36  cwalther
 *  Removing changing select_clk each time as it was crashing the MCE
 *
 *  Revision 1.33  2010/10/21 02:01:28  cwalther
 *  Removed optimal-sq1bias-sq2fb from CONFIG_HARD it is in scratch now
 *
 *  Revision 1.32  2010/10/12 22:23:31  cwalther
 *  Changed MCE command zero to flx_quanta to match UBC documents
 *
 *  Revision 1.31  2010/10/07 23:32:15  cwalther
 *  Trying to make the log file more useable
 *
 *  Revision 1.30  2010/09/30 20:32:46  cwalther
 *   Marked locations where i_clamp_val changes are needed
 *
 *  Revision 1.29  2010/09/17 19:45:16  cwalther
 *  Found that I should be looking for sq2lock-bias in scratch
 *
 *  Revision 1.28  2010/08/03 22:16:10  cwalther
 *  Changes to make sc2_setup work with SC2SCRATCH
 *
 *  Revision 1.1.1.1  2007/05/16 08:26:59  dkelly
 *  first insertion
 *  
 *
 */


// global parameter
SERVO_INFO  servoInfo;

// no of data used for average in derivative calculation and 
// check start point from the beginning of the channel data array
// change AVER_NO to 5 after seeing sq1plot  27/11/2006 x. Gao
#define AVER_NO    10
#define CHKSTART   20

// the threshold value( need to test)  for derivative if 
// there is a modulation
//#define DER_THDSSA    20
//#define DER_THDSSA    30    // for SG
#define DER_THDSSA      100   // use 20 data LSF  
#define DER_THDSQ2      10
#define DER_THDSQ1      10
#define DER_THDSQ1OPEN  10
#define DER_THDSQ2OPEN4P  10

// if peaktopeak < PEAK_THD, no modulation
#define PEAK_THDSSA       416  // This is multiplied by the number of samples, if DIVBY_N=1 in mceinitialise-old 
#define PEAK_THDSQ2      2000  
#define PEAK_THDSQ1      1000  
#define PEAK_THDSQ1OPEN   1000  
#define PEAK_THDSQ2OPEN4P 1000

// the thread hold value( need to test)  for flat 
#define FLAT_THDSSA      2   // change for using slope on both side
#define FLAT_THDSQ2      3
#define FLAT_THDSQ1      3
#define FLAT_THDSQ1OPEN  3

// average period for seraching flat
#define FLAT_AVER_NO 5 

// debug flag
#define DEBFLAT  
#define DEBFLUX  
#define DEBFLUX1  // print fluxperiod

// see each channel
//#define DEBSEREACH     

// which ch to use 
#define TEST_CH    0
#define TEST_ROW   1
#define TEST_PIXEL   TEST_ROW*COL_NUM+TEST_CH  //42
//#define PRINT_CH_FLAT  // print msg for flatDeriv 

// log midVal, cabadjTHd adjVal to logfile
//#define DEBMEANVAL 1

//#define PRINT_CH_DATA 
// if use findlockpoints1 to test, then define
//#define FINDLOCKPOINT_PROG  1
 
// use DEBSER3 to see the search for max-point
//#define DEBSER3(x...)  //printf(x)

// debug to see the search Start-stop of Midpoint
//#define DEB_MIDPOINT     

static char msg[FILE_LEN*2];


//=================sc2daa ===============================
//=======================================================

/**
 * \fn void sc2daapplyFilter(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     int whichCh, int  *biasData, int *filtedbiasData, StatusType *status)
 *
 * \brief function
 *  apply filter  to biasData and save result to filtedbiasdata 
 * 
 * \param myInfo     dasInfo structure pointer  
 * \param  setup    ARRAYSET structure pointer  
 * \param  whichCh  which channel is looked at
 * \param  biasData data pointer for each bias level 
 * \param  filtedbiasData filteddata pointer for each bias level 
 * \param  status   StatusType pointer  
 *
 * note: (N-1)/2 delay
 * 
 */
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
)
{
  int        i,feedbk;
  int        M=(FILTER_ORDER-1)/2;
  int        filtedStart, filtedEnd;
  double     xData,conData;

  if (!StatusOkP(status)) return;

  // http://www.cppreference.com/stdstring/memset.html
  // initial all val =0
  memset( setup->convolPtr, '\0', sizeof(double)*FILTER_ORDER );       

  filtedStart=2*M;
  filtedEnd= setup->fdbkNo + 2*M;

  for (i=0; i< (setup->fdbkNo+2*M); i++)
  {
     if( i < M ) /* 0 M-1 is additional points for convol */
       xData=(double)biasData[whichCh+COL_NUM*0]; 
     else if( i < (setup->fdbkNo+M) )
       xData=(double)biasData[whichCh+COL_NUM*(i-M)] ;     
     else  /* fdbkNo+M is additional points for convol*/
       xData=(double)biasData[whichCh+COL_NUM*(setup->fdbkNo-1)] ;
   
     sc2dalinearConv(setup,xData, &conData, status);
  
     // take (N-1)/2 delay off plus the pending ones
     if( i>=filtedStart &&  i<filtedEnd )
     {
       feedbk=whichCh+ COL_NUM*(i-filtedStart);
       filtedbiasData[feedbk]=(int)conData;
     }
  }
}


//========================= sc2dac  =========================
//===========================================================
/**
 * \fn void sc2dachkDeriv(ARRAYSET *setup,  void *chdata, double *derivt,
 *     int ch, int flag, StatusType *status)
 *
 * \brief function
 *  Find the first place in the data where the slope is greater than
 *  the threshold, for the current servo we are running.
 *  Return that slope in parameter derivt.  The sign of that slope
 *  Tells us whether the modulation goes up or down at the beginning
 *
 *         given and return ==>G&R        
 * 
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  chdata   void pointer for channel data 
 * \param  derivt  double pointer for derivative/ch, G&R
 * \param  ch        int
 * \param  flag      int  0: int chdata; 1: double chdata 
 * \param  status   StatusType pointer.  G&R
 *
 *  default AVER_NO=5; CHKSTART=5
 *  derivative=Average(y14:y10)-Average(y9:y5)
 */
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
)
{
  double  data, slope,dclevel;
  double  thredVal=DER_THDSSA;
  int     i, j, chkNo, start,ssafitWind=20;
  double  yData[ssafitWind],xData[ssafitWind],weightData[ssafitWind];
  int      *intArray;
  double   *floatArray;

  intArray=NULL;
  floatArray=NULL;

  if (!StatusOkP(status)) return;

  // Determine the threshold to use base on servo type

  if (setup->servo==SSARAMP || setup->servo==SSALOCK)
  {
    thredVal=DER_THDSSA;
  }
  else if (setup->servo==SQ2SERVO || setup->servo==SQ2LOCK)
    thredVal=DER_THDSQ2;
  else if (setup->servo==SQ1SERVO || setup->servo ==SQ1LOCK)
    thredVal=DER_THDSQ1;
  else if (setup->servo==SQ1OPEN )
    thredVal=DER_THDSQ1OPEN;
  else if (setup->servo==SQ2OPEN4P )
    thredVal=DER_THDSQ2OPEN4P;

  // the CHKSTART is used to avoid the first few data as when I look at 
  // data plot, it is a bit weird at beginning 
  if ( flag==1)
    floatArray=(double*)chdata;
  else
    intArray=(int*)chdata;

  chkNo=(setup->fdbkNo*2/5)/ssafitWind;
  start=CHKSTART;
  
  // search through shifted window
  for (j=0;j<chkNo;j++)
  {

    // Fill array data with ssafitWind number of data points
    for (i=0; i<ssafitWind; i++)
    {
      if ( flag==1)
        data=floatArray[start+i];
      else
        data=(double)intArray[start+i];

      xData[i]=i+1;
      weightData[i]=1;
      yData[i]=data;
    }
        
    // Do a linear fit on those data, slope is the result
    sc2damath_linfit(ssafitWind, xData, yData,weightData, &slope,&dclevel,status); 

    // Diagnostic message out
    if (setup->servo==SQ2OPEN4P )
    {
      sprintf(msg,"ch_%d  deriv=%f  thredVal=%f ",ch,slope,thredVal);
      MsgOut(status,msg);
    }

    // Is the absolute value of this slope greater than the threshold?
    if ( fabs(slope)> thredVal )
    {
      *derivt=slope;
      if (setup->servo==SQ2OPEN4P )
	{
	  printf("deriv(%f)>THD at shiftWind(from %d)",slope,(CHKSTART + j*ssafitWind) );
	}
      return;
    }
    // move the data pointer
    start +=ssafitWind;
  }
}


//========================= sc2daf  =========================
//===========================================================

/**
 * \fn void sc2dafilterFunc(dasInfoStruct_t *myInfo, double fc, double sf,
*     double *impulseRep, int filterOrder, int useotherFlag,StatusType *status)
 *
 * \brief function
 *  generate filter impulse h(n) 
 * 
 * \param myInfo     dasInfo structure pointer  
 * \param fc          double lowpass cut-off freq 
 * \param sf          double sampling freq 
 * \param impulseRep  double pointer 
 * \param filterOrder int, filter order
 * \param useotherFlag int,  1(lp):2(hp) use coeffienct from thrid party
 * \param  status     StatusType pointer  
 *
 * lowpass filter, order =N
 * PI=3.14159; 
 * h(n)=sin 2*PI*(FC/FS)*(n-(N-1)/2) /PI* (n-(N-1)/2) 
 *
 *  from sq1 data data, samplingFreq=301Hz, lowpassFre=10Hz,
 *  Lowpass 31 
 *
 */

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
)
{
  int        i;
  int        M;
  double     LP1;
  double     tmp[filterOrder];

  if (!StatusOkP(status)) return;

  if (useotherFlag==1 || useotherFlag==2)
  {
    // scopeFIR windownedSinc, FC 10Hz, SF301, order=31, very good
    double lp[]=
    {
    0.000187941840750392, 0.004194990908648022, 0.008622369643936374,
    0.013382100888852951, 0.018373423505132734, 0.023485541305617180,
    0.028600744412886091, 0.033597804161601595, 0.038355532236184611,
    0.042756388064969252, 0.046690015839649227, 0.050056594013285688,
    0.052769885720116753, 0.054759888075295048, 0.055974991428704754,
    0.056383575908738812, 0.055974991428704754, 0.054759888075295048,
    0.052769885720116753, 0.050056594013285688, 0.046690015839649227,
    0.042756388064969252, 0.038355532236184611, 0.033597804161601595,
    0.028600744412886091, 0.023485541305617180, 0.018373423505132734,
    0.013382100888852951, 0.008622369643936374, 0.004194990908648022,
    0.000187941840750392
    };
 
    // scopeFIR windownedSinc, SF301, order=31, Stopband=10Hz, passband=20Hz 
    // passband ripple=1dB, stopband attenuation=20dB
    double hp[]=
    {
     0.031671143951709085, 0.009867800403131675, 0.009084418761796206,
     0.006445151203763662, 0.001649924649358013,-0.005408521197426909,
    -0.014589662963512530,-0.025574807725622557,-0.037897407775362338,
    -0.050930540610479393,-0.063903781650189229,-0.075983538946147511,
    -0.086368413213330428,-0.094357468934419167,-0.099395015878274540,
     0.898883128394097340,-0.099395015878274540,-0.094357468934419167,
    -0.086368413213330428,-0.075983538946147511,-0.063903781650189229,
    -0.050930540610479393,-0.037897407775362338,-0.025574807725622557,
    -0.014589662963512530,-0.005408521197426909, 0.001649924649358013,
     0.006445151203763662, 0.009084418761796206, 0.009867800403131675,
     0.031671143951709085
    };

    for ( i=0; i< filterOrder; i++)
    {
      if (useotherFlag==1)
        impulseRep[i]=lp[i];
      else
        impulseRep[i]=hp[i];
    }	
    return;
  }
  // lowpass filter, order =N
  // M_PI in /usr/inculde/math.h
  LP1=fc/sf;
  M=(filterOrder-1)/2;
  for ( i=0; i< filterOrder; i++)
   tmp[i]=2.*M_PI*(i-M);

  for ( i=0; i< filterOrder; i++)
  {
    if( i == M)
      impulseRep[i]=2.*LP1;
   else
      impulseRep[i]=2.*sin(tmp[i]*LP1)/tmp[i];
  }
}


/**
 * \fn void sc2dafindcableadjVal(ARRAYSET *setup, int biasInx,
 *     StatusType *status)
 *
 * \brief function
 *  first find mean Value for each channel, then calculate  
 *  the cableadjust Val  for next bias setting 
 * 
 * \param  setup    ARRAYSET structure pointer  
 * \param  biasInx  int 
 * \param  status   StatusType pointer  
 *
 */
/*+ sc2dafindcableadjVal
*/
void sc2dafindcableadjVal
(
ARRAYSET   *setup, 
int        biasInx,
StatusType *status
)
{
  int        ch, i, overTHD,minInx,maxInx;
  int        midVal,maxVal,minVal,cabadjTHD;
  int        *cableadjPtr, *meanvalPtr;;
  double     adjVal;

  if (!StatusOkP(status)) return;

  // from char * to int *
  cableadjPtr= (int*)setup->cableadjPtr;
  meanvalPtr = (int*)setup->meanvalPtr;

  // the CHKSTART is used to avoid the first few data as when I look at 
  // data plot, it is a bit weird at beginning 
 
  for (ch=0;ch<COL_NUM;ch++)
  {
    // find the max and min, then the middle point
    maxVal=setup->allData[ch + COL_NUM*CHKSTART];
    minVal=setup->allData[ch + COL_NUM*CHKSTART];;
    for (i=CHKSTART; i<(setup->fdbkNo-1); i++)
    {
      if( setup->allData[ch + COL_NUM*i] >= maxVal) 
      {
        maxVal=setup->allData[ch + COL_NUM*i];
        maxInx=i;
      }
      else if( setup->allData[ch + COL_NUM*i] <= minVal) 
      {
        minVal=setup->allData[ch + COL_NUM*i];
        minInx=i;
      }
    } 
    // find the middle value
    midVal=(maxVal+minVal)/2;
    *(meanvalPtr + ch + biasInx*COL_NUM)=midVal;
  
    cabadjTHD=setup->cableadjThd[ch];

    // alway adjust if midVal>=0 , otherwise
    // allow | midVal| > cableadjTHD to adjust  05/01/07. X. Gao 
    if ( midVal < 0 )
    {
       if ( abs(midVal) > cabadjTHD  )
         overTHD= midVal+cabadjTHD;
       else
         overTHD=0;
    }
    else  
      overTHD= midVal;  

    // avoid ch if colMask[ch]=0
    // keep the previous one and only update if overTHD >0
    if (overTHD !=0 &&  setup->colMask[ch] !=0)
    {
      // setup->cableSlope[ch] =~ -70, avoid zero
      if(setup->cableSlope[ch] !=0 ) 
      {
        // divFlag=1, setup->cableSlope[ch]=A/N 
        // but overTHd also  = xx/N , Hence, it has no effect on adjVal
        adjVal= -setup->cableadjScale[ch]*(double)overTHD/setup->cableSlope[ch];
        setup->cableadjVal[ch] += (int) adjVal;
      }
    }
    #ifdef DEBMEANVAL
      printf(" ch_%d midVal=%d cabadjTHD=%d cabadjScale=%f,cabslope=%f,cableadjVal=%d \n",
            ch, midVal,cabadjTHD,setup->cableadjScale[ch],setup->cableSlope[ch],
            setup->cableadjVal[ch]);
    #endif
    // store meanVal for later writing
    *(cableadjPtr + ch + biasInx*COL_NUM)=setup->cableadjVal[ch];
  }
}


/**
 * \fn void sc2dafindlockPts(dasInfoStruct_t *myInfo,
 *   ARRAYSET *setup, char *databuf, FILE *fplck,
 *   char *servoFile, StatusType *status)
 *
 * \brief function
 *  find the max modulation from array of bias-setting data
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data 
 * \param  fplck   FILE point for the lock result
 * \param  servoFile  prefix filename for servodata.hex 
 * \param  status StatusType.  given and return
 *
 */
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
)
{
  int        bias, ch;
  int        *servodataPtr, *tmpData;
  char       *peakinfoPtr;
  BIAS_PEAK  *totalpeakPtr,*peakperbiasPtr;
  SERVO_DATAHEAD *rampinfoPtr;

  if (!StatusOkP(status)) return;

  // use SERVO_INFO
  servoInfo.totaldataperBias=COL_NUM*setup->fdbkNo;
  servoInfo.totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
  servoInfo.headSize = sizeof(SERVO_DATAHEAD);  
  servoInfo.totalframeSize = servoInfo.totalptsSize*COL_NUM;
  servoInfo.peakSize = sizeof(BIAS_PEAK)*setup->biasNo;
  servoInfo.maxpeakSize =  sizeof(MAX_P2P);
  servoInfo.totalservoSize = servoInfo.headSize + servoInfo.totalframeSize + 
                             servoInfo.peakSize + servoInfo.maxpeakSize;

  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);
  peakinfoPtr  = (databuf+ servoInfo.headSize+ servoInfo.totalframeSize);
   
  // point to the entry of array of BIAS_PEAK struture
  totalpeakPtr =(BIAS_PEAK *)peakinfoPtr;

  // for each bias setting, a set of data are taken at different
  // feedbk points
  for (bias=0; bias<setup->biasNo; bias++)
  {
    // point to the first data for each bias
    tmpData=servodataPtr + bias*servoInfo.totaldataperBias;  

    // point to the corresponding  BIAS_PEAK  for each bias
    peakperbiasPtr = totalpeakPtr + bias;
    for( ch=0;ch<COL_NUM; ch++)
    {
      sc2dalookeachCh(setup, tmpData, peakperbiasPtr,ch, status);  
      if(!StatusOkP(status)) 
      {
        ErsRep(0,status,"sc2dafindlockPts: failed to call sc2dalookeachCh");
        return ;
      }
    }   
  }
  
  // we have all peak-peak values for each bias setting, 
  // find max. 
  sc2dafindmaxModul(myInfo,setup,databuf,servoFile,status);

  if(setup->servo ==SQ1SERVO || setup->servo ==SQ1LOCK )
    sc2dafindrefmaxModul(myInfo,setup,databuf,status);

  // save binary "fileNameservodata.hex" data so that 
  // standalone findlockpoints can use it.
  if (setup->doServo !=SQ2SERVOSSALOCK)
      sc2dasaveservoData(setup,databuf,servoFile,status);
  
  // max mudoluation may have a flat bottom, we need to avoid it
  if( setup->servo==SSARAMP)
    sc2dafindflatBottom(&servoInfo,setup,databuf,status);

  sc2dafindinitlckPts(myInfo,setup,databuf,status);
  if (setup->doServo !=SQ2SERVOSSALOCK)
  {
    // sq1biasMax=setup->biaslckPt[ch]
    // sq2fbMax=setup->maxValue.zfact[ch];
    // sq1fbMax=setup->maxValue.initFB[ch]
    // sq1biasRef=setup->biasrefPt[ch]
    // sq2fbRef=setup->revValue..zfact[ch]; 
    // sq1fbRef=setup->refValue.initFB[ch];
    sc2dasavelockPts(myInfo,setup,databuf,fplck,status);
    sc2dasave4nextStep(myInfo,setup,databuf,fplck,status);
  }
}


/**
 * \fn void sc2dafindinitlckPts(dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, char *databuf, StatusType *status)
 *
 * \brief function
 *  after found maxModulation, from maxModuation data set
 *  find initial lock points in the middle range in X-axis
 * 
 * \param  myInfo   dasInfoStruct_t pointer
 * \param  setup    ARRAYSET structure pointer  
 * \param  databuf  char pointer for servo data
 * \param  status   StatusType pointer  G&R
 *
 */
/*+ sc2dafindinitlckPts
*/
void sc2dafindinitlckPts
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup, 
char       *databuf,
StatusType *status
)
{
  int            ch;
  int            *servodataPtr, *chdataperbiasPtr;
  SERVO_DATAHEAD *rampinfoPtr;

  if (*status != STATUS__OK) return;
 
  // get different pointers
  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);
 
  // check for each channel
  for (ch=0; ch<COL_NUM; ch++)
  {
    // point to the max-modulation bias point. but, first check
    if(setup->maxValue.peakVal[ch][0]!=NONMODULATION)
    {
      chdataperbiasPtr=
         servodataPtr +
         setup->maxValue.biasStep[ch][0]*servoInfo.totaldataperBias + 
         ch;

      // get ssaBias (SSARAMP) sq2Bias(SQ2SERVO) or sq1Bias (SQ1SERVO)
      if (setup->servo==SSARAMP || setup->servo ==SQ2SERVO || 
          setup->servo==SQ1SERVO )
      {         
        setup->biaslckPt[ch]=
             setup->maxValue.biasStep[ch][0]*setup->stepBIAS+setup->minBIAS;
      }
     
      //    (setup->servo ==SQ2SERVO )
      // since initFB, gain and Zfact are updated in ssalock, so
      // need to change them in sq2lock, no change for ssalock
      // update in 4sq2Step()
      // or see findmaxModul to see 
      //      setup->maxValue.initFB[ch]=peakperbiasPtr->initFB[ch];
      //      setup->maxValue.zfact[ch]=peakperbiasPtr->zfact[ch];
      //      setup->maxValue.gain[ch]=peakperbiasPtr->gain[ch];
     
      // find the optimal point near middle of X-axis 
      // ( either ssafb,sq2fb,or sq1fb ),
      // setup->saoutlckVal for ssa is Zfact, for sq2servo is ssafb, 
      // for sq1servo is sq2fb
      // for SSA, need to find gain as well
      // for SQ2 and SQ1, finding gain is done in (sq2,sq1)open stage
      sc2daservofindslopGain(myInfo,setup,ch, chdataperbiasPtr,status);

      // for sq2servo, find flux-period for sq1 flux-jump
      // can do that in SQ"OPEN since they are the same
      if (setup->servo ==SQ2SERVO || setup->servo ==SQ2LOCK) 
        sc2daservofindfluxPeriod(myInfo,setup,ch, chdataperbiasPtr,status);
    }
    else
    {
      if (setup->servo !=SSALOCK && setup->servo !=SQ2OPEN)        
       setup->biaslckPt[ch]=NONMOD_BIAS;
    }
  }
}


/**
 * \fn void sc2dafindlckptsGain(dasInfoStruct_t *myInfo,
 *   ARRAYSET *setup, char *databuf, FILE *fplck,
     char *servoFile, int *pixelMask, StatusType *status)
 *
 * \brief function
 *  find gain at lock points only for:
 *  SQ2OPEN SQ2OPEN4P at  sq2fb=initFB[i] sq2fdbkOpt[i]from ssaOut(sq2fb) curve 
 *  SQ1OPEN at  Zfactor from ssaOut(sq1fb) 
 *  for slopSelect[0] row, find gain from SQ2OPEN binary data for sq2fb servo
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data 
 * \param  fplck   FILE point for the lock result
 * \param  servoFile  prefix filename for servodata.hex 
 * \param  pixelMask int pointer 
 * \param  flag      int  0: for a=standalone prog, no use servoFile, 
 *                   don't call saveservoData
 * \param  status StatusType.  given and return
 *
 */
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
)
{
  int            ch, i,sqfbStep;
  int            l, tang;
  int            *servodataPtr, *chdataperbiasPtr=0;
  int            aver1,aver2;
  int            data3, data2, data1;

  int      startFB,halfWind,slopeWind=30;;
  double   xData[slopeWind],weightData[slopeWind],yData[slopeWind];
  double   slope,dclevel;
  SERVO_DATAHEAD *rampinfoPtr;

  if (*status != STATUS__OK) return;

  // use SERVO_INFO
  servoInfo.totaldataperBias=COL_NUM*setup->fdbkNo;
  servoInfo.totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
  servoInfo.headSize = sizeof(SERVO_DATAHEAD);  
  servoInfo.totalframeSize = servoInfo.totalptsSize*COL_NUM;
  servoInfo.peakSize = sizeof(BIAS_PEAK)*setup->biasNo;
  servoInfo.maxpeakSize =  sizeof(MAX_P2P);
  servoInfo.totalservoSize = servoInfo.headSize + servoInfo.totalframeSize + 
                             servoInfo.peakSize + servoInfo.maxpeakSize;
 
  // assign different pointers
  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);

  halfWind=slopeWind/2;
  if (setup->servo==SQ2OPEN  )
  {
    for (ch=0; ch<COL_NUM; ch++)
    {
      setup->gain[ch]=0;
      if(setup->initFB[ch]!=NONMODULATION)
      {
        sqfbStep=(setup->initFB[ch]-setup->minFDBK)/setup->stepFDBK;
        chdataperbiasPtr=servodataPtr + sqfbStep*COL_NUM + ch;
      
        // this is ssaOut: at ssaOut(sq2fb) for SQ2OPEN;
        setup->saoutlckVal[ch]=*(chdataperbiasPtr);

        // use least-square-fit 
        startFB=sqfbStep-halfWind;
        if (startFB <=0)
          startFB=0;
        for (i=0;i<slopeWind;i++)
        {
          xData[i]=i+1;
          weightData[i]=1;
          yData[i]=*(servodataPtr +(i+startFB)*COL_NUM +ch);
        }
        sc2damath_linfit(slopeWind, xData, yData,weightData, &slope,&dclevel,status); 
   
        jitDebug(16,"ch_%d  sqfb=%d ssaLock=%d, slope=%f  \n",
                ch,setup->initFB[ch],setup->saoutlckVal[ch],slope);

        if( slope!=0 )
          setup->gain[ch]=(double)(setup->stepFDBK)/slope;	   
      }
    }
  }
  else  if (setup->servo==SQ1OPEN)  // it is only a dummy 
  {
    for (ch=0; ch<COL_NUM; ch++)
    {
      setup->gain[ch]=0;
      for (i=0;i<(setup->fdbkNo-1);i++)
      { 
        chdataperbiasPtr=servodataPtr +i*COL_NUM +ch;
        //at ssaOut(sq1fb) for SQ1OPEN; 
        data1=*(chdataperbiasPtr);
        data2=*(chdataperbiasPtr+1);
        if ( (data1 < setup->zfact[ch] && data2 > setup->zfact[ch]) ||
             (data1 > setup->zfact[ch] && data2 < setup->zfact[ch]) )
        {
          setup->saoutlckVal[ch]=*(chdataperbiasPtr);
          break;
        }
      }
      aver1=aver2=0;
      for (l=0; l<AVER_NO; l++)
      {
        data2=*(chdataperbiasPtr+ (l         )*COL_NUM );
        data3=*(chdataperbiasPtr+ (l-AVER_NO )*COL_NUM );

        aver1 +=data2;
        aver2 +=data3;
      }
      // gain =1/tang
      tang=(aver1-aver2);
      if( tang!=0 )
        setup->gain[ch]=(double)(AVER_NO*AVER_NO*setup->stepFDBK)/tang;
    }
  }
  else
  {
    // nothing for SQ2OPEN4P, currently, later may move from findsq2fbGain to here
  }
  if ( flag==1)
    sc2dasaveservoData(setup,databuf,servoFile,status);

  sc2dasavelockPts(myInfo,setup,databuf,fplck,status);
  sc2dasave4nextStep(myInfo,setup,databuf,fplck,status); 
 
  // now for sq2fb servo gain, and file sq2open4p-FB-GZ
  // setup->slopSelect[0] > 0 =rowNo turned off for sq2fb servo 
  if ( flag==1)
  {
    jitDebug(16,"sc2dafindsq2fbGain: slopSelect[0]=%d\n", setup->slopSelect[0]);
    if ( setup->servo==SQ2OPEN4P  )  
      sc2dafindsq2fbGain(myInfo,setup,pixelMask,status); 
  }
}


/**
 * \fn void sc2dafindflatBottom(SERVO_INFO *myservoInfo,
 *     ARRAYSET *setup, char *databuf,StatusType *status)
 *
 * \brief function
 *  look for a flat bottom. from testing result, it is found that
 *  if | AVER(-3)-AVER(+3)|/stepFDBK < thredVal
 *    flat
 *  else 
 *    not 
 *
 *  this is not ture for more than two cycles ( SSA 3-1 coils). 
 *  However, since final SSA is (1-1 coil ), we can use this
 *
 * \param  myservoInfo SERVO_INFO pointer
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data 
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dafindflatBottom
*/
void sc2dafindflatBottom
(
SERVO_INFO    *myservoInfo,
ARRAYSET       *setup,
char           *databuf,
StatusType     *status
)
{
  int        i,j,ch,fbPts,deriFlag;
  int        stillFlat,a,b,*servodataPtr;
  char       *peakinfoPtr;
  int        *chdataperbiasPtr;
  BIAS_PEAK  *totalpeakPtr,*peakperbiasPtr;
  double      deriVal,thdVal;
  SERVO_DATAHEAD *rampinfoPtr;

  if (!StatusOkP(status)) return;
 
  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);
  peakinfoPtr  = (databuf+ myservoInfo->headSize + myservoInfo->totalframeSize);
   
  // point to the entry of array of BIAS_PEAK struture
  totalpeakPtr =(BIAS_PEAK *)peakinfoPtr;

  // point to max_modulation for each channel
  for (ch=0; ch<COL_NUM; ch++)
  {
    deriFlag=1; // 1: flat; 0 not flat
    stillFlat=1;
     // first check
    if(setup->nomodulCount[ch] !=setup->biasNo )
    {
      while(stillFlat==1)
      {
        // each bias has one peakEntry
        // point to the corresponding BIAS_PEAK entry for this bias
        // find one minpeak point 
        peakperbiasPtr=totalpeakPtr+setup->maxValue.biasStep[ch][0];
        for(i=0;i<3;i++)
        {
          if ( peakperbiasPtr->peakPts[ch][i] < peakperbiasPtr->peakPts[ch][i+1] )
            break;
        }
        fbPts=setup->maxValue.fdbk[ch][i];

        chdataperbiasPtr= servodataPtr +
              setup->maxValue.biasStep[ch][0]*myservoInfo->totaldataperBias + 
              ch;
   
        sc2dafindflatDeriv(setup,chdataperbiasPtr,fbPts,&deriFlag,&deriVal,&thdVal,
                           ch,status);
        if(deriFlag==0)
          stillFlat=0;
        // don't do if it is already the last one
        else if( (setup->maxValue.biasStep[ch][0]+1) >= setup->biasNo)
          stillFlat=0;
        else
        {
          // it is flat, it will be less flat in next bias
          setup->maxValue.biasStep[ch][0]++;

          // update the maxValue too
          peakperbiasPtr=totalpeakPtr+setup->maxValue.biasStep[ch][0];
          for(j=0; j<3; j++)
          {
            setup->maxValue.peakVal[ch][j]=peakperbiasPtr->p2pValue[ch][j];
            setup->maxValue.fdbk[ch][j]=peakperbiasPtr->peakInx[ch][j];
   
            a=peakperbiasPtr->peakPts[ch][j];
            b=peakperbiasPtr->peakPts[ch][j+1];          
            setup->maxValue.saOut[ch][j]=(a+b)/2;
          }
          setup->maxValue.fdbk[ch][3]=peakperbiasPtr->peakInx[ch][3];
        }
      }
    }
  }            
}


/**
 * \fn void sc2dafindflatDeriv(ARRAYSET *setup,  int *chdata, int fbStart,
 *     int *derivt, double *deriVal, double *thdVal, int ch, StatusType *status)
 *
 * \brief function
 *  find the derivative within shifted window period
 *         given and return ==>G&R        
 * 
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  chdata   int pointer for channel data 
 * \param  fbStart  int  Fdbk step No where minPeak is  
 * \param  derivt  int pointer for derivative flag, G&R
 * \param  deriVal   double pointer, pass deriVal
 * \param  thdVal    double pointer, pass thdVal
 * \param  ch        int, which ch
 * \param  status   StatusType pointer.  G&R
 *
 *  from testing result, it is found that
 *  if | AVER(-3)-AVER(+3)|/stepFDBK < thredVal
 *    flat 
 *  else 
 *   not  ( return *derivt=0)
 *
 */
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
)
{
  //double  aver1, aver2;
  //int     j,l,data1, data2;
  int     *chdataPtr;
  int     i;
  int     slopeWind=8;
  double   tmpDeriv[3];
  double  thredVal=FLAT_THDSSA;
  double  xData[slopeWind],weightData[slopeWind],yData[slopeWind];
  double  slope,dclevel;
  

  if (!StatusOkP(status)) return;

  if (setup->servo==SSARAMP)
  {
    thredVal=FLAT_THDSSA;

  }
  else if (setup->servo==SQ2SERVO || setup->servo==SQ2LOCK)
    thredVal=FLAT_THDSQ2;
  else if (setup->servo==SQ1SERVO || setup->servo==SQ1LOCK)
    thredVal=FLAT_THDSQ1;

  *thdVal=thredVal;
  
  // point to the minPeak point in fdbk ramping data for this ch
  chdataPtr=chdata + fbStart*COL_NUM;
  #ifdef PRINT_CH_FLAT
  {
    if (ch==TEST_CH)
     printf("ch-%d: flatthdVal=%f, fbStart=%d, firtsData=%d\n",
            ch, thredVal,fbStart,*chdataPtr);
  }
  #endif   

  // try to use this, get slope from both side of the minPeak, 
  // [-8,-1]
  for (i=0; i<slopeWind; i++)
  {
    xData[i]=i+1;
    weightData[i]=1;
    yData[i]=*(chdataPtr + (i-slopeWind)*COL_NUM);
  }
  sc2damath_linfit(slopeWind, xData, yData,weightData, &slope,&dclevel,status); 
  tmpDeriv[0] =slope/setup->stepFDBK;

  for (i=0; i<slopeWind; i++)
  {
    xData[i]=i+1;
    weightData[i]=1;
    yData[i]=*(chdataPtr+ (i+1)*COL_NUM);
  }
  sc2damath_linfit(slopeWind, xData, yData,weightData, &slope,&dclevel,status); 
  tmpDeriv[1] =slope/setup->stepFDBK;

  if ( fabs(tmpDeriv[0]) < thredVal ||  fabs(tmpDeriv[1]) < thredVal )
      *derivt=1;
   else
      *derivt=0;

  #ifdef PRINT_CH_FLAT
  {
    if (ch==TEST_CH)
    {
      printf("tmpDeriv[0]=%f tmpDeriv[1]=%f deriFlag=%d\n",
         tmpDeriv[0],tmpDeriv[1],*derivt);
    }
  }
  #endif   

}


/**
 * \fn void sc2dafindframelckptsGain(dasInfoStruct_t *myInfo,  ARRAYSET *setup, 
 *     int * pixelMask, StatusType *status)
 *
 * \brief function
 *  find gain at lock points only for:
 *  SQ1OPEN at ssaout=Zfactor  from ssaOut(sq1fb) curve
 *  from frame data 
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  pixelMask int pointer 
 * \param  status   StatusType pointer  G&R
 *
 */
/*+ sc2dafindframelckptsGain
*/
void sc2dafindframelckptsGain
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int             *pixelMask,
StatusType      *status
)
{
  char   *frameData,tmpFile[FILE_LEN];
  int    *framedataPtr,*framePtr;
  int    fileLen;
  int    headNo=FRAMEHEADER_NUM;
  int    fdbk,frameSize,row,col,pixel;
  BIAS_PEAK  peakInfo;  

  if (!StatusOkP(status)) return;

  if (setup->sampleNo ==0)
  {
    *status=DITS__APP_ERROR;
    sprintf(msg,"sc2dafindframelckptGain: sampleNo=0 ");
    ErsRep(0,status,msg);
    return;
  }

  // open the pixel file, but no appending
  sprintf(tmpFile,"%s-binary",myInfo->dataFile); 
  if((myInfo->fpSq1 = fopen64(tmpFile,"r")) == NULL)
    {
      *status = DITS__APP_ERROR;
      sprintf(msg,"sc2dafindframelckptGain: failed to open %s",tmpFile);
      ErsRep(0,status,msg);
      return;
    }
  sc2dareadbinaryfile2Mem(myInfo,myInfo->fpSq1, &frameData, &fileLen,status);
  if (!StatusOkP(status)) 
  {
    sprintf(msg,"sc2dafindframelckptGain: failed to call sc2dareadbinaryfile2Mem ");
    ErsRep(0,status,msg);
    return;
  }
  fclose(myInfo->fpSq1);

  if ( setup->slopSelect[15] >=0 )
  {
    fprintf(myInfo->fpLog,"use filtered data (flag=%d)\n",setup->slopSelect[15]);
  }

  // myInfo->bufsize is in byte, including CHKSUM
  frameSize=myInfo->bufSize/4;
  framePtr= (int *)frameData; 
  for (row=0;row<ROW_NUM;row++)
  { 
    for (col=0;col<COL_NUM;col++)
    {
      pixel=row*COL_NUM+col;
      setup->pixeLock.ssalockPt[pixel]=0;
      setup->pixeLock.sq1FB[pixel]=0;
      setup->pixeLock.gain[pixel]=0;
      setup->pixeLock.iVal[pixel]=0;

      // first check, we use colMask and rowMask, we need sq1fblck from sq1servo
      // it is sq1fblck[ROW_NUM*COL_NUM]
      if(setup->colMask[col] !=0 && setup->rowMask[row]!=0)
      {
       #ifndef FINDLOCKPOINT_PROG
        if ( (setup->pixeLock.sq1initFB[pixel] !=NONMODULATION ) &&
             (setup->sq2fdbkOpt[col] !=NONMODULATION) &&
             (setup->biasoptPt[row] !=NONMODULATION)   )
       #endif
        {
          // collect the pixel data from all frame for each fdbk
          for ( fdbk=0;fdbk<setup->fdbkNo;fdbk++)
          {
            // move to the right entry for each fdbk
            // setup->allData for hold all pixel's data 
            framedataPtr=framePtr + headNo + frameSize*fdbk;
            setup->allData[fdbk]=*(framedataPtr+pixel);
          }
          sc2dafindpixellckptsGain2(myInfo,setup,row,col,&peakInfo,pixelMask,status);
          if (!StatusOkP(status)) 
          {
            sprintf(msg,
              "sc2dafindframelckptGain2: failed to call sc2dafindpixellckptsGain ");
            ErsRep(0,status,msg);
            return;
          }
        }
      }
    }
  }
}



/**
 * \fn void sc2dafindpixellckptsGain2(dasInfoStruct_t *myInfo,  ARRAYSET *setup, 
 *     int row, int col, BIAS_PEAK *peakInfo, int *pixelMask, StatusType *status)
 *
 * \brief function
 *  find gain only for:SQ1OPEN around the middle point after first peak (either max or min) 
 *  from ssaOut(sq1fb) curve from pixel data, 
 *  peakInfo is not used, only for allowing lookabsMax to work
 *
 * \param myInfo    dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  row      int
 * \param  col      int
 * \param  peakInfo BIAS_PEAK pointer for the peak, not used
 * \param  pixelMask int pointer
 * \param  status   StatusType pointer  G&R
 *
 */
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
)
{
  int     pixel, chkmaxFlag=0;
  int     sq1fbStep, j,startFB,slopeWind=30, halfWind=0; 
  int     data=0, data1=0, midVal=0,start, end;
  double  xData[slopeWind],weightData[slopeWind],yData[slopeWind];
  double  slope,dclevel,ival;
  double  ivalScale=1;

  if (*status != STATUS__OK) return;
  
  pixel=row*COL_NUM + col;
  setup->pixeLock.gain[pixel]=0;
  setup->pixeLock.sq1FB[pixel]=0;
  setup->pixeLock.ssalockPt[pixel]=0;
  setup->pixeLock.fluxPeriod[pixel]=0;
  setup->pixeLock.peakInx[pixel][0]=0;
  setup->pixeLock.peakInx[pixel][1]=0;
  ivalScale=(double)setup->slopSelect[6]/(double)setup->slopSelect[7];

  if (!StatusOkP(status)) 
  {
    sprintf(msg,"sc2dafindpixellckptsGain:Error- status!=OK ");
    ErsRep(0,status,msg);
    return;
  }
  // check if ABS(max-min)>=PEAK_THD.. allData[] is int*
  // get peakInfo->peakInx[col][0]=minInx, -1]=maxInx
  sc2dalookabsMax(setup->allData,setup,col,&chkmaxFlag,0,peakInfo,status);
  // it is >=PEAK_THD, we start to look the lock point
  if ( chkmaxFlag !=0)
  {
    // find the middle point 
    midVal=( peakInfo->highest[col] + peakInfo->lowest[col] )/2;
    
    // search the ssalock from 
    if ( peakInfo->peakInx[col][0] < peakInfo->peakInx[col][1] )
    {
      start=peakInfo->peakInx[col][0];
      end  =peakInfo->peakInx[col][1];
    }
    else
    {
      start=peakInfo->peakInx[col][1];
      end  =peakInfo->peakInx[col][0];
    }
    // now serach for middle
    for (j=start; j<end; j++)
    {
      data=setup->allData[j];
      data1=setup->allData[j+1];
      if ( (data < midVal &&  data1> midVal) || 
           (data > midVal &&  data1< midVal) )
      {
        break;
      }
    }
    sq1fbStep=j;
   
    setup->pixeLock.sq1FB[pixel]=sq1fbStep*setup->stepFDBK + setup->minFDBK;

    if (setup->servo ==SQ2OPEN  || setup->servo ==SQ2OPEN4P)  
      setup->pixeLock.ssalockPt[pixel]=data;
    else  // this is normalised by sampleNo
      setup->pixeLock.ssalockPt[pixel]=data/setup->sampleNo;

    // for finding pixelfluxperiod
    setup->pixeLock.midInx[pixel]=sq1fbStep;
    setup->pixeLock.midVal[pixel]=data;

    if (setup->servo ==SQ1OPEN || (setup->servo ==SQ2OPEN4P && setup->slopSelect[5]==0) )
    {  
      // now calculate gain, use least square fit. the initial start point:
      halfWind=slopeWind/2;    startFB=sq1fbStep-halfWind;
      if (startFB<=0)          startFB=0;

      for (j=0; j<slopeWind; j++)
      {
        xData[j]=j+1;
        weightData[j]=1;
        yData[j]=setup->allData[startFB+j];
      }
      sc2damath_linfit(slopeWind, xData, yData,weightData, &slope,&dclevel,status); 

      if( slope!=0 )
      {
        if (setup->servo ==SQ2OPEN || setup->servo==SQ2OPEN4P  )  
          setup->pixeLock.gain[pixel]=(double)(setup->stepFDBK)/slope;
        else   // normalised by sampleNo, pixeLock.gain[pixel] initialVal=0;
        {
          setup->pixeLock.gain[pixel]=setup->sampleNo*(double)(setup->stepFDBK)/slope;

           // left shift 12 = *4096
          ival=( fabs( (double)setup->stepFDBK/slope) * 4096) ;
	    if (slope <0)
            setup->pixeLock.iVal[pixel]= -(int)(ival*ivalScale*pixelMask[pixel]);
          else
            setup->pixeLock.iVal[pixel]= (int)(ival*ivalScale*pixelMask[pixel]);

        }
      }
    }
    // find flux period for MCE to do flux jump if required 
    sc2dafindpixelfluxPeriod(myInfo,setup,pixel,status);
  }               
}


/**
 * \fn void sc2dafindsq2fbGain(dasInfoStruct_t *myInfo,  ARRAYSET *setup, 
 *     int *pixelMask, StatusType *status)
 *
 * \brief function
 *  find gain at middle point for SQ1-Bias turned off row from
 *  SQ2OPEN  if slopSelect[0]> 0  
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  pixelMask int pointer 
 * \param  status   StatusType pointer  G&R
 *
 * setup->slopSelect[5]=0: find lock point from middlepoint.  
 *                     =1: find lock point from sq2fbOPT 
 */
/*+ sc2dafindsq2fbGain
*/
void sc2dafindsq2fbGain
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int             *pixelMask,
StatusType      *status
)
{
  char   *frameData;
  char   tmpFile[FILE_LEN];
  int    *framedataPtr,*framePtr;
  int    fileLen, foundFB=0;
  int    headNo=FRAMEHEADER_NUM;
  int    fdbk,frameSize,row,col,pixel;
  FILE   *fpsq2fb;
  BIAS_PEAK  peakInfo;  
  int    i,slopeWind=20,halfWind, startFB;
  double xData[slopeWind],weightData[slopeWind],yData[slopeWind];
  double slope, dclevel;

  if (!StatusOkP(status)) return;

  // open the binary file,
  sprintf(tmpFile,"%s-binary",myInfo->dataFile); 
  if((myInfo->fpSq1 = fopen64(tmpFile,"r")) == NULL)
  {
    *status = DITS__APP_ERROR;
    sprintf(msg,"sc2dafindsq2fbGain: failed to open %s",tmpFile);
     ErsRep(0,status,msg);
     return;
  }
  sc2dareadbinaryfile2Mem(myInfo,myInfo->fpSq1, &frameData, &fileLen,status);
  if (!StatusOkP(status)) 
  {
    sprintf(msg,"sc2dafindsq2fbGain: failed to call sc2dareadbinaryfile2Mem ");
    ErsRep(0,status,msg);
    return;
  }
  fclose(myInfo->fpSq1);
  // myInfo->bufsize is in byte, including CHKSUM
  frameSize=myInfo->bufSize/4;
  row=setup->slopSelect[0];
  framePtr= (int *)frameData;
  
   
  // save for other package or Execl to plot
  sprintf (tmpFile, "%s-open4p", myInfo->dataFile ); 
  if((fpsq2fb = fopen(tmpFile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dafindsq2fbGain 2: Error- failed to open file %s", tmpFile); 
      return;
    }
  fprintf(fpsq2fb,"# this is column data from SQ2OPEN4P Row(%d) for sq2fb servo\n",
          setup->slopSelect[0]);
  fprintf(fpsq2fb,"# first col is fdbk \n");
  fprintf(fpsq2fb,"# second col is col_0 \n");
  fprintf(fpsq2fb,"# ................\n");
  fprintf(fpsq2fb,"# last col is col_31 \n\n");

  for ( fdbk=0; fdbk<setup->fdbkNo; fdbk++)
  {
    fprintf(fpsq2fb,"%8d ", setup->minFDBK + fdbk*setup->stepFDBK);
    for (col=0; col<COL_NUM; col++)
    {
      pixel=row*COL_NUM + col;
      framedataPtr=framePtr + headNo + frameSize*fdbk;
      fprintf(fpsq2fb,"%8d ",*(framedataPtr + pixel) );
    } 
    fprintf(fpsq2fb,"\n");
  }
  fclose(fpsq2fb);

  for (col=0; col<COL_NUM; col++)
  {
    pixel=row*COL_NUM + col;
    setup->saoutlckVal[col]=0;
    setup->initFB[col]=0;
    setup->gain[col]=0;
    setup->fluxPeriod[col]=0;

    // collect the pixel data from all frame for each fdbk
    for ( fdbk=0; fdbk<setup->fdbkNo; fdbk++)
    {
      framedataPtr=framePtr + headNo + frameSize*fdbk;
      setup->allData[fdbk]=*(framedataPtr + pixel);
    }
 
    if(setup->colMask[col] !=0 )
    {
      if(setup->slopSelect[5] ==1)
      {
        if(setup->sq2fdbkOpt[col]!=NONMODULATION)
        {
          // find the lock point at sq2fbOpt
          fdbk=(setup->sq2fdbkOpt[col]-setup->minFDBK)/setup->stepFDBK;
          foundFB=setup->minFDBK + fdbk*setup->stepFDBK;

         // use least-square-fit 
	    halfWind=slopeWind/2;
          startFB=fdbk-halfWind;
	
          if (startFB <=0)   startFB=0;
	
          for (i=0;i<slopeWind;i++)
          {
            xData[i]=i+1;
            weightData[i]=1;
            yData[i]=setup->allData[i+startFB];
          }
          sc2damath_linfit(slopeWind, xData, yData,weightData, &slope,&dclevel,status); 
          jitDebug(16,"ch_%d  sq2fdbkOpt=%d ssaLock=%d, slope=%f  \n",
                  col,setup->sq2fdbkOpt[col],setup->saoutlckVal[col],slope);
          if( slope!=0 )
          {
            setup->initFB[col]=setup->sq2fdbkOpt[col];
            setup->gain[col]=(double)(setup->stepFDBK)/slope;	   

            // sq2fbOPT[col] may not = N*stepFDBK
            setup->saoutlckVal[col]=setup->allData[fdbk]+ ((setup->initFB[col]-foundFB)*slope)/setup->stepFDBK;
	    
            sprintf(msg,"col_%2d:d[%3d]=%8d [%3d]=%8d  lckFB=%5d  foundFB=%5d slope=%6.3f lock=%8d step=%3d",
                   col,fdbk,setup->allData[fdbk], fdbk+1, setup->allData[fdbk+1],
		   setup->initFB[col],foundFB,(slope/(double)setup->stepFDBK),setup->saoutlckVal[col], setup->stepFDBK); 
            MsgOut(status, msg);
          }
	}  
      }
      // use it for flux and for lockpoints if slopSelet[5]=0
      sc2dafindpixellckptsGain2 (myInfo,setup,row,col,&peakInfo,pixelMask,status);
      setup->fluxPeriod[col]=setup->pixeLock.fluxPeriod[pixel];

      // find lockpoint from middle point
      if(setup->slopSelect[5] ==0)
      {
        setup->saoutlckVal[col]=setup->pixeLock.ssalockPt[pixel];
        setup->initFB[col]=setup->pixeLock.sq1FB[pixel];
        setup->gain[col]=setup->pixeLock.gain[pixel];  
      }
    }
  }
  jitDebug(16,"sc2dafindsq2fbGain: slopSelect[0]=%d\n", setup->slopSelect[0]);
  sc2dasavesq2open4pfbZG (myInfo, setup,  status);
}




/**
 * \fn void sc2dafindpixelfluxPeriod(dasInfoStruct_t *myInfo,  ARRAYSET *setup, 
 *     int pixel, StatusType *status)
 *
 * \brief function
 *  find fluxperiod from ssaout(sq1fb) curve, for each pixel from sq1open
 *  from SG test, it seems that there are more than one period  
 *         given and return ==>G&R        
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  pixel    int
 * \param  status         StatusType pointer  G&R
 *
 */
/*+ sc2dafindpixelfluxPeriod
*/
void sc2dafindpixelfluxPeriod
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup,
int             pixel,
StatusType      *status
)
{
  int    j,i,l,fstEnd,flag=0,startPt=5;
  int    fbStart,fbMiddle=0,fbEnd;
  int     data0,data2=0, data3=0,data4=0, data5=0;
  int     ssafitWind=20;     
  double  yData[ssafitWind],xData[ssafitWind],weightData[ssafitWind];
  double  slope, dclevel;

  if (*status != STATUS__OK) return;
  
  fbStart=startPt;
  fstEnd=setup->fdbkNo-2;

  // first find the data<> midVal
  // we start from fifth point from ssaout(sq2fb)
  data0=setup->pixeLock.midVal[pixel];

  for (j=fbStart;j<fstEnd;j++)
  {
    l=j+5;
    data2=setup->allData[j];
    if ( (l+1)== fstEnd )
      l=fstEnd;
    data3=setup->allData[l];
    if( ( (data2 <= data0)  &&  ( data3 > data0) ) ||
        ( (data2 >= data0)  &&  ( data3 < data0) ) 
      )
      break;
  }

  fbStart=j;
  while (flag==0)
  {
    data0=setup->allData[fbStart];

    // find the slope here
    for (i=0; i<ssafitWind; i++)
    {
      xData[i]=i+1;
      weightData[i]=1;
      yData[i]=setup->allData[fbStart+i];
     }
     sc2damath_linfit(ssafitWind, xData, yData,weightData, &slope,&dclevel,status); 

    //start to look for ssaOut >=data0 or <= data0 depending on
    // the slop   first half period
    fbMiddle=fbStart;
  
    for (j=(fbStart+10);j<fstEnd;j++)
    {
      l=j+5;
      data2=setup->allData[j];
      if ( (l+1)== fstEnd )
        l=fstEnd;
	
      data3=setup->allData[l];
      if ( slope > 0 )
      {
        if( (data2 <= data0)  &&  ( data3 <= data0)) 
        { 
          fbMiddle=j;	  flag=1;
          break;
        }
      }
      else //if (slope < 0 )
      {
        if ( (data2 >= data0) && ( data3 >= data0) )
        {
          fbMiddle=j;	  flag=1;
          break;
        }
      }
    }
    if (flag==1)
      break;
    else
    {
      fbStart++;
      if (fbStart > (startPt+5) )
        flag=2;
    }  
  }
   
  fbEnd=setup->fdbkNo;
  if (flag==1)
  {
    // second half period
    for (j=fbMiddle+5;j<setup->fdbkNo-1;j++)
    {
      data4=setup->allData[j];
      data5=setup->allData[j+1];
      if (slope > 0 )
      {
        if( ( data4 >= data0) &&  (data5 >= data0 ) )
        { 
          fbEnd=j;
          break;
        }
      }
      else //if (slope < 0 )
      {
        if ( (data4 <= data0)  && (data5 <= data0) )
        {
          fbEnd=j;
          break;
        }
      }
    }
    setup->pixeLock.fluxPeriod[pixel]=(fbEnd-fbStart)*setup->stepFDBK;
  }
  else
    setup->pixeLock.fluxPeriod[pixel]=(fbEnd-fbStart)*setup->stepFDBK;

  if (pixel==TEST_PIXEL)
  {    
    #ifdef DEBFLUX1
    {
      printf("==== Fluxperiod[pixel=%d]=%d     fbStart[%d],fbEnd[%d]  flag= %d ==pixel\n",
           pixel, setup->pixeLock.fluxPeriod[pixel],fbStart,fbEnd,flag);
      printf("  slope(%f):data0[%d] middle[%d]:={[%d],[%d]} end[%d]:={[%d],[%d]}\n\n",
	       slope, data0,fbMiddle,data2, data3,fbEnd, data4, data5 );    
    }
    #endif
    if (myInfo->fpLog !=NULL)
    {	              
      fprintf(myInfo->fpLog,"==== Fluxperiod[ch=%d]=%d     fbStart[%d],fbEnd[%d]  flag= %d ==pixel\n",
           pixel, setup->pixeLock.fluxPeriod[pixel],fbStart,fbEnd,flag);
      /* fprintf(myInfo->fpLog,"  slope(%f):data0[%d] middle[%d]:={[%d],[%d]} end[%d]:={[%d],[%d]}\n\n",
	 slope, data0,fbMiddle,data2, data3,fbEnd, data4, data5 ); */          
    }
  }
}



/**
 * \fn void sc2dafindmaxMin(double *chdata, BIAS_PEAK *peakInfo,
 *  ARRAYSET *setup,int ch, StatusType *status)
 *
 * \brief function
 *  find two peaks: first max, then min
 *         given and return ==>G&R        
 * 
 * \param  chdata   double pointer for channel data 
 * \param  peakInfo  BIAS_PEAK pointer for peak points/ch, G&R
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  ch       which channel is looked at
 * \param  status   StatusType pointer.  G&R
 *
 */
/*+ sc2dafindmaxMin
*/
void sc2dafindmaxMin
(
double     *chdata, 
BIAS_PEAK  *peakInfo,
ARRAYSET   *setup,
int        ch,
StatusType *status
)
{
  int      localInx,order,maxInx;
  int      peakTHD=PEAK_THDSSA;

  if (!StatusOkP(status)) return;

  // move away from the first few data
  localInx=CHKSTART;

  // Choose a threshold based on the servo type
  if(setup->servo==SQ1SERVO || setup->servo==SQ1LOCK)
    peakTHD=PEAK_THDSQ1;
  else if ( setup->servo==SQ2SERVO || setup->servo==SQ2LOCK)
    peakTHD=PEAK_THDSQ2;
  else if  ( setup->servo==SSARAMP || setup->servo==SSALOCK )
  {
    peakTHD=PEAK_THDSSA;
    if (setup->divFlag==1)
      peakTHD = PEAK_THDSSA * setup->sampleNo;
  }


  while (1)
  {
    // Find the first maximum
    order=0;
    sc2dalookMax(chdata,peakInfo,setup,ch,&localInx,order,status);
    maxInx = localInx;

    // Find the next minimum
    order ++;
    sc2dalookMin(chdata,peakInfo,setup,ch,&localInx,order,status);
    peakInfo->p2pValue[ch][0]= peakInfo->peakPts[ch][0]
                          -peakInfo->peakPts[ch][1];


    // Test if the absolute value of the peak to peak modulation is greater than our threshold
    // Absolute is used because we allow the ptop value to be negative 
    if ( abs(peakInfo->p2pValue[ch][0])>peakTHD ) 
    {
      if (setup->servo==SSARAMP || setup->servo==SSALOCK )
      {
	// Find the next maximum
        order++;
        sc2dalookMax(chdata,peakInfo,setup,ch,&localInx,order,status);
        peakInfo->p2pValue[ch][1]= peakInfo->peakPts[ch][1]
                          -peakInfo->peakPts[ch][2];

	// Find the next minimum
        order ++;
        sc2dalookMin(chdata,peakInfo,setup,ch,&localInx,order,status);
        peakInfo->p2pValue[ch][2]= peakInfo->peakPts[ch][2]
                            -peakInfo->peakPts[ch][3];
      }
    }
    else
    {
      fprintf(myFpLog, "sc2dafindmaxMin ch: %d set to NONMODULATION p2p %f thd: %d  ",
	      ch,peakInfo->p2pValue[ch][0],peakTHD);
      fprintf(myFpLog, "max: %f at %d  min: %f at %d\n",
	      peakInfo->peakPts[ch][0], maxInx, peakInfo->peakPts[ch][1], localInx);

      peakInfo->p2pValue[ch][0]=(double)NONMODULATION;
      peakInfo->p2pValue[ch][1]=(double)NONMODULATION;
      peakInfo->p2pValue[ch][2]=(double)NONMODULATION;
    }
    if (  setup->servo==SSARAMP || setup->servo==SSALOCK)
    {
      // carry on to next four if we have more than two cycles 
      if ( (peakInfo->peakInx[ch][3]*setup->stepFDBK+setup->minFDBK) > 30000 ||
          peakInfo->p2pValue[ch][0] ==NONMODULATION )
        break;
    }
    else
       break;
  }
}


/**
 * \fn void sc2dafindmeanVal(dasInfoStruct_t *myInfo,
 *   ARRAYSET *setup, char *databuf, FILE *fplck,
     char *servoFile, StatusType *status)
 *
 * \brief function
 *  find the mean Value from array of data
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data 
 * \param  fplck   FILE point for the lock result
 * \param  servoFile  prefix filename for servodata.hex 
 * \param  status StatusType.  given and return
 *
 */
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
)
{

  int        bias, ch;
  int        *servodataPtr, *tmpData;
  char       *peakinfoPtr;
  BIAS_PEAK  *totalpeakPtr,*peakperbiasPtr;
  SERVO_DATAHEAD *rampinfoPtr;


  if (!StatusOkP(status)) return;

  // use SERVO_INFO
  servoInfo.totaldataperBias=COL_NUM*setup->fdbkNo;
  servoInfo.totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
  servoInfo.headSize = sizeof(SERVO_DATAHEAD);  
  servoInfo.totalframeSize =servoInfo.totalptsSize*COL_NUM;
  servoInfo.peakSize = sizeof(BIAS_PEAK)*setup->biasNo;
  servoInfo.maxpeakSize =  sizeof(MAX_P2P);
  servoInfo.totalservoSize = servoInfo.headSize + servoInfo.totalframeSize + 
                             servoInfo.peakSize + servoInfo.maxpeakSize;


  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);
  peakinfoPtr  = (databuf + servoInfo.headSize + servoInfo.totalframeSize);
   
  // point to the entry of array of BIAS_PEAK struture
  totalpeakPtr =(BIAS_PEAK *)peakinfoPtr;

  for (bias=0; bias<setup->biasNo; bias++)
  {
    // point to the first data for each bias ( or cableoffset)
    tmpData=servodataPtr + bias*servoInfo.totaldataperBias;  

    // point to the corresponding  BIAS_PEAK  for each bias
    peakperbiasPtr=totalpeakPtr+bias;

    for( ch=0;ch<COL_NUM; ch++)
    {
      sc2dalookmeaneachCh(setup, tmpData, peakperbiasPtr,ch, status);  
      if(!StatusOkP(status)) 
      {
        ErsRep(0,status,"sc2dafindmeanVal: failed to call sc2dalookmeaneachCh");
        return ;
      }
    }   
  }
}


/**
 * \fn void sc2dafindminMax(double *chdata, BIAS_PEAK *peakInfo,
 *  ARRAYSET *setup,int ch,StatusType *status)
 *
 * \brief function
 *  find two peaks: first min, then max
 *         given and return ==>G&R        
 * 
 * \param  chdata   double pointer for channel data 
 * \param  peakInfo  BIAS_PEAK pointer for the peak points/ch, G&R
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  ch       which channel is looked at
 * \param  status   StatusType pointer.  G&R
 *
 */
/*+ sc2dafindminMax
*/
void sc2dafindminMax
(
double     *chdata, 
BIAS_PEAK  *peakInfo,
ARRAYSET   *setup,
int        ch,
StatusType *status
)
{
  int  localInx, order, minInx;
  int  peakTHD=PEAK_THDSSA;

  if (!StatusOkP(status)) return;

  // move away from the first few data
  localInx=CHKSTART;

  // Select a threshold based on the type of servo we are running

  if(setup->servo==SQ1SERVO || setup->servo==SQ1LOCK)
    peakTHD=PEAK_THDSQ1;
  else if ( setup->servo==SQ2SERVO || setup->servo==SQ2LOCK)
    peakTHD=PEAK_THDSQ2;
  else if  ( setup->servo==SSARAMP || setup->servo==SSALOCK )
  {
    peakTHD=PEAK_THDSSA;
    if (setup->divFlag==1)
      peakTHD = PEAK_THDSSA * setup->sampleNo;
  }

  if (ch==TEST_CH)
  {
    #ifdef DEBSEREACH
      printf("findminmax: ch_%d thdVal=%d\n ",ch,peakTHD);
    #endif
  }

  while (1)
  {
    // Find the first minumum in the data
    order=0;
    sc2dalookMin(chdata,peakInfo,setup,ch,&localInx,order,status); 
    minInx = localInx;

    // Then find the next maximum
    order ++;
    sc2dalookMax(chdata,peakInfo,setup,ch,&localInx,order,status);
    peakInfo->p2pValue[ch][0]= peakInfo->peakPts[ch][0]
                            -peakInfo->peakPts[ch][1];


    // Check that the absolute peak to peak value is greater than our threshold
    // Note that abs is used here because the max was subtracted from the min above   
    if ( abs(peakInfo->p2pValue[ch][0])>peakTHD ) 
    {
      if ( setup->servo==SSARAMP || setup->servo==SSALOCK )
      {
	// Continue on finding mins and maxes and computing PtoP values 
        order ++;
        sc2dalookMin(chdata,peakInfo,setup,ch,&localInx,order,status);
        peakInfo->p2pValue[ch][1]= peakInfo->peakPts[ch][1]
                            -peakInfo->peakPts[ch][2];
        order++;
        sc2dalookMax(chdata,peakInfo,setup,ch,&localInx,order,status);
        peakInfo->p2pValue[ch][2]= peakInfo->peakPts[ch][2]
                             -peakInfo->peakPts[ch][3];
      }
    }
    else
    {
      fprintf(myFpLog, "sc2dafindminMax ch: %d set to NONMODULATION p2p %f thd: %d  ",
	      ch,peakInfo->p2pValue[ch][0],peakTHD);
      fprintf(myFpLog, "max: %f at %d  min: %f at %d\n",
	      peakInfo->peakPts[ch][1], localInx, peakInfo->peakPts[ch][0], minInx);

      peakInfo->p2pValue[ch][0]=(double)NONMODULATION;
      peakInfo->p2pValue[ch][1]=(double)NONMODULATION;
      peakInfo->p2pValue[ch][2]=(double)NONMODULATION;
    }
    if (  setup->servo==SSARAMP || setup->servo==SSALOCK)
    {
      // carry on to next four if we have more than two cycle 
     if ( (peakInfo->peakInx[ch][3]*setup->stepFDBK+setup->minFDBK) > 30000 ||
          peakInfo->p2pValue[ch][2] ==NONMODULATION )
       break;
    }
    else
      break;
  }
}


/**
 * \fn void sc2dafindmaxModul(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     char *databuf, char *fileName,StatusType *status)
 *
 * \brief function
 *  find maximum modulation from all bias setting data 
 *         given and return ==>G&R        
 *
 * \param  myInfo   dasInfoStruct_t pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  databuf  char pointer for servo data  
 * \param  fileName char pointer for first part of servodata.hex  
 * \param  status   StatusType pointer  G&R
 *
 */
/*+ sc2dafindmaxModul
*/
void sc2dafindmaxModul
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup, 
char       *databuf,
char       *fileName,
StatusType *status
)
{
  int            j,ch,bias, others;
  double         a,b;
  int            *servodataPtr;
  char           *peakinfoPtr, *maxpeakinfoPtr;  
  BIAS_PEAK      *peakInfo,*peakperbiasPtr;
  SERVO_DATAHEAD *rampinfoPtr;
  MAX_P2P        *maxInfo;

  if (*status != STATUS__OK) return;

  // initial all to NONMODULATION
  for (ch=0;ch<COL_NUM;ch++)
  {
    setup->saoutlckVal[ch]=0;
    setup->gain[ch]=(double)0;
    if (setup->servo ==SSARAMP || setup->servo ==SQ2SERVO || setup->servo ==SQ1SERVO )
    {
      setup->biaslckPt[ch]=NONMODULATION;
      setup->initFB[ch]=NONMODULATION;
    }
  }
  // get different pointers
  rampinfoPtr  =  (SERVO_DATAHEAD *)databuf;
  servodataPtr =  (int *)(rampinfoPtr +1);
  peakinfoPtr  =  (databuf + servoInfo.headSize + servoInfo.totalframeSize);
  maxpeakinfoPtr =(peakinfoPtr + servoInfo.peakSize);

  // point to array of biasNo  BIAS_PEAK struture
  peakInfo =(BIAS_PEAK *)peakinfoPtr;

  // point to memory area for holding each ch's peak data
  maxInfo=(MAX_P2P *)maxpeakinfoPtr;

  // check the three p2pValue
  for (ch=0; ch<COL_NUM; ch++)
  {
    // for all these, if biaslckPt[] ==NONMOD_BIAS, skip and set no modulation
    if ( ( setup->servo==SSALOCK ||setup->servo==SQ2OPEN || setup->servo==SQ1OPEN ||
           setup->servo==SQ2LOCK || setup->servo==SQ1LOCK ) &&  
           (setup->biaslckPt[ch] ==NONMOD_BIAS)  
       )	 
      others=1;
    else if (setup->nomodulCount[ch]==setup->biasNo)
      others=1;
    else 
      others=0;

    // found that all three |p2pValue| are similar, so, only check one    
    // always starts from the peakInfo
    // we need to readjust, p2pVal must have > 0, so use 0
    if  (peakInfo->p2pValue[ch][0]==NONMODULATION)
      maxInfo->maxp2p[ch][0]=0;
    else
      maxInfo->maxp2p[ch][0]=peakInfo->p2pValue[ch][0];
   
    for (bias=0; bias<setup->biasNo; bias++)
    {
      // point to the corresponding BIAS_PEAK entry for each bias
      // if biasNo=1, there is only one, which is max for SSALOCK, SQ2LOCK, SQ1LOCK
      peakperbiasPtr=peakInfo+bias;

      // skip this bias setting, set p2p=0, later set to NONMODULATION
      // SSALOCK if biaslckPt=NOMOD_BIAS, 
      if( others==1)
      {
        maxInfo->maxp2p[ch][0]=0;
        maxInfo->maxinx0[ch][0]=0;
      }
      else
      {
        // skip NONMODULATION data
        if ( peakperbiasPtr->p2pValue[ch][0] !=NONMODULATION)
        {
          if( fabs(peakperbiasPtr->p2pValue[ch][0]) >= fabs(maxInfo->maxp2p[ch][0]) )
          {
	      for (j=0;j<3;j++)
	      {
              maxInfo->maxp2p[ch][j]=peakperbiasPtr->p2pValue[ch][j];
              maxInfo->maxinx1[ch][j]=peakperbiasPtr->peakInx[ch][j];
              maxInfo->maxinx0[ch][j]=bias;

	      //init saOut Value is half way betweeen MAx-MIN or MIN-MAX
              a=peakperbiasPtr->peakPts[ch][j];
              b=peakperbiasPtr->peakPts[ch][j+1];          
              setup->maxValue.saOut[ch][j]=(a+b)/2;
            }
            maxInfo->maxinx1[ch][3]=peakperbiasPtr->peakInx[ch][3]; 
  
            // save here for SQ2SERVOSSALOCK, we are only interested in 
	    // zfact and gain, initFB will be picked in initlckPts() 
            //setup->maxValue.initFB[ch]=peakperbiasPtr->initFB[ch];
            setup->maxValue.zfact[ch]=peakperbiasPtr->zfact[ch];
            setup->maxValue.gain[ch]=peakperbiasPtr->gain[ch];
          }
        }
      }
    }
    // have 3 max peaks for each channel, copy them to setup structure   
    for(j=0; j<3; j++)
    {
      if (maxInfo->maxp2p[ch][0]!=0)
      {
        setup->maxValue.peakVal[ch][j]=maxInfo->maxp2p[ch][j];
 
        //ssabiasMax, sq2biasmax, sq1biasMax
        setup->maxValue.biasStep[ch][j]=maxInfo->maxinx0[ch][j];
 
        // ssafbMax, sq2fbMax or sq1fbMax
        setup->maxValue.fdbk[ch][j]=maxInfo->maxinx1[ch][j];
      }
      else
      {
        setup->maxValue.peakVal[ch][j] =NONMODULATION ;
        setup->maxValue.biasStep[ch][j]=NONMODULATION ;
        setup->maxValue.fdbk[ch][j]=NONMODULATION ;
      }
    }
    if (maxInfo->maxp2p[ch][0]!=0 )
    {
      setup->maxValue.fdbk[ch][3]=maxInfo->maxinx1[ch][3];
    }
    else
    {
      setup->maxValue.fdbk[ch][3]=NONMODULATION ;
    }
  }
}



/**
 * \fn void sc2dafindrefmaxModul(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     char *databuf, StatusType *status)
 *
 * \brief function
 *  get the reference points at REF_STEP times of max_modulation 
 *  here, we always use the last one of maxbiasNo
 *
 * \param  myInfo  dasInfoStruct_t pointer
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data 
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dafindrefmaxModul
*/
void sc2dafindrefmaxModul
(
dasInfoStruct_t *myInfo,
ARRAYSET       *setup,
char           *databuf,
StatusType     *status
)
{
  int        j,ch, refBias;
  int        maxBias;
  int        maxStep,a,b,*servodataPtr;
  int        *chdataperbiasPtr;
  char       *peakinfoPtr;
  BIAS_PEAK  *totalpeakPtr,*peakperbiasPtr;
  SERVO_DATAHEAD *rampinfoPtr;

  if (!StatusOkP(status)) return;

  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);
  peakinfoPtr  = (databuf+ servoInfo.headSize + servoInfo.totalframeSize);
   
  // point to the entry of array of BIAS_PEAK struture
  totalpeakPtr =(BIAS_PEAK *)peakinfoPtr;

  // point to max_modulation for each channel
  for (ch=0; ch<COL_NUM; ch++)
  {
    // first check
    setup->biasrefPt[ch]=NONMODULATION;         //sq1refbias
    setup->refValue.zfact[ch]=NONMODULATION;    //sq2reffb
    setup->refValue.initFB[ch]=NONMODULATION;   //sq1reffb

    if(setup->maxValue.peakVal[ch][0]!=NONMODULATION)
    {
      maxStep=setup->maxValue.biasStep[ch][0];
      //get Icmax 
      maxBias= maxStep*setup->stepBIAS + setup->minBIAS;
      // #define REF_STEP 0.2 get 1.2 Icmax as Icref
      refBias=(1+REF_STEP)*maxBias;

      // always use the last one as sq1biasRef
      setup->refValue.biasStep[ch][0]=setup->biasNo-1;
      setup->biasrefPt[ch]=
          setup->refValue.biasStep[ch][0]*setup->stepBIAS + setup->minBIAS;

      if(setup->biasrefPt[ch]== maxBias)
      {
         //we treat this as bad now, SO THAT in initlckPt
         //  setup->biaslckPt[ch]=NONMOD_BIAS;
         setup->maxValue.peakVal[ch][0] =NONMODULATION ;
      }
      else
      {
        // update the ref maxValue too. [0][1] is min,max peak
        peakperbiasPtr=totalpeakPtr + setup->refValue.biasStep[ch][0];
        for(j=0; j<3; j++)
        {
          setup->refValue.peakVal[ch][j]=peakperbiasPtr->p2pValue[ch][j];
          // sq1fbRef     
          setup->refValue.fdbk[ch][j]=peakperbiasPtr->peakInx[ch][j];

          a=peakperbiasPtr->peakPts[ch][j];
          b=peakperbiasPtr->peakPts[ch][j+1];          
          setup->refValue.saOut[ch][j]=(a+b)/2;
        }
        //  
        chdataperbiasPtr=servodataPtr +
                     setup->refValue.biasStep[ch][0]*servoInfo.totaldataperBias + 
                     ch;

        sc2dafindrefPoint(myInfo,setup,ch,chdataperbiasPtr,status);
      }
    }
  } 
}


/**
 * \fn void sc2dafindrefPoint(dasInfoStruct_t *myInfo, ARRAYSET *setup, int ch, 
 *      int *chdatarefbiasPtr, StatusType *status)
 *
 * \brief function
 *  find initial reference points        
 *
 * \param  myInfo         dasInfoStruct_t pointer
 * \param  setup          ARRAYSET structure pointer  G&R
 * \param  ch             int  which channel is looked at
 * \param  chdatarefbiasPtr  int pointer for max modulation ch  
 * \param  status         StatusType pointer  G&R
 *
 */
/*+ sc2dafindrefPoint
*/
void sc2dafindrefPoint
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup,
int        ch, 
int        *chdatarefbiasPtr,
StatusType *status
)
{
  int            j;
  int            fbStart,fbStop,data0,data1;
  int            halfPt;

  if (*status != STATUS__OK) return;

  halfPt=setup->refValue.saOut[ch][0];
  // assume there are only two peaks ( min; max)
  fbStart=setup->refValue.fdbk[ch][0];
  fbStop =setup->refValue.fdbk[ch][1];
  
  for (j=fbStart;j<fbStop;j++)
  {
    //start to look for saout~=halfPt
    // each chdatmaxbiasPtr has COL_NUM*setup->fdbkNo
    // each ch's data is at j*COL_NUM for fdbk
    data0=*(chdatarefbiasPtr+    j*COL_NUM);
    data1=*(chdatarefbiasPtr+(j+1)*COL_NUM);
 
    if(  (data0 <= halfPt && data1 > halfPt) ||    // positive
         (data0 > halfPt && data1 <= halfPt)       // negative
       )
    {
      // use refValue.zfact for sq2fbRef 
      //refValue.initFB      for sq1fbRef
      setup->refValue.zfact[ch]=data0;
      setup->refValue.initFB[ch]=j*setup->stepFDBK+setup->minFDBK;
      break;
    }
  }
}



//========================= sc2dag  =========================
//===========================================================



//========================= sc2dal  =========================
//===========================================================
/**
 * \fn void sc2daleastsquareFit(ARRAYSET *setup, char *databuf,
 *     StatusType *status)
 *
 * \brief function
 *  find least square fit coefficient
 * 
 * \param  setup    ARRAYSET structure pointer  
 * \param  databuf  char pointer for servo data
 * \param  status   StatusType pointer  G&R
 *
 */
/*+ sc2daleastsquareFit
*/
void sc2daleastsquareFit
(
ARRAYSET   *setup, 
char       *databuf,
StatusType *status
)
{

  int        bias,ch;
  long       frameSize, headSize, totalptsSize;
  char       *peakinfoPtr;
  double     *meanData,*meanX, *meanWeight, a;
  BIAS_PEAK  *peakPtr,*peakperbiasPtr=0;
  SERVO_DATAHEAD *rampinfoPtr;

  if (!StatusOkP(status)) return;

  // allocate memory to hold each ch's mean Value 
  meanData=(double *)calloc(setup->biasNo,sizeof(double));
  if (meanData==NULL)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
        "sc2daleastsquareFit: failed to allocate space for meanData");
    return;
  }
  meanX=(double *)calloc(setup->biasNo,sizeof(double));
  if (meanX==NULL)
  {
    free(meanData);
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
        "sc2daleastsquareFit: failed to allocate space for meanX");
    return;
  }
  meanWeight=(double *)calloc(setup->biasNo,sizeof(double));
  if (meanWeight==NULL)
  {
    free(meanData);
    free(meanX);
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
        "sc2daleastsquareFit: failed to allocate space for meanWeight");
    return;
  }

  totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
  headSize = sizeof(SERVO_DATAHEAD); 
  frameSize = totalptsSize*COL_NUM;

  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  peakinfoPtr  = (databuf+ headSize+frameSize);
   
  // point to the entry of array of BIAS_PEAK struture
  peakPtr =(BIAS_PEAK *)peakinfoPtr;

  for( ch=0;ch<COL_NUM; ch++)
  {
    // for each bias/cable setting, a set of mean Value
    for (bias=0; bias<setup->biasNo; bias++)
    {
      // point to the corresponding  BIAS_PEAK  for each bias
      peakperbiasPtr=peakPtr+bias;
      meanData[bias]=peakperbiasPtr->p2pValue[ch][0];
      meanX[bias]=setup->minBIAS+bias*setup->stepBIAS;
      meanWeight[bias]=1;
    }
    // dcoffset=setup->maxValue.peakVal[ch][2]
    // slope=setup->maxValue.peakVal[ch][1]
    sc2damath_linfit(setup->biasNo, meanX, meanData,meanWeight,
      &setup->maxValue.peakVal[ch][1],&setup->maxValue.peakVal[ch][2],status); 

    #ifdef PRINT_CH_DATA
    {
      if (ch==TEST_CH)
        printf("setup->maxValue.peakVal[%d][1]=%7.3f\n",
          ch,setup->maxValue.peakVal[ch][1]);    
    }
    #endif
	   
    if(setup->servo==CABLECAL)
    {
      // cableOffset (init) Value, we need to decide if
      // setup->maxValue.peakVal[ch][1] ==0, setup->maxValue.peakVal[ch][0]=?  
      // if SSA_BOTTOM_LINE=0, some times a>0, which results in cableOffset(int)<0
      // MCE dose not take that. we have to force it to 0 
      a= (-SSA_BOTTOM_LINE -setup->maxValue.peakVal[ch][2]);
      if ( a >=0)
        setup->maxValue.peakVal[ch][0]=0;
      else
      {
        if (setup->maxValue.peakVal[ch][1] >0 || 
            setup->maxValue.peakVal[ch][1] <0 )
          setup->maxValue.peakVal[ch][0]=a/setup->maxValue.peakVal[ch][1];
        else 
          setup->maxValue.peakVal[ch][0]=0;
      }
    }
    else if(setup->servo==SSARAMP1)
    {
      // cableScale Value
      if( setup->cableSlope[ch] !=0 )
      {
        setup->maxValue.peakVal[ch][0]=
            setup->maxValue.peakVal[ch][1]/setup->cableSlope[ch];
      }
      else 
        setup->maxValue.peakVal[ch][0]=0;

      #ifdef PRINT_CH_DATA
      {
        if (ch==TEST_CH)
          printf("setup->maxValue.peakVal[%d][0]=%7.3f cableSlope[%d]=%7.3f\n",
            ch,setup->maxValue.peakVal[ch][0],
            ch,setup->cableSlope[ch]);    
      }
      #endif
    }
  }    
  free(meanData);
  free(meanX);
  free(meanWeight);
}

/**
 * \fn void sc2dalinearConv(ARRAYSET *setup, double xData,
 *  double *filteData, StatusType *status)
 *
 * \brief function
 *  apply linear convolution and save result to filteData 
 * 
 * \param  setup    ARRAYSET structure pointer  
 * \param  xData    double 
 * \param  filteData double pointer 
 * \param  status   StatusType pointer  
 *
 * direct convolution to filte data, conArray[] shifting input
 * y(n)=convolution ( k=0, N-1) h(k)x(n-k)
 *
 * note: (N-1)/2 delay
 */
/*+ sc2dalinearConv
*/
void sc2dalinearConv
(
ARRAYSET   *setup,
double     xData,
double     *filteData,
StatusType *status
)
{
  int    n;
  double temp=0;

  if (*status != STATUS__OK) return;
  
  setup->convolPtr[0]=xData;
  for( n=0; n<FILTER_ORDER; n++)
    temp +=setup->impulsePtr[n]*setup->convolPtr[n];
  *filteData= temp;

  // shift for next data, from high inx to low inx
  for( n=FILTER_ORDER-1;  n >0; n--)
    setup->convolPtr[n]=setup->convolPtr[n-1];
}


/**
 * \fn void sc2dalookeachCh(ARRAYSET *setup, int *biasData,
 *     BIAS_PEAK *peakInfo,int whichCh, StatusType *status)
 *
 * \brief function
 *  find peak info for each channel 
 * 
 * \param  setup    ARRAYSET structure pointer  
 * \param  biasData data pointer for each bias level 
 * \param  peakInfo BIAS_PEAK poiter     
 * \param  whichCh  which channel is looked at
 * \param  status   StatusType pointer  
 *
 */
/*+ sc2dalookeachCh
*/
void sc2dalookeachCh
(
ARRAYSET   *setup, 
int     *biasData,
BIAS_PEAK  *peakInfo,
int         whichCh,
StatusType *status
)
{
  int        i,feedbk;
  double     *chData;
  double     deriVt,thredVal=DER_THDSSA;
  int     chkmaxFlag;

  if (!StatusOkP(status)) return;

  for (i=0;i<4;i++)
  {
    peakInfo->p2pValue[whichCh][i]=(double)0;
    peakInfo->peakInx[whichCh][i]=0;
    peakInfo->peakPts[whichCh][i]=0;
  }

  // allocate memory to hold each ch's data
  chData=(double *)calloc(setup->fdbkNo,sizeof(double));
  if (chData==NULL)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
        "sc2dalookeachCh: failed to allocate space for chData");
    return;
  }

  // get whichCH's data for the whole feed back ramping
  for (feedbk=0; feedbk<setup->fdbkNo; feedbk++)
  {
    chData[feedbk]=
      (double) biasData[whichCh+COL_NUM*feedbk]*setup->colMask[whichCh];
  }

  // If the column mask is zero there is no need to work with this column
  if (setup->colMask[whichCh]==0)
  {
    setup->nomodulCount[whichCh]=setup->biasNo;
    peakInfo->p2pValue[whichCh][0]=(double)NONMODULATION;
    peakInfo->p2pValue[whichCh][1]=(double)NONMODULATION;
    peakInfo->p2pValue[whichCh][2]=(double)NONMODULATION;

    /* Just to keep the size of the log file down I will only report this in SQ2LOCK */
    if(setup->servo == SQ2LOCK)
      {
	fprintf(myFpLog, "sc2dalookeachCh colMask caused NONMODULATION ch: %d\n",whichCh);
      }

    free(chData);
    return;
  }

  // For SQ2SERVO and SQ1SERVO check that the bias lock point is not marked as NONMOD_BIAS
  if (setup->servo==SQ2SERVO || setup->servo==SQ1SERVO)
  {
    if (setup->biaslckPt[whichCh]==NONMOD_BIAS)
    {
      setup->nomodulCount[whichCh]=setup->biasNo;
      peakInfo->p2pValue[whichCh][0]=(double)NONMODULATION;
      peakInfo->p2pValue[whichCh][1]=(double)NONMODULATION;
      peakInfo->p2pValue[whichCh][2]=(double)NONMODULATION;
      fprintf(myFpLog, "sc2dalookeachCh biaslckPt caused NONMODULATION ch: %d\n",whichCh);
      free(chData);
      return;
    }  
  }

  // check if ABS(max-min) >= PEAK_THD, chkmaxFlag will be returned as zero if it is not
  sc2dalookabsMax(chData,setup,whichCh,&chkmaxFlag,1,peakInfo,status);

  // If the ptop is OK, then check that the slope is also OK
  if ( chkmaxFlag != 0)
  {

    // check if there is a modulation, initial deriVt
    deriVt=0;
    sc2dachkDeriv(setup, chData, &deriVt, whichCh,1,status);
    if (!StatusOkP(status)) 
    {
      ErsRep(0,status,"sc2dalookeachCh: sc2dachkDeriv failed ");
      free(chData);
      return;
    }
  
    // Select the threshold to use based on the servo type
    if (setup->servo==SSARAMP || setup->servo==SSALOCK)
    {
      thredVal=DER_THDSSA;
    }
    else if (setup->servo==SQ2SERVO || setup->servo==SQ2LOCK)
      thredVal=DER_THDSQ2;
    else if (setup->servo==SQ1SERVO || setup->servo==SQ1LOCK)
      thredVal=DER_THDSQ1;
    else if (setup->servo==SQ1OPEN )
      thredVal=DER_THDSQ1OPEN;

    // now, get the peak-peak values (if the begining slope (deriVt) is greater than our threshold)
    // If slope is negative look for a minimum first, otherwise look for a maximum first
    //p2pValue[ch][0] peak0-peak1
    //p2pValue[ch][1] peak1-peak2
    //p2pValue[ch][2] peak2-peak3
    if ( deriVt > thredVal )
      sc2dafindmaxMin(chData,peakInfo,setup,whichCh,status);
    else if ( deriVt < (-thredVal) )
      sc2dafindminMax(chData,peakInfo,setup,whichCh,status);
    else
    {
      setup->nomodulCount[whichCh] ++;
      peakInfo->p2pValue[whichCh][0]=(double)NONMODULATION;
      peakInfo->p2pValue[whichCh][1]=(double)NONMODULATION;
      peakInfo->p2pValue[whichCh][2]=(double)NONMODULATION;
      fprintf(myFpLog, "sc2dalookeachCh: deriVt caused NONMODULATION ch: %d deriVt: %f thredVal: %f\n",
	      whichCh, deriVt, thredVal);
    }
  }
  else //chkmaxFlag == 0
  {
    setup->nomodulCount[whichCh] ++;
    peakInfo->p2pValue[whichCh][0]=(double)NONMODULATION;
    peakInfo->p2pValue[whichCh][1]=(double)NONMODULATION;
    peakInfo->p2pValue[whichCh][2]=(double)NONMODULATION;
    fprintf(myFpLog, "sc2dalookeachCh chkmaxFlag=0 caused NONMODULATION ch: %d\n",whichCh);
  }
  free(chData);
}


/**
 * \fn void sc2dalookabsMax(void *chdata, ARRAYSET *setup, 
 *  int ch, int *thdFlag, int flag,BIAS_PEAK *peakInfo, StatusType *status)
 *
 * \brief function
 *  find the abs (highest- lowest) in whole fdbkNo period, start from CHKSTART
 *  record them in peakInfo->lowest[ch], highest[ch]     
 * 
 * \param  chdata   void pointer for channel data 
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  ch       int which channel is looked at
 * \param  thdFlag  int pointer for abs(Max-min)>THD  G&R
 * \param  flag     int   0: int array, 1: double array
 * \param  peakInfo BIAS_PEAK pointer for the peak points/ch, G&R
 * \param  status   StatusType pointer.  G&R
 *
 */
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
)
{
  int      i,maxInx=0, minInx=0;
  int      peakTHD=PEAK_THDSSA;
  double   minVal, maxVal, data;
  int      *intArray;
  double   *floatArray;

  if (!StatusOkP(status)) return; 

  intArray=NULL;
  floatArray=NULL;

  // Select the threshold to use base on the type of servo
  if(setup->servo==SQ1SERVO || setup->servo==SQ1LOCK )
  {
    peakTHD=PEAK_THDSQ1;
  }
  else if(setup->servo==SQ1OPEN)
  {
    peakTHD=PEAK_THDSQ1OPEN;
  }
  else if ( setup->servo==SQ2SERVO || setup->servo==SQ2LOCK ||
            setup->servo==SQ2OPEN )
  {
    peakTHD=PEAK_THDSQ2;

  }
  else if  ( setup->servo==SSARAMP || setup->servo==SSALOCK )
  {
    peakTHD=PEAK_THDSSA;
    if  (setup->divFlag==1)
      peakTHD = PEAK_THDSSA * setup->sampleNo;
  }
  else if(setup->servo==SQ2OPEN4P)
  {
    peakTHD=PEAK_THDSQ2OPEN4P;
  }

  else
  {
    *status=DITS__APP_ERROR;
    sprintf(msg,"sc2dalookabsMax: should not do this for this servo %s ",setup->servoName);
    ErsRep(0,status,msg);
    return;
  }

  // flag tells us if we have double or int data
  if ( flag==1)
  { 
    floatArray=(double*)chdata;
    maxVal=floatArray[CHKSTART];
    minVal=floatArray[CHKSTART];
  }
  else
  { 
    intArray=(int*)chdata;
    maxVal=(double)intArray[CHKSTART];
    minVal=(double)intArray[CHKSTART];
  }

  // Go through the entire data set looking for the maximum and minimum values
  // Save the index of those values 
  for (i=CHKSTART; i<(setup->fdbkNo-1); i++)
  {
    if ( flag==1)
      data=floatArray[i];
    else
      data=(double)intArray[i];

    if( data >= maxVal) 
    {
      maxVal=data;
      maxInx=i;
    }
    else if( data <= minVal) 
    {
      minVal=data;
      minInx=i;
    }
  }
 
  // Check to see if the total overall modulation is greater than our threshold
  if ( fabs(maxVal-minVal) >=peakTHD )
    {
      *thdFlag=1;
      peakInfo->highest[ch]=(int) maxVal;
      peakInfo->lowest[ch]=(int) minVal;
      peakInfo->peakInx[ch][0]=minInx;
      peakInfo->peakInx[ch][1]=maxInx;
    }
  else // It is not greater than the threshold
    {
      // No need to give this message if the column was flatlined by something else previously
      if((maxVal + minVal) != 0.0)
	{
	  fprintf(myFpLog, "sc2dalookabsMax setting thdFlag to zero ch %d max: %f  at %d min: %f at %d thd: %d\n",
		  ch, maxVal, maxInx, minVal, minInx, peakTHD);
	}
      *thdFlag=0;
    }
}


/**
 * \fn void sc2dalookMax(double *chdata, BIAS_PEAK *peakInfo,
 *  ARRAYSET *setup, int ch, int *dataInx,int order,StatusType *status)
 *
 * \brief function
 *  Find the next maximum peak within the data.
 *  Scrolls down the data looking for a maximum. It stops
 *  when it finds a point that is greater than all of the following
 *  points within the search period. It saves the value of the maximum
 *  And returns the index of the maximum in dataInx so that the next
 *  search for a minimum can start there.
 *
 *         given and return ==>G&R        
 * 
 * \param  chdata   double pointer for channel data 
 * \param  peakInfo BIAS_PEAK pointer for the peak points/ch, G&R
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  ch       which channel is looked at
 * \param  dataInx  pointer for where to start the search inside data  G&R
 * \param  order    Where to store the results in the p2p data set
 * \param  status   StatusType pointer.  G&R
 *
 */
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
)
{
  int      i,search,searchPerd=SEARCHPERD;
  int      peakTHD=PEAK_THDSSA;
  double   minVal, absDiff;;

  if (!StatusOkP(status)) return;

  // Select a threshold based on the servo type we use
  // Also the search period (window length) we use is servo dependent
  // This is because some servos have different numbers of data points per modulation cycle
  // The search period is in number of data points and should be just less than
  // a modulation cycle for this method to work well
  // For example an SSA modulation period is is about 33-35 samples long so the SSA SERCHPERD is 30   
  if(setup->servo==SQ1SERVO || setup->servo==SQ1LOCK)
  {
    peakTHD=PEAK_THDSQ1;
    searchPerd=setup->fdbkNo/4+20;  //SEARCHPERD_SQ1;
  }
  else if ( setup->servo==SQ2SERVO || setup->servo==SQ2LOCK)
  {
    peakTHD=PEAK_THDSQ2;
    searchPerd=setup->fdbkNo/4+20;  //SEARCHPERD_SQ2;
  }
  else if  ( setup->servo==SSARAMP || setup->servo==SSALOCK )
  {
    peakTHD=PEAK_THDSSA;
    if (setup->divFlag==1)
      peakTHD = PEAK_THDSSA * setup->sampleNo;
    searchPerd=SEARCHPERD;
  }

  // dataInx is where to start the search. If this search is made for maximums 
  // later in the data set it will not be close to zero
  while(1)
  {
    // Set some defaults to test against
    peakInfo->peakPts[ch][order]=*(chdata + *dataInx);
    minVal=(double) peakInfo->lowest[ch];
    peakInfo->peakInx[ch][order]=*dataInx;
    search=0;

    // search from dataInx to the end of the data set
    for (i=*dataInx; i<(setup->fdbkNo-1); i++)
    {
      search++;
      // If search is > searchPerd we have been looking for other 
      // maximums but have not found them in searchPerd points
      // So call this our maximum
      if( search >searchPerd ) 
        break;
      
      // Is this point greater than the last we found?
      if( chdata[i] >= peakInfo->peakPts[ch][order] )
      {
        peakInfo->peakPts[ch][order]=chdata[i];
        peakInfo->peakInx[ch][order]=i;
        search=0;
      }
    }
    *dataInx=i;

    // If the difference between the maximum we found and a minimum we found 
    // earlier is greater than the threshold stop looking
    absDiff=abs(peakInfo->peakPts[ch][order]-minVal);
    if (  absDiff>peakTHD ||  i==(setup->fdbkNo-1)  )
      break;
  }

  // Set the search point for subsequent searches to the index of the maximum we found
  *dataInx=  peakInfo->peakInx[ch][order];
}


/**
 * \fn void sc2dalookmeaneachCh(ARRAYSET *setup, int *biasData,
 *     BIAS_PEAK *peakInfo,int whichCh, StatusType *status)
 *
 * \brief function
 *  find mean Value info for each channel 
 * 
 * \param  setup    ARRAYSET structure pointer  
 * \param  biasData data pointer for each bias level 
 * \param  peakInfo BIAS_PEAK poiter     
 * \param  whichCh  which channel is looked at
 * \param  status   StatusType pointer  
 *
 */
/*+ sc2dalookmeaneachCh
*/
void sc2dalookmeaneachCh
(
ARRAYSET   *setup, 
int     *biasData,
BIAS_PEAK  *peakInfo,
int         whichCh,
StatusType *status
)
{
  int        i,feedbk,iStart,cableAver;
  double     *chData;

  if (!StatusOkP(status)) return;

  // allocate memory to hold each ch's data
  chData=(double *)calloc(setup->fdbkNo,sizeof(double));
  if (chData==NULL)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
        "sc2dalookmeaneachCh: failed to allocate space for chData");
    return;
  }
  // get whichCH's data for the whole feed back ramping
  for (feedbk=0; feedbk<setup->fdbkNo; feedbk++)
  {
    chData[feedbk]=(double) biasData[whichCh+COL_NUM*feedbk];
  }

  // the CHKSTART is used to avoid the first few data as when I look at 
  // data plot, it is a bit weird at beginning 
  // use peakInfo to save mean value
  for (i=0;i<4;i++)
  {
    peakInfo->p2pValue[whichCh][i]=0.0;
    peakInfo->peakInx[whichCh][i]=0;
    peakInfo->peakPts[whichCh][i]=0;
  }
  
  if(  (setup->fdbkNo-CHKSTART) > (setup->fdbkNo*2)/5 )
  {
    cableAver=(setup->fdbkNo*2)/5;
    iStart=CHKSTART;
  }
  else if( setup->fdbkNo<=CHKSTART )
  {
    cableAver=setup->fdbkNo;
    iStart=0;
  }
  else 
  {
    cableAver=setup->fdbkNo-CHKSTART;
    iStart=CHKSTART;
  }  
  
  for (i=0; i<cableAver; i++)
  {
   peakInfo->p2pValue[whichCh][0] +=chData[iStart+i];
   //printf("%7.1f ",chData[iStart+i]);
  }
  peakInfo->p2pValue[whichCh][0] /= cableAver;
  free(chData);
}


/**
 * \fn void sc2dalookMin(double *chdata, BIAS_PEAK *peakInfo,
 *  ARRAYSET *setup, int ch, int *dataInx,int order,StatusType *status)
 *
 * \brief function
 *  Find the next minimum peak within the data.
 *  Scrolls down the data looking for a minimum. It stops
 *  when it finds a point that is less than than all of the following
 *  points within the search period. It saves the value of the minimum
 *  and returns the index of the minimum in dataInx so that the next
 *  search for a maximum can start there.
 *
 *         given and return ==>G&R        
 * 
 * \param  chdata   double pointer for channel data 
 * \param  peakInfo BIAS_PEAK pointer for the peak points/ch, G&R
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  ch       which channel is looked at
 * \param  dataInx  pointer for where to start the search inside data  G&R
 * \param  order    Where to store the results in the p2p data set
 * \param  status   StatusType pointer.  G&R
 *
 */
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
)
{
  int  i,searchPerd=SEARCHPERD;
  int  search;

  if (!StatusOkP(status)) return;

  // Select the search period (window length) we used based on the servo type
  // This is because some servos have different numbers of data points per modulation cycle
  // The search period is in number of data points and should be just less than
  // a modulation cycle for this method to work well
  // For example an SSA modulation period is is about 33-35 samples long so the SSA SERCHPERD is 30   

  if(setup->servo==SQ1SERVO || setup->servo==SQ1LOCK)
  {
    searchPerd=setup->fdbkNo/4+20;  //SEARCHPERD_SQ1;
  }
  else if ( setup->servo==SQ2SERVO || setup->servo==SQ2LOCK)
  {
    searchPerd=SEARCHPERD_SQ2;
  }
  else if  ( setup->servo==SSARAMP || setup->servo==SSALOCK )
  {
    searchPerd=SEARCHPERD;
  }

  // Initialize before starting the search
  peakInfo->peakPts[ch][order]=*(chdata + *dataInx);
  peakInfo->peakInx[ch][order]=*dataInx;
  search=0;

  // Search from dataInx to end of data for a minimum
  for (i=*dataInx; i<(setup->fdbkNo-1); i++)
  {
    search++;

    // If search is > searchPerd we have been looking for other 
    // minimums but have not found them in searchPerd points
    // So call this our minimum
    if( search >searchPerd )
      break;

    // Is this a minimum?
    if( chdata[i] <= peakInfo->peakPts[ch][order])
    {
      peakInfo->peakPts[ch][order]=chdata[i];
      peakInfo->peakInx[ch][order]=i;
      search=0;
    }
  }

  // Set dataInx to the index where we found the minimum so that the next search for a maximum start there
  *dataInx=  peakInfo->peakInx[ch][order];
}


//========================= sc2dar  =========================
//===========================================================

/**
 * \fn void sc2dareadbinaryfile2Mem(dasInfoStruct_t *myInfo,
 *   FILE *fp, char **memPtr, int *fileLen,StatusType *status)
 *
 * \brief function:
 *   read a binary file into memory, return the memory pointer
 *
 * \param myInfo     dasInfo structure pointer
 * \param fp          FILE pointer
 * \param memPtr      double char pointer, for file memory address
 * \param fileLen     int pointer
 * \param *status       StatusType.  given and return
 *
 */
/*+ sc2dareadbinaryfile2Mem*/
void sc2dareadbinaryfile2Mem
(
dasInfoStruct_t *myInfo,      
FILE            *fp, 
char            **memPtr,       
int             *fileLen,             
StatusType      *status
)
{
  int   fileLength,dataSize;  
  char  *memAddr;

  if (*status != STATUS__OK) return;

  // find out how big the file is 
  if (fseek(fp, SEEK_SET, SEEK_END)) 
  {
    *status = DITS__APP_ERROR;
    sprintf(msg,"sc2dareadbinaryfile2Mem: Strange problems with this file");
    ErsRep(0,status,msg);
    return;
  }
  fileLength = ftell(fp);
  if (fileLength < 0 ) 
  {    
    // -1L is error return 
    *status = DITS__APP_ERROR;
    sprintf(msg,"sc2dareadbinaryfile2Mem:Error getting file length (%d) ",fileLength);
    ErsRep(0,status,msg);
    return;
  }
   *fileLen=fileLength;

  // fseek to the beginning of the file 
  if (fseek(fp, 0, SEEK_SET)) 
  {
    *status = DITS__APP_ERROR;
    sprintf(msg,"sc2dareadbinaryfile2Mem: Error- Strange problems with this file");
    ErsRep(0,status,msg);
    return;
  }
  // get the Allocate storage for the file 
  if (  (memAddr = (char*)malloc(fileLength)) == NULL) 
  {
    *status = DITS__APP_ERROR;
    sprintf(msg,"sc2dareadbinaryfile2Mem: Insufficient memory for storing this file");
    ErsRep(0,status,msg);
    return;
  }
  // read the file into memory 
  dataSize = fread(memAddr, 1, fileLength, fp);
  if (dataSize== 0) 
  {
    *status = DITS__APP_ERROR;
    sprintf(msg,"sc2dareadbinaryfile2Mem: Unable to read the file");
    ErsRep(0,status,msg);
    free(memAddr); 
    return;
  }
  // pass the memory pointer back
  *memPtr=memAddr;
}


//========================= sc2das  =========================
//===========================================================

/**
 * \fn void sc2dasaveservoData(ARRAYSET *setup, char *databuf, 
 * char *fileName,StatusType *status)
 *
 * \brief function
 *  save the servo data in binary 
 *         given and return ==>G&R        
 *
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  databuf  char pointer for servo data  
 * \param  fileName char pointer for first part of servodata.hex  
 * \param  status   StatusType pointer  G&R
 *
 */
/*+ sc2dasaveservoData
*/
void sc2dasaveservoData
(
ARRAYSET   *setup, 
char       *databuf,
char       *fileName,
StatusType *status
)
{
  int           frameSize, headSize, totalptsSize;
  int           peakSize, maxpeakSize, totalSize;
  char           tmpfile[FILE_LEN];
  FILE           *serdataFp;

  if (*status != STATUS__OK) return;

  totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
  headSize = sizeof(SERVO_DATAHEAD); 
  frameSize = totalptsSize*COL_NUM;
  peakSize = sizeof(BIAS_PEAK)*setup->biasNo;
  maxpeakSize = sizeof(MAX_P2P);
  totalSize = headSize + frameSize + peakSize + maxpeakSize;

  sprintf (tmpfile, "%sservodata.hex",fileName);
  if((serdataFp = fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasaveservoData: Error- failed to open file %s", tmpfile); 
      return;
    }
  fwrite(databuf, 1, totalSize,serdataFp);
  fclose(serdataFp);
}


/**
 * \fn void sc2dasavelockPts(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, char *databuf, FILE *fp, StatusType *status)
 *
 * \brief function
 *  save all points searched 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data (max COL*(ROW-1))
 * \param  fp      FILE pointer for result file
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavelockPts
*/
void sc2dasavelockPts
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
)
{
  int            i,ch, avefbPeriod[COL_NUM];
  int            *servodataPtr;
  char           *peakinfoPtr;  
  int             *cableadjPtr;
  BIAS_PEAK      *totalpeakPtr,*peakperbiasPtr;
  SERVO_DATAHEAD *rampinfoPtr;


  if (!StatusOkP(status)) return;

  jitDebug(16,"sc2dasavelockPts: dasInfo.mcecmdFile %s\n",myInfo->cmdrepFile);
 
  // assign different pointers
  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);
  peakinfoPtr  = (databuf+ servoInfo.headSize + servoInfo.totalframeSize);

 // point to arry of biasNo for BIAS_PEAK struture
  totalpeakPtr =(BIAS_PEAK *)peakinfoPtr;
  
  if( setup->servo!=CABLECAL && setup->servo!=SSARAMP1)
  { 
    for (i=0; i<setup->biasNo; i++)
    {
     // move to the corresponding  BIAS_PEAK  for each bias
      peakperbiasPtr=totalpeakPtr+i;
      fprintf(fp,"bias=%d \n",i*setup->stepBIAS+setup->minBIAS);
      for (ch=0;ch<COL_NUM; ch++)
      {  
        fprintf(fp,"ch_%d ",ch);
	 
        fprintf(fp,"fdVal[%7d %7d %7d %7d] ", 
         peakperbiasPtr->peakInx[ch][0]*setup->stepFDBK+setup->minFDBK,
	 peakperbiasPtr->peakInx[ch][1]*setup->stepFDBK+setup->minFDBK, 
         peakperbiasPtr->peakInx[ch][2]*setup->stepFDBK+setup->minFDBK,
	 peakperbiasPtr->peakInx[ch][3]*setup->stepFDBK+setup->minFDBK);
	
        fprintf(fp,"pkVal[%9.1f %9.1f %9.1f %9.1f] ",
         peakperbiasPtr->peakPts[ch][0],peakperbiasPtr->peakPts[ch][1], 
         peakperbiasPtr->peakPts[ch][2],peakperbiasPtr->peakPts[ch][3]);
    
        fprintf(fp,"p2pVal[%9.1f %9.1f %9.1f]\n", 
         peakperbiasPtr->p2pValue[ch][0],peakperbiasPtr->p2pValue[ch][1], 
         peakperbiasPtr->p2pValue[ch][2]);
      }
    }
    fprintf(fp,"\n");
    for (ch=0;ch<COL_NUM; ch++)
    {
      fprintf(fp,"max_Ch_%d ",ch);
      if(setup->maxValue.peakVal[ch][0]==NONMODULATION)
        fprintf(fp,"pkVal[non modulation]\n");
      else
      { 
        if(setup->servo==SSALOCK || setup->servo==SQ2OPEN ||
         setup->servo==SQ1OPEN || setup->servo == SQ2LOCK   )
        {  
          fprintf(fp,"biasVal[%7d ] ",   setup->biaslckPt[ch]);
        }
        else
        { 
          fprintf(fp,"biasVal[%7d %7d %7d] ",
            setup->maxValue.biasStep[ch][0]*setup->stepBIAS+setup->minBIAS,
            setup->maxValue.biasStep[ch][1]*setup->stepBIAS+setup->minBIAS,
            setup->maxValue.biasStep[ch][2]*setup->stepBIAS+setup->minBIAS);
        }	     
        fprintf(fp,"fdVal[%7d %7d %7d %7d] ",
          setup->maxValue.fdbk[ch][0]*setup->stepFDBK+setup->minFDBK,
          setup->maxValue.fdbk[ch][1]*setup->stepFDBK+setup->minFDBK,
          setup->maxValue.fdbk[ch][2]*setup->stepFDBK+setup->minFDBK,
          setup->maxValue.fdbk[ch][3]*setup->stepFDBK+setup->minFDBK);

       if (setup->servo ==SSARAMP || setup->servo ==SSALOCK)
       {  
          avefbPeriod[ch]= setup->maxValue.fdbk[ch][3]*setup->stepFDBK-
                       setup->maxValue.fdbk[ch][1]*setup->stepFDBK;
          avefbPeriod[ch]+= setup->maxValue.fdbk[ch][2]*setup->stepFDBK-
                         setup->maxValue.fdbk[ch][0]*setup->stepFDBK;
          avefbPeriod[ch] /=2;
       }
       else if (setup->servo ==SQ2SERVO || setup->servo ==SQ2LOCK)
       {
         avefbPeriod[ch]=setup->fluxPeriod[ch];
       }

       else
       {  
         avefbPeriod[ch]= setup->maxValue.fdbk[ch][1]*setup->stepFDBK-
                          setup->maxValue.fdbk[ch][0]*setup->stepFDBK;
          avefbPeriod[ch] *=2;
        }
        fprintf(fp,"avePeriod[%7d] ",avefbPeriod[ch]);

        fprintf(fp,"p2pVal[%9.1f %9.1f %9.1f]\n",
             setup->maxValue.peakVal[ch][0],
	     setup->maxValue.peakVal[ch][1],
	     setup->maxValue.peakVal[ch][2]);
      }
    }
    fprintf(fp,"\n");
    for (ch=0;ch<COL_NUM; ch++)
    {
      if(setup->maxValue.peakVal[ch][0]!=NONMODULATION)
      {
        fprintf(fp,"initLock for ch_%d ",ch);
        fprintf(fp,"biasVal[%7d]",setup->biaslckPt[ch]);  
        fprintf(fp,"fdVal[%7d] ",setup->initFB[ch]);  
        if(setup->servo==SQ1SERVO || setup->servo==SQ1LOCK )    
          fprintf(fp,"sq2fb[%7d] ",setup->saoutlckVal[ch]);
        else if(setup->servo==SQ2SERVO || setup->servo==SQ2LOCK)    
          fprintf(fp,"safb[%7d] Flxuperiod[%7d] ",setup->saoutlckVal[ch],
                  setup->fluxPeriod[ch]);    
        else    
          fprintf(fp,"saOut[%7d] ",setup->saoutlckVal[ch]);

        fprintf(fp,"gain[%9.4f]\n", setup->gain[ch]);
      }
    }

    fprintf(fp,"\n# BIASLCK_N: non-modulation ones, set bias to NONMOD_BIAS (%d)",
             NONMOD_BIAS);
    if ( setup->servo==SSARAMP || setup->servo==SSALOCK)
    {
      fprintf(fp,"\n# ssabias Val at max modulation");
      cableadjPtr=(int*)setup->cableadjPtr;

      // for simple test 
      for (i=0;i<setup->totalRC2use;i++)
      {
        fprintf (fp, "\n# wb rc%d sa_bias",setup->whichRC[i] ); 
        for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
          fprintf(fp,"%7d ", setup->biaslckPt[ch]);
        fprintf (fp, "\n# wb rc%d  offset ",setup->whichRC[i] ); 
        for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
         fprintf( fp,"%7d ", 
                   setup->cableOffset[ch]
                   -(int)(setup->biaslckPt[ch]*setup->cableScale[ch])
                   + *(cableadjPtr+ch +setup->maxValue.biasStep[ch][0]*COL_NUM) 
                 );
      }
    }
    else if ( setup->servo==SQ2SERVO || setup->servo==SQ2OPEN ||
              setup->servo==SQ2LOCK || setup->servo==SQ2OPEN4P)
    {
      fprintf(fp,"\n# sq2bias Val at optimal point");
      // SQ2 bias 
      fprintf (fp,"\n# wb %s flux_fb <\n",setup->sq2biasCard );
      sc2daservolckPoints(setup,fp," ", setup->biaslckPt,1,status); 
    }
    else if ( setup->servo==SQ1SERVO || setup->servo==SQ1OPEN || 
              setup->servo==SQ1LOCK 
            )
    {
      fprintf(fp,"\n# sq1bias Val at optimal point");
      fprintf(fp,"\n# wb %s on_bias <\n",setup->sq1biasCard ); 
      sc2daservolckPoints(setup,fp," ", setup->biaslckPt,1,status);
    }
    else
    {
    }

    fprintf(fp,"\n");
    sc2daservolckPoints1(setup,fp,"biaslckPt=",setup->biaslckPt,status);

    fprintf(fp,"\n# INITFB_N   non-modulation ones, set intFB to (%d)\n",
             NONMOD_BIAS);
      if ( setup->servo==SSARAMP || setup->servo==SSALOCK)
    {
      fprintf(fp,"# initial ssafb lock Value\n");
      fprintf(fp,"# wb %s flux_fb <\n", setup->safbCard);
    }
    else if ( setup->servo==SQ2SERVO || setup->servo==SQ2OPEN ||
              setup->servo==SQ2LOCK)
    {
      fprintf(fp,"# initial sq2fb lock Value\n");
      fprintf(fp,"# wb %s flux_fb <\n",setup->sq2fbCard ); 
    }
    else
    {
    }
    // for simple test 
    sc2daservolckPoints(setup,fp," ",setup->initFB,1,status);
    sc2daservolckPoints1(setup,fp,"intFB= ",setup->initFB,status);

    fprintf(fp,"\n# ZFACTOR_N:  ");
    if (setup->servo==SSARAMP || setup->servo==SSALOCK)
    {
      fprintf(fp," ssaDAOut Value \n");
    }
    else if (setup->servo==SQ2SERVO || setup->servo==SQ2LOCK)
    {
      fprintf(fp," for SQ2SERVO/LOCK, this is ssafb value at lock points used in SQ2OPEN\n");
      fprintf(fp, "# wb %s flux_fb <\n" , setup->safbCard);
    }
    else if (setup->servo==SQ1SERVO || setup->servo==SQ1LOCK)
    {
      fprintf(fp,"# for SQ1SERVO/LOCK, this is sq2fb value at lock points used in SQ1OPEN\n");
    }
    else
    {
      fprintf(fp,"  \n");
    }

    if (setup->servo==SQ2LOCK)
    {
      sc2daservolckPoints(setup,fp," ",setup->saoutlckVal,1,status);
      sc2daservolckPoints1(setup,fp,"zFactor= ",setup->saoutlckVal,status);
      fprintf(fp,"\n# GAIN_N for all channels\n");
      sc2daservolckgainPoints1(setup,fp,"Gain= ",setup->gain,status);
      fprintf(fp,"\n");
    }

    if (setup->servo==SQ2SERVO)
    {
      sc2daservolckPoints(setup,fp," ",setup->maxValue.zfact,1,status);
      sc2daservolckPoints1(setup,fp,"zFactor= ",setup->maxValue.zfact,status);

      fprintf(fp,"\n# GAIN_N for all channels\n");
      sc2daservolckgainPoints1(setup,fp,"Gain= ",setup->maxValue.gain,status);

      sc2daservolckPoints(setup,fp,
             "for SQ2SERVO, ssafb value at lock points from sq2ssalock",
             setup->maxValue.initFB,1,status);
      sc2daservolckPoints(setup,fp,
             "for SQ2SERVO, ssaout value at lock points from sq2ssalock",
             setup->maxValue.zfact,1,status);
      sc2daservolckgainPoints(setup,fp,
             "for SQ2SERVO, gain value at lock points from sq2ssalock",
             setup->maxValue.gain,1,status);
    }
    if ( setup->servo==SQ1SERVO || setup->servo==SQ1LOCK )
    {
      // sq1baisMax=setup->biaslckPt[ch]
      // sq2fbMax=setup->maxValue.zfact[ch];
      // sq1fbMax=setup->maxValue.initFB[ch]
      // sq1baisRef=setup->biasrefPt[ch]
      // sq2fbRef=setup->revValue..zfact[ch]; 
      // sq1fbRef=setup->refValue.initFB[ch];

      fprintf(fp,"\n# SQ1SERVO1 initial optimal and reference points\n");
      sc2dasq1maxrefPoints(setup,fp,"sq1bMax[]=",setup->biaslckPt,status);
      sc2dasq1maxrefPoints(setup,fp,"sq1fbMax[]=",setup->maxValue.initFB,status);
      sc2dasq1maxrefPoints(setup,fp,"sq2fbMax[]=",setup->maxValue.zfact,status);
      sc2dasq1maxrefPoints(setup,fp,"sq1bRef[]=",setup->biasrefPt,status);
      sc2dasq1maxrefPoints(setup,fp,"sq1fbRef[]=",setup->refValue.initFB,status);
      sc2dasq1maxrefPoints(setup,fp,"sq2fbRef[]=",setup->refValue.zfact,status);
    }
  }
  else
  {
    for (i=0; i<setup->biasNo; i++)
    {
     // move to the corresponding  BIAS_PEAK  for each bias
      peakperbiasPtr=totalpeakPtr+i;
      if ( setup->servo ==CABLECAL)
         fprintf(fp,"\ncableoffset=");
      else
         fprintf(fp,"\nssabias=");

      fprintf(fp,"%06d ",i*setup->stepBIAS+setup->minBIAS);
      for (ch=0;ch<COL_NUM; ch++)
      {  
        fprintf(fp," %4.3E",peakperbiasPtr->p2pValue[ch][0]*setup->colMask[ch]);
      }
    }
    fprintf(fp,"\ncableSlope=0000     ");
    for (ch=0;ch<COL_NUM; ch++)
      fprintf(fp," %4.3E",setup->maxValue.peakVal[ch][1]*setup->colMask[ch]);

    fprintf(fp,"\noffsetConstant=0000 ");
    for (ch=0;ch<COL_NUM; ch++)
      fprintf(fp," %4.3E",setup->maxValue.peakVal[ch][2]*setup->colMask[ch]);
 
    
    if( setup->servo ==CABLECAL)
      fprintf(fp,"\nCableOffset(init)=0");
    else 
      fprintf(fp,"\ncableScale=0000     ");
      
    for (ch=0;ch<COL_NUM; ch++)
        fprintf(fp," %4.3E",setup->maxValue.peakVal[ch][0]*setup->colMask[ch]);
    
    fprintf(fp,"\n");
  }
}


/**
 * \fn void sc2dasaveframelockPts(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save all pixel points from sq1open 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasaveframelockPts
*/
void sc2dasaveframelockPts
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char tmpfile[FILE_LEN];
  FILE *fpgain;

  if (*status != STATUS__OK) return;
 
  // Call the routine that writes out the bad pixel map and sets ival to zero if greater than maxGainI 
  sc2dalibsetup_writepixelInx(myInfo,setup,myInfo->dataFile,setup->pixeLock.iVal,
                     1,status);

  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ),MCEGAIN);
  if((fpgain = fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasaveframelockPts: Error- failed to open file %s", tmpfile); 
      return;
    }

  // 0: int array, 1 double array
  fprintf(fpgain,"\n# I-VAL for all pixels, it is (-1/slope *4096) \n");
  sc2daservopixellckpts4MCE(myInfo,setup, fpgain,0, "I_VAL",setup->pixeLock.iVal,status);

  fprintf(fpgain,"\n# LOCK_N for all pixels, it is normalised by samplNo\n");
  sc2daservopixellckpts4MCE(myInfo,setup,fpgain,0, "ADC_OFFSET",setup->pixeLock.ssalockPt,status);

  fprintf(fpgain,"\n# FLUX period for all pixels \n");
  sc2daservopixellckpts4MCE(myInfo,setup,fpgain,0, "FLUXPERIOD",setup->pixeLock.fluxPeriod,status);

  fprintf(fpgain,"\n# SQ1FDBK_N for all pixels \n");
  sc2daservopixellckpts4MCE(myInfo,setup,fpgain,0, "SQ1FDBK",setup->pixeLock.sq1FB,status);

  fclose(fpgain);

  // still keep original one
  sprintf (tmpfile, "%s/%s-org", getenv ( "SC2SCRATCH" ),MCEGAIN);
  if((fpgain = fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasaveframelockPts 2: Error- failed to open file %s", tmpfile); 
      return;
    }

  // sc2daservopixellckPoints() 0: int array, 1 double array
  fprintf(fpgain,"\n# GAIN_N for all pixels, it is normalised by samplNo\n");
  sc2daservopixellckPoints(myInfo,setup, fpgain,1,setup->pixeLock.gain,status);

  fprintf(fpgain,"\n\n# IVAL_N for all pixels, it is normalised by samplNo\n");
  sc2daservopixellckPoints(myInfo,setup, fpgain,0,setup->pixeLock.iVal,status);

  fprintf(fpgain,"\n# LOCK_N for all pixels, it is normalised by samplNo\n");
  sc2daservopixellckPoints(myInfo,setup,fpgain,0,setup->pixeLock.ssalockPt,status);

  fprintf(fpgain,"\n# SQ1FDBK_N for all pixels \n");
  sc2daservopixellckPoints(myInfo,setup,fpgain,0,setup->pixeLock.sq1FB,status);

  fprintf(fpgain,"\n# FLUX period for all pixels \n");
  sc2daservopixellckPoints(myInfo,setup,fpgain,0,setup->pixeLock.fluxPeriod,status);

  fclose(fpgain);
}

/**
 * \fn void sc2daservocopydataBuf(dasInfoStruct_t *myInfo, 
 *   ARRAYSET *setup, char *dataBuf, StatusType *status)
 *
 * \brief function 
 *  tmp, copy dataBuf to servodataPtr
 *
 * \param myInfo    dasInfo structure pointer
 * \param setup   ARRAYSET structure pointer
 * \param dataBuf  char pointer for servo data
 * \param status     StatusType     
 *
 * 
 *  sc2daservocopyData(&dasInfo,&setup,dataBuf,status); 
 */
/*+ sc2daservocopydataBuf
*/
void sc2daservocopydataBuf
(
dasInfoStruct_t       *myInfo,    
ARRAYSET              *setup,
char                  *dataBuf,
StatusType            *status
)
{
  int        frameSize, headSize, totalptsSize;
  int        peakSize, maxpeakSize, subtotalSize,totalSize;
  int        cableadjSize, meanvalSize;
  
  if (!StatusOkP(status)) return;
  
  if( setup->servo==CABLECAL)
  {
    totalptsSize= setup->cableNo*setup->fdbkNo*sizeof(int);
    peakSize = sizeof(BIAS_PEAK)*setup->cableNo;
  }
  else 
  {
    totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
    peakSize = sizeof(BIAS_PEAK)*setup->biasNo;
  }
  headSize = sizeof(SERVO_DATAHEAD); 
  frameSize = totalptsSize*COL_NUM;
  maxpeakSize = sizeof(MAX_P2P);
  subtotalSize=  headSize + frameSize + peakSize + maxpeakSize;

  cableadjSize =setup->biasNo*COL_NUM*sizeof(int);
  meanvalSize =  cableadjSize;
  totalSize =  subtotalSize + cableadjSize + meanvalSize;
  memcpy(setup->servodataPtr,dataBuf,totalSize);

  // copy headinfo into filtedData so that findlock and other function calls
  // can get correct info  
  memcpy(setup->filtedPtr,setup->servodataPtr,headSize);      
}


/**
 * \fn void sc2daservogetfiltedData(dasInfoStruct_t *myInfo, 
 *   ARRAYSET *setup, FILE *fp, StatusType *status)
 *
 * \brief function 
 *  apply LP filter to servo data and save the filted data
 *
 * \param myInfo    dasInfo structure pointer
 * \param setup   ARRAYSET structure pointer
 * \param fp       file pointer for servo/filtered data
 * \param status     StatusType     
 *
 */
/*+ sc2daservogetfiltedData
*/
void sc2daservogetfiltedData
(
dasInfoStruct_t       *myInfo,    
ARRAYSET              *setup,
FILE                  *fp,
StatusType            *status
)
{
  int        bias, ch;
  int        totaldataperBias;
  long       frameSize, headSize, totalptsSize;
  long       peakSize, maxpeakSize, totalSize;

  int            *servodataPtr, *tmpData;
  SERVO_DATAHEAD *rampinfoPtr;

  int            *filterservodataPtr, *filtertmpData;
  SERVO_DATAHEAD *filterrampinfoPtr;

  if (!StatusOkP(status)) return;
  
  // move away from infomation head and get the right entry pointer
  // dataBuf: SERVO_DATAHEAD string
  //          data of all Channels

  totaldataperBias=COL_NUM*setup->fdbkNo;
  totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
  headSize = sizeof(SERVO_DATAHEAD); 
  frameSize = totalptsSize*COL_NUM;
  peakSize = sizeof(BIAS_PEAK)*setup->biasNo;
  maxpeakSize = sizeof(MAX_P2P);
  totalSize = headSize + frameSize + peakSize + maxpeakSize;

  rampinfoPtr  = (SERVO_DATAHEAD *)setup->servodataPtr;
  servodataPtr = (int *)(rampinfoPtr +1);

  filterrampinfoPtr  = (SERVO_DATAHEAD *)setup->filtedPtr;
  filterservodataPtr = (int *)(filterrampinfoPtr +1);


  // for each bias setting, a set of data are taken at different
  // feedbk points
  for (bias=0; bias<setup->biasNo; bias++)
  {
    // point to the first data for each bias
    tmpData=servodataPtr + bias*totaldataperBias;  
    
    // point to the first filted data for each bias
    filtertmpData=filterservodataPtr + bias*totaldataperBias;  
    
    for( ch=0;ch<COL_NUM; ch++)
    {
      // apply filter
      sc2daapplyFilter(myInfo,setup,ch,tmpData,filtertmpData,status);

      if(!StatusOkP(status)) 
      {
        //ErsRep(0,status,"sc2dafindlockPts: failed to call sc2dalookeachCh");
        printf("servogetfiltedData: failed to call sc2daapplyfilter \n");
        return ;
      }
    }
    sc2daservosaveFilted(myInfo,setup,tmpData,filtertmpData,bias,fp,status);   
  }
}

/**
 * \fn void sc2daservolckPoints(ARRAYSET *setup, FILE *fp, char *name,
 *  int *paramarray, int flag, StatusType *status)
 *
 * \brief function:
 *  save servo lock points for each channel
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray int pointer for setup-parameter array
 * \param  flag       int  1:ad # for comment, 0: no #
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  int   ch,i;

  if (*status != STATUS__OK) return;

  if( strcmp(name," ") !=0)
      fprintf(fp,"# %s\n",name);

  for (i=0;i<4;i++)
  {
    if (flag==1)
      fprintf(fp,"#");
    for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
      fprintf(fp,"%7d ", paramarray[ch]*setup->colMask[ch]);
    if(i<3)
      fprintf(fp,"<\n");
    else
      fprintf(fp,"\n");
   }
}



/**
 * \fn void sc2dasavessalckDivbyn(ARRAYSET *setup, FILE *fp, char *name,
 *  int *paramarray, int flag, StatusType *status)
 *
 * \brief function:
 *  save ssa lock points divided by sampleNo for MCE firmware
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray int pointer for setup-parameter array
 * \param  flag       int  1:ad # for comment, 0: no #
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  int   ch,i;

  if (*status != STATUS__OK) return;

  if( strcmp(name," ") !=0)
      fprintf(fp,"# %s\n",name);

  for (i=0;i<4;i++)
  {
    if (flag==1)
      fprintf(fp,"#");
    for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
      fprintf(fp,"%7d ", (paramarray[ch]*setup->colMask[ch])/setup->sampleNo);
    if(i<3)
      fprintf(fp,"<\n");
    else
      fprintf(fp,"\n");
   }
}


/**
 * \fn void sc2daservofindfluxPeriod(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     int ch, int *chdatamaxbiasPtr, StatusType *status)
 *
 * \brief function
 *  find fluxperiod from ssafb(sq2fb) curve, for sq2servo
 *  from SG test, it seems that sq2fb has to be 0 to 65535 to have
 * a full period of ssafb change 
 *         given and return ==>G&R        
 *
 * \param  myInfo         dasInfoStruct_t pointer
 * \param  setup          ARRAYSET structure pointer  G&R
 * \param  ch             int  which channel is looked at
 * \param  chdatamaxbiasPtr  int pointer for max modulation ch  
 * \param  status         StatusType pointer  G&R
 *
 */
/*+ sc2daservofindfluxPeriod
*/
void sc2daservofindfluxPeriod
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int        ch, 
int        *chdatamaxbiasPtr,
StatusType *status
)
{
  int    j,i,l,fstEnd,flag=0,startPt=5;
  int    fbStart,fbMiddle,fbEnd;
  int     data0,data2=0, data3=0,data4=0, data5=0;
  int     ssafitWind=10;     
  double  yData[ssafitWind],xData[ssafitWind],weightData[ssafitWind];
  double  slope, dclevel;

  if (*status != STATUS__OK) return;
  
  fbStart=startPt;
  while (flag==0)
  {
    // we start from fifth point from ssafb(sq2fb)
    data0=*(chdatamaxbiasPtr+ (fbStart)*COL_NUM);

    for (i=0; i<ssafitWind; i++)
    {
      xData[i]=i+1;
      weightData[i]=1;
      yData[i]=*(chdatamaxbiasPtr+ (fbStart+i)*COL_NUM);
     }
     sc2damath_linfit(ssafitWind, xData, yData,weightData, &slope,&dclevel,status); 

    //start to look for ssafb >=data0 or <= data0 depending on
    // the slop  
    // first half period
    fstEnd=setup->fdbkNo-2;
    fbMiddle=fbStart;
  
    for (j=(fbStart+10);j<fstEnd;j++)
    {
       l=j+5;
      // each chdatmaxbiasPtr has COL_NUM*setup->fdbkNo
      // each ch's data is at j*COL_NUM for fdbk
      data2=*(chdatamaxbiasPtr+    j*COL_NUM);
      if ( (l+1)== fstEnd )
        l=fstEnd;
	
      data3=*(chdatamaxbiasPtr+    (l)*COL_NUM);
    
      if ( slope > 0 )
      {
        if( (data2 <= data0)  &&  ( data3 <= data0)) 
        { 
          fbMiddle=j;
	  flag=1;
          //DEBFLUX("slope(%f) [ch=%d] data2[%d] data3[%d]\n", slope,ch, data2,data3);
          break;
        }
      }
      else //if (slope < 0 )
      {
        if ( (data2 >= data0) && ( data3 >= data0) )
        {
          fbMiddle=j;
	  flag=1;
	  //DEBFLUX("slope(%f) [ch=%d] data2[%d] data3[%d]\n", slope,ch, data2,data3);
          break;
        }
      }
    }
    if (flag==1)
      break;
    else
    {
      fbStart++;
      if (fbStart > (startPt+5) )
        flag=2;
    }  
  }
   
  fbEnd=setup->fdbkNo;
  if (flag==1)
  {
    // second half period
    for (j=fbMiddle+5;j<setup->fdbkNo-4;j++)
    {
      data4=*(chdatamaxbiasPtr+    j*COL_NUM);
      data5=*(chdatamaxbiasPtr+    (j+4)*COL_NUM);
      if (slope > 0 )
      {
        if( ( data4 >= data0) &&  (data5 >= data0 ) )
        { 
          fbEnd=j;
          break;
        }
      }
      else //if (slope < 0 )
      {
        if ( (data4 <= data0)  && (data5 <= data0) )
        {
          fbEnd=j;
          break;
        }
      }
    }
    setup->fluxPeriod[ch]=(fbEnd-fbStart)*setup->stepFDBK;
  }
  else
    setup->fluxPeriod[ch]=(fbEnd-fbStart)*setup->stepFDBK;
    
  #ifdef DEBFLUX
  {
    printf("==== Fluxperiod[ch=%d]=%d     fbStart[%d],fbEnd[%d]  flag= %d ==servo\n",
           ch, setup->fluxPeriod[ch],fbStart,fbEnd,flag);
    printf("  slope(%f):data0[%d] middle[%d]:={[%d],[%d]} end[%d]:={[%d],[%d]}\n\n",
	       slope, data0,fbMiddle,data2, data3,fbEnd, data4, data5 );    
  }
  #endif

  if (myInfo->fpLog !=NULL)
  {	              
    fprintf(myInfo->fpLog,"==== Fluxperiod[ch=%d]=%d     fbStart[%d],fbEnd[%d]  flag= %d ==servo\n",
           ch, setup->fluxPeriod[ch],fbStart,fbEnd,flag);
    fprintf(myInfo->fpLog,"  slope(%f):data0[%d] middle[%d]:={[%d],[%d]} end[%d]:={[%d],[%d]}\n\n",
	       slope, data0,fbMiddle,data2, data3,fbEnd, data4, data5 );           
  }
}


/**
 * \fn void sc2daservosq1lckPoints(ARRAYSET *setup, FILE *fp, char *name,
 *  int *paramarray, int flag, StatusType *status)
 *
 * \brief function:
 *  save sq1 servo lock points for each channel
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray int pointer for setup-parameter array
 * \param  flag       int  1:ad # for comment, 0: no #
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  int   i;

  if (*status != STATUS__OK) return;

  if( strcmp(name," ") !=0)
      fprintf(fp,"%s= ",name);

 
  for (i=0;i<COL_NUM;i++)
  {
    fprintf(fp,"%07d ", paramarray[i]*setup->colMask[i]);
  }
  fprintf(fp,"\n");
}


/**
 * \fn void sc2daservolckgainPoints(ARRAYSET *setup, FILE *fp, char *name,
 *  double *paramarray, int flag, StatusType *status)
 *
 * \brief function:
 *  save servo lock points gain for each channel
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray double pointer for setup-parameter array
 * \param  flag       int  1:ad # for comment, 0: no #
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  int   ch,i;

  if (*status != STATUS__OK) return;

  if( strcmp(name," ") !=0)
    fprintf(fp,"# %s\n",name);

  for (i=0;i<4;i++)
  {
    if (flag==1)
      fprintf(fp,"#");
    for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
      fprintf(fp,"%7.4f ", paramarray[ch]);
    if(i<3)
      fprintf(fp,"<\n");
    else
      fprintf(fp,"\n");
   }
}



/**
 * \fn void sc2daservolckPoints1(ARRAYSET *setup, FILE *fp, char *name,
 *  int *paramarray, StatusType *status)
 *
 * \brief function:
 *  save servo lock points for each channel
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray int pointer for setup-parameter array
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2daservolckPoints1
*/
void sc2daservolckPoints1
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
int        *paramarray,
StatusType *status
)
{
  int   ch,i;

  if (*status != STATUS__OK) return;

  for (i=0;i<4;i++)
  {
    fprintf(fp,"#RC%d\n",i+1);
    fprintf(fp,"%s",name);
    for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
    {
      if ( setup->servo==SQ2FBTRACK)
        fprintf(fp,"%7d ", paramarray[ch]);
      else
        fprintf(fp,"%7d ", paramarray[ch]*setup->colMask[ch]);
    }
    fprintf(fp,"\n");
  }  
}


/**
 * \fn void sc2daservolckgainPoints1(ARRAYSET *setup, FILE *fp, char *name,
 *  double *paramarray, StatusType *status)
 *
 * \brief function:
 *  save servo lock points gain for each channel
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray double pointer for setup-parameter array
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2daservolckgainPoints1
*/
void sc2daservolckgainPoints1
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
double     *paramarray,
StatusType *status
)
{
  int   ch,i;

  if (*status != STATUS__OK) return;

  for (i=0;i<4;i++)
  {
    fprintf(fp,"#RC%d\n",i+1);
    fprintf(fp,"%s",name);
    for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
      fprintf(fp,"%7.4f ", paramarray[ch]*setup->colMask[ch]);
    fprintf(fp,"\n");
  }
}



/**
 * \fn void sc2daservopixellckpts4MCE(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *  FILE *fp, int flag, char *name, void *paramarray, StatusType *status)
 *
 * \brief function:
 *  save servo lock points for each pixel
 *
 * \param  myInfo dasInfoStruct_t structure pointer
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  flag      int  0: paramarray is *int, else is *double
 * \param  name       char point for the param name
 * \param  paramarray void pointer for setup-parameter array
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  int   row, col, pixel;
  int   card;
  int   *intArray;
  double *floatArray;

  if (*status != STATUS__OK) return;

  intArray=NULL;
  floatArray=NULL;

  if( flag==0)       
    intArray=(int*)paramarray;
  else
    floatArray=(double*)paramarray;
	
  for(card=1; card<=4; card++)
  {
    for (col=0;col<ONECARDCHS;col++)
    {
      if( strcmp(name,"I_VAL")==0)
         fprintf(fp,"wb rc%d  gaini%d <\n",card,col);
      else if( strcmp(name,"ADC_OFFSET")==0)
         fprintf(fp,"wb rc%d adc_offset%d <\n",card,col);
      else if( strcmp(name,"FLUXPERIOD")==0)
         fprintf(fp,"wb rc%d flx_quanta%d <\n ",card,col);
      else 
         fprintf(fp,"#wb rc%d xxx%d <\n# ",card,col);

      for (row=0;row<ROW_NUM;row++)
      {
         pixel=row*COL_NUM + col+ (card-1)*ONECARDCHS;

         // set the chosenRow's i=0; for sq2fb servo xg 07-11-08
         //setup->slopSelect[0]<0 for no row turned off
         if( strcmp(name,"I_VAL")==0 && row==setup->slopSelect[0])
            intArray[pixel]=0;

         if( flag==0)       
           fprintf(fp,"%d ", intArray[pixel]);
         else
           fprintf(fp,"%f ", floatArray[pixel]);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  // send pval 
  if( strcmp(name,"I_VAL")==0 && setup->slopSelect[0] > 0 )
  {
    for(card=1; card<=4; card++)
    {
      for (col=0;col<ONECARDCHS;col++)
      {
        fprintf(fp,"wb rc%d  gainp%d ",card,col);
        for (row=0;row<ROW_NUM;row++)
        {
          if ( row==setup->slopSelect[0])
            fprintf(fp,"%d ", setup->slopSelect[1]);
          else 
            fprintf(fp,"0 ");
        }
        fprintf(fp,"\n");
      }
      fprintf(fp,"\n");
    }
  }
}


/**
 * \fn void sc2daservopixellckPoints(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     FILE *fp, int flag, void *paramarray, StatusType *status)
 *
 * \brief function:
 *  save servo lock points for each channel
 *
 * \param  myInfo  dasInfoStruct_t pointer
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  flag    int 0: int array, else double
 * \param  paramarray void pointer for setup-parameter array
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  int   row, col, pixel;
  int   *intArray;
  double *floatArray;

  if (*status != STATUS__OK) return;

  intArray=NULL;
  floatArray=NULL;

  if( flag==0)       
    intArray=(int*)paramarray;
  else
    floatArray=(double*)paramarray;
	
  for (row=0;row<ROW_NUM;row++)
  {
    for (col=0;col<COL_NUM;col++)
    {
       pixel=row*COL_NUM + col;

       if( flag==0)       
         fprintf(fp,"%8d ", intArray[pixel]);
       else
         fprintf(fp,"%8.2f ", floatArray[pixel]);
    }
    if(row < (ROW_NUM-1) )
      fprintf(fp," <\n");
  }
  fprintf(fp,"\n");
}


/**
 * \fn void sc2dasq1maxrefPoints(ARRAYSET *setup, FILE *fp, char *name,
 *  int *paramarray, StatusType *status)
 *
 * \brief function:
 *  save sq1 max ref points for each channel
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray int pointer for setup-parameter array
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dasq1maxrefPoints
*/
void sc2dasq1maxrefPoints
(
ARRAYSET   *setup,
FILE       *fp,
char        *name,
int        *paramarray,
StatusType *status
)
{
  int   ch;

  if (*status != STATUS__OK) return;
 
  fprintf(fp,"\n%s ",name);
  for (ch=0;ch<COL_NUM; ch++)
  {
    if(setup->maxValue.peakVal[ch][0]==NONMODULATION)
      fprintf(fp,"NON ");
    else
      fprintf(fp,"%d ",paramarray[ch]);
  }
  fprintf(fp,"\n");
}


/**
 * \fn void sc2daservosaveFilted(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     int  *biasData, int *filtedbiasData, int bias, FILE *fp,
 *     StatusType *status)
 *
 * \brief function
 *  save biasData and filtedbiasdata  to a file
 * 
 * \param myInfo     dasInfo structure pointer  
 * \param  setup    ARRAYSET structure pointer  
 * \param  biasData data pointer for each bias level 
 * \param  filtedbiasData filteddata pointer for each bias level 
 * \param  bias     which bias is looked at
 * \param  fp       FILE pointer for saving bias/filter
 * \param  status   StatusType pointer  
 *
 * 
 */
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
)
{
  int        ch, i,feedbk;
  
  if (!StatusOkP(status)) return;

  // save to fp
  for (i=0; i< setup->fdbkNo; i++)
  {
    fprintf(fp,"%d %d ", bias*setup->stepBIAS +setup->minBIAS,
             i*setup->stepFDBK +setup->minFDBK);

    for (ch=0; ch < COL_NUM; ch++)
    {
      feedbk=ch + COL_NUM*i;
      fprintf(fp,"%d %d ",biasData[feedbk],filtedbiasData[feedbk] );
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
}  


/**
 * \fn void sc2dasave4nextStep(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, char *databuf, FILE *fp, StatusType *status)
 *
 * \brief function
 *  save different files for next step 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data (max COL*(ROW-1))
 * \param  fp      FILE pointer for result file
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasave4nextStep
*/
void sc2dasave4nextStep
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
)
{

  if (!StatusOkP(status)) return;

  if ( setup->servo==SSARAMP1 || setup->servo==CABLECAL)
    sc2dasave4cableStep(myInfo,setup,databuf,fp,status);
  else if ( setup->servo==SSARAMP || setup->servo==SSALOCK)
    sc2dasave4ssaStep(myInfo,setup,databuf,fp,status);
  else  if ( setup->servo==SQ2SERVO || setup->servo==SQ2OPEN || 
             setup->servo==SQ2LOCK  )
    sc2dasave4sq2Step(myInfo,setup,databuf,fp,status);
  else  if ( setup->servo==SQ1SERVO || setup->servo==SQ1OPEN || 
             setup->servo==SQ1LOCK )
    sc2dasave4sq1Step(myInfo,setup,databuf,fp,status);
    
}


/**
 * \fn void sc2dasave4cableStep(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, char *databuf,FILE *fp, StatusType *status)
 *
 * \brief function
 *  save result files for next step 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data (max COL*(ROW-1))
 * \param  fp      FILE pointer for result file
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasave4cableStep
*/
void sc2dasave4cableStep
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
)
{
  int      i;
  char     tmpfile1[FILE_LEN];
  int      cableOffset[COL_NUM];
  double   cable[COL_NUM];
  FILE     *fpssa1;

  if (*status != STATUS__OK) return;

  if ( setup->servo==CABLECAL )
  {
    for (i=0;i<COL_NUM;i++)
    {
      cableOffset[i]=(int)setup->maxValue.peakVal[i][0];
      cable[i] =setup->maxValue.peakVal[i][1];
    }
    sprintf (tmpfile1, "%s/%s", getenv ( "SC2SCRATCH" ), CABLE_OFFSETSLOPE );
    if((fpssa1 = fopen(tmpfile1,"w")) == NULL )
      {
	*status = DITS__APP_ERROR;
	ErsRep (0, status, "sc2dasave4cableStep: Error- failed to open file %s", tmpfile1); 
	return;
      }

    fprintf(fpssa1,"\n#define CABLEOFFSET_N  only apply to SSA\n");
    sc2daservolckPoints1(setup,fpssa1,"cableOffset=",cableOffset,status);
 
    fprintf(fpssa1,
      "\n#define CABLESLOPE_N  only apply to cable calibration\n");
    sc2daservolckgainPoints1(setup,fpssa1,"cableSlope=",cable,status);
    fclose(fpssa1);
  }
  else // setup->servo==SSARAMP1
  {
    for (i=0;i<COL_NUM;i++)
      cable[i]=setup->maxValue.peakVal[i][0];

    sprintf (tmpfile1, "%s/%s", getenv("SC2SCRATCH"), CABLE_SCALE );       
    if((fpssa1 = fopen(tmpfile1,"w")) == NULL )
      {
	*status = DITS__APP_ERROR;
	ErsRep (0, status, "sc2dasave4cableStep: Error- failed to open file %s", tmpfile1); 
	return;
      }
    fprintf(fpssa1,"\n#define CABLESCALE_N  only apply to SSA\n");
    sc2daservolckgainPoints1(setup,fpssa1,"cableScale=",cable,status);
    fclose(fpssa1);
  }
}


/**
 * \fn void sc2dasave4ssaStep(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, char *databuf, FILE *fp, StatusType *status)
 *
 * \brief function
 *  save different files for next step 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data (max COL*(ROW-1))
 * \param  fp      FILE pointer for result file
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasave4ssaStep
*/
void sc2dasave4ssaStep
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
)
{
  if (*status != STATUS__OK) return;

  if ( setup->servo==SSARAMP )
  {
    // save ssa bias cable offset to ssabiaslock.txt batch file
    sc2dassabiasoffset4batchfile(myInfo,setup,status);

    // append ssabiaslock.txt to mcelock-ssa.txt  
    sc2dassabias4mcelockSSA(myInfo,setup,status);
    // save cable adj Val for cable-correction and ssalock
    sc2dasavecableadjInit(myInfo,setup,status);

    // save ssaout mid Val
    sc2dasavessaoutmidVal(myInfo,setup,status);

    // save ssa bias and cable offset for ssalock
    sc2dasavessabiasOffset(myInfo,setup,status);
  }

  // save ssa bias,initFb,Zfactor and gain from 
  // ssa and ssalock for sq2servo, also ssabiasfdbklock.txt XG 20090507
  sc2dasavessabiasfbZG(myInfo,setup,status);

  if (setup->servo==SSALOCK)
  {
    // sq2lock need updated initFB, gain and Zfact, only after sq2servo, we can 
    // do ssalock. In case, we only run ssalock, we need to update 
    // SQ2LOCK-BIAS-FB-ZG
    // append initFB, Zfactor and Gain from ssalock to SQ2LOCK-BIAS for sq2lock 
    sc2dassafbZGAppd4SQ2lock(myInfo,setup,status);
  }
}



/**
 * \fn void sc2dasave4sq2Step(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, char *databuf,FILE *fp, StatusType *status)
 *
 * \brief function
 *  save different files for next step 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data (max COL*(ROW-1))
 * \param  fp      FILE pointer for result file
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasave4sq2Step
*/
void sc2dasave4sq2Step
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
)
{
  if (*status != STATUS__OK) return;

  if ( setup->servo==SQ2SERVO)
  {
    // save sq2 bias to batch file
    if (setup->stepBIAS !=0)
      sc2dasq2bias4batchfile(myInfo,setup,status);

    //  save bias and fbdk, Zfactor and gain for sq2lock
    sc2dasavesq2biasfbZG(myInfo,setup,status); 

    // save template1 for ssalock to use
    sc2dasq2bias4ssaLock(myInfo,setup,status);
  }

  if ( setup->servo==SQ2SERVO  || setup->servo==SQ2LOCK)
  {
    // save ssa-fb to a batch file,
    sc2dassafb4batchfile(myInfo,setup,status);

    // append mcelock-ssa. txt ssafblock.txt sq2biaslock.txt 
    // to mcelock-sq2.txt
    sc2dassafbsq2bias4mcelockSQ2(myInfo,setup,status);

    // sq2servo and sq2lock save bias fbdk points for sq2open
    sc2dasavesq2biasFB(myInfo,setup,status); 
    sc2dasavesq2Flux(myInfo,setup,status); 
  }
  else // if ( setup->servo==SQ2OPEN )
  {
    sc2dasavesq2openbiasfbZG (myInfo, setup,  status);

    if(setup->servo==SQ2OPEN)
    {
      // this is last stage where Zfact is re-find
      // the Zfact applied to all rows, so, not difference between rows???
      //save separately for MCE firmware, zfact need to be divided by N
      sc2dasq2openZ4SQ1open(myInfo,setup,status);
    }
  }
}


/**
 * \fn void sc2dasave4sq1Step(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, char *databuf, FILE *fp, StatusType *status)
 *
 * \brief function
 *  save different files for next step 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data (max COL*(ROW-1))
 * \param  fp      FILE pointer for result file
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasave4sq1Step
*/
void sc2dasave4sq1Step
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *databuf,
FILE            *fp,
StatusType      *status
)
{
  char   tmpfile[FILE_LEN];
  int    i,fdbk;
  FILE   *fpopt;

  if (*status != STATUS__OK) return;

  if ( setup->servo==SQ1SERVO || setup->servo==SQ1LOCK )
  {
    // =========== original ones  ================================
    sprintf (tmpfile, "%s/%s-org", getenv ( "CURRENTDATADIR" ),SQ1OPTPTS ); 
    // sq1servo script makes copy to  dataFile-SQ1OPTPTS-org
    // open as appending one for all rows
    if((fpopt = fopen(tmpfile,"a")) == NULL )
      {
	*status = DITS__APP_ERROR;
	ErsRep (0, status, "sc2dasave4sq1Step: Error- failed to open file %s", tmpfile); 
	return;
      }
    fprintf(fpopt,"sel_row=%d\n",setup->selRow);
    // Icmax
    sc2daservosq1lckPoints(setup,fpopt,"sq1bias   ",setup->biaslckPt,0,status);
    // sq1fb
    sc2daservosq1lckPoints(setup,fpopt,"sq2fdbk   ",setup->saoutlckVal,0,status);
    // ref sq1bias, sq2fdbk
    sc2daservosq1lckPoints(setup,fpopt,"sq1refbias", setup->biasrefPt,0,status);
    sc2daservosq1lckPoints(setup,fpopt,"sq2reffdbk", setup->refValue.zfact,0,status);
    fclose(fpopt);

    // ============ adjust for each row  =======================
    sprintf (tmpfile, "%s/%s", getenv ( "CURRENTDATADIR" ),SQ1OPTPTS );       
    // sq1servo script makes copy to  dataFile-SQ1OPTPTS
    // open as appending one for all rows
    if((fpopt = fopen(tmpfile,"a")) == NULL )
      {
	*status = DITS__APP_ERROR;
	ErsRep (0, status, "sc2dasave4sq1Step 2: Error- failed to open file %s", tmpfile); 
	return;
      }

    fprintf(fpopt,"sel_row=%d\n",setup->selRow);
    // we need to change if setup->saoutlckVal[i]<=0  or >setup->fluxPeriod[i]
    for (i=0;i<COL_NUM;i++)
    {
      // in case over two periods
      // if setup->biaslckPt=NONMODULATION, then skip, as sq2fb will be 0 too
      // sq2fb flux jump will not affect sq1bias,
      if ( (setup->fluxPeriod[i] !=0) && (setup->biaslckPt[i] !=NONMODULATION ) )
      {
        fdbk=setup->saoutlckVal[i];
        while  (fdbk <=0)  
        {
          fdbk =fdbk+setup->fluxPeriod[i]; 
          setup->refValue.zfact[i] +=setup->fluxPeriod[i];
          setup->saoutlckVal[i]=fdbk;             
        }
        while ( fdbk > setup->fluxPeriod[i])
        {         
          fdbk =fdbk - setup->fluxPeriod[i]; 
          setup->refValue.zfact[i] -=setup->fluxPeriod[i];
          setup->saoutlckVal[i] =fdbk;
        }
      }
    }
    // Icmax
    sc2daservosq1lckPoints(setup,fpopt,"sq1bias   ",setup->biaslckPt,0,status);
    // sq1fb
    sc2daservosq1lckPoints(setup,fpopt,"sq2fdbk   ",setup->saoutlckVal,0,status);
    // ref sq1bias, sq2fdbk
    sc2daservosq1lckPoints(setup,fpopt,"sq1refbias", setup->biasrefPt,0,status);
    sc2daservosq1lckPoints(setup,fpopt,"sq2reffdbk", setup->refValue.zfact,0,status);
    fclose(fpopt);


    sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ), SQ1FBLOCKFILE );       
    // open as appending one for all rows 
    if((fpopt = fopen(tmpfile,"a")) == NULL )
      {
	*status = DITS__APP_ERROR;
	ErsRep (0, status, "sc2dasave4sq1Step 3: Error- failed to open file %s", tmpfile); 
	return;
      }
    fprintf(fpopt,"sel_row=%d\n",setup->selRow);
    // save sq1fb for sq1open to find lock point and gain
    // but may be this is not better than justb use Zfact, as sq2fbopt is no
    // longer the same, so does sq1fb lock points. However, we can use it as 
    // initial start point in rearch for zfact
    sc2daservosq1lckPoints(setup,fpopt,"sq1initfdbk ",setup->initFB,0,status);
    fclose(fpopt);
  }
  else if (setup->servo==SQ1OPEN)// only do it in saveframelckpts 
  {
  }
}


//=======function for re-organise servo setup file =====//
//======================================================//
/**
 * \fn void sc2dassabiasoffset4batchfile(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save ssa bias and cable offset to a batch file, 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dassabiasoffset4batchfile
*/
void sc2dassabiasoffset4batchfile
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile1[FILE_LEN];
  FILE   *fpbias;
  int    i, ch;
  int    *cableadjPtr;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile1, "%s/%s", getenv ( "SC2SCRATCH" ),SSABIASLCK );
  if((fpbias = fopen(tmpfile1,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2download2PCI: Error- failed to open file %s", tmpfile1); 
      return;
    }

  fprintf(fpbias,"# this mce setup file is extracted after ssaramp \n");
  fprintf(fpbias,"# for non-modulation ones, set bias to NONMOD_BIAS (%d)\n",
              NONMOD_BIAS);
  fprintf(fpbias,"# offset =cableOffset[ch]-(int)(biaslckPt[ch]*cableScale[])\n");
  fprintf(fpbias,"#    +  *(cableadjPtr+ch +setup->maxValue.biasStep[ch][0]*COL_NUM) \n");

  // for simple test 
  cableadjPtr=(int *)setup->cableadjPtr;
  for (i=0;i<setup->totalRC2use;i++)
  {
    fprintf (fpbias,"\nwb rc%d sa_bias",setup->whichRC[i] ); 
    for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
    {
       fprintf(fpbias,"%7d ", setup->biaslckPt[ch]);
    }
    fprintf (fpbias,"\nwb rc%d  offset ",setup->whichRC[i] ); 
    for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
    {
      fprintf( fpbias,"%7d ", 
               setup->cableOffset[ch]
               -(int)(setup->biaslckPt[ch]*setup->cableScale[ch])
               + *(cableadjPtr+ch +setup->maxValue.biasStep[ch][0]*COL_NUM) 
             );
    }
  }
  fclose(fpbias);

}


/**
 * \fn void sc2dasavecableadjInit(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save cableadj Init Val for cable correction 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavecableadjInit
*/
void sc2dasavecableadjInit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile1[FILE_LEN];
  FILE   *fp;
  int    i;
  int    *cableadjPtr;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile1, "%s/%s", getenv ( "SC2SCRATCH" ),CABLE_ADJINIT );
  if((fp = fopen(tmpfile1, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavecableadjInit: Error- failed to open file %s", tmpfile1); 
      return;
    }

  fprintf(fp,"# this is extracted after ssaramp at max Modulation\n");
  fprintf(fp,"# cableinitVal= *(cableadjPtr+ch +setup->maxValue.biasStep[ch][0]*COL_NUM) \n");

  // for simple test 
  cableadjPtr=(int *)setup->cableadjPtr;
  for (i=0;i<COL_NUM;i++)
    setup->cableadjInit[i]= *(cableadjPtr+ i +setup->maxValue.biasStep[i][0]*COL_NUM);
  
  sc2daservolckPoints1(setup,fp,"cableadjInit= ",setup->cableadjInit,status); 

  fclose(fp);
}

/**
 * \fn void sc2dassabias4mcelockSSA(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  append ssabiaslock.txt  to mcelock-ssa.txt  
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dassabias4mcelockSSA
*/
void sc2dassabias4mcelockSSA
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   shellCmd[FILE_LEN];

  if (*status != STATUS__OK) return;

  // append ssabiaslock.txt to mcelock-ssa.txt
  // first make a copy from mcelock-template, current it is empty file,
  sprintf(shellCmd,"cp %s/%s %s/%s", getenv ( "CONFIG_ALL" ), 
         MCELOCKTEMPLATE, getenv ( "SC2SCRATCH" ),MCELOCKSSA);     
  system( shellCmd); 

  sprintf(shellCmd,"cat %s/%s >> %s/%s", getenv ( "SC2SCRATCH" ), 
          SSABIASLCK, getenv ( "SC2SCRATCH" ),MCELOCKSSA);     
  system( shellCmd); 
}


/**
 * \fn void sc2dasavessaoutmidVal(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save ssa out mid Val, 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavessaoutmidVal
*/
void sc2dasavessaoutmidVal
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile1[FILE_LEN];
  FILE   *fpmean;
  int    i, ch;
  int    *meanvalPtr;
  int    *cableadjPtr;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile1, "%s-%s", myInfo->dataFile,SSAOUTMIDVAL );       
  if((fpmean = fopen(tmpfile1,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavessaoutmidVal: Error- failed to open file %s", tmpfile1); 
      return;
    }

  fprintf(fpmean,"#ssaout mid val for later analyse, bias midVal.. \n");

  // for simple test 
  cableadjPtr=(int *)setup->cableadjPtr;
  meanvalPtr=(int*)setup->meanvalPtr;

  for( i=0;i<setup->biasNo;i++)
  {
    fprintf(fpmean, "%d_ ", i*setup->stepBIAS+setup->minBIAS);
    for ( ch=0; ch<COL_NUM; ch++ )
      fprintf( fpmean,"%7d ", *(meanvalPtr+ch +i*COL_NUM) );
    fprintf(fpmean,"\n");
  }

  fprintf(fpmean,"\n#cableadj val at each Bias point\n");
  for( i=0;i<setup->biasNo;i++)
  {
    fprintf(fpmean, "%d_ ", i*setup->stepBIAS+setup->minBIAS);
    for ( ch=0; ch<COL_NUM; ch++ )
      fprintf( fpmean,"%7d ", *(cableadjPtr+ch +i*COL_NUM) ); 
    fprintf(fpmean,"\n");
  }
  
  fprintf(fpmean,"\n#cableadj val at max-Bias-point, ch bias adjVal adjScale adjSlope\n");
  for ( ch=0; ch<COL_NUM; ch++ )
  {
    fprintf(fpmean, "%d_  %d_ ", ch, 
      (setup->maxValue.biasStep[ch][0]*setup->stepBIAS+setup->minBIAS) );
    fprintf( fpmean,"%7d ",
          *(cableadjPtr+ch +setup->maxValue.biasStep[ch][0]*COL_NUM) ); 
    fprintf( fpmean,"%f %f \n",setup->cableadjScale[ch],setup->cableSlope[ch]); 
  }
  fclose(fpmean);
}


/**
 * \fn void sc2dasavessabiasOffset(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 * save ssa bias to ssabiasland cable offset for ssalock and others
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavessabiasOffset
*/
void sc2dasavessabiasOffset
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile2[FILE_LEN];
  FILE   *fplck;

  if (*status != STATUS__OK) return;

  // SSALOCK-BIAS-CABLEOFFSET
  sprintf(tmpfile2, "%s/%s", getenv ( "SC2SCRATCH" ),SSALOCK_BIAS_CABLEOFFSET );
  if((fplck = fopen(tmpfile2, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavessabiasOffset: Error- failed to open file %s", tmpfile2); 
      return;
    }
  fprintf(fplck,"\n# for non-modulation ones, set bias to NONMOD_BIAS (%d)\n",
              NONMOD_BIAS);
  fprintf(fplck,"# SSABIASLCK_N  ssabias Val at max modulation\n");
  sc2daservolckPoints1(setup,fplck,"ssabiaslckPt= ",setup->biaslckPt,status);
   
  fprintf(fplck,"\n#define CABLEOFFSET_N  only apply to SSA \n");
  sc2daservolckPoints1(setup,fplck,"cableOffset= ",setup->cableOffset,status);
  fclose(fplck);

  // SSALOCK-BIAS
  sprintf (tmpfile2, "%s/%s", getenv ( "SC2SCRATCH" ),SSALOCK_BIAS );
  if((fplck = fopen(tmpfile2, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavessabiasOffset 2: Error- failed to open file %s", tmpfile2); 
      return;
    }

  fprintf(fplck,"\n# for non-modulation ones, set bias to NONMOD_BIAS (%d)\n",
              NONMOD_BIAS);
  fprintf(fplck,"# BIASLCK_N  ssabias Val at max modulation\n");
  sc2daservolckPoints1(setup,fplck,"biaslckPt= ",setup->biaslckPt,status);
  fclose(fplck);
  
   // SSALOCK-SSABIAS
  sprintf (tmpfile2, "%s/%s", getenv ( "SC2SCRATCH" ),SSALOCK_SSABIAS );
  if((fplck = fopen(tmpfile2, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavessabiasOffset 3: Error- failed to open file %s", tmpfile2); 
      return;
    }
  fprintf(fplck,"\n# for non-modulation ones, set bias to NONMOD_BIAS (%d)\n",
              NONMOD_BIAS);
  fprintf(fplck,"# SSABIASLCK_N  ssabias Val at max modulation\n");
  sc2daservolckPoints1(setup,fplck,"ssabiaslckPt= ",setup->biaslckPt,status);
  fclose(fplck);

}

 
/**
 * \fn void sc2dasavessabiasfbZG(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save biaslck,initFB, Zfactor and Gain from ssa stage for next step 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavessabiasfbZG
*/
void sc2dasavessabiasfbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile3[FILE_LEN];
  FILE   *fpservo;
  char   shellCmd[FILE_LEN+FILE_LEN];

  if (*status != STATUS__OK) return;

  sprintf (tmpfile3, "%s/%s", getenv ( "SC2SCRATCH" ),SSALOCK_BIAS_FB_ZG );
  if((fpservo = fopen(tmpfile3, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavessabiasfbZG: Error- failed to open file %s", tmpfile3); 
      return;
    }
    
  fprintf(fpservo,
   "\n# for non-modulation ones, set bias to NONMOD_BIAS*setup->colMask (%d)\n",
          NONMOD_BIAS);
  fprintf(fpservo,"# BIASLCK_N  ssabias Val at max modulation\n");
  sc2daservolckPoints1(setup,fpservo,"biaslckPt= ",setup->biaslckPt,status);

  fprintf(fpservo,"\n# INITFB_N   initial ssafb lock Value \n");
  fprintf(fpservo,"#  non-modulation ones, set intFB*setup->colMask to (%d)\n",
             NONMODULATION);
  sc2daservolckPoints1(setup,fpservo,"intFB= ",setup->initFB,status);
       
  fprintf(fpservo,"\n# ZFACTOR_N: ssaDAOut Value \n");
  sc2daservolckPoints1(setup,fpservo,"zFactor= ",setup->saoutlckVal,status);
    
  fprintf(fpservo,"\n# GAIN_N for all channels\n");
  sc2daservolckgainPoints1(setup,fpservo,"Gain= ",setup->gain,status);
  fclose(fpservo);

  // add ssabiasfdbklock.txt for setting ssa bias fdbk lock points
  sprintf(shellCmd,"cat $SC2SCRATCH/%s > $SC2SCRATCH/ssabiasfdbklock.txt",SSABIASLCK);
  if( system( shellCmd) !=0)
  {
     *status=DITS__APP_ERROR;
     return ;
  }  
  sprintf (tmpfile3, "%s/ssabiasfdbklock.txt", getenv ( "SC2SCRATCH" ));
  if((fpservo = fopen(tmpfile3, "a")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavessabiasfbZG 2: Error- failed to open file %s", tmpfile3); 
      return;
    }
  fprintf(fpservo,"\n\n# ssa-fdbk lock point initially from SSARAMP or SSALOCK \n");
  fprintf(fpservo,"wb %s flux_fb <\n" , setup->safbCard);
  sc2daservolckPoints(setup,fpservo," ",setup->initFB,0,status); 
  fclose(fpservo);

  // re write ssafblock.txt as it is needed after cablecorrection
  sprintf (tmpfile3, "%s/%s", getenv ( "SC2SCRATCH" ), SSAFBLCK );
  if((fpservo = fopen(tmpfile3, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavessabiasfbZG 3: Error- failed to open file %s", tmpfile3); 
      return;
    }
  fprintf(fpservo,"#this is ssafb value at lock points from SSARAMP/SSALOCK \n");
  fprintf(fpservo,"wb %s flux_fb <\n" , setup->safbCard);
  sc2daservolckPoints(setup,fpservo," ",setup->initFB,0,status); 
  fclose(fpservo);
}


/**
 * \fn void  sc2dassafbZGAppd4SQ2lock(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  append initFB, Zfactor and Gain from ssalock to
 *  SQ2LOCK-BIAS to SQ2LOCK-BIAS-FB-ZG for sq2lock 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+  sc2dassafbZGAppd4SQ2lock
*/
void  sc2dassafbZGAppd4SQ2lock
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile3[FILE_LEN];
  char   shellCmd[FILE_LEN];

  FILE   *fpservo;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile3, "%s/%s", getenv("SC2SCRATCH"),SQ2LOCK_BIAS_FB_ZG );
  sprintf(shellCmd,"cat %s/%s > %s", getenv ( "SC2SCRATCH" ), SQ2LOCK_BIAS, tmpfile3);     
  //printf("This is the shell command I will make %s\n", shellCmd);
  system( shellCmd); 
  if((fpservo = fopen(tmpfile3, "a")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dassafbZGAppd4SQ2lock: Error- failed to open file %s", tmpfile3); 
      return;
    }
 
  fprintf(fpservo,"\n# INITFB_N   initial ssafb lock Value upadted in ssalock\n");
  fprintf(fpservo,"#  non-modulation ones, set intFB*setup->colMask to (%d)\n",
             NONMODULATION);
  sc2daservolckPoints1(setup,fpservo,"intFB= ",setup->initFB,status);
       
  fprintf(fpservo,"\n# ZFACTOR_N: ssaDAOut Value upadted in ssalock\n");
  sc2daservolckPoints1(setup,fpservo,"zFactor= ",setup->saoutlckVal,status);
    
  fprintf(fpservo,"\n# GAIN_N for all channels upadted in ssalock\n");
  sc2daservolckgainPoints1(setup,fpservo,"Gain= ",setup->gain,status);
  fclose(fpservo);
}


/**
 * \fn void sc2dasq2bias4batchfile(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save sq2 bias  to a batch file, 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasq2bias4batchfile
*/
void sc2dasq2bias4batchfile
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile1[FILE_LEN];
  FILE   *fpbias;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile1, "%s/%s", getenv ( "SC2SCRATCH" ), SQ2BIASLCK );  
  if((fpbias = fopen(tmpfile1, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasq2bias4batchfile: Error- failed to open file %s", tmpfile1); 
      return;
    }
  fprintf(fpbias,"# this mce setup file is extracted after sq2servo \n");
  fprintf(fpbias,"# for non-modulation ones, set sq2 bias to NONMOD_BIAS (%d)\n",
            NONMOD_BIAS);
  fprintf(fpbias,"# BIASLCK_N  sq2bias Val at optimal point\n");
  // SQ2 bias 
  fprintf (fpbias,"wb %s flux_fb <\n",setup->sq2biasCard );
  sc2daservolckPoints(setup,fpbias," ",setup->biaslckPt,0,status); 
  fclose(fpbias);
}


/**
 * \fn void sc2dasavesq2biasfbZG(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save biaslck, initFB, Zfactor and Gain from sq2  for next step 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavesq2biasfbZG
*/
void sc2dasavesq2biasfbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile1[FILE_LEN];
  FILE   *fpservo;

  if (*status != STATUS__OK) return;
 
  sprintf (tmpfile1, "%s/%s", getenv ("SC2SCRATCH"), SQ2LOCK_BIAS_FB_ZG );
  if((fpservo = fopen(tmpfile1, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavesq2biasfbZG: Error- failed to open file %s", tmpfile1); 
      return;
    }
  fprintf(fpservo,
      "\n# for non-modulation ones, set bias to NONMOD_BIAS*setup->colMask (%d)\n",
           NONMOD_BIAS);
  fprintf(fpservo,"# BIASLCK_N  sq2bias Val at max modulation, not used\n");
  sc2daservolckPoints1(setup,fpservo,"biaslckPt= ",setup->biaslckPt,status);

  fprintf(fpservo,"\n# INITFB_N   initial ssafb lock Value updated in sq2servo\n");
  fprintf(fpservo,"#  non-modulation ones, set intFB*setup->colMask to (%d)\n",
           NONMODULATION);
  sc2daservolckPoints1(setup,fpservo,"intFB= ",setup->saoutlckVal,status);

  // zfact and gain only updated in sq2servo(as sq2ssalock) in findmaxModul
  //not sq2lock  
  fprintf(fpservo,"\n# ZFACTOR_N: ssaDAOut Value updated in sq2servo \n");
  sc2daservolckPoints1(setup,fpservo,"zFactor= ",setup->maxValue.zfact,status);
    
  fprintf(fpservo,"\n# GAIN_N for all channels updated in sq2servo\n");
  sc2daservolckgainPoints1(setup,fpservo,"Gain= ",setup->maxValue.gain,status);

  fclose(fpservo);
}


/**
 * \fn void sc2dasq2bias4ssaLock(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save sq2 biaslck only for ssalock use 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasq2bias4ssaLock
*/
void sc2dasq2bias4ssaLock
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile1[FILE_LEN];
  FILE   *fpservo;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile1, "%s/%s", getenv ( "SC2SCRATCH" ), SQ2LOCK_BIAS );
  if((fpservo = fopen(tmpfile1, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasq2bias4ssaLock: Error- failed to open file %s", tmpfile1); 
      return;
    }
    
  fprintf(fpservo,
    "\n# for non-modulation ones, set bias to NONMOD_BIAS*setup->colMask (%d)\n",
         NONMOD_BIAS);
  fprintf(fpservo,"# BIASLCK_N  sq2bias Val at max modulation, not used\n");
  sc2daservolckPoints1(setup,fpservo,"biaslckPt= ",setup->biaslckPt,status);
  fclose(fpservo);
}


/**
 * \fn void sc2dassafb4batchfile(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save ssa febdk  to a batch file, 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dassafb4batchfile
*/
void sc2dassafb4batchfile
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile[FILE_LEN];
  FILE   *fpservo;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ), SSAFBLCK );
  if((fpservo = fopen(tmpfile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dassafb4batchfile: Error- failed to open file %s", tmpfile); 
      return;
    }

  // for sq2servo, the curve is ssafb/sq2fb, so, here saoutlckVal is ssafb
  fprintf(fpservo,"#this is ssafb value at lock points from SQ2SERVO/SQ2LOCK used in SQ2OPEN\n");
  fprintf(fpservo,"wb %s flux_fb <\n" , setup->safbCard);
  sc2daservolckPoints(setup,fpservo," ",setup->saoutlckVal,0,status); 
  fclose(fpservo);

  sprintf (tmpfile, "%s/sq2fbinitlock.txt", getenv ( "SC2SCRATCH" ) );
  if((fpservo = fopen(tmpfile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dassafb4batchfile 2: Error- failed to open file %s", tmpfile); 
      return;
    }
  fprintf(fpservo, "#this is sq2servo/sq2lock sq2fb init val \n");
  fprintf(fpservo,"wb bc2 flux_fb < \n");
  sc2dalibsetup_savesq2fboptPts(fpservo,"intFB= ",setup->initFB,1,status);
  fprintf(fpservo,"\n");
  fclose(fpservo);

}


/**
 * \fn void sc2dassafbsq2bias4mcelockSQ2(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  append mcelock-ssa.txt sq2biaslock.txt ssafblock to mcelock-sq2.txt  
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dassafbsq2bias4mcelockSQ2
*/
void sc2dassafbsq2bias4mcelockSQ2
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   shellCmd[FILE_LEN];

  if (*status != STATUS__OK) return;

  // make tmp before appending
  sprintf(shellCmd,"cat %s/%s > %s/tmp", getenv ( "SC2SCRATCH" ), 
           MCELOCKSSA, getenv ( "SC2SCRATCH" ) );     
  //printf("shellcmd=%s\n", shellCmd);
  system( shellCmd); 

  // append sq2biaslock.txt to mcelock-ssa.txt
  sprintf(shellCmd,"cat %s/%s >> %s/tmp", getenv ( "SC2SCRATCH" ), 
          SQ2BIASLCK, getenv ( "SC2SCRATCH" ) );   
  //printf("shellcmd=%s\n", shellCmd);
  system( shellCmd); 

  // append ssafblock.txt to mcelock-ssa.txt 
  sprintf(shellCmd,"cat %s/%s >> %s/tmp", getenv ( "SC2SCRATCH" ), 
          SSAFBLCK, getenv ( "SC2SCRATCH" ) );     
  //printf("shellcmd=%s\n", shellCmd);
  system( shellCmd); 

  //  make a copy to mcelock-sq2,
  sprintf(shellCmd,"cp %s/tmp %s/%s", getenv ( "SC2SCRATCH" ), 
          getenv ( "SC2SCRATCH" ),MCELOCKSQ2);     
  //printf("shellcmd=%s\n", shellCmd);
  system( shellCmd); 
}


/**
 * \fn void sc2dasavesq2biasFB(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save biaslck, initFB, from sq2 and sq2lock for next step 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavesq2biasFB
*/
void sc2dasavesq2biasFB
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile3[FILE_LEN];
  FILE   *fpservo;

  if (*status != STATUS__OK) return;
 
  sprintf (tmpfile3, "%s/%s", getenv ( "SC2SCRATCH" ), SQ2LOCK_BIAS_FB );
  if((fpservo = fopen(tmpfile3, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavesq2biasFB: Error- failed to open file %s", tmpfile3); 
      return;
    }

  fprintf(fpservo,"\n# INITFB_N this is sq2fb value obtained from sq2servo/sq2lock\n");
  fprintf(fpservo,"# it is needed for finding the open look gain\n");
  sc2daservolckPoints1(setup,fpservo,"intFB= ",setup->initFB,status);

  fprintf(fpservo,"\n# for non-modulation ones, set bias to NONMOD_BIAS (%d)\n",
            NONMOD_BIAS);
  fprintf(fpservo,"# BIASLCK_N  sq2bias Val at optimal point, need for\n");
  fprintf(fpservo,"# writing result file in sq2open\n");
  sc2daservolckPoints1(setup,fpservo,"biaslckPt= ",setup->biaslckPt,status);
  fclose(fpservo);   
}


/**
 * \fn void sc2dasavesq2Flux(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save fluxperiod from sq2 and sq2lock for next step 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavesq2Flux
*/
void sc2dasavesq2Flux
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile3[FILE_LEN];
  FILE   *fpservo;

  if (*status != STATUS__OK) return;
 
  // save FLUXPERIOD from sq2  
  sprintf (tmpfile3, "%s/%s", getenv ( "SC2SCRATCH" ),SQ2LOCK_FLUXPERD );
  if((fpservo = fopen(tmpfile3, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavesq2Flux: Error- failed to open file %s", tmpfile3); 
      return;
    }

  fprintf(fpservo,"\n# FLUXPERIOD from sq2 only for sq1 stage \n");
  sc2daservolckPoints1(setup,fpservo,"fluxPeriod= ",setup->fluxPeriod,status);  
  fclose(fpservo);
}


/**
 * \fn void sc2dasavesq2open4pfbZG(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save initFB, Zfactor and Gain fluxperiod from sq2open for sq2fb servo
 * 
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavesq2open4pfbZG
*/
void sc2dasavesq2open4pfbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile[FILE_LEN];
  FILE   *fpservo;

  if (*status != STATUS__OK) return;

  if ( setup->servo ==SQ2OPEN)
      sprintf (tmpfile, "%s/%s-0", getenv ( "SC2SCRATCH" ),SQ2OPEN4P_FB_ZG );
  else  // setup->servo ==SQ2OPEN4P
    sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ),SQ2OPEN4P_FB_ZG );

  if((fpservo = fopen(tmpfile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavesq2open4pfbZG: Error- failed to open file %s", tmpfile); 
      return;
    }

  fprintf(fpservo,"\n# INITFB_N  initial sq2fb lock Value \n");
  fprintf(fpservo,"#  non-modulation ones, set intFB to (%d)\n",NONMODULATION);
  sc2daservolckPoints1(setup,fpservo,"intFB= ",setup->initFB,status);

  fprintf(fpservo,"\n# ZFACTOR_N: ssaDAOut Value \n");
  sc2daservolckPoints1(setup,fpservo,"zFactor= ",setup->saoutlckVal,status);

  fprintf(fpservo,"\n# GAIN_N for all channels\n");
  sc2daservolckgainPoints1(setup,fpservo,"Gain= ",setup->gain,status);

  fprintf(fpservo,"\n# FLUXPERIOD from sq2 only for sq2fb servo \n");
  sc2daservolckPoints1(setup,fpservo,"fluxPeriod= ",setup->fluxPeriod,status);  
  fclose(fpservo);
  
  // save sq2fbinit to sq2open4pfbint.txt
  sprintf (tmpfile, "%s/sq2fbinitlock.txt", getenv ( "SC2SCRATCH" ) );
  if((fpservo = fopen(tmpfile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavesq2open4pfbZG 2: Error- failed to open file %s", tmpfile); 
      return;
    }
  fprintf(fpservo, "#this is sq2open4p sq2fb init val \n");
  fprintf(fpservo,"wb bc2 flux_fb < \n");
  sc2dalibsetup_savesq2fboptPts(fpservo,"intFB= ",setup->initFB,1,status);
  fprintf(fpservo,"\n");
  fclose(fpservo);
  
}



/**
 * \fn void sc2dasavesq2openbiasfbZG(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save biaslck, initFB, Zfactor and Gain from sq2open for next step 
 * also fluxperiod from sq2open4p
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasavesq2openbiasfbZG
*/
void sc2dasavesq2openbiasfbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile[FILE_LEN];
  FILE   *fpservo;

  if (*status != STATUS__OK) return;

  if ( setup->servo ==SQ2OPEN)
    sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ),SQ2OPEN_BIAS_FB_ZG );

  jitDebug(16,"savesq2openbiasfbGZ: %s\n",tmpfile);
  if((fpservo = fopen(tmpfile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavesq2openbiasfbZG: Error- failed to open file %s", tmpfile); 
      return;
    }

  if ( setup->servo ==SQ2OPEN)
  {
    fprintf(fpservo,"\n# for non-modulation ones, set bias to NONMOD_BIAS (%d)\n",
           NONMOD_BIAS);
    fprintf(fpservo,"# BIASLCK_N  sq2bias Val at lock point\n");
    sc2daservolckPoints1(setup,fpservo,"biaslckPt= ",setup->biaslckPt,status);
  }


  fprintf(fpservo,"\n# INITFB_N   initial sq2fb lock Value \n");
  fprintf(fpservo,"#  non-modulation ones, set intFB to (%d)\n",
           NONMODULATION);
  sc2daservolckPoints1(setup,fpservo,"intFB= ",setup->initFB,status);

  /************* sq2open re-find the ssa out at intFB[ch] point
    // shall we use this zfact, instead of what from sq2servo?
    // ssaout: at sq2fb  SQ2OPEN;  
    //setup->saoutlckVal[ch]=*(chdataperbiasPtr);
  **************/

  fprintf(fpservo,"\n# ZFACTOR_N: ssaDAOut Value \n");
  sc2daservolckPoints1(setup,fpservo,"zFactor= ",setup->saoutlckVal,status);

  fprintf(fpservo,"\n# GAIN_N for all channels\n");
  sc2daservolckgainPoints1(setup,fpservo,"Gain= ",setup->gain,status);

  fclose(fpservo);
  //save sq2fb init to sq2openfbint.txt
  sprintf (tmpfile, "%s/sq2fbinit.txt", getenv ( "SC2SCRATCH" ) );
  if((fpservo = fopen(tmpfile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasavesq2openbiasfbZG 2: Error- failed to open file %s", tmpfile); 
      return;
    }
  fprintf(fpservo, "#this is sq2open sq2fb init val \n");
  fprintf(fpservo,"wb bc2 flux_fb < \n");
  sc2dalibsetup_savesq2fboptPts(fpservo,"intFB= ",setup->initFB,1,status);
  fprintf(fpservo,"\n");
  fclose(fpservo);
}


/**
 * \fn void sc2dasq2openZ4SQ1open(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save  Zfactor from sq2open for sq1open 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dasq2openZ4SQ1open
*/
void sc2dasq2openZ4SQ1open
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile1[FILE_LEN];
  FILE   *fpservo,*fpopen;

  if (*status != STATUS__OK) return;

  // this is last stage where Zfact is re-find
  // the Zfact applied to all rows, so, not difference between rows???
  //save separately for MCE firmware, zfact need to be divided by N
    
  sprintf (tmpfile1, "%s/%s", getenv ( "SC2SCRATCH" ),SSALCKDIVBYN );       
  if((fpservo = fopen(tmpfile1, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasq2openZ4SQ1open: Error- failed to open file %s", tmpfile1); 
      return;
    }
  fprintf(fpservo,"\n# ZFACTOR_N: ssaDAOut/samplNo Value updated in sq2open \n");
  sc2dasavessalckDivbyn(setup,fpservo,"zFactor= ",setup->saoutlckVal,0,status);
  fclose(fpservo);

  // also save zfactor to sq1open
  sprintf (tmpfile1, "%s/%s", getenv ( "SC2SCRATCH" ),SQ2OPEN_Z );
  if((fpopen = fopen(tmpfile1, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dasq2openZ4SQ1open 2: Error- failed to open file %s", tmpfile1); 
      return;
    }

  fprintf(fpopen,"\n# ZFACTOR_N: ssaDAOut Value \n");
  sc2daservolckPoints1(setup,fpopen,"zFactor= ",setup->saoutlckVal,status);
  fclose(fpopen);
}


//================================================================//

/**
 * \fn void sc2daservofindslopGain(dasInfoStruct_t *myInfo,ARRAYSET *setup, 
 *     int ch, int *chdatamaxbiasPtr, StatusType *status)
 *
 * \brief function
 *  find lock points and gain only for SSA
 *         given and return ==>G&R        
 *
 * \param myInfo          dasInfoStruct_t pointer
 * \param  setup          ARRAYSET structure pointer  G&R
 * \param  ch             int  which channel is looked at
 * \param  chdatamaxbiasPtr  int pointer for max modulation ch  
 * \param  status         StatusType pointer  G&R
 *
 */
/*+ sc2daservofindslopGain
*/
void sc2daservofindslopGain
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup,
int        ch, 
int        *chdatamaxbiasPtr,
StatusType *status
)
{
  int            j,l;
  int            fbStart,fbStop,data0,data1;
  int            halfPt,aver1,aver2;
  int            data3, data2;

  if (*status != STATUS__OK) return;

  setup->gain[ch]=0;
  setup->saoutlckVal[ch]=0;
  setup->initFB[ch]=NONMODULATION;

  // from prototype data, SSA has two period, use second fdbk,
  // others has one period, so, use the first one
  // is it still ture for SG?
  if ( setup->servo==SSARAMP  || setup->servo==SSALOCK)
  {
    halfPt=setup->maxValue.saOut[ch][1];
    fbStart=setup->maxValue.fdbk[ch][1];
    fbStop =setup->maxValue.fdbk[ch][2];
  }
  else
  {
    halfPt=setup->maxValue.saOut[ch][0];
    fbStart=setup->maxValue.fdbk[ch][0];
    fbStop =setup->maxValue.fdbk[ch][1];
  }
  if (ch==TEST_CH)
  {
    #ifdef DEBSEREACH
      printf("halfpt=%d fbstart=%d fbstop=%d\n",
           halfPt,fbStart, fbStop);
     #endif
  }
  for (j=fbStart;j<fbStop;j++)
  {
    //start to look for saout~=halfPt
    // each chdatmaxbiasPtr has COL_NUM*setup->fdbkNo
    // each ch's data is at j*COL_NUM for fdbk
    data0=*(chdatamaxbiasPtr+    j*COL_NUM);
    data1=*(chdatamaxbiasPtr+(j+1)*COL_NUM);
 
    if (ch==TEST_CH)
    {
      #ifdef DEBSEREACH
        printf("d0=%d d1=%d\n",data0,data1);
      #endif
    }
    if(  (data0 <= halfPt && data1 > halfPt) ||    // positive
         (data0 > halfPt && data1 <= halfPt)       // negative
       )
    {
      if (ch==TEST_CH)
      {
         #ifdef DEBSEREACH
           printf("found the middle point, calculate gain now\n");
         #endif
      }
      if ( setup->servo==SQ1SERVO || setup->servo==SQ1LOCK )
      {
        setup->maxValue.zfact[ch]=data0;
        setup->maxValue.initFB[ch]=j*setup->stepFDBK+setup->minFDBK;
      }

      // saoutlckVal for ssa is Zfact, for sq2servo is ssafb, 
      // for sq1servo is sq2fb 
      setup->saoutlckVal[ch]=data0;
      setup->initFB[ch]=j*setup->stepFDBK+setup->minFDBK;

      // for sq1servo, check the flux jump, 
      // this will be saved into sq1opt-points
      // we are not going to change here, just take what it is for optimised 
      // math calculation to work
      if ( setup->servo==SSARAMP  || setup->servo==SSALOCK)
      {
        // obtain the gain, SQ1SERVO SQ2SERVO do need gain here
        aver1=aver2=0;
        for (l=0; l<AVER_NO; l++)
        {
          data2=*(chdatamaxbiasPtr+ (l         +j)*COL_NUM );
          data3=*(chdatamaxbiasPtr+ (l-AVER_NO +j)*COL_NUM );
          aver1 +=data2;
          aver2 +=data3;
          if (ch==TEST_CH)
          {
            #ifdef DEBSEREACH
              printf("d2=%d d3=%d\n",data2,data3);
            #endif
          }
        }
        // gain =1/tang
        if( aver1==0 || aver2==0)
          setup->gain[ch]=0;
        else
         setup->gain[ch]=
          (double)(AVER_NO*AVER_NO*setup->stepFDBK)/(aver1-aver2);
        if (ch==TEST_CH)
        {
          #ifdef DEBSEREACH
            printf("gain=%f\n",setup->gain[ch]);
          #endif
        }
      }
      break;
    }
  }
}


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
)                                       
/* Description :
    Simple least-squares fit of a straight line. In general the weights
    are expected to be 1.0 or 0.0, allowing points to be ignored.

   History :
    20Oct2004 : original (bdk)
    13Mar2006 : copied from map.c (agg)
    09Oct2006 : check if ren==0, *grad=0. (xg)
*/

{
   double xm;           /* X mean */
   double ym;           /* Y mean */
   double rp;           /* number of values used (total weight) */
   double rnum;         /* XY factor */
   double rden;         /* X squared factor */
   double ximxm;        /* X differences */
   int i;

   if ( !StatusOkP(status) ) return;

   xm = 0.0;
   ym = 0.0;
   rp = 0.0;

   for ( i=0; i<np; i++ )
   {
      xm += x[i] * wt[i];
      ym += y[i] * wt[i];
      rp += wt[i];
   }
   xm = xm / rp;
   ym = ym / rp;

   rnum = 0.0;
   rden = 0.0;

   for ( i=0; i<np; i++ )
   {
      ximxm = ( x[i] - xm ) * wt[i];
      rnum += ximxm * y[i];
      rden += ximxm * ximxm;
   }
   if(rden > 0  || rden < 0)
     *grad = rnum / rden;
   else
     *grad = 0;
   *cons = ym - (*grad) * xm;
}



/**
 * \fn void sc2dafindpixeltwoPeaks(dasInfoStruct_t *myInfo,  ARRAYSET *setup, 
 *     int pixel,int flag,StatusType *status)
 *
 * \brief function
 *  find two peaks only for:SQ1OPEN 
 *  middle point after first peak ( either max or min) from ssaOut(sq1fb) curve
 *  from pixel data, 
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  pixel    int
 * \param  flag     int   flag=0: max then min, =1; min then max
 * \param  status   StatusType pointer  G&R
 *
 */
/*+ sc2dafindpixeltwoPeaks
*/
void sc2dafindpixeltwoPeaks
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int        pixel,
int         flag,
StatusType *status
)
{
  int     startFB,searchPerd;

  if (*status != STATUS__OK) return;

  // move away from the first few data
  startFB=5;
  searchPerd=setup->fdbkNo/4;

  if( flag==0)
  {
    //look for maxPeak
    setup->pixeLock.peakVal[pixel][0]= setup->allData[startFB];
    sc2dafindpixelmaxPeak(myInfo,setup,pixel,0,startFB,searchPerd,status);
    if( setup->pixeLock.peakInx[pixel][0]==0)
      return;
    startFB=setup->pixeLock.peakInx[pixel][0];
    setup->pixeLock.peakVal[pixel][1]=setup->pixeLock.peakVal[pixel][0];
    sc2dafindpixelminPeak(myInfo,setup,pixel,1,startFB,searchPerd,status);
    if( setup->pixeLock.peakInx[pixel][1]==0)
      return;
  }
  else
  { 
   //look for minPeak
    setup->pixeLock.peakVal[pixel][0]= setup->allData[startFB];
    sc2dafindpixelminPeak(myInfo,setup,pixel,0,startFB,searchPerd,status);
    if( setup->pixeLock.peakInx[pixel][0]==0)
      return;
    startFB=setup->pixeLock.peakInx[pixel][0];
    setup->pixeLock.peakVal[pixel][1]=setup->pixeLock.peakVal[pixel][0];
    sc2dafindpixelmaxPeak(myInfo,setup,pixel,1,startFB,searchPerd,status);
    if( setup->pixeLock.peakInx[pixel][1]==0)
      return;
  }
}


/**
 * \fn void sc2dafindpixelmaxPeak(dasInfoStruct_t *myInfo,  ARRAYSET *setup, 
 *     int pixel, int inx, int startFB, int searchPerd, StatusType *status)
 *
 * \brief function
 *  find max peak only for:SQ1OPEN 
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  pixel    int
 * \param  inx     int   either max or min as [pixel][inx]
 * \param  startFB  int   search start
 * \param  searchPerd int   search period
 * \param  status   StatusType pointer  G&R
 *
 */
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
)
{
  int     i,search=0,end;

  if (*status != STATUS__OK) return;

  end=(setup->fdbkNo-1);
  for (i=startFB; i<end; i++)
  {
    search++;
    if( search >searchPerd )
      break;
    if( setup->allData[i] >= setup->pixeLock.peakVal[pixel][inx] )
    {
      setup->pixeLock.peakVal[pixel][inx]= setup->allData[i];
      setup->pixeLock.peakInx[pixel][inx]=i;
      search=0;
    }
  }
  if (i==end) // not see two peaks, not good
    setup->pixeLock.peakInx[pixel][inx]=0;
}
 

/**
 * \fn void sc2dafindpixelminPeak(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     int pixel, int inx, int startFB, int searchPerd, StatusType *status)
 *
 * \brief function
 *  find min peak only for:SQ1OPEN 
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  pixel    int
 * \param  inx     int   either max or min as [pixel][inx]
 * \param  startFB  int   search start
 * \param  searchPerd int   search period
 * \param  status   StatusType pointer  G&R
 *
 */
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
)
{
  int     i,search=0,end;

  if (*status != STATUS__OK) return;
  
  end=(setup->fdbkNo-1);
  for (i=startFB; i<end; i++)
  {
    search++;
    if( search >searchPerd )
      break;
    if( setup->allData[i] <= setup->pixeLock.peakVal[pixel][inx] )
    {
      setup->pixeLock.peakVal[pixel][inx]= setup->allData[i];
      setup->pixeLock.peakInx[pixel][inx]=i;
      search=0;
    }
  }
  if (i==end) // not see two peaks, not good
    setup->pixeLock.peakInx[pixel][inx]=0;
}


// ================ new added  20090902  X. Gao  =============================//
// ==== re-define some function for using filter, not tested with SQs
// ==========================================================================//
/**
 * \fn void sc2daliblckpts_applycolFilter(double *impulse, int  filterOrder,
 *  int  ydim, int *chData, int *filtedData, double *convolPtr, int *status)
 *
 * \brief function
 *  apply filter to single channel Data and save result to filtedData 
 * 
 * \param impluse     double pointer filter impulse coefficient  
 * \param filterOrder int
 * \param ydim        int dimension of data 
 * \param chData      int pointer for single channel data 
 * \param filtedData  int pointer for filteddata  
 * \param convolPtr   double pointer for convolution data 
 * \param status      int pointer  
 *
 * note: (N-1)/2 delay
 * 
 */
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
)
{
  int        i,inx,M, end;
  int        filtedStart, filtedMid,filtedEnd;
  double     xData,conData;

  if (!StatusOkP(status)) return;
 
  M=(filterOrder-1)/2;
  // initial all val =0
  memset( convolPtr, '\0', sizeof(double)*filterOrder );       

  filtedStart=2*M;
  filtedMid  =ydim + M;
  filtedEnd  =ydim + 2*M;
  end        =ydim-1;

  for (i=0; i< filtedEnd; i++)
  {
     if( i < M ) /* 0 M-1 is additional points for convol */
       xData=(double)chData[0]; 
     else if( i < filtedMid )
       xData=(double)chData[(i-M)] ;     
     else   /* ydim + M is additional points for convol*/
       xData=(double)chData[end] ;
   
     sc2daliblckpts_linearConv(impulse,convolPtr, filterOrder, xData, 
                               &conData, status);
  
     // take (N-1)/2 delay off plus the pending ones
     if( i>=filtedStart &&  i<filtedEnd )
     {
       inx=(i-filtedStart);
       filtedData[inx]=(int)conData;
     }
  }
}


/**
 * \fn void sc2daliblckpts_applyFilter(double *impulse, int  filterOrder,
 *  int  ydim, int whichCh, int totalCh, int *orgData, int *filtedData,
 *  double *convolPtr, int *status)
 *
 * \brief function
 *  apply filter to orgData[ydim][totalCh] and save result to filtedData 
 * 
 * \param impluse     double pointer filter impulse coefficient  
 * \param filterOrder int
 * \param ydim         int dimension of data 
 * \param whichCh     int which channel is looked at
 * \param totalCh     int total channel 
 * \param orgData     int pointer for original data 
 * \param filtedData  int pointer for filteddata  
 * \param convolPtr   double pointer for convolution data 
 * \param status      int pointer  
 *
 * note: (N-1)/2 delay
 * 
 */
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
)
{
  int        i,inx,M, end;
  int        filtedStart, filtedMid,filtedEnd;
  double     xData,conData;

  if (!StatusOkP(status)) return;
 
  M=(filterOrder-1)/2;
  // http://www.cppreference.com/stdstring/memset.html
  // initial all val =0
  memset( convolPtr, '\0', sizeof(double)*filterOrder );       


  filtedStart=2*M;
  filtedMid  =ydim + M;
  filtedEnd  =ydim + 2*M;
  end        =ydim -1;

  for (i=0; i< filtedEnd; i++)
  {
     if( i < M ) /* 0 M-1 is additional points for convol */
       xData=(double)orgData[whichCh + totalCh*0]; 
     else if( i < filtedMid )
       xData=(double)orgData[whichCh + totalCh*(i-M)] ;     
     else   /* ydim + M is additional points for convol*/
       xData=(double)orgData[whichCh + totalCh*end] ;
   
     sc2daliblckpts_linearConv(impulse,convolPtr, filterOrder, xData, 
                               &conData, status);
  
     // take (N-1)/2 delay off plus the pending ones
     if( i>=filtedStart &&  i<filtedEnd )
     {
       inx=whichCh+ totalCh*(i-filtedStart);
       filtedData[inx]=(int)conData;
     }
  }
}

/**
 * \fn void sc2daliblckpts_filterFunc(double *impulseRep, int filterOrder, 
 *  int useFlag,int *status)
 *
 * \brief function
 *  get back filter impulse h(n) 
 * 
 * \param impulseRep   double pointer filter impulse coefficient
 * \param filterOrder  int, filter order or tap
 * \param useFlag      int,use coeffienct from thrid party
 * \param  status      int pointer  
 *
 *        = 0: lp windowSinc,SF301, corner=30, order=11
 *        = 1: lp windownedSinc, SF537, order=11, NonewindowType, corner=20
 *        = 2: lp:windownedSinc, SF537  order=31,NonewindowType, corner=20
 *        = 3: bp windownedSinc, SF301, order=31, Stopband=10Hz, passband=20Hz 
 *              passband ripple=1dB, stopband attenuation=20dB
 *
 */

/*+ sc2daliblckpts_filterFunc
*/
void sc2daliblckpts_filterFunc
(
double   *impulse, 
int      filterOrder, 
int      useFlag,
StatusType *status
)
{
  int        i;
 
  if (!StatusOkP(status)) return;

     // scopeFIR windownedSinc, SF301, order=11, NonewindowType, corner=30
    double lpwindsinc1[]=
    {
     0.000566469181072386, 0.040334540202808779, 0.086204541793040809,
     0.128870591221586370, 0.159051937715447480, 0.169943839772088300,
     0.159051937715447480, 0.128870591221586370, 0.086204541793040809,
     0.040334540202808779, 0.000566469181072386
    };

    // scopeFIR windownedSinc, SF537, order=11, NonewindowType, corner=20
    // WB--GO--WB take 1.86ms ~~ 537HZ during servo 
    double lpwindsinc2[]=
    {
     0.078352707236929567, 0.085649640812763120, 0.091586005952151717,
     0.095969996022816334, 0.098659023262730100, 0.099565253425218325,
     0.098659023262730100, 0.095969996022816334, 0.091586005952151717,
     0.085649640812763120, 0.078352707236929567
    };

    // scopeFIR windownedSinc, SF537, order=31, NonewindowType, corner=30
    // WB--GO--WB take 1.86ms ~~ 537HZ during servo 
    double lpwindsinc3[]=
    {
    -0.006602118308907399,-0.002633904069647637, 0.002099427286019177,
     0.007497667287356416, 0.013431540748811949, 0.019746797656457895,
     0.026269413668498057, 0.032811699458903713, 0.039179077192742676,
     0.045177250742100133, 0.050619476344904678, 0.055333633258657818,
     0.059168799972672347, 0.062001060522486648, 0.063738296612492015,
     0.064323763252903010, 0.063738296612492015, 0.062001060522486648,
     0.059168799972672347, 0.055333633258657818, 0.050619476344904678,
     0.045177250742100133, 0.039179077192742676, 0.032811699458903713,
     0.026269413668498057, 0.019746797656457895, 0.013431540748811949,
     0.007497667287356416, 0.002099427286019177,-0.002633904069647637,
    -0.006602118308907399
    };

    // scopeFIR windownedSinc, SF301, order=31, Stopband=10Hz, passband=20Hz 
    // passband ripple=1dB, stopband attenuation=20dB
    double bp[]=
    {
     0.031671143951709085, 0.009867800403131675, 0.009084418761796206,
     0.006445151203763662, 0.001649924649358013,-0.005408521197426909,
    -0.014589662963512530,-0.025574807725622557,-0.037897407775362338,
    -0.050930540610479393,-0.063903781650189229,-0.075983538946147511,
    -0.086368413213330428,-0.094357468934419167,-0.099395015878274540,
     0.898883128394097340,-0.099395015878274540,-0.094357468934419167,
    -0.086368413213330428,-0.075983538946147511,-0.063903781650189229,
    -0.050930540610479393,-0.037897407775362338,-0.025574807725622557,
    -0.014589662963512530,-0.005408521197426909, 0.001649924649358013,
     0.006445151203763662, 0.009084418761796206, 0.009867800403131675,
     0.031671143951709085
    };

    for ( i=0; i< filterOrder; i++)
    {
      if (useFlag==0)
        impulse[i]=lpwindsinc1[i];
      else if (useFlag==1)
        impulse[i]=lpwindsinc2[i];
      else if (useFlag==2)
        impulse[i]=lpwindsinc3[i];
      else
        impulse[i]=bp[i];
    }	
    return;
}


/**
 * \fn void sc2daliblckpts_linearConv(double *impulsePtr, double *convolPtr,
 *  int filterOrder,double xData, double *filtedData, int *status)
 *
 * \brief function
 *  apply linear convolution and save result to filteData 
 * 
 * \param  impulsePtr  double pointer filter impulse coefficient 
 * \param  convolPtr   double pointer for convolution data
 * \param  filterOrder int  
 * \param  xData       double singal data 
 * \param  filtedData  double pointer 
 * \param  status      int pointer  
 *
 * direct convolution to filte data, conArray[] shifting input
 * y(n)=convolution ( k=0, N-1) h(k)x(n-k)
 *
 * note: (N-1)/2 delay
 */
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
)
{
  int    n;
  double temp=0;

  if (*status != STATUS__OK) return;
  
  convolPtr[0]=xData;
  for( n=0; n<filterOrder; n++)
    temp +=impulsePtr[n]*convolPtr[n];
  *filtedData= temp;

  // shift for next data, from high inx to low inx
  for( n=filterOrder-1;  n >0; n--)
    convolPtr[n]=convolPtr[n-1];
}



/**
 * \fn void sc2daliblckpts_finddataSlope(int *data, int start,
 *  double *slope,int *status)
 *
 * \brief function
 *  find the slope of the data  near the start point 
 * 
 * \param data       int pointer for original data 
 * \param start      int the start for search 
 * \param ydim        int the dimension of array 
 * \param slope      double pointer
 * \param status      int pointer  
 *
 */
/*+ sc2daliblckpts_finddataSlope
*/
void sc2daliblckpts_finddataSlope
(
int    *data,
int    start,
double *slope,
StatusType *status
)
{
  int     i;
  int     fitWindow=20;
  double  yData[fitWindow];
  double  xData[fitWindow];
  double  weightData[fitWindow];  
  double  dcLevel;

  if (!StatusOkP(status)) return;

  for (i=0; i<fitWindow; i++)
  {
    xData[i]=i+1;
    weightData[i]=1;
    yData[i]=(double)data[i + start];
    //printf ("%d ",data[i + start] );
  }
  sc2damath_linfit(fitWindow, xData, yData,weightData,slope,&dcLevel,status); 
}




/**
 * \fn void sc2daliblckpts_findFiltered(ARRAYSET *setup,int *chData, 
 *  int *chfiltedData, int *status)
 *
 * \brief function
 *  find filtered data and return in chfiltedData
 *
 * \param  setup          ARRAYSET structure pointer
 * \param  chData         int  channel is looked at
 * \param  chfiltedDatar  int pointer for filtered  
 * \param  status         int pointer  G&R
 *
 */
/*+ sc2daliblckpts_findFiltered
*/
void sc2daliblckpts_findFiltered
(
ARRAYSET   *setup,
int        *chData,
int        *chfiltedData,
StatusType *status
)
{
  int    filterOrder;
  int    filterFlag;

  if (!StatusOkP(status)) return;

  filterFlag=setup->slopSelect[15];
  if ( filterFlag <2) filterOrder=11;
  else                filterOrder=31;
  
  if (setup->modFlag ==0)
  {
    // we only need to call it once, so, use setup->modFlag
    //  max dim of impulsePtr and  convolPtr is 31 , double * 
    sc2daliblckpts_filterFunc(setup->impulsePtr, filterOrder,filterFlag,status);
    setup->modFlag ++;
  }
  sc2daliblckpts_applycolFilter(setup->impulsePtr,filterOrder,setup->fdbkNo,
                                chData,chfiltedData,setup->convolPtr,status);
}


#define DEB_CHK_PEAK_CL   200 // >31 no debug display
/**
 * \fn void sc2daliblckpts_chkPeaks(int *data, int ydim, 
 *  int *peakNo, int servoFlag, StatusType *status)
 *
 * \brief function
 *  find hwo many peaks (min or max)  during [1,ydim]
 *  this is used for SQ2 as some of col has less then a flux period
 *
 * \param  data     int pointer for max modulation ch or pixel data  
 * \param  ydim      int   dimension of the array
 * \param  peakNo    int  pointer
 * \param  status    int pointer  G&R
 *
 */
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
)
{
  int   st, end, noPeak;
  int   i, j, searchPd;
  int   max, maxInx=0, min, minInx=0;

  if (*status != STATUS__OK) return;

  noPeak=0;
  searchPd=ydim/3;
  for ( j=0; j <3; j++ )
  {
    maxInx=minInx=0;

    st=j*searchPd; 
    end =st + searchPd;

    if ( data[st] > data[end] )
    {
       max= data[st]; min=data[end];
    }
    else
    {
       max= data[end]; min=data[st];
    }   
    for ( i=st; i <end; i++ )
    {
       if ( data[i] > max ) { max=data[i]; maxInx=i;}
       if ( data[i] < min ) { min=data[i]; minInx=i;}
    }
    if ( maxInx !=0 || minInx !=0 )
    {
       noPeak ++;
      #ifdef DEB_CHK_PEAK_CL 
       if ( ch == DEB_CHK_PEAK_CL )
         printf("\n%3d_ch = %d-th search,noPeak (%d) maxInx (%d) minInx (%d)\n", 
                ch, j, noPeak, maxInx, minInx);
       #endif
    }
  }
  *peakNo=noPeak;
}


/**
 * \fn void sc2daliblckpts_findsearchPt(int *chData, int ydim, 
 *  int *fbStart, int *meanVal,int *p2pval,int servoFlag, StatusType *status)
 *
 * \brief function
 *  find mean value from two end or minPeak maxPeak , p2p value
 *
 *
 * \param  chData     int pointer for max modulation ch or pixel data  
 * \param  ydim       int   dimension of the array
 * \param  fbStart    int pointer for  fdbk Start
 * \param  meanVal    int  
 * \param  p2pVal     int pointer   p2p value  
 * \param  ch         int  which col  
 * \param  status     StatusType pointer  G&R
 *
 */
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
)
{
  int    i, search0, search1,found;
  int    tryNo, max, min, noPeaks;
  int    mid;

  if (*status != STATUS__OK) return;

  max=min=chData[1]; 
  for  ( i=2; i<ydim; i++)
  {
     if ( chData[i] > max )  max=chData[i] ;
     if ( chData[i] < min )  min=chData[i] ;
  }
  *p2pVal=max-min;

  if ( servoFlag ==SSARAMP || servoFlag ==SSALOCK  || 
       servoFlag ==SQ1OPEN || servoFlag ==SQ1SERVO ||servoFlag ==SQ1LOCK)
  {
    mid=(max+min)/2;
  }
  else
  {
    mid=(chData[1]+ chData[ydim-1])/2;

    // check if we have a less a fluxperiod
    sc2daliblckpts_chkPeaks(chData, ydim, &noPeaks, servoFlag, ch,status);
    #ifdef DEB_CHK_PEAK_CL
      if ( ch == DEB_CHK_PEAK_CL )
         printf("%3d_ch: total Peaks= (%d)\n",ch, noPeaks);
    #endif
    if ( noPeaks ==1)
    {
      *value=chData[5];
      *fbStart=5;
      printf("%3d_ch: this ch only has one peak, search from data[5]\n",ch);
      return;
    }      

  }

  tryNo=0;

  // try to get to the mid value point if the curve is symmetric
  search0=1;
  search1= search0 + 20;
  found=0;

  while (found==0)
  {
    for (i=search0; i<search1; i++)
    {
      if ( (chData[i] >=mid  && chData[i+1] < mid) ||
           (chData[i] < mid && chData[i+1] >= mid) )
      {
        *value=chData[i];
        *fbStart=i;
        found=1;
        break;
      }
    }
    if ( found ==0 )
    {
      tryNo++;
      search0=search1;
      search1= search0 + 20;
    }
  }
}
 


/**
 * \fn void sc2daliblckpts_findfluxPeriod(dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, int ch, int *data, StatusType *status)
 *
 * \brief function
 *  find fluxperiod from SQ2, SQ1OPEN, SQ2OPEN4P
 *  from SG test, it seems that sq2fb has to be 0 to 65535 to have
 *  a full period of ssafb change, 
 *  it only need to do if SQ2SERVO or SQ2LOCK, SQ1OPEN, SQ2OPEN4P
 *
 *
 * \param  myInfo         dasInfoStruct_t pointer
 * \param  setup          ARRAYSET structure pointer  G&R
 * \param  ch             int  which channel is looked at
 * \param  data         int pointer for max modulation ch or pixel data ( org or filtered)  
 * \param  status         StatusType pointer  G&R
 *
 */
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
)
{
//#define TEST_COL  200 // >31 no debug display

  int    j,l,fstEnd,flag=0,startPt=5;
  int    fbStart,fbMiddle,fbEnd;
  int    data0,data2=0, data3=0,data4=0, data5=0;
  int    pixel,startFB, endFB, wdStart;
  double  slope;

  if (*status != STATUS__OK) return;

  fbStart=start;
  fstEnd=setup->fdbkNo-2;
  data0=searchVal;
  pixel=setup->minFlag*COL_NUM + ch;
 

  // start from first data0  
  while (flag==0)
  {
    // sometimes, we may find that fbStart < 10, so, need to make
    // sure finddataSlope has real data
    if (fbStart < 12 )
      wdStart=3;
    else if ( fbStart < (setup->fdbkNo-12) )
      wdStart=fbStart-10;
    else
      wdStart=setup->fdbkNo-12;

    sc2daliblckpts_finddataSlope(data,wdStart,&slope, status); 

    // Now start looking for two points five samples apart that are either
    //  >=data0 or <=data0 depending on the slope we just found
    // (This is the place where the wave passes bak through data0)
    fbMiddle=fbStart;

    #ifdef TEST_COL
    {
      if ( ch == TEST_COL)
        printf ("%3d_ch: fbStart(%d) data0(%d) slope (%3.2f)", ch, fbStart, data0, slope);
    }
    #endif

    for (j=(fbStart+10);j<fstEnd;j++)
    {
      l=j+5;
      data2=data[j];
      if ( (l+1) >= fstEnd )
        l=fstEnd;

      data3=data[l];

      if ( slope > 0 )
      {
        // Positive slope at data0 means we want to find when 
        // the two points are less than data 0
        if( (data2 <= data0)  &&  ( data3 <= data0)) 
        { 
          fbMiddle=j;	  flag=1;
          break;
        }
      }
      else //if (slope < 0 )
      {
        if ( (data2 >= data0) && ( data3 >= data0) )
        {
          fbMiddle=j;	  flag=1;
          break;
        }
      }
    }
    if (flag==1)
      break;
    else
    {
      fbStart++;
      if (fbStart > (startPt+5) )
        flag=2;
    }  
  }

  #ifdef TEST_COL
  {
    if ( ch == TEST_COL)
      printf ("fbMiddle (%d)",fbMiddle);
  }
  #endif
  // Now find the place where the wave yet again passes through data0

  fbEnd=setup->fdbkNo-1;
  if (flag==1)
  {
    // second half period
    for (j=fbMiddle+5;j<setup->fdbkNo-1;j++)
    {
      data4=data[j];
      data5=data[j+1];
      if (slope > 0 )
      {
        if( ( data4 >= data0) &&  (data5 >= data0 ) )
        { 
          fbEnd=j;
          break;
        }
      }
      else //if (slope < 0 )
      {
        if ( (data4 <= data0)  && (data5 <= data0) )
        {
          fbEnd=j;
          break;
        }
      }
    }
  }
  #ifdef TEST_COL
  {
    if ( ch == TEST_COL)
      printf ("fbEnd (%d)\n",fbEnd);
  }
  #endif

  endFB  =  fbEnd*setup->stepFDBK;
  startFB=fbStart*setup->stepFDBK;

  if ( setup->servo !=SQ1OPEN )
    printf("%3d_ch fluxPeriod =%d start(%d) end(%d)\n",
            ch,(endFB-startFB),startFB, endFB);

  setup->fluxPeriod[ch]=(endFB-startFB);
  if ( (setup->servo ==SQ1OPEN && pixel ==TEST_PIXEL) ||
       (setup->servo !=SQ1OPEN )  )
  {

    /* This is really out of place, but I know this code gets 
       executed in setup and only setup so I am going to clear the heater tracking 
       has failed flag here becasue a setup should clear the problem */
    myInfo->heatTrackFailed = 0;

    if (myInfo->fpLog != NULL)
    {	              
      fprintf(myInfo->fpLog,
       "==== Fluxperiod[ch=%d]=%d\tfbStart[%d],fbEnd[%d]  flag= %d =======\n",
           ch, setup->fluxPeriod[ch],startFB,endFB,flag);
      fprintf(myInfo->fpLog,
       "  slope(%f):data0[%d] initStart[%d],middle[%d]:={[%d],[%d]} end[%d]:={[%d],[%d]}\n\n",
	       slope, data0,start, fbMiddle,data2, data3,fbEnd, data4, data5 );           
    }
  }
}



/**
 * \fn void sc2daliblckpts_findframelckptsGain(dasInfoStruct_t *myInfo,  ARRAYSET *setup, 
 *     int * pixelMask, StatusType *status)
 *
 * \brief function
 *  find gain at lock points only for:
 *  SQ1OPEN at ssaout=Zfactor  from ssaOut(sq1fb) curve
 *  from frame data 
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  pixelMask int pointer 
 * \param  status   StatusType pointer  G&R
 *
 */
/*+ sc2daliblckpts_findframelckptsGain
*/
void sc2daliblckpts_findframelckptsGain
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int             *pixelMask,
StatusType      *status
)
{
  char   *frameData=NULL;
  int    *framedataPtr,*framePtr;
  int    headNo=FRAMEHEADER_NUM;
  int    fdbk,frameSize,row,col,pixel;
  BIAS_PEAK  peakInfo;  

  if (!StatusOkP(status)) return;

  if (setup->sampleNo ==0)
  {
    *status=DITS__APP_ERROR;
    sprintf(msg,"sc2daliblckpts_findframelckptGain: sampleNo=0 ");
    ErsRep(0,status,msg);
    return;
  }

  sc2daliblckpts_readBinary (myInfo, &frameData, status);
  if (!StatusOkP(status)) 
  {
    sprintf(msg,
          "sc2daliblckpts_findframelckptGain: _readBinary failed ");
    ErsRep(0,status,msg);
    //  need to  free frameData;
    if ( frameData !=NULL)   free (frameData);
    return;
  }

  // myInfo->bufsize is in byte, including CHKSUM
  frameSize=myInfo->bufSize/4;
  framePtr= (int *)frameData; 

  // set modFlag for _findFiltered only call filterFunc(loading h(n)) once
  setup->modFlag=0;
  for (row=0;row<ROW_NUM;row++)
  { 
    setup->minFlag=row;
    for (col=0;col<COL_NUM;col++)
    {
      pixel=row*COL_NUM+col;
      setup->pixeLock.ssalockPt[pixel]=0;
      setup->pixeLock.sq1FB[pixel]=0;
      setup->pixeLock.gain[pixel]=0;
      setup->pixeLock.iVal[pixel]=0;
      setup->pixeLock.fluxPeriod[pixel]=0;
      setup->pixeLock.peakInx[pixel][0]=0;
      setup->pixeLock.peakInx[pixel][1]=0;
 
      // first check, we use colMask and rowMask, we need sq1fblck from sq1servo
      // it is sq1fblck[ROW_NUM*COL_NUM]
      if(setup->colMask[col] !=0 && setup->rowMask[row]!=0)
      {

       #ifndef FINDLOCKPOINT_PROG
        if ( (setup->pixeLock.sq1initFB[pixel] !=NONMODULATION ) &&
             (setup->sq2fdbkOpt[col] !=NONMODULATION) &&
             (setup->biasoptPt[row] !=NONMODULATION)   )
       #endif
        {
          // collect the pixel data from all frame for each fdbk
          for ( fdbk=0;fdbk<setup->fdbkNo;fdbk++)
          {
            // move to the right entry for each fdbk
            // setup->allData for hold all pixel's data 
            framedataPtr=framePtr + headNo + frameSize*fdbk;
            setup->allData[fdbk]=*(framedataPtr+pixel);
          }
          sc2daliblckpts_findpixellckptsGain(myInfo,setup,row,col,&peakInfo,pixelMask,status);
          if (!StatusOkP(status)) 
          {
            sprintf(msg,
              "sc2daliblckpts_findframelckptGain: _findpixellckptsGain or fluxPeriod failed ");
            ErsRep(0,status,msg);
            //  need to  free frameData;
            if ( frameData !=NULL)   free (frameData);
            return;
          }
        }
      }
    }
  }
  //  need to  free frameData;
  if ( frameData !=NULL)   free (frameData);
}


/**
 * \fn void sc2daliblckpts_findinitlckPts(dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, char *databuf, StatusType *status)
 *
 * \brief function
 *  after found maxModulation, from maxModuation data set
 *  find initial lock points in the middle range in X-axis
 * 
 * \param  myInfo   dasInfoStruct_t pointer
 * \param  setup    ARRAYSET structure pointer  
 * \param  databuf  char pointer for servo data
 * \param  status   StatusType pointer  G&R
 *
 */
/*+ sc2daliblckpts_findinitlckPts
*/
void sc2daliblckpts_findinitlckPts
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup, 
char       *databuf,
StatusType *status
)
{
  int            ch;
  int            *servodataPtr, *chdataperbiasPtr=NULL;
  SERVO_DATAHEAD *rampinfoPtr;
  
if (*status != STATUS__OK) return;
 
  // get different pointers
  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);
 
  // set modFlag for _findFiltered only call filterFunc(loading h(n)) once
  setup->modFlag=0;

  // check for each channel
  if ( setup->slopSelect[15] >=0 )
  {
     fprintf(myInfo->fpLog,"use filtered data (flag=%d)\n",setup->slopSelect[15]);
  }
   
  for (ch=0; ch<COL_NUM; ch++)
  {
    // for  sq1lock, biaslckpt=biasopt[selRow] X.Gao 2009.10.20   
    // apply rowMask[j] to sq1SERVO and sq1LOCK too. X. Gao 2009.10.23
    if (setup->servo==SQ1LOCK ) 
    {
      setup->biaslckPt[ch]=setup->sq1biasOpt[setup->selRow]*setup->rowMask[setup->selRow];
    }

    // point to the max-modulation bias point. but, first check
    if(setup->maxValue.peakVal[ch][0]!=NONMODULATION || setup->servo==SQ1LOCK )
    {
      chdataperbiasPtr=
         servodataPtr +
         setup->maxValue.biasStep[ch][0]*servoInfo.totaldataperBias + 
         ch;

      // get ssaBias (SSARAMP) sq2Bias(SQ2SERVO) or sq1Bias (SQ1SERVO)
      if (setup->servo==SSARAMP || setup->servo ==SQ2SERVO || 
          setup->servo==SQ1SERVO )
      {         
        setup->biaslckPt[ch]=
             setup->maxValue.biasStep[ch][0]*setup->stepBIAS+setup->minBIAS;

        if ( setup->servo==SQ1SERVO )
              setup->biaslckPt[ch] *=setup->rowMask[setup->selRow] ;
      }

      // here we can apply filter to ecah ch if slopSelect[15]=>0
      // the fluxperiod also found either use filtered or original data
      sc2daliblckpts_findslpGain(myInfo,setup,ch, chdataperbiasPtr,status);
    }
    else if (setup->servo !=SSALOCK && setup->servo !=SQ2OPEN )
    {        
       setup->biaslckPt[ch]=NONMOD_BIAS;
    }
  }
}



/**
 * \fn void sc2daliblckpts_findlckptsGain(dasInfoStruct_t *myInfo,
 *   ARRAYSET *setup, char *databuf, FILE *fplck,
     char *servoFile, int *pixelMask, int *status)
 *
 * \brief function
 *  find gain at lock points only for:
 *  SQ2OPEN SQ2OPEN4P at  sq2fb=initFB[i] sq2fdbkOpt[i]from ssaOut(sq2fb) curve 
 *  for slopSelect[0] row, find gain from SQ2OPEN4P binary data for sq2fb servo
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data 
 * \param  fplck   FILE point for the lock result
 * \param  servoFile  prefix filename for servodata.hex 
 * \param  pixelMask int pointer 
 * \param  flag      int  0: for a standalone prog, not use servoFile, 
 *                   don't call saveservoData
 * \param  status  int  given and return
 *
 */
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
)
{
  int            ch;
  int            *servodataPtr, *chdataperbiasPtr=0;
  SERVO_DATAHEAD *rampinfoPtr;

  if (!StatusOkP(status)) return;

  // use SERVO_INFO
  servoInfo.totaldataperBias=COL_NUM*setup->fdbkNo;
  servoInfo.totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
  servoInfo.headSize = sizeof(SERVO_DATAHEAD);  
  servoInfo.totalframeSize = servoInfo.totalptsSize*COL_NUM;
  servoInfo.peakSize = sizeof(BIAS_PEAK)*setup->biasNo;
  servoInfo.maxpeakSize =  sizeof(MAX_P2P);
  servoInfo.totalservoSize = servoInfo.headSize + servoInfo.totalframeSize + 
                             servoInfo.peakSize + servoInfo.maxpeakSize;
 
  // assign different pointers
  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);
 
  // set modFlag for _findFiltered only call filterFunc(loading h(n)) once
  setup->modFlag=0;
  if ( setup->slopSelect[15] >=0 )
  {
     fprintf(myInfo->fpLog,"use filtered data (flag=%d)\n",setup->slopSelect[15]);
  }

  if (setup->servo==SQ2OPEN  )
  {
    for (ch=0; ch<COL_NUM; ch++)
    {
      if(setup->initFB[ch]!=NONMODULATION)
      {
        chdataperbiasPtr=servodataPtr + ch;
        //  pass
        //  setup->saoutlckVal[col]=midVal;
        //  setup->initFB[ch]=midInx*setup->stepFDBK + setup->minFDBK;
        //  setup->gain[ch]=(double)(setup->stepFDBK)/slope;	      
        sc2daliblckpts_findslpGain(myInfo,setup,ch, chdataperbiasPtr,status);
      }
    }
  }

  if ( flag==1)
    sc2dasaveservoData(setup,databuf,servoFile,status);

  sc2dasavelockPts(myInfo,setup,databuf,fplck,status);
  sc2dasave4nextStep(myInfo,setup,databuf,fplck,status); 
 
  // now for sq2fb servo gain, and file sq2open4p-FB-GZ
  // setup->slopSelect[0] > 0 =rowNo turned off for sq2fb servo 
  if ( flag==1)
  {
    jitDebug(16,"sc2daliblckpts_findlckptsGain: slopSelect[0]=%d\n", setup->slopSelect[0]);
    if ( setup->servo==SQ2OPEN4P  )  
      sc2daliblckpts_findsq2fbGain(myInfo,setup,pixelMask,status); 
  }
}


/**
 * \fn void sc2daliblckpts_findlockPts(dasInfoStruct_t *myInfo,
 *   ARRAYSET *setup, char *databuf, FILE *fplck,
 *   char *servoFile, StatusType *status)
 *
 * \brief function
 *  find the max modulation from array of bias-setting data
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  databuf char pointer for servo data 
 * \param  fplck   FILE point for the lock result
 * \param  servoFile  prefix filename for servodata.hex 
 * \param  status StatusType.  given and return
 *
 */
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
)
{
  int        bias, ch;
  int        *servodataPtr, *tmpData;
  char       *peakinfoPtr;
  BIAS_PEAK  *totalpeakPtr,*peakperbiasPtr;
  SERVO_DATAHEAD *rampinfoPtr;

  if (!StatusOkP(status)) return;

  // use SERVO_INFO
  servoInfo.totaldataperBias=COL_NUM*setup->fdbkNo;
  servoInfo.totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
  servoInfo.headSize = sizeof(SERVO_DATAHEAD);  
  servoInfo.totalframeSize = servoInfo.totalptsSize*COL_NUM;
  servoInfo.peakSize = sizeof(BIAS_PEAK)*setup->biasNo;
  servoInfo.maxpeakSize =  sizeof(MAX_P2P);
  servoInfo.totalservoSize = servoInfo.headSize + servoInfo.totalframeSize + 
                             servoInfo.peakSize + servoInfo.maxpeakSize;

  rampinfoPtr  = (SERVO_DATAHEAD *)databuf;
  servodataPtr = (int *)(rampinfoPtr +1);
  peakinfoPtr  = (databuf+ servoInfo.headSize+ servoInfo.totalframeSize);
   
  // point to the entry of array of BIAS_PEAK struture
  totalpeakPtr =(BIAS_PEAK *)peakinfoPtr;

  // for each bias setting, a set of data are taken at different
  // feedbk points
  for (bias=0; bias<setup->biasNo; bias++)
  {
    // point to the first data for each bias
    tmpData=servodataPtr + bias*servoInfo.totaldataperBias;  

    // point to the corresponding  BIAS_PEAK  for each bias
    peakperbiasPtr = totalpeakPtr + bias;
    for( ch=0;ch<COL_NUM; ch++)
    {
      sc2dalookeachCh(setup, tmpData, peakperbiasPtr,ch, status);  
      if(!StatusOkP(status)) 
      {
        ErsRep(0,status,"_findlockPts: failed to call sc2dalookeachCh");
        return ;
      }
    }   
  }
  
  // we have all peak-peak values for each bias setting, 
  // find max. 
  sc2dafindmaxModul(myInfo,setup,databuf,servoFile,status);

  if(setup->servo ==SQ1SERVO || setup->servo ==SQ1LOCK )
    sc2dafindrefmaxModul(myInfo,setup,databuf,status);

  // save binary "fileNameservodata.hex" data so that 
  // standalone findlockpoints can use it.
  if (setup->doServo !=SQ2SERVOSSALOCK)
      sc2dasaveservoData(setup,databuf,servoFile,status);
  
  // max mudoluation may have a flat bottom, we need to avoid it
  if( setup->servo==SSARAMP)
    sc2dafindflatBottom(&servoInfo,setup,databuf,status);

  sc2daliblckpts_findinitlckPts(myInfo,setup,databuf,status);
  if (setup->doServo !=SQ2SERVOSSALOCK)
  {
    // sq1biasMax=setup->biaslckPt[ch]
    // sq2fbMax=setup->maxValue.zfact[ch];
    // sq1fbMax=setup->maxValue.initFB[ch]
    // sq1biasRef=setup->biasrefPt[ch]
    // sq2fbRef=setup->revValue..zfact[ch]; 
    // sq1fbRef=setup->refValue.initFB[ch];
    sc2dasavelockPts(myInfo,setup,databuf,fplck,status);
    sc2dasave4nextStep(myInfo,setup,databuf,fplck,status);
  }
}


/**
 * \fn void sc2daliblckpts_findmidVal(int *data, int *midVal,int *midInx,
 *  int start,int stop,int *status)
 *
 * \brief function
 *  find the middle value from the data bewteen start and stop points 
 * 
 * \param data       int pointer for data
 * \param midVal     int pointer   return  
 * \param midInx     int pointer   return  
 * \param start      int the start for search 
 * \param stop       int the end for search 
 * \param status     int pointer  
 *
 */
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
)
{
  int  i, minInx=0, maxInx=0;
  int  min, max, mid;
  

  if (!StatusOkP(status)) return;

  min=max=data[start];
  //printf("at bejinging: min/max(%d) start(%d) stop(%d)\n", min, start, stop);
  for (i=start; i<stop; i++)
  {
    if( data[i] <= min )
    {
      min=data[i]; minInx=i;
    }
    if( data[i] >= max ) 
    {
      max=data[i]; maxInx=i;
    }
  }
  mid=(max + min)/2;
  *midInx=-1;
  
  // try to get to the mid value point if the curve is symmetric

  for (i=start; i<stop; i++)
  {
    if ( (data[i] >= mid && data[i+1] < mid) ||
         (data[i] < mid && data[i+1] >= mid) )
    {
      *midVal=data[i];
      *midInx=i;
      break;
    }
  }
}


/**
 * \fn void sc2daliblckpts_findmidvalslpFlux(dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, int ch,int *chData, int *status)
 *
 * \brief function
 *  find the middle value and slope around the middle, flux period 
 * 
 * \param myInfo     dasInfoStruct_t pointer
 * \param setup      ARRAYSET structure pointer
 * \param ch         int   
 * \param chData     int pointer for data either org or filtered
 * \param status     int pointer  
 *
 */
  
/*+ sc2daliblckpts_findmidvalslpFlux
*/
void sc2daliblckpts_findmidvalslpFlux
(
dasInfoStruct_t *myInfo,    
ARRAYSET   *setup,
int        ch, 
int        *chData,
StatusType *status
)
{
//#define DEB_CH_LOCKPT  31   // display this ch's lock point
#define SAVE_TEST_CH  1   // save this channel data to $CURRENTDATA/

  int     maxVal[2], minVal[2];
  int     maxInx[2], minInx[2];
  int     fbstartStp[2];
  int     fbStart,fbStop, fdbk, foundFB=0;
  int     midVal, midInx, start=0,data0,p2pVal,peakVal;
  double  slope, thredVal;

  if (!StatusOkP(status)) return;

#ifdef  SAVE_TEST_CH
  int   i; 
  char   tmpFile[FILE_LEN];
  FILE  *fp;

  if ( ch == SAVE_TEST_CH  && setup->servo ==SQ2LOCK)
  {
    sprintf(tmpFile,"%s/saved-sq2lcok-test-ch-%d.txt",getenv ( "CURRENTDATADIR" ),ch); 
    if((fp = fopen(tmpFile, "w")) == NULL)
      {
	*status = DITS__APP_ERROR;
	ErsRep (0, status, "sc2daliblckpts_findmidvalslpFlux: Error- failed to open file %s", tmpFile); 
	return;
      }
    for ( i=0; i<setup->fdbkNo; i++)
      fprintf(fp,"%d\n",chData[i]);
    fclose(fp);
  }
#endif

  setup->gain[ch]=0;
  setup->saoutlckVal[ch]=0;

  sc2daliblckpts_getThreshold(&thredVal, &peakVal,setup->servo,status); 

  sc2daliblckpts_findsearchPt(chData,setup->fdbkNo,&fbStart,&data0,&p2pVal,setup->servo,ch,status);
  if ( abs(p2pVal) <peakVal )
  {
    // If the column is flatlined, then this message does us no good
    if(p2pVal != 0 )
      {
	fprintf (myInfo->fpLog,"_findsearchPt: %3d_ch has p2pValue(%d) < PEAKVAL (%d)\n",ch, p2pVal, peakVal);
      }
    return;
  }
 
  // find the flux period 
  if ( setup->servo ==SQ2SERVO   || setup->servo==SQ2LOCK  ||
       setup->servo ==SQ2OPEN4P || setup->servo ==SQ1OPEN )
  {
    sc2daliblckpts_findfluxPeriod(myInfo,setup,ch,chData,data0,fbStart,status);
  }

  // from the data0, start look for peak,depeanding on the slope at data0
  sc2daliblckpts_findMaxMin(myInfo,chData, setup->fdbkNo, maxVal, maxInx, minVal, minInx, 
                  fbstartStp, fbStart,setup->servo,status);
  if ( minInx[0] <0 ||  maxInx[0] <0 )
  {
    fprintf(myInfo->fpLog,"_findmidvalslpFlux: %3d_ch has either ZERO values or flat or p2p<PEAKVAL\n",ch);
    //MsgOut(status,msg);
    return;
  }
  fbStart=fbstartStp[0];
  fbStop =fbstartStp[1];

  // find mid from [fbStart,fbStop]
  sc2daliblckpts_findmidVal(chData,&midVal,&midInx,fbStart,fbStop,status);
  if ( midInx < 0 )
  {
    sprintf(msg,"liblckpts_findmidvalSlp:%d_ch failed to find mid point\n",ch);  
    // *status=DITS__APP_ERROR;
    ErsRep(0,status,msg);  
    return;
  }     
  if ( setup->servo==SQ1SERVO || setup->servo==SQ1LOCK )
  {    
    setup->maxValue.initFB[ch]=midInx*setup->stepFDBK+setup->minFDBK;
  }
  // saoutlckVal for ssa is Zfact, for sq2servo is ssafb, 
  // for sq1servo is sq2fb 
  setup->maxValue.zfact[ch]=midVal;
  setup->saoutlckVal[ch]=midVal;
  setup->initFB[ch]=midInx*setup->stepFDBK + setup->minFDBK;

  // find the slope / gain
  if ( setup->servo !=SQ2SERVO && setup->servo !=SQ1SERVO  &&
       setup->servo !=SQ2LOCK && setup->servo != SQ1LOCK )
  {
    if ( setup->servo ==SQ2OPEN4P && setup->slopSelect[5] ==1 )    
    {
      // find the lock point at sq2fbOpt
      midInx=(setup->sq2fdbkOpt[ch]-setup->minFDBK)/setup->stepFDBK;
    }
    //if SSARAMP SSALOCK SQ2OPEN SQ1OPEN 
    // (SQ2OPEN4P && setup->slopSelect[5] ==0)    
  
    start=midInx-10; // the fit-window is 20 points

    sc2daliblckpts_finddataSlope(chData,start,&slope,status);
    if( slope!=0 )
    {
      setup->gain[ch]=(double)(setup->stepFDBK)/slope;	   
      setup->maxValue.gain[ch]=setup->gain[ch];

      if (setup->servo ==SQ2OPEN4P && setup->slopSelect[5] ==1 )
      {
        setup->initFB[ch]=setup->sq2fdbkOpt[ch];
        // sq2fbOPT[col] may not = N*stepFDBK
        // find the lock point at sq2fbOpt
        fdbk=(setup->sq2fdbkOpt[ch]-setup->minFDBK)/setup->stepFDBK;
        foundFB=setup->minFDBK + fdbk*setup->stepFDBK;

        setup->saoutlckVal[ch]=chData[fdbk]+ 
                 (int)((double)(foundFB-setup->initFB[ch])/setup->gain[ch]);
      }
    }
  }
  #ifdef  DEB_CH_LOCKPT
  {
    if ( ch == DEB_CH_LOCKPT )
    {
      printf("ch-%d: saoutlckVal(%d) initFB(%d) gain(%f) \n",
        ch, setup->saoutlckVal[ch], setup->initFB[ch],setup->gain[ch]);
    }
  }
  #endif            
}


/**
 * \fn void sc2daliblckpts_findslpGain(dasInfoStruct_t *myInfo,ARRAYSET *setup, 
 *     int ch, int *chdatamaxbiasPtr, int *status)
 *
 * \brief function
 *  find lock points,  gain only for SSA, SQ2OPEN SQ1OPEN, SQ2OPEN4P
 *
 * \param myInfo          dasInfoStruct_t pointer
 * \param  setup          ARRAYSET structure pointer  G&R
 * \param  ch             int  which channel is looked at
 * \param  chdatamaxbiasPtr  int pointer for max modulation ch  
 * \param  status         int pointer  G&R
 *
 */
/*+ sc2daliblckpts_findslpGain
*/
void sc2daliblckpts_findslpGain
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup,
int        ch, 
int        *chdatamaxbiasPtr,
StatusType *status
)
{
#define  DEB_SQ1LOCK_CH 200 // >31 no debug display
  int     j,inx;
  int     *data;
  int     filterFlag;
  int     *chfiltedData=NULL;
  int     *filterservodataPtr;
  SERVO_DATAHEAD *filterrampinfoPtr;

  if (!StatusOkP(status)) return;

  // get the channel data into setup->allData, 
  // if it has not been done, chdatamaxbiasPtr already points to the
  // beginning of this ch
  if ( setup->servo !=SQ2OPEN4P && setup->servo !=SQ1OPEN ) 
  {
    for (j=0; j<setup->fdbkNo; j++)
    {  
      inx=j*COL_NUM ;
      setup->allData[j]=chdatamaxbiasPtr[inx];
    }
  }

  data=setup->allData;
  filterFlag=setup->slopSelect[15];

  // decide if we are going to use filtered data
  if ( filterFlag >= 0 )
  {
    chfiltedData =(int*) calloc(setup->fdbkNo, sizeof(int));
    if ( chfiltedData== NULL)
    {
      sprintf(msg,"sc2daliblckpts_findslpGain:failed to calloc chfiltedData");
      *status =DITS__APP_ERROR;
      ErsRep(0,status,msg);
      return ;
    }
    sc2daliblckpts_findFiltered(setup,setup->allData, chfiltedData,status);

    // store filtered data in setup->filtedPtr, skip the SERVO_DATAHEAD
    filterrampinfoPtr  = (SERVO_DATAHEAD *)setup->filtedPtr;
    filterservodataPtr = (int *)(filterrampinfoPtr +1);

    // for SQ1OPEN, only the setup->selRow filtered data will be saved    
    if (  (setup->servo ==SQ1OPEN && setup->selRow == setup->minFlag) ||
           setup->servo !=SQ1OPEN ) 
    {
      for (j=0; j<setup->fdbkNo; j++)
      {  
        inx= j*COL_NUM + ch;
        filterservodataPtr[inx]=chfiltedData[j];
      }
    }
    // point data to filtered one; 
    data=chfiltedData;
  }

  // find the mid point and slope, pass them to 
  // setup->gain[ch], setup->saoutlckVal[ch], setup->initFB[ch]
  // also flux period
  sc2daliblckpts_findmidvalslpFlux(myInfo,setup,ch,data,status);

  if ( setup->servo ==SQ1LOCK && ch == DEB_SQ1LOCK_CH && setup->selRow == setup->minFlag )
    printf("sq1lock:row_%d:ch_%d sq2fb( %d),initFB( %d)\n",
           setup->selRow, ch, setup->saoutlckVal[ch], setup->initFB[ch]);
  free (chfiltedData);
}



/**
 * \fn void sc2daliblckpts_findsq2fbGain(dasInfoStruct_t *myInfo,  ARRAYSET *setup, 
 *     int *pixelMask, StatusType *status)
 *
 * \brief function
 *  find gain at middle point for SQ1-Bias turned off row from
 *  SQ2OPEN  if slopSelect[0]> 0, use filtered data if slopSelect[15]>=0  
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  pixelMask int pointer 
 * \param  status   StatusType pointer  G&R
 *
 * setup->slopSelect[5]=0: find lock point from middlepoint.  
 *                     =1: find lock point from sq2fbOPT 
 */
/*+ sc2daliblckpts_findsq2fbGain
*/
void sc2daliblckpts_findsq2fbGain
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
int             *pixelMask,
StatusType      *status
)
{
  char   *frameData=NULL, tmpFile[FILE_LEN];
  int    *framedataPtr,*framePtr;
  int    headNo=FRAMEHEADER_NUM;
  int    fdbk,frameSize,row,col,pixel;
  FILE   *fpsq2fb;
 
  if (!StatusOkP(status)) return;

  sc2daliblckpts_readBinary (myInfo, &frameData, status);
  if (!StatusOkP(status)) 
  {
    sprintf(msg,
          "sc2daliblckpts__findsq2fbGain: _readBinary failed ");
    ErsRep(0,status,msg);
    //  need to  free frameData;
    if ( frameData !=NULL)   free (frameData);
    return;
  }

  // myInfo->bufsize is in byte, including CHKSUM
  frameSize=myInfo->bufSize/4;
  row=setup->slopSelect[0];
  framePtr= (int *)frameData;
     
  // save THIS-ROW setup->slopSelect[0] for other package or Execl to plot
  sprintf (tmpFile, "%s-open4p", myInfo->dataFile ); 
  if((fpsq2fb = fopen(tmpFile, "w")) == NULL)
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2daliblckpts_findsq2fbGain: Error- failed to open file %s", tmpFile); 
      return;
    }
  fprintf(fpsq2fb,"# this is column data from SQ2OPEN4P Row(%d) for sq2fb servo\n",setup->slopSelect[0]);
  fprintf(fpsq2fb,"# first col is fdbk \n");
  fprintf(fpsq2fb,"# second col is col_0 \n");
  fprintf(fpsq2fb,"# ................\n");
  fprintf(fpsq2fb,"# last col is col_31 \n\n");

  for ( fdbk=0; fdbk<setup->fdbkNo; fdbk++)
  {
    fprintf(fpsq2fb,"%8d ", setup->minFDBK + fdbk*setup->stepFDBK);
    for (col=0; col<COL_NUM; col++)
    {
      pixel=row*COL_NUM + col;
      framedataPtr=framePtr + headNo + frameSize*fdbk;
      fprintf(fpsq2fb,"%8d ",*(framedataPtr + pixel) );
    } 
    fprintf(fpsq2fb,"\n");
  }
  fclose(fpsq2fb);

  for (col=0; col<COL_NUM; col++)
  {
    pixel=row*COL_NUM + col;
    setup->saoutlckVal[col]=0;
    setup->initFB[col]=0;
    setup->gain[col]=0;
    setup->fluxPeriod[col]=0;

    // collect the pixel data from frame for each fdbk
    // only for setup->slopSelect[0] row
    for ( fdbk=0; fdbk<setup->fdbkNo; fdbk++)
    {
      framedataPtr=framePtr + headNo + frameSize*fdbk;
      setup->allData[fdbk]=*(framedataPtr + pixel);
    }
 
    if(setup->colMask[col] !=0 )
    {
      if(setup->slopSelect[5] ==1)
      {
        if(setup->sq2fdbkOpt[col]!=NONMODULATION)
        {
          // already have setall->allData, pass
         // need to find minInx, maxInx lockpoints at sq2fbopt
          // midInx=(setup->sq2fdbkOpt[col]-setup->minFDBK)/setup->stepFDBK;
          // setup->gain[ch]=(double)(setup->stepFDBK)/slope;	   
          // setup->initFB[col]=setup->sq2fdbkOpt[col];
          //  setup->saoutlckVal[col]=cdData[fdbk]+ 
          //       (int)((double)(foundFB-setup->initFB[col])/setup->gain[col]);
          // also find flux period
          sc2daliblckpts_findslpGain(myInfo,setup,col, NULL,status);
        }
      }
      else
      {
        sc2daliblckpts_findslpGain(myInfo,setup,col, NULL,status);
      }
    }
  }
  jitDebug(16,"sc2daliblckpts_findsq2fbGain: slopSelect[0]=%d\n", setup->slopSelect[0]);
 
  //  need to  free frameData;
  if ( frameData !=NULL)   free (frameData);
  
  sc2dasavesq2open4pfbZG (myInfo,setup, status);
}





/**
 * \fn void sc2daliblckpts_findpixellckptsGain(dasInfoStruct_t *myInfo,  ARRAYSET *setup, 
 *     int row, int col, BIAS_PEAK *peakInfo, int *pixelMask, StatusType *status)
 *
 * \brief function
 *  find gain for:SQ1OPEN  SQ2OPEN4P&&slopSelect[5]=0 around the middle point after first peak 
 *  (either max or min) from ssaOut(sq1fb/sq2fb) curve from pixel data, filter applied if 
 *  slopSelect[15] !=0  peakInfo is not used, only for allowing lookabsMax to work
 *
 * \param myInfo    dasInfo structure pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  row      int
 * \param  col      int
 * \param  peakInfo BIAS_PEAK pointer for the peak, not used
 * \param  pixelMask int pointer
 * \param  status   StatusType pointer  G&R
 *
 */
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
)
{
  int     pixel, chkmaxFlag=0;
  double  ival,ivalScale=1;

  if (!StatusOkP(status)) return;
  
  pixel=row*COL_NUM + col;
  ivalScale=(double)setup->slopSelect[6]/(double)setup->slopSelect[7];

  // check if ABS(max-min)>=PEAK_THD.. allData[] is int*
  // get peakInfo->peakInx[col][0]=minInx, [1]=maxInx
  sc2dalookabsMax(setup->allData,setup,col,&chkmaxFlag,0,peakInfo,status);

  // it is >=PEAK_THD, we start to look the lock point
  if ( chkmaxFlag !=0)
  {
    // already have setup->allData,   
    // setup->saoutlckVal[col]=midVal;
    //  setup->initFB[ch]=midInx*setup->stepFDBK + setup->minFDBK;
    //  setup->gain[ch]=(double)(setup->stepFDBK)/slope;
    // also find flux period
    sc2daliblckpts_findslpGain(myInfo,setup,col, NULL,status);
 
    if ( setup->servo ==SQ1OPEN )  
    {
      setup->pixeLock.ssalockPt[pixel]=
          setup->saoutlckVal[col]/setup->sampleNo;

      setup->pixeLock.gain[pixel]=setup->sampleNo*setup->gain[col];
      // left shift 12 = *4096
      ival= fabs( setup->gain[col] * 4096) ;
      if ( setup->gain[col]  <0)
        setup->pixeLock.iVal[pixel]= -(int)(ival*ivalScale*pixelMask[pixel]);
      else
        setup->pixeLock.iVal[pixel]= (int)(ival*ivalScale*pixelMask[pixel]);

      // midInx  from _findmidvalSlp()
      setup->pixeLock.sq1FB[pixel]=setup->initFB[col];
      setup->pixeLock.fluxPeriod[pixel]=setup->fluxPeriod[col];

    }
  }
}


/**
 * \fn void sc2daliblckpts_readBinary(dasInfoStruct_t *myInfo, 
 *     char **frameData, int *status)
 *
 * \brief function
 *  find lock points,  gain only for SSA, SQ2OPEN SQ1OPEN, SQ2OPEN4P
 *
 * \param myInfo          dasInfoStruct_t pointer
 * \param frameData     char   address of pointer for  frame data  
 * \param  status         int pointer  G&R
 *
 */
/*+ sc2daliblckpts_readBinary
*/
void sc2daliblckpts_readBinary
(
dasInfoStruct_t *myInfo,
char            **frameData,
StatusType      *status
)
{
  char   tmpFile[FILE_LEN];
  int    fileLen;

  if (*status != STATUS__OK) return;

  // open the pixel file, but no appending
  sprintf(tmpFile,"%s-binary",myInfo->dataFile); 
  if((myInfo->fpSq1 = fopen64(tmpFile, "r")) == NULL)
  {
    *status = DITS__APP_ERROR;
    sprintf(msg,"sc2daliblckpts_readBinary: failed to open %s", myInfo->dataFile);
    ErsRep(0, status, msg);
    return;
  }
  sc2dareadbinaryfile2Mem(myInfo,myInfo->fpSq1,frameData, &fileLen,status);
  fclose(myInfo->fpSq1);

  if (!StatusOkP(status)) 
  {
    sprintf(msg,"sc2daliblckpts_findframelckptGain: failed to call sc2dareadbinaryfile2Mem ");
    ErsRep(0,status,msg);
    return;
  }
}

/**
 * \fn void sc2daliblckpts_savefilteredData(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     StatusType *status)
 *
 * \brief function
 *  save filtered data to a file
 * 
 * \param myInfo     dasInfo structure pointer  
 * \param  setup    ARRAYSET structure pointer  
 * \param  status   StatusType pointer  
 * 
 */
/*+ sc2daliblckpts_savefilteredData
*/
void sc2daliblckpts_savefilteredData
(
dasInfoStruct_t *myInfo,    
ARRAYSET        *setup, 
StatusType      *status
)
{
  int        ch,j,inx,fdbk;
  char       file[FILE_LEN];
  FILE       *fp;
  int        *filterservodataPtr;
  SERVO_DATAHEAD *filterrampinfoPtr;
  
  if (!StatusOkP(status)) return;

  // create filtered for later check-up
  sprintf(file,"%s-filtered.txt",myInfo->dataFile);
  if((fp = fopen(file, "w")) == NULL)
  { 
    *status = DITS__APP_ERROR;
    sprintf(msg,"sc2daliblckpts_savefilteredData: failed to open %s", file);
    ErsRep(0, status, msg);
    return;
  }

  // point to setup->filtedPtr, skip the SERVO_DATAHEAD
  filterrampinfoPtr  = (SERVO_DATAHEAD *)setup->filtedPtr;
  filterservodataPtr = (int *)(filterrampinfoPtr +1);

  for (j=0; j<setup->fdbkNo; j++)
  {  
    fdbk=j*setup->stepFDBK + setup->minFDBK,
    fprintf(fp,"%d\t", fdbk); 
    for (ch=0; ch<COL_NUM; ch++)
    {  
      inx= j*COL_NUM + ch;
      fprintf(fp,"%d\t",filterservodataPtr[inx] );
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}  


/**
 * \fn void sc2daliblckpts_findMax(int *data, int ydim,int *maxVal, int *maxInx,
 *  int *startOffset,int ith,int *status)
 *
 * \brief function
 *  find maximum from data and pass the maxVal and maxInx back 
 * 
 * \param data        int pointer for original data 
 * \param ydim         int  array dimension 
 * \param maxVal      int pointer  
 * \param maxInx      int pointer  
 * \param startOffset int pointer the start offset for search,  
 * \param searchPerd  int  search period
 * \param ith         int  ith max
 * \param status      int pointer  
 *
 */
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
)
{
  int i, start, found=0;;
  int search=0;

  if (*status != STATUS__OK) return;
 
  start=*startOffset;
#ifdef DEBUG_MSG
  printf ("%d_findMax: start(%d,%d)\n", ith,start,maxVal[ith]);
#endif
  
  for (i=start; i< (ydim-1); i++)
  {
    search++;
    if( search >searchPerd )
    { 
      break;
    }
    if( data[i] >= maxVal[ith] )
    {
      maxVal[ith]=data[i];  maxInx[ith]=i; search=0;
      *startOffset= maxInx[ith];
      found=1;
    }
  }
  if ( found ==0)
   printf("%d_findMax: not found\n",ith);
}  
 

/**
 * \fn void sc2daliblckpts_findMin(int *data, int ydim,int *minVal, int *minInx,
 *  int *startOffset,int ith,int *status)
 *
 * \brief function
 *  find minimum from data and pass the minVal and minInx back 
 * 
 * \param data       int pointer for original data 
 * \param ydim         int  array dimension 
 * \param minVal      int pointer  
 * \param minInx      int pointer  
 * \param startOffset int pointer the start offset for search, 
 * \param searchPerd  int  search period
 * \param ith         int ith min
 * \param status      int pointer  
 *
 */
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
)
{
  int i,start, found=0;
  int search=0;

  if (*status != STATUS__OK) return;
  
  start=*startOffset;
#ifdef DEBUG_MSG
  printf ("%d_findMin: start(%d,%d)\n", ith,start,minVal[ith]);
#endif
  for (i=start; i< (ydim-1); i++)
  {
    search++;
    if( search >searchPerd )
    { 
      break;
    }
    if( data[i] <= minVal[ith] )
    {
      minVal[ith]=data[i];  minInx[ith]=i; search=0;
      found=1;
      *startOffset= minInx[ith];     
    }
  }
  if ( found ==0)
    printf("%d_findMin: not found\n",ith);
}  


/**
 * \fn void sc2daliblckpts_getThreshold (double *thredVal, int *peakVal,
 *  int servoFlag,int *status)
 *
 * \brief function
 *  get threshold value (slope, peak)
 * 
 * \param thredval    double pointer   slope threshold
 * \param peakval     int pointer   peak threshold
 * \param servoFlag   int 
 * \param status      int pointer  
 *
 */
/*+ sc2daliblckpts_getThreshold
*/
void sc2daliblckpts_getThreshold
(
double *thredVal,
int    *peakVal, 
int    servoFlag,
int    *status
)
{

  if ( *status != STATUS__OK) return;

  if (servoFlag==SSARAMP || servoFlag==SSALOCK)
  {
     *thredVal=DER_THDSSA;
     *peakVal=PEAK_THDSSA;
  }
  else if (servoFlag==SQ2SERVO || servoFlag==SQ2LOCK)
  {
    *thredVal=DER_THDSQ2;
    *peakVal= PEAK_THDSQ2;
  }
  else if (servoFlag==SQ1SERVO || servoFlag ==SQ1LOCK)
  {
    *thredVal=DER_THDSQ1;
    *peakVal= PEAK_THDSQ1;
  }
  else if (servoFlag==SQ1OPEN )
  {
    *thredVal=DER_THDSQ1OPEN;
    *peakVal= PEAK_THDSQ1OPEN;
  }
  else if (servoFlag==SQ2OPEN4P )
  {
    *thredVal=DER_THDSQ2OPEN4P;
    *peakVal= PEAK_THDSQ2OPEN4P;
  }
}


/**
 * \fn void sc2daliblckpts_findMaxMin(dasInfoStruct_t *myInfo,int *data, int ydim,int *maxVal, 
 *  int *maxInx, int *minVal,int *minInx, int *fbstartStp, int startOffset,
 *  int servoFlag,int *status)
 *
 * \brief function
 *  find maximum or minimum value from data and pass the maxVal/minVal 
 *  and maxInd/minInx back 
 * 
 * \param myInfo          dasInfoStruct_t pointer
 * \param data       int pointer for original data 
 * \param ydim         int  array dimension 
 * \param maxVal      int pointer  
 * \param maxInx      int pointer  
 * \param minVal      int pointer  
 * \param minInx      int pointer  
 * \param fbstartStp  int pointer  
 * \param startOffset int the start offset for search, 
 * \param servoFlag   int 
 * \param status      int pointer  
 *
 */
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
)
{
  int    j;
  int    searchPerd;
  int    startPt,wdStart;
  int    peak, peakVal=PEAK_THDSSA;
  double slope, thredVal=DER_THDSSA;

  if (!StatusOkP(status)) return;

  sc2daliblckpts_getThreshold(&thredVal, &peakVal,servoFlag, status); 
 
  startPt=startOffset;
  searchPerd=ydim/4;
  maxVal[0]=maxVal[1]=data[startPt];
  maxInx[0]=maxInx[1]=-1;
  minVal[0]=minVal[1]=data[startPt];
  minInx[0]=minInx[1]=-1;

  // sometimes, we may find that fbStart < 10, so, need to make
  // sure finddataSlope has real data
  if (startPt < 12 )
    wdStart=3;
  else
    wdStart=startPt-10;

  sc2daliblckpts_finddataSlope(data,startPt,&slope, status); 
  if (abs(slope) < thredVal )
  {
    fprintf(myInfo->fpLog,"_findMaxMin: slope (%f) < threshold (%f )\n",slope, thredVal);
    return;
  }
  
  if ( slope > 0 )
  {
    // look for the max first 
    j=0;
    sc2daliblckpts_findMax(data,ydim, maxVal, maxInx, &startPt,searchPerd,
                           j, status);
    minVal[0]=maxVal[0];
    sc2daliblckpts_findMin(data,ydim, minVal, minInx, &startPt,searchPerd,
                           j, status);

    if ( startPt <(ydim-1) )
    {
      // SSA use the second peak
      if (servoFlag==SSARAMP || servoFlag==SSALOCK) 
      {  
        j=1;
        maxVal[1]=minVal[0];
        sc2daliblckpts_findMax(data,ydim, maxVal, maxInx, &startPt,searchPerd,
                               j, status);
      }			       
    }
    if (servoFlag==SSARAMP || servoFlag==SSALOCK)
    {    
      fbstartStp[0]=minInx[0];  fbstartStp[1]=maxInx[1];
    }
    else
    {    
      fbstartStp[0]=maxInx[0];  fbstartStp[1]=minInx[0];
    }
  }
  else
  {
    // look for the min first 
    j=0;
    sc2daliblckpts_findMin(data,ydim, minVal, minInx, &startPt,searchPerd,
                           j, status);
    maxVal[0]=minVal[0];
    sc2daliblckpts_findMax(data,ydim, maxVal, maxInx, &startPt,searchPerd,
                           j, status);
    if ( startPt <(ydim-1) )
    {
      if (servoFlag==SSARAMP || servoFlag==SSALOCK)
      {  
        j=1;
        minVal[1]=maxVal[0];
        sc2daliblckpts_findMin(data,ydim, minVal, minInx, &startPt,searchPerd,
                               j, status);
      }			      
    }
    if (servoFlag==SSARAMP || servoFlag==SSALOCK)
    {    
      fbstartStp[0]=maxInx[0];  fbstartStp[1]=minInx[1];
    }
    else
    {    
      fbstartStp[0]=minInx[0];  fbstartStp[1]=maxInx[0];
    }
  }
  peak=abs(data[fbstartStp[0]] - data[fbstartStp[1]]);  
  if (  peak< peakVal )
  {
     minInx[0]=maxInx[0]=-1;
     fprintf(myInfo->fpLog,"_findMaxMin: peak (%d) < threshold (%d )\n",peak, peakVal);      
  }

  // at here, we at least have maxVal[0], minVal[0], if not found then Inx[]=-1;
  //printf ("min (%d,%d) max(%d,%d)\n", minVal[0],minInx[0],maxVal[0],maxInx[0]);
  //printf ("min (%d,%d) max(%d,%d)\n", minVal[1],minInx[1],maxVal[1],maxInx[1]);
  //printf ("fbstart(%d,%d) fbstop(%d,%d) peak(%d) slope(%0.4f) \n", data[fbstartStp[0]],
  //         fbstartStp[0], data[fbstartStp[1]],fbstartStp[1],peak,slope);

}

/**
 * \fn void sc2daliblckpts_findmaxWanted(FILE, *fp, int *p2pData, 
 *  int whichCh,int ydim, int *maxp2pVal,int *maxp2pInx,
 *  int servoFlag, StatusType *status)
 *
 * \brief function
 *  find maxPeak we like to be from differetnt p2p values of SQ V-Phi  
 *
 *
 * \param  fp         FILE pointer
 * \param  p2pData    int pointer for p2pVal for all bias setting  
 * \param  whichCh    int
 * \param  ydim       int   dimension of the array
 * \param  bias       int pointer  all bias setting   
 * \param  maxp2pVal  int pointer  
 * \param  maxp2pInx  int pointer     
 * \param  thresHold  int threshold for (Max-x/Max)*100      
 * \param  servoFlag  int      
 * \param  status     StatusType pointer  G&R
 *
 */
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
)
{
  int     i;
  int     max, inx=0;
  double  fraction;

  if (!StatusOkP(status)) return;

  max=p2pData[0]; 
  for  ( i=0; i<ydim; i++)
  {
    if ( p2pData[i] > max )
    {
       max=p2pData[i] ;
       inx=i;
    }
  }
  if (max==0)
  {
    fprintf (fp,"%3d_ch: max=0\n",whichCh);
    *maxp2pVal=0;
    return;
  }
  *maxp2pVal=max;
  *maxp2pInx=inx;

  if (thresHold==1) return;

  for  ( i=0; i<ydim; i++)
  {
    fraction=(double)(max-p2pData[i])/max;
    if ( fraction < thresHold )
    {
      *maxp2pVal=p2pData[i];
      *maxp2pInx=i;
      fprintf (fp,"%3d\t%d\t%d\t%d\t%d\t%3.2E\n",
      whichCh,max,bias[inx],*maxp2pVal,bias[i],fraction);
      break;
    }
  }
}


/**
 * \fn void sc2daliblckpts_findmaxModul(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     char *databuf, char *fileName,StatusType *status)
 *
 * \brief function
 *  find maximum modulation from all bias setting data 
 *  allow user to select the wanted MaxModulation as sometime, we see a flat 
 *  p2pVal[] after Icmax and we don't want to set bais too high
 *  slopSelect[13]/slopSelect[14]( ie. 7/1000) defines the thresHold 
 *  if (max-x[i])/max < thresHold, then maxWanted found
 *
 * \param  myInfo   dasInfoStruct_t pointer
 * \param  setup    ARRAYSET structure pointer  G&R
 * \param  databuf  char pointer for servo data  
 * \param  fileName char pointer for first part of servodata.hex  
 * \param  status   StatusType pointer  G&R
 *
 */
/*+ sc2daliblckpts_findmaxModul
*/
void sc2daliblckpts_findmaxModul
(
dasInfoStruct_t *myInfo,
ARRAYSET   *setup, 
char       *databuf,
char       *fileName,
StatusType *status
)
{
  int            j,ch,bias, others;
  double         a,b;
  int            *servodataPtr;
  char           *peakinfoPtr, *maxpeakinfoPtr;  
  BIAS_PEAK      *peakInfo,*peakperbiasPtr;
  SERVO_DATAHEAD *rampinfoPtr;
  MAX_P2P        *maxInfo;
  int            *biasSet, *chp2pVal;
  int            maxWanted, maxInx;
  double         threshold;

  if (!StatusOkP(status)) return;

  // initial all to NONMODULATION
  for (ch=0;ch<COL_NUM;ch++)
  {
    setup->saoutlckVal[ch]=0;
    setup->gain[ch]=(double)0;
    if (setup->servo ==SSARAMP || setup->servo ==SQ2SERVO || setup->servo ==SQ1SERVO )
    {
      setup->biaslckPt[ch]=NONMODULATION;
      setup->initFB[ch]=NONMODULATION;
    }
  }
  // get different pointers
  rampinfoPtr  =  (SERVO_DATAHEAD *)databuf;
  servodataPtr =  (int *)(rampinfoPtr +1);
  peakinfoPtr  =  (databuf + servoInfo.headSize + servoInfo.totalframeSize);
  maxpeakinfoPtr =(peakinfoPtr + servoInfo.peakSize);

  // point to array of biasNo  BIAS_PEAK struture
  peakInfo =(BIAS_PEAK *)peakinfoPtr;

  // point to memory area for holding each ch's peak data
  maxInfo=(MAX_P2P *)maxpeakinfoPtr;

  // allocate for tmp usage 
  chp2pVal=(int*)calloc(setup->biasNo, sizeof(int));
  biasSet=(int*)calloc(setup->biasNo, sizeof(int));

  threshold=(double)setup->slopSelect[13]/setup->slopSelect[14];
  fprintf (myInfo->fpLog,
     "\nch      Max     bias   Wanted   bias   (max-X[i]/max) threshold=%3.2E\n",threshold);

  // check the three p2pValue
  for (ch=0; ch<COL_NUM; ch++)
  {
    // for all these, if biaslckPt[] ==NONMOD_BIAS, skip and set no modulation
    if ( ( setup->servo==SSALOCK ||setup->servo==SQ2OPEN || setup->servo==SQ1OPEN ||
           setup->servo==SQ2LOCK || setup->servo==SQ1LOCK ) &&  
           (setup->biaslckPt[ch] ==NONMOD_BIAS)  
       )	 
      others=1;
    else if (setup->nomodulCount[ch]==setup->biasNo)
      others=1;
    else 
      others=0;

    // found that all three |p2pValue| are similar, so, only check one    
    // always starts from the peakInfo
    // we need to readjust, p2pVal must have > 0, so use 0
    if  (peakInfo->p2pValue[ch][0]==NONMODULATION)
      maxInfo->maxp2p[ch][0]=0;
    else
      maxInfo->maxp2p[ch][0]=peakInfo->p2pValue[ch][0];

    // skip this bias setting, set p2p=0, later set to NONMODULATION
    // SSALOCK if biaslckPt=NOMOD_BIAS, 
    if( others==1)
    {
      maxInfo->maxp2p[ch][0]=0;
      maxInfo->maxinx0[ch][0]=0;
    }

    for (bias=0; bias<setup->biasNo; bias++)
    {
      // point to the corresponding BIAS_PEAK entry for each bias
      // if biasNo=1, there is only one, which is max for SSALOCK, SQ2LOCK SQ1LOCK
      biasSet[bias]=setup->minBIAS + bias*setup->stepBIAS;
      peakperbiasPtr=peakInfo + bias;

      if ( peakperbiasPtr->p2pValue[ch][0] ==NONMODULATION)
        chp2pVal[bias]=0;
      else
        chp2pVal[bias]=fabs(peakperbiasPtr->p2pValue[ch][0]);
    }

    // if (max-x[i])/max < thresHold, then maxWanted found
    sc2daliblckpts_findmaxWanted(myInfo->fpLog,chp2pVal,ch,setup->biasNo,
           biasSet,&maxWanted,&maxInx, threshold,setup->servo,status);

    if (threshold==1)
        fprintf (myInfo->fpLog,"%3d\t%d\t%d\n",ch,maxWanted,biasSet[maxInx]);
     
    if (maxWanted !=0)    
    {
      peakperbiasPtr=peakInfo + maxInx;
      for (j=0;j<3;j++)
      {
        maxInfo->maxp2p[ch][j]=peakperbiasPtr->p2pValue[ch][j];
        maxInfo->maxinx1[ch][j]=peakperbiasPtr->peakInx[ch][j];
        maxInfo->maxinx0[ch][j]=bias;

        //init saOut Value is half way betweeen MAx-MIN or MIN-MAX
        a=peakperbiasPtr->peakPts[ch][j];
        b=peakperbiasPtr->peakPts[ch][j+1];          
        setup->maxValue.saOut[ch][j]=(a+b)/2;
      }
      maxInfo->maxinx1[ch][3]=peakperbiasPtr->peakInx[ch][3]; 
  
      // save here for SQ2SERVOSSALOCK, we are only interested in 
      // zfact and gain, initFB will be picked in initlckPts() 
      //setup->maxValue.initFB[ch]=peakperbiasPtr->initFB[ch];
      setup->maxValue.zfact[ch]=peakperbiasPtr->zfact[ch];
      setup->maxValue.gain[ch]=peakperbiasPtr->gain[ch];
    }
    // have 3 max peaks for each channel, copy them to setup structure   
    for(j=0; j<3; j++)
    {
      if (maxInfo->maxp2p[ch][0]!=0)
      {
        setup->maxValue.peakVal[ch][j]=maxInfo->maxp2p[ch][j];
 
        //ssabiasMax, sq2biasmax, sq1biasMax
        setup->maxValue.biasStep[ch][j]=maxInfo->maxinx0[ch][j];
 
        // ssafbMax, sq2fbMax or sq1fbMax
        setup->maxValue.fdbk[ch][j]=maxInfo->maxinx1[ch][j];
      }
      else
      {
        setup->maxValue.peakVal[ch][j] =NONMODULATION ;
        setup->maxValue.biasStep[ch][j]=NONMODULATION ;
        setup->maxValue.fdbk[ch][j]=NONMODULATION ;
      }
    }
    if (maxInfo->maxp2p[ch][0]!=0 )
    {
      setup->maxValue.fdbk[ch][3]=maxInfo->maxinx1[ch][3];
    }
    else
    {
      setup->maxValue.fdbk[ch][3]=NONMODULATION ;
    }
  }
  free(chp2pVal); free(biasSet);
}





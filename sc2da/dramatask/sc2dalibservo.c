/* 
 * \file sc2dalibservo.c
 *
 * \brief collection of sub-functions for sc2da.c  
 *        servo algorithm
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
 *  $Log: sc2dalibservo.c,v $
 *  Revision 1.30  2011/09/02 19:16:21  bgorges
 *  Standardized entry checks on all fuctions/actions that should quietly return if entry status is bad.
 *
 *  Revision 1.29  2011/06/02 21:31:22  bgorges
 *  Checking all of fopen and fopen64 to make sure a file is actually opened.
 *
 *  Revision 1.28  2011/04/22 00:25:47  cwalther
 *  I put a considerable number of comments in this code and it has the changes made to normalize thresholds when the number of samples changes
 *
 *  Revision 1.27  2010/10/21 02:01:28  cwalther
 *  Removed optimal-sq1bias-sq2fb from CONFIG_HARD it is in scratch now
 *
 *  Revision 1.26  2010/10/07 23:32:15  cwalther
 *  Trying to make the log file more useable
 *
 *  Revision 1.25  2010/08/03 22:16:10  cwalther
 *  Changes to make sc2_setup work with SC2SCRATCH
 *
 *  Revision 1.24  2010/01/26 22:21:49  cwalther
 *  Commented out some regularly occuring messages
 *
 *  Revision 1.1.1.1  2007/05/16 08:27:00  dkelly
 *  first insertion
 *
 *
 */

// some functions defined in here
#include "sc2dalib.h"
#include "sc2dalibsetup.h"

/* global  info */
static FILE  *servoFp, *fluxjmpFp;
 
// use SDSU Timing board
//#define SDSU

#define   CH_TEST   3
#define   SQ2FBSTEPDOWN 19   // for sq2open4p

///////////////////=== servo special =============////////
///////////////////===============================///////

/**
 * \fn void sc2daservoAlgrm(SDSU_CONTEXT *con, char *byte, 
 *  dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *  int *flagPtr,int fdbkInx, char *databuf,StatusType *status)
 *
 * \brief function: 
 *  if( setup->servo ==SQ2SERVO/lock  || setup->servo ==SQ1SERVO lock)
 *      re-adjust feedback value using servo algorithm
 *  else
 *    do not apply servo algorithm to feedback 
 *
 *  doServo<0    save the whole frame, also apply servo Algorithm to Sq2 and Sq1
 * 
 *  for result data file
 *   setup->servo==SSARAMP,or SSALOCK,or SQ2BIASING or SSARAMP1 or CABLECAL
 *                                           save setup->fdbk
 *   setup->servo==SQ1SERVO lock or SQ2SERVO  lock     save setup->initFB[i]
 *   setup->servo==SQ1OPEN   or SQ2OPEN      save data[i]
 *
 *
 * for channel data
 *   save initFB[] if setup->servo==SQ1SERVO lock or SQ2SERVO lock
 *   always save data[] for the rest,  
 *
 * \param  con    SDSU_CONTEXT pointer
 * \param  byte   pointer for current raw frame data buffer
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  flagPtr int pointer for   
 * \param  fdbkInx int: feedback index in ramping
 * \param  databuf  char pointer for servo data 
 *                 ( save only needed CHs for later process)
 * \param  status  StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * the data is always placed in RC1RC3RC4 and zeros are inserted for 
 * missd RCs. see dhtask.  data structure 
 *  RC1=[0: 7]; RC2=[8:15]  RC3=[16:23]; RC4=[24:31]   
 *
 */
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
)
{
  int         *longword;
  int         dataInx,i,j,read8;
  int         data[COL_NUM];  //COL_NUM=32 
  double      gain;
  SERVO_DATAHEAD *rampinfoPtr;

  // keep the pointer to memory holding servo data
  static int  *storePtr1;     
  static int  *storePtr2;     
  static int   showDebug;
   
  if (!StatusOkP(status)) return;

  // initial data[]
  for (i=0;i<COL_NUM;i++)
    data[i]=0;

  // skip the header
  longword = (int *)( byte+HEADERSET);  

  // initialise the storePtr to the beginning of databuf
  if(*flagPtr ==0)  
  {
    // move away from infomation head and get the right entry pointer
    rampinfoPtr=(SERVO_DATAHEAD *)databuf;
    if( setup->doServo !=SQ2SERVOSSALOCK)
      storePtr1=(int *)(rampinfoPtr +1);
    else
      storePtr2=(int *)(rampinfoPtr +1);
    *flagPtr=1;
    showDebug=1;
  }
  
  // row_40 is dark row. do not apply to the missing cards
  dataInx=(setup->selRow)*COL_NUM;
  if (showDebug==1 )
    jitDebug(16,"_servoAlg:dataInx=%d\n",dataInx);

  for(i=0;i<setup->totalRC2use;i++)
  {
    read8=(setup->whichRC[i]-1)*ONECARDCHS;
    for(j=read8; j<(read8+ONECARDCHS); j++ )
    {
      data[j]=*(longword+dataInx+j);
      if ( (j==0 || j==8 || j==16 || j==24 ) && showDebug==1 )
        jitDebug(16,"data[%d]=%d ",j,data[j]);
 
      if( setup->servo ==SQ2SERVO     || setup->servo ==SQ2LOCK ||
          setup->servo ==SQ1SERVO     || setup->servo ==SQ1LOCK ||
          setup->servo ==TESBIASSERVO || setup->servo ==SQ1BIASSERVO 
        )
      {
        // the equation is based on positive slop, 
        // inittFB : ssafb for SQ2SERVO, SQ2LOCK
        // inittFB : sq2fb for SQ1SERVO  sq1lock
	  gain= setup->gain[j];	    

        if ( setup->servo ==SQ1SERVO || setup->servo ==SQ1BIASSERVO ||
	     setup->servo ==SQ1LOCK  || setup->servo ==TESBIASSERVO )
        {     
           gain *=setup->sq1gainScale[j];
        }

        setup->initFB[j] -= gain*(data[j]-setup->zfact[j]);
        setup->initFB[j] *=setup->colMask[j];
      }
    }
  }
  if (showDebug==1 )
    jitDebug(16,"\n");

  showDebug ++;
  showDebug = showDebug%100;
  
  // store channel data to the dataduf for later to find the optimal points.
  // for SQ2 and SQ1, save initFB[] to channel data if 
  // setup->servo ==SQ2SERVO lock || setup->servo ==SQ1SERVO lock
  // for the rest, always save data[] to channel data
  // the missing cards' val kept the same as in setup file
  for (j=0;j<COL_NUM;j++)
  {
    if( setup->servo ==SQ2SERVO  || setup->servo ==SQ2LOCK  ||
        setup->servo ==SQ1SERVO  || setup->servo ==SQ1LOCK  ||  
        setup->servo ==TESBIASSERVO
       )
    {
      *(storePtr1+j)=setup->initFB[j];
    }
    else 
    {
      if( setup->doServo !=SQ2SERVOSSALOCK)
        *(storePtr1+j)=data[j]*setup->colMask[j];     
      else
         *(storePtr2+j)=data[j]*setup->colMask[j];
    }
    // save the ch data now to *allData for finding the meanSSA
    if ( setup->servo ==SSARAMP || setup->servo ==CABLECORET)
        setup->allData[j + COL_NUM*fdbkInx] = data[j]*setup->colMask[j]; 
  }
  // advance storePtr by
  if( setup->doServo !=SQ2SERVOSSALOCK)
    storePtr1 +=COL_NUM;
  else
    storePtr2 +=COL_NUM;

  if(setup->doServo !=SQ2SERVOSSALOCK)
  {
    // for using idl script to do data processing/display
    if (myInfo->filesvFlag!=0)
      sc2daservosaveData(con,myInfo,setup,byte,data,fdbkInx,status);
  }
}


/**
 * \fn void sc2daservodataheadInit(ARRAYSET *setup, char *servoData,
 *  StatusType *status)
 *
 * \brief function 
 *  populate the servo data head with all informaion
 *
 *  SERVO_DATAHEAD:string information for ramping bias, fb, card, intfb .. etc
 *
 * \param setup   ARRAYSET structure pointer
 * \param servoData  pointer for servo data
 * \param status     StatusType     
 *
 */
/*+ sc2daservodataheadInit
*/
void sc2daservodataheadInit
(
ARRAYSET              *setup,
char                  *servoData,
StatusType            *status
)
{
  int    j;
  char   tmp[CMD_LEN];
  char   *rampinfoPtr;

  if (*status != STATUS__OK) return;
 
  // store ramping setting in servoData buf
  rampinfoPtr=(char *)servoData;

  sprintf(rampinfoPtr,"#servo=%s\n", setup->servoName); 
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#whichCard=");
  for ( j=0; j<setup->totalRC2use; j++ )
  {
    sprintf (tmp, " %d", setup->whichRC[j] );
    strcat(rampinfoPtr,tmp); 
  }
  strcat(rampinfoPtr,"\n"); 
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#bias=%d\n", setup->minBIAS); 
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#bstep=%d\n",setup->stepBIAS); 
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#nbias=%d\n",  setup->biasNo);  
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#feed=%d\n",setup->minFDBK); 
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#fstep=%d\n",setup->stepFDBK); 
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#nfeed=%d\n",setup->fdbkNo);   
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#row=%d\n",setup->selRow);   
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#doservo=%d\n",setup->doServo);   
  rampinfoPtr += SERV_STR_LEN;

  sprintf(rampinfoPtr,"#safbcard=%s\n",setup->safbCard);   
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#sq2fbcard=%s\n",setup->sq2fbCard);   
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#sq2biascard=%s\n",setup->sq2biasCard);   
  rampinfoPtr += SERV_STR_LEN;
  sprintf(rampinfoPtr,"#sq1biascard=%s\n",setup->sq1biasCard);   

  rampinfoPtr += SERV_STR_LEN;
  sc2daservodataheadarrayVal(rampinfoPtr,setup->initFB,"#intFB=",
                             COL_NUM,INT_NUMBER,status);

  rampinfoPtr +=(SERV_STR_LEN*4);
  sc2daservodataheadarrayVal(rampinfoPtr,setup->cableOffset,"#cableOff=",
                             COL_NUM,INT_NUMBER,status);

  rampinfoPtr +=(SERV_STR_LEN*4);
  sc2daservodataheadarrayVal(rampinfoPtr,setup->biaslckPt,"#biasLck=",
                             COL_NUM,INT_NUMBER,status);

  rampinfoPtr +=(SERV_STR_LEN*4);
  sc2daservodataheadarrayVal(rampinfoPtr,setup->slopSelect,"#slopSelect=",
                             COL_NUM,INT_NUMBER,status);

  rampinfoPtr +=(SERV_STR_LEN*4);
  sc2daservodataheadarrayVal(rampinfoPtr,setup->colMask,"#colMask=",
                             COL_NUM,INT_NUMBER,status);

  rampinfoPtr +=(SERV_STR_LEN*4);
  sc2daservodataheadarrayVal(rampinfoPtr,setup->rowMask,"#rowMask=",
                             ROW_NUM,INT_NUMBER,status);

  rampinfoPtr +=(SERV_STR_LEN*4);
  sc2daservodataheadarrayVal(rampinfoPtr,setup->cableSlope,"#cableSlope=",
                             COL_NUM,FLOAT_NUMBER,status);

  rampinfoPtr +=(SERV_STR_LEN*4);
  sc2daservodataheadarrayVal(rampinfoPtr,setup->cableScale,"#cableScale=",
                             COL_NUM,FLOAT_NUMBER,status);
}



/**
 * \fn void sc2daservodataheadarrayVal(char *infoPtr, void *paramarray,
 *  char *name, int howMany, int flag, StatusType *status)
 *
 * \brief function:
 *  populate array struture
 *
 * \param  infoPtr     char pointer for servoData address
 * \param  paramarray  void pointer for setup-parameter array
 * \param  name       char point for the param name
 * \param  howMany    int how many param need to read in this call
 * \param  flag       int  FLOAT_NUMBER, INT_NUMBER 
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  int    j, *intArray;
  double  *floatArray;
  char   tmp[CMD_LEN];

  if (*status != STATUS__OK) return;

  sprintf(infoPtr,"%s",name);
  if (flag==INT_NUMBER)
  {
     intArray=(int*)paramarray;
    for ( j=0; j<howMany; j++ )
    {
      sprintf (tmp, " %d", intArray[j] );
      strcat(infoPtr,tmp); 
    }
  }
  else
  {
    floatArray=(double*)paramarray;
    for ( j=0; j<howMany; j++ )
    {
      sprintf (tmp, " %4.3E", floatArray[j] );
      strcat(infoPtr,tmp); 
    }
  }
  strcat(infoPtr,"\n");
}



/**
 * \fn void sc2daservodataheadarrayfloatVal(char *infoPtr, double *paramarray,
 *  char *name, int howMany, StatusType *status)
 *
 * \brief function:
 *  populate array struture
 *
 * \param  infoPtr     char pointer for servoData address
 * \param  paramarray double pointer for setup-parameter array
 * \param  name       char point for the param name
 * \param  howMany    int how many param need to read in this call
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2daservodataheadarrayfloatVal
*/
void sc2daservodataheadarrayfloatVal
(
char       *infoPtr,
double     *paramarray,
char        *name,
int         howMany,
StatusType *status
)
{
  int    j;
  char   tmp[CMD_LEN];

  if (*status != STATUS__OK) return;

  sprintf(infoPtr,"%s",name);
  for ( j=0; j<howMany; j++ )
  {
    sprintf (tmp, " %f", paramarray[j] );
    strcat(infoPtr,tmp); 
  }
  strcat(infoPtr,"\n");
}



/**
 * \fn void sc2daservodatamemInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  ARRAYSET *arrayset, char **dataBuf, StatusType *status)
 *
 * \brief function 
 *  allocate the maximum memory for holding the frame data for finding max 
 *  modulation 
 *  not include MCE frame header, dark row 
 *  SERVO_DATAHEAD:string information for ramping bias, fb, card, intfb .. etc
 *  setting
 *   tmpdataBuf:  SERVO_DATAHEAD string
 *                data of all Channels or all frames
 *                BIAS_PEAK*biasNo
 *                MAX_P2P 
 *
 * \param con       SDSU context structure
 * \param myInfo    dasInfo structure pointer
 * \param arrayset  ARRAYSET structure pointer
 * \param dataBuf   address of pointer for servo data
 * \param status    StatusType     
 *
 */
/*+ sc2daservodatamemInit
*/
void sc2daservodatamemInit
(
SDSU_CONTEXT     *con,      
dasInfoStruct_t  *myInfo,    
ARRAYSET         *arrayset,
char             **dataBuf,
StatusType       *status
)
{
  int        frameSize, headSize, totalptsSize;
  int        peakSize, maxpeakSize, subtotalSize,totalSize;
  int        cableadjSize, meanvalSize;
  int        heatModulate,pixels;
  int        servobindataSize=0;

  if (*status != STATUS__OK) return;

  myInfo->servobindataPtr=NULL;
  if  (myInfo->dataFormat==MCE_BINARY_FORM)
  {
     // allocate memory  
     //     SERVOBIN_DATAHEAD
     //     data_0*COL_NUM     |    |    |
     //     xxxx_0*COL_NUM     | 4  |    |
     //     data_1*COL_NUM     |    |    |
     //     xxxx_1*COL_NUM     |    |    |
     //     .........               |    |
     //     .........               |    |
     //     .........             fdbkNo |
     //     .........               |    |
     //     .........               |    |
     //     data_0*COL_NUM     |    |    |
     //     xxxx_0*COL_NUM     | 4  |    |
     //     data_1*COL_NUM     |    |    |
     //     xxxx_1*COL_NUM     |    |    |
     //=============================  biasNo
     //     data_0*COL_NUM     |    |    |
     //     xxxx_0*COL_NUM     | 4  |    |
     //     data_1*COL_NUM     |    |    |
     //     xxxx_1*COL_NUM     |    |    |
     //     .........               |    |
     //     .........               |    |
     //     .........             fdbkNo |
     //     data_0*COL_NUM     |    |    |
     //     xxxx_0*COL_NUM     | 4  |    |
     //     data_1*COL_NUM     |    |    |
     //     xxxx_1*COL_NUM     |    |    |

     servobindataSize=sizeof(SERVOBIN_DATAHEAD) + 
           sizeof(int)*arrayset->fdbkNo*COL_NUM*4*arrayset->biasNo;
     myInfo->servobindataPtr=calloc(servobindataSize, 1);
  }

  heatModulate=arrayset->fdbkNo*4+1;  
  arrayset->totalPixel=ROW_NUM*COL_NUM;
  
  if( arrayset->servo==CABLECAL)
  {
    totalptsSize= arrayset->cableNo*arrayset->fdbkNo*sizeof(int);
    peakSize = sizeof(BIAS_PEAK)*arrayset->cableNo;
  }
  else if (arrayset->servo==TESTRANSIT)
  {
    // for tesTransit, we need all row rather than one row

    totalptsSize= arrayset->biasNo*ROW_NUM*heatModulate*sizeof(int);
    peakSize = sizeof(BIAS_PEAK)*arrayset->biasNo;
  }
  else 
  {
    totalptsSize= arrayset->biasNo*arrayset->fdbkNo*sizeof(int);
    peakSize = sizeof(BIAS_PEAK)*arrayset->biasNo;
  }
  headSize = sizeof(SERVO_DATAHEAD); 
  frameSize = totalptsSize*COL_NUM;
  maxpeakSize = sizeof(MAX_P2P);
  subtotalSize=  headSize + frameSize + peakSize + maxpeakSize;

  cableadjSize =arrayset->biasNo*COL_NUM*sizeof(int);
  meanvalSize =  cableadjSize;
  totalSize =  subtotalSize + cableadjSize + meanvalSize;

  arrayset->servodataPtr=calloc(totalSize, 1);
  if(arrayset->servodataPtr ==NULL)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
        "sc2daservodatamemInit: failed to allocate space for servodataPtr");
    return;
  }

  *dataBuf=arrayset->servodataPtr;

  arrayset->meanvalPtr =arrayset->servodataPtr + (totalSize -meanvalSize);
  arrayset->cableadjPtr=arrayset->servodataPtr + (totalSize -cableadjSize-meanvalSize);
 
  // allocate memory to hold all ch's data per bias, use the largest value 
  // to allocate if (arrayset->servo ==TESTRANSIT) 
  arrayset->allData=(int *)calloc(heatModulate*COL_NUM,sizeof(int));
  if (arrayset->allData==NULL)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
        "sc2daservodatamemInit: failed to allocate space for allData");
    return;
  }
  // allocate memory to hold each ch's data per bias
  //if (arrayset->servo ==TESTRANSIT)
  arrayset->chData=(double *)calloc(heatModulate*COL_NUM,sizeof(double));
  if (arrayset->chData==NULL)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
        "sc2daservodatamemInit: failed to allocate space for chData");
    return;
  }

  // check if TESTRANSIT has larger value
  // COL_NUM*heatModulate always <subtotalSize )
  arrayset->filtedPtr=calloc(subtotalSize, 1);
  if(arrayset->filtedPtr ==NULL)
  {
    *status=DITS__APP_ERROR;
    return;
  }
  // allocate memory to hold (folded-shifted) each ch's data
  arrayset->convolPtr=(double *)calloc( FILTER_ORDER,sizeof(double));
  if (arrayset->convolPtr==NULL)
  {
    *status=DITS__APP_ERROR;
    return;
  }
  // allocate memory to hold filter impulse
  arrayset->impulsePtr=(double *)calloc( FILTER_ORDER,sizeof(double));
  if (arrayset->impulsePtr==NULL)
  {
    *status=DITS__APP_ERROR;
    return;
  }

  if (arrayset->servo ==TESTRANSIT) 
  {
    pixels=arrayset->biasNo*arrayset->totalPixel;
    // allocate memory to hold each pixel's transitVal
    arrayset->transitVal=(double *)calloc(pixels,sizeof(double));
    if (arrayset->transitVal==NULL)
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,
          "sc2daservodatamemInit: failed to allocate space for transitVal");
      return ;
    }

    // allocate memory to hold each heater Value for the modulation
    arrayset->heater=(int *)calloc(heatModulate,sizeof(int));
    if (arrayset->heater==NULL)
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,
          "sc2daservodatamemInit: failed to allocate space for heaterValue");
      return ;
    }

    // allocate memory to hold each pixel's transitFlag
    arrayset->transitFlag=(int *)calloc(pixels,sizeof(int));
    if (arrayset->transitFlag==NULL)
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,
          "sc2daservodatamemInit: failed to allocate space for transitFlag");
      return ;
    }

    // allocate memory to hold each pixel's lockFlag
    arrayset->lockFlag=(int *)calloc(pixels,sizeof(int));
    if (arrayset->lockFlag==NULL)
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,
          "sc2daservodatamemInit: failed to allocate space for lockFlag");
      return;
    }
  }
}


/**
 * \fn void sc2daservobindataoutInit(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *     int colrow,StatusType *status)
 *
 * \brief function
 *  write all servo settings into memory for other program/scripts
 *  to process data
 *
 * \param  myInfo    dasInfoStruct_t poiter
 * \param  setup      ARRAYSET structure pointer
 * \param  colrow     int if sq1biasservo, pass column if testransit pass row
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2daservobindataoutInit
*/
void sc2daservobindataoutInit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,  
int             colrow,   
StatusType      *status
)
{
  time_t   tm;
  int      i;
  char     msg[FILE_LEN];
  SERVOBIN_DATAHEAD *binHead;

  if (!StatusOkP(status)) return;

  // change the datafile name to *.bin, it is already opened 
  if ( myInfo->fpData !=NULL)
  {
     fclose(myInfo->fpData);
     sprintf(msg,"rm -f %s",myInfo->dataFile);
     system(msg);   

     sprintf(msg,"%s.bin",myInfo->dataFile);
     if((myInfo->fpData=fopen64(msg,"w")) == NULL )
       {
	 *status = DITS__APP_ERROR;
	 ErsRep (0, status, "sc2daservobindataoutInit: Error- failed to open file %s", msg); 
	 return;
       }
  }

  binHead =(SERVOBIN_DATAHEAD *)myInfo->servobindataPtr;

  tm = time(NULL);
  sprintf (binHead->started,"%s\n",asctime(localtime(&tm)) );
  sprintf (binHead->servo, "%s\n", setup->servoName); 
  for (i=0; i<setup->totalRC2use; i++)
  { 
    sprintf (msg, "%d ",setup->whichRC[i]);
    strcat(binHead->card,msg); 
  }
  strcat (binHead->card, "\n"); 

  binHead->bias=setup->minBIAS; 
  binHead->bstep=setup->stepBIAS; 
  binHead->nbias= setup->biasNo; 
  binHead->feed=setup->minFDBK; 
  binHead->fstep=setup->stepFDBK; 
  binHead->nfeed=setup->fdbkNo; 
  binHead->nrow=setup->selRow;
  binHead->doservo=setup->doServo;
  for (i=0;i<COL_NUM;i++)
  {
    binHead->cableoffset[i]=setup->cableOffset[i];
    binHead->cableScale[i]=setup->cableScale[i];
    binHead->cableSlope[i]=setup->cableSlope[i];
    binHead->cableadjThd[i]=setup->cableadjThd[i];
    binHead->cableadjScale[i]=setup->cableadjScale[i];
    binHead->slopselect[i]=setup->slopSelect[i];
    binHead->initFB[i]=setup->initFB[i];
    binHead->sabiaslck[i]=setup->biaslckPt[i];
    binHead->gain[i]=setup->gain[i];
    binHead->zfact[i]=setup->zfact[i];
    binHead->fluxPeriod[i]=setup->fluxPeriod[i];
    binHead->colMask[i]=setup->colMask[i];
  }
  for (i=0;i<ROW_NUM;i++)
  {
    binHead->rowMask[i]=setup->rowMask[i]; 
  }
  binHead->row=colrow;
}



/**
 * \fn void sc2daservoFstGO(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  int startSeq, int endSeq, struct mcexml_struct  *mceInxpt, 
 *  char *dateTime, StatusType *status)
 *
 * \brief function: 
 *  send a single GO to MCE at the beginning of SERVO action
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param startSeq  int: 
 * \param endSeq    int:
 * \param mceInxpt  struct mcexml_struct pointer
 * \param dateTime  dateTime string pointer         
 * \param status    StatusType.  given and return
 *
 */
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
)
{
  if (!StatusOkP(status)) return;

  sc2dalib_frametakeInit(con,myInfo,startSeq,endSeq,mceInxpt,dateTime,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status, "sc2daservoFstGO: sc2dalib_frametakeInit failed");
    return;
  }
  con->process.framesetup=DA_MCE_SINGLEDATA;  
}


/**
 * \fn void sc2daservoEnd(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, char *dataBuf, char *ssalckdataBuf, int flag,
*   StatusType *status)
 *
 * \brief function:
 *  end servo action
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param setup    ARRAYSET structure pointer
 * \param dataBuf   char pointer for the servo data
 * \param ssalckdataBuf   char pointer for ssalck data
 * \param flag      int  1: check IN_SEQ, 0, don't
 * \param status    StatusType.  given and return
 *
 */
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
)
{
  jitDebug(16,"sc2daservoEnd \n");
  if (myInfo->fpData !=NULL)
    fprintf (myInfo->fpData,"</DATA>\n" );

  if(dataBuf!=NULL)
  {  free(dataBuf);
    jitDebug(16,"sc2daservoEnd free dataBuf\n");
  }

  if(ssalckdataBuf!=NULL)
    free(ssalckdataBuf);
  jitDebug(16,"sc2daservoEnd free sssalckdataBuf\n");
  
  if(setup->allData !=NULL)
    free(setup->allData);
  jitDebug(16,"sc2daservoEnd free setup->alldata\n");
  
  if(setup->chData !=NULL)
    free(setup->chData);
  jitDebug(16,"sc2daservoEnd free setup->chData\n");

  if(setup->filtedPtr !=NULL)
     free(setup->filtedPtr);
  jitDebug(16,"sc2daservoEnd free setup->filtedPtr\n");  

  if(setup->convolPtr !=NULL)
     free(setup->convolPtr);
  jitDebug(16,"sc2daservoEnd free setup->convolPtr\n");
    
  if(setup->impulsePtr !=NULL)
     free(setup->impulsePtr);
  jitDebug(16,"sc2daservoEnd free setup->impulsePtr\n");
 
 if(setup->transitFlag !=NULL)
     free(setup->transitFlag);
  jitDebug(16,"sc2daservoEnd free setup->transitFlag\n");
  
 if(setup->heater !=NULL)
     free(setup->heater);
  jitDebug(16,"sc2daservoEnd free setup->heater\n");
 
 if(setup->transitVal !=NULL)
     free(setup->transitVal);
  jitDebug(16,"sc2daservoEnd free setup->transitVal\n");

 if(setup->lockFlag !=NULL)
     free(setup->lockFlag);
  jitDebug(16,"sc2daservoEnd free setup->lockFlag\n");

  sc2dalib_actionfileEnd(con,myInfo,flag,status);
}



/**
 * \fn void sc2daservoInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  dasCmdInfo_t *myCmd, struct mcexml_struct *mceinxPtr, ARRAYSET *arrayset, 
 *  char *dateTime,  char **dataBuf, ARRAYSET *ssalckset,  char **ssalckBuf,
 *  StatusType *status)
 *
 * \brief function: 
 *  read in all args from the setup batch file for servo.
 *  allocate memory for holding the frame data (max 32*40), 
 *  open files and setup  initial values  for servo routine
 *
 * \param con           SDSU context structure
 * \param myInfo       dasInfo structure pointer
 * \param myCmd         dasCmdInfo_t pointer
 * \param mceinxPtr     mcexml_struct pointer 
 * \param arrayset      ARRAYSET structure pointer
 * \param dateTime      string pointer to dateTime string
 * \param dataBuf       address of pointer for servo data
 * \param ssalckset     ARRAYSET structure pointer
 * \param ssalckBuf     address of pointer for ssa data
 * \param status        StatusType     
 *
 */
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
)
{
  long       in_sequence,configured;
  long       range[]={0,40};
  char       servoFile[FILE_LEN];
  // the cmdselRow will overwrite setup.selRow, only for sq1
  long        cmdselRow, sq1Flag;     

  DitsArgType  argId;

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  sc2dalib_updateDebug(con,myInfo, status);

  SdpGeti("CONFIGURED", &configured, status);
  if (configured==0 )
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2daservoInit: the DA is not configured" ); 
    return;
  }
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2daservoInit: %s has not completed",
             seqStatus[myInfo->actionFlag]);
    return;
  } 
  myInfo->actionFlag=SERVOACTION;

  argId = DitsGetArgument();
  jitArgGetI( argId, "SVFILE_FLAG", 1, range, 0, 0,&myInfo->filesvFlag, status );
  if ( !StatusOkP(status) )
  {
    ErsOut(0,status, "sc2daservoInit: SVFILE_FLAGC is not given, use 2 default");
    myInfo->filesvFlag=2;
  }

  // filesvFlag=1 save cmdreply 
  // filesvFlag=2 don't save cmdreply, but save data
  // for servo, we only save data
  if(myInfo->filesvFlag==1)       myInfo->filesvFlag=2;

  jitArgGetS(argId,"DATA_FILE",2,NULL,"data.txt",0, FILE_LEN,
            myInfo->dataFile, NULL,status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2daservoInit: failed to get DATA_FILEE"); 
    return;
  }
  jitArgGetS (argId, "SETUPBATCH_FILE", 3,NULL, "sq2setup.txt", 0, FILE_LEN,
              myInfo->batchFile,  NULL,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2daservoInit: no SETUPBATCH_FILE ");
    return;
  }
  jitArgGetS (argId, "STRCHART_FILE", 4,NULL, "sq2str.txt", 0, FILE_LEN, 
              myInfo->strchartFile,  NULL,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2daservoInit: no STRCHART_FILE ");
    return;
  }

  // use cmdrepFile to save the lock points, 
  // not use it now, but use a fixed name: $datafile-lck
  jitArgGetS (argId, "LOCKPTS_FILE", 5,NULL, "lckpoints.txt", 0, FILE_LEN, 
              myInfo->cmdrepFile,  NULL,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2daservoInit: no LOCKPTS_FILE ");
    return;
  }
  sprintf( myInfo->cmdrepFile,"%s-lck", myInfo->dataFile);

  jitArgGetI( argId, "SQ1_FLAG", 6, range, 0, 0, &sq1Flag, status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status, "sc2daservoInit: SQ1_FLAG is not given");
    return;
  }
  jitArgGetI( argId, "SEL_ROW", 7, range, 0, 0, &cmdselRow, status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status, "sc2daservoInit: SEL_ROW is not given");
    return;
  }

  arrayset->allData=NULL;
  arrayset->chData=NULL;
  arrayset->servodataPtr=NULL;
  arrayset->filtedPtr=NULL;
  arrayset->convolPtr=NULL;
  arrayset->impulsePtr=NULL;
  arrayset->transitVal=NULL;
  arrayset->heater=NULL;
  arrayset->transitFlag=NULL;
  arrayset->lockFlag =NULL;
  myInfo->fpSq1=NULL;
  fluxjmpFp=NULL;

  sc2daservosetupInit(myInfo,arrayset,status);


  // readback the sample_num to re-scale the data and dataMode 
  sc2dareadmceparVal(con,myInfo,arrayset,mceinxPtr,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status,"sc2daservoInit: failed to call sc2dareadmceparVal");
    return;
  }


  //expand all include file into one single file $SC2SCRATCH/tmp
  // assign it to myInfo->batchFile
  sc2dalibsetup_servoreadsetupWrap(myInfo,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2daservoInit: sc2dalibsetup_servoreadsetupWrap failed"); 
    return;
  }


  sc2dalib_openFiles(myInfo,DATFILENOAPPEND,OTHERFILE,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2daservoInit: sc2dalib_openFiles failed"); 
    return;
  }

  jitDebug(16,"sc2daservoInit: result datafile (%s), setupfile (%s) \n",
      myInfo->dataFile, myInfo->batchFile);
      
  fprintf(myInfo->fpLog,"\n<%s> CMD from sc2da_Servo <%s>\n",
          dateTime,myInfo->dataFile);
  fflush(myInfo->fpLog);

  // read all setiings from myInfo->batchFile 
  sc2dalibsetup_servoreadSetup(myInfo,arrayset,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0, status,"sc2daservoInit:sc2dalibsetup_servoreadSetup failed ");
    return;
  }

  if(arrayset->totalRC2use==1)  // only one card
  {
    sprintf(myInfo->goCmd, "GO rc%d ret_dat 1",arrayset->whichRC[0]);
  }
  else  if(arrayset->totalRC2use==4)          /// four cards
  {
    sprintf(myInfo->goCmd, "GO rcs ret_dat 1");
  }
  else
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"sc2daservoInit: there are %d RC cards used, not four",
          arrayset->totalRC2use);
    ErsRep(0,status,"  need more info about how to send command to MCE");
    return;     
  }
  // over-write the row number for sq1
  if (sq1Flag)    arrayset->selRow=cmdselRow;

  
  if (arrayset->servo==SQ1BIASSERVO)
  {
    // cheating now
    arrayset->fdbkNo=arrayset->biasNo;
    arrayset->stepFDBK=arrayset->stepBIAS;
    arrayset->minFDBK=arrayset->minBIAS;
    arrayset->biasNo=1;
  }

  // add to read in sq1biasoptimal for sq1lock for later to set
  // biaslckPt[0..31]=sq1biasOpt[arrayset->selRow] in findinitlckPts()
  // X.Gao 2009.10.20   
  if (arrayset->servo==SQ1LOCK)
  {
    // read sq1biasoptimal in ->sq1biasOpt
    jitDebug(16,"_readOPT sq1biasoptimal.txt \n");
    sc2dalib_readOPT(myInfo,arrayset->sq1biasOpt,ROW_NUM,"sq1biasoptimal.txt",status);
    if ( !StatusOkP(status) )
      return;
  }

  jitDebug(16,"sc2daservoInit: sc2dalibsetup_servodataoutInit \n");
  if(myInfo->filesvFlag!=0 &&  myInfo->dataFormat !=MCE_BINARY_FORM)
     sc2dalibsetup_servodataoutInit(myInfo,arrayset,0,status);
 
  myInfo->msgwrtPt=0; 
  myInfo->msgreadPt=0;
  if(arrayset->servo==CABLECAL)
  {
    arrayset->biasNo=arrayset->cableNo;
    arrayset->stepBIAS=arrayset->stepCABLE;
    arrayset->minBIAS=arrayset->minCABLE;
  }
  if(arrayset->servo==BLACKBODY)
  {
    myInfo->cmdFlag=0;
  }

  // alloc all pointers 
  jitDebug(16,"sc2daservoInit: sc2daservodatamemInit \n");
  sc2daservodatamemInit(con,myInfo,arrayset,dataBuf,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,
        "sc2daservoInit:failed to allocate memory for holding frame data");
    return;
  }
  // this required by Mike to binary servo data for idl,this is default 
  if(myInfo->filesvFlag!=0 &&  myInfo->dataFormat ==MCE_BINARY_FORM)
     sc2daservobindataoutInit(myInfo,arrayset,0,status);

  strcpy(myCmd->mceCmd,myInfo->goCmd);
  sc2dalib_getcmdBuf(myInfo,myCmd,mceinxPtr,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2daservoInit: sc2dagetCmdBuf failed"); 
    return;
  }
 
  jitDebug(16,"sc2daservoInit: sc2daservoinitValue \n");
  sc2daservoinitValue(con,myInfo,arrayset,mceinxPtr,*dataBuf,dateTime,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0, status,"sc2daservoInit:sc2daservoinitValue failed ");
    return;
  }

  // need to pass divFlag to ssalckset too, so, move it here after initValue
  if(arrayset->servo==SQ2SERVO)
  {
    sc2daservossalckInit(con,myInfo,arrayset,ssalckset,ssalckBuf, status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2daservoInit: sc2daservossalckInit failed"); 
      return;
    }
  }
 
  if (arrayset->servo==SQ1OPEN || arrayset->servo==TESTRANSIT  ||
      arrayset->servo==SQ2OPEN || arrayset->servo==SQ2OPEN4P )
  {
    if (arrayset->servo==SQ1OPEN)
    {
      // read in all sq1 init fdbk int setup->pixeLock[ROW_NUM*COL_NUM]
      sc2dalibsetup_servosq1fdbkreadSetup(myInfo,arrayset,status);

      if ( !StatusOkP(status) )
      {
        ErsRep(0,status,"sc2daservoInit:sc2dalibsetup_servosq1fdbkreadSetup failed ");
        return;
      }
      jitDebug(16,"sc2daservoInit:sc2dalibsetup_servosq1fdbkreadSetup OK\n");
    }
    else if ( arrayset->servo==TESTRANSIT )
    {
      // when load sc2da task, stripchFlag=0;
      if(arrayset->slopSelect[10]==0)        arrayset->waveFlag=TRI_UP;
      else if(arrayset->slopSelect[10]==1)   arrayset->waveFlag=TRI_DOWN;
      else                                   arrayset->waveFlag=TRI_MID;
    }
    // also create datafile-binary for later usage
    sprintf(servoFile,"%s-binary",myInfo->dataFile);
    if((myInfo->fpSq1 = fopen64(servoFile,"a")) == NULL)
    { 
      *status=DITS__APP_ERROR;
      ErsRep(0, status,"sc2daservoInit:failed to open %s ", servoFile);
      return;
    }
    jitDebug(16,"save frame in Binary (%s) too\n",servoFile);
  }
  SdpPuti("IN_SEQUENCE",DA_MCE_SINGLEDATA,status);

  // servoFP is global pointer for storing servo data in ASCII format
  servoFp=NULL;
  sprintf (servoFile, "%sservodata.txt",myInfo->dataFile);
  if((servoFp = fopen(servoFile,"w")) == NULL)
  { 
    *status=DITS__APP_ERROR;
    ErsRep(0, status,"sc2daservoInit 2:failed to open %s ", servoFile);
    return;
  }


}


/**
 * \fn void sc2daservosetupInit(dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,  StatusType *status)
 *
 * \brief function: 
 *  initial all setting to 0
 *
 * \param myInfo       dasInfo structure pointer
 * \param setup        ARRAYSET structure pointer
 * \param status        StatusType     
 *
 */
/*+ sc2daservosetupInit
*/
void sc2daservosetupInit
(
dasInfoStruct_t       *myInfo,    
ARRAYSET              *setup,
StatusType            *status
)
{
  int i;

  if (*status != STATUS__OK) return;

  setup->minBIAS=0;
  setup->stepBIAS=0; 
  setup->biasNo=1; 
  setup->minFDBK=0; 
  setup->stepFDBK=0; 
  setup->fdbkNo=1; 
  setup->selRow=0;
  setup->doServo=0;
  for (i=0;i<COL_NUM;i++)
  {
    setup->cableOffset[i]=0;
    setup->cableScale[i]=0;
    setup->cableSlope[i]=0;
    setup->cableadjThd[i]=0;
    setup->cableadjScale[i]=0;
    setup->slopSelect[i]=0;
    setup->initFB[i]=0;
    setup->biaslckPt[i]=0;
    setup->gain[i]=0;
    setup->zfact[i]=0;
    setup->fluxPeriod[i]=0;
    setup->colMask[i]=0;
  }
  for (i=0;i<ROW_NUM;i++)
  {
    setup->rowMask[i]=0; 
  }
}


/**
 * \fn void sc2daservossalckInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *   ARRAYSET *arrayset, ARRAYSET *ssalckset, char **ssalckBuf,
 *   StatusType *status)
 *
 * \brief function: 
 *  set initial values for doing ssalock during sq2servo
 *
 * \param con           SDSU context structure
 * \param myInfo       dasInfo structure pointer
 * \param arrayset      ARRAYSET structure pointer
 * \param ssalckset     ARRAYSET structure pointer
 * \param ssalckBuf     address of pointer for ssa data
 * \param status        StatusType     
 *
 */
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
)
{
  char     tmpfile[128];

  if (*status != STATUS__OK) return;

  // copy all then change some, see findlockpoints for example
  memcpy(ssalckset,arrayset,sizeof(ARRAYSET));
  ssalckset->biasNo=1;
  ssalckset->stepBIAS=0;
  ssalckset->servo=SSALOCK;
  strcpy(ssalckset->servoName,"sq2ssalock");
  ssalckset->doServo=SQ2SERVOSSALOCK;
  // open fpOtheruse for ssaopen lock points during sq2servo
  sprintf (tmpfile, "%s-sq2ssalck",myInfo->dataFile );       
  if((myInfo->fpOtheruse=fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2daservossalckInit: Error- failed to open file %s", tmpfile); 
      return;
    }
  sc2daservodatamemInit(con,myInfo,ssalckset,ssalckBuf,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,
          "sc2daservossalckInit:failed to allocate memory for holding frame data");
    return;
  }
  sc2daservodataheadInit(ssalckset,*ssalckBuf,status);
}



/**
 * \fn void sc2daservoinitValue(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, struct mcexml_struct *mceinxPtr, char *servoData, 
 *  char *dateTime, StatusType *status)
 *
 * \brief function: 
 *  setup initial values for ramping at SSA, SQ2 and SQ1 
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param setup     ARRAYSET structure pointer
 * \param mceinxPtr  mcexml_struct pointer
 * \param servoData char pointer for servo data
 * \param dateTime  char pointer for date Time
 * \param status    StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 */
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
)
{
  int  i;

  if (!StatusOkP(status)) return;
  
  // copy divbynFlag to setup
  setup->divFlag=(int)myInfo->divbynFlag;

  fprintf(myInfo->fpLog,"\n<%s> %s servoInit \n",dateTime, setup->servoName);
 
  for (i=0;i<COL_NUM;i++)
  {
    // test to see if it work, cableinitVal from SSA max Modulation
     if( setup->servo==CABLECORET || setup->servo== SSALOCK) 
        setup->cableadjVal[i]=setup->cableadjInit[i];
     else
        setup->cableadjVal[i]=0;

     if (setup->servo==SSALOCK)
       setup->ssabiaslckPt[i]=setup->biaslckPt[i];

     setup->nomodulCount[i]=0;
  }
  
  if (setup->servo==SSARAMP  || setup->servo==SSALOCK    ||
      setup->servo==SSARAMP1 || setup->servo==SQ2BIASING || 
      setup->servo==CABLECAL || setup->servo==CABLECORET 
     ) 
  {
    jitDebug(2," call sc2daservoinitSSSA \n");
    sc2daservoinitSSA(con,myInfo,setup,mceinxPtr,servoData,dateTime, status);
    if (!StatusOkP(status)) 
    {
      ErsRep(0,status,"sc2daservoinitValue: failed to call sc2daservoinitSSA");
      return;
    }
  }
  else if (setup->servo==SQ2SERVO || setup->servo==SQ2OPEN ||
           setup->servo==SQ2LOCK  || setup->servo==SQ2OPEN4P )
  {
    jitDebug(16," call sc2daservoinitSQ2\n"); 
    sc2daservoinitSQ2(con,myInfo,setup,mceinxPtr, servoData,dateTime,status);
    if (!StatusOkP(status)) 
    {
      ErsRep(0,status,"sc2daservoinitValue: failed to call sc2daservoinitSQ2");
      return;
    }
  }
  else if (setup->servo==SQ1SERVO      || setup->servo==SQ1OPEN      ||
           setup->servo ==SQ1LOCK      || setup->servo ==HEATERSERVO ||  
           setup->servo ==TESBIASSERVO || setup->servo==SQ1BIASSERVO || 
           setup->servo==TESTRANSIT 
          )
  {   
    #ifndef SDSU
    jitDebug(16," call sc2daservoinitSQ1 \n");
    sc2daservoinitSQ1(con,myInfo,setup,mceinxPtr, servoData,dateTime,status);
    if (!StatusOkP(status)) 
    {
      ErsRep(0,status,"sc2daservoinitValue: failed to call sc2daservoinitSQ1");
      return;
    }
    #else
    MsgOut(status,"set hardware as SDSU, don't call sc2daservoinitSQ1");
    #endif
  }
  else
  {
     MsgOut(status,"sc2daservoinitValue: do nothing for %s",setup->servoName);
  }
  sc2daservodataheadInit(setup,servoData,status);

}




/**
 * \fn void sc2dareadmceparVal(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, struct mcexml_struct *mceinxPtr, StatusType *status)
 *
 * \brief function: 
 *  read back MCE sample_num and data_mode parameters 
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param setup     ARRAYSET structure pointer
 * \param mceinxPtr  mcexml_struct pointer
 * \param status    StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 */
/*+ sc2dareadmceparVal
*/
void sc2dareadmceparVal
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
ARRAYSET              *setup,
struct mcexml_struct  *mceinxPtr,
StatusType            *status
)
{
  // 12 Dec 2006 X.Gao
  // MCE firmware only allows rb rc1, rc2 rcs,  not rb sys
  static char  sampleCmd[]="rb rc1 sample_num 1 ";
  static char  datamodeCmd[]="rb rcs data_mode 1 ";

  if (!StatusOkP(status)) return;

  sc2dalib_readmceVal(con,myInfo,mceinxPtr,sampleCmd,&setup->sampleNo,1,status);
  if (!StatusOkP(status)) 
    return;
  else
  {
    sprintf(errmsg,"MCE's sample_num= %d", setup->sampleNo);   
    //MsgOut(status,errmsg);
    //printf("%s\n",errmsg);
  }

  if (setup->sampleNo ==0)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"sc2dareadmceparVal: sampleNo=0 ");
    return;
  }
  sc2dalib_readmceVal(con,myInfo,mceinxPtr,datamodeCmd,&setup->dataMode,1,status);
  if (!StatusOkP(status)) 
    return;
  else
  {
    sprintf(errmsg," MCE's data_mode= %d", setup->dataMode);   
    //MsgOut(status,errmsg);
    //printf("%s\n",errmsg);
  }
}



/**
 * \fn void sc2daservoinitSSA(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, struct mcexml_struct *mceInxpt,
 *  char *servoData,char *dateTime, StatusType *status)
 *
 * \brief function
 *  setup initial values for starting SSARAMP, SSALOCK, SQ2BIASING 
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param setup     ARRAYSET structure pointer
 * \param mceInxpt  mcexml_struct pointer
 * \param servoData char pointer for servo data
 * \param dateTime  char pointer for date Time
 * \param status    StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{

  if (!StatusOkP(status)) return;

  setup->fdbk=setup->minFDBK;

  if (setup->servo==SSARAMP || setup->servo==SSARAMP1 || setup->servo==CABLECAL )
  {
    setup->bias=setup->minBIAS;
    setup->biasFlag=BIASCHANGED;
  }
  else
    setup->biasFlag=BIASNOCHANGE;
  
  if  (setup->servo==SSALOCK  || setup->servo==CABLECORET)
    setup->biasNo=1;

  sc2dasetbiasfdbkSSA(con,myInfo,setup,mceInxpt,dateTime,0,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status,"sc2daservoinitSSA: failed to call sc2dasetbiasfdbkSSA");
    return;
  }
}


/**
 * \fn void sc2daservoinitSQ2(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, struct mcexml_struct *mceInxpt,
 *  char *servoData, char *dateTime,StatusType *status)
 *
 * \brief function
 *  setup initial values for starting SQ2SERVO SQ2OPEN Sq2LOCK
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param setup     ARRAYSET structure pointer
 * \param mceInxpt  mcexml_struct pointer
 * \param servoData char pointer for servo data
 * \param dateTime  char pointer for date Time
 * \param status    StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  int    i;
  char   tmp[CMD_LEN];
  dasCmdInfo_t  fbCmd;
  
  if (!StatusOkP(status)) return;

  setup->fdbk=setup->minFDBK;

  if (setup->servo==SQ2SERVO  && setup->stepBIAS != 0 )
  {
    setup->bias=setup->minBIAS;
    setup->biasFlag=BIASCHANGED;
  }
  else
    setup->biasFlag=BIASNOCHANGE;

  if  (setup->servo==SQ2LOCK || setup->servo==SQ2OPEN4P || setup->servo==SQ2OPEN)
    setup->biasNo=1;

  for (i=0;i<COL_NUM;i++) 
    setup->fluxPeriod[i]=0;
    
  if ( setup->servo==SQ2OPEN4P )
  {
    //sq2fb =0 and wait a second
    sprintf (fbCmd.mceCmd, "wb %s flux_fb",setup->sq2fbCard ); 
    for ( i=0; i<COL_NUM; i++ )
    {
      sprintf ( tmp, " %d", setup->fdbk );
      strcat(fbCmd.mceCmd,tmp); 
    }
    sc2dalib_sendCmd(con,myInfo,&fbCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,"sc2daservoinitSQ2: _sendCmd failed"); 
      return;
    }   
    sleep(WAIT_TIME_AFTER_BIAS * 2);
  }
  else
  {
    if (setup->servo==SQ2SERVO)
      jitDebug(2,"sc2daservoinitSQ2: ready for ramp SSA FB for new lock point \n");
    sc2dasetbiasfdbkSQ2(con,myInfo,setup,mceInxpt,dateTime,0,status);
    if (!StatusOkP(status)) 
    {
      ErsRep(0,status,"sc2daservoinitSQ2: failed to call sc2dasetbiasfdbkSQ2");
      return;
    }
  }
} 


/**
 * \fn void sc2daservoinitSQ1(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  ARRAYSET *setup, struct mcexml_struct *mceInxpt,
 *  char *servoData,char *dateTime,StatusType *status)
 *
 * \brief function:
 *  setup initial values for starting SQ1 
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param setup     ARRAYSET structure pointer
 * \param mceInxpt  mcexml_struct pointer
 * \param servoData char pointer for servo data
 * \param dateTime  char pointer for date Time
 * \param status    StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  if (!StatusOkP(status)) return;

  setup->fdbk=setup->minFDBK;
  setup->biasFlag=BIASNOCHANGE;

  if (setup->servo==TESTRANSIT )
  {
    // only set heater as setup->fdbk=setup->minFDBK
    setup->bias=setup->minBIAS;
    sc2datransitTES(con,myInfo,setup,mceInxpt,dateTime,status);
    if (!StatusOkP(status)) 
    {
      ErsRep(0,status, "sc2daservoinitSQ1: failed to call sc2datransitTES");
    }
    return;
  }
  // the rest  
  if (setup->servo==SQ1OPEN || setup->servo ==SQ1LOCK )
  {
      setup->biasNo=1;  
  }
  else if( setup->servo ==SQ1SERVO || setup->servo ==SQ1BIASSERVO )
  {
    setup->biasFlag=BIASCHANGED;

    if (setup->servo==SQ1SERVO)
      setup->bias=setup->minBIAS;

    // for SQ1BIASSERVO, we have already swaped fdbkNo(biasNo)
    // stepFDBK(stepBIAS), and we ramp fdbk, but still do initial  
  }

  sc2dasetbiasfdbkSQ1(con,myInfo,setup,mceInxpt,dateTime,0,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status,"sc2daservoinitSQ1: failed to call sc2dasetbiasfdbkSQ1");
    return;
  }
 }


/**
 * \fn void sc2daservoprintMsg(ARRAYSET *setup, FILE *fp, int flag,
 *  void *paramarray, char *name, StatusType *status)
 *
 * \brief function
 *  print some meaasge from array structure
 *
 * \param  setup    ARRAYSET structure pointer 
 * \param  fp       FILE pointer
 * \param  flag      int  FLOAT_NUMBER, INT_NUMBER 
 * \param  paramarray void pointer for setup-parameter array
 * \param  name       char point for the param name
 * \param  status   StatusType pointer  G&R
 *
 */
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
)
{
  int i, ch;
  int   *intArray;
  double *floatArray;

  fprintf(fp,"%s \n",name);
  if(flag==INT_NUMBER)
  {
    intArray=(int*)paramarray;
    for (i=0;i<4;i++)
    {
      for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
        fprintf(fp,"%11d ",  intArray[ch]);
      fprintf(fp,"\n");
    }
  }
  else
  {
    floatArray=(double*)paramarray;
    for (i=0;i<4;i++)
    {
      for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
        fprintf(fp,"%11.4f ", floatArray[ch]);
      fprintf(fp,"\n");
    }
  }                                                             
} 


/**
 * \fn void sc2daservoFunc(SDSU_CONTEXT *con, char * byte, 
 *  dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *  mcexml_struct *mceInxpt, int *isEnd, char *dataBuf, char *dateTime,
 *  ARRAYSET *ssalckset,  char *ssalckBuf, StatusType *status)
 *
 * \brief function
 *  start to ramp feedback and bias, apply servo algorithm
 *  depending on servo=xxxx
 *
 * \param con      SDSU_CONTEXT pointer
 * \param byte     current raw frame data pointer
 * \param myInfo  dasInfoStruct_t poiter
 * \param setup    ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param isEnd    end flag     given & return
 * \param dataBuf  char pointer for servo data
 * \param dateTime char pointer for date Time string
 * \param ssalckset     ARRAYSET structure pointer
 * \param ssalckBuf     char pointer for ssa data
 * \param status   StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.

 * in sq2servo, for each sq2bias,ssalckFunc first finds initFB[ch], zfact[ch],
 * and gain[ch] and update them for sq2servo to use
 *
 * the question is the latest ssafb found from ssafb/sq2fb curve 
 * obtained from sq2servo data will be different from ssalock.
 * do we need to do a separate ssalock to update sq2open.txt? 
 *
 */
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
)
{
  char        msg[FILE_LEN];
  char       *peakinfoPtr;
  int         i;
  int         frameSize, headSize, totalptsSize;
  BIAS_PEAK   *peakperbiasPtr;

  static int  biasInx,fdbkInx;
  static int  ptrFlag;       // flag for initialise  pointer 
  static int  sq2ssalckFlag;       // flag for initialise  pointer 
  static BIAS_PEAK    *peakbiasPtr;
  
  if (!StatusOkP(status)) return;

  if(*isEnd==1)
  {
    *isEnd=0;
    sq2ssalckFlag=1;
    biasInx=fdbkInx=ptrFlag=0;
    setup->biasFlag=BIASCHANGED;

    headSize = sizeof(SERVO_DATAHEAD); 
    totalptsSize= setup->biasNo*setup->fdbkNo*sizeof(int);
    frameSize = totalptsSize*COL_NUM;
    peakinfoPtr  = (dataBuf+ headSize+frameSize);

    // point to entry of array of BIAS_PEAK struture
    // for sc2daservossalckFunc to use
    peakbiasPtr =(BIAS_PEAK *)peakinfoPtr;
    peakperbiasPtr=peakbiasPtr + biasInx;
    for(i=0;i<COL_NUM;i++)
    {
      peakperbiasPtr->initFB[i]=setup->initFB[i];
      peakperbiasPtr->zfact[i]=setup->zfact[i];
      peakperbiasPtr->gain[i]=setup->gain[i];
    }
    fprintf(myInfo->fpLog," %s: servoFunc \n",setup->servoName);

  }

  // during servoInit  ssalckset->doServo =SQ2SERVOSSALOCK, 
  // inside ssalckFunc, sq2bias is set
  if(setup->servo==SQ2SERVO && ssalckset->doServo==SQ2SERVOSSALOCK)
  {
     //peakperbiasPtr=peakbiasPtr + biasInx;
    // re adjust ssa initFB Zfactor and Gain mainly due to sq2bias change
    sc2daservossalckFunc(con,byte,myInfo,setup,ssalckset,mceInxpt,ssalckBuf,
             dateTime,&sq2ssalckFlag,peakbiasPtr,biasInx,dataBuf,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,"sc2daservoFunc: sc2daservossalckFunc failed"); 
      return;
    }
  }
  else
  {
    while (biasInx < setup->biasNo)
    {
      if( fdbkInx==0)
      {
        sprintf(msg,"%s biasInx=<%d>.....",setup->servoName,biasInx);
        //  =0, use MsgOut, =1 use ErsOut, =2 ErsRep  3, printf, >3 no print
        sc2dalib_msgprintSave(myInfo,"servo:%s",msg,NO_PRINT,status);
        if (setup->servo==SQ2SERVO)
        {
          sc2daservoprintMsg(setup, myInfo->fpLog, INT_NUMBER,setup->initFB,"initFB",status);
          sc2daservoprintMsg(setup, myInfo->fpLog, INT_NUMBER,setup->zfact,"Zfact",status);
          sc2daservoprintMsg(setup, myInfo->fpLog, FLOAT_NUMBER,setup->gain,"Gain",status);
        }
      }
      if( setup->servo==CABLECAL)
      {
        for ( i=0; i<COL_NUM; i++ )
          setup->cableOffset[i]=setup->minCABLE + biasInx*setup->stepCABLE; 
      }

      // don't apply bias for these servo, sq1biasservo use fdbk as bias, so, 
      if ( setup->servo==SSALOCK      || setup->servo==SQ1OPEN      || 
           setup->servo==SQ2OPEN      || setup->servo==SQ1LOCK      ||
	   setup->servo==SQ2LOCK      || setup->servo==SQ1BIASSERVO ||
           setup->servo==SQ2OPEN4P    || setup->servo==CABLECORET   ||  
           setup->servo==TESBIASSERVO || setup->servo==HEATERSERVO  ||
          (setup->servo==SQ2SERVO && setup->stepBIAS==0) )
      {
         setup->biasFlag=BIASNOCHANGE;
      } 
      while( fdbkInx < setup->fdbkNo )
      {  
        *isEnd=0;      
        // if setup->servo ==SQ2SERVO,SQ2LOCK or SQ1SERVO lock SQ1BIASSERVO
        //     save sa_fb value in sa_fb row (for SQ2)
        //     save sq2_fb value in sq2_fb row (for SQ1) 
        // for SQ1BIASSERVO, apply servo algorithm, use fdbk as bias
        sc2daservoAlgrm(con,byte,myInfo,setup,&ptrFlag,fdbkInx,dataBuf,status);
        if ( !StatusOkP(status) )
        {
          ErsRep (0, status,"sc2daservoFunc: sc2daservoAlgrm failed"); 
          return;
        } 
	
	// set FB  for next data no change of Bias now
        fdbkInx++;
        setup->fdbk =setup->minFDBK + fdbkInx*setup->stepFDBK;
        setup->biasFlag=BIASNOCHANGE; 

        if (setup->servo !=BLACKBODY)
        { 
          if(fdbkInx==setup->fdbkNo)
          {
            // fdbk from 0 to fdbkNo-1, so up to here, set fdbkInx=0
            fdbkInx=0; 
            setup->fdbk =setup->minFDBK;
            break;  // jump out the feedback loop
          }
          else
          {
            sc2dasetBiasFdbk(con,myInfo,setup,mceInxpt,dateTime,fdbkInx,status);
            if ( !StatusOkP(status) )
            {
              ErsRep (0, status,"sc2daservoFunc: setBaisFdbk failed"); 
              return;
            }
          }
        }
        else // blackbody delay 
        {
          sleep(setup->slopSelect[0]);
        } 
        return;  //trigging data for next round
      }
      
      fprintf(servoFp,"\n"); 
      if(myInfo->filesvFlag!=0)
        fflush(myInfo->fpData);

      // now allow ssaramp and cablecorrect to adjust
      if (setup->servo==SSARAMP ||setup->servo==CABLECORET )
          sc2dafindcableadjVal(setup,biasInx,status);

      // if it is cable correction, then correct the offset  if
      // setup->slopSelect[9]==1
      if ( (setup->servo==CABLECORET)  &&  setup->slopSelect[9]==1  )
        sc2dacableCorrection(con,myInfo,setup,mceInxpt,dateTime,status);

      biasInx ++;
      if(biasInx < setup->biasNo)
      {
        // set FB/BIAS in SQUID to minFDBK and next bias for next data
	setup->biasFlag=BIASCHANGED;
        setup->bias = setup->minBIAS + biasInx*setup->stepBIAS;
        jitDebug(2,"sc2daservoFunc: ready for next BIAS \n");
        sc2dasetBiasFdbk(con,myInfo,setup,mceInxpt,dateTime,fdbkInx,status);		         
        if(setup->servo==SQ2SERVO)
        {
          ssalckset->doServo=SQ2SERVOSSALOCK;
          sq2ssalckFlag=1;
          jitDebug(2,"sc2daservoFunc: ready for next ramp SSA FB for new lock point \n");
        }
        return;  //trig data for next round       
      }
      else
      {
        sprintf(msg,"end biasInx=<%d> =====",biasInx);
        sc2dalib_msgprintSave(myInfo,"sc2daservoFunc: %s",msg,NO_PRINT,status);

        biasInx=0;            *isEnd=1;
        if (servoFp !=NULL)   fclose(servoFp);

        if (setup->servo==SQ1SERVO  && fluxjmpFp !=NULL)
        {
          jitDebug(16,"sc2daservoFunc: call fclose(fluxjmpFp)\n");
          fclose(fluxjmpFp);
        }

        if (setup->servo==SQ1OPEN || setup->servo==TESTRANSIT || 
            setup->servo==SQ2OPEN || setup->servo==SQ2OPEN4P )
        {
          if( myInfo->fpSq1 !=NULL)
          {
            jitDebug(16,"sc2daservoFunc: call fclose(fpSq1)\n");
            fclose(myInfo->fpSq1);
          }
        }  
        if (myInfo->dataFormat==MCE_BINARY_FORM)
          sc2daservosavebindata2File(myInfo,setup,status);

        break;
      }
    }
  }
}



/**
 * \fn void sc2daservoFunc2(SDSU_CONTEXT *con, char * byte, 
 *  dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *  mcexml_struct *mceInxpt, int *isEnd, char *dataBuf, char *dateTime,
 *  ARRAYSET *ssalckset,  char *ssalckBuf, StatusType *status)
 *
 * \brief function
 *  start to ramp feedback and bias, apply servo algorithm, but ecah channel or row
 * is different depending on biaslckPt[],biasoptPt[] and sq2fdbkopt[]
 * currently only for SQ2OPEN4P
 *  
 *
 * \param con      SDSU_CONTEXT pointer
 * \param byte     current raw frame data pointer
 * \param myInfo  dasInfoStruct_t poiter
 * \param setup    ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param isEnd    end flag     given & return
 * \param dataBuf  char pointer for servo data
 * \param dateTime char pointer for date Time string
 * \param ssalckset     ARRAYSET structure pointer
 * \param ssalckBuf     char pointer for ssa data
 * \param status   StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.

 * in sq2servo, for each sq2bias,ssalckFunc first finds initFB[ch], zfact[ch],
 * and gain[ch] and update them for sq2servo to use
 *
 * the question is the latest ssafb found from ssafb/sq2fb curve 
 * obtained from sq2servo data will be different from ssalock.
 * do we need to do a separate ssalock to update sq2open.txt? 
 *
 */
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
)
{
  static int  fdbkInx,totalStep;
  static int  ptrFlag;       // flag for initialise  pointer 
  
  if (!StatusOkP(status)) return;

  if(*isEnd==1)
  {
    *isEnd=0;
    totalStep=2*SQ2FBSTEPDOWN;
    fdbkInx=ptrFlag=0;
  }
  else
  {
    while(fdbkInx<totalStep )
    {  
      *isEnd=0;
      sc2daservoAlgrm(con,byte,myInfo,setup,&ptrFlag,fdbkInx,dataBuf,status);
      if ( !StatusOkP(status) )
      {
        ErsRep (0, status,"sc2daservoFunc: sc2daservoAlgrm failed"); 
        return;
      }   
      sc2dasetFDBK(con,myInfo,setup,mceInxpt,dateTime,fdbkInx,status);		       
      if ( !StatusOkP(status) )
      {
        ErsRep (0, status,"sc2daservoFunc2: setfdbk2 failed"); 
        return;
      }
      fdbkInx++;
      if(fdbkInx==totalStep)
      {
        fdbkInx=0;
        break;  // jump out the feedback loop
      }
      return;  //trigging data for next round
    }

    fprintf(servoFp,"\n"); 
    if(myInfo->filesvFlag!=0)
      fflush(myInfo->fpData);
  
    *isEnd=1;
    if (servoFp !=NULL)   fclose(servoFp);
    if( myInfo->fpSq1 !=NULL)
    {
      jitDebug(16,"sc2daservoFunc2: call fclose(fpSq1)\n");
      fclose(myInfo->fpSq1);
    }
    if (myInfo->dataFormat==MCE_BINARY_FORM)
      sc2daservosavebindata2File(myInfo,setup,status);
  }
}


/**
 * \fn void sc2daservossalckFunc(SDSU_CONTEXT *con, char * byte, 
 *  dasInfoStruct_t *myInfo, ARRAYSET *setup, ARRAYSET *ssalckset,
 *  mcexml_struct *mceInxpt, char *dataBuf, char *dateTime,
 *  int *initFlag, BIAS_PEAK *peakbiasPtr,int biasInx, 
 *  char *servoData, StatusType *status)
 *
 * \brief function
 *  start to do ssalock servo to update initFB , Zfact and gain
 *  to eleminate the sq2bais change effect
 *
 * \param con       SDSU_CONTEXT pointer
 * \param byte      current raw frame data pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param setup     ARRAYSET structure pointer
 * \param ssalckset ARRAYSET structure pointer
 * \param mceInxpt  mcexml_struct pointer
 * \param dataBuf   char pointer for ssalck data
 * \param dateTime  char pointer for date Time string
 * \param initFlag  int: flag for initialize some para
 * \param peakbiasPtr BIAS_PEAK structure pointer
 * \param biasInx   int the bias step index for main sq2
 * \param servoData char pointer for sq2servo data
 * \param status    StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  int        i;
   int        j;
  char       tmp[CMD_LEN];
  dasCmdInfo_t  biasCmd;
  BIAS_PEAK      *peakperbiasPtr;
  static int  fdbkInx;
  static int  ptrFlag;       // flag for initialise  pointer 
  static int  resetSSA;
  static int  sq2bias;

  if (!StatusOkP(status)) return;
 
  if( *initFlag==1)
  {
    sc2dalib_msgprintSave(myInfo,"sc2daservossalckFunc","",NO_PRINT,status);
    sq2bias = setup->minBIAS + biasInx*setup->stepBIAS;
    fdbkInx=ptrFlag=0;
    resetSSA=0;
    /*****  no need now as servoFunc set sq2bias already 
    **/
  }

  if(resetSSA==0)
  {
    // in ssalck, ssaBias ia already applied, so, don't do bias
    ssalckset->biasFlag=BIASNOCHANGE;
    // here we only do once for each ssafb, not twice
    while(fdbkInx<ssalckset->fdbkNo )
    {  
      *initFlag=0;
      // ssalckset->doServo ==SQ2SERVOSSALOCK
      sc2daservoAlgrm(con,byte,myInfo,ssalckset,&ptrFlag,fdbkInx,dataBuf,status);
      if ( !StatusOkP(status) )
      {
        ErsRep (0, status,"sc2daservossalckFunc: sc2daservoAlgrm failed"); 
        return;
      }   
 
      fdbkInx++;
      if(fdbkInx < ssalckset->fdbkNo)
      {
        ssalckset->fdbk =setup->minFDBK + fdbkInx*setup->stepFDBK;
        // send MCE command to change ssafb 
        sc2dasetBiasFdbk(con,myInfo,ssalckset,mceInxpt,dateTime,fdbkInx,status);		       
        if ( !StatusOkP(status) )
        {
          ErsRep (0, status,"sc2daservossalckFunc: setBaisFdbk failed"); 
          return;
        }
      }
      else
      {
        //sc2dalib_msgprintSave(myInfo,"sc2daservossalckFunc","",0,status);

        // findlock point from ssalckdataBuf
        // tmp save to cmdrepFile
        // ssalckset->saoutlckVal[ch]=data0;
        // ssalckset->initFB[ch]=j*ssalckset->stepFDBK+ssalckset->minFDBK;

        sprintf (tmp, "%sssalck-%d", myInfo->dataFile,biasInx );  
        sc2dafindlockPts(myInfo,ssalckset,dataBuf,myInfo->fpOtheruse,tmp,status);
     
        // ssalck never uses minBIAS, ssalckset->biasNo=1;
        ssalckset->minBIAS=sq2bias;
        sc2dasavelockPts(myInfo,ssalckset,dataBuf,myInfo->fpOtheruse,status);

        peakperbiasPtr=peakbiasPtr + biasInx;
        jitDebug(16,"sq2-ssalock:peakperbiasPtr=%0X\n", 
             ((char*)peakperbiasPtr-servoData));
        jitDebug(16,"sq2-ssalock:peakperbiasPtr->initFB[0] =%0X\n", 
             ((char*)&peakperbiasPtr->initFB[0]-servoData));

        // update initFB , Zfactor and Gain, if no peak (maxValue.peakval[ch][0]=0,
        // then initFB no changed from reading
        for( i=0;i<COL_NUM;i++)
        {
          setup->initFB[i]=ssalckset->initFB[i];
          setup->gain[i]=ssalckset->gain[i];
          setup->zfact[i]=ssalckset->saoutlckVal[i];
#ifndef SDSU
          peakperbiasPtr->initFB[i]=setup->initFB[i];
          peakperbiasPtr->zfact[i]=setup->zfact[i];
          peakperbiasPtr->gain[i]=setup->gain[i];
#else
          peakperbiasPtr->initFB[i]=setup->initFB[i]+10;
          peakperbiasPtr->zfact[i]=setup->zfact[i]+0x41;
          peakperbiasPtr->gain[i]=setup->gain[i]+0x42;
#endif
        }
        // need to set initFB and send GO again for sq2servo to work
        jitDebug(2,"ssalckFunc: End of sqBias(%d), set SSA-fdbk for servoFunc\n",sq2bias);  
        strcpy(biasCmd.mceCmd,"");
        sprintf (biasCmd.mceCmd, "wb %s flux_fb", setup->safbCard);
        for ( j=0; j<COL_NUM; j++ )
        {
          sprintf (tmp, " %d", setup->initFB[j]*setup->colMask[j] );
          strcat(biasCmd.mceCmd,tmp); 
        }
        jitDebug(2,"SSA-fdbk "); 
        sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
        if ( !StatusOkP(status) )
        {
          ErsRep (0, status,
            "sc2daservossalckFunc: sc2dalib_sendCmd failed SSA FDBK (%s)",
           biasCmd.mceCmd); 
         return;
        }
        resetSSA=1;
        break;  // jump out the feedback loop
      }
      break;
    }
  }
  else
  {
    ssalckset->doServo=NONESERVO;
  }
}



/**
 * \fn void sc2dasetBiasFdbk(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,mcexml_struct *mceInxpt, char * dateTime, int fdbkInx,
 *  StatusType *status)
 *
 * \brief function
 *  send out bias, feed back values to MCE
 *
 * \param con    SDSU_CONTEXT pointer
 * \param myInfo dasInfoStruct_t poiter
 * \param setup   ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param dateTime char pointer for date Time string
 * \param fdbkInx  int   index for feed back
 * \param status   StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  if (!StatusOkP(status)) return;

  if (setup->servo==SSARAMP  || setup->servo==SSALOCK    || 
      setup->servo==SSARAMP1 || setup->servo==SQ2BIASING || 
      setup->servo==CABLECAL || setup->servo==CABLECORET
     )
  {
    sc2dasetbiasfdbkSSA(con,myInfo,setup,mceInxpt,dateTime,fdbkInx,status);
  }
  else if (setup->servo==SQ2SERVO || setup->servo==SQ2OPEN ||
           setup->servo==SQ2LOCK  || setup->servo==SQ2OPEN4P
          )
  {
    sc2dasetbiasfdbkSQ2(con,myInfo,setup,mceInxpt,dateTime,fdbkInx,status);
  }
  else if (setup->servo==SQ1SERVO || setup->servo==SQ1OPEN ||
           setup->servo==HEATERSERVO || setup->servo==TESBIASSERVO ||
	   setup->servo ==SQ1BIASSERVO || setup->servo ==SQ1LOCK )
  {
#ifndef SDSU
    sc2dasetbiasfdbkSQ1(con,myInfo,setup,mceInxpt,dateTime,fdbkInx,status);
#endif
  }  
}


/**
 * \fn void sc2daservoresetFdbk(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,mcexml_struct *mceInxpt, char * dateTime,
 *  int *biasFlag, StatusType *status)
 *
 * \brief function: 
 *  send out 0 to feed back values 
 *  not use yet
 *
 * \param con    SDSU_CONTEXT pointer
 * \param myInfo dasInfoStruct_t poiter
 * \param setup   ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param dateTime char pointer for date Time string
 * \param biasFlag int pointer for biasFlag
 * \param status   StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  int     i,j;  
  char   tmp[CMD_LEN];
  dasCmdInfo_t  myCmd;

  if (!StatusOkP(status)) return;

  setup->fdbk=0;
  // SSA feedback 
  strcpy(myCmd.mceCmd,"");
  sprintf (myCmd.mceCmd, "wb %s flux_fb", setup->safbCard);
  for ( j=0; j<COL_NUM; j++ )
  {
    sprintf (tmp, " %d", setup->fdbk );
    strcat(myCmd.mceCmd,tmp); 
  }
  sc2dalib_sendCmd(con,myInfo,&myCmd,mceInxpt,dateTime,status);
  if ( !StatusOkP(status) )
  {
    ErsRep (0, status,
      "sc2dasetFDBK: sc2dalib_sendCmd failed SSA FDBK (%s)",
           myCmd.mceCmd); 
    return;
  }
  // Sq2 feedback 
  strcpy(myCmd.mceCmd,"");
  sprintf (myCmd.mceCmd, "wb %s flux_fb", setup->sq2fbCard);
  for ( j=0; j<COL_NUM; j++ )
  {
    sprintf (tmp, " %d", setup->fdbk );
    strcat(myCmd.mceCmd,tmp); 
  }
  sc2dalib_sendCmd(con,myInfo,&myCmd,mceInxpt,dateTime,status);
  if ( !StatusOkP(status) )
  {
    ErsRep (0, status,
      "sc2dasetFDBK: sc2dalib_sendCmd failed SQ2 FDBK (%s)",
       myCmd.mceCmd); 
    return;
  }
  //SQ1 feedback 
  strcpy(myCmd.mceCmd,"");
  for(i=0;i<setup->totalRC2use;i++)
  {
    sprintf ( myCmd.mceCmd, "wb rc%d fb_const",setup->whichRC[i] ); 
    for ( j=0; j<ONECARDCHS; j++ )
    {
      sprintf ( tmp, " %d", setup->fdbk );
      strcat(myCmd.mceCmd,tmp); 
    }
    sc2dalib_sendCmd(con,myInfo,&myCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,
       "sc2dasetFDBK: sc2dalib_sendCmd failed SQ1 FDBK (%s)",
        myCmd.mceCmd); 
      return;
    }
  }
}


/**
 * \fn void sc2dasetbiasfdbkSSA(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,mcexml_struct *mceInxpt, char * dateTime, int fdbkInx,
 *  StatusType *status)
 *
 * \brief function: 
 *  send out ssa bias, ssa feed back values to MCE 
 *
 * \param con    SDSU_CONTEXT pointer
 * \param myInfo dasInfoStruct_t poiter
 * \param setup   ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param dateTime char pointer for date Time string
 * \param status   StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * reduce cable offset change by factor 2 (bdk 13June2006) 
 * reduce cable offset change by factor CABLE_SCALE (xg 28June2006) 
 * calculate cable offset xg 25sept2006
 *
 */
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
)
{
  int    j,i,read8, a;
  char   tmp[CMD_LEN];
  dasCmdInfo_t  biasCmd;

  if (!StatusOkP(status)) return;

  strcpy(biasCmd.mceCmd,"");

  // if servo==SSALOCK, *biasFlag !=BIASCHNAGED)
  if(setup->biasFlag==BIASCHANGED)
  {
    //SSA bias  COL_NUM (32)
    if (setup->servo ==SSARAMP || setup->servo==CABLECAL  ||
        setup->servo ==SSARAMP1 )
    {       
      if (setup->servo==CABLECAL)
        setup->bias=0;

      for(i=0; i<setup->totalRC2use;i++)
      {
        read8=(setup->whichRC[i]-1)*ONECARDCHS;
        sprintf (biasCmd.mceCmd, "wb rc%d sa_bias",setup->whichRC[i] ); 
        for ( j=read8; j<(read8+ONECARDCHS); j++ )
        {
          sprintf ( tmp, " %d", setup->bias *setup->colMask[j]);
          strcat(biasCmd.mceCmd,tmp); 
        }
        sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
        if ( !StatusOkP(status) )
        {
          ErsRep (0, status,
            "sc2dasetbiasfdbkSSA: sc2dalib_sendCmd failed SSA BIAS (%s)",
             biasCmd.mceCmd); 
          return;
        }
        //SSA offset due to cable impendance 
        sprintf (biasCmd.mceCmd, "wb rc%d offset ",setup->whichRC[i] ); 
        for ( j=read8; j<(read8+ONECARDCHS); j++ )
        {
          if(setup->servo !=SSARAMP)
            setup->cableScale[j]=0;

          if ( setup->doServo ==3)
            a=0;
          else
           // cableScale <0 normally, cableadjVal >0 see sc2dafindcableadjVal()
           // see servoInit
           // setup->cableadjVal[i]=setup->cableadjInit[i]; or=0;

           a=setup->cableOffset[j]- (int)(setup->bias*setup->cableScale[j])
             +setup->cableadjVal[j];
 
          sprintf ( tmp, " %d", a*setup->colMask[j]);
          strcat(biasCmd.mceCmd,tmp); 
        }
        sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
        if ( !StatusOkP(status) )
        {
           ErsRep (0, status,
            "sc2dasetbiasfdbkSSA: sc2dalib_sendCmd failed SSA OFFSET (%s)",
             biasCmd.mceCmd); 
           return;
        }
      }
    }
    else if (setup->servo ==SQ2BIASING)
    {
      // SQ2 bias 
      sprintf (biasCmd.mceCmd, "wb %s flux_fb",setup->sq2biasCard ); 
      for ( j=0; j<COL_NUM; j++ )
      {
        sprintf ( tmp, " %d", setup->bias*setup->colMask[j] );
        strcat(biasCmd.mceCmd,tmp); 
      }
    
      sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
      if ( !StatusOkP(status) )
      {
        ErsRep (0, status,
          "sc2dasetbiasfdbkSSA: sc2dalib_sendCmd failed SQ2 BIAS (%s)",
           biasCmd.mceCmd); 
        return;
      }
    }
    // when ssabias and sq2bias change, we need to give time for circuit 
    // to settle down   in sc2da_par.h
    sleep(WAIT_TIME_AFTER_BIAS);;
  }
  // SSA feedback 
  sprintf (biasCmd.mceCmd, "wb %s flux_fb", setup->safbCard);
  for ( j=0; j<COL_NUM; j++ )
  {
     sprintf (tmp, " %d", setup->fdbk*setup->colMask[j] );
     strcat(biasCmd.mceCmd,tmp); 
  }
  jitDebug(2,"%3d-SSA-fdbk ",fdbkInx);
  sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
  if ( !StatusOkP(status) )
  {
    ErsRep (0, status,
      "sc2dasetbiasfdbkSSA: sc2dalib_sendCmd failed SSA FDBK (%s)",
       biasCmd.mceCmd); 
    return;
  }
}



/**
 * \fn void sc2dacableCorrection(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,mcexml_struct *mceInxpt, char * dateTime,
 *  StatusType *status)
 *
 * \brief function: 
 *  appply the cableadjVal to the cable offset to compensate any 
 *  cable drift
 *
 * \param con    SDSU_CONTEXT pointer
 * \param myInfo dasInfoStruct_t poiter
 * \param setup   ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param dateTime char pointer for date Time string
 * \param status   StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  int    j,i,read8, a[COL_NUM];
  char   tmp[FILE_LEN];
  char   tmpfile[FILE_LEN];
  FILE   *fp;

  dasCmdInfo_t  biasCmd;

  if (!StatusOkP(status)) return;

  strcpy(biasCmd.mceCmd,"");

  fprintf(myInfo->fpLog,"cableCorection \n");
  sprintf(tmpfile,"%s-log",myInfo->dataFile);
  if((fp = fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dacableCorrection: Error- failed to open file %s", tmpfile); 
      return;
    }

  fprintf(fp,"cableOffset  ssabiasklckPt cableScale    cableadjVal cableadjInit \n");   
  fprintf(fp,"rc* offset=Offset[j]- (int)(ssabiaslckPt[j]*Scale[j]) +adjVal[j]\n");

  for(i=0; i<setup->totalRC2use;i++)
  {
    read8=(setup->whichRC[i]-1)*ONECARDCHS;
    //SSA offset due to cable drift 
    sprintf (biasCmd.mceCmd, "wb rc%d offset ",setup->whichRC[i] ); 

    for ( j=read8; j<(read8+ONECARDCHS); j++ )
    {
      // cableScale <0 normally, cableadjVal >0 see sc2dafindcableadjVal()
      a[j]=setup->cableOffset[j]
        - (int)(setup->ssabiaslckPt[j]*setup->cableScale[j])
        +setup->cableadjVal[j];

      sprintf ( tmp, " %d", a[j]*setup->colMask[j]);
      strcat(biasCmd.mceCmd,tmp); 
    }
    fprintf(fp,"rc%d\n",setup->whichRC[i]);

    for ( j=read8; j<(read8+ONECARDCHS); j++ )
    {
      fprintf(fp,"%7d %11d %12.2f %13d %13d\n",setup->cableOffset[j],
        setup->ssabiaslckPt[j], setup->cableScale[j],setup->cableadjVal[j],
        setup->cableadjInit[j]);
    }
 
    sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,
           "sc2dacableCorrection: sc2dalib_sendCmd failed SSA OFFSET (%s)",
            biasCmd.mceCmd); 
       return;
    }
    fprintf(fp,"%s \n",biasCmd.mceCmd);
  }
  fclose(fp);

  // update ssabiaslock batch file
  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ),SSABIASLCK );
  if((fp = fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dacableCorrection 2: Error- failed to open file %s", tmpfile); 
      return;
    }

  fprintf(fp,"# this mce setup file is extracted after cable correction \n");
  fprintf(fp,"# for non-modulation ones, set bias to NONMOD_BIAS (%d)\n",
              NONMOD_BIAS);
  fprintf(fp,"# offset =cableOffset[ch]-(int)(ssabiaslckPt[ch]*cableScale[])\n");
  fprintf(fp,"#         +cableadjVal[j]\n");
  
  for (i=0;i<setup->totalRC2use;i++)
  {
    fprintf (fp,"\nwb rc%d sa_bias",setup->whichRC[i] ); 
    for (j=i*ONECARDCHS; j<(i+1)*ONECARDCHS;j++)
    {
       fprintf(fp,"%7d ", setup->ssabiaslckPt[j]);
    }
    fprintf (fp,"\nwb rc%d  offset ",setup->whichRC[i] ); 
    for (j=i*ONECARDCHS; j<(i+1)*ONECARDCHS;j++)
    {
      fprintf( fp,"%7d ", a[j]*setup->colMask[j]);
    }
  }
  fprintf (fp,"\n");
  fclose(fp);
  fflush(myInfo->fpLog);
  }



/**
 * \fn void sc2dasetbiasfdbkSQ2(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,mcexml_struct *mceInxpt, char * dateTime,int fdbkInx,
 *  StatusType *status)
 *
 * \brief function
 *  send out ssafb from initFB, sq2bias, sq2feedback values to MCE
 *   servo=SQ2SERVO,SQ2LOCK. only sq2 feed back values to MCE 
 *   servo=SQ2OPEN, SQ2OPEN4P
 *
 * \param con    SDSU_CONTEXT pointer
 * \param myInfo dasInfoStruct_t poiter
 * \param setup   ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param dateTime char pointer for date Time string
 * \param fdbkInx   int feedback index
 * \param status   StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  int    j, wait, ssafb;
  char   tmp[CMD_LEN];
  dasCmdInfo_t  biasCmd;

  strcpy(biasCmd.mceCmd,"");
    
  if (!StatusOkP(status)) return;

  if( setup->servo==SQ2SERVO || setup->servo==SQ2LOCK )
  {    
    // SSA feedback 
    //if(setup->biasFlag==BIASCHANGED )
    //  jitDebug(2,"initial SSA-fdbk=0 "); 
    //else
    //  jitDebug(2,"update SSA-fdbk ");

    sprintf (biasCmd.mceCmd, "wb %s flux_fb", setup->safbCard);
    for ( j=0; j<COL_NUM; j++ )
    {
      //if(setup->biasFlag==BIASCHANGED )
      //  ssafb=0;
      //else
      ssafb=setup->initFB[j]*setup->colMask[j];

      sprintf (tmp, " %d", ssafb );
      strcat(biasCmd.mceCmd,tmp); 
    }
    sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0,status,"sc2dasetBaisFdbk: sc2dalib_sendCmd failed SSA FDBK (%s)",
         biasCmd.mceCmd); 
      return;
    }
    if(setup->biasFlag==BIASCHANGED)
      fprintf(myInfo->fpLog,"%s,\n",biasCmd.mceCmd);   

  }

  if(setup->biasFlag==BIASCHANGED)
  {
    // SQ2 bias 
    sprintf (biasCmd.mceCmd, "wb %s flux_fb",setup->sq2biasCard ); 
    for ( j=0; j<COL_NUM; j++ )
    {
      sprintf ( tmp, " %d", setup->bias*setup->colMask[j] );
      strcat(biasCmd.mceCmd,tmp); 
    }
    jitDebug(2,"SQ2-bias ");
    sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2dasetBaisFdbk: sc2dalib_sendCmd failed SQ2 BIAS (%s)",
         biasCmd.mceCmd); 
      return;
    }
    fprintf(myInfo->fpLog,"%s,\n",biasCmd.mceCmd);   
    // when ssabias and sq2bias change, we need to give time for circuit 
    // to settle down
    sleep(WAIT_TIME_AFTER_BIAS);;
  }

  //SQ2 feedback 
  sprintf (biasCmd.mceCmd, "wb %s flux_fb",setup->sq2fbCard ); 
  for ( j=0; j<COL_NUM; j++ )
  {
    sprintf ( tmp, " %d", setup->fdbk*setup->colMask[j] );
    strcat(biasCmd.mceCmd,tmp); 
  }
  jitDebug(2,"%3d-SQ2-fdbk ",fdbkInx);
  sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
  if ( !StatusOkP(status) )
  {
    ErsRep (0, status,"sc2dasetBaisFdbk: sc2dalib_sendCmd failed SQ2 FDBK (%s)",
       biasCmd.mceCmd); 
    return;
  }
  if(setup->biasFlag==BIASCHANGED)
    fprintf(myInfo->fpLog,"%s,\n",biasCmd.mceCmd);   
  
  if (setup->servo ==SQ2OPEN4P)
  {
    //add delay here  11/12/08 xg
    // allow  fraction of second if  slopSelect[7] > 0
    if( setup->slopSelect[7] > 0 )
    {
      wait=1000000*setup->slopSelect[6]/setup->slopSelect[7]; 
      usleep(wait);
    }
  }
}


/**
 * \fn void sc2dasetbiasfdbkSQ1(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,mcexml_struct *mceInxpt, char * dateTime, int fdbkInx,
 *  StatusType *status)
 *
 * \brief function
 *  send out sq2fb from initFB, sq1 bias, feed back values to MCE 
 *  servo=SQ1SERVO,  only sq1 feed back values to MCE servo=SQ1OPEN
 *
 * \param con    SDSU_CONTEXT pointer
 * \param myInfo dasInfoStruct_t poiter
 * \param setup   ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param dateTime char pointer for date Time string
 * \param status   StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  int    j,i,read8,fdbk,flag;
  char   tmp[CMD_LEN];
  dasCmdInfo_t  biasCmd;

  strcpy(biasCmd.mceCmd,"");
    
  if (!StatusOkP(status)) return;

  // for sq1biasservo, change sq2fb
  if( setup->servo != SQ1OPEN  && setup->servo != HEATERSERVO )
  {
    // Sq2 feedback   apply *setup->colMask[j]
    sprintf (biasCmd.mceCmd, "wb %s flux_fb", setup->sq2fbCard);
    for ( j=0; j<COL_NUM; j++ )
    {
      // here we need to implement flux jump for MCE, but the sq2fb 
      // in software is kept unchanged. setup->fluxPeriod[j]=0 if NONMODULATION     
      // we also need to consider if fdkb<0, MCE will not take that
      fdbk=setup->initFB[j]*setup->colMask[j] ;
      flag=0;
      if (setup->fluxPeriod[j] !=0)
      {
        if (fdbk>=0) 
        { 
          // the fluxperiod seem ok now, it gives right answer
          while( fdbk > setup->fluxPeriod[j])
          {
           fdbk -=setup->fluxPeriod[j];
           flag=1;
          }
        }
        else
        {
          while( fdbk < 0 )
          {
            fdbk +=setup->fluxPeriod[j];
            flag=1;
          }
        }
      }
      sprintf (tmp, " %d", fdbk);
      strcat(biasCmd.mceCmd,tmp);
      if( flag ==1)
      {  
      }
    }
    jitDebug(2,"update SQ2 FB ");
    sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,"sc2dasetbiasfdbkSQ1: sc2dalib_sendCmd failed SQ2 FDBK (%s)",
        biasCmd.mceCmd); 
      return;
    }

    if(setup->biasFlag==BIASCHANGED)
      fprintf(myInfo->fpLog,"%s,\n",biasCmd.mceCmd);   

  }

  // now for bias or sq1fb
  if( setup->servo==HEATERSERVO || setup->servo==TESBIASSERVO )
  {
    // ramp heater Val
    if (setup->servo==HEATERSERVO)
      sprintf ( biasCmd.mceCmd, "wb bc1 bias %d",setup->fdbk ); 
    else // ramp TES bias Val 
      sprintf ( biasCmd.mceCmd, "wb bc2 bias %d",setup->fdbk );
 
    sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,
       "sc2dasetbiasfdbkSQ1: sc2dalib_sendCmd failed for (%s)",
          biasCmd.mceCmd); 
        return;
    }
  }
  else if ( setup->servo ==SQ1BIASSERVO )
  {
    // ramp SQ1 bias and take data for sqibiasservo
    // we have swaped fdbk with bias,  so here we use fdbk 
    sprintf ( biasCmd.mceCmd, "wb %s on_bias",setup->sq1biasCard ); 
    // zero value to masked out rows
    for ( j=0; j<ROW_NUM; j++ )
    {
      sprintf (tmp, " %d", setup->fdbk*setup->rowMask[j]);
      strcat(biasCmd.mceCmd,tmp); 
    }
    sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,
        "sc2dasetbiasfdbkSQ1: sc2dalib_sendCmd failed SQ1 BIAS (%s)",
         biasCmd.mceCmd); 
      return;
    }
  }
  else
  {
    //SQ1 feedback 
    for(i=0;i<setup->totalRC2use;i++)
    {
      sprintf ( biasCmd.mceCmd, "wb rc%d fb_const",setup->whichRC[i] ); 
      read8=(setup->whichRC[i]-1)*ONECARDCHS;
      for(j=read8; j<(read8+ONECARDCHS); j++ )
      {
        if (setup->colMask[j]==0)
          sprintf ( tmp, " -8192 ");
        else 
          sprintf ( tmp, " %d", setup->fdbk*setup->colMask[j]);
        strcat(biasCmd.mceCmd,tmp); 
      }
      jitDebug(2,"%3d-SQ1-fdbk ",fdbkInx);
      sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
      if ( !StatusOkP(status) )
      {
        ErsRep (0, status,"sc2dasetbiasfdbkSQ1: sc2dalib_sendCmd failed SQ1 FDBK (%s)",
          biasCmd.mceCmd); 
        return;
      } 
      if(setup->biasFlag==BIASCHANGED)
        fprintf(myInfo->fpLog,"%s,\n",biasCmd.mceCmd);   
    }
    // SQ1 bias 
    if(setup->biasFlag==BIASCHANGED)
    {
      sprintf ( biasCmd.mceCmd, "wb %s on_bias",setup->sq1biasCard ); 
      for ( j=0; j<ROW_NUM; j++ )
      {
        sprintf (tmp, " %d", setup->bias*setup->rowMask[j]);
        strcat(biasCmd.mceCmd,tmp); 
      }
      jitDebug(2,"SQ1-bias");
      sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
      if ( !StatusOkP(status) )
      {
        ErsRep(0,status,"sc2dasetbiasfdbkSQ1: sc2dalib_sendCmd failed SQ1 BIAS (%s)",
           biasCmd.mceCmd); 
        return;
      }
      fprintf(myInfo->fpLog,"%s,\n",biasCmd.mceCmd);   
    }
  }
}


/**
 * \fn void sc2dasetFDBK(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,mcexml_struct *mceInxpt, char * dateTime, int Step,
 *  StatusType *status)
 *
 * \brief function
 *  step up sq2fb from initFBK[col] by step only for SQ2OPEN4P
 *
 * \param con    SDSU_CONTEXT pointer
 * \param myInfo dasInfoStruct_t poiter
 * \param setup   ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param dateTime char pointer for date Time string
 * \param step     int
 * \param status   StatusType pointer.  given and return
 *
 */
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
)
{
  int    j;
  char   tmp[CMD_LEN];
  dasCmdInfo_t  biasCmd;

  if (!StatusOkP(status)) return;

  strcpy(biasCmd.mceCmd,"");

  //SQ2 feedback 
  sprintf (biasCmd.mceCmd, "wb %s flux_fb",setup->sq2fbCard ); 

  for ( j=0; j<COL_NUM; j++ )
  {
    setup->fdbk=setup->initFB[j]*setup->colMask[j];
    if ( setup->fdbk > 0 )
      setup->fdbk  += setup->initFB[j] + setup->stepFDBK*step;
     sprintf ( tmp, " %d", setup->fdbk );
     strcat(biasCmd.mceCmd,tmp); 
  }
  sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
  if ( !StatusOkP(status) )
  {
    ErsRep (0, status,"sc2dasetFDBK: sc2dalib_sendCmd failed SQ2 FDBK (%s)",
       biasCmd.mceCmd); 
    return;
  }
}


/**
 * \fn void sc2dasetstepdownFB(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,mcexml_struct *mceInxpt, char * dateTime, int manySteps,
 *  StatusType *status)
 *
 * \brief function
 *  step down sq2fb from sq2optPt[col] by manySteps only for SQ2OPEN4P
 *
 * \param con    SDSU_CONTEXT pointer
 * \param myInfo dasInfoStruct_t poiter
 * \param setup   ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param dateTime char pointer for date Time string
 * \param  manySteps int
 * \param status   StatusType pointer.  given and return
 *
 *
 */
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
)
{
  int    j, down;
  char   tmp[CMD_LEN];
  dasCmdInfo_t  biasCmd;

  if (!StatusOkP(status)) return;

  strcpy(biasCmd.mceCmd,"");

  //SQ2 feedback 
  sprintf (biasCmd.mceCmd, "wb %s flux_fb",setup->sq2fbCard ); 

  for (down=0; down<manySteps; down++) 
  {
    for ( j=0; j<COL_NUM; j++ )
    {
      setup->fdbk=setup->sq2fdbkOpt[j]*setup->colMask[j];
      setup->fdbk  += setup->sq2fdbkOpt[j]- setup->stepFDBK*down;
      if ( setup->fdbk < 0 )
        setup->fdbk=0;
      sprintf ( tmp, " %d", setup->fdbk );
      strcat(biasCmd.mceCmd,tmp); 
      // temps stored in for last servoFunc2 to use in setfdbk
      setup->initFB[j]=setup->fdbk; 
    }
    sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,"sc2dasetstepdownFB: sc2dalib_sendCmd failed SQ2 FDBK (%s)",
         biasCmd.mceCmd); 
      return;
    }
  }
}



/**
 * \fn void sc2daservosaveData(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,char *byte, int *data, int fdbkInx,StatusType *status)
 *
 * \brief function: 
 *  save servo data into datafile and stripchart data file 
 *
 *   setup->servo==SSARAMP,or SSALOCK,or SQ2BIASING save setup->fdbk
 *   setup->servo==SQ1SERVO or SQ2SERVO             save setup->initFB[i]
 *   setup->servo==SQ1OPEN   or SQ2OPEN              save data[i]
 *
 * for channel data
 *   save initFB[] if setup->servo==SQ1SERVO lock or SQ2SERVO lock, 
 *   the rest servo, always save data[] 
 *
 * \param  con       SDSU_CONTEXT pointer
 * \param  myInfo   dasInfoStruct_t poiter
 * \param  setup     ARRAYSET structure pointer
 * \param  byte      char pointer for frame data array
 * \param  data      int pointer for row data array
 * \param fdbkInx int for fedbk Index
 * \param status     StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  int    j, fdbk;
  int    fdbksubInx;
   
  if (!StatusOkP(status)) return;

  // Mike asks when doServo<0, save whole frame
  //myInfo->bufSize ( in byte) =dramamsg.bufsize in sc2da.c;
  if( setup->doServo <0 )
  {
    //save in datafile-binary
    if (setup->servo==SQ1OPEN || setup->servo==SQ2OPEN || setup->servo==SQ2OPEN4P)   
   	fwrite(byte,1,myInfo->bufSize,myInfo->fpSq1);
    else                         // depending on data-format 
      sc2dalib_saveframeData(con,myInfo,myInfo->fpData,byte,myInfo->bufSize/4,status);
  }

  if ( myInfo->dataFormat==MCE_BINARY_FORM)
  {
      if ( setup->doServo >=0 )
     sc2daservosavebinData(myInfo,setup,data,fdbkInx,status);
  }
  else
  {
    // for keeping idl script
    for ( fdbksubInx=0; fdbksubInx <2; fdbksubInx ++)
    {
      // allow sq1open to have a row data still be albe to use idl script
      if ( setup->doServo >=0 || setup->servo==SQ1OPEN || 
           setup->servo==SQ2OPEN  || setup->servo==SQ2OPEN4P 
         )
      {
        fprintf(myInfo->fpData, "data read %d: ",fdbksubInx); 

        for ( j=0; j<COL_NUM; j++ )
          fprintf ( myInfo->fpData, "%11d ",data[j] );
        fprintf(myInfo->fpData, "\n");

        if (setup->servo==SQ1SERVO || setup->servo==SQ1BIASSERVO ||
             setup->servo ==SQ1LOCK 
           )
          fprintf(myInfo->fpData, "s2_fb set %d: ",fdbksubInx);
        else if (setup->servo==SQ2SERVO || setup->servo ==SQ2LOCK)
          fprintf(myInfo->fpData, "sa_fb set %d: ",fdbksubInx);
        else if (setup->servo==SQ1OPEN || setup->servo==SQ2OPEN)
          fprintf(myInfo->fpData, "sa_out    %d: ",fdbksubInx);
        else if (setup->servo==HEATERSERVO)
          fprintf(myInfo->fpData, "heater    %d: ",fdbksubInx);
        else if (setup->servo==TESBIASSERVO)
          fprintf(myInfo->fpData, "tesbias   %d: ",fdbksubInx);
        else 
          // setup->servo==SSARAMP,or SSALOCK,or SQ2BIASING  
          // SSARAMP1, CABLECAL 
          fprintf(myInfo->fpData, "sa_fb set %d: ",fdbksubInx);
  
        for ( j=0; j<COL_NUM; j++ )
        { 
          if( setup->servo ==SQ2SERVO     || setup->servo ==SQ2LOCK ||
              setup->servo ==SQ1SERVO     || setup->servo ==SQ1LOCK ||
              setup->servo ==TESBIASSERVO || setup->servo ==SQ1BIASSERVO 
            )
          {
            // initFB[] =sa_fb for SQ2, =sq2_fb for SQ1
            fprintf(myInfo->fpData, "%11d ",setup->initFB[j]);
          }
          else if (setup->servo==SQ1OPEN   || setup->servo==SQ2OPEN ||
                   setup->servo==SQ2OPEN4P || setup->servo ==HEATERSERVO 
                   )
            fprintf(myInfo->fpData, "%11d ",data[j]);
          else
            // setup->servo==SSARAMP,or SSALOCK,or SQ2BIASING
            // SSARAMP1, CABLECAL 
            fprintf(myInfo->fpData, "%11d ",setup->fdbk);
        }  
        fprintf(myInfo->fpData, "\n");
      }
    }
    fflush (myInfo->fpData);    
  }
  // ************* for xxxservodata.txt *********************
  // ********************************************************
  // so that we can use EXCEL to plot this file
  // change to use TAB for LabView too X.Gao 20090529
  if ( setup->servo==SQ2OPEN4P && setup->slopSelect[7]==100)
  {
    for ( j=0; j<COL_NUM; j++ )
    {
      fdbk=setup->initFB[j]*setup->colMask[j];
      if ( fdbk > 0 )
        fdbk  += setup->initFB[j] + setup->stepFDBK*fdbkInx;

      fprintf(servoFp,"%d\t", fdbk);
      fprintf(servoFp,"%d\t", data[j]);
    } 
  }
  else
  {
    fprintf(servoFp, "%d\t", setup->bias);
    fprintf(servoFp, "%d\t ", setup->fdbk);
    
    for ( j=0; j<COL_NUM; j++ )
    {
      if( setup->servo ==SQ2SERVO     || setup->servo ==SQ2LOCK ||
          setup->servo ==SQ1SERVO     || setup->servo ==SQ1LOCK ||
          setup->servo ==TESBIASSERVO || setup->servo ==SQ1BIASSERVO 
        )
        fprintf(servoFp,"%d\t",setup->initFB[j]);
      else 
        fprintf(servoFp,"%d\t",data[j]);
    }
  }
  fprintf(servoFp,"\n");
  fflush (servoFp);
 
// ***************** for stripchart ***********************
// ********************************************************
// the data format for SSA  is
// col0   COL1   COL2   COL3   COL4   COL5   COL6   COL7   COL8   COL9 
// 1     fdbk0  data0  data1  data2  data3  data4  data5  data6  data7  
// 2     fdbk1  data0  data1  data2  data3  data4  data5  data6  data7  
//   ...................
//  
// the data format for Sq2(Sq1) is
// col0 COL1   COL2   COL3   COL4   COL5   COL6   COL7   COL8   COL9    COL10  COL11
// 1    fdbk0  data0  data1  data2  data3  data4  data5  data6  data7   fb0    fb1  
// 2    fdbk1  data0  data1  data2  data3  data4  data5  data6  data7   fb0    fb1 
//   ...................
// only do black body, add others if we want later
  if(setup->servo ==BLACKBODY)
  {
    fprintf(myInfo->fpStrchart, "%7d ",stripchFlag);
    fprintf(myInfo->fpStrchart, "%7d ",setup->fdbk);
    for ( j=0; j<COL_NUM; j++ )
    {
      if (setup->colMask[j] !=0) 
        fprintf ( myInfo->fpStrchart, "%11d ",data[j] );
    }	
    fprintf(myInfo->fpStrchart, "\n");
    fflush (myInfo->fpStrchart);    
    stripchFlag++;
  }
}


/**
 * \fn void sc2daservosavebinData(dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,int *data, int fdbkInx,StatusType *status)
 *
 * \brief function: 
 *  save servo data in BINARY into memory 
 *
 *   setup->servo==SSARAMP,or SSALOCK,or SQ2BIASING save setup->fdbk
 *   setup->servo==SQ1SERVO lock or SQ2SERVO lock   save setup->initFB[i]
 *   setup->servo==SQ1OPEN   or SQ2OPEN              save data[i]
 *
 * for channel data
 *   save initFB[] if setup->servo==SQ1SERVO or SQ2SERVO, 
 *   the rest servo, always save data[] 
 *
 * \param  myInfo   dasInfoStruct_t poiter
 * \param  setup     ARRAYSET structure pointer
 * \param  data      int pointer for row data array
 * \param fdbkInx int for fedbk Index
 * \param status     StatusType pointer.  given and return
 *
 */
/*+ sc2daservosavebinData
*/
void sc2daservosavebinData
(
dasInfoStruct_t  *myInfo,
ARRAYSET         *setup,
int              *data,
int              fdbkInx,
StatusType       *status
)
{
  int               j,fdbksubInx;
  SERVOBIN_DATAHEAD *binHead;
  int               *servodataPtr;
  static int        i;  //this is location index

  if (!StatusOkP(status)) return;

  binHead =(SERVOBIN_DATAHEAD *)myInfo->servobindataPtr;
  servodataPtr=(int *)(binHead +1);

  // only initialise at the beginning
  if (fdbkInx==0)    i=0;

  // allow sq1open to have a row data still be albe to use idl script
  // for keeping idl script
  if ( setup->doServo >=0 || setup->servo==SQ1OPEN)
  {
    for ( fdbksubInx=0; fdbksubInx <2; fdbksubInx ++)
    {
      for ( j=0; j<COL_NUM; j++ )
      {
        servodataPtr[i]=data[j];
        i++;
      }
      for ( j=0; j<COL_NUM; j++ )
      { 
        if( setup->servo ==SQ2SERVO     || setup->servo ==SQ2LOCK ||
            setup->servo ==SQ1SERVO     || setup->servo ==SQ1LOCK ||
            setup->servo ==TESBIASSERVO || setup->servo ==SQ1BIASSERVO 
          )
        {
          // initFB[] =sa_fb for SQ2, =sq2_fb for SQ1
          servodataPtr[i]=setup->initFB[j];
        }
        else if (setup->servo==SQ1OPEN   || setup->servo==SQ2OPEN ||
                 setup->servo==SQ2OPEN4P ||setup->servo ==HEATERSERVO 
                )
        {
          servodataPtr[i]=data[j];
        }
        else
        {  // setup->servo==SSARAMP,or SSALOCK,or SQ2BIASING
           // SSARAMP1, CABLECAL 
          servodataPtr[i]=setup->fdbk;
        }
        i++;
      }
    }  
  }
}


/**
 * \fn void sc2daservosavebindata2File(dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup,StatusType *status)
 *
 * \brief function: 
 *  save servo data from memory into datafile 
 *
 *
 * \param  myInfo   dasInfoStruct_t poiter
 * \param  setup     ARRAYSET structure pointer
 * \param status     StatusType pointer.  given and return
 *
 */
/*+ sc2daservosavebindata2File
*/
void sc2daservosavebindata2File
(
dasInfoStruct_t  *myInfo,
ARRAYSET         *setup,
StatusType       *status
)
{
  int  servobindataSize;
  
  if (!StatusOkP(status)) return;

  servobindataSize=sizeof(SERVOBIN_DATAHEAD) + 
         sizeof(int)*setup->fdbkNo*COL_NUM*4*setup->biasNo;

  fwrite(myInfo->servobindataPtr, 1, servobindataSize,myInfo->fpData);
  free( myInfo->servobindataPtr);
 }


/**
 * \fn void sc2datransitFunc(SDSU_CONTEXT *con, char * byte, 
 *  dasInfoStruct_t *myInfo, ARRAYSET *setup, mcexml_struct *mceInxpt,
 *  int *isEnd, char *dataBuf, char *dateTime,StatusType *status)
 *
 * \brief function
 *  generate a small triangle modulation for each heaterLevel 
 *
 * \param con      SDSU_CONTEXT pointer
 * \param byte     current raw frame data pointer
 * \param myInfo  dasInfoStruct_t poiter
 * \param setup    ARRAYSET structure pointer
 * \param mceInxpt mcexml_struct pointer
 * \param isEnd    end flag     given & return
 * \param dataBuf  char pointer for servo data
 * \param dateTime char pointer for date Time string
 * \param status   StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * in sc2dafindTransit(&dasInfo,&setup,dataBuf,status);
 *  save the servo data so that findlockpoints can use it.        
 */
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
)
{
  static int  heatModulate;
  static char msg[FILE_LEN];
  static int  biasInx,heaterInx;
  static int  ptrFlag,waveFlag;      // flag for initialise  pointer 
  static int  heatStart;

  if (!StatusOkP(status)) return;

  if(*isEnd==1)
  {
    *isEnd=0;
    ptrFlag=waveFlag=0;
    biasInx=0;
    heatStart=setup->minFDBK;
    setup->biasFlag=BIASNOCHANGE;
    setup->totalPixel=ROW_NUM*COL_NUM;
    heatModulate=setup->fdbkNo*4+1; 
    fprintf(myInfo->fpLog,"\n<%s>: transitFunc\n", dateTime);
  }

  while (biasInx<setup->biasNo)
  {
    // use biasNo to stepdown heater
    // for each heatLevel, there are total of (fdbk*4+1) data
    // to cover both end of triangle points: data[0~fdbk*4]
    if( waveFlag==0)
    {
      sprintf(msg,"heaterStepdown=<%d>.....",biasInx);
      sc2dalib_msgprintSave(myInfo,"sc2datransitFunc: %s",msg,USE_PRINT,status);
      // 1 use tri-angle get heater modulation. setup->minFDBK
      sc2dalibsetup_heaterGen(setup,1,status);
      if (biasInx ==0)
        heaterInx=1;
      else
        heaterInx=0;
    }
    while(heaterInx <heatModulate )
    {
      // we already have a GO in sc2da_servo with heater[0]
      waveFlag=1;  
      // save row*col pixels into dataBuf, also as full frame into fpSq1
      //"%s-binary",myInfo->dataFile
      sc2dawavegetData(byte,myInfo,setup,&ptrFlag,dataBuf,status);
      if ( !StatusOkP(status) )
      {
        ErsRep (0, status,"sc2datransitFunc: sc2dawavegetData failed"); 
        return;
      }
     
      setup->fdbk=setup->heater[heaterInx];
      sc2datransitTES(con,myInfo,setup,mceInxpt,dateTime,status);		       
      if ( !StatusOkP(status) )
      {
        ErsRep (0, status,"sc2datransitFunc: sc2datransitTES failed"); 
        return;
      }
      heaterInx ++;
    
      if (heaterInx == heatModulate)
        break;
       return;  //trigging data for next round
    } 
    waveFlag=0;
    biasInx ++;
    setup->minFDBK -=setup->stepBIAS;
    return;
  }  
  // pick up the last one
  //  (biasInx==setup->biasNo)
  sc2dawavegetData(byte,myInfo,setup,&ptrFlag,dataBuf,status);
  if ( !StatusOkP(status) )
  {
    ErsRep (0, status,"sc2datransitFunc: sc2dawavegetData failed"); 
     return;
  }   
  // a section of the triangle wave is completed
  fprintf(myInfo->fpLog,"%03d: %d \n",heatModulate*biasInx,setup->fdbk);
                                                     
  setup->minFDBK=heatStart;
  biasInx=0;
  *isEnd=1;
  if (servoFp !=NULL);
    fclose(servoFp);

  if ( myInfo->fpSq1 !=NULL)  //$datafile-binary
  {
    jitDebug(16,"sc2datransitFunc: call fclose(fpSq1)\n");
    fclose(myInfo->fpSq1);
  }
}


/**
 * \fn void sc2dawavegetData(char *byte, dasInfoStruct_t *myInfo, 
 *  ARRAYSET *setup, int *flagPtr, char *databuf,StatusType *status)
 *
 * \brief function: 
 *  doServo<0    save the whole frame, also save the selected row data
 * 
 * \param  byte   pointer for current raw frame data buffer
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  flagPtr int pointer for   
 * \param  databuf  char pointer for servo data (save CHs for later process)
 * \param  status  StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * the data is always placed in RC1RC3RC4 and zeros are inserted for 
 * missd RCs. see dhtask.  data structure 
 *  RC1=[0: 7]; RC2=[8:15]  RC3=[16:23]; RC4=[24:31]   
 *
 */
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
)
{
  int         *word, data1;
  int         row, ch,j;
  int         data[ROW_NUM*COL_NUM];
  SERVO_DATAHEAD *rampinfoPtr;
  // keep the pointer to memory holding servo data
  static int  *storePtr1;     

     
  if (!StatusOkP(status)) return;

  // initial data[] setup->totalPixel=ROW_NUM*COL_NUM set in transitFunc
  for (j=0; j<setup->totalPixel; j++)
    data[j]=0;

  // skip the header
  word = (int *)( byte+HEADERSET);  

  // initialise the storePtr to the beginning of databuf
  if(*flagPtr ==0)  
  {
    // move away from infomation head and get the right entry pointer
    rampinfoPtr=(SERVO_DATAHEAD *)databuf;
    storePtr1=(int *)(rampinfoPtr +1);
    *flagPtr=1;
  }
  // get all the framedata (no headset) into dataBuf  
  for( row=0; row<ROW_NUM;row++)
  {
    for(ch=0; ch<COL_NUM; ch++ )
    {
      j=row*COL_NUM +ch;
      data[j]=(*(word +j))*setup->colMask[ch];
      *(storePtr1+j)=data[j]; 
    }
  }
  // advance storePtr by
  storePtr1 +=setup->totalPixel;
  //save in binary mode for each heater point
  fwrite(byte,1,myInfo->bufSize,myInfo->fpSq1);

  // ************* for stripchart for this selRow ****************
  // *************************************************************
  // the data format for TESTRANSIT  is
  // col0   COL1   COL2   COL3   COL4   COL5   COL6   COL7   COL8   COL9 
  // 1     fdbk0  data0  data1  data2  data3  data4  data5  data6  data7  
  // 2     fdbk1  data0  data1  data2  data3  data4  data5  data6  data7  
  //   ...................
  fprintf(myInfo->fpStrchart, "%7d ",stripchFlag);
  fprintf(myInfo->fpStrchart, "%7d ",setup->fdbk);

   // the data has total 32 cols, always use the selRow data [0:40]
  // FDBKSHIFT(14) FDBKMASK(0x3FFFF) defined in sc2da_par.h 
  for(ch=0; ch< COL_NUM;ch++)
  {
    j=setup->selRow*COL_NUM +ch;
    if (setup->dataMode==4)
    {
      data1  = ( data[j] >> FDBKSHIFT) & FDBKMASK;
      sc2dalibsetup_signExt(myInfo,&data1,SIGNBIT_18BIT,DATAMASK_17BIT,
                   SIGN_EXT_18,status);
    }
    fprintf ( myInfo->fpStrchart, "%11d ",data1);
  }	
  fprintf(myInfo->fpStrchart, "\n");
  fflush (myInfo->fpStrchart);
  stripchFlag++;
}


/**
 * \fn void sc2datransitTES(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *     ARRAYSET *setup,mcexml_struct *mceInxpt,
 *     char * dateTime, StatusType *status)
 *
 * \brief function
 *  send out heater and tes bias (if BIASCHANGED) values to MCE  
 *
 * \param  con    SDSU_CONTEXT pointer
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  mceInxpt mcexml_struct pointer
 * \param dateTime char pointer for date Time string
 * \param status   StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
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
)
{
  dasCmdInfo_t  biasCmd;
  int           wait;

  if (!StatusOkP(status))

  strcpy(biasCmd.mceCmd,"");

  // ramp heater Val, 
  // set heater
  sprintf ( biasCmd.mceCmd, "wb bc1 bias %d",setup->fdbk ); 
  sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
  if ( !StatusOkP(status) )
  {
    ErsRep (0, status,
      "sc2dasetbiasfdbkSQ1: sc2dalib_sendCmd failed HEATER Val (%s)",
       biasCmd.mceCmd); 
    return;
  }
  // set tesBias
  if(setup->biasFlag==BIASCHANGED)
  {
    sprintf ( biasCmd.mceCmd, "wb bc2 bias %d",setup->bias ); 
    sc2dalib_sendCmd(con,myInfo,&biasCmd,mceInxpt,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,
       "sc2dasetbiasfdbkSQ1: sc2dalib_sendCmd failed TESBIAS Val (%s)",
        biasCmd.mceCmd); 
      return;
    }
  }
  // let heater settle-down if  slopSelect[12]!=0
 if( setup->slopSelect[12] > 0 )
 {
   wait=1000000*setup->slopSelect[11]/setup->slopSelect[12]; 
   usleep(wait);
   //printf("wait (%d)\n",wait);
 }
}



/**
 * \fn void sc2dafindTransit(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *  char *tesData, StatusType *status)
 *
 * \brief funcation:
 *  calculate transit value, use non-filted data
 *
 * \param myInfo  dasInfo structure pointer
 * \param setup   ARRAYSET structure pointer
 * \param tesData  char pointer for testransit data buffer
 * \param status  StatusType pointer.  given and return
 *
 */
/*+ sc2dafindTransit
*/
void sc2dafindTransit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
char            *tesData,
StatusType      *status
)
{
  int    i,heatModulate,heatLvl; 
  int    framesizeperheatLvl;
  char   *data,*dataPtr;
  int    heatStart; 
  SERVO_DATAHEAD *rampinfoPtr;
  FILE   *fp;   
  static  char *dataperheatLvl;

  if (!StatusOkP(status)) return;

  // the triangle wave modulation steps on heat Level
  dataperheatLvl=NULL;
  heatLvl=setup->biasNo;
  heatStart=setup->minFDBK;
  heatModulate=setup->fdbkNo*4+1; 
  setup->totalPixel=ROW_NUM*COL_NUM;

  // point to the framedata inside tesData buffer
  rampinfoPtr=(SERVO_DATAHEAD *)tesData;
  data=(char *)(rampinfoPtr +1);
  
  // use the full frame size now, the data has no frameState, no CHKSUM
  framesizeperheatLvl=setup->totalPixel*heatModulate*sizeof(int);
  dataperheatLvl=(char*)calloc(framesizeperheatLvl,1);
  if (dataperheatLvl==NULL)
  {
    *status=DITS__APP_ERROR;
    ErsRep (0, status,"failed to calloc dataperheatLvl\n");
    return;
  }
 
  // open a file for tmp: print triangle Val
  sprintf(msg,"%s-pixel",myInfo->dataFile);
  if((fp = fopen(msg, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dafindTransit: Error- failed to open file %s", msg); 
      return;
    }

  // step through heaterLevel
  for (i=0;i<heatLvl;i++)
  {
    //printf("call sc2dalibsetup_heaterGen\n");
    sc2dalibsetup_heaterGen(setup,1,status);
  
    dataPtr =  data + i*framesizeperheatLvl;
    memcpy(dataperheatLvl,dataPtr,framesizeperheatLvl);
  
     // use the pixelsize not full frameSize,  save the selRow
    //sc2dapixelTransit(myInfo,setup,dataperheatLvl,setup->selRow,
    //                     i,0,fp,status);
    if (*status !=STATUS__OK) 
    {
      ErsRep (0, status,"failed to call sc2dapixelTransit\n");
      free(dataperheatLvl);
      fclose(fp);
      return ; 
    } 
    setup->minFDBK -=setup->stepBIAS;
  }
  
  setup->minFDBK=heatStart;
  sc2dalibsetup_transitResult(myInfo, setup,fp,status);
  printf("complete  sc2dalibsetup_transitResult\n");
  fclose(fp);   
  free(dataperheatLvl);
  printf("free fp, dataperheatLvl\n");
}


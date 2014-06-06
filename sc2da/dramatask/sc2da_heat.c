/**
 * \file sc2dalib.c
 *
 * \brief collection heater associated functions
 *
 * \author (copier) Bryan Gorges
 *
 * \version 1.0.0
 */

//#include "sc2math.h"   // for sc2math-linfit
// #include "math.h"      // for pow
//#include "sc2da_par.h"
//#include "sc2da_struct.h"
// #include "jitXML.h"
//#include "sc2sqopt_par.h"
//#include "sc2sqopt.h"
//#include "star/atl.h"
// #include "sc2dalib.h"
#include "sc2da_heat.h"
// #include "sc2headman.h"
// #include "sc2dalibsetup.h"

static int stairCounter;

/**
 * \fn void heat_init_bias(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  int *flag, StatusType *status)
 *
 * \brief function
 *  Based on values of drcontrol and pixelHeat we are either going to do
 *  A standard SEQUENCE where the heater current is not changed, a slow
 *  fixed flat field where the current is set to a single fixed value, 
 *  a fast flat field where the heater is changed within the SEQUENCE, a
 *  heater ramp where the heater is changed within SEQUENCE, a bias sawtooth
 *  or TES bias ramp where the TES bias is changed within SEQUENCE 
 *
 *  Determine the style of SEQUENCE and handle the heater 
 *  current and TES Bias accordingly
 *
 * \param con          SDSU context structure
 * \param myInfo       dasInfo structure pointer
 * \param flag         False if nothing changes within sequence > 0 otherwise
 * \param status       StatusType     
 *
 * Used to be sc2dalib_initHeatBiasHandling, may not be used.
 */
/*+ heat_init_bias
*/
void heat_init_bias
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,   
int                   *flag,
StatusType            *status
)
{
  float ampsPerCount = 3.7842374e-10; //24.8 microamps for 65535 counts of the DtoA
  float nominal, currentWatts;
  float offset, desired;
  double squared;
  int result;
  char       mceCmd[FILE_LEN];
  char       heatVal[]="rb bc1 bias 1";
  char       biasVal[]="rb bc2 bias 1";
  int stairHeightCnts;
  int currentSetHeat, currentSetBias;

  if (*status != STATUS__OK) return;

  *flag = 0;
  stairHeightCnts = 0;
  stairCounter = 0;

  /* Read the current heater setting */
  sc2dalib_readmceVal(con, myInfo, mceInxpt, heatVal, &currentSetHeat, 1, status);
  if (!StatusOkP(status)) 
    {
      ErsRep(0,status, "heat_init_bias: sc2dalib_readmceVal(1) failed to read heater value"); 
      return;
    }

  /* Read the current TES bias setting */
  sc2dalib_readmceVal(con, myInfo, mceInxpt, biasVal, &currentSetBias, 1, status);
  if (!StatusOkP(status)) 
    {
      ErsRep(0,status, "heat_init_bias: sc2dalib_readmceVal(1) failed to read bias value"); 
      return;
    }
  myInfo->detbias = currentSetBias;

  /* Test to see if it is a fast flat field (heater sawtooth) */
  if( (myInfo->drcontrol & FLATFIELD_BIT) != 0)
    {
      *flag = 1;

      /* The stairNum we were given is the total count, but for heater sawtooths 
         it really needs to be the number of steps in a quarter cycle */
      myInfo->stairNum /= 4;

      /* Use the current heater setting and stairHeight in femto watts 
	 to calculate the step in D/A counts */

      /* Calculate:
	 currentWatts    -> The current Pixel Heater setting in watts
	 offset          -> requested offset in watts (stairHeight is in femto watts)
	 desired         -> The new desired setting in watts
	 result          -> The new D/A setting that will yield the desired one 
	 stairHeightCnts -> The different in D/A counts between result and currentSetHeat 

	 This assumes a 2 ohm resistor */

      currentWatts = (((float)currentSetHeat * ampsPerCount)*
		      ((float)currentSetHeat * ampsPerCount)) * 2.0;
      offset = (float)myInfo->stairHeight * 1.0e-15;
      desired = currentWatts + offset;
      squared = desired / 2.0;
      if(squared > 0)
	{
	  result = sqrt(squared) / ampsPerCount;
	}
      else
	{
	  result = 1000;
	}

      stairHeightCnts = result - currentSetHeat;
      if(stairHeightCnts < 1) stairHeightCnts = 1;
      jitDebug(2 ,"currentWatts = %e offset %e desired %e squared %e result %d currentSetHeat %d height %d",
	       currentWatts, offset, desired, squared, result, currentSetHeat, stairHeightCnts);

      /* Now set the heater current to the result */

      jitDebug(2," Set heater: wb bc1 bias %d \n", result);

      sprintf(mceCmd, "wb bc1 bias %d", result);
      // send cmd and get reply
      sc2dalib_setmceVal(con,myInfo,mceInxpt,mceCmd,status);
      if ( !StatusOkP(status) )
	{
	  ErsRep(0,status,"heat_init_bias: sc2dalib_setmceval(1) %s failed",mceCmd); 
	  return;
	}
      /* fprintf(myInfo->fpLog,"heat_init_bias FLATFIELD %s\n",mceCmd); */

      myInfo->stairQCycleCount = 0;
      myInfo->stairNumCount = 0;
      myInfo->stairPresentValue = result;
      myInfo->stairHeightCnts = stairHeightCnts;
      myInfo->pixelHeat = currentSetHeat;

    } /* End of fast flat field initialization */
  
  /* Is it a heater ramp SEQUENCE? */
  else if( (myInfo->drcontrol & HEATRAMP_BIT) != 0)
    {
      *flag = 2;

      /* The stairNum we were given is the total count, but for heater ramps 
         it really needs to be the number of steps in a half cycle */
      myInfo->stairNum /= 2;

      /* Heater ramps are done in DAC counts staring at stairStart stepping down stairNum/2
         steps and the stepping back up again */

      myInfo->stairHeightCnts = myInfo->stairHeight;
      stairHeightCnts = myInfo->stairHeight;
      myInfo->pixelHeat = myInfo->stairStart;
      currentSetHeat = myInfo->stairStart;

      /* Now set the heater current to this value */

      jitDebug(2," Set heater: wb bc1 bias %d \n", currentSetHeat);

      sprintf(mceCmd, "wb bc1 bias %d", currentSetHeat);
      // send cmd and get reply
      sc2dalib_setmceVal(con,myInfo,mceInxpt,mceCmd,status);
      if ( !StatusOkP(status) )
	{
	  ErsRep(0,status,"heat_init_bias: sc2dalib_setmceval(2) %s failed",mceCmd); 
	  return;
	}
      /* fprintf(myInfo->fpLog,"heat_init_bias HEATRAMP %s\n",mceCmd); */

      myInfo->stairHCycleCount = 0;
      myInfo->stairNumCount = 0;
      myInfo->stairPresentValue = currentSetHeat;

    }

  /* Is it a TES Bias sawtooth SEQUENCE? */
  else if( (myInfo->drcontrol & BIASSAW_BIT) != 0)
    {
      *flag = 3;

      /* Bias sawtooths are done in DAC counts starting at the current value (plus one step)
         and then sawtoothing around the current value and ending up at it */

      myInfo->stairHeightCnts = myInfo->stairHeight;
      stairHeightCnts = myInfo->stairHeight;

      /* The stairNum we were given is the total count, but for bias sawtooths 
         it really needs to be the number of steps in a quarter cycle */
      myInfo->stairNum /= 4;

      /* Using the current bias setting and 
         stairHeight in DAC units calculate the first setting */

      result = currentSetBias + stairHeightCnts;

      /* Now set the TES bias to this value */

      jitDebug(2," Set TES Bias: wb bc2 bias %d \n", result);

      sprintf(mceCmd, "wb bc2 bias %d", result);
      // send cmd and get reply
      fprintf(myInfo->fpLog,"heat_init_bias BIASSAW %s\n",mceCmd);
      sc2dalib_setmceVal(con,myInfo,mceInxpt,mceCmd,status);
      if ( !StatusOkP(status) )
	{
	  ErsRep(0,status,"heat_init_bias: sc2dalib_setmceval(3) %s failed",mceCmd); 
	  return;
	}

      myInfo->stairQCycleCount = 0;
      myInfo->stairNumCount = 0;
      myInfo->stairPresentValue = result;
      myInfo->pixelHeat = currentSetHeat;
    }

  /* Is it bias ramp SEQUENCE? */
  else if( (myInfo->drcontrol & BIASRAMP_BIT) != 0)
    {
      *flag = 4;

      /* The stairNum we were given is the total count, but for bias ramps 
         it really needs to be the number of steps in a half cycle */
      myInfo->stairNum /= 2;

      /* Bias ramps are done in DAC counts staring at stairStart stepping down stairNum/2
         steps and then stepping back up again */

      myInfo->stairHeightCnts = myInfo->stairHeight;
      stairHeightCnts = myInfo->stairHeight;
      myInfo->detbias = myInfo->stairStart;
      currentSetBias = myInfo->stairStart;

      /* Now set the TES Bias to this value */

      jitDebug(2," Set TES Bias: wb bc2 bias %d \n", currentSetBias);

      sprintf(mceCmd, "wb bc2 bias %d", currentSetBias);
      // send cmd and get reply
      fprintf(myInfo->fpLog,"heat_init_bias BIASRAMP %s\n",mceCmd);
      sc2dalib_setmceVal(con,myInfo,mceInxpt,mceCmd,status);
      if ( !StatusOkP(status) )
	{
	  ErsRep(0,status,"heat_init_bias: sc2dalib_setmceval(4) %s failed",mceCmd); 
	  return;
	}

      myInfo->stairHCycleCount = 0;
      myInfo->stairNumCount = 0;
      myInfo->stairPresentValue = currentSetBias;
      myInfo->pixelHeat = currentSetHeat;
    }

  else /* It will either be a standard SEQUENCE, or a slow flat field */
    {
      /* test to see if it is a standard sequence */

      if ( myInfo->pixelHeat == -99999 )
	{
	  myInfo->pixelHeat = currentSetHeat; 
	}
      else /* This is going to be a slow flat field */
	{

	  /* Calculate:
	     nominal -> The nominal Pixel Heater setting in watts
	     offset  -> requested offset in watts (pixelHeat is in femto watts)
	     desired -> The new desired setting in watts
	     result  -> The new D/A setting that will yield the desired one */

	  nominal = (((float)myInfo->nominalPixelHeat * ampsPerCount)*
		     ((float)myInfo->nominalPixelHeat * ampsPerCount)) * 2.0;
	  offset = (float)myInfo->pixelHeat * 1.0e-15;
	  desired = nominal + offset;
	  squared = desired / 2.0;
	  if(squared > 0)
	    {
	      result = sqrt(squared) / ampsPerCount;
	    }
	  else
	    {
	      result = 1000;
	    }
	  currentSetHeat = myInfo->pixelHeat = result;

	  /* Now set the heater current to this value */

	  jitDebug(2," Set heater: wb bc1 bias %d \n", myInfo->pixelHeat); 
	  sprintf(mceCmd, "wb bc1 bias %d", myInfo->pixelHeat);
	  // send cmd and get reply
	  /* fprintf(myInfo->fpLog,"heat_init_bias slow flat %s\n",mceCmd);*/
	  sc2dalib_setmceVal(con,myInfo,mceInxpt,mceCmd,status); 
	  if ( !StatusOkP(status) )
	    {
	      ErsRep(0,status,"heat_init_bias: sc2dalib_setmceval(2) %s failed",mceCmd); 
	      return;
	    }
	  // sleep for 10ms for heater to settle down
	  usleep(10000);
	}
    }

  sc2headman_initHeatBiasHandling(*flag, currentSetHeat, currentSetBias, stairHeightCnts, myInfo->stairNum,
				myInfo->stairWidth, status);
}

/**
 * \fn void heat_init_slope(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  char *dateTime, StatusType *status)
 *
 * \brief function: 
 *  read in all args from the heater setup .
 *
 * \param con           SDSU context structure
 * \param myInfo       dasInfo structure pointer
 * \param dateTime      string pointer to dateTime string
 * \param status        StatusType     
 *
 * Used to be sc2dalib_heaterslopeInit
 */
/*+ heat_init_slope
*/
void heat_init_slope
(
SDSU_CONTEXT      *con,
dasInfoStruct_t   *myInfo, 
char              *dateTime,
StatusType         *status
)
{
  long       in_sequence;
  long       range[]={0,3};
  DitsArgType  argId;

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  sc2dalib_updateDebug(con,myInfo,status);
 
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"heat_init_slope: %s has not completed",
             seqStatus[myInfo->actionFlag]);
    return;
  } 

  myInfo->actionFlag=HEATERSLOPE;

  argId = DitsGetArgument();
  jitArgGetI(argId, "SVFILE_FLAG", 1, range, 2, 0,
        &myInfo->filesvFlag, status );
  //default 2 ?
  jitArgGetS(argId,"DATA_FILE",2,NULL,"data.txt",0, FILE_LEN,
            myInfo->dataFile, NULL,status );
  jitArgGetS (argId, "SETUP_FILE", 3,NULL, "setup-heaterslope", 0, FILE_LEN,
              myInfo->batchFile,  NULL,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"heat_init_slope: failed to get variables.");
    return;
  }
  sc2dalibsetup_servoreadsetupWrap(myInfo,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"heat_init_slope: sc2dalibsetup_servoreadsetupWrap failed"); 
    return;
  }
  my_fclose(&(myInfo->fpBatch));
  if((myInfo->fpBatch = fopen(myInfo->batchFile,"r")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "Error- heat_init_slope: failed to open file %s", myInfo->batchFile); 
      return;
    }
  jitDebug(16,"heat_init_slope: read %s\n",myInfo->batchFile);

  sc2dalibsetup_readheaterSetup(myInfo, status);

  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"heat_init_slope:sc2dalibsetup_readheaterSetup failed");
    // in sc2dalib__heatslope call sc2dalib_actionfileEnd to close file
    return;
  }
  my_fclose(&(myInfo->fpData));
  if((myInfo->fpData = fopen( myInfo->dataFile, "a" )) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "Error- heat_init_slope:2 failed to open file %s", myInfo->dataFile); 
      return;
    }
  // use batchFile for storing the mean pixel data 
  sprintf(myInfo->batchFile,"%s-row%d-col%d",myInfo->dataFile, 
          myInfo->heatSlp.row, myInfo->heatSlp.col);

  my_fclose(&(myInfo->fpBatch));
  if((myInfo->fpBatch = fopen(myInfo->batchFile, "w")) == NULL)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
     "heat_init_slope: failed to open %s, check permission",
           myInfo->batchFile); 
    // in sc2dalib__heatslope call sc2dalib_actionfileEnd to close files
    return;
  }
  // use strchartFile for storing the pixel data ( every point) 
  if(myInfo->filesvFlag >0 )
  {
    sprintf(myInfo->strchartFile,"%s-row%d-col%d-pixel",myInfo->dataFile, 
          myInfo->heatSlp.row, myInfo->heatSlp.col);
    my_fclose(&(myInfo->fpStrchart));
    if((myInfo->fpStrchart = fopen(myInfo->strchartFile, "w")) == NULL)
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,
       "heat_init_slope: failed to open %s, check permission",
           myInfo->strchartFile);
      // in sc2dalib__heatslope call sc2dalib_actionfileEnd to close file
      return;
    }
  }
}

/**
 * \fn void heat_slope(dasInfoStruct_t *myInfo, 
 *  int *heaterMask, double *heaterSlope, StatusType *status)
 *
 * \brief function: 
 *  find the pixel heater slope from the heater data.
 *
 * \param myInfo       dasInfo structure pointer
 * \param heaterMask int pointer, get slope if pixelMask=1,
 * \param  heaterSlope double pointer
 * \param status        StatusType     
 *
 * Used to be: sc2dalib_heaterSlope
 */
/*+ heat_slope
*/
void heat_slope
(
dasInfoStruct_t  *myInfo,
int              *heaterMask,
double           *heaterSlope,
StatusType       *status
)
{
  int     howmany, pixelData,onepixel;
  int     i, nStep,pixel, totalPixel;
  int     svpixelFlag;
  double  *chData, *ithchData;
  double  *pixelPtr,*powerPtr,*weight;
  double  heat,dcLvl,slope;
  double  fullDC=65535, pixelResist=2, current=24.8, c;

  if (*status != STATUS__OK) return;

  nStep=myInfo->heatSlp.nStep+1;

  /* Current is full scale microamps
     fullDC is full scale DtoA counts
     pixelResist is nominal pixel heater resistor in ohms
     So current/fullDC is microamps per count
     then c equals: R * I^2 or microwatts per DtoA count */

  c= pixelResist*pow(current/fullDC,2);

  chData=(double *)calloc( nStep*ROW_NUM*COL_NUM, sizeof(double));
  if (chData==NULL)
  {  
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"heat_slope: failed to calloc for nStep*chdata\n"); 
    return;
  }

  pixelPtr=(double *)calloc( nStep, sizeof(double));
  if (pixelPtr==NULL)
  {  
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"heat_slope: failed to calloc for pixel data\n");
    free(chData); 
    return;
  }

  powerPtr=(double *)calloc( nStep, sizeof(double));
  if (powerPtr==NULL)
  {  
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"heat_slope: failed to calloc for powerVal\n"); 
    free(pixelPtr);
    free(chData);
    return;
  } 

  weight=(double *)calloc( nStep, sizeof(double));
  if (weight==NULL)
  {  
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"heat_slope:failed to calloc for weight data\n"); 
    free(pixelPtr);
    free(powerPtr);
    free(chData); 
   return;
  }

  svpixelFlag=myInfo->filesvFlag;
  if ( myInfo->filesvFlag >0 )
    myInfo->filesvFlag=1;

  // loop through all nStep
  for (i=0; i<nStep; i++)
  {

    printf("heat_slope working on step: %d\n",i);

    // myInfo->svfileFlag=3 save only individual pixel from the first 
    // heat file, myInfo->svfileFlag=1 save all
    if (svpixelFlag==3 && i==1)
      myInfo->filesvFlag=0;
      
    // use cmdrepFile for the heater data file
    sprintf(myInfo->cmdrepFile,"%s/%s-%d-eng.txt",getenv("CURRENTDATADIR"),
            myInfo->heatSlp.base, i);
    my_fclose(&(myInfo->fpOtheruse));
    if ((myInfo->fpOtheruse = fopen( myInfo->cmdrepFile, "r" )) == NULL)
    {  
      *status=DITS__APP_ERROR;
      ErsRep(0,status,
        "heat_slope: failed to open %s, check file permission",
        myInfo->cmdrepFile); 
      // in sc2dalib__heatslope call sc2dalib_actionfileEnd to close file
      free(pixelPtr);
      free(powerPtr);
      free(chData);
      free(weight);
      return ;
    }

    /* heat is starting setting (in DtoA counts) minus the size 
       of the step times i (loop counter) */
    heat=myInfo->heatSlp.stVal - myInfo->heatSlp.step*i;
    /* If we now multiply the heater DtoA setting by microwatts per DtoA count
       we get the microwatts being put into the pixel by the heater */
    heat= c * pow(heat, 2);
    fprintf(myInfo->fpData,"%11.2f", heat);
    fprintf(myInfo->fpBatch,"%11.2f",heat);

    ithchData=chData+ i*ROW_NUM*COL_NUM;

    // return the averaged frame in ithchData, a single pixel data in pixelData
    // return howmany-frame for heaterStep=i
    sc2dalibsetup_readheaterData(myInfo, &howmany, ithchData, &pixelData, heaterMask,status);
    if ( !StatusOkP(status) ) 
    {
      free(pixelPtr);
      free(powerPtr);
      free(weight);
      free(chData); 
      // in sc2dalib__heatslope call sc2dalib_actionfileEnd to close file
      ErsRep(0,status,"heat_slope: sc2dalib_readheaterData failed"); 
      return;
    }
    /* Note that weight does not change it is a constant 1
       and the power is being written to powerPtr */

    powerPtr[i]=heat;
    weight[i]=1;

    my_fclose(&(myInfo->fpOtheruse));
    
    MsgOut(status," %d frames average completed for (%s-%d)",
           howmany,  myInfo->heatSlp.base, i);

    jitDebug(2," %d frames from (%s) %s-%d, meanVal to %s %s\n",
      howmany,  getenv("CURRENTDATADIR"), myInfo->heatSlp.base, i,
      myInfo->dataFile, myInfo->batchFile);
  }

  //now, we have power[0..nstep], weight[0..nstep] and chData[0..nstep][ROW*COL]
  // go through pixels, but not dark row

  totalPixel=(ROW_NUM-1)*COL_NUM;
  onepixel=myInfo->heatSlp.row*myInfo->heatSlp.col;
  for (pixel=0; pixel<totalPixel; pixel++)
  {
    if (heaterMask[pixel]!=0) 
    {
      for (i=0; i<nStep; i++)
      {
        pixelPtr[i]= *(chData + i*ROW_NUM*COL_NUM + pixel);
  printf("Pixel: %d  power %f  sq1FDBK  %f \n",pixel,powerPtr[i],pixelPtr[i]);
      }
      
      /* linear fit of the Sq1FDBK units per microwatt */
      sc2math_linfit(nStep, powerPtr, pixelPtr,weight, &slope,&dcLvl,status);
      printf("Pixel: %d  slope %f  offset  %f \n",pixel,slope,dcLvl);
      heaterSlope[pixel]=slope;
    }
    else
     heaterSlope[pixel]=0;

    if (pixel == onepixel )
    {
       myInfo->heatSlp.slope=slope;
       myInfo->heatSlp.refFDBK=(int)pixelPtr[myInfo->heatSlp.ref];
    }
  }
  free(pixelPtr);
  free(powerPtr);
  free(weight);
  free(chData);
  heat_save_slope(myInfo, heaterSlope,status);
  return; 
}

/**
 * \fn void heat_read_slope(dasInfoStruct_t *myInfo, 
 *  double *pixelSlope, StatusType *status)
 *
 * \brief function
 *   Opens the file heaterslope.xml
 *   Reads the contents of that file into an SDS structure
 *   Walks the structure row by row and places the values read into pixelSlope
 *   I think this array is a global parameter called heaterSlope in sc2da.c
 *
 *  get args for HEATER_TRACK action and fill pixelSlope array 
 *
 * \param myInfo      dasInfo structure pointer
 * \param pixelSlope  double pointer for pixelSlope array
 * \param status       StatusType     
 *
 * Used to be: sc2dalib_heaterslopeRead
 */
/*+ heat_read_slope
*/
void heat_read_slope
(
dasInfoStruct_t   *myInfo,    
double            *pixelSlope,
StatusType        *status
)
{
  SdsIdType slopeId=0;   //id for HEATER_SLOPE structure 
  SdsIdType rowId;       //  id for row structure 
  SdsIdType id = 0;
  SdsCodeType code;
  
  char  name[FILE_LEN];
  char  delimiters[]= " ",*token;
  char  cardName[FILE_LEN], rc[50];

  long  ndims;            /* number of dimensions in row */
  long  rowNo;
  int   j,i, pixel,card;
  unsigned long dims[7];       
  unsigned long indims[7];

  if (*status != STATUS__OK) return;   

  // the heaterslope.xml needs to be in $SC2SCRATCH
  sprintf (name, "%s/%s.xml", getenv ( "SC2SCRATCH" ),PIXELHEATER_SLOPE );

  jitXML2Sds ( strlen(name),name, &slopeId,status );
  if (!StatusOkP(status))
  {
    ErsRep(0,status, "heat_read_slope: jitXML2Sds failed to read file: %s",name); 
   return;
  }
  //MsgOut(status," pixelSlopeXML=%s",name);
 
  //SdsList( slopeId,status);

  //slopeId is top SDS, the rest is sub sds. get sub-sds 
  dims[0] = dims[1] = 0;

  SdsFind ( slopeId, "row", &rowId, status );
  SdsInfo ( rowId, name, &code, &ndims, dims, status );  
  if (!StatusOkP(status))
  {
    ErsRep(0,status,"heat_read_slope: find row failed"); 
    return;
  }
  jitDebug(2,"rowId's name =%s ndims=%d, dims[0]=%d dims[1]=%d\n",
             name, ndims, dims[0],dims[1]);

  // slope=0: initialised during load task
  pixel=0;
  for ( j=0; j<dims[0]; j++ )
  {
    indims[0] = j + 1;

    // from rowId ( SDS array), get each item from row structure 
    SdsCell ( rowId, 1, indims, &id, status ); 

    // from each array[x], get item, 
    ArgGeti(id, "Count", &rowNo,status);
    if (StatusOkP(status))
    {
      pixel=(int)rowNo*32;
       
      for (card=1; card<=4; card++)
      {
        sprintf(cardName,"rc%d",card);
        ArgGetString(id, cardName, FILE_LEN, rc, status );
        if (StatusOkP(status))
        {
          token = strtok (rc, delimiters);
          pixelSlope[pixel]=atof(token);
          for ( i=1; i<8; i++)
          {
            pixel ++;
            token=strtok(NULL,delimiters);
            pixelSlope[pixel]=atof(token);
          }
          pixel ++;
        }
        else
          ErsRep(0,status, "heat_read_slope: failed in cardName row(%d) pixel(%d)",
                  j,pixel);
      }
    }
    if (!StatusOkP(status))
    {
      ErsRep(0,status,"heat_read_slope: failed in Count" ); 
      return;
    }
    SdsFreeId ( id, status );
  }
  SdsFreeId ( rowId, status );
  SdsFreeId ( slopeId, status );
}

/**
 * \fn void heat_save_slope(dasInfoStruct_t *myInfo,
 *  double *heaterSlope, StatusType *status)
 *
 * \brief function
 *  save heater slope and ref feed back to heaterslopeTable 
 *
 * \param  myInfo dasInfoStruct_t pointer
 * \param  heaterSlope double pointer
 * \param  status StatusType.  given and return
 *
 * Used to be sc2dalib_heaterslopeSave
 */
/*+ heat_save_slope
*/
void heat_save_slope
(
dasInfoStruct_t *myInfo,
double          *heaterSlope,
StatusType      *status

)
{
  char   tmpfile1[FILE_LEN];
  int    row, card,pixel,i;
  FILE   *fp;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile1, "%s/%s.xml", getenv ( "SC2SCRATCH" ),PIXELHEATER_SLOPE );
  if((fp=fopen(tmpfile1, "w")) == NULL )
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status, "heat_save_slope: Error- failed to open file %s", tmpfile1); 
    return;
  }
  fprintf(fp,"<!--- this is extracted after mceheaterslope    -->\n");
  fprintf(fp,"<!--- slope !=0 means the pixel's power values  -->\n");
  fprintf(fp,"<!--- is used for calculating updated heaterVal -->\n\n");
  fprintf(fp,"<HEATERSLOPE>\n");
  fprintf(fp,"  <REFHEAT>%d </REFHEAT>\n",myInfo->heatSlp.refHeat);
  fprintf(fp,"  <REFFDBK>%d </REFFDBK>\n",myInfo->heatSlp.refFDBK);
  // go through pixels
  for (row=0; row<40; row ++)
  {  
    fprintf(fp,"  <row Count=\"%d\"\n",row);
    pixel=row*COL_NUM;
    for (card=1; card<=4; card ++)
    {
      fprintf(fp,"    rc%d=\"",card);
      for ( i=0; i<8; i++)
      {
        fprintf(fp,"%11.3f ",heaterSlope[pixel]);
        pixel ++;
      }
      fprintf(fp,"\"\n"); 
    }
    fprintf(fp,"  />\n");     
  }
  fprintf(fp,"</HEATERSLOPE>\n");
  fclose(fp);
}

/**
 * \fn void heat_init_track(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *   dasCmdInfo_t *myCmd, struct mcexml_struct  *mceinxPtr,
 *   char *dateTime, StatusType *status)
 *
 * \brief function
 *  get args for MCETRKHEAT action and open file
 *
 * \param con       SDSU context structure pointer 
 * \param myInfo    dasInfoStruct_t pointer
 * \param myCmd      dasCmdInfo_t pointer
 * \param mceinxPtr   mcexml_struct pointer for parameter lookup table
 * \param dateTime   dateTime string pointer         
 * \param status    StatusType     
 *
 * Used to be sc2dalib_trkheatInit
 */
/*+ heat_init_track
*/
void heat_init_track
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
char                  *dateTime,
StatusType            *status
)
{
  long         in_sequence,configured;
  long         range[]={0,2};
  DitsArgType  argId;
  char         mceCard[FILE_LEN];
  char         tmp[FILE_LEN],tmp1[FILE_LEN];
  char         heatVal[]="rb bc1 bias 1";

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 

  sc2dalib_updateDebug(con,myInfo, status);

  // Make sure we have been configured and we are not in SEQUENCE

  SdpGeti("CONFIGURED", &configured, status);
  if (configured==0 )
    {
      *status = DITS__APP_ERROR;
      ErsRep(0,status,"heat_init_track: the DA is not configured" ); 
      return;
    }
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
    {
      *status = DITS__APP_ERROR;
      ErsRep(0,status,"heat_init_track: %s has not completed",
             seqStatus[myInfo->actionFlag]);
      return;
    } 
  myInfo->actionFlag=TRKHEATACTION;

  /* Using the argument sent to us from the SCUBA2 RTS client set the file save flag
     data file name and batch file name */

  argId = DitsGetArgument();
  jitArgGetI(argId, "SVFILE_FLAG", 1, range, 2, 0, &myInfo->filesvFlag,
       status );
  // filesvFlag=1 save cmdreply 
  // filesvFlag=2 don't save cmdreply, but save data, default? (changed ^)
  // for heater tracking, we only save data
  if(myInfo->filesvFlag==1)
    {
      myInfo->filesvFlag=2;
    }

  jitArgGetS(argId, "DATA_FILE", 2, NULL, "data.txt", 0, FILE_LEN,
       tmp, NULL, status);
  //print ?
  sprintf(myInfo->dataFile, "%s/%s",getenv("CURRENTDATADIR"), tmp);
  jitArgGetS(argId, "SETUP_FILE", 3, NULL, "setup-trkheater", 0, FILE_LEN,
       tmp1, NULL, status);
  sprintf(myInfo->batchFile, "%s/%s",getenv("CONFIG_ALL"), tmp1);

  /* Read The NUM_TIMES_THRU parameter from the argument */
  range[1] = 45; /* We really only think one or two */
  jitArgGetI( argId, "NUM_TIMES_THRU", 1, range, 0, 0, &myInfo->numTimesThru,
        status );
  if( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2dalib_trkheaterInit: failed to get variables.");
      return;
    }

  // now get the args from setup file
  // put all include=xxx in the setup file into a single file  
  // copy it to batchFile, for sc2dalibsetup_readheaterSetup

  sc2dalibsetup_servoreadsetupWrap(myInfo,status);
  if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2dalib_trkheaterInit: sc2dalibsetup_servoreadsetupWrap failed"); 
      return;
    }
  my_fclose(&(myInfo->fpBatch));
  if((myInfo->fpBatch=fopen(myInfo->batchFile,"r")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "heat_init_track: Error- failed to open file %s", myInfo->batchFile); 
      return;
    }
  jitDebug(16,"heat_init_track: read %s\n",myInfo->batchFile);

  /* Read in the heater slopes and the heater setup */

  heat_read_slope(myInfo, heaterSlope, status);

  sc2dalibsetup_readheaterSetup(myInfo,status);
  if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2dalib_trkheaterInit: sc2dalib_readheaterSetup failed"); 
      return;
    }
  jitDebug(16,"heat_init_track: finish readheaterSetup\n");
  my_fclose(&(myInfo->fpBatch));
  
  // Build up the goCmd always use rcs card
  sprintf(mceCard, "rcs");
  sprintf(myInfo->goCmd, "GO %s ret_dat 1",mceCard);

  sc2dalibsetup_whichRC(myInfo,mceCard,status);
  if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"heat_init_track: not recognised(%s) as MCE_WHICHCARD",
       mceCard);
      return;
    }
  strcpy(myCmd->mceCmd,myInfo->goCmd);

  sc2dalib_getcmdBuf(myInfo,myCmd,mceinxPtr,status);
  if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"heat_init_track: sc2dalib_getCmdBuf failed"); 
      return;
    }

  // Set stripchart file name for recording and displaying tracking results to trkheat.txt 

  if( (myInfo->heatSlp.option & 0x01) != 0 )
  {
    sprintf(myInfo->strchartFile,"%s/trkheat.txt", getenv("ORAC_DATA_OUT") );
  }

  // This deletes the batchfile called tmp which was made in sc2dalibsetup_servoreadsetupWrap
  // And then used by sc2dalibsetup_readheaterSetup
  sprintf( tmp, "rm -f %s",myInfo->batchFile);
  system ( tmp);

  sc2dalib_openFiles(myInfo,DATFILEAPPEND,NOBATCHFILE,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"heat_init_track: sc2dalib_openFiles failed"); 
    return;
  }
    
  fprintf(myInfo->fpLog,"\n<%s> CMD from sc2dalib__Trkheat\n",dateTime);
  SdpPuti("IN_SEQUENCE",DA_MCE_SINGLEDATA,status);
  
  // read the heater setting to use that as a starting point
  sc2dalib_readmceVal(con,myInfo,mceinxPtr,heatVal, &myInfo->heatSlp.refHeat,1,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status, "heat_init_track: sc2dalib_readmceVal failed"); 
   return;
  }
  MsgOut(status,"heat_init_track: initHeat=%d",myInfo->heatSlp.refHeat);
  myInfo->trkNo=0;   
}

#define TOTALPIXEL (ROW_NUM-1)*COL_NUM
#define HEAT_STEP 20
/**
 * \fn void heat_update_track(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceInxpt, int *data, double *slope, 
 *  StatusType *status)
 *
 * \brief function
 *  compare the pixel data with reference Value, adjust heater accordingly.
 *  use heatslope obtained from HEAT_SLOPE
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param mceInxpt  struct mcexml_struct pointer
 * \param data     int pointer for the data         
 * \param slope    double pointer for heaterSlope         
 * \param status StatusType.  given and return
 *
 * Used to be sc2dalib_trkheatUpdate
 */
/*+ heat_update_track
 */
void heat_update_track
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
int                  *data,
double               *slope,
StatusType           *status
)
{
  char         heatCmd[FILE_LEN];
  int          pixelData[HEAT_TRACK_MAX], indexArray[HEAT_TRACK_MAX];
  int          pixel,*frameData, medianIndex;
  int          headNo=FRAMEHEADER_NUM;
  double       heatVal[HEAT_TRACK_MAX], meanHeat, medianHeat, sigmaHeat;
  char         readheatCmd[]="rb bc1 bias 1";
  double       tai;
  long         taiDays;
  int row,col,i;
  int heatDiff;
  static int first, lastHeat, messCount;
  static long heatFlag;
  static int   heat, fdbkRef[HEAT_TRACK_MAX], numActivePixels, activePixels[HEAT_TRACK_MAX], testData;
  static double activeSlopes[HEAT_TRACK_MAX];
  StatusType  localStatus;

  if (!StatusOkP(status)) return;

  // skip the frame status 
  frameData = data+headNo;
 
  // use the first data from MCE as the REFERENCE data  
  if (myInfo->trkNo == 1)
  { /* Clear the heater tracking has failed flag and set first to true */
    myInfo->heatTrackFailed = 0;
    first = 1;
    messCount = 0;

    /* Initialize the current time to be TAI */
    sc2headman_getTai(&tai, status);
    taiDays = (int)tai;
    myInfo->heatSlp.curTime = tai - taiDays;

    /* Only do this if the UPDATE_HEATER_TRACKING_REFERENCE_BIT has 
       been set in the HTTRACK_FLAG */

    SdpGeti("HTTRACK_FLAG", &heatFlag, status);

    if((heatFlag & UPDATE_HEATER_TRACKING_REFERENCE_BIT) != 0)
    {
      /* Always clear the UPDATE_HEATER_TRACKING_REFERENCE_BIT, 
        If the REMEMBER_VALUE_IN_DARK is set */
      if((heatFlag & REMEMBER_VALUE_IN_DARK) != 0)
      {
        heatFlag &= ~UPDATE_HEATER_TRACKING_REFERENCE_BIT;
        SdpPuti("HTTRACK_FLAG", heatFlag, status);
      }

      sc2dalib_readmceVal(con,myInfo,mceInxpt,readheatCmd, &heat,1,status);
      if (!StatusOkP(status)) 
      { 
        ErsRep(0,status, "heat_update_track: sc2dalib_readmceVal failed"); 
        return;
      }
      myInfo->heatSlp.refHeat=heat;

      //get the pixel data,it is sq1-fdbk Val
      numActivePixels = 0;
      for (pixel=0; pixel<TOTALPIXEL; pixel++)
      {
        if( slope[pixel] != 0 )
        {
          testData = *(frameData + pixel);
          if(testData < 0 ) testData = testData * -1;
          if(testData > 10000)
          {
            fdbkRef[numActivePixels]= *(frameData + pixel);
            activePixels[numActivePixels]= pixel;
            activeSlopes[numActivePixels]= slope[pixel];
        
            col = pixel % COL_NUM;
            row = (floor)(pixel / COL_NUM);
            myInfo->parshmPtr->heat_track_values[numActivePixels][0] = col;
            myInfo->parshmPtr->heat_track_values[numActivePixels][1] = row;

            numActivePixels += 1;
            if(numActivePixels > (HEAT_TRACK_MAX-1)) numActivePixels = (HEAT_TRACK_MAX-1);
          }
        }
      }
      myInfo->parshmPtr->heat_track_num = numActivePixels;
    }
    /* If we are writing to a heater tracking file, print out the header */
    if ( (myInfo->heatSlp.option & 0x01 ) != 0)
    {
      fprintf(myInfo->fpStrchart,"# MultiPixel heater tracking Setup information:  TAI Day: %ld\n", taiDays);
      fprintf(myInfo->fpStrchart,"# Option: %d  maxOffset: %d  refHeat %d\n",
      myInfo->heatSlp.option, myInfo->heatSlp.maxOffset, myInfo->heatSlp.refHeat);

      for(i=0; i<numActivePixels; i++)
      {
        /* fprintf(myInfo->fpLog,"Pixel: %d fdbkRef: %d\n",activePixels[i], fdbkRef[i]); */
        fprintf(myInfo->fpStrchart,"# Pixel: %5d fdbkRef: %10d slope: %15.2f \n",
                activePixels[i], fdbkRef[i],activeSlopes[i]);
      }
      fprintf(myInfo->fpStrchart,"#  Day Fraction   DAC Setting  Slope    Next DAC    pixelNum  pixelData  pixelResult    Mean      Sigma\n");
    }
    lastHeat = heat;
  }   
  else /* It is not the first sample of this heater tracking event */
  {
    if((heatFlag & SIMULATE_HEATER_TRACKING_BIT) != 0)
    {
      return;
    }

    heat_pixel_track(myInfo, heat, frameData, fdbkRef, numActivePixels,
                     activePixels, activeSlopes, heatVal, first,
                     pixelData, status);
    if (!StatusOkP(status))
      return;
    first = 0;

    // If bit 1 of heatSlp.option is zero use the meanVal
    if ( (myInfo->heatSlp.option & 0x02 ) == 0)
    {
      sc2math_clipmean(3.0, numActivePixels, heatVal, &meanHeat, status);
      if (!StatusOkP(status)) 
      {
        ErsRep(0,status, "heat_update_track: sc2math_clipmean failed"); 
        return;
      }
      heat=(int)meanHeat;
    }
    else
    { // Use the median value
      /* Just to keep track of where the median is in all the other arrays */
      for(i=0; i<numActivePixels; i++)
      {
        indexArray[i] = i;
      }

      /* Get the stats for this set of proposed heater settings */
      sc2dalib_getStats(heatVal, indexArray, &medianHeat, &medianIndex, numActivePixels, &meanHeat, &sigmaHeat);
      heat = (int)medianHeat;
    }

    heatDiff = heat - lastHeat;
    if(heatDiff < 0) heatDiff = heatDiff * (-1);
    if(heatDiff > HEAT_STEP)
    {
      if((heat - lastHeat) > HEAT_STEP) 
      {
        heat = lastHeat + HEAT_STEP;
      }
      else if((heat - lastHeat) < HEAT_STEP) 
      {
        heat = lastHeat - HEAT_STEP;
      }
    }
    lastHeat = heat;

    /* If we are writing a stripchart file and using median put the stats in the log file */
    if ( (myInfo->heatSlp.option & 0x03 ) == 0x03)
    {
      fprintf(myInfo->fpStrchart,"%11d %11d %11d %11.4f %11.4f %11.4f\n",
              heat, activePixels[medianIndex], pixelData[medianIndex],
              medianHeat, meanHeat, sigmaHeat);
    }
    else if((myInfo->heatSlp.option & 0x03 ) == 0x01)
    {
      fprintf(myInfo->fpStrchart,"%11d %11.4f\n", heat, meanHeat);
    }

    // Out of bounds checking
    if ( heat > 65535 )
    {
      heat=65535;
      myInfo->heatTrackFailed = 1;
      if(messCount < 3)
      {
        localStatus = DITS__APP_ERROR;
        ErsRep(0, &localStatus,"_trkheatpixelUpdate heater tracking heatVal has hit upper rail");
      }
      messCount++;
    }
    else if( heat < 0 ) // Allow zero because it might happen on the blackbody source
    {
      heat=0;
      myInfo->heatTrackFailed = 1;
      if(messCount < 3)
      {
        localStatus = DITS__APP_ERROR;
        ErsRep(0, &localStatus,"_trkheatpixelUpdate heater tracking heatVal is below lower rail");
      }
      messCount++;
    }

    sprintf(heatCmd, "wb bc1 bias %d", heat);
    sc2dalib_setmceVal(con,myInfo,mceInxpt,heatCmd,status);
    /* fprintf(myInfo->fpLog,"_trkheatUpdate %s\n",heatCmd); */ 
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"heat_update_track: sc2dalib_setmceVal %s failed",
             heatCmd); 
      return;
    }
  }
}

#define ONE_HUNDRETH_SEC_MJD 0.000000115
/**
 * \fn void heat_pixel_track(dasInfoStruct_t *myInfo,
 *  int heat, int *data, double *slope, int *fdbkRef, double *heatVal,
 *  int first, StatusType *status)
 *
 * \brief function
 *  compare the pixel data with reference Value, adjust heater accordingly.
 *  use heatslope obtained from HEAT_SLOPE or READHEAT_SLOPE
 *
 * \param myInfo   dasInfoStruct_t     pointer
 * \param heat            int          previous heater setting
 * \param data            int          pointer for the data at the real data position          
 * \param fbdkRef         int          pointer for initial sq1 fdbk (what we try to keep the pixel reading)
 * \param numActivePixels int          Number of active pixels
 * \param activePixels    int          The pixel numbers for each active pixel
 * \param slope           double       The pre-measured heater slopes for each active pixel
 * \param heatVal         double       pointer for heat val at selected pixel (what each pixel thinks the next heater setting should be)
 * \param first           int          True if this is first time routine has been called (this event)
 * \param pixelData       int          The actual data read from each active pixel
 * \param status          StatusType.  given and return
 *
 * Used to be: sc2dalib_trkheatpixelUpdate
 */
/*+ heat_pixel_track
 */
void heat_pixel_track
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
)
{
  int          diffVal,pixel,i,diff;
  double       adjVal, adjPwr, adjCurrent, tempD;
  double       fullDC=65535, pixelResist=2, current=24.8, c;
  static       int  lastHeat, heatSlope;
  static       double  lastValues[HEAT_TRACK_MAX];
  int caught[HEAT_TRACK_MAX], numcaught=0;
  StatusType           localStatus;

  if (!StatusOkP(status)) return;

  if(first)
    {
      lastHeat=heat;
      for(i=0; i<HEAT_TRACK_MAX; i++) 
  lastValues[i] = heat;
    }

  /* full range of heater is 24.8 microamps, full range of DtoA is 65535 
     so c is the microamps per count */

  c = current/fullDC;

  if((heat < 0) || (heat > fullDC))
    {
      myInfo->heatTrackFailed = 1;
      localStatus = DITS__APP_ERROR;
      ErsRep(0,&localStatus,"Invalid heater reading passed to _trkheatpixelUpdate  heat: %d",heat);
      return;
    }

  heatSlope = heat - lastHeat;
  lastHeat = heat;

  if ( (myInfo->heatSlp.option & 0x01 ) != 0)
    fprintf(myInfo->fpStrchart,"%15.10f %11d  %6d", 
      myInfo->heatSlp.curTime, heat, heatSlope);

  myInfo->heatSlp.curTime += ONE_HUNDRETH_SEC_MJD;
 
  //calculate the heater Val for each active pixel
  for (i=0; i<numActivePixels; i++)
  {
    caught[i] = 0;
    pixel = activePixels[i];

    pixelData[i] = *(data+pixel);
    diffVal= pixelData[i] - fdbkRef[i];

    /* "diffVal" is the difference of what we would like the SQ1FDBK to be and what it 
       actually is. It is in SQ1FDBK units which are linearly proportional to power
       received by the pixel. "slope" is in SQ1FDBK units per microwatt, so
       adjPwr is simply the microwatts we need to change the power by */
  
    adjPwr= (double)diffVal/slope[i];

    /* The actual heater current is:
           heat * c
       because "heat" is in DtoA counts and c is in microamps per DtoA count. 
       Then the power that is going to the pixel due to that current is:
          pixelResist * (heat * c)^2
       We want to change that power by "adjPwr", so new power we want going to
       the pixel from the heater is:
         pixelResist * (heat * c)^2 - adjPwr
       But that power is caused by the current through "pixelResist"
          P = adjCurrent^2 * pixelResist
       If we solve for "adjCurrent" we get:
          adjCurrent^2 * pixelResist = pixelResist * (heat * c)^2 - adjPwr
          adjCurrent^2 = (heat * c)^2 - adjPwr/pixelResist
          adjCurrent = sqrt((heat * c)^2 - adjPwr/pixelResist) */

    tempD = pow((double)heat * c, 2) - (adjPwr/pixelResist);
    if(tempD > 0.0)
      {
  adjCurrent = sqrt( tempD );
      }
    else
      {
  adjCurrent = lastValues[i] * c; /* Better to guess last setting if we are having number issues */
      }

    /* Since adjCurrent is in microamps and c is in microamps per DtoA count
       We can convert to to new heater setting with "adjCurrent/c" */

    adjVal=adjCurrent/c;
    heatVal[i]=adjVal;

    /* If the requested change is too much, only change 
       by as much as we changed last time (eliminates noise spikes) */

    diff =  heatVal[i] - heat;
    if(diff < 0) diff = (-1)*diff;
    if(diff > myInfo->heatSlp.maxOffset) 
      {
  heatVal[i]=lastValues[i] + heatSlope;
  caught[i] = 1;
  numcaught++;
      }
    lastValues[i] = heatVal[i];

    /* Send to the dhtask to write into the data file */
    myInfo->parshmPtr->heat_track_values[i][2] =  heatVal[i];

    /* if we are writing a stripchart file then log this value
       if ( ((myInfo->heatSlp.option & 0x01 )!= 0) && (i < 25))
       fprintf(myInfo->fpStrchart," %11d %11.4f",pixelData[i],heatVal[i]); */

  } // end of for (i=0; i<numActivePixels; i++)

  /* Just do this stuff if we are writing to a stripchart file */

      if ((myInfo->heatSlp.option & 0x01 ) != 0)
  {
    if(numcaught > 0)
      {
        fprintf(myInfo->fpStrchart,"  #Caught: ");

        for (i=0; i<numActivePixels; i++)
    {
      if(caught[i] != 0) fprintf(myInfo->fpStrchart," %2d",i);
    }
      }
  }

  stripchFlag++;
}

/**
 * \fn void heat_step_current(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  struct mcexml_struct  *mceInxpt, StatusType *status)
 *
 * \brief function
 *  Check to see if it is a heater sawtooth or a heater ramp
 *  Based on which (half or quarter) cycle we are in change the heater current 
 *  D/A setting by +/- stairHeightCnts
 *
 *
 * \param con          SDSU context structure
 * \param myInfo       dasInfo structure pointer
 * \param mceInxpt     mcexml_struct pointer for parameter lookup table   
 * \param status       StatusType     
 *
 * Used to be sc2dalib_stepHeaterCurrent
 */
/*+ heat_step_current
*/
void heat_step_current
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,   
StatusType            *status
)
{

  char       heatCmd[FILE_LEN];
  int        currentSet;

  if (*status != STATUS__OK) return;

  currentSet = myInfo->stairPresentValue;

  myInfo->stairNumCount += 1;
  stairCounter += 1;

  /* Do not do the very last step of the SEQUENCE */
  if(stairCounter == totalStairCount) return;

  /* Is it a Heater sawtooth (fast flat field) */
  if( (myInfo->drcontrol & FLATFIELD_BIT) != 0)
  {
    if((myInfo->stairQCycleCount == 0) || (myInfo->stairQCycleCount == 3))
    {
      currentSet += myInfo->stairHeightCnts;
    }
    else
    {
      currentSet -= myInfo->stairHeightCnts;
    }
    if(myInfo->stairNumCount >=  myInfo->stairNum)
    {
      myInfo->stairNumCount = 0;
      myInfo->stairQCycleCount += 1;
      if(myInfo->stairQCycleCount >= 4)
      {
        myInfo->stairQCycleCount = 0;
      }
    }
  }

  /* Or is it a Heater ramp */
  if( (myInfo->drcontrol & HEATRAMP_BIT) != 0)
  {
    if(myInfo->stairHCycleCount == 0)
    {
      currentSet -= myInfo->stairHeightCnts;
    }
    else
    {
      currentSet += myInfo->stairHeightCnts;
    }
    if(myInfo->stairNumCount >=  myInfo->stairNum)
    {
      myInfo->stairNumCount = 0;
      myInfo->stairHCycleCount += 1;
      if(myInfo->stairHCycleCount == 1)
      {
        currentSet += myInfo->stairHeightCnts;
      }
      else
      {
        currentSet -= myInfo->stairHeightCnts;
        myInfo->stairHCycleCount = 0;
      }
    }
  }

  /* Now set the heater current to this value */

  myInfo->stairPresentValue = currentSet;
  sprintf(heatCmd, "wb bc1 bias %d", currentSet);
  /* fprintf(myInfo->fpLog,"_stepHeaterCurrent %s\n",heatCmd); */

  // send cmd and get reply
  sc2dalib_setmceVal(con, myInfo, mceInxpt, heatCmd, status); 
  if ( !StatusOkP(status) )
  {
    /* ErsRep(0,status,"heat_step_current: sc2dalib_setmceval(1) %s failed",heatCmd); */ 
    return;
  }  
}

/**
 * \fn void heat_step_TES(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  struct mcexml_struct  *mceInxpt, StatusType *status)
 *
 * \brief function
 *  Check to see if it is a bias sawtooth or a bias ramp
 *  Based on which (half or quarter) cycle we are in change the bias  
 *  D/A setting by +/- stairHeightCnts
 *
 *
 * \param con          SDSU context structure
 * \param myInfo       dasInfo structure pointer
 * \param mceInxpt     mcexml_struct pointer for parameter lookup table   
 * \param status       StatusType     
 *
 * Used to be: sc2dalib_stepTESBias
 */
/*+ heat_step_TES
*/
void heat_step_TES
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,   
StatusType            *status
)
{

  char       biasCmd[FILE_LEN];
  int        currentSet;

  if (*status != STATUS__OK) return;

  currentSet = myInfo->stairPresentValue;

  myInfo->stairNumCount += 1;
  stairCounter += 1;

  /* Do not do the very last step of the SEQUENCE */
  if(stairCounter == totalStairCount) return;

  /* Is it a Bias sawtooth */
  if( (myInfo->drcontrol & BIASSAW_BIT) != 0)
  {
    if((myInfo->stairQCycleCount == 0) || (myInfo->stairQCycleCount == 3))
    {
      currentSet += myInfo->stairHeightCnts;
    }
    else
    {
      currentSet -= myInfo->stairHeightCnts;
    }
    if(myInfo->stairNumCount >=  myInfo->stairNum)
    {
      myInfo->stairNumCount = 0;
      myInfo->stairQCycleCount += 1;
      if(myInfo->stairQCycleCount >= 4)
      {
        myInfo->stairQCycleCount = 0;
      }
    }
  }

  /* Or is it a Bias ramp */
  if( (myInfo->drcontrol & BIASRAMP_BIT) != 0)
  {
    if(myInfo->stairHCycleCount == 0)
    {
      currentSet -= myInfo->stairHeightCnts;
    }
    else
    {
      currentSet += myInfo->stairHeightCnts;
    }
    if(myInfo->stairNumCount >=  myInfo->stairNum)
    {
      myInfo->stairNumCount = 0;
      myInfo->stairHCycleCount += 1;
      if(myInfo->stairHCycleCount == 1)
      {
        currentSet += myInfo->stairHeightCnts;
      }
      else
      {
        currentSet -= myInfo->stairHeightCnts;
        myInfo->stairHCycleCount = 0;
      }
    }
  }

  /* Now set the bias to this value */
  myInfo->stairPresentValue = currentSet;
  sprintf(biasCmd, "wb bc2 bias %d", currentSet);

  /* fprintf(myInfo->fpLog,"_stepTESBias BIASSAW %s count %d\n",biasCmd,stairCounter); */

  // send cmd and get reply
  sc2dalib_setmceVal(con, myInfo, mceInxpt, biasCmd, status); 
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"heat_step_TES: sc2dalib_setmceval(1) %s failed",biasCmd); 
    return;
  } 
}

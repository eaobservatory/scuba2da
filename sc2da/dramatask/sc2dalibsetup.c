/** 
 * \file sc2dalibsetup.c
 *
 * \brief collection of sub-functions for reading parameters
 *        for servo mostly 
 *
 * \author Xiaofeng Gao, UKATC (xg@roe.ac.uk)
 *
 * \version 
 *
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
 *  $Log: sc2dalibsetup.c,v $
 *  Revision 1.31  2011/12/20 23:23:14  cwalther
 *  Implemented a style of heater tracking that ends after a fixed number of iterations
 *
 *  Revision 1.30  2011/09/02 19:16:21  bgorges
 *  Standardized entry checks on all fuctions/actions that should quietly return if entry status is bad.
 *
 *  Revision 1.29  2011/06/02 21:31:22  bgorges
 *  Checking all of fopen and fopen64 to make sure a file is actually opened.
 *
 *  Revision 1.28  2011/03/01 23:47:28  cwalther
 *  Added test for maximum gaini, set to zero if greater than maximum
 *
 *  Revision 1.27  2010/12/04 00:29:36  cwalther
 *  Removing changing select_clk each time as it was crashing the MCE
 *
 *  Revision 1.26  2010/10/11 23:55:12  cwalther
 *  Changed full scale heater current to 24.8 instead of 20 microamps
 *
 *  Revision 1.25  2010/08/03 22:16:10  cwalther
 *  Changes to make sc2_setup work with SC2SCRATCH
 *
 *  Revision 1.24  2010/03/04 01:32:54  cwalther
 *  Changes for heater tracking memory and setting back to defaults in the dark
 *
 *  Revision 1.1.1.1  2007/05/16 08:26:59  dkelly
 *  first insertion

 */ 
#include <math.h>
#include "sc2math.h"

//#define PRINT_FILTED     // fprint(myInfo->fpBatch) the filted data for row=/whichCh 
//#define PRINT_CONVOLUTION // fprint the convolution data for row=/whichCh

//#define DO_TEST_COL   // print out data[row] in sc2dalibsetup_sq2fbOutlier
//#define TEST_COL  1

static  char   errmsg[FILE_LEN];
 
/**
 * \fn void sc2dalibsetup_applychFilter(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     double *chData, int whichCh, int row, int whichRow, , int dataLen, 
 *     int *filtedchData, StatusType *status)
 *
 * \brief function
 *  apply filter to chData and save result to filtedchData 
 * 
 * \param  myInfo   dasInfo structure pointer  
 * \param  setup    ARRAYSET structure pointer 
 * \param  chData   double pointer for all channel data/heat or bias lvl 
 * \param  whichCh  int which channel is looked at
 * \param  row      int row
 * \param  whichRow int  which row is looked at
 * \param  dataLen  int  data length for each channel/column
 * \param  filtedchData int pointer for filteddata /heat or bias level 
 * \param  status   StatusType pointer  
 *
 * note: (N-1)/2 delay
 * 
 */
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
)
{
  int        i,feedbk;
  int        M=(FILTER_ORDER-1)/2;
  int        filtedStart, filtedEnd;
  double     xData,filtedData;

  if (!StatusOkP(status)) return;

  // http://www.cppreference.com/stdstring/memset.html
  // initial all val =0
  memset( setup->convolPtr, '\0', sizeof(double)*FILTER_ORDER );

  filtedStart=2*M;
  filtedEnd= dataLen + 2*M;

  for (i=0; i< filtedEnd; i++)
  {
     if( i < M ) /* 0 M-1 is additional points for convol */
       // use second point not the fisrt one as it appears to follow
       // last heatLvl
       xData=chData[whichCh+COL_NUM*1]; 
     else if( i < (dataLen+M) )
       xData=chData[whichCh+COL_NUM*(i-M)] ;     
     else  /* dataLen+M is additional points for convol, use the last-2 */
       xData=chData[whichCh+COL_NUM*(dataLen-2)] ;

     sc2dalibsetup_chlinearConv(myInfo,setup,xData, &filtedData, row,whichRow, whichCh,
                       i, filtedStart, filtedEnd,status);

     // take (N-1)/2 delay off plus the pending ones
     if( i>=filtedStart &&  i<filtedEnd )
     {
       feedbk=whichCh+ COL_NUM*(i-filtedStart);
       filtedchData[feedbk]=(int)(filtedData+0.5);
       #ifdef  PRINT_FILTED
         if(row==whichRow && whichCh==1)
           fprintf(myInfo->fpBatch,"returned [%6.4f][%6d]\n",filtedData, filtedchData[feedbk]);
       #endif
     }
  }
}


/**
 * \fn void sc2dalibsetup_chlinearConv(dasInfoStruct_t *myInfo, ARRAYSET *setup, 
 *     double xData, double *filteData, int row,int whichRow, int whichCh,
 *     int ith, int filtedBeg, int filtedEnd, StatusType *status)
 *
 * \brief function
 *  apply linear convolution and save result to filteData 
 * 
 * \param myInfo   dasInfo structure pointer  
 * \param setup    ARRAYSET structure pointer  
 * \param xData     double  pointer
 * \param filteData double pointer 
 * \param row       int  in seq for rows
 * \param whichRow  int which row is looked at
 * \param whichCh   int
 * \param ith       int the ith filted data,
 * \param filtedBeg int beginning of the filted data
 * \param filtedEnd int end of the filted data 
 * \param  status   StatusType pointer  
 *
 * direct convolution to filte data, shiting input
 * y(n)=convolution ( k=0, N-1) h(k)x(n-k)
 *
 * note: (N-1)/2 delay
 */
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
)
{
  int    n;
  double temp=0;

  if (*status != STATUS__OK) return;

  setup->convolPtr[0]=xData;
  // don't need to do the convolution if ith<filtedBeg
  if( ith >= filtedBeg && ith <filtedEnd )
  {
    for( n=0; n<FILTER_ORDER; n++)
    {
      temp +=setup->impulsePtr[n]*setup->convolPtr[n];
      #ifdef  PRINT_CONVOLUTION
        if(row==whichRow && whichCh==1)
        {  
          fprintf(myInfo->fpBatch,"[%6.4f]*[%6.4f] ",
               setup->impulsePtr[n],setup->convolPtr[n]);
        }
      #endif
    }
    *filteData= temp;

    #ifdef  PRINT_CONVOLUTION
      if(row==whichRow && whichCh==1)
        fprintf(myInfo->fpBatch,"filtedData[%6.4f]",temp);
    #endif
  }

  // shift for next data, from high inx to low inx
  for( n=FILTER_ORDER-1;  n >0; n--)
    setup->convolPtr[n]=setup->convolPtr[n-1];
}


/**
 * \fn void sc2dalibsetup_sq2fbeachcolOutlier(ARRAYSET *setup, OPT_STRUCT *sq2fbHist,
 *  int *data,  int col, StatusType *status)
 *
 * \brief find outliers before optimal routine
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  sq2fbHist  OPT_STRUCT structure pointer
 * \param  data     int array for sq2fdbk each row
 * \param  col     int 
 * \param status  StatusType pointer.  given and return
 *
 *  get minimum and maximum value from the sq2fb data.  
 *  use maximum value to malloc a buffer for sorted data.   
 *  calculate the histogram to find out the outlier
 *
 */
/*+ sc2dalibsetup_sq2fbeachcolOutlier
*/
void sc2dalibsetup_sq2fbeachcolOutlier
(
ARRAYSET   *setup,
OPT_STRUCT *sq2fbHist,
int        *data,
int        col,
StatusType *status
)
{
  // finds the smallest and largest members of a dataset 
  int min = data[0];
  int max = data[0];
  int binInx[ROW_NUM];
  int i, period, bin, whichBin=0,pixel;
  int peakVal,halfFlux;

  if (*status != STATUS__OK) return;

  // if data[i]==0, do not count, data[i]>=0
  for (i = 0; i < ROW_NUM; i++)
  {
    if( min==0)
      min=data[i];    
    else if ( data[i] !=0 && data[i] < min)
      min = data[i];
    if (data[i] > max)     
      max = data[i];
  }
  sq2fbHist->minOut[col] = min ;
  sq2fbHist->maxOut[col] = max ;

  if (min==0 && max==0)
   return;

  period=(max-min);
  bin=period/BINDIV;
  sq2fbHist->bin[col]=bin;
  
  if (min <0) 
  {  
    *status=DITS__APP_ERROR;
    sprintf(errmsg, "sc2dalibsetup_sq2fbeachcolOutliers: MUST NOT LET min (%d) <0 \n",min ); 
    ErsRep(0,status,errmsg);
    return; 
  }

  if (bin ==0) 
  {  
    // *status=DITS__APP_ERROR;
    sprintf(errmsg, "sc2dalibsetup_sq2fbeachcolOutliers: bin 0 max(%d) min(%d) \n",
              max,min ); 
    MsgOut(status,errmsg);
    return; 
  }

  for (i=0; i< ROW_NUM; i++)
  { 
    // whichBin start from 0, so, it is ok
    if (data[i]!=0)
    {
      whichBin=(data[i]-min)/bin;
      binInx[i]=whichBin;
      sq2fbHist->freqPtr[whichBin] ++;  
    }
  }

  // search for the whole freq
  max=sq2fbHist->freqPtr[0];
  for (i = 0; i <=BINDIV; i++)
  {
    if (sq2fbHist->freqPtr[i] >= max)    
    {
      max = sq2fbHist->freqPtr[i];
      whichBin=i;
    }
  }
  sq2fbHist->peakVal[col]=whichBin*bin + min + bin/2;
  sq2fbHist->peakInx[col]=whichBin;
  halfFlux=setup->fluxPeriod[col]/2;
  peakVal=sq2fbHist->peakVal[col];
  
  // take out anyone who is 15 bins away
  // consider if we can fold both end to each other
  // fold bottom end to top, this is OK,

  if (whichBin > BINTOP)
  {
    for (i=0; i<ROW_NUM; i++)
    { 
      pixel=i*COL_NUM +col;
      
      if ( binInx[i] < BINBOT )
      {
        data[i]+= setup->fluxPeriod[col];
	sq2fbHist->changeInx[pixel]=1;
      }
      else if ( abs(binInx[i]-whichBin) > SQ2FDBK_WIDTH )
        sc2dalibsetup_sq2fbhalfFlux(setup,data,i,peakVal,halfFlux,
	         &sq2fbHist->changeInx[pixel],status);
    }
  }
  else if (whichBin < BINBOT)
  {
    // fold top to bottom end, unless we take negative
    // this is not, here we fold the negative to positive
    for (i=0; i<ROW_NUM; i++)
    { 
      pixel=i*COL_NUM +col;
      if ( binInx[i] > BINTOP )
      {
        data[i] -= setup->fluxPeriod[col];
	sq2fbHist->changeInx[pixel]=-1;
      }
      else if ( abs(binInx[i]-whichBin) > SQ2FDBK_WIDTH )
        sc2dalibsetup_sq2fbhalfFlux(setup,data,i,peakVal,halfFlux,
	         &sq2fbHist->changeInx[pixel],status);
    }
  }
  else
  {
    for (i=0; i<ROW_NUM; i++)
    { 
      pixel=i*COL_NUM +col;
      if ( abs(binInx[i]-whichBin) > SQ2FDBK_WIDTH )
        sc2dalibsetup_sq2fbhalfFlux(setup,data,i,peakVal,halfFlux,
	         &sq2fbHist->changeInx[pixel],status);
    }
  }
  return;
}



/**
 * \fn void sc2dalibsetup_sq2fbhalfFlux(ARRAYSET *setup,int *data, int which, int peakVal,
 *   int halfFlux, int *changeInx, StatusType *status)
 *
 * \brief find those having half flux distance from the peak
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  data     int pointer
 * \param  which    int
 * \param  peakVal  int peak in histogram
 * \param  halfFlux int half fluxperiod
 * \param  changeInx int pointer  for recording the change
 * \param status  StatusType pointer.  given and return
 *
 */
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
)
{
  int   apart;

  if (*status != STATUS__OK) return;

  // consider fold half-flux:

  apart=abs(peakVal-data[which]);
  
  // THD_FLUX  1000 as initial threadhold for half-fluxperiod fold
  // if setup->slopSelect[12] set <0, then we don't do half-fold 
  if ( abs(apart-halfFlux) <setup->slopSelect[12] )
  {
    if( data[which] > peakVal )
    {
      data[which]-=halfFlux;
      *changeInx=-2;
    }
    else
    {
      data[which]+=halfFlux;
      *changeInx=2;
    }
  }
  else
    data[which]=0;
}     



/**
 * \fn void sc2dalibsetup_sq2fbOutlier(dasInfoStruct_t  *myInfo,   ARRAYSET *setup,
 *  char *filename, int *sq2feedback, int *sq2reffeedback, StatusType *status)
 *
 * \brief find outliers before optimal routine
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  filename    char point for the file name
 * \param  sq2feedback     int pointer
 * \param  sq2reffeedback  int pointer
 * \param status  StatusType pointer.  given and return
 *
 *  get minimum and maximum value from the sq2fb data.  
 *  use maximum value to malloc a buffer for sorted data. 
 *  calculate the histogram to find out the outlier
 * 
 *
 */
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
)
{
  char  tmpfile[FILE_LEN];
  FILE  *fp, *fpFreq;
  int   row, col,pixel, i;
  int   data[ROW_NUM], badRow[ROW_NUM];
  static OPT_STRUCT  sq2fbHist;

  if (*status != STATUS__OK) return;

  if((fp = fopen (filename, "w")) == NULL)
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalibsetup_sq2fbOutlier: Error- failed to open file %s", filename); 
      return;
    }
  sprintf(tmpfile,"%s-freq",filename);
  if((fpFreq = fopen(tmpfile, "w")) == NULL)
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalibsetup_sq2fbOutlier 2: Error- failed to open file %s", tmpfile); 
      return;
    }

  fprintf(fp,"sc2dalibsetup_sq2fbOutliers: Histogram for each col\n");
  fprintf(fpFreq,"# Histogram for each col, ");
  fprintf(fpFreq,"first line is sq2fbVal, secondline is the frequency\n");

  for ( col=0; col<COL_NUM; col++ ) 
  {
    memset(sq2fbHist.freqPtr,'\0', sizeof(int)*(BINDIV+1)); 
    sq2fbHist.peakVal[col]=sq2fbHist.peakInx[col]=0;
    sq2fbHist.maxOut[col]=sq2fbHist.minOut[col]=0;
    sq2fbHist.bin[col]=0; 
    #ifdef DO_TEST_COL
    {
      if (col==TEST_COL)
        printf("sq2fbdkcol_%d :\n",col);
    }
    #endif 
    for ( row=0; row<ROW_NUM; row++ )
    {
      badRow[row]=0;
      pixel = row * COL_NUM +col;
      if ( sq2feedback[pixel] ==VAL__BADI )
      {
        data[row]=0; badRow[row]=1; 
      }
      else 
       data[row]=sq2feedback[pixel];
       
      sq2fbHist.changeInx[pixel]=0; 
      #ifdef DO_TEST_COL
      {
        if (col==TEST_COL)
          printf("%d ",data[row]);
      }
      #endif 
    }

    #ifdef DO_TEST_COL
    {
      if (col==TEST_COL)
       printf(" pixel=(%d), now call eachcolOutlier\n",pixel);
    }
    #endif 
    
    // if dead col, just return, now data is changed after the call
    sc2dalibsetup_sq2fbeachcolOutlier(setup,&sq2fbHist,data,col,status);
    if ( *status != STATUS__OK ) 
    {
      printf(
      "sc2dalibsetup_sq2fbOutlier: col=%d, failed in sc2dalibsetup_sq2fbeachcolOutlier\n",col);
      fclose(fp);  fclose(fpFreq);
      return ;
    }

    fprintf(fp,
    "col_%3d  bin[%7d]  pealVal[~%7d] peakInx[%4d] min[%7d] max[%7d]\n",
        col,sq2fbHist.bin[col],sq2fbHist.peakVal[col],sq2fbHist.peakInx[col],
	sq2fbHist.minOut[col],sq2fbHist.maxOut[col]);
    // take out the outliers
    for ( row=0; row<ROW_NUM; row++ )
    {
      pixel = row * COL_NUM + col;
      if( sq2fbHist.bin[col]==0 )
      {
        if (sq2fbHist.maxOut[col]==0 )
        {
          fprintf(fp,"row_%3d [%7d] is outlier, set BAD \n",
	        row,sq2feedback[pixel]);
          sq2feedback[pixel]=VAL__BADI;
        }
        else        
        {
          fprintf(fp,"row_%3d [%7d] are all the same, keep  \n",
	        row,sq2feedback[pixel]);
        }
      }
      else if ( data[row]==0  &&  badRow[row]==1)
      {
        fprintf(fp,"row_%3d [%7d] is outlier, set BAD \n",
	        row,sq2feedback[pixel]);
        sq2feedback[pixel]=VAL__BADI;        
      }
      else if (sq2fbHist.changeInx[pixel]!=0)
      {
        if (abs(sq2fbHist.changeInx[pixel])==1)
          fprintf(fp,"row_%3d is folded from %7d to %7d  Flux[%d] \n",
	        row,sq2feedback[pixel],data[row], setup->fluxPeriod[col]);
        else 
          fprintf(fp,"row_%3d is half folded from %7d to %7d  Flux[%d] \n",
	        row,sq2feedback[pixel],data[row], setup->fluxPeriod[col]);
  
        sq2feedback[pixel]=data[row];
        // change ref too
        if (sq2fbHist.changeInx[pixel]==1)
          sq2reffeedback[pixel]+=setup->fluxPeriod[col];
        else if (sq2fbHist.changeInx[pixel]==-1)
          sq2reffeedback[pixel]-=setup->fluxPeriod[col];
        else if (sq2fbHist.changeInx[pixel]==2)
          sq2reffeedback[pixel]+=setup->fluxPeriod[col]/2;
        else if (sq2fbHist.changeInx[pixel]==-2)
          sq2reffeedback[pixel]-=setup->fluxPeriod[col]/2;
      }      
    }
    fprintf(fpFreq,"col_%02d ",col);
    for (i=0;i<=BINDIV;i++)
     fprintf(fpFreq,"%07d ",(sq2fbHist.bin[col]*(i+1)+sq2fbHist.minOut[col]) );
    fprintf(fpFreq,"\n");

    fprintf(fpFreq,"col_%02d ",col);
    for (i=0;i<=BINDIV;i++)
     fprintf(fpFreq,"%07d ",sq2fbHist.freqPtr[i]);
    fprintf(fpFreq,"\n");
  }
  fclose(fp);
  fclose(fpFreq);
}



/**
 * \fn void sc2dalibsetup_writepixelInx(dasInfoStruct_t  *myInfo,   
 *  ARRAYSET *setup, char *filename, void *array, int flag,StatusType *status)
 *
 * \brief function:
 *  write out pixel Inx map
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  filename    char point for the base file name
 * \param  array      void pointer for the array
 * \param  flag       int  flag=0: bef opt, =1: aft opt  
 * \param status  StatusType pointer.  given and return
 *
 */
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
)
{
  int   pixel,row,col,totalPixel;
  int   goodOne, badOne; 
  int    *intData,data;
  FILE  *fp;
  char  tmp[FILE_LEN];

  if (*status != STATUS__OK) return;

  intData=NULL;

  goodOne=badOne=0;

  if(flag==1)
    sprintf(tmp,"%s-inx-aftopt",filename);
  else 
    sprintf(tmp,"%s-inx-befopt",filename);

  if((fp = fopen ( tmp, "w" )) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalibsetup_writepixelInx: Error- failed to open file %s", tmp); 
      return;
    }

  intData=(int*)array;

  if(flag==1)
  {
    for ( row=0; row<ROW_NUM; row++ )
    {
      for ( col=0; col<COL_NUM; col++ ) 
      {
        pixel= row * COL_NUM + col;
        data=intData[pixel];
        // bad ival gain >  maxGainI
        if ( data ==0 ||  abs(data)> myInfo->maxGainI )
        {
          fprintf ( fp, "x ");
          badOne++;
	  /* The gain is greater than maximum allowable so set it to zero */
	  intData[pixel] = 0;
        } 
        else
        {
           fprintf ( fp, "1 ");
           goodOne++;
        }
      }
      fprintf ( fp, "\n");   
    }
  }
  else
  {
#define SQUAL_BAD  1   // from sc2sqopt_getSQ1scales
    for ( row=0; row<ROW_NUM; row++ )
    {
      for ( col=0; col<COL_NUM; col++ ) 
      {
        pixel= row * COL_NUM + col;
        if ( intData[pixel] ==SQUAL_BAD )
        {
          fprintf ( fp, "x ");
          badOne++;
        } 
        else
        {
           fprintf ( fp, "1 ");
           goodOne++;
        }
      }
      fprintf ( fp, "\n");   
    }
  }
  totalPixel=COL_NUM*ROW_NUM;
  fprintf ( fp, "\nTotal goodOnes= %d ratio=%f\n",goodOne,(float)goodOne/totalPixel);
  fprintf ( fp, "\nTotal  badOnes= %d ratio=%f\n",badOne,(float)badOne/totalPixel);
  fclose ( fp );   
}



/**
 * \fn void sc2dalibsetup_savesq1biasoptPts(FILE *fp, char *name,
 *  int *paramarray, int flag, StatusType *status)
 *
 * \brief function:
 *  save optimal points for the array
  *
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray int pointer for setup-parameter array
 * \param flag        0: for setupfile 1: for batch file
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_savesq1biasoptPts
*/
void sc2dalibsetup_savesq1biasoptPts
(
FILE       *fp,
char       *name,
int        *paramarray,
int        flag,
StatusType *status
)
{
  int   i, j,l, m,section=8;

  if (*status != STATUS__OK) return;
  
  if (flag==0)
  {
    fprintf(fp,"%s",name);
    for (j=0;j<ROW_NUM;j++)
    {
      // trap also for those >MAX_SQ1BIAS, asked Dan  21-09-07
      if(  paramarray[j]==VAL__BADI ||  paramarray[j] >= MAX_SQ1BIAS )
       fprintf(fp,"0 ");
      else 
       fprintf(fp,"%7d ", paramarray[j]);
    }
    fprintf(fp," \n");
  } 
  else
  {
    for (i=0;i<5;i++)
    {  
      l=i*section;
      m=l+section; 
      for (j=l;j<m;j++)
      {
        if(  paramarray[j]==VAL__BADI ||  paramarray[j] >= MAX_SQ1BIAS )
         fprintf(fp,"0 ");
        else
         fprintf(fp,"%d ", paramarray[j]);
      }
      fprintf(fp,"<\n");
    }
    if(  paramarray[40]==VAL__BADI ||  paramarray[j] >= MAX_SQ1BIAS )
      fprintf(fp,"0 ");
    else
      fprintf(fp,"%d \n", paramarray[40]);    
  }
}


/**
 * \fn void sc2dalibsetup_savesq2fboptPts(FILE *fp, char *name,
 *  int *paramarray, int flag, StatusType *status)
 *
 * \brief function:
 *  save optimal points for the array
  *
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray int pointer for setup-parameter array
 * \param flag        0: for setupfile 1: for batch file
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_savesq2fboptPts
*/
void sc2dalibsetup_savesq2fboptPts
(
FILE       *fp,
char       *name,
int        *paramarray,
int        flag,
StatusType *status
)
{
  int   ch,i;

  if (*status != STATUS__OK) return;

  if (flag==0)
  {
    for (i=0;i<4;i++)
    {
      fprintf(fp,"#RC%d\n",i+1);
      fprintf(fp,"%s",name);
      for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
      {
	if(  paramarray[ch]==VAL__BADI )
	  fprintf(fp,"      0 ");
	else
          fprintf(fp,"%7d ", paramarray[ch]);
      }
      fprintf(fp,"\n");
    }
  }
  else
  {
    for (i=0;i<4;i++)
    {
      for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
      {
	if(  paramarray[ch]==VAL__BADI )
	  fprintf(fp,"0 ");
	else
          fprintf(fp,"%d ", paramarray[ch]);
      }
      if (i<3)
        fprintf(fp,"<\n");
      else
        fprintf(fp,"\n");
    }
  }
}


/**
 * \fn void sc2dalibsetup_readheaterSetup(dasInfoStruct_t *myInfo,
 *  StatusType *status)
 *
 * \brief function:
 *  parse the setup file, set initial values for find heater slope
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param status  StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_readheaterSetup
*/
void sc2dalibsetup_readheaterSetup
(
dasInfoStruct_t *myInfo,
StatusType      *status
)
{
  char  delimiters[]= " =\n", msg1[180],msg2[180];
  char  setupfile[SERV_STR_LEN*4], *token;   
  int   noInx=0, expected; 
      
  if (!StatusOkP(status)) return;

  //sprintf(msg2,"sc2dalibsetup_readheaterSetup: setupfile=%s",myInfo->batchFile); 
  //MsgOut(status,msg2);

  while(  fgets(setupfile,SERV_STR_LEN*4,myInfo->fpBatch) !=NULL )
  {
    /* The end of the file contains XML we cannot parse */
    if(strncmp(setupfile,"<!---",5) == 0) break;

    // skip the comment and empty line
    if ( (strchr(setupfile,'#')==NULL) && (strcmp(setupfile,"\n")!=0) )
    {
      // the setupfile shall look like
      // # comments
      // startVal=10000
      // step=100
      // nStep=10
      token = strtok (setupfile, delimiters);
      sprintf(msg1,"%s=",token); 

      if (strcmp("startVal",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.stVal=atoi(token);
        sprintf(msg2,"%d",myInfo->heatSlp.stVal); 
        strcat (msg1,msg2);
        MsgOut(status,msg1);
        noInx++;
      }
      else if (strcmp("step",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.step=atoi(token);
        sprintf(msg2,"%d",myInfo->heatSlp.step); 
        strcat (msg1,msg2);
        MsgOut(status,msg1);
        noInx++;
      } 
      else if (strcmp("nStep",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.nStep=atoi(token);
        sprintf(msg2,"%d",myInfo->heatSlp.nStep); 
        strcat (msg1,msg2);
        MsgOut(status,msg1);
        noInx++;
      } 
      else if (strcmp("fileBase",token)==0)
      {
        token = strtok (NULL, delimiters);
        strcpy(myInfo->heatSlp.base,token);
        sprintf(msg2,"%s",myInfo->heatSlp.base); 
        strcat (msg1,msg2);
        MsgOut(status,msg1);
        noInx++;
      } 
      else if (strcmp("option",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.option=atoi(token);
        sprintf(msg2,"%d",myInfo->heatSlp.option); 
        strcat (msg1,msg2);
        //MsgOut(status,msg1);
        noInx++;
      }
      else if (strcmp("powerFlag",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.flag=atoi(token);
        sprintf(msg2,"%d",myInfo->heatSlp.flag); 
        strcat (msg1,msg2);
        //MsgOut(status,msg1);
        noInx++;
      } 
      else if (strcmp("pGain",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.pGain = atof(token);
        sprintf(msg2,"%f",myInfo->heatSlp.pGain); 
        strcat (msg1,msg2);
        //MsgOut(status,msg1);
        noInx++;
      } 
      else if (strcmp("iGain",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.iGain = atof(token);
        sprintf(msg2,"%f",myInfo->heatSlp.iGain); 
        strcat (msg1,msg2);
        //MsgOut(status,msg1);
        noInx++;
      } 
      else if (strcmp("heatSlope",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.slope=atof(token);
        sprintf(msg2,"%f",myInfo->heatSlp.slope); 
        strcat (msg1,msg2);
        MsgOut(status,msg1);
        noInx++;
      } 
      else if (strcmp("refFDBK",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.refFDBK=atoi(token);
        sprintf(msg2,"%d",myInfo->heatSlp.refFDBK); 
        strcat (msg1,msg2);
        MsgOut(status,msg1);
        noInx++;
      }
      else if (strcmp("refHeat",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.refHeat=atoi(token);
        sprintf(msg2,"%d",myInfo->heatSlp.refHeat); 
        strcat (msg1,msg2);
        //MsgOut(status,msg1);
        noInx++;
      } 
      else if (strcmp("maxOffset",token)==0)
      {
        token = strtok (NULL, delimiters);
        myInfo->heatSlp.maxOffset=atoi(token);
        sprintf(msg2,"%d",myInfo->heatSlp.maxOffset); 
        strcat (msg1,msg2);
        //MsgOut(status,msg1);
        noInx++;
      } 
      else
      {
        *status=DITS__APP_ERROR;
        sprintf(errmsg,
          "sc2dalibsetup_readheaterSetup: %s in %s failed to match predefined name",
            token,myInfo->batchFile);
        ErsRep(0,status,errmsg);
        return;
      }
    }
  }
  if (myInfo->actionFlag==HEATERSLOPE)
     expected=HEATERSLOPE_NO;      
  else if (myInfo->actionFlag==TRKHEATACTION)
     expected=TRKHEAT_NO ; 
  else
  {
    *status=DITS__APP_ERROR;
    sprintf(errmsg,"sc2dalibsetup_readheaterSetup: unexpected action");
    ErsRep(0,status,errmsg);
    return;
  }  
  if (noInx !=expected) 
  {
    *status=DITS__APP_ERROR;
    sprintf(errmsg, 
     "sc2dalibsetup_readheaterSetup: the myInfo paras(%d) != expectedNo(%d)",
     noInx,expected);
    ErsRep(0,status, errmsg);
    return ;
  }
  //MsgOut(status,"sc2dalibsetup_readheaterSetup:OK");
}


/**
 * \fn void sc2dalibsetup_readheaterData(dasInfoStruct_t *myInfo,
 *  int *frameNo, double *chData,int *pixelData, int* heaterMask,
 *  StatusType *status)
 *
 * \brief function 
 *  read heater data in and get mean value for each pixel
 *
 * \param myInfo   dasInfoStruct_t pointer
 * \param frameNo    int pointer
 * \param chData    double pointer  for the ch data
 * \param pixelData   int pointer for the pixel (row,col)
 * \param heaterMask  int pointer for the heaterMask, not used yet

 * \param status     StatusType     
 * 
 * use myInfo->fpOtheruse to read data
 */
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
)
{
  // for reading the heater data file ( fstatus is very long)
  #define MAXLINE 613      
  char   *fStatus, *token;
  char   line[MAXLINE];
  int    totalPixel, totalrow;
  int    row, ch, pixel,firstEnt=1;
  int    howmanyFrame; 
  char   delimiters[]= " \n";
  double data;

  if (*status != STATUS__OK) return;
 
  // initial to zero
  totalPixel=ROW_NUM*COL_NUM;
  for (pixel=0; pixel<totalPixel; pixel++)
      chData[pixel] =0;

  howmanyFrame=0;
  // include the chksum
  totalrow=ROW_NUM+1;

  while (1)
  {
    if (firstEnt==1)
    {
      // frame data has a single line for all header inf.
      fStatus = fgets (line, MAXLINE, myInfo->fpOtheruse);
      if ( fStatus == NULL)
      {
        *status=DITS__APP_ERROR;
        sprintf(errmsg,
           "fgets failed. searched for %d frames\n",howmanyFrame);
        ErsRep(0,status, errmsg);
        return;
      } 
      firstEnt=0;
    }
    else
    { 
      for (row=0; row< totalrow; row++)
      {
        fStatus = fgets (line, MAXLINE, myInfo->fpOtheruse);
        if ( fStatus == NULL)
        {
          printf ("_readData: searched for %d frames\n",howmanyFrame);
          printf ("_readData: either end of file or fgets failed \n");
          return ;
        }
	  // last one is the checksum 
        if (row < ROW_NUM)
        {  
          for (ch=0; ch<COL_NUM; ch++)
          {
            if (ch==0) 
              token = strtok (line, delimiters);
            else
              token = strtok (NULL, delimiters);

            data=atof(token);
            pixel=row*COL_NUM + ch;
   
            chData[pixel] +=data;
	    
            if ( row==myInfo->heatSlp.row && ch==myInfo->heatSlp.col  &&
		     myInfo->filesvFlag ==1 )
            {
	      fprintf(myInfo->fpStrchart,"%11d \n",(int)data);
	    } 
          }
        }
      }
      howmanyFrame ++;
      //check if end of the FILE, or read next frame status line
      fStatus = fgets (line, MAXLINE, myInfo->fpOtheruse);
      if ( fStatus == NULL)
      {
        break;
      }
    }
  }
  
  for (row=0;row<ROW_NUM;row++)
  {
    for (ch=0;ch<COL_NUM;ch++)
    {
      pixel=row*COL_NUM + ch;
      chData[pixel] /=howmanyFrame;
      if (row==myInfo->heatSlp.row )
      {
        fprintf(myInfo->fpData,"%11d",(int)chData[pixel]);
        if (ch==myInfo->heatSlp.col )
        {
          *pixelData=(int)chData[pixel];
          fprintf(myInfo->fpBatch,"%11d",*pixelData);
        }
      }
    }
  }
  fprintf(myInfo->fpData,"\n");
  fprintf(myInfo->fpBatch,"\n");
  *frameNo=howmanyFrame;
}


/**
 * \fn void sc2dalibsetup_servoreadsetupWrap(dasInfoStruct_t *myInfo,
 *  StatusType *status)
 *
 * \brief function 
 *  Opens the file called myInfo->batchFile
 *  This file should contain a list of other files listed like include=xxx
 *  This contents of each of these files is sequentially copied to 
 *  a file called tmp
 *  Changes the name in myInfo->batchFile to just tmp
 *
 *  expand all include=xxx in the setup file into a single file for 
 *  sc2dalibsetup_servoreadSetup
 *
 * \param myInfo   dasInfoStruct_t pointer
 * \param status     StatusType     
 *
 */
/*+ sc2dalibsetup_servoreadsetupWrap
*/
void sc2dalibsetup_servoreadsetupWrap
(
dasInfoStruct_t       *myInfo,
StatusType            *status
)
{
#define  MAX_SETUP_INC   8

  FILE *fp;   /* pointer for file*/
  char oneLine[FILE_LEN];
  char shellCmd[FILE_LEN+FILE_LEN],msg[FILE_LEN+FILE_LEN];
  char baseFile[MAX_SETUP_INC][FILE_LEN];
  char delimiters[]= " =\n", *token;
  int  i, noInc; 

  if (*status != STATUS__OK) return;

  if ( (fp=fopen(myInfo->batchFile,"r"))==NULL)
    {
      *status=DITS__APP_ERROR;
      sprintf(msg,"sc2dalibsetup_servoreadsetupWrap: failed to open %s",myInfo->batchFile);
      ErsRep(0,status,msg);
      return ;
    }
  // search for include string 
  noInc=0; 
  while(  fgets(oneLine,SERV_STR_LEN*4,fp) !=NULL )
    {
      // skip the comment and empty line
      if ( (strchr(oneLine,'#')==NULL) && (strcmp(oneLine,"\n")!=0) )
	{
	  // the setupfile shall look like
	  // # comments
	  // include=basefile
	  //..
	  token = strtok (oneLine, delimiters);
	  if (strcmp("include",token)==0)
	    {  
	      token = strtok (NULL, delimiters);
	      strcpy(baseFile[noInc],token);
	      sprintf(msg,"include= %s",baseFile[noInc]); 
	      //printf("_readsetupWrap(%d):%s\n",noInc,msg);
	      jitDebug(16,"_readsetupWrap:%s\n",msg);
	      //MsgOut(status,msg);
	      noInc++;
	    }
	}
    } 
  fclose(fp);
  if ( noInc > MAX_SETUP_INC  || noInc ==0)
  {
     *status=DITS__APP_ERROR;
     sprintf(msg,"sc2dalibsetup_servoreadsetupWrap: no_of_include(%d) > MAX_SETUP_INC(%d)",noInc, MAX_SETUP_INC);
     ErsRep(0,status,msg);
     return ;
  }

  // make tmp before appending
  sprintf(shellCmd,"cat %s > %s/tmp",baseFile[0], getenv ( "SC2SCRATCH" ));
  //printf("_readsetupWrap:Shellcmd(%s)\n",shellCmd);
  sprintf(msg,"shellcmd 0 =%s", shellCmd);
  jitDebug(16,"_readsetupWrap:%s\n",msg);
  // MsgOut(status,msg);
  if( system( shellCmd) !=0)
  {
     *status=DITS__APP_ERROR;
     return ;
  } 
  if (noInc >0 && noInc <= MAX_SETUP_INC)
  {  
    for (i=1;i<noInc;i++)
    {     
      sprintf(shellCmd,"cat %s >> %s/tmp",baseFile[i], getenv ( "SC2SCRATCH" ) );

      sprintf(msg,"shellcmd %d =%s", i, shellCmd);  
      // MsgOut(status, msg);

      if( system(shellCmd) !=0)
      {
        *status=DITS__APP_ERROR;
        return ;
      }
    }
  }
  // re-assign the batchFile to the newer single setupfile
  sprintf(myInfo->batchFile,"%s/tmp", getenv ( "SC2SCRATCH" ));
}


 /**
 * \fn void sc2dalibsetup_servoreadSetup(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *      StatusType *status)
 *
 * \brief function:
 *  parse the setup file, set initial values
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param status  StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_servoreadSetup
*/
void sc2dalibsetup_servoreadSetup
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char        delimiters[]= " =\n";
  char        setupfile[SERV_STR_LEN*4], *token;   
  int         paramNo=0,biasoptInx,cableadjInx,adjscaleInx;
  int         cableInx,gainInx,zfactInx,scaleInx,cableslopeInx;
  int         sq1gainscaleInx,fluxInx;
  int         initfbInx,biaslckInx,slopInx,maskInx,rowmaskInx;
  int         sq2fdbkoptInx,ssabiaslckInx,cableinitInx;
      
  if (!StatusOkP(status)) return;

  cableInx=gainInx=zfactInx=scaleInx=cableslopeInx=0;
  initfbInx=biaslckInx=slopInx=maskInx=rowmaskInx=0;
  biasoptInx=cableadjInx=adjscaleInx=sq1gainscaleInx=fluxInx=0;
  sq2fdbkoptInx=0;
  ssabiaslckInx=cableinitInx=0;
  
  setup->servo=NONESERVO;
  setup->totalRC2use=0;

  while(  fgets(setupfile,SERV_STR_LEN*4,myInfo->fpBatch) !=NULL )
  {
    // skip the comment and empty line
    if ( (strchr(setupfile,'#')==NULL) && (strcmp(setupfile,"\n")!=0) )
    {
      // the setupfile shall look like
      // # comments
      // biasStart=100
      // bStep=10
      // biasNo= 300
      token = strtok (setupfile, delimiters);

      // define OTHERS_N    4
      if (strcmp("whichServo",token)==0)
        sc2dalibsetup_servosetCard(delimiters,setup->servoName,&paramNo,status);
      else if (strcmp("whichRow",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->selRow,&paramNo,status);
      else if (strcmp("doServo",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->doServo,&paramNo,status);
      else if (strcmp("whichRC",token)==0)
      {
        sc2dalibsetup_servosetwhichRC(myInfo,delimiters,setup,&paramNo,status);
        if (!StatusOkP(status)) 
           return;
      }
      //#define CARDS_N     4
      else if (strcmp("safb_card",token)==0)
        sc2dalibsetup_servosetCard(delimiters,setup->safbCard,&paramNo,status);
      else if (strcmp("sq2fb_card",token)==0)
        sc2dalibsetup_servosetCard(delimiters,setup->sq2fbCard,&paramNo,status);
      else if (strcmp("sq2bias_card",token)==0)
        sc2dalibsetup_servosetCard(delimiters,setup->sq2biasCard,&paramNo,status);
      else if (strcmp("sq1bias_card",token)==0)
        sc2dalibsetup_servosetCard(delimiters,setup->sq1biasCard,&paramNo,status);

      //#define BIAS_N      3
      else if (strcmp("biasStart",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->minBIAS,&paramNo,status);
      else if (strcmp("bStep",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->stepBIAS,&paramNo,status);        
      else if (strcmp("biasNo",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->biasNo,&paramNo,status);
     
      //#define FEEDBK_N    3
      else if (strcmp("fbStart",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->minFDBK,&paramNo,status);
      else if (strcmp("fbStep",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->stepFDBK,&paramNo,status);
      else if (strcmp("fbNo",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->fdbkNo,&paramNo,status);

      //#define SLOPSELECT_N 4
      else if (strcmp("slopSelect",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->slopSelect,
                              "slopSelect",&slopInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }
      //#define COLMASK_N 4
      else if (strcmp("colMask",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->colMask,
                              "colMask",&maskInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }
      //#define ROWMASK_N 1
      else if (strcmp("rowMask",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->rowMask,
                              "rowMask",&rowmaskInx,&paramNo,ROW_NUM,status);
        if (!StatusOkP(status)) 
           return;
      }

      //#define BIASOPT_N 1
      else if (strcmp("biasoptPt",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->biasoptPt,
                              "biasoptPt",&biasoptInx,&paramNo,ROW_NUM,status);
        if (!StatusOkP(status)) 
           return;
      }

      //#define CABLE_N    3
      else if (strcmp("cableStart",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->minCABLE,&paramNo,status);
      else if (strcmp("cableStep",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->stepCABLE,&paramNo,status);
      else if (strcmp("cableNo",token)==0)
        sc2dalibsetup_servosetparamVal(delimiters,&setup->cableNo,&paramNo,status);

      //#define CABLEOFFSET_N 4
      else if (strcmp("cableOffset",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->cableOffset,
                              "cableOffset",&cableInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }

      //#define CABLEADJ_THD_N 4
      else if (strcmp("cableAdjthd",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->cableadjThd,
                              "cableAdjthd",&cableadjInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }

      //#define CABLEADJ_SCALE_N 4
      else if (strcmp("cableAdjscale",token)==0)
      {
        sc2dalibsetup_servosetparfloatVal(delimiters,setup->cableadjScale,
                              "cableAdjscale",&adjscaleInx,&paramNo,status);
        if (!StatusOkP(status)) 
           return;
      }
      //#define CABLEADJ_INIT_N 4
      else if (strcmp("cableadjInit",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->cableadjInit,
                              "cableadjInit",&cableinitInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }
      // #define GAIN_N      4
      else if (strcmp("Gain",token)==0)
      {
        sc2dalibsetup_servosetparfloatVal(delimiters,setup->gain,
                              "gain",&gainInx,&paramNo,status);
        if (!StatusOkP(status)) 
           return;
      }

      // #define GAINSCALE_N      4
      else if (strcmp("sq1gainScale",token)==0)
      {
        sc2dalibsetup_servosetparfloatVal(delimiters,setup->sq1gainScale,
                              "sq1gainScale",&sq1gainscaleInx,&paramNo,status);
        if (!StatusOkP(status)) 
           return;
      }

      //#define ZFACTOR_N   4
      else if (strcmp("zFactor",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->zfact,
                              "zfactor",&zfactInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }

      // #define SSABIASLCK_N   4
      else if (strcmp("ssabiaslckPt",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->ssabiaslckPt,
                              "ssabiaslckPt",&ssabiaslckInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }
      // #define BIASLCK_N   4
      else if (strcmp("biaslckPt",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->biaslckPt,
                              "biaslckPt",&biaslckInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }

      // #define FLUXPERIOD_N   4
      else if (strcmp("fluxPeriod",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->fluxPeriod,
                              "fluxPeriod",&fluxInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }
      // #define INITFB_N    4
      else if (strcmp("intFB",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->initFB,
                              "initFB",&initfbInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }
      // #define SQ2FDBKOPT_N    4
      else if (strcmp("sq2fdbkopt",token)==0)
      {
        sc2dalibsetup_servosetparamarrayVal(delimiters,setup->sq2fdbkOpt,
                              "sq2fdbkopt",&sq2fdbkoptInx,&paramNo,8,status);
        if (!StatusOkP(status)) 
           return;
      }

      //#define CALBESCALE_N   4
      else if (strcmp("cableScale",token)==0)
      {
        sc2dalibsetup_servosetparfloatVal(delimiters,setup->cableScale,
                              "cableScale",&scaleInx,&paramNo,status);
        if (!StatusOkP(status)) 
           return;
      }
      //#define CALBESLOPE_N   4
      else if (strcmp("cableSlope",token)==0)
      {
        sc2dalibsetup_servosetparfloatVal(delimiters,setup->cableSlope,
                              "cableSlope",&cableslopeInx,&paramNo,status);
        if (!StatusOkP(status)) 
           return;
      }

      else
      {
        *status=DITS__APP_ERROR;
        sprintf(errmsg,
          "sc2dalibsetup_servoreadSetup: %s in %s failed to match predefined name",
            token,myInfo->batchFile);
        ErsRep(0,status,errmsg);
        return;
      }
    }
  }
  jitDebug(16,"servoName=%s\n",setup->servoName);
  if (myInfo->actionFlag !=TRKSQ2FBACTION)
    sc2dalibsetup_servochkSetup(myInfo,setup,paramNo,status);
}


/**
 * \fn void sc2dalibsetup_servochkSetup(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *      int paramNo, StatusType *status)
 *
 * \brief function:
 *  check the setup file if anything is not right
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  paramNo int  no of parameter read
 * \param status  StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_servochkSetup
*/
void sc2dalibsetup_servochkSetup
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
int             paramNo,
StatusType      *status
)
{
  int         expectedNo=0;

  if (*status != STATUS__OK) return;
  
  if(strcmp(setup->servoName,"sq2open")==0)
  {
    setup->servo=SQ2OPEN;
  }
  else if(strcmp(setup->servoName,"sq2servo")==0)
  {
    setup->servo=SQ2SERVO;
  }
  else if(strcmp(setup->servoName,"ssalock")==0)
  {
     setup->servo=SSALOCK;
  }
  else if(strcmp(setup->servoName,"sq2biasing")==0)
  {
     setup->servo=SQ2BIASING;
  }
  else if(strcmp(setup->servoName,"sq1open")==0)
  {
     setup->servo=SQ1OPEN;
  }
  else if(strcmp(setup->servoName,"sq1servo")==0)
  {
     setup->servo=SQ1SERVO;
  }
  else if(strcmp(setup->servoName,"ssaramp")==0)
  {
    setup->servo=SSARAMP;
  }
  else if(strcmp(setup->servoName,"ssaramp1")==0)
  {
    setup->servo=SSARAMP1;
  }
  else if(strcmp(setup->servoName,"cablecal")==0)
  {
    setup->servo=CABLECAL;
  }
  else if(strcmp(setup->servoName,"blackbody")==0)
  {
    setup->servo=BLACKBODY;
  }
  else if(strcmp(setup->servoName,"sq1lock")==0)
  {
    setup->servo=SQ1LOCK;
  }
  else if(strcmp(setup->servoName,"heaterservo")==0)
  {
    setup->servo=HEATERSERVO;
  }
  else if(strcmp(setup->servoName,"tesbiasservo")==0)
  {
    setup->servo=TESBIASSERVO;
  }
  else if(strcmp(setup->servoName,"pixelmon")==0)
  {
    setup->servo=PIXELMON;
  }
  else if(strcmp(setup->servoName,"sq2lock")==0)
  {
    setup->servo=SQ2LOCK;
  }
  else if(strcmp(setup->servoName,"sq1biasservo")==0)
  {
    setup->servo=SQ1BIASSERVO;
  }
  else if(strcmp(setup->servoName,"testransit")==0)
  {
    setup->servo=TESTRANSIT;
  }
  else if(strcmp(setup->servoName,"cablecorrect")==0)
  {
    setup->servo=CABLECORET;
  }
  else if(strcmp(setup->servoName,"sq2open4p")==0)
  {
    setup->servo=SQ2OPEN4P;
  }
  else if(strcmp(setup->servoName,"sq2fbtrack")==0)
  {
    setup->servo=SQ2FBTRACK;
  }
  else if(strcmp(setup->servoName,"sq2fbgetpara")==0)
  {
    setup->servo=SQ2FBGETPARA;
  }

  else
  {
     *status=DITS__APP_ERROR;
     sprintf(errmsg, 
      "sc2dalibsetup_servochkSetup: unrecognised servoName (%s) ",
      setup->servoName);
     ErsRep(0,status,errmsg);
     return ;
  }
  *status=STATUS__OK;
  if (setup->servo==SSARAMP )    
  {
    if (paramNo !=SETUPPAR_SSARAMP)
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SSARAMP;
    }
  }
  if (setup->servo==SSARAMP1 )    
  {
    if (paramNo !=SETUPPAR_SSARAMP1)
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SSARAMP1;
    }
  }
  else if (setup->servo==SSALOCK )    
  {
    if (paramNo !=SETUPPAR_SSALOCK)
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SSALOCK;
    }
  }
  else if (setup->servo==SQ2BIASING) 
  {
    if ( paramNo !=SETUPPAR )    
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR;
    }
  }
  else if (setup->servo==BLACKBODY)
  { 
    if ( paramNo !=SETUPPAR )
       *status=DITS__APP_ERROR;
  } 
  else if (setup->servo==SQ2SERVO) 
  {
    if (paramNo !=SETUPPAR_SQ2) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SQ2;
    }
  }
  else if (setup->servo==SQ2OPEN ) 
  {
    if (paramNo !=SETUPPAR_SQ2OPEN) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SQ2OPEN;
    }
  }
  else if (setup->servo==SQ2OPEN4P ) 
  {
    if (paramNo !=SETUPPAR_SQ2OPEN4P) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SQ2OPEN4P;
    }
  }
  else if (setup->servo==SQ2LOCK ) 
  {
    if (paramNo !=SETUPPAR_SQ2LOCK) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SQ2LOCK;
    }
  }
  else if (setup->servo==SQ1SERVO) 
  {
    if (paramNo !=SETUPPAR_SQ1) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SQ1;
    }
  }
  else if (setup->servo==SQ1OPEN) 
  {
    if (paramNo !=SETUPPAR_SQ1OPEN) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SQ1OPEN;
    }
  }
  else if (setup->servo==SQ1LOCK) 
  {
    if (paramNo !=SETUPPAR_SQ1) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SQ1;
    }
  }
  else if (setup->servo==CABLECAL) 
  {
    if (paramNo !=SETUPPAR_CABLE) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_CABLE;
    }
  }
  else if (setup->servo==HEATERSERVO) 
  {
    if (paramNo !=SETUPPAR_HEAT) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_HEAT;
    }
  }
  else if (setup->servo==TESBIASSERVO) 
  {
    if (paramNo !=SETUPPAR_TES) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_TES;
    }
  }
  else if (setup->servo==PIXELMON) 
  {
    if (paramNo !=PIXELMON_PAR) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=PIXELMON_PAR;
    }
  }
  else if (setup->servo==SQ1BIASSERVO) 
  {
    if (paramNo !=SETUPPAR_SQ1BIAS) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SQ1BIAS;
    }
  }
  else if (setup->servo==TESTRANSIT) 
  {
    if (paramNo !=SETUPPAR_TRANSIT) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_TRANSIT;
    }
  }
  else if (setup->servo==CABLECORET) 
  {
    if (paramNo !=SETUPPAR_CABLECORET) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_CABLECORET;
    }
  }
  else if (setup->servo==SQ2FBTRACK || setup->servo==SQ2FBGETPARA ) 
  {
    if (paramNo !=SETUPPAR_SQ2FBTRACK) 
    {
      *status=DITS__APP_ERROR;
      expectedNo=SETUPPAR_SQ2FBTRACK;
    }
  }
  if (*status !=STATUS__OK)
  {   
    sprintf(errmsg, 
            "sc2dalibsetup_servochkSetup: the setup paras(%d) != expectedNo(%d)",
            paramNo,expectedNo);
    ErsRep(0,status, errmsg);
    return ;
  }
  if(setup->totalRC2use==0)
  {
     *status=DITS__APP_ERROR;
     sprintf(errmsg,"sc2dalibsetup_servochkSetup: whichRC has no readCard assigned");
     ErsRep(0,status,errmsg);
     return ;
  }
}
 

 
/**
 * \fn void sc2dalibsetup_servosetCard(char *delimiters, char *name,
 *     int *count,  StatusType *status)
 *
 * \brief function:
 *  read servo settings into struture
 *
 * \param  delimiters  char pointer for delimiters
 * \param  name      char pointer for setup-card or other name
 * \param  count     int pointer for param counting
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_servosetCard
*/
void sc2dalibsetup_servosetCard
(
char       *delimiters,
char       *name,
int        *count,
StatusType *status
)
{
  char *token;

  if (*status != STATUS__OK) return;
  
  // token => "the number after the word" 
  token = strtok (NULL, delimiters);      
  strcpy(name,token);
  jitDebug(16,"paramNo[%d] param= %s \n",*count,name);
  *count +=1;
}


/**
 * \fn void sc2dalibsetup_servosetparamVal(char *delimiters, int *param,
 *     int *count,  StatusType *status)
 *
 * \brief function:
 *  read servo settings into struture
 *
 * \param  delimiters  char pointer for delimiters
 * \param  param      int pointer for setup-parameter
 * \param  count     int pointer for param counting
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_servosetparamVal
*/
void sc2dalibsetup_servosetparamVal
(
char       *delimiters,
int        *param,
int        *count,
StatusType *status
)
{
  char *token;

  if (*status != STATUS__OK) return;
  
  // token => "the number after the word" 
  token = strtok (NULL, delimiters);      
  *param=atoi(token);
  jitDebug(16,"paramNo[%d] param= %d \n",*count,*param);
  *count +=1;
}


/**
 * \fn void sc2dalibsetup_servosetparamarrayVal(char *delimiters, int *paramarray,
 *     char *name, int *index, int *count, int howMany, StatusType *status)
 *
 * \brief function:
 *  read servo settings into array struture
 *
 * \param  delimiters  char pointer for delimiters
 * \param  paramarray int pointer for setup-parameter array
 * \param  name       char point for the param name
 * \param  index      int pointer for array index
 * \param  count      int pointer for param counting
 * \param  howMany    int how many param need to read in this call
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  char *token;
  int   read8=1,inx,chkInx=0;

  if (*status != STATUS__OK) return;
 
  jitDebug(16,"paramNo[%d]  entryInx[%d] %s[]= ",*count,*index,name);
  inx=*index;
  while(read8)
  {
    if( (token=strtok(NULL,delimiters))==NULL)      
    {
      *status=DITS__APP_ERROR;
      sprintf(errmsg,
        "sc2dalibsetup_servosetparamarrayVal:%s's paramamters < %d",name,howMany);
      ErsRep(0,status, errmsg);
      return ;
    }
    else
    {
      chkInx++;
      paramarray[inx]= atoi(token);
      jitDebug(16,"[%d]=%d ",inx,paramarray[inx]);
      inx++;
      if ( chkInx==howMany)
      {
        read8=0; *index=inx, jitDebug(16,"\n");
      }
    }
  }
  *count +=1;
}


/**
 * \fn void sc2dalibsetup_servosetparfloatVal(char *delimiters, double *paramarray,
 *     char *name, int *index, int *count,  StatusType *status)
 *
 * \brief function:
 *  read servo settings (float) into array struture
 *
 * \param  delimiters  char pointer for delimiters
 * \param  paramarray double pointer for setup-parameter array
 * \param  name       char point for the param name
 * \param  index      int pointer for array index
 * \param  count      int pointer for param counting
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  char *token;
  int   read8=1,inx;

  if (*status != STATUS__OK) return;
 
  jitDebug(16,"paramNo[%d]  entryInx[%d] %s[]= ",*count,*index,name);
  inx=*index;
  while(read8)
  {
    if( (token=strtok(NULL,delimiters))==NULL)      
    {
      *status=DITS__APP_ERROR;
      sprintf(errmsg,
        "sc2dalibsetup_servosetparamarrayflaotVal:%s's paramamters < 8",name);
      ErsRep(0,status,errmsg);
      return ;
    }
    else
    {
      paramarray[inx]= atof(token);
      jitDebug(16,"%f ",paramarray[inx]);
      inx++;
      if ( (inx & 0x07)==0)
      {
        read8=0; *index=inx, jitDebug(16,"\n");
      }
    }
  }
  *count +=1;
}


/**
 * \fn void sc2dalibsetup_servosetwhichRC(dasInfoStruct_t *myInfo,
 *  char *delimiters, ARRAYSET *setup,
 *  int *count,  StatusType *status)
 *
 * \brief function:
 *  read servo settings into whichRC struture
 *             also for parameter sharedMemory                       
 *
 * \param  myInfo    dasInfoStruct_t poiter
 * \param  delimiters  char pointer for delimiters
 * \param  setup      ARRAYSET structure pointer
 * \param  count      int pointer for param counting
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_servosetwhichRC
*/
void sc2dalibsetup_servosetwhichRC
(
dasInfoStruct_t *myInfo,
char       *delimiters,
ARRAYSET   *setup,
int        *count,
StatusType *status
)
{
  char *token;
  char  cardName[30];
  int   i,j=0,readRC=1;

  if (*status != STATUS__OK) return;

  jitDebug(16,"paramNo[%d] param=",*count);
  for (i=0;i<4;i++)
    setup->whichRC[i]=0;
  while (readRC )
  {
    if( (token=strtok(NULL,delimiters))==NULL)      
      readRC=0;
    else
    {
      setup->whichRC[j]= atoi (token);
      jitDebug(16," %d ",setup->whichRC[j]);
      if (setup->whichRC[j] <=0)
      {
        *status=DITS__APP_ERROR;
        sprintf(errmsg, 
         "sc2dalibsetup_servosetwhichRC: whichRC[%d]=(%d, must >=1)",
          j,setup->whichRC[j]);
        ErsRep(0,status,errmsg);
        return;
      }
      j++;
    }
  }
  jitDebug(16,"\n");
  setup->totalRC2use=j;

  // set up for parameter shared memory
  if (j==4)
    sc2dalibsetup_whichRC(myInfo,"rcs",status);
  else 
  {
    // there will only be either four or one RC
    for (i=0;i<4;i++)
    {
      if(setup->whichRC[i]!=0)
      {
        sprintf(cardName,"rc%d",setup->whichRC[i]);
        break;
      }
    }
    sc2dalibsetup_whichRC(myInfo,cardName,status);
  }
  *count +=1;
}



/**
 * \fn void sc2dalibsetup_servosetGain(char *delimiters, ARRAYSET *setup,
 *     int * index, int *count,  StatusType *status)
 *
 * \brief function:
 *  read servo settings into gain struture
 *
 * \param  delimiters  char pointer for delimiters
 * \param  setup      ARRAYSET structure pointer
 * \param  index      int pointer for array index
 * \param  count      int pointer for param counting
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_servosetGain
*/
void sc2dalibsetup_servosetGain
(
char       *delimiters,
ARRAYSET   *setup,
int         *index,
int        *count,
StatusType *status
)
{
  char *token;
  int   read8=1, inx;

  if (*status != STATUS__OK) return;

  jitDebug(16,"paramNo[%d] entryInx[%d]  Gain[]= \n",*count,*index);
  inx=*index;

  while(read8)
  {
    if( (token=strtok(NULL,delimiters))==NULL)      
    {
      *status=DITS__APP_ERROR;
      sprintf(errmsg,
        "sc2dalibsetup_servosetGain:gain paramamters < 8");
      ErsRep(0,status,errmsg);
      return ;
    }
    else
    {
      setup->gain[inx]= atof(token);
      jitDebug(16,"%f ",setup->gain[inx]);
      inx++;
      if ( (inx & 0x07)==0)
      {
        read8=0; *index=inx,jitDebug(16,"\n");
      }
    }
  }
  *count +=1;
}



/**
 * \fn void sc2dalibsetup_servosq1fdbkreadSetup(dasInfoStruct_t *myInfo, 
 *    ARRAYSET *setup, StatusType *status)
 *
 * \brief function:
 *  get sq1 fdbk initial Values from sq1servo for sq1open 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param status  StatusType pointer.  given and return
 *
 */

/*+ sc2dalibsetup_servosq1fdbkreadSetup
*/
void sc2dalibsetup_servosq1fdbkreadSetup
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char       delimiters[]= " =\n";
  char       readLine[2*FILE_LEN], tmpFile[FILE_LEN],*token;   
  int        rowNo=0,i;
  int        rowMap[32];

  if (!StatusOkP(status)) return;

  sprintf(tmpFile,"%s/%s", getenv ( "SC2SCRATCH" ),SQ1FBLOCKFILE  );
  if ((myInfo->fpOther2=fopen ( tmpFile, "r" ))==NULL)
  {
    *status=DITS__APP_ERROR;
    sprintf(errmsg, "sc2dalibsetup_servosq1fdbkreadSetup: failed to open %s, check permission", tmpFile); 
    ErsRep(0,status,errmsg);
    return;
  }

  for (i=0;i<ROW_NUM;i++)
     rowMap[i]=0;
  i=1;
  while( fgets(readLine,2*FILE_LEN,myInfo->fpOther2) !=NULL )
  {
    jitDebug(16,"\n %dth entry \n",i);
    if ( (strchr(readLine,'#')==NULL) && (strcmp(readLine,"\n")!=0) )
    {
      jitDebug(16,"%s \n",readLine);
      token = strtok (readLine, delimiters);
      if (strcmp("sel_row",token)==0)
      {
        token = strtok (NULL, delimiters);   
        rowNo= atoi (token);
        rowMap[rowNo]=1;
        jitDebug(16,"row =%d \n",rowNo);
      }  
      else if (strcmp("sq1initfdbk",token)==0)
      {      
        sc2dalibsetup_sq1biasfdbkPoints(delimiters,"sq1initfdbk",
                          setup->pixeLock.sq1initFB,  rowNo,NONMODULATION,status); 
        if (!StatusOkP(status)) 
        {
          fclose(myInfo->fpOther2);
          return;
        }
      }
      i++;
    }
  }
  for (i=0;i<ROW_NUM;i++)
  {
    if(rowMap[i]==0)
    {
      jitDebug(16,"row_%d not used, place NONMODULATION\n",i);
      sc2dalibsetup_sq1biasfdbk4notUsed("sq1initfdbk",setup->pixeLock.sq1initFB,i,
          NONMODULATION,status);
    }
  }
  fclose(myInfo->fpOther2);
}
  

/**
 * \fn void sc2dalibsetup_sq1biasfdbk4notUsed (char *name,
 *  int *paramarray, int row, int Val, StatusType *status)
 *
 * \brief function:
 *  put VAL__BADI into sq1 bias fdbk paramarray
 *
 * \param  name       char point for the param name
 * \param  paramarray int pointer for parameter array
 * \param  row        int  the row_No
 * \param  Val        int  set the value
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_sq1biasfdbk4notUsed
*/
void sc2dalibsetup_sq1biasfdbk4notUsed
(
char        *name,
int        *paramarray,
int         row,
int         Val,
StatusType *status
)
{

  int   ch; 
  int   arrayInx;

  if (*status != STATUS__OK) return;
 
  arrayInx = row*COL_NUM;

  for (ch=0;ch<COL_NUM;ch++)
      paramarray[arrayInx+ch]=Val;
}


/**
 * \fn void sc2dalibsetup_sq1biasfdbkPoints( char* delimiters,char *name,
 *  int *paramarray, int row, StatusType *status)
 *
 * \brief function:
 *  save sq1 bias fdbk points to paramarray
 *
 * \param  delimiter  char pointer for delimiters
 * \param  name       char point for the param name
 * \param  paramarray int pointer for parameter array
 * \param  row        int  the row_No
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{

  char *token;
  int   readCh=1,chkInx=0; 
  int   arrayInx=0;
  int   val;

  if (*status != STATUS__OK) return;
 
  arrayInx = row*COL_NUM;

  //printf("\n%s ",name);
  
  while(readCh)
  {
    if( (token=strtok(NULL,delimiters))==NULL)      
    {
      *status=DITS__APP_ERROR;
      printf(
        "sc2dalibsetup_sq1biasfdbkPoints: %s's dataNo < %d",name,COL_NUM);
      return ;
    }
    else
    {
      val= atoi(token);
      if (val==0)
        val=setVal;
      paramarray[arrayInx]=val;
      jitDebug(16,"%d ",paramarray[arrayInx]);
      if (arrayInx >( ROW_NUM*COL_NUM-3) )
         jitDebug(16,"Inx=%d \n",arrayInx);
      arrayInx++;
      chkInx++;
      if ( chkInx==COL_NUM)
      {
        readCh=0;
      }
    }
  }
}

/**
 * \fn void sc2dalibsetup_chksinglePixel(ARRAYSET *setup,  int *sq1comp, 
 *  int *sq2comp, int *sq1, int *sq2, int *flag, StatusType *status)
 *
 * \brief function:
 *  check if it is single pixel, if it is then pass sq1bias nad sq2fdbk 
 *  to sq1comp, sq2fcomp  
 *
 * \param  setup    ARRAYSET structure pointer
 * \param  sq1comp     int pointer, sq1bias optimal given,return; 
 * \param  sq2comp     int pointer, sq2fdbk optimal given,return
 * \param  sq1         int pointer, sq1bias given
 * \param  sq2         int pointer, sq2fdbk given
 * \param  flag        int pointer, given, return  
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{
  int row, col, biasInx=256, fdbkInx=256;
  int summy=0, sqnum;

  if (*status != STATUS__OK) return;

  for ( row=0; row<ROW_NUM; row++ )
  {
    if ( setup->rowMask[row] !=0)
    {
      summy ++;
      biasInx=row;
    }
  }
  for ( col=0; col<COL_NUM; col++ )
  {
    if ( setup->colMask[col] !=0)
    {
      summy ++;
      fdbkInx=col;
    }
  }

  if (summy >2)
  {
    *flag=0;
    return;
  }
  else
  {
    printf ("there is only one pixel(%d,%d), skip optimal\n",
             biasInx,fdbkInx);
    *flag=1;
    for ( row=0; row<ROW_NUM; row++ )
    {
      sq1comp[row]=VAL__BADI;
      for ( col=0; col<COL_NUM; col++ ) 
      {
         sqnum = row * COL_NUM + col;      
         if ( col ==fdbkInx && row ==biasInx)
           sq1comp[row]=sq1[sqnum];
      }
    }
    for ( col=0; col<COL_NUM; col++ ) 
    {
      sq2comp[col]=VAL__BADI;
      for ( row=0; row<ROW_NUM; row++ )
      {
        sqnum = row * COL_NUM + col;      
        if ( col ==fdbkInx && row ==biasInx)
          sq2comp[col]=sq2[sqnum];
      }
    }
  }
}


/**
 * \fn void sc2dalibsetup_sq1biassq2fdbkClipmean(int *sq1bias, int *sq1boptimal, int *sq2fdbk,
 *  int *sq2fdbk, int*sq2fdbkoptimal, int colNo, int rowNo, StatusType *status)
 *
 * \brief function:
 *  find sq1bias and sq2fdbk using clipmean 
 *
 * \param  sq1bias     int pointer  input sq1bias[colNo*rowNo]
 * \param  sq1boptimal int pointer  return sq1boptimal[rowNo]
 * \param  sq2fdbk     int pointer  input  sq2fdbk[colNo*rowNo]
 * \param  sq2foptimal int pointer  return sq2foptimal[colNo]
 * \param  colNo      int 
 * \param  rowNo      int 
 * \param status      StatusType pointer.  G&R
 *
 */
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
)
{
//#define DEB_OPT_ROW 0
//#define DEB_OPT_COL 1
  int     i, j, inx;
  double  data[ROW_NUM], meanVal;
  
  if ( !StatusOkP(status) ) return;

#ifdef DEB_OPT_ROW
  printf("SQ1 bias optiomal\n"); 
#endif

  for ( i=0; i<rowNo; i++ )
  {
     for ( j=0; j<colNo; j++ )
     {
        inx=i*colNo + j;
        if (sq1bias[inx] ==VAL__BADI)
          data[j]=VAL__BADD;
        else
          data[j]=(double)sq1bias[inx];
#ifdef DEB_OPT_ROW
        if ( i== DEB_OPT_ROW)
        {
          if (data[j] ==VAL__BADD)  printf("bad ");
          else                      printf("%3.0f ",data[j]);
        }
#endif
     }
     sc2math_clipmean(3.0, colNo, data, &meanVal,status);
     sq1boptimal[i]=(int)meanVal;
#ifdef DEB_OPT_ROW
     if ( i== DEB_OPT_ROW)
       printf("\n");
     printf("%d ",sq1boptimal[i]);
#endif
  }
#ifdef DEB_OPT_COL
  printf("\n"); 
  printf("SQ2 fdbk optiomal colNo(%d)rowNo(%d)\n",colNo,rowNo); 
#endif
  sleep(1);
 
  for ( i=0; i<colNo; i++ )
  {
     for ( j=0; j<rowNo; j++ )
     {
        inx = j*colNo + i;
        if ( sq2fdbk[inx]==VAL__BADI)
          data[j]=VAL__BADD; 
        else
          data[j]=(double)sq2fdbk[inx];
#ifdef DEB_OPT_COL
        if ( i== DEB_OPT_COL)
        {
          if (j==0) printf("\nsq2fdbk[0]=%d,following are for ch_%d:\n",sq2fdbk[0],i);
          if ( data[j]==VAL__BADD)   printf("bad ");
          else                       printf("%3.0f ",data[j]);
        }
#endif
     }
     sc2math_clipmean(3.0, rowNo, data, &meanVal,status);
     sq2foptimal[i]=(int)meanVal;
#ifdef DEB_OPT_COL
     if ( i== DEB_OPT_COL) 
       printf("\nsq2foptimal[0]=%d\n",sq2foptimal[0]);

     if( i==0) printf("sq2foptimal=:");
     printf("%d ",sq2foptimal[i]);
#endif
  }
  printf("\n");
}

/**
 * \fn void sc2dalibsetup_writesq1biasfdbkOut(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *  char *filename, int *sq1bias, int *sq2feedback, int *sq1refbias, 
 *  int *sq2reffeedback, int flag, StatusType *status)
 *
 * \brief function:
 *  save sq1 bias fdbk points to a text file 
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup   ARRAYSET structure pointer
 * \param  filename   char point for the file name
 * \param  sq1bias     int pointer 
 * \param  sq2feedback int pointer
 * \param  sq1refbias  int pointer
 * \param  sq2reffeedback int pointer
 * \param   flag      int 
 * \param status      StatusType pointer.  given and return
 *
 */
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
)
{

   FILE *fd, *fd1;    /* file descriptor for output */
   int col;           /* column index */
   int row;           /* row index */
   int sqnum;         /* sequential SQ1 count */
   char tmp[FILE_LEN];

   if ( !StatusOkP(status) ) return;
   
   if((fd = fopen ( filename, "w" )) == NULL )
     {
       *status = DITS__APP_ERROR;
       ErsRep (0, status, "sc2dalibsetup_writesq1biasfdbkOut: Error- failed to open file %s", filename); 
       return;
     }
   
   for ( col=0; col<COL_NUM; col++ ) 
   {
      for ( row=0; row<ROW_NUM; row++ )
      {
         sqnum = row * COL_NUM + col;
	 if( (sq2feedback[sqnum] <0 || sq2reffeedback[sqnum] <0) && flag==0)
	 {
	    sq1bias[sqnum]=VAL__BADI;
	    sq1refbias[sqnum]=VAL__BADI;
    	    sq2feedback[sqnum]=VAL__BADI;
    	    sq2reffeedback[sqnum]=VAL__BADI;
	 }	 
         fprintf ( fd, "%d  %d   %d  %d  %d  %d\n", 
           row+1, col+1, sq1bias[sqnum], sq2feedback[sqnum],
           sq1refbias[sqnum], sq2reffeedback[sqnum] );
      }  
   }
   fclose ( fd );   

/***** now we use squal from sc2sqopt_getSQ1scales  x.gao 08-nov-2007
  //sc2dalibsetup.c  bad pixel if sq2feebk=BAD_IVAL 
  //sc2dalibsetup_writepixelInx(myInfo,setup,filename,sq2feedback,
  //                   0,status);
***********/

   sprintf(tmp,"%s-sq1biassq2fb",filename);
   if((fd1 = fopen ( tmp, "w" )) == NULL )
     {
       *status = DITS__APP_ERROR;
       ErsRep (0, status, "sc2dalibsetup_writesq1biasfdbkOut: Error- failed to open file %s", tmp); 
       return;
     }
   fprintf ( fd1, "sq1bias= \n"); 
   for ( row=0; row<ROW_NUM; row++ )
   {
     fprintf ( fd1, "%2d_row ",row); 
     for ( col=0; col<COL_NUM; col++ ) 
     {
       sqnum = row * COL_NUM + col;
       if (sq1bias[sqnum]==VAL__BADI)
         fprintf ( fd1, "badval ");
       else 
         fprintf ( fd1, "%6d ",sq1bias[sqnum]); 
     }
     fprintf ( fd1, "\n"); 
   }

   fprintf ( fd1, "\nsq2fdbk= \n"); 
   for ( row=0; row<ROW_NUM; row++ )
   {
     fprintf ( fd1, "%2d_row ",row);
     for ( col=0; col<COL_NUM; col++ ) 
     {
       sqnum = row * COL_NUM + col;
       if (sq2feedback[sqnum]==VAL__BADI)
         fprintf ( fd1, "badval ");
       else 
         fprintf ( fd1, "%6d ",sq2feedback[sqnum]); 
     }
     fprintf ( fd1, "\n"); 
   }
   fclose ( fd1 );   
}

/**
 * \fn void sc2dalibsetup_sc2writequalOut( FILE *fd, double *sq1scale, int *squal,
 *  int *fqual, int *bqual, StatusType *status)
 *
 * \brief function:
 *  save middle result from sc2sqopt_getquals to a text file 
 *
 * \param  fd          FILE point for the file name
 * \param  sq1scale    double pointer 
 * \param  squal       int pointer
 * \param  fqual       int pointer, number of dad sq2fb/col
 * \param  bqual       int pointer
 * \param status       StatusType pointer.  given and return
 *
 */
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
)
{
  int    i,j,pixel;

  if (*status != STATUS__OK) return;

  fprintf ( fd, "# bqual\n");
  for ( i=0; i<ROW_NUM; i++ )
  {
    if(bqual[i]==VAL__BADI)
      fprintf ( fd, "badval " );
    else
      fprintf ( fd, "%6d ",bqual[i] );
  }
  fprintf ( fd, "\n# fqual\n");

  for ( i=0; i<COL_NUM; i++ )
  {
    if(fqual[i]==VAL__BADI)
      fprintf ( fd, "badval " );
    else
       fprintf ( fd, "%6d ",fqual[i] );
  }

  fprintf ( fd, "\n# squal\n");
  for ( i=0; i<ROW_NUM; i++ )
  {   
    for ( j=0; j<COL_NUM; j++ )
    {
      pixel=i*COL_NUM + j;    
      if(squal[pixel]==VAL__BADI)
        fprintf ( fd, "badval " );
      else
       fprintf ( fd, "%4d ",squal[pixel] );
    }
    fprintf( fd, "\n");
  }

  fprintf ( fd, "\n# sq1scale\n");

  for ( i=0; i<ROW_NUM; i++ )
  { 
    fprintf( fd, "%2d_row ",i);  
    for ( j=0; j<COL_NUM; j++ )
    {
      pixel=i*COL_NUM + j;    
      if(sq1scale[pixel]==VAL__BADD)
        fprintf ( fd, "badval " );
      else
        fprintf ( fd, "%6.3f ", sq1scale[pixel] );
    }
    fprintf( fd, "\n");
  }
}


/**
 * \fn void sc2dalibsetup_findalter4Opt(dasInfoStruct_t  *myInfo,    
 *   ARRAYSET *setup, int *optArray, int *orgArray, int flag,  
 *   int which,StatusType *status)
 *
 * \brief function:
 *  find mean value from the orgArray for the BAD sq1bias and sq2fb 
 *
 * \param  myInfo  dasInfoStruct_t poiter
 * \param  setup    ARRAYSET structure pointer
 * \param  optArray int pointer for optimal array
 * \param  orgArray int pointer for orginal array
 * \param  flag     int  0: col, 1: row
 * \param  which    int  which row or col 
 * \param status  StatusType pointer.  given and return
 *
 */
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
)
{
#define  TEST_ROW_4OPT 1000  // 0~40 do some printf, >40 not printf
  int  pixel, aver, i;
  int  data;

  if ( !StatusOkP(status) ) return;
  
  if(flag==1)  // row
  {
    if (which==TEST_ROW_4OPT)
      printf("check_row_%d \n",which);
    data=0, aver=0;
    for (i=0;i<COL_NUM;i++)
    {
      pixel=which*COL_NUM+i;
      if (which==TEST_ROW_4OPT)
        printf("pixel(%d)=%d ",pixel,orgArray[pixel]);

      if( orgArray[pixel] !=VAL__BADI )
      {
        data +=orgArray[pixel];
        aver ++;
      }
    }
    if (aver >0 )
    {
      optArray[which]=data/aver;
      // trap any value over the range
      if (optArray[which] < 0  || optArray[which] >=MAX_SQ1BIAS)
      {
        printf("row_%d: AVER value %d over the MCE range, set BAD\n",
               which,optArray[which]);
  
        optArray[which]=VAL__BADI;
      }
    }
    else
    {
      printf("row_%d: all BAD_VAL set BAD\n", which); 
      optArray[which]=VAL__BADI;
    }

    if (which==TEST_ROW_4OPT)
      printf("\n");
  }
  else   //col
  {
    data=0, aver=0;
    for (i=0;i<ROW_NUM;i++)
    {
      pixel=which +COL_NUM*i;
      if( orgArray[pixel] !=VAL__BADI )
      {
        data +=orgArray[pixel];
        aver ++;
      }
    }
    if (aver!=0)
      optArray[which]=data/aver;
  }
}


/**
 * \fn void sc2dalibsetup_writecompOut( FILE *fd, int *sq1bcomp,
 *  int *sq2fcomp, StatusType *status)
 *
 * \brief function:
 *  save middle result from sc2sqopt_fit to a text file 
 *
 * \param  fd          FILE point for the file name
 * \param  sq1bcomp    int pointer, 
 * \param  sq2fcomp    int pointer
 * \param status       StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_writesq1biasfdbkOut
*/
void sc2dalibsetup_writecompOut
(
FILE       *fd,
int        *sq1bcomp,
int        *sq2fcomp,
StatusType *status
)
{
  int   i;

  if (*status != STATUS__OK) return;

  fprintf ( fd, "# sq1bcomp\n");
  for ( i=0; i<ROW_NUM; i++ )
  {
    if (sq1bcomp[i]==VAL__BADI)
        fprintf ( fd, "badval ");
    else
     fprintf ( fd, "%6d ",sq1bcomp[i] );
  }

  fprintf ( fd, "\n# sq2fcomp \n" );
  for ( i=0; i<COL_NUM; i++ )
  {
    if (sq2fcomp[i]==VAL__BADI)
        fprintf ( fd, "badval ");
    else
     fprintf ( fd, "%6d ",sq2fcomp[i] );
  }
  fprintf ( fd, "\n" );
}

/**
 * \fn void sc2dalibsetup_readsq1Data(FILE *fp, int *bias,int *fdbk,
 *  int *refbias, int *reffdbk,StatusType *status)
 *
 * \brief funcation: read sq1 servo data
 *
 * \param  fp        FILE pointer
 * \param  bias      int pointer
 * \param  fbdk      int pointer
 * \param  refbias   int pointer
 * \param  reffdbk   int pointer
 * \param status  StatusType pointer.  given and return
 *
 */

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
)
{
//#define DEB_READSQ1SQ2_OPT_ROW 1
  char       delimiters[]= " =\n";
  char       readLine[FILE_LEN], *token;   
  char       lineCopy[FILE_LEN];
  int        rowNo=-1,i;
  int        rowMap[32];

  if (!StatusOkP(status)) return;

  for (i=0;i<ROW_NUM;i++)
     rowMap[i]=0;

  while(1)
  {
    if(  fgets(readLine,FILE_LEN,fp) !=NULL )
    {
      sprintf (lineCopy,"%s", readLine);
      token = strtok (readLine, delimiters);
      if (strcmp("sel_row",token)==0)
      {
        token = strtok (NULL, delimiters);   
        rowNo= atoi (token);
        rowMap[rowNo]=1;
        #ifdef DEB_READSQ1SQ2_OPT_ROW
          printf("\nrow =%d ",rowNo);
       #endif
      }  
      else if (strcmp("sq1bias",token)==0)  //sc2dalibsetup_sq1biasfdbkPoints from sc2dalibsetup.c
      {
        #ifdef DEB_READSQ1SQ2_OPT_ROW
         if ( rowNo == DEB_READSQ1SQ2_OPT_ROW )
          printf("sq1bias:row_%d\n%s\n",rowNo,lineCopy);
       #endif
        
        sc2dalibsetup_sq1biasfdbkPoints(delimiters,"sq1bias",bias,rowNo,VAL__BADI,status);                             
        if (!StatusOkP(status)) 
          return;
      }  
      else if (strcmp("sq2fdbk",token)==0)
      {
        sc2dalibsetup_sq1biasfdbkPoints(delimiters,"sq2fdbk",fdbk,rowNo,VAL__BADI,status);                             
        if (!StatusOkP(status)) 
          return;
      }
      else if (strcmp("sq1refbias",token)==0)
      {
        sc2dalibsetup_sq1biasfdbkPoints(delimiters,"sq1refbias",refbias,rowNo,VAL__BADI,status);                             
        if (!StatusOkP(status)) 
          return;
      }
      else if (strcmp("sq2reffdbk",token)==0)
      {
        sc2dalibsetup_sq1biasfdbkPoints(delimiters,"sq2reffdbk",reffdbk,rowNo,VAL__BADI,status);                             
        if (!StatusOkP(status)) 
          return;
      }
    }
    else
      break;
  }
  //printf("\n");
  for (i=0;i<ROW_NUM;i++)
  {
    if(rowMap[i]==0) //sc2dalibsetup_sq1biasfdbk4notUsed from sc2dalibsetup.c
    {
      printf("row_%d not used, place VAL__BADI \n",i);
      sc2dalibsetup_sq1biasfdbk4notUsed("sq1bias",bias,i,VAL__BADI,status);                             
      sc2dalibsetup_sq1biasfdbk4notUsed("sq2fdbk",fdbk,i,VAL__BADI,status);                             
      sc2dalibsetup_sq1biasfdbk4notUsed("sq1refbias",refbias,i,VAL__BADI,status);                             
      sc2dalibsetup_sq1biasfdbk4notUsed("sq2reffdbk",reffdbk,i,VAL__BADI,status);                             
    }
  }
}


/**
 * \fn void sc2dalibsetup_whichRC(dasInfoStruct_t *myInfo,    
 *  char *cardName,  StatusType *status)
 *
 * \brief function
 *  set whichRC[] in parameter shared memory 
 *
 * \param myInfo    dasInfoStruct_t pointer
 * \param cardName  char string
 * \param status     StatusType pointer
 * 
 */
/*+ sc2dalibsetup_whichRC
*/
void sc2dalibsetup_whichRC
(
dasInfoStruct_t  *myInfo,    
char             *cardName,  
StatusType       *status
)
{
  int i;
  PAR_SHARED *parentPtr;

  if (*status != STATUS__OK) return;

  // parameter shared Mempry
  parentPtr=(PAR_SHARED *) myInfo->parShm;
  for (i=0;i<4;i++)
    parentPtr->whichRC[i]=0;
 
  if ( strcmp (cardName,"rcs")==0)
  {
    for (i=0;i<4;i++)
      parentPtr->whichRC[i]=i+1;
  }
  else if ( strcmp (cardName,"rc1")==0)
    parentPtr->whichRC[0]=1;
  else if ( strcmp (cardName,"rc2")==0)
    parentPtr->whichRC[1]=2;
  else if ( strcmp (cardName,"rc3")==0)
    parentPtr->whichRC[2]=3;
  else if ( strcmp (cardName,"rc4")==0)
    parentPtr->whichRC[3]=4;
  else
    *status=DITS__APP_ERROR;
}




/**
 * \fn void sc2dalibsetup_readFile(dasInfoStruct_t *myInfo,
 *  char **filememPtr, int *fileLen, StatusType *status)
 *
 * \brief function:
 *   read  a file into memory, return the memory pointer
 *
 * \param myInfo     dasInfo structure pointer
 * \param filememPtr  double char pointer, for file memory address
 * \param fileLen     int pointer
 * \param *status       StatusType.  given and return
 *
 */
/*+ sc2dalibsetup_readFile*/
void sc2dalibsetup_readFile
(
dasInfoStruct_t *myInfo,
char            **filememPtr, 
int             *fileLen,
StatusType      *status
)
{
  int   fileLength,dataSize;  
  char  *memAddr;
  char   msg[FILE_LEN];

  if (*status != STATUS__OK) return;

  // find out how big the file is 
  if (fseek(myInfo->fpBatch, SEEK_SET, SEEK_END)) 
  {
    *status = DITS__APP_ERROR;
    sprintf(msg, "sc2dalibsetup_readFile: Strange problems with %s",
            myInfo->batchFile);
    ErsRep (0, status, msg);
    return;
  }
  fileLength = ftell(myInfo->fpBatch);
  if (fileLength < 0 ) 
  {    
    // -1L is error return 
    *status = DITS__APP_ERROR;
    sprintf(msg, 
            "sc2dalibsetup_readFile: Error getting file length (%d) of %s",
            fileLength,myInfo->batchFile);
    ErsRep (0, status, msg);
    return;
  }
   *fileLen=fileLength;

  // fseek to the beginning of the file 
  if (fseek(myInfo->fpBatch, 0, SEEK_SET)) 
  {
    *status = DITS__APP_ERROR;
    sprintf(msg, 
            "sc2dalibsetup_fpgareadFile: Error- Strange problems with %s",
            myInfo->batchFile);
    ErsRep (0, status, msg);
    return;
  }
  // get the Allocate storage for the file 
  if (  (memAddr = (char*)malloc(fileLength)) == NULL) 
  {
    *status = DITS__APP_ERROR;
    sprintf(msg, 
            "sc2dalibsetup_readFile: Insufficient memory for storing %s",
            myInfo->batchFile);
    ErsRep (0, status, msg);
    return;
  }
  // read the file into memory 
  dataSize = fread(memAddr, 1, fileLength, myInfo->fpBatch);
  if (dataSize== 0) 
  {
    *status = DITS__APP_ERROR;
    sprintf(msg, 
            "sc2dalibsetup_fpgareadFile: Unable to read %s",myInfo->batchFile);
    free(memAddr);
    ErsRep (0, status, msg); 
    return;
  }
  // pass the memory pointer back
  *filememPtr=memAddr;
}


/**
 * \fn void sc2dalibsetup_servodataoutInit(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *     int colrow,StatusType *status)
 *
 * \brief function
 *  write all servo settings into data file for other program/scripts
 *  to process data
 *
 * \param  myInfo    dasInfoStruct_t poiter
 * \param  setup      ARRAYSET structure pointer
 * \param  colrow     int if sq1biasservo, pass column if testransit pass row
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_servodataoutInit
*/
void sc2dalibsetup_servodataoutInit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,  
int             colrow,   
StatusType      *status
)
{
  time_t   tm;
  int      i;
   
  if (!StatusOkP(status))  return;

  tm = time(NULL);
  fprintf (myInfo->fpData,"<CONTROL>\n");
  fprintf (myInfo->fpData," started = %s",asctime(localtime(&tm)) );
  fprintf (myInfo->fpData," rc      = ");
  for (i=0;i<setup->totalRC2use;i++)
  { 
    fprintf (myInfo->fpData, "%d ",setup->whichRC[i]);
  }
  fprintf (myInfo->fpData, "\n"); 

  if (setup->servo==CABLECAL)
  {
    fprintf (myInfo->fpData, " servo   = %s: \n",
             setup->servoName); 
    fprintf (myInfo->fpData, " # bias is for cableOffset Val \n");
    fprintf (myInfo->fpData, " bias    = %d\n", setup->minCABLE); 
    fprintf (myInfo->fpData, " bstep   = %d\n",setup->stepCABLE); 
    fprintf (myInfo->fpData, " nbias   = %d\n",  setup->cableNo); 
  }
  else
  {
    fprintf (myInfo->fpData, " servo   = %s\n", setup->servoName); 
    fprintf (myInfo->fpData, " bias    = %d\n", setup->minBIAS); 
    fprintf (myInfo->fpData, " bstep   = %d\n",setup->stepBIAS); 
    fprintf (myInfo->fpData, " nbias   = %d\n",  setup->biasNo); 
  }
  fprintf (myInfo->fpData, " feed    = %d\n",setup->minFDBK); 
  fprintf (myInfo->fpData, " fstep   = %d\n",setup->stepFDBK); 
  fprintf (myInfo->fpData, " nfeed   = %d\n",setup->fdbkNo); 
  fprintf (myInfo->fpData, " nrow    = %d\n",setup->selRow);
  fprintf (myInfo->fpData, " doservo = %d\n",setup->doServo);

  if (setup->servo==SSARAMP  || setup->servo==SSARAMP1 )
  {
    fprintf (myInfo->fpData," cableoffset = ");
    for (i=0;i<32;i++)
       fprintf (myInfo->fpData, "%d ",setup->cableOffset[i]);
    fprintf (myInfo->fpData, "\n"); 

    if (setup->servo==SSARAMP )
    {
      fprintf (myInfo->fpData," cableScale    = ");
      for (i=0;i<32;i++)
        fprintf (myInfo->fpData, "%4.3E ",setup->cableScale[i]);
      fprintf (myInfo->fpData, "\n"); 

      fprintf (myInfo->fpData," cableadjThd    = ");
      for (i=0;i<32;i++)
        fprintf (myInfo->fpData, "%d ",setup->cableadjThd[i]);
      fprintf (myInfo->fpData, "\n"); 

      fprintf (myInfo->fpData," cableadjScale    = ");
      for (i=0;i<32;i++)
        fprintf (myInfo->fpData, "%4.3E ",setup->cableadjScale[i]);
      fprintf (myInfo->fpData, "\n"); 
    }
  }

  fprintf (myInfo->fpData," slopselect  = ");
  for (i=0;i<32;i++)
  { 
    fprintf (myInfo->fpData, "%d ",setup->slopSelect[i]);
  }
  fprintf (myInfo->fpData, "\n"); 

  fprintf (myInfo->fpData," gain        = ");
  for (i=0;i<32;i++)
  { 
    fprintf (myInfo->fpData, "%4.3E ",setup->gain[i]);
  }
  fprintf (myInfo->fpData, "\n"); 

  fprintf (myInfo->fpData," zfact       = ");
  for (i=0;i<32;i++)
  { 
    fprintf (myInfo->fpData, "%d ",setup->zfact[i]);
  }
  fprintf (myInfo->fpData, "\n"); 

  fprintf (myInfo->fpData," sabiaslck   = ");
  for (i=0;i<32;i++)
  { 
    fprintf (myInfo->fpData, "%d ",setup->biaslckPt[i]);
  }
  fprintf (myInfo->fpData, "\n"); 

  fprintf (myInfo->fpData," initfb       = ");
  for (i=0;i<32;i++)
  { 
    fprintf (myInfo->fpData, "%d ",setup->initFB[i]);
  }
  fprintf (myInfo->fpData, "\n"); 

  fprintf (myInfo->fpData," colMask       = ");
  for (i=0;i<32;i++)
  { 
    fprintf (myInfo->fpData, "%d ",setup->colMask[i]);
  }
  fprintf (myInfo->fpData, "\n"); 

  if (setup->servo==SQ1SERVO || setup->servo==SQ1OPEN ||
      setup->servo ==HEATERSERVO || setup->servo ==TESBIASSERVO ||
      setup->servo ==SQ1BIASSERVO)
  {
    fprintf (myInfo->fpData," rowMask       = ");
    for (i=0;i<ROW_NUM;i++)
    { 
      fprintf (myInfo->fpData, "%d ",setup->rowMask[i]);
    }
    fprintf (myInfo->fpData, "\n"); 
    fprintf (myInfo->fpData," fluxPeriod     = ");
    for (i=0;i<32;i++)
    { 
      fprintf (myInfo->fpData, "%d ",setup->fluxPeriod[i]);
    }

    fprintf (myInfo->fpData, "\n"); 
  }
  if( setup->servo ==SQ1BIASSERVO)
    fprintf (myInfo->fpData, " col     = %d\n",  colrow); 

  if( setup->servo ==TESTRANSIT)
  {
    fprintf (myInfo->fpData, " row     = %d\n",  colrow); 
    fprintf (myInfo->fpData, 
     " >: col0 is heatVal, col1 and col2 are org and filted data for COL0 and so on\n"); 
  }
  fprintf (myInfo->fpData,"</CONTROL>\n" );
  fprintf (myInfo->fpData,"\n<DATA>\n");
}



/********** transit point related *****************/
/*************************************************/

/**
 * \fn void sc2dalibsetup_heaterGen(ARRAYSET *setup,
 *  int funcFlag,StatusType *status)
 *
 * \brief function: 
 *  generate a heater modulation function 
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  funcFlag int   if 1: generate triangle, 0: not defined
 * \param  status  StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
/*+ sc2dalibsetup_heaterGen
*/
void sc2dalibsetup_heaterGen
(
ARRAYSET        *setup,
int             funcFlag,
StatusType      *status
)
{
  int     halfInx,i,heatModulate;
  int     inx=0, step, peakVal;

  if (*status != STATUS__OK) return;

  heatModulate=setup->fdbkNo*4+1; 

  if (funcFlag==1) // for tri-angle
  {
    if(setup->waveFlag==TRI_UP)  
    {
      step=setup->stepFDBK;
      peakVal=setup->minFDBK + setup->stepFDBK*setup->fdbkNo;
    }
    else 
    {
      step=-setup->stepFDBK;
      peakVal=setup->minFDBK - setup->stepFDBK*setup->fdbkNo;
    }

    setup->heater[0]=setup->minFDBK ;

    for(halfInx=0; halfInx<4; halfInx++)
    {
      for(i=1; i<=setup->fdbkNo; i++)
      {
        inx=halfInx*setup->fdbkNo +i;
        if( halfInx ==0 || halfInx==2)
          setup->heater[inx]=setup->heater[0] + step*i;
        else 
          setup->heater[inx]=peakVal - step*i;  
      }
    }
  }
}


/**
 * \fn void sc2dalibsetup_pixelTransit(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *     char *data,int whichRow, int heatLvl,int flag,FILE *fp, StatusType *status)
 *
 * \brief function: 
 *  check if the pixels on this row are in transition 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  data    char pointer for dataperheatlevel
 * \param  whichRow int 
 * \param  heatLvl  int 
 * \param  flag     int   1: use full frameSize, 0, use only pixelSize
 * \param  fp       FILE pointer for save the coefficient
 * \param  status  StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 * 
 */
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
)
{
  int    heat,heatModulate,data2; 
  int    *word,data1,col,row,pixel;
  double  coeff;
  double *transVal, *heatPower, *current;
  int    *lockFlag; 
  
  if (!StatusOkP(status)) return;

  heatModulate=setup->fdbkNo*4+1; 

  transVal=setup->transitVal+setup->totalPixel*heatLvl;
  lockFlag=setup->lockFlag+setup->totalPixel*heatLvl;

  heatPower=(double *)calloc(heatModulate,sizeof(double));
  current=(double *)calloc(heatModulate,sizeof(double));

  sc2dalibsetup_heater2Power(myInfo,setup,heatPower,status);
  fprintf(fp,"correlation coefficient\n");
 
  for(row=0; row<ROW_NUM; row++)
  {
    for (col=0; col<COL_NUM; col++)
    {
      pixel=row*COL_NUM +col;

      if (setup->colMask[col] ==0)
      {
        lockFlag[pixel]=0;
        transVal[pixel]=0;
      }
      else
      {
        // point to the pixel data on the first frame
        // data is in bytes, HEADERSET is in bytes too
        if (flag==1)
          word = (int *)(data + pixel*4 + HEADERSET);
        else
          word = (int *)(data + pixel*4);
  
        for (heat=0; heat<heatModulate; heat++)
        { 
          // collect all pixel per heat step in heat modulation
          // chData =(double *) COL_NUM*heatModulate   
          // allData =(int *) COL_NUM*heatModulate 
          // FRAMESIZE_INT is in INT
          if ( flag==1)
            data1=*(word + heat*FRAMESIZE_INT)*setup->colMask[col];
          else
            data1=*(word + heat*setup->totalPixel)*setup->colMask[col];

          data2=data1;
          // arrange allData[] heat0: col0,col1,....
          //                   heat1: col0,col1,....
          setup->allData[heat*heatModulate +col]=data1;
         
          if (setup->dataMode==4)
          {
            data1 = (data2 >> FDBKSHIFT) & FDBKMASK;
            sc2dalibsetup_signExt(myInfo,&data1,SIGNBIT_18BIT,DATAMASK_17BIT,SIGN_EXT_18,status);
          }
          setup->chData[heat]=(double)data1;
        }
        // we have all data for the heat modulation for this col,row
        sc2dalibsetup_correlateCoeff(heatModulate,setup->chData,heatPower,&coeff,status);
        transVal[pixel]=coeff;
      }
      fprintf(fp,"%6.2f ",coeff);
    }
    fprintf(fp,"\n");
    sc2dalibsetup_saveTransit(myInfo,setup,row,whichRow,status);
  }
  free(heatPower);
  free(current);
} 


/**
 * \fn void sc2dalibsetup_heater2Power(dasInfoStruct_t *myInfo, 
 *  ARRAYSET  *setup, double *pwer, StatusType *status)
 *
 * \brief function: 
 *  convert heater to Power   in pWater
 *
 * \param  myInfo dasInfoStruct_t pointer
 * \param setup   ARRAYSET      pointer
 * \param  pwer   double pointer 
 * \param  status  StatusType pointer.  given and return
 *  
 * P=r*I^2
 *  ...... 
 */
/*+ sc2dalibsetup_heater2Power
*/
void sc2dalibsetup_heater2Power
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
double          *pwer,
StatusType      *status
)
{
  double     fullDC=65535, current=24.8;
  double     pixelResist=2,c;
  int        i,heatModulation;

  if (*status != STATUS__OK) return;
  
  // I =c*heat;
  c= current/fullDC;
  heatModulation=setup->fdbkNo*4+1;

  for (i=0; i<heatModulation; i++)
    pwer[i]=pixelResist* pow( (c*setup->heater[i]),2);
}



/**
 * \fn void sc2dalibsetup_sq1fb2Current(dasInfoStruct_t *myInfo, 
 *  ARRAYSET  *setup, double *current, StatusType *status)
 *
 * \brief function: 
 *  convert s1fb to current 
 *
 * \param  myInfo  dasInfoStruct_t pointer
 * \param  setup   ARRAYSET  pointer
 * \param  current double pointer 
 * \param  status  StatusType pointer.  given and return
 *  
 */
/*+ sc2dalibsetup_heater2Power
*/
void sc2dalibsetup_sq1fb2Current
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
double          *current,
StatusType      *status
)
{
}



/**
 * \fn void sc2dalibsetup_signExt(dasInfoStruct_t *myInfo,
 *  int *input, int signBit, int  dataMask,
 *  int  signextMask,StatusType *status)
 *
 * \brief function: 
 *  extend the data sign if it is negative
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  input  int pointer the data need to sign extend
 * \param  signBit int  the bit mask used for sign 
 * \param  dataMask int  used to get data without sign 
 * \param  signextMask int used to extend the sign
 * \param  status  StatusType pointer.  given and return
 *
 */
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
)
{
  int  tmp;

  tmp=*input;

  //printf("signBit (%#0X) dataMask (%#0X) signextMask(%#0X) \n",
  //        signBit,dataMask,signextMask);
  if ( (tmp & signBit) > 0)
  {
    tmp = tmp & dataMask;
    tmp = tmp | signextMask;
  }
  *input=tmp;
}



/**
 * \fn void sc2dalibsetup_correlateCoeff(int np, double *x, double *y, double *coeff, 
 *     StatusType *status)
 *
 * \brief function: 
 *  calculated correlated coeffieient between heat modulate and sq1fb
 *
 * \param  np      int number of points
 * \param  x       double pointer X data
 * \param  y       double pointer y data
 * \param  coeff   double pointer correlated coefficient
 * \param  status  StatusType pointer.  given and return
 *
 */
/*+ sc2dalibsetup_correlateCoeff
*/
void sc2dalibsetup_correlateCoeff
(
int np,               /* number of points (given) */
double *x,           /* X data (given) */
double *y,           /* Y data (given) */
double *coeff,       /* correlation coefficient (returned) */
StatusType *status   /* global status (given and returned) */
) 
{
  double xm;           /* X mean */
  double ym;           /* Y mean */
  double xynum;         /* XY factor */
  double xden;         /* X squared factor */
  double yden;         /* Y squared factor */
  double xyden;         /* XY squared factor */
  double ximxm;        /* X differences */
  double yimxm;        /* y differences */
  int i; 

  if ( !StatusOkP(status) ) return;

  xm = 0.0;
  ym = 0.0;
   
  for ( i=0; i<np; i++ )
   {
      xm += x[i];
      ym += y[i] ;
   }
   xm = xm / np;
   ym = ym / np;

   xynum = 0.0;
   xden = yden =0.0;

   for ( i=0; i<np; i++ )
   {
      ximxm = ( x[i] - xm ) ;
      yimxm = ( y[i] - ym ) ;
      xynum += ximxm * yimxm;
      xden  += ximxm * ximxm;
      yden  += yimxm * yimxm;
   }
   xyden = sqrt(xden*yden);
   if(xyden > 0 )
     *coeff = xynum / xyden;
   else
     *coeff = 0;
}



/**
 * \fn void sc2dalibsetup_transitResult(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *     FILE *fp, StatusType *status)
 *
 * \brief funcation: 
 *  write the transit result for other software to use
 *
 * \param myInfo    dasInfo structure pointer
 * \param setup    ARRAYSET structure pointer
 * \param fp       FILE pointer for result file    
 * \param status  StatusType pointer.  given and return
 *
 */

/*+ sc2dalibsetup_transitResult
*/
void sc2dalibsetup_transitResult
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
FILE            *fp,
StatusType      *status
)
{
  int    i,row,col,pixel; 
  double *transVal=NULL;
  int     *lockFlag;
  int     heatStart; 

  if (!StatusOkP(status)) return;

  heatStart=setup->minFDBK;

  //setup->minFDBK is heatLvl,setup->stepBIAS is haeatLvlStep
  for ( row=0; row<ROW_NUM; row++ )
  {
    fprintf (fp,"%02d_row PowerVal \n",row); 
    setup->minFDBK=heatStart;
    // step through heaterLevel
    for (i=0;i<setup->biasNo;i++)
    {
      transVal=setup->transitVal+setup->totalPixel*i;

      //  write out pixel testransit per heatLvl
      fprintf (fp,"heatVal=%d  ",setup->minFDBK); 
      for ( col=0; col<COL_NUM; col++ ) 
      {
        pixel = row * COL_NUM + col;
        fprintf (fp,"%8.4E ",transVal[pixel]); 
      }
      fprintf ( fp, "\n");
      setup->minFDBK -=setup->stepBIAS;
    }
  }

  if (setup->slopSelect[13] > 0)
    fprintf (fp,"lockFlag \n"); 
  else
    fprintf (fp,"meanerrVal\n");

  for ( row=0; row<ROW_NUM; row++ )
  {
    fprintf (fp,"%02d_row \n",row); 
    setup->minFDBK=heatStart;
    // step through heaterLevel
    for (i=0;i<setup->biasNo;i++)
    {
      lockFlag=setup->lockFlag+setup->totalPixel*i;
      //  write out pixel testransit per heatLvl
      fprintf (fp,"heatVal=%d  ",setup->minFDBK); 
      for ( col=0; col<COL_NUM; col++ ) 
      {
        pixel = row * COL_NUM + col;
        fprintf (fp,"%6d ",lockFlag[pixel]); 
      }
      fprintf ( fp, "\n");
      setup->minFDBK -=setup->stepBIAS;
    }
  }

  // the product of error*powerLevel
  for ( row=0; row<ROW_NUM; row++ )
  {
    fprintf (fp,"%02d_row error*powerLvl \n",row); 
    setup->minFDBK=heatStart;
    // step through heaterLevel
    for (i=0;i<setup->biasNo;i++)
    {
      lockFlag=setup->lockFlag+setup->totalPixel*i;
      //  write out pixel testransit per heatLvl
      fprintf (fp,"heatVal=%d  ",setup->minFDBK); 
      for ( col=0; col<COL_NUM; col++ ) 
      {
        pixel = row * COL_NUM + col;
        fprintf (fp,"%8.4E ",lockFlag[pixel]*transVal[pixel]); 
      }
      fprintf ( fp, "\n");
      setup->minFDBK -=setup->stepBIAS;
    }
  }
}


/**
 * \fn void sc2dalibsetup_saveTransit(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *    int row ,int whichRow, StatusType *status)
 *
 * \brief function: 
 *  save the orginal data for whichRow 
 *
 * \param  myInfo   dasInfoStruct_t poiter
 * \param  setup    ARRAYSET structure pointer
 * \param   row      int 
 * \param  whichRow int  whichRow for special
 * \param  status   StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 */
/*+ sc2dalibsetup_saveTransit
*/
void sc2dalibsetup_saveTransit
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
int             row,
int             whichRow,
StatusType      *status
)
{ 
  int   col,heat,heatModulate,colHeat; 


  if (!StatusOkP(status)) return;

  if( row==whichRow)
  {
    heatModulate=setup->fdbkNo*4+1; 
    fprintf(myInfo->fpData,"\n");

    for (heat=0; heat<heatModulate; heat++)
    { 
      colHeat=heatModulate*heat;
      fprintf(myInfo->fpData,"\n%7d ",setup->heater[heat]);
      for (col=0; col<COL_NUM; col++ )  
      {
        fprintf(myInfo->fpData,"%7d ",setup->allData[colHeat+col]);
      }
    }
    fprintf(myInfo->fpData,"\n");
  }
}


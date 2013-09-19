/**
 * \file findsq1optpoints.c
 *
 * \brief stand alone application for finding sq1 optimised points
 *
 * see Dennis' sqsim.c
 *  sqsim.c - write a simulated file of SQ1 optimation measurements
 *   Method :
 *   Four numbers are generated for each first-stage SQUID in a 41x32
 *   subarray and written to an XML file to provide test data for the
 *   sqsolve program.
 *
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
 *  $Log: findsq1optpoints.c,v $
 *  Revision 1.16  2011/06/02 21:31:22  bgorges
 *  Checking all of fopen and fopen64 to make sure a file is actually opened.
 *
 *  Revision 1.15  2010/08/03 22:16:10  cwalther
 *  Changes to make sc2_setup work with SC2SCRATCH
 *
 *  Revision 1.14  2009/10/15 03:17:58  xgao
 *  take out printf for result files
 *
 *  Revision 1.13  2009/10/15 03:12:49  xgao
 *  move functions to sc2dalibsetup.c
 *
 *  Revision 1.12  2009/10/15 03:01:55  xgao
 *  filtered version (use slopSelect[15]>=0, default 1) + clipmean for sq1bias(ROW) sq2fdbk(COL) as optimalVal
 *
 *  Revision 1.11  2009/09/02 00:33:06  cwalther
 *  Added code for overriding the optical code reader, fixed some items Tim wanted fixed in the Data file
 *
 *  Revision 1.10  2008/07/02 21:28:25  xgao
 *  change function names as sc2dalib_ sc2dalibsetup_ for sc2dalib.c sc2dalibsetup.c
 *
 *  Revision 1.9  2008/06/05 00:43:42  cwalther
 *  Changes to get things working at JAC
 *
 *  Revision 1.8  2008/03/07 12:47:03  xgao
 *  tidy up after doxygen 0p44
 *
 *  Revision 1.7  2007/12/07 16:01:19  xgao
 *
 *  some modifications after SG array test, shutter only third open, new
 *  initialise.xml ( flatfile, INITFLAG) and configure.xml (DATA_MODE).
 *  negative pixHeat for loadv in SETU_SEQ. Add findheaterspike for processing
 *  heaterservo data for G/Tc.
 *
 *  Revision 1.6  2007/11/02 12:14:08  xgao
 *
 *  update to 0p33-021107, worked with JOC for INIT,CONFIG, SETUPSEQ, SEQ
 *
 *  Revision 1.5  2007/09/05 16:33:47  xgao
 *
 *  add heaterslope, trkHeater tested Ok, sq1Refbias and sq2Reffb always use the last one
 *  always do manual adjustment after optimal, fits4all is now in parShm, bigPhy=512M now
 *  0p21, ndf is turned on, tested with startmce-rtsc flatfield-seg
 *
 *  Revision 1.4  2007/08/06 15:57:40  xgao
 *  re-organise servo setupfile, use include=. 0p16-nondf
 *
 *  Revision 1.3  2007/08/02 21:09:40  xgao
 *  fix three bugs after autoset run,0p15-nondf
 *
 *  Revision 1.2  2007/05/28 17:16:15  xgao
 *  change glbInfo to myInfo move sc2dareadmceval from libservo to lib
 *
 *  Revision 1.1.1.1  2007/05/16 08:26:56  dkelly
 *  first insertion
 *
 *
 */
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "jit.h"
#include "jitXML.h"
#include "DitsParam.h"

#include "Dits_Err.h"
#include "Ers.h"
#include "DitsMsgOut.h"     /* for MsgOut           */

#include "sae_par.h"
#include "prm_par.h"
#include "dat_par.h"
#include "star/hds.h"

#include "dream_par.h"
#include "dream.h"
#include "sc2math.h"
#include "sc2sqopt_par.h"
#include "sc2sqopt.h"

#include <sdsudriver_par.h>  // COL_NUM, ROW_NUM
#include <sdsudriver_struct.h>
#include <interface.h>
#include <mcexml_par.h>
#include <mcexml_struct.h>
#include <mcexml.h>

#include "sc2da_par.h"  // need some define
#include "sc2da_struct.h"  // need structure
#include "sc2dalibsetup.h"  // need readSetup

#include "sc2dalibsetup.c"  // need sc2dalibsetup_*



int main
(
int argc,
char **argv
)
{
//#define  TEST_READIN_DATA

   int i, flag, singleFlag=0;
   int whichRow=0;
   int sq1bias[COL_NUM*ROW_NUM];
   int sq2feedback[COL_NUM*ROW_NUM]; 
   int sq1biascp[COL_NUM*ROW_NUM];
   int sq2feedbackcp[COL_NUM*ROW_NUM]; 

   int sq1refbias[COL_NUM*ROW_NUM];
   int sq2reffeedback[COL_NUM*ROW_NUM];

   /* scale factor relating SQ1 bias and SQ2 feedback */
   double sq1scale[ROW_NUM*COL_NUM];
   
   int bqual[ROW_NUM];         /* quality for each SQ1 bias row */
   int fqual[COL_NUM];         /* quality for each SQ2 feedback column */
   int squal[ROW_NUM*COL_NUM]; /* quality for each SQ1 */

   int sq1bcomp[ROW_NUM];      /* compromise SQ1 bias for each row */
   int sq2fcomp[COL_NUM];     /* compromise SQ2 feedback for each column */

   FILE *df;                /* pointer to data file*/
   FILE *fd;                /* pointer to output file*/
   FILE *fpopen,*fpfb;                
   char datafile[FILE_LEN]; /* datafile written by DAS */
   char outfile[FILE_LEN]; /* datafile written by DAS */
   char outfileXml[FILE_LEN]; /* xml file */
   char tmpfile[FILE_LEN];
   StatusType mystatus;
   StatusType *status = &mystatus;
   
   dasInfoStruct_t  dasInfo;    
   ARRAYSET         arrayset;

   *status = STATUS__OK;
   
  // check command-line arguments 
  if ( argc < 4 )
  {
    printf ("\nusage: findsq1optpoints sq1servofile outfile flag whichRow\n\n" );
    printf (" sq1servofile = $CURRENTDATADIR/sq1opt-points  written by DAS  \n" );
    printf (" outfile = name of file for optimised result\n");
    printf (" flag = 1 takeout outlier of sq2fb set outlier=BAD;  flag=0 don't take out\n");
    printf ("        3 just use sc2math_clipmean for sq1bias (each row), sq2fdbk (each column)\n");
    printf (" whichRow: option if TEST_READIN_DATA is defined, print thisROW reading\n");
    printf (" (sq1servofile need full path or in current dir \n");
    printf (" (outfile is in $CURRENTDATADIR dir), it reads $CONFIG_ALL/SQ1SERVOFILE too\n\n"); 
    return -1;
  }

  strcpy ( datafile, argv[1] );
  sprintf( outfile,"%s/%s",getenv("CURRENTDATADIR"), argv[2] );
  flag=atoi(argv[3]);
  sprintf(outfileXml,"%s-xml",outfile);
 
  if ( argc > 4)
    whichRow=atoi(argv[4]);
  
  if((df = fopen(datafile,"r")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status,"Error- findsq1points.c main()1 failed to open file %s",datafile); 
      return 1;
    }

  printf("call ...sc2dalibsetup_readsq1Data\n");
  // read the datafile into memory buffer servoData
  sc2dalibsetup_readsq1Data(df,sq1bias, sq2feedback, sq1refbias, sq2reffeedback,status);
  fclose (df);
  if (*status !=  STATUS__OK) 
  {
    return 1;
  }
  #ifdef TEST_READIN_DATA
  {
    int  *data[2],j, pixel;
  
    data[0]=sq1bias;
    data[1]=sq2feedback;

    for (i=0; i<2;i++)
    {
      if (i==0)
        printf("sq1bias =\n");
      else if (i==1)
        printf("sq2fdbk =\n");

      for ( j=0; j<COL_NUM; j++ )
      {
        pixel=whichRow*COL_NUM + j;
        printf("%d ",data[i][pixel]);
      }
      printf("\n");
    }
  }
  #endif

  if (flag ==3)
  {
    printf("call ...sc2dalibsetup_sq1biassq2fdbkClipmean\n");
    sc2dalibsetup_sq1biassq2fdbkClipmean(sq1bias, sq1bcomp, sq2feedback, sq2fcomp,COL_NUM,ROW_NUM, status);
  }
  else
  {
    //printf("call ...sc2dalibsetup_writesq1biasfdbkOut\n");
    sprintf(tmpfile, "%s-pixel",outfile);
    //printf("findsq1optpoints: %s\n",tmpfile);

    sc2dalibsetup_writesq1biasfdbkOut (&dasInfo,&arrayset, tmpfile, sq1bias, sq2feedback, sq1refbias,
       sq2reffeedback, 1,status);

    //printf("call ...sc2dalibsetup_servoreadSetup\n");
    sprintf(dasInfo.batchFile, "%s/%s",getenv("CONFIG_ALL"),SQ1SERVOFILE);
    sc2dalibsetup_servoreadsetupWrap(&dasInfo, status);

    if((dasInfo.fpBatch = fopen(dasInfo.batchFile, "r")) == NULL )
      {
	*status = DITS__APP_ERROR;
	ErsRep (0, status, "Error- findsq1optpoints.c main()2 failed to open file %s", dasInfo.batchFile); 
	return 1;
      }
    //printf("read: %s\n",dasInfo.batchFile);
    // read in SQ1SERVOFILE so that we have fluxperiod[col]
    // have to set  myInfo->parShm, as sc2dasetwhichRC use it
    dasInfo.parShm=malloc(1000);
    sc2dalibsetup_servoreadSetup(&dasInfo,&arrayset, status);
    fclose(dasInfo.fpBatch);
    free(dasInfo.parShm);

    if ( *status != STATUS__OK )
    {
      printf("sc2dalibsetup_servoreadSetup failed ");
      return 1;
    }
    // check if it is single pixel, if it is then pass sq1bias and sq2fdbk to sq1comp, sq2fcomp 
    sc2dalibsetup_chksinglePixel(&arrayset,sq1bcomp, sq2fcomp,sq1bias,sq2feedback, &singleFlag, status);
    if (singleFlag ==0)
    {
      if (flag==1)
      {
        printf("call ...sc2dalibsetup_sq2fbOutliers\n");
        sprintf(tmpfile, "%s-outlier",outfile);
        printf("4nextStep: %s\n",tmpfile);
        sc2dalibsetup_sq2fbOutlier( &dasInfo,&arrayset,tmpfile, sq2feedback,
                       sq2reffeedback, status);
        if (*status != STATUS__OK ) 
          return 1;
      }
      // make copies of sq1bias and sq2feedback, it seems that  sc2sqopt_putSQ1 changes them
      memcpy(sq1biascp, sq1bias,sizeof(int)*COL_NUM*ROW_NUM);
      memcpy(sq2feedbackcp, sq2feedback,sizeof(int)*COL_NUM*ROW_NUM);
  
      printf("call ...sc2sqopt_put\n");
      sc2sqopt_putSQ1 ( outfileXml, sq1bias, sq2feedback, sq1refbias, sq2reffeedback, status );
      if(*status != STATUS__OK)
      {
        printf("sc2sqopt_putSQ1 return bad status\n");
        return 1;
      }

      printf("call ...sc2sqopt_getSq1scales\n");
      //change from sc2sqopt_getquals to sc2sqopt_getSQ1scales  07-Nov-2007 x.gao
      sc2sqopt_getSQ1scales ( sq1bias, sq2feedback, sq1refbias,sq2reffeedback,sq1scale, squal, fqual,bqual, status );
      if( *status != STATUS__OK)
      {
        printf("sc2sqopt_getquals return bad status\n");
        return 1;
      }

      // write out
      sprintf(tmpfile, "%s-pixel",outfile);
      // sc2dalibsetup.c  bad pixel if squal=1 
      sc2dalibsetup_writepixelInx(&dasInfo,&arrayset,tmpfile,squal,0, status);

      if((fd = fopen(outfile, "w")) == NULL )
	{
	  *status = DITS__APP_ERROR;
	  ErsRep (0, status, "Error- findsq1optpoints.c main()3 failed to open file %s", outfile); 
	  return 1;
	}
      // status is not changed 
      sc2dalibsetup_sc2writequalOut(fd,sq1scale, squal, fqual, bqual, status );

      /* Perform the fit */
      printf("call ...sc2sqopt_fit\n");
      sc2sqopt_fit ( sq1bias, sq2feedback, sq1scale, squal, fqual,bqual, sq1bcomp, sq2fcomp, status );

      /* Write the solution, status no change */
      sc2dalibsetup_writecompOut(fd,sq1bcomp,sq2fcomp, status );

      // for whatever reason, we will do the manual set
      if(*status == STATUS__OK)
        printf("we are going to get AVER for those rejected by optimal function\n");
      else
      {
        printf("optimal function failed, set sq1bcomp/sq2fcomp=BAD, going to get AVER \n");
        for ( i=0; i<ROW_NUM; i++ )
          sq1bcomp[i]=VAL__BADI;
        for ( i=0; i<COL_NUM; i++ )
          sq2fcomp[i]=VAL__BADI;
      }
      // now manually chang the BAD ones
      fprintf ( fd, "sq1bias manually changed\n" );
      for ( i=0; i<ROW_NUM; i++ )
      {
        if (sq1bcomp[i]==VAL__BADI)
        {
          if (arrayset.rowMask[i] !=0) // not masked out
            sc2dalibsetup_findalter4Opt(&dasInfo,&arrayset,sq1bcomp,sq1biascp,1,i, status);
          fprintf ( fd, "row_%d %6d rowMak=%d \n",i,sq1bcomp[i],arrayset.rowMask[i] );
        }
      }
      fprintf ( fd, "sq2fdbk manually changed\n" );
      for ( i=0; i<COL_NUM; i++ )
      {
        if (sq2fcomp[i]==VAL__BADI)
        {
          if( arrayset.colMask[i] !=0)
            sc2dalibsetup_findalter4Opt(&dasInfo,&arrayset,sq2fcomp,sq2feedbackcp,0,i,status);
          fprintf ( fd, "col_%d %6d  colMask=%d \n",i, sq2fcomp[i],arrayset.colMask[i] );
        }
      }
      fprintf ( fd, "\n" );
      fclose ( fd );
    }
  }
  // save sq1bcomp to sq1biasoptimal.txt
  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ),SQ1BIASOPTFILE );

  if((fpopen = fopen(tmpfile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "Error- findsq1optpoints.c main()4 failed to open file %s", tmpfile); 
      return 1;
    }
  fprintf(fpopen,"# BIASLCK_N  sq1bias Val for optimal, set 0 for VAL__BADI \n");
  fprintf ( fpopen, "wb ac on_bias <\n" );  
  sc2dalibsetup_savesq1biasoptPts(fpopen,"biasoptPt= ",sq1bcomp,1, status);
  fclose(fpopen);

  // save sq2fcomp to sq2fboptimal.txt
  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ), SQ2FBOPTFILE );

  if((fpfb = fopen(tmpfile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "Error- findsq1optpoints.c main()5 failed to open file %s", tmpfile); 
      return 1;
    }
  fprintf(fpfb, "#this is sq2fb optimal val \n");
  fprintf(fpfb,"wb bc2 flux_fb < \n");
  sc2dalibsetup_savesq2fboptPts(fpfb,"sq2fbopt= ",sq2fcomp,1, status);
  fprintf(fpfb,"\n");
  fclose(fpfb);

  // save sq1bcomp to sq1open.txt
  sprintf (tmpfile, "%s/%s", getenv ( "CONFIG_HARD" ), OPTIMAL_SQ1B_SQ2FB );

  if((fpopen = fopen(tmpfile, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "Error- findsq1optpoints.c main()6 failed to open file %s", tmpfile); 
      return 1;
    }
  fprintf(fpopen,"\n# BIASLCK_N  sq1bias Val for optimal, set 0 for VAL__BADI \n");
  sc2dalibsetup_savesq1biasoptPts(fpopen,"biasoptPt= ",sq1bcomp,0, status);
 
  fprintf(fpopen, "\n#this is sq2fdbk optimal val, BAD =0 \n");
  sc2dalibsetup_savesq2fboptPts(fpopen,"sq2fdbkopt= ",sq2fcomp,0, status);
  fclose(fpopen);

  return *status;
}



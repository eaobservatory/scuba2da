#ifndef sc2da_struct_h
#define sc2da_struct_h
/**
 * \file sc2da_struct.h
 *
 * \brief structure definition for SCUBA2 Data Acquisition Software (DAS)  
 *
 * \author Xiaofeng Gao, UKATC (xg@roe.ac.uk)
 *
 * \version 
 *
  *
 * Copyright (c) 2005.  UK Astronomy Technology Centre (UK ATC), 
 * An establishment of the Particle Physics and Astronomy Research Council
 * PPARC).
 *
 * Web: www.roe.ac.uk
 *
 *  $Log: sc2da_struct.h,v $
 *  Revision 1.52  2012/10/19 00:56:46  cwalther
 *  changes to make sure all arrays have the exact same OBSID
 *
 *  Revision 1.51  2011/12/20 23:23:14  cwalther
 *  Implemented a style of heater tracking that ends after a fixed number of iterations
 *
 *  Revision 1.50  2011/09/27 21:14:22  cwalther
 *  Moved the heater tracking header write so that I would print when you track to the value you remember in the dark, some changes to make remembering in dark work better
 *
 *  Revision 1.49  2011/05/05 19:51:32  cwalther
 *  Fixed the zillion of messages problem when heater tracking fails
 *
 *  Revision 1.48  2011/03/01 23:47:27  cwalther
 *  Added test for maximum gaini, set to zero if greater than maximum
 *
 *  Revision 1.1.1.1  2007/05/16 08:27:00  dkelly
 *  first insertion
 *
 */

// #ifndef FILE_LEN
//  #define FILE_LEN 100
// #endif

// #ifndef CMD_LEN
//  #define CMD_LEN 290
// #endif
// standard C

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <sys/ipc.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <features.h>
#include <sys/shm.h>
#include <sys/mman.h>

#include <sdsudriver_err.h> 
#include <sdsudriver_errstring.h> 
#include <sdsudriver_par.h> 
#include <sdsudriver_struct.h>
#include <interface.h>
#include <mcexml_par.h>
#include <mcexml_struct.h>
#include <mcexml.h>

#include "Dits_Err.h"       /* for error codes      */
#include "DitsTypes.h"      /* For basic Dits types */
#include "DitsFix.h"        /* Fixed part details   */
#include "DitsInteraction.h"/* For DitsTrigger      */
#include "Ers.h"            /* For Ers routines     */
#include "mess.h"           /* For Mess routines    */
#include "arg.h"            /* For Arg routines     */
#include "DitsMsgOut.h"     /* for MsgOut           */
#include "DitsSys.h"        /* For Dits system routines         */
#include "DitsUtil.h"       /* For DitsErrorText                */
#include "Sdp.h"            /*  for the parameter system */
#include "Git.h"
#include "sds.h"            /* for SDS stuff */
#include "jit.h"            /* For jit stuff */
#include "DitsParam.h"
#include "DitsSignal.h"
#include "DitsOrphan.h"

//parameter name
typedef struct  
{
  char  logfileDir[FILE_LEN];
  char  logfileName[FILE_LEN];
  char  mcexmlFile[FILE_LEN];
} MCE_PARAM ;


// structure for DSP *.lod indentifier label 
typedef struct  
{
  char      start[10];	// _START 
  char      data[10];   // _DATA  
  char      end[10];    // _END   
} SDSU_DSP_LOD_IND ;

// structure for HEALTH Status 
typedef struct  
{
  int      part;
  uint32   mask;
  char     name[20];
  uint32   memAdr;
} PCI_HEALTH ;

typedef struct 
{
  char name[10];
  int  arg1;
} MCE_CARD_DEF ;

typedef struct 
{
  char id[10];
  char chipid[40];
  char band[10];
  char health[10];
  char task[40];
  char flatfile[FILE_LEN];
  char weightfile[FILE_LEN];
  int  mceport;
  int  x;
  int  y;
  int  override;
  int  darkHeaterI;
  int  maxGainI;
} ARRAY_DEF ;


typedef struct 
{
  char band[10];
  char label[10];
  char units[10];
  double centre;
  double width;
} WAVE_BAND ;

typedef struct 
{
  char    mceCmd[CMD_LEN*2];      // for passing mce CMD 
  char    cmdBuf[MCEBLK_SIZE];
  int     cmdBufSize;
  int     replySize;
  MCE_REPLY  reply;
} dasCmdInfo_t;

typedef struct 
{
  int    stVal;                       // heater startValue   
  int    step;                        // heater step value
  int    nStep;                       // number of step 
  int    row;                         //  
  int    col;                         // [row][col] is the pixel we are look at
  int    option;                      // option 0=>No strip chart file, 1=>strip chart file
  int    flag;                        // option 0=>Simple PI controller, 1=>Gao's method
  int    ref;                         // the heater setting point where we take refFDBK
  int    refFDBK;                     // the reference sq1feedback Value in heater tracking
  int    refHeat;                     // the heater setting at the ref
  int    maxOffset;                   // maxoffset from the initHeat
  double curTime;                     // Current time in TAI
  double pGain;                       // Gain for the proportional control section
  double iGain;                       // Gain for the integral control section
  double slope;                       // the heater slope obtained from the heater data
  char   base[FILE_LEN];              // the base used for heater datafile
  char   name[FILE_LEN];              // the base used for heater datafile
} HEATTEST;

// structure for shared memory between Linux user space processors
// NDF header structure sc2head from sc2store_struct.h
// QL header and data
// dark row is not included in IMAGE data
typedef struct
{
 char    fits[MAXQLFITS*80+2];  //
 double  data[(ROW_NUM-1)*COL_NUM]; //< the max required data   [x][y]
} QL_IMAGE ;

typedef struct
{
 char    fits[MAXQLFITS*80+2];  //
 double  data[(ROW_NUM-1+4)*(COL_NUM+4)]; 
} RECONSTR_IMAGE ;


typedef struct
{
 int      framenum;   // 
 double   timeStamp;  //
 char     fits[MAXQLFITS*80+2];  //
 double   data[(ROW_NUM-1)*COL_NUM]; //< the max required data   [x][y]
 int      subStart;
 int      subEnd;
 int      nrec;      // fits records
} QL_STRUCT ;

typedef struct
{
 int      framenum;   // 
 double   timeStamp;  //
 char    fits[MAXQLFITS*80+2];  //
 double  data[(ROW_NUM-1+4)*(COL_NUM+4)]; 
 int      subStart;
 int      subEnd;
 int      nrec;      // fits records
} RECONSTR_STRUCT ;

typedef struct
{
 int      framenum;       // 
 double   timeStamp;  //
 char     filename[FILE_LEN];
 int      subStart;
 int      subEnd;
} QL_SCAN ;


typedef struct
{
  int minOut[COL_NUM];
  int maxOut[COL_NUM];
  int peakVal[COL_NUM];
  int peakInx[COL_NUM];
  int bin[COL_NUM];
  int changeInx[ROW_NUM*COL_NUM];
  int  freqPtr[BINDIV+1];  // An array store the frequency
} OPT_STRUCT ;


// structure for parameter-shared-Memory
typedef struct 
{
int   fSize;          //= glbInfo->bufSize/4;
int   colxRow;        //= fSize-FRAMEHEADER_NUM-CHKSUM_NUM;
int   realCol;        // = colxRow/ROW_NUM;
int   realRow;        //=(ROW_NUM-1);
int   cardNo;         //= realCol/ONECARDCHS;
int   fstatSize;      //=FRAMEHEADER_NUM*sizeof(int);  
int   onecardrowSize; //=ONECARDCHS*sizeof(int);
int   onecarddataSize;    //=realRow*onecardrowSize;
int   chdataSize;     //=realCol*realRow*sizeof(int);
int   darkSize;       //=cardNo*ONECARDCHS*sizeof(int);
int   whichRC[4];     // []=0 no RC, otherwise [0]=1; [1]=2; [2]=3; [3]=4 
int   frameNum;       // the freamenum 
int   qlframeNum;      // the freamenum for QL, reset =0 in CONFIGURE
int   engFlag;        //   RTSC_MODE or   RTSC_MODE; 
int   obsMode;          
int   obsNo;
#define OBSIDSS_SIZE 30
char  obsidss[OBSIDSS_SIZE];  // The observation ID for special subsystem
int   qlNum;          // the real number for how many COADD, RECONSTR images
int   startSeq;       // start number of seq
int   procNo;         // number of frame used to do stare/dream/scan
int   dreamRound;
int   smuPattern;         // SMU steps for doing dream
int   dreamDim[2];     // reconstructed dream image Dim
int   subscanNo;
char  arrayName[40];
char  baseFile[40];
int   load;                     // sky, dark, hot
int   stairWidth;              // Number of samples per stair in any sawtooth or ramp
int   sawtoothRampFlag;       // True if we are doing sawtooths or ramps 1=heater saw, 2=heater ramp 3=bias saw 4=bias ramp
size_t actfits;
int   heat_track_num;         // Number of bolometers used in heater tracking
int   heat_track_values[HEAT_TRACK_MAX][3]; //3xn heater tracking (col,row,heaterSetting)
size_t   frameHeader;             // number of framestatus per frame 
char  fitshd[MAXFITS*80+1];     //FITS head strings
char  flatfile[FILE_LEN];       //Name of flatfield file
char  configFile[FILE_LEN];       //Name of configure.xml file
char  dreamweightFile[FILE_LEN];  //Name of dream weight file (sdf)
double dut1;                      //DUT1 = UT1 - UTC (sec)
double tel_latdeg;                //Telescope latitude (deg)
double tel_longdeg;               //Telescope longitude (deg)
double instap_x;                  //Instrument aperture x offset (arcsec)
double instap_y;                  //Instrument aperture y offset (arcsec)
int    sq2fbparaFlag;             // flag 1: do SQ2FB update during SEQ, 0:don't
int    chdataCount;     // count in DH for how many ch data being copied into memory
int    fstateCount;     // count in DH for how many mceframeState being copied into memory
int    darkrowCount;    // count in DH for how many darkrow being copied into memory
int    coaddscandreamCount; // count in DH for how many coaddframe/reconstr/scan beiing
                            // copied into sharedmemory

}  PAR_SHARED;


// structure for parameter-shared-Memory
typedef struct 
{
  char   file[FILE_LEN];      //Name of flatfield file
  size_t colsize;             /* number of pixels in column  */
  size_t rowsize;            /* number of pixels in row  */
  size_t nflat;              /* number of flat coeffs per bol  */
  char   flatname[FILE_LEN]; /* name of flatfield algorithm  */
  double *flatcal;           /* pointer to flatfield calibration  */
  double *flatpar;           /* pointer to flatfield parameters  */
  double refres;             /* Reference resistor used to create flatfield */
}  FLAT_FIELD;


#define MAX_INBEAM 150

typedef struct 
{ 
  int      msgwrtPt;          // the message write index for intermag[]
  int      msgreadPt;          // the message read index for intermag[]
  int      sharedmId;         // shared memory Segment id 
  int      sharedmparId;      // parameter shared memory Segment id 
  int      doneReadxml;          // flag for having read MCE.xml OK  
  int      lastFrame;            // flag for lastframe   
  int      frameStopped;        
  char     devName[FILE_LEN];    // for device name 
  char     cardId[CARDID_LEN];   // for cardID name 
  char     logFile[FILE_LEN];    // for log file name  
  char     logfileName[FILE_LEN];
  char     drvVersion[FILE_LEN];
  char     givenarrayName[10];
  char     okfileName[FILE_LEN];  // for completion notification
  char     baseFile[FILE_LEN];    //  date used for raw file 
  char     Date[20];             // for holding the date:ie. 14012005 
  char     stopCmd[CMD_LEN];
  char     goCmd[CMD_LEN];
  int      actionFlag;
  int      logfileFlag;
  int      cmdFlag;               // mainly for onthefly cmd 
  int      batchFlag;
  int      scriptFlag;
  int      initFlag;             // flag set in initialise.xml INITFLAGS
  int      dataFormat;
  int      engFlag;               // flag to indicate if it is ENG or RTSC mode
  int      gotsq2fbparaFlag;      // 1 if SQ2FB_GETPARA has called, 0 not
  int      rowLength;     
  int      numRows;       
  int      obsNo;                 // obs No
  int      subscanNo;             // subscan No
  char     utcshort[OBSIDSS_SIZE]; // The UT date in YYYYMMDDTHHMMSS format
  long     actIndex;              //Drama index of active action 
  long     debuglvl; 
  long     bufSize;
  char     initFile[FILE_LEN];    // initialise.xml file  
  char     batchFile[FILE_LEN];   // by user, for batchfile         
  char     cmdrepFile[FILE_LEN];  // by user, for cmdreply to/from MCE  
  char     dataFile[FILE_LEN];    // by user, for data received from MCE
  char     strchartFile[FILE_LEN];// by user, for data saved for STRIPCHARTE
  long     filesvFlag;            // 1: data from MCE/SDSU_TIMING board 
  long     simFlag;              // flag for simulation 0: no 
  long     datasplitFlag ;        // split data from error signal
  long     divbynFlag;            // if set, data read from MCE =data/sampleNo 
  char     *sharedShm;           // shared memory segment pointer
  char     *parShm;              // parameter shared memory segment pointer
  FILE     *fpMcecmd;
  FILE     *fpData;
  FILE     *fpLog;             
  FILE     *fpBatch;
  FILE     *fpStrchart;          
  FILE     *fpOtheruse;
  FILE     *fpSq1;
  FILE     *fpOther2;          
  int      trkNo; 
  int      procNum;
  int      sharedmemSize;
  int      numFrame;
  int      lastframeNo;
  int      qlNum;           // how many number for QL 
  int      ithseqWait;      // wait for ith sequence 
  int      taskchkFlag;     // 1 if all monitored tasks STATE are collected
  int      headtaskchkFlag; // 1 if it is from monitored task
  int      firstQLFlag;     // 1 if it is firstQL trig
  int      headersDone;     // 1 if all monitored headers completed
  uint32   glbCount;
  uint32   bufpostOffset;        
  SdsIdType     parsysId;
  SdsIdType     statId;
  SdsIdType     qlId;
  SdsIdType     seqId;
  SdsIdType     timeId;
  SdsIdType     filenameId;
  SdsIdType     imageId;
  SdsIdType     qldataId;
  SdsIdType     fitsId;
  caddr_t       bufAdr;            // Data buffer address  */
  char          taskList[10][20];     // for store the task to find path
  ARRAY_DEF     subArray[MAX_SUBARRAY];
  WAVE_BAND     waveBand[2];
  char          chipId[40];
  char          filter[10];
  double        wavelen;           // in meters i.e 8.5*10^-7
  int           darkHeaterI;       // in DtoA counts default heater current in the dark
  int           heatTrkStyle;      // Style of heater tracking to perform
  int           heatTrackFailed;   // The last heater tracking has failed flag
  long          numTimesThru;      // Number of times to heater track when fixed length tracking (0=until kicked)
  int           maxGainI;          // Maximum magnitude of gaini for this array
  char         loadv[FILE_LEN];
  char         load[FILE_LEN];
  int          nominalPixelHeat;
  int          pixelHeat;
  int          detbias;
  int          datamode;
  int          bbHeat;
  int          en_fb_jump;
  float        shutterFraction;
  int          seqcount;
  int          drcontrol;
  int          stairStart;
  int          stairHeight;
  int          stairHeightCnts;
  int          stairWidth;
  int          stairNum;
  int          stairQCycleCount;
  int          stairHCycleCount;
  int          stairNumCount;
  int          stairPresentValue;
  int          sawtoothRampFlag;
  char         inbeam[MAX_INBEAM];
  int          astMapState;
  HEATTEST     heatSlp;
  int         *qlPtr;
  int         *lkupflagentPtr;
  QL_STRUCT        *sharemqlPtr;
  RECONSTR_STRUCT  *sharemrecnstrPtr;
  QL_SCAN          *sharemscanPtr;
  PAR_SHARED       *parshmPtr;
  FLAT_FIELD       flatField;
  struct JCMTState *jcmtheadEntry;
  int         headerwriteNo; /* write offset for the lookupTable entry */
  int         headerreadNo;  /* reda offset for the lookupTable entry */
  int         askSDF;        /* flag =1 if DRAMA waits for SCAN file completed */
  char        *servobindataPtr; // pointer for save idl servo data in BINARY
} dasInfoStruct_t;


// need for passing more arg to DitsAltInAdd which only accepts one
typedef struct 
{
  SDSU_CONTEXT    *con;      //  
  dasInfoStruct_t *myInfo;
  DRAMA_INNERMSG  *intermsg;
} USED4DITS;


// need to create memory space for coadded frame for QL.
// not used channels are set to zero 
// fits header for coadded data
typedef struct 
{
  char   fits[MAXQLFITS*80+2];
} FITS4COADD;


typedef struct 
{
  int        framenum;      //  
  double     time;
  FITS4COADD coaddfits;
  double     maxdata[COL_NUM][ROW_NUM-1];
} USED4COADD;


typedef struct 
{
  char   servo[SERV_STR_LEN];   // servo name
  char    card[SERV_STR_LEN];  // which readout card 
  char    bias[SERV_STR_LEN];  // bias start value 
  char   bstep[SERV_STR_LEN];  // bias step value
  char   nbias[SERV_STR_LEN];  // number of bias points
  char  fdback[SERV_STR_LEN];  // feed back start value
  char   fstep[SERV_STR_LEN];  // feed back step value
  char nfdback[SERV_STR_LEN];  // number of feed back points 
  char     row[SERV_STR_LEN];  // which row, only used in sq1 
  char doservo[SERV_STR_LEN];  // choice for servo algorithm
  char safbCard[SERV_STR_LEN];
  char sq2fbCard[SERV_STR_LEN];
  char sq2biasCard[SERV_STR_LEN];
  char sq1biasCard[SERV_STR_LEN];
  char initFB[SERV_STR_LEN*4];  // sq2fb or sq1fb for SQ2OPEN SQ1OPEN 
  char cableOff[SERV_STR_LEN*4]; // cable offset fo SSA
  char biasLck[SERV_STR_LEN*4];  // biaslck (sa/sq2/sq1) for SSALOCK,SQ2OPEN SQ1OPEN 
  char slopSelect[SERV_STR_LEN*4];  // which slop (1:positive, -1 negative) 
  char   colMask[SERV_STR_LEN*4];  // col 0: not used  
  char   rowMask[SERV_STR_LEN*4];  // row:0 not used 
  char   cableSlop[SERV_STR_LEN*4];  // 
  char   cableScale[SERV_STR_LEN*4];  // 
} SERVO_DATAHEAD;

typedef struct 
{
  char  started[64];   // the date when the servo file runs
  char  servo[64];  // servo name
  char  card[64];   // which readout card 
  int   bias;       // bias start value 
  int   bstep;      // bias step value
  int   nbias;      // number of bias points
  int   feed;      // feed back start value
  int   fstep;      // feed back step value
  int   nfeed;     // number of feed back points 
  int   nrow;        // which row, only used in sq1 
  int   doservo;    // choice for servo algorithm
  int   cableoffset[COL_NUM];   // cable offset fo SSA
  float cableScale[COL_NUM];    // cable Scale fo SSA  
  float cableSlope[COL_NUM];    // cable slope fo SSA  
  int   cableadjThd[COL_NUM];   // thread for CABLECAL,SSA,SSALOCK, CABLECORRET 
  float cableadjScale[COL_NUM]; // adj Scale thread for CABLECAL,SSA,SSALOCK, CABLECORRET 
  int   slopselect[COL_NUM];    // multi use 
  int   initFB[COL_NUM];        // sq2fb or sq1fb for SQ2OPEN SQ1OPEN 
  int   sabiaslck[COL_NUM];    // biaslck (sa/sq2/sq1) for SSALOCK,SQ2OPEN SQ1OPEN 
  float gain[COL_NUM];
  int   zfact[COL_NUM];
  int   fluxPeriod[COL_NUM];
  int   colMask[COL_NUM];  // col 0: not used  
  int   rowMask[ROW_NUM];  // row:0 not used 
  int   row;            // only for SQ1BIASSERVO or TESTRANSIT ( std-alone programe)
} SERVOBIN_DATAHEAD;

//COL_NUM, ROW_NUM are defined in sdsu_driver_struct.h
typedef struct 
{
 // column data in case the checksum is included
  uint32     dataCol[COL_NUM+HEADERS+1];  
} DATACOL;

/// peak infomation per bias set 
typedef struct 
{
  double peakPts[COL_NUM][4];    // peak points in the feebbk ramp 
  int    peakInx[COL_NUM][4];    // feedbk inedx for corresponding 
                                 //peak points  
  double p2pValue[COL_NUM][4];   // peak to peak value for each channel  
  // place here for recording at optimal point in sq2servo,
  // the ssafb and zfact, gain values  
  int    initFB[COL_NUM];         // fb for ssafb (sq2servo)   
  int    zfact[COL_NUM];          // ssa out 
  double gain[COL_NUM];           // gain 
  int    lowest[COL_NUM];         // the lowest point in the ch data
  int    highest[COL_NUM];        // the highest point in the ch data
} BIAS_PEAK;


/// maxpeak infomation from all bias set 
typedef struct 
{
  double maxp2p[COL_NUM][3];  // maximum peak in each channel for MAX-MIN, MIN-MAX,
                              // MAX-MIN
  int   maxinx0[COL_NUM][3];  // biasStep where maxp2p occurs
  int   maxinx1[COL_NUM][4];  // feedbk inedx for corresponding 
                              //peak points 
} MAX_P2P;


typedef struct 
{
  int    saOut[COL_NUM][3];     // ssa output value   
  double peakVal[COL_NUM][3];   // max modulation value  
  int    biasStep[COL_NUM][3];  // the bias step when peakVal occurs 
  int    fdbk[COL_NUM][4];      // feedbk inedx for corresponding 
                               //peak points 
  int    initFB[COL_NUM];     // fb for ssafb (sq2servo)   
  int    zfact[COL_NUM];      // ssa out 
  double gain[COL_NUM];       // gain 

} CHKPOINT;


typedef struct 
{
  int    ssalockPt[ROW_NUM*COL_NUM];   // ssa output value   
  int    sq1FB[ROW_NUM*COL_NUM];       // sq1 feedback at lock point 
  int    sq1initFB[ROW_NUM*COL_NUM];   // sq1 inital fdbk at lock point
                                       // from sq1servo 
  int    iVal[ROW_NUM*COL_NUM];        // I-Val for MCE (gaini)   
  int    fluxPeriod[ROW_NUM*COL_NUM];  // fluxperiod (flx_quanta)
  int    midVal[ROW_NUM*COL_NUM];      //  the first mid
  int    midInx[ROW_NUM*COL_NUM];      // index where mid occurs   
  int    peakVal[ROW_NUM*COL_NUM][2];  //  min or max peak
  int    peakInx[ROW_NUM*COL_NUM][2];  // index where min or max peak occurs   
  double gain[ROW_NUM*COL_NUM];        // gain from open loop  
} PIXEL;

 
typedef struct 
{
  char      safbCard[30];
  char      sq2fbCard[30];
  char      sq2biasCard[30];
  char      sq1fbCard[30];
  char      sq1biasCard[30];
  char      servoName[30];
  int       totalRC2use;   /* how many readCards to use*/ 
  int       servo;
  int       selRow;        /* wich row's data to use in the frame */ 
  int       whichRC[4];   /* wich readCards */ 
  int       doServo;      /* 0 not to apply servo algorithm,others do */ 
  int       bias;         /* instant value */ 
  int       fdbk;         /* instant value */ 
  int       biasCount;    /* count for bias ramp function */
  int       fdbkCount;    /* count1 for feedback ramp function */
  int       fdCount;      /* count2 for feedback ramp function */
  int       minBIAS;      /* minimum bias Value for ramp function  */
  int       minFDBK;      /* minimum feedback Value for ramp function  */
  int       minCABLE;      /* minimum cableoffset Value for cableoffset cal */
  int       stepBIAS;     /* steps for bias ramp function  */
  int       stepFDBK;     /* steps for feedback ramp function  */
  int       stepCABLE;     /* steps for cable cal  */
  int       biasNo;       /* number of steps for bias ramp function  */
  int       fdbkNo;       /* number of steps for feedback ramp function  */
  int       cableNo;       /* number of steps for cable cal */
  int       fbFlag;       /* flag for setting ssaFb in SQ2 or sq2Fb in SQ1 SERVO  */
  int       biasFlag;     /* flag for bias setting */
  int       adjFlag;      /* flag for adjust SSA FB setting in SQ2*/
  int       modFlag;      /* flag for modulation or derivative */
  int       minFlag;      /* flag for minimum modulation value found */
  int       maxFlag;      /* flag for maximum modulation value found */
  int       divFlag;      /* flag for divide MCE data by sampleNo */
  int       waveFlag;      /* flag for different tri-angle heater wave */
  int       totalPixel;    /* total pixels for array */
  int    cableOffset[COL_NUM];     /* circuits related paramater */
  double cableSlope[COL_NUM];      /* cable slope in ssaOut(cableOffset) */ 
  double cableScale[COL_NUM];      /* cable scale factor */ 
  double gain[COL_NUM];            /* lock loop parameter 1/tang*/
  double sq1gainScale[COL_NUM];    /* */
  int    zfact[COL_NUM];   /* lock loop parameter */
  int    ssabiaslckPt[COL_NUM];  /*  only record ssa bias lock point */
  int    biaslckPt[COL_NUM];  /*  bias lock point */
  int    biasrefPt[COL_NUM];  /*  bias reference point, only for sq1 */
  int    initFB[COL_NUM]; /* initial lock point for sq2=sa_fb  sq1=sq2_fb*/
  int    lastinitFB[COL_NUM]; /* upadted initial lock point for sq1=sq2_fb ar each bias*/
  int    slopSelect[COL_NUM]; /*  select positive (1) or negative (-1) slop */
  int    colMask[COL_NUM]; /*  aply 0 to mask=0, others apply biasVal */
  int    rowMask[ROW_NUM]; /*  aply 0 to mask=0, others apply biasVal */
  int    biasoptPt[ROW_NUM];  /*  sq1 bias optimal point */
  int    saoutlckVal[COL_NUM]; /* ssa DA out Value at lock point */
  int    cableadjVal[COL_NUM]; /* used for next bias cable offset, initial =0 or =cableadjInit */
  int    cableadjThd[COL_NUM]; /* used for next bias cable offset,  */
  int    cableadjInit[COL_NUM]; /* used for cable correction  */
  int    fluxPeriod[COL_NUM];  /* used for flux jump in sq1 */
  int    nomodulCount[COL_NUM];  /* used for counting no-modulation for all bias */
  int    sq1biasOpt[ROW_NUM];  /*  sq1 bias optimal point */
  int    sq2fdbkOpt[COL_NUM];  /* used for sq2fdbk optimal Val in sq1 */
  double cableadjScale[COL_NUM]; /* cable adjust scale factor */ 
  int    sampleNo;     /* record the MCE sample_num for re_scale data */
  int    dataMode;     /* record the MCE data_mode for processing data */
  int    *allData;     /*  pointer used for all channale data perbias */
  double *chData;      /*  pointer used for each channale data perbias */
  char   *cableadjPtr; /* pointer used for all cableadjVal (TOTACHS*biasNo) */
  char   *meanvalPtr;  /* pointer used for all meanVal (TOTACHS*biasNo) */
  char   *servodataPtr;/* pointer used for all servo data (TOTACHS*biasNo*fdbkNo) */
  char   *filtedPtr;   /* pointer used for all filteded data (TOTACHS*biasNo*fdbkNo) */
  double *convolPtr;   /* pointer used for convolution data (FILT_ORDER) */
  double *impulsePtr;  /* pointer used for filter impulse data (FILT_ORDER) */
  int    *heater;       // address holding heater value for modulation   
  int    *transitFlag;  // address holding flag indicate if in transition or not   
  double *transitVal;   // address holing Val for chk if it is in transition
  int    *lockFlag;     //  addres holding flag indicate if in lock or not
  CHKPOINT  minValue;          /* minimum point   */
  CHKPOINT  maxValue;          /*maximum point  */
  CHKPOINT  refValue;          /* ref maximum point  */
  PIXEL     pixeLock;         /* sq1open stage ssalock, sq1fb and gain */
} ARRAYSET;

// structure for servo Info
typedef struct 
{
int  totalptsSize;     //= setup->biasNo*setup->fdbkNo*sizeof(int);
int  headSize;         //= sizeof(SERVO_DATAHEAD); 
int  totalframeSize;   //= totalptsSize*TOTALCH;
int  peakSize;         //= sizeof(BIAS_PEAK)*setup->biasNo;
int  maxpeakSize;      //= sizeof(MAX_P2P)*TOTALCH
int  totalservoSize;   //= headSize + frameSize + peakSize + maxpeakSize;
int  totaldataperBias; //=COL_NUM*setup->fdbkNo;
}  SERVO_INFO;

#endif

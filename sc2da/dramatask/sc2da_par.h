/**
 * \file sc2da_par.h
 *
 * \brief parameters define for SCUBA2 Data Acquisition Software (DAS)  
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
 *  $Log: sc2da_par.h,v $
 *  Revision 1.43  2011/05/11 21:58:54  cwalther
 *  Added ability to use the median in heater tracking
 *
 *  Revision 1.42  2011/04/22 00:25:47  cwalther
 *  I put a considerable number of comments in this code and it has the changes made to normalize thresholds when the number of samples changes
 *
 *  Revision 1.41  2011/03/01 23:47:27  cwalther
 *  Added test for maximum gaini, set to zero if greater than maximum
 *
 *  Revision 1.40  2010/11/16 03:10:04  cwalther
 *  Added a simulate heater tracking bit to the heater tracking style
 *
 *  Revision 1.39  2010/09/30 20:09:42  cwalther
 *  Changes made while for engineering with science mode
 *
 *  Revision 1.38  2010/08/03 22:16:10  cwalther
 *  Changes to make sc2_setup work with SC2SCRATCH
 *
 *  Revision 1.37  2010/03/04 01:32:54  cwalther
 *  Changes for heater tracking memory and setting back to defaults in the dark
 *
 *  Revision 1.36  2010/02/19 00:29:10  cwalther
 *  Changes to make fast flat fields work
 *
 *  Revision 1.35  2010/01/26 22:24:42  cwalther
 *  corrected count comments in long enum structure
 *
 *  Revision 1.34  2009/12/15 23:31:17  cwalther
 *  Added more detailed status for malloc failures
 *
 *  Revision 1.33  2009/11/24 22:02:03  cwalther
 *  Code for displaying the counters for each STATE data event comming in
 *
 *  Revision 1.32  2009/11/18 22:22:26  cwalther
 *  Updated the comments so they counted the enum correctly - removed the unused (and commented out) modes
 *
 *  Revision 1.31  2009/10/21 01:48:40  xgao
 *  add changes for sq1lock and fluxperiod less than one ( mainly SQ2), 1p5
 *
 *  Revision 1.30  2009/09/22 00:04:22  cwalther
 *  Changes to flip rows and columns
 *
 *  Revision 1.29  2009/08/11 15:10:27  xgao
 *  update to UPDATE_SQ2FB SQ1LOCK and tidyup scripts
 *
 *  Revision 1.28  2009/01/27 03:19:14  cwalther
 *  Removed all references to modes other than SCAN, STARE and DREAM. Also started publising QL filenames for DARKs during SCANs
 *
 *  Revision 1.27  2009/01/07 19:39:16  cwalther
 *  Modifications for reading in the Integrator and proportional gains
 *
 *  Revision 1.26  2008/12/31 00:08:16  cwalther
 *  Changes to make Gao's mutipixel heater tracking work
 *
 *  Revision 1.25  2008/12/12 22:27:19  xgao
 *  sq2fb servo test
 *
 *  Revision 1.24  2008/08/14 02:56:09  cwalther
 *  These changes from Gao snuck in but I am going to commit them as they seem to do no harm
 *
 *  Revision 1.23  2008/08/07 01:49:35  cwalther
 *  Lots for changes for headers from Craig. Changes in handling events from Gao
 *
 *  Revision 1.22  2008/08/01 15:39:03  xgao
 *  add PIXELHEATER_SLOPE
 *
 *  Revision 1.21  2008/06/28 20:11:54  xgao
 *
 *  during initialise, call jitPathGet if SC2ENGMODE ==0 and 1. otherwise,don't call
 *  check if USEDV==1,    seq: before GO cmd, send dvCmd to MCE
 *                    trkheat: before GO cmd, send notdvCmd to MCE
 *                      servo: before GO cmd, send notdvCmd to MCE
 *  seq:    move DitsTriger after sc2headman_seqStart
 *  trkheat move DitsTriger after receiving first data frame
 *  add pixelMask action, after task loaded, pixelMask[]=1;
 *                        after pixelmask cmd, pixelMask[] is updated by pixelmask.xml
 *                        pixelMak[] is applied to I-Val for MCE
 *  take HELP off as jit will handle it
 *
 *  Revision 1.20  2008/06/05 00:43:43  cwalther
 *  Changes to get things working at JAC
 *
 *  Revision 1.19  2008/05/14 15:10:54  xgao
 *  update version as 0p48 after windowed DREAM
 *
 *  Revision 1.18  2008/05/02 13:37:45  xgao
 *   put right das version
 *
 *  Revision 1.17  2008/05/02 13:26:26  xgao
 *
 *  tested STARE/SCAN/DREAM without restart DAS. for DREAM, reconstruct Iamge
 *  from frames (don't call _flatten)  0p47
 *
 *  Revision 1.16  2008/04/24 15:12:02  xgao
 *  add Dream bef-BDK, tidy sc2da_struct.h
 *
 *  Revision 1.15  2008/03/28 16:40:08  xgao
 *  update for sc2headman 0p45
 *
 *  Revision 1.14  2008/03/07 15:30:04  xgao
 *  tidy up after doxygen 0p44, and merge from atc and cardiff
 *
 *  Revision 1.13  2008/02/20 15:08:42  xgao
 *  JAC AT tested STARE
 *
 *  Revision 1.12  2007/12/07 16:01:19  xgao
 *
 *  some modifications after SG array test, shutter only third open, new
 *  initialise.xml ( flatfile, INITFLAG) and configure.xml (DATA_MODE).
 *  negative pixHeat for loadv in SETU_SEQ. Add findheaterspike for processing
 *  heaterservo data for G/Tc.
 *
 *  Revision 1.11  2007/11/02 12:14:08  xgao
 *
 *  update to 0p33-021107, worked with JOC for INIT,CONFIG, SETUPSEQ, SEQ
 *
 *  Revision 1.10  2007/09/05 16:33:48  xgao
 *
 *  add heaterslope, trkHeater tested Ok, sq1Refbias and sq2Reffb always use the last one
 *  always do manual adjustment after optimal, fits4all is now in parShm, bigPhy=512M now
 *  0p21, ndf is turned on, tested with startmce-rtsc flatfield-seg
 *
 *  Revision 1.9  2007/08/10 16:38:54  xgao
 *
 *  minor change from 0p16-nondf, cablecorrect updates ssabiaslock
 *  and is controlled by slopSletect[9]. 0p17-nondf
 *
 *  Revision 1.8  2007/08/06 15:57:40  xgao
 *  re-organise servo setupfile, use include=. 0p16-nondf
 *
 *  Revision 1.7  2007/08/02 21:09:40  xgao
 *  fix three bugs after autoset run,0p15-nondf
 *
 *  Revision 1.6  2007/07/23 21:10:46  xgao
 *  second update to match ICD
 *
 *  Revision 1.5  2007/07/09 16:03:47  xgao
 *  consider eng and rtsc both firstStep
 *
 *  Revision 1.4  2007/06/18 23:20:08  xgao
 *  add functions for SCRIPTDAS
 *
 *  Revision 1.3  2007/06/14 17:58:17  xgao
 *  add filtedfunc for tesTransit, tidyup
 *
 *  Revision 1.2  2007/05/23 18:08:19  xgao
 *  use JCMTState, add a few actionstate. xg
 *
 *  Revision 1.1.1.1  2007/05/16 08:27:01  dkelly
 *  first insertion
 *
 */


// This is the message buffer size, see DitsInit.
 
#define DBUFSIZE 2000000

/* Default message buffer sizes.
 * See DitsAppInit and DitsPathGet routine descriptions.
 */
#define MESSAGEBYTES 400
#define MAXMESSAGES 5
#define REPLYBYTES 800
#define MAXREPLIES 12
#define TIMEOUT 20             /* Timeout for message operations */
#define NANODELAY  1000000000L /* delay for trkheater  nanoSecond*/


//drama task name
#define DASTASKNAME "SC2DA"                     
#define FDHTASKNAME "sc2daFDHDRAMA"                     
#define RTSCTASKNAME "sc2RTSCDRAMA"                     

//parameter names
#define SC2DALOGFILE "LOG_FILE"              // Name for log file
#define SC2DAMCEXMLFILE "MCE_XML_FILE"       // Name for mce.xml 
#define SC2DADATAFORMAT "DATA_FORMAT"        // Name for data format for writing 

#define DAS_VERSION "1p6-20091208"

// for shared memeory
#define OBS_SHARED1  "/etc/passwd"
#define OBS_SHARED2   'N'

#define PAR_SHARED1  "/etc/passwd"
#define PAR_SHARED2   'S'

// used in mce.xml for recording external setting
// like shutter open,ec
#define EXTERNAL_CARD "header"

// maximum sub-array 
#define MAX_SUBARRAY  8

// for internal message buffer re-cycle
#define  MSGRECYC_NO   16
#define  MSGRECYC_HEX  0x10000

//mce stuff
#define MCEBLK_SIZE   256
#define CARDID_LEN    20
#define SSAPARAMNO  10              // for ssa servo

// for setup the databuffer,see sdsudriver_par.h
// #define FRAMEHEADER_NUM  43
#define HEADERS       FRAMEHEADER_NUM
#define FRAME_SIZE    DAS_FRAME_SIZE
#define FRAMESIZE_INT COL_NUM*ROW_NUM+FRAMEHEADER_NUM+CHKSUM_NUM


#define MAXFILESIZE 0x100000000

// define which shared memory
typedef enum 
{
  SHAREDM_PAR =0,  
  SHAREDM_OBS
} SHAREM_TYPE;

// define which shared memory
typedef enum 
{
  INT_NUMBER =0,  
  FLOAT_NUMBER
} NUMBER_TYPE;

typedef enum 
{
  CHKXMEM =0,  
  CHKYMEM,
  CHKPMEM,
  CHKFIBRE
} PCI_HEALTHCHK;

typedef enum 
{
  USE_MSGOUT =0,  
  USE_ERSOUT,
  USE_ERSREP,
  USE_PRINT,
  NO_PRINT
} PRINT_MSG;


typedef enum 
{
   ANY     =0x00,    //any Card
   PSC     =0x01,    //Power Supplier Card
   CC      =0x02,    //Clock Card
   RC1     =0x03,    //Readout Card1
   RC2     =0x04,    //Readout Card2
   RC3     =0x05,    //Readout Card3
   RC4     =0x06,    //Readout Card4
   BC2     =0x07,    //Bias  Card2
   BC3     =0x08,    //Bias  Card3
   BC1     =0x09,    //Bias  Card1
   AC      =0x0A,    //Address Card
   RCS     =0x0B,
   SYS     =0x0E,    //System  
}  MCE_CARDID;

/* Please Note! The count in this enum must match 
   the count in the array of strings called seqStatus 
   found in sc2da_actionstate.h */

// PROCESS  status, the first =0, others are numbered sequentially
typedef enum 
{
  /// frameSetup
  // from 0
  DA_MCE_NONE=0,         ///< DA not in any data taking mode
  DA_MCE_SEQ,             ///< DA in sequence frame data taking mode   
  DA_MCE_SINGLEDATA,      ///< DA in single data taking mode
  DA_MCE_TRKHEATER,       ///< DA in tracking heater data taking mode   
  DA_MCE_FPGA,            ///< DA in donloadFPGA data taking mode   
  DA_MCE_AUTOSET,         ///< DA in AUTOSET ARRARY data taking mode",
  DA_MCE_DREAM,           ///< DA in DREAM data taking mode   
  DA_MCE_SSASET,          ///< DA in SSASET ARRARY data taking mode   
  DA_MCE_SQ2SET,          ///< DA in SQ2SET ARRARY data taking mode   
  DA_MCE_SQ1SET,          ///< DA in SQ1SET ARRARY data taking mode

  //SEQUENCE
  // from 10
  
  SEQ_NOACTION,        ///< no sequence action
  SEQ_INACTION,       ///=1,< still in sequence action
  SEQ_FINISHED,       ///=2,< sequence action finished OK  
  SEQ_STOPPED,         ///=3,< sequence action finished by kick 
  SEQ_ERROR,           ///< sequence action finished with Error
 
  //actionFlag
  // from 15
  NONEACTION,        ///< none action ever called
  ABORTACTION,       ///< abort action called
  BATCHACTION,       ///< batch action called
  BACWRACTION,       ///< bac write/read action called
  CONFIGACTION,      ///< configuare action called
  CHKSUMACTION,      ///< test checksum action called
  VERSIONACTION,      ///< check version action called
  DEBUGACTION,       ///< debug action called
  DISPACTION,        ///< dispInfo action called
  DSPACTION,         ///< downloadDSP action called
  EXITACTION,        ///< exit action called
  HELPACTION,        ///< help action called
  FPGAACTION,        ///< downloadFPGA action called
  INITACTION,        ///< initial action called
  MCECMDACTION,      ///< mcecmd action called
  MCESTATUSACTION,    ///< mcestatus action called
  PCICMDACTION,      ///< pcicmd action called
  PINGACTION,        ///< ping  action called
  SEQACTION,         ///< sequence action called
  SETSEQACTION,      ///< setup sequence action called
  SEQKICK,           ///< sequence kick action called
  SERVOACTION,       ///< servo action called
  SERVOKICK,          ///<  servo kick action called
  TESTACTION,        ///< tets action called
  TRKHEATACTION,     ///< track heat action called
  TRKKICK,           ///< track kick action called
  UPDATEACTION,      ///< updatelogname action called
  FDH_START,          ///< in action for FDH_START
  FDH_KICK,           ///< stoped by  kick  
  PIXELMONACTION,     ///< monitor pixel action called     
  PIXELMONKICK,        ///<monitor pixel kick action called     
  SCRIPTDASACTION,     ///<scriptdas action called
  SCRIPTDASKICK,       ///<scriptdas kick action called
  HEATERSLOPE,         ///<heaterslope action called
  SETENGDATA,          ///<set data format and save-engData action called
  TRKSQ2FBACTION,      ///<track sq2fb action called
  SQ1OPTACTION,        ///<find sq1 optimal points
  GETSQ2FBPARAACTION,  ///< get sq2fb para for doing sq2fb track in SEQ
 
  //batchFlag
  //from 53
  BATCH_GO,           ///< GO from batch file
  BATCH_CMD,          ///< nomal cmd from batch file
  BATCH_CMNT,         ///< comment from batch file
  BATCH_END,          ///< no more cmd in batch file
  BATCH_CMDFAIL,      ///< send cmd from batch file failed
  BATCH_TRANSFAIL,    ///< translate cmd from batch file failed


  /// frame finishing reason
  // from 59
  FRAME_WAITTIMEOUT,      ///< wait for frame timeout 
  FRAME_WAITING,          ///< idle or wait for first frame 
  FRAME_COMPLETION,       ///< frame reading is finished
  FRAME_ERRGETPASSBUF,    ///< frame finished with error in sdsu_getpass buf    
  FRAME_FIRST,            ///< first frame 
  FRAME_MISSING,          ///< frame finished with frame missing
  FRAME_STOPPED,          ///< frame finished by stop/kick 
  FRAME_MORE,             ///< frame finished with more frame 
  FRAME_CHKSUMWRONG,      ///< frmae checksum wrong    
  FRAME_TRIGQL,           ///< coadd frame is ready for QL    
  FRAME_FILEFAILED,       ///< open data file failed, check permision setting
  FRAME_MALLOCFAILED,      ///< malloc for holding frame datafailed
  SHAREMEM_OBSFAILED,      ///< fail to connect shared Memory 
  FRAME_SCANREADY,        ///< scan is ready for QL if not dark
  FRAME_SQ2FBCMP,         ///< sq2fb compansation
  FRAME_MALLOCFAILED1,    ///< malloc for  holding windowed data
  FRAME_MALLOCFAILED2,    ///< malloc for  holding  temp dream
  FRAME_MALLOCFAILED3,    ///< malloc for  holding re-organised frame
  FRAME_MALLOCFAILED4,    ///< malloc for  holding  ch data
  FRAME_MALLOCFAILED5,    ///< malloc for  holding  frame status
  FRAME_MALLOCFAILED6,    ///< malloc for  holding   dark row
  FRAME_SAW_RAMP,         ///< Step heater current or bias for all sawtooths and ramps

  /// OBS mode
  // from 81
  OBS_STARE,              ///< STARE mode
  OBS_DREAM,              ///< DREAM mode
  OBS_SCAN,               ///< SCAN mode

  /// Msg FIFO error status
  // from 84
  ERR_NONE,               ///< idle or non error in put msg into msg FIFO
  ERR_SDSU_WRITMSG,       ///< writing USER_USER_MSG FIFO error
  ERR_SDSU_WRITFRAME,     ///< writing USER_USER_FRAME FIFO error    
  ERR_SDSU_WRITFERR,      ///< writing USER_USER_ERR FIFO error    
 
  /// operation on FIFOs
  // from 88
  FIFO_ERR_NONE,          ///< idle or non error in writing/opening FIFO
  FIFO_ERR_EM2EM,         ///< error in opening USER_USER_EM USER_DATA_EM FIFO
  FIFO_ERR_DRAMMSG,       ///< error in opening USER_USER_MSG FIFO
  FIFO_ERR_FRAM,          ///< error in opening USER_USER_FRAME FIFO
  FIFO_ERR_ERR,           ///< error in opening USER_USER_ERR FIFO

  ///whereabout the process
  // from 93
  WAIT_OBEY,      ///< wait for obey command
  Dits_GetSeq0,   ///< wait in DitGetSeq()=0 
  Dits_GetSeqn,   ///< wait in DitGetSeq()=n 
  WAIT_sem,       ///< wait for intask semaphore
  WAIT_emptybuf,  ///< wait for getpass_buf,ie,bufAdr in USER_USER_EMP

   ///  setup array
   // from 98
  SSASETUP,                    ///<it is for SSA",
  SQ2SETUP,                    ///<it is for SQ2",
  SQ1SETUP,                    ///<it is for SQ1",
  DO_MODULATION,          ///<it is for modulation",
  DO_DERIVATIVE,          ///<it is for modulation",
  DO_BIASRAMP,            ///<it is for bias ramp function",
  NONE_BIASRAMP,          ///<non bias ramp function",
  FOUND_MIN,              ///<found minimum modulation", 
  FOUND_MAX,              ///<found maximum modulation",

  // load 
  // from 107
  LOAD_DARK,    // shutter closed
  LOAD_HOT,     // shutter open to hot source
  LOAD_SKY,     // shutter open to sky
}  DA_PROCESS;   

typedef enum 
{
   DAS_NO_ACTION     =0x00,
   DAS_MCE_CONFG     =0x10,   
   DAS_MCE_CMD       =0x20,   
   DAS_MCE_DWLDFPGA  =0x30,   
   DAS_MCE_EXIT      =0x40,   
   DAS_MCE_GO        =0x50,   
   DAS_MCE_STOP      =0x60,   
   DAS_MCE_TEST      =0x70,   
   DAS_PCI_CMD       =0x80,   
   DAS_PCI_DWLDDSP   =0x90,
   DAS_MCE_WRBLK     =0xA0,
}  DAS_CMD;

// data save format
typedef enum 
{
   MCE_BINARY_FORM=0, //< save in binary format
   MCE_TEXT_FORM,     //< save in text format 8 data in a row
   MCE_TEXT_FORM2,    //< save in text format 1 data in a row
   DATFILENOAPPEND,   // data file no appending
   DATFILEAPPEND,     //appending for data file
   NODATAFILE,
   NOBATCHFILE,
   BATCHFILE,
   OTHERFILE
}  DATA_FORM;


/// command status
typedef enum 
{
   MCECMD_OK=0, 
   MCEXML_ERROR,           
   MCE_TIMEOUT,          
}  CMD_STATUS;

/// date time
typedef enum 
{
   DAS_DATE =0, 
   DAS_DATETIME,          
}  DATETIME_STATUS;

/// flag in savedata function 
typedef enum 
{
   NOHEADER =0, 
   HASHEADER,          
}  SAVEDATA_FLAG;

/// flag in setup-seq 
typedef enum 
{
   RTSC_MODE =0, 
   ENG_MODE,          
}  DAS_MODE;


/// flag in fillNDF 
typedef enum 
{
   AT_END_SEQ =0, 
   SEQ_NOTEND,          
}  WHARE_IN_SEQ;




// for array setup

// FITS header for QL
#define FITSSIZE 80
// FITS header for SEQ, only for test
#define FITSALLSIZE 81

#define MAXFITS 1000  // main fits head
#define MAXQLFITS 50  // sub fits head

// set to small one for DREAM, get real one 
// from reading weightsfile in DH-task
#define SMU_PATTERN   64

typedef enum 
{
 QL_WAITING =0x54494157,  //'WAIT'  big/little endian, 
 QL_READY   =0x59444552,  //'REDY',
 QL_DONE    =0x454E4F44,  //'DONE',
} QL_LOOKUP;


// from rtsDClient
// Number of tasks to be monitored 
#define SC2_MONITOR_TASKS 7  
#define MAX_INPUT_STR 15

// Constants for handling monitoring 
// from Craig's monTEst.c 
typedef enum 
{
  WAITING_FOR_ID          =1,
  WAITING_FOR_PUBLISHED_ID=2,
  ID_PUBLISHED            =3,
  TASK_NOT_DEFINED        =4,
} MONITOR_STATUS;

/// servo
typedef enum 
{
NONESERVO=0,
CABLECAL, 
CABLECORET,
HEATERSERVO,
SQ1BIASSERVO,
SQ1OPEN,
SQ1SERVO,
SQ2LOCK,
SQ2OPEN,
SQ2SERVO,
SSALOCK,
SSARAMP,
SSARAMP1,
TESBIASSERVO,
TESTRANSIT,
PIXELMON,
TRKHEATER,
TRKHEATERINIT,
SQ2FBTRACK,
SQ2OPEN4P,
SQ2BIASING,  
SQ1LOCK,
BIASCHANGED,
BIASNOCHANGE,
SQ2SERVOSSALOCK,
BLACKBODY,
HEAT_SLOPE,
SQ2FBGETPARA,
TRI_UP,
TRI_DOWN,
TRI_MID,
}  SERVO_STATUS;


// apply ssaFb in sq2servo or sq2Fb in sq1servo
typedef enum 
{
   NOT_APPLY_FB=0, 
   APPLY_FB,
} SSASQ2_FB;
 

#define COADDFILE_DAS "coadd-das.hex"
#define COADDFILE_DH  "coadd-dh.hex"
#define CHFILE_DH     "chdata-dh.hex"
#define FSTATFILE_DH  "fstat-dh.hex"
#define DARKFILE_DH   "dark-dh.hex"


// not modulation value
#define NONMODULATION   0  // 0xffffffff // VAL__BADI
#define NONMOD_BIAS     0  // 0xffffffff  // VAL__BADI

// search period for finding peak inside the channel data
#define SEARCHPERD  30
#define SEARCHPERD_SQ2  200
#define SEARCHPERD_SQ1  200
 
#define HEADERSET  HEADERS*4    // shift away from the header byte to real data 
#define TOTALROWS  ROW_NUM      // include the drak row
#define TOTALCHS   COL_NUM      // chs for four cards
#define ONECARDCHS  8          // chs per card


// for arrary setup use

#define  SERV_STR_LEN   100  // string length in servo data buffer head for 
                             // newer than 27/11/06

// settle down time (second) after ssa bias or sq2 bias
#define  WAIT_TIME_AFTER_BIAS  1


// define NO of param for the servo setup file
#define OTHERS_N    4
#define CARDS_N     4
#define BIAS_N      3
#define FEEDBK_N    3
#define CABLE_N     3
#define COLMASK_N 4
#define ROWMASK_N 1
#define BIASOPT_N 1
#define SLOPSELECT_N 4
#define CABLEOFFSET_N 4
#define CABLESCALE_N 4
#define CABLESLOPE_N 4
#define CABLEADJ_THD_N 4
#define CABLEADJ_SCALE_N 4
#define CABLEADJ_INIT_N 4
#define GAIN_N      4
#define ZFACTOR_N   4
#define BIASLCK_N   4
#define SSABIASLCK_N   4
#define INITFB_N    4
#define GAINSCALE_N 4
#define FLUXPERIOD_N 4
#define SQ2FDBKOPT_N 4

#define SETUPPAR          OTHERS_N  + CARDS_N + BIAS_N + FEEDBK_N + SLOPSELECT_N + COLMASK_N 
#define CABLEPAR1         CABLEOFFSET_N + CABLESLOPE_N + CABLESCALE_N + CABLEADJ_THD_N +CABLEADJ_SCALE_N

#define SETUPPAR_CABLE    SETUPPAR  + CABLE_N

#define SETUPPAR_SSARAMP1   SETUPPAR + CABLESLOPE_N +CABLEOFFSET_N 
#define SETUPPAR_SSARAMP    SETUPPAR + CABLEPAR1

#define SETUPPAR_SQ2FBTRACK OTHERS_N + CARDS_N + INITFB_N + COLMASK_N + GAIN_N + FLUXPERIOD_N + SLOPSELECT_N + GAINSCALE_N +ZFACTOR_N 

#define SETUPPAR_SQ2OPEN    SETUPPAR + INITFB_N + BIASLCK_N
#define SETUPPAR_SQ2LOCK    SETUPPAR + GAIN_N + ZFACTOR_N + INITFB_N + BIASLCK_N
#define SETUPPAR_SQ2OPEN4P  SETUPPAR + BIASOPT_N + SQ2FDBKOPT_N

#define SETUPPAR_SQ1OPEN    SETUPPAR + ROWMASK_N + ZFACTOR_N + BIASOPT_N + SQ2FDBKOPT_N
#define SETUPPAR_SQ1LOCK    SETUPPAR + ROWMASK_N + INITFB_N + BIASOPT_N

#define SETUPPAR_SQ2        SETUPPAR + GAIN_N + ZFACTOR_N + INITFB_N + BIASLCK_N
#define SETUPPAR_SQ1        SETUPPAR + ROWMASK_N + GAIN_N + ZFACTOR_N + INITFB_N + BIASLCK_N + GAINSCALE_N +FLUXPERIOD_N
#define SETUPPAR_HEAT       SETUPPAR + ROWMASK_N + GAIN_N + ZFACTOR_N + INITFB_N + BIASLCK_N + GAINSCALE_N +FLUXPERIOD_N
#define SETUPPAR_TES        SETUPPAR + ROWMASK_N + GAIN_N + ZFACTOR_N + INITFB_N + BIASLCK_N + GAINSCALE_N +FLUXPERIOD_N
#define PIXELMON_PAR        OTHERS_N +  SLOPSELECT_N + COLMASK_N + ROWMASK_N
#define SETUPPAR_SQ1BIAS    SETUPPAR + ROWMASK_N + GAIN_N + ZFACTOR_N + INITFB_N + BIASLCK_N+ GAINSCALE_N +FLUXPERIOD_N
#define SETUPPAR_TRANSIT    SETUPPAR

#define SETUPPAR_CABLECORET SETUPPAR + CABLEPAR1 +BIASLCK_N + CABLEADJ_INIT_N 
#define SETUPPAR_SSALOCK    SETUPPAR + CABLEPAR1 +BIASLCK_N + CABLEADJ_INIT_N 

// for findlockpint routine, include 
// servoName, whichRC[], bias, bstep, nbias, feed, fstep, nfeed,  ==>8
// whichRow, doServo, safbcard, sq2fbcard, sq2biascard, sq1biascard, ==>6
#define  SERVO_SETUP_NO1 14
// initFB[], cableOffset[], biasLck[], slopselect[]   ==>4
// cableSlope[], cableScale[], colMask[], rowMask[] ==>4
#define SERVO_SETUP_NO  SERVO_SETUP_NO1+8

#define INITLCK1_N  4
#define INITLCK2_N  4

#define MINVALUE    10
// the sliding window for check min and max value of SSAOUT 
#define  CHKNUM     50       

// the step advance for reference Icmax (=(1+REF_STEP)*Icmax)
#define  REF_STEP   0.2       

// define NO_SAMPLES used in MCE firmware 
#define  MCE_SAMPLE_NO   48

// define the SSA out bottom line value
//#define SSA_BOTTOM_LINE  8192*MCE_SAMPLE_NO*0.9 
// change it to Zero line as the bottom line
#define SSA_BOTTOM_LINE  0

// for lowpass filter
#define FILTER_ORDER      31
#define FC                10
#define SF                301


// all setup files are in $CONFIG_HARD or $CONFIG_ALL
// FB:      initFB,  
// ZG:      Zfactor and Gain
//FLUXPERD: flux period   
#define CABLE_CAL          "setup-cablecal"
#define CABLE_CORRECT      "setup-cablecorrect"
#define SSA_RAMP1          "setup-ssaramp1"
#define SSA_RAMP           "setup-ssaramp"
#define SSALCK             "setup-ssalock"

#define SQ2SERVOFILE       "setup-sq2servo"
#define SQ2OPENFILE        "setup-sq2open"
#define SQ2LOCKFILE        "setup-sq2lock"

#define SQ1SERVOFILE       "setup-sq1servo"
#define SQ1OPENFILE        "setup-sq1open"

#define HEATSERVOFILE      "setup-heaterservo"
#define TESBIASFILE        "setup-tesbiasservo"
#define SQ1BIASFILE        "setup-sq1biasservo"
#define TESTRANSITFILE     "setup-testransit"

// mid setup files
#define CABLE_OFFSETSLOPE  "cable-offsetSlope"
#define CABLE_SCALE        "cable-Scale"
#define CABLE_ADJINIT      "cable-adjInit"
#define SSALOCK_BIAS_FB_ZG "ssalock-bias-FB-ZG"   
#define SSALOCK_FB_ZG      "ssalock-FB-ZG"   
#define SSALOCK_BIAS       "ssalock-bias" 
#define SSALOCK_SSABIAS    "ssalock-ssabias" 
#define SSALOCK_BIAS_CABLEOFFSET "ssalock-bias-cableOffset" 
#define SQ2LOCK_BIAS       "sq2lock-bias"     
#define SQ2LOCK_BIAS_FB    "sq2lock-bias-FB"     
#define SQ2LOCK_BIAS_FB_ZG "sq2lock-bias-FB-ZG"
#define SQ2LOCK_FLUXPERD    "sq2lock-fluxperd"          
#define SQ2OPEN_BIAS_FB_ZG  "sq2open-bias-FB-ZG"         
#define SQ2OPEN_Z           "sq2open-Z"
#define SQ2OPEN4P_FB_ZG     "sq2open4p-FB-ZG"         

#define OPTIMAL_SQ1B_SQ2FB  "optimal-sq1bias-sq2fb"

#define SSABIASLCK         "ssabiaslock.txt"
#define SSAFBLCK           "ssafblock.txt"
#define SSALCKDIVBYN       "ssalockdivbyn.txt"

#define SQ2BIASLCK         "sq2biaslock.txt"
#define SQ2FBLCK           "sq2fblock.txt"
#define SQ1FBLOCKFILE      "sq1fblock.txt"

#define SQ1OPENTEMPLATE    "sq1open-template.txt"
#define SQ1OPENTEMPLATE1   "sq1open-template1.txt"

#define SQ1OPTPTS          "sq1opt-points"

#define SSAOUTMIDVAL       "ssaoutmidVal"
#define MCEGAIN            "mce-gain-lock-sq1fb"
#define MCELOCKTEMPLATE    "mcelock-template1.txt"
#define MCELOCKSSA         "mcelock-ssa.txt"
#define MCELOCKSQ2         "mcelock-sq2.txt"
#define SQ1BIASOPTFILE     "sq1biasoptimal.txt"
#define SQ2FBOPTFILE       "sq2fboptimal.txt"

#define  BINDIV 60
#define  SQ2FDBK_WIDTH  15
#define  BINTOP        0.9*BINDIV
#define  BINBOT        0.1*BINDIV

// data_mode==4  18fb+14err?   14 bit or 16 bit? check with UBC
#define ERRMASK   0x03FFF
#define FDBKMASK  0x3FFFF
#define FDBKSHIFT 14 
#define SIGNBIT_14BIT     0x02000
#define DATAMASK_13BIT    0x01FFF
#define SIGNBIT_18BIT     0x20000
#define DATAMASK_17BIT    0x1FFFF
#define SIGN_EXT_14    0xFFFFE000 
#define SIGN_EXT_18    0xFFFE0000 

// for heater tracking
#define TRKHEAT_NO            5
#define HEATERSLOPE_NO        6
#define HEATER_TRACK         "setup-trckheater"
#define HEATER_SLOPE         "setup-heaterslope"
#define HEATSLOPE_RESULT     "heaterslopeTable"
#define PIXELHEATER_SLOPE     "heaterslope"
#define HEAT_TRACK_MAX        100
#define WALK_HEATER_CURRENT_BIT               1
#define REMEMBER_VALUE_IN_DARK                2
#define READ_VALUE_EVERY_TIME                 4
#define UPDATE_HEATER_TRACKING_REFERENCE_BIT  8
#define SIMULATE_HEATER_TRACKING_BIT          16

// trap sq1bias
#define    MAX_SQ1BIAS       16383

//pixelmon and sq2fb servo
#define SQ2FB_UPDATE_WAIT 100

FILE *myFpLog;

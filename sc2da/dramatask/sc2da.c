/**
 * \file sc2da.c
 *
 * \brief drama application for SCUBA2 Data Acquisition Software (DAS)  
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
 *sc2da_Abort
 *
 * The DRAMA task handles all DRAMA communications and actions between MCE and
 * RTS  Client. 
 *
 * The sc2da task includes flags indicating the current state of the 
 * application -is it initialised, has it been configured or  setup, is there a
 * sequence currently being executed. The values of these flags can modify the 
 * behaviour of the actions. Some of the actions are responsible for reading
 * and parsing XML files, and storing the results in data structures. Key
 * actions are  ABORT and EXIT, which corresponds to a KICK of the SEQUENCE 
 * action. These operations must be carried out immediately, especially when a
 * SEQUENCE is actually in progress.   
 *
 * All hardware related parameters and names are defined in a XML file supplied 
 * by MCE designers.  The sc2da task parses the XML files using the "expat" 
 * C library
 *
 * When started, it establishes connection to the shared memory  (SDSU_CONTEXT
 * structure) allocated by sdsu_driver module in  kernel space, obtains a 
 * pointer of the share memory  as *con and allocates  data buffers for storing 
 * the acquisition data.
 *
 * There are 15 FIFOs created by the sdsu_driver. Four of them:
 * USER_USER_MESG, DAS_FDH_MESG, DRV_USER_MESG and USER_DRV_MESG FIFOs
 * are used for passing message,test data and error from other process
 * or driver to drama task, and vice versa. the message will trig DRAMA action
 * to reschedule or driver or other process to do things
 *
 * During operation the DRAMA task receives setup information from the RTS 
 * Client, interprets it and puts it into  local shared variables. When a 
 * SEQUENCE command is received from the RTS Client, It sets up variables in the 
 * CON struture, allocates share memory for information between drama task and 
 * data handle task, sends a GO command to the MCE to start sequencing operations. 
 * After receiving reply for GO from  the MCE,  the SEQUENCE action first sends 
 * an acknowledge (DitsTrigger) to RTS Client, then reschedules to await for message 
 * from other process. either for qlook or completion.
 * 
 * The other process(sc2dadh task) first obtains a pointer from SDSU_CONTEXT structure:CON, 
 * attachs FIFOs, opens  SDSU_DEV_FILE  and establish a mapping between the process' user address 
 * space and another shared memory (big physical memory) represented by  the device file. it then 
 * continues by entering  an infinite loop. 
 * 
 *
 *  $Log: sc2da.c,v $
 *  Revision 1.99  2013/09/20 02:22:24  ryanb
 *  fixes for kick/abort sequence
 *
 *  Revision 1.98  2013/05/22 20:16:43  ryanb
 *  alternate ABORTACTION handling attempt
 *
 *  Revision 1.95  2012/10/19 00:56:46  cwalther
 *  changes to make sure all arrays have the exact same OBSID
 *
 *  Revision 1.94  2012/06/05 23:20:52  cwalther
 *  Turned the classic MCE error report red, added timing marks throughout INITIALISE they can be removed, they are marked with XXXX
 *
 *  Revision 1.93  2012/04/18 00:10:04  ryanb
 *  fixed heater tracking double-free/corruption bug
 *
 *  Revision 1.92  2012/03/28 02:39:38  cwalther
 *  Changes to handle INBEAM correctly (now shows shutter when in the dark)
 *
 *  Revision 1.91  2011/12/20 23:23:14  cwalther
 *  Implemented a style of heater tracking that ends after a fixed number of iterations
 *
 *  Revision 1.90  2011/09/27 21:14:22  cwalther
 *  Moved the heater tracking header write so that I would print when you track to the value you remember in the dark, some changes to make remembering in dark work better
 *
 *  Revision 1.89  2011/09/02 19:16:20  bgorges
 *  Standardized entry checks on all fuctions/actions that should quietly return if entry status is bad.
 *
 *  Revision 1.88  2011/08/31 18:38:28  cwalther
 *  added ability to stay alive after the classic MCE error
 *
 *  Revision 1.87  2011/06/02 21:31:22  bgorges
 *  Checking all of fopen and fopen64 to make sure a file is actually opened.
 *
 *  Revision 1.86  2011/05/20 21:57:31  cwalther
 *  Moved time loging at end-of-SEQUENCE out of errors only into always, added MsgOut and end of CONFIGURE
 *
 *  Revision 1.85  2011/05/17 20:18:29  cwalther
 *  These are all changes to cut down on the number of messages that have to go through the SCUBA2 task
 *
 *  Revision 1.84  2011/05/11 21:58:54  cwalther
 *  Added ability to use the median in heater tracking
 *
 *  Revision 1.83  2011/05/05 19:51:32  cwalther
 *  Fixed the zillion of messages problem when heater tracking fails
 *
 *  Revision 1.82  2011/04/22 00:25:47  cwalther
 *  I put a considerable number of comments in this code and it has the changes made to normalize thresholds when the number of samples changes
 *
 *  Revision 1.81  2011/04/16 02:36:06  cwalther
 *  Took hard return out of MsgOut becasue DRAMA adds one anyway
 *
 *  Revision 1.80  2010/12/04 00:29:36  cwalther
 *  Removing changing select_clk each time as it was crashing the MCE
 *
 *  Revision 1.79  2010/12/03 02:31:41  cwalther
 *  Took out some messages that were swamping the RTS Client with 8 arrays
 *
 *  Revision 1.78  2010/10/21 02:01:28  cwalther
 *  Removed optimal-sq1bias-sq2fb from CONFIG_HARD it is in scratch now
 *
 *  Revision 1.77  2010/10/07 23:32:14  cwalther
 *  Trying to make the log file more useable
 *
 *  Revision 1.76  2010/10/06 19:10:56  cwalther
 *  removing row_dly caused a loop to go through too many times
 *
 *  Revision 1.75  2010/10/06 02:04:45  cwalther
 *  Removed row_dly so array got one shorter
 *
 *  Revision 1.74  2010/10/06 00:00:35  cwalther
 *  Return errors on problems while heater tracking
 *
 *  Revision 1.73  2010/09/30 20:09:42  cwalther
 *  Changes made while for engineering with science mode
 *
 *  Revision 1.72  2010/09/23 18:55:56  cwalther
 *  Changes for reading mce state and writing it to a file during CONFIGURE
 *
 *  Revision 1.71  2010/08/31 02:40:37  cwalther
 *  Changes for handling the clock strickly as UBC has told us
 *
 *  Revision 1.70  2010/08/30 19:32:47  cwalther
 *  Changes for reading complete MCE status and writing it out as a stringified AST object
 *
 *  Revision 1.69  2010/08/03 22:16:10  cwalther
 *  Changes to make sc2_setup work with SC2SCRATCH
 *
 *  Revision 1.68  2010/07/19 20:00:19  cwalther
 *  Changes to make the clocks exclusively external or internal and not mixed
 *
 *  Revision 1.67  2010/04/19 22:03:30  cwalther
 *  changes to put dark heater current into FITS headers
 *
 *  Revision 1.66  2010/03/04 01:32:54  cwalther
 *  Changes for heater tracking memory and setting back to defaults in the dark
 *
 *  Revision 1.65  2010/02/26 18:27:01  cwalther
 *  Changed some MsgOuts into jitDebugs
 *
 *  Revision 1.64  2010/02/19 00:29:10  cwalther
 *  Changes to make fast flat fields work
 *
 *  Revision 1.63  2010/01/26 22:27:24  cwalther
 *  Changes for getting kicks of SEQUENCE to be benign - a fastflatfield comment
 *
 *  Revision 1.1.1.1  2007/05/16 08:26:57  dkelly
 *  first insertion
 *
 */


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


/* copy from rtsMaster
 */
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


#include "expat.h"       /* For the expat XML parser */

/*************** tmp out */ 
// for using sc2store see sc2store.c
#include "ast.h"
#include "star/hds.h"
#include "sae_par.h"
#include "prm_par.h"
#include "dat_par.h"
#include "ndf.h"
#include "mers.h"
#include "jcmt/state.h"
#include "star/hds_types.h"

// in INSTALL_INCLUDE 
// from dream and sc2store sc2headman
#include <dream_par.h>
#include <sc2store.h>
#include <sc2headman.h>
//include from scuba2Rtsc
#include <sc2rtsc_par.h>

/*** where it defines 
   #define SC2RTSC__DATUMSH      1 
   #define SC2RTSC__DATUMHOT     2
   #define SC2RTSC__CLOSESH      4
   #define SC2RTSC__HOTIN        8
   #define SC2RTSC__MCERESET     32
*******/


/* in {INCDIR} */  
#include <sdsudriver_err.h> 
#include <sdsudriver_errstring.h> 
#include <sdsudriver_par.h> 
#include <sdsudriver_struct.h>
#include <interface.h>
#include <mcexml_par.h>
#include <mcexml_struct.h>
#include <mcexml.h>

//
#include "sc2da_par.h"
#include "sc2da_actionstate.h"
#include "sc2da_struct.h"
#include "sc2da_err.h"
#include "sc2dalib.h"     
#include "sc2dalibsetup.h"      
#include "sc2dalibservo.h"           
#include "sc2daliblckpts.h"           
#include "sc2da.h"           

/* global struct info */
static SDSU_CONTEXT         *con;       
static dasInfoStruct_t      dasInfo;
static dasCmdInfo_t         onflyCmd;
static DRAMA_INNERMSG       glbMsg[16]; // recycle for 16
static struct mcexml_struct glbMceInx;
static char                 taskName[30]; 
static int                  pixelMask[41*32];
static int                  heaterMask[40*32];
static double               heaterSlope[40*32];
static char                 dvCmdEx[]="wb cc use_dv 2";   
static char                 dvCmdInt[]="wb cc use_dv 0";
static char                 syncCmdEx[]="wb cc use_sync 2";      
static char                 syncCmdInt[]="wb cc use_sync 0";      
// static char                 selClkCmdEx[]="wb cc select_clk 1";      
// static char                 selClkCmdInt[]="wb cc select_clk 0";      
static char                 domainName[30];
static ARRAYSET             setup0;


// define options for sq1
// if defined, servo routine will find the optimal points 
// if not defined, use findlockpoints to find optimal points
#define LOOK4OPTIMUM

#define HAVESQ2SSALOCK


/*===============  see              ================ */
/*===============function call routines ============= */
#include "sc2dalib.c"
#include "sc2dalibsetup.c"
#include "sc2daliblckpts.c"
#include "sc2dalibservo.c"


/* Here are the meanings of the bits for debug levels

1     Print out error messages
2     print out general message
4     Print out on each successful exit of the action (like a trace)
8     print out special message
16    print out mesgg for servo
32    print out batch command
64    print out message for seq  
128   
*/


/**
 * \fn int main(int argc, char **argv)
 *
 * \brief the main drama application 
 *
 * \param argc     number of argument
 * \param **argv   argument structure pointer 
 * 
 * The function main first establishes connection to the shared memory 
 * (SDSU_CONTEXT structure) allocated by sdsu_driver module in kernel space 
 * (the driver has been loaded during system startup),  obtains a pointer of
 * the  share memory as con and allocates data buffers  for storing the
 * acquisition  data.  The semaphore for signalling the fsvTask is created. 
 *
 * The main function then declares the DRAMA actions and callbacks, and 
 * initialises DRAMA using the  DITS_M_X_COMPATIBLE flag. It opens the FIFO
 * for  communications  coming back from the fsvTask task and driver and gets 
 * descriptor for reading from the FIFO (con->fifo[USER_USER_MESG].fdr), 
 * con->fifo[DRV_USER_MESG].fdr). The input descriptors are added to the 
 * DRAMA  communications system specifying sc2da_Seqmsg(), drivermsg() as the
 * callback routines..
 *
 * The DRAMA message handling loop is then invoked using DitsAltInLoop, so
 * that messages coming through the FIFOs is waited for at the same time as 
 * DRAMA messages.
 *
 * use the arguments Number_of_Cols, Rows and default debuglevel here only for 
 * engineering software.
 */

/*+sc2da- main routine  extern int main
*/
JIT_MAIN (sc2da)
( 
int argc, 
char **argv
)
{ 
  //int              dramabufSize = DBUFSIZE;   // Message buffer size   
  long             Cols,Rows;              // temp used  for sc2da command as args 
  int              errNo=-1;
  DitsAltInType    ditsAltDesc;
  StatusType       mystatus;
  StatusType       *status = &mystatus;
  static USED4DITS used4Dits;
  char             NodeName[64];   /* test for IMP machine name */

  char * tmpName = (argc > 1 ? argv[1]: DASTASKNAME);

  *status = STATUS__OK;

  strcpy(taskName, tmpName);
  if ( argc > 2)
    strcpy(dasInfo.givenarrayName,argv[2]);
  else 
    strcpy(dasInfo.givenarrayName,"s8a");
 
  // ActionMap associates action names with routines.  ActionMapSize is the
  //  size of the array.
  static DitsActionDetailsType Actions[] =   
  {
    { sc2da_Abort,                   0, 0, 0, 0, 0, "ABORT"},
    { sc2da_Config,   sc2da_ConfigKick, 0, 0, 0, 0, "CONFIGURE"},
    { sc2da_Exit,                    0, 0, 0, 0, 0, "EXIT" },
    { sc2da_Init,       sc2da_initKick, 0, 0, 0, 0, "INITIALISE"},
    { sc2da_Ping,                    0, 0, 0, 0, 0, "PING" },
    { sc2da_SetSeq,   sc2da_SetSeqKick, 0, 0, 0, 0, "SETUP_SEQUENCE"},
    { sc2da_Seq,         sc2da_SeqKick, 0, 0, 0, 0, "SEQUENCE"},
    { sc2da_EndObs,      sc2da_ObsKick, 0, 0, 0, 0, "END_OBSERVATION"},
    
    { sc2da_Batch,                    0, 0, 0, 0, 0, "MCEBATCHGO"},  
    { sc2da_Downld2PCI,               0, 0, 0, 0, 0, "DWLOADDSP" },  // only PCI
    { sc2da_Down2FPGA,                0, 0, 0, 0, 0, "DWLOADFPGA"},
    { sc2da_Dispinfo,                 0, 0, 0, 0, 0, "DISPINFO" },  
    { sc2da_Mcecmd,                  0, 0, 0, 0, 0, "MCECMD"}, 
    { sc2da_PCIcmd,                   0, 0, 0, 0, 0, "PCICMD"},
    { sc2da_PCIblk,                   0, 0, 0, 0, 0, "PCIBLK"},
    { sc2da_UpdateLogname,            0, 0, 0, 0, 0, "UPDATELOGNAME" },
    { sc2da_Servo,      sc2da_servoKick, 0, 0, 0, 0, "SERVO"},
    { sc2da_Trkheat,      sc2da_trkKick, 0, 0, 0, 0, "HEATER_TRACK"},
    { sc2da_Trksq2fb,     sc2da_trkKick, 0, 0, 0, 0, "SQ2FB_TRACK"},
    { sc2da_Getsq2fbPara,             0, 0, 0, 0, 0, "SQ2FB_GETPARA"},
    { sc2da_Heatslope,                0, 0, 0, 0, 0, "HEAT_SLOPE"},
    { sc2da_Version,                  0,  0, 0, 0, 0,"DASVERSION"},
    { sc2da_MceStatus,                0,  0, 0, 0, 0,"MCESTATUS"},
    { sc2da_Mceonflycmd,              0,  0, 0, 0, 0,"MCEONFLYCMD"},
    { sc2da_Pixelmon,   sc2da_pixelKick, 0, 0, 0, 0, "PIXELMONITOR"},
    { sc2da_SetengData,                0, 0, 0, 0, 0, "SETENGDATA"},
    { sc2da_GetarrayName,             0, 0, 0, 0, 0, "GETARRAYNAME"},
    { sc2da_Clear,                     0, 0, 0, 0, 0, "CLEAR_STATE"},
    { sc2da_Pixelmask,                 0, 0, 0, 0, 0, "PIXEL_MASK"},
    { sc2da_Heatermask,                 0, 0, 0, 0, 0, "HEATER_MASK"},
    { sc2da_Heatersloperead,           0, 0, 0, 0, 0, "HEATERSLOPE_READ"},
    { sc2da_SQ1optpts,                 0, 0, 0, 0, 0, "SQ1OPT_POINTS"},
    { mceTest,                        0, 0, 0, 0, 0, "TEST"},
    { mceTryChkSum,                   0, 0, 0, 0, 0, "MCECHKSUM"}
   };

  //The Rows and Cols can be optionally specified on the command line. This 
  //allows us to run  this task with different frame size in simulation using
  // SDSU Timing board.

  // now, they are always set to full array size as default
  Cols = 32;
  Rows = 41;
  dasInfo.debuglvl=0;
  errNo--;
   
  ImpNodeName(NodeName,64, status);
  printf ( "NodeName = %s\n", NodeName );
  sprintf(domainName,"%s",getenv("DOMAINNAME"));

  // Get SDSU context to global con
  if ( (*status = sdsu_create(&con)) < 0) 
  {
    printf ("sc2da: failed to attach the shared memory\n"); 
    exit ( errNo);
  }
  errNo--;
  
  // Initialise global flags, after the shareDmem is established
  sc2dalib_variablesInit(con,&dasInfo,Cols,Rows,pixelMask, status);

  // Allocate buffers for storing the acquisition data
  *status = (StatusType) sdsu_allocate_buffers(con,dasInfo.bufSize);
  if (*status < 0) 
  {
    printf("sc2da: sdsu_allocate_buffers failed\n"); 
    exit(errNo) ;
  }
  errNo--;
  if( (*status=open_mastermsgfifo(con))<0 )
  {
    printf ("sc2da: failed to open dramamsgfifo\n"); 
    exit ( errNo);
  }
  errNo--;

  // Initialise DRAMA
  // use jit default if set to 0
  //jitSetDefaults( DITS_M_X_COMPATIBLE, 0.0, dramabufSize, 100000, 100000, 20000, &status );
  jitSetDefaults( DITS_M_X_COMPATIBLE, 0.0, 0, 0, 0, 0, status );
  jitAppInit ( taskName, status);
  if(*status != STATUS__OK)
  {
    printf("sc2da: Failed to initialize with jitsAppInt - exiting\n");
    exit(errNo);
  }
  //printf("sc2da: jitAppInit OK\n");
  errNo--;

  jitDebugSet( dasInfo.debuglvl, status );

  //  Initialise parameter system and create parameters  and msg queue
  sc2dalib_createSDP(&dasInfo, status);
  if( *status !=STATUS__OK )
  {
    printf ("sc2da: Failed to create parameters/queue \n"); 
    exit (errNo);
  }  
  errNo--;

  used4Dits.con=con;
  used4Dits.myInfo=&dasInfo;
  used4Dits.intermsg=glbMsg;

  // provide two call-back routines for trigging action reshedule
  DitsAltInClear ( &ditsAltDesc, status );

  DitsAltInAdd ( &ditsAltDesc, con->fifo[USER_USER_MESG].fdr, DITS_M_READ_MASK,
                 (DitsInputCallbackRoutineType) sc2dalib_callbkMsg, &used4Dits, status);
  if( *status !=STATUS__OK )
  {
    printf("sc2da: Failed to add [USER_USER_MESG] and sc2damsg\n"); 
    exit (errNo);
  }  
  errNo--;
  DitsAltInAdd ( &ditsAltDesc, con->fifo[DRV_USER_MESG].fdr, DITS_M_READ_MASK,
                 (DitsInputCallbackRoutineType) sc2dalib_callbkdrvMsg, &used4Dits, status);
  if( *status !=STATUS__OK )
  {
    printf ("sc2da: Failed to add [DRV_USER_MESG], drivermsg\n"); 
    exit (errNo);
  }  
  errNo--;

  sc2headman_makeaction (status);
  if(*status != STATUS__OK)
  {
    printf("sc2da: Failed after call to sc2headman_makeaction -exiting\n");
    exit(errNo);
  }
  errNo--;
 
  DitsPutActions(DitsNumber(Actions), Actions, status);
  if(*status != STATUS__OK)
  {
    printf("sc2da: DitsPutActions returns BAD status\n");
    exit(errNo);
  }
  errNo--;

  DitsPutOrphanHandler(sc2dalib_callorphanHandler, status);
  if( *status!=STATUS__OK)
  {
    printf("sc2da: DitsPutOrphanHandler returns BAD status\n");
    exit(errNo);
  }
  errNo--;

    
 // printf("sc2da: DitsAltInloop starts\n");
 // Enter DRAMA message handling loop, including messages from USER_USER_MESG
  DitsAltInLoop ( &ditsAltDesc, status );
  if( *status!=STATUS__OK)
  {
    ErsRep(0,status,"sc2da: DitsAltInloop returns BAD status");
    ErsFlush ( status );
    perror(0);
    exit(errNo);
  }
  //Task finishes 
  return(jitStop ( taskName, status ) );
}




/* ================== DRAMA    actions ==================== */
/* ==================                  ==================== */
/* ==================                  ==================== */
/* ================== DRAMA    actions ==================== */

/**
 * \fn void sc2da_Abort(StatusType *status)
 *
 * \brief drama action
 *  aborts the sequence.
 *
 * \param *status StatusType.  given and return
 *
 * it sends ST command to the MCE and the MCE will send back the last frame
 * with stop and last bits set
 *
 * ditscmd SC2DA ABORT
 */

/*+ sc2da_Abort - 
*/
void sc2da_Abort
(
StatusType *status
)
{
  int          rval;
  char         dateTime[40];
  long               in_sequence;

  if (*status != STATUS__OK) return;
  errno=0;   
  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  // update debug, check if it is in DA_MCE_SEQ setup actionFlag
  //sc2dalib_abortInit(con,&dasInfo,dateTime,status);
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_SEQ)
  {
    *status=DITS__APP_ERROR;	
    ErsRep(0,status, "sc2da_Abort: DA is not in SEQUENCE application!!" );
    return;
  }
  dasInfo.actionFlag=ABORTACTION;  
  fprintf(dasInfo.fpLog,"\n<%s> CMD for sc2da_Abort\n",dateTime);
  //End abortInit
  sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);  
  if (*status != STATUS__OK)  
  {
    ErsRep(0,status, "sc2da_Abort: sc2dalib_stopFramefailed");
  }
  else
   jitDebug(4,"sc2da_Abort: the action has completed\n");     
}


/**
 * \fn void sc2da_Batch(StatusType *status)
 *
 * \brief drama action
 *  execute a series of commands from a bacthfile
 *  allows users to setup command sequence to complete a test  
 *
 * \param *status StatusType.  given and return
 *
 * it sends several commands including single GOs in the batchFile, the data
 * frame from the MCE is saved in rtsltFile if filesvFlag=1. This allows users
 * to do some simple tests manually by constructing their own commands in the
 * batch file, such as setup single pixel, setup SSA or
 * some verification.
 *
 * ditscmd SC2DA  MCEBATCHGO \ 
 *         BATCH_FILE=xxx  \
 *         DATA_FILE=xxx \
 *         CMDREP_FILE=xxx \
 *         SVFILE_FLAG=1
 */
/*+ sc2da_Batch -
*/
void sc2da_Batch
(
StatusType *status
)
{
  int                     rval;
  DitsDeltaTimeType       timeout;
  char                    *glbmsgPtr,*lclmsgPtr; 
  static DRAMA_INNERMSG   dramamsg;
  static char             dateTime[40];
  static dasCmdInfo_t     batchCmd;
  static char             *byte;
  static int              mcedataWords;
  static long             batchDelay;

  if (*status != STATUS__OK) return;

  //check if it is start or completion.
  if(DitsGetSeq()==0)
  {
    rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
    // update debug, check if it is in DA_MCE_NONE setup actionFlag
    // read in args from cmd, set msgwrtPt, msgreadPT=0; 
    sc2dalib_batchInit(con,&dasInfo,dateTime,status);
    if (*status != STATUS__OK)
    {
      ErsRep(0,status,"sc2da_Batch: sc2dalib_batchInit failed");
      sc2dalib_actionfileEnd(con,&dasInfo,1,status);
      return;
    } 
    // get batchDelay (second) for delay between each GO
    SdpGeti("BATCH_DELAY", &batchDelay, status);
    byte=(char*)&batchCmd.reply;
    // send out all other cmds but not GO 
    sc2dalib_batchParser(con,&dasInfo,&batchCmd,dateTime,&glbMceInx,status);
    if(dasInfo.batchFlag == BATCH_GO)
    {
      // we have a first GO
      jitDebug(32,"sc2da_Batch: the %ldth: first GO =%s\n", 
               dasInfo.trkNo,batchCmd.mceCmd);
      usleep(batchDelay);
      
      sc2dalib_Cmd(con,&dasInfo,&batchCmd,dateTime,status);
      if (*status != STATUS__OK)
      {
        ErsRep(0,status,"sc2da_Batch: sc2dalib_Cmd failed");
        sc2dalib_actionfileEnd(con,&dasInfo,0,status);
        return;
      }
      // put timeout if read_data has not time out
      DitsDeltaTime(1, 0, &timeout);
      DitsPutDelay(&timeout, status);     
      dasInfo.actIndex=DitsGetActIndex();
      DitsPutRequest (DITS_REQ_SLEEP, status );
      post_sem(con->intask_sem,0);
    }
    else
    { 
      if (*status != STATUS__OK)
      {
        ErsRep(0,status,"sc2da_Batch: sc2dalib_BatchPareser failed");
      }
      // in case BATCH_END, BATCH_TRANSFAIL and BATCH_CMDFAIL:
      if(dasInfo.batchFlag == BATCH_END)
      {
        jitDebug(32,"sc2da_Batch: Total(%ld)commands,but no GO.\n",
                 dasInfo.trkNo);
        jitDebug(4,"sc2da_Batch: The action has completed\n");

	/* Just if we were building an AST map of the commands responses
           tell the process that we are done */
	if (dasInfo.filesvFlag == 5)
	  {
	    dasInfo.astMapState=2;
	    sc2dalib_Cmd(con,&dasInfo,&batchCmd,dateTime,status);
	  }
      }
      sc2dalib_actionfileEnd(con,&dasInfo,0,status);
    }
  }
  else
  {
    con->process.whereabout=Dits_GetSeqn;
    // get the glbmsg, we can use it diretly later
    glbmsgPtr=(char*)&glbMsg[dasInfo.msgreadPt];
    lclmsgPtr=(char*)&dramamsg;
    memcpy(lclmsgPtr,glbmsgPtr,sizeof(DRAMA_INNERMSG));       
    // recycle=16
    dasInfo.msgreadPt++;
    dasInfo.msgreadPt=dasInfo.msgreadPt & 0x000F ; 

    byte=(char*)dramamsg.data;
    mcedataWords=(int)(dramamsg.bufsize/4);
    jitDebug(32,"sc2da_Batch:framePacketSize =<%d>\n",mcedataWords);

    if (dasInfo.filesvFlag != 0)
    {
      sc2dalib_saveframeData(con,&dasInfo,dasInfo.fpData,byte,mcedataWords,status);
    }

    if(dramamsg.reason==FRAME_COMPLETION)
    {
      //check if we have more commands
      sc2dalib_batchParser(con,&dasInfo,&batchCmd,dateTime,&glbMceInx,status);
      if(dasInfo.batchFlag== BATCH_GO)
      {
        // we have another GO, 
        // driver code copy the datacount again to a local copy
        jitDebug(32,"sc2da_Batch: the %ldth: GO =%s\n",
                   dasInfo.trkNo,batchCmd.mceCmd);
        usleep(batchDelay);

        sc2dalib_Cmd(con,&dasInfo,&batchCmd,dateTime,status);
        if (*status != STATUS__OK)
        {
          ErsRep(0,status,"sc2da_Batch: sc2dalib_Cmd failed");
          sc2dalib_actionfileEnd(con,&dasInfo,0,status);
          return;
        }  
        dasInfo.actIndex=DitsGetActIndex();
        DitsDeltaTime(6, 0, &timeout);
        DitsPutDelay(&timeout, status);     
        DitsPutRequest (DITS_REQ_SLEEP, status );
        post_sem(con->intask_sem,0);
      }
      else
      { 
        if (*status != STATUS__OK)
        {
          ErsRep(0,status,"sc2da_Batch: mceBatchPareser failed");
        }
        // in case BATCH_END, BATCH_TRANSFAIL and BATCH_CMDFAIL:
        else
        {
          if(dasInfo.batchFlag== BATCH_END)
          {
            jitDebug(4,
              "sc2da_Batch: Total(%ld)commands. The action has completed\n",
               dasInfo.trkNo);
             con->process.seqstatus=SEQ_FINISHED;
          }
          else
            con->process.seqstatus=SEQ_ERROR;
        }
        sc2dalib_actionfileEnd(con,&dasInfo,0,status);
      }
    }
    else
    {
      *status=DITS__APP_ERROR;
      {
        sc2dalib_msgprintSave(&dasInfo,
           "sc2da_Batch: ERROR:%s",dramamsg.errRep,USE_ERSREP,status);
      }
      con->process.seqstatus=SEQ_ERROR;
      sc2dalib_actionfileEnd(con,&dasInfo,0,status);
    }
  }
}


/**
 * \fn void sc2da_Clear(StatusType *status)
 *
 * \brief drama action
 *  clear state flag so that other command can be executed.
 *
 * \param *status StatusType.  given and return
 *
 * ditscmd SC2DA CLEAR_STATE
 */

/*+ sc2da_Clear -
*/
void sc2da_Clear
(
StatusType *status
)
{
  if (*status != STATUS__OK) return;
  
  SdpPuti("IN_SEQUENCE", DA_MCE_NONE, status);
  con->process.framesetup =DA_MCE_NONE;
}



/**
 * \fn void sc2da_Config(StatusType *status)
 *
 * \brief drama action
 *  setup some configurations
 *
 * \param *status StatusType.  given and return
 *
 *
 * The CONFIGURE action checks the flags to ensure that the SCUBA2 DAS is 
 * initialised and that it is not currently executing a sequence. If the 
 * checks are satisfactory, it parses the configuration file and stores 
 * the information, the configured flag is set. 
 *
 * Example of configure file:
 *   Dont know yet
 *
 * ditscmd SC2DA CONFIGURE \
 *         OBS_MODE=DREAM \
 *         FRAME_NO2PROC=100 \
 *         CONFIGURATION=mce.xml
 */

/*+ sc2da_Config - 
*/
void sc2da_Config
(
StatusType *status
)
{   	        
  int          rval;
  char  dateTime[40];
  char  obsMode[FILE_LEN],configFile[FILE_LEN];
  char  shellcmd[]="checkForSc2dadh";
  SdsIdType argId;
  char  filename[FILE_LEN];
  char  myTaskName[10];
  double timeout;

  if (*status != STATUS__OK) return;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> start of sc2da_Config XXXX",dateTime);

  //If we were previously initialised, set all other flags to zero 

  SdpPuti("CONFIGURED",0,status);
  SdpPuti("SETUP",0,status);
  SdpPuti("IN_SEQUENCE",DA_MCE_NONE,status);

  // Check for the Data handling task being alive, if not throw an error

 if ( system ( shellcmd ) != 0 )
   {
     *status = DITS__APP_ERROR;
     ErsRep(0,status,"!!!!! sc2dadh is dead you must do an UNLOAD_INST - LOAD_INST cycle ");
     return;
   }


  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> sc2da_Config after checking for dhtask XXXX",dateTime);

  // update debug, check if it is in DA_MCE_NONE setup actionFlag
  // read in args from cmd 
  sc2dalib_configInit(con,&dasInfo,&glbMceInx,dateTime,obsMode, configFile,status);
  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_Config: sc2dalib_configInit failed");
    sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
    return;
  }
  jitDebug(2,"sc2da_Config: CONFIGURATION<%s>\n", configFile);
  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> sc2da_Config after sc2dalib_configInit XXXX",dateTime);


  // Only call the sc2headman routines if we are in RTS client mode
  if (dasInfo.engFlag==RTSC_MODE)
  {
    jitDebug(2,"sc2da_Config: call sc2headman_config\n ");
    sc2headman_config(dasInfo.obsNo,status);
    if (*status != STATUS__OK)
      {
	ErsRep(0,status,"sc2da_Config: sc2headman_config failed");
	return;
      }

    rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
    fprintf(dasInfo.fpLog,"\n<%s> sc2da_Config after sc2headman_config XXXX",dateTime);


    sc2headman_getmode(FILE_LEN,obsMode,status);
    jitDebug(2," RTSC: DATA_MODE=%s",obsMode); 

    // Read the current state of the MCE and write it to the 
    // files SC2SCRATCH/setupLog and SC2SCRATCH/astKeyMapFile
    ArgNew(&argId, status);
    ArgPuti(argId, "SVFILE_FLAG", 5, status); // Write the file
    ArgPuti(argId, "BATCH_DELAY", 0, status); // No delay
    sprintf(filename,"%s/logSetup.txt", getenv( "CONFIG_ALL" ));
    ArgPutString(argId, "BATCH_FILE", filename, status);
    sprintf(filename,"%s/setupLog", getenv( "SC2SCRATCH" ));
    ArgPutString(argId, "CMDREP_FILE", filename, status);
    sprintf(filename,"%s/setpLogData", getenv( "DATADIR" ));
    ArgPutString(argId, "DATA_FILE", filename, status);
    sprintf(myTaskName,"%s",getenv( "SC2DRAMATASK" ));
    timeout = 10.0;
    jitObeyWait(myTaskName, "MCEBATCHGO", argId, 0, timeout, status);
    ArgDelete(argId, status);

    if(*status != STATUS__OK)
      {
	ErsRep(0,status,"sc2da_Config: Failed when reading MCE state. status: 0x%x", *status);
	return;
      }

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> sc2da_Config after reading the MCE configuration XXXX",dateTime);


  }
  else
    MsgOut(status," ATC ENGMODE: DATA_MODE=%s, don't call sc2headman_xx",obsMode); 

  jitDebug(2,"sc2da_Config: OBS_MODE <%s> FRAME_NO2PROC<%d>\n",obsMode,
       dasInfo.parshmPtr->procNo); 

  // pass parameter got from CONFIGURE to sharedMem
  sc2dalib_configInitSet(&dasInfo,dateTime,obsMode, configFile,status);
  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> sc2da_Config after sc2dalib_configInitSet XXXX",dateTime);
 

  jitDebug(2,"sc2da_Congi: domainname=%s\n",domainName);
  if ( dasInfo.engFlag==RTSC_MODE  )
  {
    if ( strcmp(domainName,"ROE.LINUX") ==0 )
      MsgOut(status,"sc2da_Config: It is not at JAC, don't call sc2headman_startenviro");
    else
    {
      // start ENV monitor
      sc2headman_startenviro(taskName,status);
      if (*status != STATUS__OK)
        ErsRep(0,status,"sc2da_Config: sc2headman_startenviro failed");
    }
  }

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> sc2da_Config after sc2headman_startenviro  XXXX",dateTime);

  MsgOut(status,"Ending CONFIGURE action");
  jitDebug(4,"sc2da_Config: the action has completed\n");
  sc2dalib_actionfileEnd(con,&dasInfo,0,status);  
}


/**
 * \fn void sc2da_ConfigKick(StatusType *status)
 *
 * \brief drama action:
 *  kick CONFIGURATION. 
 *
 * \param *status StatusType.  given and return
 *
 *
 * ditscmd SC2DA -k CONFIGURATION 
 *
 */
/*+ sc2da_ConfigKick - Handles a "kick" of the configure action 
*/
void sc2da_ConfigKick
( 
StatusType *status    /* global status (given and returned) */
)
{
}



/**
 * \fn void sc2da_Downld2PCI(StatusType *status)
 *
 * \brief drama action
 *  downloads DSP *.lod file into PCI DSP
 *
 * \param *status StatusType.  given and return
 *
 *
 * ditscmd SC2DA DWLOADDSP filename
 */

/*+ sc2da_Downld2PCI
*/
void sc2da_Downld2PCI
(
StatusType *status
)
{
  int  rval;
  char dspfile[FILE_LEN];
  char dateTime[40];
  
   
  if (*status != STATUS__OK) return;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);

  // update debug, check if it is in DA_MCE_NONE setup actionFlag
  // read in args from cmd 
  sc2dalib_downld2pciInit(con,&dasInfo,dspfile,dateTime,status);
  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_Downld2PCI: sc2dadwonld2pciInit failed");
    sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
    return;
  }  
  sc2dalib_download2PCI(con,dspfile,status);
  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_Downld2PCI: failed after call to sc2dadownload2PCI"); 
  }
  else
    jitDebug(4, "sc2da_Downld2PCI: the action has completed" );
  sc2dalib_actionfileEnd(con,&dasInfo,0,status);
}



/**
 * \fn void sc2da_Down2FPGA(StatusType *status)
 *
 * \brief drama action
 *  downloads *.jbc file into MCE FPGA
 *
 * \param *status StatusType.  given and return
 *
 *
 * It reads xxx.jbc file into memory buffer, divides the buffer into 
 * multi-wbSize (=232 bytes) + left over, copies wbSize or the left-over 
 * bytes into WriteBlock command structure, insets checksum, sends each 
 * command to the MCE and saves the data sent out into filename if 
 * filesvFlag=1.
 *
 *ditscmd SC2DA DOWNLOADFPGA \
 *        cardId sram1(or 2) \
 *        xxx.jbc \
 *        filesvFlag \
 *        filename
 */

/*+ sc2da_Down2FPGA -
*/
void sc2da_Down2FPGA
(                    
StatusType *status
)
{
   if (*status != STATUS__OK) return;

//  temporarily take out to SSAOUT
}  



/**
 * \fn void sc2da_Dispinfo(StatusType *status)
 *
 * \brief drama action 
 *  disply some internal information
 *
 * \param *status StatusType.  given and return
 *
 * >ditscmd SC2DA   DISPINFO
 *
 */
/*+ sc2da_Dispinfo - 
*/
void sc2da_Dispinfo
(
StatusType *status    
)
{
  long    in_sequence;

  if (*status != STATUS__OK) return;

  SdpGeti("IN_SEQUENCE",&in_sequence,status);
  sc2dalib_dispInfo(con, &dasInfo,in_sequence,status);
 
  if (in_sequence != DA_MCE_NONE)
    return;
  else
    sc2dalib_pcistatusRep(con,&dasInfo,status);
}



/**
 * \fn void sc2da_EndObs(StatusType *status)
 *
 * \brief drama action
 *  end observation, tidy-up
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * ditscmd SC2DA END_OBSERVATION
 */

/*+ sc2da_EndObs -
*/
void sc2da_EndObs 
( 
StatusType *status    /* global status (given and returned) */
)
{
  long  setup;

  if (*status != STATUS__OK) return;

  SdpGeti("SETUP", &setup, status);
  if(setup != 0)
  {

    if (dasInfo.engFlag==RTSC_MODE)
       sc2headman_endobs(status);
    // MsgOut ( status, " END_OBSERVATION called");
  }
}


/**
 * \fn void sc2da_ObsKick(StatusType *status)
 *
 * \brief drama action
 *  kick end observation action
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * ditscmd SC2DA END_OBSERVATION
 */

/*+ sc2da_ObsKick -
*/
void sc2da_ObsKick 
( 
StatusType *status    /* global status (given and returned) */
)
{
  if (*status != STATUS__OK) return;

  MsgOut ( status, " Kick END_OBSERVATION called");
}



/**
 * \fn void sc2da_Exit(StatusType *status)
 *
 * \brief drama action
 *  free buffers and SDSU_CONTEXT structure before
 *  end the  drama task
 *
 * \param *status StatusType.  given and return
 *
 *
 * It tries to tell the data handle task to exit by setting a flag and using 
 * the semaphore. If the sc2da_Seq is in the middle of executing a sequence, 
 * sc2da_Exit() will use sc2da_Abort to stop the sequence. 
 *
 * In either case, sc2da_Exit() then tells DRAMA that the application is to 
 * exit.
 *
 * >ditscmd SC2DA EXIT
 */
/*+ sc2da_Exit 
*/
void sc2da_Exit 
(
StatusType *status 
)
{
  char              dateTime[40];
  int               rval;
  long              init=0,in_sequence;
  DitsDeltaTimeType timeout;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  SdpGeti("INITIALISED",&init,status);
  if ( init ) 
     fprintf(dasInfo.fpLog,"\n<%s> CMD from sc2da_Exit\n",
          dateTime);
	  
  SdpGeti("IN_SEQUENCE",&in_sequence,status);
  
  if (DitsGetSeq()==0 && in_sequence==DA_MCE_SEQ)
  {        
    if (*status != STATUS__OK) 
    {
      ErsRep(0,status,"sc2da_Exit: Error- status!=OK ");
      return;
    }
    // if in seq, abort and exit
    if ( in_sequence==DA_MCE_SEQ)  
    {
      MsgOut(status,"sc2da_Exit: Sequence is INACTION, call sc2da_Abort.");
      MsgOut(status,"         reschedule and wait for exiting ....");	
	sc2da_Abort(status);
      DitsDeltaTime(3, 0, &timeout);
      // Reschedule to wait for  timer out.
      DitsPutDelay(&timeout, status);
      DitsPutRequest(DITS_REQ_SLEEP, status);
    }
  }
  else
  {
    jitDebug(2,"sc2da_Exit: not in DA_MCE_SEQ\n");
    if ( in_sequence !=DA_MCE_NONE)
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, 
        "Exit: %s has not completed",seqStatus[dasInfo.actionFlag]);
      return;
    }
    MsgOut(status," Exiting now ..........");
    dasInfo.actionFlag=EXITACTION;
    
    SdpPuti("IN_SEQUENCE",DA_MCE_NONE,status);
    con->process.exit=1;
    if ( con->fifo[USER_USER_MESG].fdw !=-1)
    {
      // release data handle task to exit, ie close_buffer,
      //: unmapping and close device-file
      post_sem(con->intask_sem,0); 
      usleep(2500);
      //sdsu_free wait_sem until data handle task exit. but what if data
      // handle stucks? it shall not
      *status =(StatusType )sdsu_free_buffers(con);
      if( *status <0 ) 
      {
        ErsRep(0,status,"sc2da_Exit: sdsu_free_buffers failed");
        *status=STATUS__OK;
      }
    }
    else //as sdsu_free is waiting for the sem
      post_sem(con->sync, 0);

    if( sdsu_free(con)!=SDSU_OK)	 
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,"sc2da_Exit: Error- sdsu_free failed");
    }

    if (dasInfo.fpLog!=NULL)
       fclose(dasInfo.fpLog);
    dasInfo.fpLog = NULL;

    jitDebug(2,"sc2da_Exit: call sc2dalib_closesharedMem\n");
    // move to here in case no data cmd is issued before
    sc2dalib_closesharedMem(&dasInfo,SHAREDM_PAR,status);
    if (*status != STATUS__OK) 
    {
      ErsRep(0,status,"sc2da_Exit: sc2dalib_closesharedMem failed");
    }

    jitDebug(2,"sc2da_Exit: call DitsPutRequest ( DITS_REQ_EXIT) \n");
    // Tell DRAMA this task is to exit
    //  DITS_REQ_END, drama task still there
    // DITS_REQ_EXIT drama task exits 
    DitsPutRequest ( DITS_REQ_EXIT, status );
    jitDebug(2,"sc2da_Exit: called DitsPutRequest ( DITS_REQ_EXIT) \n");
  }
}


/**
 * \fn void sc2da_GetarrayName(StatusType *status)
 *
 * \brief drama action
 *  get MCE port number and find the name from the table 
 *
 * \param *status StatusType.  given and return
 *
 * ditscmd SC2DAx GETARRAYNAME 
 * 
 */
/*+ sc2da_Getarrayname -
*/
void sc2da_GetarrayName
(
StatusType *status
)
{
  int  mcePort,i;
  char arrayID[]="rb cc array_id 1 ";

  if (*status != STATUS__OK) return;
 
  if (dasInfo.doneReadxml )
  {
    sc2dalib_readmceVal(con,&dasInfo,&glbMceInx,arrayID, &mcePort,1,status);
    if (*status != STATUS__OK)
      return;
    for (i=0; i<MAX_SUBARRAY; i++)
    {
      if(dasInfo.subArray[i].mceport==mcePort)
        break;
    }
    if( i==MAX_SUBARRAY )
    {
       *status=DITS__APP_ERROR;
      ErsRep(0,status, "sc2da_getarrayname: mceport is wrong");
      return;
    }
    MsgOut(status," MCE PortNo=%d, arrayName=%s",mcePort,dasInfo.subArray[i].id);
  }
  else
    MsgOut(status," DA is not initialised");
}


/**
 * \fn void sc2da_Heatslope(StatusType *status)
 *
 * \brief drama action
 *  find heater slope for heater track
 * 
 *
 *  setdebug 16 to display some message related to function calls
 *
 * \param *status StatusType.  given and return
 *
 *
 * ditscmd SC2DA HEAT_SLOPE \
 *               DATA_FILE=$CURRENTDATADIR/$datafile \
 *               SETUP_FILE=<correct directory>/$setupfile 
 */

/*+ sc2da_Heatslope
*/
void sc2da_Heatslope
( 
StatusType *status    /* global status (given and returned) */
)
{ 
  int         rval;
  static char dateTime[40];

  if (*status != STATUS__OK) return;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  jitDebug(16,"sc2da_Heatslope: call sc2dalib_heaterslopeInit\n");     
  sc2dalib_heaterslopeInit(con, &dasInfo,dateTime, status);
  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_Heatslope: sc2daheaterslopeInit failed"); 
    sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
    return;
  }
  sc2dalib_heaterSlope(&dasInfo,heaterMask,heaterSlope,status);
  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_Heatslope: sc2daheaterSlope failed"); 
  }
  sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
}


/**
 * \fn void sc2da_Heatermask(StatusType *status)
 *
 * \brief drama action
 *  read heater Mask to heaterMask use it in heaterslope
 *
 * \param *status StatusType.  given and return
 *
 * ditscmd SC2DAx PIXEL_MASK
 * 
 */
/*+ sc2da_Heatermask
*/
void sc2da_Heatermask
(
StatusType *status
)
{
  int          i, j, pixel;
  SdsIdType    argId=0;
  SdsIdType    maskId=0;
  char         name[FILE_LEN];

  if (*status != STATUS__OK) return;

  // the pixealmask.xml need to be in $CONFIG_HARD
  argId = DitsGetArgument();

  jitArgGetXML ( argId, "HEATER_MASK", 1, NULL, sizeof(name), name,
                 &maskId, status );
  if (*status != STATUS__OK)
  {
   ErsRep(0,status, "sc2da_Heatermask: jitArgGetXML failed");
   return;
  }
  MsgOut(status," pixelMaskXML=%s",name);

  sc2dalib_readpixelMask(&dasInfo,heaterMask,maskId,status);
  SdsFreeId ( maskId, status );
  SdsFreeId ( argId, status );

  if (*status != STATUS__OK)    return;

  if (dasInfo.fpLog != NULL)
  {
    fprintf(dasInfo.fpLog,"===== HEATER_PIXEL_MASK ARRAY ====\n");
    for ( i=0; i<40; i++ )
    {
      for ( j=0; j<32; j++ )
      {
        pixel= i*32 +j;
        fprintf(dasInfo.fpLog,"%3d", heaterMask[pixel]);
      }    
      fprintf(dasInfo.fpLog,"\n");
    }
    fprintf(dasInfo.fpLog,"===== HEATER_PIXEL_MASK ARRAY ====\n");
    fflush(dasInfo.fpLog);
    fflush(dasInfo.fpLog);
  }
  MsgOut(status," sc2da_Heatermask OK");
}


/**
 * \fn void sc2da_Heatersloperead(StatusType *status)
 *
 * \brief daram action
 *  read heater slope and store in heaterSlope for tracking
 *
 * \param *status StatusType.  given and return
 *
 */
/*+ sc2da_Heatersloperead 
*/
void sc2da_Heatersloperead
( 
StatusType *status    
)
{
  int   i,j,pixel;


  if (*status != STATUS__OK) return;

  sc2dalib_heaterslopeRead(&dasInfo,heaterSlope,status);
  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_Heatersloperead: sc2dalib_heaterslopeRead failed");
    return;
  }  
  if (dasInfo.fpLog != NULL)
  {
    fprintf(dasInfo.fpLog,"===== HEATER_SLOPE ARRAY ====\n");
    for ( i=0; i<40; i++ )
    {
      for ( j=0; j<32; j++ )
      {
        pixel= i*32 +j;
        fprintf(dasInfo.fpLog,"%11.3f", heaterSlope[pixel]);
      }    
      fprintf(dasInfo.fpLog,"\n");
    }
    fprintf(dasInfo.fpLog,"===== HEATER_SLOPE ARRAY ====\n");
    fflush(dasInfo.fpLog);
    fflush(dasInfo.fpLog);
  }
  MsgOut(status," sc2da_Heatersloperead OK");
}


/**
 * \fn void sc2da_Init(StatusType *status)
 *
 * \brief drama action
 *  initialise the task, setup MCE comamnd/parameter 
 *  lookup table by reading in mce.xml  
 *
 * \param *status StatusType.  given and return
 *
 * RTSC_MODE
 * ditscmd SC2DAx INITIALISE \
 *         INITIALISE=/jac_sw/itsroot/install/scuba2Da/data/initialise.xml \
 *         SIMULATE=0 
 * 
 * ENG_MODE
 * ditscmd SC2DAx INITIALISE \
 *         LOG_FILE=$LOGFILEDIR/logcmd \
 *         MCE_XML_FILE=mce.xml \
 *         HARDWARE=mce \
 *         DATA_FORMAT=BINARY 
 */
/*+ sc2da_Init -
*/
void sc2da_Init
(
StatusType *status
)
{  
  int        rval, mcePort, i;
  char       xmlfile[FILE_LEN];
  char       arrayName[40];
  SdsIdType  id;


  if (*status != STATUS__OK) return;

  rval=sc2dalib_finddateTime(DAS_DATE,dasInfo.Date); 

  //If we were previously initialised, set all other flags to zero 
  SdpPuti("CONFIGURED",0,status);
  SdpPuti("SETUP",0,status);
  SdpPuti("IN_SEQUENCE",DA_MCE_NONE,status);

  // Clear the heater tracking has failed flag
  dasInfo.heatTrackFailed = 0;
 
  // make up wavelen, filter now 
  sc2dalib_initInit(con,&dasInfo,xmlfile,status);
  fflush(dasInfo.fpLog);

  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_Init: failed after call to sc2dalib_initInit");
    sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
    return;
  }   
  errno=0;
  mcexml_readXML (xmlfile,(int *)status);
  if ( *status ==DITS__APP_ERROR)
  {    
    ErsRep(0,status,"sc2da_Init:: failed after call to mcexml_readXML"); 
    return;
  }
  mcexml_returnXML( &glbMceInx,(int*)status);
  if ( *status ==DITS__APP_ERROR)
  {
    ErsRep(0,status,"sc2da_Init: failed to mcexml_returnXML");
    return;
  }
  dasInfo.doneReadxml=1;

  // Check to see if we need to over ride the optical code reader
  mcePort = 0xFF;  // tells sc2dalib_initgetArrayID to read port from MCE

  for(i=0; i<MAX_SUBARRAY; i++)
    {
      if(strcmp(dasInfo.givenarrayName, dasInfo.subArray[i].id) == 0 )break;
    }
  if((i < MAX_SUBARRAY) && dasInfo.subArray[i].override) mcePort = dasInfo.subArray[i].mceport;

  // only use logfile to record errors
  sc2dalib_initgetArrayID(con,&dasInfo,&glbMceInx,arrayName,xmlfile,&mcePort,status);
  if ( *status ==DITS__APP_ERROR)
   return;

  /*  MsgOut (status,"sc2da_Init: read MCE portNo=%d  arrayName (initialise.xml) is %s",
      mcePort,arrayName); */

  /* Here is the test for the optical bar code reader working */
  if (strcmp (dasInfo.givenarrayName, arrayName) != 0)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"sc2da_Init: arrayName defined in initilaise.xml != giveName(%s)",
           dasInfo.givenarrayName);
    return;
  }
  //depends on INITFLAG in initialise.xml to Reset the MCE and the PCI
  if( (dasInfo.initFlag & SC2RTSC__MCERESET) !=0 )
  {
    sc2dalib_pcisendCmd(con,"RESETMCE",&dasInfo,0,0,"y",status);
    sc2dalib_pcisendCmd(con,"RESET",&dasInfo,0,0,"y",status);
    if( *status != STATUS__OK )
    {
      ErsRep(0,status,"sc2da_Init: failed after call to sc2dalib_pcisendCmd");
      return;
    }
    usleep(250);  
    // reset pci frame count=0
    sc2dalib_pcisendCmd(con,"WRITE",&dasInfo,1,0,"X",status);
    if( *status != STATUS__OK )
    {
      ErsRep(0,status,"sc2da_Init: failed after call to sc2dalib_pcisendCmd X(1)=0");
      return;
    }
    MsgOut(status," RESETMCE RESETPCI");
  }
  /*  else
      MsgOut(status," Don't RESETMCE RESETPCI"); */

  fprintf(dasInfo.fpLog,"\n===== subarray: %s ==",arrayName);
  fflush(dasInfo.fpLog);

  // create return argument
  ArgNew ( &id, status );
  ArgPuti ( id, "ROW_LEN", dasInfo.rowLength, status );
  ArgPuti ( id, "NUM_ROWS", dasInfo.numRows, status );
  ArgPutString ( id, "SUBNAME", arrayName, status );
  ArgPutString ( id, "TASKNAME", taskName, status );

  DitsPutArgument ( id, DITS_ARG_DELETE, status );  
  if ( *status ==DITS__APP_ERROR)
  {    
    ErsRep(0,status,"sc2da_Init:: failed after call to DitsPutArgument (DITS_ARG_DELETE)"); 
    return;
  }
 
  // now call the headman,
  if (dasInfo.engFlag==RTSC_MODE)
  {
    MsgOut(status,"sc2da_Init:RTSC MODE, call sc2headman_init");
    sc2headman_init(arrayName,dasInfo.chipId, dasInfo.filter,
                    dasInfo.wavelen,status);
    if ( *status ==DITS__APP_ERROR)
    {    
      ErsRep(0,status,"sc2da_Init:: failed after call to sc2headman_init"); 
      return;
    }
  }
  else
    MsgOut(status,"sc2da_Init:in ATC ENGMODE, don't call sc2headman_int");
 
  jitDebug(2,"subarray id =%s  chipid =%s filter=%s wavelen=%2.4e\n",
             arrayName, dasInfo.chipId, dasInfo.filter,dasInfo.wavelen);
  SdpPuti("INITIALISED",1, status);
}                                                                       



/**
 * \fn void sc2da_initKick(StatusType *status)
 *
 * \brief drama action
 *  kick INITIALISE. 
 *
 * \param *status StatusType.  given and return
 *
 *
 * ditscmd SC2DAx -k INITIALISE
 *
 */
/*+ sc2da_initKick - Handles a "kick" of the initialise action 
*/
void sc2da_initKick
( 
StatusType *status    /* global status (given and returned) */
)
{
}


/**
 * \fn void sc2da_Mcecmd(StatusType *status)
 *
 * \brief drama action
 *  send a command to the MCE. 
 *
 * \param *status StatusType.  given and return
 *
 *
 * it constructs 64 32-bit words for the command, insert checksum as last
 * word and send to the MCE. the reply from the MCE is saved in rsltFile.
 *
 * ditscmd SC2DAx  MCECMD \
 *         MCE_CMD="$cmd" \
 *         CMDREP_FILE=$CMDREPDIR/$resultfile \
 *         SVFILE_FLAG=$svfileFlag
 */

/*+ sc2da_Mcecmd - 
*/
void sc2da_Mcecmd
(
StatusType *status
)
{
  int           rval;
  dasCmdInfo_t  mymceCmd;
  char          *token,*dupcmd;
  char          dateTime[40];
   
  if (*status != STATUS__OK) return;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);

  // update debug, check if it is in DA_MCE_NONE setup actionFlag
  // read in args from cmd 
  jitDebug(1,"sc2da_Mcecmd: call sc2dalib_mcecmdInit\n");
  
  sc2dalib_mcecmdInit(con,&dasInfo,&mymceCmd,dateTime,status);
  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_Mcecmd: sc2dalib_mcecmdInit failed");
    sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
    return;
  }  
  
  //block GO if it is issued by user
  dupcmd=strdupa(mymceCmd.mceCmd);  
  token=strtok(dupcmd," ");
  if(strcmp("GO",token)==0 || strcmp("go",token)==0)
  {
    *status=DITS__APP_ERROR; 
    ErsRep (0, status,"sc2da_Mcecmd: GO shall not be used by this action"); 
    sc2dalib_actionfileEnd(con,&dasInfo,0,status);   
    return;
  }

  // check if card type is EXTERNAL
  // the cmd is not be sent in sendCmd 
 jitDebug(1,"sc2da_Mcecmd: call sc2dalib_sendCmd\n"); 
 sc2dalib_sendCmd(con,&dasInfo,&mymceCmd,&glbMceInx,dateTime,status);
  if (*status != STATUS__OK)   
  {
    ErsRep(0,status,"sc2da_Mcecmd: failed after call to sc2dalib_sendCmd");
    sc2dalib_actionfileEnd(con,&dasInfo,0,status);
    return;
  }

  //print msg for EXTTERNAL
  dupcmd=strdupa(mymceCmd.mceCmd);  
  token=strtok(dupcmd," ");
  token = strtok (NULL, " ");
  if (strcmp(token, EXTERNAL_CARD) ==0)
  {
    // the EXTERNAL one only for store setting in mceInxpt
     // don't send out
    MsgOut(status," Don't send \"%s\" to MCE",
           mymceCmd.mceCmd);
  }
  else
  {
    // all errors message have been displayed in mceErrRep
    sc2dalib_dispResults(&dasInfo,&mymceCmd,status);
    sc2dalib_chkChecksum(&dasInfo,&mymceCmd,status);
    if (*status != STATUS__OK)
    {
      ErsRep(0,status,"sc2da_Mcecmd: the action completed with error"); 
      sc2dalib_actionfileEnd(con,&dasInfo,0,status);
      return;
    }
  }
  jitDebug(4,"sc2da_Mcecmd: the action has completed\n");
  sc2dalib_actionfileEnd(con,&dasInfo,0,status);
}



/**
 * \fn void sc2da_Mceonflycmd(StatusType *status)
 *
 * \brief drama action
 *  poplulate onthflycmd strut and set flag. 
 *
 * \param *status StatusType.  given and return
 *
 *
 * it constructs 64 32-bit words for the command, insert checksum as last
 * word 
 *
 * ditscmd SC2DAx  MCEONFLYCMD \
 *         MCE_CMD="$cmd" 
 */

/*+ sc2da_Mceonflycmd - 
*/
void sc2da_Mceonflycmd
(
StatusType *status
)
{
  int           rval;
  char          *token,*dupcmd;
  char          dateTime[40];
   
  if (*status != STATUS__OK) return;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);

  // update debug, check if it is in DA_MCE_NONE setup actionFlag
  // read in args from cmd 
  sc2dalib_mceonflycmdInit(con,&dasInfo,&onflyCmd,dateTime,status);
  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_Mceonflycmd: sc2dalib_mceonflycmdInit failed");
    return;
  }  
  MsgOut(status, "onflyCmd=%s",onflyCmd.mceCmd);
    
  //block GO if it is issued by user
  dupcmd=strdupa(onflyCmd.mceCmd);  
  token=strtok(dupcmd," ");
  if(strcmp("GO",token)==0 || strcmp("go",token)==0)
  {
    *status=DITS__APP_ERROR; 
    ErsRep (0, status,"sc2da_Mceonflycmd: GO shall not be used by this action"); 
    return;
  }
   // set flag to tell Pixemon action there is a onthefly cmd
  dasInfo.cmdFlag=1;
  jitDebug(4,"sc2da_Mceonflycmd: the action has completed\n");
}


/**
 * \fn void sc2da_MceStatus(StatusType *status)
 *
 * \brief drama action
 *  get all MCE status from mcexml_struct structure.
 *
 * \param *status StatusType.  given and return
 *
 * ditscmd SC2DA MCESTATUS
 */

/*+ sc2da_MceStatus - 
*/
void sc2da_MceStatus
(
StatusType *status
)
{
  int         i, rval;
  static char   dateTime[40];
  static char   *notRBCMD[]=
{
"RB sys card_type",  "RB sys slot_id", "RB sys dip",    "RB sys sample_dly",  //4
"RB sys sample_num", "RB sys fb_dly",  "RB sys flx_lp_init",                  //7
"RB sys sdsutime_errtest", "RB sys sdsutime_parablk",                         //9

"RB psc card_type",  "RB psc slot_id",  "RB psc dip",    "RB psc sample_dly", //13
"RB psc sample_num", "RB psc fb_dly",   "RB psc flx_lp_init",                 //16
"RB psc sdsutime_errtest","RB psc sdsutime_parablk",                          //18
"RB psc pow_ctrl",   "RB psc card_temp","RB psc card_id","RB psc fw_rev",     //22
"RB psc led",                                                                 //23

"RB cc dip",    "RB cc sample_dly",  "RB cc sample_num", "RB cc fb_dly",      //27
"RB cc flx_lp_init", "RB cc upload_fw",                                       //29
"RB cc sdsutime_errtest","RB cc sdsutime_parablk",                            //31

"RB ac card_type",  "RB ac slot_id",    "RB ac dip",                          //34
"RB ac sdsutime_errtest", "RB ac sdsutime_parablk",                           //36

"RB rc1 col_map",   "RB rc1 card_type", "RB rc1 dip",                         //39
"RB rc1 sdsutime_errtest", "RB rc1 sdsutime_parablk",                         //41

"RB rc2 col_map",   "RB rc2 card_type", "RB rc2 dip",                         //44
"RB rc2 sdsutime_errtest", "RB rc2 sdsutime_parablk",                         //46

"RB rc3 col_map",   "RB rc3 card_type", "RB rc3 dip",                         //49
"RB rc3 sdsutime_errtest", "RB rc3 sdsutime_parablk",                         //51

"RB rc4 col_map",   "RB rc4 card_type", "RB rc4 dip",                         //54
"RB rc4 sdsutime_errtest", "RB rc4 sdsutime_parablk",                         //56

"RB bc1 sa_htr0",   "RB bc1 sa_htr1",   "RB bc1 card_type",    "RB bc1 dip",  //60
"RB bc1 sdsutime_errtest", "RB bc1 sdsutime_parablk",                         //62

"RB bc2 sa_htr0",   "RB bc2 sa_htr1",   "RB bc2 card_type",    "RB bc2 dip",  //66
"RB bc2 sdsutime_errtest", "RB bc2 sdsutime_parablk",                         //68

"RB bc3 sa_htr0",   "RB bc3 sa_htr1",   "RB bc3 card_type",    "RB bc3 dip",  //72
"RB bc3 sdsutime_errtest", "RB bc3 sdsutime_parablk",                         //74

"RB cc card_type",                                                            //75
};

  if (*status != STATUS__OK) return;

  // Check whether this is the start or completion  
  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);   
  // update debug, check if it is in DA_MCE_SEQ setup actionFlag
  sc2dalib_mcestatusInit(con,&dasInfo,dateTime,status);
  if (*status != STATUS__OK)
  {
    ErsRep(0,status,"sc2da_MceStatus: sc2dalib_mcestatusInit failed");
    sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
    return;
  }  
  jitDebug(4,"sc2da_MceStatus: sc2dalib_mcestatusInit OK\n");
  //get all setting values from the glbMceInx
  if (dasInfo.batchFlag !=1)  
  {
    mcexml_status( dasInfo.fpMcecmd,glbMceInx,(int*)status);
   jitDebug(4,"sc2da_MceStatus: mcexml_status OK\n");
  }
  else  // send RB commands
  {
    /********** try read   */
    dasCmdInfo_t  mymceCmd;
    int done = 0;

    mcexml_initenq (glbMceInx,(int*)status);
    for ( ; ; )
    {
      if ( *status != STATUS__OK )
      {
        break;
      }
      mcexml_getenq ( glbMceInx, mymceCmd.mceCmd, &done, (int*)status );
      if ( done == 1 )
      {
        break;
      }
      jitDebug(2,"cmd: %s\n", mymceCmd.mceCmd);
      rval=0;
      for (i=0; i<75; i++) 
	{
	  if ( strcmp(mymceCmd.mceCmd, notRBCMD[i]) == 0)
	    {
	      rval=1;
	      break;
	    } 
	}
      if( rval==0)
      {
        sc2dalib_sendCmd(con,&dasInfo,&mymceCmd,&glbMceInx,dateTime,status);
        *status=STATUS__OK;
/**********
        if ( !StatusOkP(status))   
        {
          fprintf(dasInfo.fpMcecmd,"\n</HEADER>\n\n");
          ErsRep(0,status,"sc2da_MceStatus: failed after call to sc2dalib_sendCmd");
          sc2dalib_actionfileEnd(con,&dasInfo,0,status);
          return;
        }
***********/
      }
    }
  }
  fprintf(dasInfo.fpMcecmd,"\n</HEADER>\n\n");
  sc2dalib_actionfileEnd(con,&dasInfo,0,status);
  jitDebug(4,"sc2da_MceStatus: the action has completed\n");
}



/**
 * \fn void sc2da_PCIcmd(StatusType *status)
 *
 * \brief drama action
 *  send a command to PCI DSP
 *
 * \param *status StatusType.  given and return
 *
 *
 * >ditscmd SC2DAx  PCICMD \
 *          PCI_CMD=WRITE \
 *          PCI_MEMTYPE=y \
 *          PCI_MEMADDR=xxx\
 *          PCI_MEMVALUE=xxx
 */

/*+  sc2da_PCIcmd- 
*/
void sc2da_PCIcmd
(
StatusType *status
)
{   	        
  long        memAddr;
  long        memValue;
  char        cmd[CMD_LEN],memType[10];
  char        dateTime[40];
  int         rval;
  
  if (*status != STATUS__OK) return;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  sc2dalib_pcicmdInit(con,&dasInfo,cmd,memType,&memAddr,&memValue,
                  dateTime,status);
  if( *status!=STATUS__OK)
  {
    ErsRep(0,status,"sc2da_PCIcmd: failed after call to sc2dalib_pcicmdInit");
    return;
  }
  sc2dalib_pcisendCmd(con,cmd,&dasInfo,memAddr,memValue,memType,status);
  if( *status!=STATUS__OK)
    ErsRep(0,status,"sc2da_PCIcmd: faile after call to sc2dalib_pcisendCmd");
  else
     jitDebug(4,"sc2da_PCIcmd: the action has completed\n");
}



/**
 * \fn void sc2da_PCIblk(StatusType *status)
 *
 * \brief drama action: 
 *  send read command to read block memory from  PCI DSP
 *
 * \param *status StatusType.  given and return
 *
 *
 * >ditscmd SC2DA PCIBLK \
 *          PCI_MEMTYPE=y \
 *          PCI_STARTADDR=xxx \
 *          PCI_BLKSIZE=64 
 * 
 */

/*+  sc2da_PCIblk- 
*/
void sc2da_PCIblk
(
StatusType *status
)
{   	        
  long        startAddr=100;
  long        blkSize=64;
  char        cmd[CMD_LEN],memType[10];
  char        dateTime[40];
  int         rval;
  
  if (*status != STATUS__OK) return;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  strcpy(cmd,"READ");
  sc2dalib_pciblkInit(con,&dasInfo,cmd,memType,&startAddr,&blkSize,
                  dateTime,status);
  if( *status!=STATUS__OK)
  {
    ErsRep(0,status,"sc2da_PCIblk: failed after call to sc2dalib_pciblkInit");
    return;
  }
  sc2dalib_pcisendblkCmd(con,cmd,&dasInfo,startAddr,blkSize,memType,status);
  if( *status!=STATUS__OK)
    ErsRep(0,status,"sc2da_PCIblk: faile after call to sc2dalib_pcisendblkCmd");
  else
     jitDebug(4,"sc2da_PCIblk: the action has completed\n");
}

 

/**
 * \fn void sc2da_Ping(StatusType *status)
 *
 * \brief drama action
 *  send a response message to the user interface to  confirm that the 
 *  task is stil alive.
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * ditscmd SC2DA PING
 */

/*+ sc2da_Ping
*/
void sc2da_Ping
( 
StatusType *status  
)
{
  if (*status != STATUS__OK) return;
  MsgOut(status," Say hello to Ping (())(())--");  
}


/**
 * \fn void sc2da_Report(StatusType *status)
 *
 * \brief drama action
 *  print out current status of the tasks
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * ditscmd SC2DAx REPORT
 */

/*+ sc2da_Report 
*/
void sc2da_Report
( 
StatusType *status   
)
{
  if (*status != STATUS__OK) return;

  MsgOut(status," Action %s ",seqStatus[dasInfo.actionFlag]);
}



/**
 * \fn void sc2da_Pixelmon(StatusType *status)
 *
 * \brief daram action
 *  monitor pixel continuosly until being kicked, 
 *  it reads a single frame data back, write to file, check if any onflyCmd 
 *  is recevived, if so, sends the cmd.   repeats this 
 *  the dataslit_flag (=1) is used to slit error signal from data
 *
 * \param *status StatusType.  given and return
 *
 * 
 * ditscmd SC2DAx PIXELMONITOR \
 *               DATA_FILE=xxx \
 *               SETUPBATCH_FILE=pixelmon.txt\
 *               DATASLIT_FLAG=x \ 
 *               SQ2OPTFILE=xxxx
 *
 */
/*+ sc2da_Pixelmon 
*/
void sc2da_Pixelmon
( 
StatusType *status    
)
{
  int           rval, wait;
  uint32        *longword;
  char           *glbmsgPtr,*lclmsgPtr; 
  static DRAMA_INNERMSG    dramamsg;
  static uint32 mceBufsize;
  static char   *byte;
  static dasCmdInfo_t pixelCmd;
  static ARRAYSET  pixelSet;
  static char    dateTime[40];


  if (*status != STATUS__OK) return;

  if(DitsGetSeq()==0)
  {
    rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
    // initial dasInfo.cmdFlag=0;
    sc2dalib_pixelmonInit(con,&dasInfo,&pixelCmd,&glbMceInx,dateTime,&pixelSet,status);
    if (*status != STATUS__OK)
    {
      ErsRep(0,status,"sc2da_Pixelmon: sc2dalib_pixelmonInit failed");
      sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
      return;
    }  
  
    /*  here comes acknowledge to RTSC
    */
     DitsTrigger(0,status);
    if (*status != STATUS__OK) 
    {
      ErsRep(0,status, "sc2da_Pixelmon: DitsTrigger failed");
      sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
      return;
    }
    MsgOut(status,
     " Start monitoring ( use stoppixelmonitor to stop it) ............");
     

    // dasInfo.msgwrtPt=0; dasInfo.msgreadPt=0; dasInfo.trkNo=0;
    // after call sc2dalib_frametakeInit
    sc2dalib_frametakeInit(con,&dasInfo,1,1,&glbMceInx,dateTime,status);
    if (*status != STATUS__OK) 
    {
      ErsRep(0,status, "sc2da_Pixelmon: sc2dalib_frametakeInit failed");
      sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
      return;
    }
    con->process.framesetup=DA_MCE_SINGLEDATA;  
    // for data frame now
    pixelSet.fbFlag=0;
    dasInfo.actIndex=DitsGetActIndex();
    DitsPutRequest ( DITS_REQ_SLEEP, status );
    post_sem(con->intask_sem,0);
  }
  else
  {
    // get the glbmsg, we can use it directly later
    glbmsgPtr=(char*)&glbMsg[dasInfo.msgreadPt];
    lclmsgPtr=(char*)&dramamsg;
    memcpy(lclmsgPtr,glbmsgPtr,sizeof(DRAMA_INNERMSG));       
    // recycle=16
    dasInfo.msgreadPt++;
    dasInfo.msgreadPt=dasInfo.msgreadPt & 0x000F ; 

    if(dramamsg.reason==FRAME_COMPLETION)   
    {
      //data frame finishes, no more expected 
      byte=(char*)dramamsg.data;
      longword = (uint32 *)byte;
      mceBufsize=dramamsg.bufsize;
      jitDebug(2,"sc2da_Pixelmon: PacketSize from MCE=<%ld>",mceBufsize); 

      if (dasInfo.filesvFlag)
      {
        // save in trkNo  data1 data2 data2  ....
        //           1     xxx   xxx   xxxx
        //           2     xxx   xxx   xxxx
        //
        sc2dalib_savepixelData(con,&dasInfo,byte,mceBufsize/4,&pixelSet,status);
        if ( *status !=SDSU_OK)
        {
          ErsRep(0,status,"sc2da_Pixelmon: sc2dalib_savepixelData failed");
          sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
          return;
        }
      }
      
      dasInfo.trkNo++;        
      if( dasInfo.actionFlag !=PIXELMONKICK)  
      { 
        // reschedule to read again
        // check if we receive ontheFly cmd
        if ( dasInfo.cmdFlag ==1)
        {
          // send the cmd, then reset cmdflag
          sc2dalib_sendCmd(con,&dasInfo,&onflyCmd,&glbMceInx,dateTime,status);
          dasInfo.cmdFlag=0;
          if ( !StatusOkP(status) )
          { 
            ErsRep(0,status,"sc2da_Pixelmon: sc2dalib_sendCmd %s failed",
                   onflyCmd.mceCmd); 
            sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
            return;
          }
        }

        if ( dasInfo.trkNo >=SQ2FB_UPDATE_WAIT && pixelSet.slopSelect[5]>0 )
        {
          sc2dalib_changeSQVal(con,&dasInfo,&glbMceInx,dateTime,&pixelSet,status);
          if ( !StatusOkP(status) )
          { 
            ErsRep(0,status,"sc2da_Pixelmon: sc2dalib_changeSQVal failed"); 
            sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
            return;
          }
        }

        // only go is needed here, add delay here  14/02/07 xg
        // allow  fraction of second if  slopSelect[7] > 0
        if( pixelSet.slopSelect[7] > 0 )
        {
          wait=1000000*pixelSet.slopSelect[6]/pixelSet.slopSelect[7]; 
          usleep(wait);
        }


        *status=sdsu_command_mce(con,pixelCmd.cmdBuf,&pixelCmd.reply);
        if ( *status !=SDSU_OK)
        {
          ErsRep(0,status,"sc2da_Pixelmon: command_mce failed");
          sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
          return;
        }
        else
        {
          DitsPutRequest ( DITS_REQ_SLEEP, status );
          post_sem(con->intask_sem,0); 
        }    
      }
      else
      {
        if(dasInfo.debuglvl==2)    
        {
          sc2da_Dispinfo(status);
        }
        sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
        jitDebug(4,"sc2da_Pixelmon: the action has completed");
      }
    }
    else //
    {
      con->process.seqstatus=SEQ_ERROR; 
      *status=DITS__APP_ERROR;
      sc2dalib_msgprintSave(&dasInfo,
               "sc2da_Pixelmon: Error: %s",dramamsg.errRep,USE_ERSREP,status);
      sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
    }
  }
}



/**
 * \fn void sc2da_pixelKick(StatusType *status)
 *
 * \brief drama action
 *  kick sc2da_Pixelmon,set dasInfo.actionFlag=PIXELMONKICK  
 *  and use DitsSignalByIndex to end sc2da_Pixelmon
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * >ditscmd SC2DAx -k PIXELMONITOR
 */
/*+ sc2da_pixelKick - 
*/
void sc2da_pixelKick
( 
StatusType *status    
)
{
  char    dateTime[40];
  int     rval;
  int     msgreadPt;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> CMD from sc2da_pixelKick \n",
          dateTime);

  dasInfo.actionFlag=PIXELMONKICK; 

  // get the glbmsg, we need to use msgreadPt-- 
  if (dasInfo.msgreadPt !=0)
    msgreadPt=dasInfo.msgreadPt-1;
  else
    msgreadPt=15;

  if (glbMsg[msgreadPt].reason !=FRAME_COMPLETION)
  {
    if ( dasInfo.actIndex != -1 )
    {
      DitsSignalByIndex ( dasInfo.actIndex, 0, status );
    }
    else
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,"sc2da_pixelKick: not valid dasInfo.actIndex");
    }  
  }
  jitDebug(4,"sc2da_pixelKick: Kick has completed");
}


/**
 * \fn void sc2da_Pixelmask(StatusType *status)
 *
 * \brief drama action
 *  read pixel Mask and apply it to mce-gain-lock-sq1fb 
 *
 * \param *status StatusType.  given and return
 *
 * ditscmd SC2DAx PIXEL_MASK
 * 
 */
/*+ sc2da_Pixelmask -
*/
void sc2da_Pixelmask
(
StatusType *status
)
{
  int          i, j, pixel;
  SdsIdType    argId=0;
  SdsIdType    maskId=0;
  char         name[FILE_LEN];

  if (!StatusOkP(status)) return;

  // the pixealmask.xml need to be in $CONFIG_HARD
  argId = DitsGetArgument();

  jitArgGetXML ( argId, "PIXEL_MASK", 1, NULL, sizeof(name), name,
                 &maskId, status );
  if (!StatusOkP(status))
  {
   ErsRep(0,status, "sc2da_Pixelmask: jitArgGetXML failed");
   return;
  }
  MsgOut(status," pixelMaskXML=%s",name);

  sc2dalib_readpixelMask(&dasInfo,pixelMask,maskId,status);

  SdsFreeId ( maskId, status );
  SdsFreeId ( argId, status );

  if (!StatusOkP(status))    return;

  if (dasInfo.fpLog != NULL)
  {
    fprintf(dasInfo.fpLog,"===== PIXEL_MASK ARRAY ====\n");
    for ( i=0; i<41; i++ )
    {
      for ( j=0; j<32; j++ )
      {
        pixel= i*32 +j;
        fprintf(dasInfo.fpLog,"%3d", pixelMask[pixel]);
      }    
      fprintf(dasInfo.fpLog,"\n");
    }
    fprintf(dasInfo.fpLog,"===== PIXEL_MASK ARRAY ====\n");
    fflush(dasInfo.fpLog);
    fflush(dasInfo.fpLog);
  }
  MsgOut(status," sc2da_Pixelmask OK");
}


/**
 * \fn void sc2da_SetSeq(StatusType *status)
 *
 * \brief drama action
 *  setup SEQUENCE for DA
 *
 * \param *status StatusType.  given and return
 *
 *
 * The SETUP_SEQUENCE action checks the flags to ensure that the 
 * SCUBA2 DAS is initialised and configured, and it is not currently
 * executing a sequence. If the checks are satisfactory, it parses 
 * the state table file. The flag recording that a setup sequence 
 * is done is set.
 *
 *  ditscmd $SC2DRAMATASK SETUP_SEQUENCE  \
 *             GROUP=1234 \
 *             TASKS="RTS bcd cde" \
 *             SOURCE=SCIENCE \
 *             LOAD=$load \
 *             BB_TEMP=xx HEAT_CUR=  SHUT_FRAC   \
 *
 */

/*+ sc2da_SetSeq     
*/
void sc2da_SetSeq
(
StatusType *status
)
{

  char       dateTime[FILE_LEN];
  int        rval;
  static struct sc2headman_par  headPar;
  SC2STORETelpar telpar;
  PAR_SHARED *parshmPtr;

  if (!StatusOkP(status)) return;

  /* If the last time we heater tracked it failed, then retrun right away with an error */
  if (dasInfo.heatTrackFailed == 1)
    {
      *status = DITS__APP_ERROR;
      ErsRep(0, status, "sc2da_SetSeq: Heater tracking has failed please setup before continuing");
      return;
    }

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  sc2dalib_setseqInit(con,&dasInfo,&glbMceInx,dateTime,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2da_SetSeq: failed after call to sc2dalib_setseqInit");
    return;
  }
  if (dasInfo.engFlag==RTSC_MODE)
  {
     if ( strcmp(domainName,"ROE.LINUX") ==0 )
       MsgOut(status,"sc2da_SetSeq: It is not at JAC, don't call sc2headman_setup");
     else
     {
       jitDebug(2,"sc2da_setSeq: call sc2headman_setup\n");
       sc2headman_setup(status);
     }
     if ( !StatusOkP(status) )
       ErsRep(0,status,"sc2da_SetSeq: sc2headman_setup failed");
     else
     {
       headPar.pixheat=dasInfo.pixelHeat;
       if ( strcmp(dasInfo.load, "HOT") ==0 )
       {
         dasInfo.parshmPtr->load=LOAD_HOT;
       }
       else
       {
         if ( strcmp(dasInfo.load,"SKY") ==0)
           dasInfo.parshmPtr->load=LOAD_SKY;
         else
           dasInfo.parshmPtr->load=LOAD_DARK;
       }

       headPar.darkheat = dasInfo.darkHeaterI;
       headPar.shutter = dasInfo.shutterFraction;
       headPar.seqcount = dasInfo.seqcount;
       headPar.drcontrol = dasInfo.drcontrol;
       headPar.detbias = dasInfo.detbias;
       headPar.datamode = dasInfo.datamode;
       strcpy( headPar.status,"NORMAL");
       strncpy(headPar.inbeam, dasInfo.inbeam, MAX_INBEAM);

       jitDebug(2,"sc2da_setSeq: call sc2headman_putsc2par\n");

       sc2headman_putsc2par(headPar,status);
       if ( !StatusOkP(status) )
         ErsRep(0,status,"sc2da_SetSeq: sc2headman_putsc2par failed");

      /* get telescope parameters from headman and put into shared memory */


       sc2headman_gettelpar(&telpar, status);
       parshmPtr=dasInfo.parshmPtr; 
       parshmPtr->dut1 = telpar.dut1;
       parshmPtr->tel_latdeg = telpar.latdeg;
       parshmPtr->tel_longdeg = telpar.longdeg;
       parshmPtr->instap_x = telpar.instap_x;
       parshmPtr->instap_y = telpar.instap_y;

    }
  }
  SdpPuti("SETUP", 1, status);
  sc2dalib_endAction(con,&dasInfo,status);
  // reset pci frame count=0
  sc2dalib_pcisendCmd(con,"WRITE",&dasInfo,1,0,"X",status);
  if( *status!=STATUS__OK )
  {
    ErsRep(0,status,"sc2da_SetSeq: failed after call to sc2dalib_pcisendCmd X(1)=0");
    return;
  }
  if (dasInfo.engFlag==RTSC_MODE)
  {
    sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,syncCmdEx,status); 
    sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,dvCmdEx,status); 
    // sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,selClkCmdEx,status); 
    if ( !StatusOkP(status))   
      {
        ErsRep(0,status,"sc2da_SetSeq: failed after calls to sc2dalib_setmceVal to set clock source external");
        return;
      }
    jitDebug(8 ,"sc2da_Setseq: MCE to use external clocks during SEQ ");
  }
  else // For engineering mode set the clocking all internal
    {
      sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,syncCmdInt,status); 
      sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,dvCmdInt,status); 
      // sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,selClkCmdInt,status); 
      if ( !StatusOkP(status))   
	{
	  ErsRep(0,status,"sc2da_SetSeq: failed after calls to sc2dalib_setmceVal to set clock source internal");
	  return;
	}
      jitDebug(8 ,"sc2da_Setseq: MCE to use internal clocks during SEQ ");
    }

  jitDebug(2,"sc2da_SetSeq: completed\n");
}


/**
 * \fn void sc2da_SetSeqKick(StatusType *status)
 *
 * \brief drama action
 *  kick SETUP_SEQUENCE. 
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * >ditscmd SC2DA -k SETUP_SEQUENCE 
 *
 */
/*+ sc2da_SetSeqKick 
*/
void sc2da_SetSeqKick
( 
StatusType *status    /* global status (given and returned) */
)
{
}


/**
 * \fn void sc2da_Seq(StatusType *status)
 *
 * \brief drama action
 *  start SEQUENCE
 *  send a GO command to the MCE to start sequencing operations. 
 *
 * \param *status StatusType.  given and return
 *
 *
 * The SEQUENCE action checks the flags to ensure that the SCUBA2 
 * DAS is setup. If the check is satisfactory, it deduces the number 
 * of sequence steps to be executed from the arguments and stores it
 * in shared variables. The semaphore is used to tell the dhTask to
 * do data frame reading and processing, the flag recording that a 
 * sequence is in progress is set, and an acknowledgement is sent 
 * to RTS Client before DRAMA is informed that the action expects to 
 * reschedule. When a message comes through a FIFO from the dhTask, 
 * the sc2da_Seq checks if the message is for error or completion or 
 * QL and unsets the "sequence in progress" flag before exiting.
 * 
 *  ditscmd $SC2DRAMATASK SETUP_SEQUENCE  \
 *              GROUP=1234 \
 *             TASKS="RTS" \
 *             SOURCE=SCIENCE \
 *             LOAD=$load \
 *             LOADV=" "
 *
 */
/*+ sc2da_Seq - 
*/
void sc2da_Seq
(
StatusType *status
)
{ 
  int               rval, *tmpData;
  double            seq_time;          
  char              *glbmsgPtr,*lclmsgPtr; 
  char              dateTime[40];
  DitsReasonType    entReason;
  char entPathName[DITS_C_NAMELEN];
  static time_t     tloc, newtloc;   //generally time_t is a long 
  static long       frameRate, seqStart, seqEnd;
  static int        dataDone;     /* flag for data taking completed */
  static int        trigQL;
  static int        *lookupTable; /* store endsubSeq locally from shared Memory */
  static DRAMA_INNERMSG    dramamsg;
  static char obsidss[OBSIDSS_SIZE];
  static int errorDetected;     /* To keep from hammering the MCE after the first problem is detected */
  StatusType localStatus; 

  if (!StatusOkP(status)) return;

  // Check whether it is the start or completion of a sequence 
  if(DitsGetSeq()==0)
  {
    errorDetected = 0;
    dataDone = 0;
    dasInfo.headersDone = 1;
    if (dasInfo.engFlag==RTSC_MODE) 
      dasInfo.headersDone= 0;
    jitDebug(2,"sc2da_Seq: call sc2dalib_finddateTime\n");  
    rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);

    // set the ith sequence NO we are waiting for in dasInfo.ithseqWait
    jitDebug(2,"sc2da_Seq: call sc2dalib_seqInit\n");  
    sc2dalib_seqInit(con,&dasInfo,dateTime,&lookupTable,status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2da_Seq: sc2dalib_seqInit failed"); 
      sc2dalib_seqEnd(con,&dasInfo, lookupTable,status);
      return;
    }
    // get frameRate (Hz) for later to set timeOut
    SdpGeti("FRAME_RATE", &frameRate, status);
    if(frameRate==0)
    { 
      MsgOut(status,"sc2da_Seq: frameRate =0, force it=1"); 
      frameRate=1;
    }
    
    SdpGeti("SEQ_START", &seqStart, status);
    SdpGeti("SEQ_END", &seqEnd, status);
    MsgOut(status ,"sc2da_Seq: seqStart = %ld seqEnd = %ld",seqStart, seqEnd); 
      // call sc2headman_seqstart
    if (dasInfo.engFlag==RTSC_MODE)
    {
       sc2dalib_seqcallsc2headmanseqStart(&dasInfo, seqStart, seqEnd,status);
       if ( !StatusOkP(status) )
       {
         sc2dalib_seqEnd(con,&dasInfo,lookupTable, status);
         return;
       }

       sc2headman_getobsidss(OBSIDSS_SIZE, obsidss, status);
       strncpy(dasInfo.parshmPtr->obsidss,obsidss,OBSIDSS_SIZE); 
    }
    //send No_frames parameter to MCE and a GO to MCE
    //  dasInfo.msgwrtPt=0; dasInfo.msgreadPt=0;con->datacount=end-start+1
    sc2dalib_frametakeInit(con,&dasInfo,seqStart,seqEnd,&glbMceInx,dateTime,status);
    if (!StatusOkP(status)) 
    {
      ErsRep(0,status, "sc2da_Seq: sc2dalib_frametakeInit failed");
      sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);  
      sc2dalib_seqEnd(con,&dasInfo,lookupTable,status);
      return;
    }

    SdpPuti("IN_SEQUENCE",DA_MCE_SEQ,status);
    con->process.seqstatus=SEQ_INACTION;
    con->process.framesetup=DA_MCE_SEQ;
    trigQL=0; 

    // after send startSeq, endSeq and a GO, acknowledge to RTSC
    DitsTrigger(0,status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2da_Seq: DitTrigger failed"); 
      sc2dalib_seqEnd(con,&dasInfo,lookupTable, status);
      sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);  
      return;
    }

    //get time in seconds since EPOCH 
    time(&tloc);         
    dasInfo.actIndex=DitsGetActIndex();
    DitsPutRequest ( DITS_REQ_SLEEP, status );
    post_sem(con->intask_sem,0);
  }
  else
  {

    //DITS_REA_OBEY Obey message received
    //DITS_REA_KICK Kick message received
    //DITS_REA_COMPLETE Message completed
    //DITS_REA_RESCHED Action reschedule by timer expiry
    //DITS_REA_TRIGGER Action triggered by subsidiary action
    //DITS_REA_ASTINT Action triggered by signal handler
    entReason=DitsGetEntReason();
    switch( entReason )
    {    
      // rescheduled by call-back routine for FRAME, 
      // check dramamsg.reason   
      case DITS_REA_ASTINT:
      {
        // get the glbmsg, we can use it directly later
        glbmsgPtr=(char*)&glbMsg[dasInfo.msgreadPt];
        lclmsgPtr=(char*)&dramamsg;
        memcpy(lclmsgPtr,glbmsgPtr,sizeof(DRAMA_INNERMSG));       

        // recycle=16
        dasInfo.msgreadPt++;
        dasInfo.msgreadPt=dasInfo.msgreadPt & 0x000F ; 

        if( dramamsg.reason==FRAME_TRIGQL  || dramamsg.reason==FRAME_SCANREADY)
        {
          trigQL++;
          jitDebug(2,"sc2da_Seq : qlframeNum= %d\n", dasInfo.parshmPtr->qlframeNum);

          // tell seqchkQL that frame ready for parameters check or QL
          dasInfo.headtaskchkFlag=1; 
          jitDebug(2,"sc2da_Seq: received %d-th QL trig\n",trigQL);
          sc2dalib_seqchkQL(con,&dramamsg,&dasInfo,entReason,lookupTable,status);
          if ( !StatusOkP(status) )
          {
             ErsRep(0,status,"sc2da_Seq: sc2dalib_seqchkQL failed\n");
             ErsFlush(status);
             sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);
             return;
          }
          DitsPutRequest(DITS_REQ_SLEEP, status);
          return;
        }
        else if (dramamsg.reason==FRAME_CHKSUMWRONG)
        { 
          // carry on if checksum is wrong 1 use ErsOut,*status=STATUS__OK;  
          sc2dalib_msgprintSave(&dasInfo,"CHKSUMWRONG %s",dramamsg.errRep,USE_ERSOUT,status);
          *status=STATUS__OK;  
          // dhtask is going to use read_data_timed after first frame,
          // this is OK now without timeout
          DitsPutRequest(DITS_REQ_SLEEP, status);
          return;
        }
	/* We are finished with this step of heater/bias ramp/sawtooth so take another step */
	else if (dramamsg.reason==FRAME_SAW_RAMP && dasInfo.sawtoothRampFlag > 0)
	  {

	    /* If I have already had an error just go back and wait for the next message */
	    if(errorDetected)
	      {
		DitsPutRequest(DITS_REQ_SLEEP, status);
		return;
	      }

	    /* Is it a heater current something? */
	    if(dasInfo.sawtoothRampFlag < 3)
	      {
		sc2dalib_stepHeaterCurrent(con, &dasInfo, &glbMceInx, status);
	      }
	    else /* It is a bias something */
	      {
		sc2dalib_stepTESBias(con, &dasInfo, &glbMceInx, status);
	      }

	    /* If we got an MCE error, do about what we do when we get kicked */
	    if(*status != STATUS__OK)
	      {

		errorDetected = 1;

		/* First annul the error because nothing else can happen with bad error status */
		ErsAnnul( status );

		localStatus = DITS__APP_ERROR;
		ErsRep(0, &localStatus,"Classic MCE error detected - will ignore for now"); 

		/* Now abort the sequence
		   sc2da_Abort(status); */

		/* Request to be re-run again this time with abort status
		DitsPutRequest ( DITS_REQ_STAGE, status );
		return; */
	      }

	    DitsPutRequest(DITS_REQ_SLEEP, status);
	    return;
	  }
        else if (dramamsg.reason==FRAME_SQ2FBCMP && dasInfo.gotsq2fbparaFlag==1)
        { 
          // dasInfo.gotsq2fbparaFlag and myInfo->parshmPtr->sq2fbparaFlag are set in setseqInit
          // 
          // in RTSC mode, SQ2FB_GETPARA action will have myInfo->parshmPtr->sq2fbparaFlag=1 for SEQ
          tmpData=(int*)dramamsg.data;
          jitDebug(2,"sc2da_Seq : updatesq2fb frameNum= %d\n", tmpData);

           // tmpData[0] is frameNo, inside changesq2fb, if we use "param-read bc2 flux_fb 32"
          // it will tell us if all change to sq2fb is right, as it writes 
          // for (j=0;j<32;j++) val=frameNo+j ; to sq2fb
      
          sc2dalib_changesq2fbVal(con,&dasInfo,&glbMceInx,dateTime,tmpData,status);
/*
          //  need to use dasInfo.trkNo(=0 when DitGetSeq()=0)
          //  to set zfact[i]=data[i] if trkNo=0;
          sc2dalib_updatesq2fbVal(con,&dasInfo,&glbMceInx,dateTime,dramamsg.data, &setup0,
                                  setup0.sq2fdbkOpt,status);
*/
          if ( !StatusOkP(status) )
          { 
            ErsRep(0,status,"sc2da_seq: sc2dalib_updatesq2fbVal failed"); 

           // status is saved and returned
           //sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);  

            //  check and display msg, postsem(intask)
            *status=STATUS__OK;
            sc2dalib_seqchkEnd(con,&dramamsg,&dasInfo,status);
            wait_sem(con->inner_sem,0);          
            sc2dalib_seqEnd(con,&dasInfo,lookupTable,status);
            return;
          }
          dasInfo.trkNo ++;
          DitsPutRequest(DITS_REQ_SLEEP, status);
          return;
        }

        else
        {
	  dataDone = 1;
          jitDebug(2,"sc2da_Seq: received FRAME trig\n");
          if  ( dramamsg.reason ==FRAME_COMPLETION )
          {
            // if it is monitor task, headerDone is set during seqchkQL
            if (dasInfo.headersDone !=1)  // wait for header collection
            {
              MsgOut(status,"sc2da_Seq: wait for headersDone,FRAME_COMPLETION");
              DitsPutRequest(DITS_REQ_SLEEP, status);
              return;
            }
          }
          else
          {
             // also check other .reason and display msg, postsem(intask)
             sc2dalib_seqchkEnd(con,&dramamsg,&dasInfo,status);
             wait_sem(con->inner_sem,0);
             // get time in seconds since EPOCH 
             time(&tloc);
             seq_time=difftime(tloc,newtloc); 
 
             if(dramamsg.reason !=FRAME_STOPPED)
               sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);  
  
             sc2dalib_seqEnd(con,&dasInfo,lookupTable,status);
             MsgOut(status,"sc2da_seq: seq ended");
             return;
          }
        }
        break;
      }
    case DITS_REA_RESCHED: 
      {
	if( dasInfo.actionFlag == ABORTACTION) /* I was kicked and now I want to end the action */
	  {
	    int endsubSeq;
	    MsgOut(status,"sc2da_seq: ABORTACTION is TRUE - ending SEQUENCE");

	    /* Giving this semaphore lies to the dhtask telling it all of the headers are ready */
	    post_sem(con->intask_sem,0);

	    /* Now lie to ourselves in case dhtask gives FRAME_COMPLETION */
	    /*dataDone = 1; set by FRAME trigger*/
	    dasInfo.headersDone = 1;
	    dasInfo.ithseqWait = (int)seqEnd;

            /* NOTE:  The method below seems to work okay if a single task messed up,
               but leaves the DA tasks hanging if the user hits ABORT.
            */

//             /* Lie to sc2headman so it thinks all task headers are ready */
//             endsubSeq = *(lookupTable + dasInfo.headerreadNo);
//             sc2headman_settaskseq(endsubSeq, status);
// 
//             /* Lie about entry reason to seqchkQL, let it handle post_sem, dasInfo, etc. */
//             entReason = DITS_REA_ASTINT;
//             dramamsg.reason = FRAME_TRIGQL;
//             dasInfo.headtaskchkFlag=1;
//             sc2dalib_seqchkQL(con,&dramamsg,&dasInfo,entReason,lookupTable,status);
//             if ( !StatusOkP(status) )
//             {
//               ErsRep(0,status,"sc2da_Seq: sc2dalib_seqchkQL failed\n");
//               ErsFlush(status);
//               /* post_sem(con->intask_sem,0) ??? */
//               return;
//             }

	    /* Go back to sleep and wait for another trigger */
	    DitsPutRequest(DITS_REQ_SLEEP, status);
	    return;
	  }
	else
	  {
	    *status = DITS__APP_TIMEOUT;
	    ErsRep(0, status, "Timeout waiting for frame ");
	    fprintf(dasInfo.fpLog,"sc2da_Seq: Timeout waiting for frame\n");
	    sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);  
	    sc2dalib_seqEnd(con,&dasInfo,lookupTable,status);
	    return;
	  }
	  break;
      }
      case DITS_REA_TRIGGER:  // trigged by para monitor
      {
        if(dasInfo.engFlag==RTSC_MODE)
	{ 
           jitDebug(2,"sc2da_Seq: call sc2headman_trighandle =%d\n",
	          DITS_REA_TRIGGER );
           sc2headman_trighandle(status);
           if ( !StatusOkP(status) )  
           {
             ErsRep(0,status,"sc2da_Seq: sc2headman_trighandle failed");
             sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);  
             return;  
	   } 
           // check if all monitored tasks STATEs are ready and do QL
           // use sc2headman_checktasks,  also check if last parameter is ready 
           jitDebug(2,"sc2da_Seq: check if QL headers ready, trigged by paramonitor\n" );
           sc2dalib_seqchkQL(con,&dramamsg,&dasInfo,entReason,lookupTable,status); 
           if ( !StatusOkP(status) )  
           {
             ErsRep(0,status,"sc2da_Seq: sc2dalib_seqchkQL failed parMonitTrig");
             sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);
             return;  
	   } 
	 
           //also consider if no monitor-task
           if (dasInfo.headersDone !=1  ||  dasInfo.ithseqWait !=(int)seqEnd ||  dataDone !=1  ) 
	   {
              DitsPutRequest(DITS_REQ_SLEEP, status);
              return;
           }
           MsgOut(status,"sc2da_seq: this is last REA_TRIGGER, ithseq=%d",dasInfo.ithseqWait);
        }
        break;
      }
      default:
      {
        MsgOut(status,"sc2da_seq: this is default in switch");
        // Something strange has happened.
	/* Try to get the name of the task from which this message was sent */
	if((int)DitsGetEntPath() > 0)
	  DitsTaskFromPath(DitsGetEntPath(), DITS_C_NAMELEN, entPathName, status);
	else
	{
	  if((int)DitsGetParentPath() > 0)
	    DitsTaskFromPath(DitsGetParentPath(), DITS_C_NAMELEN, entPathName, status);
	  else
	    strcpy(entPathName,"Path 0");
	}
        DitsPrintReason(entReason,DitsGetEntStatus(),status);
        *status = DitsGetEntStatus();
        ErsRep(0, status, "Unexpected entry: msg from task %s entry status %s",entPathName, DitsErrorText(*status));
        fprintf(dasInfo.fpLog,
                "sc2da_Seq:Unexpected entry: msg from task %s entry status %s\n",entPathName, DitsErrorText(*status));
        sc2dalib_stopFrame(con,&dasInfo,&glbMceInx,dateTime,status);  
        sc2dalib_seqEnd(con,&dasInfo,lookupTable,status);
        return;
      }
    }
    // in case par-monitor is later than dhtask, also consider if no monitor-task
    if (dasInfo.headersDone==1 && dataDone==1 &&  dasInfo.ithseqWait==(int)seqEnd ) 
    {
 
       jitDebug(2,"sc2da_seq: headersDone=%d dataDine= %d, ithseqWait=%d\n",
               dasInfo.headersDone, dataDone, dasInfo.ithseqWait);

       jitDebug(8,"sc2da_seq: msg.reason=%s, msgreadPt=%d\n",
               seqStatus[dramamsg.reason],(int)dasInfo.msgreadPt);
 
       if(*status==STATUS__OK)
       {
         // get time in seconds since EPOCH 
         time(&newtloc);
         seq_time=difftime(newtloc,tloc); 
         jitDebug(2,"sc2da_Seq: The seq_time ~= <%f>s \n",seq_time );

         if (dasInfo.engFlag==RTSC_MODE)
         { // no need for SCAN, it is done inside seqchkQL
           if (dasInfo.parshmPtr->obsMode !=OBS_SCAN)
           {
              sc2headman_mainhead(dasInfo.parshmPtr->subscanNo,(int)seqStart, 
				  (int)seqEnd, dasInfo.utcshort, MAXFITS, 
				  &dasInfo.parshmPtr->actfits, 
				  dasInfo.parshmPtr->fitshd,status);

             if ( !StatusOkP(status) )
                 ErsRep(0,status,"sc2da_seq: sc2headman_mainhead failed"); 

             MsgOut(status,"sc2da_seq: wait for writing SDF file to complete...(^.^)\n");
           }
         }
         // also check other .reason and display msg, postsem(intask) for 
         // dhtask to end current frame taking
         sc2dalib_seqchkEnd(con,&dramamsg,&dasInfo,status);
         wait_sem(con->inner_sem,0);
         // get time in seconds since EPOCH 
         time(&tloc);
         seq_time=difftime(tloc,newtloc); 
         //printf("sc2da_Seq: writing NDF_time ~= <%f>s \n",seq_time );
         sc2dalib_seqEnd(con,&dasInfo,lookupTable,status);
         MsgOut(status,"sc2da_seq: seq completed"); 
         return;
      }
      else
      {
         ErsRep(0,status, "sc2da_seq:seq completed with bad status");
         sc2dalib_seqEnd(con,&dasInfo,lookupTable,status);
         return;
      }
    }
    else
    {
      jitDebug(8,"sc2da_seq: hedaersDone=%d  dataDone= %d, ithseqwait=%d\n",
         dasInfo.headersDone, dataDone, dasInfo.ithseqWait);
    }
  }
}




/**
 * \fn void sc2da_SeqKick(StatusType *status)
 *
 * \brief darama action
 *  kick SEQUENCE
 *  send a STOP command to the MCE to stop sequencing operations. 
 *
 * \param *status StatusType.  given and return
 *
 *
 * kick of the SEQUENCE action by ending the action. The MCE will send
 * back the last frame with stop and last bits set
 *
 *
 * ditscmd SC2DA -k SEQUENCE 
 *
 */
/*+ sc2da_SeqKick -
 */
void sc2da_SeqKick
(
StatusType *status
)
{  
  char    dateTime[40];
  int     rval;

  if (!StatusOkP(status)) return;

  dasInfo.actionFlag=SEQKICK; 

  /* Tag the data as having a problem
  strcpy( headPar.status,"ERROR");
  sc2headman_putsc2par(headPar,status); */

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> CMD from sc2da_SeqKick\n",dateTime);
  sc2da_Abort(status);  
  if ( !StatusOkP(status) )  
  {
    ErsRep(0,status,"sc2da_SeqKick: sc2da_Abort failed");
    return;  
  }
  jitDebug(4,"sc2da_SeqKick: Kick has completed\n");

  DitsPutRequest ( DITS_REQ_STAGE, status );

}


/**
 * \fn void sc2da_Servo(StatusType *status)
 *
 * \brief drama action
 *  for array setup
 *
 *  setdebug 16 to display some message related to function calls
 *
 * \param *status StatusType.  given and return
 *
 *
 * ditscmd SC2DA SERVO \
 *               SVFILE_FLAG=$svfileFlag \
 *               DATA_FILE=$CURRENTDATADIR/$datafile \
 *               SETUPBATCH_FILE=$setupfile \
 *               STRCHART_FILE=$ORAC_DATA_OUT/$strchartfile \
 *               LOCKPTS_FILE=$CURRENTDATADIR/$lockptsfile \      
 *               SQ1_FLAG=x \
 *               SEL_ROW=x
 */

/*+ sc2da_Servo - 
*/
void sc2da_Servo
( 
StatusType *status    /* global status (given and returned) */
)
{ 
  int                   rval, seq;
  char                  *byte;
  double                seq_time;
  char                  tmp[FILE_LEN];
  char                  *glbmsgPtr,*lclmsgPtr; 
  DRAMA_INNERMSG        dramamsg;
  static dasCmdInfo_t   servoCmd;
  static ARRAYSET       setup;
  static char           dateTime[40];
  static char           *dataBuf;  //memory buffer for servodata
  static ARRAYSET       ssalckset;
  static char           *ssalckdataBuf;  //memory buffer for ssalock
  static time_t          tloc, newtloc;     
  static int             isEnd;   //indicator used for sc2daservoFunc
  long                  heatFlag;

  if (!StatusOkP(status)) return;

  //  check if it is start or completion.
  if(DitsGetSeq()==0)
  {
    dataBuf=NULL;
    ssalckdataBuf=NULL;
    rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
    jitDebug(16,"sc2da_Servo: call sc2daservoInit\n");     
    sc2daservoInit(con,&dasInfo,&servoCmd,&glbMceInx,&setup,
                   dateTime,&dataBuf,&ssalckset,&ssalckdataBuf, status);
    if ( !StatusOkP(status) )
    {
      jitDebug(16,"sc2da_Servo: sc2daservoInit failed:\n");     
      ErsRep(0,status,"sc2da_Servo: sc2daservoInit failed"); 
      // fpData may not be opened yet,
      sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,1,status);
      return;
    }

    /* This routine will NOT work on external clocks because there is
       no DV from the RTS when it is being used, so set the clocks all
       internal */

    sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,syncCmdInt,status); 
    sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,dvCmdInt,status); 
    // sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,selClkCmdInt,status); 
    if ( !StatusOkP(status))   
      {
        ErsRep(0,status,"sc2da_Servo: failed after calls to sc2dalib_setmceVal");
        sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,1,status);
        return;
      }

    /* After every setup we want to update the heater tracking reference values.
       Therefore I am going to set the UPDATE_HEATER_TRACKING_REFERENCE_BIT in
       the HTTRACK_FLAG. I chose SQ1OPEN because we are pretty sure it is used
       in every type of setup we currently do */
    if (setup.servo == SQ1OPEN)
      {
	SdpGeti("HTTRACK_FLAG", &heatFlag, status);
	heatFlag = UPDATE_HEATER_TRACKING_REFERENCE_BIT | heatFlag;
	SdpPuti("HTTRACK_FLAG", heatFlag, status);
      }

    if (setup.servo ==SQ2OPEN4P)
    {
      //read sq2fboptimal.txt in and store in setup.sq2fdbkOpt
      sprintf(tmp, "sq2fboptimal.txt");
      jitDebug(16,"_readOPT file=%s\n",tmp);
      sc2dalib_readOPT(&dasInfo,setup.sq2fdbkOpt,COL_NUM,tmp,status);
      if ( !StatusOkP(status))   
      {
        ErsRep(0,status,"sc2da_Servo: failed after call to sc2dalib_readOPT");
        sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,1,status);
        return;
      }
    }
    jitDebug(16,"isEnd=%d setup.servo=%s \n",isEnd,servoStatus[setup.servo]); 
    if (setup.servo !=TESTRANSIT)
      MsgOut(status," Expecting %d reschedules .....",  setup.biasNo*setup.fdbkNo);
    else
      MsgOut(status," Expecting %d reschedules .....",  setup.biasNo*(setup.fdbkNo*4+1));
    isEnd=1;     

    //get time in seconds since EPOCH 
    time(&tloc);         

    // send a single GO to start with
    sc2daservoFstGO(con,&dasInfo,1,1,&glbMceInx,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,"sc2da_Servo: sc2daservoFstGO failed");    
      sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,0,status);
      return;
    }
    dasInfo.actIndex=DitsGetActIndex();
    DitsPutRequest ( DITS_REQ_SLEEP, status );
    post_sem(con->intask_sem,0);
  }
  else /* This is not the first time this action has been activated */
  {
    // get the glbmsg, we can use it diretly later
    glbmsgPtr=(char*)&glbMsg[dasInfo.msgreadPt];
    lclmsgPtr=(char*)&dramamsg;
    memcpy(lclmsgPtr,glbmsgPtr,sizeof(DRAMA_INNERMSG));       
    // recycle=16
    dasInfo.msgreadPt++;
    dasInfo.msgreadPt=dasInfo.msgreadPt & 0x000F ; 

    seq=DitsGetSeq();
    if(dramamsg.reason==FRAME_COMPLETION)   
    {
      byte=(char*)dramamsg.data;
      dasInfo.bufSize=dramamsg.bufsize;  // byte buffersize
      con->process.whereabout=Dits_GetSeqn;
      if(setup.servo != TESTRANSIT)
      {
        sc2daservoFunc(con,byte,&dasInfo,&setup,&glbMceInx,
                     &isEnd,dataBuf,dateTime,&ssalckset,ssalckdataBuf,status);
        if ( !StatusOkP(status) )
        {
          ErsRep(0,status,"sc2da_Servo: sc2daservoFunc failed"); 
          sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,0,status);
          return;
        }
      }
      else
      {
        sc2datransitFunc(con,byte,&dasInfo,&setup,&glbMceInx,
                     &isEnd,dataBuf,dateTime,status);
        if ( !StatusOkP(status) )
        {
          ErsRep(0,status,"sc2da_Servo: sc2datransitFunc failed"); 
          sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,0,status);
          return;
        }
      }
      //check if we have more to do or be kicked off
      if( dasInfo.actionFlag !=SERVOKICK)  
      { 
        //check if we have more to do
        if ( seq % 20==0)
           jitDebug(16,"(seq-%d) isEnd=%d \n",seq,isEnd); 
  
        if(isEnd==0)
        {
          *status=sdsu_command_mce (con,servoCmd.cmdBuf,&servoCmd.reply);
          if( *status !=SDSU_OK )
          { 
            sc2dalib_mceerrRep(&dasInfo,&servoCmd,status);
            *status=DITS__APP_ERROR;         
            ErsRep(0, status,"sc2da_Servo: sc2dalib_mceerrRep failed"); 
            sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,0,status);          
            return;
          }
          if ( seq % 20 ==0)
            jitDebug(16,"send another GO cmd and wait for data\n"); 
          dasInfo.actIndex=DitsGetActIndex();
          DitsPutRequest ( DITS_REQ_SLEEP, status );
          post_sem(con->intask_sem,0);
        }
        else /* isEnd is not == 0 */
        {
          // add  HEATERSERVO and TESBIASSERVO SQ1BIASSERVO but
          // no findlockPts for them
          // SQ1BIASSERVO 11-01-2007 . X Gao

          if(setup.servo==SSARAMP || setup.servo ==SQ2SERVO || setup.servo ==SQ1LOCK ||
             setup.servo==SSALOCK || setup.servo ==SQ2LOCK  || setup.servo ==SQ1SERVO )
          {
            // look for the max modulation/ lock points  
            // already called findflatbottom if servo==SSARAMP
            // servo data is saved in findlockPts
            // fpMcecmd is used to save the lock points ( file=dataFile-lck)
            // dataFile is used for dataFileservodata.hex in servosavedata()
            jitDebug(16,"_findlockPts\n"); 
            sc2daliblckpts_findlockPts(&dasInfo,&setup,dataBuf,dasInfo.fpMcecmd,dasInfo.dataFile,status);
            if ( !StatusOkP(status) )
            {
              ErsRep (0, status,"sc2da_Servo: sc2daliblckpts_findlockPts failed");    
              sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,0,status);
              return;
            }
          }
          else if(setup.servo==SQ2OPEN ||  setup.servo ==SQ2OPEN4P)
          {
            // find lock points, 
            jitDebug(16,"_findlckptsGain\n");
            // find gain at the sq2fb/sq1fb from SQ2OPEN  
            sc2daliblckpts_findlckptsGain(&dasInfo,&setup,dataBuf,dasInfo.fpMcecmd,dasInfo.dataFile,pixelMask,1,status);
          }
          else if(setup.servo ==SQ1OPEN )
          {
            // find lock points,this is different from others, we have to go
            // through frame for each pixel 
            jitDebug(16,"_findframelckptsGian\n");
            // find gain at the sq2fb/sq1fb from SQ1OPEN
            sc2daliblckpts_findframelckptsGain(&dasInfo,&setup,pixelMask,status);
            jitDebug(16,"sc2dasaveframelockPts\n");
            sc2dasaveframelockPts(&dasInfo,&setup,status);
          }
          else if(setup.servo==SQ2BIASING )
          {
            // save the servo data so that findlockpoints can use it.
            sc2dasaveservoData(&setup,dataBuf,dasInfo.dataFile,status);
            MsgOut(status,
             "sc2da_Servo: sq2biasing doesn't use sc2dafindlock, max, initlockPts"); 
          }
          else if(setup.servo==SSARAMP1 || setup.servo ==CABLECAL)
          { 
            // look for mean value and use these mean value to do least square 
            // fit. take the followings off, using findlockpoints 
            // save the servo data so that findlockpoints can use it.

            sc2dasaveservoData(&setup,dataBuf,dasInfo.dataFile,status);

            jitDebug(16,"sc2dafindmeanVal\n"); 
            sc2dafindmeanVal(&dasInfo,&setup,dataBuf,dasInfo.fpMcecmd,
                           dasInfo.dataFile,status);
            if ( !StatusOkP(status) )
            {
              ErsRep (0, status,"sc2da_Servo: sc2dafindmeanVal failed");    
              sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,0,status);
              return;
            }
            jitDebug(16,"sc2daleastsquareFit\n"); 
            sc2daleastsquareFit(&setup,dataBuf,status);
            jitDebug(16,"sc2dasavelockPts\n");
            sc2dasavelockPts(&dasInfo,&setup,dataBuf,dasInfo.fpMcecmd,status);
            jitDebug(16,"sc2dasave4nextStep\n");
            sc2dasave4nextStep(&dasInfo,&setup,dataBuf,dasInfo.fpMcecmd,status);
          }    
          else if(setup.servo==TESTRANSIT )
          {
            // save the servo data so that findlockpoints can use it.
            sc2dafindTransit(&dasInfo,&setup,dataBuf,status);
            MsgOut(status,
              "sc2da_Servo: test testransit, use sc2dafindTransit"); 
          }
          else if (setup.servo==CABLECORET) 
          {
             // append ssabiaslock.txt to mcelock-ssa.txt  
             sc2dassabias4mcelockSSA(&dasInfo,&setup,status);
          }
          if ( setup.slopSelect[15] >=0 &&
               (setup.servo==SSARAMP || setup.servo ==SQ2SERVO || setup.servo ==SQ1LOCK ||
                setup.servo==SSALOCK || setup.servo ==SQ2LOCK  || setup.servo ==SQ1SERVO||
                setup.servo==SQ2OPEN || setup.servo ==SQ2OPEN4P || setup.servo ==SQ1OPEN ) 
             )
          { 
            // need to save  filetered data
             sc2daliblckpts_savefilteredData (&dasInfo,&setup,status);
          } 

          // get time in seconds since EPOCH 
          time(&newtloc);
          seq_time=difftime(newtloc,tloc); 
          rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);   
          fprintf(dasInfo.fpLog,"\n<%s> CMD from sc2da_Servo end (1)\n",dateTime);
	  /*          fprintf(dasInfo.fpLog,"Totaltime ~= <%f>s  DitsGetSeq()=%ld\n",
		      seq_time,DitsGetSeq()); */
          sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,0,status);
          return;
        }
      }
      else
      {
        MsgOut(status,"sc2da_Servo: the action is kicked off"); 
        // get time in seconds since EPOCH 
        time(&newtloc);
        seq_time=difftime(newtloc,tloc); 
        rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);   
        fprintf(dasInfo.fpLog,"\n<%s> CMD from sc2da_Servo end (2)\n",dateTime);
        /* fprintf(dasInfo.fpLog,"Totaltime ~= <%f>s  DitsGetSeq()=%ld\n",
	   seq_time,DitsGetSeq()); */
        sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,0,status);
        // try to force SERVO return bad status so that servo script loop 
        // can be terminated ???
        *status=DITS__APP_ERROR; 
        return;
      }
    } /* dramamsg.reason is NOT == FRAME_COMPLETION */   
    else
    {
      *status=DITS__APP_ERROR;
      if(dramamsg.reason==FRAME_ERRGETPASSBUF)
      {
        sc2dalib_msgprintSave(&dasInfo,
          "sc2da_Servo: Ended with TIMEOUT for data buffer","",USE_ERSREP,status);
      }
      else if (dramamsg.reason==FRAME_CHKSUMWRONG)
      {
        sc2dalib_msgprintSave(&dasInfo,
        "sc2da_Servo: Ended with CHKSUMWRONG %s",dramamsg.errRep,USE_ERSREP,status);
      }
      else 
        sc2dalib_msgprintSave(&dasInfo,
	    "sc2da_Servo: %s",dramamsg.errRep,USE_ERSREP,status);

      sc2daservoEnd(con,&dasInfo,&setup,dataBuf,ssalckdataBuf,0,status);
      con->process.seqstatus=SEQ_ERROR; 
      return;
    }
  }
}


/**
 * \fn void sc2da_servoKick(StatusType *status)
 *
 * \brief drama action
 *  kick sc2da_servo,set dasInfo.actionFlag=SERVOKICK  
 *  and use DitsSignalByIndex to end sc2da_servo
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * >ditscmd SC2DA -k SERVO 
 */
/*+ sc2da_servoKick - 
*/
void sc2da_servoKick
( 
StatusType *status    
)
{
  char    dateTime[40];
  int     rval;
  int     msgreadPt;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> CMD from sc2da_servoKick \n",
          dateTime);
  dasInfo.actionFlag=SERVOKICK; 

  // get the glbmsg, we need to use msgreadPt-- 
  if (dasInfo.msgreadPt !=0)
    msgreadPt=dasInfo.msgreadPt-1;
  else
    msgreadPt=15;

  if (glbMsg[msgreadPt].reason !=FRAME_COMPLETION)
  {
    if ( dasInfo.actIndex != -1 )
    {
      DitsSignalByIndex ( dasInfo.actIndex, 0, status );
    }
    else
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,"sc2da_servoKick: not valid dasInfo.actIndex");
    }  
  }
  jitDebug(4,"sc2da_servoKick: Kick has completed");
}

/**
 * \fn void sc2da_SetengData( StatusType *status)
 *
 * \brief function
 *  get args for SETENGDATA action and set data format and save-engdata flag values 
 *
 * \param status       StatusType     
 *
 * ditscmd SC2DA  SETENGDATA  \
 *         DATA_FORMAT=BINARY \
 *         FILESAVE_FLAG=xx
 */
/*+ sc2da_SetengData
*/
void sc2da_SetengData
(
StatusType            *status
)
{
  long       in_sequence, range[] = { 0, 3};
  long       svflag;
  char       dataform[FILE_LEN];
  SdsIdType argId;

  if (*status != STATUS__OK) return;

  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status, 
      "sc2dasetengData: %s has not completed",seqStatus[dasInfo.actionFlag]);
    return;
  } 
  dasInfo.actionFlag=SETENGDATA;

  argId = DitsGetArgument();

  jitArgGetS(argId,
          SC2DADATAFORMAT,1, NULL, "BINARY",GIT_M_ARG_UPPER, FILE_LEN, dataform, 
          NULL,status );
  if ( !StatusOkP(status) )
  {
    ErsOut(0,status, 
       "sc2da_SetengData: DATA_FORMAT_NAME is not given, use BINARY as default");
    strcpy(dataform,"BINARY");
  }
  // inform sc2dadh task to sava data
  jitArgGetI( argId, "FILESAVE_FLAG", 2, range, 0, 0, &svflag, status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2da_SetengData: failed to get FILESAVE_FLAG");
    return;
  }
  con->process.svflag=(int)svflag;

  if ( strcmp(dataform,"TEXT")==0 )
  {
    dasInfo.dataFormat=MCE_TEXT_FORM;
    con->dataform=MCE_TEXT_FORM;
  }
  else if ( strcmp(dataform,"TEXT2")==0 )
  {
    dasInfo.dataFormat=MCE_TEXT_FORM2;
    con->dataform=MCE_TEXT_FORM2;
  }
  else
  {
    dasInfo.dataFormat=con->dataform=MCE_BINARY_FORM;
  }
  SdpPutString(SC2DADATAFORMAT,dataform,status);
 }

/**
 * \fn void sc2da_SQ1optpts( StatusType *status)
 *
 * \brief function
 *  find sq1 optimal points
 *
 * \param status       StatusType     
 *
 * ditscmd SC2DA  SQ1OPT_POINTS  \
 *                SQ1OPTPTS=sq1opt-points \
 *                RESULT_FILE=xx   \
 *                FLAG=xx
 */
/*+ sc2da_SQ1optpts
*/
void sc2da_SQ1optpts
(
StatusType *status
)
{                  
  int  flag;
  int sq1bias[COL_NUM*ROW_NUM];
  int sq2feedback[COL_NUM*ROW_NUM]; 
  int sq1refbias[COL_NUM*ROW_NUM];
  int sq2reffeedback[COL_NUM*ROW_NUM];

  /* scale factor relating SQ1 bias and SQ2 feedback */
  double sq1scale[ROW_NUM*COL_NUM];
  int bqual[ROW_NUM];         /* quality for each SQ1 bias row */
  int fqual[COL_NUM];         /* quality for each SQ2 feedback column */
  int squal[ROW_NUM*COL_NUM]; /* quality for each SQ1 */
  int sq1bcomp[ROW_NUM];      /* compromise SQ1 bias for each row */
  int sq2fcomp[COL_NUM];     /* compromise SQ2 feedback for each column */
  static ARRAYSET     setup;

  if (!StatusOkP(status)) return;

  // inputfile ==>dasInfo.batch, resultfile==>dasInfo.datafile
  sc2dalib_sq1optptsInit(&dasInfo, &flag, status);
  if (!StatusOkP(status))
  {
     ErsRep(0,status,"sc2da_SQ1optpts: _sq1optptsInput failed");
     return;
  }
  MsgOut(status, "sc2da_SQ1optpts:  _sq1optptInit OK");
 
  sc2dalib_sq1optptsreadInput(dasInfo.fpBatch, sq1bias,sq2feedback,sq1refbias,sq2reffeedback,status);
  if (!StatusOkP(status))
  {
     ErsRep(0,status,"sc2da_SQ1optpts: _sq1optptsreadInput failed");
     return;
  }
  if(dasInfo.fpBatch)
    fclose (dasInfo.fpBatch);
  dasInfo.fpBatch=NULL;
  MsgOut(status,"sc2da_SQ1optpts: _sq1optptsreadInput OK");

  sc2dalib_sq1optwritesq1bsq2fOut(&dasInfo,&setup, sq1bias, sq2feedback, sq1refbias, sq2reffeedback,1,status);
  if (!StatusOkP(status))
  {
    ErsRep(0,status,"sc2da_SQ1optpts: _sq1optwritesq1bsq2fOut failed");
    return;
  }
  
  // flag=1 takeout outlier of sq2fb set outlier=BAD;  flag=0 don't take out
  // flag =3 use simple clipmean for row and column
  if ( flag==3)
    sc2dalib_sq1biassq2fdbkClipmean(sq1bias, sq1bcomp, sq2feedback, sq2fcomp,COL_NUM,ROW_NUM, status);
  else
    sc2dalib_sq1optAlgrm(&dasInfo,&setup, sq1bcomp, sq2fcomp,sq1bias, sq2feedback, sq1refbias, 
                       sq2reffeedback,bqual,fqual,squal,sq1scale,flag,status);

  if (!StatusOkP(status))
  {
    ErsRep(0,status,"sc2da_SQ1optpts:_sq1biassq2fdbkClipmean or _sq1optAlgrm failed");
    return;
  }

  sc2dalib_sq1optsave2File(&dasInfo,&setup,sq1bcomp, sq2fcomp,status);
  if (!StatusOkP(status))
  {
    ErsRep(0,status,"sc2da_SQ1optpts: _sq1optsave2File failed");
    return;
  }

}



/**
 * \fn void sc2da_Trksq2fb(StatusType *status)
 *
 * \brief daram action
 *  track sq2fb, read a single frame data back and send updated sq2fb 
 *  values
 *
 * \param *status StatusType.  given and return
 *
 *
 * It is needed to adjust the sq2fb to compensate earth's magnetic pick up. 
 * this can be a single update (OPTION=0) or continuous update (OPTION!=0)
 * until told to stop. If it is continuous,
 *
 * It sends RTS Client an acknowledgement before rescheduling to wait for
 * a single data frame from the MCE.  After receiving a data frame, it 
 * calculates the required change to keep ssa lock point, sends 
 * updated values to sq2 feed back and the DRAMA is informed that the action
 * expects to reschedule.  It continues doing this until is told to end by 
 * the RTS Client. The data frame from the MCE is saved in 
 * $ORAC_DATA_OUT/sq2fbtrk.txt 
 * 
 * ditscmd SC2DA SQ2FB_TRACK \
 *               OPTION=(0 or 1) \
 *               DATA_FILE=xxx \
 *               SETU_FILE=xxx 
 *               SQ2OPT_FILE=xx
 *               SQ1OPT_FILE=xx
 *
 * SQ2FB_TRACK set sq2fb=initFB , sq1bias=0  at beginning
 *
 */
/*+ sc2da_Trksq2fb  
*/
void sc2da_Trksq2fb
( 
StatusType *status    
)
{
  int             rval,i,wait;
  uint32          mceBufsize;
  char            tmpPara[FILE_LEN]="";
  char            *byte;
  static FILE     *fp;
  static int      skipDo;
  static char     dateTime[40];
  static ARRAYSET setup;
  static long     option;
  static dasCmdInfo_t trksq2fbCmd;
  static char    renametrksq2fbFile[]="rename-trksq2fbfile";      
  static char    sq2fbCmd[]="wb bc2 flux_fb";

  if (!StatusOkP(status)) return;

  if(DitsGetSeq()==0)
  {
    rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);

    /* This routine will NOT work with external clocks because there is
       no DV from the RTS when it is being used so set the clocking to
       all internal here */

    sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,syncCmdInt,status); 
    sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,dvCmdInt,status); 
    // sc2dalib_setmceVal(con,&dasInfo,&glbMceInx,selClkCmdInt,status); 
    if ( !StatusOkP(status))   
      {
        ErsRep(0,status,"sc2da_TrkSq2Fb: failed after calls to sc2dalib_setmceVal");
        return;
      }

    //dasInfo.trkNo set to 0/ sq2fb set to setup->initFB
    sc2dalib_trksq2fbInit(con,&dasInfo,&trksq2fbCmd,&glbMceInx,dateTime,&setup,
                          &option,setup.sq2fdbkOpt,setup.sq1biasOpt,status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2da_Trksq2fb: sc2dalib_Trksq2fbInit failed");
      sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
      return;
    }
    // only for test   
    if (dasInfo.filesvFlag !=0)
    {
      MsgOut(status,"we are going to save frame data %s, can be very big!!",
              dasInfo.dataFile);
    }
    else
      MsgOut(status,"we are not going to save frame data ");

    for(i=0;i<COL_NUM;i++)
    {
      if ( (setup.colMask[i] !=0) && (setup.initFB[i] !=0) )
      {
        MsgOut(status,"sq1gainScale[%d]=%6.4E",i,setup.sq1gainScale[i]);
      }
    }
    
    // dasInfo.msgwrtPt=0; dasInfo.msgreadPt=0;
    // after call sc2dalib_frametakeInit, it set trkNo=0;
    sc2dalib_frametakeInit(con,&dasInfo,1,1,&glbMceInx,dateTime,status);
    if (!StatusOkP(status)) 
    {
      ErsRep(0,status, "sc2da_Trksq2fb: sc2dalib_frametakeInit failed");
      sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
      return;
    }
   
    // doServo: howmany frame/updatesq2FB 
    if (option ==0)
       setup.doServo=1;
    
    skipDo=setup.doServo;

    // open another file for appending, mostly used for taking sq2fb=sq2fbOPT
    if (setup.slopSelect[5]==0)
    {
      fp=NULL;  
      sprintf(tmpPara,"%s-sq2fbopt",dasInfo.dataFile);
      if((fp = fopen64(tmpPara, "a")) == NULL )
	{
	  *status = DITS__APP_ERROR;
	  ErsRep (0, status, "Error- sc2da.c sc2da_Trksq2fb() failed to open file %s", tmpPara); 
	  return;
	}
    }
    
    // for data frame now
    dasInfo.actIndex=DitsGetActIndex();
    con->process.framesetup=DA_MCE_SINGLEDATA;  
    DitsPutRequest ( DITS_REQ_SLEEP, status );
    post_sem(con->intask_sem,0);
  }
  else
  {
    if(glbMsg[dasInfo.msgreadPt].reason ==FRAME_COMPLETION)   
    {
      //data frame finishes, no more expected 
      byte=(char*)glbMsg[dasInfo.msgreadPt].data;
      mceBufsize=glbMsg[dasInfo.msgreadPt].bufsize;
      // recycle=16
      dasInfo.msgreadPt++;
      dasInfo.msgreadPt=dasInfo.msgreadPt & 0x000F ; 

      // save the raw frame data
      if (dasInfo.filesvFlag !=0)
      {
        if (setup.slopSelect[5]!=0)
          sc2dalib_saveframeData(con,&dasInfo,dasInfo.fpData,byte,mceBufsize/4,status);
        else
        {
          if ( setup.fbFlag ==1 )
            sc2dalib_saveframeData(con,&dasInfo,fp,byte,mceBufsize/4,status);
          else if ( setup.fbFlag ==0 )
            sc2dalib_saveframeData(con,&dasInfo,dasInfo.fpData,byte,mceBufsize/4,status);
        }
      }
      // the first 15 data frame don't aplly servo algorithm
      dasInfo.trkNo ++;          
      if (option ==1 && dasInfo.trkNo ==1)     
      { 
         /*  here comes acknowledge to RTSC */
        DitsTrigger(0,status);
        if (!StatusOkP(status)) 
        {
          ErsRep(0,status, "sc2da_Trksq2fb: DitsTrigger failed");
          sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
          return;
         }
         MsgOut(status," Start sq2fb tracking ............");
      }
      jitDebug(16,"sc2da_Trksq2fb: call trksq2FB skipDo=%d\n", skipDo);
 
      if (setup.slopSelect[5]!=0)
      {
        skipDo --;
        sc2dalib_trksq2FB(con,&dasInfo,&glbMceInx,dateTime,byte,&setup,
                        setup.sq2fdbkOpt,option,skipDo,status);
      }
      else
      {
         if ( setup.fbFlag ==0 )
         {
           // we calculate sq2fb delta, but send sq2fbOPT
           // the first 15 data frame don't aplly servo algorithm
           skipDo --;
           sc2dalib_trksq2FB(con,&dasInfo,&glbMceInx,dateTime,byte,&setup,
                        setup.sq2fdbkOpt,option,skipDo,status);
           setup.fbFlag =1;
          jitDebug(4,"send sq2fbOPT \n");
         }
         else
         {
           // in case of kick, don't set sq2fb=initFB
	   if (dasInfo.actionFlag !=TRKKICK ) 
	   {
             // we set initFB for next data
             sc2dalib_sendarrayCmd(con,&dasInfo,&glbMceInx,dateTime,sq2fbCmd,
                                  setup.initFB,COL_NUM,status);
             setup.fbFlag=0;
             jitDebug(4,"send sq2FB \n");
	   }
         }
      }    
      if (!StatusOkP(status)) 
      {
        ErsRep(0,status, "sc2da_Trksq2fb: _trksq2FB or _sendarrayCmd failed");
        sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
        if (setup.slopSelect[5]==0)
        {
          if(fp)
            fclose(fp);
          fp=NULL;
        }
        return;
      }

      if ( skipDo ==0) 
        skipDo=setup.doServo;
     
      if (option ==0)
      {
        sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
        MsgOut(status,"sc2da_Trksq2fb: single update completed");
        return;
      }

      if(  dasInfo.actionFlag !=TRKKICK )
      {
        // allow  fraction of second if  slopSelect[7]>0
        if( setup.slopSelect[7] > 0 )
        {
          wait=(1000000*setup.slopSelect[6])/setup.slopSelect[7];
          usleep(wait);
        }
        *status=sdsu_command_mce(con,trksq2fbCmd.cmdBuf,&trksq2fbCmd.reply);
        if ( *status !=SDSU_OK)
        {
          ErsRep(0,status,"sc2da_Trksq2fb: sdsu_command_mce failed");
          sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
          if (setup.slopSelect[5]==0)
          {
             if(fp)
               fclose(fp);
             fp=NULL;
          }
          return;
        }
        else
        {
          DitsPutRequest ( DITS_REQ_SLEEP, status );
          post_sem(con->intask_sem,0); 
        }
      }
      else
      {
        system (renametrksq2fbFile);
        sc2dalib_trksq2fbsavefbZG(&dasInfo,&setup,status);
        sc2dalib_trksq2fbsavesq2fbOPT(&dasInfo,setup.sq2fdbkOpt,status);

        sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
        if (setup.slopSelect[5]==0)
        {
          if(fp)
            fclose(fp);
          fp=NULL;
        }
        MsgOut(status,"sc2da_Trksq2fb: Kick msg received");
        jitDebug(16,"sc2da_Trksq2fb: the action has completed");
      }
    }
    else //
    {
      con->process.seqstatus=SEQ_ERROR; 
      *status=DITS__APP_ERROR;
      sc2dalib_msgprintSave(&dasInfo,
               "sc2da_Trksq2fb: Error: %s",glbMsg[dasInfo.msgreadPt].errRep,USE_ERSREP,status);
      sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
      if (setup.slopSelect[5]==0)
      {
        if(fp)
          fclose(fp);
        fp=NULL;
      }
    }
  }
}


/**
 * \fn void sc2da_Trkheat(StatusType *status)
 *
 * \brief drama action
 *  track heaters, read a single frame data back and send updated values
 *  to all heaters
 *
 * \param *status StatusType.  given and return
 *
 *
 * It is needed to adjust the heater to keep the TES array in a working 
 * condition when the shutter is opened or closed. 
 *
 * It sends RTS Client an acknowledgement before rescheduling to wait for
 * a single data frame from the MCE.  After receiving a data frame, it 
 * compares the data with some lookup table, does some calculations, sends 
 * updated values to all heaters and the DRAMA is informed that the action
 * expects to reschedule.  It continues doing this until is told to end by 
 * the RTS Client. The data frame from the MCE is saved in rtsltFile if 
 * filesvFlag=1.
 * 
 * ditscmd SC2DA HEATER_TRACK \
 *               SVFILE_FLAG=(0 or 1) \
 *               DATA_FILE=xxx \
 *               SETU_FILE=xxx 
 *
 */
/*+ sc2da_Trkheat - 
*/
void sc2da_Trkheat
( 
StatusType *status    
)
{
  int           rval;
  int           *word32;
  int           keepTracking;
  char          *glbmsgPtr,*lclmsgPtr; 
  char          heatVal[]="rb bc1 bias 1";
  long          heatFlag;
  char          shutterState[DITS_C_NAMELEN];
  char          heatCmd[20];
  static DRAMA_INNERMSG    dramamsg;
  static uint32 mceBufsize;
  static dasCmdInfo_t trkheatCmd;
  static char   *byte;
  static char    dateTime[40];
  static char    renametrkheatFile[]="rename-trkheatfile";      
  
  if (!StatusOkP(status))
  {
    MsgOut(status, "sc2da_Trkheat: called with bad status");
    return;
  }

  if(DitsGetSeq()==0)
  {
    /* DEBUG for abort testing */
    MsgOut(status, "sc2da_Trkheat: entered DitsGetSeq()=0");

    rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);

    /* This routine will NOT work with external clocks because there is no DV
       from the RTS when it is being used so set all of the clocking to internal */

    sc2dalib_setmceVal(con,&dasInfo, &glbMceInx, syncCmdInt, status); 
    sc2dalib_setmceVal(con,&dasInfo, &glbMceInx, dvCmdInt, status); 
    // sc2dalib_setmceVal(con,&dasInfo, &glbMceInx, selClkCmdInt, status); 
    if ( !StatusOkP(status))   
      {
        ErsRep(0,status,"sc2da_Trkheat: failed after calls to sc2dalib_setmceVal");
        return;
      }

    /* DEBUG for abort testing */
    MsgOut(status, "sc2da_Trkheat: starting sc2dalib_trkheatInit()");

    //Initialise for heater tracking

    sc2dalib_trkheatInit(con,&dasInfo,&trkheatCmd,&glbMceInx,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2da_Trkheat: sc2dalib_trkheatInit failed");
      sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
      return;
    }  

    // dasInfo.msgwrtPt=0; dasInfo.msgreadPt=0;
    // after call sc2dalib_frametakeInit
    sc2dalib_frametakeInit(con,&dasInfo,1,1,&glbMceInx,dateTime,status);
    if (!StatusOkP(status)) 
    {
      ErsRep(0,status, "sc2da_Trkheat: sc2dalib_frametakeInit failed");
      sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
      return;
    }
    con->process.framesetup=DA_MCE_SINGLEDATA;  
    // for data frame now
    dasInfo.actIndex=DitsGetActIndex();
    DitsPutRequest ( DITS_REQ_SLEEP, status );
    post_sem(con->intask_sem,0);
  }
  else // The entry to this action was not caused by an OBEY
  {
    // get the glbmsg, we can use it directly later
    glbmsgPtr=(char*)&glbMsg[dasInfo.msgreadPt];
    lclmsgPtr=(char*)&dramamsg;
    memcpy(lclmsgPtr,glbmsgPtr,sizeof(DRAMA_INNERMSG));       
    // recycle=16
    dasInfo.msgreadPt++;
    dasInfo.msgreadPt=dasInfo.msgreadPt & 0x000F ; 

    if(dramamsg.reason==FRAME_COMPLETION)   
    {
      //data frame finishes, no more expected 
      byte=(char*)dramamsg.data;
      word32 = (int *)byte;
      mceBufsize=dramamsg.bufsize;
      //jitDebug(2,"sc2da_Trkheat: PacketSize from MCE=<%ld>",mceBufsize); 
      if (dasInfo.filesvFlag)
      {
        sc2dalib_saveframeData(con,&dasInfo,dasInfo.fpData,byte,mceBufsize/4,status);
      }
      dasInfo.trkNo++;   
      if ( dasInfo.trkNo == 1) // Is this the second time I have entered this action, if so trigger the RTSC     
      { 
         /*  here comes acknowledge to RTSC */
        DitsTrigger(0,status);
        if (!StatusOkP(status)) 
        {
          ErsRep(0,status, "sc2da_Trkheat: DitsTrigger failed");
          sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
          return;
         }
         MsgOut(status," Start tracking ............");
      }

      /* If we have not been kicked, or reached the requested number of times through, then heater track */
      /* If dasInfo.numTimesThru is zero then we are requested to heater track until we are kicked */
      keepTracking = 1;
      if(dasInfo.numTimesThru > 0)
	{
	  if(dasInfo.trkNo > dasInfo.numTimesThru)
	    {
	      keepTracking = 0;
	    }
	}
 
      if( (dasInfo.actionFlag != TRKKICK) && keepTracking)  
      { 
        // reschedule to read again
        /*  here comes the heater lookup table
            reads the table, determines the new heat setting, 
            sends the new value to heater DAC
        */
        // it will always be a full size of the data
        sc2dalib_trkheatUpdate(con,&dasInfo,&glbMceInx,word32,heaterSlope,status);
        if (!StatusOkP(status)) 
        {
          ErsRep(0,status, "sc2da_Trkheat: sc2dalib_trkheatUpdate failed");
          sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
          return;
        }
         usleep(10000);
        /// only go is needed here, can use sc2dalib_Cmd, which saves all cmd 
        // and reply but takes time and will result in a huge logfile
        *status=sdsu_command_mce(con,trkheatCmd.cmdBuf,&trkheatCmd.reply);
        if ( *status != SDSU_OK)
        {
          ErsRep(0,status,"sc2da_Trkheat: sdsu_command_mce failed");
          sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
          return;
        }
        else
        {
          DitsPutRequest ( DITS_REQ_SLEEP, status );
          post_sem(con->intask_sem,0); 
        }    
      }
      else // dasInfo.actionFlag == TRKKICK (we got kicked)  
      {
        if(dasInfo.debuglvl==2)    
        {
          sc2da_Dispinfo(status);
        }

	/* If the WALK_HEATER_CURRENT_BIT is set in the HTTRACK_FLAG and we
           have been told that the Shutter is closed by the SCUBA2 RTS client
           then walk the heater current back to its default value */

	SdpGeti("HTTRACK_FLAG", &heatFlag, status);
	SdpGetString("SHUTTER_STATUS", DITS_C_NAMELEN, shutterState, status);
	if(((heatFlag & WALK_HEATER_CURRENT_BIT) != 0) &&
	   (strncmp(shutterState,"CLOSED",6) == 0))
	  {
	    sprintf(heatCmd, "wb bc1 bias %d", dasInfo.darkHeaterI);
	    sc2dalib_setmceVal(con, &dasInfo, &glbMceInx, heatCmd, status);
            /* fprintf(dasInfo.fpLog, "_Trkheat %s",heatCmd); */ 
	    if ( !StatusOkP(status) )
	      {
		ErsRep(0,status,"sc2dalib_Trkheat: sc2dalib_setmceval(1) %s failed",heatCmd); 
		return;
	      }
	  }

        /* Execute shell script to make a timestamped copy of trkheat.txt */
        system (renametrkheatFile);

	// read the current heater setting and store as the nominal

	sc2dalib_readmceVal(con, &dasInfo, &glbMceInx, heatVal, &dasInfo.nominalPixelHeat,1,status);
	if (!StatusOkP(status)) 
	  {
	    ErsRep(0,status, "sc2da_Trkheat: sc2dalib_readmceVal failed to read heater value"); 
	    return;
	  }
	MsgOut(status, "sc2da_Trkheat: The nominal heater current is: %d",dasInfo.nominalPixelHeat);

	/* Set all the clocks to external again */

	sc2dalib_setmceVal(con,&dasInfo, &glbMceInx, syncCmdEx, status); 
	sc2dalib_setmceVal(con,&dasInfo, &glbMceInx, dvCmdEx, status); 
	// sc2dalib_setmceVal(con,&dasInfo, &glbMceInx, selClkCmdEx, status); 
	if ( !StatusOkP(status))   
	  {
	    ErsRep(0,status,"sc2da_Trkheat: failed after calls to sc2dalib_setmceVal (external clocks)");
	    return;
	  }

	/* Tidy up files for this action */
        sc2dalib_actionfileEnd(con,&dasInfo,0,status); 
        jitDebug(4,"sc2da_Trkheat: the action has completed");

        /* Execute shell script to make a timestamped copy of trkheat.txt */
        system (renametrkheatFile);

      } /* End of things action does when kicked */
    }
    else 
    {

      MsgOut(status,"sc2da_Trkheat: The dramamsg.reason is %d\n",(int)dramamsg.reason);
      con->process.seqstatus=SEQ_ERROR; 
      *status=DITS__APP_ERROR;
      sc2dalib_msgprintSave(&dasInfo,
			    "sc2da_Trkheat: Error: %s",dramamsg.errRep,USE_ERSREP,status);
      sc2dalib_actionfileEnd(con,&dasInfo,0,status);

     }
  }
}


/**
 * \fn void sc2da_trkKick(StatusType *status)
 *
 * \brief drama action
 *  kick sc2da_Trkheat,set dasInfo.actionFlag=TRKKICK  
 *  and use DitsSignalByIndex to end sc2da_Trkheat
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * >ditscmd SC2DA -k MCETRKHEATER 
 */
/*+ sc2da_trkKick - 
*/
void sc2da_trkKick
( 
StatusType *status    
)
{
  char    dateTime[40];
  int     rval;
  int     msgreadPt;
  SdsIdType arg;
  char shutterStatus[DITS_C_NAMELEN];
  int valIndex;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  fprintf(dasInfo.fpLog,"\n<%s> CMD from sc2da_trkKick \n",
          dateTime);
  dasInfo.actionFlag=TRKKICK; 

  /* read the shutter status from the SCUBA2 RTS Client */
  arg = DitsGetArgument();
  jitArgGetS(arg, "SHUTTER_STATUS", 0, NULL, "UNKNOWN", GIT_M_ARG_UPPER, DITS_C_NAMELEN,
	     shutterStatus, &valIndex, status);
  SdpPutString("SHUTTER_STATUS", shutterStatus, status);
  /* MsgOut(status, "sc2da_trkKick: shutter_status: %s",shutterStatus); */

  // get the glbmsg, we need to use msgreadPt-- 
  if (dasInfo.msgreadPt !=0)
    msgreadPt=dasInfo.msgreadPt-1;
  else
    msgreadPt=15;


  if (glbMsg[msgreadPt].reason != FRAME_COMPLETION)
  {
    if ( dasInfo.actIndex != -1 )
    {
      DitsSignalByIndex ( dasInfo.actIndex, 0, status );
    }
    else
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,"sc2da_trkKick: not valid dasInfo.actIndex");
    }  
  }
  jitDebug(4,"sc2da_trkKick: Kick has completed");
}


/**
 * \fn void sc2da_UpdateLogname(StatusType *status)
 *
 * \brief drama action
 *  update the logfile name if the task is left for long time
 *  or over night.
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * update the dasInfo.Date and logfileName, reset MCE_COUNT_COUNT=0
 *
 * >ditscmd SC2DA UPDATELOGNAME
 */

/*+ sc2da_UpdateLogname - 
*/
void sc2da_UpdateLogname
(
StatusType *status
)
{ 
  int  rval;
  char dateTime[40];
 
  if (!StatusOkP(status)) return;

  if ( con->process.framesetup !=DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status, 
      "sc2da_UpdateLogname: %s has not completed",
       seqStatus[dasInfo.actionFlag]);
    return;
  } 
  rval=sc2dalib_finddateTime(DAS_DATE,dateTime);
  if( strcmp(dateTime,dasInfo.Date)!=0 )
  {
    strcpy(dasInfo.logfileName,dasInfo.logFile);
    strcat(dasInfo.logfileName,dasInfo.Date);
    my_fclose(&(dasInfo.fpLog));
    if ((dasInfo.fpLog = fopen64(dasInfo.logfileName,"a")) == NULL)
    {
      *status = DITS__APP_ERROR;
      ErsRep(0,status, "sc2da_UpdateLogname: failed to open %s",
	     dasInfo.logfileName );
      return;
    }

    // myFpLog is a global so all routines can write in the log file
    myFpLog = dasInfo.fpLog;

    dasInfo.glbCount=0;
    my_fclose(&(dasInfo.fpLog));
  }
  jitDebug(4,"mceFindDate: the action has completed, the logfileName=%s\n",
            dasInfo.logfileName );
}


/**
 * \fn void sc2da_Version(StatusType *status)
 *
 * \brief drama action
 *  get driver, das and pci version number.
 *
 * \param *status StatusType.  given and return
 *
 * >ditscmd SC2DA DASVERSION
 */

/*+ sc2da_Version - 
*/
void sc2da_Version
(
StatusType *status
)
{
  int         rval;
  static char   dateTime[40];
  static  DRV_USER_STRUCT userdrvMsg;
  
  if (!StatusOkP(status)) return;

  // Check whether this is the start or completion  
  if(DitsGetSeq()==0)
  {
    rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);   
    // update debug, check if it is in DA_MCE_SEQ setup actionFlag
    sc2dalib_versionInit(con,&dasInfo,dateTime,status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2da_Version: sc2dalib_versionInit failed");
      return;
    }
    userdrvMsg.action=GET_VERSION;
    strcpy(userdrvMsg.msg,taskName);
    strcpy(userdrvMsg.error,"request driver version");
    //ask driver to send back
    rval=master_write_msg2drv(con, userdrvMsg);
    if( rval !=SDSU_OK )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status,"sc2da_Version: master_write_msg2drv failed");
      return;
    }
    dasInfo.actIndex=DitsGetActIndex();
    DitsPutRequest ( DITS_REQ_SLEEP, status );
  } 
  else
  {
    //the messag from driver is back now, but stored in global driverVersion 
    sc2dalib_versionInfo(con,&dasInfo,status);
  }
}



/**
 * \fn void mceTryChkSum(StatusType *status)
 *
 * \brief drama action
 *  it is similar to sc2da_Mcecmd but the checksum is replaced by 
 *  the user input value not calculated one 
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * >ditscmd SC2DA MCECHKSUM  "WB cc sram1_strt 10" checksum
 */
/*+ mceTryChkSum - 
*/
void mceTryChkSum
(
StatusType *status
)
{ 
  int         rval;  	        
  char        *insertPtr1;
  long      checksum, *insertPtr;
  DitsArgType argId;
  dasCmdInfo_t  chksumCmd;
  char         dateTime[40];

  static struct mcexml_struct inxtmp;
     
  if (!StatusOkP(status)) return;

  if (dasInfo.doneReadxml==0 )
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status, "TryChkSum: the MCE.xml has not been read" ); 
    return;
  }
  if ( con->process.framesetup !=DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status, 
      "TryChkSum: %s has not completed",seqStatus[dasInfo.actionFlag]);
    return;
  } 

  dasInfo.actionFlag=CHKSUMACTION;
  inxtmp=glbMceInx;
  dasInfo.cmdFlag=0;
  argId = DitsGetArgument();
  ArgGetString ( argId, "Argument1", CMD_LEN, chksumCmd.mceCmd, status );
  ArgGeti      ( argId, "Argument2", &checksum, status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"TryChkSum: failed to get args");
    return;
  }
  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  dasInfo.filesvFlag=0;//no save to fpMcecmd
  fprintf(dasInfo.fpLog,"\n<%s> CMD from mceTryChkSum\n",dateTime);

  sc2dalib_getcmdBuf(&dasInfo,&chksumCmd,&glbMceInx,status);
  if ( *status !=STATUS__OK)
  {     
    ErsRep(0,status,"TryChkSum: Failed after call to sc2dalib_getcmdBuf");  
    return;
  }
  insertPtr1=(char*)&chksumCmd.cmdBuf+252;    // for checksum
  insertPtr=(long*)insertPtr1;    
  *insertPtr=checksum;
  con->testflag=1;
  
  sc2dalib_Cmd(con,&dasInfo,&chksumCmd,dateTime,status);
  if ( *status !=STATUS__OK)
  {     
    ErsRep(0,status,"TryChkSum: Failed after call to sc2dalib_Cmd");  
    return;
  }
  // all errors message have been displayed in mceErrRep
  sc2dalib_dispResults(&dasInfo,&chksumCmd,status);
  sc2dalib_chkChecksum(&dasInfo,&chksumCmd,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"TryChkSum: the action completed with error"); 
    return;
  }
  MsgOut(status,"TryChkSum: the action has completed");
}



/**
 * \fn void mceTest(StatusType *status)
 *
 * \brief drama action
 *  do some repeat tests to check time 
 *
 * \param *status StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * ditscmd SC2DA TEST \
 *         "WB cc sram1_strt 10" \
 *         number_Of_repeat_times \
 *         inctranslateCall
 */
/*+ mceTest - 
*/
void mceTest
(
StatusType *status
)
{  
  int         rval;
  long      repeatNo, i,inctranslateCall=0;
  double      seq_time;
  char         dateTime[40];
  DitsArgType argId;
  static time_t   tloc, newtloc;        
  static  dasCmdInfo_t  testCmd;
   
  if (!StatusOkP(status)) return;

  if (dasInfo.doneReadxml==0 )
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"mceTest: the MCE.xml has not been read" );
    return;
  }
  if ( con->process.framesetup !=DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status, 
      "mceTest: %s has not completed",seqStatus[dasInfo.actionFlag]);
    return;
  } 
  dasInfo.actionFlag=TESTACTION;
 
  argId = DitsGetArgument();
  ArgGetString ( argId, "Argument1", CMD_LEN, testCmd.mceCmd, status );
  ArgGeti      ( argId, "Argument2", &repeatNo, status );
  ArgGeti      ( argId, "Argument3", &inctranslateCall, status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"mceTest: failed to get args"); 
    return;
  }
  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  dasInfo.filesvFlag=0;//no save to fpMcecmd
  fprintf(dasInfo.fpLog,"\n<%s> CMD from mceTest\n",dateTime);
  sc2dalib_getcmdBuf(&dasInfo,&testCmd,&glbMceInx,status);
  if ( *status !=STATUS__OK)
  {     
    ErsRep(0,status,"mceTest: Failed after call to sc2dalib_getcmdBuf");  
    return;
  }
 
  //get time in seconds since EPOCH 
  time(&tloc);

  for (i=1;i<=repeatNo;i++)
  {
    if(inctranslateCall)
    {
      sc2dalib_getcmdBuf(&dasInfo,&testCmd,&glbMceInx,status);
      if ( *status !=STATUS__OK)
      {     
         ErsRep(0,status,"mceTest: Failed after call to sc2dalib_getcmdBuf");  
         return;
      }
    }
    *status=sdsu_command_mce(con,testCmd.cmdBuf,&testCmd.reply);
    if( *status !=SDSU_OK)
    { 
      sc2dalib_mceerrRep(&dasInfo,&testCmd,status); 
      dasInfo.cmdFlag=1;
      *status=DITS__APP_ERROR;
      ErsRep (0, status,"mceTest: the action ended with error"); 
      return;
    }
    sc2dalib_chkChecksum(&dasInfo,&testCmd,status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"mceTest: the action completed with error"); 
      return;
    }
  }
  // get time in seconds since EPOCH 
  time(&newtloc);
  seq_time=difftime(newtloc,tloc); 
  MsgOut(status,
   "mceTest: The total_time of (%ld cmds) ~= <%f>s, rate= <%5.1f>us/cmd",
    repeatNo,seq_time,(seq_time*1000000.0)/((double)repeatNo) );   
  if(inctranslateCall)
  {
    MsgOut(status,
      "mceTest: The loop includes{mceXMLtranslate; send/get cmd/reply; sc2dalib_chkChecksum} " );    
  }
  else
  {
    MsgOut(status, 
      "mceTest: The loop includes {send/get cmd/reply; sc2dalib_chkChecksum} " );
  }
  MsgOut(status,"mceTest: the action has completed");
}


/**
 * \fn void sc2da_Getsq2fbPara(StatusType *status)
 *
 * \brief daram action
 *  get args to ARRAYSET *setup0  for sq2fb, ready for _changesq2fbval to use in SEQ
 *
 * \param *status StatusType.  given and return
 *
 * ditscmd SC2DA SQ2FB_GETPARA \
 *               SETUP_FILE=xxx 
 *
 * SQ2FB_GETPARA set sq2fb=initFB , sq1bias=0  at beginning
 *
 */
/*+ sc2da_Getsq2fbPara  
*/
void sc2da_Getsq2fbPara
( 
StatusType *status    
)
{
  int             rval,i;
  static char     dateTime[40];
  static dasCmdInfo_t getsq2fbParaCmd;
 
  if (!StatusOkP(status)) return;

  rval=sc2dalib_finddateTime(DAS_DATETIME,dateTime);
  //  sq2fb set to setup->initFB
  sc2dalib_getsq2fbparaInit(con,&dasInfo,&getsq2fbParaCmd,&glbMceInx,dateTime,&setup0,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2da_getsq2fbPara: _getsq2fbparaInit failed");
    sc2dalib_actionfileEnd(con,&dasInfo,1,status); 
    return;
  }
  for(i=0;i<COL_NUM;i++)
  {
    if ( (setup0.colMask[i] !=0) && (setup0.initFB[i] !=0) )
      MsgOut(status,"sq1gainScale[%d]=%6.4E",i,setup0.sq1gainScale[i]);
  }
}

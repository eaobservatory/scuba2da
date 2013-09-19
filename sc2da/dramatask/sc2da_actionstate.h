/**
 * \file sc2da_actionstate.h
 *
 * \brief action state string for SCUBA2 Data Acquisition Software (DAS)  
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
 *  $Log: sc2da_actionstate.h,v $
 *  Revision 1.14  2010/02/19 00:29:10  cwalther
 *  Changes to make fast flat fields work
 *
 *  Revision 1.13  2009/12/15 23:31:17  cwalther
 *  Added more detailed status for malloc failures
 *
 *  Revision 1.12  2009/11/24 22:02:03  cwalther
 *  Code for displaying the counters for each STATE data event comming in
 *
 *  Revision 1.11  2009/08/11 15:10:27  xgao
 *  update to UPDATE_SQ2FB SQ1LOCK and tidyup scripts
 *
 *  Revision 1.10  2008/12/12 22:27:19  xgao
 *  sq2fb servo test
 *
 *  Revision 1.9  2008/08/14 02:56:08  cwalther
 *  These changes from Gao snuck in but I am going to commit them as they seem to do no harm
 *
 *  Revision 1.8  2008/08/07 01:49:34  cwalther
 *  Lots for changes for headers from Craig. Changes in handling events from Gao
 *
 *  Revision 1.7  2008/06/28 20:11:54  xgao
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
 *  Revision 1.6  2007/12/07 16:01:19  xgao
 *
 *  some modifications after SG array test, shutter only third open, new
 *  initialise.xml ( flatfile, INITFLAG) and configure.xml (DATA_MODE).
 *  negative pixHeat for loadv in SETU_SEQ. Add findheaterspike for processing
 *  heaterservo data for G/Tc.
 *
 *  Revision 1.5  2007/11/02 12:14:08  xgao
 *
 *  update to 0p33-021107, worked with JOC for INIT,CONFIG, SETUPSEQ, SEQ
 *
 *  Revision 1.4  2007/09/05 16:33:48  xgao
 *
 *  add heaterslope, trkHeater tested Ok, sq1Refbias and sq2Reffb always use the last one
 *  always do manual adjustment after optimal, fits4all is now in parShm, bigPhy=512M now
 *  0p21, ndf is turned on, tested with startmce-rtsc flatfield-seg
 *
 *  Revision 1.3  2007/06/18 23:22:32  xgao
 *  add functions for SCRIPTDAS
 *
 *  Revision 1.2  2007/05/23 18:09:02  xgao
 *  use JCMTState, add a few actionstate. xg
 *
 *  Revision 1.1.1.1  2007/05/16 08:27:01  dkelly
 *  first insertion
 *
 */

/* Please Note! the count in this array of strings must
   be the same as the count in the enum called DA_PROCESS
   found in sc2da_par.h */

static char *seqStatus[]=
{
  /// frameSetup
  // from 0
"DA not in any data taking mode",
"DA in sequence frame data taking mode",   
"DA in single data taking mode",
"DA in tracking heater data taking mode",   
"DA in donloadFPGA data taking mode",   
"DA in AUTOSET ARRARY data taking mode",
"DA in DREAM data taking mode",   
"DA in SSASET ARRARY data taking mode",   
"DA in SQ2SET ARRARY data taking mode",   
"DA in SQ1SET ARRARY data taking mode",

 //SEQUENCE
  // from 10
"no sequence action",
"still in sequence action",
"sequence action finished OK",  
"sequence action finished by kick", 
"sequence action finished with Error",
 
  //actionFlag
  // from 15
"none action ever called",
"abort action called",
"batch action called",
"bac write/read action called",
"configuare action called",
"test checksum action called",
"check version action called",
"debug action called",
"dispInfo action called",
"downloadDSP action called",
"exit action called",
"help action called",
"downloadFPGA action called",
"initial action called",
"mcecmd action called",
"mcestatus action called",
"pcicmd action called",
"ping  action called",
"sequence action called",
"setup sequence action called",
"sequence kick action called",
"servo action called",
"servo kick action called",
"tets action called",
"track heat action called",
"track kick action called",
"updatelogname action called",
"in action for FDH_START",
"stoped by  kick ", 
"monitor pixel action called",     
"monitor pixel kick action called",
"scriptdas action called",
"scriptdas kick action called",
"heaterslope action called",
"set data format and save-engData action called",
"track sq2fb action called",
"find sq1 optimal points action called",
"get sq2fb para for sq2fb tracking in SEQ",

  //batchFlag
  //from 53
"GO from batch file",
"nomal cmd from batch file",
"comment from batch file",
"no more cmd in batch file",
"send cmd from batch file failed",
"translate cmd from batch file failed",


  /// frame finishing reason
  // from 59
"wait for frame timeout", 
"idle or wait for first frame", 
"total frames have finished",
"frame finished with error in sdsu_getpass buf",    
"first frame",
"frame finished with less expected", 
"frame finished by stop/kick", 
"frame finished with no lastframe bit set", 
"frame checksum wrong",    
"frame is ready for headers and QL if not DARK",   
"open data file failed, check permission setting",
"malloc for holding frame data failed",
"fail to connect shared Memory",
"scan is ready for QL if not DARK",
"sq2fb compansation",
"malloc for  holding windowed data failed",
"malloc for  holding temp dream failed",
"malloc for  holding re-organised frame failed",
"malloc for  holding ch data failed",
"malloc for  holding frame status failed",
"malloc for  holding dark row failed",
"Update heater current during fast flat field",
 
 /// OBS mode
// from 81
"STARE mode",
"DREAM mode",
"SCAN mode",

  /// Msg FIFO error status
  // from 84
"idle or non error in put msg into msg FIFO",
"writing USER_USER_MSG FIFO error",
"writing USER_USER_FRAME FIFO error",    
"writing USER_USER_ERR FIFO error",    
 
  /// operation on FIFOs
  // from 88
"idle or non error in writing/opening FIFO",
"error in opening USER_USER_EM USER_DATA_EM FIFO",
"error in opening USER_USER_MSG FIFO",
"error in opening USER_USER_FRAME FIFO",
"error in opening USER_USER_ERR FIFO",

  ///whereabout the process
  // from 93
"wait for obey command",
"wait in DitGetSeq()=0", 
"wait in DitGetSeq()=n", 
"wait for intask semaphore",
"wait for getpass_buf",

   ///  setup array
  // from 98 
"it is for SSA",
"it is for SQ2",
"it is for SQ1",
"it is for modulation",
"it is for modulation",
"it is for bias ramp function",
"non bias ramp function",
"found minimum modulation", 
"found maximum modulation",

// load
// from 107
"shutter closed",
"shutter open to hot source",
"shutter open to sky",
 " "
};

static char *servoStatus[]=
{
"NONESERVO",
"CABLECAL", 
"CABLECORET",
"HEATERSERVO",
"SQ1BIASSERVO",
"SQ1OPEN",
"SQ1SERVO",
"SQ2LOCK",
"SQ2OPEN",
"SQ2SERVO",
"SSALOCK",
"SSARAMP",
"SSARAMP1",
"TESBIASSERVO",
"TESTRANSIT",
"PIXELMON",
"TRKHEATER",
"TRKHEATERINIT",
"SQ2FBTRACK",
"SQ2OPEN4P",
"TOTALBASE",
"SQ2BIASING",  
"SQ1LOCK",
"BIASCHANGED",
"BIASNOCHANGE",
"SQ2SERVOSSALOCK",
"BLACKBODY",
"HEAT_SLOPE",
"SQ2FBGETPARA",
"TRI_UP",
"TRI_DOWN",
"TRI_MID",
" "
};


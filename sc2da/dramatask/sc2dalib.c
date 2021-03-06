/**
 * \file sc2dalib.c
 *
 * \brief collection of sub-functions for sc2da.c  
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
 *  $Log: sc2dalib.c,v $
 *  Revision 1.110  2013/08/09 21:11:19  bgorges
 *  Set SVFILE_FLAG default to 2. Removed unnecessary status ok checks.
 *
 *  Revision 1.109  2013/07/16 21:56:50  cwalther
 *  Took the heater tracking print of active pixels out of log file, added a log of the file name going into the QL structure, changed sigma for clipmean when getting squid-2 feedback
 *
 *  Revision 1.108  2012/10/19 00:56:46  cwalther
 *  changes to make sure all arrays have the exact same OBSID
 *
 *  Revision 1.107  2012/07/06 22:52:02  cwalther
 *  Finnaly figured out how to simplify the modetable so I could remove the huge if statements in the DRAMA task add modes to the table for pointing/focus with POL-2 and FTS-2 in the beam
 *
 *  Revision 1.106  2012/06/29 23:56:20  cwalther
 *  Changed the format of the heater tracking files, moved out of bounds checking to on the value we are commanding, detect non-working heater tracking pixels and throw them out
 *
 *  Revision 1.105  2012/06/05 23:24:54  cwalther
 *  Changed heater tracking failure messages to red
 *
 *  Revision 1.104  2012/05/16 19:15:14  cwalther
 *  Fixed a bug in the heater tracking stuff
 *
 *  Revision 1.103  2012/04/18 00:10:04  ryanb
 *  fixed heater tracking double-free/corruption bug
 *
 *  Revision 1.102  2012/03/28 02:39:38  cwalther
 *  Changes to handle INBEAM correctly (now shows shutter when in the dark)
 *
 *  Revision 1.101  2011/12/20 23:23:14  cwalther
 *  Implemented a style of heater tracking that ends after a fixed number of iterations
 *
 *  Revision 1.100  2011/11/22 18:50:35  cwalther
 *  Fix for screwup where the minimum heater tracking step became 20
 *
 *  Revision 1.99  2011/09/27 21:14:22  cwalther
 *  Moved the heater tracking header write so that I would print when you track to the value you remember in the dark, some changes to make remembering in dark work better
 *
 *  Revision 1.98  2011/09/02 19:16:20  bgorges
 *  Standardized entry checks on all fuctions/actions that should quietly return if entry status is bad.
 *
 *  Revision 1.97  2011/06/28 01:27:47  bgorges
 *  A commit to open changes.
 *
 *  Revision 1.96  2011/06/02 21:31:22  bgorges
 *  Checking all of fopen and fopen64 to make sure a file is actually opened.
 *
 *  Revision 1.95  2011/05/20 21:57:31  cwalther
 *  Moved time loging at end-of-SEQUENCE out of errors only into always, added MsgOut and end of CONFIGURE
 *
 *  Revision 1.94  2011/05/18 21:45:24  cwalther
 *  put end-of-line characters in all of the prints to the logfile that used to be MsgOuts
 *
 *  Revision 1.93  2011/05/17 20:18:29  cwalther
 *  These are all changes to cut down on the number of messages that have to go through the SCUBA2 task
 *
 *  Revision 1.92  2011/05/11 21:58:54  cwalther
 *  Added ability to use the median in heater tracking
 *
 *  Revision 1.91  2011/05/05 19:51:32  cwalther
 *  Fixed the zillion of messages problem when heater tracking fails
 *
 *  Revision 1.90  2011/04/29 21:13:13  cwalther
 *  Removed a lot of extraneous prints on MCE errors
 *
 *  Revision 1.89  2011/04/22 00:25:47  cwalther
 *  I put a considerable number of comments in this code and it has the changes made 
 *  to normalize thresholds when the number of samples changes
 *
 *  Revision 1.88  2011/03/21 20:45:14  cwalther
 *  Had to copy maxGainI from the right arry to the global version
 *
 *  Revision 1.87  2011/03/01 23:47:28  cwalther
 *  Added test for maximum gaini, set to zero if greater than maximum
 *
 *  Revision 1.1.1.1  2007/05/16 08:26:58  dkelly
 *  first insertion
 *
 */

#include "sc2math.h"   // for sc2math-linfit
#include "sc2headman.h"
#include "math.h"      // for pow
#include "jitXML.h"
#include "sc2sqopt_par.h"
#include "sc2sqopt.h"
#include "star/atl.h"

static char engMode[10];

//#define SAVE_ENG_DATA
#define LOG_TRKHEAT_PIXEL
#define LOG_TRKSQ2FB
  
//------------ List of PCI CHECK ------------------------
//-------------------------------------------------------
static  PCI_HEALTH healthStrt[]=
{
  {CHKXMEM, 0x00000001,"X Memory test:",0},
  {CHKYMEM, 0x00000002,"Y Memory test:",0},
  {CHKPMEM, 0x00000004,"P Memory test:",0},
  {CHKFIBRE,0x00000078,"Fibre test:",   0},

}; 
// PCI_HEALTH is in sc2da_struct.h
#define HEALTH_CHK_NO (int)( sizeof(healthStrt)/sizeof(PCI_HEALTH) )


static MCE_CARD_DEF  mcecardStrt[]=
{
  { "PSC",0},
  { "CC", 1},
  { "RC4",2},
  { "RC3",3},
  { "RC2",4},
  { "RC1",5},
  { "BC3",6},
  { "BC2",7},
  { "BC1",8},
  { "AC",9},
};

#define MCE_CARD_NO (int)( sizeof(mcecardStrt)/sizeof(MCE_CARD_DEF) )


// =======sc2dalib_a*******
//====================//

/**
 * \fn void sc2dalib_actionfileEnd(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *     int flag, StatusType *status)
 *
 * \brief function
 *  end Drama action and close file 
 *
 * \param con     SDSU_CONTXT struct pointer
 * \param myInfo dasInfoStruct_t struct pointer
 * \param flag   int  1: check IN_SEQ, 0, don't
 * \param status  StatusType.  given and return
 */
/*+ sc2dalib_actionfileEnd 
*/
void sc2dalib_actionfileEnd
(
SDSU_CONTEXT    *con,           
dasInfoStruct_t *myInfo,      
int             flag,
StatusType      *status
)
{
  long    in_sequence;

  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if (flag==1 && in_sequence != DA_MCE_NONE)
  {
    return;
  }
  sc2dalib_endAction(con,myInfo,status);
  utils_close_files(myInfo);
}


// =======sc2dalib_b*******
//====================//

/**
 * \fn void sc2dalib_batchInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  char *dateTime,StatusType *status)
 *
 * \brief functio
 *  get args for BATCH action and open files
 *
 * \param con       SDSU context structure
 * \param myInfo   dasInfo structure pointer
 * \param dateTime  string pointer to dateTime string
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_batchInit
*/
void sc2dalib_batchInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *dateTime,
StatusType            *status
)
{
  long          in_sequence,initialised;
  long          batchDelay;
  long          range[]={0,200};
  SdsIdType     argId;

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  // Do not continue if not initialized or in a data taking mode 
  SdpGeti("INITIALISED", &initialised, status);
  if( initialised==0)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_batchInit: drama task is not initialised" ); 
    return;
  } 
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status,"sc2dalib_batchInit: %s has not completed",
       seqStatus[myInfo->actionFlag]);
    return;
  }
  myInfo->actionFlag=BATCHACTION;
  argId = ( SdsIdType)DitsGetArgument();
  jitArgGetS(argId,"BATCH_FILE",1,NULL,"batch.txt",0,
	     FILE_LEN,myInfo->batchFile,NULL,status );
  jitArgGetS(argId,"DATA_FILE",2,NULL,"data.txt",0,
	     FILE_LEN,myInfo->dataFile,NULL,status );
  jitArgGetS(argId,"CMDREP_FILE",3,NULL,"cmdrep.txt",0,
	     FILE_LEN,myInfo->cmdrepFile,NULL,status );
  jitArgGetI(argId, "SVFILE_FLAG", 4, range, 2, 0, &myInfo->filesvFlag,
	     status );
  jitArgGetI(argId, "BATCH_DELAY", 5, range, 0, 0, &batchDelay,
	     status );
  if ( !StatusOkP(status) )
  {
    ErsOut(0,status, 
           "sc2dalib_batchInit: MCEBATCHGO arguments no good, defaults broken");
    batchDelay=0;
  }
  SdpPuti("BATCH_DELAY",batchDelay,status);

  utils_open_files(myInfo,DATFILEAPPEND,BATCHFILE,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_batchInit: utils_open_files failed.");
    return;
  }  
  fprintf(myInfo->fpLog,"\n<%s> CMD for sc2dalib__Batch (%s)\n",
           dateTime,myInfo->batchFile);
  myInfo->trkNo=0;
  myInfo->msgwrtPt=0;
  myInfo->msgreadPt=0;
  myInfo->astMapState=0;
  con->framechk=0;  
  con->pcidacount=0;    
  con->datacount=1;
  con->process.exit=0;
  SdpPuti("IN_SEQUENCE",DA_MCE_SINGLEDATA,status);
  con->process.whereabout=Dits_GetSeq0;
  con->process.framesetup= DA_MCE_SINGLEDATA;
}


/**
 * \fn void sc2dalib_batchParser(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  dasCmdInfo_t *mycmdInfo, char *dateTime,
 *  struct mcexml_struct *mceInxpt,StatusType *status)
 *
 * \brief function
 *  parse the batch file
 *
 * \param con       SDSU context structure
 * \param myInfo   dasInfo structure pointer
 * \param mycmdInfo dasCmdInfo_t pointer
 * \param dateTime  string pointer to dateTime string
 * \param mceInxpt  mcexml_struct  pointer
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_batchParser
*/
void sc2dalib_batchParser
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,   
dasCmdInfo_t          *mycmdInfo, 
char                  *dateTime,
struct mcexml_struct  *mceInxpt,
StatusType            *status
)
{
  int    j;
  int    *cmdPtr;
  size_t cmdlen;
  char   *byte,*insertPtr1,*batchptr;
  char   delimiter[]=" ", *cp, *token;
  char   inputcmd[CMD_LEN],tmpcmd[CMD_LEN*2];

  if (*status != STATUS__OK) return;

  mycmdInfo->cmdBufSize=MCEBLK_SIZE;
  batchptr=tmpcmd;
  byte=(char*)&mycmdInfo->reply;

  jitDebug(2,"sc2dalib_batchParser: dataFile(%s) cmdrepFile(%s)\n", 
             myInfo->dataFile,myInfo->cmdrepFile);
  while(1)
  {
    if(  fgets(inputcmd,CMD_LEN,myInfo->fpBatch) !=NULL )
    {
      cmdlen=strlen(inputcmd);
      if ( (strchr(inputcmd,'#') !=NULL) || (strcmp(inputcmd,"\n")==0) )
      {
        // all comments recorded in the datafile
        myInfo->batchFlag=BATCH_CMNT;
        fprintf(myInfo->fpData,"%s\n",inputcmd);
      }
      else if ( strchr(inputcmd,'<') !=NULL)
      {
        j=1;
        while(inputcmd[cmdlen-j] !='<' )
        {
          inputcmd[cmdlen-j]=' ';  j++;
        }
        inputcmd[cmdlen-j]=' ';
        batchptr=stpcpy(batchptr, inputcmd);
      }
      else
      {
        inputcmd[cmdlen-1]=' ';
        batchptr=stpcpy(batchptr, inputcmd);
        batchptr=tmpcmd;
        strcpy(mycmdInfo->mceCmd,tmpcmd);
        mcexml_translate(tmpcmd,mceInxpt,mycmdInfo->cmdBufSize,
                          mycmdInfo->cmdBuf,status);
        if ( *status ==DITS__APP_ERROR)
        {
          utils_msg(myInfo,
                    "sc2dalib_batchParser: failed to call mcexml_translate",
                    "",USE_ERSREP,status);
          utils_msg(myInfo, "sc2dalib_batchParser: the cmd is <%s>",
                    mycmdInfo->mceCmd, USE_ERSREP, status);
          myInfo->batchFlag=BATCH_TRANSFAIL;
          return;
        }
        else
        {
          myInfo->trkNo++;
          // check if the cmd is a GO
          insertPtr1=(char*)mycmdInfo->cmdBuf+8;   
          cmdPtr=(int*)insertPtr1;
          if(*cmdPtr==MCE_GO)
          {
            // set up whichRC for parameter shared Mempry
            cp = strdupa (mycmdInfo->mceCmd); // Make writable copy. 
             // token => "command" before "space" 
            token = strtok (cp, delimiter);
            token = strtok (NULL, delimiter);
            sc2dalibsetup_whichRC(myInfo,token,status);
            if ( *status ==DITS__APP_ERROR)
            {
              utils_msg(myInfo,
                 "sc2dalib_batchParser:  not recognised(%s) as MCE_WHICHCARD",
                  mycmdInfo->mceCmd,USE_ERSREP,status);
              myInfo->batchFlag=BATCH_TRANSFAIL;
            }
            else
              // we need to reschedule, cmdbuf back from return
              myInfo->batchFlag=BATCH_GO;
            // reset datacount    
            con->datacount=1;
            return;
          }
          else
          {
             jitDebug(2,"sc2dalib_batchParser: the %ldth command= %s\n", 
                    myInfo->trkNo,mycmdInfo->mceCmd);
            // check the dummy setting for external stuff
            cp = strdupa (mycmdInfo->mceCmd); // Make writable copy. 
             // token => "command" before "space" 
            token = strtok (cp, delimiter);
            token = strtok (NULL, delimiter);
            if (strcmp(token, EXTERNAL_CARD) ==0)
            {
               // the EXTERL one only for storing setting in mceInxpt
               // don't send out
            }
            else
            {
              myInfo->batchFlag=BATCH_CMD;
              jitDebug(2,"sc2dalib_batchParser: the %ldth command= %s\n", 
                    myInfo->trkNo,mycmdInfo->mceCmd);
 
              mce_cmd(con,myInfo,mycmdInfo,dateTime,status);
              if (*status != STATUS__OK)
              {
                myInfo->batchFlag=BATCH_CMDFAIL;
                utils_msg(myInfo,
                 "sc2dalib_batchParser: failed to call mce_cmd","",
                 USE_ERSREP,status);
                return;
              }
            }
          }
        }
      }
    }
    else
    {
      myInfo->batchFlag=BATCH_END;   
      return;
    }
  }
}


// =======sc2dalib_c*******
//====================//

/**
 * \fn void sc2dalib_callbkdrvMsg(USED4DITS *argptr, StatusType *status)
 *
 * \brief function
 *  call-back routine: there is a message from DRV_USER_MESG FIFO
 *
 * \param argptr  USED4DITS pointer 
 * \param status  StatusType pointer. 
 *
 * if  *status != STATUS__OK, report error.
 *
 * This routine is called from the DRAMA altInLoop when there is something 
 * in the  DRV_USER_STRUCT FIFO  Read the message and display or do 
 * something.   message from driver
 * 
 */

/*+ sc2dalib_callbkdrvMsg 
*/
void sc2dalib_callbkdrvMsg
(
USED4DITS      *argptr,    
StatusType     *status
)
{

  time_t   tm;
  int      rval;
  static  DRV_USER_STRUCT drvuserMsg;
 
  if (!StatusOkP(status)) return;

  rval=master_read_drvmsg(argptr->con, &drvuserMsg);
  if( rval <0 )
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status,"sc2dalib_callbkdrvMsg: master_read_drvmsg failed");
    return;
  }
  if(drvuserMsg.action !=GET_VERSION)
  {
    tm = time(NULL);
    fprintf (argptr->myInfo->fpLog,"%s",asctime(localtime(&tm)) );

    MsgOut(status,
      "\nsc2dalib_callbkdrvMsg:========= SDSU_DRIVER INFORMATION ===========");
    utils_msg(argptr->myInfo,
       "==============driver ERROR ========================","",USE_PRINT, status);    
    if( drvuserMsg.action ==ERR_RESTART)
    {
      utils_msg(argptr->myInfo,"  %s",drvuserMsg.error,USE_MSGOUT, status);
      utils_msg(argptr->myInfo,"  %s",drvuserMsg.msg,USE_MSGOUT,status);
    }
    else
    {
      utils_msg(argptr->myInfo," %s",drvuserMsg.error,USE_MSGOUT, status);
      utils_msg(argptr->myInfo," %s",drvuserMsg.msg,USE_MSGOUT, status);
      utils_msg(argptr->myInfo," %s",drvErrors[drvuserMsg.action],
                        USE_MSGOUT,status);
    }
    if( drvuserMsg.action >=ERR_SEMTMO)
    {
      utils_msg(argptr->myInfo,
         "Fibre link down, or partial packets received","",USE_MSGOUT,status);
    }
  }
  else
  { // copy the version info into the global one
    strcpy(argptr->myInfo->drvVersion,drvuserMsg.msg);
    if ( argptr->myInfo->actIndex != -1 )
    {
      DitsSignalByIndex ( argptr->myInfo->actIndex, 0, status );
    }
    else
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,
             "sc2dalib_callbkdrvMsg: not valid argptr->myInfo->actIndex");
    }  
  }
}


/**
 * \fn void sc2dalib_callbkMsg(USED4DITS *argptr, StatusType *status)
 *
 * \brief function
 *  call-back routine: there is a message in USER_USER_MSG FIFO
 *
 * \param argptr  USED4DITS pointer 
 * \param status  StatusType pointer. 
 *
 * if  *status != STATUS__OK, report error.
 *
 * This routine is called from the DRAMA altInLoop when there is something 
 * in the USER_USER_MSG FIFO buffer. Read the values from the FIFO and
 * ask  DRAMA  to reschedule the associated action. message from dhandler
 */

/*+ sc2dalib_callbkMsg 
*/
void sc2dalib_callbkMsg
(
USED4DITS    *argptr,    
StatusType   *status
)
{
  int            rval;
  DRAMA_INNERMSG *intermsg;
  PAR_SHARED     *parentPtr;
  

#ifdef SAVE_ENG_DATA
  FILE          *coaddFp;
  char          coaddfileName[FILE_LEN];
#endif 

  // call-back can not send Msgout/Ers, for it is not an action
  if (!StatusOkP(status)) return;

  parentPtr=(PAR_SHARED *)argptr->myInfo->parShm;

  // myInfo is a pointer in argptr, so use ->actIndex
  // the same for intermsg
  intermsg=&argptr->intermsg[argptr->myInfo->msgwrtPt];

  if ( argptr->myInfo->actIndex != -1 )
  {
    // read  information from fifo
    if( (rval=master_read_othermsg(argptr->con, intermsg)) <0 )
    {
       printf("sc2dalib_callbkMsg: Error- sdsu_read_msg failed");
       return;
    }
    else
    {
      if (argptr->myInfo->actionFlag==SEQACTION)
      {
        jitDebug(8,"sc2dalib_callbkMsg: msg.reason=%s, msgwrtPt=%d\n",
               seqStatus[intermsg->reason],(int)argptr->myInfo->msgwrtPt);
      }
      // recycle=16, only change the msgwrtPt if we reschedule
      argptr->myInfo->msgwrtPt++;
      argptr->myInfo->msgwrtPt=argptr->myInfo->msgwrtPt & 0x000F ;
      // reschedule SEQ action
      DitsSignalByIndex ( argptr->myInfo->actIndex, 0, status );
    } 
  }
}


/**
 * \fn void sc2dalib_changeSQVal(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceInxpt, char *dateTime, ARRAYSET *setup,
 *  StatusType *status)
 *
 * \brief function
 *  depending on setup->slopSelect[5], send array cmd (ramping)
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param mceInxpt  struct mcexml_struct pointer
 * \param dateTime  dateTime string pointer   
 * \param setup    ARRAYSET structure pointer      
 * \param status StatusType.  given and return
 *
 * setup->slopSelect[5] >0  selsect  which sq to ramp
 * setup->slopSelect[4]     initial Val  // 0
 * setup->slopSelect[3]     step         // 160
 * setup->slopSelect[2]     NO of Step   // 409
 * setup->slopSelect[8]    >=0 at which step apply sq2fbopt
 *
 */
/*+ sc2dalib_changeSQVal
 */
void sc2dalib_changeSQVal
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                  *dateTime,
ARRAYSET              *setup,
StatusType            *status
)
{
  int      j, stepVal;
  char     tmp[30];
  char     ssafbCmd[]  ="wb bc1 flux_fb";
  char     sq2fbCmd[]  ="wb bc2 flux_fb";
  char     sq2biasCmd[]="wb bc3 flux_fb";
  char     sq1biasCmd[]="wb ac on_bias";
  //  char     ssabiasCmd[]="wb rc1 sa_bias";

  static dasCmdInfo_t myCmd;

  static int  inx, insertInx, stepUPDW;

  if (*status != STATUS__OK) return;

  if (myInfo->trkNo ==SQ2FB_UPDATE_WAIT )
  {
    inx=0;  insertInx=0;stepUPDW=1;
  }
  
  if ( setup->slopSelect[5] ==1 )
    sprintf (myCmd.mceCmd, "%s ",ssafbCmd ); 
  else if( setup->slopSelect[5] ==2)  
    sprintf (myCmd.mceCmd, "%s ",sq2fbCmd ); 
  else if ( setup->slopSelect[5] ==3)  
    sprintf (myCmd.mceCmd, "%s ",sq2biasCmd ); 
  else if ( setup->slopSelect[5] ==4)  
    sprintf (myCmd.mceCmd, "%s ",sq1biasCmd ); 
  else
  {
    *status=DITS__APP_ERROR;
    ErsRep (0, status,"slopSelect5]=%d not assigned",setup->slopSelect[5]);
    return;
  }
  
  if ( inx == setup->slopSelect[8]  && setup->slopSelect[8] >=0)
  {
    if (insertInx ==0)
    {
      MsgOut(status," apply sq2fbOpt now");
      sc2dalib_sendarrayCmd(con,myInfo,mceInxpt,dateTime,myCmd.mceCmd,setup->sq2fdbkOpt,COL_NUM,status);
    }
    insertInx ++;
    if (insertInx ==SQ2FB_UPDATE_WAIT)
      inx ++;
    return;
  }
  // change inx back to it was  
  if (insertInx ==SQ2FB_UPDATE_WAIT)
  {
    inx --;
    insertInx ++;
  }
         
  if ( stepUPDW ==1 )
    stepVal= setup->slopSelect[4] + setup->slopSelect[3] * inx;
  else
    stepVal= setup->slopSelect[4] + setup->slopSelect[3] * (setup->slopSelect[2]- inx);

  if ( inx %20 == 0)
    MsgOut(status,"itrkNo=%d inx=%d stepVal=%d %s",myInfo->trkNo,inx, stepVal,
           myCmd.mceCmd);

  for ( j=0; j<COL_NUM; j++ )
  {
    sprintf ( tmp, " %d", stepVal );
    strcat(myCmd.mceCmd,tmp); 
  }
  sc2dalib_sendCmd(con,myInfo,&myCmd,mceInxpt,dateTime,status);
  if ( setup->slopSelect[5] ==3)
  {
    // when ssabias and sq2bias change, we need to give time for circuit 
    // to settle down
    sleep(WAIT_TIME_AFTER_BIAS);;
  }

  inx ++;
  if ( (inx == setup->slopSelect[2]) && (stepUPDW==0) )
  {
    if ( setup->slopSelect[5] ==1 )
      sprintf (myCmd.mceCmd, "%s ",ssafbCmd ); 
    else if( setup->slopSelect[5] ==2)  
      sprintf (myCmd.mceCmd, "%s ",sq2fbCmd ); 
    else if ( setup->slopSelect[5] ==3)  
      sprintf (myCmd.mceCmd, "%s ",sq2biasCmd ); 
    else if ( setup->slopSelect[5] ==4)  
      sprintf (myCmd.mceCmd, "%s ",sq1biasCmd ); 
    
    // stop calling this function again from sc2da.c
    setup->slopSelect[5]=-1;
    MsgOut(status," apply sq2fbOpt now");
    sc2dalib_sendarrayCmd(con,myInfo,mceInxpt,dateTime,myCmd.mceCmd,setup->sq2fdbkOpt,COL_NUM,status);
  }
  else if (inx == setup->slopSelect[2] )
  {
     inx=0; stepUPDW=0; insertInx=0;
  }
}


/**
 * \fn void sc2dalib_chksdsID(dasInfoStruct_t *myInfo,
 *     StatusType *status)
 *
 * \brief function
 *  check DSD ids for QL, if !=NULL delete and free
 *  Removes the unneeded data structures from the QL structure
 *
 * \param   myInfo    dasInfoStruct_t structure pointer
 * \param  status     StatusType.  given and return
 *
 */
/*+ sc2dalib_chksdsID 
*/
void sc2dalib_chksdsID
(
dasInfoStruct_t  *myInfo,   
StatusType       *status
)
{
  if (!StatusOkP(status)) return;

  /* Remove DATA_ARRAY, FITS and IMAGE */

  if ( myInfo->qldataId !=0)
  {
     SdsDelete(myInfo->qldataId,status);
     SdsFreeId(myInfo->qldataId,status);
     myInfo->qldataId = 0;
  }
  if ( myInfo->fitsId !=0)
  {
     SdsDelete(myInfo->fitsId,status);
     SdsFreeId(myInfo->fitsId,status);
     myInfo->fitsId =0;
  }
  if ( myInfo->imageId !=0)
  {
     SdsDelete(myInfo->imageId,status);
     SdsFreeId(myInfo->imageId,status);
     myInfo->imageId = 0;
  }

  /* Remove the filename */

  if ( myInfo->filenameId !=0)
  {
     SdsDelete(myInfo->filenameId,status);
     SdsFreeId(myInfo->filenameId,status);
     myInfo->filenameId = 0;
  }
}


/**
 * \fn void sc2dalib_coadddatawriteDisp(SDSU_CONTEXT *con, 
 *  dasInfoStruct_t *myInfo, char *scanfile, char *coaddPtr,
 *  double *imagedataPtr,  StatusType *status)
 *
 * \brief function
 *  write coadd data into a file and if debug flag set, display some QL data
 *
 * \param con           SDSU context structure pointer 
 * \param myInfo       dasInfoStruct_t pointer
 * \param scanfile      char pointer for scan file name if SCAN
 * \param coaddPtr      char pointer for ith coadd/reconstr entry if STARE/DREAM
 * \param imagedataPtr  double pointer for IMAGE entry
 * \param status        StatusType.  given and return
 *
 */
/*+ sc2dalib_coadddatawriteDisp   
 */
void sc2dalib_coadddatawriteDisp
(
SDSU_CONTEXT      *con,
dasInfoStruct_t   *myInfo,
char              *scanfile,
char              *coaddPtr,
double            *imagedataPtr,    
StatusType        *status
)
{
  int    *qlentPtr;
  double *timeS; 
  int    obsMode;
  char   coaddfileName[FILE_LEN];

#ifdef SAVE_ENG_DATA  
  FILE  *coaddFp;
#endif
  
  if (!StatusOkP(status)) return;

  obsMode=myInfo->parshmPtr->obsMode;
  qlentPtr=(int *)coaddPtr;
  timeS=(double *)(qlentPtr+1);

  // for writing coadd to a temp file
  sprintf (coaddfileName, "%s/%s", getenv ( "DATADIR" ),
                  COADDFILE_DAS  ); 
#ifdef SAVE_ENG_DATA  
  if((coaddFp = fopen(coaddfileName, "a")) == NULL )
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status, "Error- sc2dalib_coadddatawriteDisp: failed to open file %s", coaddfileName); 
    return;
  }
#endif

  if( obsMode==OBS_STARE || obsMode==OBS_DREAM )
  {
#ifdef SAVE_ENG_DATA  
    if( obsMode !=OBS_DREAM) 
      fwrite(coaddPtr, 1, sizeof(QL_STRUCT), coaddFp);
    else
      fwrite(coaddPtr, 1, sizeof(RECONSTR_STRUCT), coaddFp);

    fclose (coaddFp);
#endif    
    if (myInfo->debuglvl==2) 
    { 
      int   j,i,numcard2disp=4;

      jitDebug(2,
         "sc2dalib_coadddatawriteDisp: the framenum and timestamp  in QL_STRUCT\n");
      jitDebug(2,"%11d %11.1f \n", *qlentPtr,*timeS);

      jitDebug(2,"the first 2 IMAGE data from each card in QL_STRUCT\n");
      for (j=0;j<numcard2disp;j++)
      {
        for ( i=0; i<2; i++) 
         jitDebug(2,"%11.1f ", *(imagedataPtr + j*8*40 +i) );
        jitDebug(2,"\n");  
      } 
    }
  }
  else if (myInfo->parshmPtr->obsMode ==OBS_SCAN)
  {
#ifdef SAVE_ENG_DATA  
    fwrite(coaddPtr, 1, sizeof(QL_SCAN), coaddFp);
    fclose (coaddFp);
#endif
    jitDebug(2,"sc2dalib_coadddatawriteDisp: the scanfilename =%s\n",
             scanfile);
  }
}


/**
 * \fn void sc2dalib_configInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo,struct mcexml_struct *mceInxpt,  
 *  char *dateTime, char *mode, char * config, StatusType *status)
 *
 * \brief function
 *  get args for CONFIGURE action
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfo structure pointer
 * \param dateTime  char pointer for dateTime string
 * \param mode      char pointer for observation mode
 * \param config    char pointer for configfile string
 * \param status   StatusType pointer    
 *
 */
/*+ sc2dalib_configInit
*/
void sc2dalib_configInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
struct mcexml_struct  *mceInxpt, 
char                  *dateTime,
char                  *mode,
char                  *config,
StatusType            *status
)
{
  long           in_sequence,initialised;
  long           procno,range[]={1,150000};
  SdsIdType      argId, jitId;
  char           file[FILE_LEN]=" ";
  char           heatVal[]="rb bc1 bias 1";

  if (*status != STATUS__OK) return;
  
  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  // Do not continue if not initialized or in a data taking mode 
  SdpGeti("INITIALISED", &initialised, status);
  if( initialised==0)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_configInit: drama task is not initialised" ); 
    return;
  } 
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status,"sc2dalib_configInit: %s has not completed",
       seqStatus[myInfo->actionFlag]);
    return;
  }
  SdpPuti("CONFIGURED", 0, status);
  myInfo->actionFlag =CONFIGACTION;

  argId = DitsGetArgument();

  /* Copy the SDS structure containing the parsed XML into the DRAMA parameter
   system */
  jitArgGetXML ( argId, "CONFIGURATION", 1, NULL, FILE_LEN, file, &jitId, status );
  SdpPutStruct ( "CONFIGURATION", jitId, 0, 1, status);   
  if (!StatusOkP(status))
  {
    ErsRep (0, status,"sc2dalib_configInit: SdpPutStruct CONFIGURATION failed");
    return;
  }

  // if it can not find FRAME_NO2PROC,it will take any thing from arg40, 
  // so we use a big arg to fool ?
  jitArgGetI( argId, "FRAME_NO2PROC", 40, range, 0, 0, &procno, status );
  if (!StatusOkP(status))
    return;

  if ( procno==0)
  { // it is for RTSC 
    sc2dalib_configInitRTSC(myInfo,config,status);
    ArgPutString(myInfo->statId,"DATE", myInfo->baseFile, status);
    jitDebug(2,"sc2dalib_configInit: DATE =%s\n", myInfo->baseFile);
  }
  else
  { // it is for ENG 
    sc2dalib_configInitENG(myInfo,mode,config,status);
  }

  // read the current heater setting and store as the nominal
  mce_read(con, myInfo, mceInxpt, heatVal, &myInfo->nominalPixelHeat,1, status);
  if (!StatusOkP(status)) 
    {
      ErsRep(0,status, "sc2dalib_configInit: mce_read failed to read heater value"); 
      return;
    }

}


/**
 * \fn void sc2dalib_configInitENG(dasInfoStruct_t *myInfo, 
 *  char *mode, char * config, StatusType *status)
 *
 * \brief function
 *  get args for CONFIGURE action for ENG
 *
 * \param myInfo   dasInfo structure pointer
 * \param mode      char pointer for observation mode
 * \param config    char pointer for configfile string
 * \param status   StatusType pointer    
 *
 */
/*+ sc2dalib_configInitENG
*/
void sc2dalib_configInitENG
(
dasInfoStruct_t       *myInfo,    
char                  *mode,
char                  *config,
StatusType            *status
)
{
  long               procno,range[]={1,150000};
  DitsArgType          argId;

  if (*status != STATUS__OK) return;

  // get observation mode. Can be one of DREAM, STARE, SCAN, FLATFIELD, 
  // DARK,SETUP or SKYDIP.  pipelineICD p9

  argId = DitsGetArgument();

  jitArgGetS(argId, "CONFIGURATION", 1, NULL, "config.xml", 0,
	     FILE_LEN, config, NULL, status );
  jitArgGetS(argId, "OBS_MODE", 2, NULL, "DARK", 0, FILE_LEN,
	     mode, NULL, status );
  jitArgGetI( argId, "FRAME_NO2PROC", 3, range, 200, 0, &procno, status );
  // 200 as default ?
  SdpPuti("FRAME_NO2PROC",procno,status);
  myInfo->parshmPtr->procNo=(int)procno;

  // make up things for ENG to write NDF file
  myInfo->obsNo=00001;
  sprintf(myInfo->baseFile,"%s",getenv("ARRAYDATADIR"));
}


/**
 * \fn void sc2dalib_configInitRTSC(dasInfoStruct_t *myInfo, 
 *     char *config,  StatusType *status)
 *
 * \brief function
 *  get args for CONFIGURE action from RTSC
 *
 * \param myInfo   dasInfo structure pointer
 * \param config    char pointer for configfile string
 * \param status   StatusType pointer    
 *
 */
/*+ sc2dalib_configInitRTSC
*/
void sc2dalib_configInitRTSC
(
dasInfoStruct_t       *myInfo,    
char                  *config,
StatusType            *status
)
{
  SdsIdType   argId;
  SdsIdType   initId;
  char obsNo[10];

  if (*status != STATUS__OK) return;

  argId = DitsGetArgument();
 
  jitArgGetXML ( argId, "CONFIGURATION", 1, NULL, FILE_LEN, config, &initId, status );
  jitDebug(4," configXML=%s",config); 
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_configInitRTSC: failed to get configure.xml");
    return;
  }
  jitArgGetS(argId, "DATE", 2, NULL, "20007", 0, FILE_LEN,
	     myInfo->baseFile, NULL,status);
  jitArgGetS(argId, "OBSNUM", 3, NULL, "00001", 0, 10, obsNo, NULL,status);
  myInfo->obsNo=atoi(obsNo);

  jitArgGetS(argId, "UTCSHORT", 4, NULL, "unknown", 0, 
	     OBSIDSS_SIZE, myInfo->utcshort, NULL, status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_configInitRTSC: failed to get variables."); 
    return;
  }

  // reset 
  myInfo->parshmPtr->qlframeNum=0;
  myInfo->parshmPtr->heat_track_num=0;
  SdsFreeId ( initId, status );
}


/**
 * \fn void sc2dalib_configInitSet(dasInfoStruct_t *myInfo, 
 *  char *dateTime, char *mode, char * config, StatusType *status)
 *
 * \brief function
 *  setup some parameter from CONFIGURE action for sharedMem
 *
 * \param myInfo   dasInfo structure pointer
 * \param dateTime  char pointer for dateTime string
 * \param mode      char pointer for observation mode
 * \param config    char pointer for configfile string
 * \param status   StatusType pointer    
 *
 */
/*+ sc2dalib_configInitSet
*/
void sc2dalib_configInitSet
(
dasInfoStruct_t       *myInfo,  
char                  *dateTime,  
char                  *mode,
char                  *config,
StatusType            *status
)
{
  long     procno;
  FILE     *tmpFp;
  char      arrayName[40],cmd[FILE_LEN];

  if (!StatusOkP(status)) return;

  // get procno first, then change according to other setting
  if ( myInfo->engFlag==RTSC_MODE)
    ArgGeti(myInfo->statId,"STAREADD", &procno,status);

  // myInfo->parShmPtr is   PAR_SHARED * 
  if ( strcmp(mode,"DREAM") == 0 )
  {
    myInfo->parshmPtr->obsMode=OBS_DREAM;
    if ( myInfo->engFlag==RTSC_MODE)
    { 
      // only initial smuPattern, will update at beginning of seq
      procno=myInfo->parshmPtr->dreamRound*myInfo->parshmPtr->smuPattern;  
    }
  }
  else if ( strcmp(mode,"STARE") == 0)
  {
    myInfo->parshmPtr->obsMode=OBS_STARE;
  }
  else if ( strcmp(mode,"SCAN") == 0 )
  {
    myInfo->parshmPtr->obsMode=OBS_SCAN;
    if ( myInfo->engFlag==RTSC_MODE)
      ArgGeti(myInfo->statId,"SCANWRITE", &procno,status);
  }
  else
  {
    *status=DITS__APP_ERROR; 
    ErsRep(0,status,"sc2dalib_configInitSet: OBS_MODE (%s) not recognized",
           mode); 
     return;
  }
  if ( myInfo->engFlag==RTSC_MODE)
    myInfo->parshmPtr->procNo=(int)procno;

  myInfo->procNum=myInfo->parshmPtr->procNo;

  jitDebug(2,"sc2dalib_configInitSet: set for %s OBS\n",
             seqStatus[myInfo->parshmPtr->obsMode]);

  ArgPutString(myInfo->statId,"OBS_MODE",mode, status);
  ArgPuti(myInfo->statId,"OBS_NO",myInfo->obsNo, status);
  
  ArgGetString(myInfo->statId,"ARRAY_CNNTED", 40, arrayName, status);

  // have ok and data dir ready
  //if ( parsharedEnt->engFlag==RTSC_MODE)
  {
    sprintf(cmd, "mcesetdatadir-rtsc %s %s %05d",
           arrayName,myInfo->baseFile,myInfo->obsNo);
    jitDebug(2,"sc2dalib_configInitSet:script-cmd=%s\n",cmd);
    system(cmd);
  }
  // get ok filename
  sprintf(myInfo->okfileName,"%s%s_%05d.ok",
     arrayName,myInfo->baseFile, myInfo->obsNo);

  // pass related info to DH task for writing data
  // initialise subscanNo=1, qlframeNum=0 pass to dh, so that, each time 
  // a file is written, it will be increased by one
  jitDebug(2,"sc2dalib_configInitSet: configXML=%s\n",config); 
  
  myInfo->parshmPtr->obsNo=myInfo->obsNo;
  myInfo->parshmPtr->subscanNo=myInfo->subscanNo=1;
  myInfo->parshmPtr->frameHeader=FRAMEHEADER_NUM;
  strcpy(myInfo->parshmPtr->arrayName,arrayName);
  strcpy(myInfo->parshmPtr->baseFile,myInfo->baseFile);
  strcpy(myInfo->parshmPtr->flatfile,myInfo->flatField.file);
  strcpy(myInfo->parshmPtr->configFile,config);

  // fresh /$OKROOT/tmp.ok for dhtask to append each subscanNo
  jitDebug(2,"sc2dalib_configInitSet: opening tmp.OK\n"); 
  sprintf(cmd,"%s/tmp.ok",getenv("OKROOT"));
  if((tmpFp = fopen(cmd, "w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "Error- sc2dalib_configInitSet: failed to open file %s", cmd); 
      return;
    }
  fclose(tmpFp);

  jitDebug(2,"sc2dalib_configInitSet: calling createQLDS\n"); 
  utils_create_QLSDS(myInfo, status );
  if ( !StatusOkP(status) )
    return;

  fprintf(myInfo->fpLog,"\n<%s> CMD for sc2dalib__Config",dateTime);
  fprintf(myInfo->fpLog,
          "\nOBS_MODE <%s> FRAME_NO2PROC<%d> CONFIGURATION<%s>\n",
          mode, myInfo->parshmPtr->procNo,config);

  SdpPuti("SETUP",0,status);
  SdpPuti("IN_SEQUENCE",DA_MCE_NONE,status);
  SdpPuti("CONFIGURED", 1, status);
}

 
/**
 * \fn void sc2dalib_cnvtbatch2Cmd(char *setFile,char *cmd, StatusType *status)
 *
 * \brief function
 *  convert batch file into a single cmd 
 *
 * \param  myInfo   dasInfoStruct_t poiter
 * \param  setFile  char string for the file
 * \param  cmd     char string, return cmd
 * \param  status StatusType.  given and returned
 *
 */
/*+ sc2dalib_cnvtbatch2Cmd
*/
void sc2dalib_cnvtbatch2Cmd
(
char            *setFile,
char            *cmd,
StatusType      *status
)
{
  FILE   *fp;
  char   tmpfile[FILE_LEN];
  char   input[FILE_LEN];
  size_t cmdlen;
  char   *batchptr;
  int    iscmdLine=1;
  int    j;

  if (*status != STATUS__OK) return;

  batchptr=cmd;
  fp=NULL;
  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ), setFile );
  if((fp = fopen(tmpfile,"r")) == NULL)
  { 
    *status=DITS__APP_ERROR;
    ErsRep(0, status,"sc2dalib_cnvtbatch2Cmd: failed to open %s",tmpfile);
    return;
  }

  while(1)
  {
    if(  fgets(input,CMD_LEN,fp) !=NULL )
    {
      cmdlen=strlen(input);
      // skip all comments 
      if ( (strchr(input,'#') !=NULL) || (strcmp(input,"\n")==0) )
      {
        //printf("commentline=%s\n",input);
      }
      // any line has character <
      else if ( strchr(input,'<') !=NULL)
      {
        j=1;
        // search until '<' change to space
        while(input[cmdlen-j] !='<' )
        {
          input[cmdlen-j]=' ';  j++;
        }
        input[cmdlen-j]=' ';
        
        if ( iscmdLine ==1)
        {
          iscmdLine=0;
          batchptr=stpcpy(batchptr, input);
          //printf("cmdWord=%s\n",input);
        }
        else
        {
          batchptr=stpcpy(batchptr, input);
          //printf("argWord=%s\n",input);
        }
      }
      else  // the last values
      {
        input[cmdlen-1]=' ';
        batchptr=stpcpy(batchptr, input);
        batchptr=cmd;
        break;
      }
    }
    else
      break;
  }
  fclose(fp);
}


// =======sc2dalib_d*******
//====================//

/**
 * \fn void sc2dalib_dispInfo(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo,
 *  int inSeq, StatusType *status)
 *
 * \brief functio
 *  display das information
 *
 * \param con     SDSU_CONTEXT  structure pointer
 * \param myInfo dasInfoStruct_t structure pointer
 * \param inSeq  int number indication for which state
 * \param status  StatusType.  given and return
 *
 */
/*+ sc2dalib_dispInfo
*/
void sc2dalib_dispInfo
(
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
int             inSeq,
StatusType      *status
)
{
  int rtsC, smuC, ptcsC, scuba2C, fts2C, pol2C;
  char  shellcmd[]="checkForSc2dadh";

  if (*status != STATUS__OK) return;

 //  display con->process parameter

  if(myFpLog != NULL)
    {
      fprintf(myFpLog," .actionFlag =(%s)\n", seqStatus[myInfo->actionFlag]);
      fprintf(myFpLog," .IN_SEQ     =(%s)\n", seqStatus[inSeq]);
      fprintf(myFpLog," .reason     =(%s)\n", seqStatus[con->process.reason]);
      fprintf(myFpLog," .seqstatus  =(%s)\n", seqStatus[con->process.seqstatus]);
      fprintf(myFpLog," .obsmode    =(%s)\n", seqStatus[myInfo->parshmPtr->obsMode]);
      fprintf(myFpLog," .actIndex   =(%ld)\n", myInfo->actIndex);
      fprintf(myFpLog," mceTrkheaterNo(%d)\n", myInfo->trkNo);
      fprintf(myFpLog," frameSeenby(DH)(%ld)\n", con->framechk);
      fprintf(myFpLog," frameReceived by(driver/datatask)(%ld)\n",  con->pcidacount);
      fprintf(myFpLog," frameCountSet in(driver)(%ld)\n", con->datacount);
      fprintf(myFpLog," chdataCount (DH)(%d)\n", myInfo->parshmPtr->chdataCount);
      fprintf(myFpLog," fstateCount (DH)(%d)\n", myInfo->parshmPtr->fstateCount);
      fprintf(myFpLog," darkrowCount(DH)(%d)\n", myInfo->parshmPtr->darkrowCount);
      fprintf(myFpLog," coaddscandreamCount(DH)(%d)\n", myInfo->parshmPtr->coaddscandreamCount);
      sc2headman_getcounts(&rtsC, &smuC, &ptcsC, &scuba2C, &fts2C, &pol2C);
      fprintf(myFpLog," RTS count=   %7d\n", rtsC);
      fprintf(myFpLog," SMU count=   %7d\n", smuC);
      fprintf(myFpLog," PTCS count=  %7d\n", ptcsC);
      fprintf(myFpLog," SCUBA2 count=%7d\n", scuba2C);
      fprintf(myFpLog," FTS2 count=  %7d\n", fts2C);
      fprintf(myFpLog," POL2 count=  %7d\n", pol2C);
    }

 if ( system ( shellcmd ) != 0 )
   {
     MsgOut(status,"!!!!! sc2dadh died you must do an UNLOAD_INST - LOAD_INST cycle ");
     fprintf(myFpLog,"sc2dadh died you must do an UNLOAD_INST - LOAD_INST cycle ");
   }
}


/**
 * \fn void sc2dalib_downld2pciInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *   char *file, char *dateTime, StatusType *status)
 *
 * \brief function
 *  get args for DWLOADDSP action and open files
 *
 * \param con      SDSU context structure
 * \param myInfo   dasInfo structure pointer
 * \param file      char pointer for the DSP *.lod file
 * \param dateTime  string pointer to dateTime string
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_downld2pciInit
*/
void sc2dalib_downld2pciInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *file,
char                  *dateTime,
StatusType            *status
)
{
  long               in_sequence;
  DitsArgType          argId;

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  // Do not continue if not initialized or in a data taking mode 
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
     ErsRep (0, status, "sc2dalib_downld2pciInit: %s has not completed",
	     seqStatus[myInfo->actionFlag]);
     return;
  }
  myInfo->actionFlag =DSPACTION;
  argId = DitsGetArgument();
  jitArgGetS(argId, "DSP_FILE", 1, NULL, "dsp.lod", 0,
	     FILE_LEN, file, NULL, status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_downld2pciInit: failed to get DSP_FILE"); 
    return;
  }
  fprintf(myInfo->fpLog,"\n<%s> CMD from sc2dalib__Download2PCI <%s>\n",
          dateTime,file);
 con->process.framesetup= DA_MCE_NONE;
}


/**
 * \fn void sc2dalib_download2PCI(SDSU_CONTEXT *con, char *lodfile,
 *  StatusType *status)
 *
 * \brief function
 *  download DSP *.lod file into PCI DSP
 *
 * \param con      SDSU_CONTEXT structure pointer
 * \param lodfile  char string for the *.lod file
 * \param status   StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
/*+ sc2dalib_download2PCI
*/
void sc2dalib_download2PCI
(
SDSU_CONTEXT *con, 
char        *lodfile,
StatusType *status
)
{
  FILE     *fp;
  long     memAddr;
  long     memValue, positOne;
  int      loadFinish,dataFinish,i;
  char     memTypeInt='Y', *cmdPCI="WRITE";
  char     idName[CMD_LEN],memType[CMD_LEN];
  static int      displayNo=0;
  static PCI_CMD  pci_cmd;

  if (*status != STATUS__OK) return;
   
  if((fp = fopen(lodfile,"r")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep(0, status, "sc2download2PCI: Error- failed to open file %s",
	     lodfile); 
      return;
    }
  loadFinish=0;
  while (!loadFinish)
  {
    fscanf(fp, "%s", idName);
    if ( strcmp("_START",idName) == 0)
    {
      // skip the file and version number
      for (i=0;i<6;i++) fscanf(fp, "%s", idName);
    }
    else if ( strcmp("_DATA",idName) == 0)
    { 
      displayNo=0;
      fscanf(fp, "%s", memType);
      if (strcmp(memType,"X")==0 )            memTypeInt='X';
      else if (strcmp(memType,"Y")==0 )       memTypeInt='Y';
      else if (strcmp(memType,"P")==0 )       memTypeInt='P';

      fscanf(fp, "%s", idName);
      memAddr=(long)strtoul ( idName,0,16);

      if( ( strcmp(memType,"P")==0 && (memAddr < 0x800) ) )
        MsgOut(status,"MemType=%s, MemAddr=%#lX skipping",memType,memAddr);
      else
        MsgOut(status,"MemType=%s, MemAddr=%#lX downloading.......",
               memType,memAddr); 

      // go through the current _DATA section until the next _DATA or
      // _SYMBOL
      dataFinish=0;
      if(  ( strcmp(memType,"P")==0 && (memAddr >= 0x800) ) ||
             strcmp(memType,"Y")==0 || strcmp(memType,"X")==0
        )
      {
        while(!dataFinish)
        {  
          positOne=ftell(fp);
          fscanf(fp, "%s", idName);
          *status=STATUS__OK;
          if( strcmp("_DATA",idName) == 0 || strcmp("_SYMBOL",idName) == 0)
          {
            #ifdef DEBUGDW
            if(displayNo%16==0)
              MsgOut(status,"f-positBF=%ld ",positOne);
            #endif

            dataFinish=1;     
            fseek(fp,positOne,SEEK_SET);
          }
          else
          {
            memValue=(long)strtoul ( idName,0,16);
            *status=
               sdsu_set_pcicmd(&pci_cmd,cmdPCI, memTypeInt,memAddr,&memValue);
            if( *status !=SDSU_OK)
            { 
              *status=DITS__APP_ERROR;
              ErsRep (0, status, 
                      "sc2dalib_download2PCI: Error-sdsu_set_pcicmd failed"); 
              fclose(fp);
              return;
            }  
            if( (*status =sdsu_command_pci(con,&pci_cmd)) !=SDSU_OK)
            {
              *status=DITS__APP_ERROR;         
              sc2dalib_pcierrRep(&pci_cmd,cmdPCI,status);
              ErsRep (0, status, 
                      "sc2dalib_download2PCI: Error-sdsu_command_pci failed");
              fclose(fp);
              return;
            }   
            memAddr++;
            displayNo++;
          }
        }
      }
    }
    else if ( strcmp("_SYMBOL",idName) == 0)
    {
    }
    else if ( strcmp("_END",idName) == 0)
    {  
      loadFinish=1;
    }
  }
  fclose(fp);
}


// =======sc2dalib_e*******
//====================//
/**
 * \fn void sc2dalib_endAction(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  StatusType *status)
 *
 * \brief function
 *  set args in con and myInfo for action end 
 *
 * \param con     SDSU_CONTXT structure pointer
 * \param myInfo dasInfoStruct_t structure pointer
 * \param status  StatusType pointer.  given and return
 *
 */
/*+ sc2dalib_endAction - 
*/
void sc2dalib_endAction
(
SDSU_CONTEXT    *con,           
dasInfoStruct_t *myInfo,       
StatusType      *status
)
{
  StatusType     tmp;

  tmp=STATUS__OK;
  con->datacount=0;
  fflush(myInfo->fpLog);

  con->process.whereabout=WAIT_OBEY;
  SdpPuti("IN_SEQUENCE",DA_MCE_NONE,&tmp);  
  con->process.framesetup=DA_MCE_NONE;
 
  myInfo->filesvFlag=0;
  myInfo->actIndex = -1;
  if( myInfo->debuglvl==2)
    sc2dalib_pcistatusRep(con,myInfo,&tmp);
}


// =======sc2dalib_f*******
//====================//
/**
 * \fn void sc2dalib_frametakeInit(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  int startSeq, int endSeq, struct mcexml_struct *mceInxpt, 
 *  char *dateTime, StatusType *status)
 *
 * \brief function
 *  send start and end frame Numbers to MCE, send GO to MCE.
 *  inform driver to expect data
 *
 * \param con       SDSU context structure pointer 
 * \param myInfo    dasInfoStruct_t pointer
 * \param startSeq  int: start frame 
 * \param endSeq    int: end frame
 * \param mceInxpt  struct mcexml_struct pointer
 * \param dateTime  dateTime string pointer         
 * \param status    StatusType.  given and return
 *
 */
/*+ sc2dalib_frametakeInit   
 */
void sc2dalib_frametakeInit
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
int                   startSeq,
int                   endSeq,
struct mcexml_struct  *mceInxpt,
char                  *dateTime,
StatusType            *status
)
{
  dasCmdInfo_t frameCmd;
  if (!StatusOkP(status)) return;

  errno=0;
  con->datacount=0;  // inform driver NO of data expected
  myInfo->trkNo=0;  
  myInfo->msgwrtPt=0;
  myInfo->msgreadPt=0;
  con->framechk=0;    
  con->process.exit=0;
  con->process.seqstatus=SEQ_NOACTION;  
  con->process.whereabout=Dits_GetSeq0;

  con->datacount=endSeq-startSeq + 1;   
  myInfo->parshmPtr->frameNum=myInfo->numFrame;
  myInfo->numFrame=con->datacount;   
  jitDebug(2,"sc2dalib_frametakeInit: dataCount(frameNum)=%d\n",
          myInfo->numFrame);

  sprintf(frameCmd.mceCmd, "WB cc ret_dat_s %d %d",startSeq, endSeq);
  // send cmd and get reply, save them to logfile donot save to other
  // cmdreply file specified by user
  sc2dalib_sendCmd(con,myInfo,&frameCmd,mceInxpt,dateTime,status); 
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_frametakeInit(1): sc2dalib_sendCmd %s failed",
           frameCmd.mceCmd); 
    return;
  }

  con->pcidacount=0;
  strcpy(frameCmd.mceCmd, myInfo->goCmd);
  sc2dalib_sendCmd(con,myInfo,&frameCmd,mceInxpt,dateTime,status); 
  if ( !StatusOkP(status) )
  { 
     ErsRep(0,status,"sc2dalib_frametakeInit(2): sc2dalib_sendCmd %s failed",
            frameCmd.mceCmd); 
     return;
  }
  if( myInfo->debuglvl==2)
    mce_results(myInfo,&frameCmd,status);
  con->process.reason=FRAME_WAITING;
}

/**
 * \fn void sc2dalib_fluxJump(dasInfoStruct_t *myInfo,
 *  ARRAYSET *arrayset, int *array,  int ch, StatusType *status)
 *
 * \brief function
 *  do fluxi jump 
  *
 * \param myInfo  dasInfoStruct_t pointer
 * \param setup   ARRAYSET structure pointer
 * \param array   int array, 
 * \param ch      int 
 * \param status StatusType.  given and return
 *
 */
/*+ sc2dalib_fluxJump
 */
void sc2dalib_fluxJump
(
dasInfoStruct_t  *myInfo,
ARRAYSET         *setup,
int              *array,
int              ch,
StatusType       *status
)
{
  if (*status != STATUS__OK) return;

  if (array[ch] !=0)  // masked one, kept =0
  {
    if (setup->fluxPeriod[ch] !=0)
    {
      if ( array[ch] >=0) 
      { 
         while( array[ch] > setup->fluxPeriod[ch] )
           array[ch] -=setup->fluxPeriod[ch];
      }
      else
      {
        while( array[ch] < 0 )
          array[ch] +=setup->fluxPeriod[ch];  
      }
    }
  }
}     


// =======sc2dalib_g*******
//====================//

/**
 * \fn void sc2dalib_getcmdBuf(dasInfoStruct_t  *myInfo,
 *   dasCmdInfo_t *mycmdInfo, struct mcexml_struct  *mceInxpt,
 *   StatusType *status)
 *
 * \brief functaion
 *  translate MCE command into command struture 
 *
 * \param  myInfo   dasInfoStruct_t  pointer
 * \param mycmdInfo  dasCmdInfo structure pointer
 * \param mceInxpt   mcexml_struct struct pointer
 * \param status     StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
/*+ sc2dalib_getcmdBuf 
 */
void sc2dalib_getcmdBuf
(
dasInfoStruct_t        *myInfo,
dasCmdInfo_t           *mycmdInfo,
struct mcexml_struct   *mceInxpt,
StatusType             *status
)
{
  char  cmd[CMD_LEN];

  if (!StatusOkP(status)) return;

  mycmdInfo->cmdBufSize=MCEBLK_SIZE;
  strcpy(cmd,mycmdInfo->mceCmd);
  mcexml_translate(cmd,mceInxpt,mycmdInfo->cmdBufSize,
                                mycmdInfo->cmdBuf, status);
  if ( *status ==DITS__APP_ERROR)
  {     
    utils_msg(myInfo,
       "sc2dalib_getcmdBuf: Error- failed to  call mcexml_translate","",
       USE_ERSREP,status);
    utils_msg(myInfo,
       "sc2dalib_getcmdBuf: the cmd is <%s>",mycmdInfo->mceCmd,USE_ERSREP,status);
    return ;
  }
}


// =======sc2dalib_h*******
//====================//
/**
 * \fn void sc2dalib_healthStatus(int which, int health, StatusType *status)
 *
 * \brief function
 *  test the Health bits after START PCI application 
 *
 * \param which   int: flag for which part to check
 * \param health  int:health status word return by PCI
 * \param status  StatusType.  given and return
 *
 *
 */
/*+  sc2dalib_healthStatus
 */
void sc2dalib_healthStatus
(
int        which,
int        health,
StatusType *status
)
{
  if (!StatusOkP(status)) return;

  if (which <3)
  {
    if( (health & healthStrt[which].mask) >0)
    {
      healthStrt[which].memAdr=( health >>8 ) & 0x0000FFFF;
      MsgOut(status,"%s failed at address <%#lX>",healthStrt[which].name,
             healthStrt[which].memAdr);
    }
    else 
      MsgOut(status,"%s passed !!",healthStrt[which].name);
  }
  else
  {
    if( (health & healthStrt[which].mask) >0)
    {
      if( (health & 0x00000008) >0)
      {
        MsgOut(status,
            "%s failed - no bytes receicved, check fibre connection",
               healthStrt[which].name);
        MsgOut(status,
          "(PCI card's transmitter should be looped back to PCI card's receiver)!");
      }
      if( (health & 0x00000010) >0)
        MsgOut(status,
          "%s failed - too many bytes receicved!",healthStrt[which].name);
      if( (health & 0x00000020) >0)
        MsgOut(status,
          "%s failed - not enough bytes receicved!",healthStrt[which].name);
      if( (health & 0x00000040) >0)
        MsgOut(status,
          "%s failed - corrupt data receicved!",healthStrt[which].name);
    }
    else 
      MsgOut(status,"%s passed !!",healthStrt[which].name);
  }
}


// =======sc2dalib_i*******
//====================//

/**
 * \fn void sc2dalib_initInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  char *xmlfilename, StatusType *status)
 *
 * \brief function
 *  get args for INITIALISE action and set default values 
 *
 * \param con          SDSU context structure
 * \param myInfo      dasInfo structure pointer
 * \param xmlfilename  char pointer for the mce xml file
 * \param status       StatusType     
 *
 */
/*+ sc2dalib_initInit
*/
void sc2dalib_initInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *xmlfilename,
StatusType            *status
)
{
  int        i;
  long       in_sequence;
  SdsIdType  argId;
  SdsIdType  jitId;
  char       file[FILE_LEN]=" ";

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);
  
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status, 
      "sc2dalib_initInit: %s has not completed",seqStatus[myInfo->actionFlag]);
    return;
  } 
  myInfo->actionFlag=INITACTION;

  
  argId = DitsGetArgument();

  /* Copy the SDS structure containing the parsed XML into the DRAMA parameter
     system */
  jitArgGetXML ( argId, "INITIALISE", 1, NULL, FILE_LEN, file, &jitId, status );
  SdpPutStruct ( "INITIALISATION", jitId, 0, 1, status);   
  if (!StatusOkP(status))
  {
    ErsRep (0, status,"sc2dalib_initInit: SdpPutStruct INITIALISATION failed");
    return;
  }

  // The environment variable SC2ENGMODE will be 0 or 1 for RTS mode
  // and 2 or greater is always engineering mode

  if ( (strcmp(getenv ( "SC2ENGMODE" ),"0") == 0) || (strcmp(getenv ( "SC2ENGMODE" ),"1") == 0))
  { // it is RTSC mode
    sc2dalib_initInitRTSC(con,myInfo,xmlfilename,status);
    myInfo->parshmPtr->engFlag=myInfo->engFlag=RTSC_MODE;
  }
  else
  { // it is ENG  mode
    sc2dalib_initInitENG(con,myInfo,xmlfilename,status);
    myInfo->parshmPtr->engFlag=myInfo->engFlag=ENG_MODE;
  }
  if (!StatusOkP(status))
    return;
  

  // temporarily use /tmp/user/logcmd. after arrayID, it changes back to
  // $SC2DIR/array/logfile/logcmddate 
  sprintf (myInfo->logfileName, "/tmp/%s/logcmd",getenv("USER"));

  strcat(myInfo->logfileName,myInfo->Date);
  utils_fclose(&(myInfo->fpLog));
  if((myInfo->fpLog = fopen64(myInfo->logfileName,"w")) == NULL )
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status, "Error- sc2dalib_initInit: failed to open file %s",
	    myInfo->logfileName); 
    return;
  }
  // myFpLog is a global so all routines can write in the log file
  myFpLog = myInfo->fpLog;

  jitDebug(2, "sc2dalib_initInit: the sc2dalib_logfile is %s\n",
           myInfo->logfileName);

  for (i=0;i<4;i++)
    myInfo->parshmPtr->whichRC[i]=0;

  // reset stripchFlag here as time variable for stripchart
  stripchFlag=0;
}


/**
 * \fn void sc2dalib_initInitENG(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  char *xmlfilename, StatusType *status)
 *
 * \brief function
 *  get args for INITIALISE action and set default values for ENG
 *
 * \param con          SDSU context structure
 * \param myInfo      dasInfo structure pointer
 * \param xmlfilename  char pointer for the mce xml file
 * \param status       StatusType     
 *
 */
/*+ sc2dalib_initInitENG
*/
void sc2dalib_initInitENG
(
SDSU_CONTEXT          *con, 
dasInfoStruct_t       *myInfo, 
char                  *xmlfilename,
StatusType            *status
)
{
  long         range[]={0,3};
  char         protcl[FILE_LEN]; 
  char         dataform[FILE_LEN]; 
  DitsArgType  argId;

  if (*status != STATUS__OK) return;

  argId = DitsGetArgument();
  
  //get parameters and store them in the paramater system   

  jitArgGetI(argId, "SIMULATE",2,range, 0, 0, &myInfo->simFlag, status );
  SdpPuti("SIMULATE", myInfo->simFlag, status);
  jitArgGetS(argId, SC2DADATAFORMAT, 3, NULL, "BINARY",0,
	     FILE_LEN, dataform, NULL, status );
  if ( !StatusOkP(status) )
  {
    ErsOut(0,status, 
       "sc2dalib_initInitENG: DATA_FORMAT_NAME is not given, use BINARY as default");
    SdpPutString(SC2DADATAFORMAT,"BINARY",status);
  }
  else
  {
    if ( strcmp(dataform,"TEXT")==0 )
    {
      myInfo->dataFormat=MCE_TEXT_FORM;
      con->dataform=MCE_TEXT_FORM;
      jitDebug(2,
       "sc2dalib_initInitENG: data will be saved in TEXT format 32 data/row\n");
    }
    else if ( strcmp(dataform,"TEXT2")==0 )
    {
      myInfo->dataFormat=MCE_TEXT_FORM2;
      con->dataform=MCE_TEXT_FORM2;
      jitDebug(2,
       "sc2dalib_initInitENG: data will be saved in TEXT format one data/row\n");
    }
    else
      jitDebug(2,"sc2dalib_initInitENG: data will be saved in BINARY format\n");

   jitDebug(2,
       "sc2dalib_initInitENG: dataformat:%s\n",dataform);

    SdpPutString(SC2DADATAFORMAT,dataform,status);
  }

  jitArgGetS(argId,"PROTOCOL",4, NULL, "NEW",0, FILE_LEN, protcl, 
          NULL,status );
  if ( strcmp(protcl,"NEW")==0 )
    {
      con->flag.protoCOL=NEWPROTO;
      jitDebug(2,"sc2dalib_initInitENG: use new protocol\n");
    }
  else 
    {
      con->flag.protoCOL=OLDPROTO;
      jitDebug(2,"sc2dalib_initInitENG: use old protocol\n");
    }
  jitArgGetI(argId, "DIVBY_N",5,range, 0, 0, &myInfo->divbynFlag, status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status, "sc2dalib_initInitENG: DIVBY_N is not given");
    return;
  }
  sc2dalib_initInitReadXML(con,myInfo,xmlfilename,status); 
  if (!StatusOkP(status))
    ErsRep(0,status, "sc2dalib_intInitRTSC: failed to call sc2dalib_initInitReadXML"); 
}


/**
 * \fn void sc2dalib_initInitRTSC(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  char *xmlfilename, StatusType *status)
 *
 * \brief function
 *  get args for INITIALISE action and set default values for RTSC
 *
 * \param con          SDSU context structure
 * \param myInfo      dasInfo structure pointer
 * \param xmlfilename  char pointer for the mce xml file
 * \param status       StatusType     
 *
 */
/*+ sc2dalib_initInitRTSC
*/
void sc2dalib_initInitRTSC
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *xmlfilename,
StatusType            *status
)
{
  SdsIdType argId;
  long      range[]={0,3};

  if (*status != STATUS__OK) return;

  // fix a few for ENG use
  myInfo->divbynFlag=0; 
  con->flag.protoCOL=OLDPROTO;
 
  jitDebug(2,"sc2dalib_initInitRTSC: Hey, it is RTSC  now !\n");

  argId = DitsGetArgument();
  jitArgGetI(argId, "SIMULATE",2,range, 0, 0, &myInfo->simFlag, status );
  if (!StatusOkP(status))
  {
   ErsRep(0,status, "sc2dalib_intInitRTSC: error in getting SIMULATE"); 
   return;
  }     
  SdpPuti("SIMULATE", myInfo->simFlag, status);

  sc2dalib_initInitReadXML(con,myInfo,xmlfilename,status); 
  if (!StatusOkP(status))
    ErsRep(0,status, "sc2dalib_intInitRTSC: failed to call sc2dalib_initInitReadXML"); 
}


/**
 * \fn void sc2dalib_initInitReadXML(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  char *xmlfilename, StatusType *status)
 *
 * \brief function
 *  get XML for INITIALISE action and set default values 
 *
 * \param con          SDSU context structure
 * \param myInfo      dasInfo structure pointer
 * \param xmlfilename  char pointer for the mce xml file
 * \param status       StatusType     
 *
 */
/*+ sc2dalib_initInitReadXML
*/
void sc2dalib_initInitReadXML
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *xmlfilename,
StatusType            *status
)
{
  SdsIdType argId;
  SdsIdType initId;         //id for INITIALISE structure 
  SdsIdType sc2initId;      //id for SCUBA2_INIT structure 
  SdsIdType instrumentId;  //  id for instrument structure 
  SdsIdType sc2extraId = 0;
  SdsIdType sc2taskId = 0;
  SdsIdType subarrayId = 0;
  SdsIdType wavebandId = 0;
  SdsIdType id = 0;
  SdsCodeType code;

  DitsPathInfoType  info;
  DitsPathType      path;
  double            timeOut;
  
  int   j,k,passInx=0;
  long  ndims;            /* number of dimensions in subarray */
  long  port,x,y, initFlag;
  long  stareAdd, scanWrite,dreamRound, along;
  unsigned long dims[7];       
  unsigned long indims[7];
  long darkHeaterI;
  long maxGainI;
  long anInt;

  char  name[FILE_LEN],health[10];
  char  wavelength[FILE_LEN];
  char  arrayName[10], overRide[10];

  if (*status != STATUS__OK) return;

  // the initialise.xml need to be in JITXML_DIR
  argId = DitsGetArgument();

  jitArgGetXML ( argId, "INITIALISE", 1, NULL, sizeof(name), name,
                 &initId, status );
  if (!StatusOkP(status))
  {
   ErsRep(0,status, "sc2dalib_intInitReadXML: jitArgGetXML failed"); 
   return;
  }
  strcpy(myInfo->initFile,name);
  MsgOut(status," initialXML=%s",name);

  // us SdsList(  initId,status) to print the initId SDS;
  // initID is top SDS, the rest is sub sds. get sub-sds 
  SdsFind ( initId, "SCUBA2_INIT", &sc2initId, status );
  SdsFind ( sc2initId, "INSTRUMENT", &instrumentId, status );

  dims[0] = dims[1] = 0;
  passInx=0;
  
  // from instrumentId, get subArray SDS array 
  SdsFind ( instrumentId, "subArray", &subarrayId, status );
  SdsInfo ( subarrayId, name, &code, &ndims, dims, status );  
  if (!StatusOkP(status))
  {
    ErsRep(0,status, 
    "sc2dalib_intInitReadXML: find subArray failed in passInx(%d)",passInx); 
    return;
  }
  
  // initialised with big invalid value
  for ( j=0; j<MAX_SUBARRAY; j++ )
    myInfo->subArray[j].mceport = 0xff;
 
  // we only have one dimension of subArray SDS, ie. ndims==1
  // we have 8 subArray ie. dims[0]=8, dims[x]=0
  // search each element of subArray, extract the values 
  for ( j=0; j<dims[0]; j++ )
  {
    indims[0] = j + 1;
    passInx=0;

    if (!StatusOkP(status)) break;

    // from subarrayId ( SDS array), get each  array[x] 
    SdsCell ( subarrayId, 1, indims, &id, status ); 
    // from each array[x], get item, id is the arrayName
    ArgGetString(id, "id",10,myInfo->subArray[j].id,status);

    if (StatusOkP(status))
    {
      passInx++;
      ArgGeti(id, "mceport", &port, status );
      myInfo->subArray[j].mceport=(int)port;
      if (StatusOkP(status))
      {
        passInx++;
        ArgGetString(id, "chipId",40,myInfo->subArray[j].chipid,status);
        if (StatusOkP(status))
        {
          passInx++;
          ArgGetString(id, "band",10,myInfo->subArray[j].band,status);
          if (StatusOkP(status))
          {
            passInx++;
            ArgGetString(id, "health",10,myInfo->subArray[j].health,status);
            if (StatusOkP(status))
            {
              ArgGeti(id, "x", &x, status );
              myInfo->subArray[j].x=(int)x;
              if (StatusOkP(status))
              {
                passInx++;
                ArgGeti(id, "y", &y, status );
                myInfo->subArray[j].y=(int)y;
                if (StatusOkP(status))
                {
                  passInx++;
                  ArgGetString(id, "task",40,myInfo->subArray[j].task,status);
                  if (StatusOkP(status))
                  {
                    passInx++;
                    ArgGetString(id, "flatfile",
                                 FILE_LEN,myInfo->subArray[j].flatfile,status);
                    if (StatusOkP(status))
                    {
                      passInx++;
                      ArgGetString(id, "dreamweightfile",
                                  FILE_LEN,myInfo->subArray[j].weightfile,status);
                    }
                  }
                }
              }
            }
	  }
	}
      }
      else
        break;
    }
     
    if (!StatusOkP(status))
    {
      ErsRep(0,status, 
      "sc2dalib_intInitReadXML: find subArray[%d] item failed in passInx(%d)",
       j,passInx);
      return;
    }
    SdsFreeId ( id, status );       
  }
    SdsFreeId ( subarrayId, status );

  SdsFind ( sc2initId, "SCUBA2_EXTRAS", &sc2extraId, status );
  /* printf("\nReading SCUBA2 extras\n"); */

  // from instrumentId, get indivalual item
  ArgGetString(sc2extraId, "MCEXML", FILE_LEN,xmlfilename,status);
  if (StatusOkP(status))
  {
    jitDebug(2,"sc2dalib_initInitReadXML: mce.XML=%s\n",xmlfilename);
    passInx++;
    ArgGeti(sc2extraId, "SCANWRITE", &scanWrite, status );
    if (StatusOkP(status))
    {
      passInx++;
      ArgGeti(sc2extraId, "DREAMROUND", &dreamRound, status );
      if (StatusOkP(status))
      {
        passInx++;
        ArgGeti(sc2extraId, "STAREADD", &stareAdd, status );
	if (StatusOkP(status))
	  {
	    passInx++;
	    ArgGeti(sc2extraId, "HEATTRK_STYLE", &along, status );
	    myInfo->heatTrkStyle = along;
	    if (StatusOkP(status))
	      {
		passInx++;
		ArgGeti(sc2extraId, "INITFLAGS", &initFlag, status );
		myInfo->initFlag=(int)initFlag;
	      }
	  }
      }
    }
  }

  if (!StatusOkP(status))
  {
    ErsRep(0,status, 
          "sc2dalib_intInitReadXML: %d_th ArgGet** from sc2extraId failed",
          passInx); 
    return;
  }

  /* Always make the first heater track read and save the reference values */
  anInt = UPDATE_HEATER_TRACKING_REFERENCE_BIT | myInfo->heatTrkStyle;
  SdpPuti("HTTRACK_FLAG", anInt, status);

  // initialize all of the override flags to false to start with
  for ( j=0; j<MAX_SUBARRAY; j++ )
    {
      myInfo->subArray[j].override = 0;
    }

  // from sc2extraId, get subArray SDS ID, then ask for info about it
  SdsFind ( sc2extraId, "subArray", &subarrayId, status );
  SdsInfo ( subarrayId, name, &code, &ndims, dims, status );
  if (!StatusOkP(status))
    {
      ErsRep(0,status, 
	     "sc2dalib_intInitReadXML: finding subArray failed in scuba2Extras"); 
      return;
    }

  // search each element of subArray, extract the values 
  for ( j=0; j<dims[0]; j++ )
    {
      indims[0] = j + 1;
      passInx=0;

      // from subarrayId get an id for the current cell  
      SdsCell (subarrayId, 1, indims, &id, status );
 
      if (StatusOkP(status))
	{
	  passInx++;
	  ArgGetString(id, "id",10,arrayName,status);
	  if (StatusOkP(status))
	    {
	      passInx++;
	      // overRide will be either "TRUE" or "FALSE"
	      ArgGetString(id, "codeRdrOveride", 10, overRide, status );
	      if (StatusOkP(status))
		{
		  passInx++;
		  // Read the dark heater value 
		  ArgGeti(id, "darkHeaterI", &darkHeaterI, status );
		  if (StatusOkP(status))
		    {
		      passInx++;
		      // Read the magnitude of maximum allowable integrator gain (in the MCE PPL) 
		      ArgGeti(id, "maxGainI", &maxGainI, status );
		    }
		}
	    }
	}
      if (!StatusOkP(status))
	{
	  ErsRep(0, status, "sc2dalib_intInitReadXML: %d_th ArgGet** from sc2extras-subArray failed",
		 passInx);
	  return;
	}

      // Now search the subArray name list to find this particular array
      for(k=0; k<MAX_SUBARRAY; k++)
	{

	  // Is this the subArray I have the info for?
	  if(strcmp(arrayName, myInfo->subArray[k].id) == 0)
	    {
	      // If override is TRUE change its value to 1 
	      if (strcmp(overRide, "TRUE") == 0)
		{
		  // MsgOut(status, "Over ride optical code reader for array %s", arrayName);
		  myInfo->subArray[k].override = 1;
		}
	      // Set this particular array's dark heater current to that which was read
	      myInfo->subArray[k].darkHeaterI = darkHeaterI;
	      // Set this particular array's maximim integrator gain to that which was read
	      myInfo->subArray[k].maxGainI = maxGainI;
	    }
	}
      // free this cell's ID
      SdsFreeId(id, status);
    }
  // Free the subArray ID
  SdsFreeId ( subarrayId, status );

  /* Find and decode the waveBand information in the instrument XML */
  
  SdsFind ( instrumentId, "waveBand", &wavebandId, status );
  SdsInfo ( wavebandId, name, &code, &ndims, dims, status );
  if (!StatusOkP(status))
  {
    ErsRep(0,status,"sc2dalib_intInitReadXML: find waveband failed ");
    return;
  }

  //SdsList(  wavebandId,status);
  // we only have one dimension of waveband SDS, ie. ndims==1
  // we have 2 waveband ie. dims[0]=2, dims[x]=0
  for ( j=0; j<dims[0]; j++ )
  {
    indims[0] = j + 1;
    passInx=0;
    if (!StatusOkP(status)) break;

    SdsCell ( wavebandId, 1, indims, &id, status ); 

    ArgGetString(id, "band",10,myInfo->waveBand[j].band,status);
    if (StatusOkP(status))
    {
      passInx++;
      ArgGetString(id, "label",10,myInfo->waveBand[j].label,status);
      if (StatusOkP(status))
      {
        passInx++;
        ArgGetString(id, "units",10,myInfo->waveBand[j].units,status);
        if (StatusOkP(status))
        {
          passInx++;
          ArgGetd(id, "centre",&myInfo->waveBand[j].centre,status);
          if (StatusOkP(status))
          {
            passInx++;
            ArgGetd(id, "width",&myInfo->waveBand[j].width,status);
          }
        }
      }
    }
    SdsFreeId ( id, status );     
  }
  SdsFreeId ( wavebandId, status );
  if (!StatusOkP(status))
  {
    ErsRep(0,status, 
      "sc2dalib_intInitReadXML: find waveband[%d] item failed in passInx(%d)",
       j,passInx); 
     return;
  }
 

  // =====  get SC2TASK,  save to myInfo and jitPathget if engMose =0,1
  sprintf(engMode,"%s",getenv("SC2ENGMODE"));
  if ( strcmp (engMode, "0") == 0  || strcmp (engMode, "1") == 0 )
  { 
    SdsFind ( sc2extraId, "SC2TASK", &sc2taskId, status );
    SdsInfo ( sc2taskId, name, &code, &ndims, dims, status );
    if (!StatusOkP(status))
    {
      ErsRep(0,status,"sc2dalib_intInitReadXML: find SC2TASK failed ");
      return;
    }
    info.MessageBytes = 500;
    info.MaxMessages = 10;
    info.ReplyBytes = 500;
    info.MaxReplies = 10;
    timeOut = 10;

    for ( j=0; j<dims[0]; j++ )
    {
      indims[0] = j + 1;
      SdsCell ( sc2taskId, 1, indims, &id, status ); 
      ArgGetString(id, "id",20,myInfo->taskList[j],status);
      if (StatusOkP(status))
      {
        ArgGetString(id, "health",10,health,status);
        if (StatusOkP(status))
        {
          if ( strcmp( health,"ON") == 0 )
          {
            /* MsgOut(status,"sc2dalib_intInitReadXML: find Path for %s",
	       myInfo->taskList[j]); */
            jitPathGet(myInfo->taskList[j], &info, timeOut,&path, status);
          }
        }
      }
      SdsFreeId ( id, status );
      if (!StatusOkP(status)) break;
    }

    SdsFreeId ( sc2taskId, status );
    if (!StatusOkP(status))
    {
       ErsRep(0,status,"sc2dalib_intInitReadXML: jitPathGet (%s) status 0x%x",
              myInfo->taskList[j], (int)*status);
    }
  }
  SdsFreeId ( sc2extraId, status );
  SdsFreeId ( instrumentId, status );
  SdsFreeId ( sc2initId, status );
  SdsFreeId ( initId, status );

  jitDebug(2,"STAREDADD=%d\n",stareAdd);
  jitDebug(2, "SCANWRITE=%d\n",scanWrite );
  jitDebug(2, "DREAMROUND=%d\n",dreamRound);
  jitDebug(2, "INITFLAGS=%d\n",initFlag);
  jitDebug(2, "MCEXML=%s\n",xmlfilename); 
 
  ArgPuti(myInfo->statId, "STAREADD",stareAdd, status );
  ArgPuti(myInfo->statId, "SCANWRITE",scanWrite, status );
  ArgPuti(myInfo->statId, "DREAMROUND",dreamRound, status );
  ArgPuti(myInfo->statId, "INITFLAGS",initFlag, status );
  ArgPutString(myInfo->statId, "MCEXML",xmlfilename,status); 
  ArgPutString(myInfo->statId, "WAVELENGTH",wavelength,status); 
  SdpPutString(SC2DAMCEXMLFILE, xmlfilename, status);
  myInfo->parshmPtr->dreamRound=dreamRound;
  myInfo->parshmPtr->smuPattern=SMU_PATTERN; // initial to 64
}



/**
 * \fn void sc2dalib_initgetArrayID(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceinxPtr, char *name,  char *mcexml, int *mcePort,
 *  StatusType *status)
 *
 * \brief function: 
 *  read back MCE parameters Val and assign ARRAY ID 
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param mceinxPtr  mcexml_struct pointer
 * \param name      char  pointer for array   G$R
 * \param mcexml    char  pointer for mcexml file
 * \param mcePort   int  pointer for mcePort No
 * \param status    StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 */
/*+ sc2dalib_initgetArrayID
*/
void sc2dalib_initgetArrayID
(
SDSU_CONTEXT          *con,
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceinxPtr,
char                  *name,
char                  *mcexml,
int                   *mcePort,
StatusType            *status
)
{
  char dateTime[40];
  char dataform[40];

  if (*status != STATUS__OK) return;

  long tmp;
  int  i,j;
  char  arrayID[]="rb cc array_id 1 ";

  // it needs a existing logfile for writing error if any
  // doesnt need cmdrepfile as it is NULL,if this is ENG_MODE
  // then logfile is $LOGFILEDIR/logcmddate
  //but for RTSC_MODE, at moment it is /tmp/user/logcmddate,

  /* If *mcePort is 0xFF, read the port number from the MCE, else use the port number given */
 
  if(*mcePort == 0xFF)
    {
      mce_read(con,myInfo,mceinxPtr,arrayID, mcePort,1,status);
      if (!StatusOkP(status)) 
	return;
    }

  for (i=0; i<MAX_SUBARRAY; i++)
    if(myInfo->subArray[i].mceport==*mcePort)
      break;
  if( i==MAX_SUBARRAY )
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status, "sc2dalib_initgetArrayID: mceport is wrong"); 
      return;
    }
    
  myInfo->darkHeaterI = myInfo->subArray[i].darkHeaterI;
  myInfo->maxGainI = myInfo->subArray[i].maxGainI;
  printf("This is my dark heater current: %d\n",myInfo->darkHeaterI);
  strcpy(name, myInfo->subArray[i].id);
  strcpy(myInfo->chipId,myInfo->subArray[i].chipid);
  strcpy(myInfo->flatField.file,myInfo->subArray[i].flatfile);  
  strcpy(myInfo->parshmPtr->dreamweightFile,myInfo->subArray[i].weightfile);  
  for (j=0; j<2; j++)
  {
    if (strcmp (myInfo->subArray[i].band, myInfo->waveBand[j].band) ==0 )
    {
      myInfo->wavelen=myInfo->waveBand[j].centre;
      strcpy(myInfo->filter, myInfo->waveBand[j].label);
    }
  }
  MsgOut(status, " subArray[%d].id=%s wavelength=%4.1e filter=%s",
        i,myInfo->subArray[i].id,myInfo->wavelen,myInfo->filter);
  
  // use this to get row_len, num_rows 
  mce_read_frame_rate(con,myInfo,mceinxPtr,&tmp,status);
  if ( *status ==DITS__APP_ERROR)
    return;

  utils_get_time(DAS_DATETIME,dateTime);

  ArgPutString(myInfo->statId, "ARRAY_CNNTED",name, status);
  // change logfile to array related
  sprintf (myInfo->logfileName, "%s/%s/logfile/logcmd", 
          getenv ("SC2LOGDIR"),name);
  strcpy(myInfo->logFile,myInfo->logfileName);
  SdpPutString(SC2DALOGFILE,myInfo->logFile, status);
  
  strcat(myInfo->logfileName,myInfo->Date);
  utils_fclose(&(myInfo->fpLog));
  if((myInfo->fpLog = fopen64(myInfo->logfileName, "a")) == NULL)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"Error- sc2dalib_initgetArrayID: failed to open %s",
           myInfo->logfileName );
    return;
  }
  // myFpLog is a global so all routines can write in the log file
  myFpLog = myInfo->fpLog;

  myInfo->logfileFlag=1; 

  jitDebug(2, "sc2dalib_initgetarrayId: the sc2dalib_logfile is %s\n",
           myInfo->logfileName);

 if (con->dataform ==MCE_BINARY_FORM)
     sprintf(dataform,"BINARY");
  else
     sprintf(dataform,"TEXT");

  fprintf(myInfo->fpLog,"\n<%s> CMD from sc2dalib__Init, <%s>", 
          dateTime,myInfo->initFile);
  fprintf(myInfo->fpLog,"\n<%s> DATA_FORMAT <%s>\n", mcexml,dataform);
  fprintf(myInfo->fpLog,"====> DASVERSION sourced= <%s> <====\n", getenv ( "SC2VERSION" ));
}



/**
 * \fn void sc2dalib_isFile(const char *fname,StatusType *status)
 *
 * \brief functaion
 *  check if file name is a existing regular file 
 *
 * \param fname   char pointer for the filename to be checked
 * \param status  StatusType.  given and return
 * 
 */
/*+ sc2dalib_isFile 
*/
void sc2dalib_isFile
(
const char *fname,
StatusType *status
) 
{
  struct stat sbuf;

  if (*status != STATUS__OK) return;

  //check if the file exists (=0, not exist=-1)
  if (lstat(fname, &sbuf) !=0) 
  {
    *status=DITS__APP_ERROR;
  }
  else
  { 
      *status=STATUS__OK;
  }
}

// =======sc2dalib_m*******
//====================//
/**
 * \fn void sc2dalib_mcecmdInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *   dasCmdInfo_t  *mymceCmd, char *dateTime, StatusType *status)
 *
 * \brief function
 *  get args for MCECMD action and open files
 *
 * \param con       SDSU context structure
 * \param myInfo   dasInfo structure pointer
 * \param mymceCmd dasCmdInfo_t pointer
 * \param dateTime  string pointer to dateTime string
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_mcecmdInit
*/
void sc2dalib_mcecmdInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *mymceCmd,
char                  *dateTime,
StatusType            *status
)
{
  long               in_sequence;
  long               range[]={0,2};
  DitsArgType          argId;

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  // Do not continue if MCE-xml is not read or in a data taking mode 
  if( myInfo->doneReadxml==0)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_mcecmdInit: MCE-xml is not read" ); 
    return;
  }
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  myInfo->actionFlag=MCECMDACTION;

  argId = DitsGetArgument();
  jitArgGetS(argId, "MCE_CMD", 1, NULL, "wb cc sram1_strt 15",
	     0, FILE_LEN, mymceCmd->mceCmd, NULL, status);
  jitArgGetS(argId, "CMDREP_FILE", 2, NULL, "cmdrep.txt",
	     0, FILE_LEN, myInfo->cmdrepFile, NULL, status);
  jitArgGetI(argId, "SVFILE_FLAG", 3, range, 0, 0, 
             &myInfo->filesvFlag, status);
  //sv default 1?
  if ( !StatusOkP(status) )
  {
    ErsOut(0,status, 
       "sc2dalib_mcecmdInit: Variables and defaults not good.");
    myInfo->filesvFlag=1;
  }

  // use filesvFlag ==2 to send cmd while in_seq
  if( myInfo->filesvFlag !=2)
  {
    if(in_sequence != DA_MCE_NONE)
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status,"sc2dalib_mcecmdInit: %s has not completed",
              seqStatus[myInfo->actionFlag]);
      return;
    }
  } 
  utils_open_files(myInfo,NODATAFILE,NOBATCHFILE,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_mcecmdInit: utils_open_files failed.");
    return;
  }  
  fprintf(myInfo->fpLog,"\n<%s> CMD for sc2dalib__Mcecmd <%s>\n",
          dateTime,mymceCmd->mceCmd);

  /* if(myInfo->filesvFlag !=0)
     fprintf(myInfo->fpLog,"<%s>\n",myInfo->cmdrepFile); */

  con->testflag=0; // calculated checksum in command_mce    
  con->framechk=0;  
  con->pcidacount=0;    
  con->datacount=0;
  con->process.exit=0;
  SdpPuti("IN_SEQUENCE",DA_MCE_NONE,status);
  con->process.whereabout=Dits_GetSeq0;
  con->process.framesetup= DA_MCE_NONE;
}



/**
 * \fn void sc2dalib_mceonflycmdInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *   dasCmdInfo_t  *mymceCmd, char *dateTime, StatusType *status)
 *
 * \brief function
 *  get args for MCECMD action and open files
 *
 * \param con       SDSU context structure
 * \param myInfo   dasInfo structure pointer
 * \param mymceCmd dasCmdInfo_t pointer
 * \param dateTime  string pointer to dateTime string
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_mceonflycmdInit
*/
void sc2dalib_mceonflycmdInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *mymceCmd,
char                  *dateTime,
StatusType            *status
)
{
  long               initialised;
  char               setfile[FILE_LEN];
  DitsArgType          argId;

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  // Do not continue if not initialized or in a data taking mode 
  SdpGeti("INITIALISED", &initialised, status);
  if( initialised==0)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_mcecmdInit: drama task is not initialised" ); 
    return;
  } 

  argId = DitsGetArgument();
  jitArgGetS(argId,"SETFILE",1,NULL,"sq2fboptimal.txt",0, FILE_LEN,
            setfile,NULL,status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_mcecmdInit: failed to get MCE_CMD"); 
    return;
  }
  sc2dalib_cnvtbatch2Cmd(setfile,mymceCmd->mceCmd, status);
}


/**
 * \fn void sc2dalib_mceerrnewRep(dasInfoStruct_t *myInfo,
 *  dasCmdInfo_t *mycmdInfo, StatusType *status)
 *
 * \brief functaion
 *  check reply status from MCE and report it for new protocol
 *
 * \param myInfo    dasInfoStruct_t  pointer
 * \param mycmdInfo  dasCmdInfo_t  pointer
 * \param status     StatusType.  given and return(STATUS__OK)
 *
 * report the Hardware error
 */
/*+ sc2dalib_mceerrnewRep - 
*/
void sc2dalib_mceerrnewRep
(
dasInfoStruct_t *myInfo,
dasCmdInfo_t    *mycmdInfo,
StatusType      *status
)
{
  char cmd[10]="";
  int  i, errNo,carderrNo;
  char localmsg[FILE_LEN];

  // message print and save to logfile, 
  // flag=0, use MsgOut, 
  // flag=1, use ErsOut, 
  // flag=2 use ErsRep
 
  strncpy(cmd,mycmdInfo->mceCmd,2);

  // *status != OK from entry, *status not changed after call
  // utils_msg as 2 is used
  utils_msg(myInfo,"sc2dalib_mceerrnewRep: executing command (%s)",
          mycmdInfo->mceCmd,USE_ERSREP,status);

  if (mycmdInfo->reply.status==MCE_NONREPLY)   
  {
    utils_msg(myInfo,
      "sc2dalib_mceerrnewRep: NO REPLY from driver when sdsu_command_mce ",
      "",USE_ERSOUT,status);
    utils_msg(myInfo,
      "sc2dalib_mceerrnewRep: mceTask exits?", "",USE_ERSREP,status);
  }
  else if (mycmdInfo->reply.status==DRV_MCEERR)   
  {
    if( mycmdInfo->reply.data[0]==PCI2MCE_TMO || 
        mycmdInfo->reply.data[0]==INT2MCE_TMO )
    { 
      utils_msg(myInfo,
        "sc2dalib_mceerrnewRep: %s",drvErrors[mycmdInfo->reply.data[0]],
        USE_ERSREP,status);
    }
    else 
    {
      utils_msg(myInfo,
        "sc2dalib_mceerrnewRep: %s",drvErrors[mycmdInfo->reply.data[0]],
        USE_ERSREP,status);
        sprintf(localmsg,"%#lX",mycmdInfo->reply.data[1]);
      utils_msg(myInfo,"sc2dalib_mceerrnewRep: %s",localmsg,
        USE_ERSREP,status);
    }
  }
  else                                
  {
    sprintf(localmsg,"%s:the error word(%08lx)",
              cmd,mycmdInfo->reply.status);
    utils_msg(myInfo,"sc2dalib_mceerrnewRep:%s",localmsg, 
        USE_ERSREP,status);  
    
    errNo= mycmdInfo->reply.status;
    for (i=0;i<MCE_CARD_NO;i++)
    {
      carderrNo=errNo & 0x07;
      sc2dalib_mceerreachCard(myInfo,cmd, mcecardStrt[i].name,carderrNo,status);
      errNo= errNo>>3;
    }
  }
}


/**
 * \fn void sc2dalib_mceerreachCard(dasInfoStruct_t *myInfo,
 *  char *mceCmd, char *cardName, int errNo, StatusType *status)
 *
 * \brief functaion
 *  check error for each card from MCE and report it for new protocol
 *
 * \param myInfo    dasInfoStruct_t  pointer
 * \param mceCmd     char  pointer
 * \param cardName   char  pointer
 * \param errNo      int   
 * \param status     StatusType pointer  
 *
 */
/*+ sc2dalib_mceerreachCard 
*/
void sc2dalib_mceerreachCard
(
dasInfoStruct_t *myInfo,
char            *mceCmd,
char            *cardName,
int             errNo,
StatusType      *status
)
{
  int  mceErr;
  char localmsg[FILE_LEN];

  // *status != OK from entry, *status not changed after call
  // utils_msg as 2 is used
  if ( errNo!=0)
  {
    if (errNo==1)
      mceErr=NO_CARD_ERR;
    else if (errNo==2)
      mceErr=CRC_ERR;
    else if (errNo==4)
      mceErr=READONLY_ERR;
    else 
    {
      sprintf(localmsg,"%s(%3s, %#X), wrong errno",mceCmd,cardName, errNo);
      utils_msg(myInfo,"sc2dalib_mceerreachCard: %s",localmsg,
         USE_ERSREP,status);
      return;
    }
    sprintf(localmsg,"%s(%3s, %#X), %s",mceCmd,cardName, errNo,mceErrors[mceErr]);
    utils_msg(myInfo,"sc2dalib_mceerreachCard: %s",localmsg,
         USE_ERSREP,status);
  }
}



/**
 * \fn void sc2dalib_mcestatusInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *   char *dateTime, StatusType *status)
 *
 * \brief function
 *  get args for MCESTATUS action and open file
 *
 * \param con       SDSU context structure
 * \param myInfo   dasInfo structure pointer
 * \param dateTime  string pointer to dateTime string
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_mcestatusInit
*/
void sc2dalib_mcestatusInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *dateTime,
StatusType            *status
)
{
  long            in_sequence,initialised;
  long            range[]={0,2}, flag;    
  DitsArgType     argId;

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  // Do not continue if not initialized or in a data taking mode 
  SdpGeti("INITIALISED", &initialised, status);
  if( initialised==0)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_mcestatusInit: drama task is not initialised" ); 
    return;
  } 
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status,"sc2dalib_mcestatusInit: %s has not completed",
       seqStatus[myInfo->actionFlag]);
    return;
  }
  myInfo->actionFlag=MCESTATUSACTION;

  argId = DitsGetArgument();
  jitArgGetS(argId, "STATUS_FILE", 1, NULL, "mcestatus.txt",
	     0, FILE_LEN, myInfo->cmdrepFile, NULL, status);
  jitArgGetI(argId, "BATCH_FLAG", 2, range, 0, 0, &flag, status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_mcestatusInit: variables not good.");
    return;
  }       
  myInfo->batchFlag=(int)flag;
  // always save ( appending ) the result
  myInfo->filesvFlag=1;

  utils_open_files(myInfo,NODATAFILE,NOBATCHFILE,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_mcestatusInit: utils_open_files failed");
    return;
  }  
  fprintf(myInfo->fpLog,"\n<%s> CMD for sc2dalib__MceStatus ",dateTime);
  fprintf(myInfo->fpLog,"<%s>\n",myInfo->cmdrepFile);
  fprintf(myInfo->fpMcecmd,"<HEADER>\n");
  SdpPuti("IN_SEQUENCE",DA_MCE_NONE,status);
  con->process.framesetup= DA_MCE_NONE;
}


// =======sc2dalib_p*******
//====================//

/**
 * \fn void sc2dalib_pcicmdInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  char *pciCmd, char *memType, long *memAddr, long *memValue,
 *  char  *dateTime, StatusType *status)
 *
 * \brief function
 *  get args for PCICMD action 
 *
 * \param con       SDSU context structure
 * \param myInfo   dasInfo structure pointer
 * \param pciCmd    pointer for pciCmd string
 * \param memType   pointer for memType string
 * \param memAddr   long pointer for  memAddr
 * \param memValue  long pointer for  memValue
 * \param dateTime  pointer to dateTime string
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_pcicmdInit
*/
void sc2dalib_pcicmdInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *pciCmd,
char                  *memType,
long                  *memAddr,
long                  *memValue,
char                  *dateTime,
StatusType            *status
)
{
  long        in_sequence,initialised;
  long        range[]={0,2000};
  long        rangeVal[]={0,2000000};
  DitsArgType argId;
  
  if (!StatusOkP(status)) return;
  
  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  // Do not continue if not initialized or in a data taking mode 
  SdpGeti("INITIALISED", &initialised, status);
/********
  if (myInfo->engFlag==RTSC_MODE)
  {
    if( initialised==0)
    {
      *status = DITS__APP_ERROR;
      ErsRep(0,status,"sc2dalib_pcicmdInit: drama task is not initialised" ); 
      return;
    }
  } 
*********/
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status, 
      "sc2dalib_pcicmdInit: %s has not completed ",seqStatus[myInfo->actionFlag]);
    return;
  } 

  myInfo->actionFlag=PCICMDACTION;

  argId = DitsGetArgument();
  jitArgGetS(argId, "PCI_CMD", 1, NULL, "READ", 0, CMD_LEN,
	     pciCmd, NULL, status );
  jitArgGetS(argId, "PCI_MEMTYPE", 2, NULL, "y", 0, 10,
	     memType, NULL, status );
  jitArgGetI( argId, "PCI_MEMADDR", 3, range, 0, 0, memAddr, status);
  jitArgGetI( argId, "PCI_MEMVALUE", 4, rangeVal, 0, 0, memValue, status);
  if ( !StatusOkP(status) )
  {
    ErsOut(0,status, "sc2dalib__Pcicmd: Variables or defaults not good.");
    return;
  }

  // pci cmd is still included in the counting 
  myInfo->glbCount++;
  fprintf(myInfo->fpLog,"\n<%s> CMD(%ld) for sc2dalib__PCIcmd <%s %s %ld %ld>\n",
          dateTime,myInfo->glbCount,pciCmd,memType,*memAddr,*memValue );
  fflush(myInfo->fpLog);
}
 


/**
 * \fn void sc2dalib_pciblkInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  char *pciCmd, char *memType, long *startAddr, long *blkSize,
 *  char  *dateTime, StatusType *status)
 *
 * \brief function
 *  get args for PCIBLK action 
 *
 * \param con       SDSU context structure
 * \param myInfo   dasInfo structure pointer
 * \param pciCmd    pointer for pciCmd string
 * \param memType   pointer for memType string
 * \param startAddr long pointer for  startAddr
 * \param blkSize   long pointer for  blkSize
 * \param dateTime  pointer to dateTime string
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_pciblkInit
*/
void sc2dalib_pciblkInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *pciCmd,
char                  *memType,
long                  *startAddr,
long                  *blkSize,
char                  *dateTime,
StatusType            *status
)
{
  long        in_sequence;
  long        range[]={0,200000};
  DitsArgType argId;
  
  if (!StatusOkP(status)) return;
  
  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status, 
      "sc2dalib_pciblkInit: %s has not completed ",seqStatus[myInfo->actionFlag]);
    return;
  } 

  myInfo->actionFlag=PCICMDACTION;
  // make sure don't save to fpMcecmd
  myInfo->filesvFlag=0;

  argId = DitsGetArgument();
  jitArgGetS(argId, "PCI_MEMTYPE", 1, NULL, "y", 0, 10,
	     memType, NULL, status);
  jitArgGetI( argId, "PCI_STARTADDR", 2, range, 0, 0, startAddr, status);
  jitArgGetI( argId, "PCI_BLKSIZE", 3, range, 0, 0, blkSize, status );
  if ( !StatusOkP(status) )
  {
    ErsOut(0,status, "sc2dalib__Pciblk: failed to get variables");
    return;
  }

  // pci cmd is still included in the counting 
  myInfo->glbCount++;
  fprintf(myInfo->fpLog,"\n<%s> CMD(%ld) for sc2dalib__PCIblk <%s %s %ld %ld>\n",
          dateTime,myInfo->glbCount,pciCmd,memType,*startAddr,*blkSize );

 jitDebug(2,"sc2dalib__Pciblk: memType(%s), startAddr(%ld %#0X>, blkSize(%ld)\n",
             memType,*startAddr,*startAddr,*blkSize);     
  fflush(myInfo->fpLog);
}
 


/**
 * \fn void sc2dalib_pcierrRep(PCI_CMD *pciptr, char *cmd,StatusType *status)
 *
 * \brief function
 *  check reply status from PCI and report it
 *
 * \param  pciptr  PCI_CMD structure pointer,
 * \param  cmd     char pointer for PCI CMD
 * \param  status  StatusType.  given and return
 */

/*+ sc2dalib_pcierrRep -
*/
void sc2dalib_pcierrRep
(
PCI_CMD *pciptr,
char * cmd,
StatusType *status
)
{
  // pass back the reply from command_pci
  // pciptr->arg1=pci_reply.status;
  // pciptr->arg2=pci_reply.data;
  // *ststua!=STATUS__OK when called
  // these are inside driver
  if( pciptr->arg1==PCI_ERR )
  {
     ErsRep(0,status,"sc2dalib_pcierrRep:cmd <%s>- %s",
         cmd, drvErrors[pciptr->arg2]);
  }
  else if( (pciptr->arg1 &0x00FFFFFF)==DRV_PCITMO)
  {
    ErsRep(0,status, "sc2pciErrRep:cmd <%s> timeout wait for reply",
      cmd);
  }
  // these error msg from PCI
  else if( (pciptr->arg2 &0x00FFFFFF)==PCI_NAL)
  { 
    ErsRep(0,status,"sc2dalib_pcierrRep:cmd <%s>- No Application Loaded",
          cmd);
  }
  else if( (pciptr->arg2 &0x00FFFFFF)==PCI_CNE)
  { 
    ErsRep(0,status,"sc2dalib_pcierrRep:cmd <%s> Command Name Error",
          cmd);
  }
  else if( (pciptr->arg2 &0x00FFFFFF)==PCI_MTE)
  { 
    ErsRep(0,status,
     "sc2dalib_pcierrRep:cmd <%s> Memory Type Error - memory type not valid",
      cmd);
  }
  else
  { 
    ErsRep (0,status,
     "something is Wrong for cmd <%s> with reply.status <%#lx> arg<%#lx>",
      cmd, pciptr->arg1,pciptr->arg2);
  }
}




/**
 * \fn void sc2dalib_pcistatusRep(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo,
 *  StatusType *status)
 *
 * \brief function
 *  read from PCI and report some status
 *
 * \param con     SDSU_CONTEXT  structure pointer
 * \param myInfo dasInfoStruct_t structure pointer
 * \param status  StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
/*+ sc2dalib_pcistatusRep - 
*/
void sc2dalib_pcistatusRep
(
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
StatusType      *status
)
{
  int             i;
  char            memTypeInt='X';
  long            memAddr[5],memValue[5]; 
  PCI_CMD  pci_cmd;
  char     localmsg[FILE_LEN];

  *status=STATUS__OK;

  for (i=0;i<3;i++)
  {  
    memValue[i]=0; memAddr[i]=i;
    *status=
      sdsu_set_pcicmd(&pci_cmd,"READ", memTypeInt,memAddr[i],&memValue[i]);
    if( *status !=SDSU_OK)
    { 
      ErsRep (0,status, "pcistatusRep: Error-sdsu_set_pcicmd failed");  
      *status=DITS__APP_ERROR;
      return;
    }
    if ((*status=sdsu_command_pci(con,&pci_cmd)) !=SDSU_OK)
    {
      ErsRep (0, status, 
              "sc2dalib_pcistatusRep: Error-sdsu_command_pci failed"); 
      *status=DITS__APP_ERROR;
      return;
    }
    memValue[i]=pci_cmd.arg2;
  } 
  sprintf(localmsg,
     "pcistatusRep:PCI(X:0)=%#lX (X:2)=%#lx PCI FrameCount(X:1)=%#lX(%ld)", 
      memValue[0], memValue[2], memValue[1], memValue[1]);
  utils_msg(myInfo,localmsg,"",USE_MSGOUT,status);
  return;
}
  

/**
 * \fn void  sc2dalib_pcisendCmd(SDSU_CONTEXT *con,char *cmd,
 *  dasInfoStruct_t *myInfo, long  memAddr, long  memValue, 
 *  char *memType, StatusType *status)
 *
 * \brief function
 *  send a command to PCI DSP
 *
 * \param con      SDSU_CONTEXT  structure pointer
 * \param cmd   char pointer for pcicmd
 * \param myInfo  dasInfoStruct_t structure pointer
 * \param memAddr  long memory addresss
 * \param memValue long: memory value
 * \param memType  char pointer for memory Type
 * \param *status  StatusType.  given and return
 *
 */
/*+   sc2dalib_pcisendCmd
*/
void  sc2dalib_pcisendCmd
(
SDSU_CONTEXT    *con,
char            *cmd,
dasInfoStruct_t *myInfo,
long            memAddr,
long            memValue,
char            *memType,
StatusType      *status
)
{   	        
  long            repeatLoop,i;
  char            memTypeInt;
  static PCI_CMD  pcicmdStr;
 
  if ( !StatusOkP(status) ) return;
  
  if(strcmp(memType,"X")==0 || strcmp(memType,"x")==0 )
    memTypeInt='X'; 
  else if(strcmp(memType,"Y")==0 || strcmp(memType,"y")==0) 
    memTypeInt='Y';
  else if(strcmp(memType,"P")==0 || strcmp(memType,"p")==0)
    memTypeInt='P'; 
  else 
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_pcisendCmd: -unknown memory Type(P:X:Y?)"); 
    return;
  }
 
  // allow PCI cmd to execute even if doneInit!=1  
  // the reply from command_pci is passed back to pci_cmd
  // pci_cmd->arg1=pci_reply.status;
  // pci_cmd->arg2=pci_reply.data;
  if( strcmp("WRITE",cmd) == 0 || strcmp("READ",cmd) == 0  || 
      strcmp("START",cmd) == 0 || strcmp("STOP",cmd) == 0  || 
      strcmp("RESET",cmd) == 0 || strcmp("RESETMCE",cmd) == 0 ||
      strcmp("REPEAT",cmd) == 0 )
  {
    *status=sdsu_set_pcicmd(&pcicmdStr,cmd, memTypeInt,memAddr,&memValue);   
    if( *status !=SDSU_OK )
    { 
      *status=DITS__APP_ERROR;
      ErsRep(0,status, "sc2dalib_pcisendCmd: sdsu_set_pcicmd failed");   
      return;
    }  
    if ( strcmp("REPEAT",cmd) == 0)
    {             
      repeatLoop=memValue;
      jitDebug(2,"sc2dalib_pcisendCmd: Start the repeat of writing memValue");      
      for (i=0;i<repeatLoop;i++)
      {
        sdsu_set_pciwritecmd(&pcicmdStr,0);
        *status=sdsu_command_pci(con,&pcicmdStr);
        if( *status !=SDSU_OK )
        {  
          *status=DITS__APP_ERROR;         
          sc2dalib_pcierrRep(&pcicmdStr,cmd,status);
          ErsRep(0,status, "sc2dalib_pcisendCmd: sdsu_command_pci failed"); 
          break;
        }
      }
      if (*status ==SDSU_OK)
      {
        jitDebug(2, "sc2dalib_pcisendCmd: OK! stop at %ld repeat CMD\n",i);
      }
      else
        return;
    }
    else
    {      
      *status=sdsu_command_pci(con,&pcicmdStr);
      if( *status !=SDSU_OK )
      {
        *status=DITS__APP_ERROR;         
        sc2dalib_pcierrRep(&pcicmdStr,cmd,status);
        ErsRep (0, status, "sc2dalib_pcisendCmd: sdsu_command_pci failed"); 
        if ( strcmp("RESET",cmd) == 0)
        {
          ErsRep(0,status,
         "interrupts over PCI or pci write to host memory (bus slave): failed ");
          ErsRep(0,status," switch off/on the PC and repeat this test!");
        }
        return;
      }
      if ( strcmp("WRITE",cmd) == 0)
      {  
        jitDebug(2,"sc2dalib_pcisendCmd: %s Value(%#lX) to %s MemAddr(%#lX)\n",
                  cmd,memValue,memType,memAddr);    
      }
      else if ( strcmp("READ",cmd) == 0)
      {   
        MsgOut(status,"sc2dalib_pcisendCmd: %s Value(%#lX: %ld) from %s MemAddr(%#lX)",
               cmd,pcicmdStr.arg2,pcicmdStr.arg2,memType,memAddr); 
      }
      else if ( strcmp("START",cmd) == 0)
      {  
        MsgOut(status,"sc2dalib_pcisendCmd: %s (application No=%ld) is Okay!",
               cmd,memValue);
        //application 0
        if(memValue==0)
        { 
          MsgOut(status,"====== Health status(%#lX)===========",
                 pcicmdStr.arg2);

          for (i=0;i<HEALTH_CHK_NO;i++)
            sc2dalib_healthStatus(i,pcicmdStr.arg2,status);

          MsgOut(status,"====== Health status  end ============");
        }
        else if(memValue==1)
        { 
          MsgOut(status,
           "===== Application %ld now running......        ======",memValue);  
          MsgOut(status,
           "===== PCI card echos back any word received on ======"); 
          MsgOut(status,
           "===== fibre and saves in Y memory (from Y:0).  ======");
        }
        else
        { 
          MsgOut(status,
            "======no associated application number(%ld)===============",
              memValue);
          MsgOut(status,
             "======(downloaded application is running) ================");
        }
      }
      else 
      {
        jitDebug(2,"sc2dalib_pcisendCmd: %s Okay!\n", cmd);    
      }
    }
  }
  else
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status, "sc2dalib_pcisendCmd: unknow PCI command");
    return;
  }
}


/**
 * \fn void  sc2dalib_pcisendblkCmd(SDSU_CONTEXT *con,char *cmd,
 *  dasInfoStruct_t *myInfo, long startAddr, long  blkSize, 
 *  char *memType, StatusType *status)
 *
 * \brief function
 *  read blk memory from PCI DSP
 *
 * \param con      SDSU_CONTEXT  structure pointer
 * \param cmd   char pointer for pcicmd
 * \param myInfo  dasInfoStruct_t structure pointer
 * \param startAddr  long start memory addresss
 * \param blkSize   long: block size ( 32bits words)
 * \param memType  char pointer for memory Type
 * \param *status  StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 */
/*+   sc2dalib_pcisendblkCmd
*/
void  sc2dalib_pcisendblkCmd
(
SDSU_CONTEXT    *con,
char            *cmd,
dasInfoStruct_t *myInfo,
long            startAddr,
long            blkSize,
char            *memType,
StatusType      *status
)
{
  int             i, multiNo,remainNo;
  long            memAddr=0, *retVal, tmpVal=0;
  char            memTypeInt;
  char            remainData[FILE_LEN], remainData1[FILE_LEN]="";
  char            tmpData[2][FILE_LEN];
  char            localmsg[FILE_LEN];
  static PCI_CMD  pcicmdStr;
 
  if ( !StatusOkP(status) ) return;
  
  if(strcmp(memType,"X")==0 || strcmp(memType,"x")==0 )
    memTypeInt='X'; 
  else if(strcmp(memType,"Y")==0 || strcmp(memType,"y")==0) 
    memTypeInt='Y';
  else if(strcmp(memType,"P")==0 || strcmp(memType,"p")==0)
    memTypeInt='P'; 
  else 
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_pcisendblkCmd: -unknown memory Type(P:X:Y?)"); 
    return;
  }
  // allocate  buffer for the return val
  if( (retVal=(long*)calloc(blkSize, sizeof(long))) == NULL)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_pcisendblkCmd: -failed to allocate buffer for retVals"); 
    return;
  }
  // allow PCI cmd to execute even if doneInit!=1  
  // the reply from command_pci is passed back to pci_cmd
  // pci_cmd->arg1=pci_reply.status;
  // pci_cmd->arg2=pci_reply.data;
  for (i=0;i<blkSize;i++)
  { 
     memAddr=startAddr+i;
    *status=sdsu_set_pcicmd(&pcicmdStr,cmd, memTypeInt,memAddr,&tmpVal);   
    if( *status !=SDSU_OK )
    { 
      *status=DITS__APP_ERROR;
      ErsRep(0,status, "sc2dFapcisendblkCmd: sdsu_set_pcicmd failed");   
      return;
    }  
    *status=sdsu_command_pci(con,&pcicmdStr);
    if( *status !=SDSU_OK )
    {
      *status=DITS__APP_ERROR;         
      sc2dalib_pcierrRep(&pcicmdStr,cmd,status);
      ErsRep (0, status, "sc2dalib_pcisendblkCmd: sdsu_command_pci failed"); 
      return;
    }
    *(retVal+i)=(pcicmdStr.arg2 & 0xFFFFFF);
  }
  
  sprintf(localmsg,
   " %ld Val read from %s:%#lX (for fiber-related CMD: only (16bit) in little Endian LSB MSB)",
         blkSize, memType,startAddr); 
  utils_msg(myInfo,localmsg,"",USE_MSGOUT,status);
  multiNo  =blkSize/16;
  remainNo =blkSize- multiNo*16;
  for (i=0;i<multiNo;i++)  
  {
    sprintf(tmpData[0],"%06lX %06lX %06lX %06lX %06lX %06lX %06lX %06lX ",
         retVal[i*16+0],retVal[i*16+1],retVal[i*16+2],retVal[i*16+3],
         retVal[i*16+4],retVal[i*16+5],retVal[i*16+6],retVal[i*16+7]);
 
    sprintf(tmpData[1],"%06lX %06lX %06lX %06lX %06lX %06lX %06lX %06lX",
         retVal[i*16+8],retVal[i*16+9],retVal[i*16+10],retVal[i*16+11],
         retVal[i*16+12],retVal[i*16+13],retVal[i*16+14],retVal[i*16+15]);
 
    strcat(tmpData[0],tmpData[1]);
    utils_msg(myInfo,tmpData[0],"",USE_MSGOUT,status);
  }
  if(remainNo !=0)
  {
    for (i=0;i<remainNo;i++)  
    {
      sprintf(remainData,"%s %06lX", remainData1,retVal[multiNo*16+i]);
      strcpy(remainData1,remainData);
    } 
    utils_msg(myInfo,remainData,"",USE_MSGOUT,status);
  }
}


/**
 * \fn void sc2dalib_pixelmonInit(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *   dasCmdInfo_t *myCmd, struct mcexml_struct  *mceinxPtr,
 *   char *dateTime, ARRAYSET *pixelset, StatusType *status)
 *
 * \brief function
 *  get args for PIXELMONITOR action and open file
 *
 * \param con       SDSU context structure pointer 
 * \param myInfo    dasInfoStruct_t pointer
 * \param myCmd      dasCmdInfo_t pointer
 * \param mceinxPtr   mcexml_struct pointer for parameter lookup table
 * \param dateTime   dateTime string pointer         
 * \param pixelset   ARRAYSET structure pointer for pixelMon parameters
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_pixelmonInit
*/
void sc2dalib_pixelmonInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
char                  *dateTime,
ARRAYSET              *pixelset,
StatusType            *status
)
{
  long         range[]={0,1};
  long         in_sequence,configured;
  char          tmp[FILE_LEN];
  DitsArgType  argId;

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  SdpGeti("CONFIGURED", &configured, status);
  if (configured==0 )
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_pixelmonInit: the DA is not configured" ); 
    return;
  }
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_pixelmonInit: %s has not completed",
             seqStatus[myInfo->actionFlag]);
    return;
  } 
  myInfo->actionFlag=PIXELMONACTION;

  argId = DitsGetArgument();

  // filesvFlag=1 save cmdreply 
  // filesvFlag=2 don't save cmdreply, but save data
  // for servo, we only save data
  myInfo->filesvFlag=2;

  jitArgGetS(argId, "DATA_FILE", 1, NULL, "data.txt", 0, FILE_LEN,
	     myInfo->dataFile, NULL, status);
  jitArgGetS(argId, "SETUPBATCH_FILE", 2, NULL, "pixelsetup.txt", 0, FILE_LEN,
	     myInfo->batchFile, NULL, status);
  jitArgGetI(argId, "DATASPLIT_FLAG", 3, range, 0, 0,
	     &myInfo->datasplitFlag, status );
  if ( !StatusOkP(status) )
  {
    ErsOut(0,status, "sc2dalib_pixelmonInit: failed to get Variables");
    return;
  }

  fprintf(myInfo->fpLog,"\n<%s> CMD from sc2dalib__pixelmon <%s>\n",
          dateTime,myInfo->batchFile);
  fflush(myInfo->fpLog);

  jitArgGetS (argId, "SQ2OPT_FILE", 4,NULL, "sq2fboptimal.txt", 0, FILE_LEN,
              tmp,  NULL,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_pixelmonInit: failed to get SQ2OPT_FILE ");
    return;
  }
  // read sq2fboptimal in and store in global sq2fboptInit
  sc2dalib_readOPT(myInfo,pixelset->sq2fdbkOpt,COL_NUM,tmp,status);
  if ( !StatusOkP(status) )
    return;
  MsgOut(status,"c2dalib_pixelmonInit: read %s",tmp);
  
  sc2dalibsetup_servoreadsetupWrap(myInfo,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_pixelmonInit: sc2dalibsetup_servoreadsetupWrap failed"); 
    return;
  }

  // open datafile "a" for stripchart to look, 
  // open otherfile "r" for setup
  utils_open_files(myInfo,DATFILEAPPEND,OTHERFILE,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_pixelmonInit: utils_open_files failed"); 
    return;
  }
  // use this to read pixel monitor setting
  sc2dalibsetup_servoreadSetup(myInfo,pixelset,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0, status,"sc2dalib_pixelmonInit:sc2dalibsetup_servoreadSetup failed ");
    return;
  }
  if(pixelset->totalRC2use==1)  // only one card
  {
    sprintf(myInfo->goCmd, "GO rc%d ret_dat 1",pixelset->whichRC[0]);
  }
  else  if(pixelset->totalRC2use==4)          /// four cards
  {
    sprintf(myInfo->goCmd, "GO rcs ret_dat 1");
  }
  else
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,
          "sc2dalib_pixelmonInit: there are %d RC cards used, not four",
          pixelset->totalRC2use);
    ErsRep(0,status,
          "                I need more information about how to");
    ErsRep(0,status,
          "                send command to MCE");
    return;     
  }
  strcpy(myCmd->mceCmd,myInfo->goCmd);
  sc2dalib_getcmdBuf(myInfo,myCmd,mceinxPtr,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_pixelmonInit: sc2dalib_getCmdBuf failed"); 
    return;
  }
  myInfo->cmdFlag=0;  // noflycmd
}


// =======sc2dalib_r*******
//====================//
/**
 * \fn void sc2dalib_readpixelMask(dasInfoStruct_t *myInfo, 
 *  int *pixelMask, SdsIdType maskId, StatusType *status)
 *
 * \brief function
 *  get args for PIXEL_MASK action and fill pixelMask array 
 *
 * \param myInfo      dasInfo structure pointer
 * \param pixelMask   int pointer for pixel mask array
 * \param maskId      SdsIdType for he mask
 * \param status       StatusType     
 *
 */
/*+ sc2dalib_readpixelMask
*/
void sc2dalib_readpixelMask
(
dasInfoStruct_t   *myInfo,    
int               *pixelMask,
SdsIdType         maskId,
StatusType        *status
)
{
  SdsIdType pixelmaskId; //id for PIXELMASK structure 
  SdsIdType rowId;        //  id for row structure 
  SdsIdType id = 0;
  SdsCodeType code;
  
  char  name[FILE_LEN];
  char  delimiters[]= " ",*token;
  char  cardName[40], rc[50];

  long  ndims;            /* number of dimensions in row */
  long  rowNo;
  int   j,i, pixel,card;
  unsigned long dims[7];       
  unsigned long indims[7];    


  // initID is top SDS, the rest is sub sds. get sub-sds 
  SdsFind ( maskId, "PIXELMASK", &pixelmaskId, status );

  dims[0] = dims[1] = 0;

   // from maskId, get row 
  SdsFind ( pixelmaskId, "row", &rowId, status );
  if (!StatusOkP(status))
  {
    ErsRep(0,status,"sc2dalib_readpixelMask: find row failed");
    SdsFreeId ( pixelmaskId, status ); 
    return;
  }
  SdsInfo ( rowId, name, &code, &ndims, dims, status );
  jitDebug(2,"rowId's name =%s ndims=%d, dims[0]=%d dims[1]=%d\n",
             name, ndims, dims[0],dims[1]);

  // we only have one dimension of row SDS, ie. ndims==1
  // we has dim[0]=40 row SDS, dims[0]=0
  // search each element of row, extract the values 
  //  mask=1: initialised during load task
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
        ArgGetString(id, cardName, 40, rc, status );
        if (StatusOkP(status))
        {
          token = strtok (rc, delimiters);
          pixelMask[pixel]=atoi(token);
          for ( i=1; i<8; i++)
          {
            pixel ++;
            token=strtok(NULL,delimiters);
            pixelMask[pixel]=atoi(token);
          }
          pixel ++;
        }
      }
    }    
    SdsFreeId ( id, status );
    if (!StatusOkP(status))
    {
      ErsRep(0,status, 
      "sc2dalib_readpixelMask: row(%d) item failed in pixel(%d)",
       j,pixel); 
      break;
    }
  }
  jitDebug (2,"sc2dalib_readpixelMask: get rowCount=%d\n",j);
  SdsFreeId ( rowId, status );
  SdsFreeId ( pixelmaskId, status );
}


// =======sc2dalib_s*******
//====================//

/**
 * \fn void sc2dalib_saveframeData(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  FILE *fp, char * byte, int datasize,StatusType *status)
 *
 * \brief functaion
 *  save frame data in TEXT or BINARY format
 *
 * \param con      SDSU context structure
 * \param myInfo  dasInfo structure pointer
 * \param byte     pointer to current frame data buffer 
 * \param datasize frame data size in 32bits
 * \param fp       FILE pointer
 * \param status   StatusType, given and retutrn: 
 *
 */
/*+ sc2dalib_saveframeData 
*/
void sc2dalib_saveframeData
(
SDSU_CONTEXT    *con,       
dasInfoStruct_t *myInfo,
FILE            *fp,    
char            *byte,
int             datasize,
StatusType      *status
)
{
  int     *longword;
  int     j, wordsize,l,toploop,midloop, remain,headNo;

  if (*status != STATUS__OK) return;

  midloop=32;
  headNo=FRAMEHEADER_NUM;
  wordsize=datasize-headNo;
  toploop=(wordsize)/midloop;
  remain =(wordsize)%midloop;
  longword = (int *)byte;
   
  jitDebug(8,"sc2dalib_saveframeData: FrameCount in frame_status words =%d\n",
           *(longword+1));
  
  if (myInfo->dataFormat==MCE_TEXT_FORM)
  {  
    //request by UBC to have a single line for all header inf.
    for (j=0;j<headNo;j++)
      fprintf(fp,"%11d ", *(longword +j));
    fprintf(fp,"\n");
      longword +=headNo;

    for (j=0;j<toploop;j++)
    {
      for(l=0;l<midloop;l++)
        fprintf(fp,"%11d ",*(longword+ l+midloop*j) );
      fprintf(fp,"\n");
    }
    toploop =toploop*midloop;
    for(l=0;l<remain;l++)
      fprintf(fp,"%11d ", *(longword+ l+ toploop) );
    if(remain)
      fprintf(fp,"\n");
  }
  else if (myInfo->dataFormat==MCE_TEXT_FORM2)
  {
    for(j=0;j<wordsize;j++)  
      fprintf(fp,"%d \n", *(longword+j) );
    fprintf(fp,"\n" );
  }
  else
  {
    fwrite(byte, 1, datasize*4, fp);
  }
  return;
}


/**
 * \fn void sc2dalib_savepixelData(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *  char *byte, int datasize, ARRAYSET *pixelset, StatusType *status)
 *
 * \brief functaion
 *  save pixel data for stripchart
 *
 * \param con      SDSU context structure
 * \param myInfo  dasInfo structure pointer
 * \param byte     pointer to current frame data buffer 
 * \param datasize frame data size in 32bits
 * \param pixelset  ARRAYSET    pointer
 * \param status   StatusType, given and retutrn: 
 *
 */
/*+ sc2dalib_savepixelData 
*/
void sc2dalib_savepixelData
(
SDSU_CONTEXT    *con,       
dasInfoStruct_t *myInfo,    
char            *byte,
int             datasize,
ARRAYSET        *pixelset,
StatusType      *status
)
{
  int     *pixelPtr,*longword;
  int      i, j;
  int      errorData,data,data_1;

  if (*status != STATUS__OK) return;

  // point to frame data as INT
  longword = (int *)byte;
  jitDebug(8,"sc2dalib_savepixelData: FrameCount in frame_status words =%d\n",
           *(longword+1));

  // skip the header
  pixelPtr=longword+FRAMEHEADER_NUM; 
  
  // save data in 
  // [row_0][col_0] ~ [row_0][col_M] [row_1][col_0] ~ [row_1][col_M]...
  // save data  if dataslitFlag=1 
  // [row_0][col_0] [row_0][col_0 error] ..[row_0][col_M] [row_0][col_M error] ...

  fprintf(myInfo->fpData,"%7d ",myInfo->trkNo);
  for (i=0;i<ROW_NUM;i++)
  {
    if( pixelset->rowMask[i] >0)
    {    
      for (j=0;j<COL_NUM;j++)
      {
        if( pixelset->colMask[j] !=0)
        {
          jitDebug(8,"row=%d, col=%d datasplitFlag=%d\n",i,j,myInfo->datasplitFlag);
          data = *(pixelPtr + j + i*COL_NUM);
          if (myInfo->datasplitFlag ==0)
            fprintf(myInfo->fpData,"%7d ",data );
          else
          {
             data_1 = (data >> 16) & 0xFFFF;
             errorData= data & 0x0000FFFF;
             fprintf(myInfo->fpData,"%7d %7d",data_1, errorData);
          }
        }
      }
    }
  }
  fprintf(myInfo->fpData,"\n");
  fflush(myInfo->fpData);
}


/**
 * \fn void sc2dalib_savemcerepData(dasInfoStruct_t *myInfo, FILE *fpout,
 *  int mxbuf, char * byte)
 *
 * \brief function
 *  save MCE reply data into a file, used for mcestatus
 *
 * \param myInfo    dasInfoStruct_t pointer
 * \param fpout  FILE pointer for store the reply 
 * \param mxbuf  int: the buffer size in 32bits
 * \param byte   pointer to the reply buffer
 *
 */
/*+ sc2dalib_savemcerepData 
*/
void sc2dalib_savemcerepData
(
dasInfoStruct_t  *myInfo,
FILE             *fpout,   //store the replied data  
int              mxbuf,
char             *byte    // pointer to reply Buffer 
)
{
  int   *longword;
  int   j,toploop;

  // we do not want the first two RBOK and CardId_paraId
  // and last chksum
  toploop=(mxbuf-3);
  longword = (int *)byte;

  for (j=0;j<toploop;j++)
    fprintf(fpout,"%08d ",*(longword+2 + j ) );
  fprintf(fpout,"\n");
}


/**
 * \fn void sc2dalib_sendCmd(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  dasCmdInfo_t *mycmdInfo, struct mcexml_struct *mceInxpt,
 * char *dateTime,StatusType *status)
 *
 * \brief function
 *  send a cmd to MCE
 *
 * \param con       SDSU context structure pointer 
 * \param myInfo    dasInfoStruct_t pointer
 * \param mycmdInfo  dasCmdInfo_t pointer
 * \param mceInxpt   mcexml_struct pointer for parameter lookup table
 * \param dateTime   dateTime string pointer         
 * \param status     StatusType.  given and return
 *
 * if  *status != STATUS__OK, report error.
 *
 * translate MCE command into command struture, send it to MCE
 * save the comamnd buffer and reply or error in logfile 
 *
 */
/*+ sc2dalib_sendCmd   
 */
void sc2dalib_sendCmd
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
dasCmdInfo_t          *mycmdInfo,
struct mcexml_struct  *mceInxpt,
char                  *dateTime,
StatusType            *status
)
{
  char         delimiters[]= " ",*token;;
  char  cp[FILE_LEN];

  if (!StatusOkP(status)) return;

  sc2dalib_getcmdBuf(myInfo,mycmdInfo,mceInxpt,status);
  if ( *status ==DITS__APP_ERROR)
  {     
    ErsRep (0, status,"sc2dalib_sendCmd: faile to call sc2dalib_getcmdBufi");
    return ;
  }
  
  // check for the dummy setting for external stuff
  // token => "command" before "space"
  strcpy(cp,mycmdInfo->mceCmd); 
  token = strtok (cp, delimiters);
  token = strtok (NULL,delimiters);
  //jitDebug(2,"sc2dalib_sendCmd: token (%s)\n",token);
  
  if (strcmp(token, EXTERNAL_CARD) ==0)
  {
    // the EXTERL one only for store setting in mceInxpt
    // don't send out
  }
  else
  {
    jitDebug(2,"sc2dalib_sendCmd:(%s)\n",mycmdInfo->mceCmd);
    mce_cmd(con,myInfo,mycmdInfo,dateTime,status);  
    if ( *status ==DITS__APP_ERROR)
    {     
      ErsRep (0, status,"sc2dalib_sendCmd: failed to call mce_cmd");
      return;
    }
    if(myInfo->debuglvl==DISP_GEN_MSG)     
    {
      MsgOut(status,
                  "sc2dalib_sendCmd: PacketSize from MCE=<%d>",
                  mycmdInfo->replySize); 
    }
  }
}


/**
 * \fn void sc2dalib_setseqInit(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceInxpt, char *dateTime,StatusType *status)
 *
 * \brief function
 *  get args for SETUP_SEQUENCE action,
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param mceInxpt   mcexml_struct pointer for parameter lookup table
 * \param dateTime  dateTime string pointer         
 * \param status    StatusType.  given and return
 *
 */
/*+ sc2dalib_setseqInit   
 */
void sc2dalib_setseqInit
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,
char                  *dateTime,
StatusType            *status
)
{
  char  mceCard[FILE_LEN];
  long  configured, in_sequence;
  SdsIdType    argId;
  long     frameRate; 

  if (*status != STATUS__OK) return;
 
  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  // Do not continue if not initialized, not configured or in data
  // takeing mode 
  SdpGeti("CONFIGURED", &configured, status);
  if(configured == 0 && myInfo->engFlag==RTSC_MODE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_setseqInit: drama task is not configured" );
    return;
  }
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence !=DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_setseqInit: %s has not completed",
           seqStatus[myInfo->actionFlag]);
    return;
  } 
  myInfo->actionFlag=SETSEQACTION;

  argId = (SdsIdType)DitsGetArgument();

  jitArgGetS (argId, "MCE_WHICHCARD", 40,NULL, " ",0, FILE_LEN, mceCard,
              NULL,status);
   if (!StatusOkP(status))
    return;

  //SdsList( argId, status);

  jitDebug(2,"sc2dalib_setseqInit: MCE_WHICHCARD= %s\n", mceCard);
  // basefile is passed from RTSC as yyyymmdd ( date )

  if ( strcmp(mceCard," ")==0)
  { 
    // it is for RTSC  myInfo->baseFile is set in configure
    // con->process.svfile is not updated
    if( myInfo->engFlag==RTSC_MODE)
    {
      sc2dalib_setseqInitRTSC(con,myInfo,mceInxpt,dateTime,mceCard,status);
      jitDebug(2,"sc2dalib_setseqInit: MCE_WHICHCARD =%s\n", mceCard);
      // pass baseFile in case svFlag>0
      strncpy(con->process.svfile,myInfo->baseFile,sizeof(con->process.svfile));
    }
    else
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,"sc2dalib_setseqInit: mceconfig-rtsc not set");
      return;
    }
  }
  else
  { 
    // it is for ENG myInfo->baseFile is set in configure
    // con->process.svfile is updated  from SETUP_SEQ arg
    if( myInfo->engFlag==ENG_MODE)
      sc2dalib_setseqInitENG(con,myInfo,mceCard,status);
    else
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,"sc2dalib_setseqInit: mceconfig-stare (scan/dream) not set");
      return;
    }
  }

  if (!StatusOkP(status))
    return;

  sprintf(myInfo->goCmd, "GO %s ret_dat 1",mceCard);
  sprintf(myInfo->stopCmd, "ST %s ret_dat",mceCard);

  SdpPutString("MCE_WHICHCARD",mceCard, status);

  sc2dalibsetup_whichRC(myInfo,mceCard,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_setseqInit: not recognised(%s) as MCE_WHICHCARD",
           mceCard);
    return;
  }

  mce_read_frame_rate(con,myInfo,mceInxpt,&frameRate,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status,"sc2dalib_setseqInit failed to call mce_read_frame_rate");
    return;

  }
  SdpPuti("FRAME_RATE", frameRate, status);

  fprintf(myInfo->fpLog,"\n<%s> CMD for sc2dalib__SetSeq \n",dateTime);
  /* fprintf(myInfo->fpLog,"FFRAME_FILENAME <%s> MCE_WHICHCARD <%s>\n",
     con->process.svfile, mceCard); */
}


/**
 * \fn void sc2dalib_setseqInitENG(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  char *whichCard,  StatusType *status)
 *
 * \brief function
 *  get args for SETUP_SEQUENCE action for ENG
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo    dasInfoStruct_t pointer
 * \param whichCard string pointer         
 * \param status    StatusType.  given and return
 *
 */
/*+ sc2dalib_setseqInitENG
 */
void sc2dalib_setseqInitENG
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
char                  *whichCard,
StatusType            *status
)
{
  long       svflag;
  long       range[] = { 0, 3};
  SdsIdType  argId;

  if (*status != STATUS__OK) return;

  argId = (SdsIdType)DitsGetArgument();

  jitArgGetS(argId, "FRAME_FILENAME", 1, NULL, "20007", 0, FILE_LEN,
	     con->process.svfile, NULL, status); 
  jitArgGetS(argId, "MCE_WHICHCARD", 2, NULL, "rc1", 0, FILE_LEN,
	     whichCard, NULL, status);
  // inform sc2dalib_dh task to sava data
  jitArgGetI(argId, "FILESAVE_FLAG", 3, range, 0, 0, &svflag, status);
  con->process.svflag=(int)svflag;
  jitArgGetI( argId, "UPDATE_SQ2FB", 4, range, 0, 0, &svflag, status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_setseqInitENG: Failed to get Variables");
    return;
  }
  myInfo->gotsq2fbparaFlag=(int)svflag;
  // pass on to DH task
  myInfo->parshmPtr->sq2fbparaFlag=(int)svflag;
}


/**
 * \fn void sc2dalib_setseqInitRTSC(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceInxpt, char *dateTime,char * whichCard,  
 *  StatusType *status)
 *
 * \brief function
 *  get args for SETUP_SEQUENCE action for RTSC
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param mceInxpt   mcexml_struct pointer for parameter lookup table
 * \param dateTime  dateTime string pointer         
 * \param whichCard string pointer         
 * \param status    StatusType.  given and return
 *
 */
/*+ sc2dalib_setseqInitRTSC
 */
void sc2dalib_setseqInitRTSC
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,
char                  *dateTime,
char                  *whichCard,
StatusType            *status
)
{
  SdsIdType  argId;
  long       group=0,heat, bbheat, seqcount, obsend, drcontrol;
  long       bb_inbeam, fts_inbeam, pol_inbeam,sq2fbflag;
  long       stairNum, stairHeight, stairWidth, stairStart;
  long       range[] = { -99999, 32766};  // from rtsDClient
  double     rangeD[]={0.0, 1},shutter;
  int        val_indx;
  char       source[MAX_INPUT_STR];
  char       inBeamArray[MAX_INBEAM];
  const char *source_values[] = {"SCIENCE", "REFERENCE", 0};
  const char *load_values[] = {"SKY", "LOAD2", "HALFHOT", "LINE", "DARK", "HOT",0};
  char       datamodeVal[]="rb rcs data_mode 1";
  int        sawtoothRampFlag;

  if (*status != STATUS__OK) return;

  // for RTSC, always use rcs?? no save ENG-DATA
  sprintf(whichCard, "rcs");

  argId = (SdsIdType)DitsGetArgument();

  /* GROUP */
  jitArgGetI( argId, "GROUP", 0, range, 0, 0, &group, status);
  jitDebug(2,"sc2dalib_setseqInitRTSC: group=%ld \n",group);

  /* TASKS  DAS doesn't need to read tasks */

  /* SOURCE */
  jitArgGetS( argId, "SOURCE", 0, source_values, "SCIENCE", GIT_M_ARG_UPPER, 
		  MAX_INPUT_STR, source, &val_indx, status);

  jitDebug(2,"sc2dalib_setseqInitRTSC: source=%s \n",source);

  /* LOAD */
  jitArgGetS (argId, "LOAD", 0,load_values, "SKY",GIT_M_ARG_UPPER, FILE_LEN, 
              myInfo->load, &val_indx,status);

  jitDebug(2,"sc2dalib_setseqInitRTSC: load=%s \n",myInfo->load);

  /* other arg */
  jitArgGetI( argId, "HEAT_CUR", 0, range, 0, 0, &heat, status);
  jitArgGetI( argId, "BB_TEMP", 0, range, 0, 0, &bbheat, status);
  range[0] = 0;
  jitArgGetI( argId, "SEQCOUNT", 0, range, 0, 0, &seqcount, status);
  jitArgGetI( argId, "DRCONTROL", 0, range, 0, 0, &drcontrol, status);
  jitArgGetD( argId, "SHUT_FRAC", 0, rangeD, 0.0, 0, &shutter, status);
  range[1] = 1;
  jitArgGetI( argId, "OBSEND", 0, range, 0, 0, &obsend, status);

  /* The stair handling information */
  range[1] = 65535;
  jitArgGetI( argId, "STAIR_NUM", 0, range, 0, 0, &stairNum, status);
  jitArgGetI( argId, "STAIR_WIDTH", 0, range, 0, 0, &stairWidth, status);
  jitArgGetI( argId, "STAIR_HEIGHT", 0, range, 0, 0, &stairHeight, status);
  jitArgGetI( argId, "STAIR_START", 0, range, 0, 0, &stairStart, status);
  myInfo->stairHeight = stairHeight;
  myInfo->stairWidth = stairWidth;
  myInfo->stairNum = stairNum;
  myInfo->stairStart = stairStart;

  /* Grab BB_INBEAM from the rtsClient, then call the fillInbeam routine to
     build up the entire inBeam array */
  range[1] = 7; /* there is more than one part of the Polarimeter that can be in the beam at any one time */
  jitArgGetI( argId, "BB_INBEAM", 0, range, 0, 0, &bb_inbeam, status);
  jitArgGetI( argId, "FTS_INBEAM", 0, range, 0, 0, &fts_inbeam, status);
  jitArgGetI( argId, "POL_INBEAM", 0, range, 0, 0, &pol_inbeam, status);
  sc2dalib_fillInBeam(bb_inbeam, fts_inbeam, pol_inbeam, shutter, inBeamArray, status);
  strcpy(myInfo->inbeam,inBeamArray);

  // add sq2fb flag for updating 
  jitArgGetI( argId, "UPDATE_SQ2FB", 0, range, 0, 0, &sq2fbflag, status );
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_setseqInitRTSC: failed to get UPDATE_SQ2FB");
    return;
  }
  myInfo->gotsq2fbparaFlag=(int)sq2fbflag;
  // pass on to DH task
  myInfo->parshmPtr->sq2fbparaFlag=(int)sq2fbflag;
 

  myInfo->pixelHeat=(int)heat;
  myInfo->bbHeat   =(int)bbheat;
  myInfo->seqcount   =(int)seqcount;
  myInfo->drcontrol   =(int)drcontrol;
  myInfo->shutterFraction=(float)shutter;
  jitDebug(2,"sc2dalib_setseqInitRTSC: HEAT_CUR(%d) BB_TEMP(%d) SHUT_FRAC(%4.2f)\n",
              myInfo->pixelHeat, myInfo->bbHeat, myInfo->shutterFraction);

  // Handle initialization of heater control
  heat_init_bias(con, myInfo, mceInxpt, &sawtoothRampFlag, status);
  if (!StatusOkP(status)) 
    return;
  // pass on to DH task
  myInfo->parshmPtr->sawtoothRampFlag = sawtoothRampFlag;
  myInfo->sawtoothRampFlag = sawtoothRampFlag;
  myInfo->parshmPtr->stairWidth = myInfo->stairWidth;

  // read the data_mode setting from the MCE
  mce_read(con,myInfo,mceInxpt, datamodeVal, &myInfo->datamode,1,status);
  if (!StatusOkP(status)) 
    {
      ErsRep(0,status, "sc2dalib_setseqInitRTSC: mce_read failed to read data_mode"); 
      return;
    }


  /* Set the value of OBSEND based on what we got from the rtsc */
  if(obsend)
    sc2headman_setobsend(1,status);
  else
    sc2headman_setobsend(0,status);

  // make sure it does not save ENG_DATA
  con->process.svflag=0; 
  jitDebug(2,"sc2dalib_setseqInitRTSC:OK\n");
}


/**
 * \fn void sc2dalib_fillInBeam(long bb_inbeam, long fts_inbeam, long pol_inbeam, double shutter 
 *    char * inBeamArray, StatusType *status)
 *
 * \brief function
 *  Create the INBEAM string for the FITS headers
 *
 * \param bb_inbeam    int         TRUE or FALSE for Blackbody in beam
 * \param fts_inbeam   int         TRUE or FALSE for FTS2 in beam
 * \param pol_inbeam   int         Bits for various pieces of POL2 in beam
 * \param shutter      double      0.0 closed 1.0 Open, any thing else is a little in beam
 * \param inBeamArray  char*       Where to write the string
 * \param status       StatusType  given and return
 *
 */
/*+ sc2dalib_fillInBeam   
 */
void sc2dalib_fillInBeam
(
long              bb_inbeam,         
long              fts_inbeam,         
long              pol_inbeam,
double            shutter,         
char             *inBeamArray,
StatusType       *status
)
{

  if ( *status != STATUS__OK ) return;

  strcpy(inBeamArray, "");

  /*  create a space separated list of INBEAM values
      from SCUBA2, POL2 and FTS2 */
  
  if(bb_inbeam)
    {
      strcat(inBeamArray, "blackbody ");
    }

  /* FTS2 */

  if(fts_inbeam)
    {
      strcat(inBeamArray, "fts2 ");
    }

  /* POL2 */

#define CALIBRATOR 1
#define ANALYSER   2
#define WAVEPLATE  4

  /* The standard configuration for POL-2 is the waveplate and analyser in the beam and the
     calibrator out.  In this case we just put POL in inbeam otherwise spell them out separately */

  if(((pol_inbeam & WAVEPLATE) == WAVEPLATE) && ((pol_inbeam & ANALYSER) == ANALYSER) &&
     ((pol_inbeam & CALIBRATOR) != CALIBRATOR))
    {
      strcat(inBeamArray, "pol ");
    }
  else
    {
      if((pol_inbeam & WAVEPLATE) == WAVEPLATE)
	{
	  strcat(inBeamArray, "pol2_wave ");
	}

      if((pol_inbeam & CALIBRATOR) == CALIBRATOR)
	{
	  strcat(inBeamArray, "pol2_cal ");
	}

      if((pol_inbeam & ANALYSER) == ANALYSER)
	{
	  strcat(inBeamArray, "pol2_ana ");
	}
    }

  /* Shutter */
  if(shutter != 1.0)
   {
      strcat(inBeamArray, "shutter ");
    }
}



/**
 * \fn void sc2dalib_seqInit(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  char *dateTime,int **lookupTable,StatusType *status)
 *
 * \brief function
 *  get args for SEQUENCE action,allocate shared memory,
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param dateTime  dateTime string pointer         
 * \param lookupTable address of int pointer         
 * \param status    StatusType.  given and return
 *
 */
/*+ sc2dalib_seqInit   
 */
void sc2dalib_seqInit
(
SDSU_CONTEXT          *con,         
dasInfoStruct_t       *myInfo,
char                  *dateTime,
int                   **lookupTable,
StatusType            *status
)
{
  long          setup,in_sequence;  	
  long          range[] = { 1, 0x7FFFFFFE};
  long          startSeq, endSeq,dWell;
  DitsArgType   argId;
  int           *table=NULL;
  char          obsMode[40];  

  if (!StatusOkP(status)) return;

  jitDebug(2,"sc2dalib_seqInit: call utils_update_debug\n");  
  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);
  SdpGeti("SETUP", &setup, status);
  if(setup == 0)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_seqInit: the sequence has not been setup" ); 
    return;
  }
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if ( in_sequence !=DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_seqInit: %s has not completed",
         seqStatus[myInfo->actionFlag]);
    return;
  } 
  myInfo->actionFlag=SEQACTION;

  jitDebug(2,"sc2dalib_seqInit: call DitsGetArgument()\n");  
  // get the argument, report if error.
  argId = DitsGetArgument();
  jitArgGetI(argId, "START", 1, range, 0, 0, &startSeq, status);      
  jitArgGetI( argId, "END", 2, range, 10, 0, &endSeq, status);
  jitArgGetI(argId, "DWELL", 3, range, 0, 0, &dWell, status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_seqInit: failed to get variables");
    return;
  }       

  SdpPuti("SEQ_START", startSeq, status);
  SdpPuti("SEQ_END", endSeq, status);
  SdpPuti("SEQ_DWELL", dWell, status);

  ArgPuti(myInfo->statId,"SEQ_START", startSeq,status);
  ArgPuti(myInfo->statId,"SEQ_END", endSeq,status);
  ArgPuti(myInfo->statId,"SEQ_DWELL", dWell,status);

  fprintf(myInfo->fpLog,"\n<%s> CMD for sc2dalib__Seq ",dateTime);
  fprintf(myInfo->fpLog,"SEQ_START <%ld> SEQ_END <%ld> DWELL <%ld>\n",
          startSeq, endSeq, dWell);
	  
  jitDebug(2,"sc2dalib_seqInit: call utils_alloc_memory\n");

  /* Calculate the number of stair steps that will be taken if this is a 
     fast flat field SEQUENCE. Protect against stairWidth being = zero */
  if(myInfo->stairWidth < 1)  myInfo->stairWidth = 1;
  totalStairCount =  (endSeq - startSeq + 1 - extraSteps) / myInfo->stairWidth;
  // fprintf(myInfo->fpLog,"totalStairCount: %d\n",totalStairCount);

  // allocate shared memory (dasdrama<=> dhtask) here. 
  // since we know how big it is now
  utils_alloc_memory(con, myInfo, startSeq, endSeq, status);
  if (!StatusOkP(status)) 
    return;
  jitDebug(2,"sc2dalib_seqInit: utils_alloc_memory OK\n");

  // initial the ith seqNum for checking parameter
  myInfo->firstQLFlag=1;
  myInfo->headtaskchkFlag=0;
  myInfo->taskchkFlag=0;
  myInfo->headerwriteNo=0;
  myInfo->headerreadNo=0;
  myInfo->askSDF=0;
  
  myInfo->parshmPtr->startSeq=startSeq; 
  myInfo->ithseqWait=startSeq + myInfo->parshmPtr->procNo -1;

  // since data is saved in dhtask, don't need to open file
  // don't save cmd to fpMcecmd
  myInfo->filesvFlag=0;

  table=(int*)calloc(myInfo->qlNum,sizeof(int));
  *lookupTable=table;

  ArgGetString(myInfo->statId,"OBS_MODE",40, obsMode, status);
}



/**
 * \fn void sc2dalib_seqcallsc2headmanseqStart(dasInfoStruct_t *myInfo,
 *  long startSeq, long endSeq,StatusType *status)
 *
 * \brief function
 *  call sc2headman-seqstart
 *
 * \param myInfo   dasInfoStruct_t pointer
 * \param startSeq long
 * \param endSeq   long
 * \param status    StatusType.  given and return
 *
 */
/*+ sc2dalib_seqcallsc2headmanseqStart   
 */
void sc2dalib_seqcallsc2headmanseqStart
(
dasInfoStruct_t       *myInfo,
long                  startSeq,
long                  endSeq,
StatusType            *status
)
{
  if (!StatusOkP(status)) return;
  jitDebug(2,
  "sc2dalib_seqcallsc2headmanseqStart: call sc2headman_seqstart startSeq(%d) EndSeq(%d)\n",
     (int)startSeq, (int)endSeq); 
  sc2headman_seqstart((int)startSeq, (int)endSeq,myInfo->jcmtheadEntry, status);
  if ( !StatusOkP(status) )
    ErsRep(0,status,"sc2dalib_seqcallsc2headmanseqStart: sc2headmanseqStart failed"); 
} 



/**
 * \fn void sc2dalib_seqcallsc2headmanmainHead(dasInfoStruct_t *myInfo,
 *  long startSeq, long endSeq,StatusType *status)
 *
 * \brief function
 *  call sc2headman-mainhead,Fill-in FITS main header records at start and end times
 *
 * \param myInfo   dasInfoStruct_t pointer
 * \param startSeq long
 * \param endSeq   long
 * \param status    StatusType.  given and return
 *
 */
/*+ sc2dalib_seqcallsc2headmanmainHead   
 */
void sc2dalib_seqcallsc2headmanmainHead
(
dasInfoStruct_t       *myInfo,
long                  startSeq,
long                  endSeq,
StatusType            *status
)
{
  size_t        actfits;         /* number of FITS records */

  if ( !StatusOkP(status) )     return;

   // shared par memory 
  //printf("sc2dalib_seqcallsc2headmanmainHead: pass MAXFITS=%d\n", MAXFITS);
  sc2headman_mainhead(myInfo->parshmPtr->subscanNo,(int)startSeq, (int)endSeq, 
		      myInfo->utcshort, MAXFITS, &actfits, myInfo->parshmPtr->fitshd, status);
  if ( !StatusOkP(status) )
    ErsRep(0,status,"sc2dalib_seqcallsc2headmanmainHead: sc2headmanmainHead failed"); 
  
  //printf("sc2dalib_seqcallsc2headmanmainHead: actfits returned=%d\n", actfits);
  myInfo->parshmPtr->actfits=actfits;
  jitDebug(2,"sc2dalib_seqcallsc2headmanmainHead: actfits=%d\n",
           myInfo->parshmPtr->actfits);
} 


/**
 * \fn void sc2dalib_seqchkQL(SDSU_CONTEXT *con, DRAMA_INNERMSG *dhmsg,
 *  dasInfoStruct_t *myInfo, DitsReasonType entReason,
 *  int *localTable,StatusType *status)
 *
 * \brief function
 *  updata QL if all drama parameters are received 
 *
 * \param con     SDSU context structure pointer
 * \param dhmsg    DRAMA_INNERMSG  struct pointer
 * \param myInfo  dasInfoStruct_t pointer
 * \param entReason DitsReasonType 
 * \param localTable  int pointer for local look up table 

 * \param status   StatusType.  given and return
 *
 * each time, DHTASK puts subendseq in a look-up-table (beginning at myInfo->lkupflagentPtr)
 * when a block of frames is ready and signal DRAMATASK to check if parameters are ready.
 *
 * For SCAN, DHTASK signals DRAMATASK twice for each subscan in a SEQ, Firstly when 
 * the block of subscan frames is ready and secondly when the subscan file is written.
 *
 * this call checks if all parameters are collected up to subendseq from the look-up-table
 * after the first block of frames ready from DHTASK.(no need to check before that).
 *
 * The look-up-table is always reset=0 during seqInit, the headerreadNo is always 
 * reset=0 for SEQ, when first block of frames are ready from DHTASK 
 *
 * if parameters up to the current subendseq are collected, then headrreadNo is increased 
 * by one for next subendseq, otherwise, it stays the same for next call to check. If 
 * headerreadNo = maxNo of Ix imag for STARE/DREAM or maxNo of subcan for SCAN,it increases
 * no more.
 *
 */
/*+ sc2dalib_seqchkQL   
 */
void sc2dalib_seqchkQL
(
SDSU_CONTEXT          *con,
DRAMA_INNERMSG        *dhmsg,         
dasInfoStruct_t       *myInfo,
DitsReasonType        entReason,
int                   *localTable,
StatusType            *status
)
{
  int     *ithscannumPtr=NULL;
  int     *start=NULL, *end=NULL; 
  double  *timePtr=NULL;
  char    *scanfileName=NULL;
  static long  startsubSeq, endsubSeq, endSeq;

  if (!StatusOkP(status)) return;

  jitDebug(8,"_seqchkQL: headerreadNo= %d  headerwriteNo= %d askSDF= %d at the entry\n",
           myInfo->headerreadNo,myInfo->headerwriteNo, myInfo->askSDF); 

  if ( myInfo->headtaskchkFlag==0)
  {
    // no block of frames are ready from dhtask, so, don't check parameters
    // this is where par-monotoring comes before first block frame ready
    return;
  }

  // firstFlag rest=0 after first seqsendQL 
  if (myInfo->firstQLFlag==1)
  {
    SdpGeti("SEQ_START", &startsubSeq, status);  
    SdpGeti("SEQ_END", &endSeq, status);
  }
 
  if ( entReason==DITS_REA_ASTINT  && dhmsg->reason == FRAME_TRIGQL)
  {
    //  copy lookupTable to localTable
    // headerwriteNO is reset =0 at the beginning of SEQ
    *(localTable + myInfo->headerwriteNo) = *( myInfo->lkupflagentPtr + myInfo->headerwriteNo);
    if ( (myInfo->headerwriteNo+1) < myInfo->qlNum)
    {
        myInfo->headerwriteNo ++;
    }
  }

  endsubSeq=*(localTable + myInfo->headerreadNo);
  myInfo->taskchkFlag=0;

  // paraMonitor can get here before frame ready, as headerreadNo is changed
  // after last seqsendQL. so check if frame ready
  if(endsubSeq == 0)
  {
    // frame not ready
    jitDebug(8, "sc2dalib_seqchkQL: frames not ready \n");
    return;
  }
  else 
  {
   jitDebug(8, "sc2dalib_seqchkQL: endsubSeq=%d \n",endsubSeq);
  }
  // used for scan QL ready 
  if ( entReason==DITS_REA_ASTINT && dhmsg->reason == FRAME_SCANREADY ) 
  {
    // scan already has sdf file written, so sendQL
    // only seqsendQL if it is signaled by dhtask
    jitDebug(8,"_seqchkQL: SDF done. seqsendQL *ithlkupflag %#0X (YDER=>REDY)\n",
             *(myInfo->lkupflagentPtr + myInfo->headerreadNo) );
    sc2dalib_seqsendQL(con,dhmsg,myInfo,status);
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"sc2dalib_seqchkQL: for scan: sc2dalib_seqsendQL failed");
      return;
    }
    if ( ( myInfo->headerreadNo+1) < myInfo->qlNum)
    {
      myInfo->headerreadNo++;
      myInfo->askSDF=0;
      startsubSeq=endsubSeq + 1; 
      if ( startsubSeq < endSeq)
      {
        jitDebug(8,"sc2dalib_seqchkQL: SCAN next startsubSeq=%ld\n",startsubSeq);
      }
    }
    return;
  }
  else
  {
    if ( myInfo->engFlag==RTSC_MODE )
    {
      // during waiting for SDF to be written in DHTASK, no check for task parameters
      if ( (myInfo->parshmPtr->obsMode ==OBS_SCAN) && (myInfo->askSDF ==1) )
      {
        jitDebug(8,"_seqchkQL: wait for SDF file to be written in DH\n");
        return;
      }

      sc2headman_checktasks(endsubSeq, &myInfo->taskchkFlag,status);
      if ( !StatusOkP(status) )  
      { 
        ErsRep(0,status,
            "sc2dalib_seqchkQL:sc2headman_checktasks(endsubSeq = %ld ) failed",
             endsubSeq);
        return;
      }
      // if all monitored tasks STATE are collected, call sendQL, 
      // if no monitoring task at all, taskchkFlag==1 
      if (myInfo->taskchkFlag==1)
      {
        // frame is ready update the ithseqWait,
        myInfo->ithseqWait =endsubSeq;
        jitDebug(8,"sc2dalib_seqchkQL: RTSmode ithseqWait(%d) headerreadNo(%d)\n",
                    myInfo->ithseqWait, myInfo->headerreadNo); 

        // also check endSeq here too. otherwise, we can stuck if no monitor task
        // headersDone==1 if no monitoring task at all
        sc2headman_checktasks(endSeq, &myInfo->headersDone,status);
        if ( !StatusOkP(status) )  
        { 
          ErsRep(0,status,
               "sc2dalib_seqchkQL:sc2headman_checktasks(endSeq= %ld ) failed",
               endSeq);
          return;
        }
        // not scan mode, so, sendQL 
        if (myInfo->parshmPtr->obsMode !=OBS_SCAN) 
        { 
          sc2dalib_seqsendQL(con,dhmsg,myInfo,status);
          if ( !StatusOkP(status) )  
          { 
            ErsRep(0,status,"sc2dalib_seqchkQL sc2dalib_seqsendQL failed");
            return;
          }

          if ( ( myInfo->headerreadNo+1) < myInfo->qlNum)
            myInfo->headerreadNo++;

          jitDebug(8,"_seqchkQL: headerreadNo= %d after seqsendQL myInfo->qlNum=%d\n",
                     myInfo->headerreadNo, myInfo->qlNum);
        }
        else  
        {
          // scan mode, tell dhtask that all FITs headers are filled and  
          // ask dhtask to write sdf file, myInfo->parshmPtr->fitshd is used here
          if ( startsubSeq > endSeq )
          {
            *status=DITS__APP_ERROR;
            ErsRep(0,status,"sc2dalib_seqchkQL: SCAN, startsubSeq > endSeq");
            post_sem(con->intask_sem,0); // allow dhtask to get out wait_sem
            return;
          }
          sc2headman_mainhead(myInfo->parshmPtr->subscanNo,(int)startsubSeq, 
			      (int)endsubSeq, myInfo->utcshort, MAXFITS, 
			      &myInfo->parshmPtr->actfits, myInfo->parshmPtr->fitshd,status);

          // whatever, have to allow dhtask to get out wait_sem
          post_sem(con->intask_sem,0);
          if (!StatusOkP(status)) 
          {
            ErsRep(0,status,
               "sc2dalib_seqchkQL: sc2headman_mainhead failed");
            return;
          } 
          ithscannumPtr = (int*)(myInfo->sharemscanPtr + myInfo->headerreadNo );
          timePtr =(double *)(ithscannumPtr+1);
          scanfileName=(char*)(timePtr+1);
          start       = (int *)(scanfileName + FILE_LEN);
          end    = start +1; 
          *start=startsubSeq;
          *end  =endsubSeq;
          myInfo->askSDF=1;
        }
      }
    }
    else
    {
       // frame is ready update the ithseqWait,
       myInfo->ithseqWait =endsubSeq;
      jitDebug(8,"sc2dalib_seqchkQL: eng_mode ithseqWait %d\n",myInfo->ithseqWait);
      if (myInfo->parshmPtr->obsMode !=OBS_SCAN)
      { 

        sc2dalib_seqsendQL(con,dhmsg,myInfo,status);
        if ( !StatusOkP(status) )  
        { 
           ErsRep(0,status,"sc2dalib_seqchkQL sc2dalib_seqsendQL failed");
           return;
        }
        myInfo->headerreadNo++;
      }
      else 
      {
        post_sem(con->intask_sem,0);
      }
    }
  }
}


/**
 * \fn void sc2dalib_seqsendQL(SDSU_CONTEXT *con, DRAMA_INNERMSG *dhmsg,
 *  dasInfoStruct_t *myInfo, StatusType *status)
 *
 * \brief function
 *  updata QL if all drama parameters are received 
 *
 * \param con     SDSU context structure pointer
 * \param dhmsg    DRAMA_INNERMSG  struct pointer
 * \param myInfo  dasInfoStruct_t pointer
 * \param status   StatusType.  given and return
 *
 */
/*+ sc2dalib_seqsendQL   
 */
void sc2dalib_seqsendQL
(
SDSU_CONTEXT          *con,
DRAMA_INNERMSG        *dhmsg,         
dasInfoStruct_t       *myInfo,
StatusType            *status
)
{
  size_t          actfits;
  int             *nrec=NULL;
  int             *start=NULL,*end=NULL;
  int             obsMode;  
  int             *ithcoaddnumPtr=NULL;
  int             *ithlkupflagPtr=NULL;
  double          *timePtr=NULL, *dataPtr=NULL;
  char            *scanfileName=" ";
  char            *fitshd=NULL;
  FITS4COADD      *fitsPtr=NULL;
  QL_IMAGE        *coaddImage=NULL;
  RECONSTR_IMAGE  *reconstrImage=NULL;
  
  long             seqStart=0;
  static  int      qlNum;
  static  uint32   dataDim[2],fitsDim[2];
  static  long     startsubSeq, endsubSeq,endSeq;

  if (!StatusOkP(status)) return;

  jitDebug(8,"sc2dalib_seqsendQL: qlRecod= %d-th QL trig\n",
           *myInfo->qlPtr+1);

  obsMode=myInfo->parshmPtr->obsMode;
 
  // check if it is the first QL trig  
  // point to the right entry at shared memory, static value
  if (myInfo->firstQLFlag >=1)
  {
    qlNum=0;

    SdpGeti("SEQ_START", &startsubSeq, status);
    SdpGeti("SEQ_END", &endSeq, status);

    if ( obsMode !=OBS_DREAM)
    {
      dataDim[SC2STORE__COL_INDEX]=COL_NUM;   dataDim[SC2STORE__ROW_INDEX]=(ROW_NUM-1);
    }
    else
    {
      dataDim[1]=myInfo->parshmPtr->dreamDim[1];   
      dataDim[0]=myInfo->parshmPtr->dreamDim[0];
    }
    fitsDim[0]=FITSSIZE;
    myInfo->firstQLFlag=0;
  }
  seqStart=startsubSeq;

  jitDebug(8,"sc2dalib_seqsendQL: myInfo->firstQLFlag=%d\n",
           myInfo->firstQLFlag );

  // update the endsubseq
  // (*ithlkupflagPtr=count4QL+parshmptr->startSeq-1 in dhtask)

  ithlkupflagPtr = myInfo->lkupflagentPtr + qlNum;
  endsubSeq=myInfo->ithseqWait;

  if ( obsMode==OBS_STARE || obsMode==OBS_DREAM)
  {
    jitDebug(8,"sc2dalib_seqsendQL: stare/dream lkupflag= %d\n", 
             *ithlkupflagPtr);
    if ( obsMode==OBS_DREAM )
      ithcoaddnumPtr = (int*)(myInfo->sharemrecnstrPtr + qlNum );
    else
      ithcoaddnumPtr = (int*)(myInfo->sharemqlPtr + qlNum );

    timePtr =(double *)(ithcoaddnumPtr+1);
    fitsPtr=(FITS4COADD *) (timePtr +1);
    dataPtr=(double*) (fitsPtr+ 1);
    fitshd=(char*)fitsPtr;
    if ( obsMode==OBS_DREAM )
    {
      reconstrImage= (RECONSTR_IMAGE *)(timePtr +1);
      start  =(int *) (reconstrImage +1);
    }
    else
    {
      coaddImage= (QL_IMAGE *)(timePtr +1);
      start  =(int *) (coaddImage +1);
    }
    end    = start +1; 
    nrec   = end +1;

    jitDebug(8, "sc2dalib_seqsendQL: seqStart (%d) endsubSeq(%d)\n",
            (int)seqStart, (int)endsubSeq);

    if ( startsubSeq > endSeq )
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,"sc2dalib_seqsendQL: STARE/DREAM, startsubSeq > endSeq");
      return;
    }

    if (myInfo->engFlag==RTSC_MODE )
    {
      sc2headman_subhead((int)seqStart, (int)endsubSeq, MAXQLFITS,&actfits,
                       fitshd,status);
      if (!StatusOkP(status)) 
      {
        ErsRep(0,status,"sc2dalib_seqsendQL: sc2headman_subhead failed");
        return;
      }
      //fitsDim[1]=(uint32)actfits;
      *nrec =(int)actfits;
      *start=startsubSeq;
      *end  =endsubSeq;
    }
    startsubSeq=endsubSeq + 1;
  }
  else if ( obsMode ==OBS_SCAN )
  {
    ithcoaddnumPtr = (int*)(myInfo->sharemscanPtr + qlNum );
    timePtr =(double *)(ithcoaddnumPtr+1);
    scanfileName=(char*)(timePtr+1);
    start       = (int *)(scanfileName + FILE_LEN);
    end    = start +1; 
    *start=startsubSeq;
    *end  =endsubSeq;
    startsubSeq=endsubSeq + 1;
  }

  /*  sc2headman_qlhead fails anytime the telescope is not participating
      in the current SEQUENCE.  This happens whenever we are not in RTS 
      mode, in all DARKs, and in SEQUENCEs where the rtsClient is in simulation
      like FLATFIELD. All of these types of SEQUENCEs now have obsMode == SCAN
      and since SCAN does not use the FITS stuff in the QL structure anyway, I 
      am just not going to call sc2headman_qlhead for any SCAN type observations */


  /* If not in RTS mode do not do anything with the QL structure */
  if ( myInfo->parshmPtr->engFlag==RTSC_MODE)
  {
    if(obsMode != OBS_SCAN && (myInfo->parshmPtr->load != LOAD_DARK) )
  	{
  	  jitDebug(8,
         "sc2dalib_seqsendQL: call sc2headman_qlhead startSeq (%d) endSeq(%d)\n",
  		   (int)seqStart, (int)endsubSeq);

  	  sc2headman_qlhead(myInfo->parshmPtr->subscanNo,(int)seqStart,
                        (int)endsubSeq, myInfo->utcshort, MAXFITS, &actfits,
                        myInfo->parshmPtr->fitshd, status);
  	  if ( !StatusOkP(status) )
	    {
	      ErsRep(0,status,"sc2dalib_seqsendQL: call sc2headman_qlHead failed"); 
	      return;
	    } 
  	}

	  fitsDim[1]=actfits;  
	  utils_update_QLSDS(con, myInfo, ithcoaddnumPtr, timePtr,
                       myInfo->parshmPtr->fitshd, fitsDim, dataPtr,
                       dataDim, scanfileName, status);
	  if (!StatusOkP(status))
	    return;
  }

  // just put a value for lookup flag table, increase qlNum
  *ithlkupflagPtr= QL_DONE;
  qlNum ++;
  *myInfo->qlPtr=qlNum; 

  // write single coadd struct into a temp file and display if setdebug=2
  sc2dalib_coadddatawriteDisp(con,myInfo,scanfileName,(char*)ithcoaddnumPtr,
			      dataPtr,status); 
}


/**
 * \fn void sc2dalib_seqchkEnd(SDSU_CONTEXT *con, DRAMA_INNERMSG *dramamsg,
 *  dasInfoStruct_t *myInfo, StatusType *status)
 *
 * \brief function
 *  check reason for ending sequence action and report. inform driver not
 *  to expect any data
 *
 * \param con     SDSU context structure pointer
 * \param dramamsg DRAMA_INNERMSG  struct pointer
 * \param myInfo  dasInfoStruct_t pointer
 * \param status   StatusType.  given and return
 *
 */
/*+ sc2dalib_seqchkEnd   
 */
void sc2dalib_seqchkEnd
(
SDSU_CONTEXT          *con,
DRAMA_INNERMSG        *dramamsg,         
dasInfoStruct_t       *myInfo,
StatusType            *status
)
{
  time_t  tm;

  if (!StatusOkP(status)) return;

#ifdef PRINT_LOG_MESG
  fprintf(dasInfo.fpLog,"sc2dalib__Seq: received FRAME trig\n");
#endif

  //jitDebug(2,"sc2dalib__Seq: received FRAME trig\n");

  con->datacount=0;  // inform driver non data expected
  con->process.whereabout=Dits_GetSeqn;

  // Write the time into the logfile
  tm = time(NULL);
  fprintf (myInfo->fpLog,"%s",asctime(localtime(&tm)) );

  // use 0 MsgOut,3  printf; 4 no print at all, all save to log
  utils_msg(myInfo,"sc2dalib_seqchkEnd: %s",
                    dramamsg->endMsg,NO_PRINT,status);
 
  if (dramamsg->reason !=FRAME_COMPLETION)
  {
    if(dramamsg->reason==FRAME_STOPPED)
    {
      utils_msg(myInfo,
        "sc2dalib_seqchkEnd: the action is stopped by kick","",USE_MSGOUT,status);
      con->process.seqstatus=SEQ_STOPPED;
    }
    else
    {
      *status=DITS__APP_ERROR;
      con->process.seqstatus=SEQ_ERROR;    

      utils_msg(myInfo,"sc2dalib_seqchkEnd: -----ERROR-------------",
                        "",USE_ERSREP,status); 
      if (dramamsg->reason==FRAME_ERRGETPASSBUF)
      {         
        utils_msg(myInfo," Ended with waiting data buffer TIMEOUT",
                          "",USE_ERSREP,status);
      }
      else if (dramamsg->reason==FRAME_MISSING)
      {  
        utils_msg(myInfo," Ended with less frames then asked",
                          "",USE_ERSREP,status);
      }
      else if (dramamsg->reason==FRAME_MORE)
      {  
        utils_msg(myInfo, 
          "either corrupted frame status or too many frames","",USE_ERSREP,status);
      }
      else if (dramamsg->reason==FRAME_MALLOCFAILED)
      {  
        utils_msg(myInfo," %s",dramamsg->errRep,USE_ERSREP,status); 
      }
      else 
      {
       utils_msg(myInfo," reason unknown or","",USE_ERSREP,status);
       utils_msg(myInfo," %s",dramamsg->errRep,USE_ERSREP,status);
      }
    }
    // in case of error, need to unblock dh task, but not wait for 
    // parameters
    post_sem(con->intask_sem,0);
  }
  else
  {
    con->process.seqstatus=SEQ_FINISHED;
    // read in all parameters and placed them in the sharedMem
    // post sem to dh task
     post_sem(con->intask_sem,0);
  }
}


/**                       
 * \fn void sc2dalib_seqEnd (SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *    int *lookupTable,StatusType *status)
 *
 * \brief function
 *  free shared memory and end seq Drama action 
 *
 * \param con      SDSU_CONTXT structure pointer
 * \param myInfo  dasInfoStruct_t structure pointer 
 * \param lookupTable int pointer
 * \param status   StatusType.  given and return
 */
/*+ sc2dalib_seqEnd - 
*/
void sc2dalib_seqEnd
(
SDSU_CONTEXT    *con,           
dasInfoStruct_t *myInfo,       
int             *lookupTable,
StatusType      *status
)
{
  StatusType   inStatus=*status;
  *status=STATUS__OK;

  if ( lookupTable !=NULL)
     free (lookupTable);

  if (myInfo->engFlag==RTSC_MODE)
    sc2headman_endseq(status);

  utils_close_memory(myInfo,SHAREDM_OBS,status);
  sc2dalib_endAction(con,myInfo,status);
  *status=inStatus;
}


/**
 * \fn void sc2dalib_stopFrame(SDSU_CONTEXT    *con,                    
 *   dasInfoStruct_t *myInfo, struct mcexml_struct  *mceinxPtr,
 *   char *dateTime, StatusType *status)
 *
 * \brief function
 *   send ST to MCE to stop frame taking
 *
 * \param con         SDSU context structure
 * \param myInfo     dasInfo structure pointer 
 * \param mceinxPtr   mcexml_struct pointer for parameter lookup table
 * \param dateTime   dateTime string pointer         
 * \param status      StatusType pointer
 * 
 */
/*+ sc2dalib_stopFrame
*/
void sc2dalib_stopFrame
(                  
SDSU_CONTEXT         *con,
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceinxPtr,
char                 *dateTime,
StatusType           *status
)
{
  dasCmdInfo_t abortCmd;
  StatusType   inStatus=*status;

  *status=STATUS__OK;

  strcpy(abortCmd.mceCmd,myInfo->stopCmd);
  myInfo->filesvFlag=0;
  // send cmd and get reply, save them to logfile but not save to other
  // cmdreply file specified by user
  sc2dalib_sendCmd(con,myInfo,&abortCmd,mceinxPtr,dateTime,status);               
  if ( !StatusOkP(status) )  
    {
      ErsRep(0,status, "sc2dalib_stopFrame: sc2dalib_sendCmd failed");
      return;
    }
  else
    {
      *status=inStatus;
    }
  
}


/**
 * \fn int sc2dalib_getStats( const void *a, const void *b)
 *
 *  double    theArray[]     The proposed new heater setting for each pixel (given)
 *  int       index[]        Just an array to get sorted the same way as theArray to retain the index
                             Needs to be filled with 0,1,2,3, ...
 *  double    *theMedian;    The median value in theArray
 *  int       *theIndex;     The starting index of the median value (where it was in the original theArray) 
 *  int       n,             The number of values to sort
 *  double    *avg,          The average of theArray (throw out +/- %5)
 *  double    *sigma;        The standard deviation of theArray
 *
 * Sorts theArray based on value and sorts index the same way a theArray
 * Returns the median, the index of the median, the average and the standard deviation
 *
 */
/*+ sc2dalib_getStats
 */
void sc2dalib_getStats
(
 double theArray[],
 int index[],
 double *theMedian,
 int *theIndex,
 int n,
 double *avg,
 double *sigma
)
{
  double temp;
  int i, j, tempIndex;
  float sum, sumSqu, sigmaSqu;

  /* No use messing around if we only use 2 pixels */
  *theMedian = theArray[0];
  *theIndex = index[0];
  *avg = 0.0;
  *sigma = 0.0;
  if(n < 3) return;

  /* Sort the array by magnitude, do the same sort on index */
  for(i=(n-2); i >= 0; i--)
    {
      for(j=0;j<=i;j++)
	{
	  if(theArray[j]>=theArray[j+1])
	    {
	      temp=theArray[j];
	      tempIndex = index[j];
	      theArray[j]=theArray[j+1];
	      index[j]=index[j+1];
	      theArray[j+1]=temp;
	      index[j+1] = tempIndex;
	    }
	}
    }
  *theMedian = theArray[n/2];
  *theIndex = index[n/2];

  /* Now throw out the top %5 and bottom 5% of the outliers
     and calculate an average and standard deviation */
  i = 0.05 * n;
  tempIndex = 0;
  sum = 0.0;
  sumSqu = 0.0;

  for(j=i; j<(n-i); j++)
    {
      sum = sum + theArray[j];
      sumSqu = sumSqu + (theArray[j] * theArray[j]);
      tempIndex += 1;
    }

  temp = tempIndex;
  *avg = sum / temp;
  sigmaSqu = (temp * sumSqu - (sum * sum)) / (temp * (temp - 1));
  if(sigmaSqu > 0.0)
    {
      *sigma = sqrt(sigmaSqu);
    }

}

/**
 * \fn void sc2dalib_getsq2fbparaInit(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *   dasCmdInfo_t *myCmd, struct mcexml_struct  *mceinxPtr,
 *   char *dateTime, ARRAYSET *arrayset,StatusType *status)
 *
 * \brief function
 *  get args for drama action SQ2FBGETPARA and open stripchart file
 *
 * \param con       SDSU context structure pointer 
 * \param myInfo    dasInfoStruct_t pointer
 * \param myCmd      dasCmdInfo_t pointer, return go cmdBuf
 * \param mceinxPtr   mcexml_struct pointer for parameter lookup table
 * \param dateTime   dateTime string pointer         
 * \param arrayset      ARRAYSET structure pointer ( point to global one)
 * \param status    StatusType     
 *
 * hard wire $CURRENTDATADIR and $CONFIG_HARD in sc2dalib.c so that
 * RTSC can call it 
 */
/*+ sc2dalib_getsq2fbparaInit
*/
void sc2dalib_getsq2fbparaInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
char                  *dateTime,
ARRAYSET              *arrayset,
StatusType            *status
)
{
  int           i;
  long         in_sequence,configured;
  DitsArgType  argId;
  char         tmp[FILE_LEN],tmp1[FILE_LEN];
  char         readsq2fbCmd[]="rb bc2 flux_fb 32";

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  SdpGeti("CONFIGURED", &configured, status);
  if (configured==0 )
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_getsq2fbparaInit: the DA is not configured" ); 
    return;
  }
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_getsq2fbparaInit: %s has not completed",
             seqStatus[myInfo->actionFlag]);
    return;
  } 
  myInfo->actionFlag=GETSQ2FBPARAACTION;

  // readback the sample_num to re-scale the data and dataMode 
  sc2dareadmceparVal(con,myInfo,arrayset,mceinxPtr,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status,"sc2dalib_getsq2fbparaInit: failed to call sc2dareadmceparVal");
    return;
  }

  argId = DitsGetArgument();
  jitArgGetS(argId, "SETUP_FILE", 1,NULL, "setup-sq2fbgetpara", 0, FILE_LEN,
	     tmp1, NULL, status);
  //print ?
  sprintf(myInfo->batchFile, "%s/%s",getenv("SC2SCRATCH"), tmp1);
  myInfo->filesvFlag=0;
  jitArgGetS(argId, "SQ2OPT_FILE", 2,NULL, "sq2fboptimal.txt", 0, FILE_LEN,
	     tmp, NULL, status);
  jitDebug(16,"_readOPT %s \n",tmp);
  // read sq2fboptimal in ->sq2fdbkOpt
  sc2dalib_readOPT(myInfo, arrayset->sq2fdbkOpt, COL_NUM, tmp, status);
  if ( !StatusOkP(status) )
    return;
  jitArgGetS(argId, "SQ1OPT_FILE", 3,NULL, "sq1biasoptimal.txt", 0, FILE_LEN,
	     tmp,  NULL,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_getsq2fbparaInit: failed to get variables.");
    return;
  }
  // read sq1biasoptimal in ->sq1biasOpt
  jitDebug(16,"_readOPT %s \n",tmp);
  sc2dalib_readOPT(myInfo,arrayset->sq1biasOpt,ROW_NUM,tmp,status);
  if ( !StatusOkP(status) )
    return;
  
  // now get the args from setup file
  // put all include=xxx in the setup file into a single file  
  // copy it to batchFile, for sc2dalibsetup_readsq2fbSetup
  sc2dalibsetup_servoreadsetupWrap(myInfo,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,
        "sc2dalib_getsq2fbparaInit: sc2dalibsetup_servoreadsetupWrap failed"); 
    return;
  }
  utils_fclose(&(myInfo->fpBatch));
  if((myInfo->fpBatch=fopen(myInfo->batchFile,"r")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_getsq2fbparaInit: Error- failed to open file %s", myInfo->batchFile); 
      return;
    }
  jitDebug(16,"sc2dalib_getsq2fbparaInit read: %s\n",myInfo->batchFile);

  sc2dalibsetup_servoreadSetup(myInfo,arrayset,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_getsq2fbparaInit: sc2dalib_servoreadSetup failed"); 
    return;
  }
  jitDebug(16,"sc2dalib_getsq2fbparaInit finish servoreadSetup\n");
  utils_fclose(&(myInfo->fpBatch));

  if (arrayset->slopSelect[0] <0 )
  {
    *status=DITS__APP_ERROR;	
    ErsRep(0,status,"sc2dalib_getsq2fbparaInit: No row selected in SQ2FBTRACK slopSelrct[0]"); 
    return;
  }
  if (arrayset->doServo <=0 )
  {
    *status=DITS__APP_ERROR;	
    ErsRep(0,status,"sc2dalib_getsq2fbparaInit: sq2fbgetpara-base:doServoNo must > 0"); 
    return;
  }
  MsgOut(status, "_getsq2fbparaInit: doServo=%d, slopSelect[0] [5]=%d %d",
         arrayset->doServo,arrayset->slopSelect[0],arrayset->slopSelect[5]);

  // always remove previous tmp setupfile
  sprintf( tmp, "rm -f %s",myInfo->batchFile);
  system ( tmp);
 
  // fix stripchart name for disply tracking result, 
  //inside it, the col are: == No sq2fb delta sq2fbopt ...(32 col)  ==
  sprintf(myInfo->strchartFile,"%s/trksq2fb.txt", getenv("ORAC_DATA_OUT") );
  utils_fclose(&(myInfo->fpStrchart));
  if ((myInfo->fpStrchart=fopen64(myInfo->strchartFile,"a"))==NULL)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"_getsq2fbparaInit: failed to open %s",myInfo->strchartFile);
    ErsRep(0,status,"_getsq2fbparaInit: %s",strerror(errno));
    return;
  }
  fprintf(myInfo->fpLog,"\n<%s> CMD from sc2da_Getsq2fbPara\n",dateTime);
  myInfo->gotsq2fbparaFlag=1;
  myInfo->parshmPtr->sq2fbparaFlag=1;
  stripchFlag=0;

  // read sq2fb values from the MCE to sq2fdbkOPt
  mce_read(con,myInfo,mceinxPtr,readsq2fbCmd,arrayset->sq2fdbkOpt,32,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status, "sc2dalib_getsq2fbparaInit:mce_read failed to read sq2fb"); 
    return;
  }
  fprintf(myInfo->fpLog,"sq2fdbkOpt readback===== \n");
  for (i=0;i<32;i++)
    fprintf(myInfo->fpLog,"%d ",arrayset->sq2fdbkOpt[i]);
  fprintf(myInfo->fpLog,"\n ===== \n");
}


/**
 * \fn void sc2dalib_trksq2fbInit(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *   dasCmdInfo_t *myCmd, struct mcexml_struct  *mceinxPtr,
 *   char *dateTime, ARRAYSET *arrayset,long *option, int *sq2fboptInit,
 *   int *sq1biasOpt, StatusType *status)
 *
 * \brief function
 *  get args for drama action and open file
 *
 * \param con       SDSU context structure pointer 
 * \param myInfo    dasInfoStruct_t pointer
 * \param myCmd      dasCmdInfo_t pointer, return go cmdBuf
 * \param mceinxPtr   mcexml_struct pointer for parameter lookup table
 * \param dateTime   dateTime string pointer         
 * \param arrayset      ARRAYSET structure pointer ( point to global one)
 * \param option     long pointer
 * \param sq2fboptInit  int pointer,  sq2fbopt when sq1bias !=0
 * \param sq1biasOpt    int pointer,  sq1biasopt when sq1bias!=0
 * \param status    StatusType     
 *
 * hard wire $CURRENTDATADIR and $CONFIG_HARD in sc2dalib.c so that
 * RTSC can call it 
 */
/*+ sc2dalib_trksq2fbInit
*/
void sc2dalib_trksq2fbInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
char                  *dateTime,
ARRAYSET              *arrayset,
long                  *option,
int                   *sq2fboptInit,
int                   *sq1biasOpt,
StatusType            *status
)
{
  int          i;
  long         in_sequence,configured;
  long         range[]={0,2};
  DitsArgType  argId;
  int          wait;
  char        sq2fbCmd[]="wb bc2 flux_fb";
  char        readsq2fbCmd[]="rb bc2 flux_fb 32";

  char  mceCard[FILE_LEN];
  char  tmp[FILE_LEN],tmp1[FILE_LEN];

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);

  SdpGeti("CONFIGURED", &configured, status);
  if (configured==0 )
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_trksq2fbInit: the DA is not configured" ); 
    return;
  }
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"sc2dalib_trksq2fbInit: %s has not completed",
             seqStatus[myInfo->actionFlag]);
    return;
  } 
  myInfo->actionFlag=TRKSQ2FBACTION;

  // readback the sample_num to re-scale the data and dataMode 
  sc2dareadmceparVal(con,myInfo,arrayset,mceinxPtr,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status,"sc2dalib_trksq2fbInit: failed to call sc2dareadmceparVal");
    return;
  }

  argId = DitsGetArgument();
  jitArgGetI(argId, "OPTION", 1, range, 0, 0, option, status );
  //option default 0?
  myInfo->filesvFlag=2;
  jitArgGetS(argId, "DATA_FILE", 2, NULL, "data.txt", 0, FILE_LEN,
	     tmp, NULL, status);
  sprintf(myInfo->dataFile, "%s/%s",getenv("CURRENTDATADIR"), tmp);
  //print ?
  jitArgGetS(argId, "SQ2OPT_FILE", 4, NULL, "sq2fboptimal.txt", 0, FILE_LEN,
	     tmp, NULL, status);
  jitDebug(16,"_trksq2fbreadOPT %s \n",tmp);
  // read sq2fboptimal in and store in global sq2fboptInit
  sc2dalib_readOPT(myInfo,sq2fboptInit,COL_NUM,tmp,status);
  if ( !StatusOkP(status) )
    return;
  jitArgGetS(argId, "SQ1OPT_FILE", 5, NULL, "sq1biasoptimal.txt", 0, FILE_LEN,
	     tmp, NULL, status);
  // read sq1biasoptimal in and store in global sq1biasOpt
  jitDebug(16,"_trksq2fbreadOPT %s \n",tmp);
  sc2dalib_readOPT(myInfo,sq1biasOpt,ROW_NUM,tmp,status);
  if ( !StatusOkP(status) )
    return;
  jitArgGetS(argId, "SETUP_FILE", 3, NULL, "setup-trksq2fb", 0, FILE_LEN,
	     tmp1, NULL, status);
  sprintf(myInfo->batchFile, "%s/%s",getenv("SC2SCRATCH"), tmp1);
  //print ? sv default 0?
  jitArgGetI(argId, "SVFILE_FLAG", 4, range, 0, 0, &myInfo->filesvFlag, status);
  if ( !StatusOkP(status) )
  {
    ErsOut(0,status, "sc2dalib_trksq2fbInit: failed to get variables.");
    return;
  }
  
  // now get the args from setup file
  // put all include=xxx in the setup file into a single file  
  // copy it to batchFile, for sc2dalibsetup_readsq2fbSetup ?
  sc2dalibsetup_servoreadsetupWrap(myInfo,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,
        "sc2dalib_trksq2fbInit: sc2dalibsetup_servoreadsetupWrap failed"); 
    return;
  }
  utils_fclose(&(myInfo->fpBatch));
  if((myInfo->fpBatch=fopen(myInfo->batchFile,"r")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_trksq2fbInit: Error- failed to open file %s", myInfo->batchFile); 
      return;
    }
  jitDebug(16,"sc2dalib_trksq2fbInit read: %s\n",myInfo->batchFile);

  sc2dalibsetup_servoreadSetup(myInfo,arrayset,status);
  jitDebug(16,"sc2dalib_trksq2fbInit finish servoreadSetup\n");

  if (arrayset->slopSelect[0] <0 )
  {
    *status=DITS__APP_ERROR;	
    ErsRep(0,status,"sc2dalib_trksq2fbInit: No row selected in SQ2FBTRACK slopSelrct[0]"); 
    return;
  }
  if (arrayset->doServo <=0 )
  {
    *status=DITS__APP_ERROR;	
    ErsRep(0,status,"sc2dalib_trksq2fbInit: sq2fbtrack-base:doServoNo must > 0"); 
    return;
  }
  MsgOut(status, "_trksq2fbInit: doServo=%d, slopSelect[0] [5]=%d %d",
         arrayset->doServo,arrayset->slopSelect[0],arrayset->slopSelect[5]);

  // always remove previous tmp setupfile
  sprintf( tmp, "rm -f %s",myInfo->batchFile);
  system ( tmp);
  utils_fclose(&(myInfo->fpBatch));

  //populate cmd buff myCmd.cmdBuf,always use rcs now
  jitDebug(16,"populate cmd buff\n");
  sprintf(mceCard, "rcs");
  sprintf(myInfo->goCmd, "GO %s ret_dat 1",mceCard);
  sc2dalibsetup_whichRC(myInfo,mceCard,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_trksq2fbInit: not recognised(%s) as MCE_WHICHCARD",
           mceCard);
    return;
  }
  strcpy(myCmd->mceCmd,myInfo->goCmd);
  sc2dalib_getcmdBuf(myInfo,myCmd,mceinxPtr,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_trksq2fbInit: sc2dalib_getCmdBuf failed"); 
    return;
  }

   // fix stripchart name for disply tracking result, 
  //inside it, the col are: == No sq2fb delta sq2fbopt ...(32 col)  ==
  sprintf(myInfo->strchartFile,"%s/trksq2fb.txt", getenv("ORAC_DATA_OUT") );

  jitDebug(16,"_openFiles\n");
  utils_open_files(myInfo,DATFILEAPPEND,NOBATCHFILE,status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"sc2dalib_trksq2fbInit: utils_open_files failed"); 
    return;
  }
 
  fprintf(myInfo->fpLog,"\n<%s> CMD from sc2dalib__Trksq2fb\n",dateTime);
  SdpPuti("IN_SEQUENCE",DA_MCE_SINGLEDATA,status);

  // read sq2fb values from the MCE to sq2fdbkOPt
  mce_read(con,myInfo,mceinxPtr,readsq2fbCmd,arrayset->sq2fdbkOpt,32,status);
  if (!StatusOkP(status)) 
  {
    ErsRep(0,status, "sc2dalib_trksq2fbInit:mce_read failed to read sq2fb"); 
    return;
  }
  fprintf(myInfo->fpLog,"sq2fdbkOpt readback ===== \n");
  for (i=0;i<32;i++)
    fprintf(myInfo->fpLog,"%d ",arrayset->sq2fdbkOpt[i]);
  fprintf(myInfo->fpLog,"\n ===== \n");

  if (arrayset->slopSelect[5]==0)  // fbFlag=0 take data sq2fb=sq2fbINIT
  {
    arrayset->fbFlag=0;
    sc2dalib_sendarrayCmd(con,myInfo,mceinxPtr,dateTime,sq2fbCmd,arrayset->initFB,COL_NUM,status);
    if ( !StatusOkP(status) )
    {
      ErsRep (0, status,"sc2dalib_trksq2fbInit: _sendarrayCmd failed SQ2 FDBK (%s)",sq2fbCmd); 
    }
    if( arrayset->slopSelect[7] > 0 )
    {
      wait=(1000000*arrayset->slopSelect[6])/arrayset->slopSelect[7];
      usleep(wait);
    }
  }
}



/**
 * \fn void sc2dalib_trksq2FB(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceInxpt, char *dateTime,char *data, 
 *  ARRAYSET *arrayset, int *sq2fbOpt,  int option, int skipDo, StatusType *status)
 *
 * \brief function
 *  calculate delta sq2fb using pixel data and ssalock, gain obtained 
 *  from SQ2OPEN
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param mceInxpt  struct mcexml_struct pointer
 * \param dateTime  dateTime string pointer         
 * \param data     char pointer for the data         
 * \param setup    ARRAYSET structure pointer
 * \param sq2fbOpt int array, sq2fb optimal values
 * \param option   int 0: single, 1 continuous
 * \param skipDo   int  0: update sq2fb, otherwise skip upadte but still write trksq2fb.txt
 * \param status StatusType.  given and return
 *
 */
/*+ sc2dalib_trksq2FB
 */
void sc2dalib_trksq2FB
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                  *dateTime,
char                  *data,
ARRAYSET              *setup,
int                   *sq2fbOpt,
int                   option,
int                   skipDo,
StatusType            *status
)
{
  int            *frameData,i;
  int            dataInx, delta,refRow;
  char           sq2fbCmd[]="wb bc2 flux_fb";
  static double  gain[COL_NUM];
  static  int    coaddData[COL_NUM];
  static int     flag[COL_NUM];

  if (*status != STATUS__OK) return;

  refRow=setup->slopSelect[0];

  // skip the header and point to the refRow 
  frameData = (int *)( data+HEADERSET);  
  dataInx=refRow*COL_NUM;

  if (myInfo->trkNo==1)
  { 
    for(i=0;i<COL_NUM;i++)
    {
      if ( setup->slopSelect[5]==1 )
        gain[i]=setup->sq1gainScale[i]*setup->gain[i];
      else
        gain[i]=setup->gain[i];

      coaddData[i]=0;

      // use the first data as zfact, x Gao 2009.07.09
      setup->zfact[i] = *(frameData+dataInx + i); 

      if (setup->colMask[i] !=0 )
        MsgOut(status,"Gain[%d]=(%f)%f", i,setup->gain[i],gain[i]);
    }  
  }
 
 
  for(i=0;i<COL_NUM;i++)
     coaddData[i] += *(frameData+dataInx + i); 

  if (skipDo ==0)
  {
    for(i=0;i<COL_NUM;i++)
      coaddData[i] /=setup->doServo; 
  }    
  
  fprintf(myInfo->fpStrchart,"%7d ",stripchFlag);

  jitDebug(16,"sc2dalib_trksq2FB starts \n");
  for(i=0;i<COL_NUM;i++)
  {
    delta=0;
    if (setup->colMask[i] !=0)
      fprintf(myInfo->fpStrchart,"%7d ",setup->initFB[i]);
      
    if ( myInfo->trkNo >= SQ2FB_UPDATE_WAIT ) // allow some initial data do nothing 100
    {
      if (skipDo ==0 )
      {
        if (setup->initFB[i] !=0)  // masked one, kept =0
        {
          delta= ( coaddData[i] - setup->zfact[i] )*gain[i];
          setup->initFB[i] -= delta; 
        }
 
        if (sq2fbOpt[i] !=0)  // masked one, kept =0
          sq2fbOpt[i] -=delta;
 
       //  flux jump ???
        sc2dalib_fluxJump(myInfo,setup,setup->initFB,i,status);   
        sc2dalib_fluxJump(myInfo,setup,sq2fbOpt,i,status);   
      
        if ( setup->slopSelect[5]==1 && setup->slopSelect[8]==1 )
          sc2dalib_trksq2fbupdateGain(myInfo,setup,i,delta, gain, flag,status);   
      }
    }     
    if (setup->colMask[i] !=0)
    {
      fprintf(myInfo->fpStrchart,"%5d %5d ",coaddData[i],setup->zfact[i]);
      fprintf(myInfo->fpStrchart,"%5d %7d ",delta, setup->fluxPeriod[i]);
      fprintf(myInfo->fpStrchart,"%6.4E %3d %5d",gain[i],flag[i],sq2fbOpt[i]);
    }
    fprintf(myInfo->fpLog,"%7d ",setup->initFB[i]);
  }
  fprintf(myInfo->fpLog,"\n");
  fflush(myInfo->fpLog);
  fprintf(myInfo->fpStrchart,"\n");
  fflush(myInfo->fpStrchart);
  stripchFlag++;
  if (option ==0)
  {
    sc2dalib_trksq2fbsavefbZG(myInfo,setup,status);
    sc2dalib_trksq2fbsavesq2fbOPT(myInfo,sq2fbOpt,status);
  }
  else
  {
    if ( skipDo ==0 )
    {
      for(i=0;i<COL_NUM;i++)
        coaddData[i]=0;
      
      //update SQ2 feedback for continuous tracking
      if ( setup->slopSelect[5]==1 ) 
        sc2dalib_sendarrayCmd(con,myInfo,mceInxpt,dateTime,sq2fbCmd,sq2fbOpt,COL_NUM,status);
      else
      {   
        if ( setup->fbFlag==0 && setup->slopSelect[4]==1) 
          sc2dalib_sendarrayCmd(con,myInfo,mceInxpt,dateTime,sq2fbCmd,sq2fbOpt,COL_NUM,status);
      }
      if ( !StatusOkP(status) )
      {
        ErsRep (0, status,
        "sc2dalib_trksq2FB: sc2dalib_sendarrayCmd failed SQ2 FDBK (%s)",sq2fbCmd); 
        return;
      }
    }
  }
}


/**
 * \fn void sc2dalib_updatesq2fbVal(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceInxpt, char *dateTime,char *data, 
 *  ARRAYSET *arrayset, int *sq2fbOpt,  StatusType *status)
 *
 * \brief function
 *  calculate delta sq2fb using pixel data,lockpoint and gain obtained from SQ2OPEN
 *  and update sq2fbOPT
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param mceInxpt  struct mcexml_struct pointer
 * \param dateTime  dateTime string pointer         
 * \param data     char pointer for the data         
 * \param setup    ARRAYSET structure pointer
 * \param sq2fbOpt int array, sq2fb optimal values
 * \param status StatusType.  given and return
 *
 */
/*+ sc2dalib_updatesq2fbVal
 */
void sc2dalib_updatesq2fbVal
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                  *dateTime,
char                  *data,
ARRAYSET              *setup,
int                   *sq2fbOpt,
StatusType            *status
)
{
  int            *frameData,i;
  int            dataInx, delta,refRow;
  char           sq2fbCmd[]="wb bc2 flux_fb";
  static double  gain[COL_NUM];
  static  int    array[COL_NUM];
  static int     flag[COL_NUM];

  if (*status != STATUS__OK) return;

  refRow=setup->slopSelect[0];

  // skip the header and point to the refRow 
  frameData = (int *)( data+HEADERSET);  
  dataInx=refRow*COL_NUM;

  for(i=0;i<COL_NUM;i++)
  {
    if ( setup->slopSelect[5]==1 )
      gain[i]=setup->sq1gainScale[i]*setup->gain[i];
    else
      gain[i]=setup->gain[i];

    // use the first frame as zfact[i]
    if ( myInfo->trkNo < 1 )
      setup->zfact[i]=*(frameData+dataInx + i);
  }  
 
  for(i=0;i<COL_NUM;i++)
     array[i] += *(frameData+dataInx + i); 

  fprintf(myInfo->fpStrchart,"%7d ",stripchFlag);
  for(i=0;i<COL_NUM;i++)
  {
    delta=0;
    if (setup->colMask[i] !=0)
      fprintf(myInfo->fpStrchart,"%d\t",setup->initFB[i]);

    if (setup->initFB[i] !=0)  // masked one, kept =0
    {
      delta= ( array[i] - setup->zfact[i] )*gain[i];
      setup->initFB[i] -= delta; 
    }

    if (sq2fbOpt[i] !=0)  // masked one, kept =0
      sq2fbOpt[i] -=delta;

    //  flux jump ???
    sc2dalib_fluxJump(myInfo,setup,setup->initFB,i,status);   
    sc2dalib_fluxJump(myInfo,setup,sq2fbOpt,i,status);   

    if ( setup->slopSelect[5]==1 && setup->slopSelect[8]==1 )
       sc2dalib_trksq2fbupdateGain(myInfo,setup,i,delta, gain, flag,status);   

    if (setup->colMask[i] !=0)
    {
      fprintf(myInfo->fpStrchart,"%d\t%d\t",array[i],setup->zfact[i]);
      fprintf(myInfo->fpStrchart,"%d\t%d\t ",delta, setup->fluxPeriod[i]);
      fprintf(myInfo->fpStrchart,"%f\t%d\t%d\t",gain[i],flag[i],sq2fbOpt[i]);
    }
    //fprintf(myInfo->fpLog,"%7d ",setup->initFB[i]);
  }
  fprintf(myInfo->fpStrchart,"\n");
  //fprintf(myInfo->fpLog,"\n");
  fflush(myInfo->fpStrchart);
  //fflush(myInfo->fpLog);
  stripchFlag++;

  //update SQ2 feedback for continuous tracking
  if ( setup->slopSelect[5]==1 ) 
    sc2dalib_sendarrayCmd(con,myInfo,mceInxpt,dateTime,sq2fbCmd,sq2fbOpt,COL_NUM,status);
  else
  {   
    if ( setup->fbFlag==0 && setup->slopSelect[4]==1) 
        sc2dalib_sendarrayCmd(con,myInfo,mceInxpt,dateTime,sq2fbCmd,sq2fbOpt,COL_NUM,status);
  }
  if ( !StatusOkP(status) )
  {
    ErsRep (0,status,"_trksq2FB: sc2dalib_sendarrayCmd failed SQ2 FDBK (%s)",sq2fbCmd); 
    return;
  }
}


/**
 * \fn void sc2dalib_changesq2fbVal(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceInxpt, char *dateTime, int *valArray,
 *  StatusType *status)
 *
 * \brief function
 *   send wb bc2 flux_fb xxxx xxxx  to sq2fb 
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param mceInxpt  struct mcexml_struct pointer
 * \param dateTime  dateTime string pointer   
 * \param valArray   int array pointer      
 * \param status StatusType.  given and return
 *
 */
/*+ sc2dalib_changesq2fbVal
 */
void sc2dalib_changesq2fbVal
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                 *dateTime,
int                  *valArray,
StatusType           *status
)
{
  int       j; 
  char    tmp[30];
  char    sq2fbCmd[]  ="wb bc2 flux_fb";
  dasCmdInfo_t  myCmd;

  if (*status != STATUS__OK) return;

  sprintf (myCmd.mceCmd, "%s ",sq2fbCmd ); 
  // valArray[0[ is frameCount
  for ( j=0; j<COL_NUM; j++ )
  {
    // test firstto see if on-fly-cmd work
    sprintf (tmp, " %d", j + valArray[0]);

    //sprintf (tmp, " %d", valArray[j] );
    strcat(myCmd.mceCmd,tmp); 
  }
  fprintf(myInfo->fpLog,"%d ",valArray[0]);
  sc2dalib_sendCmd(con,myInfo,&myCmd,mceInxpt,dateTime,status);
}


/**
 * \fn void sc2dalib_sendarrayCmd(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceInxpt, char *dateTime,char *mceCmd, 
 *  int *array,  int howMany, StatusType *status)
 *
 * \brief function
 *  send array cmd
 *
 * \param con      SDSU context structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param mceInxpt  struct mcexml_struct pointer
 * \param dateTime  dateTime string pointer         
 * \param mceCmd    char pointer          
 * \param array    int array, 
 * \param howMany   int for the array
 * \param status StatusType.  given and return
 *
 */
/*+ sc2dalib_sendarrayCmd
 */
void sc2dalib_sendarrayCmd
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
char                  *dateTime,
char                  *mceCmd,
int                   *array,
int                   howMany,
StatusType            *status
)
{
  int          i;
  dasCmdInfo_t myCmd;
  char         tmp[15];

  if (*status != STATUS__OK) return;

  sprintf (myCmd.mceCmd, "%s",mceCmd ); 
  for ( i=0; i<howMany; i++ )
  {
    sprintf ( tmp, " %d", array[i] );
    strcat(myCmd.mceCmd,tmp); 
  }
  sc2dalib_sendCmd(con,myInfo,&myCmd,mceInxpt,dateTime,status);
}  

 

/**
 * \fn void sc2dalib_trksq2fbupdateGain(dasInfoStruct_t *myInfo,
 *  ARRAYSET *arrayset, int ch,  int delta, double *gain, 
 *  int *flag, StatusType *status)
 *
 * \brief function
 *  update gain if to higher or too smaller, 
 *
 * \param myInfo   dasInfoStruct_t pointer
 * \param setup    ARRAYSET structure pointer
 * \param ch      int 
 * \param delta   int sq2fb delta
 * \param gain     double pointer
 * \param flag    int  counter
 * \param status StatusType.  given and return
 *
 */
/*+ sc2dalib_trksq2fbupdateGain
 */
void sc2dalib_trksq2fbupdateGain
(
dasInfoStruct_t      *myInfo,
ARRAYSET              *setup,
int                   ch,
int                   delta,
double                *gain,
int                   *flag,
StatusType            *status
)
{
  int          i,j, aveErr;
  static int   integFlag;
  static int   delta0[COL_NUM];
  static int   positFlag[COL_NUM];
  static int   negatFlag[COL_NUM];
  static int   positNo[COL_NUM];
  static int   negatNo[COL_NUM];
  static int   integErr[COL_NUM][6];
  static double   initGain[COL_NUM];

  if (*status != STATUS__OK) return;
 
  if (myInfo->trkNo==SQ2FB_UPDATE_WAIT)
  { 
    integFlag=0;
    for(i=0;i<COL_NUM;i++)
    {
      flag[i]=0;
      delta0[i]=-1;
      positFlag[i]=negatFlag[i]=0;
      positNo[i]=negatNo[i]=0;
      initGain[i]=gain[i];
      
      for (j=0;j<6; j++ )
       integErr[i][j]=0;
    }
  }
  
  // adapt gain 
  if (flag[ch] < 6 )
  {
    if ( delta > 1 &&  delta >= delta0[ch] )
    {
      positFlag[ch] ++; negatFlag[ch] --; positNo[ch] ++;
    }
    else if ( delta < -1   &&  delta <= delta0[ch] )
    {
      positFlag[ch] --; negatFlag[ch] ++; negatNo[ch] ++;
    }
    
    delta0[ch]=delta;

    if ( positFlag[ch] >2 || negatFlag[ch] >2 )
    {
      // go either direction, so increase gain
      gain[ch] += initGain[ch]*1.5;
      flag[ch]=0;
    }
    else if ( positNo[ch] >2 || negatNo[ch] >2 )
    {
      // oscillate. so, reduce gain
      gain[ch] -= (0.3*gain[ch]);
      flag[ch]=0;
    }
    else
      flag[ch] ++;
  }
  else
    flag[ch]=0;
  
  if (flag[ch]== 0 )
  {
    positFlag[ch]=negatFlag[ch]=0;
    positNo[ch]=negatNo[ch]=0;
  }
  
  // don't do any for this
  integFlag ++;
  if ( integFlag==6 )
  {
     integFlag=0;
     aveErr=0;
     for (j=0; j<6; j++)
       aveErr +=integErr[ch][j];
     
     aveErr /=6;
  }          
}     

 
/**
 * \fn void sc2dalib_readOPT(dasInfoStruct_t *myInfo,
 *     int *optArray, int dim, char *optFile,StatusType *status)
 *
 * \brief function
 *  read sq2fb (sq1bias) optimal Value and store in global array 
 *
 * \param  myInfo   dasInfoStruct_t poiter
 * \param  optArray int pointer, sq2fb (sq1bias) optimal
 * \param  dim      int array dimension
 * \param  optFile  char string for the file
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dalib_readOPT
*/
void sc2dalib_readOPT
(
dasInfoStruct_t *myInfo,
int             *optArray,
int             dim,
char            *optFile,
StatusType      *status
)
{
  FILE   *fp;
  char   tmpfile[FILE_LEN];
  char   input[CMD_LEN];
  size_t cmdlen;
  char   *batchptr;
  int    iscmdLine=1;
  char   delimiter[]=" ", *token;
  char   tmpcmd[CMD_LEN*2];
  int    i,j;

  if (*status != STATUS__OK) return;

  batchptr=tmpcmd;
  fp=NULL;
  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ), optFile );
  if ((fp = fopen(tmpfile, "r")) == NULL)
  { 
    *status=DITS__APP_ERROR;
    ErsRep(0, status,"sc2dalib_readOPT: failed to open optFile: %s",tmpfile);
    return;
  }

  jitDebug(16,"setupFile=%s\n",tmpfile);
  while(1)
  {
    if(  fgets(input,CMD_LEN,fp) !=NULL )
    {
      cmdlen=strlen(input);
      // skip all comments 
      if ( (strchr(input,'#') !=NULL) || (strcmp(input,"\n")==0) )
      {
        jitDebug(16,"commentline=%s\n",input);
      }
      // any line has character <
      else if ( strchr(input,'<') !=NULL)
      {
        j=1;
        // search until '<' change to space
        while(input[cmdlen-j] !='<' )
        {
          input[cmdlen-j]=' ';  j++;
        }
        input[cmdlen-j]=' ';
        
        if ( iscmdLine ==1)
        {
          iscmdLine=0;
          jitDebug(16,"commentWord=%s\n",input);
        }
        else
        {
          batchptr=stpcpy(batchptr, input);
          jitDebug(16,"commentWord=%s\n",input);
        }
      }
      else  // the last values
      {
        input[cmdlen-1]=' ';
        batchptr=stpcpy(batchptr, input);
        batchptr=tmpcmd;
        break;
      }
    }
    else
      break;
  }
  fprintf(myInfo->fpLog,"sc2dalib_readOPT read %s\n",tmpfile);
  jitDebug(16,"cmdline=%s\n",tmpcmd);
  token = strtok (tmpcmd, delimiter);
  optArray[0]=atoi(token);
  jitDebug(16,"optArray(0)=%d\n",optArray[0]);
  /* fprintf(myInfo->fpLog,"%d ",optArray[0]); */
  for ( i=1; i<dim; i++)
  {
    token=strtok(NULL,delimiter);
    optArray[i]=atoi(token);
    /* fprintf(myInfo->fpLog,"%d ",optArray[i]); */
  }
  /* fprintf(myInfo->fpLog,"\n"); */
  fclose(fp);

}


/**
 * \fn void sc2dalib_trksq2fbsavefbZG(dasInfoStruct_t *myInfo,
 *     ARRAYSET *setup, StatusType *status)
 *
 * \brief function
 *  save Zfactor and Gain( unchanged) initFB from trksq2fb for next round 
 *
 * \param  myInfo dasInfoStruct_t poiter
 * \param  setup   ARRAYSET structure pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dalib_trksq2fbsavefbZG
*/
void sc2dalib_trksq2fbsavefbZG
(
dasInfoStruct_t *myInfo,
ARRAYSET        *setup,
StatusType      *status
)
{
  char   tmpfile[FILE_LEN];
  FILE   *fpservo;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ),SQ2OPEN4P_FB_ZG );
  
  if((fpservo=fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_trksq2fbsavefbZG: Error- failed to open file %s", tmpfile); 
      return;
    }

  fprintf(fpservo,"\n# INITFB_N  initial sq2fb lock Value \n");
  fprintf(fpservo,"#  non-modulation ones, set intFB to (%d)\n",
           NONMODULATION);
  sc2dalib_savelckPoints1(setup,fpservo,"intFB= ",setup->initFB,0,status);

  fprintf(fpservo,"\n# ZFACTOR_N: ssaDAOut Value \n");
  sc2dalib_savelckPoints1(setup,fpservo,"zFactor= ",setup->zfact,0,status);

  fprintf(fpservo,"\n# GAIN_N for all channels\n");
  sc2dalib_savelckPoints1(setup,fpservo,"Gain= ",setup->gain,1,status);

  fprintf(fpservo,"\n# FLUXPERIOD from sq2 only for sq2fb servo \n");
  sc2dalib_savelckPoints1(setup,fpservo,"fluxPeriod= ",setup->fluxPeriod,0,status);  

  fclose(fpservo);
}



/**
 * \fn void sc2dalib_savelckPoints1(ARRAYSET *setup, FILE *fp, char *name,
 *  void *paramarray, int flag,StatusType *status)
 *
 * \brief function:
 *  save servo lock points for each channel
 *
 * \param  setup   ARRAYSET structure pointer
 * \param  fp      FILE pointer for result file
 * \param  name       char point for the param name
 * \param  paramarray void pointer for setup-parameter array
 * \param  flag      int  0: paramarray is *int, 1: *double
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalib_savelckPoints1
*/
void sc2dalib_savelckPoints1
(
ARRAYSET   *setup,
FILE       *fp,
char       *name,
void       *paramarray,
int        flag,
StatusType *status
)
{
  int   ch,i;
  int   *intArray;
  double *floatArray;

  if (*status != STATUS__OK) return;

  intArray=NULL;
  floatArray=NULL;

  if( flag==0)       
    intArray=(int*)paramarray;
  else
    floatArray=(double*)paramarray; 

  for (i=0;i<4;i++)
  {
    fprintf(fp,"#RC%d\n",i+1);
    fprintf(fp,"%s",name);
    for (ch=i*ONECARDCHS; ch<(i+1)*ONECARDCHS;ch++)
    {
      if( flag==0)
      {
        fprintf(fp,"%7d ", intArray[ch]);
      }
      else
      {
        fprintf(fp,"%7.4f ", floatArray[ch]);
      }
    }
    fprintf(fp,"\n");
  }  
}



/**
 * \fn void sc2dalib_trksq2fbsavesq2fbOPT(dasInfoStruct_t *myInfo,
 *     int sq2fbOPT, StatusType *status)
 *
 * \brief function
 *  save updated sq2fb optimal value 
 *
 * \param  myInfo   dasInfoStruct_t poiter
 * \param  sq2fbOpt int pointer
 * \param  status StatusType.  given and return
 *
 */
/*+ sc2dalib_trksq2fbsavesq2fbOPT
*/
void sc2dalib_trksq2fbsavesq2fbOPT
(
dasInfoStruct_t *myInfo,
int             *sq2fbOpt,
StatusType      *status
)
{
  char   tmpfile[FILE_LEN];
  FILE   *fp;

  if (*status != STATUS__OK) return;

  sprintf (tmpfile, "%s/sq2fboptimal.txt", getenv ( "SC2SCRATCH" ) );

  if((fp=fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_trksq2fbsavesq2fbOPT: Error- failed to open file %s", tmpfile); 
      return;
    }

  fprintf(fp,"#this is sq2fb optimal value, updated after sq2fbtracking\n");
  fprintf(fp,"wb bc2 flux_fb < \n");
  sc2dalibsetup_savesq2fboptPts(fp,"sq2fbopt= ",sq2fbOpt,1,status);
  fprintf(fp,"\n");
  fclose(fp);
}


// =======sc2dalib_v*******
//====================//
/**
 * \fn void sc2dalib_versionInit(SDSU_CONTEXT *con,dasInfoStruct_t *myInfo, 
 *  char *dateTime,StatusType *status)
 *
 * \brief fountion
 *  initial set up for VERSION action 
 *
 * \param con      SDSU context structure
 * \param myInfo  dasInfo structure pointer
 * \param dateTime  string pointer to dateTime string
 * \param status  StatusType     
 *
 */
/*+ sc2dalib_versionInit
*/
void sc2dalib_versionInit
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
char                  *dateTime,
StatusType            *status
)
{
  long               in_sequence;

  if (*status != STATUS__OK) return;

  // Update debug flag, in case it has changed 
  utils_update_debug(con,myInfo, status);
  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0, status,"sc2dalib_versionInit: %s has not completed",
       seqStatus[myInfo->actionFlag]);
    return;
  }
  if ( myInfo->logfileFlag==0)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0, status,"sc2dalib_versionInit: the logfile is not valid");
    return;
  }
  myInfo->actionFlag=VERSIONACTION;
  // the logfile already opened in sc2dalib__Init
  fprintf(myInfo->fpLog,"\n<%s> CMD for sc2dalib__Version\n",dateTime);
}


/**
 * \fn void sc2dalib_versionInfo(SDSU_CONTEXT *localcon, 
 *   dasInfoStruct_t *myInfo, StatusType *status)
 *
 * \brief function
 *  get version infomation, display and save.
 *
 * \param localcon   SDSU_CONTEXT pointer
 * \param myInfo    dasInfoStaruct pointer
 * \param status     StatusType.  given and return
 *
 */
/*+ sc2dalib_versionInfo 
*/
void sc2dalib_versionInfo
(
SDSU_CONTEXT    *localcon,
dasInfoStruct_t *myInfo,
StatusType      *status
)
{
  int            rval,i;
  long           memAddr;
  long           pciVer[2];
  char           memTypeInt,versionType[10];
  PCI_CMD pcicmd;
  char    localmsg[FILE_LEN];

  if (*status != STATUS__OK) return;

  //now read PCI version, stored in X:3,4 
  memTypeInt='X'; 
  memAddr=3;
  for (i=0;i<2;i++)
  {
    memAddr=memAddr+i;
    jitDebug(2,"dasVerion: X:%ld is read \n",memAddr);
    rval=sdsu_set_pcicmd(&pcicmd,"READ", memTypeInt,memAddr,&pciVer[i]);
    if( rval!=SDSU_OK )
    { 
      *status = DITS__APP_ERROR;
      ErsRep(0,status,"sc2dalib_versionInfo: sdsu_set_pcicmd failed");   
      return;
    }
    rval=sdsu_command_pci(localcon,&pcicmd);
    if( rval !=SDSU_OK)
    {
     *status=DITS__APP_ERROR; 
      sc2dalib_pcierrRep(&pcicmd,"READ",status);        
      ErsRep(0,status,"sc2dalib_versionInfo: sdsu_command_pci failed"); 
      return;
    }
    pciVer[i]=pcicmd.arg2; 
  }  
  memTypeInt=(char)(pciVer[0]>>16);
  if( memTypeInt=='A')
     strcpy(versionType,"A");
  else if( memTypeInt=='B')
     strcpy(versionType,"B");
  else if( memTypeInt=='X')
     strcpy(versionType,"X");
  else
  {
     utils_msg(myInfo,
        "sc2dalib_versionInfo: wrong PCI version Letter, it shall be A or B or X",
        "",USE_MSGOUT,status);
     strcpy(versionType,"NOTVALID --");
  }
  memAddr=(pciVer[0] & 0x0000FFFF);
  utils_msg(myInfo,
                    "===========   VERSION INFORMATION   ===========",
                    "",USE_MSGOUT,status);
  utils_msg(myInfo,"sc2dalib_versionInfo: driver version= %s",
                    myInfo->drvVersion,USE_MSGOUT,status);
  utils_msg(myInfo,"sc2dalib_versionInfo: DAS version=%s",
                     DAS_VERSION,0,status);
  sprintf(localmsg," PCI DSP version=%s%04lX (X:03)",versionType,memAddr);
  utils_msg(myInfo,"sc2dalib_versionInfo:%s",localmsg,USE_MSGOUT,status);
  sprintf(localmsg," PCI DSP date=%06lX (X:04)",pciVer[1]);
  utils_msg(myInfo,"sc2dalib_versionInfo:%s",localmsg,USE_MSGOUT,status);
  utils_msg(myInfo,
                    "===============================================",
                    "",USE_MSGOUT,status);
}


/**
 * \fn void sc2dalib_callorphanHandler(StatusType *status)
 *
 * \brief function
 *   stop orphan messsage
 * \param status     StatusType.  given and return
 *
 */
/*+ sc2dalib_callorphanHandler
*/
void sc2dalib_callorphanHandler
(
StatusType      *status
)
{
}



//===============  now place sq1 optimization into dramatask action  =====
//
//================ X. Gao 20091022  ======================================

/**
 * \fn void sc2dalib_sq1optptsInit(dasInfoStruct_t *myInfo, int *option,StatusType *status)
 *
 * \brief function
 *  get args for drama action and open file
 *
 * \param myInfo    dasInfoStruct_t pointer
 * \param option        int pointer
 * \param status    StatusType     
 *
 */
/*+ sc2dalib_sq1optptsInit
*/
void sc2dalib_sq1optptsInit
(
dasInfoStruct_t       *myInfo,
int                    *option,
StatusType            *status
)
{
  long       in_sequence, range[] = { 0, 3};
  long       flag;
  char       tmp[FILE_LEN];
  SdsIdType argId;

  if (!StatusOkP(status)) return;

  SdpGeti("IN_SEQUENCE", &in_sequence, status);
  if(in_sequence != DA_MCE_NONE)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"_sq1optptsInit: %s has not completed",seqStatus[dasInfo.actionFlag]);
    return;
  } 
  myInfo->actionFlag=SQ1OPTACTION;

  argId = DitsGetArgument();
  jitArgGetS(argId, "SQ1OPTPTS", 1, NULL, "sq1opt-points", 0, FILE_LEN,
	     tmp, NULL, status);
  sprintf( myInfo->batchFile,"%s",tmp);
  //print ?
  jitArgGetS(argId, "RESULT_FILE", 2, NULL, "sq1opt-xxx", 0, FILE_LEN,
	     tmp, NULL, status);
  sprintf(myInfo->dataFile,"%s",tmp);
  jitArgGetI(argId, "FLAG", 3, range, 0, 0, &flag, status);
  if ( !StatusOkP(status) )
  {
    ErsRep(0,status,"_sq1optptsInit: failed to get variables.");
    return;
  }
  *option=flag;
  jitDebug(2,"_sq1optptsInit:sq1optPoint=%s\n",myInfo->batchFile);
  utils_fclose(&(myInfo->fpBatch));
  if ((myInfo->fpBatch = fopen(myInfo->batchFile,"r")) == NULL)
  {
    *status = DITS__APP_ERROR;
    ErsRep(0,status,"_sq1optptsInit: failed to open %s",myInfo->batchFile);
    return;
  } 
}



/**
 * \fn void sc2dalib_sq1optptsreadInput(FILE *fp, int *bias,int *fdbk,
 *  int *refbias, int *reffdbk,StatusType *status)
 *
 * \brief funcation: read sq1 optimal input file
 *
 * \param  fp        FILE pointer
 * \param  bias      int pointer
 * \param  fbdk      int pointer
 * \param  refbias   int pointer
 * \param  reffdbk   int pointer
 * \param status  StatusType pointer.  given and return
 *
 */
/*+ sc2dalib_sq1optptsreadInput
*/
void sc2dalib_sq1optptsreadInput
(
FILE       *fp,
int        *bias,
int        *fdbk,
int        *refbias,
int        *reffdbk,
StatusType *status
)
{
  char       delimiters[]= " =\n";
  char       readLine[FILE_LEN], *token;   
  int        rowNo=0,i;
  int        rowMap[32];

  if (!StatusOkP(status)) return;

  for (i=0;i<ROW_NUM;i++)
     rowMap[i]=0;

  while(1)
  {
    if(  fgets(readLine,FILE_LEN,fp) !=NULL )
    {
      token = strtok (readLine, delimiters);
      if (strcmp("sel_row",token)==0)
      {
        token = strtok (NULL, delimiters);   
        rowNo= atoi (token);
        rowMap[rowNo]=1;
        //printf("\nrow =%d ",rowNo);
      }  
      else if (strcmp("sq1bias",token)==0)
      {
        sc2dalib_sq1biasfdbkPoints(delimiters,"sq1bias",bias,rowNo,VAL__BADI,status);
        if (!StatusOkP(status)) 
          return;
      }  
      else if (strcmp("sq2fdbk",token)==0)
      {
        sc2dalib_sq1biasfdbkPoints(delimiters,"sq2fdbk",fdbk,rowNo,VAL__BADI,status);
        if (!StatusOkP(status)) 
          return;
      }
      else if (strcmp("sq1refbias",token)==0)
      {
        sc2dalib_sq1biasfdbkPoints(delimiters,"sq1refbias",refbias,rowNo,VAL__BADI,status);
        if (!StatusOkP(status)) 
          return;
      }
      else if (strcmp("sq2reffdbk",token)==0)
      {
        sc2dalib_sq1biasfdbkPoints(delimiters,"sq2reffdbk",reffdbk,rowNo,VAL__BADI,status); 
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
    if(rowMap[i]==0)
    {
      MsgOut(status,"_sq1optptsreadInput:row_%d not used, place VAL__BADI",i);
      sc2dalib_sq1biasfdbk4notUsed("sq1bias",bias,i,VAL__BADI,status); 
      sc2dalib_sq1biasfdbk4notUsed("sq2fdbk",fdbk,i,VAL__BADI,status); 
      sc2dalib_sq1biasfdbk4notUsed("sq1refbias",refbias,i,VAL__BADI,status);
      sc2dalib_sq1biasfdbk4notUsed("sq2reffdbk",reffdbk,i,VAL__BADI,status); 
    }
  }
}



/**
 * \fn void sc2dalib_sq1biasfdbkPoints( char* delimiters,char *name,
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
/*+ sc2dalib_sq1biasfdbkPoints
*/
void sc2dalib_sq1biasfdbkPoints
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

  if (!StatusOkP(status)) return;
 
  arrayInx = row*COL_NUM;
  //printf("\n%s ",name);
  while(readCh)
  {
    if( (token=strtok(NULL,delimiters))==NULL)      
    {
      *status=DITS__APP_ERROR;
      ErsRep(0,status,"_sq1biasfdbkPoints: %s's dataNo < %d",name,COL_NUM);
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
 * \fn void sc2dalib_sq1biasfdbk4notUsed (char *name,
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
/*+ sc2dalib_sq1biasfdbk4notUsed
*/
void sc2dalib_sq1biasfdbk4notUsed
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

  if (!StatusOkP(status)) return;
 
  arrayInx = row*COL_NUM;

  for (ch=0;ch<COL_NUM;ch++)
      paramarray[arrayInx+ch]=Val;
}




/**
 * \fn void sc2dalib_sq1optwritesq1bsq2fOut(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *  int *sq1bias, int *sq2feedback, int *sq1refbias,  int *sq2reffeedback, int flag, StatusType *status)
 *
 * \brief function:
 *  save sq1 bias sq2 fdbk points to a text file 
 *
 * \param myInfo  dasInfo structure pointer
 * \param  setup   ARRAYSET structure pointer
 * \param  sq1bias     int pointer 
 * \param  sq2feedback int pointer
 * \param  sq1refbias  int pointer
 * \param  sq2reffeedback int pointer
 * \param   flag      int 
 * \param status      StatusType pointer.  given and return
 *
 */
/*+ sc2dalib_sq1optwritesq1bsq2fOut
*/
void sc2dalib_sq1optwritesq1bsq2fOut
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
int              *sq1bias,       //  bias for each SQ1 (given) 
int              *sq2feedback,    //SQ2 feedback at lock point for sq1bias(given) 
int              *sq1refbias,     // reference bias for each SQ1 (given) 
int              *sq2reffeedback, // SQ2 feedback at lock point for sq1refbias (given)
int               flag,            // = 1 allow sq2fb< 0;  =0 set sq2fb<0 pixel =BAD
StatusType *status
)
{
  FILE *fd, *fd1;    /* file descriptor for output */
  int col;           /* column index */
  int row;           /* row index */
  int sqnum;         /* sequential SQ1 count */
  char tmp[FILE_LEN];
  char file[FILE_LEN];
 
  if (!StatusOkP(status)) return;
 
  sprintf(file, "%s-pixel",myInfo->dataFile);
  jitDebug(2,"_sq1optwritesq1bsq2fOut: %s\n",file);

  if((fd = fopen( file, "w" )) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_sq1optwritesq1bsq2fOut: Error- failed to open file %s", file); 
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
      fprintf (fd, "%d  %d   %d  %d  %d  %d\n", row+1, col+1, sq1bias[sqnum], 
               sq2feedback[sqnum], sq1refbias[sqnum], sq2reffeedback[sqnum] );
    }
  }
  fclose (fd ); 

  sprintf(tmp,"%s-sq1biassq2fb",file);
  if((fd1 = fopen( tmp, "w" )) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_sq1optwritesq1bsq2fOut: 2 Error- failed to open file %s", tmp); 
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
 * \fn void sc2dalib_chksinglePixel(ARRAYSET *setup,  int *sq1comp, 
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
/*+ sc2dalib_chksinglePixel
*/
void sc2dalib_chksinglePixel
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

  if (!StatusOkP(status)) return;

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
    MsgOut(status,"_chksinglePixel:there is only one pixel(%d,%d), skip optimal",
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
 * \fn void sc2dalib_sq2fbeachcolOutlier(ARRAYSET *setup, OPT_STRUCT *sq2fbHist,
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
/*+ sc2dalib_sq2fbeachcolOutlier
*/
void sc2dalib_sq2fbeachcolOutlier
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
  char errmsg[FILE_LEN];

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
    sprintf(errmsg, "_sq2fbeachcolOutliers: MUST NOT LET min (%d) <0 \n",min ); 
    ErsRep(0,status,errmsg);
    return; 
  }

  if (bin ==0) 
  {  
    // *status=DITS__APP_ERROR;
    sprintf(errmsg, "_sq2fbeachcolOutliers: bin 0 max(%d) min(%d) \n",max,min ); 
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
        sc2dalib_sq2fbhalfFlux(setup,data,i,peakVal,halfFlux,
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
        sc2dalib_sq2fbhalfFlux(setup,data,i,peakVal,halfFlux,
	         &sq2fbHist->changeInx[pixel],status);
    }
  }
  else
  {
    for (i=0; i<ROW_NUM; i++)
    { 
      pixel=i*COL_NUM +col;
      if ( abs(binInx[i]-whichBin) > SQ2FDBK_WIDTH )
        sc2dalib_sq2fbhalfFlux(setup,data,i,peakVal,halfFlux,
	         &sq2fbHist->changeInx[pixel],status);
    }
  }
  return;
}


/**
 * \fn void sc2dalib_sq2fbhalfFlux(ARRAYSET *setup,int *data, int which, int peakVal,
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
/*+ sc2dalib_sq2fbhalfFlux
*/
void sc2dalib_sq2fbhalfFlux
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
  // if setup->select[12] set <0, then we don't do half-fold 
  if ( abs(apart-halfFlux) < setup->slopSelect[12] )
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
 * \fn void sc2dalib_sq2fbOutlier(dasInfoStruct_t  *myInfo,   ARRAYSET *setup,
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
/*+ sc2dalib_sq2fbOutlier
*/
void sc2dalib_sq2fbOutlier
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

  if((fp = fopen( filename, "w" )) == NULL )
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status, "sc2dalib_sq2fbOutlier: Error- failed to open file %s", filename); 
    return;
  }
  sprintf(tmpfile,"%s-freq",filename);
  if((fpFreq = fopen( tmpfile, "w" )) == NULL )
  {
    *status = DITS__APP_ERROR;
    ErsRep (0, status, "sc2dalib_sq2fbOutlier: 2 Error- failed to open file %s", tmpfile); 
    return;
  }

  fprintf(fp,"sc2dalib_sq2fbOutliers: Histogram for each col\n");
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
    sc2dalib_sq2fbeachcolOutlier(setup,&sq2fbHist,data,col,status);
    if ( *status != STATUS__OK ) 
    {
      ErsRep(0,status,"sc2dalibsetup_sq2fbOutlier: col=%d, failed in sc2dalibsetup_sq2fbeachcolOutlier\n",col);
      fclose(fp);  fclose(fpFreq);
      return ;
    }

    fprintf(fp,"col_%3d  bin[%7d]  pealVal[~%7d] peakInx[%4d] min[%7d] max[%7d]\n",
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
          fprintf(fp,"row_%3d [%7d] is outlier, set BAD \n",row,sq2feedback[pixel]);
          sq2feedback[pixel]=VAL__BADI;
        }
        else        
        {
          fprintf(fp,"row_%3d [%7d] are all the same, keep \n",row,sq2feedback[pixel]);
        }
      }
      else if ( data[row]==0  &&  badRow[row]==1)
      {
        fprintf(fp,"row_%3d [%7d] is outlier, set BAD \n",row,sq2feedback[pixel]);
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
 * \fn void sc2dalib_sq1biassq2fdbkClipmean(int *sq1bias, int *sq1boptimal, int *sq2fdbk,
 *  int *sq2fdbk, int*sq2fdbkoptimal, int colNo, int rowNo, StatusType *status)
 *
 * \brief function:
 *  find sq1bias and sq2fdbk using clipmean 
 *  only used with drama, replace sc2dalibsetup_sq1biassq2fdbkClipmean 
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
/*+ sc2dalib_sq1biassq2fdbkClipmean
*/
void sc2dalib_sq1biassq2fdbkClipmean
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
     //Changing for new sigma was 3.0
     sc2math_clipmean(2.0, colNo, data, &meanVal,status);
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
     //Changing for new sigma was 3.0
     sc2math_clipmean(2.0, rowNo, data, &meanVal,status);
     sq2foptimal[i]=(int)meanVal;
#ifdef DEB_OPT_COL
     if ( i== DEB_OPT_COL) 
       printf("\nsq2foptimal[0]=%d\n",sq2foptimal[0]);

     if( i==0) printf("sq2foptimal=:");
     printf("%d ",sq2foptimal[i]);
#endif
  }
  //printf("\n");
}


/**
 * \fn void sc2dalib_sq1optAlgrm(dasInfoStruct_t *myInfo, ARRAYSET *setup, int *sq1bcomp,
 *  int *sq2fcomp, int *sq1bias, int *sq2feedback, int *sq1refbias,  int *sq2reffeedback, 
 *  int *bqual, int *fqual int *squal, double *sq1scale, int flag, StatusType *status)
 *
 * \brief function:
 *  find optimal sq1bias and sq2feedback 
 *
 * \param  myInfo         dasInfo structure pointer
 * \param  setup          ARRAYSET structure pointer
 * \param  sq1bcomp       int pointer   compromise SQ1 bias for each row    (return)
 * \param  sq2fcomp       int pointer   compromise SQ2 feedback for each column (return)
 * \param  sq1bias        int pointer   bias for each SQ1 (given) 
 * \param  sq2feedback    int pointer   SQ2 feedback at lock point for sq1bias
 * \param  sq1refbias     int pointer   reference bias for each SQ1 
 * \param  sq2reffeedback int pointer   SQ2 feedback at lock point for sq1refbias (given)
 * \param  bqual          int pointer   quality for each SQ1 bias row 
 * \param  fqual          int pointer   quality for each SQ2 feedback column 
 * \param  squal          int pointer   quality for each SQ1 
 * \param  sq1scale       double pointer   scale factor relating SQ1 bias and SQ2 feedback 
 * \param  flag           int    1 takeout outlier of sq2fb set outlier=BAD;  flag=0 don't take out
 * \param  status          StatusType pointer.  given and return
 *
 */
/*+ sc2dalib_sq1optAlgrm
*/
void sc2dalib_sq1optAlgrm
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
int              *sq1bcomp,
int              *sq2fcomp,
int              *sq1bias,         
int              *sq2feedback,    
int              *sq1refbias,     
int              *sq2reffeedback, 
int              *bqual,         
int              *fqual,         
int              *squal,         
double           *sq1scale,      
int              flag,           
StatusType *status
)
{
  int  i,singleFlag=0;
  char outfileXml[FILE_LEN];
  char tmpfile[FILE_LEN];
  int  sq1biascp[COL_NUM*ROW_NUM];
  int  sq2feedbackcp[COL_NUM*ROW_NUM]; 

  if (!StatusOkP(status)) return;

  sprintf(outfileXml,"%s-xml",myInfo->dataFile);

  // check if it is single pixel, if it is then pass sq1bias and sq2fdbk to sq1comp, sq2fcomp 
  jitDebug(2,"_sq1optAlgrm: _chksinglePixel\n");
  sc2dalib_chksinglePixel(setup,sq1bcomp, sq2fcomp,sq1bias,sq2feedback, &singleFlag, status);
  if (singleFlag ==0)
  {
    if (flag==1)
    {
      jitDebug(2,"_sq1optAlgrm: sq2fbOutliers\n");
      sprintf(tmpfile, "%s-outlier",myInfo->dataFile);
      jitDebug(2,"_sq1optAlgrm::4nextStep: %s\n",tmpfile);
      // find outliers before optimal routine  and save into dasInfo.datafile-outlier
      sc2dalib_sq2fbOutlier( myInfo,setup,tmpfile,sq2feedback,sq2reffeedback,status);
      if ( !StatusOkP(status) )
        return;
    }
    // make copies of sq1bias and sq2feedback, it seems that sc2sqopt_putSQ1 changes them
    memcpy(sq1biascp, sq1bias,sizeof(int)*COL_NUM*ROW_NUM);
    memcpy(sq2feedbackcp, sq2feedback,sizeof(int)*COL_NUM*ROW_NUM);
  
    jitDebug(2,"_sq1optAlgrm: acllsc2sqopt_putSQ1\n");
    sc2sqopt_putSQ1 ( outfileXml, sq1bias, sq2feedback,sq1refbias, sq2reffeedback, status );
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"_sq1optAlgrm: sc2sqopt_putSQ1 bfailed");
      return;
    }

    jitDebug(2,"_sq1optAlgrm: call sc2sqopt_getSq1scales\n");
    //change from sc2sqopt_getquals to sc2sqopt_getSQ1scales  07-Nov-2007 x.gao
    sc2sqopt_getSQ1scales(sq1bias,sq2feedback,sq1refbias,sq2reffeedback,sq1scale,squal,fqual,bqual,status );
    if ( !StatusOkP(status) )
    {
      ErsRep(0,status,"_sq1optAlgrm: sc2sqopt_getSQ1scales failed");
      return;
    }
   
    // write out
    sprintf(tmpfile, "%s-pixel",myInfo->dataFile);
    //  bad pixel if squal=1 
    sc2dalib_writepixelInx(myInfo,setup,tmpfile,squal,0,status);
    utils_fclose(&(myInfo->fpData));
    if((myInfo->fpData = fopen( myInfo->dataFile, "w" )) == NULL )
      {
	*status = DITS__APP_ERROR;
	ErsRep (0, status, "sc2dalib_sq1optAlgrm: Error- failed to open file %s", myInfo->dataFile); 
	return;
  }
    sc2dalib_writequalOut(myInfo->fpData,sq1scale,squal,fqual,bqual,status );

    /* Perform the fit */
    MsgOut(status,"_sq1optAlgrm: ...sc2sqopt_fit\n");
    sc2sqopt_fit (sq1bias,sq2feedback,sq1scale,squal,fqual,bqual,sq1bcomp,sq2fcomp,status);

    /* Write the solution, status no change */
    sc2dalib_writecompOut(myInfo->fpData,sq1bcomp,sq2fcomp,status );

    // for whatever reason, we will do the manual set
    if( status == STATUS__OK)
      MsgOut(status,"_sq1optAlgrm:we are going to get AVER for those rejected by optimal");
    else
    {
      *status = STATUS__OK;
       MsgOut(status,"_sq1optAlgrm: optimal failed, set sq1bcomp/sq2fcomp=BAD, going to get AVER");
      for ( i=0; i<ROW_NUM; i++ )
        sq1bcomp[i]=VAL__BADI;
      for ( i=0; i<COL_NUM; i++ )
        sq2fcomp[i]=VAL__BADI;
    }
    // no w manually chang the BAD ones
    fprintf (myInfo->fpData, "sq1bias manually changed\n" );
    for ( i=0; i<ROW_NUM; i++ )
    {
      if (sq1bcomp[i]==VAL__BADI)
      {
        if (setup->rowMask[i] !=0) // not masked out
          sc2dalib_sc2dafindalter4Opt(myInfo,setup,sq1bcomp,sq1biascp,1,i,status);
        fprintf ( myInfo->fpData, "row_%d %6d rowMak=%d \n",i,sq1bcomp[i],setup->rowMask[i] );
      }
    }
    fprintf ( myInfo->fpData, "sq2fdbk manually changed\n" );
    for ( i=0; i<COL_NUM; i++ )
    {
      if (sq2fcomp[i]==VAL__BADI)
      {
        if( setup->colMask[i] !=0)
          sc2dalib_sc2dafindalter4Opt(myInfo,setup,sq2fcomp,sq2feedbackcp,0,i,status);
        fprintf ( myInfo->fpData, "col_%d %6d  colMask=%d \n",i, sq2fcomp[i],setup->colMask[i] );
      }
    }
    fprintf ( myInfo->fpData, "\n" );
    utils_fclose(&(myInfo->fpData));
  }
}




/**
 * \fn void sc2dalib_writepixelInx(dasInfoStruct_t  *myInfo,   
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
/*+ sc2dalib_writepixelInx
*/
void sc2dalib_writepixelInx
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
  
  if((fp = fopen( tmp, "w" )) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_writepixelInx: Error- failed to open file %s", tmp); 
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
	  /* The calculated integrator gain was too much so set it to zero */
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
  fprintf ( fp, "\nTotal badOnes= %d  ratio=%f\n",badOne,(float)badOne/totalPixel);
  fclose ( fp );   
}



/**
 * \fn void sc2dalib_writequalOut( FILE *fd, double *sq1scale, int *squal,
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
/*+ writesq1biasfdbkOut
*/
void sc2dalib_writequalOut
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

  if (!StatusOkP(status)) return;

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
 * \fn void sc2dalib_writecompOut( FILE *fd, int *sq1bcomp,
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
/*+ writesq1biasfdbkOut
*/
void sc2dalib_writecompOut
(
FILE       *fd,
int        *sq1bcomp,
int        *sq2fcomp,
StatusType *status
)
{
  int   i;

  if (!StatusOkP(status)) return;

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
 * \fn void sc2dalib_sc2dafindalter4Opt(dasInfoStruct_t  *myInfo,    
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
/*+ sc2dalib_sc2dafindalter4Opt
*/
void sc2dalib_sc2dafindalter4Opt
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
  int  pixel, aver, i;
  int  data;

  if (!StatusOkP(status)) return;
  
  if(flag==1)  // row
  {
    data=0, aver=0;
    for (i=0;i<COL_NUM;i++)
    {
      pixel=which*COL_NUM+i;

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
        MsgOut(status,"row_%d: AVER value %d over the MCE range, set BAD", which,optArray[which]);
        optArray[which]=VAL__BADI;
      }
    }
    else
    {
      MsgOut(status,"row_%d: all BAD_VAL set BAD", which); 
      optArray[which]=VAL__BADI;
    }
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
 * \fn void sc2dalib_sq1optsave2File(dasInfoStruct_t *myInfo, ARRAYSET *setup,
 *  int *sq1bcomp, int *sq2fcomp, StatusType *status)
 *
 * \brief function:
 *  save optimal sq1bias and sq2feedback 
 *
 * \param  myInfo         dasInfo structure pointer
 * \param  setup          ARRAYSET structure pointer
 * \param  sq1bcomp       int pointer   compromise SQ1 bias for each row   
 * \param  sq2fcomp       int pointer   compromise SQ2 feedback for each column 
 * \param  status          StatusType pointer.  given and return
 *
 */
/*+ sc2dalib_sq1optsave2File
*/
void sc2dalib_sq1optsave2File
(
dasInfoStruct_t  *myInfo,    
ARRAYSET         *setup,
int              *sq1bcomp,
int              *sq2fcomp,
StatusType *status
)
{
  FILE *fpopen;
  char  tmpfile[FILE_LEN];

  if (!StatusOkP(status)) return;

  // save sq1bcomp to sq1biasoptimal.txt
  sprintf (tmpfile, "%s/%s", getenv ("SC2SCRATCH"),SQ1BIASOPTFILE );

  if((fpopen=fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_sq1optsave2File: Error- failed to open file %s", tmpfile); 
      return;
    }
  fprintf(fpopen,"# BIASLCK_N  sq1bias Val for optimal, set 0 for VAL__BADI \n");
  fprintf ( fpopen, "wb ac on_bias <\n" );  
  sc2dalib_sq1optsavesq1bComp(fpopen,"biasoptPt= ",sq1bcomp,1,status);
  fclose(fpopen);

  // save sq2fcomp to sq2fboptimal.txt
  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ), SQ2FBOPTFILE );
  if((fpopen=fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_sq1optsave2File 2: Error- failed to open file %s", tmpfile); 
      return;
    }
  fprintf(fpopen, "#this is sq2fb optimal val \n");
  fprintf(fpopen,"wb bc2 flux_fb < \n");
  sc2dalib_sq1optsavesq2fComp(fpopen,"sq2fbopt= ",sq2fcomp,1,status);
  fprintf(fpopen,"\n");
  fclose(fpopen);

  // save sq1bcomp sq2fvomp to 
  sprintf (tmpfile, "%s/%s", getenv ( "SC2SCRATCH" ), OPTIMAL_SQ1B_SQ2FB );
  if((fpopen=fopen(tmpfile,"w")) == NULL )
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "sc2dalib_sq1optsave2File 3: Error- failed to open file %s", tmpfile); 
      return;
    }

  fprintf(fpopen,"# BIASLCK_N  sq1bias Val for optimal, set 0 for VAL__BADI\n");
  sc2dalib_sq1optsavesq1bComp(fpopen,"biasoptPt= ",sq1bcomp,0,status);

  fprintf(fpopen, "\n\n#this is sq2fdbk optimal val, BAD =0\n");
  utils_array_to_file(fpopen,"sq2fdbkopt= ",(void *)sq2fcomp, INT_NUMBER,status);

 
  fclose(fpopen);
}


/**
 * \fn void sc2dalib_sq1optsavesq1bComp(FILE *fp, char *name,
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
/*+ sc2dalib_sq1optsavesq1bComp
*/
void sc2dalib_sq1optsavesq1bComp
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
    fprintf(fp,"    %s ",name);
    for (j=0;j<ROW_NUM;j++)
    {
      // trap also for those >MAX_SQ1BIAS, asked Dan  21-09-07
      if(  paramarray[j]==VAL__BADI ||  paramarray[j] >= MAX_SQ1BIAS )
       fprintf(fp,"0 ");
      else 
       fprintf(fp,"%7d ", paramarray[j]);
    }
    fprintf(fp,"\n");
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
 * \fn void sc2dalib_sq1optsavesq2fComp(FILE *fp, char *name,
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
/*+ sc2dalib_sq1optsavesq2fComp
*/
void sc2dalib_sq1optsavesq2fComp
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


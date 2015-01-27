
#include "sc2da_mce.h"

/**
 * \fn void mce_read(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceinxPtr, char *cmd,  int *val, int howMany,
 *  StatusType *status)
 *
 * \brief function: 
 *  read back MCE parameters Val 
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param mceinxPtr  mcexml_struct pointer
 * \param cmd        char pointer for command
 * \param  val      int pointer for returned Val 
 * \param howMany      int howmany to return 
 * \param status    StatusType pointer.  given and return
 *
 * if  *status != STATUS__OK, report error.
 * 
 * Used to be sc2dalib_readmceVal
 */
/*+ mce_read
*/
void mce_read
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceinxPtr,
char                  *cmd,
int                   *val,
int                   howMany,
StatusType            *status
)
{
  int           i;
  char         errmsg[FILE_LEN];
  dasCmdInfo_t  mycmdInfo;

  if (*status != STATUS__OK) return;

  mycmdInfo.cmdBufSize=MCEBLK_SIZE;
  strcpy(mycmdInfo.mceCmd,cmd);  
  mcexml_translate(mycmdInfo.mceCmd, mceinxPtr,mycmdInfo.cmdBufSize,
                          mycmdInfo.cmdBuf,status);
  if ( *status ==DITS__APP_ERROR)
  {
    sprintf(errmsg,"mce_read: failed to call mcexml_translate");  
    ErsRep (0, status, errmsg);
    return;
  }
  if( (*status=sdsu_command_mce (con,mycmdInfo.cmdBuf,&mycmdInfo.reply)) <0)
  { 
     strcpy(mycmdInfo.mceCmd,cmd);
     mce_error_reply(myInfo,&mycmdInfo,status); 
     sprintf(errmsg,"mce_read: failed to call sdsu_command_mce");
     ErsRep (0, status, errmsg);
     return;
  }
  if( mycmdInfo.reply.status == MCE_RBOK )
  {
    if ( howMany==1) 
      *val=(int)mycmdInfo.reply.data[1];
    else
    {
      for (i=0; i<howMany; i++)
        val[i]=(int)mycmdInfo.reply.data[i+1];
    }
  }
  else
  {
    *status=DITS__APP_ERROR;
    sprintf(errmsg,"mce_read: MCE reply(%#0lx) != RBOK",
            mycmdInfo.reply.status );  
    ErsRep (0, status, errmsg);
    return;
  }    
}


/**
 * \fn void mce_set(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceinxPtr, char *cmd,  StatusType *status)
 *
 * \brief function: 
 *  set MCE parameters Val 
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param mceinxPtr  mcexml_struct pointer
 * \param cmd        char pointer for command
 * \param status    StatusType pointer.  given and return
 *
 * Used to be sc2dalib_setmceVal
 */
/*+ mce_set
*/
void mce_set
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceinxPtr,
char                  *cmd,
StatusType            *status
)
{
  char         errmsg[FILE_LEN];
  dasCmdInfo_t mycmdInfo;

  if(*status != STATUS__OK) return;

  mycmdInfo.cmdBufSize=MCEBLK_SIZE;
  strcpy(mycmdInfo.mceCmd,cmd);  
  mcexml_translate(mycmdInfo.mceCmd, mceinxPtr,mycmdInfo.cmdBufSize,
                          mycmdInfo.cmdBuf,status);
  if ( *status != STATUS__OK)
  {
    sprintf(errmsg,
         "mce_set: mcexml_translate returned with bad status");  
    ErsRep (0, status, errmsg);
    return;
  }
  if( (*status=sdsu_command_mce (con,mycmdInfo.cmdBuf,&mycmdInfo.reply)) < 0)
  { 
     strcpy(mycmdInfo.mceCmd,cmd);
     mce_error_reply(myInfo,&mycmdInfo,status); 
     sprintf(errmsg,"mce_set: sdsu_command_mce returned bad status");
     ErsRep (0, status, errmsg);
     return;
  }
}

/**
 * \fn void mce_read_frame_rate(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  struct mcexml_struct *mceinxPtr, long *rate, 
 *  StatusType *status)
 *
 * \brief function: 
 *  read back MCE parameters Val and calculate frame rate 
 *
 * \param con       SDSU context structure pointer
 * \param myInfo   dasInfoStruct_t poiter
 * \param mceinxPtr  mcexml_struct pointer
 * \param rate       long pointer for frame rate
 * \param status    StatusType pointer.  given and return
 *
 * Used to be sc2dalib_readframeRate
 */
/*+ mce_read_frame_rate
*/
void mce_read_frame_rate
(
SDSU_CONTEXT          *con,
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceinxPtr,
long                  *rate,
StatusType            *status
)
{
  if (*status != STATUS__OK) return;

  int dataRate;
  // MCE firmware only allows rb rc1, rc2 rcs,  not rb sys
  char  rowlengthCmd[]="rb sys row_len 1 ";
  char    numrowCmd[] ="rb sys num_rows 1 ";
  char  datarateCmd[] ="rb cc  data_rate 1";
  char  fbJumpCmd[] ="rb rcs en_fb_jump 1";

  mce_read(con,myInfo,mceinxPtr,rowlengthCmd,&myInfo->rowLength,1,status);
  if (!StatusOkP(status)) 
    return;

  mce_read(con,myInfo,mceinxPtr,numrowCmd,&myInfo->numRows,1,status);
  if (!StatusOkP(status)) 
    return;

  mce_read(con,myInfo,mceinxPtr,datarateCmd,&dataRate,1,status);
  if (!StatusOkP(status)) 
    return;

  mce_read(con,myInfo,mceinxPtr,fbJumpCmd,&myInfo->en_fb_jump,1,status);
  if (!StatusOkP(status)) 
    return;

  *rate= (float)50000000/(myInfo->rowLength*myInfo->numRows*dataRate)+0.5;
  jitDebug(2,"=============================\n");
  jitDebug(2,"MCE's rowLength(=%d) numRows(=%d),dataRate(=%d)\n",
            myInfo->rowLength, myInfo->numRows, dataRate);
  jitDebug(2,"The calculated frameRate =%ld Hz\n",*rate); 
  jitDebug(2,"============================\n");
}


/**
 * \fn void mce_checksum(dasInfoStruct_t *myInfo,
 *      dasCmdInfo_t *mycmdInfo,   StatusType *status)
 *
 * \brief function
 *  check reply checksum and report it
 *
 * \param   myInfo    dasInfoStruct_t structure pointer
 * \param  mycmdInfo  dasCmdInfo_t structure pointer
 * \param  status     StatusType.  given and return
 *
 * Used to be: sc2dalib_chkChecksum
 */
/*+ mce_checksum 
*/
void mce_checksum
(
dasInfoStruct_t  *myInfo,   
dasCmdInfo_t     *mycmdInfo,
StatusType       *status
)
{
  int         chksum=0,chksumRec;
  int         chksumNo,j;
  static char chksumerr[FILE_LEN];
 
  if (!StatusOkP(status)) return;

  // my reply struct has total 64 words
  //{ uint32  status  ;    //<cmd and status, i.e "WMOK" or "WMER" or errNo for NEW PROTOCOL 
  //  uint32  data[62];    //< Data argument 
  //  uint32  checksum;    //< word64 checksum 
  //}
  //
  // MCE's reply is (old protocol) 
  // { status; 
  //   cardId+parId;
  //   58 data (max);
  //   chksum
  // }
  // mycmdInfo->replySize <=61
  // ( 1(status) + 1(cardId+parId) + 58(No_data) + 1(chksum) )
  //
  // MCE's reply is (new protocol) 
  // {  errno; 
  //   cardId+parId;
  //   58 data (max);
  //   chksum
  // }
  // mycmdInfo->replySize <=61
  // ( 1(errno) + 1(cardId+parId) +58(No_data) + 1(chksum) )
  //
  // hence, it never fills the reply structure upto reply.checksum

  chksumNo=mycmdInfo->replySize - 2;
  chksum ^= mycmdInfo->reply.status;

  for (j=0;j<chksumNo;j++)     
    chksum ^=mycmdInfo->reply.data[j];

  // checksum from MCE is alway in reply.data[] 
  chksumRec=mycmdInfo->reply.data[chksumNo];
  jitDebug(2, "mce_checksum: RTL's CHECKSUM <%08X>, MCE's CHECKSUM  <%08X>\n",
           chksum,chksumRec);

  if( chksumRec !=chksum)
  {
    sprintf(chksumerr, "the calculated chksum(%#X)!=mcechksum(%#X)",
            chksum, chksumRec);
    utils_msg(myInfo, " mce_checksum: %s", chksumerr,
              USE_ERSOUT, status);
    //only place myInfo is used.
    *status=DITS__APP_ERROR;
  }
}


/**
 * \fn void mce_results(dasInfoStruct_t *myInfo, 
 *  dasCmdInfo_t *mycmdInfo,  StatusType *status)
 *
 * \brief functaion
 *  display reply from MCE 
 *
 * \param myInfo   dasInfoStruct_t structure pointer
 * \param mycmdInfo dasCmdInfo_t structure pointer 
 * \param status    StatusType.  given and return
 *
 * Used to be sc2dalib_dispResults
 */
/*+ mce_results - 
*/
void mce_results
(
dasInfoStruct_t       *myInfo,   
dasCmdInfo_t          *mycmdInfo, 
StatusType            *status
)
{
  int   dataNo,i,j;
  int  multi,remain;
  char  remainData[FILE_LEN], remainData1[FILE_LEN]="";
  char  localmsg[FILE_LEN];
  char  fwrevString[]="fw_rev";
  char  heatString[]="bc1 bias";
  char  biasString[]="bc2 bias";
  char  arrayString[]="cc array_id";


  if (!StatusOkP(status)) return;

  sprintf(localmsg,"CardIdParam=%08lX Word3=%08lX",
             mycmdInfo->reply.data[0],mycmdInfo->reply.data[1]);
  if (mycmdInfo->reply.status == MCE_WBOK )
    jitDebug(2,"WBOK %s\n",localmsg);  
  else if (mycmdInfo->reply.status == MCE_GOOK )
    jitDebug(2,"GOOK %s\n",localmsg);
  else if (mycmdInfo->reply.status == MCE_STOK )
    jitDebug(2,"STOK %s\n",localmsg);
  else if (mycmdInfo->reply.status == MCE_RSOK )
    jitDebug(2,"RSOK %s\n",localmsg);

  // all error msgs have been displayed in mce_error_reply
  // only use reply.data[]:1,.... 
  // not include reply.status,data[0]=(cardId+ParaId) and chechsum
  dataNo=mycmdInfo->replySize-3;  
  jitDebug(2,"dataNO (%d) from RBOK)\n",dataNo);
 
  // firmware revision
  if( (strstr(mycmdInfo->mceCmd,fwrevString)!= NULL) && 
             (mycmdInfo->reply.status == MCE_RBOK) 
    )  
  {
    int buildNo, minorNo, majorNo;
    //print revision 
    buildNo=(int)( (mycmdInfo->reply.data[1]    )&0x0000FFFF );
    minorNo=(int)( (mycmdInfo->reply.data[1]>>16)&0x000000FF );
    majorNo=(int)( (mycmdInfo->reply.data[1]>>24)&0x000000FF );
    MsgOut(status,"%s => Revision %d .%d build %d", 
     mycmdInfo->mceCmd,majorNo,minorNo,buildNo);
  }
  else if( ( strstr(mycmdInfo->mceCmd,heatString)!= NULL ||
             strstr(mycmdInfo->mceCmd,arrayString)!= NULL ||
       strstr(mycmdInfo->mceCmd,biasString)!= NULL ) && 
                   (mycmdInfo->reply.status == MCE_RBOK) 
         )  
  {
    mce_reply(mycmdInfo,1,status);
  }
  else if( ( myInfo->filesvFlag > 1) && ( myInfo->actionFlag==MCECMDACTION) && 
           (mycmdInfo->reply.status == MCE_RBOK) 
         )  
  {
    mce_reply(mycmdInfo,dataNo,status);
  }
   
  if(myInfo->debuglvl==2)     
  {
    if(mycmdInfo->reply.status == MCE_RBOK )
    {
      multi=dataNo/8; remain=dataNo%8;
      jitDebug(2,"Data from ReadBlock\n");
      for (j=0;j<multi;j++)  
      {
        jitDebug(2,"  %8ld %8ld %8ld %8ld %8ld %8ld %8ld %8ld\n",
         mycmdInfo->reply.data[j*8+1],mycmdInfo->reply.data[j*8+2],
         mycmdInfo->reply.data[j*8+3],mycmdInfo->reply.data[j*8+4],
         mycmdInfo->reply.data[j*8+5],mycmdInfo->reply.data[j*8+6],
         mycmdInfo->reply.data[j*8+7],mycmdInfo->reply.data[j*8+8]);
      }
      if(remain)
      {
        for (i=1;i<=remain;i++)  
        {
          sprintf(remainData,"%s %8ld", remainData1,
                   mycmdInfo->reply.data[j*8+i]);
          strcpy(remainData1,remainData);
        } 
        jitDebug(2," %s\n", remainData);
      }
    }
  }
}


/**
 * \fn void mce_reply(dasCmdInfo_t *mycmdInfo, int howMany,
 *  StatusType *status)
 *
 * \brief functaion
 *  display reply 
 *
 *\param mycmdInfo  dasCmdInfo_t structure pointer 
 *\param howMany    int
 *\param status     StatusType.  given and return
 *
 * Used to be sc2dalib_dispReply
 */
/*+ mce_reply - 
*/
void mce_reply
(
dasCmdInfo_t  *mycmdInfo,
int           howMany, 
StatusType    *status
)
{
  int   i; 
  char  msg[FILE_LEN],tmp[FILE_LEN]="";
 
  if (!StatusOkP(status)) return;
  
  for ( i=0; i<howMany; i++)
  {
    sprintf(msg,"%d ",(int)mycmdInfo->reply.data[i+1]);
    strcat(tmp,msg);
  }
  MsgOut(status," %s",tmp);   
}


/**
 * \fn void mce_error_reply(dasInfoStruct_t *myInfo,
 *  dasCmdInfo_t *mycmdInfo, StatusType *status)
 *
 * \brief functaion
 *  check reply status from MCE and report it
 *
 * \param myInfo    dasInfoStruct_t  pointer
 * \param mycmdInfo  dasCmdInfo_t  pointer
 * \param status     StatusType.  given and return(STATUS__OK)
 *
 * report the Hardware error
 * Used to be sc2dalib_mceerrRep
 */
/*+ mce_error_reply - 
*/
void mce_error_reply
(
dasInfoStruct_t *myInfo,
dasCmdInfo_t    *mycmdInfo,
StatusType      *status
)
{
  time_t   tm;
  char localmsg[FILE_LEN];

  // message print and save to logfile, 
  // flag=0, use MsgOut, 
  // flag=1, use ErsOut, 
  // flag=2 use ErsRep
  tm = time(NULL);
  fprintf (myInfo->fpLog,"%s",asctime(localtime(&tm)) );

  // *status != OK from entry, *status not changed after call
  // utils_msg as 2 is used
  utils_msg(myInfo,"mce_error_reply: executing command (%s)",
            mycmdInfo->mceCmd, USE_ERSREP, status);
  if (mycmdInfo->reply.status==MCE_NONREPLY)   
  {
    utils_msg(myInfo,
      "mce_error_reply: NO REPLY from driver when sdsu_command_mce ",
      "", USE_ERSOUT, status);
    utils_msg(myInfo, "mce_error_reply: mceTask exits?",
              "", USE_ERSREP, status);
  }
  else if (mycmdInfo->reply.status==DRV_MCEERR)   
  {
    if( mycmdInfo->reply.data[0]==PCI2MCE_TMO || 
        mycmdInfo->reply.data[0]==INT2MCE_TMO )
    { 
      utils_msg(myInfo, "mce_error_reply: %s",
                drvErrors[mycmdInfo->reply.data[0]], USE_ERSREP, status);
    }
    else 
    {
      utils_msg(myInfo, "mce_error_reply: %s",
                drvErrors[mycmdInfo->reply.data[0]], USE_ERSREP, status);
      sprintf(localmsg, "%#lX", mycmdInfo->reply.data[1]);
      utils_msg(myInfo,"mce_error_reply: %s",localmsg,USE_ERSREP,status);
    }
  }
  else                                
  {
    sprintf(localmsg,"CardIdParam=%08lX ErrNo=%08lX",
             mycmdInfo->reply.data[0],mycmdInfo->reply.data[1]);
    if (mycmdInfo->reply.status == MCE_WBER ) 
      utils_msg(myInfo,
                       "mce_error_reply: WBER %s",localmsg,USE_ERSREP,status); 
    else if (mycmdInfo->reply.status == MCE_RBER )
      utils_msg(myInfo,
                       "mce_error_reply: RBER %s",localmsg,USE_ERSREP,status); 
    else if (mycmdInfo->reply.status == MCE_GOER )
     utils_msg(myInfo,
                       "mce_error_reply: GOER %s",localmsg,USE_ERSREP,status); 
    else if (mycmdInfo->reply.status == MCE_STER )
      utils_msg(myInfo,
                       "mce_error_reply: STER %s",localmsg,USE_ERSREP,status); 
    else if (mycmdInfo->reply.status == MCE_RSER )
      utils_msg(myInfo,
                       "mce_error_reply: RSER %s",localmsg,USE_ERSREP,status);
    else if (mycmdInfo->reply.status == MCE_RSER )
      utils_msg(myInfo,
                       "mce_error_reply: RSER %s",localmsg,USE_ERSREP,status);
    else 
    {
      sprintf(localmsg,"the error Status word(%08lx) is wrong",
                mycmdInfo->reply.status);
      utils_msg(myInfo,"mce_error_reply:%s",localmsg, USE_ERSREP,status);  
    }
  } 
}

/**
 * \fn void mce_cmd(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo,
 *  dasCmdInfo_t *mycmdInfo, char *dateTime,StatusType *status)
 *
 * \brief function
 *  send a cmd to MCE
 *
 * \param con        SDSU context structure pointer 
 * \param myInfo    dasInfoStruct_t pointer
 * \param mycmdInfo  dasCmdInfo_t pointer
 * \param dateTime   dateTime string pointer         
 * \param status     StatusType.  given and return
 *
 * send the command buffer to MCE
 * save the command buffer and reply or error in logfile 
 *
 * Used to be sc2dalib_Cmd
 */
/*+ mce_cmd 
 */
void mce_cmd
(
SDSU_CONTEXT      *con,
dasInfoStruct_t   *myInfo,
dasCmdInfo_t      *mycmdInfo,
char              *dateTime,
StatusType        *status
)
{
  char   *byte;
  int   *longword;
  
  if (!StatusOkP(status)) return;
  
  // initial reply as non-reply
  mycmdInfo->reply.status = MCE_NONREPLY;
  longword = (int *)mycmdInfo->cmdBuf;
  byte = (char*)&mycmdInfo->reply;
  myInfo->glbCount++;

  /* In this case we have been building an AST keymap, but have reached the end 
     of the file so we just want to write the key map out and end */
  if ((myInfo->filesvFlag == 5) && (myInfo->astMapState == 2))
    {
      mce_reply_hash(myInfo,myInfo->fpMcecmd,mycmdInfo->replySize,
        byte,mycmdInfo->mceCmd);
      return;
    }
#ifdef PRINT_LOG_MESG
  if (myInfo->actionFlag == MCESTATUSACTION || myInfo->actionFlag == SERVOACTION)
  {
    fprintf(myInfo->fpLog,"\nCMD(%ld): <%s> sent to the MCE\n",
          myInfo->glbCount,mycmdInfo->mceCmd);
  }
  else
  {
    fprintf(myInfo->fpLog,"\n<%s>CMD(%ld): <%s> sent to the MCE\n",
          dateTime,myInfo->glbCount,mycmdInfo->mceCmd);
  }
#endif

  if (myInfo->filesvFlag == 1 )
  {
    if (myInfo->actionFlag == MCESTATUSACTION)
    {
      fprintf(myInfo->fpMcecmd, "<%s> ",mycmdInfo->mceCmd);
    }
    else
    {
      // also save into the specified file, 
      fprintf(myInfo->fpMcecmd,"\n<%s>CMD(%d): <%s> sent to MCE \n",
            dateTime,myInfo->trkNo,mycmdInfo->mceCmd);
    }
  }

  // return SDSU_OK if WBOK RBOK GOOK, STOK or RSOK
  // or errNo for new protocol
  *status=sdsu_command_mce (con,mycmdInfo->cmdBuf,&mycmdInfo->reply);

  // save cmd and reply first, then check *status, so that we have checksum saved 
  // don't save cmd buffer to the logfile if it is a batch or
  // servo cmd or track heat action
  if( con->process.framesetup != DA_MCE_SINGLEDATA)
  { 
#ifdef PRINT_LOG_MESG
    if ( myInfo->actionFlag !=MCESTATUSACTION && 
        myInfo->actionFlag !=SERVOACTION )
      mce_cmd_save(myInfo->fpLog,mycmdInfo->cmdBufSize,mycmdInfo->cmdBuf);
#endif
    if (myInfo->filesvFlag==1 && myInfo->actionFlag !=MCESTATUSACTION)
      mce_cmd_save(myInfo->fpMcecmd,mycmdInfo->cmdBufSize,
                      mycmdInfo->cmdBuf);
  }
  mycmdInfo->replySize=con->packsize;
#ifdef PRINT_LOG_MESG
  mce_reply_save(myInfo->fpLog,mycmdInfo->replySize, 
                                    byte,myInfo->glbCount);
#endif
  if (myInfo->filesvFlag == 1)
  {
    if ( myInfo->actionFlag != MCESTATUSACTION )
      mce_reply_save(myInfo->fpMcecmd,mycmdInfo->replySize,byte,
                                      myInfo->trkNo);
    else
      sc2dalib_savemcerepData(myInfo,myInfo->fpMcecmd,mycmdInfo->replySize,byte);
  }
  else  if (myInfo->filesvFlag == 5)
    {
      mce_reply_hash(myInfo,myInfo->fpMcecmd,mycmdInfo->replySize,
        byte,mycmdInfo->mceCmd);
    }
  
  if( *status !=SDSU_OK)
  { 
    *status=DITS__APP_ERROR;
    if (con->flag.protoCOL==OLDPROTO)
      mce_error_reply(myInfo,mycmdInfo,status); 
    else
      // need to decide if it is fatal or non fatal error
      sc2dalib_mceerrnewRep(myInfo,mycmdInfo,status); 
    return;
  }
}

/**
 * \fn void mce_cmd_save(FILE *fpout,int mxbuf,char buffer[mxbuf])
 *
 * \brief functaion
 *  save command buffer into a file
 *
 * \param fpout    FILE pointer for store the cmdBuf sent out 
 * \param mxbuf    int the buffer size
 * \param buffer[] char array for the command buffer
 *
 * Used to be sc2dalib_savecmdBuf
 */
/*+ mce_cmd_save - save command buffer into a file
*/
void mce_cmd_save
(
FILE *fpout,          
int   mxbuf,
char  buffer[mxbuf]  
)
{
   int *longword;
   int  j,l,toploop,midloop,remain;

   midloop = 8;
   mxbuf /= 4;
   toploop = (mxbuf) / midloop;
   remain = (mxbuf) % midloop;
   longword = (int *)buffer;

   for (j=0;j<toploop;j++)
   {
      for(l=0;l<midloop;l++)
         fprintf(fpout,"%08X ",*(longword+ l+midloop*j) );
      fprintf(fpout,"\n");
   }
   toploop =toploop*midloop;
   for(l=0;l<remain;l++)
      fprintf(fpout,"%08X ", *(longword+ l+ toploop) );
   if(remain)
      fprintf(fpout,"\n");
}

/**
 * \fn void mce_reply_save(FILE *fpout,int mxbuf,char * byte,
 *   uint32 count)
 *
 * \brief function
 *  save MCE reply into a file
 *
 * \param fpout  FILE pointer for store the reply 
 * \param mxbuf  int: the buffer size in 32bits
 * \param byte   pointer to the reply buffer
 * \param count  the count_th repy since dasDrama starts
 *
 * Used to be sc2dalib_savemceReply
 */
/*+ mce_reply_save -save MCE reply into a file
*/
void mce_reply_save
(
FILE *fpout,    //store the cmdBuf sent to PCI, in byte order  
int   mxbuf,
char  *byte,    // pointer to reply Buffer 
uint32 count    // the count_Th reply 
)
{
  int   *longword;
  int   j,l,toploop,midloop, remain;

  midloop=8;
  toploop=(mxbuf)/midloop;
  remain =(mxbuf)%midloop;
  longword = (int *)byte;
 
  fprintf(fpout,"CMD(%ld)'s reply is here \n",count);
  for (j=0;j<toploop;j++)
  {
    for(l=0;l<midloop;l++)
      fprintf(fpout,"%08X ",*(longword+ l+midloop*j) );
    fprintf(fpout,"\n");
  }
  toploop =toploop*midloop;
  for(l=0;l<remain;l++)
    fprintf(fpout,"%08X ", *(longword+ l+ toploop) );
  if(remain)
    fprintf(fpout,"\n");
}


/**
 * \fn void mce_reply_hash(dasInfoStruct_t *myInfo, FILE *fpout,
 *  int mxbuf, char * byte, char *mycmd)
 *
 * \brief function
 *  save the MCE reply data into a file and put it into an AST hash
 *
 * \param myInfo     dasInfoStruct_t pointer
 * \param fpout      FILE pointer of file where to store the reply 
 * \param mxbuf      int: the buffer size in 32bits
 * \param byte       pointer to the reply buffer
 * \param mycmd      Pointer to the actual command
 *
 * Used to be sc2dalib_saveMceRepToHash
 */
/*+ mce_reply_hash
*/
void mce_reply_hash
(
dasInfoStruct_t  *myInfo,
FILE             *fpout,      // file pointer where to store the replied data  
int              mxbuf,       // How many words came back from the MCE
char             *byte,       // Pointer to reply Buffer
char            *mycmd       // Pointer to command which caused the data
)
{
  int   *longword;
  int   j, expected;
  char *myPntr;
  char delim[6];
  char myBuff[80];
  char card[8];
  char cmd[30];
  static AstKeyMap *theMap;
  char key[40];
  int values[50];
  int status;
  int lenCmd, lenPossible;

  /* myInfo->astMapState will be 0 the first time this routine is called,
     1 on all subsequent calls except when we are done then it will be 2 */
  if(myInfo->astMapState == 0)
  {
    theMap = NULL;
    theMap = astKeyMap(" ");
  }

  else if(myInfo->astMapState == 2)
  {
    /* Write map out to a file called astKeyMapFile */
    status = SAI__OK;
    sprintf(myBuff,"%s/astKeyMapFile", getenv ("SC2SCRATCH"));
    atlShow((AstObject *)theMap, myBuff, "", &status );
    astAnnul(theMap);
    return;
  }

  /* Until the very end we want the state equal to 1 */
  myInfo->astMapState = 1;

  /* The command will be something like this:rb bc2 bias 1
     Where the second thing (bc2) is the card number and the
     third thing is the command and the fourth thing is how
     many words I was expecting to get back, so parse it */
  delim[0] = ' ';
  delim[1] = 0;
  myBuff[18]=0;
  strncpy(myBuff, mycmd, 40);
  myPntr = strtok(myBuff, delim); /* This should be pointing at the rb */
  myPntr = strtok(NULL, delim);   /* The card number */
  strncpy(card,myPntr,8);
  myPntr = strtok(NULL, delim);   /* The command itself */
  strncpy(cmd,myPntr,20);
  myPntr = strtok(NULL, delim);   /* The number I was expecting */
  expected = atoi(myPntr);

  /* Chop the end off the command so that it plus the card string, plus one <= 15 characters */
  lenPossible = 14 - strlen(card);
  lenCmd = strlen(cmd);
  if(lenCmd > lenPossible) cmd[lenPossible-1] = 0;
  
  fprintf(fpout,"%s_%s:  ",card, cmd);

  if(expected > (mxbuf-3)) expected = mxbuf-3;

  // we do not want the first two RBOK and CardId_paraId
  // and last chksum

  longword = (int *)byte;

  for (j=0; j<expected; j++)
  {
    fprintf(fpout,"%8d ",*(longword+2 + j ) );
    values[j] = *(longword+2 + j );
  }
  fprintf(fpout,"\n");

  /* build up key word and add to key map */
  strcpy(key,cmd);
  strcat(key,"_");
  strcat(key,card);
  astMapPut1I(theMap, key, expected, values, " ");
}

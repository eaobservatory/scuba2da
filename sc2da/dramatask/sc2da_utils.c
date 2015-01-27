#include "sc2da_utils.h"

static long qlFileCounter;  

/**
 * \fn void utils_msg(dasInfoStruct_t *myInfo, char *string, 
 *   char *string2, int flag, StatusType *status)
 *
 * \brief function
 *  report error and save it into logfile
 *
 * \param myInfo  dasInfoStruct_t  pointer
 * \param string   char pointer  the string
 * \param string2  char pointer  the second string
 * \param flag     int  =0, use MsgOut, =1 use ErsOut, =2 ErsRep  3, printf, >3 no print
 * \param status   StatusType.  given and return(STATUS__OK)
 *
 * Used to be sc2dalib_msgprintSave
 */
/*+ utils_msg
*/
void utils_msg
(
dasInfoStruct_t *myInfo, 
char            *string, 
char            *string2,
int             flag,
StatusType      *status
)
{
  if(flag == USE_MSGOUT)
  	MsgOut(status, string,string2);
  else if (flag == USE_ERSOUT)
    ErsOut(0,status,string,string2);
  else if  (flag == USE_ERSREP)
    ErsRep(0,status,string,string2);
  else if  (flag == USE_PRINT)
  {
    printf(string,string2); printf("\n");
  }

  fprintf(myInfo->fpLog, string,string2);
  fprintf(myInfo->fpLog, "\n");
  fflush(myInfo->fpLog);
  if (myInfo->filesvFlag==1)
  {
    // also save into the specified file
    if(myInfo->fpMcecmd !=NULL)
    {  
      fprintf(myInfo->fpMcecmd,string,string2);
      fprintf(myInfo->fpMcecmd, "\n");
      fflush(myInfo->fpMcecmd);
    }
  }
}


/**
 * \fn utils_fclose(FILE **fp)
 * 
 * \brief function
 *  Checks if file is closed, should be used before fopen
 *
 * \param **fp file pointer
 * Stat: used 25 times
 * Used to be my_fclose
 */
/*+ utils_fclose
*/
void utils_fclose(FILE **fp)
{
  if(fp == NULL || *fp == NULL)
    {
      return;
    }
  fclose(*fp);
  *fp = NULL;
}


/**
 * \fn utils_close_files(dasInfoStruct_t *myInfo)
 * 
 * \brief function
 *  Closes myInfo-> (fpData, fpMcecmd, fpBatch, Strchart, fpOtheruse)
 *  Flushes log file.
 *
 * \param myInfo
 *  Stat: used 2 times
 * Used to be my_closeFiles
 */
/*+ utils_close_files
*/
void utils_close_files(dasInfoStruct_t *myInfo)
{
  fflush(myInfo->fpLog);
  utils_fclose(&(myInfo->fpData));
  jitDebug(2,"utils_close_files: closed fpData\n"); 

  utils_fclose(&(myInfo->fpMcecmd));
  jitDebug(2,"utils_close_files: closed fpMcecmd\n"); 

  utils_fclose(&(myInfo->fpBatch));
  jitDebug(2,"utils_close_files: closed fpBatch\n"); 

  utils_fclose(&(myInfo->fpStrchart));
  jitDebug(2,"utils_close_files: closed fpStrchart\n"); 

  utils_fclose(&(myInfo->fpOtheruse));
  jitDebug(2,"utils_close_files: closed fpOtheruse\n"); 
}


/**
 * \fn void utils_open_files(dasInfoStruct_t *myInfo,
 * int getData, int getBatch, StatusType *status)
 *
 * \brief function
 *  open all files, set myInfo->trkNo=0 if no cmdrepFile exits
 *
 * \param myInfo     dasInfo structure pointer
 * \param getData    int
 * \param getBatch   int
 * \param status     StatusType.  given and return
 *
 * getData=
 *   DATFILENOAPPEND,   // data file no appending
 *   DATFILEAPPEND,     //appending for data file
 *   NODATAFILE,
 * getBatch=
 *   NOBATCHFILE,
 *   BATCHFILE,          // it is batch
 *   OTHERFILE           // it is other file,i.e setup file
 *
 * Used to be sc2dalib_openFiles
 * Subsumed sc2da_isFile
 */
/*+ utils_open_files
*/
void utils_open_files
(
dasInfoStruct_t *myInfo,       
int             getData,
int             getBatch,
StatusType      *status
)
{
  struct stat sbuf;
  if (!StatusOkP(status)) return;

  utils_close_files(myInfo);

  //filesvFlag !=1, don't save cmd-reply to fpMcecmd

  jitDebug(8,
      "utils_open_files: cmdBuf to MCE and Reply from MCE stored in %s\n", 
      myInfo->logfileName);

  if ((myInfo->filesvFlag==1) || (myInfo->filesvFlag==5)) 
  {
    // if no file exists,
    // reset the trkNo, otherwise, leave trkNo unchanged
    if (lstat(myInfo->cmdrepFile, &sbuf) !=0) 
  	{
      myInfo->trkNo=0;
    }
    else
      myInfo->trkNo++;
       
    if ((myInfo->fpMcecmd=fopen64(myInfo->cmdrepFile,"w")) == NULL)
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "utils_open_files:cmdrepFile failed to open %s",
              myInfo->cmdrepFile);
      return;
    }
    jitDebug(8,
       "utils_open_files: the CMD-REPLY results are also stored in %s\n",
       myInfo->cmdrepFile);
  }
  else if (myInfo->actionFlag==SERVOACTION)
  {
    if((myInfo->fpMcecmd = fopen64(myInfo->cmdrepFile,"w")) == NULL)
    {
      *status = DITS__APP_ERROR;
      ErsRep (0, status, "utils_open_files:cmdrepFile 1  failed to open %s",
              myInfo->cmdrepFile);
      return;
    }
    jitDebug(16,
       "utils_open_files: the lockpoint results are also stored in %s\n",
       myInfo->cmdrepFile);
  }

  if( myInfo->filesvFlag !=0 ) 
  {
    if ( getData !=NODATAFILE )
    { 
      if ( getData==DATFILEAPPEND ) 
      {
        // open another file for appending, mostly used for taking 
        // frame data
	if((myInfo->fpData = fopen64(myInfo->dataFile,"a")) == NULL )
	  {
	    *status = DITS__APP_ERROR;
	    ErsRep (0, status, "utils_open_files: Error- failed to open file %s", myInfo->dataFile); 
	    return;
	  }
      }
      else  //( getData==DATFILENOAPPEND ) 
      {
        // open another file, but no appending
	if((myInfo->fpData = fopen64(myInfo->dataFile,"w")) == NULL )
	  {
	    *status = DITS__APP_ERROR;
	    ErsRep (0, status, "utils_open_files: Error- failed to open file %s", myInfo->dataFile); 
	    return;
	  }
      }
      if ( myInfo->fpData==NULL)
      {
        *status = DITS__APP_ERROR;
        ErsRep(0,status,"utils_open_files: dataFile failed to open %s",
               myInfo->dataFile);
        ErsRep(0,status, "utils_open_files: %s",strerror(errno));
        return;
      }
      jitDebug(8,"utils_open_files: the data result stored in %s\n",
              myInfo->dataFile);
    }
  }
  // open another file for stripchart, do not append so file will not grow 
  if(   ( (myInfo->actionFlag==SERVOACTION) && (myInfo->filesvFlag !=0 ) ) ||
        ( (myInfo->actionFlag==TRKHEATACTION) && ( (myInfo->heatSlp.option & 0x01)!=0 ) ) ||
        (myInfo->actionFlag==TRKSQ2FBACTION) 
    )
  {
    if ((myInfo->fpStrchart=fopen64(myInfo->strchartFile,"w"))==NULL)
    {
      *status = DITS__APP_ERROR;
      ErsRep(0,status,"utils_open_files: strchartFile failed to open %s",
             myInfo->strchartFile);
      ErsRep(0,status, "utils_open_files: %s",strerror(errno));
      return;
    }
  }
  
  // open another file, for batch or setup
  if(getBatch==BATCHFILE || getBatch==OTHERFILE)  
  {
   jitDebug(8,"utils_open_files: the specified file:%s\n",
              myInfo->batchFile);
    if ((myInfo->fpBatch = fopen(myInfo->batchFile,"r")) == NULL)
    {
      jitDebug(8,"utils_open_files: batchFile failed to open: %s\n",
              myInfo->batchFile);     
      *status=DITS__APP_ERROR;
      ErsRep(0,status,"utils_open_files: batchFile failed to open: %s\n",
             myInfo->batchFile);
      return;
    }
    jitDebug(8,"utils_open_files: the batch or specified file:%s\n",
             myInfo->batchFile);     
  }
}

/**
 * \fn void utils_update_debug(SDSU_CONTEXT *con, dasInfoStruct_t  *myInfo, 
    StatusType *status)
 *
 * \brief function
 *  get DEBUG parameter and update debug level 
 *
 * \param con     SDSU_CONTEXT pointer 
 * \param myInfo dasInfoStruct_t pointer
 * \param status  StatusType pointer.  given and return
 *
 * debuglvl:bit  7-0   for DAS
 * debuglvl:bit  15-8   for DATAHANDLE
 * debuglvl:bit  23-16   for INTERFACE
 * debuglvl:bit 31-24  for DRIVER
 *
 * debug level: 
 *  NO_GEN_MSG=0,  diplay very limited message, no filenames
 *  DISP_GEN_MSG=1,  display general message, no filenames
 *  DISP_FILE_MSG=3, display only filenames
 *  the followings are not used currently
 *  DATA_BEG_MSG=2,  display part of the beg. of data ,no filenames
 *  DISP_GEN2_MSG=4,  disply special message
 *  DISP_MCETASK=5,  disply general message in mce task
 *  DISP_PCITASK=6, disply general message in pci task
 *  DISP_DATATASK=7, disply general message in data task
 *
 * Used to be sc2dalib_updateDebug
 */
/*+ utils_update_debug 
*/
void utils_update_debug
(
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
StatusType      *status    
)
{
  long debuglvl;
  
  if (!StatusOkP(status)) return;

  SdpGeti("DEBUG", &debuglvl, status);
  con->process.debuglvl[DAS]       =(char)((debuglvl    )& 0x000000FF);
  con->process.debuglvl[DATAHANDLE]=(char)((debuglvl>> 8)& 0x000000FF);
  con->process.debuglvl[INTERFACE] =(char)((debuglvl>> 16)& 0x000000FF);
  con->process.debuglvl[DRIVER]    =(char)((debuglvl>> 24)& 0x000000FF);
  myInfo->debuglvl=(long)con->process.debuglvl[DAS];
}


/**
 * \fn void utils_alloc_memory(SDSU_CONTEXT *con, 
 * dasInfoStruct_t *myInfo,  int startSeq, int endSeq, StatusType *status)
 *
 * \brief function
 *  allocate share memory for information between drama 
 *  task and data handle task
 *
 * \param con      SDSU_CONTEXT structure pointer 
 * \param myInfo   dasInfoStruct_t pointer
 * \param startSeq  start Sequence Num
 * \param endSeq    end Sequence Num
 * \param status    StatusType.  given and return
 *
 *  // the whole shared memory looks like 
 * //  
 * // |  int        ithqlrecord |
 * // |-------------------------|
 * // |  qlNum*       QL struct |
 * // |               ......... |
 * // |               QL struct |
 * // |-------------------------|
 * // |  qlNum* int   lkupflage |
 * // |               ........  |   
 * // |               lkupflag  |
 * // |-------------------------|
 * // |  myInfo->numFrame *  sc2head struct |
 * // |                     ............... |
 * // |                     sc2head struct  |
 * // |-------------------------------------|
 * // in SCAN mode, the QL_STRUCT is different
 *
 * Used to be sc2dalib_allctsharedMem
 */
/*+ utils_alloc_memory   
 */
void utils_alloc_memory
(
SDSU_CONTEXT       *con,         
dasInfoStruct_t    *myInfo,
int                startSeq,
int                endSeq,
StatusType         *status
)
{
  int     i, *lkupPtr;
  int     sharememSize=0;
  int     *scannumPtr;
  int     obsMode;
  double  *timePtr;  
  char    *scanfileName;
  QL_SCAN *scanPtr; 
  int     X,Y;

  if (!StatusOkP(status)) return;

  myInfo->numFrame=(endSeq-startSeq+1);
  obsMode=myInfo->parshmPtr->obsMode;

  myInfo->qlNum=myInfo->numFrame/myInfo->procNum;
  // in case we have no-integer number of procNum
  // for DREAM, the qlNum may not be right,
  //(it is maximum, as default SMU_PATTERN=64)
  // need to re-calculated after reading weightsfile in DH-task
  if ( (myInfo->numFrame - myInfo->qlNum*myInfo->procNum) >0 )
    myInfo->qlNum++;
  
  myInfo->parshmPtr->qlNum=myInfo->qlNum;

  jitDebug(2,"utils_alloc_memory: numfram(%d),procNum(%d), qlNum(%d)\n",
            myInfo->numFrame, myInfo->procNum, myInfo->qlNum);

  sharememSize = sizeof(int); // which QL is updated

  if ( obsMode==OBS_STARE )
   sharememSize += sizeof(QL_STRUCT)*myInfo->qlNum; 
  else if ( obsMode==OBS_DREAM )
   sharememSize += sizeof(RECONSTR_STRUCT)*myInfo->qlNum; 
  else if ( myInfo->parshmPtr->obsMode ==OBS_SCAN )
    sharememSize += sizeof(QL_SCAN)*myInfo->qlNum;   

  sharememSize += sizeof(JCMTState)*myInfo->numFrame;  // for all sc2head
  sharememSize += sizeof(int)*myInfo->qlNum;      // for all lookupflag table 
  
  jitDebug(2, "sc2da_seq: sizeof(JCMTState)= %d bytes, required sharedMem = %d MBs",
              sizeof(JCMTState),sharememSize/(1024*1024));
 
  utils_set_memory(con, myInfo, SHAREDM_OBS, sharememSize, status);
  if (!StatusOkP(status)) 
     return;  

  myInfo->qlPtr=(int*)myInfo->sharedShm;

  if ( obsMode==OBS_STARE ) 
  {
    if(SC2STORE__COL_INDEX == 0)
      {
	X = COL_NUM;
	Y = ROW_NUM-1;
      }
    else
      {
	X = ROW_NUM-1;
	Y = COL_NUM;
      }

    int s1=sizeof(QL_STRUCT);
    int s2=sizeof(QL_IMAGE); 

    jitDebug(2,"  sizeof(double)(%d) doube[Y*X] (%d) \n",
            sizeof(double), sizeof(double)*X*Y); 
    jitDebug(2," sizeof(QL) (%d) sizeof(IMAGE) (%d) \n",s1,s2 ); 

    myInfo->sharemqlPtr=(QL_STRUCT *)(myInfo->qlPtr+1);
    myInfo->lkupflagentPtr= (int *)(myInfo->sharemqlPtr + myInfo->qlNum);   
  }
  else if ( obsMode==OBS_DREAM ) 
  {
    if(SC2STORE__COL_INDEX == 0)
      {
	X = COL_NUM+4;
	Y = ROW_NUM-1+4;
      }
    else
      {
	X = ROW_NUM-1+4;
	Y = COL_NUM+4;
      }
   
    int s1=sizeof(RECONSTR_STRUCT);
    int s2=sizeof(RECONSTR_IMAGE); 

    jitDebug(2,"  sizeof(double)(%d) doube[Y*X] (%d) \n",
            sizeof(double), sizeof(double)*X*Y); 
    jitDebug(2," sizeof(RESCONSTR) (%d) sizeof(IMAGE) (%d) \n",s1,s2 ); 

    //ithqlPtr+1: skip the qlrecord 
    myInfo->sharemrecnstrPtr=(RECONSTR_STRUCT *)(myInfo->qlPtr+1);
    myInfo->lkupflagentPtr= (int *)(myInfo->sharemrecnstrPtr + myInfo->qlNum);  
  }  
  else if (obsMode==OBS_SCAN) 
  {
    myInfo->sharemscanPtr=(QL_SCAN *)(myInfo->qlPtr+1);
    myInfo->lkupflagentPtr= (int *)(myInfo->sharemscanPtr +myInfo->qlNum);

    // initialised scanName=NULL, lkupflag=0
    for (i=0; i<myInfo->qlNum; i++)
    {
      scanPtr=myInfo->sharemscanPtr+i;
      scannumPtr = (int *)scanPtr;
      timePtr =(double *)(scannumPtr+1);
      scanfileName=(char*) (timePtr+1);
      strcpy(scanfileName,"NULL\n");
    }
  }
  else
  {
    *status=DITS__APP_ERROR;  
    ErsRep(0,status, "utils_alloc_memory: don't know how to do for this OBSMODE");
    return;
  }
  // move qlNum away  
  myInfo->jcmtheadEntry=(JCMTState *)(myInfo->lkupflagentPtr + myInfo->qlNum); 

  for (i=0; i<myInfo->qlNum; i++)
  {
    lkupPtr=myInfo->lkupflagentPtr +i;
    *lkupPtr=0;
  }
}


/**
 * \fn void utils_set_memory(SDSU_CONTEXT    *con,                    
 *  dasInfoStruct_t *myInfo, int whichShared, int sharedmemSize,
 *   StatusType *status)
 *
 * \brief function
 *  get shared memory 
 *
 * \param con         SDSU context structure
 * \param myInfo     dasInfo structure pointer 
 * \param whichShared  int 
 * \param sharedmemSize int 
 * \param status      StatusType pointer
 * 
 * Used to be sc2dalib_setsharedMem
 */
/*+ utils_set_memory
*/
void utils_set_memory
(                  
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
int             whichShared,
int             sharedmemSize,  
StatusType      *status
)
{
  int   key;
  int   sharedmId;
  char *sharedmPtr;
  char  whichName[30];

  if (*status != STATUS__OK) return;

  errno=0;
  if ( whichShared==SHAREDM_PAR)
  {
    key = ftok(PAR_SHARED1, PAR_SHARED2);   // Get a key 
    strcpy(whichName,"PAR");
  }
  else
  {
    key = ftok(OBS_SHARED1, OBS_SHARED2);   // Get a key 
    strcpy(whichName,"OBS");
  }
  if (key ==(key_t)-1 )
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"utils_set_memory:failed to get key ");
    con->process.reason=FRAME_MALLOCFAILED;
    return;
  }
  jitDebug(2,"utils_set_memory: (%s)  shareMemSize=%d byte\n",
         whichName,sharedmemSize);

  // if a shared-memory segment exists, get it; otherwise, create one 
  sharedmId = shmget(key, sharedmemSize, 0666 | IPC_CREAT);
  if ( sharedmId < 0)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"utils_set_memory:failed to get shared memory");
    ErsRep(0,status,"utils_set_memory: request MemSize(%f)M   SHMMAX=(?M)",
               (double)sharedmemSize/(1024*1024));
    return;
  }
  if ( whichShared==SHAREDM_PAR)
    myInfo->sharedmparId=sharedmId;
  else
    myInfo->sharedmId=sharedmId;

  // Attach segment to process. Use an attach address of zero to
  // let the system find a correct virtual address to attach.
  sharedmPtr= shmat(sharedmId, 0, 0644);
  if (sharedmPtr == (char *) -1)
  {
    *status=DITS__APP_ERROR;
    ErsRep(0,status,"utils_set_memory:failed to attach to shared memory");
    con->process.reason=FRAME_MALLOCFAILED;
    return;
  }
  // this shall set initial value for the shared memory
  memset(sharedmPtr, '0', sharedmemSize);

  if ( whichShared==SHAREDM_PAR)
  {  
     myInfo->parShm=sharedmPtr;    
     myInfo->parshmPtr=(PAR_SHARED *)myInfo->parShm;
  }
  else
    myInfo->sharedShm=sharedmPtr;
}

/**
 * \fn void utils_close_memory(dasInfoStruct_t *myInfo,
 *   int whichShared, StatusType *status)
 *
 * \brief founction
 *  close the shared memory (QL) used for drama task and data handle task
 *
 * \param myInfo       dasInfoStruct_t pointer
 * \param whichShared  int 
 * \param status       StatusType. pointer  given and return
 *
 * Used to be sc2dalib_closesharedMem
 */
/*+ utils_close_memory   
 */
void utils_close_memory
(
dasInfoStruct_t    *myInfo,
int                whichShared,
StatusType         *status
)
{
  int   sharedmId;
  char *sharedmPtr;

  struct shmid_ds     stShmId;

  if (*status != STATUS__OK) return;
  
  if ( whichShared==SHAREDM_PAR)
    sharedmPtr=myInfo->parShm;
  else
    sharedmPtr=myInfo->sharedShm;
  errno=0;
  // check if the shared memory segment attached to process.
  if (sharedmPtr != (char *) -1)
  {
    if( shmdt(sharedmPtr)==-1)
    {
      *status = DITS__APP_ERROR;
      if ( whichShared==SHAREDM_PAR)
      {
        ErsRep (0, status,"utils_close_memory: (PAR)shmdt %s",
                strerror(errno));
      }
      else
      {
        ErsRep (0, status,"utils_close_memory: (OBS)shmdt %s",
                strerror(errno));
      }
      return;
    }
  }
  // check if the shared memory segment exists.
  if ( whichShared==SHAREDM_PAR)
    sharedmId=myInfo->sharedmparId;
  else
    sharedmId=myInfo->sharedmId;

  errno=0;
  if ( sharedmId >0 )
  {
    if (shmctl(sharedmId, IPC_RMID, &stShmId))
    {
      *status = DITS__APP_ERROR;
      if ( whichShared==SHAREDM_PAR)
        ErsRep (0, status,"utils_close_memory: (PAR)shmctl %s",
                strerror(errno));
      else
        ErsRep (0, status,"utils_close_memory: (OBS)shmctl %s",
                strerror(errno));
      return;
    }
  }
}

/**
 * \fn void utils_init_variables(SDSU_CONTEXT *con, dasInfoStruct_t *myInfo, 
 *     long cols, long rows, int *pixelMask, StatusType *status)
 *
 * \brief function
 *  set some initial values for con and myInfo 
 *
 * \param con      SDSU_CONTEXT pointer
 * \param myInfo  dasInfoStaruct pointer
 * \param cols     long: initial column number
 * \param rows     long: initial row number
 * \param pixelMask int pointer
 * \param status   StatusType pointer. 
 *
 * if  *status != STATUS__OK, report error.
 *
 * Used to be sc2dalib_variablesInit
 */
/*+ utils_init_variables
*/
void utils_init_variables
(
SDSU_CONTEXT    *con,
dasInfoStruct_t *myInfo,
long            cols,
long            rows,
int             *pixelMask,
StatusType      *status
)
{
  int    i, j, pixel;

  if (*status != STATUS__OK) return;

  // use the stripchFlag here as time variable for stripchart
  stripchFlag=0;

  for (j=0; j<41; j++ )
  {
    for ( i=0; i<32; i++)
    {
      pixel=j*32 +i;
      pixelMask[pixel]=1;
    }
  }

  //FRAMEHEADER_NUM, CHKSUM_NUM defined in dasDrama_par.h
  myInfo->bufSize = (cols*rows + FRAMEHEADER_NUM + CHKSUM_NUM)*4;
  con->process.framebufsize = myInfo->bufSize;

  con->process.debuglvl[DAS]       =(char)( myInfo->debuglvl     &0x000000FF);
  con->process.debuglvl[DATAHANDLE]=(char)((myInfo->debuglvl>> 8)&0x000000FF);
  con->process.debuglvl[INTERFACE] =(char)((myInfo->debuglvl>>16)&0x000000FF);
  con->process.debuglvl[DRIVER]    =(char)((myInfo->debuglvl>>24)&0x000000FF);

  myInfo->debuglvl= (long)con->process.debuglvl[DAS];
  myInfo->doneReadxml=0;
  myInfo->fpLog=NULL;
  myInfo->fpMcecmd=NULL;
  myInfo->fpData=NULL;
  myInfo->fpBatch=NULL;
  myInfo->fpStrchart=NULL;
  myInfo->logfileFlag=0;
  myInfo->actIndex = -1;
  myInfo->actionFlag=NONEACTION;
  myInfo->glbCount=0;
  myInfo->trkNo=0;
  // SdsIdType is long
  myInfo->qlId=0;
  myInfo->seqId=0;
  myInfo->timeId=0;
  myInfo->filenameId=0;
  myInfo->imageId=0;
  myInfo->fitsId=0;
  myInfo->qldataId=0;
  myInfo->dataFormat=MCE_BINARY_FORM;

  con->gofailflag=SDSU_FALSE;
  con->dataform=MCE_BINARY_FORM;
  con->process.seqstatus=SEQ_NOACTION;
  con->process.framesetup=DA_MCE_NONE;
  con->process.reason=FRAME_WAITING;

  con->process.exit=0;
  con->process.svflag=0;
  con->process.procnum=200;
  con->process.totalframe=0;
  con->process.msgFIFO=ERR_NONE;
  con->process.wrFIFO=FIFO_ERR_NONE;
  con->process.whereabout=WAIT_OBEY;
  
  // parameter shared memory
  utils_set_memory(con,myInfo,SHAREDM_PAR,sizeof(PAR_SHARED),status);   
  if (!StatusOkP(status))
    return;
  // parameter shared Mempry
  myInfo->parshmPtr->obsMode=OBS_STARE;  // default mode
  myInfo->parshmPtr->procNo = 200;       // default
}

/**
 * \fn void utils_get_time(short which, char *dateArray)
 *
 * \brief function
 *  get the time date in DayMonthYear:hour:minute:second 
 *  (DAS_DATETIME) or DayMonthYear (DAS_DATE)format and return 
 *
 * \param which     short: either date only or data and time
 * \param dateArray char: pointer for the date string
 *
 * \retval 1: DAS_DATE 0: otherwise
 *
 * Used to be sc2dalib_finddateTime
 */
/*+ int utils_get_time 
 */
void utils_get_time
(
short         which,
char          *dateArray
)
{
  time_t tm;
  struct tm *ptr;
  char   date[40]="";
  
  tm = time(NULL);
  ptr = localtime(&tm);

  if(which==DAS_DATE)
  {
    strftime(date,40,"%Y%m%d",ptr);
    strcpy(dateArray,date);
  }
  else
  {
    strftime(date,40,"%Y%m%d:%X",ptr);
    strcpy(dateArray,date); 
  }
  return;
}



/**
 * \fn void utils_array_to_file(FILE *fp,  char *itemName, void *array,
 *    int flag, StatusType *status)
 *
 * \brief function
 *  write item  to middle setup file 
 *
 * \param fp      FILE pointer
 * \param item    string
 * \param array   void pointer
 * \param flag    int  FLOAT_NUMBER, INT_NUMBER 
 * \param dim     int  array dim
 * \param status   StatusType pointer.  given and return
 *
 * Used to be sc2dalib_writearray2setupFile
 */
/*+ utils_array_to_file
*/
void utils_array_to_file
(
FILE       *fp,
char       *item,
void       *array,
int        flag,
StatusType *status
)
{
  int   col, card,i;
  int   *intArray=NULL;
  float *floatArray=NULL;

  if ( !StatusOkP(status) ) return;

  if (flag==INT_NUMBER)
    intArray=(int*)array;
  else
    floatArray=(float*)array;

  col=0;
  for (card=1; card<=4; card++)
  {
    fprintf(fp,"#RC%1d\n%s",card,item);
    for ( i=0; i<8; i++)
    {
      if (flag==INT_NUMBER)
      { 
        if(  intArray[col]==VAL__BADI )
          fprintf(fp,"0  ");
        else
          fprintf(fp,"%4d ", intArray[col]);
      }
      else
      { 
        if(  floatArray[col]==VAL__BADI )
          fprintf(fp,"0  ");
        else
         fprintf(fp,"%4.2f ",floatArray[col]);
      }
      col ++;
    }
    fprintf(fp,"\n");
  }
}

/**
 * \fn void utils_create_SDP(dasInfoStruct_t  *myInfo, StatusType *status)
 *
 * \brief functaion
 *  create parameter SDP struct to store local info used by this task.
 *
 * \param myInfo  dasInfoStruct_t pointer
 * \param status   StatusType pointer.  given and return
 *
 * Used to be sc2dalib_createSDP
 */
/*+ utils_create_SDP
 */
void utils_create_SDP
(
dasInfoStruct_t  *myInfo,
StatusType       *status
)
{   
  int   allzeros=0;
  char  noname[]=" ";

  //have to be static
  static long           zero=0;
  static double         dzero=0.0;  

  
  if (!StatusOkP(status)) return;

  // parameter supported by this task 
  static SdpParDefType  daParDefs[]=
  {
   { "STATE", &zero,  SDS_STRUCT  },
   { "CONFIGURE_ID", &zero,  SDS_INT  },
   { "SETUP_SEQ_ID", &zero,  SDS_INT  },
   { "SEQUENCE_ID", &zero,  SDS_INT  },
   { "TAI", &dzero,  SDS_DOUBLE  },
   { "INITIALISED", &zero, SDS_INT},
   { "CONFIGURED", &zero, SDS_INT},
   { "SETUP", &zero, SDS_INT},
   { "IN_SEQUENCE", &zero, SDS_INT},
   { "SEQ_START", &zero, SDS_INT},
   { "SEQ_END", &zero, SDS_INT},
   { "SEQ_DWELL", &zero, SDS_INT},
   { "FRAME_RATE", &zero, SDS_INT},
   { "BATCH_DELAY", &zero, SDS_INT},
   { "QL", &zero, SDS_STRUCT},
   { "DATA_FILE", "A long filename", ARG_STRING},
   { "FRAME_NO2PROC", &zero, SDS_INT},
   { "FRAME_TOTALNO", &zero, SDS_INT},
   { "FRAME_BUFSIZE", &zero, SDS_INT},
   { "ARRAY_CNNTED", "array name", ARG_STRING}, 
   { "WAVELENGTH", "850", ARG_STRING}, 
   { "OBS_NO", &zero, SDS_INT},
   { "SUBSCAN_NO", &zero, SDS_INT},
   { "DATE", "A long filename", ARG_STRING},
   { "FRAME_FILENAME", "A long filename", ARG_STRING},
   { "MCE_WHICHCARD", "A card name", ARG_STRING},
   { "SHUTTER_STATUS", "UNKNOWN", ARG_STRING},
   { "OBS_MODE", "DARK", SDP_STRING},
   { "STAREADD", &zero, SDS_INT},
   { "SCANWRITE", &zero, SDS_INT},
   { "DREAMROUND", &zero, SDS_INT},
   { "HTTRACK_FLAG", &zero, SDS_INT},
   { "INITFLAGS", &zero, SDS_INT},
   { SC2DALOGFILE, "A long filename", ARG_STRING},
   { SC2DAMCEXMLFILE, "A long filename", ARG_STRING},
   { SC2DADATAFORMAT, " format", ARG_STRING}
  };
 
  myInfo->parsysId = (SdsIdType)DitsGetParId();
  SdpCreate(myInfo->parsysId,DitsNumber(daParDefs), daParDefs, status);
  if( *status != STATUS__OK)
  {
    printf("utils_create_SDP: SdpCreate failed \n");
    return;
  }

  // put all the parameters into the STATE structure */
  SdpGetSds("STATE", &myInfo->statId, status);
  if( *status != STATUS__OK)
  {
    printf("utils_create_SDP:could not initialized STATE sds structure \n");
    return;
  }
  ArgPuti(myInfo->statId, "NUMBER",allzeros, status);
  ArgPuti(myInfo->statId, "SEQ_START",allzeros, status);
  ArgPuti(myInfo->statId, "SEQ_END",allzeros,status);
  ArgPuti(myInfo->statId, "SEQ_DWELL",allzeros,status);
  ArgPuti(myInfo->statId, "FRAME_RATE",1,status);
  ArgPuti(myInfo->statId, "BATCH_DELAY",0,status);
  ArgPutString(myInfo->statId, "OBS_MODE",noname, status);
  ArgPutString(myInfo->statId, "MCE_WHICHCARD",noname, status);
  ArgPutString(myInfo->statId, "ARRAY_CNNTED",noname, status);
  ArgPutString(myInfo->statId, "DATE",noname, status);
  ArgPutString(myInfo->statId, "WAVELENGTH",noname, status);
  ArgPuti(myInfo->statId, "OBS_NO",allzeros, status);
  ArgPuti(myInfo->statId, "SUBSCAN_NO",allzeros, status);
  ArgPuti(myInfo->statId, "STAREADD",200, status);
  ArgPuti(myInfo->statId, "SCANWRITE",40000, status);
  ArgPuti(myInfo->statId, "DREAMROUND",3, status);
  // leave create QL SDS struct in configuration since it depends on 
  // obs_mode
}
   

/**
 * \fn void utils_create_QLSDS(dasInfoStruct_t  *myInfo,
 *  StatusType *status)
 *
 * \brief functaion
 *  create QL parameter structure depending on obsmode 
 *
 * \param myInfo  dasInfoStruct_t pointer
 * \param status   StatusType pointer.  given and return
 *
 * Used to be sc2dalib_createQLSDS
 */
/*+ utils_create_QLSDS
 */
void utils_create_QLSDS
(
dasInfoStruct_t  *myInfo,
StatusType       *status
)
{
  int   obsMode;
  char  scanfile[]="                     ";
  long tLong;
  double tDouble;

  if (*status != STATUS__OK) return;

  /* Get the ID of the QL structure */

  SdpGetSds("QL", &myInfo->qlId, status);
  if( *status != STATUS__OK)
  {
    ErsRep(0,status,"utils_create_QLSDS: SdpGetSds QL failed \n");
    return;
  }

  /* Mark the QL structure as not containing valid data */
  ArgPutString(myInfo->qlId,"DATAVALID", "NO", status);

  // delete and free unsed sub sds data from qlId

  sc2dalib_chksdsID(myInfo,status);
  if( *status != STATUS__OK)
  {
    ErsRep(0,status,"utils_create_QLSDS: sc2dalib_chksdsID failed \n");
    return;
  }

  // Place the parameters that are shared into the QL structure 
  tLong =  0;
  ArgPuti(myInfo->qlId,"FRAMENUM", tLong, status);
  SdsFind(myInfo->qlId,"FRAMENUM", &myInfo->seqId, status);
  if( *status != STATUS__OK)
  {
    ErsRep(0,status,"utils_create_QLSDS: ArgPuti (framenum) failed \n");
    return;
  } 
  tDouble = 0.0;
  ArgPutd(myInfo->qlId,"TIMESTAMP", tDouble, status);
  SdsFind(myInfo->qlId,"TIMESTAMP", &myInfo->timeId, status);
  if( *status != STATUS__OK)
  {
    ErsRep(0,status,"utils_create_QLSDS: ArgPutd (timestamp) failed \n");
    return;
  } 

  obsMode=myInfo->parshmPtr->obsMode;

  /* The FILENAME is only used in OBS_SCAN */  
  if( obsMode == OBS_SCAN )
  {
    ArgPutString(myInfo->qlId,"FILENAME",scanfile,status);
    SdsFind(myInfo->qlId,"FILENAME", &myInfo->filenameId, status);
    if( *status != STATUS__OK)
    {
      ErsRep(0,status,"utils_create_QLSDS: argputstring (scanfile) failed");
      return;
    }

    qlFileCounter=1;
  }
  else if( obsMode==OBS_STARE || obsMode==OBS_DREAM)
  {
    unsigned long  datadims[2];
    unsigned long  fitsdims[2];

    // Create a x*y array; y=datadims[0]=COL_NUM; x=datadims[1]=ROW_NUM-1;
    if ( obsMode !=OBS_DREAM )
    {
      datadims[SC2STORE__COL_INDEX] = COL_NUM; 
      datadims[SC2STORE__ROW_INDEX] = ROW_NUM-1;
    }
    else
    {
      if(SC2STORE__COL_INDEX == 0)
      {
        datadims[0] = COL_NUM+4;
        datadims[1] = ROW_NUM-1+4;
      }
      else
      {
        datadims[0] = ROW_NUM-1+4;
        datadims[1] = COL_NUM+4;
      }
    }
   
    // Both the data array and the FITS headers go in the structure called "IMAGE"
    SdsNew(myInfo->qlId,"IMAGE", 0, NULL, SDS_STRUCT, 0, NULL, 
           &myInfo->imageId,status);
    if( *status != STATUS__OK)
    {
      ErsRep(0,status,"utils_create_QLSDS: SdsNew (image) failed");
      return;
    } 
    // Stick the space for the FITS headers into IMAGE 
    fitsdims[0] = FITSSIZE;
    fitsdims[1] = MAXFITS;
    SdsNew(myInfo->imageId, "FITS", 0, NULL, SDS_CHAR, 2, fitsdims, 
           &myInfo->fitsId,status);
    if( *status != STATUS__OK)
    {
      ErsRep(0,status,"utils_create_QLSDS: SdsNew (fits) failed ");
      return;
    }

    // Stick the DATA_ARRAY into IMAGE
    SdsNew(myInfo->imageId,"DATA_ARRAY", 0, NULL, SDS_DOUBLE, 2, datadims, 
           &myInfo->qldataId,status);
    if( *status != STATUS__OK)
    {
      ErsRep(0,status,"utils_create_QLSDS: SdsNew (data-array)failed ");
      return;
    }
  }
}  


/**
 * \fn void utils_update_QLSDS(SDSU_CONTEXT *con, dasInfoStruct_t  *myInfo,
 * int *coaddnumPtr, double *timePtr, char  *fitsPtr, uint32 *fitsDim,
 * double *coadddataPtr, uint32 *dataDim, char  *scanfileName,
 *  StatusType *status)
 *
 * \brief functaion
 *  resize QL  structure depending on realdata 
 *
 * \param con          SDSU context structure pointer
 * \param myInfo      dasInfoStruct_t pointer
 * \param coaddnumPtr  int pointer for ith codd Number,
 * \param timePtr      double pointer for time
 * \param fitsPtr      char pointer for FITS 
 * \param fitsDim      uint32 pointer, fits dimension [0] size_n, [1] num_n
 * \param coadddataPtr  double pointer for coadd data entry,
 * \param dataDim      uint32 pointer,  data dimension [0] col_n, [1] row_n
 * \param scanfileName char pointer for scanfile name 
 * \param status   StatusType pointer.  given and return
 *
 * Used to be sc2dalib_updateQLSDS
 */
/*+ utils_update_QLSDS
 */
void utils_update_QLSDS
(
SDSU_CONTEXT     *con,
dasInfoStruct_t  *myInfo,
int              *coaddnumPtr,
double           *timePtr,
char             *fitsPtr,
uint32           *fitsDim,    // fits dimension [0] size_n, [1] num_n
double           *coadddataPtr,
uint32           *dataDim,      // data dimension [0] col_n, [1] row_n
char             *scanfileName,
StatusType       *status
)
{   
  int   bufSize;
  int   obsMode;
  static char dummyFileName[FILE_LEN + 4];

  if (*status != STATUS__OK) return;

  obsMode=myInfo->parshmPtr->obsMode;

  SdsPut(myInfo->timeId,sizeof(double),0,timePtr,status);
  ArgPuti(myInfo->qlId,"OBSNUM", (long)myInfo->obsNo, status);
  
  if(obsMode==OBS_STARE || obsMode==OBS_DREAM)
  {
    if(myInfo->parshmPtr->load != LOAD_DARK)
      {

  // Put the QL frame number in FRAMENUM
  SdsPut(myInfo->seqId, sizeof(int),0,coaddnumPtr,status);

  // Create a dataDim[0]*[1] array within dasInfo.qlId + others 
  bufSize=dataDim[0]*dataDim[1]*sizeof(double);
  // update now
  SdsResize(myInfo->fitsId,2,fitsDim,status);
  SdsPut(myInfo->fitsId,fitsDim[0]*fitsDim[1],0,fitsPtr,status);
  // the coaddaataPtr is the coadd frame   
  SdsPut(myInfo->qldataId,bufSize,0,coadddataPtr,status);

  /* Mark the QL structure as containing valid data
     and then update the QL structure */
  
  ArgPutString(myInfo->qlId,"DATAVALID", "YES", status);
  SdpUpdate(myInfo->qlId,status);

      }

  }
  else if( obsMode==OBS_SCAN )
  {
    /* Add .sdf to the end of the scanfileName */
    strncpy(dummyFileName, scanfileName, FILE_LEN);
    strcat(dummyFileName,".sdf");
    ArgPutString(myInfo->qlId,"FILENAME",dummyFileName, status);

    /* Log the filename to the log file */
    fprintf(myInfo->fpLog,
     "utils_update_QLSDS: Putting this file name into the QL structure %s\n",
     dummyFileName);

    /* Increment counter so it increments with each new file (and put in FRAMENUM)*/

    ArgPuti(myInfo->qlId, "FRAMENUM", qlFileCounter, status);

  /* Mark the QL structure as containing valid data
     and then update the QL structure */
  
    ArgPutString(myInfo->qlId,"DATAVALID", "YES", status);
    SdpUpdate(myInfo->qlId,status);

    qlFileCounter++;

  }
} 

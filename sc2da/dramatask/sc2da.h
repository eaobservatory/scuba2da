#ifndef HEADGEN____sc2da_h
#define HEADGEN____sc2da_h 
 
 
/*+sc2da- main routine  extern int main
*/
JIT_MAIN (sc2da)
( 
int argc, 
char **argv
);

/*+ sc2da_Abort - 
*/
void sc2da_Abort
(
StatusType *status
);

/*+ sc2da_Batch -
*/
void sc2da_Batch
(
StatusType *status
);

/*+ sc2da_Clear -
*/
void sc2da_Clear
(
StatusType *status
);

/*+ sc2da_Config - 
*/
void sc2da_Config
(
StatusType *status
);

/*+ sc2da_ConfigKick - Handles a "kick" of the configure action 
*/
void sc2da_ConfigKick
( 
StatusType *status    /* global status (given and returned) */
);

/*+ sc2da_Downld2PCI
*/
void sc2da_Downld2PCI
(
StatusType *status
);

/*+ sc2da_Down2FPGA -
*/
void sc2da_Down2FPGA
(                    
StatusType *status
);

/*+ sc2da_Dispinfo - 
*/
void sc2da_Dispinfo
(
StatusType *status    
);

/*+ sc2da_EndObs -
*/
void sc2da_EndObs 
( 
StatusType *status    /* global status (given and returned) */
);

/*+ sc2da_ObsKick -
*/
void sc2da_ObsKick 
( 
StatusType *status    /* global status (given and returned) */
);

/*+ sc2da_Exit 
*/
void sc2da_Exit 
(
StatusType *status 
);

/*+ sc2da_Getarrayname -
*/
void sc2da_GetarrayName
(
StatusType *status
);

/*+ sc2da_Heatslope
*/
void sc2da_Heatslope
( 
StatusType *status    /* global status (given and returned) */
);

/*+ sc2da_Heatermask
*/
void sc2da_Heatermask
(
StatusType *status
);

/*+ sc2da_Heatersloperead 
*/
void sc2da_Heatersloperead
( 
StatusType *status    
);

/*+ sc2da_Init -
*/
void sc2da_Init
(
StatusType *status
);

/*+ sc2da_initKick - Handles a "kick" of the initialise action 
*/
void sc2da_initKick
( 
StatusType *status    /* global status (given and returned) */
);

/*+ sc2da_Mcecmd - 
*/
void sc2da_Mcecmd
(
StatusType *status
);

/*+ sc2da_Mceonflycmd - 
*/
void sc2da_Mceonflycmd
(
StatusType *status
);

/*+ sc2da_MceStatus - 
*/
void sc2da_MceStatus
(
StatusType *status
);

/*+  sc2da_PCIcmd- 
*/
void sc2da_PCIcmd
(
StatusType *status
);

/*+  sc2da_PCIblk- 
*/
void sc2da_PCIblk
(
StatusType *status
);

/*+ sc2da_Ping
*/
void sc2da_Ping
( 
StatusType *status  
);

/*+ sc2da_Report 
*/
void sc2da_Report
( 
StatusType *status   
);

/*+ sc2da_Pixelmon 
*/
void sc2da_Pixelmon
( 
StatusType *status    
);

/*+ sc2da_pixelKick - 
*/
void sc2da_pixelKick
( 
StatusType *status    
);

/*+ sc2da_Pixelmask -
*/
void sc2da_Pixelmask
(
StatusType *status
);

/*+ sc2da_SetSeq     
*/
void sc2da_SetSeq
(
StatusType *status
);

/*+ sc2da_SetSeqKick 
*/
void sc2da_SetSeqKick
( 
StatusType *status    /* global status (given and returned) */
);

/*+ sc2da_Seq - 
*/
void sc2da_Seq
(
StatusType *status
);

/*+ sc2da_SeqKick -
 */
void sc2da_SeqKick
(
StatusType *status
);

/*+ sc2da_Servo - 
*/
void sc2da_Servo
( 
StatusType *status    /* global status (given and returned) */
);

/*+ sc2da_servoKick - 
*/
void sc2da_servoKick
( 
StatusType *status    
);

/*+ sc2da_SetengData
*/
void sc2da_SetengData
(
StatusType            *status
);

/*+ sc2da_SQ1optpts
*/
void sc2da_SQ1optpts
(
StatusType *status
);

/*+ sc2da_Trksq2fb  
*/
void sc2da_Trksq2fb
( 
StatusType *status    
);

/*+ sc2da_Trkheat - 
*/
void sc2da_Trkheat
( 
StatusType *status    
);

/*+ sc2da_trkKick - 
*/
void sc2da_trkKick
( 
StatusType *status    
);

/*+ sc2da_UpdateLogname - 
*/
void sc2da_UpdateLogname
(
StatusType *status
);

/*+ sc2da_Version - 
*/
void sc2da_Version
(
StatusType *status
);

/*+ mceTryChkSum - 
*/
void mceTryChkSum
(
StatusType *status
);

/*+ mceTest - 
*/
void mceTest
(
StatusType *status
);

/*+ sc2da_Getsq2fbPara  
*/
void sc2da_Getsq2fbPara
( 
StatusType *status    
);

 
 
#endif

/**
 * \file sc2da_err.h
 *
 * \brief error message define for SCUBA2 Data Acquisition Software (DAS)  
 *
 * \author Xiaofeng Gao, UKATC (xg@roe.ac.uk)
 *
 * \version 
 *
 * \date 22/09/2005
 *
 * Copyright (c) 2005.  UK Astronomy Technology Centre (UK ATC), 
 * An establishment of the Particle Physics and Astronomy Research Council
 * PPARC).
 *
 * Web: www.roe.ac.uk
 *
 *  $Log: sc2da_err.h,v $
 *  Revision 1.1.1.1  2007/05/16 08:27:01  dkelly
 *  first insertion
 *
 *
 */


/* MCE error numbers */
typedef enum 
{
   NON_ERROR      =0,    ///=0 
   CHKSUM_ERROR ,         ///=1,< fiber checksum error (command packet)
   NO_CARD_ERR,
   CRC_ERR,
   READONLY_ERR
}MCEERROR;


char *mceErrors[] =
{
"NON_ERROR       ///<no error",
"CHKSUM_ERROR    ///< fiber checksum error (command packet)",
"Card not present in the subrack",                                                                                                                                                                                              
"CRC error is detected, or waiting reply times out",                                                                                                                                                                                        
"a read-only slave has been written to",        
" "   
};


#define MCE__OK            0
#define MCE__ERROR        -200
#define MCE__SDSUCMD      -201
#define MCE__LASRBITSET   -202
#define MCE__STOPBITSET   -203
#define MCE__FREEBUF      -204
#define MCE__TANSLATE     -205


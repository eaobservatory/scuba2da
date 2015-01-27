#ifndef sc2da_mce_h
#define sc2da_mce_h


#include "sc2da_inc.h"

void mce_read
(
	SDSU_CONTEXT          *con,      
	dasInfoStruct_t       *myInfo,
	struct mcexml_struct  *mceinxPtr,
	char                  *cmd,
	int                   *val,
	int                   howMany,
	StatusType            *status
);

void mce_set
(
	SDSU_CONTEXT          *con,      
	dasInfoStruct_t       *myInfo,
	struct mcexml_struct  *mceinxPtr,
	char                  *cmd,
	StatusType            *status
);

void mce_read_frame_rate
(
	SDSU_CONTEXT          *con,
	dasInfoStruct_t       *myInfo,
	struct mcexml_struct  *mceinxPtr,
	long                  *rate,
	StatusType            *status
);

void mce_checksum
(
	dasInfoStruct_t  *myInfo,   
	dasCmdInfo_t     *mycmdInfo,
	StatusType       *status
);

void mce_results
(
	dasInfoStruct_t       *myInfo,   
	dasCmdInfo_t          *mycmdInfo, 
	StatusType            *status
);

void mce_reply
(
	dasCmdInfo_t  *mycmdInfo,
	int           howMany, 
	StatusType    *status
);

void mce_error_reply
(
dasInfoStruct_t *myInfo,
dasCmdInfo_t    *mycmdInfo,
StatusType      *status
);

void mce_cmd
(
SDSU_CONTEXT      *con,
dasInfoStruct_t   *myInfo,
dasCmdInfo_t      *mycmdInfo,
char              *dateTime,
StatusType        *status
);

#endif

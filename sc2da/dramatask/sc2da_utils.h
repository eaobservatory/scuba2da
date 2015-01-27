#ifndef sc2da_utils_h
#define sc2da_utils_h

#include "sc2da_inc.h"

void utils_msg
(
	dasInfoStruct_t *myInfo, 
	char            *string, 
	char            *string2,
	int             flag,
	StatusType      *status
);

void utils_fclose(FILE **fp);

void utils_close_files(dasInfoStruct_t *myInfo);

void utils_open_files
(
	dasInfoStruct_t *myInfo,       
	int             getData,
	int             getBatch,
	StatusType      *status
);

void utils_update_debug
(
	SDSU_CONTEXT    *con,
	dasInfoStruct_t *myInfo,
	StatusType      *status    
);

void utils_alloc_memory
(
	SDSU_CONTEXT       *con,         
	dasInfoStruct_t    *myInfo,
	int                startSeq,
	int                endSeq,
	StatusType         *status
);

void utils_set_memory
(                  
	SDSU_CONTEXT    *con,
	dasInfoStruct_t *myInfo,
	int             whichShared,
	int             sharedmemSize,  
	StatusType      *status
);

void utils_close_memory
(
	dasInfoStruct_t    *myInfo,
	int                whichShared,
	StatusType         *status
);

void utils_init_variables
(
	SDSU_CONTEXT    *con,
	dasInfoStruct_t *myInfo,
	long            cols,
	long            rows,
	int             *pixelMask,
	StatusType      *status
);

void utils_get_time
(
short         which,
char          *dateArray
);

void utils_array_to_file
(
FILE       *fp,
char       *item,
void       *array,
int        flag,
StatusType *status
);

void utils_create_SDP
(
dasInfoStruct_t  *myInfo,
StatusType       *status
);

void utils_create_QLSDS
(
dasInfoStruct_t  *myInfo,
StatusType       *status
);

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
);

#endif

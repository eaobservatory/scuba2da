#ifndef HEADGEN____sc2da_heat_h
#define HEADGEN____sc2da_heat_h 

#include "sc2da_par.h"
#include "sc2da_struct.h"
#include "sc2math.h"   // for sc2math-linfit
#include "sc2headman.h"
#include "math.h"      // for pow
#include "jitXML.h"
#include "sc2sqopt_par.h"
#include "sc2sqopt.h"
#include "star/atl.h"
#include "sc2dalib.h"
#include "sc2headman.h"
#include "sc2dalibsetup.h"

/*+ heat_init_bias
*/
void heat_init_bias
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,   
int                   *flag,
StatusType            *status
);

/*+ heat_init_slope
*/
void heat_init_slope
(
SDSU_CONTEXT      *con,
dasInfoStruct_t   *myInfo, 
char              *dateTime,
StatusType         *status
);

/*+ heat_slope
*/
void heat_slope
(
dasInfoStruct_t  *myInfo,
int              *heaterMask,
double           *heaterSlope,
StatusType       *status
);

/*+ heat_read_slope
*/
void heat_read_slope
(
dasInfoStruct_t   *myInfo,    
double            *pixelSlope,
StatusType        *status
);

/*+ heat_save_slope
*/
void heat_save_slope
(
dasInfoStruct_t *myInfo,
double          *heaterSlope,
StatusType      *status

);

/*+ heat_init_track
*/
void heat_init_track
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,    
dasCmdInfo_t          *myCmd,
struct mcexml_struct  *mceinxPtr,
char                  *dateTime,
StatusType            *status
);

/*+ heat_update_track
 */
void heat_update_track
(
SDSU_CONTEXT         *con,         
dasInfoStruct_t      *myInfo,
struct mcexml_struct *mceInxpt,
int                  *data,
double               *slope,
StatusType           *status
);

/*+ heat_pixel_track
 */
void heat_pixel_track
(
dasInfoStruct_t      *myInfo,
int                  heat,
int                  *data,
int                  *fdbkRef,
int                  numActivePixels,
int                  *activePixels,
double               *slope,
double               *heatVal,
int                  first,
int                  *pixelData,
StatusType           *status
);

/*+ heat_step_current
*/
void heat_step_current
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,   
StatusType            *status
);

/*+ heat_step_TES
*/
void heat_step_TES
(
SDSU_CONTEXT          *con,      
dasInfoStruct_t       *myInfo,
struct mcexml_struct  *mceInxpt,   
StatusType            *status
);
 
#endif

#ifndef MAIN_SIMULINK_cap_host_h__
#define MAIN_SIMULINK_cap_host_h__
#ifdef HOST_CAPI_BUILD
#include "rtw_capi.h"
#include "rtw_modelmap_simtarget.h"
typedef struct { rtwCAPI_ModelMappingInfo mmi ; }
MAIN_SIMULINK_host_DataMapInfo_T ;
#ifdef __cplusplus
extern "C" {
#endif
void MAIN_SIMULINK_host_InitializeDataMapInfo ( MAIN_SIMULINK_host_DataMapInfo_T
* dataMap , const char * path ) ;
#ifdef __cplusplus
}
#endif
#endif
#endif

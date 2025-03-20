#ifndef Test_SIMULINK_cap_host_h__
#define Test_SIMULINK_cap_host_h__
#ifdef HOST_CAPI_BUILD
#include "rtw_capi.h"
#include "rtw_modelmap_simtarget.h"
typedef struct { rtwCAPI_ModelMappingInfo mmi ; }
Test_SIMULINK_host_DataMapInfo_T ;
#ifdef __cplusplus
extern "C" {
#endif
void Test_SIMULINK_host_InitializeDataMapInfo ( Test_SIMULINK_host_DataMapInfo_T
* dataMap , const char * path ) ;
#ifdef __cplusplus
}
#endif
#endif
#endif

#include "Test_SIMULINK_capi_host.h"
static Test_SIMULINK_host_DataMapInfo_T root;
static int initialized = 0;
__declspec( dllexport ) rtwCAPI_ModelMappingInfo *getRootMappingInfo()
{
    if (initialized == 0) {
        initialized = 1;
        Test_SIMULINK_host_InitializeDataMapInfo(&(root), "Test_SIMULINK");
    }
    return &root.mmi;
}

rtwCAPI_ModelMappingInfo *mexFunction(){return(getRootMappingInfo());}

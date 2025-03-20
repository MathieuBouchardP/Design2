#ifndef MAIN_SIMULINK_h_
#define MAIN_SIMULINK_h_
#ifndef MAIN_SIMULINK_COMMON_INCLUDES_
#define MAIN_SIMULINK_COMMON_INCLUDES_
#include <stdio.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "sigstream_rtw.h"
#include "simtarget/slSimTgtSigstreamRTW.h"
#include "simtarget/slSimTgtSlioCoreRTW.h"
#include "simtarget/slSimTgtSlioClientsRTW.h"
#include "simtarget/slSimTgtSlioSdiRTW.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "raccel.h"
#include "slsv_diagnostic_codegen_c_api.h"
#include "rt_logging_simtarget.h"
#include "rt_nonfinite.h"
#include "math.h"
#include "dt_info.h"
#include "ext_work.h"
#endif
#include "MAIN_SIMULINK_types.h"
#include <stddef.h>
#include "rtGetInf.h"
#include "rtw_modelmap_simtarget.h"
#include "rt_defines.h"
#include <string.h>
#define MODEL_NAME MAIN_SIMULINK
#define NSAMPLE_TIMES (6) 
#define NINPUTS (0)       
#define NOUTPUTS (0)     
#define NBLOCKIO (31) 
#define NUM_ZC_EVENTS (0) 
#ifndef NCSTATES
#define NCSTATES (6)   
#elif NCSTATES != 6
#error Invalid specification of NCSTATES defined in compiler command
#endif
#ifndef rtmGetDataMapInfo
#define rtmGetDataMapInfo(rtm) (*rt_dataMapInfoPtr)
#endif
#ifndef rtmSetDataMapInfo
#define rtmSetDataMapInfo(rtm, val) (rt_dataMapInfoPtr = &val)
#endif
#ifndef IN_RACCEL_MAIN
#endif
typedef struct { real_T bu3amqfjr4 ; } hnt03bld54 ; typedef struct { int32_T
bwn2zcbnfw ; boolean_T j5neklhggr ; } egc1zchle4 ; typedef struct { real_T
p0it2epuk4 ; } ib0sbtg3ox ; typedef struct { int32_T kvftsnjgid ; boolean_T
apuvtbdrc3 ; } n1bfl5llqh ; typedef struct { real_T ordo20z2um ; } m4pxvochto
; typedef struct { int32_T pgvjqvnmcd ; boolean_T kbykrgfalh ; } f1zatxozzs ;
typedef struct { real_T fsu23y5ts3 ; real_T ccwcs4orvs ; real_T exsilx1jdl ;
real_T jfqdwkeqgu ; real_T oobv21my4e ; real_T h13eabb5l1 ; real_T k1g1lcxbvk
; real_T hssf4p0lhx ; real_T eelznhvfta ; real_T fqgdvg40sk ; real_T
mzld5yp0dz ; real_T lma1ydcnif ; real_T db4m1bc3gg ; real_T olqwozfrl3 ;
real_T jl4tzlqkvr ; real_T lijzxwwv3g ; real_T eblcjzhc4e ; real_T jccqxfucp1
; real_T hxdb1eoutr ; real_T op004jo4lt ; real_T mxyvw51mkl ; real_T
hkgktvl34f ; real_T lbhnidlzjv ; real_T p2k4zdl445 ; real_T ay1xeqbo54 ;
m4pxvochto mfc24avgl1 ; ib0sbtg3ox n55rhng3cw ; hnt03bld54 arhszp1yvs ;
m4pxvochto cdh5o32feg ; ib0sbtg3ox nlg3idhhmk ; hnt03bld54 d355u2baya ; } B ;
typedef struct { real_T crbbuehctw [ 2 ] ; real_T iswnzmgxkk ; real_T
coltpdu5lu ; real_T ova1iuexhr ; real_T eblypyylx2 ; struct { real_T
modelTStart ; } pxrw0rbfcb ; struct { real_T modelTStart ; } geadk4kxpi ;
struct { void * TUbufferPtrs [ 2 ] ; } je1znfcfuw ; struct { void *
LoggedData ; } awoseuipwg ; struct { void * TUbufferPtrs [ 2 ] ; } b0bf5yxgs5
; struct { int_T Tail ; int_T Head ; int_T Last ; int_T CircularBufSize ;
int_T MaxNewBufSize ; } msak0ukcbx ; struct { int_T Tail ; int_T Head ; int_T
Last ; int_T CircularBufSize ; int_T MaxNewBufSize ; } cylxtx5seq ; int_T
anvved14ni ; int_T d0equhukvj ; boolean_T ajdhv40vzi ; boolean_T g5c5uvs55y ;
boolean_T pwovuz3v0u ; boolean_T dlr30m0qeo ; f1zatxozzs mfc24avgl1 ;
n1bfl5llqh n55rhng3cw ; egc1zchle4 arhszp1yvs ; f1zatxozzs cdh5o32feg ;
n1bfl5llqh nlg3idhhmk ; egc1zchle4 d355u2baya ; } DW ; typedef struct {
real_T h5fxknj5eo ; real_T c5tdg2xm33 [ 3 ] ; real_T fa04rowvce ; real_T
aqpnssxlaz ; } X ; typedef struct { real_T h5fxknj5eo ; real_T c5tdg2xm33 [ 3
] ; real_T fa04rowvce ; real_T aqpnssxlaz ; } XDot ; typedef struct {
boolean_T h5fxknj5eo ; boolean_T c5tdg2xm33 [ 3 ] ; boolean_T fa04rowvce ;
boolean_T aqpnssxlaz ; } XDis ; typedef struct { real_T h5fxknj5eo ; real_T
c5tdg2xm33 [ 3 ] ; real_T fa04rowvce ; real_T aqpnssxlaz ; } CStateAbsTol ;
typedef struct { real_T h5fxknj5eo ; real_T c5tdg2xm33 [ 3 ] ; real_T
fa04rowvce ; real_T aqpnssxlaz ; } CXPtMin ; typedef struct { real_T
h5fxknj5eo ; real_T c5tdg2xm33 [ 3 ] ; real_T fa04rowvce ; real_T aqpnssxlaz
; } CXPtMax ; typedef struct { real_T er4q3gqmjl ; real_T hyejaapxwj ; } ZCV
; typedef struct { rtwCAPI_ModelMappingInfo mmi ; } DataMapInfo ; struct P_ {
real_T TEC ; real_T PWM_Period ; real_T DerivativeKickDF_X0 ; real_T
DerivativeKickDF_X0_pde1xmny1t ; real_T DerivativeKickDF_X0_pcnom2kezu ;
real_T DiscreteTransferFcnwithinitialstates2_X0 ; real_T
DiscreteTransferFcnwithinitialstates1_X0 ; real_T s_Gain ; real_T
Consigner_Time ; real_T Consigner_Y0 ; real_T Consigner_YFinal ; real_T
Dynamique2_A ; real_T Dynamique2_C ; real_T gain_in_Gain ; real_T _Gain ;
real_T _Gain_j0borcplw2 ; real_T _Gain_lwe4ykerbv ; real_T
Bruitdemesure_Amplitude ; real_T Bruitdemesure_Frequency ; real_T
DiscreteStateSpace_A [ 2 ] ; real_T DiscreteStateSpace_C [ 2 ] ; real_T
DiscreteStateSpace_D ; real_T Poids1_Gain ; real_T Retard2_Delay ; real_T
Retard2_InitOutput ; real_T gain_in_Gain_noqrl1m1xd ; real_T _Gain_gnigm4iwxo
; real_T _Gain_lsczmwd1h2 ; real_T _Gain_j1zbr4lzxq ; real_T
Bruitdemesure_Amplitude_hbiskpfbn4 ; real_T
Bruitdemesure_Frequency_jbbh4d4sn1 ; real_T DiscreteStateSpace_A_mopcmmoov2 ;
real_T DiscreteStateSpace_C_j3dzy0w2di ; real_T
DiscreteStateSpace_D_fbqt1av0ft ; real_T Poids2_Gain ; real_T
FiltretripleRC_A [ 3 ] ; real_T FiltretripleRC_C [ 3 ] ; real_T
DiscreteStateSpace_A_asx1t5v32r ; real_T DiscreteStateSpace_C_jykma24tmg ;
real_T DiscreteStateSpace_D_iep130ejoo ; real_T u_op_Time ; real_T u_op_Y0 ;
real_T u_op_YFinal ; real_T DiscreteStateSpace_D_mvsskepuks ; real_T
DiscreteStateSpace_A_jfbtcvit52 ; real_T DiscreteStateSpace_C_j3s1yy5c2n ;
real_T Saturation_UpperSat ; real_T Saturation_LowerSat ; real_T TEC2_Gain ;
real_T TEC1_Gain ; real_T Dynamique1_A ; real_T Dynamique1_C ; real_T
Dynamique3_A ; real_T Dynamique3_C ; real_T Retard3_Delay ; real_T
Retard3_InitOutput ; real_T baseline_temp4_Value ; real_T offset_Value ;
real_T offset1_Value ; real_T baseline_temp1_Value ; real_T
baseline_temp2_Value ; real_T offset_Value_jz2p0lxylj ; real_T
offset1_Value_nyr0ia51bl ; real_T baseline_temp1_Value_lk15cqnntl ; real_T
baseline_temp2_Value_jauwwxrr5h ; real_T baseline_temp3_Value ; real_T
baseline_temp_Value ; uint8_T ManualSwitch1_CurrentSetting ; uint8_T
ManualSwitch_CurrentSetting ; } ; extern const char_T *
RT_MEMORY_ALLOCATION_ERROR ; extern B rtB ; extern X rtX ; extern DW rtDW ;
extern P rtP ; extern mxArray * mr_MAIN_SIMULINK_GetDWork ( ) ; extern void
mr_MAIN_SIMULINK_SetDWork ( const mxArray * ssDW ) ; extern mxArray *
mr_MAIN_SIMULINK_GetSimStateDisallowedBlocks ( ) ; extern const
rtwCAPI_ModelMappingStaticInfo * MAIN_SIMULINK_GetCAPIStaticMap ( void ) ;
extern SimStruct * const rtS ; extern DataMapInfo * rt_dataMapInfoPtr ;
extern rtwCAPI_ModelMappingInfo * rt_modelMapInfoPtr ; void MdlOutputs ( int_T
tid ) ; void MdlOutputsParameterSampleTime ( int_T tid ) ; void MdlUpdate ( int_T tid ) ; void MdlTerminate ( void ) ; void MdlInitializeSizes ( void ) ; void MdlInitializeSampleTimes ( void ) ; SimStruct * raccel_register_model ( ssExecutionInfo * executionInfo ) ;
#endif

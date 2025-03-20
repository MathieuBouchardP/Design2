#ifndef Test_SIMULINK_h_
#define Test_SIMULINK_h_
#ifndef Test_SIMULINK_COMMON_INCLUDES_
#define Test_SIMULINK_COMMON_INCLUDES_
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
#include "Test_SIMULINK_types.h"
#include <stddef.h>
#include "rtGetInf.h"
#include "rtw_modelmap_simtarget.h"
#include "rt_defines.h"
#include <string.h>
#define MODEL_NAME Test_SIMULINK
#define NSAMPLE_TIMES (5) 
#define NINPUTS (0)       
#define NOUTPUTS (0)     
#define NBLOCKIO (25) 
#define NUM_ZC_EVENTS (0) 
#ifndef NCSTATES
#define NCSTATES (5)   
#elif NCSTATES != 5
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
typedef struct { real_T hkl5tp0nal ; } mtm33kqezg ; typedef struct { int32_T
hqhfxzgi1o ; boolean_T l2dv0gcr51 ; } ktdbkejzpq ; typedef struct { real_T
on2zrfzugk ; real_T k0yhgn3yum ; real_T bxgnxcai3w ; real_T pueubsjsy3 ;
real_T ajnbcnsq0m ; real_T iqglqsdfl0 ; real_T e4mzvctw3d ; real_T ceikr4ok1o
; real_T mpv4x2aipu ; real_T ong5hipbvc ; real_T a2fqeu1reh ; real_T
oe0vvwikyz ; real_T lomy52t04g ; real_T erlql5weea ; real_T eepwwoikht ;
real_T oq3b1dbzbl ; real_T hk5rehgenf ; real_T f3wtreuwqi ; real_T gzpepbf35g
; real_T l25jt0va01 ; real_T jbhastp02r ; mtm33kqezg j41fhxein4 ; mtm33kqezg
oz433gll2h ; } B ; typedef struct { real_T edgpmmcdxz [ 2 ] ; real_T
g3bm2udrh3 ; real_T mirhypy3qp ; real_T lw4cj301yx ; real_T mxttvbicvy ;
real_T lfqr0kdyo1 ; real_T cybiqktxai ; uint64_T deqpvmkrvg ; struct { real_T
modelTStart ; } ambxv44yhm ; struct { void * TUbufferPtrs [ 2 ] ; }
gfeejxrons ; struct { void * LoggedData ; } llzi04jmpj ; int32_T jj2clsenbw ;
int32_T cgiaeud0qw ; int32_T puo5ij3u42 ; int32_T kvlgkps2mn ; struct { int_T
Tail ; int_T Head ; int_T Last ; int_T CircularBufSize ; int_T MaxNewBufSize
; } e2ywx2u2km ; int_T pxo55lx3k5 ; int_T c1cwvs2f2k ; boolean_T p5uek20fdm ;
boolean_T k5rlxihffg ; boolean_T pyirayv2ej ; boolean_T om2ukmvpgs ;
boolean_T mfrbedwhfy ; boolean_T pjbp3mhdwx ; boolean_T pasxll2tbi ;
boolean_T gcjvvbwzzj ; ktdbkejzpq j41fhxein4 ; ktdbkejzpq oz433gll2h ; } DW ;
typedef struct { real_T fb2c0z2l0l ; real_T fz0tkjcd02 [ 3 ] ; real_T
cfwi3jz4pa ; } X ; typedef struct { real_T fb2c0z2l0l ; real_T fz0tkjcd02 [ 3
] ; real_T cfwi3jz4pa ; } XDot ; typedef struct { boolean_T fb2c0z2l0l ;
boolean_T fz0tkjcd02 [ 3 ] ; boolean_T cfwi3jz4pa ; } XDis ; typedef struct {
real_T fb2c0z2l0l ; real_T fz0tkjcd02 [ 3 ] ; real_T cfwi3jz4pa ; }
CStateAbsTol ; typedef struct { real_T fb2c0z2l0l ; real_T fz0tkjcd02 [ 3 ] ;
real_T cfwi3jz4pa ; } CXPtMin ; typedef struct { real_T fb2c0z2l0l ; real_T
fz0tkjcd02 [ 3 ] ; real_T cfwi3jz4pa ; } CXPtMax ; typedef struct { real_T
bosj54yo2g ; real_T c5wh512a2k ; } ZCV ; typedef struct {
rtwCAPI_ModelMappingInfo mmi ; } DataMapInfo ; struct P_ { real_T TEC ;
real_T PWM_Period ; real_T DerivativeKickDF1_X0 ; real_T DerivativeKickDF2_X0
; real_T DerivativeKickDF_X0 ; real_T
DiscreteTransferFcnwithinitialstates2_X0 ; real_T
DiscreteTransferFcnwithinitialstates1_X0 ; real_T Consigner_Time ; real_T
Consigner_Y0 ; real_T Consigner_YFinal ; real_T Dynamique2_A ; real_T
Dynamique2_C ; real_T gain_in_Gain ; real_T _Gain ; real_T _Gain_lwe4ykerbv ;
real_T Bruitdemesure_Amplitude ; real_T Bruitdemesure_Frequency ; real_T
DiscreteStateSpace_A [ 2 ] ; real_T DiscreteStateSpace_C [ 2 ] ; real_T
DiscreteStateSpace_D ; real_T Poids1_Gain ; real_T Retard2_Delay ; real_T
Retard2_InitOutput ; real_T gain_in_Gain_noqrl1m1xd ; real_T _Gain_gnigm4iwxo
; real_T _Gain_lsczmwd1h2 ; real_T _Gain_j1zbr4lzxq ; real_T
Bruitdemesure_Amplitude_hbiskpfbn4 ; real_T
Bruitdemesure_Frequency_jbbh4d4sn1 ; real_T DiscreteStateSpace_A_mopcmmoov2 ;
real_T DiscreteStateSpace_C_j3dzy0w2di ; real_T
DiscreteStateSpace_D_fbqt1av0ft ; real_T Poids2_Gain ; real_T
FiltretripleRC_A [ 3 ] ; real_T FiltretripleRC_C [ 3 ] ; real_T s_Gain ;
real_T DiscreteStateSpace_A_asx1t5v32r ; real_T
DiscreteStateSpace_C_jykma24tmg ; real_T DiscreteStateSpace_D_iep130ejoo ;
real_T Dynamique1_A ; real_T Dynamique1_C ; real_T u_op_Time ; real_T u_op_Y0
; real_T u_op_YFinal ; real_T DiscreteStateSpace_D_mvsskepuks ; real_T
DiscreteStateSpace_A_jfbtcvit52 ; real_T DiscreteStateSpace_C_j3s1yy5c2n ;
real_T Saturation_UpperSat ; real_T Saturation_LowerSat ; real_T TEC2_Gain ;
real_T TEC1_Gain ; real_T baseline_temp4_Value ; real_T baseline_temp1_Value
; real_T baseline_temp2_Value ; real_T offset_Value ; real_T offset1_Value ;
real_T offset_Value_jz2p0lxylj ; real_T offset1_Value_nyr0ia51bl ; real_T
baseline_temp1_Value_lk15cqnntl ; real_T baseline_temp2_Value_jauwwxrr5h ;
real_T baseline_temp3_Value ; real_T baseline_temp_Value ; } ; extern const
char_T * RT_MEMORY_ALLOCATION_ERROR ; extern B rtB ; extern X rtX ; extern DW
rtDW ; extern P rtP ; extern mxArray * mr_Test_SIMULINK_GetDWork ( ) ; extern
void mr_Test_SIMULINK_SetDWork ( const mxArray * ssDW ) ; extern mxArray *
mr_Test_SIMULINK_GetSimStateDisallowedBlocks ( ) ; extern const
rtwCAPI_ModelMappingStaticInfo * Test_SIMULINK_GetCAPIStaticMap ( void ) ;
extern SimStruct * const rtS ; extern DataMapInfo * rt_dataMapInfoPtr ;
extern rtwCAPI_ModelMappingInfo * rt_modelMapInfoPtr ; void MdlOutputs ( int_T
tid ) ; void MdlOutputsParameterSampleTime ( int_T tid ) ; void MdlUpdate ( int_T tid ) ; void MdlTerminate ( void ) ; void MdlInitializeSizes ( void ) ; void MdlInitializeSampleTimes ( void ) ; SimStruct * raccel_register_model ( ssExecutionInfo * executionInfo ) ;
#endif

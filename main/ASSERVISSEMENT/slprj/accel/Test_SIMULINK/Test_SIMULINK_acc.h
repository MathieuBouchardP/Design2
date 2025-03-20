#ifndef Test_SIMULINK_acc_h_
#define Test_SIMULINK_acc_h_
#ifndef Test_SIMULINK_acc_COMMON_INCLUDES_
#define Test_SIMULINK_acc_COMMON_INCLUDES_
#include <stdlib.h>
#define S_FUNCTION_NAME simulink_only_sfcn
#define S_FUNCTION_LEVEL 2
#ifndef RTW_GENERATED_S_FUNCTION
#define RTW_GENERATED_S_FUNCTION
#endif
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_nonfinite.h"
#include "math.h"
#endif
#include "Test_SIMULINK_acc_types.h"
#include <stddef.h>
#include <string.h>
#include "rtGetInf.h"
#include "rt_defines.h"
typedef struct { real_T B_6_0_0 ; real_T B_6_1_8 ; real_T B_6_2_16 ; real_T
B_6_3_24 ; real_T B_6_4_32 ; real_T B_6_5_40 ; real_T B_6_6_48 ; real_T
B_6_7_56 ; real_T B_6_8_64 ; real_T B_6_9_72 ; real_T B_6_10_80 ; real_T
B_6_11_88 ; real_T B_6_12_96 ; real_T B_6_13_104 ; real_T B_6_14_112 ; real_T
B_6_15_120 ; real_T B_6_16_128 ; real_T B_6_17_136 ; real_T B_6_18_144 ;
real_T B_6_19_152 ; real_T B_6_20_160 ; real_T B_6_21_168 ; real_T B_6_22_176
; real_T B_6_23_184 ; real_T B_6_24_192 ; real_T B_6_25_200 ; real_T
B_6_26_208 ; real_T B_6_27_216 ; real_T B_6_28_224 ; real_T B_6_29_232 ;
real_T B_6_30_240 ; real_T B_6_31_248 ; real_T B_6_32_256 ; real_T B_6_33_264
; real_T B_6_34_272 ; real_T B_5_35_280 ; real_T B_4_36_288 ; real_T
B_3_37_296 ; real_T B_2_38_304 ; real_T B_1_39_312 ; real_T B_0_40_320 ; }
B_Test_SIMULINK_T ; typedef struct { real_T DiscreteStateSpace_DSTATE [ 2 ] ;
real_T DiscreteStateSpace_DSTATE_i ; real_T DiscreteStateSpace_DSTATE_j ;
real_T DiscreteStateSpace_DSTATE_c ; real_T startTimeOfNextCycle ; real_T
timeWhenPeriodLastChanged ; real_T previousPeriod ; uint64_T
nCyclesWithSamePeriod ; struct { real_T modelTStart ; } Retard2_RWORK ;
struct { void * TUbufferPtrs [ 2 ] ; } Retard2_PWORK ; void *
Consigneretsortiey4_PWORK ; int32_T Thermistance_sysIdxToRun ; int32_T
Converttotemp_sysIdxToRun ; int32_T ConvRsistancetension_sysIdxToRun ;
int32_T Thermistance_sysIdxToRun_i ; int32_T Converttotemp_sysIdxToRun_n ;
int32_T ConvRsistancetension_sysIdxToRun_a ; struct { int_T Tail ; int_T Head
; int_T Last ; int_T CircularBufSize ; int_T MaxNewBufSize ; } Retard2_IWORK
; int_T Consigner_MODE ; int_T u_op_MODE ; boolean_T nextOutput ; boolean_T
isStartOfNextCycle ; boolean_T isFirstWarningDCGreaterThanOne ; boolean_T
isFirstWarningDCLessThanZero ; } DW_Test_SIMULINK_T ; typedef struct { real_T
Dynamique2_CSTATE ; real_T FiltretripleRC_CSTATE [ 3 ] ; real_T
Dynamique1_CSTATE ; } X_Test_SIMULINK_T ; typedef struct { real_T
Dynamique2_CSTATE ; real_T FiltretripleRC_CSTATE [ 3 ] ; real_T
Dynamique1_CSTATE ; } XDot_Test_SIMULINK_T ; typedef struct { boolean_T
Dynamique2_CSTATE ; boolean_T FiltretripleRC_CSTATE [ 3 ] ; boolean_T
Dynamique1_CSTATE ; } XDis_Test_SIMULINK_T ; typedef struct { real_T
Dynamique2_CSTATE ; real_T FiltretripleRC_CSTATE [ 3 ] ; real_T
Dynamique1_CSTATE ; } CStateAbsTol_Test_SIMULINK_T ; typedef struct { real_T
Dynamique2_CSTATE ; real_T FiltretripleRC_CSTATE [ 3 ] ; real_T
Dynamique1_CSTATE ; } CXPtMin_Test_SIMULINK_T ; typedef struct { real_T
Dynamique2_CSTATE ; real_T FiltretripleRC_CSTATE [ 3 ] ; real_T
Dynamique1_CSTATE ; } CXPtMax_Test_SIMULINK_T ; typedef struct { real_T
Consigner_StepTime_ZC ; real_T u_op_StepTime_ZC ; } ZCV_Test_SIMULINK_T ;
typedef struct { ZCSigState Consigner_StepTime_ZCE ; ZCSigState
u_op_StepTime_ZCE ; } PrevZCX_Test_SIMULINK_T ; struct P_Test_SIMULINK_T_ {
real_T P_0 ; real_T P_1 ; real_T P_2 ; real_T P_3 ; real_T P_4 ; real_T P_5 ;
real_T P_6 ; real_T P_7 ; real_T P_8 ; real_T P_9 ; real_T P_10 [ 3 ] ;
real_T P_11 ; real_T P_12 [ 2 ] ; real_T P_13 ; real_T P_14 ; real_T P_15 ;
real_T P_16 ; real_T P_17 ; real_T P_18 ; real_T P_19 ; real_T P_20 ; real_T
P_21 ; real_T P_22 ; real_T P_23 ; real_T P_24 ; real_T P_25 ; real_T P_26 ;
real_T P_27 ; real_T P_28 ; real_T P_29 ; real_T P_30 [ 3 ] ; real_T P_31 [ 3
] ; real_T P_32 ; real_T P_33 ; real_T P_34 ; real_T P_35 ; real_T P_36 ;
real_T P_37 ; real_T P_38 ; real_T P_39 ; real_T P_40 ; real_T P_41 ; real_T
P_42 ; real_T P_46 ; real_T P_47 ; real_T P_48 ; real_T P_49 ; real_T P_50 ;
real_T P_52 ; real_T P_53 ; real_T P_54 ; real_T P_55 ; real_T P_56 ; real_T
P_57 ; real_T P_58 ; real_T P_59 ; real_T P_60 ; real_T P_61 ; real_T P_62 ;
real_T P_63 ; real_T P_64 ; real_T P_65 ; real_T P_66 ; real_T P_67 ; real_T
P_68 ; real_T P_69 ; } ; extern P_Test_SIMULINK_T Test_SIMULINK_rtDefaultP ;
#endif

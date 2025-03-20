#include "rtw_capi.h"
#ifdef HOST_CAPI_BUILD
#include "MAIN_SIMULINK_capi_host.h"
#define sizeof(...) ((size_t)(0xFFFF))
#undef rt_offsetof
#define rt_offsetof(s,el) ((uint16_T)(0xFFFF))
#define TARGET_CONST
#define TARGET_STRING(s) (s)
#ifndef SS_UINT64
#define SS_UINT64 17
#endif
#ifndef SS_INT64
#define SS_INT64 18
#endif
#else
#include "builtin_typeid_types.h"
#include "MAIN_SIMULINK.h"
#include "MAIN_SIMULINK_capi.h"
#include "MAIN_SIMULINK_private.h"
#ifdef LIGHT_WEIGHT_CAPI
#define TARGET_CONST
#define TARGET_STRING(s)               ((NULL))
#else
#define TARGET_CONST                   const
#define TARGET_STRING(s)               (s)
#endif
#endif
static const rtwCAPI_Signals rtBlockSignals [ ] = { { 0 , 0 , TARGET_STRING ( "MAIN_SIMULINK/TEC" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 1 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Sum" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 2 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Sum1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 3 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Sum2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 4 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Manual Switch1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 5 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 6 , 1 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Conv. Résistance à tension" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 7 , 2 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Convert to temp" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 8 , 3 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Thermistance" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 9 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Sum3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 10 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Sum8" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 11 , 4 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Conv. Résistance à tension" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 12 , 5 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Convert to temp" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 13 , 6 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Thermistance" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 14 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Sum3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 15 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Sum8" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 16 , 0 , TARGET_STRING ( "MAIN_SIMULINK/PWM/TEC1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 3 } , { 17 , 0 , TARGET_STRING ( "MAIN_SIMULINK/PWM/TEC2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 18 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Saturation" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 19 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Sum4" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 20 , 0 , TARGET_STRING ( "MAIN_SIMULINK/T1/Dynamique2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 21 , 0 , TARGET_STRING ( "MAIN_SIMULINK/T2/Dynamique1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 22 , 0 , TARGET_STRING ( "MAIN_SIMULINK/T2/Retard2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 23 , 0 , TARGET_STRING ( "MAIN_SIMULINK/T3/Dynamique3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 24 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/ADC/Zero-Order Hold" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 25 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Conditionnement 2/gain_in" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 26 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/ADC/Zero-Order Hold" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 27 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Conditionnement 2/gain_in" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 28 , 0 , TARGET_STRING ( "MAIN_SIMULINK/PWM/PWM/Variable Pulse Generator" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 3 } , { 29 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)1/Discrete State Space" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 30 , 0 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)2/Discrete State Space" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 } } ; static const rtwCAPI_BlockParameters rtBlockParameters [ ] = { { 31 , TARGET_STRING ( "MAIN_SIMULINK/Derivative Kick (DF)" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 32 , TARGET_STRING ( "MAIN_SIMULINK/baseline_temp" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 33 , TARGET_STRING ( "MAIN_SIMULINK/Poids 1" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 34 , TARGET_STRING ( "MAIN_SIMULINK/Poids 2" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 35 , TARGET_STRING ( "MAIN_SIMULINK/Consigne r" ) , TARGET_STRING ( "Time" ) , 0 , 0 , 0 } , { 36 , TARGET_STRING ( "MAIN_SIMULINK/Consigne r" ) , TARGET_STRING ( "Before" ) , 0 , 0 , 0 } , { 37 , TARGET_STRING ( "MAIN_SIMULINK/Consigne r" ) , TARGET_STRING ( "After" ) , 0 , 0 , 0 } , { 38 , TARGET_STRING ( "MAIN_SIMULINK/Manual Switch" ) , TARGET_STRING ( "CurrentSetting" ) , 1 , 0 , 0 } , { 39 , TARGET_STRING ( "MAIN_SIMULINK/Manual Switch1" ) , TARGET_STRING ( "CurrentSetting" ) , 1 , 0 , 0 } , { 40 , TARGET_STRING ( "MAIN_SIMULINK/Conditionnement 1/baseline_temp4" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 41 , TARGET_STRING ( "MAIN_SIMULINK/Conditionnement 1/s" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 42 , TARGET_STRING ( "MAIN_SIMULINK/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 43 , TARGET_STRING ( "MAIN_SIMULINK/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 44 , TARGET_STRING ( "MAIN_SIMULINK/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "D" ) , 0 , 0 , 0 } , { 45 , TARGET_STRING ( "MAIN_SIMULINK/Filtre trible RC/Filtre triple RC" ) , TARGET_STRING ( "A" ) , 0 , 1 , 0 } , { 46 , TARGET_STRING ( "MAIN_SIMULINK/Filtre trible RC/Filtre triple RC" ) , TARGET_STRING ( "C" ) , 0 , 2 , 0 } , { 47 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/baseline_temp1" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 48 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/baseline_temp2" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 49 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Bruit de mesure" ) , TARGET_STRING ( "Amplitude" ) , 0 , 0 , 0 } , { 50 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Bruit de mesure" ) , TARGET_STRING ( "Frequency" ) , 0 , 0 , 0 } , { 51 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/baseline_temp1" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 52 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/baseline_temp2" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 53 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Bruit de mesure" ) , TARGET_STRING ( "Amplitude" ) , 0 , 0 , 0 } , { 54 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Bruit de mesure" ) , TARGET_STRING ( "Frequency" ) , 0 , 0 , 0 } , { 55 , TARGET_STRING ( "MAIN_SIMULINK/PWM/PWM" ) , TARGET_STRING ( "Period" ) , 0 , 0 , 0 } , { 56 , TARGET_STRING ( "MAIN_SIMULINK/PWM/baseline_temp3" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 57 , TARGET_STRING ( "MAIN_SIMULINK/PWM/TEC1" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 58 , TARGET_STRING ( "MAIN_SIMULINK/PWM/TEC2" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 59 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)1" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 60 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)2" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 61 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Saturation" ) , TARGET_STRING ( "UpperLimit" ) , 0 , 0 , 0 } , { 62 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Saturation" ) , TARGET_STRING ( "LowerLimit" ) , 0 , 0 , 0 } , { 63 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/u_op" ) , TARGET_STRING ( "Time" ) , 0 , 0 , 0 } , { 64 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/u_op" ) , TARGET_STRING ( "Before" ) , 0 , 0 , 0 } , { 65 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/u_op" ) , TARGET_STRING ( "After" ) , 0 , 0 , 0 } , { 66 , TARGET_STRING ( "MAIN_SIMULINK/T1/Dynamique2" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 67 , TARGET_STRING ( "MAIN_SIMULINK/T1/Dynamique2" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 68 , TARGET_STRING ( "MAIN_SIMULINK/T1->T3/Derivative Kick (DF)" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 69 , TARGET_STRING ( "MAIN_SIMULINK/T2/Dynamique1" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 70 , TARGET_STRING ( "MAIN_SIMULINK/T2/Dynamique1" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 71 , TARGET_STRING ( "MAIN_SIMULINK/T2/Retard2" ) , TARGET_STRING ( "DelayTime" ) , 0 , 0 , 0 } , { 72 , TARGET_STRING ( "MAIN_SIMULINK/T2/Retard2" ) , TARGET_STRING ( "InitialOutput" ) , 0 , 0 , 0 } , { 73 , TARGET_STRING ( "MAIN_SIMULINK/T2->T3/Derivative Kick (DF)" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 74 , TARGET_STRING ( "MAIN_SIMULINK/T3/Dynamique3" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 75 , TARGET_STRING ( "MAIN_SIMULINK/T3/Dynamique3" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 76 , TARGET_STRING ( "MAIN_SIMULINK/T3/Retard3" ) , TARGET_STRING ( "DelayTime" ) , 0 , 0 , 0 } , { 77 , TARGET_STRING ( "MAIN_SIMULINK/T3/Retard3" ) , TARGET_STRING ( "InitialOutput" ) , 0 , 0 , 0 } , { 78 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/ADC/ " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 79 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Conditionnement 2/offset" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 80 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Conditionnement 2/gain_in" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 81 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Conditionnement 3/offset1" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 82 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Conditionnement 3/  " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 83 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T1/Conditionnement 3/   " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 84 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/ADC/ " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 85 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Conditionnement 2/offset" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 86 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Conditionnement 2/gain_in" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 87 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Conditionnement 3/offset1" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 88 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Conditionnement 3/  " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 89 , TARGET_STRING ( "MAIN_SIMULINK/Obtenir température T2/Conditionnement 3/   " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 90 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)1/Discrete State Space" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 91 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)1/Discrete State Space" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 92 , TARGET_STRING ( "MAIN_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)2/Discrete State Space" ) , TARGET_STRING ( "D" ) , 0 , 0 , 0 } , { 93 , TARGET_STRING ( "MAIN_SIMULINK/T1->T3/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "A" ) , 0 , 3 , 0 } , { 94 , TARGET_STRING ( "MAIN_SIMULINK/T1->T3/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "C" ) , 0 , 4 , 0 } , { 95 , TARGET_STRING ( "MAIN_SIMULINK/T1->T3/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "D" ) , 0 , 0 , 0 } , { 96 , TARGET_STRING ( "MAIN_SIMULINK/T2->T3/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 97 , TARGET_STRING ( "MAIN_SIMULINK/T2->T3/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 98 , TARGET_STRING ( "MAIN_SIMULINK/T2->T3/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "D" ) , 0 , 0 , 0 } , { 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 } } ; static int_T rt_LoggedStateIdxList [ ] = { - 1 } ; static const rtwCAPI_Signals rtRootInputs [ ] = { { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 } } ; static const rtwCAPI_Signals rtRootOutputs [ ] = { { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 } } ; static const rtwCAPI_ModelParameters rtModelParameters [ ] = { { 99 , TARGET_STRING ( "TEC" ) , 0 , 0 , 0 } , { 0 , ( NULL ) , 0 , 0 , 0 } } ;
#ifndef HOST_CAPI_BUILD
static void * rtDataAddrMap [ ] = { & rtB . ay1xeqbo54 , & rtB . jl4tzlqkvr ,
& rtB . lma1ydcnif , & rtB . fsu23y5ts3 , & rtB . db4m1bc3gg , & rtB .
olqwozfrl3 , & rtB . d355u2baya . bu3amqfjr4 , & rtB . nlg3idhhmk .
p0it2epuk4 , & rtB . cdh5o32feg . ordo20z2um , & rtB . h13eabb5l1 , & rtB .
oobv21my4e , & rtB . arhszp1yvs . bu3amqfjr4 , & rtB . n55rhng3cw .
p0it2epuk4 , & rtB . mfc24avgl1 . ordo20z2um , & rtB . mzld5yp0dz , & rtB .
fqgdvg40sk , & rtB . mxyvw51mkl , & rtB . hxdb1eoutr , & rtB . jccqxfucp1 , &
rtB . hkgktvl34f , & rtB . ccwcs4orvs , & rtB . lbhnidlzjv , & rtB .
k1g1lcxbvk , & rtB . p2k4zdl445 , & rtB . jfqdwkeqgu , & rtB . exsilx1jdl , &
rtB . eelznhvfta , & rtB . hssf4p0lhx , & rtB . op004jo4lt , & rtB .
eblcjzhc4e , & rtB . lijzxwwv3g , & rtP . DerivativeKickDF_X0_pcnom2kezu , &
rtP . baseline_temp_Value , & rtP . Poids1_Gain , & rtP . Poids2_Gain , & rtP
. Consigner_Time , & rtP . Consigner_Y0 , & rtP . Consigner_YFinal , & rtP .
ManualSwitch_CurrentSetting , & rtP . ManualSwitch1_CurrentSetting , & rtP .
baseline_temp4_Value , & rtP . s_Gain , & rtP .
DiscreteStateSpace_A_asx1t5v32r , & rtP . DiscreteStateSpace_C_jykma24tmg , &
rtP . DiscreteStateSpace_D_iep130ejoo , & rtP . FiltretripleRC_A [ 0 ] , &
rtP . FiltretripleRC_C [ 0 ] , & rtP . baseline_temp1_Value , & rtP .
baseline_temp2_Value , & rtP . Bruitdemesure_Amplitude , & rtP .
Bruitdemesure_Frequency , & rtP . baseline_temp1_Value_lk15cqnntl , & rtP .
baseline_temp2_Value_jauwwxrr5h , & rtP . Bruitdemesure_Amplitude_hbiskpfbn4
, & rtP . Bruitdemesure_Frequency_jbbh4d4sn1 , & rtP . PWM_Period , & rtP .
baseline_temp3_Value , & rtP . TEC1_Gain , & rtP . TEC2_Gain , & rtP .
DiscreteTransferFcnwithinitialstates1_X0 , & rtP .
DiscreteTransferFcnwithinitialstates2_X0 , & rtP . Saturation_UpperSat , &
rtP . Saturation_LowerSat , & rtP . u_op_Time , & rtP . u_op_Y0 , & rtP .
u_op_YFinal , & rtP . Dynamique2_A , & rtP . Dynamique2_C , & rtP .
DerivativeKickDF_X0 , & rtP . Dynamique1_A , & rtP . Dynamique1_C , & rtP .
Retard2_Delay , & rtP . Retard2_InitOutput , & rtP .
DerivativeKickDF_X0_pde1xmny1t , & rtP . Dynamique3_A , & rtP . Dynamique3_C
, & rtP . Retard3_Delay , & rtP . Retard3_InitOutput , & rtP . _Gain , & rtP
. offset_Value , & rtP . gain_in_Gain , & rtP . offset1_Value , & rtP .
_Gain_j0borcplw2 , & rtP . _Gain_lwe4ykerbv , & rtP . _Gain_gnigm4iwxo , &
rtP . offset_Value_jz2p0lxylj , & rtP . gain_in_Gain_noqrl1m1xd , & rtP .
offset1_Value_nyr0ia51bl , & rtP . _Gain_lsczmwd1h2 , & rtP .
_Gain_j1zbr4lzxq , & rtP . DiscreteStateSpace_A_jfbtcvit52 , & rtP .
DiscreteStateSpace_C_j3s1yy5c2n , & rtP . DiscreteStateSpace_D_mvsskepuks , &
rtP . DiscreteStateSpace_A [ 0 ] , & rtP . DiscreteStateSpace_C [ 0 ] , & rtP
. DiscreteStateSpace_D , & rtP . DiscreteStateSpace_A_mopcmmoov2 , & rtP .
DiscreteStateSpace_C_j3dzy0w2di , & rtP . DiscreteStateSpace_D_fbqt1av0ft , &
rtP . TEC , } ; static int32_T * rtVarDimsAddrMap [ ] = { ( NULL ) } ;
#endif
static TARGET_CONST rtwCAPI_DataTypeMap rtDataTypeMap [ ] = { { "double" ,
"real_T" , 0 , 0 , sizeof ( real_T ) , ( uint8_T ) SS_DOUBLE , 0 , 0 , 0 } ,
{ "unsigned char" , "uint8_T" , 0 , 0 , sizeof ( uint8_T ) , ( uint8_T )
SS_UINT8 , 0 , 0 , 0 } } ;
#ifdef HOST_CAPI_BUILD
#undef sizeof
#endif
static TARGET_CONST rtwCAPI_ElementMap rtElementMap [ ] = { { ( NULL ) , 0 ,
0 , 0 , 0 } , } ; static const rtwCAPI_DimensionMap rtDimensionMap [ ] = { {
rtwCAPI_SCALAR , 0 , 2 , 0 } , { rtwCAPI_VECTOR , 2 , 2 , 0 } , {
rtwCAPI_VECTOR , 4 , 2 , 0 } , { rtwCAPI_VECTOR , 6 , 2 , 0 } , {
rtwCAPI_VECTOR , 8 , 2 , 0 } } ; static const uint_T rtDimensionArray [ ] = {
1 , 1 , 3 , 1 , 1 , 3 , 2 , 1 , 1 , 2 } ; static const real_T
rtcapiStoredFloats [ ] = { 0.0 , 1.0 , 5.0 , 0.1 } ; static const
rtwCAPI_FixPtMap rtFixPtMap [ ] = { { ( NULL ) , ( NULL ) ,
rtwCAPI_FIX_RESERVED , 0 , 0 , ( boolean_T ) 0 } , } ; static const
rtwCAPI_SampleTimeMap rtSampleTimeMap [ ] = { { ( const void * ) &
rtcapiStoredFloats [ 0 ] , ( const void * ) & rtcapiStoredFloats [ 0 ] , ( int8_T ) 0 , ( uint8_T ) 0 } , { ( const void * ) & rtcapiStoredFloats [ 0 ] , ( const void * ) & rtcapiStoredFloats [ 1 ] , ( int8_T ) 1 , ( uint8_T ) 0 } , { ( const void * ) & rtcapiStoredFloats [ 2 ] , ( const void * ) & rtcapiStoredFloats [ 0 ] , ( int8_T ) 3 , ( uint8_T ) 0 } , { ( const void * ) & rtcapiStoredFloats [ 3 ] , ( const void * ) & rtcapiStoredFloats [ 0 ] , ( int8_T ) 2 , ( uint8_T ) 0 } } ; static rtwCAPI_ModelMappingStaticInfo mmiStatic = { { rtBlockSignals , 31 , rtRootInputs , 0 , rtRootOutputs , 0 } , { rtBlockParameters , 68 , rtModelParameters , 1 } , { ( NULL ) , 0 } , { rtDataTypeMap , rtDimensionMap , rtFixPtMap , rtElementMap , rtSampleTimeMap , rtDimensionArray } , "float" , { 839251671U , 2922706996U , 573706206U , 3100317870U } , ( NULL ) , 0 , ( boolean_T ) 0 , rt_LoggedStateIdxList } ; const rtwCAPI_ModelMappingStaticInfo * MAIN_SIMULINK_GetCAPIStaticMap ( void ) { return & mmiStatic ; }
#ifndef HOST_CAPI_BUILD
void MAIN_SIMULINK_InitializeDataMapInfo ( void ) { rtwCAPI_SetVersion ( ( *
rt_dataMapInfoPtr ) . mmi , 1 ) ; rtwCAPI_SetStaticMap ( ( *
rt_dataMapInfoPtr ) . mmi , & mmiStatic ) ; rtwCAPI_SetLoggingStaticMap ( ( *
rt_dataMapInfoPtr ) . mmi , ( NULL ) ) ; rtwCAPI_SetDataAddressMap ( ( *
rt_dataMapInfoPtr ) . mmi , rtDataAddrMap ) ; rtwCAPI_SetVarDimsAddressMap ( ( *
rt_dataMapInfoPtr ) . mmi , rtVarDimsAddrMap ) ;
rtwCAPI_SetInstanceLoggingInfo ( ( * rt_dataMapInfoPtr ) . mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArray ( ( * rt_dataMapInfoPtr ) . mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArrayLen ( ( * rt_dataMapInfoPtr ) . mmi , 0 ) ; }
#else
#ifdef __cplusplus
extern "C" {
#endif
void MAIN_SIMULINK_host_InitializeDataMapInfo ( MAIN_SIMULINK_host_DataMapInfo_T
* dataMap , const char * path ) { rtwCAPI_SetVersion ( dataMap -> mmi , 1 ) ;
rtwCAPI_SetStaticMap ( dataMap -> mmi , & mmiStatic ) ;
rtwCAPI_SetDataAddressMap ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetVarDimsAddressMap ( dataMap -> mmi , ( NULL ) ) ; rtwCAPI_SetPath
( dataMap -> mmi , path ) ; rtwCAPI_SetFullPath ( dataMap -> mmi , ( NULL ) )
; rtwCAPI_SetChildMMIArray ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArrayLen ( dataMap -> mmi , 0 ) ; }
#ifdef __cplusplus
}
#endif
#endif

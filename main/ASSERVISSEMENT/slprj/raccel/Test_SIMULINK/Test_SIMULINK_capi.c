#include "rtw_capi.h"
#ifdef HOST_CAPI_BUILD
#include "Test_SIMULINK_capi_host.h"
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
#include "Test_SIMULINK.h"
#include "Test_SIMULINK_capi.h"
#include "Test_SIMULINK_private.h"
#ifdef LIGHT_WEIGHT_CAPI
#define TARGET_CONST
#define TARGET_STRING(s)               ((NULL))
#else
#define TARGET_CONST                   const
#define TARGET_STRING(s)               (s)
#endif
#endif
static const rtwCAPI_Signals rtBlockSignals [ ] = { { 0 , 0 , TARGET_STRING ( "Test_SIMULINK/TEC" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 1 , 0 , TARGET_STRING ( "Test_SIMULINK/Sum" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 2 , 0 , TARGET_STRING ( "Test_SIMULINK/Sum1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 3 , 0 , TARGET_STRING ( "Test_SIMULINK/Sum2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 4 , 0 , TARGET_STRING ( "Test_SIMULINK/Dynamique1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 5 , 0 , TARGET_STRING ( "Test_SIMULINK/Dynamique2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 6 , 0 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 7 , 1 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/Conv. Résistance à tension" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 8 , 3 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/Thermistance" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 9 , 0 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/gain_in" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 10 , 0 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/Sum3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 11 , 0 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/Sum8" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 12 , 4 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Conv. Résistance à tension" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 13 , 6 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Thermistance" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 14 , 0 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Sum3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 15 , 0 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Sum8" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 16 , 0 , TARGET_STRING ( "Test_SIMULINK/PWM/TEC1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 17 , 0 , TARGET_STRING ( "Test_SIMULINK/PWM/TEC2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 18 , 0 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Sum4" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 19 , 0 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Conditionnement 2/gain_in" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 20 , 0 , TARGET_STRING ( "Test_SIMULINK/PWM/PWM/Variable Pulse Generator" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 21 , 0 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)1/Discrete State Space" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 22 , 0 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)2/Discrete State Space" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 2 } , { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 } } ; static const rtwCAPI_BlockParameters rtBlockParameters [ ] = { { 23 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 24 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)1" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 25 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)2" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 26 , TARGET_STRING ( "Test_SIMULINK/baseline_temp" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 27 , TARGET_STRING ( "Test_SIMULINK/Poids 1" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 28 , TARGET_STRING ( "Test_SIMULINK/Poids 2" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 29 , TARGET_STRING ( "Test_SIMULINK/Consigne r" ) , TARGET_STRING ( "Time" ) , 0 , 0 , 0 } , { 30 , TARGET_STRING ( "Test_SIMULINK/Consigne r" ) , TARGET_STRING ( "Before" ) , 0 , 0 , 0 } , { 31 , TARGET_STRING ( "Test_SIMULINK/Consigne r" ) , TARGET_STRING ( "After" ) , 0 , 0 , 0 } , { 32 , TARGET_STRING ( "Test_SIMULINK/Dynamique1" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 33 , TARGET_STRING ( "Test_SIMULINK/Dynamique1" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 34 , TARGET_STRING ( "Test_SIMULINK/Dynamique2" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 35 , TARGET_STRING ( "Test_SIMULINK/Dynamique2" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 36 , TARGET_STRING ( "Test_SIMULINK/Retard2" ) , TARGET_STRING ( "DelayTime" ) , 0 , 0 , 0 } , { 37 , TARGET_STRING ( "Test_SIMULINK/Retard2" ) , TARGET_STRING ( "InitialOutput" ) , 0 , 0 , 0 } , { 38 , TARGET_STRING ( "Test_SIMULINK/Conditionnement 1/baseline_temp4" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 39 , TARGET_STRING ( "Test_SIMULINK/Conditionnement 1/s" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 40 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 41 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 42 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)/Discrete State Space" ) , TARGET_STRING ( "D" ) , 0 , 0 , 0 } , { 43 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)1/Discrete State Space" ) , TARGET_STRING ( "A" ) , 0 , 1 , 0 } , { 44 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)1/Discrete State Space" ) , TARGET_STRING ( "C" ) , 0 , 2 , 0 } , { 45 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)1/Discrete State Space" ) , TARGET_STRING ( "D" ) , 0 , 0 , 0 } , { 46 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)2/Discrete State Space" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 47 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)2/Discrete State Space" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 48 , TARGET_STRING ( "Test_SIMULINK/Derivative Kick (DF)2/Discrete State Space" ) , TARGET_STRING ( "D" ) , 0 , 0 , 0 } , { 49 , TARGET_STRING ( "Test_SIMULINK/Filtre trible RC/Filtre triple RC" ) , TARGET_STRING ( "A" ) , 0 , 3 , 0 } , { 50 , TARGET_STRING ( "Test_SIMULINK/Filtre trible RC/Filtre triple RC" ) , TARGET_STRING ( "C" ) , 0 , 4 , 0 } , { 51 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/baseline_temp1" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 52 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/baseline_temp2" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 53 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/offset" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 54 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/offset1" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 55 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/   " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 56 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/gain_in" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 57 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/Bruit de mesure" ) , TARGET_STRING ( "Amplitude" ) , 0 , 0 , 0 } , { 58 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/Bruit de mesure" ) , TARGET_STRING ( "Frequency" ) , 0 , 0 , 0 } , { 59 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/baseline_temp1" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 60 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/baseline_temp2" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 61 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Bruit de mesure" ) , TARGET_STRING ( "Amplitude" ) , 0 , 0 , 0 } , { 62 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Bruit de mesure" ) , TARGET_STRING ( "Frequency" ) , 0 , 0 , 0 } , { 63 , TARGET_STRING ( "Test_SIMULINK/PWM/PWM" ) , TARGET_STRING ( "Period" ) , 0 , 0 , 0 } , { 64 , TARGET_STRING ( "Test_SIMULINK/PWM/baseline_temp3" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 65 , TARGET_STRING ( "Test_SIMULINK/PWM/TEC1" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 66 , TARGET_STRING ( "Test_SIMULINK/PWM/TEC2" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 67 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)1" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 68 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)2" ) , TARGET_STRING ( "X0" ) , 0 , 0 , 0 } , { 69 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Saturation" ) , TARGET_STRING ( "UpperLimit" ) , 0 , 0 , 0 } , { 70 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Saturation" ) , TARGET_STRING ( "LowerLimit" ) , 0 , 0 , 0 } , { 71 , TARGET_STRING ( "Test_SIMULINK/Régulateur/u_op" ) , TARGET_STRING ( "Time" ) , 0 , 0 , 0 } , { 72 , TARGET_STRING ( "Test_SIMULINK/Régulateur/u_op" ) , TARGET_STRING ( "Before" ) , 0 , 0 , 0 } , { 73 , TARGET_STRING ( "Test_SIMULINK/Régulateur/u_op" ) , TARGET_STRING ( "After" ) , 0 , 0 , 0 } , { 74 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T1/ADC/ " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 75 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/ADC/ " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 76 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Conditionnement 2/offset" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 77 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Conditionnement 2/gain_in" ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 78 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Conditionnement 3/offset1" ) , TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 79 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Conditionnement 3/  " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 80 , TARGET_STRING ( "Test_SIMULINK/Obtenir température T2/Conditionnement 3/   " ) , TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 81 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)1/Discrete State Space" ) , TARGET_STRING ( "A" ) , 0 , 0 , 0 } , { 82 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)1/Discrete State Space" ) , TARGET_STRING ( "C" ) , 0 , 0 , 0 } , { 83 , TARGET_STRING ( "Test_SIMULINK/Régulateur/Discrete Transfer Fcn (with initial states)2/Discrete State Space" ) , TARGET_STRING ( "D" ) , 0 , 0 , 0 } , { 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 } } ; static int_T rt_LoggedStateIdxList [ ] = { - 1 } ; static const rtwCAPI_Signals rtRootInputs [ ] = { { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 } } ; static const rtwCAPI_Signals rtRootOutputs [ ] = { { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 } } ; static const rtwCAPI_ModelParameters rtModelParameters [ ] = { { 84 , TARGET_STRING ( "TEC" ) , 0 , 0 , 0 } , { 0 , ( NULL ) , 0 , 0 , 0 } } ;
#ifndef HOST_CAPI_BUILD
static void * rtDataAddrMap [ ] = { & rtB . gzpepbf35g , & rtB . oe0vvwikyz ,
& rtB . mpv4x2aipu , & rtB . on2zrfzugk , & rtB . a2fqeu1reh , & rtB .
k0yhgn3yum , & rtB . ong5hipbvc , & rtB . oz433gll2h . hkl5tp0nal , & rtB .
jbhastp02r , & rtB . bxgnxcai3w , & rtB . ajnbcnsq0m , & rtB . pueubsjsy3 , &
rtB . j41fhxein4 . hkl5tp0nal , & rtB . l25jt0va01 , & rtB . ceikr4ok1o , &
rtB . e4mzvctw3d , & rtB . hk5rehgenf , & rtB . eepwwoikht , & rtB .
f3wtreuwqi , & rtB . iqglqsdfl0 , & rtB . oq3b1dbzbl , & rtB . erlql5weea , &
rtB . lomy52t04g , & rtP . DerivativeKickDF_X0 , & rtP . DerivativeKickDF1_X0
, & rtP . DerivativeKickDF2_X0 , & rtP . baseline_temp_Value , & rtP .
Poids1_Gain , & rtP . Poids2_Gain , & rtP . Consigner_Time , & rtP .
Consigner_Y0 , & rtP . Consigner_YFinal , & rtP . Dynamique1_A , & rtP .
Dynamique1_C , & rtP . Dynamique2_A , & rtP . Dynamique2_C , & rtP .
Retard2_Delay , & rtP . Retard2_InitOutput , & rtP . baseline_temp4_Value , &
rtP . s_Gain , & rtP . DiscreteStateSpace_A_asx1t5v32r , & rtP .
DiscreteStateSpace_C_jykma24tmg , & rtP . DiscreteStateSpace_D_iep130ejoo , &
rtP . DiscreteStateSpace_A [ 0 ] , & rtP . DiscreteStateSpace_C [ 0 ] , & rtP
. DiscreteStateSpace_D , & rtP . DiscreteStateSpace_A_mopcmmoov2 , & rtP .
DiscreteStateSpace_C_j3dzy0w2di , & rtP . DiscreteStateSpace_D_fbqt1av0ft , &
rtP . FiltretripleRC_A [ 0 ] , & rtP . FiltretripleRC_C [ 0 ] , & rtP .
baseline_temp1_Value , & rtP . baseline_temp2_Value , & rtP . offset_Value ,
& rtP . offset1_Value , & rtP . _Gain_lwe4ykerbv , & rtP . gain_in_Gain , &
rtP . Bruitdemesure_Amplitude , & rtP . Bruitdemesure_Frequency , & rtP .
baseline_temp1_Value_lk15cqnntl , & rtP . baseline_temp2_Value_jauwwxrr5h , &
rtP . Bruitdemesure_Amplitude_hbiskpfbn4 , & rtP .
Bruitdemesure_Frequency_jbbh4d4sn1 , & rtP . PWM_Period , & rtP .
baseline_temp3_Value , & rtP . TEC1_Gain , & rtP . TEC2_Gain , & rtP .
DiscreteTransferFcnwithinitialstates1_X0 , & rtP .
DiscreteTransferFcnwithinitialstates2_X0 , & rtP . Saturation_UpperSat , &
rtP . Saturation_LowerSat , & rtP . u_op_Time , & rtP . u_op_Y0 , & rtP .
u_op_YFinal , & rtP . _Gain , & rtP . _Gain_gnigm4iwxo , & rtP .
offset_Value_jz2p0lxylj , & rtP . gain_in_Gain_noqrl1m1xd , & rtP .
offset1_Value_nyr0ia51bl , & rtP . _Gain_lsczmwd1h2 , & rtP .
_Gain_j1zbr4lzxq , & rtP . DiscreteStateSpace_A_jfbtcvit52 , & rtP .
DiscreteStateSpace_C_j3s1yy5c2n , & rtP . DiscreteStateSpace_D_mvsskepuks , &
rtP . TEC , } ; static int32_T * rtVarDimsAddrMap [ ] = { ( NULL ) } ;
#endif
static TARGET_CONST rtwCAPI_DataTypeMap rtDataTypeMap [ ] = { { "double" ,
"real_T" , 0 , 0 , sizeof ( real_T ) , ( uint8_T ) SS_DOUBLE , 0 , 0 , 0 } }
;
#ifdef HOST_CAPI_BUILD
#undef sizeof
#endif
static TARGET_CONST rtwCAPI_ElementMap rtElementMap [ ] = { { ( NULL ) , 0 ,
0 , 0 , 0 } , } ; static const rtwCAPI_DimensionMap rtDimensionMap [ ] = { {
rtwCAPI_SCALAR , 0 , 2 , 0 } , { rtwCAPI_VECTOR , 2 , 2 , 0 } , {
rtwCAPI_VECTOR , 4 , 2 , 0 } , { rtwCAPI_VECTOR , 6 , 2 , 0 } , {
rtwCAPI_VECTOR , 8 , 2 , 0 } } ; static const uint_T rtDimensionArray [ ] = {
1 , 1 , 2 , 1 , 1 , 2 , 3 , 1 , 1 , 3 } ; static const real_T
rtcapiStoredFloats [ ] = { 0.0 , 1.0 , 5.0 } ; static const rtwCAPI_FixPtMap
rtFixPtMap [ ] = { { ( NULL ) , ( NULL ) , rtwCAPI_FIX_RESERVED , 0 , 0 , ( boolean_T ) 0 } , } ; static const rtwCAPI_SampleTimeMap rtSampleTimeMap [ ] = { { ( const void * ) & rtcapiStoredFloats [ 0 ] , ( const void * ) & rtcapiStoredFloats [ 0 ] , ( int8_T ) 0 , ( uint8_T ) 0 } , { ( const void * ) & rtcapiStoredFloats [ 0 ] , ( const void * ) & rtcapiStoredFloats [ 1 ] , ( int8_T ) 1 , ( uint8_T ) 0 } , { ( const void * ) & rtcapiStoredFloats [ 2 ] , ( const void * ) & rtcapiStoredFloats [ 0 ] , ( int8_T ) 2 , ( uint8_T ) 0 } } ; static rtwCAPI_ModelMappingStaticInfo mmiStatic = { { rtBlockSignals , 23 , rtRootInputs , 0 , rtRootOutputs , 0 } , { rtBlockParameters , 61 , rtModelParameters , 1 } , { ( NULL ) , 0 } , { rtDataTypeMap , rtDimensionMap , rtFixPtMap , rtElementMap , rtSampleTimeMap , rtDimensionArray } , "float" , { 979203763U , 2237083764U , 3775564064U , 168080820U } , ( NULL ) , 0 , ( boolean_T ) 0 , rt_LoggedStateIdxList } ; const rtwCAPI_ModelMappingStaticInfo * Test_SIMULINK_GetCAPIStaticMap ( void ) { return & mmiStatic ; }
#ifndef HOST_CAPI_BUILD
void Test_SIMULINK_InitializeDataMapInfo ( void ) { rtwCAPI_SetVersion ( ( *
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
void Test_SIMULINK_host_InitializeDataMapInfo ( Test_SIMULINK_host_DataMapInfo_T
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

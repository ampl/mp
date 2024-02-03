#ifndef NLW2_BASIC_DEFS_C_H
#define NLW2_BASIC_DEFS_C_H

#ifdef __cplusplus
extern "C" {
#endif

const int NLW2_ObjSenseMinimize = 0;
const int NLW2_ObjSenseMaximize = 1;

const int NLW2_VarTypeContinuous = 0;
const int NLW2_VarTypeInteger = 1;

//const int NLW2_VarTypeSemiContinuous = 2;
//const int NLW2_VarTypeSemiInteger = 3;
//const int NLW2_VarTypeImplicitInteger = 4;

const int NLW2_MatrixFormatColwise = 1;
const int NLW2_MatrixFormatRowwise = 2;

const int NLW2_HessianFormatTriangular = 1;
const int NLW2_HessianFormatSquare = 2;

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLW2_BASIC_DEFS_C_H

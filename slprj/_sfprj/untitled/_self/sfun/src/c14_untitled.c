/* Include files */

#include <stddef.h>
#include "blas.h"
#include "untitled_sfun.h"
#include "c14_untitled.h"
#include <math.h>
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "untitled_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c14_debug_family_names[4] = { "nargin", "nargout", "A",
  "exponent_A" };

/* Function Declarations */
static void initialize_c14_untitled(SFc14_untitledInstanceStruct *chartInstance);
static void initialize_params_c14_untitled(SFc14_untitledInstanceStruct
  *chartInstance);
static void enable_c14_untitled(SFc14_untitledInstanceStruct *chartInstance);
static void disable_c14_untitled(SFc14_untitledInstanceStruct *chartInstance);
static void c14_update_debugger_state_c14_untitled(SFc14_untitledInstanceStruct *
  chartInstance);
static const mxArray *get_sim_state_c14_untitled(SFc14_untitledInstanceStruct
  *chartInstance);
static void set_sim_state_c14_untitled(SFc14_untitledInstanceStruct
  *chartInstance, const mxArray *c14_st);
static void finalize_c14_untitled(SFc14_untitledInstanceStruct *chartInstance);
static void sf_gateway_c14_untitled(SFc14_untitledInstanceStruct *chartInstance);
static void initSimStructsc14_untitled(SFc14_untitledInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c14_machineNumber, uint32_T
  c14_chartNumber, uint32_T c14_instanceNumber);
static const mxArray *c14_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static void c14_emlrt_marshallIn(SFc14_untitledInstanceStruct *chartInstance,
  const mxArray *c14_exponent_A, const char_T *c14_identifier, real_T c14_y[16]);
static void c14_b_emlrt_marshallIn(SFc14_untitledInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId, real_T c14_y[16]);
static void c14_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static const mxArray *c14_b_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static real_T c14_c_emlrt_marshallIn(SFc14_untitledInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId);
static void c14_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static void c14_info_helper(const mxArray **c14_info);
static const mxArray *c14_emlrt_marshallOut(const char * c14_u);
static const mxArray *c14_b_emlrt_marshallOut(const uint32_T c14_u);
static void c14_b_info_helper(const mxArray **c14_info);
static void c14_c_info_helper(const mxArray **c14_info);
static void c14_expm(SFc14_untitledInstanceStruct *chartInstance, real_T c14_A
                     [16], real_T c14_F[16]);
static void c14_PadeApproximantOfDegree(SFc14_untitledInstanceStruct
  *chartInstance, real_T c14_A[16], real_T c14_m, real_T c14_F[16]);
static void c14_eml_scalar_eg(SFc14_untitledInstanceStruct *chartInstance);
static void c14_eml_xgemm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16], real_T c14_C[16], real_T c14_b_C[16]);
static void c14_eml_switch_helper(SFc14_untitledInstanceStruct *chartInstance);
static void c14_eml_matlab_zgetrf(SFc14_untitledInstanceStruct *chartInstance,
  real_T c14_A[16], real_T c14_b_A[16], int32_T c14_ipiv[4], int32_T *c14_info);
static int32_T c14_eml_ixamax(SFc14_untitledInstanceStruct *chartInstance,
  int32_T c14_n, real_T c14_x[16], int32_T c14_ix0);
static void c14_check_forloop_overflow_error(SFc14_untitledInstanceStruct
  *chartInstance, boolean_T c14_overflow);
static void c14_threshold(SFc14_untitledInstanceStruct *chartInstance);
static void c14_eml_xgeru(SFc14_untitledInstanceStruct *chartInstance, int32_T
  c14_m, int32_T c14_n, real_T c14_alpha1, int32_T c14_ix0, int32_T c14_iy0,
  real_T c14_A[16], int32_T c14_ia0, real_T c14_b_A[16]);
static void c14_b_eml_scalar_eg(SFc14_untitledInstanceStruct *chartInstance);
static void c14_eml_warning(SFc14_untitledInstanceStruct *chartInstance);
static void c14_eml_xtrsm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16], real_T c14_b_B[16]);
static void c14_b_threshold(SFc14_untitledInstanceStruct *chartInstance);
static void c14_scalarEg(SFc14_untitledInstanceStruct *chartInstance);
static void c14_b_eml_xtrsm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16], real_T c14_b_B[16]);
static const mxArray *c14_c_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static int32_T c14_d_emlrt_marshallIn(SFc14_untitledInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId);
static void c14_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static uint8_T c14_e_emlrt_marshallIn(SFc14_untitledInstanceStruct
  *chartInstance, const mxArray *c14_b_is_active_c14_untitled, const char_T
  *c14_identifier);
static uint8_T c14_f_emlrt_marshallIn(SFc14_untitledInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId);
static void c14_b_eml_xgemm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16], real_T c14_C[16]);
static void c14_b_eml_matlab_zgetrf(SFc14_untitledInstanceStruct *chartInstance,
  real_T c14_A[16], int32_T c14_ipiv[4], int32_T *c14_info);
static void c14_b_eml_xgeru(SFc14_untitledInstanceStruct *chartInstance, int32_T
  c14_m, int32_T c14_n, real_T c14_alpha1, int32_T c14_ix0, int32_T c14_iy0,
  real_T c14_A[16], int32_T c14_ia0);
static void c14_c_eml_xtrsm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16]);
static void c14_d_eml_xtrsm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16]);
static void init_dsm_address_info(SFc14_untitledInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c14_untitled(SFc14_untitledInstanceStruct *chartInstance)
{
  chartInstance->c14_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c14_is_active_c14_untitled = 0U;
}

static void initialize_params_c14_untitled(SFc14_untitledInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c14_untitled(SFc14_untitledInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c14_untitled(SFc14_untitledInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c14_update_debugger_state_c14_untitled(SFc14_untitledInstanceStruct *
  chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c14_untitled(SFc14_untitledInstanceStruct
  *chartInstance)
{
  const mxArray *c14_st;
  const mxArray *c14_y = NULL;
  int32_T c14_i0;
  real_T c14_u[16];
  const mxArray *c14_b_y = NULL;
  uint8_T c14_hoistedGlobal;
  uint8_T c14_b_u;
  const mxArray *c14_c_y = NULL;
  real_T (*c14_exponent_A)[16];
  c14_exponent_A = (real_T (*)[16])ssGetOutputPortSignal(chartInstance->S, 1);
  c14_st = NULL;
  c14_st = NULL;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_createcellmatrix(2, 1), false);
  for (c14_i0 = 0; c14_i0 < 16; c14_i0++) {
    c14_u[c14_i0] = (*c14_exponent_A)[c14_i0];
  }

  c14_b_y = NULL;
  sf_mex_assign(&c14_b_y, sf_mex_create("y", c14_u, 0, 0U, 1U, 0U, 2, 4, 4),
                false);
  sf_mex_setcell(c14_y, 0, c14_b_y);
  c14_hoistedGlobal = chartInstance->c14_is_active_c14_untitled;
  c14_b_u = c14_hoistedGlobal;
  c14_c_y = NULL;
  sf_mex_assign(&c14_c_y, sf_mex_create("y", &c14_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c14_y, 1, c14_c_y);
  sf_mex_assign(&c14_st, c14_y, false);
  return c14_st;
}

static void set_sim_state_c14_untitled(SFc14_untitledInstanceStruct
  *chartInstance, const mxArray *c14_st)
{
  const mxArray *c14_u;
  real_T c14_dv0[16];
  int32_T c14_i1;
  real_T (*c14_exponent_A)[16];
  c14_exponent_A = (real_T (*)[16])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c14_doneDoubleBufferReInit = true;
  c14_u = sf_mex_dup(c14_st);
  c14_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c14_u, 0)),
                       "exponent_A", c14_dv0);
  for (c14_i1 = 0; c14_i1 < 16; c14_i1++) {
    (*c14_exponent_A)[c14_i1] = c14_dv0[c14_i1];
  }

  chartInstance->c14_is_active_c14_untitled = c14_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c14_u, 1)),
     "is_active_c14_untitled");
  sf_mex_destroy(&c14_u);
  c14_update_debugger_state_c14_untitled(chartInstance);
  sf_mex_destroy(&c14_st);
}

static void finalize_c14_untitled(SFc14_untitledInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c14_untitled(SFc14_untitledInstanceStruct *chartInstance)
{
  int32_T c14_i2;
  int32_T c14_i3;
  real_T c14_A[16];
  uint32_T c14_debug_family_var_map[4];
  real_T c14_nargin = 1.0;
  real_T c14_nargout = 1.0;
  real_T c14_exponent_A[16];
  int32_T c14_i4;
  real_T c14_b_A[16];
  real_T c14_dv1[16];
  int32_T c14_i5;
  int32_T c14_i6;
  int32_T c14_i7;
  real_T (*c14_b_exponent_A)[16];
  real_T (*c14_c_A)[16];
  c14_b_exponent_A = (real_T (*)[16])ssGetOutputPortSignal(chartInstance->S, 1);
  c14_c_A = (real_T (*)[16])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c14_sfEvent);
  for (c14_i2 = 0; c14_i2 < 16; c14_i2++) {
    _SFD_DATA_RANGE_CHECK((*c14_c_A)[c14_i2], 0U);
  }

  chartInstance->c14_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c14_sfEvent);
  for (c14_i3 = 0; c14_i3 < 16; c14_i3++) {
    c14_A[c14_i3] = (*c14_c_A)[c14_i3];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c14_debug_family_names,
    c14_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_nargin, 0U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_nargout, 1U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c14_A, 2U, c14_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c14_exponent_A, 3U, c14_sf_marshallOut,
    c14_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 3);
  for (c14_i4 = 0; c14_i4 < 16; c14_i4++) {
    c14_b_A[c14_i4] = c14_A[c14_i4];
  }

  c14_expm(chartInstance, c14_b_A, c14_dv1);
  for (c14_i5 = 0; c14_i5 < 16; c14_i5++) {
    c14_exponent_A[c14_i5] = c14_dv1[c14_i5];
  }

  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, -3);
  _SFD_SYMBOL_SCOPE_POP();
  for (c14_i6 = 0; c14_i6 < 16; c14_i6++) {
    (*c14_b_exponent_A)[c14_i6] = c14_exponent_A[c14_i6];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c14_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_untitledMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c14_i7 = 0; c14_i7 < 16; c14_i7++) {
    _SFD_DATA_RANGE_CHECK((*c14_b_exponent_A)[c14_i7], 1U);
  }
}

static void initSimStructsc14_untitled(SFc14_untitledInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c14_machineNumber, uint32_T
  c14_chartNumber, uint32_T c14_instanceNumber)
{
  (void)c14_machineNumber;
  (void)c14_chartNumber;
  (void)c14_instanceNumber;
}

static const mxArray *c14_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  int32_T c14_i8;
  int32_T c14_i9;
  int32_T c14_i10;
  real_T c14_b_inData[16];
  int32_T c14_i11;
  int32_T c14_i12;
  int32_T c14_i13;
  real_T c14_u[16];
  const mxArray *c14_y = NULL;
  SFc14_untitledInstanceStruct *chartInstance;
  chartInstance = (SFc14_untitledInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  c14_i8 = 0;
  for (c14_i9 = 0; c14_i9 < 4; c14_i9++) {
    for (c14_i10 = 0; c14_i10 < 4; c14_i10++) {
      c14_b_inData[c14_i10 + c14_i8] = (*(real_T (*)[16])c14_inData)[c14_i10 +
        c14_i8];
    }

    c14_i8 += 4;
  }

  c14_i11 = 0;
  for (c14_i12 = 0; c14_i12 < 4; c14_i12++) {
    for (c14_i13 = 0; c14_i13 < 4; c14_i13++) {
      c14_u[c14_i13 + c14_i11] = c14_b_inData[c14_i13 + c14_i11];
    }

    c14_i11 += 4;
  }

  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", c14_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, false);
  return c14_mxArrayOutData;
}

static void c14_emlrt_marshallIn(SFc14_untitledInstanceStruct *chartInstance,
  const mxArray *c14_exponent_A, const char_T *c14_identifier, real_T c14_y[16])
{
  emlrtMsgIdentifier c14_thisId;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_exponent_A), &c14_thisId,
    c14_y);
  sf_mex_destroy(&c14_exponent_A);
}

static void c14_b_emlrt_marshallIn(SFc14_untitledInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId, real_T c14_y[16])
{
  real_T c14_dv2[16];
  int32_T c14_i14;
  (void)chartInstance;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), c14_dv2, 1, 0, 0U, 1, 0U, 2, 4,
                4);
  for (c14_i14 = 0; c14_i14 < 16; c14_i14++) {
    c14_y[c14_i14] = c14_dv2[c14_i14];
  }

  sf_mex_destroy(&c14_u);
}

static void c14_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_exponent_A;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  real_T c14_y[16];
  int32_T c14_i15;
  int32_T c14_i16;
  int32_T c14_i17;
  SFc14_untitledInstanceStruct *chartInstance;
  chartInstance = (SFc14_untitledInstanceStruct *)chartInstanceVoid;
  c14_exponent_A = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_exponent_A), &c14_thisId,
    c14_y);
  sf_mex_destroy(&c14_exponent_A);
  c14_i15 = 0;
  for (c14_i16 = 0; c14_i16 < 4; c14_i16++) {
    for (c14_i17 = 0; c14_i17 < 4; c14_i17++) {
      (*(real_T (*)[16])c14_outData)[c14_i17 + c14_i15] = c14_y[c14_i17 +
        c14_i15];
    }

    c14_i15 += 4;
  }

  sf_mex_destroy(&c14_mxArrayInData);
}

static const mxArray *c14_b_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  real_T c14_u;
  const mxArray *c14_y = NULL;
  SFc14_untitledInstanceStruct *chartInstance;
  chartInstance = (SFc14_untitledInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  c14_u = *(real_T *)c14_inData;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", &c14_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, false);
  return c14_mxArrayOutData;
}

static real_T c14_c_emlrt_marshallIn(SFc14_untitledInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId)
{
  real_T c14_y;
  real_T c14_d0;
  (void)chartInstance;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), &c14_d0, 1, 0, 0U, 0, 0U, 0);
  c14_y = c14_d0;
  sf_mex_destroy(&c14_u);
  return c14_y;
}

static void c14_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_nargout;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  real_T c14_y;
  SFc14_untitledInstanceStruct *chartInstance;
  chartInstance = (SFc14_untitledInstanceStruct *)chartInstanceVoid;
  c14_nargout = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_y = c14_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_nargout),
    &c14_thisId);
  sf_mex_destroy(&c14_nargout);
  *(real_T *)c14_outData = c14_y;
  sf_mex_destroy(&c14_mxArrayInData);
}

const mxArray *sf_c14_untitled_get_eml_resolved_functions_info(void)
{
  const mxArray *c14_nameCaptureInfo = NULL;
  c14_nameCaptureInfo = NULL;
  sf_mex_assign(&c14_nameCaptureInfo, sf_mex_createstruct("structure", 2, 164, 1),
                false);
  c14_info_helper(&c14_nameCaptureInfo);
  c14_b_info_helper(&c14_nameCaptureInfo);
  c14_c_info_helper(&c14_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c14_nameCaptureInfo);
  return c14_nameCaptureInfo;
}

static void c14_info_helper(const mxArray **c14_info)
{
  const mxArray *c14_rhs0 = NULL;
  const mxArray *c14_lhs0 = NULL;
  const mxArray *c14_rhs1 = NULL;
  const mxArray *c14_lhs1 = NULL;
  const mxArray *c14_rhs2 = NULL;
  const mxArray *c14_lhs2 = NULL;
  const mxArray *c14_rhs3 = NULL;
  const mxArray *c14_lhs3 = NULL;
  const mxArray *c14_rhs4 = NULL;
  const mxArray *c14_lhs4 = NULL;
  const mxArray *c14_rhs5 = NULL;
  const mxArray *c14_lhs5 = NULL;
  const mxArray *c14_rhs6 = NULL;
  const mxArray *c14_lhs6 = NULL;
  const mxArray *c14_rhs7 = NULL;
  const mxArray *c14_lhs7 = NULL;
  const mxArray *c14_rhs8 = NULL;
  const mxArray *c14_lhs8 = NULL;
  const mxArray *c14_rhs9 = NULL;
  const mxArray *c14_lhs9 = NULL;
  const mxArray *c14_rhs10 = NULL;
  const mxArray *c14_lhs10 = NULL;
  const mxArray *c14_rhs11 = NULL;
  const mxArray *c14_lhs11 = NULL;
  const mxArray *c14_rhs12 = NULL;
  const mxArray *c14_lhs12 = NULL;
  const mxArray *c14_rhs13 = NULL;
  const mxArray *c14_lhs13 = NULL;
  const mxArray *c14_rhs14 = NULL;
  const mxArray *c14_lhs14 = NULL;
  const mxArray *c14_rhs15 = NULL;
  const mxArray *c14_lhs15 = NULL;
  const mxArray *c14_rhs16 = NULL;
  const mxArray *c14_lhs16 = NULL;
  const mxArray *c14_rhs17 = NULL;
  const mxArray *c14_lhs17 = NULL;
  const mxArray *c14_rhs18 = NULL;
  const mxArray *c14_lhs18 = NULL;
  const mxArray *c14_rhs19 = NULL;
  const mxArray *c14_lhs19 = NULL;
  const mxArray *c14_rhs20 = NULL;
  const mxArray *c14_lhs20 = NULL;
  const mxArray *c14_rhs21 = NULL;
  const mxArray *c14_lhs21 = NULL;
  const mxArray *c14_rhs22 = NULL;
  const mxArray *c14_lhs22 = NULL;
  const mxArray *c14_rhs23 = NULL;
  const mxArray *c14_lhs23 = NULL;
  const mxArray *c14_rhs24 = NULL;
  const mxArray *c14_lhs24 = NULL;
  const mxArray *c14_rhs25 = NULL;
  const mxArray *c14_lhs25 = NULL;
  const mxArray *c14_rhs26 = NULL;
  const mxArray *c14_lhs26 = NULL;
  const mxArray *c14_rhs27 = NULL;
  const mxArray *c14_lhs27 = NULL;
  const mxArray *c14_rhs28 = NULL;
  const mxArray *c14_lhs28 = NULL;
  const mxArray *c14_rhs29 = NULL;
  const mxArray *c14_lhs29 = NULL;
  const mxArray *c14_rhs30 = NULL;
  const mxArray *c14_lhs30 = NULL;
  const mxArray *c14_rhs31 = NULL;
  const mxArray *c14_lhs31 = NULL;
  const mxArray *c14_rhs32 = NULL;
  const mxArray *c14_lhs32 = NULL;
  const mxArray *c14_rhs33 = NULL;
  const mxArray *c14_lhs33 = NULL;
  const mxArray *c14_rhs34 = NULL;
  const mxArray *c14_lhs34 = NULL;
  const mxArray *c14_rhs35 = NULL;
  const mxArray *c14_lhs35 = NULL;
  const mxArray *c14_rhs36 = NULL;
  const mxArray *c14_lhs36 = NULL;
  const mxArray *c14_rhs37 = NULL;
  const mxArray *c14_lhs37 = NULL;
  const mxArray *c14_rhs38 = NULL;
  const mxArray *c14_lhs38 = NULL;
  const mxArray *c14_rhs39 = NULL;
  const mxArray *c14_lhs39 = NULL;
  const mxArray *c14_rhs40 = NULL;
  const mxArray *c14_lhs40 = NULL;
  const mxArray *c14_rhs41 = NULL;
  const mxArray *c14_lhs41 = NULL;
  const mxArray *c14_rhs42 = NULL;
  const mxArray *c14_lhs42 = NULL;
  const mxArray *c14_rhs43 = NULL;
  const mxArray *c14_lhs43 = NULL;
  const mxArray *c14_rhs44 = NULL;
  const mxArray *c14_lhs44 = NULL;
  const mxArray *c14_rhs45 = NULL;
  const mxArray *c14_lhs45 = NULL;
  const mxArray *c14_rhs46 = NULL;
  const mxArray *c14_lhs46 = NULL;
  const mxArray *c14_rhs47 = NULL;
  const mxArray *c14_lhs47 = NULL;
  const mxArray *c14_rhs48 = NULL;
  const mxArray *c14_lhs48 = NULL;
  const mxArray *c14_rhs49 = NULL;
  const mxArray *c14_lhs49 = NULL;
  const mxArray *c14_rhs50 = NULL;
  const mxArray *c14_lhs50 = NULL;
  const mxArray *c14_rhs51 = NULL;
  const mxArray *c14_lhs51 = NULL;
  const mxArray *c14_rhs52 = NULL;
  const mxArray *c14_lhs52 = NULL;
  const mxArray *c14_rhs53 = NULL;
  const mxArray *c14_lhs53 = NULL;
  const mxArray *c14_rhs54 = NULL;
  const mxArray *c14_lhs54 = NULL;
  const mxArray *c14_rhs55 = NULL;
  const mxArray *c14_lhs55 = NULL;
  const mxArray *c14_rhs56 = NULL;
  const mxArray *c14_lhs56 = NULL;
  const mxArray *c14_rhs57 = NULL;
  const mxArray *c14_lhs57 = NULL;
  const mxArray *c14_rhs58 = NULL;
  const mxArray *c14_lhs58 = NULL;
  const mxArray *c14_rhs59 = NULL;
  const mxArray *c14_lhs59 = NULL;
  const mxArray *c14_rhs60 = NULL;
  const mxArray *c14_lhs60 = NULL;
  const mxArray *c14_rhs61 = NULL;
  const mxArray *c14_lhs61 = NULL;
  const mxArray *c14_rhs62 = NULL;
  const mxArray *c14_lhs62 = NULL;
  const mxArray *c14_rhs63 = NULL;
  const mxArray *c14_lhs63 = NULL;
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("expm"), "name", "name", 0);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1381857504U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c14_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 1);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c14_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("ismatrix"), "name", "name",
                  2);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c14_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("norm"), "name", "name", 3);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c14_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c14_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 5);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("abs"), "name", "name", 5);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c14_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 6);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c14_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 7);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c14_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 8);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("isnan"), "name", "name", 8);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c14_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 9);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c14_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 10);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 10);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c14_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "context", "context", 11);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 11);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c14_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m!PadeApproximantOfDegree"),
                  "context", "context", 12);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 12);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c14_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 13);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 13);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c14_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 14);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 14);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c14_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 15);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 15);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c14_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 16);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c14_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  17);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c14_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 18);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c14_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 19);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c14_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 20);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 20);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c14_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 21);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 21);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c14_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 22);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 22);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c14_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 23);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 23);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c14_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 24);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.refblas.xgemm"), "name", "name", 24);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c14_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m!PadeApproximantOfDegree"),
                  "context", "context", 25);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("mldivide"), "name", "name",
                  25);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1319737166U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c14_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "context",
                  "context", 26);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_lusolve"), "name",
                  "name", 26);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1370017086U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c14_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 27);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_xgetrf"), "name", "name",
                  27);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286826006U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c14_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "context", "context", 28);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_lapack_xgetrf"), "name",
                  "name", 28);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286826010U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c14_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "context", "context", 29);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_matlab_zgetrf"), "name",
                  "name", 29);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1302696194U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c14_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("realmin"), "name", "name",
                  30);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 30);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c14_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_realmin"), "name",
                  "name", 31);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 31);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1307658444U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c14_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 32);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c14_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 33);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eps"), "name", "name", 33);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c14_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 34);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 34);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c14_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_eps"), "name", "name",
                  35);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 35);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c14_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 36);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 36);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c14_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 37);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("min"), "name", "name", 37);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 37);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 37);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1311262518U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c14_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 38);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1378303184U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c14_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 39);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 39);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 39);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c14_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 40);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 40);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c14_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 41);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 41);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 41);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c14_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 42);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 42);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 42);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c14_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 43);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 43);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c14_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 44);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 44);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 44);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c14_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 45);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 45);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 45);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c14_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 46);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("colon"), "name", "name", 46);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1378303188U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c14_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 47);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("colon"), "name", "name", 47);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 47);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1378303188U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c14_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 48);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 48);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c14_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 49);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 49);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 49);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c14_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 50);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("floor"), "name", "name", 50);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 50);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c14_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 51);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 51);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c14_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 52);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 52);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c14_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 53);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("intmin"), "name", "name", 53);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 53);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c14_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 54);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 54);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c14_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 55);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("intmax"), "name", "name", 55);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 55);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c14_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 56);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 56);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c14_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 57);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("intmin"), "name", "name", 57);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 57);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c14_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 58);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("intmax"), "name", "name", 58);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 58);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c14_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 59);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_isa_uint"), "name",
                  "name", 59);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 59);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 59);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c14_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "context",
                  "context", 60);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.isaUint"),
                  "name", "name", 60);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 60);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/isaUint.p"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c14_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 61);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_unsigned_class"), "name",
                  "name", 61);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c14_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "context", "context", 62);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.unsignedClass"), "name", "name", 62);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c14_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "context", "context", 63);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 63);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c14_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c14_rhs0);
  sf_mex_destroy(&c14_lhs0);
  sf_mex_destroy(&c14_rhs1);
  sf_mex_destroy(&c14_lhs1);
  sf_mex_destroy(&c14_rhs2);
  sf_mex_destroy(&c14_lhs2);
  sf_mex_destroy(&c14_rhs3);
  sf_mex_destroy(&c14_lhs3);
  sf_mex_destroy(&c14_rhs4);
  sf_mex_destroy(&c14_lhs4);
  sf_mex_destroy(&c14_rhs5);
  sf_mex_destroy(&c14_lhs5);
  sf_mex_destroy(&c14_rhs6);
  sf_mex_destroy(&c14_lhs6);
  sf_mex_destroy(&c14_rhs7);
  sf_mex_destroy(&c14_lhs7);
  sf_mex_destroy(&c14_rhs8);
  sf_mex_destroy(&c14_lhs8);
  sf_mex_destroy(&c14_rhs9);
  sf_mex_destroy(&c14_lhs9);
  sf_mex_destroy(&c14_rhs10);
  sf_mex_destroy(&c14_lhs10);
  sf_mex_destroy(&c14_rhs11);
  sf_mex_destroy(&c14_lhs11);
  sf_mex_destroy(&c14_rhs12);
  sf_mex_destroy(&c14_lhs12);
  sf_mex_destroy(&c14_rhs13);
  sf_mex_destroy(&c14_lhs13);
  sf_mex_destroy(&c14_rhs14);
  sf_mex_destroy(&c14_lhs14);
  sf_mex_destroy(&c14_rhs15);
  sf_mex_destroy(&c14_lhs15);
  sf_mex_destroy(&c14_rhs16);
  sf_mex_destroy(&c14_lhs16);
  sf_mex_destroy(&c14_rhs17);
  sf_mex_destroy(&c14_lhs17);
  sf_mex_destroy(&c14_rhs18);
  sf_mex_destroy(&c14_lhs18);
  sf_mex_destroy(&c14_rhs19);
  sf_mex_destroy(&c14_lhs19);
  sf_mex_destroy(&c14_rhs20);
  sf_mex_destroy(&c14_lhs20);
  sf_mex_destroy(&c14_rhs21);
  sf_mex_destroy(&c14_lhs21);
  sf_mex_destroy(&c14_rhs22);
  sf_mex_destroy(&c14_lhs22);
  sf_mex_destroy(&c14_rhs23);
  sf_mex_destroy(&c14_lhs23);
  sf_mex_destroy(&c14_rhs24);
  sf_mex_destroy(&c14_lhs24);
  sf_mex_destroy(&c14_rhs25);
  sf_mex_destroy(&c14_lhs25);
  sf_mex_destroy(&c14_rhs26);
  sf_mex_destroy(&c14_lhs26);
  sf_mex_destroy(&c14_rhs27);
  sf_mex_destroy(&c14_lhs27);
  sf_mex_destroy(&c14_rhs28);
  sf_mex_destroy(&c14_lhs28);
  sf_mex_destroy(&c14_rhs29);
  sf_mex_destroy(&c14_lhs29);
  sf_mex_destroy(&c14_rhs30);
  sf_mex_destroy(&c14_lhs30);
  sf_mex_destroy(&c14_rhs31);
  sf_mex_destroy(&c14_lhs31);
  sf_mex_destroy(&c14_rhs32);
  sf_mex_destroy(&c14_lhs32);
  sf_mex_destroy(&c14_rhs33);
  sf_mex_destroy(&c14_lhs33);
  sf_mex_destroy(&c14_rhs34);
  sf_mex_destroy(&c14_lhs34);
  sf_mex_destroy(&c14_rhs35);
  sf_mex_destroy(&c14_lhs35);
  sf_mex_destroy(&c14_rhs36);
  sf_mex_destroy(&c14_lhs36);
  sf_mex_destroy(&c14_rhs37);
  sf_mex_destroy(&c14_lhs37);
  sf_mex_destroy(&c14_rhs38);
  sf_mex_destroy(&c14_lhs38);
  sf_mex_destroy(&c14_rhs39);
  sf_mex_destroy(&c14_lhs39);
  sf_mex_destroy(&c14_rhs40);
  sf_mex_destroy(&c14_lhs40);
  sf_mex_destroy(&c14_rhs41);
  sf_mex_destroy(&c14_lhs41);
  sf_mex_destroy(&c14_rhs42);
  sf_mex_destroy(&c14_lhs42);
  sf_mex_destroy(&c14_rhs43);
  sf_mex_destroy(&c14_lhs43);
  sf_mex_destroy(&c14_rhs44);
  sf_mex_destroy(&c14_lhs44);
  sf_mex_destroy(&c14_rhs45);
  sf_mex_destroy(&c14_lhs45);
  sf_mex_destroy(&c14_rhs46);
  sf_mex_destroy(&c14_lhs46);
  sf_mex_destroy(&c14_rhs47);
  sf_mex_destroy(&c14_lhs47);
  sf_mex_destroy(&c14_rhs48);
  sf_mex_destroy(&c14_lhs48);
  sf_mex_destroy(&c14_rhs49);
  sf_mex_destroy(&c14_lhs49);
  sf_mex_destroy(&c14_rhs50);
  sf_mex_destroy(&c14_lhs50);
  sf_mex_destroy(&c14_rhs51);
  sf_mex_destroy(&c14_lhs51);
  sf_mex_destroy(&c14_rhs52);
  sf_mex_destroy(&c14_lhs52);
  sf_mex_destroy(&c14_rhs53);
  sf_mex_destroy(&c14_lhs53);
  sf_mex_destroy(&c14_rhs54);
  sf_mex_destroy(&c14_lhs54);
  sf_mex_destroy(&c14_rhs55);
  sf_mex_destroy(&c14_lhs55);
  sf_mex_destroy(&c14_rhs56);
  sf_mex_destroy(&c14_lhs56);
  sf_mex_destroy(&c14_rhs57);
  sf_mex_destroy(&c14_lhs57);
  sf_mex_destroy(&c14_rhs58);
  sf_mex_destroy(&c14_lhs58);
  sf_mex_destroy(&c14_rhs59);
  sf_mex_destroy(&c14_lhs59);
  sf_mex_destroy(&c14_rhs60);
  sf_mex_destroy(&c14_lhs60);
  sf_mex_destroy(&c14_rhs61);
  sf_mex_destroy(&c14_lhs61);
  sf_mex_destroy(&c14_rhs62);
  sf_mex_destroy(&c14_lhs62);
  sf_mex_destroy(&c14_rhs63);
  sf_mex_destroy(&c14_lhs63);
}

static const mxArray *c14_emlrt_marshallOut(const char * c14_u)
{
  const mxArray *c14_y = NULL;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", c14_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c14_u)), false);
  return c14_y;
}

static const mxArray *c14_b_emlrt_marshallOut(const uint32_T c14_u)
{
  const mxArray *c14_y = NULL;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", &c14_u, 7, 0U, 0U, 0U, 0), false);
  return c14_y;
}

static void c14_b_info_helper(const mxArray **c14_info)
{
  const mxArray *c14_rhs64 = NULL;
  const mxArray *c14_lhs64 = NULL;
  const mxArray *c14_rhs65 = NULL;
  const mxArray *c14_lhs65 = NULL;
  const mxArray *c14_rhs66 = NULL;
  const mxArray *c14_lhs66 = NULL;
  const mxArray *c14_rhs67 = NULL;
  const mxArray *c14_lhs67 = NULL;
  const mxArray *c14_rhs68 = NULL;
  const mxArray *c14_lhs68 = NULL;
  const mxArray *c14_rhs69 = NULL;
  const mxArray *c14_lhs69 = NULL;
  const mxArray *c14_rhs70 = NULL;
  const mxArray *c14_lhs70 = NULL;
  const mxArray *c14_rhs71 = NULL;
  const mxArray *c14_lhs71 = NULL;
  const mxArray *c14_rhs72 = NULL;
  const mxArray *c14_lhs72 = NULL;
  const mxArray *c14_rhs73 = NULL;
  const mxArray *c14_lhs73 = NULL;
  const mxArray *c14_rhs74 = NULL;
  const mxArray *c14_lhs74 = NULL;
  const mxArray *c14_rhs75 = NULL;
  const mxArray *c14_lhs75 = NULL;
  const mxArray *c14_rhs76 = NULL;
  const mxArray *c14_lhs76 = NULL;
  const mxArray *c14_rhs77 = NULL;
  const mxArray *c14_lhs77 = NULL;
  const mxArray *c14_rhs78 = NULL;
  const mxArray *c14_lhs78 = NULL;
  const mxArray *c14_rhs79 = NULL;
  const mxArray *c14_lhs79 = NULL;
  const mxArray *c14_rhs80 = NULL;
  const mxArray *c14_lhs80 = NULL;
  const mxArray *c14_rhs81 = NULL;
  const mxArray *c14_lhs81 = NULL;
  const mxArray *c14_rhs82 = NULL;
  const mxArray *c14_lhs82 = NULL;
  const mxArray *c14_rhs83 = NULL;
  const mxArray *c14_lhs83 = NULL;
  const mxArray *c14_rhs84 = NULL;
  const mxArray *c14_lhs84 = NULL;
  const mxArray *c14_rhs85 = NULL;
  const mxArray *c14_lhs85 = NULL;
  const mxArray *c14_rhs86 = NULL;
  const mxArray *c14_lhs86 = NULL;
  const mxArray *c14_rhs87 = NULL;
  const mxArray *c14_lhs87 = NULL;
  const mxArray *c14_rhs88 = NULL;
  const mxArray *c14_lhs88 = NULL;
  const mxArray *c14_rhs89 = NULL;
  const mxArray *c14_lhs89 = NULL;
  const mxArray *c14_rhs90 = NULL;
  const mxArray *c14_lhs90 = NULL;
  const mxArray *c14_rhs91 = NULL;
  const mxArray *c14_lhs91 = NULL;
  const mxArray *c14_rhs92 = NULL;
  const mxArray *c14_lhs92 = NULL;
  const mxArray *c14_rhs93 = NULL;
  const mxArray *c14_lhs93 = NULL;
  const mxArray *c14_rhs94 = NULL;
  const mxArray *c14_lhs94 = NULL;
  const mxArray *c14_rhs95 = NULL;
  const mxArray *c14_lhs95 = NULL;
  const mxArray *c14_rhs96 = NULL;
  const mxArray *c14_lhs96 = NULL;
  const mxArray *c14_rhs97 = NULL;
  const mxArray *c14_lhs97 = NULL;
  const mxArray *c14_rhs98 = NULL;
  const mxArray *c14_lhs98 = NULL;
  const mxArray *c14_rhs99 = NULL;
  const mxArray *c14_lhs99 = NULL;
  const mxArray *c14_rhs100 = NULL;
  const mxArray *c14_lhs100 = NULL;
  const mxArray *c14_rhs101 = NULL;
  const mxArray *c14_lhs101 = NULL;
  const mxArray *c14_rhs102 = NULL;
  const mxArray *c14_lhs102 = NULL;
  const mxArray *c14_rhs103 = NULL;
  const mxArray *c14_lhs103 = NULL;
  const mxArray *c14_rhs104 = NULL;
  const mxArray *c14_lhs104 = NULL;
  const mxArray *c14_rhs105 = NULL;
  const mxArray *c14_lhs105 = NULL;
  const mxArray *c14_rhs106 = NULL;
  const mxArray *c14_lhs106 = NULL;
  const mxArray *c14_rhs107 = NULL;
  const mxArray *c14_lhs107 = NULL;
  const mxArray *c14_rhs108 = NULL;
  const mxArray *c14_lhs108 = NULL;
  const mxArray *c14_rhs109 = NULL;
  const mxArray *c14_lhs109 = NULL;
  const mxArray *c14_rhs110 = NULL;
  const mxArray *c14_lhs110 = NULL;
  const mxArray *c14_rhs111 = NULL;
  const mxArray *c14_lhs111 = NULL;
  const mxArray *c14_rhs112 = NULL;
  const mxArray *c14_lhs112 = NULL;
  const mxArray *c14_rhs113 = NULL;
  const mxArray *c14_lhs113 = NULL;
  const mxArray *c14_rhs114 = NULL;
  const mxArray *c14_lhs114 = NULL;
  const mxArray *c14_rhs115 = NULL;
  const mxArray *c14_lhs115 = NULL;
  const mxArray *c14_rhs116 = NULL;
  const mxArray *c14_lhs116 = NULL;
  const mxArray *c14_rhs117 = NULL;
  const mxArray *c14_lhs117 = NULL;
  const mxArray *c14_rhs118 = NULL;
  const mxArray *c14_lhs118 = NULL;
  const mxArray *c14_rhs119 = NULL;
  const mxArray *c14_lhs119 = NULL;
  const mxArray *c14_rhs120 = NULL;
  const mxArray *c14_lhs120 = NULL;
  const mxArray *c14_rhs121 = NULL;
  const mxArray *c14_lhs121 = NULL;
  const mxArray *c14_rhs122 = NULL;
  const mxArray *c14_lhs122 = NULL;
  const mxArray *c14_rhs123 = NULL;
  const mxArray *c14_lhs123 = NULL;
  const mxArray *c14_rhs124 = NULL;
  const mxArray *c14_lhs124 = NULL;
  const mxArray *c14_rhs125 = NULL;
  const mxArray *c14_lhs125 = NULL;
  const mxArray *c14_rhs126 = NULL;
  const mxArray *c14_lhs126 = NULL;
  const mxArray *c14_rhs127 = NULL;
  const mxArray *c14_lhs127 = NULL;
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "context", "context", 64);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 64);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c14_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 65);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 65);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c14_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 66);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("intmax"), "name", "name", 66);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 66);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c14_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 67);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_isa_uint"), "name",
                  "name", 67);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 67);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 67);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c14_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 68);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 68);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c14_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 69);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 69);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 69);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c14_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_signed_integer_colon"),
                  "context", "context", 70);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 70);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c14_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 71);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("intmax"), "name", "name", 71);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 71);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c14_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 72);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 72);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c14_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 73);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 73);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c14_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 74);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 74);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c14_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 75);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 75);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 75);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c14_rhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 76);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 76);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 76);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c14_rhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 77);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 77);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 77);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 77);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c14_rhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 78);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 78);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c14_rhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 79);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 79);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 79);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c14_rhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 80);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 80);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 80);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 80);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c14_rhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 81);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 81);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 81);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 81);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c14_rhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 82);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 82);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 82);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 82);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c14_rhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 83);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_ixamax"), "name", "name",
                  83);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "resolved", "resolved", 83);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c14_rhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 84);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 84);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 84);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c14_rhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 85);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.ixamax"),
                  "name", "name", 85);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c14_rhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "context", "context", 86);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 86);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 86);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c14_rhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p!below_threshold"),
                  "context", "context", 87);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 87);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c14_rhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p!below_threshold"),
                  "context", "context", 88);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("length"), "name", "name", 88);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 88);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1303153406U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c14_rhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength"),
                  "context", "context", 89);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 89);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c14_rhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "context", "context", 90);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.refblas.ixamax"), "name", "name", 90);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "resolved", "resolved", 90);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c14_rhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 91);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.refblas.xcabs1"), "name", "name", 91);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c14_rhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "context", "context", 92);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("abs"), "name", "name", 92);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 92);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c14_rhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 93);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 93);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 93);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c14_rhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 94);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 94);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 94);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 94);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c14_rhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 95);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_xswap"), "name", "name",
                  95);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 95);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"),
                  "resolved", "resolved", 95);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c14_rhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 96);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 96);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 96);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 96);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c14_rhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 97);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.xswap"),
                  "name", "name", 97);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 97);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "resolved", "resolved", 97);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c14_rhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 98);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 98);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 98);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 98);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c14_rhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p!below_threshold"),
                  "context", "context", 99);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 99);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 99);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c14_rhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 100);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.refblas.xswap"), "name", "name", 100);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "resolved", "resolved", 100);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c14_rhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs100), "rhs",
                  "rhs", 100);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs100), "lhs",
                  "lhs", 100);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 101);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("abs"), "name", "name", 101);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 101);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 101);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c14_rhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs101), "rhs",
                  "rhs", 101);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs101), "lhs",
                  "lhs", 101);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 102);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 102);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 102);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 102);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c14_rhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs102), "rhs",
                  "rhs", 102);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs102), "lhs",
                  "lhs", 102);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 103);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 103);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 103);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 103);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c14_rhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs103), "rhs",
                  "rhs", 103);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs103), "lhs",
                  "lhs", 103);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 104);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 104);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 104);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c14_rhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs104), "rhs",
                  "rhs", 104);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs104), "lhs",
                  "lhs", 104);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 105);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 105);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 105);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 105);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c14_rhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs105), "rhs",
                  "rhs", 105);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs105), "lhs",
                  "lhs", 105);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 106);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_div"), "name", "name",
                  106);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 106);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 106);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c14_rhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs106), "rhs",
                  "rhs", 106);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs106), "lhs",
                  "lhs", 106);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 107);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 107);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 107);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 107);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c14_rhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs107), "rhs",
                  "rhs", 107);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs107), "lhs",
                  "lhs", 107);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 108);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_xgeru"), "name", "name",
                  108);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 108);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c14_rhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs108), "rhs",
                  "rhs", 108);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs108), "lhs",
                  "lhs", 108);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 109);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 109);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 109);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 109);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c14_rhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs109), "rhs",
                  "rhs", 109);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs109), "lhs",
                  "lhs", 109);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 110);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.xgeru"),
                  "name", "name", 110);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgeru.p"),
                  "resolved", "resolved", 110);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c14_rhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs110), "rhs",
                  "rhs", 110);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs110), "lhs",
                  "lhs", 110);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgeru.p"),
                  "context", "context", 111);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.xger"),
                  "name", "name", 111);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 111);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "resolved", "resolved", 111);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c14_rhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs111), "rhs",
                  "rhs", 111);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs111), "lhs",
                  "lhs", 111);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "context", "context", 112);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 112);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 112);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c14_rhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs112), "rhs",
                  "rhs", 112);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs112), "lhs",
                  "lhs", 112);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 113);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 113);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c14_rhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs113), "rhs",
                  "rhs", 113);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs113), "lhs",
                  "lhs", 113);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 114);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.int"),
                  "name", "name", 114);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/int.p"),
                  "resolved", "resolved", 114);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c14_rhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs114), "rhs",
                  "rhs", 114);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs114), "lhs",
                  "lhs", 114);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 115);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("intmax"), "name", "name",
                  115);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 115);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c14_rhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs115), "rhs",
                  "rhs", 115);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs115), "lhs",
                  "lhs", 115);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 116);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("min"), "name", "name", 116);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 116);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 116);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1311262518U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c14_rhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs116), "rhs",
                  "rhs", 116);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs116), "lhs",
                  "lhs", 116);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 117);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 117);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 117);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 117);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c14_rhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs117), "rhs",
                  "rhs", 117);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs117), "lhs",
                  "lhs", 117);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 118);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 118);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 118);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 118);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c14_rhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs118), "rhs",
                  "rhs", 118);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs118), "lhs",
                  "lhs", 118);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 119);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 119);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 119);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 119);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c14_rhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs119), "rhs",
                  "rhs", 119);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs119), "lhs",
                  "lhs", 119);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 120);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 120);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 120);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 120);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c14_rhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs120), "rhs",
                  "rhs", 120);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs120), "lhs",
                  "lhs", 120);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 121);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 121);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 121);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 121);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c14_rhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs121), "rhs",
                  "rhs", 121);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs121), "lhs",
                  "lhs", 121);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "context", "context", 122);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.refblas.xger"),
                  "name", "name", 122);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xger.p"),
                  "resolved", "resolved", 122);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c14_rhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs122), "rhs",
                  "rhs", 122);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs122), "lhs",
                  "lhs", 122);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xger.p"),
                  "context", "context", 123);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.refblas.xgerx"), "name", "name", 123);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 123);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "resolved", "resolved", 123);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c14_rhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs123), "rhs",
                  "rhs", 123);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs123), "lhs",
                  "lhs", 123);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 124);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("abs"), "name", "name", 124);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 124);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 124);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c14_rhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs124), "rhs",
                  "rhs", 124);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs124), "lhs",
                  "lhs", 124);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 125);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 125);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 125);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 125);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c14_rhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs125), "rhs",
                  "rhs", 125);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs125), "lhs",
                  "lhs", 125);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 126);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 126);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 126);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c14_rhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs126), "rhs",
                  "rhs", 126);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs126), "lhs",
                  "lhs", 126);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 127);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 127);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 127);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c14_rhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs127), "rhs",
                  "rhs", 127);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs127), "lhs",
                  "lhs", 127);
  sf_mex_destroy(&c14_rhs64);
  sf_mex_destroy(&c14_lhs64);
  sf_mex_destroy(&c14_rhs65);
  sf_mex_destroy(&c14_lhs65);
  sf_mex_destroy(&c14_rhs66);
  sf_mex_destroy(&c14_lhs66);
  sf_mex_destroy(&c14_rhs67);
  sf_mex_destroy(&c14_lhs67);
  sf_mex_destroy(&c14_rhs68);
  sf_mex_destroy(&c14_lhs68);
  sf_mex_destroy(&c14_rhs69);
  sf_mex_destroy(&c14_lhs69);
  sf_mex_destroy(&c14_rhs70);
  sf_mex_destroy(&c14_lhs70);
  sf_mex_destroy(&c14_rhs71);
  sf_mex_destroy(&c14_lhs71);
  sf_mex_destroy(&c14_rhs72);
  sf_mex_destroy(&c14_lhs72);
  sf_mex_destroy(&c14_rhs73);
  sf_mex_destroy(&c14_lhs73);
  sf_mex_destroy(&c14_rhs74);
  sf_mex_destroy(&c14_lhs74);
  sf_mex_destroy(&c14_rhs75);
  sf_mex_destroy(&c14_lhs75);
  sf_mex_destroy(&c14_rhs76);
  sf_mex_destroy(&c14_lhs76);
  sf_mex_destroy(&c14_rhs77);
  sf_mex_destroy(&c14_lhs77);
  sf_mex_destroy(&c14_rhs78);
  sf_mex_destroy(&c14_lhs78);
  sf_mex_destroy(&c14_rhs79);
  sf_mex_destroy(&c14_lhs79);
  sf_mex_destroy(&c14_rhs80);
  sf_mex_destroy(&c14_lhs80);
  sf_mex_destroy(&c14_rhs81);
  sf_mex_destroy(&c14_lhs81);
  sf_mex_destroy(&c14_rhs82);
  sf_mex_destroy(&c14_lhs82);
  sf_mex_destroy(&c14_rhs83);
  sf_mex_destroy(&c14_lhs83);
  sf_mex_destroy(&c14_rhs84);
  sf_mex_destroy(&c14_lhs84);
  sf_mex_destroy(&c14_rhs85);
  sf_mex_destroy(&c14_lhs85);
  sf_mex_destroy(&c14_rhs86);
  sf_mex_destroy(&c14_lhs86);
  sf_mex_destroy(&c14_rhs87);
  sf_mex_destroy(&c14_lhs87);
  sf_mex_destroy(&c14_rhs88);
  sf_mex_destroy(&c14_lhs88);
  sf_mex_destroy(&c14_rhs89);
  sf_mex_destroy(&c14_lhs89);
  sf_mex_destroy(&c14_rhs90);
  sf_mex_destroy(&c14_lhs90);
  sf_mex_destroy(&c14_rhs91);
  sf_mex_destroy(&c14_lhs91);
  sf_mex_destroy(&c14_rhs92);
  sf_mex_destroy(&c14_lhs92);
  sf_mex_destroy(&c14_rhs93);
  sf_mex_destroy(&c14_lhs93);
  sf_mex_destroy(&c14_rhs94);
  sf_mex_destroy(&c14_lhs94);
  sf_mex_destroy(&c14_rhs95);
  sf_mex_destroy(&c14_lhs95);
  sf_mex_destroy(&c14_rhs96);
  sf_mex_destroy(&c14_lhs96);
  sf_mex_destroy(&c14_rhs97);
  sf_mex_destroy(&c14_lhs97);
  sf_mex_destroy(&c14_rhs98);
  sf_mex_destroy(&c14_lhs98);
  sf_mex_destroy(&c14_rhs99);
  sf_mex_destroy(&c14_lhs99);
  sf_mex_destroy(&c14_rhs100);
  sf_mex_destroy(&c14_lhs100);
  sf_mex_destroy(&c14_rhs101);
  sf_mex_destroy(&c14_lhs101);
  sf_mex_destroy(&c14_rhs102);
  sf_mex_destroy(&c14_lhs102);
  sf_mex_destroy(&c14_rhs103);
  sf_mex_destroy(&c14_lhs103);
  sf_mex_destroy(&c14_rhs104);
  sf_mex_destroy(&c14_lhs104);
  sf_mex_destroy(&c14_rhs105);
  sf_mex_destroy(&c14_lhs105);
  sf_mex_destroy(&c14_rhs106);
  sf_mex_destroy(&c14_lhs106);
  sf_mex_destroy(&c14_rhs107);
  sf_mex_destroy(&c14_lhs107);
  sf_mex_destroy(&c14_rhs108);
  sf_mex_destroy(&c14_lhs108);
  sf_mex_destroy(&c14_rhs109);
  sf_mex_destroy(&c14_lhs109);
  sf_mex_destroy(&c14_rhs110);
  sf_mex_destroy(&c14_lhs110);
  sf_mex_destroy(&c14_rhs111);
  sf_mex_destroy(&c14_lhs111);
  sf_mex_destroy(&c14_rhs112);
  sf_mex_destroy(&c14_lhs112);
  sf_mex_destroy(&c14_rhs113);
  sf_mex_destroy(&c14_lhs113);
  sf_mex_destroy(&c14_rhs114);
  sf_mex_destroy(&c14_lhs114);
  sf_mex_destroy(&c14_rhs115);
  sf_mex_destroy(&c14_lhs115);
  sf_mex_destroy(&c14_rhs116);
  sf_mex_destroy(&c14_lhs116);
  sf_mex_destroy(&c14_rhs117);
  sf_mex_destroy(&c14_lhs117);
  sf_mex_destroy(&c14_rhs118);
  sf_mex_destroy(&c14_lhs118);
  sf_mex_destroy(&c14_rhs119);
  sf_mex_destroy(&c14_lhs119);
  sf_mex_destroy(&c14_rhs120);
  sf_mex_destroy(&c14_lhs120);
  sf_mex_destroy(&c14_rhs121);
  sf_mex_destroy(&c14_lhs121);
  sf_mex_destroy(&c14_rhs122);
  sf_mex_destroy(&c14_lhs122);
  sf_mex_destroy(&c14_rhs123);
  sf_mex_destroy(&c14_lhs123);
  sf_mex_destroy(&c14_rhs124);
  sf_mex_destroy(&c14_lhs124);
  sf_mex_destroy(&c14_rhs125);
  sf_mex_destroy(&c14_lhs125);
  sf_mex_destroy(&c14_rhs126);
  sf_mex_destroy(&c14_lhs126);
  sf_mex_destroy(&c14_rhs127);
  sf_mex_destroy(&c14_lhs127);
}

static void c14_c_info_helper(const mxArray **c14_info)
{
  const mxArray *c14_rhs128 = NULL;
  const mxArray *c14_lhs128 = NULL;
  const mxArray *c14_rhs129 = NULL;
  const mxArray *c14_lhs129 = NULL;
  const mxArray *c14_rhs130 = NULL;
  const mxArray *c14_lhs130 = NULL;
  const mxArray *c14_rhs131 = NULL;
  const mxArray *c14_lhs131 = NULL;
  const mxArray *c14_rhs132 = NULL;
  const mxArray *c14_lhs132 = NULL;
  const mxArray *c14_rhs133 = NULL;
  const mxArray *c14_lhs133 = NULL;
  const mxArray *c14_rhs134 = NULL;
  const mxArray *c14_lhs134 = NULL;
  const mxArray *c14_rhs135 = NULL;
  const mxArray *c14_lhs135 = NULL;
  const mxArray *c14_rhs136 = NULL;
  const mxArray *c14_lhs136 = NULL;
  const mxArray *c14_rhs137 = NULL;
  const mxArray *c14_lhs137 = NULL;
  const mxArray *c14_rhs138 = NULL;
  const mxArray *c14_lhs138 = NULL;
  const mxArray *c14_rhs139 = NULL;
  const mxArray *c14_lhs139 = NULL;
  const mxArray *c14_rhs140 = NULL;
  const mxArray *c14_lhs140 = NULL;
  const mxArray *c14_rhs141 = NULL;
  const mxArray *c14_lhs141 = NULL;
  const mxArray *c14_rhs142 = NULL;
  const mxArray *c14_lhs142 = NULL;
  const mxArray *c14_rhs143 = NULL;
  const mxArray *c14_lhs143 = NULL;
  const mxArray *c14_rhs144 = NULL;
  const mxArray *c14_lhs144 = NULL;
  const mxArray *c14_rhs145 = NULL;
  const mxArray *c14_lhs145 = NULL;
  const mxArray *c14_rhs146 = NULL;
  const mxArray *c14_lhs146 = NULL;
  const mxArray *c14_rhs147 = NULL;
  const mxArray *c14_lhs147 = NULL;
  const mxArray *c14_rhs148 = NULL;
  const mxArray *c14_lhs148 = NULL;
  const mxArray *c14_rhs149 = NULL;
  const mxArray *c14_lhs149 = NULL;
  const mxArray *c14_rhs150 = NULL;
  const mxArray *c14_lhs150 = NULL;
  const mxArray *c14_rhs151 = NULL;
  const mxArray *c14_lhs151 = NULL;
  const mxArray *c14_rhs152 = NULL;
  const mxArray *c14_lhs152 = NULL;
  const mxArray *c14_rhs153 = NULL;
  const mxArray *c14_lhs153 = NULL;
  const mxArray *c14_rhs154 = NULL;
  const mxArray *c14_lhs154 = NULL;
  const mxArray *c14_rhs155 = NULL;
  const mxArray *c14_lhs155 = NULL;
  const mxArray *c14_rhs156 = NULL;
  const mxArray *c14_lhs156 = NULL;
  const mxArray *c14_rhs157 = NULL;
  const mxArray *c14_lhs157 = NULL;
  const mxArray *c14_rhs158 = NULL;
  const mxArray *c14_lhs158 = NULL;
  const mxArray *c14_rhs159 = NULL;
  const mxArray *c14_lhs159 = NULL;
  const mxArray *c14_rhs160 = NULL;
  const mxArray *c14_lhs160 = NULL;
  const mxArray *c14_rhs161 = NULL;
  const mxArray *c14_lhs161 = NULL;
  const mxArray *c14_rhs162 = NULL;
  const mxArray *c14_lhs162 = NULL;
  const mxArray *c14_rhs163 = NULL;
  const mxArray *c14_lhs163 = NULL;
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 128);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 128);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 128);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 128);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c14_rhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs128), "rhs",
                  "rhs", 128);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs128), "lhs",
                  "lhs", 128);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular"),
                  "context", "context", 129);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_warning"), "name",
                  "name", 129);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 129);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 129);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c14_rhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs129), "rhs",
                  "rhs", 129);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs129), "lhs",
                  "lhs", 129);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 130);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 130);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 130);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 130);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c14_rhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs130), "rhs",
                  "rhs", 130);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs130), "lhs",
                  "lhs", 130);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 131);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 131);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 131);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c14_rhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs131), "rhs",
                  "rhs", 131);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs131), "lhs",
                  "lhs", 131);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 132);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_xtrsm"), "name", "name",
                  132);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 132);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"),
                  "resolved", "resolved", 132);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c14_rhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs132), "rhs",
                  "rhs", 132);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs132), "lhs",
                  "lhs", 132);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 133);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 133);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 133);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 133);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c14_rhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs133), "rhs",
                  "rhs", 133);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs133), "lhs",
                  "lhs", 133);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 134);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.blas.xtrsm"),
                  "name", "name", 134);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "resolved", "resolved", 134);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c14_rhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs134), "rhs",
                  "rhs", 134);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs134), "lhs",
                  "lhs", 134);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 135);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 135);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 135);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c14_rhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs135), "rhs",
                  "rhs", 135);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs135), "lhs",
                  "lhs", 135);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p!below_threshold"),
                  "context", "context", 136);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 136);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 136);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c14_rhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs136), "rhs",
                  "rhs", 136);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs136), "lhs",
                  "lhs", 136);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 137);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 137);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 137);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c14_rhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs137), "rhs",
                  "rhs", 137);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs137), "lhs",
                  "lhs", 137);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 138);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.refblas.xtrsm"), "name", "name", 138);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "resolved", "resolved", 138);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c14_rhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs138), "rhs",
                  "rhs", 138);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs138), "lhs",
                  "lhs", 138);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 139);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 139);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 139);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c14_rhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs139), "rhs",
                  "rhs", 139);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs139), "lhs",
                  "lhs", 139);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 140);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 140);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 140);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c14_rhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs140), "rhs",
                  "rhs", 140);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs140), "lhs",
                  "lhs", 140);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 141);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 141);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 141);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 141);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c14_rhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs141), "rhs",
                  "rhs", 141);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs141), "lhs",
                  "lhs", 141);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 142);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("intmin"), "name", "name",
                  142);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 142);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c14_rhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs142), "rhs",
                  "rhs", 142);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs142), "lhs",
                  "lhs", 142);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 143);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("rdivide"), "name", "name",
                  143);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 143);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 143);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c14_rhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs143), "rhs",
                  "rhs", 143);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs143), "lhs",
                  "lhs", 143);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 144);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 144);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 144);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c14_rhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs144), "rhs",
                  "rhs", 144);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs144), "lhs",
                  "lhs", 144);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 145);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 145);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 145);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c14_rhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs145), "rhs",
                  "rhs", 145);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs145), "lhs",
                  "lhs", 145);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 146);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_div"), "name", "name",
                  146);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 146);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c14_rhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs146), "rhs",
                  "rhs", 146);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs146), "lhs",
                  "lhs", 146);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 147);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("log2"), "name", "name", 147);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log2.m"), "resolved",
                  "resolved", 147);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1343837582U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c14_rhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs147), "rhs",
                  "rhs", 147);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs147), "lhs",
                  "lhs", 147);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log2.m!scalar_frexp"),
                  "context", "context", 148);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("isfinite"), "name", "name",
                  148);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "resolved",
                  "resolved", 148);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c14_rhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs148), "rhs",
                  "rhs", 148);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs148), "lhs",
                  "lhs", 148);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 149);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 149);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 149);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c14_rhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs149), "rhs",
                  "rhs", 149);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs149), "lhs",
                  "lhs", 149);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 150);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("isinf"), "name", "name", 150);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 150);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c14_rhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs150), "rhs",
                  "rhs", 150);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs150), "lhs",
                  "lhs", 150);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 151);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 151);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 151);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c14_rhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs151), "rhs",
                  "rhs", 151);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs151), "lhs",
                  "lhs", 151);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 152);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("isnan"), "name", "name", 152);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 152);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c14_rhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs152), "rhs",
                  "rhs", 152);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs152), "lhs",
                  "lhs", 152);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 153);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("pow2"), "name", "name", 153);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/pow2.m"), "resolved",
                  "resolved", 153);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1343837582U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c14_rhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs153), "rhs",
                  "rhs", 153);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs153), "lhs",
                  "lhs", 153);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/pow2.m"), "context",
                  "context", 154);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_pow2"), "name",
                  "name", 154);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 154);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_pow2.m"),
                  "resolved", "resolved", 154);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286825932U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c14_rhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs154), "rhs",
                  "rhs", 154);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs154), "lhs",
                  "lhs", 154);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_pow2.m"),
                  "context", "context", 155);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("power"), "name", "name", 155);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 155);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c14_rhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs155), "rhs",
                  "rhs", 155);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs155), "lhs",
                  "lhs", 155);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 156);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 156);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 156);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c14_rhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs156), "rhs",
                  "rhs", 156);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs156), "lhs",
                  "lhs", 156);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 157);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 157);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 157);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c14_rhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs157), "rhs",
                  "rhs", 157);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs157), "lhs",
                  "lhs", 157);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 158);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 158);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 158);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 158);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c14_rhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs158), "rhs",
                  "rhs", 158);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs158), "lhs",
                  "lhs", 158);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 159);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("floor"), "name", "name", 159);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 159);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 159);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c14_rhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs159), "rhs",
                  "rhs", 159);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs159), "lhs",
                  "lhs", 159);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 160);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_error"), "name", "name",
                  160);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 160);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 160);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c14_rhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs160), "rhs",
                  "rhs", 160);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs160), "lhs",
                  "lhs", 160);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 161);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 161);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 161);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 161);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c14_rhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs161), "rhs",
                  "rhs", 161);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs161), "lhs",
                  "lhs", 161);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 162);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_div"), "name", "name",
                  162);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 162);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 162);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c14_rhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs162), "rhs",
                  "rhs", 162);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs162), "lhs",
                  "lhs", 162);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 163);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 163);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 163);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c14_rhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c14_lhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs163), "rhs",
                  "rhs", 163);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs163), "lhs",
                  "lhs", 163);
  sf_mex_destroy(&c14_rhs128);
  sf_mex_destroy(&c14_lhs128);
  sf_mex_destroy(&c14_rhs129);
  sf_mex_destroy(&c14_lhs129);
  sf_mex_destroy(&c14_rhs130);
  sf_mex_destroy(&c14_lhs130);
  sf_mex_destroy(&c14_rhs131);
  sf_mex_destroy(&c14_lhs131);
  sf_mex_destroy(&c14_rhs132);
  sf_mex_destroy(&c14_lhs132);
  sf_mex_destroy(&c14_rhs133);
  sf_mex_destroy(&c14_lhs133);
  sf_mex_destroy(&c14_rhs134);
  sf_mex_destroy(&c14_lhs134);
  sf_mex_destroy(&c14_rhs135);
  sf_mex_destroy(&c14_lhs135);
  sf_mex_destroy(&c14_rhs136);
  sf_mex_destroy(&c14_lhs136);
  sf_mex_destroy(&c14_rhs137);
  sf_mex_destroy(&c14_lhs137);
  sf_mex_destroy(&c14_rhs138);
  sf_mex_destroy(&c14_lhs138);
  sf_mex_destroy(&c14_rhs139);
  sf_mex_destroy(&c14_lhs139);
  sf_mex_destroy(&c14_rhs140);
  sf_mex_destroy(&c14_lhs140);
  sf_mex_destroy(&c14_rhs141);
  sf_mex_destroy(&c14_lhs141);
  sf_mex_destroy(&c14_rhs142);
  sf_mex_destroy(&c14_lhs142);
  sf_mex_destroy(&c14_rhs143);
  sf_mex_destroy(&c14_lhs143);
  sf_mex_destroy(&c14_rhs144);
  sf_mex_destroy(&c14_lhs144);
  sf_mex_destroy(&c14_rhs145);
  sf_mex_destroy(&c14_lhs145);
  sf_mex_destroy(&c14_rhs146);
  sf_mex_destroy(&c14_lhs146);
  sf_mex_destroy(&c14_rhs147);
  sf_mex_destroy(&c14_lhs147);
  sf_mex_destroy(&c14_rhs148);
  sf_mex_destroy(&c14_lhs148);
  sf_mex_destroy(&c14_rhs149);
  sf_mex_destroy(&c14_lhs149);
  sf_mex_destroy(&c14_rhs150);
  sf_mex_destroy(&c14_lhs150);
  sf_mex_destroy(&c14_rhs151);
  sf_mex_destroy(&c14_lhs151);
  sf_mex_destroy(&c14_rhs152);
  sf_mex_destroy(&c14_lhs152);
  sf_mex_destroy(&c14_rhs153);
  sf_mex_destroy(&c14_lhs153);
  sf_mex_destroy(&c14_rhs154);
  sf_mex_destroy(&c14_lhs154);
  sf_mex_destroy(&c14_rhs155);
  sf_mex_destroy(&c14_lhs155);
  sf_mex_destroy(&c14_rhs156);
  sf_mex_destroy(&c14_lhs156);
  sf_mex_destroy(&c14_rhs157);
  sf_mex_destroy(&c14_lhs157);
  sf_mex_destroy(&c14_rhs158);
  sf_mex_destroy(&c14_lhs158);
  sf_mex_destroy(&c14_rhs159);
  sf_mex_destroy(&c14_lhs159);
  sf_mex_destroy(&c14_rhs160);
  sf_mex_destroy(&c14_lhs160);
  sf_mex_destroy(&c14_rhs161);
  sf_mex_destroy(&c14_lhs161);
  sf_mex_destroy(&c14_rhs162);
  sf_mex_destroy(&c14_lhs162);
  sf_mex_destroy(&c14_rhs163);
  sf_mex_destroy(&c14_lhs163);
}

static void c14_expm(SFc14_untitledInstanceStruct *chartInstance, real_T c14_A
                     [16], real_T c14_F[16])
{
  real_T c14_normA;
  int32_T c14_j;
  real_T c14_b_j;
  real_T c14_s;
  int32_T c14_i;
  real_T c14_b_i;
  real_T c14_x;
  real_T c14_b_x;
  real_T c14_y;
  real_T c14_c_x;
  boolean_T c14_b;
  int32_T c14_c_i;
  real_T c14_d_i;
  static real_T c14_theta[5] = { 0.01495585217958292, 0.253939833006323,
    0.95041789961629319, 2.097847961257068, 5.3719203511481517 };

  int32_T c14_i18;
  real_T c14_b_A[16];
  static real_T c14_dv3[5] = { 3.0, 5.0, 7.0, 9.0, 13.0 };

  real_T c14_d_x;
  real_T c14_e_x;
  real_T c14_f_x;
  real_T c14_g_x;
  boolean_T c14_b_b;
  boolean_T c14_b0;
  real_T c14_h_x;
  boolean_T c14_c_b;
  boolean_T c14_b1;
  boolean_T c14_d_b;
  int32_T c14_eint;
  real_T c14_fdbl;
  int32_T c14_b_eint;
  real_T c14_b_fdbl;
  int32_T c14_c_eint;
  real_T c14_d1;
  real_T c14_d2;
  real_T c14_t;
  real_T c14_b_s;
  real_T c14_b_t;
  real_T c14_c_s;
  real_T c14_a;
  real_T c14_b_a;
  real_T c14_e_b;
  real_T c14_f_b;
  real_T c14_bk;
  real_T c14_g_b;
  real_T c14_br;
  real_T c14_b_y;
  real_T c14_c_y;
  real_T c14_d_y;
  int32_T c14_i19;
  int32_T c14_i20;
  real_T c14_c_A[16];
  real_T c14_d_s;
  int32_T c14_i21;
  int32_T c14_c_j;
  int32_T c14_i22;
  real_T c14_c_a[16];
  int32_T c14_i23;
  int32_T c14_i24;
  real_T c14_d_a[16];
  int32_T c14_i25;
  real_T c14_e_a[16];
  boolean_T exitg1;
  boolean_T exitg2;
  c14_normA = 0.0;
  c14_j = 0;
  exitg2 = false;
  while ((exitg2 == false) && (c14_j < 4)) {
    c14_b_j = 1.0 + (real_T)c14_j;
    c14_s = 0.0;
    for (c14_i = 0; c14_i < 4; c14_i++) {
      c14_b_i = 1.0 + (real_T)c14_i;
      c14_x = c14_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c14_b_i), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK
                       ("", (int32_T)_SFD_INTEGER_CHECK("", c14_b_j), 1, 4, 2, 0)
        - 1) << 2)) - 1];
      c14_b_x = c14_x;
      c14_y = muDoubleScalarAbs(c14_b_x);
      c14_s += c14_y;
    }

    c14_c_x = c14_s;
    c14_b = muDoubleScalarIsNaN(c14_c_x);
    if (c14_b) {
      c14_normA = rtNaN;
      exitg2 = true;
    } else {
      if (c14_s > c14_normA) {
        c14_normA = c14_s;
      }

      c14_j++;
    }
  }

  if (c14_normA <= 5.3719203511481517) {
    c14_c_i = 0;
    exitg1 = false;
    while ((exitg1 == false) && (c14_c_i < 5)) {
      c14_d_i = 1.0 + (real_T)c14_c_i;
      if (c14_normA <= c14_theta[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
           _SFD_INTEGER_CHECK("", c14_d_i), 1, 5, 1, 0) - 1]) {
        for (c14_i18 = 0; c14_i18 < 16; c14_i18++) {
          c14_b_A[c14_i18] = c14_A[c14_i18];
        }

        c14_PadeApproximantOfDegree(chartInstance, c14_b_A,
          c14_dv3[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c14_d_i), 1, 5, 1, 0) - 1], c14_F);
        exitg1 = true;
      } else {
        c14_c_i++;
      }
    }
  } else {
    c14_d_x = c14_normA / 5.3719203511481517;
    c14_e_x = c14_d_x;
    c14_f_x = c14_e_x;
    c14_g_x = c14_f_x;
    c14_b_b = muDoubleScalarIsInf(c14_g_x);
    c14_b0 = !c14_b_b;
    c14_h_x = c14_f_x;
    c14_c_b = muDoubleScalarIsNaN(c14_h_x);
    c14_b1 = !c14_c_b;
    c14_d_b = (c14_b0 && c14_b1);
    if (c14_d_b) {
      c14_fdbl = frexp(c14_e_x, &c14_eint);
      c14_b_eint = c14_eint;
      c14_b_fdbl = c14_fdbl;
      c14_c_eint = c14_b_eint;
      c14_d1 = c14_b_fdbl;
      c14_d2 = (real_T)c14_c_eint;
    } else {
      c14_d1 = c14_e_x;
      c14_d2 = 0.0;
    }

    c14_t = c14_d1;
    c14_b_s = c14_d2;
    c14_b_t = c14_t;
    c14_c_s = c14_b_s;
    if (c14_b_t == 0.5) {
      c14_c_s--;
    }

    c14_a = c14_c_s;
    c14_b_a = c14_a;
    c14_e_b = c14_b_a;
    c14_f_b = c14_e_b;
    c14_b_eml_scalar_eg(chartInstance);
    c14_bk = c14_f_b;
    c14_g_b = c14_bk;
    c14_b_eml_scalar_eg(chartInstance);
    c14_br = c14_g_b;
    c14_b_y = muDoubleScalarPower(2.0, c14_br);
    c14_c_y = c14_b_y;
    c14_d_y = c14_c_y;
    for (c14_i19 = 0; c14_i19 < 16; c14_i19++) {
      c14_A[c14_i19] /= c14_d_y;
    }

    for (c14_i20 = 0; c14_i20 < 16; c14_i20++) {
      c14_c_A[c14_i20] = c14_A[c14_i20];
    }

    c14_PadeApproximantOfDegree(chartInstance, c14_c_A, 13.0, c14_F);
    c14_d_s = c14_c_s;
    c14_i21 = (int32_T)c14_d_s;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c14_d_s, mxDOUBLE_CLASS, c14_i21);
    for (c14_c_j = 0; c14_c_j < c14_i21; c14_c_j++) {
      for (c14_i22 = 0; c14_i22 < 16; c14_i22++) {
        c14_c_a[c14_i22] = c14_F[c14_i22];
      }

      c14_eml_scalar_eg(chartInstance);
      c14_eml_scalar_eg(chartInstance);
      for (c14_i23 = 0; c14_i23 < 16; c14_i23++) {
        c14_F[c14_i23] = 0.0;
      }

      for (c14_i24 = 0; c14_i24 < 16; c14_i24++) {
        c14_d_a[c14_i24] = c14_c_a[c14_i24];
      }

      for (c14_i25 = 0; c14_i25 < 16; c14_i25++) {
        c14_e_a[c14_i25] = c14_c_a[c14_i25];
      }

      c14_b_eml_xgemm(chartInstance, c14_d_a, c14_e_a, c14_F);
    }
  }
}

static void c14_PadeApproximantOfDegree(SFc14_untitledInstanceStruct
  *chartInstance, real_T c14_A[16], real_T c14_m, real_T c14_F[16])
{
  int32_T c14_i26;
  real_T c14_A2[16];
  int32_T c14_i27;
  real_T c14_b_A[16];
  int32_T c14_i28;
  real_T c14_c_A[16];
  int32_T c14_i29;
  real_T c14_U[16];
  int32_T c14_k;
  real_T c14_b_k;
  int32_T c14_i30;
  real_T c14_y[16];
  int32_T c14_i31;
  int32_T c14_i32;
  real_T c14_d_A[16];
  int32_T c14_i33;
  real_T c14_b_y[16];
  int32_T c14_i34;
  real_T c14_d;
  int32_T c14_i35;
  real_T c14_A3[16];
  int32_T c14_i36;
  real_T c14_b_A2[16];
  int32_T c14_i37;
  real_T c14_c_A2[16];
  int32_T c14_i38;
  int32_T c14_i39;
  int32_T c14_c_k;
  int32_T c14_i40;
  int32_T c14_i41;
  int32_T c14_i42;
  real_T c14_e_A[16];
  int32_T c14_i43;
  real_T c14_c_y[16];
  int32_T c14_i44;
  int32_T c14_i45;
  int32_T c14_i46;
  int32_T c14_i47;
  real_T c14_A4[16];
  int32_T c14_i48;
  real_T c14_b_A3[16];
  int32_T c14_i49;
  real_T c14_d_A2[16];
  int32_T c14_i50;
  int32_T c14_i51;
  real_T c14_d_y[16];
  int32_T c14_i52;
  int32_T c14_d_k;
  int32_T c14_i53;
  int32_T c14_i54;
  int32_T c14_i55;
  real_T c14_f_A[16];
  int32_T c14_i56;
  real_T c14_e_y[16];
  int32_T c14_i57;
  int32_T c14_i58;
  int32_T c14_i59;
  int32_T c14_i60;
  int32_T c14_i61;
  int32_T c14_i62;
  real_T c14_b_A4[16];
  int32_T c14_i63;
  real_T c14_e_A2[16];
  int32_T c14_i64;
  int32_T c14_i65;
  int32_T c14_i66;
  real_T c14_f_y[16];
  int32_T c14_i67;
  int32_T c14_e_k;
  int32_T c14_i68;
  int32_T c14_i69;
  int32_T c14_i70;
  real_T c14_g_A[16];
  int32_T c14_i71;
  real_T c14_g_y[16];
  int32_T c14_i72;
  int32_T c14_i73;
  int32_T c14_i74;
  int32_T c14_i75;
  int32_T c14_i76;
  int32_T c14_i77;
  int32_T c14_i78;
  int32_T c14_i79;
  int32_T c14_i80;
  int32_T c14_f_k;
  int32_T c14_i81;
  int32_T c14_i82;
  int32_T c14_i83;
  int32_T c14_i84;
  int32_T c14_i85;
  real_T c14_c_A4[16];
  int32_T c14_i86;
  real_T c14_h_y[16];
  int32_T c14_i87;
  int32_T c14_i88;
  int32_T c14_i89;
  real_T c14_h_A[16];
  int32_T c14_i90;
  real_T c14_i_y[16];
  int32_T c14_i91;
  int32_T c14_i92;
  int32_T c14_i93;
  int32_T c14_i94;
  int32_T c14_i95;
  int32_T c14_i96;
  real_T c14_d_A4[16];
  int32_T c14_i97;
  real_T c14_j_y[16];
  int32_T c14_i98;
  int32_T c14_i99;
  int32_T c14_i100;
  int32_T c14_i101;
  int32_T c14_g_k;
  int32_T c14_h_k;
  real_T c14_uk;
  int32_T c14_info;
  int32_T c14_ipiv[4];
  int32_T c14_b_info;
  int32_T c14_c_info;
  int32_T c14_d_info;
  int32_T c14_xi;
  int32_T c14_b_xi;
  int32_T c14_ip;
  int32_T c14_xj;
  int32_T c14_b_xj;
  real_T c14_temp;
  int32_T c14_i102;
  real_T c14_b_U[16];
  int32_T c14_i103;
  real_T c14_c_U[16];
  c14_eml_scalar_eg(chartInstance);
  c14_eml_scalar_eg(chartInstance);
  for (c14_i26 = 0; c14_i26 < 16; c14_i26++) {
    c14_A2[c14_i26] = 0.0;
  }

  for (c14_i27 = 0; c14_i27 < 16; c14_i27++) {
    c14_b_A[c14_i27] = c14_A[c14_i27];
  }

  for (c14_i28 = 0; c14_i28 < 16; c14_i28++) {
    c14_c_A[c14_i28] = c14_A[c14_i28];
  }

  c14_b_eml_xgemm(chartInstance, c14_b_A, c14_c_A, c14_A2);
  if (c14_m == 3.0) {
    for (c14_i29 = 0; c14_i29 < 16; c14_i29++) {
      c14_U[c14_i29] = c14_A2[c14_i29];
    }

    for (c14_k = 0; c14_k < 4; c14_k++) {
      c14_b_k = 1.0 + (real_T)c14_k;
      c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
               c14_b_k), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 2, 0) - 1) << 2))
        - 1] = c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 1, 0) +
                      ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 2, 0) - 1) << 2)) - 1] + 60.0;
    }

    for (c14_i30 = 0; c14_i30 < 16; c14_i30++) {
      c14_y[c14_i30] = c14_U[c14_i30];
    }

    c14_eml_scalar_eg(chartInstance);
    c14_eml_scalar_eg(chartInstance);
    for (c14_i31 = 0; c14_i31 < 16; c14_i31++) {
      c14_U[c14_i31] = 0.0;
    }

    for (c14_i32 = 0; c14_i32 < 16; c14_i32++) {
      c14_d_A[c14_i32] = c14_A[c14_i32];
    }

    for (c14_i33 = 0; c14_i33 < 16; c14_i33++) {
      c14_b_y[c14_i33] = c14_y[c14_i33];
    }

    c14_b_eml_xgemm(chartInstance, c14_d_A, c14_b_y, c14_U);
    for (c14_i34 = 0; c14_i34 < 16; c14_i34++) {
      c14_F[c14_i34] = 12.0 * c14_A2[c14_i34];
    }

    c14_d = 120.0;
  } else {
    c14_eml_scalar_eg(chartInstance);
    c14_eml_scalar_eg(chartInstance);
    for (c14_i35 = 0; c14_i35 < 16; c14_i35++) {
      c14_A3[c14_i35] = 0.0;
    }

    for (c14_i36 = 0; c14_i36 < 16; c14_i36++) {
      c14_b_A2[c14_i36] = c14_A2[c14_i36];
    }

    for (c14_i37 = 0; c14_i37 < 16; c14_i37++) {
      c14_c_A2[c14_i37] = c14_A2[c14_i37];
    }

    c14_b_eml_xgemm(chartInstance, c14_b_A2, c14_c_A2, c14_A3);
    if (c14_m == 5.0) {
      for (c14_i38 = 0; c14_i38 < 16; c14_i38++) {
        c14_U[c14_i38] = 420.0 * c14_A2[c14_i38];
      }

      for (c14_i39 = 0; c14_i39 < 16; c14_i39++) {
        c14_U[c14_i39] += c14_A3[c14_i39];
      }

      for (c14_c_k = 0; c14_c_k < 4; c14_c_k++) {
        c14_b_k = 1.0 + (real_T)c14_c_k;
        c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 c14_b_k), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK("",
                  (int32_T)_SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 2, 0) - 1) <<
                2)) - 1] = c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 1, 0) +
          ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c14_b_k), 1, 4, 2, 0) - 1) << 2)) - 1] + 15120.0;
      }

      for (c14_i40 = 0; c14_i40 < 16; c14_i40++) {
        c14_y[c14_i40] = c14_U[c14_i40];
      }

      c14_eml_scalar_eg(chartInstance);
      c14_eml_scalar_eg(chartInstance);
      for (c14_i41 = 0; c14_i41 < 16; c14_i41++) {
        c14_U[c14_i41] = 0.0;
      }

      for (c14_i42 = 0; c14_i42 < 16; c14_i42++) {
        c14_e_A[c14_i42] = c14_A[c14_i42];
      }

      for (c14_i43 = 0; c14_i43 < 16; c14_i43++) {
        c14_c_y[c14_i43] = c14_y[c14_i43];
      }

      c14_b_eml_xgemm(chartInstance, c14_e_A, c14_c_y, c14_U);
      for (c14_i44 = 0; c14_i44 < 16; c14_i44++) {
        c14_A3[c14_i44] *= 30.0;
      }

      for (c14_i45 = 0; c14_i45 < 16; c14_i45++) {
        c14_A2[c14_i45] *= 3360.0;
      }

      for (c14_i46 = 0; c14_i46 < 16; c14_i46++) {
        c14_F[c14_i46] = c14_A3[c14_i46] + c14_A2[c14_i46];
      }

      c14_d = 30240.0;
    } else {
      c14_eml_scalar_eg(chartInstance);
      c14_eml_scalar_eg(chartInstance);
      for (c14_i47 = 0; c14_i47 < 16; c14_i47++) {
        c14_A4[c14_i47] = 0.0;
      }

      for (c14_i48 = 0; c14_i48 < 16; c14_i48++) {
        c14_b_A3[c14_i48] = c14_A3[c14_i48];
      }

      for (c14_i49 = 0; c14_i49 < 16; c14_i49++) {
        c14_d_A2[c14_i49] = c14_A2[c14_i49];
      }

      c14_b_eml_xgemm(chartInstance, c14_b_A3, c14_d_A2, c14_A4);
      if (c14_m == 7.0) {
        for (c14_i50 = 0; c14_i50 < 16; c14_i50++) {
          c14_U[c14_i50] = 1512.0 * c14_A3[c14_i50];
        }

        for (c14_i51 = 0; c14_i51 < 16; c14_i51++) {
          c14_d_y[c14_i51] = 277200.0 * c14_A2[c14_i51];
        }

        for (c14_i52 = 0; c14_i52 < 16; c14_i52++) {
          c14_U[c14_i52] = (c14_A4[c14_i52] + c14_U[c14_i52]) + c14_d_y[c14_i52];
        }

        for (c14_d_k = 0; c14_d_k < 4; c14_d_k++) {
          c14_b_k = 1.0 + (real_T)c14_d_k;
          c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   c14_b_k), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 2, 0) - 1) <<
                  2)) - 1] = c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 1, 0) +
            ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            c14_b_k), 1, 4, 2, 0) - 1) << 2)) - 1] + 8.64864E+6;
        }

        for (c14_i53 = 0; c14_i53 < 16; c14_i53++) {
          c14_y[c14_i53] = c14_U[c14_i53];
        }

        c14_eml_scalar_eg(chartInstance);
        c14_eml_scalar_eg(chartInstance);
        for (c14_i54 = 0; c14_i54 < 16; c14_i54++) {
          c14_U[c14_i54] = 0.0;
        }

        for (c14_i55 = 0; c14_i55 < 16; c14_i55++) {
          c14_f_A[c14_i55] = c14_A[c14_i55];
        }

        for (c14_i56 = 0; c14_i56 < 16; c14_i56++) {
          c14_e_y[c14_i56] = c14_y[c14_i56];
        }

        c14_b_eml_xgemm(chartInstance, c14_f_A, c14_e_y, c14_U);
        for (c14_i57 = 0; c14_i57 < 16; c14_i57++) {
          c14_A4[c14_i57] *= 56.0;
        }

        for (c14_i58 = 0; c14_i58 < 16; c14_i58++) {
          c14_A3[c14_i58] *= 25200.0;
        }

        for (c14_i59 = 0; c14_i59 < 16; c14_i59++) {
          c14_A2[c14_i59] *= 1.99584E+6;
        }

        for (c14_i60 = 0; c14_i60 < 16; c14_i60++) {
          c14_F[c14_i60] = (c14_A4[c14_i60] + c14_A3[c14_i60]) + c14_A2[c14_i60];
        }

        c14_d = 1.729728E+7;
      } else if (c14_m == 9.0) {
        c14_eml_scalar_eg(chartInstance);
        c14_eml_scalar_eg(chartInstance);
        for (c14_i61 = 0; c14_i61 < 16; c14_i61++) {
          c14_F[c14_i61] = 0.0;
        }

        for (c14_i62 = 0; c14_i62 < 16; c14_i62++) {
          c14_b_A4[c14_i62] = c14_A4[c14_i62];
        }

        for (c14_i63 = 0; c14_i63 < 16; c14_i63++) {
          c14_e_A2[c14_i63] = c14_A2[c14_i63];
        }

        c14_b_eml_xgemm(chartInstance, c14_b_A4, c14_e_A2, c14_F);
        for (c14_i64 = 0; c14_i64 < 16; c14_i64++) {
          c14_U[c14_i64] = 3960.0 * c14_A4[c14_i64];
        }

        for (c14_i65 = 0; c14_i65 < 16; c14_i65++) {
          c14_d_y[c14_i65] = 2.16216E+6 * c14_A3[c14_i65];
        }

        for (c14_i66 = 0; c14_i66 < 16; c14_i66++) {
          c14_f_y[c14_i66] = 3.027024E+8 * c14_A2[c14_i66];
        }

        for (c14_i67 = 0; c14_i67 < 16; c14_i67++) {
          c14_U[c14_i67] = ((c14_F[c14_i67] + c14_U[c14_i67]) + c14_d_y[c14_i67])
            + c14_f_y[c14_i67];
        }

        for (c14_e_k = 0; c14_e_k < 4; c14_e_k++) {
          c14_b_k = 1.0 + (real_T)c14_e_k;
          c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   c14_b_k), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 2, 0) - 1) <<
                  2)) - 1] = c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 1, 0) +
            ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            c14_b_k), 1, 4, 2, 0) - 1) << 2)) - 1] + 8.8216128E+9;
        }

        for (c14_i68 = 0; c14_i68 < 16; c14_i68++) {
          c14_y[c14_i68] = c14_U[c14_i68];
        }

        c14_eml_scalar_eg(chartInstance);
        c14_eml_scalar_eg(chartInstance);
        for (c14_i69 = 0; c14_i69 < 16; c14_i69++) {
          c14_U[c14_i69] = 0.0;
        }

        for (c14_i70 = 0; c14_i70 < 16; c14_i70++) {
          c14_g_A[c14_i70] = c14_A[c14_i70];
        }

        for (c14_i71 = 0; c14_i71 < 16; c14_i71++) {
          c14_g_y[c14_i71] = c14_y[c14_i71];
        }

        c14_b_eml_xgemm(chartInstance, c14_g_A, c14_g_y, c14_U);
        for (c14_i72 = 0; c14_i72 < 16; c14_i72++) {
          c14_F[c14_i72] *= 90.0;
        }

        for (c14_i73 = 0; c14_i73 < 16; c14_i73++) {
          c14_A4[c14_i73] *= 110880.0;
        }

        for (c14_i74 = 0; c14_i74 < 16; c14_i74++) {
          c14_A3[c14_i74] *= 3.027024E+7;
        }

        for (c14_i75 = 0; c14_i75 < 16; c14_i75++) {
          c14_A2[c14_i75] *= 2.0756736E+9;
        }

        for (c14_i76 = 0; c14_i76 < 16; c14_i76++) {
          c14_F[c14_i76] = ((c14_F[c14_i76] + c14_A4[c14_i76]) + c14_A3[c14_i76])
            + c14_A2[c14_i76];
        }

        c14_d = 1.76432256E+10;
      } else {
        for (c14_i77 = 0; c14_i77 < 16; c14_i77++) {
          c14_U[c14_i77] = 3.352212864E+10 * c14_A4[c14_i77];
        }

        for (c14_i78 = 0; c14_i78 < 16; c14_i78++) {
          c14_d_y[c14_i78] = 1.05594705216E+13 * c14_A3[c14_i78];
        }

        for (c14_i79 = 0; c14_i79 < 16; c14_i79++) {
          c14_f_y[c14_i79] = 1.1873537964288E+15 * c14_A2[c14_i79];
        }

        for (c14_i80 = 0; c14_i80 < 16; c14_i80++) {
          c14_U[c14_i80] = (c14_U[c14_i80] + c14_d_y[c14_i80]) + c14_f_y[c14_i80];
        }

        for (c14_f_k = 0; c14_f_k < 4; c14_f_k++) {
          c14_b_k = 1.0 + (real_T)c14_f_k;
          c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   c14_b_k), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 2, 0) - 1) <<
                  2)) - 1] = c14_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 1, 0) +
            ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            c14_b_k), 1, 4, 2, 0) - 1) << 2)) - 1] + 3.238237626624E+16;
        }

        for (c14_i81 = 0; c14_i81 < 16; c14_i81++) {
          c14_d_y[c14_i81] = 16380.0 * c14_A3[c14_i81];
        }

        for (c14_i82 = 0; c14_i82 < 16; c14_i82++) {
          c14_f_y[c14_i82] = 4.08408E+7 * c14_A2[c14_i82];
        }

        for (c14_i83 = 0; c14_i83 < 16; c14_i83++) {
          c14_d_y[c14_i83] = (c14_A4[c14_i83] + c14_d_y[c14_i83]) +
            c14_f_y[c14_i83];
        }

        c14_eml_scalar_eg(chartInstance);
        c14_eml_scalar_eg(chartInstance);
        for (c14_i84 = 0; c14_i84 < 16; c14_i84++) {
          c14_f_y[c14_i84] = 0.0;
        }

        for (c14_i85 = 0; c14_i85 < 16; c14_i85++) {
          c14_c_A4[c14_i85] = c14_A4[c14_i85];
        }

        for (c14_i86 = 0; c14_i86 < 16; c14_i86++) {
          c14_h_y[c14_i86] = c14_d_y[c14_i86];
        }

        c14_b_eml_xgemm(chartInstance, c14_c_A4, c14_h_y, c14_f_y);
        for (c14_i87 = 0; c14_i87 < 16; c14_i87++) {
          c14_f_y[c14_i87] += c14_U[c14_i87];
        }

        c14_eml_scalar_eg(chartInstance);
        c14_eml_scalar_eg(chartInstance);
        for (c14_i88 = 0; c14_i88 < 16; c14_i88++) {
          c14_U[c14_i88] = 0.0;
        }

        for (c14_i89 = 0; c14_i89 < 16; c14_i89++) {
          c14_h_A[c14_i89] = c14_A[c14_i89];
        }

        for (c14_i90 = 0; c14_i90 < 16; c14_i90++) {
          c14_i_y[c14_i90] = c14_f_y[c14_i90];
        }

        c14_b_eml_xgemm(chartInstance, c14_h_A, c14_i_y, c14_U);
        for (c14_i91 = 0; c14_i91 < 16; c14_i91++) {
          c14_d_y[c14_i91] = 182.0 * c14_A4[c14_i91];
        }

        for (c14_i92 = 0; c14_i92 < 16; c14_i92++) {
          c14_f_y[c14_i92] = 960960.0 * c14_A3[c14_i92];
        }

        for (c14_i93 = 0; c14_i93 < 16; c14_i93++) {
          c14_y[c14_i93] = 1.32324192E+9 * c14_A2[c14_i93];
        }

        for (c14_i94 = 0; c14_i94 < 16; c14_i94++) {
          c14_d_y[c14_i94] = (c14_d_y[c14_i94] + c14_f_y[c14_i94]) +
            c14_y[c14_i94];
        }

        c14_eml_scalar_eg(chartInstance);
        c14_eml_scalar_eg(chartInstance);
        for (c14_i95 = 0; c14_i95 < 16; c14_i95++) {
          c14_F[c14_i95] = 0.0;
        }

        for (c14_i96 = 0; c14_i96 < 16; c14_i96++) {
          c14_d_A4[c14_i96] = c14_A4[c14_i96];
        }

        for (c14_i97 = 0; c14_i97 < 16; c14_i97++) {
          c14_j_y[c14_i97] = c14_d_y[c14_i97];
        }

        c14_b_eml_xgemm(chartInstance, c14_d_A4, c14_j_y, c14_F);
        for (c14_i98 = 0; c14_i98 < 16; c14_i98++) {
          c14_A4[c14_i98] *= 6.704425728E+11;
        }

        for (c14_i99 = 0; c14_i99 < 16; c14_i99++) {
          c14_A3[c14_i99] *= 1.29060195264E+14;
        }

        for (c14_i100 = 0; c14_i100 < 16; c14_i100++) {
          c14_A2[c14_i100] *= 7.7717703038976E+15;
        }

        for (c14_i101 = 0; c14_i101 < 16; c14_i101++) {
          c14_F[c14_i101] = ((c14_F[c14_i101] + c14_A4[c14_i101]) +
                             c14_A3[c14_i101]) + c14_A2[c14_i101];
        }

        c14_d = 6.476475253248E+16;
      }
    }
  }

  for (c14_g_k = 0; c14_g_k < 4; c14_g_k++) {
    c14_b_k = 1.0 + (real_T)c14_g_k;
    c14_F[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
             c14_b_k), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 2, 0) - 1) << 2)) - 1] =
      c14_F[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
               c14_b_k), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c14_b_k), 1, 4, 2, 0) - 1) << 2))
      - 1] + c14_d;
  }

  for (c14_h_k = 0; c14_h_k < 16; c14_h_k++) {
    c14_b_k = 1.0 + (real_T)c14_h_k;
    c14_uk = c14_U[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", c14_b_k), 1, 16, 1, 0) - 1];
    c14_U[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c14_b_k), 1, 16, 1, 0) - 1] = c14_F[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", c14_b_k), 1, 16, 1, 0) - 1] - c14_uk;
    c14_F[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c14_b_k), 1, 16, 1, 0) - 1] = c14_F[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", c14_b_k), 1, 16, 1, 0) - 1] + c14_uk;
  }

  c14_b_eml_matlab_zgetrf(chartInstance, c14_U, c14_ipiv, &c14_info);
  c14_b_info = c14_info;
  c14_c_info = c14_b_info;
  c14_d_info = c14_c_info;
  if (c14_d_info > 0) {
    c14_eml_warning(chartInstance);
  }

  c14_eml_scalar_eg(chartInstance);
  for (c14_xi = 1; c14_xi < 4; c14_xi++) {
    c14_b_xi = c14_xi;
    if (c14_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c14_b_xi), 1, 4, 1, 0) - 1] != c14_b_xi) {
      c14_ip = c14_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c14_b_xi), 1, 4, 1, 0) - 1];
      for (c14_xj = 1; c14_xj < 5; c14_xj++) {
        c14_b_xj = c14_xj;
        c14_temp = c14_F[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c14_b_xi), 1, 4, 1, 0) +
                          ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c14_b_xj), 1, 4, 2, 0) - 1) << 2)) - 1];
        c14_F[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c14_b_xi), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK(
                  "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c14_b_xj), 1, 4, 2,
                  0) - 1) << 2)) - 1] = c14_F[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)_SFD_INTEGER_CHECK("", (real_T)c14_ip), 1, 4, 1, 0) +
          ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c14_b_xj), 1, 4, 2, 0) - 1) << 2)) - 1];
        c14_F[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c14_ip), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK("",
                  (int32_T)_SFD_INTEGER_CHECK("", (real_T)c14_b_xj), 1, 4, 2, 0)
                 - 1) << 2)) - 1] = c14_temp;
      }
    }
  }

  for (c14_i102 = 0; c14_i102 < 16; c14_i102++) {
    c14_b_U[c14_i102] = c14_U[c14_i102];
  }

  c14_c_eml_xtrsm(chartInstance, c14_b_U, c14_F);
  for (c14_i103 = 0; c14_i103 < 16; c14_i103++) {
    c14_c_U[c14_i103] = c14_U[c14_i103];
  }

  c14_d_eml_xtrsm(chartInstance, c14_c_U, c14_F);
}

static void c14_eml_scalar_eg(SFc14_untitledInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c14_eml_xgemm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16], real_T c14_C[16], real_T c14_b_C[16])
{
  int32_T c14_i104;
  int32_T c14_i105;
  real_T c14_b_A[16];
  int32_T c14_i106;
  real_T c14_b_B[16];
  for (c14_i104 = 0; c14_i104 < 16; c14_i104++) {
    c14_b_C[c14_i104] = c14_C[c14_i104];
  }

  for (c14_i105 = 0; c14_i105 < 16; c14_i105++) {
    c14_b_A[c14_i105] = c14_A[c14_i105];
  }

  for (c14_i106 = 0; c14_i106 < 16; c14_i106++) {
    c14_b_B[c14_i106] = c14_B[c14_i106];
  }

  c14_b_eml_xgemm(chartInstance, c14_b_A, c14_b_B, c14_b_C);
}

static void c14_eml_switch_helper(SFc14_untitledInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c14_eml_matlab_zgetrf(SFc14_untitledInstanceStruct *chartInstance,
  real_T c14_A[16], real_T c14_b_A[16], int32_T c14_ipiv[4], int32_T *c14_info)
{
  int32_T c14_i107;
  for (c14_i107 = 0; c14_i107 < 16; c14_i107++) {
    c14_b_A[c14_i107] = c14_A[c14_i107];
  }

  c14_b_eml_matlab_zgetrf(chartInstance, c14_b_A, c14_ipiv, c14_info);
}

static int32_T c14_eml_ixamax(SFc14_untitledInstanceStruct *chartInstance,
  int32_T c14_n, real_T c14_x[16], int32_T c14_ix0)
{
  int32_T c14_idxmax;
  int32_T c14_b_n;
  int32_T c14_b_ix0;
  int32_T c14_c_n;
  int32_T c14_c_ix0;
  int32_T c14_ix;
  real_T c14_b_x;
  real_T c14_c_x;
  real_T c14_d_x;
  real_T c14_y;
  real_T c14_e_x;
  real_T c14_f_x;
  real_T c14_b_y;
  real_T c14_smax;
  int32_T c14_d_n;
  int32_T c14_b;
  int32_T c14_b_b;
  boolean_T c14_overflow;
  int32_T c14_k;
  int32_T c14_b_k;
  int32_T c14_a;
  real_T c14_g_x;
  real_T c14_h_x;
  real_T c14_i_x;
  real_T c14_c_y;
  real_T c14_j_x;
  real_T c14_k_x;
  real_T c14_d_y;
  real_T c14_s;
  c14_b_n = c14_n;
  c14_b_ix0 = c14_ix0;
  c14_c_n = c14_b_n;
  c14_c_ix0 = c14_b_ix0;
  if (c14_c_n < 1) {
    c14_idxmax = 0;
  } else {
    c14_idxmax = 1;
    if (c14_c_n > 1) {
      c14_ix = c14_c_ix0;
      c14_b_x = c14_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c14_ix), 1, 16, 1, 0) - 1];
      c14_c_x = c14_b_x;
      c14_d_x = c14_c_x;
      c14_y = muDoubleScalarAbs(c14_d_x);
      c14_e_x = 0.0;
      c14_f_x = c14_e_x;
      c14_b_y = muDoubleScalarAbs(c14_f_x);
      c14_smax = c14_y + c14_b_y;
      c14_d_n = c14_c_n;
      c14_b = c14_d_n;
      c14_b_b = c14_b;
      if (2 > c14_b_b) {
        c14_overflow = false;
      } else {
        c14_eml_switch_helper(chartInstance);
        c14_overflow = (c14_b_b > 2147483646);
      }

      if (c14_overflow) {
        c14_check_forloop_overflow_error(chartInstance, c14_overflow);
      }

      for (c14_k = 2; c14_k <= c14_d_n; c14_k++) {
        c14_b_k = c14_k;
        c14_a = c14_ix + 1;
        c14_ix = c14_a;
        c14_g_x = c14_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c14_ix), 1, 16, 1, 0) - 1];
        c14_h_x = c14_g_x;
        c14_i_x = c14_h_x;
        c14_c_y = muDoubleScalarAbs(c14_i_x);
        c14_j_x = 0.0;
        c14_k_x = c14_j_x;
        c14_d_y = muDoubleScalarAbs(c14_k_x);
        c14_s = c14_c_y + c14_d_y;
        if (c14_s > c14_smax) {
          c14_idxmax = c14_b_k;
          c14_smax = c14_s;
        }
      }
    }
  }

  return c14_idxmax;
}

static void c14_check_forloop_overflow_error(SFc14_untitledInstanceStruct
  *chartInstance, boolean_T c14_overflow)
{
  int32_T c14_i108;
  static char_T c14_cv0[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c14_u[34];
  const mxArray *c14_y = NULL;
  int32_T c14_i109;
  static char_T c14_cv1[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c14_b_u[23];
  const mxArray *c14_b_y = NULL;
  (void)chartInstance;
  if (!c14_overflow) {
  } else {
    for (c14_i108 = 0; c14_i108 < 34; c14_i108++) {
      c14_u[c14_i108] = c14_cv0[c14_i108];
    }

    c14_y = NULL;
    sf_mex_assign(&c14_y, sf_mex_create("y", c14_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  false);
    for (c14_i109 = 0; c14_i109 < 23; c14_i109++) {
      c14_b_u[c14_i109] = c14_cv1[c14_i109];
    }

    c14_b_y = NULL;
    sf_mex_assign(&c14_b_y, sf_mex_create("y", c14_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 2U, 14, c14_y, 14, c14_b_y));
  }
}

static void c14_threshold(SFc14_untitledInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c14_eml_xgeru(SFc14_untitledInstanceStruct *chartInstance, int32_T
  c14_m, int32_T c14_n, real_T c14_alpha1, int32_T c14_ix0, int32_T c14_iy0,
  real_T c14_A[16], int32_T c14_ia0, real_T c14_b_A[16])
{
  int32_T c14_i110;
  for (c14_i110 = 0; c14_i110 < 16; c14_i110++) {
    c14_b_A[c14_i110] = c14_A[c14_i110];
  }

  c14_b_eml_xgeru(chartInstance, c14_m, c14_n, c14_alpha1, c14_ix0, c14_iy0,
                  c14_b_A, c14_ia0);
}

static void c14_b_eml_scalar_eg(SFc14_untitledInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c14_eml_warning(SFc14_untitledInstanceStruct *chartInstance)
{
  int32_T c14_i111;
  static char_T c14_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c14_u[27];
  const mxArray *c14_y = NULL;
  (void)chartInstance;
  for (c14_i111 = 0; c14_i111 < 27; c14_i111++) {
    c14_u[c14_i111] = c14_varargin_1[c14_i111];
  }

  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", c14_u, 10, 0U, 1U, 0U, 2, 1, 27),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c14_y));
}

static void c14_eml_xtrsm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16], real_T c14_b_B[16])
{
  int32_T c14_i112;
  int32_T c14_i113;
  real_T c14_b_A[16];
  for (c14_i112 = 0; c14_i112 < 16; c14_i112++) {
    c14_b_B[c14_i112] = c14_B[c14_i112];
  }

  for (c14_i113 = 0; c14_i113 < 16; c14_i113++) {
    c14_b_A[c14_i113] = c14_A[c14_i113];
  }

  c14_c_eml_xtrsm(chartInstance, c14_b_A, c14_b_B);
}

static void c14_b_threshold(SFc14_untitledInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c14_scalarEg(SFc14_untitledInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c14_b_eml_xtrsm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16], real_T c14_b_B[16])
{
  int32_T c14_i114;
  int32_T c14_i115;
  real_T c14_b_A[16];
  for (c14_i114 = 0; c14_i114 < 16; c14_i114++) {
    c14_b_B[c14_i114] = c14_B[c14_i114];
  }

  for (c14_i115 = 0; c14_i115 < 16; c14_i115++) {
    c14_b_A[c14_i115] = c14_A[c14_i115];
  }

  c14_d_eml_xtrsm(chartInstance, c14_b_A, c14_b_B);
}

static const mxArray *c14_c_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  int32_T c14_u;
  const mxArray *c14_y = NULL;
  SFc14_untitledInstanceStruct *chartInstance;
  chartInstance = (SFc14_untitledInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  c14_u = *(int32_T *)c14_inData;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", &c14_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, false);
  return c14_mxArrayOutData;
}

static int32_T c14_d_emlrt_marshallIn(SFc14_untitledInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId)
{
  int32_T c14_y;
  int32_T c14_i116;
  (void)chartInstance;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), &c14_i116, 1, 6, 0U, 0, 0U, 0);
  c14_y = c14_i116;
  sf_mex_destroy(&c14_u);
  return c14_y;
}

static void c14_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_b_sfEvent;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  int32_T c14_y;
  SFc14_untitledInstanceStruct *chartInstance;
  chartInstance = (SFc14_untitledInstanceStruct *)chartInstanceVoid;
  c14_b_sfEvent = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_y = c14_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_b_sfEvent),
    &c14_thisId);
  sf_mex_destroy(&c14_b_sfEvent);
  *(int32_T *)c14_outData = c14_y;
  sf_mex_destroy(&c14_mxArrayInData);
}

static uint8_T c14_e_emlrt_marshallIn(SFc14_untitledInstanceStruct
  *chartInstance, const mxArray *c14_b_is_active_c14_untitled, const char_T
  *c14_identifier)
{
  uint8_T c14_y;
  emlrtMsgIdentifier c14_thisId;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_y = c14_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c14_b_is_active_c14_untitled), &c14_thisId);
  sf_mex_destroy(&c14_b_is_active_c14_untitled);
  return c14_y;
}

static uint8_T c14_f_emlrt_marshallIn(SFc14_untitledInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId)
{
  uint8_T c14_y;
  uint8_T c14_u0;
  (void)chartInstance;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), &c14_u0, 1, 3, 0U, 0, 0U, 0);
  c14_y = c14_u0;
  sf_mex_destroy(&c14_u);
  return c14_y;
}

static void c14_b_eml_xgemm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16], real_T c14_C[16])
{
  int32_T c14_i117;
  int32_T c14_i118;
  int32_T c14_i119;
  int32_T c14_i120;
  int32_T c14_i121;
  (void)chartInstance;
  for (c14_i117 = 0; c14_i117 < 4; c14_i117++) {
    c14_i118 = 0;
    for (c14_i119 = 0; c14_i119 < 4; c14_i119++) {
      c14_C[c14_i118 + c14_i117] = 0.0;
      c14_i120 = 0;
      for (c14_i121 = 0; c14_i121 < 4; c14_i121++) {
        c14_C[c14_i118 + c14_i117] += c14_A[c14_i120 + c14_i117] *
          c14_B[c14_i121 + c14_i118];
        c14_i120 += 4;
      }

      c14_i118 += 4;
    }
  }
}

static void c14_b_eml_matlab_zgetrf(SFc14_untitledInstanceStruct *chartInstance,
  real_T c14_A[16], int32_T c14_ipiv[4], int32_T *c14_info)
{
  int32_T c14_i122;
  int32_T c14_j;
  int32_T c14_b_j;
  int32_T c14_a;
  int32_T c14_b_a;
  int32_T c14_jm1;
  int32_T c14_b;
  int32_T c14_b_b;
  int32_T c14_mmj;
  int32_T c14_c_a;
  int32_T c14_d_a;
  int32_T c14_c;
  int32_T c14_c_b;
  int32_T c14_d_b;
  int32_T c14_jj;
  int32_T c14_e_a;
  int32_T c14_f_a;
  int32_T c14_jp1j;
  int32_T c14_g_a;
  int32_T c14_h_a;
  int32_T c14_b_c;
  int32_T c14_i123;
  int32_T c14_i124;
  int32_T c14_i125;
  real_T c14_b_A[16];
  int32_T c14_i_a;
  int32_T c14_j_a;
  int32_T c14_jpiv_offset;
  int32_T c14_k_a;
  int32_T c14_e_b;
  int32_T c14_l_a;
  int32_T c14_f_b;
  int32_T c14_jpiv;
  int32_T c14_m_a;
  int32_T c14_g_b;
  int32_T c14_n_a;
  int32_T c14_h_b;
  int32_T c14_c_c;
  int32_T c14_i_b;
  int32_T c14_j_b;
  int32_T c14_jrow;
  int32_T c14_o_a;
  int32_T c14_k_b;
  int32_T c14_p_a;
  int32_T c14_l_b;
  int32_T c14_jprow;
  int32_T c14_ix0;
  int32_T c14_iy0;
  int32_T c14_b_ix0;
  int32_T c14_b_iy0;
  int32_T c14_c_ix0;
  int32_T c14_c_iy0;
  int32_T c14_ix;
  int32_T c14_iy;
  int32_T c14_k;
  real_T c14_temp;
  int32_T c14_q_a;
  int32_T c14_r_a;
  int32_T c14_b_jp1j;
  int32_T c14_s_a;
  int32_T c14_t_a;
  int32_T c14_d_c;
  int32_T c14_u_a;
  int32_T c14_m_b;
  int32_T c14_v_a;
  int32_T c14_n_b;
  int32_T c14_i126;
  int32_T c14_w_a;
  int32_T c14_o_b;
  int32_T c14_x_a;
  int32_T c14_p_b;
  boolean_T c14_overflow;
  int32_T c14_i;
  int32_T c14_b_i;
  real_T c14_x;
  real_T c14_y;
  real_T c14_b_x;
  real_T c14_b_y;
  real_T c14_z;
  int32_T c14_q_b;
  int32_T c14_r_b;
  int32_T c14_e_c;
  int32_T c14_y_a;
  int32_T c14_ab_a;
  int32_T c14_f_c;
  int32_T c14_bb_a;
  int32_T c14_cb_a;
  int32_T c14_g_c;
  real_T c14_d3;
  for (c14_i122 = 0; c14_i122 < 4; c14_i122++) {
    c14_ipiv[c14_i122] = 1 + c14_i122;
  }

  *c14_info = 0;
  for (c14_j = 1; c14_j < 4; c14_j++) {
    c14_b_j = c14_j;
    c14_a = c14_b_j;
    c14_b_a = c14_a - 1;
    c14_jm1 = c14_b_a;
    c14_b = c14_b_j;
    c14_b_b = c14_b;
    c14_mmj = 4 - c14_b_b;
    c14_c_a = c14_jm1;
    c14_d_a = c14_c_a;
    c14_c = c14_d_a * 5;
    c14_c_b = c14_c;
    c14_d_b = c14_c_b + 1;
    c14_jj = c14_d_b;
    c14_e_a = c14_jj;
    c14_f_a = c14_e_a + 1;
    c14_jp1j = c14_f_a;
    c14_g_a = c14_mmj;
    c14_h_a = c14_g_a;
    c14_b_c = c14_h_a;
    c14_i123 = 0;
    for (c14_i124 = 0; c14_i124 < 4; c14_i124++) {
      for (c14_i125 = 0; c14_i125 < 4; c14_i125++) {
        c14_b_A[c14_i125 + c14_i123] = c14_A[c14_i125 + c14_i123];
      }

      c14_i123 += 4;
    }

    c14_i_a = c14_eml_ixamax(chartInstance, c14_b_c + 1, c14_b_A, c14_jj);
    c14_j_a = c14_i_a - 1;
    c14_jpiv_offset = c14_j_a;
    c14_k_a = c14_jj;
    c14_e_b = c14_jpiv_offset;
    c14_l_a = c14_k_a;
    c14_f_b = c14_e_b;
    c14_jpiv = c14_l_a + c14_f_b;
    if (c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c14_jpiv), 1, 16, 1, 0) - 1] != 0.0) {
      if (c14_jpiv_offset != 0) {
        c14_m_a = c14_b_j;
        c14_g_b = c14_jpiv_offset;
        c14_n_a = c14_m_a;
        c14_h_b = c14_g_b;
        c14_c_c = c14_n_a + c14_h_b;
        c14_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c14_b_j), 1, 4, 1, 0) - 1] = c14_c_c;
        c14_i_b = c14_jm1;
        c14_j_b = c14_i_b + 1;
        c14_jrow = c14_j_b;
        c14_o_a = c14_jrow;
        c14_k_b = c14_jpiv_offset;
        c14_p_a = c14_o_a;
        c14_l_b = c14_k_b;
        c14_jprow = c14_p_a + c14_l_b;
        c14_ix0 = c14_jrow;
        c14_iy0 = c14_jprow;
        c14_b_ix0 = c14_ix0;
        c14_b_iy0 = c14_iy0;
        c14_threshold(chartInstance);
        c14_c_ix0 = c14_b_ix0;
        c14_c_iy0 = c14_b_iy0;
        c14_ix = c14_c_ix0;
        c14_iy = c14_c_iy0;
        for (c14_k = 1; c14_k < 5; c14_k++) {
          c14_temp = c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c14_ix), 1, 16, 1, 0) - 1];
          c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c14_ix), 1, 16, 1, 0) - 1] =
            c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c14_iy), 1, 16, 1, 0) - 1];
          c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c14_iy), 1, 16, 1, 0) - 1] = c14_temp;
          c14_q_a = c14_ix + 4;
          c14_ix = c14_q_a;
          c14_r_a = c14_iy + 4;
          c14_iy = c14_r_a;
        }
      }

      c14_b_jp1j = c14_jp1j;
      c14_s_a = c14_mmj;
      c14_t_a = c14_s_a;
      c14_d_c = c14_t_a;
      c14_u_a = c14_jp1j;
      c14_m_b = c14_d_c - 1;
      c14_v_a = c14_u_a;
      c14_n_b = c14_m_b;
      c14_i126 = c14_v_a + c14_n_b;
      c14_w_a = c14_b_jp1j;
      c14_o_b = c14_i126;
      c14_x_a = c14_w_a;
      c14_p_b = c14_o_b;
      if (c14_x_a > c14_p_b) {
        c14_overflow = false;
      } else {
        c14_eml_switch_helper(chartInstance);
        c14_overflow = (c14_p_b > 2147483646);
      }

      if (c14_overflow) {
        c14_check_forloop_overflow_error(chartInstance, c14_overflow);
      }

      for (c14_i = c14_b_jp1j; c14_i <= c14_i126; c14_i++) {
        c14_b_i = c14_i;
        c14_x = c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c14_b_i), 1, 16, 1, 0) - 1];
        c14_y = c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c14_jj), 1, 16, 1, 0) - 1];
        c14_b_x = c14_x;
        c14_b_y = c14_y;
        c14_z = c14_b_x / c14_b_y;
        c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c14_b_i), 1, 16, 1, 0) - 1] = c14_z;
      }
    } else {
      *c14_info = c14_b_j;
    }

    c14_q_b = c14_b_j;
    c14_r_b = c14_q_b;
    c14_e_c = 4 - c14_r_b;
    c14_y_a = c14_jj;
    c14_ab_a = c14_y_a;
    c14_f_c = c14_ab_a;
    c14_bb_a = c14_jj;
    c14_cb_a = c14_bb_a;
    c14_g_c = c14_cb_a;
    c14_d3 = -1.0;
    c14_b_eml_xgeru(chartInstance, c14_mmj, c14_e_c, c14_d3, c14_jp1j, c14_f_c +
                    4, c14_A, c14_g_c + 5);
  }

  if (*c14_info == 0) {
    if (!(c14_A[15] != 0.0)) {
      *c14_info = 4;
    }
  }
}

static void c14_b_eml_xgeru(SFc14_untitledInstanceStruct *chartInstance, int32_T
  c14_m, int32_T c14_n, real_T c14_alpha1, int32_T c14_ix0, int32_T c14_iy0,
  real_T c14_A[16], int32_T c14_ia0)
{
  int32_T c14_b_m;
  int32_T c14_b_n;
  real_T c14_b_alpha1;
  int32_T c14_b_ix0;
  int32_T c14_b_iy0;
  int32_T c14_b_ia0;
  int32_T c14_c_m;
  int32_T c14_c_n;
  real_T c14_c_alpha1;
  int32_T c14_c_ix0;
  int32_T c14_c_iy0;
  int32_T c14_c_ia0;
  int32_T c14_d_m;
  int32_T c14_d_n;
  real_T c14_d_alpha1;
  int32_T c14_d_ix0;
  int32_T c14_d_iy0;
  int32_T c14_d_ia0;
  int32_T c14_e_m;
  int32_T c14_e_n;
  real_T c14_e_alpha1;
  int32_T c14_e_ix0;
  int32_T c14_e_iy0;
  int32_T c14_e_ia0;
  int32_T c14_ixstart;
  int32_T c14_a;
  int32_T c14_jA;
  int32_T c14_jy;
  int32_T c14_f_n;
  int32_T c14_b;
  int32_T c14_b_b;
  boolean_T c14_overflow;
  int32_T c14_j;
  real_T c14_yjy;
  real_T c14_temp;
  int32_T c14_ix;
  int32_T c14_c_b;
  int32_T c14_i127;
  int32_T c14_b_a;
  int32_T c14_d_b;
  int32_T c14_i128;
  int32_T c14_c_a;
  int32_T c14_e_b;
  int32_T c14_d_a;
  int32_T c14_f_b;
  boolean_T c14_b_overflow;
  int32_T c14_ijA;
  int32_T c14_b_ijA;
  int32_T c14_e_a;
  int32_T c14_f_a;
  int32_T c14_g_a;
  c14_b_m = c14_m;
  c14_b_n = c14_n;
  c14_b_alpha1 = c14_alpha1;
  c14_b_ix0 = c14_ix0;
  c14_b_iy0 = c14_iy0;
  c14_b_ia0 = c14_ia0;
  c14_c_m = c14_b_m;
  c14_c_n = c14_b_n;
  c14_c_alpha1 = c14_b_alpha1;
  c14_c_ix0 = c14_b_ix0;
  c14_c_iy0 = c14_b_iy0;
  c14_c_ia0 = c14_b_ia0;
  c14_d_m = c14_c_m;
  c14_d_n = c14_c_n;
  c14_d_alpha1 = c14_c_alpha1;
  c14_d_ix0 = c14_c_ix0;
  c14_d_iy0 = c14_c_iy0;
  c14_d_ia0 = c14_c_ia0;
  c14_e_m = c14_d_m;
  c14_e_n = c14_d_n;
  c14_e_alpha1 = c14_d_alpha1;
  c14_e_ix0 = c14_d_ix0;
  c14_e_iy0 = c14_d_iy0;
  c14_e_ia0 = c14_d_ia0;
  if (c14_e_alpha1 == 0.0) {
  } else {
    c14_ixstart = c14_e_ix0;
    c14_a = c14_e_ia0 - 1;
    c14_jA = c14_a;
    c14_jy = c14_e_iy0;
    c14_f_n = c14_e_n;
    c14_b = c14_f_n;
    c14_b_b = c14_b;
    if (1 > c14_b_b) {
      c14_overflow = false;
    } else {
      c14_eml_switch_helper(chartInstance);
      c14_overflow = (c14_b_b > 2147483646);
    }

    if (c14_overflow) {
      c14_check_forloop_overflow_error(chartInstance, c14_overflow);
    }

    for (c14_j = 1; c14_j <= c14_f_n; c14_j++) {
      c14_yjy = c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c14_jy), 1, 16, 1, 0) - 1];
      if (c14_yjy != 0.0) {
        c14_temp = c14_yjy * c14_e_alpha1;
        c14_ix = c14_ixstart;
        c14_c_b = c14_jA + 1;
        c14_i127 = c14_c_b;
        c14_b_a = c14_e_m;
        c14_d_b = c14_jA;
        c14_i128 = c14_b_a + c14_d_b;
        c14_c_a = c14_i127;
        c14_e_b = c14_i128;
        c14_d_a = c14_c_a;
        c14_f_b = c14_e_b;
        if (c14_d_a > c14_f_b) {
          c14_b_overflow = false;
        } else {
          c14_eml_switch_helper(chartInstance);
          c14_b_overflow = (c14_f_b > 2147483646);
        }

        if (c14_b_overflow) {
          c14_check_forloop_overflow_error(chartInstance, c14_b_overflow);
        }

        for (c14_ijA = c14_i127; c14_ijA <= c14_i128; c14_ijA++) {
          c14_b_ijA = c14_ijA;
          c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c14_b_ijA), 1, 16, 1, 0) - 1] =
            c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c14_b_ijA), 1, 16, 1, 0) - 1] +
            c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c14_ix), 1, 16, 1, 0) - 1] * c14_temp;
          c14_e_a = c14_ix + 1;
          c14_ix = c14_e_a;
        }
      }

      c14_f_a = c14_jy + 4;
      c14_jy = c14_f_a;
      c14_g_a = c14_jA + 4;
      c14_jA = c14_g_a;
    }
  }
}

static void c14_c_eml_xtrsm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16])
{
  int32_T c14_j;
  int32_T c14_b_j;
  int32_T c14_jBcol;
  int32_T c14_k;
  int32_T c14_b_k;
  int32_T c14_kAcol;
  int32_T c14_i129;
  int32_T c14_a;
  int32_T c14_b_a;
  boolean_T c14_overflow;
  int32_T c14_i;
  int32_T c14_b_i;
  c14_b_threshold(chartInstance);
  c14_scalarEg(chartInstance);
  for (c14_j = 1; c14_j < 5; c14_j++) {
    c14_b_j = c14_j - 1;
    c14_jBcol = c14_b_j << 2;
    for (c14_k = 1; c14_k < 5; c14_k++) {
      c14_b_k = c14_k;
      c14_kAcol = (c14_b_k - 1) << 2;
      if (c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_k + c14_jBcol)), 1, 16, 1, 0) - 1] != 0.0) {
        c14_i129 = c14_b_k + 1;
        c14_a = c14_i129;
        c14_b_a = c14_a;
        if (c14_b_a > 4) {
          c14_overflow = false;
        } else {
          c14_eml_switch_helper(chartInstance);
          c14_overflow = false;
        }

        if (c14_overflow) {
          c14_check_forloop_overflow_error(chartInstance, c14_overflow);
        }

        for (c14_i = c14_i129; c14_i < 5; c14_i++) {
          c14_b_i = c14_i;
          c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_i + c14_jBcol)), 1, 16, 1, 0) - 1] =
            c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_i + c14_jBcol)), 1, 16, 1, 0) - 1] -
            c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_k + c14_jBcol)), 1, 16, 1, 0) - 1] *
            c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_i + c14_kAcol)), 1, 16, 1, 0) - 1];
        }
      }
    }
  }
}

static void c14_d_eml_xtrsm(SFc14_untitledInstanceStruct *chartInstance, real_T
  c14_A[16], real_T c14_B[16])
{
  int32_T c14_j;
  int32_T c14_b_j;
  int32_T c14_jBcol;
  int32_T c14_k;
  int32_T c14_b_k;
  int32_T c14_kAcol;
  real_T c14_x;
  real_T c14_y;
  real_T c14_b_x;
  real_T c14_b_y;
  real_T c14_c_x;
  real_T c14_c_y;
  real_T c14_z;
  int32_T c14_i130;
  int32_T c14_b;
  int32_T c14_b_b;
  boolean_T c14_overflow;
  int32_T c14_i;
  int32_T c14_b_i;
  c14_b_threshold(chartInstance);
  c14_scalarEg(chartInstance);
  for (c14_j = 1; c14_j < 5; c14_j++) {
    c14_b_j = c14_j - 1;
    c14_jBcol = c14_b_j << 2;
    for (c14_k = 4; c14_k > 0; c14_k--) {
      c14_b_k = c14_k;
      c14_kAcol = (c14_b_k - 1) << 2;
      if (c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_k + c14_jBcol)), 1, 16, 1, 0) - 1] != 0.0) {
        c14_x = c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)(c14_b_k + c14_jBcol)), 1, 16, 1, 0) -
          1];
        c14_y = c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)(c14_b_k + c14_kAcol)), 1, 16, 1, 0) -
          1];
        c14_b_x = c14_x;
        c14_b_y = c14_y;
        c14_c_x = c14_b_x;
        c14_c_y = c14_b_y;
        c14_z = c14_c_x / c14_c_y;
        c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)(c14_b_k + c14_jBcol)), 1, 16, 1, 0) - 1] = c14_z;
        c14_i130 = c14_b_k - 1;
        c14_b = c14_i130;
        c14_b_b = c14_b;
        if (1 > c14_b_b) {
          c14_overflow = false;
        } else {
          c14_eml_switch_helper(chartInstance);
          c14_overflow = (c14_b_b > 2147483646);
        }

        if (c14_overflow) {
          c14_check_forloop_overflow_error(chartInstance, c14_overflow);
        }

        for (c14_i = 1; c14_i <= c14_i130; c14_i++) {
          c14_b_i = c14_i;
          c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_i + c14_jBcol)), 1, 16, 1, 0) - 1] =
            c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_i + c14_jBcol)), 1, 16, 1, 0) - 1] -
            c14_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_k + c14_jBcol)), 1, 16, 1, 0) - 1] *
            c14_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c14_b_i + c14_kAcol)), 1, 16, 1, 0) - 1];
        }
      }
    }
  }
}

static void init_dsm_address_info(SFc14_untitledInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c14_untitled_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1262591729U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(540335555U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(980997694U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(54445756U);
}

mxArray *sf_c14_untitled_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("UvElzSK2mfRURATcmJQU0B");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(4);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(4);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c14_untitled_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c14_untitled_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c14_untitled(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"exponent_A\",},{M[8],M[0],T\"is_active_c14_untitled\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c14_untitled_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc14_untitledInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc14_untitledInstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _untitledMachineNumber_,
           14,
           1,
           1,
           0,
           2,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation(_untitledMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_untitledMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _untitledMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"A");
          _SFD_SET_DATA_PROPS(1,2,0,1,"exponent_A");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,68);

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_sf_marshallOut,(MexInFcnForType)
            c14_sf_marshallIn);
        }

        {
          real_T (*c14_A)[16];
          real_T (*c14_exponent_A)[16];
          c14_exponent_A = (real_T (*)[16])ssGetOutputPortSignal
            (chartInstance->S, 1);
          c14_A = (real_T (*)[16])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c14_A);
          _SFD_SET_DATA_VALUE_PTR(1U, *c14_exponent_A);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _untitledMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "qASixlSJ5ukrAvn9UDoRYE";
}

static void sf_opaque_initialize_c14_untitled(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc14_untitledInstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c14_untitled((SFc14_untitledInstanceStruct*)
    chartInstanceVar);
  initialize_c14_untitled((SFc14_untitledInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c14_untitled(void *chartInstanceVar)
{
  enable_c14_untitled((SFc14_untitledInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c14_untitled(void *chartInstanceVar)
{
  disable_c14_untitled((SFc14_untitledInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c14_untitled(void *chartInstanceVar)
{
  sf_gateway_c14_untitled((SFc14_untitledInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c14_untitled(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c14_untitled((SFc14_untitledInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c14_untitled();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c14_untitled(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c14_untitled();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c14_untitled((SFc14_untitledInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c14_untitled(SimStruct* S)
{
  return sf_internal_get_sim_state_c14_untitled(S);
}

static void sf_opaque_set_sim_state_c14_untitled(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c14_untitled(S, st);
}

static void sf_opaque_terminate_c14_untitled(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc14_untitledInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_untitled_optimization_info();
    }

    finalize_c14_untitled((SFc14_untitledInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc14_untitled((SFc14_untitledInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c14_untitled(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c14_untitled((SFc14_untitledInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c14_untitled(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_untitled_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      14);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,14,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,14,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,14);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,14,1);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,14,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 1; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,14);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(145181960U));
  ssSetChecksum1(S,(2538254662U));
  ssSetChecksum2(S,(51937714U));
  ssSetChecksum3(S,(315814594U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c14_untitled(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c14_untitled(SimStruct *S)
{
  SFc14_untitledInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc14_untitledInstanceStruct *)utMalloc(sizeof
    (SFc14_untitledInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc14_untitledInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c14_untitled;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c14_untitled;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c14_untitled;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c14_untitled;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c14_untitled;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c14_untitled;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c14_untitled;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c14_untitled;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c14_untitled;
  chartInstance->chartInfo.mdlStart = mdlStart_c14_untitled;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c14_untitled;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c14_untitled_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c14_untitled(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c14_untitled(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c14_untitled(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c14_untitled_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

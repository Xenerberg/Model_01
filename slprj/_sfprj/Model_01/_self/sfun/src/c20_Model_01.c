/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_01_sfun.h"
#include "c20_Model_01.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "Model_01_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c20_debug_family_names[8] = { "del_ita_0", "del_ita",
  "Tensor_del_ita", "nargin", "nargout", "del_ita_v", "nominal_ita", "ita" };

static const char * c20_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

static const char * c20_c_debug_family_names[8] = { "q_v", "q_0", "cross_q_v",
  "nargin", "nargout", "q", "flag", "CrossTensor" };

/* Function Declarations */
static void initialize_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance);
static void initialize_params_c20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance);
static void enable_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance);
static void disable_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance);
static void c20_update_debugger_state_c20_Model_01(SFc20_Model_01InstanceStruct *
  chartInstance);
static const mxArray *get_sim_state_c20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance);
static void set_sim_state_c20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance, const mxArray *c20_st);
static void finalize_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance);
static void sf_gateway_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance);
static void c20_chartstep_c20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance);
static void initSimStructsc20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c20_machineNumber, uint32_T
  c20_chartNumber, uint32_T c20_instanceNumber);
static const mxArray *c20_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static void c20_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_ita, const char_T *c20_identifier, real_T c20_y[4]);
static void c20_b_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[4]);
static void c20_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static const mxArray *c20_b_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static const mxArray *c20_c_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static real_T c20_c_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId);
static void c20_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static const mxArray *c20_d_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static void c20_d_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[16]);
static void c20_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static const mxArray *c20_e_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static void c20_e_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[9]);
static void c20_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static void c20_f_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[3]);
static void c20_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static void c20_info_helper(const mxArray **c20_info);
static const mxArray *c20_emlrt_marshallOut(const char * c20_u);
static const mxArray *c20_b_emlrt_marshallOut(const uint32_T c20_u);
static void c20_b_info_helper(const mxArray **c20_info);
static real_T c20_eml_xnrm2(SFc20_Model_01InstanceStruct *chartInstance, real_T
  c20_x[3]);
static void c20_eml_scalar_eg(SFc20_Model_01InstanceStruct *chartInstance);
static void c20_eml_error(SFc20_Model_01InstanceStruct *chartInstance);
static void c20_b_eml_scalar_eg(SFc20_Model_01InstanceStruct *chartInstance);
static void c20_eml_xgemm(SFc20_Model_01InstanceStruct *chartInstance, real_T
  c20_A[16], real_T c20_B[4], real_T c20_C[4], real_T c20_b_C[4]);
static const mxArray *c20_f_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static int32_T c20_g_emlrt_marshallIn(SFc20_Model_01InstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId);
static void c20_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static uint8_T c20_h_emlrt_marshallIn(SFc20_Model_01InstanceStruct
  *chartInstance, const mxArray *c20_b_is_active_c20_Model_01, const char_T
  *c20_identifier);
static uint8_T c20_i_emlrt_marshallIn(SFc20_Model_01InstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId);
static void c20_b_eml_xgemm(SFc20_Model_01InstanceStruct *chartInstance, real_T
  c20_A[16], real_T c20_B[4], real_T c20_C[4]);
static void init_dsm_address_info(SFc20_Model_01InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance)
{
  chartInstance->c20_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c20_is_active_c20_Model_01 = 0U;
}

static void initialize_params_c20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c20_update_debugger_state_c20_Model_01(SFc20_Model_01InstanceStruct *
  chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance)
{
  const mxArray *c20_st;
  const mxArray *c20_y = NULL;
  int32_T c20_i0;
  real_T c20_u[4];
  const mxArray *c20_b_y = NULL;
  uint8_T c20_hoistedGlobal;
  uint8_T c20_b_u;
  const mxArray *c20_c_y = NULL;
  real_T (*c20_ita)[4];
  c20_ita = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c20_st = NULL;
  c20_st = NULL;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_createcellmatrix(2, 1), false);
  for (c20_i0 = 0; c20_i0 < 4; c20_i0++) {
    c20_u[c20_i0] = (*c20_ita)[c20_i0];
  }

  c20_b_y = NULL;
  sf_mex_assign(&c20_b_y, sf_mex_create("y", c20_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_setcell(c20_y, 0, c20_b_y);
  c20_hoistedGlobal = chartInstance->c20_is_active_c20_Model_01;
  c20_b_u = c20_hoistedGlobal;
  c20_c_y = NULL;
  sf_mex_assign(&c20_c_y, sf_mex_create("y", &c20_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c20_y, 1, c20_c_y);
  sf_mex_assign(&c20_st, c20_y, false);
  return c20_st;
}

static void set_sim_state_c20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance, const mxArray *c20_st)
{
  const mxArray *c20_u;
  real_T c20_dv0[4];
  int32_T c20_i1;
  real_T (*c20_ita)[4];
  c20_ita = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c20_doneDoubleBufferReInit = true;
  c20_u = sf_mex_dup(c20_st);
  c20_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c20_u, 0)),
                       "ita", c20_dv0);
  for (c20_i1 = 0; c20_i1 < 4; c20_i1++) {
    (*c20_ita)[c20_i1] = c20_dv0[c20_i1];
  }

  chartInstance->c20_is_active_c20_Model_01 = c20_h_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c20_u, 1)),
     "is_active_c20_Model_01");
  sf_mex_destroy(&c20_u);
  c20_update_debugger_state_c20_Model_01(chartInstance);
  sf_mex_destroy(&c20_st);
}

static void finalize_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c20_Model_01(SFc20_Model_01InstanceStruct *chartInstance)
{
  int32_T c20_i2;
  int32_T c20_i3;
  int32_T c20_i4;
  real_T (*c20_nominal_ita)[4];
  real_T (*c20_ita)[4];
  real_T (*c20_del_ita_v)[3];
  c20_nominal_ita = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 1);
  c20_ita = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c20_del_ita_v = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 19U, chartInstance->c20_sfEvent);
  for (c20_i2 = 0; c20_i2 < 3; c20_i2++) {
    _SFD_DATA_RANGE_CHECK((*c20_del_ita_v)[c20_i2], 0U);
  }

  chartInstance->c20_sfEvent = CALL_EVENT;
  c20_chartstep_c20_Model_01(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_01MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c20_i3 = 0; c20_i3 < 4; c20_i3++) {
    _SFD_DATA_RANGE_CHECK((*c20_ita)[c20_i3], 1U);
  }

  for (c20_i4 = 0; c20_i4 < 4; c20_i4++) {
    _SFD_DATA_RANGE_CHECK((*c20_nominal_ita)[c20_i4], 2U);
  }
}

static void c20_chartstep_c20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance)
{
  int32_T c20_i5;
  real_T c20_del_ita_v[3];
  int32_T c20_i6;
  real_T c20_nominal_ita[4];
  uint32_T c20_debug_family_var_map[8];
  real_T c20_del_ita_0;
  real_T c20_del_ita[4];
  real_T c20_Tensor_del_ita[16];
  real_T c20_nargin = 2.0;
  real_T c20_nargout = 1.0;
  real_T c20_ita[4];
  int32_T c20_i7;
  real_T c20_x[3];
  int32_T c20_i8;
  real_T c20_b_x[3];
  real_T c20_y;
  int32_T c20_i9;
  int32_T c20_i10;
  real_T c20_c_x[3];
  real_T c20_b_y;
  int32_T c20_i11;
  real_T c20_B;
  real_T c20_c_y;
  real_T c20_d_y;
  real_T c20_e_y;
  int32_T c20_i12;
  int32_T c20_i13;
  int32_T c20_i14;
  real_T c20_d_x[3];
  real_T c20_f_y;
  real_T c20_a;
  real_T c20_b_a;
  real_T c20_c_a;
  real_T c20_ak;
  real_T c20_d_a;
  real_T c20_c;
  real_T c20_e_x;
  real_T c20_f_x;
  int32_T c20_i15;
  int32_T c20_i16;
  real_T c20_q[4];
  real_T c20_flag;
  real_T c20_q_v[3];
  real_T c20_q_0;
  real_T c20_cross_q_v[9];
  real_T c20_b_nargin = 2.0;
  real_T c20_b_nargout = 1.0;
  int32_T c20_i17;
  int32_T c20_i18;
  int32_T c20_i19;
  real_T c20_v[3];
  uint32_T c20_b_debug_family_var_map[4];
  real_T c20_c_nargin = 1.0;
  real_T c20_c_nargout = 1.0;
  int32_T c20_i20;
  real_T c20_e_a;
  int32_T c20_i21;
  static real_T c20_b[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c20_g_y[9];
  int32_T c20_i22;
  int32_T c20_i23;
  int32_T c20_i24;
  int32_T c20_i25;
  int32_T c20_i26;
  int32_T c20_i27;
  int32_T c20_i28;
  real_T c20_f_a;
  int32_T c20_i29;
  int32_T c20_i30;
  int32_T c20_i31;
  int32_T c20_i32;
  int32_T c20_i33;
  int32_T c20_i34;
  int32_T c20_i35;
  int32_T c20_i36;
  int32_T c20_i37;
  real_T c20_g_a[16];
  int32_T c20_i38;
  real_T c20_b_b[4];
  int32_T c20_i39;
  int32_T c20_i40;
  int32_T c20_i41;
  real_T c20_dv1[16];
  int32_T c20_i42;
  real_T c20_dv2[4];
  int32_T c20_i43;
  real_T c20_dv3[16];
  int32_T c20_i44;
  real_T c20_dv4[4];
  int32_T c20_i45;
  real_T (*c20_b_ita)[4];
  real_T (*c20_b_nominal_ita)[4];
  real_T (*c20_b_del_ita_v)[3];
  c20_b_nominal_ita = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 1);
  c20_b_ita = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c20_b_del_ita_v = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 19U, chartInstance->c20_sfEvent);
  for (c20_i5 = 0; c20_i5 < 3; c20_i5++) {
    c20_del_ita_v[c20_i5] = (*c20_b_del_ita_v)[c20_i5];
  }

  for (c20_i6 = 0; c20_i6 < 4; c20_i6++) {
    c20_nominal_ita[c20_i6] = (*c20_b_nominal_ita)[c20_i6];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c20_debug_family_names,
    c20_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_del_ita_0, 0U, c20_c_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_del_ita, 1U, c20_sf_marshallOut,
    c20_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_Tensor_del_ita, 2U,
    c20_d_sf_marshallOut, c20_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_nargin, 3U, c20_c_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_nargout, 4U, c20_c_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c20_del_ita_v, 5U, c20_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c20_nominal_ita, 6U, c20_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_ita, 7U, c20_sf_marshallOut,
    c20_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, 3);
  for (c20_i7 = 0; c20_i7 < 3; c20_i7++) {
    c20_x[c20_i7] = c20_del_ita_v[c20_i7];
  }

  for (c20_i8 = 0; c20_i8 < 3; c20_i8++) {
    c20_b_x[c20_i8] = c20_x[c20_i8];
  }

  c20_y = c20_eml_xnrm2(chartInstance, c20_b_x);
  if (CV_EML_IF(0, 1, 0, c20_y > 1.0)) {
    _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, 4);
    for (c20_i9 = 0; c20_i9 < 3; c20_i9++) {
      c20_x[c20_i9] = c20_del_ita_v[c20_i9];
    }

    for (c20_i10 = 0; c20_i10 < 3; c20_i10++) {
      c20_c_x[c20_i10] = c20_x[c20_i10];
    }

    c20_b_y = c20_eml_xnrm2(chartInstance, c20_c_x);
    for (c20_i11 = 0; c20_i11 < 3; c20_i11++) {
      c20_x[c20_i11] = c20_del_ita_v[c20_i11];
    }

    c20_B = c20_b_y;
    c20_c_y = c20_B;
    c20_d_y = c20_c_y;
    c20_e_y = c20_d_y;
    for (c20_i12 = 0; c20_i12 < 3; c20_i12++) {
      c20_del_ita_v[c20_i12] = c20_x[c20_i12] / c20_e_y;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, 6);
  for (c20_i13 = 0; c20_i13 < 3; c20_i13++) {
    c20_x[c20_i13] = c20_del_ita_v[c20_i13];
  }

  for (c20_i14 = 0; c20_i14 < 3; c20_i14++) {
    c20_d_x[c20_i14] = c20_x[c20_i14];
  }

  c20_f_y = c20_eml_xnrm2(chartInstance, c20_d_x);
  c20_a = c20_f_y;
  c20_b_a = c20_a;
  c20_c_a = c20_b_a;
  c20_eml_scalar_eg(chartInstance);
  c20_ak = c20_c_a;
  c20_d_a = c20_ak;
  c20_eml_scalar_eg(chartInstance);
  c20_c = c20_d_a * c20_d_a;
  c20_e_x = 1.0 - c20_c;
  c20_del_ita_0 = c20_e_x;
  if (c20_del_ita_0 < 0.0) {
    c20_eml_error(chartInstance);
  }

  c20_f_x = c20_del_ita_0;
  c20_del_ita_0 = c20_f_x;
  c20_del_ita_0 = muDoubleScalarSqrt(c20_del_ita_0);
  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, 7);
  for (c20_i15 = 0; c20_i15 < 3; c20_i15++) {
    c20_del_ita[c20_i15] = c20_del_ita_v[c20_i15];
  }

  c20_del_ita[3] = c20_del_ita_0;
  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, 8);
  for (c20_i16 = 0; c20_i16 < 4; c20_i16++) {
    c20_q[c20_i16] = c20_del_ita[c20_i16];
  }

  c20_flag = 1.0;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c20_c_debug_family_names,
    c20_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_q_v, 0U, c20_b_sf_marshallOut,
    c20_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_q_0, 1U, c20_c_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_cross_q_v, 2U, c20_e_sf_marshallOut,
    c20_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_b_nargin, 3U, c20_c_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_b_nargout, 4U, c20_c_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_q, 5U, c20_sf_marshallOut,
    c20_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_flag, 6U, c20_c_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_Tensor_del_ita, 7U,
    c20_d_sf_marshallOut, c20_c_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c20_sfEvent, 4);
  for (c20_i17 = 0; c20_i17 < 16; c20_i17++) {
    c20_Tensor_del_ita[c20_i17] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c20_sfEvent, 5);
  for (c20_i18 = 0; c20_i18 < 3; c20_i18++) {
    c20_q_v[c20_i18] = c20_q[c20_i18];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c20_sfEvent, 6);
  c20_q_0 = c20_q[3];
  _SFD_SCRIPT_CALL(0U, chartInstance->c20_sfEvent, 7);
  for (c20_i19 = 0; c20_i19 < 3; c20_i19++) {
    c20_v[c20_i19] = c20_q_v[c20_i19];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c20_b_debug_family_names,
    c20_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_c_nargin, 0U, c20_c_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_c_nargout, 1U, c20_c_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_v, 2U, c20_b_sf_marshallOut,
    c20_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_cross_q_v, 3U, c20_e_sf_marshallOut,
    c20_d_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 2);
  for (c20_i20 = 0; c20_i20 < 9; c20_i20++) {
    c20_cross_q_v[c20_i20] = 0.0;
  }

  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 3);
  c20_cross_q_v[0] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 4);
  c20_cross_q_v[4] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 5);
  c20_cross_q_v[8] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 6);
  c20_cross_q_v[3] = -c20_v[2];
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 7);
  c20_cross_q_v[6] = c20_v[1];
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 8);
  c20_cross_q_v[7] = -c20_v[0];
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 9);
  c20_cross_q_v[1] = c20_v[2];
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 10);
  c20_cross_q_v[2] = -c20_v[1];
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, 11);
  c20_cross_q_v[5] = c20_v[0];
  _SFD_SCRIPT_CALL(1U, chartInstance->c20_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_SCRIPT_CALL(0U, chartInstance->c20_sfEvent, 8);
  switch ((int32_T)_SFD_INTEGER_CHECK("flag", c20_flag)) {
   case 0:
    CV_SCRIPT_SWITCH(0, 0, 1);
    _SFD_SCRIPT_CALL(0U, chartInstance->c20_sfEvent, 10);
    c20_e_a = c20_q_0;
    for (c20_i21 = 0; c20_i21 < 9; c20_i21++) {
      c20_g_y[c20_i21] = c20_e_a * c20_b[c20_i21];
    }

    c20_i22 = 0;
    c20_i23 = 0;
    for (c20_i24 = 0; c20_i24 < 3; c20_i24++) {
      for (c20_i25 = 0; c20_i25 < 3; c20_i25++) {
        c20_Tensor_del_ita[c20_i25 + c20_i22] = -c20_cross_q_v[c20_i25 + c20_i23]
          + c20_g_y[c20_i25 + c20_i23];
      }

      c20_i22 += 4;
      c20_i23 += 3;
    }

    for (c20_i26 = 0; c20_i26 < 3; c20_i26++) {
      c20_Tensor_del_ita[c20_i26 + 12] = c20_q_v[c20_i26];
    }

    c20_i27 = 0;
    for (c20_i28 = 0; c20_i28 < 3; c20_i28++) {
      c20_Tensor_del_ita[c20_i27 + 3] = -c20_q_v[c20_i28];
      c20_i27 += 4;
    }

    c20_Tensor_del_ita[15] = c20_q_0;
    break;

   case 1:
    CV_SCRIPT_SWITCH(0, 0, 2);
    _SFD_SCRIPT_CALL(0U, chartInstance->c20_sfEvent, 12);
    c20_f_a = c20_q_0;
    for (c20_i29 = 0; c20_i29 < 9; c20_i29++) {
      c20_g_y[c20_i29] = c20_f_a * c20_b[c20_i29];
    }

    c20_i30 = 0;
    c20_i31 = 0;
    for (c20_i32 = 0; c20_i32 < 3; c20_i32++) {
      for (c20_i33 = 0; c20_i33 < 3; c20_i33++) {
        c20_Tensor_del_ita[c20_i33 + c20_i30] = c20_cross_q_v[c20_i33 + c20_i31]
          + c20_g_y[c20_i33 + c20_i31];
      }

      c20_i30 += 4;
      c20_i31 += 3;
    }

    for (c20_i34 = 0; c20_i34 < 3; c20_i34++) {
      c20_Tensor_del_ita[c20_i34 + 12] = c20_q_v[c20_i34];
    }

    c20_i35 = 0;
    for (c20_i36 = 0; c20_i36 < 3; c20_i36++) {
      c20_Tensor_del_ita[c20_i35 + 3] = -c20_q_v[c20_i36];
      c20_i35 += 4;
    }

    c20_Tensor_del_ita[15] = c20_q_0;
    break;

   default:
    CV_SCRIPT_SWITCH(0, 0, 0);
    break;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c20_sfEvent, -12);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, 9);
  for (c20_i37 = 0; c20_i37 < 16; c20_i37++) {
    c20_g_a[c20_i37] = c20_Tensor_del_ita[c20_i37];
  }

  for (c20_i38 = 0; c20_i38 < 4; c20_i38++) {
    c20_b_b[c20_i38] = c20_nominal_ita[c20_i38];
  }

  c20_b_eml_scalar_eg(chartInstance);
  c20_b_eml_scalar_eg(chartInstance);
  for (c20_i39 = 0; c20_i39 < 4; c20_i39++) {
    c20_ita[c20_i39] = 0.0;
  }

  for (c20_i40 = 0; c20_i40 < 4; c20_i40++) {
    c20_ita[c20_i40] = 0.0;
  }

  for (c20_i41 = 0; c20_i41 < 16; c20_i41++) {
    c20_dv1[c20_i41] = c20_g_a[c20_i41];
  }

  for (c20_i42 = 0; c20_i42 < 4; c20_i42++) {
    c20_dv2[c20_i42] = c20_b_b[c20_i42];
  }

  for (c20_i43 = 0; c20_i43 < 16; c20_i43++) {
    c20_dv3[c20_i43] = c20_dv1[c20_i43];
  }

  for (c20_i44 = 0; c20_i44 < 4; c20_i44++) {
    c20_dv4[c20_i44] = c20_dv2[c20_i44];
  }

  c20_b_eml_xgemm(chartInstance, c20_dv3, c20_dv4, c20_ita);
  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, -9);
  _SFD_SYMBOL_SCOPE_POP();
  for (c20_i45 = 0; c20_i45 < 4; c20_i45++) {
    (*c20_b_ita)[c20_i45] = c20_ita[c20_i45];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 19U, chartInstance->c20_sfEvent);
}

static void initSimStructsc20_Model_01(SFc20_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c20_machineNumber, uint32_T
  c20_chartNumber, uint32_T c20_instanceNumber)
{
  (void)c20_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c20_chartNumber, c20_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg\\Documents\\MATLAB\\Model_01\\fn_CrossTensor.m"));
  _SFD_SCRIPT_TRANSLATION(c20_chartNumber, c20_instanceNumber, 1U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c20_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  int32_T c20_i46;
  real_T c20_b_inData[4];
  int32_T c20_i47;
  real_T c20_u[4];
  const mxArray *c20_y = NULL;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  for (c20_i46 = 0; c20_i46 < 4; c20_i46++) {
    c20_b_inData[c20_i46] = (*(real_T (*)[4])c20_inData)[c20_i46];
  }

  for (c20_i47 = 0; c20_i47 < 4; c20_i47++) {
    c20_u[c20_i47] = c20_b_inData[c20_i47];
  }

  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, false);
  return c20_mxArrayOutData;
}

static void c20_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_ita, const char_T *c20_identifier, real_T c20_y[4])
{
  emlrtMsgIdentifier c20_thisId;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_ita), &c20_thisId, c20_y);
  sf_mex_destroy(&c20_ita);
}

static void c20_b_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[4])
{
  real_T c20_dv5[4];
  int32_T c20_i48;
  (void)chartInstance;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), c20_dv5, 1, 0, 0U, 1, 0U, 1, 4);
  for (c20_i48 = 0; c20_i48 < 4; c20_i48++) {
    c20_y[c20_i48] = c20_dv5[c20_i48];
  }

  sf_mex_destroy(&c20_u);
}

static void c20_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_ita;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  real_T c20_y[4];
  int32_T c20_i49;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_ita = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_ita), &c20_thisId, c20_y);
  sf_mex_destroy(&c20_ita);
  for (c20_i49 = 0; c20_i49 < 4; c20_i49++) {
    (*(real_T (*)[4])c20_outData)[c20_i49] = c20_y[c20_i49];
  }

  sf_mex_destroy(&c20_mxArrayInData);
}

static const mxArray *c20_b_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  int32_T c20_i50;
  real_T c20_b_inData[3];
  int32_T c20_i51;
  real_T c20_u[3];
  const mxArray *c20_y = NULL;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  for (c20_i50 = 0; c20_i50 < 3; c20_i50++) {
    c20_b_inData[c20_i50] = (*(real_T (*)[3])c20_inData)[c20_i50];
  }

  for (c20_i51 = 0; c20_i51 < 3; c20_i51++) {
    c20_u[c20_i51] = c20_b_inData[c20_i51];
  }

  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, false);
  return c20_mxArrayOutData;
}

static const mxArray *c20_c_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  real_T c20_u;
  const mxArray *c20_y = NULL;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  c20_u = *(real_T *)c20_inData;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", &c20_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, false);
  return c20_mxArrayOutData;
}

static real_T c20_c_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId)
{
  real_T c20_y;
  real_T c20_d0;
  (void)chartInstance;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), &c20_d0, 1, 0, 0U, 0, 0U, 0);
  c20_y = c20_d0;
  sf_mex_destroy(&c20_u);
  return c20_y;
}

static void c20_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_nargout;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  real_T c20_y;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_nargout = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_y = c20_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_nargout),
    &c20_thisId);
  sf_mex_destroy(&c20_nargout);
  *(real_T *)c20_outData = c20_y;
  sf_mex_destroy(&c20_mxArrayInData);
}

static const mxArray *c20_d_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  int32_T c20_i52;
  int32_T c20_i53;
  int32_T c20_i54;
  real_T c20_b_inData[16];
  int32_T c20_i55;
  int32_T c20_i56;
  int32_T c20_i57;
  real_T c20_u[16];
  const mxArray *c20_y = NULL;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  c20_i52 = 0;
  for (c20_i53 = 0; c20_i53 < 4; c20_i53++) {
    for (c20_i54 = 0; c20_i54 < 4; c20_i54++) {
      c20_b_inData[c20_i54 + c20_i52] = (*(real_T (*)[16])c20_inData)[c20_i54 +
        c20_i52];
    }

    c20_i52 += 4;
  }

  c20_i55 = 0;
  for (c20_i56 = 0; c20_i56 < 4; c20_i56++) {
    for (c20_i57 = 0; c20_i57 < 4; c20_i57++) {
      c20_u[c20_i57 + c20_i55] = c20_b_inData[c20_i57 + c20_i55];
    }

    c20_i55 += 4;
  }

  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, false);
  return c20_mxArrayOutData;
}

static void c20_d_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[16])
{
  real_T c20_dv6[16];
  int32_T c20_i58;
  (void)chartInstance;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), c20_dv6, 1, 0, 0U, 1, 0U, 2, 4,
                4);
  for (c20_i58 = 0; c20_i58 < 16; c20_i58++) {
    c20_y[c20_i58] = c20_dv6[c20_i58];
  }

  sf_mex_destroy(&c20_u);
}

static void c20_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_Tensor_del_ita;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  real_T c20_y[16];
  int32_T c20_i59;
  int32_T c20_i60;
  int32_T c20_i61;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_Tensor_del_ita = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_Tensor_del_ita),
    &c20_thisId, c20_y);
  sf_mex_destroy(&c20_Tensor_del_ita);
  c20_i59 = 0;
  for (c20_i60 = 0; c20_i60 < 4; c20_i60++) {
    for (c20_i61 = 0; c20_i61 < 4; c20_i61++) {
      (*(real_T (*)[16])c20_outData)[c20_i61 + c20_i59] = c20_y[c20_i61 +
        c20_i59];
    }

    c20_i59 += 4;
  }

  sf_mex_destroy(&c20_mxArrayInData);
}

static const mxArray *c20_e_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  int32_T c20_i62;
  int32_T c20_i63;
  int32_T c20_i64;
  real_T c20_b_inData[9];
  int32_T c20_i65;
  int32_T c20_i66;
  int32_T c20_i67;
  real_T c20_u[9];
  const mxArray *c20_y = NULL;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  c20_i62 = 0;
  for (c20_i63 = 0; c20_i63 < 3; c20_i63++) {
    for (c20_i64 = 0; c20_i64 < 3; c20_i64++) {
      c20_b_inData[c20_i64 + c20_i62] = (*(real_T (*)[9])c20_inData)[c20_i64 +
        c20_i62];
    }

    c20_i62 += 3;
  }

  c20_i65 = 0;
  for (c20_i66 = 0; c20_i66 < 3; c20_i66++) {
    for (c20_i67 = 0; c20_i67 < 3; c20_i67++) {
      c20_u[c20_i67 + c20_i65] = c20_b_inData[c20_i67 + c20_i65];
    }

    c20_i65 += 3;
  }

  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, false);
  return c20_mxArrayOutData;
}

static void c20_e_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[9])
{
  real_T c20_dv7[9];
  int32_T c20_i68;
  (void)chartInstance;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), c20_dv7, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c20_i68 = 0; c20_i68 < 9; c20_i68++) {
    c20_y[c20_i68] = c20_dv7[c20_i68];
  }

  sf_mex_destroy(&c20_u);
}

static void c20_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_SkewSymmetricTensor;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  real_T c20_y[9];
  int32_T c20_i69;
  int32_T c20_i70;
  int32_T c20_i71;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_SkewSymmetricTensor = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_SkewSymmetricTensor),
    &c20_thisId, c20_y);
  sf_mex_destroy(&c20_SkewSymmetricTensor);
  c20_i69 = 0;
  for (c20_i70 = 0; c20_i70 < 3; c20_i70++) {
    for (c20_i71 = 0; c20_i71 < 3; c20_i71++) {
      (*(real_T (*)[9])c20_outData)[c20_i71 + c20_i69] = c20_y[c20_i71 + c20_i69];
    }

    c20_i69 += 3;
  }

  sf_mex_destroy(&c20_mxArrayInData);
}

static void c20_f_emlrt_marshallIn(SFc20_Model_01InstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[3])
{
  real_T c20_dv8[3];
  int32_T c20_i72;
  (void)chartInstance;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), c20_dv8, 1, 0, 0U, 1, 0U, 1, 3);
  for (c20_i72 = 0; c20_i72 < 3; c20_i72++) {
    c20_y[c20_i72] = c20_dv8[c20_i72];
  }

  sf_mex_destroy(&c20_u);
}

static void c20_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_v;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  real_T c20_y[3];
  int32_T c20_i73;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_v = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_v), &c20_thisId, c20_y);
  sf_mex_destroy(&c20_v);
  for (c20_i73 = 0; c20_i73 < 3; c20_i73++) {
    (*(real_T (*)[3])c20_outData)[c20_i73] = c20_y[c20_i73];
  }

  sf_mex_destroy(&c20_mxArrayInData);
}

const mxArray *sf_c20_Model_01_get_eml_resolved_functions_info(void)
{
  const mxArray *c20_nameCaptureInfo = NULL;
  c20_nameCaptureInfo = NULL;
  sf_mex_assign(&c20_nameCaptureInfo, sf_mex_createstruct("structure", 2, 75, 1),
                false);
  c20_info_helper(&c20_nameCaptureInfo);
  c20_b_info_helper(&c20_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c20_nameCaptureInfo);
  return c20_nameCaptureInfo;
}

static void c20_info_helper(const mxArray **c20_info)
{
  const mxArray *c20_rhs0 = NULL;
  const mxArray *c20_lhs0 = NULL;
  const mxArray *c20_rhs1 = NULL;
  const mxArray *c20_lhs1 = NULL;
  const mxArray *c20_rhs2 = NULL;
  const mxArray *c20_lhs2 = NULL;
  const mxArray *c20_rhs3 = NULL;
  const mxArray *c20_lhs3 = NULL;
  const mxArray *c20_rhs4 = NULL;
  const mxArray *c20_lhs4 = NULL;
  const mxArray *c20_rhs5 = NULL;
  const mxArray *c20_lhs5 = NULL;
  const mxArray *c20_rhs6 = NULL;
  const mxArray *c20_lhs6 = NULL;
  const mxArray *c20_rhs7 = NULL;
  const mxArray *c20_lhs7 = NULL;
  const mxArray *c20_rhs8 = NULL;
  const mxArray *c20_lhs8 = NULL;
  const mxArray *c20_rhs9 = NULL;
  const mxArray *c20_lhs9 = NULL;
  const mxArray *c20_rhs10 = NULL;
  const mxArray *c20_lhs10 = NULL;
  const mxArray *c20_rhs11 = NULL;
  const mxArray *c20_lhs11 = NULL;
  const mxArray *c20_rhs12 = NULL;
  const mxArray *c20_lhs12 = NULL;
  const mxArray *c20_rhs13 = NULL;
  const mxArray *c20_lhs13 = NULL;
  const mxArray *c20_rhs14 = NULL;
  const mxArray *c20_lhs14 = NULL;
  const mxArray *c20_rhs15 = NULL;
  const mxArray *c20_lhs15 = NULL;
  const mxArray *c20_rhs16 = NULL;
  const mxArray *c20_lhs16 = NULL;
  const mxArray *c20_rhs17 = NULL;
  const mxArray *c20_lhs17 = NULL;
  const mxArray *c20_rhs18 = NULL;
  const mxArray *c20_lhs18 = NULL;
  const mxArray *c20_rhs19 = NULL;
  const mxArray *c20_lhs19 = NULL;
  const mxArray *c20_rhs20 = NULL;
  const mxArray *c20_lhs20 = NULL;
  const mxArray *c20_rhs21 = NULL;
  const mxArray *c20_lhs21 = NULL;
  const mxArray *c20_rhs22 = NULL;
  const mxArray *c20_lhs22 = NULL;
  const mxArray *c20_rhs23 = NULL;
  const mxArray *c20_lhs23 = NULL;
  const mxArray *c20_rhs24 = NULL;
  const mxArray *c20_lhs24 = NULL;
  const mxArray *c20_rhs25 = NULL;
  const mxArray *c20_lhs25 = NULL;
  const mxArray *c20_rhs26 = NULL;
  const mxArray *c20_lhs26 = NULL;
  const mxArray *c20_rhs27 = NULL;
  const mxArray *c20_lhs27 = NULL;
  const mxArray *c20_rhs28 = NULL;
  const mxArray *c20_lhs28 = NULL;
  const mxArray *c20_rhs29 = NULL;
  const mxArray *c20_lhs29 = NULL;
  const mxArray *c20_rhs30 = NULL;
  const mxArray *c20_lhs30 = NULL;
  const mxArray *c20_rhs31 = NULL;
  const mxArray *c20_lhs31 = NULL;
  const mxArray *c20_rhs32 = NULL;
  const mxArray *c20_lhs32 = NULL;
  const mxArray *c20_rhs33 = NULL;
  const mxArray *c20_lhs33 = NULL;
  const mxArray *c20_rhs34 = NULL;
  const mxArray *c20_lhs34 = NULL;
  const mxArray *c20_rhs35 = NULL;
  const mxArray *c20_lhs35 = NULL;
  const mxArray *c20_rhs36 = NULL;
  const mxArray *c20_lhs36 = NULL;
  const mxArray *c20_rhs37 = NULL;
  const mxArray *c20_lhs37 = NULL;
  const mxArray *c20_rhs38 = NULL;
  const mxArray *c20_lhs38 = NULL;
  const mxArray *c20_rhs39 = NULL;
  const mxArray *c20_lhs39 = NULL;
  const mxArray *c20_rhs40 = NULL;
  const mxArray *c20_lhs40 = NULL;
  const mxArray *c20_rhs41 = NULL;
  const mxArray *c20_lhs41 = NULL;
  const mxArray *c20_rhs42 = NULL;
  const mxArray *c20_lhs42 = NULL;
  const mxArray *c20_rhs43 = NULL;
  const mxArray *c20_lhs43 = NULL;
  const mxArray *c20_rhs44 = NULL;
  const mxArray *c20_lhs44 = NULL;
  const mxArray *c20_rhs45 = NULL;
  const mxArray *c20_lhs45 = NULL;
  const mxArray *c20_rhs46 = NULL;
  const mxArray *c20_lhs46 = NULL;
  const mxArray *c20_rhs47 = NULL;
  const mxArray *c20_lhs47 = NULL;
  const mxArray *c20_rhs48 = NULL;
  const mxArray *c20_lhs48 = NULL;
  const mxArray *c20_rhs49 = NULL;
  const mxArray *c20_lhs49 = NULL;
  const mxArray *c20_rhs50 = NULL;
  const mxArray *c20_lhs50 = NULL;
  const mxArray *c20_rhs51 = NULL;
  const mxArray *c20_lhs51 = NULL;
  const mxArray *c20_rhs52 = NULL;
  const mxArray *c20_lhs52 = NULL;
  const mxArray *c20_rhs53 = NULL;
  const mxArray *c20_lhs53 = NULL;
  const mxArray *c20_rhs54 = NULL;
  const mxArray *c20_lhs54 = NULL;
  const mxArray *c20_rhs55 = NULL;
  const mxArray *c20_lhs55 = NULL;
  const mxArray *c20_rhs56 = NULL;
  const mxArray *c20_lhs56 = NULL;
  const mxArray *c20_rhs57 = NULL;
  const mxArray *c20_lhs57 = NULL;
  const mxArray *c20_rhs58 = NULL;
  const mxArray *c20_lhs58 = NULL;
  const mxArray *c20_rhs59 = NULL;
  const mxArray *c20_lhs59 = NULL;
  const mxArray *c20_rhs60 = NULL;
  const mxArray *c20_lhs60 = NULL;
  const mxArray *c20_rhs61 = NULL;
  const mxArray *c20_lhs61 = NULL;
  const mxArray *c20_rhs62 = NULL;
  const mxArray *c20_lhs62 = NULL;
  const mxArray *c20_rhs63 = NULL;
  const mxArray *c20_lhs63 = NULL;
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("norm"), "name", "name", 0);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c20_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 1);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 1);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c20_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 2);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 2);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c20_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 3);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  3);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c20_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 4);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c20_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.blas.xnrm2"),
                  "name", "name", 5);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c20_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 6);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 6);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c20_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p!below_threshold"),
                  "context", "context", 7);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 7);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c20_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 8);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 8);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c20_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 9);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.refblas.xnrm2"), "name", "name", 9);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c20_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 10);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("realmin"), "name", "name",
                  10);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c20_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_realmin"), "name",
                  "name", 11);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1307658444U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c20_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 12);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c20_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 13);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 13);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c20_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 14);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 14);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 14);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c20_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 15);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 15);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 15);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c20_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 16);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 16);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c20_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 17);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("intmax"), "name", "name", 17);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c20_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 18);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c20_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 19);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("abs"), "name", "name", 19);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c20_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 20);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c20_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 21);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c20_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 22);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("mrdivide"), "name", "name",
                  22);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c20_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 23);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 23);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c20_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 24);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("rdivide"), "name", "name",
                  24);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c20_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 25);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c20_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 26);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c20_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_div"), "name", "name",
                  27);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c20_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 28);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c20_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 29);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("mpower"), "name", "name", 29);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c20_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 30);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 30);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c20_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("ismatrix"), "name", "name",
                  31);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 31);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c20_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("power"), "name", "name", 32);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c20_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 33);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 33);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c20_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 34);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 34);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 34);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c20_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 35);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c20_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 36);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 36);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c20_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 37);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 37);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c20_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 38);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("floor"), "name", "name", 38);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 38);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c20_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 39);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 39);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c20_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 40);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c20_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 41);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 41);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 41);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c20_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 42);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("sqrt"), "name", "name", 42);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c20_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 43);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_error"), "name", "name",
                  43);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 43);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c20_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 44);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 44);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286825938U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c20_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 45);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("fn_CrossTensor"), "name",
                  "name", 45);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1450348648U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c20_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_CrossTensor.m"), "context",
                  "context", 46);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 46);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1447321639U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c20_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_CrossTensor.m"), "context",
                  "context", 47);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eye"), "name", "name", 47);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 47);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c20_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 48);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 48);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c20_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 49);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 49);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c20_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 50);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("isinf"), "name", "name", 50);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 50);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c20_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 51);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 51);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c20_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 52);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_is_integer_class"),
                  "name", "name", 52);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c20_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 53);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("intmax"), "name", "name", 53);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 53);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c20_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 54);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("intmin"), "name", "name", 54);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 54);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c20_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 55);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 55);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c20_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 56);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.indexIntRelop"), "name", "name", 56);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c20_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 57);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 57);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c20_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 58);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 58);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c20_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 59);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("intmin"), "name", "name", 59);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 59);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c20_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 60);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 60);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c20_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 61);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("intmax"), "name", "name", 61);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 61);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c20_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 62);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 62);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c20_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_CrossTensor.m"), "context",
                  "context", 63);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 63);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c20_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c20_rhs0);
  sf_mex_destroy(&c20_lhs0);
  sf_mex_destroy(&c20_rhs1);
  sf_mex_destroy(&c20_lhs1);
  sf_mex_destroy(&c20_rhs2);
  sf_mex_destroy(&c20_lhs2);
  sf_mex_destroy(&c20_rhs3);
  sf_mex_destroy(&c20_lhs3);
  sf_mex_destroy(&c20_rhs4);
  sf_mex_destroy(&c20_lhs4);
  sf_mex_destroy(&c20_rhs5);
  sf_mex_destroy(&c20_lhs5);
  sf_mex_destroy(&c20_rhs6);
  sf_mex_destroy(&c20_lhs6);
  sf_mex_destroy(&c20_rhs7);
  sf_mex_destroy(&c20_lhs7);
  sf_mex_destroy(&c20_rhs8);
  sf_mex_destroy(&c20_lhs8);
  sf_mex_destroy(&c20_rhs9);
  sf_mex_destroy(&c20_lhs9);
  sf_mex_destroy(&c20_rhs10);
  sf_mex_destroy(&c20_lhs10);
  sf_mex_destroy(&c20_rhs11);
  sf_mex_destroy(&c20_lhs11);
  sf_mex_destroy(&c20_rhs12);
  sf_mex_destroy(&c20_lhs12);
  sf_mex_destroy(&c20_rhs13);
  sf_mex_destroy(&c20_lhs13);
  sf_mex_destroy(&c20_rhs14);
  sf_mex_destroy(&c20_lhs14);
  sf_mex_destroy(&c20_rhs15);
  sf_mex_destroy(&c20_lhs15);
  sf_mex_destroy(&c20_rhs16);
  sf_mex_destroy(&c20_lhs16);
  sf_mex_destroy(&c20_rhs17);
  sf_mex_destroy(&c20_lhs17);
  sf_mex_destroy(&c20_rhs18);
  sf_mex_destroy(&c20_lhs18);
  sf_mex_destroy(&c20_rhs19);
  sf_mex_destroy(&c20_lhs19);
  sf_mex_destroy(&c20_rhs20);
  sf_mex_destroy(&c20_lhs20);
  sf_mex_destroy(&c20_rhs21);
  sf_mex_destroy(&c20_lhs21);
  sf_mex_destroy(&c20_rhs22);
  sf_mex_destroy(&c20_lhs22);
  sf_mex_destroy(&c20_rhs23);
  sf_mex_destroy(&c20_lhs23);
  sf_mex_destroy(&c20_rhs24);
  sf_mex_destroy(&c20_lhs24);
  sf_mex_destroy(&c20_rhs25);
  sf_mex_destroy(&c20_lhs25);
  sf_mex_destroy(&c20_rhs26);
  sf_mex_destroy(&c20_lhs26);
  sf_mex_destroy(&c20_rhs27);
  sf_mex_destroy(&c20_lhs27);
  sf_mex_destroy(&c20_rhs28);
  sf_mex_destroy(&c20_lhs28);
  sf_mex_destroy(&c20_rhs29);
  sf_mex_destroy(&c20_lhs29);
  sf_mex_destroy(&c20_rhs30);
  sf_mex_destroy(&c20_lhs30);
  sf_mex_destroy(&c20_rhs31);
  sf_mex_destroy(&c20_lhs31);
  sf_mex_destroy(&c20_rhs32);
  sf_mex_destroy(&c20_lhs32);
  sf_mex_destroy(&c20_rhs33);
  sf_mex_destroy(&c20_lhs33);
  sf_mex_destroy(&c20_rhs34);
  sf_mex_destroy(&c20_lhs34);
  sf_mex_destroy(&c20_rhs35);
  sf_mex_destroy(&c20_lhs35);
  sf_mex_destroy(&c20_rhs36);
  sf_mex_destroy(&c20_lhs36);
  sf_mex_destroy(&c20_rhs37);
  sf_mex_destroy(&c20_lhs37);
  sf_mex_destroy(&c20_rhs38);
  sf_mex_destroy(&c20_lhs38);
  sf_mex_destroy(&c20_rhs39);
  sf_mex_destroy(&c20_lhs39);
  sf_mex_destroy(&c20_rhs40);
  sf_mex_destroy(&c20_lhs40);
  sf_mex_destroy(&c20_rhs41);
  sf_mex_destroy(&c20_lhs41);
  sf_mex_destroy(&c20_rhs42);
  sf_mex_destroy(&c20_lhs42);
  sf_mex_destroy(&c20_rhs43);
  sf_mex_destroy(&c20_lhs43);
  sf_mex_destroy(&c20_rhs44);
  sf_mex_destroy(&c20_lhs44);
  sf_mex_destroy(&c20_rhs45);
  sf_mex_destroy(&c20_lhs45);
  sf_mex_destroy(&c20_rhs46);
  sf_mex_destroy(&c20_lhs46);
  sf_mex_destroy(&c20_rhs47);
  sf_mex_destroy(&c20_lhs47);
  sf_mex_destroy(&c20_rhs48);
  sf_mex_destroy(&c20_lhs48);
  sf_mex_destroy(&c20_rhs49);
  sf_mex_destroy(&c20_lhs49);
  sf_mex_destroy(&c20_rhs50);
  sf_mex_destroy(&c20_lhs50);
  sf_mex_destroy(&c20_rhs51);
  sf_mex_destroy(&c20_lhs51);
  sf_mex_destroy(&c20_rhs52);
  sf_mex_destroy(&c20_lhs52);
  sf_mex_destroy(&c20_rhs53);
  sf_mex_destroy(&c20_lhs53);
  sf_mex_destroy(&c20_rhs54);
  sf_mex_destroy(&c20_lhs54);
  sf_mex_destroy(&c20_rhs55);
  sf_mex_destroy(&c20_lhs55);
  sf_mex_destroy(&c20_rhs56);
  sf_mex_destroy(&c20_lhs56);
  sf_mex_destroy(&c20_rhs57);
  sf_mex_destroy(&c20_lhs57);
  sf_mex_destroy(&c20_rhs58);
  sf_mex_destroy(&c20_lhs58);
  sf_mex_destroy(&c20_rhs59);
  sf_mex_destroy(&c20_lhs59);
  sf_mex_destroy(&c20_rhs60);
  sf_mex_destroy(&c20_lhs60);
  sf_mex_destroy(&c20_rhs61);
  sf_mex_destroy(&c20_lhs61);
  sf_mex_destroy(&c20_rhs62);
  sf_mex_destroy(&c20_lhs62);
  sf_mex_destroy(&c20_rhs63);
  sf_mex_destroy(&c20_lhs63);
}

static const mxArray *c20_emlrt_marshallOut(const char * c20_u)
{
  const mxArray *c20_y = NULL;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c20_u)), false);
  return c20_y;
}

static const mxArray *c20_b_emlrt_marshallOut(const uint32_T c20_u)
{
  const mxArray *c20_y = NULL;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", &c20_u, 7, 0U, 0U, 0U, 0), false);
  return c20_y;
}

static void c20_b_info_helper(const mxArray **c20_info)
{
  const mxArray *c20_rhs64 = NULL;
  const mxArray *c20_lhs64 = NULL;
  const mxArray *c20_rhs65 = NULL;
  const mxArray *c20_lhs65 = NULL;
  const mxArray *c20_rhs66 = NULL;
  const mxArray *c20_lhs66 = NULL;
  const mxArray *c20_rhs67 = NULL;
  const mxArray *c20_lhs67 = NULL;
  const mxArray *c20_rhs68 = NULL;
  const mxArray *c20_lhs68 = NULL;
  const mxArray *c20_rhs69 = NULL;
  const mxArray *c20_lhs69 = NULL;
  const mxArray *c20_rhs70 = NULL;
  const mxArray *c20_lhs70 = NULL;
  const mxArray *c20_rhs71 = NULL;
  const mxArray *c20_lhs71 = NULL;
  const mxArray *c20_rhs72 = NULL;
  const mxArray *c20_lhs72 = NULL;
  const mxArray *c20_rhs73 = NULL;
  const mxArray *c20_lhs73 = NULL;
  const mxArray *c20_rhs74 = NULL;
  const mxArray *c20_lhs74 = NULL;
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 64);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 64);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c20_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 65);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 65);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c20_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 66);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 66);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 66);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c20_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 67);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 67);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 67);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c20_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 68);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  68);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c20_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 69);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 69);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 69);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c20_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 70);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 70);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c20_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 71);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 71);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c20_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 72);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 72);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c20_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 73);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 73);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c20_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 74);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.refblas.xgemm"), "name", "name", 74);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c20_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c20_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs74), "lhs", "lhs",
                  74);
  sf_mex_destroy(&c20_rhs64);
  sf_mex_destroy(&c20_lhs64);
  sf_mex_destroy(&c20_rhs65);
  sf_mex_destroy(&c20_lhs65);
  sf_mex_destroy(&c20_rhs66);
  sf_mex_destroy(&c20_lhs66);
  sf_mex_destroy(&c20_rhs67);
  sf_mex_destroy(&c20_lhs67);
  sf_mex_destroy(&c20_rhs68);
  sf_mex_destroy(&c20_lhs68);
  sf_mex_destroy(&c20_rhs69);
  sf_mex_destroy(&c20_lhs69);
  sf_mex_destroy(&c20_rhs70);
  sf_mex_destroy(&c20_lhs70);
  sf_mex_destroy(&c20_rhs71);
  sf_mex_destroy(&c20_lhs71);
  sf_mex_destroy(&c20_rhs72);
  sf_mex_destroy(&c20_lhs72);
  sf_mex_destroy(&c20_rhs73);
  sf_mex_destroy(&c20_lhs73);
  sf_mex_destroy(&c20_rhs74);
  sf_mex_destroy(&c20_lhs74);
}

static real_T c20_eml_xnrm2(SFc20_Model_01InstanceStruct *chartInstance, real_T
  c20_x[3])
{
  real_T c20_y;
  real_T c20_scale;
  int32_T c20_k;
  int32_T c20_b_k;
  real_T c20_b_x;
  real_T c20_c_x;
  real_T c20_absxk;
  real_T c20_t;
  (void)chartInstance;
  c20_y = 0.0;
  c20_scale = 2.2250738585072014E-308;
  for (c20_k = 1; c20_k < 4; c20_k++) {
    c20_b_k = c20_k;
    c20_b_x = c20_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c20_b_k), 1, 3, 1, 0) - 1];
    c20_c_x = c20_b_x;
    c20_absxk = muDoubleScalarAbs(c20_c_x);
    if (c20_absxk > c20_scale) {
      c20_t = c20_scale / c20_absxk;
      c20_y = 1.0 + c20_y * c20_t * c20_t;
      c20_scale = c20_absxk;
    } else {
      c20_t = c20_absxk / c20_scale;
      c20_y += c20_t * c20_t;
    }
  }

  return c20_scale * muDoubleScalarSqrt(c20_y);
}

static void c20_eml_scalar_eg(SFc20_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c20_eml_error(SFc20_Model_01InstanceStruct *chartInstance)
{
  int32_T c20_i74;
  static char_T c20_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c20_u[30];
  const mxArray *c20_y = NULL;
  int32_T c20_i75;
  static char_T c20_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c20_b_u[4];
  const mxArray *c20_b_y = NULL;
  (void)chartInstance;
  for (c20_i74 = 0; c20_i74 < 30; c20_i74++) {
    c20_u[c20_i74] = c20_cv0[c20_i74];
  }

  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 10, 0U, 1U, 0U, 2, 1, 30),
                false);
  for (c20_i75 = 0; c20_i75 < 4; c20_i75++) {
    c20_b_u[c20_i75] = c20_cv1[c20_i75];
  }

  c20_b_y = NULL;
  sf_mex_assign(&c20_b_y, sf_mex_create("y", c20_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c20_y, 14, c20_b_y));
}

static void c20_b_eml_scalar_eg(SFc20_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c20_eml_xgemm(SFc20_Model_01InstanceStruct *chartInstance, real_T
  c20_A[16], real_T c20_B[4], real_T c20_C[4], real_T c20_b_C[4])
{
  int32_T c20_i76;
  int32_T c20_i77;
  real_T c20_b_A[16];
  int32_T c20_i78;
  real_T c20_b_B[4];
  for (c20_i76 = 0; c20_i76 < 4; c20_i76++) {
    c20_b_C[c20_i76] = c20_C[c20_i76];
  }

  for (c20_i77 = 0; c20_i77 < 16; c20_i77++) {
    c20_b_A[c20_i77] = c20_A[c20_i77];
  }

  for (c20_i78 = 0; c20_i78 < 4; c20_i78++) {
    c20_b_B[c20_i78] = c20_B[c20_i78];
  }

  c20_b_eml_xgemm(chartInstance, c20_b_A, c20_b_B, c20_b_C);
}

static const mxArray *c20_f_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  int32_T c20_u;
  const mxArray *c20_y = NULL;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  c20_u = *(int32_T *)c20_inData;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", &c20_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, false);
  return c20_mxArrayOutData;
}

static int32_T c20_g_emlrt_marshallIn(SFc20_Model_01InstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId)
{
  int32_T c20_y;
  int32_T c20_i79;
  (void)chartInstance;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), &c20_i79, 1, 6, 0U, 0, 0U, 0);
  c20_y = c20_i79;
  sf_mex_destroy(&c20_u);
  return c20_y;
}

static void c20_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_b_sfEvent;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  int32_T c20_y;
  SFc20_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc20_Model_01InstanceStruct *)chartInstanceVoid;
  c20_b_sfEvent = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_y = c20_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_b_sfEvent),
    &c20_thisId);
  sf_mex_destroy(&c20_b_sfEvent);
  *(int32_T *)c20_outData = c20_y;
  sf_mex_destroy(&c20_mxArrayInData);
}

static uint8_T c20_h_emlrt_marshallIn(SFc20_Model_01InstanceStruct
  *chartInstance, const mxArray *c20_b_is_active_c20_Model_01, const char_T
  *c20_identifier)
{
  uint8_T c20_y;
  emlrtMsgIdentifier c20_thisId;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_y = c20_i_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c20_b_is_active_c20_Model_01), &c20_thisId);
  sf_mex_destroy(&c20_b_is_active_c20_Model_01);
  return c20_y;
}

static uint8_T c20_i_emlrt_marshallIn(SFc20_Model_01InstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId)
{
  uint8_T c20_y;
  uint8_T c20_u0;
  (void)chartInstance;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), &c20_u0, 1, 3, 0U, 0, 0U, 0);
  c20_y = c20_u0;
  sf_mex_destroy(&c20_u);
  return c20_y;
}

static void c20_b_eml_xgemm(SFc20_Model_01InstanceStruct *chartInstance, real_T
  c20_A[16], real_T c20_B[4], real_T c20_C[4])
{
  int32_T c20_i80;
  int32_T c20_i81;
  int32_T c20_i82;
  (void)chartInstance;
  for (c20_i80 = 0; c20_i80 < 4; c20_i80++) {
    c20_C[c20_i80] = 0.0;
    c20_i81 = 0;
    for (c20_i82 = 0; c20_i82 < 4; c20_i82++) {
      c20_C[c20_i80] += c20_A[c20_i81 + c20_i80] * c20_B[c20_i82];
      c20_i81 += 4;
    }
  }
}

static void init_dsm_address_info(SFc20_Model_01InstanceStruct *chartInstance)
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

void sf_c20_Model_01_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(556751772U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1712430319U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2361977879U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1164671196U);
}

mxArray *sf_c20_Model_01_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("FymDh5Efx9OfAaYAG4DJPH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
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
      pr[1] = (double)(1);
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

mxArray *sf_c20_Model_01_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c20_Model_01_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c20_Model_01(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"ita\",},{M[8],M[0],T\"is_active_c20_Model_01\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c20_Model_01_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc20_Model_01InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc20_Model_01InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_01MachineNumber_,
           20,
           1,
           1,
           0,
           3,
           0,
           0,
           0,
           0,
           2,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation(_Model_01MachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_Model_01MachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _Model_01MachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"del_ita_v");
          _SFD_SET_DATA_PROPS(1,2,0,1,"ita");
          _SFD_SET_DATA_PROPS(2,1,1,0,"nominal_ita");
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
        _SFD_CV_INIT_EML(0,1,1,1,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,323);
        _SFD_CV_INIT_EML_IF(0,1,0,66,90,-1,149);
        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,1,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"fn_CrossTensor",59,-1,457);

        {
          static int caseStart[] = { 423, 254, 339 };

          static int caseExprEnd[] = { 432, 260, 345 };

          _SFD_CV_INIT_SCRIPT_SWITCH(0,0,232,245,453,3,&(caseStart[0]),
            &(caseExprEnd[0]));
        }

        _SFD_CV_INIT_SCRIPT(1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(1,0,"fn_VectorToSkewSymmetricTensor",0,-1,433);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c20_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c20_sf_marshallOut,(MexInFcnForType)
            c20_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c20_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T (*c20_del_ita_v)[3];
          real_T (*c20_ita)[4];
          real_T (*c20_nominal_ita)[4];
          c20_nominal_ita = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S,
            1);
          c20_ita = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
          c20_del_ita_v = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S,
            0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c20_del_ita_v);
          _SFD_SET_DATA_VALUE_PTR(1U, *c20_ita);
          _SFD_SET_DATA_VALUE_PTR(2U, *c20_nominal_ita);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _Model_01MachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "VcgptJofnyDTRtHEmw9adE";
}

static void sf_opaque_initialize_c20_Model_01(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc20_Model_01InstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c20_Model_01((SFc20_Model_01InstanceStruct*)
    chartInstanceVar);
  initialize_c20_Model_01((SFc20_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c20_Model_01(void *chartInstanceVar)
{
  enable_c20_Model_01((SFc20_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c20_Model_01(void *chartInstanceVar)
{
  disable_c20_Model_01((SFc20_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c20_Model_01(void *chartInstanceVar)
{
  sf_gateway_c20_Model_01((SFc20_Model_01InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c20_Model_01(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c20_Model_01((SFc20_Model_01InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c20_Model_01();/* state var info */
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

extern void sf_internal_set_sim_state_c20_Model_01(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c20_Model_01();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c20_Model_01((SFc20_Model_01InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c20_Model_01(SimStruct* S)
{
  return sf_internal_get_sim_state_c20_Model_01(S);
}

static void sf_opaque_set_sim_state_c20_Model_01(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c20_Model_01(S, st);
}

static void sf_opaque_terminate_c20_Model_01(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc20_Model_01InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_01_optimization_info();
    }

    finalize_c20_Model_01((SFc20_Model_01InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc20_Model_01((SFc20_Model_01InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c20_Model_01(SimStruct *S)
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
    initialize_params_c20_Model_01((SFc20_Model_01InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c20_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_01_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      20);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,20,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,20,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,20);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,20,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,20,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,20);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3278941852U));
  ssSetChecksum1(S,(2357191565U));
  ssSetChecksum2(S,(1822098290U));
  ssSetChecksum3(S,(2279141086U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c20_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c20_Model_01(SimStruct *S)
{
  SFc20_Model_01InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc20_Model_01InstanceStruct *)utMalloc(sizeof
    (SFc20_Model_01InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc20_Model_01InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c20_Model_01;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c20_Model_01;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c20_Model_01;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c20_Model_01;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c20_Model_01;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c20_Model_01;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c20_Model_01;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c20_Model_01;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c20_Model_01;
  chartInstance->chartInfo.mdlStart = mdlStart_c20_Model_01;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c20_Model_01;
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

void c20_Model_01_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c20_Model_01(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c20_Model_01(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c20_Model_01(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c20_Model_01_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

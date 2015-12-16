/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_01_sfun.h"
#include "c8_Model_01.h"
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
static const char * c8_debug_family_names[11] = { "v_n", "q_n", "q_om",
  "Tensor_om", "Tensor_n", "nargin", "nargout", "om", "n", "q", "q_dot" };

static const char * c8_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

static const char * c8_c_debug_family_names[8] = { "q_v", "q_0", "cross_q_v",
  "nargin", "nargout", "q", "flag", "CrossTensor" };

/* Function Declarations */
static void initialize_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance);
static void initialize_params_c8_Model_01(SFc8_Model_01InstanceStruct
  *chartInstance);
static void enable_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance);
static void disable_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance);
static void c8_update_debugger_state_c8_Model_01(SFc8_Model_01InstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c8_Model_01(SFc8_Model_01InstanceStruct
  *chartInstance);
static void set_sim_state_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_st);
static void finalize_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance);
static void sf_gateway_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance);
static void initSimStructsc8_Model_01(SFc8_Model_01InstanceStruct *chartInstance);
static void c8_fn_CrossTensor(SFc8_Model_01InstanceStruct *chartInstance, real_T
  c8_q[4], real_T c8_flag, real_T c8_CrossTensor[16]);
static void init_script_number_translation(uint32_T c8_machineNumber, uint32_T
  c8_chartNumber, uint32_T c8_instanceNumber);
static const mxArray *c8_sf_marshallOut(void *chartInstanceVoid, void *c8_inData);
static void c8_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_q_dot, const char_T *c8_identifier, real_T c8_y[4]);
static void c8_b_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId, real_T c8_y[4]);
static void c8_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData);
static const mxArray *c8_b_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData);
static const mxArray *c8_c_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData);
static real_T c8_c_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId);
static void c8_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData);
static const mxArray *c8_d_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData);
static void c8_d_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId, real_T c8_y[16]);
static void c8_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData);
static void c8_e_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId, real_T c8_y[3]);
static void c8_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData);
static const mxArray *c8_e_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData);
static void c8_f_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId, real_T c8_y[9]);
static void c8_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData);
static void c8_info_helper(const mxArray **c8_info);
static const mxArray *c8_emlrt_marshallOut(const char * c8_u);
static const mxArray *c8_b_emlrt_marshallOut(const uint32_T c8_u);
static void c8_eml_scalar_eg(SFc8_Model_01InstanceStruct *chartInstance);
static void c8_eml_xgemm(SFc8_Model_01InstanceStruct *chartInstance, real_T
  c8_A[16], real_T c8_B[4], real_T c8_C[4], real_T c8_b_C[4]);
static const mxArray *c8_f_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData);
static int32_T c8_g_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId);
static void c8_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData);
static uint8_T c8_h_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_b_is_active_c8_Model_01, const char_T *c8_identifier);
static uint8_T c8_i_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId);
static void c8_b_eml_xgemm(SFc8_Model_01InstanceStruct *chartInstance, real_T
  c8_A[16], real_T c8_B[4], real_T c8_C[4]);
static void init_dsm_address_info(SFc8_Model_01InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance)
{
  chartInstance->c8_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c8_is_active_c8_Model_01 = 0U;
}

static void initialize_params_c8_Model_01(SFc8_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c8_update_debugger_state_c8_Model_01(SFc8_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c8_Model_01(SFc8_Model_01InstanceStruct
  *chartInstance)
{
  const mxArray *c8_st;
  const mxArray *c8_y = NULL;
  int32_T c8_i0;
  real_T c8_u[4];
  const mxArray *c8_b_y = NULL;
  uint8_T c8_hoistedGlobal;
  uint8_T c8_b_u;
  const mxArray *c8_c_y = NULL;
  real_T (*c8_q_dot)[4];
  c8_q_dot = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c8_st = NULL;
  c8_st = NULL;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_createcellmatrix(2, 1), false);
  for (c8_i0 = 0; c8_i0 < 4; c8_i0++) {
    c8_u[c8_i0] = (*c8_q_dot)[c8_i0];
  }

  c8_b_y = NULL;
  sf_mex_assign(&c8_b_y, sf_mex_create("y", c8_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_setcell(c8_y, 0, c8_b_y);
  c8_hoistedGlobal = chartInstance->c8_is_active_c8_Model_01;
  c8_b_u = c8_hoistedGlobal;
  c8_c_y = NULL;
  sf_mex_assign(&c8_c_y, sf_mex_create("y", &c8_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c8_y, 1, c8_c_y);
  sf_mex_assign(&c8_st, c8_y, false);
  return c8_st;
}

static void set_sim_state_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_st)
{
  const mxArray *c8_u;
  real_T c8_dv0[4];
  int32_T c8_i1;
  real_T (*c8_q_dot)[4];
  c8_q_dot = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c8_doneDoubleBufferReInit = true;
  c8_u = sf_mex_dup(c8_st);
  c8_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c8_u, 0)),
                      "q_dot", c8_dv0);
  for (c8_i1 = 0; c8_i1 < 4; c8_i1++) {
    (*c8_q_dot)[c8_i1] = c8_dv0[c8_i1];
  }

  chartInstance->c8_is_active_c8_Model_01 = c8_h_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c8_u, 1)), "is_active_c8_Model_01");
  sf_mex_destroy(&c8_u);
  c8_update_debugger_state_c8_Model_01(chartInstance);
  sf_mex_destroy(&c8_st);
}

static void finalize_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c8_Model_01(SFc8_Model_01InstanceStruct *chartInstance)
{
  int32_T c8_i2;
  real_T c8_hoistedGlobal;
  int32_T c8_i3;
  real_T c8_om[3];
  real_T c8_n;
  int32_T c8_i4;
  real_T c8_q[4];
  uint32_T c8_debug_family_var_map[11];
  real_T c8_v_n[3];
  real_T c8_q_n[4];
  real_T c8_q_om[4];
  real_T c8_Tensor_om[16];
  real_T c8_Tensor_n[16];
  real_T c8_nargin = 3.0;
  real_T c8_nargout = 1.0;
  real_T c8_q_dot[4];
  int32_T c8_i5;
  int32_T c8_i6;
  int32_T c8_i7;
  real_T c8_b_q_om[4];
  real_T c8_dv1[16];
  int32_T c8_i8;
  int32_T c8_i9;
  real_T c8_b_q_n[4];
  real_T c8_dv2[16];
  int32_T c8_i10;
  int32_T c8_i11;
  real_T c8_b[16];
  int32_T c8_i12;
  int32_T c8_i13;
  real_T c8_b_b[4];
  int32_T c8_i14;
  int32_T c8_i15;
  int32_T c8_i16;
  real_T c8_dv3[16];
  int32_T c8_i17;
  real_T c8_dv4[4];
  int32_T c8_i18;
  real_T c8_dv5[16];
  int32_T c8_i19;
  real_T c8_dv6[4];
  int32_T c8_i20;
  int32_T c8_i21;
  int32_T c8_i22;
  real_T *c8_b_n;
  real_T (*c8_b_q_dot)[4];
  real_T (*c8_b_q)[4];
  real_T (*c8_b_om)[3];
  c8_b_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 2);
  c8_b_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c8_b_q_dot = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c8_b_om = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 7U, chartInstance->c8_sfEvent);
  for (c8_i2 = 0; c8_i2 < 3; c8_i2++) {
    _SFD_DATA_RANGE_CHECK((*c8_b_om)[c8_i2], 0U);
  }

  chartInstance->c8_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 7U, chartInstance->c8_sfEvent);
  c8_hoistedGlobal = *c8_b_n;
  for (c8_i3 = 0; c8_i3 < 3; c8_i3++) {
    c8_om[c8_i3] = (*c8_b_om)[c8_i3];
  }

  c8_n = c8_hoistedGlobal;
  for (c8_i4 = 0; c8_i4 < 4; c8_i4++) {
    c8_q[c8_i4] = (*c8_b_q)[c8_i4];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 11U, 11U, c8_debug_family_names,
    c8_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_v_n, 0U, c8_c_sf_marshallOut,
    c8_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_q_n, 1U, c8_sf_marshallOut,
    c8_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_q_om, 2U, c8_sf_marshallOut,
    c8_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_Tensor_om, 3U, c8_d_sf_marshallOut,
    c8_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_Tensor_n, 4U, c8_d_sf_marshallOut,
    c8_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_nargin, 5U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_nargout, 6U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c8_om, 7U, c8_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_n, 8U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c8_q, 9U, c8_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_q_dot, 10U, c8_sf_marshallOut,
    c8_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 2);
  c8_v_n[0] = 0.0;
  c8_v_n[1] = 0.0;
  c8_v_n[2] = c8_n;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 3);
  for (c8_i5 = 0; c8_i5 < 3; c8_i5++) {
    c8_q_n[c8_i5] = c8_v_n[c8_i5];
  }

  c8_q_n[3] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 4);
  for (c8_i6 = 0; c8_i6 < 3; c8_i6++) {
    c8_q_om[c8_i6] = c8_om[c8_i6];
  }

  c8_q_om[3] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 5);
  for (c8_i7 = 0; c8_i7 < 4; c8_i7++) {
    c8_b_q_om[c8_i7] = c8_q_om[c8_i7];
  }

  c8_fn_CrossTensor(chartInstance, c8_b_q_om, 0.0, c8_dv1);
  for (c8_i8 = 0; c8_i8 < 16; c8_i8++) {
    c8_Tensor_om[c8_i8] = c8_dv1[c8_i8];
  }

  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 6);
  for (c8_i9 = 0; c8_i9 < 4; c8_i9++) {
    c8_b_q_n[c8_i9] = c8_q_n[c8_i9];
  }

  c8_fn_CrossTensor(chartInstance, c8_b_q_n, 1.0, c8_dv2);
  for (c8_i10 = 0; c8_i10 < 16; c8_i10++) {
    c8_Tensor_n[c8_i10] = c8_dv2[c8_i10];
  }

  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 7);
  for (c8_i11 = 0; c8_i11 < 16; c8_i11++) {
    c8_b[c8_i11] = c8_Tensor_om[c8_i11] - c8_Tensor_n[c8_i11];
  }

  for (c8_i12 = 0; c8_i12 < 16; c8_i12++) {
    c8_b[c8_i12] *= 0.5;
  }

  for (c8_i13 = 0; c8_i13 < 4; c8_i13++) {
    c8_b_b[c8_i13] = c8_q[c8_i13];
  }

  c8_eml_scalar_eg(chartInstance);
  c8_eml_scalar_eg(chartInstance);
  for (c8_i14 = 0; c8_i14 < 4; c8_i14++) {
    c8_q_dot[c8_i14] = 0.0;
  }

  for (c8_i15 = 0; c8_i15 < 4; c8_i15++) {
    c8_q_dot[c8_i15] = 0.0;
  }

  for (c8_i16 = 0; c8_i16 < 16; c8_i16++) {
    c8_dv3[c8_i16] = c8_b[c8_i16];
  }

  for (c8_i17 = 0; c8_i17 < 4; c8_i17++) {
    c8_dv4[c8_i17] = c8_b_b[c8_i17];
  }

  for (c8_i18 = 0; c8_i18 < 16; c8_i18++) {
    c8_dv5[c8_i18] = c8_dv3[c8_i18];
  }

  for (c8_i19 = 0; c8_i19 < 4; c8_i19++) {
    c8_dv6[c8_i19] = c8_dv4[c8_i19];
  }

  c8_b_eml_xgemm(chartInstance, c8_dv5, c8_dv6, c8_q_dot);
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, -7);
  _SFD_SYMBOL_SCOPE_POP();
  for (c8_i20 = 0; c8_i20 < 4; c8_i20++) {
    (*c8_b_q_dot)[c8_i20] = c8_q_dot[c8_i20];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 7U, chartInstance->c8_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_01MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c8_i21 = 0; c8_i21 < 4; c8_i21++) {
    _SFD_DATA_RANGE_CHECK((*c8_b_q_dot)[c8_i21], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c8_b_n, 2U);
  for (c8_i22 = 0; c8_i22 < 4; c8_i22++) {
    _SFD_DATA_RANGE_CHECK((*c8_b_q)[c8_i22], 3U);
  }
}

static void initSimStructsc8_Model_01(SFc8_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c8_fn_CrossTensor(SFc8_Model_01InstanceStruct *chartInstance, real_T
  c8_q[4], real_T c8_flag, real_T c8_CrossTensor[16])
{
  uint32_T c8_debug_family_var_map[8];
  real_T c8_q_v[3];
  real_T c8_q_0;
  real_T c8_cross_q_v[9];
  real_T c8_nargin = 2.0;
  real_T c8_nargout = 1.0;
  int32_T c8_i23;
  int32_T c8_i24;
  int32_T c8_i25;
  real_T c8_v[3];
  uint32_T c8_b_debug_family_var_map[4];
  real_T c8_b_nargin = 1.0;
  real_T c8_b_nargout = 1.0;
  int32_T c8_i26;
  real_T c8_a;
  int32_T c8_i27;
  static real_T c8_b[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c8_y[9];
  int32_T c8_i28;
  int32_T c8_i29;
  int32_T c8_i30;
  int32_T c8_i31;
  int32_T c8_i32;
  int32_T c8_i33;
  int32_T c8_i34;
  real_T c8_b_a;
  int32_T c8_i35;
  int32_T c8_i36;
  int32_T c8_i37;
  int32_T c8_i38;
  int32_T c8_i39;
  int32_T c8_i40;
  int32_T c8_i41;
  int32_T c8_i42;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c8_c_debug_family_names,
    c8_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_q_v, 0U, c8_c_sf_marshallOut,
    c8_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_q_0, 1U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_cross_q_v, 2U, c8_e_sf_marshallOut,
    c8_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_nargin, 3U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_nargout, 4U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_q, 5U, c8_sf_marshallOut,
    c8_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_flag, 6U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_CrossTensor, 7U, c8_d_sf_marshallOut,
    c8_c_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c8_sfEvent, 4);
  for (c8_i23 = 0; c8_i23 < 16; c8_i23++) {
    c8_CrossTensor[c8_i23] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c8_sfEvent, 5);
  for (c8_i24 = 0; c8_i24 < 3; c8_i24++) {
    c8_q_v[c8_i24] = c8_q[c8_i24];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c8_sfEvent, 6);
  c8_q_0 = c8_q[3];
  _SFD_SCRIPT_CALL(0U, chartInstance->c8_sfEvent, 7);
  for (c8_i25 = 0; c8_i25 < 3; c8_i25++) {
    c8_v[c8_i25] = c8_q_v[c8_i25];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c8_b_debug_family_names,
    c8_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_b_nargin, 0U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_b_nargout, 1U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_v, 2U, c8_c_sf_marshallOut,
    c8_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_cross_q_v, 3U, c8_e_sf_marshallOut,
    c8_e_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 2);
  for (c8_i26 = 0; c8_i26 < 9; c8_i26++) {
    c8_cross_q_v[c8_i26] = 0.0;
  }

  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 3);
  c8_cross_q_v[0] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 4);
  c8_cross_q_v[4] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 5);
  c8_cross_q_v[8] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 6);
  c8_cross_q_v[3] = -c8_v[2];
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 7);
  c8_cross_q_v[6] = c8_v[1];
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 8);
  c8_cross_q_v[7] = -c8_v[0];
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 9);
  c8_cross_q_v[1] = c8_v[2];
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 10);
  c8_cross_q_v[2] = -c8_v[1];
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, 11);
  c8_cross_q_v[5] = c8_v[0];
  _SFD_SCRIPT_CALL(1U, chartInstance->c8_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_SCRIPT_CALL(0U, chartInstance->c8_sfEvent, 8);
  switch ((int32_T)_SFD_INTEGER_CHECK("flag", c8_flag)) {
   case 0:
    CV_SCRIPT_SWITCH(0, 0, 1);
    _SFD_SCRIPT_CALL(0U, chartInstance->c8_sfEvent, 10);
    c8_a = c8_q_0;
    for (c8_i27 = 0; c8_i27 < 9; c8_i27++) {
      c8_y[c8_i27] = c8_a * c8_b[c8_i27];
    }

    c8_i28 = 0;
    c8_i29 = 0;
    for (c8_i30 = 0; c8_i30 < 3; c8_i30++) {
      for (c8_i31 = 0; c8_i31 < 3; c8_i31++) {
        c8_CrossTensor[c8_i31 + c8_i28] = -c8_cross_q_v[c8_i31 + c8_i29] +
          c8_y[c8_i31 + c8_i29];
      }

      c8_i28 += 4;
      c8_i29 += 3;
    }

    for (c8_i32 = 0; c8_i32 < 3; c8_i32++) {
      c8_CrossTensor[c8_i32 + 12] = c8_q_v[c8_i32];
    }

    c8_i33 = 0;
    for (c8_i34 = 0; c8_i34 < 3; c8_i34++) {
      c8_CrossTensor[c8_i33 + 3] = -c8_q_v[c8_i34];
      c8_i33 += 4;
    }

    c8_CrossTensor[15] = c8_q_0;
    break;

   case 1:
    CV_SCRIPT_SWITCH(0, 0, 2);
    _SFD_SCRIPT_CALL(0U, chartInstance->c8_sfEvent, 12);
    c8_b_a = c8_q_0;
    for (c8_i35 = 0; c8_i35 < 9; c8_i35++) {
      c8_y[c8_i35] = c8_b_a * c8_b[c8_i35];
    }

    c8_i36 = 0;
    c8_i37 = 0;
    for (c8_i38 = 0; c8_i38 < 3; c8_i38++) {
      for (c8_i39 = 0; c8_i39 < 3; c8_i39++) {
        c8_CrossTensor[c8_i39 + c8_i36] = c8_cross_q_v[c8_i39 + c8_i37] +
          c8_y[c8_i39 + c8_i37];
      }

      c8_i36 += 4;
      c8_i37 += 3;
    }

    for (c8_i40 = 0; c8_i40 < 3; c8_i40++) {
      c8_CrossTensor[c8_i40 + 12] = c8_q_v[c8_i40];
    }

    c8_i41 = 0;
    for (c8_i42 = 0; c8_i42 < 3; c8_i42++) {
      c8_CrossTensor[c8_i41 + 3] = -c8_q_v[c8_i42];
      c8_i41 += 4;
    }

    c8_CrossTensor[15] = c8_q_0;
    break;

   default:
    CV_SCRIPT_SWITCH(0, 0, 0);
    break;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c8_sfEvent, -12);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c8_machineNumber, uint32_T
  c8_chartNumber, uint32_T c8_instanceNumber)
{
  (void)c8_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c8_chartNumber, c8_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg-2\\Documents\\MATLAB\\Model_01\\fn_CrossTensor.m"));
  _SFD_SCRIPT_TRANSLATION(c8_chartNumber, c8_instanceNumber, 1U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg-2\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c8_sf_marshallOut(void *chartInstanceVoid, void *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  int32_T c8_i43;
  real_T c8_b_inData[4];
  int32_T c8_i44;
  real_T c8_u[4];
  const mxArray *c8_y = NULL;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  for (c8_i43 = 0; c8_i43 < 4; c8_i43++) {
    c8_b_inData[c8_i43] = (*(real_T (*)[4])c8_inData)[c8_i43];
  }

  for (c8_i44 = 0; c8_i44 < 4; c8_i44++) {
    c8_u[c8_i44] = c8_b_inData[c8_i44];
  }

  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", c8_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static void c8_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_q_dot, const char_T *c8_identifier, real_T c8_y[4])
{
  emlrtMsgIdentifier c8_thisId;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_q_dot), &c8_thisId, c8_y);
  sf_mex_destroy(&c8_q_dot);
}

static void c8_b_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId, real_T c8_y[4])
{
  real_T c8_dv7[4];
  int32_T c8_i45;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), c8_dv7, 1, 0, 0U, 1, 0U, 1, 4);
  for (c8_i45 = 0; c8_i45 < 4; c8_i45++) {
    c8_y[c8_i45] = c8_dv7[c8_i45];
  }

  sf_mex_destroy(&c8_u);
}

static void c8_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData)
{
  const mxArray *c8_q_dot;
  const char_T *c8_identifier;
  emlrtMsgIdentifier c8_thisId;
  real_T c8_y[4];
  int32_T c8_i46;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_q_dot = sf_mex_dup(c8_mxArrayInData);
  c8_identifier = c8_varName;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_q_dot), &c8_thisId, c8_y);
  sf_mex_destroy(&c8_q_dot);
  for (c8_i46 = 0; c8_i46 < 4; c8_i46++) {
    (*(real_T (*)[4])c8_outData)[c8_i46] = c8_y[c8_i46];
  }

  sf_mex_destroy(&c8_mxArrayInData);
}

static const mxArray *c8_b_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  real_T c8_u;
  const mxArray *c8_y = NULL;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  c8_u = *(real_T *)c8_inData;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", &c8_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static const mxArray *c8_c_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  int32_T c8_i47;
  real_T c8_b_inData[3];
  int32_T c8_i48;
  real_T c8_u[3];
  const mxArray *c8_y = NULL;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  for (c8_i47 = 0; c8_i47 < 3; c8_i47++) {
    c8_b_inData[c8_i47] = (*(real_T (*)[3])c8_inData)[c8_i47];
  }

  for (c8_i48 = 0; c8_i48 < 3; c8_i48++) {
    c8_u[c8_i48] = c8_b_inData[c8_i48];
  }

  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", c8_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static real_T c8_c_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId)
{
  real_T c8_y;
  real_T c8_d0;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), &c8_d0, 1, 0, 0U, 0, 0U, 0);
  c8_y = c8_d0;
  sf_mex_destroy(&c8_u);
  return c8_y;
}

static void c8_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData)
{
  const mxArray *c8_nargout;
  const char_T *c8_identifier;
  emlrtMsgIdentifier c8_thisId;
  real_T c8_y;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_nargout = sf_mex_dup(c8_mxArrayInData);
  c8_identifier = c8_varName;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_y = c8_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_nargout), &c8_thisId);
  sf_mex_destroy(&c8_nargout);
  *(real_T *)c8_outData = c8_y;
  sf_mex_destroy(&c8_mxArrayInData);
}

static const mxArray *c8_d_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  int32_T c8_i49;
  int32_T c8_i50;
  int32_T c8_i51;
  real_T c8_b_inData[16];
  int32_T c8_i52;
  int32_T c8_i53;
  int32_T c8_i54;
  real_T c8_u[16];
  const mxArray *c8_y = NULL;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  c8_i49 = 0;
  for (c8_i50 = 0; c8_i50 < 4; c8_i50++) {
    for (c8_i51 = 0; c8_i51 < 4; c8_i51++) {
      c8_b_inData[c8_i51 + c8_i49] = (*(real_T (*)[16])c8_inData)[c8_i51 +
        c8_i49];
    }

    c8_i49 += 4;
  }

  c8_i52 = 0;
  for (c8_i53 = 0; c8_i53 < 4; c8_i53++) {
    for (c8_i54 = 0; c8_i54 < 4; c8_i54++) {
      c8_u[c8_i54 + c8_i52] = c8_b_inData[c8_i54 + c8_i52];
    }

    c8_i52 += 4;
  }

  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", c8_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static void c8_d_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId, real_T c8_y[16])
{
  real_T c8_dv8[16];
  int32_T c8_i55;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), c8_dv8, 1, 0, 0U, 1, 0U, 2, 4, 4);
  for (c8_i55 = 0; c8_i55 < 16; c8_i55++) {
    c8_y[c8_i55] = c8_dv8[c8_i55];
  }

  sf_mex_destroy(&c8_u);
}

static void c8_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData)
{
  const mxArray *c8_Tensor_n;
  const char_T *c8_identifier;
  emlrtMsgIdentifier c8_thisId;
  real_T c8_y[16];
  int32_T c8_i56;
  int32_T c8_i57;
  int32_T c8_i58;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_Tensor_n = sf_mex_dup(c8_mxArrayInData);
  c8_identifier = c8_varName;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_Tensor_n), &c8_thisId, c8_y);
  sf_mex_destroy(&c8_Tensor_n);
  c8_i56 = 0;
  for (c8_i57 = 0; c8_i57 < 4; c8_i57++) {
    for (c8_i58 = 0; c8_i58 < 4; c8_i58++) {
      (*(real_T (*)[16])c8_outData)[c8_i58 + c8_i56] = c8_y[c8_i58 + c8_i56];
    }

    c8_i56 += 4;
  }

  sf_mex_destroy(&c8_mxArrayInData);
}

static void c8_e_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId, real_T c8_y[3])
{
  real_T c8_dv9[3];
  int32_T c8_i59;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), c8_dv9, 1, 0, 0U, 1, 0U, 1, 3);
  for (c8_i59 = 0; c8_i59 < 3; c8_i59++) {
    c8_y[c8_i59] = c8_dv9[c8_i59];
  }

  sf_mex_destroy(&c8_u);
}

static void c8_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData)
{
  const mxArray *c8_v_n;
  const char_T *c8_identifier;
  emlrtMsgIdentifier c8_thisId;
  real_T c8_y[3];
  int32_T c8_i60;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_v_n = sf_mex_dup(c8_mxArrayInData);
  c8_identifier = c8_varName;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_v_n), &c8_thisId, c8_y);
  sf_mex_destroy(&c8_v_n);
  for (c8_i60 = 0; c8_i60 < 3; c8_i60++) {
    (*(real_T (*)[3])c8_outData)[c8_i60] = c8_y[c8_i60];
  }

  sf_mex_destroy(&c8_mxArrayInData);
}

static const mxArray *c8_e_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  int32_T c8_i61;
  int32_T c8_i62;
  int32_T c8_i63;
  real_T c8_b_inData[9];
  int32_T c8_i64;
  int32_T c8_i65;
  int32_T c8_i66;
  real_T c8_u[9];
  const mxArray *c8_y = NULL;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  c8_i61 = 0;
  for (c8_i62 = 0; c8_i62 < 3; c8_i62++) {
    for (c8_i63 = 0; c8_i63 < 3; c8_i63++) {
      c8_b_inData[c8_i63 + c8_i61] = (*(real_T (*)[9])c8_inData)[c8_i63 + c8_i61];
    }

    c8_i61 += 3;
  }

  c8_i64 = 0;
  for (c8_i65 = 0; c8_i65 < 3; c8_i65++) {
    for (c8_i66 = 0; c8_i66 < 3; c8_i66++) {
      c8_u[c8_i66 + c8_i64] = c8_b_inData[c8_i66 + c8_i64];
    }

    c8_i64 += 3;
  }

  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", c8_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static void c8_f_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId, real_T c8_y[9])
{
  real_T c8_dv10[9];
  int32_T c8_i67;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), c8_dv10, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c8_i67 = 0; c8_i67 < 9; c8_i67++) {
    c8_y[c8_i67] = c8_dv10[c8_i67];
  }

  sf_mex_destroy(&c8_u);
}

static void c8_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData)
{
  const mxArray *c8_SkewSymmetricTensor;
  const char_T *c8_identifier;
  emlrtMsgIdentifier c8_thisId;
  real_T c8_y[9];
  int32_T c8_i68;
  int32_T c8_i69;
  int32_T c8_i70;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_SkewSymmetricTensor = sf_mex_dup(c8_mxArrayInData);
  c8_identifier = c8_varName;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_SkewSymmetricTensor),
                        &c8_thisId, c8_y);
  sf_mex_destroy(&c8_SkewSymmetricTensor);
  c8_i68 = 0;
  for (c8_i69 = 0; c8_i69 < 3; c8_i69++) {
    for (c8_i70 = 0; c8_i70 < 3; c8_i70++) {
      (*(real_T (*)[9])c8_outData)[c8_i70 + c8_i68] = c8_y[c8_i70 + c8_i68];
    }

    c8_i68 += 3;
  }

  sf_mex_destroy(&c8_mxArrayInData);
}

const mxArray *sf_c8_Model_01_get_eml_resolved_functions_info(void)
{
  const mxArray *c8_nameCaptureInfo = NULL;
  c8_nameCaptureInfo = NULL;
  sf_mex_assign(&c8_nameCaptureInfo, sf_mex_createstruct("structure", 2, 41, 1),
                false);
  c8_info_helper(&c8_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c8_nameCaptureInfo);
  return c8_nameCaptureInfo;
}

static void c8_info_helper(const mxArray **c8_info)
{
  const mxArray *c8_rhs0 = NULL;
  const mxArray *c8_lhs0 = NULL;
  const mxArray *c8_rhs1 = NULL;
  const mxArray *c8_lhs1 = NULL;
  const mxArray *c8_rhs2 = NULL;
  const mxArray *c8_lhs2 = NULL;
  const mxArray *c8_rhs3 = NULL;
  const mxArray *c8_lhs3 = NULL;
  const mxArray *c8_rhs4 = NULL;
  const mxArray *c8_lhs4 = NULL;
  const mxArray *c8_rhs5 = NULL;
  const mxArray *c8_lhs5 = NULL;
  const mxArray *c8_rhs6 = NULL;
  const mxArray *c8_lhs6 = NULL;
  const mxArray *c8_rhs7 = NULL;
  const mxArray *c8_lhs7 = NULL;
  const mxArray *c8_rhs8 = NULL;
  const mxArray *c8_lhs8 = NULL;
  const mxArray *c8_rhs9 = NULL;
  const mxArray *c8_lhs9 = NULL;
  const mxArray *c8_rhs10 = NULL;
  const mxArray *c8_lhs10 = NULL;
  const mxArray *c8_rhs11 = NULL;
  const mxArray *c8_lhs11 = NULL;
  const mxArray *c8_rhs12 = NULL;
  const mxArray *c8_lhs12 = NULL;
  const mxArray *c8_rhs13 = NULL;
  const mxArray *c8_lhs13 = NULL;
  const mxArray *c8_rhs14 = NULL;
  const mxArray *c8_lhs14 = NULL;
  const mxArray *c8_rhs15 = NULL;
  const mxArray *c8_lhs15 = NULL;
  const mxArray *c8_rhs16 = NULL;
  const mxArray *c8_lhs16 = NULL;
  const mxArray *c8_rhs17 = NULL;
  const mxArray *c8_lhs17 = NULL;
  const mxArray *c8_rhs18 = NULL;
  const mxArray *c8_lhs18 = NULL;
  const mxArray *c8_rhs19 = NULL;
  const mxArray *c8_lhs19 = NULL;
  const mxArray *c8_rhs20 = NULL;
  const mxArray *c8_lhs20 = NULL;
  const mxArray *c8_rhs21 = NULL;
  const mxArray *c8_lhs21 = NULL;
  const mxArray *c8_rhs22 = NULL;
  const mxArray *c8_lhs22 = NULL;
  const mxArray *c8_rhs23 = NULL;
  const mxArray *c8_lhs23 = NULL;
  const mxArray *c8_rhs24 = NULL;
  const mxArray *c8_lhs24 = NULL;
  const mxArray *c8_rhs25 = NULL;
  const mxArray *c8_lhs25 = NULL;
  const mxArray *c8_rhs26 = NULL;
  const mxArray *c8_lhs26 = NULL;
  const mxArray *c8_rhs27 = NULL;
  const mxArray *c8_lhs27 = NULL;
  const mxArray *c8_rhs28 = NULL;
  const mxArray *c8_lhs28 = NULL;
  const mxArray *c8_rhs29 = NULL;
  const mxArray *c8_lhs29 = NULL;
  const mxArray *c8_rhs30 = NULL;
  const mxArray *c8_lhs30 = NULL;
  const mxArray *c8_rhs31 = NULL;
  const mxArray *c8_lhs31 = NULL;
  const mxArray *c8_rhs32 = NULL;
  const mxArray *c8_lhs32 = NULL;
  const mxArray *c8_rhs33 = NULL;
  const mxArray *c8_lhs33 = NULL;
  const mxArray *c8_rhs34 = NULL;
  const mxArray *c8_lhs34 = NULL;
  const mxArray *c8_rhs35 = NULL;
  const mxArray *c8_lhs35 = NULL;
  const mxArray *c8_rhs36 = NULL;
  const mxArray *c8_lhs36 = NULL;
  const mxArray *c8_rhs37 = NULL;
  const mxArray *c8_lhs37 = NULL;
  const mxArray *c8_rhs38 = NULL;
  const mxArray *c8_lhs38 = NULL;
  const mxArray *c8_rhs39 = NULL;
  const mxArray *c8_lhs39 = NULL;
  const mxArray *c8_rhs40 = NULL;
  const mxArray *c8_lhs40 = NULL;
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("fn_CrossTensor"), "name",
                  "name", 0);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1450227411U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c8_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "context", "context", 1);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 1);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1450040424U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c8_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "context", "context", 2);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eye"), "name", "name", 2);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c8_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 3);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c8_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 4);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c8_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 5);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("isinf"), "name", "name", 5);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c8_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 6);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c8_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 7);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 7);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c8_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 8);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("intmax"), "name", "name", 8);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c8_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 9);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c8_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 10);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("intmin"), "name", "name", 10);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c8_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 11);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c8_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 12);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 12);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c8_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 13);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 13);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c8_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 14);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 14);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c8_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 15);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("intmin"), "name", "name", 15);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c8_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 16);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c8_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("intmax"), "name", "name", 17);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c8_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 18);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c8_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 19);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("intmax"), "name", "name", 19);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c8_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "context", "context", 20);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 20);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c8_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 21);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 21);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c8_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "context", "context", 22);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("mrdivide"), "name", "name", 22);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c8_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 23);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 23);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c8_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 24);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("rdivide"), "name", "name", 24);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c8_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 25);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c8_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 26);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c8_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_div"), "name", "name", 27);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c8_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 28);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c8_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "context", "context", 29);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 29);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c8_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 30);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c8_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 31);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 31);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 31);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c8_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 32);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c8_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 33);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  33);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c8_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 34);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 34);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c8_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 35);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c8_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 36);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 36);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c8_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 37);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 37);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c8_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 38);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 38);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c8_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 39);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 39);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c8_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 40);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 40);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c8_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs40), "lhs", "lhs",
                  40);
  sf_mex_destroy(&c8_rhs0);
  sf_mex_destroy(&c8_lhs0);
  sf_mex_destroy(&c8_rhs1);
  sf_mex_destroy(&c8_lhs1);
  sf_mex_destroy(&c8_rhs2);
  sf_mex_destroy(&c8_lhs2);
  sf_mex_destroy(&c8_rhs3);
  sf_mex_destroy(&c8_lhs3);
  sf_mex_destroy(&c8_rhs4);
  sf_mex_destroy(&c8_lhs4);
  sf_mex_destroy(&c8_rhs5);
  sf_mex_destroy(&c8_lhs5);
  sf_mex_destroy(&c8_rhs6);
  sf_mex_destroy(&c8_lhs6);
  sf_mex_destroy(&c8_rhs7);
  sf_mex_destroy(&c8_lhs7);
  sf_mex_destroy(&c8_rhs8);
  sf_mex_destroy(&c8_lhs8);
  sf_mex_destroy(&c8_rhs9);
  sf_mex_destroy(&c8_lhs9);
  sf_mex_destroy(&c8_rhs10);
  sf_mex_destroy(&c8_lhs10);
  sf_mex_destroy(&c8_rhs11);
  sf_mex_destroy(&c8_lhs11);
  sf_mex_destroy(&c8_rhs12);
  sf_mex_destroy(&c8_lhs12);
  sf_mex_destroy(&c8_rhs13);
  sf_mex_destroy(&c8_lhs13);
  sf_mex_destroy(&c8_rhs14);
  sf_mex_destroy(&c8_lhs14);
  sf_mex_destroy(&c8_rhs15);
  sf_mex_destroy(&c8_lhs15);
  sf_mex_destroy(&c8_rhs16);
  sf_mex_destroy(&c8_lhs16);
  sf_mex_destroy(&c8_rhs17);
  sf_mex_destroy(&c8_lhs17);
  sf_mex_destroy(&c8_rhs18);
  sf_mex_destroy(&c8_lhs18);
  sf_mex_destroy(&c8_rhs19);
  sf_mex_destroy(&c8_lhs19);
  sf_mex_destroy(&c8_rhs20);
  sf_mex_destroy(&c8_lhs20);
  sf_mex_destroy(&c8_rhs21);
  sf_mex_destroy(&c8_lhs21);
  sf_mex_destroy(&c8_rhs22);
  sf_mex_destroy(&c8_lhs22);
  sf_mex_destroy(&c8_rhs23);
  sf_mex_destroy(&c8_lhs23);
  sf_mex_destroy(&c8_rhs24);
  sf_mex_destroy(&c8_lhs24);
  sf_mex_destroy(&c8_rhs25);
  sf_mex_destroy(&c8_lhs25);
  sf_mex_destroy(&c8_rhs26);
  sf_mex_destroy(&c8_lhs26);
  sf_mex_destroy(&c8_rhs27);
  sf_mex_destroy(&c8_lhs27);
  sf_mex_destroy(&c8_rhs28);
  sf_mex_destroy(&c8_lhs28);
  sf_mex_destroy(&c8_rhs29);
  sf_mex_destroy(&c8_lhs29);
  sf_mex_destroy(&c8_rhs30);
  sf_mex_destroy(&c8_lhs30);
  sf_mex_destroy(&c8_rhs31);
  sf_mex_destroy(&c8_lhs31);
  sf_mex_destroy(&c8_rhs32);
  sf_mex_destroy(&c8_lhs32);
  sf_mex_destroy(&c8_rhs33);
  sf_mex_destroy(&c8_lhs33);
  sf_mex_destroy(&c8_rhs34);
  sf_mex_destroy(&c8_lhs34);
  sf_mex_destroy(&c8_rhs35);
  sf_mex_destroy(&c8_lhs35);
  sf_mex_destroy(&c8_rhs36);
  sf_mex_destroy(&c8_lhs36);
  sf_mex_destroy(&c8_rhs37);
  sf_mex_destroy(&c8_lhs37);
  sf_mex_destroy(&c8_rhs38);
  sf_mex_destroy(&c8_lhs38);
  sf_mex_destroy(&c8_rhs39);
  sf_mex_destroy(&c8_lhs39);
  sf_mex_destroy(&c8_rhs40);
  sf_mex_destroy(&c8_lhs40);
}

static const mxArray *c8_emlrt_marshallOut(const char * c8_u)
{
  const mxArray *c8_y = NULL;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", c8_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c8_u)), false);
  return c8_y;
}

static const mxArray *c8_b_emlrt_marshallOut(const uint32_T c8_u)
{
  const mxArray *c8_y = NULL;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", &c8_u, 7, 0U, 0U, 0U, 0), false);
  return c8_y;
}

static void c8_eml_scalar_eg(SFc8_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c8_eml_xgemm(SFc8_Model_01InstanceStruct *chartInstance, real_T
  c8_A[16], real_T c8_B[4], real_T c8_C[4], real_T c8_b_C[4])
{
  int32_T c8_i71;
  int32_T c8_i72;
  real_T c8_b_A[16];
  int32_T c8_i73;
  real_T c8_b_B[4];
  for (c8_i71 = 0; c8_i71 < 4; c8_i71++) {
    c8_b_C[c8_i71] = c8_C[c8_i71];
  }

  for (c8_i72 = 0; c8_i72 < 16; c8_i72++) {
    c8_b_A[c8_i72] = c8_A[c8_i72];
  }

  for (c8_i73 = 0; c8_i73 < 4; c8_i73++) {
    c8_b_B[c8_i73] = c8_B[c8_i73];
  }

  c8_b_eml_xgemm(chartInstance, c8_b_A, c8_b_B, c8_b_C);
}

static const mxArray *c8_f_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  int32_T c8_u;
  const mxArray *c8_y = NULL;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  c8_u = *(int32_T *)c8_inData;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", &c8_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static int32_T c8_g_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId)
{
  int32_T c8_y;
  int32_T c8_i74;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), &c8_i74, 1, 6, 0U, 0, 0U, 0);
  c8_y = c8_i74;
  sf_mex_destroy(&c8_u);
  return c8_y;
}

static void c8_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData)
{
  const mxArray *c8_b_sfEvent;
  const char_T *c8_identifier;
  emlrtMsgIdentifier c8_thisId;
  int32_T c8_y;
  SFc8_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc8_Model_01InstanceStruct *)chartInstanceVoid;
  c8_b_sfEvent = sf_mex_dup(c8_mxArrayInData);
  c8_identifier = c8_varName;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_y = c8_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_b_sfEvent),
    &c8_thisId);
  sf_mex_destroy(&c8_b_sfEvent);
  *(int32_T *)c8_outData = c8_y;
  sf_mex_destroy(&c8_mxArrayInData);
}

static uint8_T c8_h_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_b_is_active_c8_Model_01, const char_T *c8_identifier)
{
  uint8_T c8_y;
  emlrtMsgIdentifier c8_thisId;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_y = c8_i_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c8_b_is_active_c8_Model_01), &c8_thisId);
  sf_mex_destroy(&c8_b_is_active_c8_Model_01);
  return c8_y;
}

static uint8_T c8_i_emlrt_marshallIn(SFc8_Model_01InstanceStruct *chartInstance,
  const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId)
{
  uint8_T c8_y;
  uint8_T c8_u0;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), &c8_u0, 1, 3, 0U, 0, 0U, 0);
  c8_y = c8_u0;
  sf_mex_destroy(&c8_u);
  return c8_y;
}

static void c8_b_eml_xgemm(SFc8_Model_01InstanceStruct *chartInstance, real_T
  c8_A[16], real_T c8_B[4], real_T c8_C[4])
{
  int32_T c8_i75;
  int32_T c8_i76;
  int32_T c8_i77;
  (void)chartInstance;
  for (c8_i75 = 0; c8_i75 < 4; c8_i75++) {
    c8_C[c8_i75] = 0.0;
    c8_i76 = 0;
    for (c8_i77 = 0; c8_i77 < 4; c8_i77++) {
      c8_C[c8_i75] += c8_A[c8_i76 + c8_i75] * c8_B[c8_i77];
      c8_i76 += 4;
    }
  }
}

static void init_dsm_address_info(SFc8_Model_01InstanceStruct *chartInstance)
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

void sf_c8_Model_01_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3113265144U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3536496537U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3085676960U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3810426170U);
}

mxArray *sf_c8_Model_01_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("StOLOaIbyTPtcQv8FGc3PC");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

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
      pr[0] = (double)(1);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
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

mxArray *sf_c8_Model_01_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c8_Model_01_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c8_Model_01(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"q_dot\",},{M[8],M[0],T\"is_active_c8_Model_01\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c8_Model_01_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc8_Model_01InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc8_Model_01InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_01MachineNumber_,
           8,
           1,
           1,
           0,
           4,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"om");
          _SFD_SET_DATA_PROPS(1,2,0,1,"q_dot");
          _SFD_SET_DATA_PROPS(2,1,1,0,"n");
          _SFD_SET_DATA_PROPS(3,1,1,0,"q");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,225);
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
            1.0,0,0,(MexFcnForType)c8_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c8_sf_marshallOut,(MexInFcnForType)
            c8_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c8_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c8_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c8_n;
          real_T (*c8_om)[3];
          real_T (*c8_q_dot)[4];
          real_T (*c8_q)[4];
          c8_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 2);
          c8_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c8_q_dot = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
          c8_om = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c8_om);
          _SFD_SET_DATA_VALUE_PTR(1U, *c8_q_dot);
          _SFD_SET_DATA_VALUE_PTR(2U, c8_n);
          _SFD_SET_DATA_VALUE_PTR(3U, *c8_q);
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
  return "aWEZbcySXBVyXqdTVPBGpG";
}

static void sf_opaque_initialize_c8_Model_01(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc8_Model_01InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c8_Model_01((SFc8_Model_01InstanceStruct*) chartInstanceVar);
  initialize_c8_Model_01((SFc8_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c8_Model_01(void *chartInstanceVar)
{
  enable_c8_Model_01((SFc8_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c8_Model_01(void *chartInstanceVar)
{
  disable_c8_Model_01((SFc8_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c8_Model_01(void *chartInstanceVar)
{
  sf_gateway_c8_Model_01((SFc8_Model_01InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c8_Model_01(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c8_Model_01((SFc8_Model_01InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c8_Model_01();/* state var info */
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

extern void sf_internal_set_sim_state_c8_Model_01(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c8_Model_01();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c8_Model_01((SFc8_Model_01InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c8_Model_01(SimStruct* S)
{
  return sf_internal_get_sim_state_c8_Model_01(S);
}

static void sf_opaque_set_sim_state_c8_Model_01(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c8_Model_01(S, st);
}

static void sf_opaque_terminate_c8_Model_01(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc8_Model_01InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_01_optimization_info();
    }

    finalize_c8_Model_01((SFc8_Model_01InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc8_Model_01((SFc8_Model_01InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c8_Model_01(SimStruct *S)
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
    initialize_params_c8_Model_01((SFc8_Model_01InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c8_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_01_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,8);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,8,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,8,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,8);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,8,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,8,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 3; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,8);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1631386465U));
  ssSetChecksum1(S,(3803898727U));
  ssSetChecksum2(S,(766063664U));
  ssSetChecksum3(S,(1303678544U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c8_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c8_Model_01(SimStruct *S)
{
  SFc8_Model_01InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc8_Model_01InstanceStruct *)utMalloc(sizeof
    (SFc8_Model_01InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc8_Model_01InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c8_Model_01;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c8_Model_01;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c8_Model_01;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c8_Model_01;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c8_Model_01;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c8_Model_01;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c8_Model_01;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c8_Model_01;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c8_Model_01;
  chartInstance->chartInfo.mdlStart = mdlStart_c8_Model_01;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c8_Model_01;
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

void c8_Model_01_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c8_Model_01(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c8_Model_01(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c8_Model_01(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c8_Model_01_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

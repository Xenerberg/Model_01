/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_01_sfun.h"
#include "c9_Model_01.h"
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
static const char * c9_debug_family_names[15] = { "parameter_gravitation", "re",
  "v_n", "Tensor_n", "Term_1", "Term_2", "Term_3", "Term_4", "nargin", "nargout",
  "n", "r", "e_force", "r_dot", "r_ddot" };

static const char * c9_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

/* Function Declarations */
static void initialize_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance);
static void initialize_params_c9_Model_01(SFc9_Model_01InstanceStruct
  *chartInstance);
static void enable_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance);
static void disable_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance);
static void c9_update_debugger_state_c9_Model_01(SFc9_Model_01InstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c9_Model_01(SFc9_Model_01InstanceStruct
  *chartInstance);
static void set_sim_state_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_st);
static void finalize_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance);
static void sf_gateway_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance);
static void c9_chartstep_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance);
static void initSimStructsc9_Model_01(SFc9_Model_01InstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c9_machineNumber, uint32_T
  c9_chartNumber, uint32_T c9_instanceNumber);
static const mxArray *c9_sf_marshallOut(void *chartInstanceVoid, void *c9_inData);
static void c9_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_r_ddot, const char_T *c9_identifier, real_T c9_y[3]);
static void c9_b_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, real_T c9_y[3]);
static void c9_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData);
static const mxArray *c9_b_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData);
static real_T c9_c_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId);
static void c9_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData);
static const mxArray *c9_c_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData);
static void c9_d_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, real_T c9_y[9]);
static void c9_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData);
static void c9_info_helper(const mxArray **c9_info);
static const mxArray *c9_emlrt_marshallOut(const char * c9_u);
static const mxArray *c9_b_emlrt_marshallOut(const uint32_T c9_u);
static real_T c9_mpower(SFc9_Model_01InstanceStruct *chartInstance, real_T c9_a);
static void c9_eml_scalar_eg(SFc9_Model_01InstanceStruct *chartInstance);
static real_T c9_b_mpower(SFc9_Model_01InstanceStruct *chartInstance, real_T
  c9_a);
static void c9_eml_error(SFc9_Model_01InstanceStruct *chartInstance);
static void c9_b_eml_scalar_eg(SFc9_Model_01InstanceStruct *chartInstance);
static void c9_eml_xgemm(SFc9_Model_01InstanceStruct *chartInstance, real_T
  c9_A[9], real_T c9_B[3], real_T c9_C[3], real_T c9_b_C[3]);
static real_T c9_eml_xnrm2(SFc9_Model_01InstanceStruct *chartInstance, real_T
  c9_x[3]);
static const mxArray *c9_d_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData);
static int32_T c9_e_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId);
static void c9_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData);
static uint8_T c9_f_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_b_is_active_c9_Model_01, const char_T *c9_identifier);
static uint8_T c9_g_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId);
static void c9_b_eml_xgemm(SFc9_Model_01InstanceStruct *chartInstance, real_T
  c9_A[9], real_T c9_B[3], real_T c9_C[3]);
static void init_dsm_address_info(SFc9_Model_01InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance)
{
  chartInstance->c9_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c9_is_active_c9_Model_01 = 0U;
}

static void initialize_params_c9_Model_01(SFc9_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c9_update_debugger_state_c9_Model_01(SFc9_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c9_Model_01(SFc9_Model_01InstanceStruct
  *chartInstance)
{
  const mxArray *c9_st;
  const mxArray *c9_y = NULL;
  int32_T c9_i0;
  real_T c9_u[3];
  const mxArray *c9_b_y = NULL;
  uint8_T c9_hoistedGlobal;
  uint8_T c9_b_u;
  const mxArray *c9_c_y = NULL;
  real_T (*c9_r_ddot)[3];
  c9_r_ddot = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c9_st = NULL;
  c9_st = NULL;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_createcellmatrix(2, 1), false);
  for (c9_i0 = 0; c9_i0 < 3; c9_i0++) {
    c9_u[c9_i0] = (*c9_r_ddot)[c9_i0];
  }

  c9_b_y = NULL;
  sf_mex_assign(&c9_b_y, sf_mex_create("y", c9_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c9_y, 0, c9_b_y);
  c9_hoistedGlobal = chartInstance->c9_is_active_c9_Model_01;
  c9_b_u = c9_hoistedGlobal;
  c9_c_y = NULL;
  sf_mex_assign(&c9_c_y, sf_mex_create("y", &c9_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c9_y, 1, c9_c_y);
  sf_mex_assign(&c9_st, c9_y, false);
  return c9_st;
}

static void set_sim_state_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_st)
{
  const mxArray *c9_u;
  real_T c9_dv0[3];
  int32_T c9_i1;
  real_T (*c9_r_ddot)[3];
  c9_r_ddot = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c9_doneDoubleBufferReInit = true;
  c9_u = sf_mex_dup(c9_st);
  c9_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c9_u, 0)),
                      "r_ddot", c9_dv0);
  for (c9_i1 = 0; c9_i1 < 3; c9_i1++) {
    (*c9_r_ddot)[c9_i1] = c9_dv0[c9_i1];
  }

  chartInstance->c9_is_active_c9_Model_01 = c9_f_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c9_u, 1)), "is_active_c9_Model_01");
  sf_mex_destroy(&c9_u);
  c9_update_debugger_state_c9_Model_01(chartInstance);
  sf_mex_destroy(&c9_st);
}

static void finalize_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance)
{
  int32_T c9_i2;
  int32_T c9_i3;
  int32_T c9_i4;
  int32_T c9_i5;
  real_T *c9_n;
  real_T (*c9_r_dot)[3];
  real_T (*c9_e_force)[3];
  real_T (*c9_r)[3];
  real_T (*c9_r_ddot)[3];
  c9_r_dot = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 3);
  c9_e_force = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 2);
  c9_r = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c9_r_ddot = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c9_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 8U, chartInstance->c9_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c9_n, 0U);
  chartInstance->c9_sfEvent = CALL_EVENT;
  c9_chartstep_c9_Model_01(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_01MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c9_i2 = 0; c9_i2 < 3; c9_i2++) {
    _SFD_DATA_RANGE_CHECK((*c9_r_ddot)[c9_i2], 1U);
  }

  for (c9_i3 = 0; c9_i3 < 3; c9_i3++) {
    _SFD_DATA_RANGE_CHECK((*c9_r)[c9_i3], 2U);
  }

  for (c9_i4 = 0; c9_i4 < 3; c9_i4++) {
    _SFD_DATA_RANGE_CHECK((*c9_e_force)[c9_i4], 3U);
  }

  for (c9_i5 = 0; c9_i5 < 3; c9_i5++) {
    _SFD_DATA_RANGE_CHECK((*c9_r_dot)[c9_i5], 4U);
  }
}

static void c9_chartstep_c9_Model_01(SFc9_Model_01InstanceStruct *chartInstance)
{
  real_T c9_hoistedGlobal;
  real_T c9_n;
  int32_T c9_i6;
  real_T c9_r[3];
  int32_T c9_i7;
  real_T c9_e_force[3];
  int32_T c9_i8;
  real_T c9_r_dot[3];
  uint32_T c9_debug_family_var_map[15];
  real_T c9_parameter_gravitation;
  real_T c9_re[3];
  real_T c9_v_n[3];
  real_T c9_Tensor_n[9];
  real_T c9_Term_1[3];
  real_T c9_Term_2[3];
  real_T c9_Term_3[3];
  real_T c9_Term_4[3];
  real_T c9_nargin = 4.0;
  real_T c9_nargout = 1.0;
  real_T c9_r_ddot[3];
  real_T c9_B;
  real_T c9_y;
  real_T c9_b_y;
  real_T c9_c_y;
  real_T c9_d_y;
  real_T c9_d0;
  int32_T c9_i9;
  real_T c9_v[3];
  uint32_T c9_b_debug_family_var_map[4];
  real_T c9_b_nargin = 1.0;
  real_T c9_b_nargout = 1.0;
  int32_T c9_i10;
  int32_T c9_i11;
  real_T c9_b[9];
  int32_T c9_i12;
  int32_T c9_i13;
  real_T c9_b_b[3];
  int32_T c9_i14;
  int32_T c9_i15;
  int32_T c9_i16;
  real_T c9_dv1[9];
  int32_T c9_i17;
  real_T c9_dv2[3];
  int32_T c9_i18;
  real_T c9_dv3[9];
  int32_T c9_i19;
  real_T c9_dv4[3];
  int32_T c9_i20;
  int32_T c9_i21;
  int32_T c9_i22;
  real_T c9_e_y[3];
  int32_T c9_i23;
  real_T c9_c_b[9];
  int32_T c9_i24;
  real_T c9_d_b[3];
  int32_T c9_i25;
  int32_T c9_i26;
  int32_T c9_i27;
  int32_T c9_i28;
  real_T c9_dv5[9];
  int32_T c9_i29;
  real_T c9_dv6[3];
  int32_T c9_i30;
  real_T c9_dv7[9];
  int32_T c9_i31;
  real_T c9_dv8[3];
  int32_T c9_i32;
  int32_T c9_i33;
  int32_T c9_i34;
  int32_T c9_i35;
  real_T c9_f_y[3];
  real_T c9_g_y;
  real_T c9_a;
  real_T c9_b_a;
  real_T c9_c_a;
  real_T c9_ak;
  real_T c9_d_a;
  real_T c9_ar;
  real_T c9_c;
  real_T c9_b_B;
  real_T c9_h_y;
  real_T c9_i_y;
  real_T c9_j_y;
  int32_T c9_i36;
  real_T c9_e_a;
  int32_T c9_i37;
  int32_T c9_i38;
  int32_T c9_i39;
  int32_T c9_i40;
  real_T *c9_b_n;
  real_T (*c9_b_r_ddot)[3];
  real_T (*c9_b_r_dot)[3];
  real_T (*c9_b_e_force)[3];
  real_T (*c9_b_r)[3];
  c9_b_r_dot = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 3);
  c9_b_e_force = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 2);
  c9_b_r = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c9_b_r_ddot = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c9_b_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 8U, chartInstance->c9_sfEvent);
  c9_hoistedGlobal = *c9_b_n;
  c9_n = c9_hoistedGlobal;
  for (c9_i6 = 0; c9_i6 < 3; c9_i6++) {
    c9_r[c9_i6] = (*c9_b_r)[c9_i6];
  }

  for (c9_i7 = 0; c9_i7 < 3; c9_i7++) {
    c9_e_force[c9_i7] = (*c9_b_e_force)[c9_i7];
  }

  for (c9_i8 = 0; c9_i8 < 3; c9_i8++) {
    c9_r_dot[c9_i8] = (*c9_b_r_dot)[c9_i8];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 15U, 15U, c9_debug_family_names,
    c9_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c9_parameter_gravitation, 0U, c9_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_re, 1U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_v_n, 2U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_Tensor_n, 3U, c9_c_sf_marshallOut,
    c9_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_Term_1, 4U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_Term_2, 5U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_Term_3, 6U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_Term_4, 7U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_nargin, 8U, c9_b_sf_marshallOut,
    c9_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_nargout, 9U, c9_b_sf_marshallOut,
    c9_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c9_n, 10U, c9_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c9_r, 11U, c9_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c9_e_force, 12U, c9_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c9_r_dot, 13U, c9_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_r_ddot, 14U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 3);
  c9_parameter_gravitation = 3.986005E+14;
  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 4);
  c9_B = c9_mpower(chartInstance, c9_n);
  c9_y = c9_B;
  c9_b_y = c9_y;
  c9_c_y = c9_b_y;
  c9_d_y = 3.986005E+14 / c9_c_y;
  c9_d0 = c9_b_mpower(chartInstance, c9_d_y);
  c9_re[0] = c9_d0;
  c9_re[1] = 0.0;
  c9_re[2] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 5);
  c9_v_n[0] = 0.0;
  c9_v_n[1] = 0.0;
  c9_v_n[2] = c9_n;
  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 6);
  for (c9_i9 = 0; c9_i9 < 3; c9_i9++) {
    c9_v[c9_i9] = c9_v_n[c9_i9];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c9_b_debug_family_names,
    c9_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_b_nargin, 0U, c9_b_sf_marshallOut,
    c9_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_b_nargout, 1U, c9_b_sf_marshallOut,
    c9_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_v, 2U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c9_Tensor_n, 3U, c9_c_sf_marshallOut,
    c9_c_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 2);
  for (c9_i10 = 0; c9_i10 < 9; c9_i10++) {
    c9_Tensor_n[c9_i10] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 3);
  c9_Tensor_n[0] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 4);
  c9_Tensor_n[4] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 5);
  c9_Tensor_n[8] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 6);
  c9_Tensor_n[3] = -c9_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 7);
  c9_Tensor_n[6] = c9_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 8);
  c9_Tensor_n[7] = -c9_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 9);
  c9_Tensor_n[1] = c9_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 10);
  c9_Tensor_n[2] = -c9_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, 11);
  c9_Tensor_n[5] = c9_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c9_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 7);
  for (c9_i11 = 0; c9_i11 < 9; c9_i11++) {
    c9_b[c9_i11] = c9_Tensor_n[c9_i11];
  }

  for (c9_i12 = 0; c9_i12 < 9; c9_i12++) {
    c9_b[c9_i12] *= -2.0;
  }

  for (c9_i13 = 0; c9_i13 < 3; c9_i13++) {
    c9_b_b[c9_i13] = c9_r_dot[c9_i13];
  }

  c9_b_eml_scalar_eg(chartInstance);
  c9_b_eml_scalar_eg(chartInstance);
  for (c9_i14 = 0; c9_i14 < 3; c9_i14++) {
    c9_Term_1[c9_i14] = 0.0;
  }

  for (c9_i15 = 0; c9_i15 < 3; c9_i15++) {
    c9_Term_1[c9_i15] = 0.0;
  }

  for (c9_i16 = 0; c9_i16 < 9; c9_i16++) {
    c9_dv1[c9_i16] = c9_b[c9_i16];
  }

  for (c9_i17 = 0; c9_i17 < 3; c9_i17++) {
    c9_dv2[c9_i17] = c9_b_b[c9_i17];
  }

  for (c9_i18 = 0; c9_i18 < 9; c9_i18++) {
    c9_dv3[c9_i18] = c9_dv1[c9_i18];
  }

  for (c9_i19 = 0; c9_i19 < 3; c9_i19++) {
    c9_dv4[c9_i19] = c9_dv2[c9_i19];
  }

  c9_b_eml_xgemm(chartInstance, c9_dv3, c9_dv4, c9_Term_1);
  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 8);
  for (c9_i20 = 0; c9_i20 < 9; c9_i20++) {
    c9_b[c9_i20] = c9_Tensor_n[c9_i20];
  }

  for (c9_i21 = 0; c9_i21 < 3; c9_i21++) {
    c9_b_b[c9_i21] = c9_r[c9_i21];
  }

  c9_b_eml_scalar_eg(chartInstance);
  c9_b_eml_scalar_eg(chartInstance);
  for (c9_i22 = 0; c9_i22 < 3; c9_i22++) {
    c9_e_y[c9_i22] = 0.0;
  }

  for (c9_i23 = 0; c9_i23 < 9; c9_i23++) {
    c9_c_b[c9_i23] = c9_b[c9_i23];
  }

  for (c9_i24 = 0; c9_i24 < 3; c9_i24++) {
    c9_d_b[c9_i24] = c9_b_b[c9_i24];
  }

  c9_b_eml_xgemm(chartInstance, c9_c_b, c9_d_b, c9_e_y);
  for (c9_i25 = 0; c9_i25 < 9; c9_i25++) {
    c9_b[c9_i25] = c9_Tensor_n[c9_i25];
  }

  c9_b_eml_scalar_eg(chartInstance);
  c9_b_eml_scalar_eg(chartInstance);
  for (c9_i26 = 0; c9_i26 < 3; c9_i26++) {
    c9_Term_2[c9_i26] = 0.0;
  }

  for (c9_i27 = 0; c9_i27 < 3; c9_i27++) {
    c9_Term_2[c9_i27] = 0.0;
  }

  for (c9_i28 = 0; c9_i28 < 9; c9_i28++) {
    c9_dv5[c9_i28] = c9_b[c9_i28];
  }

  for (c9_i29 = 0; c9_i29 < 3; c9_i29++) {
    c9_dv6[c9_i29] = c9_e_y[c9_i29];
  }

  for (c9_i30 = 0; c9_i30 < 9; c9_i30++) {
    c9_dv7[c9_i30] = c9_dv5[c9_i30];
  }

  for (c9_i31 = 0; c9_i31 < 3; c9_i31++) {
    c9_dv8[c9_i31] = c9_dv6[c9_i31];
  }

  c9_b_eml_xgemm(chartInstance, c9_dv7, c9_dv8, c9_Term_2);
  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 9);
  for (c9_i32 = 0; c9_i32 < 3; c9_i32++) {
    c9_b_b[c9_i32] = c9_re[c9_i32] + c9_r[c9_i32];
  }

  for (c9_i33 = 0; c9_i33 < 3; c9_i33++) {
    c9_b_b[c9_i33] *= 3.986005E+14;
  }

  for (c9_i34 = 0; c9_i34 < 3; c9_i34++) {
    c9_e_y[c9_i34] = c9_re[c9_i34] + c9_r[c9_i34];
  }

  for (c9_i35 = 0; c9_i35 < 3; c9_i35++) {
    c9_f_y[c9_i35] = c9_e_y[c9_i35];
  }

  c9_g_y = c9_eml_xnrm2(chartInstance, c9_f_y);
  c9_a = c9_g_y;
  c9_b_a = c9_a;
  c9_c_a = c9_b_a;
  c9_eml_scalar_eg(chartInstance);
  c9_ak = c9_c_a;
  c9_d_a = c9_ak;
  c9_eml_scalar_eg(chartInstance);
  c9_ar = c9_d_a;
  c9_c = muDoubleScalarPower(c9_ar, 3.0);
  c9_b_B = c9_c;
  c9_h_y = c9_b_B;
  c9_i_y = c9_h_y;
  c9_j_y = c9_i_y;
  for (c9_i36 = 0; c9_i36 < 3; c9_i36++) {
    c9_Term_3[c9_i36] = c9_b_b[c9_i36] / c9_j_y;
  }

  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 10);
  c9_e_a = c9_mpower(chartInstance, c9_n);
  for (c9_i37 = 0; c9_i37 < 3; c9_i37++) {
    c9_b_b[c9_i37] = c9_re[c9_i37];
  }

  for (c9_i38 = 0; c9_i38 < 3; c9_i38++) {
    c9_Term_4[c9_i38] = c9_e_a * c9_b_b[c9_i38];
  }

  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 12);
  for (c9_i39 = 0; c9_i39 < 3; c9_i39++) {
    c9_r_ddot[c9_i39] = (((c9_Term_1[c9_i39] - c9_Term_2[c9_i39]) -
                          c9_Term_3[c9_i39]) + c9_Term_4[c9_i39]) +
      c9_e_force[c9_i39];
  }

  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, -12);
  _SFD_SYMBOL_SCOPE_POP();
  for (c9_i40 = 0; c9_i40 < 3; c9_i40++) {
    (*c9_b_r_ddot)[c9_i40] = c9_r_ddot[c9_i40];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 8U, chartInstance->c9_sfEvent);
}

static void initSimStructsc9_Model_01(SFc9_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c9_machineNumber, uint32_T
  c9_chartNumber, uint32_T c9_instanceNumber)
{
  (void)c9_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c9_chartNumber, c9_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c9_sf_marshallOut(void *chartInstanceVoid, void *c9_inData)
{
  const mxArray *c9_mxArrayOutData = NULL;
  int32_T c9_i41;
  real_T c9_b_inData[3];
  int32_T c9_i42;
  real_T c9_u[3];
  const mxArray *c9_y = NULL;
  SFc9_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc9_Model_01InstanceStruct *)chartInstanceVoid;
  c9_mxArrayOutData = NULL;
  for (c9_i41 = 0; c9_i41 < 3; c9_i41++) {
    c9_b_inData[c9_i41] = (*(real_T (*)[3])c9_inData)[c9_i41];
  }

  for (c9_i42 = 0; c9_i42 < 3; c9_i42++) {
    c9_u[c9_i42] = c9_b_inData[c9_i42];
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c9_mxArrayOutData, c9_y, false);
  return c9_mxArrayOutData;
}

static void c9_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_r_ddot, const char_T *c9_identifier, real_T c9_y[3])
{
  emlrtMsgIdentifier c9_thisId;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_r_ddot), &c9_thisId, c9_y);
  sf_mex_destroy(&c9_r_ddot);
}

static void c9_b_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, real_T c9_y[3])
{
  real_T c9_dv9[3];
  int32_T c9_i43;
  (void)chartInstance;
  sf_mex_import(c9_parentId, sf_mex_dup(c9_u), c9_dv9, 1, 0, 0U, 1, 0U, 1, 3);
  for (c9_i43 = 0; c9_i43 < 3; c9_i43++) {
    c9_y[c9_i43] = c9_dv9[c9_i43];
  }

  sf_mex_destroy(&c9_u);
}

static void c9_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData)
{
  const mxArray *c9_r_ddot;
  const char_T *c9_identifier;
  emlrtMsgIdentifier c9_thisId;
  real_T c9_y[3];
  int32_T c9_i44;
  SFc9_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc9_Model_01InstanceStruct *)chartInstanceVoid;
  c9_r_ddot = sf_mex_dup(c9_mxArrayInData);
  c9_identifier = c9_varName;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_r_ddot), &c9_thisId, c9_y);
  sf_mex_destroy(&c9_r_ddot);
  for (c9_i44 = 0; c9_i44 < 3; c9_i44++) {
    (*(real_T (*)[3])c9_outData)[c9_i44] = c9_y[c9_i44];
  }

  sf_mex_destroy(&c9_mxArrayInData);
}

static const mxArray *c9_b_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData)
{
  const mxArray *c9_mxArrayOutData = NULL;
  real_T c9_u;
  const mxArray *c9_y = NULL;
  SFc9_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc9_Model_01InstanceStruct *)chartInstanceVoid;
  c9_mxArrayOutData = NULL;
  c9_u = *(real_T *)c9_inData;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", &c9_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c9_mxArrayOutData, c9_y, false);
  return c9_mxArrayOutData;
}

static real_T c9_c_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId)
{
  real_T c9_y;
  real_T c9_d1;
  (void)chartInstance;
  sf_mex_import(c9_parentId, sf_mex_dup(c9_u), &c9_d1, 1, 0, 0U, 0, 0U, 0);
  c9_y = c9_d1;
  sf_mex_destroy(&c9_u);
  return c9_y;
}

static void c9_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData)
{
  const mxArray *c9_nargout;
  const char_T *c9_identifier;
  emlrtMsgIdentifier c9_thisId;
  real_T c9_y;
  SFc9_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc9_Model_01InstanceStruct *)chartInstanceVoid;
  c9_nargout = sf_mex_dup(c9_mxArrayInData);
  c9_identifier = c9_varName;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_y = c9_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_nargout), &c9_thisId);
  sf_mex_destroy(&c9_nargout);
  *(real_T *)c9_outData = c9_y;
  sf_mex_destroy(&c9_mxArrayInData);
}

static const mxArray *c9_c_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData)
{
  const mxArray *c9_mxArrayOutData = NULL;
  int32_T c9_i45;
  int32_T c9_i46;
  int32_T c9_i47;
  real_T c9_b_inData[9];
  int32_T c9_i48;
  int32_T c9_i49;
  int32_T c9_i50;
  real_T c9_u[9];
  const mxArray *c9_y = NULL;
  SFc9_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc9_Model_01InstanceStruct *)chartInstanceVoid;
  c9_mxArrayOutData = NULL;
  c9_i45 = 0;
  for (c9_i46 = 0; c9_i46 < 3; c9_i46++) {
    for (c9_i47 = 0; c9_i47 < 3; c9_i47++) {
      c9_b_inData[c9_i47 + c9_i45] = (*(real_T (*)[9])c9_inData)[c9_i47 + c9_i45];
    }

    c9_i45 += 3;
  }

  c9_i48 = 0;
  for (c9_i49 = 0; c9_i49 < 3; c9_i49++) {
    for (c9_i50 = 0; c9_i50 < 3; c9_i50++) {
      c9_u[c9_i50 + c9_i48] = c9_b_inData[c9_i50 + c9_i48];
    }

    c9_i48 += 3;
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c9_mxArrayOutData, c9_y, false);
  return c9_mxArrayOutData;
}

static void c9_d_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, real_T c9_y[9])
{
  real_T c9_dv10[9];
  int32_T c9_i51;
  (void)chartInstance;
  sf_mex_import(c9_parentId, sf_mex_dup(c9_u), c9_dv10, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c9_i51 = 0; c9_i51 < 9; c9_i51++) {
    c9_y[c9_i51] = c9_dv10[c9_i51];
  }

  sf_mex_destroy(&c9_u);
}

static void c9_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData)
{
  const mxArray *c9_Tensor_n;
  const char_T *c9_identifier;
  emlrtMsgIdentifier c9_thisId;
  real_T c9_y[9];
  int32_T c9_i52;
  int32_T c9_i53;
  int32_T c9_i54;
  SFc9_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc9_Model_01InstanceStruct *)chartInstanceVoid;
  c9_Tensor_n = sf_mex_dup(c9_mxArrayInData);
  c9_identifier = c9_varName;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_Tensor_n), &c9_thisId, c9_y);
  sf_mex_destroy(&c9_Tensor_n);
  c9_i52 = 0;
  for (c9_i53 = 0; c9_i53 < 3; c9_i53++) {
    for (c9_i54 = 0; c9_i54 < 3; c9_i54++) {
      (*(real_T (*)[9])c9_outData)[c9_i54 + c9_i52] = c9_y[c9_i54 + c9_i52];
    }

    c9_i52 += 3;
  }

  sf_mex_destroy(&c9_mxArrayInData);
}

const mxArray *sf_c9_Model_01_get_eml_resolved_functions_info(void)
{
  const mxArray *c9_nameCaptureInfo = NULL;
  c9_nameCaptureInfo = NULL;
  sf_mex_assign(&c9_nameCaptureInfo, sf_mex_createstruct("structure", 2, 55, 1),
                false);
  c9_info_helper(&c9_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c9_nameCaptureInfo);
  return c9_nameCaptureInfo;
}

static void c9_info_helper(const mxArray **c9_info)
{
  const mxArray *c9_rhs0 = NULL;
  const mxArray *c9_lhs0 = NULL;
  const mxArray *c9_rhs1 = NULL;
  const mxArray *c9_lhs1 = NULL;
  const mxArray *c9_rhs2 = NULL;
  const mxArray *c9_lhs2 = NULL;
  const mxArray *c9_rhs3 = NULL;
  const mxArray *c9_lhs3 = NULL;
  const mxArray *c9_rhs4 = NULL;
  const mxArray *c9_lhs4 = NULL;
  const mxArray *c9_rhs5 = NULL;
  const mxArray *c9_lhs5 = NULL;
  const mxArray *c9_rhs6 = NULL;
  const mxArray *c9_lhs6 = NULL;
  const mxArray *c9_rhs7 = NULL;
  const mxArray *c9_lhs7 = NULL;
  const mxArray *c9_rhs8 = NULL;
  const mxArray *c9_lhs8 = NULL;
  const mxArray *c9_rhs9 = NULL;
  const mxArray *c9_lhs9 = NULL;
  const mxArray *c9_rhs10 = NULL;
  const mxArray *c9_lhs10 = NULL;
  const mxArray *c9_rhs11 = NULL;
  const mxArray *c9_lhs11 = NULL;
  const mxArray *c9_rhs12 = NULL;
  const mxArray *c9_lhs12 = NULL;
  const mxArray *c9_rhs13 = NULL;
  const mxArray *c9_lhs13 = NULL;
  const mxArray *c9_rhs14 = NULL;
  const mxArray *c9_lhs14 = NULL;
  const mxArray *c9_rhs15 = NULL;
  const mxArray *c9_lhs15 = NULL;
  const mxArray *c9_rhs16 = NULL;
  const mxArray *c9_lhs16 = NULL;
  const mxArray *c9_rhs17 = NULL;
  const mxArray *c9_lhs17 = NULL;
  const mxArray *c9_rhs18 = NULL;
  const mxArray *c9_lhs18 = NULL;
  const mxArray *c9_rhs19 = NULL;
  const mxArray *c9_lhs19 = NULL;
  const mxArray *c9_rhs20 = NULL;
  const mxArray *c9_lhs20 = NULL;
  const mxArray *c9_rhs21 = NULL;
  const mxArray *c9_lhs21 = NULL;
  const mxArray *c9_rhs22 = NULL;
  const mxArray *c9_lhs22 = NULL;
  const mxArray *c9_rhs23 = NULL;
  const mxArray *c9_lhs23 = NULL;
  const mxArray *c9_rhs24 = NULL;
  const mxArray *c9_lhs24 = NULL;
  const mxArray *c9_rhs25 = NULL;
  const mxArray *c9_lhs25 = NULL;
  const mxArray *c9_rhs26 = NULL;
  const mxArray *c9_lhs26 = NULL;
  const mxArray *c9_rhs27 = NULL;
  const mxArray *c9_lhs27 = NULL;
  const mxArray *c9_rhs28 = NULL;
  const mxArray *c9_lhs28 = NULL;
  const mxArray *c9_rhs29 = NULL;
  const mxArray *c9_lhs29 = NULL;
  const mxArray *c9_rhs30 = NULL;
  const mxArray *c9_lhs30 = NULL;
  const mxArray *c9_rhs31 = NULL;
  const mxArray *c9_lhs31 = NULL;
  const mxArray *c9_rhs32 = NULL;
  const mxArray *c9_lhs32 = NULL;
  const mxArray *c9_rhs33 = NULL;
  const mxArray *c9_lhs33 = NULL;
  const mxArray *c9_rhs34 = NULL;
  const mxArray *c9_lhs34 = NULL;
  const mxArray *c9_rhs35 = NULL;
  const mxArray *c9_lhs35 = NULL;
  const mxArray *c9_rhs36 = NULL;
  const mxArray *c9_lhs36 = NULL;
  const mxArray *c9_rhs37 = NULL;
  const mxArray *c9_lhs37 = NULL;
  const mxArray *c9_rhs38 = NULL;
  const mxArray *c9_lhs38 = NULL;
  const mxArray *c9_rhs39 = NULL;
  const mxArray *c9_lhs39 = NULL;
  const mxArray *c9_rhs40 = NULL;
  const mxArray *c9_lhs40 = NULL;
  const mxArray *c9_rhs41 = NULL;
  const mxArray *c9_lhs41 = NULL;
  const mxArray *c9_rhs42 = NULL;
  const mxArray *c9_lhs42 = NULL;
  const mxArray *c9_rhs43 = NULL;
  const mxArray *c9_lhs43 = NULL;
  const mxArray *c9_rhs44 = NULL;
  const mxArray *c9_lhs44 = NULL;
  const mxArray *c9_rhs45 = NULL;
  const mxArray *c9_lhs45 = NULL;
  const mxArray *c9_rhs46 = NULL;
  const mxArray *c9_lhs46 = NULL;
  const mxArray *c9_rhs47 = NULL;
  const mxArray *c9_lhs47 = NULL;
  const mxArray *c9_rhs48 = NULL;
  const mxArray *c9_lhs48 = NULL;
  const mxArray *c9_rhs49 = NULL;
  const mxArray *c9_lhs49 = NULL;
  const mxArray *c9_rhs50 = NULL;
  const mxArray *c9_lhs50 = NULL;
  const mxArray *c9_rhs51 = NULL;
  const mxArray *c9_lhs51 = NULL;
  const mxArray *c9_rhs52 = NULL;
  const mxArray *c9_lhs52 = NULL;
  const mxArray *c9_rhs53 = NULL;
  const mxArray *c9_lhs53 = NULL;
  const mxArray *c9_rhs54 = NULL;
  const mxArray *c9_lhs54 = NULL;
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("mpower"), "name", "name", 0);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c9_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c9_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("ismatrix"), "name", "name", 2);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c9_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("power"), "name", "name", 3);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c9_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c9_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 5);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 5);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c9_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 6);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c9_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 7);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 7);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c9_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 8);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 8);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c9_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 9);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("floor"), "name", "name", 9);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c9_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 10);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c9_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 11);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c9_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 12);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c9_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 13);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("mrdivide"), "name", "name", 13);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c9_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 14);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 14);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c9_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 15);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("rdivide"), "name", "name", 15);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c9_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 16);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c9_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 17);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c9_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_div"), "name", "name", 18);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c9_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 19);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c9_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 20);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_error"), "name", "name",
                  20);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c9_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 21);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 21);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1447321639U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c9_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 22);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 22);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c9_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 23);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 23);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c9_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 24);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 24);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c9_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 25);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 25);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c9_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 26);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  26);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c9_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 27);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c9_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 28);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c9_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 29);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 29);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c9_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 30);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 30);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c9_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 31);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 31);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c9_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 32);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 32);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c9_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 33);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 33);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c9_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 34);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("norm"), "name", "name", 34);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 34);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c9_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 35);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 35);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c9_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 36);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 36);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c9_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 37);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  37);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c9_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 38);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c9_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 39);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.blas.xnrm2"),
                  "name", "name", 39);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c9_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 40);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 40);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c9_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p!below_threshold"),
                  "context", "context", 41);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 41);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c9_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 42);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.refblas.xnrm2"),
                  "name", "name", 42);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c9_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 43);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("realmin"), "name", "name", 43);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 43);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c9_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 44);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_realmin"), "name", "name",
                  44);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1307658444U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c9_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 45);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 45);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c9_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 46);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 46);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c9_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 47);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 47);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 47);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c9_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 48);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 48);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 48);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c9_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 49);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 49);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c9_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 50);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmax"), "name", "name", 50);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 50);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c9_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 51);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 51);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c9_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 52);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("abs"), "name", "name", 52);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 52);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c9_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 53);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 53);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c9_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 54);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 54);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c9_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c9_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs54), "lhs", "lhs",
                  54);
  sf_mex_destroy(&c9_rhs0);
  sf_mex_destroy(&c9_lhs0);
  sf_mex_destroy(&c9_rhs1);
  sf_mex_destroy(&c9_lhs1);
  sf_mex_destroy(&c9_rhs2);
  sf_mex_destroy(&c9_lhs2);
  sf_mex_destroy(&c9_rhs3);
  sf_mex_destroy(&c9_lhs3);
  sf_mex_destroy(&c9_rhs4);
  sf_mex_destroy(&c9_lhs4);
  sf_mex_destroy(&c9_rhs5);
  sf_mex_destroy(&c9_lhs5);
  sf_mex_destroy(&c9_rhs6);
  sf_mex_destroy(&c9_lhs6);
  sf_mex_destroy(&c9_rhs7);
  sf_mex_destroy(&c9_lhs7);
  sf_mex_destroy(&c9_rhs8);
  sf_mex_destroy(&c9_lhs8);
  sf_mex_destroy(&c9_rhs9);
  sf_mex_destroy(&c9_lhs9);
  sf_mex_destroy(&c9_rhs10);
  sf_mex_destroy(&c9_lhs10);
  sf_mex_destroy(&c9_rhs11);
  sf_mex_destroy(&c9_lhs11);
  sf_mex_destroy(&c9_rhs12);
  sf_mex_destroy(&c9_lhs12);
  sf_mex_destroy(&c9_rhs13);
  sf_mex_destroy(&c9_lhs13);
  sf_mex_destroy(&c9_rhs14);
  sf_mex_destroy(&c9_lhs14);
  sf_mex_destroy(&c9_rhs15);
  sf_mex_destroy(&c9_lhs15);
  sf_mex_destroy(&c9_rhs16);
  sf_mex_destroy(&c9_lhs16);
  sf_mex_destroy(&c9_rhs17);
  sf_mex_destroy(&c9_lhs17);
  sf_mex_destroy(&c9_rhs18);
  sf_mex_destroy(&c9_lhs18);
  sf_mex_destroy(&c9_rhs19);
  sf_mex_destroy(&c9_lhs19);
  sf_mex_destroy(&c9_rhs20);
  sf_mex_destroy(&c9_lhs20);
  sf_mex_destroy(&c9_rhs21);
  sf_mex_destroy(&c9_lhs21);
  sf_mex_destroy(&c9_rhs22);
  sf_mex_destroy(&c9_lhs22);
  sf_mex_destroy(&c9_rhs23);
  sf_mex_destroy(&c9_lhs23);
  sf_mex_destroy(&c9_rhs24);
  sf_mex_destroy(&c9_lhs24);
  sf_mex_destroy(&c9_rhs25);
  sf_mex_destroy(&c9_lhs25);
  sf_mex_destroy(&c9_rhs26);
  sf_mex_destroy(&c9_lhs26);
  sf_mex_destroy(&c9_rhs27);
  sf_mex_destroy(&c9_lhs27);
  sf_mex_destroy(&c9_rhs28);
  sf_mex_destroy(&c9_lhs28);
  sf_mex_destroy(&c9_rhs29);
  sf_mex_destroy(&c9_lhs29);
  sf_mex_destroy(&c9_rhs30);
  sf_mex_destroy(&c9_lhs30);
  sf_mex_destroy(&c9_rhs31);
  sf_mex_destroy(&c9_lhs31);
  sf_mex_destroy(&c9_rhs32);
  sf_mex_destroy(&c9_lhs32);
  sf_mex_destroy(&c9_rhs33);
  sf_mex_destroy(&c9_lhs33);
  sf_mex_destroy(&c9_rhs34);
  sf_mex_destroy(&c9_lhs34);
  sf_mex_destroy(&c9_rhs35);
  sf_mex_destroy(&c9_lhs35);
  sf_mex_destroy(&c9_rhs36);
  sf_mex_destroy(&c9_lhs36);
  sf_mex_destroy(&c9_rhs37);
  sf_mex_destroy(&c9_lhs37);
  sf_mex_destroy(&c9_rhs38);
  sf_mex_destroy(&c9_lhs38);
  sf_mex_destroy(&c9_rhs39);
  sf_mex_destroy(&c9_lhs39);
  sf_mex_destroy(&c9_rhs40);
  sf_mex_destroy(&c9_lhs40);
  sf_mex_destroy(&c9_rhs41);
  sf_mex_destroy(&c9_lhs41);
  sf_mex_destroy(&c9_rhs42);
  sf_mex_destroy(&c9_lhs42);
  sf_mex_destroy(&c9_rhs43);
  sf_mex_destroy(&c9_lhs43);
  sf_mex_destroy(&c9_rhs44);
  sf_mex_destroy(&c9_lhs44);
  sf_mex_destroy(&c9_rhs45);
  sf_mex_destroy(&c9_lhs45);
  sf_mex_destroy(&c9_rhs46);
  sf_mex_destroy(&c9_lhs46);
  sf_mex_destroy(&c9_rhs47);
  sf_mex_destroy(&c9_lhs47);
  sf_mex_destroy(&c9_rhs48);
  sf_mex_destroy(&c9_lhs48);
  sf_mex_destroy(&c9_rhs49);
  sf_mex_destroy(&c9_lhs49);
  sf_mex_destroy(&c9_rhs50);
  sf_mex_destroy(&c9_lhs50);
  sf_mex_destroy(&c9_rhs51);
  sf_mex_destroy(&c9_lhs51);
  sf_mex_destroy(&c9_rhs52);
  sf_mex_destroy(&c9_lhs52);
  sf_mex_destroy(&c9_rhs53);
  sf_mex_destroy(&c9_lhs53);
  sf_mex_destroy(&c9_rhs54);
  sf_mex_destroy(&c9_lhs54);
}

static const mxArray *c9_emlrt_marshallOut(const char * c9_u)
{
  const mxArray *c9_y = NULL;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c9_u)), false);
  return c9_y;
}

static const mxArray *c9_b_emlrt_marshallOut(const uint32_T c9_u)
{
  const mxArray *c9_y = NULL;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", &c9_u, 7, 0U, 0U, 0U, 0), false);
  return c9_y;
}

static real_T c9_mpower(SFc9_Model_01InstanceStruct *chartInstance, real_T c9_a)
{
  real_T c9_b_a;
  real_T c9_c_a;
  real_T c9_ak;
  real_T c9_d_a;
  c9_b_a = c9_a;
  c9_c_a = c9_b_a;
  c9_eml_scalar_eg(chartInstance);
  c9_ak = c9_c_a;
  c9_d_a = c9_ak;
  c9_eml_scalar_eg(chartInstance);
  return c9_d_a * c9_d_a;
}

static void c9_eml_scalar_eg(SFc9_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c9_b_mpower(SFc9_Model_01InstanceStruct *chartInstance, real_T
  c9_a)
{
  real_T c9_b_a;
  real_T c9_c_a;
  real_T c9_ak;
  real_T c9_d_a;
  real_T c9_ar;
  c9_b_a = c9_a;
  c9_c_a = c9_b_a;
  c9_eml_scalar_eg(chartInstance);
  c9_ak = c9_c_a;
  if (c9_ak < 0.0) {
    c9_eml_error(chartInstance);
  }

  c9_d_a = c9_ak;
  c9_eml_scalar_eg(chartInstance);
  c9_ar = c9_d_a;
  return muDoubleScalarPower(c9_ar, 0.33333333333333331);
}

static void c9_eml_error(SFc9_Model_01InstanceStruct *chartInstance)
{
  int32_T c9_i55;
  static char_T c9_cv0[31] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'p', 'o', 'w', 'e', 'r', '_', 'd', 'o', 'm', 'a', 'i',
    'n', 'E', 'r', 'r', 'o', 'r' };

  char_T c9_u[31];
  const mxArray *c9_y = NULL;
  (void)chartInstance;
  for (c9_i55 = 0; c9_i55 < 31; c9_i55++) {
    c9_u[c9_i55] = c9_cv0[c9_i55];
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 10, 0U, 1U, 0U, 2, 1, 31), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c9_y));
}

static void c9_b_eml_scalar_eg(SFc9_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c9_eml_xgemm(SFc9_Model_01InstanceStruct *chartInstance, real_T
  c9_A[9], real_T c9_B[3], real_T c9_C[3], real_T c9_b_C[3])
{
  int32_T c9_i56;
  int32_T c9_i57;
  real_T c9_b_A[9];
  int32_T c9_i58;
  real_T c9_b_B[3];
  for (c9_i56 = 0; c9_i56 < 3; c9_i56++) {
    c9_b_C[c9_i56] = c9_C[c9_i56];
  }

  for (c9_i57 = 0; c9_i57 < 9; c9_i57++) {
    c9_b_A[c9_i57] = c9_A[c9_i57];
  }

  for (c9_i58 = 0; c9_i58 < 3; c9_i58++) {
    c9_b_B[c9_i58] = c9_B[c9_i58];
  }

  c9_b_eml_xgemm(chartInstance, c9_b_A, c9_b_B, c9_b_C);
}

static real_T c9_eml_xnrm2(SFc9_Model_01InstanceStruct *chartInstance, real_T
  c9_x[3])
{
  real_T c9_y;
  real_T c9_scale;
  int32_T c9_k;
  int32_T c9_b_k;
  real_T c9_b_x;
  real_T c9_c_x;
  real_T c9_absxk;
  real_T c9_t;
  (void)chartInstance;
  c9_y = 0.0;
  c9_scale = 2.2250738585072014E-308;
  for (c9_k = 1; c9_k < 4; c9_k++) {
    c9_b_k = c9_k;
    c9_b_x = c9_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c9_b_k), 1, 3, 1, 0) - 1];
    c9_c_x = c9_b_x;
    c9_absxk = muDoubleScalarAbs(c9_c_x);
    if (c9_absxk > c9_scale) {
      c9_t = c9_scale / c9_absxk;
      c9_y = 1.0 + c9_y * c9_t * c9_t;
      c9_scale = c9_absxk;
    } else {
      c9_t = c9_absxk / c9_scale;
      c9_y += c9_t * c9_t;
    }
  }

  return c9_scale * muDoubleScalarSqrt(c9_y);
}

static const mxArray *c9_d_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData)
{
  const mxArray *c9_mxArrayOutData = NULL;
  int32_T c9_u;
  const mxArray *c9_y = NULL;
  SFc9_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc9_Model_01InstanceStruct *)chartInstanceVoid;
  c9_mxArrayOutData = NULL;
  c9_u = *(int32_T *)c9_inData;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", &c9_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c9_mxArrayOutData, c9_y, false);
  return c9_mxArrayOutData;
}

static int32_T c9_e_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId)
{
  int32_T c9_y;
  int32_T c9_i59;
  (void)chartInstance;
  sf_mex_import(c9_parentId, sf_mex_dup(c9_u), &c9_i59, 1, 6, 0U, 0, 0U, 0);
  c9_y = c9_i59;
  sf_mex_destroy(&c9_u);
  return c9_y;
}

static void c9_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData)
{
  const mxArray *c9_b_sfEvent;
  const char_T *c9_identifier;
  emlrtMsgIdentifier c9_thisId;
  int32_T c9_y;
  SFc9_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc9_Model_01InstanceStruct *)chartInstanceVoid;
  c9_b_sfEvent = sf_mex_dup(c9_mxArrayInData);
  c9_identifier = c9_varName;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_y = c9_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_b_sfEvent),
    &c9_thisId);
  sf_mex_destroy(&c9_b_sfEvent);
  *(int32_T *)c9_outData = c9_y;
  sf_mex_destroy(&c9_mxArrayInData);
}

static uint8_T c9_f_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_b_is_active_c9_Model_01, const char_T *c9_identifier)
{
  uint8_T c9_y;
  emlrtMsgIdentifier c9_thisId;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_y = c9_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c9_b_is_active_c9_Model_01), &c9_thisId);
  sf_mex_destroy(&c9_b_is_active_c9_Model_01);
  return c9_y;
}

static uint8_T c9_g_emlrt_marshallIn(SFc9_Model_01InstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId)
{
  uint8_T c9_y;
  uint8_T c9_u0;
  (void)chartInstance;
  sf_mex_import(c9_parentId, sf_mex_dup(c9_u), &c9_u0, 1, 3, 0U, 0, 0U, 0);
  c9_y = c9_u0;
  sf_mex_destroy(&c9_u);
  return c9_y;
}

static void c9_b_eml_xgemm(SFc9_Model_01InstanceStruct *chartInstance, real_T
  c9_A[9], real_T c9_B[3], real_T c9_C[3])
{
  int32_T c9_i60;
  int32_T c9_i61;
  int32_T c9_i62;
  (void)chartInstance;
  for (c9_i60 = 0; c9_i60 < 3; c9_i60++) {
    c9_C[c9_i60] = 0.0;
    c9_i61 = 0;
    for (c9_i62 = 0; c9_i62 < 3; c9_i62++) {
      c9_C[c9_i60] += c9_A[c9_i61 + c9_i60] * c9_B[c9_i62];
      c9_i61 += 3;
    }
  }
}

static void init_dsm_address_info(SFc9_Model_01InstanceStruct *chartInstance)
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

void sf_c9_Model_01_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4284250812U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1167420104U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3797578767U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3695247764U);
}

mxArray *sf_c9_Model_01_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("Sq0CQBHL4wQCOelA2JdeiH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
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
      pr[0] = (double)(3);
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
      pr[0] = (double)(3);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c9_Model_01_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c9_Model_01_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c9_Model_01(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"r_ddot\",},{M[8],M[0],T\"is_active_c9_Model_01\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c9_Model_01_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc9_Model_01InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc9_Model_01InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_01MachineNumber_,
           9,
           1,
           1,
           0,
           5,
           0,
           0,
           0,
           0,
           1,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"n");
          _SFD_SET_DATA_PROPS(1,2,0,1,"r_ddot");
          _SFD_SET_DATA_PROPS(2,1,1,0,"r");
          _SFD_SET_DATA_PROPS(3,1,1,0,"e_force");
          _SFD_SET_DATA_PROPS(4,1,1,0,"r_dot");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,444);
        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"fn_VectorToSkewSymmetricTensor",0,-1,433);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c9_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c9_sf_marshallOut,(MexInFcnForType)
            c9_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c9_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c9_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c9_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c9_n;
          real_T (*c9_r_ddot)[3];
          real_T (*c9_r)[3];
          real_T (*c9_e_force)[3];
          real_T (*c9_r_dot)[3];
          c9_r_dot = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 3);
          c9_e_force = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 2);
          c9_r = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
          c9_r_ddot = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
          c9_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c9_n);
          _SFD_SET_DATA_VALUE_PTR(1U, *c9_r_ddot);
          _SFD_SET_DATA_VALUE_PTR(2U, *c9_r);
          _SFD_SET_DATA_VALUE_PTR(3U, *c9_e_force);
          _SFD_SET_DATA_VALUE_PTR(4U, *c9_r_dot);
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
  return "23i2VtcMhW0jI67I50N18C";
}

static void sf_opaque_initialize_c9_Model_01(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc9_Model_01InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c9_Model_01((SFc9_Model_01InstanceStruct*) chartInstanceVar);
  initialize_c9_Model_01((SFc9_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c9_Model_01(void *chartInstanceVar)
{
  enable_c9_Model_01((SFc9_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c9_Model_01(void *chartInstanceVar)
{
  disable_c9_Model_01((SFc9_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c9_Model_01(void *chartInstanceVar)
{
  sf_gateway_c9_Model_01((SFc9_Model_01InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c9_Model_01(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c9_Model_01((SFc9_Model_01InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c9_Model_01();/* state var info */
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

extern void sf_internal_set_sim_state_c9_Model_01(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c9_Model_01();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c9_Model_01((SFc9_Model_01InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c9_Model_01(SimStruct* S)
{
  return sf_internal_get_sim_state_c9_Model_01(S);
}

static void sf_opaque_set_sim_state_c9_Model_01(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c9_Model_01(S, st);
}

static void sf_opaque_terminate_c9_Model_01(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc9_Model_01InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_01_optimization_info();
    }

    finalize_c9_Model_01((SFc9_Model_01InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc9_Model_01((SFc9_Model_01InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c9_Model_01(SimStruct *S)
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
    initialize_params_c9_Model_01((SFc9_Model_01InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c9_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_01_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,9);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,9,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,9,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,9);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,9,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,9,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 4; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,9);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2147268498U));
  ssSetChecksum1(S,(4270389843U));
  ssSetChecksum2(S,(425553704U));
  ssSetChecksum3(S,(218522884U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c9_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c9_Model_01(SimStruct *S)
{
  SFc9_Model_01InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc9_Model_01InstanceStruct *)utMalloc(sizeof
    (SFc9_Model_01InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc9_Model_01InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c9_Model_01;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c9_Model_01;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c9_Model_01;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c9_Model_01;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c9_Model_01;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c9_Model_01;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c9_Model_01;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c9_Model_01;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c9_Model_01;
  chartInstance->chartInfo.mdlStart = mdlStart_c9_Model_01;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c9_Model_01;
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

void c9_Model_01_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c9_Model_01(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c9_Model_01(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c9_Model_01(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c9_Model_01_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

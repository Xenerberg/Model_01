/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_01_sfun.h"
#include "c21_Model_01.h"
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
static const char * c21_debug_family_names[16] = { "omega_norm", "q_omega", "c",
  "s", "exponential_first", "exponential_second", "n_v", "q_n", "exponential_om",
  "nargin", "nargout", "omega", "t_delta", "n", "q", "q_nominal" };

static const char * c21_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

static const char * c21_c_debug_family_names[8] = { "q_v", "q_0", "cross_q_v",
  "nargin", "nargout", "q", "flag", "CrossTensor" };

/* Function Declarations */
static void initialize_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance);
static void initialize_params_c21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance);
static void enable_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance);
static void disable_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance);
static void c21_update_debugger_state_c21_Model_01(SFc21_Model_01InstanceStruct *
  chartInstance);
static const mxArray *get_sim_state_c21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance);
static void set_sim_state_c21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance, const mxArray *c21_st);
static void finalize_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance);
static void sf_gateway_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance);
static void c21_chartstep_c21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance);
static void initSimStructsc21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance);
static void c21_fn_CrossTensor(SFc21_Model_01InstanceStruct *chartInstance,
  real_T c21_q[4], real_T c21_flag, real_T c21_CrossTensor[16]);
static void init_script_number_translation(uint32_T c21_machineNumber, uint32_T
  c21_chartNumber, uint32_T c21_instanceNumber);
static const mxArray *c21_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData);
static void c21_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_q_nominal, const char_T *c21_identifier, real_T c21_y[4]);
static void c21_b_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[4]);
static void c21_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData);
static const mxArray *c21_b_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData);
static const mxArray *c21_c_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData);
static real_T c21_c_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId);
static void c21_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData);
static const mxArray *c21_d_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData);
static void c21_d_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[16]);
static void c21_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData);
static void c21_e_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[3]);
static void c21_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData);
static const mxArray *c21_e_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData);
static void c21_f_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[9]);
static void c21_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData);
static void c21_info_helper(const mxArray **c21_info);
static const mxArray *c21_emlrt_marshallOut(const char * c21_u);
static const mxArray *c21_b_emlrt_marshallOut(const uint32_T c21_u);
static void c21_b_info_helper(const mxArray **c21_info);
static real_T c21_norm(SFc21_Model_01InstanceStruct *chartInstance, real_T
  c21_x[3]);
static void c21_threshold(SFc21_Model_01InstanceStruct *chartInstance);
static void c21_realmin(SFc21_Model_01InstanceStruct *chartInstance);
static void c21_eml_scalar_eg(SFc21_Model_01InstanceStruct *chartInstance);
static void c21_eml_xgemm(SFc21_Model_01InstanceStruct *chartInstance, real_T
  c21_A[16], real_T c21_B[4], real_T c21_C[4], real_T c21_b_C[4]);
static const mxArray *c21_f_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData);
static int32_T c21_g_emlrt_marshallIn(SFc21_Model_01InstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId);
static void c21_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData);
static uint8_T c21_h_emlrt_marshallIn(SFc21_Model_01InstanceStruct
  *chartInstance, const mxArray *c21_b_is_active_c21_Model_01, const char_T
  *c21_identifier);
static uint8_T c21_i_emlrt_marshallIn(SFc21_Model_01InstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId);
static void c21_b_eml_xgemm(SFc21_Model_01InstanceStruct *chartInstance, real_T
  c21_A[16], real_T c21_B[4], real_T c21_C[4]);
static void init_dsm_address_info(SFc21_Model_01InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance)
{
  chartInstance->c21_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c21_is_active_c21_Model_01 = 0U;
}

static void initialize_params_c21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c21_update_debugger_state_c21_Model_01(SFc21_Model_01InstanceStruct *
  chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance)
{
  const mxArray *c21_st;
  const mxArray *c21_y = NULL;
  int32_T c21_i0;
  real_T c21_u[4];
  const mxArray *c21_b_y = NULL;
  uint8_T c21_hoistedGlobal;
  uint8_T c21_b_u;
  const mxArray *c21_c_y = NULL;
  real_T (*c21_q_nominal)[4];
  c21_q_nominal = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c21_st = NULL;
  c21_st = NULL;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_createcellmatrix(2, 1), false);
  for (c21_i0 = 0; c21_i0 < 4; c21_i0++) {
    c21_u[c21_i0] = (*c21_q_nominal)[c21_i0];
  }

  c21_b_y = NULL;
  sf_mex_assign(&c21_b_y, sf_mex_create("y", c21_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_setcell(c21_y, 0, c21_b_y);
  c21_hoistedGlobal = chartInstance->c21_is_active_c21_Model_01;
  c21_b_u = c21_hoistedGlobal;
  c21_c_y = NULL;
  sf_mex_assign(&c21_c_y, sf_mex_create("y", &c21_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c21_y, 1, c21_c_y);
  sf_mex_assign(&c21_st, c21_y, false);
  return c21_st;
}

static void set_sim_state_c21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance, const mxArray *c21_st)
{
  const mxArray *c21_u;
  real_T c21_dv0[4];
  int32_T c21_i1;
  real_T (*c21_q_nominal)[4];
  c21_q_nominal = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c21_doneDoubleBufferReInit = true;
  c21_u = sf_mex_dup(c21_st);
  c21_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c21_u, 0)),
                       "q_nominal", c21_dv0);
  for (c21_i1 = 0; c21_i1 < 4; c21_i1++) {
    (*c21_q_nominal)[c21_i1] = c21_dv0[c21_i1];
  }

  chartInstance->c21_is_active_c21_Model_01 = c21_h_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c21_u, 1)),
     "is_active_c21_Model_01");
  sf_mex_destroy(&c21_u);
  c21_update_debugger_state_c21_Model_01(chartInstance);
  sf_mex_destroy(&c21_st);
}

static void finalize_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c21_Model_01(SFc21_Model_01InstanceStruct *chartInstance)
{
  int32_T c21_i2;
  int32_T c21_i3;
  int32_T c21_i4;
  real_T *c21_t_delta;
  real_T *c21_n;
  real_T (*c21_q)[4];
  real_T (*c21_q_nominal)[4];
  real_T (*c21_omega)[3];
  c21_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 3);
  c21_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c21_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c21_q_nominal = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c21_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 20U, chartInstance->c21_sfEvent);
  for (c21_i2 = 0; c21_i2 < 3; c21_i2++) {
    _SFD_DATA_RANGE_CHECK((*c21_omega)[c21_i2], 0U);
  }

  chartInstance->c21_sfEvent = CALL_EVENT;
  c21_chartstep_c21_Model_01(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_01MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c21_i3 = 0; c21_i3 < 4; c21_i3++) {
    _SFD_DATA_RANGE_CHECK((*c21_q_nominal)[c21_i3], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c21_t_delta, 2U);
  _SFD_DATA_RANGE_CHECK(*c21_n, 3U);
  for (c21_i4 = 0; c21_i4 < 4; c21_i4++) {
    _SFD_DATA_RANGE_CHECK((*c21_q)[c21_i4], 4U);
  }
}

static void c21_chartstep_c21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance)
{
  real_T c21_hoistedGlobal;
  real_T c21_b_hoistedGlobal;
  int32_T c21_i5;
  real_T c21_omega[3];
  real_T c21_t_delta;
  real_T c21_n;
  int32_T c21_i6;
  real_T c21_q[4];
  uint32_T c21_debug_family_var_map[16];
  real_T c21_omega_norm;
  real_T c21_q_omega[4];
  real_T c21_c;
  real_T c21_s;
  real_T c21_exponential_first[16];
  real_T c21_exponential_second[16];
  real_T c21_n_v[3];
  real_T c21_q_n[4];
  real_T c21_exponential_om[16];
  real_T c21_nargin = 4.0;
  real_T c21_nargout = 1.0;
  real_T c21_q_nominal[4];
  int32_T c21_i7;
  int32_T c21_i8;
  real_T c21_b_omega[3];
  int32_T c21_i9;
  real_T c21_A;
  real_T c21_x;
  real_T c21_b_x;
  real_T c21_c_x;
  real_T c21_y;
  real_T c21_d_x;
  real_T c21_e_x;
  real_T c21_b_A;
  real_T c21_f_x;
  real_T c21_g_x;
  real_T c21_h_x;
  real_T c21_b_y;
  real_T c21_i_x;
  real_T c21_j_x;
  real_T c21_a;
  int32_T c21_i10;
  static real_T c21_b[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c21_c_A;
  real_T c21_B;
  real_T c21_k_x;
  real_T c21_c_y;
  real_T c21_l_x;
  real_T c21_d_y;
  real_T c21_m_x;
  real_T c21_e_y;
  real_T c21_f_y;
  real_T c21_b_a;
  int32_T c21_i11;
  real_T c21_b_q_omega[4];
  real_T c21_b_b[16];
  int32_T c21_i12;
  int32_T c21_i13;
  int32_T c21_i14;
  boolean_T c21_c_b[16];
  int32_T c21_k;
  int32_T c21_i;
  int32_T c21_b_i;
  int32_T c21_c_a;
  int32_T c21_d_a;
  const mxArray *c21_g_y = NULL;
  int32_T c21_tmp_sizes;
  int32_T c21_loop_ub;
  int32_T c21_i15;
  int32_T c21_tmp_data[16];
  int32_T c21_b_tmp_sizes;
  int32_T c21_j;
  int32_T c21_c_i;
  int32_T c21_d_i;
  int32_T c21_b_tmp_data[16];
  int32_T c21_e_a;
  int32_T c21_f_a;
  int32_T c21_b_loop_ub;
  int32_T c21_i16;
  int32_T c21_i17;
  int32_T c21_i18;
  int32_T c21_i19;
  int32_T c21_i20;
  real_T c21_d_b[4];
  int32_T c21_i21;
  int32_T c21_i22;
  int32_T c21_i23;
  real_T c21_dv1[16];
  int32_T c21_i24;
  real_T c21_dv2[4];
  int32_T c21_i25;
  real_T c21_dv3[16];
  int32_T c21_i26;
  real_T c21_dv4[4];
  int32_T c21_i27;
  real_T c21_h_y;
  real_T c21_scale;
  int32_T c21_b_k;
  int32_T c21_c_k;
  real_T c21_n_x;
  real_T c21_o_x;
  real_T c21_absxk;
  real_T c21_t;
  int32_T c21_i28;
  real_T c21_b_B;
  real_T c21_i_y;
  real_T c21_j_y;
  real_T c21_k_y;
  int32_T c21_i29;
  int32_T c21_i30;
  real_T *c21_b_t_delta;
  real_T *c21_b_n;
  real_T (*c21_b_q_nominal)[4];
  real_T (*c21_b_q)[4];
  real_T (*c21_c_omega)[3];
  c21_b_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 3);
  c21_b_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c21_b_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c21_b_q_nominal = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c21_c_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 20U, chartInstance->c21_sfEvent);
  c21_hoistedGlobal = *c21_b_t_delta;
  c21_b_hoistedGlobal = *c21_b_n;
  for (c21_i5 = 0; c21_i5 < 3; c21_i5++) {
    c21_omega[c21_i5] = (*c21_c_omega)[c21_i5];
  }

  c21_t_delta = c21_hoistedGlobal;
  c21_n = c21_b_hoistedGlobal;
  for (c21_i6 = 0; c21_i6 < 4; c21_i6++) {
    c21_q[c21_i6] = (*c21_b_q)[c21_i6];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 16U, 16U, c21_debug_family_names,
    c21_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_omega_norm, 0U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_q_omega, 1U, c21_sf_marshallOut,
    c21_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_c, 2U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_s, 3U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_exponential_first, 4U,
    c21_d_sf_marshallOut, c21_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_exponential_second, 5U,
    c21_d_sf_marshallOut, c21_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_n_v, 6U, c21_c_sf_marshallOut,
    c21_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_q_n, 7U, c21_sf_marshallOut,
    c21_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_exponential_om, 8U,
    c21_d_sf_marshallOut, c21_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_nargin, 9U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_nargout, 10U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c21_omega, 11U, c21_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c21_t_delta, 12U, c21_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c21_n, 13U, c21_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c21_q, 14U, c21_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_q_nominal, 15U, c21_sf_marshallOut,
    c21_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 3);
  for (c21_i7 = 0; c21_i7 < 4; c21_i7++) {
    c21_q_nominal[c21_i7] = c21_q[c21_i7];
  }

  _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 4);
  for (c21_i8 = 0; c21_i8 < 3; c21_i8++) {
    c21_b_omega[c21_i8] = c21_omega[c21_i8];
  }

  c21_omega_norm = c21_norm(chartInstance, c21_b_omega);
  _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 5);
  if (CV_EML_IF(0, 1, 0, c21_omega_norm != 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 6);
    for (c21_i9 = 0; c21_i9 < 3; c21_i9++) {
      c21_q_omega[c21_i9] = c21_omega[c21_i9];
    }

    c21_q_omega[3] = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 7);
    c21_A = c21_omega_norm * c21_t_delta;
    c21_x = c21_A;
    c21_b_x = c21_x;
    c21_c_x = c21_b_x;
    c21_y = c21_c_x / 2.0;
    c21_d_x = c21_y;
    c21_c = c21_d_x;
    c21_e_x = c21_c;
    c21_c = c21_e_x;
    c21_c = muDoubleScalarCos(c21_c);
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 8);
    c21_b_A = c21_omega_norm * c21_t_delta;
    c21_f_x = c21_b_A;
    c21_g_x = c21_f_x;
    c21_h_x = c21_g_x;
    c21_b_y = c21_h_x / 2.0;
    c21_i_x = c21_b_y;
    c21_s = c21_i_x;
    c21_j_x = c21_s;
    c21_s = c21_j_x;
    c21_s = muDoubleScalarSin(c21_s);
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 9);
    c21_a = c21_c + c21_s;
    for (c21_i10 = 0; c21_i10 < 16; c21_i10++) {
      c21_exponential_first[c21_i10] = c21_a * c21_b[c21_i10];
    }

    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 10);
    c21_c_A = c21_omega_norm * c21_t_delta * c21_c - 4.0 * c21_s;
    c21_B = 2.0 * c21_omega_norm;
    c21_k_x = c21_c_A;
    c21_c_y = c21_B;
    c21_l_x = c21_k_x;
    c21_d_y = c21_c_y;
    c21_m_x = c21_l_x;
    c21_e_y = c21_d_y;
    c21_f_y = c21_m_x / c21_e_y;
    c21_b_a = c21_f_y;
    for (c21_i11 = 0; c21_i11 < 4; c21_i11++) {
      c21_b_q_omega[c21_i11] = c21_q_omega[c21_i11];
    }

    c21_fn_CrossTensor(chartInstance, c21_b_q_omega, 0.0, c21_b_b);
    for (c21_i12 = 0; c21_i12 < 16; c21_i12++) {
      c21_exponential_second[c21_i12] = c21_b_a * c21_b_b[c21_i12];
    }

    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 11);
    for (c21_i13 = 0; c21_i13 < 16; c21_i13++) {
      c21_b_b[c21_i13] = c21_exponential_second[c21_i13];
    }

    for (c21_i14 = 0; c21_i14 < 16; c21_i14++) {
      c21_c_b[c21_i14] = muDoubleScalarIsNaN(c21_b_b[c21_i14]);
    }

    c21_k = 0;
    for (c21_i = 1; c21_i < 17; c21_i++) {
      c21_b_i = c21_i - 1;
      if (c21_c_b[c21_b_i]) {
        c21_c_a = c21_k;
        c21_d_a = c21_c_a + 1;
        c21_k = c21_d_a;
      }
    }

    if (c21_k <= 16) {
    } else {
      c21_g_y = NULL;
      sf_mex_assign(&c21_g_y, sf_mex_create("y", "Assertion failed.", 15, 0U, 0U,
        0U, 2, 1, strlen("Assertion failed.")), false);
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                        c21_g_y);
    }

    c21_tmp_sizes = (int32_T)_SFD_NON_NEGATIVE_CHECK("", (real_T)c21_k);
    c21_loop_ub = (int32_T)_SFD_NON_NEGATIVE_CHECK("", (real_T)c21_k) - 1;
    for (c21_i15 = 0; c21_i15 <= c21_loop_ub; c21_i15++) {
      c21_tmp_data[c21_i15] = 0;
    }

    c21_b_tmp_sizes = c21_tmp_sizes;
    c21_j = 1;
    for (c21_c_i = 1; c21_c_i < 17; c21_c_i++) {
      c21_d_i = c21_c_i;
      if (c21_c_b[c21_d_i - 1]) {
        c21_b_tmp_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c21_j, 1, c21_b_tmp_sizes,
          1, 0) - 1] = c21_d_i;
        c21_e_a = c21_j;
        c21_f_a = c21_e_a + 1;
        c21_j = c21_f_a;
      }
    }

    c21_b_loop_ub = c21_b_tmp_sizes - 1;
    for (c21_i16 = 0; c21_i16 <= c21_b_loop_ub; c21_i16++) {
      c21_exponential_second[c21_b_tmp_data[c21_i16] - 1] = 0.0;
    }

    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 14);
    c21_n_v[0] = 0.0;
    c21_n_v[1] = 0.0;
    c21_n_v[2] = c21_n;
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 15);
    for (c21_i17 = 0; c21_i17 < 3; c21_i17++) {
      c21_q_n[c21_i17] = c21_n_v[c21_i17];
    }

    c21_q_n[3] = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 17);
    for (c21_i18 = 0; c21_i18 < 16; c21_i18++) {
      c21_exponential_om[c21_i18] = c21_exponential_first[c21_i18] -
        c21_exponential_second[c21_i18];
    }

    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 21);
    for (c21_i19 = 0; c21_i19 < 16; c21_i19++) {
      c21_b_b[c21_i19] = c21_exponential_om[c21_i19];
    }

    for (c21_i20 = 0; c21_i20 < 4; c21_i20++) {
      c21_d_b[c21_i20] = c21_q[c21_i20];
    }

    c21_eml_scalar_eg(chartInstance);
    c21_eml_scalar_eg(chartInstance);
    for (c21_i21 = 0; c21_i21 < 4; c21_i21++) {
      c21_q_nominal[c21_i21] = 0.0;
    }

    for (c21_i22 = 0; c21_i22 < 4; c21_i22++) {
      c21_q_nominal[c21_i22] = 0.0;
    }

    for (c21_i23 = 0; c21_i23 < 16; c21_i23++) {
      c21_dv1[c21_i23] = c21_b_b[c21_i23];
    }

    for (c21_i24 = 0; c21_i24 < 4; c21_i24++) {
      c21_dv2[c21_i24] = c21_d_b[c21_i24];
    }

    for (c21_i25 = 0; c21_i25 < 16; c21_i25++) {
      c21_dv3[c21_i25] = c21_dv1[c21_i25];
    }

    for (c21_i26 = 0; c21_i26 < 4; c21_i26++) {
      c21_dv4[c21_i26] = c21_dv2[c21_i26];
    }

    c21_b_eml_xgemm(chartInstance, c21_dv3, c21_dv4, c21_q_nominal);
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 22);
    for (c21_i27 = 0; c21_i27 < 4; c21_i27++) {
      c21_d_b[c21_i27] = c21_q_nominal[c21_i27];
    }

    c21_threshold(chartInstance);
    c21_h_y = 0.0;
    c21_realmin(chartInstance);
    c21_scale = 2.2250738585072014E-308;
    for (c21_b_k = 1; c21_b_k < 5; c21_b_k++) {
      c21_c_k = c21_b_k - 1;
      c21_n_x = c21_d_b[c21_c_k];
      c21_o_x = c21_n_x;
      c21_absxk = muDoubleScalarAbs(c21_o_x);
      if (c21_absxk > c21_scale) {
        c21_t = c21_scale / c21_absxk;
        c21_h_y = 1.0 + c21_h_y * c21_t * c21_t;
        c21_scale = c21_absxk;
      } else {
        c21_t = c21_absxk / c21_scale;
        c21_h_y += c21_t * c21_t;
      }
    }

    c21_h_y = c21_scale * muDoubleScalarSqrt(c21_h_y);
    for (c21_i28 = 0; c21_i28 < 4; c21_i28++) {
      c21_d_b[c21_i28] = c21_q_nominal[c21_i28];
    }

    c21_b_B = c21_h_y;
    c21_i_y = c21_b_B;
    c21_j_y = c21_i_y;
    c21_k_y = c21_j_y;
    for (c21_i29 = 0; c21_i29 < 4; c21_i29++) {
      c21_q_nominal[c21_i29] = c21_d_b[c21_i29] / c21_k_y;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, -22);
  _SFD_SYMBOL_SCOPE_POP();
  for (c21_i30 = 0; c21_i30 < 4; c21_i30++) {
    (*c21_b_q_nominal)[c21_i30] = c21_q_nominal[c21_i30];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 20U, chartInstance->c21_sfEvent);
}

static void initSimStructsc21_Model_01(SFc21_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c21_fn_CrossTensor(SFc21_Model_01InstanceStruct *chartInstance,
  real_T c21_q[4], real_T c21_flag, real_T c21_CrossTensor[16])
{
  uint32_T c21_debug_family_var_map[8];
  real_T c21_q_v[3];
  real_T c21_q_0;
  real_T c21_cross_q_v[9];
  real_T c21_nargin = 2.0;
  real_T c21_nargout = 1.0;
  int32_T c21_i31;
  int32_T c21_i32;
  int32_T c21_i33;
  real_T c21_v[3];
  uint32_T c21_b_debug_family_var_map[4];
  real_T c21_b_nargin = 1.0;
  real_T c21_b_nargout = 1.0;
  int32_T c21_i34;
  int32_T c21_i35;
  int32_T c21_i36;
  int32_T c21_i37;
  int32_T c21_i38;
  int32_T c21_i39;
  int32_T c21_i40;
  int32_T c21_i41;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c21_c_debug_family_names,
    c21_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_q_v, 0U, c21_c_sf_marshallOut,
    c21_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_q_0, 1U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_cross_q_v, 2U, c21_e_sf_marshallOut,
    c21_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_nargin, 3U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_nargout, 4U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_q, 5U, c21_sf_marshallOut,
    c21_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_flag, 6U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_CrossTensor, 7U, c21_d_sf_marshallOut,
    c21_c_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c21_sfEvent, 4);
  for (c21_i31 = 0; c21_i31 < 16; c21_i31++) {
    c21_CrossTensor[c21_i31] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c21_sfEvent, 5);
  for (c21_i32 = 0; c21_i32 < 3; c21_i32++) {
    c21_q_v[c21_i32] = c21_q[c21_i32];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c21_sfEvent, 6);
  c21_q_0 = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c21_sfEvent, 7);
  for (c21_i33 = 0; c21_i33 < 3; c21_i33++) {
    c21_v[c21_i33] = c21_q_v[c21_i33];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c21_b_debug_family_names,
    c21_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_b_nargin, 0U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_b_nargout, 1U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_v, 2U, c21_c_sf_marshallOut,
    c21_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_cross_q_v, 3U, c21_e_sf_marshallOut,
    c21_e_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 2);
  for (c21_i34 = 0; c21_i34 < 9; c21_i34++) {
    c21_cross_q_v[c21_i34] = 0.0;
  }

  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 3);
  c21_cross_q_v[0] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 4);
  c21_cross_q_v[4] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 5);
  c21_cross_q_v[8] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 6);
  c21_cross_q_v[3] = -c21_v[2];
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 7);
  c21_cross_q_v[6] = c21_v[1];
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 8);
  c21_cross_q_v[7] = -c21_v[0];
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 9);
  c21_cross_q_v[1] = c21_v[2];
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 10);
  c21_cross_q_v[2] = -c21_v[1];
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, 11);
  c21_cross_q_v[5] = c21_v[0];
  _SFD_SCRIPT_CALL(1U, chartInstance->c21_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_SCRIPT_CALL(0U, chartInstance->c21_sfEvent, 8);
  CV_SCRIPT_SWITCH(0, 0, 1);
  _SFD_SCRIPT_CALL(0U, chartInstance->c21_sfEvent, 10);
  c21_i35 = 0;
  c21_i36 = 0;
  for (c21_i37 = 0; c21_i37 < 3; c21_i37++) {
    for (c21_i38 = 0; c21_i38 < 3; c21_i38++) {
      c21_CrossTensor[c21_i38 + c21_i35] = -c21_cross_q_v[c21_i38 + c21_i36];
    }

    c21_i35 += 4;
    c21_i36 += 3;
  }

  for (c21_i39 = 0; c21_i39 < 3; c21_i39++) {
    c21_CrossTensor[c21_i39 + 12] = c21_q_v[c21_i39];
  }

  c21_i40 = 0;
  for (c21_i41 = 0; c21_i41 < 3; c21_i41++) {
    c21_CrossTensor[c21_i40 + 3] = -c21_q_v[c21_i41];
    c21_i40 += 4;
  }

  c21_CrossTensor[15] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c21_sfEvent, -12);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c21_machineNumber, uint32_T
  c21_chartNumber, uint32_T c21_instanceNumber)
{
  (void)c21_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c21_chartNumber, c21_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg\\Documents\\MATLAB\\Model_01\\fn_CrossTensor.m"));
  _SFD_SCRIPT_TRANSLATION(c21_chartNumber, c21_instanceNumber, 1U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c21_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData)
{
  const mxArray *c21_mxArrayOutData = NULL;
  int32_T c21_i42;
  real_T c21_b_inData[4];
  int32_T c21_i43;
  real_T c21_u[4];
  const mxArray *c21_y = NULL;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_mxArrayOutData = NULL;
  for (c21_i42 = 0; c21_i42 < 4; c21_i42++) {
    c21_b_inData[c21_i42] = (*(real_T (*)[4])c21_inData)[c21_i42];
  }

  for (c21_i43 = 0; c21_i43 < 4; c21_i43++) {
    c21_u[c21_i43] = c21_b_inData[c21_i43];
  }

  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", c21_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c21_mxArrayOutData, c21_y, false);
  return c21_mxArrayOutData;
}

static void c21_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_q_nominal, const char_T *c21_identifier, real_T c21_y[4])
{
  emlrtMsgIdentifier c21_thisId;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_q_nominal), &c21_thisId,
    c21_y);
  sf_mex_destroy(&c21_q_nominal);
}

static void c21_b_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[4])
{
  real_T c21_dv5[4];
  int32_T c21_i44;
  (void)chartInstance;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), c21_dv5, 1, 0, 0U, 1, 0U, 1, 4);
  for (c21_i44 = 0; c21_i44 < 4; c21_i44++) {
    c21_y[c21_i44] = c21_dv5[c21_i44];
  }

  sf_mex_destroy(&c21_u);
}

static void c21_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData)
{
  const mxArray *c21_q_nominal;
  const char_T *c21_identifier;
  emlrtMsgIdentifier c21_thisId;
  real_T c21_y[4];
  int32_T c21_i45;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_q_nominal = sf_mex_dup(c21_mxArrayInData);
  c21_identifier = c21_varName;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_q_nominal), &c21_thisId,
    c21_y);
  sf_mex_destroy(&c21_q_nominal);
  for (c21_i45 = 0; c21_i45 < 4; c21_i45++) {
    (*(real_T (*)[4])c21_outData)[c21_i45] = c21_y[c21_i45];
  }

  sf_mex_destroy(&c21_mxArrayInData);
}

static const mxArray *c21_b_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData)
{
  const mxArray *c21_mxArrayOutData = NULL;
  real_T c21_u;
  const mxArray *c21_y = NULL;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_mxArrayOutData = NULL;
  c21_u = *(real_T *)c21_inData;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", &c21_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c21_mxArrayOutData, c21_y, false);
  return c21_mxArrayOutData;
}

static const mxArray *c21_c_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData)
{
  const mxArray *c21_mxArrayOutData = NULL;
  int32_T c21_i46;
  real_T c21_b_inData[3];
  int32_T c21_i47;
  real_T c21_u[3];
  const mxArray *c21_y = NULL;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_mxArrayOutData = NULL;
  for (c21_i46 = 0; c21_i46 < 3; c21_i46++) {
    c21_b_inData[c21_i46] = (*(real_T (*)[3])c21_inData)[c21_i46];
  }

  for (c21_i47 = 0; c21_i47 < 3; c21_i47++) {
    c21_u[c21_i47] = c21_b_inData[c21_i47];
  }

  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", c21_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c21_mxArrayOutData, c21_y, false);
  return c21_mxArrayOutData;
}

static real_T c21_c_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId)
{
  real_T c21_y;
  real_T c21_d0;
  (void)chartInstance;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), &c21_d0, 1, 0, 0U, 0, 0U, 0);
  c21_y = c21_d0;
  sf_mex_destroy(&c21_u);
  return c21_y;
}

static void c21_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData)
{
  const mxArray *c21_nargout;
  const char_T *c21_identifier;
  emlrtMsgIdentifier c21_thisId;
  real_T c21_y;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_nargout = sf_mex_dup(c21_mxArrayInData);
  c21_identifier = c21_varName;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_y = c21_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_nargout),
    &c21_thisId);
  sf_mex_destroy(&c21_nargout);
  *(real_T *)c21_outData = c21_y;
  sf_mex_destroy(&c21_mxArrayInData);
}

static const mxArray *c21_d_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData)
{
  const mxArray *c21_mxArrayOutData = NULL;
  int32_T c21_i48;
  int32_T c21_i49;
  int32_T c21_i50;
  real_T c21_b_inData[16];
  int32_T c21_i51;
  int32_T c21_i52;
  int32_T c21_i53;
  real_T c21_u[16];
  const mxArray *c21_y = NULL;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_mxArrayOutData = NULL;
  c21_i48 = 0;
  for (c21_i49 = 0; c21_i49 < 4; c21_i49++) {
    for (c21_i50 = 0; c21_i50 < 4; c21_i50++) {
      c21_b_inData[c21_i50 + c21_i48] = (*(real_T (*)[16])c21_inData)[c21_i50 +
        c21_i48];
    }

    c21_i48 += 4;
  }

  c21_i51 = 0;
  for (c21_i52 = 0; c21_i52 < 4; c21_i52++) {
    for (c21_i53 = 0; c21_i53 < 4; c21_i53++) {
      c21_u[c21_i53 + c21_i51] = c21_b_inData[c21_i53 + c21_i51];
    }

    c21_i51 += 4;
  }

  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", c21_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c21_mxArrayOutData, c21_y, false);
  return c21_mxArrayOutData;
}

static void c21_d_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[16])
{
  real_T c21_dv6[16];
  int32_T c21_i54;
  (void)chartInstance;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), c21_dv6, 1, 0, 0U, 1, 0U, 2, 4,
                4);
  for (c21_i54 = 0; c21_i54 < 16; c21_i54++) {
    c21_y[c21_i54] = c21_dv6[c21_i54];
  }

  sf_mex_destroy(&c21_u);
}

static void c21_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData)
{
  const mxArray *c21_exponential_om;
  const char_T *c21_identifier;
  emlrtMsgIdentifier c21_thisId;
  real_T c21_y[16];
  int32_T c21_i55;
  int32_T c21_i56;
  int32_T c21_i57;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_exponential_om = sf_mex_dup(c21_mxArrayInData);
  c21_identifier = c21_varName;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_exponential_om),
    &c21_thisId, c21_y);
  sf_mex_destroy(&c21_exponential_om);
  c21_i55 = 0;
  for (c21_i56 = 0; c21_i56 < 4; c21_i56++) {
    for (c21_i57 = 0; c21_i57 < 4; c21_i57++) {
      (*(real_T (*)[16])c21_outData)[c21_i57 + c21_i55] = c21_y[c21_i57 +
        c21_i55];
    }

    c21_i55 += 4;
  }

  sf_mex_destroy(&c21_mxArrayInData);
}

static void c21_e_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[3])
{
  real_T c21_dv7[3];
  int32_T c21_i58;
  (void)chartInstance;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), c21_dv7, 1, 0, 0U, 1, 0U, 1, 3);
  for (c21_i58 = 0; c21_i58 < 3; c21_i58++) {
    c21_y[c21_i58] = c21_dv7[c21_i58];
  }

  sf_mex_destroy(&c21_u);
}

static void c21_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData)
{
  const mxArray *c21_n_v;
  const char_T *c21_identifier;
  emlrtMsgIdentifier c21_thisId;
  real_T c21_y[3];
  int32_T c21_i59;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_n_v = sf_mex_dup(c21_mxArrayInData);
  c21_identifier = c21_varName;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_n_v), &c21_thisId, c21_y);
  sf_mex_destroy(&c21_n_v);
  for (c21_i59 = 0; c21_i59 < 3; c21_i59++) {
    (*(real_T (*)[3])c21_outData)[c21_i59] = c21_y[c21_i59];
  }

  sf_mex_destroy(&c21_mxArrayInData);
}

static const mxArray *c21_e_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData)
{
  const mxArray *c21_mxArrayOutData = NULL;
  int32_T c21_i60;
  int32_T c21_i61;
  int32_T c21_i62;
  real_T c21_b_inData[9];
  int32_T c21_i63;
  int32_T c21_i64;
  int32_T c21_i65;
  real_T c21_u[9];
  const mxArray *c21_y = NULL;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_mxArrayOutData = NULL;
  c21_i60 = 0;
  for (c21_i61 = 0; c21_i61 < 3; c21_i61++) {
    for (c21_i62 = 0; c21_i62 < 3; c21_i62++) {
      c21_b_inData[c21_i62 + c21_i60] = (*(real_T (*)[9])c21_inData)[c21_i62 +
        c21_i60];
    }

    c21_i60 += 3;
  }

  c21_i63 = 0;
  for (c21_i64 = 0; c21_i64 < 3; c21_i64++) {
    for (c21_i65 = 0; c21_i65 < 3; c21_i65++) {
      c21_u[c21_i65 + c21_i63] = c21_b_inData[c21_i65 + c21_i63];
    }

    c21_i63 += 3;
  }

  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", c21_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c21_mxArrayOutData, c21_y, false);
  return c21_mxArrayOutData;
}

static void c21_f_emlrt_marshallIn(SFc21_Model_01InstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[9])
{
  real_T c21_dv8[9];
  int32_T c21_i66;
  (void)chartInstance;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), c21_dv8, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c21_i66 = 0; c21_i66 < 9; c21_i66++) {
    c21_y[c21_i66] = c21_dv8[c21_i66];
  }

  sf_mex_destroy(&c21_u);
}

static void c21_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData)
{
  const mxArray *c21_SkewSymmetricTensor;
  const char_T *c21_identifier;
  emlrtMsgIdentifier c21_thisId;
  real_T c21_y[9];
  int32_T c21_i67;
  int32_T c21_i68;
  int32_T c21_i69;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_SkewSymmetricTensor = sf_mex_dup(c21_mxArrayInData);
  c21_identifier = c21_varName;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_SkewSymmetricTensor),
    &c21_thisId, c21_y);
  sf_mex_destroy(&c21_SkewSymmetricTensor);
  c21_i67 = 0;
  for (c21_i68 = 0; c21_i68 < 3; c21_i68++) {
    for (c21_i69 = 0; c21_i69 < 3; c21_i69++) {
      (*(real_T (*)[9])c21_outData)[c21_i69 + c21_i67] = c21_y[c21_i69 + c21_i67];
    }

    c21_i67 += 3;
  }

  sf_mex_destroy(&c21_mxArrayInData);
}

const mxArray *sf_c21_Model_01_get_eml_resolved_functions_info(void)
{
  const mxArray *c21_nameCaptureInfo = NULL;
  c21_nameCaptureInfo = NULL;
  sf_mex_assign(&c21_nameCaptureInfo, sf_mex_createstruct("structure", 2, 75, 1),
                false);
  c21_info_helper(&c21_nameCaptureInfo);
  c21_b_info_helper(&c21_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c21_nameCaptureInfo);
  return c21_nameCaptureInfo;
}

static void c21_info_helper(const mxArray **c21_info)
{
  const mxArray *c21_rhs0 = NULL;
  const mxArray *c21_lhs0 = NULL;
  const mxArray *c21_rhs1 = NULL;
  const mxArray *c21_lhs1 = NULL;
  const mxArray *c21_rhs2 = NULL;
  const mxArray *c21_lhs2 = NULL;
  const mxArray *c21_rhs3 = NULL;
  const mxArray *c21_lhs3 = NULL;
  const mxArray *c21_rhs4 = NULL;
  const mxArray *c21_lhs4 = NULL;
  const mxArray *c21_rhs5 = NULL;
  const mxArray *c21_lhs5 = NULL;
  const mxArray *c21_rhs6 = NULL;
  const mxArray *c21_lhs6 = NULL;
  const mxArray *c21_rhs7 = NULL;
  const mxArray *c21_lhs7 = NULL;
  const mxArray *c21_rhs8 = NULL;
  const mxArray *c21_lhs8 = NULL;
  const mxArray *c21_rhs9 = NULL;
  const mxArray *c21_lhs9 = NULL;
  const mxArray *c21_rhs10 = NULL;
  const mxArray *c21_lhs10 = NULL;
  const mxArray *c21_rhs11 = NULL;
  const mxArray *c21_lhs11 = NULL;
  const mxArray *c21_rhs12 = NULL;
  const mxArray *c21_lhs12 = NULL;
  const mxArray *c21_rhs13 = NULL;
  const mxArray *c21_lhs13 = NULL;
  const mxArray *c21_rhs14 = NULL;
  const mxArray *c21_lhs14 = NULL;
  const mxArray *c21_rhs15 = NULL;
  const mxArray *c21_lhs15 = NULL;
  const mxArray *c21_rhs16 = NULL;
  const mxArray *c21_lhs16 = NULL;
  const mxArray *c21_rhs17 = NULL;
  const mxArray *c21_lhs17 = NULL;
  const mxArray *c21_rhs18 = NULL;
  const mxArray *c21_lhs18 = NULL;
  const mxArray *c21_rhs19 = NULL;
  const mxArray *c21_lhs19 = NULL;
  const mxArray *c21_rhs20 = NULL;
  const mxArray *c21_lhs20 = NULL;
  const mxArray *c21_rhs21 = NULL;
  const mxArray *c21_lhs21 = NULL;
  const mxArray *c21_rhs22 = NULL;
  const mxArray *c21_lhs22 = NULL;
  const mxArray *c21_rhs23 = NULL;
  const mxArray *c21_lhs23 = NULL;
  const mxArray *c21_rhs24 = NULL;
  const mxArray *c21_lhs24 = NULL;
  const mxArray *c21_rhs25 = NULL;
  const mxArray *c21_lhs25 = NULL;
  const mxArray *c21_rhs26 = NULL;
  const mxArray *c21_lhs26 = NULL;
  const mxArray *c21_rhs27 = NULL;
  const mxArray *c21_lhs27 = NULL;
  const mxArray *c21_rhs28 = NULL;
  const mxArray *c21_lhs28 = NULL;
  const mxArray *c21_rhs29 = NULL;
  const mxArray *c21_lhs29 = NULL;
  const mxArray *c21_rhs30 = NULL;
  const mxArray *c21_lhs30 = NULL;
  const mxArray *c21_rhs31 = NULL;
  const mxArray *c21_lhs31 = NULL;
  const mxArray *c21_rhs32 = NULL;
  const mxArray *c21_lhs32 = NULL;
  const mxArray *c21_rhs33 = NULL;
  const mxArray *c21_lhs33 = NULL;
  const mxArray *c21_rhs34 = NULL;
  const mxArray *c21_lhs34 = NULL;
  const mxArray *c21_rhs35 = NULL;
  const mxArray *c21_lhs35 = NULL;
  const mxArray *c21_rhs36 = NULL;
  const mxArray *c21_lhs36 = NULL;
  const mxArray *c21_rhs37 = NULL;
  const mxArray *c21_lhs37 = NULL;
  const mxArray *c21_rhs38 = NULL;
  const mxArray *c21_lhs38 = NULL;
  const mxArray *c21_rhs39 = NULL;
  const mxArray *c21_lhs39 = NULL;
  const mxArray *c21_rhs40 = NULL;
  const mxArray *c21_lhs40 = NULL;
  const mxArray *c21_rhs41 = NULL;
  const mxArray *c21_lhs41 = NULL;
  const mxArray *c21_rhs42 = NULL;
  const mxArray *c21_lhs42 = NULL;
  const mxArray *c21_rhs43 = NULL;
  const mxArray *c21_lhs43 = NULL;
  const mxArray *c21_rhs44 = NULL;
  const mxArray *c21_lhs44 = NULL;
  const mxArray *c21_rhs45 = NULL;
  const mxArray *c21_lhs45 = NULL;
  const mxArray *c21_rhs46 = NULL;
  const mxArray *c21_lhs46 = NULL;
  const mxArray *c21_rhs47 = NULL;
  const mxArray *c21_lhs47 = NULL;
  const mxArray *c21_rhs48 = NULL;
  const mxArray *c21_lhs48 = NULL;
  const mxArray *c21_rhs49 = NULL;
  const mxArray *c21_lhs49 = NULL;
  const mxArray *c21_rhs50 = NULL;
  const mxArray *c21_lhs50 = NULL;
  const mxArray *c21_rhs51 = NULL;
  const mxArray *c21_lhs51 = NULL;
  const mxArray *c21_rhs52 = NULL;
  const mxArray *c21_lhs52 = NULL;
  const mxArray *c21_rhs53 = NULL;
  const mxArray *c21_lhs53 = NULL;
  const mxArray *c21_rhs54 = NULL;
  const mxArray *c21_lhs54 = NULL;
  const mxArray *c21_rhs55 = NULL;
  const mxArray *c21_lhs55 = NULL;
  const mxArray *c21_rhs56 = NULL;
  const mxArray *c21_lhs56 = NULL;
  const mxArray *c21_rhs57 = NULL;
  const mxArray *c21_lhs57 = NULL;
  const mxArray *c21_rhs58 = NULL;
  const mxArray *c21_lhs58 = NULL;
  const mxArray *c21_rhs59 = NULL;
  const mxArray *c21_lhs59 = NULL;
  const mxArray *c21_rhs60 = NULL;
  const mxArray *c21_lhs60 = NULL;
  const mxArray *c21_rhs61 = NULL;
  const mxArray *c21_lhs61 = NULL;
  const mxArray *c21_rhs62 = NULL;
  const mxArray *c21_lhs62 = NULL;
  const mxArray *c21_rhs63 = NULL;
  const mxArray *c21_lhs63 = NULL;
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("norm"), "name", "name", 0);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c21_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 1);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 1);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c21_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 2);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 2);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c21_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 3);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  3);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c21_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 4);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c21_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.blas.xnrm2"),
                  "name", "name", 5);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c21_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 6);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 6);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c21_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p!below_threshold"),
                  "context", "context", 7);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 7);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c21_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 8);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 8);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c21_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 9);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.refblas.xnrm2"), "name", "name", 9);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c21_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 10);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("realmin"), "name", "name",
                  10);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c21_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_realmin"), "name",
                  "name", 11);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1307658444U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c21_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 12);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c21_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 13);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 13);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c21_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 14);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 14);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 14);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c21_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 15);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 15);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 15);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c21_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 16);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 16);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c21_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 17);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("intmax"), "name", "name", 17);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c21_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 18);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c21_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 19);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("abs"), "name", "name", 19);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c21_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 20);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c21_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 21);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c21_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 22);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("mrdivide"), "name", "name",
                  22);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c21_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 23);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 23);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c21_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 24);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("rdivide"), "name", "name",
                  24);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c21_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 25);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c21_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 26);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c21_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_div"), "name", "name",
                  27);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c21_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 28);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c21_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 29);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("cos"), "name", "name", 29);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1343837572U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c21_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 30);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 30);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1286825922U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c21_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 31);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("sin"), "name", "name", 31);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 31);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c21_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 32);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1286825936U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c21_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 33);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eye"), "name", "name", 33);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c21_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 34);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 34);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c21_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 35);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c21_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 36);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("isinf"), "name", "name", 36);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 36);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c21_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 37);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 37);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c21_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 38);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_is_integer_class"),
                  "name", "name", 38);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c21_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 39);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("intmax"), "name", "name", 39);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c21_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 40);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("intmin"), "name", "name", 40);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 40);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c21_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 41);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c21_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 42);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.indexIntRelop"), "name", "name", 42);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c21_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 43);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 43);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c21_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 44);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 44);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c21_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 45);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("intmin"), "name", "name", 45);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 45);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c21_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 46);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 46);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c21_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 47);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("intmax"), "name", "name", 47);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 47);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c21_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 48);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 48);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c21_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 49);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 49);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c21_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 50);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 50);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c21_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 51);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("fn_CrossTensor"), "name",
                  "name", 51);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1450348648U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c21_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_CrossTensor.m"), "context",
                  "context", 52);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 52);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1447321639U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c21_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_CrossTensor.m"), "context",
                  "context", 53);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eye"), "name", "name", 53);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 53);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c21_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_CrossTensor.m"), "context",
                  "context", 54);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 54);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c21_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 55);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("isnan"), "name", "name", 55);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 55);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c21_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 56);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 56);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c21_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 57);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_li_find"), "name",
                  "name", 57);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m"), "resolved",
                  "resolved", 57);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1286825986U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c21_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m"), "context",
                  "context", 58);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 58);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c21_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m!compute_nones"),
                  "context", "context", 59);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 59);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c21_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m!compute_nones"),
                  "context", "context", 60);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 60);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c21_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m!compute_nones"),
                  "context", "context", 61);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 61);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c21_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 62);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 62);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c21_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m"), "context",
                  "context", 63);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 63);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c21_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c21_rhs0);
  sf_mex_destroy(&c21_lhs0);
  sf_mex_destroy(&c21_rhs1);
  sf_mex_destroy(&c21_lhs1);
  sf_mex_destroy(&c21_rhs2);
  sf_mex_destroy(&c21_lhs2);
  sf_mex_destroy(&c21_rhs3);
  sf_mex_destroy(&c21_lhs3);
  sf_mex_destroy(&c21_rhs4);
  sf_mex_destroy(&c21_lhs4);
  sf_mex_destroy(&c21_rhs5);
  sf_mex_destroy(&c21_lhs5);
  sf_mex_destroy(&c21_rhs6);
  sf_mex_destroy(&c21_lhs6);
  sf_mex_destroy(&c21_rhs7);
  sf_mex_destroy(&c21_lhs7);
  sf_mex_destroy(&c21_rhs8);
  sf_mex_destroy(&c21_lhs8);
  sf_mex_destroy(&c21_rhs9);
  sf_mex_destroy(&c21_lhs9);
  sf_mex_destroy(&c21_rhs10);
  sf_mex_destroy(&c21_lhs10);
  sf_mex_destroy(&c21_rhs11);
  sf_mex_destroy(&c21_lhs11);
  sf_mex_destroy(&c21_rhs12);
  sf_mex_destroy(&c21_lhs12);
  sf_mex_destroy(&c21_rhs13);
  sf_mex_destroy(&c21_lhs13);
  sf_mex_destroy(&c21_rhs14);
  sf_mex_destroy(&c21_lhs14);
  sf_mex_destroy(&c21_rhs15);
  sf_mex_destroy(&c21_lhs15);
  sf_mex_destroy(&c21_rhs16);
  sf_mex_destroy(&c21_lhs16);
  sf_mex_destroy(&c21_rhs17);
  sf_mex_destroy(&c21_lhs17);
  sf_mex_destroy(&c21_rhs18);
  sf_mex_destroy(&c21_lhs18);
  sf_mex_destroy(&c21_rhs19);
  sf_mex_destroy(&c21_lhs19);
  sf_mex_destroy(&c21_rhs20);
  sf_mex_destroy(&c21_lhs20);
  sf_mex_destroy(&c21_rhs21);
  sf_mex_destroy(&c21_lhs21);
  sf_mex_destroy(&c21_rhs22);
  sf_mex_destroy(&c21_lhs22);
  sf_mex_destroy(&c21_rhs23);
  sf_mex_destroy(&c21_lhs23);
  sf_mex_destroy(&c21_rhs24);
  sf_mex_destroy(&c21_lhs24);
  sf_mex_destroy(&c21_rhs25);
  sf_mex_destroy(&c21_lhs25);
  sf_mex_destroy(&c21_rhs26);
  sf_mex_destroy(&c21_lhs26);
  sf_mex_destroy(&c21_rhs27);
  sf_mex_destroy(&c21_lhs27);
  sf_mex_destroy(&c21_rhs28);
  sf_mex_destroy(&c21_lhs28);
  sf_mex_destroy(&c21_rhs29);
  sf_mex_destroy(&c21_lhs29);
  sf_mex_destroy(&c21_rhs30);
  sf_mex_destroy(&c21_lhs30);
  sf_mex_destroy(&c21_rhs31);
  sf_mex_destroy(&c21_lhs31);
  sf_mex_destroy(&c21_rhs32);
  sf_mex_destroy(&c21_lhs32);
  sf_mex_destroy(&c21_rhs33);
  sf_mex_destroy(&c21_lhs33);
  sf_mex_destroy(&c21_rhs34);
  sf_mex_destroy(&c21_lhs34);
  sf_mex_destroy(&c21_rhs35);
  sf_mex_destroy(&c21_lhs35);
  sf_mex_destroy(&c21_rhs36);
  sf_mex_destroy(&c21_lhs36);
  sf_mex_destroy(&c21_rhs37);
  sf_mex_destroy(&c21_lhs37);
  sf_mex_destroy(&c21_rhs38);
  sf_mex_destroy(&c21_lhs38);
  sf_mex_destroy(&c21_rhs39);
  sf_mex_destroy(&c21_lhs39);
  sf_mex_destroy(&c21_rhs40);
  sf_mex_destroy(&c21_lhs40);
  sf_mex_destroy(&c21_rhs41);
  sf_mex_destroy(&c21_lhs41);
  sf_mex_destroy(&c21_rhs42);
  sf_mex_destroy(&c21_lhs42);
  sf_mex_destroy(&c21_rhs43);
  sf_mex_destroy(&c21_lhs43);
  sf_mex_destroy(&c21_rhs44);
  sf_mex_destroy(&c21_lhs44);
  sf_mex_destroy(&c21_rhs45);
  sf_mex_destroy(&c21_lhs45);
  sf_mex_destroy(&c21_rhs46);
  sf_mex_destroy(&c21_lhs46);
  sf_mex_destroy(&c21_rhs47);
  sf_mex_destroy(&c21_lhs47);
  sf_mex_destroy(&c21_rhs48);
  sf_mex_destroy(&c21_lhs48);
  sf_mex_destroy(&c21_rhs49);
  sf_mex_destroy(&c21_lhs49);
  sf_mex_destroy(&c21_rhs50);
  sf_mex_destroy(&c21_lhs50);
  sf_mex_destroy(&c21_rhs51);
  sf_mex_destroy(&c21_lhs51);
  sf_mex_destroy(&c21_rhs52);
  sf_mex_destroy(&c21_lhs52);
  sf_mex_destroy(&c21_rhs53);
  sf_mex_destroy(&c21_lhs53);
  sf_mex_destroy(&c21_rhs54);
  sf_mex_destroy(&c21_lhs54);
  sf_mex_destroy(&c21_rhs55);
  sf_mex_destroy(&c21_lhs55);
  sf_mex_destroy(&c21_rhs56);
  sf_mex_destroy(&c21_lhs56);
  sf_mex_destroy(&c21_rhs57);
  sf_mex_destroy(&c21_lhs57);
  sf_mex_destroy(&c21_rhs58);
  sf_mex_destroy(&c21_lhs58);
  sf_mex_destroy(&c21_rhs59);
  sf_mex_destroy(&c21_lhs59);
  sf_mex_destroy(&c21_rhs60);
  sf_mex_destroy(&c21_lhs60);
  sf_mex_destroy(&c21_rhs61);
  sf_mex_destroy(&c21_lhs61);
  sf_mex_destroy(&c21_rhs62);
  sf_mex_destroy(&c21_lhs62);
  sf_mex_destroy(&c21_rhs63);
  sf_mex_destroy(&c21_lhs63);
}

static const mxArray *c21_emlrt_marshallOut(const char * c21_u)
{
  const mxArray *c21_y = NULL;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", c21_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c21_u)), false);
  return c21_y;
}

static const mxArray *c21_b_emlrt_marshallOut(const uint32_T c21_u)
{
  const mxArray *c21_y = NULL;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", &c21_u, 7, 0U, 0U, 0U, 0), false);
  return c21_y;
}

static void c21_b_info_helper(const mxArray **c21_info)
{
  const mxArray *c21_rhs64 = NULL;
  const mxArray *c21_lhs64 = NULL;
  const mxArray *c21_rhs65 = NULL;
  const mxArray *c21_lhs65 = NULL;
  const mxArray *c21_rhs66 = NULL;
  const mxArray *c21_lhs66 = NULL;
  const mxArray *c21_rhs67 = NULL;
  const mxArray *c21_lhs67 = NULL;
  const mxArray *c21_rhs68 = NULL;
  const mxArray *c21_lhs68 = NULL;
  const mxArray *c21_rhs69 = NULL;
  const mxArray *c21_lhs69 = NULL;
  const mxArray *c21_rhs70 = NULL;
  const mxArray *c21_lhs70 = NULL;
  const mxArray *c21_rhs71 = NULL;
  const mxArray *c21_lhs71 = NULL;
  const mxArray *c21_rhs72 = NULL;
  const mxArray *c21_lhs72 = NULL;
  const mxArray *c21_rhs73 = NULL;
  const mxArray *c21_lhs73 = NULL;
  const mxArray *c21_rhs74 = NULL;
  const mxArray *c21_lhs74 = NULL;
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m"), "context",
                  "context", 64);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 64);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c21_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 65);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 65);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c21_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 66);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 66);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 66);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c21_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 67);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 67);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 67);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c21_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 68);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  68);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c21_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 69);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 69);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 69);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c21_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 70);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 70);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c21_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 71);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 71);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c21_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 72);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 72);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c21_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 73);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 73);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c21_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 74);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.refblas.xgemm"), "name", "name", 74);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c21_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c21_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs74), "lhs", "lhs",
                  74);
  sf_mex_destroy(&c21_rhs64);
  sf_mex_destroy(&c21_lhs64);
  sf_mex_destroy(&c21_rhs65);
  sf_mex_destroy(&c21_lhs65);
  sf_mex_destroy(&c21_rhs66);
  sf_mex_destroy(&c21_lhs66);
  sf_mex_destroy(&c21_rhs67);
  sf_mex_destroy(&c21_lhs67);
  sf_mex_destroy(&c21_rhs68);
  sf_mex_destroy(&c21_lhs68);
  sf_mex_destroy(&c21_rhs69);
  sf_mex_destroy(&c21_lhs69);
  sf_mex_destroy(&c21_rhs70);
  sf_mex_destroy(&c21_lhs70);
  sf_mex_destroy(&c21_rhs71);
  sf_mex_destroy(&c21_lhs71);
  sf_mex_destroy(&c21_rhs72);
  sf_mex_destroy(&c21_lhs72);
  sf_mex_destroy(&c21_rhs73);
  sf_mex_destroy(&c21_lhs73);
  sf_mex_destroy(&c21_rhs74);
  sf_mex_destroy(&c21_lhs74);
}

static real_T c21_norm(SFc21_Model_01InstanceStruct *chartInstance, real_T
  c21_x[3])
{
  real_T c21_y;
  real_T c21_scale;
  int32_T c21_k;
  int32_T c21_b_k;
  real_T c21_b_x;
  real_T c21_c_x;
  real_T c21_absxk;
  real_T c21_t;
  c21_threshold(chartInstance);
  c21_y = 0.0;
  c21_realmin(chartInstance);
  c21_scale = 2.2250738585072014E-308;
  for (c21_k = 1; c21_k < 4; c21_k++) {
    c21_b_k = c21_k - 1;
    c21_b_x = c21_x[c21_b_k];
    c21_c_x = c21_b_x;
    c21_absxk = muDoubleScalarAbs(c21_c_x);
    if (c21_absxk > c21_scale) {
      c21_t = c21_scale / c21_absxk;
      c21_y = 1.0 + c21_y * c21_t * c21_t;
      c21_scale = c21_absxk;
    } else {
      c21_t = c21_absxk / c21_scale;
      c21_y += c21_t * c21_t;
    }
  }

  return c21_scale * muDoubleScalarSqrt(c21_y);
}

static void c21_threshold(SFc21_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c21_realmin(SFc21_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c21_eml_scalar_eg(SFc21_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c21_eml_xgemm(SFc21_Model_01InstanceStruct *chartInstance, real_T
  c21_A[16], real_T c21_B[4], real_T c21_C[4], real_T c21_b_C[4])
{
  int32_T c21_i70;
  int32_T c21_i71;
  real_T c21_b_A[16];
  int32_T c21_i72;
  real_T c21_b_B[4];
  for (c21_i70 = 0; c21_i70 < 4; c21_i70++) {
    c21_b_C[c21_i70] = c21_C[c21_i70];
  }

  for (c21_i71 = 0; c21_i71 < 16; c21_i71++) {
    c21_b_A[c21_i71] = c21_A[c21_i71];
  }

  for (c21_i72 = 0; c21_i72 < 4; c21_i72++) {
    c21_b_B[c21_i72] = c21_B[c21_i72];
  }

  c21_b_eml_xgemm(chartInstance, c21_b_A, c21_b_B, c21_b_C);
}

static const mxArray *c21_f_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData)
{
  const mxArray *c21_mxArrayOutData = NULL;
  int32_T c21_u;
  const mxArray *c21_y = NULL;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_mxArrayOutData = NULL;
  c21_u = *(int32_T *)c21_inData;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", &c21_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c21_mxArrayOutData, c21_y, false);
  return c21_mxArrayOutData;
}

static int32_T c21_g_emlrt_marshallIn(SFc21_Model_01InstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId)
{
  int32_T c21_y;
  int32_T c21_i73;
  (void)chartInstance;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), &c21_i73, 1, 6, 0U, 0, 0U, 0);
  c21_y = c21_i73;
  sf_mex_destroy(&c21_u);
  return c21_y;
}

static void c21_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData)
{
  const mxArray *c21_b_sfEvent;
  const char_T *c21_identifier;
  emlrtMsgIdentifier c21_thisId;
  int32_T c21_y;
  SFc21_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc21_Model_01InstanceStruct *)chartInstanceVoid;
  c21_b_sfEvent = sf_mex_dup(c21_mxArrayInData);
  c21_identifier = c21_varName;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_y = c21_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_b_sfEvent),
    &c21_thisId);
  sf_mex_destroy(&c21_b_sfEvent);
  *(int32_T *)c21_outData = c21_y;
  sf_mex_destroy(&c21_mxArrayInData);
}

static uint8_T c21_h_emlrt_marshallIn(SFc21_Model_01InstanceStruct
  *chartInstance, const mxArray *c21_b_is_active_c21_Model_01, const char_T
  *c21_identifier)
{
  uint8_T c21_y;
  emlrtMsgIdentifier c21_thisId;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_y = c21_i_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c21_b_is_active_c21_Model_01), &c21_thisId);
  sf_mex_destroy(&c21_b_is_active_c21_Model_01);
  return c21_y;
}

static uint8_T c21_i_emlrt_marshallIn(SFc21_Model_01InstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId)
{
  uint8_T c21_y;
  uint8_T c21_u0;
  (void)chartInstance;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), &c21_u0, 1, 3, 0U, 0, 0U, 0);
  c21_y = c21_u0;
  sf_mex_destroy(&c21_u);
  return c21_y;
}

static void c21_b_eml_xgemm(SFc21_Model_01InstanceStruct *chartInstance, real_T
  c21_A[16], real_T c21_B[4], real_T c21_C[4])
{
  int32_T c21_i74;
  int32_T c21_i75;
  int32_T c21_i76;
  (void)chartInstance;
  for (c21_i74 = 0; c21_i74 < 4; c21_i74++) {
    c21_C[c21_i74] = 0.0;
    c21_i75 = 0;
    for (c21_i76 = 0; c21_i76 < 4; c21_i76++) {
      c21_C[c21_i74] += c21_A[c21_i75 + c21_i74] * c21_B[c21_i76];
      c21_i75 += 4;
    }
  }
}

static void init_dsm_address_info(SFc21_Model_01InstanceStruct *chartInstance)
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

void sf_c21_Model_01_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2930957761U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3519179060U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1539375366U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(17274774U);
}

mxArray *sf_c21_Model_01_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("M0FAt4GVCmVCI61eKyXeBF");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

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
      pr[0] = (double)(1);
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
      pr[0] = (double)(4);
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

mxArray *sf_c21_Model_01_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c21_Model_01_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c21_Model_01(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"q_nominal\",},{M[8],M[0],T\"is_active_c21_Model_01\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c21_Model_01_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc21_Model_01InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc21_Model_01InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_01MachineNumber_,
           21,
           1,
           1,
           0,
           5,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"omega");
          _SFD_SET_DATA_PROPS(1,2,0,1,"q_nominal");
          _SFD_SET_DATA_PROPS(2,1,1,0,"t_delta");
          _SFD_SET_DATA_PROPS(3,1,1,0,"n");
          _SFD_SET_DATA_PROPS(4,1,1,0,"q");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,868);
        _SFD_CV_INIT_EML_IF(0,1,0,122,142,-1,863);
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
            1.0,0,0,(MexFcnForType)c21_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c21_sf_marshallOut,(MexInFcnForType)
            c21_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c21_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c21_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c21_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c21_t_delta;
          real_T *c21_n;
          real_T (*c21_omega)[3];
          real_T (*c21_q_nominal)[4];
          real_T (*c21_q)[4];
          c21_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 3);
          c21_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c21_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c21_q_nominal = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S,
            1);
          c21_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c21_omega);
          _SFD_SET_DATA_VALUE_PTR(1U, *c21_q_nominal);
          _SFD_SET_DATA_VALUE_PTR(2U, c21_t_delta);
          _SFD_SET_DATA_VALUE_PTR(3U, c21_n);
          _SFD_SET_DATA_VALUE_PTR(4U, *c21_q);
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
  return "ptcXT9rvDNI6oi3hVogNqF";
}

static void sf_opaque_initialize_c21_Model_01(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc21_Model_01InstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c21_Model_01((SFc21_Model_01InstanceStruct*)
    chartInstanceVar);
  initialize_c21_Model_01((SFc21_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c21_Model_01(void *chartInstanceVar)
{
  enable_c21_Model_01((SFc21_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c21_Model_01(void *chartInstanceVar)
{
  disable_c21_Model_01((SFc21_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c21_Model_01(void *chartInstanceVar)
{
  sf_gateway_c21_Model_01((SFc21_Model_01InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c21_Model_01(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c21_Model_01((SFc21_Model_01InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c21_Model_01();/* state var info */
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

extern void sf_internal_set_sim_state_c21_Model_01(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c21_Model_01();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c21_Model_01((SFc21_Model_01InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c21_Model_01(SimStruct* S)
{
  return sf_internal_get_sim_state_c21_Model_01(S);
}

static void sf_opaque_set_sim_state_c21_Model_01(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c21_Model_01(S, st);
}

static void sf_opaque_terminate_c21_Model_01(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc21_Model_01InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_01_optimization_info();
    }

    finalize_c21_Model_01((SFc21_Model_01InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc21_Model_01((SFc21_Model_01InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c21_Model_01(SimStruct *S)
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
    initialize_params_c21_Model_01((SFc21_Model_01InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c21_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_01_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      21);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,21,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,21,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,21);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,21,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,21,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,21);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1775602629U));
  ssSetChecksum1(S,(347188542U));
  ssSetChecksum2(S,(905647487U));
  ssSetChecksum3(S,(293222638U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c21_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c21_Model_01(SimStruct *S)
{
  SFc21_Model_01InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc21_Model_01InstanceStruct *)utMalloc(sizeof
    (SFc21_Model_01InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc21_Model_01InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c21_Model_01;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c21_Model_01;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c21_Model_01;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c21_Model_01;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c21_Model_01;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c21_Model_01;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c21_Model_01;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c21_Model_01;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c21_Model_01;
  chartInstance->chartInfo.mdlStart = mdlStart_c21_Model_01;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c21_Model_01;
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

void c21_Model_01_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c21_Model_01(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c21_Model_01(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c21_Model_01(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c21_Model_01_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_01_sfun.h"
#include "c25_Model_01.h"
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
static const char * c25_debug_family_names[12] = { "mu", "rc", "ita_star",
  "q_star", "v", "nargin", "nargout", "signal", "rho_c", "q_k", "ita_k", "zk" };

static const char * c25_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

static const char * c25_c_debug_family_names[8] = { "q_v", "q_0", "cross_q_v",
  "nargin", "nargout", "q", "flag", "CrossTensor" };

/* Function Declarations */
static void initialize_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance);
static void initialize_params_c25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance);
static void enable_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance);
static void disable_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance);
static void c25_update_debugger_state_c25_Model_01(SFc25_Model_01InstanceStruct *
  chartInstance);
static const mxArray *get_sim_state_c25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance);
static void set_sim_state_c25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance, const mxArray *c25_st);
static void finalize_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance);
static void sf_gateway_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance);
static void c25_chartstep_c25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance);
static void initSimStructsc25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance);
static void c25_fn_CrossTensor(SFc25_Model_01InstanceStruct *chartInstance,
  real_T c25_q[4], real_T c25_flag, real_T c25_CrossTensor[16]);
static void init_script_number_translation(uint32_T c25_machineNumber, uint32_T
  c25_chartNumber, uint32_T c25_instanceNumber);
static const mxArray *c25_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData);
static void c25_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_zk, const char_T *c25_identifier, real_T c25_y[6]);
static void c25_b_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[6]);
static void c25_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData);
static const mxArray *c25_b_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData);
static const mxArray *c25_c_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData);
static const mxArray *c25_d_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData);
static const mxArray *c25_e_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData);
static real_T c25_c_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId);
static void c25_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData);
static void c25_d_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[4]);
static void c25_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData);
static void c25_e_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[3]);
static void c25_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData);
static const mxArray *c25_f_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData);
static void c25_f_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[9]);
static void c25_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData);
static const mxArray *c25_g_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData);
static void c25_g_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[16]);
static void c25_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData);
static void c25_info_helper(const mxArray **c25_info);
static const mxArray *c25_emlrt_marshallOut(const char * c25_u);
static const mxArray *c25_b_emlrt_marshallOut(const uint32_T c25_u);
static void c25_eml_scalar_eg(SFc25_Model_01InstanceStruct *chartInstance);
static void c25_threshold(SFc25_Model_01InstanceStruct *chartInstance);
static void c25_b_eml_scalar_eg(SFc25_Model_01InstanceStruct *chartInstance);
static const mxArray *c25_h_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData);
static int32_T c25_h_emlrt_marshallIn(SFc25_Model_01InstanceStruct
  *chartInstance, const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId);
static void c25_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData);
static uint8_T c25_i_emlrt_marshallIn(SFc25_Model_01InstanceStruct
  *chartInstance, const mxArray *c25_b_is_active_c25_Model_01, const char_T
  *c25_identifier);
static uint8_T c25_j_emlrt_marshallIn(SFc25_Model_01InstanceStruct
  *chartInstance, const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId);
static void init_dsm_address_info(SFc25_Model_01InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance)
{
  chartInstance->c25_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c25_is_active_c25_Model_01 = 0U;
}

static void initialize_params_c25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c25_update_debugger_state_c25_Model_01(SFc25_Model_01InstanceStruct *
  chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance)
{
  const mxArray *c25_st;
  const mxArray *c25_y = NULL;
  int32_T c25_i0;
  real_T c25_u[6];
  const mxArray *c25_b_y = NULL;
  uint8_T c25_hoistedGlobal;
  uint8_T c25_b_u;
  const mxArray *c25_c_y = NULL;
  real_T (*c25_zk)[6];
  c25_zk = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c25_st = NULL;
  c25_st = NULL;
  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_createcellmatrix(2, 1), false);
  for (c25_i0 = 0; c25_i0 < 6; c25_i0++) {
    c25_u[c25_i0] = (*c25_zk)[c25_i0];
  }

  c25_b_y = NULL;
  sf_mex_assign(&c25_b_y, sf_mex_create("y", c25_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_setcell(c25_y, 0, c25_b_y);
  c25_hoistedGlobal = chartInstance->c25_is_active_c25_Model_01;
  c25_b_u = c25_hoistedGlobal;
  c25_c_y = NULL;
  sf_mex_assign(&c25_c_y, sf_mex_create("y", &c25_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c25_y, 1, c25_c_y);
  sf_mex_assign(&c25_st, c25_y, false);
  return c25_st;
}

static void set_sim_state_c25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance, const mxArray *c25_st)
{
  const mxArray *c25_u;
  real_T c25_dv0[6];
  int32_T c25_i1;
  real_T (*c25_zk)[6];
  c25_zk = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c25_doneDoubleBufferReInit = true;
  c25_u = sf_mex_dup(c25_st);
  c25_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c25_u, 0)), "zk",
                       c25_dv0);
  for (c25_i1 = 0; c25_i1 < 6; c25_i1++) {
    (*c25_zk)[c25_i1] = c25_dv0[c25_i1];
  }

  chartInstance->c25_is_active_c25_Model_01 = c25_i_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c25_u, 1)),
     "is_active_c25_Model_01");
  sf_mex_destroy(&c25_u);
  c25_update_debugger_state_c25_Model_01(chartInstance);
  sf_mex_destroy(&c25_st);
}

static void finalize_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c25_Model_01(SFc25_Model_01InstanceStruct *chartInstance)
{
  int32_T c25_i2;
  int32_T c25_i3;
  int32_T c25_i4;
  int32_T c25_i5;
  int32_T c25_i6;
  real_T (*c25_ita_k)[4];
  real_T (*c25_q_k)[4];
  real_T (*c25_rho_c)[3];
  real_T (*c25_zk)[6];
  real_T (*c25_signal)[7];
  c25_ita_k = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 3);
  c25_q_k = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 2);
  c25_rho_c = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c25_zk = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c25_signal = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 24U, chartInstance->c25_sfEvent);
  for (c25_i2 = 0; c25_i2 < 7; c25_i2++) {
    _SFD_DATA_RANGE_CHECK((*c25_signal)[c25_i2], 0U);
  }

  chartInstance->c25_sfEvent = CALL_EVENT;
  c25_chartstep_c25_Model_01(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_01MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c25_i3 = 0; c25_i3 < 6; c25_i3++) {
    _SFD_DATA_RANGE_CHECK((*c25_zk)[c25_i3], 1U);
  }

  for (c25_i4 = 0; c25_i4 < 3; c25_i4++) {
    _SFD_DATA_RANGE_CHECK((*c25_rho_c)[c25_i4], 2U);
  }

  for (c25_i5 = 0; c25_i5 < 4; c25_i5++) {
    _SFD_DATA_RANGE_CHECK((*c25_q_k)[c25_i5], 3U);
  }

  for (c25_i6 = 0; c25_i6 < 4; c25_i6++) {
    _SFD_DATA_RANGE_CHECK((*c25_ita_k)[c25_i6], 4U);
  }
}

static void c25_chartstep_c25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance)
{
  int32_T c25_i7;
  real_T c25_signal[7];
  int32_T c25_i8;
  real_T c25_rho_c[3];
  int32_T c25_i9;
  real_T c25_q_k[4];
  int32_T c25_i10;
  real_T c25_ita_k[4];
  uint32_T c25_debug_family_var_map[12];
  real_T c25_mu[4];
  real_T c25_rc[3];
  real_T c25_ita_star[4];
  real_T c25_q_star[4];
  real_T c25_v[4];
  real_T c25_nargin = 4.0;
  real_T c25_nargout = 1.0;
  real_T c25_zk[6];
  int32_T c25_i11;
  int32_T c25_i12;
  int32_T c25_i13;
  int32_T c25_i14;
  int32_T c25_i15;
  int32_T c25_i16;
  int32_T c25_i17;
  real_T c25_b_ita_star[4];
  real_T c25_a[16];
  int32_T c25_i18;
  real_T c25_b_mu[4];
  real_T c25_b[16];
  int32_T c25_i19;
  int32_T c25_i20;
  int32_T c25_i21;
  real_T c25_y[16];
  int32_T c25_i22;
  int32_T c25_i23;
  int32_T c25_i24;
  real_T c25_b_b[4];
  int32_T c25_i25;
  int32_T c25_i26;
  int32_T c25_i27;
  real_T c25_C[4];
  int32_T c25_i28;
  int32_T c25_i29;
  int32_T c25_i30;
  int32_T c25_i31;
  int32_T c25_i32;
  int32_T c25_i33;
  int32_T c25_i34;
  int32_T c25_i35;
  real_T (*c25_b_zk)[6];
  real_T (*c25_b_ita_k)[4];
  real_T (*c25_b_q_k)[4];
  real_T (*c25_b_rho_c)[3];
  real_T (*c25_b_signal)[7];
  c25_b_ita_k = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 3);
  c25_b_q_k = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 2);
  c25_b_rho_c = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c25_b_zk = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c25_b_signal = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 24U, chartInstance->c25_sfEvent);
  for (c25_i7 = 0; c25_i7 < 7; c25_i7++) {
    c25_signal[c25_i7] = (*c25_b_signal)[c25_i7];
  }

  for (c25_i8 = 0; c25_i8 < 3; c25_i8++) {
    c25_rho_c[c25_i8] = (*c25_b_rho_c)[c25_i8];
  }

  for (c25_i9 = 0; c25_i9 < 4; c25_i9++) {
    c25_q_k[c25_i9] = (*c25_b_q_k)[c25_i9];
  }

  for (c25_i10 = 0; c25_i10 < 4; c25_i10++) {
    c25_ita_k[c25_i10] = (*c25_b_ita_k)[c25_i10];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 12U, 12U, c25_debug_family_names,
    c25_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_mu, 0U, c25_b_sf_marshallOut,
    c25_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_rc, 1U, c25_c_sf_marshallOut,
    c25_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_ita_star, 2U, c25_b_sf_marshallOut,
    c25_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_q_star, 3U, c25_b_sf_marshallOut,
    c25_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_v, 4U, c25_b_sf_marshallOut,
    c25_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c25_nargin, 5U, c25_e_sf_marshallOut,
    c25_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c25_nargout, 6U, c25_e_sf_marshallOut,
    c25_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c25_signal, 7U, c25_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c25_rho_c, 8U, c25_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c25_q_k, 9U, c25_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c25_ita_k, 10U, c25_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_zk, 11U, c25_sf_marshallOut,
    c25_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c25_sfEvent, 3);
  for (c25_i11 = 0; c25_i11 < 6; c25_i11++) {
    c25_zk[c25_i11] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c25_sfEvent, 4);
  for (c25_i12 = 0; c25_i12 < 4; c25_i12++) {
    c25_mu[c25_i12] = c25_signal[c25_i12];
  }

  _SFD_EML_CALL(0U, chartInstance->c25_sfEvent, 5);
  for (c25_i13 = 0; c25_i13 < 3; c25_i13++) {
    c25_rc[c25_i13] = c25_signal[c25_i13 + 4];
  }

  _SFD_EML_CALL(0U, chartInstance->c25_sfEvent, 7);
  for (c25_i14 = 0; c25_i14 < 3; c25_i14++) {
    c25_ita_star[c25_i14] = -c25_ita_k[c25_i14];
  }

  c25_ita_star[3] = c25_ita_k[3];
  _SFD_EML_CALL(0U, chartInstance->c25_sfEvent, 8);
  for (c25_i15 = 0; c25_i15 < 3; c25_i15++) {
    c25_q_star[c25_i15] = -c25_q_k[c25_i15];
  }

  c25_q_star[3] = c25_q_k[3];
  _SFD_EML_CALL(0U, chartInstance->c25_sfEvent, 9);
  for (c25_i16 = 0; c25_i16 < 3; c25_i16++) {
    c25_zk[c25_i16] = c25_rc[c25_i16] - c25_rho_c[c25_i16];
  }

  _SFD_EML_CALL(0U, chartInstance->c25_sfEvent, 10);
  for (c25_i17 = 0; c25_i17 < 4; c25_i17++) {
    c25_b_ita_star[c25_i17] = c25_ita_star[c25_i17];
  }

  c25_fn_CrossTensor(chartInstance, c25_b_ita_star, 0.0, c25_a);
  for (c25_i18 = 0; c25_i18 < 4; c25_i18++) {
    c25_b_mu[c25_i18] = c25_mu[c25_i18];
  }

  c25_fn_CrossTensor(chartInstance, c25_b_mu, 0.0, c25_b);
  c25_eml_scalar_eg(chartInstance);
  c25_eml_scalar_eg(chartInstance);
  c25_threshold(chartInstance);
  for (c25_i19 = 0; c25_i19 < 4; c25_i19++) {
    c25_i20 = 0;
    for (c25_i21 = 0; c25_i21 < 4; c25_i21++) {
      c25_y[c25_i20 + c25_i19] = 0.0;
      c25_i22 = 0;
      for (c25_i23 = 0; c25_i23 < 4; c25_i23++) {
        c25_y[c25_i20 + c25_i19] += c25_a[c25_i22 + c25_i19] * c25_b[c25_i23 +
          c25_i20];
        c25_i22 += 4;
      }

      c25_i20 += 4;
    }
  }

  for (c25_i24 = 0; c25_i24 < 4; c25_i24++) {
    c25_b_b[c25_i24] = c25_q_star[c25_i24];
  }

  c25_b_eml_scalar_eg(chartInstance);
  c25_b_eml_scalar_eg(chartInstance);
  for (c25_i25 = 0; c25_i25 < 4; c25_i25++) {
    c25_v[c25_i25] = 0.0;
  }

  for (c25_i26 = 0; c25_i26 < 4; c25_i26++) {
    c25_v[c25_i26] = 0.0;
  }

  for (c25_i27 = 0; c25_i27 < 4; c25_i27++) {
    c25_C[c25_i27] = c25_v[c25_i27];
  }

  for (c25_i28 = 0; c25_i28 < 4; c25_i28++) {
    c25_v[c25_i28] = c25_C[c25_i28];
  }

  c25_threshold(chartInstance);
  for (c25_i29 = 0; c25_i29 < 4; c25_i29++) {
    c25_C[c25_i29] = c25_v[c25_i29];
  }

  for (c25_i30 = 0; c25_i30 < 4; c25_i30++) {
    c25_v[c25_i30] = c25_C[c25_i30];
  }

  for (c25_i31 = 0; c25_i31 < 4; c25_i31++) {
    c25_v[c25_i31] = 0.0;
    c25_i32 = 0;
    for (c25_i33 = 0; c25_i33 < 4; c25_i33++) {
      c25_v[c25_i31] += c25_y[c25_i32 + c25_i31] * c25_b_b[c25_i33];
      c25_i32 += 4;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c25_sfEvent, 11);
  for (c25_i34 = 0; c25_i34 < 3; c25_i34++) {
    c25_zk[c25_i34 + 3] = c25_v[c25_i34];
  }

  _SFD_EML_CALL(0U, chartInstance->c25_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  for (c25_i35 = 0; c25_i35 < 6; c25_i35++) {
    (*c25_b_zk)[c25_i35] = c25_zk[c25_i35];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 24U, chartInstance->c25_sfEvent);
}

static void initSimStructsc25_Model_01(SFc25_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c25_fn_CrossTensor(SFc25_Model_01InstanceStruct *chartInstance,
  real_T c25_q[4], real_T c25_flag, real_T c25_CrossTensor[16])
{
  uint32_T c25_debug_family_var_map[8];
  real_T c25_q_v[3];
  real_T c25_q_0;
  real_T c25_cross_q_v[9];
  real_T c25_nargin = 2.0;
  real_T c25_nargout = 1.0;
  int32_T c25_i36;
  int32_T c25_i37;
  int32_T c25_i38;
  real_T c25_v[3];
  uint32_T c25_b_debug_family_var_map[4];
  real_T c25_b_nargin = 1.0;
  real_T c25_b_nargout = 1.0;
  int32_T c25_i39;
  real_T c25_a;
  int32_T c25_i40;
  static real_T c25_b[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c25_y[9];
  int32_T c25_i41;
  int32_T c25_i42;
  int32_T c25_i43;
  int32_T c25_i44;
  int32_T c25_i45;
  int32_T c25_i46;
  int32_T c25_i47;
  real_T c25_b_a;
  int32_T c25_i48;
  int32_T c25_i49;
  int32_T c25_i50;
  int32_T c25_i51;
  int32_T c25_i52;
  int32_T c25_i53;
  int32_T c25_i54;
  int32_T c25_i55;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c25_c_debug_family_names,
    c25_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_q_v, 0U, c25_c_sf_marshallOut,
    c25_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c25_q_0, 1U, c25_e_sf_marshallOut,
    c25_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_cross_q_v, 2U, c25_f_sf_marshallOut,
    c25_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c25_nargin, 3U, c25_e_sf_marshallOut,
    c25_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c25_nargout, 4U, c25_e_sf_marshallOut,
    c25_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_q, 5U, c25_b_sf_marshallOut,
    c25_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c25_flag, 6U, c25_e_sf_marshallOut,
    c25_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_CrossTensor, 7U, c25_g_sf_marshallOut,
    c25_f_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c25_sfEvent, 4);
  for (c25_i36 = 0; c25_i36 < 16; c25_i36++) {
    c25_CrossTensor[c25_i36] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c25_sfEvent, 5);
  for (c25_i37 = 0; c25_i37 < 3; c25_i37++) {
    c25_q_v[c25_i37] = c25_q[c25_i37];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c25_sfEvent, 6);
  c25_q_0 = c25_q[3];
  _SFD_SCRIPT_CALL(0U, chartInstance->c25_sfEvent, 7);
  for (c25_i38 = 0; c25_i38 < 3; c25_i38++) {
    c25_v[c25_i38] = c25_q_v[c25_i38];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c25_b_debug_family_names,
    c25_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c25_b_nargin, 0U, c25_e_sf_marshallOut,
    c25_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c25_b_nargout, 1U, c25_e_sf_marshallOut,
    c25_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_v, 2U, c25_c_sf_marshallOut,
    c25_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c25_cross_q_v, 3U, c25_f_sf_marshallOut,
    c25_e_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 2);
  for (c25_i39 = 0; c25_i39 < 9; c25_i39++) {
    c25_cross_q_v[c25_i39] = 0.0;
  }

  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 3);
  c25_cross_q_v[0] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 4);
  c25_cross_q_v[4] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 5);
  c25_cross_q_v[8] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 6);
  c25_cross_q_v[3] = -c25_v[2];
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 7);
  c25_cross_q_v[6] = c25_v[1];
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 8);
  c25_cross_q_v[7] = -c25_v[0];
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 9);
  c25_cross_q_v[1] = c25_v[2];
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 10);
  c25_cross_q_v[2] = -c25_v[1];
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, 11);
  c25_cross_q_v[5] = c25_v[0];
  _SFD_SCRIPT_CALL(1U, chartInstance->c25_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_SCRIPT_CALL(0U, chartInstance->c25_sfEvent, 8);
  switch ((int32_T)_SFD_INTEGER_CHECK("flag", c25_flag)) {
   case 0:
    CV_SCRIPT_SWITCH(0, 0, 1);
    _SFD_SCRIPT_CALL(0U, chartInstance->c25_sfEvent, 10);
    c25_a = c25_q_0;
    for (c25_i40 = 0; c25_i40 < 9; c25_i40++) {
      c25_y[c25_i40] = c25_a * c25_b[c25_i40];
    }

    c25_i41 = 0;
    c25_i42 = 0;
    for (c25_i43 = 0; c25_i43 < 3; c25_i43++) {
      for (c25_i44 = 0; c25_i44 < 3; c25_i44++) {
        c25_CrossTensor[c25_i44 + c25_i41] = -c25_cross_q_v[c25_i44 + c25_i42] +
          c25_y[c25_i44 + c25_i42];
      }

      c25_i41 += 4;
      c25_i42 += 3;
    }

    for (c25_i45 = 0; c25_i45 < 3; c25_i45++) {
      c25_CrossTensor[c25_i45 + 12] = c25_q_v[c25_i45];
    }

    c25_i46 = 0;
    for (c25_i47 = 0; c25_i47 < 3; c25_i47++) {
      c25_CrossTensor[c25_i46 + 3] = -c25_q_v[c25_i47];
      c25_i46 += 4;
    }

    c25_CrossTensor[15] = c25_q_0;
    break;

   case 1:
    CV_SCRIPT_SWITCH(0, 0, 2);
    _SFD_SCRIPT_CALL(0U, chartInstance->c25_sfEvent, 12);
    c25_b_a = c25_q_0;
    for (c25_i48 = 0; c25_i48 < 9; c25_i48++) {
      c25_y[c25_i48] = c25_b_a * c25_b[c25_i48];
    }

    c25_i49 = 0;
    c25_i50 = 0;
    for (c25_i51 = 0; c25_i51 < 3; c25_i51++) {
      for (c25_i52 = 0; c25_i52 < 3; c25_i52++) {
        c25_CrossTensor[c25_i52 + c25_i49] = c25_cross_q_v[c25_i52 + c25_i50] +
          c25_y[c25_i52 + c25_i50];
      }

      c25_i49 += 4;
      c25_i50 += 3;
    }

    for (c25_i53 = 0; c25_i53 < 3; c25_i53++) {
      c25_CrossTensor[c25_i53 + 12] = c25_q_v[c25_i53];
    }

    c25_i54 = 0;
    for (c25_i55 = 0; c25_i55 < 3; c25_i55++) {
      c25_CrossTensor[c25_i54 + 3] = -c25_q_v[c25_i55];
      c25_i54 += 4;
    }

    c25_CrossTensor[15] = c25_q_0;
    break;

   default:
    CV_SCRIPT_SWITCH(0, 0, 0);
    break;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c25_sfEvent, -12);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c25_machineNumber, uint32_T
  c25_chartNumber, uint32_T c25_instanceNumber)
{
  (void)c25_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c25_chartNumber, c25_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg-2\\Documents\\MATLAB\\Model_01\\fn_CrossTensor.m"));
  _SFD_SCRIPT_TRANSLATION(c25_chartNumber, c25_instanceNumber, 1U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg-2\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c25_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData)
{
  const mxArray *c25_mxArrayOutData = NULL;
  int32_T c25_i56;
  real_T c25_b_inData[6];
  int32_T c25_i57;
  real_T c25_u[6];
  const mxArray *c25_y = NULL;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_mxArrayOutData = NULL;
  for (c25_i56 = 0; c25_i56 < 6; c25_i56++) {
    c25_b_inData[c25_i56] = (*(real_T (*)[6])c25_inData)[c25_i56];
  }

  for (c25_i57 = 0; c25_i57 < 6; c25_i57++) {
    c25_u[c25_i57] = c25_b_inData[c25_i57];
  }

  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", c25_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_assign(&c25_mxArrayOutData, c25_y, false);
  return c25_mxArrayOutData;
}

static void c25_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_zk, const char_T *c25_identifier, real_T c25_y[6])
{
  emlrtMsgIdentifier c25_thisId;
  c25_thisId.fIdentifier = c25_identifier;
  c25_thisId.fParent = NULL;
  c25_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c25_zk), &c25_thisId, c25_y);
  sf_mex_destroy(&c25_zk);
}

static void c25_b_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[6])
{
  real_T c25_dv1[6];
  int32_T c25_i58;
  (void)chartInstance;
  sf_mex_import(c25_parentId, sf_mex_dup(c25_u), c25_dv1, 1, 0, 0U, 1, 0U, 1, 6);
  for (c25_i58 = 0; c25_i58 < 6; c25_i58++) {
    c25_y[c25_i58] = c25_dv1[c25_i58];
  }

  sf_mex_destroy(&c25_u);
}

static void c25_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData)
{
  const mxArray *c25_zk;
  const char_T *c25_identifier;
  emlrtMsgIdentifier c25_thisId;
  real_T c25_y[6];
  int32_T c25_i59;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_zk = sf_mex_dup(c25_mxArrayInData);
  c25_identifier = c25_varName;
  c25_thisId.fIdentifier = c25_identifier;
  c25_thisId.fParent = NULL;
  c25_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c25_zk), &c25_thisId, c25_y);
  sf_mex_destroy(&c25_zk);
  for (c25_i59 = 0; c25_i59 < 6; c25_i59++) {
    (*(real_T (*)[6])c25_outData)[c25_i59] = c25_y[c25_i59];
  }

  sf_mex_destroy(&c25_mxArrayInData);
}

static const mxArray *c25_b_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData)
{
  const mxArray *c25_mxArrayOutData = NULL;
  int32_T c25_i60;
  real_T c25_b_inData[4];
  int32_T c25_i61;
  real_T c25_u[4];
  const mxArray *c25_y = NULL;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_mxArrayOutData = NULL;
  for (c25_i60 = 0; c25_i60 < 4; c25_i60++) {
    c25_b_inData[c25_i60] = (*(real_T (*)[4])c25_inData)[c25_i60];
  }

  for (c25_i61 = 0; c25_i61 < 4; c25_i61++) {
    c25_u[c25_i61] = c25_b_inData[c25_i61];
  }

  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", c25_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c25_mxArrayOutData, c25_y, false);
  return c25_mxArrayOutData;
}

static const mxArray *c25_c_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData)
{
  const mxArray *c25_mxArrayOutData = NULL;
  int32_T c25_i62;
  real_T c25_b_inData[3];
  int32_T c25_i63;
  real_T c25_u[3];
  const mxArray *c25_y = NULL;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_mxArrayOutData = NULL;
  for (c25_i62 = 0; c25_i62 < 3; c25_i62++) {
    c25_b_inData[c25_i62] = (*(real_T (*)[3])c25_inData)[c25_i62];
  }

  for (c25_i63 = 0; c25_i63 < 3; c25_i63++) {
    c25_u[c25_i63] = c25_b_inData[c25_i63];
  }

  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", c25_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c25_mxArrayOutData, c25_y, false);
  return c25_mxArrayOutData;
}

static const mxArray *c25_d_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData)
{
  const mxArray *c25_mxArrayOutData = NULL;
  int32_T c25_i64;
  real_T c25_b_inData[7];
  int32_T c25_i65;
  real_T c25_u[7];
  const mxArray *c25_y = NULL;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_mxArrayOutData = NULL;
  for (c25_i64 = 0; c25_i64 < 7; c25_i64++) {
    c25_b_inData[c25_i64] = (*(real_T (*)[7])c25_inData)[c25_i64];
  }

  for (c25_i65 = 0; c25_i65 < 7; c25_i65++) {
    c25_u[c25_i65] = c25_b_inData[c25_i65];
  }

  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", c25_u, 0, 0U, 1U, 0U, 2, 7, 1), false);
  sf_mex_assign(&c25_mxArrayOutData, c25_y, false);
  return c25_mxArrayOutData;
}

static const mxArray *c25_e_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData)
{
  const mxArray *c25_mxArrayOutData = NULL;
  real_T c25_u;
  const mxArray *c25_y = NULL;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_mxArrayOutData = NULL;
  c25_u = *(real_T *)c25_inData;
  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", &c25_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c25_mxArrayOutData, c25_y, false);
  return c25_mxArrayOutData;
}

static real_T c25_c_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId)
{
  real_T c25_y;
  real_T c25_d0;
  (void)chartInstance;
  sf_mex_import(c25_parentId, sf_mex_dup(c25_u), &c25_d0, 1, 0, 0U, 0, 0U, 0);
  c25_y = c25_d0;
  sf_mex_destroy(&c25_u);
  return c25_y;
}

static void c25_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData)
{
  const mxArray *c25_nargout;
  const char_T *c25_identifier;
  emlrtMsgIdentifier c25_thisId;
  real_T c25_y;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_nargout = sf_mex_dup(c25_mxArrayInData);
  c25_identifier = c25_varName;
  c25_thisId.fIdentifier = c25_identifier;
  c25_thisId.fParent = NULL;
  c25_y = c25_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c25_nargout),
    &c25_thisId);
  sf_mex_destroy(&c25_nargout);
  *(real_T *)c25_outData = c25_y;
  sf_mex_destroy(&c25_mxArrayInData);
}

static void c25_d_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[4])
{
  real_T c25_dv2[4];
  int32_T c25_i66;
  (void)chartInstance;
  sf_mex_import(c25_parentId, sf_mex_dup(c25_u), c25_dv2, 1, 0, 0U, 1, 0U, 1, 4);
  for (c25_i66 = 0; c25_i66 < 4; c25_i66++) {
    c25_y[c25_i66] = c25_dv2[c25_i66];
  }

  sf_mex_destroy(&c25_u);
}

static void c25_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData)
{
  const mxArray *c25_v;
  const char_T *c25_identifier;
  emlrtMsgIdentifier c25_thisId;
  real_T c25_y[4];
  int32_T c25_i67;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_v = sf_mex_dup(c25_mxArrayInData);
  c25_identifier = c25_varName;
  c25_thisId.fIdentifier = c25_identifier;
  c25_thisId.fParent = NULL;
  c25_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c25_v), &c25_thisId, c25_y);
  sf_mex_destroy(&c25_v);
  for (c25_i67 = 0; c25_i67 < 4; c25_i67++) {
    (*(real_T (*)[4])c25_outData)[c25_i67] = c25_y[c25_i67];
  }

  sf_mex_destroy(&c25_mxArrayInData);
}

static void c25_e_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[3])
{
  real_T c25_dv3[3];
  int32_T c25_i68;
  (void)chartInstance;
  sf_mex_import(c25_parentId, sf_mex_dup(c25_u), c25_dv3, 1, 0, 0U, 1, 0U, 1, 3);
  for (c25_i68 = 0; c25_i68 < 3; c25_i68++) {
    c25_y[c25_i68] = c25_dv3[c25_i68];
  }

  sf_mex_destroy(&c25_u);
}

static void c25_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData)
{
  const mxArray *c25_rc;
  const char_T *c25_identifier;
  emlrtMsgIdentifier c25_thisId;
  real_T c25_y[3];
  int32_T c25_i69;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_rc = sf_mex_dup(c25_mxArrayInData);
  c25_identifier = c25_varName;
  c25_thisId.fIdentifier = c25_identifier;
  c25_thisId.fParent = NULL;
  c25_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c25_rc), &c25_thisId, c25_y);
  sf_mex_destroy(&c25_rc);
  for (c25_i69 = 0; c25_i69 < 3; c25_i69++) {
    (*(real_T (*)[3])c25_outData)[c25_i69] = c25_y[c25_i69];
  }

  sf_mex_destroy(&c25_mxArrayInData);
}

static const mxArray *c25_f_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData)
{
  const mxArray *c25_mxArrayOutData = NULL;
  int32_T c25_i70;
  int32_T c25_i71;
  int32_T c25_i72;
  real_T c25_b_inData[9];
  int32_T c25_i73;
  int32_T c25_i74;
  int32_T c25_i75;
  real_T c25_u[9];
  const mxArray *c25_y = NULL;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_mxArrayOutData = NULL;
  c25_i70 = 0;
  for (c25_i71 = 0; c25_i71 < 3; c25_i71++) {
    for (c25_i72 = 0; c25_i72 < 3; c25_i72++) {
      c25_b_inData[c25_i72 + c25_i70] = (*(real_T (*)[9])c25_inData)[c25_i72 +
        c25_i70];
    }

    c25_i70 += 3;
  }

  c25_i73 = 0;
  for (c25_i74 = 0; c25_i74 < 3; c25_i74++) {
    for (c25_i75 = 0; c25_i75 < 3; c25_i75++) {
      c25_u[c25_i75 + c25_i73] = c25_b_inData[c25_i75 + c25_i73];
    }

    c25_i73 += 3;
  }

  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", c25_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c25_mxArrayOutData, c25_y, false);
  return c25_mxArrayOutData;
}

static void c25_f_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[9])
{
  real_T c25_dv4[9];
  int32_T c25_i76;
  (void)chartInstance;
  sf_mex_import(c25_parentId, sf_mex_dup(c25_u), c25_dv4, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c25_i76 = 0; c25_i76 < 9; c25_i76++) {
    c25_y[c25_i76] = c25_dv4[c25_i76];
  }

  sf_mex_destroy(&c25_u);
}

static void c25_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData)
{
  const mxArray *c25_SkewSymmetricTensor;
  const char_T *c25_identifier;
  emlrtMsgIdentifier c25_thisId;
  real_T c25_y[9];
  int32_T c25_i77;
  int32_T c25_i78;
  int32_T c25_i79;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_SkewSymmetricTensor = sf_mex_dup(c25_mxArrayInData);
  c25_identifier = c25_varName;
  c25_thisId.fIdentifier = c25_identifier;
  c25_thisId.fParent = NULL;
  c25_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c25_SkewSymmetricTensor),
    &c25_thisId, c25_y);
  sf_mex_destroy(&c25_SkewSymmetricTensor);
  c25_i77 = 0;
  for (c25_i78 = 0; c25_i78 < 3; c25_i78++) {
    for (c25_i79 = 0; c25_i79 < 3; c25_i79++) {
      (*(real_T (*)[9])c25_outData)[c25_i79 + c25_i77] = c25_y[c25_i79 + c25_i77];
    }

    c25_i77 += 3;
  }

  sf_mex_destroy(&c25_mxArrayInData);
}

static const mxArray *c25_g_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData)
{
  const mxArray *c25_mxArrayOutData = NULL;
  int32_T c25_i80;
  int32_T c25_i81;
  int32_T c25_i82;
  real_T c25_b_inData[16];
  int32_T c25_i83;
  int32_T c25_i84;
  int32_T c25_i85;
  real_T c25_u[16];
  const mxArray *c25_y = NULL;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_mxArrayOutData = NULL;
  c25_i80 = 0;
  for (c25_i81 = 0; c25_i81 < 4; c25_i81++) {
    for (c25_i82 = 0; c25_i82 < 4; c25_i82++) {
      c25_b_inData[c25_i82 + c25_i80] = (*(real_T (*)[16])c25_inData)[c25_i82 +
        c25_i80];
    }

    c25_i80 += 4;
  }

  c25_i83 = 0;
  for (c25_i84 = 0; c25_i84 < 4; c25_i84++) {
    for (c25_i85 = 0; c25_i85 < 4; c25_i85++) {
      c25_u[c25_i85 + c25_i83] = c25_b_inData[c25_i85 + c25_i83];
    }

    c25_i83 += 4;
  }

  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", c25_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c25_mxArrayOutData, c25_y, false);
  return c25_mxArrayOutData;
}

static void c25_g_emlrt_marshallIn(SFc25_Model_01InstanceStruct *chartInstance,
  const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId, real_T c25_y[16])
{
  real_T c25_dv5[16];
  int32_T c25_i86;
  (void)chartInstance;
  sf_mex_import(c25_parentId, sf_mex_dup(c25_u), c25_dv5, 1, 0, 0U, 1, 0U, 2, 4,
                4);
  for (c25_i86 = 0; c25_i86 < 16; c25_i86++) {
    c25_y[c25_i86] = c25_dv5[c25_i86];
  }

  sf_mex_destroy(&c25_u);
}

static void c25_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData)
{
  const mxArray *c25_CrossTensor;
  const char_T *c25_identifier;
  emlrtMsgIdentifier c25_thisId;
  real_T c25_y[16];
  int32_T c25_i87;
  int32_T c25_i88;
  int32_T c25_i89;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_CrossTensor = sf_mex_dup(c25_mxArrayInData);
  c25_identifier = c25_varName;
  c25_thisId.fIdentifier = c25_identifier;
  c25_thisId.fParent = NULL;
  c25_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c25_CrossTensor), &c25_thisId,
    c25_y);
  sf_mex_destroy(&c25_CrossTensor);
  c25_i87 = 0;
  for (c25_i88 = 0; c25_i88 < 4; c25_i88++) {
    for (c25_i89 = 0; c25_i89 < 4; c25_i89++) {
      (*(real_T (*)[16])c25_outData)[c25_i89 + c25_i87] = c25_y[c25_i89 +
        c25_i87];
    }

    c25_i87 += 4;
  }

  sf_mex_destroy(&c25_mxArrayInData);
}

const mxArray *sf_c25_Model_01_get_eml_resolved_functions_info(void)
{
  const mxArray *c25_nameCaptureInfo = NULL;
  c25_nameCaptureInfo = NULL;
  sf_mex_assign(&c25_nameCaptureInfo, sf_mex_createstruct("structure", 2, 34, 1),
                false);
  c25_info_helper(&c25_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c25_nameCaptureInfo);
  return c25_nameCaptureInfo;
}

static void c25_info_helper(const mxArray **c25_info)
{
  const mxArray *c25_rhs0 = NULL;
  const mxArray *c25_lhs0 = NULL;
  const mxArray *c25_rhs1 = NULL;
  const mxArray *c25_lhs1 = NULL;
  const mxArray *c25_rhs2 = NULL;
  const mxArray *c25_lhs2 = NULL;
  const mxArray *c25_rhs3 = NULL;
  const mxArray *c25_lhs3 = NULL;
  const mxArray *c25_rhs4 = NULL;
  const mxArray *c25_lhs4 = NULL;
  const mxArray *c25_rhs5 = NULL;
  const mxArray *c25_lhs5 = NULL;
  const mxArray *c25_rhs6 = NULL;
  const mxArray *c25_lhs6 = NULL;
  const mxArray *c25_rhs7 = NULL;
  const mxArray *c25_lhs7 = NULL;
  const mxArray *c25_rhs8 = NULL;
  const mxArray *c25_lhs8 = NULL;
  const mxArray *c25_rhs9 = NULL;
  const mxArray *c25_lhs9 = NULL;
  const mxArray *c25_rhs10 = NULL;
  const mxArray *c25_lhs10 = NULL;
  const mxArray *c25_rhs11 = NULL;
  const mxArray *c25_lhs11 = NULL;
  const mxArray *c25_rhs12 = NULL;
  const mxArray *c25_lhs12 = NULL;
  const mxArray *c25_rhs13 = NULL;
  const mxArray *c25_lhs13 = NULL;
  const mxArray *c25_rhs14 = NULL;
  const mxArray *c25_lhs14 = NULL;
  const mxArray *c25_rhs15 = NULL;
  const mxArray *c25_lhs15 = NULL;
  const mxArray *c25_rhs16 = NULL;
  const mxArray *c25_lhs16 = NULL;
  const mxArray *c25_rhs17 = NULL;
  const mxArray *c25_lhs17 = NULL;
  const mxArray *c25_rhs18 = NULL;
  const mxArray *c25_lhs18 = NULL;
  const mxArray *c25_rhs19 = NULL;
  const mxArray *c25_lhs19 = NULL;
  const mxArray *c25_rhs20 = NULL;
  const mxArray *c25_lhs20 = NULL;
  const mxArray *c25_rhs21 = NULL;
  const mxArray *c25_lhs21 = NULL;
  const mxArray *c25_rhs22 = NULL;
  const mxArray *c25_lhs22 = NULL;
  const mxArray *c25_rhs23 = NULL;
  const mxArray *c25_lhs23 = NULL;
  const mxArray *c25_rhs24 = NULL;
  const mxArray *c25_lhs24 = NULL;
  const mxArray *c25_rhs25 = NULL;
  const mxArray *c25_lhs25 = NULL;
  const mxArray *c25_rhs26 = NULL;
  const mxArray *c25_lhs26 = NULL;
  const mxArray *c25_rhs27 = NULL;
  const mxArray *c25_lhs27 = NULL;
  const mxArray *c25_rhs28 = NULL;
  const mxArray *c25_lhs28 = NULL;
  const mxArray *c25_rhs29 = NULL;
  const mxArray *c25_lhs29 = NULL;
  const mxArray *c25_rhs30 = NULL;
  const mxArray *c25_lhs30 = NULL;
  const mxArray *c25_rhs31 = NULL;
  const mxArray *c25_lhs31 = NULL;
  const mxArray *c25_rhs32 = NULL;
  const mxArray *c25_lhs32 = NULL;
  const mxArray *c25_rhs33 = NULL;
  const mxArray *c25_lhs33 = NULL;
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("fn_CrossTensor"), "name",
                  "name", 0);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1450227411U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c25_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "context", "context", 1);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 1);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1450040424U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c25_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "context", "context", 2);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eye"), "name", "name", 2);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c25_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 3);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c25_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 4);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c25_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 5);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("isinf"), "name", "name", 5);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c25_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 6);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c25_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 7);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_is_integer_class"),
                  "name", "name", 7);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c25_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 8);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("intmax"), "name", "name", 8);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c25_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 9);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c25_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 10);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("intmin"), "name", "name", 10);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c25_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 11);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c25_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 12);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "coder.internal.indexIntRelop"), "name", "name", 12);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c25_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 13);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 13);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c25_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 14);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 14);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c25_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 15);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("intmin"), "name", "name", 15);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c25_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 16);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c25_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("intmax"), "name", "name", 17);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c25_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 18);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c25_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 19);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("intmax"), "name", "name", 19);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c25_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_CrossTensor.m"),
                  "context", "context", 20);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 20);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c25_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 21);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 21);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c25_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "context", "context", 22);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 22);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c25_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 23);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 23);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c25_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 24);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 24);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c25_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 25);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c25_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 26);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  26);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c25_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 27);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c25_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 28);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c25_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 29);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 29);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c25_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 30);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 30);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c25_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 31);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 31);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c25_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 32);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 32);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c25_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 33);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "coder.internal.refblas.xgemm"), "name", "name", 33);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c25_info, c25_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c25_info, c25_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c25_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c25_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c25_info, sf_mex_duplicatearraysafe(&c25_lhs33), "lhs", "lhs",
                  33);
  sf_mex_destroy(&c25_rhs0);
  sf_mex_destroy(&c25_lhs0);
  sf_mex_destroy(&c25_rhs1);
  sf_mex_destroy(&c25_lhs1);
  sf_mex_destroy(&c25_rhs2);
  sf_mex_destroy(&c25_lhs2);
  sf_mex_destroy(&c25_rhs3);
  sf_mex_destroy(&c25_lhs3);
  sf_mex_destroy(&c25_rhs4);
  sf_mex_destroy(&c25_lhs4);
  sf_mex_destroy(&c25_rhs5);
  sf_mex_destroy(&c25_lhs5);
  sf_mex_destroy(&c25_rhs6);
  sf_mex_destroy(&c25_lhs6);
  sf_mex_destroy(&c25_rhs7);
  sf_mex_destroy(&c25_lhs7);
  sf_mex_destroy(&c25_rhs8);
  sf_mex_destroy(&c25_lhs8);
  sf_mex_destroy(&c25_rhs9);
  sf_mex_destroy(&c25_lhs9);
  sf_mex_destroy(&c25_rhs10);
  sf_mex_destroy(&c25_lhs10);
  sf_mex_destroy(&c25_rhs11);
  sf_mex_destroy(&c25_lhs11);
  sf_mex_destroy(&c25_rhs12);
  sf_mex_destroy(&c25_lhs12);
  sf_mex_destroy(&c25_rhs13);
  sf_mex_destroy(&c25_lhs13);
  sf_mex_destroy(&c25_rhs14);
  sf_mex_destroy(&c25_lhs14);
  sf_mex_destroy(&c25_rhs15);
  sf_mex_destroy(&c25_lhs15);
  sf_mex_destroy(&c25_rhs16);
  sf_mex_destroy(&c25_lhs16);
  sf_mex_destroy(&c25_rhs17);
  sf_mex_destroy(&c25_lhs17);
  sf_mex_destroy(&c25_rhs18);
  sf_mex_destroy(&c25_lhs18);
  sf_mex_destroy(&c25_rhs19);
  sf_mex_destroy(&c25_lhs19);
  sf_mex_destroy(&c25_rhs20);
  sf_mex_destroy(&c25_lhs20);
  sf_mex_destroy(&c25_rhs21);
  sf_mex_destroy(&c25_lhs21);
  sf_mex_destroy(&c25_rhs22);
  sf_mex_destroy(&c25_lhs22);
  sf_mex_destroy(&c25_rhs23);
  sf_mex_destroy(&c25_lhs23);
  sf_mex_destroy(&c25_rhs24);
  sf_mex_destroy(&c25_lhs24);
  sf_mex_destroy(&c25_rhs25);
  sf_mex_destroy(&c25_lhs25);
  sf_mex_destroy(&c25_rhs26);
  sf_mex_destroy(&c25_lhs26);
  sf_mex_destroy(&c25_rhs27);
  sf_mex_destroy(&c25_lhs27);
  sf_mex_destroy(&c25_rhs28);
  sf_mex_destroy(&c25_lhs28);
  sf_mex_destroy(&c25_rhs29);
  sf_mex_destroy(&c25_lhs29);
  sf_mex_destroy(&c25_rhs30);
  sf_mex_destroy(&c25_lhs30);
  sf_mex_destroy(&c25_rhs31);
  sf_mex_destroy(&c25_lhs31);
  sf_mex_destroy(&c25_rhs32);
  sf_mex_destroy(&c25_lhs32);
  sf_mex_destroy(&c25_rhs33);
  sf_mex_destroy(&c25_lhs33);
}

static const mxArray *c25_emlrt_marshallOut(const char * c25_u)
{
  const mxArray *c25_y = NULL;
  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", c25_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c25_u)), false);
  return c25_y;
}

static const mxArray *c25_b_emlrt_marshallOut(const uint32_T c25_u)
{
  const mxArray *c25_y = NULL;
  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", &c25_u, 7, 0U, 0U, 0U, 0), false);
  return c25_y;
}

static void c25_eml_scalar_eg(SFc25_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c25_threshold(SFc25_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c25_b_eml_scalar_eg(SFc25_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *c25_h_sf_marshallOut(void *chartInstanceVoid, void
  *c25_inData)
{
  const mxArray *c25_mxArrayOutData = NULL;
  int32_T c25_u;
  const mxArray *c25_y = NULL;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_mxArrayOutData = NULL;
  c25_u = *(int32_T *)c25_inData;
  c25_y = NULL;
  sf_mex_assign(&c25_y, sf_mex_create("y", &c25_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c25_mxArrayOutData, c25_y, false);
  return c25_mxArrayOutData;
}

static int32_T c25_h_emlrt_marshallIn(SFc25_Model_01InstanceStruct
  *chartInstance, const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId)
{
  int32_T c25_y;
  int32_T c25_i90;
  (void)chartInstance;
  sf_mex_import(c25_parentId, sf_mex_dup(c25_u), &c25_i90, 1, 6, 0U, 0, 0U, 0);
  c25_y = c25_i90;
  sf_mex_destroy(&c25_u);
  return c25_y;
}

static void c25_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c25_mxArrayInData, const char_T *c25_varName, void *c25_outData)
{
  const mxArray *c25_b_sfEvent;
  const char_T *c25_identifier;
  emlrtMsgIdentifier c25_thisId;
  int32_T c25_y;
  SFc25_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc25_Model_01InstanceStruct *)chartInstanceVoid;
  c25_b_sfEvent = sf_mex_dup(c25_mxArrayInData);
  c25_identifier = c25_varName;
  c25_thisId.fIdentifier = c25_identifier;
  c25_thisId.fParent = NULL;
  c25_y = c25_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c25_b_sfEvent),
    &c25_thisId);
  sf_mex_destroy(&c25_b_sfEvent);
  *(int32_T *)c25_outData = c25_y;
  sf_mex_destroy(&c25_mxArrayInData);
}

static uint8_T c25_i_emlrt_marshallIn(SFc25_Model_01InstanceStruct
  *chartInstance, const mxArray *c25_b_is_active_c25_Model_01, const char_T
  *c25_identifier)
{
  uint8_T c25_y;
  emlrtMsgIdentifier c25_thisId;
  c25_thisId.fIdentifier = c25_identifier;
  c25_thisId.fParent = NULL;
  c25_y = c25_j_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c25_b_is_active_c25_Model_01), &c25_thisId);
  sf_mex_destroy(&c25_b_is_active_c25_Model_01);
  return c25_y;
}

static uint8_T c25_j_emlrt_marshallIn(SFc25_Model_01InstanceStruct
  *chartInstance, const mxArray *c25_u, const emlrtMsgIdentifier *c25_parentId)
{
  uint8_T c25_y;
  uint8_T c25_u0;
  (void)chartInstance;
  sf_mex_import(c25_parentId, sf_mex_dup(c25_u), &c25_u0, 1, 3, 0U, 0, 0U, 0);
  c25_y = c25_u0;
  sf_mex_destroy(&c25_u);
  return c25_y;
}

static void init_dsm_address_info(SFc25_Model_01InstanceStruct *chartInstance)
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

void sf_c25_Model_01_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1240831421U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2049415141U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(4166439370U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4051849543U);
}

mxArray *sf_c25_Model_01_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("Fl47OdJDSWCKRT01czohIG");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(7);
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
      pr[0] = (double)(6);
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

mxArray *sf_c25_Model_01_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c25_Model_01_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c25_Model_01(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"zk\",},{M[8],M[0],T\"is_active_c25_Model_01\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c25_Model_01_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc25_Model_01InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc25_Model_01InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_01MachineNumber_,
           25,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"signal");
          _SFD_SET_DATA_PROPS(1,2,0,1,"zk");
          _SFD_SET_DATA_PROPS(2,1,1,0,"rho_c");
          _SFD_SET_DATA_PROPS(3,1,1,0,"q_k");
          _SFD_SET_DATA_PROPS(4,1,1,0,"ita_k");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,322);
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
          unsigned int dimVector[2];
          dimVector[0]= 7;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c25_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c25_sf_marshallOut,(MexInFcnForType)
            c25_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c25_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c25_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c25_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T (*c25_signal)[7];
          real_T (*c25_zk)[6];
          real_T (*c25_rho_c)[3];
          real_T (*c25_q_k)[4];
          real_T (*c25_ita_k)[4];
          c25_ita_k = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 3);
          c25_q_k = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 2);
          c25_rho_c = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
          c25_zk = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
          c25_signal = (real_T (*)[7])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c25_signal);
          _SFD_SET_DATA_VALUE_PTR(1U, *c25_zk);
          _SFD_SET_DATA_VALUE_PTR(2U, *c25_rho_c);
          _SFD_SET_DATA_VALUE_PTR(3U, *c25_q_k);
          _SFD_SET_DATA_VALUE_PTR(4U, *c25_ita_k);
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
  return "dwoFzTByORxjD9xyp4fZcE";
}

static void sf_opaque_initialize_c25_Model_01(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc25_Model_01InstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c25_Model_01((SFc25_Model_01InstanceStruct*)
    chartInstanceVar);
  initialize_c25_Model_01((SFc25_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c25_Model_01(void *chartInstanceVar)
{
  enable_c25_Model_01((SFc25_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c25_Model_01(void *chartInstanceVar)
{
  disable_c25_Model_01((SFc25_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c25_Model_01(void *chartInstanceVar)
{
  sf_gateway_c25_Model_01((SFc25_Model_01InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c25_Model_01(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c25_Model_01((SFc25_Model_01InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c25_Model_01();/* state var info */
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

extern void sf_internal_set_sim_state_c25_Model_01(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c25_Model_01();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c25_Model_01((SFc25_Model_01InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c25_Model_01(SimStruct* S)
{
  return sf_internal_get_sim_state_c25_Model_01(S);
}

static void sf_opaque_set_sim_state_c25_Model_01(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c25_Model_01(S, st);
}

static void sf_opaque_terminate_c25_Model_01(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc25_Model_01InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_01_optimization_info();
    }

    finalize_c25_Model_01((SFc25_Model_01InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc25_Model_01((SFc25_Model_01InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c25_Model_01(SimStruct *S)
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
    initialize_params_c25_Model_01((SFc25_Model_01InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c25_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_01_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      25);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,25,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,25,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,25);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,25,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,25,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,25);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2974963941U));
  ssSetChecksum1(S,(3219236262U));
  ssSetChecksum2(S,(955568224U));
  ssSetChecksum3(S,(2268901172U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c25_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c25_Model_01(SimStruct *S)
{
  SFc25_Model_01InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc25_Model_01InstanceStruct *)utMalloc(sizeof
    (SFc25_Model_01InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc25_Model_01InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c25_Model_01;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c25_Model_01;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c25_Model_01;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c25_Model_01;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c25_Model_01;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c25_Model_01;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c25_Model_01;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c25_Model_01;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c25_Model_01;
  chartInstance->chartInfo.mdlStart = mdlStart_c25_Model_01;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c25_Model_01;
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

void c25_Model_01_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c25_Model_01(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c25_Model_01(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c25_Model_01(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c25_Model_01_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

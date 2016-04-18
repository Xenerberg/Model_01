/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_02_sfun.h"
#include "c7_Model_02.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "Model_02_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c7_debug_family_names[13] = { "J", "Q_r11", "Q_r12", "Q_r22",
  "nargin", "nargout", "Phi", "unused1", "t_delta", "tau", "unused", "p", "Q_r"
};

static const char * c7_b_debug_family_names[4] = { "nargin", "nargout", "p", "J"
};

static const char * c7_c_debug_family_names[7] = { "nargin", "nargout",
  "phi_r12", "J", "t_delta", "tau", "Q_r11" };

static const char * c7_d_debug_family_names[8] = { "nargin", "nargout",
  "phi_r12", "phi_r22", "J", "t_delta", "tau", "Q_r12" };

static const char * c7_e_debug_family_names[7] = { "nargin", "nargout",
  "phi_r22", "J", "t_delta", "tau", "Q_r22" };

/* Function Declarations */
static void initialize_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance);
static void initialize_params_c7_Model_02(SFc7_Model_02InstanceStruct
  *chartInstance);
static void enable_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance);
static void disable_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance);
static void c7_update_debugger_state_c7_Model_02(SFc7_Model_02InstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c7_Model_02(SFc7_Model_02InstanceStruct
  *chartInstance);
static void set_sim_state_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_st);
static void finalize_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance);
static void sf_gateway_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance);
static void c7_chartstep_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance);
static void initSimStructsc7_Model_02(SFc7_Model_02InstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c7_machineNumber, uint32_T
  c7_chartNumber, uint32_T c7_instanceNumber);
static const mxArray *c7_sf_marshallOut(void *chartInstanceVoid, void *c7_inData);
static void c7_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_Q_r, const char_T *c7_identifier, real_T c7_y[36]);
static void c7_b_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId, real_T c7_y[36]);
static void c7_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData);
static const mxArray *c7_b_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData);
static const mxArray *c7_c_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData);
static const mxArray *c7_d_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData);
static real_T c7_c_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId);
static void c7_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData);
static const mxArray *c7_e_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData);
static void c7_d_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId, real_T c7_y[9]);
static void c7_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData);
static void c7_e_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId, real_T c7_y[3]);
static void c7_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData);
static void c7_info_helper(const mxArray **c7_info);
static const mxArray *c7_emlrt_marshallOut(const char * c7_u);
static const mxArray *c7_b_emlrt_marshallOut(const uint32_T c7_u);
static void c7_fn_Create_Q_r11(SFc7_Model_02InstanceStruct *chartInstance,
  real_T c7_phi_r12[9], real_T c7_J[9], real_T c7_t_delta, real_T c7_tau, real_T
  c7_Q_r11[9]);
static void c7_mpower(SFc7_Model_02InstanceStruct *chartInstance, real_T c7_a[9],
                      real_T c7_c[9]);
static void c7_eml_scalar_eg(SFc7_Model_02InstanceStruct *chartInstance);
static void c7_eml_xgemm(SFc7_Model_02InstanceStruct *chartInstance, real_T
  c7_A[9], real_T c7_B[9], real_T c7_C[9], real_T c7_b_C[9]);
static const mxArray *c7_f_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData);
static int32_T c7_f_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId);
static void c7_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData);
static uint8_T c7_g_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_b_is_active_c7_Model_02, const char_T *c7_identifier);
static uint8_T c7_h_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId);
static void c7_b_eml_xgemm(SFc7_Model_02InstanceStruct *chartInstance, real_T
  c7_A[9], real_T c7_B[9], real_T c7_C[9]);
static void init_dsm_address_info(SFc7_Model_02InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance)
{
  chartInstance->c7_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c7_is_active_c7_Model_02 = 0U;
}

static void initialize_params_c7_Model_02(SFc7_Model_02InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c7_update_debugger_state_c7_Model_02(SFc7_Model_02InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c7_Model_02(SFc7_Model_02InstanceStruct
  *chartInstance)
{
  const mxArray *c7_st;
  const mxArray *c7_y = NULL;
  int32_T c7_i0;
  real_T c7_u[36];
  const mxArray *c7_b_y = NULL;
  uint8_T c7_hoistedGlobal;
  uint8_T c7_b_u;
  const mxArray *c7_c_y = NULL;
  real_T (*c7_Q_r)[36];
  c7_Q_r = (real_T (*)[36])ssGetOutputPortSignal(chartInstance->S, 1);
  c7_st = NULL;
  c7_st = NULL;
  c7_y = NULL;
  sf_mex_assign(&c7_y, sf_mex_createcellmatrix(2, 1), false);
  for (c7_i0 = 0; c7_i0 < 36; c7_i0++) {
    c7_u[c7_i0] = (*c7_Q_r)[c7_i0];
  }

  c7_b_y = NULL;
  sf_mex_assign(&c7_b_y, sf_mex_create("y", c7_u, 0, 0U, 1U, 0U, 2, 6, 6), false);
  sf_mex_setcell(c7_y, 0, c7_b_y);
  c7_hoistedGlobal = chartInstance->c7_is_active_c7_Model_02;
  c7_b_u = c7_hoistedGlobal;
  c7_c_y = NULL;
  sf_mex_assign(&c7_c_y, sf_mex_create("y", &c7_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c7_y, 1, c7_c_y);
  sf_mex_assign(&c7_st, c7_y, false);
  return c7_st;
}

static void set_sim_state_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_st)
{
  const mxArray *c7_u;
  real_T c7_dv0[36];
  int32_T c7_i1;
  real_T (*c7_Q_r)[36];
  c7_Q_r = (real_T (*)[36])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c7_doneDoubleBufferReInit = true;
  c7_u = sf_mex_dup(c7_st);
  c7_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c7_u, 0)), "Q_r",
                      c7_dv0);
  for (c7_i1 = 0; c7_i1 < 36; c7_i1++) {
    (*c7_Q_r)[c7_i1] = c7_dv0[c7_i1];
  }

  chartInstance->c7_is_active_c7_Model_02 = c7_g_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c7_u, 1)), "is_active_c7_Model_02");
  sf_mex_destroy(&c7_u);
  c7_update_debugger_state_c7_Model_02(chartInstance);
  sf_mex_destroy(&c7_st);
}

static void finalize_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance)
{
  int32_T c7_i2;
  int32_T c7_i3;
  int32_T c7_i4;
  real_T *c7_unused1;
  real_T *c7_t_delta;
  real_T *c7_tau;
  real_T *c7_unused;
  real_T (*c7_p)[3];
  real_T (*c7_Q_r)[36];
  real_T (*c7_Phi)[144];
  c7_p = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c7_unused = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c7_tau = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c7_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c7_unused1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c7_Q_r = (real_T (*)[36])ssGetOutputPortSignal(chartInstance->S, 1);
  c7_Phi = (real_T (*)[144])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 6U, chartInstance->c7_sfEvent);
  for (c7_i2 = 0; c7_i2 < 144; c7_i2++) {
    _SFD_DATA_RANGE_CHECK((*c7_Phi)[c7_i2], 0U);
  }

  chartInstance->c7_sfEvent = CALL_EVENT;
  c7_chartstep_c7_Model_02(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_02MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c7_i3 = 0; c7_i3 < 36; c7_i3++) {
    _SFD_DATA_RANGE_CHECK((*c7_Q_r)[c7_i3], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c7_unused1, 2U);
  _SFD_DATA_RANGE_CHECK(*c7_t_delta, 3U);
  _SFD_DATA_RANGE_CHECK(*c7_tau, 4U);
  _SFD_DATA_RANGE_CHECK(*c7_unused, 5U);
  for (c7_i4 = 0; c7_i4 < 3; c7_i4++) {
    _SFD_DATA_RANGE_CHECK((*c7_p)[c7_i4], 6U);
  }
}

static void c7_chartstep_c7_Model_02(SFc7_Model_02InstanceStruct *chartInstance)
{
  real_T c7_hoistedGlobal;
  real_T c7_b_hoistedGlobal;
  real_T c7_c_hoistedGlobal;
  real_T c7_d_hoistedGlobal;
  int32_T c7_i5;
  real_T c7_Phi[144];
  real_T c7_unused1;
  real_T c7_t_delta;
  real_T c7_tau;
  real_T c7_unused;
  int32_T c7_i6;
  real_T c7_p[3];
  uint32_T c7_debug_family_var_map[13];
  real_T c7_J[9];
  real_T c7_Q_r11[9];
  real_T c7_Q_r12[9];
  real_T c7_Q_r22[9];
  real_T c7_nargin = 6.0;
  real_T c7_nargout = 1.0;
  real_T c7_Q_r[36];
  int32_T c7_i7;
  real_T c7_b_p[3];
  uint32_T c7_b_debug_family_var_map[4];
  real_T c7_b_nargin = 1.0;
  real_T c7_b_nargout = 1.0;
  real_T c7_A;
  real_T c7_B;
  real_T c7_x;
  real_T c7_y;
  real_T c7_b_x;
  real_T c7_b_y;
  real_T c7_c_x;
  real_T c7_c_y;
  real_T c7_d_y;
  real_T c7_b_A;
  real_T c7_b_B;
  real_T c7_d_x;
  real_T c7_e_y;
  real_T c7_e_x;
  real_T c7_f_y;
  real_T c7_f_x;
  real_T c7_g_y;
  real_T c7_h_y;
  int32_T c7_i8;
  int32_T c7_i9;
  static real_T c7_dv1[3] = { 1.0, 0.0, 0.0 };

  int32_T c7_i10;
  int32_T c7_i11;
  int32_T c7_i12;
  int32_T c7_i13;
  real_T c7_b_Phi[9];
  int32_T c7_i14;
  real_T c7_b_J[9];
  real_T c7_b[9];
  int32_T c7_i15;
  int32_T c7_i16;
  int32_T c7_i17;
  int32_T c7_i18;
  int32_T c7_i19;
  real_T c7_phi_r12[9];
  int32_T c7_i20;
  int32_T c7_i21;
  int32_T c7_i22;
  int32_T c7_i23;
  real_T c7_phi_r22[9];
  int32_T c7_i24;
  real_T c7_c_J[9];
  real_T c7_b_t_delta;
  real_T c7_b_tau;
  uint32_T c7_c_debug_family_var_map[8];
  real_T c7_c_nargin = 5.0;
  real_T c7_c_nargout = 1.0;
  real_T c7_b_Q_r12[9];
  real_T c7_a;
  int32_T c7_i25;
  int32_T c7_i26;
  int32_T c7_i27;
  real_T c7_d_J[9];
  real_T c7_b_b[9];
  int32_T c7_i28;
  real_T c7_i_y[9];
  int32_T c7_i29;
  real_T c7_c_b[9];
  int32_T c7_i30;
  real_T c7_d_b[9];
  int32_T c7_i31;
  int32_T c7_i32;
  int32_T c7_i33;
  real_T c7_j_y[9];
  int32_T c7_i34;
  real_T c7_e_b[9];
  real_T c7_f_b;
  int32_T c7_i35;
  int32_T c7_i36;
  int32_T c7_i37;
  int32_T c7_i38;
  int32_T c7_i39;
  int32_T c7_i40;
  int32_T c7_i41;
  real_T c7_b_phi_r22[9];
  int32_T c7_i42;
  real_T c7_e_J[9];
  real_T c7_c_t_delta;
  real_T c7_c_tau;
  uint32_T c7_d_debug_family_var_map[7];
  real_T c7_d_nargin = 4.0;
  real_T c7_d_nargout = 1.0;
  real_T c7_b_a;
  int32_T c7_i43;
  int32_T c7_i44;
  int32_T c7_i45;
  real_T c7_f_J[9];
  int32_T c7_i46;
  int32_T c7_i47;
  real_T c7_g_b[9];
  int32_T c7_i48;
  real_T c7_h_b[9];
  int32_T c7_i49;
  int32_T c7_i50;
  int32_T c7_i51;
  real_T c7_k_y[9];
  int32_T c7_i52;
  real_T c7_i_b[9];
  real_T c7_j_b;
  int32_T c7_i53;
  int32_T c7_i54;
  int32_T c7_i55;
  int32_T c7_i56;
  int32_T c7_i57;
  int32_T c7_i58;
  int32_T c7_i59;
  int32_T c7_i60;
  int32_T c7_i61;
  int32_T c7_i62;
  int32_T c7_i63;
  int32_T c7_i64;
  int32_T c7_i65;
  int32_T c7_i66;
  int32_T c7_i67;
  int32_T c7_i68;
  int32_T c7_i69;
  int32_T c7_i70;
  real_T *c7_b_unused;
  real_T *c7_d_tau;
  real_T *c7_d_t_delta;
  real_T *c7_b_unused1;
  real_T (*c7_b_Q_r)[36];
  real_T (*c7_c_p)[3];
  real_T (*c7_c_Phi)[144];
  c7_c_p = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c7_b_unused = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c7_d_tau = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c7_d_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c7_b_unused1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c7_b_Q_r = (real_T (*)[36])ssGetOutputPortSignal(chartInstance->S, 1);
  c7_c_Phi = (real_T (*)[144])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 6U, chartInstance->c7_sfEvent);
  c7_hoistedGlobal = *c7_b_unused1;
  c7_b_hoistedGlobal = *c7_d_t_delta;
  c7_c_hoistedGlobal = *c7_d_tau;
  c7_d_hoistedGlobal = *c7_b_unused;
  for (c7_i5 = 0; c7_i5 < 144; c7_i5++) {
    c7_Phi[c7_i5] = (*c7_c_Phi)[c7_i5];
  }

  c7_unused1 = c7_hoistedGlobal;
  c7_t_delta = c7_b_hoistedGlobal;
  c7_tau = c7_c_hoistedGlobal;
  c7_unused = c7_d_hoistedGlobal;
  for (c7_i6 = 0; c7_i6 < 3; c7_i6++) {
    c7_p[c7_i6] = (*c7_c_p)[c7_i6];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 13U, 13U, c7_debug_family_names,
    c7_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_J, 0U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_Q_r11, 1U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_Q_r12, 2U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_Q_r22, 3U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_nargin, 4U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_nargout, 5U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c7_Phi, 6U, c7_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c7_unused1, 7U, c7_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c7_t_delta, 8U, c7_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c7_tau, 9U, c7_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c7_unused, 10U, c7_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c7_p, 11U, c7_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_Q_r, 12U, c7_sf_marshallOut,
    c7_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, 12);
  for (c7_i7 = 0; c7_i7 < 3; c7_i7++) {
    c7_b_p[c7_i7] = c7_p[c7_i7];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c7_b_debug_family_names,
    c7_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_b_nargin, 0U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_b_nargout, 1U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_b_p, 2U, c7_b_sf_marshallOut,
    c7_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_J, 3U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  CV_EML_FCN(0, 4);
  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, 36);
  c7_A = 1.0 - c7_b_p[1];
  c7_B = 1.0 + c7_b_p[0];
  c7_x = c7_A;
  c7_y = c7_B;
  c7_b_x = c7_x;
  c7_b_y = c7_y;
  c7_c_x = c7_b_x;
  c7_c_y = c7_b_y;
  c7_d_y = c7_c_x / c7_c_y;
  c7_b_A = 1.0 + c7_b_p[2];
  c7_b_B = 1.0 - c7_b_p[0];
  c7_d_x = c7_b_A;
  c7_e_y = c7_b_B;
  c7_e_x = c7_d_x;
  c7_f_y = c7_e_y;
  c7_f_x = c7_e_x;
  c7_g_y = c7_f_y;
  c7_h_y = c7_f_x / c7_g_y;
  c7_i8 = 0;
  for (c7_i9 = 0; c7_i9 < 3; c7_i9++) {
    c7_J[c7_i8] = c7_dv1[c7_i9];
    c7_i8 += 3;
  }

  c7_J[1] = 0.0;
  c7_J[4] = c7_d_y;
  c7_J[7] = 0.0;
  c7_J[2] = 0.0;
  c7_J[5] = 0.0;
  c7_J[8] = c7_h_y;
  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, -36);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, 13);
  c7_i10 = 0;
  c7_i11 = 0;
  for (c7_i12 = 0; c7_i12 < 3; c7_i12++) {
    for (c7_i13 = 0; c7_i13 < 3; c7_i13++) {
      c7_b_Phi[c7_i13 + c7_i10] = c7_Phi[(c7_i13 + c7_i11) + 36];
    }

    c7_i10 += 3;
    c7_i11 += 12;
  }

  for (c7_i14 = 0; c7_i14 < 9; c7_i14++) {
    c7_b_J[c7_i14] = c7_J[c7_i14];
  }

  c7_fn_Create_Q_r11(chartInstance, c7_b_Phi, c7_b_J, c7_t_delta, c7_tau, c7_b);
  for (c7_i15 = 0; c7_i15 < 9; c7_i15++) {
    c7_Q_r11[c7_i15] = 0.25 * c7_b[c7_i15];
  }

  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, 14);
  c7_i16 = 0;
  c7_i17 = 0;
  for (c7_i18 = 0; c7_i18 < 3; c7_i18++) {
    for (c7_i19 = 0; c7_i19 < 3; c7_i19++) {
      c7_phi_r12[c7_i19 + c7_i16] = c7_Phi[(c7_i19 + c7_i17) + 36];
    }

    c7_i16 += 3;
    c7_i17 += 12;
  }

  c7_i20 = 0;
  c7_i21 = 0;
  for (c7_i22 = 0; c7_i22 < 3; c7_i22++) {
    for (c7_i23 = 0; c7_i23 < 3; c7_i23++) {
      c7_phi_r22[c7_i23 + c7_i20] = c7_Phi[(c7_i23 + c7_i21) + 39];
    }

    c7_i20 += 3;
    c7_i21 += 12;
  }

  for (c7_i24 = 0; c7_i24 < 9; c7_i24++) {
    c7_c_J[c7_i24] = c7_J[c7_i24];
  }

  c7_b_t_delta = c7_t_delta;
  c7_b_tau = c7_tau;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c7_d_debug_family_names,
    c7_c_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_c_nargin, 0U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_c_nargout, 1U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_phi_r12, 2U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_phi_r22, 3U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_c_J, 4U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_b_t_delta, 5U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_b_tau, 6U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_b_Q_r12, 7U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  CV_EML_FCN(0, 2);
  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, 25);
  c7_a = c7_b_tau;
  for (c7_i25 = 0; c7_i25 < 9; c7_i25++) {
    c7_b[c7_i25] = c7_phi_r12[c7_i25];
  }

  for (c7_i26 = 0; c7_i26 < 9; c7_i26++) {
    c7_b[c7_i26] *= c7_a;
  }

  for (c7_i27 = 0; c7_i27 < 9; c7_i27++) {
    c7_d_J[c7_i27] = c7_c_J[c7_i27];
  }

  c7_mpower(chartInstance, c7_d_J, c7_b_b);
  c7_eml_scalar_eg(chartInstance);
  c7_eml_scalar_eg(chartInstance);
  for (c7_i28 = 0; c7_i28 < 9; c7_i28++) {
    c7_i_y[c7_i28] = 0.0;
  }

  for (c7_i29 = 0; c7_i29 < 9; c7_i29++) {
    c7_c_b[c7_i29] = c7_b[c7_i29];
  }

  for (c7_i30 = 0; c7_i30 < 9; c7_i30++) {
    c7_d_b[c7_i30] = c7_b_b[c7_i30];
  }

  c7_b_eml_xgemm(chartInstance, c7_c_b, c7_d_b, c7_i_y);
  for (c7_i31 = 0; c7_i31 < 9; c7_i31++) {
    c7_b[c7_i31] = c7_phi_r22[c7_i31];
  }

  c7_eml_scalar_eg(chartInstance);
  c7_eml_scalar_eg(chartInstance);
  for (c7_i32 = 0; c7_i32 < 9; c7_i32++) {
    c7_b_b[c7_i32] = 0.0;
  }

  for (c7_i33 = 0; c7_i33 < 9; c7_i33++) {
    c7_j_y[c7_i33] = c7_i_y[c7_i33];
  }

  for (c7_i34 = 0; c7_i34 < 9; c7_i34++) {
    c7_e_b[c7_i34] = c7_b[c7_i34];
  }

  c7_b_eml_xgemm(chartInstance, c7_j_y, c7_e_b, c7_b_b);
  c7_f_b = c7_b_t_delta;
  for (c7_i35 = 0; c7_i35 < 9; c7_i35++) {
    c7_b_Q_r12[c7_i35] = c7_b_b[c7_i35] * c7_f_b;
  }

  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, -25);
  _SFD_SYMBOL_SCOPE_POP();
  for (c7_i36 = 0; c7_i36 < 9; c7_i36++) {
    c7_b[c7_i36] = c7_b_Q_r12[c7_i36];
  }

  for (c7_i37 = 0; c7_i37 < 9; c7_i37++) {
    c7_Q_r12[c7_i37] = 0.5 * c7_b[c7_i37];
  }

  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, 15);
  c7_i38 = 0;
  c7_i39 = 0;
  for (c7_i40 = 0; c7_i40 < 3; c7_i40++) {
    for (c7_i41 = 0; c7_i41 < 3; c7_i41++) {
      c7_b_phi_r22[c7_i41 + c7_i38] = c7_Phi[(c7_i41 + c7_i39) + 39];
    }

    c7_i38 += 3;
    c7_i39 += 12;
  }

  for (c7_i42 = 0; c7_i42 < 9; c7_i42++) {
    c7_e_J[c7_i42] = c7_J[c7_i42];
  }

  c7_c_t_delta = c7_t_delta;
  c7_c_tau = c7_tau;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 7U, 7U, c7_e_debug_family_names,
    c7_d_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_d_nargin, 0U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_d_nargout, 1U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_b_phi_r22, 2U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_e_J, 3U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_c_t_delta, 4U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_c_tau, 5U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_Q_r22, 6U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  CV_EML_FCN(0, 3);
  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, 28);
  c7_b_a = c7_c_tau;
  for (c7_i43 = 0; c7_i43 < 9; c7_i43++) {
    c7_b[c7_i43] = c7_b_phi_r22[c7_i43];
  }

  for (c7_i44 = 0; c7_i44 < 9; c7_i44++) {
    c7_b[c7_i44] *= c7_b_a;
  }

  for (c7_i45 = 0; c7_i45 < 9; c7_i45++) {
    c7_f_J[c7_i45] = c7_e_J[c7_i45];
  }

  c7_mpower(chartInstance, c7_f_J, c7_b_b);
  c7_eml_scalar_eg(chartInstance);
  c7_eml_scalar_eg(chartInstance);
  for (c7_i46 = 0; c7_i46 < 9; c7_i46++) {
    c7_i_y[c7_i46] = 0.0;
  }

  for (c7_i47 = 0; c7_i47 < 9; c7_i47++) {
    c7_g_b[c7_i47] = c7_b[c7_i47];
  }

  for (c7_i48 = 0; c7_i48 < 9; c7_i48++) {
    c7_h_b[c7_i48] = c7_b_b[c7_i48];
  }

  c7_b_eml_xgemm(chartInstance, c7_g_b, c7_h_b, c7_i_y);
  for (c7_i49 = 0; c7_i49 < 9; c7_i49++) {
    c7_b[c7_i49] = c7_b_phi_r22[c7_i49];
  }

  c7_eml_scalar_eg(chartInstance);
  c7_eml_scalar_eg(chartInstance);
  for (c7_i50 = 0; c7_i50 < 9; c7_i50++) {
    c7_b_b[c7_i50] = 0.0;
  }

  for (c7_i51 = 0; c7_i51 < 9; c7_i51++) {
    c7_k_y[c7_i51] = c7_i_y[c7_i51];
  }

  for (c7_i52 = 0; c7_i52 < 9; c7_i52++) {
    c7_i_b[c7_i52] = c7_b[c7_i52];
  }

  c7_b_eml_xgemm(chartInstance, c7_k_y, c7_i_b, c7_b_b);
  c7_j_b = c7_c_t_delta;
  for (c7_i53 = 0; c7_i53 < 9; c7_i53++) {
    c7_Q_r22[c7_i53] = c7_b_b[c7_i53] * c7_j_b;
  }

  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, -28);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, 16);
  c7_i54 = 0;
  c7_i55 = 0;
  for (c7_i56 = 0; c7_i56 < 3; c7_i56++) {
    for (c7_i57 = 0; c7_i57 < 3; c7_i57++) {
      c7_Q_r[c7_i57 + c7_i54] = c7_Q_r11[c7_i57 + c7_i55];
    }

    c7_i54 += 6;
    c7_i55 += 3;
  }

  c7_i58 = 0;
  c7_i59 = 0;
  for (c7_i60 = 0; c7_i60 < 3; c7_i60++) {
    for (c7_i61 = 0; c7_i61 < 3; c7_i61++) {
      c7_Q_r[(c7_i61 + c7_i58) + 18] = c7_Q_r12[c7_i61 + c7_i59];
    }

    c7_i58 += 6;
    c7_i59 += 3;
  }

  c7_i62 = 0;
  for (c7_i63 = 0; c7_i63 < 3; c7_i63++) {
    c7_i64 = 0;
    for (c7_i65 = 0; c7_i65 < 3; c7_i65++) {
      c7_Q_r[(c7_i65 + c7_i62) + 3] = c7_Q_r12[c7_i64 + c7_i63];
      c7_i64 += 3;
    }

    c7_i62 += 6;
  }

  c7_i66 = 0;
  c7_i67 = 0;
  for (c7_i68 = 0; c7_i68 < 3; c7_i68++) {
    for (c7_i69 = 0; c7_i69 < 3; c7_i69++) {
      c7_Q_r[(c7_i69 + c7_i66) + 21] = c7_Q_r22[c7_i69 + c7_i67];
    }

    c7_i66 += 6;
    c7_i67 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, -16);
  _SFD_SYMBOL_SCOPE_POP();
  for (c7_i70 = 0; c7_i70 < 36; c7_i70++) {
    (*c7_b_Q_r)[c7_i70] = c7_Q_r[c7_i70];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 6U, chartInstance->c7_sfEvent);
}

static void initSimStructsc7_Model_02(SFc7_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c7_machineNumber, uint32_T
  c7_chartNumber, uint32_T c7_instanceNumber)
{
  (void)c7_machineNumber;
  (void)c7_chartNumber;
  (void)c7_instanceNumber;
}

static const mxArray *c7_sf_marshallOut(void *chartInstanceVoid, void *c7_inData)
{
  const mxArray *c7_mxArrayOutData = NULL;
  int32_T c7_i71;
  int32_T c7_i72;
  int32_T c7_i73;
  real_T c7_b_inData[36];
  int32_T c7_i74;
  int32_T c7_i75;
  int32_T c7_i76;
  real_T c7_u[36];
  const mxArray *c7_y = NULL;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_mxArrayOutData = NULL;
  c7_i71 = 0;
  for (c7_i72 = 0; c7_i72 < 6; c7_i72++) {
    for (c7_i73 = 0; c7_i73 < 6; c7_i73++) {
      c7_b_inData[c7_i73 + c7_i71] = (*(real_T (*)[36])c7_inData)[c7_i73 +
        c7_i71];
    }

    c7_i71 += 6;
  }

  c7_i74 = 0;
  for (c7_i75 = 0; c7_i75 < 6; c7_i75++) {
    for (c7_i76 = 0; c7_i76 < 6; c7_i76++) {
      c7_u[c7_i76 + c7_i74] = c7_b_inData[c7_i76 + c7_i74];
    }

    c7_i74 += 6;
  }

  c7_y = NULL;
  sf_mex_assign(&c7_y, sf_mex_create("y", c7_u, 0, 0U, 1U, 0U, 2, 6, 6), false);
  sf_mex_assign(&c7_mxArrayOutData, c7_y, false);
  return c7_mxArrayOutData;
}

static void c7_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_Q_r, const char_T *c7_identifier, real_T c7_y[36])
{
  emlrtMsgIdentifier c7_thisId;
  c7_thisId.fIdentifier = c7_identifier;
  c7_thisId.fParent = NULL;
  c7_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c7_Q_r), &c7_thisId, c7_y);
  sf_mex_destroy(&c7_Q_r);
}

static void c7_b_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId, real_T c7_y[36])
{
  real_T c7_dv2[36];
  int32_T c7_i77;
  (void)chartInstance;
  sf_mex_import(c7_parentId, sf_mex_dup(c7_u), c7_dv2, 1, 0, 0U, 1, 0U, 2, 6, 6);
  for (c7_i77 = 0; c7_i77 < 36; c7_i77++) {
    c7_y[c7_i77] = c7_dv2[c7_i77];
  }

  sf_mex_destroy(&c7_u);
}

static void c7_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData)
{
  const mxArray *c7_Q_r;
  const char_T *c7_identifier;
  emlrtMsgIdentifier c7_thisId;
  real_T c7_y[36];
  int32_T c7_i78;
  int32_T c7_i79;
  int32_T c7_i80;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_Q_r = sf_mex_dup(c7_mxArrayInData);
  c7_identifier = c7_varName;
  c7_thisId.fIdentifier = c7_identifier;
  c7_thisId.fParent = NULL;
  c7_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c7_Q_r), &c7_thisId, c7_y);
  sf_mex_destroy(&c7_Q_r);
  c7_i78 = 0;
  for (c7_i79 = 0; c7_i79 < 6; c7_i79++) {
    for (c7_i80 = 0; c7_i80 < 6; c7_i80++) {
      (*(real_T (*)[36])c7_outData)[c7_i80 + c7_i78] = c7_y[c7_i80 + c7_i78];
    }

    c7_i78 += 6;
  }

  sf_mex_destroy(&c7_mxArrayInData);
}

static const mxArray *c7_b_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData)
{
  const mxArray *c7_mxArrayOutData = NULL;
  int32_T c7_i81;
  real_T c7_b_inData[3];
  int32_T c7_i82;
  real_T c7_u[3];
  const mxArray *c7_y = NULL;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_mxArrayOutData = NULL;
  for (c7_i81 = 0; c7_i81 < 3; c7_i81++) {
    c7_b_inData[c7_i81] = (*(real_T (*)[3])c7_inData)[c7_i81];
  }

  for (c7_i82 = 0; c7_i82 < 3; c7_i82++) {
    c7_u[c7_i82] = c7_b_inData[c7_i82];
  }

  c7_y = NULL;
  sf_mex_assign(&c7_y, sf_mex_create("y", c7_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c7_mxArrayOutData, c7_y, false);
  return c7_mxArrayOutData;
}

static const mxArray *c7_c_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData)
{
  const mxArray *c7_mxArrayOutData = NULL;
  real_T c7_u;
  const mxArray *c7_y = NULL;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_mxArrayOutData = NULL;
  c7_u = *(real_T *)c7_inData;
  c7_y = NULL;
  sf_mex_assign(&c7_y, sf_mex_create("y", &c7_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c7_mxArrayOutData, c7_y, false);
  return c7_mxArrayOutData;
}

static const mxArray *c7_d_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData)
{
  const mxArray *c7_mxArrayOutData = NULL;
  int32_T c7_i83;
  int32_T c7_i84;
  int32_T c7_i85;
  real_T c7_b_inData[144];
  int32_T c7_i86;
  int32_T c7_i87;
  int32_T c7_i88;
  real_T c7_u[144];
  const mxArray *c7_y = NULL;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_mxArrayOutData = NULL;
  c7_i83 = 0;
  for (c7_i84 = 0; c7_i84 < 12; c7_i84++) {
    for (c7_i85 = 0; c7_i85 < 12; c7_i85++) {
      c7_b_inData[c7_i85 + c7_i83] = (*(real_T (*)[144])c7_inData)[c7_i85 +
        c7_i83];
    }

    c7_i83 += 12;
  }

  c7_i86 = 0;
  for (c7_i87 = 0; c7_i87 < 12; c7_i87++) {
    for (c7_i88 = 0; c7_i88 < 12; c7_i88++) {
      c7_u[c7_i88 + c7_i86] = c7_b_inData[c7_i88 + c7_i86];
    }

    c7_i86 += 12;
  }

  c7_y = NULL;
  sf_mex_assign(&c7_y, sf_mex_create("y", c7_u, 0, 0U, 1U, 0U, 2, 12, 12), false);
  sf_mex_assign(&c7_mxArrayOutData, c7_y, false);
  return c7_mxArrayOutData;
}

static real_T c7_c_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId)
{
  real_T c7_y;
  real_T c7_d0;
  (void)chartInstance;
  sf_mex_import(c7_parentId, sf_mex_dup(c7_u), &c7_d0, 1, 0, 0U, 0, 0U, 0);
  c7_y = c7_d0;
  sf_mex_destroy(&c7_u);
  return c7_y;
}

static void c7_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData)
{
  const mxArray *c7_nargout;
  const char_T *c7_identifier;
  emlrtMsgIdentifier c7_thisId;
  real_T c7_y;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_nargout = sf_mex_dup(c7_mxArrayInData);
  c7_identifier = c7_varName;
  c7_thisId.fIdentifier = c7_identifier;
  c7_thisId.fParent = NULL;
  c7_y = c7_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c7_nargout), &c7_thisId);
  sf_mex_destroy(&c7_nargout);
  *(real_T *)c7_outData = c7_y;
  sf_mex_destroy(&c7_mxArrayInData);
}

static const mxArray *c7_e_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData)
{
  const mxArray *c7_mxArrayOutData = NULL;
  int32_T c7_i89;
  int32_T c7_i90;
  int32_T c7_i91;
  real_T c7_b_inData[9];
  int32_T c7_i92;
  int32_T c7_i93;
  int32_T c7_i94;
  real_T c7_u[9];
  const mxArray *c7_y = NULL;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_mxArrayOutData = NULL;
  c7_i89 = 0;
  for (c7_i90 = 0; c7_i90 < 3; c7_i90++) {
    for (c7_i91 = 0; c7_i91 < 3; c7_i91++) {
      c7_b_inData[c7_i91 + c7_i89] = (*(real_T (*)[9])c7_inData)[c7_i91 + c7_i89];
    }

    c7_i89 += 3;
  }

  c7_i92 = 0;
  for (c7_i93 = 0; c7_i93 < 3; c7_i93++) {
    for (c7_i94 = 0; c7_i94 < 3; c7_i94++) {
      c7_u[c7_i94 + c7_i92] = c7_b_inData[c7_i94 + c7_i92];
    }

    c7_i92 += 3;
  }

  c7_y = NULL;
  sf_mex_assign(&c7_y, sf_mex_create("y", c7_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c7_mxArrayOutData, c7_y, false);
  return c7_mxArrayOutData;
}

static void c7_d_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId, real_T c7_y[9])
{
  real_T c7_dv3[9];
  int32_T c7_i95;
  (void)chartInstance;
  sf_mex_import(c7_parentId, sf_mex_dup(c7_u), c7_dv3, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c7_i95 = 0; c7_i95 < 9; c7_i95++) {
    c7_y[c7_i95] = c7_dv3[c7_i95];
  }

  sf_mex_destroy(&c7_u);
}

static void c7_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData)
{
  const mxArray *c7_Q_r22;
  const char_T *c7_identifier;
  emlrtMsgIdentifier c7_thisId;
  real_T c7_y[9];
  int32_T c7_i96;
  int32_T c7_i97;
  int32_T c7_i98;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_Q_r22 = sf_mex_dup(c7_mxArrayInData);
  c7_identifier = c7_varName;
  c7_thisId.fIdentifier = c7_identifier;
  c7_thisId.fParent = NULL;
  c7_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c7_Q_r22), &c7_thisId, c7_y);
  sf_mex_destroy(&c7_Q_r22);
  c7_i96 = 0;
  for (c7_i97 = 0; c7_i97 < 3; c7_i97++) {
    for (c7_i98 = 0; c7_i98 < 3; c7_i98++) {
      (*(real_T (*)[9])c7_outData)[c7_i98 + c7_i96] = c7_y[c7_i98 + c7_i96];
    }

    c7_i96 += 3;
  }

  sf_mex_destroy(&c7_mxArrayInData);
}

static void c7_e_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId, real_T c7_y[3])
{
  real_T c7_dv4[3];
  int32_T c7_i99;
  (void)chartInstance;
  sf_mex_import(c7_parentId, sf_mex_dup(c7_u), c7_dv4, 1, 0, 0U, 1, 0U, 1, 3);
  for (c7_i99 = 0; c7_i99 < 3; c7_i99++) {
    c7_y[c7_i99] = c7_dv4[c7_i99];
  }

  sf_mex_destroy(&c7_u);
}

static void c7_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData)
{
  const mxArray *c7_p;
  const char_T *c7_identifier;
  emlrtMsgIdentifier c7_thisId;
  real_T c7_y[3];
  int32_T c7_i100;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_p = sf_mex_dup(c7_mxArrayInData);
  c7_identifier = c7_varName;
  c7_thisId.fIdentifier = c7_identifier;
  c7_thisId.fParent = NULL;
  c7_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c7_p), &c7_thisId, c7_y);
  sf_mex_destroy(&c7_p);
  for (c7_i100 = 0; c7_i100 < 3; c7_i100++) {
    (*(real_T (*)[3])c7_outData)[c7_i100] = c7_y[c7_i100];
  }

  sf_mex_destroy(&c7_mxArrayInData);
}

const mxArray *sf_c7_Model_02_get_eml_resolved_functions_info(void)
{
  const mxArray *c7_nameCaptureInfo = NULL;
  c7_nameCaptureInfo = NULL;
  sf_mex_assign(&c7_nameCaptureInfo, sf_mex_createstruct("structure", 2, 26, 1),
                false);
  c7_info_helper(&c7_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c7_nameCaptureInfo);
  return c7_nameCaptureInfo;
}

static void c7_info_helper(const mxArray **c7_info)
{
  const mxArray *c7_rhs0 = NULL;
  const mxArray *c7_lhs0 = NULL;
  const mxArray *c7_rhs1 = NULL;
  const mxArray *c7_lhs1 = NULL;
  const mxArray *c7_rhs2 = NULL;
  const mxArray *c7_lhs2 = NULL;
  const mxArray *c7_rhs3 = NULL;
  const mxArray *c7_lhs3 = NULL;
  const mxArray *c7_rhs4 = NULL;
  const mxArray *c7_lhs4 = NULL;
  const mxArray *c7_rhs5 = NULL;
  const mxArray *c7_lhs5 = NULL;
  const mxArray *c7_rhs6 = NULL;
  const mxArray *c7_lhs6 = NULL;
  const mxArray *c7_rhs7 = NULL;
  const mxArray *c7_lhs7 = NULL;
  const mxArray *c7_rhs8 = NULL;
  const mxArray *c7_lhs8 = NULL;
  const mxArray *c7_rhs9 = NULL;
  const mxArray *c7_lhs9 = NULL;
  const mxArray *c7_rhs10 = NULL;
  const mxArray *c7_lhs10 = NULL;
  const mxArray *c7_rhs11 = NULL;
  const mxArray *c7_lhs11 = NULL;
  const mxArray *c7_rhs12 = NULL;
  const mxArray *c7_lhs12 = NULL;
  const mxArray *c7_rhs13 = NULL;
  const mxArray *c7_lhs13 = NULL;
  const mxArray *c7_rhs14 = NULL;
  const mxArray *c7_lhs14 = NULL;
  const mxArray *c7_rhs15 = NULL;
  const mxArray *c7_lhs15 = NULL;
  const mxArray *c7_rhs16 = NULL;
  const mxArray *c7_lhs16 = NULL;
  const mxArray *c7_rhs17 = NULL;
  const mxArray *c7_lhs17 = NULL;
  const mxArray *c7_rhs18 = NULL;
  const mxArray *c7_lhs18 = NULL;
  const mxArray *c7_rhs19 = NULL;
  const mxArray *c7_lhs19 = NULL;
  const mxArray *c7_rhs20 = NULL;
  const mxArray *c7_lhs20 = NULL;
  const mxArray *c7_rhs21 = NULL;
  const mxArray *c7_lhs21 = NULL;
  const mxArray *c7_rhs22 = NULL;
  const mxArray *c7_lhs22 = NULL;
  const mxArray *c7_rhs23 = NULL;
  const mxArray *c7_lhs23 = NULL;
  const mxArray *c7_rhs24 = NULL;
  const mxArray *c7_lhs24 = NULL;
  const mxArray *c7_rhs25 = NULL;
  const mxArray *c7_lhs25 = NULL;
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("mrdivide"), "name", "name", 0);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c7_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 1);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 1);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c7_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 2);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("rdivide"), "name", "name", 2);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c7_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c7_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 4);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c7_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_div"), "name", "name", 5);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c7_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 6);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c7_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(""), "context", "context", 7);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 7);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c7_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 8);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 8);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c7_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(""), "context", "context", 9);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("mpower"), "name", "name", 9);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c7_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 10);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c7_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("ismatrix"), "name", "name", 11);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c7_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 12);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c7_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 13);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 13);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c7_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 14);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c7_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 15);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 15);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c7_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 16);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c7_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 17);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c7_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 18);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  18);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c7_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 19);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c7_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 20);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c7_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 21);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 21);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c7_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 22);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 22);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c7_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 23);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 23);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c7_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 24);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 24);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c7_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 25);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 25);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c7_info, c7_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c7_info, c7_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c7_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c7_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c7_info, sf_mex_duplicatearraysafe(&c7_lhs25), "lhs", "lhs",
                  25);
  sf_mex_destroy(&c7_rhs0);
  sf_mex_destroy(&c7_lhs0);
  sf_mex_destroy(&c7_rhs1);
  sf_mex_destroy(&c7_lhs1);
  sf_mex_destroy(&c7_rhs2);
  sf_mex_destroy(&c7_lhs2);
  sf_mex_destroy(&c7_rhs3);
  sf_mex_destroy(&c7_lhs3);
  sf_mex_destroy(&c7_rhs4);
  sf_mex_destroy(&c7_lhs4);
  sf_mex_destroy(&c7_rhs5);
  sf_mex_destroy(&c7_lhs5);
  sf_mex_destroy(&c7_rhs6);
  sf_mex_destroy(&c7_lhs6);
  sf_mex_destroy(&c7_rhs7);
  sf_mex_destroy(&c7_lhs7);
  sf_mex_destroy(&c7_rhs8);
  sf_mex_destroy(&c7_lhs8);
  sf_mex_destroy(&c7_rhs9);
  sf_mex_destroy(&c7_lhs9);
  sf_mex_destroy(&c7_rhs10);
  sf_mex_destroy(&c7_lhs10);
  sf_mex_destroy(&c7_rhs11);
  sf_mex_destroy(&c7_lhs11);
  sf_mex_destroy(&c7_rhs12);
  sf_mex_destroy(&c7_lhs12);
  sf_mex_destroy(&c7_rhs13);
  sf_mex_destroy(&c7_lhs13);
  sf_mex_destroy(&c7_rhs14);
  sf_mex_destroy(&c7_lhs14);
  sf_mex_destroy(&c7_rhs15);
  sf_mex_destroy(&c7_lhs15);
  sf_mex_destroy(&c7_rhs16);
  sf_mex_destroy(&c7_lhs16);
  sf_mex_destroy(&c7_rhs17);
  sf_mex_destroy(&c7_lhs17);
  sf_mex_destroy(&c7_rhs18);
  sf_mex_destroy(&c7_lhs18);
  sf_mex_destroy(&c7_rhs19);
  sf_mex_destroy(&c7_lhs19);
  sf_mex_destroy(&c7_rhs20);
  sf_mex_destroy(&c7_lhs20);
  sf_mex_destroy(&c7_rhs21);
  sf_mex_destroy(&c7_lhs21);
  sf_mex_destroy(&c7_rhs22);
  sf_mex_destroy(&c7_lhs22);
  sf_mex_destroy(&c7_rhs23);
  sf_mex_destroy(&c7_lhs23);
  sf_mex_destroy(&c7_rhs24);
  sf_mex_destroy(&c7_lhs24);
  sf_mex_destroy(&c7_rhs25);
  sf_mex_destroy(&c7_lhs25);
}

static const mxArray *c7_emlrt_marshallOut(const char * c7_u)
{
  const mxArray *c7_y = NULL;
  c7_y = NULL;
  sf_mex_assign(&c7_y, sf_mex_create("y", c7_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c7_u)), false);
  return c7_y;
}

static const mxArray *c7_b_emlrt_marshallOut(const uint32_T c7_u)
{
  const mxArray *c7_y = NULL;
  c7_y = NULL;
  sf_mex_assign(&c7_y, sf_mex_create("y", &c7_u, 7, 0U, 0U, 0U, 0), false);
  return c7_y;
}

static void c7_fn_Create_Q_r11(SFc7_Model_02InstanceStruct *chartInstance,
  real_T c7_phi_r12[9], real_T c7_J[9], real_T c7_t_delta, real_T c7_tau, real_T
  c7_Q_r11[9])
{
  uint32_T c7_debug_family_var_map[7];
  real_T c7_nargin = 4.0;
  real_T c7_nargout = 1.0;
  real_T c7_a;
  int32_T c7_i101;
  real_T c7_b[9];
  int32_T c7_i102;
  int32_T c7_i103;
  real_T c7_b_J[9];
  real_T c7_b_b[9];
  int32_T c7_i104;
  real_T c7_y[9];
  int32_T c7_i105;
  real_T c7_c_b[9];
  int32_T c7_i106;
  real_T c7_d_b[9];
  int32_T c7_i107;
  int32_T c7_i108;
  int32_T c7_i109;
  int32_T c7_i110;
  int32_T c7_i111;
  int32_T c7_i112;
  real_T c7_b_y[9];
  int32_T c7_i113;
  real_T c7_e_b[9];
  real_T c7_f_b;
  int32_T c7_i114;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 7U, 7U, c7_c_debug_family_names,
    c7_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_nargin, 0U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_nargout, 1U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_phi_r12, 2U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_J, 3U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_t_delta, 4U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c7_tau, 5U, c7_c_sf_marshallOut,
    c7_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c7_Q_r11, 6U, c7_e_sf_marshallOut,
    c7_c_sf_marshallIn);
  CV_EML_FCN(0, 1);
  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, 21);
  c7_a = c7_tau;
  for (c7_i101 = 0; c7_i101 < 9; c7_i101++) {
    c7_b[c7_i101] = c7_phi_r12[c7_i101];
  }

  for (c7_i102 = 0; c7_i102 < 9; c7_i102++) {
    c7_b[c7_i102] *= c7_a;
  }

  for (c7_i103 = 0; c7_i103 < 9; c7_i103++) {
    c7_b_J[c7_i103] = c7_J[c7_i103];
  }

  c7_mpower(chartInstance, c7_b_J, c7_b_b);
  c7_eml_scalar_eg(chartInstance);
  c7_eml_scalar_eg(chartInstance);
  for (c7_i104 = 0; c7_i104 < 9; c7_i104++) {
    c7_y[c7_i104] = 0.0;
  }

  for (c7_i105 = 0; c7_i105 < 9; c7_i105++) {
    c7_c_b[c7_i105] = c7_b[c7_i105];
  }

  for (c7_i106 = 0; c7_i106 < 9; c7_i106++) {
    c7_d_b[c7_i106] = c7_b_b[c7_i106];
  }

  c7_b_eml_xgemm(chartInstance, c7_c_b, c7_d_b, c7_y);
  c7_i107 = 0;
  for (c7_i108 = 0; c7_i108 < 3; c7_i108++) {
    c7_i109 = 0;
    for (c7_i110 = 0; c7_i110 < 3; c7_i110++) {
      c7_b[c7_i110 + c7_i107] = c7_phi_r12[c7_i109 + c7_i108];
      c7_i109 += 3;
    }

    c7_i107 += 3;
  }

  c7_eml_scalar_eg(chartInstance);
  c7_eml_scalar_eg(chartInstance);
  for (c7_i111 = 0; c7_i111 < 9; c7_i111++) {
    c7_b_b[c7_i111] = 0.0;
  }

  for (c7_i112 = 0; c7_i112 < 9; c7_i112++) {
    c7_b_y[c7_i112] = c7_y[c7_i112];
  }

  for (c7_i113 = 0; c7_i113 < 9; c7_i113++) {
    c7_e_b[c7_i113] = c7_b[c7_i113];
  }

  c7_b_eml_xgemm(chartInstance, c7_b_y, c7_e_b, c7_b_b);
  c7_f_b = c7_t_delta;
  for (c7_i114 = 0; c7_i114 < 9; c7_i114++) {
    c7_Q_r11[c7_i114] = c7_b_b[c7_i114] * c7_f_b;
  }

  _SFD_EML_CALL(0U, chartInstance->c7_sfEvent, -21);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c7_mpower(SFc7_Model_02InstanceStruct *chartInstance, real_T c7_a[9],
                      real_T c7_c[9])
{
  int32_T c7_i115;
  int32_T c7_i116;
  real_T c7_b_a[9];
  int32_T c7_i117;
  real_T c7_c_a[9];
  c7_eml_scalar_eg(chartInstance);
  c7_eml_scalar_eg(chartInstance);
  for (c7_i115 = 0; c7_i115 < 9; c7_i115++) {
    c7_c[c7_i115] = 0.0;
  }

  for (c7_i116 = 0; c7_i116 < 9; c7_i116++) {
    c7_b_a[c7_i116] = c7_a[c7_i116];
  }

  for (c7_i117 = 0; c7_i117 < 9; c7_i117++) {
    c7_c_a[c7_i117] = c7_a[c7_i117];
  }

  c7_b_eml_xgemm(chartInstance, c7_b_a, c7_c_a, c7_c);
}

static void c7_eml_scalar_eg(SFc7_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c7_eml_xgemm(SFc7_Model_02InstanceStruct *chartInstance, real_T
  c7_A[9], real_T c7_B[9], real_T c7_C[9], real_T c7_b_C[9])
{
  int32_T c7_i118;
  int32_T c7_i119;
  real_T c7_b_A[9];
  int32_T c7_i120;
  real_T c7_b_B[9];
  for (c7_i118 = 0; c7_i118 < 9; c7_i118++) {
    c7_b_C[c7_i118] = c7_C[c7_i118];
  }

  for (c7_i119 = 0; c7_i119 < 9; c7_i119++) {
    c7_b_A[c7_i119] = c7_A[c7_i119];
  }

  for (c7_i120 = 0; c7_i120 < 9; c7_i120++) {
    c7_b_B[c7_i120] = c7_B[c7_i120];
  }

  c7_b_eml_xgemm(chartInstance, c7_b_A, c7_b_B, c7_b_C);
}

static const mxArray *c7_f_sf_marshallOut(void *chartInstanceVoid, void
  *c7_inData)
{
  const mxArray *c7_mxArrayOutData = NULL;
  int32_T c7_u;
  const mxArray *c7_y = NULL;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_mxArrayOutData = NULL;
  c7_u = *(int32_T *)c7_inData;
  c7_y = NULL;
  sf_mex_assign(&c7_y, sf_mex_create("y", &c7_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c7_mxArrayOutData, c7_y, false);
  return c7_mxArrayOutData;
}

static int32_T c7_f_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId)
{
  int32_T c7_y;
  int32_T c7_i121;
  (void)chartInstance;
  sf_mex_import(c7_parentId, sf_mex_dup(c7_u), &c7_i121, 1, 6, 0U, 0, 0U, 0);
  c7_y = c7_i121;
  sf_mex_destroy(&c7_u);
  return c7_y;
}

static void c7_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c7_mxArrayInData, const char_T *c7_varName, void *c7_outData)
{
  const mxArray *c7_b_sfEvent;
  const char_T *c7_identifier;
  emlrtMsgIdentifier c7_thisId;
  int32_T c7_y;
  SFc7_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc7_Model_02InstanceStruct *)chartInstanceVoid;
  c7_b_sfEvent = sf_mex_dup(c7_mxArrayInData);
  c7_identifier = c7_varName;
  c7_thisId.fIdentifier = c7_identifier;
  c7_thisId.fParent = NULL;
  c7_y = c7_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c7_b_sfEvent),
    &c7_thisId);
  sf_mex_destroy(&c7_b_sfEvent);
  *(int32_T *)c7_outData = c7_y;
  sf_mex_destroy(&c7_mxArrayInData);
}

static uint8_T c7_g_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_b_is_active_c7_Model_02, const char_T *c7_identifier)
{
  uint8_T c7_y;
  emlrtMsgIdentifier c7_thisId;
  c7_thisId.fIdentifier = c7_identifier;
  c7_thisId.fParent = NULL;
  c7_y = c7_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c7_b_is_active_c7_Model_02), &c7_thisId);
  sf_mex_destroy(&c7_b_is_active_c7_Model_02);
  return c7_y;
}

static uint8_T c7_h_emlrt_marshallIn(SFc7_Model_02InstanceStruct *chartInstance,
  const mxArray *c7_u, const emlrtMsgIdentifier *c7_parentId)
{
  uint8_T c7_y;
  uint8_T c7_u0;
  (void)chartInstance;
  sf_mex_import(c7_parentId, sf_mex_dup(c7_u), &c7_u0, 1, 3, 0U, 0, 0U, 0);
  c7_y = c7_u0;
  sf_mex_destroy(&c7_u);
  return c7_y;
}

static void c7_b_eml_xgemm(SFc7_Model_02InstanceStruct *chartInstance, real_T
  c7_A[9], real_T c7_B[9], real_T c7_C[9])
{
  int32_T c7_i122;
  int32_T c7_i123;
  int32_T c7_i124;
  int32_T c7_i125;
  int32_T c7_i126;
  (void)chartInstance;
  for (c7_i122 = 0; c7_i122 < 3; c7_i122++) {
    c7_i123 = 0;
    for (c7_i124 = 0; c7_i124 < 3; c7_i124++) {
      c7_C[c7_i123 + c7_i122] = 0.0;
      c7_i125 = 0;
      for (c7_i126 = 0; c7_i126 < 3; c7_i126++) {
        c7_C[c7_i123 + c7_i122] += c7_A[c7_i125 + c7_i122] * c7_B[c7_i126 +
          c7_i123];
        c7_i125 += 3;
      }

      c7_i123 += 3;
    }
  }
}

static void init_dsm_address_info(SFc7_Model_02InstanceStruct *chartInstance)
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

void sf_c7_Model_02_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(928989769U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(152734646U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3757477034U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2311851141U);
}

mxArray *sf_c7_Model_02_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("YCUUo4NLtqpoc8IBQ9KNGB");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,6,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(12);
      pr[1] = (double)(12);
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
      pr[0] = (double)(1);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));
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
      pr[1] = (double)(6);
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

mxArray *sf_c7_Model_02_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c7_Model_02_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c7_Model_02(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"Q_r\",},{M[8],M[0],T\"is_active_c7_Model_02\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c7_Model_02_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc7_Model_02InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc7_Model_02InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_02MachineNumber_,
           7,
           1,
           1,
           0,
           7,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation(_Model_02MachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_Model_02MachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _Model_02MachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"Phi");
          _SFD_SET_DATA_PROPS(1,2,0,1,"Q_r");
          _SFD_SET_DATA_PROPS(2,1,1,0,"unused1");
          _SFD_SET_DATA_PROPS(3,1,1,0,"t_delta");
          _SFD_SET_DATA_PROPS(4,1,1,0,"tau");
          _SFD_SET_DATA_PROPS(5,1,1,0,"unused");
          _SFD_SET_DATA_PROPS(6,1,1,0,"p");
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
        _SFD_CV_INIT_EML(0,1,5,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",343,-1,678);
        _SFD_CV_INIT_EML_FCN(0,1,"fn_Create_Q_r11",723,-1,832);
        _SFD_CV_INIT_EML_FCN(0,2,"fn_Create_Q_r12",834,-1,948);
        _SFD_CV_INIT_EML_FCN(0,3,"fn_Create_Q_r22",949,-1,1056);
        _SFD_CV_INIT_EML_FCN(0,4,"fn_Create_J",1212,-1,1304);

        {
          unsigned int dimVector[2];
          dimVector[0]= 12;
          dimVector[1]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c7_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c7_sf_marshallOut,(MexInFcnForType)
            c7_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c7_c_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c7_c_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c7_c_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c7_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c7_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c7_unused1;
          real_T *c7_t_delta;
          real_T *c7_tau;
          real_T *c7_unused;
          real_T (*c7_Phi)[144];
          real_T (*c7_Q_r)[36];
          real_T (*c7_p)[3];
          c7_p = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
          c7_unused = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c7_tau = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c7_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c7_unused1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c7_Q_r = (real_T (*)[36])ssGetOutputPortSignal(chartInstance->S, 1);
          c7_Phi = (real_T (*)[144])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c7_Phi);
          _SFD_SET_DATA_VALUE_PTR(1U, *c7_Q_r);
          _SFD_SET_DATA_VALUE_PTR(2U, c7_unused1);
          _SFD_SET_DATA_VALUE_PTR(3U, c7_t_delta);
          _SFD_SET_DATA_VALUE_PTR(4U, c7_tau);
          _SFD_SET_DATA_VALUE_PTR(5U, c7_unused);
          _SFD_SET_DATA_VALUE_PTR(6U, *c7_p);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _Model_02MachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "EforzPMyvPho6YHhxHakWC";
}

static void sf_opaque_initialize_c7_Model_02(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc7_Model_02InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c7_Model_02((SFc7_Model_02InstanceStruct*) chartInstanceVar);
  initialize_c7_Model_02((SFc7_Model_02InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c7_Model_02(void *chartInstanceVar)
{
  enable_c7_Model_02((SFc7_Model_02InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c7_Model_02(void *chartInstanceVar)
{
  disable_c7_Model_02((SFc7_Model_02InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c7_Model_02(void *chartInstanceVar)
{
  sf_gateway_c7_Model_02((SFc7_Model_02InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c7_Model_02(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c7_Model_02((SFc7_Model_02InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c7_Model_02();/* state var info */
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

extern void sf_internal_set_sim_state_c7_Model_02(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c7_Model_02();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c7_Model_02((SFc7_Model_02InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c7_Model_02(SimStruct* S)
{
  return sf_internal_get_sim_state_c7_Model_02(S);
}

static void sf_opaque_set_sim_state_c7_Model_02(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c7_Model_02(S, st);
}

static void sf_opaque_terminate_c7_Model_02(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc7_Model_02InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_02_optimization_info();
    }

    finalize_c7_Model_02((SFc7_Model_02InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc7_Model_02((SFc7_Model_02InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c7_Model_02(SimStruct *S)
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
    initialize_params_c7_Model_02((SFc7_Model_02InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c7_Model_02(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_02_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,7);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,7,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,7,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,7);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,7,6);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,7,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 6; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,7);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2183303508U));
  ssSetChecksum1(S,(3869545555U));
  ssSetChecksum2(S,(3849384848U));
  ssSetChecksum3(S,(620228778U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c7_Model_02(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c7_Model_02(SimStruct *S)
{
  SFc7_Model_02InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc7_Model_02InstanceStruct *)utMalloc(sizeof
    (SFc7_Model_02InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc7_Model_02InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c7_Model_02;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c7_Model_02;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c7_Model_02;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c7_Model_02;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c7_Model_02;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c7_Model_02;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c7_Model_02;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c7_Model_02;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c7_Model_02;
  chartInstance->chartInfo.mdlStart = mdlStart_c7_Model_02;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c7_Model_02;
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

void c7_Model_02_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c7_Model_02(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c7_Model_02(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c7_Model_02(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c7_Model_02_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

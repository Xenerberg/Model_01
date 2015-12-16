/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_01_sfun.h"
#include "c11_Model_01.h"
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
static const char * c11_debug_family_names[10] = { "q_0", "q_v", "Q_v", "R",
  "first", "nargin", "nargout", "q", "rho_t_k", "Matrix_1" };

static const char * c11_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

/* Function Declarations */
static void initialize_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance);
static void initialize_params_c11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance);
static void enable_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance);
static void disable_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance);
static void c11_update_debugger_state_c11_Model_01(SFc11_Model_01InstanceStruct *
  chartInstance);
static const mxArray *get_sim_state_c11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance);
static void set_sim_state_c11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance, const mxArray *c11_st);
static void finalize_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance);
static void sf_gateway_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance);
static void c11_chartstep_c11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance);
static void initSimStructsc11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance);
static void c11_fn_VectorToSkewSymmetricTensor(SFc11_Model_01InstanceStruct
  *chartInstance, real_T c11_v[3], real_T c11_SkewSymmetricTensor[9]);
static void init_script_number_translation(uint32_T c11_machineNumber, uint32_T
  c11_chartNumber, uint32_T c11_instanceNumber);
static const mxArray *c11_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData);
static void c11_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_Matrix_1, const char_T *c11_identifier, real_T c11_y[63]);
static void c11_b_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId, real_T c11_y[63]);
static void c11_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData);
static const mxArray *c11_b_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData);
static const mxArray *c11_c_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData);
static const mxArray *c11_d_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData);
static real_T c11_c_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId);
static void c11_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData);
static const mxArray *c11_e_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData);
static void c11_d_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId, real_T c11_y[9]);
static void c11_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData);
static void c11_e_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId, real_T c11_y[3]);
static void c11_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData);
static void c11_info_helper(const mxArray **c11_info);
static const mxArray *c11_emlrt_marshallOut(const char * c11_u);
static const mxArray *c11_b_emlrt_marshallOut(const uint32_T c11_u);
static void c11_eml_scalar_eg(SFc11_Model_01InstanceStruct *chartInstance);
static void c11_b_eml_scalar_eg(SFc11_Model_01InstanceStruct *chartInstance);
static void c11_eml_xgemm(SFc11_Model_01InstanceStruct *chartInstance, real_T
  c11_A[9], real_T c11_B[9], real_T c11_C[9], real_T c11_b_C[9]);
static const mxArray *c11_f_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData);
static int32_T c11_f_emlrt_marshallIn(SFc11_Model_01InstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId);
static void c11_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData);
static uint8_T c11_g_emlrt_marshallIn(SFc11_Model_01InstanceStruct
  *chartInstance, const mxArray *c11_b_is_active_c11_Model_01, const char_T
  *c11_identifier);
static uint8_T c11_h_emlrt_marshallIn(SFc11_Model_01InstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId);
static void c11_b_eml_xgemm(SFc11_Model_01InstanceStruct *chartInstance, real_T
  c11_A[9], real_T c11_B[9], real_T c11_C[9]);
static void init_dsm_address_info(SFc11_Model_01InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance)
{
  chartInstance->c11_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c11_is_active_c11_Model_01 = 0U;
}

static void initialize_params_c11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c11_update_debugger_state_c11_Model_01(SFc11_Model_01InstanceStruct *
  chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance)
{
  const mxArray *c11_st;
  const mxArray *c11_y = NULL;
  int32_T c11_i0;
  real_T c11_u[63];
  const mxArray *c11_b_y = NULL;
  uint8_T c11_hoistedGlobal;
  uint8_T c11_b_u;
  const mxArray *c11_c_y = NULL;
  real_T (*c11_Matrix_1)[63];
  c11_Matrix_1 = (real_T (*)[63])ssGetOutputPortSignal(chartInstance->S, 1);
  c11_st = NULL;
  c11_st = NULL;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_createcellmatrix(2, 1), false);
  for (c11_i0 = 0; c11_i0 < 63; c11_i0++) {
    c11_u[c11_i0] = (*c11_Matrix_1)[c11_i0];
  }

  c11_b_y = NULL;
  sf_mex_assign(&c11_b_y, sf_mex_create("y", c11_u, 0, 0U, 1U, 0U, 2, 3, 21),
                false);
  sf_mex_setcell(c11_y, 0, c11_b_y);
  c11_hoistedGlobal = chartInstance->c11_is_active_c11_Model_01;
  c11_b_u = c11_hoistedGlobal;
  c11_c_y = NULL;
  sf_mex_assign(&c11_c_y, sf_mex_create("y", &c11_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c11_y, 1, c11_c_y);
  sf_mex_assign(&c11_st, c11_y, false);
  return c11_st;
}

static void set_sim_state_c11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance, const mxArray *c11_st)
{
  const mxArray *c11_u;
  real_T c11_dv0[63];
  int32_T c11_i1;
  real_T (*c11_Matrix_1)[63];
  c11_Matrix_1 = (real_T (*)[63])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c11_doneDoubleBufferReInit = true;
  c11_u = sf_mex_dup(c11_st);
  c11_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c11_u, 0)),
                       "Matrix_1", c11_dv0);
  for (c11_i1 = 0; c11_i1 < 63; c11_i1++) {
    (*c11_Matrix_1)[c11_i1] = c11_dv0[c11_i1];
  }

  chartInstance->c11_is_active_c11_Model_01 = c11_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c11_u, 1)),
     "is_active_c11_Model_01");
  sf_mex_destroy(&c11_u);
  c11_update_debugger_state_c11_Model_01(chartInstance);
  sf_mex_destroy(&c11_st);
}

static void finalize_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c11_Model_01(SFc11_Model_01InstanceStruct *chartInstance)
{
  int32_T c11_i2;
  int32_T c11_i3;
  int32_T c11_i4;
  real_T (*c11_rho_t_k)[3];
  real_T (*c11_Matrix_1)[63];
  real_T (*c11_q)[4];
  c11_rho_t_k = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c11_Matrix_1 = (real_T (*)[63])ssGetOutputPortSignal(chartInstance->S, 1);
  c11_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 10U, chartInstance->c11_sfEvent);
  for (c11_i2 = 0; c11_i2 < 4; c11_i2++) {
    _SFD_DATA_RANGE_CHECK((*c11_q)[c11_i2], 0U);
  }

  chartInstance->c11_sfEvent = CALL_EVENT;
  c11_chartstep_c11_Model_01(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_01MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c11_i3 = 0; c11_i3 < 63; c11_i3++) {
    _SFD_DATA_RANGE_CHECK((*c11_Matrix_1)[c11_i3], 1U);
  }

  for (c11_i4 = 0; c11_i4 < 3; c11_i4++) {
    _SFD_DATA_RANGE_CHECK((*c11_rho_t_k)[c11_i4], 2U);
  }
}

static void c11_chartstep_c11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance)
{
  int32_T c11_i5;
  real_T c11_q[4];
  int32_T c11_i6;
  real_T c11_rho_t_k[3];
  uint32_T c11_debug_family_var_map[10];
  real_T c11_q_0;
  real_T c11_q_v[3];
  real_T c11_Q_v[9];
  real_T c11_R[9];
  real_T c11_first[9];
  real_T c11_nargin = 2.0;
  real_T c11_nargout = 1.0;
  real_T c11_Matrix_1[63];
  int32_T c11_i7;
  int32_T c11_i8;
  real_T c11_b_q_v[3];
  real_T c11_dv1[9];
  int32_T c11_i9;
  real_T c11_a;
  real_T c11_b_a;
  real_T c11_c_a;
  real_T c11_ak;
  real_T c11_d_a;
  real_T c11_c;
  real_T c11_e_a;
  int32_T c11_i10;
  static real_T c11_b[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c11_y[9];
  real_T c11_f_a;
  int32_T c11_i11;
  real_T c11_b_b[9];
  int32_T c11_i12;
  int32_T c11_i13;
  real_T c11_c_b[3];
  int32_T c11_i14;
  int32_T c11_i15;
  real_T c11_d_b[3];
  int32_T c11_i16;
  int32_T c11_i17;
  int32_T c11_i18;
  real_T c11_b_y[9];
  int32_T c11_i19;
  int32_T c11_i20;
  int32_T c11_i21;
  int32_T c11_i22;
  real_T c11_c_q_v[3];
  int32_T c11_i23;
  int32_T c11_i24;
  int32_T c11_i25;
  real_T c11_dv2[9];
  int32_T c11_i26;
  real_T c11_dv3[9];
  int32_T c11_i27;
  real_T c11_dv4[9];
  int32_T c11_i28;
  real_T c11_dv5[9];
  int32_T c11_i29;
  int32_T c11_k;
  int32_T c11_b_k;
  int32_T c11_i30;
  int32_T c11_i31;
  int32_T c11_i32;
  int32_T c11_i33;
  int32_T c11_i34;
  int32_T c11_i35;
  int32_T c11_i36;
  int32_T c11_i37;
  int32_T c11_i38;
  int32_T c11_i39;
  int32_T c11_i40;
  int32_T c11_i41;
  int32_T c11_i42;
  int32_T c11_i43;
  int32_T c11_i44;
  int32_T c11_i45;
  int32_T c11_i46;
  int32_T c11_i47;
  int32_T c11_i48;
  real_T (*c11_b_Matrix_1)[63];
  real_T (*c11_b_rho_t_k)[3];
  real_T (*c11_b_q)[4];
  c11_b_rho_t_k = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c11_b_Matrix_1 = (real_T (*)[63])ssGetOutputPortSignal(chartInstance->S, 1);
  c11_b_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 10U, chartInstance->c11_sfEvent);
  for (c11_i5 = 0; c11_i5 < 4; c11_i5++) {
    c11_q[c11_i5] = (*c11_b_q)[c11_i5];
  }

  for (c11_i6 = 0; c11_i6 < 3; c11_i6++) {
    c11_rho_t_k[c11_i6] = (*c11_b_rho_t_k)[c11_i6];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 10U, 10U, c11_debug_family_names,
    c11_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c11_q_0, 0U, c11_d_sf_marshallOut,
    c11_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_q_v, 1U, c11_b_sf_marshallOut,
    c11_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_Q_v, 2U, c11_e_sf_marshallOut,
    c11_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_R, 3U, c11_e_sf_marshallOut,
    c11_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_first, 4U, c11_e_sf_marshallOut,
    c11_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c11_nargin, 5U, c11_d_sf_marshallOut,
    c11_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c11_nargout, 6U, c11_d_sf_marshallOut,
    c11_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c11_q, 7U, c11_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c11_rho_t_k, 8U, c11_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_Matrix_1, 9U, c11_sf_marshallOut,
    c11_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 3);
  c11_q_0 = c11_q[3];
  _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 4);
  for (c11_i7 = 0; c11_i7 < 3; c11_i7++) {
    c11_q_v[c11_i7] = c11_q[c11_i7];
  }

  _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 5);
  for (c11_i8 = 0; c11_i8 < 3; c11_i8++) {
    c11_b_q_v[c11_i8] = c11_q_v[c11_i8];
  }

  c11_fn_VectorToSkewSymmetricTensor(chartInstance, c11_b_q_v, c11_dv1);
  for (c11_i9 = 0; c11_i9 < 9; c11_i9++) {
    c11_Q_v[c11_i9] = c11_dv1[c11_i9];
  }

  _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 6);
  c11_a = c11_q_0;
  c11_b_a = c11_a;
  c11_c_a = c11_b_a;
  c11_eml_scalar_eg(chartInstance);
  c11_ak = c11_c_a;
  c11_d_a = c11_ak;
  c11_eml_scalar_eg(chartInstance);
  c11_c = c11_d_a * c11_d_a;
  c11_e_a = 2.0 * c11_c - 1.0;
  for (c11_i10 = 0; c11_i10 < 9; c11_i10++) {
    c11_y[c11_i10] = c11_e_a * c11_b[c11_i10];
  }

  c11_f_a = 2.0 * c11_q_0;
  for (c11_i11 = 0; c11_i11 < 9; c11_i11++) {
    c11_b_b[c11_i11] = c11_Q_v[c11_i11];
  }

  for (c11_i12 = 0; c11_i12 < 9; c11_i12++) {
    c11_b_b[c11_i12] *= c11_f_a;
  }

  for (c11_i13 = 0; c11_i13 < 3; c11_i13++) {
    c11_c_b[c11_i13] = c11_q_v[c11_i13];
  }

  for (c11_i14 = 0; c11_i14 < 3; c11_i14++) {
    c11_c_b[c11_i14] *= 2.0;
  }

  for (c11_i15 = 0; c11_i15 < 3; c11_i15++) {
    c11_d_b[c11_i15] = c11_q_v[c11_i15];
  }

  for (c11_i16 = 0; c11_i16 < 3; c11_i16++) {
    c11_i17 = 0;
    for (c11_i18 = 0; c11_i18 < 3; c11_i18++) {
      c11_b_y[c11_i17 + c11_i16] = c11_c_b[c11_i16] * c11_d_b[c11_i18];
      c11_i17 += 3;
    }
  }

  for (c11_i19 = 0; c11_i19 < 9; c11_i19++) {
    c11_R[c11_i19] = (c11_y[c11_i19] + c11_b_b[c11_i19]) + c11_b_y[c11_i19];
  }

  _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 8);
  for (c11_i20 = 0; c11_i20 < 9; c11_i20++) {
    c11_b_b[c11_i20] = c11_R[c11_i20];
  }

  for (c11_i21 = 0; c11_i21 < 9; c11_i21++) {
    c11_b_b[c11_i21] *= -2.0;
  }

  for (c11_i22 = 0; c11_i22 < 3; c11_i22++) {
    c11_c_q_v[c11_i22] = c11_q_v[c11_i22];
  }

  c11_fn_VectorToSkewSymmetricTensor(chartInstance, c11_c_q_v, c11_y);
  c11_b_eml_scalar_eg(chartInstance);
  c11_b_eml_scalar_eg(chartInstance);
  for (c11_i23 = 0; c11_i23 < 9; c11_i23++) {
    c11_first[c11_i23] = 0.0;
  }

  for (c11_i24 = 0; c11_i24 < 9; c11_i24++) {
    c11_first[c11_i24] = 0.0;
  }

  for (c11_i25 = 0; c11_i25 < 9; c11_i25++) {
    c11_dv2[c11_i25] = c11_b_b[c11_i25];
  }

  for (c11_i26 = 0; c11_i26 < 9; c11_i26++) {
    c11_dv3[c11_i26] = c11_y[c11_i26];
  }

  for (c11_i27 = 0; c11_i27 < 9; c11_i27++) {
    c11_dv4[c11_i27] = c11_dv2[c11_i27];
  }

  for (c11_i28 = 0; c11_i28 < 9; c11_i28++) {
    c11_dv5[c11_i28] = c11_dv3[c11_i28];
  }

  c11_b_eml_xgemm(chartInstance, c11_dv4, c11_dv5, c11_first);
  _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 9);
  for (c11_i29 = 0; c11_i29 < 9; c11_i29++) {
    c11_y[c11_i29] = 0.0;
  }

  for (c11_k = 1; c11_k < 4; c11_k++) {
    c11_b_k = c11_k;
    c11_y[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
             (real_T)c11_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
             (int32_T)_SFD_INTEGER_CHECK("", (real_T)c11_b_k), 1, 3, 2, 0) - 1))
      - 1] = 1.0;
  }

  c11_i30 = 0;
  for (c11_i31 = 0; c11_i31 < 3; c11_i31++) {
    for (c11_i32 = 0; c11_i32 < 3; c11_i32++) {
      c11_Matrix_1[c11_i32 + c11_i30] = c11_first[c11_i32 + c11_i30];
    }

    c11_i30 += 3;
  }

  c11_i33 = 0;
  for (c11_i34 = 0; c11_i34 < 6; c11_i34++) {
    for (c11_i35 = 0; c11_i35 < 3; c11_i35++) {
      c11_Matrix_1[(c11_i35 + c11_i33) + 9] = 0.0;
    }

    c11_i33 += 3;
  }

  c11_i36 = 0;
  for (c11_i37 = 0; c11_i37 < 3; c11_i37++) {
    for (c11_i38 = 0; c11_i38 < 3; c11_i38++) {
      c11_Matrix_1[(c11_i38 + c11_i36) + 27] = c11_y[c11_i38 + c11_i36];
    }

    c11_i36 += 3;
  }

  c11_i39 = 0;
  for (c11_i40 = 0; c11_i40 < 3; c11_i40++) {
    for (c11_i41 = 0; c11_i41 < 3; c11_i41++) {
      c11_Matrix_1[(c11_i41 + c11_i39) + 36] = 0.0;
    }

    c11_i39 += 3;
  }

  c11_i42 = 0;
  for (c11_i43 = 0; c11_i43 < 3; c11_i43++) {
    for (c11_i44 = 0; c11_i44 < 3; c11_i44++) {
      c11_Matrix_1[(c11_i44 + c11_i42) + 45] = c11_R[c11_i44 + c11_i42];
    }

    c11_i42 += 3;
  }

  c11_i45 = 0;
  for (c11_i46 = 0; c11_i46 < 3; c11_i46++) {
    for (c11_i47 = 0; c11_i47 < 3; c11_i47++) {
      c11_Matrix_1[(c11_i47 + c11_i45) + 54] = 0.0;
    }

    c11_i45 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, -9);
  _SFD_SYMBOL_SCOPE_POP();
  for (c11_i48 = 0; c11_i48 < 63; c11_i48++) {
    (*c11_b_Matrix_1)[c11_i48] = c11_Matrix_1[c11_i48];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 10U, chartInstance->c11_sfEvent);
}

static void initSimStructsc11_Model_01(SFc11_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c11_fn_VectorToSkewSymmetricTensor(SFc11_Model_01InstanceStruct
  *chartInstance, real_T c11_v[3], real_T c11_SkewSymmetricTensor[9])
{
  uint32_T c11_debug_family_var_map[4];
  real_T c11_nargin = 1.0;
  real_T c11_nargout = 1.0;
  int32_T c11_i49;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c11_b_debug_family_names,
    c11_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c11_nargin, 0U, c11_d_sf_marshallOut,
    c11_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c11_nargout, 1U, c11_d_sf_marshallOut,
    c11_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_v, 2U, c11_b_sf_marshallOut,
    c11_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_SkewSymmetricTensor, 3U,
    c11_e_sf_marshallOut, c11_c_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 2);
  for (c11_i49 = 0; c11_i49 < 9; c11_i49++) {
    c11_SkewSymmetricTensor[c11_i49] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 3);
  c11_SkewSymmetricTensor[0] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 4);
  c11_SkewSymmetricTensor[4] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 5);
  c11_SkewSymmetricTensor[8] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 6);
  c11_SkewSymmetricTensor[3] = -c11_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 7);
  c11_SkewSymmetricTensor[6] = c11_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 8);
  c11_SkewSymmetricTensor[7] = -c11_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 9);
  c11_SkewSymmetricTensor[1] = c11_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 10);
  c11_SkewSymmetricTensor[2] = -c11_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, 11);
  c11_SkewSymmetricTensor[5] = c11_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c11_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c11_machineNumber, uint32_T
  c11_chartNumber, uint32_T c11_instanceNumber)
{
  (void)c11_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c11_chartNumber, c11_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg-2\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c11_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData)
{
  const mxArray *c11_mxArrayOutData = NULL;
  int32_T c11_i50;
  int32_T c11_i51;
  int32_T c11_i52;
  real_T c11_b_inData[63];
  int32_T c11_i53;
  int32_T c11_i54;
  int32_T c11_i55;
  real_T c11_u[63];
  const mxArray *c11_y = NULL;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_mxArrayOutData = NULL;
  c11_i50 = 0;
  for (c11_i51 = 0; c11_i51 < 21; c11_i51++) {
    for (c11_i52 = 0; c11_i52 < 3; c11_i52++) {
      c11_b_inData[c11_i52 + c11_i50] = (*(real_T (*)[63])c11_inData)[c11_i52 +
        c11_i50];
    }

    c11_i50 += 3;
  }

  c11_i53 = 0;
  for (c11_i54 = 0; c11_i54 < 21; c11_i54++) {
    for (c11_i55 = 0; c11_i55 < 3; c11_i55++) {
      c11_u[c11_i55 + c11_i53] = c11_b_inData[c11_i55 + c11_i53];
    }

    c11_i53 += 3;
  }

  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", c11_u, 0, 0U, 1U, 0U, 2, 3, 21),
                false);
  sf_mex_assign(&c11_mxArrayOutData, c11_y, false);
  return c11_mxArrayOutData;
}

static void c11_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_Matrix_1, const char_T *c11_identifier, real_T c11_y[63])
{
  emlrtMsgIdentifier c11_thisId;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_Matrix_1), &c11_thisId,
    c11_y);
  sf_mex_destroy(&c11_Matrix_1);
}

static void c11_b_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId, real_T c11_y[63])
{
  real_T c11_dv6[63];
  int32_T c11_i56;
  (void)chartInstance;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), c11_dv6, 1, 0, 0U, 1, 0U, 2, 3,
                21);
  for (c11_i56 = 0; c11_i56 < 63; c11_i56++) {
    c11_y[c11_i56] = c11_dv6[c11_i56];
  }

  sf_mex_destroy(&c11_u);
}

static void c11_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData)
{
  const mxArray *c11_Matrix_1;
  const char_T *c11_identifier;
  emlrtMsgIdentifier c11_thisId;
  real_T c11_y[63];
  int32_T c11_i57;
  int32_T c11_i58;
  int32_T c11_i59;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_Matrix_1 = sf_mex_dup(c11_mxArrayInData);
  c11_identifier = c11_varName;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_Matrix_1), &c11_thisId,
    c11_y);
  sf_mex_destroy(&c11_Matrix_1);
  c11_i57 = 0;
  for (c11_i58 = 0; c11_i58 < 21; c11_i58++) {
    for (c11_i59 = 0; c11_i59 < 3; c11_i59++) {
      (*(real_T (*)[63])c11_outData)[c11_i59 + c11_i57] = c11_y[c11_i59 +
        c11_i57];
    }

    c11_i57 += 3;
  }

  sf_mex_destroy(&c11_mxArrayInData);
}

static const mxArray *c11_b_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData)
{
  const mxArray *c11_mxArrayOutData = NULL;
  int32_T c11_i60;
  real_T c11_b_inData[3];
  int32_T c11_i61;
  real_T c11_u[3];
  const mxArray *c11_y = NULL;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_mxArrayOutData = NULL;
  for (c11_i60 = 0; c11_i60 < 3; c11_i60++) {
    c11_b_inData[c11_i60] = (*(real_T (*)[3])c11_inData)[c11_i60];
  }

  for (c11_i61 = 0; c11_i61 < 3; c11_i61++) {
    c11_u[c11_i61] = c11_b_inData[c11_i61];
  }

  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", c11_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c11_mxArrayOutData, c11_y, false);
  return c11_mxArrayOutData;
}

static const mxArray *c11_c_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData)
{
  const mxArray *c11_mxArrayOutData = NULL;
  int32_T c11_i62;
  real_T c11_b_inData[4];
  int32_T c11_i63;
  real_T c11_u[4];
  const mxArray *c11_y = NULL;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_mxArrayOutData = NULL;
  for (c11_i62 = 0; c11_i62 < 4; c11_i62++) {
    c11_b_inData[c11_i62] = (*(real_T (*)[4])c11_inData)[c11_i62];
  }

  for (c11_i63 = 0; c11_i63 < 4; c11_i63++) {
    c11_u[c11_i63] = c11_b_inData[c11_i63];
  }

  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", c11_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c11_mxArrayOutData, c11_y, false);
  return c11_mxArrayOutData;
}

static const mxArray *c11_d_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData)
{
  const mxArray *c11_mxArrayOutData = NULL;
  real_T c11_u;
  const mxArray *c11_y = NULL;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_mxArrayOutData = NULL;
  c11_u = *(real_T *)c11_inData;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", &c11_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c11_mxArrayOutData, c11_y, false);
  return c11_mxArrayOutData;
}

static real_T c11_c_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId)
{
  real_T c11_y;
  real_T c11_d0;
  (void)chartInstance;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), &c11_d0, 1, 0, 0U, 0, 0U, 0);
  c11_y = c11_d0;
  sf_mex_destroy(&c11_u);
  return c11_y;
}

static void c11_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData)
{
  const mxArray *c11_nargout;
  const char_T *c11_identifier;
  emlrtMsgIdentifier c11_thisId;
  real_T c11_y;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_nargout = sf_mex_dup(c11_mxArrayInData);
  c11_identifier = c11_varName;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_y = c11_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_nargout),
    &c11_thisId);
  sf_mex_destroy(&c11_nargout);
  *(real_T *)c11_outData = c11_y;
  sf_mex_destroy(&c11_mxArrayInData);
}

static const mxArray *c11_e_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData)
{
  const mxArray *c11_mxArrayOutData = NULL;
  int32_T c11_i64;
  int32_T c11_i65;
  int32_T c11_i66;
  real_T c11_b_inData[9];
  int32_T c11_i67;
  int32_T c11_i68;
  int32_T c11_i69;
  real_T c11_u[9];
  const mxArray *c11_y = NULL;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_mxArrayOutData = NULL;
  c11_i64 = 0;
  for (c11_i65 = 0; c11_i65 < 3; c11_i65++) {
    for (c11_i66 = 0; c11_i66 < 3; c11_i66++) {
      c11_b_inData[c11_i66 + c11_i64] = (*(real_T (*)[9])c11_inData)[c11_i66 +
        c11_i64];
    }

    c11_i64 += 3;
  }

  c11_i67 = 0;
  for (c11_i68 = 0; c11_i68 < 3; c11_i68++) {
    for (c11_i69 = 0; c11_i69 < 3; c11_i69++) {
      c11_u[c11_i69 + c11_i67] = c11_b_inData[c11_i69 + c11_i67];
    }

    c11_i67 += 3;
  }

  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", c11_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c11_mxArrayOutData, c11_y, false);
  return c11_mxArrayOutData;
}

static void c11_d_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId, real_T c11_y[9])
{
  real_T c11_dv7[9];
  int32_T c11_i70;
  (void)chartInstance;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), c11_dv7, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c11_i70 = 0; c11_i70 < 9; c11_i70++) {
    c11_y[c11_i70] = c11_dv7[c11_i70];
  }

  sf_mex_destroy(&c11_u);
}

static void c11_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData)
{
  const mxArray *c11_first;
  const char_T *c11_identifier;
  emlrtMsgIdentifier c11_thisId;
  real_T c11_y[9];
  int32_T c11_i71;
  int32_T c11_i72;
  int32_T c11_i73;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_first = sf_mex_dup(c11_mxArrayInData);
  c11_identifier = c11_varName;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_first), &c11_thisId,
    c11_y);
  sf_mex_destroy(&c11_first);
  c11_i71 = 0;
  for (c11_i72 = 0; c11_i72 < 3; c11_i72++) {
    for (c11_i73 = 0; c11_i73 < 3; c11_i73++) {
      (*(real_T (*)[9])c11_outData)[c11_i73 + c11_i71] = c11_y[c11_i73 + c11_i71];
    }

    c11_i71 += 3;
  }

  sf_mex_destroy(&c11_mxArrayInData);
}

static void c11_e_emlrt_marshallIn(SFc11_Model_01InstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId, real_T c11_y[3])
{
  real_T c11_dv8[3];
  int32_T c11_i74;
  (void)chartInstance;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), c11_dv8, 1, 0, 0U, 1, 0U, 1, 3);
  for (c11_i74 = 0; c11_i74 < 3; c11_i74++) {
    c11_y[c11_i74] = c11_dv8[c11_i74];
  }

  sf_mex_destroy(&c11_u);
}

static void c11_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData)
{
  const mxArray *c11_q_v;
  const char_T *c11_identifier;
  emlrtMsgIdentifier c11_thisId;
  real_T c11_y[3];
  int32_T c11_i75;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_q_v = sf_mex_dup(c11_mxArrayInData);
  c11_identifier = c11_varName;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_q_v), &c11_thisId, c11_y);
  sf_mex_destroy(&c11_q_v);
  for (c11_i75 = 0; c11_i75 < 3; c11_i75++) {
    (*(real_T (*)[3])c11_outData)[c11_i75] = c11_y[c11_i75];
  }

  sf_mex_destroy(&c11_mxArrayInData);
}

const mxArray *sf_c11_Model_01_get_eml_resolved_functions_info(void)
{
  const mxArray *c11_nameCaptureInfo = NULL;
  c11_nameCaptureInfo = NULL;
  sf_mex_assign(&c11_nameCaptureInfo, sf_mex_createstruct("structure", 2, 44, 1),
                false);
  c11_info_helper(&c11_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c11_nameCaptureInfo);
  return c11_nameCaptureInfo;
}

static void c11_info_helper(const mxArray **c11_info)
{
  const mxArray *c11_rhs0 = NULL;
  const mxArray *c11_lhs0 = NULL;
  const mxArray *c11_rhs1 = NULL;
  const mxArray *c11_lhs1 = NULL;
  const mxArray *c11_rhs2 = NULL;
  const mxArray *c11_lhs2 = NULL;
  const mxArray *c11_rhs3 = NULL;
  const mxArray *c11_lhs3 = NULL;
  const mxArray *c11_rhs4 = NULL;
  const mxArray *c11_lhs4 = NULL;
  const mxArray *c11_rhs5 = NULL;
  const mxArray *c11_lhs5 = NULL;
  const mxArray *c11_rhs6 = NULL;
  const mxArray *c11_lhs6 = NULL;
  const mxArray *c11_rhs7 = NULL;
  const mxArray *c11_lhs7 = NULL;
  const mxArray *c11_rhs8 = NULL;
  const mxArray *c11_lhs8 = NULL;
  const mxArray *c11_rhs9 = NULL;
  const mxArray *c11_lhs9 = NULL;
  const mxArray *c11_rhs10 = NULL;
  const mxArray *c11_lhs10 = NULL;
  const mxArray *c11_rhs11 = NULL;
  const mxArray *c11_lhs11 = NULL;
  const mxArray *c11_rhs12 = NULL;
  const mxArray *c11_lhs12 = NULL;
  const mxArray *c11_rhs13 = NULL;
  const mxArray *c11_lhs13 = NULL;
  const mxArray *c11_rhs14 = NULL;
  const mxArray *c11_lhs14 = NULL;
  const mxArray *c11_rhs15 = NULL;
  const mxArray *c11_lhs15 = NULL;
  const mxArray *c11_rhs16 = NULL;
  const mxArray *c11_lhs16 = NULL;
  const mxArray *c11_rhs17 = NULL;
  const mxArray *c11_lhs17 = NULL;
  const mxArray *c11_rhs18 = NULL;
  const mxArray *c11_lhs18 = NULL;
  const mxArray *c11_rhs19 = NULL;
  const mxArray *c11_lhs19 = NULL;
  const mxArray *c11_rhs20 = NULL;
  const mxArray *c11_lhs20 = NULL;
  const mxArray *c11_rhs21 = NULL;
  const mxArray *c11_lhs21 = NULL;
  const mxArray *c11_rhs22 = NULL;
  const mxArray *c11_lhs22 = NULL;
  const mxArray *c11_rhs23 = NULL;
  const mxArray *c11_lhs23 = NULL;
  const mxArray *c11_rhs24 = NULL;
  const mxArray *c11_lhs24 = NULL;
  const mxArray *c11_rhs25 = NULL;
  const mxArray *c11_lhs25 = NULL;
  const mxArray *c11_rhs26 = NULL;
  const mxArray *c11_lhs26 = NULL;
  const mxArray *c11_rhs27 = NULL;
  const mxArray *c11_lhs27 = NULL;
  const mxArray *c11_rhs28 = NULL;
  const mxArray *c11_lhs28 = NULL;
  const mxArray *c11_rhs29 = NULL;
  const mxArray *c11_lhs29 = NULL;
  const mxArray *c11_rhs30 = NULL;
  const mxArray *c11_lhs30 = NULL;
  const mxArray *c11_rhs31 = NULL;
  const mxArray *c11_lhs31 = NULL;
  const mxArray *c11_rhs32 = NULL;
  const mxArray *c11_lhs32 = NULL;
  const mxArray *c11_rhs33 = NULL;
  const mxArray *c11_lhs33 = NULL;
  const mxArray *c11_rhs34 = NULL;
  const mxArray *c11_lhs34 = NULL;
  const mxArray *c11_rhs35 = NULL;
  const mxArray *c11_lhs35 = NULL;
  const mxArray *c11_rhs36 = NULL;
  const mxArray *c11_lhs36 = NULL;
  const mxArray *c11_rhs37 = NULL;
  const mxArray *c11_lhs37 = NULL;
  const mxArray *c11_rhs38 = NULL;
  const mxArray *c11_lhs38 = NULL;
  const mxArray *c11_rhs39 = NULL;
  const mxArray *c11_lhs39 = NULL;
  const mxArray *c11_rhs40 = NULL;
  const mxArray *c11_lhs40 = NULL;
  const mxArray *c11_rhs41 = NULL;
  const mxArray *c11_lhs41 = NULL;
  const mxArray *c11_rhs42 = NULL;
  const mxArray *c11_lhs42 = NULL;
  const mxArray *c11_rhs43 = NULL;
  const mxArray *c11_lhs43 = NULL;
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 0);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1450040424U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c11_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "context", "context", 1);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("mpower"), "name", "name", 1);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c11_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 2);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c11_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("ismatrix"), "name", "name",
                  3);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c11_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("power"), "name", "name", 4);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c11_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 5);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c11_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 6);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 6);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c11_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 7);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c11_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 8);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 8);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c11_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 9);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c11_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 10);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("floor"), "name", "name", 10);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c11_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 11);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c11_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 12);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c11_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 13);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 13);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c11_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "context", "context", 14);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eye"), "name", "name", 14);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c11_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 15);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c11_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 16);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c11_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 17);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("isinf"), "name", "name", 17);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c11_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 18);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c11_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 19);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_is_integer_class"),
                  "name", "name", 19);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c11_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 20);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("intmax"), "name", "name", 20);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c11_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 21);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c11_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 22);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("intmin"), "name", "name", 22);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c11_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 23);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c11_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 24);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.indexIntRelop"), "name", "name", 24);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c11_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 25);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 25);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c11_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 26);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 26);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c11_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 27);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("intmin"), "name", "name", 27);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c11_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 28);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 28);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c11_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 29);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("intmax"), "name", "name", 29);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c11_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 30);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 30);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c11_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 31);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("intmax"), "name", "name", 31);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 31);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c11_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "context", "context", 32);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 32);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c11_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 33);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 33);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c11_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 34);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 34);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c11_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 35);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 35);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c11_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 36);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  36);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c11_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 37);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 37);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c11_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 38);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c11_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 39);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 39);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c11_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 40);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.blas.threshold"), "name", "name", 40);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c11_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 41);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 41);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c11_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 42);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 42);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c11_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 43);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.refblas.xgemm"), "name", "name", 43);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c11_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c11_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs43), "lhs", "lhs",
                  43);
  sf_mex_destroy(&c11_rhs0);
  sf_mex_destroy(&c11_lhs0);
  sf_mex_destroy(&c11_rhs1);
  sf_mex_destroy(&c11_lhs1);
  sf_mex_destroy(&c11_rhs2);
  sf_mex_destroy(&c11_lhs2);
  sf_mex_destroy(&c11_rhs3);
  sf_mex_destroy(&c11_lhs3);
  sf_mex_destroy(&c11_rhs4);
  sf_mex_destroy(&c11_lhs4);
  sf_mex_destroy(&c11_rhs5);
  sf_mex_destroy(&c11_lhs5);
  sf_mex_destroy(&c11_rhs6);
  sf_mex_destroy(&c11_lhs6);
  sf_mex_destroy(&c11_rhs7);
  sf_mex_destroy(&c11_lhs7);
  sf_mex_destroy(&c11_rhs8);
  sf_mex_destroy(&c11_lhs8);
  sf_mex_destroy(&c11_rhs9);
  sf_mex_destroy(&c11_lhs9);
  sf_mex_destroy(&c11_rhs10);
  sf_mex_destroy(&c11_lhs10);
  sf_mex_destroy(&c11_rhs11);
  sf_mex_destroy(&c11_lhs11);
  sf_mex_destroy(&c11_rhs12);
  sf_mex_destroy(&c11_lhs12);
  sf_mex_destroy(&c11_rhs13);
  sf_mex_destroy(&c11_lhs13);
  sf_mex_destroy(&c11_rhs14);
  sf_mex_destroy(&c11_lhs14);
  sf_mex_destroy(&c11_rhs15);
  sf_mex_destroy(&c11_lhs15);
  sf_mex_destroy(&c11_rhs16);
  sf_mex_destroy(&c11_lhs16);
  sf_mex_destroy(&c11_rhs17);
  sf_mex_destroy(&c11_lhs17);
  sf_mex_destroy(&c11_rhs18);
  sf_mex_destroy(&c11_lhs18);
  sf_mex_destroy(&c11_rhs19);
  sf_mex_destroy(&c11_lhs19);
  sf_mex_destroy(&c11_rhs20);
  sf_mex_destroy(&c11_lhs20);
  sf_mex_destroy(&c11_rhs21);
  sf_mex_destroy(&c11_lhs21);
  sf_mex_destroy(&c11_rhs22);
  sf_mex_destroy(&c11_lhs22);
  sf_mex_destroy(&c11_rhs23);
  sf_mex_destroy(&c11_lhs23);
  sf_mex_destroy(&c11_rhs24);
  sf_mex_destroy(&c11_lhs24);
  sf_mex_destroy(&c11_rhs25);
  sf_mex_destroy(&c11_lhs25);
  sf_mex_destroy(&c11_rhs26);
  sf_mex_destroy(&c11_lhs26);
  sf_mex_destroy(&c11_rhs27);
  sf_mex_destroy(&c11_lhs27);
  sf_mex_destroy(&c11_rhs28);
  sf_mex_destroy(&c11_lhs28);
  sf_mex_destroy(&c11_rhs29);
  sf_mex_destroy(&c11_lhs29);
  sf_mex_destroy(&c11_rhs30);
  sf_mex_destroy(&c11_lhs30);
  sf_mex_destroy(&c11_rhs31);
  sf_mex_destroy(&c11_lhs31);
  sf_mex_destroy(&c11_rhs32);
  sf_mex_destroy(&c11_lhs32);
  sf_mex_destroy(&c11_rhs33);
  sf_mex_destroy(&c11_lhs33);
  sf_mex_destroy(&c11_rhs34);
  sf_mex_destroy(&c11_lhs34);
  sf_mex_destroy(&c11_rhs35);
  sf_mex_destroy(&c11_lhs35);
  sf_mex_destroy(&c11_rhs36);
  sf_mex_destroy(&c11_lhs36);
  sf_mex_destroy(&c11_rhs37);
  sf_mex_destroy(&c11_lhs37);
  sf_mex_destroy(&c11_rhs38);
  sf_mex_destroy(&c11_lhs38);
  sf_mex_destroy(&c11_rhs39);
  sf_mex_destroy(&c11_lhs39);
  sf_mex_destroy(&c11_rhs40);
  sf_mex_destroy(&c11_lhs40);
  sf_mex_destroy(&c11_rhs41);
  sf_mex_destroy(&c11_lhs41);
  sf_mex_destroy(&c11_rhs42);
  sf_mex_destroy(&c11_lhs42);
  sf_mex_destroy(&c11_rhs43);
  sf_mex_destroy(&c11_lhs43);
}

static const mxArray *c11_emlrt_marshallOut(const char * c11_u)
{
  const mxArray *c11_y = NULL;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", c11_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c11_u)), false);
  return c11_y;
}

static const mxArray *c11_b_emlrt_marshallOut(const uint32_T c11_u)
{
  const mxArray *c11_y = NULL;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", &c11_u, 7, 0U, 0U, 0U, 0), false);
  return c11_y;
}

static void c11_eml_scalar_eg(SFc11_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c11_b_eml_scalar_eg(SFc11_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c11_eml_xgemm(SFc11_Model_01InstanceStruct *chartInstance, real_T
  c11_A[9], real_T c11_B[9], real_T c11_C[9], real_T c11_b_C[9])
{
  int32_T c11_i76;
  int32_T c11_i77;
  real_T c11_b_A[9];
  int32_T c11_i78;
  real_T c11_b_B[9];
  for (c11_i76 = 0; c11_i76 < 9; c11_i76++) {
    c11_b_C[c11_i76] = c11_C[c11_i76];
  }

  for (c11_i77 = 0; c11_i77 < 9; c11_i77++) {
    c11_b_A[c11_i77] = c11_A[c11_i77];
  }

  for (c11_i78 = 0; c11_i78 < 9; c11_i78++) {
    c11_b_B[c11_i78] = c11_B[c11_i78];
  }

  c11_b_eml_xgemm(chartInstance, c11_b_A, c11_b_B, c11_b_C);
}

static const mxArray *c11_f_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData)
{
  const mxArray *c11_mxArrayOutData = NULL;
  int32_T c11_u;
  const mxArray *c11_y = NULL;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_mxArrayOutData = NULL;
  c11_u = *(int32_T *)c11_inData;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", &c11_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c11_mxArrayOutData, c11_y, false);
  return c11_mxArrayOutData;
}

static int32_T c11_f_emlrt_marshallIn(SFc11_Model_01InstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId)
{
  int32_T c11_y;
  int32_T c11_i79;
  (void)chartInstance;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), &c11_i79, 1, 6, 0U, 0, 0U, 0);
  c11_y = c11_i79;
  sf_mex_destroy(&c11_u);
  return c11_y;
}

static void c11_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData)
{
  const mxArray *c11_b_sfEvent;
  const char_T *c11_identifier;
  emlrtMsgIdentifier c11_thisId;
  int32_T c11_y;
  SFc11_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc11_Model_01InstanceStruct *)chartInstanceVoid;
  c11_b_sfEvent = sf_mex_dup(c11_mxArrayInData);
  c11_identifier = c11_varName;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_y = c11_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_b_sfEvent),
    &c11_thisId);
  sf_mex_destroy(&c11_b_sfEvent);
  *(int32_T *)c11_outData = c11_y;
  sf_mex_destroy(&c11_mxArrayInData);
}

static uint8_T c11_g_emlrt_marshallIn(SFc11_Model_01InstanceStruct
  *chartInstance, const mxArray *c11_b_is_active_c11_Model_01, const char_T
  *c11_identifier)
{
  uint8_T c11_y;
  emlrtMsgIdentifier c11_thisId;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_y = c11_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c11_b_is_active_c11_Model_01), &c11_thisId);
  sf_mex_destroy(&c11_b_is_active_c11_Model_01);
  return c11_y;
}

static uint8_T c11_h_emlrt_marshallIn(SFc11_Model_01InstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId)
{
  uint8_T c11_y;
  uint8_T c11_u0;
  (void)chartInstance;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), &c11_u0, 1, 3, 0U, 0, 0U, 0);
  c11_y = c11_u0;
  sf_mex_destroy(&c11_u);
  return c11_y;
}

static void c11_b_eml_xgemm(SFc11_Model_01InstanceStruct *chartInstance, real_T
  c11_A[9], real_T c11_B[9], real_T c11_C[9])
{
  int32_T c11_i80;
  int32_T c11_i81;
  int32_T c11_i82;
  int32_T c11_i83;
  int32_T c11_i84;
  (void)chartInstance;
  for (c11_i80 = 0; c11_i80 < 3; c11_i80++) {
    c11_i81 = 0;
    for (c11_i82 = 0; c11_i82 < 3; c11_i82++) {
      c11_C[c11_i81 + c11_i80] = 0.0;
      c11_i83 = 0;
      for (c11_i84 = 0; c11_i84 < 3; c11_i84++) {
        c11_C[c11_i81 + c11_i80] += c11_A[c11_i83 + c11_i80] * c11_B[c11_i84 +
          c11_i81];
        c11_i83 += 3;
      }

      c11_i81 += 3;
    }
  }
}

static void init_dsm_address_info(SFc11_Model_01InstanceStruct *chartInstance)
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

void sf_c11_Model_01_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(524028496U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2372881526U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(521254949U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(573049994U);
}

mxArray *sf_c11_Model_01_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("IS5kaJJMSAHNXsjGcxQjbB");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

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
      pr[1] = (double)(21);
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

mxArray *sf_c11_Model_01_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c11_Model_01_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c11_Model_01(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"Matrix_1\",},{M[8],M[0],T\"is_active_c11_Model_01\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c11_Model_01_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc11_Model_01InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc11_Model_01InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_01MachineNumber_,
           11,
           1,
           1,
           0,
           3,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"q");
          _SFD_SET_DATA_PROPS(1,2,0,1,"Matrix_1");
          _SFD_SET_DATA_PROPS(2,1,1,0,"rho_t_k");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,327);
        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"fn_VectorToSkewSymmetricTensor",0,-1,433);

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c11_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 21;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c11_sf_marshallOut,(MexInFcnForType)
            c11_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c11_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T (*c11_q)[4];
          real_T (*c11_Matrix_1)[63];
          real_T (*c11_rho_t_k)[3];
          c11_rho_t_k = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
          c11_Matrix_1 = (real_T (*)[63])ssGetOutputPortSignal(chartInstance->S,
            1);
          c11_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c11_q);
          _SFD_SET_DATA_VALUE_PTR(1U, *c11_Matrix_1);
          _SFD_SET_DATA_VALUE_PTR(2U, *c11_rho_t_k);
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
  return "QwETwijcyK4H4WrgvmSxxG";
}

static void sf_opaque_initialize_c11_Model_01(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc11_Model_01InstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c11_Model_01((SFc11_Model_01InstanceStruct*)
    chartInstanceVar);
  initialize_c11_Model_01((SFc11_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c11_Model_01(void *chartInstanceVar)
{
  enable_c11_Model_01((SFc11_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c11_Model_01(void *chartInstanceVar)
{
  disable_c11_Model_01((SFc11_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c11_Model_01(void *chartInstanceVar)
{
  sf_gateway_c11_Model_01((SFc11_Model_01InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c11_Model_01(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c11_Model_01((SFc11_Model_01InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c11_Model_01();/* state var info */
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

extern void sf_internal_set_sim_state_c11_Model_01(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c11_Model_01();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c11_Model_01((SFc11_Model_01InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c11_Model_01(SimStruct* S)
{
  return sf_internal_get_sim_state_c11_Model_01(S);
}

static void sf_opaque_set_sim_state_c11_Model_01(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c11_Model_01(S, st);
}

static void sf_opaque_terminate_c11_Model_01(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc11_Model_01InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_01_optimization_info();
    }

    finalize_c11_Model_01((SFc11_Model_01InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc11_Model_01((SFc11_Model_01InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c11_Model_01(SimStruct *S)
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
    initialize_params_c11_Model_01((SFc11_Model_01InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c11_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_01_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      11);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,11,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,11,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,11);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,11,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,11,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,11);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3410834018U));
  ssSetChecksum1(S,(1385172485U));
  ssSetChecksum2(S,(4080443286U));
  ssSetChecksum3(S,(3777171823U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c11_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c11_Model_01(SimStruct *S)
{
  SFc11_Model_01InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc11_Model_01InstanceStruct *)utMalloc(sizeof
    (SFc11_Model_01InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc11_Model_01InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c11_Model_01;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c11_Model_01;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c11_Model_01;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c11_Model_01;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c11_Model_01;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c11_Model_01;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c11_Model_01;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c11_Model_01;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c11_Model_01;
  chartInstance->chartInfo.mdlStart = mdlStart_c11_Model_01;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c11_Model_01;
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

void c11_Model_01_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c11_Model_01(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c11_Model_01(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c11_Model_01(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c11_Model_01_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

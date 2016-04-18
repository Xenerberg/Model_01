/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_02_sfun.h"
#include "c6_Model_02.h"
#include <math.h>
#include "mwmathutil.h"
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
static const char * c6_debug_family_names[10] = { "Phi_t", "Phi_r", "nargin",
  "nargout", "X_a", "n", "t_Kalman", "p", "M", "Phi" };

static const char * c6_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

static const char * c6_c_debug_family_names[7] = { "Phi_t12", "Phi_t22",
  "nargin", "nargout", "n", "t_delta", "Phi_t" };

static const char * c6_d_debug_family_names[5] = { "nargin", "nargout", "p",
  "omega", "M" };

static const char * c6_e_debug_family_names[9] = { "omega", "M", "A", "Phi_r",
  "nargin", "nargout", "X_a", "p", "t_Kalman" };

/* Function Declarations */
static void initialize_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance);
static void initialize_params_c6_Model_02(SFc6_Model_02InstanceStruct
  *chartInstance);
static void enable_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance);
static void disable_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance);
static void c6_update_debugger_state_c6_Model_02(SFc6_Model_02InstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c6_Model_02(SFc6_Model_02InstanceStruct
  *chartInstance);
static void set_sim_state_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_st);
static void finalize_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance);
static void sf_gateway_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance);
static void initSimStructsc6_Model_02(SFc6_Model_02InstanceStruct *chartInstance);
static void c6_fn_VectorToSkewSymmetricTensor(SFc6_Model_02InstanceStruct
  *chartInstance, real_T c6_v[3], real_T c6_SkewSymmetricTensor[9]);
static void init_script_number_translation(uint32_T c6_machineNumber, uint32_T
  c6_chartNumber, uint32_T c6_instanceNumber);
static const mxArray *c6_sf_marshallOut(void *chartInstanceVoid, void *c6_inData);
static void c6_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_Phi, const char_T *c6_identifier, real_T c6_y[144]);
static void c6_b_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[144]);
static void c6_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static const mxArray *c6_b_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static const mxArray *c6_c_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static const mxArray *c6_d_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static const mxArray *c6_e_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static real_T c6_c_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId);
static void c6_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static const mxArray *c6_f_sf_marshallOut(void *chartInstanceVoid, real_T
  c6_inData_data[], int32_T c6_inData_sizes[2]);
static void c6_d_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y_data[],
  int32_T c6_y_sizes[2]);
static void c6_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, real_T c6_outData_data[], int32_T
  c6_outData_sizes[2]);
static const mxArray *c6_g_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static void c6_e_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[36]);
static void c6_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static void c6_f_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[9]);
static void c6_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static void c6_g_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[3]);
static void c6_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static void c6_h_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[13]);
static void c6_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static void c6_info_helper(const mxArray **c6_info);
static const mxArray *c6_emlrt_marshallOut(const char * c6_u);
static const mxArray *c6_b_emlrt_marshallOut(const uint32_T c6_u);
static void c6_b_info_helper(const mxArray **c6_info);
static void c6_c_info_helper(const mxArray **c6_info);
static void c6_eml_scalar_eg(SFc6_Model_02InstanceStruct *chartInstance);
static void c6_eml_switch_helper(SFc6_Model_02InstanceStruct *chartInstance);
static void c6_eye(SFc6_Model_02InstanceStruct *chartInstance, real_T c6_I[9]);
static void c6_expm(SFc6_Model_02InstanceStruct *chartInstance, real_T c6_A[36],
                    real_T c6_F[36]);
static void c6_PadeApproximantOfDegree(SFc6_Model_02InstanceStruct
  *chartInstance, real_T c6_A[36], real_T c6_m, real_T c6_F[36]);
static void c6_b_eml_scalar_eg(SFc6_Model_02InstanceStruct *chartInstance);
static void c6_eml_xgemm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36], real_T c6_C[36], real_T c6_b_C[36]);
static void c6_eml_matlab_zgetrf(SFc6_Model_02InstanceStruct *chartInstance,
  real_T c6_A[36], real_T c6_b_A[36], int32_T c6_ipiv[6], int32_T *c6_info);
static int32_T c6_eml_ixamax(SFc6_Model_02InstanceStruct *chartInstance, int32_T
  c6_n, real_T c6_x[36], int32_T c6_ix0);
static void c6_check_forloop_overflow_error(SFc6_Model_02InstanceStruct
  *chartInstance, boolean_T c6_overflow);
static void c6_threshold(SFc6_Model_02InstanceStruct *chartInstance);
static void c6_eml_xgeru(SFc6_Model_02InstanceStruct *chartInstance, int32_T
  c6_m, int32_T c6_n, real_T c6_alpha1, int32_T c6_ix0, int32_T c6_iy0, real_T
  c6_A[36], int32_T c6_ia0, real_T c6_b_A[36]);
static void c6_eml_warning(SFc6_Model_02InstanceStruct *chartInstance);
static void c6_eml_xtrsm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36], real_T c6_b_B[36]);
static void c6_b_threshold(SFc6_Model_02InstanceStruct *chartInstance);
static void c6_scalarEg(SFc6_Model_02InstanceStruct *chartInstance);
static void c6_b_eml_xtrsm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36], real_T c6_b_B[36]);
static const mxArray *c6_h_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static int32_T c6_i_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId);
static void c6_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static uint8_T c6_j_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_b_is_active_c6_Model_02, const char_T *c6_identifier);
static uint8_T c6_k_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId);
static void c6_b_eml_xgemm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36], real_T c6_C[36]);
static void c6_b_eml_matlab_zgetrf(SFc6_Model_02InstanceStruct *chartInstance,
  real_T c6_A[36], int32_T c6_ipiv[6], int32_T *c6_info);
static void c6_b_eml_xgeru(SFc6_Model_02InstanceStruct *chartInstance, int32_T
  c6_m, int32_T c6_n, real_T c6_alpha1, int32_T c6_ix0, int32_T c6_iy0, real_T
  c6_A[36], int32_T c6_ia0);
static void c6_c_eml_xtrsm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36]);
static void c6_d_eml_xtrsm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36]);
static void init_dsm_address_info(SFc6_Model_02InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance)
{
  chartInstance->c6_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c6_is_active_c6_Model_02 = 0U;
}

static void initialize_params_c6_Model_02(SFc6_Model_02InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c6_update_debugger_state_c6_Model_02(SFc6_Model_02InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c6_Model_02(SFc6_Model_02InstanceStruct
  *chartInstance)
{
  const mxArray *c6_st;
  const mxArray *c6_y = NULL;
  int32_T c6_i0;
  real_T c6_u[144];
  const mxArray *c6_b_y = NULL;
  uint8_T c6_hoistedGlobal;
  uint8_T c6_b_u;
  const mxArray *c6_c_y = NULL;
  real_T (*c6_Phi)[144];
  c6_Phi = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 1);
  c6_st = NULL;
  c6_st = NULL;
  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_createcellmatrix(2, 1), false);
  for (c6_i0 = 0; c6_i0 < 144; c6_i0++) {
    c6_u[c6_i0] = (*c6_Phi)[c6_i0];
  }

  c6_b_y = NULL;
  sf_mex_assign(&c6_b_y, sf_mex_create("y", c6_u, 0, 0U, 1U, 0U, 2, 12, 12),
                false);
  sf_mex_setcell(c6_y, 0, c6_b_y);
  c6_hoistedGlobal = chartInstance->c6_is_active_c6_Model_02;
  c6_b_u = c6_hoistedGlobal;
  c6_c_y = NULL;
  sf_mex_assign(&c6_c_y, sf_mex_create("y", &c6_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c6_y, 1, c6_c_y);
  sf_mex_assign(&c6_st, c6_y, false);
  return c6_st;
}

static void set_sim_state_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_st)
{
  const mxArray *c6_u;
  real_T c6_dv0[144];
  int32_T c6_i1;
  real_T (*c6_Phi)[144];
  c6_Phi = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c6_doneDoubleBufferReInit = true;
  c6_u = sf_mex_dup(c6_st);
  c6_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c6_u, 0)), "Phi",
                      c6_dv0);
  for (c6_i1 = 0; c6_i1 < 144; c6_i1++) {
    (*c6_Phi)[c6_i1] = c6_dv0[c6_i1];
  }

  chartInstance->c6_is_active_c6_Model_02 = c6_j_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c6_u, 1)), "is_active_c6_Model_02");
  sf_mex_destroy(&c6_u);
  c6_update_debugger_state_c6_Model_02(chartInstance);
  sf_mex_destroy(&c6_st);
}

static void finalize_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c6_Model_02(SFc6_Model_02InstanceStruct *chartInstance)
{
  int32_T c6_i2;
  real_T c6_hoistedGlobal;
  real_T c6_b_hoistedGlobal;
  int32_T c6_i3;
  real_T c6_X_a[13];
  real_T c6_n;
  real_T c6_t_Kalman;
  int32_T c6_i4;
  real_T c6_p[3];
  int32_T c6_i5;
  real_T c6_M[9];
  uint32_T c6_debug_family_var_map[10];
  real_T c6_Phi_t[36];
  int32_T c6_Phi_r_sizes[2];
  real_T c6_Phi_r_data[36];
  real_T c6_nargin = 5.0;
  real_T c6_nargout = 1.0;
  real_T c6_Phi[144];
  real_T c6_b_n;
  real_T c6_t_delta;
  uint32_T c6_b_debug_family_var_map[7];
  real_T c6_Phi_t12[9];
  real_T c6_Phi_t22[9];
  real_T c6_b_nargin = 2.0;
  real_T c6_b_nargout = 1.0;
  real_T c6_a;
  real_T c6_b_a;
  real_T c6_c_a;
  real_T c6_ak;
  real_T c6_d_a;
  real_T c6_c;
  real_T c6_dv1[3];
  real_T c6_b[9];
  int32_T c6_i6;
  real_T c6_b_b;
  int32_T c6_i7;
  real_T c6_dv2[9];
  int32_T c6_i8;
  int32_T c6_i9;
  int32_T c6_i10;
  int32_T c6_i11;
  int32_T c6_i12;
  int32_T c6_i13;
  int32_T c6_i14;
  int32_T c6_i15;
  int32_T c6_i16;
  int32_T c6_i17;
  int32_T c6_i18;
  int32_T c6_i19;
  int32_T c6_i20;
  int32_T c6_i21;
  int32_T c6_i22;
  int32_T c6_i23;
  int32_T c6_i24;
  real_T c6_b_X_a[13];
  int32_T c6_i25;
  real_T c6_b_p[3];
  real_T c6_b_t_Kalman;
  uint32_T c6_c_debug_family_var_map[9];
  real_T c6_omega[3];
  real_T c6_b_M[9];
  real_T c6_A[36];
  real_T c6_Phi_r[144];
  real_T c6_c_nargin = 3.0;
  real_T c6_c_nargout = 1.0;
  int32_T c6_i26;
  int32_T c6_i27;
  int32_T c6_i28;
  real_T c6_c_p[3];
  int32_T c6_i29;
  real_T c6_b_omega[3];
  uint32_T c6_d_debug_family_var_map[5];
  real_T c6_d_nargin = 2.0;
  real_T c6_d_nargout = 1.0;
  int32_T c6_i30;
  real_T c6_c_omega[3];
  int32_T c6_i31;
  int32_T c6_i32;
  int32_T c6_i33;
  int32_T c6_i34;
  int32_T c6_i35;
  int32_T c6_i36;
  int32_T c6_i37;
  int32_T c6_i38;
  static real_T c6_y[9] = { 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5 };

  int32_T c6_i39;
  int32_T c6_i40;
  int32_T c6_i41;
  int32_T c6_i42;
  int32_T c6_i43;
  int32_T c6_i44;
  int32_T c6_i45;
  int32_T c6_i46;
  real_T c6_e_a[36];
  real_T c6_c_b;
  int32_T c6_i47;
  int32_T c6_i48;
  real_T c6_f_a[36];
  int32_T c6_b_Phi_r;
  int32_T c6_c_Phi_r;
  int32_T c6_i49;
  int32_T c6_tmp_sizes[2];
  int32_T c6_i50;
  int32_T c6_i51;
  real_T c6_tmp_data[72];
  int32_T c6_i52;
  int32_T c6_i53;
  int32_T c6_i54;
  int32_T c6_i55;
  int32_T c6_i56;
  int32_T c6_i57;
  int32_T c6_i58;
  int32_T c6_i59;
  int32_T c6_i60;
  int32_T c6_i61;
  int32_T c6_i62;
  int32_T c6_i63;
  int32_T c6_i64;
  int32_T c6_i65;
  int32_T c6_i66;
  real_T *c6_c_n;
  real_T *c6_c_t_Kalman;
  real_T (*c6_b_Phi)[144];
  real_T (*c6_c_M)[9];
  real_T (*c6_d_p)[3];
  real_T (*c6_c_X_a)[13];
  c6_c_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 4);
  c6_d_p = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 3);
  c6_c_t_Kalman = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c6_c_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c6_b_Phi = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 1);
  c6_c_X_a = (real_T (*)[13])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 5U, chartInstance->c6_sfEvent);
  for (c6_i2 = 0; c6_i2 < 13; c6_i2++) {
    _SFD_DATA_RANGE_CHECK((*c6_c_X_a)[c6_i2], 0U);
  }

  chartInstance->c6_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 5U, chartInstance->c6_sfEvent);
  c6_hoistedGlobal = *c6_c_n;
  c6_b_hoistedGlobal = *c6_c_t_Kalman;
  for (c6_i3 = 0; c6_i3 < 13; c6_i3++) {
    c6_X_a[c6_i3] = (*c6_c_X_a)[c6_i3];
  }

  c6_n = c6_hoistedGlobal;
  c6_t_Kalman = c6_b_hoistedGlobal;
  for (c6_i4 = 0; c6_i4 < 3; c6_i4++) {
    c6_p[c6_i4] = (*c6_d_p)[c6_i4];
  }

  for (c6_i5 = 0; c6_i5 < 9; c6_i5++) {
    c6_M[c6_i5] = (*c6_c_M)[c6_i5];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 10U, 10U, c6_debug_family_names,
    c6_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_Phi_t, 0U, c6_g_sf_marshallOut,
    c6_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_IMPORTABLE(c6_Phi_r_data, (const int32_T *)
    &c6_Phi_r_sizes, NULL, 0, 1, (void *)c6_f_sf_marshallOut, (void *)
    c6_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_nargin, 2U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_nargout, 3U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c6_X_a, 4U, c6_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c6_n, 5U, c6_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c6_t_Kalman, 6U, c6_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c6_p, 7U, c6_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c6_M, 8U, c6_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_Phi, 9U, c6_sf_marshallOut,
    c6_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 9);
  c6_b_n = c6_n;
  c6_t_delta = c6_t_Kalman;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 7U, 7U, c6_c_debug_family_names,
    c6_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_Phi_t12, 0U, c6_b_sf_marshallOut,
    c6_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_Phi_t22, 1U, c6_b_sf_marshallOut,
    c6_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_b_nargin, 2U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_b_nargout, 3U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_b_n, 4U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_t_delta, 5U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_Phi_t, 6U, c6_g_sf_marshallOut,
    c6_d_sf_marshallIn);
  CV_EML_FCN(0, 4);
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 70);
  c6_a = c6_t_delta;
  c6_b_a = c6_a;
  c6_c_a = c6_b_a;
  c6_eml_scalar_eg(chartInstance);
  c6_ak = c6_c_a;
  c6_d_a = c6_ak;
  c6_eml_scalar_eg(chartInstance);
  c6_c = c6_d_a * c6_d_a;
  c6_Phi_t12[0] = c6_t_delta;
  c6_Phi_t12[3] = c6_b_n * c6_c;
  c6_Phi_t12[6] = 0.0;
  c6_Phi_t12[1] = 0.0;
  c6_Phi_t12[4] = c6_t_delta;
  c6_Phi_t12[7] = 0.0;
  c6_Phi_t12[2] = 0.0;
  c6_Phi_t12[5] = 0.0;
  c6_Phi_t12[8] = c6_t_delta;
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 71);
  c6_dv1[0] = 0.0;
  c6_dv1[1] = 0.0;
  c6_dv1[2] = c6_b_n;
  c6_fn_VectorToSkewSymmetricTensor(chartInstance, c6_dv1, c6_b);
  for (c6_i6 = 0; c6_i6 < 9; c6_i6++) {
    c6_b[c6_i6] *= 2.0;
  }

  c6_b_b = c6_t_delta;
  for (c6_i7 = 0; c6_i7 < 9; c6_i7++) {
    c6_b[c6_i7] *= c6_b_b;
  }

  c6_eye(chartInstance, c6_dv2);
  for (c6_i8 = 0; c6_i8 < 9; c6_i8++) {
    c6_Phi_t22[c6_i8] = c6_dv2[c6_i8] - c6_b[c6_i8];
  }

  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 73);
  c6_eye(chartInstance, c6_dv2);
  c6_i9 = 0;
  c6_i10 = 0;
  for (c6_i11 = 0; c6_i11 < 3; c6_i11++) {
    for (c6_i12 = 0; c6_i12 < 3; c6_i12++) {
      c6_Phi_t[c6_i12 + c6_i9] = c6_dv2[c6_i12 + c6_i10];
    }

    c6_i9 += 6;
    c6_i10 += 3;
  }

  c6_i13 = 0;
  c6_i14 = 0;
  for (c6_i15 = 0; c6_i15 < 3; c6_i15++) {
    for (c6_i16 = 0; c6_i16 < 3; c6_i16++) {
      c6_Phi_t[(c6_i16 + c6_i13) + 18] = c6_Phi_t12[c6_i16 + c6_i14];
    }

    c6_i13 += 6;
    c6_i14 += 3;
  }

  c6_i17 = 0;
  for (c6_i18 = 0; c6_i18 < 3; c6_i18++) {
    for (c6_i19 = 0; c6_i19 < 3; c6_i19++) {
      c6_Phi_t[(c6_i19 + c6_i17) + 3] = 0.0;
    }

    c6_i17 += 6;
  }

  c6_i20 = 0;
  c6_i21 = 0;
  for (c6_i22 = 0; c6_i22 < 3; c6_i22++) {
    for (c6_i23 = 0; c6_i23 < 3; c6_i23++) {
      c6_Phi_t[(c6_i23 + c6_i20) + 21] = c6_Phi_t22[c6_i23 + c6_i21];
    }

    c6_i20 += 6;
    c6_i21 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, -73);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 10);
  for (c6_i24 = 0; c6_i24 < 13; c6_i24++) {
    c6_b_X_a[c6_i24] = c6_X_a[c6_i24];
  }

  for (c6_i25 = 0; c6_i25 < 3; c6_i25++) {
    c6_b_p[c6_i25] = c6_p[c6_i25];
  }

  c6_b_t_Kalman = c6_t_Kalman;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 9U, 10U, c6_e_debug_family_names,
    c6_c_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_omega, 0U, c6_c_sf_marshallOut,
    c6_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_b_M, 1U, c6_b_sf_marshallOut,
    c6_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_A, 2U, c6_g_sf_marshallOut,
    c6_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_Phi_r, MAX_uint32_T, c6_sf_marshallOut,
    c6_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_c_nargin, 4U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_c_nargout, 5U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_b_X_a, 6U, c6_e_sf_marshallOut,
    c6_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_b_p, 7U, c6_c_sf_marshallOut,
    c6_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_b_t_Kalman, 8U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_IMPORTABLE(c6_Phi_r_data, (const int32_T *)
    &c6_Phi_r_sizes, NULL, 0, -1, (void *)c6_f_sf_marshallOut, (void *)
    c6_c_sf_marshallIn);
  CV_EML_FCN(0, 1);
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 23);
  for (c6_i26 = 0; c6_i26 < 144; c6_i26++) {
    c6_Phi_r[c6_i26] = 0.0;
  }

  _SFD_SYMBOL_SWITCH(3U, 3U);
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 24);
  for (c6_i27 = 0; c6_i27 < 3; c6_i27++) {
    c6_omega[c6_i27] = c6_b_X_a[c6_i27 + 4];
  }

  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 25);
  for (c6_i28 = 0; c6_i28 < 3; c6_i28++) {
    c6_c_p[c6_i28] = c6_b_p[c6_i28];
  }

  for (c6_i29 = 0; c6_i29 < 3; c6_i29++) {
    c6_b_omega[c6_i29] = c6_omega[c6_i29];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 5U, 5U, c6_d_debug_family_names,
    c6_d_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_d_nargin, 0U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_d_nargout, 1U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_c_p, 2U, c6_c_sf_marshallOut,
    c6_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_b_omega, 3U, c6_c_sf_marshallOut,
    c6_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_b_M, 4U, c6_b_sf_marshallOut,
    c6_e_sf_marshallIn);
  CV_EML_FCN(0, 2);
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 48);
  c6_b_M[0] = 0.0;
  c6_b_M[3] = c6_c_p[0] * c6_b_omega[2];
  c6_b_M[6] = c6_c_p[0] * c6_b_omega[1];
  c6_b_M[1] = c6_c_p[1] * c6_b_omega[2];
  c6_b_M[4] = 0.0;
  c6_b_M[7] = c6_c_p[1] * c6_b_omega[0];
  c6_b_M[2] = c6_c_p[2] * c6_b_omega[1];
  c6_b_M[5] = c6_c_p[2] * c6_b_omega[0];
  c6_b_M[8] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, -48);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 27);
  for (c6_i30 = 0; c6_i30 < 3; c6_i30++) {
    c6_c_omega[c6_i30] = c6_omega[c6_i30];
  }

  c6_fn_VectorToSkewSymmetricTensor(chartInstance, c6_c_omega, c6_dv2);
  c6_i31 = 0;
  c6_i32 = 0;
  for (c6_i33 = 0; c6_i33 < 3; c6_i33++) {
    for (c6_i34 = 0; c6_i34 < 3; c6_i34++) {
      c6_A[c6_i34 + c6_i31] = -c6_dv2[c6_i34 + c6_i32];
    }

    c6_i31 += 6;
    c6_i32 += 3;
  }

  c6_i35 = 0;
  c6_i36 = 0;
  for (c6_i37 = 0; c6_i37 < 3; c6_i37++) {
    for (c6_i38 = 0; c6_i38 < 3; c6_i38++) {
      c6_A[(c6_i38 + c6_i35) + 18] = c6_y[c6_i38 + c6_i36];
    }

    c6_i35 += 6;
    c6_i36 += 3;
  }

  c6_i39 = 0;
  for (c6_i40 = 0; c6_i40 < 3; c6_i40++) {
    for (c6_i41 = 0; c6_i41 < 3; c6_i41++) {
      c6_A[(c6_i41 + c6_i39) + 3] = 0.0;
    }

    c6_i39 += 6;
  }

  c6_i42 = 0;
  c6_i43 = 0;
  for (c6_i44 = 0; c6_i44 < 3; c6_i44++) {
    for (c6_i45 = 0; c6_i45 < 3; c6_i45++) {
      c6_A[(c6_i45 + c6_i42) + 21] = c6_b_M[c6_i45 + c6_i43];
    }

    c6_i42 += 6;
    c6_i43 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 28);
  for (c6_i46 = 0; c6_i46 < 36; c6_i46++) {
    c6_e_a[c6_i46] = c6_A[c6_i46];
  }

  c6_c_b = c6_b_t_Kalman;
  for (c6_i47 = 0; c6_i47 < 36; c6_i47++) {
    c6_e_a[c6_i47] *= c6_c_b;
  }

  for (c6_i48 = 0; c6_i48 < 36; c6_i48++) {
    c6_f_a[c6_i48] = c6_e_a[c6_i48];
  }

  c6_expm(chartInstance, c6_f_a, c6_e_a);
  c6_Phi_r_sizes[0] = 6;
  c6_Phi_r_sizes[1] = 6;
  c6_b_Phi_r = c6_Phi_r_sizes[0];
  c6_c_Phi_r = c6_Phi_r_sizes[1];
  for (c6_i49 = 0; c6_i49 < 36; c6_i49++) {
    c6_Phi_r_data[c6_i49] = c6_e_a[c6_i49];
  }

  _SFD_SYMBOL_SWITCH(3U, 9U);
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, -28);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 11);
  c6_tmp_sizes[0] = 6;
  c6_tmp_sizes[1] = 12;
  for (c6_i50 = 0; c6_i50 < 6; c6_i50++) {
    for (c6_i51 = 0; c6_i51 < 6; c6_i51++) {
      c6_tmp_data[c6_i51 + c6_tmp_sizes[0] * c6_i50] = c6_Phi_r_data[c6_i51 +
        c6_Phi_r_sizes[0] * c6_i50];
    }
  }

  for (c6_i52 = 0; c6_i52 < 6; c6_i52++) {
    for (c6_i53 = 0; c6_i53 < 6; c6_i53++) {
      c6_tmp_data[c6_i53 + c6_tmp_sizes[0] * (c6_i52 + 6)] = 0.0;
    }
  }

  for (c6_i54 = 0; c6_i54 < 12; c6_i54++) {
    for (c6_i55 = 0; c6_i55 < 6; c6_i55++) {
      c6_Phi[c6_i55 + 12 * c6_i54] = c6_tmp_data[c6_i55 + c6_tmp_sizes[0] *
        c6_i54];
    }
  }

  c6_i56 = 0;
  for (c6_i57 = 0; c6_i57 < 6; c6_i57++) {
    for (c6_i58 = 0; c6_i58 < 6; c6_i58++) {
      c6_Phi[(c6_i58 + c6_i56) + 6] = 0.0;
    }

    c6_i56 += 12;
  }

  c6_i59 = 0;
  c6_i60 = 0;
  for (c6_i61 = 0; c6_i61 < 6; c6_i61++) {
    for (c6_i62 = 0; c6_i62 < 6; c6_i62++) {
      c6_Phi[(c6_i62 + c6_i59) + 78] = c6_Phi_t[c6_i62 + c6_i60];
    }

    c6_i59 += 12;
    c6_i60 += 6;
  }

  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  for (c6_i63 = 0; c6_i63 < 144; c6_i63++) {
    (*c6_b_Phi)[c6_i63] = c6_Phi[c6_i63];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 5U, chartInstance->c6_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_02MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c6_i64 = 0; c6_i64 < 144; c6_i64++) {
    _SFD_DATA_RANGE_CHECK((*c6_b_Phi)[c6_i64], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c6_c_n, 2U);
  _SFD_DATA_RANGE_CHECK(*c6_c_t_Kalman, 3U);
  for (c6_i65 = 0; c6_i65 < 3; c6_i65++) {
    _SFD_DATA_RANGE_CHECK((*c6_d_p)[c6_i65], 4U);
  }

  for (c6_i66 = 0; c6_i66 < 9; c6_i66++) {
    _SFD_DATA_RANGE_CHECK((*c6_c_M)[c6_i66], 5U);
  }
}

static void initSimStructsc6_Model_02(SFc6_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c6_fn_VectorToSkewSymmetricTensor(SFc6_Model_02InstanceStruct
  *chartInstance, real_T c6_v[3], real_T c6_SkewSymmetricTensor[9])
{
  uint32_T c6_debug_family_var_map[4];
  real_T c6_nargin = 1.0;
  real_T c6_nargout = 1.0;
  int32_T c6_i67;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c6_b_debug_family_names,
    c6_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_nargin, 0U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_nargout, 1U, c6_d_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_v, 2U, c6_c_sf_marshallOut,
    c6_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_SkewSymmetricTensor, 3U,
    c6_b_sf_marshallOut, c6_e_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 2);
  for (c6_i67 = 0; c6_i67 < 9; c6_i67++) {
    c6_SkewSymmetricTensor[c6_i67] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 3);
  c6_SkewSymmetricTensor[0] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 4);
  c6_SkewSymmetricTensor[4] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 5);
  c6_SkewSymmetricTensor[8] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 6);
  c6_SkewSymmetricTensor[3] = -c6_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 7);
  c6_SkewSymmetricTensor[6] = c6_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 8);
  c6_SkewSymmetricTensor[7] = -c6_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 9);
  c6_SkewSymmetricTensor[1] = c6_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 10);
  c6_SkewSymmetricTensor[2] = -c6_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, 11);
  c6_SkewSymmetricTensor[5] = c6_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c6_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c6_machineNumber, uint32_T
  c6_chartNumber, uint32_T c6_instanceNumber)
{
  (void)c6_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c6_chartNumber, c6_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg-2\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c6_sf_marshallOut(void *chartInstanceVoid, void *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  int32_T c6_i68;
  int32_T c6_i69;
  int32_T c6_i70;
  real_T c6_b_inData[144];
  int32_T c6_i71;
  int32_T c6_i72;
  int32_T c6_i73;
  real_T c6_u[144];
  const mxArray *c6_y = NULL;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_i68 = 0;
  for (c6_i69 = 0; c6_i69 < 12; c6_i69++) {
    for (c6_i70 = 0; c6_i70 < 12; c6_i70++) {
      c6_b_inData[c6_i70 + c6_i68] = (*(real_T (*)[144])c6_inData)[c6_i70 +
        c6_i68];
    }

    c6_i68 += 12;
  }

  c6_i71 = 0;
  for (c6_i72 = 0; c6_i72 < 12; c6_i72++) {
    for (c6_i73 = 0; c6_i73 < 12; c6_i73++) {
      c6_u[c6_i73 + c6_i71] = c6_b_inData[c6_i73 + c6_i71];
    }

    c6_i71 += 12;
  }

  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u, 0, 0U, 1U, 0U, 2, 12, 12), false);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, false);
  return c6_mxArrayOutData;
}

static void c6_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_Phi, const char_T *c6_identifier, real_T c6_y[144])
{
  emlrtMsgIdentifier c6_thisId;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_Phi), &c6_thisId, c6_y);
  sf_mex_destroy(&c6_Phi);
}

static void c6_b_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[144])
{
  real_T c6_dv3[144];
  int32_T c6_i74;
  (void)chartInstance;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), c6_dv3, 1, 0, 0U, 1, 0U, 2, 12,
                12);
  for (c6_i74 = 0; c6_i74 < 144; c6_i74++) {
    c6_y[c6_i74] = c6_dv3[c6_i74];
  }

  sf_mex_destroy(&c6_u);
}

static void c6_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_Phi;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  real_T c6_y[144];
  int32_T c6_i75;
  int32_T c6_i76;
  int32_T c6_i77;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_Phi = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_Phi), &c6_thisId, c6_y);
  sf_mex_destroy(&c6_Phi);
  c6_i75 = 0;
  for (c6_i76 = 0; c6_i76 < 12; c6_i76++) {
    for (c6_i77 = 0; c6_i77 < 12; c6_i77++) {
      (*(real_T (*)[144])c6_outData)[c6_i77 + c6_i75] = c6_y[c6_i77 + c6_i75];
    }

    c6_i75 += 12;
  }

  sf_mex_destroy(&c6_mxArrayInData);
}

static const mxArray *c6_b_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  int32_T c6_i78;
  int32_T c6_i79;
  int32_T c6_i80;
  real_T c6_b_inData[9];
  int32_T c6_i81;
  int32_T c6_i82;
  int32_T c6_i83;
  real_T c6_u[9];
  const mxArray *c6_y = NULL;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_i78 = 0;
  for (c6_i79 = 0; c6_i79 < 3; c6_i79++) {
    for (c6_i80 = 0; c6_i80 < 3; c6_i80++) {
      c6_b_inData[c6_i80 + c6_i78] = (*(real_T (*)[9])c6_inData)[c6_i80 + c6_i78];
    }

    c6_i78 += 3;
  }

  c6_i81 = 0;
  for (c6_i82 = 0; c6_i82 < 3; c6_i82++) {
    for (c6_i83 = 0; c6_i83 < 3; c6_i83++) {
      c6_u[c6_i83 + c6_i81] = c6_b_inData[c6_i83 + c6_i81];
    }

    c6_i81 += 3;
  }

  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, false);
  return c6_mxArrayOutData;
}

static const mxArray *c6_c_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  int32_T c6_i84;
  real_T c6_b_inData[3];
  int32_T c6_i85;
  real_T c6_u[3];
  const mxArray *c6_y = NULL;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  for (c6_i84 = 0; c6_i84 < 3; c6_i84++) {
    c6_b_inData[c6_i84] = (*(real_T (*)[3])c6_inData)[c6_i84];
  }

  for (c6_i85 = 0; c6_i85 < 3; c6_i85++) {
    c6_u[c6_i85] = c6_b_inData[c6_i85];
  }

  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, false);
  return c6_mxArrayOutData;
}

static const mxArray *c6_d_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  real_T c6_u;
  const mxArray *c6_y = NULL;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_u = *(real_T *)c6_inData;
  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", &c6_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, false);
  return c6_mxArrayOutData;
}

static const mxArray *c6_e_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  int32_T c6_i86;
  real_T c6_b_inData[13];
  int32_T c6_i87;
  real_T c6_u[13];
  const mxArray *c6_y = NULL;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  for (c6_i86 = 0; c6_i86 < 13; c6_i86++) {
    c6_b_inData[c6_i86] = (*(real_T (*)[13])c6_inData)[c6_i86];
  }

  for (c6_i87 = 0; c6_i87 < 13; c6_i87++) {
    c6_u[c6_i87] = c6_b_inData[c6_i87];
  }

  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u, 0, 0U, 1U, 0U, 1, 13), false);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, false);
  return c6_mxArrayOutData;
}

static real_T c6_c_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId)
{
  real_T c6_y;
  real_T c6_d0;
  (void)chartInstance;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), &c6_d0, 1, 0, 0U, 0, 0U, 0);
  c6_y = c6_d0;
  sf_mex_destroy(&c6_u);
  return c6_y;
}

static void c6_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_nargout;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  real_T c6_y;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_nargout = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_y = c6_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_nargout), &c6_thisId);
  sf_mex_destroy(&c6_nargout);
  *(real_T *)c6_outData = c6_y;
  sf_mex_destroy(&c6_mxArrayInData);
}

static const mxArray *c6_f_sf_marshallOut(void *chartInstanceVoid, real_T
  c6_inData_data[], int32_T c6_inData_sizes[2])
{
  const mxArray *c6_mxArrayOutData = NULL;
  int32_T c6_b_inData_sizes[2];
  int32_T c6_i88;
  int32_T c6_i89;
  real_T c6_b_inData_data[36];
  int32_T c6_u_sizes[2];
  int32_T c6_i90;
  int32_T c6_i91;
  real_T c6_u_data[36];
  const mxArray *c6_y = NULL;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_b_inData_sizes[0] = 6;
  c6_b_inData_sizes[1] = 6;
  for (c6_i88 = 0; c6_i88 < 6; c6_i88++) {
    for (c6_i89 = 0; c6_i89 < 6; c6_i89++) {
      c6_b_inData_data[c6_i89 + c6_b_inData_sizes[0] * c6_i88] =
        c6_inData_data[c6_i89 + c6_inData_sizes[0] * c6_i88];
    }
  }

  c6_u_sizes[0] = 6;
  c6_u_sizes[1] = 6;
  for (c6_i90 = 0; c6_i90 < 6; c6_i90++) {
    for (c6_i91 = 0; c6_i91 < 6; c6_i91++) {
      c6_u_data[c6_i91 + c6_u_sizes[0] * c6_i90] = c6_b_inData_data[c6_i91 +
        c6_b_inData_sizes[0] * c6_i90];
    }
  }

  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u_data, 0, 0U, 1U, 0U, 2,
    c6_u_sizes[0], c6_u_sizes[1]), false);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, false);
  return c6_mxArrayOutData;
}

static void c6_d_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y_data[],
  int32_T c6_y_sizes[2])
{
  int32_T c6_i92;
  uint32_T c6_uv0[2];
  int32_T c6_i93;
  boolean_T c6_bv0[2];
  int32_T c6_tmp_sizes[2];
  real_T c6_tmp_data[36];
  int32_T c6_y;
  int32_T c6_b_y;
  int32_T c6_i94;
  (void)chartInstance;
  for (c6_i92 = 0; c6_i92 < 2; c6_i92++) {
    c6_uv0[c6_i92] = 6U;
  }

  for (c6_i93 = 0; c6_i93 < 2; c6_i93++) {
    c6_bv0[c6_i93] = false;
  }

  sf_mex_import_vs(c6_parentId, sf_mex_dup(c6_u), c6_tmp_data, 1, 0, 0U, 1, 0U,
                   2, c6_bv0, c6_uv0, c6_tmp_sizes);
  c6_y_sizes[0] = 6;
  c6_y_sizes[1] = 6;
  c6_y = c6_y_sizes[0];
  c6_b_y = c6_y_sizes[1];
  for (c6_i94 = 0; c6_i94 < 36; c6_i94++) {
    c6_y_data[c6_i94] = c6_tmp_data[c6_i94];
  }

  sf_mex_destroy(&c6_u);
}

static void c6_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, real_T c6_outData_data[], int32_T
  c6_outData_sizes[2])
{
  const mxArray *c6_Phi_r;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  int32_T c6_y_sizes[2];
  real_T c6_y_data[36];
  int32_T c6_i95;
  int32_T c6_i96;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_Phi_r = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_Phi_r), &c6_thisId,
                        c6_y_data, c6_y_sizes);
  sf_mex_destroy(&c6_Phi_r);
  c6_outData_sizes[0] = 6;
  c6_outData_sizes[1] = 6;
  for (c6_i95 = 0; c6_i95 < 6; c6_i95++) {
    for (c6_i96 = 0; c6_i96 < 6; c6_i96++) {
      c6_outData_data[c6_i96 + c6_outData_sizes[0] * c6_i95] = c6_y_data[c6_i96
        + c6_y_sizes[0] * c6_i95];
    }
  }

  sf_mex_destroy(&c6_mxArrayInData);
}

static const mxArray *c6_g_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  int32_T c6_i97;
  int32_T c6_i98;
  int32_T c6_i99;
  real_T c6_b_inData[36];
  int32_T c6_i100;
  int32_T c6_i101;
  int32_T c6_i102;
  real_T c6_u[36];
  const mxArray *c6_y = NULL;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_i97 = 0;
  for (c6_i98 = 0; c6_i98 < 6; c6_i98++) {
    for (c6_i99 = 0; c6_i99 < 6; c6_i99++) {
      c6_b_inData[c6_i99 + c6_i97] = (*(real_T (*)[36])c6_inData)[c6_i99 +
        c6_i97];
    }

    c6_i97 += 6;
  }

  c6_i100 = 0;
  for (c6_i101 = 0; c6_i101 < 6; c6_i101++) {
    for (c6_i102 = 0; c6_i102 < 6; c6_i102++) {
      c6_u[c6_i102 + c6_i100] = c6_b_inData[c6_i102 + c6_i100];
    }

    c6_i100 += 6;
  }

  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u, 0, 0U, 1U, 0U, 2, 6, 6), false);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, false);
  return c6_mxArrayOutData;
}

static void c6_e_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[36])
{
  real_T c6_dv4[36];
  int32_T c6_i103;
  (void)chartInstance;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), c6_dv4, 1, 0, 0U, 1, 0U, 2, 6, 6);
  for (c6_i103 = 0; c6_i103 < 36; c6_i103++) {
    c6_y[c6_i103] = c6_dv4[c6_i103];
  }

  sf_mex_destroy(&c6_u);
}

static void c6_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_Phi_t;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  real_T c6_y[36];
  int32_T c6_i104;
  int32_T c6_i105;
  int32_T c6_i106;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_Phi_t = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_Phi_t), &c6_thisId, c6_y);
  sf_mex_destroy(&c6_Phi_t);
  c6_i104 = 0;
  for (c6_i105 = 0; c6_i105 < 6; c6_i105++) {
    for (c6_i106 = 0; c6_i106 < 6; c6_i106++) {
      (*(real_T (*)[36])c6_outData)[c6_i106 + c6_i104] = c6_y[c6_i106 + c6_i104];
    }

    c6_i104 += 6;
  }

  sf_mex_destroy(&c6_mxArrayInData);
}

static void c6_f_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[9])
{
  real_T c6_dv5[9];
  int32_T c6_i107;
  (void)chartInstance;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), c6_dv5, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c6_i107 = 0; c6_i107 < 9; c6_i107++) {
    c6_y[c6_i107] = c6_dv5[c6_i107];
  }

  sf_mex_destroy(&c6_u);
}

static void c6_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_SkewSymmetricTensor;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  real_T c6_y[9];
  int32_T c6_i108;
  int32_T c6_i109;
  int32_T c6_i110;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_SkewSymmetricTensor = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_SkewSymmetricTensor),
                        &c6_thisId, c6_y);
  sf_mex_destroy(&c6_SkewSymmetricTensor);
  c6_i108 = 0;
  for (c6_i109 = 0; c6_i109 < 3; c6_i109++) {
    for (c6_i110 = 0; c6_i110 < 3; c6_i110++) {
      (*(real_T (*)[9])c6_outData)[c6_i110 + c6_i108] = c6_y[c6_i110 + c6_i108];
    }

    c6_i108 += 3;
  }

  sf_mex_destroy(&c6_mxArrayInData);
}

static void c6_g_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[3])
{
  real_T c6_dv6[3];
  int32_T c6_i111;
  (void)chartInstance;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), c6_dv6, 1, 0, 0U, 1, 0U, 1, 3);
  for (c6_i111 = 0; c6_i111 < 3; c6_i111++) {
    c6_y[c6_i111] = c6_dv6[c6_i111];
  }

  sf_mex_destroy(&c6_u);
}

static void c6_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_v;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  real_T c6_y[3];
  int32_T c6_i112;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_v = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_v), &c6_thisId, c6_y);
  sf_mex_destroy(&c6_v);
  for (c6_i112 = 0; c6_i112 < 3; c6_i112++) {
    (*(real_T (*)[3])c6_outData)[c6_i112] = c6_y[c6_i112];
  }

  sf_mex_destroy(&c6_mxArrayInData);
}

static void c6_h_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[13])
{
  real_T c6_dv7[13];
  int32_T c6_i113;
  (void)chartInstance;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), c6_dv7, 1, 0, 0U, 1, 0U, 1, 13);
  for (c6_i113 = 0; c6_i113 < 13; c6_i113++) {
    c6_y[c6_i113] = c6_dv7[c6_i113];
  }

  sf_mex_destroy(&c6_u);
}

static void c6_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_X_a;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  real_T c6_y[13];
  int32_T c6_i114;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_X_a = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_X_a), &c6_thisId, c6_y);
  sf_mex_destroy(&c6_X_a);
  for (c6_i114 = 0; c6_i114 < 13; c6_i114++) {
    (*(real_T (*)[13])c6_outData)[c6_i114] = c6_y[c6_i114];
  }

  sf_mex_destroy(&c6_mxArrayInData);
}

const mxArray *sf_c6_Model_02_get_eml_resolved_functions_info(void)
{
  const mxArray *c6_nameCaptureInfo = NULL;
  c6_nameCaptureInfo = NULL;
  sf_mex_assign(&c6_nameCaptureInfo, sf_mex_createstruct("structure", 2, 184, 1),
                false);
  c6_info_helper(&c6_nameCaptureInfo);
  c6_b_info_helper(&c6_nameCaptureInfo);
  c6_c_info_helper(&c6_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c6_nameCaptureInfo);
  return c6_nameCaptureInfo;
}

static void c6_info_helper(const mxArray **c6_info)
{
  const mxArray *c6_rhs0 = NULL;
  const mxArray *c6_lhs0 = NULL;
  const mxArray *c6_rhs1 = NULL;
  const mxArray *c6_lhs1 = NULL;
  const mxArray *c6_rhs2 = NULL;
  const mxArray *c6_lhs2 = NULL;
  const mxArray *c6_rhs3 = NULL;
  const mxArray *c6_lhs3 = NULL;
  const mxArray *c6_rhs4 = NULL;
  const mxArray *c6_lhs4 = NULL;
  const mxArray *c6_rhs5 = NULL;
  const mxArray *c6_lhs5 = NULL;
  const mxArray *c6_rhs6 = NULL;
  const mxArray *c6_lhs6 = NULL;
  const mxArray *c6_rhs7 = NULL;
  const mxArray *c6_lhs7 = NULL;
  const mxArray *c6_rhs8 = NULL;
  const mxArray *c6_lhs8 = NULL;
  const mxArray *c6_rhs9 = NULL;
  const mxArray *c6_lhs9 = NULL;
  const mxArray *c6_rhs10 = NULL;
  const mxArray *c6_lhs10 = NULL;
  const mxArray *c6_rhs11 = NULL;
  const mxArray *c6_lhs11 = NULL;
  const mxArray *c6_rhs12 = NULL;
  const mxArray *c6_lhs12 = NULL;
  const mxArray *c6_rhs13 = NULL;
  const mxArray *c6_lhs13 = NULL;
  const mxArray *c6_rhs14 = NULL;
  const mxArray *c6_lhs14 = NULL;
  const mxArray *c6_rhs15 = NULL;
  const mxArray *c6_lhs15 = NULL;
  const mxArray *c6_rhs16 = NULL;
  const mxArray *c6_lhs16 = NULL;
  const mxArray *c6_rhs17 = NULL;
  const mxArray *c6_lhs17 = NULL;
  const mxArray *c6_rhs18 = NULL;
  const mxArray *c6_lhs18 = NULL;
  const mxArray *c6_rhs19 = NULL;
  const mxArray *c6_lhs19 = NULL;
  const mxArray *c6_rhs20 = NULL;
  const mxArray *c6_lhs20 = NULL;
  const mxArray *c6_rhs21 = NULL;
  const mxArray *c6_lhs21 = NULL;
  const mxArray *c6_rhs22 = NULL;
  const mxArray *c6_lhs22 = NULL;
  const mxArray *c6_rhs23 = NULL;
  const mxArray *c6_lhs23 = NULL;
  const mxArray *c6_rhs24 = NULL;
  const mxArray *c6_lhs24 = NULL;
  const mxArray *c6_rhs25 = NULL;
  const mxArray *c6_lhs25 = NULL;
  const mxArray *c6_rhs26 = NULL;
  const mxArray *c6_lhs26 = NULL;
  const mxArray *c6_rhs27 = NULL;
  const mxArray *c6_lhs27 = NULL;
  const mxArray *c6_rhs28 = NULL;
  const mxArray *c6_lhs28 = NULL;
  const mxArray *c6_rhs29 = NULL;
  const mxArray *c6_lhs29 = NULL;
  const mxArray *c6_rhs30 = NULL;
  const mxArray *c6_lhs30 = NULL;
  const mxArray *c6_rhs31 = NULL;
  const mxArray *c6_lhs31 = NULL;
  const mxArray *c6_rhs32 = NULL;
  const mxArray *c6_lhs32 = NULL;
  const mxArray *c6_rhs33 = NULL;
  const mxArray *c6_lhs33 = NULL;
  const mxArray *c6_rhs34 = NULL;
  const mxArray *c6_lhs34 = NULL;
  const mxArray *c6_rhs35 = NULL;
  const mxArray *c6_lhs35 = NULL;
  const mxArray *c6_rhs36 = NULL;
  const mxArray *c6_lhs36 = NULL;
  const mxArray *c6_rhs37 = NULL;
  const mxArray *c6_lhs37 = NULL;
  const mxArray *c6_rhs38 = NULL;
  const mxArray *c6_lhs38 = NULL;
  const mxArray *c6_rhs39 = NULL;
  const mxArray *c6_lhs39 = NULL;
  const mxArray *c6_rhs40 = NULL;
  const mxArray *c6_lhs40 = NULL;
  const mxArray *c6_rhs41 = NULL;
  const mxArray *c6_lhs41 = NULL;
  const mxArray *c6_rhs42 = NULL;
  const mxArray *c6_lhs42 = NULL;
  const mxArray *c6_rhs43 = NULL;
  const mxArray *c6_lhs43 = NULL;
  const mxArray *c6_rhs44 = NULL;
  const mxArray *c6_lhs44 = NULL;
  const mxArray *c6_rhs45 = NULL;
  const mxArray *c6_lhs45 = NULL;
  const mxArray *c6_rhs46 = NULL;
  const mxArray *c6_lhs46 = NULL;
  const mxArray *c6_rhs47 = NULL;
  const mxArray *c6_lhs47 = NULL;
  const mxArray *c6_rhs48 = NULL;
  const mxArray *c6_lhs48 = NULL;
  const mxArray *c6_rhs49 = NULL;
  const mxArray *c6_lhs49 = NULL;
  const mxArray *c6_rhs50 = NULL;
  const mxArray *c6_lhs50 = NULL;
  const mxArray *c6_rhs51 = NULL;
  const mxArray *c6_lhs51 = NULL;
  const mxArray *c6_rhs52 = NULL;
  const mxArray *c6_lhs52 = NULL;
  const mxArray *c6_rhs53 = NULL;
  const mxArray *c6_lhs53 = NULL;
  const mxArray *c6_rhs54 = NULL;
  const mxArray *c6_lhs54 = NULL;
  const mxArray *c6_rhs55 = NULL;
  const mxArray *c6_lhs55 = NULL;
  const mxArray *c6_rhs56 = NULL;
  const mxArray *c6_lhs56 = NULL;
  const mxArray *c6_rhs57 = NULL;
  const mxArray *c6_lhs57 = NULL;
  const mxArray *c6_rhs58 = NULL;
  const mxArray *c6_lhs58 = NULL;
  const mxArray *c6_rhs59 = NULL;
  const mxArray *c6_lhs59 = NULL;
  const mxArray *c6_rhs60 = NULL;
  const mxArray *c6_lhs60 = NULL;
  const mxArray *c6_rhs61 = NULL;
  const mxArray *c6_lhs61 = NULL;
  const mxArray *c6_rhs62 = NULL;
  const mxArray *c6_lhs62 = NULL;
  const mxArray *c6_rhs63 = NULL;
  const mxArray *c6_lhs63 = NULL;
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("mpower"), "name", "name", 0);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c6_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c6_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("ismatrix"), "name", "name", 2);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c6_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("power"), "name", "name", 3);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c6_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c6_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 5);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 5);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c6_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 6);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c6_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 7);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 7);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c6_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 8);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 8);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c6_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 9);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("floor"), "name", "name", 9);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c6_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 10);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c6_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 11);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c6_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 12);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c6_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "context", "context", 13);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eye"), "name", "name", 13);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c6_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 14);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c6_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 15);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 15);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c6_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 16);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("isinf"), "name", "name", 16);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c6_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 17);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c6_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 18);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 18);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c6_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 19);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmax"), "name", "name", 19);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c6_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 20);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c6_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 21);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmin"), "name", "name", 21);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c6_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 22);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c6_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 23);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 23);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c6_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 24);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 24);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c6_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 25);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 25);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c6_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 26);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmin"), "name", "name", 26);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c6_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 27);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 27);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c6_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 28);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmax"), "name", "name", 28);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c6_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 29);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 29);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c6_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 30);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmax"), "name", "name", 30);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 30);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c6_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "context", "context", 31);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 31);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1450040424U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c6_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "context", "context", 32);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 32);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c6_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 33);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 33);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c6_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "context", "context", 34);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("expm"), "name", "name", 34);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "resolved",
                  "resolved", 34);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1381857504U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c6_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 35);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c6_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 36);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("ismatrix"), "name", "name", 36);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 36);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c6_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 37);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("norm"), "name", "name", 37);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 37);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c6_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 38);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c6_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 39);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("abs"), "name", "name", 39);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c6_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 40);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c6_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 41);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c6_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 42);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("isnan"), "name", "name", 42);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c6_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 43);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 43);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c6_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 44);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 44);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c6_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "context", "context", 45);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 45);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c6_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m!PadeApproximantOfDegree"),
                  "context", "context", 46);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 46);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c6_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 47);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 47);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c6_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 48);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 48);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 48);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c6_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 49);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  49);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c6_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 50);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 50);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c6_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 51);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 51);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c6_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 52);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 52);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c6_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 53);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 53);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c6_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 54);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 54);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c6_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 55);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 55);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c6_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 56);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 56);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c6_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m!PadeApproximantOfDegree"),
                  "context", "context", 57);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("mldivide"), "name", "name", 57);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "resolved",
                  "resolved", 57);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1319737166U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c6_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "context",
                  "context", 58);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_lusolve"), "name", "name",
                  58);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "resolved",
                  "resolved", 58);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1370017086U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c6_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 59);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_xgetrf"), "name", "name",
                  59);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286826006U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c6_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "context", "context", 60);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_lapack_xgetrf"), "name",
                  "name", 60);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286826010U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c6_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "context", "context", 61);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_matlab_zgetrf"), "name",
                  "name", 61);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1302696194U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c6_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 62);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("realmin"), "name", "name", 62);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 62);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c6_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 63);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_realmin"), "name", "name",
                  63);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 63);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1307658444U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c6_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c6_rhs0);
  sf_mex_destroy(&c6_lhs0);
  sf_mex_destroy(&c6_rhs1);
  sf_mex_destroy(&c6_lhs1);
  sf_mex_destroy(&c6_rhs2);
  sf_mex_destroy(&c6_lhs2);
  sf_mex_destroy(&c6_rhs3);
  sf_mex_destroy(&c6_lhs3);
  sf_mex_destroy(&c6_rhs4);
  sf_mex_destroy(&c6_lhs4);
  sf_mex_destroy(&c6_rhs5);
  sf_mex_destroy(&c6_lhs5);
  sf_mex_destroy(&c6_rhs6);
  sf_mex_destroy(&c6_lhs6);
  sf_mex_destroy(&c6_rhs7);
  sf_mex_destroy(&c6_lhs7);
  sf_mex_destroy(&c6_rhs8);
  sf_mex_destroy(&c6_lhs8);
  sf_mex_destroy(&c6_rhs9);
  sf_mex_destroy(&c6_lhs9);
  sf_mex_destroy(&c6_rhs10);
  sf_mex_destroy(&c6_lhs10);
  sf_mex_destroy(&c6_rhs11);
  sf_mex_destroy(&c6_lhs11);
  sf_mex_destroy(&c6_rhs12);
  sf_mex_destroy(&c6_lhs12);
  sf_mex_destroy(&c6_rhs13);
  sf_mex_destroy(&c6_lhs13);
  sf_mex_destroy(&c6_rhs14);
  sf_mex_destroy(&c6_lhs14);
  sf_mex_destroy(&c6_rhs15);
  sf_mex_destroy(&c6_lhs15);
  sf_mex_destroy(&c6_rhs16);
  sf_mex_destroy(&c6_lhs16);
  sf_mex_destroy(&c6_rhs17);
  sf_mex_destroy(&c6_lhs17);
  sf_mex_destroy(&c6_rhs18);
  sf_mex_destroy(&c6_lhs18);
  sf_mex_destroy(&c6_rhs19);
  sf_mex_destroy(&c6_lhs19);
  sf_mex_destroy(&c6_rhs20);
  sf_mex_destroy(&c6_lhs20);
  sf_mex_destroy(&c6_rhs21);
  sf_mex_destroy(&c6_lhs21);
  sf_mex_destroy(&c6_rhs22);
  sf_mex_destroy(&c6_lhs22);
  sf_mex_destroy(&c6_rhs23);
  sf_mex_destroy(&c6_lhs23);
  sf_mex_destroy(&c6_rhs24);
  sf_mex_destroy(&c6_lhs24);
  sf_mex_destroy(&c6_rhs25);
  sf_mex_destroy(&c6_lhs25);
  sf_mex_destroy(&c6_rhs26);
  sf_mex_destroy(&c6_lhs26);
  sf_mex_destroy(&c6_rhs27);
  sf_mex_destroy(&c6_lhs27);
  sf_mex_destroy(&c6_rhs28);
  sf_mex_destroy(&c6_lhs28);
  sf_mex_destroy(&c6_rhs29);
  sf_mex_destroy(&c6_lhs29);
  sf_mex_destroy(&c6_rhs30);
  sf_mex_destroy(&c6_lhs30);
  sf_mex_destroy(&c6_rhs31);
  sf_mex_destroy(&c6_lhs31);
  sf_mex_destroy(&c6_rhs32);
  sf_mex_destroy(&c6_lhs32);
  sf_mex_destroy(&c6_rhs33);
  sf_mex_destroy(&c6_lhs33);
  sf_mex_destroy(&c6_rhs34);
  sf_mex_destroy(&c6_lhs34);
  sf_mex_destroy(&c6_rhs35);
  sf_mex_destroy(&c6_lhs35);
  sf_mex_destroy(&c6_rhs36);
  sf_mex_destroy(&c6_lhs36);
  sf_mex_destroy(&c6_rhs37);
  sf_mex_destroy(&c6_lhs37);
  sf_mex_destroy(&c6_rhs38);
  sf_mex_destroy(&c6_lhs38);
  sf_mex_destroy(&c6_rhs39);
  sf_mex_destroy(&c6_lhs39);
  sf_mex_destroy(&c6_rhs40);
  sf_mex_destroy(&c6_lhs40);
  sf_mex_destroy(&c6_rhs41);
  sf_mex_destroy(&c6_lhs41);
  sf_mex_destroy(&c6_rhs42);
  sf_mex_destroy(&c6_lhs42);
  sf_mex_destroy(&c6_rhs43);
  sf_mex_destroy(&c6_lhs43);
  sf_mex_destroy(&c6_rhs44);
  sf_mex_destroy(&c6_lhs44);
  sf_mex_destroy(&c6_rhs45);
  sf_mex_destroy(&c6_lhs45);
  sf_mex_destroy(&c6_rhs46);
  sf_mex_destroy(&c6_lhs46);
  sf_mex_destroy(&c6_rhs47);
  sf_mex_destroy(&c6_lhs47);
  sf_mex_destroy(&c6_rhs48);
  sf_mex_destroy(&c6_lhs48);
  sf_mex_destroy(&c6_rhs49);
  sf_mex_destroy(&c6_lhs49);
  sf_mex_destroy(&c6_rhs50);
  sf_mex_destroy(&c6_lhs50);
  sf_mex_destroy(&c6_rhs51);
  sf_mex_destroy(&c6_lhs51);
  sf_mex_destroy(&c6_rhs52);
  sf_mex_destroy(&c6_lhs52);
  sf_mex_destroy(&c6_rhs53);
  sf_mex_destroy(&c6_lhs53);
  sf_mex_destroy(&c6_rhs54);
  sf_mex_destroy(&c6_lhs54);
  sf_mex_destroy(&c6_rhs55);
  sf_mex_destroy(&c6_lhs55);
  sf_mex_destroy(&c6_rhs56);
  sf_mex_destroy(&c6_lhs56);
  sf_mex_destroy(&c6_rhs57);
  sf_mex_destroy(&c6_lhs57);
  sf_mex_destroy(&c6_rhs58);
  sf_mex_destroy(&c6_lhs58);
  sf_mex_destroy(&c6_rhs59);
  sf_mex_destroy(&c6_lhs59);
  sf_mex_destroy(&c6_rhs60);
  sf_mex_destroy(&c6_lhs60);
  sf_mex_destroy(&c6_rhs61);
  sf_mex_destroy(&c6_lhs61);
  sf_mex_destroy(&c6_rhs62);
  sf_mex_destroy(&c6_lhs62);
  sf_mex_destroy(&c6_rhs63);
  sf_mex_destroy(&c6_lhs63);
}

static const mxArray *c6_emlrt_marshallOut(const char * c6_u)
{
  const mxArray *c6_y = NULL;
  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c6_u)), false);
  return c6_y;
}

static const mxArray *c6_b_emlrt_marshallOut(const uint32_T c6_u)
{
  const mxArray *c6_y = NULL;
  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", &c6_u, 7, 0U, 0U, 0U, 0), false);
  return c6_y;
}

static void c6_b_info_helper(const mxArray **c6_info)
{
  const mxArray *c6_rhs64 = NULL;
  const mxArray *c6_lhs64 = NULL;
  const mxArray *c6_rhs65 = NULL;
  const mxArray *c6_lhs65 = NULL;
  const mxArray *c6_rhs66 = NULL;
  const mxArray *c6_lhs66 = NULL;
  const mxArray *c6_rhs67 = NULL;
  const mxArray *c6_lhs67 = NULL;
  const mxArray *c6_rhs68 = NULL;
  const mxArray *c6_lhs68 = NULL;
  const mxArray *c6_rhs69 = NULL;
  const mxArray *c6_lhs69 = NULL;
  const mxArray *c6_rhs70 = NULL;
  const mxArray *c6_lhs70 = NULL;
  const mxArray *c6_rhs71 = NULL;
  const mxArray *c6_lhs71 = NULL;
  const mxArray *c6_rhs72 = NULL;
  const mxArray *c6_lhs72 = NULL;
  const mxArray *c6_rhs73 = NULL;
  const mxArray *c6_lhs73 = NULL;
  const mxArray *c6_rhs74 = NULL;
  const mxArray *c6_lhs74 = NULL;
  const mxArray *c6_rhs75 = NULL;
  const mxArray *c6_lhs75 = NULL;
  const mxArray *c6_rhs76 = NULL;
  const mxArray *c6_lhs76 = NULL;
  const mxArray *c6_rhs77 = NULL;
  const mxArray *c6_lhs77 = NULL;
  const mxArray *c6_rhs78 = NULL;
  const mxArray *c6_lhs78 = NULL;
  const mxArray *c6_rhs79 = NULL;
  const mxArray *c6_lhs79 = NULL;
  const mxArray *c6_rhs80 = NULL;
  const mxArray *c6_lhs80 = NULL;
  const mxArray *c6_rhs81 = NULL;
  const mxArray *c6_lhs81 = NULL;
  const mxArray *c6_rhs82 = NULL;
  const mxArray *c6_lhs82 = NULL;
  const mxArray *c6_rhs83 = NULL;
  const mxArray *c6_lhs83 = NULL;
  const mxArray *c6_rhs84 = NULL;
  const mxArray *c6_lhs84 = NULL;
  const mxArray *c6_rhs85 = NULL;
  const mxArray *c6_lhs85 = NULL;
  const mxArray *c6_rhs86 = NULL;
  const mxArray *c6_lhs86 = NULL;
  const mxArray *c6_rhs87 = NULL;
  const mxArray *c6_lhs87 = NULL;
  const mxArray *c6_rhs88 = NULL;
  const mxArray *c6_lhs88 = NULL;
  const mxArray *c6_rhs89 = NULL;
  const mxArray *c6_lhs89 = NULL;
  const mxArray *c6_rhs90 = NULL;
  const mxArray *c6_lhs90 = NULL;
  const mxArray *c6_rhs91 = NULL;
  const mxArray *c6_lhs91 = NULL;
  const mxArray *c6_rhs92 = NULL;
  const mxArray *c6_lhs92 = NULL;
  const mxArray *c6_rhs93 = NULL;
  const mxArray *c6_lhs93 = NULL;
  const mxArray *c6_rhs94 = NULL;
  const mxArray *c6_lhs94 = NULL;
  const mxArray *c6_rhs95 = NULL;
  const mxArray *c6_lhs95 = NULL;
  const mxArray *c6_rhs96 = NULL;
  const mxArray *c6_lhs96 = NULL;
  const mxArray *c6_rhs97 = NULL;
  const mxArray *c6_lhs97 = NULL;
  const mxArray *c6_rhs98 = NULL;
  const mxArray *c6_lhs98 = NULL;
  const mxArray *c6_rhs99 = NULL;
  const mxArray *c6_lhs99 = NULL;
  const mxArray *c6_rhs100 = NULL;
  const mxArray *c6_lhs100 = NULL;
  const mxArray *c6_rhs101 = NULL;
  const mxArray *c6_lhs101 = NULL;
  const mxArray *c6_rhs102 = NULL;
  const mxArray *c6_lhs102 = NULL;
  const mxArray *c6_rhs103 = NULL;
  const mxArray *c6_lhs103 = NULL;
  const mxArray *c6_rhs104 = NULL;
  const mxArray *c6_lhs104 = NULL;
  const mxArray *c6_rhs105 = NULL;
  const mxArray *c6_lhs105 = NULL;
  const mxArray *c6_rhs106 = NULL;
  const mxArray *c6_lhs106 = NULL;
  const mxArray *c6_rhs107 = NULL;
  const mxArray *c6_lhs107 = NULL;
  const mxArray *c6_rhs108 = NULL;
  const mxArray *c6_lhs108 = NULL;
  const mxArray *c6_rhs109 = NULL;
  const mxArray *c6_lhs109 = NULL;
  const mxArray *c6_rhs110 = NULL;
  const mxArray *c6_lhs110 = NULL;
  const mxArray *c6_rhs111 = NULL;
  const mxArray *c6_lhs111 = NULL;
  const mxArray *c6_rhs112 = NULL;
  const mxArray *c6_lhs112 = NULL;
  const mxArray *c6_rhs113 = NULL;
  const mxArray *c6_lhs113 = NULL;
  const mxArray *c6_rhs114 = NULL;
  const mxArray *c6_lhs114 = NULL;
  const mxArray *c6_rhs115 = NULL;
  const mxArray *c6_lhs115 = NULL;
  const mxArray *c6_rhs116 = NULL;
  const mxArray *c6_lhs116 = NULL;
  const mxArray *c6_rhs117 = NULL;
  const mxArray *c6_lhs117 = NULL;
  const mxArray *c6_rhs118 = NULL;
  const mxArray *c6_lhs118 = NULL;
  const mxArray *c6_rhs119 = NULL;
  const mxArray *c6_lhs119 = NULL;
  const mxArray *c6_rhs120 = NULL;
  const mxArray *c6_lhs120 = NULL;
  const mxArray *c6_rhs121 = NULL;
  const mxArray *c6_lhs121 = NULL;
  const mxArray *c6_rhs122 = NULL;
  const mxArray *c6_lhs122 = NULL;
  const mxArray *c6_rhs123 = NULL;
  const mxArray *c6_lhs123 = NULL;
  const mxArray *c6_rhs124 = NULL;
  const mxArray *c6_lhs124 = NULL;
  const mxArray *c6_rhs125 = NULL;
  const mxArray *c6_lhs125 = NULL;
  const mxArray *c6_rhs126 = NULL;
  const mxArray *c6_lhs126 = NULL;
  const mxArray *c6_rhs127 = NULL;
  const mxArray *c6_lhs127 = NULL;
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 64);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 64);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c6_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 65);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eps"), "name", "name", 65);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 65);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c6_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 66);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 66);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 66);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c6_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 67);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_eps"), "name", "name", 67);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 67);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c6_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 68);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 68);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c6_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 69);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("min"), "name", "name", 69);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 69);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 69);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1311262518U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c6_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "context",
                  "context", 70);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 70);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1378303184U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c6_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 71);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 71);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 71);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 71);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c6_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 72);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 72);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 72);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c6_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 73);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 73);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 73);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c6_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 74);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 74);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 74);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c6_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 75);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 75);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 75);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c6_rhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 76);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 76);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 76);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 76);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c6_rhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 77);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 77);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 77);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 77);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c6_rhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("colon"), "name", "name", 78);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 78);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 78);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1378303188U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c6_rhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 79);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("colon"), "name", "name", 79);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 79);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1378303188U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c6_rhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 80);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 80);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 80);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c6_rhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 81);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 81);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 81);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 81);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c6_rhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 82);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("floor"), "name", "name", 82);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 82);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c6_rhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 83);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmin"), "name", "name", 83);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 83);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c6_rhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 84);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmax"), "name", "name", 84);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 84);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c6_rhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 85);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmin"), "name", "name", 85);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 85);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c6_rhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 86);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmax"), "name", "name", 86);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 86);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c6_rhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 87);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_isa_uint"), "name", "name",
                  87);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 87);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 87);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c6_rhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "context",
                  "context", 88);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.isaUint"),
                  "name", "name", 88);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 88);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/isaUint.p"),
                  "resolved", "resolved", 88);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c6_rhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 89);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_unsigned_class"), "name",
                  "name", 89);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c6_rhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "context", "context", 90);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.unsignedClass"),
                  "name", "name", 90);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "resolved", "resolved", 90);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c6_rhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "context", "context", 91);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 91);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c6_rhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "context", "context", 92);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 92);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 92);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c6_rhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 93);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 93);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 93);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c6_rhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 94);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmax"), "name", "name", 94);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 94);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 94);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c6_rhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 95);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_isa_uint"), "name", "name",
                  95);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 95);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 95);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c6_rhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 96);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 96);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 96);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 96);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c6_rhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 97);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 97);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 97);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 97);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c6_rhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_signed_integer_colon"),
                  "context", "context", 98);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 98);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 98);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 98);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c6_rhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 99);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 99);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 99);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c6_rhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 100);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 100);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 100);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c6_rhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs100), "rhs", "rhs",
                  100);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs100), "lhs", "lhs",
                  100);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 101);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 101);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 101);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 101);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c6_rhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs101), "rhs", "rhs",
                  101);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs101), "lhs", "lhs",
                  101);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 102);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 102);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 102);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c6_rhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs102), "rhs", "rhs",
                  102);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs102), "lhs", "lhs",
                  102);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 103);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 103);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 103);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 103);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c6_rhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs103), "rhs", "rhs",
                  103);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs103), "lhs", "lhs",
                  103);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 104);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 104);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 104);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 104);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c6_rhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs104), "rhs", "rhs",
                  104);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs104), "lhs", "lhs",
                  104);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 105);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 105);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 105);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 105);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c6_rhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs105), "rhs", "rhs",
                  105);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs105), "lhs", "lhs",
                  105);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 106);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 106);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 106);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 106);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c6_rhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs106), "rhs", "rhs",
                  106);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs106), "lhs", "lhs",
                  106);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 107);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 107);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 107);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 107);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c6_rhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs107), "rhs", "rhs",
                  107);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs107), "lhs", "lhs",
                  107);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 108);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 108);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 108);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c6_rhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs108), "rhs", "rhs",
                  108);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs108), "lhs", "lhs",
                  108);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 109);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 109);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 109);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 109);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c6_rhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs109), "rhs", "rhs",
                  109);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs109), "lhs", "lhs",
                  109);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 110);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_ixamax"), "name", "name",
                  110);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "resolved", "resolved", 110);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c6_rhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs110), "rhs", "rhs",
                  110);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs110), "lhs", "lhs",
                  110);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 111);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 111);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 111);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 111);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c6_rhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs111), "rhs", "rhs",
                  111);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs111), "lhs", "lhs",
                  111);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 112);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.ixamax"),
                  "name", "name", 112);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "resolved", "resolved", 112);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c6_rhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs112), "rhs", "rhs",
                  112);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs112), "lhs", "lhs",
                  112);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "context", "context", 113);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 113);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c6_rhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs113), "rhs", "rhs",
                  113);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs113), "lhs", "lhs",
                  113);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p!below_threshold"),
                  "context", "context", 114);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 114);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 114);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c6_rhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs114), "rhs", "rhs",
                  114);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs114), "lhs", "lhs",
                  114);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p!below_threshold"),
                  "context", "context", 115);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("length"), "name", "name", 115);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 115);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1303153406U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c6_rhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs115), "rhs", "rhs",
                  115);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs115), "lhs", "lhs",
                  115);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength"),
                  "context", "context", 116);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 116);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 116);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 116);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c6_rhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs116), "rhs", "rhs",
                  116);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs116), "lhs", "lhs",
                  116);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "context", "context", 117);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.refblas.ixamax"),
                  "name", "name", 117);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 117);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "resolved", "resolved", 117);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c6_rhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs117), "rhs", "rhs",
                  117);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs117), "lhs", "lhs",
                  117);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 118);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.refblas.xcabs1"),
                  "name", "name", 118);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 118);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "resolved", "resolved", 118);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c6_rhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs118), "rhs", "rhs",
                  118);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs118), "lhs", "lhs",
                  118);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "context", "context", 119);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("abs"), "name", "name", 119);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 119);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 119);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c6_rhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs119), "rhs", "rhs",
                  119);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs119), "lhs", "lhs",
                  119);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 120);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 120);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 120);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 120);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c6_rhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs120), "rhs", "rhs",
                  120);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs120), "lhs", "lhs",
                  120);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 121);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 121);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 121);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 121);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c6_rhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs121), "rhs", "rhs",
                  121);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs121), "lhs", "lhs",
                  121);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 122);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_xswap"), "name", "name",
                  122);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"),
                  "resolved", "resolved", 122);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c6_rhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs122), "rhs", "rhs",
                  122);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs122), "lhs", "lhs",
                  122);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 123);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 123);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 123);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 123);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c6_rhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs123), "rhs", "rhs",
                  123);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs123), "lhs", "lhs",
                  123);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 124);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.xswap"),
                  "name", "name", 124);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 124);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "resolved", "resolved", 124);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c6_rhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs124), "rhs", "rhs",
                  124);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs124), "lhs", "lhs",
                  124);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 125);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 125);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 125);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 125);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c6_rhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs125), "rhs", "rhs",
                  125);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs125), "lhs", "lhs",
                  125);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p!below_threshold"),
                  "context", "context", 126);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 126);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 126);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c6_rhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs126), "rhs", "rhs",
                  126);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs126), "lhs", "lhs",
                  126);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 127);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.refblas.xswap"),
                  "name", "name", 127);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "resolved", "resolved", 127);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c6_rhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs127), "rhs", "rhs",
                  127);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs127), "lhs", "lhs",
                  127);
  sf_mex_destroy(&c6_rhs64);
  sf_mex_destroy(&c6_lhs64);
  sf_mex_destroy(&c6_rhs65);
  sf_mex_destroy(&c6_lhs65);
  sf_mex_destroy(&c6_rhs66);
  sf_mex_destroy(&c6_lhs66);
  sf_mex_destroy(&c6_rhs67);
  sf_mex_destroy(&c6_lhs67);
  sf_mex_destroy(&c6_rhs68);
  sf_mex_destroy(&c6_lhs68);
  sf_mex_destroy(&c6_rhs69);
  sf_mex_destroy(&c6_lhs69);
  sf_mex_destroy(&c6_rhs70);
  sf_mex_destroy(&c6_lhs70);
  sf_mex_destroy(&c6_rhs71);
  sf_mex_destroy(&c6_lhs71);
  sf_mex_destroy(&c6_rhs72);
  sf_mex_destroy(&c6_lhs72);
  sf_mex_destroy(&c6_rhs73);
  sf_mex_destroy(&c6_lhs73);
  sf_mex_destroy(&c6_rhs74);
  sf_mex_destroy(&c6_lhs74);
  sf_mex_destroy(&c6_rhs75);
  sf_mex_destroy(&c6_lhs75);
  sf_mex_destroy(&c6_rhs76);
  sf_mex_destroy(&c6_lhs76);
  sf_mex_destroy(&c6_rhs77);
  sf_mex_destroy(&c6_lhs77);
  sf_mex_destroy(&c6_rhs78);
  sf_mex_destroy(&c6_lhs78);
  sf_mex_destroy(&c6_rhs79);
  sf_mex_destroy(&c6_lhs79);
  sf_mex_destroy(&c6_rhs80);
  sf_mex_destroy(&c6_lhs80);
  sf_mex_destroy(&c6_rhs81);
  sf_mex_destroy(&c6_lhs81);
  sf_mex_destroy(&c6_rhs82);
  sf_mex_destroy(&c6_lhs82);
  sf_mex_destroy(&c6_rhs83);
  sf_mex_destroy(&c6_lhs83);
  sf_mex_destroy(&c6_rhs84);
  sf_mex_destroy(&c6_lhs84);
  sf_mex_destroy(&c6_rhs85);
  sf_mex_destroy(&c6_lhs85);
  sf_mex_destroy(&c6_rhs86);
  sf_mex_destroy(&c6_lhs86);
  sf_mex_destroy(&c6_rhs87);
  sf_mex_destroy(&c6_lhs87);
  sf_mex_destroy(&c6_rhs88);
  sf_mex_destroy(&c6_lhs88);
  sf_mex_destroy(&c6_rhs89);
  sf_mex_destroy(&c6_lhs89);
  sf_mex_destroy(&c6_rhs90);
  sf_mex_destroy(&c6_lhs90);
  sf_mex_destroy(&c6_rhs91);
  sf_mex_destroy(&c6_lhs91);
  sf_mex_destroy(&c6_rhs92);
  sf_mex_destroy(&c6_lhs92);
  sf_mex_destroy(&c6_rhs93);
  sf_mex_destroy(&c6_lhs93);
  sf_mex_destroy(&c6_rhs94);
  sf_mex_destroy(&c6_lhs94);
  sf_mex_destroy(&c6_rhs95);
  sf_mex_destroy(&c6_lhs95);
  sf_mex_destroy(&c6_rhs96);
  sf_mex_destroy(&c6_lhs96);
  sf_mex_destroy(&c6_rhs97);
  sf_mex_destroy(&c6_lhs97);
  sf_mex_destroy(&c6_rhs98);
  sf_mex_destroy(&c6_lhs98);
  sf_mex_destroy(&c6_rhs99);
  sf_mex_destroy(&c6_lhs99);
  sf_mex_destroy(&c6_rhs100);
  sf_mex_destroy(&c6_lhs100);
  sf_mex_destroy(&c6_rhs101);
  sf_mex_destroy(&c6_lhs101);
  sf_mex_destroy(&c6_rhs102);
  sf_mex_destroy(&c6_lhs102);
  sf_mex_destroy(&c6_rhs103);
  sf_mex_destroy(&c6_lhs103);
  sf_mex_destroy(&c6_rhs104);
  sf_mex_destroy(&c6_lhs104);
  sf_mex_destroy(&c6_rhs105);
  sf_mex_destroy(&c6_lhs105);
  sf_mex_destroy(&c6_rhs106);
  sf_mex_destroy(&c6_lhs106);
  sf_mex_destroy(&c6_rhs107);
  sf_mex_destroy(&c6_lhs107);
  sf_mex_destroy(&c6_rhs108);
  sf_mex_destroy(&c6_lhs108);
  sf_mex_destroy(&c6_rhs109);
  sf_mex_destroy(&c6_lhs109);
  sf_mex_destroy(&c6_rhs110);
  sf_mex_destroy(&c6_lhs110);
  sf_mex_destroy(&c6_rhs111);
  sf_mex_destroy(&c6_lhs111);
  sf_mex_destroy(&c6_rhs112);
  sf_mex_destroy(&c6_lhs112);
  sf_mex_destroy(&c6_rhs113);
  sf_mex_destroy(&c6_lhs113);
  sf_mex_destroy(&c6_rhs114);
  sf_mex_destroy(&c6_lhs114);
  sf_mex_destroy(&c6_rhs115);
  sf_mex_destroy(&c6_lhs115);
  sf_mex_destroy(&c6_rhs116);
  sf_mex_destroy(&c6_lhs116);
  sf_mex_destroy(&c6_rhs117);
  sf_mex_destroy(&c6_lhs117);
  sf_mex_destroy(&c6_rhs118);
  sf_mex_destroy(&c6_lhs118);
  sf_mex_destroy(&c6_rhs119);
  sf_mex_destroy(&c6_lhs119);
  sf_mex_destroy(&c6_rhs120);
  sf_mex_destroy(&c6_lhs120);
  sf_mex_destroy(&c6_rhs121);
  sf_mex_destroy(&c6_lhs121);
  sf_mex_destroy(&c6_rhs122);
  sf_mex_destroy(&c6_lhs122);
  sf_mex_destroy(&c6_rhs123);
  sf_mex_destroy(&c6_lhs123);
  sf_mex_destroy(&c6_rhs124);
  sf_mex_destroy(&c6_lhs124);
  sf_mex_destroy(&c6_rhs125);
  sf_mex_destroy(&c6_lhs125);
  sf_mex_destroy(&c6_rhs126);
  sf_mex_destroy(&c6_lhs126);
  sf_mex_destroy(&c6_rhs127);
  sf_mex_destroy(&c6_lhs127);
}

static void c6_c_info_helper(const mxArray **c6_info)
{
  const mxArray *c6_rhs128 = NULL;
  const mxArray *c6_lhs128 = NULL;
  const mxArray *c6_rhs129 = NULL;
  const mxArray *c6_lhs129 = NULL;
  const mxArray *c6_rhs130 = NULL;
  const mxArray *c6_lhs130 = NULL;
  const mxArray *c6_rhs131 = NULL;
  const mxArray *c6_lhs131 = NULL;
  const mxArray *c6_rhs132 = NULL;
  const mxArray *c6_lhs132 = NULL;
  const mxArray *c6_rhs133 = NULL;
  const mxArray *c6_lhs133 = NULL;
  const mxArray *c6_rhs134 = NULL;
  const mxArray *c6_lhs134 = NULL;
  const mxArray *c6_rhs135 = NULL;
  const mxArray *c6_lhs135 = NULL;
  const mxArray *c6_rhs136 = NULL;
  const mxArray *c6_lhs136 = NULL;
  const mxArray *c6_rhs137 = NULL;
  const mxArray *c6_lhs137 = NULL;
  const mxArray *c6_rhs138 = NULL;
  const mxArray *c6_lhs138 = NULL;
  const mxArray *c6_rhs139 = NULL;
  const mxArray *c6_lhs139 = NULL;
  const mxArray *c6_rhs140 = NULL;
  const mxArray *c6_lhs140 = NULL;
  const mxArray *c6_rhs141 = NULL;
  const mxArray *c6_lhs141 = NULL;
  const mxArray *c6_rhs142 = NULL;
  const mxArray *c6_lhs142 = NULL;
  const mxArray *c6_rhs143 = NULL;
  const mxArray *c6_lhs143 = NULL;
  const mxArray *c6_rhs144 = NULL;
  const mxArray *c6_lhs144 = NULL;
  const mxArray *c6_rhs145 = NULL;
  const mxArray *c6_lhs145 = NULL;
  const mxArray *c6_rhs146 = NULL;
  const mxArray *c6_lhs146 = NULL;
  const mxArray *c6_rhs147 = NULL;
  const mxArray *c6_lhs147 = NULL;
  const mxArray *c6_rhs148 = NULL;
  const mxArray *c6_lhs148 = NULL;
  const mxArray *c6_rhs149 = NULL;
  const mxArray *c6_lhs149 = NULL;
  const mxArray *c6_rhs150 = NULL;
  const mxArray *c6_lhs150 = NULL;
  const mxArray *c6_rhs151 = NULL;
  const mxArray *c6_lhs151 = NULL;
  const mxArray *c6_rhs152 = NULL;
  const mxArray *c6_lhs152 = NULL;
  const mxArray *c6_rhs153 = NULL;
  const mxArray *c6_lhs153 = NULL;
  const mxArray *c6_rhs154 = NULL;
  const mxArray *c6_lhs154 = NULL;
  const mxArray *c6_rhs155 = NULL;
  const mxArray *c6_lhs155 = NULL;
  const mxArray *c6_rhs156 = NULL;
  const mxArray *c6_lhs156 = NULL;
  const mxArray *c6_rhs157 = NULL;
  const mxArray *c6_lhs157 = NULL;
  const mxArray *c6_rhs158 = NULL;
  const mxArray *c6_lhs158 = NULL;
  const mxArray *c6_rhs159 = NULL;
  const mxArray *c6_lhs159 = NULL;
  const mxArray *c6_rhs160 = NULL;
  const mxArray *c6_lhs160 = NULL;
  const mxArray *c6_rhs161 = NULL;
  const mxArray *c6_lhs161 = NULL;
  const mxArray *c6_rhs162 = NULL;
  const mxArray *c6_lhs162 = NULL;
  const mxArray *c6_rhs163 = NULL;
  const mxArray *c6_lhs163 = NULL;
  const mxArray *c6_rhs164 = NULL;
  const mxArray *c6_lhs164 = NULL;
  const mxArray *c6_rhs165 = NULL;
  const mxArray *c6_lhs165 = NULL;
  const mxArray *c6_rhs166 = NULL;
  const mxArray *c6_lhs166 = NULL;
  const mxArray *c6_rhs167 = NULL;
  const mxArray *c6_lhs167 = NULL;
  const mxArray *c6_rhs168 = NULL;
  const mxArray *c6_lhs168 = NULL;
  const mxArray *c6_rhs169 = NULL;
  const mxArray *c6_lhs169 = NULL;
  const mxArray *c6_rhs170 = NULL;
  const mxArray *c6_lhs170 = NULL;
  const mxArray *c6_rhs171 = NULL;
  const mxArray *c6_lhs171 = NULL;
  const mxArray *c6_rhs172 = NULL;
  const mxArray *c6_lhs172 = NULL;
  const mxArray *c6_rhs173 = NULL;
  const mxArray *c6_lhs173 = NULL;
  const mxArray *c6_rhs174 = NULL;
  const mxArray *c6_lhs174 = NULL;
  const mxArray *c6_rhs175 = NULL;
  const mxArray *c6_lhs175 = NULL;
  const mxArray *c6_rhs176 = NULL;
  const mxArray *c6_lhs176 = NULL;
  const mxArray *c6_rhs177 = NULL;
  const mxArray *c6_lhs177 = NULL;
  const mxArray *c6_rhs178 = NULL;
  const mxArray *c6_lhs178 = NULL;
  const mxArray *c6_rhs179 = NULL;
  const mxArray *c6_lhs179 = NULL;
  const mxArray *c6_rhs180 = NULL;
  const mxArray *c6_lhs180 = NULL;
  const mxArray *c6_rhs181 = NULL;
  const mxArray *c6_lhs181 = NULL;
  const mxArray *c6_rhs182 = NULL;
  const mxArray *c6_lhs182 = NULL;
  const mxArray *c6_rhs183 = NULL;
  const mxArray *c6_lhs183 = NULL;
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 128);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("abs"), "name", "name", 128);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 128);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 128);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c6_rhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs128), "rhs", "rhs",
                  128);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs128), "lhs", "lhs",
                  128);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 129);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 129);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 129);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 129);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c6_rhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs129), "rhs", "rhs",
                  129);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs129), "lhs", "lhs",
                  129);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 130);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 130);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 130);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 130);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c6_rhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs130), "rhs", "rhs",
                  130);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs130), "lhs", "lhs",
                  130);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 131);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 131);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 131);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c6_rhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs131), "rhs", "rhs",
                  131);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs131), "lhs", "lhs",
                  131);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 132);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 132);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 132);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 132);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c6_rhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs132), "rhs", "rhs",
                  132);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs132), "lhs", "lhs",
                  132);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 133);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_div"), "name", "name", 133);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 133);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 133);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c6_rhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs133), "rhs", "rhs",
                  133);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs133), "lhs", "lhs",
                  133);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 134);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 134);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 134);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c6_rhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs134), "rhs", "rhs",
                  134);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs134), "lhs", "lhs",
                  134);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 135);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_xgeru"), "name", "name",
                  135);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"),
                  "resolved", "resolved", 135);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c6_rhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs135), "rhs", "rhs",
                  135);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs135), "lhs", "lhs",
                  135);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 136);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 136);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 136);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c6_rhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs136), "rhs", "rhs",
                  136);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs136), "lhs", "lhs",
                  136);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 137);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.xgeru"),
                  "name", "name", 137);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgeru.p"),
                  "resolved", "resolved", 137);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c6_rhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs137), "rhs", "rhs",
                  137);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs137), "lhs", "lhs",
                  137);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgeru.p"),
                  "context", "context", 138);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.xger"),
                  "name", "name", 138);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "resolved", "resolved", 138);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c6_rhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs138), "rhs", "rhs",
                  138);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs138), "lhs", "lhs",
                  138);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "context", "context", 139);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 139);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 139);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c6_rhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs139), "rhs", "rhs",
                  139);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs139), "lhs", "lhs",
                  139);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 140);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 140);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 140);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c6_rhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs140), "rhs", "rhs",
                  140);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs140), "lhs", "lhs",
                  140);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 141);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.int"),
                  "name", "name", 141);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 141);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/int.p"),
                  "resolved", "resolved", 141);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c6_rhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs141), "rhs", "rhs",
                  141);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs141), "lhs", "lhs",
                  141);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 142);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmax"), "name", "name", 142);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 142);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c6_rhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs142), "rhs", "rhs",
                  142);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs142), "lhs", "lhs",
                  142);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 143);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("min"), "name", "name", 143);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 143);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 143);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1311262518U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c6_rhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs143), "rhs", "rhs",
                  143);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs143), "lhs", "lhs",
                  143);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 144);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 144);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 144);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c6_rhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs144), "rhs", "rhs",
                  144);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs144), "lhs", "lhs",
                  144);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 145);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 145);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 145);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c6_rhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs145), "rhs", "rhs",
                  145);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs145), "lhs", "lhs",
                  145);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 146);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 146);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 146);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c6_rhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs146), "rhs", "rhs",
                  146);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs146), "lhs", "lhs",
                  146);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 147);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 147);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 147);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c6_rhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs147), "rhs", "rhs",
                  147);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs147), "lhs", "lhs",
                  147);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "context", "context", 148);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.refblas.xger"),
                  "name", "name", 148);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xger.p"),
                  "resolved", "resolved", 148);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c6_rhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs148), "rhs", "rhs",
                  148);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs148), "lhs", "lhs",
                  148);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xger.p"),
                  "context", "context", 149);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.refblas.xgerx"),
                  "name", "name", 149);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "resolved", "resolved", 149);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c6_rhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs149), "rhs", "rhs",
                  149);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs149), "lhs", "lhs",
                  149);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 150);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("abs"), "name", "name", 150);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 150);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 150);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c6_rhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs150), "rhs", "rhs",
                  150);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs150), "lhs", "lhs",
                  150);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 151);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 151);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 151);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c6_rhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs151), "rhs", "rhs",
                  151);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs151), "lhs", "lhs",
                  151);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 152);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 152);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 152);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c6_rhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs152), "rhs", "rhs",
                  152);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs152), "lhs", "lhs",
                  152);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 153);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 153);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 153);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c6_rhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs153), "rhs", "rhs",
                  153);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs153), "lhs", "lhs",
                  153);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 154);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 154);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 154);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 154);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c6_rhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs154), "rhs", "rhs",
                  154);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs154), "lhs", "lhs",
                  154);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular"),
                  "context", "context", 155);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_warning"), "name", "name",
                  155);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 155);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c6_rhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs155), "rhs", "rhs",
                  155);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs155), "lhs", "lhs",
                  155);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 156);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 156);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 156);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c6_rhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs156), "rhs", "rhs",
                  156);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs156), "lhs", "lhs",
                  156);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 157);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 157);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 157);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c6_rhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs157), "rhs", "rhs",
                  157);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs157), "lhs", "lhs",
                  157);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 158);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_xtrsm"), "name", "name",
                  158);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 158);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"),
                  "resolved", "resolved", 158);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c6_rhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs158), "rhs", "rhs",
                  158);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs158), "lhs", "lhs",
                  158);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 159);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 159);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 159);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 159);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c6_rhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs159), "rhs", "rhs",
                  159);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs159), "lhs", "lhs",
                  159);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 160);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.xtrsm"),
                  "name", "name", 160);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 160);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "resolved", "resolved", 160);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c6_rhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs160), "rhs", "rhs",
                  160);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs160), "lhs", "lhs",
                  160);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 161);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 161);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 161);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 161);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c6_rhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs161), "rhs", "rhs",
                  161);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs161), "lhs", "lhs",
                  161);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p!below_threshold"),
                  "context", "context", 162);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 162);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 162);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 162);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c6_rhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs162), "rhs", "rhs",
                  162);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs162), "lhs", "lhs",
                  162);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 163);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 163);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 163);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c6_rhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs163), "rhs", "rhs",
                  163);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs163), "lhs", "lhs",
                  163);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 164);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.refblas.xtrsm"),
                  "name", "name", 164);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 164);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "resolved", "resolved", 164);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 164);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 164);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 164);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 164);
  sf_mex_assign(&c6_rhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs164), "rhs", "rhs",
                  164);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs164), "lhs", "lhs",
                  164);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 165);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 165);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 165);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 165);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 165);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 165);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 165);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 165);
  sf_mex_assign(&c6_rhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs165), "rhs", "rhs",
                  165);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs165), "lhs", "lhs",
                  165);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 166);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 166);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 166);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 166);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 166);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 166);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 166);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 166);
  sf_mex_assign(&c6_rhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs166), "rhs", "rhs",
                  166);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs166), "lhs", "lhs",
                  166);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 167);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 167);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 167);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 167);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 167);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 167);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 167);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 167);
  sf_mex_assign(&c6_rhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs167), "rhs", "rhs",
                  167);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs167), "lhs", "lhs",
                  167);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 168);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("intmin"), "name", "name", 168);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 168);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 168);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 168);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 168);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 168);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 168);
  sf_mex_assign(&c6_rhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs168), "rhs", "rhs",
                  168);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs168), "lhs", "lhs",
                  168);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 169);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("rdivide"), "name", "name", 169);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 169);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 169);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 169);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 169);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 169);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 169);
  sf_mex_assign(&c6_rhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs169), "rhs", "rhs",
                  169);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs169), "lhs", "lhs",
                  169);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 170);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 170);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 170);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 170);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 170);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 170);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 170);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 170);
  sf_mex_assign(&c6_rhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs170), "rhs", "rhs",
                  170);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs170), "lhs", "lhs",
                  170);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 171);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 171);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 171);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 171);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 171);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 171);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 171);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 171);
  sf_mex_assign(&c6_rhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs171), "rhs", "rhs",
                  171);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs171), "lhs", "lhs",
                  171);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 172);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_div"), "name", "name", 172);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 172);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 172);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 172);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 172);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 172);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 172);
  sf_mex_assign(&c6_rhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs172), "rhs", "rhs",
                  172);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs172), "lhs", "lhs",
                  172);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 173);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("log2"), "name", "name", 173);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 173);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log2.m"), "resolved",
                  "resolved", 173);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1343837582U), "fileTimeLo",
                  "fileTimeLo", 173);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 173);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 173);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 173);
  sf_mex_assign(&c6_rhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs173), "rhs", "rhs",
                  173);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs173), "lhs", "lhs",
                  173);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log2.m!scalar_frexp"),
                  "context", "context", 174);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("isfinite"), "name", "name",
                  174);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 174);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "resolved",
                  "resolved", 174);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 174);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 174);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 174);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 174);
  sf_mex_assign(&c6_rhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs174), "rhs", "rhs",
                  174);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs174), "lhs", "lhs",
                  174);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 175);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 175);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 175);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 175);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 175);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 175);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 175);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 175);
  sf_mex_assign(&c6_rhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs175), "rhs", "rhs",
                  175);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs175), "lhs", "lhs",
                  175);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 176);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("isinf"), "name", "name", 176);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 176);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 176);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 176);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 176);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 176);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 176);
  sf_mex_assign(&c6_rhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs176), "rhs", "rhs",
                  176);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs176), "lhs", "lhs",
                  176);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 177);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("isnan"), "name", "name", 177);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 177);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 177);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 177);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 177);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 177);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 177);
  sf_mex_assign(&c6_rhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs177), "rhs", "rhs",
                  177);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs177), "lhs", "lhs",
                  177);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 178);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("pow2"), "name", "name", 178);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 178);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/pow2.m"), "resolved",
                  "resolved", 178);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1343837582U), "fileTimeLo",
                  "fileTimeLo", 178);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 178);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 178);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 178);
  sf_mex_assign(&c6_rhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs178), "rhs", "rhs",
                  178);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs178), "lhs", "lhs",
                  178);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/pow2.m"), "context",
                  "context", 179);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_scalar_pow2"), "name",
                  "name", 179);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 179);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_pow2.m"),
                  "resolved", "resolved", 179);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1286825932U), "fileTimeLo",
                  "fileTimeLo", 179);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 179);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 179);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 179);
  sf_mex_assign(&c6_rhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs179), "rhs", "rhs",
                  179);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs179), "lhs", "lhs",
                  179);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_pow2.m"),
                  "context", "context", 180);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("power"), "name", "name", 180);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 180);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 180);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 180);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 180);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 180);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 180);
  sf_mex_assign(&c6_rhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs180), "rhs", "rhs",
                  180);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs180), "lhs", "lhs",
                  180);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 181);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_error"), "name", "name",
                  181);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 181);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 181);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 181);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 181);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 181);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 181);
  sf_mex_assign(&c6_rhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs181), "rhs", "rhs",
                  181);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs181), "lhs", "lhs",
                  181);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 182);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_div"), "name", "name", 182);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 182);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 182);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 182);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 182);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 182);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 182);
  sf_mex_assign(&c6_rhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs182), "rhs", "rhs",
                  182);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs182), "lhs", "lhs",
                  182);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 183);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 183);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 183);
  sf_mex_addfield(*c6_info, c6_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 183);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 183);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 183);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 183);
  sf_mex_addfield(*c6_info, c6_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 183);
  sf_mex_assign(&c6_rhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c6_lhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_rhs183), "rhs", "rhs",
                  183);
  sf_mex_addfield(*c6_info, sf_mex_duplicatearraysafe(&c6_lhs183), "lhs", "lhs",
                  183);
  sf_mex_destroy(&c6_rhs128);
  sf_mex_destroy(&c6_lhs128);
  sf_mex_destroy(&c6_rhs129);
  sf_mex_destroy(&c6_lhs129);
  sf_mex_destroy(&c6_rhs130);
  sf_mex_destroy(&c6_lhs130);
  sf_mex_destroy(&c6_rhs131);
  sf_mex_destroy(&c6_lhs131);
  sf_mex_destroy(&c6_rhs132);
  sf_mex_destroy(&c6_lhs132);
  sf_mex_destroy(&c6_rhs133);
  sf_mex_destroy(&c6_lhs133);
  sf_mex_destroy(&c6_rhs134);
  sf_mex_destroy(&c6_lhs134);
  sf_mex_destroy(&c6_rhs135);
  sf_mex_destroy(&c6_lhs135);
  sf_mex_destroy(&c6_rhs136);
  sf_mex_destroy(&c6_lhs136);
  sf_mex_destroy(&c6_rhs137);
  sf_mex_destroy(&c6_lhs137);
  sf_mex_destroy(&c6_rhs138);
  sf_mex_destroy(&c6_lhs138);
  sf_mex_destroy(&c6_rhs139);
  sf_mex_destroy(&c6_lhs139);
  sf_mex_destroy(&c6_rhs140);
  sf_mex_destroy(&c6_lhs140);
  sf_mex_destroy(&c6_rhs141);
  sf_mex_destroy(&c6_lhs141);
  sf_mex_destroy(&c6_rhs142);
  sf_mex_destroy(&c6_lhs142);
  sf_mex_destroy(&c6_rhs143);
  sf_mex_destroy(&c6_lhs143);
  sf_mex_destroy(&c6_rhs144);
  sf_mex_destroy(&c6_lhs144);
  sf_mex_destroy(&c6_rhs145);
  sf_mex_destroy(&c6_lhs145);
  sf_mex_destroy(&c6_rhs146);
  sf_mex_destroy(&c6_lhs146);
  sf_mex_destroy(&c6_rhs147);
  sf_mex_destroy(&c6_lhs147);
  sf_mex_destroy(&c6_rhs148);
  sf_mex_destroy(&c6_lhs148);
  sf_mex_destroy(&c6_rhs149);
  sf_mex_destroy(&c6_lhs149);
  sf_mex_destroy(&c6_rhs150);
  sf_mex_destroy(&c6_lhs150);
  sf_mex_destroy(&c6_rhs151);
  sf_mex_destroy(&c6_lhs151);
  sf_mex_destroy(&c6_rhs152);
  sf_mex_destroy(&c6_lhs152);
  sf_mex_destroy(&c6_rhs153);
  sf_mex_destroy(&c6_lhs153);
  sf_mex_destroy(&c6_rhs154);
  sf_mex_destroy(&c6_lhs154);
  sf_mex_destroy(&c6_rhs155);
  sf_mex_destroy(&c6_lhs155);
  sf_mex_destroy(&c6_rhs156);
  sf_mex_destroy(&c6_lhs156);
  sf_mex_destroy(&c6_rhs157);
  sf_mex_destroy(&c6_lhs157);
  sf_mex_destroy(&c6_rhs158);
  sf_mex_destroy(&c6_lhs158);
  sf_mex_destroy(&c6_rhs159);
  sf_mex_destroy(&c6_lhs159);
  sf_mex_destroy(&c6_rhs160);
  sf_mex_destroy(&c6_lhs160);
  sf_mex_destroy(&c6_rhs161);
  sf_mex_destroy(&c6_lhs161);
  sf_mex_destroy(&c6_rhs162);
  sf_mex_destroy(&c6_lhs162);
  sf_mex_destroy(&c6_rhs163);
  sf_mex_destroy(&c6_lhs163);
  sf_mex_destroy(&c6_rhs164);
  sf_mex_destroy(&c6_lhs164);
  sf_mex_destroy(&c6_rhs165);
  sf_mex_destroy(&c6_lhs165);
  sf_mex_destroy(&c6_rhs166);
  sf_mex_destroy(&c6_lhs166);
  sf_mex_destroy(&c6_rhs167);
  sf_mex_destroy(&c6_lhs167);
  sf_mex_destroy(&c6_rhs168);
  sf_mex_destroy(&c6_lhs168);
  sf_mex_destroy(&c6_rhs169);
  sf_mex_destroy(&c6_lhs169);
  sf_mex_destroy(&c6_rhs170);
  sf_mex_destroy(&c6_lhs170);
  sf_mex_destroy(&c6_rhs171);
  sf_mex_destroy(&c6_lhs171);
  sf_mex_destroy(&c6_rhs172);
  sf_mex_destroy(&c6_lhs172);
  sf_mex_destroy(&c6_rhs173);
  sf_mex_destroy(&c6_lhs173);
  sf_mex_destroy(&c6_rhs174);
  sf_mex_destroy(&c6_lhs174);
  sf_mex_destroy(&c6_rhs175);
  sf_mex_destroy(&c6_lhs175);
  sf_mex_destroy(&c6_rhs176);
  sf_mex_destroy(&c6_lhs176);
  sf_mex_destroy(&c6_rhs177);
  sf_mex_destroy(&c6_lhs177);
  sf_mex_destroy(&c6_rhs178);
  sf_mex_destroy(&c6_lhs178);
  sf_mex_destroy(&c6_rhs179);
  sf_mex_destroy(&c6_lhs179);
  sf_mex_destroy(&c6_rhs180);
  sf_mex_destroy(&c6_lhs180);
  sf_mex_destroy(&c6_rhs181);
  sf_mex_destroy(&c6_lhs181);
  sf_mex_destroy(&c6_rhs182);
  sf_mex_destroy(&c6_lhs182);
  sf_mex_destroy(&c6_rhs183);
  sf_mex_destroy(&c6_lhs183);
}

static void c6_eml_scalar_eg(SFc6_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c6_eml_switch_helper(SFc6_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c6_eye(SFc6_Model_02InstanceStruct *chartInstance, real_T c6_I[9])
{
  int32_T c6_i115;
  int32_T c6_k;
  int32_T c6_b_k;
  (void)chartInstance;
  for (c6_i115 = 0; c6_i115 < 9; c6_i115++) {
    c6_I[c6_i115] = 0.0;
  }

  for (c6_k = 1; c6_k < 4; c6_k++) {
    c6_b_k = c6_k - 1;
    c6_I[c6_b_k + 3 * c6_b_k] = 1.0;
  }
}

static void c6_expm(SFc6_Model_02InstanceStruct *chartInstance, real_T c6_A[36],
                    real_T c6_F[36])
{
  real_T c6_normA;
  int32_T c6_j;
  real_T c6_b_j;
  real_T c6_s;
  int32_T c6_i;
  real_T c6_b_i;
  real_T c6_x;
  real_T c6_b_x;
  real_T c6_y;
  real_T c6_c_x;
  boolean_T c6_b;
  int32_T c6_c_i;
  real_T c6_d_i;
  static real_T c6_theta[5] = { 0.01495585217958292, 0.253939833006323,
    0.95041789961629319, 2.097847961257068, 5.3719203511481517 };

  int32_T c6_i116;
  real_T c6_b_A[36];
  static real_T c6_dv8[5] = { 3.0, 5.0, 7.0, 9.0, 13.0 };

  real_T c6_d_x;
  real_T c6_e_x;
  real_T c6_f_x;
  real_T c6_g_x;
  boolean_T c6_b_b;
  boolean_T c6_b0;
  real_T c6_h_x;
  boolean_T c6_c_b;
  boolean_T c6_b1;
  boolean_T c6_d_b;
  int32_T c6_eint;
  real_T c6_fdbl;
  int32_T c6_b_eint;
  real_T c6_b_fdbl;
  int32_T c6_c_eint;
  real_T c6_d1;
  real_T c6_d2;
  real_T c6_t;
  real_T c6_b_s;
  real_T c6_b_t;
  real_T c6_c_s;
  real_T c6_a;
  real_T c6_b_a;
  real_T c6_e_b;
  real_T c6_f_b;
  real_T c6_bk;
  real_T c6_g_b;
  real_T c6_br;
  real_T c6_b_y;
  real_T c6_c_y;
  real_T c6_d_y;
  int32_T c6_i117;
  int32_T c6_i118;
  real_T c6_c_A[36];
  real_T c6_d_s;
  int32_T c6_i119;
  int32_T c6_c_j;
  int32_T c6_i120;
  real_T c6_c_a[36];
  int32_T c6_i121;
  int32_T c6_i122;
  real_T c6_d_a[36];
  int32_T c6_i123;
  real_T c6_e_a[36];
  boolean_T exitg1;
  boolean_T exitg2;
  c6_normA = 0.0;
  c6_j = 0;
  exitg2 = false;
  while ((exitg2 == false) && (c6_j < 6)) {
    c6_b_j = 1.0 + (real_T)c6_j;
    c6_s = 0.0;
    for (c6_i = 0; c6_i < 6; c6_i++) {
      c6_b_i = 1.0 + (real_T)c6_i;
      c6_x = c6_A[((int32_T)c6_b_i + 6 * ((int32_T)c6_b_j - 1)) - 1];
      c6_b_x = c6_x;
      c6_y = muDoubleScalarAbs(c6_b_x);
      c6_s += c6_y;
    }

    c6_c_x = c6_s;
    c6_b = muDoubleScalarIsNaN(c6_c_x);
    if (c6_b) {
      c6_normA = rtNaN;
      exitg2 = true;
    } else {
      if (c6_s > c6_normA) {
        c6_normA = c6_s;
      }

      c6_j++;
    }
  }

  if (c6_normA <= 5.3719203511481517) {
    c6_c_i = 0;
    exitg1 = false;
    while ((exitg1 == false) && (c6_c_i < 5)) {
      c6_d_i = 1.0 + (real_T)c6_c_i;
      if (c6_normA <= c6_theta[(int32_T)c6_d_i - 1]) {
        for (c6_i116 = 0; c6_i116 < 36; c6_i116++) {
          c6_b_A[c6_i116] = c6_A[c6_i116];
        }

        c6_PadeApproximantOfDegree(chartInstance, c6_b_A, c6_dv8[(int32_T)c6_d_i
          - 1], c6_F);
        exitg1 = true;
      } else {
        c6_c_i++;
      }
    }
  } else {
    c6_d_x = c6_normA / 5.3719203511481517;
    c6_e_x = c6_d_x;
    c6_f_x = c6_e_x;
    c6_g_x = c6_f_x;
    c6_b_b = muDoubleScalarIsInf(c6_g_x);
    c6_b0 = !c6_b_b;
    c6_h_x = c6_f_x;
    c6_c_b = muDoubleScalarIsNaN(c6_h_x);
    c6_b1 = !c6_c_b;
    c6_d_b = (c6_b0 && c6_b1);
    if (c6_d_b) {
      c6_fdbl = frexp(c6_e_x, &c6_eint);
      c6_b_eint = c6_eint;
      c6_b_fdbl = c6_fdbl;
      c6_c_eint = c6_b_eint;
      c6_d1 = c6_b_fdbl;
      c6_d2 = (real_T)c6_c_eint;
    } else {
      c6_d1 = c6_e_x;
      c6_d2 = 0.0;
    }

    c6_t = c6_d1;
    c6_b_s = c6_d2;
    c6_b_t = c6_t;
    c6_c_s = c6_b_s;
    if (c6_b_t == 0.5) {
      c6_c_s--;
    }

    c6_a = c6_c_s;
    c6_b_a = c6_a;
    c6_e_b = c6_b_a;
    c6_f_b = c6_e_b;
    c6_eml_scalar_eg(chartInstance);
    c6_bk = c6_f_b;
    c6_g_b = c6_bk;
    c6_eml_scalar_eg(chartInstance);
    c6_br = c6_g_b;
    c6_b_y = muDoubleScalarPower(2.0, c6_br);
    c6_c_y = c6_b_y;
    c6_d_y = c6_c_y;
    for (c6_i117 = 0; c6_i117 < 36; c6_i117++) {
      c6_A[c6_i117] /= c6_d_y;
    }

    for (c6_i118 = 0; c6_i118 < 36; c6_i118++) {
      c6_c_A[c6_i118] = c6_A[c6_i118];
    }

    c6_PadeApproximantOfDegree(chartInstance, c6_c_A, 13.0, c6_F);
    c6_d_s = c6_c_s;
    c6_i119 = (int32_T)c6_d_s;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c6_d_s, mxDOUBLE_CLASS, c6_i119);
    for (c6_c_j = 0; c6_c_j < c6_i119; c6_c_j++) {
      for (c6_i120 = 0; c6_i120 < 36; c6_i120++) {
        c6_c_a[c6_i120] = c6_F[c6_i120];
      }

      c6_b_eml_scalar_eg(chartInstance);
      c6_b_eml_scalar_eg(chartInstance);
      for (c6_i121 = 0; c6_i121 < 36; c6_i121++) {
        c6_F[c6_i121] = 0.0;
      }

      for (c6_i122 = 0; c6_i122 < 36; c6_i122++) {
        c6_d_a[c6_i122] = c6_c_a[c6_i122];
      }

      for (c6_i123 = 0; c6_i123 < 36; c6_i123++) {
        c6_e_a[c6_i123] = c6_c_a[c6_i123];
      }

      c6_b_eml_xgemm(chartInstance, c6_d_a, c6_e_a, c6_F);
    }
  }
}

static void c6_PadeApproximantOfDegree(SFc6_Model_02InstanceStruct
  *chartInstance, real_T c6_A[36], real_T c6_m, real_T c6_F[36])
{
  int32_T c6_i124;
  real_T c6_A2[36];
  int32_T c6_i125;
  real_T c6_b_A[36];
  int32_T c6_i126;
  real_T c6_c_A[36];
  int32_T c6_i127;
  real_T c6_U[36];
  int32_T c6_k;
  real_T c6_b_k;
  int32_T c6_i128;
  real_T c6_y[36];
  int32_T c6_i129;
  int32_T c6_i130;
  real_T c6_d_A[36];
  int32_T c6_i131;
  real_T c6_b_y[36];
  int32_T c6_i132;
  real_T c6_d;
  int32_T c6_i133;
  real_T c6_A3[36];
  int32_T c6_i134;
  real_T c6_b_A2[36];
  int32_T c6_i135;
  real_T c6_c_A2[36];
  int32_T c6_i136;
  int32_T c6_i137;
  int32_T c6_c_k;
  int32_T c6_i138;
  int32_T c6_i139;
  int32_T c6_i140;
  real_T c6_e_A[36];
  int32_T c6_i141;
  real_T c6_c_y[36];
  int32_T c6_i142;
  int32_T c6_i143;
  int32_T c6_i144;
  int32_T c6_i145;
  real_T c6_A4[36];
  int32_T c6_i146;
  real_T c6_b_A3[36];
  int32_T c6_i147;
  real_T c6_d_A2[36];
  int32_T c6_i148;
  int32_T c6_i149;
  real_T c6_d_y[36];
  int32_T c6_i150;
  int32_T c6_d_k;
  int32_T c6_i151;
  int32_T c6_i152;
  int32_T c6_i153;
  real_T c6_f_A[36];
  int32_T c6_i154;
  real_T c6_e_y[36];
  int32_T c6_i155;
  int32_T c6_i156;
  int32_T c6_i157;
  int32_T c6_i158;
  int32_T c6_i159;
  int32_T c6_i160;
  real_T c6_b_A4[36];
  int32_T c6_i161;
  real_T c6_e_A2[36];
  int32_T c6_i162;
  int32_T c6_i163;
  int32_T c6_i164;
  real_T c6_f_y[36];
  int32_T c6_i165;
  int32_T c6_e_k;
  int32_T c6_i166;
  int32_T c6_i167;
  int32_T c6_i168;
  real_T c6_g_A[36];
  int32_T c6_i169;
  real_T c6_g_y[36];
  int32_T c6_i170;
  int32_T c6_i171;
  int32_T c6_i172;
  int32_T c6_i173;
  int32_T c6_i174;
  int32_T c6_i175;
  int32_T c6_i176;
  int32_T c6_i177;
  int32_T c6_i178;
  int32_T c6_f_k;
  int32_T c6_i179;
  int32_T c6_i180;
  int32_T c6_i181;
  int32_T c6_i182;
  int32_T c6_i183;
  real_T c6_c_A4[36];
  int32_T c6_i184;
  real_T c6_h_y[36];
  int32_T c6_i185;
  int32_T c6_i186;
  int32_T c6_i187;
  real_T c6_h_A[36];
  int32_T c6_i188;
  real_T c6_i_y[36];
  int32_T c6_i189;
  int32_T c6_i190;
  int32_T c6_i191;
  int32_T c6_i192;
  int32_T c6_i193;
  int32_T c6_i194;
  real_T c6_d_A4[36];
  int32_T c6_i195;
  real_T c6_j_y[36];
  int32_T c6_i196;
  int32_T c6_i197;
  int32_T c6_i198;
  int32_T c6_i199;
  int32_T c6_g_k;
  int32_T c6_h_k;
  real_T c6_uk;
  int32_T c6_info;
  int32_T c6_ipiv[6];
  int32_T c6_b_info;
  int32_T c6_c_info;
  int32_T c6_d_info;
  int32_T c6_xi;
  int32_T c6_b_xi;
  int32_T c6_ip;
  int32_T c6_xj;
  int32_T c6_b_xj;
  real_T c6_temp;
  int32_T c6_i200;
  real_T c6_b_U[36];
  int32_T c6_i201;
  real_T c6_c_U[36];
  c6_b_eml_scalar_eg(chartInstance);
  c6_b_eml_scalar_eg(chartInstance);
  for (c6_i124 = 0; c6_i124 < 36; c6_i124++) {
    c6_A2[c6_i124] = 0.0;
  }

  for (c6_i125 = 0; c6_i125 < 36; c6_i125++) {
    c6_b_A[c6_i125] = c6_A[c6_i125];
  }

  for (c6_i126 = 0; c6_i126 < 36; c6_i126++) {
    c6_c_A[c6_i126] = c6_A[c6_i126];
  }

  c6_b_eml_xgemm(chartInstance, c6_b_A, c6_c_A, c6_A2);
  if (c6_m == 3.0) {
    for (c6_i127 = 0; c6_i127 < 36; c6_i127++) {
      c6_U[c6_i127] = c6_A2[c6_i127];
    }

    for (c6_k = 0; c6_k < 6; c6_k++) {
      c6_b_k = 1.0 + (real_T)c6_k;
      c6_U[((int32_T)c6_b_k + 6 * ((int32_T)c6_b_k - 1)) - 1] += 60.0;
    }

    for (c6_i128 = 0; c6_i128 < 36; c6_i128++) {
      c6_y[c6_i128] = c6_U[c6_i128];
    }

    c6_b_eml_scalar_eg(chartInstance);
    c6_b_eml_scalar_eg(chartInstance);
    for (c6_i129 = 0; c6_i129 < 36; c6_i129++) {
      c6_U[c6_i129] = 0.0;
    }

    for (c6_i130 = 0; c6_i130 < 36; c6_i130++) {
      c6_d_A[c6_i130] = c6_A[c6_i130];
    }

    for (c6_i131 = 0; c6_i131 < 36; c6_i131++) {
      c6_b_y[c6_i131] = c6_y[c6_i131];
    }

    c6_b_eml_xgemm(chartInstance, c6_d_A, c6_b_y, c6_U);
    for (c6_i132 = 0; c6_i132 < 36; c6_i132++) {
      c6_F[c6_i132] = 12.0 * c6_A2[c6_i132];
    }

    c6_d = 120.0;
  } else {
    c6_b_eml_scalar_eg(chartInstance);
    c6_b_eml_scalar_eg(chartInstance);
    for (c6_i133 = 0; c6_i133 < 36; c6_i133++) {
      c6_A3[c6_i133] = 0.0;
    }

    for (c6_i134 = 0; c6_i134 < 36; c6_i134++) {
      c6_b_A2[c6_i134] = c6_A2[c6_i134];
    }

    for (c6_i135 = 0; c6_i135 < 36; c6_i135++) {
      c6_c_A2[c6_i135] = c6_A2[c6_i135];
    }

    c6_b_eml_xgemm(chartInstance, c6_b_A2, c6_c_A2, c6_A3);
    if (c6_m == 5.0) {
      for (c6_i136 = 0; c6_i136 < 36; c6_i136++) {
        c6_U[c6_i136] = 420.0 * c6_A2[c6_i136];
      }

      for (c6_i137 = 0; c6_i137 < 36; c6_i137++) {
        c6_U[c6_i137] += c6_A3[c6_i137];
      }

      for (c6_c_k = 0; c6_c_k < 6; c6_c_k++) {
        c6_b_k = 1.0 + (real_T)c6_c_k;
        c6_U[((int32_T)c6_b_k + 6 * ((int32_T)c6_b_k - 1)) - 1] += 15120.0;
      }

      for (c6_i138 = 0; c6_i138 < 36; c6_i138++) {
        c6_y[c6_i138] = c6_U[c6_i138];
      }

      c6_b_eml_scalar_eg(chartInstance);
      c6_b_eml_scalar_eg(chartInstance);
      for (c6_i139 = 0; c6_i139 < 36; c6_i139++) {
        c6_U[c6_i139] = 0.0;
      }

      for (c6_i140 = 0; c6_i140 < 36; c6_i140++) {
        c6_e_A[c6_i140] = c6_A[c6_i140];
      }

      for (c6_i141 = 0; c6_i141 < 36; c6_i141++) {
        c6_c_y[c6_i141] = c6_y[c6_i141];
      }

      c6_b_eml_xgemm(chartInstance, c6_e_A, c6_c_y, c6_U);
      for (c6_i142 = 0; c6_i142 < 36; c6_i142++) {
        c6_A3[c6_i142] *= 30.0;
      }

      for (c6_i143 = 0; c6_i143 < 36; c6_i143++) {
        c6_A2[c6_i143] *= 3360.0;
      }

      for (c6_i144 = 0; c6_i144 < 36; c6_i144++) {
        c6_F[c6_i144] = c6_A3[c6_i144] + c6_A2[c6_i144];
      }

      c6_d = 30240.0;
    } else {
      c6_b_eml_scalar_eg(chartInstance);
      c6_b_eml_scalar_eg(chartInstance);
      for (c6_i145 = 0; c6_i145 < 36; c6_i145++) {
        c6_A4[c6_i145] = 0.0;
      }

      for (c6_i146 = 0; c6_i146 < 36; c6_i146++) {
        c6_b_A3[c6_i146] = c6_A3[c6_i146];
      }

      for (c6_i147 = 0; c6_i147 < 36; c6_i147++) {
        c6_d_A2[c6_i147] = c6_A2[c6_i147];
      }

      c6_b_eml_xgemm(chartInstance, c6_b_A3, c6_d_A2, c6_A4);
      if (c6_m == 7.0) {
        for (c6_i148 = 0; c6_i148 < 36; c6_i148++) {
          c6_U[c6_i148] = 1512.0 * c6_A3[c6_i148];
        }

        for (c6_i149 = 0; c6_i149 < 36; c6_i149++) {
          c6_d_y[c6_i149] = 277200.0 * c6_A2[c6_i149];
        }

        for (c6_i150 = 0; c6_i150 < 36; c6_i150++) {
          c6_U[c6_i150] = (c6_A4[c6_i150] + c6_U[c6_i150]) + c6_d_y[c6_i150];
        }

        for (c6_d_k = 0; c6_d_k < 6; c6_d_k++) {
          c6_b_k = 1.0 + (real_T)c6_d_k;
          c6_U[((int32_T)c6_b_k + 6 * ((int32_T)c6_b_k - 1)) - 1] += 8.64864E+6;
        }

        for (c6_i151 = 0; c6_i151 < 36; c6_i151++) {
          c6_y[c6_i151] = c6_U[c6_i151];
        }

        c6_b_eml_scalar_eg(chartInstance);
        c6_b_eml_scalar_eg(chartInstance);
        for (c6_i152 = 0; c6_i152 < 36; c6_i152++) {
          c6_U[c6_i152] = 0.0;
        }

        for (c6_i153 = 0; c6_i153 < 36; c6_i153++) {
          c6_f_A[c6_i153] = c6_A[c6_i153];
        }

        for (c6_i154 = 0; c6_i154 < 36; c6_i154++) {
          c6_e_y[c6_i154] = c6_y[c6_i154];
        }

        c6_b_eml_xgemm(chartInstance, c6_f_A, c6_e_y, c6_U);
        for (c6_i155 = 0; c6_i155 < 36; c6_i155++) {
          c6_A4[c6_i155] *= 56.0;
        }

        for (c6_i156 = 0; c6_i156 < 36; c6_i156++) {
          c6_A3[c6_i156] *= 25200.0;
        }

        for (c6_i157 = 0; c6_i157 < 36; c6_i157++) {
          c6_A2[c6_i157] *= 1.99584E+6;
        }

        for (c6_i158 = 0; c6_i158 < 36; c6_i158++) {
          c6_F[c6_i158] = (c6_A4[c6_i158] + c6_A3[c6_i158]) + c6_A2[c6_i158];
        }

        c6_d = 1.729728E+7;
      } else if (c6_m == 9.0) {
        c6_b_eml_scalar_eg(chartInstance);
        c6_b_eml_scalar_eg(chartInstance);
        for (c6_i159 = 0; c6_i159 < 36; c6_i159++) {
          c6_F[c6_i159] = 0.0;
        }

        for (c6_i160 = 0; c6_i160 < 36; c6_i160++) {
          c6_b_A4[c6_i160] = c6_A4[c6_i160];
        }

        for (c6_i161 = 0; c6_i161 < 36; c6_i161++) {
          c6_e_A2[c6_i161] = c6_A2[c6_i161];
        }

        c6_b_eml_xgemm(chartInstance, c6_b_A4, c6_e_A2, c6_F);
        for (c6_i162 = 0; c6_i162 < 36; c6_i162++) {
          c6_U[c6_i162] = 3960.0 * c6_A4[c6_i162];
        }

        for (c6_i163 = 0; c6_i163 < 36; c6_i163++) {
          c6_d_y[c6_i163] = 2.16216E+6 * c6_A3[c6_i163];
        }

        for (c6_i164 = 0; c6_i164 < 36; c6_i164++) {
          c6_f_y[c6_i164] = 3.027024E+8 * c6_A2[c6_i164];
        }

        for (c6_i165 = 0; c6_i165 < 36; c6_i165++) {
          c6_U[c6_i165] = ((c6_F[c6_i165] + c6_U[c6_i165]) + c6_d_y[c6_i165]) +
            c6_f_y[c6_i165];
        }

        for (c6_e_k = 0; c6_e_k < 6; c6_e_k++) {
          c6_b_k = 1.0 + (real_T)c6_e_k;
          c6_U[((int32_T)c6_b_k + 6 * ((int32_T)c6_b_k - 1)) - 1] +=
            8.8216128E+9;
        }

        for (c6_i166 = 0; c6_i166 < 36; c6_i166++) {
          c6_y[c6_i166] = c6_U[c6_i166];
        }

        c6_b_eml_scalar_eg(chartInstance);
        c6_b_eml_scalar_eg(chartInstance);
        for (c6_i167 = 0; c6_i167 < 36; c6_i167++) {
          c6_U[c6_i167] = 0.0;
        }

        for (c6_i168 = 0; c6_i168 < 36; c6_i168++) {
          c6_g_A[c6_i168] = c6_A[c6_i168];
        }

        for (c6_i169 = 0; c6_i169 < 36; c6_i169++) {
          c6_g_y[c6_i169] = c6_y[c6_i169];
        }

        c6_b_eml_xgemm(chartInstance, c6_g_A, c6_g_y, c6_U);
        for (c6_i170 = 0; c6_i170 < 36; c6_i170++) {
          c6_F[c6_i170] *= 90.0;
        }

        for (c6_i171 = 0; c6_i171 < 36; c6_i171++) {
          c6_A4[c6_i171] *= 110880.0;
        }

        for (c6_i172 = 0; c6_i172 < 36; c6_i172++) {
          c6_A3[c6_i172] *= 3.027024E+7;
        }

        for (c6_i173 = 0; c6_i173 < 36; c6_i173++) {
          c6_A2[c6_i173] *= 2.0756736E+9;
        }

        for (c6_i174 = 0; c6_i174 < 36; c6_i174++) {
          c6_F[c6_i174] = ((c6_F[c6_i174] + c6_A4[c6_i174]) + c6_A3[c6_i174]) +
            c6_A2[c6_i174];
        }

        c6_d = 1.76432256E+10;
      } else {
        for (c6_i175 = 0; c6_i175 < 36; c6_i175++) {
          c6_U[c6_i175] = 3.352212864E+10 * c6_A4[c6_i175];
        }

        for (c6_i176 = 0; c6_i176 < 36; c6_i176++) {
          c6_d_y[c6_i176] = 1.05594705216E+13 * c6_A3[c6_i176];
        }

        for (c6_i177 = 0; c6_i177 < 36; c6_i177++) {
          c6_f_y[c6_i177] = 1.1873537964288E+15 * c6_A2[c6_i177];
        }

        for (c6_i178 = 0; c6_i178 < 36; c6_i178++) {
          c6_U[c6_i178] = (c6_U[c6_i178] + c6_d_y[c6_i178]) + c6_f_y[c6_i178];
        }

        for (c6_f_k = 0; c6_f_k < 6; c6_f_k++) {
          c6_b_k = 1.0 + (real_T)c6_f_k;
          c6_U[((int32_T)c6_b_k + 6 * ((int32_T)c6_b_k - 1)) - 1] +=
            3.238237626624E+16;
        }

        for (c6_i179 = 0; c6_i179 < 36; c6_i179++) {
          c6_d_y[c6_i179] = 16380.0 * c6_A3[c6_i179];
        }

        for (c6_i180 = 0; c6_i180 < 36; c6_i180++) {
          c6_f_y[c6_i180] = 4.08408E+7 * c6_A2[c6_i180];
        }

        for (c6_i181 = 0; c6_i181 < 36; c6_i181++) {
          c6_d_y[c6_i181] = (c6_A4[c6_i181] + c6_d_y[c6_i181]) + c6_f_y[c6_i181];
        }

        c6_b_eml_scalar_eg(chartInstance);
        c6_b_eml_scalar_eg(chartInstance);
        for (c6_i182 = 0; c6_i182 < 36; c6_i182++) {
          c6_f_y[c6_i182] = 0.0;
        }

        for (c6_i183 = 0; c6_i183 < 36; c6_i183++) {
          c6_c_A4[c6_i183] = c6_A4[c6_i183];
        }

        for (c6_i184 = 0; c6_i184 < 36; c6_i184++) {
          c6_h_y[c6_i184] = c6_d_y[c6_i184];
        }

        c6_b_eml_xgemm(chartInstance, c6_c_A4, c6_h_y, c6_f_y);
        for (c6_i185 = 0; c6_i185 < 36; c6_i185++) {
          c6_f_y[c6_i185] += c6_U[c6_i185];
        }

        c6_b_eml_scalar_eg(chartInstance);
        c6_b_eml_scalar_eg(chartInstance);
        for (c6_i186 = 0; c6_i186 < 36; c6_i186++) {
          c6_U[c6_i186] = 0.0;
        }

        for (c6_i187 = 0; c6_i187 < 36; c6_i187++) {
          c6_h_A[c6_i187] = c6_A[c6_i187];
        }

        for (c6_i188 = 0; c6_i188 < 36; c6_i188++) {
          c6_i_y[c6_i188] = c6_f_y[c6_i188];
        }

        c6_b_eml_xgemm(chartInstance, c6_h_A, c6_i_y, c6_U);
        for (c6_i189 = 0; c6_i189 < 36; c6_i189++) {
          c6_d_y[c6_i189] = 182.0 * c6_A4[c6_i189];
        }

        for (c6_i190 = 0; c6_i190 < 36; c6_i190++) {
          c6_f_y[c6_i190] = 960960.0 * c6_A3[c6_i190];
        }

        for (c6_i191 = 0; c6_i191 < 36; c6_i191++) {
          c6_y[c6_i191] = 1.32324192E+9 * c6_A2[c6_i191];
        }

        for (c6_i192 = 0; c6_i192 < 36; c6_i192++) {
          c6_d_y[c6_i192] = (c6_d_y[c6_i192] + c6_f_y[c6_i192]) + c6_y[c6_i192];
        }

        c6_b_eml_scalar_eg(chartInstance);
        c6_b_eml_scalar_eg(chartInstance);
        for (c6_i193 = 0; c6_i193 < 36; c6_i193++) {
          c6_F[c6_i193] = 0.0;
        }

        for (c6_i194 = 0; c6_i194 < 36; c6_i194++) {
          c6_d_A4[c6_i194] = c6_A4[c6_i194];
        }

        for (c6_i195 = 0; c6_i195 < 36; c6_i195++) {
          c6_j_y[c6_i195] = c6_d_y[c6_i195];
        }

        c6_b_eml_xgemm(chartInstance, c6_d_A4, c6_j_y, c6_F);
        for (c6_i196 = 0; c6_i196 < 36; c6_i196++) {
          c6_A4[c6_i196] *= 6.704425728E+11;
        }

        for (c6_i197 = 0; c6_i197 < 36; c6_i197++) {
          c6_A3[c6_i197] *= 1.29060195264E+14;
        }

        for (c6_i198 = 0; c6_i198 < 36; c6_i198++) {
          c6_A2[c6_i198] *= 7.7717703038976E+15;
        }

        for (c6_i199 = 0; c6_i199 < 36; c6_i199++) {
          c6_F[c6_i199] = ((c6_F[c6_i199] + c6_A4[c6_i199]) + c6_A3[c6_i199]) +
            c6_A2[c6_i199];
        }

        c6_d = 6.476475253248E+16;
      }
    }
  }

  for (c6_g_k = 0; c6_g_k < 6; c6_g_k++) {
    c6_b_k = 1.0 + (real_T)c6_g_k;
    c6_F[((int32_T)c6_b_k + 6 * ((int32_T)c6_b_k - 1)) - 1] += c6_d;
  }

  for (c6_h_k = 0; c6_h_k < 36; c6_h_k++) {
    c6_b_k = 1.0 + (real_T)c6_h_k;
    c6_uk = c6_U[(int32_T)c6_b_k - 1];
    c6_U[(int32_T)c6_b_k - 1] = c6_F[(int32_T)c6_b_k - 1] - c6_uk;
    c6_F[(int32_T)c6_b_k - 1] += c6_uk;
  }

  c6_b_eml_matlab_zgetrf(chartInstance, c6_U, c6_ipiv, &c6_info);
  c6_b_info = c6_info;
  c6_c_info = c6_b_info;
  c6_d_info = c6_c_info;
  if (c6_d_info > 0) {
    c6_eml_warning(chartInstance);
  }

  c6_b_eml_scalar_eg(chartInstance);
  for (c6_xi = 1; c6_xi < 6; c6_xi++) {
    c6_b_xi = c6_xi - 1;
    if (c6_ipiv[c6_b_xi] != c6_b_xi + 1) {
      c6_ip = c6_ipiv[c6_b_xi];
      for (c6_xj = 1; c6_xj < 7; c6_xj++) {
        c6_b_xj = c6_xj - 1;
        c6_temp = c6_F[c6_b_xi + 6 * c6_b_xj];
        c6_F[c6_b_xi + 6 * c6_b_xj] = c6_F[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
          c6_ip, 1, 6, 1, 0) + 6 * c6_b_xj) - 1];
        c6_F[(c6_ip + 6 * c6_b_xj) - 1] = c6_temp;
      }
    }
  }

  for (c6_i200 = 0; c6_i200 < 36; c6_i200++) {
    c6_b_U[c6_i200] = c6_U[c6_i200];
  }

  c6_c_eml_xtrsm(chartInstance, c6_b_U, c6_F);
  for (c6_i201 = 0; c6_i201 < 36; c6_i201++) {
    c6_c_U[c6_i201] = c6_U[c6_i201];
  }

  c6_d_eml_xtrsm(chartInstance, c6_c_U, c6_F);
}

static void c6_b_eml_scalar_eg(SFc6_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c6_eml_xgemm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36], real_T c6_C[36], real_T c6_b_C[36])
{
  int32_T c6_i202;
  int32_T c6_i203;
  real_T c6_b_A[36];
  int32_T c6_i204;
  real_T c6_b_B[36];
  for (c6_i202 = 0; c6_i202 < 36; c6_i202++) {
    c6_b_C[c6_i202] = c6_C[c6_i202];
  }

  for (c6_i203 = 0; c6_i203 < 36; c6_i203++) {
    c6_b_A[c6_i203] = c6_A[c6_i203];
  }

  for (c6_i204 = 0; c6_i204 < 36; c6_i204++) {
    c6_b_B[c6_i204] = c6_B[c6_i204];
  }

  c6_b_eml_xgemm(chartInstance, c6_b_A, c6_b_B, c6_b_C);
}

static void c6_eml_matlab_zgetrf(SFc6_Model_02InstanceStruct *chartInstance,
  real_T c6_A[36], real_T c6_b_A[36], int32_T c6_ipiv[6], int32_T *c6_info)
{
  int32_T c6_i205;
  for (c6_i205 = 0; c6_i205 < 36; c6_i205++) {
    c6_b_A[c6_i205] = c6_A[c6_i205];
  }

  c6_b_eml_matlab_zgetrf(chartInstance, c6_b_A, c6_ipiv, c6_info);
}

static int32_T c6_eml_ixamax(SFc6_Model_02InstanceStruct *chartInstance, int32_T
  c6_n, real_T c6_x[36], int32_T c6_ix0)
{
  int32_T c6_idxmax;
  int32_T c6_b_n;
  int32_T c6_b_ix0;
  int32_T c6_c_n;
  int32_T c6_c_ix0;
  int32_T c6_ix;
  real_T c6_b_x;
  real_T c6_c_x;
  real_T c6_d_x;
  real_T c6_y;
  real_T c6_e_x;
  real_T c6_f_x;
  real_T c6_b_y;
  real_T c6_smax;
  int32_T c6_d_n;
  int32_T c6_b;
  int32_T c6_b_b;
  boolean_T c6_overflow;
  int32_T c6_k;
  int32_T c6_b_k;
  int32_T c6_a;
  real_T c6_g_x;
  real_T c6_h_x;
  real_T c6_i_x;
  real_T c6_c_y;
  real_T c6_j_x;
  real_T c6_k_x;
  real_T c6_d_y;
  real_T c6_s;
  c6_b_n = c6_n;
  c6_b_ix0 = c6_ix0;
  c6_c_n = c6_b_n;
  c6_c_ix0 = c6_b_ix0;
  if (c6_c_n < 1) {
    c6_idxmax = 0;
  } else {
    c6_idxmax = 1;
    if (c6_c_n > 1) {
      c6_ix = c6_c_ix0;
      c6_b_x = c6_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_ix, 1, 36, 1, 0) - 1];
      c6_c_x = c6_b_x;
      c6_d_x = c6_c_x;
      c6_y = muDoubleScalarAbs(c6_d_x);
      c6_e_x = 0.0;
      c6_f_x = c6_e_x;
      c6_b_y = muDoubleScalarAbs(c6_f_x);
      c6_smax = c6_y + c6_b_y;
      c6_d_n = c6_c_n;
      c6_b = c6_d_n;
      c6_b_b = c6_b;
      c6_eml_switch_helper(chartInstance);
      c6_overflow = (c6_b_b > 2147483646);
      if (c6_overflow) {
        c6_check_forloop_overflow_error(chartInstance, true);
      }

      for (c6_k = 2; c6_k <= c6_d_n; c6_k++) {
        c6_b_k = c6_k;
        c6_a = c6_ix + 1;
        c6_ix = c6_a;
        c6_g_x = c6_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_ix, 1, 36, 1, 0) - 1];
        c6_h_x = c6_g_x;
        c6_i_x = c6_h_x;
        c6_c_y = muDoubleScalarAbs(c6_i_x);
        c6_j_x = 0.0;
        c6_k_x = c6_j_x;
        c6_d_y = muDoubleScalarAbs(c6_k_x);
        c6_s = c6_c_y + c6_d_y;
        if (c6_s > c6_smax) {
          c6_idxmax = c6_b_k;
          c6_smax = c6_s;
        }
      }
    }
  }

  return c6_idxmax;
}

static void c6_check_forloop_overflow_error(SFc6_Model_02InstanceStruct
  *chartInstance, boolean_T c6_overflow)
{
  int32_T c6_i206;
  static char_T c6_cv0[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c6_u[34];
  const mxArray *c6_y = NULL;
  int32_T c6_i207;
  static char_T c6_cv1[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c6_b_u[23];
  const mxArray *c6_b_y = NULL;
  (void)chartInstance;
  (void)c6_overflow;
  for (c6_i206 = 0; c6_i206 < 34; c6_i206++) {
    c6_u[c6_i206] = c6_cv0[c6_i206];
  }

  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u, 10, 0U, 1U, 0U, 2, 1, 34), false);
  for (c6_i207 = 0; c6_i207 < 23; c6_i207++) {
    c6_b_u[c6_i207] = c6_cv1[c6_i207];
  }

  c6_b_y = NULL;
  sf_mex_assign(&c6_b_y, sf_mex_create("y", c6_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c6_y, 14, c6_b_y));
}

static void c6_threshold(SFc6_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c6_eml_xgeru(SFc6_Model_02InstanceStruct *chartInstance, int32_T
  c6_m, int32_T c6_n, real_T c6_alpha1, int32_T c6_ix0, int32_T c6_iy0, real_T
  c6_A[36], int32_T c6_ia0, real_T c6_b_A[36])
{
  int32_T c6_i208;
  for (c6_i208 = 0; c6_i208 < 36; c6_i208++) {
    c6_b_A[c6_i208] = c6_A[c6_i208];
  }

  c6_b_eml_xgeru(chartInstance, c6_m, c6_n, c6_alpha1, c6_ix0, c6_iy0, c6_b_A,
                 c6_ia0);
}

static void c6_eml_warning(SFc6_Model_02InstanceStruct *chartInstance)
{
  int32_T c6_i209;
  static char_T c6_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c6_u[27];
  const mxArray *c6_y = NULL;
  (void)chartInstance;
  for (c6_i209 = 0; c6_i209 < 27; c6_i209++) {
    c6_u[c6_i209] = c6_varargin_1[c6_i209];
  }

  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u, 10, 0U, 1U, 0U, 2, 1, 27), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c6_y));
}

static void c6_eml_xtrsm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36], real_T c6_b_B[36])
{
  int32_T c6_i210;
  int32_T c6_i211;
  real_T c6_b_A[36];
  for (c6_i210 = 0; c6_i210 < 36; c6_i210++) {
    c6_b_B[c6_i210] = c6_B[c6_i210];
  }

  for (c6_i211 = 0; c6_i211 < 36; c6_i211++) {
    c6_b_A[c6_i211] = c6_A[c6_i211];
  }

  c6_c_eml_xtrsm(chartInstance, c6_b_A, c6_b_B);
}

static void c6_b_threshold(SFc6_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c6_scalarEg(SFc6_Model_02InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c6_b_eml_xtrsm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36], real_T c6_b_B[36])
{
  int32_T c6_i212;
  int32_T c6_i213;
  real_T c6_b_A[36];
  for (c6_i212 = 0; c6_i212 < 36; c6_i212++) {
    c6_b_B[c6_i212] = c6_B[c6_i212];
  }

  for (c6_i213 = 0; c6_i213 < 36; c6_i213++) {
    c6_b_A[c6_i213] = c6_A[c6_i213];
  }

  c6_d_eml_xtrsm(chartInstance, c6_b_A, c6_b_B);
}

static const mxArray *c6_h_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  int32_T c6_u;
  const mxArray *c6_y = NULL;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_u = *(int32_T *)c6_inData;
  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", &c6_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, false);
  return c6_mxArrayOutData;
}

static int32_T c6_i_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId)
{
  int32_T c6_y;
  int32_T c6_i214;
  (void)chartInstance;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), &c6_i214, 1, 6, 0U, 0, 0U, 0);
  c6_y = c6_i214;
  sf_mex_destroy(&c6_u);
  return c6_y;
}

static void c6_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_b_sfEvent;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  int32_T c6_y;
  SFc6_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc6_Model_02InstanceStruct *)chartInstanceVoid;
  c6_b_sfEvent = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_y = c6_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_b_sfEvent),
    &c6_thisId);
  sf_mex_destroy(&c6_b_sfEvent);
  *(int32_T *)c6_outData = c6_y;
  sf_mex_destroy(&c6_mxArrayInData);
}

static uint8_T c6_j_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_b_is_active_c6_Model_02, const char_T *c6_identifier)
{
  uint8_T c6_y;
  emlrtMsgIdentifier c6_thisId;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_y = c6_k_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c6_b_is_active_c6_Model_02), &c6_thisId);
  sf_mex_destroy(&c6_b_is_active_c6_Model_02);
  return c6_y;
}

static uint8_T c6_k_emlrt_marshallIn(SFc6_Model_02InstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId)
{
  uint8_T c6_y;
  uint8_T c6_u0;
  (void)chartInstance;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), &c6_u0, 1, 3, 0U, 0, 0U, 0);
  c6_y = c6_u0;
  sf_mex_destroy(&c6_u);
  return c6_y;
}

static void c6_b_eml_xgemm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36], real_T c6_C[36])
{
  int32_T c6_i215;
  int32_T c6_i216;
  int32_T c6_i217;
  int32_T c6_i218;
  int32_T c6_i219;
  (void)chartInstance;
  for (c6_i215 = 0; c6_i215 < 6; c6_i215++) {
    c6_i216 = 0;
    for (c6_i217 = 0; c6_i217 < 6; c6_i217++) {
      c6_C[c6_i216 + c6_i215] = 0.0;
      c6_i218 = 0;
      for (c6_i219 = 0; c6_i219 < 6; c6_i219++) {
        c6_C[c6_i216 + c6_i215] += c6_A[c6_i218 + c6_i215] * c6_B[c6_i219 +
          c6_i216];
        c6_i218 += 6;
      }

      c6_i216 += 6;
    }
  }
}

static void c6_b_eml_matlab_zgetrf(SFc6_Model_02InstanceStruct *chartInstance,
  real_T c6_A[36], int32_T c6_ipiv[6], int32_T *c6_info)
{
  int32_T c6_i220;
  int32_T c6_j;
  int32_T c6_b_j;
  int32_T c6_a;
  int32_T c6_b_a;
  int32_T c6_jm1;
  int32_T c6_b;
  int32_T c6_b_b;
  int32_T c6_mmj;
  int32_T c6_c_a;
  int32_T c6_d_a;
  int32_T c6_c;
  int32_T c6_c_b;
  int32_T c6_d_b;
  int32_T c6_jj;
  int32_T c6_e_a;
  int32_T c6_f_a;
  int32_T c6_jp1j;
  int32_T c6_g_a;
  int32_T c6_h_a;
  int32_T c6_b_c;
  int32_T c6_i221;
  int32_T c6_i222;
  int32_T c6_i223;
  real_T c6_b_A[36];
  int32_T c6_i_a;
  int32_T c6_j_a;
  int32_T c6_jpiv_offset;
  int32_T c6_k_a;
  int32_T c6_e_b;
  int32_T c6_l_a;
  int32_T c6_f_b;
  int32_T c6_jpiv;
  int32_T c6_m_a;
  int32_T c6_g_b;
  int32_T c6_n_a;
  int32_T c6_h_b;
  int32_T c6_c_c;
  int32_T c6_i_b;
  int32_T c6_j_b;
  int32_T c6_jrow;
  int32_T c6_o_a;
  int32_T c6_k_b;
  int32_T c6_p_a;
  int32_T c6_l_b;
  int32_T c6_jprow;
  int32_T c6_ix0;
  int32_T c6_iy0;
  int32_T c6_b_ix0;
  int32_T c6_b_iy0;
  int32_T c6_c_ix0;
  int32_T c6_c_iy0;
  int32_T c6_ix;
  int32_T c6_iy;
  int32_T c6_k;
  real_T c6_temp;
  int32_T c6_q_a;
  int32_T c6_r_a;
  int32_T c6_b_jp1j;
  int32_T c6_s_a;
  int32_T c6_t_a;
  int32_T c6_d_c;
  int32_T c6_u_a;
  int32_T c6_m_b;
  int32_T c6_v_a;
  int32_T c6_n_b;
  int32_T c6_i224;
  int32_T c6_w_a;
  int32_T c6_o_b;
  int32_T c6_x_a;
  int32_T c6_p_b;
  boolean_T c6_overflow;
  int32_T c6_i;
  int32_T c6_b_i;
  real_T c6_x;
  real_T c6_y;
  real_T c6_b_x;
  real_T c6_b_y;
  real_T c6_z;
  int32_T c6_q_b;
  int32_T c6_r_b;
  int32_T c6_e_c;
  int32_T c6_y_a;
  int32_T c6_ab_a;
  int32_T c6_f_c;
  int32_T c6_bb_a;
  int32_T c6_cb_a;
  int32_T c6_g_c;
  real_T c6_d3;
  for (c6_i220 = 0; c6_i220 < 6; c6_i220++) {
    c6_ipiv[c6_i220] = 1 + c6_i220;
  }

  *c6_info = 0;
  for (c6_j = 1; c6_j < 6; c6_j++) {
    c6_b_j = c6_j;
    c6_a = c6_b_j;
    c6_b_a = c6_a - 1;
    c6_jm1 = c6_b_a;
    c6_b = c6_b_j;
    c6_b_b = c6_b;
    c6_mmj = 6 - c6_b_b;
    c6_c_a = c6_jm1;
    c6_d_a = c6_c_a;
    c6_c = c6_d_a * 7;
    c6_c_b = c6_c;
    c6_d_b = c6_c_b + 1;
    c6_jj = c6_d_b;
    c6_e_a = c6_jj;
    c6_f_a = c6_e_a + 1;
    c6_jp1j = c6_f_a;
    c6_g_a = c6_mmj;
    c6_h_a = c6_g_a;
    c6_b_c = c6_h_a;
    c6_i221 = 0;
    for (c6_i222 = 0; c6_i222 < 6; c6_i222++) {
      for (c6_i223 = 0; c6_i223 < 6; c6_i223++) {
        c6_b_A[c6_i223 + c6_i221] = c6_A[c6_i223 + c6_i221];
      }

      c6_i221 += 6;
    }

    c6_i_a = c6_eml_ixamax(chartInstance, c6_b_c + 1, c6_b_A, c6_jj);
    c6_j_a = c6_i_a - 1;
    c6_jpiv_offset = c6_j_a;
    c6_k_a = c6_jj;
    c6_e_b = c6_jpiv_offset;
    c6_l_a = c6_k_a;
    c6_f_b = c6_e_b;
    c6_jpiv = c6_l_a + c6_f_b;
    if (c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_jpiv, 1, 36, 1, 0) - 1] != 0.0)
    {
      if (c6_jpiv_offset != 0) {
        c6_m_a = c6_b_j;
        c6_g_b = c6_jpiv_offset;
        c6_n_a = c6_m_a;
        c6_h_b = c6_g_b;
        c6_c_c = c6_n_a + c6_h_b;
        c6_ipiv[c6_b_j - 1] = c6_c_c;
        c6_i_b = c6_jm1;
        c6_j_b = c6_i_b + 1;
        c6_jrow = c6_j_b;
        c6_o_a = c6_jrow;
        c6_k_b = c6_jpiv_offset;
        c6_p_a = c6_o_a;
        c6_l_b = c6_k_b;
        c6_jprow = c6_p_a + c6_l_b;
        c6_ix0 = c6_jrow;
        c6_iy0 = c6_jprow;
        c6_b_ix0 = c6_ix0;
        c6_b_iy0 = c6_iy0;
        c6_threshold(chartInstance);
        c6_c_ix0 = c6_b_ix0;
        c6_c_iy0 = c6_b_iy0;
        c6_ix = c6_c_ix0;
        c6_iy = c6_c_iy0;
        for (c6_k = 1; c6_k < 7; c6_k++) {
          c6_temp = c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_ix, 1, 36, 1, 0) - 1];
          c6_A[c6_ix - 1] = c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_iy, 1, 36, 1,
            0) - 1];
          c6_A[c6_iy - 1] = c6_temp;
          c6_q_a = c6_ix + 6;
          c6_ix = c6_q_a;
          c6_r_a = c6_iy + 6;
          c6_iy = c6_r_a;
        }
      }

      c6_b_jp1j = c6_jp1j;
      c6_s_a = c6_mmj;
      c6_t_a = c6_s_a;
      c6_d_c = c6_t_a;
      c6_u_a = c6_jp1j;
      c6_m_b = c6_d_c - 1;
      c6_v_a = c6_u_a;
      c6_n_b = c6_m_b;
      c6_i224 = c6_v_a + c6_n_b;
      c6_w_a = c6_b_jp1j;
      c6_o_b = c6_i224;
      c6_x_a = c6_w_a;
      c6_p_b = c6_o_b;
      if (c6_x_a > c6_p_b) {
        c6_overflow = false;
      } else {
        c6_eml_switch_helper(chartInstance);
        c6_overflow = (c6_p_b > 2147483646);
      }

      if (c6_overflow) {
        c6_check_forloop_overflow_error(chartInstance, true);
      }

      for (c6_i = c6_b_jp1j; c6_i <= c6_i224; c6_i++) {
        c6_b_i = c6_i;
        c6_x = c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_b_i, 1, 36, 1, 0) - 1];
        c6_y = c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_jj, 1, 36, 1, 0) - 1];
        c6_b_x = c6_x;
        c6_b_y = c6_y;
        c6_z = c6_b_x / c6_b_y;
        c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_b_i, 1, 36, 1, 0) - 1] = c6_z;
      }
    } else {
      *c6_info = c6_b_j;
    }

    c6_q_b = c6_b_j;
    c6_r_b = c6_q_b;
    c6_e_c = 6 - c6_r_b;
    c6_y_a = c6_jj;
    c6_ab_a = c6_y_a;
    c6_f_c = c6_ab_a;
    c6_bb_a = c6_jj;
    c6_cb_a = c6_bb_a;
    c6_g_c = c6_cb_a;
    c6_d3 = -1.0;
    c6_b_eml_xgeru(chartInstance, c6_mmj, c6_e_c, c6_d3, c6_jp1j, c6_f_c + 6,
                   c6_A, c6_g_c + 7);
  }

  if (*c6_info == 0) {
    if (!(c6_A[35] != 0.0)) {
      *c6_info = 6;
    }
  }
}

static void c6_b_eml_xgeru(SFc6_Model_02InstanceStruct *chartInstance, int32_T
  c6_m, int32_T c6_n, real_T c6_alpha1, int32_T c6_ix0, int32_T c6_iy0, real_T
  c6_A[36], int32_T c6_ia0)
{
  int32_T c6_b_m;
  int32_T c6_b_n;
  int32_T c6_b_ix0;
  int32_T c6_b_iy0;
  int32_T c6_b_ia0;
  int32_T c6_c_m;
  int32_T c6_c_n;
  int32_T c6_c_ix0;
  int32_T c6_c_iy0;
  int32_T c6_c_ia0;
  int32_T c6_d_m;
  int32_T c6_d_n;
  int32_T c6_d_ix0;
  int32_T c6_d_iy0;
  int32_T c6_d_ia0;
  int32_T c6_e_m;
  int32_T c6_e_n;
  int32_T c6_e_ix0;
  int32_T c6_e_iy0;
  int32_T c6_e_ia0;
  int32_T c6_ixstart;
  int32_T c6_a;
  int32_T c6_jA;
  int32_T c6_jy;
  int32_T c6_f_n;
  int32_T c6_b;
  int32_T c6_b_b;
  int32_T c6_j;
  real_T c6_yjy;
  real_T c6_temp;
  int32_T c6_ix;
  int32_T c6_c_b;
  int32_T c6_i225;
  int32_T c6_b_a;
  int32_T c6_d_b;
  int32_T c6_i226;
  int32_T c6_c_a;
  int32_T c6_e_b;
  int32_T c6_d_a;
  int32_T c6_f_b;
  boolean_T c6_overflow;
  int32_T c6_ijA;
  int32_T c6_b_ijA;
  int32_T c6_e_a;
  int32_T c6_f_a;
  int32_T c6_g_a;
  (void)c6_alpha1;
  c6_b_m = c6_m;
  c6_b_n = c6_n;
  c6_b_ix0 = c6_ix0;
  c6_b_iy0 = c6_iy0;
  c6_b_ia0 = c6_ia0;
  c6_c_m = c6_b_m;
  c6_c_n = c6_b_n;
  c6_c_ix0 = c6_b_ix0;
  c6_c_iy0 = c6_b_iy0;
  c6_c_ia0 = c6_b_ia0;
  c6_d_m = c6_c_m;
  c6_d_n = c6_c_n;
  c6_d_ix0 = c6_c_ix0;
  c6_d_iy0 = c6_c_iy0;
  c6_d_ia0 = c6_c_ia0;
  c6_e_m = c6_d_m;
  c6_e_n = c6_d_n;
  c6_e_ix0 = c6_d_ix0;
  c6_e_iy0 = c6_d_iy0;
  c6_e_ia0 = c6_d_ia0;
  c6_ixstart = c6_e_ix0;
  c6_a = c6_e_ia0 - 1;
  c6_jA = c6_a;
  c6_jy = c6_e_iy0;
  c6_f_n = c6_e_n;
  c6_b = c6_f_n;
  c6_b_b = c6_b;
  if (1 > c6_b_b) {
  } else {
    c6_eml_switch_helper(chartInstance);
  }

  for (c6_j = 1; c6_j <= c6_f_n; c6_j++) {
    c6_yjy = c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_jy, 1, 36, 1, 0) - 1];
    if (c6_yjy != 0.0) {
      c6_temp = -c6_yjy;
      c6_ix = c6_ixstart;
      c6_c_b = c6_jA + 1;
      c6_i225 = c6_c_b;
      c6_b_a = c6_e_m;
      c6_d_b = c6_jA;
      c6_i226 = c6_b_a + c6_d_b;
      c6_c_a = c6_i225;
      c6_e_b = c6_i226;
      c6_d_a = c6_c_a;
      c6_f_b = c6_e_b;
      if (c6_d_a > c6_f_b) {
        c6_overflow = false;
      } else {
        c6_eml_switch_helper(chartInstance);
        c6_overflow = (c6_f_b > 2147483646);
      }

      if (c6_overflow) {
        c6_check_forloop_overflow_error(chartInstance, true);
      }

      for (c6_ijA = c6_i225; c6_ijA <= c6_i226; c6_ijA++) {
        c6_b_ijA = c6_ijA;
        c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_b_ijA, 1, 36, 1, 0) - 1] =
          c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_b_ijA, 1, 36, 1, 0) - 1] +
          c6_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c6_ix, 1, 36, 1, 0) - 1] *
          c6_temp;
        c6_e_a = c6_ix + 1;
        c6_ix = c6_e_a;
      }
    }

    c6_f_a = c6_jy + 6;
    c6_jy = c6_f_a;
    c6_g_a = c6_jA + 6;
    c6_jA = c6_g_a;
  }
}

static void c6_c_eml_xtrsm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36])
{
  int32_T c6_j;
  int32_T c6_b_j;
  int32_T c6_jBcol;
  int32_T c6_k;
  int32_T c6_b_k;
  int32_T c6_kAcol;
  int32_T c6_i227;
  int32_T c6_a;
  int32_T c6_b_a;
  int32_T c6_i;
  int32_T c6_b_i;
  c6_b_threshold(chartInstance);
  c6_scalarEg(chartInstance);
  for (c6_j = 1; c6_j < 7; c6_j++) {
    c6_b_j = c6_j - 1;
    c6_jBcol = 6 * c6_b_j - 1;
    for (c6_k = 1; c6_k < 7; c6_k++) {
      c6_b_k = c6_k;
      c6_kAcol = 6 * (c6_b_k - 1) - 1;
      if (c6_B[c6_b_k + c6_jBcol] != 0.0) {
        c6_i227 = c6_b_k + 1;
        c6_a = c6_i227;
        c6_b_a = c6_a;
        if (c6_b_a > 6) {
        } else {
          c6_eml_switch_helper(chartInstance);
        }

        for (c6_i = c6_i227; c6_i < 7; c6_i++) {
          c6_b_i = c6_i;
          c6_B[c6_b_i + c6_jBcol] -= c6_B[c6_b_k + c6_jBcol] * c6_A[c6_b_i +
            c6_kAcol];
        }
      }
    }
  }
}

static void c6_d_eml_xtrsm(SFc6_Model_02InstanceStruct *chartInstance, real_T
  c6_A[36], real_T c6_B[36])
{
  int32_T c6_j;
  int32_T c6_b_j;
  int32_T c6_jBcol;
  int32_T c6_k;
  int32_T c6_b_k;
  int32_T c6_kAcol;
  real_T c6_x;
  real_T c6_y;
  real_T c6_b_x;
  real_T c6_b_y;
  real_T c6_c_x;
  real_T c6_c_y;
  real_T c6_z;
  int32_T c6_i228;
  int32_T c6_b;
  int32_T c6_b_b;
  int32_T c6_i;
  int32_T c6_b_i;
  c6_b_threshold(chartInstance);
  c6_scalarEg(chartInstance);
  for (c6_j = 1; c6_j < 7; c6_j++) {
    c6_b_j = c6_j - 1;
    c6_jBcol = 6 * c6_b_j - 1;
    for (c6_k = 6; c6_k > 0; c6_k--) {
      c6_b_k = c6_k;
      c6_kAcol = 6 * (c6_b_k - 1) - 1;
      if (c6_B[c6_b_k + c6_jBcol] != 0.0) {
        c6_x = c6_B[c6_b_k + c6_jBcol];
        c6_y = c6_A[c6_b_k + c6_kAcol];
        c6_b_x = c6_x;
        c6_b_y = c6_y;
        c6_c_x = c6_b_x;
        c6_c_y = c6_b_y;
        c6_z = c6_c_x / c6_c_y;
        c6_B[c6_b_k + c6_jBcol] = c6_z;
        c6_i228 = c6_b_k - 1;
        c6_b = c6_i228;
        c6_b_b = c6_b;
        if (1 > c6_b_b) {
        } else {
          c6_eml_switch_helper(chartInstance);
        }

        for (c6_i = 1; c6_i <= c6_i228; c6_i++) {
          c6_b_i = c6_i;
          c6_B[c6_b_i + c6_jBcol] -= c6_B[c6_b_k + c6_jBcol] * c6_A[c6_b_i +
            c6_kAcol];
        }
      }
    }
  }
}

static void init_dsm_address_info(SFc6_Model_02InstanceStruct *chartInstance)
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

void sf_c6_Model_02_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3173891208U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2316637846U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1250221759U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2242581754U);
}

mxArray *sf_c6_Model_02_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("OkUz36wxTwoos5y8btPRvG");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(13);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(3);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c6_Model_02_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c6_Model_02_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c6_Model_02(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"Phi\",},{M[8],M[0],T\"is_active_c6_Model_02\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c6_Model_02_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc6_Model_02InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc6_Model_02InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_02MachineNumber_,
           6,
           1,
           1,
           0,
           6,
           0,
           0,
           0,
           0,
           1,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"X_a");
          _SFD_SET_DATA_PROPS(1,2,0,1,"Phi");
          _SFD_SET_DATA_PROPS(2,1,1,0,"n");
          _SFD_SET_DATA_PROPS(3,1,1,0,"t_Kalman");
          _SFD_SET_DATA_PROPS(4,1,1,0,"p");
          _SFD_SET_DATA_PROPS(5,1,1,0,"M");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",260,-1,449);
        _SFD_CV_INIT_EML_FCN(0,1,"fn_Create_Phi_r",736,-1,1354);
        _SFD_CV_INIT_EML_FCN(0,2,"fn_Create_M",1612,-1,1775);
        _SFD_CV_INIT_EML_FCN(0,3,"fn_Create_N",2003,-1,2125);
        _SFD_CV_INIT_EML_FCN(0,4,"fn_Create_Phi_t",2377,-1,2629);
        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"fn_VectorToSkewSymmetricTensor",0,-1,433);

        {
          unsigned int dimVector[1];
          dimVector[0]= 13;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c6_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 12;
          dimVector[1]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c6_sf_marshallOut,(MexInFcnForType)
            c6_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c6_d_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c6_d_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c6_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c6_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c6_n;
          real_T *c6_t_Kalman;
          real_T (*c6_X_a)[13];
          real_T (*c6_Phi)[144];
          real_T (*c6_p)[3];
          real_T (*c6_M)[9];
          c6_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 4);
          c6_p = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 3);
          c6_t_Kalman = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c6_n = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c6_Phi = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 1);
          c6_X_a = (real_T (*)[13])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c6_X_a);
          _SFD_SET_DATA_VALUE_PTR(1U, *c6_Phi);
          _SFD_SET_DATA_VALUE_PTR(2U, c6_n);
          _SFD_SET_DATA_VALUE_PTR(3U, c6_t_Kalman);
          _SFD_SET_DATA_VALUE_PTR(4U, *c6_p);
          _SFD_SET_DATA_VALUE_PTR(5U, *c6_M);
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
  return "2bLvARcyFLfB8CMWpJ0uOE";
}

static void sf_opaque_initialize_c6_Model_02(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc6_Model_02InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c6_Model_02((SFc6_Model_02InstanceStruct*) chartInstanceVar);
  initialize_c6_Model_02((SFc6_Model_02InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c6_Model_02(void *chartInstanceVar)
{
  enable_c6_Model_02((SFc6_Model_02InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c6_Model_02(void *chartInstanceVar)
{
  disable_c6_Model_02((SFc6_Model_02InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c6_Model_02(void *chartInstanceVar)
{
  sf_gateway_c6_Model_02((SFc6_Model_02InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c6_Model_02(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c6_Model_02((SFc6_Model_02InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c6_Model_02();/* state var info */
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

extern void sf_internal_set_sim_state_c6_Model_02(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c6_Model_02();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c6_Model_02((SFc6_Model_02InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c6_Model_02(SimStruct* S)
{
  return sf_internal_get_sim_state_c6_Model_02(S);
}

static void sf_opaque_set_sim_state_c6_Model_02(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c6_Model_02(S, st);
}

static void sf_opaque_terminate_c6_Model_02(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc6_Model_02InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_02_optimization_info();
    }

    finalize_c6_Model_02((SFc6_Model_02InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc6_Model_02((SFc6_Model_02InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c6_Model_02(SimStruct *S)
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
    initialize_params_c6_Model_02((SFc6_Model_02InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c6_Model_02(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_02_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,6);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,6,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,6,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,6);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,6,5);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,6,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 5; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,6);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3984997255U));
  ssSetChecksum1(S,(476736259U));
  ssSetChecksum2(S,(3004143715U));
  ssSetChecksum3(S,(626950274U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c6_Model_02(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c6_Model_02(SimStruct *S)
{
  SFc6_Model_02InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc6_Model_02InstanceStruct *)utMalloc(sizeof
    (SFc6_Model_02InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc6_Model_02InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c6_Model_02;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c6_Model_02;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c6_Model_02;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c6_Model_02;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c6_Model_02;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c6_Model_02;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c6_Model_02;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c6_Model_02;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c6_Model_02;
  chartInstance->chartInfo.mdlStart = mdlStart_c6_Model_02;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c6_Model_02;
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

void c6_Model_02_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c6_Model_02(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c6_Model_02(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c6_Model_02(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c6_Model_02_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

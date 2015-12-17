/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_01_sfun.h"
#include "c5_Model_01.h"
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
static const char * c5_debug_family_names[15] = { "Omega_Tensor", "omega_norm",
  "v_lambda", "A", "Gamma", "i", "j", "k", "nargin", "nargout", "omega",
  "t_delta", "M", "N", "Phi_r23" };

static const char * c5_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

static const char * c5_c_debug_family_names[6] = { "nargin", "nargout", "lambda",
  "omega_norm", "t_delta", "phi_j3" };

static const char * c5_d_debug_family_names[6] = { "nargin", "nargout", "lambda",
  "omega_norm", "t_delta", "phi_j2" };

static const char * c5_e_debug_family_names[5] = { "nargin", "nargout", "lambda",
  "t_delta", "phi_j1" };

static const char * c5_f_debug_family_names[7] = { "nargin", "nargout", "k",
  "lambda", "omega_norm", "t_delta", "phi_jk" };

/* Function Declarations */
static void initialize_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance);
static void initialize_params_c5_Model_01(SFc5_Model_01InstanceStruct
  *chartInstance);
static void enable_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance);
static void disable_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_update_debugger_state_c5_Model_01(SFc5_Model_01InstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c5_Model_01(SFc5_Model_01InstanceStruct
  *chartInstance);
static void set_sim_state_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_st);
static void finalize_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance);
static void sf_gateway_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_chartstep_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance);
static void initSimStructsc5_Model_01(SFc5_Model_01InstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber, uint32_T c5_instanceNumber);
static const mxArray *c5_sf_marshallOut(void *chartInstanceVoid, void *c5_inData);
static void c5_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_Phi_r23, const char_T *c5_identifier, real_T c5_y[9]);
static void c5_b_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[9]);
static void c5_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static const mxArray *c5_b_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static const mxArray *c5_c_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static real_T c5_c_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static void c5_d_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[3]);
static void c5_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static void c5_info_helper(const mxArray **c5_info);
static const mxArray *c5_emlrt_marshallOut(const char * c5_u);
static const mxArray *c5_b_emlrt_marshallOut(const uint32_T c5_u);
static void c5_b_info_helper(const mxArray **c5_info);
static void c5_c_info_helper(const mxArray **c5_info);
static void c5_d_info_helper(const mxArray **c5_info);
static void c5_e_info_helper(const mxArray **c5_info);
static real_T c5_norm(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x[3]);
static void c5_eml_switch_helper(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_realmin(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_b_eml_switch_helper(SFc5_Model_01InstanceStruct *chartInstance);
static real_T c5_abs(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x);
static void c5_eig(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_A[9],
                   creal_T c5_V[3]);
static void c5_eml_error(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_eps(SFc5_Model_01InstanceStruct *chartInstance);
static real_T c5_eml_matlab_zlangeM(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_x[9]);
static real_T c5_b_abs(SFc5_Model_01InstanceStruct *chartInstance, creal_T c5_x);
static boolean_T c5_isfinite(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_x);
static void c5_eml_matlab_zlascl(SFc5_Model_01InstanceStruct *chartInstance,
  real_T c5_cfrom, real_T c5_cto, creal_T c5_A[9], creal_T c5_b_A[9]);
static void c5_eml_matlab_zggbal(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], creal_T c5_b_A[9], int32_T *c5_ilo, int32_T *c5_ihi, int32_T
  c5_rscale[3]);
static void c5_check_forloop_overflow_error(SFc5_Model_01InstanceStruct
  *chartInstance, boolean_T c5_overflow);
static void c5_eml_matlab_zlartg(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_f, creal_T c5_g, real_T *c5_cs, creal_T *c5_sn, creal_T *c5_r);
static void c5_eml_scalar_eg(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_eml_matlab_zhgeqz(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T c5_ilo, int32_T c5_ihi, real_T *c5_info, creal_T
  c5_alpha1[3], creal_T c5_beta1[3]);
static real_T c5_eml_matlab_zlanhs(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T c5_ilo, int32_T c5_ihi);
static int32_T c5_mod(SFc5_Model_01InstanceStruct *chartInstance, int32_T c5_x);
static creal_T c5_eml_div(SFc5_Model_01InstanceStruct *chartInstance, creal_T
  c5_x, real_T c5_y);
static void c5_scalarEg(SFc5_Model_01InstanceStruct *chartInstance);
static creal_T c5_sqrt(SFc5_Model_01InstanceStruct *chartInstance, creal_T c5_x);
static void c5_realmax(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_b_eml_matlab_zlartg(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_f, creal_T c5_g, real_T *c5_cs, creal_T *c5_sn);
static void c5_b_eml_matlab_zlascl(SFc5_Model_01InstanceStruct *chartInstance,
  real_T c5_cfrom, real_T c5_cto, creal_T c5_A[3], creal_T c5_b_A[3]);
static void c5_b_eml_div(SFc5_Model_01InstanceStruct *chartInstance, creal_T
  c5_x[3], creal_T c5_y[3], creal_T c5_z[3]);
static void c5_eml_warning(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_b_eml_warning(SFc5_Model_01InstanceStruct *chartInstance);
static real_T c5_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_a);
static void c5_inv(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x[9],
                   real_T c5_y[9]);
static void c5_inv3x3(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x[9],
                      real_T c5_y[9]);
static real_T c5_b_norm(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x
  [9]);
static void c5_c_eml_warning(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_d_eml_warning(SFc5_Model_01InstanceStruct *chartInstance, char_T
  c5_varargin_2[14]);
static real_T c5_b_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_a);
static real_T c5_c_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_a);
static real_T c5_d_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_a);
static real_T c5_e_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_a);
static void c5_f_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_a
  [9], real_T c5_b, real_T c5_c[9]);
static void c5_b_eml_scalar_eg(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_c_eml_scalar_eg(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_eml_xgemm(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_A[9], real_T c5_B[9], real_T c5_C[9], real_T c5_b_C[9]);
static void c5_b_eml_error(SFc5_Model_01InstanceStruct *chartInstance);
static void c5_b_eig(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_A[9],
                     creal_T c5_V[9], creal_T c5_D[9]);
static void c5_eml_matlab_zggev(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], real_T *c5_info, creal_T c5_alpha1[3], creal_T c5_beta1[3],
  creal_T c5_V[9]);
static void c5_eml_matlab_zgghrd(SFc5_Model_01InstanceStruct *chartInstance,
  int32_T c5_ilo, int32_T c5_ihi, creal_T c5_A[9], creal_T c5_b_A[9], creal_T
  c5_Z[9]);
static void c5_b_eml_matlab_zhgeqz(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T c5_ilo, int32_T c5_ihi, creal_T c5_Z[9], real_T
  *c5_info, creal_T c5_alpha1[3], creal_T c5_beta1[3], creal_T c5_b_A[9],
  creal_T c5_b_Z[9]);
static void c5_eml_matlab_ztgevc(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], creal_T c5_V[9], creal_T c5_b_V[9]);
static creal_T c5_rdivide(SFc5_Model_01InstanceStruct *chartInstance, creal_T
  c5_x, creal_T c5_y);
static creal_T c5_power(SFc5_Model_01InstanceStruct *chartInstance, creal_T c5_a,
  real_T c5_b);
static void c5_eml_lusolve(SFc5_Model_01InstanceStruct *chartInstance, creal_T
  c5_A[9], creal_T c5_B[9], creal_T c5_X[9]);
static void c5_e_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_sprintf, const char_T *c5_identifier, char_T c5_y[14]);
static void c5_f_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, char_T c5_y[14]);
static const mxArray *c5_d_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static int32_T c5_g_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static uint8_T c5_h_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_b_is_active_c5_Model_01, const char_T *c5_identifier);
static uint8_T c5_i_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_c_eml_matlab_zlascl(SFc5_Model_01InstanceStruct *chartInstance,
  real_T c5_cfrom, real_T c5_cto, creal_T c5_A[9]);
static void c5_b_eml_matlab_zggbal(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T *c5_ilo, int32_T *c5_ihi, int32_T c5_rscale[3]);
static void c5_b_sqrt(SFc5_Model_01InstanceStruct *chartInstance, creal_T *c5_x);
static void c5_d_eml_matlab_zlascl(SFc5_Model_01InstanceStruct *chartInstance,
  real_T c5_cfrom, real_T c5_cto, creal_T c5_A[3]);
static void c5_b_eml_xgemm(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_A[9], real_T c5_B[9], real_T c5_C[9]);
static void c5_b_eml_matlab_zgghrd(SFc5_Model_01InstanceStruct *chartInstance,
  int32_T c5_ilo, int32_T c5_ihi, creal_T c5_A[9], creal_T c5_Z[9]);
static void c5_c_eml_matlab_zhgeqz(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T c5_ilo, int32_T c5_ihi, creal_T c5_Z[9], real_T
  *c5_info, creal_T c5_alpha1[3], creal_T c5_beta1[3]);
static void c5_b_eml_matlab_ztgevc(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], creal_T c5_V[9]);
static int32_T c5_div_s32(SFc5_Model_01InstanceStruct *chartInstance, int32_T
  c5_numerator, int32_T c5_denominator);
static void init_dsm_address_info(SFc5_Model_01InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance)
{
  chartInstance->c5_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c5_is_active_c5_Model_01 = 0U;
}

static void initialize_params_c5_Model_01(SFc5_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c5_update_debugger_state_c5_Model_01(SFc5_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c5_Model_01(SFc5_Model_01InstanceStruct
  *chartInstance)
{
  const mxArray *c5_st;
  const mxArray *c5_y = NULL;
  int32_T c5_i0;
  real_T c5_u[9];
  const mxArray *c5_b_y = NULL;
  uint8_T c5_hoistedGlobal;
  uint8_T c5_b_u;
  const mxArray *c5_c_y = NULL;
  real_T (*c5_Phi_r23)[9];
  c5_Phi_r23 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c5_st = NULL;
  c5_st = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_createcellmatrix(2, 1), false);
  for (c5_i0 = 0; c5_i0 < 9; c5_i0++) {
    c5_u[c5_i0] = (*c5_Phi_r23)[c5_i0];
  }

  c5_b_y = NULL;
  sf_mex_assign(&c5_b_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_setcell(c5_y, 0, c5_b_y);
  c5_hoistedGlobal = chartInstance->c5_is_active_c5_Model_01;
  c5_b_u = c5_hoistedGlobal;
  c5_c_y = NULL;
  sf_mex_assign(&c5_c_y, sf_mex_create("y", &c5_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c5_y, 1, c5_c_y);
  sf_mex_assign(&c5_st, c5_y, false);
  return c5_st;
}

static void set_sim_state_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_st)
{
  const mxArray *c5_u;
  real_T c5_dv0[9];
  int32_T c5_i1;
  real_T (*c5_Phi_r23)[9];
  c5_Phi_r23 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c5_doneDoubleBufferReInit = true;
  c5_u = sf_mex_dup(c5_st);
  c5_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 0)),
                      "Phi_r23", c5_dv0);
  for (c5_i1 = 0; c5_i1 < 9; c5_i1++) {
    (*c5_Phi_r23)[c5_i1] = c5_dv0[c5_i1];
  }

  chartInstance->c5_is_active_c5_Model_01 = c5_h_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c5_u, 1)), "is_active_c5_Model_01");
  sf_mex_destroy(&c5_u);
  c5_update_debugger_state_c5_Model_01(chartInstance);
  sf_mex_destroy(&c5_st);
}

static void finalize_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance)
{
  int32_T c5_i2;
  int32_T c5_i3;
  int32_T c5_i4;
  int32_T c5_i5;
  real_T *c5_t_delta;
  real_T (*c5_N)[9];
  real_T (*c5_M)[9];
  real_T (*c5_Phi_r23)[9];
  real_T (*c5_omega)[3];
  c5_N = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
  c5_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 2);
  c5_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c5_Phi_r23 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c5_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
  for (c5_i2 = 0; c5_i2 < 3; c5_i2++) {
    _SFD_DATA_RANGE_CHECK((*c5_omega)[c5_i2], 0U);
  }

  chartInstance->c5_sfEvent = CALL_EVENT;
  c5_chartstep_c5_Model_01(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_01MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c5_i3 = 0; c5_i3 < 9; c5_i3++) {
    _SFD_DATA_RANGE_CHECK((*c5_Phi_r23)[c5_i3], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c5_t_delta, 2U);
  for (c5_i4 = 0; c5_i4 < 9; c5_i4++) {
    _SFD_DATA_RANGE_CHECK((*c5_M)[c5_i4], 3U);
  }

  for (c5_i5 = 0; c5_i5 < 9; c5_i5++) {
    _SFD_DATA_RANGE_CHECK((*c5_N)[c5_i5], 4U);
  }
}

static void c5_chartstep_c5_Model_01(SFc5_Model_01InstanceStruct *chartInstance)
{
  real_T c5_hoistedGlobal;
  int32_T c5_i6;
  real_T c5_omega[3];
  real_T c5_t_delta;
  int32_T c5_i7;
  real_T c5_M[9];
  int32_T c5_i8;
  real_T c5_N[9];
  uint32_T c5_debug_family_var_map[15];
  real_T c5_Omega_Tensor[9];
  real_T c5_omega_norm;
  real_T c5_v_lambda[3];
  real_T c5_A[9];
  real_T c5_Gamma[9];
  real_T c5_i;
  real_T c5_j;
  real_T c5_k;
  real_T c5_nargin = 4.0;
  real_T c5_nargout = 1.0;
  real_T c5_Phi_r23[9];
  int32_T c5_i9;
  int32_T c5_i10;
  real_T c5_v[3];
  uint32_T c5_b_debug_family_var_map[4];
  real_T c5_b_nargin = 1.0;
  real_T c5_b_nargout = 1.0;
  int32_T c5_i11;
  int32_T c5_i12;
  real_T c5_b_omega[3];
  int32_T c5_i13;
  real_T c5_b_M[9];
  creal_T c5_dcv0[3];
  int32_T c5_i14;
  real_T c5_d0;
  real_T c5_d1;
  real_T c5_d2;
  int32_T c5_i15;
  real_T c5_b_A[9];
  real_T c5_dv1[9];
  int32_T c5_i16;
  int32_T c5_b_i;
  int32_T c5_b_j;
  int32_T c5_b_k;
  real_T c5_c_k;
  real_T c5_lambda;
  real_T c5_b_omega_norm;
  real_T c5_b_t_delta;
  uint32_T c5_c_debug_family_var_map[7];
  real_T c5_c_nargin = 4.0;
  real_T c5_c_nargout = 1.0;
  real_T c5_phi_jk;
  real_T c5_b_lambda;
  real_T c5_c_t_delta;
  uint32_T c5_d_debug_family_var_map[5];
  real_T c5_d_nargin = 2.0;
  real_T c5_d_nargout = 1.0;
  real_T c5_x;
  real_T c5_b_x;
  real_T c5_c_lambda;
  real_T c5_c_omega_norm;
  real_T c5_d_t_delta;
  uint32_T c5_e_debug_family_var_map[6];
  real_T c5_e_nargin = 3.0;
  real_T c5_e_nargout = 1.0;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_e_x;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_h_x;
  real_T c5_d_lambda;
  real_T c5_d_omega_norm;
  real_T c5_e_t_delta;
  real_T c5_f_nargin = 3.0;
  real_T c5_f_nargout = 1.0;
  real_T c5_a;
  real_T c5_b_a;
  real_T c5_c_a;
  real_T c5_ak;
  real_T c5_d_a;
  real_T c5_ar;
  real_T c5_c;
  real_T c5_i_x;
  real_T c5_j_x;
  real_T c5_k_x;
  real_T c5_l_x;
  real_T c5_m_x;
  real_T c5_n_x;
  real_T c5_e_a;
  int32_T c5_i17;
  real_T c5_b_Omega_Tensor[9];
  real_T c5_b[9];
  int32_T c5_i18;
  int32_T c5_i19;
  real_T c5_c_M[9];
  real_T c5_b_b[9];
  int32_T c5_i20;
  real_T c5_y[9];
  int32_T c5_i21;
  real_T c5_c_b[9];
  int32_T c5_i22;
  real_T c5_d_b[9];
  int32_T c5_i23;
  int32_T c5_i24;
  int32_T c5_i25;
  real_T c5_b_y[9];
  int32_T c5_i26;
  real_T c5_e_b[9];
  int32_T c5_i27;
  int32_T c5_i28;
  real_T *c5_f_t_delta;
  real_T (*c5_b_Phi_r23)[9];
  real_T (*c5_b_N)[9];
  real_T (*c5_d_M)[9];
  real_T (*c5_c_omega)[3];
  c5_b_N = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
  c5_d_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 2);
  c5_f_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c5_b_Phi_r23 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c5_c_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
  c5_hoistedGlobal = *c5_f_t_delta;
  for (c5_i6 = 0; c5_i6 < 3; c5_i6++) {
    c5_omega[c5_i6] = (*c5_c_omega)[c5_i6];
  }

  c5_t_delta = c5_hoistedGlobal;
  for (c5_i7 = 0; c5_i7 < 9; c5_i7++) {
    c5_M[c5_i7] = (*c5_d_M)[c5_i7];
  }

  for (c5_i8 = 0; c5_i8 < 9; c5_i8++) {
    c5_N[c5_i8] = (*c5_b_N)[c5_i8];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 15U, 15U, c5_debug_family_names,
    c5_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_Omega_Tensor, 0U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_omega_norm, 1U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_v_lambda, 2U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_A, 3U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_Gamma, 4U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_i, 5U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_j, 6U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_k, 7U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_nargin, 8U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_nargout, 9U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c5_omega, 10U, c5_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_t_delta, 11U, c5_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c5_M, 12U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c5_N, 13U, c5_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_Phi_r23, 14U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 3);
  for (c5_i9 = 0; c5_i9 < 9; c5_i9++) {
    c5_Phi_r23[c5_i9] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 4);
  for (c5_i10 = 0; c5_i10 < 3; c5_i10++) {
    c5_v[c5_i10] = c5_omega[c5_i10];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c5_b_debug_family_names,
    c5_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_b_nargin, 0U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_b_nargout, 1U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_v, 2U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_Omega_Tensor, 3U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 2);
  for (c5_i11 = 0; c5_i11 < 9; c5_i11++) {
    c5_Omega_Tensor[c5_i11] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 3);
  c5_Omega_Tensor[0] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 4);
  c5_Omega_Tensor[4] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 5);
  c5_Omega_Tensor[8] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 6);
  c5_Omega_Tensor[3] = -c5_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 7);
  c5_Omega_Tensor[6] = c5_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 8);
  c5_Omega_Tensor[7] = -c5_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 9);
  c5_Omega_Tensor[1] = c5_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 10);
  c5_Omega_Tensor[2] = -c5_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, 11);
  c5_Omega_Tensor[5] = c5_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c5_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 5);
  for (c5_i12 = 0; c5_i12 < 3; c5_i12++) {
    c5_b_omega[c5_i12] = c5_omega[c5_i12];
  }

  c5_omega_norm = c5_norm(chartInstance, c5_b_omega);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 6);
  for (c5_i13 = 0; c5_i13 < 9; c5_i13++) {
    c5_b_M[c5_i13] = c5_M[c5_i13];
  }

  c5_eig(chartInstance, c5_b_M, c5_dcv0);
  for (c5_i14 = 0; c5_i14 < 3; c5_i14++) {
    c5_v_lambda[c5_i14] = c5_dcv0[c5_i14].re;
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 7);
  c5_d0 = c5_mpower(chartInstance, c5_v_lambda[0]);
  c5_d1 = c5_mpower(chartInstance, c5_v_lambda[1]);
  c5_d2 = c5_mpower(chartInstance, c5_v_lambda[2]);
  c5_A[0] = 1.0;
  c5_A[3] = c5_v_lambda[0];
  c5_A[6] = c5_d0;
  c5_A[1] = 1.0;
  c5_A[4] = c5_v_lambda[1];
  c5_A[7] = c5_d1;
  c5_A[2] = 1.0;
  c5_A[5] = c5_v_lambda[2];
  c5_A[8] = c5_d2;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 8);
  for (c5_i15 = 0; c5_i15 < 9; c5_i15++) {
    c5_b_A[c5_i15] = c5_A[c5_i15];
  }

  c5_inv(chartInstance, c5_b_A, c5_dv1);
  for (c5_i16 = 0; c5_i16 < 9; c5_i16++) {
    c5_Gamma[c5_i16] = c5_dv1[c5_i16];
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 9);
  c5_i = 1.0;
  c5_b_i = 0;
  while (c5_b_i < 3) {
    c5_i = 1.0 + (real_T)c5_b_i;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 10);
    c5_j = 1.0;
    c5_b_j = 0;
    while (c5_b_j < 3) {
      c5_j = 1.0 + (real_T)c5_b_j;
      CV_EML_FOR(0, 1, 1, 1);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 11);
      c5_k = 1.0;
      c5_b_k = 0;
      while (c5_b_k < 3) {
        c5_k = 1.0 + (real_T)c5_b_k;
        CV_EML_FOR(0, 1, 2, 1);
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 12);
        c5_c_k = c5_k;
        c5_lambda = c5_v_lambda[_SFD_EML_ARRAY_BOUNDS_CHECK("v_lambda", (int32_T)
          _SFD_INTEGER_CHECK("j", c5_j), 1, 3, 1, 0) - 1];
        c5_b_omega_norm = c5_omega_norm;
        c5_b_t_delta = c5_t_delta;
        _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 7U, 7U, c5_f_debug_family_names,
          c5_c_debug_family_var_map);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_c_nargin, 0U,
          c5_b_sf_marshallOut, c5_b_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_c_nargout, 1U,
          c5_b_sf_marshallOut, c5_b_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_c_k, 2U, c5_b_sf_marshallOut,
          c5_b_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_lambda, 3U, c5_b_sf_marshallOut,
          c5_b_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_b_omega_norm, 4U,
          c5_b_sf_marshallOut, c5_b_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_b_t_delta, 5U,
          c5_b_sf_marshallOut, c5_b_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_phi_jk, 6U, c5_b_sf_marshallOut,
          c5_b_sf_marshallIn);
        CV_EML_FCN(0, 1);
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 18);
        switch ((int32_T)_SFD_INTEGER_CHECK("k", c5_c_k)) {
         case 1:
          CV_EML_SWITCH(0, 1, 0, 1);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 20);
          c5_b_lambda = c5_lambda;
          c5_c_t_delta = c5_b_t_delta;
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 5U, 5U, c5_e_debug_family_names,
            c5_d_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_d_nargin, 0U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_d_nargout, 1U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_b_lambda, 2U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_c_t_delta, 3U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_phi_jk, 4U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          CV_EML_FCN(0, 2);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 31);
          c5_x = c5_b_lambda * c5_c_t_delta;
          c5_b_x = c5_x;
          c5_b_x = muDoubleScalarExp(c5_b_x);
          c5_phi_jk = c5_b_mpower(chartInstance, c5_b_lambda) * (1.0 - c5_b_x) -
            c5_e_mpower(chartInstance, c5_b_lambda) * c5_c_t_delta;
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, -31);
          _SFD_SYMBOL_SCOPE_POP();
          break;

         case 2:
          CV_EML_SWITCH(0, 1, 0, 2);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 22);
          c5_c_lambda = c5_lambda;
          c5_c_omega_norm = c5_b_omega_norm;
          c5_d_t_delta = c5_b_t_delta;
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c5_d_debug_family_names,
            c5_e_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_e_nargin, 0U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_e_nargout, 1U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_c_lambda, 2U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_c_omega_norm, 3U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_d_t_delta, 4U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_phi_jk, 5U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          CV_EML_FCN(0, 3);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 35);
          c5_c_x = c5_c_omega_norm * c5_d_t_delta;
          c5_d_x = c5_c_x;
          c5_d_x = muDoubleScalarSin(c5_d_x);
          c5_e_x = c5_c_omega_norm * c5_d_t_delta;
          c5_f_x = c5_e_x;
          c5_f_x = muDoubleScalarCos(c5_f_x);
          c5_g_x = c5_c_omega_norm * c5_d_t_delta;
          c5_h_x = c5_g_x;
          c5_h_x = muDoubleScalarExp(c5_h_x);
          c5_phi_jk = (c5_d_mpower(chartInstance, c5_c_lambda) * c5_mpower
                       (chartInstance, c5_c_omega_norm) + c5_e_mpower
                       (chartInstance, c5_c_lambda * c5_c_mpower(chartInstance,
            c5_c_omega_norm))) * ((((c5_c_omega_norm * c5_c_lambda * c5_d_x +
            c5_mpower(chartInstance, c5_c_lambda) * c5_f_x) + c5_mpower
            (chartInstance, c5_c_omega_norm) * c5_h_x) - c5_mpower(chartInstance,
            c5_c_lambda)) - c5_mpower(chartInstance, c5_c_omega_norm));
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, -35);
          _SFD_SYMBOL_SCOPE_POP();
          break;

         case 3:
          CV_EML_SWITCH(0, 1, 0, 3);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 24);
          c5_d_lambda = c5_lambda;
          c5_d_omega_norm = c5_b_omega_norm;
          c5_e_t_delta = c5_b_t_delta;
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c5_c_debug_family_names,
            c5_e_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_f_nargin, 0U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_f_nargout, 1U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_d_lambda, 2U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_d_omega_norm, 3U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_e_t_delta, 4U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_phi_jk, 5U,
            c5_b_sf_marshallOut, c5_b_sf_marshallIn);
          CV_EML_FCN(0, 4);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 39);
          c5_a = c5_d_omega_norm;
          c5_b_a = c5_a;
          c5_c_a = c5_b_a;
          c5_eml_scalar_eg(chartInstance);
          c5_ak = c5_c_a;
          c5_d_a = c5_ak;
          c5_eml_scalar_eg(chartInstance);
          c5_ar = c5_d_a;
          c5_c = muDoubleScalarPower(c5_ar, 5.0);
          c5_i_x = c5_d_omega_norm * c5_e_t_delta;
          c5_j_x = c5_i_x;
          c5_j_x = muDoubleScalarSin(c5_j_x);
          c5_k_x = c5_d_omega_norm * c5_e_t_delta;
          c5_l_x = c5_k_x;
          c5_l_x = muDoubleScalarCos(c5_l_x);
          c5_m_x = c5_d_omega_norm * c5_e_t_delta;
          c5_n_x = c5_m_x;
          c5_n_x = muDoubleScalarExp(c5_n_x);
          c5_phi_jk = c5_b_mpower(chartInstance, c5_d_omega_norm) * c5_b_mpower
            (chartInstance, c5_d_lambda) + c5_e_mpower(chartInstance,
            c5_c_mpower(chartInstance, c5_d_lambda) * c5_d_mpower(chartInstance,
            c5_d_omega_norm) + c5_mpower(chartInstance, c5_d_lambda) * c5_c) *
            (((c5_d_mpower(chartInstance, c5_d_lambda) * c5_j_x -
               c5_d_omega_norm * c5_mpower(chartInstance, c5_d_lambda) * c5_l_x)
              - c5_d_mpower(chartInstance, c5_d_omega_norm) * c5_n_x) -
             (c5_d_mpower(chartInstance, c5_d_omega_norm) * c5_d_lambda +
              c5_d_mpower(chartInstance, c5_d_lambda) * c5_d_omega_norm) *
             c5_e_t_delta);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, -39);
          _SFD_SYMBOL_SCOPE_POP();
          break;

         default:
          CV_EML_SWITCH(0, 1, 0, 0);
          _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 26);
          c5_phi_jk = 0.0;
          break;
        }

        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, -26);
        _SFD_SYMBOL_SCOPE_POP();
        c5_e_a = c5_Gamma[(_SFD_EML_ARRAY_BOUNDS_CHECK("Gamma", (int32_T)
          _SFD_INTEGER_CHECK("i", c5_i), 1, 3, 1, 0) + 3 *
                           (_SFD_EML_ARRAY_BOUNDS_CHECK("Gamma", (int32_T)
          _SFD_INTEGER_CHECK("j", c5_j), 1, 3, 2, 0) - 1)) - 1] * c5_phi_jk;
        for (c5_i17 = 0; c5_i17 < 9; c5_i17++) {
          c5_b_Omega_Tensor[c5_i17] = c5_Omega_Tensor[c5_i17];
        }

        c5_f_mpower(chartInstance, c5_b_Omega_Tensor, c5_k - 1.0, c5_b);
        for (c5_i18 = 0; c5_i18 < 9; c5_i18++) {
          c5_b[c5_i18] *= c5_e_a;
        }

        for (c5_i19 = 0; c5_i19 < 9; c5_i19++) {
          c5_c_M[c5_i19] = c5_M[c5_i19];
        }

        c5_f_mpower(chartInstance, c5_c_M, c5_i - 1.0, c5_b_b);
        c5_c_eml_scalar_eg(chartInstance);
        c5_c_eml_scalar_eg(chartInstance);
        for (c5_i20 = 0; c5_i20 < 9; c5_i20++) {
          c5_y[c5_i20] = 0.0;
        }

        for (c5_i21 = 0; c5_i21 < 9; c5_i21++) {
          c5_c_b[c5_i21] = c5_b[c5_i21];
        }

        for (c5_i22 = 0; c5_i22 < 9; c5_i22++) {
          c5_d_b[c5_i22] = c5_b_b[c5_i22];
        }

        c5_b_eml_xgemm(chartInstance, c5_c_b, c5_d_b, c5_y);
        for (c5_i23 = 0; c5_i23 < 9; c5_i23++) {
          c5_b[c5_i23] = c5_N[c5_i23];
        }

        c5_c_eml_scalar_eg(chartInstance);
        c5_c_eml_scalar_eg(chartInstance);
        for (c5_i24 = 0; c5_i24 < 9; c5_i24++) {
          c5_b_b[c5_i24] = 0.0;
        }

        for (c5_i25 = 0; c5_i25 < 9; c5_i25++) {
          c5_b_y[c5_i25] = c5_y[c5_i25];
        }

        for (c5_i26 = 0; c5_i26 < 9; c5_i26++) {
          c5_e_b[c5_i26] = c5_b[c5_i26];
        }

        c5_b_eml_xgemm(chartInstance, c5_b_y, c5_e_b, c5_b_b);
        for (c5_i27 = 0; c5_i27 < 9; c5_i27++) {
          c5_Phi_r23[c5_i27] += c5_b_b[c5_i27];
        }

        c5_b_k++;
        _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
      }

      CV_EML_FOR(0, 1, 2, 0);
      c5_b_j++;
      _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
    }

    CV_EML_FOR(0, 1, 1, 0);
    c5_b_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, -12);
  _SFD_SYMBOL_SCOPE_POP();
  for (c5_i28 = 0; c5_i28 < 9; c5_i28++) {
    (*c5_b_Phi_r23)[c5_i28] = c5_Phi_r23[c5_i28];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
}

static void initSimStructsc5_Model_01(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber, uint32_T c5_instanceNumber)
{
  (void)c5_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c5_chartNumber, c5_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c5_sf_marshallOut(void *chartInstanceVoid, void *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_i29;
  int32_T c5_i30;
  int32_T c5_i31;
  real_T c5_b_inData[9];
  int32_T c5_i32;
  int32_T c5_i33;
  int32_T c5_i34;
  real_T c5_u[9];
  const mxArray *c5_y = NULL;
  SFc5_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc5_Model_01InstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_i29 = 0;
  for (c5_i30 = 0; c5_i30 < 3; c5_i30++) {
    for (c5_i31 = 0; c5_i31 < 3; c5_i31++) {
      c5_b_inData[c5_i31 + c5_i29] = (*(real_T (*)[9])c5_inData)[c5_i31 + c5_i29];
    }

    c5_i29 += 3;
  }

  c5_i32 = 0;
  for (c5_i33 = 0; c5_i33 < 3; c5_i33++) {
    for (c5_i34 = 0; c5_i34 < 3; c5_i34++) {
      c5_u[c5_i34 + c5_i32] = c5_b_inData[c5_i34 + c5_i32];
    }

    c5_i32 += 3;
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static void c5_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_Phi_r23, const char_T *c5_identifier, real_T c5_y[9])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_Phi_r23), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_Phi_r23);
}

static void c5_b_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[9])
{
  real_T c5_dv2[9];
  int32_T c5_i35;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv2, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c5_i35 = 0; c5_i35 < 9; c5_i35++) {
    c5_y[c5_i35] = c5_dv2[c5_i35];
  }

  sf_mex_destroy(&c5_u);
}

static void c5_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_Phi_r23;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y[9];
  int32_T c5_i36;
  int32_T c5_i37;
  int32_T c5_i38;
  SFc5_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc5_Model_01InstanceStruct *)chartInstanceVoid;
  c5_Phi_r23 = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_Phi_r23), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_Phi_r23);
  c5_i36 = 0;
  for (c5_i37 = 0; c5_i37 < 3; c5_i37++) {
    for (c5_i38 = 0; c5_i38 < 3; c5_i38++) {
      (*(real_T (*)[9])c5_outData)[c5_i38 + c5_i36] = c5_y[c5_i38 + c5_i36];
    }

    c5_i36 += 3;
  }

  sf_mex_destroy(&c5_mxArrayInData);
}

static const mxArray *c5_b_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  real_T c5_u;
  const mxArray *c5_y = NULL;
  SFc5_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc5_Model_01InstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_u = *(real_T *)c5_inData;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static const mxArray *c5_c_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_i39;
  real_T c5_b_inData[3];
  int32_T c5_i40;
  real_T c5_u[3];
  const mxArray *c5_y = NULL;
  SFc5_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc5_Model_01InstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  for (c5_i39 = 0; c5_i39 < 3; c5_i39++) {
    c5_b_inData[c5_i39] = (*(real_T (*)[3])c5_inData)[c5_i39];
  }

  for (c5_i40 = 0; c5_i40 < 3; c5_i40++) {
    c5_u[c5_i40] = c5_b_inData[c5_i40];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static real_T c5_c_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  real_T c5_y;
  real_T c5_d3;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_d3, 1, 0, 0U, 0, 0U, 0);
  c5_y = c5_d3;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_nargout;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y;
  SFc5_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc5_Model_01InstanceStruct *)chartInstanceVoid;
  c5_nargout = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_nargout), &c5_thisId);
  sf_mex_destroy(&c5_nargout);
  *(real_T *)c5_outData = c5_y;
  sf_mex_destroy(&c5_mxArrayInData);
}

static void c5_d_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[3])
{
  real_T c5_dv3[3];
  int32_T c5_i41;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv3, 1, 0, 0U, 1, 0U, 1, 3);
  for (c5_i41 = 0; c5_i41 < 3; c5_i41++) {
    c5_y[c5_i41] = c5_dv3[c5_i41];
  }

  sf_mex_destroy(&c5_u);
}

static void c5_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_v_lambda;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y[3];
  int32_T c5_i42;
  SFc5_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc5_Model_01InstanceStruct *)chartInstanceVoid;
  c5_v_lambda = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_v_lambda), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_v_lambda);
  for (c5_i42 = 0; c5_i42 < 3; c5_i42++) {
    (*(real_T (*)[3])c5_outData)[c5_i42] = c5_y[c5_i42];
  }

  sf_mex_destroy(&c5_mxArrayInData);
}

const mxArray *sf_c5_Model_01_get_eml_resolved_functions_info(void)
{
  const mxArray *c5_nameCaptureInfo = NULL;
  c5_nameCaptureInfo = NULL;
  sf_mex_assign(&c5_nameCaptureInfo, sf_mex_createstruct("structure", 2, 280, 1),
                false);
  c5_info_helper(&c5_nameCaptureInfo);
  c5_b_info_helper(&c5_nameCaptureInfo);
  c5_c_info_helper(&c5_nameCaptureInfo);
  c5_d_info_helper(&c5_nameCaptureInfo);
  c5_e_info_helper(&c5_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c5_nameCaptureInfo);
  return c5_nameCaptureInfo;
}

static void c5_info_helper(const mxArray **c5_info)
{
  const mxArray *c5_rhs0 = NULL;
  const mxArray *c5_lhs0 = NULL;
  const mxArray *c5_rhs1 = NULL;
  const mxArray *c5_lhs1 = NULL;
  const mxArray *c5_rhs2 = NULL;
  const mxArray *c5_lhs2 = NULL;
  const mxArray *c5_rhs3 = NULL;
  const mxArray *c5_lhs3 = NULL;
  const mxArray *c5_rhs4 = NULL;
  const mxArray *c5_lhs4 = NULL;
  const mxArray *c5_rhs5 = NULL;
  const mxArray *c5_lhs5 = NULL;
  const mxArray *c5_rhs6 = NULL;
  const mxArray *c5_lhs6 = NULL;
  const mxArray *c5_rhs7 = NULL;
  const mxArray *c5_lhs7 = NULL;
  const mxArray *c5_rhs8 = NULL;
  const mxArray *c5_lhs8 = NULL;
  const mxArray *c5_rhs9 = NULL;
  const mxArray *c5_lhs9 = NULL;
  const mxArray *c5_rhs10 = NULL;
  const mxArray *c5_lhs10 = NULL;
  const mxArray *c5_rhs11 = NULL;
  const mxArray *c5_lhs11 = NULL;
  const mxArray *c5_rhs12 = NULL;
  const mxArray *c5_lhs12 = NULL;
  const mxArray *c5_rhs13 = NULL;
  const mxArray *c5_lhs13 = NULL;
  const mxArray *c5_rhs14 = NULL;
  const mxArray *c5_lhs14 = NULL;
  const mxArray *c5_rhs15 = NULL;
  const mxArray *c5_lhs15 = NULL;
  const mxArray *c5_rhs16 = NULL;
  const mxArray *c5_lhs16 = NULL;
  const mxArray *c5_rhs17 = NULL;
  const mxArray *c5_lhs17 = NULL;
  const mxArray *c5_rhs18 = NULL;
  const mxArray *c5_lhs18 = NULL;
  const mxArray *c5_rhs19 = NULL;
  const mxArray *c5_lhs19 = NULL;
  const mxArray *c5_rhs20 = NULL;
  const mxArray *c5_lhs20 = NULL;
  const mxArray *c5_rhs21 = NULL;
  const mxArray *c5_lhs21 = NULL;
  const mxArray *c5_rhs22 = NULL;
  const mxArray *c5_lhs22 = NULL;
  const mxArray *c5_rhs23 = NULL;
  const mxArray *c5_lhs23 = NULL;
  const mxArray *c5_rhs24 = NULL;
  const mxArray *c5_lhs24 = NULL;
  const mxArray *c5_rhs25 = NULL;
  const mxArray *c5_lhs25 = NULL;
  const mxArray *c5_rhs26 = NULL;
  const mxArray *c5_lhs26 = NULL;
  const mxArray *c5_rhs27 = NULL;
  const mxArray *c5_lhs27 = NULL;
  const mxArray *c5_rhs28 = NULL;
  const mxArray *c5_lhs28 = NULL;
  const mxArray *c5_rhs29 = NULL;
  const mxArray *c5_lhs29 = NULL;
  const mxArray *c5_rhs30 = NULL;
  const mxArray *c5_lhs30 = NULL;
  const mxArray *c5_rhs31 = NULL;
  const mxArray *c5_lhs31 = NULL;
  const mxArray *c5_rhs32 = NULL;
  const mxArray *c5_lhs32 = NULL;
  const mxArray *c5_rhs33 = NULL;
  const mxArray *c5_lhs33 = NULL;
  const mxArray *c5_rhs34 = NULL;
  const mxArray *c5_lhs34 = NULL;
  const mxArray *c5_rhs35 = NULL;
  const mxArray *c5_lhs35 = NULL;
  const mxArray *c5_rhs36 = NULL;
  const mxArray *c5_lhs36 = NULL;
  const mxArray *c5_rhs37 = NULL;
  const mxArray *c5_lhs37 = NULL;
  const mxArray *c5_rhs38 = NULL;
  const mxArray *c5_lhs38 = NULL;
  const mxArray *c5_rhs39 = NULL;
  const mxArray *c5_lhs39 = NULL;
  const mxArray *c5_rhs40 = NULL;
  const mxArray *c5_lhs40 = NULL;
  const mxArray *c5_rhs41 = NULL;
  const mxArray *c5_lhs41 = NULL;
  const mxArray *c5_rhs42 = NULL;
  const mxArray *c5_lhs42 = NULL;
  const mxArray *c5_rhs43 = NULL;
  const mxArray *c5_lhs43 = NULL;
  const mxArray *c5_rhs44 = NULL;
  const mxArray *c5_lhs44 = NULL;
  const mxArray *c5_rhs45 = NULL;
  const mxArray *c5_lhs45 = NULL;
  const mxArray *c5_rhs46 = NULL;
  const mxArray *c5_lhs46 = NULL;
  const mxArray *c5_rhs47 = NULL;
  const mxArray *c5_lhs47 = NULL;
  const mxArray *c5_rhs48 = NULL;
  const mxArray *c5_lhs48 = NULL;
  const mxArray *c5_rhs49 = NULL;
  const mxArray *c5_lhs49 = NULL;
  const mxArray *c5_rhs50 = NULL;
  const mxArray *c5_lhs50 = NULL;
  const mxArray *c5_rhs51 = NULL;
  const mxArray *c5_lhs51 = NULL;
  const mxArray *c5_rhs52 = NULL;
  const mxArray *c5_lhs52 = NULL;
  const mxArray *c5_rhs53 = NULL;
  const mxArray *c5_lhs53 = NULL;
  const mxArray *c5_rhs54 = NULL;
  const mxArray *c5_lhs54 = NULL;
  const mxArray *c5_rhs55 = NULL;
  const mxArray *c5_lhs55 = NULL;
  const mxArray *c5_rhs56 = NULL;
  const mxArray *c5_lhs56 = NULL;
  const mxArray *c5_rhs57 = NULL;
  const mxArray *c5_lhs57 = NULL;
  const mxArray *c5_rhs58 = NULL;
  const mxArray *c5_lhs58 = NULL;
  const mxArray *c5_rhs59 = NULL;
  const mxArray *c5_lhs59 = NULL;
  const mxArray *c5_rhs60 = NULL;
  const mxArray *c5_lhs60 = NULL;
  const mxArray *c5_rhs61 = NULL;
  const mxArray *c5_lhs61 = NULL;
  const mxArray *c5_rhs62 = NULL;
  const mxArray *c5_lhs62 = NULL;
  const mxArray *c5_rhs63 = NULL;
  const mxArray *c5_lhs63 = NULL;
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1447321639U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c5_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("norm"), "name", "name", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c5_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c5_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c5_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_xnrm2"), "name", "name", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c5_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c5_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.blas.xnrm2"),
                  "name", "name", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c5_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c5_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p!below_threshold"),
                  "context", "context", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c5_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c5_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.refblas.xnrm2"),
                  "name", "name", 10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c5_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("realmin"), "name", "name", 11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c5_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_realmin"), "name", "name",
                  12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658444U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c5_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c5_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c5_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c5_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c5_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c5_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 18);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("intmax"), "name", "name", 18);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c5_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 19);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c5_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 20);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 20);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c5_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 21);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c5_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 22);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c5_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 23);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eig"), "name", "name", 23);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "resolved",
                  "resolved", 23);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c5_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 24);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 24);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c5_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 25);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c5_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_xgeev"), "name", "name",
                  26);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826004U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c5_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m"),
                  "context", "context", 27);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_lapack_xgeev"), "name",
                  "name", 27);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1301335668U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c5_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 28);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zggev"), "name",
                  "name", 28);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826018U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c5_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 29);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("realmin"), "name", "name", 29);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c5_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("sqrt"), "name", "name", 30);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 30);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c5_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_error"), "name", "name",
                  31);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 31);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c5_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 32);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825938U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c5_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 33);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eps"), "name", "name", 33);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c5_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 34);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 34);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c5_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_eps"), "name", "name", 35);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 35);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c5_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 36);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 36);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c5_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 37);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zlangeM"), "name",
                  "name", 37);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826020U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c5_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "context", "context", 38);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 38);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 38);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c5_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "context", "context", 39);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  39);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1350417854U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c5_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "context", "context", 40);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isnan"), "name", "name", 40);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 40);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c5_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 41);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c5_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "context", "context", 42);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 42);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c5_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "context", "context", 43);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 43);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c5_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 44);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isfinite"), "name", "name", 44);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c5_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 45);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 45);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c5_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 46);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isinf"), "name", "name", 46);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c5_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 47);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 47);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c5_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 48);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isnan"), "name", "name", 48);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 48);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c5_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 49);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 49);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c5_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 50);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zlascl"), "name",
                  "name", 50);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826022U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c5_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 51);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("realmin"), "name", "name", 51);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 51);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c5_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 52);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eps"), "name", "name", 52);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 52);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c5_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 53);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 53);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 53);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c5_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 54);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 54);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c5_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 55);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 55);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c5_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zggbal"), "name",
                  "name", 56);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826018U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c5_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "context", "context", 57);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 57);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c5_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_rows"),
                  "context", "context", 58);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 58);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c5_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_rows"),
                  "context", "context", 59);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 59);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c5_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_rows"),
                  "context", "context", 60);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 60);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 60);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c5_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 61);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 61);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 61);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c5_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_simtran"),
                  "context", "context", 62);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 62);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c5_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_simtran"),
                  "context", "context", 63);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 63);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c5_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c5_rhs0);
  sf_mex_destroy(&c5_lhs0);
  sf_mex_destroy(&c5_rhs1);
  sf_mex_destroy(&c5_lhs1);
  sf_mex_destroy(&c5_rhs2);
  sf_mex_destroy(&c5_lhs2);
  sf_mex_destroy(&c5_rhs3);
  sf_mex_destroy(&c5_lhs3);
  sf_mex_destroy(&c5_rhs4);
  sf_mex_destroy(&c5_lhs4);
  sf_mex_destroy(&c5_rhs5);
  sf_mex_destroy(&c5_lhs5);
  sf_mex_destroy(&c5_rhs6);
  sf_mex_destroy(&c5_lhs6);
  sf_mex_destroy(&c5_rhs7);
  sf_mex_destroy(&c5_lhs7);
  sf_mex_destroy(&c5_rhs8);
  sf_mex_destroy(&c5_lhs8);
  sf_mex_destroy(&c5_rhs9);
  sf_mex_destroy(&c5_lhs9);
  sf_mex_destroy(&c5_rhs10);
  sf_mex_destroy(&c5_lhs10);
  sf_mex_destroy(&c5_rhs11);
  sf_mex_destroy(&c5_lhs11);
  sf_mex_destroy(&c5_rhs12);
  sf_mex_destroy(&c5_lhs12);
  sf_mex_destroy(&c5_rhs13);
  sf_mex_destroy(&c5_lhs13);
  sf_mex_destroy(&c5_rhs14);
  sf_mex_destroy(&c5_lhs14);
  sf_mex_destroy(&c5_rhs15);
  sf_mex_destroy(&c5_lhs15);
  sf_mex_destroy(&c5_rhs16);
  sf_mex_destroy(&c5_lhs16);
  sf_mex_destroy(&c5_rhs17);
  sf_mex_destroy(&c5_lhs17);
  sf_mex_destroy(&c5_rhs18);
  sf_mex_destroy(&c5_lhs18);
  sf_mex_destroy(&c5_rhs19);
  sf_mex_destroy(&c5_lhs19);
  sf_mex_destroy(&c5_rhs20);
  sf_mex_destroy(&c5_lhs20);
  sf_mex_destroy(&c5_rhs21);
  sf_mex_destroy(&c5_lhs21);
  sf_mex_destroy(&c5_rhs22);
  sf_mex_destroy(&c5_lhs22);
  sf_mex_destroy(&c5_rhs23);
  sf_mex_destroy(&c5_lhs23);
  sf_mex_destroy(&c5_rhs24);
  sf_mex_destroy(&c5_lhs24);
  sf_mex_destroy(&c5_rhs25);
  sf_mex_destroy(&c5_lhs25);
  sf_mex_destroy(&c5_rhs26);
  sf_mex_destroy(&c5_lhs26);
  sf_mex_destroy(&c5_rhs27);
  sf_mex_destroy(&c5_lhs27);
  sf_mex_destroy(&c5_rhs28);
  sf_mex_destroy(&c5_lhs28);
  sf_mex_destroy(&c5_rhs29);
  sf_mex_destroy(&c5_lhs29);
  sf_mex_destroy(&c5_rhs30);
  sf_mex_destroy(&c5_lhs30);
  sf_mex_destroy(&c5_rhs31);
  sf_mex_destroy(&c5_lhs31);
  sf_mex_destroy(&c5_rhs32);
  sf_mex_destroy(&c5_lhs32);
  sf_mex_destroy(&c5_rhs33);
  sf_mex_destroy(&c5_lhs33);
  sf_mex_destroy(&c5_rhs34);
  sf_mex_destroy(&c5_lhs34);
  sf_mex_destroy(&c5_rhs35);
  sf_mex_destroy(&c5_lhs35);
  sf_mex_destroy(&c5_rhs36);
  sf_mex_destroy(&c5_lhs36);
  sf_mex_destroy(&c5_rhs37);
  sf_mex_destroy(&c5_lhs37);
  sf_mex_destroy(&c5_rhs38);
  sf_mex_destroy(&c5_lhs38);
  sf_mex_destroy(&c5_rhs39);
  sf_mex_destroy(&c5_lhs39);
  sf_mex_destroy(&c5_rhs40);
  sf_mex_destroy(&c5_lhs40);
  sf_mex_destroy(&c5_rhs41);
  sf_mex_destroy(&c5_lhs41);
  sf_mex_destroy(&c5_rhs42);
  sf_mex_destroy(&c5_lhs42);
  sf_mex_destroy(&c5_rhs43);
  sf_mex_destroy(&c5_lhs43);
  sf_mex_destroy(&c5_rhs44);
  sf_mex_destroy(&c5_lhs44);
  sf_mex_destroy(&c5_rhs45);
  sf_mex_destroy(&c5_lhs45);
  sf_mex_destroy(&c5_rhs46);
  sf_mex_destroy(&c5_lhs46);
  sf_mex_destroy(&c5_rhs47);
  sf_mex_destroy(&c5_lhs47);
  sf_mex_destroy(&c5_rhs48);
  sf_mex_destroy(&c5_lhs48);
  sf_mex_destroy(&c5_rhs49);
  sf_mex_destroy(&c5_lhs49);
  sf_mex_destroy(&c5_rhs50);
  sf_mex_destroy(&c5_lhs50);
  sf_mex_destroy(&c5_rhs51);
  sf_mex_destroy(&c5_lhs51);
  sf_mex_destroy(&c5_rhs52);
  sf_mex_destroy(&c5_lhs52);
  sf_mex_destroy(&c5_rhs53);
  sf_mex_destroy(&c5_lhs53);
  sf_mex_destroy(&c5_rhs54);
  sf_mex_destroy(&c5_lhs54);
  sf_mex_destroy(&c5_rhs55);
  sf_mex_destroy(&c5_lhs55);
  sf_mex_destroy(&c5_rhs56);
  sf_mex_destroy(&c5_lhs56);
  sf_mex_destroy(&c5_rhs57);
  sf_mex_destroy(&c5_lhs57);
  sf_mex_destroy(&c5_rhs58);
  sf_mex_destroy(&c5_lhs58);
  sf_mex_destroy(&c5_rhs59);
  sf_mex_destroy(&c5_lhs59);
  sf_mex_destroy(&c5_rhs60);
  sf_mex_destroy(&c5_lhs60);
  sf_mex_destroy(&c5_rhs61);
  sf_mex_destroy(&c5_lhs61);
  sf_mex_destroy(&c5_rhs62);
  sf_mex_destroy(&c5_lhs62);
  sf_mex_destroy(&c5_rhs63);
  sf_mex_destroy(&c5_lhs63);
}

static const mxArray *c5_emlrt_marshallOut(const char * c5_u)
{
  const mxArray *c5_y = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c5_u)), false);
  return c5_y;
}

static const mxArray *c5_b_emlrt_marshallOut(const uint32_T c5_u)
{
  const mxArray *c5_y = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 7, 0U, 0U, 0U, 0), false);
  return c5_y;
}

static void c5_b_info_helper(const mxArray **c5_info)
{
  const mxArray *c5_rhs64 = NULL;
  const mxArray *c5_lhs64 = NULL;
  const mxArray *c5_rhs65 = NULL;
  const mxArray *c5_lhs65 = NULL;
  const mxArray *c5_rhs66 = NULL;
  const mxArray *c5_lhs66 = NULL;
  const mxArray *c5_rhs67 = NULL;
  const mxArray *c5_lhs67 = NULL;
  const mxArray *c5_rhs68 = NULL;
  const mxArray *c5_lhs68 = NULL;
  const mxArray *c5_rhs69 = NULL;
  const mxArray *c5_lhs69 = NULL;
  const mxArray *c5_rhs70 = NULL;
  const mxArray *c5_lhs70 = NULL;
  const mxArray *c5_rhs71 = NULL;
  const mxArray *c5_lhs71 = NULL;
  const mxArray *c5_rhs72 = NULL;
  const mxArray *c5_lhs72 = NULL;
  const mxArray *c5_rhs73 = NULL;
  const mxArray *c5_lhs73 = NULL;
  const mxArray *c5_rhs74 = NULL;
  const mxArray *c5_lhs74 = NULL;
  const mxArray *c5_rhs75 = NULL;
  const mxArray *c5_lhs75 = NULL;
  const mxArray *c5_rhs76 = NULL;
  const mxArray *c5_lhs76 = NULL;
  const mxArray *c5_rhs77 = NULL;
  const mxArray *c5_lhs77 = NULL;
  const mxArray *c5_rhs78 = NULL;
  const mxArray *c5_lhs78 = NULL;
  const mxArray *c5_rhs79 = NULL;
  const mxArray *c5_lhs79 = NULL;
  const mxArray *c5_rhs80 = NULL;
  const mxArray *c5_lhs80 = NULL;
  const mxArray *c5_rhs81 = NULL;
  const mxArray *c5_lhs81 = NULL;
  const mxArray *c5_rhs82 = NULL;
  const mxArray *c5_lhs82 = NULL;
  const mxArray *c5_rhs83 = NULL;
  const mxArray *c5_lhs83 = NULL;
  const mxArray *c5_rhs84 = NULL;
  const mxArray *c5_lhs84 = NULL;
  const mxArray *c5_rhs85 = NULL;
  const mxArray *c5_lhs85 = NULL;
  const mxArray *c5_rhs86 = NULL;
  const mxArray *c5_lhs86 = NULL;
  const mxArray *c5_rhs87 = NULL;
  const mxArray *c5_lhs87 = NULL;
  const mxArray *c5_rhs88 = NULL;
  const mxArray *c5_lhs88 = NULL;
  const mxArray *c5_rhs89 = NULL;
  const mxArray *c5_lhs89 = NULL;
  const mxArray *c5_rhs90 = NULL;
  const mxArray *c5_lhs90 = NULL;
  const mxArray *c5_rhs91 = NULL;
  const mxArray *c5_lhs91 = NULL;
  const mxArray *c5_rhs92 = NULL;
  const mxArray *c5_lhs92 = NULL;
  const mxArray *c5_rhs93 = NULL;
  const mxArray *c5_lhs93 = NULL;
  const mxArray *c5_rhs94 = NULL;
  const mxArray *c5_lhs94 = NULL;
  const mxArray *c5_rhs95 = NULL;
  const mxArray *c5_lhs95 = NULL;
  const mxArray *c5_rhs96 = NULL;
  const mxArray *c5_lhs96 = NULL;
  const mxArray *c5_rhs97 = NULL;
  const mxArray *c5_lhs97 = NULL;
  const mxArray *c5_rhs98 = NULL;
  const mxArray *c5_lhs98 = NULL;
  const mxArray *c5_rhs99 = NULL;
  const mxArray *c5_lhs99 = NULL;
  const mxArray *c5_rhs100 = NULL;
  const mxArray *c5_lhs100 = NULL;
  const mxArray *c5_rhs101 = NULL;
  const mxArray *c5_lhs101 = NULL;
  const mxArray *c5_rhs102 = NULL;
  const mxArray *c5_lhs102 = NULL;
  const mxArray *c5_rhs103 = NULL;
  const mxArray *c5_lhs103 = NULL;
  const mxArray *c5_rhs104 = NULL;
  const mxArray *c5_lhs104 = NULL;
  const mxArray *c5_rhs105 = NULL;
  const mxArray *c5_lhs105 = NULL;
  const mxArray *c5_rhs106 = NULL;
  const mxArray *c5_lhs106 = NULL;
  const mxArray *c5_rhs107 = NULL;
  const mxArray *c5_lhs107 = NULL;
  const mxArray *c5_rhs108 = NULL;
  const mxArray *c5_lhs108 = NULL;
  const mxArray *c5_rhs109 = NULL;
  const mxArray *c5_lhs109 = NULL;
  const mxArray *c5_rhs110 = NULL;
  const mxArray *c5_lhs110 = NULL;
  const mxArray *c5_rhs111 = NULL;
  const mxArray *c5_lhs111 = NULL;
  const mxArray *c5_rhs112 = NULL;
  const mxArray *c5_lhs112 = NULL;
  const mxArray *c5_rhs113 = NULL;
  const mxArray *c5_lhs113 = NULL;
  const mxArray *c5_rhs114 = NULL;
  const mxArray *c5_lhs114 = NULL;
  const mxArray *c5_rhs115 = NULL;
  const mxArray *c5_lhs115 = NULL;
  const mxArray *c5_rhs116 = NULL;
  const mxArray *c5_lhs116 = NULL;
  const mxArray *c5_rhs117 = NULL;
  const mxArray *c5_lhs117 = NULL;
  const mxArray *c5_rhs118 = NULL;
  const mxArray *c5_lhs118 = NULL;
  const mxArray *c5_rhs119 = NULL;
  const mxArray *c5_lhs119 = NULL;
  const mxArray *c5_rhs120 = NULL;
  const mxArray *c5_lhs120 = NULL;
  const mxArray *c5_rhs121 = NULL;
  const mxArray *c5_lhs121 = NULL;
  const mxArray *c5_rhs122 = NULL;
  const mxArray *c5_lhs122 = NULL;
  const mxArray *c5_rhs123 = NULL;
  const mxArray *c5_lhs123 = NULL;
  const mxArray *c5_rhs124 = NULL;
  const mxArray *c5_lhs124 = NULL;
  const mxArray *c5_rhs125 = NULL;
  const mxArray *c5_lhs125 = NULL;
  const mxArray *c5_rhs126 = NULL;
  const mxArray *c5_lhs126 = NULL;
  const mxArray *c5_rhs127 = NULL;
  const mxArray *c5_lhs127 = NULL;
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "context", "context", 64);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 64);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 64);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c5_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_cols"),
                  "context", "context", 65);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 65);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c5_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_cols"),
                  "context", "context", 66);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 66);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 66);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c5_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "context", "context", 67);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 67);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 67);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 67);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c5_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 68);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 68);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 68);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c5_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 69);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zgghrd"), "name",
                  "name", 69);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "resolved", "resolved", 69);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826020U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c5_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 70);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 70);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c5_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 71);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 71);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c5_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 72);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 72);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c5_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 73);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 73);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 73);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c5_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 74);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 74);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 74);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c5_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 75);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zlartg"), "name",
                  "name", 75);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "resolved", "resolved", 75);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826022U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c5_rhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 76);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("realmin"), "name", "name", 76);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 76);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c5_rhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 77);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eps"), "name", "name", 77);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 77);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 77);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c5_rhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("fix"), "name", "name", 78);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 78);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m"), "resolved",
                  "resolved", 78);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c5_rhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m"), "context",
                  "context", 79);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 79);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c5_rhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m"), "context",
                  "context", 80);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_fix"), "name",
                  "name", 80);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_fix.m"),
                  "resolved", "resolved", 80);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658438U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c5_rhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 81);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mpower"), "name", "name", 81);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 81);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c5_rhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 82);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 82);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 82);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c5_rhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 83);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("ismatrix"), "name", "name", 83);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 83);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c5_rhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 84);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("power"), "name", "name", 84);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 84);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c5_rhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 85);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 85);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c5_rhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 86);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 86);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 86);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c5_rhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 87);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 87);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c5_rhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 88);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 88);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 88);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c5_rhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 89);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("floor"), "name", "name", 89);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 89);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c5_rhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 90);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 90);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 90);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c5_rhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 91);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 91);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c5_rhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 92);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 92);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 92);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c5_rhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 93);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 93);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 93);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c5_rhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m!absinf"),
                  "context", "context", 94);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 94);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 94);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 94);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c5_rhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 95);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 95);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 95);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 95);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c5_rhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 96);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 96);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 96);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 96);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c5_rhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 97);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  97);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 97);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 97);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1350417854U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c5_rhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 98);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("sqrt"), "name", "name", 98);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 98);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 98);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c5_rhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 99);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 99);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 99);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c5_rhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 100);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_zrot_rows"), "name",
                  "name", 100);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "resolved", "resolved", 100);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c5_rhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs100), "rhs", "rhs",
                  100);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs100), "lhs", "lhs",
                  100);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "context", "context", 101);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 101);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 101);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 101);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c5_rhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs101), "rhs", "rhs",
                  101);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs101), "lhs", "lhs",
                  101);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "context", "context", 102);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 102);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 102);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c5_rhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs102), "rhs", "rhs",
                  102);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs102), "lhs", "lhs",
                  102);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "context", "context", 103);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.conjtimes"),
                  "name", "name", 103);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 103);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/conjtimes.m"),
                  "resolved", "resolved", 103);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1360286186U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c5_rhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs103), "rhs", "rhs",
                  103);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs103), "lhs", "lhs",
                  103);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 104);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_zrot_cols"), "name",
                  "name", 104);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "resolved", "resolved", 104);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c5_rhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs104), "rhs", "rhs",
                  104);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs104), "lhs", "lhs",
                  104);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "context", "context", 105);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 105);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 105);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 105);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c5_rhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs105), "rhs", "rhs",
                  105);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs105), "lhs", "lhs",
                  105);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "context", "context", 106);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 106);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 106);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 106);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c5_rhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs106), "rhs", "rhs",
                  106);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs106), "lhs", "lhs",
                  106);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "context", "context", 107);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.conjtimes"),
                  "name", "name", 107);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 107);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/conjtimes.m"),
                  "resolved", "resolved", 107);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1360286186U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c5_rhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs107), "rhs", "rhs",
                  107);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs107), "lhs", "lhs",
                  107);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 108);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zhgeqz"), "name",
                  "name", 108);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 108);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1368190232U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c5_rhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs108), "rhs", "rhs",
                  108);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs108), "lhs", "lhs",
                  108);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 109);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 109);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 109);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 109);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c5_rhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs109), "rhs", "rhs",
                  109);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs109), "lhs", "lhs",
                  109);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 110);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eps"), "name", "name", 110);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 110);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c5_rhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs110), "rhs", "rhs",
                  110);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs110), "lhs", "lhs",
                  110);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 111);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("realmin"), "name", "name", 111);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 111);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 111);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c5_rhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs111), "rhs", "rhs",
                  111);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs111), "lhs", "lhs",
                  111);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 112);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zlanhs"), "name",
                  "name", 112);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "resolved", "resolved", 112);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826020U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c5_rhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs112), "rhs", "rhs",
                  112);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs112), "lhs", "lhs",
                  112);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 113);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 113);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c5_rhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs113), "rhs", "rhs",
                  113);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs113), "lhs", "lhs",
                  113);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 114);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 114);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 114);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c5_rhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs114), "rhs", "rhs",
                  114);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs114), "lhs", "lhs",
                  114);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 115);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 115);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 115);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 115);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c5_rhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs115), "rhs", "rhs",
                  115);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs115), "lhs", "lhs",
                  115);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 116);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 116);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 116);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 116);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c5_rhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs116), "rhs", "rhs",
                  116);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs116), "lhs", "lhs",
                  116);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 117);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("sqrt"), "name", "name", 117);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 117);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 117);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c5_rhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs117), "rhs", "rhs",
                  117);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs117), "lhs", "lhs",
                  117);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 118);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 118);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 118);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 118);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c5_rhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs118), "rhs", "rhs",
                  118);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs118), "lhs", "lhs",
                  118);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 119);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 119);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 119);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 119);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c5_rhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs119), "rhs", "rhs",
                  119);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs119), "lhs", "lhs",
                  119);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 120);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 120);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 120);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 120);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c5_rhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs120), "rhs", "rhs",
                  120);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs120), "lhs", "lhs",
                  120);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 121);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 121);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 121);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 121);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c5_rhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs121), "rhs", "rhs",
                  121);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs121), "lhs", "lhs",
                  121);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 122);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 122);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 122);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c5_rhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs122), "rhs", "rhs",
                  122);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs122), "lhs", "lhs",
                  122);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m!abs1"),
                  "context", "context", 123);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 123);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 123);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 123);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c5_rhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs123), "rhs", "rhs",
                  123);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs123), "lhs", "lhs",
                  123);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 124);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 124);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 124);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 124);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c5_rhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs124), "rhs", "rhs",
                  124);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs124), "lhs", "lhs",
                  124);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 125);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zlartg"), "name",
                  "name", 125);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 125);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "resolved", "resolved", 125);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826022U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c5_rhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs125), "rhs", "rhs",
                  125);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs125), "lhs", "lhs",
                  125);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 126);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_zrot_cols"), "name",
                  "name", 126);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "resolved", "resolved", 126);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c5_rhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs126), "rhs", "rhs",
                  126);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs126), "lhs", "lhs",
                  126);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 127);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mod"), "name", "name", 127);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "resolved",
                  "resolved", 127);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c5_rhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs127), "rhs", "rhs",
                  127);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs127), "lhs", "lhs",
                  127);
  sf_mex_destroy(&c5_rhs64);
  sf_mex_destroy(&c5_lhs64);
  sf_mex_destroy(&c5_rhs65);
  sf_mex_destroy(&c5_lhs65);
  sf_mex_destroy(&c5_rhs66);
  sf_mex_destroy(&c5_lhs66);
  sf_mex_destroy(&c5_rhs67);
  sf_mex_destroy(&c5_lhs67);
  sf_mex_destroy(&c5_rhs68);
  sf_mex_destroy(&c5_lhs68);
  sf_mex_destroy(&c5_rhs69);
  sf_mex_destroy(&c5_lhs69);
  sf_mex_destroy(&c5_rhs70);
  sf_mex_destroy(&c5_lhs70);
  sf_mex_destroy(&c5_rhs71);
  sf_mex_destroy(&c5_lhs71);
  sf_mex_destroy(&c5_rhs72);
  sf_mex_destroy(&c5_lhs72);
  sf_mex_destroy(&c5_rhs73);
  sf_mex_destroy(&c5_lhs73);
  sf_mex_destroy(&c5_rhs74);
  sf_mex_destroy(&c5_lhs74);
  sf_mex_destroy(&c5_rhs75);
  sf_mex_destroy(&c5_lhs75);
  sf_mex_destroy(&c5_rhs76);
  sf_mex_destroy(&c5_lhs76);
  sf_mex_destroy(&c5_rhs77);
  sf_mex_destroy(&c5_lhs77);
  sf_mex_destroy(&c5_rhs78);
  sf_mex_destroy(&c5_lhs78);
  sf_mex_destroy(&c5_rhs79);
  sf_mex_destroy(&c5_lhs79);
  sf_mex_destroy(&c5_rhs80);
  sf_mex_destroy(&c5_lhs80);
  sf_mex_destroy(&c5_rhs81);
  sf_mex_destroy(&c5_lhs81);
  sf_mex_destroy(&c5_rhs82);
  sf_mex_destroy(&c5_lhs82);
  sf_mex_destroy(&c5_rhs83);
  sf_mex_destroy(&c5_lhs83);
  sf_mex_destroy(&c5_rhs84);
  sf_mex_destroy(&c5_lhs84);
  sf_mex_destroy(&c5_rhs85);
  sf_mex_destroy(&c5_lhs85);
  sf_mex_destroy(&c5_rhs86);
  sf_mex_destroy(&c5_lhs86);
  sf_mex_destroy(&c5_rhs87);
  sf_mex_destroy(&c5_lhs87);
  sf_mex_destroy(&c5_rhs88);
  sf_mex_destroy(&c5_lhs88);
  sf_mex_destroy(&c5_rhs89);
  sf_mex_destroy(&c5_lhs89);
  sf_mex_destroy(&c5_rhs90);
  sf_mex_destroy(&c5_lhs90);
  sf_mex_destroy(&c5_rhs91);
  sf_mex_destroy(&c5_lhs91);
  sf_mex_destroy(&c5_rhs92);
  sf_mex_destroy(&c5_lhs92);
  sf_mex_destroy(&c5_rhs93);
  sf_mex_destroy(&c5_lhs93);
  sf_mex_destroy(&c5_rhs94);
  sf_mex_destroy(&c5_lhs94);
  sf_mex_destroy(&c5_rhs95);
  sf_mex_destroy(&c5_lhs95);
  sf_mex_destroy(&c5_rhs96);
  sf_mex_destroy(&c5_lhs96);
  sf_mex_destroy(&c5_rhs97);
  sf_mex_destroy(&c5_lhs97);
  sf_mex_destroy(&c5_rhs98);
  sf_mex_destroy(&c5_lhs98);
  sf_mex_destroy(&c5_rhs99);
  sf_mex_destroy(&c5_lhs99);
  sf_mex_destroy(&c5_rhs100);
  sf_mex_destroy(&c5_lhs100);
  sf_mex_destroy(&c5_rhs101);
  sf_mex_destroy(&c5_lhs101);
  sf_mex_destroy(&c5_rhs102);
  sf_mex_destroy(&c5_lhs102);
  sf_mex_destroy(&c5_rhs103);
  sf_mex_destroy(&c5_lhs103);
  sf_mex_destroy(&c5_rhs104);
  sf_mex_destroy(&c5_lhs104);
  sf_mex_destroy(&c5_rhs105);
  sf_mex_destroy(&c5_lhs105);
  sf_mex_destroy(&c5_rhs106);
  sf_mex_destroy(&c5_lhs106);
  sf_mex_destroy(&c5_rhs107);
  sf_mex_destroy(&c5_lhs107);
  sf_mex_destroy(&c5_rhs108);
  sf_mex_destroy(&c5_lhs108);
  sf_mex_destroy(&c5_rhs109);
  sf_mex_destroy(&c5_lhs109);
  sf_mex_destroy(&c5_rhs110);
  sf_mex_destroy(&c5_lhs110);
  sf_mex_destroy(&c5_rhs111);
  sf_mex_destroy(&c5_lhs111);
  sf_mex_destroy(&c5_rhs112);
  sf_mex_destroy(&c5_lhs112);
  sf_mex_destroy(&c5_rhs113);
  sf_mex_destroy(&c5_lhs113);
  sf_mex_destroy(&c5_rhs114);
  sf_mex_destroy(&c5_lhs114);
  sf_mex_destroy(&c5_rhs115);
  sf_mex_destroy(&c5_lhs115);
  sf_mex_destroy(&c5_rhs116);
  sf_mex_destroy(&c5_lhs116);
  sf_mex_destroy(&c5_rhs117);
  sf_mex_destroy(&c5_lhs117);
  sf_mex_destroy(&c5_rhs118);
  sf_mex_destroy(&c5_lhs118);
  sf_mex_destroy(&c5_rhs119);
  sf_mex_destroy(&c5_lhs119);
  sf_mex_destroy(&c5_rhs120);
  sf_mex_destroy(&c5_lhs120);
  sf_mex_destroy(&c5_rhs121);
  sf_mex_destroy(&c5_lhs121);
  sf_mex_destroy(&c5_rhs122);
  sf_mex_destroy(&c5_lhs122);
  sf_mex_destroy(&c5_rhs123);
  sf_mex_destroy(&c5_lhs123);
  sf_mex_destroy(&c5_rhs124);
  sf_mex_destroy(&c5_lhs124);
  sf_mex_destroy(&c5_rhs125);
  sf_mex_destroy(&c5_lhs125);
  sf_mex_destroy(&c5_rhs126);
  sf_mex_destroy(&c5_lhs126);
  sf_mex_destroy(&c5_rhs127);
  sf_mex_destroy(&c5_lhs127);
}

static void c5_c_info_helper(const mxArray **c5_info)
{
  const mxArray *c5_rhs128 = NULL;
  const mxArray *c5_lhs128 = NULL;
  const mxArray *c5_rhs129 = NULL;
  const mxArray *c5_lhs129 = NULL;
  const mxArray *c5_rhs130 = NULL;
  const mxArray *c5_lhs130 = NULL;
  const mxArray *c5_rhs131 = NULL;
  const mxArray *c5_lhs131 = NULL;
  const mxArray *c5_rhs132 = NULL;
  const mxArray *c5_lhs132 = NULL;
  const mxArray *c5_rhs133 = NULL;
  const mxArray *c5_lhs133 = NULL;
  const mxArray *c5_rhs134 = NULL;
  const mxArray *c5_lhs134 = NULL;
  const mxArray *c5_rhs135 = NULL;
  const mxArray *c5_lhs135 = NULL;
  const mxArray *c5_rhs136 = NULL;
  const mxArray *c5_lhs136 = NULL;
  const mxArray *c5_rhs137 = NULL;
  const mxArray *c5_lhs137 = NULL;
  const mxArray *c5_rhs138 = NULL;
  const mxArray *c5_lhs138 = NULL;
  const mxArray *c5_rhs139 = NULL;
  const mxArray *c5_lhs139 = NULL;
  const mxArray *c5_rhs140 = NULL;
  const mxArray *c5_lhs140 = NULL;
  const mxArray *c5_rhs141 = NULL;
  const mxArray *c5_lhs141 = NULL;
  const mxArray *c5_rhs142 = NULL;
  const mxArray *c5_lhs142 = NULL;
  const mxArray *c5_rhs143 = NULL;
  const mxArray *c5_lhs143 = NULL;
  const mxArray *c5_rhs144 = NULL;
  const mxArray *c5_lhs144 = NULL;
  const mxArray *c5_rhs145 = NULL;
  const mxArray *c5_lhs145 = NULL;
  const mxArray *c5_rhs146 = NULL;
  const mxArray *c5_lhs146 = NULL;
  const mxArray *c5_rhs147 = NULL;
  const mxArray *c5_lhs147 = NULL;
  const mxArray *c5_rhs148 = NULL;
  const mxArray *c5_lhs148 = NULL;
  const mxArray *c5_rhs149 = NULL;
  const mxArray *c5_lhs149 = NULL;
  const mxArray *c5_rhs150 = NULL;
  const mxArray *c5_lhs150 = NULL;
  const mxArray *c5_rhs151 = NULL;
  const mxArray *c5_lhs151 = NULL;
  const mxArray *c5_rhs152 = NULL;
  const mxArray *c5_lhs152 = NULL;
  const mxArray *c5_rhs153 = NULL;
  const mxArray *c5_lhs153 = NULL;
  const mxArray *c5_rhs154 = NULL;
  const mxArray *c5_lhs154 = NULL;
  const mxArray *c5_rhs155 = NULL;
  const mxArray *c5_lhs155 = NULL;
  const mxArray *c5_rhs156 = NULL;
  const mxArray *c5_lhs156 = NULL;
  const mxArray *c5_rhs157 = NULL;
  const mxArray *c5_lhs157 = NULL;
  const mxArray *c5_rhs158 = NULL;
  const mxArray *c5_lhs158 = NULL;
  const mxArray *c5_rhs159 = NULL;
  const mxArray *c5_lhs159 = NULL;
  const mxArray *c5_rhs160 = NULL;
  const mxArray *c5_lhs160 = NULL;
  const mxArray *c5_rhs161 = NULL;
  const mxArray *c5_lhs161 = NULL;
  const mxArray *c5_rhs162 = NULL;
  const mxArray *c5_lhs162 = NULL;
  const mxArray *c5_rhs163 = NULL;
  const mxArray *c5_lhs163 = NULL;
  const mxArray *c5_rhs164 = NULL;
  const mxArray *c5_lhs164 = NULL;
  const mxArray *c5_rhs165 = NULL;
  const mxArray *c5_lhs165 = NULL;
  const mxArray *c5_rhs166 = NULL;
  const mxArray *c5_lhs166 = NULL;
  const mxArray *c5_rhs167 = NULL;
  const mxArray *c5_lhs167 = NULL;
  const mxArray *c5_rhs168 = NULL;
  const mxArray *c5_lhs168 = NULL;
  const mxArray *c5_rhs169 = NULL;
  const mxArray *c5_lhs169 = NULL;
  const mxArray *c5_rhs170 = NULL;
  const mxArray *c5_lhs170 = NULL;
  const mxArray *c5_rhs171 = NULL;
  const mxArray *c5_lhs171 = NULL;
  const mxArray *c5_rhs172 = NULL;
  const mxArray *c5_lhs172 = NULL;
  const mxArray *c5_rhs173 = NULL;
  const mxArray *c5_lhs173 = NULL;
  const mxArray *c5_rhs174 = NULL;
  const mxArray *c5_lhs174 = NULL;
  const mxArray *c5_rhs175 = NULL;
  const mxArray *c5_lhs175 = NULL;
  const mxArray *c5_rhs176 = NULL;
  const mxArray *c5_lhs176 = NULL;
  const mxArray *c5_rhs177 = NULL;
  const mxArray *c5_lhs177 = NULL;
  const mxArray *c5_rhs178 = NULL;
  const mxArray *c5_lhs178 = NULL;
  const mxArray *c5_rhs179 = NULL;
  const mxArray *c5_lhs179 = NULL;
  const mxArray *c5_rhs180 = NULL;
  const mxArray *c5_lhs180 = NULL;
  const mxArray *c5_rhs181 = NULL;
  const mxArray *c5_lhs181 = NULL;
  const mxArray *c5_rhs182 = NULL;
  const mxArray *c5_lhs182 = NULL;
  const mxArray *c5_rhs183 = NULL;
  const mxArray *c5_lhs183 = NULL;
  const mxArray *c5_rhs184 = NULL;
  const mxArray *c5_lhs184 = NULL;
  const mxArray *c5_rhs185 = NULL;
  const mxArray *c5_lhs185 = NULL;
  const mxArray *c5_rhs186 = NULL;
  const mxArray *c5_lhs186 = NULL;
  const mxArray *c5_rhs187 = NULL;
  const mxArray *c5_lhs187 = NULL;
  const mxArray *c5_rhs188 = NULL;
  const mxArray *c5_lhs188 = NULL;
  const mxArray *c5_rhs189 = NULL;
  const mxArray *c5_lhs189 = NULL;
  const mxArray *c5_rhs190 = NULL;
  const mxArray *c5_lhs190 = NULL;
  const mxArray *c5_rhs191 = NULL;
  const mxArray *c5_lhs191 = NULL;
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 128);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 128);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 128);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 128);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c5_rhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs128), "rhs", "rhs",
                  128);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs128), "lhs", "lhs",
                  128);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 129);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 129);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 129);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 129);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c5_rhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs129), "rhs", "rhs",
                  129);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs129), "lhs", "lhs",
                  129);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 130);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 130);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 130);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 130);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c5_rhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs130), "rhs", "rhs",
                  130);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs130), "lhs", "lhs",
                  130);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 131);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 131);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 131);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c5_rhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs131), "rhs", "rhs",
                  131);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs131), "lhs", "lhs",
                  131);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 132);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 132);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 132);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 132);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c5_rhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs132), "rhs", "rhs",
                  132);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs132), "lhs", "lhs",
                  132);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!intmod"), "context",
                  "context", 133);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 133);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 133);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 133);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c5_rhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs133), "rhs", "rhs",
                  133);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs133), "lhs", "lhs",
                  133);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 134);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 134);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 134);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 134);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c5_rhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs134), "rhs", "rhs",
                  134);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs134), "lhs", "lhs",
                  134);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 135);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_div"), "name", "name", 135);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 135);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c5_rhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs135), "rhs", "rhs",
                  135);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs135), "lhs", "lhs",
                  135);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 136);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 136);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 136);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c5_rhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs136), "rhs", "rhs",
                  136);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs136), "lhs", "lhs",
                  136);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p!eml_fldiv"),
                  "context", "context", 137);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 137);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 137);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c5_rhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs137), "rhs", "rhs",
                  137);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs137), "lhs", "lhs",
                  137);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p!eml_fldiv"),
                  "context", "context", 138);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 138);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 138);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c5_rhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs138), "rhs", "rhs",
                  138);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs138), "lhs", "lhs",
                  138);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p!eml_fldiv"),
                  "context", "context", 139);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 139);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 139);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c5_rhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs139), "rhs", "rhs",
                  139);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs139), "lhs", "lhs",
                  139);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 140);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("sqrt"), "name", "name", 140);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 140);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c5_rhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs140), "rhs", "rhs",
                  140);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs140), "lhs", "lhs",
                  140);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 141);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("rdivide"), "name", "name", 141);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 141);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 141);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c5_rhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs141), "rhs", "rhs",
                  141);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs141), "lhs", "lhs",
                  141);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 142);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 142);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 142);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c5_rhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs142), "rhs", "rhs",
                  142);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs142), "lhs", "lhs",
                  142);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 143);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 143);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 143);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 143);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c5_rhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs143), "rhs", "rhs",
                  143);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs143), "lhs", "lhs",
                  143);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 144);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_div"), "name", "name", 144);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 144);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c5_rhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs144), "rhs", "rhs",
                  144);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs144), "lhs", "lhs",
                  144);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 145);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isnan"), "name", "name", 145);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 145);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c5_rhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs145), "rhs", "rhs",
                  145);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs145), "lhs", "lhs",
                  145);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 146);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 146);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 146);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c5_rhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs146), "rhs", "rhs",
                  146);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs146), "lhs", "lhs",
                  146);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 147);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isinf"), "name", "name", 147);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 147);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c5_rhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs147), "rhs", "rhs",
                  147);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs147), "lhs", "lhs",
                  147);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 148);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_guarded_inf"), "name",
                  "name", 148);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_inf.m"),
                  "resolved", "resolved", 148);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c5_rhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs148), "rhs", "rhs",
                  148);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs148), "lhs", "lhs",
                  148);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_inf.m"),
                  "context", "context", 149);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 149);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 149);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c5_rhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs149), "rhs", "rhs",
                  149);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs149), "lhs", "lhs",
                  149);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 150);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("realmax"), "name", "name", 150);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m"), "resolved",
                  "resolved", 150);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c5_rhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs150), "rhs", "rhs",
                  150);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs150), "lhs", "lhs",
                  150);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m"), "context",
                  "context", 151);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_realmax"), "name", "name",
                  151);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmax.m"), "resolved",
                  "resolved", 151);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c5_rhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs151), "rhs", "rhs",
                  151);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs151), "lhs", "lhs",
                  151);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmax.m"), "context",
                  "context", 152);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 152);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 152);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c5_rhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs152), "rhs", "rhs",
                  152);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs152), "lhs", "lhs",
                  152);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 153);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mrdivide"), "name", "name",
                  153);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 153);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c5_rhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs153), "rhs", "rhs",
                  153);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs153), "lhs", "lhs",
                  153);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 154);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 154);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 154);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 154);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c5_rhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs154), "rhs", "rhs",
                  154);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs154), "lhs", "lhs",
                  154);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 155);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("rdivide"), "name", "name", 155);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 155);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c5_rhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs155), "rhs", "rhs",
                  155);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs155), "lhs", "lhs",
                  155);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 156);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_hypot"), "name",
                  "name", 156);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "resolved", "resolved", 156);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c5_rhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs156), "rhs", "rhs",
                  156);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs156), "lhs", "lhs",
                  156);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "context", "context", 157);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 157);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 157);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c5_rhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs157), "rhs", "rhs",
                  157);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs157), "lhs", "lhs",
                  157);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "context", "context", 158);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  158);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 158);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 158);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1350417854U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c5_rhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs158), "rhs", "rhs",
                  158);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs158), "lhs", "lhs",
                  158);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 159);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 159);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 159);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 159);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c5_rhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs159), "rhs", "rhs",
                  159);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs159), "lhs", "lhs",
                  159);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 160);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_zrot_rows"), "name",
                  "name", 160);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 160);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "resolved", "resolved", 160);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c5_rhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs160), "rhs", "rhs",
                  160);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs160), "lhs", "lhs",
                  160);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 161);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_div"), "name", "name", 161);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 161);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 161);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c5_rhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs161), "rhs", "rhs",
                  161);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs161), "lhs", "lhs",
                  161);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p!equalsize"),
                  "context", "context", 162);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("max"), "name", "name", 162);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 162);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "resolved",
                  "resolved", 162);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1311262516U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c5_rhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs162), "rhs", "rhs",
                  162);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs162), "lhs", "lhs",
                  162);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "context",
                  "context", 163);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 163);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 163);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1378303184U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c5_rhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs163), "rhs", "rhs",
                  163);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs163), "lhs", "lhs",
                  163);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 164);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 164);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 164);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 164);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 164);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 164);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 164);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 164);
  sf_mex_assign(&c5_rhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs164), "rhs", "rhs",
                  164);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs164), "lhs", "lhs",
                  164);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 165);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 165);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 165);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 165);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 165);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 165);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 165);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 165);
  sf_mex_assign(&c5_rhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs165), "rhs", "rhs",
                  165);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs165), "lhs", "lhs",
                  165);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 166);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 166);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 166);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 166);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 166);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 166);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 166);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 166);
  sf_mex_assign(&c5_rhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs166), "rhs", "rhs",
                  166);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs166), "lhs", "lhs",
                  166);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 167);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 167);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 167);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 167);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 167);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 167);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 167);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 167);
  sf_mex_assign(&c5_rhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs167), "rhs", "rhs",
                  167);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs167), "lhs", "lhs",
                  167);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 168);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 168);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 168);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 168);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 168);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 168);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 168);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 168);
  sf_mex_assign(&c5_rhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs168), "rhs", "rhs",
                  168);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs168), "lhs", "lhs",
                  168);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 169);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_warning"), "name", "name",
                  169);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 169);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 169);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 169);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 169);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 169);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 169);
  sf_mex_assign(&c5_rhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs169), "rhs", "rhs",
                  169);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs169), "lhs", "lhs",
                  169);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 170);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mpower"), "name", "name", 170);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 170);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 170);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 170);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 170);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 170);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 170);
  sf_mex_assign(&c5_rhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs170), "rhs", "rhs",
                  170);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs170), "lhs", "lhs",
                  170);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 171);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("inv"), "name", "name", 171);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 171);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 171);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 171);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 171);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 171);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 171);
  sf_mex_assign(&c5_rhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs171), "rhs", "rhs",
                  171);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs171), "lhs", "lhs",
                  171);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 172);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 172);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 172);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 172);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 172);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 172);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 172);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 172);
  sf_mex_assign(&c5_rhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs172), "rhs", "rhs",
                  172);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs172), "lhs", "lhs",
                  172);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 173);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 173);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 173);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 173);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 173);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 173);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 173);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 173);
  sf_mex_assign(&c5_rhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs173), "rhs", "rhs",
                  173);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs173), "lhs", "lhs",
                  173);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 174);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_div"), "name", "name", 174);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 174);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 174);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 174);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 174);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 174);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 174);
  sf_mex_assign(&c5_rhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs174), "rhs", "rhs",
                  174);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs174), "lhs", "lhs",
                  174);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 175);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 175);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 175);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 175);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 175);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 175);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 175);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 175);
  sf_mex_assign(&c5_rhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs175), "rhs", "rhs",
                  175);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs175), "lhs", "lhs",
                  175);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 176);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("norm"), "name", "name", 176);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 176);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 176);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 176);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 176);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 176);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 176);
  sf_mex_assign(&c5_rhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs176), "rhs", "rhs",
                  176);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs176), "lhs", "lhs",
                  176);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 177);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 177);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 177);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 177);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 177);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 177);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 177);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 177);
  sf_mex_assign(&c5_rhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs177), "rhs", "rhs",
                  177);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs177), "lhs", "lhs",
                  177);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 178);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 178);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 178);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 178);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 178);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 178);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 178);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 178);
  sf_mex_assign(&c5_rhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs178), "rhs", "rhs",
                  178);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs178), "lhs", "lhs",
                  178);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 179);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isnan"), "name", "name", 179);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 179);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 179);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 179);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 179);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 179);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 179);
  sf_mex_assign(&c5_rhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs179), "rhs", "rhs",
                  179);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs179), "lhs", "lhs",
                  179);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 180);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 180);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 180);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 180);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 180);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 180);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 180);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 180);
  sf_mex_assign(&c5_rhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs180), "rhs", "rhs",
                  180);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs180), "lhs", "lhs",
                  180);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 181);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_warning"), "name", "name",
                  181);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 181);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 181);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 181);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 181);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 181);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 181);
  sf_mex_assign(&c5_rhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs181), "rhs", "rhs",
                  181);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs181), "lhs", "lhs",
                  181);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 182);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isnan"), "name", "name", 182);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 182);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 182);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 182);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 182);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 182);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 182);
  sf_mex_assign(&c5_rhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs182), "rhs", "rhs",
                  182);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs182), "lhs", "lhs",
                  182);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 183);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eps"), "name", "name", 183);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 183);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 183);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 183);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 183);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 183);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 183);
  sf_mex_assign(&c5_rhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs183), "rhs", "rhs",
                  183);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs183), "lhs", "lhs",
                  183);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 184);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_flt2str"), "name", "name",
                  184);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 184);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "resolved",
                  "resolved", 184);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 184);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 184);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 184);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 184);
  sf_mex_assign(&c5_rhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs184), "rhs", "rhs",
                  184);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs184), "lhs", "lhs",
                  184);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "context",
                  "context", 185);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "name", "name", 185);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 185);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m"), "resolved",
                  "resolved", 185);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1319737168U), "fileTimeLo",
                  "fileTimeLo", 185);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 185);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 185);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 185);
  sf_mex_assign(&c5_rhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs185), "rhs", "rhs",
                  185);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs185), "lhs", "lhs",
                  185);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 186);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mrdivide"), "name", "name",
                  186);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 186);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 186);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 186);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 186);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 186);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 186);
  sf_mex_assign(&c5_rhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs186), "rhs", "rhs",
                  186);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs186), "lhs", "lhs",
                  186);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 187);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("sin"), "name", "name", 187);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 187);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 187);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 187);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 187);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 187);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 187);
  sf_mex_assign(&c5_rhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs187), "rhs", "rhs",
                  187);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs187), "lhs", "lhs",
                  187);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 188);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 188);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 188);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 188);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825936U), "fileTimeLo",
                  "fileTimeLo", 188);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 188);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 188);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 188);
  sf_mex_assign(&c5_rhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs188), "rhs", "rhs",
                  188);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs188), "lhs", "lhs",
                  188);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 189);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("cos"), "name", "name", 189);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 189);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 189);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837572U), "fileTimeLo",
                  "fileTimeLo", 189);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 189);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 189);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 189);
  sf_mex_assign(&c5_rhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs189), "rhs", "rhs",
                  189);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs189), "lhs", "lhs",
                  189);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 190);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 190);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 190);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 190);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825922U), "fileTimeLo",
                  "fileTimeLo", 190);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 190);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 190);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 190);
  sf_mex_assign(&c5_rhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs190), "rhs", "rhs",
                  190);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs190), "lhs", "lhs",
                  190);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 191);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("exp"), "name", "name", 191);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 191);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m"), "resolved",
                  "resolved", 191);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837580U), "fileTimeLo",
                  "fileTimeLo", 191);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 191);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 191);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 191);
  sf_mex_assign(&c5_rhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs191), "rhs", "rhs",
                  191);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs191), "lhs", "lhs",
                  191);
  sf_mex_destroy(&c5_rhs128);
  sf_mex_destroy(&c5_lhs128);
  sf_mex_destroy(&c5_rhs129);
  sf_mex_destroy(&c5_lhs129);
  sf_mex_destroy(&c5_rhs130);
  sf_mex_destroy(&c5_lhs130);
  sf_mex_destroy(&c5_rhs131);
  sf_mex_destroy(&c5_lhs131);
  sf_mex_destroy(&c5_rhs132);
  sf_mex_destroy(&c5_lhs132);
  sf_mex_destroy(&c5_rhs133);
  sf_mex_destroy(&c5_lhs133);
  sf_mex_destroy(&c5_rhs134);
  sf_mex_destroy(&c5_lhs134);
  sf_mex_destroy(&c5_rhs135);
  sf_mex_destroy(&c5_lhs135);
  sf_mex_destroy(&c5_rhs136);
  sf_mex_destroy(&c5_lhs136);
  sf_mex_destroy(&c5_rhs137);
  sf_mex_destroy(&c5_lhs137);
  sf_mex_destroy(&c5_rhs138);
  sf_mex_destroy(&c5_lhs138);
  sf_mex_destroy(&c5_rhs139);
  sf_mex_destroy(&c5_lhs139);
  sf_mex_destroy(&c5_rhs140);
  sf_mex_destroy(&c5_lhs140);
  sf_mex_destroy(&c5_rhs141);
  sf_mex_destroy(&c5_lhs141);
  sf_mex_destroy(&c5_rhs142);
  sf_mex_destroy(&c5_lhs142);
  sf_mex_destroy(&c5_rhs143);
  sf_mex_destroy(&c5_lhs143);
  sf_mex_destroy(&c5_rhs144);
  sf_mex_destroy(&c5_lhs144);
  sf_mex_destroy(&c5_rhs145);
  sf_mex_destroy(&c5_lhs145);
  sf_mex_destroy(&c5_rhs146);
  sf_mex_destroy(&c5_lhs146);
  sf_mex_destroy(&c5_rhs147);
  sf_mex_destroy(&c5_lhs147);
  sf_mex_destroy(&c5_rhs148);
  sf_mex_destroy(&c5_lhs148);
  sf_mex_destroy(&c5_rhs149);
  sf_mex_destroy(&c5_lhs149);
  sf_mex_destroy(&c5_rhs150);
  sf_mex_destroy(&c5_lhs150);
  sf_mex_destroy(&c5_rhs151);
  sf_mex_destroy(&c5_lhs151);
  sf_mex_destroy(&c5_rhs152);
  sf_mex_destroy(&c5_lhs152);
  sf_mex_destroy(&c5_rhs153);
  sf_mex_destroy(&c5_lhs153);
  sf_mex_destroy(&c5_rhs154);
  sf_mex_destroy(&c5_lhs154);
  sf_mex_destroy(&c5_rhs155);
  sf_mex_destroy(&c5_lhs155);
  sf_mex_destroy(&c5_rhs156);
  sf_mex_destroy(&c5_lhs156);
  sf_mex_destroy(&c5_rhs157);
  sf_mex_destroy(&c5_lhs157);
  sf_mex_destroy(&c5_rhs158);
  sf_mex_destroy(&c5_lhs158);
  sf_mex_destroy(&c5_rhs159);
  sf_mex_destroy(&c5_lhs159);
  sf_mex_destroy(&c5_rhs160);
  sf_mex_destroy(&c5_lhs160);
  sf_mex_destroy(&c5_rhs161);
  sf_mex_destroy(&c5_lhs161);
  sf_mex_destroy(&c5_rhs162);
  sf_mex_destroy(&c5_lhs162);
  sf_mex_destroy(&c5_rhs163);
  sf_mex_destroy(&c5_lhs163);
  sf_mex_destroy(&c5_rhs164);
  sf_mex_destroy(&c5_lhs164);
  sf_mex_destroy(&c5_rhs165);
  sf_mex_destroy(&c5_lhs165);
  sf_mex_destroy(&c5_rhs166);
  sf_mex_destroy(&c5_lhs166);
  sf_mex_destroy(&c5_rhs167);
  sf_mex_destroy(&c5_lhs167);
  sf_mex_destroy(&c5_rhs168);
  sf_mex_destroy(&c5_lhs168);
  sf_mex_destroy(&c5_rhs169);
  sf_mex_destroy(&c5_lhs169);
  sf_mex_destroy(&c5_rhs170);
  sf_mex_destroy(&c5_lhs170);
  sf_mex_destroy(&c5_rhs171);
  sf_mex_destroy(&c5_lhs171);
  sf_mex_destroy(&c5_rhs172);
  sf_mex_destroy(&c5_lhs172);
  sf_mex_destroy(&c5_rhs173);
  sf_mex_destroy(&c5_lhs173);
  sf_mex_destroy(&c5_rhs174);
  sf_mex_destroy(&c5_lhs174);
  sf_mex_destroy(&c5_rhs175);
  sf_mex_destroy(&c5_lhs175);
  sf_mex_destroy(&c5_rhs176);
  sf_mex_destroy(&c5_lhs176);
  sf_mex_destroy(&c5_rhs177);
  sf_mex_destroy(&c5_lhs177);
  sf_mex_destroy(&c5_rhs178);
  sf_mex_destroy(&c5_lhs178);
  sf_mex_destroy(&c5_rhs179);
  sf_mex_destroy(&c5_lhs179);
  sf_mex_destroy(&c5_rhs180);
  sf_mex_destroy(&c5_lhs180);
  sf_mex_destroy(&c5_rhs181);
  sf_mex_destroy(&c5_lhs181);
  sf_mex_destroy(&c5_rhs182);
  sf_mex_destroy(&c5_lhs182);
  sf_mex_destroy(&c5_rhs183);
  sf_mex_destroy(&c5_lhs183);
  sf_mex_destroy(&c5_rhs184);
  sf_mex_destroy(&c5_lhs184);
  sf_mex_destroy(&c5_rhs185);
  sf_mex_destroy(&c5_lhs185);
  sf_mex_destroy(&c5_rhs186);
  sf_mex_destroy(&c5_lhs186);
  sf_mex_destroy(&c5_rhs187);
  sf_mex_destroy(&c5_lhs187);
  sf_mex_destroy(&c5_rhs188);
  sf_mex_destroy(&c5_lhs188);
  sf_mex_destroy(&c5_rhs189);
  sf_mex_destroy(&c5_lhs189);
  sf_mex_destroy(&c5_rhs190);
  sf_mex_destroy(&c5_lhs190);
  sf_mex_destroy(&c5_rhs191);
  sf_mex_destroy(&c5_lhs191);
}

static void c5_d_info_helper(const mxArray **c5_info)
{
  const mxArray *c5_rhs192 = NULL;
  const mxArray *c5_lhs192 = NULL;
  const mxArray *c5_rhs193 = NULL;
  const mxArray *c5_lhs193 = NULL;
  const mxArray *c5_rhs194 = NULL;
  const mxArray *c5_lhs194 = NULL;
  const mxArray *c5_rhs195 = NULL;
  const mxArray *c5_lhs195 = NULL;
  const mxArray *c5_rhs196 = NULL;
  const mxArray *c5_lhs196 = NULL;
  const mxArray *c5_rhs197 = NULL;
  const mxArray *c5_lhs197 = NULL;
  const mxArray *c5_rhs198 = NULL;
  const mxArray *c5_lhs198 = NULL;
  const mxArray *c5_rhs199 = NULL;
  const mxArray *c5_lhs199 = NULL;
  const mxArray *c5_rhs200 = NULL;
  const mxArray *c5_lhs200 = NULL;
  const mxArray *c5_rhs201 = NULL;
  const mxArray *c5_lhs201 = NULL;
  const mxArray *c5_rhs202 = NULL;
  const mxArray *c5_lhs202 = NULL;
  const mxArray *c5_rhs203 = NULL;
  const mxArray *c5_lhs203 = NULL;
  const mxArray *c5_rhs204 = NULL;
  const mxArray *c5_lhs204 = NULL;
  const mxArray *c5_rhs205 = NULL;
  const mxArray *c5_lhs205 = NULL;
  const mxArray *c5_rhs206 = NULL;
  const mxArray *c5_lhs206 = NULL;
  const mxArray *c5_rhs207 = NULL;
  const mxArray *c5_lhs207 = NULL;
  const mxArray *c5_rhs208 = NULL;
  const mxArray *c5_lhs208 = NULL;
  const mxArray *c5_rhs209 = NULL;
  const mxArray *c5_lhs209 = NULL;
  const mxArray *c5_rhs210 = NULL;
  const mxArray *c5_lhs210 = NULL;
  const mxArray *c5_rhs211 = NULL;
  const mxArray *c5_lhs211 = NULL;
  const mxArray *c5_rhs212 = NULL;
  const mxArray *c5_lhs212 = NULL;
  const mxArray *c5_rhs213 = NULL;
  const mxArray *c5_lhs213 = NULL;
  const mxArray *c5_rhs214 = NULL;
  const mxArray *c5_lhs214 = NULL;
  const mxArray *c5_rhs215 = NULL;
  const mxArray *c5_lhs215 = NULL;
  const mxArray *c5_rhs216 = NULL;
  const mxArray *c5_lhs216 = NULL;
  const mxArray *c5_rhs217 = NULL;
  const mxArray *c5_lhs217 = NULL;
  const mxArray *c5_rhs218 = NULL;
  const mxArray *c5_lhs218 = NULL;
  const mxArray *c5_rhs219 = NULL;
  const mxArray *c5_lhs219 = NULL;
  const mxArray *c5_rhs220 = NULL;
  const mxArray *c5_lhs220 = NULL;
  const mxArray *c5_rhs221 = NULL;
  const mxArray *c5_lhs221 = NULL;
  const mxArray *c5_rhs222 = NULL;
  const mxArray *c5_lhs222 = NULL;
  const mxArray *c5_rhs223 = NULL;
  const mxArray *c5_lhs223 = NULL;
  const mxArray *c5_rhs224 = NULL;
  const mxArray *c5_lhs224 = NULL;
  const mxArray *c5_rhs225 = NULL;
  const mxArray *c5_lhs225 = NULL;
  const mxArray *c5_rhs226 = NULL;
  const mxArray *c5_lhs226 = NULL;
  const mxArray *c5_rhs227 = NULL;
  const mxArray *c5_lhs227 = NULL;
  const mxArray *c5_rhs228 = NULL;
  const mxArray *c5_lhs228 = NULL;
  const mxArray *c5_rhs229 = NULL;
  const mxArray *c5_lhs229 = NULL;
  const mxArray *c5_rhs230 = NULL;
  const mxArray *c5_lhs230 = NULL;
  const mxArray *c5_rhs231 = NULL;
  const mxArray *c5_lhs231 = NULL;
  const mxArray *c5_rhs232 = NULL;
  const mxArray *c5_lhs232 = NULL;
  const mxArray *c5_rhs233 = NULL;
  const mxArray *c5_lhs233 = NULL;
  const mxArray *c5_rhs234 = NULL;
  const mxArray *c5_lhs234 = NULL;
  const mxArray *c5_rhs235 = NULL;
  const mxArray *c5_lhs235 = NULL;
  const mxArray *c5_rhs236 = NULL;
  const mxArray *c5_lhs236 = NULL;
  const mxArray *c5_rhs237 = NULL;
  const mxArray *c5_lhs237 = NULL;
  const mxArray *c5_rhs238 = NULL;
  const mxArray *c5_lhs238 = NULL;
  const mxArray *c5_rhs239 = NULL;
  const mxArray *c5_lhs239 = NULL;
  const mxArray *c5_rhs240 = NULL;
  const mxArray *c5_lhs240 = NULL;
  const mxArray *c5_rhs241 = NULL;
  const mxArray *c5_lhs241 = NULL;
  const mxArray *c5_rhs242 = NULL;
  const mxArray *c5_lhs242 = NULL;
  const mxArray *c5_rhs243 = NULL;
  const mxArray *c5_lhs243 = NULL;
  const mxArray *c5_rhs244 = NULL;
  const mxArray *c5_lhs244 = NULL;
  const mxArray *c5_rhs245 = NULL;
  const mxArray *c5_lhs245 = NULL;
  const mxArray *c5_rhs246 = NULL;
  const mxArray *c5_lhs246 = NULL;
  const mxArray *c5_rhs247 = NULL;
  const mxArray *c5_lhs247 = NULL;
  const mxArray *c5_rhs248 = NULL;
  const mxArray *c5_lhs248 = NULL;
  const mxArray *c5_rhs249 = NULL;
  const mxArray *c5_lhs249 = NULL;
  const mxArray *c5_rhs250 = NULL;
  const mxArray *c5_lhs250 = NULL;
  const mxArray *c5_rhs251 = NULL;
  const mxArray *c5_lhs251 = NULL;
  const mxArray *c5_rhs252 = NULL;
  const mxArray *c5_lhs252 = NULL;
  const mxArray *c5_rhs253 = NULL;
  const mxArray *c5_lhs253 = NULL;
  const mxArray *c5_rhs254 = NULL;
  const mxArray *c5_lhs254 = NULL;
  const mxArray *c5_rhs255 = NULL;
  const mxArray *c5_lhs255 = NULL;
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m"), "context",
                  "context", 192);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_exp"), "name",
                  "name", 192);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 192);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_exp.m"),
                  "resolved", "resolved", 192);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1301335664U), "fileTimeLo",
                  "fileTimeLo", 192);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 192);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 192);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 192);
  sf_mex_assign(&c5_rhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs192), "rhs", "rhs",
                  192);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs192), "lhs", "lhs",
                  192);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 193);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 193);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 193);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 193);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 193);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 193);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 193);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 193);
  sf_mex_assign(&c5_rhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs193), "rhs", "rhs",
                  193);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs193), "lhs", "lhs",
                  193);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 194);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 194);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 194);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 194);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 194);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 194);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 194);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 194);
  sf_mex_assign(&c5_rhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs194), "rhs", "rhs",
                  194);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs194), "lhs", "lhs",
                  194);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 195);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 195);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 195);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 195);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 195);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 195);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 195);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 195);
  sf_mex_assign(&c5_rhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs195), "rhs", "rhs",
                  195);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs195), "lhs", "lhs",
                  195);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 196);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 196);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 196);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 196);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 196);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 196);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 196);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 196);
  sf_mex_assign(&c5_rhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs196), "rhs", "rhs",
                  196);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs196), "lhs", "lhs",
                  196);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 197);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 197);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 197);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 197);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 197);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 197);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 197);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 197);
  sf_mex_assign(&c5_rhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs197), "rhs", "rhs",
                  197);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs197), "lhs", "lhs",
                  197);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 198);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 198);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 198);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 198);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 198);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 198);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 198);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 198);
  sf_mex_assign(&c5_rhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs198), "rhs", "rhs",
                  198);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs198), "lhs", "lhs",
                  198);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 199);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 199);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 199);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 199);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 199);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 199);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 199);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 199);
  sf_mex_assign(&c5_rhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs199), "rhs", "rhs",
                  199);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs199), "lhs", "lhs",
                  199);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 200);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  200);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 200);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 200);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 200);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 200);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 200);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 200);
  sf_mex_assign(&c5_rhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs200), "rhs", "rhs",
                  200);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs200), "lhs", "lhs",
                  200);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 201);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 201);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 201);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 201);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 201);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 201);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 201);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 201);
  sf_mex_assign(&c5_rhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs201), "rhs", "rhs",
                  201);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs201), "lhs", "lhs",
                  201);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 202);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 202);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 202);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 202);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 202);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 202);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 202);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 202);
  sf_mex_assign(&c5_rhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs202), "rhs", "rhs",
                  202);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs202), "lhs", "lhs",
                  202);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 203);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 203);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 203);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 203);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 203);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 203);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 203);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 203);
  sf_mex_assign(&c5_rhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs203), "rhs", "rhs",
                  203);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs203), "lhs", "lhs",
                  203);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 204);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 204);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 204);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 204);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 204);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 204);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 204);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 204);
  sf_mex_assign(&c5_rhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs204), "rhs", "rhs",
                  204);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs204), "lhs", "lhs",
                  204);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 205);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 205);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 205);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 205);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 205);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 205);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 205);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 205);
  sf_mex_assign(&c5_rhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs205), "rhs", "rhs",
                  205);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs205), "lhs", "lhs",
                  205);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 206);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 206);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 206);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 206);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 206);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 206);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 206);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 206);
  sf_mex_assign(&c5_rhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs206), "rhs", "rhs",
                  206);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs206), "lhs", "lhs",
                  206);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 207);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("inv"), "name", "name", 207);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 207);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 207);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 207);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 207);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 207);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 207);
  sf_mex_assign(&c5_rhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs207), "rhs", "rhs",
                  207);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs207), "lhs", "lhs",
                  207);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 208);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 208);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 208);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 208);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 208);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 208);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 208);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 208);
  sf_mex_assign(&c5_rhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs208), "rhs", "rhs",
                  208);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs208), "lhs", "lhs",
                  208);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 209);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_error"), "name", "name",
                  209);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 209);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 209);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 209);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 209);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 209);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 209);
  sf_mex_assign(&c5_rhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs209), "rhs", "rhs",
                  209);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs209), "lhs", "lhs",
                  209);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 210);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eig"), "name", "name", 210);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 210);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "resolved",
                  "resolved", 210);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 210);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 210);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 210);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 210);
  sf_mex_assign(&c5_rhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs210), "rhs", "rhs",
                  210);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs210), "lhs", "lhs",
                  210);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 211);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 211);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 211);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 211);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 211);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 211);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 211);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 211);
  sf_mex_assign(&c5_rhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs211), "rhs", "rhs",
                  211);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs211), "lhs", "lhs",
                  211);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 212);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eye"), "name", "name", 212);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 212);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 212);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 212);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 212);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 212);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 212);
  sf_mex_assign(&c5_rhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs212), "rhs", "rhs",
                  212);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs212), "lhs", "lhs",
                  212);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 213);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 213);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 213);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 213);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 213);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 213);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 213);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 213);
  sf_mex_assign(&c5_rhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs213), "rhs", "rhs",
                  213);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs213), "lhs", "lhs",
                  213);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 214);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 214);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 214);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 214);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 214);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 214);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 214);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 214);
  sf_mex_assign(&c5_rhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs214), "rhs", "rhs",
                  214);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs214), "lhs", "lhs",
                  214);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 215);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 215);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 215);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 215);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 215);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 215);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 215);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 215);
  sf_mex_assign(&c5_rhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs215), "rhs", "rhs",
                  215);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs215), "lhs", "lhs",
                  215);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 216);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("intmax"), "name", "name", 216);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 216);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 216);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 216);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 216);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 216);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 216);
  sf_mex_assign(&c5_rhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs216), "rhs", "rhs",
                  216);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs216), "lhs", "lhs",
                  216);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 217);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("intmin"), "name", "name", 217);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 217);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 217);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 217);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 217);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 217);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 217);
  sf_mex_assign(&c5_rhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs217), "rhs", "rhs",
                  217);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs217), "lhs", "lhs",
                  217);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 218);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 218);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 218);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 218);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 218);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 218);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 218);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 218);
  sf_mex_assign(&c5_rhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs218), "rhs", "rhs",
                  218);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs218), "lhs", "lhs",
                  218);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 219);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 219);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 219);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 219);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 219);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 219);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 219);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 219);
  sf_mex_assign(&c5_rhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs219), "rhs", "rhs",
                  219);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs219), "lhs", "lhs",
                  219);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 220);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("intmax"), "name", "name", 220);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 220);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 220);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 220);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 220);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 220);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 220);
  sf_mex_assign(&c5_rhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs220), "rhs", "rhs",
                  220);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs220), "lhs", "lhs",
                  220);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 221);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 221);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 221);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 221);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 221);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 221);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 221);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 221);
  sf_mex_assign(&c5_rhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs221), "rhs", "rhs",
                  221);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs221), "lhs", "lhs",
                  221);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 222);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 222);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 222);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 222);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 222);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 222);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 222);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 222);
  sf_mex_assign(&c5_rhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs222), "rhs", "rhs",
                  222);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs222), "lhs", "lhs",
                  222);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 223);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_ztgevc"), "name",
                  "name", 223);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 223);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "resolved", "resolved", 223);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826024U), "fileTimeLo",
                  "fileTimeLo", 223);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 223);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 223);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 223);
  sf_mex_assign(&c5_rhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs223), "rhs", "rhs",
                  223);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs223), "lhs", "lhs",
                  223);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 224);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eps"), "name", "name", 224);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 224);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 224);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 224);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 224);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 224);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 224);
  sf_mex_assign(&c5_rhs224, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs224, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs224), "rhs", "rhs",
                  224);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs224), "lhs", "lhs",
                  224);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 225);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("realmin"), "name", "name", 225);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 225);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 225);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 225);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 225);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 225);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 225);
  sf_mex_assign(&c5_rhs225, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs225, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs225), "rhs", "rhs",
                  225);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs225), "lhs", "lhs",
                  225);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m!abs1"),
                  "context", "context", 226);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 226);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 226);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 226);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 226);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 226);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 226);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 226);
  sf_mex_assign(&c5_rhs226, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs226, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs226), "rhs", "rhs",
                  226);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs226), "lhs", "lhs",
                  226);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 227);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 227);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 227);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 227);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 227);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 227);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 227);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 227);
  sf_mex_assign(&c5_rhs227, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs227, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs227), "rhs", "rhs",
                  227);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs227), "lhs", "lhs",
                  227);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 228);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 228);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 228);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 228);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 228);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 228);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 228);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 228);
  sf_mex_assign(&c5_rhs228, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs228, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs228), "rhs", "rhs",
                  228);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs228), "lhs", "lhs",
                  228);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 229);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mrdivide"), "name", "name",
                  229);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 229);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 229);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 229);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 229);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 229);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 229);
  sf_mex_assign(&c5_rhs229, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs229, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs229), "rhs", "rhs",
                  229);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs229), "lhs", "lhs",
                  229);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 230);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_matlab_zggbak"), "name",
                  "name", 230);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 230);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "resolved", "resolved", 230);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826018U), "fileTimeLo",
                  "fileTimeLo", 230);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 230);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 230);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 230);
  sf_mex_assign(&c5_rhs230, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs230, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs230), "rhs", "rhs",
                  230);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs230), "lhs", "lhs",
                  230);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 231);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 231);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 231);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 231);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 231);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 231);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 231);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 231);
  sf_mex_assign(&c5_rhs231, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs231, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs231), "rhs", "rhs",
                  231);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs231), "lhs", "lhs",
                  231);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 232);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 232);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 232);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 232);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 232);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 232);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 232);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 232);
  sf_mex_assign(&c5_rhs232, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs232, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs232), "rhs", "rhs",
                  232);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs232), "lhs", "lhs",
                  232);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 233);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 233);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 233);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 233);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 233);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 233);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 233);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 233);
  sf_mex_assign(&c5_rhs233, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs233, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs233), "rhs", "rhs",
                  233);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs233), "lhs", "lhs",
                  233);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 234);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 234);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 234);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 234);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 234);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 234);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 234);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 234);
  sf_mex_assign(&c5_rhs234, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs234, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs234), "rhs", "rhs",
                  234);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs234), "lhs", "lhs",
                  234);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m!abs1"),
                  "context", "context", 235);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 235);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 235);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 235);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 235);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 235);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 235);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 235);
  sf_mex_assign(&c5_rhs235, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs235, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs235), "rhs", "rhs",
                  235);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs235), "lhs", "lhs",
                  235);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 236);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 236);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 236);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 236);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 236);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 236);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 236);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 236);
  sf_mex_assign(&c5_rhs236, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs236, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs236), "rhs", "rhs",
                  236);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs236), "lhs", "lhs",
                  236);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 237);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 237);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 237);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 237);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 237);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 237);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 237);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 237);
  sf_mex_assign(&c5_rhs237, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs237, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs237), "rhs", "rhs",
                  237);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs237), "lhs", "lhs",
                  237);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 238);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 238);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 238);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 238);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 238);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 238);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 238);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 238);
  sf_mex_assign(&c5_rhs238, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs238, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs238), "rhs", "rhs",
                  238);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs238), "lhs", "lhs",
                  238);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 239);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 239);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 239);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 239);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 239);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 239);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 239);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 239);
  sf_mex_assign(&c5_rhs239, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs239, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs239), "rhs", "rhs",
                  239);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs239), "lhs", "lhs",
                  239);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 240);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 240);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 240);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 240);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 240);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 240);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 240);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 240);
  sf_mex_assign(&c5_rhs240, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs240, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs240), "rhs", "rhs",
                  240);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs240), "lhs", "lhs",
                  240);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 241);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 241);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 241);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 241);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 241);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 241);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 241);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 241);
  sf_mex_assign(&c5_rhs241, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs241, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs241), "rhs", "rhs",
                  241);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs241), "lhs", "lhs",
                  241);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 242);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 242);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 242);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 242);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 242);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 242);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 242);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 242);
  sf_mex_assign(&c5_rhs242, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs242, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs242), "rhs", "rhs",
                  242);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs242), "lhs", "lhs",
                  242);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 243);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 243);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 243);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 243);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 243);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 243);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 243);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 243);
  sf_mex_assign(&c5_rhs243, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs243, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs243), "rhs", "rhs",
                  243);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs243), "lhs", "lhs",
                  243);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 244);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  244);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 244);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 244);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 244);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 244);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 244);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 244);
  sf_mex_assign(&c5_rhs244, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs244, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs244), "rhs", "rhs",
                  244);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs244), "lhs", "lhs",
                  244);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 245);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mrdivide"), "name", "name",
                  245);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 245);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 245);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 245);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 245);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 245);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 245);
  sf_mex_assign(&c5_rhs245, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs245, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs245), "rhs", "rhs",
                  245);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs245), "lhs", "lhs",
                  245);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 246);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("diag"), "name", "name", 246);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 246);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "resolved",
                  "resolved", 246);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 246);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 246);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 246);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 246);
  sf_mex_assign(&c5_rhs246, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs246, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs246), "rhs", "rhs",
                  246);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs246), "lhs", "lhs",
                  246);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 247);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("ismatrix"), "name", "name",
                  247);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 247);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 247);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 247);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 247);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 247);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 247);
  sf_mex_assign(&c5_rhs247, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs247, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs247), "rhs", "rhs",
                  247);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs247), "lhs", "lhs",
                  247);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 248);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 248);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 248);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 248);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 248);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 248);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 248);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 248);
  sf_mex_assign(&c5_rhs248, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs248, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs248), "rhs", "rhs",
                  248);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs248), "lhs", "lhs",
                  248);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 249);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 249);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 249);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 249);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 249);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 249);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 249);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 249);
  sf_mex_assign(&c5_rhs249, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs249, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs249), "rhs", "rhs",
                  249);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs249), "lhs", "lhs",
                  249);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 250);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 250);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 250);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 250);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 250);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 250);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 250);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 250);
  sf_mex_assign(&c5_rhs250, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs250, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs250), "rhs", "rhs",
                  250);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs250), "lhs", "lhs",
                  250);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 251);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("power"), "name", "name", 251);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 251);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 251);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 251);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 251);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 251);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 251);
  sf_mex_assign(&c5_rhs251, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs251, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs251), "rhs", "rhs",
                  251);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs251), "lhs", "lhs",
                  251);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 252);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("flintmax"), "name", "name",
                  252);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 252);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/flintmax.m"), "resolved",
                  "resolved", 252);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1348199116U), "fileTimeLo",
                  "fileTimeLo", 252);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 252);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 252);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 252);
  sf_mex_assign(&c5_rhs252, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs252, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs252), "rhs", "rhs",
                  252);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs252), "lhs", "lhs",
                  252);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/flintmax.m"), "context",
                  "context", 253);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 253);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 253);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 253);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 253);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 253);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 253);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 253);
  sf_mex_assign(&c5_rhs253, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs253, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs253), "rhs", "rhs",
                  253);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs253), "lhs", "lhs",
                  253);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 254);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mod"), "name", "name", 254);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 254);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "resolved",
                  "resolved", 254);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 254);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 254);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 254);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 254);
  sf_mex_assign(&c5_rhs254, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs254, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs254), "rhs", "rhs",
                  254);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs254), "lhs", "lhs",
                  254);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 255);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 255);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 255);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 255);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 255);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 255);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 255);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 255);
  sf_mex_assign(&c5_rhs255, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs255, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs255), "rhs", "rhs",
                  255);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs255), "lhs", "lhs",
                  255);
  sf_mex_destroy(&c5_rhs192);
  sf_mex_destroy(&c5_lhs192);
  sf_mex_destroy(&c5_rhs193);
  sf_mex_destroy(&c5_lhs193);
  sf_mex_destroy(&c5_rhs194);
  sf_mex_destroy(&c5_lhs194);
  sf_mex_destroy(&c5_rhs195);
  sf_mex_destroy(&c5_lhs195);
  sf_mex_destroy(&c5_rhs196);
  sf_mex_destroy(&c5_lhs196);
  sf_mex_destroy(&c5_rhs197);
  sf_mex_destroy(&c5_lhs197);
  sf_mex_destroy(&c5_rhs198);
  sf_mex_destroy(&c5_lhs198);
  sf_mex_destroy(&c5_rhs199);
  sf_mex_destroy(&c5_lhs199);
  sf_mex_destroy(&c5_rhs200);
  sf_mex_destroy(&c5_lhs200);
  sf_mex_destroy(&c5_rhs201);
  sf_mex_destroy(&c5_lhs201);
  sf_mex_destroy(&c5_rhs202);
  sf_mex_destroy(&c5_lhs202);
  sf_mex_destroy(&c5_rhs203);
  sf_mex_destroy(&c5_lhs203);
  sf_mex_destroy(&c5_rhs204);
  sf_mex_destroy(&c5_lhs204);
  sf_mex_destroy(&c5_rhs205);
  sf_mex_destroy(&c5_lhs205);
  sf_mex_destroy(&c5_rhs206);
  sf_mex_destroy(&c5_lhs206);
  sf_mex_destroy(&c5_rhs207);
  sf_mex_destroy(&c5_lhs207);
  sf_mex_destroy(&c5_rhs208);
  sf_mex_destroy(&c5_lhs208);
  sf_mex_destroy(&c5_rhs209);
  sf_mex_destroy(&c5_lhs209);
  sf_mex_destroy(&c5_rhs210);
  sf_mex_destroy(&c5_lhs210);
  sf_mex_destroy(&c5_rhs211);
  sf_mex_destroy(&c5_lhs211);
  sf_mex_destroy(&c5_rhs212);
  sf_mex_destroy(&c5_lhs212);
  sf_mex_destroy(&c5_rhs213);
  sf_mex_destroy(&c5_lhs213);
  sf_mex_destroy(&c5_rhs214);
  sf_mex_destroy(&c5_lhs214);
  sf_mex_destroy(&c5_rhs215);
  sf_mex_destroy(&c5_lhs215);
  sf_mex_destroy(&c5_rhs216);
  sf_mex_destroy(&c5_lhs216);
  sf_mex_destroy(&c5_rhs217);
  sf_mex_destroy(&c5_lhs217);
  sf_mex_destroy(&c5_rhs218);
  sf_mex_destroy(&c5_lhs218);
  sf_mex_destroy(&c5_rhs219);
  sf_mex_destroy(&c5_lhs219);
  sf_mex_destroy(&c5_rhs220);
  sf_mex_destroy(&c5_lhs220);
  sf_mex_destroy(&c5_rhs221);
  sf_mex_destroy(&c5_lhs221);
  sf_mex_destroy(&c5_rhs222);
  sf_mex_destroy(&c5_lhs222);
  sf_mex_destroy(&c5_rhs223);
  sf_mex_destroy(&c5_lhs223);
  sf_mex_destroy(&c5_rhs224);
  sf_mex_destroy(&c5_lhs224);
  sf_mex_destroy(&c5_rhs225);
  sf_mex_destroy(&c5_lhs225);
  sf_mex_destroy(&c5_rhs226);
  sf_mex_destroy(&c5_lhs226);
  sf_mex_destroy(&c5_rhs227);
  sf_mex_destroy(&c5_lhs227);
  sf_mex_destroy(&c5_rhs228);
  sf_mex_destroy(&c5_lhs228);
  sf_mex_destroy(&c5_rhs229);
  sf_mex_destroy(&c5_lhs229);
  sf_mex_destroy(&c5_rhs230);
  sf_mex_destroy(&c5_lhs230);
  sf_mex_destroy(&c5_rhs231);
  sf_mex_destroy(&c5_lhs231);
  sf_mex_destroy(&c5_rhs232);
  sf_mex_destroy(&c5_lhs232);
  sf_mex_destroy(&c5_rhs233);
  sf_mex_destroy(&c5_lhs233);
  sf_mex_destroy(&c5_rhs234);
  sf_mex_destroy(&c5_lhs234);
  sf_mex_destroy(&c5_rhs235);
  sf_mex_destroy(&c5_lhs235);
  sf_mex_destroy(&c5_rhs236);
  sf_mex_destroy(&c5_lhs236);
  sf_mex_destroy(&c5_rhs237);
  sf_mex_destroy(&c5_lhs237);
  sf_mex_destroy(&c5_rhs238);
  sf_mex_destroy(&c5_lhs238);
  sf_mex_destroy(&c5_rhs239);
  sf_mex_destroy(&c5_lhs239);
  sf_mex_destroy(&c5_rhs240);
  sf_mex_destroy(&c5_lhs240);
  sf_mex_destroy(&c5_rhs241);
  sf_mex_destroy(&c5_lhs241);
  sf_mex_destroy(&c5_rhs242);
  sf_mex_destroy(&c5_lhs242);
  sf_mex_destroy(&c5_rhs243);
  sf_mex_destroy(&c5_lhs243);
  sf_mex_destroy(&c5_rhs244);
  sf_mex_destroy(&c5_lhs244);
  sf_mex_destroy(&c5_rhs245);
  sf_mex_destroy(&c5_lhs245);
  sf_mex_destroy(&c5_rhs246);
  sf_mex_destroy(&c5_lhs246);
  sf_mex_destroy(&c5_rhs247);
  sf_mex_destroy(&c5_lhs247);
  sf_mex_destroy(&c5_rhs248);
  sf_mex_destroy(&c5_lhs248);
  sf_mex_destroy(&c5_rhs249);
  sf_mex_destroy(&c5_lhs249);
  sf_mex_destroy(&c5_rhs250);
  sf_mex_destroy(&c5_lhs250);
  sf_mex_destroy(&c5_rhs251);
  sf_mex_destroy(&c5_lhs251);
  sf_mex_destroy(&c5_rhs252);
  sf_mex_destroy(&c5_lhs252);
  sf_mex_destroy(&c5_rhs253);
  sf_mex_destroy(&c5_lhs253);
  sf_mex_destroy(&c5_rhs254);
  sf_mex_destroy(&c5_lhs254);
  sf_mex_destroy(&c5_rhs255);
  sf_mex_destroy(&c5_lhs255);
}

static void c5_e_info_helper(const mxArray **c5_info)
{
  const mxArray *c5_rhs256 = NULL;
  const mxArray *c5_lhs256 = NULL;
  const mxArray *c5_rhs257 = NULL;
  const mxArray *c5_lhs257 = NULL;
  const mxArray *c5_rhs258 = NULL;
  const mxArray *c5_lhs258 = NULL;
  const mxArray *c5_rhs259 = NULL;
  const mxArray *c5_lhs259 = NULL;
  const mxArray *c5_rhs260 = NULL;
  const mxArray *c5_lhs260 = NULL;
  const mxArray *c5_rhs261 = NULL;
  const mxArray *c5_lhs261 = NULL;
  const mxArray *c5_rhs262 = NULL;
  const mxArray *c5_lhs262 = NULL;
  const mxArray *c5_rhs263 = NULL;
  const mxArray *c5_lhs263 = NULL;
  const mxArray *c5_rhs264 = NULL;
  const mxArray *c5_lhs264 = NULL;
  const mxArray *c5_rhs265 = NULL;
  const mxArray *c5_lhs265 = NULL;
  const mxArray *c5_rhs266 = NULL;
  const mxArray *c5_lhs266 = NULL;
  const mxArray *c5_rhs267 = NULL;
  const mxArray *c5_lhs267 = NULL;
  const mxArray *c5_rhs268 = NULL;
  const mxArray *c5_lhs268 = NULL;
  const mxArray *c5_rhs269 = NULL;
  const mxArray *c5_lhs269 = NULL;
  const mxArray *c5_rhs270 = NULL;
  const mxArray *c5_lhs270 = NULL;
  const mxArray *c5_rhs271 = NULL;
  const mxArray *c5_lhs271 = NULL;
  const mxArray *c5_rhs272 = NULL;
  const mxArray *c5_lhs272 = NULL;
  const mxArray *c5_rhs273 = NULL;
  const mxArray *c5_lhs273 = NULL;
  const mxArray *c5_rhs274 = NULL;
  const mxArray *c5_lhs274 = NULL;
  const mxArray *c5_rhs275 = NULL;
  const mxArray *c5_lhs275 = NULL;
  const mxArray *c5_rhs276 = NULL;
  const mxArray *c5_lhs276 = NULL;
  const mxArray *c5_rhs277 = NULL;
  const mxArray *c5_lhs277 = NULL;
  const mxArray *c5_rhs278 = NULL;
  const mxArray *c5_lhs278 = NULL;
  const mxArray *c5_rhs279 = NULL;
  const mxArray *c5_lhs279 = NULL;
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 256);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 256);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 256);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 256);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 256);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 256);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 256);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 256);
  sf_mex_assign(&c5_rhs256, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs256, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs256), "rhs", "rhs",
                  256);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs256), "lhs", "lhs",
                  256);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 257);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_round"), "name",
                  "name", 257);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 257);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_round.m"),
                  "resolved", "resolved", 257);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658438U), "fileTimeLo",
                  "fileTimeLo", 257);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 257);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 257);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 257);
  sf_mex_assign(&c5_rhs257, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs257, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs257), "rhs", "rhs",
                  257);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs257), "lhs", "lhs",
                  257);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 258);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 258);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 258);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 258);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 258);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 258);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 258);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 258);
  sf_mex_assign(&c5_rhs258, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs258, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs258), "rhs", "rhs",
                  258);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs258), "lhs", "lhs",
                  258);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 259);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eps"), "name", "name", 259);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 259);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 259);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 259);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 259);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 259);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 259);
  sf_mex_assign(&c5_rhs259, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs259, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs259), "rhs", "rhs",
                  259);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs259), "lhs", "lhs",
                  259);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 260);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("log"), "name", "name", 260);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 260);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m"), "resolved",
                  "resolved", 260);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1343837580U), "fileTimeLo",
                  "fileTimeLo", 260);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 260);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 260);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 260);
  sf_mex_assign(&c5_rhs260, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs260, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs260), "rhs", "rhs",
                  260);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs260), "lhs", "lhs",
                  260);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m"), "context",
                  "context", 261);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_log"), "name",
                  "name", 261);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 261);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "resolved", "resolved", 261);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825928U), "fileTimeLo",
                  "fileTimeLo", 261);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 261);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 261);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 261);
  sf_mex_assign(&c5_rhs261, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs261, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs261), "rhs", "rhs",
                  261);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs261), "lhs", "lhs",
                  261);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 262);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("realmax"), "name", "name", 262);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 262);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m"), "resolved",
                  "resolved", 262);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 262);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 262);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 262);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 262);
  sf_mex_assign(&c5_rhs262, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs262, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs262), "rhs", "rhs",
                  262);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs262), "lhs", "lhs",
                  262);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 263);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mrdivide"), "name", "name",
                  263);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 263);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 263);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 263);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 263);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 263);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 263);
  sf_mex_assign(&c5_rhs263, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs263, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs263), "rhs", "rhs",
                  263);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs263), "lhs", "lhs",
                  263);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 264);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("isnan"), "name", "name", 264);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 264);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 264);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 264);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 264);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 264);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 264);
  sf_mex_assign(&c5_rhs264, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs264, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs264), "rhs", "rhs",
                  264);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs264), "lhs", "lhs",
                  264);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 265);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_hypot"), "name",
                  "name", 265);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 265);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "resolved", "resolved", 265);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 265);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 265);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 265);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 265);
  sf_mex_assign(&c5_rhs265, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs265, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs265), "rhs", "rhs",
                  265);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs265), "lhs", "lhs",
                  265);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 266);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_atan2"), "name",
                  "name", 266);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 266);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m"),
                  "resolved", "resolved", 266);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825920U), "fileTimeLo",
                  "fileTimeLo", 266);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 266);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 266);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 266);
  sf_mex_assign(&c5_rhs266, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs266, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs266), "rhs", "rhs",
                  266);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs266), "lhs", "lhs",
                  266);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 267);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 267);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 267);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 267);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 267);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 267);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 267);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 267);
  sf_mex_assign(&c5_rhs267, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs267, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs267), "rhs", "rhs",
                  267);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs267), "lhs", "lhs",
                  267);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 268);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 268);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 268);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 268);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 268);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 268);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 268);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 268);
  sf_mex_assign(&c5_rhs268, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs268, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs268), "rhs", "rhs",
                  268);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs268), "lhs", "lhs",
                  268);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 269);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("mrdivide"), "name", "name",
                  269);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 269);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 269);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 269);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 269);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 269);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 269);
  sf_mex_assign(&c5_rhs269, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs269, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs269), "rhs", "rhs",
                  269);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs269), "lhs", "lhs",
                  269);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 270);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("ismatrix"), "name", "name",
                  270);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 270);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 270);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 270);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 270);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 270);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 270);
  sf_mex_assign(&c5_rhs270, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs270, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs270), "rhs", "rhs",
                  270);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs270), "lhs", "lhs",
                  270);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 271);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_lusolve"), "name", "name",
                  271);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 271);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "resolved",
                  "resolved", 271);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1370017086U), "fileTimeLo",
                  "fileTimeLo", 271);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 271);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 271);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 271);
  sf_mex_assign(&c5_rhs271, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs271, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs271), "rhs", "rhs",
                  271);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs271), "lhs", "lhs",
                  271);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 272);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_xcabs1"), "name", "name",
                  272);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 272);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "resolved", "resolved", 272);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 272);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 272);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 272);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 272);
  sf_mex_assign(&c5_rhs272, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs272, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs272), "rhs", "rhs",
                  272);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs272), "lhs", "lhs",
                  272);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "context", "context", 273);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("coder.internal.refblas.xcabs1"),
                  "name", "name", 273);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 273);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "resolved", "resolved", 273);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 273);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 273);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 273);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 273);
  sf_mex_assign(&c5_rhs273, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs273, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs273), "rhs", "rhs",
                  273);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs273), "lhs", "lhs",
                  273);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "context", "context", 274);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("abs"), "name", "name", 274);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 274);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 274);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 274);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 274);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 274);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 274);
  sf_mex_assign(&c5_rhs274, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs274, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs274), "rhs", "rhs",
                  274);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs274), "lhs", "lhs",
                  274);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 275);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("rdivide"), "name", "name", 275);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 275);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 275);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 275);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 275);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 275);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 275);
  sf_mex_assign(&c5_rhs275, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs275, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs275), "rhs", "rhs",
                  275);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs275), "lhs", "lhs",
                  275);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular"),
                  "context", "context", 276);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_warning"), "name", "name",
                  276);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 276);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 276);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 276);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 276);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 276);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 276);
  sf_mex_assign(&c5_rhs276, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs276, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs276), "rhs", "rhs",
                  276);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs276), "lhs", "lhs",
                  276);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 277);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 277);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 277);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 277);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 277);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 277);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 277);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 277);
  sf_mex_assign(&c5_rhs277, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs277, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs277), "rhs", "rhs",
                  277);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs277), "lhs", "lhs",
                  277);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 278);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 278);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 278);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 278);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 278);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 278);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 278);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 278);
  sf_mex_assign(&c5_rhs278, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs278, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs278), "rhs", "rhs",
                  278);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs278), "lhs", "lhs",
                  278);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "context", "context", 279);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 279);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 279);
  sf_mex_addfield(*c5_info, c5_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 279);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 279);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 279);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 279);
  sf_mex_addfield(*c5_info, c5_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 279);
  sf_mex_assign(&c5_rhs279, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c5_lhs279, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_rhs279), "rhs", "rhs",
                  279);
  sf_mex_addfield(*c5_info, sf_mex_duplicatearraysafe(&c5_lhs279), "lhs", "lhs",
                  279);
  sf_mex_destroy(&c5_rhs256);
  sf_mex_destroy(&c5_lhs256);
  sf_mex_destroy(&c5_rhs257);
  sf_mex_destroy(&c5_lhs257);
  sf_mex_destroy(&c5_rhs258);
  sf_mex_destroy(&c5_lhs258);
  sf_mex_destroy(&c5_rhs259);
  sf_mex_destroy(&c5_lhs259);
  sf_mex_destroy(&c5_rhs260);
  sf_mex_destroy(&c5_lhs260);
  sf_mex_destroy(&c5_rhs261);
  sf_mex_destroy(&c5_lhs261);
  sf_mex_destroy(&c5_rhs262);
  sf_mex_destroy(&c5_lhs262);
  sf_mex_destroy(&c5_rhs263);
  sf_mex_destroy(&c5_lhs263);
  sf_mex_destroy(&c5_rhs264);
  sf_mex_destroy(&c5_lhs264);
  sf_mex_destroy(&c5_rhs265);
  sf_mex_destroy(&c5_lhs265);
  sf_mex_destroy(&c5_rhs266);
  sf_mex_destroy(&c5_lhs266);
  sf_mex_destroy(&c5_rhs267);
  sf_mex_destroy(&c5_lhs267);
  sf_mex_destroy(&c5_rhs268);
  sf_mex_destroy(&c5_lhs268);
  sf_mex_destroy(&c5_rhs269);
  sf_mex_destroy(&c5_lhs269);
  sf_mex_destroy(&c5_rhs270);
  sf_mex_destroy(&c5_lhs270);
  sf_mex_destroy(&c5_rhs271);
  sf_mex_destroy(&c5_lhs271);
  sf_mex_destroy(&c5_rhs272);
  sf_mex_destroy(&c5_lhs272);
  sf_mex_destroy(&c5_rhs273);
  sf_mex_destroy(&c5_lhs273);
  sf_mex_destroy(&c5_rhs274);
  sf_mex_destroy(&c5_lhs274);
  sf_mex_destroy(&c5_rhs275);
  sf_mex_destroy(&c5_lhs275);
  sf_mex_destroy(&c5_rhs276);
  sf_mex_destroy(&c5_lhs276);
  sf_mex_destroy(&c5_rhs277);
  sf_mex_destroy(&c5_lhs277);
  sf_mex_destroy(&c5_rhs278);
  sf_mex_destroy(&c5_lhs278);
  sf_mex_destroy(&c5_rhs279);
  sf_mex_destroy(&c5_lhs279);
}

static real_T c5_norm(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x[3])
{
  real_T c5_y;
  real_T c5_scale;
  int32_T c5_k;
  int32_T c5_b_k;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_absxk;
  real_T c5_t;
  c5_eml_switch_helper(chartInstance);
  c5_y = 0.0;
  c5_realmin(chartInstance);
  c5_scale = 2.2250738585072014E-308;
  for (c5_k = 1; c5_k < 4; c5_k++) {
    c5_b_k = c5_k;
    c5_b_x = c5_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_b_k), 1, 3, 1, 0) - 1];
    c5_c_x = c5_b_x;
    c5_absxk = muDoubleScalarAbs(c5_c_x);
    if (c5_absxk > c5_scale) {
      c5_t = c5_scale / c5_absxk;
      c5_y = 1.0 + c5_y * c5_t * c5_t;
      c5_scale = c5_absxk;
    } else {
      c5_t = c5_absxk / c5_scale;
      c5_y += c5_t * c5_t;
    }
  }

  return c5_scale * muDoubleScalarSqrt(c5_y);
}

static void c5_eml_switch_helper(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c5_realmin(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c5_b_eml_switch_helper(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c5_abs(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x)
{
  real_T c5_b_x;
  (void)chartInstance;
  c5_b_x = c5_x;
  return muDoubleScalarAbs(c5_b_x);
}

static void c5_eig(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_A[9],
                   creal_T c5_V[3])
{
  int32_T c5_i43;
  static creal_T c5_dc0 = { 0.0, 0.0 };

  creal_T c5_b_A[9];
  real_T c5_info;
  int32_T c5_i44;
  creal_T c5_c_A[9];
  real_T c5_anrm;
  int32_T c5_i45;
  creal_T c5_alpha1[3];
  int32_T c5_i46;
  creal_T c5_beta1[3];
  boolean_T c5_ilascl;
  real_T c5_anrmto;
  int32_T c5_rscale[3];
  int32_T c5_ihi;
  int32_T c5_ilo;
  int32_T c5_b_ilo;
  int32_T c5_b_ihi;
  int32_T c5_c_ilo;
  int32_T c5_c_ihi;
  int32_T c5_a;
  int32_T c5_b_a;
  int32_T c5_c;
  int32_T c5_c_a;
  int32_T c5_d_a;
  int32_T c5_ihim1;
  int32_T c5_jcol;
  int32_T c5_e_a;
  int32_T c5_f_a;
  int32_T c5_jcolp1;
  int32_T c5_jrow;
  int32_T c5_g_a;
  int32_T c5_h_a;
  int32_T c5_jrowm1;
  creal_T c5_d_A;
  creal_T c5_e_A;
  creal_T c5_b;
  creal_T c5_s;
  real_T c5_b_c;
  real_T c5_c_c;
  real_T c5_d_c;
  int32_T c5_xrow;
  int32_T c5_yrow;
  int32_T c5_jlo;
  int32_T c5_jhi;
  int32_T c5_b_jlo;
  int32_T c5_b_jhi;
  int32_T c5_i_a;
  int32_T c5_b_b;
  int32_T c5_j_a;
  int32_T c5_c_b;
  boolean_T c5_overflow;
  int32_T c5_j;
  int32_T c5_b_j;
  real_T c5_k_a;
  creal_T c5_y;
  creal_T c5_b_s;
  creal_T c5_stemp;
  real_T c5_l_a;
  creal_T c5_d_b;
  creal_T c5_e_b;
  creal_T c5_f_b;
  creal_T c5_g_b;
  real_T c5_e_c;
  int32_T c5_xcol;
  int32_T c5_ycol;
  int32_T c5_d_ilo;
  int32_T c5_d_ihi;
  int32_T c5_e_ilo;
  int32_T c5_e_ihi;
  int32_T c5_m_a;
  int32_T c5_h_b;
  int32_T c5_n_a;
  int32_T c5_i_b;
  boolean_T c5_b_overflow;
  int32_T c5_i;
  int32_T c5_b_i;
  real_T c5_o_a;
  creal_T c5_c_s;
  real_T c5_p_a;
  creal_T c5_j_b;
  creal_T c5_k_b;
  creal_T c5_l_b;
  creal_T c5_m_b;
  int32_T c5_i47;
  creal_T c5_f_A[9];
  real_T c5_b_info;
  real_T c5_c_info;
  real_T c5_d_info;
  real_T c5_e_info;
  int32_T c5_i48;
  creal_T c5_b_alpha1[3];
  int32_T c5_i49;
  creal_T c5_b_beta1[3];
  boolean_T guard1 = false;
  for (c5_i43 = 0; c5_i43 < 9; c5_i43++) {
    c5_b_A[c5_i43].re = c5_A[c5_i43] + c5_dc0.re;
    c5_b_A[c5_i43].im = c5_dc0.im;
  }

  c5_info = 0.0;
  c5_realmin(chartInstance);
  c5_eps(chartInstance);
  for (c5_i44 = 0; c5_i44 < 9; c5_i44++) {
    c5_c_A[c5_i44] = c5_b_A[c5_i44];
  }

  c5_anrm = c5_eml_matlab_zlangeM(chartInstance, c5_c_A);
  if (!c5_isfinite(chartInstance, c5_anrm)) {
    for (c5_i45 = 0; c5_i45 < 3; c5_i45++) {
      c5_alpha1[c5_i45].re = rtNaN;
      c5_alpha1[c5_i45].im = 0.0;
    }

    for (c5_i46 = 0; c5_i46 < 3; c5_i46++) {
      c5_beta1[c5_i46].re = rtNaN;
      c5_beta1[c5_i46].im = 0.0;
    }
  } else {
    c5_ilascl = false;
    c5_anrmto = c5_anrm;
    guard1 = false;
    if (c5_anrm > 0.0) {
      if (c5_anrm < 6.7178761075670888E-139) {
        c5_anrmto = 6.7178761075670888E-139;
        c5_ilascl = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      if (c5_anrm > 1.4885657073574029E+138) {
        c5_anrmto = 1.4885657073574029E+138;
        c5_ilascl = true;
      }
    }

    if (c5_ilascl) {
      c5_c_eml_matlab_zlascl(chartInstance, c5_anrm, c5_anrmto, c5_b_A);
    }

    c5_b_eml_matlab_zggbal(chartInstance, c5_b_A, &c5_ilo, &c5_ihi, c5_rscale);
    c5_b_ilo = c5_ilo;
    c5_b_ihi = c5_ihi;
    c5_c_ilo = c5_b_ilo;
    c5_c_ihi = c5_b_ihi;
    c5_a = c5_c_ilo;
    c5_b_a = c5_a + 2;
    c5_c = c5_b_a;
    if (c5_c_ihi < c5_c) {
    } else {
      c5_c_a = c5_c_ihi;
      c5_d_a = c5_c_a - 1;
      c5_ihim1 = c5_d_a;
      c5_jcol = c5_c_ilo;
      while (c5_jcol < c5_ihim1) {
        c5_e_a = c5_jcol;
        c5_f_a = c5_e_a + 1;
        c5_jcolp1 = c5_f_a;
        c5_jrow = c5_c_ihi;
        while (c5_jrow > c5_jcolp1) {
          c5_g_a = c5_jrow;
          c5_h_a = c5_g_a - 1;
          c5_jrowm1 = c5_h_a;
          c5_d_A.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_jrowm1), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].re;
          c5_d_A.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_jrowm1), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].im;
          c5_e_A.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_jrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].re;
          c5_e_A.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_jrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].im;
          c5_eml_matlab_zlartg(chartInstance, c5_d_A, c5_e_A, &c5_b_c, &c5_s,
                               &c5_b);
          c5_c_c = c5_b_c;
          c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_jrowm1), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].re = c5_b.re;
          c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_jrowm1), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].im = c5_b.im;
          c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_jrow), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].re = c5_dc0.re;
          c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_jrow), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].im = c5_dc0.im;
          c5_d_c = c5_c_c;
          c5_xrow = c5_jrowm1;
          c5_yrow = c5_jrow;
          c5_jlo = c5_jcolp1;
          c5_jhi = c5_c_ihi;
          c5_b_jlo = c5_jlo;
          c5_b_jhi = c5_jhi;
          c5_i_a = c5_b_jlo;
          c5_b_b = c5_b_jhi;
          c5_j_a = c5_i_a;
          c5_c_b = c5_b_b;
          if (c5_j_a > c5_c_b) {
            c5_overflow = false;
          } else {
            c5_b_eml_switch_helper(chartInstance);
            c5_overflow = (c5_c_b > 2147483646);
          }

          if (c5_overflow) {
            c5_check_forloop_overflow_error(chartInstance, c5_overflow);
          }

          for (c5_j = c5_b_jlo; c5_j <= c5_b_jhi; c5_j++) {
            c5_b_j = c5_j;
            c5_k_a = c5_d_c;
            c5_b.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c5_b.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c5_y.re = c5_k_a * c5_b.re;
            c5_y.im = c5_k_a * c5_b.im;
            c5_b_s.re = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re - c5_s.im * c5_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c5_b_s.im = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im + c5_s.im * c5_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c5_stemp.re = c5_y.re + c5_b_s.re;
            c5_stemp.im = c5_y.im + c5_b_s.im;
            c5_l_a = c5_d_c;
            c5_b.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c5_b.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c5_y.re = c5_l_a * c5_b.re;
            c5_y.im = c5_l_a * c5_b.im;
            c5_b.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c5_b.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c5_d_b = c5_b;
            c5_e_b = c5_b;
            c5_f_b = c5_b;
            c5_g_b = c5_b;
            c5_b.re = c5_s.re * c5_d_b.re + c5_s.im * c5_e_b.im;
            c5_b.im = c5_s.re * c5_f_b.im - c5_s.im * c5_g_b.re;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re = c5_y.re
              - c5_b.re;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im = c5_y.im
              - c5_b.im;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re =
              c5_stemp.re;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im =
              c5_stemp.im;
          }

          c5_s.re = -c5_s.re;
          c5_s.im = -c5_s.im;
          c5_e_c = c5_c_c;
          c5_xcol = c5_jrow;
          c5_ycol = c5_jrowm1;
          c5_d_ilo = c5_c_ilo;
          c5_d_ihi = c5_c_ihi;
          c5_e_ilo = c5_d_ilo;
          c5_e_ihi = c5_d_ihi;
          c5_m_a = c5_e_ilo;
          c5_h_b = c5_e_ihi;
          c5_n_a = c5_m_a;
          c5_i_b = c5_h_b;
          if (c5_n_a > c5_i_b) {
            c5_b_overflow = false;
          } else {
            c5_b_eml_switch_helper(chartInstance);
            c5_b_overflow = (c5_i_b > 2147483646);
          }

          if (c5_b_overflow) {
            c5_check_forloop_overflow_error(chartInstance, c5_b_overflow);
          }

          for (c5_i = c5_e_ilo; c5_i <= c5_e_ihi; c5_i++) {
            c5_b_i = c5_i;
            c5_o_a = c5_e_c;
            c5_b.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].re;
            c5_b.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].im;
            c5_y.re = c5_o_a * c5_b.re;
            c5_y.im = c5_o_a * c5_b.im;
            c5_c_s.re = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re - c5_s.im * c5_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im;
            c5_c_s.im = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im + c5_s.im * c5_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re;
            c5_stemp.re = c5_y.re + c5_c_s.re;
            c5_stemp.im = c5_y.im + c5_c_s.im;
            c5_p_a = c5_e_c;
            c5_b.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re;
            c5_b.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im;
            c5_y.re = c5_p_a * c5_b.re;
            c5_y.im = c5_p_a * c5_b.im;
            c5_b.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].re;
            c5_b.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].im;
            c5_j_b = c5_b;
            c5_k_b = c5_b;
            c5_l_b = c5_b;
            c5_m_b = c5_b;
            c5_b.re = c5_s.re * c5_j_b.re + c5_s.im * c5_k_b.im;
            c5_b.im = c5_s.re * c5_l_b.im - c5_s.im * c5_m_b.re;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re = c5_y.re
              - c5_b.re;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im = c5_y.im
              - c5_b.im;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].re =
              c5_stemp.re;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].im =
              c5_stemp.im;
          }

          c5_jrow = c5_jrowm1;
        }

        c5_jcol = c5_jcolp1;
      }
    }

    for (c5_i47 = 0; c5_i47 < 9; c5_i47++) {
      c5_f_A[c5_i47] = c5_b_A[c5_i47];
    }

    c5_eml_matlab_zhgeqz(chartInstance, c5_f_A, c5_b_ilo, c5_b_ihi, &c5_b_info,
                         c5_alpha1, c5_beta1);
    c5_info = c5_b_info;
    if (c5_info != 0.0) {
    } else {
      if (c5_ilascl) {
        c5_d_eml_matlab_zlascl(chartInstance, c5_anrmto, c5_anrm, c5_alpha1);
      }
    }
  }

  c5_c_info = c5_info;
  c5_d_info = c5_c_info;
  c5_e_info = c5_d_info;
  for (c5_i48 = 0; c5_i48 < 3; c5_i48++) {
    c5_b_alpha1[c5_i48] = c5_alpha1[c5_i48];
  }

  for (c5_i49 = 0; c5_i49 < 3; c5_i49++) {
    c5_b_beta1[c5_i49] = c5_beta1[c5_i49];
  }

  c5_b_eml_div(chartInstance, c5_b_alpha1, c5_b_beta1, c5_V);
  if (c5_e_info < 0.0) {
    c5_eml_warning(chartInstance);
  } else {
    if (c5_e_info > 0.0) {
      c5_b_eml_warning(chartInstance);
    }
  }
}

static void c5_eml_error(SFc5_Model_01InstanceStruct *chartInstance)
{
  int32_T c5_i50;
  static char_T c5_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c5_u[30];
  const mxArray *c5_y = NULL;
  int32_T c5_i51;
  static char_T c5_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c5_b_u[4];
  const mxArray *c5_b_y = NULL;
  (void)chartInstance;
  for (c5_i50 = 0; c5_i50 < 30; c5_i50++) {
    c5_u[c5_i50] = c5_cv0[c5_i50];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  for (c5_i51 = 0; c5_i51 < 4; c5_i51++) {
    c5_b_u[c5_i51] = c5_cv1[c5_i51];
  }

  c5_b_y = NULL;
  sf_mex_assign(&c5_b_y, sf_mex_create("y", c5_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c5_y, 14, c5_b_y));
}

static void c5_eps(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c5_eml_matlab_zlangeM(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_x[9])
{
  real_T c5_y;
  int32_T c5_k;
  real_T c5_b_k;
  creal_T c5_b_x;
  real_T c5_x1;
  real_T c5_x2;
  real_T c5_a;
  real_T c5_b;
  real_T c5_absxk;
  real_T c5_c_x;
  boolean_T c5_b_b;
  boolean_T exitg1;
  (void)chartInstance;
  c5_y = 0.0;
  c5_k = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c5_k < 9)) {
    c5_b_k = 1.0 + (real_T)c5_k;
    c5_b_x.re = c5_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", c5_b_k), 1, 9, 1, 0) - 1].re;
    c5_b_x.im = c5_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", c5_b_k), 1, 9, 1, 0) - 1].im;
    c5_x1 = c5_b_x.re;
    c5_x2 = c5_b_x.im;
    c5_a = c5_x1;
    c5_b = c5_x2;
    c5_absxk = muDoubleScalarHypot(c5_a, c5_b);
    c5_c_x = c5_absxk;
    c5_b_b = muDoubleScalarIsNaN(c5_c_x);
    if (c5_b_b) {
      c5_y = rtNaN;
      exitg1 = true;
    } else {
      if (c5_absxk > c5_y) {
        c5_y = c5_absxk;
      }

      c5_k++;
    }
  }

  return c5_y;
}

static real_T c5_b_abs(SFc5_Model_01InstanceStruct *chartInstance, creal_T c5_x)
{
  real_T c5_x1;
  real_T c5_x2;
  real_T c5_a;
  real_T c5_b;
  (void)chartInstance;
  c5_x1 = c5_x.re;
  c5_x2 = c5_x.im;
  c5_a = c5_x1;
  c5_b = c5_x2;
  return muDoubleScalarHypot(c5_a, c5_b);
}

static boolean_T c5_isfinite(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_x)
{
  real_T c5_b_x;
  boolean_T c5_b_b;
  boolean_T c5_b0;
  real_T c5_c_x;
  boolean_T c5_c_b;
  boolean_T c5_b1;
  (void)chartInstance;
  c5_b_x = c5_x;
  c5_b_b = muDoubleScalarIsInf(c5_b_x);
  c5_b0 = !c5_b_b;
  c5_c_x = c5_x;
  c5_c_b = muDoubleScalarIsNaN(c5_c_x);
  c5_b1 = !c5_c_b;
  return c5_b0 && c5_b1;
}

static void c5_eml_matlab_zlascl(SFc5_Model_01InstanceStruct *chartInstance,
  real_T c5_cfrom, real_T c5_cto, creal_T c5_A[9], creal_T c5_b_A[9])
{
  int32_T c5_i52;
  for (c5_i52 = 0; c5_i52 < 9; c5_i52++) {
    c5_b_A[c5_i52] = c5_A[c5_i52];
  }

  c5_c_eml_matlab_zlascl(chartInstance, c5_cfrom, c5_cto, c5_b_A);
}

static void c5_eml_matlab_zggbal(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], creal_T c5_b_A[9], int32_T *c5_ilo, int32_T *c5_ihi, int32_T
  c5_rscale[3])
{
  int32_T c5_i53;
  for (c5_i53 = 0; c5_i53 < 9; c5_i53++) {
    c5_b_A[c5_i53] = c5_A[c5_i53];
  }

  c5_b_eml_matlab_zggbal(chartInstance, c5_b_A, c5_ilo, c5_ihi, c5_rscale);
}

static void c5_check_forloop_overflow_error(SFc5_Model_01InstanceStruct
  *chartInstance, boolean_T c5_overflow)
{
  int32_T c5_i54;
  static char_T c5_cv2[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c5_u[34];
  const mxArray *c5_y = NULL;
  int32_T c5_i55;
  static char_T c5_cv3[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c5_b_u[23];
  const mxArray *c5_b_y = NULL;
  (void)chartInstance;
  if (!c5_overflow) {
  } else {
    for (c5_i54 = 0; c5_i54 < 34; c5_i54++) {
      c5_u[c5_i54] = c5_cv2[c5_i54];
    }

    c5_y = NULL;
    sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  false);
    for (c5_i55 = 0; c5_i55 < 23; c5_i55++) {
      c5_b_u[c5_i55] = c5_cv3[c5_i55];
    }

    c5_b_y = NULL;
    sf_mex_assign(&c5_b_y, sf_mex_create("y", c5_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 2U, 14, c5_y, 14, c5_b_y));
  }
}

static void c5_eml_matlab_zlartg(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_f, creal_T c5_g, real_T *c5_cs, creal_T *c5_sn, creal_T *c5_r)
{
  real_T c5_x;
  real_T c5_b_x;
  real_T c5_y;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_b_y;
  real_T c5_e_x;
  real_T c5_c_y;
  real_T c5_d_y;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_e_y;
  real_T c5_h_x;
  real_T c5_i_x;
  real_T c5_f_y;
  real_T c5_j_x;
  real_T c5_g_y;
  real_T c5_h_y;
  real_T c5_k_x;
  real_T c5_i_y;
  real_T c5_scale;
  creal_T c5_fs;
  creal_T c5_gs;
  int32_T c5_count;
  real_T c5_rescaledir;
  int32_T c5_a;
  int32_T c5_b_a;
  static creal_T c5_dc1 = { 0.0, 0.0 };

  boolean_T c5_b_g;
  int32_T c5_c_a;
  int32_T c5_d_a;
  real_T c5_f2;
  real_T c5_g2;
  real_T c5_l_x;
  real_T c5_m_x;
  boolean_T c5_b_f;
  real_T c5_x1;
  real_T c5_x2;
  real_T c5_e_a;
  real_T c5_b;
  real_T c5_j_y;
  real_T c5_b_x1;
  real_T c5_b_x2;
  real_T c5_f_a;
  real_T c5_b_b;
  real_T c5_d;
  real_T c5_c_x1;
  real_T c5_c_x2;
  real_T c5_g_a;
  real_T c5_c_b;
  real_T c5_f2s;
  real_T c5_n_x;
  real_T c5_g2s;
  real_T c5_o_x;
  real_T c5_p_x;
  real_T c5_k_y;
  real_T c5_q_x;
  real_T c5_r_x;
  real_T c5_l_y;
  real_T c5_s_x;
  real_T c5_m_y;
  real_T c5_n_y;
  real_T c5_d_x1;
  real_T c5_d_x2;
  real_T c5_h_a;
  real_T c5_d_b;
  real_T c5_dr;
  real_T c5_di;
  real_T c5_e_x1;
  real_T c5_e_x2;
  real_T c5_i_a;
  real_T c5_e_b;
  creal_T c5_b_gs;
  real_T c5_j_a;
  creal_T c5_b_sn;
  real_T c5_t_x;
  creal_T c5_c_gs;
  creal_T c5_c_sn;
  int32_T c5_b_count;
  int32_T c5_f_b;
  int32_T c5_g_b;
  boolean_T c5_overflow;
  int32_T c5_i;
  int32_T c5_c_count;
  int32_T c5_h_b;
  int32_T c5_i_b;
  boolean_T c5_b_overflow;
  int32_T c5_b_i;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  c5_realmin(chartInstance);
  c5_eps(chartInstance);
  c5_eps(chartInstance);
  c5_x = c5_f.re;
  c5_b_x = c5_x;
  c5_y = muDoubleScalarAbs(c5_b_x);
  c5_c_x = c5_f.im;
  c5_d_x = c5_c_x;
  c5_b_y = muDoubleScalarAbs(c5_d_x);
  c5_e_x = c5_y;
  c5_c_y = c5_b_y;
  c5_d_y = c5_e_x;
  if (c5_c_y > c5_d_y) {
    c5_d_y = c5_c_y;
  }

  c5_f_x = c5_g.re;
  c5_g_x = c5_f_x;
  c5_e_y = muDoubleScalarAbs(c5_g_x);
  c5_h_x = c5_g.im;
  c5_i_x = c5_h_x;
  c5_f_y = muDoubleScalarAbs(c5_i_x);
  c5_j_x = c5_e_y;
  c5_g_y = c5_f_y;
  c5_h_y = c5_j_x;
  if (c5_g_y > c5_h_y) {
    c5_h_y = c5_g_y;
  }

  c5_k_x = c5_d_y;
  c5_i_y = c5_h_y;
  c5_scale = c5_k_x;
  if (c5_i_y > c5_scale) {
    c5_scale = c5_i_y;
  }

  c5_fs = c5_f;
  c5_gs = c5_g;
  c5_count = 0;
  c5_rescaledir = 0.0;
  guard1 = false;
  guard2 = false;
  if (c5_scale >= 7.4428285367870146E+137) {
    do {
      c5_a = c5_count;
      c5_b_a = c5_a + 1;
      c5_count = c5_b_a;
      c5_fs.re *= 1.3435752215134178E-138;
      c5_fs.im *= 1.3435752215134178E-138;
      c5_gs.re *= 1.3435752215134178E-138;
      c5_gs.im *= 1.3435752215134178E-138;
      c5_scale *= 1.3435752215134178E-138;
    } while (!(c5_scale < 7.4428285367870146E+137));

    c5_rescaledir = 1.0;
    guard1 = true;
  } else if (c5_scale <= 1.3435752215134178E-138) {
    c5_b_g = ((c5_g.re == c5_dc1.re) && (c5_g.im == c5_dc1.im));
    if (c5_b_g) {
      *c5_cs = 1.0;
      *c5_sn = c5_dc1;
      *c5_r = c5_f;
    } else {
      do {
        c5_c_a = c5_count;
        c5_d_a = c5_c_a + 1;
        c5_count = c5_d_a;
        c5_fs.re *= 7.4428285367870146E+137;
        c5_fs.im *= 7.4428285367870146E+137;
        c5_gs.re *= 7.4428285367870146E+137;
        c5_gs.im *= 7.4428285367870146E+137;
        c5_scale *= 7.4428285367870146E+137;
      } while (!(c5_scale > 1.3435752215134178E-138));

      c5_rescaledir = -1.0;
      guard2 = true;
    }
  } else {
    guard2 = true;
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c5_f2 = c5_fs.re * c5_fs.re + c5_fs.im * c5_fs.im;
    c5_g2 = c5_gs.re * c5_gs.re + c5_gs.im * c5_gs.im;
    c5_l_x = c5_g2;
    c5_m_x = c5_l_x;
    if (1.0 > c5_m_x) {
      c5_m_x = 1.0;
    }

    if (c5_f2 <= c5_m_x * 2.0041683600089728E-292) {
      c5_b_f = ((c5_f.re == c5_dc1.re) && (c5_f.im == c5_dc1.im));
      if (c5_b_f) {
        *c5_cs = 0.0;
        c5_x1 = c5_g.re;
        c5_x2 = c5_g.im;
        c5_e_a = c5_x1;
        c5_b = c5_x2;
        c5_j_y = muDoubleScalarHypot(c5_e_a, c5_b);
        c5_r->re = c5_j_y;
        c5_r->im = 0.0;
        c5_b_x1 = c5_gs.re;
        c5_b_x2 = c5_gs.im;
        c5_f_a = c5_b_x1;
        c5_b_b = c5_b_x2;
        c5_d = muDoubleScalarHypot(c5_f_a, c5_b_b);
        c5_sn->re = c5_gs.re / c5_d;
        c5_sn->im = -c5_gs.im / c5_d;
      } else {
        c5_c_x1 = c5_fs.re;
        c5_c_x2 = c5_fs.im;
        c5_g_a = c5_c_x1;
        c5_c_b = c5_c_x2;
        c5_f2s = muDoubleScalarHypot(c5_g_a, c5_c_b);
        c5_n_x = c5_g2;
        c5_g2s = c5_n_x;
        if (c5_g2s < 0.0) {
          c5_eml_error(chartInstance);
        }

        c5_g2s = muDoubleScalarSqrt(c5_g2s);
        *c5_cs = c5_f2s / c5_g2s;
        c5_o_x = c5_f.re;
        c5_p_x = c5_o_x;
        c5_k_y = muDoubleScalarAbs(c5_p_x);
        c5_q_x = c5_f.im;
        c5_r_x = c5_q_x;
        c5_l_y = muDoubleScalarAbs(c5_r_x);
        c5_s_x = c5_k_y;
        c5_m_y = c5_l_y;
        c5_n_y = c5_s_x;
        if (c5_m_y > c5_n_y) {
          c5_n_y = c5_m_y;
        }

        if (c5_n_y > 1.0) {
          c5_d_x1 = c5_f.re;
          c5_d_x2 = c5_f.im;
          c5_h_a = c5_d_x1;
          c5_d_b = c5_d_x2;
          c5_d = muDoubleScalarHypot(c5_h_a, c5_d_b);
          c5_fs.re = c5_f.re / c5_d;
          c5_fs.im = c5_f.im / c5_d;
        } else {
          c5_dr = 7.4428285367870146E+137 * c5_f.re;
          c5_di = 7.4428285367870146E+137 * c5_f.im;
          c5_e_x1 = c5_dr;
          c5_e_x2 = c5_di;
          c5_i_a = c5_e_x1;
          c5_e_b = c5_e_x2;
          c5_d = muDoubleScalarHypot(c5_i_a, c5_e_b);
          c5_fs.re = c5_dr / c5_d;
          c5_fs.im = c5_di / c5_d;
        }

        c5_b_gs.re = c5_gs.re / c5_g2s;
        c5_b_gs.im = -c5_gs.im / c5_g2s;
        c5_sn->re = c5_fs.re * c5_b_gs.re - c5_fs.im * c5_b_gs.im;
        c5_sn->im = c5_fs.re * c5_b_gs.im + c5_fs.im * c5_b_gs.re;
        c5_j_a = *c5_cs;
        c5_fs.re = c5_j_a * c5_f.re;
        c5_fs.im = c5_j_a * c5_f.im;
        c5_b_sn.re = c5_sn->re * c5_g.re - c5_sn->im * c5_g.im;
        c5_b_sn.im = c5_sn->re * c5_g.im + c5_sn->im * c5_g.re;
        c5_r->re = c5_fs.re + c5_b_sn.re;
        c5_r->im = c5_fs.im + c5_b_sn.im;
      }
    } else {
      c5_t_x = 1.0 + c5_g2 / c5_f2;
      c5_f2s = c5_t_x;
      if (c5_f2s < 0.0) {
        c5_eml_error(chartInstance);
      }

      c5_f2s = muDoubleScalarSqrt(c5_f2s);
      c5_r->re = c5_f2s * c5_fs.re;
      c5_r->im = c5_f2s * c5_fs.im;
      *c5_cs = 1.0 / c5_f2s;
      c5_d = c5_f2 + c5_g2;
      c5_sn->re = c5_r->re / c5_d;
      c5_sn->im = c5_r->im / c5_d;
      c5_c_gs.re = c5_gs.re;
      c5_c_gs.im = -c5_gs.im;
      c5_c_sn = *c5_sn;
      c5_sn->re = c5_c_sn.re * c5_c_gs.re - c5_c_sn.im * c5_c_gs.im;
      c5_sn->im = c5_c_sn.re * c5_c_gs.im + c5_c_sn.im * c5_c_gs.re;
      if (c5_rescaledir > 0.0) {
        c5_b_count = c5_count;
        c5_f_b = c5_b_count;
        c5_g_b = c5_f_b;
        if (1 > c5_g_b) {
          c5_overflow = false;
        } else {
          c5_b_eml_switch_helper(chartInstance);
          c5_overflow = (c5_g_b > 2147483646);
        }

        if (c5_overflow) {
          c5_check_forloop_overflow_error(chartInstance, c5_overflow);
        }

        for (c5_i = 1; c5_i <= c5_b_count; c5_i++) {
          c5_r->re *= 7.4428285367870146E+137;
          c5_r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (c5_rescaledir < 0.0) {
          c5_c_count = c5_count;
          c5_h_b = c5_c_count;
          c5_i_b = c5_h_b;
          if (1 > c5_i_b) {
            c5_b_overflow = false;
          } else {
            c5_b_eml_switch_helper(chartInstance);
            c5_b_overflow = (c5_i_b > 2147483646);
          }

          if (c5_b_overflow) {
            c5_check_forloop_overflow_error(chartInstance, c5_b_overflow);
          }

          for (c5_b_i = 1; c5_b_i <= c5_c_count; c5_b_i++) {
            c5_r->re *= 1.3435752215134178E-138;
            c5_r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

static void c5_eml_scalar_eg(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c5_eml_matlab_zhgeqz(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T c5_ilo, int32_T c5_ihi, real_T *c5_info, creal_T
  c5_alpha1[3], creal_T c5_beta1[3])
{
  static creal_T c5_dc2 = { 0.0, 0.0 };

  int32_T c5_i56;
  creal_T c5_b_A[9];
  int32_T c5_i57;
  int32_T c5_i58;
  static creal_T c5_dc3 = { 0.0, 0.0 };

  creal_T c5_eshift;
  creal_T c5_ctemp;
  creal_T c5_rho;
  int32_T c5_i59;
  creal_T c5_c_A[9];
  real_T c5_anorm;
  real_T c5_y;
  real_T c5_atol;
  real_T c5_b_y;
  real_T c5_x;
  real_T c5_ascale;
  boolean_T c5_failed;
  int32_T c5_a;
  int32_T c5_b_a;
  int32_T c5_i60;
  int32_T c5_c_a;
  int32_T c5_d_a;
  boolean_T c5_overflow;
  int32_T c5_j;
  int32_T c5_b_j;
  int32_T c5_ifirst;
  int32_T c5_istart;
  int32_T c5_ilast;
  int32_T c5_e_a;
  int32_T c5_f_a;
  int32_T c5_ilastm1;
  int32_T c5_ifrstm;
  int32_T c5_ilastm;
  int32_T c5_iiter;
  int32_T c5_g_a;
  int32_T c5_b;
  int32_T c5_h_a;
  int32_T c5_b_b;
  int32_T c5_c;
  int32_T c5_i_a;
  int32_T c5_j_a;
  int32_T c5_b_c;
  int32_T c5_c_b;
  int32_T c5_d_b;
  int32_T c5_maxit;
  boolean_T c5_goto50;
  boolean_T c5_goto60;
  boolean_T c5_goto70;
  boolean_T c5_goto90;
  int32_T c5_b_maxit;
  int32_T c5_e_b;
  int32_T c5_f_b;
  boolean_T c5_b_overflow;
  int32_T c5_jiter;
  creal_T c5_a22;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_c_y;
  real_T c5_d_x;
  real_T c5_e_x;
  real_T c5_d_y;
  real_T c5_e_y;
  int32_T c5_k_a;
  int32_T c5_l_a;
  int32_T c5_jm1;
  boolean_T c5_ilazro;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_f_y;
  real_T c5_h_x;
  real_T c5_i_x;
  real_T c5_g_y;
  real_T c5_h_y;
  boolean_T c5_b2;
  int32_T c5_i61;
  int32_T c5_i62;
  creal_T c5_d_A;
  creal_T c5_e_A;
  creal_T c5_s;
  real_T c5_c_c;
  real_T c5_d_c;
  real_T c5_e_c;
  int32_T c5_xcol;
  int32_T c5_ycol;
  int32_T c5_b_ilo;
  int32_T c5_b_ihi;
  int32_T c5_c_ilo;
  int32_T c5_c_ihi;
  int32_T c5_m_a;
  int32_T c5_g_b;
  int32_T c5_n_a;
  int32_T c5_h_b;
  boolean_T c5_c_overflow;
  int32_T c5_i;
  int32_T c5_b_i;
  real_T c5_o_a;
  creal_T c5_a12;
  creal_T c5_b_s;
  creal_T c5_a21;
  real_T c5_p_a;
  creal_T c5_b_a22;
  creal_T c5_c_a22;
  creal_T c5_d_a22;
  creal_T c5_e_a22;
  int32_T c5_q_a;
  int32_T c5_r_a;
  int32_T c5_s_a;
  int32_T c5_t_a;
  creal_T c5_r2;
  creal_T c5_f_a22;
  creal_T c5_b_rho;
  creal_T c5_b_a12;
  creal_T c5_c_a12;
  creal_T c5_b_a21;
  real_T c5_d4;
  real_T c5_d5;
  int32_T c5_u_a;
  int32_T c5_v_a;
  int32_T c5_jp1;
  int32_T c5_w_a;
  int32_T c5_x_a;
  real_T c5_j_x;
  real_T c5_k_x;
  real_T c5_i_y;
  real_T c5_l_x;
  real_T c5_m_x;
  real_T c5_j_y;
  real_T c5_k_y;
  real_T c5_temp;
  real_T c5_n_x;
  real_T c5_o_x;
  real_T c5_l_y;
  real_T c5_p_x;
  real_T c5_q_x;
  real_T c5_m_y;
  real_T c5_n_y;
  real_T c5_temp2;
  real_T c5_r_x;
  real_T c5_o_y;
  real_T c5_tempr;
  real_T c5_s_x;
  real_T c5_t_x;
  real_T c5_p_y;
  real_T c5_u_x;
  real_T c5_v_x;
  real_T c5_q_y;
  real_T c5_r_y;
  int32_T c5_y_a;
  int32_T c5_ab_a;
  int32_T c5_f_c;
  real_T c5_g_c;
  int32_T c5_bb_a;
  int32_T c5_cb_a;
  int32_T c5_db_a;
  int32_T c5_eb_a;
  creal_T c5_f_A;
  creal_T c5_g_A;
  real_T c5_h_c;
  real_T c5_i_c;
  int32_T c5_xrow;
  int32_T c5_yrow;
  int32_T c5_jlo;
  int32_T c5_jhi;
  int32_T c5_b_jlo;
  int32_T c5_b_jhi;
  int32_T c5_fb_a;
  int32_T c5_i_b;
  int32_T c5_gb_a;
  int32_T c5_j_b;
  boolean_T c5_d_overflow;
  int32_T c5_c_j;
  int32_T c5_d_j;
  real_T c5_hb_a;
  creal_T c5_c_s;
  real_T c5_ib_a;
  creal_T c5_g_a22;
  creal_T c5_h_a22;
  creal_T c5_i_a22;
  creal_T c5_j_a22;
  int32_T c5_jb_a;
  int32_T c5_kb_a;
  int32_T c5_j_c;
  int32_T c5_w_x;
  int32_T c5_s_y;
  int32_T c5_x_x;
  real_T c5_k_c;
  int32_T c5_b_xcol;
  int32_T c5_b_ycol;
  int32_T c5_d_ilo;
  int32_T c5_d_ihi;
  int32_T c5_e_ilo;
  int32_T c5_e_ihi;
  int32_T c5_lb_a;
  int32_T c5_k_b;
  int32_T c5_mb_a;
  int32_T c5_l_b;
  boolean_T c5_e_overflow;
  int32_T c5_c_i;
  int32_T c5_d_i;
  real_T c5_nb_a;
  creal_T c5_d_s;
  real_T c5_ob_a;
  creal_T c5_k_a22;
  creal_T c5_l_a22;
  creal_T c5_m_a22;
  creal_T c5_n_a22;
  int32_T c5_b_ilast;
  int32_T c5_m_b;
  int32_T c5_n_b;
  boolean_T c5_f_overflow;
  int32_T c5_k;
  int32_T c5_b_k;
  int32_T c5_pb_a;
  int32_T c5_qb_a;
  int32_T c5_i63;
  int32_T c5_o_b;
  int32_T c5_p_b;
  boolean_T c5_g_overflow;
  int32_T c5_e_j;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  int32_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T guard11 = false;
  c5_dc2.re = rtNaN;
  for (c5_i56 = 0; c5_i56 < 9; c5_i56++) {
    c5_b_A[c5_i56] = c5_A[c5_i56];
  }

  for (c5_i57 = 0; c5_i57 < 3; c5_i57++) {
    c5_alpha1[c5_i57].re = 0.0;
    c5_alpha1[c5_i57].im = 0.0;
  }

  for (c5_i58 = 0; c5_i58 < 3; c5_i58++) {
    c5_beta1[c5_i58].re = 1.0;
    c5_beta1[c5_i58].im = 0.0;
  }

  c5_eps(chartInstance);
  c5_realmin(chartInstance);
  c5_eshift = c5_dc3;
  c5_ctemp = c5_dc3;
  c5_rho = c5_dc3;
  for (c5_i59 = 0; c5_i59 < 9; c5_i59++) {
    c5_c_A[c5_i59] = c5_b_A[c5_i59];
  }

  c5_anorm = c5_eml_matlab_zlanhs(chartInstance, c5_c_A, c5_ilo, c5_ihi);
  c5_y = 2.2204460492503131E-16 * c5_anorm;
  c5_atol = 2.2250738585072014E-308;
  if (c5_y > 2.2250738585072014E-308) {
    c5_atol = c5_y;
  }

  c5_b_y = c5_anorm;
  c5_x = 2.2250738585072014E-308;
  if (c5_b_y > 2.2250738585072014E-308) {
    c5_x = c5_b_y;
  }

  c5_ascale = 1.0 / c5_x;
  c5_failed = true;
  c5_a = c5_ihi;
  c5_b_a = c5_a + 1;
  c5_i60 = c5_b_a;
  c5_c_a = c5_i60;
  c5_d_a = c5_c_a;
  if (c5_d_a > 3) {
    c5_overflow = false;
  } else {
    c5_b_eml_switch_helper(chartInstance);
    c5_overflow = false;
  }

  if (c5_overflow) {
    c5_check_forloop_overflow_error(chartInstance, c5_overflow);
  }

  for (c5_j = c5_i60; c5_j < 4; c5_j++) {
    c5_b_j = c5_j;
    c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_b_j), 1, 3, 1, 0) - 1].re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
      ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
    c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_b_j), 1, 3, 1, 0) - 1].im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
      ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
  }

  guard1 = false;
  guard2 = false;
  if (c5_ihi >= c5_ilo) {
    c5_ifirst = c5_ilo;
    c5_istart = c5_ilo;
    c5_ilast = c5_ihi;
    c5_e_a = c5_ilast;
    c5_f_a = c5_e_a - 1;
    c5_ilastm1 = c5_f_a;
    c5_ifrstm = c5_ilo;
    c5_ilastm = c5_ihi;
    c5_iiter = 0;
    c5_g_a = c5_ihi;
    c5_b = c5_ilo;
    c5_h_a = c5_g_a;
    c5_b_b = c5_b;
    c5_c = c5_h_a - c5_b_b;
    c5_i_a = c5_c;
    c5_j_a = c5_i_a;
    c5_b_c = c5_j_a;
    c5_c_b = c5_b_c + 1;
    c5_d_b = c5_c_b;
    c5_maxit = 30 * c5_d_b;
    c5_goto50 = false;
    c5_goto60 = false;
    c5_goto70 = false;
    c5_goto90 = false;
    c5_b_maxit = c5_maxit;
    c5_e_b = c5_b_maxit;
    c5_f_b = c5_e_b;
    if (1 > c5_f_b) {
      c5_b_overflow = false;
    } else {
      c5_b_eml_switch_helper(chartInstance);
      c5_b_overflow = (c5_f_b > 2147483646);
    }

    if (c5_b_overflow) {
      c5_check_forloop_overflow_error(chartInstance, c5_b_overflow);
    }

    c5_jiter = 1;
    do {
      exitg1 = 0;
      if (c5_jiter <= c5_b_maxit) {
        if (c5_ilast == c5_ilo) {
          c5_goto60 = true;
        } else {
          c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            re;
          c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            im;
          c5_b_x = c5_a22.re;
          c5_c_x = c5_b_x;
          c5_c_y = muDoubleScalarAbs(c5_c_x);
          c5_d_x = c5_a22.im;
          c5_e_x = c5_d_x;
          c5_d_y = muDoubleScalarAbs(c5_e_x);
          c5_e_y = c5_c_y + c5_d_y;
          if (c5_e_y <= c5_atol) {
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c5_dc3.re;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c5_dc3.im;
            c5_goto60 = true;
          } else {
            c5_b_j = c5_ilastm1;
            exitg3 = false;
            while ((exitg3 == false) && (c5_b_j >= c5_ilo)) {
              c5_k_a = c5_b_j;
              c5_l_a = c5_k_a - 1;
              c5_jm1 = c5_l_a;
              if (c5_b_j == c5_ilo) {
                c5_ilazro = true;
              } else {
                c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c5_f_x = c5_a22.re;
                c5_g_x = c5_f_x;
                c5_f_y = muDoubleScalarAbs(c5_g_x);
                c5_h_x = c5_a22.im;
                c5_i_x = c5_h_x;
                c5_g_y = muDoubleScalarAbs(c5_i_x);
                c5_h_y = c5_f_y + c5_g_y;
                if (c5_h_y <= c5_atol) {
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0)
                               - 1)) - 1].re = c5_dc3.re;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0)
                               - 1)) - 1].im = c5_dc3.im;
                  c5_ilazro = true;
                } else {
                  c5_ilazro = false;
                }
              }

              if (c5_ilazro) {
                c5_ifirst = c5_b_j;
                c5_goto70 = true;
                exitg3 = true;
              } else {
                c5_b_j = c5_jm1;
              }
            }
          }
        }

        guard3 = false;
        guard4 = false;
        if (c5_goto50) {
          guard4 = true;
        } else if (c5_goto60) {
          guard4 = true;
        } else if (c5_goto70) {
          guard3 = true;
        } else {
          c5_b2 = false;
        }

        if (guard4 == true) {
          guard3 = true;
        }

        if (guard3 == true) {
          c5_b2 = true;
        }

        if (!c5_b2) {
          for (c5_i61 = 0; c5_i61 < 3; c5_i61++) {
            c5_alpha1[c5_i61] = c5_dc2;
          }

          for (c5_i62 = 0; c5_i62 < 3; c5_i62++) {
            c5_beta1[c5_i62] = c5_dc2;
          }

          *c5_info = -1.0;
          exitg1 = 1;
        } else {
          if (c5_goto50) {
            c5_goto50 = false;
            c5_d_A.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].
              re;
            c5_d_A.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].
              im;
            c5_e_A.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .re;
            c5_e_A.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .im;
            c5_eml_matlab_zlartg(chartInstance, c5_d_A, c5_e_A, &c5_c_c, &c5_s,
                                 &c5_a22);
            c5_d_c = c5_c_c;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].re =
              c5_a22.re;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].im =
              c5_a22.im;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c5_dc3.re;
            c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c5_dc3.im;
            c5_e_c = c5_d_c;
            c5_xcol = c5_ilast;
            c5_ycol = c5_ilastm1;
            c5_b_ilo = c5_ifrstm;
            c5_b_ihi = c5_ilastm1;
            c5_c_ilo = c5_b_ilo;
            c5_c_ihi = c5_b_ihi;
            c5_m_a = c5_c_ilo;
            c5_g_b = c5_c_ihi;
            c5_n_a = c5_m_a;
            c5_h_b = c5_g_b;
            if (c5_n_a > c5_h_b) {
              c5_c_overflow = false;
            } else {
              c5_b_eml_switch_helper(chartInstance);
              c5_c_overflow = (c5_h_b > 2147483646);
            }

            if (c5_c_overflow) {
              c5_check_forloop_overflow_error(chartInstance, c5_c_overflow);
            }

            for (c5_i = c5_c_ilo; c5_i <= c5_c_ihi; c5_i++) {
              c5_b_i = c5_i;
              c5_o_a = c5_e_c;
              c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c5_a12.re = c5_o_a * c5_a22.re;
              c5_a12.im = c5_o_a * c5_a22.im;
              c5_b_s.re = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re - c5_s.im *
                c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) -
                           1)) - 1].im;
              c5_b_s.im = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im + c5_s.im *
                c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) -
                           1)) - 1].re;
              c5_a21.re = c5_a12.re + c5_b_s.re;
              c5_a21.im = c5_a12.im + c5_b_s.im;
              c5_p_a = c5_e_c;
              c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c5_a12.re = c5_p_a * c5_a22.re;
              c5_a12.im = c5_p_a * c5_a22.im;
              c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c5_b_a22 = c5_a22;
              c5_c_a22 = c5_a22;
              c5_d_a22 = c5_a22;
              c5_e_a22 = c5_a22;
              c5_a22.re = c5_s.re * c5_b_a22.re + c5_s.im * c5_c_a22.im;
              c5_a22.im = c5_s.re * c5_d_a22.im - c5_s.im * c5_e_a22.re;
              c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1))
                - 1].re = c5_a12.re - c5_a22.re;
              c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1))
                - 1].im = c5_a12.im - c5_a22.im;
              c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1))
                - 1].re = c5_a21.re;
              c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1))
                - 1].im = c5_a21.im;
            }

            c5_goto60 = true;
          }

          guard11 = false;
          if (c5_goto60) {
            c5_goto60 = false;
            c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) - 1].re =
              c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3
                      * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) -
                         1)) - 1].re;
            c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) - 1].im =
              c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3
                      * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) -
                         1)) - 1].im;
            c5_ilast = c5_ilastm1;
            c5_q_a = c5_ilast;
            c5_r_a = c5_q_a - 1;
            c5_ilastm1 = c5_r_a;
            if (c5_ilast < c5_ilo) {
              c5_failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              c5_iiter = 0;
              c5_eshift = c5_dc3;
              c5_ilastm = c5_ilast;
              if (c5_ifrstm > c5_ilast) {
                c5_ifrstm = c5_ilo;
              }

              guard11 = true;
            }
          } else {
            if (c5_goto70) {
              c5_goto70 = false;
              c5_s_a = c5_iiter;
              c5_t_a = c5_s_a + 1;
              c5_iiter = c5_t_a;
              c5_ifrstm = c5_ifirst;
              if (c5_mod(chartInstance, c5_iiter) != 0) {
                c5_s.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c5_s.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c5_r2.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                   (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c5_r2.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                   (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c5_a12.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c5_a12.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c5_a21.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c5_a21.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c5_a22.re = c5_r2.re - c5_s.re;
                c5_a22.im = c5_r2.im - c5_s.im;
                c5_f_a22.re = -c5_a22.re;
                c5_f_a22.im = -c5_a22.im;
                c5_rho = c5_eml_div(chartInstance, c5_f_a22, 2.0);
                c5_b_rho.re = c5_rho.re * c5_rho.re - c5_rho.im * c5_rho.im;
                c5_b_rho.im = c5_rho.re * c5_rho.im + c5_rho.im * c5_rho.re;
                c5_b_a12.re = c5_a12.re * c5_a21.re - c5_a12.im * c5_a21.im;
                c5_b_a12.im = c5_a12.re * c5_a21.im + c5_a12.im * c5_a21.re;
                c5_a22.re = c5_b_rho.re + c5_b_a12.re;
                c5_a22.im = c5_b_rho.im + c5_b_a12.im;
                c5_b_sqrt(chartInstance, &c5_a22);
                c5_a12.re = c5_s.re - (c5_rho.re - c5_a22.re);
                c5_a12.im = c5_s.im - (c5_rho.im - c5_a22.im);
                c5_a21.re = c5_s.re - (c5_rho.re + c5_a22.re);
                c5_a21.im = c5_s.im - (c5_rho.im + c5_a22.im);
                c5_c_a12.re = c5_a12.re - c5_r2.re;
                c5_c_a12.im = c5_a12.im - c5_r2.im;
                c5_b_a21.re = c5_a21.re - c5_r2.re;
                c5_b_a21.im = c5_a21.im - c5_r2.im;
                c5_d4 = c5_b_abs(chartInstance, c5_c_a12);
                c5_d5 = c5_b_abs(chartInstance, c5_b_a21);
                if (c5_d4 <= c5_d5) {
                  c5_a21 = c5_a12;
                  c5_rho.re -= c5_a22.re;
                  c5_rho.im -= c5_a22.im;
                } else {
                  c5_rho.re += c5_a22.re;
                  c5_rho.im += c5_a22.im;
                }
              } else {
                c5_eshift.re += c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                  "", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].re;
                c5_eshift.im += c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                  "", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].im;
                c5_a21 = c5_eshift;
              }

              c5_b_j = c5_ilastm1;
              c5_u_a = c5_b_j;
              c5_v_a = c5_u_a + 1;
              c5_jp1 = c5_v_a;
              exitg2 = false;
              while ((exitg2 == false) && (c5_b_j > c5_ifirst)) {
                c5_w_a = c5_b_j;
                c5_x_a = c5_w_a - 1;
                c5_jm1 = c5_x_a;
                c5_istart = c5_b_j;
                c5_ctemp.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re - c5_a21.re;
                c5_ctemp.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im - c5_a21.im;
                c5_j_x = c5_ctemp.re;
                c5_k_x = c5_j_x;
                c5_i_y = muDoubleScalarAbs(c5_k_x);
                c5_l_x = c5_ctemp.im;
                c5_m_x = c5_l_x;
                c5_j_y = muDoubleScalarAbs(c5_m_x);
                c5_k_y = c5_i_y + c5_j_y;
                c5_temp = c5_ascale * c5_k_y;
                c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c5_n_x = c5_a22.re;
                c5_o_x = c5_n_x;
                c5_l_y = muDoubleScalarAbs(c5_o_x);
                c5_p_x = c5_a22.im;
                c5_q_x = c5_p_x;
                c5_m_y = muDoubleScalarAbs(c5_q_x);
                c5_n_y = c5_l_y + c5_m_y;
                c5_temp2 = c5_ascale * c5_n_y;
                c5_r_x = c5_temp;
                c5_o_y = c5_temp2;
                c5_tempr = c5_r_x;
                if (c5_o_y > c5_tempr) {
                  c5_tempr = c5_o_y;
                }

                if (c5_tempr < 1.0) {
                  if (c5_tempr != 0.0) {
                    c5_temp /= c5_tempr;
                    c5_temp2 /= c5_tempr;
                  }
                }

                c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c5_s_x = c5_a22.re;
                c5_t_x = c5_s_x;
                c5_p_y = muDoubleScalarAbs(c5_t_x);
                c5_u_x = c5_a22.im;
                c5_v_x = c5_u_x;
                c5_q_y = muDoubleScalarAbs(c5_v_x);
                c5_r_y = c5_p_y + c5_q_y;
                if (c5_r_y * c5_temp2 <= c5_temp * c5_atol) {
                  c5_goto90 = true;
                  exitg2 = true;
                } else {
                  c5_jp1 = c5_b_j;
                  c5_b_j = c5_jm1;
                }
              }

              if (!c5_goto90) {
                c5_istart = c5_ifirst;
                if (c5_istart == c5_ilastm1) {
                  c5_ctemp = c5_rho;
                } else {
                  c5_ctemp.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 1, 0) + 3 *
                                        (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 2,
                    0) - 1)) - 1].re - c5_a21.re;
                  c5_ctemp.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 1, 0) + 3 *
                                        (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 2,
                    0) - 1)) - 1].im - c5_a21.im;
                }

                c5_goto90 = true;
              }
            }

            if (c5_goto90) {
              c5_goto90 = false;
              c5_y_a = c5_istart;
              c5_ab_a = c5_y_a;
              c5_f_c = c5_ab_a;
              c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c5_f_c + 1)), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c5_f_c + 1)), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c5_b_eml_matlab_zlartg(chartInstance, c5_ctemp, c5_a22, &c5_g_c,
                &c5_s);
              c5_d_c = c5_g_c;
              c5_b_j = c5_istart;
              c5_bb_a = c5_b_j;
              c5_cb_a = c5_bb_a - 1;
              c5_jm1 = c5_cb_a;
              while (c5_b_j < c5_ilast) {
                c5_db_a = c5_b_j;
                c5_eb_a = c5_db_a + 1;
                c5_jp1 = c5_eb_a;
                if (c5_b_j > c5_istart) {
                  c5_f_A.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_f_A.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_g_A.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_g_A.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_eml_matlab_zlartg(chartInstance, c5_f_A, c5_g_A, &c5_h_c,
                                       &c5_s, &c5_a22);
                  c5_d_c = c5_h_c;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0)
                               - 1)) - 1].re = c5_a22.re;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0)
                               - 1)) - 1].im = c5_a22.im;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0)
                               - 1)) - 1].re = c5_dc3.re;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0)
                               - 1)) - 1].im = c5_dc3.im;
                }

                c5_i_c = c5_d_c;
                c5_xrow = c5_b_j;
                c5_yrow = c5_jp1;
                c5_jlo = c5_b_j;
                c5_jhi = c5_ilastm;
                c5_b_jlo = c5_jlo;
                c5_b_jhi = c5_jhi;
                c5_fb_a = c5_b_jlo;
                c5_i_b = c5_b_jhi;
                c5_gb_a = c5_fb_a;
                c5_j_b = c5_i_b;
                if (c5_gb_a > c5_j_b) {
                  c5_d_overflow = false;
                } else {
                  c5_b_eml_switch_helper(chartInstance);
                  c5_d_overflow = (c5_j_b > 2147483646);
                }

                if (c5_d_overflow) {
                  c5_check_forloop_overflow_error(chartInstance, c5_d_overflow);
                }

                for (c5_c_j = c5_b_jlo; c5_c_j <= c5_b_jhi; c5_c_j++) {
                  c5_d_j = c5_c_j;
                  c5_hb_a = c5_i_c;
                  c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_a12.re = c5_hb_a * c5_a22.re;
                  c5_a12.im = c5_hb_a * c5_a22.im;
                  c5_c_s.re = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re - c5_s.im * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_c_s.im = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im + c5_s.im * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_a21.re = c5_a12.re + c5_c_s.re;
                  c5_a21.im = c5_a12.im + c5_c_s.im;
                  c5_ib_a = c5_i_c;
                  c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_a12.re = c5_ib_a * c5_a22.re;
                  c5_a12.im = c5_ib_a * c5_a22.im;
                  c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_g_a22 = c5_a22;
                  c5_h_a22 = c5_a22;
                  c5_i_a22 = c5_a22;
                  c5_j_a22 = c5_a22;
                  c5_a22.re = c5_s.re * c5_g_a22.re + c5_s.im * c5_h_a22.im;
                  c5_a22.im = c5_s.re * c5_i_a22.im - c5_s.im * c5_j_a22.re;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].re = c5_a12.re - c5_a22.re;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].im = c5_a12.im - c5_a22.im;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].re = c5_a21.re;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].im = c5_a21.im;
                }

                c5_s.re = -c5_s.re;
                c5_s.im = -c5_s.im;
                c5_jb_a = c5_jp1;
                c5_kb_a = c5_jb_a;
                c5_j_c = c5_kb_a;
                c5_w_x = c5_j_c + 1;
                c5_s_y = c5_ilast;
                c5_x_x = c5_w_x;
                if (c5_s_y < c5_x_x) {
                  c5_x_x = c5_s_y;
                }

                c5_k_c = c5_d_c;
                c5_b_xcol = c5_jp1;
                c5_b_ycol = c5_b_j;
                c5_d_ilo = c5_ifrstm;
                c5_d_ihi = c5_x_x;
                c5_e_ilo = c5_d_ilo;
                c5_e_ihi = c5_d_ihi;
                c5_lb_a = c5_e_ilo;
                c5_k_b = c5_e_ihi;
                c5_mb_a = c5_lb_a;
                c5_l_b = c5_k_b;
                if (c5_mb_a > c5_l_b) {
                  c5_e_overflow = false;
                } else {
                  c5_b_eml_switch_helper(chartInstance);
                  c5_e_overflow = (c5_l_b > 2147483646);
                }

                if (c5_e_overflow) {
                  c5_check_forloop_overflow_error(chartInstance, c5_e_overflow);
                }

                for (c5_c_i = c5_e_ilo; c5_c_i <= c5_e_ihi; c5_c_i++) {
                  c5_d_i = c5_c_i;
                  c5_nb_a = c5_k_c;
                  c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_a12.re = c5_nb_a * c5_a22.re;
                  c5_a12.im = c5_nb_a * c5_a22.im;
                  c5_d_s.re = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].re - c5_s.im * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_d_s.im = c5_s.re * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].im + c5_s.im * c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a21.re = c5_a12.re + c5_d_s.re;
                  c5_a21.im = c5_a12.im + c5_d_s.im;
                  c5_ob_a = c5_k_c;
                  c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_a12.re = c5_ob_a * c5_a22.re;
                  c5_a12.im = c5_ob_a * c5_a22.im;
                  c5_a22.re = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a22.im = c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_k_a22 = c5_a22;
                  c5_l_a22 = c5_a22;
                  c5_m_a22 = c5_a22;
                  c5_n_a22 = c5_a22;
                  c5_a22.re = c5_s.re * c5_k_a22.re + c5_s.im * c5_l_a22.im;
                  c5_a22.im = c5_s.re * c5_m_a22.im - c5_s.im * c5_n_a22.re;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2,
                            0) - 1)) - 1].re = c5_a12.re - c5_a22.re;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2,
                            0) - 1)) - 1].im = c5_a12.im - c5_a22.im;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2,
                            0) - 1)) - 1].re = c5_a21.re;
                  c5_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2,
                            0) - 1)) - 1].im = c5_a21.im;
                }

                c5_jm1 = c5_b_j;
                c5_b_j = c5_jp1;
              }
            }

            guard11 = true;
          }

          if (guard11 == true) {
            c5_jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2 == true) {
    if (c5_failed) {
      *c5_info = (real_T)c5_ilast;
      c5_b_ilast = c5_ilast;
      c5_m_b = c5_b_ilast;
      c5_n_b = c5_m_b;
      if (1 > c5_n_b) {
        c5_f_overflow = false;
      } else {
        c5_b_eml_switch_helper(chartInstance);
        c5_f_overflow = (c5_n_b > 2147483646);
      }

      if (c5_f_overflow) {
        c5_check_forloop_overflow_error(chartInstance, c5_f_overflow);
      }

      for (c5_k = 1; c5_k <= c5_b_ilast; c5_k++) {
        c5_b_k = c5_k;
        c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c5_b_k), 1, 3, 1, 0) - 1].re = c5_dc2.re;
        c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c5_b_k), 1, 3, 1, 0) - 1].im = c5_dc2.im;
        c5_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c5_b_k), 1, 3, 1, 0) - 1].re = c5_dc2.re;
        c5_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c5_b_k), 1, 3, 1, 0) - 1].im = c5_dc2.im;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1 == true) {
    c5_pb_a = c5_ilo;
    c5_qb_a = c5_pb_a - 1;
    c5_i63 = c5_qb_a;
    c5_o_b = c5_i63;
    c5_p_b = c5_o_b;
    if (1 > c5_p_b) {
      c5_g_overflow = false;
    } else {
      c5_b_eml_switch_helper(chartInstance);
      c5_g_overflow = (c5_p_b > 2147483646);
    }

    if (c5_g_overflow) {
      c5_check_forloop_overflow_error(chartInstance, c5_g_overflow);
    }

    for (c5_e_j = 1; c5_e_j <= c5_i63; c5_e_j++) {
      c5_b_j = c5_e_j;
      c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c5_b_j), 1, 3, 1, 0) - 1].re = c5_b_A
        [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) -
        1].re;
      c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c5_b_j), 1, 3, 1, 0) - 1].im = c5_b_A
        [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) -
        1].im;
    }

    *c5_info = 0.0;
  }
}

static real_T c5_eml_matlab_zlanhs(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T c5_ilo, int32_T c5_ihi)
{
  real_T c5_f;
  real_T c5_scale;
  real_T c5_sumsq;
  boolean_T c5_firstNonZero;
  int32_T c5_b_ilo;
  int32_T c5_b_ihi;
  int32_T c5_a;
  int32_T c5_b;
  int32_T c5_b_a;
  int32_T c5_b_b;
  boolean_T c5_overflow;
  int32_T c5_j;
  int32_T c5_b_j;
  int32_T c5_c_ilo;
  int32_T c5_c_a;
  int32_T c5_d_a;
  int32_T c5_c;
  int32_T c5_x;
  int32_T c5_y;
  int32_T c5_i64;
  int32_T c5_e_a;
  int32_T c5_c_b;
  int32_T c5_f_a;
  int32_T c5_d_b;
  boolean_T c5_b_overflow;
  int32_T c5_i;
  int32_T c5_b_i;
  creal_T c5_Aij;
  real_T c5_reAij;
  real_T c5_imAij;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_temp1;
  real_T c5_temp2;
  real_T c5_d_x;
  real_T c5_e_x;
  real_T c5_f_x;
  real_T c5_g_x;
  c5_f = 0.0;
  if (c5_ilo > c5_ihi) {
  } else {
    c5_scale = 0.0;
    c5_sumsq = 0.0;
    c5_firstNonZero = true;
    c5_b_ilo = c5_ilo;
    c5_b_ihi = c5_ihi;
    c5_a = c5_b_ilo;
    c5_b = c5_b_ihi;
    c5_b_a = c5_a;
    c5_b_b = c5_b;
    if (c5_b_a > c5_b_b) {
      c5_overflow = false;
    } else {
      c5_b_eml_switch_helper(chartInstance);
      c5_overflow = (c5_b_b > 2147483646);
    }

    if (c5_overflow) {
      c5_check_forloop_overflow_error(chartInstance, c5_overflow);
    }

    for (c5_j = c5_b_ilo; c5_j <= c5_b_ihi; c5_j++) {
      c5_b_j = c5_j;
      c5_c_ilo = c5_ilo;
      c5_c_a = c5_b_j;
      c5_d_a = c5_c_a;
      c5_c = c5_d_a;
      c5_x = c5_c + 1;
      c5_y = c5_ihi;
      c5_i64 = c5_x;
      if (c5_y < c5_i64) {
        c5_i64 = c5_y;
      }

      c5_e_a = c5_c_ilo;
      c5_c_b = c5_i64;
      c5_f_a = c5_e_a;
      c5_d_b = c5_c_b;
      if (c5_f_a > c5_d_b) {
        c5_b_overflow = false;
      } else {
        c5_b_eml_switch_helper(chartInstance);
        c5_b_overflow = (c5_d_b > 2147483646);
      }

      if (c5_b_overflow) {
        c5_check_forloop_overflow_error(chartInstance, c5_b_overflow);
      }

      for (c5_i = c5_c_ilo; c5_i <= c5_i64; c5_i++) {
        c5_b_i = c5_i;
        c5_Aij.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
        c5_Aij.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
        c5_reAij = c5_Aij.re;
        c5_imAij = c5_Aij.im;
        if (c5_reAij != 0.0) {
          c5_b_x = c5_reAij;
          c5_c_x = c5_b_x;
          c5_temp1 = muDoubleScalarAbs(c5_c_x);
          if (c5_firstNonZero) {
            c5_sumsq = 1.0;
            c5_scale = c5_temp1;
            c5_firstNonZero = false;
          } else if (c5_scale < c5_temp1) {
            c5_temp2 = c5_scale / c5_temp1;
            c5_sumsq = 1.0 + c5_sumsq * c5_temp2 * c5_temp2;
            c5_scale = c5_temp1;
          } else {
            c5_temp2 = c5_temp1 / c5_scale;
            c5_sumsq += c5_temp2 * c5_temp2;
          }
        }

        if (c5_imAij != 0.0) {
          c5_d_x = c5_imAij;
          c5_e_x = c5_d_x;
          c5_temp1 = muDoubleScalarAbs(c5_e_x);
          if (c5_firstNonZero) {
            c5_sumsq = 1.0;
            c5_scale = c5_temp1;
            c5_firstNonZero = false;
          } else if (c5_scale < c5_temp1) {
            c5_temp2 = c5_scale / c5_temp1;
            c5_sumsq = 1.0 + c5_sumsq * c5_temp2 * c5_temp2;
            c5_scale = c5_temp1;
          } else {
            c5_temp2 = c5_temp1 / c5_scale;
            c5_sumsq += c5_temp2 * c5_temp2;
          }
        }
      }
    }

    c5_f_x = c5_sumsq;
    c5_g_x = c5_f_x;
    if (c5_g_x < 0.0) {
      c5_eml_error(chartInstance);
    }

    c5_g_x = muDoubleScalarSqrt(c5_g_x);
    c5_f = c5_scale * c5_g_x;
  }

  return c5_f;
}

static int32_T c5_mod(SFc5_Model_01InstanceStruct *chartInstance, int32_T c5_x)
{
  int32_T c5_b_x;
  int32_T c5_t;
  c5_b_x = c5_x;
  c5_t = c5_div_s32(chartInstance, c5_b_x, 10);
  c5_t *= 10;
  return c5_b_x - c5_t;
}

static creal_T c5_eml_div(SFc5_Model_01InstanceStruct *chartInstance, creal_T
  c5_x, real_T c5_y)
{
  creal_T c5_z;
  real_T c5_b_y;
  real_T c5_ar;
  real_T c5_ai;
  real_T c5_br;
  real_T c5_bi;
  real_T c5_brm;
  real_T c5_bim;
  real_T c5_s;
  real_T c5_d;
  real_T c5_nr;
  real_T c5_ni;
  real_T c5_sgnbr;
  real_T c5_sgnbi;
  (void)chartInstance;
  c5_b_y = c5_y;
  c5_ar = c5_x.re;
  c5_ai = c5_x.im;
  c5_br = c5_b_y;
  c5_bi = 0.0;
  if (c5_bi == 0.0) {
    if (c5_ai == 0.0) {
      c5_z.re = c5_ar / c5_br;
      c5_z.im = 0.0;
    } else if (c5_ar == 0.0) {
      c5_z.re = 0.0;
      c5_z.im = c5_ai / c5_br;
    } else {
      c5_z.re = c5_ar / c5_br;
      c5_z.im = c5_ai / c5_br;
    }
  } else if (c5_br == 0.0) {
    if (c5_ar == 0.0) {
      c5_z.re = c5_ai / c5_bi;
      c5_z.im = 0.0;
    } else if (c5_ai == 0.0) {
      c5_z.re = 0.0;
      c5_z.im = -(c5_ar / c5_bi);
    } else {
      c5_z.re = c5_ai / c5_bi;
      c5_z.im = -(c5_ar / c5_bi);
    }
  } else {
    c5_brm = muDoubleScalarAbs(c5_br);
    c5_bim = muDoubleScalarAbs(c5_bi);
    if (c5_brm > c5_bim) {
      c5_s = c5_bi / c5_br;
      c5_d = c5_br + c5_s * c5_bi;
      c5_nr = c5_ar + c5_s * c5_ai;
      c5_ni = c5_ai - c5_s * c5_ar;
      c5_z.re = c5_nr / c5_d;
      c5_z.im = c5_ni / c5_d;
    } else if (c5_bim == c5_brm) {
      if (c5_br > 0.0) {
        c5_sgnbr = 0.5;
      } else {
        c5_sgnbr = -0.5;
      }

      if (c5_bi > 0.0) {
        c5_sgnbi = 0.5;
      } else {
        c5_sgnbi = -0.5;
      }

      c5_nr = c5_ar * c5_sgnbr + c5_ai * c5_sgnbi;
      c5_ni = c5_ai * c5_sgnbr - c5_ar * c5_sgnbi;
      c5_z.re = c5_nr / c5_brm;
      c5_z.im = c5_ni / c5_brm;
    } else {
      c5_s = c5_br / c5_bi;
      c5_d = c5_bi + c5_s * c5_br;
      c5_nr = c5_s * c5_ar + c5_ai;
      c5_ni = c5_s * c5_ai - c5_ar;
      c5_z.re = c5_nr / c5_d;
      c5_z.im = c5_ni / c5_d;
    }
  }

  return c5_z;
}

static void c5_scalarEg(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static creal_T c5_sqrt(SFc5_Model_01InstanceStruct *chartInstance, creal_T c5_x)
{
  creal_T c5_b_x;
  c5_b_x = c5_x;
  c5_b_sqrt(chartInstance, &c5_b_x);
  return c5_b_x;
}

static void c5_realmax(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c5_b_eml_matlab_zlartg(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_f, creal_T c5_g, real_T *c5_cs, creal_T *c5_sn)
{
  real_T c5_x;
  real_T c5_b_x;
  real_T c5_y;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_b_y;
  real_T c5_e_x;
  real_T c5_c_y;
  real_T c5_d_y;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_e_y;
  real_T c5_h_x;
  real_T c5_i_x;
  real_T c5_f_y;
  real_T c5_j_x;
  real_T c5_g_y;
  real_T c5_h_y;
  real_T c5_k_x;
  real_T c5_i_y;
  real_T c5_scale;
  creal_T c5_fs;
  creal_T c5_gs;
  int32_T c5_count;
  real_T c5_rescaledir;
  int32_T c5_a;
  int32_T c5_b_a;
  static creal_T c5_dc4 = { 0.0, 0.0 };

  boolean_T c5_b_g;
  int32_T c5_c_a;
  int32_T c5_d_a;
  real_T c5_f2;
  real_T c5_g2;
  real_T c5_l_x;
  real_T c5_m_x;
  boolean_T c5_b_f;
  real_T c5_x1;
  real_T c5_x2;
  real_T c5_e_a;
  real_T c5_b;
  real_T c5_d;
  real_T c5_b_x1;
  real_T c5_b_x2;
  real_T c5_f_a;
  real_T c5_b_b;
  real_T c5_f2s;
  real_T c5_n_x;
  real_T c5_g2s;
  real_T c5_o_x;
  real_T c5_p_x;
  real_T c5_j_y;
  real_T c5_q_x;
  real_T c5_r_x;
  real_T c5_k_y;
  real_T c5_s_x;
  real_T c5_l_y;
  real_T c5_m_y;
  real_T c5_c_x1;
  real_T c5_c_x2;
  real_T c5_g_a;
  real_T c5_c_b;
  real_T c5_dr;
  real_T c5_di;
  real_T c5_d_x1;
  real_T c5_d_x2;
  real_T c5_h_a;
  real_T c5_d_b;
  creal_T c5_b_gs;
  real_T c5_t_x;
  creal_T c5_b_fs;
  creal_T c5_c_fs;
  creal_T c5_c_gs;
  creal_T c5_b_sn;
  int32_T c5_b_count;
  int32_T c5_e_b;
  int32_T c5_f_b;
  boolean_T c5_overflow;
  int32_T c5_c_count;
  int32_T c5_g_b;
  int32_T c5_h_b;
  boolean_T c5_b_overflow;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  c5_realmin(chartInstance);
  c5_eps(chartInstance);
  c5_eps(chartInstance);
  c5_x = c5_f.re;
  c5_b_x = c5_x;
  c5_y = muDoubleScalarAbs(c5_b_x);
  c5_c_x = c5_f.im;
  c5_d_x = c5_c_x;
  c5_b_y = muDoubleScalarAbs(c5_d_x);
  c5_e_x = c5_y;
  c5_c_y = c5_b_y;
  c5_d_y = c5_e_x;
  if (c5_c_y > c5_d_y) {
    c5_d_y = c5_c_y;
  }

  c5_f_x = c5_g.re;
  c5_g_x = c5_f_x;
  c5_e_y = muDoubleScalarAbs(c5_g_x);
  c5_h_x = c5_g.im;
  c5_i_x = c5_h_x;
  c5_f_y = muDoubleScalarAbs(c5_i_x);
  c5_j_x = c5_e_y;
  c5_g_y = c5_f_y;
  c5_h_y = c5_j_x;
  if (c5_g_y > c5_h_y) {
    c5_h_y = c5_g_y;
  }

  c5_k_x = c5_d_y;
  c5_i_y = c5_h_y;
  c5_scale = c5_k_x;
  if (c5_i_y > c5_scale) {
    c5_scale = c5_i_y;
  }

  c5_fs = c5_f;
  c5_gs = c5_g;
  c5_count = 0;
  c5_rescaledir = 0.0;
  guard1 = false;
  guard2 = false;
  if (c5_scale >= 7.4428285367870146E+137) {
    do {
      c5_a = c5_count;
      c5_b_a = c5_a + 1;
      c5_count = c5_b_a;
      c5_fs.re *= 1.3435752215134178E-138;
      c5_fs.im *= 1.3435752215134178E-138;
      c5_gs.re *= 1.3435752215134178E-138;
      c5_gs.im *= 1.3435752215134178E-138;
      c5_scale *= 1.3435752215134178E-138;
    } while (!(c5_scale < 7.4428285367870146E+137));

    c5_rescaledir = 1.0;
    guard1 = true;
  } else if (c5_scale <= 1.3435752215134178E-138) {
    c5_b_g = ((c5_g.re == c5_dc4.re) && (c5_g.im == c5_dc4.im));
    if (c5_b_g) {
      *c5_cs = 1.0;
      *c5_sn = c5_dc4;
    } else {
      do {
        c5_c_a = c5_count;
        c5_d_a = c5_c_a + 1;
        c5_count = c5_d_a;
        c5_fs.re *= 7.4428285367870146E+137;
        c5_fs.im *= 7.4428285367870146E+137;
        c5_gs.re *= 7.4428285367870146E+137;
        c5_gs.im *= 7.4428285367870146E+137;
        c5_scale *= 7.4428285367870146E+137;
      } while (!(c5_scale > 1.3435752215134178E-138));

      c5_rescaledir = -1.0;
      guard2 = true;
    }
  } else {
    guard2 = true;
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c5_f2 = c5_fs.re * c5_fs.re + c5_fs.im * c5_fs.im;
    c5_g2 = c5_gs.re * c5_gs.re + c5_gs.im * c5_gs.im;
    c5_l_x = c5_g2;
    c5_m_x = c5_l_x;
    if (1.0 > c5_m_x) {
      c5_m_x = 1.0;
    }

    if (c5_f2 <= c5_m_x * 2.0041683600089728E-292) {
      c5_b_f = ((c5_f.re == c5_dc4.re) && (c5_f.im == c5_dc4.im));
      if (c5_b_f) {
        *c5_cs = 0.0;
        c5_x1 = c5_gs.re;
        c5_x2 = c5_gs.im;
        c5_e_a = c5_x1;
        c5_b = c5_x2;
        c5_d = muDoubleScalarHypot(c5_e_a, c5_b);
        c5_sn->re = c5_gs.re / c5_d;
        c5_sn->im = -c5_gs.im / c5_d;
      } else {
        c5_b_x1 = c5_fs.re;
        c5_b_x2 = c5_fs.im;
        c5_f_a = c5_b_x1;
        c5_b_b = c5_b_x2;
        c5_f2s = muDoubleScalarHypot(c5_f_a, c5_b_b);
        c5_n_x = c5_g2;
        c5_g2s = c5_n_x;
        if (c5_g2s < 0.0) {
          c5_eml_error(chartInstance);
        }

        c5_g2s = muDoubleScalarSqrt(c5_g2s);
        *c5_cs = c5_f2s / c5_g2s;
        c5_o_x = c5_f.re;
        c5_p_x = c5_o_x;
        c5_j_y = muDoubleScalarAbs(c5_p_x);
        c5_q_x = c5_f.im;
        c5_r_x = c5_q_x;
        c5_k_y = muDoubleScalarAbs(c5_r_x);
        c5_s_x = c5_j_y;
        c5_l_y = c5_k_y;
        c5_m_y = c5_s_x;
        if (c5_l_y > c5_m_y) {
          c5_m_y = c5_l_y;
        }

        if (c5_m_y > 1.0) {
          c5_c_x1 = c5_f.re;
          c5_c_x2 = c5_f.im;
          c5_g_a = c5_c_x1;
          c5_c_b = c5_c_x2;
          c5_d = muDoubleScalarHypot(c5_g_a, c5_c_b);
          c5_fs.re = c5_f.re / c5_d;
          c5_fs.im = c5_f.im / c5_d;
        } else {
          c5_dr = 7.4428285367870146E+137 * c5_f.re;
          c5_di = 7.4428285367870146E+137 * c5_f.im;
          c5_d_x1 = c5_dr;
          c5_d_x2 = c5_di;
          c5_h_a = c5_d_x1;
          c5_d_b = c5_d_x2;
          c5_d = muDoubleScalarHypot(c5_h_a, c5_d_b);
          c5_fs.re = c5_dr / c5_d;
          c5_fs.im = c5_di / c5_d;
        }

        c5_b_gs.re = c5_gs.re / c5_g2s;
        c5_b_gs.im = -c5_gs.im / c5_g2s;
        c5_sn->re = c5_fs.re * c5_b_gs.re - c5_fs.im * c5_b_gs.im;
        c5_sn->im = c5_fs.re * c5_b_gs.im + c5_fs.im * c5_b_gs.re;
      }
    } else {
      c5_t_x = 1.0 + c5_g2 / c5_f2;
      c5_f2s = c5_t_x;
      if (c5_f2s < 0.0) {
        c5_eml_error(chartInstance);
      }

      c5_f2s = muDoubleScalarSqrt(c5_f2s);
      c5_b_fs = c5_fs;
      c5_c_fs = c5_fs;
      c5_fs.re = c5_f2s * c5_b_fs.re;
      c5_fs.im = c5_f2s * c5_c_fs.im;
      *c5_cs = 1.0 / c5_f2s;
      c5_d = c5_f2 + c5_g2;
      c5_sn->re = c5_fs.re / c5_d;
      c5_sn->im = c5_fs.im / c5_d;
      c5_c_gs.re = c5_gs.re;
      c5_c_gs.im = -c5_gs.im;
      c5_b_sn = *c5_sn;
      c5_sn->re = c5_b_sn.re * c5_c_gs.re - c5_b_sn.im * c5_c_gs.im;
      c5_sn->im = c5_b_sn.re * c5_c_gs.im + c5_b_sn.im * c5_c_gs.re;
      if (c5_rescaledir > 0.0) {
        c5_b_count = c5_count;
        c5_e_b = c5_b_count;
        c5_f_b = c5_e_b;
        if (1 > c5_f_b) {
          c5_overflow = false;
        } else {
          c5_b_eml_switch_helper(chartInstance);
          c5_overflow = (c5_f_b > 2147483646);
        }

        if (c5_overflow) {
          c5_check_forloop_overflow_error(chartInstance, c5_overflow);
        }
      } else {
        if (c5_rescaledir < 0.0) {
          c5_c_count = c5_count;
          c5_g_b = c5_c_count;
          c5_h_b = c5_g_b;
          if (1 > c5_h_b) {
            c5_b_overflow = false;
          } else {
            c5_b_eml_switch_helper(chartInstance);
            c5_b_overflow = (c5_h_b > 2147483646);
          }

          if (c5_b_overflow) {
            c5_check_forloop_overflow_error(chartInstance, c5_b_overflow);
          }
        }
      }
    }
  }
}

static void c5_b_eml_matlab_zlascl(SFc5_Model_01InstanceStruct *chartInstance,
  real_T c5_cfrom, real_T c5_cto, creal_T c5_A[3], creal_T c5_b_A[3])
{
  int32_T c5_i65;
  for (c5_i65 = 0; c5_i65 < 3; c5_i65++) {
    c5_b_A[c5_i65] = c5_A[c5_i65];
  }

  c5_d_eml_matlab_zlascl(chartInstance, c5_cfrom, c5_cto, c5_b_A);
}

static void c5_b_eml_div(SFc5_Model_01InstanceStruct *chartInstance, creal_T
  c5_x[3], creal_T c5_y[3], creal_T c5_z[3])
{
  int32_T c5_i66;
  real_T c5_ar;
  real_T c5_ai;
  real_T c5_br;
  real_T c5_bi;
  real_T c5_brm;
  real_T c5_bim;
  real_T c5_s;
  real_T c5_d;
  real_T c5_nr;
  real_T c5_ni;
  real_T c5_sgnbr;
  real_T c5_sgnbi;
  (void)chartInstance;
  for (c5_i66 = 0; c5_i66 < 3; c5_i66++) {
    c5_ar = c5_x[c5_i66].re;
    c5_ai = c5_x[c5_i66].im;
    c5_br = c5_y[c5_i66].re;
    c5_bi = c5_y[c5_i66].im;
    if (c5_bi == 0.0) {
      if (c5_ai == 0.0) {
        c5_z[c5_i66].re = c5_ar / c5_br;
        c5_z[c5_i66].im = 0.0;
      } else if (c5_ar == 0.0) {
        c5_z[c5_i66].re = 0.0;
        c5_z[c5_i66].im = c5_ai / c5_br;
      } else {
        c5_z[c5_i66].re = c5_ar / c5_br;
        c5_z[c5_i66].im = c5_ai / c5_br;
      }
    } else if (c5_br == 0.0) {
      if (c5_ar == 0.0) {
        c5_z[c5_i66].re = c5_ai / c5_bi;
        c5_z[c5_i66].im = 0.0;
      } else if (c5_ai == 0.0) {
        c5_z[c5_i66].re = 0.0;
        c5_z[c5_i66].im = -(c5_ar / c5_bi);
      } else {
        c5_z[c5_i66].re = c5_ai / c5_bi;
        c5_z[c5_i66].im = -(c5_ar / c5_bi);
      }
    } else {
      c5_brm = muDoubleScalarAbs(c5_br);
      c5_bim = muDoubleScalarAbs(c5_bi);
      if (c5_brm > c5_bim) {
        c5_s = c5_bi / c5_br;
        c5_d = c5_br + c5_s * c5_bi;
        c5_nr = c5_ar + c5_s * c5_ai;
        c5_ni = c5_ai - c5_s * c5_ar;
        c5_z[c5_i66].re = c5_nr / c5_d;
        c5_z[c5_i66].im = c5_ni / c5_d;
      } else if (c5_bim == c5_brm) {
        if (c5_br > 0.0) {
          c5_sgnbr = 0.5;
        } else {
          c5_sgnbr = -0.5;
        }

        if (c5_bi > 0.0) {
          c5_sgnbi = 0.5;
        } else {
          c5_sgnbi = -0.5;
        }

        c5_nr = c5_ar * c5_sgnbr + c5_ai * c5_sgnbi;
        c5_ni = c5_ai * c5_sgnbr - c5_ar * c5_sgnbi;
        c5_z[c5_i66].re = c5_nr / c5_brm;
        c5_z[c5_i66].im = c5_ni / c5_brm;
      } else {
        c5_s = c5_br / c5_bi;
        c5_d = c5_bi + c5_s * c5_br;
        c5_nr = c5_s * c5_ar + c5_ai;
        c5_ni = c5_s * c5_ai - c5_ar;
        c5_z[c5_i66].re = c5_nr / c5_d;
        c5_z[c5_i66].im = c5_ni / c5_d;
      }
    }
  }
}

static void c5_eml_warning(SFc5_Model_01InstanceStruct *chartInstance)
{
  int32_T c5_i67;
  static char_T c5_varargin_1[26] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'f', 'a', 'i',
    'l', 'e', 'd' };

  char_T c5_u[26];
  const mxArray *c5_y = NULL;
  (void)chartInstance;
  for (c5_i67 = 0; c5_i67 < 26; c5_i67++) {
    c5_u[c5_i67] = c5_varargin_1[c5_i67];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 26), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c5_y));
}

static void c5_b_eml_warning(SFc5_Model_01InstanceStruct *chartInstance)
{
  int32_T c5_i68;
  static char_T c5_varargin_1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'n', 'o', 'n',
    'c', 'o', 'n', 'v', 'e', 'r', 'g', 'e', 'n', 'c', 'e' };

  char_T c5_u[34];
  const mxArray *c5_y = NULL;
  (void)chartInstance;
  for (c5_i68 = 0; c5_i68 < 34; c5_i68++) {
    c5_u[c5_i68] = c5_varargin_1[c5_i68];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 34), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c5_y));
}

static real_T c5_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_a)
{
  real_T c5_b_a;
  real_T c5_c_a;
  real_T c5_ak;
  real_T c5_d_a;
  c5_b_a = c5_a;
  c5_c_a = c5_b_a;
  c5_eml_scalar_eg(chartInstance);
  c5_ak = c5_c_a;
  c5_d_a = c5_ak;
  c5_eml_scalar_eg(chartInstance);
  return c5_d_a * c5_d_a;
}

static void c5_inv(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x[9],
                   real_T c5_y[9])
{
  int32_T c5_i69;
  real_T c5_b_x[9];
  int32_T c5_i70;
  real_T c5_c_x[9];
  real_T c5_n1x;
  int32_T c5_i71;
  real_T c5_b_y[9];
  real_T c5_n1xinv;
  real_T c5_rc;
  real_T c5_d_x;
  boolean_T c5_b;
  real_T c5_e_x;
  int32_T c5_i72;
  static char_T c5_cv4[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T c5_u[8];
  const mxArray *c5_c_y = NULL;
  real_T c5_b_u;
  const mxArray *c5_d_y = NULL;
  real_T c5_c_u;
  const mxArray *c5_e_y = NULL;
  real_T c5_d_u;
  const mxArray *c5_f_y = NULL;
  char_T c5_str[14];
  int32_T c5_i73;
  char_T c5_b_str[14];
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  for (c5_i69 = 0; c5_i69 < 9; c5_i69++) {
    c5_b_x[c5_i69] = c5_x[c5_i69];
  }

  c5_inv3x3(chartInstance, c5_b_x, c5_y);
  for (c5_i70 = 0; c5_i70 < 9; c5_i70++) {
    c5_c_x[c5_i70] = c5_x[c5_i70];
  }

  c5_n1x = c5_b_norm(chartInstance, c5_c_x);
  for (c5_i71 = 0; c5_i71 < 9; c5_i71++) {
    c5_b_y[c5_i71] = c5_y[c5_i71];
  }

  c5_n1xinv = c5_b_norm(chartInstance, c5_b_y);
  c5_rc = 1.0 / (c5_n1x * c5_n1xinv);
  guard1 = false;
  guard2 = false;
  if (c5_n1x == 0.0) {
    guard2 = true;
  } else if (c5_n1xinv == 0.0) {
    guard2 = true;
  } else if (c5_rc == 0.0) {
    guard1 = true;
  } else {
    c5_d_x = c5_rc;
    c5_b = muDoubleScalarIsNaN(c5_d_x);
    guard3 = false;
    if (c5_b) {
      guard3 = true;
    } else {
      c5_eps(chartInstance);
      if (c5_rc < 2.2204460492503131E-16) {
        guard3 = true;
      }
    }

    if (guard3 == true) {
      c5_e_x = c5_rc;
      for (c5_i72 = 0; c5_i72 < 8; c5_i72++) {
        c5_u[c5_i72] = c5_cv4[c5_i72];
      }

      c5_c_y = NULL;
      sf_mex_assign(&c5_c_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 8),
                    false);
      c5_b_u = 14.0;
      c5_d_y = NULL;
      sf_mex_assign(&c5_d_y, sf_mex_create("y", &c5_b_u, 0, 0U, 0U, 0U, 0),
                    false);
      c5_c_u = 6.0;
      c5_e_y = NULL;
      sf_mex_assign(&c5_e_y, sf_mex_create("y", &c5_c_u, 0, 0U, 0U, 0U, 0),
                    false);
      c5_d_u = c5_e_x;
      c5_f_y = NULL;
      sf_mex_assign(&c5_f_y, sf_mex_create("y", &c5_d_u, 0, 0U, 0U, 0U, 0),
                    false);
      c5_e_emlrt_marshallIn(chartInstance, sf_mex_call_debug
                            (sfGlobalDebugInstanceStruct, "sprintf", 1U, 2U, 14,
        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "sprintf", 1U, 3U, 14,
                          c5_c_y, 14, c5_d_y, 14, c5_e_y), 14, c5_f_y),
                            "sprintf", c5_str);
      for (c5_i73 = 0; c5_i73 < 14; c5_i73++) {
        c5_b_str[c5_i73] = c5_str[c5_i73];
      }

      c5_d_eml_warning(chartInstance, c5_b_str);
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c5_c_eml_warning(chartInstance);
  }
}

static void c5_inv3x3(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x[9],
                      real_T c5_y[9])
{
  int32_T c5_p1;
  int32_T c5_p2;
  int32_T c5_p3;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_absx11;
  real_T c5_d_x;
  real_T c5_e_x;
  real_T c5_absx21;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_absx31;
  real_T c5_t1;
  real_T c5_h_x;
  real_T c5_b_y;
  real_T c5_i_x;
  real_T c5_c_y;
  real_T c5_z;
  real_T c5_j_x;
  real_T c5_d_y;
  real_T c5_k_x;
  real_T c5_e_y;
  real_T c5_b_z;
  real_T c5_l_x;
  real_T c5_m_x;
  real_T c5_f_y;
  real_T c5_n_x;
  real_T c5_o_x;
  real_T c5_g_y;
  int32_T c5_itmp;
  real_T c5_p_x;
  real_T c5_h_y;
  real_T c5_q_x;
  real_T c5_i_y;
  real_T c5_c_z;
  real_T c5_r_x;
  real_T c5_j_y;
  real_T c5_s_x;
  real_T c5_k_y;
  real_T c5_t3;
  real_T c5_t_x;
  real_T c5_l_y;
  real_T c5_u_x;
  real_T c5_m_y;
  real_T c5_t2;
  int32_T c5_a;
  int32_T c5_b_a;
  int32_T c5_c;
  real_T c5_v_x;
  real_T c5_n_y;
  real_T c5_w_x;
  real_T c5_o_y;
  real_T c5_d_z;
  int32_T c5_c_a;
  int32_T c5_d_a;
  int32_T c5_b_c;
  int32_T c5_e_a;
  int32_T c5_f_a;
  int32_T c5_c_c;
  real_T c5_x_x;
  real_T c5_p_y;
  real_T c5_y_x;
  real_T c5_q_y;
  real_T c5_ab_x;
  real_T c5_r_y;
  real_T c5_bb_x;
  real_T c5_s_y;
  int32_T c5_g_a;
  int32_T c5_h_a;
  int32_T c5_d_c;
  real_T c5_cb_x;
  real_T c5_t_y;
  real_T c5_db_x;
  real_T c5_u_y;
  real_T c5_e_z;
  int32_T c5_i_a;
  int32_T c5_j_a;
  int32_T c5_e_c;
  int32_T c5_k_a;
  int32_T c5_l_a;
  int32_T c5_f_c;
  real_T c5_v_y;
  real_T c5_w_y;
  real_T c5_eb_x;
  real_T c5_x_y;
  real_T c5_fb_x;
  real_T c5_y_y;
  int32_T c5_m_a;
  int32_T c5_n_a;
  int32_T c5_g_c;
  real_T c5_gb_x;
  real_T c5_ab_y;
  real_T c5_hb_x;
  real_T c5_bb_y;
  real_T c5_f_z;
  int32_T c5_o_a;
  int32_T c5_p_a;
  int32_T c5_h_c;
  int32_T c5_q_a;
  int32_T c5_r_a;
  int32_T c5_i_c;
  boolean_T guard1 = false;
  (void)chartInstance;
  c5_p1 = 0;
  c5_p2 = 3;
  c5_p3 = 6;
  c5_b_x = c5_x[0];
  c5_c_x = c5_b_x;
  c5_absx11 = muDoubleScalarAbs(c5_c_x);
  c5_d_x = c5_x[1];
  c5_e_x = c5_d_x;
  c5_absx21 = muDoubleScalarAbs(c5_e_x);
  c5_f_x = c5_x[2];
  c5_g_x = c5_f_x;
  c5_absx31 = muDoubleScalarAbs(c5_g_x);
  guard1 = false;
  if (c5_absx21 > c5_absx11) {
    if (c5_absx21 > c5_absx31) {
      c5_p1 = 3;
      c5_p2 = 0;
      c5_t1 = c5_x[0];
      c5_x[0] = c5_x[1];
      c5_x[1] = c5_t1;
      c5_t1 = c5_x[3];
      c5_x[3] = c5_x[4];
      c5_x[4] = c5_t1;
      c5_t1 = c5_x[6];
      c5_x[6] = c5_x[7];
      c5_x[7] = c5_t1;
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1 == true) {
    if (c5_absx31 > c5_absx11) {
      c5_p1 = 6;
      c5_p3 = 0;
      c5_t1 = c5_x[0];
      c5_x[0] = c5_x[2];
      c5_x[2] = c5_t1;
      c5_t1 = c5_x[3];
      c5_x[3] = c5_x[5];
      c5_x[5] = c5_t1;
      c5_t1 = c5_x[6];
      c5_x[6] = c5_x[8];
      c5_x[8] = c5_t1;
    }
  }

  c5_h_x = c5_x[1];
  c5_b_y = c5_x[0];
  c5_i_x = c5_h_x;
  c5_c_y = c5_b_y;
  c5_z = c5_i_x / c5_c_y;
  c5_x[1] = c5_z;
  c5_j_x = c5_x[2];
  c5_d_y = c5_x[0];
  c5_k_x = c5_j_x;
  c5_e_y = c5_d_y;
  c5_b_z = c5_k_x / c5_e_y;
  c5_x[2] = c5_b_z;
  c5_x[4] -= c5_x[1] * c5_x[3];
  c5_x[5] -= c5_x[2] * c5_x[3];
  c5_x[7] -= c5_x[1] * c5_x[6];
  c5_x[8] -= c5_x[2] * c5_x[6];
  c5_l_x = c5_x[5];
  c5_m_x = c5_l_x;
  c5_f_y = muDoubleScalarAbs(c5_m_x);
  c5_n_x = c5_x[4];
  c5_o_x = c5_n_x;
  c5_g_y = muDoubleScalarAbs(c5_o_x);
  if (c5_f_y > c5_g_y) {
    c5_itmp = c5_p2;
    c5_p2 = c5_p3;
    c5_p3 = c5_itmp;
    c5_t1 = c5_x[1];
    c5_x[1] = c5_x[2];
    c5_x[2] = c5_t1;
    c5_t1 = c5_x[4];
    c5_x[4] = c5_x[5];
    c5_x[5] = c5_t1;
    c5_t1 = c5_x[7];
    c5_x[7] = c5_x[8];
    c5_x[8] = c5_t1;
  }

  c5_p_x = c5_x[5];
  c5_h_y = c5_x[4];
  c5_q_x = c5_p_x;
  c5_i_y = c5_h_y;
  c5_c_z = c5_q_x / c5_i_y;
  c5_x[5] = c5_c_z;
  c5_x[8] -= c5_x[5] * c5_x[7];
  c5_r_x = c5_x[5] * c5_x[1] - c5_x[2];
  c5_j_y = c5_x[8];
  c5_s_x = c5_r_x;
  c5_k_y = c5_j_y;
  c5_t3 = c5_s_x / c5_k_y;
  c5_t_x = -(c5_x[1] + c5_x[7] * c5_t3);
  c5_l_y = c5_x[4];
  c5_u_x = c5_t_x;
  c5_m_y = c5_l_y;
  c5_t2 = c5_u_x / c5_m_y;
  c5_a = c5_p1;
  c5_b_a = c5_a + 1;
  c5_c = c5_b_a;
  c5_v_x = (1.0 - c5_x[3] * c5_t2) - c5_x[6] * c5_t3;
  c5_n_y = c5_x[0];
  c5_w_x = c5_v_x;
  c5_o_y = c5_n_y;
  c5_d_z = c5_w_x / c5_o_y;
  c5_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_c), 1, 9, 1, 0) - 1] = c5_d_z;
  c5_c_a = c5_p1;
  c5_d_a = c5_c_a + 2;
  c5_b_c = c5_d_a;
  c5_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_b_c), 1, 9, 1, 0) - 1] = c5_t2;
  c5_e_a = c5_p1;
  c5_f_a = c5_e_a + 3;
  c5_c_c = c5_f_a;
  c5_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_c_c), 1, 9, 1, 0) - 1] = c5_t3;
  c5_x_x = -c5_x[5];
  c5_p_y = c5_x[8];
  c5_y_x = c5_x_x;
  c5_q_y = c5_p_y;
  c5_t3 = c5_y_x / c5_q_y;
  c5_ab_x = 1.0 - c5_x[7] * c5_t3;
  c5_r_y = c5_x[4];
  c5_bb_x = c5_ab_x;
  c5_s_y = c5_r_y;
  c5_t2 = c5_bb_x / c5_s_y;
  c5_g_a = c5_p2;
  c5_h_a = c5_g_a + 1;
  c5_d_c = c5_h_a;
  c5_cb_x = -(c5_x[3] * c5_t2 + c5_x[6] * c5_t3);
  c5_t_y = c5_x[0];
  c5_db_x = c5_cb_x;
  c5_u_y = c5_t_y;
  c5_e_z = c5_db_x / c5_u_y;
  c5_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_d_c), 1, 9, 1, 0) - 1] = c5_e_z;
  c5_i_a = c5_p2;
  c5_j_a = c5_i_a + 2;
  c5_e_c = c5_j_a;
  c5_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_e_c), 1, 9, 1, 0) - 1] = c5_t2;
  c5_k_a = c5_p2;
  c5_l_a = c5_k_a + 3;
  c5_f_c = c5_l_a;
  c5_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_f_c), 1, 9, 1, 0) - 1] = c5_t3;
  c5_v_y = c5_x[8];
  c5_w_y = c5_v_y;
  c5_t3 = 1.0 / c5_w_y;
  c5_eb_x = -c5_x[7] * c5_t3;
  c5_x_y = c5_x[4];
  c5_fb_x = c5_eb_x;
  c5_y_y = c5_x_y;
  c5_t2 = c5_fb_x / c5_y_y;
  c5_m_a = c5_p3;
  c5_n_a = c5_m_a + 1;
  c5_g_c = c5_n_a;
  c5_gb_x = -(c5_x[3] * c5_t2 + c5_x[6] * c5_t3);
  c5_ab_y = c5_x[0];
  c5_hb_x = c5_gb_x;
  c5_bb_y = c5_ab_y;
  c5_f_z = c5_hb_x / c5_bb_y;
  c5_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_g_c), 1, 9, 1, 0) - 1] = c5_f_z;
  c5_o_a = c5_p3;
  c5_p_a = c5_o_a + 2;
  c5_h_c = c5_p_a;
  c5_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_h_c), 1, 9, 1, 0) - 1] = c5_t2;
  c5_q_a = c5_p3;
  c5_r_a = c5_q_a + 3;
  c5_i_c = c5_r_a;
  c5_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_i_c), 1, 9, 1, 0) - 1] = c5_t3;
}

static real_T c5_b_norm(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_x
  [9])
{
  real_T c5_y;
  int32_T c5_j;
  real_T c5_b_j;
  real_T c5_s;
  int32_T c5_i;
  real_T c5_b_i;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_b_y;
  real_T c5_d_x;
  boolean_T c5_b;
  boolean_T exitg1;
  (void)chartInstance;
  c5_y = 0.0;
  c5_j = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c5_j < 3)) {
    c5_b_j = 1.0 + (real_T)c5_j;
    c5_s = 0.0;
    for (c5_i = 0; c5_i < 3; c5_i++) {
      c5_b_i = 1.0 + (real_T)c5_i;
      c5_b_x = c5_x[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c5_b_i), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1];
      c5_c_x = c5_b_x;
      c5_b_y = muDoubleScalarAbs(c5_c_x);
      c5_s += c5_b_y;
    }

    c5_d_x = c5_s;
    c5_b = muDoubleScalarIsNaN(c5_d_x);
    if (c5_b) {
      c5_y = rtNaN;
      exitg1 = true;
    } else {
      if (c5_s > c5_y) {
        c5_y = c5_s;
      }

      c5_j++;
    }
  }

  return c5_y;
}

static void c5_c_eml_warning(SFc5_Model_01InstanceStruct *chartInstance)
{
  int32_T c5_i74;
  static char_T c5_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c5_u[27];
  const mxArray *c5_y = NULL;
  (void)chartInstance;
  for (c5_i74 = 0; c5_i74 < 27; c5_i74++) {
    c5_u[c5_i74] = c5_varargin_1[c5_i74];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 27), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c5_y));
}

static void c5_d_eml_warning(SFc5_Model_01InstanceStruct *chartInstance, char_T
  c5_varargin_2[14])
{
  int32_T c5_i75;
  static char_T c5_varargin_1[33] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 'i', 'l', 'l', 'C', 'o', 'n', 'd', 'i', 't', 'i',
    'o', 'n', 'e', 'd', 'M', 'a', 't', 'r', 'i', 'x' };

  char_T c5_u[33];
  const mxArray *c5_y = NULL;
  int32_T c5_i76;
  char_T c5_b_u[14];
  const mxArray *c5_b_y = NULL;
  (void)chartInstance;
  for (c5_i75 = 0; c5_i75 < 33; c5_i75++) {
    c5_u[c5_i75] = c5_varargin_1[c5_i75];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 33), false);
  for (c5_i76 = 0; c5_i76 < 14; c5_i76++) {
    c5_b_u[c5_i76] = c5_varargin_2[c5_i76];
  }

  c5_b_y = NULL;
  sf_mex_assign(&c5_b_y, sf_mex_create("y", c5_b_u, 10, 0U, 1U, 0U, 2, 1, 14),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c5_y, 14, c5_b_y));
}

static real_T c5_b_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_a)
{
  real_T c5_b_a;
  real_T c5_c_a;
  real_T c5_ak;
  real_T c5_d_a;
  real_T c5_ar;
  c5_b_a = c5_a;
  c5_c_a = c5_b_a;
  c5_eml_scalar_eg(chartInstance);
  c5_ak = c5_c_a;
  c5_d_a = c5_ak;
  c5_eml_scalar_eg(chartInstance);
  c5_ar = c5_d_a;
  return muDoubleScalarPower(c5_ar, -2.0);
}

static real_T c5_c_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_a)
{
  real_T c5_b_a;
  real_T c5_c_a;
  real_T c5_ak;
  real_T c5_d_a;
  real_T c5_ar;
  c5_b_a = c5_a;
  c5_c_a = c5_b_a;
  c5_eml_scalar_eg(chartInstance);
  c5_ak = c5_c_a;
  c5_d_a = c5_ak;
  c5_eml_scalar_eg(chartInstance);
  c5_ar = c5_d_a;
  return muDoubleScalarPower(c5_ar, 4.0);
}

static real_T c5_d_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_a)
{
  real_T c5_b_a;
  real_T c5_c_a;
  real_T c5_ak;
  real_T c5_d_a;
  real_T c5_ar;
  c5_b_a = c5_a;
  c5_c_a = c5_b_a;
  c5_eml_scalar_eg(chartInstance);
  c5_ak = c5_c_a;
  c5_d_a = c5_ak;
  c5_eml_scalar_eg(chartInstance);
  c5_ar = c5_d_a;
  return muDoubleScalarPower(c5_ar, 3.0);
}

static real_T c5_e_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_a)
{
  real_T c5_b_a;
  real_T c5_c_a;
  real_T c5_ak;
  real_T c5_d_a;
  real_T c5_B;
  real_T c5_y;
  real_T c5_b_y;
  real_T c5_c_y;
  c5_b_a = c5_a;
  c5_c_a = c5_b_a;
  c5_eml_scalar_eg(chartInstance);
  c5_ak = c5_c_a;
  c5_d_a = c5_ak;
  c5_eml_scalar_eg(chartInstance);
  c5_B = c5_d_a;
  c5_y = c5_B;
  c5_b_y = c5_y;
  c5_c_y = c5_b_y;
  return 1.0 / c5_c_y;
}

static void c5_f_mpower(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_a
  [9], real_T c5_b, real_T c5_c[9])
{
  real_T c5_x;
  real_T c5_b_x;
  int32_T c5_i77;
  real_T c5_b_a[9];
  real_T c5_b_b;
  int32_T c5_i78;
  int32_T c5_k;
  real_T c5_b_k;
  real_T c5_c_x;
  real_T c5_e;
  boolean_T c5_firstmult;
  real_T c5_d_x;
  real_T c5_ed2;
  int32_T c5_i79;
  int32_T c5_i80;
  real_T c5_c_a[9];
  int32_T c5_i81;
  int32_T c5_i82;
  real_T c5_d_a[9];
  int32_T c5_i83;
  real_T c5_e_a[9];
  int32_T c5_i84;
  real_T c5_b_c[9];
  int32_T c5_i85;
  int32_T c5_i86;
  int32_T c5_i87;
  real_T c5_f_a[9];
  int32_T c5_i88;
  real_T c5_g_a[9];
  real_T c5_c_b;
  int32_T c5_i89;
  real_T c5_h_a[9];
  creal_T c5_D[9];
  creal_T c5_V[9];
  int32_T c5_j;
  real_T c5_b_j;
  creal_T c5_b_D;
  creal_T c5_r;
  int32_T c5_i;
  real_T c5_b_i;
  int32_T c5_i90;
  creal_T c5_b_V[9];
  int32_T c5_i91;
  creal_T c5_c_D[9];
  int32_T c5_i92;
  int32_T exitg1;
  c5_x = c5_b;
  c5_b_x = c5_x;
  c5_b_x = muDoubleScalarFloor(c5_b_x);
  if (c5_b_x == c5_b) {
    for (c5_i77 = 0; c5_i77 < 9; c5_i77++) {
      c5_b_a[c5_i77] = c5_a[c5_i77];
    }

    c5_b_b = c5_b;
    c5_b_eml_scalar_eg(chartInstance);
    if (c5_b_b == 0.0) {
      for (c5_i78 = 0; c5_i78 < 9; c5_i78++) {
        c5_c[c5_i78] = 0.0;
      }

      for (c5_k = 0; c5_k < 3; c5_k++) {
        c5_b_k = 1.0 + (real_T)c5_k;
        c5_c[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c5_b_k), 1, 3, 2, 0) - 1)) - 1] =
          1.0;
      }
    } else {
      c5_c_x = c5_b_b;
      c5_e = muDoubleScalarAbs(c5_c_x);
      c5_firstmult = true;
      do {
        exitg1 = 0;
        c5_d_x = c5_e / 2.0;
        c5_ed2 = c5_d_x;
        c5_ed2 = muDoubleScalarFloor(c5_ed2);
        if (2.0 * c5_ed2 != c5_e) {
          if (c5_firstmult) {
            for (c5_i79 = 0; c5_i79 < 9; c5_i79++) {
              c5_c[c5_i79] = c5_b_a[c5_i79];
            }

            c5_firstmult = false;
          } else {
            for (c5_i80 = 0; c5_i80 < 9; c5_i80++) {
              c5_c_a[c5_i80] = c5_c[c5_i80];
            }

            c5_c_eml_scalar_eg(chartInstance);
            c5_c_eml_scalar_eg(chartInstance);
            for (c5_i81 = 0; c5_i81 < 9; c5_i81++) {
              c5_c[c5_i81] = 0.0;
            }

            for (c5_i82 = 0; c5_i82 < 9; c5_i82++) {
              c5_d_a[c5_i82] = c5_c_a[c5_i82];
            }

            for (c5_i83 = 0; c5_i83 < 9; c5_i83++) {
              c5_e_a[c5_i83] = c5_b_a[c5_i83];
            }

            c5_b_eml_xgemm(chartInstance, c5_d_a, c5_e_a, c5_c);
          }
        }

        if (c5_ed2 == 0.0) {
          exitg1 = 1;
        } else {
          c5_e = c5_ed2;
          for (c5_i85 = 0; c5_i85 < 9; c5_i85++) {
            c5_c_a[c5_i85] = c5_b_a[c5_i85];
          }

          c5_c_eml_scalar_eg(chartInstance);
          c5_c_eml_scalar_eg(chartInstance);
          for (c5_i86 = 0; c5_i86 < 9; c5_i86++) {
            c5_b_a[c5_i86] = 0.0;
          }

          for (c5_i87 = 0; c5_i87 < 9; c5_i87++) {
            c5_f_a[c5_i87] = c5_c_a[c5_i87];
          }

          for (c5_i88 = 0; c5_i88 < 9; c5_i88++) {
            c5_g_a[c5_i88] = c5_c_a[c5_i88];
          }

          c5_b_eml_xgemm(chartInstance, c5_f_a, c5_g_a, c5_b_a);
        }
      } while (exitg1 == 0);

      if (c5_b_b < 0.0) {
        for (c5_i84 = 0; c5_i84 < 9; c5_i84++) {
          c5_b_c[c5_i84] = c5_c[c5_i84];
        }

        c5_inv(chartInstance, c5_b_c, c5_c);
      }
    }
  } else {
    c5_c_b = c5_b;
    c5_b_eml_scalar_eg(chartInstance);
    c5_b_eml_error(chartInstance);
    for (c5_i89 = 0; c5_i89 < 9; c5_i89++) {
      c5_h_a[c5_i89] = c5_a[c5_i89];
    }

    c5_b_eig(chartInstance, c5_h_a, c5_V, c5_D);
    for (c5_j = 0; c5_j < 3; c5_j++) {
      c5_b_j = 1.0 + (real_T)c5_j;
      c5_b_D.re = c5_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) + 3 *
                        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c5_b_D.im = c5_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) + 3 *
                        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
      c5_r = c5_power(chartInstance, c5_b_D, c5_c_b);
      for (c5_i = 0; c5_i < 3; c5_i++) {
        c5_b_i = 1.0 + (real_T)c5_i;
        c5_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c5_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].
          re = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c5_b_i), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].re * c5_r.re -
          c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  c5_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                  (int32_T)_SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1]
          .im * c5_r.im;
        c5_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c5_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].
          im = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c5_b_i), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].re * c5_r.im +
          c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  c5_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                  (int32_T)_SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1]
          .im * c5_r.re;
      }
    }

    for (c5_i90 = 0; c5_i90 < 9; c5_i90++) {
      c5_b_V[c5_i90] = c5_V[c5_i90];
    }

    for (c5_i91 = 0; c5_i91 < 9; c5_i91++) {
      c5_c_D[c5_i91] = c5_D[c5_i91];
    }

    c5_eml_lusolve(chartInstance, c5_b_V, c5_c_D, c5_D);
    for (c5_i92 = 0; c5_i92 < 9; c5_i92++) {
      c5_c[c5_i92] = c5_D[c5_i92].re;
    }
  }
}

static void c5_b_eml_scalar_eg(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c5_c_eml_scalar_eg(SFc5_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c5_eml_xgemm(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_A[9], real_T c5_B[9], real_T c5_C[9], real_T c5_b_C[9])
{
  int32_T c5_i93;
  int32_T c5_i94;
  real_T c5_b_A[9];
  int32_T c5_i95;
  real_T c5_b_B[9];
  for (c5_i93 = 0; c5_i93 < 9; c5_i93++) {
    c5_b_C[c5_i93] = c5_C[c5_i93];
  }

  for (c5_i94 = 0; c5_i94 < 9; c5_i94++) {
    c5_b_A[c5_i94] = c5_A[c5_i94];
  }

  for (c5_i95 = 0; c5_i95 < 9; c5_i95++) {
    c5_b_B[c5_i95] = c5_B[c5_i95];
  }

  c5_b_eml_xgemm(chartInstance, c5_b_A, c5_b_B, c5_b_C);
}

static void c5_b_eml_error(SFc5_Model_01InstanceStruct *chartInstance)
{
  int32_T c5_i96;
  static char_T c5_cv5[37] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'm', 'p', 'o', 'w', 'e', 'r', '_', 'n', 'e', 'e', 'd',
    'C', 'o', 'm', 'p', 'l', 'e', 'x', 'I', 'n', 'p', 'u', 't' };

  char_T c5_u[37];
  const mxArray *c5_y = NULL;
  (void)chartInstance;
  for (c5_i96 = 0; c5_i96 < 37; c5_i96++) {
    c5_u[c5_i96] = c5_cv5[c5_i96];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 37), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c5_y));
}

static void c5_b_eig(SFc5_Model_01InstanceStruct *chartInstance, real_T c5_A[9],
                     creal_T c5_V[9], creal_T c5_D[9])
{
  int32_T c5_i97;
  static creal_T c5_dc5 = { 0.0, 0.0 };

  creal_T c5_b_A[9];
  int32_T c5_i98;
  creal_T c5_c_A[9];
  creal_T c5_beta1[3];
  creal_T c5_alpha1[3];
  real_T c5_info;
  real_T c5_b_info;
  int32_T c5_coltop;
  int32_T c5_b_coltop;
  int32_T c5_ix0;
  int32_T c5_b_ix0;
  int32_T c5_c_ix0;
  real_T c5_colnorm;
  real_T c5_scale;
  int32_T c5_kstart;
  int32_T c5_a;
  int32_T c5_kend;
  int32_T c5_b_kstart;
  int32_T c5_b_kend;
  int32_T c5_b_a;
  int32_T c5_b;
  int32_T c5_c_a;
  int32_T c5_b_b;
  boolean_T c5_overflow;
  int32_T c5_k;
  int32_T c5_b_k;
  real_T c5_x;
  real_T c5_b_x;
  real_T c5_absxk;
  real_T c5_t;
  real_T c5_c_x;
  real_T c5_d_x;
  int32_T c5_c_coltop;
  int32_T c5_d_a;
  int32_T c5_e_a;
  int32_T c5_i99;
  int32_T c5_f_a;
  int32_T c5_c_b;
  int32_T c5_g_a;
  int32_T c5_d_b;
  boolean_T c5_b_overflow;
  int32_T c5_j;
  int32_T c5_b_j;
  creal_T c5_d_A;
  real_T c5_B;
  real_T c5_y;
  real_T c5_c_info;
  real_T c5_d_info;
  int32_T c5_i100;
  creal_T c5_b_alpha1[3];
  int32_T c5_i101;
  creal_T c5_b_beta1[3];
  int32_T c5_i102;
  int32_T c5_c_j;
  int32_T c5_d_j;
  for (c5_i97 = 0; c5_i97 < 9; c5_i97++) {
    c5_b_A[c5_i97].re = c5_A[c5_i97] + c5_dc5.re;
    c5_b_A[c5_i97].im = c5_dc5.im;
  }

  for (c5_i98 = 0; c5_i98 < 9; c5_i98++) {
    c5_c_A[c5_i98] = c5_b_A[c5_i98];
  }

  c5_eml_matlab_zggev(chartInstance, c5_c_A, &c5_info, c5_alpha1, c5_beta1, c5_V);
  c5_b_info = c5_info;
  for (c5_coltop = 1; c5_coltop < 8; c5_coltop += 3) {
    c5_b_coltop = c5_coltop;
    c5_ix0 = c5_b_coltop;
    c5_b_ix0 = c5_ix0;
    c5_eml_switch_helper(chartInstance);
    c5_c_ix0 = c5_b_ix0;
    c5_colnorm = 0.0;
    c5_realmin(chartInstance);
    c5_scale = 2.2250738585072014E-308;
    c5_kstart = c5_c_ix0;
    c5_a = c5_kstart;
    c5_kend = c5_a;
    c5_b_kstart = c5_kstart;
    c5_b_kend = c5_kend + 2;
    c5_b_a = c5_b_kstart;
    c5_b = c5_b_kend;
    c5_c_a = c5_b_a;
    c5_b_b = c5_b;
    if (c5_c_a > c5_b_b) {
      c5_overflow = false;
    } else {
      c5_b_eml_switch_helper(chartInstance);
      c5_overflow = (c5_b_b > 2147483646);
    }

    if (c5_overflow) {
      c5_check_forloop_overflow_error(chartInstance, c5_overflow);
    }

    for (c5_k = c5_b_kstart; c5_k <= c5_b_kend; c5_k++) {
      c5_b_k = c5_k;
      c5_x = c5_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c5_b_k), 1, 9, 1, 0) - 1].re;
      c5_b_x = c5_x;
      c5_absxk = muDoubleScalarAbs(c5_b_x);
      if (c5_absxk > c5_scale) {
        c5_t = c5_scale / c5_absxk;
        c5_colnorm = 1.0 + c5_colnorm * c5_t * c5_t;
        c5_scale = c5_absxk;
      } else {
        c5_t = c5_absxk / c5_scale;
        c5_colnorm += c5_t * c5_t;
      }

      c5_c_x = c5_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c5_b_k), 1, 9, 1, 0) - 1].im;
      c5_d_x = c5_c_x;
      c5_absxk = muDoubleScalarAbs(c5_d_x);
      if (c5_absxk > c5_scale) {
        c5_t = c5_scale / c5_absxk;
        c5_colnorm = 1.0 + c5_colnorm * c5_t * c5_t;
        c5_scale = c5_absxk;
      } else {
        c5_t = c5_absxk / c5_scale;
        c5_colnorm += c5_t * c5_t;
      }
    }

    c5_colnorm = c5_scale * muDoubleScalarSqrt(c5_colnorm);
    c5_c_coltop = c5_b_coltop;
    c5_d_a = c5_b_coltop;
    c5_e_a = c5_d_a + 2;
    c5_i99 = c5_e_a;
    c5_f_a = c5_c_coltop;
    c5_c_b = c5_i99;
    c5_g_a = c5_f_a;
    c5_d_b = c5_c_b;
    if (c5_g_a > c5_d_b) {
      c5_b_overflow = false;
    } else {
      c5_b_eml_switch_helper(chartInstance);
      c5_b_overflow = (c5_d_b > 2147483646);
    }

    if (c5_b_overflow) {
      c5_check_forloop_overflow_error(chartInstance, c5_b_overflow);
    }

    for (c5_j = c5_c_coltop; c5_j <= c5_i99; c5_j++) {
      c5_b_j = c5_j;
      c5_d_A.re = c5_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 9, 1, 0) - 1].re;
      c5_d_A.im = c5_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 9, 1, 0) - 1].im;
      c5_B = c5_colnorm;
      c5_y = c5_B;
      c5_d_A = c5_eml_div(chartInstance, c5_d_A, c5_y);
      c5_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c5_b_j), 1, 9, 1, 0) - 1].re = c5_d_A.re;
      c5_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c5_b_j), 1, 9, 1, 0) - 1].im = c5_d_A.im;
    }
  }

  c5_c_info = c5_b_info;
  c5_d_info = c5_c_info;
  for (c5_i100 = 0; c5_i100 < 3; c5_i100++) {
    c5_b_alpha1[c5_i100] = c5_alpha1[c5_i100];
  }

  for (c5_i101 = 0; c5_i101 < 3; c5_i101++) {
    c5_b_beta1[c5_i101] = c5_beta1[c5_i101];
  }

  c5_b_eml_div(chartInstance, c5_b_alpha1, c5_b_beta1, c5_alpha1);
  for (c5_i102 = 0; c5_i102 < 9; c5_i102++) {
    c5_D[c5_i102] = c5_dc5;
  }

  for (c5_c_j = 1; c5_c_j < 4; c5_c_j++) {
    c5_d_j = c5_c_j;
    c5_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_d_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
      1].re = c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 1, 0) - 1].re;
    c5_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_d_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
      1].im = c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 1, 0) - 1].im;
  }

  if (c5_d_info < 0.0) {
    c5_eml_warning(chartInstance);
  } else {
    if (c5_d_info > 0.0) {
      c5_b_eml_warning(chartInstance);
    }
  }
}

static void c5_eml_matlab_zggev(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], real_T *c5_info, creal_T c5_alpha1[3], creal_T c5_beta1[3],
  creal_T c5_V[9])
{
  int32_T c5_i103;
  creal_T c5_b_A[9];
  real_T c5_anrm;
  int32_T c5_i104;
  int32_T c5_i105;
  int32_T c5_i106;
  boolean_T c5_ilascl;
  real_T c5_anrmto;
  int32_T c5_rscale[3];
  int32_T c5_ihi;
  int32_T c5_ilo;
  int32_T c5_b_ilo;
  int32_T c5_b_ihi;
  real_T c5_b_info;
  int32_T c5_i107;
  creal_T c5_c_A[9];
  int32_T c5_c_ilo;
  int32_T c5_c_ihi;
  int32_T c5_a;
  int32_T c5_b_a;
  int32_T c5_i;
  int32_T c5_k;
  int32_T c5_j;
  int32_T c5_b_j;
  creal_T c5_tmp;
  int32_T c5_c_a;
  int32_T c5_d_a;
  int32_T c5_e_a;
  int32_T c5_f_a;
  int32_T c5_i108;
  int32_T c5_g_a;
  int32_T c5_h_a;
  boolean_T c5_overflow;
  int32_T c5_b_i;
  int32_T c5_c_j;
  int32_T c5_jc;
  real_T c5_b_jc;
  real_T c5_x;
  real_T c5_b_x;
  real_T c5_y;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_b_y;
  real_T c5_vtemp;
  int32_T c5_jr;
  real_T c5_b_jr;
  real_T c5_e_x;
  real_T c5_f_x;
  real_T c5_c_y;
  real_T c5_g_x;
  real_T c5_h_x;
  real_T c5_d_y;
  real_T c5_e_y;
  real_T c5_f_y;
  int32_T c5_c_jr;
  real_T c5_b;
  boolean_T guard1 = false;
  *c5_info = 0.0;
  c5_realmin(chartInstance);
  c5_eps(chartInstance);
  for (c5_i103 = 0; c5_i103 < 9; c5_i103++) {
    c5_b_A[c5_i103] = c5_A[c5_i103];
  }

  c5_anrm = c5_eml_matlab_zlangeM(chartInstance, c5_b_A);
  if (!c5_isfinite(chartInstance, c5_anrm)) {
    for (c5_i104 = 0; c5_i104 < 3; c5_i104++) {
      c5_alpha1[c5_i104].re = rtNaN;
      c5_alpha1[c5_i104].im = 0.0;
    }

    for (c5_i105 = 0; c5_i105 < 3; c5_i105++) {
      c5_beta1[c5_i105].re = rtNaN;
      c5_beta1[c5_i105].im = 0.0;
    }

    for (c5_i106 = 0; c5_i106 < 9; c5_i106++) {
      c5_V[c5_i106].re = rtNaN;
      c5_V[c5_i106].im = 0.0;
    }
  } else {
    c5_ilascl = false;
    c5_anrmto = c5_anrm;
    guard1 = false;
    if (c5_anrm > 0.0) {
      if (c5_anrm < 6.7178761075670888E-139) {
        c5_anrmto = 6.7178761075670888E-139;
        c5_ilascl = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      if (c5_anrm > 1.4885657073574029E+138) {
        c5_anrmto = 1.4885657073574029E+138;
        c5_ilascl = true;
      }
    }

    if (c5_ilascl) {
      c5_c_eml_matlab_zlascl(chartInstance, c5_anrm, c5_anrmto, c5_A);
    }

    c5_b_eml_matlab_zggbal(chartInstance, c5_A, &c5_ilo, &c5_ihi, c5_rscale);
    c5_b_ilo = c5_ilo;
    c5_b_ihi = c5_ihi;
    c5_b_eml_matlab_zgghrd(chartInstance, c5_b_ilo, c5_b_ihi, c5_A, c5_V);
    c5_c_eml_matlab_zhgeqz(chartInstance, c5_A, c5_b_ilo, c5_b_ihi, c5_V,
      &c5_b_info, c5_alpha1, c5_beta1);
    *c5_info = c5_b_info;
    if (*c5_info != 0.0) {
    } else {
      for (c5_i107 = 0; c5_i107 < 9; c5_i107++) {
        c5_c_A[c5_i107] = c5_A[c5_i107];
      }

      c5_b_eml_matlab_ztgevc(chartInstance, c5_c_A, c5_V);
      c5_c_ilo = c5_b_ilo;
      c5_c_ihi = c5_b_ihi;
      if (c5_c_ilo > 1) {
        c5_a = c5_c_ilo;
        c5_b_a = c5_a - 1;
        c5_i = c5_b_a;
        while (c5_i >= 1) {
          c5_k = c5_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_i), 1, 3, 1, 0) - 1];
          if (c5_k != c5_i) {
            for (c5_j = 1; c5_j < 4; c5_j++) {
              c5_b_j = c5_j;
              c5_tmp.re = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].
                re;
              c5_tmp.im = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].
                im;
              c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re = c5_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
              c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im = c5_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
              c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re =
                c5_tmp.re;
              c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im =
                c5_tmp.im;
            }
          }

          c5_c_a = c5_i;
          c5_d_a = c5_c_a - 1;
          c5_i = c5_d_a;
        }
      }

      if (c5_c_ihi < 3) {
        c5_e_a = c5_c_ihi;
        c5_f_a = c5_e_a + 1;
        c5_i108 = c5_f_a;
        c5_g_a = c5_i108;
        c5_h_a = c5_g_a;
        if (c5_h_a > 3) {
          c5_overflow = false;
        } else {
          c5_b_eml_switch_helper(chartInstance);
          c5_overflow = false;
        }

        if (c5_overflow) {
          c5_check_forloop_overflow_error(chartInstance, c5_overflow);
        }

        for (c5_b_i = c5_i108; c5_b_i < 4; c5_b_i++) {
          c5_i = c5_b_i;
          c5_k = c5_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_i), 1, 3, 1, 0) - 1];
          if (c5_k != c5_i) {
            for (c5_c_j = 1; c5_c_j < 4; c5_c_j++) {
              c5_b_j = c5_c_j;
              c5_tmp.re = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].
                re;
              c5_tmp.im = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].
                im;
              c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re = c5_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
              c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im = c5_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
              c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re =
                c5_tmp.re;
              c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im =
                c5_tmp.im;
            }
          }
        }
      }

      for (c5_jc = 0; c5_jc < 3; c5_jc++) {
        c5_b_jc = 1.0 + (real_T)c5_jc;
        c5_tmp.re = c5_V[3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1)].re;
        c5_tmp.im = c5_V[3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1)].im;
        c5_x = c5_tmp.re;
        c5_b_x = c5_x;
        c5_y = muDoubleScalarAbs(c5_b_x);
        c5_c_x = c5_tmp.im;
        c5_d_x = c5_c_x;
        c5_b_y = muDoubleScalarAbs(c5_d_x);
        c5_vtemp = c5_y + c5_b_y;
        for (c5_jr = 0; c5_jr < 2; c5_jr++) {
          c5_b_jr = 2.0 + (real_T)c5_jr;
          c5_tmp.re = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1)) - 1].re;
          c5_tmp.im = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1)) - 1].im;
          c5_e_x = c5_tmp.re;
          c5_f_x = c5_e_x;
          c5_c_y = muDoubleScalarAbs(c5_f_x);
          c5_g_x = c5_tmp.im;
          c5_h_x = c5_g_x;
          c5_d_y = muDoubleScalarAbs(c5_h_x);
          c5_e_y = c5_c_y + c5_d_y;
          c5_f_y = c5_e_y;
          if (c5_f_y > c5_vtemp) {
            c5_vtemp = c5_f_y;
          }
        }

        if (c5_vtemp >= 6.7178761075670888E-139) {
          c5_vtemp = 1.0 / c5_vtemp;
          for (c5_c_jr = 0; c5_c_jr < 3; c5_c_jr++) {
            c5_b_jr = 1.0 + (real_T)c5_c_jr;
            c5_tmp.re = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1)) - 1].re;
            c5_tmp.im = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1)) - 1].im;
            c5_b = c5_vtemp;
            c5_tmp.re *= c5_b;
            c5_tmp.im *= c5_b;
            c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    c5_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1))
              - 1].re = c5_tmp.re;
            c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    c5_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1))
              - 1].im = c5_tmp.im;
          }
        }
      }

      if (c5_ilascl) {
        c5_d_eml_matlab_zlascl(chartInstance, c5_anrmto, c5_anrm, c5_alpha1);
      }
    }
  }
}

static void c5_eml_matlab_zgghrd(SFc5_Model_01InstanceStruct *chartInstance,
  int32_T c5_ilo, int32_T c5_ihi, creal_T c5_A[9], creal_T c5_b_A[9], creal_T
  c5_Z[9])
{
  int32_T c5_i109;
  for (c5_i109 = 0; c5_i109 < 9; c5_i109++) {
    c5_b_A[c5_i109] = c5_A[c5_i109];
  }

  c5_b_eml_matlab_zgghrd(chartInstance, c5_ilo, c5_ihi, c5_b_A, c5_Z);
}

static void c5_b_eml_matlab_zhgeqz(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T c5_ilo, int32_T c5_ihi, creal_T c5_Z[9], real_T
  *c5_info, creal_T c5_alpha1[3], creal_T c5_beta1[3], creal_T c5_b_A[9],
  creal_T c5_b_Z[9])
{
  int32_T c5_i110;
  int32_T c5_i111;
  for (c5_i110 = 0; c5_i110 < 9; c5_i110++) {
    c5_b_A[c5_i110] = c5_A[c5_i110];
  }

  for (c5_i111 = 0; c5_i111 < 9; c5_i111++) {
    c5_b_Z[c5_i111] = c5_Z[c5_i111];
  }

  c5_c_eml_matlab_zhgeqz(chartInstance, c5_b_A, c5_ilo, c5_ihi, c5_b_Z, c5_info,
    c5_alpha1, c5_beta1);
}

static void c5_eml_matlab_ztgevc(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], creal_T c5_V[9], creal_T c5_b_V[9])
{
  int32_T c5_i112;
  int32_T c5_i113;
  creal_T c5_b_A[9];
  for (c5_i112 = 0; c5_i112 < 9; c5_i112++) {
    c5_b_V[c5_i112] = c5_V[c5_i112];
  }

  for (c5_i113 = 0; c5_i113 < 9; c5_i113++) {
    c5_b_A[c5_i113] = c5_A[c5_i113];
  }

  c5_b_eml_matlab_ztgevc(chartInstance, c5_b_A, c5_b_V);
}

static creal_T c5_rdivide(SFc5_Model_01InstanceStruct *chartInstance, creal_T
  c5_x, creal_T c5_y)
{
  creal_T c5_z;
  real_T c5_ar;
  real_T c5_ai;
  real_T c5_br;
  real_T c5_bi;
  real_T c5_brm;
  real_T c5_bim;
  real_T c5_s;
  real_T c5_d;
  real_T c5_nr;
  real_T c5_ni;
  real_T c5_sgnbr;
  real_T c5_sgnbi;
  (void)chartInstance;
  c5_ar = c5_x.re;
  c5_ai = c5_x.im;
  c5_br = c5_y.re;
  c5_bi = c5_y.im;
  if (c5_bi == 0.0) {
    if (c5_ai == 0.0) {
      c5_z.re = c5_ar / c5_br;
      c5_z.im = 0.0;
    } else if (c5_ar == 0.0) {
      c5_z.re = 0.0;
      c5_z.im = c5_ai / c5_br;
    } else {
      c5_z.re = c5_ar / c5_br;
      c5_z.im = c5_ai / c5_br;
    }
  } else if (c5_br == 0.0) {
    if (c5_ar == 0.0) {
      c5_z.re = c5_ai / c5_bi;
      c5_z.im = 0.0;
    } else if (c5_ai == 0.0) {
      c5_z.re = 0.0;
      c5_z.im = -(c5_ar / c5_bi);
    } else {
      c5_z.re = c5_ai / c5_bi;
      c5_z.im = -(c5_ar / c5_bi);
    }
  } else {
    c5_brm = muDoubleScalarAbs(c5_br);
    c5_bim = muDoubleScalarAbs(c5_bi);
    if (c5_brm > c5_bim) {
      c5_s = c5_bi / c5_br;
      c5_d = c5_br + c5_s * c5_bi;
      c5_nr = c5_ar + c5_s * c5_ai;
      c5_ni = c5_ai - c5_s * c5_ar;
      c5_z.re = c5_nr / c5_d;
      c5_z.im = c5_ni / c5_d;
    } else if (c5_bim == c5_brm) {
      if (c5_br > 0.0) {
        c5_sgnbr = 0.5;
      } else {
        c5_sgnbr = -0.5;
      }

      if (c5_bi > 0.0) {
        c5_sgnbi = 0.5;
      } else {
        c5_sgnbi = -0.5;
      }

      c5_nr = c5_ar * c5_sgnbr + c5_ai * c5_sgnbi;
      c5_ni = c5_ai * c5_sgnbr - c5_ar * c5_sgnbi;
      c5_z.re = c5_nr / c5_brm;
      c5_z.im = c5_ni / c5_brm;
    } else {
      c5_s = c5_br / c5_bi;
      c5_d = c5_bi + c5_s * c5_br;
      c5_nr = c5_s * c5_ar + c5_ai;
      c5_ni = c5_s * c5_ai - c5_ar;
      c5_z.re = c5_nr / c5_d;
      c5_z.im = c5_ni / c5_d;
    }
  }

  return c5_z;
}

static creal_T c5_power(SFc5_Model_01InstanceStruct *chartInstance, creal_T c5_a,
  real_T c5_b)
{
  creal_T c5_y;
  real_T c5_b_b;
  real_T c5_bk;
  real_T c5_c_b;
  real_T c5_ar;
  real_T c5_ai;
  real_T c5_br;
  real_T c5_bi;
  real_T c5_ytmp;
  real_T c5_x;
  real_T c5_xk;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_r;
  real_T c5_d6;
  int8_T c5_i114;
  int8_T c5_b_r;
  creal_T c5_t;
  real_T c5_e_x;
  boolean_T c5_d_b;
  real_T c5_A;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_h_x;
  real_T c5_b_y;
  real_T c5_b_A;
  real_T c5_i_x;
  real_T c5_j_x;
  real_T c5_k_x;
  real_T c5_c_y;
  real_T c5_l_x;
  real_T c5_d_y;
  real_T c5_x1;
  real_T c5_x2;
  real_T c5_b_a;
  real_T c5_e_b;
  real_T c5_z;
  real_T c5_e_y;
  real_T c5_m_x;
  real_T c5_c_r;
  real_T c5_b_x1;
  real_T c5_b_x2;
  real_T c5_c_a;
  real_T c5_f_b;
  real_T c5_f_y;
  real_T c5_g_y;
  real_T c5_n_x;
  real_T c5_d_r;
  real_T c5_d_a;
  real_T c5_tr;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T guard5 = false;
  boolean_T guard6 = false;
  boolean_T guard7 = false;
  c5_b_b = c5_b;
  c5_scalarEg(chartInstance);
  c5_bk = c5_b_b;
  c5_c_b = c5_bk;
  c5_scalarEg(chartInstance);
  c5_ar = c5_a.re;
  c5_ai = c5_a.im;
  c5_br = c5_c_b;
  c5_bi = 0.0;
  guard1 = false;
  guard2 = false;
  if (c5_ai == 0.0) {
    if (c5_bi == 0.0) {
      if (c5_ar >= 0.0) {
        c5_y.re = muDoubleScalarPower(c5_ar, c5_br);
        c5_y.im = 0.0;
      } else {
        guard1 = true;
      }
    } else {
      guard2 = true;
    }
  } else {
    guard2 = true;
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    guard3 = false;
    guard4 = false;
    guard5 = false;
    if (c5_ar == 0.0) {
      if (c5_bi == 0.0) {
        if (muDoubleScalarFloor(c5_br) == c5_br) {
          if (muDoubleScalarAbs(c5_br) <= 9.007199254740992E+15) {
            c5_ytmp = muDoubleScalarPower(c5_ai, c5_br);
            c5_x = c5_br;
            c5_eml_scalar_eg(chartInstance);
            c5_xk = c5_x;
            c5_b_x = c5_xk;
            c5_eml_scalar_eg(chartInstance);
            c5_c_x = c5_b_x / 4.0;
            c5_d_x = c5_c_x;
            c5_d_x = muDoubleScalarFloor(c5_d_x);
            c5_r = c5_b_x - c5_d_x * 4.0;
            c5_d6 = muDoubleScalarRound(c5_r);
            if (c5_d6 < 128.0) {
              if (c5_d6 >= -128.0) {
                c5_i114 = (int8_T)c5_d6;
              } else {
                c5_i114 = MIN_int8_T;
              }
            } else if (c5_d6 >= 128.0) {
              c5_i114 = MAX_int8_T;
            } else {
              c5_i114 = 0;
            }

            c5_b_r = c5_i114;
            if ((real_T)c5_b_r == 3.0) {
              c5_y.re = 0.0;
              c5_y.im = -c5_ytmp;
            } else if ((real_T)c5_b_r == 2.0) {
              c5_y.re = -c5_ytmp;
              c5_y.im = 0.0;
            } else if ((real_T)c5_b_r == 1.0) {
              c5_y.re = 0.0;
              c5_y.im = c5_ytmp;
            } else {
              c5_y.re = c5_ytmp;
              c5_y.im = 0.0;
            }
          } else {
            guard3 = true;
          }
        } else {
          guard4 = true;
        }
      } else {
        guard5 = true;
      }
    } else {
      guard5 = true;
    }

    if (guard5 == true) {
      guard4 = true;
    }

    if (guard4 == true) {
      guard3 = true;
    }

    if (guard3 == true) {
      c5_t = c5_a;
      c5_realmax(chartInstance);
      guard6 = false;
      if (c5_t.im == 0.0) {
        c5_e_x = c5_t.re;
        c5_d_b = muDoubleScalarIsNaN(c5_e_x);
        if (c5_d_b) {
        } else {
          guard6 = true;
        }
      } else {
        guard6 = true;
      }

      if (guard6 == true) {
        guard7 = false;
        if (muDoubleScalarAbs(c5_t.re) > 8.9884656743115785E+307) {
          guard7 = true;
        } else if (muDoubleScalarAbs(c5_t.im) > 8.9884656743115785E+307) {
          guard7 = true;
        } else {
          c5_b_x1 = c5_t.re;
          c5_b_x2 = c5_t.im;
          c5_c_a = c5_b_x1;
          c5_f_b = c5_b_x2;
          c5_f_y = muDoubleScalarHypot(c5_c_a, c5_f_b);
          c5_g_y = c5_t.im;
          c5_n_x = c5_t.re;
          c5_d_r = muDoubleScalarAtan2(c5_g_y, c5_n_x);
          c5_t.re = muDoubleScalarLog(c5_f_y);
          c5_t.im = c5_d_r;
        }

        if (guard7 == true) {
          c5_A = c5_t.re;
          c5_f_x = c5_A;
          c5_g_x = c5_f_x;
          c5_h_x = c5_g_x;
          c5_b_y = c5_h_x / 2.0;
          c5_b_A = c5_t.im;
          c5_i_x = c5_b_A;
          c5_j_x = c5_i_x;
          c5_k_x = c5_j_x;
          c5_c_y = c5_k_x / 2.0;
          c5_l_x = c5_b_y;
          c5_d_y = c5_c_y;
          c5_eml_scalar_eg(chartInstance);
          c5_x1 = c5_l_x;
          c5_x2 = c5_d_y;
          c5_b_a = c5_x1;
          c5_e_b = c5_x2;
          c5_z = muDoubleScalarHypot(c5_b_a, c5_e_b);
          c5_e_y = c5_t.im;
          c5_m_x = c5_t.re;
          c5_c_r = muDoubleScalarAtan2(c5_e_y, c5_m_x);
          c5_t.re = muDoubleScalarLog(c5_z) + 0.69314718055994529;
          c5_t.im = c5_c_r;
        }
      }

      c5_d_a = c5_c_b;
      c5_t.re *= c5_d_a;
      c5_t.im *= c5_d_a;
      c5_tr = muDoubleScalarExp(c5_t.re);
      c5_y.re = c5_tr * muDoubleScalarCos(c5_t.im);
      c5_y.im = c5_tr * muDoubleScalarSin(c5_t.im);
    }
  }

  return c5_y;
}

static void c5_eml_lusolve(SFc5_Model_01InstanceStruct *chartInstance, creal_T
  c5_A[9], creal_T c5_B[9], creal_T c5_X[9])
{
  int32_T c5_i115;
  creal_T c5_b_A[9];
  int32_T c5_r1;
  int32_T c5_r2;
  int32_T c5_r3;
  creal_T c5_x;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_y;
  real_T c5_d_x;
  real_T c5_e_x;
  real_T c5_b_y;
  real_T c5_maxval;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_c_y;
  real_T c5_h_x;
  real_T c5_i_x;
  real_T c5_d_y;
  real_T c5_a21;
  real_T c5_j_x;
  real_T c5_k_x;
  real_T c5_e_y;
  real_T c5_l_x;
  real_T c5_m_x;
  real_T c5_f_y;
  real_T c5_d;
  creal_T c5_c_A;
  creal_T c5_d_A;
  creal_T c5_e_A;
  creal_T c5_f_A;
  creal_T c5_g_A;
  creal_T c5_h_A;
  creal_T c5_i_A;
  creal_T c5_j_A;
  real_T c5_n_x;
  real_T c5_o_x;
  real_T c5_g_y;
  real_T c5_p_x;
  real_T c5_q_x;
  real_T c5_h_y;
  real_T c5_b_d;
  real_T c5_r_x;
  real_T c5_s_x;
  real_T c5_i_y;
  real_T c5_t_x;
  real_T c5_u_x;
  real_T c5_j_y;
  real_T c5_c_d;
  int32_T c5_rtemp;
  creal_T c5_k_A;
  creal_T c5_l_A;
  creal_T c5_m_A;
  static creal_T c5_dc6 = { 0.0, 0.0 };

  boolean_T c5_n_A;
  boolean_T c5_o_A;
  boolean_T c5_p_A;
  int32_T c5_k;
  int32_T c5_b_k;
  creal_T c5_b_B;
  creal_T c5_q_A;
  creal_T c5_b_X;
  creal_T c5_c_X;
  creal_T c5_d_X;
  creal_T c5_r_A;
  creal_T c5_e_X;
  creal_T c5_f_X;
  creal_T c5_s_A;
  creal_T c5_g_X;
  creal_T c5_h_X;
  creal_T c5_i_X;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  for (c5_i115 = 0; c5_i115 < 9; c5_i115++) {
    c5_b_A[c5_i115] = c5_A[c5_i115];
  }

  c5_r1 = 1;
  c5_r2 = 2;
  c5_r3 = 3;
  c5_x = c5_b_A[0];
  c5_b_x = c5_x.re;
  c5_c_x = c5_b_x;
  c5_y = muDoubleScalarAbs(c5_c_x);
  c5_d_x = c5_x.im;
  c5_e_x = c5_d_x;
  c5_b_y = muDoubleScalarAbs(c5_e_x);
  c5_maxval = c5_y + c5_b_y;
  c5_x = c5_b_A[1];
  c5_f_x = c5_x.re;
  c5_g_x = c5_f_x;
  c5_c_y = muDoubleScalarAbs(c5_g_x);
  c5_h_x = c5_x.im;
  c5_i_x = c5_h_x;
  c5_d_y = muDoubleScalarAbs(c5_i_x);
  c5_a21 = c5_c_y + c5_d_y;
  if (c5_a21 > c5_maxval) {
    c5_maxval = c5_a21;
    c5_r1 = 2;
    c5_r2 = 1;
  }

  c5_x = c5_b_A[2];
  c5_j_x = c5_x.re;
  c5_k_x = c5_j_x;
  c5_e_y = muDoubleScalarAbs(c5_k_x);
  c5_l_x = c5_x.im;
  c5_m_x = c5_l_x;
  c5_f_y = muDoubleScalarAbs(c5_m_x);
  c5_d = c5_e_y + c5_f_y;
  if (c5_d > c5_maxval) {
    c5_r1 = 3;
    c5_r2 = 2;
    c5_r3 = 1;
  }

  c5_c_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r2), 1, 3, 1, 0) - 1].re;
  c5_c_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r2), 1, 3, 1, 0) - 1].im;
  c5_d_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r1), 1, 3, 1, 0) - 1].re;
  c5_d_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r1), 1, 3, 1, 0) - 1].im;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r2), 1, 3, 1, 0) - 1] = c5_rdivide(chartInstance, c5_c_A, c5_d_A);
  c5_e_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) - 1].re;
  c5_e_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) - 1].im;
  c5_f_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r1), 1, 3, 1, 0) - 1].re;
  c5_f_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r1), 1, 3, 1, 0) - 1].im;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r3), 1, 3, 1, 0) - 1] = c5_rdivide(chartInstance, c5_e_A, c5_f_A);
  c5_g_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r2), 1, 3, 1, 0) - 1].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 2].re - c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) - 1].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 2].im;
  c5_g_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r2), 1, 3, 1, 0) - 1].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 2].im + c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) - 1].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 2].re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r2), 1, 3, 1, 0) + 2].re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 2].re -
    c5_g_A.re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r2), 1, 3, 1, 0) + 2].im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 2].im -
    c5_g_A.im;
  c5_h_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) - 1].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 2].re - c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) - 1].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 2].im;
  c5_h_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) - 1].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 2].im + c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) - 1].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 2].re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r3), 1, 3, 1, 0) + 2].re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 2].re -
    c5_h_A.re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r3), 1, 3, 1, 0) + 2].im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 2].im -
    c5_h_A.im;
  c5_i_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r2), 1, 3, 1, 0) - 1].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 5].re - c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) - 1].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 5].im;
  c5_i_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r2), 1, 3, 1, 0) - 1].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 5].im + c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) - 1].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 5].re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r2), 1, 3, 1, 0) + 5].re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 5].re -
    c5_i_A.re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r2), 1, 3, 1, 0) + 5].im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 5].im -
    c5_i_A.im;
  c5_j_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) - 1].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 5].re - c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) - 1].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 5].im;
  c5_j_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) - 1].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 5].im + c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) - 1].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r1), 1, 3, 1, 0) + 5].re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r3), 1, 3, 1, 0) + 5].re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 5].re -
    c5_j_A.re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r3), 1, 3, 1, 0) + 5].im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 5].im -
    c5_j_A.im;
  c5_x.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c5_r3), 1, 3, 1, 0) + 2].re;
  c5_x.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c5_r3), 1, 3, 1, 0) + 2].im;
  c5_n_x = c5_x.re;
  c5_o_x = c5_n_x;
  c5_g_y = muDoubleScalarAbs(c5_o_x);
  c5_p_x = c5_x.im;
  c5_q_x = c5_p_x;
  c5_h_y = muDoubleScalarAbs(c5_q_x);
  c5_b_d = c5_g_y + c5_h_y;
  c5_x.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c5_r2), 1, 3, 1, 0) + 2].re;
  c5_x.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c5_r2), 1, 3, 1, 0) + 2].im;
  c5_r_x = c5_x.re;
  c5_s_x = c5_r_x;
  c5_i_y = muDoubleScalarAbs(c5_s_x);
  c5_t_x = c5_x.im;
  c5_u_x = c5_t_x;
  c5_j_y = muDoubleScalarAbs(c5_u_x);
  c5_c_d = c5_i_y + c5_j_y;
  if (c5_b_d > c5_c_d) {
    c5_rtemp = c5_r2;
    c5_r2 = c5_r3;
    c5_r3 = c5_rtemp;
  }

  c5_k_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) + 2].re;
  c5_k_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) + 2].im;
  c5_l_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r2), 1, 3, 1, 0) + 2].re;
  c5_l_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r2), 1, 3, 1, 0) + 2].im;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r3), 1, 3, 1, 0) + 2] = c5_rdivide(chartInstance, c5_k_A, c5_l_A);
  c5_m_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) + 2].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r2), 1, 3, 1, 0) + 5].re - c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 2].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r2), 1, 3, 1, 0) + 5].im;
  c5_m_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c5_r3), 1, 3, 1, 0) + 2].re *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r2), 1, 3, 1, 0) + 5].im + c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 2].im *
    c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c5_r2), 1, 3, 1, 0) + 5].re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r3), 1, 3, 1, 0) + 5].re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 5].re -
    c5_m_A.re;
  c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c5_r3), 1, 3, 1, 0) + 5].im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 5].im -
    c5_m_A.im;
  c5_n_A = ((c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c5_r1), 1, 3, 1, 0) - 1].re == c5_dc6.re) &&
            (c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c5_r1), 1, 3, 1, 0) - 1].im == c5_dc6.im));
  guard1 = false;
  guard2 = false;
  if (c5_n_A) {
    guard2 = true;
  } else {
    c5_o_A = ((c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 2].re ==
               c5_dc6.re) && (c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 2].im ==
               c5_dc6.im));
    if (c5_o_A) {
      guard2 = true;
    } else {
      c5_p_A = ((c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 5].re ==
                 c5_dc6.re) && (c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 5].im ==
                 c5_dc6.im));
      if (c5_p_A) {
        guard1 = true;
      }
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c5_c_eml_warning(chartInstance);
  }

  for (c5_k = 1; c5_k < 4; c5_k++) {
    c5_b_k = c5_k;
    c5_b_B.re = c5_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", (real_T)c5_b_k), 1, 3, 1, 0) - 1].re;
    c5_b_B.im = c5_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", (real_T)c5_b_k), 1, 3, 1, 0) - 1].im;
    c5_q_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 1, 0) - 1].re;
    c5_q_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 1, 0) - 1].im;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) -
      1] = c5_rdivide(chartInstance, c5_b_B, c5_q_A);
    c5_b_X.re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r1), 1, 3, 1, 0) + 2].re - c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r1), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 1, 0) + 2].im;
    c5_b_X.im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r1), 1, 3, 1, 0) + 2].im + c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r1), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 1, 0) + 2].re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) -
      1].re = c5_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 2].re - c5_b_X.re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) -
      1].im = c5_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 2].im - c5_b_X.im;
    c5_c_X.re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r1), 1, 3, 1, 0) + 5].re - c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r1), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 1, 0) + 5].im;
    c5_c_X.im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r1), 1, 3, 1, 0) + 5].im + c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r1), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 1, 0) + 5].re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) -
      1].re = c5_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 5].re - c5_c_X.re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) -
      1].im = c5_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 5].im - c5_c_X.im;
    c5_d_X.re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) - 1].re;
    c5_d_X.im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) - 1].im;
    c5_r_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 2].re;
    c5_r_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 2].im;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) -
      1] = c5_rdivide(chartInstance, c5_d_X, c5_r_A);
    c5_e_X.re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r2), 1, 3, 1, 0) + 5].re - c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r2), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 5].im;
    c5_e_X.im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r2), 1, 3, 1, 0) + 5].im + c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r2), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) + 5].re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) -
      1].re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) - 1].re
      - c5_e_X.re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) -
      1].im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) - 1].im
      - c5_e_X.im;
    c5_f_X.re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) - 1].re;
    c5_f_X.im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) - 1].im;
    c5_s_A.re = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 5].re;
    c5_s_A.im = c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 5].im;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) -
      1] = c5_rdivide(chartInstance, c5_f_X, c5_s_A);
    c5_g_X.re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r3), 1, 3, 1, 0) + 2].re - c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r3), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 2].im;
    c5_g_X.im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r3), 1, 3, 1, 0) + 2].im + c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r3), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) + 2].re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) -
      1].re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) - 1].re
      - c5_g_X.re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) -
      1].im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) - 1].im
      - c5_g_X.im;
    c5_h_X.re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r3), 1, 3, 1, 0) - 1].re - c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r3), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) - 1].im;
    c5_h_X.im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r3), 1, 3, 1, 0) - 1].im + c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r3), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r3), 1, 3, 1, 0) - 1].re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) -
      1].re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) - 1].re
      - c5_h_X.re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) -
      1].im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) - 1].im
      - c5_h_X.im;
    c5_i_X.re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r2), 1, 3, 1, 0) - 1].re - c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r2), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) - 1].im;
    c5_i_X.im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_r2), 1, 3, 1, 0) - 1].im + c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_r2), 1, 3, 2, 0) - 1)) - 1].im * c5_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r2), 1, 3, 1, 0) - 1].re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) -
      1].re = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) - 1].re
      - c5_i_X.re;
    c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) -
      1].im = c5_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c5_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_r1), 1, 3, 2, 0) - 1)) - 1].im
      - c5_i_X.im;
  }
}

static void c5_e_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_sprintf, const char_T *c5_identifier, char_T c5_y[14])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_sprintf), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_sprintf);
}

static void c5_f_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, char_T c5_y[14])
{
  char_T c5_cv6[14];
  int32_T c5_i116;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_cv6, 1, 10, 0U, 1, 0U, 2, 1,
                14);
  for (c5_i116 = 0; c5_i116 < 14; c5_i116++) {
    c5_y[c5_i116] = c5_cv6[c5_i116];
  }

  sf_mex_destroy(&c5_u);
}

static const mxArray *c5_d_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_u;
  const mxArray *c5_y = NULL;
  SFc5_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc5_Model_01InstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_u = *(int32_T *)c5_inData;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, false);
  return c5_mxArrayOutData;
}

static int32_T c5_g_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  int32_T c5_y;
  int32_T c5_i117;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_i117, 1, 6, 0U, 0, 0U, 0);
  c5_y = c5_i117;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_b_sfEvent;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  int32_T c5_y;
  SFc5_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc5_Model_01InstanceStruct *)chartInstanceVoid;
  c5_b_sfEvent = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_sfEvent),
    &c5_thisId);
  sf_mex_destroy(&c5_b_sfEvent);
  *(int32_T *)c5_outData = c5_y;
  sf_mex_destroy(&c5_mxArrayInData);
}

static uint8_T c5_h_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_b_is_active_c5_Model_01, const char_T *c5_identifier)
{
  uint8_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_i_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c5_b_is_active_c5_Model_01), &c5_thisId);
  sf_mex_destroy(&c5_b_is_active_c5_Model_01);
  return c5_y;
}

static uint8_T c5_i_emlrt_marshallIn(SFc5_Model_01InstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  uint8_T c5_y;
  uint8_T c5_u0;
  (void)chartInstance;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_u0, 1, 3, 0U, 0, 0U, 0);
  c5_y = c5_u0;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_c_eml_matlab_zlascl(SFc5_Model_01InstanceStruct *chartInstance,
  real_T c5_cfrom, real_T c5_cto, creal_T c5_A[9])
{
  real_T c5_cfromc;
  real_T c5_ctoc;
  boolean_T c5_notdone;
  real_T c5_cfrom1;
  real_T c5_cto1;
  real_T c5_x;
  real_T c5_b_x;
  real_T c5_y;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_b_y;
  real_T c5_mul;
  real_T c5_e_x;
  real_T c5_f_x;
  real_T c5_c_y;
  real_T c5_g_x;
  real_T c5_h_x;
  real_T c5_d_y;
  real_T c5_a;
  int32_T c5_i118;
  int32_T c5_i119;
  int32_T c5_i120;
  boolean_T guard1 = false;
  c5_realmin(chartInstance);
  c5_eps(chartInstance);
  c5_cfromc = c5_cfrom;
  c5_ctoc = c5_cto;
  c5_notdone = true;
  while (c5_notdone) {
    c5_cfrom1 = c5_cfromc * 2.0041683600089728E-292;
    c5_cto1 = c5_ctoc / 4.9896007738368E+291;
    c5_x = c5_cfrom1;
    c5_b_x = c5_x;
    c5_y = muDoubleScalarAbs(c5_b_x);
    c5_c_x = c5_ctoc;
    c5_d_x = c5_c_x;
    c5_b_y = muDoubleScalarAbs(c5_d_x);
    guard1 = false;
    if (c5_y > c5_b_y) {
      if (c5_ctoc != 0.0) {
        c5_mul = 2.0041683600089728E-292;
        c5_notdone = true;
        c5_cfromc = c5_cfrom1;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      c5_e_x = c5_cto1;
      c5_f_x = c5_e_x;
      c5_c_y = muDoubleScalarAbs(c5_f_x);
      c5_g_x = c5_cfromc;
      c5_h_x = c5_g_x;
      c5_d_y = muDoubleScalarAbs(c5_h_x);
      if (c5_c_y > c5_d_y) {
        c5_mul = 4.9896007738368E+291;
        c5_notdone = true;
        c5_ctoc = c5_cto1;
      } else {
        c5_mul = c5_ctoc / c5_cfromc;
        c5_notdone = false;
      }
    }

    c5_a = c5_mul;
    c5_i118 = 0;
    for (c5_i119 = 0; c5_i119 < 3; c5_i119++) {
      for (c5_i120 = 0; c5_i120 < 3; c5_i120++) {
        c5_A[c5_i120 + c5_i118].re *= c5_a;
        c5_A[c5_i120 + c5_i118].im *= c5_a;
      }

      c5_i118 += 3;
    }
  }
}

static void c5_b_eml_matlab_zggbal(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T *c5_ilo, int32_T *c5_ihi, int32_T c5_rscale[3])
{
  int32_T c5_i121;
  int32_T c5_b_ihi;
  int32_T c5_i;
  int32_T c5_j;
  boolean_T c5_found;
  int32_T c5_ii;
  real_T c5_nzcount;
  int32_T c5_c_ihi;
  int32_T c5_b;
  int32_T c5_b_b;
  boolean_T c5_overflow;
  int32_T c5_jj;
  int32_T c5_b_jj;
  static creal_T c5_dc7 = { 0.0, 0.0 };

  boolean_T c5_b_A;
  int32_T c5_a;
  int32_T c5_b_a;
  int32_T c5_b_i;
  int32_T c5_b_j;
  boolean_T c5_b_found;
  int32_T c5_b_ilo;
  int32_T c5_d_ihi;
  int32_T c5_c_i;
  int32_T c5_c_j;
  boolean_T c5_c_found;
  int32_T c5_c_ilo;
  int32_T c5_e_ihi;
  int32_T c5_c_a;
  int32_T c5_c_b;
  int32_T c5_d_a;
  int32_T c5_d_b;
  boolean_T c5_b_overflow;
  int32_T c5_c_jj;
  int32_T c5_d_jj;
  real_T c5_b_nzcount;
  int32_T c5_d_ilo;
  int32_T c5_f_ihi;
  int32_T c5_e_a;
  int32_T c5_e_b;
  int32_T c5_f_a;
  int32_T c5_f_b;
  boolean_T c5_c_overflow;
  int32_T c5_b_ii;
  int32_T c5_c_ii;
  boolean_T c5_c_A;
  int32_T c5_m;
  int32_T c5_d_i;
  int32_T c5_d_j;
  int32_T c5_e_ilo;
  int32_T c5_g_ihi;
  int32_T c5_f_ilo;
  int32_T c5_g_a;
  int32_T c5_h_a;
  boolean_T c5_d_overflow;
  int32_T c5_k;
  int32_T c5_b_k;
  creal_T c5_atmp;
  int32_T c5_h_ihi;
  int32_T c5_g_b;
  int32_T c5_h_b;
  boolean_T c5_e_overflow;
  int32_T c5_c_k;
  int32_T c5_i_a;
  int32_T c5_j_a;
  int32_T c5_b_m;
  int32_T c5_e_i;
  int32_T c5_e_j;
  int32_T c5_i_ihi;
  int32_T c5_d_k;
  int32_T c5_e_k;
  int32_T c5_j_ihi;
  int32_T c5_i_b;
  int32_T c5_j_b;
  boolean_T c5_f_overflow;
  int32_T c5_f_k;
  int32_T c5_k_a;
  int32_T c5_l_a;
  int32_T exitg1;
  int32_T exitg2;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T exitg5;
  boolean_T exitg6;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  for (c5_i121 = 0; c5_i121 < 3; c5_i121++) {
    c5_rscale[c5_i121] = 0;
  }

  *c5_ilo = 1;
  *c5_ihi = 3;
  do {
    exitg2 = 0;
    c5_b_ihi = *c5_ihi;
    c5_i = 0;
    c5_j = 0;
    c5_found = false;
    c5_ii = c5_b_ihi;
    exitg5 = false;
    while ((exitg5 == false) && (c5_ii > 0)) {
      c5_nzcount = 0.0;
      c5_i = c5_ii;
      c5_j = c5_b_ihi;
      c5_c_ihi = c5_b_ihi;
      c5_b = c5_c_ihi;
      c5_b_b = c5_b;
      if (1 > c5_b_b) {
        c5_overflow = false;
      } else {
        c5_b_eml_switch_helper(chartInstance);
        c5_overflow = (c5_b_b > 2147483646);
      }

      if (c5_overflow) {
        c5_check_forloop_overflow_error(chartInstance, c5_overflow);
      }

      c5_jj = 1;
      exitg6 = false;
      while ((exitg6 == false) && (c5_jj <= c5_c_ihi)) {
        c5_b_jj = c5_jj;
        c5_b_A = ((c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_ii), 1, 3, 1, 0) + 3 *
                         (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_b_jj), 1, 3, 2, 0) - 1)) - 1].re !=
                   c5_dc7.re) || (c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_ii), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_jj), 1, 3, 2, 0) - 1)) - 1].im !=
                   c5_dc7.im));
        guard3 = false;
        guard4 = false;
        if (c5_b_A) {
          guard4 = true;
        } else if (c5_ii == c5_b_jj) {
          guard4 = true;
        } else {
          guard3 = true;
        }

        if (guard4 == true) {
          if (c5_nzcount == 0.0) {
            c5_j = c5_b_jj;
            c5_nzcount = 1.0;
            guard3 = true;
          } else {
            c5_nzcount = 2.0;
            exitg6 = true;
          }
        }

        if (guard3 == true) {
          c5_jj++;
        }
      }

      if (c5_nzcount < 2.0) {
        c5_found = true;
        exitg5 = true;
      } else {
        c5_a = c5_ii;
        c5_b_a = c5_a - 1;
        c5_ii = c5_b_a;
      }
    }

    c5_b_i = c5_i;
    c5_b_j = c5_j;
    c5_b_found = c5_found;
    if (!c5_b_found) {
      exitg2 = 2;
    } else {
      c5_b_m = *c5_ihi;
      c5_e_i = c5_b_i;
      c5_e_j = c5_b_j;
      c5_i_ihi = *c5_ihi;
      if (c5_e_i != c5_b_m) {
        c5_b_eml_switch_helper(chartInstance);
        for (c5_d_k = 1; c5_d_k < 4; c5_d_k++) {
          c5_e_k = c5_d_k;
          c5_atmp.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_e_i), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_e_k), 1, 3, 2, 0) - 1)) - 1].re;
          c5_atmp.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_e_i), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_e_k), 1, 3, 2, 0) - 1)) - 1].im;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_e_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_e_k), 1, 3, 2, 0) - 1)) - 1].re = c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_b_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_e_k), 1, 3, 2, 0)
               - 1)) - 1].re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_e_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_e_k), 1, 3, 2, 0) - 1)) - 1].im = c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_b_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_e_k), 1, 3, 2, 0)
               - 1)) - 1].im;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_m), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_e_k), 1, 3, 2, 0) - 1)) - 1].re = c5_atmp.re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_m), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_e_k), 1, 3, 2, 0) - 1)) - 1].im = c5_atmp.im;
        }
      }

      if (c5_e_j != c5_b_m) {
        c5_j_ihi = c5_i_ihi;
        c5_i_b = c5_j_ihi;
        c5_j_b = c5_i_b;
        if (1 > c5_j_b) {
          c5_f_overflow = false;
        } else {
          c5_b_eml_switch_helper(chartInstance);
          c5_f_overflow = (c5_j_b > 2147483646);
        }

        if (c5_f_overflow) {
          c5_check_forloop_overflow_error(chartInstance, c5_f_overflow);
        }

        for (c5_f_k = 1; c5_f_k <= c5_j_ihi; c5_f_k++) {
          c5_e_k = c5_f_k;
          c5_atmp.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_e_k), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_e_j), 1, 3, 2, 0) - 1)) - 1].re;
          c5_atmp.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_e_k), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_e_j), 1, 3, 2, 0) - 1)) - 1].im;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_e_j), 1, 3, 2, 0) - 1)) - 1].re = c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_e_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_m), 1, 3, 2, 0)
               - 1)) - 1].re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_e_j), 1, 3, 2, 0) - 1)) - 1].im = c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_e_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_m), 1, 3, 2, 0)
               - 1)) - 1].im;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_m), 1, 3, 2, 0) - 1)) - 1].re = c5_atmp.re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_m), 1, 3, 2, 0) - 1)) - 1].im = c5_atmp.im;
        }
      }

      c5_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)*c5_ihi), 1, 3, 1, 0) - 1] = c5_b_j;
      c5_k_a = *c5_ihi;
      c5_l_a = c5_k_a - 1;
      *c5_ihi = c5_l_a;
      if (*c5_ihi == 1) {
        c5_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)*c5_ihi), 1, 3, 1, 0) - 1] = *c5_ihi;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    do {
      exitg1 = 0;
      c5_b_ilo = *c5_ilo;
      c5_d_ihi = *c5_ihi;
      c5_c_i = 0;
      c5_c_j = 0;
      c5_c_found = false;
      c5_c_ilo = c5_b_ilo;
      c5_e_ihi = c5_d_ihi;
      c5_c_a = c5_c_ilo;
      c5_c_b = c5_e_ihi;
      c5_d_a = c5_c_a;
      c5_d_b = c5_c_b;
      if (c5_d_a > c5_d_b) {
        c5_b_overflow = false;
      } else {
        c5_b_eml_switch_helper(chartInstance);
        c5_b_overflow = (c5_d_b > 2147483646);
      }

      if (c5_b_overflow) {
        c5_check_forloop_overflow_error(chartInstance, c5_b_overflow);
      }

      c5_c_jj = c5_c_ilo;
      exitg3 = false;
      while ((exitg3 == false) && (c5_c_jj <= c5_e_ihi)) {
        c5_d_jj = c5_c_jj;
        c5_b_nzcount = 0.0;
        c5_c_i = c5_d_ihi;
        c5_c_j = c5_d_jj;
        c5_d_ilo = c5_b_ilo;
        c5_f_ihi = c5_d_ihi;
        c5_e_a = c5_d_ilo;
        c5_e_b = c5_f_ihi;
        c5_f_a = c5_e_a;
        c5_f_b = c5_e_b;
        if (c5_f_a > c5_f_b) {
          c5_c_overflow = false;
        } else {
          c5_b_eml_switch_helper(chartInstance);
          c5_c_overflow = (c5_f_b > 2147483646);
        }

        if (c5_c_overflow) {
          c5_check_forloop_overflow_error(chartInstance, c5_c_overflow);
        }

        c5_b_ii = c5_d_ilo;
        exitg4 = false;
        while ((exitg4 == false) && (c5_b_ii <= c5_f_ihi)) {
          c5_c_ii = c5_b_ii;
          c5_c_A = ((c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_c_ii), 1, 3, 1, 0) + 3 *
                           (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_jj), 1, 3, 2, 0) - 1)) - 1].re
                     != c5_dc7.re) || (c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_c_ii), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_jj), 1, 3, 2, 0) - 1)) - 1].im
                     != c5_dc7.im));
          guard1 = false;
          guard2 = false;
          if (c5_c_A) {
            guard2 = true;
          } else if (c5_c_ii == c5_d_jj) {
            guard2 = true;
          } else {
            guard1 = true;
          }

          if (guard2 == true) {
            if (c5_b_nzcount == 0.0) {
              c5_c_i = c5_c_ii;
              c5_b_nzcount = 1.0;
              guard1 = true;
            } else {
              c5_b_nzcount = 2.0;
              exitg4 = true;
            }
          }

          if (guard1 == true) {
            c5_b_ii++;
          }
        }

        if (c5_b_nzcount < 2.0) {
          c5_c_found = true;
          exitg3 = true;
        } else {
          c5_c_jj++;
        }
      }

      c5_b_i = c5_c_i;
      c5_b_j = c5_c_j;
      c5_b_found = c5_c_found;
      if (!c5_b_found) {
        exitg1 = 1;
      } else {
        c5_m = *c5_ilo;
        c5_d_i = c5_b_i;
        c5_d_j = c5_b_j;
        c5_e_ilo = *c5_ilo;
        c5_g_ihi = *c5_ihi;
        if (c5_d_i != c5_m) {
          c5_f_ilo = c5_e_ilo;
          c5_g_a = c5_f_ilo;
          c5_h_a = c5_g_a;
          if (c5_h_a > 3) {
            c5_d_overflow = false;
          } else {
            c5_b_eml_switch_helper(chartInstance);
            c5_d_overflow = false;
          }

          if (c5_d_overflow) {
            c5_check_forloop_overflow_error(chartInstance, c5_d_overflow);
          }

          for (c5_k = c5_f_ilo; c5_k < 4; c5_k++) {
            c5_b_k = c5_k;
            c5_atmp.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 2, 0) - 1)) - 1].re;
            c5_atmp.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 2, 0) - 1)) - 1].im;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_b_k), 1, 3, 2, 0) - 1)) - 1].re = c5_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                  "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 2,
                  0) - 1)) - 1].re;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_b_k), 1, 3, 2, 0) - 1)) - 1].im = c5_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                  "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 2,
                  0) - 1)) - 1].im;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_m), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_b_k), 1, 3, 2, 0) - 1)) - 1].re = c5_atmp.re;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_m), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_b_k), 1, 3, 2, 0) - 1)) - 1].im = c5_atmp.im;
          }
        }

        if (c5_d_j != c5_m) {
          c5_h_ihi = c5_g_ihi;
          c5_g_b = c5_h_ihi;
          c5_h_b = c5_g_b;
          if (1 > c5_h_b) {
            c5_e_overflow = false;
          } else {
            c5_b_eml_switch_helper(chartInstance);
            c5_e_overflow = (c5_h_b > 2147483646);
          }

          if (c5_e_overflow) {
            c5_check_forloop_overflow_error(chartInstance, c5_e_overflow);
          }

          for (c5_c_k = 1; c5_c_k <= c5_h_ihi; c5_c_k++) {
            c5_b_k = c5_c_k;
            c5_atmp.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) - 1].re;
            c5_atmp.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) - 1].im;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) - 1].re = c5_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_m), 1, 3, 2, 0) - 1)) - 1].re;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) - 1].im = c5_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_m), 1, 3, 2, 0) - 1)) - 1].im;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_m), 1, 3, 2, 0) - 1)) - 1].re = c5_atmp.re;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_m), 1, 3, 2, 0) - 1)) - 1].im = c5_atmp.im;
          }
        }

        c5_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)*c5_ilo), 1, 3, 1, 0) - 1] = c5_b_j;
        c5_i_a = *c5_ilo;
        c5_j_a = c5_i_a + 1;
        *c5_ilo = c5_j_a;
        if (*c5_ilo == *c5_ihi) {
          c5_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
            "", (real_T)*c5_ilo), 1, 3, 1, 0) - 1] = *c5_ilo;
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

static void c5_b_sqrt(SFc5_Model_01InstanceStruct *chartInstance, creal_T *c5_x)
{
  real_T c5_yr;
  real_T c5_yi;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_z;
  real_T c5_e_x;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_b_z;
  boolean_T c5_b3;
  boolean_T c5_b4;
  boolean_T c5_b;
  real_T c5_h_x;
  boolean_T c5_b_b;
  real_T c5_i_x;
  boolean_T c5_c_b;
  real_T c5_absxr;
  real_T c5_absxi;
  real_T c5_j_x;
  real_T c5_y;
  real_T c5_x1;
  real_T c5_x2;
  real_T c5_a;
  real_T c5_d_b;
  real_T c5_absxd2;
  real_T c5_k_x;
  real_T c5_b_y;
  real_T c5_l_x;
  real_T c5_c_y;
  real_T c5_m_x;
  real_T c5_d_y;
  real_T c5_c_z;
  real_T c5_n_x;
  real_T c5_e_y;
  real_T c5_b_x1;
  real_T c5_b_x2;
  real_T c5_b_a;
  real_T c5_e_b;
  real_T c5_d_z;
  real_T c5_o_x;
  real_T c5_f_y;
  real_T c5_p_x;
  real_T c5_g_y;
  real_T c5_q_x;
  real_T c5_h_y;
  real_T c5_e_z;
  real_T c5_r_x;
  real_T c5_i_y;
  real_T c5_s_x;
  real_T c5_j_y;
  real_T c5_t_x;
  real_T c5_k_y;
  real_T c5_f_z;
  boolean_T guard1 = false;
  if (c5_x->im == 0.0) {
    if (c5_x->re < 0.0) {
      c5_yr = 0.0;
      c5_yi = muDoubleScalarSqrt(muDoubleScalarAbs(c5_x->re));
    } else {
      c5_yr = muDoubleScalarSqrt(c5_x->re);
      c5_yi = 0.0;
    }
  } else if (c5_x->re == 0.0) {
    if (c5_x->im < 0.0) {
      c5_b_x = -c5_x->im;
      c5_c_x = c5_b_x;
      c5_d_x = c5_c_x;
      c5_z = c5_d_x / 2.0;
      c5_yr = muDoubleScalarSqrt(c5_z);
      c5_yi = -c5_yr;
    } else {
      c5_e_x = c5_x->im;
      c5_f_x = c5_e_x;
      c5_g_x = c5_f_x;
      c5_b_z = c5_g_x / 2.0;
      c5_yr = muDoubleScalarSqrt(c5_b_z);
      c5_yi = c5_yr;
    }
  } else {
    c5_b3 = muDoubleScalarIsNaN(c5_x->re);
    c5_b4 = muDoubleScalarIsNaN(c5_x->im);
    c5_b = (c5_b3 || c5_b4);
    if (c5_b) {
      c5_yr = rtNaN;
      c5_yi = rtNaN;
    } else {
      c5_h_x = c5_x->im;
      c5_b_b = muDoubleScalarIsInf(c5_h_x);
      if (c5_b_b) {
        c5_yr = rtInf;
        c5_yi = c5_x->im;
      } else {
        c5_i_x = c5_x->re;
        c5_c_b = muDoubleScalarIsInf(c5_i_x);
        if (c5_c_b) {
          if (c5_x->re < 0.0) {
            c5_yr = 0.0;
            c5_yi = rtInf;
          } else {
            c5_yr = rtInf;
            c5_yi = 0.0;
          }
        } else {
          c5_absxr = muDoubleScalarAbs(c5_x->re);
          c5_absxi = muDoubleScalarAbs(c5_x->im);
          c5_realmax(chartInstance);
          guard1 = false;
          if (c5_absxr > 4.4942328371557893E+307) {
            guard1 = true;
          } else {
            c5_realmax(chartInstance);
            if (c5_absxi > 4.4942328371557893E+307) {
              guard1 = true;
            } else {
              c5_n_x = c5_absxr;
              c5_e_y = c5_absxi;
              c5_eml_scalar_eg(chartInstance);
              c5_b_x1 = c5_n_x;
              c5_b_x2 = c5_e_y;
              c5_b_a = c5_b_x1;
              c5_e_b = c5_b_x2;
              c5_d_z = muDoubleScalarHypot(c5_b_a, c5_e_b);
              c5_yr = muDoubleScalarSqrt((c5_d_z + c5_absxr) * 0.5);
            }
          }

          if (guard1 == true) {
            c5_absxr *= 0.5;
            c5_absxi *= 0.5;
            c5_j_x = c5_absxr;
            c5_y = c5_absxi;
            c5_eml_scalar_eg(chartInstance);
            c5_x1 = c5_j_x;
            c5_x2 = c5_y;
            c5_a = c5_x1;
            c5_d_b = c5_x2;
            c5_absxd2 = muDoubleScalarHypot(c5_a, c5_d_b);
            if (c5_absxd2 > c5_absxr) {
              c5_k_x = c5_absxr;
              c5_b_y = c5_absxd2;
              c5_l_x = c5_k_x;
              c5_c_y = c5_b_y;
              c5_m_x = c5_l_x;
              c5_d_y = c5_c_y;
              c5_c_z = c5_m_x / c5_d_y;
              c5_yr = muDoubleScalarSqrt(c5_absxd2) * muDoubleScalarSqrt(1.0 +
                c5_c_z);
            } else {
              c5_yr = muDoubleScalarSqrt(c5_absxd2) * 1.4142135623730951;
            }
          }

          if (c5_x->re > 0.0) {
            c5_o_x = c5_x->im;
            c5_f_y = c5_yr;
            c5_p_x = c5_o_x;
            c5_g_y = c5_f_y;
            c5_q_x = c5_p_x;
            c5_h_y = c5_g_y;
            c5_e_z = c5_q_x / c5_h_y;
            c5_yi = 0.5 * c5_e_z;
          } else {
            if (c5_x->im < 0.0) {
              c5_yi = -c5_yr;
            } else {
              c5_yi = c5_yr;
            }

            c5_r_x = c5_x->im;
            c5_i_y = c5_yi;
            c5_s_x = c5_r_x;
            c5_j_y = c5_i_y;
            c5_t_x = c5_s_x;
            c5_k_y = c5_j_y;
            c5_f_z = c5_t_x / c5_k_y;
            c5_yr = 0.5 * c5_f_z;
          }
        }
      }
    }
  }

  c5_x->re = c5_yr;
  c5_x->im = c5_yi;
}

static void c5_d_eml_matlab_zlascl(SFc5_Model_01InstanceStruct *chartInstance,
  real_T c5_cfrom, real_T c5_cto, creal_T c5_A[3])
{
  real_T c5_cfromc;
  real_T c5_ctoc;
  boolean_T c5_notdone;
  real_T c5_cfrom1;
  real_T c5_cto1;
  real_T c5_x;
  real_T c5_b_x;
  real_T c5_y;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_b_y;
  real_T c5_mul;
  real_T c5_e_x;
  real_T c5_f_x;
  real_T c5_c_y;
  real_T c5_g_x;
  real_T c5_h_x;
  real_T c5_d_y;
  real_T c5_a;
  int32_T c5_i122;
  boolean_T guard1 = false;
  c5_realmin(chartInstance);
  c5_eps(chartInstance);
  c5_cfromc = c5_cfrom;
  c5_ctoc = c5_cto;
  c5_notdone = true;
  while (c5_notdone) {
    c5_cfrom1 = c5_cfromc * 2.0041683600089728E-292;
    c5_cto1 = c5_ctoc / 4.9896007738368E+291;
    c5_x = c5_cfrom1;
    c5_b_x = c5_x;
    c5_y = muDoubleScalarAbs(c5_b_x);
    c5_c_x = c5_ctoc;
    c5_d_x = c5_c_x;
    c5_b_y = muDoubleScalarAbs(c5_d_x);
    guard1 = false;
    if (c5_y > c5_b_y) {
      if (c5_ctoc != 0.0) {
        c5_mul = 2.0041683600089728E-292;
        c5_notdone = true;
        c5_cfromc = c5_cfrom1;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      c5_e_x = c5_cto1;
      c5_f_x = c5_e_x;
      c5_c_y = muDoubleScalarAbs(c5_f_x);
      c5_g_x = c5_cfromc;
      c5_h_x = c5_g_x;
      c5_d_y = muDoubleScalarAbs(c5_h_x);
      if (c5_c_y > c5_d_y) {
        c5_mul = 4.9896007738368E+291;
        c5_notdone = true;
        c5_ctoc = c5_cto1;
      } else {
        c5_mul = c5_ctoc / c5_cfromc;
        c5_notdone = false;
      }
    }

    c5_a = c5_mul;
    for (c5_i122 = 0; c5_i122 < 3; c5_i122++) {
      c5_A[c5_i122].re *= c5_a;
      c5_A[c5_i122].im *= c5_a;
    }
  }
}

static void c5_b_eml_xgemm(SFc5_Model_01InstanceStruct *chartInstance, real_T
  c5_A[9], real_T c5_B[9], real_T c5_C[9])
{
  int32_T c5_i123;
  int32_T c5_i124;
  int32_T c5_i125;
  int32_T c5_i126;
  int32_T c5_i127;
  (void)chartInstance;
  for (c5_i123 = 0; c5_i123 < 3; c5_i123++) {
    c5_i124 = 0;
    for (c5_i125 = 0; c5_i125 < 3; c5_i125++) {
      c5_C[c5_i124 + c5_i123] = 0.0;
      c5_i126 = 0;
      for (c5_i127 = 0; c5_i127 < 3; c5_i127++) {
        c5_C[c5_i124 + c5_i123] += c5_A[c5_i126 + c5_i123] * c5_B[c5_i127 +
          c5_i124];
        c5_i126 += 3;
      }

      c5_i124 += 3;
    }
  }
}

static void c5_b_eml_matlab_zgghrd(SFc5_Model_01InstanceStruct *chartInstance,
  int32_T c5_ilo, int32_T c5_ihi, creal_T c5_A[9], creal_T c5_Z[9])
{
  int32_T c5_i128;
  static real_T c5_dv4[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  int32_T c5_a;
  int32_T c5_b_a;
  int32_T c5_c;
  int32_T c5_c_a;
  int32_T c5_d_a;
  int32_T c5_ihim1;
  int32_T c5_jcol;
  int32_T c5_e_a;
  int32_T c5_f_a;
  int32_T c5_jcolp1;
  int32_T c5_jrow;
  int32_T c5_g_a;
  int32_T c5_h_a;
  int32_T c5_jrowm1;
  creal_T c5_b_A;
  creal_T c5_c_A;
  creal_T c5_b;
  creal_T c5_s;
  real_T c5_b_c;
  real_T c5_c_c;
  static creal_T c5_dc8 = { 0.0, 0.0 };

  real_T c5_d_c;
  int32_T c5_xrow;
  int32_T c5_yrow;
  int32_T c5_jlo;
  int32_T c5_jhi;
  int32_T c5_b_jlo;
  int32_T c5_b_jhi;
  int32_T c5_i_a;
  int32_T c5_b_b;
  int32_T c5_j_a;
  int32_T c5_c_b;
  boolean_T c5_overflow;
  int32_T c5_j;
  int32_T c5_b_j;
  real_T c5_k_a;
  creal_T c5_y;
  creal_T c5_b_s;
  creal_T c5_stemp;
  real_T c5_l_a;
  creal_T c5_d_b;
  creal_T c5_e_b;
  creal_T c5_f_b;
  creal_T c5_g_b;
  real_T c5_e_c;
  int32_T c5_xcol;
  int32_T c5_ycol;
  int32_T c5_b_ilo;
  int32_T c5_b_ihi;
  int32_T c5_c_ilo;
  int32_T c5_c_ihi;
  int32_T c5_m_a;
  int32_T c5_h_b;
  int32_T c5_n_a;
  int32_T c5_i_b;
  boolean_T c5_b_overflow;
  int32_T c5_i;
  int32_T c5_b_i;
  real_T c5_o_a;
  creal_T c5_c_s;
  real_T c5_p_a;
  creal_T c5_j_b;
  creal_T c5_k_b;
  creal_T c5_l_b;
  creal_T c5_m_b;
  real_T c5_f_c;
  int32_T c5_b_xcol;
  int32_T c5_b_ycol;
  int32_T c5_c_i;
  int32_T c5_d_i;
  real_T c5_q_a;
  creal_T c5_d_s;
  real_T c5_r_a;
  creal_T c5_n_b;
  creal_T c5_o_b;
  creal_T c5_p_b;
  creal_T c5_q_b;
  for (c5_i128 = 0; c5_i128 < 9; c5_i128++) {
    c5_Z[c5_i128].re = c5_dv4[c5_i128];
    c5_Z[c5_i128].im = 0.0;
  }

  c5_a = c5_ilo;
  c5_b_a = c5_a;
  c5_c = c5_b_a;
  if (c5_ihi < c5_c + 2) {
  } else {
    c5_c_a = c5_ihi;
    c5_d_a = c5_c_a;
    c5_ihim1 = c5_d_a;
    c5_jcol = c5_ilo;
    while (c5_jcol < c5_ihim1 - 1) {
      c5_e_a = c5_jcol;
      c5_f_a = c5_e_a + 1;
      c5_jcolp1 = c5_f_a;
      c5_jrow = c5_ihi;
      while (c5_jrow > c5_jcolp1) {
        c5_g_a = c5_jrow;
        c5_h_a = c5_g_a - 1;
        c5_jrowm1 = c5_h_a;
        c5_b_A.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_jrowm1), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].re;
        c5_b_A.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_jrowm1), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].im;
        c5_c_A.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_jrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].re;
        c5_c_A.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_jrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].im;
        c5_eml_matlab_zlartg(chartInstance, c5_b_A, c5_c_A, &c5_b_c, &c5_s,
                             &c5_b);
        c5_c_c = c5_b_c;
        c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_jrowm1), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].re = c5_b.re;
        c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_jrowm1), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c5_jcol), 1, 3, 2, 0) - 1)) - 1].im = c5_b.im;
        c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_jrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0)
               - 1)) - 1].re = c5_dc8.re;
        c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_jrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_jcol), 1, 3, 2, 0)
               - 1)) - 1].im = c5_dc8.im;
        c5_d_c = c5_c_c;
        c5_xrow = c5_jrowm1;
        c5_yrow = c5_jrow;
        c5_jlo = c5_jcolp1;
        c5_jhi = c5_ihi;
        c5_b_jlo = c5_jlo;
        c5_b_jhi = c5_jhi;
        c5_i_a = c5_b_jlo;
        c5_b_b = c5_b_jhi;
        c5_j_a = c5_i_a;
        c5_c_b = c5_b_b;
        if (c5_j_a > c5_c_b) {
          c5_overflow = false;
        } else {
          c5_b_eml_switch_helper(chartInstance);
          c5_overflow = (c5_c_b > 2147483646);
        }

        if (c5_overflow) {
          c5_check_forloop_overflow_error(chartInstance, c5_overflow);
        }

        for (c5_j = c5_b_jlo; c5_j <= c5_b_jhi; c5_j++) {
          c5_b_j = c5_j;
          c5_k_a = c5_d_c;
          c5_b.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
          c5_b.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
          c5_y.re = c5_k_a * c5_b.re;
          c5_y.im = c5_k_a * c5_b.im;
          c5_b_s.re = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re - c5_s.im * c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_yrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0)
               - 1)) - 1].im;
          c5_b_s.im = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im + c5_s.im * c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_yrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0)
               - 1)) - 1].re;
          c5_stemp.re = c5_y.re + c5_b_s.re;
          c5_stemp.im = c5_y.im + c5_b_s.im;
          c5_l_a = c5_d_c;
          c5_b.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
          c5_b.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
          c5_y.re = c5_l_a * c5_b.re;
          c5_y.im = c5_l_a * c5_b.im;
          c5_b.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
          c5_b.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
          c5_d_b = c5_b;
          c5_e_b = c5_b;
          c5_f_b = c5_b;
          c5_g_b = c5_b;
          c5_b.re = c5_s.re * c5_d_b.re + c5_s.im * c5_e_b.im;
          c5_b.im = c5_s.re * c5_f_b.im - c5_s.im * c5_g_b.re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re = c5_y.re -
            c5_b.re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im = c5_y.im -
            c5_b.im;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].re = c5_stemp.re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1].im = c5_stemp.im;
        }

        c5_s.re = -c5_s.re;
        c5_s.im = -c5_s.im;
        c5_e_c = c5_c_c;
        c5_xcol = c5_jrow;
        c5_ycol = c5_jrowm1;
        c5_b_ilo = c5_ilo;
        c5_b_ihi = c5_ihi;
        c5_c_ilo = c5_b_ilo;
        c5_c_ihi = c5_b_ihi;
        c5_m_a = c5_c_ilo;
        c5_h_b = c5_c_ihi;
        c5_n_a = c5_m_a;
        c5_i_b = c5_h_b;
        if (c5_n_a > c5_i_b) {
          c5_b_overflow = false;
        } else {
          c5_b_eml_switch_helper(chartInstance);
          c5_b_overflow = (c5_i_b > 2147483646);
        }

        if (c5_b_overflow) {
          c5_check_forloop_overflow_error(chartInstance, c5_b_overflow);
        }

        for (c5_i = c5_c_ilo; c5_i <= c5_c_ihi; c5_i++) {
          c5_b_i = c5_i;
          c5_o_a = c5_e_c;
          c5_b.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c5_b.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c5_y.re = c5_o_a * c5_b.re;
          c5_y.im = c5_o_a * c5_b.im;
          c5_c_s.re = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re - c5_s.im * c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0)
               - 1)) - 1].im;
          c5_c_s.im = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im + c5_s.im * c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0)
               - 1)) - 1].re;
          c5_stemp.re = c5_y.re + c5_c_s.re;
          c5_stemp.im = c5_y.im + c5_c_s.im;
          c5_p_a = c5_e_c;
          c5_b.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re;
          c5_b.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im;
          c5_y.re = c5_p_a * c5_b.re;
          c5_y.im = c5_p_a * c5_b.im;
          c5_b.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c5_b.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c5_j_b = c5_b;
          c5_k_b = c5_b;
          c5_l_b = c5_b;
          c5_m_b = c5_b;
          c5_b.re = c5_s.re * c5_j_b.re + c5_s.im * c5_k_b.im;
          c5_b.im = c5_s.re * c5_l_b.im - c5_s.im * c5_m_b.re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re = c5_y.re -
            c5_b.re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im = c5_y.im -
            c5_b.im;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].re = c5_stemp.re;
          c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].im = c5_stemp.im;
        }

        c5_f_c = c5_c_c;
        c5_b_xcol = c5_jrow;
        c5_b_ycol = c5_jrowm1;
        c5_b_eml_switch_helper(chartInstance);
        for (c5_c_i = 1; c5_c_i < 4; c5_c_i++) {
          c5_d_i = c5_c_i;
          c5_q_a = c5_f_c;
          c5_b.re = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c5_b.im = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c5_y.re = c5_q_a * c5_b.re;
          c5_y.im = c5_q_a * c5_b.im;
          c5_d_s.re = c5_s.re * c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].re - c5_s.im * c5_Z
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_d_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2,
                0) - 1)) - 1].im;
          c5_d_s.im = c5_s.re * c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].im + c5_s.im * c5_Z
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c5_d_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2,
                0) - 1)) - 1].re;
          c5_stemp.re = c5_y.re + c5_d_s.re;
          c5_stemp.im = c5_y.im + c5_d_s.im;
          c5_r_a = c5_f_c;
          c5_b.re = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].re;
          c5_b.im = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].im;
          c5_y.re = c5_r_a * c5_b.re;
          c5_y.im = c5_r_a * c5_b.im;
          c5_b.re = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c5_b.im = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c5_n_b = c5_b;
          c5_o_b = c5_b;
          c5_p_b = c5_b;
          c5_q_b = c5_b;
          c5_b.re = c5_s.re * c5_n_b.re + c5_s.im * c5_o_b.im;
          c5_b.im = c5_s.re * c5_p_b.im - c5_s.im * c5_q_b.re;
          c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].re = c5_y.re -
            c5_b.re;
          c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].im = c5_y.im -
            c5_b.im;
          c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1].re = c5_stemp.re;
          c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1].im = c5_stemp.im;
        }

        c5_jrow = c5_jrowm1;
      }

      c5_jcol = c5_jcolp1;
    }
  }
}

static void c5_c_eml_matlab_zhgeqz(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], int32_T c5_ilo, int32_T c5_ihi, creal_T c5_Z[9], real_T
  *c5_info, creal_T c5_alpha1[3], creal_T c5_beta1[3])
{
  static creal_T c5_dc9 = { 0.0, 0.0 };

  int32_T c5_i129;
  int32_T c5_i130;
  static creal_T c5_dc10 = { 0.0, 0.0 };

  creal_T c5_eshift;
  creal_T c5_ctemp;
  creal_T c5_rho;
  int32_T c5_i131;
  int32_T c5_i132;
  int32_T c5_i133;
  creal_T c5_b_A[9];
  real_T c5_anorm;
  real_T c5_y;
  real_T c5_atol;
  real_T c5_b_y;
  real_T c5_x;
  real_T c5_ascale;
  boolean_T c5_failed;
  int32_T c5_a;
  int32_T c5_b_a;
  int32_T c5_i134;
  int32_T c5_c_a;
  int32_T c5_d_a;
  boolean_T c5_overflow;
  int32_T c5_j;
  int32_T c5_b_j;
  int32_T c5_ifirst;
  int32_T c5_istart;
  int32_T c5_ilast;
  int32_T c5_e_a;
  int32_T c5_f_a;
  int32_T c5_ilastm1;
  int32_T c5_iiter;
  int32_T c5_g_a;
  int32_T c5_b;
  int32_T c5_h_a;
  int32_T c5_b_b;
  int32_T c5_c;
  int32_T c5_i_a;
  int32_T c5_j_a;
  int32_T c5_b_c;
  int32_T c5_c_b;
  int32_T c5_d_b;
  int32_T c5_maxit;
  boolean_T c5_goto50;
  boolean_T c5_goto60;
  boolean_T c5_goto70;
  boolean_T c5_goto90;
  int32_T c5_b_maxit;
  int32_T c5_e_b;
  int32_T c5_f_b;
  boolean_T c5_b_overflow;
  int32_T c5_jiter;
  creal_T c5_a22;
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_c_y;
  real_T c5_d_x;
  real_T c5_e_x;
  real_T c5_d_y;
  real_T c5_e_y;
  int32_T c5_k_a;
  int32_T c5_l_a;
  int32_T c5_jm1;
  boolean_T c5_ilazro;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_f_y;
  real_T c5_h_x;
  real_T c5_i_x;
  real_T c5_g_y;
  real_T c5_h_y;
  boolean_T c5_b5;
  int32_T c5_i135;
  int32_T c5_i136;
  int32_T c5_i137;
  creal_T c5_c_A;
  creal_T c5_d_A;
  creal_T c5_s;
  real_T c5_c_c;
  real_T c5_d_c;
  real_T c5_e_c;
  int32_T c5_xcol;
  int32_T c5_ycol;
  int32_T c5_b_ihi;
  int32_T c5_c_ihi;
  int32_T c5_g_b;
  int32_T c5_h_b;
  boolean_T c5_c_overflow;
  int32_T c5_i;
  int32_T c5_b_i;
  real_T c5_m_a;
  creal_T c5_a12;
  creal_T c5_b_s;
  creal_T c5_a21;
  real_T c5_n_a;
  creal_T c5_b_a22;
  creal_T c5_c_a22;
  creal_T c5_d_a22;
  creal_T c5_e_a22;
  real_T c5_f_c;
  int32_T c5_b_xcol;
  int32_T c5_b_ycol;
  int32_T c5_c_i;
  int32_T c5_d_i;
  real_T c5_o_a;
  creal_T c5_c_s;
  real_T c5_p_a;
  creal_T c5_f_a22;
  creal_T c5_g_a22;
  creal_T c5_h_a22;
  creal_T c5_i_a22;
  int32_T c5_q_a;
  int32_T c5_r_a;
  int32_T c5_s_a;
  int32_T c5_t_a;
  creal_T c5_r2;
  creal_T c5_j_a22;
  creal_T c5_b_rho;
  creal_T c5_b_a12;
  creal_T c5_c_a12;
  creal_T c5_b_a21;
  real_T c5_d7;
  real_T c5_d8;
  int32_T c5_u_a;
  int32_T c5_v_a;
  int32_T c5_jp1;
  int32_T c5_w_a;
  int32_T c5_x_a;
  real_T c5_j_x;
  real_T c5_k_x;
  real_T c5_i_y;
  real_T c5_l_x;
  real_T c5_m_x;
  real_T c5_j_y;
  real_T c5_k_y;
  real_T c5_temp;
  real_T c5_n_x;
  real_T c5_o_x;
  real_T c5_l_y;
  real_T c5_p_x;
  real_T c5_q_x;
  real_T c5_m_y;
  real_T c5_n_y;
  real_T c5_temp2;
  real_T c5_r_x;
  real_T c5_o_y;
  real_T c5_tempr;
  real_T c5_s_x;
  real_T c5_t_x;
  real_T c5_p_y;
  real_T c5_u_x;
  real_T c5_v_x;
  real_T c5_q_y;
  real_T c5_r_y;
  int32_T c5_y_a;
  int32_T c5_ab_a;
  int32_T c5_g_c;
  real_T c5_h_c;
  int32_T c5_bb_a;
  int32_T c5_cb_a;
  int32_T c5_db_a;
  int32_T c5_eb_a;
  creal_T c5_e_A;
  creal_T c5_f_A;
  real_T c5_i_c;
  real_T c5_j_c;
  int32_T c5_xrow;
  int32_T c5_yrow;
  int32_T c5_jlo;
  int32_T c5_b_jlo;
  int32_T c5_fb_a;
  int32_T c5_gb_a;
  boolean_T c5_d_overflow;
  int32_T c5_c_j;
  int32_T c5_d_j;
  real_T c5_hb_a;
  creal_T c5_d_s;
  real_T c5_ib_a;
  creal_T c5_k_a22;
  creal_T c5_l_a22;
  creal_T c5_m_a22;
  creal_T c5_n_a22;
  int32_T c5_jb_a;
  int32_T c5_kb_a;
  int32_T c5_k_c;
  int32_T c5_w_x;
  int32_T c5_s_y;
  int32_T c5_x_x;
  real_T c5_l_c;
  int32_T c5_c_xcol;
  int32_T c5_c_ycol;
  int32_T c5_d_ihi;
  int32_T c5_e_ihi;
  int32_T c5_i_b;
  int32_T c5_j_b;
  boolean_T c5_e_overflow;
  int32_T c5_e_i;
  int32_T c5_f_i;
  real_T c5_lb_a;
  creal_T c5_e_s;
  real_T c5_mb_a;
  creal_T c5_o_a22;
  creal_T c5_p_a22;
  creal_T c5_q_a22;
  creal_T c5_r_a22;
  real_T c5_m_c;
  int32_T c5_d_xcol;
  int32_T c5_d_ycol;
  int32_T c5_g_i;
  int32_T c5_h_i;
  real_T c5_nb_a;
  creal_T c5_f_s;
  real_T c5_ob_a;
  creal_T c5_s_a22;
  creal_T c5_t_a22;
  creal_T c5_u_a22;
  creal_T c5_v_a22;
  int32_T c5_b_ilast;
  int32_T c5_k_b;
  int32_T c5_l_b;
  boolean_T c5_f_overflow;
  int32_T c5_k;
  int32_T c5_b_k;
  int32_T c5_i138;
  int32_T c5_pb_a;
  int32_T c5_qb_a;
  int32_T c5_i139;
  int32_T c5_m_b;
  int32_T c5_n_b;
  boolean_T c5_g_overflow;
  int32_T c5_e_j;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  int32_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T guard11 = false;
  c5_dc9.re = rtNaN;
  for (c5_i129 = 0; c5_i129 < 3; c5_i129++) {
    c5_alpha1[c5_i129].re = 0.0;
    c5_alpha1[c5_i129].im = 0.0;
  }

  for (c5_i130 = 0; c5_i130 < 3; c5_i130++) {
    c5_beta1[c5_i130].re = 1.0;
    c5_beta1[c5_i130].im = 0.0;
  }

  c5_eps(chartInstance);
  c5_realmin(chartInstance);
  c5_eshift = c5_dc10;
  c5_ctemp = c5_dc10;
  c5_rho = c5_dc10;
  c5_i131 = 0;
  for (c5_i132 = 0; c5_i132 < 3; c5_i132++) {
    for (c5_i133 = 0; c5_i133 < 3; c5_i133++) {
      c5_b_A[c5_i133 + c5_i131] = c5_A[c5_i133 + c5_i131];
    }

    c5_i131 += 3;
  }

  c5_anorm = c5_eml_matlab_zlanhs(chartInstance, c5_b_A, c5_ilo, c5_ihi);
  c5_y = 2.2204460492503131E-16 * c5_anorm;
  c5_atol = 2.2250738585072014E-308;
  if (c5_y > 2.2250738585072014E-308) {
    c5_atol = c5_y;
  }

  c5_b_y = c5_anorm;
  c5_x = 2.2250738585072014E-308;
  if (c5_b_y > 2.2250738585072014E-308) {
    c5_x = c5_b_y;
  }

  c5_ascale = 1.0 / c5_x;
  c5_failed = true;
  c5_a = c5_ihi;
  c5_b_a = c5_a + 1;
  c5_i134 = c5_b_a;
  c5_c_a = c5_i134;
  c5_d_a = c5_c_a;
  if (c5_d_a > 3) {
    c5_overflow = false;
  } else {
    c5_b_eml_switch_helper(chartInstance);
    c5_overflow = false;
  }

  if (c5_overflow) {
    c5_check_forloop_overflow_error(chartInstance, c5_overflow);
  }

  for (c5_j = c5_i134; c5_j < 4; c5_j++) {
    c5_b_j = c5_j;
    c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_b_j), 1, 3, 1, 0) - 1].re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK(
      "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
    c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c5_b_j), 1, 3, 1, 0) - 1].im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK(
      "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
  }

  guard1 = false;
  guard2 = false;
  if (c5_ihi >= c5_ilo) {
    c5_ifirst = c5_ilo;
    c5_istart = c5_ilo;
    c5_ilast = c5_ihi;
    c5_e_a = c5_ilast;
    c5_f_a = c5_e_a - 1;
    c5_ilastm1 = c5_f_a;
    c5_iiter = 0;
    c5_g_a = c5_ihi;
    c5_b = c5_ilo;
    c5_h_a = c5_g_a;
    c5_b_b = c5_b;
    c5_c = c5_h_a - c5_b_b;
    c5_i_a = c5_c;
    c5_j_a = c5_i_a;
    c5_b_c = c5_j_a;
    c5_c_b = c5_b_c + 1;
    c5_d_b = c5_c_b;
    c5_maxit = 30 * c5_d_b;
    c5_goto50 = false;
    c5_goto60 = false;
    c5_goto70 = false;
    c5_goto90 = false;
    c5_b_maxit = c5_maxit;
    c5_e_b = c5_b_maxit;
    c5_f_b = c5_e_b;
    if (1 > c5_f_b) {
      c5_b_overflow = false;
    } else {
      c5_b_eml_switch_helper(chartInstance);
      c5_b_overflow = (c5_f_b > 2147483646);
    }

    if (c5_b_overflow) {
      c5_check_forloop_overflow_error(chartInstance, c5_b_overflow);
    }

    c5_jiter = 1;
    do {
      exitg1 = 0;
      if (c5_jiter <= c5_b_maxit) {
        if (c5_ilast == c5_ilo) {
          c5_goto60 = true;
        } else {
          c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            re;
          c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            im;
          c5_b_x = c5_a22.re;
          c5_c_x = c5_b_x;
          c5_c_y = muDoubleScalarAbs(c5_c_x);
          c5_d_x = c5_a22.im;
          c5_e_x = c5_d_x;
          c5_d_y = muDoubleScalarAbs(c5_e_x);
          c5_e_y = c5_c_y + c5_d_y;
          if (c5_e_y <= c5_atol) {
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c5_dc10.re;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c5_dc10.im;
            c5_goto60 = true;
          } else {
            c5_b_j = c5_ilastm1;
            exitg3 = false;
            while ((exitg3 == false) && (c5_b_j >= c5_ilo)) {
              c5_k_a = c5_b_j;
              c5_l_a = c5_k_a - 1;
              c5_jm1 = c5_l_a;
              if (c5_b_j == c5_ilo) {
                c5_ilazro = true;
              } else {
                c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c5_f_x = c5_a22.re;
                c5_g_x = c5_f_x;
                c5_f_y = muDoubleScalarAbs(c5_g_x);
                c5_h_x = c5_a22.im;
                c5_i_x = c5_h_x;
                c5_g_y = muDoubleScalarAbs(c5_i_x);
                c5_h_y = c5_f_y + c5_g_y;
                if (c5_h_y <= c5_atol) {
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) -
                           1)) - 1].re = c5_dc10.re;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) -
                           1)) - 1].im = c5_dc10.im;
                  c5_ilazro = true;
                } else {
                  c5_ilazro = false;
                }
              }

              if (c5_ilazro) {
                c5_ifirst = c5_b_j;
                c5_goto70 = true;
                exitg3 = true;
              } else {
                c5_b_j = c5_jm1;
              }
            }
          }
        }

        guard3 = false;
        guard4 = false;
        if (c5_goto50) {
          guard4 = true;
        } else if (c5_goto60) {
          guard4 = true;
        } else if (c5_goto70) {
          guard3 = true;
        } else {
          c5_b5 = false;
        }

        if (guard4 == true) {
          guard3 = true;
        }

        if (guard3 == true) {
          c5_b5 = true;
        }

        if (!c5_b5) {
          for (c5_i135 = 0; c5_i135 < 3; c5_i135++) {
            c5_alpha1[c5_i135] = c5_dc9;
          }

          for (c5_i136 = 0; c5_i136 < 3; c5_i136++) {
            c5_beta1[c5_i136] = c5_dc9;
          }

          for (c5_i137 = 0; c5_i137 < 9; c5_i137++) {
            c5_Z[c5_i137] = c5_dc9;
          }

          *c5_info = -1.0;
          exitg1 = 1;
        } else {
          if (c5_goto50) {
            c5_goto50 = false;
            c5_c_A.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].
              re;
            c5_c_A.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].
              im;
            c5_d_A.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .re;
            c5_d_A.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .im;
            c5_eml_matlab_zlartg(chartInstance, c5_c_A, c5_d_A, &c5_c_c, &c5_s,
                                 &c5_a22);
            c5_d_c = c5_c_c;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].re =
              c5_a22.re;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].im =
              c5_a22.im;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c5_dc10.re;
            c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c5_dc10.im;
            c5_e_c = c5_d_c;
            c5_xcol = c5_ilast;
            c5_ycol = c5_ilastm1;
            c5_b_ihi = c5_ilastm1;
            c5_c_ihi = c5_b_ihi;
            c5_g_b = c5_c_ihi;
            c5_h_b = c5_g_b;
            if (1 > c5_h_b) {
              c5_c_overflow = false;
            } else {
              c5_b_eml_switch_helper(chartInstance);
              c5_c_overflow = (c5_h_b > 2147483646);
            }

            if (c5_c_overflow) {
              c5_check_forloop_overflow_error(chartInstance, c5_c_overflow);
            }

            for (c5_i = 1; c5_i <= c5_c_ihi; c5_i++) {
              c5_b_i = c5_i;
              c5_m_a = c5_e_c;
              c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c5_a12.re = c5_m_a * c5_a22.re;
              c5_a12.im = c5_m_a * c5_a22.im;
              c5_b_s.re = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re - c5_s.im *
                c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1))
                - 1].im;
              c5_b_s.im = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im + c5_s.im *
                c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1))
                - 1].re;
              c5_a21.re = c5_a12.re + c5_b_s.re;
              c5_a21.im = c5_a12.im + c5_b_s.im;
              c5_n_a = c5_e_c;
              c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c5_a12.re = c5_n_a * c5_a22.re;
              c5_a12.im = c5_n_a * c5_a22.im;
              c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c5_b_a22 = c5_a22;
              c5_c_a22 = c5_a22;
              c5_d_a22 = c5_a22;
              c5_e_a22 = c5_a22;
              c5_a22.re = c5_s.re * c5_b_a22.re + c5_s.im * c5_c_a22.im;
              c5_a22.im = c5_s.re * c5_d_a22.im - c5_s.im * c5_e_a22.re;
              c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].re =
                c5_a12.re - c5_a22.re;
              c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ycol), 1, 3, 2, 0) - 1)) - 1].im =
                c5_a12.im - c5_a22.im;
              c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].re =
                c5_a21.re;
              c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_xcol), 1, 3, 2, 0) - 1)) - 1].im =
                c5_a21.im;
            }

            c5_f_c = c5_d_c;
            c5_b_xcol = c5_ilast;
            c5_b_ycol = c5_ilastm1;
            c5_b_eml_switch_helper(chartInstance);
            for (c5_c_i = 1; c5_c_i < 4; c5_c_i++) {
              c5_d_i = c5_c_i;
              c5_o_a = c5_f_c;
              c5_a22.re = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c5_a22.im = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c5_a12.re = c5_o_a * c5_a22.re;
              c5_a12.im = c5_o_a * c5_a22.im;
              c5_c_s.re = c5_s.re * c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].re - c5_s.im *
                c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) -
                       1)) - 1].im;
              c5_c_s.im = c5_s.re * c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].im + c5_s.im *
                c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) -
                       1)) - 1].re;
              c5_a21.re = c5_a12.re + c5_c_s.re;
              c5_a21.im = c5_a12.im + c5_c_s.im;
              c5_p_a = c5_f_c;
              c5_a22.re = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c5_a22.im = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c5_a12.re = c5_p_a * c5_a22.re;
              c5_a12.im = c5_p_a * c5_a22.im;
              c5_a22.re = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c5_a22.im = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c5_f_a22 = c5_a22;
              c5_g_a22 = c5_a22;
              c5_h_a22 = c5_a22;
              c5_i_a22 = c5_a22;
              c5_a22.re = c5_s.re * c5_f_a22.re + c5_s.im * c5_g_a22.im;
              c5_a22.im = c5_s.re * c5_h_a22.im - c5_s.im * c5_i_a22.re;
              c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].re =
                c5_a12.re - c5_a22.re;
              c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_ycol), 1, 3, 2, 0) - 1)) - 1].im =
                c5_a12.im - c5_a22.im;
              c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1].re =
                c5_a21.re;
              c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_b_xcol), 1, 3, 2, 0) - 1)) - 1].im =
                c5_a21.im;
            }

            c5_goto60 = true;
          }

          guard11 = false;
          if (c5_goto60) {
            c5_goto60 = false;
            c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) - 1].re =
              c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].re;
            c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) - 1].im =
              c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) - 1].im;
            c5_ilast = c5_ilastm1;
            c5_q_a = c5_ilast;
            c5_r_a = c5_q_a - 1;
            c5_ilastm1 = c5_r_a;
            if (c5_ilast < c5_ilo) {
              c5_failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              c5_iiter = 0;
              c5_eshift = c5_dc10;
              guard11 = true;
            }
          } else {
            if (c5_goto70) {
              c5_goto70 = false;
              c5_s_a = c5_iiter;
              c5_t_a = c5_s_a + 1;
              c5_iiter = c5_t_a;
              if (c5_mod(chartInstance, c5_iiter) != 0) {
                c5_s.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c5_s.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c5_r2.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c5_r2.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c5_a12.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c5_a12.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c5_a21.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c5_a21.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c5_a22.re = c5_r2.re - c5_s.re;
                c5_a22.im = c5_r2.im - c5_s.im;
                c5_j_a22.re = -c5_a22.re;
                c5_j_a22.im = -c5_a22.im;
                c5_rho = c5_eml_div(chartInstance, c5_j_a22, 2.0);
                c5_b_rho.re = c5_rho.re * c5_rho.re - c5_rho.im * c5_rho.im;
                c5_b_rho.im = c5_rho.re * c5_rho.im + c5_rho.im * c5_rho.re;
                c5_b_a12.re = c5_a12.re * c5_a21.re - c5_a12.im * c5_a21.im;
                c5_b_a12.im = c5_a12.re * c5_a21.im + c5_a12.im * c5_a21.re;
                c5_a22.re = c5_b_rho.re + c5_b_a12.re;
                c5_a22.im = c5_b_rho.im + c5_b_a12.im;
                c5_b_sqrt(chartInstance, &c5_a22);
                c5_a12.re = c5_s.re - (c5_rho.re - c5_a22.re);
                c5_a12.im = c5_s.im - (c5_rho.im - c5_a22.im);
                c5_a21.re = c5_s.re - (c5_rho.re + c5_a22.re);
                c5_a21.im = c5_s.im - (c5_rho.im + c5_a22.im);
                c5_c_a12.re = c5_a12.re - c5_r2.re;
                c5_c_a12.im = c5_a12.im - c5_r2.im;
                c5_b_a21.re = c5_a21.re - c5_r2.re;
                c5_b_a21.im = c5_a21.im - c5_r2.im;
                c5_d7 = c5_b_abs(chartInstance, c5_c_a12);
                c5_d8 = c5_b_abs(chartInstance, c5_b_a21);
                if (c5_d7 <= c5_d8) {
                  c5_a21 = c5_a12;
                  c5_rho.re -= c5_a22.re;
                  c5_rho.im -= c5_a22.im;
                } else {
                  c5_rho.re += c5_a22.re;
                  c5_rho.im += c5_a22.im;
                }
              } else {
                c5_eshift.re += c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c5_eshift.im += c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilast), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c5_a21 = c5_eshift;
              }

              c5_b_j = c5_ilastm1;
              c5_u_a = c5_b_j;
              c5_v_a = c5_u_a + 1;
              c5_jp1 = c5_v_a;
              exitg2 = false;
              while ((exitg2 == false) && (c5_b_j > c5_ifirst)) {
                c5_w_a = c5_b_j;
                c5_x_a = c5_w_a - 1;
                c5_jm1 = c5_x_a;
                c5_istart = c5_b_j;
                c5_ctemp.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re - c5_a21.re;
                c5_ctemp.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im - c5_a21.im;
                c5_j_x = c5_ctemp.re;
                c5_k_x = c5_j_x;
                c5_i_y = muDoubleScalarAbs(c5_k_x);
                c5_l_x = c5_ctemp.im;
                c5_m_x = c5_l_x;
                c5_j_y = muDoubleScalarAbs(c5_m_x);
                c5_k_y = c5_i_y + c5_j_y;
                c5_temp = c5_ascale * c5_k_y;
                c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c5_n_x = c5_a22.re;
                c5_o_x = c5_n_x;
                c5_l_y = muDoubleScalarAbs(c5_o_x);
                c5_p_x = c5_a22.im;
                c5_q_x = c5_p_x;
                c5_m_y = muDoubleScalarAbs(c5_q_x);
                c5_n_y = c5_l_y + c5_m_y;
                c5_temp2 = c5_ascale * c5_n_y;
                c5_r_x = c5_temp;
                c5_o_y = c5_temp2;
                c5_tempr = c5_r_x;
                if (c5_o_y > c5_tempr) {
                  c5_tempr = c5_o_y;
                }

                if (c5_tempr < 1.0) {
                  if (c5_tempr != 0.0) {
                    c5_temp /= c5_tempr;
                    c5_temp2 /= c5_tempr;
                  }
                }

                c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c5_s_x = c5_a22.re;
                c5_t_x = c5_s_x;
                c5_p_y = muDoubleScalarAbs(c5_t_x);
                c5_u_x = c5_a22.im;
                c5_v_x = c5_u_x;
                c5_q_y = muDoubleScalarAbs(c5_v_x);
                c5_r_y = c5_p_y + c5_q_y;
                if (c5_r_y * c5_temp2 <= c5_temp * c5_atol) {
                  c5_goto90 = true;
                  exitg2 = true;
                } else {
                  c5_jp1 = c5_b_j;
                  c5_b_j = c5_jm1;
                }
              }

              if (!c5_goto90) {
                c5_istart = c5_ifirst;
                if (c5_istart == c5_ilastm1) {
                  c5_ctemp = c5_rho;
                } else {
                  c5_ctemp.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 2, 0) - 1))
                    - 1].re - c5_a21.re;
                  c5_ctemp.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 2, 0) - 1))
                    - 1].im - c5_a21.im;
                }

                c5_goto90 = true;
              }
            }

            if (c5_goto90) {
              c5_goto90 = false;
              c5_y_a = c5_istart;
              c5_ab_a = c5_y_a;
              c5_g_c = c5_ab_a;
              c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c5_g_c + 1)), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c5_g_c + 1)), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c5_istart), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c5_b_eml_matlab_zlartg(chartInstance, c5_ctemp, c5_a22, &c5_h_c,
                &c5_s);
              c5_d_c = c5_h_c;
              c5_b_j = c5_istart;
              c5_bb_a = c5_b_j;
              c5_cb_a = c5_bb_a - 1;
              c5_jm1 = c5_cb_a;
              while (c5_b_j < c5_ilast) {
                c5_db_a = c5_b_j;
                c5_eb_a = c5_db_a + 1;
                c5_jp1 = c5_eb_a;
                if (c5_b_j > c5_istart) {
                  c5_e_A.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_e_A.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_f_A.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_f_A.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_eml_matlab_zlartg(chartInstance, c5_e_A, c5_f_A, &c5_i_c,
                                       &c5_s, &c5_a22);
                  c5_d_c = c5_i_c;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) -
                           1)) - 1].re = c5_a22.re;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) -
                           1)) - 1].im = c5_a22.im;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) -
                           1)) - 1].re = c5_dc10.re;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_jp1), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_jm1), 1, 3, 2, 0) -
                           1)) - 1].im = c5_dc10.im;
                }

                c5_j_c = c5_d_c;
                c5_xrow = c5_b_j;
                c5_yrow = c5_jp1;
                c5_jlo = c5_b_j;
                c5_b_jlo = c5_jlo;
                c5_fb_a = c5_b_jlo;
                c5_gb_a = c5_fb_a;
                if (c5_gb_a > 3) {
                  c5_d_overflow = false;
                } else {
                  c5_b_eml_switch_helper(chartInstance);
                  c5_d_overflow = false;
                }

                if (c5_d_overflow) {
                  c5_check_forloop_overflow_error(chartInstance, c5_d_overflow);
                }

                for (c5_c_j = c5_b_jlo; c5_c_j < 4; c5_c_j++) {
                  c5_d_j = c5_c_j;
                  c5_hb_a = c5_j_c;
                  c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_a12.re = c5_hb_a * c5_a22.re;
                  c5_a12.im = c5_hb_a * c5_a22.im;
                  c5_d_s.re = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re - c5_s.im * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_d_s.im = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im + c5_s.im * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_a21.re = c5_a12.re + c5_d_s.re;
                  c5_a21.im = c5_a12.im + c5_d_s.im;
                  c5_ib_a = c5_j_c;
                  c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_a12.re = c5_ib_a * c5_a22.re;
                  c5_a12.im = c5_ib_a * c5_a22.im;
                  c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c5_k_a22 = c5_a22;
                  c5_l_a22 = c5_a22;
                  c5_m_a22 = c5_a22;
                  c5_n_a22 = c5_a22;
                  c5_a22.re = c5_s.re * c5_k_a22.re + c5_s.im * c5_l_a22.im;
                  c5_a22.im = c5_s.re * c5_m_a22.im - c5_s.im * c5_n_a22.re;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) -
                             1)) - 1].re = c5_a12.re - c5_a22.re;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_yrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) -
                             1)) - 1].im = c5_a12.im - c5_a22.im;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) -
                             1)) - 1].re = c5_a21.re;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_xrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_d_j), 1, 3, 2, 0) -
                             1)) - 1].im = c5_a21.im;
                }

                c5_s.re = -c5_s.re;
                c5_s.im = -c5_s.im;
                c5_jb_a = c5_jp1;
                c5_kb_a = c5_jb_a;
                c5_k_c = c5_kb_a;
                c5_w_x = c5_k_c + 1;
                c5_s_y = c5_ilast;
                c5_x_x = c5_w_x;
                if (c5_s_y < c5_x_x) {
                  c5_x_x = c5_s_y;
                }

                c5_l_c = c5_d_c;
                c5_c_xcol = c5_jp1;
                c5_c_ycol = c5_b_j;
                c5_d_ihi = c5_x_x;
                c5_e_ihi = c5_d_ihi;
                c5_i_b = c5_e_ihi;
                c5_j_b = c5_i_b;
                if (1 > c5_j_b) {
                  c5_e_overflow = false;
                } else {
                  c5_b_eml_switch_helper(chartInstance);
                  c5_e_overflow = (c5_j_b > 2147483646);
                }

                if (c5_e_overflow) {
                  c5_check_forloop_overflow_error(chartInstance, c5_e_overflow);
                }

                for (c5_e_i = 1; c5_e_i <= c5_e_ihi; c5_e_i++) {
                  c5_f_i = c5_e_i;
                  c5_lb_a = c5_l_c;
                  c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_a12.re = c5_lb_a * c5_a22.re;
                  c5_a12.im = c5_lb_a * c5_a22.im;
                  c5_e_s.re = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].re - c5_s.im * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_e_s.im = c5_s.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].im + c5_s.im * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a21.re = c5_a12.re + c5_e_s.re;
                  c5_a21.im = c5_a12.im + c5_e_s.im;
                  c5_mb_a = c5_l_c;
                  c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_a12.re = c5_mb_a * c5_a22.re;
                  c5_a12.im = c5_mb_a * c5_a22.im;
                  c5_a22.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a22.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_o_a22 = c5_a22;
                  c5_p_a22 = c5_a22;
                  c5_q_a22 = c5_a22;
                  c5_r_a22 = c5_a22;
                  c5_a22.re = c5_s.re * c5_o_a22.re + c5_s.im * c5_p_a22.im;
                  c5_a22.im = c5_s.re * c5_q_a22.im - c5_s.im * c5_r_a22.re;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_c_ycol), 1, 3, 2, 0)
                           - 1)) - 1].re = c5_a12.re - c5_a22.re;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_c_ycol), 1, 3, 2, 0)
                           - 1)) - 1].im = c5_a12.im - c5_a22.im;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_c_xcol), 1, 3, 2, 0)
                           - 1)) - 1].re = c5_a21.re;
                  c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_c_xcol), 1, 3, 2, 0)
                           - 1)) - 1].im = c5_a21.im;
                }

                c5_m_c = c5_d_c;
                c5_d_xcol = c5_jp1;
                c5_d_ycol = c5_b_j;
                c5_b_eml_switch_helper(chartInstance);
                for (c5_g_i = 1; c5_g_i < 4; c5_g_i++) {
                  c5_h_i = c5_g_i;
                  c5_nb_a = c5_m_c;
                  c5_a22.re = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a22.im = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_a12.re = c5_nb_a * c5_a22.re;
                  c5_a12.im = c5_nb_a * c5_a22.im;
                  c5_f_s.re = c5_s.re * c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].re - c5_s.im * c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_f_s.im = c5_s.re * c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].im + c5_s.im * c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a21.re = c5_a12.re + c5_f_s.re;
                  c5_a21.im = c5_a12.im + c5_f_s.im;
                  c5_ob_a = c5_m_c;
                  c5_a22.re = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a22.im = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_a12.re = c5_ob_a * c5_a22.re;
                  c5_a12.im = c5_ob_a * c5_a22.im;
                  c5_a22.re = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c5_a22.im = c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c5_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c5_s_a22 = c5_a22;
                  c5_t_a22 = c5_a22;
                  c5_u_a22 = c5_a22;
                  c5_v_a22 = c5_a22;
                  c5_a22.re = c5_s.re * c5_s_a22.re + c5_s.im * c5_t_a22.im;
                  c5_a22.im = c5_s.re * c5_u_a22.im - c5_s.im * c5_v_a22.re;
                  c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_d_ycol), 1, 3, 2, 0)
                           - 1)) - 1].re = c5_a12.re - c5_a22.re;
                  c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_d_ycol), 1, 3, 2, 0)
                           - 1)) - 1].im = c5_a12.im - c5_a22.im;
                  c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_d_xcol), 1, 3, 2, 0)
                           - 1)) - 1].re = c5_a21.re;
                  c5_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c5_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c5_d_xcol), 1, 3, 2, 0)
                           - 1)) - 1].im = c5_a21.im;
                }

                c5_jm1 = c5_b_j;
                c5_b_j = c5_jp1;
              }
            }

            guard11 = true;
          }

          if (guard11 == true) {
            c5_jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2 == true) {
    if (c5_failed) {
      *c5_info = (real_T)c5_ilast;
      c5_b_ilast = c5_ilast;
      c5_k_b = c5_b_ilast;
      c5_l_b = c5_k_b;
      if (1 > c5_l_b) {
        c5_f_overflow = false;
      } else {
        c5_b_eml_switch_helper(chartInstance);
        c5_f_overflow = (c5_l_b > 2147483646);
      }

      if (c5_f_overflow) {
        c5_check_forloop_overflow_error(chartInstance, c5_f_overflow);
      }

      for (c5_k = 1; c5_k <= c5_b_ilast; c5_k++) {
        c5_b_k = c5_k;
        c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c5_b_k), 1, 3, 1, 0) - 1].re = c5_dc9.re;
        c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c5_b_k), 1, 3, 1, 0) - 1].im = c5_dc9.im;
        c5_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c5_b_k), 1, 3, 1, 0) - 1].re = c5_dc9.re;
        c5_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c5_b_k), 1, 3, 1, 0) - 1].im = c5_dc9.im;
      }

      for (c5_i138 = 0; c5_i138 < 9; c5_i138++) {
        c5_Z[c5_i138] = c5_dc9;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1 == true) {
    c5_pb_a = c5_ilo;
    c5_qb_a = c5_pb_a - 1;
    c5_i139 = c5_qb_a;
    c5_m_b = c5_i139;
    c5_n_b = c5_m_b;
    if (1 > c5_n_b) {
      c5_g_overflow = false;
    } else {
      c5_b_eml_switch_helper(chartInstance);
      c5_g_overflow = (c5_n_b > 2147483646);
    }

    if (c5_g_overflow) {
      c5_check_forloop_overflow_error(chartInstance, c5_g_overflow);
    }

    for (c5_e_j = 1; c5_e_j <= c5_i139; c5_e_j++) {
      c5_b_j = c5_e_j;
      c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c5_b_j), 1, 3, 1, 0) - 1].re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
        ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
        c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c5_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c5_b_j), 1, 3, 1, 0) - 1].im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
        ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c5_b_j), 1, 3, 1, 0) + 3 *
        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
        c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
    }

    *c5_info = 0.0;
  }
}

static void c5_b_eml_matlab_ztgevc(SFc5_Model_01InstanceStruct *chartInstance,
  creal_T c5_A[9], creal_T c5_V[9])
{
  int32_T c5_i140;
  creal_T c5_work1[3];
  int32_T c5_i141;
  creal_T c5_work2[3];
  int32_T c5_i142;
  real_T c5_rworka[3];
  creal_T c5_ca;
  real_T c5_x;
  real_T c5_b_x;
  real_T c5_y;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_b_y;
  real_T c5_anorm;
  int32_T c5_j;
  real_T c5_b_j;
  real_T c5_d9;
  int32_T c5_i143;
  int32_T c5_i;
  real_T c5_b_i;
  real_T c5_e_x;
  real_T c5_f_x;
  real_T c5_c_y;
  real_T c5_g_x;
  real_T c5_h_x;
  real_T c5_d_y;
  real_T c5_e_y;
  real_T c5_i_x;
  real_T c5_j_x;
  real_T c5_f_y;
  real_T c5_k_x;
  real_T c5_l_x;
  real_T c5_g_y;
  real_T c5_h_y;
  real_T c5_i_y;
  real_T c5_m_x;
  real_T c5_n_x;
  real_T c5_ascale;
  int32_T c5_je;
  real_T c5_b_je;
  real_T c5_ieig;
  real_T c5_o_x;
  real_T c5_p_x;
  real_T c5_j_y;
  real_T c5_q_x;
  real_T c5_r_x;
  real_T c5_k_y;
  real_T c5_l_y;
  real_T c5_s_x;
  real_T c5_t_x;
  real_T c5_temp;
  real_T c5_sbeta;
  real_T c5_a;
  real_T c5_b;
  creal_T c5_salpha;
  real_T c5_acoeff;
  boolean_T c5_b6;
  boolean_T c5_lscalea;
  real_T c5_u_x;
  real_T c5_v_x;
  real_T c5_m_y;
  real_T c5_w_x;
  real_T c5_x_x;
  real_T c5_n_y;
  real_T c5_o_y;
  real_T c5_y_x;
  real_T c5_ab_x;
  real_T c5_p_y;
  real_T c5_bb_x;
  real_T c5_cb_x;
  real_T c5_q_y;
  real_T c5_r_y;
  boolean_T c5_b7;
  boolean_T c5_lscaleb;
  real_T c5_scale;
  real_T c5_db_x;
  real_T c5_eb_x;
  real_T c5_fb_x;
  real_T c5_gb_x;
  real_T c5_s_y;
  real_T c5_hb_x;
  real_T c5_ib_x;
  real_T c5_t_y;
  real_T c5_u_y;
  real_T c5_v_y;
  real_T c5_jb_x;
  real_T c5_kb_x;
  real_T c5_w_y;
  real_T c5_lb_x;
  real_T c5_mb_x;
  real_T c5_x_y;
  real_T c5_y_y;
  real_T c5_nb_x;
  real_T c5_z;
  real_T c5_ob_x;
  real_T c5_ab_y;
  real_T c5_b_a;
  real_T c5_c_a;
  real_T c5_acoefa;
  real_T c5_pb_x;
  real_T c5_qb_x;
  real_T c5_bb_y;
  real_T c5_rb_x;
  real_T c5_sb_x;
  real_T c5_cb_y;
  real_T c5_bcoefa;
  int32_T c5_jr;
  real_T c5_b_jr;
  static creal_T c5_dc11 = { 0.0, 0.0 };

  static creal_T c5_dc12 = { 1.0, 0.0 };

  real_T c5_tb_x;
  real_T c5_db_y;
  real_T c5_dmin;
  real_T c5_d10;
  int32_T c5_i144;
  int32_T c5_c_jr;
  real_T c5_d_a;
  real_T c5_d11;
  int32_T c5_i145;
  int32_T c5_c_j;
  real_T c5_e_a;
  creal_T c5_d;
  real_T c5_ub_x;
  real_T c5_vb_x;
  real_T c5_eb_y;
  real_T c5_wb_x;
  real_T c5_xb_x;
  real_T c5_fb_y;
  real_T c5_gb_y;
  real_T c5_yb_x;
  real_T c5_ac_x;
  real_T c5_hb_y;
  real_T c5_bc_x;
  real_T c5_cc_x;
  real_T c5_ib_y;
  real_T c5_jb_y;
  real_T c5_dc_x;
  real_T c5_ec_x;
  real_T c5_kb_y;
  real_T c5_fc_x;
  real_T c5_gc_x;
  real_T c5_lb_y;
  real_T c5_mb_y;
  real_T c5_hc_x;
  real_T c5_ic_x;
  real_T c5_nb_y;
  real_T c5_jc_x;
  real_T c5_kc_x;
  real_T c5_ob_y;
  real_T c5_pb_y;
  real_T c5_lc_x;
  real_T c5_mc_x;
  real_T c5_qb_y;
  real_T c5_nc_x;
  real_T c5_oc_x;
  real_T c5_rb_y;
  real_T c5_sb_y;
  real_T c5_c_je;
  int32_T c5_i146;
  int32_T c5_d_jr;
  real_T c5_f_a;
  real_T c5_pc_x;
  real_T c5_qc_x;
  real_T c5_tb_y;
  real_T c5_rc_x;
  real_T c5_sc_x;
  real_T c5_ub_y;
  real_T c5_vb_y;
  real_T c5_tc_x;
  real_T c5_uc_x;
  real_T c5_wb_y;
  real_T c5_vc_x;
  real_T c5_wc_x;
  real_T c5_xb_y;
  real_T c5_yb_y;
  real_T c5_d_je;
  int32_T c5_i147;
  int32_T c5_e_jr;
  real_T c5_g_a;
  real_T c5_h_a;
  real_T c5_d12;
  int32_T c5_i148;
  int32_T c5_f_jr;
  creal_T c5_b_ca;
  int32_T c5_g_jr;
  real_T c5_e_je;
  int32_T c5_i149;
  int32_T c5_jc;
  real_T c5_b_jc;
  int32_T c5_h_jr;
  creal_T c5_b_V;
  real_T c5_xc_x;
  real_T c5_yc_x;
  real_T c5_ac_y;
  real_T c5_ad_x;
  real_T c5_bd_x;
  real_T c5_bc_y;
  real_T c5_xmx;
  int32_T c5_i_jr;
  real_T c5_cd_x;
  real_T c5_dd_x;
  real_T c5_cc_y;
  real_T c5_ed_x;
  real_T c5_fd_x;
  real_T c5_dc_y;
  real_T c5_ec_y;
  real_T c5_fc_y;
  int32_T c5_j_jr;
  real_T c5_i_a;
  int32_T c5_k_jr;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  for (c5_i140 = 0; c5_i140 < 3; c5_i140++) {
    c5_work1[c5_i140].re = 0.0;
    c5_work1[c5_i140].im = 0.0;
  }

  for (c5_i141 = 0; c5_i141 < 3; c5_i141++) {
    c5_work2[c5_i141].re = 0.0;
    c5_work2[c5_i141].im = 0.0;
  }

  c5_eps(chartInstance);
  c5_realmin(chartInstance);
  for (c5_i142 = 0; c5_i142 < 3; c5_i142++) {
    c5_rworka[c5_i142] = 0.0;
  }

  c5_ca = c5_A[0];
  c5_x = c5_ca.re;
  c5_b_x = c5_x;
  c5_y = muDoubleScalarAbs(c5_b_x);
  c5_c_x = c5_ca.im;
  c5_d_x = c5_c_x;
  c5_b_y = muDoubleScalarAbs(c5_d_x);
  c5_anorm = c5_y + c5_b_y;
  for (c5_j = 0; c5_j < 2; c5_j++) {
    c5_b_j = 2.0 + (real_T)c5_j;
    c5_d9 = c5_b_j - 1.0;
    c5_i143 = (int32_T)c5_d9;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c5_d9, mxDOUBLE_CLASS, c5_i143);
    for (c5_i = 0; c5_i < c5_i143; c5_i++) {
      c5_b_i = 1.0 + (real_T)c5_i;
      c5_ca.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_i), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c5_ca.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_i), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
      c5_e_x = c5_ca.re;
      c5_f_x = c5_e_x;
      c5_c_y = muDoubleScalarAbs(c5_f_x);
      c5_g_x = c5_ca.im;
      c5_h_x = c5_g_x;
      c5_d_y = muDoubleScalarAbs(c5_h_x);
      c5_e_y = c5_c_y + c5_d_y;
      c5_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_j), 1, 3, 1, 0) - 1] = c5_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)_SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1] + c5_e_y;
    }

    c5_ca.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c5_b_j), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
    c5_ca.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c5_b_j), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
    c5_i_x = c5_ca.re;
    c5_j_x = c5_i_x;
    c5_f_y = muDoubleScalarAbs(c5_j_x);
    c5_k_x = c5_ca.im;
    c5_l_x = c5_k_x;
    c5_g_y = muDoubleScalarAbs(c5_l_x);
    c5_h_y = c5_f_y + c5_g_y;
    c5_i_y = c5_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1] + c5_h_y;
    if (c5_i_y > c5_anorm) {
      c5_anorm = c5_i_y;
    }
  }

  c5_m_x = c5_anorm;
  c5_n_x = c5_m_x;
  if (2.2250738585072014E-308 > c5_n_x) {
    c5_n_x = 2.2250738585072014E-308;
  }

  c5_ascale = 1.0 / c5_n_x;
  for (c5_je = 0; c5_je < 3; c5_je++) {
    c5_b_je = 3.0 + -(real_T)c5_je;
    c5_ieig = c5_b_je;
    c5_ca.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c5_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c5_b_je), 1, 3, 2, 0) - 1)) - 1].re;
    c5_ca.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c5_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c5_b_je), 1, 3, 2, 0) - 1)) - 1].im;
    c5_o_x = c5_ca.re;
    c5_p_x = c5_o_x;
    c5_j_y = muDoubleScalarAbs(c5_p_x);
    c5_q_x = c5_ca.im;
    c5_r_x = c5_q_x;
    c5_k_y = muDoubleScalarAbs(c5_r_x);
    c5_l_y = c5_j_y + c5_k_y;
    c5_s_x = c5_l_y * c5_ascale;
    c5_t_x = c5_s_x;
    if (1.0 > c5_t_x) {
      c5_t_x = 1.0;
    }

    c5_temp = 1.0 / c5_t_x;
    c5_sbeta = c5_temp;
    c5_a = c5_temp;
    c5_ca.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c5_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c5_b_je), 1, 3, 2, 0) - 1)) - 1].re;
    c5_ca.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c5_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c5_b_je), 1, 3, 2, 0) - 1)) - 1].im;
    c5_ca.re *= c5_a;
    c5_ca.im *= c5_a;
    c5_b = c5_ascale;
    c5_salpha.re = c5_b * c5_ca.re;
    c5_salpha.im = c5_b * c5_ca.im;
    c5_acoeff = c5_sbeta * c5_ascale;
    guard3 = false;
    if (c5_abs(chartInstance, c5_sbeta) >= 2.2250738585072014E-308) {
      if (c5_abs(chartInstance, c5_acoeff) < 3.0062525400134592E-292) {
        c5_b6 = true;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }

    if (guard3 == true) {
      c5_b6 = false;
    }

    c5_lscalea = c5_b6;
    c5_u_x = c5_salpha.re;
    c5_v_x = c5_u_x;
    c5_m_y = muDoubleScalarAbs(c5_v_x);
    c5_w_x = c5_salpha.im;
    c5_x_x = c5_w_x;
    c5_n_y = muDoubleScalarAbs(c5_x_x);
    c5_o_y = c5_m_y + c5_n_y;
    guard2 = false;
    if (c5_o_y >= 2.2250738585072014E-308) {
      c5_y_x = c5_salpha.re;
      c5_ab_x = c5_y_x;
      c5_p_y = muDoubleScalarAbs(c5_ab_x);
      c5_bb_x = c5_salpha.im;
      c5_cb_x = c5_bb_x;
      c5_q_y = muDoubleScalarAbs(c5_cb_x);
      c5_r_y = c5_p_y + c5_q_y;
      if (c5_r_y < 3.0062525400134592E-292) {
        c5_b7 = true;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2 == true) {
      c5_b7 = false;
    }

    c5_lscaleb = c5_b7;
    c5_scale = 1.0;
    if (c5_lscalea) {
      c5_db_x = c5_anorm;
      c5_eb_x = c5_db_x;
      if (3.3264005158911995E+291 < c5_eb_x) {
        c5_eb_x = 3.3264005158911995E+291;
      }

      c5_scale = 3.0062525400134592E-292 / c5_abs(chartInstance, c5_sbeta) *
        c5_eb_x;
    }

    if (c5_lscaleb) {
      c5_fb_x = c5_salpha.re;
      c5_gb_x = c5_fb_x;
      c5_s_y = muDoubleScalarAbs(c5_gb_x);
      c5_hb_x = c5_salpha.im;
      c5_ib_x = c5_hb_x;
      c5_t_y = muDoubleScalarAbs(c5_ib_x);
      c5_u_y = c5_s_y + c5_t_y;
      c5_v_y = 3.0062525400134592E-292 / c5_u_y;
      if (c5_v_y > c5_scale) {
        c5_scale = c5_v_y;
      }
    }

    guard1 = false;
    if (c5_lscalea) {
      guard1 = true;
    } else {
      if (c5_lscaleb) {
        guard1 = true;
      }
    }

    if (guard1 == true) {
      c5_jb_x = c5_salpha.re;
      c5_kb_x = c5_jb_x;
      c5_w_y = muDoubleScalarAbs(c5_kb_x);
      c5_lb_x = c5_salpha.im;
      c5_mb_x = c5_lb_x;
      c5_x_y = muDoubleScalarAbs(c5_mb_x);
      c5_y_y = c5_w_y + c5_x_y;
      c5_nb_x = c5_abs(chartInstance, c5_acoeff);
      c5_z = c5_y_y;
      c5_ob_x = c5_nb_x;
      if (1.0 > c5_ob_x) {
        c5_ob_x = 1.0;
      }

      if (c5_z > c5_ob_x) {
        c5_ob_x = c5_z;
      }

      c5_ab_y = 1.0 / (2.2250738585072014E-308 * c5_ob_x);
      if (c5_ab_y < c5_scale) {
        c5_scale = c5_ab_y;
      }

      if (c5_lscalea) {
        c5_acoeff = c5_ascale * (c5_scale * c5_sbeta);
      } else {
        c5_acoeff *= c5_scale;
      }

      if (c5_lscaleb) {
        c5_b_a = c5_scale;
        c5_salpha.re *= c5_b_a;
        c5_salpha.im *= c5_b_a;
      } else {
        c5_c_a = c5_scale;
        c5_salpha.re *= c5_c_a;
        c5_salpha.im *= c5_c_a;
      }
    }

    c5_acoefa = c5_abs(chartInstance, c5_acoeff);
    c5_pb_x = c5_salpha.re;
    c5_qb_x = c5_pb_x;
    c5_bb_y = muDoubleScalarAbs(c5_qb_x);
    c5_rb_x = c5_salpha.im;
    c5_sb_x = c5_rb_x;
    c5_cb_y = muDoubleScalarAbs(c5_sb_x);
    c5_bcoefa = c5_bb_y + c5_cb_y;
    for (c5_jr = 0; c5_jr < 3; c5_jr++) {
      c5_b_jr = 1.0 + (real_T)c5_jr;
      c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_jr), 1, 3, 1, 0) - 1].re = c5_dc11.re;
      c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_jr), 1, 3, 1, 0) - 1].im = c5_dc11.im;
    }

    c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c5_b_je), 1, 3, 1, 0) - 1].re = c5_dc12.re;
    c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c5_b_je), 1, 3, 1, 0) - 1].im = c5_dc12.im;
    c5_tb_x = 2.2204460492503131E-16 * c5_acoefa * c5_anorm;
    c5_db_y = 2.2204460492503131E-16 * c5_bcoefa;
    c5_dmin = c5_tb_x;
    if (c5_db_y > c5_dmin) {
      c5_dmin = c5_db_y;
    }

    if (2.2250738585072014E-308 > c5_dmin) {
      c5_dmin = 2.2250738585072014E-308;
    }

    c5_d10 = c5_b_je - 1.0;
    c5_i144 = (int32_T)c5_d10;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c5_d10, mxDOUBLE_CLASS, c5_i144);
    for (c5_c_jr = 0; c5_c_jr < c5_i144; c5_c_jr++) {
      c5_b_jr = 1.0 + (real_T)c5_c_jr;
      c5_d_a = c5_acoeff;
      c5_ca.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_je), 1, 3, 2, 0) - 1)) - 1].re;
      c5_ca.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_je), 1, 3, 2, 0) - 1)) - 1].im;
      c5_ca.re *= c5_d_a;
      c5_ca.im *= c5_d_a;
      c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_jr), 1, 3, 1, 0) - 1].re = c5_ca.re;
      c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_jr), 1, 3, 1, 0) - 1].im = c5_ca.im;
    }

    c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c5_b_je), 1, 3, 1, 0) - 1].re = c5_dc12.re;
    c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c5_b_je), 1, 3, 1, 0) - 1].im = c5_dc12.im;
    c5_d11 = c5_b_je - 1.0;
    c5_i145 = (int32_T)-(1.0 + (-1.0 - c5_d11));
    _SFD_FOR_LOOP_VECTOR_CHECK(c5_d11, -1.0, 1.0, mxDOUBLE_CLASS, c5_i145);
    for (c5_c_j = 0; c5_c_j < c5_i145; c5_c_j++) {
      c5_b_j = c5_d11 + -(real_T)c5_c_j;
      c5_e_a = c5_acoeff;
      c5_ca.re = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c5_ca.im = c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].im;
      c5_ca.re *= c5_e_a;
      c5_ca.im *= c5_e_a;
      c5_d.re = c5_ca.re - c5_salpha.re;
      c5_d.im = c5_ca.im - c5_salpha.im;
      c5_ub_x = c5_d.re;
      c5_vb_x = c5_ub_x;
      c5_eb_y = muDoubleScalarAbs(c5_vb_x);
      c5_wb_x = c5_d.im;
      c5_xb_x = c5_wb_x;
      c5_fb_y = muDoubleScalarAbs(c5_xb_x);
      c5_gb_y = c5_eb_y + c5_fb_y;
      if (c5_gb_y <= c5_dmin) {
        c5_d.re = c5_dmin;
        c5_d.im = 0.0;
      }

      c5_yb_x = c5_d.re;
      c5_ac_x = c5_yb_x;
      c5_hb_y = muDoubleScalarAbs(c5_ac_x);
      c5_bc_x = c5_d.im;
      c5_cc_x = c5_bc_x;
      c5_ib_y = muDoubleScalarAbs(c5_cc_x);
      c5_jb_y = c5_hb_y + c5_ib_y;
      if (c5_jb_y < 1.0) {
        c5_ca.re = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].re;
        c5_ca.im = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].im;
        c5_dc_x = c5_ca.re;
        c5_ec_x = c5_dc_x;
        c5_kb_y = muDoubleScalarAbs(c5_ec_x);
        c5_fc_x = c5_ca.im;
        c5_gc_x = c5_fc_x;
        c5_lb_y = muDoubleScalarAbs(c5_gc_x);
        c5_mb_y = c5_kb_y + c5_lb_y;
        c5_hc_x = c5_d.re;
        c5_ic_x = c5_hc_x;
        c5_nb_y = muDoubleScalarAbs(c5_ic_x);
        c5_jc_x = c5_d.im;
        c5_kc_x = c5_jc_x;
        c5_ob_y = muDoubleScalarAbs(c5_kc_x);
        c5_pb_y = c5_nb_y + c5_ob_y;
        if (c5_mb_y >= 1.4980776123852632E+307 * c5_pb_y) {
          c5_ca.re = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].re;
          c5_ca.im = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].im;
          c5_lc_x = c5_ca.re;
          c5_mc_x = c5_lc_x;
          c5_qb_y = muDoubleScalarAbs(c5_mc_x);
          c5_nc_x = c5_ca.im;
          c5_oc_x = c5_nc_x;
          c5_rb_y = muDoubleScalarAbs(c5_oc_x);
          c5_sb_y = c5_qb_y + c5_rb_y;
          c5_temp = 1.0 / c5_sb_y;
          c5_c_je = c5_b_je;
          c5_i146 = (int32_T)c5_c_je;
          _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c5_c_je, mxDOUBLE_CLASS, c5_i146);
          for (c5_d_jr = 0; c5_d_jr < c5_i146; c5_d_jr++) {
            c5_b_jr = 1.0 + (real_T)c5_d_jr;
            c5_f_a = c5_temp;
            c5_ca.re = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].re;
            c5_ca.im = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].im;
            c5_ca.re *= c5_f_a;
            c5_ca.im *= c5_f_a;
            c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
              ("", c5_b_jr), 1, 3, 1, 0) - 1].re = c5_ca.re;
            c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
              ("", c5_b_jr), 1, 3, 1, 0) - 1].im = c5_ca.im;
          }
        }
      }

      c5_ca.re = -c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].re;
      c5_ca.im = -c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].im;
      c5_ca = c5_rdivide(chartInstance, c5_ca, c5_d);
      c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_j), 1, 3, 1, 0) - 1].re = c5_ca.re;
      c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_j), 1, 3, 1, 0) - 1].im = c5_ca.im;
      if (c5_b_j > 1.0) {
        c5_ca.re = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].re;
        c5_ca.im = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].im;
        c5_pc_x = c5_ca.re;
        c5_qc_x = c5_pc_x;
        c5_tb_y = muDoubleScalarAbs(c5_qc_x);
        c5_rc_x = c5_ca.im;
        c5_sc_x = c5_rc_x;
        c5_ub_y = muDoubleScalarAbs(c5_sc_x);
        c5_vb_y = c5_tb_y + c5_ub_y;
        if (c5_vb_y > 1.0) {
          c5_ca.re = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].re;
          c5_ca.im = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].im;
          c5_tc_x = c5_ca.re;
          c5_uc_x = c5_tc_x;
          c5_wb_y = muDoubleScalarAbs(c5_uc_x);
          c5_vc_x = c5_ca.im;
          c5_wc_x = c5_vc_x;
          c5_xb_y = muDoubleScalarAbs(c5_wc_x);
          c5_yb_y = c5_wb_y + c5_xb_y;
          c5_temp = 1.0 / c5_yb_y;
          if (c5_acoefa * c5_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
               _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1] >=
              1.4980776123852632E+307 * c5_temp) {
            c5_d_je = c5_b_je;
            c5_i147 = (int32_T)c5_d_je;
            _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c5_d_je, mxDOUBLE_CLASS,
              c5_i147);
            for (c5_e_jr = 0; c5_e_jr < c5_i147; c5_e_jr++) {
              c5_b_jr = 1.0 + (real_T)c5_e_jr;
              c5_g_a = c5_temp;
              c5_ca.re = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].re;
              c5_ca.im = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].im;
              c5_ca.re *= c5_g_a;
              c5_ca.im *= c5_g_a;
              c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].re = c5_ca.re;
              c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].im = c5_ca.im;
            }
          }
        }

        c5_h_a = c5_acoeff;
        c5_ca.re = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].re;
        c5_ca.im = c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 1, 0) - 1].im;
        c5_ca.re *= c5_h_a;
        c5_ca.im *= c5_h_a;
        c5_d12 = c5_b_j - 1.0;
        c5_i148 = (int32_T)c5_d12;
        _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c5_d12, mxDOUBLE_CLASS, c5_i148);
        for (c5_f_jr = 0; c5_f_jr < c5_i148; c5_f_jr++) {
          c5_b_jr = 1.0 + (real_T)c5_f_jr;
          c5_b_ca.re = c5_ca.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            c5_b_j), 1, 3, 2, 0) - 1)) - 1].re - c5_ca.im * c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c5_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].
            im;
          c5_b_ca.im = c5_ca.re * c5_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            c5_b_j), 1, 3, 2, 0) - 1)) - 1].im + c5_ca.im * c5_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c5_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c5_b_j), 1, 3, 2, 0) - 1)) - 1].
            re;
          c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
            "", c5_b_jr), 1, 3, 1, 0) - 1].re =
            c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c5_b_jr), 1, 3, 1, 0) - 1].re + c5_b_ca.re;
          c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
            "", c5_b_jr), 1, 3, 1, 0) - 1].im =
            c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c5_b_jr), 1, 3, 1, 0) - 1].im + c5_b_ca.im;
        }
      }
    }

    for (c5_g_jr = 0; c5_g_jr < 3; c5_g_jr++) {
      c5_b_jr = 1.0 + (real_T)c5_g_jr;
      c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_jr), 1, 3, 1, 0) - 1].re = c5_dc11.re;
      c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c5_b_jr), 1, 3, 1, 0) - 1].im = c5_dc11.im;
    }

    c5_e_je = c5_b_je;
    c5_i149 = (int32_T)c5_e_je;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c5_e_je, mxDOUBLE_CLASS, c5_i149);
    for (c5_jc = 0; c5_jc < c5_i149; c5_jc++) {
      c5_b_jc = 1.0 + (real_T)c5_jc;
      for (c5_h_jr = 0; c5_h_jr < 3; c5_h_jr++) {
        c5_b_jr = 1.0 + (real_T)c5_h_jr;
        c5_b_V.re = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1)) - 1].re *
          c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c5_b_jc), 1, 3, 1, 0) - 1].re - c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c5_b_jc), 1, 3, 2, 0) - 1)) - 1].im *
          c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c5_b_jc), 1, 3, 1, 0) - 1].im;
        c5_b_V.im = c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_jc), 1, 3, 2, 0) - 1)) - 1].re *
          c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c5_b_jc), 1, 3, 1, 0) - 1].im + c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) + 3 *
          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c5_b_jc), 1, 3, 2, 0) - 1)) - 1].im *
          c5_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c5_b_jc), 1, 3, 1, 0) - 1].re;
        c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c5_b_jr), 1, 3, 1, 0) - 1].re = c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].re +
          c5_b_V.re;
        c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c5_b_jr), 1, 3, 1, 0) - 1].im = c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].im +
          c5_b_V.im;
      }
    }

    c5_ca = c5_work2[0];
    c5_xc_x = c5_ca.re;
    c5_yc_x = c5_xc_x;
    c5_ac_y = muDoubleScalarAbs(c5_yc_x);
    c5_ad_x = c5_ca.im;
    c5_bd_x = c5_ad_x;
    c5_bc_y = muDoubleScalarAbs(c5_bd_x);
    c5_xmx = c5_ac_y + c5_bc_y;
    for (c5_i_jr = 0; c5_i_jr < 2; c5_i_jr++) {
      c5_b_jr = 2.0 + (real_T)c5_i_jr;
      c5_ca.re = c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].re;
      c5_ca.im = c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].im;
      c5_cd_x = c5_ca.re;
      c5_dd_x = c5_cd_x;
      c5_cc_y = muDoubleScalarAbs(c5_dd_x);
      c5_ed_x = c5_ca.im;
      c5_fd_x = c5_ed_x;
      c5_dc_y = muDoubleScalarAbs(c5_fd_x);
      c5_ec_y = c5_cc_y + c5_dc_y;
      c5_fc_y = c5_ec_y;
      if (c5_fc_y > c5_xmx) {
        c5_xmx = c5_fc_y;
      }
    }

    if (c5_xmx > 2.2250738585072014E-308) {
      c5_temp = 1.0 / c5_xmx;
      for (c5_j_jr = 0; c5_j_jr < 3; c5_j_jr++) {
        c5_b_jr = 1.0 + (real_T)c5_j_jr;
        c5_i_a = c5_temp;
        c5_ca.re = c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].re;
        c5_ca.im = c5_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c5_b_jr), 1, 3, 1, 0) - 1].im;
        c5_ca.re *= c5_i_a;
        c5_ca.im *= c5_i_a;
        c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c5_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c5_ieig), 1, 3, 2, 0) - 1)) - 1]
          .re = c5_ca.re;
        c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c5_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c5_ieig), 1, 3, 2, 0) - 1)) - 1]
          .im = c5_ca.im;
      }
    } else {
      for (c5_k_jr = 0; c5_k_jr < 3; c5_k_jr++) {
        c5_b_jr = 1.0 + (real_T)c5_k_jr;
        c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c5_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c5_ieig), 1, 3, 2, 0) - 1)) - 1]
          .re = c5_dc11.re;
        c5_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c5_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c5_ieig), 1, 3, 2, 0) - 1)) - 1]
          .im = c5_dc11.im;
      }
    }
  }
}

static int32_T c5_div_s32(SFc5_Model_01InstanceStruct *chartInstance, int32_T
  c5_numerator, int32_T c5_denominator)
{
  int32_T c5_quotient;
  uint32_T c5_absNumerator;
  uint32_T c5_absDenominator;
  boolean_T c5_quotientNeedsNegation;
  uint32_T c5_tempAbsQuotient;
  (void)chartInstance;
  if (c5_denominator == 0) {
    if (c5_numerator >= 0) {
      c5_quotient = MAX_int32_T;
    } else {
      c5_quotient = MIN_int32_T;
    }

    _SFD_OVERFLOW_DETECTION(SFDB_DIVIDE_BY_ZERO);
  } else {
    if (c5_numerator >= 0) {
      c5_absNumerator = (uint32_T)c5_numerator;
    } else {
      c5_absNumerator = (uint32_T)-c5_numerator;
    }

    if (c5_denominator >= 0) {
      c5_absDenominator = (uint32_T)c5_denominator;
    } else {
      c5_absDenominator = (uint32_T)-c5_denominator;
    }

    c5_quotientNeedsNegation = (c5_numerator < 0 != c5_denominator < 0);
    c5_tempAbsQuotient = c5_absNumerator / c5_absDenominator;
    if (c5_quotientNeedsNegation) {
      c5_quotient = -(int32_T)c5_tempAbsQuotient;
    } else {
      c5_quotient = (int32_T)c5_tempAbsQuotient;
    }
  }

  return c5_quotient;
}

static void init_dsm_address_info(SFc5_Model_01InstanceStruct *chartInstance)
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

void sf_c5_Model_01_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(705285550U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(968497424U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3109536952U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3196241693U);
}

mxArray *sf_c5_Model_01_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("1SXmiNdsrSH20fPu4M4XFB");
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
      pr[0] = (double)(3);
      pr[1] = (double)(3);
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
      pr[1] = (double)(3);
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
      pr[1] = (double)(3);
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

mxArray *sf_c5_Model_01_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c5_Model_01_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c5_Model_01(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"Phi_r23\",},{M[8],M[0],T\"is_active_c5_Model_01\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c5_Model_01_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc5_Model_01InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc5_Model_01InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_01MachineNumber_,
           5,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"omega");
          _SFD_SET_DATA_PROPS(1,2,0,1,"Phi_r23");
          _SFD_SET_DATA_PROPS(2,1,1,0,"t_delta");
          _SFD_SET_DATA_PROPS(3,1,1,0,"M");
          _SFD_SET_DATA_PROPS(4,1,1,0,"N");
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
        _SFD_CV_INIT_EML(0,1,5,0,0,0,1,3,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,556);
        _SFD_CV_INIT_EML_FCN(0,1,"fn_phi_jk",557,-1,896);
        _SFD_CV_INIT_EML_FCN(0,2,"fn_phi_j1",898,-1,1018);
        _SFD_CV_INIT_EML_FCN(0,3,"fn_phi_j2",1020,-1,1299);
        _SFD_CV_INIT_EML_FCN(0,4,"fn_phi_j3",1305,-1,1654);
        _SFD_CV_INIT_EML_FOR(0,1,0,337,348,552);
        _SFD_CV_INIT_EML_FOR(0,1,1,356,368,544);
        _SFD_CV_INIT_EML_FOR(0,1,2,380,392,532);

        {
          static int caseStart[] = { 851, 636, 700, 775 };

          static int caseExprEnd[] = { 860, 642, 706, 781 };

          _SFD_CV_INIT_EML_SWITCH(0,1,0,618,627,892,4,&(caseStart[0]),
            &(caseExprEnd[0]));
        }

        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"fn_VectorToSkewSymmetricTensor",0,-1,433);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c5_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)
            c5_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c5_t_delta;
          real_T (*c5_omega)[3];
          real_T (*c5_Phi_r23)[9];
          real_T (*c5_M)[9];
          real_T (*c5_N)[9];
          c5_N = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
          c5_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 2);
          c5_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c5_Phi_r23 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
          c5_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c5_omega);
          _SFD_SET_DATA_VALUE_PTR(1U, *c5_Phi_r23);
          _SFD_SET_DATA_VALUE_PTR(2U, c5_t_delta);
          _SFD_SET_DATA_VALUE_PTR(3U, *c5_M);
          _SFD_SET_DATA_VALUE_PTR(4U, *c5_N);
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
  return "iX9Zl6FlNLqtZng6BCwrjF";
}

static void sf_opaque_initialize_c5_Model_01(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc5_Model_01InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c5_Model_01((SFc5_Model_01InstanceStruct*) chartInstanceVar);
  initialize_c5_Model_01((SFc5_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c5_Model_01(void *chartInstanceVar)
{
  enable_c5_Model_01((SFc5_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c5_Model_01(void *chartInstanceVar)
{
  disable_c5_Model_01((SFc5_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c5_Model_01(void *chartInstanceVar)
{
  sf_gateway_c5_Model_01((SFc5_Model_01InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c5_Model_01(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c5_Model_01((SFc5_Model_01InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c5_Model_01();/* state var info */
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

extern void sf_internal_set_sim_state_c5_Model_01(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c5_Model_01();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c5_Model_01((SFc5_Model_01InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c5_Model_01(SimStruct* S)
{
  return sf_internal_get_sim_state_c5_Model_01(S);
}

static void sf_opaque_set_sim_state_c5_Model_01(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c5_Model_01(S, st);
}

static void sf_opaque_terminate_c5_Model_01(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc5_Model_01InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_01_optimization_info();
    }

    finalize_c5_Model_01((SFc5_Model_01InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc5_Model_01((SFc5_Model_01InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c5_Model_01(SimStruct *S)
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
    initialize_params_c5_Model_01((SFc5_Model_01InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c5_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_01_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,5);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,5,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,5,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,5);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,5,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,5,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,5);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3136522695U));
  ssSetChecksum1(S,(2338157068U));
  ssSetChecksum2(S,(1113484602U));
  ssSetChecksum3(S,(1360097174U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c5_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c5_Model_01(SimStruct *S)
{
  SFc5_Model_01InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc5_Model_01InstanceStruct *)utMalloc(sizeof
    (SFc5_Model_01InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc5_Model_01InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c5_Model_01;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c5_Model_01;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c5_Model_01;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c5_Model_01;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c5_Model_01;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c5_Model_01;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c5_Model_01;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c5_Model_01;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c5_Model_01;
  chartInstance->chartInfo.mdlStart = mdlStart_c5_Model_01;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c5_Model_01;
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

void c5_Model_01_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c5_Model_01(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c5_Model_01(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c5_Model_01(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c5_Model_01_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

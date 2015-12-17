/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_01_sfun.h"
#include "c1_Model_01.h"
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
static const char * c1_debug_family_names[14] = { "v_lambda", "A", "Gamma",
  "Omega_Tensor", "omega_norm", "l", "j", "k", "nargin", "nargout", "omega", "M",
  "t_delta", "Phi_r12" };

static const char * c1_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

static const char * c1_c_debug_family_names[6] = { "nargin", "nargout", "lambda",
  "omega_norm", "t_delta", "phi_j3" };

static const char * c1_d_debug_family_names[6] = { "nargin", "nargout", "lambda",
  "omega_norm", "t_delta", "phi_j2" };

static const char * c1_e_debug_family_names[5] = { "nargin", "nargout", "lambda",
  "t_delta", "phi_j1" };

static const char * c1_f_debug_family_names[7] = { "nargin", "nargout", "k",
  "lambda", "t_delta", "omega_norm", "phi_jk" };

/* Function Declarations */
static void initialize_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance);
static void initialize_params_c1_Model_01(SFc1_Model_01InstanceStruct
  *chartInstance);
static void enable_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance);
static void disable_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_Model_01(SFc1_Model_01InstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c1_Model_01(SFc1_Model_01InstanceStruct
  *chartInstance);
static void set_sim_state_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_st);
static void finalize_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance);
static void sf_gateway_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_chartstep_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance);
static void initSimStructsc1_Model_01(SFc1_Model_01InstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static void c1_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_Phi_r12, const char_T *c1_identifier, real_T c1_y[9]);
static void c1_b_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[9]);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static real_T c1_c_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_d_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3]);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_info_helper(const mxArray **c1_info);
static const mxArray *c1_emlrt_marshallOut(const char * c1_u);
static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u);
static void c1_b_info_helper(const mxArray **c1_info);
static void c1_c_info_helper(const mxArray **c1_info);
static void c1_d_info_helper(const mxArray **c1_info);
static void c1_e_info_helper(const mxArray **c1_info);
static void c1_eig(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_A[9],
                   creal_T c1_V[3]);
static void c1_realmin(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_eml_error(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_eps(SFc1_Model_01InstanceStruct *chartInstance);
static real_T c1_eml_matlab_zlangeM(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_x[9]);
static real_T c1_abs(SFc1_Model_01InstanceStruct *chartInstance, creal_T c1_x);
static boolean_T c1_isfinite(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_x);
static void c1_eml_matlab_zlascl(SFc1_Model_01InstanceStruct *chartInstance,
  real_T c1_cfrom, real_T c1_cto, creal_T c1_A[9], creal_T c1_b_A[9]);
static real_T c1_b_abs(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_x);
static void c1_eml_matlab_zggbal(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], creal_T c1_b_A[9], int32_T *c1_ilo, int32_T *c1_ihi, int32_T
  c1_rscale[3]);
static void c1_eml_switch_helper(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_check_forloop_overflow_error(SFc1_Model_01InstanceStruct
  *chartInstance, boolean_T c1_overflow);
static void c1_eml_matlab_zlartg(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_f, creal_T c1_g, real_T *c1_cs, creal_T *c1_sn, creal_T *c1_r);
static void c1_eml_scalar_eg(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_eml_matlab_zhgeqz(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T c1_ilo, int32_T c1_ihi, real_T *c1_info, creal_T
  c1_alpha1[3], creal_T c1_beta1[3]);
static real_T c1_eml_matlab_zlanhs(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T c1_ilo, int32_T c1_ihi);
static int32_T c1_mod(SFc1_Model_01InstanceStruct *chartInstance, int32_T c1_x);
static creal_T c1_eml_div(SFc1_Model_01InstanceStruct *chartInstance, creal_T
  c1_x, real_T c1_y);
static void c1_scalarEg(SFc1_Model_01InstanceStruct *chartInstance);
static creal_T c1_sqrt(SFc1_Model_01InstanceStruct *chartInstance, creal_T c1_x);
static void c1_realmax(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_b_eml_matlab_zlartg(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_f, creal_T c1_g, real_T *c1_cs, creal_T *c1_sn);
static void c1_b_eml_matlab_zlascl(SFc1_Model_01InstanceStruct *chartInstance,
  real_T c1_cfrom, real_T c1_cto, creal_T c1_A[3], creal_T c1_b_A[3]);
static void c1_b_eml_div(SFc1_Model_01InstanceStruct *chartInstance, creal_T
  c1_x[3], creal_T c1_y[3], creal_T c1_z[3]);
static void c1_eml_warning(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_b_eml_warning(SFc1_Model_01InstanceStruct *chartInstance);
static real_T c1_mpower(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_a);
static void c1_inv(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_x[9],
                   real_T c1_y[9]);
static void c1_inv3x3(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_x[9],
                      real_T c1_y[9]);
static real_T c1_norm(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_x[9]);
static void c1_c_eml_warning(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_d_eml_warning(SFc1_Model_01InstanceStruct *chartInstance, char_T
  c1_varargin_2[14]);
static void c1_b_eml_switch_helper(SFc1_Model_01InstanceStruct *chartInstance);
static real_T c1_fn_phi_jk(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_k, real_T c1_lambda, real_T c1_t_delta, real_T c1_omega_norm);
static real_T c1_b_mpower(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_a);
static real_T c1_c_mpower(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_a);
static void c1_d_mpower(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_a
  [9], real_T c1_b, real_T c1_c[9]);
static void c1_b_eml_scalar_eg(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_c_eml_scalar_eg(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_eml_xgemm(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_A[9], real_T c1_B[9], real_T c1_C[9], real_T c1_b_C[9]);
static void c1_b_eml_error(SFc1_Model_01InstanceStruct *chartInstance);
static void c1_b_eig(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_A[9],
                     creal_T c1_V[9], creal_T c1_D[9]);
static void c1_eml_matlab_zggev(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], real_T *c1_info, creal_T c1_alpha1[3], creal_T c1_beta1[3],
  creal_T c1_V[9]);
static void c1_eml_matlab_zgghrd(SFc1_Model_01InstanceStruct *chartInstance,
  int32_T c1_ilo, int32_T c1_ihi, creal_T c1_A[9], creal_T c1_b_A[9], creal_T
  c1_Z[9]);
static void c1_b_eml_matlab_zhgeqz(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T c1_ilo, int32_T c1_ihi, creal_T c1_Z[9], real_T
  *c1_info, creal_T c1_alpha1[3], creal_T c1_beta1[3], creal_T c1_b_A[9],
  creal_T c1_b_Z[9]);
static void c1_eml_matlab_ztgevc(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], creal_T c1_V[9], creal_T c1_b_V[9]);
static creal_T c1_rdivide(SFc1_Model_01InstanceStruct *chartInstance, creal_T
  c1_x, creal_T c1_y);
static creal_T c1_power(SFc1_Model_01InstanceStruct *chartInstance, creal_T c1_a,
  real_T c1_b);
static void c1_eml_lusolve(SFc1_Model_01InstanceStruct *chartInstance, creal_T
  c1_A[9], creal_T c1_B[9], creal_T c1_X[9]);
static void c1_e_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_sprintf, const char_T *c1_identifier, char_T c1_y[14]);
static void c1_f_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, char_T c1_y[14]);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_g_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_h_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_Model_01, const char_T *c1_identifier);
static uint8_T c1_i_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_c_eml_matlab_zlascl(SFc1_Model_01InstanceStruct *chartInstance,
  real_T c1_cfrom, real_T c1_cto, creal_T c1_A[9]);
static void c1_b_eml_matlab_zggbal(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T *c1_ilo, int32_T *c1_ihi, int32_T c1_rscale[3]);
static void c1_b_sqrt(SFc1_Model_01InstanceStruct *chartInstance, creal_T *c1_x);
static void c1_d_eml_matlab_zlascl(SFc1_Model_01InstanceStruct *chartInstance,
  real_T c1_cfrom, real_T c1_cto, creal_T c1_A[3]);
static void c1_b_eml_xgemm(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_A[9], real_T c1_B[9], real_T c1_C[9]);
static void c1_b_eml_matlab_zgghrd(SFc1_Model_01InstanceStruct *chartInstance,
  int32_T c1_ilo, int32_T c1_ihi, creal_T c1_A[9], creal_T c1_Z[9]);
static void c1_c_eml_matlab_zhgeqz(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T c1_ilo, int32_T c1_ihi, creal_T c1_Z[9], real_T
  *c1_info, creal_T c1_alpha1[3], creal_T c1_beta1[3]);
static void c1_b_eml_matlab_ztgevc(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], creal_T c1_V[9]);
static int32_T c1_div_s32(SFc1_Model_01InstanceStruct *chartInstance, int32_T
  c1_numerator, int32_T c1_denominator);
static void init_dsm_address_info(SFc1_Model_01InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance)
{
  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c1_is_active_c1_Model_01 = 0U;
}

static void initialize_params_c1_Model_01(SFc1_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c1_update_debugger_state_c1_Model_01(SFc1_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c1_Model_01(SFc1_Model_01InstanceStruct
  *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  int32_T c1_i0;
  real_T c1_u[9];
  const mxArray *c1_b_y = NULL;
  uint8_T c1_hoistedGlobal;
  uint8_T c1_b_u;
  const mxArray *c1_c_y = NULL;
  real_T (*c1_Phi_r12)[9];
  c1_Phi_r12 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c1_st = NULL;
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellmatrix(2, 1), false);
  for (c1_i0 = 0; c1_i0 < 9; c1_i0++) {
    c1_u[c1_i0] = (*c1_Phi_r12)[c1_i0];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  c1_hoistedGlobal = chartInstance->c1_is_active_c1_Model_01;
  c1_b_u = c1_hoistedGlobal;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  sf_mex_assign(&c1_st, c1_y, false);
  return c1_st;
}

static void set_sim_state_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_st)
{
  const mxArray *c1_u;
  real_T c1_dv0[9];
  int32_T c1_i1;
  real_T (*c1_Phi_r12)[9];
  c1_Phi_r12 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c1_doneDoubleBufferReInit = true;
  c1_u = sf_mex_dup(c1_st);
  c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 0)),
                      "Phi_r12", c1_dv0);
  for (c1_i1 = 0; c1_i1 < 9; c1_i1++) {
    (*c1_Phi_r12)[c1_i1] = c1_dv0[c1_i1];
  }

  chartInstance->c1_is_active_c1_Model_01 = c1_h_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c1_u, 1)), "is_active_c1_Model_01");
  sf_mex_destroy(&c1_u);
  c1_update_debugger_state_c1_Model_01(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance)
{
  int32_T c1_i2;
  int32_T c1_i3;
  int32_T c1_i4;
  real_T *c1_t_delta;
  real_T (*c1_M)[9];
  real_T (*c1_Phi_r12)[9];
  real_T (*c1_omega)[3];
  c1_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c1_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 1);
  c1_Phi_r12 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c1_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  for (c1_i2 = 0; c1_i2 < 3; c1_i2++) {
    _SFD_DATA_RANGE_CHECK((*c1_omega)[c1_i2], 0U);
  }

  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_Model_01(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_01MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c1_i3 = 0; c1_i3 < 9; c1_i3++) {
    _SFD_DATA_RANGE_CHECK((*c1_Phi_r12)[c1_i3], 1U);
  }

  for (c1_i4 = 0; c1_i4 < 9; c1_i4++) {
    _SFD_DATA_RANGE_CHECK((*c1_M)[c1_i4], 2U);
  }

  _SFD_DATA_RANGE_CHECK(*c1_t_delta, 3U);
}

static void c1_chartstep_c1_Model_01(SFc1_Model_01InstanceStruct *chartInstance)
{
  real_T c1_hoistedGlobal;
  int32_T c1_i5;
  real_T c1_omega[3];
  int32_T c1_i6;
  real_T c1_M[9];
  real_T c1_t_delta;
  uint32_T c1_debug_family_var_map[14];
  real_T c1_v_lambda[3];
  real_T c1_A[9];
  real_T c1_Gamma[9];
  real_T c1_Omega_Tensor[9];
  real_T c1_omega_norm;
  real_T c1_l;
  real_T c1_j;
  real_T c1_k;
  real_T c1_nargin = 3.0;
  real_T c1_nargout = 1.0;
  real_T c1_Phi_r12[9];
  int32_T c1_i7;
  int32_T c1_i8;
  real_T c1_b_M[9];
  creal_T c1_dcv0[3];
  int32_T c1_i9;
  real_T c1_d0;
  real_T c1_d1;
  real_T c1_d2;
  int32_T c1_i10;
  real_T c1_b_A[9];
  real_T c1_dv1[9];
  int32_T c1_i11;
  int32_T c1_i12;
  real_T c1_v[3];
  uint32_T c1_b_debug_family_var_map[4];
  real_T c1_b_nargin = 1.0;
  real_T c1_b_nargout = 1.0;
  int32_T c1_i13;
  int32_T c1_i14;
  real_T c1_x[3];
  real_T c1_scale;
  int32_T c1_b_k;
  int32_T c1_c_k;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_absxk;
  real_T c1_t;
  int32_T c1_b_l;
  int32_T c1_b_j;
  int32_T c1_d_k;
  real_T c1_a;
  int32_T c1_i15;
  real_T c1_b_Omega_Tensor[9];
  real_T c1_b[9];
  int32_T c1_i16;
  int32_T c1_i17;
  real_T c1_c_M[9];
  real_T c1_b_b[9];
  int32_T c1_i18;
  real_T c1_y[9];
  int32_T c1_i19;
  real_T c1_c_b[9];
  int32_T c1_i20;
  real_T c1_d_b[9];
  int32_T c1_i21;
  int32_T c1_i22;
  int32_T c1_i23;
  int32_T c1_i24;
  real_T *c1_b_t_delta;
  real_T (*c1_b_Phi_r12)[9];
  real_T (*c1_d_M)[9];
  real_T (*c1_b_omega)[3];
  c1_b_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c1_d_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 1);
  c1_b_Phi_r12 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c1_b_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  c1_hoistedGlobal = *c1_b_t_delta;
  for (c1_i5 = 0; c1_i5 < 3; c1_i5++) {
    c1_omega[c1_i5] = (*c1_b_omega)[c1_i5];
  }

  for (c1_i6 = 0; c1_i6 < 9; c1_i6++) {
    c1_M[c1_i6] = (*c1_d_M)[c1_i6];
  }

  c1_t_delta = c1_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 14U, 14U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_v_lambda, 0U, c1_c_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_A, 1U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_Gamma, 2U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_Omega_Tensor, 3U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_omega_norm, 4U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_l, 5U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_j, 6U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_k, 7U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 8U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 9U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_omega, 10U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_M, 11U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_t_delta, 12U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_Phi_r12, 13U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 3);
  for (c1_i7 = 0; c1_i7 < 9; c1_i7++) {
    c1_Phi_r12[c1_i7] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 4);
  for (c1_i8 = 0; c1_i8 < 9; c1_i8++) {
    c1_b_M[c1_i8] = c1_M[c1_i8];
  }

  c1_eig(chartInstance, c1_b_M, c1_dcv0);
  for (c1_i9 = 0; c1_i9 < 3; c1_i9++) {
    c1_v_lambda[c1_i9] = c1_dcv0[c1_i9].re;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 6);
  c1_d0 = c1_mpower(chartInstance, c1_v_lambda[0]);
  c1_d1 = c1_mpower(chartInstance, c1_v_lambda[1]);
  c1_d2 = c1_mpower(chartInstance, c1_v_lambda[2]);
  c1_A[0] = 1.0;
  c1_A[3] = c1_v_lambda[0];
  c1_A[6] = c1_d0;
  c1_A[1] = 1.0;
  c1_A[4] = c1_v_lambda[1];
  c1_A[7] = c1_d1;
  c1_A[2] = 1.0;
  c1_A[5] = c1_v_lambda[2];
  c1_A[8] = c1_d2;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 7);
  for (c1_i10 = 0; c1_i10 < 9; c1_i10++) {
    c1_b_A[c1_i10] = c1_A[c1_i10];
  }

  c1_inv(chartInstance, c1_b_A, c1_dv1);
  for (c1_i11 = 0; c1_i11 < 9; c1_i11++) {
    c1_Gamma[c1_i11] = c1_dv1[c1_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 8);
  for (c1_i12 = 0; c1_i12 < 3; c1_i12++) {
    c1_v[c1_i12] = c1_omega[c1_i12];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c1_b_debug_family_names,
    c1_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargin, 0U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargout, 1U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_v, 2U, c1_c_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_Omega_Tensor, 3U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 2);
  for (c1_i13 = 0; c1_i13 < 9; c1_i13++) {
    c1_Omega_Tensor[c1_i13] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 3);
  c1_Omega_Tensor[0] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 4);
  c1_Omega_Tensor[4] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 5);
  c1_Omega_Tensor[8] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 6);
  c1_Omega_Tensor[3] = -c1_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 7);
  c1_Omega_Tensor[6] = c1_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 8);
  c1_Omega_Tensor[7] = -c1_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 9);
  c1_Omega_Tensor[1] = c1_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 10);
  c1_Omega_Tensor[2] = -c1_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 11);
  c1_Omega_Tensor[5] = c1_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 9);
  for (c1_i14 = 0; c1_i14 < 3; c1_i14++) {
    c1_x[c1_i14] = c1_omega[c1_i14];
  }

  c1_b_eml_switch_helper(chartInstance);
  c1_omega_norm = 0.0;
  c1_realmin(chartInstance);
  c1_scale = 2.2250738585072014E-308;
  for (c1_b_k = 1; c1_b_k < 4; c1_b_k++) {
    c1_c_k = c1_b_k;
    c1_b_x = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_c_k), 1, 3, 1, 0) - 1];
    c1_c_x = c1_b_x;
    c1_absxk = muDoubleScalarAbs(c1_c_x);
    if (c1_absxk > c1_scale) {
      c1_t = c1_scale / c1_absxk;
      c1_omega_norm = 1.0 + c1_omega_norm * c1_t * c1_t;
      c1_scale = c1_absxk;
    } else {
      c1_t = c1_absxk / c1_scale;
      c1_omega_norm += c1_t * c1_t;
    }
  }

  c1_omega_norm = c1_scale * muDoubleScalarSqrt(c1_omega_norm);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 10);
  c1_l = 1.0;
  c1_b_l = 0;
  while (c1_b_l < 3) {
    c1_l = 1.0 + (real_T)c1_b_l;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 11);
    c1_j = 1.0;
    c1_b_j = 0;
    while (c1_b_j < 3) {
      c1_j = 1.0 + (real_T)c1_b_j;
      CV_EML_FOR(0, 1, 1, 1);
      _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 12);
      c1_k = 1.0;
      c1_d_k = 0;
      while (c1_d_k < 3) {
        c1_k = 1.0 + (real_T)c1_d_k;
        CV_EML_FOR(0, 1, 2, 1);
        _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 13);
        c1_a = c1_Gamma[(_SFD_EML_ARRAY_BOUNDS_CHECK("Gamma", (int32_T)
          _SFD_INTEGER_CHECK("l", c1_l), 1, 3, 1, 0) + 3 *
                         (_SFD_EML_ARRAY_BOUNDS_CHECK("Gamma", (int32_T)
          _SFD_INTEGER_CHECK("j", c1_j), 1, 3, 2, 0) - 1)) - 1] * c1_fn_phi_jk
          (chartInstance, c1_k, c1_v_lambda[_SFD_EML_ARRAY_BOUNDS_CHECK(
            "v_lambda", (int32_T)_SFD_INTEGER_CHECK("j", c1_j), 1, 3, 1, 0) - 1],
           c1_t_delta, c1_omega_norm);
        for (c1_i15 = 0; c1_i15 < 9; c1_i15++) {
          c1_b_Omega_Tensor[c1_i15] = c1_Omega_Tensor[c1_i15];
        }

        c1_d_mpower(chartInstance, c1_b_Omega_Tensor, c1_k - 1.0, c1_b);
        for (c1_i16 = 0; c1_i16 < 9; c1_i16++) {
          c1_b[c1_i16] *= c1_a;
        }

        for (c1_i17 = 0; c1_i17 < 9; c1_i17++) {
          c1_c_M[c1_i17] = c1_M[c1_i17];
        }

        c1_d_mpower(chartInstance, c1_c_M, c1_l - 1.0, c1_b_b);
        c1_c_eml_scalar_eg(chartInstance);
        c1_c_eml_scalar_eg(chartInstance);
        for (c1_i18 = 0; c1_i18 < 9; c1_i18++) {
          c1_y[c1_i18] = 0.0;
        }

        for (c1_i19 = 0; c1_i19 < 9; c1_i19++) {
          c1_c_b[c1_i19] = c1_b[c1_i19];
        }

        for (c1_i20 = 0; c1_i20 < 9; c1_i20++) {
          c1_d_b[c1_i20] = c1_b_b[c1_i20];
        }

        c1_b_eml_xgemm(chartInstance, c1_c_b, c1_d_b, c1_y);
        for (c1_i21 = 0; c1_i21 < 9; c1_i21++) {
          c1_Phi_r12[c1_i21] += c1_y[c1_i21];
        }

        c1_d_k++;
        _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
      }

      CV_EML_FOR(0, 1, 2, 0);
      c1_b_j++;
      _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
    }

    CV_EML_FOR(0, 1, 1, 0);
    c1_b_l++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 17);
  for (c1_i22 = 0; c1_i22 < 9; c1_i22++) {
    c1_b[c1_i22] = c1_Phi_r12[c1_i22];
  }

  for (c1_i23 = 0; c1_i23 < 9; c1_i23++) {
    c1_Phi_r12[c1_i23] = 0.5 * c1_b[c1_i23];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -17);
  _SFD_SYMBOL_SCOPE_POP();
  for (c1_i24 = 0; c1_i24 < 9; c1_i24++) {
    (*c1_b_Phi_r12)[c1_i24] = c1_Phi_r12[c1_i24];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
}

static void initSimStructsc1_Model_01(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber)
{
  (void)c1_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i25;
  int32_T c1_i26;
  int32_T c1_i27;
  real_T c1_b_inData[9];
  int32_T c1_i28;
  int32_T c1_i29;
  int32_T c1_i30;
  real_T c1_u[9];
  const mxArray *c1_y = NULL;
  SFc1_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc1_Model_01InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i25 = 0;
  for (c1_i26 = 0; c1_i26 < 3; c1_i26++) {
    for (c1_i27 = 0; c1_i27 < 3; c1_i27++) {
      c1_b_inData[c1_i27 + c1_i25] = (*(real_T (*)[9])c1_inData)[c1_i27 + c1_i25];
    }

    c1_i25 += 3;
  }

  c1_i28 = 0;
  for (c1_i29 = 0; c1_i29 < 3; c1_i29++) {
    for (c1_i30 = 0; c1_i30 < 3; c1_i30++) {
      c1_u[c1_i30 + c1_i28] = c1_b_inData[c1_i30 + c1_i28];
    }

    c1_i28 += 3;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_Phi_r12, const char_T *c1_identifier, real_T c1_y[9])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Phi_r12), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_Phi_r12);
}

static void c1_b_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[9])
{
  real_T c1_dv2[9];
  int32_T c1_i31;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv2, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c1_i31 = 0; c1_i31 < 9; c1_i31++) {
    c1_y[c1_i31] = c1_dv2[c1_i31];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Phi_r12;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[9];
  int32_T c1_i32;
  int32_T c1_i33;
  int32_T c1_i34;
  SFc1_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc1_Model_01InstanceStruct *)chartInstanceVoid;
  c1_Phi_r12 = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Phi_r12), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_Phi_r12);
  c1_i32 = 0;
  for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
    for (c1_i34 = 0; c1_i34 < 3; c1_i34++) {
      (*(real_T (*)[9])c1_outData)[c1_i34 + c1_i32] = c1_y[c1_i34 + c1_i32];
    }

    c1_i32 += 3;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  real_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc1_Model_01InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(real_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i35;
  real_T c1_b_inData[3];
  int32_T c1_i36;
  real_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc1_Model_01InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i35 = 0; c1_i35 < 3; c1_i35++) {
    c1_b_inData[c1_i35] = (*(real_T (*)[3])c1_inData)[c1_i35];
  }

  for (c1_i36 = 0; c1_i36 < 3; c1_i36++) {
    c1_u[c1_i36] = c1_b_inData[c1_i36];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static real_T c1_c_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d3;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d3, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d3;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_nargout;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y;
  SFc1_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc1_Model_01InstanceStruct *)chartInstanceVoid;
  c1_nargout = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_nargout), &c1_thisId);
  sf_mex_destroy(&c1_nargout);
  *(real_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static void c1_d_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3])
{
  real_T c1_dv3[3];
  int32_T c1_i37;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv3, 1, 0, 0U, 1, 0U, 1, 3);
  for (c1_i37 = 0; c1_i37 < 3; c1_i37++) {
    c1_y[c1_i37] = c1_dv3[c1_i37];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_v_lambda;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[3];
  int32_T c1_i38;
  SFc1_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc1_Model_01InstanceStruct *)chartInstanceVoid;
  c1_v_lambda = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_v_lambda), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_v_lambda);
  for (c1_i38 = 0; c1_i38 < 3; c1_i38++) {
    (*(real_T (*)[3])c1_outData)[c1_i38] = c1_y[c1_i38];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_Model_01_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  sf_mex_assign(&c1_nameCaptureInfo, sf_mex_createstruct("structure", 2, 280, 1),
                false);
  c1_info_helper(&c1_nameCaptureInfo);
  c1_b_info_helper(&c1_nameCaptureInfo);
  c1_c_info_helper(&c1_nameCaptureInfo);
  c1_d_info_helper(&c1_nameCaptureInfo);
  c1_e_info_helper(&c1_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c1_nameCaptureInfo);
  return c1_nameCaptureInfo;
}

static void c1_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs0 = NULL;
  const mxArray *c1_lhs0 = NULL;
  const mxArray *c1_rhs1 = NULL;
  const mxArray *c1_lhs1 = NULL;
  const mxArray *c1_rhs2 = NULL;
  const mxArray *c1_lhs2 = NULL;
  const mxArray *c1_rhs3 = NULL;
  const mxArray *c1_lhs3 = NULL;
  const mxArray *c1_rhs4 = NULL;
  const mxArray *c1_lhs4 = NULL;
  const mxArray *c1_rhs5 = NULL;
  const mxArray *c1_lhs5 = NULL;
  const mxArray *c1_rhs6 = NULL;
  const mxArray *c1_lhs6 = NULL;
  const mxArray *c1_rhs7 = NULL;
  const mxArray *c1_lhs7 = NULL;
  const mxArray *c1_rhs8 = NULL;
  const mxArray *c1_lhs8 = NULL;
  const mxArray *c1_rhs9 = NULL;
  const mxArray *c1_lhs9 = NULL;
  const mxArray *c1_rhs10 = NULL;
  const mxArray *c1_lhs10 = NULL;
  const mxArray *c1_rhs11 = NULL;
  const mxArray *c1_lhs11 = NULL;
  const mxArray *c1_rhs12 = NULL;
  const mxArray *c1_lhs12 = NULL;
  const mxArray *c1_rhs13 = NULL;
  const mxArray *c1_lhs13 = NULL;
  const mxArray *c1_rhs14 = NULL;
  const mxArray *c1_lhs14 = NULL;
  const mxArray *c1_rhs15 = NULL;
  const mxArray *c1_lhs15 = NULL;
  const mxArray *c1_rhs16 = NULL;
  const mxArray *c1_lhs16 = NULL;
  const mxArray *c1_rhs17 = NULL;
  const mxArray *c1_lhs17 = NULL;
  const mxArray *c1_rhs18 = NULL;
  const mxArray *c1_lhs18 = NULL;
  const mxArray *c1_rhs19 = NULL;
  const mxArray *c1_lhs19 = NULL;
  const mxArray *c1_rhs20 = NULL;
  const mxArray *c1_lhs20 = NULL;
  const mxArray *c1_rhs21 = NULL;
  const mxArray *c1_lhs21 = NULL;
  const mxArray *c1_rhs22 = NULL;
  const mxArray *c1_lhs22 = NULL;
  const mxArray *c1_rhs23 = NULL;
  const mxArray *c1_lhs23 = NULL;
  const mxArray *c1_rhs24 = NULL;
  const mxArray *c1_lhs24 = NULL;
  const mxArray *c1_rhs25 = NULL;
  const mxArray *c1_lhs25 = NULL;
  const mxArray *c1_rhs26 = NULL;
  const mxArray *c1_lhs26 = NULL;
  const mxArray *c1_rhs27 = NULL;
  const mxArray *c1_lhs27 = NULL;
  const mxArray *c1_rhs28 = NULL;
  const mxArray *c1_lhs28 = NULL;
  const mxArray *c1_rhs29 = NULL;
  const mxArray *c1_lhs29 = NULL;
  const mxArray *c1_rhs30 = NULL;
  const mxArray *c1_lhs30 = NULL;
  const mxArray *c1_rhs31 = NULL;
  const mxArray *c1_lhs31 = NULL;
  const mxArray *c1_rhs32 = NULL;
  const mxArray *c1_lhs32 = NULL;
  const mxArray *c1_rhs33 = NULL;
  const mxArray *c1_lhs33 = NULL;
  const mxArray *c1_rhs34 = NULL;
  const mxArray *c1_lhs34 = NULL;
  const mxArray *c1_rhs35 = NULL;
  const mxArray *c1_lhs35 = NULL;
  const mxArray *c1_rhs36 = NULL;
  const mxArray *c1_lhs36 = NULL;
  const mxArray *c1_rhs37 = NULL;
  const mxArray *c1_lhs37 = NULL;
  const mxArray *c1_rhs38 = NULL;
  const mxArray *c1_lhs38 = NULL;
  const mxArray *c1_rhs39 = NULL;
  const mxArray *c1_lhs39 = NULL;
  const mxArray *c1_rhs40 = NULL;
  const mxArray *c1_lhs40 = NULL;
  const mxArray *c1_rhs41 = NULL;
  const mxArray *c1_lhs41 = NULL;
  const mxArray *c1_rhs42 = NULL;
  const mxArray *c1_lhs42 = NULL;
  const mxArray *c1_rhs43 = NULL;
  const mxArray *c1_lhs43 = NULL;
  const mxArray *c1_rhs44 = NULL;
  const mxArray *c1_lhs44 = NULL;
  const mxArray *c1_rhs45 = NULL;
  const mxArray *c1_lhs45 = NULL;
  const mxArray *c1_rhs46 = NULL;
  const mxArray *c1_lhs46 = NULL;
  const mxArray *c1_rhs47 = NULL;
  const mxArray *c1_lhs47 = NULL;
  const mxArray *c1_rhs48 = NULL;
  const mxArray *c1_lhs48 = NULL;
  const mxArray *c1_rhs49 = NULL;
  const mxArray *c1_lhs49 = NULL;
  const mxArray *c1_rhs50 = NULL;
  const mxArray *c1_lhs50 = NULL;
  const mxArray *c1_rhs51 = NULL;
  const mxArray *c1_lhs51 = NULL;
  const mxArray *c1_rhs52 = NULL;
  const mxArray *c1_lhs52 = NULL;
  const mxArray *c1_rhs53 = NULL;
  const mxArray *c1_lhs53 = NULL;
  const mxArray *c1_rhs54 = NULL;
  const mxArray *c1_lhs54 = NULL;
  const mxArray *c1_rhs55 = NULL;
  const mxArray *c1_lhs55 = NULL;
  const mxArray *c1_rhs56 = NULL;
  const mxArray *c1_lhs56 = NULL;
  const mxArray *c1_rhs57 = NULL;
  const mxArray *c1_lhs57 = NULL;
  const mxArray *c1_rhs58 = NULL;
  const mxArray *c1_lhs58 = NULL;
  const mxArray *c1_rhs59 = NULL;
  const mxArray *c1_lhs59 = NULL;
  const mxArray *c1_rhs60 = NULL;
  const mxArray *c1_lhs60 = NULL;
  const mxArray *c1_rhs61 = NULL;
  const mxArray *c1_lhs61 = NULL;
  const mxArray *c1_rhs62 = NULL;
  const mxArray *c1_lhs62 = NULL;
  const mxArray *c1_rhs63 = NULL;
  const mxArray *c1_lhs63 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eig"), "name", "name", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c1_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c1_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c1_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgeev"), "name", "name", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826004U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c1_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m"),
                  "context", "context", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_lapack_xgeev"), "name",
                  "name", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1301335668U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c1_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zggev"), "name",
                  "name", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826018U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c1_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmin"), "name", "name", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c1_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_realmin"), "name", "name",
                  7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658444U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c1_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c1_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sqrt"), "name", "name", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c1_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_error"), "name", "name",
                  10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c1_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825938U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c1_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c1_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c1_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_eps"), "name", "name", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c1_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c1_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zlangeM"), "name",
                  "name", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826020U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c1_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c1_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c1_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c1_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "context", "context", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1350417854U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c1_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "context", "context", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c1_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c1_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "context", "context", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c1_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "context", "context", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c1_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isfinite"), "name", "name", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c1_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c1_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isinf"), "name", "name", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c1_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c1_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c1_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c1_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zlascl"), "name",
                  "name", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826022U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c1_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmin"), "name", "name", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c1_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c1_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c1_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c1_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c1_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zggbal"), "name",
                  "name", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826018U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c1_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "context", "context", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c1_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_rows"),
                  "context", "context", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c1_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_rows"),
                  "context", "context", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c1_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c1_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c1_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_rows"),
                  "context", "context", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c1_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c1_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_simtran"),
                  "context", "context", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c1_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_simtran"),
                  "context", "context", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c1_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "context", "context", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c1_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_cols"),
                  "context", "context", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c1_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_cols"),
                  "context", "context", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c1_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "context", "context", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c1_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c1_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zgghrd"), "name",
                  "name", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826020U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c1_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c1_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c1_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c1_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c1_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c1_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zlartg"), "name",
                  "name", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826022U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c1_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmin"), "name", "name", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c1_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c1_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("fix"), "name", "name", 61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m"), "resolved",
                  "resolved", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c1_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m"), "context",
                  "context", 62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c1_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m"), "context",
                  "context", 63);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_fix"), "name",
                  "name", 63);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_fix.m"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658438U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c1_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c1_rhs0);
  sf_mex_destroy(&c1_lhs0);
  sf_mex_destroy(&c1_rhs1);
  sf_mex_destroy(&c1_lhs1);
  sf_mex_destroy(&c1_rhs2);
  sf_mex_destroy(&c1_lhs2);
  sf_mex_destroy(&c1_rhs3);
  sf_mex_destroy(&c1_lhs3);
  sf_mex_destroy(&c1_rhs4);
  sf_mex_destroy(&c1_lhs4);
  sf_mex_destroy(&c1_rhs5);
  sf_mex_destroy(&c1_lhs5);
  sf_mex_destroy(&c1_rhs6);
  sf_mex_destroy(&c1_lhs6);
  sf_mex_destroy(&c1_rhs7);
  sf_mex_destroy(&c1_lhs7);
  sf_mex_destroy(&c1_rhs8);
  sf_mex_destroy(&c1_lhs8);
  sf_mex_destroy(&c1_rhs9);
  sf_mex_destroy(&c1_lhs9);
  sf_mex_destroy(&c1_rhs10);
  sf_mex_destroy(&c1_lhs10);
  sf_mex_destroy(&c1_rhs11);
  sf_mex_destroy(&c1_lhs11);
  sf_mex_destroy(&c1_rhs12);
  sf_mex_destroy(&c1_lhs12);
  sf_mex_destroy(&c1_rhs13);
  sf_mex_destroy(&c1_lhs13);
  sf_mex_destroy(&c1_rhs14);
  sf_mex_destroy(&c1_lhs14);
  sf_mex_destroy(&c1_rhs15);
  sf_mex_destroy(&c1_lhs15);
  sf_mex_destroy(&c1_rhs16);
  sf_mex_destroy(&c1_lhs16);
  sf_mex_destroy(&c1_rhs17);
  sf_mex_destroy(&c1_lhs17);
  sf_mex_destroy(&c1_rhs18);
  sf_mex_destroy(&c1_lhs18);
  sf_mex_destroy(&c1_rhs19);
  sf_mex_destroy(&c1_lhs19);
  sf_mex_destroy(&c1_rhs20);
  sf_mex_destroy(&c1_lhs20);
  sf_mex_destroy(&c1_rhs21);
  sf_mex_destroy(&c1_lhs21);
  sf_mex_destroy(&c1_rhs22);
  sf_mex_destroy(&c1_lhs22);
  sf_mex_destroy(&c1_rhs23);
  sf_mex_destroy(&c1_lhs23);
  sf_mex_destroy(&c1_rhs24);
  sf_mex_destroy(&c1_lhs24);
  sf_mex_destroy(&c1_rhs25);
  sf_mex_destroy(&c1_lhs25);
  sf_mex_destroy(&c1_rhs26);
  sf_mex_destroy(&c1_lhs26);
  sf_mex_destroy(&c1_rhs27);
  sf_mex_destroy(&c1_lhs27);
  sf_mex_destroy(&c1_rhs28);
  sf_mex_destroy(&c1_lhs28);
  sf_mex_destroy(&c1_rhs29);
  sf_mex_destroy(&c1_lhs29);
  sf_mex_destroy(&c1_rhs30);
  sf_mex_destroy(&c1_lhs30);
  sf_mex_destroy(&c1_rhs31);
  sf_mex_destroy(&c1_lhs31);
  sf_mex_destroy(&c1_rhs32);
  sf_mex_destroy(&c1_lhs32);
  sf_mex_destroy(&c1_rhs33);
  sf_mex_destroy(&c1_lhs33);
  sf_mex_destroy(&c1_rhs34);
  sf_mex_destroy(&c1_lhs34);
  sf_mex_destroy(&c1_rhs35);
  sf_mex_destroy(&c1_lhs35);
  sf_mex_destroy(&c1_rhs36);
  sf_mex_destroy(&c1_lhs36);
  sf_mex_destroy(&c1_rhs37);
  sf_mex_destroy(&c1_lhs37);
  sf_mex_destroy(&c1_rhs38);
  sf_mex_destroy(&c1_lhs38);
  sf_mex_destroy(&c1_rhs39);
  sf_mex_destroy(&c1_lhs39);
  sf_mex_destroy(&c1_rhs40);
  sf_mex_destroy(&c1_lhs40);
  sf_mex_destroy(&c1_rhs41);
  sf_mex_destroy(&c1_lhs41);
  sf_mex_destroy(&c1_rhs42);
  sf_mex_destroy(&c1_lhs42);
  sf_mex_destroy(&c1_rhs43);
  sf_mex_destroy(&c1_lhs43);
  sf_mex_destroy(&c1_rhs44);
  sf_mex_destroy(&c1_lhs44);
  sf_mex_destroy(&c1_rhs45);
  sf_mex_destroy(&c1_lhs45);
  sf_mex_destroy(&c1_rhs46);
  sf_mex_destroy(&c1_lhs46);
  sf_mex_destroy(&c1_rhs47);
  sf_mex_destroy(&c1_lhs47);
  sf_mex_destroy(&c1_rhs48);
  sf_mex_destroy(&c1_lhs48);
  sf_mex_destroy(&c1_rhs49);
  sf_mex_destroy(&c1_lhs49);
  sf_mex_destroy(&c1_rhs50);
  sf_mex_destroy(&c1_lhs50);
  sf_mex_destroy(&c1_rhs51);
  sf_mex_destroy(&c1_lhs51);
  sf_mex_destroy(&c1_rhs52);
  sf_mex_destroy(&c1_lhs52);
  sf_mex_destroy(&c1_rhs53);
  sf_mex_destroy(&c1_lhs53);
  sf_mex_destroy(&c1_rhs54);
  sf_mex_destroy(&c1_lhs54);
  sf_mex_destroy(&c1_rhs55);
  sf_mex_destroy(&c1_lhs55);
  sf_mex_destroy(&c1_rhs56);
  sf_mex_destroy(&c1_lhs56);
  sf_mex_destroy(&c1_rhs57);
  sf_mex_destroy(&c1_lhs57);
  sf_mex_destroy(&c1_rhs58);
  sf_mex_destroy(&c1_lhs58);
  sf_mex_destroy(&c1_rhs59);
  sf_mex_destroy(&c1_lhs59);
  sf_mex_destroy(&c1_rhs60);
  sf_mex_destroy(&c1_lhs60);
  sf_mex_destroy(&c1_rhs61);
  sf_mex_destroy(&c1_lhs61);
  sf_mex_destroy(&c1_rhs62);
  sf_mex_destroy(&c1_lhs62);
  sf_mex_destroy(&c1_rhs63);
  sf_mex_destroy(&c1_lhs63);
}

static const mxArray *c1_emlrt_marshallOut(const char * c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c1_u)), false);
  return c1_y;
}

static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 7, 0U, 0U, 0U, 0), false);
  return c1_y;
}

static void c1_b_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs64 = NULL;
  const mxArray *c1_lhs64 = NULL;
  const mxArray *c1_rhs65 = NULL;
  const mxArray *c1_lhs65 = NULL;
  const mxArray *c1_rhs66 = NULL;
  const mxArray *c1_lhs66 = NULL;
  const mxArray *c1_rhs67 = NULL;
  const mxArray *c1_lhs67 = NULL;
  const mxArray *c1_rhs68 = NULL;
  const mxArray *c1_lhs68 = NULL;
  const mxArray *c1_rhs69 = NULL;
  const mxArray *c1_lhs69 = NULL;
  const mxArray *c1_rhs70 = NULL;
  const mxArray *c1_lhs70 = NULL;
  const mxArray *c1_rhs71 = NULL;
  const mxArray *c1_lhs71 = NULL;
  const mxArray *c1_rhs72 = NULL;
  const mxArray *c1_lhs72 = NULL;
  const mxArray *c1_rhs73 = NULL;
  const mxArray *c1_lhs73 = NULL;
  const mxArray *c1_rhs74 = NULL;
  const mxArray *c1_lhs74 = NULL;
  const mxArray *c1_rhs75 = NULL;
  const mxArray *c1_lhs75 = NULL;
  const mxArray *c1_rhs76 = NULL;
  const mxArray *c1_lhs76 = NULL;
  const mxArray *c1_rhs77 = NULL;
  const mxArray *c1_lhs77 = NULL;
  const mxArray *c1_rhs78 = NULL;
  const mxArray *c1_lhs78 = NULL;
  const mxArray *c1_rhs79 = NULL;
  const mxArray *c1_lhs79 = NULL;
  const mxArray *c1_rhs80 = NULL;
  const mxArray *c1_lhs80 = NULL;
  const mxArray *c1_rhs81 = NULL;
  const mxArray *c1_lhs81 = NULL;
  const mxArray *c1_rhs82 = NULL;
  const mxArray *c1_lhs82 = NULL;
  const mxArray *c1_rhs83 = NULL;
  const mxArray *c1_lhs83 = NULL;
  const mxArray *c1_rhs84 = NULL;
  const mxArray *c1_lhs84 = NULL;
  const mxArray *c1_rhs85 = NULL;
  const mxArray *c1_lhs85 = NULL;
  const mxArray *c1_rhs86 = NULL;
  const mxArray *c1_lhs86 = NULL;
  const mxArray *c1_rhs87 = NULL;
  const mxArray *c1_lhs87 = NULL;
  const mxArray *c1_rhs88 = NULL;
  const mxArray *c1_lhs88 = NULL;
  const mxArray *c1_rhs89 = NULL;
  const mxArray *c1_lhs89 = NULL;
  const mxArray *c1_rhs90 = NULL;
  const mxArray *c1_lhs90 = NULL;
  const mxArray *c1_rhs91 = NULL;
  const mxArray *c1_lhs91 = NULL;
  const mxArray *c1_rhs92 = NULL;
  const mxArray *c1_lhs92 = NULL;
  const mxArray *c1_rhs93 = NULL;
  const mxArray *c1_lhs93 = NULL;
  const mxArray *c1_rhs94 = NULL;
  const mxArray *c1_lhs94 = NULL;
  const mxArray *c1_rhs95 = NULL;
  const mxArray *c1_lhs95 = NULL;
  const mxArray *c1_rhs96 = NULL;
  const mxArray *c1_lhs96 = NULL;
  const mxArray *c1_rhs97 = NULL;
  const mxArray *c1_lhs97 = NULL;
  const mxArray *c1_rhs98 = NULL;
  const mxArray *c1_lhs98 = NULL;
  const mxArray *c1_rhs99 = NULL;
  const mxArray *c1_lhs99 = NULL;
  const mxArray *c1_rhs100 = NULL;
  const mxArray *c1_lhs100 = NULL;
  const mxArray *c1_rhs101 = NULL;
  const mxArray *c1_lhs101 = NULL;
  const mxArray *c1_rhs102 = NULL;
  const mxArray *c1_lhs102 = NULL;
  const mxArray *c1_rhs103 = NULL;
  const mxArray *c1_lhs103 = NULL;
  const mxArray *c1_rhs104 = NULL;
  const mxArray *c1_lhs104 = NULL;
  const mxArray *c1_rhs105 = NULL;
  const mxArray *c1_lhs105 = NULL;
  const mxArray *c1_rhs106 = NULL;
  const mxArray *c1_lhs106 = NULL;
  const mxArray *c1_rhs107 = NULL;
  const mxArray *c1_lhs107 = NULL;
  const mxArray *c1_rhs108 = NULL;
  const mxArray *c1_lhs108 = NULL;
  const mxArray *c1_rhs109 = NULL;
  const mxArray *c1_lhs109 = NULL;
  const mxArray *c1_rhs110 = NULL;
  const mxArray *c1_lhs110 = NULL;
  const mxArray *c1_rhs111 = NULL;
  const mxArray *c1_lhs111 = NULL;
  const mxArray *c1_rhs112 = NULL;
  const mxArray *c1_lhs112 = NULL;
  const mxArray *c1_rhs113 = NULL;
  const mxArray *c1_lhs113 = NULL;
  const mxArray *c1_rhs114 = NULL;
  const mxArray *c1_lhs114 = NULL;
  const mxArray *c1_rhs115 = NULL;
  const mxArray *c1_lhs115 = NULL;
  const mxArray *c1_rhs116 = NULL;
  const mxArray *c1_lhs116 = NULL;
  const mxArray *c1_rhs117 = NULL;
  const mxArray *c1_lhs117 = NULL;
  const mxArray *c1_rhs118 = NULL;
  const mxArray *c1_lhs118 = NULL;
  const mxArray *c1_rhs119 = NULL;
  const mxArray *c1_lhs119 = NULL;
  const mxArray *c1_rhs120 = NULL;
  const mxArray *c1_lhs120 = NULL;
  const mxArray *c1_rhs121 = NULL;
  const mxArray *c1_lhs121 = NULL;
  const mxArray *c1_rhs122 = NULL;
  const mxArray *c1_lhs122 = NULL;
  const mxArray *c1_rhs123 = NULL;
  const mxArray *c1_lhs123 = NULL;
  const mxArray *c1_rhs124 = NULL;
  const mxArray *c1_lhs124 = NULL;
  const mxArray *c1_rhs125 = NULL;
  const mxArray *c1_lhs125 = NULL;
  const mxArray *c1_rhs126 = NULL;
  const mxArray *c1_lhs126 = NULL;
  const mxArray *c1_rhs127 = NULL;
  const mxArray *c1_lhs127 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mpower"), "name", "name", 64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c1_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c1_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("ismatrix"), "name", "name", 66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c1_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("power"), "name", "name", 67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c1_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c1_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c1_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c1_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c1_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("floor"), "name", "name", 72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c1_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c1_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c1_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c1_rhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c1_rhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m!absinf"),
                  "context", "context", 77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c1_rhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c1_rhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c1_rhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1350417854U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c1_rhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sqrt"), "name", "name", 81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c1_rhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c1_rhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_zrot_rows"), "name",
                  "name", 83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "resolved", "resolved", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c1_rhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "context", "context", 84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c1_rhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "context", "context", 85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c1_rhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "context", "context", 86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.conjtimes"),
                  "name", "name", 86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/conjtimes.m"),
                  "resolved", "resolved", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1360286186U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c1_rhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_zrot_cols"), "name",
                  "name", 87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c1_rhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "context", "context", 88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c1_rhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "context", "context", 89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c1_rhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "context", "context", 90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.conjtimes"),
                  "name", "name", 90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/conjtimes.m"),
                  "resolved", "resolved", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1360286186U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c1_rhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zhgeqz"), "name",
                  "name", 91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1368190232U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c1_rhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c1_rhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c1_rhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmin"), "name", "name", 94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c1_rhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zlanhs"), "name",
                  "name", 95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "resolved", "resolved", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826020U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c1_rhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c1_rhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c1_rhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c1_rhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c1_rhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sqrt"), "name", "name", 100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c1_rhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs100), "rhs", "rhs",
                  100);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs100), "lhs", "lhs",
                  100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c1_rhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs101), "rhs", "rhs",
                  101);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs101), "lhs", "lhs",
                  101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c1_rhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs102), "rhs", "rhs",
                  102);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs102), "lhs", "lhs",
                  102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c1_rhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs103), "rhs", "rhs",
                  103);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs103), "lhs", "lhs",
                  103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c1_rhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs104), "rhs", "rhs",
                  104);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs104), "lhs", "lhs",
                  104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c1_rhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs105), "rhs", "rhs",
                  105);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs105), "lhs", "lhs",
                  105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m!abs1"),
                  "context", "context", 106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c1_rhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs106), "rhs", "rhs",
                  106);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs106), "lhs", "lhs",
                  106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c1_rhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs107), "rhs", "rhs",
                  107);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs107), "lhs", "lhs",
                  107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zlartg"), "name",
                  "name", 108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826022U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c1_rhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs108), "rhs", "rhs",
                  108);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs108), "lhs", "lhs",
                  108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_zrot_cols"), "name",
                  "name", 109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "resolved", "resolved", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c1_rhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs109), "rhs", "rhs",
                  109);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs109), "lhs", "lhs",
                  109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mod"), "name", "name", 110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "resolved",
                  "resolved", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c1_rhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs110), "rhs", "rhs",
                  110);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs110), "lhs", "lhs",
                  110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c1_rhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs111), "rhs", "rhs",
                  111);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs111), "lhs", "lhs",
                  111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c1_rhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs112), "rhs", "rhs",
                  112);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs112), "lhs", "lhs",
                  112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c1_rhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs113), "rhs", "rhs",
                  113);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs113), "lhs", "lhs",
                  113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c1_rhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs114), "rhs", "rhs",
                  114);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs114), "lhs", "lhs",
                  114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c1_rhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs115), "rhs", "rhs",
                  115);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs115), "lhs", "lhs",
                  115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!intmod"), "context",
                  "context", 116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c1_rhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs116), "rhs", "rhs",
                  116);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs116), "lhs", "lhs",
                  116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c1_rhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs117), "rhs", "rhs",
                  117);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs117), "lhs", "lhs",
                  117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c1_rhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs118), "rhs", "rhs",
                  118);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs118), "lhs", "lhs",
                  118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c1_rhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs119), "rhs", "rhs",
                  119);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs119), "lhs", "lhs",
                  119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p!eml_fldiv"),
                  "context", "context", 120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c1_rhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs120), "rhs", "rhs",
                  120);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs120), "lhs", "lhs",
                  120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p!eml_fldiv"),
                  "context", "context", 121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c1_rhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs121), "rhs", "rhs",
                  121);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs121), "lhs", "lhs",
                  121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p!eml_fldiv"),
                  "context", "context", 122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c1_rhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs122), "rhs", "rhs",
                  122);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs122), "lhs", "lhs",
                  122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sqrt"), "name", "name", 123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c1_rhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs123), "rhs", "rhs",
                  123);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs123), "lhs", "lhs",
                  123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c1_rhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs124), "rhs", "rhs",
                  124);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs124), "lhs", "lhs",
                  124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c1_rhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs125), "rhs", "rhs",
                  125);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs125), "lhs", "lhs",
                  125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c1_rhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs126), "rhs", "rhs",
                  126);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs126), "lhs", "lhs",
                  126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 127);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 127);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c1_rhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs127), "rhs", "rhs",
                  127);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs127), "lhs", "lhs",
                  127);
  sf_mex_destroy(&c1_rhs64);
  sf_mex_destroy(&c1_lhs64);
  sf_mex_destroy(&c1_rhs65);
  sf_mex_destroy(&c1_lhs65);
  sf_mex_destroy(&c1_rhs66);
  sf_mex_destroy(&c1_lhs66);
  sf_mex_destroy(&c1_rhs67);
  sf_mex_destroy(&c1_lhs67);
  sf_mex_destroy(&c1_rhs68);
  sf_mex_destroy(&c1_lhs68);
  sf_mex_destroy(&c1_rhs69);
  sf_mex_destroy(&c1_lhs69);
  sf_mex_destroy(&c1_rhs70);
  sf_mex_destroy(&c1_lhs70);
  sf_mex_destroy(&c1_rhs71);
  sf_mex_destroy(&c1_lhs71);
  sf_mex_destroy(&c1_rhs72);
  sf_mex_destroy(&c1_lhs72);
  sf_mex_destroy(&c1_rhs73);
  sf_mex_destroy(&c1_lhs73);
  sf_mex_destroy(&c1_rhs74);
  sf_mex_destroy(&c1_lhs74);
  sf_mex_destroy(&c1_rhs75);
  sf_mex_destroy(&c1_lhs75);
  sf_mex_destroy(&c1_rhs76);
  sf_mex_destroy(&c1_lhs76);
  sf_mex_destroy(&c1_rhs77);
  sf_mex_destroy(&c1_lhs77);
  sf_mex_destroy(&c1_rhs78);
  sf_mex_destroy(&c1_lhs78);
  sf_mex_destroy(&c1_rhs79);
  sf_mex_destroy(&c1_lhs79);
  sf_mex_destroy(&c1_rhs80);
  sf_mex_destroy(&c1_lhs80);
  sf_mex_destroy(&c1_rhs81);
  sf_mex_destroy(&c1_lhs81);
  sf_mex_destroy(&c1_rhs82);
  sf_mex_destroy(&c1_lhs82);
  sf_mex_destroy(&c1_rhs83);
  sf_mex_destroy(&c1_lhs83);
  sf_mex_destroy(&c1_rhs84);
  sf_mex_destroy(&c1_lhs84);
  sf_mex_destroy(&c1_rhs85);
  sf_mex_destroy(&c1_lhs85);
  sf_mex_destroy(&c1_rhs86);
  sf_mex_destroy(&c1_lhs86);
  sf_mex_destroy(&c1_rhs87);
  sf_mex_destroy(&c1_lhs87);
  sf_mex_destroy(&c1_rhs88);
  sf_mex_destroy(&c1_lhs88);
  sf_mex_destroy(&c1_rhs89);
  sf_mex_destroy(&c1_lhs89);
  sf_mex_destroy(&c1_rhs90);
  sf_mex_destroy(&c1_lhs90);
  sf_mex_destroy(&c1_rhs91);
  sf_mex_destroy(&c1_lhs91);
  sf_mex_destroy(&c1_rhs92);
  sf_mex_destroy(&c1_lhs92);
  sf_mex_destroy(&c1_rhs93);
  sf_mex_destroy(&c1_lhs93);
  sf_mex_destroy(&c1_rhs94);
  sf_mex_destroy(&c1_lhs94);
  sf_mex_destroy(&c1_rhs95);
  sf_mex_destroy(&c1_lhs95);
  sf_mex_destroy(&c1_rhs96);
  sf_mex_destroy(&c1_lhs96);
  sf_mex_destroy(&c1_rhs97);
  sf_mex_destroy(&c1_lhs97);
  sf_mex_destroy(&c1_rhs98);
  sf_mex_destroy(&c1_lhs98);
  sf_mex_destroy(&c1_rhs99);
  sf_mex_destroy(&c1_lhs99);
  sf_mex_destroy(&c1_rhs100);
  sf_mex_destroy(&c1_lhs100);
  sf_mex_destroy(&c1_rhs101);
  sf_mex_destroy(&c1_lhs101);
  sf_mex_destroy(&c1_rhs102);
  sf_mex_destroy(&c1_lhs102);
  sf_mex_destroy(&c1_rhs103);
  sf_mex_destroy(&c1_lhs103);
  sf_mex_destroy(&c1_rhs104);
  sf_mex_destroy(&c1_lhs104);
  sf_mex_destroy(&c1_rhs105);
  sf_mex_destroy(&c1_lhs105);
  sf_mex_destroy(&c1_rhs106);
  sf_mex_destroy(&c1_lhs106);
  sf_mex_destroy(&c1_rhs107);
  sf_mex_destroy(&c1_lhs107);
  sf_mex_destroy(&c1_rhs108);
  sf_mex_destroy(&c1_lhs108);
  sf_mex_destroy(&c1_rhs109);
  sf_mex_destroy(&c1_lhs109);
  sf_mex_destroy(&c1_rhs110);
  sf_mex_destroy(&c1_lhs110);
  sf_mex_destroy(&c1_rhs111);
  sf_mex_destroy(&c1_lhs111);
  sf_mex_destroy(&c1_rhs112);
  sf_mex_destroy(&c1_lhs112);
  sf_mex_destroy(&c1_rhs113);
  sf_mex_destroy(&c1_lhs113);
  sf_mex_destroy(&c1_rhs114);
  sf_mex_destroy(&c1_lhs114);
  sf_mex_destroy(&c1_rhs115);
  sf_mex_destroy(&c1_lhs115);
  sf_mex_destroy(&c1_rhs116);
  sf_mex_destroy(&c1_lhs116);
  sf_mex_destroy(&c1_rhs117);
  sf_mex_destroy(&c1_lhs117);
  sf_mex_destroy(&c1_rhs118);
  sf_mex_destroy(&c1_lhs118);
  sf_mex_destroy(&c1_rhs119);
  sf_mex_destroy(&c1_lhs119);
  sf_mex_destroy(&c1_rhs120);
  sf_mex_destroy(&c1_lhs120);
  sf_mex_destroy(&c1_rhs121);
  sf_mex_destroy(&c1_lhs121);
  sf_mex_destroy(&c1_rhs122);
  sf_mex_destroy(&c1_lhs122);
  sf_mex_destroy(&c1_rhs123);
  sf_mex_destroy(&c1_lhs123);
  sf_mex_destroy(&c1_rhs124);
  sf_mex_destroy(&c1_lhs124);
  sf_mex_destroy(&c1_rhs125);
  sf_mex_destroy(&c1_lhs125);
  sf_mex_destroy(&c1_rhs126);
  sf_mex_destroy(&c1_lhs126);
  sf_mex_destroy(&c1_rhs127);
  sf_mex_destroy(&c1_lhs127);
}

static void c1_c_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs128 = NULL;
  const mxArray *c1_lhs128 = NULL;
  const mxArray *c1_rhs129 = NULL;
  const mxArray *c1_lhs129 = NULL;
  const mxArray *c1_rhs130 = NULL;
  const mxArray *c1_lhs130 = NULL;
  const mxArray *c1_rhs131 = NULL;
  const mxArray *c1_lhs131 = NULL;
  const mxArray *c1_rhs132 = NULL;
  const mxArray *c1_lhs132 = NULL;
  const mxArray *c1_rhs133 = NULL;
  const mxArray *c1_lhs133 = NULL;
  const mxArray *c1_rhs134 = NULL;
  const mxArray *c1_lhs134 = NULL;
  const mxArray *c1_rhs135 = NULL;
  const mxArray *c1_lhs135 = NULL;
  const mxArray *c1_rhs136 = NULL;
  const mxArray *c1_lhs136 = NULL;
  const mxArray *c1_rhs137 = NULL;
  const mxArray *c1_lhs137 = NULL;
  const mxArray *c1_rhs138 = NULL;
  const mxArray *c1_lhs138 = NULL;
  const mxArray *c1_rhs139 = NULL;
  const mxArray *c1_lhs139 = NULL;
  const mxArray *c1_rhs140 = NULL;
  const mxArray *c1_lhs140 = NULL;
  const mxArray *c1_rhs141 = NULL;
  const mxArray *c1_lhs141 = NULL;
  const mxArray *c1_rhs142 = NULL;
  const mxArray *c1_lhs142 = NULL;
  const mxArray *c1_rhs143 = NULL;
  const mxArray *c1_lhs143 = NULL;
  const mxArray *c1_rhs144 = NULL;
  const mxArray *c1_lhs144 = NULL;
  const mxArray *c1_rhs145 = NULL;
  const mxArray *c1_lhs145 = NULL;
  const mxArray *c1_rhs146 = NULL;
  const mxArray *c1_lhs146 = NULL;
  const mxArray *c1_rhs147 = NULL;
  const mxArray *c1_lhs147 = NULL;
  const mxArray *c1_rhs148 = NULL;
  const mxArray *c1_lhs148 = NULL;
  const mxArray *c1_rhs149 = NULL;
  const mxArray *c1_lhs149 = NULL;
  const mxArray *c1_rhs150 = NULL;
  const mxArray *c1_lhs150 = NULL;
  const mxArray *c1_rhs151 = NULL;
  const mxArray *c1_lhs151 = NULL;
  const mxArray *c1_rhs152 = NULL;
  const mxArray *c1_lhs152 = NULL;
  const mxArray *c1_rhs153 = NULL;
  const mxArray *c1_lhs153 = NULL;
  const mxArray *c1_rhs154 = NULL;
  const mxArray *c1_lhs154 = NULL;
  const mxArray *c1_rhs155 = NULL;
  const mxArray *c1_lhs155 = NULL;
  const mxArray *c1_rhs156 = NULL;
  const mxArray *c1_lhs156 = NULL;
  const mxArray *c1_rhs157 = NULL;
  const mxArray *c1_lhs157 = NULL;
  const mxArray *c1_rhs158 = NULL;
  const mxArray *c1_lhs158 = NULL;
  const mxArray *c1_rhs159 = NULL;
  const mxArray *c1_lhs159 = NULL;
  const mxArray *c1_rhs160 = NULL;
  const mxArray *c1_lhs160 = NULL;
  const mxArray *c1_rhs161 = NULL;
  const mxArray *c1_lhs161 = NULL;
  const mxArray *c1_rhs162 = NULL;
  const mxArray *c1_lhs162 = NULL;
  const mxArray *c1_rhs163 = NULL;
  const mxArray *c1_lhs163 = NULL;
  const mxArray *c1_rhs164 = NULL;
  const mxArray *c1_lhs164 = NULL;
  const mxArray *c1_rhs165 = NULL;
  const mxArray *c1_lhs165 = NULL;
  const mxArray *c1_rhs166 = NULL;
  const mxArray *c1_lhs166 = NULL;
  const mxArray *c1_rhs167 = NULL;
  const mxArray *c1_lhs167 = NULL;
  const mxArray *c1_rhs168 = NULL;
  const mxArray *c1_lhs168 = NULL;
  const mxArray *c1_rhs169 = NULL;
  const mxArray *c1_lhs169 = NULL;
  const mxArray *c1_rhs170 = NULL;
  const mxArray *c1_lhs170 = NULL;
  const mxArray *c1_rhs171 = NULL;
  const mxArray *c1_lhs171 = NULL;
  const mxArray *c1_rhs172 = NULL;
  const mxArray *c1_lhs172 = NULL;
  const mxArray *c1_rhs173 = NULL;
  const mxArray *c1_lhs173 = NULL;
  const mxArray *c1_rhs174 = NULL;
  const mxArray *c1_lhs174 = NULL;
  const mxArray *c1_rhs175 = NULL;
  const mxArray *c1_lhs175 = NULL;
  const mxArray *c1_rhs176 = NULL;
  const mxArray *c1_lhs176 = NULL;
  const mxArray *c1_rhs177 = NULL;
  const mxArray *c1_lhs177 = NULL;
  const mxArray *c1_rhs178 = NULL;
  const mxArray *c1_lhs178 = NULL;
  const mxArray *c1_rhs179 = NULL;
  const mxArray *c1_lhs179 = NULL;
  const mxArray *c1_rhs180 = NULL;
  const mxArray *c1_lhs180 = NULL;
  const mxArray *c1_rhs181 = NULL;
  const mxArray *c1_lhs181 = NULL;
  const mxArray *c1_rhs182 = NULL;
  const mxArray *c1_lhs182 = NULL;
  const mxArray *c1_rhs183 = NULL;
  const mxArray *c1_lhs183 = NULL;
  const mxArray *c1_rhs184 = NULL;
  const mxArray *c1_lhs184 = NULL;
  const mxArray *c1_rhs185 = NULL;
  const mxArray *c1_lhs185 = NULL;
  const mxArray *c1_rhs186 = NULL;
  const mxArray *c1_lhs186 = NULL;
  const mxArray *c1_rhs187 = NULL;
  const mxArray *c1_lhs187 = NULL;
  const mxArray *c1_rhs188 = NULL;
  const mxArray *c1_lhs188 = NULL;
  const mxArray *c1_rhs189 = NULL;
  const mxArray *c1_lhs189 = NULL;
  const mxArray *c1_rhs190 = NULL;
  const mxArray *c1_lhs190 = NULL;
  const mxArray *c1_rhs191 = NULL;
  const mxArray *c1_lhs191 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c1_rhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs128), "rhs", "rhs",
                  128);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs128), "lhs", "lhs",
                  128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c1_rhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs129), "rhs", "rhs",
                  129);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs129), "lhs", "lhs",
                  129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isinf"), "name", "name", 130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c1_rhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs130), "rhs", "rhs",
                  130);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs130), "lhs", "lhs",
                  130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_guarded_inf"), "name",
                  "name", 131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_inf.m"),
                  "resolved", "resolved", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c1_rhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs131), "rhs", "rhs",
                  131);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs131), "lhs", "lhs",
                  131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_inf.m"),
                  "context", "context", 132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c1_rhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs132), "rhs", "rhs",
                  132);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs132), "lhs", "lhs",
                  132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmax"), "name", "name", 133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m"), "resolved",
                  "resolved", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c1_rhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs133), "rhs", "rhs",
                  133);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs133), "lhs", "lhs",
                  133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m"), "context",
                  "context", 134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_realmax"), "name", "name",
                  134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmax.m"), "resolved",
                  "resolved", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c1_rhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs134), "rhs", "rhs",
                  134);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs134), "lhs", "lhs",
                  134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmax.m"), "context",
                  "context", 135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c1_rhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs135), "rhs", "rhs",
                  135);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs135), "lhs", "lhs",
                  135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name",
                  136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c1_rhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs136), "rhs", "rhs",
                  136);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs136), "lhs", "lhs",
                  136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c1_rhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs137), "rhs", "rhs",
                  137);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs137), "lhs", "lhs",
                  137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c1_rhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs138), "rhs", "rhs",
                  138);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs138), "lhs", "lhs",
                  138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_hypot"), "name",
                  "name", 139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "resolved", "resolved", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c1_rhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs139), "rhs", "rhs",
                  139);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs139), "lhs", "lhs",
                  139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "context", "context", 140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c1_rhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs140), "rhs", "rhs",
                  140);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs140), "lhs", "lhs",
                  140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "context", "context", 141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1350417854U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c1_rhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs141), "rhs", "rhs",
                  141);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs141), "lhs", "lhs",
                  141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c1_rhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs142), "rhs", "rhs",
                  142);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs142), "lhs", "lhs",
                  142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_zrot_rows"), "name",
                  "name", 143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "resolved", "resolved", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c1_rhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs143), "rhs", "rhs",
                  143);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs143), "lhs", "lhs",
                  143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c1_rhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs144), "rhs", "rhs",
                  144);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs144), "lhs", "lhs",
                  144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p!equalsize"),
                  "context", "context", 145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("max"), "name", "name", 145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "resolved",
                  "resolved", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1311262516U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c1_rhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs145), "rhs", "rhs",
                  145);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs145), "lhs", "lhs",
                  145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "context",
                  "context", 146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1378303184U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c1_rhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs146), "rhs", "rhs",
                  146);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs146), "lhs", "lhs",
                  146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c1_rhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs147), "rhs", "rhs",
                  147);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs147), "lhs", "lhs",
                  147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c1_rhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs148), "rhs", "rhs",
                  148);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs148), "lhs", "lhs",
                  148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c1_rhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs149), "rhs", "rhs",
                  149);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs149), "lhs", "lhs",
                  149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c1_rhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs150), "rhs", "rhs",
                  150);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs150), "lhs", "lhs",
                  150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c1_rhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs151), "rhs", "rhs",
                  151);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs151), "lhs", "lhs",
                  151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_warning"), "name", "name",
                  152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c1_rhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs152), "rhs", "rhs",
                  152);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs152), "lhs", "lhs",
                  152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mpower"), "name", "name", 153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c1_rhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs153), "rhs", "rhs",
                  153);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs153), "lhs", "lhs",
                  153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("inv"), "name", "name", 154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c1_rhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs154), "rhs", "rhs",
                  154);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs154), "lhs", "lhs",
                  154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c1_rhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs155), "rhs", "rhs",
                  155);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs155), "lhs", "lhs",
                  155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c1_rhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs156), "rhs", "rhs",
                  156);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs156), "lhs", "lhs",
                  156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c1_rhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs157), "rhs", "rhs",
                  157);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs157), "lhs", "lhs",
                  157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c1_rhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs158), "rhs", "rhs",
                  158);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs158), "lhs", "lhs",
                  158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("norm"), "name", "name", 159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c1_rhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs159), "rhs", "rhs",
                  159);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs159), "lhs", "lhs",
                  159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c1_rhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs160), "rhs", "rhs",
                  160);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs160), "lhs", "lhs",
                  160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c1_rhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs161), "rhs", "rhs",
                  161);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs161), "lhs", "lhs",
                  161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c1_rhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs162), "rhs", "rhs",
                  162);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs162), "lhs", "lhs",
                  162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c1_rhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs163), "rhs", "rhs",
                  163);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs163), "lhs", "lhs",
                  163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_warning"), "name", "name",
                  164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 164);
  sf_mex_assign(&c1_rhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs164), "rhs", "rhs",
                  164);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs164), "lhs", "lhs",
                  164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 165);
  sf_mex_assign(&c1_rhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs165), "rhs", "rhs",
                  165);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs165), "lhs", "lhs",
                  165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 166);
  sf_mex_assign(&c1_rhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs166), "rhs", "rhs",
                  166);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs166), "lhs", "lhs",
                  166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_flt2str"), "name", "name",
                  167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "resolved",
                  "resolved", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 167);
  sf_mex_assign(&c1_rhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs167), "rhs", "rhs",
                  167);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs167), "lhs", "lhs",
                  167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "context",
                  "context", 168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "name", "name", 168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m"), "resolved",
                  "resolved", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1319737168U), "fileTimeLo",
                  "fileTimeLo", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 168);
  sf_mex_assign(&c1_rhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs168), "rhs", "rhs",
                  168);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs168), "lhs", "lhs",
                  168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/Iseberg/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1447321639U), "fileTimeLo",
                  "fileTimeLo", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 169);
  sf_mex_assign(&c1_rhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs169), "rhs", "rhs",
                  169);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs169), "lhs", "lhs",
                  169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("norm"), "name", "name", 170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 170);
  sf_mex_assign(&c1_rhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs170), "rhs", "rhs",
                  170);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs170), "lhs", "lhs",
                  170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 171);
  sf_mex_assign(&c1_rhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs171), "rhs", "rhs",
                  171);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs171), "lhs", "lhs",
                  171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 172);
  sf_mex_assign(&c1_rhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs172), "rhs", "rhs",
                  172);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs172), "lhs", "lhs",
                  172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 173);
  sf_mex_assign(&c1_rhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs173), "rhs", "rhs",
                  173);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs173), "lhs", "lhs",
                  173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 174);
  sf_mex_assign(&c1_rhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs174), "rhs", "rhs",
                  174);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs174), "lhs", "lhs",
                  174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xnrm2"),
                  "name", "name", 175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "resolved", "resolved", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 175);
  sf_mex_assign(&c1_rhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs175), "rhs", "rhs",
                  175);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs175), "lhs", "lhs",
                  175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 176);
  sf_mex_assign(&c1_rhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs176), "rhs", "rhs",
                  176);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs176), "lhs", "lhs",
                  176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p!below_threshold"),
                  "context", "context", 177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 177);
  sf_mex_assign(&c1_rhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs177), "rhs", "rhs",
                  177);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs177), "lhs", "lhs",
                  177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 178);
  sf_mex_assign(&c1_rhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs178), "rhs", "rhs",
                  178);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs178), "lhs", "lhs",
                  178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xnrm2"),
                  "name", "name", 179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "resolved", "resolved", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 179);
  sf_mex_assign(&c1_rhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs179), "rhs", "rhs",
                  179);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs179), "lhs", "lhs",
                  179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmin"), "name", "name", 180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 180);
  sf_mex_assign(&c1_rhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs180), "rhs", "rhs",
                  180);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs180), "lhs", "lhs",
                  180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 181);
  sf_mex_assign(&c1_rhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs181), "rhs", "rhs",
                  181);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs181), "lhs", "lhs",
                  181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 182);
  sf_mex_assign(&c1_rhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs182), "rhs", "rhs",
                  182);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs182), "lhs", "lhs",
                  182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 183);
  sf_mex_assign(&c1_rhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs183), "rhs", "rhs",
                  183);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs183), "lhs", "lhs",
                  183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 184);
  sf_mex_assign(&c1_rhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs184), "rhs", "rhs",
                  184);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs184), "lhs", "lhs",
                  184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 185);
  sf_mex_assign(&c1_rhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs185), "rhs", "rhs",
                  185);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs185), "lhs", "lhs",
                  185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name",
                  186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 186);
  sf_mex_assign(&c1_rhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs186), "rhs", "rhs",
                  186);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs186), "lhs", "lhs",
                  186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cos"), "name", "name", 187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837572U), "fileTimeLo",
                  "fileTimeLo", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 187);
  sf_mex_assign(&c1_rhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs187), "rhs", "rhs",
                  187);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs187), "lhs", "lhs",
                  187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825922U), "fileTimeLo",
                  "fileTimeLo", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 188);
  sf_mex_assign(&c1_rhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs188), "rhs", "rhs",
                  188);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs188), "lhs", "lhs",
                  188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sin"), "name", "name", 189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 189);
  sf_mex_assign(&c1_rhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs189), "rhs", "rhs",
                  189);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs189), "lhs", "lhs",
                  189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825936U), "fileTimeLo",
                  "fileTimeLo", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 190);
  sf_mex_assign(&c1_rhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs190), "rhs", "rhs",
                  190);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs190), "lhs", "lhs",
                  190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 191);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("exp"), "name", "name", 191);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 191);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m"), "resolved",
                  "resolved", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837580U), "fileTimeLo",
                  "fileTimeLo", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 191);
  sf_mex_assign(&c1_rhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs191), "rhs", "rhs",
                  191);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs191), "lhs", "lhs",
                  191);
  sf_mex_destroy(&c1_rhs128);
  sf_mex_destroy(&c1_lhs128);
  sf_mex_destroy(&c1_rhs129);
  sf_mex_destroy(&c1_lhs129);
  sf_mex_destroy(&c1_rhs130);
  sf_mex_destroy(&c1_lhs130);
  sf_mex_destroy(&c1_rhs131);
  sf_mex_destroy(&c1_lhs131);
  sf_mex_destroy(&c1_rhs132);
  sf_mex_destroy(&c1_lhs132);
  sf_mex_destroy(&c1_rhs133);
  sf_mex_destroy(&c1_lhs133);
  sf_mex_destroy(&c1_rhs134);
  sf_mex_destroy(&c1_lhs134);
  sf_mex_destroy(&c1_rhs135);
  sf_mex_destroy(&c1_lhs135);
  sf_mex_destroy(&c1_rhs136);
  sf_mex_destroy(&c1_lhs136);
  sf_mex_destroy(&c1_rhs137);
  sf_mex_destroy(&c1_lhs137);
  sf_mex_destroy(&c1_rhs138);
  sf_mex_destroy(&c1_lhs138);
  sf_mex_destroy(&c1_rhs139);
  sf_mex_destroy(&c1_lhs139);
  sf_mex_destroy(&c1_rhs140);
  sf_mex_destroy(&c1_lhs140);
  sf_mex_destroy(&c1_rhs141);
  sf_mex_destroy(&c1_lhs141);
  sf_mex_destroy(&c1_rhs142);
  sf_mex_destroy(&c1_lhs142);
  sf_mex_destroy(&c1_rhs143);
  sf_mex_destroy(&c1_lhs143);
  sf_mex_destroy(&c1_rhs144);
  sf_mex_destroy(&c1_lhs144);
  sf_mex_destroy(&c1_rhs145);
  sf_mex_destroy(&c1_lhs145);
  sf_mex_destroy(&c1_rhs146);
  sf_mex_destroy(&c1_lhs146);
  sf_mex_destroy(&c1_rhs147);
  sf_mex_destroy(&c1_lhs147);
  sf_mex_destroy(&c1_rhs148);
  sf_mex_destroy(&c1_lhs148);
  sf_mex_destroy(&c1_rhs149);
  sf_mex_destroy(&c1_lhs149);
  sf_mex_destroy(&c1_rhs150);
  sf_mex_destroy(&c1_lhs150);
  sf_mex_destroy(&c1_rhs151);
  sf_mex_destroy(&c1_lhs151);
  sf_mex_destroy(&c1_rhs152);
  sf_mex_destroy(&c1_lhs152);
  sf_mex_destroy(&c1_rhs153);
  sf_mex_destroy(&c1_lhs153);
  sf_mex_destroy(&c1_rhs154);
  sf_mex_destroy(&c1_lhs154);
  sf_mex_destroy(&c1_rhs155);
  sf_mex_destroy(&c1_lhs155);
  sf_mex_destroy(&c1_rhs156);
  sf_mex_destroy(&c1_lhs156);
  sf_mex_destroy(&c1_rhs157);
  sf_mex_destroy(&c1_lhs157);
  sf_mex_destroy(&c1_rhs158);
  sf_mex_destroy(&c1_lhs158);
  sf_mex_destroy(&c1_rhs159);
  sf_mex_destroy(&c1_lhs159);
  sf_mex_destroy(&c1_rhs160);
  sf_mex_destroy(&c1_lhs160);
  sf_mex_destroy(&c1_rhs161);
  sf_mex_destroy(&c1_lhs161);
  sf_mex_destroy(&c1_rhs162);
  sf_mex_destroy(&c1_lhs162);
  sf_mex_destroy(&c1_rhs163);
  sf_mex_destroy(&c1_lhs163);
  sf_mex_destroy(&c1_rhs164);
  sf_mex_destroy(&c1_lhs164);
  sf_mex_destroy(&c1_rhs165);
  sf_mex_destroy(&c1_lhs165);
  sf_mex_destroy(&c1_rhs166);
  sf_mex_destroy(&c1_lhs166);
  sf_mex_destroy(&c1_rhs167);
  sf_mex_destroy(&c1_lhs167);
  sf_mex_destroy(&c1_rhs168);
  sf_mex_destroy(&c1_lhs168);
  sf_mex_destroy(&c1_rhs169);
  sf_mex_destroy(&c1_lhs169);
  sf_mex_destroy(&c1_rhs170);
  sf_mex_destroy(&c1_lhs170);
  sf_mex_destroy(&c1_rhs171);
  sf_mex_destroy(&c1_lhs171);
  sf_mex_destroy(&c1_rhs172);
  sf_mex_destroy(&c1_lhs172);
  sf_mex_destroy(&c1_rhs173);
  sf_mex_destroy(&c1_lhs173);
  sf_mex_destroy(&c1_rhs174);
  sf_mex_destroy(&c1_lhs174);
  sf_mex_destroy(&c1_rhs175);
  sf_mex_destroy(&c1_lhs175);
  sf_mex_destroy(&c1_rhs176);
  sf_mex_destroy(&c1_lhs176);
  sf_mex_destroy(&c1_rhs177);
  sf_mex_destroy(&c1_lhs177);
  sf_mex_destroy(&c1_rhs178);
  sf_mex_destroy(&c1_lhs178);
  sf_mex_destroy(&c1_rhs179);
  sf_mex_destroy(&c1_lhs179);
  sf_mex_destroy(&c1_rhs180);
  sf_mex_destroy(&c1_lhs180);
  sf_mex_destroy(&c1_rhs181);
  sf_mex_destroy(&c1_lhs181);
  sf_mex_destroy(&c1_rhs182);
  sf_mex_destroy(&c1_lhs182);
  sf_mex_destroy(&c1_rhs183);
  sf_mex_destroy(&c1_lhs183);
  sf_mex_destroy(&c1_rhs184);
  sf_mex_destroy(&c1_lhs184);
  sf_mex_destroy(&c1_rhs185);
  sf_mex_destroy(&c1_lhs185);
  sf_mex_destroy(&c1_rhs186);
  sf_mex_destroy(&c1_lhs186);
  sf_mex_destroy(&c1_rhs187);
  sf_mex_destroy(&c1_lhs187);
  sf_mex_destroy(&c1_rhs188);
  sf_mex_destroy(&c1_lhs188);
  sf_mex_destroy(&c1_rhs189);
  sf_mex_destroy(&c1_lhs189);
  sf_mex_destroy(&c1_rhs190);
  sf_mex_destroy(&c1_lhs190);
  sf_mex_destroy(&c1_rhs191);
  sf_mex_destroy(&c1_lhs191);
}

static void c1_d_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs192 = NULL;
  const mxArray *c1_lhs192 = NULL;
  const mxArray *c1_rhs193 = NULL;
  const mxArray *c1_lhs193 = NULL;
  const mxArray *c1_rhs194 = NULL;
  const mxArray *c1_lhs194 = NULL;
  const mxArray *c1_rhs195 = NULL;
  const mxArray *c1_lhs195 = NULL;
  const mxArray *c1_rhs196 = NULL;
  const mxArray *c1_lhs196 = NULL;
  const mxArray *c1_rhs197 = NULL;
  const mxArray *c1_lhs197 = NULL;
  const mxArray *c1_rhs198 = NULL;
  const mxArray *c1_lhs198 = NULL;
  const mxArray *c1_rhs199 = NULL;
  const mxArray *c1_lhs199 = NULL;
  const mxArray *c1_rhs200 = NULL;
  const mxArray *c1_lhs200 = NULL;
  const mxArray *c1_rhs201 = NULL;
  const mxArray *c1_lhs201 = NULL;
  const mxArray *c1_rhs202 = NULL;
  const mxArray *c1_lhs202 = NULL;
  const mxArray *c1_rhs203 = NULL;
  const mxArray *c1_lhs203 = NULL;
  const mxArray *c1_rhs204 = NULL;
  const mxArray *c1_lhs204 = NULL;
  const mxArray *c1_rhs205 = NULL;
  const mxArray *c1_lhs205 = NULL;
  const mxArray *c1_rhs206 = NULL;
  const mxArray *c1_lhs206 = NULL;
  const mxArray *c1_rhs207 = NULL;
  const mxArray *c1_lhs207 = NULL;
  const mxArray *c1_rhs208 = NULL;
  const mxArray *c1_lhs208 = NULL;
  const mxArray *c1_rhs209 = NULL;
  const mxArray *c1_lhs209 = NULL;
  const mxArray *c1_rhs210 = NULL;
  const mxArray *c1_lhs210 = NULL;
  const mxArray *c1_rhs211 = NULL;
  const mxArray *c1_lhs211 = NULL;
  const mxArray *c1_rhs212 = NULL;
  const mxArray *c1_lhs212 = NULL;
  const mxArray *c1_rhs213 = NULL;
  const mxArray *c1_lhs213 = NULL;
  const mxArray *c1_rhs214 = NULL;
  const mxArray *c1_lhs214 = NULL;
  const mxArray *c1_rhs215 = NULL;
  const mxArray *c1_lhs215 = NULL;
  const mxArray *c1_rhs216 = NULL;
  const mxArray *c1_lhs216 = NULL;
  const mxArray *c1_rhs217 = NULL;
  const mxArray *c1_lhs217 = NULL;
  const mxArray *c1_rhs218 = NULL;
  const mxArray *c1_lhs218 = NULL;
  const mxArray *c1_rhs219 = NULL;
  const mxArray *c1_lhs219 = NULL;
  const mxArray *c1_rhs220 = NULL;
  const mxArray *c1_lhs220 = NULL;
  const mxArray *c1_rhs221 = NULL;
  const mxArray *c1_lhs221 = NULL;
  const mxArray *c1_rhs222 = NULL;
  const mxArray *c1_lhs222 = NULL;
  const mxArray *c1_rhs223 = NULL;
  const mxArray *c1_lhs223 = NULL;
  const mxArray *c1_rhs224 = NULL;
  const mxArray *c1_lhs224 = NULL;
  const mxArray *c1_rhs225 = NULL;
  const mxArray *c1_lhs225 = NULL;
  const mxArray *c1_rhs226 = NULL;
  const mxArray *c1_lhs226 = NULL;
  const mxArray *c1_rhs227 = NULL;
  const mxArray *c1_lhs227 = NULL;
  const mxArray *c1_rhs228 = NULL;
  const mxArray *c1_lhs228 = NULL;
  const mxArray *c1_rhs229 = NULL;
  const mxArray *c1_lhs229 = NULL;
  const mxArray *c1_rhs230 = NULL;
  const mxArray *c1_lhs230 = NULL;
  const mxArray *c1_rhs231 = NULL;
  const mxArray *c1_lhs231 = NULL;
  const mxArray *c1_rhs232 = NULL;
  const mxArray *c1_lhs232 = NULL;
  const mxArray *c1_rhs233 = NULL;
  const mxArray *c1_lhs233 = NULL;
  const mxArray *c1_rhs234 = NULL;
  const mxArray *c1_lhs234 = NULL;
  const mxArray *c1_rhs235 = NULL;
  const mxArray *c1_lhs235 = NULL;
  const mxArray *c1_rhs236 = NULL;
  const mxArray *c1_lhs236 = NULL;
  const mxArray *c1_rhs237 = NULL;
  const mxArray *c1_lhs237 = NULL;
  const mxArray *c1_rhs238 = NULL;
  const mxArray *c1_lhs238 = NULL;
  const mxArray *c1_rhs239 = NULL;
  const mxArray *c1_lhs239 = NULL;
  const mxArray *c1_rhs240 = NULL;
  const mxArray *c1_lhs240 = NULL;
  const mxArray *c1_rhs241 = NULL;
  const mxArray *c1_lhs241 = NULL;
  const mxArray *c1_rhs242 = NULL;
  const mxArray *c1_lhs242 = NULL;
  const mxArray *c1_rhs243 = NULL;
  const mxArray *c1_lhs243 = NULL;
  const mxArray *c1_rhs244 = NULL;
  const mxArray *c1_lhs244 = NULL;
  const mxArray *c1_rhs245 = NULL;
  const mxArray *c1_lhs245 = NULL;
  const mxArray *c1_rhs246 = NULL;
  const mxArray *c1_lhs246 = NULL;
  const mxArray *c1_rhs247 = NULL;
  const mxArray *c1_lhs247 = NULL;
  const mxArray *c1_rhs248 = NULL;
  const mxArray *c1_lhs248 = NULL;
  const mxArray *c1_rhs249 = NULL;
  const mxArray *c1_lhs249 = NULL;
  const mxArray *c1_rhs250 = NULL;
  const mxArray *c1_lhs250 = NULL;
  const mxArray *c1_rhs251 = NULL;
  const mxArray *c1_lhs251 = NULL;
  const mxArray *c1_rhs252 = NULL;
  const mxArray *c1_lhs252 = NULL;
  const mxArray *c1_rhs253 = NULL;
  const mxArray *c1_lhs253 = NULL;
  const mxArray *c1_rhs254 = NULL;
  const mxArray *c1_lhs254 = NULL;
  const mxArray *c1_rhs255 = NULL;
  const mxArray *c1_lhs255 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m"), "context",
                  "context", 192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_exp"), "name",
                  "name", 192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_exp.m"),
                  "resolved", "resolved", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1301335664U), "fileTimeLo",
                  "fileTimeLo", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 192);
  sf_mex_assign(&c1_rhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs192), "rhs", "rhs",
                  192);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs192), "lhs", "lhs",
                  192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 193);
  sf_mex_assign(&c1_rhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs193), "rhs", "rhs",
                  193);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs193), "lhs", "lhs",
                  193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 194);
  sf_mex_assign(&c1_rhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs194), "rhs", "rhs",
                  194);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs194), "lhs", "lhs",
                  194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 195);
  sf_mex_assign(&c1_rhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs195), "rhs", "rhs",
                  195);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs195), "lhs", "lhs",
                  195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 196);
  sf_mex_assign(&c1_rhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs196), "rhs", "rhs",
                  196);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs196), "lhs", "lhs",
                  196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 197);
  sf_mex_assign(&c1_rhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs197), "rhs", "rhs",
                  197);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs197), "lhs", "lhs",
                  197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 198);
  sf_mex_assign(&c1_rhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs198), "rhs", "rhs",
                  198);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs198), "lhs", "lhs",
                  198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 199);
  sf_mex_assign(&c1_rhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs199), "rhs", "rhs",
                  199);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs199), "lhs", "lhs",
                  199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 200);
  sf_mex_assign(&c1_rhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs200), "rhs", "rhs",
                  200);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs200), "lhs", "lhs",
                  200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 201);
  sf_mex_assign(&c1_rhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs201), "rhs", "rhs",
                  201);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs201), "lhs", "lhs",
                  201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 202);
  sf_mex_assign(&c1_rhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs202), "rhs", "rhs",
                  202);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs202), "lhs", "lhs",
                  202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 203);
  sf_mex_assign(&c1_rhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs203), "rhs", "rhs",
                  203);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs203), "lhs", "lhs",
                  203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 204);
  sf_mex_assign(&c1_rhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs204), "rhs", "rhs",
                  204);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs204), "lhs", "lhs",
                  204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 205);
  sf_mex_assign(&c1_rhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs205), "rhs", "rhs",
                  205);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs205), "lhs", "lhs",
                  205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 206);
  sf_mex_assign(&c1_rhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs206), "rhs", "rhs",
                  206);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs206), "lhs", "lhs",
                  206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("inv"), "name", "name", 207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 207);
  sf_mex_assign(&c1_rhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs207), "rhs", "rhs",
                  207);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs207), "lhs", "lhs",
                  207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 208);
  sf_mex_assign(&c1_rhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs208), "rhs", "rhs",
                  208);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs208), "lhs", "lhs",
                  208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_error"), "name", "name",
                  209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 209);
  sf_mex_assign(&c1_rhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs209), "rhs", "rhs",
                  209);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs209), "lhs", "lhs",
                  209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 210);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eig"), "name", "name", 210);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 210);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "resolved",
                  "resolved", 210);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 210);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 210);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 210);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 210);
  sf_mex_assign(&c1_rhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs210), "rhs", "rhs",
                  210);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs210), "lhs", "lhs",
                  210);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 211);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 211);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 211);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 211);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 211);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 211);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 211);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 211);
  sf_mex_assign(&c1_rhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs211), "rhs", "rhs",
                  211);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs211), "lhs", "lhs",
                  211);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 212);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eye"), "name", "name", 212);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 212);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 212);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 212);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 212);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 212);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 212);
  sf_mex_assign(&c1_rhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs212), "rhs", "rhs",
                  212);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs212), "lhs", "lhs",
                  212);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 213);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 213);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 213);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 213);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 213);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 213);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 213);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 213);
  sf_mex_assign(&c1_rhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs213), "rhs", "rhs",
                  213);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs213), "lhs", "lhs",
                  213);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 214);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 214);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 214);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 214);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 214);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 214);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 214);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 214);
  sf_mex_assign(&c1_rhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs214), "rhs", "rhs",
                  214);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs214), "lhs", "lhs",
                  214);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 215);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 215);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 215);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 215);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 215);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 215);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 215);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 215);
  sf_mex_assign(&c1_rhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs215), "rhs", "rhs",
                  215);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs215), "lhs", "lhs",
                  215);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 216);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 216);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 216);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 216);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 216);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 216);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 216);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 216);
  sf_mex_assign(&c1_rhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs216), "rhs", "rhs",
                  216);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs216), "lhs", "lhs",
                  216);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 217);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 217);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 217);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 217);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 217);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 217);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 217);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 217);
  sf_mex_assign(&c1_rhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs217), "rhs", "rhs",
                  217);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs217), "lhs", "lhs",
                  217);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 218);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 218);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 218);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 218);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 218);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 218);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 218);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 218);
  sf_mex_assign(&c1_rhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs218), "rhs", "rhs",
                  218);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs218), "lhs", "lhs",
                  218);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 219);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 219);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 219);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 219);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 219);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 219);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 219);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 219);
  sf_mex_assign(&c1_rhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs219), "rhs", "rhs",
                  219);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs219), "lhs", "lhs",
                  219);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 220);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 220);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 220);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 220);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 220);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 220);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 220);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 220);
  sf_mex_assign(&c1_rhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs220), "rhs", "rhs",
                  220);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs220), "lhs", "lhs",
                  220);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 221);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 221);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 221);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 221);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 221);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 221);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 221);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 221);
  sf_mex_assign(&c1_rhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs221), "rhs", "rhs",
                  221);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs221), "lhs", "lhs",
                  221);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 222);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 222);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 222);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 222);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 222);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 222);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 222);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 222);
  sf_mex_assign(&c1_rhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs222), "rhs", "rhs",
                  222);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs222), "lhs", "lhs",
                  222);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 223);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_ztgevc"), "name",
                  "name", 223);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 223);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "resolved", "resolved", 223);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826024U), "fileTimeLo",
                  "fileTimeLo", 223);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 223);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 223);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 223);
  sf_mex_assign(&c1_rhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs223), "rhs", "rhs",
                  223);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs223), "lhs", "lhs",
                  223);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 224);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 224);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 224);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 224);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 224);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 224);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 224);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 224);
  sf_mex_assign(&c1_rhs224, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs224, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs224), "rhs", "rhs",
                  224);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs224), "lhs", "lhs",
                  224);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 225);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmin"), "name", "name", 225);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 225);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 225);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 225);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 225);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 225);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 225);
  sf_mex_assign(&c1_rhs225, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs225, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs225), "rhs", "rhs",
                  225);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs225), "lhs", "lhs",
                  225);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m!abs1"),
                  "context", "context", 226);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 226);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 226);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 226);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 226);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 226);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 226);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 226);
  sf_mex_assign(&c1_rhs226, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs226, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs226), "rhs", "rhs",
                  226);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs226), "lhs", "lhs",
                  226);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 227);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 227);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 227);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 227);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 227);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 227);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 227);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 227);
  sf_mex_assign(&c1_rhs227, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs227, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs227), "rhs", "rhs",
                  227);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs227), "lhs", "lhs",
                  227);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 228);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 228);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 228);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 228);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 228);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 228);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 228);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 228);
  sf_mex_assign(&c1_rhs228, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs228, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs228), "rhs", "rhs",
                  228);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs228), "lhs", "lhs",
                  228);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 229);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name",
                  229);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 229);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 229);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 229);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 229);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 229);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 229);
  sf_mex_assign(&c1_rhs229, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs229, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs229), "rhs", "rhs",
                  229);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs229), "lhs", "lhs",
                  229);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 230);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zggbak"), "name",
                  "name", 230);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 230);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "resolved", "resolved", 230);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826018U), "fileTimeLo",
                  "fileTimeLo", 230);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 230);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 230);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 230);
  sf_mex_assign(&c1_rhs230, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs230, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs230), "rhs", "rhs",
                  230);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs230), "lhs", "lhs",
                  230);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 231);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 231);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 231);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 231);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 231);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 231);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 231);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 231);
  sf_mex_assign(&c1_rhs231, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs231, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs231), "rhs", "rhs",
                  231);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs231), "lhs", "lhs",
                  231);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 232);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 232);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 232);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 232);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 232);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 232);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 232);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 232);
  sf_mex_assign(&c1_rhs232, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs232, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs232), "rhs", "rhs",
                  232);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs232), "lhs", "lhs",
                  232);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 233);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 233);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 233);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 233);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 233);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 233);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 233);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 233);
  sf_mex_assign(&c1_rhs233, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs233, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs233), "rhs", "rhs",
                  233);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs233), "lhs", "lhs",
                  233);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 234);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 234);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 234);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 234);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 234);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 234);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 234);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 234);
  sf_mex_assign(&c1_rhs234, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs234, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs234), "rhs", "rhs",
                  234);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs234), "lhs", "lhs",
                  234);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m!abs1"),
                  "context", "context", 235);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 235);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 235);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 235);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 235);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 235);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 235);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 235);
  sf_mex_assign(&c1_rhs235, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs235, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs235), "rhs", "rhs",
                  235);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs235), "lhs", "lhs",
                  235);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 236);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 236);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 236);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 236);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 236);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 236);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 236);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 236);
  sf_mex_assign(&c1_rhs236, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs236, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs236), "rhs", "rhs",
                  236);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs236), "lhs", "lhs",
                  236);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 237);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 237);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 237);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 237);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 237);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 237);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 237);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 237);
  sf_mex_assign(&c1_rhs237, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs237, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs237), "rhs", "rhs",
                  237);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs237), "lhs", "lhs",
                  237);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 238);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 238);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 238);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 238);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 238);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 238);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 238);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 238);
  sf_mex_assign(&c1_rhs238, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs238, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs238), "rhs", "rhs",
                  238);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs238), "lhs", "lhs",
                  238);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 239);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 239);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 239);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 239);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 239);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 239);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 239);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 239);
  sf_mex_assign(&c1_rhs239, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs239, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs239), "rhs", "rhs",
                  239);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs239), "lhs", "lhs",
                  239);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 240);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 240);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 240);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 240);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 240);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 240);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 240);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 240);
  sf_mex_assign(&c1_rhs240, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs240, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs240), "rhs", "rhs",
                  240);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs240), "lhs", "lhs",
                  240);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 241);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 241);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 241);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 241);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 241);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 241);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 241);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 241);
  sf_mex_assign(&c1_rhs241, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs241, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs241), "rhs", "rhs",
                  241);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs241), "lhs", "lhs",
                  241);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 242);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 242);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 242);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 242);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 242);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 242);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 242);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 242);
  sf_mex_assign(&c1_rhs242, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs242, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs242), "rhs", "rhs",
                  242);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs242), "lhs", "lhs",
                  242);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 243);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 243);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 243);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 243);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 243);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 243);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 243);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 243);
  sf_mex_assign(&c1_rhs243, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs243, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs243), "rhs", "rhs",
                  243);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs243), "lhs", "lhs",
                  243);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 244);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  244);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 244);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 244);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 244);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 244);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 244);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 244);
  sf_mex_assign(&c1_rhs244, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs244, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs244), "rhs", "rhs",
                  244);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs244), "lhs", "lhs",
                  244);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 245);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name",
                  245);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 245);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 245);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 245);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 245);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 245);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 245);
  sf_mex_assign(&c1_rhs245, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs245, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs245), "rhs", "rhs",
                  245);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs245), "lhs", "lhs",
                  245);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 246);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("diag"), "name", "name", 246);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 246);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "resolved",
                  "resolved", 246);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 246);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 246);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 246);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 246);
  sf_mex_assign(&c1_rhs246, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs246, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs246), "rhs", "rhs",
                  246);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs246), "lhs", "lhs",
                  246);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 247);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("ismatrix"), "name", "name",
                  247);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 247);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 247);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 247);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 247);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 247);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 247);
  sf_mex_assign(&c1_rhs247, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs247, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs247), "rhs", "rhs",
                  247);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs247), "lhs", "lhs",
                  247);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 248);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 248);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 248);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 248);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 248);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 248);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 248);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 248);
  sf_mex_assign(&c1_rhs248, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs248, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs248), "rhs", "rhs",
                  248);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs248), "lhs", "lhs",
                  248);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 249);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 249);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 249);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 249);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 249);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 249);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 249);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 249);
  sf_mex_assign(&c1_rhs249, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs249, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs249), "rhs", "rhs",
                  249);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs249), "lhs", "lhs",
                  249);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 250);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 250);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 250);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 250);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 250);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 250);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 250);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 250);
  sf_mex_assign(&c1_rhs250, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs250, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs250), "rhs", "rhs",
                  250);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs250), "lhs", "lhs",
                  250);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 251);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("power"), "name", "name", 251);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 251);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 251);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 251);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 251);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 251);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 251);
  sf_mex_assign(&c1_rhs251, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs251, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs251), "rhs", "rhs",
                  251);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs251), "lhs", "lhs",
                  251);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 252);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("flintmax"), "name", "name",
                  252);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 252);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/flintmax.m"), "resolved",
                  "resolved", 252);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1348199116U), "fileTimeLo",
                  "fileTimeLo", 252);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 252);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 252);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 252);
  sf_mex_assign(&c1_rhs252, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs252, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs252), "rhs", "rhs",
                  252);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs252), "lhs", "lhs",
                  252);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/flintmax.m"), "context",
                  "context", 253);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 253);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 253);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 253);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 253);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 253);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 253);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 253);
  sf_mex_assign(&c1_rhs253, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs253, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs253), "rhs", "rhs",
                  253);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs253), "lhs", "lhs",
                  253);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 254);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mod"), "name", "name", 254);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 254);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "resolved",
                  "resolved", 254);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 254);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 254);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 254);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 254);
  sf_mex_assign(&c1_rhs254, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs254, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs254), "rhs", "rhs",
                  254);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs254), "lhs", "lhs",
                  254);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 255);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 255);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 255);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 255);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 255);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 255);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 255);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 255);
  sf_mex_assign(&c1_rhs255, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs255, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs255), "rhs", "rhs",
                  255);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs255), "lhs", "lhs",
                  255);
  sf_mex_destroy(&c1_rhs192);
  sf_mex_destroy(&c1_lhs192);
  sf_mex_destroy(&c1_rhs193);
  sf_mex_destroy(&c1_lhs193);
  sf_mex_destroy(&c1_rhs194);
  sf_mex_destroy(&c1_lhs194);
  sf_mex_destroy(&c1_rhs195);
  sf_mex_destroy(&c1_lhs195);
  sf_mex_destroy(&c1_rhs196);
  sf_mex_destroy(&c1_lhs196);
  sf_mex_destroy(&c1_rhs197);
  sf_mex_destroy(&c1_lhs197);
  sf_mex_destroy(&c1_rhs198);
  sf_mex_destroy(&c1_lhs198);
  sf_mex_destroy(&c1_rhs199);
  sf_mex_destroy(&c1_lhs199);
  sf_mex_destroy(&c1_rhs200);
  sf_mex_destroy(&c1_lhs200);
  sf_mex_destroy(&c1_rhs201);
  sf_mex_destroy(&c1_lhs201);
  sf_mex_destroy(&c1_rhs202);
  sf_mex_destroy(&c1_lhs202);
  sf_mex_destroy(&c1_rhs203);
  sf_mex_destroy(&c1_lhs203);
  sf_mex_destroy(&c1_rhs204);
  sf_mex_destroy(&c1_lhs204);
  sf_mex_destroy(&c1_rhs205);
  sf_mex_destroy(&c1_lhs205);
  sf_mex_destroy(&c1_rhs206);
  sf_mex_destroy(&c1_lhs206);
  sf_mex_destroy(&c1_rhs207);
  sf_mex_destroy(&c1_lhs207);
  sf_mex_destroy(&c1_rhs208);
  sf_mex_destroy(&c1_lhs208);
  sf_mex_destroy(&c1_rhs209);
  sf_mex_destroy(&c1_lhs209);
  sf_mex_destroy(&c1_rhs210);
  sf_mex_destroy(&c1_lhs210);
  sf_mex_destroy(&c1_rhs211);
  sf_mex_destroy(&c1_lhs211);
  sf_mex_destroy(&c1_rhs212);
  sf_mex_destroy(&c1_lhs212);
  sf_mex_destroy(&c1_rhs213);
  sf_mex_destroy(&c1_lhs213);
  sf_mex_destroy(&c1_rhs214);
  sf_mex_destroy(&c1_lhs214);
  sf_mex_destroy(&c1_rhs215);
  sf_mex_destroy(&c1_lhs215);
  sf_mex_destroy(&c1_rhs216);
  sf_mex_destroy(&c1_lhs216);
  sf_mex_destroy(&c1_rhs217);
  sf_mex_destroy(&c1_lhs217);
  sf_mex_destroy(&c1_rhs218);
  sf_mex_destroy(&c1_lhs218);
  sf_mex_destroy(&c1_rhs219);
  sf_mex_destroy(&c1_lhs219);
  sf_mex_destroy(&c1_rhs220);
  sf_mex_destroy(&c1_lhs220);
  sf_mex_destroy(&c1_rhs221);
  sf_mex_destroy(&c1_lhs221);
  sf_mex_destroy(&c1_rhs222);
  sf_mex_destroy(&c1_lhs222);
  sf_mex_destroy(&c1_rhs223);
  sf_mex_destroy(&c1_lhs223);
  sf_mex_destroy(&c1_rhs224);
  sf_mex_destroy(&c1_lhs224);
  sf_mex_destroy(&c1_rhs225);
  sf_mex_destroy(&c1_lhs225);
  sf_mex_destroy(&c1_rhs226);
  sf_mex_destroy(&c1_lhs226);
  sf_mex_destroy(&c1_rhs227);
  sf_mex_destroy(&c1_lhs227);
  sf_mex_destroy(&c1_rhs228);
  sf_mex_destroy(&c1_lhs228);
  sf_mex_destroy(&c1_rhs229);
  sf_mex_destroy(&c1_lhs229);
  sf_mex_destroy(&c1_rhs230);
  sf_mex_destroy(&c1_lhs230);
  sf_mex_destroy(&c1_rhs231);
  sf_mex_destroy(&c1_lhs231);
  sf_mex_destroy(&c1_rhs232);
  sf_mex_destroy(&c1_lhs232);
  sf_mex_destroy(&c1_rhs233);
  sf_mex_destroy(&c1_lhs233);
  sf_mex_destroy(&c1_rhs234);
  sf_mex_destroy(&c1_lhs234);
  sf_mex_destroy(&c1_rhs235);
  sf_mex_destroy(&c1_lhs235);
  sf_mex_destroy(&c1_rhs236);
  sf_mex_destroy(&c1_lhs236);
  sf_mex_destroy(&c1_rhs237);
  sf_mex_destroy(&c1_lhs237);
  sf_mex_destroy(&c1_rhs238);
  sf_mex_destroy(&c1_lhs238);
  sf_mex_destroy(&c1_rhs239);
  sf_mex_destroy(&c1_lhs239);
  sf_mex_destroy(&c1_rhs240);
  sf_mex_destroy(&c1_lhs240);
  sf_mex_destroy(&c1_rhs241);
  sf_mex_destroy(&c1_lhs241);
  sf_mex_destroy(&c1_rhs242);
  sf_mex_destroy(&c1_lhs242);
  sf_mex_destroy(&c1_rhs243);
  sf_mex_destroy(&c1_lhs243);
  sf_mex_destroy(&c1_rhs244);
  sf_mex_destroy(&c1_lhs244);
  sf_mex_destroy(&c1_rhs245);
  sf_mex_destroy(&c1_lhs245);
  sf_mex_destroy(&c1_rhs246);
  sf_mex_destroy(&c1_lhs246);
  sf_mex_destroy(&c1_rhs247);
  sf_mex_destroy(&c1_lhs247);
  sf_mex_destroy(&c1_rhs248);
  sf_mex_destroy(&c1_lhs248);
  sf_mex_destroy(&c1_rhs249);
  sf_mex_destroy(&c1_lhs249);
  sf_mex_destroy(&c1_rhs250);
  sf_mex_destroy(&c1_lhs250);
  sf_mex_destroy(&c1_rhs251);
  sf_mex_destroy(&c1_lhs251);
  sf_mex_destroy(&c1_rhs252);
  sf_mex_destroy(&c1_lhs252);
  sf_mex_destroy(&c1_rhs253);
  sf_mex_destroy(&c1_lhs253);
  sf_mex_destroy(&c1_rhs254);
  sf_mex_destroy(&c1_lhs254);
  sf_mex_destroy(&c1_rhs255);
  sf_mex_destroy(&c1_lhs255);
}

static void c1_e_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs256 = NULL;
  const mxArray *c1_lhs256 = NULL;
  const mxArray *c1_rhs257 = NULL;
  const mxArray *c1_lhs257 = NULL;
  const mxArray *c1_rhs258 = NULL;
  const mxArray *c1_lhs258 = NULL;
  const mxArray *c1_rhs259 = NULL;
  const mxArray *c1_lhs259 = NULL;
  const mxArray *c1_rhs260 = NULL;
  const mxArray *c1_lhs260 = NULL;
  const mxArray *c1_rhs261 = NULL;
  const mxArray *c1_lhs261 = NULL;
  const mxArray *c1_rhs262 = NULL;
  const mxArray *c1_lhs262 = NULL;
  const mxArray *c1_rhs263 = NULL;
  const mxArray *c1_lhs263 = NULL;
  const mxArray *c1_rhs264 = NULL;
  const mxArray *c1_lhs264 = NULL;
  const mxArray *c1_rhs265 = NULL;
  const mxArray *c1_lhs265 = NULL;
  const mxArray *c1_rhs266 = NULL;
  const mxArray *c1_lhs266 = NULL;
  const mxArray *c1_rhs267 = NULL;
  const mxArray *c1_lhs267 = NULL;
  const mxArray *c1_rhs268 = NULL;
  const mxArray *c1_lhs268 = NULL;
  const mxArray *c1_rhs269 = NULL;
  const mxArray *c1_lhs269 = NULL;
  const mxArray *c1_rhs270 = NULL;
  const mxArray *c1_lhs270 = NULL;
  const mxArray *c1_rhs271 = NULL;
  const mxArray *c1_lhs271 = NULL;
  const mxArray *c1_rhs272 = NULL;
  const mxArray *c1_lhs272 = NULL;
  const mxArray *c1_rhs273 = NULL;
  const mxArray *c1_lhs273 = NULL;
  const mxArray *c1_rhs274 = NULL;
  const mxArray *c1_lhs274 = NULL;
  const mxArray *c1_rhs275 = NULL;
  const mxArray *c1_lhs275 = NULL;
  const mxArray *c1_rhs276 = NULL;
  const mxArray *c1_lhs276 = NULL;
  const mxArray *c1_rhs277 = NULL;
  const mxArray *c1_lhs277 = NULL;
  const mxArray *c1_rhs278 = NULL;
  const mxArray *c1_lhs278 = NULL;
  const mxArray *c1_rhs279 = NULL;
  const mxArray *c1_lhs279 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 256);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 256);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 256);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 256);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 256);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 256);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 256);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 256);
  sf_mex_assign(&c1_rhs256, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs256, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs256), "rhs", "rhs",
                  256);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs256), "lhs", "lhs",
                  256);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 257);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_round"), "name",
                  "name", 257);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 257);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_round.m"),
                  "resolved", "resolved", 257);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658438U), "fileTimeLo",
                  "fileTimeLo", 257);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 257);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 257);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 257);
  sf_mex_assign(&c1_rhs257, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs257, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs257), "rhs", "rhs",
                  257);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs257), "lhs", "lhs",
                  257);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 258);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 258);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 258);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 258);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 258);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 258);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 258);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 258);
  sf_mex_assign(&c1_rhs258, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs258, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs258), "rhs", "rhs",
                  258);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs258), "lhs", "lhs",
                  258);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 259);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 259);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 259);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 259);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 259);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 259);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 259);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 259);
  sf_mex_assign(&c1_rhs259, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs259, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs259), "rhs", "rhs",
                  259);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs259), "lhs", "lhs",
                  259);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 260);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("log"), "name", "name", 260);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 260);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m"), "resolved",
                  "resolved", 260);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837580U), "fileTimeLo",
                  "fileTimeLo", 260);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 260);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 260);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 260);
  sf_mex_assign(&c1_rhs260, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs260, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs260), "rhs", "rhs",
                  260);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs260), "lhs", "lhs",
                  260);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m"), "context",
                  "context", 261);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_log"), "name",
                  "name", 261);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 261);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "resolved", "resolved", 261);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825928U), "fileTimeLo",
                  "fileTimeLo", 261);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 261);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 261);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 261);
  sf_mex_assign(&c1_rhs261, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs261, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs261), "rhs", "rhs",
                  261);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs261), "lhs", "lhs",
                  261);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 262);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmax"), "name", "name", 262);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 262);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m"), "resolved",
                  "resolved", 262);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 262);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 262);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 262);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 262);
  sf_mex_assign(&c1_rhs262, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs262, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs262), "rhs", "rhs",
                  262);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs262), "lhs", "lhs",
                  262);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 263);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name",
                  263);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 263);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 263);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 263);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 263);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 263);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 263);
  sf_mex_assign(&c1_rhs263, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs263, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs263), "rhs", "rhs",
                  263);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs263), "lhs", "lhs",
                  263);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 264);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 264);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 264);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 264);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 264);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 264);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 264);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 264);
  sf_mex_assign(&c1_rhs264, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs264, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs264), "rhs", "rhs",
                  264);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs264), "lhs", "lhs",
                  264);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 265);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_hypot"), "name",
                  "name", 265);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 265);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "resolved", "resolved", 265);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 265);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 265);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 265);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 265);
  sf_mex_assign(&c1_rhs265, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs265, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs265), "rhs", "rhs",
                  265);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs265), "lhs", "lhs",
                  265);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 266);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_atan2"), "name",
                  "name", 266);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 266);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m"),
                  "resolved", "resolved", 266);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825920U), "fileTimeLo",
                  "fileTimeLo", 266);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 266);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 266);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 266);
  sf_mex_assign(&c1_rhs266, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs266, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs266), "rhs", "rhs",
                  266);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs266), "lhs", "lhs",
                  266);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 267);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 267);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 267);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 267);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 267);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 267);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 267);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 267);
  sf_mex_assign(&c1_rhs267, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs267, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs267), "rhs", "rhs",
                  267);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs267), "lhs", "lhs",
                  267);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 268);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 268);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 268);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 268);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 268);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 268);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 268);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 268);
  sf_mex_assign(&c1_rhs268, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs268, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs268), "rhs", "rhs",
                  268);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs268), "lhs", "lhs",
                  268);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 269);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name",
                  269);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 269);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 269);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 269);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 269);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 269);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 269);
  sf_mex_assign(&c1_rhs269, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs269, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs269), "rhs", "rhs",
                  269);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs269), "lhs", "lhs",
                  269);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 270);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("ismatrix"), "name", "name",
                  270);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 270);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 270);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 270);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 270);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 270);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 270);
  sf_mex_assign(&c1_rhs270, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs270, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs270), "rhs", "rhs",
                  270);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs270), "lhs", "lhs",
                  270);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 271);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_lusolve"), "name", "name",
                  271);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 271);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "resolved",
                  "resolved", 271);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "fileTimeLo",
                  "fileTimeLo", 271);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 271);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 271);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 271);
  sf_mex_assign(&c1_rhs271, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs271, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs271), "rhs", "rhs",
                  271);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs271), "lhs", "lhs",
                  271);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 272);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xcabs1"), "name", "name",
                  272);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 272);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "resolved", "resolved", 272);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 272);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 272);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 272);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 272);
  sf_mex_assign(&c1_rhs272, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs272, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs272), "rhs", "rhs",
                  272);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs272), "lhs", "lhs",
                  272);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "context", "context", 273);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xcabs1"),
                  "name", "name", 273);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 273);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "resolved", "resolved", 273);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 273);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 273);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 273);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 273);
  sf_mex_assign(&c1_rhs273, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs273, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs273), "rhs", "rhs",
                  273);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs273), "lhs", "lhs",
                  273);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "context", "context", 274);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 274);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 274);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 274);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 274);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 274);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 274);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 274);
  sf_mex_assign(&c1_rhs274, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs274, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs274), "rhs", "rhs",
                  274);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs274), "lhs", "lhs",
                  274);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 275);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 275);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 275);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 275);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 275);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 275);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 275);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 275);
  sf_mex_assign(&c1_rhs275, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs275, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs275), "rhs", "rhs",
                  275);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs275), "lhs", "lhs",
                  275);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular"),
                  "context", "context", 276);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_warning"), "name", "name",
                  276);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 276);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 276);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 276);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 276);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 276);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 276);
  sf_mex_assign(&c1_rhs276, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs276, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs276), "rhs", "rhs",
                  276);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs276), "lhs", "lhs",
                  276);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 277);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 277);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 277);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 277);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 277);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 277);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 277);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 277);
  sf_mex_assign(&c1_rhs277, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs277, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs277), "rhs", "rhs",
                  277);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs277), "lhs", "lhs",
                  277);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 278);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 278);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 278);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 278);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 278);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 278);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 278);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 278);
  sf_mex_assign(&c1_rhs278, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs278, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs278), "rhs", "rhs",
                  278);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs278), "lhs", "lhs",
                  278);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 279);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 279);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 279);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 279);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 279);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 279);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 279);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 279);
  sf_mex_assign(&c1_rhs279, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs279, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs279), "rhs", "rhs",
                  279);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs279), "lhs", "lhs",
                  279);
  sf_mex_destroy(&c1_rhs256);
  sf_mex_destroy(&c1_lhs256);
  sf_mex_destroy(&c1_rhs257);
  sf_mex_destroy(&c1_lhs257);
  sf_mex_destroy(&c1_rhs258);
  sf_mex_destroy(&c1_lhs258);
  sf_mex_destroy(&c1_rhs259);
  sf_mex_destroy(&c1_lhs259);
  sf_mex_destroy(&c1_rhs260);
  sf_mex_destroy(&c1_lhs260);
  sf_mex_destroy(&c1_rhs261);
  sf_mex_destroy(&c1_lhs261);
  sf_mex_destroy(&c1_rhs262);
  sf_mex_destroy(&c1_lhs262);
  sf_mex_destroy(&c1_rhs263);
  sf_mex_destroy(&c1_lhs263);
  sf_mex_destroy(&c1_rhs264);
  sf_mex_destroy(&c1_lhs264);
  sf_mex_destroy(&c1_rhs265);
  sf_mex_destroy(&c1_lhs265);
  sf_mex_destroy(&c1_rhs266);
  sf_mex_destroy(&c1_lhs266);
  sf_mex_destroy(&c1_rhs267);
  sf_mex_destroy(&c1_lhs267);
  sf_mex_destroy(&c1_rhs268);
  sf_mex_destroy(&c1_lhs268);
  sf_mex_destroy(&c1_rhs269);
  sf_mex_destroy(&c1_lhs269);
  sf_mex_destroy(&c1_rhs270);
  sf_mex_destroy(&c1_lhs270);
  sf_mex_destroy(&c1_rhs271);
  sf_mex_destroy(&c1_lhs271);
  sf_mex_destroy(&c1_rhs272);
  sf_mex_destroy(&c1_lhs272);
  sf_mex_destroy(&c1_rhs273);
  sf_mex_destroy(&c1_lhs273);
  sf_mex_destroy(&c1_rhs274);
  sf_mex_destroy(&c1_lhs274);
  sf_mex_destroy(&c1_rhs275);
  sf_mex_destroy(&c1_lhs275);
  sf_mex_destroy(&c1_rhs276);
  sf_mex_destroy(&c1_lhs276);
  sf_mex_destroy(&c1_rhs277);
  sf_mex_destroy(&c1_lhs277);
  sf_mex_destroy(&c1_rhs278);
  sf_mex_destroy(&c1_lhs278);
  sf_mex_destroy(&c1_rhs279);
  sf_mex_destroy(&c1_lhs279);
}

static void c1_eig(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_A[9],
                   creal_T c1_V[3])
{
  int32_T c1_i39;
  static creal_T c1_dc0 = { 0.0, 0.0 };

  creal_T c1_b_A[9];
  real_T c1_info;
  int32_T c1_i40;
  creal_T c1_c_A[9];
  real_T c1_anrm;
  int32_T c1_i41;
  creal_T c1_alpha1[3];
  int32_T c1_i42;
  creal_T c1_beta1[3];
  boolean_T c1_ilascl;
  real_T c1_anrmto;
  int32_T c1_rscale[3];
  int32_T c1_ihi;
  int32_T c1_ilo;
  int32_T c1_b_ilo;
  int32_T c1_b_ihi;
  int32_T c1_c_ilo;
  int32_T c1_c_ihi;
  int32_T c1_a;
  int32_T c1_b_a;
  int32_T c1_c;
  int32_T c1_c_a;
  int32_T c1_d_a;
  int32_T c1_ihim1;
  int32_T c1_jcol;
  int32_T c1_e_a;
  int32_T c1_f_a;
  int32_T c1_jcolp1;
  int32_T c1_jrow;
  int32_T c1_g_a;
  int32_T c1_h_a;
  int32_T c1_jrowm1;
  creal_T c1_d_A;
  creal_T c1_e_A;
  creal_T c1_b;
  creal_T c1_s;
  real_T c1_b_c;
  real_T c1_c_c;
  real_T c1_d_c;
  int32_T c1_xrow;
  int32_T c1_yrow;
  int32_T c1_jlo;
  int32_T c1_jhi;
  int32_T c1_b_jlo;
  int32_T c1_b_jhi;
  int32_T c1_i_a;
  int32_T c1_b_b;
  int32_T c1_j_a;
  int32_T c1_c_b;
  boolean_T c1_overflow;
  int32_T c1_j;
  int32_T c1_b_j;
  real_T c1_k_a;
  creal_T c1_y;
  creal_T c1_b_s;
  creal_T c1_stemp;
  real_T c1_l_a;
  creal_T c1_d_b;
  creal_T c1_e_b;
  creal_T c1_f_b;
  creal_T c1_g_b;
  real_T c1_e_c;
  int32_T c1_xcol;
  int32_T c1_ycol;
  int32_T c1_d_ilo;
  int32_T c1_d_ihi;
  int32_T c1_e_ilo;
  int32_T c1_e_ihi;
  int32_T c1_m_a;
  int32_T c1_h_b;
  int32_T c1_n_a;
  int32_T c1_i_b;
  boolean_T c1_b_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  real_T c1_o_a;
  creal_T c1_c_s;
  real_T c1_p_a;
  creal_T c1_j_b;
  creal_T c1_k_b;
  creal_T c1_l_b;
  creal_T c1_m_b;
  int32_T c1_i43;
  creal_T c1_f_A[9];
  real_T c1_b_info;
  real_T c1_c_info;
  real_T c1_d_info;
  real_T c1_e_info;
  int32_T c1_i44;
  creal_T c1_b_alpha1[3];
  int32_T c1_i45;
  creal_T c1_b_beta1[3];
  boolean_T guard1 = false;
  for (c1_i39 = 0; c1_i39 < 9; c1_i39++) {
    c1_b_A[c1_i39].re = c1_A[c1_i39] + c1_dc0.re;
    c1_b_A[c1_i39].im = c1_dc0.im;
  }

  c1_info = 0.0;
  c1_realmin(chartInstance);
  c1_eps(chartInstance);
  for (c1_i40 = 0; c1_i40 < 9; c1_i40++) {
    c1_c_A[c1_i40] = c1_b_A[c1_i40];
  }

  c1_anrm = c1_eml_matlab_zlangeM(chartInstance, c1_c_A);
  if (!c1_isfinite(chartInstance, c1_anrm)) {
    for (c1_i41 = 0; c1_i41 < 3; c1_i41++) {
      c1_alpha1[c1_i41].re = rtNaN;
      c1_alpha1[c1_i41].im = 0.0;
    }

    for (c1_i42 = 0; c1_i42 < 3; c1_i42++) {
      c1_beta1[c1_i42].re = rtNaN;
      c1_beta1[c1_i42].im = 0.0;
    }
  } else {
    c1_ilascl = false;
    c1_anrmto = c1_anrm;
    guard1 = false;
    if (c1_anrm > 0.0) {
      if (c1_anrm < 6.7178761075670888E-139) {
        c1_anrmto = 6.7178761075670888E-139;
        c1_ilascl = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      if (c1_anrm > 1.4885657073574029E+138) {
        c1_anrmto = 1.4885657073574029E+138;
        c1_ilascl = true;
      }
    }

    if (c1_ilascl) {
      c1_c_eml_matlab_zlascl(chartInstance, c1_anrm, c1_anrmto, c1_b_A);
    }

    c1_b_eml_matlab_zggbal(chartInstance, c1_b_A, &c1_ilo, &c1_ihi, c1_rscale);
    c1_b_ilo = c1_ilo;
    c1_b_ihi = c1_ihi;
    c1_c_ilo = c1_b_ilo;
    c1_c_ihi = c1_b_ihi;
    c1_a = c1_c_ilo;
    c1_b_a = c1_a + 2;
    c1_c = c1_b_a;
    if (c1_c_ihi < c1_c) {
    } else {
      c1_c_a = c1_c_ihi;
      c1_d_a = c1_c_a - 1;
      c1_ihim1 = c1_d_a;
      c1_jcol = c1_c_ilo;
      while (c1_jcol < c1_ihim1) {
        c1_e_a = c1_jcol;
        c1_f_a = c1_e_a + 1;
        c1_jcolp1 = c1_f_a;
        c1_jrow = c1_c_ihi;
        while (c1_jrow > c1_jcolp1) {
          c1_g_a = c1_jrow;
          c1_h_a = c1_g_a - 1;
          c1_jrowm1 = c1_h_a;
          c1_d_A.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_jrowm1), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].re;
          c1_d_A.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_jrowm1), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].im;
          c1_e_A.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_jrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].re;
          c1_e_A.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_jrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].im;
          c1_eml_matlab_zlartg(chartInstance, c1_d_A, c1_e_A, &c1_b_c, &c1_s,
                               &c1_b);
          c1_c_c = c1_b_c;
          c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_jrowm1), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].re = c1_b.re;
          c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_jrowm1), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].im = c1_b.im;
          c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_jrow), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].re = c1_dc0.re;
          c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_jrow), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].im = c1_dc0.im;
          c1_d_c = c1_c_c;
          c1_xrow = c1_jrowm1;
          c1_yrow = c1_jrow;
          c1_jlo = c1_jcolp1;
          c1_jhi = c1_c_ihi;
          c1_b_jlo = c1_jlo;
          c1_b_jhi = c1_jhi;
          c1_i_a = c1_b_jlo;
          c1_b_b = c1_b_jhi;
          c1_j_a = c1_i_a;
          c1_c_b = c1_b_b;
          if (c1_j_a > c1_c_b) {
            c1_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_overflow = (c1_c_b > 2147483646);
          }

          if (c1_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_overflow);
          }

          for (c1_j = c1_b_jlo; c1_j <= c1_b_jhi; c1_j++) {
            c1_b_j = c1_j;
            c1_k_a = c1_d_c;
            c1_b.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c1_b.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c1_y.re = c1_k_a * c1_b.re;
            c1_y.im = c1_k_a * c1_b.im;
            c1_b_s.re = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re - c1_s.im * c1_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c1_b_s.im = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im + c1_s.im * c1_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c1_stemp.re = c1_y.re + c1_b_s.re;
            c1_stemp.im = c1_y.im + c1_b_s.im;
            c1_l_a = c1_d_c;
            c1_b.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c1_b.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c1_y.re = c1_l_a * c1_b.re;
            c1_y.im = c1_l_a * c1_b.im;
            c1_b.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c1_b.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c1_d_b = c1_b;
            c1_e_b = c1_b;
            c1_f_b = c1_b;
            c1_g_b = c1_b;
            c1_b.re = c1_s.re * c1_d_b.re + c1_s.im * c1_e_b.im;
            c1_b.im = c1_s.re * c1_f_b.im - c1_s.im * c1_g_b.re;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re = c1_y.re
              - c1_b.re;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im = c1_y.im
              - c1_b.im;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re =
              c1_stemp.re;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im =
              c1_stemp.im;
          }

          c1_s.re = -c1_s.re;
          c1_s.im = -c1_s.im;
          c1_e_c = c1_c_c;
          c1_xcol = c1_jrow;
          c1_ycol = c1_jrowm1;
          c1_d_ilo = c1_c_ilo;
          c1_d_ihi = c1_c_ihi;
          c1_e_ilo = c1_d_ilo;
          c1_e_ihi = c1_d_ihi;
          c1_m_a = c1_e_ilo;
          c1_h_b = c1_e_ihi;
          c1_n_a = c1_m_a;
          c1_i_b = c1_h_b;
          if (c1_n_a > c1_i_b) {
            c1_b_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_b_overflow = (c1_i_b > 2147483646);
          }

          if (c1_b_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
          }

          for (c1_i = c1_e_ilo; c1_i <= c1_e_ihi; c1_i++) {
            c1_b_i = c1_i;
            c1_o_a = c1_e_c;
            c1_b.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].re;
            c1_b.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].im;
            c1_y.re = c1_o_a * c1_b.re;
            c1_y.im = c1_o_a * c1_b.im;
            c1_c_s.re = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re - c1_s.im * c1_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im;
            c1_c_s.im = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im + c1_s.im * c1_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re;
            c1_stemp.re = c1_y.re + c1_c_s.re;
            c1_stemp.im = c1_y.im + c1_c_s.im;
            c1_p_a = c1_e_c;
            c1_b.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re;
            c1_b.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im;
            c1_y.re = c1_p_a * c1_b.re;
            c1_y.im = c1_p_a * c1_b.im;
            c1_b.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].re;
            c1_b.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].im;
            c1_j_b = c1_b;
            c1_k_b = c1_b;
            c1_l_b = c1_b;
            c1_m_b = c1_b;
            c1_b.re = c1_s.re * c1_j_b.re + c1_s.im * c1_k_b.im;
            c1_b.im = c1_s.re * c1_l_b.im - c1_s.im * c1_m_b.re;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re = c1_y.re
              - c1_b.re;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im = c1_y.im
              - c1_b.im;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].re =
              c1_stemp.re;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].im =
              c1_stemp.im;
          }

          c1_jrow = c1_jrowm1;
        }

        c1_jcol = c1_jcolp1;
      }
    }

    for (c1_i43 = 0; c1_i43 < 9; c1_i43++) {
      c1_f_A[c1_i43] = c1_b_A[c1_i43];
    }

    c1_eml_matlab_zhgeqz(chartInstance, c1_f_A, c1_b_ilo, c1_b_ihi, &c1_b_info,
                         c1_alpha1, c1_beta1);
    c1_info = c1_b_info;
    if (c1_info != 0.0) {
    } else {
      if (c1_ilascl) {
        c1_d_eml_matlab_zlascl(chartInstance, c1_anrmto, c1_anrm, c1_alpha1);
      }
    }
  }

  c1_c_info = c1_info;
  c1_d_info = c1_c_info;
  c1_e_info = c1_d_info;
  for (c1_i44 = 0; c1_i44 < 3; c1_i44++) {
    c1_b_alpha1[c1_i44] = c1_alpha1[c1_i44];
  }

  for (c1_i45 = 0; c1_i45 < 3; c1_i45++) {
    c1_b_beta1[c1_i45] = c1_beta1[c1_i45];
  }

  c1_b_eml_div(chartInstance, c1_b_alpha1, c1_b_beta1, c1_V);
  if (c1_e_info < 0.0) {
    c1_eml_warning(chartInstance);
  } else {
    if (c1_e_info > 0.0) {
      c1_b_eml_warning(chartInstance);
    }
  }
}

static void c1_realmin(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_eml_error(SFc1_Model_01InstanceStruct *chartInstance)
{
  int32_T c1_i46;
  static char_T c1_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c1_u[30];
  const mxArray *c1_y = NULL;
  int32_T c1_i47;
  static char_T c1_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c1_b_u[4];
  const mxArray *c1_b_y = NULL;
  (void)chartInstance;
  for (c1_i46 = 0; c1_i46 < 30; c1_i46++) {
    c1_u[c1_i46] = c1_cv0[c1_i46];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  for (c1_i47 = 0; c1_i47 < 4; c1_i47++) {
    c1_b_u[c1_i47] = c1_cv1[c1_i47];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c1_y, 14, c1_b_y));
}

static void c1_eps(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c1_eml_matlab_zlangeM(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_x[9])
{
  real_T c1_y;
  int32_T c1_k;
  real_T c1_b_k;
  creal_T c1_b_x;
  real_T c1_x1;
  real_T c1_x2;
  real_T c1_a;
  real_T c1_b;
  real_T c1_absxk;
  real_T c1_c_x;
  boolean_T c1_b_b;
  boolean_T exitg1;
  (void)chartInstance;
  c1_y = 0.0;
  c1_k = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c1_k < 9)) {
    c1_b_k = 1.0 + (real_T)c1_k;
    c1_b_x.re = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", c1_b_k), 1, 9, 1, 0) - 1].re;
    c1_b_x.im = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", c1_b_k), 1, 9, 1, 0) - 1].im;
    c1_x1 = c1_b_x.re;
    c1_x2 = c1_b_x.im;
    c1_a = c1_x1;
    c1_b = c1_x2;
    c1_absxk = muDoubleScalarHypot(c1_a, c1_b);
    c1_c_x = c1_absxk;
    c1_b_b = muDoubleScalarIsNaN(c1_c_x);
    if (c1_b_b) {
      c1_y = rtNaN;
      exitg1 = true;
    } else {
      if (c1_absxk > c1_y) {
        c1_y = c1_absxk;
      }

      c1_k++;
    }
  }

  return c1_y;
}

static real_T c1_abs(SFc1_Model_01InstanceStruct *chartInstance, creal_T c1_x)
{
  real_T c1_x1;
  real_T c1_x2;
  real_T c1_a;
  real_T c1_b;
  (void)chartInstance;
  c1_x1 = c1_x.re;
  c1_x2 = c1_x.im;
  c1_a = c1_x1;
  c1_b = c1_x2;
  return muDoubleScalarHypot(c1_a, c1_b);
}

static boolean_T c1_isfinite(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_x)
{
  real_T c1_b_x;
  boolean_T c1_b_b;
  boolean_T c1_b0;
  real_T c1_c_x;
  boolean_T c1_c_b;
  boolean_T c1_b1;
  (void)chartInstance;
  c1_b_x = c1_x;
  c1_b_b = muDoubleScalarIsInf(c1_b_x);
  c1_b0 = !c1_b_b;
  c1_c_x = c1_x;
  c1_c_b = muDoubleScalarIsNaN(c1_c_x);
  c1_b1 = !c1_c_b;
  return c1_b0 && c1_b1;
}

static void c1_eml_matlab_zlascl(SFc1_Model_01InstanceStruct *chartInstance,
  real_T c1_cfrom, real_T c1_cto, creal_T c1_A[9], creal_T c1_b_A[9])
{
  int32_T c1_i48;
  for (c1_i48 = 0; c1_i48 < 9; c1_i48++) {
    c1_b_A[c1_i48] = c1_A[c1_i48];
  }

  c1_c_eml_matlab_zlascl(chartInstance, c1_cfrom, c1_cto, c1_b_A);
}

static real_T c1_b_abs(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_x)
{
  real_T c1_b_x;
  (void)chartInstance;
  c1_b_x = c1_x;
  return muDoubleScalarAbs(c1_b_x);
}

static void c1_eml_matlab_zggbal(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], creal_T c1_b_A[9], int32_T *c1_ilo, int32_T *c1_ihi, int32_T
  c1_rscale[3])
{
  int32_T c1_i49;
  for (c1_i49 = 0; c1_i49 < 9; c1_i49++) {
    c1_b_A[c1_i49] = c1_A[c1_i49];
  }

  c1_b_eml_matlab_zggbal(chartInstance, c1_b_A, c1_ilo, c1_ihi, c1_rscale);
}

static void c1_eml_switch_helper(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_check_forloop_overflow_error(SFc1_Model_01InstanceStruct
  *chartInstance, boolean_T c1_overflow)
{
  int32_T c1_i50;
  static char_T c1_cv2[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c1_u[34];
  const mxArray *c1_y = NULL;
  int32_T c1_i51;
  static char_T c1_cv3[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c1_b_u[23];
  const mxArray *c1_b_y = NULL;
  (void)chartInstance;
  if (!c1_overflow) {
  } else {
    for (c1_i50 = 0; c1_i50 < 34; c1_i50++) {
      c1_u[c1_i50] = c1_cv2[c1_i50];
    }

    c1_y = NULL;
    sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  false);
    for (c1_i51 = 0; c1_i51 < 23; c1_i51++) {
      c1_b_u[c1_i51] = c1_cv3[c1_i51];
    }

    c1_b_y = NULL;
    sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 2U, 14, c1_y, 14, c1_b_y));
  }
}

static void c1_eml_matlab_zlartg(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_f, creal_T c1_g, real_T *c1_cs, creal_T *c1_sn, creal_T *c1_r)
{
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_y;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_b_y;
  real_T c1_e_x;
  real_T c1_c_y;
  real_T c1_d_y;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_e_y;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_f_y;
  real_T c1_j_x;
  real_T c1_g_y;
  real_T c1_h_y;
  real_T c1_k_x;
  real_T c1_i_y;
  real_T c1_scale;
  creal_T c1_fs;
  creal_T c1_gs;
  int32_T c1_count;
  real_T c1_rescaledir;
  int32_T c1_a;
  int32_T c1_b_a;
  static creal_T c1_dc1 = { 0.0, 0.0 };

  boolean_T c1_b_g;
  int32_T c1_c_a;
  int32_T c1_d_a;
  real_T c1_f2;
  real_T c1_g2;
  real_T c1_l_x;
  real_T c1_m_x;
  boolean_T c1_b_f;
  real_T c1_x1;
  real_T c1_x2;
  real_T c1_e_a;
  real_T c1_b;
  real_T c1_j_y;
  real_T c1_b_x1;
  real_T c1_b_x2;
  real_T c1_f_a;
  real_T c1_b_b;
  real_T c1_d;
  real_T c1_c_x1;
  real_T c1_c_x2;
  real_T c1_g_a;
  real_T c1_c_b;
  real_T c1_f2s;
  real_T c1_n_x;
  real_T c1_g2s;
  real_T c1_o_x;
  real_T c1_p_x;
  real_T c1_k_y;
  real_T c1_q_x;
  real_T c1_r_x;
  real_T c1_l_y;
  real_T c1_s_x;
  real_T c1_m_y;
  real_T c1_n_y;
  real_T c1_d_x1;
  real_T c1_d_x2;
  real_T c1_h_a;
  real_T c1_d_b;
  real_T c1_dr;
  real_T c1_di;
  real_T c1_e_x1;
  real_T c1_e_x2;
  real_T c1_i_a;
  real_T c1_e_b;
  creal_T c1_b_gs;
  real_T c1_j_a;
  creal_T c1_b_sn;
  real_T c1_t_x;
  creal_T c1_c_gs;
  creal_T c1_c_sn;
  int32_T c1_b_count;
  int32_T c1_f_b;
  int32_T c1_g_b;
  boolean_T c1_overflow;
  int32_T c1_i;
  int32_T c1_c_count;
  int32_T c1_h_b;
  int32_T c1_i_b;
  boolean_T c1_b_overflow;
  int32_T c1_b_i;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  c1_realmin(chartInstance);
  c1_eps(chartInstance);
  c1_eps(chartInstance);
  c1_x = c1_f.re;
  c1_b_x = c1_x;
  c1_y = muDoubleScalarAbs(c1_b_x);
  c1_c_x = c1_f.im;
  c1_d_x = c1_c_x;
  c1_b_y = muDoubleScalarAbs(c1_d_x);
  c1_e_x = c1_y;
  c1_c_y = c1_b_y;
  c1_d_y = c1_e_x;
  if (c1_c_y > c1_d_y) {
    c1_d_y = c1_c_y;
  }

  c1_f_x = c1_g.re;
  c1_g_x = c1_f_x;
  c1_e_y = muDoubleScalarAbs(c1_g_x);
  c1_h_x = c1_g.im;
  c1_i_x = c1_h_x;
  c1_f_y = muDoubleScalarAbs(c1_i_x);
  c1_j_x = c1_e_y;
  c1_g_y = c1_f_y;
  c1_h_y = c1_j_x;
  if (c1_g_y > c1_h_y) {
    c1_h_y = c1_g_y;
  }

  c1_k_x = c1_d_y;
  c1_i_y = c1_h_y;
  c1_scale = c1_k_x;
  if (c1_i_y > c1_scale) {
    c1_scale = c1_i_y;
  }

  c1_fs = c1_f;
  c1_gs = c1_g;
  c1_count = 0;
  c1_rescaledir = 0.0;
  guard1 = false;
  guard2 = false;
  if (c1_scale >= 7.4428285367870146E+137) {
    do {
      c1_a = c1_count;
      c1_b_a = c1_a + 1;
      c1_count = c1_b_a;
      c1_fs.re *= 1.3435752215134178E-138;
      c1_fs.im *= 1.3435752215134178E-138;
      c1_gs.re *= 1.3435752215134178E-138;
      c1_gs.im *= 1.3435752215134178E-138;
      c1_scale *= 1.3435752215134178E-138;
    } while (!(c1_scale < 7.4428285367870146E+137));

    c1_rescaledir = 1.0;
    guard1 = true;
  } else if (c1_scale <= 1.3435752215134178E-138) {
    c1_b_g = ((c1_g.re == c1_dc1.re) && (c1_g.im == c1_dc1.im));
    if (c1_b_g) {
      *c1_cs = 1.0;
      *c1_sn = c1_dc1;
      *c1_r = c1_f;
    } else {
      do {
        c1_c_a = c1_count;
        c1_d_a = c1_c_a + 1;
        c1_count = c1_d_a;
        c1_fs.re *= 7.4428285367870146E+137;
        c1_fs.im *= 7.4428285367870146E+137;
        c1_gs.re *= 7.4428285367870146E+137;
        c1_gs.im *= 7.4428285367870146E+137;
        c1_scale *= 7.4428285367870146E+137;
      } while (!(c1_scale > 1.3435752215134178E-138));

      c1_rescaledir = -1.0;
      guard2 = true;
    }
  } else {
    guard2 = true;
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c1_f2 = c1_fs.re * c1_fs.re + c1_fs.im * c1_fs.im;
    c1_g2 = c1_gs.re * c1_gs.re + c1_gs.im * c1_gs.im;
    c1_l_x = c1_g2;
    c1_m_x = c1_l_x;
    if (1.0 > c1_m_x) {
      c1_m_x = 1.0;
    }

    if (c1_f2 <= c1_m_x * 2.0041683600089728E-292) {
      c1_b_f = ((c1_f.re == c1_dc1.re) && (c1_f.im == c1_dc1.im));
      if (c1_b_f) {
        *c1_cs = 0.0;
        c1_x1 = c1_g.re;
        c1_x2 = c1_g.im;
        c1_e_a = c1_x1;
        c1_b = c1_x2;
        c1_j_y = muDoubleScalarHypot(c1_e_a, c1_b);
        c1_r->re = c1_j_y;
        c1_r->im = 0.0;
        c1_b_x1 = c1_gs.re;
        c1_b_x2 = c1_gs.im;
        c1_f_a = c1_b_x1;
        c1_b_b = c1_b_x2;
        c1_d = muDoubleScalarHypot(c1_f_a, c1_b_b);
        c1_sn->re = c1_gs.re / c1_d;
        c1_sn->im = -c1_gs.im / c1_d;
      } else {
        c1_c_x1 = c1_fs.re;
        c1_c_x2 = c1_fs.im;
        c1_g_a = c1_c_x1;
        c1_c_b = c1_c_x2;
        c1_f2s = muDoubleScalarHypot(c1_g_a, c1_c_b);
        c1_n_x = c1_g2;
        c1_g2s = c1_n_x;
        if (c1_g2s < 0.0) {
          c1_eml_error(chartInstance);
        }

        c1_g2s = muDoubleScalarSqrt(c1_g2s);
        *c1_cs = c1_f2s / c1_g2s;
        c1_o_x = c1_f.re;
        c1_p_x = c1_o_x;
        c1_k_y = muDoubleScalarAbs(c1_p_x);
        c1_q_x = c1_f.im;
        c1_r_x = c1_q_x;
        c1_l_y = muDoubleScalarAbs(c1_r_x);
        c1_s_x = c1_k_y;
        c1_m_y = c1_l_y;
        c1_n_y = c1_s_x;
        if (c1_m_y > c1_n_y) {
          c1_n_y = c1_m_y;
        }

        if (c1_n_y > 1.0) {
          c1_d_x1 = c1_f.re;
          c1_d_x2 = c1_f.im;
          c1_h_a = c1_d_x1;
          c1_d_b = c1_d_x2;
          c1_d = muDoubleScalarHypot(c1_h_a, c1_d_b);
          c1_fs.re = c1_f.re / c1_d;
          c1_fs.im = c1_f.im / c1_d;
        } else {
          c1_dr = 7.4428285367870146E+137 * c1_f.re;
          c1_di = 7.4428285367870146E+137 * c1_f.im;
          c1_e_x1 = c1_dr;
          c1_e_x2 = c1_di;
          c1_i_a = c1_e_x1;
          c1_e_b = c1_e_x2;
          c1_d = muDoubleScalarHypot(c1_i_a, c1_e_b);
          c1_fs.re = c1_dr / c1_d;
          c1_fs.im = c1_di / c1_d;
        }

        c1_b_gs.re = c1_gs.re / c1_g2s;
        c1_b_gs.im = -c1_gs.im / c1_g2s;
        c1_sn->re = c1_fs.re * c1_b_gs.re - c1_fs.im * c1_b_gs.im;
        c1_sn->im = c1_fs.re * c1_b_gs.im + c1_fs.im * c1_b_gs.re;
        c1_j_a = *c1_cs;
        c1_fs.re = c1_j_a * c1_f.re;
        c1_fs.im = c1_j_a * c1_f.im;
        c1_b_sn.re = c1_sn->re * c1_g.re - c1_sn->im * c1_g.im;
        c1_b_sn.im = c1_sn->re * c1_g.im + c1_sn->im * c1_g.re;
        c1_r->re = c1_fs.re + c1_b_sn.re;
        c1_r->im = c1_fs.im + c1_b_sn.im;
      }
    } else {
      c1_t_x = 1.0 + c1_g2 / c1_f2;
      c1_f2s = c1_t_x;
      if (c1_f2s < 0.0) {
        c1_eml_error(chartInstance);
      }

      c1_f2s = muDoubleScalarSqrt(c1_f2s);
      c1_r->re = c1_f2s * c1_fs.re;
      c1_r->im = c1_f2s * c1_fs.im;
      *c1_cs = 1.0 / c1_f2s;
      c1_d = c1_f2 + c1_g2;
      c1_sn->re = c1_r->re / c1_d;
      c1_sn->im = c1_r->im / c1_d;
      c1_c_gs.re = c1_gs.re;
      c1_c_gs.im = -c1_gs.im;
      c1_c_sn = *c1_sn;
      c1_sn->re = c1_c_sn.re * c1_c_gs.re - c1_c_sn.im * c1_c_gs.im;
      c1_sn->im = c1_c_sn.re * c1_c_gs.im + c1_c_sn.im * c1_c_gs.re;
      if (c1_rescaledir > 0.0) {
        c1_b_count = c1_count;
        c1_f_b = c1_b_count;
        c1_g_b = c1_f_b;
        if (1 > c1_g_b) {
          c1_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_overflow = (c1_g_b > 2147483646);
        }

        if (c1_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_overflow);
        }

        for (c1_i = 1; c1_i <= c1_b_count; c1_i++) {
          c1_r->re *= 7.4428285367870146E+137;
          c1_r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (c1_rescaledir < 0.0) {
          c1_c_count = c1_count;
          c1_h_b = c1_c_count;
          c1_i_b = c1_h_b;
          if (1 > c1_i_b) {
            c1_b_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_b_overflow = (c1_i_b > 2147483646);
          }

          if (c1_b_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
          }

          for (c1_b_i = 1; c1_b_i <= c1_c_count; c1_b_i++) {
            c1_r->re *= 1.3435752215134178E-138;
            c1_r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

static void c1_eml_scalar_eg(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_eml_matlab_zhgeqz(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T c1_ilo, int32_T c1_ihi, real_T *c1_info, creal_T
  c1_alpha1[3], creal_T c1_beta1[3])
{
  static creal_T c1_dc2 = { 0.0, 0.0 };

  int32_T c1_i52;
  creal_T c1_b_A[9];
  int32_T c1_i53;
  int32_T c1_i54;
  static creal_T c1_dc3 = { 0.0, 0.0 };

  creal_T c1_eshift;
  creal_T c1_ctemp;
  creal_T c1_rho;
  int32_T c1_i55;
  creal_T c1_c_A[9];
  real_T c1_anorm;
  real_T c1_y;
  real_T c1_atol;
  real_T c1_b_y;
  real_T c1_x;
  real_T c1_ascale;
  boolean_T c1_failed;
  int32_T c1_a;
  int32_T c1_b_a;
  int32_T c1_i56;
  int32_T c1_c_a;
  int32_T c1_d_a;
  boolean_T c1_overflow;
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_ifirst;
  int32_T c1_istart;
  int32_T c1_ilast;
  int32_T c1_e_a;
  int32_T c1_f_a;
  int32_T c1_ilastm1;
  int32_T c1_ifrstm;
  int32_T c1_ilastm;
  int32_T c1_iiter;
  int32_T c1_g_a;
  int32_T c1_b;
  int32_T c1_h_a;
  int32_T c1_b_b;
  int32_T c1_c;
  int32_T c1_i_a;
  int32_T c1_j_a;
  int32_T c1_b_c;
  int32_T c1_c_b;
  int32_T c1_d_b;
  int32_T c1_maxit;
  boolean_T c1_goto50;
  boolean_T c1_goto60;
  boolean_T c1_goto70;
  boolean_T c1_goto90;
  int32_T c1_b_maxit;
  int32_T c1_e_b;
  int32_T c1_f_b;
  boolean_T c1_b_overflow;
  int32_T c1_jiter;
  creal_T c1_a22;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_c_y;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_d_y;
  real_T c1_e_y;
  int32_T c1_k_a;
  int32_T c1_l_a;
  int32_T c1_jm1;
  boolean_T c1_ilazro;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_f_y;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_g_y;
  real_T c1_h_y;
  boolean_T c1_b2;
  int32_T c1_i57;
  int32_T c1_i58;
  creal_T c1_d_A;
  creal_T c1_e_A;
  creal_T c1_s;
  real_T c1_c_c;
  real_T c1_d_c;
  real_T c1_e_c;
  int32_T c1_xcol;
  int32_T c1_ycol;
  int32_T c1_b_ilo;
  int32_T c1_b_ihi;
  int32_T c1_c_ilo;
  int32_T c1_c_ihi;
  int32_T c1_m_a;
  int32_T c1_g_b;
  int32_T c1_n_a;
  int32_T c1_h_b;
  boolean_T c1_c_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  real_T c1_o_a;
  creal_T c1_a12;
  creal_T c1_b_s;
  creal_T c1_a21;
  real_T c1_p_a;
  creal_T c1_b_a22;
  creal_T c1_c_a22;
  creal_T c1_d_a22;
  creal_T c1_e_a22;
  int32_T c1_q_a;
  int32_T c1_r_a;
  int32_T c1_s_a;
  int32_T c1_t_a;
  creal_T c1_r2;
  creal_T c1_f_a22;
  creal_T c1_b_rho;
  creal_T c1_b_a12;
  creal_T c1_c_a12;
  creal_T c1_b_a21;
  real_T c1_d4;
  real_T c1_d5;
  int32_T c1_u_a;
  int32_T c1_v_a;
  int32_T c1_jp1;
  int32_T c1_w_a;
  int32_T c1_x_a;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_i_y;
  real_T c1_l_x;
  real_T c1_m_x;
  real_T c1_j_y;
  real_T c1_k_y;
  real_T c1_temp;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_l_y;
  real_T c1_p_x;
  real_T c1_q_x;
  real_T c1_m_y;
  real_T c1_n_y;
  real_T c1_temp2;
  real_T c1_r_x;
  real_T c1_o_y;
  real_T c1_tempr;
  real_T c1_s_x;
  real_T c1_t_x;
  real_T c1_p_y;
  real_T c1_u_x;
  real_T c1_v_x;
  real_T c1_q_y;
  real_T c1_r_y;
  int32_T c1_y_a;
  int32_T c1_ab_a;
  int32_T c1_f_c;
  real_T c1_g_c;
  int32_T c1_bb_a;
  int32_T c1_cb_a;
  int32_T c1_db_a;
  int32_T c1_eb_a;
  creal_T c1_f_A;
  creal_T c1_g_A;
  real_T c1_h_c;
  real_T c1_i_c;
  int32_T c1_xrow;
  int32_T c1_yrow;
  int32_T c1_jlo;
  int32_T c1_jhi;
  int32_T c1_b_jlo;
  int32_T c1_b_jhi;
  int32_T c1_fb_a;
  int32_T c1_i_b;
  int32_T c1_gb_a;
  int32_T c1_j_b;
  boolean_T c1_d_overflow;
  int32_T c1_c_j;
  int32_T c1_d_j;
  real_T c1_hb_a;
  creal_T c1_c_s;
  real_T c1_ib_a;
  creal_T c1_g_a22;
  creal_T c1_h_a22;
  creal_T c1_i_a22;
  creal_T c1_j_a22;
  int32_T c1_jb_a;
  int32_T c1_kb_a;
  int32_T c1_j_c;
  int32_T c1_w_x;
  int32_T c1_s_y;
  int32_T c1_x_x;
  real_T c1_k_c;
  int32_T c1_b_xcol;
  int32_T c1_b_ycol;
  int32_T c1_d_ilo;
  int32_T c1_d_ihi;
  int32_T c1_e_ilo;
  int32_T c1_e_ihi;
  int32_T c1_lb_a;
  int32_T c1_k_b;
  int32_T c1_mb_a;
  int32_T c1_l_b;
  boolean_T c1_e_overflow;
  int32_T c1_c_i;
  int32_T c1_d_i;
  real_T c1_nb_a;
  creal_T c1_d_s;
  real_T c1_ob_a;
  creal_T c1_k_a22;
  creal_T c1_l_a22;
  creal_T c1_m_a22;
  creal_T c1_n_a22;
  int32_T c1_b_ilast;
  int32_T c1_m_b;
  int32_T c1_n_b;
  boolean_T c1_f_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  int32_T c1_pb_a;
  int32_T c1_qb_a;
  int32_T c1_i59;
  int32_T c1_o_b;
  int32_T c1_p_b;
  boolean_T c1_g_overflow;
  int32_T c1_e_j;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  int32_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T guard11 = false;
  c1_dc2.re = rtNaN;
  for (c1_i52 = 0; c1_i52 < 9; c1_i52++) {
    c1_b_A[c1_i52] = c1_A[c1_i52];
  }

  for (c1_i53 = 0; c1_i53 < 3; c1_i53++) {
    c1_alpha1[c1_i53].re = 0.0;
    c1_alpha1[c1_i53].im = 0.0;
  }

  for (c1_i54 = 0; c1_i54 < 3; c1_i54++) {
    c1_beta1[c1_i54].re = 1.0;
    c1_beta1[c1_i54].im = 0.0;
  }

  c1_eps(chartInstance);
  c1_realmin(chartInstance);
  c1_eshift = c1_dc3;
  c1_ctemp = c1_dc3;
  c1_rho = c1_dc3;
  for (c1_i55 = 0; c1_i55 < 9; c1_i55++) {
    c1_c_A[c1_i55] = c1_b_A[c1_i55];
  }

  c1_anorm = c1_eml_matlab_zlanhs(chartInstance, c1_c_A, c1_ilo, c1_ihi);
  c1_y = 2.2204460492503131E-16 * c1_anorm;
  c1_atol = 2.2250738585072014E-308;
  if (c1_y > 2.2250738585072014E-308) {
    c1_atol = c1_y;
  }

  c1_b_y = c1_anorm;
  c1_x = 2.2250738585072014E-308;
  if (c1_b_y > 2.2250738585072014E-308) {
    c1_x = c1_b_y;
  }

  c1_ascale = 1.0 / c1_x;
  c1_failed = true;
  c1_a = c1_ihi;
  c1_b_a = c1_a + 1;
  c1_i56 = c1_b_a;
  c1_c_a = c1_i56;
  c1_d_a = c1_c_a;
  if (c1_d_a > 3) {
    c1_overflow = false;
  } else {
    c1_eml_switch_helper(chartInstance);
    c1_overflow = false;
  }

  if (c1_overflow) {
    c1_check_forloop_overflow_error(chartInstance, c1_overflow);
  }

  for (c1_j = c1_i56; c1_j < 4; c1_j++) {
    c1_b_j = c1_j;
    c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_b_j), 1, 3, 1, 0) - 1].re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
      ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
    c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_b_j), 1, 3, 1, 0) - 1].im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
      ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
  }

  guard1 = false;
  guard2 = false;
  if (c1_ihi >= c1_ilo) {
    c1_ifirst = c1_ilo;
    c1_istart = c1_ilo;
    c1_ilast = c1_ihi;
    c1_e_a = c1_ilast;
    c1_f_a = c1_e_a - 1;
    c1_ilastm1 = c1_f_a;
    c1_ifrstm = c1_ilo;
    c1_ilastm = c1_ihi;
    c1_iiter = 0;
    c1_g_a = c1_ihi;
    c1_b = c1_ilo;
    c1_h_a = c1_g_a;
    c1_b_b = c1_b;
    c1_c = c1_h_a - c1_b_b;
    c1_i_a = c1_c;
    c1_j_a = c1_i_a;
    c1_b_c = c1_j_a;
    c1_c_b = c1_b_c + 1;
    c1_d_b = c1_c_b;
    c1_maxit = 30 * c1_d_b;
    c1_goto50 = false;
    c1_goto60 = false;
    c1_goto70 = false;
    c1_goto90 = false;
    c1_b_maxit = c1_maxit;
    c1_e_b = c1_b_maxit;
    c1_f_b = c1_e_b;
    if (1 > c1_f_b) {
      c1_b_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_b_overflow = (c1_f_b > 2147483646);
    }

    if (c1_b_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
    }

    c1_jiter = 1;
    do {
      exitg1 = 0;
      if (c1_jiter <= c1_b_maxit) {
        if (c1_ilast == c1_ilo) {
          c1_goto60 = true;
        } else {
          c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            re;
          c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            im;
          c1_b_x = c1_a22.re;
          c1_c_x = c1_b_x;
          c1_c_y = muDoubleScalarAbs(c1_c_x);
          c1_d_x = c1_a22.im;
          c1_e_x = c1_d_x;
          c1_d_y = muDoubleScalarAbs(c1_e_x);
          c1_e_y = c1_c_y + c1_d_y;
          if (c1_e_y <= c1_atol) {
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c1_dc3.re;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c1_dc3.im;
            c1_goto60 = true;
          } else {
            c1_b_j = c1_ilastm1;
            exitg3 = false;
            while ((exitg3 == false) && (c1_b_j >= c1_ilo)) {
              c1_k_a = c1_b_j;
              c1_l_a = c1_k_a - 1;
              c1_jm1 = c1_l_a;
              if (c1_b_j == c1_ilo) {
                c1_ilazro = true;
              } else {
                c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c1_f_x = c1_a22.re;
                c1_g_x = c1_f_x;
                c1_f_y = muDoubleScalarAbs(c1_g_x);
                c1_h_x = c1_a22.im;
                c1_i_x = c1_h_x;
                c1_g_y = muDoubleScalarAbs(c1_i_x);
                c1_h_y = c1_f_y + c1_g_y;
                if (c1_h_y <= c1_atol) {
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0)
                               - 1)) - 1].re = c1_dc3.re;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0)
                               - 1)) - 1].im = c1_dc3.im;
                  c1_ilazro = true;
                } else {
                  c1_ilazro = false;
                }
              }

              if (c1_ilazro) {
                c1_ifirst = c1_b_j;
                c1_goto70 = true;
                exitg3 = true;
              } else {
                c1_b_j = c1_jm1;
              }
            }
          }
        }

        guard3 = false;
        guard4 = false;
        if (c1_goto50) {
          guard4 = true;
        } else if (c1_goto60) {
          guard4 = true;
        } else if (c1_goto70) {
          guard3 = true;
        } else {
          c1_b2 = false;
        }

        if (guard4 == true) {
          guard3 = true;
        }

        if (guard3 == true) {
          c1_b2 = true;
        }

        if (!c1_b2) {
          for (c1_i57 = 0; c1_i57 < 3; c1_i57++) {
            c1_alpha1[c1_i57] = c1_dc2;
          }

          for (c1_i58 = 0; c1_i58 < 3; c1_i58++) {
            c1_beta1[c1_i58] = c1_dc2;
          }

          *c1_info = -1.0;
          exitg1 = 1;
        } else {
          if (c1_goto50) {
            c1_goto50 = false;
            c1_d_A.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].
              re;
            c1_d_A.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].
              im;
            c1_e_A.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .re;
            c1_e_A.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .im;
            c1_eml_matlab_zlartg(chartInstance, c1_d_A, c1_e_A, &c1_c_c, &c1_s,
                                 &c1_a22);
            c1_d_c = c1_c_c;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].re =
              c1_a22.re;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].im =
              c1_a22.im;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c1_dc3.re;
            c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c1_dc3.im;
            c1_e_c = c1_d_c;
            c1_xcol = c1_ilast;
            c1_ycol = c1_ilastm1;
            c1_b_ilo = c1_ifrstm;
            c1_b_ihi = c1_ilastm1;
            c1_c_ilo = c1_b_ilo;
            c1_c_ihi = c1_b_ihi;
            c1_m_a = c1_c_ilo;
            c1_g_b = c1_c_ihi;
            c1_n_a = c1_m_a;
            c1_h_b = c1_g_b;
            if (c1_n_a > c1_h_b) {
              c1_c_overflow = false;
            } else {
              c1_eml_switch_helper(chartInstance);
              c1_c_overflow = (c1_h_b > 2147483646);
            }

            if (c1_c_overflow) {
              c1_check_forloop_overflow_error(chartInstance, c1_c_overflow);
            }

            for (c1_i = c1_c_ilo; c1_i <= c1_c_ihi; c1_i++) {
              c1_b_i = c1_i;
              c1_o_a = c1_e_c;
              c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c1_a12.re = c1_o_a * c1_a22.re;
              c1_a12.im = c1_o_a * c1_a22.im;
              c1_b_s.re = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re - c1_s.im *
                c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) -
                           1)) - 1].im;
              c1_b_s.im = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im + c1_s.im *
                c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) -
                           1)) - 1].re;
              c1_a21.re = c1_a12.re + c1_b_s.re;
              c1_a21.im = c1_a12.im + c1_b_s.im;
              c1_p_a = c1_e_c;
              c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c1_a12.re = c1_p_a * c1_a22.re;
              c1_a12.im = c1_p_a * c1_a22.im;
              c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c1_b_a22 = c1_a22;
              c1_c_a22 = c1_a22;
              c1_d_a22 = c1_a22;
              c1_e_a22 = c1_a22;
              c1_a22.re = c1_s.re * c1_b_a22.re + c1_s.im * c1_c_a22.im;
              c1_a22.im = c1_s.re * c1_d_a22.im - c1_s.im * c1_e_a22.re;
              c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1))
                - 1].re = c1_a12.re - c1_a22.re;
              c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1))
                - 1].im = c1_a12.im - c1_a22.im;
              c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1))
                - 1].re = c1_a21.re;
              c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1))
                - 1].im = c1_a21.im;
            }

            c1_goto60 = true;
          }

          guard11 = false;
          if (c1_goto60) {
            c1_goto60 = false;
            c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) - 1].re =
              c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3
                      * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) -
                         1)) - 1].re;
            c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) - 1].im =
              c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3
                      * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) -
                         1)) - 1].im;
            c1_ilast = c1_ilastm1;
            c1_q_a = c1_ilast;
            c1_r_a = c1_q_a - 1;
            c1_ilastm1 = c1_r_a;
            if (c1_ilast < c1_ilo) {
              c1_failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              c1_iiter = 0;
              c1_eshift = c1_dc3;
              c1_ilastm = c1_ilast;
              if (c1_ifrstm > c1_ilast) {
                c1_ifrstm = c1_ilo;
              }

              guard11 = true;
            }
          } else {
            if (c1_goto70) {
              c1_goto70 = false;
              c1_s_a = c1_iiter;
              c1_t_a = c1_s_a + 1;
              c1_iiter = c1_t_a;
              c1_ifrstm = c1_ifirst;
              if (c1_mod(chartInstance, c1_iiter) != 0) {
                c1_s.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c1_s.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c1_r2.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                   (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c1_r2.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                   (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c1_a12.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c1_a12.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c1_a21.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c1_a21.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c1_a22.re = c1_r2.re - c1_s.re;
                c1_a22.im = c1_r2.im - c1_s.im;
                c1_f_a22.re = -c1_a22.re;
                c1_f_a22.im = -c1_a22.im;
                c1_rho = c1_eml_div(chartInstance, c1_f_a22, 2.0);
                c1_b_rho.re = c1_rho.re * c1_rho.re - c1_rho.im * c1_rho.im;
                c1_b_rho.im = c1_rho.re * c1_rho.im + c1_rho.im * c1_rho.re;
                c1_b_a12.re = c1_a12.re * c1_a21.re - c1_a12.im * c1_a21.im;
                c1_b_a12.im = c1_a12.re * c1_a21.im + c1_a12.im * c1_a21.re;
                c1_a22.re = c1_b_rho.re + c1_b_a12.re;
                c1_a22.im = c1_b_rho.im + c1_b_a12.im;
                c1_b_sqrt(chartInstance, &c1_a22);
                c1_a12.re = c1_s.re - (c1_rho.re - c1_a22.re);
                c1_a12.im = c1_s.im - (c1_rho.im - c1_a22.im);
                c1_a21.re = c1_s.re - (c1_rho.re + c1_a22.re);
                c1_a21.im = c1_s.im - (c1_rho.im + c1_a22.im);
                c1_c_a12.re = c1_a12.re - c1_r2.re;
                c1_c_a12.im = c1_a12.im - c1_r2.im;
                c1_b_a21.re = c1_a21.re - c1_r2.re;
                c1_b_a21.im = c1_a21.im - c1_r2.im;
                c1_d4 = c1_abs(chartInstance, c1_c_a12);
                c1_d5 = c1_abs(chartInstance, c1_b_a21);
                if (c1_d4 <= c1_d5) {
                  c1_a21 = c1_a12;
                  c1_rho.re -= c1_a22.re;
                  c1_rho.im -= c1_a22.im;
                } else {
                  c1_rho.re += c1_a22.re;
                  c1_rho.im += c1_a22.im;
                }
              } else {
                c1_eshift.re += c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                  "", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].re;
                c1_eshift.im += c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                  "", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].im;
                c1_a21 = c1_eshift;
              }

              c1_b_j = c1_ilastm1;
              c1_u_a = c1_b_j;
              c1_v_a = c1_u_a + 1;
              c1_jp1 = c1_v_a;
              exitg2 = false;
              while ((exitg2 == false) && (c1_b_j > c1_ifirst)) {
                c1_w_a = c1_b_j;
                c1_x_a = c1_w_a - 1;
                c1_jm1 = c1_x_a;
                c1_istart = c1_b_j;
                c1_ctemp.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re - c1_a21.re;
                c1_ctemp.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im - c1_a21.im;
                c1_j_x = c1_ctemp.re;
                c1_k_x = c1_j_x;
                c1_i_y = muDoubleScalarAbs(c1_k_x);
                c1_l_x = c1_ctemp.im;
                c1_m_x = c1_l_x;
                c1_j_y = muDoubleScalarAbs(c1_m_x);
                c1_k_y = c1_i_y + c1_j_y;
                c1_temp = c1_ascale * c1_k_y;
                c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c1_n_x = c1_a22.re;
                c1_o_x = c1_n_x;
                c1_l_y = muDoubleScalarAbs(c1_o_x);
                c1_p_x = c1_a22.im;
                c1_q_x = c1_p_x;
                c1_m_y = muDoubleScalarAbs(c1_q_x);
                c1_n_y = c1_l_y + c1_m_y;
                c1_temp2 = c1_ascale * c1_n_y;
                c1_r_x = c1_temp;
                c1_o_y = c1_temp2;
                c1_tempr = c1_r_x;
                if (c1_o_y > c1_tempr) {
                  c1_tempr = c1_o_y;
                }

                if (c1_tempr < 1.0) {
                  if (c1_tempr != 0.0) {
                    c1_temp /= c1_tempr;
                    c1_temp2 /= c1_tempr;
                  }
                }

                c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c1_s_x = c1_a22.re;
                c1_t_x = c1_s_x;
                c1_p_y = muDoubleScalarAbs(c1_t_x);
                c1_u_x = c1_a22.im;
                c1_v_x = c1_u_x;
                c1_q_y = muDoubleScalarAbs(c1_v_x);
                c1_r_y = c1_p_y + c1_q_y;
                if (c1_r_y * c1_temp2 <= c1_temp * c1_atol) {
                  c1_goto90 = true;
                  exitg2 = true;
                } else {
                  c1_jp1 = c1_b_j;
                  c1_b_j = c1_jm1;
                }
              }

              if (!c1_goto90) {
                c1_istart = c1_ifirst;
                if (c1_istart == c1_ilastm1) {
                  c1_ctemp = c1_rho;
                } else {
                  c1_ctemp.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 1, 0) + 3 *
                                        (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 2,
                    0) - 1)) - 1].re - c1_a21.re;
                  c1_ctemp.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 1, 0) + 3 *
                                        (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 2,
                    0) - 1)) - 1].im - c1_a21.im;
                }

                c1_goto90 = true;
              }
            }

            if (c1_goto90) {
              c1_goto90 = false;
              c1_y_a = c1_istart;
              c1_ab_a = c1_y_a;
              c1_f_c = c1_ab_a;
              c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c1_f_c + 1)), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c1_f_c + 1)), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c1_b_eml_matlab_zlartg(chartInstance, c1_ctemp, c1_a22, &c1_g_c,
                &c1_s);
              c1_d_c = c1_g_c;
              c1_b_j = c1_istart;
              c1_bb_a = c1_b_j;
              c1_cb_a = c1_bb_a - 1;
              c1_jm1 = c1_cb_a;
              while (c1_b_j < c1_ilast) {
                c1_db_a = c1_b_j;
                c1_eb_a = c1_db_a + 1;
                c1_jp1 = c1_eb_a;
                if (c1_b_j > c1_istart) {
                  c1_f_A.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_f_A.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_g_A.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_g_A.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_eml_matlab_zlartg(chartInstance, c1_f_A, c1_g_A, &c1_h_c,
                                       &c1_s, &c1_a22);
                  c1_d_c = c1_h_c;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0)
                               - 1)) - 1].re = c1_a22.re;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0)
                               - 1)) - 1].im = c1_a22.im;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0)
                               - 1)) - 1].re = c1_dc3.re;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0)
                               - 1)) - 1].im = c1_dc3.im;
                }

                c1_i_c = c1_d_c;
                c1_xrow = c1_b_j;
                c1_yrow = c1_jp1;
                c1_jlo = c1_b_j;
                c1_jhi = c1_ilastm;
                c1_b_jlo = c1_jlo;
                c1_b_jhi = c1_jhi;
                c1_fb_a = c1_b_jlo;
                c1_i_b = c1_b_jhi;
                c1_gb_a = c1_fb_a;
                c1_j_b = c1_i_b;
                if (c1_gb_a > c1_j_b) {
                  c1_d_overflow = false;
                } else {
                  c1_eml_switch_helper(chartInstance);
                  c1_d_overflow = (c1_j_b > 2147483646);
                }

                if (c1_d_overflow) {
                  c1_check_forloop_overflow_error(chartInstance, c1_d_overflow);
                }

                for (c1_c_j = c1_b_jlo; c1_c_j <= c1_b_jhi; c1_c_j++) {
                  c1_d_j = c1_c_j;
                  c1_hb_a = c1_i_c;
                  c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_a12.re = c1_hb_a * c1_a22.re;
                  c1_a12.im = c1_hb_a * c1_a22.im;
                  c1_c_s.re = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re - c1_s.im * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_c_s.im = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im + c1_s.im * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_a21.re = c1_a12.re + c1_c_s.re;
                  c1_a21.im = c1_a12.im + c1_c_s.im;
                  c1_ib_a = c1_i_c;
                  c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_a12.re = c1_ib_a * c1_a22.re;
                  c1_a12.im = c1_ib_a * c1_a22.im;
                  c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_g_a22 = c1_a22;
                  c1_h_a22 = c1_a22;
                  c1_i_a22 = c1_a22;
                  c1_j_a22 = c1_a22;
                  c1_a22.re = c1_s.re * c1_g_a22.re + c1_s.im * c1_h_a22.im;
                  c1_a22.im = c1_s.re * c1_i_a22.im - c1_s.im * c1_j_a22.re;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].re = c1_a12.re - c1_a22.re;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].im = c1_a12.im - c1_a22.im;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].re = c1_a21.re;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].im = c1_a21.im;
                }

                c1_s.re = -c1_s.re;
                c1_s.im = -c1_s.im;
                c1_jb_a = c1_jp1;
                c1_kb_a = c1_jb_a;
                c1_j_c = c1_kb_a;
                c1_w_x = c1_j_c + 1;
                c1_s_y = c1_ilast;
                c1_x_x = c1_w_x;
                if (c1_s_y < c1_x_x) {
                  c1_x_x = c1_s_y;
                }

                c1_k_c = c1_d_c;
                c1_b_xcol = c1_jp1;
                c1_b_ycol = c1_b_j;
                c1_d_ilo = c1_ifrstm;
                c1_d_ihi = c1_x_x;
                c1_e_ilo = c1_d_ilo;
                c1_e_ihi = c1_d_ihi;
                c1_lb_a = c1_e_ilo;
                c1_k_b = c1_e_ihi;
                c1_mb_a = c1_lb_a;
                c1_l_b = c1_k_b;
                if (c1_mb_a > c1_l_b) {
                  c1_e_overflow = false;
                } else {
                  c1_eml_switch_helper(chartInstance);
                  c1_e_overflow = (c1_l_b > 2147483646);
                }

                if (c1_e_overflow) {
                  c1_check_forloop_overflow_error(chartInstance, c1_e_overflow);
                }

                for (c1_c_i = c1_e_ilo; c1_c_i <= c1_e_ihi; c1_c_i++) {
                  c1_d_i = c1_c_i;
                  c1_nb_a = c1_k_c;
                  c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_a12.re = c1_nb_a * c1_a22.re;
                  c1_a12.im = c1_nb_a * c1_a22.im;
                  c1_d_s.re = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].re - c1_s.im * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_d_s.im = c1_s.re * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].im + c1_s.im * c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a21.re = c1_a12.re + c1_d_s.re;
                  c1_a21.im = c1_a12.im + c1_d_s.im;
                  c1_ob_a = c1_k_c;
                  c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_a12.re = c1_ob_a * c1_a22.re;
                  c1_a12.im = c1_ob_a * c1_a22.im;
                  c1_a22.re = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a22.im = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_k_a22 = c1_a22;
                  c1_l_a22 = c1_a22;
                  c1_m_a22 = c1_a22;
                  c1_n_a22 = c1_a22;
                  c1_a22.re = c1_s.re * c1_k_a22.re + c1_s.im * c1_l_a22.im;
                  c1_a22.im = c1_s.re * c1_m_a22.im - c1_s.im * c1_n_a22.re;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2,
                            0) - 1)) - 1].re = c1_a12.re - c1_a22.re;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2,
                            0) - 1)) - 1].im = c1_a12.im - c1_a22.im;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2,
                            0) - 1)) - 1].re = c1_a21.re;
                  c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2,
                            0) - 1)) - 1].im = c1_a21.im;
                }

                c1_jm1 = c1_b_j;
                c1_b_j = c1_jp1;
              }
            }

            guard11 = true;
          }

          if (guard11 == true) {
            c1_jiter++;
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
    if (c1_failed) {
      *c1_info = (real_T)c1_ilast;
      c1_b_ilast = c1_ilast;
      c1_m_b = c1_b_ilast;
      c1_n_b = c1_m_b;
      if (1 > c1_n_b) {
        c1_f_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_f_overflow = (c1_n_b > 2147483646);
      }

      if (c1_f_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_f_overflow);
      }

      for (c1_k = 1; c1_k <= c1_b_ilast; c1_k++) {
        c1_b_k = c1_k;
        c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_k), 1, 3, 1, 0) - 1].re = c1_dc2.re;
        c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_k), 1, 3, 1, 0) - 1].im = c1_dc2.im;
        c1_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_k), 1, 3, 1, 0) - 1].re = c1_dc2.re;
        c1_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_k), 1, 3, 1, 0) - 1].im = c1_dc2.im;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1 == true) {
    c1_pb_a = c1_ilo;
    c1_qb_a = c1_pb_a - 1;
    c1_i59 = c1_qb_a;
    c1_o_b = c1_i59;
    c1_p_b = c1_o_b;
    if (1 > c1_p_b) {
      c1_g_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_g_overflow = (c1_p_b > 2147483646);
    }

    if (c1_g_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_g_overflow);
    }

    for (c1_e_j = 1; c1_e_j <= c1_i59; c1_e_j++) {
      c1_b_j = c1_e_j;
      c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_j), 1, 3, 1, 0) - 1].re = c1_b_A
        [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) -
        1].re;
      c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_j), 1, 3, 1, 0) - 1].im = c1_b_A
        [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) -
        1].im;
    }

    *c1_info = 0.0;
  }
}

static real_T c1_eml_matlab_zlanhs(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T c1_ilo, int32_T c1_ihi)
{
  real_T c1_f;
  real_T c1_scale;
  real_T c1_sumsq;
  boolean_T c1_firstNonZero;
  int32_T c1_b_ilo;
  int32_T c1_b_ihi;
  int32_T c1_a;
  int32_T c1_b;
  int32_T c1_b_a;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_c_ilo;
  int32_T c1_c_a;
  int32_T c1_d_a;
  int32_T c1_c;
  int32_T c1_x;
  int32_T c1_y;
  int32_T c1_i60;
  int32_T c1_e_a;
  int32_T c1_c_b;
  int32_T c1_f_a;
  int32_T c1_d_b;
  boolean_T c1_b_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  creal_T c1_Aij;
  real_T c1_reAij;
  real_T c1_imAij;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_temp1;
  real_T c1_temp2;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_g_x;
  c1_f = 0.0;
  if (c1_ilo > c1_ihi) {
  } else {
    c1_scale = 0.0;
    c1_sumsq = 0.0;
    c1_firstNonZero = true;
    c1_b_ilo = c1_ilo;
    c1_b_ihi = c1_ihi;
    c1_a = c1_b_ilo;
    c1_b = c1_b_ihi;
    c1_b_a = c1_a;
    c1_b_b = c1_b;
    if (c1_b_a > c1_b_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_b_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_j = c1_b_ilo; c1_j <= c1_b_ihi; c1_j++) {
      c1_b_j = c1_j;
      c1_c_ilo = c1_ilo;
      c1_c_a = c1_b_j;
      c1_d_a = c1_c_a;
      c1_c = c1_d_a;
      c1_x = c1_c + 1;
      c1_y = c1_ihi;
      c1_i60 = c1_x;
      if (c1_y < c1_i60) {
        c1_i60 = c1_y;
      }

      c1_e_a = c1_c_ilo;
      c1_c_b = c1_i60;
      c1_f_a = c1_e_a;
      c1_d_b = c1_c_b;
      if (c1_f_a > c1_d_b) {
        c1_b_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_b_overflow = (c1_d_b > 2147483646);
      }

      if (c1_b_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
      }

      for (c1_i = c1_c_ilo; c1_i <= c1_i60; c1_i++) {
        c1_b_i = c1_i;
        c1_Aij.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
        c1_Aij.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
        c1_reAij = c1_Aij.re;
        c1_imAij = c1_Aij.im;
        if (c1_reAij != 0.0) {
          c1_b_x = c1_reAij;
          c1_c_x = c1_b_x;
          c1_temp1 = muDoubleScalarAbs(c1_c_x);
          if (c1_firstNonZero) {
            c1_sumsq = 1.0;
            c1_scale = c1_temp1;
            c1_firstNonZero = false;
          } else if (c1_scale < c1_temp1) {
            c1_temp2 = c1_scale / c1_temp1;
            c1_sumsq = 1.0 + c1_sumsq * c1_temp2 * c1_temp2;
            c1_scale = c1_temp1;
          } else {
            c1_temp2 = c1_temp1 / c1_scale;
            c1_sumsq += c1_temp2 * c1_temp2;
          }
        }

        if (c1_imAij != 0.0) {
          c1_d_x = c1_imAij;
          c1_e_x = c1_d_x;
          c1_temp1 = muDoubleScalarAbs(c1_e_x);
          if (c1_firstNonZero) {
            c1_sumsq = 1.0;
            c1_scale = c1_temp1;
            c1_firstNonZero = false;
          } else if (c1_scale < c1_temp1) {
            c1_temp2 = c1_scale / c1_temp1;
            c1_sumsq = 1.0 + c1_sumsq * c1_temp2 * c1_temp2;
            c1_scale = c1_temp1;
          } else {
            c1_temp2 = c1_temp1 / c1_scale;
            c1_sumsq += c1_temp2 * c1_temp2;
          }
        }
      }
    }

    c1_f_x = c1_sumsq;
    c1_g_x = c1_f_x;
    if (c1_g_x < 0.0) {
      c1_eml_error(chartInstance);
    }

    c1_g_x = muDoubleScalarSqrt(c1_g_x);
    c1_f = c1_scale * c1_g_x;
  }

  return c1_f;
}

static int32_T c1_mod(SFc1_Model_01InstanceStruct *chartInstance, int32_T c1_x)
{
  int32_T c1_b_x;
  int32_T c1_t;
  c1_b_x = c1_x;
  c1_t = c1_div_s32(chartInstance, c1_b_x, 10);
  c1_t *= 10;
  return c1_b_x - c1_t;
}

static creal_T c1_eml_div(SFc1_Model_01InstanceStruct *chartInstance, creal_T
  c1_x, real_T c1_y)
{
  creal_T c1_z;
  real_T c1_b_y;
  real_T c1_ar;
  real_T c1_ai;
  real_T c1_br;
  real_T c1_bi;
  real_T c1_brm;
  real_T c1_bim;
  real_T c1_s;
  real_T c1_d;
  real_T c1_nr;
  real_T c1_ni;
  real_T c1_sgnbr;
  real_T c1_sgnbi;
  (void)chartInstance;
  c1_b_y = c1_y;
  c1_ar = c1_x.re;
  c1_ai = c1_x.im;
  c1_br = c1_b_y;
  c1_bi = 0.0;
  if (c1_bi == 0.0) {
    if (c1_ai == 0.0) {
      c1_z.re = c1_ar / c1_br;
      c1_z.im = 0.0;
    } else if (c1_ar == 0.0) {
      c1_z.re = 0.0;
      c1_z.im = c1_ai / c1_br;
    } else {
      c1_z.re = c1_ar / c1_br;
      c1_z.im = c1_ai / c1_br;
    }
  } else if (c1_br == 0.0) {
    if (c1_ar == 0.0) {
      c1_z.re = c1_ai / c1_bi;
      c1_z.im = 0.0;
    } else if (c1_ai == 0.0) {
      c1_z.re = 0.0;
      c1_z.im = -(c1_ar / c1_bi);
    } else {
      c1_z.re = c1_ai / c1_bi;
      c1_z.im = -(c1_ar / c1_bi);
    }
  } else {
    c1_brm = muDoubleScalarAbs(c1_br);
    c1_bim = muDoubleScalarAbs(c1_bi);
    if (c1_brm > c1_bim) {
      c1_s = c1_bi / c1_br;
      c1_d = c1_br + c1_s * c1_bi;
      c1_nr = c1_ar + c1_s * c1_ai;
      c1_ni = c1_ai - c1_s * c1_ar;
      c1_z.re = c1_nr / c1_d;
      c1_z.im = c1_ni / c1_d;
    } else if (c1_bim == c1_brm) {
      if (c1_br > 0.0) {
        c1_sgnbr = 0.5;
      } else {
        c1_sgnbr = -0.5;
      }

      if (c1_bi > 0.0) {
        c1_sgnbi = 0.5;
      } else {
        c1_sgnbi = -0.5;
      }

      c1_nr = c1_ar * c1_sgnbr + c1_ai * c1_sgnbi;
      c1_ni = c1_ai * c1_sgnbr - c1_ar * c1_sgnbi;
      c1_z.re = c1_nr / c1_brm;
      c1_z.im = c1_ni / c1_brm;
    } else {
      c1_s = c1_br / c1_bi;
      c1_d = c1_bi + c1_s * c1_br;
      c1_nr = c1_s * c1_ar + c1_ai;
      c1_ni = c1_s * c1_ai - c1_ar;
      c1_z.re = c1_nr / c1_d;
      c1_z.im = c1_ni / c1_d;
    }
  }

  return c1_z;
}

static void c1_scalarEg(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static creal_T c1_sqrt(SFc1_Model_01InstanceStruct *chartInstance, creal_T c1_x)
{
  creal_T c1_b_x;
  c1_b_x = c1_x;
  c1_b_sqrt(chartInstance, &c1_b_x);
  return c1_b_x;
}

static void c1_realmax(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_b_eml_matlab_zlartg(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_f, creal_T c1_g, real_T *c1_cs, creal_T *c1_sn)
{
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_y;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_b_y;
  real_T c1_e_x;
  real_T c1_c_y;
  real_T c1_d_y;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_e_y;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_f_y;
  real_T c1_j_x;
  real_T c1_g_y;
  real_T c1_h_y;
  real_T c1_k_x;
  real_T c1_i_y;
  real_T c1_scale;
  creal_T c1_fs;
  creal_T c1_gs;
  int32_T c1_count;
  real_T c1_rescaledir;
  int32_T c1_a;
  int32_T c1_b_a;
  static creal_T c1_dc4 = { 0.0, 0.0 };

  boolean_T c1_b_g;
  int32_T c1_c_a;
  int32_T c1_d_a;
  real_T c1_f2;
  real_T c1_g2;
  real_T c1_l_x;
  real_T c1_m_x;
  boolean_T c1_b_f;
  real_T c1_x1;
  real_T c1_x2;
  real_T c1_e_a;
  real_T c1_b;
  real_T c1_d;
  real_T c1_b_x1;
  real_T c1_b_x2;
  real_T c1_f_a;
  real_T c1_b_b;
  real_T c1_f2s;
  real_T c1_n_x;
  real_T c1_g2s;
  real_T c1_o_x;
  real_T c1_p_x;
  real_T c1_j_y;
  real_T c1_q_x;
  real_T c1_r_x;
  real_T c1_k_y;
  real_T c1_s_x;
  real_T c1_l_y;
  real_T c1_m_y;
  real_T c1_c_x1;
  real_T c1_c_x2;
  real_T c1_g_a;
  real_T c1_c_b;
  real_T c1_dr;
  real_T c1_di;
  real_T c1_d_x1;
  real_T c1_d_x2;
  real_T c1_h_a;
  real_T c1_d_b;
  creal_T c1_b_gs;
  real_T c1_t_x;
  creal_T c1_b_fs;
  creal_T c1_c_fs;
  creal_T c1_c_gs;
  creal_T c1_b_sn;
  int32_T c1_b_count;
  int32_T c1_e_b;
  int32_T c1_f_b;
  boolean_T c1_overflow;
  int32_T c1_c_count;
  int32_T c1_g_b;
  int32_T c1_h_b;
  boolean_T c1_b_overflow;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  c1_realmin(chartInstance);
  c1_eps(chartInstance);
  c1_eps(chartInstance);
  c1_x = c1_f.re;
  c1_b_x = c1_x;
  c1_y = muDoubleScalarAbs(c1_b_x);
  c1_c_x = c1_f.im;
  c1_d_x = c1_c_x;
  c1_b_y = muDoubleScalarAbs(c1_d_x);
  c1_e_x = c1_y;
  c1_c_y = c1_b_y;
  c1_d_y = c1_e_x;
  if (c1_c_y > c1_d_y) {
    c1_d_y = c1_c_y;
  }

  c1_f_x = c1_g.re;
  c1_g_x = c1_f_x;
  c1_e_y = muDoubleScalarAbs(c1_g_x);
  c1_h_x = c1_g.im;
  c1_i_x = c1_h_x;
  c1_f_y = muDoubleScalarAbs(c1_i_x);
  c1_j_x = c1_e_y;
  c1_g_y = c1_f_y;
  c1_h_y = c1_j_x;
  if (c1_g_y > c1_h_y) {
    c1_h_y = c1_g_y;
  }

  c1_k_x = c1_d_y;
  c1_i_y = c1_h_y;
  c1_scale = c1_k_x;
  if (c1_i_y > c1_scale) {
    c1_scale = c1_i_y;
  }

  c1_fs = c1_f;
  c1_gs = c1_g;
  c1_count = 0;
  c1_rescaledir = 0.0;
  guard1 = false;
  guard2 = false;
  if (c1_scale >= 7.4428285367870146E+137) {
    do {
      c1_a = c1_count;
      c1_b_a = c1_a + 1;
      c1_count = c1_b_a;
      c1_fs.re *= 1.3435752215134178E-138;
      c1_fs.im *= 1.3435752215134178E-138;
      c1_gs.re *= 1.3435752215134178E-138;
      c1_gs.im *= 1.3435752215134178E-138;
      c1_scale *= 1.3435752215134178E-138;
    } while (!(c1_scale < 7.4428285367870146E+137));

    c1_rescaledir = 1.0;
    guard1 = true;
  } else if (c1_scale <= 1.3435752215134178E-138) {
    c1_b_g = ((c1_g.re == c1_dc4.re) && (c1_g.im == c1_dc4.im));
    if (c1_b_g) {
      *c1_cs = 1.0;
      *c1_sn = c1_dc4;
    } else {
      do {
        c1_c_a = c1_count;
        c1_d_a = c1_c_a + 1;
        c1_count = c1_d_a;
        c1_fs.re *= 7.4428285367870146E+137;
        c1_fs.im *= 7.4428285367870146E+137;
        c1_gs.re *= 7.4428285367870146E+137;
        c1_gs.im *= 7.4428285367870146E+137;
        c1_scale *= 7.4428285367870146E+137;
      } while (!(c1_scale > 1.3435752215134178E-138));

      c1_rescaledir = -1.0;
      guard2 = true;
    }
  } else {
    guard2 = true;
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c1_f2 = c1_fs.re * c1_fs.re + c1_fs.im * c1_fs.im;
    c1_g2 = c1_gs.re * c1_gs.re + c1_gs.im * c1_gs.im;
    c1_l_x = c1_g2;
    c1_m_x = c1_l_x;
    if (1.0 > c1_m_x) {
      c1_m_x = 1.0;
    }

    if (c1_f2 <= c1_m_x * 2.0041683600089728E-292) {
      c1_b_f = ((c1_f.re == c1_dc4.re) && (c1_f.im == c1_dc4.im));
      if (c1_b_f) {
        *c1_cs = 0.0;
        c1_x1 = c1_gs.re;
        c1_x2 = c1_gs.im;
        c1_e_a = c1_x1;
        c1_b = c1_x2;
        c1_d = muDoubleScalarHypot(c1_e_a, c1_b);
        c1_sn->re = c1_gs.re / c1_d;
        c1_sn->im = -c1_gs.im / c1_d;
      } else {
        c1_b_x1 = c1_fs.re;
        c1_b_x2 = c1_fs.im;
        c1_f_a = c1_b_x1;
        c1_b_b = c1_b_x2;
        c1_f2s = muDoubleScalarHypot(c1_f_a, c1_b_b);
        c1_n_x = c1_g2;
        c1_g2s = c1_n_x;
        if (c1_g2s < 0.0) {
          c1_eml_error(chartInstance);
        }

        c1_g2s = muDoubleScalarSqrt(c1_g2s);
        *c1_cs = c1_f2s / c1_g2s;
        c1_o_x = c1_f.re;
        c1_p_x = c1_o_x;
        c1_j_y = muDoubleScalarAbs(c1_p_x);
        c1_q_x = c1_f.im;
        c1_r_x = c1_q_x;
        c1_k_y = muDoubleScalarAbs(c1_r_x);
        c1_s_x = c1_j_y;
        c1_l_y = c1_k_y;
        c1_m_y = c1_s_x;
        if (c1_l_y > c1_m_y) {
          c1_m_y = c1_l_y;
        }

        if (c1_m_y > 1.0) {
          c1_c_x1 = c1_f.re;
          c1_c_x2 = c1_f.im;
          c1_g_a = c1_c_x1;
          c1_c_b = c1_c_x2;
          c1_d = muDoubleScalarHypot(c1_g_a, c1_c_b);
          c1_fs.re = c1_f.re / c1_d;
          c1_fs.im = c1_f.im / c1_d;
        } else {
          c1_dr = 7.4428285367870146E+137 * c1_f.re;
          c1_di = 7.4428285367870146E+137 * c1_f.im;
          c1_d_x1 = c1_dr;
          c1_d_x2 = c1_di;
          c1_h_a = c1_d_x1;
          c1_d_b = c1_d_x2;
          c1_d = muDoubleScalarHypot(c1_h_a, c1_d_b);
          c1_fs.re = c1_dr / c1_d;
          c1_fs.im = c1_di / c1_d;
        }

        c1_b_gs.re = c1_gs.re / c1_g2s;
        c1_b_gs.im = -c1_gs.im / c1_g2s;
        c1_sn->re = c1_fs.re * c1_b_gs.re - c1_fs.im * c1_b_gs.im;
        c1_sn->im = c1_fs.re * c1_b_gs.im + c1_fs.im * c1_b_gs.re;
      }
    } else {
      c1_t_x = 1.0 + c1_g2 / c1_f2;
      c1_f2s = c1_t_x;
      if (c1_f2s < 0.0) {
        c1_eml_error(chartInstance);
      }

      c1_f2s = muDoubleScalarSqrt(c1_f2s);
      c1_b_fs = c1_fs;
      c1_c_fs = c1_fs;
      c1_fs.re = c1_f2s * c1_b_fs.re;
      c1_fs.im = c1_f2s * c1_c_fs.im;
      *c1_cs = 1.0 / c1_f2s;
      c1_d = c1_f2 + c1_g2;
      c1_sn->re = c1_fs.re / c1_d;
      c1_sn->im = c1_fs.im / c1_d;
      c1_c_gs.re = c1_gs.re;
      c1_c_gs.im = -c1_gs.im;
      c1_b_sn = *c1_sn;
      c1_sn->re = c1_b_sn.re * c1_c_gs.re - c1_b_sn.im * c1_c_gs.im;
      c1_sn->im = c1_b_sn.re * c1_c_gs.im + c1_b_sn.im * c1_c_gs.re;
      if (c1_rescaledir > 0.0) {
        c1_b_count = c1_count;
        c1_e_b = c1_b_count;
        c1_f_b = c1_e_b;
        if (1 > c1_f_b) {
          c1_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_overflow = (c1_f_b > 2147483646);
        }

        if (c1_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_overflow);
        }
      } else {
        if (c1_rescaledir < 0.0) {
          c1_c_count = c1_count;
          c1_g_b = c1_c_count;
          c1_h_b = c1_g_b;
          if (1 > c1_h_b) {
            c1_b_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_b_overflow = (c1_h_b > 2147483646);
          }

          if (c1_b_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
          }
        }
      }
    }
  }
}

static void c1_b_eml_matlab_zlascl(SFc1_Model_01InstanceStruct *chartInstance,
  real_T c1_cfrom, real_T c1_cto, creal_T c1_A[3], creal_T c1_b_A[3])
{
  int32_T c1_i61;
  for (c1_i61 = 0; c1_i61 < 3; c1_i61++) {
    c1_b_A[c1_i61] = c1_A[c1_i61];
  }

  c1_d_eml_matlab_zlascl(chartInstance, c1_cfrom, c1_cto, c1_b_A);
}

static void c1_b_eml_div(SFc1_Model_01InstanceStruct *chartInstance, creal_T
  c1_x[3], creal_T c1_y[3], creal_T c1_z[3])
{
  int32_T c1_i62;
  real_T c1_ar;
  real_T c1_ai;
  real_T c1_br;
  real_T c1_bi;
  real_T c1_brm;
  real_T c1_bim;
  real_T c1_s;
  real_T c1_d;
  real_T c1_nr;
  real_T c1_ni;
  real_T c1_sgnbr;
  real_T c1_sgnbi;
  (void)chartInstance;
  for (c1_i62 = 0; c1_i62 < 3; c1_i62++) {
    c1_ar = c1_x[c1_i62].re;
    c1_ai = c1_x[c1_i62].im;
    c1_br = c1_y[c1_i62].re;
    c1_bi = c1_y[c1_i62].im;
    if (c1_bi == 0.0) {
      if (c1_ai == 0.0) {
        c1_z[c1_i62].re = c1_ar / c1_br;
        c1_z[c1_i62].im = 0.0;
      } else if (c1_ar == 0.0) {
        c1_z[c1_i62].re = 0.0;
        c1_z[c1_i62].im = c1_ai / c1_br;
      } else {
        c1_z[c1_i62].re = c1_ar / c1_br;
        c1_z[c1_i62].im = c1_ai / c1_br;
      }
    } else if (c1_br == 0.0) {
      if (c1_ar == 0.0) {
        c1_z[c1_i62].re = c1_ai / c1_bi;
        c1_z[c1_i62].im = 0.0;
      } else if (c1_ai == 0.0) {
        c1_z[c1_i62].re = 0.0;
        c1_z[c1_i62].im = -(c1_ar / c1_bi);
      } else {
        c1_z[c1_i62].re = c1_ai / c1_bi;
        c1_z[c1_i62].im = -(c1_ar / c1_bi);
      }
    } else {
      c1_brm = muDoubleScalarAbs(c1_br);
      c1_bim = muDoubleScalarAbs(c1_bi);
      if (c1_brm > c1_bim) {
        c1_s = c1_bi / c1_br;
        c1_d = c1_br + c1_s * c1_bi;
        c1_nr = c1_ar + c1_s * c1_ai;
        c1_ni = c1_ai - c1_s * c1_ar;
        c1_z[c1_i62].re = c1_nr / c1_d;
        c1_z[c1_i62].im = c1_ni / c1_d;
      } else if (c1_bim == c1_brm) {
        if (c1_br > 0.0) {
          c1_sgnbr = 0.5;
        } else {
          c1_sgnbr = -0.5;
        }

        if (c1_bi > 0.0) {
          c1_sgnbi = 0.5;
        } else {
          c1_sgnbi = -0.5;
        }

        c1_nr = c1_ar * c1_sgnbr + c1_ai * c1_sgnbi;
        c1_ni = c1_ai * c1_sgnbr - c1_ar * c1_sgnbi;
        c1_z[c1_i62].re = c1_nr / c1_brm;
        c1_z[c1_i62].im = c1_ni / c1_brm;
      } else {
        c1_s = c1_br / c1_bi;
        c1_d = c1_bi + c1_s * c1_br;
        c1_nr = c1_s * c1_ar + c1_ai;
        c1_ni = c1_s * c1_ai - c1_ar;
        c1_z[c1_i62].re = c1_nr / c1_d;
        c1_z[c1_i62].im = c1_ni / c1_d;
      }
    }
  }
}

static void c1_eml_warning(SFc1_Model_01InstanceStruct *chartInstance)
{
  int32_T c1_i63;
  static char_T c1_varargin_1[26] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'f', 'a', 'i',
    'l', 'e', 'd' };

  char_T c1_u[26];
  const mxArray *c1_y = NULL;
  (void)chartInstance;
  for (c1_i63 = 0; c1_i63 < 26; c1_i63++) {
    c1_u[c1_i63] = c1_varargin_1[c1_i63];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 26), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c1_y));
}

static void c1_b_eml_warning(SFc1_Model_01InstanceStruct *chartInstance)
{
  int32_T c1_i64;
  static char_T c1_varargin_1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'n', 'o', 'n',
    'c', 'o', 'n', 'v', 'e', 'r', 'g', 'e', 'n', 'c', 'e' };

  char_T c1_u[34];
  const mxArray *c1_y = NULL;
  (void)chartInstance;
  for (c1_i64 = 0; c1_i64 < 34; c1_i64++) {
    c1_u[c1_i64] = c1_varargin_1[c1_i64];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 34), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c1_y));
}

static real_T c1_mpower(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_a)
{
  real_T c1_b_a;
  real_T c1_c_a;
  real_T c1_ak;
  real_T c1_d_a;
  c1_b_a = c1_a;
  c1_c_a = c1_b_a;
  c1_eml_scalar_eg(chartInstance);
  c1_ak = c1_c_a;
  c1_d_a = c1_ak;
  c1_eml_scalar_eg(chartInstance);
  return c1_d_a * c1_d_a;
}

static void c1_inv(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_x[9],
                   real_T c1_y[9])
{
  int32_T c1_i65;
  real_T c1_b_x[9];
  int32_T c1_i66;
  real_T c1_c_x[9];
  real_T c1_n1x;
  int32_T c1_i67;
  real_T c1_b_y[9];
  real_T c1_n1xinv;
  real_T c1_rc;
  real_T c1_d_x;
  boolean_T c1_b;
  real_T c1_e_x;
  int32_T c1_i68;
  static char_T c1_cv4[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T c1_u[8];
  const mxArray *c1_c_y = NULL;
  real_T c1_b_u;
  const mxArray *c1_d_y = NULL;
  real_T c1_c_u;
  const mxArray *c1_e_y = NULL;
  real_T c1_d_u;
  const mxArray *c1_f_y = NULL;
  char_T c1_str[14];
  int32_T c1_i69;
  char_T c1_b_str[14];
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  for (c1_i65 = 0; c1_i65 < 9; c1_i65++) {
    c1_b_x[c1_i65] = c1_x[c1_i65];
  }

  c1_inv3x3(chartInstance, c1_b_x, c1_y);
  for (c1_i66 = 0; c1_i66 < 9; c1_i66++) {
    c1_c_x[c1_i66] = c1_x[c1_i66];
  }

  c1_n1x = c1_norm(chartInstance, c1_c_x);
  for (c1_i67 = 0; c1_i67 < 9; c1_i67++) {
    c1_b_y[c1_i67] = c1_y[c1_i67];
  }

  c1_n1xinv = c1_norm(chartInstance, c1_b_y);
  c1_rc = 1.0 / (c1_n1x * c1_n1xinv);
  guard1 = false;
  guard2 = false;
  if (c1_n1x == 0.0) {
    guard2 = true;
  } else if (c1_n1xinv == 0.0) {
    guard2 = true;
  } else if (c1_rc == 0.0) {
    guard1 = true;
  } else {
    c1_d_x = c1_rc;
    c1_b = muDoubleScalarIsNaN(c1_d_x);
    guard3 = false;
    if (c1_b) {
      guard3 = true;
    } else {
      c1_eps(chartInstance);
      if (c1_rc < 2.2204460492503131E-16) {
        guard3 = true;
      }
    }

    if (guard3 == true) {
      c1_e_x = c1_rc;
      for (c1_i68 = 0; c1_i68 < 8; c1_i68++) {
        c1_u[c1_i68] = c1_cv4[c1_i68];
      }

      c1_c_y = NULL;
      sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 8),
                    false);
      c1_b_u = 14.0;
      c1_d_y = NULL;
      sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0),
                    false);
      c1_c_u = 6.0;
      c1_e_y = NULL;
      sf_mex_assign(&c1_e_y, sf_mex_create("y", &c1_c_u, 0, 0U, 0U, 0U, 0),
                    false);
      c1_d_u = c1_e_x;
      c1_f_y = NULL;
      sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_d_u, 0, 0U, 0U, 0U, 0),
                    false);
      c1_e_emlrt_marshallIn(chartInstance, sf_mex_call_debug
                            (sfGlobalDebugInstanceStruct, "sprintf", 1U, 2U, 14,
        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "sprintf", 1U, 3U, 14,
                          c1_c_y, 14, c1_d_y, 14, c1_e_y), 14, c1_f_y),
                            "sprintf", c1_str);
      for (c1_i69 = 0; c1_i69 < 14; c1_i69++) {
        c1_b_str[c1_i69] = c1_str[c1_i69];
      }

      c1_d_eml_warning(chartInstance, c1_b_str);
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c1_c_eml_warning(chartInstance);
  }
}

static void c1_inv3x3(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_x[9],
                      real_T c1_y[9])
{
  int32_T c1_p1;
  int32_T c1_p2;
  int32_T c1_p3;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_absx11;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_absx21;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_absx31;
  real_T c1_t1;
  real_T c1_h_x;
  real_T c1_b_y;
  real_T c1_i_x;
  real_T c1_c_y;
  real_T c1_z;
  real_T c1_j_x;
  real_T c1_d_y;
  real_T c1_k_x;
  real_T c1_e_y;
  real_T c1_b_z;
  real_T c1_l_x;
  real_T c1_m_x;
  real_T c1_f_y;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_g_y;
  int32_T c1_itmp;
  real_T c1_p_x;
  real_T c1_h_y;
  real_T c1_q_x;
  real_T c1_i_y;
  real_T c1_c_z;
  real_T c1_r_x;
  real_T c1_j_y;
  real_T c1_s_x;
  real_T c1_k_y;
  real_T c1_t3;
  real_T c1_t_x;
  real_T c1_l_y;
  real_T c1_u_x;
  real_T c1_m_y;
  real_T c1_t2;
  int32_T c1_a;
  int32_T c1_b_a;
  int32_T c1_c;
  real_T c1_v_x;
  real_T c1_n_y;
  real_T c1_w_x;
  real_T c1_o_y;
  real_T c1_d_z;
  int32_T c1_c_a;
  int32_T c1_d_a;
  int32_T c1_b_c;
  int32_T c1_e_a;
  int32_T c1_f_a;
  int32_T c1_c_c;
  real_T c1_x_x;
  real_T c1_p_y;
  real_T c1_y_x;
  real_T c1_q_y;
  real_T c1_ab_x;
  real_T c1_r_y;
  real_T c1_bb_x;
  real_T c1_s_y;
  int32_T c1_g_a;
  int32_T c1_h_a;
  int32_T c1_d_c;
  real_T c1_cb_x;
  real_T c1_t_y;
  real_T c1_db_x;
  real_T c1_u_y;
  real_T c1_e_z;
  int32_T c1_i_a;
  int32_T c1_j_a;
  int32_T c1_e_c;
  int32_T c1_k_a;
  int32_T c1_l_a;
  int32_T c1_f_c;
  real_T c1_v_y;
  real_T c1_w_y;
  real_T c1_eb_x;
  real_T c1_x_y;
  real_T c1_fb_x;
  real_T c1_y_y;
  int32_T c1_m_a;
  int32_T c1_n_a;
  int32_T c1_g_c;
  real_T c1_gb_x;
  real_T c1_ab_y;
  real_T c1_hb_x;
  real_T c1_bb_y;
  real_T c1_f_z;
  int32_T c1_o_a;
  int32_T c1_p_a;
  int32_T c1_h_c;
  int32_T c1_q_a;
  int32_T c1_r_a;
  int32_T c1_i_c;
  boolean_T guard1 = false;
  (void)chartInstance;
  c1_p1 = 0;
  c1_p2 = 3;
  c1_p3 = 6;
  c1_b_x = c1_x[0];
  c1_c_x = c1_b_x;
  c1_absx11 = muDoubleScalarAbs(c1_c_x);
  c1_d_x = c1_x[1];
  c1_e_x = c1_d_x;
  c1_absx21 = muDoubleScalarAbs(c1_e_x);
  c1_f_x = c1_x[2];
  c1_g_x = c1_f_x;
  c1_absx31 = muDoubleScalarAbs(c1_g_x);
  guard1 = false;
  if (c1_absx21 > c1_absx11) {
    if (c1_absx21 > c1_absx31) {
      c1_p1 = 3;
      c1_p2 = 0;
      c1_t1 = c1_x[0];
      c1_x[0] = c1_x[1];
      c1_x[1] = c1_t1;
      c1_t1 = c1_x[3];
      c1_x[3] = c1_x[4];
      c1_x[4] = c1_t1;
      c1_t1 = c1_x[6];
      c1_x[6] = c1_x[7];
      c1_x[7] = c1_t1;
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1 == true) {
    if (c1_absx31 > c1_absx11) {
      c1_p1 = 6;
      c1_p3 = 0;
      c1_t1 = c1_x[0];
      c1_x[0] = c1_x[2];
      c1_x[2] = c1_t1;
      c1_t1 = c1_x[3];
      c1_x[3] = c1_x[5];
      c1_x[5] = c1_t1;
      c1_t1 = c1_x[6];
      c1_x[6] = c1_x[8];
      c1_x[8] = c1_t1;
    }
  }

  c1_h_x = c1_x[1];
  c1_b_y = c1_x[0];
  c1_i_x = c1_h_x;
  c1_c_y = c1_b_y;
  c1_z = c1_i_x / c1_c_y;
  c1_x[1] = c1_z;
  c1_j_x = c1_x[2];
  c1_d_y = c1_x[0];
  c1_k_x = c1_j_x;
  c1_e_y = c1_d_y;
  c1_b_z = c1_k_x / c1_e_y;
  c1_x[2] = c1_b_z;
  c1_x[4] -= c1_x[1] * c1_x[3];
  c1_x[5] -= c1_x[2] * c1_x[3];
  c1_x[7] -= c1_x[1] * c1_x[6];
  c1_x[8] -= c1_x[2] * c1_x[6];
  c1_l_x = c1_x[5];
  c1_m_x = c1_l_x;
  c1_f_y = muDoubleScalarAbs(c1_m_x);
  c1_n_x = c1_x[4];
  c1_o_x = c1_n_x;
  c1_g_y = muDoubleScalarAbs(c1_o_x);
  if (c1_f_y > c1_g_y) {
    c1_itmp = c1_p2;
    c1_p2 = c1_p3;
    c1_p3 = c1_itmp;
    c1_t1 = c1_x[1];
    c1_x[1] = c1_x[2];
    c1_x[2] = c1_t1;
    c1_t1 = c1_x[4];
    c1_x[4] = c1_x[5];
    c1_x[5] = c1_t1;
    c1_t1 = c1_x[7];
    c1_x[7] = c1_x[8];
    c1_x[8] = c1_t1;
  }

  c1_p_x = c1_x[5];
  c1_h_y = c1_x[4];
  c1_q_x = c1_p_x;
  c1_i_y = c1_h_y;
  c1_c_z = c1_q_x / c1_i_y;
  c1_x[5] = c1_c_z;
  c1_x[8] -= c1_x[5] * c1_x[7];
  c1_r_x = c1_x[5] * c1_x[1] - c1_x[2];
  c1_j_y = c1_x[8];
  c1_s_x = c1_r_x;
  c1_k_y = c1_j_y;
  c1_t3 = c1_s_x / c1_k_y;
  c1_t_x = -(c1_x[1] + c1_x[7] * c1_t3);
  c1_l_y = c1_x[4];
  c1_u_x = c1_t_x;
  c1_m_y = c1_l_y;
  c1_t2 = c1_u_x / c1_m_y;
  c1_a = c1_p1;
  c1_b_a = c1_a + 1;
  c1_c = c1_b_a;
  c1_v_x = (1.0 - c1_x[3] * c1_t2) - c1_x[6] * c1_t3;
  c1_n_y = c1_x[0];
  c1_w_x = c1_v_x;
  c1_o_y = c1_n_y;
  c1_d_z = c1_w_x / c1_o_y;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_c), 1, 9, 1, 0) - 1] = c1_d_z;
  c1_c_a = c1_p1;
  c1_d_a = c1_c_a + 2;
  c1_b_c = c1_d_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_b_c), 1, 9, 1, 0) - 1] = c1_t2;
  c1_e_a = c1_p1;
  c1_f_a = c1_e_a + 3;
  c1_c_c = c1_f_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_c_c), 1, 9, 1, 0) - 1] = c1_t3;
  c1_x_x = -c1_x[5];
  c1_p_y = c1_x[8];
  c1_y_x = c1_x_x;
  c1_q_y = c1_p_y;
  c1_t3 = c1_y_x / c1_q_y;
  c1_ab_x = 1.0 - c1_x[7] * c1_t3;
  c1_r_y = c1_x[4];
  c1_bb_x = c1_ab_x;
  c1_s_y = c1_r_y;
  c1_t2 = c1_bb_x / c1_s_y;
  c1_g_a = c1_p2;
  c1_h_a = c1_g_a + 1;
  c1_d_c = c1_h_a;
  c1_cb_x = -(c1_x[3] * c1_t2 + c1_x[6] * c1_t3);
  c1_t_y = c1_x[0];
  c1_db_x = c1_cb_x;
  c1_u_y = c1_t_y;
  c1_e_z = c1_db_x / c1_u_y;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_d_c), 1, 9, 1, 0) - 1] = c1_e_z;
  c1_i_a = c1_p2;
  c1_j_a = c1_i_a + 2;
  c1_e_c = c1_j_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_e_c), 1, 9, 1, 0) - 1] = c1_t2;
  c1_k_a = c1_p2;
  c1_l_a = c1_k_a + 3;
  c1_f_c = c1_l_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_f_c), 1, 9, 1, 0) - 1] = c1_t3;
  c1_v_y = c1_x[8];
  c1_w_y = c1_v_y;
  c1_t3 = 1.0 / c1_w_y;
  c1_eb_x = -c1_x[7] * c1_t3;
  c1_x_y = c1_x[4];
  c1_fb_x = c1_eb_x;
  c1_y_y = c1_x_y;
  c1_t2 = c1_fb_x / c1_y_y;
  c1_m_a = c1_p3;
  c1_n_a = c1_m_a + 1;
  c1_g_c = c1_n_a;
  c1_gb_x = -(c1_x[3] * c1_t2 + c1_x[6] * c1_t3);
  c1_ab_y = c1_x[0];
  c1_hb_x = c1_gb_x;
  c1_bb_y = c1_ab_y;
  c1_f_z = c1_hb_x / c1_bb_y;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_g_c), 1, 9, 1, 0) - 1] = c1_f_z;
  c1_o_a = c1_p3;
  c1_p_a = c1_o_a + 2;
  c1_h_c = c1_p_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_h_c), 1, 9, 1, 0) - 1] = c1_t2;
  c1_q_a = c1_p3;
  c1_r_a = c1_q_a + 3;
  c1_i_c = c1_r_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_i_c), 1, 9, 1, 0) - 1] = c1_t3;
}

static real_T c1_norm(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_x[9])
{
  real_T c1_y;
  int32_T c1_j;
  real_T c1_b_j;
  real_T c1_s;
  int32_T c1_i;
  real_T c1_b_i;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_b_y;
  real_T c1_d_x;
  boolean_T c1_b;
  boolean_T exitg1;
  (void)chartInstance;
  c1_y = 0.0;
  c1_j = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c1_j < 3)) {
    c1_b_j = 1.0 + (real_T)c1_j;
    c1_s = 0.0;
    for (c1_i = 0; c1_i < 3; c1_i++) {
      c1_b_i = 1.0 + (real_T)c1_i;
      c1_b_x = c1_x[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_i), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1];
      c1_c_x = c1_b_x;
      c1_b_y = muDoubleScalarAbs(c1_c_x);
      c1_s += c1_b_y;
    }

    c1_d_x = c1_s;
    c1_b = muDoubleScalarIsNaN(c1_d_x);
    if (c1_b) {
      c1_y = rtNaN;
      exitg1 = true;
    } else {
      if (c1_s > c1_y) {
        c1_y = c1_s;
      }

      c1_j++;
    }
  }

  return c1_y;
}

static void c1_c_eml_warning(SFc1_Model_01InstanceStruct *chartInstance)
{
  int32_T c1_i70;
  static char_T c1_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c1_u[27];
  const mxArray *c1_y = NULL;
  (void)chartInstance;
  for (c1_i70 = 0; c1_i70 < 27; c1_i70++) {
    c1_u[c1_i70] = c1_varargin_1[c1_i70];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 27), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c1_y));
}

static void c1_d_eml_warning(SFc1_Model_01InstanceStruct *chartInstance, char_T
  c1_varargin_2[14])
{
  int32_T c1_i71;
  static char_T c1_varargin_1[33] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 'i', 'l', 'l', 'C', 'o', 'n', 'd', 'i', 't', 'i',
    'o', 'n', 'e', 'd', 'M', 'a', 't', 'r', 'i', 'x' };

  char_T c1_u[33];
  const mxArray *c1_y = NULL;
  int32_T c1_i72;
  char_T c1_b_u[14];
  const mxArray *c1_b_y = NULL;
  (void)chartInstance;
  for (c1_i71 = 0; c1_i71 < 33; c1_i71++) {
    c1_u[c1_i71] = c1_varargin_1[c1_i71];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 33), false);
  for (c1_i72 = 0; c1_i72 < 14; c1_i72++) {
    c1_b_u[c1_i72] = c1_varargin_2[c1_i72];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 14),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c1_y, 14, c1_b_y));
}

static void c1_b_eml_switch_helper(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c1_fn_phi_jk(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_k, real_T c1_lambda, real_T c1_t_delta, real_T c1_omega_norm)
{
  real_T c1_phi_jk;
  uint32_T c1_debug_family_var_map[7];
  real_T c1_nargin = 4.0;
  real_T c1_nargout = 1.0;
  real_T c1_b_lambda;
  real_T c1_b_t_delta;
  uint32_T c1_b_debug_family_var_map[5];
  real_T c1_b_nargin = 2.0;
  real_T c1_b_nargout = 1.0;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_c_lambda;
  real_T c1_b_omega_norm;
  real_T c1_c_t_delta;
  uint32_T c1_c_debug_family_var_map[6];
  real_T c1_c_nargin = 3.0;
  real_T c1_c_nargout = 1.0;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_d_lambda;
  real_T c1_c_omega_norm;
  real_T c1_d_t_delta;
  real_T c1_d_nargin = 3.0;
  real_T c1_d_nargout = 1.0;
  real_T c1_a;
  real_T c1_b_a;
  real_T c1_c_a;
  real_T c1_ak;
  real_T c1_d_a;
  real_T c1_ar;
  real_T c1_c;
  real_T c1_e_a;
  real_T c1_f_a;
  real_T c1_g_a;
  real_T c1_b_ak;
  real_T c1_h_a;
  real_T c1_b_ar;
  real_T c1_b_c;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_l_x;
  real_T c1_m_x;
  real_T c1_n_x;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 7U, 7U, c1_f_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 0U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 1U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_k, 2U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_lambda, 3U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t_delta, 4U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_omega_norm, 5U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_phi_jk, 6U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  CV_EML_FCN(0, 1);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 21);
  switch ((int32_T)_SFD_INTEGER_CHECK("k", c1_k)) {
   case 1:
    CV_EML_SWITCH(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 23);
    c1_b_lambda = c1_lambda;
    c1_b_t_delta = c1_t_delta;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 5U, 5U, c1_e_debug_family_names,
      c1_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargin, 0U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargout, 1U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_lambda, 2U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_t_delta, 3U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_phi_jk, 4U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    CV_EML_FCN(0, 2);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 35);
    c1_x = c1_b_lambda * c1_b_t_delta;
    c1_b_x = c1_x;
    c1_b_x = muDoubleScalarExp(c1_b_x);
    c1_phi_jk = c1_b_mpower(chartInstance, c1_b_lambda) * (c1_b_x - 1.0);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -35);
    _SFD_SYMBOL_SCOPE_POP();
    break;

   case 2:
    CV_EML_SWITCH(0, 1, 0, 2);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 25);
    c1_c_lambda = c1_lambda;
    c1_b_omega_norm = c1_omega_norm;
    c1_c_t_delta = c1_t_delta;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c1_d_debug_family_names,
      c1_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_c_nargin, 0U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_c_nargout, 1U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_c_lambda, 2U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_omega_norm, 3U,
      c1_b_sf_marshallOut, c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_c_t_delta, 4U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_phi_jk, 5U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    CV_EML_FCN(0, 3);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 38);
    c1_c_x = c1_b_omega_norm * c1_c_t_delta;
    c1_d_x = c1_c_x;
    c1_d_x = muDoubleScalarCos(c1_d_x);
    c1_e_x = c1_b_omega_norm * c1_c_t_delta;
    c1_f_x = c1_e_x;
    c1_f_x = muDoubleScalarSin(c1_f_x);
    c1_g_x = c1_c_lambda * c1_c_t_delta;
    c1_h_x = c1_g_x;
    c1_h_x = muDoubleScalarExp(c1_h_x);
    c1_phi_jk = c1_b_mpower(chartInstance, c1_mpower(chartInstance, c1_c_lambda)
      * c1_b_omega_norm + c1_c_mpower(chartInstance, c1_b_omega_norm)) *
      ((c1_b_omega_norm * c1_d_x + c1_c_lambda * c1_f_x) - c1_h_x *
       c1_b_omega_norm);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -38);
    _SFD_SYMBOL_SCOPE_POP();
    break;

   case 3:
    CV_EML_SWITCH(0, 1, 0, 3);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 27);
    c1_d_lambda = c1_lambda;
    c1_c_omega_norm = c1_omega_norm;
    c1_d_t_delta = c1_t_delta;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c1_c_debug_family_names,
      c1_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_d_nargin, 0U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_d_nargout, 1U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_d_lambda, 2U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_c_omega_norm, 3U,
      c1_b_sf_marshallOut, c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_d_t_delta, 4U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_phi_jk, 5U, c1_b_sf_marshallOut,
      c1_b_sf_marshallIn);
    CV_EML_FCN(0, 4);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 41);
    c1_a = c1_c_omega_norm;
    c1_b_a = c1_a;
    c1_c_a = c1_b_a;
    c1_eml_scalar_eg(chartInstance);
    c1_ak = c1_c_a;
    c1_d_a = c1_ak;
    c1_eml_scalar_eg(chartInstance);
    c1_ar = c1_d_a;
    c1_c = muDoubleScalarPower(c1_ar, -2.0);
    c1_e_a = c1_c_omega_norm;
    c1_f_a = c1_e_a;
    c1_g_a = c1_f_a;
    c1_eml_scalar_eg(chartInstance);
    c1_b_ak = c1_g_a;
    c1_h_a = c1_b_ak;
    c1_eml_scalar_eg(chartInstance);
    c1_b_ar = c1_h_a;
    c1_b_c = muDoubleScalarPower(c1_b_ar, 4.0);
    c1_i_x = c1_c_omega_norm * c1_d_t_delta;
    c1_j_x = c1_i_x;
    c1_j_x = muDoubleScalarCos(c1_j_x);
    c1_k_x = c1_c_omega_norm * c1_d_t_delta;
    c1_l_x = c1_k_x;
    c1_l_x = muDoubleScalarSin(c1_l_x);
    c1_m_x = c1_c_omega_norm * c1_d_t_delta;
    c1_n_x = c1_m_x;
    c1_n_x = muDoubleScalarExp(c1_n_x);
    c1_phi_jk = c1_c * c1_b_mpower(chartInstance, c1_d_lambda) + (c1_c_mpower
      (chartInstance, c1_d_lambda) * c1_mpower(chartInstance, c1_c_omega_norm) +
      c1_d_lambda * c1_b_c) * ((c1_mpower(chartInstance, c1_d_lambda) * c1_j_x -
      c1_c_omega_norm * c1_d_lambda * c1_l_x) + c1_mpower(chartInstance,
      c1_c_omega_norm) * c1_n_x);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -41);
    _SFD_SYMBOL_SCOPE_POP();
    break;

   default:
    CV_EML_SWITCH(0, 1, 0, 0);
    _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 29);
    c1_phi_jk = 0.0;
    break;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -29);
  _SFD_SYMBOL_SCOPE_POP();
  return c1_phi_jk;
}

static real_T c1_b_mpower(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_a)
{
  real_T c1_b_a;
  real_T c1_c_a;
  real_T c1_ak;
  real_T c1_d_a;
  real_T c1_B;
  real_T c1_y;
  real_T c1_b_y;
  real_T c1_c_y;
  c1_b_a = c1_a;
  c1_c_a = c1_b_a;
  c1_eml_scalar_eg(chartInstance);
  c1_ak = c1_c_a;
  c1_d_a = c1_ak;
  c1_eml_scalar_eg(chartInstance);
  c1_B = c1_d_a;
  c1_y = c1_B;
  c1_b_y = c1_y;
  c1_c_y = c1_b_y;
  return 1.0 / c1_c_y;
}

static real_T c1_c_mpower(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_a)
{
  real_T c1_b_a;
  real_T c1_c_a;
  real_T c1_ak;
  real_T c1_d_a;
  real_T c1_ar;
  c1_b_a = c1_a;
  c1_c_a = c1_b_a;
  c1_eml_scalar_eg(chartInstance);
  c1_ak = c1_c_a;
  c1_d_a = c1_ak;
  c1_eml_scalar_eg(chartInstance);
  c1_ar = c1_d_a;
  return muDoubleScalarPower(c1_ar, 3.0);
}

static void c1_d_mpower(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_a
  [9], real_T c1_b, real_T c1_c[9])
{
  real_T c1_x;
  real_T c1_b_x;
  int32_T c1_i73;
  real_T c1_b_a[9];
  real_T c1_b_b;
  int32_T c1_i74;
  int32_T c1_k;
  real_T c1_b_k;
  real_T c1_c_x;
  real_T c1_e;
  boolean_T c1_firstmult;
  real_T c1_d_x;
  real_T c1_ed2;
  int32_T c1_i75;
  int32_T c1_i76;
  real_T c1_c_a[9];
  int32_T c1_i77;
  int32_T c1_i78;
  real_T c1_d_a[9];
  int32_T c1_i79;
  real_T c1_e_a[9];
  int32_T c1_i80;
  real_T c1_b_c[9];
  int32_T c1_i81;
  int32_T c1_i82;
  int32_T c1_i83;
  real_T c1_f_a[9];
  int32_T c1_i84;
  real_T c1_g_a[9];
  real_T c1_c_b;
  int32_T c1_i85;
  real_T c1_h_a[9];
  creal_T c1_D[9];
  creal_T c1_V[9];
  int32_T c1_j;
  real_T c1_b_j;
  creal_T c1_b_D;
  creal_T c1_r;
  int32_T c1_i;
  real_T c1_b_i;
  int32_T c1_i86;
  creal_T c1_b_V[9];
  int32_T c1_i87;
  creal_T c1_c_D[9];
  int32_T c1_i88;
  int32_T exitg1;
  c1_x = c1_b;
  c1_b_x = c1_x;
  c1_b_x = muDoubleScalarFloor(c1_b_x);
  if (c1_b_x == c1_b) {
    for (c1_i73 = 0; c1_i73 < 9; c1_i73++) {
      c1_b_a[c1_i73] = c1_a[c1_i73];
    }

    c1_b_b = c1_b;
    c1_b_eml_scalar_eg(chartInstance);
    if (c1_b_b == 0.0) {
      for (c1_i74 = 0; c1_i74 < 9; c1_i74++) {
        c1_c[c1_i74] = 0.0;
      }

      for (c1_k = 0; c1_k < 3; c1_k++) {
        c1_b_k = 1.0 + (real_T)c1_k;
        c1_c[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c1_b_k), 1, 3, 2, 0) - 1)) - 1] =
          1.0;
      }
    } else {
      c1_c_x = c1_b_b;
      c1_e = muDoubleScalarAbs(c1_c_x);
      c1_firstmult = true;
      do {
        exitg1 = 0;
        c1_d_x = c1_e / 2.0;
        c1_ed2 = c1_d_x;
        c1_ed2 = muDoubleScalarFloor(c1_ed2);
        if (2.0 * c1_ed2 != c1_e) {
          if (c1_firstmult) {
            for (c1_i75 = 0; c1_i75 < 9; c1_i75++) {
              c1_c[c1_i75] = c1_b_a[c1_i75];
            }

            c1_firstmult = false;
          } else {
            for (c1_i76 = 0; c1_i76 < 9; c1_i76++) {
              c1_c_a[c1_i76] = c1_c[c1_i76];
            }

            c1_c_eml_scalar_eg(chartInstance);
            c1_c_eml_scalar_eg(chartInstance);
            for (c1_i77 = 0; c1_i77 < 9; c1_i77++) {
              c1_c[c1_i77] = 0.0;
            }

            for (c1_i78 = 0; c1_i78 < 9; c1_i78++) {
              c1_d_a[c1_i78] = c1_c_a[c1_i78];
            }

            for (c1_i79 = 0; c1_i79 < 9; c1_i79++) {
              c1_e_a[c1_i79] = c1_b_a[c1_i79];
            }

            c1_b_eml_xgemm(chartInstance, c1_d_a, c1_e_a, c1_c);
          }
        }

        if (c1_ed2 == 0.0) {
          exitg1 = 1;
        } else {
          c1_e = c1_ed2;
          for (c1_i81 = 0; c1_i81 < 9; c1_i81++) {
            c1_c_a[c1_i81] = c1_b_a[c1_i81];
          }

          c1_c_eml_scalar_eg(chartInstance);
          c1_c_eml_scalar_eg(chartInstance);
          for (c1_i82 = 0; c1_i82 < 9; c1_i82++) {
            c1_b_a[c1_i82] = 0.0;
          }

          for (c1_i83 = 0; c1_i83 < 9; c1_i83++) {
            c1_f_a[c1_i83] = c1_c_a[c1_i83];
          }

          for (c1_i84 = 0; c1_i84 < 9; c1_i84++) {
            c1_g_a[c1_i84] = c1_c_a[c1_i84];
          }

          c1_b_eml_xgemm(chartInstance, c1_f_a, c1_g_a, c1_b_a);
        }
      } while (exitg1 == 0);

      if (c1_b_b < 0.0) {
        for (c1_i80 = 0; c1_i80 < 9; c1_i80++) {
          c1_b_c[c1_i80] = c1_c[c1_i80];
        }

        c1_inv(chartInstance, c1_b_c, c1_c);
      }
    }
  } else {
    c1_c_b = c1_b;
    c1_b_eml_scalar_eg(chartInstance);
    c1_b_eml_error(chartInstance);
    for (c1_i85 = 0; c1_i85 < 9; c1_i85++) {
      c1_h_a[c1_i85] = c1_a[c1_i85];
    }

    c1_b_eig(chartInstance, c1_h_a, c1_V, c1_D);
    for (c1_j = 0; c1_j < 3; c1_j++) {
      c1_b_j = 1.0 + (real_T)c1_j;
      c1_b_D.re = c1_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) + 3 *
                        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c1_b_D.im = c1_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) + 3 *
                        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
      c1_r = c1_power(chartInstance, c1_b_D, c1_c_b);
      for (c1_i = 0; c1_i < 3; c1_i++) {
        c1_b_i = 1.0 + (real_T)c1_i;
        c1_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c1_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].
          re = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_i), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].re * c1_r.re -
          c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  c1_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                  (int32_T)_SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1]
          .im * c1_r.im;
        c1_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c1_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].
          im = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_i), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].re * c1_r.im +
          c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  c1_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                  (int32_T)_SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1]
          .im * c1_r.re;
      }
    }

    for (c1_i86 = 0; c1_i86 < 9; c1_i86++) {
      c1_b_V[c1_i86] = c1_V[c1_i86];
    }

    for (c1_i87 = 0; c1_i87 < 9; c1_i87++) {
      c1_c_D[c1_i87] = c1_D[c1_i87];
    }

    c1_eml_lusolve(chartInstance, c1_b_V, c1_c_D, c1_D);
    for (c1_i88 = 0; c1_i88 < 9; c1_i88++) {
      c1_c[c1_i88] = c1_D[c1_i88].re;
    }
  }
}

static void c1_b_eml_scalar_eg(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_c_eml_scalar_eg(SFc1_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_eml_xgemm(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_A[9], real_T c1_B[9], real_T c1_C[9], real_T c1_b_C[9])
{
  int32_T c1_i89;
  int32_T c1_i90;
  real_T c1_b_A[9];
  int32_T c1_i91;
  real_T c1_b_B[9];
  for (c1_i89 = 0; c1_i89 < 9; c1_i89++) {
    c1_b_C[c1_i89] = c1_C[c1_i89];
  }

  for (c1_i90 = 0; c1_i90 < 9; c1_i90++) {
    c1_b_A[c1_i90] = c1_A[c1_i90];
  }

  for (c1_i91 = 0; c1_i91 < 9; c1_i91++) {
    c1_b_B[c1_i91] = c1_B[c1_i91];
  }

  c1_b_eml_xgemm(chartInstance, c1_b_A, c1_b_B, c1_b_C);
}

static void c1_b_eml_error(SFc1_Model_01InstanceStruct *chartInstance)
{
  int32_T c1_i92;
  static char_T c1_cv5[37] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'm', 'p', 'o', 'w', 'e', 'r', '_', 'n', 'e', 'e', 'd',
    'C', 'o', 'm', 'p', 'l', 'e', 'x', 'I', 'n', 'p', 'u', 't' };

  char_T c1_u[37];
  const mxArray *c1_y = NULL;
  (void)chartInstance;
  for (c1_i92 = 0; c1_i92 < 37; c1_i92++) {
    c1_u[c1_i92] = c1_cv5[c1_i92];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 37), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c1_y));
}

static void c1_b_eig(SFc1_Model_01InstanceStruct *chartInstance, real_T c1_A[9],
                     creal_T c1_V[9], creal_T c1_D[9])
{
  int32_T c1_i93;
  static creal_T c1_dc5 = { 0.0, 0.0 };

  creal_T c1_b_A[9];
  int32_T c1_i94;
  creal_T c1_c_A[9];
  creal_T c1_beta1[3];
  creal_T c1_alpha1[3];
  real_T c1_info;
  real_T c1_b_info;
  int32_T c1_coltop;
  int32_T c1_b_coltop;
  int32_T c1_ix0;
  int32_T c1_b_ix0;
  int32_T c1_c_ix0;
  real_T c1_colnorm;
  real_T c1_scale;
  int32_T c1_kstart;
  int32_T c1_a;
  int32_T c1_kend;
  int32_T c1_b_kstart;
  int32_T c1_b_kend;
  int32_T c1_b_a;
  int32_T c1_b;
  int32_T c1_c_a;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_absxk;
  real_T c1_t;
  real_T c1_c_x;
  real_T c1_d_x;
  int32_T c1_c_coltop;
  int32_T c1_d_a;
  int32_T c1_e_a;
  int32_T c1_i95;
  int32_T c1_f_a;
  int32_T c1_c_b;
  int32_T c1_g_a;
  int32_T c1_d_b;
  boolean_T c1_b_overflow;
  int32_T c1_j;
  int32_T c1_b_j;
  creal_T c1_d_A;
  real_T c1_B;
  real_T c1_y;
  real_T c1_c_info;
  real_T c1_d_info;
  int32_T c1_i96;
  creal_T c1_b_alpha1[3];
  int32_T c1_i97;
  creal_T c1_b_beta1[3];
  int32_T c1_i98;
  int32_T c1_c_j;
  int32_T c1_d_j;
  for (c1_i93 = 0; c1_i93 < 9; c1_i93++) {
    c1_b_A[c1_i93].re = c1_A[c1_i93] + c1_dc5.re;
    c1_b_A[c1_i93].im = c1_dc5.im;
  }

  for (c1_i94 = 0; c1_i94 < 9; c1_i94++) {
    c1_c_A[c1_i94] = c1_b_A[c1_i94];
  }

  c1_eml_matlab_zggev(chartInstance, c1_c_A, &c1_info, c1_alpha1, c1_beta1, c1_V);
  c1_b_info = c1_info;
  for (c1_coltop = 1; c1_coltop < 8; c1_coltop += 3) {
    c1_b_coltop = c1_coltop;
    c1_ix0 = c1_b_coltop;
    c1_b_ix0 = c1_ix0;
    c1_b_eml_switch_helper(chartInstance);
    c1_c_ix0 = c1_b_ix0;
    c1_colnorm = 0.0;
    c1_realmin(chartInstance);
    c1_scale = 2.2250738585072014E-308;
    c1_kstart = c1_c_ix0;
    c1_a = c1_kstart;
    c1_kend = c1_a;
    c1_b_kstart = c1_kstart;
    c1_b_kend = c1_kend + 2;
    c1_b_a = c1_b_kstart;
    c1_b = c1_b_kend;
    c1_c_a = c1_b_a;
    c1_b_b = c1_b;
    if (c1_c_a > c1_b_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_b_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_k = c1_b_kstart; c1_k <= c1_b_kend; c1_k++) {
      c1_b_k = c1_k;
      c1_x = c1_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_k), 1, 9, 1, 0) - 1].re;
      c1_b_x = c1_x;
      c1_absxk = muDoubleScalarAbs(c1_b_x);
      if (c1_absxk > c1_scale) {
        c1_t = c1_scale / c1_absxk;
        c1_colnorm = 1.0 + c1_colnorm * c1_t * c1_t;
        c1_scale = c1_absxk;
      } else {
        c1_t = c1_absxk / c1_scale;
        c1_colnorm += c1_t * c1_t;
      }

      c1_c_x = c1_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c1_b_k), 1, 9, 1, 0) - 1].im;
      c1_d_x = c1_c_x;
      c1_absxk = muDoubleScalarAbs(c1_d_x);
      if (c1_absxk > c1_scale) {
        c1_t = c1_scale / c1_absxk;
        c1_colnorm = 1.0 + c1_colnorm * c1_t * c1_t;
        c1_scale = c1_absxk;
      } else {
        c1_t = c1_absxk / c1_scale;
        c1_colnorm += c1_t * c1_t;
      }
    }

    c1_colnorm = c1_scale * muDoubleScalarSqrt(c1_colnorm);
    c1_c_coltop = c1_b_coltop;
    c1_d_a = c1_b_coltop;
    c1_e_a = c1_d_a + 2;
    c1_i95 = c1_e_a;
    c1_f_a = c1_c_coltop;
    c1_c_b = c1_i95;
    c1_g_a = c1_f_a;
    c1_d_b = c1_c_b;
    if (c1_g_a > c1_d_b) {
      c1_b_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_b_overflow = (c1_d_b > 2147483646);
    }

    if (c1_b_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
    }

    for (c1_j = c1_c_coltop; c1_j <= c1_i95; c1_j++) {
      c1_b_j = c1_j;
      c1_d_A.re = c1_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 9, 1, 0) - 1].re;
      c1_d_A.im = c1_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 9, 1, 0) - 1].im;
      c1_B = c1_colnorm;
      c1_y = c1_B;
      c1_d_A = c1_eml_div(chartInstance, c1_d_A, c1_y);
      c1_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_j), 1, 9, 1, 0) - 1].re = c1_d_A.re;
      c1_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_j), 1, 9, 1, 0) - 1].im = c1_d_A.im;
    }
  }

  c1_c_info = c1_b_info;
  c1_d_info = c1_c_info;
  for (c1_i96 = 0; c1_i96 < 3; c1_i96++) {
    c1_b_alpha1[c1_i96] = c1_alpha1[c1_i96];
  }

  for (c1_i97 = 0; c1_i97 < 3; c1_i97++) {
    c1_b_beta1[c1_i97] = c1_beta1[c1_i97];
  }

  c1_b_eml_div(chartInstance, c1_b_alpha1, c1_b_beta1, c1_alpha1);
  for (c1_i98 = 0; c1_i98 < 9; c1_i98++) {
    c1_D[c1_i98] = c1_dc5;
  }

  for (c1_c_j = 1; c1_c_j < 4; c1_c_j++) {
    c1_d_j = c1_c_j;
    c1_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_d_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
      1].re = c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 1, 0) - 1].re;
    c1_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_d_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
      1].im = c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 1, 0) - 1].im;
  }

  if (c1_d_info < 0.0) {
    c1_eml_warning(chartInstance);
  } else {
    if (c1_d_info > 0.0) {
      c1_b_eml_warning(chartInstance);
    }
  }
}

static void c1_eml_matlab_zggev(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], real_T *c1_info, creal_T c1_alpha1[3], creal_T c1_beta1[3],
  creal_T c1_V[9])
{
  int32_T c1_i99;
  creal_T c1_b_A[9];
  real_T c1_anrm;
  int32_T c1_i100;
  int32_T c1_i101;
  int32_T c1_i102;
  boolean_T c1_ilascl;
  real_T c1_anrmto;
  int32_T c1_rscale[3];
  int32_T c1_ihi;
  int32_T c1_ilo;
  int32_T c1_b_ilo;
  int32_T c1_b_ihi;
  real_T c1_b_info;
  int32_T c1_i103;
  creal_T c1_c_A[9];
  int32_T c1_c_ilo;
  int32_T c1_c_ihi;
  int32_T c1_a;
  int32_T c1_b_a;
  int32_T c1_i;
  int32_T c1_k;
  int32_T c1_j;
  int32_T c1_b_j;
  creal_T c1_tmp;
  int32_T c1_c_a;
  int32_T c1_d_a;
  int32_T c1_e_a;
  int32_T c1_f_a;
  int32_T c1_i104;
  int32_T c1_g_a;
  int32_T c1_h_a;
  boolean_T c1_overflow;
  int32_T c1_b_i;
  int32_T c1_c_j;
  int32_T c1_jc;
  real_T c1_b_jc;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_y;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_b_y;
  real_T c1_vtemp;
  int32_T c1_jr;
  real_T c1_b_jr;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_c_y;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_d_y;
  real_T c1_e_y;
  real_T c1_f_y;
  int32_T c1_c_jr;
  real_T c1_b;
  boolean_T guard1 = false;
  *c1_info = 0.0;
  c1_realmin(chartInstance);
  c1_eps(chartInstance);
  for (c1_i99 = 0; c1_i99 < 9; c1_i99++) {
    c1_b_A[c1_i99] = c1_A[c1_i99];
  }

  c1_anrm = c1_eml_matlab_zlangeM(chartInstance, c1_b_A);
  if (!c1_isfinite(chartInstance, c1_anrm)) {
    for (c1_i100 = 0; c1_i100 < 3; c1_i100++) {
      c1_alpha1[c1_i100].re = rtNaN;
      c1_alpha1[c1_i100].im = 0.0;
    }

    for (c1_i101 = 0; c1_i101 < 3; c1_i101++) {
      c1_beta1[c1_i101].re = rtNaN;
      c1_beta1[c1_i101].im = 0.0;
    }

    for (c1_i102 = 0; c1_i102 < 9; c1_i102++) {
      c1_V[c1_i102].re = rtNaN;
      c1_V[c1_i102].im = 0.0;
    }
  } else {
    c1_ilascl = false;
    c1_anrmto = c1_anrm;
    guard1 = false;
    if (c1_anrm > 0.0) {
      if (c1_anrm < 6.7178761075670888E-139) {
        c1_anrmto = 6.7178761075670888E-139;
        c1_ilascl = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      if (c1_anrm > 1.4885657073574029E+138) {
        c1_anrmto = 1.4885657073574029E+138;
        c1_ilascl = true;
      }
    }

    if (c1_ilascl) {
      c1_c_eml_matlab_zlascl(chartInstance, c1_anrm, c1_anrmto, c1_A);
    }

    c1_b_eml_matlab_zggbal(chartInstance, c1_A, &c1_ilo, &c1_ihi, c1_rscale);
    c1_b_ilo = c1_ilo;
    c1_b_ihi = c1_ihi;
    c1_b_eml_matlab_zgghrd(chartInstance, c1_b_ilo, c1_b_ihi, c1_A, c1_V);
    c1_c_eml_matlab_zhgeqz(chartInstance, c1_A, c1_b_ilo, c1_b_ihi, c1_V,
      &c1_b_info, c1_alpha1, c1_beta1);
    *c1_info = c1_b_info;
    if (*c1_info != 0.0) {
    } else {
      for (c1_i103 = 0; c1_i103 < 9; c1_i103++) {
        c1_c_A[c1_i103] = c1_A[c1_i103];
      }

      c1_b_eml_matlab_ztgevc(chartInstance, c1_c_A, c1_V);
      c1_c_ilo = c1_b_ilo;
      c1_c_ihi = c1_b_ihi;
      if (c1_c_ilo > 1) {
        c1_a = c1_c_ilo;
        c1_b_a = c1_a - 1;
        c1_i = c1_b_a;
        while (c1_i >= 1) {
          c1_k = c1_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_i), 1, 3, 1, 0) - 1];
          if (c1_k != c1_i) {
            for (c1_j = 1; c1_j < 4; c1_j++) {
              c1_b_j = c1_j;
              c1_tmp.re = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].
                re;
              c1_tmp.im = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].
                im;
              c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re = c1_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
              c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im = c1_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
              c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re =
                c1_tmp.re;
              c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im =
                c1_tmp.im;
            }
          }

          c1_c_a = c1_i;
          c1_d_a = c1_c_a - 1;
          c1_i = c1_d_a;
        }
      }

      if (c1_c_ihi < 3) {
        c1_e_a = c1_c_ihi;
        c1_f_a = c1_e_a + 1;
        c1_i104 = c1_f_a;
        c1_g_a = c1_i104;
        c1_h_a = c1_g_a;
        if (c1_h_a > 3) {
          c1_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_overflow = false;
        }

        if (c1_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_overflow);
        }

        for (c1_b_i = c1_i104; c1_b_i < 4; c1_b_i++) {
          c1_i = c1_b_i;
          c1_k = c1_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_i), 1, 3, 1, 0) - 1];
          if (c1_k != c1_i) {
            for (c1_c_j = 1; c1_c_j < 4; c1_c_j++) {
              c1_b_j = c1_c_j;
              c1_tmp.re = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].
                re;
              c1_tmp.im = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].
                im;
              c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re = c1_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
              c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im = c1_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
              c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re =
                c1_tmp.re;
              c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im =
                c1_tmp.im;
            }
          }
        }
      }

      for (c1_jc = 0; c1_jc < 3; c1_jc++) {
        c1_b_jc = 1.0 + (real_T)c1_jc;
        c1_tmp.re = c1_V[3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1)].re;
        c1_tmp.im = c1_V[3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1)].im;
        c1_x = c1_tmp.re;
        c1_b_x = c1_x;
        c1_y = muDoubleScalarAbs(c1_b_x);
        c1_c_x = c1_tmp.im;
        c1_d_x = c1_c_x;
        c1_b_y = muDoubleScalarAbs(c1_d_x);
        c1_vtemp = c1_y + c1_b_y;
        for (c1_jr = 0; c1_jr < 2; c1_jr++) {
          c1_b_jr = 2.0 + (real_T)c1_jr;
          c1_tmp.re = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1)) - 1].re;
          c1_tmp.im = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1)) - 1].im;
          c1_e_x = c1_tmp.re;
          c1_f_x = c1_e_x;
          c1_c_y = muDoubleScalarAbs(c1_f_x);
          c1_g_x = c1_tmp.im;
          c1_h_x = c1_g_x;
          c1_d_y = muDoubleScalarAbs(c1_h_x);
          c1_e_y = c1_c_y + c1_d_y;
          c1_f_y = c1_e_y;
          if (c1_f_y > c1_vtemp) {
            c1_vtemp = c1_f_y;
          }
        }

        if (c1_vtemp >= 6.7178761075670888E-139) {
          c1_vtemp = 1.0 / c1_vtemp;
          for (c1_c_jr = 0; c1_c_jr < 3; c1_c_jr++) {
            c1_b_jr = 1.0 + (real_T)c1_c_jr;
            c1_tmp.re = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1)) - 1].re;
            c1_tmp.im = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1)) - 1].im;
            c1_b = c1_vtemp;
            c1_tmp.re *= c1_b;
            c1_tmp.im *= c1_b;
            c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    c1_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1))
              - 1].re = c1_tmp.re;
            c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    c1_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1))
              - 1].im = c1_tmp.im;
          }
        }
      }

      if (c1_ilascl) {
        c1_d_eml_matlab_zlascl(chartInstance, c1_anrmto, c1_anrm, c1_alpha1);
      }
    }
  }
}

static void c1_eml_matlab_zgghrd(SFc1_Model_01InstanceStruct *chartInstance,
  int32_T c1_ilo, int32_T c1_ihi, creal_T c1_A[9], creal_T c1_b_A[9], creal_T
  c1_Z[9])
{
  int32_T c1_i105;
  for (c1_i105 = 0; c1_i105 < 9; c1_i105++) {
    c1_b_A[c1_i105] = c1_A[c1_i105];
  }

  c1_b_eml_matlab_zgghrd(chartInstance, c1_ilo, c1_ihi, c1_b_A, c1_Z);
}

static void c1_b_eml_matlab_zhgeqz(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T c1_ilo, int32_T c1_ihi, creal_T c1_Z[9], real_T
  *c1_info, creal_T c1_alpha1[3], creal_T c1_beta1[3], creal_T c1_b_A[9],
  creal_T c1_b_Z[9])
{
  int32_T c1_i106;
  int32_T c1_i107;
  for (c1_i106 = 0; c1_i106 < 9; c1_i106++) {
    c1_b_A[c1_i106] = c1_A[c1_i106];
  }

  for (c1_i107 = 0; c1_i107 < 9; c1_i107++) {
    c1_b_Z[c1_i107] = c1_Z[c1_i107];
  }

  c1_c_eml_matlab_zhgeqz(chartInstance, c1_b_A, c1_ilo, c1_ihi, c1_b_Z, c1_info,
    c1_alpha1, c1_beta1);
}

static void c1_eml_matlab_ztgevc(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], creal_T c1_V[9], creal_T c1_b_V[9])
{
  int32_T c1_i108;
  int32_T c1_i109;
  creal_T c1_b_A[9];
  for (c1_i108 = 0; c1_i108 < 9; c1_i108++) {
    c1_b_V[c1_i108] = c1_V[c1_i108];
  }

  for (c1_i109 = 0; c1_i109 < 9; c1_i109++) {
    c1_b_A[c1_i109] = c1_A[c1_i109];
  }

  c1_b_eml_matlab_ztgevc(chartInstance, c1_b_A, c1_b_V);
}

static creal_T c1_rdivide(SFc1_Model_01InstanceStruct *chartInstance, creal_T
  c1_x, creal_T c1_y)
{
  creal_T c1_z;
  real_T c1_ar;
  real_T c1_ai;
  real_T c1_br;
  real_T c1_bi;
  real_T c1_brm;
  real_T c1_bim;
  real_T c1_s;
  real_T c1_d;
  real_T c1_nr;
  real_T c1_ni;
  real_T c1_sgnbr;
  real_T c1_sgnbi;
  (void)chartInstance;
  c1_ar = c1_x.re;
  c1_ai = c1_x.im;
  c1_br = c1_y.re;
  c1_bi = c1_y.im;
  if (c1_bi == 0.0) {
    if (c1_ai == 0.0) {
      c1_z.re = c1_ar / c1_br;
      c1_z.im = 0.0;
    } else if (c1_ar == 0.0) {
      c1_z.re = 0.0;
      c1_z.im = c1_ai / c1_br;
    } else {
      c1_z.re = c1_ar / c1_br;
      c1_z.im = c1_ai / c1_br;
    }
  } else if (c1_br == 0.0) {
    if (c1_ar == 0.0) {
      c1_z.re = c1_ai / c1_bi;
      c1_z.im = 0.0;
    } else if (c1_ai == 0.0) {
      c1_z.re = 0.0;
      c1_z.im = -(c1_ar / c1_bi);
    } else {
      c1_z.re = c1_ai / c1_bi;
      c1_z.im = -(c1_ar / c1_bi);
    }
  } else {
    c1_brm = muDoubleScalarAbs(c1_br);
    c1_bim = muDoubleScalarAbs(c1_bi);
    if (c1_brm > c1_bim) {
      c1_s = c1_bi / c1_br;
      c1_d = c1_br + c1_s * c1_bi;
      c1_nr = c1_ar + c1_s * c1_ai;
      c1_ni = c1_ai - c1_s * c1_ar;
      c1_z.re = c1_nr / c1_d;
      c1_z.im = c1_ni / c1_d;
    } else if (c1_bim == c1_brm) {
      if (c1_br > 0.0) {
        c1_sgnbr = 0.5;
      } else {
        c1_sgnbr = -0.5;
      }

      if (c1_bi > 0.0) {
        c1_sgnbi = 0.5;
      } else {
        c1_sgnbi = -0.5;
      }

      c1_nr = c1_ar * c1_sgnbr + c1_ai * c1_sgnbi;
      c1_ni = c1_ai * c1_sgnbr - c1_ar * c1_sgnbi;
      c1_z.re = c1_nr / c1_brm;
      c1_z.im = c1_ni / c1_brm;
    } else {
      c1_s = c1_br / c1_bi;
      c1_d = c1_bi + c1_s * c1_br;
      c1_nr = c1_s * c1_ar + c1_ai;
      c1_ni = c1_s * c1_ai - c1_ar;
      c1_z.re = c1_nr / c1_d;
      c1_z.im = c1_ni / c1_d;
    }
  }

  return c1_z;
}

static creal_T c1_power(SFc1_Model_01InstanceStruct *chartInstance, creal_T c1_a,
  real_T c1_b)
{
  creal_T c1_y;
  real_T c1_b_b;
  real_T c1_bk;
  real_T c1_c_b;
  real_T c1_ar;
  real_T c1_ai;
  real_T c1_br;
  real_T c1_bi;
  real_T c1_ytmp;
  real_T c1_x;
  real_T c1_xk;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_r;
  real_T c1_d6;
  int8_T c1_i110;
  int8_T c1_b_r;
  creal_T c1_t;
  real_T c1_e_x;
  boolean_T c1_d_b;
  real_T c1_A;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_b_y;
  real_T c1_b_A;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_c_y;
  real_T c1_l_x;
  real_T c1_d_y;
  real_T c1_x1;
  real_T c1_x2;
  real_T c1_b_a;
  real_T c1_e_b;
  real_T c1_z;
  real_T c1_e_y;
  real_T c1_m_x;
  real_T c1_c_r;
  real_T c1_b_x1;
  real_T c1_b_x2;
  real_T c1_c_a;
  real_T c1_f_b;
  real_T c1_f_y;
  real_T c1_g_y;
  real_T c1_n_x;
  real_T c1_d_r;
  real_T c1_d_a;
  real_T c1_tr;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T guard5 = false;
  boolean_T guard6 = false;
  boolean_T guard7 = false;
  c1_b_b = c1_b;
  c1_scalarEg(chartInstance);
  c1_bk = c1_b_b;
  c1_c_b = c1_bk;
  c1_scalarEg(chartInstance);
  c1_ar = c1_a.re;
  c1_ai = c1_a.im;
  c1_br = c1_c_b;
  c1_bi = 0.0;
  guard1 = false;
  guard2 = false;
  if (c1_ai == 0.0) {
    if (c1_bi == 0.0) {
      if (c1_ar >= 0.0) {
        c1_y.re = muDoubleScalarPower(c1_ar, c1_br);
        c1_y.im = 0.0;
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
    if (c1_ar == 0.0) {
      if (c1_bi == 0.0) {
        if (muDoubleScalarFloor(c1_br) == c1_br) {
          if (muDoubleScalarAbs(c1_br) <= 9.007199254740992E+15) {
            c1_ytmp = muDoubleScalarPower(c1_ai, c1_br);
            c1_x = c1_br;
            c1_eml_scalar_eg(chartInstance);
            c1_xk = c1_x;
            c1_b_x = c1_xk;
            c1_eml_scalar_eg(chartInstance);
            c1_c_x = c1_b_x / 4.0;
            c1_d_x = c1_c_x;
            c1_d_x = muDoubleScalarFloor(c1_d_x);
            c1_r = c1_b_x - c1_d_x * 4.0;
            c1_d6 = muDoubleScalarRound(c1_r);
            if (c1_d6 < 128.0) {
              if (c1_d6 >= -128.0) {
                c1_i110 = (int8_T)c1_d6;
              } else {
                c1_i110 = MIN_int8_T;
              }
            } else if (c1_d6 >= 128.0) {
              c1_i110 = MAX_int8_T;
            } else {
              c1_i110 = 0;
            }

            c1_b_r = c1_i110;
            if ((real_T)c1_b_r == 3.0) {
              c1_y.re = 0.0;
              c1_y.im = -c1_ytmp;
            } else if ((real_T)c1_b_r == 2.0) {
              c1_y.re = -c1_ytmp;
              c1_y.im = 0.0;
            } else if ((real_T)c1_b_r == 1.0) {
              c1_y.re = 0.0;
              c1_y.im = c1_ytmp;
            } else {
              c1_y.re = c1_ytmp;
              c1_y.im = 0.0;
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
      c1_t = c1_a;
      c1_realmax(chartInstance);
      guard6 = false;
      if (c1_t.im == 0.0) {
        c1_e_x = c1_t.re;
        c1_d_b = muDoubleScalarIsNaN(c1_e_x);
        if (c1_d_b) {
        } else {
          guard6 = true;
        }
      } else {
        guard6 = true;
      }

      if (guard6 == true) {
        guard7 = false;
        if (muDoubleScalarAbs(c1_t.re) > 8.9884656743115785E+307) {
          guard7 = true;
        } else if (muDoubleScalarAbs(c1_t.im) > 8.9884656743115785E+307) {
          guard7 = true;
        } else {
          c1_b_x1 = c1_t.re;
          c1_b_x2 = c1_t.im;
          c1_c_a = c1_b_x1;
          c1_f_b = c1_b_x2;
          c1_f_y = muDoubleScalarHypot(c1_c_a, c1_f_b);
          c1_g_y = c1_t.im;
          c1_n_x = c1_t.re;
          c1_d_r = muDoubleScalarAtan2(c1_g_y, c1_n_x);
          c1_t.re = muDoubleScalarLog(c1_f_y);
          c1_t.im = c1_d_r;
        }

        if (guard7 == true) {
          c1_A = c1_t.re;
          c1_f_x = c1_A;
          c1_g_x = c1_f_x;
          c1_h_x = c1_g_x;
          c1_b_y = c1_h_x / 2.0;
          c1_b_A = c1_t.im;
          c1_i_x = c1_b_A;
          c1_j_x = c1_i_x;
          c1_k_x = c1_j_x;
          c1_c_y = c1_k_x / 2.0;
          c1_l_x = c1_b_y;
          c1_d_y = c1_c_y;
          c1_eml_scalar_eg(chartInstance);
          c1_x1 = c1_l_x;
          c1_x2 = c1_d_y;
          c1_b_a = c1_x1;
          c1_e_b = c1_x2;
          c1_z = muDoubleScalarHypot(c1_b_a, c1_e_b);
          c1_e_y = c1_t.im;
          c1_m_x = c1_t.re;
          c1_c_r = muDoubleScalarAtan2(c1_e_y, c1_m_x);
          c1_t.re = muDoubleScalarLog(c1_z) + 0.69314718055994529;
          c1_t.im = c1_c_r;
        }
      }

      c1_d_a = c1_c_b;
      c1_t.re *= c1_d_a;
      c1_t.im *= c1_d_a;
      c1_tr = muDoubleScalarExp(c1_t.re);
      c1_y.re = c1_tr * muDoubleScalarCos(c1_t.im);
      c1_y.im = c1_tr * muDoubleScalarSin(c1_t.im);
    }
  }

  return c1_y;
}

static void c1_eml_lusolve(SFc1_Model_01InstanceStruct *chartInstance, creal_T
  c1_A[9], creal_T c1_B[9], creal_T c1_X[9])
{
  int32_T c1_i111;
  creal_T c1_b_A[9];
  int32_T c1_r1;
  int32_T c1_r2;
  int32_T c1_r3;
  creal_T c1_x;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_y;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_b_y;
  real_T c1_maxval;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_c_y;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_d_y;
  real_T c1_a21;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_e_y;
  real_T c1_l_x;
  real_T c1_m_x;
  real_T c1_f_y;
  real_T c1_d;
  creal_T c1_c_A;
  creal_T c1_d_A;
  creal_T c1_e_A;
  creal_T c1_f_A;
  creal_T c1_g_A;
  creal_T c1_h_A;
  creal_T c1_i_A;
  creal_T c1_j_A;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_g_y;
  real_T c1_p_x;
  real_T c1_q_x;
  real_T c1_h_y;
  real_T c1_b_d;
  real_T c1_r_x;
  real_T c1_s_x;
  real_T c1_i_y;
  real_T c1_t_x;
  real_T c1_u_x;
  real_T c1_j_y;
  real_T c1_c_d;
  int32_T c1_rtemp;
  creal_T c1_k_A;
  creal_T c1_l_A;
  creal_T c1_m_A;
  static creal_T c1_dc6 = { 0.0, 0.0 };

  boolean_T c1_n_A;
  boolean_T c1_o_A;
  boolean_T c1_p_A;
  int32_T c1_k;
  int32_T c1_b_k;
  creal_T c1_b_B;
  creal_T c1_q_A;
  creal_T c1_b_X;
  creal_T c1_c_X;
  creal_T c1_d_X;
  creal_T c1_r_A;
  creal_T c1_e_X;
  creal_T c1_f_X;
  creal_T c1_s_A;
  creal_T c1_g_X;
  creal_T c1_h_X;
  creal_T c1_i_X;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  for (c1_i111 = 0; c1_i111 < 9; c1_i111++) {
    c1_b_A[c1_i111] = c1_A[c1_i111];
  }

  c1_r1 = 1;
  c1_r2 = 2;
  c1_r3 = 3;
  c1_x = c1_b_A[0];
  c1_b_x = c1_x.re;
  c1_c_x = c1_b_x;
  c1_y = muDoubleScalarAbs(c1_c_x);
  c1_d_x = c1_x.im;
  c1_e_x = c1_d_x;
  c1_b_y = muDoubleScalarAbs(c1_e_x);
  c1_maxval = c1_y + c1_b_y;
  c1_x = c1_b_A[1];
  c1_f_x = c1_x.re;
  c1_g_x = c1_f_x;
  c1_c_y = muDoubleScalarAbs(c1_g_x);
  c1_h_x = c1_x.im;
  c1_i_x = c1_h_x;
  c1_d_y = muDoubleScalarAbs(c1_i_x);
  c1_a21 = c1_c_y + c1_d_y;
  if (c1_a21 > c1_maxval) {
    c1_maxval = c1_a21;
    c1_r1 = 2;
    c1_r2 = 1;
  }

  c1_x = c1_b_A[2];
  c1_j_x = c1_x.re;
  c1_k_x = c1_j_x;
  c1_e_y = muDoubleScalarAbs(c1_k_x);
  c1_l_x = c1_x.im;
  c1_m_x = c1_l_x;
  c1_f_y = muDoubleScalarAbs(c1_m_x);
  c1_d = c1_e_y + c1_f_y;
  if (c1_d > c1_maxval) {
    c1_r1 = 3;
    c1_r2 = 2;
    c1_r3 = 1;
  }

  c1_c_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r2), 1, 3, 1, 0) - 1].re;
  c1_c_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r2), 1, 3, 1, 0) - 1].im;
  c1_d_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r1), 1, 3, 1, 0) - 1].re;
  c1_d_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r1), 1, 3, 1, 0) - 1].im;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r2), 1, 3, 1, 0) - 1] = c1_rdivide(chartInstance, c1_c_A, c1_d_A);
  c1_e_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) - 1].re;
  c1_e_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) - 1].im;
  c1_f_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r1), 1, 3, 1, 0) - 1].re;
  c1_f_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r1), 1, 3, 1, 0) - 1].im;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) - 1] = c1_rdivide(chartInstance, c1_e_A, c1_f_A);
  c1_g_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r2), 1, 3, 1, 0) - 1].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 2].re - c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) - 1].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 2].im;
  c1_g_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r2), 1, 3, 1, 0) - 1].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 2].im + c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) - 1].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 2].re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r2), 1, 3, 1, 0) + 2].re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 2].re -
    c1_g_A.re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r2), 1, 3, 1, 0) + 2].im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 2].im -
    c1_g_A.im;
  c1_h_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) - 1].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 2].re - c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) - 1].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 2].im;
  c1_h_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) - 1].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 2].im + c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) - 1].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 2].re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 2].re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 2].re -
    c1_h_A.re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 2].im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 2].im -
    c1_h_A.im;
  c1_i_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r2), 1, 3, 1, 0) - 1].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 5].re - c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) - 1].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 5].im;
  c1_i_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r2), 1, 3, 1, 0) - 1].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 5].im + c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) - 1].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 5].re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r2), 1, 3, 1, 0) + 5].re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 5].re -
    c1_i_A.re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r2), 1, 3, 1, 0) + 5].im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 5].im -
    c1_i_A.im;
  c1_j_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) - 1].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 5].re - c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) - 1].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 5].im;
  c1_j_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) - 1].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 5].im + c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) - 1].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) + 5].re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 5].re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5].re -
    c1_j_A.re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 5].im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5].im -
    c1_j_A.im;
  c1_x.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c1_r3), 1, 3, 1, 0) + 2].re;
  c1_x.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c1_r3), 1, 3, 1, 0) + 2].im;
  c1_n_x = c1_x.re;
  c1_o_x = c1_n_x;
  c1_g_y = muDoubleScalarAbs(c1_o_x);
  c1_p_x = c1_x.im;
  c1_q_x = c1_p_x;
  c1_h_y = muDoubleScalarAbs(c1_q_x);
  c1_b_d = c1_g_y + c1_h_y;
  c1_x.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c1_r2), 1, 3, 1, 0) + 2].re;
  c1_x.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c1_r2), 1, 3, 1, 0) + 2].im;
  c1_r_x = c1_x.re;
  c1_s_x = c1_r_x;
  c1_i_y = muDoubleScalarAbs(c1_s_x);
  c1_t_x = c1_x.im;
  c1_u_x = c1_t_x;
  c1_j_y = muDoubleScalarAbs(c1_u_x);
  c1_c_d = c1_i_y + c1_j_y;
  if (c1_b_d > c1_c_d) {
    c1_rtemp = c1_r2;
    c1_r2 = c1_r3;
    c1_r3 = c1_rtemp;
  }

  c1_k_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) + 2].re;
  c1_k_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) + 2].im;
  c1_l_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r2), 1, 3, 1, 0) + 2].re;
  c1_l_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r2), 1, 3, 1, 0) + 2].im;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 2] = c1_rdivide(chartInstance, c1_k_A, c1_l_A);
  c1_m_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) + 2].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) + 5].re - c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 2].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) + 5].im;
  c1_m_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c1_r3), 1, 3, 1, 0) + 2].re *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) + 5].im + c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 2].im *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) + 5].re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 5].re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5].re -
    c1_m_A.re;
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 5].im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5].im -
    c1_m_A.im;
  c1_n_A = ((c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c1_r1), 1, 3, 1, 0) - 1].re == c1_dc6.re) &&
            (c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c1_r1), 1, 3, 1, 0) - 1].im == c1_dc6.im));
  guard1 = false;
  guard2 = false;
  if (c1_n_A) {
    guard2 = true;
  } else {
    c1_o_A = ((c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 2].re ==
               c1_dc6.re) && (c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 2].im ==
               c1_dc6.im));
    if (c1_o_A) {
      guard2 = true;
    } else {
      c1_p_A = ((c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5].re ==
                 c1_dc6.re) && (c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5].im ==
                 c1_dc6.im));
      if (c1_p_A) {
        guard1 = true;
      }
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c1_c_eml_warning(chartInstance);
  }

  for (c1_k = 1; c1_k < 4; c1_k++) {
    c1_b_k = c1_k;
    c1_b_B.re = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", (real_T)c1_b_k), 1, 3, 1, 0) - 1].re;
    c1_b_B.im = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", (real_T)c1_b_k), 1, 3, 1, 0) - 1].im;
    c1_q_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) - 1].re;
    c1_q_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) - 1].im;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) -
      1] = c1_rdivide(chartInstance, c1_b_B, c1_q_A);
    c1_b_X.re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r1), 1, 3, 1, 0) + 2].re - c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r1), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 2].im;
    c1_b_X.im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r1), 1, 3, 1, 0) + 2].im + c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r1), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 2].re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) -
      1].re = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 2].re - c1_b_X.re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) -
      1].im = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 2].im - c1_b_X.im;
    c1_c_X.re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r1), 1, 3, 1, 0) + 5].re - c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r1), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 5].im;
    c1_c_X.im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r1), 1, 3, 1, 0) + 5].im + c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r1), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 5].re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) -
      1].re = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 5].re - c1_c_X.re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) -
      1].im = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 5].im - c1_c_X.im;
    c1_d_X.re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) - 1].re;
    c1_d_X.im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) - 1].im;
    c1_r_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 2].re;
    c1_r_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 2].im;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) -
      1] = c1_rdivide(chartInstance, c1_d_X, c1_r_A);
    c1_e_X.re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r2), 1, 3, 1, 0) + 5].re - c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r2), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 5].im;
    c1_e_X.im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r2), 1, 3, 1, 0) + 5].im + c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r2), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 5].re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) -
      1].re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) - 1].re
      - c1_e_X.re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) -
      1].im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) - 1].im
      - c1_e_X.im;
    c1_f_X.re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) - 1].re;
    c1_f_X.im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) - 1].im;
    c1_s_A.re = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5].re;
    c1_s_A.im = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5].im;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) -
      1] = c1_rdivide(chartInstance, c1_f_X, c1_s_A);
    c1_g_X.re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r3), 1, 3, 1, 0) + 2].re - c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r3), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 2].im;
    c1_g_X.im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r3), 1, 3, 1, 0) + 2].im + c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r3), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 2].re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) -
      1].re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) - 1].re
      - c1_g_X.re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) -
      1].im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) - 1].im
      - c1_g_X.im;
    c1_h_X.re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r3), 1, 3, 1, 0) - 1].re - c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r3), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) - 1].im;
    c1_h_X.im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r3), 1, 3, 1, 0) - 1].im + c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r3), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) - 1].re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) -
      1].re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) - 1].re
      - c1_h_X.re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) -
      1].im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) - 1].im
      - c1_h_X.im;
    c1_i_X.re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r2), 1, 3, 1, 0) - 1].re - c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r2), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) - 1].im;
    c1_i_X.im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_r2), 1, 3, 1, 0) - 1].im + c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_r2), 1, 3, 2, 0) - 1)) - 1].im * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) - 1].re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) -
      1].re = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) - 1].re
      - c1_i_X.re;
    c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) -
      1].im = c1_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 2, 0) - 1)) - 1].im
      - c1_i_X.im;
  }
}

static void c1_e_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_sprintf, const char_T *c1_identifier, char_T c1_y[14])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_sprintf), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_sprintf);
}

static void c1_f_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, char_T c1_y[14])
{
  char_T c1_cv6[14];
  int32_T c1_i112;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_cv6, 1, 10, 0U, 1, 0U, 2, 1,
                14);
  for (c1_i112 = 0; c1_i112 < 14; c1_i112++) {
    c1_y[c1_i112] = c1_cv6[c1_i112];
  }

  sf_mex_destroy(&c1_u);
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc1_Model_01InstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(int32_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static int32_T c1_g_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_y;
  int32_T c1_i113;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_i113, 1, 6, 0U, 0, 0U, 0);
  c1_y = c1_i113;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_y;
  SFc1_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc1_Model_01InstanceStruct *)chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_h_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_Model_01, const char_T *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_i_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_Model_01), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_Model_01);
  return c1_y;
}

static uint8_T c1_i_emlrt_marshallIn(SFc1_Model_01InstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_c_eml_matlab_zlascl(SFc1_Model_01InstanceStruct *chartInstance,
  real_T c1_cfrom, real_T c1_cto, creal_T c1_A[9])
{
  real_T c1_cfromc;
  real_T c1_ctoc;
  boolean_T c1_notdone;
  real_T c1_cfrom1;
  real_T c1_cto1;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_y;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_b_y;
  real_T c1_mul;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_c_y;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_d_y;
  real_T c1_a;
  int32_T c1_i114;
  int32_T c1_i115;
  int32_T c1_i116;
  boolean_T guard1 = false;
  c1_realmin(chartInstance);
  c1_eps(chartInstance);
  c1_cfromc = c1_cfrom;
  c1_ctoc = c1_cto;
  c1_notdone = true;
  while (c1_notdone) {
    c1_cfrom1 = c1_cfromc * 2.0041683600089728E-292;
    c1_cto1 = c1_ctoc / 4.9896007738368E+291;
    c1_x = c1_cfrom1;
    c1_b_x = c1_x;
    c1_y = muDoubleScalarAbs(c1_b_x);
    c1_c_x = c1_ctoc;
    c1_d_x = c1_c_x;
    c1_b_y = muDoubleScalarAbs(c1_d_x);
    guard1 = false;
    if (c1_y > c1_b_y) {
      if (c1_ctoc != 0.0) {
        c1_mul = 2.0041683600089728E-292;
        c1_notdone = true;
        c1_cfromc = c1_cfrom1;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      c1_e_x = c1_cto1;
      c1_f_x = c1_e_x;
      c1_c_y = muDoubleScalarAbs(c1_f_x);
      c1_g_x = c1_cfromc;
      c1_h_x = c1_g_x;
      c1_d_y = muDoubleScalarAbs(c1_h_x);
      if (c1_c_y > c1_d_y) {
        c1_mul = 4.9896007738368E+291;
        c1_notdone = true;
        c1_ctoc = c1_cto1;
      } else {
        c1_mul = c1_ctoc / c1_cfromc;
        c1_notdone = false;
      }
    }

    c1_a = c1_mul;
    c1_i114 = 0;
    for (c1_i115 = 0; c1_i115 < 3; c1_i115++) {
      for (c1_i116 = 0; c1_i116 < 3; c1_i116++) {
        c1_A[c1_i116 + c1_i114].re *= c1_a;
        c1_A[c1_i116 + c1_i114].im *= c1_a;
      }

      c1_i114 += 3;
    }
  }
}

static void c1_b_eml_matlab_zggbal(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T *c1_ilo, int32_T *c1_ihi, int32_T c1_rscale[3])
{
  int32_T c1_i117;
  int32_T c1_b_ihi;
  int32_T c1_i;
  int32_T c1_j;
  boolean_T c1_found;
  int32_T c1_ii;
  real_T c1_nzcount;
  int32_T c1_c_ihi;
  int32_T c1_b;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_jj;
  int32_T c1_b_jj;
  static creal_T c1_dc7 = { 0.0, 0.0 };

  boolean_T c1_b_A;
  int32_T c1_a;
  int32_T c1_b_a;
  int32_T c1_b_i;
  int32_T c1_b_j;
  boolean_T c1_b_found;
  int32_T c1_b_ilo;
  int32_T c1_d_ihi;
  int32_T c1_c_i;
  int32_T c1_c_j;
  boolean_T c1_c_found;
  int32_T c1_c_ilo;
  int32_T c1_e_ihi;
  int32_T c1_c_a;
  int32_T c1_c_b;
  int32_T c1_d_a;
  int32_T c1_d_b;
  boolean_T c1_b_overflow;
  int32_T c1_c_jj;
  int32_T c1_d_jj;
  real_T c1_b_nzcount;
  int32_T c1_d_ilo;
  int32_T c1_f_ihi;
  int32_T c1_e_a;
  int32_T c1_e_b;
  int32_T c1_f_a;
  int32_T c1_f_b;
  boolean_T c1_c_overflow;
  int32_T c1_b_ii;
  int32_T c1_c_ii;
  boolean_T c1_c_A;
  int32_T c1_m;
  int32_T c1_d_i;
  int32_T c1_d_j;
  int32_T c1_e_ilo;
  int32_T c1_g_ihi;
  int32_T c1_f_ilo;
  int32_T c1_g_a;
  int32_T c1_h_a;
  boolean_T c1_d_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  creal_T c1_atmp;
  int32_T c1_h_ihi;
  int32_T c1_g_b;
  int32_T c1_h_b;
  boolean_T c1_e_overflow;
  int32_T c1_c_k;
  int32_T c1_i_a;
  int32_T c1_j_a;
  int32_T c1_b_m;
  int32_T c1_e_i;
  int32_T c1_e_j;
  int32_T c1_i_ihi;
  int32_T c1_d_k;
  int32_T c1_e_k;
  int32_T c1_j_ihi;
  int32_T c1_i_b;
  int32_T c1_j_b;
  boolean_T c1_f_overflow;
  int32_T c1_f_k;
  int32_T c1_k_a;
  int32_T c1_l_a;
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
  for (c1_i117 = 0; c1_i117 < 3; c1_i117++) {
    c1_rscale[c1_i117] = 0;
  }

  *c1_ilo = 1;
  *c1_ihi = 3;
  do {
    exitg2 = 0;
    c1_b_ihi = *c1_ihi;
    c1_i = 0;
    c1_j = 0;
    c1_found = false;
    c1_ii = c1_b_ihi;
    exitg5 = false;
    while ((exitg5 == false) && (c1_ii > 0)) {
      c1_nzcount = 0.0;
      c1_i = c1_ii;
      c1_j = c1_b_ihi;
      c1_c_ihi = c1_b_ihi;
      c1_b = c1_c_ihi;
      c1_b_b = c1_b;
      if (1 > c1_b_b) {
        c1_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_overflow = (c1_b_b > 2147483646);
      }

      if (c1_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_overflow);
      }

      c1_jj = 1;
      exitg6 = false;
      while ((exitg6 == false) && (c1_jj <= c1_c_ihi)) {
        c1_b_jj = c1_jj;
        c1_b_A = ((c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_ii), 1, 3, 1, 0) + 3 *
                         (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_jj), 1, 3, 2, 0) - 1)) - 1].re !=
                   c1_dc7.re) || (c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_ii), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_jj), 1, 3, 2, 0) - 1)) - 1].im !=
                   c1_dc7.im));
        guard3 = false;
        guard4 = false;
        if (c1_b_A) {
          guard4 = true;
        } else if (c1_ii == c1_b_jj) {
          guard4 = true;
        } else {
          guard3 = true;
        }

        if (guard4 == true) {
          if (c1_nzcount == 0.0) {
            c1_j = c1_b_jj;
            c1_nzcount = 1.0;
            guard3 = true;
          } else {
            c1_nzcount = 2.0;
            exitg6 = true;
          }
        }

        if (guard3 == true) {
          c1_jj++;
        }
      }

      if (c1_nzcount < 2.0) {
        c1_found = true;
        exitg5 = true;
      } else {
        c1_a = c1_ii;
        c1_b_a = c1_a - 1;
        c1_ii = c1_b_a;
      }
    }

    c1_b_i = c1_i;
    c1_b_j = c1_j;
    c1_b_found = c1_found;
    if (!c1_b_found) {
      exitg2 = 2;
    } else {
      c1_b_m = *c1_ihi;
      c1_e_i = c1_b_i;
      c1_e_j = c1_b_j;
      c1_i_ihi = *c1_ihi;
      if (c1_e_i != c1_b_m) {
        c1_eml_switch_helper(chartInstance);
        for (c1_d_k = 1; c1_d_k < 4; c1_d_k++) {
          c1_e_k = c1_d_k;
          c1_atmp.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_e_i), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_e_k), 1, 3, 2, 0) - 1)) - 1].re;
          c1_atmp.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_e_i), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_e_k), 1, 3, 2, 0) - 1)) - 1].im;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_e_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_e_k), 1, 3, 2, 0) - 1)) - 1].re = c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_b_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_e_k), 1, 3, 2, 0)
               - 1)) - 1].re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_e_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_e_k), 1, 3, 2, 0) - 1)) - 1].im = c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_b_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_e_k), 1, 3, 2, 0)
               - 1)) - 1].im;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_m), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_e_k), 1, 3, 2, 0) - 1)) - 1].re = c1_atmp.re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_m), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_e_k), 1, 3, 2, 0) - 1)) - 1].im = c1_atmp.im;
        }
      }

      if (c1_e_j != c1_b_m) {
        c1_j_ihi = c1_i_ihi;
        c1_i_b = c1_j_ihi;
        c1_j_b = c1_i_b;
        if (1 > c1_j_b) {
          c1_f_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_f_overflow = (c1_j_b > 2147483646);
        }

        if (c1_f_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_f_overflow);
        }

        for (c1_f_k = 1; c1_f_k <= c1_j_ihi; c1_f_k++) {
          c1_e_k = c1_f_k;
          c1_atmp.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_e_k), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_e_j), 1, 3, 2, 0) - 1)) - 1].re;
          c1_atmp.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_e_k), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_e_j), 1, 3, 2, 0) - 1)) - 1].im;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_e_j), 1, 3, 2, 0) - 1)) - 1].re = c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_e_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_m), 1, 3, 2, 0)
               - 1)) - 1].re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_e_j), 1, 3, 2, 0) - 1)) - 1].im = c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_e_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_m), 1, 3, 2, 0)
               - 1)) - 1].im;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_m), 1, 3, 2, 0) - 1)) - 1].re = c1_atmp.re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_m), 1, 3, 2, 0) - 1)) - 1].im = c1_atmp.im;
        }
      }

      c1_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)*c1_ihi), 1, 3, 1, 0) - 1] = c1_b_j;
      c1_k_a = *c1_ihi;
      c1_l_a = c1_k_a - 1;
      *c1_ihi = c1_l_a;
      if (*c1_ihi == 1) {
        c1_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)*c1_ihi), 1, 3, 1, 0) - 1] = *c1_ihi;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    do {
      exitg1 = 0;
      c1_b_ilo = *c1_ilo;
      c1_d_ihi = *c1_ihi;
      c1_c_i = 0;
      c1_c_j = 0;
      c1_c_found = false;
      c1_c_ilo = c1_b_ilo;
      c1_e_ihi = c1_d_ihi;
      c1_c_a = c1_c_ilo;
      c1_c_b = c1_e_ihi;
      c1_d_a = c1_c_a;
      c1_d_b = c1_c_b;
      if (c1_d_a > c1_d_b) {
        c1_b_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_b_overflow = (c1_d_b > 2147483646);
      }

      if (c1_b_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
      }

      c1_c_jj = c1_c_ilo;
      exitg3 = false;
      while ((exitg3 == false) && (c1_c_jj <= c1_e_ihi)) {
        c1_d_jj = c1_c_jj;
        c1_b_nzcount = 0.0;
        c1_c_i = c1_d_ihi;
        c1_c_j = c1_d_jj;
        c1_d_ilo = c1_b_ilo;
        c1_f_ihi = c1_d_ihi;
        c1_e_a = c1_d_ilo;
        c1_e_b = c1_f_ihi;
        c1_f_a = c1_e_a;
        c1_f_b = c1_e_b;
        if (c1_f_a > c1_f_b) {
          c1_c_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_c_overflow = (c1_f_b > 2147483646);
        }

        if (c1_c_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_c_overflow);
        }

        c1_b_ii = c1_d_ilo;
        exitg4 = false;
        while ((exitg4 == false) && (c1_b_ii <= c1_f_ihi)) {
          c1_c_ii = c1_b_ii;
          c1_c_A = ((c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_c_ii), 1, 3, 1, 0) + 3 *
                           (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_jj), 1, 3, 2, 0) - 1)) - 1].re
                     != c1_dc7.re) || (c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_c_ii), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_jj), 1, 3, 2, 0) - 1)) - 1].im
                     != c1_dc7.im));
          guard1 = false;
          guard2 = false;
          if (c1_c_A) {
            guard2 = true;
          } else if (c1_c_ii == c1_d_jj) {
            guard2 = true;
          } else {
            guard1 = true;
          }

          if (guard2 == true) {
            if (c1_b_nzcount == 0.0) {
              c1_c_i = c1_c_ii;
              c1_b_nzcount = 1.0;
              guard1 = true;
            } else {
              c1_b_nzcount = 2.0;
              exitg4 = true;
            }
          }

          if (guard1 == true) {
            c1_b_ii++;
          }
        }

        if (c1_b_nzcount < 2.0) {
          c1_c_found = true;
          exitg3 = true;
        } else {
          c1_c_jj++;
        }
      }

      c1_b_i = c1_c_i;
      c1_b_j = c1_c_j;
      c1_b_found = c1_c_found;
      if (!c1_b_found) {
        exitg1 = 1;
      } else {
        c1_m = *c1_ilo;
        c1_d_i = c1_b_i;
        c1_d_j = c1_b_j;
        c1_e_ilo = *c1_ilo;
        c1_g_ihi = *c1_ihi;
        if (c1_d_i != c1_m) {
          c1_f_ilo = c1_e_ilo;
          c1_g_a = c1_f_ilo;
          c1_h_a = c1_g_a;
          if (c1_h_a > 3) {
            c1_d_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_d_overflow = false;
          }

          if (c1_d_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_d_overflow);
          }

          for (c1_k = c1_f_ilo; c1_k < 4; c1_k++) {
            c1_b_k = c1_k;
            c1_atmp.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 2, 0) - 1)) - 1].re;
            c1_atmp.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 2, 0) - 1)) - 1].im;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_b_k), 1, 3, 2, 0) - 1)) - 1].re = c1_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                  "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 2,
                  0) - 1)) - 1].re;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_b_k), 1, 3, 2, 0) - 1)) - 1].im = c1_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                  "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 2,
                  0) - 1)) - 1].im;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_m), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_b_k), 1, 3, 2, 0) - 1)) - 1].re = c1_atmp.re;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_m), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_b_k), 1, 3, 2, 0) - 1)) - 1].im = c1_atmp.im;
          }
        }

        if (c1_d_j != c1_m) {
          c1_h_ihi = c1_g_ihi;
          c1_g_b = c1_h_ihi;
          c1_h_b = c1_g_b;
          if (1 > c1_h_b) {
            c1_e_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_e_overflow = (c1_h_b > 2147483646);
          }

          if (c1_e_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_e_overflow);
          }

          for (c1_c_k = 1; c1_c_k <= c1_h_ihi; c1_c_k++) {
            c1_b_k = c1_c_k;
            c1_atmp.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) - 1].re;
            c1_atmp.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) - 1].im;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) - 1].re = c1_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_m), 1, 3, 2, 0) - 1)) - 1].re;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) - 1].im = c1_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_m), 1, 3, 2, 0) - 1)) - 1].im;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_m), 1, 3, 2, 0) - 1)) - 1].re = c1_atmp.re;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_m), 1, 3, 2, 0) - 1)) - 1].im = c1_atmp.im;
          }
        }

        c1_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)*c1_ilo), 1, 3, 1, 0) - 1] = c1_b_j;
        c1_i_a = *c1_ilo;
        c1_j_a = c1_i_a + 1;
        *c1_ilo = c1_j_a;
        if (*c1_ilo == *c1_ihi) {
          c1_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
            "", (real_T)*c1_ilo), 1, 3, 1, 0) - 1] = *c1_ilo;
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

static void c1_b_sqrt(SFc1_Model_01InstanceStruct *chartInstance, creal_T *c1_x)
{
  real_T c1_yr;
  real_T c1_yi;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_z;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_b_z;
  boolean_T c1_b3;
  boolean_T c1_b4;
  boolean_T c1_b;
  real_T c1_h_x;
  boolean_T c1_b_b;
  real_T c1_i_x;
  boolean_T c1_c_b;
  real_T c1_absxr;
  real_T c1_absxi;
  real_T c1_j_x;
  real_T c1_y;
  real_T c1_x1;
  real_T c1_x2;
  real_T c1_a;
  real_T c1_d_b;
  real_T c1_absxd2;
  real_T c1_k_x;
  real_T c1_b_y;
  real_T c1_l_x;
  real_T c1_c_y;
  real_T c1_m_x;
  real_T c1_d_y;
  real_T c1_c_z;
  real_T c1_n_x;
  real_T c1_e_y;
  real_T c1_b_x1;
  real_T c1_b_x2;
  real_T c1_b_a;
  real_T c1_e_b;
  real_T c1_d_z;
  real_T c1_o_x;
  real_T c1_f_y;
  real_T c1_p_x;
  real_T c1_g_y;
  real_T c1_q_x;
  real_T c1_h_y;
  real_T c1_e_z;
  real_T c1_r_x;
  real_T c1_i_y;
  real_T c1_s_x;
  real_T c1_j_y;
  real_T c1_t_x;
  real_T c1_k_y;
  real_T c1_f_z;
  boolean_T guard1 = false;
  if (c1_x->im == 0.0) {
    if (c1_x->re < 0.0) {
      c1_yr = 0.0;
      c1_yi = muDoubleScalarSqrt(muDoubleScalarAbs(c1_x->re));
    } else {
      c1_yr = muDoubleScalarSqrt(c1_x->re);
      c1_yi = 0.0;
    }
  } else if (c1_x->re == 0.0) {
    if (c1_x->im < 0.0) {
      c1_b_x = -c1_x->im;
      c1_c_x = c1_b_x;
      c1_d_x = c1_c_x;
      c1_z = c1_d_x / 2.0;
      c1_yr = muDoubleScalarSqrt(c1_z);
      c1_yi = -c1_yr;
    } else {
      c1_e_x = c1_x->im;
      c1_f_x = c1_e_x;
      c1_g_x = c1_f_x;
      c1_b_z = c1_g_x / 2.0;
      c1_yr = muDoubleScalarSqrt(c1_b_z);
      c1_yi = c1_yr;
    }
  } else {
    c1_b3 = muDoubleScalarIsNaN(c1_x->re);
    c1_b4 = muDoubleScalarIsNaN(c1_x->im);
    c1_b = (c1_b3 || c1_b4);
    if (c1_b) {
      c1_yr = rtNaN;
      c1_yi = rtNaN;
    } else {
      c1_h_x = c1_x->im;
      c1_b_b = muDoubleScalarIsInf(c1_h_x);
      if (c1_b_b) {
        c1_yr = rtInf;
        c1_yi = c1_x->im;
      } else {
        c1_i_x = c1_x->re;
        c1_c_b = muDoubleScalarIsInf(c1_i_x);
        if (c1_c_b) {
          if (c1_x->re < 0.0) {
            c1_yr = 0.0;
            c1_yi = rtInf;
          } else {
            c1_yr = rtInf;
            c1_yi = 0.0;
          }
        } else {
          c1_absxr = muDoubleScalarAbs(c1_x->re);
          c1_absxi = muDoubleScalarAbs(c1_x->im);
          c1_realmax(chartInstance);
          guard1 = false;
          if (c1_absxr > 4.4942328371557893E+307) {
            guard1 = true;
          } else {
            c1_realmax(chartInstance);
            if (c1_absxi > 4.4942328371557893E+307) {
              guard1 = true;
            } else {
              c1_n_x = c1_absxr;
              c1_e_y = c1_absxi;
              c1_eml_scalar_eg(chartInstance);
              c1_b_x1 = c1_n_x;
              c1_b_x2 = c1_e_y;
              c1_b_a = c1_b_x1;
              c1_e_b = c1_b_x2;
              c1_d_z = muDoubleScalarHypot(c1_b_a, c1_e_b);
              c1_yr = muDoubleScalarSqrt((c1_d_z + c1_absxr) * 0.5);
            }
          }

          if (guard1 == true) {
            c1_absxr *= 0.5;
            c1_absxi *= 0.5;
            c1_j_x = c1_absxr;
            c1_y = c1_absxi;
            c1_eml_scalar_eg(chartInstance);
            c1_x1 = c1_j_x;
            c1_x2 = c1_y;
            c1_a = c1_x1;
            c1_d_b = c1_x2;
            c1_absxd2 = muDoubleScalarHypot(c1_a, c1_d_b);
            if (c1_absxd2 > c1_absxr) {
              c1_k_x = c1_absxr;
              c1_b_y = c1_absxd2;
              c1_l_x = c1_k_x;
              c1_c_y = c1_b_y;
              c1_m_x = c1_l_x;
              c1_d_y = c1_c_y;
              c1_c_z = c1_m_x / c1_d_y;
              c1_yr = muDoubleScalarSqrt(c1_absxd2) * muDoubleScalarSqrt(1.0 +
                c1_c_z);
            } else {
              c1_yr = muDoubleScalarSqrt(c1_absxd2) * 1.4142135623730951;
            }
          }

          if (c1_x->re > 0.0) {
            c1_o_x = c1_x->im;
            c1_f_y = c1_yr;
            c1_p_x = c1_o_x;
            c1_g_y = c1_f_y;
            c1_q_x = c1_p_x;
            c1_h_y = c1_g_y;
            c1_e_z = c1_q_x / c1_h_y;
            c1_yi = 0.5 * c1_e_z;
          } else {
            if (c1_x->im < 0.0) {
              c1_yi = -c1_yr;
            } else {
              c1_yi = c1_yr;
            }

            c1_r_x = c1_x->im;
            c1_i_y = c1_yi;
            c1_s_x = c1_r_x;
            c1_j_y = c1_i_y;
            c1_t_x = c1_s_x;
            c1_k_y = c1_j_y;
            c1_f_z = c1_t_x / c1_k_y;
            c1_yr = 0.5 * c1_f_z;
          }
        }
      }
    }
  }

  c1_x->re = c1_yr;
  c1_x->im = c1_yi;
}

static void c1_d_eml_matlab_zlascl(SFc1_Model_01InstanceStruct *chartInstance,
  real_T c1_cfrom, real_T c1_cto, creal_T c1_A[3])
{
  real_T c1_cfromc;
  real_T c1_ctoc;
  boolean_T c1_notdone;
  real_T c1_cfrom1;
  real_T c1_cto1;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_y;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_b_y;
  real_T c1_mul;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_c_y;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_d_y;
  real_T c1_a;
  int32_T c1_i118;
  boolean_T guard1 = false;
  c1_realmin(chartInstance);
  c1_eps(chartInstance);
  c1_cfromc = c1_cfrom;
  c1_ctoc = c1_cto;
  c1_notdone = true;
  while (c1_notdone) {
    c1_cfrom1 = c1_cfromc * 2.0041683600089728E-292;
    c1_cto1 = c1_ctoc / 4.9896007738368E+291;
    c1_x = c1_cfrom1;
    c1_b_x = c1_x;
    c1_y = muDoubleScalarAbs(c1_b_x);
    c1_c_x = c1_ctoc;
    c1_d_x = c1_c_x;
    c1_b_y = muDoubleScalarAbs(c1_d_x);
    guard1 = false;
    if (c1_y > c1_b_y) {
      if (c1_ctoc != 0.0) {
        c1_mul = 2.0041683600089728E-292;
        c1_notdone = true;
        c1_cfromc = c1_cfrom1;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      c1_e_x = c1_cto1;
      c1_f_x = c1_e_x;
      c1_c_y = muDoubleScalarAbs(c1_f_x);
      c1_g_x = c1_cfromc;
      c1_h_x = c1_g_x;
      c1_d_y = muDoubleScalarAbs(c1_h_x);
      if (c1_c_y > c1_d_y) {
        c1_mul = 4.9896007738368E+291;
        c1_notdone = true;
        c1_ctoc = c1_cto1;
      } else {
        c1_mul = c1_ctoc / c1_cfromc;
        c1_notdone = false;
      }
    }

    c1_a = c1_mul;
    for (c1_i118 = 0; c1_i118 < 3; c1_i118++) {
      c1_A[c1_i118].re *= c1_a;
      c1_A[c1_i118].im *= c1_a;
    }
  }
}

static void c1_b_eml_xgemm(SFc1_Model_01InstanceStruct *chartInstance, real_T
  c1_A[9], real_T c1_B[9], real_T c1_C[9])
{
  int32_T c1_i119;
  int32_T c1_i120;
  int32_T c1_i121;
  int32_T c1_i122;
  int32_T c1_i123;
  (void)chartInstance;
  for (c1_i119 = 0; c1_i119 < 3; c1_i119++) {
    c1_i120 = 0;
    for (c1_i121 = 0; c1_i121 < 3; c1_i121++) {
      c1_C[c1_i120 + c1_i119] = 0.0;
      c1_i122 = 0;
      for (c1_i123 = 0; c1_i123 < 3; c1_i123++) {
        c1_C[c1_i120 + c1_i119] += c1_A[c1_i122 + c1_i119] * c1_B[c1_i123 +
          c1_i120];
        c1_i122 += 3;
      }

      c1_i120 += 3;
    }
  }
}

static void c1_b_eml_matlab_zgghrd(SFc1_Model_01InstanceStruct *chartInstance,
  int32_T c1_ilo, int32_T c1_ihi, creal_T c1_A[9], creal_T c1_Z[9])
{
  int32_T c1_i124;
  static real_T c1_dv4[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  int32_T c1_a;
  int32_T c1_b_a;
  int32_T c1_c;
  int32_T c1_c_a;
  int32_T c1_d_a;
  int32_T c1_ihim1;
  int32_T c1_jcol;
  int32_T c1_e_a;
  int32_T c1_f_a;
  int32_T c1_jcolp1;
  int32_T c1_jrow;
  int32_T c1_g_a;
  int32_T c1_h_a;
  int32_T c1_jrowm1;
  creal_T c1_b_A;
  creal_T c1_c_A;
  creal_T c1_b;
  creal_T c1_s;
  real_T c1_b_c;
  real_T c1_c_c;
  static creal_T c1_dc8 = { 0.0, 0.0 };

  real_T c1_d_c;
  int32_T c1_xrow;
  int32_T c1_yrow;
  int32_T c1_jlo;
  int32_T c1_jhi;
  int32_T c1_b_jlo;
  int32_T c1_b_jhi;
  int32_T c1_i_a;
  int32_T c1_b_b;
  int32_T c1_j_a;
  int32_T c1_c_b;
  boolean_T c1_overflow;
  int32_T c1_j;
  int32_T c1_b_j;
  real_T c1_k_a;
  creal_T c1_y;
  creal_T c1_b_s;
  creal_T c1_stemp;
  real_T c1_l_a;
  creal_T c1_d_b;
  creal_T c1_e_b;
  creal_T c1_f_b;
  creal_T c1_g_b;
  real_T c1_e_c;
  int32_T c1_xcol;
  int32_T c1_ycol;
  int32_T c1_b_ilo;
  int32_T c1_b_ihi;
  int32_T c1_c_ilo;
  int32_T c1_c_ihi;
  int32_T c1_m_a;
  int32_T c1_h_b;
  int32_T c1_n_a;
  int32_T c1_i_b;
  boolean_T c1_b_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  real_T c1_o_a;
  creal_T c1_c_s;
  real_T c1_p_a;
  creal_T c1_j_b;
  creal_T c1_k_b;
  creal_T c1_l_b;
  creal_T c1_m_b;
  real_T c1_f_c;
  int32_T c1_b_xcol;
  int32_T c1_b_ycol;
  int32_T c1_c_i;
  int32_T c1_d_i;
  real_T c1_q_a;
  creal_T c1_d_s;
  real_T c1_r_a;
  creal_T c1_n_b;
  creal_T c1_o_b;
  creal_T c1_p_b;
  creal_T c1_q_b;
  for (c1_i124 = 0; c1_i124 < 9; c1_i124++) {
    c1_Z[c1_i124].re = c1_dv4[c1_i124];
    c1_Z[c1_i124].im = 0.0;
  }

  c1_a = c1_ilo;
  c1_b_a = c1_a;
  c1_c = c1_b_a;
  if (c1_ihi < c1_c + 2) {
  } else {
    c1_c_a = c1_ihi;
    c1_d_a = c1_c_a;
    c1_ihim1 = c1_d_a;
    c1_jcol = c1_ilo;
    while (c1_jcol < c1_ihim1 - 1) {
      c1_e_a = c1_jcol;
      c1_f_a = c1_e_a + 1;
      c1_jcolp1 = c1_f_a;
      c1_jrow = c1_ihi;
      while (c1_jrow > c1_jcolp1) {
        c1_g_a = c1_jrow;
        c1_h_a = c1_g_a - 1;
        c1_jrowm1 = c1_h_a;
        c1_b_A.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_jrowm1), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].re;
        c1_b_A.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_jrowm1), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].im;
        c1_c_A.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_jrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].re;
        c1_c_A.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_jrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].im;
        c1_eml_matlab_zlartg(chartInstance, c1_b_A, c1_c_A, &c1_b_c, &c1_s,
                             &c1_b);
        c1_c_c = c1_b_c;
        c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_jrowm1), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].re = c1_b.re;
        c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_jrowm1), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c1_jcol), 1, 3, 2, 0) - 1)) - 1].im = c1_b.im;
        c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_jrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0)
               - 1)) - 1].re = c1_dc8.re;
        c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_jrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_jcol), 1, 3, 2, 0)
               - 1)) - 1].im = c1_dc8.im;
        c1_d_c = c1_c_c;
        c1_xrow = c1_jrowm1;
        c1_yrow = c1_jrow;
        c1_jlo = c1_jcolp1;
        c1_jhi = c1_ihi;
        c1_b_jlo = c1_jlo;
        c1_b_jhi = c1_jhi;
        c1_i_a = c1_b_jlo;
        c1_b_b = c1_b_jhi;
        c1_j_a = c1_i_a;
        c1_c_b = c1_b_b;
        if (c1_j_a > c1_c_b) {
          c1_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_overflow = (c1_c_b > 2147483646);
        }

        if (c1_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_overflow);
        }

        for (c1_j = c1_b_jlo; c1_j <= c1_b_jhi; c1_j++) {
          c1_b_j = c1_j;
          c1_k_a = c1_d_c;
          c1_b.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
          c1_b.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
          c1_y.re = c1_k_a * c1_b.re;
          c1_y.im = c1_k_a * c1_b.im;
          c1_b_s.re = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re - c1_s.im * c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_yrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0)
               - 1)) - 1].im;
          c1_b_s.im = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im + c1_s.im * c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_yrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0)
               - 1)) - 1].re;
          c1_stemp.re = c1_y.re + c1_b_s.re;
          c1_stemp.im = c1_y.im + c1_b_s.im;
          c1_l_a = c1_d_c;
          c1_b.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
          c1_b.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
          c1_y.re = c1_l_a * c1_b.re;
          c1_y.im = c1_l_a * c1_b.im;
          c1_b.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
          c1_b.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
          c1_d_b = c1_b;
          c1_e_b = c1_b;
          c1_f_b = c1_b;
          c1_g_b = c1_b;
          c1_b.re = c1_s.re * c1_d_b.re + c1_s.im * c1_e_b.im;
          c1_b.im = c1_s.re * c1_f_b.im - c1_s.im * c1_g_b.re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re = c1_y.re -
            c1_b.re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im = c1_y.im -
            c1_b.im;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].re = c1_stemp.re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1].im = c1_stemp.im;
        }

        c1_s.re = -c1_s.re;
        c1_s.im = -c1_s.im;
        c1_e_c = c1_c_c;
        c1_xcol = c1_jrow;
        c1_ycol = c1_jrowm1;
        c1_b_ilo = c1_ilo;
        c1_b_ihi = c1_ihi;
        c1_c_ilo = c1_b_ilo;
        c1_c_ihi = c1_b_ihi;
        c1_m_a = c1_c_ilo;
        c1_h_b = c1_c_ihi;
        c1_n_a = c1_m_a;
        c1_i_b = c1_h_b;
        if (c1_n_a > c1_i_b) {
          c1_b_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_b_overflow = (c1_i_b > 2147483646);
        }

        if (c1_b_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
        }

        for (c1_i = c1_c_ilo; c1_i <= c1_c_ihi; c1_i++) {
          c1_b_i = c1_i;
          c1_o_a = c1_e_c;
          c1_b.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c1_b.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c1_y.re = c1_o_a * c1_b.re;
          c1_y.im = c1_o_a * c1_b.im;
          c1_c_s.re = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re - c1_s.im * c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0)
               - 1)) - 1].im;
          c1_c_s.im = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im + c1_s.im * c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0)
               - 1)) - 1].re;
          c1_stemp.re = c1_y.re + c1_c_s.re;
          c1_stemp.im = c1_y.im + c1_c_s.im;
          c1_p_a = c1_e_c;
          c1_b.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re;
          c1_b.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im;
          c1_y.re = c1_p_a * c1_b.re;
          c1_y.im = c1_p_a * c1_b.im;
          c1_b.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c1_b.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c1_j_b = c1_b;
          c1_k_b = c1_b;
          c1_l_b = c1_b;
          c1_m_b = c1_b;
          c1_b.re = c1_s.re * c1_j_b.re + c1_s.im * c1_k_b.im;
          c1_b.im = c1_s.re * c1_l_b.im - c1_s.im * c1_m_b.re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re = c1_y.re -
            c1_b.re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im = c1_y.im -
            c1_b.im;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].re = c1_stemp.re;
          c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].im = c1_stemp.im;
        }

        c1_f_c = c1_c_c;
        c1_b_xcol = c1_jrow;
        c1_b_ycol = c1_jrowm1;
        c1_eml_switch_helper(chartInstance);
        for (c1_c_i = 1; c1_c_i < 4; c1_c_i++) {
          c1_d_i = c1_c_i;
          c1_q_a = c1_f_c;
          c1_b.re = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c1_b.im = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c1_y.re = c1_q_a * c1_b.re;
          c1_y.im = c1_q_a * c1_b.im;
          c1_d_s.re = c1_s.re * c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].re - c1_s.im * c1_Z
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_d_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2,
                0) - 1)) - 1].im;
          c1_d_s.im = c1_s.re * c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].im + c1_s.im * c1_Z
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_d_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2,
                0) - 1)) - 1].re;
          c1_stemp.re = c1_y.re + c1_d_s.re;
          c1_stemp.im = c1_y.im + c1_d_s.im;
          c1_r_a = c1_f_c;
          c1_b.re = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].re;
          c1_b.im = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].im;
          c1_y.re = c1_r_a * c1_b.re;
          c1_y.im = c1_r_a * c1_b.im;
          c1_b.re = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c1_b.im = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c1_n_b = c1_b;
          c1_o_b = c1_b;
          c1_p_b = c1_b;
          c1_q_b = c1_b;
          c1_b.re = c1_s.re * c1_n_b.re + c1_s.im * c1_o_b.im;
          c1_b.im = c1_s.re * c1_p_b.im - c1_s.im * c1_q_b.re;
          c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].re = c1_y.re -
            c1_b.re;
          c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].im = c1_y.im -
            c1_b.im;
          c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1].re = c1_stemp.re;
          c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1].im = c1_stemp.im;
        }

        c1_jrow = c1_jrowm1;
      }

      c1_jcol = c1_jcolp1;
    }
  }
}

static void c1_c_eml_matlab_zhgeqz(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], int32_T c1_ilo, int32_T c1_ihi, creal_T c1_Z[9], real_T
  *c1_info, creal_T c1_alpha1[3], creal_T c1_beta1[3])
{
  static creal_T c1_dc9 = { 0.0, 0.0 };

  int32_T c1_i125;
  int32_T c1_i126;
  static creal_T c1_dc10 = { 0.0, 0.0 };

  creal_T c1_eshift;
  creal_T c1_ctemp;
  creal_T c1_rho;
  int32_T c1_i127;
  int32_T c1_i128;
  int32_T c1_i129;
  creal_T c1_b_A[9];
  real_T c1_anorm;
  real_T c1_y;
  real_T c1_atol;
  real_T c1_b_y;
  real_T c1_x;
  real_T c1_ascale;
  boolean_T c1_failed;
  int32_T c1_a;
  int32_T c1_b_a;
  int32_T c1_i130;
  int32_T c1_c_a;
  int32_T c1_d_a;
  boolean_T c1_overflow;
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_ifirst;
  int32_T c1_istart;
  int32_T c1_ilast;
  int32_T c1_e_a;
  int32_T c1_f_a;
  int32_T c1_ilastm1;
  int32_T c1_iiter;
  int32_T c1_g_a;
  int32_T c1_b;
  int32_T c1_h_a;
  int32_T c1_b_b;
  int32_T c1_c;
  int32_T c1_i_a;
  int32_T c1_j_a;
  int32_T c1_b_c;
  int32_T c1_c_b;
  int32_T c1_d_b;
  int32_T c1_maxit;
  boolean_T c1_goto50;
  boolean_T c1_goto60;
  boolean_T c1_goto70;
  boolean_T c1_goto90;
  int32_T c1_b_maxit;
  int32_T c1_e_b;
  int32_T c1_f_b;
  boolean_T c1_b_overflow;
  int32_T c1_jiter;
  creal_T c1_a22;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_c_y;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_d_y;
  real_T c1_e_y;
  int32_T c1_k_a;
  int32_T c1_l_a;
  int32_T c1_jm1;
  boolean_T c1_ilazro;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_f_y;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_g_y;
  real_T c1_h_y;
  boolean_T c1_b5;
  int32_T c1_i131;
  int32_T c1_i132;
  int32_T c1_i133;
  creal_T c1_c_A;
  creal_T c1_d_A;
  creal_T c1_s;
  real_T c1_c_c;
  real_T c1_d_c;
  real_T c1_e_c;
  int32_T c1_xcol;
  int32_T c1_ycol;
  int32_T c1_b_ihi;
  int32_T c1_c_ihi;
  int32_T c1_g_b;
  int32_T c1_h_b;
  boolean_T c1_c_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  real_T c1_m_a;
  creal_T c1_a12;
  creal_T c1_b_s;
  creal_T c1_a21;
  real_T c1_n_a;
  creal_T c1_b_a22;
  creal_T c1_c_a22;
  creal_T c1_d_a22;
  creal_T c1_e_a22;
  real_T c1_f_c;
  int32_T c1_b_xcol;
  int32_T c1_b_ycol;
  int32_T c1_c_i;
  int32_T c1_d_i;
  real_T c1_o_a;
  creal_T c1_c_s;
  real_T c1_p_a;
  creal_T c1_f_a22;
  creal_T c1_g_a22;
  creal_T c1_h_a22;
  creal_T c1_i_a22;
  int32_T c1_q_a;
  int32_T c1_r_a;
  int32_T c1_s_a;
  int32_T c1_t_a;
  creal_T c1_r2;
  creal_T c1_j_a22;
  creal_T c1_b_rho;
  creal_T c1_b_a12;
  creal_T c1_c_a12;
  creal_T c1_b_a21;
  real_T c1_d7;
  real_T c1_d8;
  int32_T c1_u_a;
  int32_T c1_v_a;
  int32_T c1_jp1;
  int32_T c1_w_a;
  int32_T c1_x_a;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_i_y;
  real_T c1_l_x;
  real_T c1_m_x;
  real_T c1_j_y;
  real_T c1_k_y;
  real_T c1_temp;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_l_y;
  real_T c1_p_x;
  real_T c1_q_x;
  real_T c1_m_y;
  real_T c1_n_y;
  real_T c1_temp2;
  real_T c1_r_x;
  real_T c1_o_y;
  real_T c1_tempr;
  real_T c1_s_x;
  real_T c1_t_x;
  real_T c1_p_y;
  real_T c1_u_x;
  real_T c1_v_x;
  real_T c1_q_y;
  real_T c1_r_y;
  int32_T c1_y_a;
  int32_T c1_ab_a;
  int32_T c1_g_c;
  real_T c1_h_c;
  int32_T c1_bb_a;
  int32_T c1_cb_a;
  int32_T c1_db_a;
  int32_T c1_eb_a;
  creal_T c1_e_A;
  creal_T c1_f_A;
  real_T c1_i_c;
  real_T c1_j_c;
  int32_T c1_xrow;
  int32_T c1_yrow;
  int32_T c1_jlo;
  int32_T c1_b_jlo;
  int32_T c1_fb_a;
  int32_T c1_gb_a;
  boolean_T c1_d_overflow;
  int32_T c1_c_j;
  int32_T c1_d_j;
  real_T c1_hb_a;
  creal_T c1_d_s;
  real_T c1_ib_a;
  creal_T c1_k_a22;
  creal_T c1_l_a22;
  creal_T c1_m_a22;
  creal_T c1_n_a22;
  int32_T c1_jb_a;
  int32_T c1_kb_a;
  int32_T c1_k_c;
  int32_T c1_w_x;
  int32_T c1_s_y;
  int32_T c1_x_x;
  real_T c1_l_c;
  int32_T c1_c_xcol;
  int32_T c1_c_ycol;
  int32_T c1_d_ihi;
  int32_T c1_e_ihi;
  int32_T c1_i_b;
  int32_T c1_j_b;
  boolean_T c1_e_overflow;
  int32_T c1_e_i;
  int32_T c1_f_i;
  real_T c1_lb_a;
  creal_T c1_e_s;
  real_T c1_mb_a;
  creal_T c1_o_a22;
  creal_T c1_p_a22;
  creal_T c1_q_a22;
  creal_T c1_r_a22;
  real_T c1_m_c;
  int32_T c1_d_xcol;
  int32_T c1_d_ycol;
  int32_T c1_g_i;
  int32_T c1_h_i;
  real_T c1_nb_a;
  creal_T c1_f_s;
  real_T c1_ob_a;
  creal_T c1_s_a22;
  creal_T c1_t_a22;
  creal_T c1_u_a22;
  creal_T c1_v_a22;
  int32_T c1_b_ilast;
  int32_T c1_k_b;
  int32_T c1_l_b;
  boolean_T c1_f_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  int32_T c1_i134;
  int32_T c1_pb_a;
  int32_T c1_qb_a;
  int32_T c1_i135;
  int32_T c1_m_b;
  int32_T c1_n_b;
  boolean_T c1_g_overflow;
  int32_T c1_e_j;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  int32_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T guard11 = false;
  c1_dc9.re = rtNaN;
  for (c1_i125 = 0; c1_i125 < 3; c1_i125++) {
    c1_alpha1[c1_i125].re = 0.0;
    c1_alpha1[c1_i125].im = 0.0;
  }

  for (c1_i126 = 0; c1_i126 < 3; c1_i126++) {
    c1_beta1[c1_i126].re = 1.0;
    c1_beta1[c1_i126].im = 0.0;
  }

  c1_eps(chartInstance);
  c1_realmin(chartInstance);
  c1_eshift = c1_dc10;
  c1_ctemp = c1_dc10;
  c1_rho = c1_dc10;
  c1_i127 = 0;
  for (c1_i128 = 0; c1_i128 < 3; c1_i128++) {
    for (c1_i129 = 0; c1_i129 < 3; c1_i129++) {
      c1_b_A[c1_i129 + c1_i127] = c1_A[c1_i129 + c1_i127];
    }

    c1_i127 += 3;
  }

  c1_anorm = c1_eml_matlab_zlanhs(chartInstance, c1_b_A, c1_ilo, c1_ihi);
  c1_y = 2.2204460492503131E-16 * c1_anorm;
  c1_atol = 2.2250738585072014E-308;
  if (c1_y > 2.2250738585072014E-308) {
    c1_atol = c1_y;
  }

  c1_b_y = c1_anorm;
  c1_x = 2.2250738585072014E-308;
  if (c1_b_y > 2.2250738585072014E-308) {
    c1_x = c1_b_y;
  }

  c1_ascale = 1.0 / c1_x;
  c1_failed = true;
  c1_a = c1_ihi;
  c1_b_a = c1_a + 1;
  c1_i130 = c1_b_a;
  c1_c_a = c1_i130;
  c1_d_a = c1_c_a;
  if (c1_d_a > 3) {
    c1_overflow = false;
  } else {
    c1_eml_switch_helper(chartInstance);
    c1_overflow = false;
  }

  if (c1_overflow) {
    c1_check_forloop_overflow_error(chartInstance, c1_overflow);
  }

  for (c1_j = c1_i130; c1_j < 4; c1_j++) {
    c1_b_j = c1_j;
    c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_b_j), 1, 3, 1, 0) - 1].re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK(
      "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
    c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_b_j), 1, 3, 1, 0) - 1].im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK(
      "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
  }

  guard1 = false;
  guard2 = false;
  if (c1_ihi >= c1_ilo) {
    c1_ifirst = c1_ilo;
    c1_istart = c1_ilo;
    c1_ilast = c1_ihi;
    c1_e_a = c1_ilast;
    c1_f_a = c1_e_a - 1;
    c1_ilastm1 = c1_f_a;
    c1_iiter = 0;
    c1_g_a = c1_ihi;
    c1_b = c1_ilo;
    c1_h_a = c1_g_a;
    c1_b_b = c1_b;
    c1_c = c1_h_a - c1_b_b;
    c1_i_a = c1_c;
    c1_j_a = c1_i_a;
    c1_b_c = c1_j_a;
    c1_c_b = c1_b_c + 1;
    c1_d_b = c1_c_b;
    c1_maxit = 30 * c1_d_b;
    c1_goto50 = false;
    c1_goto60 = false;
    c1_goto70 = false;
    c1_goto90 = false;
    c1_b_maxit = c1_maxit;
    c1_e_b = c1_b_maxit;
    c1_f_b = c1_e_b;
    if (1 > c1_f_b) {
      c1_b_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_b_overflow = (c1_f_b > 2147483646);
    }

    if (c1_b_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
    }

    c1_jiter = 1;
    do {
      exitg1 = 0;
      if (c1_jiter <= c1_b_maxit) {
        if (c1_ilast == c1_ilo) {
          c1_goto60 = true;
        } else {
          c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            re;
          c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            im;
          c1_b_x = c1_a22.re;
          c1_c_x = c1_b_x;
          c1_c_y = muDoubleScalarAbs(c1_c_x);
          c1_d_x = c1_a22.im;
          c1_e_x = c1_d_x;
          c1_d_y = muDoubleScalarAbs(c1_e_x);
          c1_e_y = c1_c_y + c1_d_y;
          if (c1_e_y <= c1_atol) {
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c1_dc10.re;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c1_dc10.im;
            c1_goto60 = true;
          } else {
            c1_b_j = c1_ilastm1;
            exitg3 = false;
            while ((exitg3 == false) && (c1_b_j >= c1_ilo)) {
              c1_k_a = c1_b_j;
              c1_l_a = c1_k_a - 1;
              c1_jm1 = c1_l_a;
              if (c1_b_j == c1_ilo) {
                c1_ilazro = true;
              } else {
                c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c1_f_x = c1_a22.re;
                c1_g_x = c1_f_x;
                c1_f_y = muDoubleScalarAbs(c1_g_x);
                c1_h_x = c1_a22.im;
                c1_i_x = c1_h_x;
                c1_g_y = muDoubleScalarAbs(c1_i_x);
                c1_h_y = c1_f_y + c1_g_y;
                if (c1_h_y <= c1_atol) {
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) -
                           1)) - 1].re = c1_dc10.re;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) -
                           1)) - 1].im = c1_dc10.im;
                  c1_ilazro = true;
                } else {
                  c1_ilazro = false;
                }
              }

              if (c1_ilazro) {
                c1_ifirst = c1_b_j;
                c1_goto70 = true;
                exitg3 = true;
              } else {
                c1_b_j = c1_jm1;
              }
            }
          }
        }

        guard3 = false;
        guard4 = false;
        if (c1_goto50) {
          guard4 = true;
        } else if (c1_goto60) {
          guard4 = true;
        } else if (c1_goto70) {
          guard3 = true;
        } else {
          c1_b5 = false;
        }

        if (guard4 == true) {
          guard3 = true;
        }

        if (guard3 == true) {
          c1_b5 = true;
        }

        if (!c1_b5) {
          for (c1_i131 = 0; c1_i131 < 3; c1_i131++) {
            c1_alpha1[c1_i131] = c1_dc9;
          }

          for (c1_i132 = 0; c1_i132 < 3; c1_i132++) {
            c1_beta1[c1_i132] = c1_dc9;
          }

          for (c1_i133 = 0; c1_i133 < 9; c1_i133++) {
            c1_Z[c1_i133] = c1_dc9;
          }

          *c1_info = -1.0;
          exitg1 = 1;
        } else {
          if (c1_goto50) {
            c1_goto50 = false;
            c1_c_A.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].
              re;
            c1_c_A.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].
              im;
            c1_d_A.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .re;
            c1_d_A.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .im;
            c1_eml_matlab_zlartg(chartInstance, c1_c_A, c1_d_A, &c1_c_c, &c1_s,
                                 &c1_a22);
            c1_d_c = c1_c_c;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].re =
              c1_a22.re;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].im =
              c1_a22.im;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c1_dc10.re;
            c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c1_dc10.im;
            c1_e_c = c1_d_c;
            c1_xcol = c1_ilast;
            c1_ycol = c1_ilastm1;
            c1_b_ihi = c1_ilastm1;
            c1_c_ihi = c1_b_ihi;
            c1_g_b = c1_c_ihi;
            c1_h_b = c1_g_b;
            if (1 > c1_h_b) {
              c1_c_overflow = false;
            } else {
              c1_eml_switch_helper(chartInstance);
              c1_c_overflow = (c1_h_b > 2147483646);
            }

            if (c1_c_overflow) {
              c1_check_forloop_overflow_error(chartInstance, c1_c_overflow);
            }

            for (c1_i = 1; c1_i <= c1_c_ihi; c1_i++) {
              c1_b_i = c1_i;
              c1_m_a = c1_e_c;
              c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c1_a12.re = c1_m_a * c1_a22.re;
              c1_a12.im = c1_m_a * c1_a22.im;
              c1_b_s.re = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re - c1_s.im *
                c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1))
                - 1].im;
              c1_b_s.im = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im + c1_s.im *
                c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1))
                - 1].re;
              c1_a21.re = c1_a12.re + c1_b_s.re;
              c1_a21.im = c1_a12.im + c1_b_s.im;
              c1_n_a = c1_e_c;
              c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c1_a12.re = c1_n_a * c1_a22.re;
              c1_a12.im = c1_n_a * c1_a22.im;
              c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c1_b_a22 = c1_a22;
              c1_c_a22 = c1_a22;
              c1_d_a22 = c1_a22;
              c1_e_a22 = c1_a22;
              c1_a22.re = c1_s.re * c1_b_a22.re + c1_s.im * c1_c_a22.im;
              c1_a22.im = c1_s.re * c1_d_a22.im - c1_s.im * c1_e_a22.re;
              c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].re =
                c1_a12.re - c1_a22.re;
              c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ycol), 1, 3, 2, 0) - 1)) - 1].im =
                c1_a12.im - c1_a22.im;
              c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].re =
                c1_a21.re;
              c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_xcol), 1, 3, 2, 0) - 1)) - 1].im =
                c1_a21.im;
            }

            c1_f_c = c1_d_c;
            c1_b_xcol = c1_ilast;
            c1_b_ycol = c1_ilastm1;
            c1_eml_switch_helper(chartInstance);
            for (c1_c_i = 1; c1_c_i < 4; c1_c_i++) {
              c1_d_i = c1_c_i;
              c1_o_a = c1_f_c;
              c1_a22.re = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c1_a22.im = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c1_a12.re = c1_o_a * c1_a22.re;
              c1_a12.im = c1_o_a * c1_a22.im;
              c1_c_s.re = c1_s.re * c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].re - c1_s.im *
                c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) -
                       1)) - 1].im;
              c1_c_s.im = c1_s.re * c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].im + c1_s.im *
                c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) -
                       1)) - 1].re;
              c1_a21.re = c1_a12.re + c1_c_s.re;
              c1_a21.im = c1_a12.im + c1_c_s.im;
              c1_p_a = c1_f_c;
              c1_a22.re = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c1_a22.im = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c1_a12.re = c1_p_a * c1_a22.re;
              c1_a12.im = c1_p_a * c1_a22.im;
              c1_a22.re = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c1_a22.im = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c1_f_a22 = c1_a22;
              c1_g_a22 = c1_a22;
              c1_h_a22 = c1_a22;
              c1_i_a22 = c1_a22;
              c1_a22.re = c1_s.re * c1_f_a22.re + c1_s.im * c1_g_a22.im;
              c1_a22.im = c1_s.re * c1_h_a22.im - c1_s.im * c1_i_a22.re;
              c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].re =
                c1_a12.re - c1_a22.re;
              c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_ycol), 1, 3, 2, 0) - 1)) - 1].im =
                c1_a12.im - c1_a22.im;
              c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1].re =
                c1_a21.re;
              c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_xcol), 1, 3, 2, 0) - 1)) - 1].im =
                c1_a21.im;
            }

            c1_goto60 = true;
          }

          guard11 = false;
          if (c1_goto60) {
            c1_goto60 = false;
            c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) - 1].re =
              c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].re;
            c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) - 1].im =
              c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) - 1].im;
            c1_ilast = c1_ilastm1;
            c1_q_a = c1_ilast;
            c1_r_a = c1_q_a - 1;
            c1_ilastm1 = c1_r_a;
            if (c1_ilast < c1_ilo) {
              c1_failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              c1_iiter = 0;
              c1_eshift = c1_dc10;
              guard11 = true;
            }
          } else {
            if (c1_goto70) {
              c1_goto70 = false;
              c1_s_a = c1_iiter;
              c1_t_a = c1_s_a + 1;
              c1_iiter = c1_t_a;
              if (c1_mod(chartInstance, c1_iiter) != 0) {
                c1_s.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c1_s.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c1_r2.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c1_r2.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c1_a12.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c1_a12.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c1_a21.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c1_a21.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c1_a22.re = c1_r2.re - c1_s.re;
                c1_a22.im = c1_r2.im - c1_s.im;
                c1_j_a22.re = -c1_a22.re;
                c1_j_a22.im = -c1_a22.im;
                c1_rho = c1_eml_div(chartInstance, c1_j_a22, 2.0);
                c1_b_rho.re = c1_rho.re * c1_rho.re - c1_rho.im * c1_rho.im;
                c1_b_rho.im = c1_rho.re * c1_rho.im + c1_rho.im * c1_rho.re;
                c1_b_a12.re = c1_a12.re * c1_a21.re - c1_a12.im * c1_a21.im;
                c1_b_a12.im = c1_a12.re * c1_a21.im + c1_a12.im * c1_a21.re;
                c1_a22.re = c1_b_rho.re + c1_b_a12.re;
                c1_a22.im = c1_b_rho.im + c1_b_a12.im;
                c1_b_sqrt(chartInstance, &c1_a22);
                c1_a12.re = c1_s.re - (c1_rho.re - c1_a22.re);
                c1_a12.im = c1_s.im - (c1_rho.im - c1_a22.im);
                c1_a21.re = c1_s.re - (c1_rho.re + c1_a22.re);
                c1_a21.im = c1_s.im - (c1_rho.im + c1_a22.im);
                c1_c_a12.re = c1_a12.re - c1_r2.re;
                c1_c_a12.im = c1_a12.im - c1_r2.im;
                c1_b_a21.re = c1_a21.re - c1_r2.re;
                c1_b_a21.im = c1_a21.im - c1_r2.im;
                c1_d7 = c1_abs(chartInstance, c1_c_a12);
                c1_d8 = c1_abs(chartInstance, c1_b_a21);
                if (c1_d7 <= c1_d8) {
                  c1_a21 = c1_a12;
                  c1_rho.re -= c1_a22.re;
                  c1_rho.im -= c1_a22.im;
                } else {
                  c1_rho.re += c1_a22.re;
                  c1_rho.im += c1_a22.im;
                }
              } else {
                c1_eshift.re += c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c1_eshift.im += c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilast), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c1_a21 = c1_eshift;
              }

              c1_b_j = c1_ilastm1;
              c1_u_a = c1_b_j;
              c1_v_a = c1_u_a + 1;
              c1_jp1 = c1_v_a;
              exitg2 = false;
              while ((exitg2 == false) && (c1_b_j > c1_ifirst)) {
                c1_w_a = c1_b_j;
                c1_x_a = c1_w_a - 1;
                c1_jm1 = c1_x_a;
                c1_istart = c1_b_j;
                c1_ctemp.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re - c1_a21.re;
                c1_ctemp.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im - c1_a21.im;
                c1_j_x = c1_ctemp.re;
                c1_k_x = c1_j_x;
                c1_i_y = muDoubleScalarAbs(c1_k_x);
                c1_l_x = c1_ctemp.im;
                c1_m_x = c1_l_x;
                c1_j_y = muDoubleScalarAbs(c1_m_x);
                c1_k_y = c1_i_y + c1_j_y;
                c1_temp = c1_ascale * c1_k_y;
                c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c1_n_x = c1_a22.re;
                c1_o_x = c1_n_x;
                c1_l_y = muDoubleScalarAbs(c1_o_x);
                c1_p_x = c1_a22.im;
                c1_q_x = c1_p_x;
                c1_m_y = muDoubleScalarAbs(c1_q_x);
                c1_n_y = c1_l_y + c1_m_y;
                c1_temp2 = c1_ascale * c1_n_y;
                c1_r_x = c1_temp;
                c1_o_y = c1_temp2;
                c1_tempr = c1_r_x;
                if (c1_o_y > c1_tempr) {
                  c1_tempr = c1_o_y;
                }

                if (c1_tempr < 1.0) {
                  if (c1_tempr != 0.0) {
                    c1_temp /= c1_tempr;
                    c1_temp2 /= c1_tempr;
                  }
                }

                c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c1_s_x = c1_a22.re;
                c1_t_x = c1_s_x;
                c1_p_y = muDoubleScalarAbs(c1_t_x);
                c1_u_x = c1_a22.im;
                c1_v_x = c1_u_x;
                c1_q_y = muDoubleScalarAbs(c1_v_x);
                c1_r_y = c1_p_y + c1_q_y;
                if (c1_r_y * c1_temp2 <= c1_temp * c1_atol) {
                  c1_goto90 = true;
                  exitg2 = true;
                } else {
                  c1_jp1 = c1_b_j;
                  c1_b_j = c1_jm1;
                }
              }

              if (!c1_goto90) {
                c1_istart = c1_ifirst;
                if (c1_istart == c1_ilastm1) {
                  c1_ctemp = c1_rho;
                } else {
                  c1_ctemp.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 2, 0) - 1))
                    - 1].re - c1_a21.re;
                  c1_ctemp.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 2, 0) - 1))
                    - 1].im - c1_a21.im;
                }

                c1_goto90 = true;
              }
            }

            if (c1_goto90) {
              c1_goto90 = false;
              c1_y_a = c1_istart;
              c1_ab_a = c1_y_a;
              c1_g_c = c1_ab_a;
              c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c1_g_c + 1)), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c1_g_c + 1)), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_istart), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c1_b_eml_matlab_zlartg(chartInstance, c1_ctemp, c1_a22, &c1_h_c,
                &c1_s);
              c1_d_c = c1_h_c;
              c1_b_j = c1_istart;
              c1_bb_a = c1_b_j;
              c1_cb_a = c1_bb_a - 1;
              c1_jm1 = c1_cb_a;
              while (c1_b_j < c1_ilast) {
                c1_db_a = c1_b_j;
                c1_eb_a = c1_db_a + 1;
                c1_jp1 = c1_eb_a;
                if (c1_b_j > c1_istart) {
                  c1_e_A.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_e_A.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_f_A.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_f_A.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_eml_matlab_zlartg(chartInstance, c1_e_A, c1_f_A, &c1_i_c,
                                       &c1_s, &c1_a22);
                  c1_d_c = c1_i_c;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) -
                           1)) - 1].re = c1_a22.re;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) -
                           1)) - 1].im = c1_a22.im;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) -
                           1)) - 1].re = c1_dc10.re;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_jp1), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_jm1), 1, 3, 2, 0) -
                           1)) - 1].im = c1_dc10.im;
                }

                c1_j_c = c1_d_c;
                c1_xrow = c1_b_j;
                c1_yrow = c1_jp1;
                c1_jlo = c1_b_j;
                c1_b_jlo = c1_jlo;
                c1_fb_a = c1_b_jlo;
                c1_gb_a = c1_fb_a;
                if (c1_gb_a > 3) {
                  c1_d_overflow = false;
                } else {
                  c1_eml_switch_helper(chartInstance);
                  c1_d_overflow = false;
                }

                if (c1_d_overflow) {
                  c1_check_forloop_overflow_error(chartInstance, c1_d_overflow);
                }

                for (c1_c_j = c1_b_jlo; c1_c_j < 4; c1_c_j++) {
                  c1_d_j = c1_c_j;
                  c1_hb_a = c1_j_c;
                  c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_a12.re = c1_hb_a * c1_a22.re;
                  c1_a12.im = c1_hb_a * c1_a22.im;
                  c1_d_s.re = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re - c1_s.im * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_d_s.im = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im + c1_s.im * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_a21.re = c1_a12.re + c1_d_s.re;
                  c1_a21.im = c1_a12.im + c1_d_s.im;
                  c1_ib_a = c1_j_c;
                  c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_a12.re = c1_ib_a * c1_a22.re;
                  c1_a12.im = c1_ib_a * c1_a22.im;
                  c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c1_k_a22 = c1_a22;
                  c1_l_a22 = c1_a22;
                  c1_m_a22 = c1_a22;
                  c1_n_a22 = c1_a22;
                  c1_a22.re = c1_s.re * c1_k_a22.re + c1_s.im * c1_l_a22.im;
                  c1_a22.im = c1_s.re * c1_m_a22.im - c1_s.im * c1_n_a22.re;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) -
                             1)) - 1].re = c1_a12.re - c1_a22.re;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_yrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) -
                             1)) - 1].im = c1_a12.im - c1_a22.im;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) -
                             1)) - 1].re = c1_a21.re;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_xrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_d_j), 1, 3, 2, 0) -
                             1)) - 1].im = c1_a21.im;
                }

                c1_s.re = -c1_s.re;
                c1_s.im = -c1_s.im;
                c1_jb_a = c1_jp1;
                c1_kb_a = c1_jb_a;
                c1_k_c = c1_kb_a;
                c1_w_x = c1_k_c + 1;
                c1_s_y = c1_ilast;
                c1_x_x = c1_w_x;
                if (c1_s_y < c1_x_x) {
                  c1_x_x = c1_s_y;
                }

                c1_l_c = c1_d_c;
                c1_c_xcol = c1_jp1;
                c1_c_ycol = c1_b_j;
                c1_d_ihi = c1_x_x;
                c1_e_ihi = c1_d_ihi;
                c1_i_b = c1_e_ihi;
                c1_j_b = c1_i_b;
                if (1 > c1_j_b) {
                  c1_e_overflow = false;
                } else {
                  c1_eml_switch_helper(chartInstance);
                  c1_e_overflow = (c1_j_b > 2147483646);
                }

                if (c1_e_overflow) {
                  c1_check_forloop_overflow_error(chartInstance, c1_e_overflow);
                }

                for (c1_e_i = 1; c1_e_i <= c1_e_ihi; c1_e_i++) {
                  c1_f_i = c1_e_i;
                  c1_lb_a = c1_l_c;
                  c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_a12.re = c1_lb_a * c1_a22.re;
                  c1_a12.im = c1_lb_a * c1_a22.im;
                  c1_e_s.re = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].re - c1_s.im * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_e_s.im = c1_s.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].im + c1_s.im * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a21.re = c1_a12.re + c1_e_s.re;
                  c1_a21.im = c1_a12.im + c1_e_s.im;
                  c1_mb_a = c1_l_c;
                  c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_a12.re = c1_mb_a * c1_a22.re;
                  c1_a12.im = c1_mb_a * c1_a22.im;
                  c1_a22.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a22.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_o_a22 = c1_a22;
                  c1_p_a22 = c1_a22;
                  c1_q_a22 = c1_a22;
                  c1_r_a22 = c1_a22;
                  c1_a22.re = c1_s.re * c1_o_a22.re + c1_s.im * c1_p_a22.im;
                  c1_a22.im = c1_s.re * c1_q_a22.im - c1_s.im * c1_r_a22.re;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_c_ycol), 1, 3, 2, 0)
                           - 1)) - 1].re = c1_a12.re - c1_a22.re;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_c_ycol), 1, 3, 2, 0)
                           - 1)) - 1].im = c1_a12.im - c1_a22.im;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_c_xcol), 1, 3, 2, 0)
                           - 1)) - 1].re = c1_a21.re;
                  c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_c_xcol), 1, 3, 2, 0)
                           - 1)) - 1].im = c1_a21.im;
                }

                c1_m_c = c1_d_c;
                c1_d_xcol = c1_jp1;
                c1_d_ycol = c1_b_j;
                c1_eml_switch_helper(chartInstance);
                for (c1_g_i = 1; c1_g_i < 4; c1_g_i++) {
                  c1_h_i = c1_g_i;
                  c1_nb_a = c1_m_c;
                  c1_a22.re = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a22.im = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_a12.re = c1_nb_a * c1_a22.re;
                  c1_a12.im = c1_nb_a * c1_a22.im;
                  c1_f_s.re = c1_s.re * c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].re - c1_s.im * c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_f_s.im = c1_s.re * c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].im + c1_s.im * c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a21.re = c1_a12.re + c1_f_s.re;
                  c1_a21.im = c1_a12.im + c1_f_s.im;
                  c1_ob_a = c1_m_c;
                  c1_a22.re = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a22.im = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_a12.re = c1_ob_a * c1_a22.re;
                  c1_a12.im = c1_ob_a * c1_a22.im;
                  c1_a22.re = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c1_a22.im = c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c1_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c1_s_a22 = c1_a22;
                  c1_t_a22 = c1_a22;
                  c1_u_a22 = c1_a22;
                  c1_v_a22 = c1_a22;
                  c1_a22.re = c1_s.re * c1_s_a22.re + c1_s.im * c1_t_a22.im;
                  c1_a22.im = c1_s.re * c1_u_a22.im - c1_s.im * c1_v_a22.re;
                  c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_d_ycol), 1, 3, 2, 0)
                           - 1)) - 1].re = c1_a12.re - c1_a22.re;
                  c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_d_ycol), 1, 3, 2, 0)
                           - 1)) - 1].im = c1_a12.im - c1_a22.im;
                  c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_d_xcol), 1, 3, 2, 0)
                           - 1)) - 1].re = c1_a21.re;
                  c1_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c1_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c1_d_xcol), 1, 3, 2, 0)
                           - 1)) - 1].im = c1_a21.im;
                }

                c1_jm1 = c1_b_j;
                c1_b_j = c1_jp1;
              }
            }

            guard11 = true;
          }

          if (guard11 == true) {
            c1_jiter++;
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
    if (c1_failed) {
      *c1_info = (real_T)c1_ilast;
      c1_b_ilast = c1_ilast;
      c1_k_b = c1_b_ilast;
      c1_l_b = c1_k_b;
      if (1 > c1_l_b) {
        c1_f_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_f_overflow = (c1_l_b > 2147483646);
      }

      if (c1_f_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_f_overflow);
      }

      for (c1_k = 1; c1_k <= c1_b_ilast; c1_k++) {
        c1_b_k = c1_k;
        c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_k), 1, 3, 1, 0) - 1].re = c1_dc9.re;
        c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_k), 1, 3, 1, 0) - 1].im = c1_dc9.im;
        c1_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_k), 1, 3, 1, 0) - 1].re = c1_dc9.re;
        c1_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_k), 1, 3, 1, 0) - 1].im = c1_dc9.im;
      }

      for (c1_i134 = 0; c1_i134 < 9; c1_i134++) {
        c1_Z[c1_i134] = c1_dc9;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1 == true) {
    c1_pb_a = c1_ilo;
    c1_qb_a = c1_pb_a - 1;
    c1_i135 = c1_qb_a;
    c1_m_b = c1_i135;
    c1_n_b = c1_m_b;
    if (1 > c1_n_b) {
      c1_g_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_g_overflow = (c1_n_b > 2147483646);
    }

    if (c1_g_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_g_overflow);
    }

    for (c1_e_j = 1; c1_e_j <= c1_i135; c1_e_j++) {
      c1_b_j = c1_e_j;
      c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_j), 1, 3, 1, 0) - 1].re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
        ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
        c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c1_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_j), 1, 3, 1, 0) - 1].im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
        ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 3, 1, 0) + 3 *
        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
        c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
    }

    *c1_info = 0.0;
  }
}

static void c1_b_eml_matlab_ztgevc(SFc1_Model_01InstanceStruct *chartInstance,
  creal_T c1_A[9], creal_T c1_V[9])
{
  int32_T c1_i136;
  creal_T c1_work1[3];
  int32_T c1_i137;
  creal_T c1_work2[3];
  int32_T c1_i138;
  real_T c1_rworka[3];
  creal_T c1_ca;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_y;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_b_y;
  real_T c1_anorm;
  int32_T c1_j;
  real_T c1_b_j;
  real_T c1_d9;
  int32_T c1_i139;
  int32_T c1_i;
  real_T c1_b_i;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_c_y;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_d_y;
  real_T c1_e_y;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_f_y;
  real_T c1_k_x;
  real_T c1_l_x;
  real_T c1_g_y;
  real_T c1_h_y;
  real_T c1_i_y;
  real_T c1_m_x;
  real_T c1_n_x;
  real_T c1_ascale;
  int32_T c1_je;
  real_T c1_b_je;
  real_T c1_ieig;
  real_T c1_o_x;
  real_T c1_p_x;
  real_T c1_j_y;
  real_T c1_q_x;
  real_T c1_r_x;
  real_T c1_k_y;
  real_T c1_l_y;
  real_T c1_s_x;
  real_T c1_t_x;
  real_T c1_temp;
  real_T c1_sbeta;
  real_T c1_a;
  real_T c1_b;
  creal_T c1_salpha;
  real_T c1_acoeff;
  boolean_T c1_b6;
  boolean_T c1_lscalea;
  real_T c1_u_x;
  real_T c1_v_x;
  real_T c1_m_y;
  real_T c1_w_x;
  real_T c1_x_x;
  real_T c1_n_y;
  real_T c1_o_y;
  real_T c1_y_x;
  real_T c1_ab_x;
  real_T c1_p_y;
  real_T c1_bb_x;
  real_T c1_cb_x;
  real_T c1_q_y;
  real_T c1_r_y;
  boolean_T c1_b7;
  boolean_T c1_lscaleb;
  real_T c1_scale;
  real_T c1_db_x;
  real_T c1_eb_x;
  real_T c1_fb_x;
  real_T c1_gb_x;
  real_T c1_s_y;
  real_T c1_hb_x;
  real_T c1_ib_x;
  real_T c1_t_y;
  real_T c1_u_y;
  real_T c1_v_y;
  real_T c1_jb_x;
  real_T c1_kb_x;
  real_T c1_w_y;
  real_T c1_lb_x;
  real_T c1_mb_x;
  real_T c1_x_y;
  real_T c1_y_y;
  real_T c1_nb_x;
  real_T c1_z;
  real_T c1_ob_x;
  real_T c1_ab_y;
  real_T c1_b_a;
  real_T c1_c_a;
  real_T c1_acoefa;
  real_T c1_pb_x;
  real_T c1_qb_x;
  real_T c1_bb_y;
  real_T c1_rb_x;
  real_T c1_sb_x;
  real_T c1_cb_y;
  real_T c1_bcoefa;
  int32_T c1_jr;
  real_T c1_b_jr;
  static creal_T c1_dc11 = { 0.0, 0.0 };

  static creal_T c1_dc12 = { 1.0, 0.0 };

  real_T c1_tb_x;
  real_T c1_db_y;
  real_T c1_dmin;
  real_T c1_d10;
  int32_T c1_i140;
  int32_T c1_c_jr;
  real_T c1_d_a;
  real_T c1_d11;
  int32_T c1_i141;
  int32_T c1_c_j;
  real_T c1_e_a;
  creal_T c1_d;
  real_T c1_ub_x;
  real_T c1_vb_x;
  real_T c1_eb_y;
  real_T c1_wb_x;
  real_T c1_xb_x;
  real_T c1_fb_y;
  real_T c1_gb_y;
  real_T c1_yb_x;
  real_T c1_ac_x;
  real_T c1_hb_y;
  real_T c1_bc_x;
  real_T c1_cc_x;
  real_T c1_ib_y;
  real_T c1_jb_y;
  real_T c1_dc_x;
  real_T c1_ec_x;
  real_T c1_kb_y;
  real_T c1_fc_x;
  real_T c1_gc_x;
  real_T c1_lb_y;
  real_T c1_mb_y;
  real_T c1_hc_x;
  real_T c1_ic_x;
  real_T c1_nb_y;
  real_T c1_jc_x;
  real_T c1_kc_x;
  real_T c1_ob_y;
  real_T c1_pb_y;
  real_T c1_lc_x;
  real_T c1_mc_x;
  real_T c1_qb_y;
  real_T c1_nc_x;
  real_T c1_oc_x;
  real_T c1_rb_y;
  real_T c1_sb_y;
  real_T c1_c_je;
  int32_T c1_i142;
  int32_T c1_d_jr;
  real_T c1_f_a;
  real_T c1_pc_x;
  real_T c1_qc_x;
  real_T c1_tb_y;
  real_T c1_rc_x;
  real_T c1_sc_x;
  real_T c1_ub_y;
  real_T c1_vb_y;
  real_T c1_tc_x;
  real_T c1_uc_x;
  real_T c1_wb_y;
  real_T c1_vc_x;
  real_T c1_wc_x;
  real_T c1_xb_y;
  real_T c1_yb_y;
  real_T c1_d_je;
  int32_T c1_i143;
  int32_T c1_e_jr;
  real_T c1_g_a;
  real_T c1_h_a;
  real_T c1_d12;
  int32_T c1_i144;
  int32_T c1_f_jr;
  creal_T c1_b_ca;
  int32_T c1_g_jr;
  real_T c1_e_je;
  int32_T c1_i145;
  int32_T c1_jc;
  real_T c1_b_jc;
  int32_T c1_h_jr;
  creal_T c1_b_V;
  real_T c1_xc_x;
  real_T c1_yc_x;
  real_T c1_ac_y;
  real_T c1_ad_x;
  real_T c1_bd_x;
  real_T c1_bc_y;
  real_T c1_xmx;
  int32_T c1_i_jr;
  real_T c1_cd_x;
  real_T c1_dd_x;
  real_T c1_cc_y;
  real_T c1_ed_x;
  real_T c1_fd_x;
  real_T c1_dc_y;
  real_T c1_ec_y;
  real_T c1_fc_y;
  int32_T c1_j_jr;
  real_T c1_i_a;
  int32_T c1_k_jr;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  for (c1_i136 = 0; c1_i136 < 3; c1_i136++) {
    c1_work1[c1_i136].re = 0.0;
    c1_work1[c1_i136].im = 0.0;
  }

  for (c1_i137 = 0; c1_i137 < 3; c1_i137++) {
    c1_work2[c1_i137].re = 0.0;
    c1_work2[c1_i137].im = 0.0;
  }

  c1_eps(chartInstance);
  c1_realmin(chartInstance);
  for (c1_i138 = 0; c1_i138 < 3; c1_i138++) {
    c1_rworka[c1_i138] = 0.0;
  }

  c1_ca = c1_A[0];
  c1_x = c1_ca.re;
  c1_b_x = c1_x;
  c1_y = muDoubleScalarAbs(c1_b_x);
  c1_c_x = c1_ca.im;
  c1_d_x = c1_c_x;
  c1_b_y = muDoubleScalarAbs(c1_d_x);
  c1_anorm = c1_y + c1_b_y;
  for (c1_j = 0; c1_j < 2; c1_j++) {
    c1_b_j = 2.0 + (real_T)c1_j;
    c1_d9 = c1_b_j - 1.0;
    c1_i139 = (int32_T)c1_d9;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c1_d9, mxDOUBLE_CLASS, c1_i139);
    for (c1_i = 0; c1_i < c1_i139; c1_i++) {
      c1_b_i = 1.0 + (real_T)c1_i;
      c1_ca.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_i), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c1_ca.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_i), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
      c1_e_x = c1_ca.re;
      c1_f_x = c1_e_x;
      c1_c_y = muDoubleScalarAbs(c1_f_x);
      c1_g_x = c1_ca.im;
      c1_h_x = c1_g_x;
      c1_d_y = muDoubleScalarAbs(c1_h_x);
      c1_e_y = c1_c_y + c1_d_y;
      c1_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_j), 1, 3, 1, 0) - 1] = c1_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)_SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1] + c1_e_y;
    }

    c1_ca.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_j), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
    c1_ca.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_j), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
    c1_i_x = c1_ca.re;
    c1_j_x = c1_i_x;
    c1_f_y = muDoubleScalarAbs(c1_j_x);
    c1_k_x = c1_ca.im;
    c1_l_x = c1_k_x;
    c1_g_y = muDoubleScalarAbs(c1_l_x);
    c1_h_y = c1_f_y + c1_g_y;
    c1_i_y = c1_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1] + c1_h_y;
    if (c1_i_y > c1_anorm) {
      c1_anorm = c1_i_y;
    }
  }

  c1_m_x = c1_anorm;
  c1_n_x = c1_m_x;
  if (2.2250738585072014E-308 > c1_n_x) {
    c1_n_x = 2.2250738585072014E-308;
  }

  c1_ascale = 1.0 / c1_n_x;
  for (c1_je = 0; c1_je < 3; c1_je++) {
    c1_b_je = 3.0 + -(real_T)c1_je;
    c1_ieig = c1_b_je;
    c1_ca.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c1_b_je), 1, 3, 2, 0) - 1)) - 1].re;
    c1_ca.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c1_b_je), 1, 3, 2, 0) - 1)) - 1].im;
    c1_o_x = c1_ca.re;
    c1_p_x = c1_o_x;
    c1_j_y = muDoubleScalarAbs(c1_p_x);
    c1_q_x = c1_ca.im;
    c1_r_x = c1_q_x;
    c1_k_y = muDoubleScalarAbs(c1_r_x);
    c1_l_y = c1_j_y + c1_k_y;
    c1_s_x = c1_l_y * c1_ascale;
    c1_t_x = c1_s_x;
    if (1.0 > c1_t_x) {
      c1_t_x = 1.0;
    }

    c1_temp = 1.0 / c1_t_x;
    c1_sbeta = c1_temp;
    c1_a = c1_temp;
    c1_ca.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c1_b_je), 1, 3, 2, 0) - 1)) - 1].re;
    c1_ca.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c1_b_je), 1, 3, 2, 0) - 1)) - 1].im;
    c1_ca.re *= c1_a;
    c1_ca.im *= c1_a;
    c1_b = c1_ascale;
    c1_salpha.re = c1_b * c1_ca.re;
    c1_salpha.im = c1_b * c1_ca.im;
    c1_acoeff = c1_sbeta * c1_ascale;
    guard3 = false;
    if (c1_b_abs(chartInstance, c1_sbeta) >= 2.2250738585072014E-308) {
      if (c1_b_abs(chartInstance, c1_acoeff) < 3.0062525400134592E-292) {
        c1_b6 = true;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }

    if (guard3 == true) {
      c1_b6 = false;
    }

    c1_lscalea = c1_b6;
    c1_u_x = c1_salpha.re;
    c1_v_x = c1_u_x;
    c1_m_y = muDoubleScalarAbs(c1_v_x);
    c1_w_x = c1_salpha.im;
    c1_x_x = c1_w_x;
    c1_n_y = muDoubleScalarAbs(c1_x_x);
    c1_o_y = c1_m_y + c1_n_y;
    guard2 = false;
    if (c1_o_y >= 2.2250738585072014E-308) {
      c1_y_x = c1_salpha.re;
      c1_ab_x = c1_y_x;
      c1_p_y = muDoubleScalarAbs(c1_ab_x);
      c1_bb_x = c1_salpha.im;
      c1_cb_x = c1_bb_x;
      c1_q_y = muDoubleScalarAbs(c1_cb_x);
      c1_r_y = c1_p_y + c1_q_y;
      if (c1_r_y < 3.0062525400134592E-292) {
        c1_b7 = true;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2 == true) {
      c1_b7 = false;
    }

    c1_lscaleb = c1_b7;
    c1_scale = 1.0;
    if (c1_lscalea) {
      c1_db_x = c1_anorm;
      c1_eb_x = c1_db_x;
      if (3.3264005158911995E+291 < c1_eb_x) {
        c1_eb_x = 3.3264005158911995E+291;
      }

      c1_scale = 3.0062525400134592E-292 / c1_b_abs(chartInstance, c1_sbeta) *
        c1_eb_x;
    }

    if (c1_lscaleb) {
      c1_fb_x = c1_salpha.re;
      c1_gb_x = c1_fb_x;
      c1_s_y = muDoubleScalarAbs(c1_gb_x);
      c1_hb_x = c1_salpha.im;
      c1_ib_x = c1_hb_x;
      c1_t_y = muDoubleScalarAbs(c1_ib_x);
      c1_u_y = c1_s_y + c1_t_y;
      c1_v_y = 3.0062525400134592E-292 / c1_u_y;
      if (c1_v_y > c1_scale) {
        c1_scale = c1_v_y;
      }
    }

    guard1 = false;
    if (c1_lscalea) {
      guard1 = true;
    } else {
      if (c1_lscaleb) {
        guard1 = true;
      }
    }

    if (guard1 == true) {
      c1_jb_x = c1_salpha.re;
      c1_kb_x = c1_jb_x;
      c1_w_y = muDoubleScalarAbs(c1_kb_x);
      c1_lb_x = c1_salpha.im;
      c1_mb_x = c1_lb_x;
      c1_x_y = muDoubleScalarAbs(c1_mb_x);
      c1_y_y = c1_w_y + c1_x_y;
      c1_nb_x = c1_b_abs(chartInstance, c1_acoeff);
      c1_z = c1_y_y;
      c1_ob_x = c1_nb_x;
      if (1.0 > c1_ob_x) {
        c1_ob_x = 1.0;
      }

      if (c1_z > c1_ob_x) {
        c1_ob_x = c1_z;
      }

      c1_ab_y = 1.0 / (2.2250738585072014E-308 * c1_ob_x);
      if (c1_ab_y < c1_scale) {
        c1_scale = c1_ab_y;
      }

      if (c1_lscalea) {
        c1_acoeff = c1_ascale * (c1_scale * c1_sbeta);
      } else {
        c1_acoeff *= c1_scale;
      }

      if (c1_lscaleb) {
        c1_b_a = c1_scale;
        c1_salpha.re *= c1_b_a;
        c1_salpha.im *= c1_b_a;
      } else {
        c1_c_a = c1_scale;
        c1_salpha.re *= c1_c_a;
        c1_salpha.im *= c1_c_a;
      }
    }

    c1_acoefa = c1_b_abs(chartInstance, c1_acoeff);
    c1_pb_x = c1_salpha.re;
    c1_qb_x = c1_pb_x;
    c1_bb_y = muDoubleScalarAbs(c1_qb_x);
    c1_rb_x = c1_salpha.im;
    c1_sb_x = c1_rb_x;
    c1_cb_y = muDoubleScalarAbs(c1_sb_x);
    c1_bcoefa = c1_bb_y + c1_cb_y;
    for (c1_jr = 0; c1_jr < 3; c1_jr++) {
      c1_b_jr = 1.0 + (real_T)c1_jr;
      c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_jr), 1, 3, 1, 0) - 1].re = c1_dc11.re;
      c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_jr), 1, 3, 1, 0) - 1].im = c1_dc11.im;
    }

    c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c1_b_je), 1, 3, 1, 0) - 1].re = c1_dc12.re;
    c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c1_b_je), 1, 3, 1, 0) - 1].im = c1_dc12.im;
    c1_tb_x = 2.2204460492503131E-16 * c1_acoefa * c1_anorm;
    c1_db_y = 2.2204460492503131E-16 * c1_bcoefa;
    c1_dmin = c1_tb_x;
    if (c1_db_y > c1_dmin) {
      c1_dmin = c1_db_y;
    }

    if (2.2250738585072014E-308 > c1_dmin) {
      c1_dmin = 2.2250738585072014E-308;
    }

    c1_d10 = c1_b_je - 1.0;
    c1_i140 = (int32_T)c1_d10;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c1_d10, mxDOUBLE_CLASS, c1_i140);
    for (c1_c_jr = 0; c1_c_jr < c1_i140; c1_c_jr++) {
      c1_b_jr = 1.0 + (real_T)c1_c_jr;
      c1_d_a = c1_acoeff;
      c1_ca.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_je), 1, 3, 2, 0) - 1)) - 1].re;
      c1_ca.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_je), 1, 3, 2, 0) - 1)) - 1].im;
      c1_ca.re *= c1_d_a;
      c1_ca.im *= c1_d_a;
      c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_jr), 1, 3, 1, 0) - 1].re = c1_ca.re;
      c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_jr), 1, 3, 1, 0) - 1].im = c1_ca.im;
    }

    c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c1_b_je), 1, 3, 1, 0) - 1].re = c1_dc12.re;
    c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c1_b_je), 1, 3, 1, 0) - 1].im = c1_dc12.im;
    c1_d11 = c1_b_je - 1.0;
    c1_i141 = (int32_T)-(1.0 + (-1.0 - c1_d11));
    _SFD_FOR_LOOP_VECTOR_CHECK(c1_d11, -1.0, 1.0, mxDOUBLE_CLASS, c1_i141);
    for (c1_c_j = 0; c1_c_j < c1_i141; c1_c_j++) {
      c1_b_j = c1_d11 + -(real_T)c1_c_j;
      c1_e_a = c1_acoeff;
      c1_ca.re = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c1_ca.im = c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].im;
      c1_ca.re *= c1_e_a;
      c1_ca.im *= c1_e_a;
      c1_d.re = c1_ca.re - c1_salpha.re;
      c1_d.im = c1_ca.im - c1_salpha.im;
      c1_ub_x = c1_d.re;
      c1_vb_x = c1_ub_x;
      c1_eb_y = muDoubleScalarAbs(c1_vb_x);
      c1_wb_x = c1_d.im;
      c1_xb_x = c1_wb_x;
      c1_fb_y = muDoubleScalarAbs(c1_xb_x);
      c1_gb_y = c1_eb_y + c1_fb_y;
      if (c1_gb_y <= c1_dmin) {
        c1_d.re = c1_dmin;
        c1_d.im = 0.0;
      }

      c1_yb_x = c1_d.re;
      c1_ac_x = c1_yb_x;
      c1_hb_y = muDoubleScalarAbs(c1_ac_x);
      c1_bc_x = c1_d.im;
      c1_cc_x = c1_bc_x;
      c1_ib_y = muDoubleScalarAbs(c1_cc_x);
      c1_jb_y = c1_hb_y + c1_ib_y;
      if (c1_jb_y < 1.0) {
        c1_ca.re = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].re;
        c1_ca.im = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].im;
        c1_dc_x = c1_ca.re;
        c1_ec_x = c1_dc_x;
        c1_kb_y = muDoubleScalarAbs(c1_ec_x);
        c1_fc_x = c1_ca.im;
        c1_gc_x = c1_fc_x;
        c1_lb_y = muDoubleScalarAbs(c1_gc_x);
        c1_mb_y = c1_kb_y + c1_lb_y;
        c1_hc_x = c1_d.re;
        c1_ic_x = c1_hc_x;
        c1_nb_y = muDoubleScalarAbs(c1_ic_x);
        c1_jc_x = c1_d.im;
        c1_kc_x = c1_jc_x;
        c1_ob_y = muDoubleScalarAbs(c1_kc_x);
        c1_pb_y = c1_nb_y + c1_ob_y;
        if (c1_mb_y >= 1.4980776123852632E+307 * c1_pb_y) {
          c1_ca.re = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].re;
          c1_ca.im = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].im;
          c1_lc_x = c1_ca.re;
          c1_mc_x = c1_lc_x;
          c1_qb_y = muDoubleScalarAbs(c1_mc_x);
          c1_nc_x = c1_ca.im;
          c1_oc_x = c1_nc_x;
          c1_rb_y = muDoubleScalarAbs(c1_oc_x);
          c1_sb_y = c1_qb_y + c1_rb_y;
          c1_temp = 1.0 / c1_sb_y;
          c1_c_je = c1_b_je;
          c1_i142 = (int32_T)c1_c_je;
          _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c1_c_je, mxDOUBLE_CLASS, c1_i142);
          for (c1_d_jr = 0; c1_d_jr < c1_i142; c1_d_jr++) {
            c1_b_jr = 1.0 + (real_T)c1_d_jr;
            c1_f_a = c1_temp;
            c1_ca.re = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].re;
            c1_ca.im = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].im;
            c1_ca.re *= c1_f_a;
            c1_ca.im *= c1_f_a;
            c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
              ("", c1_b_jr), 1, 3, 1, 0) - 1].re = c1_ca.re;
            c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
              ("", c1_b_jr), 1, 3, 1, 0) - 1].im = c1_ca.im;
          }
        }
      }

      c1_ca.re = -c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].re;
      c1_ca.im = -c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].im;
      c1_ca = c1_rdivide(chartInstance, c1_ca, c1_d);
      c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_j), 1, 3, 1, 0) - 1].re = c1_ca.re;
      c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_j), 1, 3, 1, 0) - 1].im = c1_ca.im;
      if (c1_b_j > 1.0) {
        c1_ca.re = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].re;
        c1_ca.im = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].im;
        c1_pc_x = c1_ca.re;
        c1_qc_x = c1_pc_x;
        c1_tb_y = muDoubleScalarAbs(c1_qc_x);
        c1_rc_x = c1_ca.im;
        c1_sc_x = c1_rc_x;
        c1_ub_y = muDoubleScalarAbs(c1_sc_x);
        c1_vb_y = c1_tb_y + c1_ub_y;
        if (c1_vb_y > 1.0) {
          c1_ca.re = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].re;
          c1_ca.im = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].im;
          c1_tc_x = c1_ca.re;
          c1_uc_x = c1_tc_x;
          c1_wb_y = muDoubleScalarAbs(c1_uc_x);
          c1_vc_x = c1_ca.im;
          c1_wc_x = c1_vc_x;
          c1_xb_y = muDoubleScalarAbs(c1_wc_x);
          c1_yb_y = c1_wb_y + c1_xb_y;
          c1_temp = 1.0 / c1_yb_y;
          if (c1_acoefa * c1_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
               _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1] >=
              1.4980776123852632E+307 * c1_temp) {
            c1_d_je = c1_b_je;
            c1_i143 = (int32_T)c1_d_je;
            _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c1_d_je, mxDOUBLE_CLASS,
              c1_i143);
            for (c1_e_jr = 0; c1_e_jr < c1_i143; c1_e_jr++) {
              c1_b_jr = 1.0 + (real_T)c1_e_jr;
              c1_g_a = c1_temp;
              c1_ca.re = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].re;
              c1_ca.im = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].im;
              c1_ca.re *= c1_g_a;
              c1_ca.im *= c1_g_a;
              c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].re = c1_ca.re;
              c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].im = c1_ca.im;
            }
          }
        }

        c1_h_a = c1_acoeff;
        c1_ca.re = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].re;
        c1_ca.im = c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 1, 0) - 1].im;
        c1_ca.re *= c1_h_a;
        c1_ca.im *= c1_h_a;
        c1_d12 = c1_b_j - 1.0;
        c1_i144 = (int32_T)c1_d12;
        _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c1_d12, mxDOUBLE_CLASS, c1_i144);
        for (c1_f_jr = 0; c1_f_jr < c1_i144; c1_f_jr++) {
          c1_b_jr = 1.0 + (real_T)c1_f_jr;
          c1_b_ca.re = c1_ca.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            c1_b_j), 1, 3, 2, 0) - 1)) - 1].re - c1_ca.im * c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c1_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].
            im;
          c1_b_ca.im = c1_ca.re * c1_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            c1_b_j), 1, 3, 2, 0) - 1)) - 1].im + c1_ca.im * c1_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c1_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1].
            re;
          c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
            "", c1_b_jr), 1, 3, 1, 0) - 1].re =
            c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c1_b_jr), 1, 3, 1, 0) - 1].re + c1_b_ca.re;
          c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
            "", c1_b_jr), 1, 3, 1, 0) - 1].im =
            c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c1_b_jr), 1, 3, 1, 0) - 1].im + c1_b_ca.im;
        }
      }
    }

    for (c1_g_jr = 0; c1_g_jr < 3; c1_g_jr++) {
      c1_b_jr = 1.0 + (real_T)c1_g_jr;
      c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_jr), 1, 3, 1, 0) - 1].re = c1_dc11.re;
      c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_jr), 1, 3, 1, 0) - 1].im = c1_dc11.im;
    }

    c1_e_je = c1_b_je;
    c1_i145 = (int32_T)c1_e_je;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c1_e_je, mxDOUBLE_CLASS, c1_i145);
    for (c1_jc = 0; c1_jc < c1_i145; c1_jc++) {
      c1_b_jc = 1.0 + (real_T)c1_jc;
      for (c1_h_jr = 0; c1_h_jr < 3; c1_h_jr++) {
        c1_b_jr = 1.0 + (real_T)c1_h_jr;
        c1_b_V.re = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1)) - 1].re *
          c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c1_b_jc), 1, 3, 1, 0) - 1].re - c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c1_b_jc), 1, 3, 2, 0) - 1)) - 1].im *
          c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c1_b_jc), 1, 3, 1, 0) - 1].im;
        c1_b_V.im = c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_jc), 1, 3, 2, 0) - 1)) - 1].re *
          c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c1_b_jc), 1, 3, 1, 0) - 1].im + c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) + 3 *
          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c1_b_jc), 1, 3, 2, 0) - 1)) - 1].im *
          c1_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c1_b_jc), 1, 3, 1, 0) - 1].re;
        c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c1_b_jr), 1, 3, 1, 0) - 1].re = c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].re +
          c1_b_V.re;
        c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c1_b_jr), 1, 3, 1, 0) - 1].im = c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].im +
          c1_b_V.im;
      }
    }

    c1_ca = c1_work2[0];
    c1_xc_x = c1_ca.re;
    c1_yc_x = c1_xc_x;
    c1_ac_y = muDoubleScalarAbs(c1_yc_x);
    c1_ad_x = c1_ca.im;
    c1_bd_x = c1_ad_x;
    c1_bc_y = muDoubleScalarAbs(c1_bd_x);
    c1_xmx = c1_ac_y + c1_bc_y;
    for (c1_i_jr = 0; c1_i_jr < 2; c1_i_jr++) {
      c1_b_jr = 2.0 + (real_T)c1_i_jr;
      c1_ca.re = c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].re;
      c1_ca.im = c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].im;
      c1_cd_x = c1_ca.re;
      c1_dd_x = c1_cd_x;
      c1_cc_y = muDoubleScalarAbs(c1_dd_x);
      c1_ed_x = c1_ca.im;
      c1_fd_x = c1_ed_x;
      c1_dc_y = muDoubleScalarAbs(c1_fd_x);
      c1_ec_y = c1_cc_y + c1_dc_y;
      c1_fc_y = c1_ec_y;
      if (c1_fc_y > c1_xmx) {
        c1_xmx = c1_fc_y;
      }
    }

    if (c1_xmx > 2.2250738585072014E-308) {
      c1_temp = 1.0 / c1_xmx;
      for (c1_j_jr = 0; c1_j_jr < 3; c1_j_jr++) {
        c1_b_jr = 1.0 + (real_T)c1_j_jr;
        c1_i_a = c1_temp;
        c1_ca.re = c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].re;
        c1_ca.im = c1_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c1_b_jr), 1, 3, 1, 0) - 1].im;
        c1_ca.re *= c1_i_a;
        c1_ca.im *= c1_i_a;
        c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c1_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c1_ieig), 1, 3, 2, 0) - 1)) - 1]
          .re = c1_ca.re;
        c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c1_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c1_ieig), 1, 3, 2, 0) - 1)) - 1]
          .im = c1_ca.im;
      }
    } else {
      for (c1_k_jr = 0; c1_k_jr < 3; c1_k_jr++) {
        c1_b_jr = 1.0 + (real_T)c1_k_jr;
        c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c1_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c1_ieig), 1, 3, 2, 0) - 1)) - 1]
          .re = c1_dc11.re;
        c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c1_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c1_ieig), 1, 3, 2, 0) - 1)) - 1]
          .im = c1_dc11.im;
      }
    }
  }
}

static int32_T c1_div_s32(SFc1_Model_01InstanceStruct *chartInstance, int32_T
  c1_numerator, int32_T c1_denominator)
{
  int32_T c1_quotient;
  uint32_T c1_absNumerator;
  uint32_T c1_absDenominator;
  boolean_T c1_quotientNeedsNegation;
  uint32_T c1_tempAbsQuotient;
  (void)chartInstance;
  if (c1_denominator == 0) {
    if (c1_numerator >= 0) {
      c1_quotient = MAX_int32_T;
    } else {
      c1_quotient = MIN_int32_T;
    }

    _SFD_OVERFLOW_DETECTION(SFDB_DIVIDE_BY_ZERO);
  } else {
    if (c1_numerator >= 0) {
      c1_absNumerator = (uint32_T)c1_numerator;
    } else {
      c1_absNumerator = (uint32_T)-c1_numerator;
    }

    if (c1_denominator >= 0) {
      c1_absDenominator = (uint32_T)c1_denominator;
    } else {
      c1_absDenominator = (uint32_T)-c1_denominator;
    }

    c1_quotientNeedsNegation = (c1_numerator < 0 != c1_denominator < 0);
    c1_tempAbsQuotient = c1_absNumerator / c1_absDenominator;
    if (c1_quotientNeedsNegation) {
      c1_quotient = -(int32_T)c1_tempAbsQuotient;
    } else {
      c1_quotient = (int32_T)c1_tempAbsQuotient;
    }
  }

  return c1_quotient;
}

static void init_dsm_address_info(SFc1_Model_01InstanceStruct *chartInstance)
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

void sf_c1_Model_01_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(854522139U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2183244584U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(621910901U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(442532969U);
}

mxArray *sf_c1_Model_01_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("mh8d8Qhhhm8NVF2I6X4oeH");
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
      pr[0] = (double)(3);
      pr[1] = (double)(3);
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

mxArray *sf_c1_Model_01_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c1_Model_01_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c1_Model_01(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"Phi_r12\",},{M[8],M[0],T\"is_active_c1_Model_01\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_Model_01_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_Model_01InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc1_Model_01InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_01MachineNumber_,
           1,
           1,
           1,
           0,
           4,
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
          _SFD_SET_DATA_PROPS(1,2,0,1,"Phi_r12");
          _SFD_SET_DATA_PROPS(2,1,1,0,"M");
          _SFD_SET_DATA_PROPS(3,1,1,0,"t_delta");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,612);
        _SFD_CV_INIT_EML_FCN(0,1,"fn_phi_jk",617,-1,955);
        _SFD_CV_INIT_EML_FCN(0,2,"fn_phi_j1",997,-1,1100);
        _SFD_CV_INIT_EML_FCN(0,3,"fn_phi_j2",1101,-1,1320);
        _SFD_CV_INIT_EML_FCN(0,4,"fn_phi_j3",1321,-1,1593);
        _SFD_CV_INIT_EML_FOR(0,1,0,371,381,581);
        _SFD_CV_INIT_EML_FOR(0,1,1,389,399,573);
        _SFD_CV_INIT_EML_FOR(0,1,2,411,422,561);

        {
          static int caseStart[] = { 910, 696, 759, 834 };

          static int caseExprEnd[] = { 919, 702, 765, 840 };

          _SFD_CV_INIT_EML_SWITCH(0,1,0,678,687,951,4,&(caseStart[0]),
            &(caseExprEnd[0]));
        }

        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"fn_VectorToSkewSymmetricTensor",0,-1,433);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)
            c1_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c1_t_delta;
          real_T (*c1_omega)[3];
          real_T (*c1_Phi_r12)[9];
          real_T (*c1_M)[9];
          c1_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c1_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 1);
          c1_Phi_r12 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
          c1_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c1_omega);
          _SFD_SET_DATA_VALUE_PTR(1U, *c1_Phi_r12);
          _SFD_SET_DATA_VALUE_PTR(2U, *c1_M);
          _SFD_SET_DATA_VALUE_PTR(3U, c1_t_delta);
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
  return "WXBmQUuIQcrIgHsskbGIjC";
}

static void sf_opaque_initialize_c1_Model_01(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_Model_01InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c1_Model_01((SFc1_Model_01InstanceStruct*) chartInstanceVar);
  initialize_c1_Model_01((SFc1_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_Model_01(void *chartInstanceVar)
{
  enable_c1_Model_01((SFc1_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_Model_01(void *chartInstanceVar)
{
  disable_c1_Model_01((SFc1_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_Model_01(void *chartInstanceVar)
{
  sf_gateway_c1_Model_01((SFc1_Model_01InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c1_Model_01(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c1_Model_01((SFc1_Model_01InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c1_Model_01();/* state var info */
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

extern void sf_internal_set_sim_state_c1_Model_01(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c1_Model_01();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c1_Model_01((SFc1_Model_01InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c1_Model_01(SimStruct* S)
{
  return sf_internal_get_sim_state_c1_Model_01(S);
}

static void sf_opaque_set_sim_state_c1_Model_01(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c1_Model_01(S, st);
}

static void sf_opaque_terminate_c1_Model_01(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_Model_01InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_01_optimization_info();
    }

    finalize_c1_Model_01((SFc1_Model_01InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_Model_01((SFc1_Model_01InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_Model_01(SimStruct *S)
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
    initialize_params_c1_Model_01((SFc1_Model_01InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c1_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_01_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,1,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,1,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,1);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(309801711U));
  ssSetChecksum1(S,(1809273752U));
  ssSetChecksum2(S,(2982949255U));
  ssSetChecksum3(S,(2097790181U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c1_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_Model_01(SimStruct *S)
{
  SFc1_Model_01InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc1_Model_01InstanceStruct *)utMalloc(sizeof
    (SFc1_Model_01InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc1_Model_01InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c1_Model_01;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c1_Model_01;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c1_Model_01;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c1_Model_01;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c1_Model_01;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c1_Model_01;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c1_Model_01;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c1_Model_01;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_Model_01;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_Model_01;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c1_Model_01;
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

void c1_Model_01_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_Model_01(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_Model_01(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_Model_01(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_Model_01_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

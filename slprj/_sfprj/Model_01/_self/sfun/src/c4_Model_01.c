/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_01_sfun.h"
#include "c4_Model_01.h"
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
static const char * c4_debug_family_names[11] = { "v_lambda", "A", "Gamma", "i",
  "j", "nargin", "nargout", "omega", "t_delta", "M", "Phi_r22" };

/* Function Declarations */
static void initialize_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance);
static void initialize_params_c4_Model_01(SFc4_Model_01InstanceStruct
  *chartInstance);
static void enable_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance);
static void disable_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_update_debugger_state_c4_Model_01(SFc4_Model_01InstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c4_Model_01(SFc4_Model_01InstanceStruct
  *chartInstance);
static void set_sim_state_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_st);
static void finalize_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance);
static void sf_gateway_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_chartstep_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance);
static void initSimStructsc4_Model_01(SFc4_Model_01InstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber, uint32_T c4_instanceNumber);
static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData);
static void c4_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_Phi_r22, const char_T *c4_identifier, real_T c4_y[9]);
static void c4_b_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[9]);
static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static real_T c4_c_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static void c4_d_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[3]);
static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static void c4_info_helper(const mxArray **c4_info);
static const mxArray *c4_emlrt_marshallOut(const char * c4_u);
static const mxArray *c4_b_emlrt_marshallOut(const uint32_T c4_u);
static void c4_b_info_helper(const mxArray **c4_info);
static void c4_c_info_helper(const mxArray **c4_info);
static void c4_d_info_helper(const mxArray **c4_info);
static void c4_e_info_helper(const mxArray **c4_info);
static void c4_eig(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_A[9],
                   creal_T c4_V[3]);
static void c4_realmin(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_eml_error(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_eps(SFc4_Model_01InstanceStruct *chartInstance);
static real_T c4_eml_matlab_zlangeM(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_x[9]);
static real_T c4_abs(SFc4_Model_01InstanceStruct *chartInstance, creal_T c4_x);
static boolean_T c4_isfinite(SFc4_Model_01InstanceStruct *chartInstance, real_T
  c4_x);
static void c4_eml_matlab_zlascl(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_cfrom, real_T c4_cto, creal_T c4_A[9], creal_T c4_b_A[9]);
static real_T c4_b_abs(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_x);
static void c4_eml_matlab_zggbal(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], creal_T c4_b_A[9], int32_T *c4_ilo, int32_T *c4_ihi, int32_T
  c4_rscale[3]);
static void c4_eml_switch_helper(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_check_forloop_overflow_error(SFc4_Model_01InstanceStruct
  *chartInstance, boolean_T c4_overflow);
static void c4_eml_matlab_zlartg(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_f, creal_T c4_g, real_T *c4_cs, creal_T *c4_sn, creal_T *c4_r);
static void c4_eml_scalar_eg(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_eml_matlab_zhgeqz(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T c4_ilo, int32_T c4_ihi, real_T *c4_info, creal_T
  c4_alpha1[3], creal_T c4_beta1[3]);
static real_T c4_eml_matlab_zlanhs(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T c4_ilo, int32_T c4_ihi);
static int32_T c4_mod(SFc4_Model_01InstanceStruct *chartInstance, int32_T c4_x);
static creal_T c4_eml_div(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_x, real_T c4_y);
static void c4_scalarEg(SFc4_Model_01InstanceStruct *chartInstance);
static creal_T c4_sqrt(SFc4_Model_01InstanceStruct *chartInstance, creal_T c4_x);
static void c4_realmax(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_b_eml_matlab_zlartg(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_f, creal_T c4_g, real_T *c4_cs, creal_T *c4_sn);
static void c4_b_eml_matlab_zlascl(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_cfrom, real_T c4_cto, creal_T c4_A[3], creal_T c4_b_A[3]);
static void c4_b_eml_div(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_x[3], creal_T c4_y[3], creal_T c4_z[3]);
static void c4_eml_warning(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_b_eml_warning(SFc4_Model_01InstanceStruct *chartInstance);
static real_T c4_mpower(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_a);
static void c4_inv(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_x[9],
                   real_T c4_y[9]);
static void c4_inv3x3(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_x[9],
                      real_T c4_y[9]);
static real_T c4_norm(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_x[9]);
static void c4_c_eml_warning(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_d_eml_warning(SFc4_Model_01InstanceStruct *chartInstance, char_T
  c4_varargin_2[14]);
static void c4_b_eml_scalar_eg(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_c_eml_scalar_eg(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_eml_xgemm(SFc4_Model_01InstanceStruct *chartInstance, real_T
  c4_A[9], real_T c4_B[9], real_T c4_C[9], real_T c4_b_C[9]);
static void c4_matrix_to_scalar_power(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_a[9], real_T c4_b, real_T c4_c[9]);
static void c4_b_eml_error(SFc4_Model_01InstanceStruct *chartInstance);
static void c4_eml_matlab_zggev(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], real_T *c4_info, creal_T c4_alpha1[3], creal_T c4_beta1[3],
  creal_T c4_V[9]);
static void c4_eml_matlab_zgghrd(SFc4_Model_01InstanceStruct *chartInstance,
  int32_T c4_ilo, int32_T c4_ihi, creal_T c4_A[9], creal_T c4_b_A[9], creal_T
  c4_Z[9]);
static void c4_b_eml_matlab_zhgeqz(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T c4_ilo, int32_T c4_ihi, creal_T c4_Z[9], real_T
  *c4_info, creal_T c4_alpha1[3], creal_T c4_beta1[3], creal_T c4_b_A[9],
  creal_T c4_b_Z[9]);
static void c4_eml_matlab_ztgevc(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], creal_T c4_V[9], creal_T c4_b_V[9]);
static creal_T c4_rdivide(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_x, creal_T c4_y);
static real_T c4_eml_xnrm2(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_x[9], int32_T c4_ix0);
static creal_T c4_power(SFc4_Model_01InstanceStruct *chartInstance, creal_T c4_a,
  real_T c4_b);
static void c4_eml_lusolve(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_A[9], creal_T c4_B[9], creal_T c4_X[9]);
static void c4_e_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_sprintf, const char_T *c4_identifier, char_T c4_y[14]);
static void c4_f_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, char_T c4_y[14]);
static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static int32_T c4_g_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static uint8_T c4_h_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_b_is_active_c4_Model_01, const char_T *c4_identifier);
static uint8_T c4_i_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_c_eml_matlab_zlascl(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_cfrom, real_T c4_cto, creal_T c4_A[9]);
static void c4_b_eml_matlab_zggbal(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T *c4_ilo, int32_T *c4_ihi, int32_T c4_rscale[3]);
static void c4_b_sqrt(SFc4_Model_01InstanceStruct *chartInstance, creal_T *c4_x);
static void c4_d_eml_matlab_zlascl(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_cfrom, real_T c4_cto, creal_T c4_A[3]);
static void c4_b_eml_xgemm(SFc4_Model_01InstanceStruct *chartInstance, real_T
  c4_A[9], real_T c4_B[9], real_T c4_C[9]);
static void c4_b_eml_matlab_zgghrd(SFc4_Model_01InstanceStruct *chartInstance,
  int32_T c4_ilo, int32_T c4_ihi, creal_T c4_A[9], creal_T c4_Z[9]);
static void c4_c_eml_matlab_zhgeqz(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T c4_ilo, int32_T c4_ihi, creal_T c4_Z[9], real_T
  *c4_info, creal_T c4_alpha1[3], creal_T c4_beta1[3]);
static void c4_b_eml_matlab_ztgevc(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], creal_T c4_V[9]);
static int32_T c4_div_s32(SFc4_Model_01InstanceStruct *chartInstance, int32_T
  c4_numerator, int32_T c4_denominator);
static void init_dsm_address_info(SFc4_Model_01InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance)
{
  chartInstance->c4_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c4_is_active_c4_Model_01 = 0U;
}

static void initialize_params_c4_Model_01(SFc4_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c4_update_debugger_state_c4_Model_01(SFc4_Model_01InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c4_Model_01(SFc4_Model_01InstanceStruct
  *chartInstance)
{
  const mxArray *c4_st;
  const mxArray *c4_y = NULL;
  int32_T c4_i0;
  real_T c4_u[9];
  const mxArray *c4_b_y = NULL;
  uint8_T c4_hoistedGlobal;
  uint8_T c4_b_u;
  const mxArray *c4_c_y = NULL;
  real_T (*c4_Phi_r22)[9];
  c4_Phi_r22 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c4_st = NULL;
  c4_st = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_createcellmatrix(2, 1), false);
  for (c4_i0 = 0; c4_i0 < 9; c4_i0++) {
    c4_u[c4_i0] = (*c4_Phi_r22)[c4_i0];
  }

  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_setcell(c4_y, 0, c4_b_y);
  c4_hoistedGlobal = chartInstance->c4_is_active_c4_Model_01;
  c4_b_u = c4_hoistedGlobal;
  c4_c_y = NULL;
  sf_mex_assign(&c4_c_y, sf_mex_create("y", &c4_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c4_y, 1, c4_c_y);
  sf_mex_assign(&c4_st, c4_y, false);
  return c4_st;
}

static void set_sim_state_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_st)
{
  const mxArray *c4_u;
  real_T c4_dv0[9];
  int32_T c4_i1;
  real_T (*c4_Phi_r22)[9];
  c4_Phi_r22 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c4_doneDoubleBufferReInit = true;
  c4_u = sf_mex_dup(c4_st);
  c4_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c4_u, 0)),
                      "Phi_r22", c4_dv0);
  for (c4_i1 = 0; c4_i1 < 9; c4_i1++) {
    (*c4_Phi_r22)[c4_i1] = c4_dv0[c4_i1];
  }

  chartInstance->c4_is_active_c4_Model_01 = c4_h_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c4_u, 1)), "is_active_c4_Model_01");
  sf_mex_destroy(&c4_u);
  c4_update_debugger_state_c4_Model_01(chartInstance);
  sf_mex_destroy(&c4_st);
}

static void finalize_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance)
{
  int32_T c4_i2;
  int32_T c4_i3;
  int32_T c4_i4;
  real_T *c4_t_delta;
  real_T (*c4_M)[9];
  real_T (*c4_Phi_r22)[9];
  real_T (*c4_omega)[3];
  c4_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 2);
  c4_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c4_Phi_r22 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c4_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 3U, chartInstance->c4_sfEvent);
  for (c4_i2 = 0; c4_i2 < 3; c4_i2++) {
    _SFD_DATA_RANGE_CHECK((*c4_omega)[c4_i2], 0U);
  }

  chartInstance->c4_sfEvent = CALL_EVENT;
  c4_chartstep_c4_Model_01(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_01MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c4_i3 = 0; c4_i3 < 9; c4_i3++) {
    _SFD_DATA_RANGE_CHECK((*c4_Phi_r22)[c4_i3], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c4_t_delta, 2U);
  for (c4_i4 = 0; c4_i4 < 9; c4_i4++) {
    _SFD_DATA_RANGE_CHECK((*c4_M)[c4_i4], 3U);
  }
}

static void c4_chartstep_c4_Model_01(SFc4_Model_01InstanceStruct *chartInstance)
{
  real_T c4_hoistedGlobal;
  int32_T c4_i5;
  real_T c4_omega[3];
  real_T c4_t_delta;
  int32_T c4_i6;
  real_T c4_M[9];
  uint32_T c4_debug_family_var_map[11];
  real_T c4_v_lambda[3];
  real_T c4_A[9];
  real_T c4_Gamma[9];
  real_T c4_i;
  real_T c4_j;
  real_T c4_nargin = 3.0;
  real_T c4_nargout = 1.0;
  real_T c4_Phi_r22[9];
  int32_T c4_i7;
  int32_T c4_i8;
  real_T c4_b_M[9];
  creal_T c4_dcv0[3];
  int32_T c4_i9;
  real_T c4_d0;
  real_T c4_d1;
  real_T c4_d2;
  int32_T c4_i10;
  real_T c4_b_A[9];
  real_T c4_dv1[9];
  int32_T c4_i11;
  int32_T c4_b_i;
  int32_T c4_b_j;
  real_T c4_x;
  real_T c4_b_x;
  int32_T c4_i12;
  real_T c4_a[9];
  real_T c4_b;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_b_b;
  int32_T c4_i13;
  real_T c4_c[9];
  int32_T c4_k;
  real_T c4_b_k;
  real_T c4_e_x;
  real_T c4_e;
  boolean_T c4_firstmult;
  real_T c4_f_x;
  real_T c4_ed2;
  int32_T c4_i14;
  int32_T c4_i15;
  real_T c4_b_a[9];
  int32_T c4_i16;
  int32_T c4_i17;
  real_T c4_c_a[9];
  int32_T c4_i18;
  real_T c4_d_a[9];
  int32_T c4_i19;
  real_T c4_b_c[9];
  int32_T c4_i20;
  int32_T c4_i21;
  int32_T c4_i22;
  real_T c4_e_a[9];
  int32_T c4_i23;
  real_T c4_f_a[9];
  int32_T c4_i24;
  real_T c4_g_a[9];
  real_T c4_h_a;
  int32_T c4_i25;
  int32_T c4_i26;
  int32_T c4_i27;
  real_T *c4_b_t_delta;
  real_T (*c4_b_Phi_r22)[9];
  real_T (*c4_c_M)[9];
  real_T (*c4_b_omega)[3];
  int32_T exitg1;
  c4_c_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 2);
  c4_b_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c4_b_Phi_r22 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c4_b_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 3U, chartInstance->c4_sfEvent);
  c4_hoistedGlobal = *c4_b_t_delta;
  for (c4_i5 = 0; c4_i5 < 3; c4_i5++) {
    c4_omega[c4_i5] = (*c4_b_omega)[c4_i5];
  }

  c4_t_delta = c4_hoistedGlobal;
  for (c4_i6 = 0; c4_i6 < 9; c4_i6++) {
    c4_M[c4_i6] = (*c4_c_M)[c4_i6];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 11U, 11U, c4_debug_family_names,
    c4_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_v_lambda, 0U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_A, 1U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_Gamma, 2U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_i, 3U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_j, 4U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargin, 5U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargout, 6U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_omega, 7U, c4_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_t_delta, 8U, c4_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_M, 9U, c4_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_Phi_r22, 10U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 3);
  for (c4_i7 = 0; c4_i7 < 9; c4_i7++) {
    c4_Phi_r22[c4_i7] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 4);
  for (c4_i8 = 0; c4_i8 < 9; c4_i8++) {
    c4_b_M[c4_i8] = c4_M[c4_i8];
  }

  c4_eig(chartInstance, c4_b_M, c4_dcv0);
  for (c4_i9 = 0; c4_i9 < 3; c4_i9++) {
    c4_v_lambda[c4_i9] = c4_dcv0[c4_i9].re;
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 5);
  c4_d0 = c4_mpower(chartInstance, c4_v_lambda[0]);
  c4_d1 = c4_mpower(chartInstance, c4_v_lambda[1]);
  c4_d2 = c4_mpower(chartInstance, c4_v_lambda[2]);
  c4_A[0] = 1.0;
  c4_A[3] = c4_v_lambda[0];
  c4_A[6] = c4_d0;
  c4_A[1] = 1.0;
  c4_A[4] = c4_v_lambda[1];
  c4_A[7] = c4_d1;
  c4_A[2] = 1.0;
  c4_A[5] = c4_v_lambda[2];
  c4_A[8] = c4_d2;
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 6);
  for (c4_i10 = 0; c4_i10 < 9; c4_i10++) {
    c4_b_A[c4_i10] = c4_A[c4_i10];
  }

  c4_inv(chartInstance, c4_b_A, c4_dv1);
  for (c4_i11 = 0; c4_i11 < 9; c4_i11++) {
    c4_Gamma[c4_i11] = c4_dv1[c4_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 7);
  c4_i = 1.0;
  c4_b_i = 0;
  while (c4_b_i < 3) {
    c4_i = 1.0 + (real_T)c4_b_i;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 8);
    c4_j = 1.0;
    c4_b_j = 0;
    while (c4_b_j < 3) {
      c4_j = 1.0 + (real_T)c4_b_j;
      CV_EML_FOR(0, 1, 1, 1);
      _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 9);
      c4_x = c4_v_lambda[_SFD_EML_ARRAY_BOUNDS_CHECK("v_lambda", (int32_T)
        _SFD_INTEGER_CHECK("j", c4_j), 1, 3, 1, 0) - 1] * c4_t_delta;
      c4_b_x = c4_x;
      c4_b_x = muDoubleScalarExp(c4_b_x);
      for (c4_i12 = 0; c4_i12 < 9; c4_i12++) {
        c4_a[c4_i12] = c4_M[c4_i12];
      }

      c4_b = c4_i - 1.0;
      c4_c_x = c4_b;
      c4_d_x = c4_c_x;
      c4_d_x = muDoubleScalarFloor(c4_d_x);
      if (c4_d_x == c4_b) {
        c4_b_b = c4_b;
        c4_b_eml_scalar_eg(chartInstance);
        if (c4_b_b == 0.0) {
          for (c4_i13 = 0; c4_i13 < 9; c4_i13++) {
            c4_c[c4_i13] = 0.0;
          }

          for (c4_k = 0; c4_k < 3; c4_k++) {
            c4_b_k = 1.0 + (real_T)c4_k;
            c4_c[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c4_b_k), 1, 3, 2, 0) - 1)) -
              1] = 1.0;
          }
        } else {
          c4_e_x = c4_b_b;
          c4_e = muDoubleScalarAbs(c4_e_x);
          c4_firstmult = true;
          do {
            exitg1 = 0;
            c4_f_x = c4_e / 2.0;
            c4_ed2 = c4_f_x;
            c4_ed2 = muDoubleScalarFloor(c4_ed2);
            if (2.0 * c4_ed2 != c4_e) {
              if (c4_firstmult) {
                for (c4_i14 = 0; c4_i14 < 9; c4_i14++) {
                  c4_c[c4_i14] = c4_a[c4_i14];
                }

                c4_firstmult = false;
              } else {
                for (c4_i15 = 0; c4_i15 < 9; c4_i15++) {
                  c4_b_a[c4_i15] = c4_c[c4_i15];
                }

                c4_c_eml_scalar_eg(chartInstance);
                c4_c_eml_scalar_eg(chartInstance);
                for (c4_i16 = 0; c4_i16 < 9; c4_i16++) {
                  c4_c[c4_i16] = 0.0;
                }

                for (c4_i17 = 0; c4_i17 < 9; c4_i17++) {
                  c4_c_a[c4_i17] = c4_b_a[c4_i17];
                }

                for (c4_i18 = 0; c4_i18 < 9; c4_i18++) {
                  c4_d_a[c4_i18] = c4_a[c4_i18];
                }

                c4_b_eml_xgemm(chartInstance, c4_c_a, c4_d_a, c4_c);
              }
            }

            if (c4_ed2 == 0.0) {
              exitg1 = 1;
            } else {
              c4_e = c4_ed2;
              for (c4_i20 = 0; c4_i20 < 9; c4_i20++) {
                c4_b_a[c4_i20] = c4_a[c4_i20];
              }

              c4_c_eml_scalar_eg(chartInstance);
              c4_c_eml_scalar_eg(chartInstance);
              for (c4_i21 = 0; c4_i21 < 9; c4_i21++) {
                c4_a[c4_i21] = 0.0;
              }

              for (c4_i22 = 0; c4_i22 < 9; c4_i22++) {
                c4_e_a[c4_i22] = c4_b_a[c4_i22];
              }

              for (c4_i23 = 0; c4_i23 < 9; c4_i23++) {
                c4_f_a[c4_i23] = c4_b_a[c4_i23];
              }

              c4_b_eml_xgemm(chartInstance, c4_e_a, c4_f_a, c4_a);
            }
          } while (exitg1 == 0);

          if (c4_b_b < 0.0) {
            for (c4_i19 = 0; c4_i19 < 9; c4_i19++) {
              c4_b_c[c4_i19] = c4_c[c4_i19];
            }

            c4_inv(chartInstance, c4_b_c, c4_c);
          }
        }
      } else {
        for (c4_i24 = 0; c4_i24 < 9; c4_i24++) {
          c4_g_a[c4_i24] = c4_a[c4_i24];
        }

        c4_matrix_to_scalar_power(chartInstance, c4_g_a, c4_b, c4_c);
      }

      c4_h_a = c4_Gamma[(_SFD_EML_ARRAY_BOUNDS_CHECK("Gamma", (int32_T)
        _SFD_INTEGER_CHECK("i", c4_i), 1, 3, 1, 0) + 3 *
                         (_SFD_EML_ARRAY_BOUNDS_CHECK("Gamma", (int32_T)
        _SFD_INTEGER_CHECK("j", c4_j), 1, 3, 2, 0) - 1)) - 1] * c4_b_x;
      for (c4_i25 = 0; c4_i25 < 9; c4_i25++) {
        c4_c[c4_i25] *= c4_h_a;
      }

      for (c4_i26 = 0; c4_i26 < 9; c4_i26++) {
        c4_Phi_r22[c4_i26] += c4_c[c4_i26];
      }

      c4_b_j++;
      _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
    }

    CV_EML_FOR(0, 1, 1, 0);
    c4_b_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, -9);
  _SFD_SYMBOL_SCOPE_POP();
  for (c4_i27 = 0; c4_i27 < 9; c4_i27++) {
    (*c4_b_Phi_r22)[c4_i27] = c4_Phi_r22[c4_i27];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 3U, chartInstance->c4_sfEvent);
}

static void initSimStructsc4_Model_01(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber, uint32_T c4_instanceNumber)
{
  (void)c4_machineNumber;
  (void)c4_chartNumber;
  (void)c4_instanceNumber;
}

static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i28;
  int32_T c4_i29;
  int32_T c4_i30;
  real_T c4_b_inData[9];
  int32_T c4_i31;
  int32_T c4_i32;
  int32_T c4_i33;
  real_T c4_u[9];
  const mxArray *c4_y = NULL;
  SFc4_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc4_Model_01InstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_i28 = 0;
  for (c4_i29 = 0; c4_i29 < 3; c4_i29++) {
    for (c4_i30 = 0; c4_i30 < 3; c4_i30++) {
      c4_b_inData[c4_i30 + c4_i28] = (*(real_T (*)[9])c4_inData)[c4_i30 + c4_i28];
    }

    c4_i28 += 3;
  }

  c4_i31 = 0;
  for (c4_i32 = 0; c4_i32 < 3; c4_i32++) {
    for (c4_i33 = 0; c4_i33 < 3; c4_i33++) {
      c4_u[c4_i33 + c4_i31] = c4_b_inData[c4_i33 + c4_i31];
    }

    c4_i31 += 3;
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static void c4_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_Phi_r22, const char_T *c4_identifier, real_T c4_y[9])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_Phi_r22), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_Phi_r22);
}

static void c4_b_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[9])
{
  real_T c4_dv2[9];
  int32_T c4_i34;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv2, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c4_i34 = 0; c4_i34 < 9; c4_i34++) {
    c4_y[c4_i34] = c4_dv2[c4_i34];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_Phi_r22;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[9];
  int32_T c4_i35;
  int32_T c4_i36;
  int32_T c4_i37;
  SFc4_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc4_Model_01InstanceStruct *)chartInstanceVoid;
  c4_Phi_r22 = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_Phi_r22), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_Phi_r22);
  c4_i35 = 0;
  for (c4_i36 = 0; c4_i36 < 3; c4_i36++) {
    for (c4_i37 = 0; c4_i37 < 3; c4_i37++) {
      (*(real_T (*)[9])c4_outData)[c4_i37 + c4_i35] = c4_y[c4_i37 + c4_i35];
    }

    c4_i35 += 3;
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  real_T c4_u;
  const mxArray *c4_y = NULL;
  SFc4_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc4_Model_01InstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_u = *(real_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i38;
  real_T c4_b_inData[3];
  int32_T c4_i39;
  real_T c4_u[3];
  const mxArray *c4_y = NULL;
  SFc4_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc4_Model_01InstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i38 = 0; c4_i38 < 3; c4_i38++) {
    c4_b_inData[c4_i38] = (*(real_T (*)[3])c4_inData)[c4_i38];
  }

  for (c4_i39 = 0; c4_i39 < 3; c4_i39++) {
    c4_u[c4_i39] = c4_b_inData[c4_i39];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static real_T c4_c_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  real_T c4_y;
  real_T c4_d3;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_d3, 1, 0, 0U, 0, 0U, 0);
  c4_y = c4_d3;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_nargout;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y;
  SFc4_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc4_Model_01InstanceStruct *)chartInstanceVoid;
  c4_nargout = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_nargout), &c4_thisId);
  sf_mex_destroy(&c4_nargout);
  *(real_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static void c4_d_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[3])
{
  real_T c4_dv3[3];
  int32_T c4_i40;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv3, 1, 0, 0U, 1, 0U, 1, 3);
  for (c4_i40 = 0; c4_i40 < 3; c4_i40++) {
    c4_y[c4_i40] = c4_dv3[c4_i40];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_v_lambda;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[3];
  int32_T c4_i41;
  SFc4_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc4_Model_01InstanceStruct *)chartInstanceVoid;
  c4_v_lambda = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_v_lambda), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_v_lambda);
  for (c4_i41 = 0; c4_i41 < 3; c4_i41++) {
    (*(real_T (*)[3])c4_outData)[c4_i41] = c4_y[c4_i41];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

const mxArray *sf_c4_Model_01_get_eml_resolved_functions_info(void)
{
  const mxArray *c4_nameCaptureInfo = NULL;
  c4_nameCaptureInfo = NULL;
  sf_mex_assign(&c4_nameCaptureInfo, sf_mex_createstruct("structure", 2, 270, 1),
                false);
  c4_info_helper(&c4_nameCaptureInfo);
  c4_b_info_helper(&c4_nameCaptureInfo);
  c4_c_info_helper(&c4_nameCaptureInfo);
  c4_d_info_helper(&c4_nameCaptureInfo);
  c4_e_info_helper(&c4_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c4_nameCaptureInfo);
  return c4_nameCaptureInfo;
}

static void c4_info_helper(const mxArray **c4_info)
{
  const mxArray *c4_rhs0 = NULL;
  const mxArray *c4_lhs0 = NULL;
  const mxArray *c4_rhs1 = NULL;
  const mxArray *c4_lhs1 = NULL;
  const mxArray *c4_rhs2 = NULL;
  const mxArray *c4_lhs2 = NULL;
  const mxArray *c4_rhs3 = NULL;
  const mxArray *c4_lhs3 = NULL;
  const mxArray *c4_rhs4 = NULL;
  const mxArray *c4_lhs4 = NULL;
  const mxArray *c4_rhs5 = NULL;
  const mxArray *c4_lhs5 = NULL;
  const mxArray *c4_rhs6 = NULL;
  const mxArray *c4_lhs6 = NULL;
  const mxArray *c4_rhs7 = NULL;
  const mxArray *c4_lhs7 = NULL;
  const mxArray *c4_rhs8 = NULL;
  const mxArray *c4_lhs8 = NULL;
  const mxArray *c4_rhs9 = NULL;
  const mxArray *c4_lhs9 = NULL;
  const mxArray *c4_rhs10 = NULL;
  const mxArray *c4_lhs10 = NULL;
  const mxArray *c4_rhs11 = NULL;
  const mxArray *c4_lhs11 = NULL;
  const mxArray *c4_rhs12 = NULL;
  const mxArray *c4_lhs12 = NULL;
  const mxArray *c4_rhs13 = NULL;
  const mxArray *c4_lhs13 = NULL;
  const mxArray *c4_rhs14 = NULL;
  const mxArray *c4_lhs14 = NULL;
  const mxArray *c4_rhs15 = NULL;
  const mxArray *c4_lhs15 = NULL;
  const mxArray *c4_rhs16 = NULL;
  const mxArray *c4_lhs16 = NULL;
  const mxArray *c4_rhs17 = NULL;
  const mxArray *c4_lhs17 = NULL;
  const mxArray *c4_rhs18 = NULL;
  const mxArray *c4_lhs18 = NULL;
  const mxArray *c4_rhs19 = NULL;
  const mxArray *c4_lhs19 = NULL;
  const mxArray *c4_rhs20 = NULL;
  const mxArray *c4_lhs20 = NULL;
  const mxArray *c4_rhs21 = NULL;
  const mxArray *c4_lhs21 = NULL;
  const mxArray *c4_rhs22 = NULL;
  const mxArray *c4_lhs22 = NULL;
  const mxArray *c4_rhs23 = NULL;
  const mxArray *c4_lhs23 = NULL;
  const mxArray *c4_rhs24 = NULL;
  const mxArray *c4_lhs24 = NULL;
  const mxArray *c4_rhs25 = NULL;
  const mxArray *c4_lhs25 = NULL;
  const mxArray *c4_rhs26 = NULL;
  const mxArray *c4_lhs26 = NULL;
  const mxArray *c4_rhs27 = NULL;
  const mxArray *c4_lhs27 = NULL;
  const mxArray *c4_rhs28 = NULL;
  const mxArray *c4_lhs28 = NULL;
  const mxArray *c4_rhs29 = NULL;
  const mxArray *c4_lhs29 = NULL;
  const mxArray *c4_rhs30 = NULL;
  const mxArray *c4_lhs30 = NULL;
  const mxArray *c4_rhs31 = NULL;
  const mxArray *c4_lhs31 = NULL;
  const mxArray *c4_rhs32 = NULL;
  const mxArray *c4_lhs32 = NULL;
  const mxArray *c4_rhs33 = NULL;
  const mxArray *c4_lhs33 = NULL;
  const mxArray *c4_rhs34 = NULL;
  const mxArray *c4_lhs34 = NULL;
  const mxArray *c4_rhs35 = NULL;
  const mxArray *c4_lhs35 = NULL;
  const mxArray *c4_rhs36 = NULL;
  const mxArray *c4_lhs36 = NULL;
  const mxArray *c4_rhs37 = NULL;
  const mxArray *c4_lhs37 = NULL;
  const mxArray *c4_rhs38 = NULL;
  const mxArray *c4_lhs38 = NULL;
  const mxArray *c4_rhs39 = NULL;
  const mxArray *c4_lhs39 = NULL;
  const mxArray *c4_rhs40 = NULL;
  const mxArray *c4_lhs40 = NULL;
  const mxArray *c4_rhs41 = NULL;
  const mxArray *c4_lhs41 = NULL;
  const mxArray *c4_rhs42 = NULL;
  const mxArray *c4_lhs42 = NULL;
  const mxArray *c4_rhs43 = NULL;
  const mxArray *c4_lhs43 = NULL;
  const mxArray *c4_rhs44 = NULL;
  const mxArray *c4_lhs44 = NULL;
  const mxArray *c4_rhs45 = NULL;
  const mxArray *c4_lhs45 = NULL;
  const mxArray *c4_rhs46 = NULL;
  const mxArray *c4_lhs46 = NULL;
  const mxArray *c4_rhs47 = NULL;
  const mxArray *c4_lhs47 = NULL;
  const mxArray *c4_rhs48 = NULL;
  const mxArray *c4_lhs48 = NULL;
  const mxArray *c4_rhs49 = NULL;
  const mxArray *c4_lhs49 = NULL;
  const mxArray *c4_rhs50 = NULL;
  const mxArray *c4_lhs50 = NULL;
  const mxArray *c4_rhs51 = NULL;
  const mxArray *c4_lhs51 = NULL;
  const mxArray *c4_rhs52 = NULL;
  const mxArray *c4_lhs52 = NULL;
  const mxArray *c4_rhs53 = NULL;
  const mxArray *c4_lhs53 = NULL;
  const mxArray *c4_rhs54 = NULL;
  const mxArray *c4_lhs54 = NULL;
  const mxArray *c4_rhs55 = NULL;
  const mxArray *c4_lhs55 = NULL;
  const mxArray *c4_rhs56 = NULL;
  const mxArray *c4_lhs56 = NULL;
  const mxArray *c4_rhs57 = NULL;
  const mxArray *c4_lhs57 = NULL;
  const mxArray *c4_rhs58 = NULL;
  const mxArray *c4_lhs58 = NULL;
  const mxArray *c4_rhs59 = NULL;
  const mxArray *c4_lhs59 = NULL;
  const mxArray *c4_rhs60 = NULL;
  const mxArray *c4_lhs60 = NULL;
  const mxArray *c4_rhs61 = NULL;
  const mxArray *c4_lhs61 = NULL;
  const mxArray *c4_rhs62 = NULL;
  const mxArray *c4_lhs62 = NULL;
  const mxArray *c4_rhs63 = NULL;
  const mxArray *c4_lhs63 = NULL;
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eig"), "name", "name", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c4_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c4_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c4_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_xgeev"), "name", "name", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826004U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c4_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeev.m"),
                  "context", "context", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_lapack_xgeev"), "name",
                  "name", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1301335668U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c4_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zggev"), "name",
                  "name", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826018U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c4_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("realmin"), "name", "name", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c4_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_realmin"), "name", "name",
                  7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658444U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c4_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c4_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("sqrt"), "name", "name", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c4_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_error"), "name", "name",
                  10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c4_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825938U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c4_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eps"), "name", "name", 12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c4_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c4_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_eps"), "name", "name", 14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c4_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c4_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zlangeM"), "name",
                  "name", 16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826020U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c4_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c4_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c4_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c4_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "context", "context", 20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1350417854U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c4_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "context", "context", 21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isnan"), "name", "name", 21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c4_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c4_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlangeM.m"),
                  "context", "context", 23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c4_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "context", "context", 24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c4_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isfinite"), "name", "name", 25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c4_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c4_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isinf"), "name", "name", 27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c4_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c4_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isnan"), "name", "name", 29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c4_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c4_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zlascl"), "name",
                  "name", 31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826022U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c4_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("realmin"), "name", "name", 32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c4_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 33);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eps"), "name", "name", 33);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c4_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 34);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 34);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 34);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c4_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlascl.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 35);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c4_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 36);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 36);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c4_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 37);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zggbal"), "name",
                  "name", 37);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826018U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c4_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "context", "context", 38);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 38);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c4_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_rows"),
                  "context", "context", 39);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 39);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c4_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_rows"),
                  "context", "context", 40);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 40);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c4_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 41);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("intmax"), "name", "name", 41);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 41);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c4_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 42);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 42);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c4_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_rows"),
                  "context", "context", 43);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 43);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 43);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c4_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 44);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 44);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 44);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c4_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_simtran"),
                  "context", "context", 45);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 45);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c4_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_simtran"),
                  "context", "context", 46);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 46);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c4_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "context", "context", 47);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 47);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 47);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c4_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_cols"),
                  "context", "context", 48);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 48);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c4_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m!eml_zggbal_eigsearch_cols"),
                  "context", "context", 49);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 49);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c4_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbal.m"),
                  "context", "context", 50);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 50);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 50);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c4_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 51);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 51);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 51);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c4_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 52);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zgghrd"), "name",
                  "name", 52);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826020U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c4_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 53);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 53);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c4_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 54);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 54);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c4_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 55);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 55);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c4_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 56);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 56);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c4_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 57);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 57);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 57);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c4_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 58);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zlartg"), "name",
                  "name", 58);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826022U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c4_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 59);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("realmin"), "name", "name", 59);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 59);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c4_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 60);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eps"), "name", "name", 60);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 60);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c4_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 61);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("fix"), "name", "name", 61);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m"), "resolved",
                  "resolved", 61);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c4_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m"), "context",
                  "context", 62);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 62);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c4_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/fix.m"), "context",
                  "context", 63);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_fix"), "name",
                  "name", 63);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_fix.m"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658438U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c4_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c4_rhs0);
  sf_mex_destroy(&c4_lhs0);
  sf_mex_destroy(&c4_rhs1);
  sf_mex_destroy(&c4_lhs1);
  sf_mex_destroy(&c4_rhs2);
  sf_mex_destroy(&c4_lhs2);
  sf_mex_destroy(&c4_rhs3);
  sf_mex_destroy(&c4_lhs3);
  sf_mex_destroy(&c4_rhs4);
  sf_mex_destroy(&c4_lhs4);
  sf_mex_destroy(&c4_rhs5);
  sf_mex_destroy(&c4_lhs5);
  sf_mex_destroy(&c4_rhs6);
  sf_mex_destroy(&c4_lhs6);
  sf_mex_destroy(&c4_rhs7);
  sf_mex_destroy(&c4_lhs7);
  sf_mex_destroy(&c4_rhs8);
  sf_mex_destroy(&c4_lhs8);
  sf_mex_destroy(&c4_rhs9);
  sf_mex_destroy(&c4_lhs9);
  sf_mex_destroy(&c4_rhs10);
  sf_mex_destroy(&c4_lhs10);
  sf_mex_destroy(&c4_rhs11);
  sf_mex_destroy(&c4_lhs11);
  sf_mex_destroy(&c4_rhs12);
  sf_mex_destroy(&c4_lhs12);
  sf_mex_destroy(&c4_rhs13);
  sf_mex_destroy(&c4_lhs13);
  sf_mex_destroy(&c4_rhs14);
  sf_mex_destroy(&c4_lhs14);
  sf_mex_destroy(&c4_rhs15);
  sf_mex_destroy(&c4_lhs15);
  sf_mex_destroy(&c4_rhs16);
  sf_mex_destroy(&c4_lhs16);
  sf_mex_destroy(&c4_rhs17);
  sf_mex_destroy(&c4_lhs17);
  sf_mex_destroy(&c4_rhs18);
  sf_mex_destroy(&c4_lhs18);
  sf_mex_destroy(&c4_rhs19);
  sf_mex_destroy(&c4_lhs19);
  sf_mex_destroy(&c4_rhs20);
  sf_mex_destroy(&c4_lhs20);
  sf_mex_destroy(&c4_rhs21);
  sf_mex_destroy(&c4_lhs21);
  sf_mex_destroy(&c4_rhs22);
  sf_mex_destroy(&c4_lhs22);
  sf_mex_destroy(&c4_rhs23);
  sf_mex_destroy(&c4_lhs23);
  sf_mex_destroy(&c4_rhs24);
  sf_mex_destroy(&c4_lhs24);
  sf_mex_destroy(&c4_rhs25);
  sf_mex_destroy(&c4_lhs25);
  sf_mex_destroy(&c4_rhs26);
  sf_mex_destroy(&c4_lhs26);
  sf_mex_destroy(&c4_rhs27);
  sf_mex_destroy(&c4_lhs27);
  sf_mex_destroy(&c4_rhs28);
  sf_mex_destroy(&c4_lhs28);
  sf_mex_destroy(&c4_rhs29);
  sf_mex_destroy(&c4_lhs29);
  sf_mex_destroy(&c4_rhs30);
  sf_mex_destroy(&c4_lhs30);
  sf_mex_destroy(&c4_rhs31);
  sf_mex_destroy(&c4_lhs31);
  sf_mex_destroy(&c4_rhs32);
  sf_mex_destroy(&c4_lhs32);
  sf_mex_destroy(&c4_rhs33);
  sf_mex_destroy(&c4_lhs33);
  sf_mex_destroy(&c4_rhs34);
  sf_mex_destroy(&c4_lhs34);
  sf_mex_destroy(&c4_rhs35);
  sf_mex_destroy(&c4_lhs35);
  sf_mex_destroy(&c4_rhs36);
  sf_mex_destroy(&c4_lhs36);
  sf_mex_destroy(&c4_rhs37);
  sf_mex_destroy(&c4_lhs37);
  sf_mex_destroy(&c4_rhs38);
  sf_mex_destroy(&c4_lhs38);
  sf_mex_destroy(&c4_rhs39);
  sf_mex_destroy(&c4_lhs39);
  sf_mex_destroy(&c4_rhs40);
  sf_mex_destroy(&c4_lhs40);
  sf_mex_destroy(&c4_rhs41);
  sf_mex_destroy(&c4_lhs41);
  sf_mex_destroy(&c4_rhs42);
  sf_mex_destroy(&c4_lhs42);
  sf_mex_destroy(&c4_rhs43);
  sf_mex_destroy(&c4_lhs43);
  sf_mex_destroy(&c4_rhs44);
  sf_mex_destroy(&c4_lhs44);
  sf_mex_destroy(&c4_rhs45);
  sf_mex_destroy(&c4_lhs45);
  sf_mex_destroy(&c4_rhs46);
  sf_mex_destroy(&c4_lhs46);
  sf_mex_destroy(&c4_rhs47);
  sf_mex_destroy(&c4_lhs47);
  sf_mex_destroy(&c4_rhs48);
  sf_mex_destroy(&c4_lhs48);
  sf_mex_destroy(&c4_rhs49);
  sf_mex_destroy(&c4_lhs49);
  sf_mex_destroy(&c4_rhs50);
  sf_mex_destroy(&c4_lhs50);
  sf_mex_destroy(&c4_rhs51);
  sf_mex_destroy(&c4_lhs51);
  sf_mex_destroy(&c4_rhs52);
  sf_mex_destroy(&c4_lhs52);
  sf_mex_destroy(&c4_rhs53);
  sf_mex_destroy(&c4_lhs53);
  sf_mex_destroy(&c4_rhs54);
  sf_mex_destroy(&c4_lhs54);
  sf_mex_destroy(&c4_rhs55);
  sf_mex_destroy(&c4_lhs55);
  sf_mex_destroy(&c4_rhs56);
  sf_mex_destroy(&c4_lhs56);
  sf_mex_destroy(&c4_rhs57);
  sf_mex_destroy(&c4_lhs57);
  sf_mex_destroy(&c4_rhs58);
  sf_mex_destroy(&c4_lhs58);
  sf_mex_destroy(&c4_rhs59);
  sf_mex_destroy(&c4_lhs59);
  sf_mex_destroy(&c4_rhs60);
  sf_mex_destroy(&c4_lhs60);
  sf_mex_destroy(&c4_rhs61);
  sf_mex_destroy(&c4_lhs61);
  sf_mex_destroy(&c4_rhs62);
  sf_mex_destroy(&c4_lhs62);
  sf_mex_destroy(&c4_rhs63);
  sf_mex_destroy(&c4_lhs63);
}

static const mxArray *c4_emlrt_marshallOut(const char * c4_u)
{
  const mxArray *c4_y = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c4_u)), false);
  return c4_y;
}

static const mxArray *c4_b_emlrt_marshallOut(const uint32_T c4_u)
{
  const mxArray *c4_y = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 7, 0U, 0U, 0U, 0), false);
  return c4_y;
}

static void c4_b_info_helper(const mxArray **c4_info)
{
  const mxArray *c4_rhs64 = NULL;
  const mxArray *c4_lhs64 = NULL;
  const mxArray *c4_rhs65 = NULL;
  const mxArray *c4_lhs65 = NULL;
  const mxArray *c4_rhs66 = NULL;
  const mxArray *c4_lhs66 = NULL;
  const mxArray *c4_rhs67 = NULL;
  const mxArray *c4_lhs67 = NULL;
  const mxArray *c4_rhs68 = NULL;
  const mxArray *c4_lhs68 = NULL;
  const mxArray *c4_rhs69 = NULL;
  const mxArray *c4_lhs69 = NULL;
  const mxArray *c4_rhs70 = NULL;
  const mxArray *c4_lhs70 = NULL;
  const mxArray *c4_rhs71 = NULL;
  const mxArray *c4_lhs71 = NULL;
  const mxArray *c4_rhs72 = NULL;
  const mxArray *c4_lhs72 = NULL;
  const mxArray *c4_rhs73 = NULL;
  const mxArray *c4_lhs73 = NULL;
  const mxArray *c4_rhs74 = NULL;
  const mxArray *c4_lhs74 = NULL;
  const mxArray *c4_rhs75 = NULL;
  const mxArray *c4_lhs75 = NULL;
  const mxArray *c4_rhs76 = NULL;
  const mxArray *c4_lhs76 = NULL;
  const mxArray *c4_rhs77 = NULL;
  const mxArray *c4_lhs77 = NULL;
  const mxArray *c4_rhs78 = NULL;
  const mxArray *c4_lhs78 = NULL;
  const mxArray *c4_rhs79 = NULL;
  const mxArray *c4_lhs79 = NULL;
  const mxArray *c4_rhs80 = NULL;
  const mxArray *c4_lhs80 = NULL;
  const mxArray *c4_rhs81 = NULL;
  const mxArray *c4_lhs81 = NULL;
  const mxArray *c4_rhs82 = NULL;
  const mxArray *c4_lhs82 = NULL;
  const mxArray *c4_rhs83 = NULL;
  const mxArray *c4_lhs83 = NULL;
  const mxArray *c4_rhs84 = NULL;
  const mxArray *c4_lhs84 = NULL;
  const mxArray *c4_rhs85 = NULL;
  const mxArray *c4_lhs85 = NULL;
  const mxArray *c4_rhs86 = NULL;
  const mxArray *c4_lhs86 = NULL;
  const mxArray *c4_rhs87 = NULL;
  const mxArray *c4_lhs87 = NULL;
  const mxArray *c4_rhs88 = NULL;
  const mxArray *c4_lhs88 = NULL;
  const mxArray *c4_rhs89 = NULL;
  const mxArray *c4_lhs89 = NULL;
  const mxArray *c4_rhs90 = NULL;
  const mxArray *c4_lhs90 = NULL;
  const mxArray *c4_rhs91 = NULL;
  const mxArray *c4_lhs91 = NULL;
  const mxArray *c4_rhs92 = NULL;
  const mxArray *c4_lhs92 = NULL;
  const mxArray *c4_rhs93 = NULL;
  const mxArray *c4_lhs93 = NULL;
  const mxArray *c4_rhs94 = NULL;
  const mxArray *c4_lhs94 = NULL;
  const mxArray *c4_rhs95 = NULL;
  const mxArray *c4_lhs95 = NULL;
  const mxArray *c4_rhs96 = NULL;
  const mxArray *c4_lhs96 = NULL;
  const mxArray *c4_rhs97 = NULL;
  const mxArray *c4_lhs97 = NULL;
  const mxArray *c4_rhs98 = NULL;
  const mxArray *c4_lhs98 = NULL;
  const mxArray *c4_rhs99 = NULL;
  const mxArray *c4_lhs99 = NULL;
  const mxArray *c4_rhs100 = NULL;
  const mxArray *c4_lhs100 = NULL;
  const mxArray *c4_rhs101 = NULL;
  const mxArray *c4_lhs101 = NULL;
  const mxArray *c4_rhs102 = NULL;
  const mxArray *c4_lhs102 = NULL;
  const mxArray *c4_rhs103 = NULL;
  const mxArray *c4_lhs103 = NULL;
  const mxArray *c4_rhs104 = NULL;
  const mxArray *c4_lhs104 = NULL;
  const mxArray *c4_rhs105 = NULL;
  const mxArray *c4_lhs105 = NULL;
  const mxArray *c4_rhs106 = NULL;
  const mxArray *c4_lhs106 = NULL;
  const mxArray *c4_rhs107 = NULL;
  const mxArray *c4_lhs107 = NULL;
  const mxArray *c4_rhs108 = NULL;
  const mxArray *c4_lhs108 = NULL;
  const mxArray *c4_rhs109 = NULL;
  const mxArray *c4_lhs109 = NULL;
  const mxArray *c4_rhs110 = NULL;
  const mxArray *c4_lhs110 = NULL;
  const mxArray *c4_rhs111 = NULL;
  const mxArray *c4_lhs111 = NULL;
  const mxArray *c4_rhs112 = NULL;
  const mxArray *c4_lhs112 = NULL;
  const mxArray *c4_rhs113 = NULL;
  const mxArray *c4_lhs113 = NULL;
  const mxArray *c4_rhs114 = NULL;
  const mxArray *c4_lhs114 = NULL;
  const mxArray *c4_rhs115 = NULL;
  const mxArray *c4_lhs115 = NULL;
  const mxArray *c4_rhs116 = NULL;
  const mxArray *c4_lhs116 = NULL;
  const mxArray *c4_rhs117 = NULL;
  const mxArray *c4_lhs117 = NULL;
  const mxArray *c4_rhs118 = NULL;
  const mxArray *c4_lhs118 = NULL;
  const mxArray *c4_rhs119 = NULL;
  const mxArray *c4_lhs119 = NULL;
  const mxArray *c4_rhs120 = NULL;
  const mxArray *c4_lhs120 = NULL;
  const mxArray *c4_rhs121 = NULL;
  const mxArray *c4_lhs121 = NULL;
  const mxArray *c4_rhs122 = NULL;
  const mxArray *c4_lhs122 = NULL;
  const mxArray *c4_rhs123 = NULL;
  const mxArray *c4_lhs123 = NULL;
  const mxArray *c4_rhs124 = NULL;
  const mxArray *c4_lhs124 = NULL;
  const mxArray *c4_rhs125 = NULL;
  const mxArray *c4_lhs125 = NULL;
  const mxArray *c4_rhs126 = NULL;
  const mxArray *c4_lhs126 = NULL;
  const mxArray *c4_rhs127 = NULL;
  const mxArray *c4_lhs127 = NULL;
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 64);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mpower"), "name", "name", 64);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 64);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c4_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 65);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 65);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c4_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 66);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("ismatrix"), "name", "name", 66);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 66);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c4_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 67);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("power"), "name", "name", 67);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 67);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c4_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 68);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 68);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c4_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 69);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 69);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 69);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c4_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 70);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 70);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c4_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 71);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 71);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c4_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 72);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("floor"), "name", "name", 72);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 72);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c4_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 73);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 73);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c4_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 74);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 74);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c4_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 75);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 75);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 75);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c4_rhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 76);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 76);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 76);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c4_rhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m!absinf"),
                  "context", "context", 77);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 77);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 77);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 77);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c4_rhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 78);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 78);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c4_rhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 79);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 79);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c4_rhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 80);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  80);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 80);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1350417854U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c4_rhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 81);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("sqrt"), "name", "name", 81);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 81);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c4_rhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "context", "context", 82);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 82);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 82);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c4_rhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 83);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_zrot_rows"), "name",
                  "name", 83);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "resolved", "resolved", 83);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c4_rhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "context", "context", 84);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 84);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 84);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c4_rhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "context", "context", 85);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 85);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c4_rhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "context", "context", 86);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.conjtimes"),
                  "name", "name", 86);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/conjtimes.m"),
                  "resolved", "resolved", 86);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1360286186U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c4_rhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 87);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_zrot_cols"), "name",
                  "name", 87);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c4_rhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "context", "context", 88);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 88);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 88);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c4_rhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "context", "context", 89);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 89);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c4_rhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "context", "context", 90);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.conjtimes"),
                  "name", "name", 90);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/conjtimes.m"),
                  "resolved", "resolved", 90);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1360286186U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c4_rhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 91);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zhgeqz"), "name",
                  "name", 91);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1368190232U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c4_rhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 92);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 92);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 92);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c4_rhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 93);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eps"), "name", "name", 93);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 93);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c4_rhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 94);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("realmin"), "name", "name", 94);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 94);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 94);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c4_rhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 95);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zlanhs"), "name",
                  "name", 95);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 95);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "resolved", "resolved", 95);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826020U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c4_rhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 96);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 96);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 96);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 96);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c4_rhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 97);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 97);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 97);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 97);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c4_rhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 98);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 98);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 98);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 98);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c4_rhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 99);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 99);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 99);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c4_rhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlanhs.m"),
                  "context", "context", 100);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("sqrt"), "name", "name", 100);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 100);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c4_rhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs100), "rhs", "rhs",
                  100);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs100), "lhs", "lhs",
                  100);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 101);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 101);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 101);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 101);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c4_rhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs101), "rhs", "rhs",
                  101);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs101), "lhs", "lhs",
                  101);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 102);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 102);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 102);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c4_rhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs102), "rhs", "rhs",
                  102);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs102), "lhs", "lhs",
                  102);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 103);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 103);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 103);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 103);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c4_rhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs103), "rhs", "rhs",
                  103);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs103), "lhs", "lhs",
                  103);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 104);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 104);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 104);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c4_rhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs104), "rhs", "rhs",
                  104);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs104), "lhs", "lhs",
                  104);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 105);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 105);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 105);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 105);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c4_rhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs105), "rhs", "rhs",
                  105);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs105), "lhs", "lhs",
                  105);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m!abs1"),
                  "context", "context", 106);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 106);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 106);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 106);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c4_rhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs106), "rhs", "rhs",
                  106);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs106), "lhs", "lhs",
                  106);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 107);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 107);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 107);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 107);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c4_rhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs107), "rhs", "rhs",
                  107);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs107), "lhs", "lhs",
                  107);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 108);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zlartg"), "name",
                  "name", 108);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 108);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlartg.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826022U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c4_rhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs108), "rhs", "rhs",
                  108);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs108), "lhs", "lhs",
                  108);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 109);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_zrot_cols"), "name",
                  "name", 109);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 109);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_cols.m"),
                  "resolved", "resolved", 109);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c4_rhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs109), "rhs", "rhs",
                  109);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs109), "lhs", "lhs",
                  109);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 110);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mod"), "name", "name", 110);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "resolved",
                  "resolved", 110);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c4_rhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs110), "rhs", "rhs",
                  110);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs110), "lhs", "lhs",
                  110);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 111);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 111);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 111);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 111);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c4_rhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs111), "rhs", "rhs",
                  111);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs111), "lhs", "lhs",
                  111);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 112);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 112);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 112);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c4_rhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs112), "rhs", "rhs",
                  112);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs112), "lhs", "lhs",
                  112);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 113);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 113);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c4_rhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs113), "rhs", "rhs",
                  113);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs113), "lhs", "lhs",
                  113);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 114);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 114);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 114);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c4_rhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs114), "rhs", "rhs",
                  114);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs114), "lhs", "lhs",
                  114);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "context",
                  "context", 115);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 115);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 115);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c4_rhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs115), "rhs", "rhs",
                  115);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs115), "lhs", "lhs",
                  115);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!intmod"), "context",
                  "context", 116);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 116);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 116);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 116);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c4_rhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs116), "rhs", "rhs",
                  116);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs116), "lhs", "lhs",
                  116);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 117);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 117);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 117);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 117);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c4_rhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs117), "rhs", "rhs",
                  117);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs117), "lhs", "lhs",
                  117);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 118);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_div"), "name", "name", 118);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 118);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 118);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c4_rhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs118), "rhs", "rhs",
                  118);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs118), "lhs", "lhs",
                  118);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 119);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 119);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 119);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 119);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c4_rhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs119), "rhs", "rhs",
                  119);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs119), "lhs", "lhs",
                  119);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p!eml_fldiv"),
                  "context", "context", 120);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 120);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 120);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 120);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c4_rhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs120), "rhs", "rhs",
                  120);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs120), "lhs", "lhs",
                  120);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p!eml_fldiv"),
                  "context", "context", 121);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 121);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 121);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 121);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c4_rhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs121), "rhs", "rhs",
                  121);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs121), "lhs", "lhs",
                  121);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p!eml_fldiv"),
                  "context", "context", 122);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 122);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 122);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c4_rhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs122), "rhs", "rhs",
                  122);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs122), "lhs", "lhs",
                  122);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 123);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("sqrt"), "name", "name", 123);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 123);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 123);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c4_rhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs123), "rhs", "rhs",
                  123);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs123), "lhs", "lhs",
                  123);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 124);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("rdivide"), "name", "name", 124);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 124);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 124);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c4_rhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs124), "rhs", "rhs",
                  124);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs124), "lhs", "lhs",
                  124);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 125);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 125);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 125);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 125);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c4_rhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs125), "rhs", "rhs",
                  125);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs125), "lhs", "lhs",
                  125);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 126);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 126);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 126);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c4_rhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs126), "rhs", "rhs",
                  126);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs126), "lhs", "lhs",
                  126);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 127);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_div"), "name", "name", 127);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 127);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c4_rhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs127), "rhs", "rhs",
                  127);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs127), "lhs", "lhs",
                  127);
  sf_mex_destroy(&c4_rhs64);
  sf_mex_destroy(&c4_lhs64);
  sf_mex_destroy(&c4_rhs65);
  sf_mex_destroy(&c4_lhs65);
  sf_mex_destroy(&c4_rhs66);
  sf_mex_destroy(&c4_lhs66);
  sf_mex_destroy(&c4_rhs67);
  sf_mex_destroy(&c4_lhs67);
  sf_mex_destroy(&c4_rhs68);
  sf_mex_destroy(&c4_lhs68);
  sf_mex_destroy(&c4_rhs69);
  sf_mex_destroy(&c4_lhs69);
  sf_mex_destroy(&c4_rhs70);
  sf_mex_destroy(&c4_lhs70);
  sf_mex_destroy(&c4_rhs71);
  sf_mex_destroy(&c4_lhs71);
  sf_mex_destroy(&c4_rhs72);
  sf_mex_destroy(&c4_lhs72);
  sf_mex_destroy(&c4_rhs73);
  sf_mex_destroy(&c4_lhs73);
  sf_mex_destroy(&c4_rhs74);
  sf_mex_destroy(&c4_lhs74);
  sf_mex_destroy(&c4_rhs75);
  sf_mex_destroy(&c4_lhs75);
  sf_mex_destroy(&c4_rhs76);
  sf_mex_destroy(&c4_lhs76);
  sf_mex_destroy(&c4_rhs77);
  sf_mex_destroy(&c4_lhs77);
  sf_mex_destroy(&c4_rhs78);
  sf_mex_destroy(&c4_lhs78);
  sf_mex_destroy(&c4_rhs79);
  sf_mex_destroy(&c4_lhs79);
  sf_mex_destroy(&c4_rhs80);
  sf_mex_destroy(&c4_lhs80);
  sf_mex_destroy(&c4_rhs81);
  sf_mex_destroy(&c4_lhs81);
  sf_mex_destroy(&c4_rhs82);
  sf_mex_destroy(&c4_lhs82);
  sf_mex_destroy(&c4_rhs83);
  sf_mex_destroy(&c4_lhs83);
  sf_mex_destroy(&c4_rhs84);
  sf_mex_destroy(&c4_lhs84);
  sf_mex_destroy(&c4_rhs85);
  sf_mex_destroy(&c4_lhs85);
  sf_mex_destroy(&c4_rhs86);
  sf_mex_destroy(&c4_lhs86);
  sf_mex_destroy(&c4_rhs87);
  sf_mex_destroy(&c4_lhs87);
  sf_mex_destroy(&c4_rhs88);
  sf_mex_destroy(&c4_lhs88);
  sf_mex_destroy(&c4_rhs89);
  sf_mex_destroy(&c4_lhs89);
  sf_mex_destroy(&c4_rhs90);
  sf_mex_destroy(&c4_lhs90);
  sf_mex_destroy(&c4_rhs91);
  sf_mex_destroy(&c4_lhs91);
  sf_mex_destroy(&c4_rhs92);
  sf_mex_destroy(&c4_lhs92);
  sf_mex_destroy(&c4_rhs93);
  sf_mex_destroy(&c4_lhs93);
  sf_mex_destroy(&c4_rhs94);
  sf_mex_destroy(&c4_lhs94);
  sf_mex_destroy(&c4_rhs95);
  sf_mex_destroy(&c4_lhs95);
  sf_mex_destroy(&c4_rhs96);
  sf_mex_destroy(&c4_lhs96);
  sf_mex_destroy(&c4_rhs97);
  sf_mex_destroy(&c4_lhs97);
  sf_mex_destroy(&c4_rhs98);
  sf_mex_destroy(&c4_lhs98);
  sf_mex_destroy(&c4_rhs99);
  sf_mex_destroy(&c4_lhs99);
  sf_mex_destroy(&c4_rhs100);
  sf_mex_destroy(&c4_lhs100);
  sf_mex_destroy(&c4_rhs101);
  sf_mex_destroy(&c4_lhs101);
  sf_mex_destroy(&c4_rhs102);
  sf_mex_destroy(&c4_lhs102);
  sf_mex_destroy(&c4_rhs103);
  sf_mex_destroy(&c4_lhs103);
  sf_mex_destroy(&c4_rhs104);
  sf_mex_destroy(&c4_lhs104);
  sf_mex_destroy(&c4_rhs105);
  sf_mex_destroy(&c4_lhs105);
  sf_mex_destroy(&c4_rhs106);
  sf_mex_destroy(&c4_lhs106);
  sf_mex_destroy(&c4_rhs107);
  sf_mex_destroy(&c4_lhs107);
  sf_mex_destroy(&c4_rhs108);
  sf_mex_destroy(&c4_lhs108);
  sf_mex_destroy(&c4_rhs109);
  sf_mex_destroy(&c4_lhs109);
  sf_mex_destroy(&c4_rhs110);
  sf_mex_destroy(&c4_lhs110);
  sf_mex_destroy(&c4_rhs111);
  sf_mex_destroy(&c4_lhs111);
  sf_mex_destroy(&c4_rhs112);
  sf_mex_destroy(&c4_lhs112);
  sf_mex_destroy(&c4_rhs113);
  sf_mex_destroy(&c4_lhs113);
  sf_mex_destroy(&c4_rhs114);
  sf_mex_destroy(&c4_lhs114);
  sf_mex_destroy(&c4_rhs115);
  sf_mex_destroy(&c4_lhs115);
  sf_mex_destroy(&c4_rhs116);
  sf_mex_destroy(&c4_lhs116);
  sf_mex_destroy(&c4_rhs117);
  sf_mex_destroy(&c4_lhs117);
  sf_mex_destroy(&c4_rhs118);
  sf_mex_destroy(&c4_lhs118);
  sf_mex_destroy(&c4_rhs119);
  sf_mex_destroy(&c4_lhs119);
  sf_mex_destroy(&c4_rhs120);
  sf_mex_destroy(&c4_lhs120);
  sf_mex_destroy(&c4_rhs121);
  sf_mex_destroy(&c4_lhs121);
  sf_mex_destroy(&c4_rhs122);
  sf_mex_destroy(&c4_lhs122);
  sf_mex_destroy(&c4_rhs123);
  sf_mex_destroy(&c4_lhs123);
  sf_mex_destroy(&c4_rhs124);
  sf_mex_destroy(&c4_lhs124);
  sf_mex_destroy(&c4_rhs125);
  sf_mex_destroy(&c4_lhs125);
  sf_mex_destroy(&c4_rhs126);
  sf_mex_destroy(&c4_lhs126);
  sf_mex_destroy(&c4_rhs127);
  sf_mex_destroy(&c4_lhs127);
}

static void c4_c_info_helper(const mxArray **c4_info)
{
  const mxArray *c4_rhs128 = NULL;
  const mxArray *c4_lhs128 = NULL;
  const mxArray *c4_rhs129 = NULL;
  const mxArray *c4_lhs129 = NULL;
  const mxArray *c4_rhs130 = NULL;
  const mxArray *c4_lhs130 = NULL;
  const mxArray *c4_rhs131 = NULL;
  const mxArray *c4_lhs131 = NULL;
  const mxArray *c4_rhs132 = NULL;
  const mxArray *c4_lhs132 = NULL;
  const mxArray *c4_rhs133 = NULL;
  const mxArray *c4_lhs133 = NULL;
  const mxArray *c4_rhs134 = NULL;
  const mxArray *c4_lhs134 = NULL;
  const mxArray *c4_rhs135 = NULL;
  const mxArray *c4_lhs135 = NULL;
  const mxArray *c4_rhs136 = NULL;
  const mxArray *c4_lhs136 = NULL;
  const mxArray *c4_rhs137 = NULL;
  const mxArray *c4_lhs137 = NULL;
  const mxArray *c4_rhs138 = NULL;
  const mxArray *c4_lhs138 = NULL;
  const mxArray *c4_rhs139 = NULL;
  const mxArray *c4_lhs139 = NULL;
  const mxArray *c4_rhs140 = NULL;
  const mxArray *c4_lhs140 = NULL;
  const mxArray *c4_rhs141 = NULL;
  const mxArray *c4_lhs141 = NULL;
  const mxArray *c4_rhs142 = NULL;
  const mxArray *c4_lhs142 = NULL;
  const mxArray *c4_rhs143 = NULL;
  const mxArray *c4_lhs143 = NULL;
  const mxArray *c4_rhs144 = NULL;
  const mxArray *c4_lhs144 = NULL;
  const mxArray *c4_rhs145 = NULL;
  const mxArray *c4_lhs145 = NULL;
  const mxArray *c4_rhs146 = NULL;
  const mxArray *c4_lhs146 = NULL;
  const mxArray *c4_rhs147 = NULL;
  const mxArray *c4_lhs147 = NULL;
  const mxArray *c4_rhs148 = NULL;
  const mxArray *c4_lhs148 = NULL;
  const mxArray *c4_rhs149 = NULL;
  const mxArray *c4_lhs149 = NULL;
  const mxArray *c4_rhs150 = NULL;
  const mxArray *c4_lhs150 = NULL;
  const mxArray *c4_rhs151 = NULL;
  const mxArray *c4_lhs151 = NULL;
  const mxArray *c4_rhs152 = NULL;
  const mxArray *c4_lhs152 = NULL;
  const mxArray *c4_rhs153 = NULL;
  const mxArray *c4_lhs153 = NULL;
  const mxArray *c4_rhs154 = NULL;
  const mxArray *c4_lhs154 = NULL;
  const mxArray *c4_rhs155 = NULL;
  const mxArray *c4_lhs155 = NULL;
  const mxArray *c4_rhs156 = NULL;
  const mxArray *c4_lhs156 = NULL;
  const mxArray *c4_rhs157 = NULL;
  const mxArray *c4_lhs157 = NULL;
  const mxArray *c4_rhs158 = NULL;
  const mxArray *c4_lhs158 = NULL;
  const mxArray *c4_rhs159 = NULL;
  const mxArray *c4_lhs159 = NULL;
  const mxArray *c4_rhs160 = NULL;
  const mxArray *c4_lhs160 = NULL;
  const mxArray *c4_rhs161 = NULL;
  const mxArray *c4_lhs161 = NULL;
  const mxArray *c4_rhs162 = NULL;
  const mxArray *c4_lhs162 = NULL;
  const mxArray *c4_rhs163 = NULL;
  const mxArray *c4_lhs163 = NULL;
  const mxArray *c4_rhs164 = NULL;
  const mxArray *c4_lhs164 = NULL;
  const mxArray *c4_rhs165 = NULL;
  const mxArray *c4_lhs165 = NULL;
  const mxArray *c4_rhs166 = NULL;
  const mxArray *c4_lhs166 = NULL;
  const mxArray *c4_rhs167 = NULL;
  const mxArray *c4_lhs167 = NULL;
  const mxArray *c4_rhs168 = NULL;
  const mxArray *c4_lhs168 = NULL;
  const mxArray *c4_rhs169 = NULL;
  const mxArray *c4_lhs169 = NULL;
  const mxArray *c4_rhs170 = NULL;
  const mxArray *c4_lhs170 = NULL;
  const mxArray *c4_rhs171 = NULL;
  const mxArray *c4_lhs171 = NULL;
  const mxArray *c4_rhs172 = NULL;
  const mxArray *c4_lhs172 = NULL;
  const mxArray *c4_rhs173 = NULL;
  const mxArray *c4_lhs173 = NULL;
  const mxArray *c4_rhs174 = NULL;
  const mxArray *c4_lhs174 = NULL;
  const mxArray *c4_rhs175 = NULL;
  const mxArray *c4_lhs175 = NULL;
  const mxArray *c4_rhs176 = NULL;
  const mxArray *c4_lhs176 = NULL;
  const mxArray *c4_rhs177 = NULL;
  const mxArray *c4_lhs177 = NULL;
  const mxArray *c4_rhs178 = NULL;
  const mxArray *c4_lhs178 = NULL;
  const mxArray *c4_rhs179 = NULL;
  const mxArray *c4_lhs179 = NULL;
  const mxArray *c4_rhs180 = NULL;
  const mxArray *c4_lhs180 = NULL;
  const mxArray *c4_rhs181 = NULL;
  const mxArray *c4_lhs181 = NULL;
  const mxArray *c4_rhs182 = NULL;
  const mxArray *c4_lhs182 = NULL;
  const mxArray *c4_rhs183 = NULL;
  const mxArray *c4_lhs183 = NULL;
  const mxArray *c4_rhs184 = NULL;
  const mxArray *c4_lhs184 = NULL;
  const mxArray *c4_rhs185 = NULL;
  const mxArray *c4_lhs185 = NULL;
  const mxArray *c4_rhs186 = NULL;
  const mxArray *c4_lhs186 = NULL;
  const mxArray *c4_rhs187 = NULL;
  const mxArray *c4_lhs187 = NULL;
  const mxArray *c4_rhs188 = NULL;
  const mxArray *c4_lhs188 = NULL;
  const mxArray *c4_rhs189 = NULL;
  const mxArray *c4_lhs189 = NULL;
  const mxArray *c4_rhs190 = NULL;
  const mxArray *c4_lhs190 = NULL;
  const mxArray *c4_rhs191 = NULL;
  const mxArray *c4_lhs191 = NULL;
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 128);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isnan"), "name", "name", 128);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 128);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 128);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c4_rhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs128), "rhs", "rhs",
                  128);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs128), "lhs", "lhs",
                  128);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 129);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 129);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 129);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 129);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c4_rhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs129), "rhs", "rhs",
                  129);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs129), "lhs", "lhs",
                  129);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 130);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isinf"), "name", "name", 130);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 130);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 130);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c4_rhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs130), "rhs", "rhs",
                  130);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs130), "lhs", "lhs",
                  130);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 131);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_guarded_inf"), "name",
                  "name", 131);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_inf.m"),
                  "resolved", "resolved", 131);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c4_rhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs131), "rhs", "rhs",
                  131);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs131), "lhs", "lhs",
                  131);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_inf.m"),
                  "context", "context", 132);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 132);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 132);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 132);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c4_rhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs132), "rhs", "rhs",
                  132);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs132), "lhs", "lhs",
                  132);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 133);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("realmax"), "name", "name", 133);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 133);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m"), "resolved",
                  "resolved", 133);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c4_rhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs133), "rhs", "rhs",
                  133);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs133), "lhs", "lhs",
                  133);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m"), "context",
                  "context", 134);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_realmax"), "name", "name",
                  134);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmax.m"), "resolved",
                  "resolved", 134);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c4_rhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs134), "rhs", "rhs",
                  134);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs134), "lhs", "lhs",
                  134);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmax.m"), "context",
                  "context", 135);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 135);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 135);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c4_rhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs135), "rhs", "rhs",
                  135);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs135), "lhs", "lhs",
                  135);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 136);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mrdivide"), "name", "name",
                  136);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 136);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c4_rhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs136), "rhs", "rhs",
                  136);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs136), "lhs", "lhs",
                  136);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 137);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 137);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 137);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c4_rhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs137), "rhs", "rhs",
                  137);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs137), "lhs", "lhs",
                  137);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 138);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("rdivide"), "name", "name", 138);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 138);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c4_rhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs138), "rhs", "rhs",
                  138);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs138), "lhs", "lhs",
                  138);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m!eml_complex_scalar_sqrt"),
                  "context", "context", 139);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_hypot"), "name",
                  "name", 139);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "resolved", "resolved", 139);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c4_rhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs139), "rhs", "rhs",
                  139);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs139), "lhs", "lhs",
                  139);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "context", "context", 140);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 140);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 140);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c4_rhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs140), "rhs", "rhs",
                  140);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs140), "lhs", "lhs",
                  140);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "context", "context", 141);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  141);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 141);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 141);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1350417854U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c4_rhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs141), "rhs", "rhs",
                  141);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs141), "lhs", "lhs",
                  141);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 142);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 142);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 142);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c4_rhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs142), "rhs", "rhs",
                  142);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs142), "lhs", "lhs",
                  142);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zhgeqz.m"),
                  "context", "context", 143);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_zrot_rows"), "name",
                  "name", 143);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 143);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_zrot_rows.m"),
                  "resolved", "resolved", 143);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c4_rhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs143), "rhs", "rhs",
                  143);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs143), "lhs", "lhs",
                  143);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 144);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_div"), "name", "name", 144);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 144);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c4_rhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs144), "rhs", "rhs",
                  144);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs144), "lhs", "lhs",
                  144);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p!equalsize"),
                  "context", "context", 145);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("max"), "name", "name", 145);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "resolved",
                  "resolved", 145);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1311262516U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c4_rhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs145), "rhs", "rhs",
                  145);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs145), "lhs", "lhs",
                  145);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "context",
                  "context", 146);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 146);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 146);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1378303184U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c4_rhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs146), "rhs", "rhs",
                  146);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs146), "lhs", "lhs",
                  146);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 147);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 147);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 147);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c4_rhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs147), "rhs", "rhs",
                  147);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs147), "lhs", "lhs",
                  147);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 148);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 148);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 148);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c4_rhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs148), "rhs", "rhs",
                  148);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs148), "lhs", "lhs",
                  148);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 149);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 149);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 149);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c4_rhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs149), "rhs", "rhs",
                  149);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs149), "lhs", "lhs",
                  149);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 150);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 150);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 150);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c4_rhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs150), "rhs", "rhs",
                  150);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs150), "lhs", "lhs",
                  150);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 151);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 151);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 151);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c4_rhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs151), "rhs", "rhs",
                  151);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs151), "lhs", "lhs",
                  151);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 152);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_warning"), "name", "name",
                  152);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 152);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c4_rhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs152), "rhs", "rhs",
                  152);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs152), "lhs", "lhs",
                  152);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 153);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mpower"), "name", "name", 153);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 153);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c4_rhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs153), "rhs", "rhs",
                  153);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs153), "lhs", "lhs",
                  153);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 154);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("inv"), "name", "name", 154);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 154);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 154);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c4_rhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs154), "rhs", "rhs",
                  154);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs154), "lhs", "lhs",
                  154);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 155);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 155);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 155);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c4_rhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs155), "rhs", "rhs",
                  155);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs155), "lhs", "lhs",
                  155);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 156);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 156);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 156);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c4_rhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs156), "rhs", "rhs",
                  156);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs156), "lhs", "lhs",
                  156);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 157);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_div"), "name", "name", 157);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 157);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c4_rhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs157), "rhs", "rhs",
                  157);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs157), "lhs", "lhs",
                  157);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 158);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 158);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 158);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 158);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c4_rhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs158), "rhs", "rhs",
                  158);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs158), "lhs", "lhs",
                  158);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 159);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("norm"), "name", "name", 159);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 159);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 159);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c4_rhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs159), "rhs", "rhs",
                  159);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs159), "lhs", "lhs",
                  159);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 160);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 160);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 160);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 160);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c4_rhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs160), "rhs", "rhs",
                  160);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs160), "lhs", "lhs",
                  160);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 161);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 161);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 161);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 161);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c4_rhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs161), "rhs", "rhs",
                  161);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs161), "lhs", "lhs",
                  161);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 162);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isnan"), "name", "name", 162);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 162);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 162);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c4_rhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs162), "rhs", "rhs",
                  162);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs162), "lhs", "lhs",
                  162);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 163);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 163);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 163);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825976U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c4_rhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs163), "rhs", "rhs",
                  163);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs163), "lhs", "lhs",
                  163);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 164);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_warning"), "name", "name",
                  164);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 164);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 164);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 164);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 164);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 164);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 164);
  sf_mex_assign(&c4_rhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs164), "rhs", "rhs",
                  164);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs164), "lhs", "lhs",
                  164);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 165);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isnan"), "name", "name", 165);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 165);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 165);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 165);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 165);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 165);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 165);
  sf_mex_assign(&c4_rhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs165), "rhs", "rhs",
                  165);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs165), "lhs", "lhs",
                  165);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 166);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eps"), "name", "name", 166);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 166);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 166);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 166);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 166);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 166);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 166);
  sf_mex_assign(&c4_rhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs166), "rhs", "rhs",
                  166);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs166), "lhs", "lhs",
                  166);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 167);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_flt2str"), "name", "name",
                  167);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 167);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "resolved",
                  "resolved", 167);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 167);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 167);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 167);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 167);
  sf_mex_assign(&c4_rhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs167), "rhs", "rhs",
                  167);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs167), "lhs", "lhs",
                  167);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "context",
                  "context", 168);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "name", "name", 168);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 168);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m"), "resolved",
                  "resolved", 168);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1319737168U), "fileTimeLo",
                  "fileTimeLo", 168);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 168);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 168);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 168);
  sf_mex_assign(&c4_rhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs168), "rhs", "rhs",
                  168);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs168), "lhs", "lhs",
                  168);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 169);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("exp"), "name", "name", 169);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 169);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m"), "resolved",
                  "resolved", 169);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837580U), "fileTimeLo",
                  "fileTimeLo", 169);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 169);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 169);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 169);
  sf_mex_assign(&c4_rhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs169), "rhs", "rhs",
                  169);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs169), "lhs", "lhs",
                  169);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/exp.m"), "context",
                  "context", 170);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_exp"), "name",
                  "name", 170);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 170);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_exp.m"),
                  "resolved", "resolved", 170);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1301335664U), "fileTimeLo",
                  "fileTimeLo", 170);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 170);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 170);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 170);
  sf_mex_assign(&c4_rhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs170), "rhs", "rhs",
                  170);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs170), "lhs", "lhs",
                  170);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 171);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 171);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 171);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 171);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 171);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 171);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 171);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 171);
  sf_mex_assign(&c4_rhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs171), "rhs", "rhs",
                  171);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs171), "lhs", "lhs",
                  171);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 172);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 172);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 172);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 172);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 172);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 172);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 172);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 172);
  sf_mex_assign(&c4_rhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs172), "rhs", "rhs",
                  172);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs172), "lhs", "lhs",
                  172);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 173);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 173);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 173);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 173);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 173);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 173);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 173);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 173);
  sf_mex_assign(&c4_rhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs173), "rhs", "rhs",
                  173);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs173), "lhs", "lhs",
                  173);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 174);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 174);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 174);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 174);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 174);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 174);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 174);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 174);
  sf_mex_assign(&c4_rhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs174), "rhs", "rhs",
                  174);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs174), "lhs", "lhs",
                  174);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 175);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 175);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 175);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 175);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 175);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 175);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 175);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 175);
  sf_mex_assign(&c4_rhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs175), "rhs", "rhs",
                  175);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs175), "lhs", "lhs",
                  175);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 176);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 176);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 176);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 176);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 176);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 176);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 176);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 176);
  sf_mex_assign(&c4_rhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs176), "rhs", "rhs",
                  176);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs176), "lhs", "lhs",
                  176);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 177);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 177);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 177);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 177);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 177);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 177);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 177);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 177);
  sf_mex_assign(&c4_rhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs177), "rhs", "rhs",
                  177);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs177), "lhs", "lhs",
                  177);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 178);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  178);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 178);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 178);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 178);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 178);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 178);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 178);
  sf_mex_assign(&c4_rhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs178), "rhs", "rhs",
                  178);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs178), "lhs", "lhs",
                  178);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 179);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 179);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 179);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 179);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 179);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 179);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 179);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 179);
  sf_mex_assign(&c4_rhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs179), "rhs", "rhs",
                  179);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs179), "lhs", "lhs",
                  179);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 180);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 180);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 180);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 180);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 180);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 180);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 180);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 180);
  sf_mex_assign(&c4_rhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs180), "rhs", "rhs",
                  180);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs180), "lhs", "lhs",
                  180);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 181);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 181);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 181);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 181);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 181);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 181);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 181);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 181);
  sf_mex_assign(&c4_rhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs181), "rhs", "rhs",
                  181);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs181), "lhs", "lhs",
                  181);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 182);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 182);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 182);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 182);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 182);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 182);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 182);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 182);
  sf_mex_assign(&c4_rhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs182), "rhs", "rhs",
                  182);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs182), "lhs", "lhs",
                  182);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 183);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 183);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 183);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 183);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 183);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 183);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 183);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 183);
  sf_mex_assign(&c4_rhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs183), "rhs", "rhs",
                  183);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs183), "lhs", "lhs",
                  183);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 184);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 184);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 184);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 184);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 184);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 184);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 184);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 184);
  sf_mex_assign(&c4_rhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs184), "rhs", "rhs",
                  184);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs184), "lhs", "lhs",
                  184);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 185);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 185);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 185);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 185);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 185);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 185);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 185);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 185);
  sf_mex_assign(&c4_rhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs185), "rhs", "rhs",
                  185);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs185), "lhs", "lhs",
                  185);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_integer_power"),
                  "context", "context", 186);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("inv"), "name", "name", 186);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 186);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 186);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 186);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 186);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 186);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 186);
  sf_mex_assign(&c4_rhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs186), "rhs", "rhs",
                  186);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs186), "lhs", "lhs",
                  186);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 187);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 187);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 187);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 187);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 187);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 187);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 187);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 187);
  sf_mex_assign(&c4_rhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs187), "rhs", "rhs",
                  187);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs187), "lhs", "lhs",
                  187);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 188);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_error"), "name", "name",
                  188);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 188);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 188);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 188);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 188);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 188);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 188);
  sf_mex_assign(&c4_rhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs188), "rhs", "rhs",
                  188);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs188), "lhs", "lhs",
                  188);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 189);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eig"), "name", "name", 189);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 189);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "resolved",
                  "resolved", 189);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1305325200U), "fileTimeLo",
                  "fileTimeLo", 189);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 189);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 189);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 189);
  sf_mex_assign(&c4_rhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs189), "rhs", "rhs",
                  189);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs189), "lhs", "lhs",
                  189);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 190);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 190);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 190);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 190);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 190);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 190);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 190);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 190);
  sf_mex_assign(&c4_rhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs190), "rhs", "rhs",
                  190);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs190), "lhs", "lhs",
                  190);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgghrd.m"),
                  "context", "context", 191);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eye"), "name", "name", 191);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 191);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 191);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 191);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 191);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 191);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 191);
  sf_mex_assign(&c4_rhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs191), "rhs", "rhs",
                  191);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs191), "lhs", "lhs",
                  191);
  sf_mex_destroy(&c4_rhs128);
  sf_mex_destroy(&c4_lhs128);
  sf_mex_destroy(&c4_rhs129);
  sf_mex_destroy(&c4_lhs129);
  sf_mex_destroy(&c4_rhs130);
  sf_mex_destroy(&c4_lhs130);
  sf_mex_destroy(&c4_rhs131);
  sf_mex_destroy(&c4_lhs131);
  sf_mex_destroy(&c4_rhs132);
  sf_mex_destroy(&c4_lhs132);
  sf_mex_destroy(&c4_rhs133);
  sf_mex_destroy(&c4_lhs133);
  sf_mex_destroy(&c4_rhs134);
  sf_mex_destroy(&c4_lhs134);
  sf_mex_destroy(&c4_rhs135);
  sf_mex_destroy(&c4_lhs135);
  sf_mex_destroy(&c4_rhs136);
  sf_mex_destroy(&c4_lhs136);
  sf_mex_destroy(&c4_rhs137);
  sf_mex_destroy(&c4_lhs137);
  sf_mex_destroy(&c4_rhs138);
  sf_mex_destroy(&c4_lhs138);
  sf_mex_destroy(&c4_rhs139);
  sf_mex_destroy(&c4_lhs139);
  sf_mex_destroy(&c4_rhs140);
  sf_mex_destroy(&c4_lhs140);
  sf_mex_destroy(&c4_rhs141);
  sf_mex_destroy(&c4_lhs141);
  sf_mex_destroy(&c4_rhs142);
  sf_mex_destroy(&c4_lhs142);
  sf_mex_destroy(&c4_rhs143);
  sf_mex_destroy(&c4_lhs143);
  sf_mex_destroy(&c4_rhs144);
  sf_mex_destroy(&c4_lhs144);
  sf_mex_destroy(&c4_rhs145);
  sf_mex_destroy(&c4_lhs145);
  sf_mex_destroy(&c4_rhs146);
  sf_mex_destroy(&c4_lhs146);
  sf_mex_destroy(&c4_rhs147);
  sf_mex_destroy(&c4_lhs147);
  sf_mex_destroy(&c4_rhs148);
  sf_mex_destroy(&c4_lhs148);
  sf_mex_destroy(&c4_rhs149);
  sf_mex_destroy(&c4_lhs149);
  sf_mex_destroy(&c4_rhs150);
  sf_mex_destroy(&c4_lhs150);
  sf_mex_destroy(&c4_rhs151);
  sf_mex_destroy(&c4_lhs151);
  sf_mex_destroy(&c4_rhs152);
  sf_mex_destroy(&c4_lhs152);
  sf_mex_destroy(&c4_rhs153);
  sf_mex_destroy(&c4_lhs153);
  sf_mex_destroy(&c4_rhs154);
  sf_mex_destroy(&c4_lhs154);
  sf_mex_destroy(&c4_rhs155);
  sf_mex_destroy(&c4_lhs155);
  sf_mex_destroy(&c4_rhs156);
  sf_mex_destroy(&c4_lhs156);
  sf_mex_destroy(&c4_rhs157);
  sf_mex_destroy(&c4_lhs157);
  sf_mex_destroy(&c4_rhs158);
  sf_mex_destroy(&c4_lhs158);
  sf_mex_destroy(&c4_rhs159);
  sf_mex_destroy(&c4_lhs159);
  sf_mex_destroy(&c4_rhs160);
  sf_mex_destroy(&c4_lhs160);
  sf_mex_destroy(&c4_rhs161);
  sf_mex_destroy(&c4_lhs161);
  sf_mex_destroy(&c4_rhs162);
  sf_mex_destroy(&c4_lhs162);
  sf_mex_destroy(&c4_rhs163);
  sf_mex_destroy(&c4_lhs163);
  sf_mex_destroy(&c4_rhs164);
  sf_mex_destroy(&c4_lhs164);
  sf_mex_destroy(&c4_rhs165);
  sf_mex_destroy(&c4_lhs165);
  sf_mex_destroy(&c4_rhs166);
  sf_mex_destroy(&c4_lhs166);
  sf_mex_destroy(&c4_rhs167);
  sf_mex_destroy(&c4_lhs167);
  sf_mex_destroy(&c4_rhs168);
  sf_mex_destroy(&c4_lhs168);
  sf_mex_destroy(&c4_rhs169);
  sf_mex_destroy(&c4_lhs169);
  sf_mex_destroy(&c4_rhs170);
  sf_mex_destroy(&c4_lhs170);
  sf_mex_destroy(&c4_rhs171);
  sf_mex_destroy(&c4_lhs171);
  sf_mex_destroy(&c4_rhs172);
  sf_mex_destroy(&c4_lhs172);
  sf_mex_destroy(&c4_rhs173);
  sf_mex_destroy(&c4_lhs173);
  sf_mex_destroy(&c4_rhs174);
  sf_mex_destroy(&c4_lhs174);
  sf_mex_destroy(&c4_rhs175);
  sf_mex_destroy(&c4_lhs175);
  sf_mex_destroy(&c4_rhs176);
  sf_mex_destroy(&c4_lhs176);
  sf_mex_destroy(&c4_rhs177);
  sf_mex_destroy(&c4_lhs177);
  sf_mex_destroy(&c4_rhs178);
  sf_mex_destroy(&c4_lhs178);
  sf_mex_destroy(&c4_rhs179);
  sf_mex_destroy(&c4_lhs179);
  sf_mex_destroy(&c4_rhs180);
  sf_mex_destroy(&c4_lhs180);
  sf_mex_destroy(&c4_rhs181);
  sf_mex_destroy(&c4_lhs181);
  sf_mex_destroy(&c4_rhs182);
  sf_mex_destroy(&c4_lhs182);
  sf_mex_destroy(&c4_rhs183);
  sf_mex_destroy(&c4_lhs183);
  sf_mex_destroy(&c4_rhs184);
  sf_mex_destroy(&c4_lhs184);
  sf_mex_destroy(&c4_rhs185);
  sf_mex_destroy(&c4_lhs185);
  sf_mex_destroy(&c4_rhs186);
  sf_mex_destroy(&c4_lhs186);
  sf_mex_destroy(&c4_rhs187);
  sf_mex_destroy(&c4_lhs187);
  sf_mex_destroy(&c4_rhs188);
  sf_mex_destroy(&c4_lhs188);
  sf_mex_destroy(&c4_rhs189);
  sf_mex_destroy(&c4_lhs189);
  sf_mex_destroy(&c4_rhs190);
  sf_mex_destroy(&c4_lhs190);
  sf_mex_destroy(&c4_rhs191);
  sf_mex_destroy(&c4_lhs191);
}

static void c4_d_info_helper(const mxArray **c4_info)
{
  const mxArray *c4_rhs192 = NULL;
  const mxArray *c4_lhs192 = NULL;
  const mxArray *c4_rhs193 = NULL;
  const mxArray *c4_lhs193 = NULL;
  const mxArray *c4_rhs194 = NULL;
  const mxArray *c4_lhs194 = NULL;
  const mxArray *c4_rhs195 = NULL;
  const mxArray *c4_lhs195 = NULL;
  const mxArray *c4_rhs196 = NULL;
  const mxArray *c4_lhs196 = NULL;
  const mxArray *c4_rhs197 = NULL;
  const mxArray *c4_lhs197 = NULL;
  const mxArray *c4_rhs198 = NULL;
  const mxArray *c4_lhs198 = NULL;
  const mxArray *c4_rhs199 = NULL;
  const mxArray *c4_lhs199 = NULL;
  const mxArray *c4_rhs200 = NULL;
  const mxArray *c4_lhs200 = NULL;
  const mxArray *c4_rhs201 = NULL;
  const mxArray *c4_lhs201 = NULL;
  const mxArray *c4_rhs202 = NULL;
  const mxArray *c4_lhs202 = NULL;
  const mxArray *c4_rhs203 = NULL;
  const mxArray *c4_lhs203 = NULL;
  const mxArray *c4_rhs204 = NULL;
  const mxArray *c4_lhs204 = NULL;
  const mxArray *c4_rhs205 = NULL;
  const mxArray *c4_lhs205 = NULL;
  const mxArray *c4_rhs206 = NULL;
  const mxArray *c4_lhs206 = NULL;
  const mxArray *c4_rhs207 = NULL;
  const mxArray *c4_lhs207 = NULL;
  const mxArray *c4_rhs208 = NULL;
  const mxArray *c4_lhs208 = NULL;
  const mxArray *c4_rhs209 = NULL;
  const mxArray *c4_lhs209 = NULL;
  const mxArray *c4_rhs210 = NULL;
  const mxArray *c4_lhs210 = NULL;
  const mxArray *c4_rhs211 = NULL;
  const mxArray *c4_lhs211 = NULL;
  const mxArray *c4_rhs212 = NULL;
  const mxArray *c4_lhs212 = NULL;
  const mxArray *c4_rhs213 = NULL;
  const mxArray *c4_lhs213 = NULL;
  const mxArray *c4_rhs214 = NULL;
  const mxArray *c4_lhs214 = NULL;
  const mxArray *c4_rhs215 = NULL;
  const mxArray *c4_lhs215 = NULL;
  const mxArray *c4_rhs216 = NULL;
  const mxArray *c4_lhs216 = NULL;
  const mxArray *c4_rhs217 = NULL;
  const mxArray *c4_lhs217 = NULL;
  const mxArray *c4_rhs218 = NULL;
  const mxArray *c4_lhs218 = NULL;
  const mxArray *c4_rhs219 = NULL;
  const mxArray *c4_lhs219 = NULL;
  const mxArray *c4_rhs220 = NULL;
  const mxArray *c4_lhs220 = NULL;
  const mxArray *c4_rhs221 = NULL;
  const mxArray *c4_lhs221 = NULL;
  const mxArray *c4_rhs222 = NULL;
  const mxArray *c4_lhs222 = NULL;
  const mxArray *c4_rhs223 = NULL;
  const mxArray *c4_lhs223 = NULL;
  const mxArray *c4_rhs224 = NULL;
  const mxArray *c4_lhs224 = NULL;
  const mxArray *c4_rhs225 = NULL;
  const mxArray *c4_lhs225 = NULL;
  const mxArray *c4_rhs226 = NULL;
  const mxArray *c4_lhs226 = NULL;
  const mxArray *c4_rhs227 = NULL;
  const mxArray *c4_lhs227 = NULL;
  const mxArray *c4_rhs228 = NULL;
  const mxArray *c4_lhs228 = NULL;
  const mxArray *c4_rhs229 = NULL;
  const mxArray *c4_lhs229 = NULL;
  const mxArray *c4_rhs230 = NULL;
  const mxArray *c4_lhs230 = NULL;
  const mxArray *c4_rhs231 = NULL;
  const mxArray *c4_lhs231 = NULL;
  const mxArray *c4_rhs232 = NULL;
  const mxArray *c4_lhs232 = NULL;
  const mxArray *c4_rhs233 = NULL;
  const mxArray *c4_lhs233 = NULL;
  const mxArray *c4_rhs234 = NULL;
  const mxArray *c4_lhs234 = NULL;
  const mxArray *c4_rhs235 = NULL;
  const mxArray *c4_lhs235 = NULL;
  const mxArray *c4_rhs236 = NULL;
  const mxArray *c4_lhs236 = NULL;
  const mxArray *c4_rhs237 = NULL;
  const mxArray *c4_lhs237 = NULL;
  const mxArray *c4_rhs238 = NULL;
  const mxArray *c4_lhs238 = NULL;
  const mxArray *c4_rhs239 = NULL;
  const mxArray *c4_lhs239 = NULL;
  const mxArray *c4_rhs240 = NULL;
  const mxArray *c4_lhs240 = NULL;
  const mxArray *c4_rhs241 = NULL;
  const mxArray *c4_lhs241 = NULL;
  const mxArray *c4_rhs242 = NULL;
  const mxArray *c4_lhs242 = NULL;
  const mxArray *c4_rhs243 = NULL;
  const mxArray *c4_lhs243 = NULL;
  const mxArray *c4_rhs244 = NULL;
  const mxArray *c4_lhs244 = NULL;
  const mxArray *c4_rhs245 = NULL;
  const mxArray *c4_lhs245 = NULL;
  const mxArray *c4_rhs246 = NULL;
  const mxArray *c4_lhs246 = NULL;
  const mxArray *c4_rhs247 = NULL;
  const mxArray *c4_lhs247 = NULL;
  const mxArray *c4_rhs248 = NULL;
  const mxArray *c4_lhs248 = NULL;
  const mxArray *c4_rhs249 = NULL;
  const mxArray *c4_lhs249 = NULL;
  const mxArray *c4_rhs250 = NULL;
  const mxArray *c4_lhs250 = NULL;
  const mxArray *c4_rhs251 = NULL;
  const mxArray *c4_lhs251 = NULL;
  const mxArray *c4_rhs252 = NULL;
  const mxArray *c4_lhs252 = NULL;
  const mxArray *c4_rhs253 = NULL;
  const mxArray *c4_lhs253 = NULL;
  const mxArray *c4_rhs254 = NULL;
  const mxArray *c4_lhs254 = NULL;
  const mxArray *c4_rhs255 = NULL;
  const mxArray *c4_lhs255 = NULL;
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 192);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 192);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 192);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 192);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 192);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 192);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 192);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 192);
  sf_mex_assign(&c4_rhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs192), "rhs", "rhs",
                  192);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs192), "lhs", "lhs",
                  192);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 193);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 193);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 193);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 193);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 193);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 193);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 193);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 193);
  sf_mex_assign(&c4_rhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs193), "rhs", "rhs",
                  193);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs193), "lhs", "lhs",
                  193);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 194);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 194);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 194);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 194);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 194);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 194);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 194);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 194);
  sf_mex_assign(&c4_rhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs194), "rhs", "rhs",
                  194);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs194), "lhs", "lhs",
                  194);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 195);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("intmax"), "name", "name", 195);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 195);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 195);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 195);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 195);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 195);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 195);
  sf_mex_assign(&c4_rhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs195), "rhs", "rhs",
                  195);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs195), "lhs", "lhs",
                  195);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 196);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("intmin"), "name", "name", 196);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 196);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 196);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 196);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 196);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 196);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 196);
  sf_mex_assign(&c4_rhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs196), "rhs", "rhs",
                  196);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs196), "lhs", "lhs",
                  196);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 197);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 197);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 197);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 197);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 197);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 197);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 197);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 197);
  sf_mex_assign(&c4_rhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs197), "rhs", "rhs",
                  197);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs197), "lhs", "lhs",
                  197);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 198);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 198);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 198);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 198);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 198);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 198);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 198);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 198);
  sf_mex_assign(&c4_rhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs198), "rhs", "rhs",
                  198);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs198), "lhs", "lhs",
                  198);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 199);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("intmax"), "name", "name", 199);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 199);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 199);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 199);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 199);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 199);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 199);
  sf_mex_assign(&c4_rhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs199), "rhs", "rhs",
                  199);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs199), "lhs", "lhs",
                  199);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 200);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 200);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 200);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 200);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 200);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 200);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 200);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 200);
  sf_mex_assign(&c4_rhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs200), "rhs", "rhs",
                  200);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs200), "lhs", "lhs",
                  200);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 201);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 201);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 201);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 201);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 201);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 201);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 201);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 201);
  sf_mex_assign(&c4_rhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs201), "rhs", "rhs",
                  201);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs201), "lhs", "lhs",
                  201);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 202);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_ztgevc"), "name",
                  "name", 202);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 202);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "resolved", "resolved", 202);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826024U), "fileTimeLo",
                  "fileTimeLo", 202);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 202);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 202);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 202);
  sf_mex_assign(&c4_rhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs202), "rhs", "rhs",
                  202);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs202), "lhs", "lhs",
                  202);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 203);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eps"), "name", "name", 203);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 203);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 203);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 203);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 203);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 203);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 203);
  sf_mex_assign(&c4_rhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs203), "rhs", "rhs",
                  203);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs203), "lhs", "lhs",
                  203);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 204);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("realmin"), "name", "name", 204);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 204);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 204);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 204);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 204);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 204);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 204);
  sf_mex_assign(&c4_rhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs204), "rhs", "rhs",
                  204);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs204), "lhs", "lhs",
                  204);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m!abs1"),
                  "context", "context", 205);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 205);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 205);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 205);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 205);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 205);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 205);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 205);
  sf_mex_assign(&c4_rhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs205), "rhs", "rhs",
                  205);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs205), "lhs", "lhs",
                  205);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 206);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 206);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 206);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 206);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 206);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 206);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 206);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 206);
  sf_mex_assign(&c4_rhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs206), "rhs", "rhs",
                  206);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs206), "lhs", "lhs",
                  206);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 207);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 207);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 207);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 207);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 207);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 207);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 207);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 207);
  sf_mex_assign(&c4_rhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs207), "rhs", "rhs",
                  207);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs207), "lhs", "lhs",
                  207);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_ztgevc.m"),
                  "context", "context", 208);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mrdivide"), "name", "name",
                  208);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 208);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 208);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 208);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 208);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 208);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 208);
  sf_mex_assign(&c4_rhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs208), "rhs", "rhs",
                  208);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs208), "lhs", "lhs",
                  208);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 209);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_matlab_zggbak"), "name",
                  "name", 209);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 209);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "resolved", "resolved", 209);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826018U), "fileTimeLo",
                  "fileTimeLo", 209);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 209);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 209);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 209);
  sf_mex_assign(&c4_rhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs209), "rhs", "rhs",
                  209);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs209), "lhs", "lhs",
                  209);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 210);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 210);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 210);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 210);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 210);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 210);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 210);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 210);
  sf_mex_assign(&c4_rhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs210), "rhs", "rhs",
                  210);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs210), "lhs", "lhs",
                  210);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 211);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 211);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 211);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 211);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 211);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 211);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 211);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 211);
  sf_mex_assign(&c4_rhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs211), "rhs", "rhs",
                  211);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs211), "lhs", "lhs",
                  211);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 212);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 212);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 212);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 212);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 212);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 212);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 212);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 212);
  sf_mex_assign(&c4_rhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs212), "rhs", "rhs",
                  212);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs212), "lhs", "lhs",
                  212);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggbak.m"),
                  "context", "context", 213);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 213);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 213);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 213);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 213);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 213);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 213);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 213);
  sf_mex_assign(&c4_rhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs213), "rhs", "rhs",
                  213);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs213), "lhs", "lhs",
                  213);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m!abs1"),
                  "context", "context", 214);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 214);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 214);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 214);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 214);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 214);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 214);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 214);
  sf_mex_assign(&c4_rhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs214), "rhs", "rhs",
                  214);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs214), "lhs", "lhs",
                  214);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zggev.m"),
                  "context", "context", 215);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 215);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 215);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 215);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 215);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 215);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 215);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 215);
  sf_mex_assign(&c4_rhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs215), "rhs", "rhs",
                  215);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs215), "lhs", "lhs",
                  215);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 216);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 216);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 216);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 216);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 216);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 216);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 216);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 216);
  sf_mex_assign(&c4_rhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs216), "rhs", "rhs",
                  216);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs216), "lhs", "lhs",
                  216);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 217);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 217);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 217);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 217);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 217);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 217);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 217);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 217);
  sf_mex_assign(&c4_rhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs217), "rhs", "rhs",
                  217);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs217), "lhs", "lhs",
                  217);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 218);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 218);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 218);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 218);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 218);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 218);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 218);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 218);
  sf_mex_assign(&c4_rhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs218), "rhs", "rhs",
                  218);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs218), "lhs", "lhs",
                  218);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 219);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 219);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 219);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 219);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 219);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 219);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 219);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 219);
  sf_mex_assign(&c4_rhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs219), "rhs", "rhs",
                  219);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs219), "lhs", "lhs",
                  219);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 220);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 220);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 220);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 220);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 220);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 220);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 220);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 220);
  sf_mex_assign(&c4_rhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs220), "rhs", "rhs",
                  220);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs220), "lhs", "lhs",
                  220);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 221);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 221);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 221);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 221);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 221);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 221);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 221);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 221);
  sf_mex_assign(&c4_rhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs221), "rhs", "rhs",
                  221);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs221), "lhs", "lhs",
                  221);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 222);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 222);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 222);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 222);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 222);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 222);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 222);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 222);
  sf_mex_assign(&c4_rhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs222), "rhs", "rhs",
                  222);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs222), "lhs", "lhs",
                  222);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 223);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  223);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 223);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 223);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 223);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 223);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 223);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 223);
  sf_mex_assign(&c4_rhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs223), "rhs", "rhs",
                  223);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs223), "lhs", "lhs",
                  223);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 224);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 224);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 224);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 224);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 224);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 224);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 224);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 224);
  sf_mex_assign(&c4_rhs224, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs224, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs224), "rhs", "rhs",
                  224);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs224), "lhs", "lhs",
                  224);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 225);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.xnrm2"),
                  "name", "name", 225);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 225);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "resolved", "resolved", 225);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 225);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 225);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 225);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 225);
  sf_mex_assign(&c4_rhs225, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs225, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs225), "rhs", "rhs",
                  225);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs225), "lhs", "lhs",
                  225);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 226);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 226);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 226);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 226);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 226);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 226);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 226);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 226);
  sf_mex_assign(&c4_rhs226, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs226, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs226), "rhs", "rhs",
                  226);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs226), "lhs", "lhs",
                  226);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p!below_threshold"),
                  "context", "context", 227);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 227);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 227);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 227);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 227);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 227);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 227);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 227);
  sf_mex_assign(&c4_rhs227, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs227, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs227), "rhs", "rhs",
                  227);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs227), "lhs", "lhs",
                  227);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 228);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.refblas.xnrm2"),
                  "name", "name", 228);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 228);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "resolved", "resolved", 228);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 228);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 228);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 228);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 228);
  sf_mex_assign(&c4_rhs228, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs228, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs228), "rhs", "rhs",
                  228);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs228), "lhs", "lhs",
                  228);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 229);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("realmin"), "name", "name", 229);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 229);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 229);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 229);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 229);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 229);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 229);
  sf_mex_assign(&c4_rhs229, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs229, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs229), "rhs", "rhs",
                  229);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs229), "lhs", "lhs",
                  229);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 230);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 230);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 230);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 230);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 230);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 230);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 230);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 230);
  sf_mex_assign(&c4_rhs230, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs230, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs230), "rhs", "rhs",
                  230);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs230), "lhs", "lhs",
                  230);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 231);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 231);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 231);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 231);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 231);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 231);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 231);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 231);
  sf_mex_assign(&c4_rhs231, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs231, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs231), "rhs", "rhs",
                  231);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs231), "lhs", "lhs",
                  231);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 232);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 232);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 232);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 232);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 232);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 232);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 232);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 232);
  sf_mex_assign(&c4_rhs232, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs232, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs232), "rhs", "rhs",
                  232);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs232), "lhs", "lhs",
                  232);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 233);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 233);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 233);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 233);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 233);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 233);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 233);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 233);
  sf_mex_assign(&c4_rhs233, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs233, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs233), "rhs", "rhs",
                  233);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs233), "lhs", "lhs",
                  233);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 234);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 234);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 234);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 234);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 234);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 234);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 234);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 234);
  sf_mex_assign(&c4_rhs234, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs234, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs234), "rhs", "rhs",
                  234);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs234), "lhs", "lhs",
                  234);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeev.m"),
                  "context", "context", 235);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mrdivide"), "name", "name",
                  235);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 235);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 235);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 235);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 235);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 235);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 235);
  sf_mex_assign(&c4_rhs235, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs235, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs235), "rhs", "rhs",
                  235);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs235), "lhs", "lhs",
                  235);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/eig.m"), "context",
                  "context", 236);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("diag"), "name", "name", 236);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 236);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "resolved",
                  "resolved", 236);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 236);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 236);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 236);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 236);
  sf_mex_assign(&c4_rhs236, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs236, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs236), "rhs", "rhs",
                  236);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs236), "lhs", "lhs",
                  236);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 237);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("ismatrix"), "name", "name",
                  237);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 237);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 237);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 237);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 237);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 237);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 237);
  sf_mex_assign(&c4_rhs237, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs237, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs237), "rhs", "rhs",
                  237);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs237), "lhs", "lhs",
                  237);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 238);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 238);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 238);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 238);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 238);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 238);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 238);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 238);
  sf_mex_assign(&c4_rhs238, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs238, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs238), "rhs", "rhs",
                  238);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs238), "lhs", "lhs",
                  238);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 239);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 239);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 239);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 239);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 239);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 239);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 239);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 239);
  sf_mex_assign(&c4_rhs239, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs239, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs239), "rhs", "rhs",
                  239);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs239), "lhs", "lhs",
                  239);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 240);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 240);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 240);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 240);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 240);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 240);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 240);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 240);
  sf_mex_assign(&c4_rhs240, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs240, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs240), "rhs", "rhs",
                  240);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs240), "lhs", "lhs",
                  240);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 241);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("power"), "name", "name", 241);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 241);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 241);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 241);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 241);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 241);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 241);
  sf_mex_assign(&c4_rhs241, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs241, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs241), "rhs", "rhs",
                  241);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs241), "lhs", "lhs",
                  241);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 242);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("flintmax"), "name", "name",
                  242);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 242);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/flintmax.m"), "resolved",
                  "resolved", 242);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1348199116U), "fileTimeLo",
                  "fileTimeLo", 242);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 242);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 242);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 242);
  sf_mex_assign(&c4_rhs242, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs242, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs242), "rhs", "rhs",
                  242);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs242), "lhs", "lhs",
                  242);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/flintmax.m"), "context",
                  "context", 243);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 243);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 243);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 243);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 243);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 243);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 243);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 243);
  sf_mex_assign(&c4_rhs243, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs243, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs243), "rhs", "rhs",
                  243);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs243), "lhs", "lhs",
                  243);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 244);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mod"), "name", "name", 244);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 244);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m"), "resolved",
                  "resolved", 244);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 244);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 244);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 244);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 244);
  sf_mex_assign(&c4_rhs244, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs244, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs244), "rhs", "rhs",
                  244);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs244), "lhs", "lhs",
                  244);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 245);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 245);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 245);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 245);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 245);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 245);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 245);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 245);
  sf_mex_assign(&c4_rhs245, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs245, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs245), "rhs", "rhs",
                  245);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs245), "lhs", "lhs",
                  245);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 246);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 246);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 246);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 246);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 246);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 246);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 246);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 246);
  sf_mex_assign(&c4_rhs246, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs246, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs246), "rhs", "rhs",
                  246);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs246), "lhs", "lhs",
                  246);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 247);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_round"), "name",
                  "name", 247);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 247);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_round.m"),
                  "resolved", "resolved", 247);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658438U), "fileTimeLo",
                  "fileTimeLo", 247);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 247);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 247);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 247);
  sf_mex_assign(&c4_rhs247, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs247, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs247), "rhs", "rhs",
                  247);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs247), "lhs", "lhs",
                  247);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 248);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 248);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 248);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 248);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 248);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 248);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 248);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 248);
  sf_mex_assign(&c4_rhs248, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs248, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs248), "rhs", "rhs",
                  248);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs248), "lhs", "lhs",
                  248);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/mod.m!floatmod"), "context",
                  "context", 249);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eps"), "name", "name", 249);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 249);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 249);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 249);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 249);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 249);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 249);
  sf_mex_assign(&c4_rhs249, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs249, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs249), "rhs", "rhs",
                  249);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs249), "lhs", "lhs",
                  249);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 250);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("log"), "name", "name", 250);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 250);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m"), "resolved",
                  "resolved", 250);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1343837580U), "fileTimeLo",
                  "fileTimeLo", 250);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 250);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 250);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 250);
  sf_mex_assign(&c4_rhs250, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs250, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs250), "rhs", "rhs",
                  250);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs250), "lhs", "lhs",
                  250);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log.m"), "context",
                  "context", 251);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_log"), "name",
                  "name", 251);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 251);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "resolved", "resolved", 251);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825928U), "fileTimeLo",
                  "fileTimeLo", 251);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 251);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 251);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 251);
  sf_mex_assign(&c4_rhs251, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs251, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs251), "rhs", "rhs",
                  251);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs251), "lhs", "lhs",
                  251);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 252);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("realmax"), "name", "name", 252);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 252);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmax.m"), "resolved",
                  "resolved", 252);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 252);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 252);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 252);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 252);
  sf_mex_assign(&c4_rhs252, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs252, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs252), "rhs", "rhs",
                  252);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs252), "lhs", "lhs",
                  252);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 253);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mrdivide"), "name", "name",
                  253);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 253);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 253);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 253);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 253);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 253);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 253);
  sf_mex_assign(&c4_rhs253, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs253, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs253), "rhs", "rhs",
                  253);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs253), "lhs", "lhs",
                  253);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 254);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("isnan"), "name", "name", 254);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 254);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 254);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 254);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 254);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 254);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 254);
  sf_mex_assign(&c4_rhs254, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs254, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs254), "rhs", "rhs",
                  254);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs254), "lhs", "lhs",
                  254);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 255);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_hypot"), "name",
                  "name", 255);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 255);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_hypot.m"),
                  "resolved", "resolved", 255);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825926U), "fileTimeLo",
                  "fileTimeLo", 255);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 255);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 255);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 255);
  sf_mex_assign(&c4_rhs255, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs255, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs255), "rhs", "rhs",
                  255);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs255), "lhs", "lhs",
                  255);
  sf_mex_destroy(&c4_rhs192);
  sf_mex_destroy(&c4_lhs192);
  sf_mex_destroy(&c4_rhs193);
  sf_mex_destroy(&c4_lhs193);
  sf_mex_destroy(&c4_rhs194);
  sf_mex_destroy(&c4_lhs194);
  sf_mex_destroy(&c4_rhs195);
  sf_mex_destroy(&c4_lhs195);
  sf_mex_destroy(&c4_rhs196);
  sf_mex_destroy(&c4_lhs196);
  sf_mex_destroy(&c4_rhs197);
  sf_mex_destroy(&c4_lhs197);
  sf_mex_destroy(&c4_rhs198);
  sf_mex_destroy(&c4_lhs198);
  sf_mex_destroy(&c4_rhs199);
  sf_mex_destroy(&c4_lhs199);
  sf_mex_destroy(&c4_rhs200);
  sf_mex_destroy(&c4_lhs200);
  sf_mex_destroy(&c4_rhs201);
  sf_mex_destroy(&c4_lhs201);
  sf_mex_destroy(&c4_rhs202);
  sf_mex_destroy(&c4_lhs202);
  sf_mex_destroy(&c4_rhs203);
  sf_mex_destroy(&c4_lhs203);
  sf_mex_destroy(&c4_rhs204);
  sf_mex_destroy(&c4_lhs204);
  sf_mex_destroy(&c4_rhs205);
  sf_mex_destroy(&c4_lhs205);
  sf_mex_destroy(&c4_rhs206);
  sf_mex_destroy(&c4_lhs206);
  sf_mex_destroy(&c4_rhs207);
  sf_mex_destroy(&c4_lhs207);
  sf_mex_destroy(&c4_rhs208);
  sf_mex_destroy(&c4_lhs208);
  sf_mex_destroy(&c4_rhs209);
  sf_mex_destroy(&c4_lhs209);
  sf_mex_destroy(&c4_rhs210);
  sf_mex_destroy(&c4_lhs210);
  sf_mex_destroy(&c4_rhs211);
  sf_mex_destroy(&c4_lhs211);
  sf_mex_destroy(&c4_rhs212);
  sf_mex_destroy(&c4_lhs212);
  sf_mex_destroy(&c4_rhs213);
  sf_mex_destroy(&c4_lhs213);
  sf_mex_destroy(&c4_rhs214);
  sf_mex_destroy(&c4_lhs214);
  sf_mex_destroy(&c4_rhs215);
  sf_mex_destroy(&c4_lhs215);
  sf_mex_destroy(&c4_rhs216);
  sf_mex_destroy(&c4_lhs216);
  sf_mex_destroy(&c4_rhs217);
  sf_mex_destroy(&c4_lhs217);
  sf_mex_destroy(&c4_rhs218);
  sf_mex_destroy(&c4_lhs218);
  sf_mex_destroy(&c4_rhs219);
  sf_mex_destroy(&c4_lhs219);
  sf_mex_destroy(&c4_rhs220);
  sf_mex_destroy(&c4_lhs220);
  sf_mex_destroy(&c4_rhs221);
  sf_mex_destroy(&c4_lhs221);
  sf_mex_destroy(&c4_rhs222);
  sf_mex_destroy(&c4_lhs222);
  sf_mex_destroy(&c4_rhs223);
  sf_mex_destroy(&c4_lhs223);
  sf_mex_destroy(&c4_rhs224);
  sf_mex_destroy(&c4_lhs224);
  sf_mex_destroy(&c4_rhs225);
  sf_mex_destroy(&c4_lhs225);
  sf_mex_destroy(&c4_rhs226);
  sf_mex_destroy(&c4_lhs226);
  sf_mex_destroy(&c4_rhs227);
  sf_mex_destroy(&c4_lhs227);
  sf_mex_destroy(&c4_rhs228);
  sf_mex_destroy(&c4_lhs228);
  sf_mex_destroy(&c4_rhs229);
  sf_mex_destroy(&c4_lhs229);
  sf_mex_destroy(&c4_rhs230);
  sf_mex_destroy(&c4_lhs230);
  sf_mex_destroy(&c4_rhs231);
  sf_mex_destroy(&c4_lhs231);
  sf_mex_destroy(&c4_rhs232);
  sf_mex_destroy(&c4_lhs232);
  sf_mex_destroy(&c4_rhs233);
  sf_mex_destroy(&c4_lhs233);
  sf_mex_destroy(&c4_rhs234);
  sf_mex_destroy(&c4_lhs234);
  sf_mex_destroy(&c4_rhs235);
  sf_mex_destroy(&c4_lhs235);
  sf_mex_destroy(&c4_rhs236);
  sf_mex_destroy(&c4_lhs236);
  sf_mex_destroy(&c4_rhs237);
  sf_mex_destroy(&c4_lhs237);
  sf_mex_destroy(&c4_rhs238);
  sf_mex_destroy(&c4_lhs238);
  sf_mex_destroy(&c4_rhs239);
  sf_mex_destroy(&c4_lhs239);
  sf_mex_destroy(&c4_rhs240);
  sf_mex_destroy(&c4_lhs240);
  sf_mex_destroy(&c4_rhs241);
  sf_mex_destroy(&c4_lhs241);
  sf_mex_destroy(&c4_rhs242);
  sf_mex_destroy(&c4_lhs242);
  sf_mex_destroy(&c4_rhs243);
  sf_mex_destroy(&c4_lhs243);
  sf_mex_destroy(&c4_rhs244);
  sf_mex_destroy(&c4_lhs244);
  sf_mex_destroy(&c4_rhs245);
  sf_mex_destroy(&c4_lhs245);
  sf_mex_destroy(&c4_rhs246);
  sf_mex_destroy(&c4_lhs246);
  sf_mex_destroy(&c4_rhs247);
  sf_mex_destroy(&c4_lhs247);
  sf_mex_destroy(&c4_rhs248);
  sf_mex_destroy(&c4_lhs248);
  sf_mex_destroy(&c4_rhs249);
  sf_mex_destroy(&c4_lhs249);
  sf_mex_destroy(&c4_rhs250);
  sf_mex_destroy(&c4_lhs250);
  sf_mex_destroy(&c4_rhs251);
  sf_mex_destroy(&c4_lhs251);
  sf_mex_destroy(&c4_rhs252);
  sf_mex_destroy(&c4_lhs252);
  sf_mex_destroy(&c4_rhs253);
  sf_mex_destroy(&c4_lhs253);
  sf_mex_destroy(&c4_rhs254);
  sf_mex_destroy(&c4_lhs254);
  sf_mex_destroy(&c4_rhs255);
  sf_mex_destroy(&c4_lhs255);
}

static void c4_e_info_helper(const mxArray **c4_info)
{
  const mxArray *c4_rhs256 = NULL;
  const mxArray *c4_lhs256 = NULL;
  const mxArray *c4_rhs257 = NULL;
  const mxArray *c4_lhs257 = NULL;
  const mxArray *c4_rhs258 = NULL;
  const mxArray *c4_lhs258 = NULL;
  const mxArray *c4_rhs259 = NULL;
  const mxArray *c4_lhs259 = NULL;
  const mxArray *c4_rhs260 = NULL;
  const mxArray *c4_lhs260 = NULL;
  const mxArray *c4_rhs261 = NULL;
  const mxArray *c4_lhs261 = NULL;
  const mxArray *c4_rhs262 = NULL;
  const mxArray *c4_lhs262 = NULL;
  const mxArray *c4_rhs263 = NULL;
  const mxArray *c4_lhs263 = NULL;
  const mxArray *c4_rhs264 = NULL;
  const mxArray *c4_lhs264 = NULL;
  const mxArray *c4_rhs265 = NULL;
  const mxArray *c4_lhs265 = NULL;
  const mxArray *c4_rhs266 = NULL;
  const mxArray *c4_lhs266 = NULL;
  const mxArray *c4_rhs267 = NULL;
  const mxArray *c4_lhs267 = NULL;
  const mxArray *c4_rhs268 = NULL;
  const mxArray *c4_lhs268 = NULL;
  const mxArray *c4_rhs269 = NULL;
  const mxArray *c4_lhs269 = NULL;
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 256);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_atan2"), "name",
                  "name", 256);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 256);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m"),
                  "resolved", "resolved", 256);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825920U), "fileTimeLo",
                  "fileTimeLo", 256);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 256);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 256);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 256);
  sf_mex_assign(&c4_rhs256, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs256, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs256), "rhs", "rhs",
                  256);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs256), "lhs", "lhs",
                  256);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_log.m"),
                  "context", "context", 257);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 257);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 257);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 257);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 257);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 257);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 257);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 257);
  sf_mex_assign(&c4_rhs257, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs257, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs257), "rhs", "rhs",
                  257);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs257), "lhs", "lhs",
                  257);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 258);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 258);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 258);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 258);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 258);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 258);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 258);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 258);
  sf_mex_assign(&c4_rhs258, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs258, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs258), "rhs", "rhs",
                  258);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs258), "lhs", "lhs",
                  258);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m!matrix_to_scalar_power"),
                  "context", "context", 259);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mrdivide"), "name", "name",
                  259);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 259);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 259);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 259);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 259);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 259);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 259);
  sf_mex_assign(&c4_rhs259, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs259, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs259), "rhs", "rhs",
                  259);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs259), "lhs", "lhs",
                  259);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 260);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("ismatrix"), "name", "name",
                  260);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 260);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 260);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 260);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 260);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 260);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 260);
  sf_mex_assign(&c4_rhs260, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs260, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs260), "rhs", "rhs",
                  260);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs260), "lhs", "lhs",
                  260);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 261);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_lusolve"), "name", "name",
                  261);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 261);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "resolved",
                  "resolved", 261);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1370017086U), "fileTimeLo",
                  "fileTimeLo", 261);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 261);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 261);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 261);
  sf_mex_assign(&c4_rhs261, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs261, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs261), "rhs", "rhs",
                  261);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs261), "lhs", "lhs",
                  261);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 262);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_xcabs1"), "name", "name",
                  262);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 262);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "resolved", "resolved", 262);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 262);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 262);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 262);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 262);
  sf_mex_assign(&c4_rhs262, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs262, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs262), "rhs", "rhs",
                  262);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs262), "lhs", "lhs",
                  262);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "context", "context", 263);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("coder.internal.refblas.xcabs1"),
                  "name", "name", 263);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 263);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "resolved", "resolved", 263);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 263);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 263);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 263);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 263);
  sf_mex_assign(&c4_rhs263, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs263, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs263), "rhs", "rhs",
                  263);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs263), "lhs", "lhs",
                  263);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "context", "context", 264);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("abs"), "name", "name", 264);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 264);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 264);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 264);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 264);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 264);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 264);
  sf_mex_assign(&c4_rhs264, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs264, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs264), "rhs", "rhs",
                  264);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs264), "lhs", "lhs",
                  264);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 265);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("rdivide"), "name", "name", 265);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 265);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 265);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 265);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 265);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 265);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 265);
  sf_mex_assign(&c4_rhs265, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs265, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs265), "rhs", "rhs",
                  265);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs265), "lhs", "lhs",
                  265);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular"),
                  "context", "context", 266);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_warning"), "name", "name",
                  266);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 266);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 266);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1286826002U), "fileTimeLo",
                  "fileTimeLo", 266);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 266);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 266);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 266);
  sf_mex_assign(&c4_rhs266, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs266, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs266), "rhs", "rhs",
                  266);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs266), "lhs", "lhs",
                  266);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 267);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 267);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 267);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 267);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 267);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 267);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 267);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 267);
  sf_mex_assign(&c4_rhs267, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs267, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs267), "rhs", "rhs",
                  267);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs267), "lhs", "lhs",
                  267);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 268);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 268);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 268);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 268);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 268);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 268);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 268);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 268);
  sf_mex_assign(&c4_rhs268, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs268, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs268), "rhs", "rhs",
                  268);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs268), "lhs", "lhs",
                  268);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 269);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 269);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 269);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 269);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 269);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 269);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 269);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 269);
  sf_mex_assign(&c4_rhs269, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c4_lhs269, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs269), "rhs", "rhs",
                  269);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs269), "lhs", "lhs",
                  269);
  sf_mex_destroy(&c4_rhs256);
  sf_mex_destroy(&c4_lhs256);
  sf_mex_destroy(&c4_rhs257);
  sf_mex_destroy(&c4_lhs257);
  sf_mex_destroy(&c4_rhs258);
  sf_mex_destroy(&c4_lhs258);
  sf_mex_destroy(&c4_rhs259);
  sf_mex_destroy(&c4_lhs259);
  sf_mex_destroy(&c4_rhs260);
  sf_mex_destroy(&c4_lhs260);
  sf_mex_destroy(&c4_rhs261);
  sf_mex_destroy(&c4_lhs261);
  sf_mex_destroy(&c4_rhs262);
  sf_mex_destroy(&c4_lhs262);
  sf_mex_destroy(&c4_rhs263);
  sf_mex_destroy(&c4_lhs263);
  sf_mex_destroy(&c4_rhs264);
  sf_mex_destroy(&c4_lhs264);
  sf_mex_destroy(&c4_rhs265);
  sf_mex_destroy(&c4_lhs265);
  sf_mex_destroy(&c4_rhs266);
  sf_mex_destroy(&c4_lhs266);
  sf_mex_destroy(&c4_rhs267);
  sf_mex_destroy(&c4_lhs267);
  sf_mex_destroy(&c4_rhs268);
  sf_mex_destroy(&c4_lhs268);
  sf_mex_destroy(&c4_rhs269);
  sf_mex_destroy(&c4_lhs269);
}

static void c4_eig(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_A[9],
                   creal_T c4_V[3])
{
  int32_T c4_i42;
  static creal_T c4_dc0 = { 0.0, 0.0 };

  creal_T c4_b_A[9];
  real_T c4_info;
  int32_T c4_i43;
  creal_T c4_c_A[9];
  real_T c4_anrm;
  int32_T c4_i44;
  creal_T c4_alpha1[3];
  int32_T c4_i45;
  creal_T c4_beta1[3];
  boolean_T c4_ilascl;
  real_T c4_anrmto;
  int32_T c4_rscale[3];
  int32_T c4_ihi;
  int32_T c4_ilo;
  int32_T c4_b_ilo;
  int32_T c4_b_ihi;
  int32_T c4_c_ilo;
  int32_T c4_c_ihi;
  int32_T c4_a;
  int32_T c4_b_a;
  int32_T c4_c;
  int32_T c4_c_a;
  int32_T c4_d_a;
  int32_T c4_ihim1;
  int32_T c4_jcol;
  int32_T c4_e_a;
  int32_T c4_f_a;
  int32_T c4_jcolp1;
  int32_T c4_jrow;
  int32_T c4_g_a;
  int32_T c4_h_a;
  int32_T c4_jrowm1;
  creal_T c4_d_A;
  creal_T c4_e_A;
  creal_T c4_b;
  creal_T c4_s;
  real_T c4_b_c;
  real_T c4_c_c;
  real_T c4_d_c;
  int32_T c4_xrow;
  int32_T c4_yrow;
  int32_T c4_jlo;
  int32_T c4_jhi;
  int32_T c4_b_jlo;
  int32_T c4_b_jhi;
  int32_T c4_i_a;
  int32_T c4_b_b;
  int32_T c4_j_a;
  int32_T c4_c_b;
  boolean_T c4_overflow;
  int32_T c4_j;
  int32_T c4_b_j;
  real_T c4_k_a;
  creal_T c4_y;
  creal_T c4_b_s;
  creal_T c4_stemp;
  real_T c4_l_a;
  creal_T c4_d_b;
  creal_T c4_e_b;
  creal_T c4_f_b;
  creal_T c4_g_b;
  real_T c4_e_c;
  int32_T c4_xcol;
  int32_T c4_ycol;
  int32_T c4_d_ilo;
  int32_T c4_d_ihi;
  int32_T c4_e_ilo;
  int32_T c4_e_ihi;
  int32_T c4_m_a;
  int32_T c4_h_b;
  int32_T c4_n_a;
  int32_T c4_i_b;
  boolean_T c4_b_overflow;
  int32_T c4_i;
  int32_T c4_b_i;
  real_T c4_o_a;
  creal_T c4_c_s;
  real_T c4_p_a;
  creal_T c4_j_b;
  creal_T c4_k_b;
  creal_T c4_l_b;
  creal_T c4_m_b;
  int32_T c4_i46;
  creal_T c4_f_A[9];
  real_T c4_b_info;
  real_T c4_c_info;
  real_T c4_d_info;
  real_T c4_e_info;
  int32_T c4_i47;
  creal_T c4_b_alpha1[3];
  int32_T c4_i48;
  creal_T c4_b_beta1[3];
  boolean_T guard1 = false;
  for (c4_i42 = 0; c4_i42 < 9; c4_i42++) {
    c4_b_A[c4_i42].re = c4_A[c4_i42] + c4_dc0.re;
    c4_b_A[c4_i42].im = c4_dc0.im;
  }

  c4_info = 0.0;
  c4_realmin(chartInstance);
  c4_eps(chartInstance);
  for (c4_i43 = 0; c4_i43 < 9; c4_i43++) {
    c4_c_A[c4_i43] = c4_b_A[c4_i43];
  }

  c4_anrm = c4_eml_matlab_zlangeM(chartInstance, c4_c_A);
  if (!c4_isfinite(chartInstance, c4_anrm)) {
    for (c4_i44 = 0; c4_i44 < 3; c4_i44++) {
      c4_alpha1[c4_i44].re = rtNaN;
      c4_alpha1[c4_i44].im = 0.0;
    }

    for (c4_i45 = 0; c4_i45 < 3; c4_i45++) {
      c4_beta1[c4_i45].re = rtNaN;
      c4_beta1[c4_i45].im = 0.0;
    }
  } else {
    c4_ilascl = false;
    c4_anrmto = c4_anrm;
    guard1 = false;
    if (c4_anrm > 0.0) {
      if (c4_anrm < 6.7178761075670888E-139) {
        c4_anrmto = 6.7178761075670888E-139;
        c4_ilascl = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      if (c4_anrm > 1.4885657073574029E+138) {
        c4_anrmto = 1.4885657073574029E+138;
        c4_ilascl = true;
      }
    }

    if (c4_ilascl) {
      c4_c_eml_matlab_zlascl(chartInstance, c4_anrm, c4_anrmto, c4_b_A);
    }

    c4_b_eml_matlab_zggbal(chartInstance, c4_b_A, &c4_ilo, &c4_ihi, c4_rscale);
    c4_b_ilo = c4_ilo;
    c4_b_ihi = c4_ihi;
    c4_c_ilo = c4_b_ilo;
    c4_c_ihi = c4_b_ihi;
    c4_a = c4_c_ilo;
    c4_b_a = c4_a + 2;
    c4_c = c4_b_a;
    if (c4_c_ihi < c4_c) {
    } else {
      c4_c_a = c4_c_ihi;
      c4_d_a = c4_c_a - 1;
      c4_ihim1 = c4_d_a;
      c4_jcol = c4_c_ilo;
      while (c4_jcol < c4_ihim1) {
        c4_e_a = c4_jcol;
        c4_f_a = c4_e_a + 1;
        c4_jcolp1 = c4_f_a;
        c4_jrow = c4_c_ihi;
        while (c4_jrow > c4_jcolp1) {
          c4_g_a = c4_jrow;
          c4_h_a = c4_g_a - 1;
          c4_jrowm1 = c4_h_a;
          c4_d_A.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_jrowm1), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].re;
          c4_d_A.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_jrowm1), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].im;
          c4_e_A.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_jrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].re;
          c4_e_A.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_jrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].im;
          c4_eml_matlab_zlartg(chartInstance, c4_d_A, c4_e_A, &c4_b_c, &c4_s,
                               &c4_b);
          c4_c_c = c4_b_c;
          c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_jrowm1), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].re = c4_b.re;
          c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_jrowm1), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].im = c4_b.im;
          c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_jrow), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].re = c4_dc0.re;
          c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_jrow), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].im = c4_dc0.im;
          c4_d_c = c4_c_c;
          c4_xrow = c4_jrowm1;
          c4_yrow = c4_jrow;
          c4_jlo = c4_jcolp1;
          c4_jhi = c4_c_ihi;
          c4_b_jlo = c4_jlo;
          c4_b_jhi = c4_jhi;
          c4_i_a = c4_b_jlo;
          c4_b_b = c4_b_jhi;
          c4_j_a = c4_i_a;
          c4_c_b = c4_b_b;
          if (c4_j_a > c4_c_b) {
            c4_overflow = false;
          } else {
            c4_eml_switch_helper(chartInstance);
            c4_overflow = (c4_c_b > 2147483646);
          }

          if (c4_overflow) {
            c4_check_forloop_overflow_error(chartInstance, c4_overflow);
          }

          for (c4_j = c4_b_jlo; c4_j <= c4_b_jhi; c4_j++) {
            c4_b_j = c4_j;
            c4_k_a = c4_d_c;
            c4_b.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c4_b.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c4_y.re = c4_k_a * c4_b.re;
            c4_y.im = c4_k_a * c4_b.im;
            c4_b_s.re = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re - c4_s.im * c4_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c4_b_s.im = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im + c4_s.im * c4_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c4_stemp.re = c4_y.re + c4_b_s.re;
            c4_stemp.im = c4_y.im + c4_b_s.im;
            c4_l_a = c4_d_c;
            c4_b.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c4_b.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c4_y.re = c4_l_a * c4_b.re;
            c4_y.im = c4_l_a * c4_b.im;
            c4_b.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
            c4_b.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
            c4_d_b = c4_b;
            c4_e_b = c4_b;
            c4_f_b = c4_b;
            c4_g_b = c4_b;
            c4_b.re = c4_s.re * c4_d_b.re + c4_s.im * c4_e_b.im;
            c4_b.im = c4_s.re * c4_f_b.im - c4_s.im * c4_g_b.re;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re = c4_y.re
              - c4_b.re;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im = c4_y.im
              - c4_b.im;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re =
              c4_stemp.re;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im =
              c4_stemp.im;
          }

          c4_s.re = -c4_s.re;
          c4_s.im = -c4_s.im;
          c4_e_c = c4_c_c;
          c4_xcol = c4_jrow;
          c4_ycol = c4_jrowm1;
          c4_d_ilo = c4_c_ilo;
          c4_d_ihi = c4_c_ihi;
          c4_e_ilo = c4_d_ilo;
          c4_e_ihi = c4_d_ihi;
          c4_m_a = c4_e_ilo;
          c4_h_b = c4_e_ihi;
          c4_n_a = c4_m_a;
          c4_i_b = c4_h_b;
          if (c4_n_a > c4_i_b) {
            c4_b_overflow = false;
          } else {
            c4_eml_switch_helper(chartInstance);
            c4_b_overflow = (c4_i_b > 2147483646);
          }

          if (c4_b_overflow) {
            c4_check_forloop_overflow_error(chartInstance, c4_b_overflow);
          }

          for (c4_i = c4_e_ilo; c4_i <= c4_e_ihi; c4_i++) {
            c4_b_i = c4_i;
            c4_o_a = c4_e_c;
            c4_b.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].re;
            c4_b.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].im;
            c4_y.re = c4_o_a * c4_b.re;
            c4_y.im = c4_o_a * c4_b.im;
            c4_c_s.re = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re - c4_s.im * c4_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im;
            c4_c_s.im = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im + c4_s.im * c4_b_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re;
            c4_stemp.re = c4_y.re + c4_c_s.re;
            c4_stemp.im = c4_y.im + c4_c_s.im;
            c4_p_a = c4_e_c;
            c4_b.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re;
            c4_b.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im;
            c4_y.re = c4_p_a * c4_b.re;
            c4_y.im = c4_p_a * c4_b.im;
            c4_b.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].re;
            c4_b.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].im;
            c4_j_b = c4_b;
            c4_k_b = c4_b;
            c4_l_b = c4_b;
            c4_m_b = c4_b;
            c4_b.re = c4_s.re * c4_j_b.re + c4_s.im * c4_k_b.im;
            c4_b.im = c4_s.re * c4_l_b.im - c4_s.im * c4_m_b.re;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re = c4_y.re
              - c4_b.re;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im = c4_y.im
              - c4_b.im;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].re =
              c4_stemp.re;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].im =
              c4_stemp.im;
          }

          c4_jrow = c4_jrowm1;
        }

        c4_jcol = c4_jcolp1;
      }
    }

    for (c4_i46 = 0; c4_i46 < 9; c4_i46++) {
      c4_f_A[c4_i46] = c4_b_A[c4_i46];
    }

    c4_eml_matlab_zhgeqz(chartInstance, c4_f_A, c4_b_ilo, c4_b_ihi, &c4_b_info,
                         c4_alpha1, c4_beta1);
    c4_info = c4_b_info;
    if (c4_info != 0.0) {
    } else {
      if (c4_ilascl) {
        c4_d_eml_matlab_zlascl(chartInstance, c4_anrmto, c4_anrm, c4_alpha1);
      }
    }
  }

  c4_c_info = c4_info;
  c4_d_info = c4_c_info;
  c4_e_info = c4_d_info;
  for (c4_i47 = 0; c4_i47 < 3; c4_i47++) {
    c4_b_alpha1[c4_i47] = c4_alpha1[c4_i47];
  }

  for (c4_i48 = 0; c4_i48 < 3; c4_i48++) {
    c4_b_beta1[c4_i48] = c4_beta1[c4_i48];
  }

  c4_b_eml_div(chartInstance, c4_b_alpha1, c4_b_beta1, c4_V);
  if (c4_e_info < 0.0) {
    c4_eml_warning(chartInstance);
  } else {
    if (c4_e_info > 0.0) {
      c4_b_eml_warning(chartInstance);
    }
  }
}

static void c4_realmin(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c4_eml_error(SFc4_Model_01InstanceStruct *chartInstance)
{
  int32_T c4_i49;
  static char_T c4_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c4_u[30];
  const mxArray *c4_y = NULL;
  int32_T c4_i50;
  static char_T c4_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c4_b_u[4];
  const mxArray *c4_b_y = NULL;
  (void)chartInstance;
  for (c4_i49 = 0; c4_i49 < 30; c4_i49++) {
    c4_u[c4_i49] = c4_cv0[c4_i49];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  for (c4_i50 = 0; c4_i50 < 4; c4_i50++) {
    c4_b_u[c4_i50] = c4_cv1[c4_i50];
  }

  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", c4_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c4_y, 14, c4_b_y));
}

static void c4_eps(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c4_eml_matlab_zlangeM(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_x[9])
{
  real_T c4_y;
  int32_T c4_k;
  real_T c4_b_k;
  creal_T c4_b_x;
  real_T c4_x1;
  real_T c4_x2;
  real_T c4_a;
  real_T c4_b;
  real_T c4_absxk;
  real_T c4_c_x;
  boolean_T c4_b_b;
  boolean_T exitg1;
  (void)chartInstance;
  c4_y = 0.0;
  c4_k = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c4_k < 9)) {
    c4_b_k = 1.0 + (real_T)c4_k;
    c4_b_x.re = c4_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", c4_b_k), 1, 9, 1, 0) - 1].re;
    c4_b_x.im = c4_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", c4_b_k), 1, 9, 1, 0) - 1].im;
    c4_x1 = c4_b_x.re;
    c4_x2 = c4_b_x.im;
    c4_a = c4_x1;
    c4_b = c4_x2;
    c4_absxk = muDoubleScalarHypot(c4_a, c4_b);
    c4_c_x = c4_absxk;
    c4_b_b = muDoubleScalarIsNaN(c4_c_x);
    if (c4_b_b) {
      c4_y = rtNaN;
      exitg1 = true;
    } else {
      if (c4_absxk > c4_y) {
        c4_y = c4_absxk;
      }

      c4_k++;
    }
  }

  return c4_y;
}

static real_T c4_abs(SFc4_Model_01InstanceStruct *chartInstance, creal_T c4_x)
{
  real_T c4_x1;
  real_T c4_x2;
  real_T c4_a;
  real_T c4_b;
  (void)chartInstance;
  c4_x1 = c4_x.re;
  c4_x2 = c4_x.im;
  c4_a = c4_x1;
  c4_b = c4_x2;
  return muDoubleScalarHypot(c4_a, c4_b);
}

static boolean_T c4_isfinite(SFc4_Model_01InstanceStruct *chartInstance, real_T
  c4_x)
{
  real_T c4_b_x;
  boolean_T c4_b_b;
  boolean_T c4_b0;
  real_T c4_c_x;
  boolean_T c4_c_b;
  boolean_T c4_b1;
  (void)chartInstance;
  c4_b_x = c4_x;
  c4_b_b = muDoubleScalarIsInf(c4_b_x);
  c4_b0 = !c4_b_b;
  c4_c_x = c4_x;
  c4_c_b = muDoubleScalarIsNaN(c4_c_x);
  c4_b1 = !c4_c_b;
  return c4_b0 && c4_b1;
}

static void c4_eml_matlab_zlascl(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_cfrom, real_T c4_cto, creal_T c4_A[9], creal_T c4_b_A[9])
{
  int32_T c4_i51;
  for (c4_i51 = 0; c4_i51 < 9; c4_i51++) {
    c4_b_A[c4_i51] = c4_A[c4_i51];
  }

  c4_c_eml_matlab_zlascl(chartInstance, c4_cfrom, c4_cto, c4_b_A);
}

static real_T c4_b_abs(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_x)
{
  real_T c4_b_x;
  (void)chartInstance;
  c4_b_x = c4_x;
  return muDoubleScalarAbs(c4_b_x);
}

static void c4_eml_matlab_zggbal(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], creal_T c4_b_A[9], int32_T *c4_ilo, int32_T *c4_ihi, int32_T
  c4_rscale[3])
{
  int32_T c4_i52;
  for (c4_i52 = 0; c4_i52 < 9; c4_i52++) {
    c4_b_A[c4_i52] = c4_A[c4_i52];
  }

  c4_b_eml_matlab_zggbal(chartInstance, c4_b_A, c4_ilo, c4_ihi, c4_rscale);
}

static void c4_eml_switch_helper(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c4_check_forloop_overflow_error(SFc4_Model_01InstanceStruct
  *chartInstance, boolean_T c4_overflow)
{
  int32_T c4_i53;
  static char_T c4_cv2[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c4_u[34];
  const mxArray *c4_y = NULL;
  int32_T c4_i54;
  static char_T c4_cv3[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c4_b_u[23];
  const mxArray *c4_b_y = NULL;
  (void)chartInstance;
  if (!c4_overflow) {
  } else {
    for (c4_i53 = 0; c4_i53 < 34; c4_i53++) {
      c4_u[c4_i53] = c4_cv2[c4_i53];
    }

    c4_y = NULL;
    sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  false);
    for (c4_i54 = 0; c4_i54 < 23; c4_i54++) {
      c4_b_u[c4_i54] = c4_cv3[c4_i54];
    }

    c4_b_y = NULL;
    sf_mex_assign(&c4_b_y, sf_mex_create("y", c4_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 2U, 14, c4_y, 14, c4_b_y));
  }
}

static void c4_eml_matlab_zlartg(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_f, creal_T c4_g, real_T *c4_cs, creal_T *c4_sn, creal_T *c4_r)
{
  real_T c4_x;
  real_T c4_b_x;
  real_T c4_y;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_b_y;
  real_T c4_e_x;
  real_T c4_c_y;
  real_T c4_d_y;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_e_y;
  real_T c4_h_x;
  real_T c4_i_x;
  real_T c4_f_y;
  real_T c4_j_x;
  real_T c4_g_y;
  real_T c4_h_y;
  real_T c4_k_x;
  real_T c4_i_y;
  real_T c4_scale;
  creal_T c4_fs;
  creal_T c4_gs;
  int32_T c4_count;
  real_T c4_rescaledir;
  int32_T c4_a;
  int32_T c4_b_a;
  static creal_T c4_dc1 = { 0.0, 0.0 };

  boolean_T c4_b_g;
  int32_T c4_c_a;
  int32_T c4_d_a;
  real_T c4_f2;
  real_T c4_g2;
  real_T c4_l_x;
  real_T c4_m_x;
  boolean_T c4_b_f;
  real_T c4_x1;
  real_T c4_x2;
  real_T c4_e_a;
  real_T c4_b;
  real_T c4_j_y;
  real_T c4_b_x1;
  real_T c4_b_x2;
  real_T c4_f_a;
  real_T c4_b_b;
  real_T c4_d;
  real_T c4_c_x1;
  real_T c4_c_x2;
  real_T c4_g_a;
  real_T c4_c_b;
  real_T c4_f2s;
  real_T c4_n_x;
  real_T c4_g2s;
  real_T c4_o_x;
  real_T c4_p_x;
  real_T c4_k_y;
  real_T c4_q_x;
  real_T c4_r_x;
  real_T c4_l_y;
  real_T c4_s_x;
  real_T c4_m_y;
  real_T c4_n_y;
  real_T c4_d_x1;
  real_T c4_d_x2;
  real_T c4_h_a;
  real_T c4_d_b;
  real_T c4_dr;
  real_T c4_di;
  real_T c4_e_x1;
  real_T c4_e_x2;
  real_T c4_i_a;
  real_T c4_e_b;
  creal_T c4_b_gs;
  real_T c4_j_a;
  creal_T c4_b_sn;
  real_T c4_t_x;
  creal_T c4_c_gs;
  creal_T c4_c_sn;
  int32_T c4_b_count;
  int32_T c4_f_b;
  int32_T c4_g_b;
  boolean_T c4_overflow;
  int32_T c4_i;
  int32_T c4_c_count;
  int32_T c4_h_b;
  int32_T c4_i_b;
  boolean_T c4_b_overflow;
  int32_T c4_b_i;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  c4_realmin(chartInstance);
  c4_eps(chartInstance);
  c4_eps(chartInstance);
  c4_x = c4_f.re;
  c4_b_x = c4_x;
  c4_y = muDoubleScalarAbs(c4_b_x);
  c4_c_x = c4_f.im;
  c4_d_x = c4_c_x;
  c4_b_y = muDoubleScalarAbs(c4_d_x);
  c4_e_x = c4_y;
  c4_c_y = c4_b_y;
  c4_d_y = c4_e_x;
  if (c4_c_y > c4_d_y) {
    c4_d_y = c4_c_y;
  }

  c4_f_x = c4_g.re;
  c4_g_x = c4_f_x;
  c4_e_y = muDoubleScalarAbs(c4_g_x);
  c4_h_x = c4_g.im;
  c4_i_x = c4_h_x;
  c4_f_y = muDoubleScalarAbs(c4_i_x);
  c4_j_x = c4_e_y;
  c4_g_y = c4_f_y;
  c4_h_y = c4_j_x;
  if (c4_g_y > c4_h_y) {
    c4_h_y = c4_g_y;
  }

  c4_k_x = c4_d_y;
  c4_i_y = c4_h_y;
  c4_scale = c4_k_x;
  if (c4_i_y > c4_scale) {
    c4_scale = c4_i_y;
  }

  c4_fs = c4_f;
  c4_gs = c4_g;
  c4_count = 0;
  c4_rescaledir = 0.0;
  guard1 = false;
  guard2 = false;
  if (c4_scale >= 7.4428285367870146E+137) {
    do {
      c4_a = c4_count;
      c4_b_a = c4_a + 1;
      c4_count = c4_b_a;
      c4_fs.re *= 1.3435752215134178E-138;
      c4_fs.im *= 1.3435752215134178E-138;
      c4_gs.re *= 1.3435752215134178E-138;
      c4_gs.im *= 1.3435752215134178E-138;
      c4_scale *= 1.3435752215134178E-138;
    } while (!(c4_scale < 7.4428285367870146E+137));

    c4_rescaledir = 1.0;
    guard1 = true;
  } else if (c4_scale <= 1.3435752215134178E-138) {
    c4_b_g = ((c4_g.re == c4_dc1.re) && (c4_g.im == c4_dc1.im));
    if (c4_b_g) {
      *c4_cs = 1.0;
      *c4_sn = c4_dc1;
      *c4_r = c4_f;
    } else {
      do {
        c4_c_a = c4_count;
        c4_d_a = c4_c_a + 1;
        c4_count = c4_d_a;
        c4_fs.re *= 7.4428285367870146E+137;
        c4_fs.im *= 7.4428285367870146E+137;
        c4_gs.re *= 7.4428285367870146E+137;
        c4_gs.im *= 7.4428285367870146E+137;
        c4_scale *= 7.4428285367870146E+137;
      } while (!(c4_scale > 1.3435752215134178E-138));

      c4_rescaledir = -1.0;
      guard2 = true;
    }
  } else {
    guard2 = true;
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c4_f2 = c4_fs.re * c4_fs.re + c4_fs.im * c4_fs.im;
    c4_g2 = c4_gs.re * c4_gs.re + c4_gs.im * c4_gs.im;
    c4_l_x = c4_g2;
    c4_m_x = c4_l_x;
    if (1.0 > c4_m_x) {
      c4_m_x = 1.0;
    }

    if (c4_f2 <= c4_m_x * 2.0041683600089728E-292) {
      c4_b_f = ((c4_f.re == c4_dc1.re) && (c4_f.im == c4_dc1.im));
      if (c4_b_f) {
        *c4_cs = 0.0;
        c4_x1 = c4_g.re;
        c4_x2 = c4_g.im;
        c4_e_a = c4_x1;
        c4_b = c4_x2;
        c4_j_y = muDoubleScalarHypot(c4_e_a, c4_b);
        c4_r->re = c4_j_y;
        c4_r->im = 0.0;
        c4_b_x1 = c4_gs.re;
        c4_b_x2 = c4_gs.im;
        c4_f_a = c4_b_x1;
        c4_b_b = c4_b_x2;
        c4_d = muDoubleScalarHypot(c4_f_a, c4_b_b);
        c4_sn->re = c4_gs.re / c4_d;
        c4_sn->im = -c4_gs.im / c4_d;
      } else {
        c4_c_x1 = c4_fs.re;
        c4_c_x2 = c4_fs.im;
        c4_g_a = c4_c_x1;
        c4_c_b = c4_c_x2;
        c4_f2s = muDoubleScalarHypot(c4_g_a, c4_c_b);
        c4_n_x = c4_g2;
        c4_g2s = c4_n_x;
        if (c4_g2s < 0.0) {
          c4_eml_error(chartInstance);
        }

        c4_g2s = muDoubleScalarSqrt(c4_g2s);
        *c4_cs = c4_f2s / c4_g2s;
        c4_o_x = c4_f.re;
        c4_p_x = c4_o_x;
        c4_k_y = muDoubleScalarAbs(c4_p_x);
        c4_q_x = c4_f.im;
        c4_r_x = c4_q_x;
        c4_l_y = muDoubleScalarAbs(c4_r_x);
        c4_s_x = c4_k_y;
        c4_m_y = c4_l_y;
        c4_n_y = c4_s_x;
        if (c4_m_y > c4_n_y) {
          c4_n_y = c4_m_y;
        }

        if (c4_n_y > 1.0) {
          c4_d_x1 = c4_f.re;
          c4_d_x2 = c4_f.im;
          c4_h_a = c4_d_x1;
          c4_d_b = c4_d_x2;
          c4_d = muDoubleScalarHypot(c4_h_a, c4_d_b);
          c4_fs.re = c4_f.re / c4_d;
          c4_fs.im = c4_f.im / c4_d;
        } else {
          c4_dr = 7.4428285367870146E+137 * c4_f.re;
          c4_di = 7.4428285367870146E+137 * c4_f.im;
          c4_e_x1 = c4_dr;
          c4_e_x2 = c4_di;
          c4_i_a = c4_e_x1;
          c4_e_b = c4_e_x2;
          c4_d = muDoubleScalarHypot(c4_i_a, c4_e_b);
          c4_fs.re = c4_dr / c4_d;
          c4_fs.im = c4_di / c4_d;
        }

        c4_b_gs.re = c4_gs.re / c4_g2s;
        c4_b_gs.im = -c4_gs.im / c4_g2s;
        c4_sn->re = c4_fs.re * c4_b_gs.re - c4_fs.im * c4_b_gs.im;
        c4_sn->im = c4_fs.re * c4_b_gs.im + c4_fs.im * c4_b_gs.re;
        c4_j_a = *c4_cs;
        c4_fs.re = c4_j_a * c4_f.re;
        c4_fs.im = c4_j_a * c4_f.im;
        c4_b_sn.re = c4_sn->re * c4_g.re - c4_sn->im * c4_g.im;
        c4_b_sn.im = c4_sn->re * c4_g.im + c4_sn->im * c4_g.re;
        c4_r->re = c4_fs.re + c4_b_sn.re;
        c4_r->im = c4_fs.im + c4_b_sn.im;
      }
    } else {
      c4_t_x = 1.0 + c4_g2 / c4_f2;
      c4_f2s = c4_t_x;
      if (c4_f2s < 0.0) {
        c4_eml_error(chartInstance);
      }

      c4_f2s = muDoubleScalarSqrt(c4_f2s);
      c4_r->re = c4_f2s * c4_fs.re;
      c4_r->im = c4_f2s * c4_fs.im;
      *c4_cs = 1.0 / c4_f2s;
      c4_d = c4_f2 + c4_g2;
      c4_sn->re = c4_r->re / c4_d;
      c4_sn->im = c4_r->im / c4_d;
      c4_c_gs.re = c4_gs.re;
      c4_c_gs.im = -c4_gs.im;
      c4_c_sn = *c4_sn;
      c4_sn->re = c4_c_sn.re * c4_c_gs.re - c4_c_sn.im * c4_c_gs.im;
      c4_sn->im = c4_c_sn.re * c4_c_gs.im + c4_c_sn.im * c4_c_gs.re;
      if (c4_rescaledir > 0.0) {
        c4_b_count = c4_count;
        c4_f_b = c4_b_count;
        c4_g_b = c4_f_b;
        if (1 > c4_g_b) {
          c4_overflow = false;
        } else {
          c4_eml_switch_helper(chartInstance);
          c4_overflow = (c4_g_b > 2147483646);
        }

        if (c4_overflow) {
          c4_check_forloop_overflow_error(chartInstance, c4_overflow);
        }

        for (c4_i = 1; c4_i <= c4_b_count; c4_i++) {
          c4_r->re *= 7.4428285367870146E+137;
          c4_r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (c4_rescaledir < 0.0) {
          c4_c_count = c4_count;
          c4_h_b = c4_c_count;
          c4_i_b = c4_h_b;
          if (1 > c4_i_b) {
            c4_b_overflow = false;
          } else {
            c4_eml_switch_helper(chartInstance);
            c4_b_overflow = (c4_i_b > 2147483646);
          }

          if (c4_b_overflow) {
            c4_check_forloop_overflow_error(chartInstance, c4_b_overflow);
          }

          for (c4_b_i = 1; c4_b_i <= c4_c_count; c4_b_i++) {
            c4_r->re *= 1.3435752215134178E-138;
            c4_r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

static void c4_eml_scalar_eg(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c4_eml_matlab_zhgeqz(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T c4_ilo, int32_T c4_ihi, real_T *c4_info, creal_T
  c4_alpha1[3], creal_T c4_beta1[3])
{
  static creal_T c4_dc2 = { 0.0, 0.0 };

  int32_T c4_i55;
  creal_T c4_b_A[9];
  int32_T c4_i56;
  int32_T c4_i57;
  static creal_T c4_dc3 = { 0.0, 0.0 };

  creal_T c4_eshift;
  creal_T c4_ctemp;
  creal_T c4_rho;
  int32_T c4_i58;
  creal_T c4_c_A[9];
  real_T c4_anorm;
  real_T c4_y;
  real_T c4_atol;
  real_T c4_b_y;
  real_T c4_x;
  real_T c4_ascale;
  boolean_T c4_failed;
  int32_T c4_a;
  int32_T c4_b_a;
  int32_T c4_i59;
  int32_T c4_c_a;
  int32_T c4_d_a;
  boolean_T c4_overflow;
  int32_T c4_j;
  int32_T c4_b_j;
  int32_T c4_ifirst;
  int32_T c4_istart;
  int32_T c4_ilast;
  int32_T c4_e_a;
  int32_T c4_f_a;
  int32_T c4_ilastm1;
  int32_T c4_ifrstm;
  int32_T c4_ilastm;
  int32_T c4_iiter;
  int32_T c4_g_a;
  int32_T c4_b;
  int32_T c4_h_a;
  int32_T c4_b_b;
  int32_T c4_c;
  int32_T c4_i_a;
  int32_T c4_j_a;
  int32_T c4_b_c;
  int32_T c4_c_b;
  int32_T c4_d_b;
  int32_T c4_maxit;
  boolean_T c4_goto50;
  boolean_T c4_goto60;
  boolean_T c4_goto70;
  boolean_T c4_goto90;
  int32_T c4_b_maxit;
  int32_T c4_e_b;
  int32_T c4_f_b;
  boolean_T c4_b_overflow;
  int32_T c4_jiter;
  creal_T c4_a22;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_c_y;
  real_T c4_d_x;
  real_T c4_e_x;
  real_T c4_d_y;
  real_T c4_e_y;
  int32_T c4_k_a;
  int32_T c4_l_a;
  int32_T c4_jm1;
  boolean_T c4_ilazro;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_f_y;
  real_T c4_h_x;
  real_T c4_i_x;
  real_T c4_g_y;
  real_T c4_h_y;
  boolean_T c4_b2;
  int32_T c4_i60;
  int32_T c4_i61;
  creal_T c4_d_A;
  creal_T c4_e_A;
  creal_T c4_s;
  real_T c4_c_c;
  real_T c4_d_c;
  real_T c4_e_c;
  int32_T c4_xcol;
  int32_T c4_ycol;
  int32_T c4_b_ilo;
  int32_T c4_b_ihi;
  int32_T c4_c_ilo;
  int32_T c4_c_ihi;
  int32_T c4_m_a;
  int32_T c4_g_b;
  int32_T c4_n_a;
  int32_T c4_h_b;
  boolean_T c4_c_overflow;
  int32_T c4_i;
  int32_T c4_b_i;
  real_T c4_o_a;
  creal_T c4_a12;
  creal_T c4_b_s;
  creal_T c4_a21;
  real_T c4_p_a;
  creal_T c4_b_a22;
  creal_T c4_c_a22;
  creal_T c4_d_a22;
  creal_T c4_e_a22;
  int32_T c4_q_a;
  int32_T c4_r_a;
  int32_T c4_s_a;
  int32_T c4_t_a;
  creal_T c4_r2;
  creal_T c4_f_a22;
  creal_T c4_b_rho;
  creal_T c4_b_a12;
  creal_T c4_c_a12;
  creal_T c4_b_a21;
  real_T c4_d4;
  real_T c4_d5;
  int32_T c4_u_a;
  int32_T c4_v_a;
  int32_T c4_jp1;
  int32_T c4_w_a;
  int32_T c4_x_a;
  real_T c4_j_x;
  real_T c4_k_x;
  real_T c4_i_y;
  real_T c4_l_x;
  real_T c4_m_x;
  real_T c4_j_y;
  real_T c4_k_y;
  real_T c4_temp;
  real_T c4_n_x;
  real_T c4_o_x;
  real_T c4_l_y;
  real_T c4_p_x;
  real_T c4_q_x;
  real_T c4_m_y;
  real_T c4_n_y;
  real_T c4_temp2;
  real_T c4_r_x;
  real_T c4_o_y;
  real_T c4_tempr;
  real_T c4_s_x;
  real_T c4_t_x;
  real_T c4_p_y;
  real_T c4_u_x;
  real_T c4_v_x;
  real_T c4_q_y;
  real_T c4_r_y;
  int32_T c4_y_a;
  int32_T c4_ab_a;
  int32_T c4_f_c;
  real_T c4_g_c;
  int32_T c4_bb_a;
  int32_T c4_cb_a;
  int32_T c4_db_a;
  int32_T c4_eb_a;
  creal_T c4_f_A;
  creal_T c4_g_A;
  real_T c4_h_c;
  real_T c4_i_c;
  int32_T c4_xrow;
  int32_T c4_yrow;
  int32_T c4_jlo;
  int32_T c4_jhi;
  int32_T c4_b_jlo;
  int32_T c4_b_jhi;
  int32_T c4_fb_a;
  int32_T c4_i_b;
  int32_T c4_gb_a;
  int32_T c4_j_b;
  boolean_T c4_d_overflow;
  int32_T c4_c_j;
  int32_T c4_d_j;
  real_T c4_hb_a;
  creal_T c4_c_s;
  real_T c4_ib_a;
  creal_T c4_g_a22;
  creal_T c4_h_a22;
  creal_T c4_i_a22;
  creal_T c4_j_a22;
  int32_T c4_jb_a;
  int32_T c4_kb_a;
  int32_T c4_j_c;
  int32_T c4_w_x;
  int32_T c4_s_y;
  int32_T c4_x_x;
  real_T c4_k_c;
  int32_T c4_b_xcol;
  int32_T c4_b_ycol;
  int32_T c4_d_ilo;
  int32_T c4_d_ihi;
  int32_T c4_e_ilo;
  int32_T c4_e_ihi;
  int32_T c4_lb_a;
  int32_T c4_k_b;
  int32_T c4_mb_a;
  int32_T c4_l_b;
  boolean_T c4_e_overflow;
  int32_T c4_c_i;
  int32_T c4_d_i;
  real_T c4_nb_a;
  creal_T c4_d_s;
  real_T c4_ob_a;
  creal_T c4_k_a22;
  creal_T c4_l_a22;
  creal_T c4_m_a22;
  creal_T c4_n_a22;
  int32_T c4_b_ilast;
  int32_T c4_m_b;
  int32_T c4_n_b;
  boolean_T c4_f_overflow;
  int32_T c4_k;
  int32_T c4_b_k;
  int32_T c4_pb_a;
  int32_T c4_qb_a;
  int32_T c4_i62;
  int32_T c4_o_b;
  int32_T c4_p_b;
  boolean_T c4_g_overflow;
  int32_T c4_e_j;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  int32_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T guard11 = false;
  c4_dc2.re = rtNaN;
  for (c4_i55 = 0; c4_i55 < 9; c4_i55++) {
    c4_b_A[c4_i55] = c4_A[c4_i55];
  }

  for (c4_i56 = 0; c4_i56 < 3; c4_i56++) {
    c4_alpha1[c4_i56].re = 0.0;
    c4_alpha1[c4_i56].im = 0.0;
  }

  for (c4_i57 = 0; c4_i57 < 3; c4_i57++) {
    c4_beta1[c4_i57].re = 1.0;
    c4_beta1[c4_i57].im = 0.0;
  }

  c4_eps(chartInstance);
  c4_realmin(chartInstance);
  c4_eshift = c4_dc3;
  c4_ctemp = c4_dc3;
  c4_rho = c4_dc3;
  for (c4_i58 = 0; c4_i58 < 9; c4_i58++) {
    c4_c_A[c4_i58] = c4_b_A[c4_i58];
  }

  c4_anorm = c4_eml_matlab_zlanhs(chartInstance, c4_c_A, c4_ilo, c4_ihi);
  c4_y = 2.2204460492503131E-16 * c4_anorm;
  c4_atol = 2.2250738585072014E-308;
  if (c4_y > 2.2250738585072014E-308) {
    c4_atol = c4_y;
  }

  c4_b_y = c4_anorm;
  c4_x = 2.2250738585072014E-308;
  if (c4_b_y > 2.2250738585072014E-308) {
    c4_x = c4_b_y;
  }

  c4_ascale = 1.0 / c4_x;
  c4_failed = true;
  c4_a = c4_ihi;
  c4_b_a = c4_a + 1;
  c4_i59 = c4_b_a;
  c4_c_a = c4_i59;
  c4_d_a = c4_c_a;
  if (c4_d_a > 3) {
    c4_overflow = false;
  } else {
    c4_eml_switch_helper(chartInstance);
    c4_overflow = false;
  }

  if (c4_overflow) {
    c4_check_forloop_overflow_error(chartInstance, c4_overflow);
  }

  for (c4_j = c4_i59; c4_j < 4; c4_j++) {
    c4_b_j = c4_j;
    c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_b_j), 1, 3, 1, 0) - 1].re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
      ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
    c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_b_j), 1, 3, 1, 0) - 1].im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
      ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
  }

  guard1 = false;
  guard2 = false;
  if (c4_ihi >= c4_ilo) {
    c4_ifirst = c4_ilo;
    c4_istart = c4_ilo;
    c4_ilast = c4_ihi;
    c4_e_a = c4_ilast;
    c4_f_a = c4_e_a - 1;
    c4_ilastm1 = c4_f_a;
    c4_ifrstm = c4_ilo;
    c4_ilastm = c4_ihi;
    c4_iiter = 0;
    c4_g_a = c4_ihi;
    c4_b = c4_ilo;
    c4_h_a = c4_g_a;
    c4_b_b = c4_b;
    c4_c = c4_h_a - c4_b_b;
    c4_i_a = c4_c;
    c4_j_a = c4_i_a;
    c4_b_c = c4_j_a;
    c4_c_b = c4_b_c + 1;
    c4_d_b = c4_c_b;
    c4_maxit = 30 * c4_d_b;
    c4_goto50 = false;
    c4_goto60 = false;
    c4_goto70 = false;
    c4_goto90 = false;
    c4_b_maxit = c4_maxit;
    c4_e_b = c4_b_maxit;
    c4_f_b = c4_e_b;
    if (1 > c4_f_b) {
      c4_b_overflow = false;
    } else {
      c4_eml_switch_helper(chartInstance);
      c4_b_overflow = (c4_f_b > 2147483646);
    }

    if (c4_b_overflow) {
      c4_check_forloop_overflow_error(chartInstance, c4_b_overflow);
    }

    c4_jiter = 1;
    do {
      exitg1 = 0;
      if (c4_jiter <= c4_b_maxit) {
        if (c4_ilast == c4_ilo) {
          c4_goto60 = true;
        } else {
          c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            re;
          c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            im;
          c4_b_x = c4_a22.re;
          c4_c_x = c4_b_x;
          c4_c_y = muDoubleScalarAbs(c4_c_x);
          c4_d_x = c4_a22.im;
          c4_e_x = c4_d_x;
          c4_d_y = muDoubleScalarAbs(c4_e_x);
          c4_e_y = c4_c_y + c4_d_y;
          if (c4_e_y <= c4_atol) {
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c4_dc3.re;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c4_dc3.im;
            c4_goto60 = true;
          } else {
            c4_b_j = c4_ilastm1;
            exitg3 = false;
            while ((exitg3 == false) && (c4_b_j >= c4_ilo)) {
              c4_k_a = c4_b_j;
              c4_l_a = c4_k_a - 1;
              c4_jm1 = c4_l_a;
              if (c4_b_j == c4_ilo) {
                c4_ilazro = true;
              } else {
                c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c4_f_x = c4_a22.re;
                c4_g_x = c4_f_x;
                c4_f_y = muDoubleScalarAbs(c4_g_x);
                c4_h_x = c4_a22.im;
                c4_i_x = c4_h_x;
                c4_g_y = muDoubleScalarAbs(c4_i_x);
                c4_h_y = c4_f_y + c4_g_y;
                if (c4_h_y <= c4_atol) {
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0)
                               - 1)) - 1].re = c4_dc3.re;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0)
                               - 1)) - 1].im = c4_dc3.im;
                  c4_ilazro = true;
                } else {
                  c4_ilazro = false;
                }
              }

              if (c4_ilazro) {
                c4_ifirst = c4_b_j;
                c4_goto70 = true;
                exitg3 = true;
              } else {
                c4_b_j = c4_jm1;
              }
            }
          }
        }

        guard3 = false;
        guard4 = false;
        if (c4_goto50) {
          guard4 = true;
        } else if (c4_goto60) {
          guard4 = true;
        } else if (c4_goto70) {
          guard3 = true;
        } else {
          c4_b2 = false;
        }

        if (guard4 == true) {
          guard3 = true;
        }

        if (guard3 == true) {
          c4_b2 = true;
        }

        if (!c4_b2) {
          for (c4_i60 = 0; c4_i60 < 3; c4_i60++) {
            c4_alpha1[c4_i60] = c4_dc2;
          }

          for (c4_i61 = 0; c4_i61 < 3; c4_i61++) {
            c4_beta1[c4_i61] = c4_dc2;
          }

          *c4_info = -1.0;
          exitg1 = 1;
        } else {
          if (c4_goto50) {
            c4_goto50 = false;
            c4_d_A.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].
              re;
            c4_d_A.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].
              im;
            c4_e_A.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .re;
            c4_e_A.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .im;
            c4_eml_matlab_zlartg(chartInstance, c4_d_A, c4_e_A, &c4_c_c, &c4_s,
                                 &c4_a22);
            c4_d_c = c4_c_c;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].re =
              c4_a22.re;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].im =
              c4_a22.im;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c4_dc3.re;
            c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c4_dc3.im;
            c4_e_c = c4_d_c;
            c4_xcol = c4_ilast;
            c4_ycol = c4_ilastm1;
            c4_b_ilo = c4_ifrstm;
            c4_b_ihi = c4_ilastm1;
            c4_c_ilo = c4_b_ilo;
            c4_c_ihi = c4_b_ihi;
            c4_m_a = c4_c_ilo;
            c4_g_b = c4_c_ihi;
            c4_n_a = c4_m_a;
            c4_h_b = c4_g_b;
            if (c4_n_a > c4_h_b) {
              c4_c_overflow = false;
            } else {
              c4_eml_switch_helper(chartInstance);
              c4_c_overflow = (c4_h_b > 2147483646);
            }

            if (c4_c_overflow) {
              c4_check_forloop_overflow_error(chartInstance, c4_c_overflow);
            }

            for (c4_i = c4_c_ilo; c4_i <= c4_c_ihi; c4_i++) {
              c4_b_i = c4_i;
              c4_o_a = c4_e_c;
              c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c4_a12.re = c4_o_a * c4_a22.re;
              c4_a12.im = c4_o_a * c4_a22.im;
              c4_b_s.re = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re - c4_s.im *
                c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) -
                           1)) - 1].im;
              c4_b_s.im = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im + c4_s.im *
                c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) -
                           1)) - 1].re;
              c4_a21.re = c4_a12.re + c4_b_s.re;
              c4_a21.im = c4_a12.im + c4_b_s.im;
              c4_p_a = c4_e_c;
              c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c4_a12.re = c4_p_a * c4_a22.re;
              c4_a12.im = c4_p_a * c4_a22.im;
              c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c4_b_a22 = c4_a22;
              c4_c_a22 = c4_a22;
              c4_d_a22 = c4_a22;
              c4_e_a22 = c4_a22;
              c4_a22.re = c4_s.re * c4_b_a22.re + c4_s.im * c4_c_a22.im;
              c4_a22.im = c4_s.re * c4_d_a22.im - c4_s.im * c4_e_a22.re;
              c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1))
                - 1].re = c4_a12.re - c4_a22.re;
              c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1))
                - 1].im = c4_a12.im - c4_a22.im;
              c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1))
                - 1].re = c4_a21.re;
              c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1))
                - 1].im = c4_a21.im;
            }

            c4_goto60 = true;
          }

          guard11 = false;
          if (c4_goto60) {
            c4_goto60 = false;
            c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) - 1].re =
              c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3
                      * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) -
                         1)) - 1].re;
            c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) - 1].im =
              c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3
                      * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) -
                         1)) - 1].im;
            c4_ilast = c4_ilastm1;
            c4_q_a = c4_ilast;
            c4_r_a = c4_q_a - 1;
            c4_ilastm1 = c4_r_a;
            if (c4_ilast < c4_ilo) {
              c4_failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              c4_iiter = 0;
              c4_eshift = c4_dc3;
              c4_ilastm = c4_ilast;
              if (c4_ifrstm > c4_ilast) {
                c4_ifrstm = c4_ilo;
              }

              guard11 = true;
            }
          } else {
            if (c4_goto70) {
              c4_goto70 = false;
              c4_s_a = c4_iiter;
              c4_t_a = c4_s_a + 1;
              c4_iiter = c4_t_a;
              c4_ifrstm = c4_ifirst;
              if (c4_mod(chartInstance, c4_iiter) != 0) {
                c4_s.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c4_s.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c4_r2.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                   (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c4_r2.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                   (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c4_a12.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c4_a12.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c4_a21.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c4_a21.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c4_a22.re = c4_r2.re - c4_s.re;
                c4_a22.im = c4_r2.im - c4_s.im;
                c4_f_a22.re = -c4_a22.re;
                c4_f_a22.im = -c4_a22.im;
                c4_rho = c4_eml_div(chartInstance, c4_f_a22, 2.0);
                c4_b_rho.re = c4_rho.re * c4_rho.re - c4_rho.im * c4_rho.im;
                c4_b_rho.im = c4_rho.re * c4_rho.im + c4_rho.im * c4_rho.re;
                c4_b_a12.re = c4_a12.re * c4_a21.re - c4_a12.im * c4_a21.im;
                c4_b_a12.im = c4_a12.re * c4_a21.im + c4_a12.im * c4_a21.re;
                c4_a22.re = c4_b_rho.re + c4_b_a12.re;
                c4_a22.im = c4_b_rho.im + c4_b_a12.im;
                c4_b_sqrt(chartInstance, &c4_a22);
                c4_a12.re = c4_s.re - (c4_rho.re - c4_a22.re);
                c4_a12.im = c4_s.im - (c4_rho.im - c4_a22.im);
                c4_a21.re = c4_s.re - (c4_rho.re + c4_a22.re);
                c4_a21.im = c4_s.im - (c4_rho.im + c4_a22.im);
                c4_c_a12.re = c4_a12.re - c4_r2.re;
                c4_c_a12.im = c4_a12.im - c4_r2.im;
                c4_b_a21.re = c4_a21.re - c4_r2.re;
                c4_b_a21.im = c4_a21.im - c4_r2.im;
                c4_d4 = c4_abs(chartInstance, c4_c_a12);
                c4_d5 = c4_abs(chartInstance, c4_b_a21);
                if (c4_d4 <= c4_d5) {
                  c4_a21 = c4_a12;
                  c4_rho.re -= c4_a22.re;
                  c4_rho.im -= c4_a22.im;
                } else {
                  c4_rho.re += c4_a22.re;
                  c4_rho.im += c4_a22.im;
                }
              } else {
                c4_eshift.re += c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                  "", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].re;
                c4_eshift.im += c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                  "", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].im;
                c4_a21 = c4_eshift;
              }

              c4_b_j = c4_ilastm1;
              c4_u_a = c4_b_j;
              c4_v_a = c4_u_a + 1;
              c4_jp1 = c4_v_a;
              exitg2 = false;
              while ((exitg2 == false) && (c4_b_j > c4_ifirst)) {
                c4_w_a = c4_b_j;
                c4_x_a = c4_w_a - 1;
                c4_jm1 = c4_x_a;
                c4_istart = c4_b_j;
                c4_ctemp.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re - c4_a21.re;
                c4_ctemp.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im - c4_a21.im;
                c4_j_x = c4_ctemp.re;
                c4_k_x = c4_j_x;
                c4_i_y = muDoubleScalarAbs(c4_k_x);
                c4_l_x = c4_ctemp.im;
                c4_m_x = c4_l_x;
                c4_j_y = muDoubleScalarAbs(c4_m_x);
                c4_k_y = c4_i_y + c4_j_y;
                c4_temp = c4_ascale * c4_k_y;
                c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c4_n_x = c4_a22.re;
                c4_o_x = c4_n_x;
                c4_l_y = muDoubleScalarAbs(c4_o_x);
                c4_p_x = c4_a22.im;
                c4_q_x = c4_p_x;
                c4_m_y = muDoubleScalarAbs(c4_q_x);
                c4_n_y = c4_l_y + c4_m_y;
                c4_temp2 = c4_ascale * c4_n_y;
                c4_r_x = c4_temp;
                c4_o_y = c4_temp2;
                c4_tempr = c4_r_x;
                if (c4_o_y > c4_tempr) {
                  c4_tempr = c4_o_y;
                }

                if (c4_tempr < 1.0) {
                  if (c4_tempr != 0.0) {
                    c4_temp /= c4_tempr;
                    c4_temp2 /= c4_tempr;
                  }
                }

                c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c4_s_x = c4_a22.re;
                c4_t_x = c4_s_x;
                c4_p_y = muDoubleScalarAbs(c4_t_x);
                c4_u_x = c4_a22.im;
                c4_v_x = c4_u_x;
                c4_q_y = muDoubleScalarAbs(c4_v_x);
                c4_r_y = c4_p_y + c4_q_y;
                if (c4_r_y * c4_temp2 <= c4_temp * c4_atol) {
                  c4_goto90 = true;
                  exitg2 = true;
                } else {
                  c4_jp1 = c4_b_j;
                  c4_b_j = c4_jm1;
                }
              }

              if (!c4_goto90) {
                c4_istart = c4_ifirst;
                if (c4_istart == c4_ilastm1) {
                  c4_ctemp = c4_rho;
                } else {
                  c4_ctemp.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 1, 0) + 3 *
                                        (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 2,
                    0) - 1)) - 1].re - c4_a21.re;
                  c4_ctemp.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 1, 0) + 3 *
                                        (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 2,
                    0) - 1)) - 1].im - c4_a21.im;
                }

                c4_goto90 = true;
              }
            }

            if (c4_goto90) {
              c4_goto90 = false;
              c4_y_a = c4_istart;
              c4_ab_a = c4_y_a;
              c4_f_c = c4_ab_a;
              c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c4_f_c + 1)), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c4_f_c + 1)), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c4_b_eml_matlab_zlartg(chartInstance, c4_ctemp, c4_a22, &c4_g_c,
                &c4_s);
              c4_d_c = c4_g_c;
              c4_b_j = c4_istart;
              c4_bb_a = c4_b_j;
              c4_cb_a = c4_bb_a - 1;
              c4_jm1 = c4_cb_a;
              while (c4_b_j < c4_ilast) {
                c4_db_a = c4_b_j;
                c4_eb_a = c4_db_a + 1;
                c4_jp1 = c4_eb_a;
                if (c4_b_j > c4_istart) {
                  c4_f_A.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_f_A.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_g_A.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_g_A.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_eml_matlab_zlartg(chartInstance, c4_f_A, c4_g_A, &c4_h_c,
                                       &c4_s, &c4_a22);
                  c4_d_c = c4_h_c;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0)
                               - 1)) - 1].re = c4_a22.re;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0)
                               - 1)) - 1].im = c4_a22.im;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0)
                               - 1)) - 1].re = c4_dc3.re;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0)
                               - 1)) - 1].im = c4_dc3.im;
                }

                c4_i_c = c4_d_c;
                c4_xrow = c4_b_j;
                c4_yrow = c4_jp1;
                c4_jlo = c4_b_j;
                c4_jhi = c4_ilastm;
                c4_b_jlo = c4_jlo;
                c4_b_jhi = c4_jhi;
                c4_fb_a = c4_b_jlo;
                c4_i_b = c4_b_jhi;
                c4_gb_a = c4_fb_a;
                c4_j_b = c4_i_b;
                if (c4_gb_a > c4_j_b) {
                  c4_d_overflow = false;
                } else {
                  c4_eml_switch_helper(chartInstance);
                  c4_d_overflow = (c4_j_b > 2147483646);
                }

                if (c4_d_overflow) {
                  c4_check_forloop_overflow_error(chartInstance, c4_d_overflow);
                }

                for (c4_c_j = c4_b_jlo; c4_c_j <= c4_b_jhi; c4_c_j++) {
                  c4_d_j = c4_c_j;
                  c4_hb_a = c4_i_c;
                  c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_a12.re = c4_hb_a * c4_a22.re;
                  c4_a12.im = c4_hb_a * c4_a22.im;
                  c4_c_s.re = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re - c4_s.im * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_c_s.im = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im + c4_s.im * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_a21.re = c4_a12.re + c4_c_s.re;
                  c4_a21.im = c4_a12.im + c4_c_s.im;
                  c4_ib_a = c4_i_c;
                  c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_a12.re = c4_ib_a * c4_a22.re;
                  c4_a12.im = c4_ib_a * c4_a22.im;
                  c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_g_a22 = c4_a22;
                  c4_h_a22 = c4_a22;
                  c4_i_a22 = c4_a22;
                  c4_j_a22 = c4_a22;
                  c4_a22.re = c4_s.re * c4_g_a22.re + c4_s.im * c4_h_a22.im;
                  c4_a22.im = c4_s.re * c4_i_a22.im - c4_s.im * c4_j_a22.re;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].re = c4_a12.re - c4_a22.re;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].im = c4_a12.im - c4_a22.im;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].re = c4_a21.re;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0)
                          + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0)
                                 - 1)) - 1].im = c4_a21.im;
                }

                c4_s.re = -c4_s.re;
                c4_s.im = -c4_s.im;
                c4_jb_a = c4_jp1;
                c4_kb_a = c4_jb_a;
                c4_j_c = c4_kb_a;
                c4_w_x = c4_j_c + 1;
                c4_s_y = c4_ilast;
                c4_x_x = c4_w_x;
                if (c4_s_y < c4_x_x) {
                  c4_x_x = c4_s_y;
                }

                c4_k_c = c4_d_c;
                c4_b_xcol = c4_jp1;
                c4_b_ycol = c4_b_j;
                c4_d_ilo = c4_ifrstm;
                c4_d_ihi = c4_x_x;
                c4_e_ilo = c4_d_ilo;
                c4_e_ihi = c4_d_ihi;
                c4_lb_a = c4_e_ilo;
                c4_k_b = c4_e_ihi;
                c4_mb_a = c4_lb_a;
                c4_l_b = c4_k_b;
                if (c4_mb_a > c4_l_b) {
                  c4_e_overflow = false;
                } else {
                  c4_eml_switch_helper(chartInstance);
                  c4_e_overflow = (c4_l_b > 2147483646);
                }

                if (c4_e_overflow) {
                  c4_check_forloop_overflow_error(chartInstance, c4_e_overflow);
                }

                for (c4_c_i = c4_e_ilo; c4_c_i <= c4_e_ihi; c4_c_i++) {
                  c4_d_i = c4_c_i;
                  c4_nb_a = c4_k_c;
                  c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_a12.re = c4_nb_a * c4_a22.re;
                  c4_a12.im = c4_nb_a * c4_a22.im;
                  c4_d_s.re = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].re - c4_s.im * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_d_s.im = c4_s.re * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].im + c4_s.im * c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a21.re = c4_a12.re + c4_d_s.re;
                  c4_a21.im = c4_a12.im + c4_d_s.im;
                  c4_ob_a = c4_k_c;
                  c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_a12.re = c4_ob_a * c4_a22.re;
                  c4_a12.im = c4_ob_a * c4_a22.im;
                  c4_a22.re = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a22.im = c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_k_a22 = c4_a22;
                  c4_l_a22 = c4_a22;
                  c4_m_a22 = c4_a22;
                  c4_n_a22 = c4_a22;
                  c4_a22.re = c4_s.re * c4_k_a22.re + c4_s.im * c4_l_a22.im;
                  c4_a22.im = c4_s.re * c4_m_a22.im - c4_s.im * c4_n_a22.re;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2,
                            0) - 1)) - 1].re = c4_a12.re - c4_a22.re;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2,
                            0) - 1)) - 1].im = c4_a12.im - c4_a22.im;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2,
                            0) - 1)) - 1].re = c4_a21.re;
                  c4_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                           _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) +
                          3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                            _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2,
                            0) - 1)) - 1].im = c4_a21.im;
                }

                c4_jm1 = c4_b_j;
                c4_b_j = c4_jp1;
              }
            }

            guard11 = true;
          }

          if (guard11 == true) {
            c4_jiter++;
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
    if (c4_failed) {
      *c4_info = (real_T)c4_ilast;
      c4_b_ilast = c4_ilast;
      c4_m_b = c4_b_ilast;
      c4_n_b = c4_m_b;
      if (1 > c4_n_b) {
        c4_f_overflow = false;
      } else {
        c4_eml_switch_helper(chartInstance);
        c4_f_overflow = (c4_n_b > 2147483646);
      }

      if (c4_f_overflow) {
        c4_check_forloop_overflow_error(chartInstance, c4_f_overflow);
      }

      for (c4_k = 1; c4_k <= c4_b_ilast; c4_k++) {
        c4_b_k = c4_k;
        c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c4_b_k), 1, 3, 1, 0) - 1].re = c4_dc2.re;
        c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c4_b_k), 1, 3, 1, 0) - 1].im = c4_dc2.im;
        c4_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c4_b_k), 1, 3, 1, 0) - 1].re = c4_dc2.re;
        c4_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c4_b_k), 1, 3, 1, 0) - 1].im = c4_dc2.im;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1 == true) {
    c4_pb_a = c4_ilo;
    c4_qb_a = c4_pb_a - 1;
    c4_i62 = c4_qb_a;
    c4_o_b = c4_i62;
    c4_p_b = c4_o_b;
    if (1 > c4_p_b) {
      c4_g_overflow = false;
    } else {
      c4_eml_switch_helper(chartInstance);
      c4_g_overflow = (c4_p_b > 2147483646);
    }

    if (c4_g_overflow) {
      c4_check_forloop_overflow_error(chartInstance, c4_g_overflow);
    }

    for (c4_e_j = 1; c4_e_j <= c4_i62; c4_e_j++) {
      c4_b_j = c4_e_j;
      c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c4_b_j), 1, 3, 1, 0) - 1].re = c4_b_A
        [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) -
        1].re;
      c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c4_b_j), 1, 3, 1, 0) - 1].im = c4_b_A
        [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) -
        1].im;
    }

    *c4_info = 0.0;
  }
}

static real_T c4_eml_matlab_zlanhs(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T c4_ilo, int32_T c4_ihi)
{
  real_T c4_f;
  real_T c4_scale;
  real_T c4_sumsq;
  boolean_T c4_firstNonZero;
  int32_T c4_b_ilo;
  int32_T c4_b_ihi;
  int32_T c4_a;
  int32_T c4_b;
  int32_T c4_b_a;
  int32_T c4_b_b;
  boolean_T c4_overflow;
  int32_T c4_j;
  int32_T c4_b_j;
  int32_T c4_c_ilo;
  int32_T c4_c_a;
  int32_T c4_d_a;
  int32_T c4_c;
  int32_T c4_x;
  int32_T c4_y;
  int32_T c4_i63;
  int32_T c4_e_a;
  int32_T c4_c_b;
  int32_T c4_f_a;
  int32_T c4_d_b;
  boolean_T c4_b_overflow;
  int32_T c4_i;
  int32_T c4_b_i;
  creal_T c4_Aij;
  real_T c4_reAij;
  real_T c4_imAij;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_temp1;
  real_T c4_temp2;
  real_T c4_d_x;
  real_T c4_e_x;
  real_T c4_f_x;
  real_T c4_g_x;
  c4_f = 0.0;
  if (c4_ilo > c4_ihi) {
  } else {
    c4_scale = 0.0;
    c4_sumsq = 0.0;
    c4_firstNonZero = true;
    c4_b_ilo = c4_ilo;
    c4_b_ihi = c4_ihi;
    c4_a = c4_b_ilo;
    c4_b = c4_b_ihi;
    c4_b_a = c4_a;
    c4_b_b = c4_b;
    if (c4_b_a > c4_b_b) {
      c4_overflow = false;
    } else {
      c4_eml_switch_helper(chartInstance);
      c4_overflow = (c4_b_b > 2147483646);
    }

    if (c4_overflow) {
      c4_check_forloop_overflow_error(chartInstance, c4_overflow);
    }

    for (c4_j = c4_b_ilo; c4_j <= c4_b_ihi; c4_j++) {
      c4_b_j = c4_j;
      c4_c_ilo = c4_ilo;
      c4_c_a = c4_b_j;
      c4_d_a = c4_c_a;
      c4_c = c4_d_a;
      c4_x = c4_c + 1;
      c4_y = c4_ihi;
      c4_i63 = c4_x;
      if (c4_y < c4_i63) {
        c4_i63 = c4_y;
      }

      c4_e_a = c4_c_ilo;
      c4_c_b = c4_i63;
      c4_f_a = c4_e_a;
      c4_d_b = c4_c_b;
      if (c4_f_a > c4_d_b) {
        c4_b_overflow = false;
      } else {
        c4_eml_switch_helper(chartInstance);
        c4_b_overflow = (c4_d_b > 2147483646);
      }

      if (c4_b_overflow) {
        c4_check_forloop_overflow_error(chartInstance, c4_b_overflow);
      }

      for (c4_i = c4_c_ilo; c4_i <= c4_i63; c4_i++) {
        c4_b_i = c4_i;
        c4_Aij.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
        c4_Aij.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
        c4_reAij = c4_Aij.re;
        c4_imAij = c4_Aij.im;
        if (c4_reAij != 0.0) {
          c4_b_x = c4_reAij;
          c4_c_x = c4_b_x;
          c4_temp1 = muDoubleScalarAbs(c4_c_x);
          if (c4_firstNonZero) {
            c4_sumsq = 1.0;
            c4_scale = c4_temp1;
            c4_firstNonZero = false;
          } else if (c4_scale < c4_temp1) {
            c4_temp2 = c4_scale / c4_temp1;
            c4_sumsq = 1.0 + c4_sumsq * c4_temp2 * c4_temp2;
            c4_scale = c4_temp1;
          } else {
            c4_temp2 = c4_temp1 / c4_scale;
            c4_sumsq += c4_temp2 * c4_temp2;
          }
        }

        if (c4_imAij != 0.0) {
          c4_d_x = c4_imAij;
          c4_e_x = c4_d_x;
          c4_temp1 = muDoubleScalarAbs(c4_e_x);
          if (c4_firstNonZero) {
            c4_sumsq = 1.0;
            c4_scale = c4_temp1;
            c4_firstNonZero = false;
          } else if (c4_scale < c4_temp1) {
            c4_temp2 = c4_scale / c4_temp1;
            c4_sumsq = 1.0 + c4_sumsq * c4_temp2 * c4_temp2;
            c4_scale = c4_temp1;
          } else {
            c4_temp2 = c4_temp1 / c4_scale;
            c4_sumsq += c4_temp2 * c4_temp2;
          }
        }
      }
    }

    c4_f_x = c4_sumsq;
    c4_g_x = c4_f_x;
    if (c4_g_x < 0.0) {
      c4_eml_error(chartInstance);
    }

    c4_g_x = muDoubleScalarSqrt(c4_g_x);
    c4_f = c4_scale * c4_g_x;
  }

  return c4_f;
}

static int32_T c4_mod(SFc4_Model_01InstanceStruct *chartInstance, int32_T c4_x)
{
  int32_T c4_b_x;
  int32_T c4_t;
  c4_b_x = c4_x;
  c4_t = c4_div_s32(chartInstance, c4_b_x, 10);
  c4_t *= 10;
  return c4_b_x - c4_t;
}

static creal_T c4_eml_div(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_x, real_T c4_y)
{
  creal_T c4_z;
  real_T c4_b_y;
  real_T c4_ar;
  real_T c4_ai;
  real_T c4_br;
  real_T c4_bi;
  real_T c4_brm;
  real_T c4_bim;
  real_T c4_s;
  real_T c4_d;
  real_T c4_nr;
  real_T c4_ni;
  real_T c4_sgnbr;
  real_T c4_sgnbi;
  (void)chartInstance;
  c4_b_y = c4_y;
  c4_ar = c4_x.re;
  c4_ai = c4_x.im;
  c4_br = c4_b_y;
  c4_bi = 0.0;
  if (c4_bi == 0.0) {
    if (c4_ai == 0.0) {
      c4_z.re = c4_ar / c4_br;
      c4_z.im = 0.0;
    } else if (c4_ar == 0.0) {
      c4_z.re = 0.0;
      c4_z.im = c4_ai / c4_br;
    } else {
      c4_z.re = c4_ar / c4_br;
      c4_z.im = c4_ai / c4_br;
    }
  } else if (c4_br == 0.0) {
    if (c4_ar == 0.0) {
      c4_z.re = c4_ai / c4_bi;
      c4_z.im = 0.0;
    } else if (c4_ai == 0.0) {
      c4_z.re = 0.0;
      c4_z.im = -(c4_ar / c4_bi);
    } else {
      c4_z.re = c4_ai / c4_bi;
      c4_z.im = -(c4_ar / c4_bi);
    }
  } else {
    c4_brm = muDoubleScalarAbs(c4_br);
    c4_bim = muDoubleScalarAbs(c4_bi);
    if (c4_brm > c4_bim) {
      c4_s = c4_bi / c4_br;
      c4_d = c4_br + c4_s * c4_bi;
      c4_nr = c4_ar + c4_s * c4_ai;
      c4_ni = c4_ai - c4_s * c4_ar;
      c4_z.re = c4_nr / c4_d;
      c4_z.im = c4_ni / c4_d;
    } else if (c4_bim == c4_brm) {
      if (c4_br > 0.0) {
        c4_sgnbr = 0.5;
      } else {
        c4_sgnbr = -0.5;
      }

      if (c4_bi > 0.0) {
        c4_sgnbi = 0.5;
      } else {
        c4_sgnbi = -0.5;
      }

      c4_nr = c4_ar * c4_sgnbr + c4_ai * c4_sgnbi;
      c4_ni = c4_ai * c4_sgnbr - c4_ar * c4_sgnbi;
      c4_z.re = c4_nr / c4_brm;
      c4_z.im = c4_ni / c4_brm;
    } else {
      c4_s = c4_br / c4_bi;
      c4_d = c4_bi + c4_s * c4_br;
      c4_nr = c4_s * c4_ar + c4_ai;
      c4_ni = c4_s * c4_ai - c4_ar;
      c4_z.re = c4_nr / c4_d;
      c4_z.im = c4_ni / c4_d;
    }
  }

  return c4_z;
}

static void c4_scalarEg(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static creal_T c4_sqrt(SFc4_Model_01InstanceStruct *chartInstance, creal_T c4_x)
{
  creal_T c4_b_x;
  c4_b_x = c4_x;
  c4_b_sqrt(chartInstance, &c4_b_x);
  return c4_b_x;
}

static void c4_realmax(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c4_b_eml_matlab_zlartg(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_f, creal_T c4_g, real_T *c4_cs, creal_T *c4_sn)
{
  real_T c4_x;
  real_T c4_b_x;
  real_T c4_y;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_b_y;
  real_T c4_e_x;
  real_T c4_c_y;
  real_T c4_d_y;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_e_y;
  real_T c4_h_x;
  real_T c4_i_x;
  real_T c4_f_y;
  real_T c4_j_x;
  real_T c4_g_y;
  real_T c4_h_y;
  real_T c4_k_x;
  real_T c4_i_y;
  real_T c4_scale;
  creal_T c4_fs;
  creal_T c4_gs;
  int32_T c4_count;
  real_T c4_rescaledir;
  int32_T c4_a;
  int32_T c4_b_a;
  static creal_T c4_dc4 = { 0.0, 0.0 };

  boolean_T c4_b_g;
  int32_T c4_c_a;
  int32_T c4_d_a;
  real_T c4_f2;
  real_T c4_g2;
  real_T c4_l_x;
  real_T c4_m_x;
  boolean_T c4_b_f;
  real_T c4_x1;
  real_T c4_x2;
  real_T c4_e_a;
  real_T c4_b;
  real_T c4_d;
  real_T c4_b_x1;
  real_T c4_b_x2;
  real_T c4_f_a;
  real_T c4_b_b;
  real_T c4_f2s;
  real_T c4_n_x;
  real_T c4_g2s;
  real_T c4_o_x;
  real_T c4_p_x;
  real_T c4_j_y;
  real_T c4_q_x;
  real_T c4_r_x;
  real_T c4_k_y;
  real_T c4_s_x;
  real_T c4_l_y;
  real_T c4_m_y;
  real_T c4_c_x1;
  real_T c4_c_x2;
  real_T c4_g_a;
  real_T c4_c_b;
  real_T c4_dr;
  real_T c4_di;
  real_T c4_d_x1;
  real_T c4_d_x2;
  real_T c4_h_a;
  real_T c4_d_b;
  creal_T c4_b_gs;
  real_T c4_t_x;
  creal_T c4_b_fs;
  creal_T c4_c_fs;
  creal_T c4_c_gs;
  creal_T c4_b_sn;
  int32_T c4_b_count;
  int32_T c4_e_b;
  int32_T c4_f_b;
  boolean_T c4_overflow;
  int32_T c4_c_count;
  int32_T c4_g_b;
  int32_T c4_h_b;
  boolean_T c4_b_overflow;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  c4_realmin(chartInstance);
  c4_eps(chartInstance);
  c4_eps(chartInstance);
  c4_x = c4_f.re;
  c4_b_x = c4_x;
  c4_y = muDoubleScalarAbs(c4_b_x);
  c4_c_x = c4_f.im;
  c4_d_x = c4_c_x;
  c4_b_y = muDoubleScalarAbs(c4_d_x);
  c4_e_x = c4_y;
  c4_c_y = c4_b_y;
  c4_d_y = c4_e_x;
  if (c4_c_y > c4_d_y) {
    c4_d_y = c4_c_y;
  }

  c4_f_x = c4_g.re;
  c4_g_x = c4_f_x;
  c4_e_y = muDoubleScalarAbs(c4_g_x);
  c4_h_x = c4_g.im;
  c4_i_x = c4_h_x;
  c4_f_y = muDoubleScalarAbs(c4_i_x);
  c4_j_x = c4_e_y;
  c4_g_y = c4_f_y;
  c4_h_y = c4_j_x;
  if (c4_g_y > c4_h_y) {
    c4_h_y = c4_g_y;
  }

  c4_k_x = c4_d_y;
  c4_i_y = c4_h_y;
  c4_scale = c4_k_x;
  if (c4_i_y > c4_scale) {
    c4_scale = c4_i_y;
  }

  c4_fs = c4_f;
  c4_gs = c4_g;
  c4_count = 0;
  c4_rescaledir = 0.0;
  guard1 = false;
  guard2 = false;
  if (c4_scale >= 7.4428285367870146E+137) {
    do {
      c4_a = c4_count;
      c4_b_a = c4_a + 1;
      c4_count = c4_b_a;
      c4_fs.re *= 1.3435752215134178E-138;
      c4_fs.im *= 1.3435752215134178E-138;
      c4_gs.re *= 1.3435752215134178E-138;
      c4_gs.im *= 1.3435752215134178E-138;
      c4_scale *= 1.3435752215134178E-138;
    } while (!(c4_scale < 7.4428285367870146E+137));

    c4_rescaledir = 1.0;
    guard1 = true;
  } else if (c4_scale <= 1.3435752215134178E-138) {
    c4_b_g = ((c4_g.re == c4_dc4.re) && (c4_g.im == c4_dc4.im));
    if (c4_b_g) {
      *c4_cs = 1.0;
      *c4_sn = c4_dc4;
    } else {
      do {
        c4_c_a = c4_count;
        c4_d_a = c4_c_a + 1;
        c4_count = c4_d_a;
        c4_fs.re *= 7.4428285367870146E+137;
        c4_fs.im *= 7.4428285367870146E+137;
        c4_gs.re *= 7.4428285367870146E+137;
        c4_gs.im *= 7.4428285367870146E+137;
        c4_scale *= 7.4428285367870146E+137;
      } while (!(c4_scale > 1.3435752215134178E-138));

      c4_rescaledir = -1.0;
      guard2 = true;
    }
  } else {
    guard2 = true;
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c4_f2 = c4_fs.re * c4_fs.re + c4_fs.im * c4_fs.im;
    c4_g2 = c4_gs.re * c4_gs.re + c4_gs.im * c4_gs.im;
    c4_l_x = c4_g2;
    c4_m_x = c4_l_x;
    if (1.0 > c4_m_x) {
      c4_m_x = 1.0;
    }

    if (c4_f2 <= c4_m_x * 2.0041683600089728E-292) {
      c4_b_f = ((c4_f.re == c4_dc4.re) && (c4_f.im == c4_dc4.im));
      if (c4_b_f) {
        *c4_cs = 0.0;
        c4_x1 = c4_gs.re;
        c4_x2 = c4_gs.im;
        c4_e_a = c4_x1;
        c4_b = c4_x2;
        c4_d = muDoubleScalarHypot(c4_e_a, c4_b);
        c4_sn->re = c4_gs.re / c4_d;
        c4_sn->im = -c4_gs.im / c4_d;
      } else {
        c4_b_x1 = c4_fs.re;
        c4_b_x2 = c4_fs.im;
        c4_f_a = c4_b_x1;
        c4_b_b = c4_b_x2;
        c4_f2s = muDoubleScalarHypot(c4_f_a, c4_b_b);
        c4_n_x = c4_g2;
        c4_g2s = c4_n_x;
        if (c4_g2s < 0.0) {
          c4_eml_error(chartInstance);
        }

        c4_g2s = muDoubleScalarSqrt(c4_g2s);
        *c4_cs = c4_f2s / c4_g2s;
        c4_o_x = c4_f.re;
        c4_p_x = c4_o_x;
        c4_j_y = muDoubleScalarAbs(c4_p_x);
        c4_q_x = c4_f.im;
        c4_r_x = c4_q_x;
        c4_k_y = muDoubleScalarAbs(c4_r_x);
        c4_s_x = c4_j_y;
        c4_l_y = c4_k_y;
        c4_m_y = c4_s_x;
        if (c4_l_y > c4_m_y) {
          c4_m_y = c4_l_y;
        }

        if (c4_m_y > 1.0) {
          c4_c_x1 = c4_f.re;
          c4_c_x2 = c4_f.im;
          c4_g_a = c4_c_x1;
          c4_c_b = c4_c_x2;
          c4_d = muDoubleScalarHypot(c4_g_a, c4_c_b);
          c4_fs.re = c4_f.re / c4_d;
          c4_fs.im = c4_f.im / c4_d;
        } else {
          c4_dr = 7.4428285367870146E+137 * c4_f.re;
          c4_di = 7.4428285367870146E+137 * c4_f.im;
          c4_d_x1 = c4_dr;
          c4_d_x2 = c4_di;
          c4_h_a = c4_d_x1;
          c4_d_b = c4_d_x2;
          c4_d = muDoubleScalarHypot(c4_h_a, c4_d_b);
          c4_fs.re = c4_dr / c4_d;
          c4_fs.im = c4_di / c4_d;
        }

        c4_b_gs.re = c4_gs.re / c4_g2s;
        c4_b_gs.im = -c4_gs.im / c4_g2s;
        c4_sn->re = c4_fs.re * c4_b_gs.re - c4_fs.im * c4_b_gs.im;
        c4_sn->im = c4_fs.re * c4_b_gs.im + c4_fs.im * c4_b_gs.re;
      }
    } else {
      c4_t_x = 1.0 + c4_g2 / c4_f2;
      c4_f2s = c4_t_x;
      if (c4_f2s < 0.0) {
        c4_eml_error(chartInstance);
      }

      c4_f2s = muDoubleScalarSqrt(c4_f2s);
      c4_b_fs = c4_fs;
      c4_c_fs = c4_fs;
      c4_fs.re = c4_f2s * c4_b_fs.re;
      c4_fs.im = c4_f2s * c4_c_fs.im;
      *c4_cs = 1.0 / c4_f2s;
      c4_d = c4_f2 + c4_g2;
      c4_sn->re = c4_fs.re / c4_d;
      c4_sn->im = c4_fs.im / c4_d;
      c4_c_gs.re = c4_gs.re;
      c4_c_gs.im = -c4_gs.im;
      c4_b_sn = *c4_sn;
      c4_sn->re = c4_b_sn.re * c4_c_gs.re - c4_b_sn.im * c4_c_gs.im;
      c4_sn->im = c4_b_sn.re * c4_c_gs.im + c4_b_sn.im * c4_c_gs.re;
      if (c4_rescaledir > 0.0) {
        c4_b_count = c4_count;
        c4_e_b = c4_b_count;
        c4_f_b = c4_e_b;
        if (1 > c4_f_b) {
          c4_overflow = false;
        } else {
          c4_eml_switch_helper(chartInstance);
          c4_overflow = (c4_f_b > 2147483646);
        }

        if (c4_overflow) {
          c4_check_forloop_overflow_error(chartInstance, c4_overflow);
        }
      } else {
        if (c4_rescaledir < 0.0) {
          c4_c_count = c4_count;
          c4_g_b = c4_c_count;
          c4_h_b = c4_g_b;
          if (1 > c4_h_b) {
            c4_b_overflow = false;
          } else {
            c4_eml_switch_helper(chartInstance);
            c4_b_overflow = (c4_h_b > 2147483646);
          }

          if (c4_b_overflow) {
            c4_check_forloop_overflow_error(chartInstance, c4_b_overflow);
          }
        }
      }
    }
  }
}

static void c4_b_eml_matlab_zlascl(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_cfrom, real_T c4_cto, creal_T c4_A[3], creal_T c4_b_A[3])
{
  int32_T c4_i64;
  for (c4_i64 = 0; c4_i64 < 3; c4_i64++) {
    c4_b_A[c4_i64] = c4_A[c4_i64];
  }

  c4_d_eml_matlab_zlascl(chartInstance, c4_cfrom, c4_cto, c4_b_A);
}

static void c4_b_eml_div(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_x[3], creal_T c4_y[3], creal_T c4_z[3])
{
  int32_T c4_i65;
  real_T c4_ar;
  real_T c4_ai;
  real_T c4_br;
  real_T c4_bi;
  real_T c4_brm;
  real_T c4_bim;
  real_T c4_s;
  real_T c4_d;
  real_T c4_nr;
  real_T c4_ni;
  real_T c4_sgnbr;
  real_T c4_sgnbi;
  (void)chartInstance;
  for (c4_i65 = 0; c4_i65 < 3; c4_i65++) {
    c4_ar = c4_x[c4_i65].re;
    c4_ai = c4_x[c4_i65].im;
    c4_br = c4_y[c4_i65].re;
    c4_bi = c4_y[c4_i65].im;
    if (c4_bi == 0.0) {
      if (c4_ai == 0.0) {
        c4_z[c4_i65].re = c4_ar / c4_br;
        c4_z[c4_i65].im = 0.0;
      } else if (c4_ar == 0.0) {
        c4_z[c4_i65].re = 0.0;
        c4_z[c4_i65].im = c4_ai / c4_br;
      } else {
        c4_z[c4_i65].re = c4_ar / c4_br;
        c4_z[c4_i65].im = c4_ai / c4_br;
      }
    } else if (c4_br == 0.0) {
      if (c4_ar == 0.0) {
        c4_z[c4_i65].re = c4_ai / c4_bi;
        c4_z[c4_i65].im = 0.0;
      } else if (c4_ai == 0.0) {
        c4_z[c4_i65].re = 0.0;
        c4_z[c4_i65].im = -(c4_ar / c4_bi);
      } else {
        c4_z[c4_i65].re = c4_ai / c4_bi;
        c4_z[c4_i65].im = -(c4_ar / c4_bi);
      }
    } else {
      c4_brm = muDoubleScalarAbs(c4_br);
      c4_bim = muDoubleScalarAbs(c4_bi);
      if (c4_brm > c4_bim) {
        c4_s = c4_bi / c4_br;
        c4_d = c4_br + c4_s * c4_bi;
        c4_nr = c4_ar + c4_s * c4_ai;
        c4_ni = c4_ai - c4_s * c4_ar;
        c4_z[c4_i65].re = c4_nr / c4_d;
        c4_z[c4_i65].im = c4_ni / c4_d;
      } else if (c4_bim == c4_brm) {
        if (c4_br > 0.0) {
          c4_sgnbr = 0.5;
        } else {
          c4_sgnbr = -0.5;
        }

        if (c4_bi > 0.0) {
          c4_sgnbi = 0.5;
        } else {
          c4_sgnbi = -0.5;
        }

        c4_nr = c4_ar * c4_sgnbr + c4_ai * c4_sgnbi;
        c4_ni = c4_ai * c4_sgnbr - c4_ar * c4_sgnbi;
        c4_z[c4_i65].re = c4_nr / c4_brm;
        c4_z[c4_i65].im = c4_ni / c4_brm;
      } else {
        c4_s = c4_br / c4_bi;
        c4_d = c4_bi + c4_s * c4_br;
        c4_nr = c4_s * c4_ar + c4_ai;
        c4_ni = c4_s * c4_ai - c4_ar;
        c4_z[c4_i65].re = c4_nr / c4_d;
        c4_z[c4_i65].im = c4_ni / c4_d;
      }
    }
  }
}

static void c4_eml_warning(SFc4_Model_01InstanceStruct *chartInstance)
{
  int32_T c4_i66;
  static char_T c4_varargin_1[26] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'f', 'a', 'i',
    'l', 'e', 'd' };

  char_T c4_u[26];
  const mxArray *c4_y = NULL;
  (void)chartInstance;
  for (c4_i66 = 0; c4_i66 < 26; c4_i66++) {
    c4_u[c4_i66] = c4_varargin_1[c4_i66];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 26), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c4_y));
}

static void c4_b_eml_warning(SFc4_Model_01InstanceStruct *chartInstance)
{
  int32_T c4_i67;
  static char_T c4_varargin_1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'e', 'i', 'g', '_', 'Q', 'Z', 'n', 'o', 'n',
    'c', 'o', 'n', 'v', 'e', 'r', 'g', 'e', 'n', 'c', 'e' };

  char_T c4_u[34];
  const mxArray *c4_y = NULL;
  (void)chartInstance;
  for (c4_i67 = 0; c4_i67 < 34; c4_i67++) {
    c4_u[c4_i67] = c4_varargin_1[c4_i67];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 34), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c4_y));
}

static real_T c4_mpower(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_a)
{
  real_T c4_b_a;
  real_T c4_c_a;
  real_T c4_ak;
  real_T c4_d_a;
  c4_b_a = c4_a;
  c4_c_a = c4_b_a;
  c4_eml_scalar_eg(chartInstance);
  c4_ak = c4_c_a;
  c4_d_a = c4_ak;
  c4_eml_scalar_eg(chartInstance);
  return c4_d_a * c4_d_a;
}

static void c4_inv(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_x[9],
                   real_T c4_y[9])
{
  int32_T c4_i68;
  real_T c4_b_x[9];
  int32_T c4_i69;
  real_T c4_c_x[9];
  real_T c4_n1x;
  int32_T c4_i70;
  real_T c4_b_y[9];
  real_T c4_n1xinv;
  real_T c4_rc;
  real_T c4_d_x;
  boolean_T c4_b;
  real_T c4_e_x;
  int32_T c4_i71;
  static char_T c4_cv4[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T c4_u[8];
  const mxArray *c4_c_y = NULL;
  real_T c4_b_u;
  const mxArray *c4_d_y = NULL;
  real_T c4_c_u;
  const mxArray *c4_e_y = NULL;
  real_T c4_d_u;
  const mxArray *c4_f_y = NULL;
  char_T c4_str[14];
  int32_T c4_i72;
  char_T c4_b_str[14];
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  for (c4_i68 = 0; c4_i68 < 9; c4_i68++) {
    c4_b_x[c4_i68] = c4_x[c4_i68];
  }

  c4_inv3x3(chartInstance, c4_b_x, c4_y);
  for (c4_i69 = 0; c4_i69 < 9; c4_i69++) {
    c4_c_x[c4_i69] = c4_x[c4_i69];
  }

  c4_n1x = c4_norm(chartInstance, c4_c_x);
  for (c4_i70 = 0; c4_i70 < 9; c4_i70++) {
    c4_b_y[c4_i70] = c4_y[c4_i70];
  }

  c4_n1xinv = c4_norm(chartInstance, c4_b_y);
  c4_rc = 1.0 / (c4_n1x * c4_n1xinv);
  guard1 = false;
  guard2 = false;
  if (c4_n1x == 0.0) {
    guard2 = true;
  } else if (c4_n1xinv == 0.0) {
    guard2 = true;
  } else if (c4_rc == 0.0) {
    guard1 = true;
  } else {
    c4_d_x = c4_rc;
    c4_b = muDoubleScalarIsNaN(c4_d_x);
    guard3 = false;
    if (c4_b) {
      guard3 = true;
    } else {
      c4_eps(chartInstance);
      if (c4_rc < 2.2204460492503131E-16) {
        guard3 = true;
      }
    }

    if (guard3 == true) {
      c4_e_x = c4_rc;
      for (c4_i71 = 0; c4_i71 < 8; c4_i71++) {
        c4_u[c4_i71] = c4_cv4[c4_i71];
      }

      c4_c_y = NULL;
      sf_mex_assign(&c4_c_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 8),
                    false);
      c4_b_u = 14.0;
      c4_d_y = NULL;
      sf_mex_assign(&c4_d_y, sf_mex_create("y", &c4_b_u, 0, 0U, 0U, 0U, 0),
                    false);
      c4_c_u = 6.0;
      c4_e_y = NULL;
      sf_mex_assign(&c4_e_y, sf_mex_create("y", &c4_c_u, 0, 0U, 0U, 0U, 0),
                    false);
      c4_d_u = c4_e_x;
      c4_f_y = NULL;
      sf_mex_assign(&c4_f_y, sf_mex_create("y", &c4_d_u, 0, 0U, 0U, 0U, 0),
                    false);
      c4_e_emlrt_marshallIn(chartInstance, sf_mex_call_debug
                            (sfGlobalDebugInstanceStruct, "sprintf", 1U, 2U, 14,
        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "sprintf", 1U, 3U, 14,
                          c4_c_y, 14, c4_d_y, 14, c4_e_y), 14, c4_f_y),
                            "sprintf", c4_str);
      for (c4_i72 = 0; c4_i72 < 14; c4_i72++) {
        c4_b_str[c4_i72] = c4_str[c4_i72];
      }

      c4_d_eml_warning(chartInstance, c4_b_str);
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c4_c_eml_warning(chartInstance);
  }
}

static void c4_inv3x3(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_x[9],
                      real_T c4_y[9])
{
  int32_T c4_p1;
  int32_T c4_p2;
  int32_T c4_p3;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_absx11;
  real_T c4_d_x;
  real_T c4_e_x;
  real_T c4_absx21;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_absx31;
  real_T c4_t1;
  real_T c4_h_x;
  real_T c4_b_y;
  real_T c4_i_x;
  real_T c4_c_y;
  real_T c4_z;
  real_T c4_j_x;
  real_T c4_d_y;
  real_T c4_k_x;
  real_T c4_e_y;
  real_T c4_b_z;
  real_T c4_l_x;
  real_T c4_m_x;
  real_T c4_f_y;
  real_T c4_n_x;
  real_T c4_o_x;
  real_T c4_g_y;
  int32_T c4_itmp;
  real_T c4_p_x;
  real_T c4_h_y;
  real_T c4_q_x;
  real_T c4_i_y;
  real_T c4_c_z;
  real_T c4_r_x;
  real_T c4_j_y;
  real_T c4_s_x;
  real_T c4_k_y;
  real_T c4_t3;
  real_T c4_t_x;
  real_T c4_l_y;
  real_T c4_u_x;
  real_T c4_m_y;
  real_T c4_t2;
  int32_T c4_a;
  int32_T c4_b_a;
  int32_T c4_c;
  real_T c4_v_x;
  real_T c4_n_y;
  real_T c4_w_x;
  real_T c4_o_y;
  real_T c4_d_z;
  int32_T c4_c_a;
  int32_T c4_d_a;
  int32_T c4_b_c;
  int32_T c4_e_a;
  int32_T c4_f_a;
  int32_T c4_c_c;
  real_T c4_x_x;
  real_T c4_p_y;
  real_T c4_y_x;
  real_T c4_q_y;
  real_T c4_ab_x;
  real_T c4_r_y;
  real_T c4_bb_x;
  real_T c4_s_y;
  int32_T c4_g_a;
  int32_T c4_h_a;
  int32_T c4_d_c;
  real_T c4_cb_x;
  real_T c4_t_y;
  real_T c4_db_x;
  real_T c4_u_y;
  real_T c4_e_z;
  int32_T c4_i_a;
  int32_T c4_j_a;
  int32_T c4_e_c;
  int32_T c4_k_a;
  int32_T c4_l_a;
  int32_T c4_f_c;
  real_T c4_v_y;
  real_T c4_w_y;
  real_T c4_eb_x;
  real_T c4_x_y;
  real_T c4_fb_x;
  real_T c4_y_y;
  int32_T c4_m_a;
  int32_T c4_n_a;
  int32_T c4_g_c;
  real_T c4_gb_x;
  real_T c4_ab_y;
  real_T c4_hb_x;
  real_T c4_bb_y;
  real_T c4_f_z;
  int32_T c4_o_a;
  int32_T c4_p_a;
  int32_T c4_h_c;
  int32_T c4_q_a;
  int32_T c4_r_a;
  int32_T c4_i_c;
  boolean_T guard1 = false;
  (void)chartInstance;
  c4_p1 = 0;
  c4_p2 = 3;
  c4_p3 = 6;
  c4_b_x = c4_x[0];
  c4_c_x = c4_b_x;
  c4_absx11 = muDoubleScalarAbs(c4_c_x);
  c4_d_x = c4_x[1];
  c4_e_x = c4_d_x;
  c4_absx21 = muDoubleScalarAbs(c4_e_x);
  c4_f_x = c4_x[2];
  c4_g_x = c4_f_x;
  c4_absx31 = muDoubleScalarAbs(c4_g_x);
  guard1 = false;
  if (c4_absx21 > c4_absx11) {
    if (c4_absx21 > c4_absx31) {
      c4_p1 = 3;
      c4_p2 = 0;
      c4_t1 = c4_x[0];
      c4_x[0] = c4_x[1];
      c4_x[1] = c4_t1;
      c4_t1 = c4_x[3];
      c4_x[3] = c4_x[4];
      c4_x[4] = c4_t1;
      c4_t1 = c4_x[6];
      c4_x[6] = c4_x[7];
      c4_x[7] = c4_t1;
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1 == true) {
    if (c4_absx31 > c4_absx11) {
      c4_p1 = 6;
      c4_p3 = 0;
      c4_t1 = c4_x[0];
      c4_x[0] = c4_x[2];
      c4_x[2] = c4_t1;
      c4_t1 = c4_x[3];
      c4_x[3] = c4_x[5];
      c4_x[5] = c4_t1;
      c4_t1 = c4_x[6];
      c4_x[6] = c4_x[8];
      c4_x[8] = c4_t1;
    }
  }

  c4_h_x = c4_x[1];
  c4_b_y = c4_x[0];
  c4_i_x = c4_h_x;
  c4_c_y = c4_b_y;
  c4_z = c4_i_x / c4_c_y;
  c4_x[1] = c4_z;
  c4_j_x = c4_x[2];
  c4_d_y = c4_x[0];
  c4_k_x = c4_j_x;
  c4_e_y = c4_d_y;
  c4_b_z = c4_k_x / c4_e_y;
  c4_x[2] = c4_b_z;
  c4_x[4] -= c4_x[1] * c4_x[3];
  c4_x[5] -= c4_x[2] * c4_x[3];
  c4_x[7] -= c4_x[1] * c4_x[6];
  c4_x[8] -= c4_x[2] * c4_x[6];
  c4_l_x = c4_x[5];
  c4_m_x = c4_l_x;
  c4_f_y = muDoubleScalarAbs(c4_m_x);
  c4_n_x = c4_x[4];
  c4_o_x = c4_n_x;
  c4_g_y = muDoubleScalarAbs(c4_o_x);
  if (c4_f_y > c4_g_y) {
    c4_itmp = c4_p2;
    c4_p2 = c4_p3;
    c4_p3 = c4_itmp;
    c4_t1 = c4_x[1];
    c4_x[1] = c4_x[2];
    c4_x[2] = c4_t1;
    c4_t1 = c4_x[4];
    c4_x[4] = c4_x[5];
    c4_x[5] = c4_t1;
    c4_t1 = c4_x[7];
    c4_x[7] = c4_x[8];
    c4_x[8] = c4_t1;
  }

  c4_p_x = c4_x[5];
  c4_h_y = c4_x[4];
  c4_q_x = c4_p_x;
  c4_i_y = c4_h_y;
  c4_c_z = c4_q_x / c4_i_y;
  c4_x[5] = c4_c_z;
  c4_x[8] -= c4_x[5] * c4_x[7];
  c4_r_x = c4_x[5] * c4_x[1] - c4_x[2];
  c4_j_y = c4_x[8];
  c4_s_x = c4_r_x;
  c4_k_y = c4_j_y;
  c4_t3 = c4_s_x / c4_k_y;
  c4_t_x = -(c4_x[1] + c4_x[7] * c4_t3);
  c4_l_y = c4_x[4];
  c4_u_x = c4_t_x;
  c4_m_y = c4_l_y;
  c4_t2 = c4_u_x / c4_m_y;
  c4_a = c4_p1;
  c4_b_a = c4_a + 1;
  c4_c = c4_b_a;
  c4_v_x = (1.0 - c4_x[3] * c4_t2) - c4_x[6] * c4_t3;
  c4_n_y = c4_x[0];
  c4_w_x = c4_v_x;
  c4_o_y = c4_n_y;
  c4_d_z = c4_w_x / c4_o_y;
  c4_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_c), 1, 9, 1, 0) - 1] = c4_d_z;
  c4_c_a = c4_p1;
  c4_d_a = c4_c_a + 2;
  c4_b_c = c4_d_a;
  c4_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_b_c), 1, 9, 1, 0) - 1] = c4_t2;
  c4_e_a = c4_p1;
  c4_f_a = c4_e_a + 3;
  c4_c_c = c4_f_a;
  c4_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_c_c), 1, 9, 1, 0) - 1] = c4_t3;
  c4_x_x = -c4_x[5];
  c4_p_y = c4_x[8];
  c4_y_x = c4_x_x;
  c4_q_y = c4_p_y;
  c4_t3 = c4_y_x / c4_q_y;
  c4_ab_x = 1.0 - c4_x[7] * c4_t3;
  c4_r_y = c4_x[4];
  c4_bb_x = c4_ab_x;
  c4_s_y = c4_r_y;
  c4_t2 = c4_bb_x / c4_s_y;
  c4_g_a = c4_p2;
  c4_h_a = c4_g_a + 1;
  c4_d_c = c4_h_a;
  c4_cb_x = -(c4_x[3] * c4_t2 + c4_x[6] * c4_t3);
  c4_t_y = c4_x[0];
  c4_db_x = c4_cb_x;
  c4_u_y = c4_t_y;
  c4_e_z = c4_db_x / c4_u_y;
  c4_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_d_c), 1, 9, 1, 0) - 1] = c4_e_z;
  c4_i_a = c4_p2;
  c4_j_a = c4_i_a + 2;
  c4_e_c = c4_j_a;
  c4_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_e_c), 1, 9, 1, 0) - 1] = c4_t2;
  c4_k_a = c4_p2;
  c4_l_a = c4_k_a + 3;
  c4_f_c = c4_l_a;
  c4_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_f_c), 1, 9, 1, 0) - 1] = c4_t3;
  c4_v_y = c4_x[8];
  c4_w_y = c4_v_y;
  c4_t3 = 1.0 / c4_w_y;
  c4_eb_x = -c4_x[7] * c4_t3;
  c4_x_y = c4_x[4];
  c4_fb_x = c4_eb_x;
  c4_y_y = c4_x_y;
  c4_t2 = c4_fb_x / c4_y_y;
  c4_m_a = c4_p3;
  c4_n_a = c4_m_a + 1;
  c4_g_c = c4_n_a;
  c4_gb_x = -(c4_x[3] * c4_t2 + c4_x[6] * c4_t3);
  c4_ab_y = c4_x[0];
  c4_hb_x = c4_gb_x;
  c4_bb_y = c4_ab_y;
  c4_f_z = c4_hb_x / c4_bb_y;
  c4_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_g_c), 1, 9, 1, 0) - 1] = c4_f_z;
  c4_o_a = c4_p3;
  c4_p_a = c4_o_a + 2;
  c4_h_c = c4_p_a;
  c4_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_h_c), 1, 9, 1, 0) - 1] = c4_t2;
  c4_q_a = c4_p3;
  c4_r_a = c4_q_a + 3;
  c4_i_c = c4_r_a;
  c4_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_i_c), 1, 9, 1, 0) - 1] = c4_t3;
}

static real_T c4_norm(SFc4_Model_01InstanceStruct *chartInstance, real_T c4_x[9])
{
  real_T c4_y;
  int32_T c4_j;
  real_T c4_b_j;
  real_T c4_s;
  int32_T c4_i;
  real_T c4_b_i;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_b_y;
  real_T c4_d_x;
  boolean_T c4_b;
  boolean_T exitg1;
  (void)chartInstance;
  c4_y = 0.0;
  c4_j = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c4_j < 3)) {
    c4_b_j = 1.0 + (real_T)c4_j;
    c4_s = 0.0;
    for (c4_i = 0; c4_i < 3; c4_i++) {
      c4_b_i = 1.0 + (real_T)c4_i;
      c4_b_x = c4_x[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c4_b_i), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 2, 0) - 1)) - 1];
      c4_c_x = c4_b_x;
      c4_b_y = muDoubleScalarAbs(c4_c_x);
      c4_s += c4_b_y;
    }

    c4_d_x = c4_s;
    c4_b = muDoubleScalarIsNaN(c4_d_x);
    if (c4_b) {
      c4_y = rtNaN;
      exitg1 = true;
    } else {
      if (c4_s > c4_y) {
        c4_y = c4_s;
      }

      c4_j++;
    }
  }

  return c4_y;
}

static void c4_c_eml_warning(SFc4_Model_01InstanceStruct *chartInstance)
{
  int32_T c4_i73;
  static char_T c4_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c4_u[27];
  const mxArray *c4_y = NULL;
  (void)chartInstance;
  for (c4_i73 = 0; c4_i73 < 27; c4_i73++) {
    c4_u[c4_i73] = c4_varargin_1[c4_i73];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 27), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c4_y));
}

static void c4_d_eml_warning(SFc4_Model_01InstanceStruct *chartInstance, char_T
  c4_varargin_2[14])
{
  int32_T c4_i74;
  static char_T c4_varargin_1[33] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 'i', 'l', 'l', 'C', 'o', 'n', 'd', 'i', 't', 'i',
    'o', 'n', 'e', 'd', 'M', 'a', 't', 'r', 'i', 'x' };

  char_T c4_u[33];
  const mxArray *c4_y = NULL;
  int32_T c4_i75;
  char_T c4_b_u[14];
  const mxArray *c4_b_y = NULL;
  (void)chartInstance;
  for (c4_i74 = 0; c4_i74 < 33; c4_i74++) {
    c4_u[c4_i74] = c4_varargin_1[c4_i74];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 33), false);
  for (c4_i75 = 0; c4_i75 < 14; c4_i75++) {
    c4_b_u[c4_i75] = c4_varargin_2[c4_i75];
  }

  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", c4_b_u, 10, 0U, 1U, 0U, 2, 1, 14),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c4_y, 14, c4_b_y));
}

static void c4_b_eml_scalar_eg(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c4_c_eml_scalar_eg(SFc4_Model_01InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c4_eml_xgemm(SFc4_Model_01InstanceStruct *chartInstance, real_T
  c4_A[9], real_T c4_B[9], real_T c4_C[9], real_T c4_b_C[9])
{
  int32_T c4_i76;
  int32_T c4_i77;
  real_T c4_b_A[9];
  int32_T c4_i78;
  real_T c4_b_B[9];
  for (c4_i76 = 0; c4_i76 < 9; c4_i76++) {
    c4_b_C[c4_i76] = c4_C[c4_i76];
  }

  for (c4_i77 = 0; c4_i77 < 9; c4_i77++) {
    c4_b_A[c4_i77] = c4_A[c4_i77];
  }

  for (c4_i78 = 0; c4_i78 < 9; c4_i78++) {
    c4_b_B[c4_i78] = c4_B[c4_i78];
  }

  c4_b_eml_xgemm(chartInstance, c4_b_A, c4_b_B, c4_b_C);
}

static void c4_matrix_to_scalar_power(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_a[9], real_T c4_b, real_T c4_c[9])
{
  int32_T c4_i79;
  static creal_T c4_dc5 = { 0.0, 0.0 };

  creal_T c4_V[9];
  int32_T c4_i80;
  creal_T c4_b_V[9];
  creal_T c4_beta1[3];
  creal_T c4_alpha1[3];
  real_T c4_info;
  real_T c4_b_info;
  int32_T c4_coltop;
  int32_T c4_b_coltop;
  int32_T c4_i81;
  creal_T c4_c_V[9];
  real_T c4_colnorm;
  int32_T c4_c_coltop;
  int32_T c4_b_a;
  int32_T c4_c_a;
  int32_T c4_i82;
  int32_T c4_d_a;
  int32_T c4_b_b;
  int32_T c4_e_a;
  int32_T c4_c_b;
  boolean_T c4_overflow;
  int32_T c4_j;
  int32_T c4_b_j;
  creal_T c4_r;
  real_T c4_B;
  real_T c4_y;
  real_T c4_c_info;
  real_T c4_d_info;
  int32_T c4_i83;
  creal_T c4_b_alpha1[3];
  int32_T c4_i84;
  creal_T c4_b_beta1[3];
  int32_T c4_i85;
  creal_T c4_D[9];
  int32_T c4_c_j;
  int32_T c4_d_j;
  int32_T c4_e_j;
  real_T c4_f_j;
  creal_T c4_b_D;
  int32_T c4_i;
  real_T c4_b_i;
  int32_T c4_i86;
  creal_T c4_d_V[9];
  int32_T c4_i87;
  creal_T c4_c_D[9];
  int32_T c4_i88;
  c4_b_eml_scalar_eg(chartInstance);
  c4_b_eml_error(chartInstance);
  for (c4_i79 = 0; c4_i79 < 9; c4_i79++) {
    c4_V[c4_i79].re = c4_a[c4_i79] + c4_dc5.re;
    c4_V[c4_i79].im = c4_dc5.im;
  }

  for (c4_i80 = 0; c4_i80 < 9; c4_i80++) {
    c4_b_V[c4_i80] = c4_V[c4_i80];
  }

  c4_eml_matlab_zggev(chartInstance, c4_b_V, &c4_info, c4_alpha1, c4_beta1, c4_V);
  c4_b_info = c4_info;
  for (c4_coltop = 1; c4_coltop < 8; c4_coltop += 3) {
    c4_b_coltop = c4_coltop;
    for (c4_i81 = 0; c4_i81 < 9; c4_i81++) {
      c4_c_V[c4_i81] = c4_V[c4_i81];
    }

    c4_colnorm = c4_eml_xnrm2(chartInstance, c4_c_V, c4_b_coltop);
    c4_c_coltop = c4_b_coltop;
    c4_b_a = c4_b_coltop;
    c4_c_a = c4_b_a + 2;
    c4_i82 = c4_c_a;
    c4_d_a = c4_c_coltop;
    c4_b_b = c4_i82;
    c4_e_a = c4_d_a;
    c4_c_b = c4_b_b;
    if (c4_e_a > c4_c_b) {
      c4_overflow = false;
    } else {
      c4_eml_switch_helper(chartInstance);
      c4_overflow = (c4_c_b > 2147483646);
    }

    if (c4_overflow) {
      c4_check_forloop_overflow_error(chartInstance, c4_overflow);
    }

    for (c4_j = c4_c_coltop; c4_j <= c4_i82; c4_j++) {
      c4_b_j = c4_j;
      c4_r.re = c4_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
        ("", (real_T)c4_b_j), 1, 9, 1, 0) - 1].re;
      c4_r.im = c4_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
        ("", (real_T)c4_b_j), 1, 9, 1, 0) - 1].im;
      c4_B = c4_colnorm;
      c4_y = c4_B;
      c4_r = c4_eml_div(chartInstance, c4_r, c4_y);
      c4_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c4_b_j), 1, 9, 1, 0) - 1].re = c4_r.re;
      c4_V[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c4_b_j), 1, 9, 1, 0) - 1].im = c4_r.im;
    }
  }

  c4_c_info = c4_b_info;
  c4_d_info = c4_c_info;
  for (c4_i83 = 0; c4_i83 < 3; c4_i83++) {
    c4_b_alpha1[c4_i83] = c4_alpha1[c4_i83];
  }

  for (c4_i84 = 0; c4_i84 < 3; c4_i84++) {
    c4_b_beta1[c4_i84] = c4_beta1[c4_i84];
  }

  c4_b_eml_div(chartInstance, c4_b_alpha1, c4_b_beta1, c4_alpha1);
  for (c4_i85 = 0; c4_i85 < 9; c4_i85++) {
    c4_D[c4_i85] = c4_dc5;
  }

  for (c4_c_j = 1; c4_c_j < 4; c4_c_j++) {
    c4_d_j = c4_c_j;
    c4_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_d_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
      1].re = c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 1, 0) - 1].re;
    c4_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_d_j), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
      1].im = c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 1, 0) - 1].im;
  }

  if (c4_d_info < 0.0) {
    c4_eml_warning(chartInstance);
  } else {
    if (c4_d_info > 0.0) {
      c4_b_eml_warning(chartInstance);
    }
  }

  for (c4_e_j = 0; c4_e_j < 3; c4_e_j++) {
    c4_f_j = 1.0 + (real_T)c4_e_j;
    c4_b_D.re = c4_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_f_j), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_f_j), 1, 3, 2, 0) - 1)) - 1].re;
    c4_b_D.im = c4_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_f_j), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_f_j), 1, 3, 2, 0) - 1)) - 1].im;
    c4_r = c4_power(chartInstance, c4_b_D, c4_b);
    for (c4_i = 0; c4_i < 3; c4_i++) {
      c4_b_i = 1.0 + (real_T)c4_i;
      c4_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              c4_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", c4_f_j), 1, 3, 2, 0) - 1)) - 1].re
        = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  c4_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                  (int32_T)_SFD_INTEGER_CHECK("", c4_f_j), 1, 3, 2, 0) - 1)) - 1]
        .re * c4_r.re - c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_i), 1, 3, 1, 0) + 3 *
        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", c4_f_j),
        1, 3, 2, 0) - 1)) - 1].im * c4_r.im;
      c4_D[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              c4_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", c4_f_j), 1, 3, 2, 0) - 1)) - 1].im
        = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  c4_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                  (int32_T)_SFD_INTEGER_CHECK("", c4_f_j), 1, 3, 2, 0) - 1)) - 1]
        .re * c4_r.im + c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_i), 1, 3, 1, 0) + 3 *
        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", c4_f_j),
        1, 3, 2, 0) - 1)) - 1].im * c4_r.re;
    }
  }

  for (c4_i86 = 0; c4_i86 < 9; c4_i86++) {
    c4_d_V[c4_i86] = c4_V[c4_i86];
  }

  for (c4_i87 = 0; c4_i87 < 9; c4_i87++) {
    c4_c_D[c4_i87] = c4_D[c4_i87];
  }

  c4_eml_lusolve(chartInstance, c4_d_V, c4_c_D, c4_D);
  for (c4_i88 = 0; c4_i88 < 9; c4_i88++) {
    c4_c[c4_i88] = c4_D[c4_i88].re;
  }
}

static void c4_b_eml_error(SFc4_Model_01InstanceStruct *chartInstance)
{
  int32_T c4_i89;
  static char_T c4_cv5[37] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'm', 'p', 'o', 'w', 'e', 'r', '_', 'n', 'e', 'e', 'd',
    'C', 'o', 'm', 'p', 'l', 'e', 'x', 'I', 'n', 'p', 'u', 't' };

  char_T c4_u[37];
  const mxArray *c4_y = NULL;
  (void)chartInstance;
  for (c4_i89 = 0; c4_i89 < 37; c4_i89++) {
    c4_u[c4_i89] = c4_cv5[c4_i89];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 37), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c4_y));
}

static void c4_eml_matlab_zggev(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], real_T *c4_info, creal_T c4_alpha1[3], creal_T c4_beta1[3],
  creal_T c4_V[9])
{
  int32_T c4_i90;
  creal_T c4_b_A[9];
  real_T c4_anrm;
  int32_T c4_i91;
  int32_T c4_i92;
  int32_T c4_i93;
  boolean_T c4_ilascl;
  real_T c4_anrmto;
  int32_T c4_rscale[3];
  int32_T c4_ihi;
  int32_T c4_ilo;
  int32_T c4_b_ilo;
  int32_T c4_b_ihi;
  real_T c4_b_info;
  int32_T c4_i94;
  creal_T c4_c_A[9];
  int32_T c4_c_ilo;
  int32_T c4_c_ihi;
  int32_T c4_a;
  int32_T c4_b_a;
  int32_T c4_i;
  int32_T c4_k;
  int32_T c4_j;
  int32_T c4_b_j;
  creal_T c4_tmp;
  int32_T c4_c_a;
  int32_T c4_d_a;
  int32_T c4_e_a;
  int32_T c4_f_a;
  int32_T c4_i95;
  int32_T c4_g_a;
  int32_T c4_h_a;
  boolean_T c4_overflow;
  int32_T c4_b_i;
  int32_T c4_c_j;
  int32_T c4_jc;
  real_T c4_b_jc;
  real_T c4_x;
  real_T c4_b_x;
  real_T c4_y;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_b_y;
  real_T c4_vtemp;
  int32_T c4_jr;
  real_T c4_b_jr;
  real_T c4_e_x;
  real_T c4_f_x;
  real_T c4_c_y;
  real_T c4_g_x;
  real_T c4_h_x;
  real_T c4_d_y;
  real_T c4_e_y;
  real_T c4_f_y;
  int32_T c4_c_jr;
  real_T c4_b;
  boolean_T guard1 = false;
  *c4_info = 0.0;
  c4_realmin(chartInstance);
  c4_eps(chartInstance);
  for (c4_i90 = 0; c4_i90 < 9; c4_i90++) {
    c4_b_A[c4_i90] = c4_A[c4_i90];
  }

  c4_anrm = c4_eml_matlab_zlangeM(chartInstance, c4_b_A);
  if (!c4_isfinite(chartInstance, c4_anrm)) {
    for (c4_i91 = 0; c4_i91 < 3; c4_i91++) {
      c4_alpha1[c4_i91].re = rtNaN;
      c4_alpha1[c4_i91].im = 0.0;
    }

    for (c4_i92 = 0; c4_i92 < 3; c4_i92++) {
      c4_beta1[c4_i92].re = rtNaN;
      c4_beta1[c4_i92].im = 0.0;
    }

    for (c4_i93 = 0; c4_i93 < 9; c4_i93++) {
      c4_V[c4_i93].re = rtNaN;
      c4_V[c4_i93].im = 0.0;
    }
  } else {
    c4_ilascl = false;
    c4_anrmto = c4_anrm;
    guard1 = false;
    if (c4_anrm > 0.0) {
      if (c4_anrm < 6.7178761075670888E-139) {
        c4_anrmto = 6.7178761075670888E-139;
        c4_ilascl = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      if (c4_anrm > 1.4885657073574029E+138) {
        c4_anrmto = 1.4885657073574029E+138;
        c4_ilascl = true;
      }
    }

    if (c4_ilascl) {
      c4_c_eml_matlab_zlascl(chartInstance, c4_anrm, c4_anrmto, c4_A);
    }

    c4_b_eml_matlab_zggbal(chartInstance, c4_A, &c4_ilo, &c4_ihi, c4_rscale);
    c4_b_ilo = c4_ilo;
    c4_b_ihi = c4_ihi;
    c4_b_eml_matlab_zgghrd(chartInstance, c4_b_ilo, c4_b_ihi, c4_A, c4_V);
    c4_c_eml_matlab_zhgeqz(chartInstance, c4_A, c4_b_ilo, c4_b_ihi, c4_V,
      &c4_b_info, c4_alpha1, c4_beta1);
    *c4_info = c4_b_info;
    if (*c4_info != 0.0) {
    } else {
      for (c4_i94 = 0; c4_i94 < 9; c4_i94++) {
        c4_c_A[c4_i94] = c4_A[c4_i94];
      }

      c4_b_eml_matlab_ztgevc(chartInstance, c4_c_A, c4_V);
      c4_c_ilo = c4_b_ilo;
      c4_c_ihi = c4_b_ihi;
      if (c4_c_ilo > 1) {
        c4_a = c4_c_ilo;
        c4_b_a = c4_a - 1;
        c4_i = c4_b_a;
        while (c4_i >= 1) {
          c4_k = c4_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_i), 1, 3, 1, 0) - 1];
          if (c4_k != c4_i) {
            for (c4_j = 1; c4_j < 4; c4_j++) {
              c4_b_j = c4_j;
              c4_tmp.re = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].
                re;
              c4_tmp.im = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].
                im;
              c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re = c4_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
              c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im = c4_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
              c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re =
                c4_tmp.re;
              c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im =
                c4_tmp.im;
            }
          }

          c4_c_a = c4_i;
          c4_d_a = c4_c_a - 1;
          c4_i = c4_d_a;
        }
      }

      if (c4_c_ihi < 3) {
        c4_e_a = c4_c_ihi;
        c4_f_a = c4_e_a + 1;
        c4_i95 = c4_f_a;
        c4_g_a = c4_i95;
        c4_h_a = c4_g_a;
        if (c4_h_a > 3) {
          c4_overflow = false;
        } else {
          c4_eml_switch_helper(chartInstance);
          c4_overflow = false;
        }

        if (c4_overflow) {
          c4_check_forloop_overflow_error(chartInstance, c4_overflow);
        }

        for (c4_b_i = c4_i95; c4_b_i < 4; c4_b_i++) {
          c4_i = c4_b_i;
          c4_k = c4_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_i), 1, 3, 1, 0) - 1];
          if (c4_k != c4_i) {
            for (c4_c_j = 1; c4_c_j < 4; c4_c_j++) {
              c4_b_j = c4_c_j;
              c4_tmp.re = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].
                re;
              c4_tmp.im = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].
                im;
              c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re = c4_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
              c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im = c4_V
                [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
              c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re =
                c4_tmp.re;
              c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_k), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im =
                c4_tmp.im;
            }
          }
        }
      }

      for (c4_jc = 0; c4_jc < 3; c4_jc++) {
        c4_b_jc = 1.0 + (real_T)c4_jc;
        c4_tmp.re = c4_V[3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1)].re;
        c4_tmp.im = c4_V[3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1)].im;
        c4_x = c4_tmp.re;
        c4_b_x = c4_x;
        c4_y = muDoubleScalarAbs(c4_b_x);
        c4_c_x = c4_tmp.im;
        c4_d_x = c4_c_x;
        c4_b_y = muDoubleScalarAbs(c4_d_x);
        c4_vtemp = c4_y + c4_b_y;
        for (c4_jr = 0; c4_jr < 2; c4_jr++) {
          c4_b_jr = 2.0 + (real_T)c4_jr;
          c4_tmp.re = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1)) - 1].re;
          c4_tmp.im = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1)) - 1].im;
          c4_e_x = c4_tmp.re;
          c4_f_x = c4_e_x;
          c4_c_y = muDoubleScalarAbs(c4_f_x);
          c4_g_x = c4_tmp.im;
          c4_h_x = c4_g_x;
          c4_d_y = muDoubleScalarAbs(c4_h_x);
          c4_e_y = c4_c_y + c4_d_y;
          c4_f_y = c4_e_y;
          if (c4_f_y > c4_vtemp) {
            c4_vtemp = c4_f_y;
          }
        }

        if (c4_vtemp >= 6.7178761075670888E-139) {
          c4_vtemp = 1.0 / c4_vtemp;
          for (c4_c_jr = 0; c4_c_jr < 3; c4_c_jr++) {
            c4_b_jr = 1.0 + (real_T)c4_c_jr;
            c4_tmp.re = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1)) - 1].re;
            c4_tmp.im = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1)) - 1].im;
            c4_b = c4_vtemp;
            c4_tmp.re *= c4_b;
            c4_tmp.im *= c4_b;
            c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    c4_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1))
              - 1].re = c4_tmp.re;
            c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    c4_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1))
              - 1].im = c4_tmp.im;
          }
        }
      }

      if (c4_ilascl) {
        c4_d_eml_matlab_zlascl(chartInstance, c4_anrmto, c4_anrm, c4_alpha1);
      }
    }
  }
}

static void c4_eml_matlab_zgghrd(SFc4_Model_01InstanceStruct *chartInstance,
  int32_T c4_ilo, int32_T c4_ihi, creal_T c4_A[9], creal_T c4_b_A[9], creal_T
  c4_Z[9])
{
  int32_T c4_i96;
  for (c4_i96 = 0; c4_i96 < 9; c4_i96++) {
    c4_b_A[c4_i96] = c4_A[c4_i96];
  }

  c4_b_eml_matlab_zgghrd(chartInstance, c4_ilo, c4_ihi, c4_b_A, c4_Z);
}

static void c4_b_eml_matlab_zhgeqz(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T c4_ilo, int32_T c4_ihi, creal_T c4_Z[9], real_T
  *c4_info, creal_T c4_alpha1[3], creal_T c4_beta1[3], creal_T c4_b_A[9],
  creal_T c4_b_Z[9])
{
  int32_T c4_i97;
  int32_T c4_i98;
  for (c4_i97 = 0; c4_i97 < 9; c4_i97++) {
    c4_b_A[c4_i97] = c4_A[c4_i97];
  }

  for (c4_i98 = 0; c4_i98 < 9; c4_i98++) {
    c4_b_Z[c4_i98] = c4_Z[c4_i98];
  }

  c4_c_eml_matlab_zhgeqz(chartInstance, c4_b_A, c4_ilo, c4_ihi, c4_b_Z, c4_info,
    c4_alpha1, c4_beta1);
}

static void c4_eml_matlab_ztgevc(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], creal_T c4_V[9], creal_T c4_b_V[9])
{
  int32_T c4_i99;
  int32_T c4_i100;
  creal_T c4_b_A[9];
  for (c4_i99 = 0; c4_i99 < 9; c4_i99++) {
    c4_b_V[c4_i99] = c4_V[c4_i99];
  }

  for (c4_i100 = 0; c4_i100 < 9; c4_i100++) {
    c4_b_A[c4_i100] = c4_A[c4_i100];
  }

  c4_b_eml_matlab_ztgevc(chartInstance, c4_b_A, c4_b_V);
}

static creal_T c4_rdivide(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_x, creal_T c4_y)
{
  creal_T c4_z;
  real_T c4_ar;
  real_T c4_ai;
  real_T c4_br;
  real_T c4_bi;
  real_T c4_brm;
  real_T c4_bim;
  real_T c4_s;
  real_T c4_d;
  real_T c4_nr;
  real_T c4_ni;
  real_T c4_sgnbr;
  real_T c4_sgnbi;
  (void)chartInstance;
  c4_ar = c4_x.re;
  c4_ai = c4_x.im;
  c4_br = c4_y.re;
  c4_bi = c4_y.im;
  if (c4_bi == 0.0) {
    if (c4_ai == 0.0) {
      c4_z.re = c4_ar / c4_br;
      c4_z.im = 0.0;
    } else if (c4_ar == 0.0) {
      c4_z.re = 0.0;
      c4_z.im = c4_ai / c4_br;
    } else {
      c4_z.re = c4_ar / c4_br;
      c4_z.im = c4_ai / c4_br;
    }
  } else if (c4_br == 0.0) {
    if (c4_ar == 0.0) {
      c4_z.re = c4_ai / c4_bi;
      c4_z.im = 0.0;
    } else if (c4_ai == 0.0) {
      c4_z.re = 0.0;
      c4_z.im = -(c4_ar / c4_bi);
    } else {
      c4_z.re = c4_ai / c4_bi;
      c4_z.im = -(c4_ar / c4_bi);
    }
  } else {
    c4_brm = muDoubleScalarAbs(c4_br);
    c4_bim = muDoubleScalarAbs(c4_bi);
    if (c4_brm > c4_bim) {
      c4_s = c4_bi / c4_br;
      c4_d = c4_br + c4_s * c4_bi;
      c4_nr = c4_ar + c4_s * c4_ai;
      c4_ni = c4_ai - c4_s * c4_ar;
      c4_z.re = c4_nr / c4_d;
      c4_z.im = c4_ni / c4_d;
    } else if (c4_bim == c4_brm) {
      if (c4_br > 0.0) {
        c4_sgnbr = 0.5;
      } else {
        c4_sgnbr = -0.5;
      }

      if (c4_bi > 0.0) {
        c4_sgnbi = 0.5;
      } else {
        c4_sgnbi = -0.5;
      }

      c4_nr = c4_ar * c4_sgnbr + c4_ai * c4_sgnbi;
      c4_ni = c4_ai * c4_sgnbr - c4_ar * c4_sgnbi;
      c4_z.re = c4_nr / c4_brm;
      c4_z.im = c4_ni / c4_brm;
    } else {
      c4_s = c4_br / c4_bi;
      c4_d = c4_bi + c4_s * c4_br;
      c4_nr = c4_s * c4_ar + c4_ai;
      c4_ni = c4_s * c4_ai - c4_ar;
      c4_z.re = c4_nr / c4_d;
      c4_z.im = c4_ni / c4_d;
    }
  }

  return c4_z;
}

static real_T c4_eml_xnrm2(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_x[9], int32_T c4_ix0)
{
  real_T c4_y;
  int32_T c4_b_ix0;
  int32_T c4_c_ix0;
  real_T c4_scale;
  int32_T c4_kstart;
  int32_T c4_a;
  int32_T c4_kend;
  int32_T c4_b_kstart;
  int32_T c4_b_kend;
  int32_T c4_b_a;
  int32_T c4_b;
  int32_T c4_c_a;
  int32_T c4_b_b;
  boolean_T c4_overflow;
  int32_T c4_k;
  int32_T c4_b_k;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_absxk;
  real_T c4_t;
  real_T c4_d_x;
  real_T c4_e_x;
  c4_b_ix0 = c4_ix0;
  c4_c_ix0 = c4_b_ix0;
  c4_y = 0.0;
  c4_realmin(chartInstance);
  c4_scale = 2.2250738585072014E-308;
  c4_kstart = c4_c_ix0;
  c4_a = c4_kstart;
  c4_kend = c4_a;
  c4_b_kstart = c4_kstart;
  c4_b_kend = c4_kend + 2;
  c4_b_a = c4_b_kstart;
  c4_b = c4_b_kend;
  c4_c_a = c4_b_a;
  c4_b_b = c4_b;
  if (c4_c_a > c4_b_b) {
    c4_overflow = false;
  } else {
    c4_eml_switch_helper(chartInstance);
    c4_overflow = (c4_b_b > 2147483646);
  }

  if (c4_overflow) {
    c4_check_forloop_overflow_error(chartInstance, c4_overflow);
  }

  for (c4_k = c4_b_kstart; c4_k <= c4_b_kend; c4_k++) {
    c4_b_k = c4_k;
    c4_b_x = c4_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_b_k), 1, 9, 1, 0) - 1].re;
    c4_c_x = c4_b_x;
    c4_absxk = muDoubleScalarAbs(c4_c_x);
    if (c4_absxk > c4_scale) {
      c4_t = c4_scale / c4_absxk;
      c4_y = 1.0 + c4_y * c4_t * c4_t;
      c4_scale = c4_absxk;
    } else {
      c4_t = c4_absxk / c4_scale;
      c4_y += c4_t * c4_t;
    }

    c4_d_x = c4_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_b_k), 1, 9, 1, 0) - 1].im;
    c4_e_x = c4_d_x;
    c4_absxk = muDoubleScalarAbs(c4_e_x);
    if (c4_absxk > c4_scale) {
      c4_t = c4_scale / c4_absxk;
      c4_y = 1.0 + c4_y * c4_t * c4_t;
      c4_scale = c4_absxk;
    } else {
      c4_t = c4_absxk / c4_scale;
      c4_y += c4_t * c4_t;
    }
  }

  return c4_scale * muDoubleScalarSqrt(c4_y);
}

static creal_T c4_power(SFc4_Model_01InstanceStruct *chartInstance, creal_T c4_a,
  real_T c4_b)
{
  creal_T c4_y;
  real_T c4_b_b;
  real_T c4_bk;
  real_T c4_c_b;
  real_T c4_ar;
  real_T c4_ai;
  real_T c4_br;
  real_T c4_bi;
  real_T c4_ytmp;
  real_T c4_x;
  real_T c4_xk;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_r;
  real_T c4_d6;
  int8_T c4_i101;
  int8_T c4_b_r;
  creal_T c4_t;
  real_T c4_e_x;
  boolean_T c4_d_b;
  real_T c4_A;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_h_x;
  real_T c4_b_y;
  real_T c4_b_A;
  real_T c4_i_x;
  real_T c4_j_x;
  real_T c4_k_x;
  real_T c4_c_y;
  real_T c4_l_x;
  real_T c4_d_y;
  real_T c4_x1;
  real_T c4_x2;
  real_T c4_b_a;
  real_T c4_e_b;
  real_T c4_z;
  real_T c4_e_y;
  real_T c4_m_x;
  real_T c4_c_r;
  real_T c4_b_x1;
  real_T c4_b_x2;
  real_T c4_c_a;
  real_T c4_f_b;
  real_T c4_f_y;
  real_T c4_g_y;
  real_T c4_n_x;
  real_T c4_d_r;
  real_T c4_d_a;
  real_T c4_tr;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T guard5 = false;
  boolean_T guard6 = false;
  boolean_T guard7 = false;
  c4_b_b = c4_b;
  c4_scalarEg(chartInstance);
  c4_bk = c4_b_b;
  c4_c_b = c4_bk;
  c4_scalarEg(chartInstance);
  c4_ar = c4_a.re;
  c4_ai = c4_a.im;
  c4_br = c4_c_b;
  c4_bi = 0.0;
  guard1 = false;
  guard2 = false;
  if (c4_ai == 0.0) {
    if (c4_bi == 0.0) {
      if (c4_ar >= 0.0) {
        c4_y.re = muDoubleScalarPower(c4_ar, c4_br);
        c4_y.im = 0.0;
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
    if (c4_ar == 0.0) {
      if (c4_bi == 0.0) {
        if (muDoubleScalarFloor(c4_br) == c4_br) {
          if (muDoubleScalarAbs(c4_br) <= 9.007199254740992E+15) {
            c4_ytmp = muDoubleScalarPower(c4_ai, c4_br);
            c4_x = c4_br;
            c4_eml_scalar_eg(chartInstance);
            c4_xk = c4_x;
            c4_b_x = c4_xk;
            c4_eml_scalar_eg(chartInstance);
            c4_c_x = c4_b_x / 4.0;
            c4_d_x = c4_c_x;
            c4_d_x = muDoubleScalarFloor(c4_d_x);
            c4_r = c4_b_x - c4_d_x * 4.0;
            c4_d6 = muDoubleScalarRound(c4_r);
            if (c4_d6 < 128.0) {
              if (c4_d6 >= -128.0) {
                c4_i101 = (int8_T)c4_d6;
              } else {
                c4_i101 = MIN_int8_T;
              }
            } else if (c4_d6 >= 128.0) {
              c4_i101 = MAX_int8_T;
            } else {
              c4_i101 = 0;
            }

            c4_b_r = c4_i101;
            if ((real_T)c4_b_r == 3.0) {
              c4_y.re = 0.0;
              c4_y.im = -c4_ytmp;
            } else if ((real_T)c4_b_r == 2.0) {
              c4_y.re = -c4_ytmp;
              c4_y.im = 0.0;
            } else if ((real_T)c4_b_r == 1.0) {
              c4_y.re = 0.0;
              c4_y.im = c4_ytmp;
            } else {
              c4_y.re = c4_ytmp;
              c4_y.im = 0.0;
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
      c4_t = c4_a;
      c4_realmax(chartInstance);
      guard6 = false;
      if (c4_t.im == 0.0) {
        c4_e_x = c4_t.re;
        c4_d_b = muDoubleScalarIsNaN(c4_e_x);
        if (c4_d_b) {
        } else {
          guard6 = true;
        }
      } else {
        guard6 = true;
      }

      if (guard6 == true) {
        guard7 = false;
        if (muDoubleScalarAbs(c4_t.re) > 8.9884656743115785E+307) {
          guard7 = true;
        } else if (muDoubleScalarAbs(c4_t.im) > 8.9884656743115785E+307) {
          guard7 = true;
        } else {
          c4_b_x1 = c4_t.re;
          c4_b_x2 = c4_t.im;
          c4_c_a = c4_b_x1;
          c4_f_b = c4_b_x2;
          c4_f_y = muDoubleScalarHypot(c4_c_a, c4_f_b);
          c4_g_y = c4_t.im;
          c4_n_x = c4_t.re;
          c4_d_r = muDoubleScalarAtan2(c4_g_y, c4_n_x);
          c4_t.re = muDoubleScalarLog(c4_f_y);
          c4_t.im = c4_d_r;
        }

        if (guard7 == true) {
          c4_A = c4_t.re;
          c4_f_x = c4_A;
          c4_g_x = c4_f_x;
          c4_h_x = c4_g_x;
          c4_b_y = c4_h_x / 2.0;
          c4_b_A = c4_t.im;
          c4_i_x = c4_b_A;
          c4_j_x = c4_i_x;
          c4_k_x = c4_j_x;
          c4_c_y = c4_k_x / 2.0;
          c4_l_x = c4_b_y;
          c4_d_y = c4_c_y;
          c4_eml_scalar_eg(chartInstance);
          c4_x1 = c4_l_x;
          c4_x2 = c4_d_y;
          c4_b_a = c4_x1;
          c4_e_b = c4_x2;
          c4_z = muDoubleScalarHypot(c4_b_a, c4_e_b);
          c4_e_y = c4_t.im;
          c4_m_x = c4_t.re;
          c4_c_r = muDoubleScalarAtan2(c4_e_y, c4_m_x);
          c4_t.re = muDoubleScalarLog(c4_z) + 0.69314718055994529;
          c4_t.im = c4_c_r;
        }
      }

      c4_d_a = c4_c_b;
      c4_t.re *= c4_d_a;
      c4_t.im *= c4_d_a;
      c4_tr = muDoubleScalarExp(c4_t.re);
      c4_y.re = c4_tr * muDoubleScalarCos(c4_t.im);
      c4_y.im = c4_tr * muDoubleScalarSin(c4_t.im);
    }
  }

  return c4_y;
}

static void c4_eml_lusolve(SFc4_Model_01InstanceStruct *chartInstance, creal_T
  c4_A[9], creal_T c4_B[9], creal_T c4_X[9])
{
  int32_T c4_i102;
  creal_T c4_b_A[9];
  int32_T c4_r1;
  int32_T c4_r2;
  int32_T c4_r3;
  creal_T c4_x;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_y;
  real_T c4_d_x;
  real_T c4_e_x;
  real_T c4_b_y;
  real_T c4_maxval;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_c_y;
  real_T c4_h_x;
  real_T c4_i_x;
  real_T c4_d_y;
  real_T c4_a21;
  real_T c4_j_x;
  real_T c4_k_x;
  real_T c4_e_y;
  real_T c4_l_x;
  real_T c4_m_x;
  real_T c4_f_y;
  real_T c4_d;
  creal_T c4_c_A;
  creal_T c4_d_A;
  creal_T c4_e_A;
  creal_T c4_f_A;
  creal_T c4_g_A;
  creal_T c4_h_A;
  creal_T c4_i_A;
  creal_T c4_j_A;
  real_T c4_n_x;
  real_T c4_o_x;
  real_T c4_g_y;
  real_T c4_p_x;
  real_T c4_q_x;
  real_T c4_h_y;
  real_T c4_b_d;
  real_T c4_r_x;
  real_T c4_s_x;
  real_T c4_i_y;
  real_T c4_t_x;
  real_T c4_u_x;
  real_T c4_j_y;
  real_T c4_c_d;
  int32_T c4_rtemp;
  creal_T c4_k_A;
  creal_T c4_l_A;
  creal_T c4_m_A;
  static creal_T c4_dc6 = { 0.0, 0.0 };

  boolean_T c4_n_A;
  boolean_T c4_o_A;
  boolean_T c4_p_A;
  int32_T c4_k;
  int32_T c4_b_k;
  creal_T c4_b_B;
  creal_T c4_q_A;
  creal_T c4_b_X;
  creal_T c4_c_X;
  creal_T c4_d_X;
  creal_T c4_r_A;
  creal_T c4_e_X;
  creal_T c4_f_X;
  creal_T c4_s_A;
  creal_T c4_g_X;
  creal_T c4_h_X;
  creal_T c4_i_X;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  for (c4_i102 = 0; c4_i102 < 9; c4_i102++) {
    c4_b_A[c4_i102] = c4_A[c4_i102];
  }

  c4_r1 = 1;
  c4_r2 = 2;
  c4_r3 = 3;
  c4_x = c4_b_A[0];
  c4_b_x = c4_x.re;
  c4_c_x = c4_b_x;
  c4_y = muDoubleScalarAbs(c4_c_x);
  c4_d_x = c4_x.im;
  c4_e_x = c4_d_x;
  c4_b_y = muDoubleScalarAbs(c4_e_x);
  c4_maxval = c4_y + c4_b_y;
  c4_x = c4_b_A[1];
  c4_f_x = c4_x.re;
  c4_g_x = c4_f_x;
  c4_c_y = muDoubleScalarAbs(c4_g_x);
  c4_h_x = c4_x.im;
  c4_i_x = c4_h_x;
  c4_d_y = muDoubleScalarAbs(c4_i_x);
  c4_a21 = c4_c_y + c4_d_y;
  if (c4_a21 > c4_maxval) {
    c4_maxval = c4_a21;
    c4_r1 = 2;
    c4_r2 = 1;
  }

  c4_x = c4_b_A[2];
  c4_j_x = c4_x.re;
  c4_k_x = c4_j_x;
  c4_e_y = muDoubleScalarAbs(c4_k_x);
  c4_l_x = c4_x.im;
  c4_m_x = c4_l_x;
  c4_f_y = muDoubleScalarAbs(c4_m_x);
  c4_d = c4_e_y + c4_f_y;
  if (c4_d > c4_maxval) {
    c4_r1 = 3;
    c4_r2 = 2;
    c4_r3 = 1;
  }

  c4_c_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r2), 1, 3, 1, 0) - 1].re;
  c4_c_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r2), 1, 3, 1, 0) - 1].im;
  c4_d_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r1), 1, 3, 1, 0) - 1].re;
  c4_d_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r1), 1, 3, 1, 0) - 1].im;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r2), 1, 3, 1, 0) - 1] = c4_rdivide(chartInstance, c4_c_A, c4_d_A);
  c4_e_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) - 1].re;
  c4_e_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) - 1].im;
  c4_f_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r1), 1, 3, 1, 0) - 1].re;
  c4_f_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r1), 1, 3, 1, 0) - 1].im;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r3), 1, 3, 1, 0) - 1] = c4_rdivide(chartInstance, c4_e_A, c4_f_A);
  c4_g_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r2), 1, 3, 1, 0) - 1].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 2].re - c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) - 1].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 2].im;
  c4_g_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r2), 1, 3, 1, 0) - 1].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 2].im + c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) - 1].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 2].re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r2), 1, 3, 1, 0) + 2].re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 2].re -
    c4_g_A.re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r2), 1, 3, 1, 0) + 2].im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 2].im -
    c4_g_A.im;
  c4_h_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) - 1].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 2].re - c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) - 1].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 2].im;
  c4_h_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) - 1].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 2].im + c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) - 1].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 2].re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r3), 1, 3, 1, 0) + 2].re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 2].re -
    c4_h_A.re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r3), 1, 3, 1, 0) + 2].im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 2].im -
    c4_h_A.im;
  c4_i_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r2), 1, 3, 1, 0) - 1].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 5].re - c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) - 1].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 5].im;
  c4_i_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r2), 1, 3, 1, 0) - 1].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 5].im + c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) - 1].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 5].re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r2), 1, 3, 1, 0) + 5].re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 5].re -
    c4_i_A.re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r2), 1, 3, 1, 0) + 5].im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 5].im -
    c4_i_A.im;
  c4_j_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) - 1].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 5].re - c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) - 1].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 5].im;
  c4_j_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) - 1].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 5].im + c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) - 1].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r1), 1, 3, 1, 0) + 5].re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r3), 1, 3, 1, 0) + 5].re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 5].re -
    c4_j_A.re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r3), 1, 3, 1, 0) + 5].im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 5].im -
    c4_j_A.im;
  c4_x.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c4_r3), 1, 3, 1, 0) + 2].re;
  c4_x.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c4_r3), 1, 3, 1, 0) + 2].im;
  c4_n_x = c4_x.re;
  c4_o_x = c4_n_x;
  c4_g_y = muDoubleScalarAbs(c4_o_x);
  c4_p_x = c4_x.im;
  c4_q_x = c4_p_x;
  c4_h_y = muDoubleScalarAbs(c4_q_x);
  c4_b_d = c4_g_y + c4_h_y;
  c4_x.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c4_r2), 1, 3, 1, 0) + 2].re;
  c4_x.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c4_r2), 1, 3, 1, 0) + 2].im;
  c4_r_x = c4_x.re;
  c4_s_x = c4_r_x;
  c4_i_y = muDoubleScalarAbs(c4_s_x);
  c4_t_x = c4_x.im;
  c4_u_x = c4_t_x;
  c4_j_y = muDoubleScalarAbs(c4_u_x);
  c4_c_d = c4_i_y + c4_j_y;
  if (c4_b_d > c4_c_d) {
    c4_rtemp = c4_r2;
    c4_r2 = c4_r3;
    c4_r3 = c4_rtemp;
  }

  c4_k_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) + 2].re;
  c4_k_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) + 2].im;
  c4_l_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r2), 1, 3, 1, 0) + 2].re;
  c4_l_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r2), 1, 3, 1, 0) + 2].im;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r3), 1, 3, 1, 0) + 2] = c4_rdivide(chartInstance, c4_k_A, c4_l_A);
  c4_m_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) + 2].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r2), 1, 3, 1, 0) + 5].re - c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 2].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r2), 1, 3, 1, 0) + 5].im;
  c4_m_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
    ("", (real_T)c4_r3), 1, 3, 1, 0) + 2].re *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r2), 1, 3, 1, 0) + 5].im + c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 2].im *
    c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c4_r2), 1, 3, 1, 0) + 5].re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r3), 1, 3, 1, 0) + 5].re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 5].re -
    c4_m_A.re;
  c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c4_r3), 1, 3, 1, 0) + 5].im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 5].im -
    c4_m_A.im;
  c4_n_A = ((c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c4_r1), 1, 3, 1, 0) - 1].re == c4_dc6.re) &&
            (c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c4_r1), 1, 3, 1, 0) - 1].im == c4_dc6.im));
  guard1 = false;
  guard2 = false;
  if (c4_n_A) {
    guard2 = true;
  } else {
    c4_o_A = ((c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 2].re ==
               c4_dc6.re) && (c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 2].im ==
               c4_dc6.im));
    if (c4_o_A) {
      guard2 = true;
    } else {
      c4_p_A = ((c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 5].re ==
                 c4_dc6.re) && (c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 5].im ==
                 c4_dc6.im));
      if (c4_p_A) {
        guard1 = true;
      }
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c4_c_eml_warning(chartInstance);
  }

  for (c4_k = 1; c4_k < 4; c4_k++) {
    c4_b_k = c4_k;
    c4_b_B.re = c4_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", (real_T)c4_b_k), 1, 3, 1, 0) - 1].re;
    c4_b_B.im = c4_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", (real_T)c4_b_k), 1, 3, 1, 0) - 1].im;
    c4_q_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 1, 0) - 1].re;
    c4_q_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 1, 0) - 1].im;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) -
      1] = c4_rdivide(chartInstance, c4_b_B, c4_q_A);
    c4_b_X.re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r1), 1, 3, 1, 0) + 2].re - c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r1), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 1, 0) + 2].im;
    c4_b_X.im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r1), 1, 3, 1, 0) + 2].im + c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r1), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 1, 0) + 2].re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) -
      1].re = c4_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 2].re - c4_b_X.re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) -
      1].im = c4_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 2].im - c4_b_X.im;
    c4_c_X.re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r1), 1, 3, 1, 0) + 5].re - c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r1), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 1, 0) + 5].im;
    c4_c_X.im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r1), 1, 3, 1, 0) + 5].im + c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r1), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 1, 0) + 5].re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) -
      1].re = c4_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 5].re - c4_c_X.re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) -
      1].im = c4_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 5].im - c4_c_X.im;
    c4_d_X.re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) - 1].re;
    c4_d_X.im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) - 1].im;
    c4_r_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 2].re;
    c4_r_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 2].im;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) -
      1] = c4_rdivide(chartInstance, c4_d_X, c4_r_A);
    c4_e_X.re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r2), 1, 3, 1, 0) + 5].re - c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r2), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 5].im;
    c4_e_X.im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r2), 1, 3, 1, 0) + 5].im + c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r2), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) + 5].re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) -
      1].re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) - 1].re
      - c4_e_X.re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) -
      1].im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) - 1].im
      - c4_e_X.im;
    c4_f_X.re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) - 1].re;
    c4_f_X.im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) - 1].im;
    c4_s_A.re = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 5].re;
    c4_s_A.im = c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 5].im;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) -
      1] = c4_rdivide(chartInstance, c4_f_X, c4_s_A);
    c4_g_X.re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r3), 1, 3, 1, 0) + 2].re - c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r3), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 2].im;
    c4_g_X.im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r3), 1, 3, 1, 0) + 2].im + c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r3), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) + 2].re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) -
      1].re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) - 1].re
      - c4_g_X.re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) -
      1].im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) - 1].im
      - c4_g_X.im;
    c4_h_X.re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r3), 1, 3, 1, 0) - 1].re - c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r3), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) - 1].im;
    c4_h_X.im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r3), 1, 3, 1, 0) - 1].im + c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r3), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r3), 1, 3, 1, 0) - 1].re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) -
      1].re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) - 1].re
      - c4_h_X.re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) -
      1].im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) - 1].im
      - c4_h_X.im;
    c4_i_X.re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r2), 1, 3, 1, 0) - 1].re - c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r2), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) - 1].im;
    c4_i_X.im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 2, 0) - 1)) - 1].re *
      c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_r2), 1, 3, 1, 0) - 1].im + c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_r2), 1, 3, 2, 0) - 1)) - 1].im * c4_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r2), 1, 3, 1, 0) - 1].re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) -
      1].re = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) - 1].re
      - c4_i_X.re;
    c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) -
      1].im = c4_X[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c4_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_r1), 1, 3, 2, 0) - 1)) - 1].im
      - c4_i_X.im;
  }
}

static void c4_e_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_sprintf, const char_T *c4_identifier, char_T c4_y[14])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_sprintf), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_sprintf);
}

static void c4_f_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, char_T c4_y[14])
{
  char_T c4_cv6[14];
  int32_T c4_i103;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_cv6, 1, 10, 0U, 1, 0U, 2, 1,
                14);
  for (c4_i103 = 0; c4_i103 < 14; c4_i103++) {
    c4_y[c4_i103] = c4_cv6[c4_i103];
  }

  sf_mex_destroy(&c4_u);
}

static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_u;
  const mxArray *c4_y = NULL;
  SFc4_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc4_Model_01InstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_u = *(int32_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, false);
  return c4_mxArrayOutData;
}

static int32_T c4_g_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  int32_T c4_y;
  int32_T c4_i104;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_i104, 1, 6, 0U, 0, 0U, 0);
  c4_y = c4_i104;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_b_sfEvent;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  int32_T c4_y;
  SFc4_Model_01InstanceStruct *chartInstance;
  chartInstance = (SFc4_Model_01InstanceStruct *)chartInstanceVoid;
  c4_b_sfEvent = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_sfEvent),
    &c4_thisId);
  sf_mex_destroy(&c4_b_sfEvent);
  *(int32_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static uint8_T c4_h_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_b_is_active_c4_Model_01, const char_T *c4_identifier)
{
  uint8_T c4_y;
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_i_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c4_b_is_active_c4_Model_01), &c4_thisId);
  sf_mex_destroy(&c4_b_is_active_c4_Model_01);
  return c4_y;
}

static uint8_T c4_i_emlrt_marshallIn(SFc4_Model_01InstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  uint8_T c4_y;
  uint8_T c4_u0;
  (void)chartInstance;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_u0, 1, 3, 0U, 0, 0U, 0);
  c4_y = c4_u0;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_c_eml_matlab_zlascl(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_cfrom, real_T c4_cto, creal_T c4_A[9])
{
  real_T c4_cfromc;
  real_T c4_ctoc;
  boolean_T c4_notdone;
  real_T c4_cfrom1;
  real_T c4_cto1;
  real_T c4_x;
  real_T c4_b_x;
  real_T c4_y;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_b_y;
  real_T c4_mul;
  real_T c4_e_x;
  real_T c4_f_x;
  real_T c4_c_y;
  real_T c4_g_x;
  real_T c4_h_x;
  real_T c4_d_y;
  real_T c4_a;
  int32_T c4_i105;
  int32_T c4_i106;
  int32_T c4_i107;
  boolean_T guard1 = false;
  c4_realmin(chartInstance);
  c4_eps(chartInstance);
  c4_cfromc = c4_cfrom;
  c4_ctoc = c4_cto;
  c4_notdone = true;
  while (c4_notdone) {
    c4_cfrom1 = c4_cfromc * 2.0041683600089728E-292;
    c4_cto1 = c4_ctoc / 4.9896007738368E+291;
    c4_x = c4_cfrom1;
    c4_b_x = c4_x;
    c4_y = muDoubleScalarAbs(c4_b_x);
    c4_c_x = c4_ctoc;
    c4_d_x = c4_c_x;
    c4_b_y = muDoubleScalarAbs(c4_d_x);
    guard1 = false;
    if (c4_y > c4_b_y) {
      if (c4_ctoc != 0.0) {
        c4_mul = 2.0041683600089728E-292;
        c4_notdone = true;
        c4_cfromc = c4_cfrom1;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      c4_e_x = c4_cto1;
      c4_f_x = c4_e_x;
      c4_c_y = muDoubleScalarAbs(c4_f_x);
      c4_g_x = c4_cfromc;
      c4_h_x = c4_g_x;
      c4_d_y = muDoubleScalarAbs(c4_h_x);
      if (c4_c_y > c4_d_y) {
        c4_mul = 4.9896007738368E+291;
        c4_notdone = true;
        c4_ctoc = c4_cto1;
      } else {
        c4_mul = c4_ctoc / c4_cfromc;
        c4_notdone = false;
      }
    }

    c4_a = c4_mul;
    c4_i105 = 0;
    for (c4_i106 = 0; c4_i106 < 3; c4_i106++) {
      for (c4_i107 = 0; c4_i107 < 3; c4_i107++) {
        c4_A[c4_i107 + c4_i105].re *= c4_a;
        c4_A[c4_i107 + c4_i105].im *= c4_a;
      }

      c4_i105 += 3;
    }
  }
}

static void c4_b_eml_matlab_zggbal(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T *c4_ilo, int32_T *c4_ihi, int32_T c4_rscale[3])
{
  int32_T c4_i108;
  int32_T c4_b_ihi;
  int32_T c4_i;
  int32_T c4_j;
  boolean_T c4_found;
  int32_T c4_ii;
  real_T c4_nzcount;
  int32_T c4_c_ihi;
  int32_T c4_b;
  int32_T c4_b_b;
  boolean_T c4_overflow;
  int32_T c4_jj;
  int32_T c4_b_jj;
  static creal_T c4_dc7 = { 0.0, 0.0 };

  boolean_T c4_b_A;
  int32_T c4_a;
  int32_T c4_b_a;
  int32_T c4_b_i;
  int32_T c4_b_j;
  boolean_T c4_b_found;
  int32_T c4_b_ilo;
  int32_T c4_d_ihi;
  int32_T c4_c_i;
  int32_T c4_c_j;
  boolean_T c4_c_found;
  int32_T c4_c_ilo;
  int32_T c4_e_ihi;
  int32_T c4_c_a;
  int32_T c4_c_b;
  int32_T c4_d_a;
  int32_T c4_d_b;
  boolean_T c4_b_overflow;
  int32_T c4_c_jj;
  int32_T c4_d_jj;
  real_T c4_b_nzcount;
  int32_T c4_d_ilo;
  int32_T c4_f_ihi;
  int32_T c4_e_a;
  int32_T c4_e_b;
  int32_T c4_f_a;
  int32_T c4_f_b;
  boolean_T c4_c_overflow;
  int32_T c4_b_ii;
  int32_T c4_c_ii;
  boolean_T c4_c_A;
  int32_T c4_m;
  int32_T c4_d_i;
  int32_T c4_d_j;
  int32_T c4_e_ilo;
  int32_T c4_g_ihi;
  int32_T c4_f_ilo;
  int32_T c4_g_a;
  int32_T c4_h_a;
  boolean_T c4_d_overflow;
  int32_T c4_k;
  int32_T c4_b_k;
  creal_T c4_atmp;
  int32_T c4_h_ihi;
  int32_T c4_g_b;
  int32_T c4_h_b;
  boolean_T c4_e_overflow;
  int32_T c4_c_k;
  int32_T c4_i_a;
  int32_T c4_j_a;
  int32_T c4_b_m;
  int32_T c4_e_i;
  int32_T c4_e_j;
  int32_T c4_i_ihi;
  int32_T c4_d_k;
  int32_T c4_e_k;
  int32_T c4_j_ihi;
  int32_T c4_i_b;
  int32_T c4_j_b;
  boolean_T c4_f_overflow;
  int32_T c4_f_k;
  int32_T c4_k_a;
  int32_T c4_l_a;
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
  for (c4_i108 = 0; c4_i108 < 3; c4_i108++) {
    c4_rscale[c4_i108] = 0;
  }

  *c4_ilo = 1;
  *c4_ihi = 3;
  do {
    exitg2 = 0;
    c4_b_ihi = *c4_ihi;
    c4_i = 0;
    c4_j = 0;
    c4_found = false;
    c4_ii = c4_b_ihi;
    exitg5 = false;
    while ((exitg5 == false) && (c4_ii > 0)) {
      c4_nzcount = 0.0;
      c4_i = c4_ii;
      c4_j = c4_b_ihi;
      c4_c_ihi = c4_b_ihi;
      c4_b = c4_c_ihi;
      c4_b_b = c4_b;
      if (1 > c4_b_b) {
        c4_overflow = false;
      } else {
        c4_eml_switch_helper(chartInstance);
        c4_overflow = (c4_b_b > 2147483646);
      }

      if (c4_overflow) {
        c4_check_forloop_overflow_error(chartInstance, c4_overflow);
      }

      c4_jj = 1;
      exitg6 = false;
      while ((exitg6 == false) && (c4_jj <= c4_c_ihi)) {
        c4_b_jj = c4_jj;
        c4_b_A = ((c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_ii), 1, 3, 1, 0) + 3 *
                         (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_b_jj), 1, 3, 2, 0) - 1)) - 1].re !=
                   c4_dc7.re) || (c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_ii), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_jj), 1, 3, 2, 0) - 1)) - 1].im !=
                   c4_dc7.im));
        guard3 = false;
        guard4 = false;
        if (c4_b_A) {
          guard4 = true;
        } else if (c4_ii == c4_b_jj) {
          guard4 = true;
        } else {
          guard3 = true;
        }

        if (guard4 == true) {
          if (c4_nzcount == 0.0) {
            c4_j = c4_b_jj;
            c4_nzcount = 1.0;
            guard3 = true;
          } else {
            c4_nzcount = 2.0;
            exitg6 = true;
          }
        }

        if (guard3 == true) {
          c4_jj++;
        }
      }

      if (c4_nzcount < 2.0) {
        c4_found = true;
        exitg5 = true;
      } else {
        c4_a = c4_ii;
        c4_b_a = c4_a - 1;
        c4_ii = c4_b_a;
      }
    }

    c4_b_i = c4_i;
    c4_b_j = c4_j;
    c4_b_found = c4_found;
    if (!c4_b_found) {
      exitg2 = 2;
    } else {
      c4_b_m = *c4_ihi;
      c4_e_i = c4_b_i;
      c4_e_j = c4_b_j;
      c4_i_ihi = *c4_ihi;
      if (c4_e_i != c4_b_m) {
        c4_eml_switch_helper(chartInstance);
        for (c4_d_k = 1; c4_d_k < 4; c4_d_k++) {
          c4_e_k = c4_d_k;
          c4_atmp.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_e_i), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_e_k), 1, 3, 2, 0) - 1)) - 1].re;
          c4_atmp.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_e_i), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_e_k), 1, 3, 2, 0) - 1)) - 1].im;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_e_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_e_k), 1, 3, 2, 0) - 1)) - 1].re = c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_b_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_e_k), 1, 3, 2, 0)
               - 1)) - 1].re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_e_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_e_k), 1, 3, 2, 0) - 1)) - 1].im = c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_b_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_e_k), 1, 3, 2, 0)
               - 1)) - 1].im;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_m), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_e_k), 1, 3, 2, 0) - 1)) - 1].re = c4_atmp.re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_m), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_e_k), 1, 3, 2, 0) - 1)) - 1].im = c4_atmp.im;
        }
      }

      if (c4_e_j != c4_b_m) {
        c4_j_ihi = c4_i_ihi;
        c4_i_b = c4_j_ihi;
        c4_j_b = c4_i_b;
        if (1 > c4_j_b) {
          c4_f_overflow = false;
        } else {
          c4_eml_switch_helper(chartInstance);
          c4_f_overflow = (c4_j_b > 2147483646);
        }

        if (c4_f_overflow) {
          c4_check_forloop_overflow_error(chartInstance, c4_f_overflow);
        }

        for (c4_f_k = 1; c4_f_k <= c4_j_ihi; c4_f_k++) {
          c4_e_k = c4_f_k;
          c4_atmp.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_e_k), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_e_j), 1, 3, 2, 0) - 1)) - 1].re;
          c4_atmp.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_e_k), 1, 3, 1, 0) + 3 *
                             (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_e_j), 1, 3, 2, 0) - 1)) - 1].im;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_e_j), 1, 3, 2, 0) - 1)) - 1].re = c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_e_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_m), 1, 3, 2, 0)
               - 1)) - 1].re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_e_j), 1, 3, 2, 0) - 1)) - 1].im = c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_e_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_m), 1, 3, 2, 0)
               - 1)) - 1].im;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_m), 1, 3, 2, 0) - 1)) - 1].re = c4_atmp.re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_e_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_m), 1, 3, 2, 0) - 1)) - 1].im = c4_atmp.im;
        }
      }

      c4_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)*c4_ihi), 1, 3, 1, 0) - 1] = c4_b_j;
      c4_k_a = *c4_ihi;
      c4_l_a = c4_k_a - 1;
      *c4_ihi = c4_l_a;
      if (*c4_ihi == 1) {
        c4_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)*c4_ihi), 1, 3, 1, 0) - 1] = *c4_ihi;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    do {
      exitg1 = 0;
      c4_b_ilo = *c4_ilo;
      c4_d_ihi = *c4_ihi;
      c4_c_i = 0;
      c4_c_j = 0;
      c4_c_found = false;
      c4_c_ilo = c4_b_ilo;
      c4_e_ihi = c4_d_ihi;
      c4_c_a = c4_c_ilo;
      c4_c_b = c4_e_ihi;
      c4_d_a = c4_c_a;
      c4_d_b = c4_c_b;
      if (c4_d_a > c4_d_b) {
        c4_b_overflow = false;
      } else {
        c4_eml_switch_helper(chartInstance);
        c4_b_overflow = (c4_d_b > 2147483646);
      }

      if (c4_b_overflow) {
        c4_check_forloop_overflow_error(chartInstance, c4_b_overflow);
      }

      c4_c_jj = c4_c_ilo;
      exitg3 = false;
      while ((exitg3 == false) && (c4_c_jj <= c4_e_ihi)) {
        c4_d_jj = c4_c_jj;
        c4_b_nzcount = 0.0;
        c4_c_i = c4_d_ihi;
        c4_c_j = c4_d_jj;
        c4_d_ilo = c4_b_ilo;
        c4_f_ihi = c4_d_ihi;
        c4_e_a = c4_d_ilo;
        c4_e_b = c4_f_ihi;
        c4_f_a = c4_e_a;
        c4_f_b = c4_e_b;
        if (c4_f_a > c4_f_b) {
          c4_c_overflow = false;
        } else {
          c4_eml_switch_helper(chartInstance);
          c4_c_overflow = (c4_f_b > 2147483646);
        }

        if (c4_c_overflow) {
          c4_check_forloop_overflow_error(chartInstance, c4_c_overflow);
        }

        c4_b_ii = c4_d_ilo;
        exitg4 = false;
        while ((exitg4 == false) && (c4_b_ii <= c4_f_ihi)) {
          c4_c_ii = c4_b_ii;
          c4_c_A = ((c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_c_ii), 1, 3, 1, 0) + 3 *
                           (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_jj), 1, 3, 2, 0) - 1)) - 1].re
                     != c4_dc7.re) || (c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_c_ii), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_jj), 1, 3, 2, 0) - 1)) - 1].im
                     != c4_dc7.im));
          guard1 = false;
          guard2 = false;
          if (c4_c_A) {
            guard2 = true;
          } else if (c4_c_ii == c4_d_jj) {
            guard2 = true;
          } else {
            guard1 = true;
          }

          if (guard2 == true) {
            if (c4_b_nzcount == 0.0) {
              c4_c_i = c4_c_ii;
              c4_b_nzcount = 1.0;
              guard1 = true;
            } else {
              c4_b_nzcount = 2.0;
              exitg4 = true;
            }
          }

          if (guard1 == true) {
            c4_b_ii++;
          }
        }

        if (c4_b_nzcount < 2.0) {
          c4_c_found = true;
          exitg3 = true;
        } else {
          c4_c_jj++;
        }
      }

      c4_b_i = c4_c_i;
      c4_b_j = c4_c_j;
      c4_b_found = c4_c_found;
      if (!c4_b_found) {
        exitg1 = 1;
      } else {
        c4_m = *c4_ilo;
        c4_d_i = c4_b_i;
        c4_d_j = c4_b_j;
        c4_e_ilo = *c4_ilo;
        c4_g_ihi = *c4_ihi;
        if (c4_d_i != c4_m) {
          c4_f_ilo = c4_e_ilo;
          c4_g_a = c4_f_ilo;
          c4_h_a = c4_g_a;
          if (c4_h_a > 3) {
            c4_d_overflow = false;
          } else {
            c4_eml_switch_helper(chartInstance);
            c4_d_overflow = false;
          }

          if (c4_d_overflow) {
            c4_check_forloop_overflow_error(chartInstance, c4_d_overflow);
          }

          for (c4_k = c4_f_ilo; c4_k < 4; c4_k++) {
            c4_b_k = c4_k;
            c4_atmp.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 2, 0) - 1)) - 1].re;
            c4_atmp.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 2, 0) - 1)) - 1].im;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_b_k), 1, 3, 2, 0) - 1)) - 1].re = c4_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                  "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 2,
                  0) - 1)) - 1].re;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_b_k), 1, 3, 2, 0) - 1)) - 1].im = c4_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_m), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                  "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 2,
                  0) - 1)) - 1].im;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_m), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_b_k), 1, 3, 2, 0) - 1)) - 1].re = c4_atmp.re;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_m), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_b_k), 1, 3, 2, 0) - 1)) - 1].im = c4_atmp.im;
          }
        }

        if (c4_d_j != c4_m) {
          c4_h_ihi = c4_g_ihi;
          c4_g_b = c4_h_ihi;
          c4_h_b = c4_g_b;
          if (1 > c4_h_b) {
            c4_e_overflow = false;
          } else {
            c4_eml_switch_helper(chartInstance);
            c4_e_overflow = (c4_h_b > 2147483646);
          }

          if (c4_e_overflow) {
            c4_check_forloop_overflow_error(chartInstance, c4_e_overflow);
          }

          for (c4_c_k = 1; c4_c_k <= c4_h_ihi; c4_c_k++) {
            c4_b_k = c4_c_k;
            c4_atmp.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) - 1].re;
            c4_atmp.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) - 1].im;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) - 1].re = c4_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_m), 1, 3, 2, 0) - 1)) - 1].re;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) - 1].im = c4_A
              [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_m), 1, 3, 2, 0) - 1)) - 1].im;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_m), 1, 3, 2, 0) - 1)) - 1].re = c4_atmp.re;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_b_k), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_m), 1, 3, 2, 0) - 1)) - 1].im = c4_atmp.im;
          }
        }

        c4_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)*c4_ilo), 1, 3, 1, 0) - 1] = c4_b_j;
        c4_i_a = *c4_ilo;
        c4_j_a = c4_i_a + 1;
        *c4_ilo = c4_j_a;
        if (*c4_ilo == *c4_ihi) {
          c4_rscale[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
            "", (real_T)*c4_ilo), 1, 3, 1, 0) - 1] = *c4_ilo;
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

static void c4_b_sqrt(SFc4_Model_01InstanceStruct *chartInstance, creal_T *c4_x)
{
  real_T c4_yr;
  real_T c4_yi;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_z;
  real_T c4_e_x;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_b_z;
  boolean_T c4_b3;
  boolean_T c4_b4;
  boolean_T c4_b;
  real_T c4_h_x;
  boolean_T c4_b_b;
  real_T c4_i_x;
  boolean_T c4_c_b;
  real_T c4_absxr;
  real_T c4_absxi;
  real_T c4_j_x;
  real_T c4_y;
  real_T c4_x1;
  real_T c4_x2;
  real_T c4_a;
  real_T c4_d_b;
  real_T c4_absxd2;
  real_T c4_k_x;
  real_T c4_b_y;
  real_T c4_l_x;
  real_T c4_c_y;
  real_T c4_m_x;
  real_T c4_d_y;
  real_T c4_c_z;
  real_T c4_n_x;
  real_T c4_e_y;
  real_T c4_b_x1;
  real_T c4_b_x2;
  real_T c4_b_a;
  real_T c4_e_b;
  real_T c4_d_z;
  real_T c4_o_x;
  real_T c4_f_y;
  real_T c4_p_x;
  real_T c4_g_y;
  real_T c4_q_x;
  real_T c4_h_y;
  real_T c4_e_z;
  real_T c4_r_x;
  real_T c4_i_y;
  real_T c4_s_x;
  real_T c4_j_y;
  real_T c4_t_x;
  real_T c4_k_y;
  real_T c4_f_z;
  boolean_T guard1 = false;
  if (c4_x->im == 0.0) {
    if (c4_x->re < 0.0) {
      c4_yr = 0.0;
      c4_yi = muDoubleScalarSqrt(muDoubleScalarAbs(c4_x->re));
    } else {
      c4_yr = muDoubleScalarSqrt(c4_x->re);
      c4_yi = 0.0;
    }
  } else if (c4_x->re == 0.0) {
    if (c4_x->im < 0.0) {
      c4_b_x = -c4_x->im;
      c4_c_x = c4_b_x;
      c4_d_x = c4_c_x;
      c4_z = c4_d_x / 2.0;
      c4_yr = muDoubleScalarSqrt(c4_z);
      c4_yi = -c4_yr;
    } else {
      c4_e_x = c4_x->im;
      c4_f_x = c4_e_x;
      c4_g_x = c4_f_x;
      c4_b_z = c4_g_x / 2.0;
      c4_yr = muDoubleScalarSqrt(c4_b_z);
      c4_yi = c4_yr;
    }
  } else {
    c4_b3 = muDoubleScalarIsNaN(c4_x->re);
    c4_b4 = muDoubleScalarIsNaN(c4_x->im);
    c4_b = (c4_b3 || c4_b4);
    if (c4_b) {
      c4_yr = rtNaN;
      c4_yi = rtNaN;
    } else {
      c4_h_x = c4_x->im;
      c4_b_b = muDoubleScalarIsInf(c4_h_x);
      if (c4_b_b) {
        c4_yr = rtInf;
        c4_yi = c4_x->im;
      } else {
        c4_i_x = c4_x->re;
        c4_c_b = muDoubleScalarIsInf(c4_i_x);
        if (c4_c_b) {
          if (c4_x->re < 0.0) {
            c4_yr = 0.0;
            c4_yi = rtInf;
          } else {
            c4_yr = rtInf;
            c4_yi = 0.0;
          }
        } else {
          c4_absxr = muDoubleScalarAbs(c4_x->re);
          c4_absxi = muDoubleScalarAbs(c4_x->im);
          c4_realmax(chartInstance);
          guard1 = false;
          if (c4_absxr > 4.4942328371557893E+307) {
            guard1 = true;
          } else {
            c4_realmax(chartInstance);
            if (c4_absxi > 4.4942328371557893E+307) {
              guard1 = true;
            } else {
              c4_n_x = c4_absxr;
              c4_e_y = c4_absxi;
              c4_eml_scalar_eg(chartInstance);
              c4_b_x1 = c4_n_x;
              c4_b_x2 = c4_e_y;
              c4_b_a = c4_b_x1;
              c4_e_b = c4_b_x2;
              c4_d_z = muDoubleScalarHypot(c4_b_a, c4_e_b);
              c4_yr = muDoubleScalarSqrt((c4_d_z + c4_absxr) * 0.5);
            }
          }

          if (guard1 == true) {
            c4_absxr *= 0.5;
            c4_absxi *= 0.5;
            c4_j_x = c4_absxr;
            c4_y = c4_absxi;
            c4_eml_scalar_eg(chartInstance);
            c4_x1 = c4_j_x;
            c4_x2 = c4_y;
            c4_a = c4_x1;
            c4_d_b = c4_x2;
            c4_absxd2 = muDoubleScalarHypot(c4_a, c4_d_b);
            if (c4_absxd2 > c4_absxr) {
              c4_k_x = c4_absxr;
              c4_b_y = c4_absxd2;
              c4_l_x = c4_k_x;
              c4_c_y = c4_b_y;
              c4_m_x = c4_l_x;
              c4_d_y = c4_c_y;
              c4_c_z = c4_m_x / c4_d_y;
              c4_yr = muDoubleScalarSqrt(c4_absxd2) * muDoubleScalarSqrt(1.0 +
                c4_c_z);
            } else {
              c4_yr = muDoubleScalarSqrt(c4_absxd2) * 1.4142135623730951;
            }
          }

          if (c4_x->re > 0.0) {
            c4_o_x = c4_x->im;
            c4_f_y = c4_yr;
            c4_p_x = c4_o_x;
            c4_g_y = c4_f_y;
            c4_q_x = c4_p_x;
            c4_h_y = c4_g_y;
            c4_e_z = c4_q_x / c4_h_y;
            c4_yi = 0.5 * c4_e_z;
          } else {
            if (c4_x->im < 0.0) {
              c4_yi = -c4_yr;
            } else {
              c4_yi = c4_yr;
            }

            c4_r_x = c4_x->im;
            c4_i_y = c4_yi;
            c4_s_x = c4_r_x;
            c4_j_y = c4_i_y;
            c4_t_x = c4_s_x;
            c4_k_y = c4_j_y;
            c4_f_z = c4_t_x / c4_k_y;
            c4_yr = 0.5 * c4_f_z;
          }
        }
      }
    }
  }

  c4_x->re = c4_yr;
  c4_x->im = c4_yi;
}

static void c4_d_eml_matlab_zlascl(SFc4_Model_01InstanceStruct *chartInstance,
  real_T c4_cfrom, real_T c4_cto, creal_T c4_A[3])
{
  real_T c4_cfromc;
  real_T c4_ctoc;
  boolean_T c4_notdone;
  real_T c4_cfrom1;
  real_T c4_cto1;
  real_T c4_x;
  real_T c4_b_x;
  real_T c4_y;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_b_y;
  real_T c4_mul;
  real_T c4_e_x;
  real_T c4_f_x;
  real_T c4_c_y;
  real_T c4_g_x;
  real_T c4_h_x;
  real_T c4_d_y;
  real_T c4_a;
  int32_T c4_i109;
  boolean_T guard1 = false;
  c4_realmin(chartInstance);
  c4_eps(chartInstance);
  c4_cfromc = c4_cfrom;
  c4_ctoc = c4_cto;
  c4_notdone = true;
  while (c4_notdone) {
    c4_cfrom1 = c4_cfromc * 2.0041683600089728E-292;
    c4_cto1 = c4_ctoc / 4.9896007738368E+291;
    c4_x = c4_cfrom1;
    c4_b_x = c4_x;
    c4_y = muDoubleScalarAbs(c4_b_x);
    c4_c_x = c4_ctoc;
    c4_d_x = c4_c_x;
    c4_b_y = muDoubleScalarAbs(c4_d_x);
    guard1 = false;
    if (c4_y > c4_b_y) {
      if (c4_ctoc != 0.0) {
        c4_mul = 2.0041683600089728E-292;
        c4_notdone = true;
        c4_cfromc = c4_cfrom1;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      c4_e_x = c4_cto1;
      c4_f_x = c4_e_x;
      c4_c_y = muDoubleScalarAbs(c4_f_x);
      c4_g_x = c4_cfromc;
      c4_h_x = c4_g_x;
      c4_d_y = muDoubleScalarAbs(c4_h_x);
      if (c4_c_y > c4_d_y) {
        c4_mul = 4.9896007738368E+291;
        c4_notdone = true;
        c4_ctoc = c4_cto1;
      } else {
        c4_mul = c4_ctoc / c4_cfromc;
        c4_notdone = false;
      }
    }

    c4_a = c4_mul;
    for (c4_i109 = 0; c4_i109 < 3; c4_i109++) {
      c4_A[c4_i109].re *= c4_a;
      c4_A[c4_i109].im *= c4_a;
    }
  }
}

static void c4_b_eml_xgemm(SFc4_Model_01InstanceStruct *chartInstance, real_T
  c4_A[9], real_T c4_B[9], real_T c4_C[9])
{
  int32_T c4_i110;
  int32_T c4_i111;
  int32_T c4_i112;
  int32_T c4_i113;
  int32_T c4_i114;
  (void)chartInstance;
  for (c4_i110 = 0; c4_i110 < 3; c4_i110++) {
    c4_i111 = 0;
    for (c4_i112 = 0; c4_i112 < 3; c4_i112++) {
      c4_C[c4_i111 + c4_i110] = 0.0;
      c4_i113 = 0;
      for (c4_i114 = 0; c4_i114 < 3; c4_i114++) {
        c4_C[c4_i111 + c4_i110] += c4_A[c4_i113 + c4_i110] * c4_B[c4_i114 +
          c4_i111];
        c4_i113 += 3;
      }

      c4_i111 += 3;
    }
  }
}

static void c4_b_eml_matlab_zgghrd(SFc4_Model_01InstanceStruct *chartInstance,
  int32_T c4_ilo, int32_T c4_ihi, creal_T c4_A[9], creal_T c4_Z[9])
{
  int32_T c4_i115;
  static real_T c4_dv4[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  int32_T c4_a;
  int32_T c4_b_a;
  int32_T c4_c;
  int32_T c4_c_a;
  int32_T c4_d_a;
  int32_T c4_ihim1;
  int32_T c4_jcol;
  int32_T c4_e_a;
  int32_T c4_f_a;
  int32_T c4_jcolp1;
  int32_T c4_jrow;
  int32_T c4_g_a;
  int32_T c4_h_a;
  int32_T c4_jrowm1;
  creal_T c4_b_A;
  creal_T c4_c_A;
  creal_T c4_b;
  creal_T c4_s;
  real_T c4_b_c;
  real_T c4_c_c;
  static creal_T c4_dc8 = { 0.0, 0.0 };

  real_T c4_d_c;
  int32_T c4_xrow;
  int32_T c4_yrow;
  int32_T c4_jlo;
  int32_T c4_jhi;
  int32_T c4_b_jlo;
  int32_T c4_b_jhi;
  int32_T c4_i_a;
  int32_T c4_b_b;
  int32_T c4_j_a;
  int32_T c4_c_b;
  boolean_T c4_overflow;
  int32_T c4_j;
  int32_T c4_b_j;
  real_T c4_k_a;
  creal_T c4_y;
  creal_T c4_b_s;
  creal_T c4_stemp;
  real_T c4_l_a;
  creal_T c4_d_b;
  creal_T c4_e_b;
  creal_T c4_f_b;
  creal_T c4_g_b;
  real_T c4_e_c;
  int32_T c4_xcol;
  int32_T c4_ycol;
  int32_T c4_b_ilo;
  int32_T c4_b_ihi;
  int32_T c4_c_ilo;
  int32_T c4_c_ihi;
  int32_T c4_m_a;
  int32_T c4_h_b;
  int32_T c4_n_a;
  int32_T c4_i_b;
  boolean_T c4_b_overflow;
  int32_T c4_i;
  int32_T c4_b_i;
  real_T c4_o_a;
  creal_T c4_c_s;
  real_T c4_p_a;
  creal_T c4_j_b;
  creal_T c4_k_b;
  creal_T c4_l_b;
  creal_T c4_m_b;
  real_T c4_f_c;
  int32_T c4_b_xcol;
  int32_T c4_b_ycol;
  int32_T c4_c_i;
  int32_T c4_d_i;
  real_T c4_q_a;
  creal_T c4_d_s;
  real_T c4_r_a;
  creal_T c4_n_b;
  creal_T c4_o_b;
  creal_T c4_p_b;
  creal_T c4_q_b;
  for (c4_i115 = 0; c4_i115 < 9; c4_i115++) {
    c4_Z[c4_i115].re = c4_dv4[c4_i115];
    c4_Z[c4_i115].im = 0.0;
  }

  c4_a = c4_ilo;
  c4_b_a = c4_a;
  c4_c = c4_b_a;
  if (c4_ihi < c4_c + 2) {
  } else {
    c4_c_a = c4_ihi;
    c4_d_a = c4_c_a;
    c4_ihim1 = c4_d_a;
    c4_jcol = c4_ilo;
    while (c4_jcol < c4_ihim1 - 1) {
      c4_e_a = c4_jcol;
      c4_f_a = c4_e_a + 1;
      c4_jcolp1 = c4_f_a;
      c4_jrow = c4_ihi;
      while (c4_jrow > c4_jcolp1) {
        c4_g_a = c4_jrow;
        c4_h_a = c4_g_a - 1;
        c4_jrowm1 = c4_h_a;
        c4_b_A.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_jrowm1), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].re;
        c4_b_A.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_jrowm1), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].im;
        c4_c_A.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_jrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].re;
        c4_c_A.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_jrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].im;
        c4_eml_matlab_zlartg(chartInstance, c4_b_A, c4_c_A, &c4_b_c, &c4_s,
                             &c4_b);
        c4_c_c = c4_b_c;
        c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_jrowm1), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].re = c4_b.re;
        c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_jrowm1), 1, 3, 1, 0) + 3 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c4_jcol), 1, 3, 2, 0) - 1)) - 1].im = c4_b.im;
        c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_jrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0)
               - 1)) - 1].re = c4_dc8.re;
        c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_jrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_jcol), 1, 3, 2, 0)
               - 1)) - 1].im = c4_dc8.im;
        c4_d_c = c4_c_c;
        c4_xrow = c4_jrowm1;
        c4_yrow = c4_jrow;
        c4_jlo = c4_jcolp1;
        c4_jhi = c4_ihi;
        c4_b_jlo = c4_jlo;
        c4_b_jhi = c4_jhi;
        c4_i_a = c4_b_jlo;
        c4_b_b = c4_b_jhi;
        c4_j_a = c4_i_a;
        c4_c_b = c4_b_b;
        if (c4_j_a > c4_c_b) {
          c4_overflow = false;
        } else {
          c4_eml_switch_helper(chartInstance);
          c4_overflow = (c4_c_b > 2147483646);
        }

        if (c4_overflow) {
          c4_check_forloop_overflow_error(chartInstance, c4_overflow);
        }

        for (c4_j = c4_b_jlo; c4_j <= c4_b_jhi; c4_j++) {
          c4_b_j = c4_j;
          c4_k_a = c4_d_c;
          c4_b.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
          c4_b.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
          c4_y.re = c4_k_a * c4_b.re;
          c4_y.im = c4_k_a * c4_b.im;
          c4_b_s.re = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re - c4_s.im * c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_yrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0)
               - 1)) - 1].im;
          c4_b_s.im = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im + c4_s.im * c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_yrow), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0)
               - 1)) - 1].re;
          c4_stemp.re = c4_y.re + c4_b_s.re;
          c4_stemp.im = c4_y.im + c4_b_s.im;
          c4_l_a = c4_d_c;
          c4_b.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
          c4_b.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
          c4_y.re = c4_l_a * c4_b.re;
          c4_y.im = c4_l_a * c4_b.im;
          c4_b.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
          c4_b.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
          c4_d_b = c4_b;
          c4_e_b = c4_b;
          c4_f_b = c4_b;
          c4_g_b = c4_b;
          c4_b.re = c4_s.re * c4_d_b.re + c4_s.im * c4_e_b.im;
          c4_b.im = c4_s.re * c4_f_b.im - c4_s.im * c4_g_b.re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re = c4_y.re -
            c4_b.re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im = c4_y.im -
            c4_b.im;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].re = c4_stemp.re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1].im = c4_stemp.im;
        }

        c4_s.re = -c4_s.re;
        c4_s.im = -c4_s.im;
        c4_e_c = c4_c_c;
        c4_xcol = c4_jrow;
        c4_ycol = c4_jrowm1;
        c4_b_ilo = c4_ilo;
        c4_b_ihi = c4_ihi;
        c4_c_ilo = c4_b_ilo;
        c4_c_ihi = c4_b_ihi;
        c4_m_a = c4_c_ilo;
        c4_h_b = c4_c_ihi;
        c4_n_a = c4_m_a;
        c4_i_b = c4_h_b;
        if (c4_n_a > c4_i_b) {
          c4_b_overflow = false;
        } else {
          c4_eml_switch_helper(chartInstance);
          c4_b_overflow = (c4_i_b > 2147483646);
        }

        if (c4_b_overflow) {
          c4_check_forloop_overflow_error(chartInstance, c4_b_overflow);
        }

        for (c4_i = c4_c_ilo; c4_i <= c4_c_ihi; c4_i++) {
          c4_b_i = c4_i;
          c4_o_a = c4_e_c;
          c4_b.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c4_b.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c4_y.re = c4_o_a * c4_b.re;
          c4_y.im = c4_o_a * c4_b.im;
          c4_c_s.re = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re - c4_s.im * c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0)
               - 1)) - 1].im;
          c4_c_s.im = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im + c4_s.im * c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_b_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0)
               - 1)) - 1].re;
          c4_stemp.re = c4_y.re + c4_c_s.re;
          c4_stemp.im = c4_y.im + c4_c_s.im;
          c4_p_a = c4_e_c;
          c4_b.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re;
          c4_b.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im;
          c4_y.re = c4_p_a * c4_b.re;
          c4_y.im = c4_p_a * c4_b.im;
          c4_b.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c4_b.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c4_j_b = c4_b;
          c4_k_b = c4_b;
          c4_l_b = c4_b;
          c4_m_b = c4_b;
          c4_b.re = c4_s.re * c4_j_b.re + c4_s.im * c4_k_b.im;
          c4_b.im = c4_s.re * c4_l_b.im - c4_s.im * c4_m_b.re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re = c4_y.re -
            c4_b.re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im = c4_y.im -
            c4_b.im;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].re = c4_stemp.re;
          c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].im = c4_stemp.im;
        }

        c4_f_c = c4_c_c;
        c4_b_xcol = c4_jrow;
        c4_b_ycol = c4_jrowm1;
        c4_eml_switch_helper(chartInstance);
        for (c4_c_i = 1; c4_c_i < 4; c4_c_i++) {
          c4_d_i = c4_c_i;
          c4_q_a = c4_f_c;
          c4_b.re = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c4_b.im = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c4_y.re = c4_q_a * c4_b.re;
          c4_y.im = c4_q_a * c4_b.im;
          c4_d_s.re = c4_s.re * c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].re - c4_s.im * c4_Z
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_d_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2,
                0) - 1)) - 1].im;
          c4_d_s.im = c4_s.re * c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].im + c4_s.im * c4_Z
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c4_d_i), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2,
                0) - 1)) - 1].re;
          c4_stemp.re = c4_y.re + c4_d_s.re;
          c4_stemp.im = c4_y.im + c4_d_s.im;
          c4_r_a = c4_f_c;
          c4_b.re = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].re;
          c4_b.im = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].im;
          c4_y.re = c4_r_a * c4_b.re;
          c4_y.im = c4_r_a * c4_b.im;
          c4_b.re = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1].re;
          c4_b.im = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1].im;
          c4_n_b = c4_b;
          c4_o_b = c4_b;
          c4_p_b = c4_b;
          c4_q_b = c4_b;
          c4_b.re = c4_s.re * c4_n_b.re + c4_s.im * c4_o_b.im;
          c4_b.im = c4_s.re * c4_p_b.im - c4_s.im * c4_q_b.re;
          c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].re = c4_y.re -
            c4_b.re;
          c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].im = c4_y.im -
            c4_b.im;
          c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1].re = c4_stemp.re;
          c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1].im = c4_stemp.im;
        }

        c4_jrow = c4_jrowm1;
      }

      c4_jcol = c4_jcolp1;
    }
  }
}

static void c4_c_eml_matlab_zhgeqz(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], int32_T c4_ilo, int32_T c4_ihi, creal_T c4_Z[9], real_T
  *c4_info, creal_T c4_alpha1[3], creal_T c4_beta1[3])
{
  static creal_T c4_dc9 = { 0.0, 0.0 };

  int32_T c4_i116;
  int32_T c4_i117;
  static creal_T c4_dc10 = { 0.0, 0.0 };

  creal_T c4_eshift;
  creal_T c4_ctemp;
  creal_T c4_rho;
  int32_T c4_i118;
  int32_T c4_i119;
  int32_T c4_i120;
  creal_T c4_b_A[9];
  real_T c4_anorm;
  real_T c4_y;
  real_T c4_atol;
  real_T c4_b_y;
  real_T c4_x;
  real_T c4_ascale;
  boolean_T c4_failed;
  int32_T c4_a;
  int32_T c4_b_a;
  int32_T c4_i121;
  int32_T c4_c_a;
  int32_T c4_d_a;
  boolean_T c4_overflow;
  int32_T c4_j;
  int32_T c4_b_j;
  int32_T c4_ifirst;
  int32_T c4_istart;
  int32_T c4_ilast;
  int32_T c4_e_a;
  int32_T c4_f_a;
  int32_T c4_ilastm1;
  int32_T c4_iiter;
  int32_T c4_g_a;
  int32_T c4_b;
  int32_T c4_h_a;
  int32_T c4_b_b;
  int32_T c4_c;
  int32_T c4_i_a;
  int32_T c4_j_a;
  int32_T c4_b_c;
  int32_T c4_c_b;
  int32_T c4_d_b;
  int32_T c4_maxit;
  boolean_T c4_goto50;
  boolean_T c4_goto60;
  boolean_T c4_goto70;
  boolean_T c4_goto90;
  int32_T c4_b_maxit;
  int32_T c4_e_b;
  int32_T c4_f_b;
  boolean_T c4_b_overflow;
  int32_T c4_jiter;
  creal_T c4_a22;
  real_T c4_b_x;
  real_T c4_c_x;
  real_T c4_c_y;
  real_T c4_d_x;
  real_T c4_e_x;
  real_T c4_d_y;
  real_T c4_e_y;
  int32_T c4_k_a;
  int32_T c4_l_a;
  int32_T c4_jm1;
  boolean_T c4_ilazro;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_f_y;
  real_T c4_h_x;
  real_T c4_i_x;
  real_T c4_g_y;
  real_T c4_h_y;
  boolean_T c4_b5;
  int32_T c4_i122;
  int32_T c4_i123;
  int32_T c4_i124;
  creal_T c4_c_A;
  creal_T c4_d_A;
  creal_T c4_s;
  real_T c4_c_c;
  real_T c4_d_c;
  real_T c4_e_c;
  int32_T c4_xcol;
  int32_T c4_ycol;
  int32_T c4_b_ihi;
  int32_T c4_c_ihi;
  int32_T c4_g_b;
  int32_T c4_h_b;
  boolean_T c4_c_overflow;
  int32_T c4_i;
  int32_T c4_b_i;
  real_T c4_m_a;
  creal_T c4_a12;
  creal_T c4_b_s;
  creal_T c4_a21;
  real_T c4_n_a;
  creal_T c4_b_a22;
  creal_T c4_c_a22;
  creal_T c4_d_a22;
  creal_T c4_e_a22;
  real_T c4_f_c;
  int32_T c4_b_xcol;
  int32_T c4_b_ycol;
  int32_T c4_c_i;
  int32_T c4_d_i;
  real_T c4_o_a;
  creal_T c4_c_s;
  real_T c4_p_a;
  creal_T c4_f_a22;
  creal_T c4_g_a22;
  creal_T c4_h_a22;
  creal_T c4_i_a22;
  int32_T c4_q_a;
  int32_T c4_r_a;
  int32_T c4_s_a;
  int32_T c4_t_a;
  creal_T c4_r2;
  creal_T c4_j_a22;
  creal_T c4_b_rho;
  creal_T c4_b_a12;
  creal_T c4_c_a12;
  creal_T c4_b_a21;
  real_T c4_d7;
  real_T c4_d8;
  int32_T c4_u_a;
  int32_T c4_v_a;
  int32_T c4_jp1;
  int32_T c4_w_a;
  int32_T c4_x_a;
  real_T c4_j_x;
  real_T c4_k_x;
  real_T c4_i_y;
  real_T c4_l_x;
  real_T c4_m_x;
  real_T c4_j_y;
  real_T c4_k_y;
  real_T c4_temp;
  real_T c4_n_x;
  real_T c4_o_x;
  real_T c4_l_y;
  real_T c4_p_x;
  real_T c4_q_x;
  real_T c4_m_y;
  real_T c4_n_y;
  real_T c4_temp2;
  real_T c4_r_x;
  real_T c4_o_y;
  real_T c4_tempr;
  real_T c4_s_x;
  real_T c4_t_x;
  real_T c4_p_y;
  real_T c4_u_x;
  real_T c4_v_x;
  real_T c4_q_y;
  real_T c4_r_y;
  int32_T c4_y_a;
  int32_T c4_ab_a;
  int32_T c4_g_c;
  real_T c4_h_c;
  int32_T c4_bb_a;
  int32_T c4_cb_a;
  int32_T c4_db_a;
  int32_T c4_eb_a;
  creal_T c4_e_A;
  creal_T c4_f_A;
  real_T c4_i_c;
  real_T c4_j_c;
  int32_T c4_xrow;
  int32_T c4_yrow;
  int32_T c4_jlo;
  int32_T c4_b_jlo;
  int32_T c4_fb_a;
  int32_T c4_gb_a;
  boolean_T c4_d_overflow;
  int32_T c4_c_j;
  int32_T c4_d_j;
  real_T c4_hb_a;
  creal_T c4_d_s;
  real_T c4_ib_a;
  creal_T c4_k_a22;
  creal_T c4_l_a22;
  creal_T c4_m_a22;
  creal_T c4_n_a22;
  int32_T c4_jb_a;
  int32_T c4_kb_a;
  int32_T c4_k_c;
  int32_T c4_w_x;
  int32_T c4_s_y;
  int32_T c4_x_x;
  real_T c4_l_c;
  int32_T c4_c_xcol;
  int32_T c4_c_ycol;
  int32_T c4_d_ihi;
  int32_T c4_e_ihi;
  int32_T c4_i_b;
  int32_T c4_j_b;
  boolean_T c4_e_overflow;
  int32_T c4_e_i;
  int32_T c4_f_i;
  real_T c4_lb_a;
  creal_T c4_e_s;
  real_T c4_mb_a;
  creal_T c4_o_a22;
  creal_T c4_p_a22;
  creal_T c4_q_a22;
  creal_T c4_r_a22;
  real_T c4_m_c;
  int32_T c4_d_xcol;
  int32_T c4_d_ycol;
  int32_T c4_g_i;
  int32_T c4_h_i;
  real_T c4_nb_a;
  creal_T c4_f_s;
  real_T c4_ob_a;
  creal_T c4_s_a22;
  creal_T c4_t_a22;
  creal_T c4_u_a22;
  creal_T c4_v_a22;
  int32_T c4_b_ilast;
  int32_T c4_k_b;
  int32_T c4_l_b;
  boolean_T c4_f_overflow;
  int32_T c4_k;
  int32_T c4_b_k;
  int32_T c4_i125;
  int32_T c4_pb_a;
  int32_T c4_qb_a;
  int32_T c4_i126;
  int32_T c4_m_b;
  int32_T c4_n_b;
  boolean_T c4_g_overflow;
  int32_T c4_e_j;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  int32_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T guard11 = false;
  c4_dc9.re = rtNaN;
  for (c4_i116 = 0; c4_i116 < 3; c4_i116++) {
    c4_alpha1[c4_i116].re = 0.0;
    c4_alpha1[c4_i116].im = 0.0;
  }

  for (c4_i117 = 0; c4_i117 < 3; c4_i117++) {
    c4_beta1[c4_i117].re = 1.0;
    c4_beta1[c4_i117].im = 0.0;
  }

  c4_eps(chartInstance);
  c4_realmin(chartInstance);
  c4_eshift = c4_dc10;
  c4_ctemp = c4_dc10;
  c4_rho = c4_dc10;
  c4_i118 = 0;
  for (c4_i119 = 0; c4_i119 < 3; c4_i119++) {
    for (c4_i120 = 0; c4_i120 < 3; c4_i120++) {
      c4_b_A[c4_i120 + c4_i118] = c4_A[c4_i120 + c4_i118];
    }

    c4_i118 += 3;
  }

  c4_anorm = c4_eml_matlab_zlanhs(chartInstance, c4_b_A, c4_ilo, c4_ihi);
  c4_y = 2.2204460492503131E-16 * c4_anorm;
  c4_atol = 2.2250738585072014E-308;
  if (c4_y > 2.2250738585072014E-308) {
    c4_atol = c4_y;
  }

  c4_b_y = c4_anorm;
  c4_x = 2.2250738585072014E-308;
  if (c4_b_y > 2.2250738585072014E-308) {
    c4_x = c4_b_y;
  }

  c4_ascale = 1.0 / c4_x;
  c4_failed = true;
  c4_a = c4_ihi;
  c4_b_a = c4_a + 1;
  c4_i121 = c4_b_a;
  c4_c_a = c4_i121;
  c4_d_a = c4_c_a;
  if (c4_d_a > 3) {
    c4_overflow = false;
  } else {
    c4_eml_switch_helper(chartInstance);
    c4_overflow = false;
  }

  if (c4_overflow) {
    c4_check_forloop_overflow_error(chartInstance, c4_overflow);
  }

  for (c4_j = c4_i121; c4_j < 4; c4_j++) {
    c4_b_j = c4_j;
    c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_b_j), 1, 3, 1, 0) - 1].re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK(
      "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
    c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_b_j), 1, 3, 1, 0) - 1].im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK(
      "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
  }

  guard1 = false;
  guard2 = false;
  if (c4_ihi >= c4_ilo) {
    c4_ifirst = c4_ilo;
    c4_istart = c4_ilo;
    c4_ilast = c4_ihi;
    c4_e_a = c4_ilast;
    c4_f_a = c4_e_a - 1;
    c4_ilastm1 = c4_f_a;
    c4_iiter = 0;
    c4_g_a = c4_ihi;
    c4_b = c4_ilo;
    c4_h_a = c4_g_a;
    c4_b_b = c4_b;
    c4_c = c4_h_a - c4_b_b;
    c4_i_a = c4_c;
    c4_j_a = c4_i_a;
    c4_b_c = c4_j_a;
    c4_c_b = c4_b_c + 1;
    c4_d_b = c4_c_b;
    c4_maxit = 30 * c4_d_b;
    c4_goto50 = false;
    c4_goto60 = false;
    c4_goto70 = false;
    c4_goto90 = false;
    c4_b_maxit = c4_maxit;
    c4_e_b = c4_b_maxit;
    c4_f_b = c4_e_b;
    if (1 > c4_f_b) {
      c4_b_overflow = false;
    } else {
      c4_eml_switch_helper(chartInstance);
      c4_b_overflow = (c4_f_b > 2147483646);
    }

    if (c4_b_overflow) {
      c4_check_forloop_overflow_error(chartInstance, c4_b_overflow);
    }

    c4_jiter = 1;
    do {
      exitg1 = 0;
      if (c4_jiter <= c4_b_maxit) {
        if (c4_ilast == c4_ilo) {
          c4_goto60 = true;
        } else {
          c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            re;
          c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].
            im;
          c4_b_x = c4_a22.re;
          c4_c_x = c4_b_x;
          c4_c_y = muDoubleScalarAbs(c4_c_x);
          c4_d_x = c4_a22.im;
          c4_e_x = c4_d_x;
          c4_d_y = muDoubleScalarAbs(c4_e_x);
          c4_e_y = c4_c_y + c4_d_y;
          if (c4_e_y <= c4_atol) {
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c4_dc10.re;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c4_dc10.im;
            c4_goto60 = true;
          } else {
            c4_b_j = c4_ilastm1;
            exitg3 = false;
            while ((exitg3 == false) && (c4_b_j >= c4_ilo)) {
              c4_k_a = c4_b_j;
              c4_l_a = c4_k_a - 1;
              c4_jm1 = c4_l_a;
              if (c4_b_j == c4_ilo) {
                c4_ilazro = true;
              } else {
                c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c4_f_x = c4_a22.re;
                c4_g_x = c4_f_x;
                c4_f_y = muDoubleScalarAbs(c4_g_x);
                c4_h_x = c4_a22.im;
                c4_i_x = c4_h_x;
                c4_g_y = muDoubleScalarAbs(c4_i_x);
                c4_h_y = c4_f_y + c4_g_y;
                if (c4_h_y <= c4_atol) {
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) -
                           1)) - 1].re = c4_dc10.re;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) -
                           1)) - 1].im = c4_dc10.im;
                  c4_ilazro = true;
                } else {
                  c4_ilazro = false;
                }
              }

              if (c4_ilazro) {
                c4_ifirst = c4_b_j;
                c4_goto70 = true;
                exitg3 = true;
              } else {
                c4_b_j = c4_jm1;
              }
            }
          }
        }

        guard3 = false;
        guard4 = false;
        if (c4_goto50) {
          guard4 = true;
        } else if (c4_goto60) {
          guard4 = true;
        } else if (c4_goto70) {
          guard3 = true;
        } else {
          c4_b5 = false;
        }

        if (guard4 == true) {
          guard3 = true;
        }

        if (guard3 == true) {
          c4_b5 = true;
        }

        if (!c4_b5) {
          for (c4_i122 = 0; c4_i122 < 3; c4_i122++) {
            c4_alpha1[c4_i122] = c4_dc9;
          }

          for (c4_i123 = 0; c4_i123 < 3; c4_i123++) {
            c4_beta1[c4_i123] = c4_dc9;
          }

          for (c4_i124 = 0; c4_i124 < 9; c4_i124++) {
            c4_Z[c4_i124] = c4_dc9;
          }

          *c4_info = -1.0;
          exitg1 = 1;
        } else {
          if (c4_goto50) {
            c4_goto50 = false;
            c4_c_A.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].
              re;
            c4_c_A.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].
              im;
            c4_d_A.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .re;
            c4_d_A.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1]
              .im;
            c4_eml_matlab_zlartg(chartInstance, c4_c_A, c4_d_A, &c4_c_c, &c4_s,
                                 &c4_a22);
            c4_d_c = c4_c_c;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].re =
              c4_a22.re;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].im =
              c4_a22.im;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].re =
              c4_dc10.re;
            c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                     "", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1)) - 1].im =
              c4_dc10.im;
            c4_e_c = c4_d_c;
            c4_xcol = c4_ilast;
            c4_ycol = c4_ilastm1;
            c4_b_ihi = c4_ilastm1;
            c4_c_ihi = c4_b_ihi;
            c4_g_b = c4_c_ihi;
            c4_h_b = c4_g_b;
            if (1 > c4_h_b) {
              c4_c_overflow = false;
            } else {
              c4_eml_switch_helper(chartInstance);
              c4_c_overflow = (c4_h_b > 2147483646);
            }

            if (c4_c_overflow) {
              c4_check_forloop_overflow_error(chartInstance, c4_c_overflow);
            }

            for (c4_i = 1; c4_i <= c4_c_ihi; c4_i++) {
              c4_b_i = c4_i;
              c4_m_a = c4_e_c;
              c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c4_a12.re = c4_m_a * c4_a22.re;
              c4_a12.im = c4_m_a * c4_a22.im;
              c4_b_s.re = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re - c4_s.im *
                c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1))
                - 1].im;
              c4_b_s.im = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im + c4_s.im *
                c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1))
                - 1].re;
              c4_a21.re = c4_a12.re + c4_b_s.re;
              c4_a21.im = c4_a12.im + c4_b_s.im;
              c4_n_a = c4_e_c;
              c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c4_a12.re = c4_n_a * c4_a22.re;
              c4_a12.im = c4_n_a * c4_a22.im;
              c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].
                re;
              c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].
                im;
              c4_b_a22 = c4_a22;
              c4_c_a22 = c4_a22;
              c4_d_a22 = c4_a22;
              c4_e_a22 = c4_a22;
              c4_a22.re = c4_s.re * c4_b_a22.re + c4_s.im * c4_c_a22.im;
              c4_a22.im = c4_s.re * c4_d_a22.im - c4_s.im * c4_e_a22.re;
              c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].re =
                c4_a12.re - c4_a22.re;
              c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ycol), 1, 3, 2, 0) - 1)) - 1].im =
                c4_a12.im - c4_a22.im;
              c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].re =
                c4_a21.re;
              c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_b_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_xcol), 1, 3, 2, 0) - 1)) - 1].im =
                c4_a21.im;
            }

            c4_f_c = c4_d_c;
            c4_b_xcol = c4_ilast;
            c4_b_ycol = c4_ilastm1;
            c4_eml_switch_helper(chartInstance);
            for (c4_c_i = 1; c4_c_i < 4; c4_c_i++) {
              c4_d_i = c4_c_i;
              c4_o_a = c4_f_c;
              c4_a22.re = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c4_a22.im = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c4_a12.re = c4_o_a * c4_a22.re;
              c4_a12.im = c4_o_a * c4_a22.im;
              c4_c_s.re = c4_s.re * c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].re - c4_s.im *
                c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) -
                       1)) - 1].im;
              c4_c_s.im = c4_s.re * c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3
                * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].im + c4_s.im *
                c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                       _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) -
                       1)) - 1].re;
              c4_a21.re = c4_a12.re + c4_c_s.re;
              c4_a21.im = c4_a12.im + c4_c_s.im;
              c4_p_a = c4_f_c;
              c4_a22.re = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c4_a22.im = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c4_a12.re = c4_p_a * c4_a22.re;
              c4_a12.im = c4_p_a * c4_a22.im;
              c4_a22.re = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c4_a22.im = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c4_f_a22 = c4_a22;
              c4_g_a22 = c4_a22;
              c4_h_a22 = c4_a22;
              c4_i_a22 = c4_a22;
              c4_a22.re = c4_s.re * c4_f_a22.re + c4_s.im * c4_g_a22.im;
              c4_a22.im = c4_s.re * c4_h_a22.im - c4_s.im * c4_i_a22.re;
              c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].re =
                c4_a12.re - c4_a22.re;
              c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_ycol), 1, 3, 2, 0) - 1)) - 1].im =
                c4_a12.im - c4_a22.im;
              c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1].re =
                c4_a21.re;
              c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_d_i), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_b_xcol), 1, 3, 2, 0) - 1)) - 1].im =
                c4_a21.im;
            }

            c4_goto60 = true;
          }

          guard11 = false;
          if (c4_goto60) {
            c4_goto60 = false;
            c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) - 1].re =
              c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].re;
            c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) - 1].im =
              c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) - 1].im;
            c4_ilast = c4_ilastm1;
            c4_q_a = c4_ilast;
            c4_r_a = c4_q_a - 1;
            c4_ilastm1 = c4_r_a;
            if (c4_ilast < c4_ilo) {
              c4_failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              c4_iiter = 0;
              c4_eshift = c4_dc10;
              guard11 = true;
            }
          } else {
            if (c4_goto70) {
              c4_goto70 = false;
              c4_s_a = c4_iiter;
              c4_t_a = c4_s_a + 1;
              c4_iiter = c4_t_a;
              if (c4_mod(chartInstance, c4_iiter) != 0) {
                c4_s.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c4_s.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c4_r2.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c4_r2.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c4_a12.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) -
                  1].re;
                c4_a12.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 2, 0) - 1)) -
                  1].im;
                c4_a21.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c4_a21.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c4_a22.re = c4_r2.re - c4_s.re;
                c4_a22.im = c4_r2.im - c4_s.im;
                c4_j_a22.re = -c4_a22.re;
                c4_j_a22.im = -c4_a22.im;
                c4_rho = c4_eml_div(chartInstance, c4_j_a22, 2.0);
                c4_b_rho.re = c4_rho.re * c4_rho.re - c4_rho.im * c4_rho.im;
                c4_b_rho.im = c4_rho.re * c4_rho.im + c4_rho.im * c4_rho.re;
                c4_b_a12.re = c4_a12.re * c4_a21.re - c4_a12.im * c4_a21.im;
                c4_b_a12.im = c4_a12.re * c4_a21.im + c4_a12.im * c4_a21.re;
                c4_a22.re = c4_b_rho.re + c4_b_a12.re;
                c4_a22.im = c4_b_rho.im + c4_b_a12.im;
                c4_b_sqrt(chartInstance, &c4_a22);
                c4_a12.re = c4_s.re - (c4_rho.re - c4_a22.re);
                c4_a12.im = c4_s.im - (c4_rho.im - c4_a22.im);
                c4_a21.re = c4_s.re - (c4_rho.re + c4_a22.re);
                c4_a21.im = c4_s.im - (c4_rho.im + c4_a22.im);
                c4_c_a12.re = c4_a12.re - c4_r2.re;
                c4_c_a12.im = c4_a12.im - c4_r2.im;
                c4_b_a21.re = c4_a21.re - c4_r2.re;
                c4_b_a21.im = c4_a21.im - c4_r2.im;
                c4_d7 = c4_abs(chartInstance, c4_c_a12);
                c4_d8 = c4_abs(chartInstance, c4_b_a21);
                if (c4_d7 <= c4_d8) {
                  c4_a21 = c4_a12;
                  c4_rho.re -= c4_a22.re;
                  c4_rho.im -= c4_a22.im;
                } else {
                  c4_rho.re += c4_a22.re;
                  c4_rho.im += c4_a22.im;
                }
              } else {
                c4_eshift.re += c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].re;
                c4_eshift.im += c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilast), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_ilastm1), 1, 3, 2, 0) - 1))
                  - 1].im;
                c4_a21 = c4_eshift;
              }

              c4_b_j = c4_ilastm1;
              c4_u_a = c4_b_j;
              c4_v_a = c4_u_a + 1;
              c4_jp1 = c4_v_a;
              exitg2 = false;
              while ((exitg2 == false) && (c4_b_j > c4_ifirst)) {
                c4_w_a = c4_b_j;
                c4_x_a = c4_w_a - 1;
                c4_jm1 = c4_x_a;
                c4_istart = c4_b_j;
                c4_ctemp.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re - c4_a21.re;
                c4_ctemp.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im - c4_a21.im;
                c4_j_x = c4_ctemp.re;
                c4_k_x = c4_j_x;
                c4_i_y = muDoubleScalarAbs(c4_k_x);
                c4_l_x = c4_ctemp.im;
                c4_m_x = c4_l_x;
                c4_j_y = muDoubleScalarAbs(c4_m_x);
                c4_k_y = c4_i_y + c4_j_y;
                c4_temp = c4_ascale * c4_k_y;
                c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c4_n_x = c4_a22.re;
                c4_o_x = c4_n_x;
                c4_l_y = muDoubleScalarAbs(c4_o_x);
                c4_p_x = c4_a22.im;
                c4_q_x = c4_p_x;
                c4_m_y = muDoubleScalarAbs(c4_q_x);
                c4_n_y = c4_l_y + c4_m_y;
                c4_temp2 = c4_ascale * c4_n_y;
                c4_r_x = c4_temp;
                c4_o_y = c4_temp2;
                c4_tempr = c4_r_x;
                if (c4_o_y > c4_tempr) {
                  c4_tempr = c4_o_y;
                }

                if (c4_tempr < 1.0) {
                  if (c4_tempr != 0.0) {
                    c4_temp /= c4_tempr;
                    c4_temp2 /= c4_tempr;
                  }
                }

                c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .re;
                c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                  (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                  _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) - 1]
                  .im;
                c4_s_x = c4_a22.re;
                c4_t_x = c4_s_x;
                c4_p_y = muDoubleScalarAbs(c4_t_x);
                c4_u_x = c4_a22.im;
                c4_v_x = c4_u_x;
                c4_q_y = muDoubleScalarAbs(c4_v_x);
                c4_r_y = c4_p_y + c4_q_y;
                if (c4_r_y * c4_temp2 <= c4_temp * c4_atol) {
                  c4_goto90 = true;
                  exitg2 = true;
                } else {
                  c4_jp1 = c4_b_j;
                  c4_b_j = c4_jm1;
                }
              }

              if (!c4_goto90) {
                c4_istart = c4_ifirst;
                if (c4_istart == c4_ilastm1) {
                  c4_ctemp = c4_rho;
                } else {
                  c4_ctemp.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 2, 0) - 1))
                    - 1].re - c4_a21.re;
                  c4_ctemp.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 1, 0) + 3 *
                                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 2, 0) - 1))
                    - 1].im - c4_a21.im;
                }

                c4_goto90 = true;
              }
            }

            if (c4_goto90) {
              c4_goto90 = false;
              c4_y_a = c4_istart;
              c4_ab_a = c4_y_a;
              c4_g_c = c4_ab_a;
              c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c4_g_c + 1)), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 2, 0) - 1)) - 1]
                .re;
              c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)(c4_g_c + 1)), 1, 3, 1, 0) + 3 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c4_istart), 1, 3, 2, 0) - 1)) - 1]
                .im;
              c4_b_eml_matlab_zlartg(chartInstance, c4_ctemp, c4_a22, &c4_h_c,
                &c4_s);
              c4_d_c = c4_h_c;
              c4_b_j = c4_istart;
              c4_bb_a = c4_b_j;
              c4_cb_a = c4_bb_a - 1;
              c4_jm1 = c4_cb_a;
              while (c4_b_j < c4_ilast) {
                c4_db_a = c4_b_j;
                c4_eb_a = c4_db_a + 1;
                c4_jp1 = c4_eb_a;
                if (c4_b_j > c4_istart) {
                  c4_e_A.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_e_A.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_f_A.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_f_A.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_eml_matlab_zlartg(chartInstance, c4_e_A, c4_f_A, &c4_i_c,
                                       &c4_s, &c4_a22);
                  c4_d_c = c4_i_c;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) -
                           1)) - 1].re = c4_a22.re;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) -
                           1)) - 1].im = c4_a22.im;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) -
                           1)) - 1].re = c4_dc10.re;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_jp1), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_jm1), 1, 3, 2, 0) -
                           1)) - 1].im = c4_dc10.im;
                }

                c4_j_c = c4_d_c;
                c4_xrow = c4_b_j;
                c4_yrow = c4_jp1;
                c4_jlo = c4_b_j;
                c4_b_jlo = c4_jlo;
                c4_fb_a = c4_b_jlo;
                c4_gb_a = c4_fb_a;
                if (c4_gb_a > 3) {
                  c4_d_overflow = false;
                } else {
                  c4_eml_switch_helper(chartInstance);
                  c4_d_overflow = false;
                }

                if (c4_d_overflow) {
                  c4_check_forloop_overflow_error(chartInstance, c4_d_overflow);
                }

                for (c4_c_j = c4_b_jlo; c4_c_j < 4; c4_c_j++) {
                  c4_d_j = c4_c_j;
                  c4_hb_a = c4_j_c;
                  c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_a12.re = c4_hb_a * c4_a22.re;
                  c4_a12.im = c4_hb_a * c4_a22.im;
                  c4_d_s.re = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re - c4_s.im * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_d_s.im = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im + c4_s.im * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_a21.re = c4_a12.re + c4_d_s.re;
                  c4_a21.im = c4_a12.im + c4_d_s.im;
                  c4_ib_a = c4_j_c;
                  c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_a12.re = c4_ib_a * c4_a22.re;
                  c4_a12.im = c4_ib_a * c4_a22.im;
                  c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].re;
                  c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) - 1)) -
                    1].im;
                  c4_k_a22 = c4_a22;
                  c4_l_a22 = c4_a22;
                  c4_m_a22 = c4_a22;
                  c4_n_a22 = c4_a22;
                  c4_a22.re = c4_s.re * c4_k_a22.re + c4_s.im * c4_l_a22.im;
                  c4_a22.im = c4_s.re * c4_m_a22.im - c4_s.im * c4_n_a22.re;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) -
                             1)) - 1].re = c4_a12.re - c4_a22.re;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_yrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) -
                             1)) - 1].im = c4_a12.im - c4_a22.im;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) -
                             1)) - 1].re = c4_a21.re;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_xrow), 1, 3, 1, 0) +
                        3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_d_j), 1, 3, 2, 0) -
                             1)) - 1].im = c4_a21.im;
                }

                c4_s.re = -c4_s.re;
                c4_s.im = -c4_s.im;
                c4_jb_a = c4_jp1;
                c4_kb_a = c4_jb_a;
                c4_k_c = c4_kb_a;
                c4_w_x = c4_k_c + 1;
                c4_s_y = c4_ilast;
                c4_x_x = c4_w_x;
                if (c4_s_y < c4_x_x) {
                  c4_x_x = c4_s_y;
                }

                c4_l_c = c4_d_c;
                c4_c_xcol = c4_jp1;
                c4_c_ycol = c4_b_j;
                c4_d_ihi = c4_x_x;
                c4_e_ihi = c4_d_ihi;
                c4_i_b = c4_e_ihi;
                c4_j_b = c4_i_b;
                if (1 > c4_j_b) {
                  c4_e_overflow = false;
                } else {
                  c4_eml_switch_helper(chartInstance);
                  c4_e_overflow = (c4_j_b > 2147483646);
                }

                if (c4_e_overflow) {
                  c4_check_forloop_overflow_error(chartInstance, c4_e_overflow);
                }

                for (c4_e_i = 1; c4_e_i <= c4_e_ihi; c4_e_i++) {
                  c4_f_i = c4_e_i;
                  c4_lb_a = c4_l_c;
                  c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_a12.re = c4_lb_a * c4_a22.re;
                  c4_a12.im = c4_lb_a * c4_a22.im;
                  c4_e_s.re = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].re - c4_s.im * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_e_s.im = c4_s.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].im + c4_s.im * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a21.re = c4_a12.re + c4_e_s.re;
                  c4_a21.im = c4_a12.im + c4_e_s.im;
                  c4_mb_a = c4_l_c;
                  c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_a12.re = c4_mb_a * c4_a22.re;
                  c4_a12.im = c4_mb_a * c4_a22.im;
                  c4_a22.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a22.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_c_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_o_a22 = c4_a22;
                  c4_p_a22 = c4_a22;
                  c4_q_a22 = c4_a22;
                  c4_r_a22 = c4_a22;
                  c4_a22.re = c4_s.re * c4_o_a22.re + c4_s.im * c4_p_a22.im;
                  c4_a22.im = c4_s.re * c4_q_a22.im - c4_s.im * c4_r_a22.re;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_c_ycol), 1, 3, 2, 0)
                           - 1)) - 1].re = c4_a12.re - c4_a22.re;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_c_ycol), 1, 3, 2, 0)
                           - 1)) - 1].im = c4_a12.im - c4_a22.im;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_c_xcol), 1, 3, 2, 0)
                           - 1)) - 1].re = c4_a21.re;
                  c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_f_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_c_xcol), 1, 3, 2, 0)
                           - 1)) - 1].im = c4_a21.im;
                }

                c4_m_c = c4_d_c;
                c4_d_xcol = c4_jp1;
                c4_d_ycol = c4_b_j;
                c4_eml_switch_helper(chartInstance);
                for (c4_g_i = 1; c4_g_i < 4; c4_g_i++) {
                  c4_h_i = c4_g_i;
                  c4_nb_a = c4_m_c;
                  c4_a22.re = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a22.im = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_a12.re = c4_nb_a * c4_a22.re;
                  c4_a12.im = c4_nb_a * c4_a22.im;
                  c4_f_s.re = c4_s.re * c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].re - c4_s.im * c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_f_s.im = c4_s.re * c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].im + c4_s.im * c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
                    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0)
                    + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a21.re = c4_a12.re + c4_f_s.re;
                  c4_a21.im = c4_a12.im + c4_f_s.im;
                  c4_ob_a = c4_m_c;
                  c4_a22.re = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a22.im = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_ycol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_a12.re = c4_ob_a * c4_a22.re;
                  c4_a12.im = c4_ob_a * c4_a22.im;
                  c4_a22.re = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].re;
                  c4_a22.im = c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3 *
                                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                    _SFD_INTEGER_CHECK("", (real_T)c4_d_xcol), 1, 3, 2, 0) - 1))
                    - 1].im;
                  c4_s_a22 = c4_a22;
                  c4_t_a22 = c4_a22;
                  c4_u_a22 = c4_a22;
                  c4_v_a22 = c4_a22;
                  c4_a22.re = c4_s.re * c4_s_a22.re + c4_s.im * c4_t_a22.im;
                  c4_a22.im = c4_s.re * c4_u_a22.im - c4_s.im * c4_v_a22.re;
                  c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_d_ycol), 1, 3, 2, 0)
                           - 1)) - 1].re = c4_a12.re - c4_a22.re;
                  c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_d_ycol), 1, 3, 2, 0)
                           - 1)) - 1].im = c4_a12.im - c4_a22.im;
                  c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_d_xcol), 1, 3, 2, 0)
                           - 1)) - 1].re = c4_a21.re;
                  c4_Z[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                         _SFD_INTEGER_CHECK("", (real_T)c4_h_i), 1, 3, 1, 0) + 3
                        * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                          _SFD_INTEGER_CHECK("", (real_T)c4_d_xcol), 1, 3, 2, 0)
                           - 1)) - 1].im = c4_a21.im;
                }

                c4_jm1 = c4_b_j;
                c4_b_j = c4_jp1;
              }
            }

            guard11 = true;
          }

          if (guard11 == true) {
            c4_jiter++;
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
    if (c4_failed) {
      *c4_info = (real_T)c4_ilast;
      c4_b_ilast = c4_ilast;
      c4_k_b = c4_b_ilast;
      c4_l_b = c4_k_b;
      if (1 > c4_l_b) {
        c4_f_overflow = false;
      } else {
        c4_eml_switch_helper(chartInstance);
        c4_f_overflow = (c4_l_b > 2147483646);
      }

      if (c4_f_overflow) {
        c4_check_forloop_overflow_error(chartInstance, c4_f_overflow);
      }

      for (c4_k = 1; c4_k <= c4_b_ilast; c4_k++) {
        c4_b_k = c4_k;
        c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c4_b_k), 1, 3, 1, 0) - 1].re = c4_dc9.re;
        c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c4_b_k), 1, 3, 1, 0) - 1].im = c4_dc9.im;
        c4_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c4_b_k), 1, 3, 1, 0) - 1].re = c4_dc9.re;
        c4_beta1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c4_b_k), 1, 3, 1, 0) - 1].im = c4_dc9.im;
      }

      for (c4_i125 = 0; c4_i125 < 9; c4_i125++) {
        c4_Z[c4_i125] = c4_dc9;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1 == true) {
    c4_pb_a = c4_ilo;
    c4_qb_a = c4_pb_a - 1;
    c4_i126 = c4_qb_a;
    c4_m_b = c4_i126;
    c4_n_b = c4_m_b;
    if (1 > c4_n_b) {
      c4_g_overflow = false;
    } else {
      c4_eml_switch_helper(chartInstance);
      c4_g_overflow = (c4_n_b > 2147483646);
    }

    if (c4_g_overflow) {
      c4_check_forloop_overflow_error(chartInstance, c4_g_overflow);
    }

    for (c4_e_j = 1; c4_e_j <= c4_i126; c4_e_j++) {
      c4_b_j = c4_e_j;
      c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c4_b_j), 1, 3, 1, 0) - 1].re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
        ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
        c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c4_alpha1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c4_b_j), 1, 3, 1, 0) - 1].im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK
        ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c4_b_j), 1, 3, 1, 0) + 3 *
        (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
        c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
    }

    *c4_info = 0.0;
  }
}

static void c4_b_eml_matlab_ztgevc(SFc4_Model_01InstanceStruct *chartInstance,
  creal_T c4_A[9], creal_T c4_V[9])
{
  int32_T c4_i127;
  creal_T c4_work1[3];
  int32_T c4_i128;
  creal_T c4_work2[3];
  int32_T c4_i129;
  real_T c4_rworka[3];
  creal_T c4_ca;
  real_T c4_x;
  real_T c4_b_x;
  real_T c4_y;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_b_y;
  real_T c4_anorm;
  int32_T c4_j;
  real_T c4_b_j;
  real_T c4_d9;
  int32_T c4_i130;
  int32_T c4_i;
  real_T c4_b_i;
  real_T c4_e_x;
  real_T c4_f_x;
  real_T c4_c_y;
  real_T c4_g_x;
  real_T c4_h_x;
  real_T c4_d_y;
  real_T c4_e_y;
  real_T c4_i_x;
  real_T c4_j_x;
  real_T c4_f_y;
  real_T c4_k_x;
  real_T c4_l_x;
  real_T c4_g_y;
  real_T c4_h_y;
  real_T c4_i_y;
  real_T c4_m_x;
  real_T c4_n_x;
  real_T c4_ascale;
  int32_T c4_je;
  real_T c4_b_je;
  real_T c4_ieig;
  real_T c4_o_x;
  real_T c4_p_x;
  real_T c4_j_y;
  real_T c4_q_x;
  real_T c4_r_x;
  real_T c4_k_y;
  real_T c4_l_y;
  real_T c4_s_x;
  real_T c4_t_x;
  real_T c4_temp;
  real_T c4_sbeta;
  real_T c4_a;
  real_T c4_b;
  creal_T c4_salpha;
  real_T c4_acoeff;
  boolean_T c4_b6;
  boolean_T c4_lscalea;
  real_T c4_u_x;
  real_T c4_v_x;
  real_T c4_m_y;
  real_T c4_w_x;
  real_T c4_x_x;
  real_T c4_n_y;
  real_T c4_o_y;
  real_T c4_y_x;
  real_T c4_ab_x;
  real_T c4_p_y;
  real_T c4_bb_x;
  real_T c4_cb_x;
  real_T c4_q_y;
  real_T c4_r_y;
  boolean_T c4_b7;
  boolean_T c4_lscaleb;
  real_T c4_scale;
  real_T c4_db_x;
  real_T c4_eb_x;
  real_T c4_fb_x;
  real_T c4_gb_x;
  real_T c4_s_y;
  real_T c4_hb_x;
  real_T c4_ib_x;
  real_T c4_t_y;
  real_T c4_u_y;
  real_T c4_v_y;
  real_T c4_jb_x;
  real_T c4_kb_x;
  real_T c4_w_y;
  real_T c4_lb_x;
  real_T c4_mb_x;
  real_T c4_x_y;
  real_T c4_y_y;
  real_T c4_nb_x;
  real_T c4_z;
  real_T c4_ob_x;
  real_T c4_ab_y;
  real_T c4_b_a;
  real_T c4_c_a;
  real_T c4_acoefa;
  real_T c4_pb_x;
  real_T c4_qb_x;
  real_T c4_bb_y;
  real_T c4_rb_x;
  real_T c4_sb_x;
  real_T c4_cb_y;
  real_T c4_bcoefa;
  int32_T c4_jr;
  real_T c4_b_jr;
  static creal_T c4_dc11 = { 0.0, 0.0 };

  static creal_T c4_dc12 = { 1.0, 0.0 };

  real_T c4_tb_x;
  real_T c4_db_y;
  real_T c4_dmin;
  real_T c4_d10;
  int32_T c4_i131;
  int32_T c4_c_jr;
  real_T c4_d_a;
  real_T c4_d11;
  int32_T c4_i132;
  int32_T c4_c_j;
  real_T c4_e_a;
  creal_T c4_d;
  real_T c4_ub_x;
  real_T c4_vb_x;
  real_T c4_eb_y;
  real_T c4_wb_x;
  real_T c4_xb_x;
  real_T c4_fb_y;
  real_T c4_gb_y;
  real_T c4_yb_x;
  real_T c4_ac_x;
  real_T c4_hb_y;
  real_T c4_bc_x;
  real_T c4_cc_x;
  real_T c4_ib_y;
  real_T c4_jb_y;
  real_T c4_dc_x;
  real_T c4_ec_x;
  real_T c4_kb_y;
  real_T c4_fc_x;
  real_T c4_gc_x;
  real_T c4_lb_y;
  real_T c4_mb_y;
  real_T c4_hc_x;
  real_T c4_ic_x;
  real_T c4_nb_y;
  real_T c4_jc_x;
  real_T c4_kc_x;
  real_T c4_ob_y;
  real_T c4_pb_y;
  real_T c4_lc_x;
  real_T c4_mc_x;
  real_T c4_qb_y;
  real_T c4_nc_x;
  real_T c4_oc_x;
  real_T c4_rb_y;
  real_T c4_sb_y;
  real_T c4_c_je;
  int32_T c4_i133;
  int32_T c4_d_jr;
  real_T c4_f_a;
  real_T c4_pc_x;
  real_T c4_qc_x;
  real_T c4_tb_y;
  real_T c4_rc_x;
  real_T c4_sc_x;
  real_T c4_ub_y;
  real_T c4_vb_y;
  real_T c4_tc_x;
  real_T c4_uc_x;
  real_T c4_wb_y;
  real_T c4_vc_x;
  real_T c4_wc_x;
  real_T c4_xb_y;
  real_T c4_yb_y;
  real_T c4_d_je;
  int32_T c4_i134;
  int32_T c4_e_jr;
  real_T c4_g_a;
  real_T c4_h_a;
  real_T c4_d12;
  int32_T c4_i135;
  int32_T c4_f_jr;
  creal_T c4_b_ca;
  int32_T c4_g_jr;
  real_T c4_e_je;
  int32_T c4_i136;
  int32_T c4_jc;
  real_T c4_b_jc;
  int32_T c4_h_jr;
  creal_T c4_b_V;
  real_T c4_xc_x;
  real_T c4_yc_x;
  real_T c4_ac_y;
  real_T c4_ad_x;
  real_T c4_bd_x;
  real_T c4_bc_y;
  real_T c4_xmx;
  int32_T c4_i_jr;
  real_T c4_cd_x;
  real_T c4_dd_x;
  real_T c4_cc_y;
  real_T c4_ed_x;
  real_T c4_fd_x;
  real_T c4_dc_y;
  real_T c4_ec_y;
  real_T c4_fc_y;
  int32_T c4_j_jr;
  real_T c4_i_a;
  int32_T c4_k_jr;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  for (c4_i127 = 0; c4_i127 < 3; c4_i127++) {
    c4_work1[c4_i127].re = 0.0;
    c4_work1[c4_i127].im = 0.0;
  }

  for (c4_i128 = 0; c4_i128 < 3; c4_i128++) {
    c4_work2[c4_i128].re = 0.0;
    c4_work2[c4_i128].im = 0.0;
  }

  c4_eps(chartInstance);
  c4_realmin(chartInstance);
  for (c4_i129 = 0; c4_i129 < 3; c4_i129++) {
    c4_rworka[c4_i129] = 0.0;
  }

  c4_ca = c4_A[0];
  c4_x = c4_ca.re;
  c4_b_x = c4_x;
  c4_y = muDoubleScalarAbs(c4_b_x);
  c4_c_x = c4_ca.im;
  c4_d_x = c4_c_x;
  c4_b_y = muDoubleScalarAbs(c4_d_x);
  c4_anorm = c4_y + c4_b_y;
  for (c4_j = 0; c4_j < 2; c4_j++) {
    c4_b_j = 2.0 + (real_T)c4_j;
    c4_d9 = c4_b_j - 1.0;
    c4_i130 = (int32_T)c4_d9;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c4_d9, mxDOUBLE_CLASS, c4_i130);
    for (c4_i = 0; c4_i < c4_i130; c4_i++) {
      c4_b_i = 1.0 + (real_T)c4_i;
      c4_ca.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_i), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c4_ca.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_i), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
      c4_e_x = c4_ca.re;
      c4_f_x = c4_e_x;
      c4_c_y = muDoubleScalarAbs(c4_f_x);
      c4_g_x = c4_ca.im;
      c4_h_x = c4_g_x;
      c4_d_y = muDoubleScalarAbs(c4_h_x);
      c4_e_y = c4_c_y + c4_d_y;
      c4_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_j), 1, 3, 1, 0) - 1] = c4_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)_SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1] + c4_e_y;
    }

    c4_ca.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c4_b_j), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
    c4_ca.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c4_b_j), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
    c4_i_x = c4_ca.re;
    c4_j_x = c4_i_x;
    c4_f_y = muDoubleScalarAbs(c4_j_x);
    c4_k_x = c4_ca.im;
    c4_l_x = c4_k_x;
    c4_g_y = muDoubleScalarAbs(c4_l_x);
    c4_h_y = c4_f_y + c4_g_y;
    c4_i_y = c4_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1] + c4_h_y;
    if (c4_i_y > c4_anorm) {
      c4_anorm = c4_i_y;
    }
  }

  c4_m_x = c4_anorm;
  c4_n_x = c4_m_x;
  if (2.2250738585072014E-308 > c4_n_x) {
    c4_n_x = 2.2250738585072014E-308;
  }

  c4_ascale = 1.0 / c4_n_x;
  for (c4_je = 0; c4_je < 3; c4_je++) {
    c4_b_je = 3.0 + -(real_T)c4_je;
    c4_ieig = c4_b_je;
    c4_ca.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c4_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_b_je), 1, 3, 2, 0) - 1)) - 1].re;
    c4_ca.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c4_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_b_je), 1, 3, 2, 0) - 1)) - 1].im;
    c4_o_x = c4_ca.re;
    c4_p_x = c4_o_x;
    c4_j_y = muDoubleScalarAbs(c4_p_x);
    c4_q_x = c4_ca.im;
    c4_r_x = c4_q_x;
    c4_k_y = muDoubleScalarAbs(c4_r_x);
    c4_l_y = c4_j_y + c4_k_y;
    c4_s_x = c4_l_y * c4_ascale;
    c4_t_x = c4_s_x;
    if (1.0 > c4_t_x) {
      c4_t_x = 1.0;
    }

    c4_temp = 1.0 / c4_t_x;
    c4_sbeta = c4_temp;
    c4_a = c4_temp;
    c4_ca.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c4_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_b_je), 1, 3, 2, 0) - 1)) - 1].re;
    c4_ca.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c4_b_je), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c4_b_je), 1, 3, 2, 0) - 1)) - 1].im;
    c4_ca.re *= c4_a;
    c4_ca.im *= c4_a;
    c4_b = c4_ascale;
    c4_salpha.re = c4_b * c4_ca.re;
    c4_salpha.im = c4_b * c4_ca.im;
    c4_acoeff = c4_sbeta * c4_ascale;
    guard3 = false;
    if (c4_b_abs(chartInstance, c4_sbeta) >= 2.2250738585072014E-308) {
      if (c4_b_abs(chartInstance, c4_acoeff) < 3.0062525400134592E-292) {
        c4_b6 = true;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }

    if (guard3 == true) {
      c4_b6 = false;
    }

    c4_lscalea = c4_b6;
    c4_u_x = c4_salpha.re;
    c4_v_x = c4_u_x;
    c4_m_y = muDoubleScalarAbs(c4_v_x);
    c4_w_x = c4_salpha.im;
    c4_x_x = c4_w_x;
    c4_n_y = muDoubleScalarAbs(c4_x_x);
    c4_o_y = c4_m_y + c4_n_y;
    guard2 = false;
    if (c4_o_y >= 2.2250738585072014E-308) {
      c4_y_x = c4_salpha.re;
      c4_ab_x = c4_y_x;
      c4_p_y = muDoubleScalarAbs(c4_ab_x);
      c4_bb_x = c4_salpha.im;
      c4_cb_x = c4_bb_x;
      c4_q_y = muDoubleScalarAbs(c4_cb_x);
      c4_r_y = c4_p_y + c4_q_y;
      if (c4_r_y < 3.0062525400134592E-292) {
        c4_b7 = true;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2 == true) {
      c4_b7 = false;
    }

    c4_lscaleb = c4_b7;
    c4_scale = 1.0;
    if (c4_lscalea) {
      c4_db_x = c4_anorm;
      c4_eb_x = c4_db_x;
      if (3.3264005158911995E+291 < c4_eb_x) {
        c4_eb_x = 3.3264005158911995E+291;
      }

      c4_scale = 3.0062525400134592E-292 / c4_b_abs(chartInstance, c4_sbeta) *
        c4_eb_x;
    }

    if (c4_lscaleb) {
      c4_fb_x = c4_salpha.re;
      c4_gb_x = c4_fb_x;
      c4_s_y = muDoubleScalarAbs(c4_gb_x);
      c4_hb_x = c4_salpha.im;
      c4_ib_x = c4_hb_x;
      c4_t_y = muDoubleScalarAbs(c4_ib_x);
      c4_u_y = c4_s_y + c4_t_y;
      c4_v_y = 3.0062525400134592E-292 / c4_u_y;
      if (c4_v_y > c4_scale) {
        c4_scale = c4_v_y;
      }
    }

    guard1 = false;
    if (c4_lscalea) {
      guard1 = true;
    } else {
      if (c4_lscaleb) {
        guard1 = true;
      }
    }

    if (guard1 == true) {
      c4_jb_x = c4_salpha.re;
      c4_kb_x = c4_jb_x;
      c4_w_y = muDoubleScalarAbs(c4_kb_x);
      c4_lb_x = c4_salpha.im;
      c4_mb_x = c4_lb_x;
      c4_x_y = muDoubleScalarAbs(c4_mb_x);
      c4_y_y = c4_w_y + c4_x_y;
      c4_nb_x = c4_b_abs(chartInstance, c4_acoeff);
      c4_z = c4_y_y;
      c4_ob_x = c4_nb_x;
      if (1.0 > c4_ob_x) {
        c4_ob_x = 1.0;
      }

      if (c4_z > c4_ob_x) {
        c4_ob_x = c4_z;
      }

      c4_ab_y = 1.0 / (2.2250738585072014E-308 * c4_ob_x);
      if (c4_ab_y < c4_scale) {
        c4_scale = c4_ab_y;
      }

      if (c4_lscalea) {
        c4_acoeff = c4_ascale * (c4_scale * c4_sbeta);
      } else {
        c4_acoeff *= c4_scale;
      }

      if (c4_lscaleb) {
        c4_b_a = c4_scale;
        c4_salpha.re *= c4_b_a;
        c4_salpha.im *= c4_b_a;
      } else {
        c4_c_a = c4_scale;
        c4_salpha.re *= c4_c_a;
        c4_salpha.im *= c4_c_a;
      }
    }

    c4_acoefa = c4_b_abs(chartInstance, c4_acoeff);
    c4_pb_x = c4_salpha.re;
    c4_qb_x = c4_pb_x;
    c4_bb_y = muDoubleScalarAbs(c4_qb_x);
    c4_rb_x = c4_salpha.im;
    c4_sb_x = c4_rb_x;
    c4_cb_y = muDoubleScalarAbs(c4_sb_x);
    c4_bcoefa = c4_bb_y + c4_cb_y;
    for (c4_jr = 0; c4_jr < 3; c4_jr++) {
      c4_b_jr = 1.0 + (real_T)c4_jr;
      c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_jr), 1, 3, 1, 0) - 1].re = c4_dc11.re;
      c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_jr), 1, 3, 1, 0) - 1].im = c4_dc11.im;
    }

    c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c4_b_je), 1, 3, 1, 0) - 1].re = c4_dc12.re;
    c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c4_b_je), 1, 3, 1, 0) - 1].im = c4_dc12.im;
    c4_tb_x = 2.2204460492503131E-16 * c4_acoefa * c4_anorm;
    c4_db_y = 2.2204460492503131E-16 * c4_bcoefa;
    c4_dmin = c4_tb_x;
    if (c4_db_y > c4_dmin) {
      c4_dmin = c4_db_y;
    }

    if (2.2250738585072014E-308 > c4_dmin) {
      c4_dmin = 2.2250738585072014E-308;
    }

    c4_d10 = c4_b_je - 1.0;
    c4_i131 = (int32_T)c4_d10;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c4_d10, mxDOUBLE_CLASS, c4_i131);
    for (c4_c_jr = 0; c4_c_jr < c4_i131; c4_c_jr++) {
      c4_b_jr = 1.0 + (real_T)c4_c_jr;
      c4_d_a = c4_acoeff;
      c4_ca.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_je), 1, 3, 2, 0) - 1)) - 1].re;
      c4_ca.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_je), 1, 3, 2, 0) - 1)) - 1].im;
      c4_ca.re *= c4_d_a;
      c4_ca.im *= c4_d_a;
      c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_jr), 1, 3, 1, 0) - 1].re = c4_ca.re;
      c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_jr), 1, 3, 1, 0) - 1].im = c4_ca.im;
    }

    c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c4_b_je), 1, 3, 1, 0) - 1].re = c4_dc12.re;
    c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c4_b_je), 1, 3, 1, 0) - 1].im = c4_dc12.im;
    c4_d11 = c4_b_je - 1.0;
    c4_i132 = (int32_T)-(1.0 + (-1.0 - c4_d11));
    _SFD_FOR_LOOP_VECTOR_CHECK(c4_d11, -1.0, 1.0, mxDOUBLE_CLASS, c4_i132);
    for (c4_c_j = 0; c4_c_j < c4_i132; c4_c_j++) {
      c4_b_j = c4_d11 + -(real_T)c4_c_j;
      c4_e_a = c4_acoeff;
      c4_ca.re = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 2, 0) - 1)) - 1].re;
      c4_ca.im = c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) + 3 *
                       (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 2, 0) - 1)) - 1].im;
      c4_ca.re *= c4_e_a;
      c4_ca.im *= c4_e_a;
      c4_d.re = c4_ca.re - c4_salpha.re;
      c4_d.im = c4_ca.im - c4_salpha.im;
      c4_ub_x = c4_d.re;
      c4_vb_x = c4_ub_x;
      c4_eb_y = muDoubleScalarAbs(c4_vb_x);
      c4_wb_x = c4_d.im;
      c4_xb_x = c4_wb_x;
      c4_fb_y = muDoubleScalarAbs(c4_xb_x);
      c4_gb_y = c4_eb_y + c4_fb_y;
      if (c4_gb_y <= c4_dmin) {
        c4_d.re = c4_dmin;
        c4_d.im = 0.0;
      }

      c4_yb_x = c4_d.re;
      c4_ac_x = c4_yb_x;
      c4_hb_y = muDoubleScalarAbs(c4_ac_x);
      c4_bc_x = c4_d.im;
      c4_cc_x = c4_bc_x;
      c4_ib_y = muDoubleScalarAbs(c4_cc_x);
      c4_jb_y = c4_hb_y + c4_ib_y;
      if (c4_jb_y < 1.0) {
        c4_ca.re = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].re;
        c4_ca.im = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].im;
        c4_dc_x = c4_ca.re;
        c4_ec_x = c4_dc_x;
        c4_kb_y = muDoubleScalarAbs(c4_ec_x);
        c4_fc_x = c4_ca.im;
        c4_gc_x = c4_fc_x;
        c4_lb_y = muDoubleScalarAbs(c4_gc_x);
        c4_mb_y = c4_kb_y + c4_lb_y;
        c4_hc_x = c4_d.re;
        c4_ic_x = c4_hc_x;
        c4_nb_y = muDoubleScalarAbs(c4_ic_x);
        c4_jc_x = c4_d.im;
        c4_kc_x = c4_jc_x;
        c4_ob_y = muDoubleScalarAbs(c4_kc_x);
        c4_pb_y = c4_nb_y + c4_ob_y;
        if (c4_mb_y >= 1.4980776123852632E+307 * c4_pb_y) {
          c4_ca.re = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].re;
          c4_ca.im = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].im;
          c4_lc_x = c4_ca.re;
          c4_mc_x = c4_lc_x;
          c4_qb_y = muDoubleScalarAbs(c4_mc_x);
          c4_nc_x = c4_ca.im;
          c4_oc_x = c4_nc_x;
          c4_rb_y = muDoubleScalarAbs(c4_oc_x);
          c4_sb_y = c4_qb_y + c4_rb_y;
          c4_temp = 1.0 / c4_sb_y;
          c4_c_je = c4_b_je;
          c4_i133 = (int32_T)c4_c_je;
          _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c4_c_je, mxDOUBLE_CLASS, c4_i133);
          for (c4_d_jr = 0; c4_d_jr < c4_i133; c4_d_jr++) {
            c4_b_jr = 1.0 + (real_T)c4_d_jr;
            c4_f_a = c4_temp;
            c4_ca.re = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].re;
            c4_ca.im = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].im;
            c4_ca.re *= c4_f_a;
            c4_ca.im *= c4_f_a;
            c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
              ("", c4_b_jr), 1, 3, 1, 0) - 1].re = c4_ca.re;
            c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
              ("", c4_b_jr), 1, 3, 1, 0) - 1].im = c4_ca.im;
          }
        }
      }

      c4_ca.re = -c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].re;
      c4_ca.im = -c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].im;
      c4_ca = c4_rdivide(chartInstance, c4_ca, c4_d);
      c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_j), 1, 3, 1, 0) - 1].re = c4_ca.re;
      c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_j), 1, 3, 1, 0) - 1].im = c4_ca.im;
      if (c4_b_j > 1.0) {
        c4_ca.re = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].re;
        c4_ca.im = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].im;
        c4_pc_x = c4_ca.re;
        c4_qc_x = c4_pc_x;
        c4_tb_y = muDoubleScalarAbs(c4_qc_x);
        c4_rc_x = c4_ca.im;
        c4_sc_x = c4_rc_x;
        c4_ub_y = muDoubleScalarAbs(c4_sc_x);
        c4_vb_y = c4_tb_y + c4_ub_y;
        if (c4_vb_y > 1.0) {
          c4_ca.re = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].re;
          c4_ca.im = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].im;
          c4_tc_x = c4_ca.re;
          c4_uc_x = c4_tc_x;
          c4_wb_y = muDoubleScalarAbs(c4_uc_x);
          c4_vc_x = c4_ca.im;
          c4_wc_x = c4_vc_x;
          c4_xb_y = muDoubleScalarAbs(c4_wc_x);
          c4_yb_y = c4_wb_y + c4_xb_y;
          c4_temp = 1.0 / c4_yb_y;
          if (c4_acoefa * c4_rworka[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
               _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1] >=
              1.4980776123852632E+307 * c4_temp) {
            c4_d_je = c4_b_je;
            c4_i134 = (int32_T)c4_d_je;
            _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c4_d_je, mxDOUBLE_CLASS,
              c4_i134);
            for (c4_e_jr = 0; c4_e_jr < c4_i134; c4_e_jr++) {
              c4_b_jr = 1.0 + (real_T)c4_e_jr;
              c4_g_a = c4_temp;
              c4_ca.re = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].re;
              c4_ca.im = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].im;
              c4_ca.re *= c4_g_a;
              c4_ca.im *= c4_g_a;
              c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].re = c4_ca.re;
              c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].im = c4_ca.im;
            }
          }
        }

        c4_h_a = c4_acoeff;
        c4_ca.re = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].re;
        c4_ca.im = c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 1, 0) - 1].im;
        c4_ca.re *= c4_h_a;
        c4_ca.im *= c4_h_a;
        c4_d12 = c4_b_j - 1.0;
        c4_i135 = (int32_T)c4_d12;
        _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c4_d12, mxDOUBLE_CLASS, c4_i135);
        for (c4_f_jr = 0; c4_f_jr < c4_i135; c4_f_jr++) {
          c4_b_jr = 1.0 + (real_T)c4_f_jr;
          c4_b_ca.re = c4_ca.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            c4_b_j), 1, 3, 2, 0) - 1)) - 1].re - c4_ca.im * c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c4_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 2, 0) - 1)) - 1].
            im;
          c4_b_ca.im = c4_ca.re * c4_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
            (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            c4_b_j), 1, 3, 2, 0) - 1)) - 1].im + c4_ca.im * c4_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c4_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c4_b_j), 1, 3, 2, 0) - 1)) - 1].
            re;
          c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
            "", c4_b_jr), 1, 3, 1, 0) - 1].re =
            c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c4_b_jr), 1, 3, 1, 0) - 1].re + c4_b_ca.re;
          c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
            "", c4_b_jr), 1, 3, 1, 0) - 1].im =
            c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
            ("", c4_b_jr), 1, 3, 1, 0) - 1].im + c4_b_ca.im;
        }
      }
    }

    for (c4_g_jr = 0; c4_g_jr < 3; c4_g_jr++) {
      c4_b_jr = 1.0 + (real_T)c4_g_jr;
      c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_jr), 1, 3, 1, 0) - 1].re = c4_dc11.re;
      c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c4_b_jr), 1, 3, 1, 0) - 1].im = c4_dc11.im;
    }

    c4_e_je = c4_b_je;
    c4_i136 = (int32_T)c4_e_je;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c4_e_je, mxDOUBLE_CLASS, c4_i136);
    for (c4_jc = 0; c4_jc < c4_i136; c4_jc++) {
      c4_b_jc = 1.0 + (real_T)c4_jc;
      for (c4_h_jr = 0; c4_h_jr < 3; c4_h_jr++) {
        c4_b_jr = 1.0 + (real_T)c4_h_jr;
        c4_b_V.re = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1)) - 1].re *
          c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c4_b_jc), 1, 3, 1, 0) - 1].re - c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c4_b_jc), 1, 3, 2, 0) - 1)) - 1].im *
          c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c4_b_jc), 1, 3, 1, 0) - 1].im;
        c4_b_V.im = c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_jc), 1, 3, 2, 0) - 1)) - 1].re *
          c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c4_b_jc), 1, 3, 1, 0) - 1].im + c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) + 3 *
          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c4_b_jc), 1, 3, 2, 0) - 1)) - 1].im *
          c4_work1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", c4_b_jc), 1, 3, 1, 0) - 1].re;
        c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c4_b_jr), 1, 3, 1, 0) - 1].re = c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].re +
          c4_b_V.re;
        c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          c4_b_jr), 1, 3, 1, 0) - 1].im = c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].im +
          c4_b_V.im;
      }
    }

    c4_ca = c4_work2[0];
    c4_xc_x = c4_ca.re;
    c4_yc_x = c4_xc_x;
    c4_ac_y = muDoubleScalarAbs(c4_yc_x);
    c4_ad_x = c4_ca.im;
    c4_bd_x = c4_ad_x;
    c4_bc_y = muDoubleScalarAbs(c4_bd_x);
    c4_xmx = c4_ac_y + c4_bc_y;
    for (c4_i_jr = 0; c4_i_jr < 2; c4_i_jr++) {
      c4_b_jr = 2.0 + (real_T)c4_i_jr;
      c4_ca.re = c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].re;
      c4_ca.im = c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].im;
      c4_cd_x = c4_ca.re;
      c4_dd_x = c4_cd_x;
      c4_cc_y = muDoubleScalarAbs(c4_dd_x);
      c4_ed_x = c4_ca.im;
      c4_fd_x = c4_ed_x;
      c4_dc_y = muDoubleScalarAbs(c4_fd_x);
      c4_ec_y = c4_cc_y + c4_dc_y;
      c4_fc_y = c4_ec_y;
      if (c4_fc_y > c4_xmx) {
        c4_xmx = c4_fc_y;
      }
    }

    if (c4_xmx > 2.2250738585072014E-308) {
      c4_temp = 1.0 / c4_xmx;
      for (c4_j_jr = 0; c4_j_jr < 3; c4_j_jr++) {
        c4_b_jr = 1.0 + (real_T)c4_j_jr;
        c4_i_a = c4_temp;
        c4_ca.re = c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].re;
        c4_ca.im = c4_work2[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c4_b_jr), 1, 3, 1, 0) - 1].im;
        c4_ca.re *= c4_i_a;
        c4_ca.im *= c4_i_a;
        c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c4_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c4_ieig), 1, 3, 2, 0) - 1)) - 1]
          .re = c4_ca.re;
        c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c4_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c4_ieig), 1, 3, 2, 0) - 1)) - 1]
          .im = c4_ca.im;
      }
    } else {
      for (c4_k_jr = 0; c4_k_jr < 3; c4_k_jr++) {
        c4_b_jr = 1.0 + (real_T)c4_k_jr;
        c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c4_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c4_ieig), 1, 3, 2, 0) - 1)) - 1]
          .re = c4_dc11.re;
        c4_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                c4_b_jr), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                (int32_T)_SFD_INTEGER_CHECK("", c4_ieig), 1, 3, 2, 0) - 1)) - 1]
          .im = c4_dc11.im;
      }
    }
  }
}

static int32_T c4_div_s32(SFc4_Model_01InstanceStruct *chartInstance, int32_T
  c4_numerator, int32_T c4_denominator)
{
  int32_T c4_quotient;
  uint32_T c4_absNumerator;
  uint32_T c4_absDenominator;
  boolean_T c4_quotientNeedsNegation;
  uint32_T c4_tempAbsQuotient;
  (void)chartInstance;
  if (c4_denominator == 0) {
    if (c4_numerator >= 0) {
      c4_quotient = MAX_int32_T;
    } else {
      c4_quotient = MIN_int32_T;
    }

    _SFD_OVERFLOW_DETECTION(SFDB_DIVIDE_BY_ZERO);
  } else {
    if (c4_numerator >= 0) {
      c4_absNumerator = (uint32_T)c4_numerator;
    } else {
      c4_absNumerator = (uint32_T)-c4_numerator;
    }

    if (c4_denominator >= 0) {
      c4_absDenominator = (uint32_T)c4_denominator;
    } else {
      c4_absDenominator = (uint32_T)-c4_denominator;
    }

    c4_quotientNeedsNegation = (c4_numerator < 0 != c4_denominator < 0);
    c4_tempAbsQuotient = c4_absNumerator / c4_absDenominator;
    if (c4_quotientNeedsNegation) {
      c4_quotient = -(int32_T)c4_tempAbsQuotient;
    } else {
      c4_quotient = (int32_T)c4_tempAbsQuotient;
    }
  }

  return c4_quotient;
}

static void init_dsm_address_info(SFc4_Model_01InstanceStruct *chartInstance)
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

void sf_c4_Model_01_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(409655994U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1935394810U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2779466918U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4058392699U);
}

mxArray *sf_c4_Model_01_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("BSLVtnxiujwcYccoMPm0xD");
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

mxArray *sf_c4_Model_01_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c4_Model_01_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c4_Model_01(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"Phi_r22\",},{M[8],M[0],T\"is_active_c4_Model_01\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c4_Model_01_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc4_Model_01InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc4_Model_01InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_01MachineNumber_,
           4,
           1,
           1,
           0,
           4,
           0,
           0,
           0,
           0,
           0,
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
          _SFD_SET_DATA_PROPS(1,2,0,1,"Phi_r22");
          _SFD_SET_DATA_PROPS(2,1,1,0,"t_delta");
          _SFD_SET_DATA_PROPS(3,1,1,0,"M");
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
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,2,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,376);
        _SFD_CV_INIT_EML_FOR(0,1,0,244,256,372);
        _SFD_CV_INIT_EML_FOR(0,1,1,264,276,364);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)
            c4_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c4_t_delta;
          real_T (*c4_omega)[3];
          real_T (*c4_Phi_r22)[9];
          real_T (*c4_M)[9];
          c4_M = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 2);
          c4_t_delta = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c4_Phi_r22 = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
          c4_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c4_omega);
          _SFD_SET_DATA_VALUE_PTR(1U, *c4_Phi_r22);
          _SFD_SET_DATA_VALUE_PTR(2U, c4_t_delta);
          _SFD_SET_DATA_VALUE_PTR(3U, *c4_M);
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
  return "tGF1mbJPGUsxezq4V9CaKB";
}

static void sf_opaque_initialize_c4_Model_01(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc4_Model_01InstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c4_Model_01((SFc4_Model_01InstanceStruct*) chartInstanceVar);
  initialize_c4_Model_01((SFc4_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c4_Model_01(void *chartInstanceVar)
{
  enable_c4_Model_01((SFc4_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c4_Model_01(void *chartInstanceVar)
{
  disable_c4_Model_01((SFc4_Model_01InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c4_Model_01(void *chartInstanceVar)
{
  sf_gateway_c4_Model_01((SFc4_Model_01InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c4_Model_01(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c4_Model_01((SFc4_Model_01InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c4_Model_01();/* state var info */
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

extern void sf_internal_set_sim_state_c4_Model_01(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c4_Model_01();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c4_Model_01((SFc4_Model_01InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c4_Model_01(SimStruct* S)
{
  return sf_internal_get_sim_state_c4_Model_01(S);
}

static void sf_opaque_set_sim_state_c4_Model_01(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c4_Model_01(S, st);
}

static void sf_opaque_terminate_c4_Model_01(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc4_Model_01InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_01_optimization_info();
    }

    finalize_c4_Model_01((SFc4_Model_01InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc4_Model_01((SFc4_Model_01InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c4_Model_01(SimStruct *S)
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
    initialize_params_c4_Model_01((SFc4_Model_01InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c4_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_01_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,4);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,4,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,4,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,4);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,4,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,4,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,4);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1910036532U));
  ssSetChecksum1(S,(2597560996U));
  ssSetChecksum2(S,(964861719U));
  ssSetChecksum3(S,(1184992409U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c4_Model_01(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c4_Model_01(SimStruct *S)
{
  SFc4_Model_01InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc4_Model_01InstanceStruct *)utMalloc(sizeof
    (SFc4_Model_01InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc4_Model_01InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c4_Model_01;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c4_Model_01;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c4_Model_01;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c4_Model_01;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c4_Model_01;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c4_Model_01;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c4_Model_01;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c4_Model_01;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c4_Model_01;
  chartInstance->chartInfo.mdlStart = mdlStart_c4_Model_01;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c4_Model_01;
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

void c4_Model_01_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c4_Model_01(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c4_Model_01(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c4_Model_01(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c4_Model_01_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

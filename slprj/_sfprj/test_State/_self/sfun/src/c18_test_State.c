/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_State_sfun.h"
#include "c18_test_State.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "test_State_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c18_debug_family_names[8] = { "q_v", "q_0", "cross_q_v",
  "nargin", "nargout", "q", "flag", "CrossTensor" };

static const char * c18_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

/* Function Declarations */
static void initialize_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance);
static void initialize_params_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance);
static void enable_c18_test_State(SFc18_test_StateInstanceStruct *chartInstance);
static void disable_c18_test_State(SFc18_test_StateInstanceStruct *chartInstance);
static void c18_update_debugger_state_c18_test_State
  (SFc18_test_StateInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c18_test_State
  (SFc18_test_StateInstanceStruct *chartInstance);
static void set_sim_state_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_st);
static void finalize_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance);
static void sf_gateway_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance);
static void initSimStructsc18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c18_machineNumber, uint32_T
  c18_chartNumber, uint32_T c18_instanceNumber);
static const mxArray *c18_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData);
static void c18_emlrt_marshallIn(SFc18_test_StateInstanceStruct *chartInstance,
  const mxArray *c18_CrossTensor, const char_T *c18_identifier, real_T c18_y[16]);
static void c18_b_emlrt_marshallIn(SFc18_test_StateInstanceStruct *chartInstance,
  const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId, real_T c18_y[16]);
static void c18_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData);
static const mxArray *c18_b_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData);
static const mxArray *c18_c_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData);
static real_T c18_c_emlrt_marshallIn(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId);
static void c18_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData);
static const mxArray *c18_d_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData);
static void c18_d_emlrt_marshallIn(SFc18_test_StateInstanceStruct *chartInstance,
  const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId, real_T c18_y[9]);
static void c18_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData);
static const mxArray *c18_e_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData);
static void c18_e_emlrt_marshallIn(SFc18_test_StateInstanceStruct *chartInstance,
  const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId, real_T c18_y[3]);
static void c18_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData);
static void c18_info_helper(const mxArray **c18_info);
static const mxArray *c18_emlrt_marshallOut(const char * c18_u);
static const mxArray *c18_b_emlrt_marshallOut(const uint32_T c18_u);
static const mxArray *c18_f_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData);
static int32_T c18_f_emlrt_marshallIn(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId);
static void c18_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData);
static uint8_T c18_g_emlrt_marshallIn(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_b_is_active_c18_test_State, const char_T
  *c18_identifier);
static uint8_T c18_h_emlrt_marshallIn(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId);
static void init_dsm_address_info(SFc18_test_StateInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance)
{
  chartInstance->c18_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c18_is_active_c18_test_State = 0U;
}

static void initialize_params_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c18_test_State(SFc18_test_StateInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c18_test_State(SFc18_test_StateInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c18_update_debugger_state_c18_test_State
  (SFc18_test_StateInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c18_test_State
  (SFc18_test_StateInstanceStruct *chartInstance)
{
  const mxArray *c18_st;
  const mxArray *c18_y = NULL;
  int32_T c18_i0;
  real_T c18_u[16];
  const mxArray *c18_b_y = NULL;
  uint8_T c18_hoistedGlobal;
  uint8_T c18_b_u;
  const mxArray *c18_c_y = NULL;
  real_T (*c18_CrossTensor)[16];
  c18_CrossTensor = (real_T (*)[16])ssGetOutputPortSignal(chartInstance->S, 1);
  c18_st = NULL;
  c18_st = NULL;
  c18_y = NULL;
  sf_mex_assign(&c18_y, sf_mex_createcellmatrix(2, 1), false);
  for (c18_i0 = 0; c18_i0 < 16; c18_i0++) {
    c18_u[c18_i0] = (*c18_CrossTensor)[c18_i0];
  }

  c18_b_y = NULL;
  sf_mex_assign(&c18_b_y, sf_mex_create("y", c18_u, 0, 0U, 1U, 0U, 2, 4, 4),
                false);
  sf_mex_setcell(c18_y, 0, c18_b_y);
  c18_hoistedGlobal = chartInstance->c18_is_active_c18_test_State;
  c18_b_u = c18_hoistedGlobal;
  c18_c_y = NULL;
  sf_mex_assign(&c18_c_y, sf_mex_create("y", &c18_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c18_y, 1, c18_c_y);
  sf_mex_assign(&c18_st, c18_y, false);
  return c18_st;
}

static void set_sim_state_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_st)
{
  const mxArray *c18_u;
  real_T c18_dv0[16];
  int32_T c18_i1;
  real_T (*c18_CrossTensor)[16];
  c18_CrossTensor = (real_T (*)[16])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c18_doneDoubleBufferReInit = true;
  c18_u = sf_mex_dup(c18_st);
  c18_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c18_u, 0)),
                       "CrossTensor", c18_dv0);
  for (c18_i1 = 0; c18_i1 < 16; c18_i1++) {
    (*c18_CrossTensor)[c18_i1] = c18_dv0[c18_i1];
  }

  chartInstance->c18_is_active_c18_test_State = c18_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c18_u, 1)),
     "is_active_c18_test_State");
  sf_mex_destroy(&c18_u);
  c18_update_debugger_state_c18_test_State(chartInstance);
  sf_mex_destroy(&c18_st);
}

static void finalize_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance)
{
  int32_T c18_i2;
  real_T c18_hoistedGlobal;
  int32_T c18_i3;
  real_T c18_q[4];
  real_T c18_flag;
  uint32_T c18_debug_family_var_map[8];
  real_T c18_q_v[3];
  real_T c18_q_0;
  real_T c18_cross_q_v[9];
  real_T c18_nargin = 2.0;
  real_T c18_nargout = 1.0;
  real_T c18_CrossTensor[16];
  int32_T c18_i4;
  int32_T c18_i5;
  int32_T c18_i6;
  real_T c18_v[3];
  uint32_T c18_b_debug_family_var_map[4];
  real_T c18_b_nargin = 1.0;
  real_T c18_b_nargout = 1.0;
  int32_T c18_i7;
  real_T c18_a;
  int32_T c18_i8;
  static real_T c18_b[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c18_y[9];
  int32_T c18_i9;
  int32_T c18_i10;
  int32_T c18_i11;
  int32_T c18_i12;
  int32_T c18_i13;
  int32_T c18_i14;
  int32_T c18_i15;
  real_T c18_b_a;
  int32_T c18_i16;
  int32_T c18_i17;
  int32_T c18_i18;
  int32_T c18_i19;
  int32_T c18_i20;
  int32_T c18_i21;
  int32_T c18_i22;
  int32_T c18_i23;
  int32_T c18_i24;
  int32_T c18_i25;
  real_T *c18_b_flag;
  real_T (*c18_b_CrossTensor)[16];
  real_T (*c18_b_q)[4];
  c18_b_flag = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c18_b_CrossTensor = (real_T (*)[16])ssGetOutputPortSignal(chartInstance->S, 1);
  c18_b_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 2U, chartInstance->c18_sfEvent);
  for (c18_i2 = 0; c18_i2 < 4; c18_i2++) {
    _SFD_DATA_RANGE_CHECK((*c18_b_q)[c18_i2], 0U);
  }

  chartInstance->c18_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 2U, chartInstance->c18_sfEvent);
  c18_hoistedGlobal = *c18_b_flag;
  for (c18_i3 = 0; c18_i3 < 4; c18_i3++) {
    c18_q[c18_i3] = (*c18_b_q)[c18_i3];
  }

  c18_flag = c18_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c18_debug_family_names,
    c18_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c18_q_v, 0U, c18_e_sf_marshallOut,
    c18_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c18_q_0, 1U, c18_b_sf_marshallOut,
    c18_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c18_cross_q_v, 2U, c18_d_sf_marshallOut,
    c18_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c18_nargin, 3U, c18_b_sf_marshallOut,
    c18_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c18_nargout, 4U, c18_b_sf_marshallOut,
    c18_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c18_q, 5U, c18_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c18_flag, 6U, c18_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c18_CrossTensor, 7U, c18_sf_marshallOut,
    c18_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c18_sfEvent, 5);
  for (c18_i4 = 0; c18_i4 < 16; c18_i4++) {
    c18_CrossTensor[c18_i4] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c18_sfEvent, 6);
  for (c18_i5 = 0; c18_i5 < 3; c18_i5++) {
    c18_q_v[c18_i5] = c18_q[c18_i5];
  }

  _SFD_EML_CALL(0U, chartInstance->c18_sfEvent, 7);
  c18_q_0 = c18_q[3];
  _SFD_EML_CALL(0U, chartInstance->c18_sfEvent, 8);
  for (c18_i6 = 0; c18_i6 < 3; c18_i6++) {
    c18_v[c18_i6] = c18_q_v[c18_i6];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c18_b_debug_family_names,
    c18_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c18_b_nargin, 0U, c18_b_sf_marshallOut,
    c18_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c18_b_nargout, 1U, c18_b_sf_marshallOut,
    c18_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c18_v, 2U, c18_e_sf_marshallOut,
    c18_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c18_cross_q_v, 3U, c18_d_sf_marshallOut,
    c18_c_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 2);
  for (c18_i7 = 0; c18_i7 < 9; c18_i7++) {
    c18_cross_q_v[c18_i7] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 3);
  c18_cross_q_v[0] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 4);
  c18_cross_q_v[4] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 5);
  c18_cross_q_v[8] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 6);
  c18_cross_q_v[3] = -c18_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 7);
  c18_cross_q_v[6] = c18_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 8);
  c18_cross_q_v[7] = -c18_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 9);
  c18_cross_q_v[1] = c18_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 10);
  c18_cross_q_v[2] = -c18_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, 11);
  c18_cross_q_v[5] = c18_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c18_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c18_sfEvent, 9);
  switch ((int32_T)_SFD_INTEGER_CHECK("flag", c18_flag)) {
   case 0:
    CV_EML_SWITCH(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c18_sfEvent, 11);
    c18_a = c18_q_0;
    for (c18_i8 = 0; c18_i8 < 9; c18_i8++) {
      c18_y[c18_i8] = c18_a * c18_b[c18_i8];
    }

    c18_i9 = 0;
    c18_i10 = 0;
    for (c18_i11 = 0; c18_i11 < 3; c18_i11++) {
      for (c18_i12 = 0; c18_i12 < 3; c18_i12++) {
        c18_CrossTensor[c18_i12 + c18_i9] = -c18_cross_q_v[c18_i12 + c18_i10] +
          c18_y[c18_i12 + c18_i10];
      }

      c18_i9 += 4;
      c18_i10 += 3;
    }

    for (c18_i13 = 0; c18_i13 < 3; c18_i13++) {
      c18_CrossTensor[c18_i13 + 12] = c18_q_v[c18_i13];
    }

    c18_i14 = 0;
    for (c18_i15 = 0; c18_i15 < 3; c18_i15++) {
      c18_CrossTensor[c18_i14 + 3] = -c18_q_v[c18_i15];
      c18_i14 += 4;
    }

    c18_CrossTensor[15] = c18_q_0;
    break;

   case 1:
    CV_EML_SWITCH(0, 1, 0, 2);
    _SFD_EML_CALL(0U, chartInstance->c18_sfEvent, 13);
    c18_b_a = c18_q_0;
    for (c18_i16 = 0; c18_i16 < 9; c18_i16++) {
      c18_y[c18_i16] = c18_b_a * c18_b[c18_i16];
    }

    c18_i17 = 0;
    c18_i18 = 0;
    for (c18_i19 = 0; c18_i19 < 3; c18_i19++) {
      for (c18_i20 = 0; c18_i20 < 3; c18_i20++) {
        c18_CrossTensor[c18_i20 + c18_i17] = c18_cross_q_v[c18_i20 + c18_i18] +
          c18_y[c18_i20 + c18_i18];
      }

      c18_i17 += 4;
      c18_i18 += 3;
    }

    for (c18_i21 = 0; c18_i21 < 3; c18_i21++) {
      c18_CrossTensor[c18_i21 + 12] = c18_q_v[c18_i21];
    }

    c18_i22 = 0;
    for (c18_i23 = 0; c18_i23 < 3; c18_i23++) {
      c18_CrossTensor[c18_i22 + 3] = -c18_q_v[c18_i23];
      c18_i22 += 4;
    }

    c18_CrossTensor[15] = c18_q_0;
    break;

   default:
    CV_EML_SWITCH(0, 1, 0, 0);
    break;
  }

  _SFD_EML_CALL(0U, chartInstance->c18_sfEvent, -13);
  _SFD_SYMBOL_SCOPE_POP();
  for (c18_i24 = 0; c18_i24 < 16; c18_i24++) {
    (*c18_b_CrossTensor)[c18_i24] = c18_CrossTensor[c18_i24];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 2U, chartInstance->c18_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_test_StateMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c18_i25 = 0; c18_i25 < 16; c18_i25++) {
    _SFD_DATA_RANGE_CHECK((*c18_b_CrossTensor)[c18_i25], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c18_b_flag, 2U);
}

static void initSimStructsc18_test_State(SFc18_test_StateInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c18_machineNumber, uint32_T
  c18_chartNumber, uint32_T c18_instanceNumber)
{
  (void)c18_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c18_chartNumber, c18_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg-2\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c18_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData)
{
  const mxArray *c18_mxArrayOutData = NULL;
  int32_T c18_i26;
  int32_T c18_i27;
  int32_T c18_i28;
  real_T c18_b_inData[16];
  int32_T c18_i29;
  int32_T c18_i30;
  int32_T c18_i31;
  real_T c18_u[16];
  const mxArray *c18_y = NULL;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_mxArrayOutData = NULL;
  c18_i26 = 0;
  for (c18_i27 = 0; c18_i27 < 4; c18_i27++) {
    for (c18_i28 = 0; c18_i28 < 4; c18_i28++) {
      c18_b_inData[c18_i28 + c18_i26] = (*(real_T (*)[16])c18_inData)[c18_i28 +
        c18_i26];
    }

    c18_i26 += 4;
  }

  c18_i29 = 0;
  for (c18_i30 = 0; c18_i30 < 4; c18_i30++) {
    for (c18_i31 = 0; c18_i31 < 4; c18_i31++) {
      c18_u[c18_i31 + c18_i29] = c18_b_inData[c18_i31 + c18_i29];
    }

    c18_i29 += 4;
  }

  c18_y = NULL;
  sf_mex_assign(&c18_y, sf_mex_create("y", c18_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c18_mxArrayOutData, c18_y, false);
  return c18_mxArrayOutData;
}

static void c18_emlrt_marshallIn(SFc18_test_StateInstanceStruct *chartInstance,
  const mxArray *c18_CrossTensor, const char_T *c18_identifier, real_T c18_y[16])
{
  emlrtMsgIdentifier c18_thisId;
  c18_thisId.fIdentifier = c18_identifier;
  c18_thisId.fParent = NULL;
  c18_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c18_CrossTensor), &c18_thisId,
    c18_y);
  sf_mex_destroy(&c18_CrossTensor);
}

static void c18_b_emlrt_marshallIn(SFc18_test_StateInstanceStruct *chartInstance,
  const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId, real_T c18_y[16])
{
  real_T c18_dv1[16];
  int32_T c18_i32;
  (void)chartInstance;
  sf_mex_import(c18_parentId, sf_mex_dup(c18_u), c18_dv1, 1, 0, 0U, 1, 0U, 2, 4,
                4);
  for (c18_i32 = 0; c18_i32 < 16; c18_i32++) {
    c18_y[c18_i32] = c18_dv1[c18_i32];
  }

  sf_mex_destroy(&c18_u);
}

static void c18_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData)
{
  const mxArray *c18_CrossTensor;
  const char_T *c18_identifier;
  emlrtMsgIdentifier c18_thisId;
  real_T c18_y[16];
  int32_T c18_i33;
  int32_T c18_i34;
  int32_T c18_i35;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_CrossTensor = sf_mex_dup(c18_mxArrayInData);
  c18_identifier = c18_varName;
  c18_thisId.fIdentifier = c18_identifier;
  c18_thisId.fParent = NULL;
  c18_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c18_CrossTensor), &c18_thisId,
    c18_y);
  sf_mex_destroy(&c18_CrossTensor);
  c18_i33 = 0;
  for (c18_i34 = 0; c18_i34 < 4; c18_i34++) {
    for (c18_i35 = 0; c18_i35 < 4; c18_i35++) {
      (*(real_T (*)[16])c18_outData)[c18_i35 + c18_i33] = c18_y[c18_i35 +
        c18_i33];
    }

    c18_i33 += 4;
  }

  sf_mex_destroy(&c18_mxArrayInData);
}

static const mxArray *c18_b_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData)
{
  const mxArray *c18_mxArrayOutData = NULL;
  real_T c18_u;
  const mxArray *c18_y = NULL;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_mxArrayOutData = NULL;
  c18_u = *(real_T *)c18_inData;
  c18_y = NULL;
  sf_mex_assign(&c18_y, sf_mex_create("y", &c18_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c18_mxArrayOutData, c18_y, false);
  return c18_mxArrayOutData;
}

static const mxArray *c18_c_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData)
{
  const mxArray *c18_mxArrayOutData = NULL;
  int32_T c18_i36;
  real_T c18_b_inData[4];
  int32_T c18_i37;
  real_T c18_u[4];
  const mxArray *c18_y = NULL;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_mxArrayOutData = NULL;
  for (c18_i36 = 0; c18_i36 < 4; c18_i36++) {
    c18_b_inData[c18_i36] = (*(real_T (*)[4])c18_inData)[c18_i36];
  }

  for (c18_i37 = 0; c18_i37 < 4; c18_i37++) {
    c18_u[c18_i37] = c18_b_inData[c18_i37];
  }

  c18_y = NULL;
  sf_mex_assign(&c18_y, sf_mex_create("y", c18_u, 0, 0U, 1U, 0U, 2, 4, 1), false);
  sf_mex_assign(&c18_mxArrayOutData, c18_y, false);
  return c18_mxArrayOutData;
}

static real_T c18_c_emlrt_marshallIn(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId)
{
  real_T c18_y;
  real_T c18_d0;
  (void)chartInstance;
  sf_mex_import(c18_parentId, sf_mex_dup(c18_u), &c18_d0, 1, 0, 0U, 0, 0U, 0);
  c18_y = c18_d0;
  sf_mex_destroy(&c18_u);
  return c18_y;
}

static void c18_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData)
{
  const mxArray *c18_nargout;
  const char_T *c18_identifier;
  emlrtMsgIdentifier c18_thisId;
  real_T c18_y;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_nargout = sf_mex_dup(c18_mxArrayInData);
  c18_identifier = c18_varName;
  c18_thisId.fIdentifier = c18_identifier;
  c18_thisId.fParent = NULL;
  c18_y = c18_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c18_nargout),
    &c18_thisId);
  sf_mex_destroy(&c18_nargout);
  *(real_T *)c18_outData = c18_y;
  sf_mex_destroy(&c18_mxArrayInData);
}

static const mxArray *c18_d_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData)
{
  const mxArray *c18_mxArrayOutData = NULL;
  int32_T c18_i38;
  int32_T c18_i39;
  int32_T c18_i40;
  real_T c18_b_inData[9];
  int32_T c18_i41;
  int32_T c18_i42;
  int32_T c18_i43;
  real_T c18_u[9];
  const mxArray *c18_y = NULL;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_mxArrayOutData = NULL;
  c18_i38 = 0;
  for (c18_i39 = 0; c18_i39 < 3; c18_i39++) {
    for (c18_i40 = 0; c18_i40 < 3; c18_i40++) {
      c18_b_inData[c18_i40 + c18_i38] = (*(real_T (*)[9])c18_inData)[c18_i40 +
        c18_i38];
    }

    c18_i38 += 3;
  }

  c18_i41 = 0;
  for (c18_i42 = 0; c18_i42 < 3; c18_i42++) {
    for (c18_i43 = 0; c18_i43 < 3; c18_i43++) {
      c18_u[c18_i43 + c18_i41] = c18_b_inData[c18_i43 + c18_i41];
    }

    c18_i41 += 3;
  }

  c18_y = NULL;
  sf_mex_assign(&c18_y, sf_mex_create("y", c18_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c18_mxArrayOutData, c18_y, false);
  return c18_mxArrayOutData;
}

static void c18_d_emlrt_marshallIn(SFc18_test_StateInstanceStruct *chartInstance,
  const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId, real_T c18_y[9])
{
  real_T c18_dv2[9];
  int32_T c18_i44;
  (void)chartInstance;
  sf_mex_import(c18_parentId, sf_mex_dup(c18_u), c18_dv2, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c18_i44 = 0; c18_i44 < 9; c18_i44++) {
    c18_y[c18_i44] = c18_dv2[c18_i44];
  }

  sf_mex_destroy(&c18_u);
}

static void c18_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData)
{
  const mxArray *c18_cross_q_v;
  const char_T *c18_identifier;
  emlrtMsgIdentifier c18_thisId;
  real_T c18_y[9];
  int32_T c18_i45;
  int32_T c18_i46;
  int32_T c18_i47;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_cross_q_v = sf_mex_dup(c18_mxArrayInData);
  c18_identifier = c18_varName;
  c18_thisId.fIdentifier = c18_identifier;
  c18_thisId.fParent = NULL;
  c18_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c18_cross_q_v), &c18_thisId,
    c18_y);
  sf_mex_destroy(&c18_cross_q_v);
  c18_i45 = 0;
  for (c18_i46 = 0; c18_i46 < 3; c18_i46++) {
    for (c18_i47 = 0; c18_i47 < 3; c18_i47++) {
      (*(real_T (*)[9])c18_outData)[c18_i47 + c18_i45] = c18_y[c18_i47 + c18_i45];
    }

    c18_i45 += 3;
  }

  sf_mex_destroy(&c18_mxArrayInData);
}

static const mxArray *c18_e_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData)
{
  const mxArray *c18_mxArrayOutData = NULL;
  int32_T c18_i48;
  real_T c18_b_inData[3];
  int32_T c18_i49;
  real_T c18_u[3];
  const mxArray *c18_y = NULL;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_mxArrayOutData = NULL;
  for (c18_i48 = 0; c18_i48 < 3; c18_i48++) {
    c18_b_inData[c18_i48] = (*(real_T (*)[3])c18_inData)[c18_i48];
  }

  for (c18_i49 = 0; c18_i49 < 3; c18_i49++) {
    c18_u[c18_i49] = c18_b_inData[c18_i49];
  }

  c18_y = NULL;
  sf_mex_assign(&c18_y, sf_mex_create("y", c18_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c18_mxArrayOutData, c18_y, false);
  return c18_mxArrayOutData;
}

static void c18_e_emlrt_marshallIn(SFc18_test_StateInstanceStruct *chartInstance,
  const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId, real_T c18_y[3])
{
  real_T c18_dv3[3];
  int32_T c18_i50;
  (void)chartInstance;
  sf_mex_import(c18_parentId, sf_mex_dup(c18_u), c18_dv3, 1, 0, 0U, 1, 0U, 1, 3);
  for (c18_i50 = 0; c18_i50 < 3; c18_i50++) {
    c18_y[c18_i50] = c18_dv3[c18_i50];
  }

  sf_mex_destroy(&c18_u);
}

static void c18_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData)
{
  const mxArray *c18_q_v;
  const char_T *c18_identifier;
  emlrtMsgIdentifier c18_thisId;
  real_T c18_y[3];
  int32_T c18_i51;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_q_v = sf_mex_dup(c18_mxArrayInData);
  c18_identifier = c18_varName;
  c18_thisId.fIdentifier = c18_identifier;
  c18_thisId.fParent = NULL;
  c18_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c18_q_v), &c18_thisId, c18_y);
  sf_mex_destroy(&c18_q_v);
  for (c18_i51 = 0; c18_i51 < 3; c18_i51++) {
    (*(real_T (*)[3])c18_outData)[c18_i51] = c18_y[c18_i51];
  }

  sf_mex_destroy(&c18_mxArrayInData);
}

const mxArray *sf_c18_test_State_get_eml_resolved_functions_info(void)
{
  const mxArray *c18_nameCaptureInfo = NULL;
  c18_nameCaptureInfo = NULL;
  sf_mex_assign(&c18_nameCaptureInfo, sf_mex_createstruct("structure", 2, 21, 1),
                false);
  c18_info_helper(&c18_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c18_nameCaptureInfo);
  return c18_nameCaptureInfo;
}

static void c18_info_helper(const mxArray **c18_info)
{
  const mxArray *c18_rhs0 = NULL;
  const mxArray *c18_lhs0 = NULL;
  const mxArray *c18_rhs1 = NULL;
  const mxArray *c18_lhs1 = NULL;
  const mxArray *c18_rhs2 = NULL;
  const mxArray *c18_lhs2 = NULL;
  const mxArray *c18_rhs3 = NULL;
  const mxArray *c18_lhs3 = NULL;
  const mxArray *c18_rhs4 = NULL;
  const mxArray *c18_lhs4 = NULL;
  const mxArray *c18_rhs5 = NULL;
  const mxArray *c18_lhs5 = NULL;
  const mxArray *c18_rhs6 = NULL;
  const mxArray *c18_lhs6 = NULL;
  const mxArray *c18_rhs7 = NULL;
  const mxArray *c18_lhs7 = NULL;
  const mxArray *c18_rhs8 = NULL;
  const mxArray *c18_lhs8 = NULL;
  const mxArray *c18_rhs9 = NULL;
  const mxArray *c18_lhs9 = NULL;
  const mxArray *c18_rhs10 = NULL;
  const mxArray *c18_lhs10 = NULL;
  const mxArray *c18_rhs11 = NULL;
  const mxArray *c18_lhs11 = NULL;
  const mxArray *c18_rhs12 = NULL;
  const mxArray *c18_lhs12 = NULL;
  const mxArray *c18_rhs13 = NULL;
  const mxArray *c18_lhs13 = NULL;
  const mxArray *c18_rhs14 = NULL;
  const mxArray *c18_lhs14 = NULL;
  const mxArray *c18_rhs15 = NULL;
  const mxArray *c18_lhs15 = NULL;
  const mxArray *c18_rhs16 = NULL;
  const mxArray *c18_lhs16 = NULL;
  const mxArray *c18_rhs17 = NULL;
  const mxArray *c18_lhs17 = NULL;
  const mxArray *c18_rhs18 = NULL;
  const mxArray *c18_lhs18 = NULL;
  const mxArray *c18_rhs19 = NULL;
  const mxArray *c18_lhs19 = NULL;
  const mxArray *c18_rhs20 = NULL;
  const mxArray *c18_lhs20 = NULL;
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 0);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1450040424U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c18_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "context", "context", 1);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("eye"), "name", "name", 1);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c18_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 2);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c18_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 3);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c18_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 4);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("isinf"), "name", "name", 4);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c18_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 5);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c18_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 6);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("eml_is_integer_class"),
                  "name", "name", 6);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c18_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 7);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("intmax"), "name", "name", 7);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c18_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 8);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c18_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 9);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("intmin"), "name", "name", 9);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c18_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 10);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c18_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 11);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "coder.internal.indexIntRelop"), "name", "name", 11);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c18_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 12);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 12);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c18_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 13);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 13);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c18_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 14);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("intmin"), "name", "name", 14);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c18_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 15);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 15);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c18_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("intmax"), "name", "name", 16);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c18_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 17);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c18_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 18);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("intmax"), "name", "name", 18);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c18_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "context", "context", 19);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 19);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c18_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 20);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 20);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c18_info, c18_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c18_info, c18_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c18_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c18_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c18_info, sf_mex_duplicatearraysafe(&c18_lhs20), "lhs", "lhs",
                  20);
  sf_mex_destroy(&c18_rhs0);
  sf_mex_destroy(&c18_lhs0);
  sf_mex_destroy(&c18_rhs1);
  sf_mex_destroy(&c18_lhs1);
  sf_mex_destroy(&c18_rhs2);
  sf_mex_destroy(&c18_lhs2);
  sf_mex_destroy(&c18_rhs3);
  sf_mex_destroy(&c18_lhs3);
  sf_mex_destroy(&c18_rhs4);
  sf_mex_destroy(&c18_lhs4);
  sf_mex_destroy(&c18_rhs5);
  sf_mex_destroy(&c18_lhs5);
  sf_mex_destroy(&c18_rhs6);
  sf_mex_destroy(&c18_lhs6);
  sf_mex_destroy(&c18_rhs7);
  sf_mex_destroy(&c18_lhs7);
  sf_mex_destroy(&c18_rhs8);
  sf_mex_destroy(&c18_lhs8);
  sf_mex_destroy(&c18_rhs9);
  sf_mex_destroy(&c18_lhs9);
  sf_mex_destroy(&c18_rhs10);
  sf_mex_destroy(&c18_lhs10);
  sf_mex_destroy(&c18_rhs11);
  sf_mex_destroy(&c18_lhs11);
  sf_mex_destroy(&c18_rhs12);
  sf_mex_destroy(&c18_lhs12);
  sf_mex_destroy(&c18_rhs13);
  sf_mex_destroy(&c18_lhs13);
  sf_mex_destroy(&c18_rhs14);
  sf_mex_destroy(&c18_lhs14);
  sf_mex_destroy(&c18_rhs15);
  sf_mex_destroy(&c18_lhs15);
  sf_mex_destroy(&c18_rhs16);
  sf_mex_destroy(&c18_lhs16);
  sf_mex_destroy(&c18_rhs17);
  sf_mex_destroy(&c18_lhs17);
  sf_mex_destroy(&c18_rhs18);
  sf_mex_destroy(&c18_lhs18);
  sf_mex_destroy(&c18_rhs19);
  sf_mex_destroy(&c18_lhs19);
  sf_mex_destroy(&c18_rhs20);
  sf_mex_destroy(&c18_lhs20);
}

static const mxArray *c18_emlrt_marshallOut(const char * c18_u)
{
  const mxArray *c18_y = NULL;
  c18_y = NULL;
  sf_mex_assign(&c18_y, sf_mex_create("y", c18_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c18_u)), false);
  return c18_y;
}

static const mxArray *c18_b_emlrt_marshallOut(const uint32_T c18_u)
{
  const mxArray *c18_y = NULL;
  c18_y = NULL;
  sf_mex_assign(&c18_y, sf_mex_create("y", &c18_u, 7, 0U, 0U, 0U, 0), false);
  return c18_y;
}

static const mxArray *c18_f_sf_marshallOut(void *chartInstanceVoid, void
  *c18_inData)
{
  const mxArray *c18_mxArrayOutData = NULL;
  int32_T c18_u;
  const mxArray *c18_y = NULL;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_mxArrayOutData = NULL;
  c18_u = *(int32_T *)c18_inData;
  c18_y = NULL;
  sf_mex_assign(&c18_y, sf_mex_create("y", &c18_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c18_mxArrayOutData, c18_y, false);
  return c18_mxArrayOutData;
}

static int32_T c18_f_emlrt_marshallIn(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId)
{
  int32_T c18_y;
  int32_T c18_i52;
  (void)chartInstance;
  sf_mex_import(c18_parentId, sf_mex_dup(c18_u), &c18_i52, 1, 6, 0U, 0, 0U, 0);
  c18_y = c18_i52;
  sf_mex_destroy(&c18_u);
  return c18_y;
}

static void c18_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c18_mxArrayInData, const char_T *c18_varName, void *c18_outData)
{
  const mxArray *c18_b_sfEvent;
  const char_T *c18_identifier;
  emlrtMsgIdentifier c18_thisId;
  int32_T c18_y;
  SFc18_test_StateInstanceStruct *chartInstance;
  chartInstance = (SFc18_test_StateInstanceStruct *)chartInstanceVoid;
  c18_b_sfEvent = sf_mex_dup(c18_mxArrayInData);
  c18_identifier = c18_varName;
  c18_thisId.fIdentifier = c18_identifier;
  c18_thisId.fParent = NULL;
  c18_y = c18_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c18_b_sfEvent),
    &c18_thisId);
  sf_mex_destroy(&c18_b_sfEvent);
  *(int32_T *)c18_outData = c18_y;
  sf_mex_destroy(&c18_mxArrayInData);
}

static uint8_T c18_g_emlrt_marshallIn(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_b_is_active_c18_test_State, const char_T
  *c18_identifier)
{
  uint8_T c18_y;
  emlrtMsgIdentifier c18_thisId;
  c18_thisId.fIdentifier = c18_identifier;
  c18_thisId.fParent = NULL;
  c18_y = c18_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c18_b_is_active_c18_test_State), &c18_thisId);
  sf_mex_destroy(&c18_b_is_active_c18_test_State);
  return c18_y;
}

static uint8_T c18_h_emlrt_marshallIn(SFc18_test_StateInstanceStruct
  *chartInstance, const mxArray *c18_u, const emlrtMsgIdentifier *c18_parentId)
{
  uint8_T c18_y;
  uint8_T c18_u0;
  (void)chartInstance;
  sf_mex_import(c18_parentId, sf_mex_dup(c18_u), &c18_u0, 1, 3, 0U, 0, 0U, 0);
  c18_y = c18_u0;
  sf_mex_destroy(&c18_u);
  return c18_y;
}

static void init_dsm_address_info(SFc18_test_StateInstanceStruct *chartInstance)
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

void sf_c18_test_State_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2177456119U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3164679274U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1063462902U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2828855279U);
}

mxArray *sf_c18_test_State_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("dAwER4BQVW7WyQwhwHKjN");
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
      pr[1] = (double)(4);
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

mxArray *sf_c18_test_State_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c18_test_State_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c18_test_State(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"CrossTensor\",},{M[8],M[0],T\"is_active_c18_test_State\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c18_test_State_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc18_test_StateInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc18_test_StateInstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _test_StateMachineNumber_,
           18,
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
        init_script_number_translation(_test_StateMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_test_StateMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _test_StateMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"q");
          _SFD_SET_DATA_PROPS(1,2,0,1,"CrossTensor");
          _SFD_SET_DATA_PROPS(2,1,1,0,"flag");
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
        _SFD_CV_INIT_EML(0,1,1,0,0,0,1,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",59,-1,470);

        {
          static int caseStart[] = { 436, 267, 352 };

          static int caseExprEnd[] = { 445, 273, 358 };

          _SFD_CV_INIT_EML_SWITCH(0,1,0,245,258,466,3,&(caseStart[0]),
            &(caseExprEnd[0]));
        }

        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"fn_VectorToSkewSymmetricTensor",0,-1,433);

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c18_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c18_sf_marshallOut,(MexInFcnForType)
            c18_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c18_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c18_flag;
          real_T (*c18_q)[4];
          real_T (*c18_CrossTensor)[16];
          c18_flag = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c18_CrossTensor = (real_T (*)[16])ssGetOutputPortSignal
            (chartInstance->S, 1);
          c18_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c18_q);
          _SFD_SET_DATA_VALUE_PTR(1U, *c18_CrossTensor);
          _SFD_SET_DATA_VALUE_PTR(2U, c18_flag);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _test_StateMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "CA1Hoy753gVYtifiO7STHH";
}

static void sf_opaque_initialize_c18_test_State(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc18_test_StateInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c18_test_State((SFc18_test_StateInstanceStruct*)
    chartInstanceVar);
  initialize_c18_test_State((SFc18_test_StateInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c18_test_State(void *chartInstanceVar)
{
  enable_c18_test_State((SFc18_test_StateInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c18_test_State(void *chartInstanceVar)
{
  disable_c18_test_State((SFc18_test_StateInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c18_test_State(void *chartInstanceVar)
{
  sf_gateway_c18_test_State((SFc18_test_StateInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c18_test_State(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c18_test_State
    ((SFc18_test_StateInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c18_test_State();/* state var info */
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

extern void sf_internal_set_sim_state_c18_test_State(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c18_test_State();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c18_test_State((SFc18_test_StateInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c18_test_State(SimStruct* S)
{
  return sf_internal_get_sim_state_c18_test_State(S);
}

static void sf_opaque_set_sim_state_c18_test_State(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c18_test_State(S, st);
}

static void sf_opaque_terminate_c18_test_State(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc18_test_StateInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_State_optimization_info();
    }

    finalize_c18_test_State((SFc18_test_StateInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc18_test_State((SFc18_test_StateInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c18_test_State(SimStruct *S)
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
    initialize_params_c18_test_State((SFc18_test_StateInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c18_test_State(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_State_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      18);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,18,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,18,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,18);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,18,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,18,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,18);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2098792550U));
  ssSetChecksum1(S,(1254947437U));
  ssSetChecksum2(S,(27937861U));
  ssSetChecksum3(S,(1584991101U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c18_test_State(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c18_test_State(SimStruct *S)
{
  SFc18_test_StateInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc18_test_StateInstanceStruct *)utMalloc(sizeof
    (SFc18_test_StateInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc18_test_StateInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c18_test_State;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c18_test_State;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c18_test_State;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c18_test_State;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c18_test_State;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c18_test_State;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c18_test_State;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c18_test_State;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c18_test_State;
  chartInstance->chartInfo.mdlStart = mdlStart_c18_test_State;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c18_test_State;
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

void c18_test_State_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c18_test_State(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c18_test_State(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c18_test_State(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c18_test_State_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

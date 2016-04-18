/* Include files */

#include <stddef.h>
#include "blas.h"
#include "untitled1_sfun.h"
#include "c17_untitled1.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "untitled1_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c17_debug_family_names[8] = { "q_v", "q_0", "cross_q_v",
  "nargin", "nargout", "q", "flag", "CrossTensor" };

static const char * c17_b_debug_family_names[4] = { "nargin", "nargout", "v",
  "SkewSymmetricTensor" };

/* Function Declarations */
static void initialize_c17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance);
static void initialize_params_c17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance);
static void enable_c17_untitled1(SFc17_untitled1InstanceStruct *chartInstance);
static void disable_c17_untitled1(SFc17_untitled1InstanceStruct *chartInstance);
static void c17_update_debugger_state_c17_untitled1
  (SFc17_untitled1InstanceStruct *chartInstance);
static const mxArray *get_sim_state_c17_untitled1(SFc17_untitled1InstanceStruct *
  chartInstance);
static void set_sim_state_c17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_st);
static void finalize_c17_untitled1(SFc17_untitled1InstanceStruct *chartInstance);
static void sf_gateway_c17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance);
static void initSimStructsc17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c17_machineNumber, uint32_T
  c17_chartNumber, uint32_T c17_instanceNumber);
static const mxArray *c17_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static void c17_emlrt_marshallIn(SFc17_untitled1InstanceStruct *chartInstance,
  const mxArray *c17_CrossTensor, const char_T *c17_identifier, real_T c17_y[16]);
static void c17_b_emlrt_marshallIn(SFc17_untitled1InstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[16]);
static void c17_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static const mxArray *c17_b_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static const mxArray *c17_c_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static real_T c17_c_emlrt_marshallIn(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId);
static void c17_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static const mxArray *c17_d_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static void c17_d_emlrt_marshallIn(SFc17_untitled1InstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[9]);
static void c17_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static const mxArray *c17_e_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static void c17_e_emlrt_marshallIn(SFc17_untitled1InstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[3]);
static void c17_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static void c17_info_helper(const mxArray **c17_info);
static const mxArray *c17_emlrt_marshallOut(const char * c17_u);
static const mxArray *c17_b_emlrt_marshallOut(const uint32_T c17_u);
static const mxArray *c17_f_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static int32_T c17_f_emlrt_marshallIn(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId);
static void c17_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static uint8_T c17_g_emlrt_marshallIn(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_b_is_active_c17_untitled1, const char_T
  *c17_identifier);
static uint8_T c17_h_emlrt_marshallIn(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId);
static void init_dsm_address_info(SFc17_untitled1InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance)
{
  chartInstance->c17_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c17_is_active_c17_untitled1 = 0U;
}

static void initialize_params_c17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c17_untitled1(SFc17_untitled1InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c17_untitled1(SFc17_untitled1InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c17_update_debugger_state_c17_untitled1
  (SFc17_untitled1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c17_untitled1(SFc17_untitled1InstanceStruct *
  chartInstance)
{
  const mxArray *c17_st;
  const mxArray *c17_y = NULL;
  int32_T c17_i0;
  real_T c17_u[16];
  const mxArray *c17_b_y = NULL;
  uint8_T c17_hoistedGlobal;
  uint8_T c17_b_u;
  const mxArray *c17_c_y = NULL;
  real_T (*c17_CrossTensor)[16];
  c17_CrossTensor = (real_T (*)[16])ssGetOutputPortSignal(chartInstance->S, 1);
  c17_st = NULL;
  c17_st = NULL;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_createcellmatrix(2, 1), false);
  for (c17_i0 = 0; c17_i0 < 16; c17_i0++) {
    c17_u[c17_i0] = (*c17_CrossTensor)[c17_i0];
  }

  c17_b_y = NULL;
  sf_mex_assign(&c17_b_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 2, 4, 4),
                false);
  sf_mex_setcell(c17_y, 0, c17_b_y);
  c17_hoistedGlobal = chartInstance->c17_is_active_c17_untitled1;
  c17_b_u = c17_hoistedGlobal;
  c17_c_y = NULL;
  sf_mex_assign(&c17_c_y, sf_mex_create("y", &c17_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c17_y, 1, c17_c_y);
  sf_mex_assign(&c17_st, c17_y, false);
  return c17_st;
}

static void set_sim_state_c17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_st)
{
  const mxArray *c17_u;
  real_T c17_dv0[16];
  int32_T c17_i1;
  real_T (*c17_CrossTensor)[16];
  c17_CrossTensor = (real_T (*)[16])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c17_doneDoubleBufferReInit = true;
  c17_u = sf_mex_dup(c17_st);
  c17_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c17_u, 0)),
                       "CrossTensor", c17_dv0);
  for (c17_i1 = 0; c17_i1 < 16; c17_i1++) {
    (*c17_CrossTensor)[c17_i1] = c17_dv0[c17_i1];
  }

  chartInstance->c17_is_active_c17_untitled1 = c17_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c17_u, 1)),
     "is_active_c17_untitled1");
  sf_mex_destroy(&c17_u);
  c17_update_debugger_state_c17_untitled1(chartInstance);
  sf_mex_destroy(&c17_st);
}

static void finalize_c17_untitled1(SFc17_untitled1InstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance)
{
  int32_T c17_i2;
  real_T c17_hoistedGlobal;
  int32_T c17_i3;
  real_T c17_q[4];
  real_T c17_flag;
  uint32_T c17_debug_family_var_map[8];
  real_T c17_q_v[3];
  real_T c17_q_0;
  real_T c17_cross_q_v[9];
  real_T c17_nargin = 2.0;
  real_T c17_nargout = 1.0;
  real_T c17_CrossTensor[16];
  int32_T c17_i4;
  int32_T c17_i5;
  int32_T c17_i6;
  real_T c17_v[3];
  uint32_T c17_b_debug_family_var_map[4];
  real_T c17_b_nargin = 1.0;
  real_T c17_b_nargout = 1.0;
  int32_T c17_i7;
  real_T c17_a;
  int32_T c17_i8;
  static real_T c17_b[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c17_y[9];
  int32_T c17_i9;
  int32_T c17_i10;
  int32_T c17_i11;
  int32_T c17_i12;
  int32_T c17_i13;
  int32_T c17_i14;
  int32_T c17_i15;
  real_T c17_b_a;
  int32_T c17_i16;
  int32_T c17_i17;
  int32_T c17_i18;
  int32_T c17_i19;
  int32_T c17_i20;
  int32_T c17_i21;
  int32_T c17_i22;
  int32_T c17_i23;
  int32_T c17_i24;
  int32_T c17_i25;
  real_T *c17_b_flag;
  real_T (*c17_b_CrossTensor)[16];
  real_T (*c17_b_q)[4];
  c17_b_flag = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c17_b_CrossTensor = (real_T (*)[16])ssGetOutputPortSignal(chartInstance->S, 1);
  c17_b_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c17_sfEvent);
  for (c17_i2 = 0; c17_i2 < 4; c17_i2++) {
    _SFD_DATA_RANGE_CHECK((*c17_b_q)[c17_i2], 0U);
  }

  chartInstance->c17_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c17_sfEvent);
  c17_hoistedGlobal = *c17_b_flag;
  for (c17_i3 = 0; c17_i3 < 4; c17_i3++) {
    c17_q[c17_i3] = (*c17_b_q)[c17_i3];
  }

  c17_flag = c17_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c17_debug_family_names,
    c17_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_q_v, 0U, c17_e_sf_marshallOut,
    c17_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_q_0, 1U, c17_b_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_cross_q_v, 2U, c17_d_sf_marshallOut,
    c17_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_nargin, 3U, c17_b_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_nargout, 4U, c17_b_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c17_q, 5U, c17_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c17_flag, 6U, c17_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_CrossTensor, 7U, c17_sf_marshallOut,
    c17_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 5);
  for (c17_i4 = 0; c17_i4 < 16; c17_i4++) {
    c17_CrossTensor[c17_i4] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 6);
  for (c17_i5 = 0; c17_i5 < 3; c17_i5++) {
    c17_q_v[c17_i5] = c17_q[c17_i5];
  }

  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 7);
  c17_q_0 = c17_q[3];
  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 8);
  for (c17_i6 = 0; c17_i6 < 3; c17_i6++) {
    c17_v[c17_i6] = c17_q_v[c17_i6];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c17_b_debug_family_names,
    c17_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_b_nargin, 0U, c17_b_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_b_nargout, 1U, c17_b_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_v, 2U, c17_e_sf_marshallOut,
    c17_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_cross_q_v, 3U, c17_d_sf_marshallOut,
    c17_c_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 2);
  for (c17_i7 = 0; c17_i7 < 9; c17_i7++) {
    c17_cross_q_v[c17_i7] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 3);
  c17_cross_q_v[0] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 4);
  c17_cross_q_v[4] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 5);
  c17_cross_q_v[8] = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 6);
  c17_cross_q_v[3] = -c17_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 7);
  c17_cross_q_v[6] = c17_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 8);
  c17_cross_q_v[7] = -c17_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 9);
  c17_cross_q_v[1] = c17_v[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 10);
  c17_cross_q_v[2] = -c17_v[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, 11);
  c17_cross_q_v[5] = c17_v[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c17_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 9);
  switch ((int32_T)_SFD_INTEGER_CHECK("flag", c17_flag)) {
   case 0:
    CV_EML_SWITCH(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 11);
    c17_a = c17_q_0;
    for (c17_i8 = 0; c17_i8 < 9; c17_i8++) {
      c17_y[c17_i8] = c17_a * c17_b[c17_i8];
    }

    c17_i9 = 0;
    c17_i10 = 0;
    for (c17_i11 = 0; c17_i11 < 3; c17_i11++) {
      for (c17_i12 = 0; c17_i12 < 3; c17_i12++) {
        c17_CrossTensor[c17_i12 + c17_i9] = -c17_cross_q_v[c17_i12 + c17_i10] +
          c17_y[c17_i12 + c17_i10];
      }

      c17_i9 += 4;
      c17_i10 += 3;
    }

    for (c17_i13 = 0; c17_i13 < 3; c17_i13++) {
      c17_CrossTensor[c17_i13 + 12] = c17_q_v[c17_i13];
    }

    c17_i14 = 0;
    for (c17_i15 = 0; c17_i15 < 3; c17_i15++) {
      c17_CrossTensor[c17_i14 + 3] = -c17_q_v[c17_i15];
      c17_i14 += 4;
    }

    c17_CrossTensor[15] = c17_q_0;
    break;

   case 1:
    CV_EML_SWITCH(0, 1, 0, 2);
    _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 13);
    c17_b_a = c17_q_0;
    for (c17_i16 = 0; c17_i16 < 9; c17_i16++) {
      c17_y[c17_i16] = c17_b_a * c17_b[c17_i16];
    }

    c17_i17 = 0;
    c17_i18 = 0;
    for (c17_i19 = 0; c17_i19 < 3; c17_i19++) {
      for (c17_i20 = 0; c17_i20 < 3; c17_i20++) {
        c17_CrossTensor[c17_i20 + c17_i17] = c17_cross_q_v[c17_i20 + c17_i18] +
          c17_y[c17_i20 + c17_i18];
      }

      c17_i17 += 4;
      c17_i18 += 3;
    }

    for (c17_i21 = 0; c17_i21 < 3; c17_i21++) {
      c17_CrossTensor[c17_i21 + 12] = c17_q_v[c17_i21];
    }

    c17_i22 = 0;
    for (c17_i23 = 0; c17_i23 < 3; c17_i23++) {
      c17_CrossTensor[c17_i22 + 3] = -c17_q_v[c17_i23];
      c17_i22 += 4;
    }

    c17_CrossTensor[15] = c17_q_0;
    break;

   default:
    CV_EML_SWITCH(0, 1, 0, 0);
    break;
  }

  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, -13);
  _SFD_SYMBOL_SCOPE_POP();
  for (c17_i24 = 0; c17_i24 < 16; c17_i24++) {
    (*c17_b_CrossTensor)[c17_i24] = c17_CrossTensor[c17_i24];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c17_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_untitled1MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c17_i25 = 0; c17_i25 < 16; c17_i25++) {
    _SFD_DATA_RANGE_CHECK((*c17_b_CrossTensor)[c17_i25], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c17_b_flag, 2U);
}

static void initSimStructsc17_untitled1(SFc17_untitled1InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c17_machineNumber, uint32_T
  c17_chartNumber, uint32_T c17_instanceNumber)
{
  (void)c17_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c17_chartNumber, c17_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Iseberg-2\\Documents\\MATLAB\\Model_01\\fn_VectorToSkewSymmetricTensor.m"));
}

static const mxArray *c17_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_i26;
  int32_T c17_i27;
  int32_T c17_i28;
  real_T c17_b_inData[16];
  int32_T c17_i29;
  int32_T c17_i30;
  int32_T c17_i31;
  real_T c17_u[16];
  const mxArray *c17_y = NULL;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  c17_i26 = 0;
  for (c17_i27 = 0; c17_i27 < 4; c17_i27++) {
    for (c17_i28 = 0; c17_i28 < 4; c17_i28++) {
      c17_b_inData[c17_i28 + c17_i26] = (*(real_T (*)[16])c17_inData)[c17_i28 +
        c17_i26];
    }

    c17_i26 += 4;
  }

  c17_i29 = 0;
  for (c17_i30 = 0; c17_i30 < 4; c17_i30++) {
    for (c17_i31 = 0; c17_i31 < 4; c17_i31++) {
      c17_u[c17_i31 + c17_i29] = c17_b_inData[c17_i31 + c17_i29];
    }

    c17_i29 += 4;
  }

  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static void c17_emlrt_marshallIn(SFc17_untitled1InstanceStruct *chartInstance,
  const mxArray *c17_CrossTensor, const char_T *c17_identifier, real_T c17_y[16])
{
  emlrtMsgIdentifier c17_thisId;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_CrossTensor), &c17_thisId,
    c17_y);
  sf_mex_destroy(&c17_CrossTensor);
}

static void c17_b_emlrt_marshallIn(SFc17_untitled1InstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[16])
{
  real_T c17_dv1[16];
  int32_T c17_i32;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), c17_dv1, 1, 0, 0U, 1, 0U, 2, 4,
                4);
  for (c17_i32 = 0; c17_i32 < 16; c17_i32++) {
    c17_y[c17_i32] = c17_dv1[c17_i32];
  }

  sf_mex_destroy(&c17_u);
}

static void c17_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_CrossTensor;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  real_T c17_y[16];
  int32_T c17_i33;
  int32_T c17_i34;
  int32_T c17_i35;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_CrossTensor = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_CrossTensor), &c17_thisId,
    c17_y);
  sf_mex_destroy(&c17_CrossTensor);
  c17_i33 = 0;
  for (c17_i34 = 0; c17_i34 < 4; c17_i34++) {
    for (c17_i35 = 0; c17_i35 < 4; c17_i35++) {
      (*(real_T (*)[16])c17_outData)[c17_i35 + c17_i33] = c17_y[c17_i35 +
        c17_i33];
    }

    c17_i33 += 4;
  }

  sf_mex_destroy(&c17_mxArrayInData);
}

static const mxArray *c17_b_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  real_T c17_u;
  const mxArray *c17_y = NULL;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  c17_u = *(real_T *)c17_inData;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", &c17_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static const mxArray *c17_c_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_i36;
  real_T c17_b_inData[4];
  int32_T c17_i37;
  real_T c17_u[4];
  const mxArray *c17_y = NULL;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  for (c17_i36 = 0; c17_i36 < 4; c17_i36++) {
    c17_b_inData[c17_i36] = (*(real_T (*)[4])c17_inData)[c17_i36];
  }

  for (c17_i37 = 0; c17_i37 < 4; c17_i37++) {
    c17_u[c17_i37] = c17_b_inData[c17_i37];
  }

  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 2, 4, 1), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static real_T c17_c_emlrt_marshallIn(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId)
{
  real_T c17_y;
  real_T c17_d0;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), &c17_d0, 1, 0, 0U, 0, 0U, 0);
  c17_y = c17_d0;
  sf_mex_destroy(&c17_u);
  return c17_y;
}

static void c17_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_nargout;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  real_T c17_y;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_nargout = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_y = c17_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_nargout),
    &c17_thisId);
  sf_mex_destroy(&c17_nargout);
  *(real_T *)c17_outData = c17_y;
  sf_mex_destroy(&c17_mxArrayInData);
}

static const mxArray *c17_d_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_i38;
  int32_T c17_i39;
  int32_T c17_i40;
  real_T c17_b_inData[9];
  int32_T c17_i41;
  int32_T c17_i42;
  int32_T c17_i43;
  real_T c17_u[9];
  const mxArray *c17_y = NULL;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  c17_i38 = 0;
  for (c17_i39 = 0; c17_i39 < 3; c17_i39++) {
    for (c17_i40 = 0; c17_i40 < 3; c17_i40++) {
      c17_b_inData[c17_i40 + c17_i38] = (*(real_T (*)[9])c17_inData)[c17_i40 +
        c17_i38];
    }

    c17_i38 += 3;
  }

  c17_i41 = 0;
  for (c17_i42 = 0; c17_i42 < 3; c17_i42++) {
    for (c17_i43 = 0; c17_i43 < 3; c17_i43++) {
      c17_u[c17_i43 + c17_i41] = c17_b_inData[c17_i43 + c17_i41];
    }

    c17_i41 += 3;
  }

  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static void c17_d_emlrt_marshallIn(SFc17_untitled1InstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[9])
{
  real_T c17_dv2[9];
  int32_T c17_i44;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), c17_dv2, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c17_i44 = 0; c17_i44 < 9; c17_i44++) {
    c17_y[c17_i44] = c17_dv2[c17_i44];
  }

  sf_mex_destroy(&c17_u);
}

static void c17_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_cross_q_v;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  real_T c17_y[9];
  int32_T c17_i45;
  int32_T c17_i46;
  int32_T c17_i47;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_cross_q_v = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_cross_q_v), &c17_thisId,
    c17_y);
  sf_mex_destroy(&c17_cross_q_v);
  c17_i45 = 0;
  for (c17_i46 = 0; c17_i46 < 3; c17_i46++) {
    for (c17_i47 = 0; c17_i47 < 3; c17_i47++) {
      (*(real_T (*)[9])c17_outData)[c17_i47 + c17_i45] = c17_y[c17_i47 + c17_i45];
    }

    c17_i45 += 3;
  }

  sf_mex_destroy(&c17_mxArrayInData);
}

static const mxArray *c17_e_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_i48;
  real_T c17_b_inData[3];
  int32_T c17_i49;
  real_T c17_u[3];
  const mxArray *c17_y = NULL;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  for (c17_i48 = 0; c17_i48 < 3; c17_i48++) {
    c17_b_inData[c17_i48] = (*(real_T (*)[3])c17_inData)[c17_i48];
  }

  for (c17_i49 = 0; c17_i49 < 3; c17_i49++) {
    c17_u[c17_i49] = c17_b_inData[c17_i49];
  }

  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static void c17_e_emlrt_marshallIn(SFc17_untitled1InstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[3])
{
  real_T c17_dv3[3];
  int32_T c17_i50;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), c17_dv3, 1, 0, 0U, 1, 0U, 1, 3);
  for (c17_i50 = 0; c17_i50 < 3; c17_i50++) {
    c17_y[c17_i50] = c17_dv3[c17_i50];
  }

  sf_mex_destroy(&c17_u);
}

static void c17_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_q_v;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  real_T c17_y[3];
  int32_T c17_i51;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_q_v = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_q_v), &c17_thisId, c17_y);
  sf_mex_destroy(&c17_q_v);
  for (c17_i51 = 0; c17_i51 < 3; c17_i51++) {
    (*(real_T (*)[3])c17_outData)[c17_i51] = c17_y[c17_i51];
  }

  sf_mex_destroy(&c17_mxArrayInData);
}

const mxArray *sf_c17_untitled1_get_eml_resolved_functions_info(void)
{
  const mxArray *c17_nameCaptureInfo = NULL;
  c17_nameCaptureInfo = NULL;
  sf_mex_assign(&c17_nameCaptureInfo, sf_mex_createstruct("structure", 2, 21, 1),
                false);
  c17_info_helper(&c17_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c17_nameCaptureInfo);
  return c17_nameCaptureInfo;
}

static void c17_info_helper(const mxArray **c17_info)
{
  const mxArray *c17_rhs0 = NULL;
  const mxArray *c17_lhs0 = NULL;
  const mxArray *c17_rhs1 = NULL;
  const mxArray *c17_lhs1 = NULL;
  const mxArray *c17_rhs2 = NULL;
  const mxArray *c17_lhs2 = NULL;
  const mxArray *c17_rhs3 = NULL;
  const mxArray *c17_lhs3 = NULL;
  const mxArray *c17_rhs4 = NULL;
  const mxArray *c17_lhs4 = NULL;
  const mxArray *c17_rhs5 = NULL;
  const mxArray *c17_lhs5 = NULL;
  const mxArray *c17_rhs6 = NULL;
  const mxArray *c17_lhs6 = NULL;
  const mxArray *c17_rhs7 = NULL;
  const mxArray *c17_lhs7 = NULL;
  const mxArray *c17_rhs8 = NULL;
  const mxArray *c17_lhs8 = NULL;
  const mxArray *c17_rhs9 = NULL;
  const mxArray *c17_lhs9 = NULL;
  const mxArray *c17_rhs10 = NULL;
  const mxArray *c17_lhs10 = NULL;
  const mxArray *c17_rhs11 = NULL;
  const mxArray *c17_lhs11 = NULL;
  const mxArray *c17_rhs12 = NULL;
  const mxArray *c17_lhs12 = NULL;
  const mxArray *c17_rhs13 = NULL;
  const mxArray *c17_lhs13 = NULL;
  const mxArray *c17_rhs14 = NULL;
  const mxArray *c17_lhs14 = NULL;
  const mxArray *c17_rhs15 = NULL;
  const mxArray *c17_lhs15 = NULL;
  const mxArray *c17_rhs16 = NULL;
  const mxArray *c17_lhs16 = NULL;
  const mxArray *c17_rhs17 = NULL;
  const mxArray *c17_lhs17 = NULL;
  const mxArray *c17_rhs18 = NULL;
  const mxArray *c17_lhs18 = NULL;
  const mxArray *c17_rhs19 = NULL;
  const mxArray *c17_lhs19 = NULL;
  const mxArray *c17_rhs20 = NULL;
  const mxArray *c17_lhs20 = NULL;
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "fn_VectorToSkewSymmetricTensor"), "name", "name", 0);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[E]C:/Users/Iseberg-2/Documents/MATLAB/Model_01/fn_VectorToSkewSymmetricTensor.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1450040424U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c17_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "context", "context", 1);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("eye"), "name", "name", 1);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c17_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 2);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c17_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 3);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c17_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 4);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("isinf"), "name", "name", 4);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c17_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 5);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c17_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 6);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("eml_is_integer_class"),
                  "name", "name", 6);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c17_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 7);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("intmax"), "name", "name", 7);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c17_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 8);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c17_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 9);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("intmin"), "name", "name", 9);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c17_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 10);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c17_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 11);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "coder.internal.indexIntRelop"), "name", "name", 11);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c17_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 12);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 12);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c17_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 13);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 13);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c17_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 14);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("intmin"), "name", "name", 14);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c17_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 15);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 15);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c17_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("intmax"), "name", "name", 16);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c17_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 17);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c17_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 18);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("intmax"), "name", "name", 18);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c17_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "context", "context", 19);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 19);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c17_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 20);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 20);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c17_info, c17_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c17_info, c17_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c17_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c17_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c17_info, sf_mex_duplicatearraysafe(&c17_lhs20), "lhs", "lhs",
                  20);
  sf_mex_destroy(&c17_rhs0);
  sf_mex_destroy(&c17_lhs0);
  sf_mex_destroy(&c17_rhs1);
  sf_mex_destroy(&c17_lhs1);
  sf_mex_destroy(&c17_rhs2);
  sf_mex_destroy(&c17_lhs2);
  sf_mex_destroy(&c17_rhs3);
  sf_mex_destroy(&c17_lhs3);
  sf_mex_destroy(&c17_rhs4);
  sf_mex_destroy(&c17_lhs4);
  sf_mex_destroy(&c17_rhs5);
  sf_mex_destroy(&c17_lhs5);
  sf_mex_destroy(&c17_rhs6);
  sf_mex_destroy(&c17_lhs6);
  sf_mex_destroy(&c17_rhs7);
  sf_mex_destroy(&c17_lhs7);
  sf_mex_destroy(&c17_rhs8);
  sf_mex_destroy(&c17_lhs8);
  sf_mex_destroy(&c17_rhs9);
  sf_mex_destroy(&c17_lhs9);
  sf_mex_destroy(&c17_rhs10);
  sf_mex_destroy(&c17_lhs10);
  sf_mex_destroy(&c17_rhs11);
  sf_mex_destroy(&c17_lhs11);
  sf_mex_destroy(&c17_rhs12);
  sf_mex_destroy(&c17_lhs12);
  sf_mex_destroy(&c17_rhs13);
  sf_mex_destroy(&c17_lhs13);
  sf_mex_destroy(&c17_rhs14);
  sf_mex_destroy(&c17_lhs14);
  sf_mex_destroy(&c17_rhs15);
  sf_mex_destroy(&c17_lhs15);
  sf_mex_destroy(&c17_rhs16);
  sf_mex_destroy(&c17_lhs16);
  sf_mex_destroy(&c17_rhs17);
  sf_mex_destroy(&c17_lhs17);
  sf_mex_destroy(&c17_rhs18);
  sf_mex_destroy(&c17_lhs18);
  sf_mex_destroy(&c17_rhs19);
  sf_mex_destroy(&c17_lhs19);
  sf_mex_destroy(&c17_rhs20);
  sf_mex_destroy(&c17_lhs20);
}

static const mxArray *c17_emlrt_marshallOut(const char * c17_u)
{
  const mxArray *c17_y = NULL;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", c17_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c17_u)), false);
  return c17_y;
}

static const mxArray *c17_b_emlrt_marshallOut(const uint32_T c17_u)
{
  const mxArray *c17_y = NULL;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", &c17_u, 7, 0U, 0U, 0U, 0), false);
  return c17_y;
}

static const mxArray *c17_f_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_u;
  const mxArray *c17_y = NULL;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  c17_u = *(int32_T *)c17_inData;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", &c17_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static int32_T c17_f_emlrt_marshallIn(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId)
{
  int32_T c17_y;
  int32_T c17_i52;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), &c17_i52, 1, 6, 0U, 0, 0U, 0);
  c17_y = c17_i52;
  sf_mex_destroy(&c17_u);
  return c17_y;
}

static void c17_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_b_sfEvent;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  int32_T c17_y;
  SFc17_untitled1InstanceStruct *chartInstance;
  chartInstance = (SFc17_untitled1InstanceStruct *)chartInstanceVoid;
  c17_b_sfEvent = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_y = c17_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_b_sfEvent),
    &c17_thisId);
  sf_mex_destroy(&c17_b_sfEvent);
  *(int32_T *)c17_outData = c17_y;
  sf_mex_destroy(&c17_mxArrayInData);
}

static uint8_T c17_g_emlrt_marshallIn(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_b_is_active_c17_untitled1, const char_T
  *c17_identifier)
{
  uint8_T c17_y;
  emlrtMsgIdentifier c17_thisId;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_y = c17_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c17_b_is_active_c17_untitled1), &c17_thisId);
  sf_mex_destroy(&c17_b_is_active_c17_untitled1);
  return c17_y;
}

static uint8_T c17_h_emlrt_marshallIn(SFc17_untitled1InstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId)
{
  uint8_T c17_y;
  uint8_T c17_u0;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), &c17_u0, 1, 3, 0U, 0, 0U, 0);
  c17_y = c17_u0;
  sf_mex_destroy(&c17_u);
  return c17_y;
}

static void init_dsm_address_info(SFc17_untitled1InstanceStruct *chartInstance)
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

void sf_c17_untitled1_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2177456119U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3164679274U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1063462902U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2828855279U);
}

mxArray *sf_c17_untitled1_get_autoinheritance_info(void)
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

mxArray *sf_c17_untitled1_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c17_untitled1_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c17_untitled1(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"CrossTensor\",},{M[8],M[0],T\"is_active_c17_untitled1\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c17_untitled1_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc17_untitled1InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc17_untitled1InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _untitled1MachineNumber_,
           17,
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
        init_script_number_translation(_untitled1MachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_untitled1MachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _untitled1MachineNumber_,
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
            1.0,0,0,(MexFcnForType)c17_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c17_sf_marshallOut,(MexInFcnForType)
            c17_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c17_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c17_flag;
          real_T (*c17_q)[4];
          real_T (*c17_CrossTensor)[16];
          c17_flag = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c17_CrossTensor = (real_T (*)[16])ssGetOutputPortSignal
            (chartInstance->S, 1);
          c17_q = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c17_q);
          _SFD_SET_DATA_VALUE_PTR(1U, *c17_CrossTensor);
          _SFD_SET_DATA_VALUE_PTR(2U, c17_flag);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _untitled1MachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "CA1Hoy753gVYtifiO7STHH";
}

static void sf_opaque_initialize_c17_untitled1(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc17_untitled1InstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c17_untitled1((SFc17_untitled1InstanceStruct*)
    chartInstanceVar);
  initialize_c17_untitled1((SFc17_untitled1InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c17_untitled1(void *chartInstanceVar)
{
  enable_c17_untitled1((SFc17_untitled1InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c17_untitled1(void *chartInstanceVar)
{
  disable_c17_untitled1((SFc17_untitled1InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c17_untitled1(void *chartInstanceVar)
{
  sf_gateway_c17_untitled1((SFc17_untitled1InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c17_untitled1(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c17_untitled1
    ((SFc17_untitled1InstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c17_untitled1();/* state var info */
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

extern void sf_internal_set_sim_state_c17_untitled1(SimStruct* S, const mxArray *
  st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c17_untitled1();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c17_untitled1((SFc17_untitled1InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c17_untitled1(SimStruct* S)
{
  return sf_internal_get_sim_state_c17_untitled1(S);
}

static void sf_opaque_set_sim_state_c17_untitled1(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c17_untitled1(S, st);
}

static void sf_opaque_terminate_c17_untitled1(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc17_untitled1InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_untitled1_optimization_info();
    }

    finalize_c17_untitled1((SFc17_untitled1InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc17_untitled1((SFc17_untitled1InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c17_untitled1(SimStruct *S)
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
    initialize_params_c17_untitled1((SFc17_untitled1InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c17_untitled1(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_untitled1_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      17);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,17,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,17,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,17);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,17,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,17,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,17);
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

static void mdlRTW_c17_untitled1(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c17_untitled1(SimStruct *S)
{
  SFc17_untitled1InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc17_untitled1InstanceStruct *)utMalloc(sizeof
    (SFc17_untitled1InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc17_untitled1InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c17_untitled1;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c17_untitled1;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c17_untitled1;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c17_untitled1;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c17_untitled1;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c17_untitled1;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c17_untitled1;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c17_untitled1;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c17_untitled1;
  chartInstance->chartInfo.mdlStart = mdlStart_c17_untitled1;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c17_untitled1;
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

void c17_untitled1_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c17_untitled1(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c17_untitled1(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c17_untitled1(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c17_untitled1_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

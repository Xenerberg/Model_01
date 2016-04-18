/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Model_02_sfun.h"
#include "c19_Model_02.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "Model_02_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)
#define c19_IN_NO_ACTIVE_CHILD         ((uint8_T)0U)
#define c19_IN_S                       ((uint8_T)1U)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c19_debug_family_names[2] = { "nargin", "nargout" };

static const char * c19_b_debug_family_names[2] = { "nargin", "nargout" };

static const char * c19_c_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static const char * c19_d_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static boolean_T c19_dataWrittenToVector[1];

/* Function Declarations */
static void initialize_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance);
static void initialize_params_c19_Model_02(SFc19_Model_02InstanceStruct
  *chartInstance);
static void enable_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance);
static void disable_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance);
static void c19_update_debugger_state_c19_Model_02(SFc19_Model_02InstanceStruct *
  chartInstance);
static const mxArray *get_sim_state_c19_Model_02(SFc19_Model_02InstanceStruct
  *chartInstance);
static void set_sim_state_c19_Model_02(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_st);
static void c19_set_sim_state_side_effects_c19_Model_02
  (SFc19_Model_02InstanceStruct *chartInstance);
static void finalize_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance);
static void sf_gateway_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance);
static void initSimStructsc19_Model_02(SFc19_Model_02InstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c19_machineNumber, uint32_T
  c19_chartNumber, uint32_T c19_instanceNumber);
static const mxArray *c19_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData);
static real_T c19_emlrt_marshallIn(SFc19_Model_02InstanceStruct *chartInstance,
  const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId);
static void c19_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData);
static const mxArray *c19_b_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData);
static boolean_T c19_b_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId);
static void c19_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData);
static const mxArray *c19_c_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData);
static int32_T c19_c_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId);
static void c19_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData);
static const mxArray *c19_d_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData);
static uint8_T c19_d_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_b_tp_S, const char_T *c19_identifier);
static uint8_T c19_e_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId);
static void c19_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData);
static const mxArray *c19_f_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_b_setSimStateSideEffectsInfo, const char_T *
  c19_identifier);
static const mxArray *c19_g_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId);
static void init_dsm_address_info(SFc19_Model_02InstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance)
{
  chartInstance->c19_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c19_doSetSimStateSideEffects = 0U;
  chartInstance->c19_setSimStateSideEffectsInfo = NULL;
  chartInstance->c19_tp_S = 0U;
  chartInstance->c19_is_active_c19_Model_02 = 0U;
  chartInstance->c19_is_c19_Model_02 = c19_IN_NO_ACTIVE_CHILD;
}

static void initialize_params_c19_Model_02(SFc19_Model_02InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
  sf_call_output_fcn_enable(chartInstance->S, 0, "event_out", 0);
  if (chartInstance->c19_is_c19_Model_02 == c19_IN_S) {
    sf_call_output_fcn_enable(chartInstance->S, 1, "Integ", 0);
  }
}

static void disable_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
  sf_call_output_fcn_disable(chartInstance->S, 0, "event_out", 0);
  if (chartInstance->c19_is_c19_Model_02 == c19_IN_S) {
    sf_call_output_fcn_disable(chartInstance->S, 1, "Integ", 0);
  }
}

static void c19_update_debugger_state_c19_Model_02(SFc19_Model_02InstanceStruct *
  chartInstance)
{
  uint32_T c19_prevAniVal;
  c19_prevAniVal = _SFD_GET_ANIMATION();
  _SFD_SET_ANIMATION(0U);
  _SFD_SET_HONOR_BREAKPOINTS(0U);
  if (chartInstance->c19_is_active_c19_Model_02 == 1U) {
    _SFD_CC_CALL(CHART_ACTIVE_TAG, 16U, chartInstance->c19_sfEvent);
  }

  if (chartInstance->c19_is_c19_Model_02 == c19_IN_S) {
    _SFD_CS_CALL(STATE_ACTIVE_TAG, 0U, chartInstance->c19_sfEvent);
  } else {
    _SFD_CS_CALL(STATE_INACTIVE_TAG, 0U, chartInstance->c19_sfEvent);
  }

  _SFD_SET_ANIMATION(c19_prevAniVal);
  _SFD_SET_HONOR_BREAKPOINTS(1U);
  _SFD_ANIMATE();
}

static const mxArray *get_sim_state_c19_Model_02(SFc19_Model_02InstanceStruct
  *chartInstance)
{
  const mxArray *c19_st;
  const mxArray *c19_y = NULL;
  uint8_T c19_hoistedGlobal;
  uint8_T c19_u;
  const mxArray *c19_b_y = NULL;
  uint8_T c19_b_hoistedGlobal;
  uint8_T c19_b_u;
  const mxArray *c19_c_y = NULL;
  c19_st = NULL;
  c19_st = NULL;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_createcellmatrix(2, 1), false);
  c19_hoistedGlobal = chartInstance->c19_is_active_c19_Model_02;
  c19_u = c19_hoistedGlobal;
  c19_b_y = NULL;
  sf_mex_assign(&c19_b_y, sf_mex_create("y", &c19_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c19_y, 0, c19_b_y);
  c19_b_hoistedGlobal = chartInstance->c19_is_c19_Model_02;
  c19_b_u = c19_b_hoistedGlobal;
  c19_c_y = NULL;
  sf_mex_assign(&c19_c_y, sf_mex_create("y", &c19_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c19_y, 1, c19_c_y);
  sf_mex_assign(&c19_st, c19_y, false);
  return c19_st;
}

static void set_sim_state_c19_Model_02(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_st)
{
  const mxArray *c19_u;
  c19_u = sf_mex_dup(c19_st);
  chartInstance->c19_is_active_c19_Model_02 = c19_d_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c19_u, 0)),
     "is_active_c19_Model_02");
  chartInstance->c19_is_c19_Model_02 = c19_d_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c19_u, 1)), "is_c19_Model_02");
  sf_mex_assign(&chartInstance->c19_setSimStateSideEffectsInfo,
                c19_f_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c19_u, 2)), "setSimStateSideEffectsInfo"), true);
  sf_mex_destroy(&c19_u);
  chartInstance->c19_doSetSimStateSideEffects = 1U;
  c19_update_debugger_state_c19_Model_02(chartInstance);
  sf_mex_destroy(&c19_st);
}

static void c19_set_sim_state_side_effects_c19_Model_02
  (SFc19_Model_02InstanceStruct *chartInstance)
{
  if (chartInstance->c19_doSetSimStateSideEffects != 0) {
    if (chartInstance->c19_is_c19_Model_02 == c19_IN_S) {
      chartInstance->c19_tp_S = 1U;
      if (sf_mex_sub(chartInstance->c19_setSimStateSideEffectsInfo,
                     "setSimStateSideEffectsInfo", 1, 2) == 0.0) {
        sf_call_output_fcn_enable(chartInstance->S, 1, "Integ", 0);
      }
    } else {
      chartInstance->c19_tp_S = 0U;
      if (sf_mex_sub(chartInstance->c19_setSimStateSideEffectsInfo,
                     "setSimStateSideEffectsInfo", 1, 2) > 0.0) {
        sf_call_output_fcn_disable(chartInstance->S, 1, "Integ", 0);
      }
    }

    chartInstance->c19_doSetSimStateSideEffects = 0U;
  }
}

static void finalize_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance)
{
  sf_mex_destroy(&chartInstance->c19_setSimStateSideEffectsInfo);
}

static void sf_gateway_c19_Model_02(SFc19_Model_02InstanceStruct *chartInstance)
{
  uint32_T c19_debug_family_var_map[2];
  real_T c19_nargin = 0.0;
  real_T c19_nargout = 0.0;
  uint32_T c19_b_debug_family_var_map[3];
  real_T c19_b_nargin = 0.0;
  real_T c19_b_nargout = 1.0;
  boolean_T c19_out;
  real_T c19_c_nargin = 0.0;
  real_T c19_c_nargout = 0.0;
  real_T c19_d_nargin = 0.0;
  real_T c19_d_nargout = 0.0;
  real_T *c19_u;
  c19_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c19_set_sim_state_side_effects_c19_Model_02(chartInstance);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 16U, chartInstance->c19_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c19_u, 0U);
  chartInstance->c19_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 16U, chartInstance->c19_sfEvent);
  if (chartInstance->c19_is_active_c19_Model_02 == 0U) {
    _SFD_CC_CALL(CHART_ENTER_ENTRY_FUNCTION_TAG, 16U, chartInstance->c19_sfEvent);
    chartInstance->c19_is_active_c19_Model_02 = 1U;
    _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 16U, chartInstance->c19_sfEvent);
    _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 0U, chartInstance->c19_sfEvent);
    chartInstance->c19_is_c19_Model_02 = c19_IN_S;
    _SFD_CS_CALL(STATE_ACTIVE_TAG, 0U, chartInstance->c19_sfEvent);
    chartInstance->c19_tp_S = 1U;
    sf_call_output_fcn_enable(chartInstance->S, 1, "Integ", 0);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c19_debug_family_names,
      c19_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_nargin, 0U, c19_sf_marshallOut,
      c19_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_nargout, 1U, c19_sf_marshallOut,
      c19_sf_marshallIn);
    sf_call_output_fcn_call(chartInstance->S, 0, "event_out", 0);
    _SFD_SYMBOL_SCOPE_POP();
  } else {
    _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 1U,
                 chartInstance->c19_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c19_c_debug_family_names,
      c19_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_nargin, 0U, c19_sf_marshallOut,
      c19_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_nargout, 1U, c19_sf_marshallOut,
      c19_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_out, 2U, c19_b_sf_marshallOut,
      c19_b_sf_marshallIn);
    c19_out = CV_EML_IF(1, 0, 0, *c19_u != 0.0);
    _SFD_SYMBOL_SCOPE_POP();
    if (c19_out) {
      _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 1U, chartInstance->c19_sfEvent);
      chartInstance->c19_tp_S = 0U;
      _SFD_CS_CALL(STATE_ENTER_EXIT_FUNCTION_TAG, 0U, chartInstance->c19_sfEvent);
      sf_call_output_fcn_disable(chartInstance->S, 1, "Integ", 0);
      _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c19_sfEvent);
      _SFD_CS_CALL(STATE_INACTIVE_TAG, 0U, chartInstance->c19_sfEvent);
      chartInstance->c19_is_c19_Model_02 = c19_IN_S;
      _SFD_CS_CALL(STATE_ACTIVE_TAG, 0U, chartInstance->c19_sfEvent);
      chartInstance->c19_tp_S = 1U;
      sf_call_output_fcn_enable(chartInstance->S, 1, "Integ", 0);
      _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c19_debug_family_names,
        c19_debug_family_var_map);
      _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_c_nargin, 0U, c19_sf_marshallOut,
        c19_sf_marshallIn);
      _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_c_nargout, 1U,
        c19_sf_marshallOut, c19_sf_marshallIn);
      sf_call_output_fcn_call(chartInstance->S, 0, "event_out", 0);
      _SFD_SYMBOL_SCOPE_POP();
    } else {
      _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 0U,
                   chartInstance->c19_sfEvent);
      _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c19_b_debug_family_names,
        c19_debug_family_var_map);
      _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_d_nargin, 0U, c19_sf_marshallOut,
        c19_sf_marshallIn);
      _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_d_nargout, 1U,
        c19_sf_marshallOut, c19_sf_marshallIn);
      sf_call_output_fcn_call(chartInstance->S, 0, "event_out", 0);
      _SFD_SYMBOL_SCOPE_POP();
    }

    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c19_sfEvent);
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 16U, chartInstance->c19_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Model_02MachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc19_Model_02(SFc19_Model_02InstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c19_machineNumber, uint32_T
  c19_chartNumber, uint32_T c19_instanceNumber)
{
  (void)c19_machineNumber;
  (void)c19_chartNumber;
  (void)c19_instanceNumber;
}

static const mxArray *c19_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData)
{
  const mxArray *c19_mxArrayOutData = NULL;
  real_T c19_u;
  const mxArray *c19_y = NULL;
  SFc19_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc19_Model_02InstanceStruct *)chartInstanceVoid;
  c19_mxArrayOutData = NULL;
  c19_u = *(real_T *)c19_inData;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_create("y", &c19_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c19_mxArrayOutData, c19_y, false);
  return c19_mxArrayOutData;
}

static real_T c19_emlrt_marshallIn(SFc19_Model_02InstanceStruct *chartInstance,
  const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId)
{
  real_T c19_y;
  real_T c19_d0;
  (void)chartInstance;
  sf_mex_import(c19_parentId, sf_mex_dup(c19_u), &c19_d0, 1, 0, 0U, 0, 0U, 0);
  c19_y = c19_d0;
  sf_mex_destroy(&c19_u);
  return c19_y;
}

static void c19_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData)
{
  const mxArray *c19_nargout;
  const char_T *c19_identifier;
  emlrtMsgIdentifier c19_thisId;
  real_T c19_y;
  SFc19_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc19_Model_02InstanceStruct *)chartInstanceVoid;
  c19_nargout = sf_mex_dup(c19_mxArrayInData);
  c19_identifier = c19_varName;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_y = c19_emlrt_marshallIn(chartInstance, sf_mex_dup(c19_nargout),
    &c19_thisId);
  sf_mex_destroy(&c19_nargout);
  *(real_T *)c19_outData = c19_y;
  sf_mex_destroy(&c19_mxArrayInData);
}

static const mxArray *c19_b_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData)
{
  const mxArray *c19_mxArrayOutData = NULL;
  boolean_T c19_u;
  const mxArray *c19_y = NULL;
  SFc19_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc19_Model_02InstanceStruct *)chartInstanceVoid;
  c19_mxArrayOutData = NULL;
  c19_u = *(boolean_T *)c19_inData;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_create("y", &c19_u, 11, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c19_mxArrayOutData, c19_y, false);
  return c19_mxArrayOutData;
}

static boolean_T c19_b_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId)
{
  boolean_T c19_y;
  boolean_T c19_b0;
  (void)chartInstance;
  sf_mex_import(c19_parentId, sf_mex_dup(c19_u), &c19_b0, 1, 11, 0U, 0, 0U, 0);
  c19_y = c19_b0;
  sf_mex_destroy(&c19_u);
  return c19_y;
}

static void c19_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData)
{
  const mxArray *c19_sf_internal_predicateOutput;
  const char_T *c19_identifier;
  emlrtMsgIdentifier c19_thisId;
  boolean_T c19_y;
  SFc19_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc19_Model_02InstanceStruct *)chartInstanceVoid;
  c19_sf_internal_predicateOutput = sf_mex_dup(c19_mxArrayInData);
  c19_identifier = c19_varName;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_y = c19_b_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c19_sf_internal_predicateOutput), &c19_thisId);
  sf_mex_destroy(&c19_sf_internal_predicateOutput);
  *(boolean_T *)c19_outData = c19_y;
  sf_mex_destroy(&c19_mxArrayInData);
}

const mxArray *sf_c19_Model_02_get_eml_resolved_functions_info(void)
{
  const mxArray *c19_nameCaptureInfo = NULL;
  c19_nameCaptureInfo = NULL;
  sf_mex_assign(&c19_nameCaptureInfo, sf_mex_create("nameCaptureInfo", NULL, 0,
    0U, 1U, 0U, 2, 0, 1), false);
  return c19_nameCaptureInfo;
}

static const mxArray *c19_c_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData)
{
  const mxArray *c19_mxArrayOutData = NULL;
  int32_T c19_u;
  const mxArray *c19_y = NULL;
  SFc19_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc19_Model_02InstanceStruct *)chartInstanceVoid;
  c19_mxArrayOutData = NULL;
  c19_u = *(int32_T *)c19_inData;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_create("y", &c19_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c19_mxArrayOutData, c19_y, false);
  return c19_mxArrayOutData;
}

static int32_T c19_c_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId)
{
  int32_T c19_y;
  int32_T c19_i0;
  (void)chartInstance;
  sf_mex_import(c19_parentId, sf_mex_dup(c19_u), &c19_i0, 1, 6, 0U, 0, 0U, 0);
  c19_y = c19_i0;
  sf_mex_destroy(&c19_u);
  return c19_y;
}

static void c19_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData)
{
  const mxArray *c19_b_sfEvent;
  const char_T *c19_identifier;
  emlrtMsgIdentifier c19_thisId;
  int32_T c19_y;
  SFc19_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc19_Model_02InstanceStruct *)chartInstanceVoid;
  c19_b_sfEvent = sf_mex_dup(c19_mxArrayInData);
  c19_identifier = c19_varName;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_y = c19_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c19_b_sfEvent),
    &c19_thisId);
  sf_mex_destroy(&c19_b_sfEvent);
  *(int32_T *)c19_outData = c19_y;
  sf_mex_destroy(&c19_mxArrayInData);
}

static const mxArray *c19_d_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData)
{
  const mxArray *c19_mxArrayOutData = NULL;
  uint8_T c19_u;
  const mxArray *c19_y = NULL;
  SFc19_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc19_Model_02InstanceStruct *)chartInstanceVoid;
  c19_mxArrayOutData = NULL;
  c19_u = *(uint8_T *)c19_inData;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_create("y", &c19_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c19_mxArrayOutData, c19_y, false);
  return c19_mxArrayOutData;
}

static uint8_T c19_d_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_b_tp_S, const char_T *c19_identifier)
{
  uint8_T c19_y;
  emlrtMsgIdentifier c19_thisId;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_y = c19_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c19_b_tp_S),
    &c19_thisId);
  sf_mex_destroy(&c19_b_tp_S);
  return c19_y;
}

static uint8_T c19_e_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId)
{
  uint8_T c19_y;
  uint8_T c19_u0;
  (void)chartInstance;
  sf_mex_import(c19_parentId, sf_mex_dup(c19_u), &c19_u0, 1, 3, 0U, 0, 0U, 0);
  c19_y = c19_u0;
  sf_mex_destroy(&c19_u);
  return c19_y;
}

static void c19_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData)
{
  const mxArray *c19_b_tp_S;
  const char_T *c19_identifier;
  emlrtMsgIdentifier c19_thisId;
  uint8_T c19_y;
  SFc19_Model_02InstanceStruct *chartInstance;
  chartInstance = (SFc19_Model_02InstanceStruct *)chartInstanceVoid;
  c19_b_tp_S = sf_mex_dup(c19_mxArrayInData);
  c19_identifier = c19_varName;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_y = c19_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c19_b_tp_S),
    &c19_thisId);
  sf_mex_destroy(&c19_b_tp_S);
  *(uint8_T *)c19_outData = c19_y;
  sf_mex_destroy(&c19_mxArrayInData);
}

static const mxArray *c19_f_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_b_setSimStateSideEffectsInfo, const char_T *
  c19_identifier)
{
  const mxArray *c19_y = NULL;
  emlrtMsgIdentifier c19_thisId;
  c19_y = NULL;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  sf_mex_assign(&c19_y, c19_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c19_b_setSimStateSideEffectsInfo), &c19_thisId), false);
  sf_mex_destroy(&c19_b_setSimStateSideEffectsInfo);
  return c19_y;
}

static const mxArray *c19_g_emlrt_marshallIn(SFc19_Model_02InstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId)
{
  const mxArray *c19_y = NULL;
  (void)chartInstance;
  (void)c19_parentId;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_duplicatearraysafe(&c19_u), false);
  sf_mex_destroy(&c19_u);
  return c19_y;
}

static void init_dsm_address_info(SFc19_Model_02InstanceStruct *chartInstance)
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

void sf_c19_Model_02_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4219072931U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2511915074U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(539507835U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1695700737U);
}

mxArray *sf_c19_Model_02_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("s0EZFnP2VGeUkpkD0ItnXE");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c19_Model_02_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c19_Model_02_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c19_Model_02(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[8],M[0],T\"is_active_c19_Model_02\",},{M[9],M[0],T\"is_c19_Model_02\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c19_Model_02_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc19_Model_02InstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc19_Model_02InstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Model_02MachineNumber_,
           19,
           2,
           2,
           0,
           5,
           1,
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
            1,
            1,
            1);
          _SFD_SET_DATA_PROPS(0,1,1,0,"u");
          _SFD_SET_DATA_PROPS(1,8,0,0,"");
          _SFD_SET_DATA_PROPS(2,8,0,0,"");
          _SFD_SET_DATA_PROPS(3,8,0,0,"");
          _SFD_SET_DATA_PROPS(4,9,0,0,"");
          _SFD_EVENT_SCOPE(0,2);
          _SFD_STATE_INFO(0,0,0);
          _SFD_STATE_INFO(1,0,2);
          _SFD_CH_SUBSTATE_COUNT(1);
          _SFD_CH_SUBSTATE_DECOMP(0);
          _SFD_CH_SUBSTATE_INDEX(0,0);
          _SFD_ST_SUBSTATE_COUNT(0,0);
        }

        _SFD_CV_INIT_CHART(1,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        {
          _SFD_CV_INIT_STATE(1,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(1,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,0,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML(1,0,0,1,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_IF(1,0,0,1,2,1,2);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_sf_marshallOut,(MexInFcnForType)c19_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_sf_marshallOut,(MexInFcnForType)c19_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_sf_marshallOut,(MexInFcnForType)c19_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_sf_marshallOut,(MexInFcnForType)c19_sf_marshallIn);
        _SFD_SET_DATA_VALUE_PTR(1,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(2,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(3,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(4,(void *)(NULL));

        {
          real_T *c19_u;
          c19_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c19_u);
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
  return "LKwCT1JRUAa5pWDBhuS50F";
}

static void sf_opaque_initialize_c19_Model_02(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc19_Model_02InstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c19_Model_02((SFc19_Model_02InstanceStruct*)
    chartInstanceVar);
  initialize_c19_Model_02((SFc19_Model_02InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c19_Model_02(void *chartInstanceVar)
{
  enable_c19_Model_02((SFc19_Model_02InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c19_Model_02(void *chartInstanceVar)
{
  disable_c19_Model_02((SFc19_Model_02InstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c19_Model_02(void *chartInstanceVar)
{
  sf_gateway_c19_Model_02((SFc19_Model_02InstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c19_Model_02(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c19_Model_02((SFc19_Model_02InstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c19_Model_02();/* state var info */
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

extern void sf_internal_set_sim_state_c19_Model_02(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c19_Model_02();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c19_Model_02((SFc19_Model_02InstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c19_Model_02(SimStruct* S)
{
  return sf_internal_get_sim_state_c19_Model_02(S);
}

static void sf_opaque_set_sim_state_c19_Model_02(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c19_Model_02(S, st);
}

static void sf_opaque_terminate_c19_Model_02(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc19_Model_02InstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Model_02_optimization_info();
    }

    finalize_c19_Model_02((SFc19_Model_02InstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc19_Model_02((SFc19_Model_02InstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c19_Model_02(SimStruct *S)
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
    initialize_params_c19_Model_02((SFc19_Model_02InstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c19_Model_02(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Model_02_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      19);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,19,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,19,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,19);
    sf_mark_output_events_with_multiple_callers(S,sf_get_instance_specialization
      (),infoStruct,19,2);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,19,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=3; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,19);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3216940482U));
  ssSetChecksum1(S,(3972562521U));
  ssSetChecksum2(S,(2059278380U));
  ssSetChecksum3(S,(3258672081U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c19_Model_02(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Stateflow");
  }
}

static void mdlStart_c19_Model_02(SimStruct *S)
{
  SFc19_Model_02InstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc19_Model_02InstanceStruct *)utMalloc(sizeof
    (SFc19_Model_02InstanceStruct));
  memset(chartInstance, 0, sizeof(SFc19_Model_02InstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 0;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c19_Model_02;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c19_Model_02;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c19_Model_02;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c19_Model_02;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c19_Model_02;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c19_Model_02;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c19_Model_02;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c19_Model_02;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c19_Model_02;
  chartInstance->chartInfo.mdlStart = mdlStart_c19_Model_02;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c19_Model_02;
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

void c19_Model_02_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c19_Model_02(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c19_Model_02(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c19_Model_02(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c19_Model_02_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}

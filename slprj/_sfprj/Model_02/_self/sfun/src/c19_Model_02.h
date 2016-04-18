#ifndef __c19_Model_02_h__
#define __c19_Model_02_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc19_Model_02InstanceStruct
#define typedef_SFc19_Model_02InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c19_sfEvent;
  uint8_T c19_tp_S;
  boolean_T c19_isStable;
  uint8_T c19_is_active_c19_Model_02;
  uint8_T c19_is_c19_Model_02;
  uint8_T c19_doSetSimStateSideEffects;
  const mxArray *c19_setSimStateSideEffectsInfo;
} SFc19_Model_02InstanceStruct;

#endif                                 /*typedef_SFc19_Model_02InstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c19_Model_02_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c19_Model_02_get_check_sum(mxArray *plhs[]);
extern void c19_Model_02_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif

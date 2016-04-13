#ifndef __c8_Model_02_h__
#define __c8_Model_02_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc8_Model_02InstanceStruct
#define typedef_SFc8_Model_02InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c8_sfEvent;
  boolean_T c8_isStable;
  boolean_T c8_doneDoubleBufferReInit;
  uint8_T c8_is_active_c8_Model_02;
} SFc8_Model_02InstanceStruct;

#endif                                 /*typedef_SFc8_Model_02InstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c8_Model_02_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c8_Model_02_get_check_sum(mxArray *plhs[]);
extern void c8_Model_02_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif

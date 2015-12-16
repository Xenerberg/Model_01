#ifndef __c17_Model_01_h__
#define __c17_Model_01_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc17_Model_01InstanceStruct
#define typedef_SFc17_Model_01InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c17_sfEvent;
  boolean_T c17_isStable;
  boolean_T c17_doneDoubleBufferReInit;
  uint8_T c17_is_active_c17_Model_01;
} SFc17_Model_01InstanceStruct;

#endif                                 /*typedef_SFc17_Model_01InstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c17_Model_01_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c17_Model_01_get_check_sum(mxArray *plhs[]);
extern void c17_Model_01_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif

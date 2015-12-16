#ifndef __c16_Model_01_h__
#define __c16_Model_01_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc16_Model_01InstanceStruct
#define typedef_SFc16_Model_01InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c16_sfEvent;
  boolean_T c16_isStable;
  boolean_T c16_doneDoubleBufferReInit;
  uint8_T c16_is_active_c16_Model_01;
} SFc16_Model_01InstanceStruct;

#endif                                 /*typedef_SFc16_Model_01InstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c16_Model_01_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c16_Model_01_get_check_sum(mxArray *plhs[]);
extern void c16_Model_01_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif

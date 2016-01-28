#ifndef __c25_Model_01_h__
#define __c25_Model_01_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc25_Model_01InstanceStruct
#define typedef_SFc25_Model_01InstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c25_sfEvent;
  boolean_T c25_isStable;
  boolean_T c25_doneDoubleBufferReInit;
  uint8_T c25_is_active_c25_Model_01;
} SFc25_Model_01InstanceStruct;

#endif                                 /*typedef_SFc25_Model_01InstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c25_Model_01_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c25_Model_01_get_check_sum(mxArray *plhs[]);
extern void c25_Model_01_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif

/* Include files */

#include "Model_02_sfun.h"
#include "Model_02_sfun_debug_macros.h"
#include "c1_Model_02.h"
#include "c2_Model_02.h"
#include "c3_Model_02.h"
#include "c4_Model_02.h"
#include "c5_Model_02.h"
#include "c6_Model_02.h"
#include "c7_Model_02.h"
#include "c8_Model_02.h"
#include "c9_Model_02.h"
#include "c10_Model_02.h"
#include "c11_Model_02.h"
#include "c12_Model_02.h"
#include "c13_Model_02.h"
#include "c14_Model_02.h"
#include "c15_Model_02.h"
#include "c18_Model_02.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _Model_02MachineNumber_;

/* Function Declarations */

/* Function Definitions */
void Model_02_initializer(void)
{
}

void Model_02_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_Model_02_method_dispatcher(SimStruct *simstructPtr, unsigned int
  chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==1) {
    c1_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==2) {
    c2_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==3) {
    c3_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==4) {
    c4_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==5) {
    c5_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==6) {
    c6_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==7) {
    c7_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==8) {
    c8_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==9) {
    c9_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==10) {
    c10_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==11) {
    c11_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==12) {
    c12_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==13) {
    c13_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==14) {
    c14_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==15) {
    c15_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==18) {
    c18_Model_02_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_Model_02_process_check_sum_call( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(288660444U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1573660797U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2981196288U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3208649436U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3852544024U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1877877891U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2928356134U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2408195997U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 1:
        {
          extern void sf_c1_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c1_Model_02_get_check_sum(plhs);
          break;
        }

       case 2:
        {
          extern void sf_c2_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c2_Model_02_get_check_sum(plhs);
          break;
        }

       case 3:
        {
          extern void sf_c3_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c3_Model_02_get_check_sum(plhs);
          break;
        }

       case 4:
        {
          extern void sf_c4_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c4_Model_02_get_check_sum(plhs);
          break;
        }

       case 5:
        {
          extern void sf_c5_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c5_Model_02_get_check_sum(plhs);
          break;
        }

       case 6:
        {
          extern void sf_c6_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c6_Model_02_get_check_sum(plhs);
          break;
        }

       case 7:
        {
          extern void sf_c7_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c7_Model_02_get_check_sum(plhs);
          break;
        }

       case 8:
        {
          extern void sf_c8_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c8_Model_02_get_check_sum(plhs);
          break;
        }

       case 9:
        {
          extern void sf_c9_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c9_Model_02_get_check_sum(plhs);
          break;
        }

       case 10:
        {
          extern void sf_c10_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c10_Model_02_get_check_sum(plhs);
          break;
        }

       case 11:
        {
          extern void sf_c11_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c11_Model_02_get_check_sum(plhs);
          break;
        }

       case 12:
        {
          extern void sf_c12_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c12_Model_02_get_check_sum(plhs);
          break;
        }

       case 13:
        {
          extern void sf_c13_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c13_Model_02_get_check_sum(plhs);
          break;
        }

       case 14:
        {
          extern void sf_c14_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c14_Model_02_get_check_sum(plhs);
          break;
        }

       case 15:
        {
          extern void sf_c15_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c15_Model_02_get_check_sum(plhs);
          break;
        }

       case 18:
        {
          extern void sf_c18_Model_02_get_check_sum(mxArray *plhs[]);
          sf_c18_Model_02_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3031367619U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(4001028638U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3978939492U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(838979348U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3841119648U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2635076831U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(927393108U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1177194732U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Model_02_autoinheritance_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char aiChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;
  mxGetString(prhs[2], aiChksum,sizeof(aiChksum)/sizeof(char));
  aiChksum[(sizeof(aiChksum)/sizeof(char)-1)] = '\0';

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(aiChksum, "s0EZFnP2VGeUkpkD0ItnXE") == 0) {
          extern mxArray *sf_c1_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c1_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 2:
      {
        if (strcmp(aiChksum, "dAwER4BQVW7WyQwhwHKjN") == 0) {
          extern mxArray *sf_c2_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c2_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 3:
      {
        if (strcmp(aiChksum, "d4vouQUsiG0gKw23NeUxLE") == 0) {
          extern mxArray *sf_c3_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c3_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 4:
      {
        if (strcmp(aiChksum, "d4vouQUsiG0gKw23NeUxLE") == 0) {
          extern mxArray *sf_c4_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c4_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 5:
      {
        if (strcmp(aiChksum, "d4vouQUsiG0gKw23NeUxLE") == 0) {
          extern mxArray *sf_c5_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c5_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 6:
      {
        if (strcmp(aiChksum, "OkUz36wxTwoos5y8btPRvG") == 0) {
          extern mxArray *sf_c6_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c6_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 7:
      {
        if (strcmp(aiChksum, "YCUUo4NLtqpoc8IBQ9KNGB") == 0) {
          extern mxArray *sf_c7_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c7_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 8:
      {
        if (strcmp(aiChksum, "yW53vFYTTYynsrcdU0R6wE") == 0) {
          extern mxArray *sf_c8_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c8_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 9:
      {
        if (strcmp(aiChksum, "d4vouQUsiG0gKw23NeUxLE") == 0) {
          extern mxArray *sf_c9_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c9_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 10:
      {
        if (strcmp(aiChksum, "9Su5bOd41H9LAwJwHG1sfD") == 0) {
          extern mxArray *sf_c10_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c10_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 11:
      {
        if (strcmp(aiChksum, "d4vouQUsiG0gKw23NeUxLE") == 0) {
          extern mxArray *sf_c11_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c11_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 12:
      {
        if (strcmp(aiChksum, "d4vouQUsiG0gKw23NeUxLE") == 0) {
          extern mxArray *sf_c12_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c12_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 13:
      {
        if (strcmp(aiChksum, "d4vouQUsiG0gKw23NeUxLE") == 0) {
          extern mxArray *sf_c13_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c13_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 14:
      {
        if (strcmp(aiChksum, "UvElzSK2mfRURATcmJQU0B") == 0) {
          extern mxArray *sf_c14_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c14_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 15:
      {
        if (strcmp(aiChksum, "d4vouQUsiG0gKw23NeUxLE") == 0) {
          extern mxArray *sf_c15_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c15_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 18:
      {
        if (strcmp(aiChksum, "dAwER4BQVW7WyQwhwHKjN") == 0) {
          extern mxArray *sf_c18_Model_02_get_autoinheritance_info(void);
          plhs[0] = sf_c18_Model_02_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Model_02_get_eml_resolved_functions_info( int nlhs, mxArray *
  plhs[], int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        extern const mxArray *sf_c1_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c1_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 2:
      {
        extern const mxArray *sf_c2_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 3:
      {
        extern const mxArray *sf_c3_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c3_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 4:
      {
        extern const mxArray *sf_c4_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c4_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 5:
      {
        extern const mxArray *sf_c5_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c5_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 6:
      {
        extern const mxArray *sf_c6_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c6_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 7:
      {
        extern const mxArray *sf_c7_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c7_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 8:
      {
        extern const mxArray *sf_c8_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c8_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 9:
      {
        extern const mxArray *sf_c9_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c9_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 10:
      {
        extern const mxArray *sf_c10_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c10_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 11:
      {
        extern const mxArray *sf_c11_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c11_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 12:
      {
        extern const mxArray *sf_c12_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c12_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 13:
      {
        extern const mxArray *sf_c13_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c13_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 14:
      {
        extern const mxArray *sf_c14_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c14_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 15:
      {
        extern const mxArray *sf_c15_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c15_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 18:
      {
        extern const mxArray *sf_c18_Model_02_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c18_Model_02_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Model_02_third_party_uses_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the third_party_uses_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_third_party_uses_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "LKwCT1JRUAa5pWDBhuS50F") == 0) {
          extern mxArray *sf_c1_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c1_Model_02_third_party_uses_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "CA1Hoy753gVYtifiO7STHH") == 0) {
          extern mxArray *sf_c2_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c2_Model_02_third_party_uses_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c3_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c3_Model_02_third_party_uses_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c4_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c4_Model_02_third_party_uses_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c5_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c5_Model_02_third_party_uses_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "2bLvARcyFLfB8CMWpJ0uOE") == 0) {
          extern mxArray *sf_c6_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c6_Model_02_third_party_uses_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "EforzPMyvPho6YHhxHakWC") == 0) {
          extern mxArray *sf_c7_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c7_Model_02_third_party_uses_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "kSiq8AVWPXEug3LVXJKZGH") == 0) {
          extern mxArray *sf_c8_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c8_Model_02_third_party_uses_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c9_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c9_Model_02_third_party_uses_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "cIgJO9Fg2ZUlYMATVlKSqD") == 0) {
          extern mxArray *sf_c10_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c10_Model_02_third_party_uses_info();
          break;
        }
      }

     case 11:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c11_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c11_Model_02_third_party_uses_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c12_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c12_Model_02_third_party_uses_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c13_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c13_Model_02_third_party_uses_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "qASixlSJ5ukrAvn9UDoRYE") == 0) {
          extern mxArray *sf_c14_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c14_Model_02_third_party_uses_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c15_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c15_Model_02_third_party_uses_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "CA1Hoy753gVYtifiO7STHH") == 0) {
          extern mxArray *sf_c18_Model_02_third_party_uses_info(void);
          plhs[0] = sf_c18_Model_02_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_Model_02_updateBuildInfo_args_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the updateBuildInfo_args_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_updateBuildInfo_args_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "LKwCT1JRUAa5pWDBhuS50F") == 0) {
          extern mxArray *sf_c1_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c1_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "CA1Hoy753gVYtifiO7STHH") == 0) {
          extern mxArray *sf_c2_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c3_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c3_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c4_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c4_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c5_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c5_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "2bLvARcyFLfB8CMWpJ0uOE") == 0) {
          extern mxArray *sf_c6_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c6_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "EforzPMyvPho6YHhxHakWC") == 0) {
          extern mxArray *sf_c7_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c7_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "kSiq8AVWPXEug3LVXJKZGH") == 0) {
          extern mxArray *sf_c8_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c8_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c9_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c9_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "cIgJO9Fg2ZUlYMATVlKSqD") == 0) {
          extern mxArray *sf_c10_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c10_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 11:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c11_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c11_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c12_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c12_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c13_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c13_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "qASixlSJ5ukrAvn9UDoRYE") == 0) {
          extern mxArray *sf_c14_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c14_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "cIBiYNm46H62nSvGIWwbME") == 0) {
          extern mxArray *sf_c15_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c15_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "CA1Hoy753gVYtifiO7STHH") == 0) {
          extern mxArray *sf_c18_Model_02_updateBuildInfo_args_info(void);
          plhs[0] = sf_c18_Model_02_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void Model_02_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _Model_02MachineNumber_ = sf_debug_initialize_machine(debugInstance,"Model_02",
    "sfun",0,16,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_Model_02MachineNumber_,0,
    0);
  sf_debug_set_machine_data_thresholds(debugInstance,_Model_02MachineNumber_,0);
}

void Model_02_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_Model_02_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("Model_02",
      "Model_02");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_Model_02_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}

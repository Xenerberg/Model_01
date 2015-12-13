/* Include files */

#include "Model_01_sfun.h"
#include "Model_01_sfun_debug_macros.h"
#include "c1_Model_01.h"
#include "c2_Model_01.h"
#include "c3_Model_01.h"
#include "c4_Model_01.h"
#include "c5_Model_01.h"
#include "c6_Model_01.h"
#include "c7_Model_01.h"
#include "c8_Model_01.h"
#include "c9_Model_01.h"
#include "c10_Model_01.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _Model_01MachineNumber_;

/* Function Declarations */

/* Function Definitions */
void Model_01_initializer(void)
{
}

void Model_01_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_Model_01_method_dispatcher(SimStruct *simstructPtr, unsigned int
  chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==1) {
    c1_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==2) {
    c2_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==3) {
    c3_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==4) {
    c4_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==5) {
    c5_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==6) {
    c6_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==7) {
    c7_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==8) {
    c8_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==9) {
    c9_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==10) {
    c10_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_Model_01_process_check_sum_call( int nlhs, mxArray * plhs[], int
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
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2395467678U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2655813948U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1733208674U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3001864846U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1378275903U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2746056357U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1013953898U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(268299227U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 1:
        {
          extern void sf_c1_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c1_Model_01_get_check_sum(plhs);
          break;
        }

       case 2:
        {
          extern void sf_c2_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c2_Model_01_get_check_sum(plhs);
          break;
        }

       case 3:
        {
          extern void sf_c3_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c3_Model_01_get_check_sum(plhs);
          break;
        }

       case 4:
        {
          extern void sf_c4_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c4_Model_01_get_check_sum(plhs);
          break;
        }

       case 5:
        {
          extern void sf_c5_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c5_Model_01_get_check_sum(plhs);
          break;
        }

       case 6:
        {
          extern void sf_c6_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c6_Model_01_get_check_sum(plhs);
          break;
        }

       case 7:
        {
          extern void sf_c7_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c7_Model_01_get_check_sum(plhs);
          break;
        }

       case 8:
        {
          extern void sf_c8_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c8_Model_01_get_check_sum(plhs);
          break;
        }

       case 9:
        {
          extern void sf_c9_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c9_Model_01_get_check_sum(plhs);
          break;
        }

       case 10:
        {
          extern void sf_c10_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c10_Model_01_get_check_sum(plhs);
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
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(39944841U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3539119247U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3907438040U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(345630706U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Model_01_autoinheritance_info( int nlhs, mxArray * plhs[], int
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
        if (strcmp(aiChksum, "mh8d8Qhhhm8NVF2I6X4oeH") == 0) {
          extern mxArray *sf_c1_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c1_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 2:
      {
        if (strcmp(aiChksum, "B1tON1hdGmzSnYOY95SP4B") == 0) {
          extern mxArray *sf_c2_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c2_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 3:
      {
        if (strcmp(aiChksum, "EGwcflfkiYxPecBCz3KLTH") == 0) {
          extern mxArray *sf_c3_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c3_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 4:
      {
        if (strcmp(aiChksum, "BSLVtnxiujwcYccoMPm0xD") == 0) {
          extern mxArray *sf_c4_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c4_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 5:
      {
        if (strcmp(aiChksum, "1SXmiNdsrSH20fPu4M4XFB") == 0) {
          extern mxArray *sf_c5_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c5_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 6:
      {
        if (strcmp(aiChksum, "v74Dva2dvAgckH566yBRCB") == 0) {
          extern mxArray *sf_c6_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c6_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 7:
      {
        if (strcmp(aiChksum, "u2fKBbGJm8l67sYcBYBc4G") == 0) {
          extern mxArray *sf_c7_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c7_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 8:
      {
        if (strcmp(aiChksum, "StOLOaIbyTPtcQv8FGc3PC") == 0) {
          extern mxArray *sf_c8_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c8_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 9:
      {
        if (strcmp(aiChksum, "Sq0CQBHL4wQCOelA2JdeiH") == 0) {
          extern mxArray *sf_c9_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c9_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 10:
      {
        if (strcmp(aiChksum, "JNHPi3lXTvdiOQ1pPckTcH") == 0) {
          extern mxArray *sf_c10_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c10_Model_01_get_autoinheritance_info();
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

unsigned int sf_Model_01_get_eml_resolved_functions_info( int nlhs, mxArray *
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
        extern const mxArray *sf_c1_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c1_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 2:
      {
        extern const mxArray *sf_c2_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 3:
      {
        extern const mxArray *sf_c3_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c3_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 4:
      {
        extern const mxArray *sf_c4_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c4_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 5:
      {
        extern const mxArray *sf_c5_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c5_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 6:
      {
        extern const mxArray *sf_c6_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c6_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 7:
      {
        extern const mxArray *sf_c7_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c7_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 8:
      {
        extern const mxArray *sf_c8_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c8_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 9:
      {
        extern const mxArray *sf_c9_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c9_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 10:
      {
        extern const mxArray *sf_c10_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c10_Model_01_get_eml_resolved_functions_info();
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

unsigned int sf_Model_01_third_party_uses_info( int nlhs, mxArray * plhs[], int
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
        if (strcmp(tpChksum, "WXBmQUuIQcrIgHsskbGIjC") == 0) {
          extern mxArray *sf_c1_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c1_Model_01_third_party_uses_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "bSvlT79syXJgZtTFVgmsyB") == 0) {
          extern mxArray *sf_c2_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c2_Model_01_third_party_uses_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "zpdI4lpeEhcXCoU3wpKKiF") == 0) {
          extern mxArray *sf_c3_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c3_Model_01_third_party_uses_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "tGF1mbJPGUsxezq4V9CaKB") == 0) {
          extern mxArray *sf_c4_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c4_Model_01_third_party_uses_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "iX9Zl6FlNLqtZng6BCwrjF") == 0) {
          extern mxArray *sf_c5_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c5_Model_01_third_party_uses_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "Dl2Q4AOdZkowJfWOlWgDe") == 0) {
          extern mxArray *sf_c6_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c6_Model_01_third_party_uses_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "BYHJKIDcMcOyp0LDH6Zm3D") == 0) {
          extern mxArray *sf_c7_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c7_Model_01_third_party_uses_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "aWEZbcySXBVyXqdTVPBGpG") == 0) {
          extern mxArray *sf_c8_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c8_Model_01_third_party_uses_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "23i2VtcMhW0jI67I50N18C") == 0) {
          extern mxArray *sf_c9_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c9_Model_01_third_party_uses_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "5eb6GG5yWF4WdG34nQsttE") == 0) {
          extern mxArray *sf_c10_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c10_Model_01_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_Model_01_updateBuildInfo_args_info( int nlhs, mxArray * plhs[],
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
        if (strcmp(tpChksum, "WXBmQUuIQcrIgHsskbGIjC") == 0) {
          extern mxArray *sf_c1_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c1_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "bSvlT79syXJgZtTFVgmsyB") == 0) {
          extern mxArray *sf_c2_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "zpdI4lpeEhcXCoU3wpKKiF") == 0) {
          extern mxArray *sf_c3_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c3_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "tGF1mbJPGUsxezq4V9CaKB") == 0) {
          extern mxArray *sf_c4_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c4_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "iX9Zl6FlNLqtZng6BCwrjF") == 0) {
          extern mxArray *sf_c5_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c5_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "Dl2Q4AOdZkowJfWOlWgDe") == 0) {
          extern mxArray *sf_c6_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c6_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "BYHJKIDcMcOyp0LDH6Zm3D") == 0) {
          extern mxArray *sf_c7_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c7_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "aWEZbcySXBVyXqdTVPBGpG") == 0) {
          extern mxArray *sf_c8_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c8_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "23i2VtcMhW0jI67I50N18C") == 0) {
          extern mxArray *sf_c9_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c9_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "5eb6GG5yWF4WdG34nQsttE") == 0) {
          extern mxArray *sf_c10_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c10_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void Model_01_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _Model_01MachineNumber_ = sf_debug_initialize_machine(debugInstance,"Model_01",
    "sfun",0,10,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_Model_01MachineNumber_,0,
    0);
  sf_debug_set_machine_data_thresholds(debugInstance,_Model_01MachineNumber_,0);
}

void Model_01_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_Model_01_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("Model_01",
      "Model_01");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_Model_01_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}

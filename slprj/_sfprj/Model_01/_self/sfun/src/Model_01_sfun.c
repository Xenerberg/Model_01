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
#include "c11_Model_01.h"
#include "c12_Model_01.h"
#include "c13_Model_01.h"
#include "c14_Model_01.h"
#include "c15_Model_01.h"
#include "c16_Model_01.h"
#include "c17_Model_01.h"
#include "c18_Model_01.h"
#include "c19_Model_01.h"
#include "c20_Model_01.h"
#include "c21_Model_01.h"
#include "c22_Model_01.h"
#include "c23_Model_01.h"
#include "c24_Model_01.h"
#include "c25_Model_01.h"

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

  if (chartFileNumber==11) {
    c11_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==12) {
    c12_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==13) {
    c13_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==14) {
    c14_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==15) {
    c15_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==16) {
    c16_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==17) {
    c17_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==18) {
    c18_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==19) {
    c19_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==20) {
    c20_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==21) {
    c21_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==22) {
    c22_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==23) {
    c23_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==24) {
    c24_Model_01_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==25) {
    c25_Model_01_method_dispatcher(simstructPtr, method, data);
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
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1467923453U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2112315684U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(892916413U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1268429609U);
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

       case 11:
        {
          extern void sf_c11_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c11_Model_01_get_check_sum(plhs);
          break;
        }

       case 12:
        {
          extern void sf_c12_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c12_Model_01_get_check_sum(plhs);
          break;
        }

       case 13:
        {
          extern void sf_c13_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c13_Model_01_get_check_sum(plhs);
          break;
        }

       case 14:
        {
          extern void sf_c14_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c14_Model_01_get_check_sum(plhs);
          break;
        }

       case 15:
        {
          extern void sf_c15_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c15_Model_01_get_check_sum(plhs);
          break;
        }

       case 16:
        {
          extern void sf_c16_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c16_Model_01_get_check_sum(plhs);
          break;
        }

       case 17:
        {
          extern void sf_c17_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c17_Model_01_get_check_sum(plhs);
          break;
        }

       case 18:
        {
          extern void sf_c18_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c18_Model_01_get_check_sum(plhs);
          break;
        }

       case 19:
        {
          extern void sf_c19_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c19_Model_01_get_check_sum(plhs);
          break;
        }

       case 20:
        {
          extern void sf_c20_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c20_Model_01_get_check_sum(plhs);
          break;
        }

       case 21:
        {
          extern void sf_c21_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c21_Model_01_get_check_sum(plhs);
          break;
        }

       case 22:
        {
          extern void sf_c22_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c22_Model_01_get_check_sum(plhs);
          break;
        }

       case 23:
        {
          extern void sf_c23_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c23_Model_01_get_check_sum(plhs);
          break;
        }

       case 24:
        {
          extern void sf_c24_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c24_Model_01_get_check_sum(plhs);
          break;
        }

       case 25:
        {
          extern void sf_c25_Model_01_get_check_sum(mxArray *plhs[]);
          sf_c25_Model_01_get_check_sum(plhs);
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
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3102475941U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1259861897U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(268821384U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3536486563U);
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

     case 11:
      {
        if (strcmp(aiChksum, "IS5kaJJMSAHNXsjGcxQjbB") == 0) {
          extern mxArray *sf_c11_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c11_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 12:
      {
        if (strcmp(aiChksum, "Z5zcJhuDQMmn7USpJLM42F") == 0) {
          extern mxArray *sf_c12_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c12_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 13:
      {
        if (strcmp(aiChksum, "HxYR2xACudhglmitq1hOGC") == 0) {
          extern mxArray *sf_c13_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c13_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 14:
      {
        if (strcmp(aiChksum, "i1F2hgPs1THBF49GEkUOIF") == 0) {
          extern mxArray *sf_c14_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c14_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 15:
      {
        if (strcmp(aiChksum, "yfela4PNL1hjBN6seFvK2B") == 0) {
          extern mxArray *sf_c15_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c15_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 16:
      {
        if (strcmp(aiChksum, "Y8MIbm9zgLrDSY06AU6nOB") == 0) {
          extern mxArray *sf_c16_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c16_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 17:
      {
        if (strcmp(aiChksum, "UeYxRy69v8HjCdEadZZsSD") == 0) {
          extern mxArray *sf_c17_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c17_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 18:
      {
        if (strcmp(aiChksum, "06BKYiKMrCeapG4pb5meFD") == 0) {
          extern mxArray *sf_c18_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c18_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 19:
      {
        if (strcmp(aiChksum, "QV4wCykUN7fTOLJMPe6HZB") == 0) {
          extern mxArray *sf_c19_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c19_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 20:
      {
        if (strcmp(aiChksum, "5Gr85UYM9WW8xen0JoONaG") == 0) {
          extern mxArray *sf_c20_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c20_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 21:
      {
        if (strcmp(aiChksum, "wGlXd845UI3k1QllG41hKG") == 0) {
          extern mxArray *sf_c21_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c21_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 22:
      {
        if (strcmp(aiChksum, "3znuVi0YRCKeQiaOhIYOfB") == 0) {
          extern mxArray *sf_c22_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c22_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 23:
      {
        if (strcmp(aiChksum, "5giLtXLjUUdWL7aCBaeNOF") == 0) {
          extern mxArray *sf_c23_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c23_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 24:
      {
        if (strcmp(aiChksum, "WqHAg9XJ7ZlFUzRTrkt8ND") == 0) {
          extern mxArray *sf_c24_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c24_Model_01_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 25:
      {
        if (strcmp(aiChksum, "Fl47OdJDSWCKRT01czohIG") == 0) {
          extern mxArray *sf_c25_Model_01_get_autoinheritance_info(void);
          plhs[0] = sf_c25_Model_01_get_autoinheritance_info();
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

     case 11:
      {
        extern const mxArray *sf_c11_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c11_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 12:
      {
        extern const mxArray *sf_c12_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c12_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 13:
      {
        extern const mxArray *sf_c13_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c13_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 14:
      {
        extern const mxArray *sf_c14_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c14_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 15:
      {
        extern const mxArray *sf_c15_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c15_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 16:
      {
        extern const mxArray *sf_c16_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c16_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 17:
      {
        extern const mxArray *sf_c17_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c17_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 18:
      {
        extern const mxArray *sf_c18_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c18_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 19:
      {
        extern const mxArray *sf_c19_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c19_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 20:
      {
        extern const mxArray *sf_c20_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c20_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 21:
      {
        extern const mxArray *sf_c21_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c21_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 22:
      {
        extern const mxArray *sf_c22_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c22_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 23:
      {
        extern const mxArray *sf_c23_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c23_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 24:
      {
        extern const mxArray *sf_c24_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c24_Model_01_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 25:
      {
        extern const mxArray *sf_c25_Model_01_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c25_Model_01_get_eml_resolved_functions_info();
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

     case 11:
      {
        if (strcmp(tpChksum, "QwETwijcyK4H4WrgvmSxxG") == 0) {
          extern mxArray *sf_c11_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c11_Model_01_third_party_uses_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "DE5yHmseDeJz1ulxRRysRE") == 0) {
          extern mxArray *sf_c12_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c12_Model_01_third_party_uses_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "H7mISRpE6SQo3UP9VXgsDE") == 0) {
          extern mxArray *sf_c13_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c13_Model_01_third_party_uses_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "85DGOMCP8MjTH9hhiMHlSE") == 0) {
          extern mxArray *sf_c14_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c14_Model_01_third_party_uses_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "7DZY7rufJy9vvlu0fijQQH") == 0) {
          extern mxArray *sf_c15_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c15_Model_01_third_party_uses_info();
          break;
        }
      }

     case 16:
      {
        if (strcmp(tpChksum, "IuSSsJZIfFJ6iHjFwXSqXF") == 0) {
          extern mxArray *sf_c16_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c16_Model_01_third_party_uses_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "KQBb6GGKaZ9GEqWQrEZmKH") == 0) {
          extern mxArray *sf_c17_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c17_Model_01_third_party_uses_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "9ED3OSWUWBvLguXkiXLqRC") == 0) {
          extern mxArray *sf_c18_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c18_Model_01_third_party_uses_info();
          break;
        }
      }

     case 19:
      {
        if (strcmp(tpChksum, "Tl9sDoIkjP3E007JlGDTCE") == 0) {
          extern mxArray *sf_c19_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c19_Model_01_third_party_uses_info();
          break;
        }
      }

     case 20:
      {
        if (strcmp(tpChksum, "wphbv33LDzi0MTnqKdWGDE") == 0) {
          extern mxArray *sf_c20_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c20_Model_01_third_party_uses_info();
          break;
        }
      }

     case 21:
      {
        if (strcmp(tpChksum, "vg45KOUeguIO0CLb0w70oB") == 0) {
          extern mxArray *sf_c21_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c21_Model_01_third_party_uses_info();
          break;
        }
      }

     case 22:
      {
        if (strcmp(tpChksum, "MnGZIoF4qkyXhgrzEsOubG") == 0) {
          extern mxArray *sf_c22_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c22_Model_01_third_party_uses_info();
          break;
        }
      }

     case 23:
      {
        if (strcmp(tpChksum, "J3FURByDxB3py7ir0RMGL") == 0) {
          extern mxArray *sf_c23_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c23_Model_01_third_party_uses_info();
          break;
        }
      }

     case 24:
      {
        if (strcmp(tpChksum, "N7Jw8eSBiPMnv467LAHSXF") == 0) {
          extern mxArray *sf_c24_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c24_Model_01_third_party_uses_info();
          break;
        }
      }

     case 25:
      {
        if (strcmp(tpChksum, "dwoFzTByORxjD9xyp4fZcE") == 0) {
          extern mxArray *sf_c25_Model_01_third_party_uses_info(void);
          plhs[0] = sf_c25_Model_01_third_party_uses_info();
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

     case 11:
      {
        if (strcmp(tpChksum, "QwETwijcyK4H4WrgvmSxxG") == 0) {
          extern mxArray *sf_c11_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c11_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "DE5yHmseDeJz1ulxRRysRE") == 0) {
          extern mxArray *sf_c12_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c12_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "H7mISRpE6SQo3UP9VXgsDE") == 0) {
          extern mxArray *sf_c13_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c13_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "85DGOMCP8MjTH9hhiMHlSE") == 0) {
          extern mxArray *sf_c14_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c14_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "7DZY7rufJy9vvlu0fijQQH") == 0) {
          extern mxArray *sf_c15_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c15_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 16:
      {
        if (strcmp(tpChksum, "IuSSsJZIfFJ6iHjFwXSqXF") == 0) {
          extern mxArray *sf_c16_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c16_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "KQBb6GGKaZ9GEqWQrEZmKH") == 0) {
          extern mxArray *sf_c17_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c17_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "9ED3OSWUWBvLguXkiXLqRC") == 0) {
          extern mxArray *sf_c18_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c18_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 19:
      {
        if (strcmp(tpChksum, "Tl9sDoIkjP3E007JlGDTCE") == 0) {
          extern mxArray *sf_c19_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c19_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 20:
      {
        if (strcmp(tpChksum, "wphbv33LDzi0MTnqKdWGDE") == 0) {
          extern mxArray *sf_c20_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c20_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 21:
      {
        if (strcmp(tpChksum, "vg45KOUeguIO0CLb0w70oB") == 0) {
          extern mxArray *sf_c21_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c21_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 22:
      {
        if (strcmp(tpChksum, "MnGZIoF4qkyXhgrzEsOubG") == 0) {
          extern mxArray *sf_c22_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c22_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 23:
      {
        if (strcmp(tpChksum, "J3FURByDxB3py7ir0RMGL") == 0) {
          extern mxArray *sf_c23_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c23_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 24:
      {
        if (strcmp(tpChksum, "N7Jw8eSBiPMnv467LAHSXF") == 0) {
          extern mxArray *sf_c24_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c24_Model_01_updateBuildInfo_args_info();
          break;
        }
      }

     case 25:
      {
        if (strcmp(tpChksum, "dwoFzTByORxjD9xyp4fZcE") == 0) {
          extern mxArray *sf_c25_Model_01_updateBuildInfo_args_info(void);
          plhs[0] = sf_c25_Model_01_updateBuildInfo_args_info();
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
    "sfun",0,25,0,0,0);
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

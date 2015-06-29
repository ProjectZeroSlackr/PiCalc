#include "pi.h"
#include "bigint.h"
#include "bigmul.h"
#include "agm.h"
#include <time.h>

/*
** This file is a 'debugging' / 'tuning' type of file that does a bunch
** of timings for the various functions.  You can use it too, just to
** make sure that everything works okay on your system.
**
** These routines are activated by passing a NEGATIVE number of digits
** to the pi program.  As long as this routine is compiled to do the
** testings, and not to abort with an error, of course.
*/

static volatile clock_t TestTime,ElapTime;

static void
StartTimings(char *Str)
{
ClearTimings();
fprintf(stdout,"%s: ",Str);fflush(stdout);
/*fprintf(stderr,"%s: ",Str);fflush(stderr);*/

ElapTime=clock();do {TestTime=clock();} while (ElapTime==TestTime);
}

static void
EndTimings(void)
{
ElapTime=clock()-TestTime;
fprintf(stdout,"Total time=%.0f ticks\n",(double)ElapTime);
/*fprintf(stderr,"Total time=%.0f ticks\n",(double)ElapTime);*/
/*
DumpTimings(ElapTime/CLOCKS_PER_SEC);
*/
}


void
TestMath(size_t Len)
{
  BigInt AGM_A, AGM_B, AGM_C, AGMWork;

  ClearTimings();
  fprintf(stderr,"\nPerform some timings and self checks on numbers that\n");
  fprintf(stderr,"are %lu digits long.  The self checks, multiplies,\n",
          (ULINT)Len*RawIntDigits);
  fprintf(stderr,"square root and reciprocation may take a while.\n\n");

  AGM_A = CreateBigInt(Len);
  AGM_B = CreateBigInt(Len);
  AGM_C = CreateBigInt(Len);
  AGMWork = CreateBigInt(Len*2);

  Num1IsCached = Num2IsCached = 0;

  fprintf(stderr, "Initialization\n");
  SetNum(AGM_A,Len,BI_One,0);
  SetNum(AGM_C,Len,BI_One,0);
  ClearBigInt(AGMWork,Len);

  fprintf(stderr,"There are %d clock ticks per second.\n",(int)CLOCKS_PER_SEC);

#if 0
  StartTimings("SetNum");
  SetNum(AGM_C,Len,BI_One,0);
  EndTimings();

  StartTimings("ClearBigInt");
  ClearBigInt(AGMWork,Len);
  EndTimings();

  StartTimings("Add one var");
  Add(AGMWork,AGM_A,AGM_A,Len);
  EndTimings();

  StartTimings("Add two vars");
  Add(AGMWork,AGM_A,AGM_B,Len);
  EndTimings();

  StartTimings("Copy");
  Copy(AGMWork,AGM_A,Len);
  EndTimings();

  {size_t x,y;
  StartTimings("CountZeros");
  CountZeros(AGM_A,&x,&y,Len);
  EndTimings();
  }

  StartTimings("DivBy 2");
  DivBy(AGMWork,AGMWork,2,Len);
  EndTimings();

  StartTimings("DivBy 5");
  DivBy(AGMWork,AGMWork,5,Len);
  EndTimings();

  StartTimings("MulBy 2");
  MulBy(AGMWork,AGMWork,2,Len);
  EndTimings();

  StartTimings("MulBy 10");
  MulBy(AGMWork,AGMWork,10,Len);
  EndTimings();

  StartTimings("MulBy float");
  MulByFloat(AGMWork,17.0,Len);
  EndTimings();

  StartTimings("Negate");
  Negate(AGM_A,Len);
  EndTimings();

  SetNum(AGM_A,Len,One,0);
  StartTimings("RevSubInt 5-1");
  RevSubInt(5,AGM_A,Len);
  EndTimings();

  SetNum(AGM_A,Len,5*BI_One,0);
  StartTimings("RevSubInt 1-5");
  RevSubInt(1,AGM_A,Len);
  EndTimings();

  StartTimings("Sub");
  Sub(AGM_A,AGM_A,AGM_C,Len);
  EndTimings();

  StartTimings("SpecialAGMFunc1");
  SpecialAGMFunc1(AGM_A,AGM_B,AGM_C,Len);
  EndTimings();

  StartTimings("HalfDiff");
  HalfDiff(AGM_A,AGM_B,AGM_C,Len);
  EndTimings();
#endif

  SetNum(AGM_A,Len,BI_One,0);
  SetNum(AGM_C,Len,BI_One,0);
  ClearBigInt(AGMWork,Len);

#if 1
/* Set AGM_A to 99999999...99999 */
  ClearBigInt(AGM_A,Len);
  ClearBigInt(AGM_C,Len);SetBigIntDigit(AGM_C,Len-1,1);
  Sub(AGM_A,AGM_A,AGM_C,Len);
  StartTimings("Testing square FFT multiply");
/*  FullMul(AGMWork,AGM_A,AGM_A,Len);*/
  FFTMul(AGMWork,AGM_A,AGM_A,Len,Len*2,0);
  EndTimings();
/* Check it against 999...98000...01 */
  ClearBigInt(AGM_A,Len);
  ClearBigInt(AGM_C,Len);SetBigIntDigit(AGM_C,Len-1,2);
  Sub(AGM_A,AGM_A,AGM_C,Len);
  Sub(AGM_B,AGM_A,AGMWork,Len);
  if (!IsZero(AGM_B,Len)) FatalError("Sqr mul failed first check.\n");
  if ( (GetBigIntDigit(AGMWork,Len*2-1) != 1) ||
       (!IsZero(AGMWork+Len,Len-1)))
    FatalError("Sqr mul failed second check.\n");
/*  fprintf(stderr,"Sqr passed.\n");*/
#endif

#if 0
/* Set AGM_A and AGM_C to 99999999...99999 */
  ClearBigInt(AGM_A,Len);
  ClearBigInt(AGM_C,Len);SetBigIntDigit(AGM_C,Len-1,1);
  Sub(AGM_A,AGM_A,AGM_C,Len);
  Copy(AGM_C,AGM_A,Len);
  StartTimings("Testing non-square FFT multiply");
/*  FullMul(AGMWork,AGM_A,AGM_C,Len);*/
  FFTMul(AGMWork,AGM_A,AGM_C,Len,Len*2,0);
  EndTimings();
/* Check it against 999...98000...01 */
  ClearBigInt(AGM_A,Len);
  ClearBigInt(AGM_C,Len);SetBigIntDigit(AGM_C,Len-1,2);
  Sub(AGM_A,AGM_A,AGM_C,Len);
  Sub(AGM_B,AGM_A,AGMWork,Len);
  if (!IsZero(AGM_B,Len)) FatalError("Non-Sqr mul failed first check.\n");
  if ( (GetBigIntDigit(AGMWork,Len*2-1) != 1) ||
       (!IsZero(AGMWork+Len,Len-1)))
    FatalError("Non-Sqr mul failed second check.\n");
/*  fprintf(stderr,"Non-Sqr passed.\n");*/
#endif

#if 0
  StartTimings("Testing square root 0.5");
  Sqrt05(AGM_B,Len);
  EndTimings();

  StartTimings("Testing square root of sqrt(0.5)");
  AGMSqrt(AGM_A,AGM_B,Len,0);
  EndTimings();

  SetNum(AGM_C,Len,9138931,62088927);
  SetNum(AGM_B,Len,BI_Two,0);
  StartTimings("Divide");
  AGMDivide(AGM_A,AGM_B,AGM_C,Len,AGMWork);
  EndTimings();
#endif

  fprintf(stderr,"Testing done.\n");
  DeInitBigIntPkg();

  ExitPrg(EXIT_SUCCESS);
}



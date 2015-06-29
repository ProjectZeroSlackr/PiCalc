#include "pi.h"
#include "agm.h"
#include "bigint.h"
#include "cache.h"
#include "bigmul.h"


/*
** d = a/b by computing the reciprocal of b and then multiplying
** that by a.
*/
void
AGMDivide(BigInt R, BigInt Num1, BigInt Num2, size_t Len, BigInt Work)
{ int Sign;
  size_t SubLen=2;
  size_t Redo=REDO_LEN;

  if      (NumIs(Num2, 9138931,62088927)) SetNum(R,Len,10942198, 7613238);
  else if (NumIs(Num2,12300397,35639667)) SetNum(R,Len, 8129818,66378456);
  else
    {
     DumpBigInt("Unknown Divisor: ",Num2,4);
     ExitPrg(EXIT_FAILURE);
    }

  if (SubLen >= Redo) Redo=0;
  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
  while (SubLen < Len/2)
    {
      SubLen *= 2;if (SubLen > Len) SubLen = Len;
      if (!Cfg.Macintosh)
        fprintf(stderr,"Recip: %4s",Num2Str(SubLen*RawIntDigits));
      FlushFFTCache(0);

/* Perform safety check */
      if ( (strcmp(GetCheckStr(R),"1094219807613238")!=0) &&
           (strcmp(GetCheckStr(R),"0812981866378456")!=0) )
        fprintf(stderr,"** WARNING **\a\nAGMDivide may be failing.\n%s\n",GetCheckStr(R));

      SaveNum1FFT = 10;
      if (SubLen==Len/2) SaveNum2FFT=11;
      ClearBigInt(R+SubLen/2,SubLen/2);
      if (!Cfg.Macintosh) fputc('|',stderr);
      N1R0Mul(DSWork, R, Num2, Work,SubLen);

      Sign = RevSubInt(BI_One,DSWork,SubLen);

      Num1IsCached=10;
      if (!Cfg.Macintosh) fputc('.',stderr);
      HalfMul(DSWork,R,DSWork,SubLen);
      if (Sign) Sub(R,R,DSWork,SubLen);
      else      Add(R,R,DSWork,SubLen);

      if (Redo == ULONG_MAX) Redo=0;
      if (SubLen == Redo) {SubLen/=2;Redo=ULONG_MAX;}
      if (!Cfg.Macintosh) BackSpace(13);
    }

  if (!Cfg.Macintosh) fprintf(stderr,"Recip: %4s",Num2Str(Len*RawIntDigits));
  FlushFFTCache(11);
  ClearBigInt(R+Len/2,Len/2);
  SaveNum1FFT = 10;
  if (!Cfg.Macintosh) fputc('.',stderr);
  HalfMul(OldRoot,R,Num1,Len);

  ClearBigInt(OldRoot+Len/2,Len/2);
  Num2IsCached=11;
  if (!Cfg.Macintosh) fputc('|',stderr);
  N1R0Mul(DSWork,OldRoot,Num2,Work,Len);

  Sign=Sub(DSWork,Num1,DSWork,Len);
  if (Sign) Negate(DSWork,Len);

  Num1IsCached=10;
  if (!Cfg.Macintosh) fputc('.',stderr);
  HalfMul(DSWork,R,DSWork,Len);
  if (Sign) Sub(R,OldRoot,DSWork,Len);
  else      Add(R,OldRoot,DSWork,Len);

  if (!Cfg.Macintosh) BackSpace(14);
  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
}

/*
** d = a/b by computing the reciprocal of b and then multiplying
** that by a.
**
** r = r*(2-br)
** d = a * r
*/
void
ClassicDivide(BigInt D, BigInt Num, BigInt Denom, size_t Len)
{ size_t SubLen;
  size_t Redo=REDO_LEN;

  {INT32 PreSet;
   double Prep;
   Prep=(1.0*BI_One)/GetBigIntDigit(Denom,0);
   PreSet=Prep*BI_One;
   SubLen=1;
   SetNum(D,1,PreSet,0);
   ClearBigInt(D+SubLen,Len-SubLen);
  }

  while (SubLen < Len)
    {
      SubLen *= 2;if (SubLen > Len) SubLen = Len;
      if (!Cfg.Macintosh) fprintf(stderr,"Div: %4s",Num2Str(SubLen*RawIntDigits));

      FullMul(DSWork, D, Denom, SubLen);
      RevSubInt(BI_Two,DSWork,SubLen);
      FullMul(D, D, DSWork, SubLen);
      if (SubLen == Redo) {SubLen/=2;Redo=0;}
      if (!Cfg.Macintosh) BackSpace(9);
    }

  if (!Cfg.Macintosh) fprintf(stderr,"Div: Mul");
  FullMul(D, Num, D, Len);
  if (!Cfg.Macintosh) BackSpace(8);
}


/*
** Computes the square root of 'n' by computing its reciprocal square
** root using the Newton's method, and then multiplying by the original
** number.  This is faster than doing the normal Newton 'divide & average'
** method.
*/
/* static */ void
AGMSqrt(BigInt Root, BigInt Num, size_t Len, size_t SubLen)
{
  int Sign;
  size_t Redo=REDO_LEN;

  if (SubLen <= 0) SubLen = 2;
  if (SubLen > Len) SubLen = Len;

  if      (NumIs(Num, 7071067,81186547)) SetNum(OldRoot,2,11892071,15002721);
  else if (NumIs(Num, 7177499,86377537)) SetNum(OldRoot,2,11803570,64195417);
  else if (NumIs(Num, 7177700,10976296)) SetNum(OldRoot,2,11803405,99073515);
  else if (NumIs(Num, 7177700,11046129)) SetNum(OldRoot,2,11803405,99016096);
/* The second AGM */
  else if (NumIs(Num, 9659258,26289068)) SetNum(OldRoot,2,10174852,23681446);
  else if (NumIs(Num, 9660709,46551927)) SetNum(OldRoot,2,10174087,99031044);
  else if (NumIs(Num, 9660709,49277277)) SetNum(OldRoot,2,10174087,97595956);

  else
    {
     DumpBigInt("Unknown Sqrt: ",Num,4);
     ExitPrg(EXIT_FAILURE);
    }
  ClearBigInt(OldRoot+SubLen,Len-SubLen);

  if (SubLen >= Redo) Redo=0;
  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
  while (SubLen < Len/2)
    {
      SubLen *= 2;if (SubLen > Len) SubLen = Len;
      if (!Cfg.Macintosh)
        fprintf(stderr,"Sqrt: %4s",Num2Str(SubLen*RawIntDigits));
      FlushFFTCache(0);

/* Perform safety check */
      {char *Str=GetCheckStr(OldRoot);
       if ( (strcmp(Str,"1189207115002721")!=0) &&
            (strcmp(Str,"1180357064195417")!=0) &&
            (strcmp(Str,"1180340599073515")!=0) &&
            (strcmp(Str,"1180340599016096")!=0) &&

            (strcmp(Str,"1017485223681446")!=0) &&
            (strcmp(Str,"1017408799031044")!=0) &&
            (strcmp(Str,"1017408797595956")!=0)
          )
         fprintf(stderr,"** WARNING **\a\nAGMSqrt may be failing.\n%s\n",Str);
      }

      ClearBigInt(OldRoot+SubLen/2,SubLen/2);
      if (!Cfg.Macintosh) fputc('.',stderr);
      SaveNum1FFT = 20;
      HalfMul(DSWork, OldRoot, OldRoot, SubLen);
      if (SubLen == Len/2) SaveNum1FFT = 21;
      if (!Cfg.Macintosh) fputc('|',stderr);
      FullMul(DSWork, Num, DSWork, SubLen);

      Sign = RevSubInt(BI_One,DSWork,SubLen);

      Num1IsCached=20;
      if (!Cfg.Macintosh) fputc('.',stderr);
      HalfMul2(DSWork,OldRoot,DSWork,SubLen);

      if (Sign) Sub(OldRoot,OldRoot,DSWork,SubLen);
      else      Add(OldRoot,OldRoot,DSWork,SubLen);

      if (Redo == ULONG_MAX) Redo = 0;
      if (SubLen == Redo) {SubLen/=2;Redo=ULONG_MAX;}
      if (!Cfg.Macintosh) BackSpace(13);
    }

  if (!Cfg.Macintosh) fprintf(stderr,"Sqrt: %4s",Num2Str(Len*RawIntDigits));
  FlushFFTCache(21);
  if (!Cfg.Macintosh) fputc('.',stderr);
  SaveNum2FFT=20;
  Num1IsCached=21;
  HalfMul(Root, Num, OldRoot, Len);
  ClearBigInt(Root+Len/2,Len/2);
  if (!Cfg.Macintosh) fputc('.',stderr);
  HalfMul(DSWork, Root, Root, Len);

  Sign=Sub(DSWork, Num, DSWork, Len);
  if (Sign) Negate(DSWork,Len);

  Num1IsCached = 20;
  if (!Cfg.Macintosh) fputc('.',stderr);
  HalfMul2(DSWork,OldRoot,DSWork,Len);

  if (Sign) Sub(Root, Root, DSWork, Len);
  else      Add(Root, Root, DSWork, Len);

  if (!Cfg.Macintosh) BackSpace(13);
  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
}

/*
** Computes the square root of 'n' by computing the reciprocal and then
** multiplying that by the original number.
** r = r*(3-n*r^2)/2
*/
void
ClassicSqrt(BigInt Root, BigInt Num, size_t Len)
{
  size_t SubLen=1;
  size_t Redo=REDO_LEN;

  {INT32 PreSet;
   double Prep;
   Prep=GetBigIntDigit(Num,0)/((double)BI_One);
   Prep=sqrt(1.0/Prep);
   PreSet=Prep*BI_One;
   SubLen=1;
   SetNum(Root,1,PreSet,0);
   ClearBigInt(Root+SubLen,Len-SubLen);
  }

  while (SubLen < Len)
    {
      SubLen *= 2;if (SubLen > Len) SubLen = Len;
      if (!Cfg.Macintosh)
        fprintf(stderr,"Sqrt: %4s",Num2Str(SubLen*RawIntDigits));

      FullMul(DSWork, Root, Root, SubLen);
      FullMul(DSWork, Num, DSWork, SubLen);
      RevSubInt(BI_Three,DSWork,SubLen);
      FullMul(Root,Root,DSWork,SubLen);
      DivBy(Root, Root, 2, SubLen);
      if (SubLen == Redo) {SubLen/=2;Redo=0;}
      if (!Cfg.Macintosh) BackSpace(10);
    }

  if (!Cfg.Macintosh) fprintf(stderr,"Sqrt: Mul");
  FullMul(Root,Num,Root,Len);
  if (!Cfg.Macintosh) BackSpace(9);
}

static void
InvSqrtInt(BigInt Root, INT32 INum,size_t Len)
/*
** Generic routine to compute the inverse of the square root of
** an integer.  If you need the regular form, you can easily
** do a MulBy(x,x,INum,Len); yourself
*/
{
  int Sign;
  size_t Redo=REDO_LEN;
  size_t SubLen=1;
  char Check[128];

  {double Prep;INT32 PreSet;
/* Yuck.  I hate calculating the starting value. */
   Prep=sqrt(1.0/INum);
   Prep=Prep*10000000; PreSet=Prep;
   SetNum(Root,1,PreSet,0);
   ClearBigInt(Root+SubLen,Len-SubLen);
   sprintf(Check,"%08d",PreSet);
  }

  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
  while (SubLen < Len)
    {
      SubLen *= 2;if (SubLen > Len) SubLen = Len;
      if (!Cfg.Macintosh)
        fprintf(stderr,"Sqrt%d: %4s",INum,Num2Str(SubLen*RawIntDigits));
      FlushFFTCache(0);

/* Perform safety check */
      if (strncmp(GetCheckStr(Root),Check,8)!=0)
        fprintf(stderr,"** WARNING **\a\nInvSqrtInt %d may be failing.\n%s vs. %s\n",
                INum,Check,GetCheckStr(Root));

      ClearBigInt(Root+SubLen/2,SubLen/2);
      SaveNum1FFT = 30;
      if (!Cfg.Macintosh) fputc('.',stderr);
      HalfMul(DSWork, Root, Root, SubLen);
      MulBy(DSWork,DSWork,INum,SubLen);

      Sign = RevSubInt(BI_One,DSWork,SubLen);

      Num1IsCached=30;
      if (!Cfg.Macintosh) fputc('.',stderr);
      HalfMul(DSWork,Root,DSWork,SubLen);
      DivBy(DSWork,DSWork,2,SubLen);

      if (Sign) Sub(Root,Root,DSWork,SubLen);
      else      Add(Root,Root,DSWork,SubLen);

      if (Redo == ULONG_MAX) Redo = 0;
      if (SubLen == Redo) {SubLen/=2;Redo=ULONG_MAX;}
      if (!Cfg.Macintosh) BackSpace(13);
    }
 FlushFFTCache(0);
}


/* static */ void
Sqrt05(BigInt Root, size_t Len)
{
InvSqrtInt(Root,2,Len);
}

/* static */ void
Sqrt20(BigInt Root, size_t Len)
{
InvSqrtInt(Root,2,Len);
MulBy(Root,Root,2,Len);
}

/* static */ void
Sqrt30(BigInt Root, size_t Len)
{
InvSqrtInt(Root,3,Len);
MulBy(Root,Root,3,Len);
}

/* static */ void
Sqrt60(BigInt Root, size_t Len)
{
InvSqrtInt(Root,6,Len);
MulBy(Root,Root,6,Len);
}

static void
CorruptVar(BigInt Var,size_t Len)
/*
** Deliberately corrupt a variable so we can test to see if the
** self checking will catch it.
*/
{size_t P;
 INT32 Dig;

if (Len < 1024) return;

srand((unsigned int)time(NULL));
P=rand()%(Len-200);
Dig=GetBigIntDigit(Var,P);
Dig=(INT32)(Dig+(rand()%9998)+1);
Dig=(INT32)(Dig % 99999999);
SetBigIntDigit(Var,P,Dig);
printf("Corrupting position %lu with digit %d\n",(ULINT)P,(INT32)Dig);
}

static void
DoCheck(BigInt Scratch, BigInt Num,int Which, size_t Len)
{size_t z;
if (Len < 1024) return;

if (Sub(Scratch,Scratch,Num,Len)) Negate(Scratch,Len);
z=FindFirstNonZero(Scratch,Len);
if (z < Len-200/RawIntDigits)
  {
   fprintf(stderr,"The AGM self check #%d failed. Exp: %lu Got: %lu\n",
         Which,(ULINT)Len-200/RawIntDigits,(ULINT)z);
/*
     ExitPrg(EXIT_FAILURE);
*/
  }
}


/*
**      The AGM itself
**
** A[0]=1 B[0]=1/sqrt(2) Sum=1
**
** n=1..inf
** A[n] = (A[n-1] + B[n-1])/2
** B[n] = Sqrt(A[n-1]*B[n-1])
** C[n] = (A[n-1]-B[n-1])/2    or:
** C[n]^2 = A[n]^2 - B[n]^2    or:
** C[n]^2 = 4A[n+1]*C[n+1]
** Sum  = Sum - C[n]^2*(2^(n+1))
** PI[n] = 4A[n+1]^2 / Sum
**
** However, it's not implemented that way.  We can save a full
** sized multiplication (with two numbers) by rearranging it,
** and keeping the squares of the A and B.  Also, we can do the
** first pass a little differently since A and A^2 will both be 1.
**
** A[0]=1           A[0]^2=1
** B[0]=1/sqrt(2)   B[0]^2=0.5
**
** First pass:
**
** A[1]     = (A[0] + B[0])/2
** B[1]^2   = A[0]*B[0]  ; Since A[0]==1, B[1]^2=B[0]
** C[1]^2   = ((A[0]^2+B[0]^2)/2-B[1]^2)/2
**
** Remainging passes:
**
** C[n]     = (A[n-1]-B[n-1])/2; C[n] is actually temp usage of C[n]^2
** A[n]     = A[n-1] - C[n]
** C[n]^2   = C[n]*C[n]
** B[n]^2   = (A[n-1]^2+B[n-1]^2-4C[n]^2)/2
**
** Then the rest of the formula is done the same:
**
** A[n]^2   = C[n]^2 + B[n]^2
** B[n]     = sqrt(B[n]^2)
** Sum      = Sum - C[n]^2*(2^(n+1))
**
** It is an unusual arrangment, but they are all derived directly
** from the standard AGM formulas as given by Salamin.
**
*/
static int
ComputeFastAGM(size_t Len, size_t MaxPasses)
{
  int Sign;
  size_t Pass;
  double Pow2;
  clock_t LoopTime,EndTime,StartTime;
  BigInt AGM_A, AGM_B, AGM_Sum, AGM_C2;
  BigInt AGM_A2, AGM_B2; /* Squares */
  int Done; /* boolean */

  if ((Cfg.PiFormulaToUse != 2) && (Cfg.PiFormulaToUse != 3))
    FatalError("FastAGM was somehow called with pi formula %d\n",
               Cfg.PiFormulaToUse);

  if (!Cfg.Macintosh) fprintf(stderr,"Creating vars");
  AGM_A   = CreateBigInt(Len);
  AGM_A2  = CreateBigInt(Len);
  AGM_B   = CreateBigInt(Len);
  AGM_B2  = CreateBigInt(Len);
  AGM_C2  = CreateBigInt(Len);
  AGM_Sum = CreateBigInt(Len);
  if (!Cfg.Macintosh) BackSpace(13);

  Num1IsCached = Num2IsCached = 0;
  StartTime = clock();
  Pow2 = 4.0;
  Pass = 0;
  if (!LoadData(AGM_A,AGM_A2,AGM_B,AGM_B2,AGM_Sum,OldRoot,
    &StartTime,&Pow2,&Pass,Len,Cfg.PiFormulaToUse))
    {
     fprintf(stderr, "Init     : ");
     LoopTime = clock();
     Pow2 = 4.0;
     Pass = 0;

     if (Cfg.PiFormulaToUse==2)
       {
        if (!Cfg.Macintosh) fprintf(stderr,"Setting vars");
        SetNum(AGM_A,Len,  BI_One,0);
        SetNum(AGM_A2,Len, BI_One,0);
        SetNum(AGM_B2,Len, BI_OneHalf,0);
        SetNum(AGM_Sum,Len,BI_One,0);
        ClearBigInt(OldRoot,Len);
        if (!Cfg.Macintosh) BackSpace(12);
        Sqrt05(AGM_B,Len);
        if ((Cfg.AGMSelfCheck == 3) || (Cfg.AGMSelfCheck == 4))
          {
         /* We can use OldRoot and DSWork as scratch. */
           if (!Cfg.Macintosh) fprintf(stderr,"Self Check");
           SpecialSquare(OldRoot,AGM_B,Len,DSWork);
           DoCheck(OldRoot,AGM_B2,0,Len);
           ClearBigInt(OldRoot,Len);
           if (!Cfg.Macintosh) BackSpace(10);
          }
       }
     else
       {
        Sqrt20(AGM_A,Len);
        Sqrt60(AGM_A2,Len);
        if (!Cfg.Macintosh) fprintf(stderr,"Square");
        Add(AGM_B,AGM_A,AGM_A2,Len);
        DivBy(AGM_B,AGM_B,4,Len);
        SpecialSquare(AGM_B2,AGM_B,Len,AGM_A);
        if (!Cfg.Macintosh) BackSpace(6);
        if (!Cfg.Macintosh) fprintf(stderr,"Setting vars");
        ClearBigInt(OldRoot,Len);
        SetNum(AGM_A,Len,  BI_One,0);
        SetNum(AGM_A2,Len, BI_One,0);
        SetNum(AGM_Sum,Len,BI_One,0);
        if (!Cfg.Macintosh) BackSpace(12);
        Pass=1;
       }

     EndTime = clock();
     fprintf(stderr, "Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     DumpTimings(((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     if (Cfg.AlwaysSaveData)
        SaveData(AGM_A,AGM_A2,AGM_B,AGM_B2,AGM_Sum,OldRoot,
                 ((double)EndTime-(double)StartTime)/CLOCKS_PER_SEC,
	  Pow2,Pass,Len,Cfg.PiFormulaToUse);
    }

/*
  DumpBigInt("AGM_A  ",AGM_A,Len);
  DumpBigInt("AGM_A2 ",AGM_A2,Len);
  DumpBigInt("AGM_B  ",AGM_B,Len);
  DumpBigInt("AGM_B2 ",AGM_B2,Len);
  DumpBigInt("AGM_C2 ",AGM_C2,Len);
  DumpBigInt("AGM_Sum",AGM_Sum,Len);
*/

  Done=0;
  MaxPasses+=Pass;
  while ((!Done) && (Pass < MaxPasses))
    {size_t FNZ;
     int Key=0;
/* DJGPP stuff to test the self checking.
** Comment out the keyboard break statement at the end of the loop.
     extern int getkey(void);
     if (TestKeyboard()) {Key=getkey();srand((unsigned int)time(NULL));}
*/

      fprintf(stderr, "Pass %4s: ", Num2Str(2<<(++Pass)));
      LoopTime = clock();

      if ((Pass==1) &&
          (strcmp(GetCheckStr(AGM_A),"1000000000000000")==0))
        {
/*
** Since AGM_A==1, we can do the first pass a little differently
** and avoide the multiplication completely.
*/
         if (!Cfg.Macintosh) fprintf(stderr,"AGM Part1");

         /* Compute a=(a+b)/2 */
         if (Key=='1') CorruptVar(AGM_A,Len);
         Add(AGM_A, AGM_A, AGM_B, Len);DivBy(AGM_A,AGM_A, 2, Len);
         if (Cfg.AGMSelfCheck >= 1)
           {
            DivBy(DSWork,AGM_B,2,Len);AddInt(DSWork,BI_OneHalf);
            DoCheck(DSWork,AGM_A,1,Len);
           }

         /* Compute C^2 and B^2 */
         Add(AGM_C2,AGM_A2,AGM_B2,Len);DivBy(AGM_C2,AGM_C2,2,Len);
         if (Cfg.AGMSelfCheck >= 1) Add(DSWork,AGM_A2,AGM_B2,Len);
         Copy(AGM_B2,AGM_B,Len); /* AGM_A == 1.0, so we can copy */
         Sign=Sub(AGM_C2,AGM_C2,AGM_B2,Len);DivBy(AGM_C2,AGM_C2,2,Len);
         if (Key=='4') CorruptVar(AGM_B2,Len);
         if (Key=='5') CorruptVar(AGM_C2,Len);
         if (Sign) FatalError("AGM_C2 should never be negative.\n");
         if (!Cfg.Macintosh) BackSpace(9);
        }
      else
        {
         if (!Cfg.Macintosh) fprintf(stderr,"AGM Part1");

         /* Compute C */
/*
         Sign=Sub(AGM_C2, AGM_A, AGM_B, Len);DivBy(AGM_C2,AGM_C2,2,Len);
*/
         Sign=HalfDiff(AGM_C2,AGM_A,AGM_B,Len);
         if (Sign) FatalError("AGM_C2 should never be negative.\n");

         /* Compute A.  AGM_C2 is still only C at the moment. */
         if (Cfg.AGMSelfCheck >= 1)
           {Add(DSWork,AGM_A,AGM_B,Len);DivBy(DSWork,DSWork,2,Len);}
         if (Sub(AGM_A,AGM_A,AGM_C2,Len))
            FatalError("AGM_A should never be negative.\n");
         if (Key=='1') CorruptVar(AGM_A,Len);
         if (Cfg.AGMSelfCheck >= 1) DoCheck(DSWork,AGM_A,2,Len);

         if (!Cfg.Macintosh) BackSpace(9);
         if (!Cfg.Macintosh) fprintf(stderr,"Sqr");
         /* Compute C^2 */
         SpecialSquare(AGM_C2,AGM_C2,Len,AGM_B);
         if (Key=='5') CorruptVar(AGM_C2,Len);
         if (!Cfg.Macintosh) BackSpace(3);

         if (!Cfg.Macintosh) fprintf(stderr,"AGM Part2");
         if (Cfg.AGMSelfCheck >= 1) Add(DSWork,AGM_A2,AGM_B2,Len);
         /* Compute B^2 */
         Sign=SpecialAGMFunc1(AGM_B2,AGM_A2,AGM_C2,Len);
         if (Key=='4') CorruptVar(AGM_B2,Len);
/*
         Add(AGM_B2,AGM_B2,AGM_A2,Len);
         MulBy(AGM_A2,AGM_C2,4,Len);
         Sign=Sub(AGM_B2,AGM_B2,AGM_A2,Len);
         DivBy(AGM_B2,AGM_B2,2,Len);
*/

         if (Sign) FatalError("AGM_B2 should never be negative.\n");
         if (!Cfg.Macintosh) BackSpace(9);
        }

      if (!Cfg.Macintosh) fprintf(stderr,"AGM Part3");
      /* Compute A^2 */
      if (Cfg.AGMSelfCheck >= 1)
        {
         Add(DSWork,DSWork,AGM_B2,Len);
         Add(DSWork,DSWork,AGM_B2,Len);
         DivBy(DSWork,DSWork,4,Len);
        }
      Add(AGM_A2,AGM_C2,AGM_B2,Len);
      if (Key=='2') CorruptVar(AGM_A2,Len);
      if (Cfg.AGMSelfCheck >= 1) DoCheck(DSWork,AGM_A2,3,Len);
      if (!Cfg.Macintosh) BackSpace(9);

/*
** Do some self checking.  The estimate formula predicts the number of
** leading zeros.  It's based on the regular AGM accuracy formula.  It
** predicts just a couple of digits less than what is actually there,
** so if anything happens, this should fail.  But, since it is so
** close, it just might fail for the wrong reason.
*/
      if (!Cfg.Macintosh) fprintf(stderr,"Self Check");
      FNZ=FindFirstNonZero(AGM_C2,Len);
      if ((Pass > 3) && (FNZ != Len))
        {size_t P=Pass-2;size_t Est;
         /* formula based on the AGM 'number of digits right' formula */
         if (Cfg.PiFormulaToUse==2)
            Est=(1.364376354*(1<<(P))-P*.301029995-2.342434419)/2;
         else /* other agm does 1 less pass, so pass is one higher, so div by 4 */
            Est=(sqrt(3.0)*1.364376354*(1<<(P))-P*.301029995-2.690531961)/4;

         if (Est >= FNZ)
           {
            fprintf(stderr,"\aAGM self check fails.  Predicted: %lu Actual: %lu\n",
                    (ULINT)Est,(ULINT)FNZ);
            fprintf(stderr,"Beyond what I've tested, this may not be an actual failure\n");
            fprintf(stderr,"but an inaccuracy in the self check formula itself.\n");
           }
        }

      /* See if it's done.  (ie: more than half zeros) */
      if (FNZ > Len/2+4) Done=1;

      if ((Cfg.AGMSelfCheck == 2) || (Cfg.AGMSelfCheck == 4))
        {
      /* We can use AGM_B and DSWork as scratch. */
         SpecialSquare(AGM_B,AGM_A,Len,DSWork);
         DoCheck(AGM_B,AGM_A2,4,Len);
        }
      if (!Cfg.Macintosh) BackSpace(10);

      /* Compute B=sqrt(B^2) */
      if (!Done)
        AGMSqrt(AGM_B, AGM_B2, Len,  (Pass>5) ? (2<<(Pass-5)) : 0);
      if (Key=='3') CorruptVar(AGM_B,Len);

      if (!Cfg.Macintosh) fprintf(stderr,"AGM Part4");
      /* Sum = Sum - C2 * 2^(J+1) */
      MulByFloat(AGM_C2,Pow2,Len);Pow2 *= 2.0;
      if (Sub(AGM_Sum,AGM_Sum,AGM_C2,Len))
        FatalError("AGM_Sum should never be negative.\n");

      EndTime=clock();
/*
      DumpBigInt("AGM_A  ",AGM_A,Len);
      DumpBigInt("AGM_A2 ",AGM_A2,Len);
      DumpBigInt("AGM_B  ",AGM_B,Len);
      DumpBigInt("AGM_B2 ",AGM_B2,Len);
      DumpBigInt("AGM_C2 ",AGM_C2,Len);
      DumpBigInt("AGM_Sum",AGM_Sum,Len);
*/

      if (!Cfg.Macintosh) BackSpace(9);
      if (!Done && ((Cfg.AGMSelfCheck == 3) || (Cfg.AGMSelfCheck == 4)))
        {
      /* We can use AGM_C2 and DSWork as scratch. */
         if (!Cfg.Macintosh) fprintf(stderr,"Self Check");
         SpecialSquare(AGM_C2,AGM_B,Len,DSWork);
         DoCheck(AGM_C2,AGM_B2,5,Len);
         if (!Cfg.Macintosh) BackSpace(10);
        }

      DumpDebug("Pass %u took:", Pass);
      fprintf(stderr,"Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
      DumpTimings(((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);

      if (Cfg.AlwaysSaveData)
        SaveData(AGM_A,AGM_A2,AGM_B,AGM_B2,AGM_Sum,OldRoot,
                 ((double)EndTime-(double)StartTime)
	     /CLOCKS_PER_SEC,Pow2,Pass,Len,Cfg.PiFormulaToUse);
      if (TestKeyboard()) break;
    }


/* We've done the AGM part, now do the final calculation. */
  LoopTime=clock();
  if (Done)
    {
     fprintf(stderr, "Final    : ");
/*
** Since AGM_A and AGM_B have converged, I can compute the AGM_A2
** by calculating AGM_B2.  That can be done like up above, except
** we don't need the 4*AGM_C2 because that will be zero.
*/
     if (!Cfg.Macintosh) fprintf(stderr,"Part 1");
     Add(AGM_B2,AGM_B2,AGM_A2,Len);
     MulBy(AGM_B2,AGM_B2,2,Len);
     if (!Cfg.Macintosh) BackSpace(6);

     if (Cfg.PiFormulaToUse==3)
       {
        if (!Cfg.Macintosh) fprintf(stderr,"Part 2");
        if (!Cfg.Macintosh) BackSpace(6);
        Sqrt30(AGM_A2,Len);
        if (!Cfg.Macintosh) fprintf(stderr,"Part 3");
        SpecialFullMul(AGM_A,AGM_Sum,AGM_A2,Len,AGM_B);
        SubInt(AGM_A, BI_OneHalf);
        Copy(AGM_Sum,AGM_A,Len);
        if (!Cfg.Macintosh) BackSpace(6);
       }

     AGMDivide(AGM_B, AGM_B2, AGM_Sum, Len, AGM_A2);

     EndTime=clock();
     fprintf(stderr, "Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     DumpTimings(((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);

     if (Cfg.PiFormulaToUse==2)
        PrintFormattedPi("Fast 1/sqrt(2) AGM",((double)EndTime-(double)StartTime)
	     /CLOCKS_PER_SEC, AGM_B, Len);
     else
        PrintFormattedPi("Fast sqrt(3) AGM",((double)EndTime-(double)StartTime)
	     /CLOCKS_PER_SEC, AGM_B, Len);

     DeleteSaveFile();
    }
  else if (Pass >= FatalPasses)
       FatalError("The AGM didn't converge.\n");
  else SaveData(AGM_A,AGM_A2,AGM_B,AGM_B2,AGM_Sum,OldRoot,
                ((double)clock()-(double)StartTime)
	     /CLOCKS_PER_SEC,Pow2,Pass,Len,Cfg.PiFormulaToUse);

  fprintf(stderr,"\nTotal Execution time: %0.2f\n\n",((double)clock()-(double)StartTime)
	     /CLOCKS_PER_SEC);
  if (!Done) DumpTimings(((double)clock()-(double)StartTime)
	     /CLOCKS_PER_SEC);

  return (Done);
}

/*
*****************************************************
**      The AGM itself                             **
**                                                 **
** A[n] = (A[n-1] + B[n-1])/2                      **
** B[n] = Sqrt(A[n-1]*B[n-1])                      **
** C[n] = (A[n-1]-B[n-1])/2                        **
**                                                 **
**                          n                      **
** PI[n] = 4A[n+1]^2 / (1-(Sum (2^(j+1))*C[j]^2))  **
**                        j = 1                    **
*****************************************************
*/
static int
ComputeSlowAGM(size_t Len, size_t MaxPasses)
{
  int Sign;
  size_t Pass;
  double Pow2;
  clock_t LoopTime,StartTime,EndTime;
  BigInt AGM_A, AGM_B, AGM_C, AGMWork;
  int Done; /* boolean */

  if ((Cfg.PiFormulaToUse != 2) && (Cfg.PiFormulaToUse != 3))
    FatalError("SlowAGM was somehow called with pi formula %d\n",
               Cfg.PiFormulaToUse);

  AGM_A = CreateBigInt(Len);
  AGM_B = CreateBigInt(Len);
  AGM_C = CreateBigInt(Len);
  AGMWork = CreateBigInt(Len);

  Num1IsCached = Num2IsCached = 0;
  StartTime = clock();
  Pow2 = 4.0;
  Pass = 0;

  if (!LoadData(AGM_A,AGM_B,AGM_C,OldRoot,NO_NUM,NO_NUM,
                &StartTime,&Pow2,&Pass,Len,Cfg.PiFormulaToUse))
    {
     LoopTime = clock();
     fprintf(stderr, "Init     : ");

     if (Cfg.PiFormulaToUse==2)
       {
        SetNum(AGM_A,Len,BI_One,0);
        SetNum(AGM_C,Len,BI_One,0);
        ClearBigInt(OldRoot,Len);
        Sqrt05(AGM_B,Len);
       }
     else
       {
        Sqrt20(AGM_A,Len);
        Sqrt60(AGM_C,Len);
        if (!Cfg.Macintosh) fprintf(stderr,"Square");
        Add(AGM_B,AGM_A,AGM_C,Len);
        DivBy(AGM_B,AGM_B,4,Len);
        if (!Cfg.Macintosh) BackSpace(6);
        if (!Cfg.Macintosh) fprintf(stderr,"Setting vars");
        ClearBigInt(OldRoot,Len);
        SetNum(AGM_A,Len,BI_One,0);
        SetNum(AGM_C,Len,BI_One,0);
        if (!Cfg.Macintosh) BackSpace(12);
        Pass=1;
       }

     EndTime = clock();
     fprintf(stderr,"Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     DumpTimings(((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     if (Cfg.AlwaysSaveData)
        SaveData(AGM_A,AGM_B,AGM_C,OldRoot,NO_NUM,NO_NUM,
                 ((double)EndTime-(double)StartTime)
     /CLOCKS_PER_SEC,Pow2,Pass,Len,Cfg.PiFormulaToUse);
    }

/*
  DumpBigInt("AGM_A",AGM_A,Len);
  DumpBigInt("AGM_B",AGM_B,Len);
  DumpBigInt("AGM_C",AGM_C,Len);
*/

  Done=0;
  MaxPasses+=Pass;
  while ((!Done) && (Pass < MaxPasses))
    {
      fprintf(stderr, "Pass %4s: ", Num2Str(1<<(++Pass)));
      LoopTime = clock();

      if (!Cfg.Macintosh) fprintf(stderr,"AGM Part1");
      /* w = (a-b)/2 */
      Sign = Sub(AGMWork, AGM_A, AGM_B, Len);
      if (Sign) Negate(AGMWork,Len);
      DivBy(AGMWork,AGMWork, 2, Len);

      {size_t Nz;
       Nz=FindFirstNonZero(AGMWork,Len);
       if ((Pass > 4) && (Nz != Len))
         if ((1<<(Pass-4)) > Nz)
           fprintf(stderr,"AGMCheck1 fails.  Predicted: %lu Actual: %lu\a\n",
                   (unsigned long)1<<(Pass-4),(ULINT)Nz);
      }
      if (!Cfg.Macintosh) BackSpace(9);

      /* m = w*w */
      if (!Cfg.Macintosh) fprintf(stderr,"Sqr1");
      if (IsZero(AGMWork,Len/2)) Done=1;
      if (Done) ClearBigInt(AGMWork,Len);
      else FullMul(AGMWork, AGMWork, AGMWork, Len);
      if (!Cfg.Macintosh) BackSpace(4);

      if (!Cfg.Macintosh) fprintf(stderr,"AGM Part2");
      {size_t Nz;
      Nz=FindFirstNonZero(AGMWork,Len);
      if ((Pass > 3) && (Nz != Len))
        if ((1<<(Pass-3)) > Nz)
          fprintf(stderr,"AGMCheck2 fails.  Predicted: %lu Actual: %lu\a\n",
                  (unsigned long)1<<(Pass-3),(ULINT)Nz);
      }

      /* m = m* w^(J+1) */
      MulByFloat(AGMWork,Pow2,Len);
      Pow2 *= 2.0;

      /* c = c - m */
      if (Sign) Add(AGM_C, AGM_C, AGMWork, Len);
      else      Sub(AGM_C, AGM_C, AGMWork, Len);

      /* See if it's done */
      if (IsZero(AGMWork,Len-(Len/16))) Done=1;
      if (!Cfg.Macintosh) BackSpace(9);

      /* m = a*b */
      if (!Cfg.Macintosh) fprintf(stderr,"Sqr2");
      if (Pass==1) Copy(AGMWork,AGM_B,Len); /* first pass, AGM_A = 1.0 */
      else if (!Done)   /* If last iter, we can skip it. */
         FullMul(AGMWork, AGM_A, AGM_B, Len);
      if (!Cfg.Macintosh) BackSpace(4);

      /* a = (a+b)/2 */
      Add(AGM_A, AGM_A, AGM_B, Len);
      DivBy(AGM_A,AGM_A, 2, Len);

      /* b = sqrt(a*b) */
      if (!Done) /* Optimization */
          AGMSqrt(AGM_B, AGMWork, Len, (Pass>5) ? (2<<(Pass-5)) : 0 );

      EndTime=clock();

/*
      DumpBigInt("AGM_A",AGM_A,Len);
      DumpBigInt("AGM_B",AGM_B,Len);
      DumpBigInt("AGM_C",AGM_C,Len);
*/
      DumpDebug("Pass %u took:", Pass);
      fprintf(stderr, "Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
     /CLOCKS_PER_SEC);
      DumpTimings(((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);

      if (TestKeyboard()) break;
     if (Cfg.AlwaysSaveData)
        SaveData(AGM_A,AGM_B,AGM_C,OldRoot,NO_NUM,NO_NUM,
                 ((double)EndTime-(double)StartTime)
	     /CLOCKS_PER_SEC,Pow2,Pass,Len,Cfg.PiFormulaToUse);
    }

  LoopTime=clock();
  if (Done)
    {
     fprintf(stderr,"Final    : ");
     if (!Cfg.Macintosh) fprintf(stderr,"Sqr");
     SpecialSquare(AGM_A, AGM_A, Len, AGMWork);
     MulBy(AGM_A,AGM_A,4,Len);
     if (!Cfg.Macintosh) BackSpace(3);

     if (Cfg.PiFormulaToUse==3)
       {
        Sqrt30(AGM_B,Len);
        if (!Cfg.Macintosh) fprintf(stderr,"Part 2");
        FullMul(AGM_C,AGM_C,AGM_B,Len);
        SubInt(AGM_C, BI_OneHalf);
        if (!Cfg.Macintosh) BackSpace(6);
       }

     AGMDivide(AGM_B, AGM_A, AGM_C, Len, AGMWork);
     EndTime=clock();
     fprintf(stderr, "Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     DumpTimings(((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     if (Cfg.PiFormulaToUse==2)
        PrintFormattedPi("Slow 1/sqrt(2) AGM",((double)EndTime-(double)StartTime)
	     /CLOCKS_PER_SEC, AGM_B, Len);
     else
        PrintFormattedPi("Slow sqrt(3) AGM",((double)EndTime-(double)StartTime)
	     /CLOCKS_PER_SEC, AGM_B, Len);
     DeleteSaveFile();
    }
  else if (Pass >= FatalPasses)
       FatalError("The AGM didn't converge.\n");
  else SaveData(AGM_A,AGM_B,AGM_C,OldRoot,NO_NUM,NO_NUM,
                ((double)clock()-(double)StartTime)
	     /CLOCKS_PER_SEC,Pow2,Pass,Len,Cfg.PiFormulaToUse);

  fprintf(stderr,"\nTotal Execution time: %0.2f seconds.\n\n",((double)clock()-(double)LoopTime)
	     /CLOCKS_PER_SEC);
  if (!Done) DumpTimings(((double)clock()-(double)LoopTime)
	     /CLOCKS_PER_SEC);

  return (Done);
}

static int
ComputeClassicAGM(size_t Len, size_t MaxPasses)
{
  int Sign;
  size_t Pass;
  double Pow2;
  clock_t LoopTime,EndTime,StartTime;
  BigInt AGM_A, AGM_B, AGM_C, AGMWork;
  int Done; /* boolean */

  AGM_A = CreateBigInt(Len);
  AGM_B = CreateBigInt(Len);
  AGM_C = CreateBigInt(Len);
  AGMWork = OldRoot;

  StartTime = clock();
  Pow2 = 4.0;
  Pass = 0;

  if (!LoadData(AGM_A,AGM_B,AGM_C,NO_NUM,NO_NUM,NO_NUM,
                &StartTime,&Pow2,&Pass,Len,Cfg.PiFormulaToUse))
    {
     fprintf(stderr, "Init     : ");
     LoopTime = clock();
     SetNum(AGM_A,Len,BI_One,0);
     SetNum(AGM_C,Len,BI_OneHalf,0);
     ClassicSqrt(AGM_B,AGM_C,Len);
     SetNum(AGM_C,Len,BI_One,0);
     EndTime = clock();
     fprintf(stderr,"Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     DumpTimings(((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     if (Cfg.AlwaysSaveData)
        SaveData(AGM_A,AGM_B,AGM_C,NO_NUM,NO_NUM,NO_NUM,
                 ((double)EndTime-(double)StartTime)
	     /CLOCKS_PER_SEC,Pow2,Pass,Len,Cfg.PiFormulaToUse);
    }

/*
  DumpBigInt("AGM_A",AGM_A,Len);
  DumpBigInt("AGM_B",AGM_B,Len);
  DumpBigInt("AGM_C",AGM_C,Len);
*/

  Done=0;
  MaxPasses+=Pass;
  while ((!Done) && (Pass < MaxPasses))
    {
      fprintf(stderr, "Pass %4s: ", Num2Str(1<<(++Pass)));
      LoopTime = clock();

      if (!Cfg.Macintosh) fprintf(stderr,"AGM Part1");
      /* w = (a-b)/2 */
      Sign = Sub(AGMWork, AGM_A, AGM_B, Len);
      if (Sign) Negate(AGMWork,Len);
      DivBy(AGMWork,AGMWork, 2, Len);
      if (!Cfg.Macintosh) BackSpace(9);

      /* m = w*w */
      if (!Cfg.Macintosh) fprintf(stderr,"Sqr1");
      if (IsZero(AGMWork,Len/2)) Done=1;
      if (Done) ClearBigInt(AGMWork,Len);
      else FullMul(AGMWork, AGMWork, AGMWork, Len);
      if (!Cfg.Macintosh) BackSpace(4);

      if (!Cfg.Macintosh) fprintf(stderr,"AGM Part2");

      /* m = m* w^(J+1) */
      MulByFloat(AGMWork,Pow2,Len);
      Pow2 *= 2.0;

      /* c = c - m */
      if (Sign) Add(AGM_C, AGM_C, AGMWork, Len);
      else      Sub(AGM_C, AGM_C, AGMWork, Len);

      /* See if it's done */
      if (IsZero(AGMWork,Len-(Len/16))) Done=1;
      if (!Cfg.Macintosh) BackSpace(9);

      /* m = a*b */
      if (!Cfg.Macintosh) fprintf(stderr,"Sqr2");
      if (Pass==1) Copy(AGMWork,AGM_B,Len); /* first pass, AGM_A = 1.0 */
      else if (!Done)   /* If last iter, we can skip it. */
         FullMul(AGMWork, AGM_A, AGM_B, Len);

      /* a = (a+b)/2 */
      Add(AGM_A, AGM_A, AGM_B, Len);
      DivBy(AGM_A,AGM_A, 2, Len);
      if (!Cfg.Macintosh) BackSpace(4);

      /* b = sqrt(a*b) */
      if (!Done) /* Optimization */
          ClassicSqrt(AGM_B, AGMWork, Len);

      EndTime=clock();

/*
      DumpBigInt("AGM_A",AGM_A,Len);
      DumpBigInt("AGM_B",AGM_B,Len);
      DumpBigInt("AGM_C",AGM_C,Len);
*/
      DumpDebug("Pass %u took:", Pass);
      fprintf(stderr, "Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
      DumpTimings(((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);

      if (TestKeyboard()) break;
     if (Cfg.AlwaysSaveData)
        SaveData(AGM_A,AGM_B,AGM_C,NO_NUM,NO_NUM,NO_NUM,
                 ((double)EndTime-(double)StartTime)
	     /CLOCKS_PER_SEC,Pow2,Pass,Len,Cfg.PiFormulaToUse);
    }

  LoopTime=clock();
  if (Done)
    {
     fprintf(stderr,"Final    : ");
     if (!Cfg.Macintosh) fprintf(stderr,"Sqr");
     FullMul(AGM_A, AGM_A, AGM_A, Len);
     MulBy(AGM_A,AGM_A,4,Len);
     if (!Cfg.Macintosh) BackSpace(3);

     ClassicDivide(AGM_B, AGM_A, AGM_C,Len);
     EndTime=clock();
     fprintf(stderr, "Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     DumpTimings(((double)EndTime-(double)LoopTime)
	     /CLOCKS_PER_SEC);
     PrintFormattedPi("Classic AGM",((double)EndTime-(double)StartTime)
	     /CLOCKS_PER_SEC, AGM_B, Len);
     DeleteSaveFile();
    }
  else if (Pass >= FatalPasses)
       FatalError("The AGM didn't converge.\n");
  else SaveData(AGM_A,AGM_B,AGM_C,NO_NUM,NO_NUM,NO_NUM,
                ((double)clock()-(double)StartTime)
	     /CLOCKS_PER_SEC,Pow2,Pass,Len,Cfg.PiFormulaToUse);

  fprintf(stderr,"\nTotal Execution time: %0.2f seconds.\n\n",((double)clock()-(double)StartTime)
	     /CLOCKS_PER_SEC);
  if (!Done) DumpTimings(((double)clock()-(double)StartTime)
	     /CLOCKS_PER_SEC);

  return (Done);
}



int
ComputeAGM(size_t Len, size_t MaxPasses)
{
if (Cfg.PiFormulaToUse==1)
                    return ComputeClassicAGM(Len,MaxPasses);

if (Cfg.UseFastAGM) return ComputeFastAGM(Len,MaxPasses);
else                return ComputeSlowAGM(Len,MaxPasses);
}


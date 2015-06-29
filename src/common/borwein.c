#include "pi.h"
#include "borwein.h"
#include "bigint.h"
#include "cache.h"
#include "bigmul.h"

/*
** d = a/b by computing the reciprocal of b and then multiplying
** that by a.
*/
static void
BorDivide(BigInt R, BigInt Num1, BigInt Num2, BigInt Work,
          size_t Len, size_t SubLen)
{ int Sign;
  size_t Redo=REDO_LEN;
  int FirstPass=0;

  if (SubLen <= 0) SubLen = 2;
  if (SubLen > Len) SubLen = Len;

/* There are only a few divisors, so I can cheat. */
  if      (NumIs(Num2,20074977,74254721)) SetNum(R,Len, 4981325,57268337);
  else if (NumIs(Num2,20000000,   48646)) SetNum(R,Len, 4999999,99987838);
  else if (NumIs(Num2,20000000,       0))
    {
     SetNum(R,Len, BI_OneHalf,0);
     FirstPass=1;
    }
  else
    {
     DumpBigInt("Unknown Divisor: ",Num2,4);
     ExitPrg(EXIT_FAILURE);
    }

  if (SubLen >= Redo) Redo=0;
  if (SubLen*2 == Redo) Redo=0;
  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
  while (SubLen < Len/2)
    {
      SubLen *= 2;if (SubLen > Len) SubLen = Len;
      if (!Cfg.Macintosh) 
        fprintf(stderr,"Div: %4s",Num2Str(SubLen*RawIntDigits));
      FlushFFTCache(0);

/* Perform safety check */
      {char *Str=GetCheckStr(R);
       if ( (strcmp(Str,"0498132557268337")!=0) &&
            (strcmp(Str,"0499999999987838")!=0) &&
            (strcmp(Str,"0499999999999999")!=0) &&
            (strcmp(Str,"0500000000000000")!=0) )
         fprintf(stderr,"** WARNING **\a\nBorDivide may be failing.\n%s\n",Str);
      }

      ClearBigInt(R+SubLen/2,SubLen/2);
      if (!Cfg.Macintosh) fputc('|',stderr);
      if (FirstPass) DivBy(DSWork,Num2,2,SubLen);
      else {
            SaveNum1FFT = 50;
            if (SubLen==Len/2) SaveNum2FFT=51;
            N1R0Mul(DSWork, R, Num2, Work, SubLen);
           }

      Sign = RevSubInt(BI_One,DSWork,SubLen);

      if (!Cfg.Macintosh) fputc('.',stderr);
      if (FirstPass) {DivBy(DSWork,DSWork,2,SubLen);}
      else {
            Num1IsCached=50;
            HalfMul(DSWork,R,DSWork,SubLen);
           }
      if (Sign) Sub(R,R,DSWork,SubLen);
      else      Add(R,R,DSWork,SubLen);

      if (Redo == ULONG_MAX) Redo=0;
      if (SubLen == Redo) {SubLen/=2;Redo=ULONG_MAX;}
      FirstPass=0;
      if (!Cfg.Macintosh) BackSpace(11);
    }

  if (!Cfg.Macintosh) fprintf(stderr,"Div: %4s",Num2Str(Len*RawIntDigits));
  FlushFFTCache(51);
  ClearBigInt(R+Len/2,Len/2);
  if (!Cfg.Macintosh) fputc('.',stderr);
  if (FirstPass) DivBy(OldRoot,Num1,2,Len);
  else {
        SaveNum1FFT = 50;
        HalfMul(OldRoot,R,Num1,Len);
       }

  ClearBigInt(OldRoot+Len/2,Len/2);
  if (!Cfg.Macintosh) fputc('|',stderr);
  if (FirstPass) {Add(DSWork,OldRoot,OldRoot,Len);}
  else {
        Num2IsCached=51;
        N1R0Mul(DSWork,OldRoot,Num2,Work,Len);
       }

  Sign=Sub(DSWork,Num1,DSWork,Len);
  if (Sign) Negate(DSWork,Len);

  if (!Cfg.Macintosh) fputc('.',stderr);
  if (FirstPass) {DivBy(DSWork,DSWork,2,Len);}
  else {
        Num1IsCached=50;
        HalfMul(DSWork,R,DSWork,Len);
       }
  if (Sign) Sub(R,OldRoot,DSWork,Len);
  else      Add(R,OldRoot,DSWork,Len);
  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
  if (!Cfg.Macintosh) BackSpace(12);
}

static void
BorReciprocal(BigInt R, BigInt Num2, BigInt Work, size_t Len)
{ int Sign;
  size_t Redo=REDO_LEN;
  size_t SubLen=2;

/* There are only a few divisors, so I can cheat. */
  if (NumIs(Num2, 3183098,86183790)) SetNum(R,Len,31415926,53589793);
  else
    {
     DumpBigInt("Unknown Recip: ",Num2,4);
     ExitPrg(EXIT_FAILURE);
    }

  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
  while (SubLen < Len/2)
    {
      SubLen *= 2;if (SubLen > Len) SubLen = Len;
      FlushFFTCache(0);
      if (!Cfg.Macintosh) fprintf(stderr,"Recip: %4s",Num2Str(SubLen*RawIntDigits));

/* Perform safety check */
      if (strcmp(GetCheckStr(R),"3141592653589793")!=0)
        fprintf(stderr,"** WARNING **\a\nBorRecip may be failing.\n%s\n",GetCheckStr(R));

      ClearBigInt(R+SubLen/2,SubLen/2);
      SaveNum1FFT = 60;
      if (SubLen==Len/2) SaveNum2FFT=61;
      if (!Cfg.Macintosh) fputc('|',stderr);
      N1R0Mul(DSWork, R, Num2, Work, SubLen);

      Sign = RevSubInt(BI_One,DSWork,SubLen);

      Num1IsCached=60;
      if (!Cfg.Macintosh) fputc('.',stderr);
      HalfMul(DSWork,R,DSWork,SubLen);
      if (Sign) Sub(R,R,DSWork,SubLen);
      else      Add(R,R,DSWork,SubLen);

      if (Redo == ULONG_MAX) Redo=0;
      if (SubLen == Redo) {SubLen/=2;Redo=ULONG_MAX;}
      if (!Cfg.Macintosh) BackSpace(13);
    }

  if (!Cfg.Macintosh) fprintf(stderr,"Recip: %4s",Num2Str(Len*RawIntDigits));
  FlushFFTCache(61);
  ClearBigInt(R+Len/2,Len/2);
  Copy(OldRoot,R,Len);

  if (!Cfg.Macintosh) fputc('|',stderr);
  ClearBigInt(OldRoot+Len/2,Len/2);
  Num2IsCached=61;
  N1R0Mul(DSWork,OldRoot,Num2,Work,Len);

  Sign = RevSubInt(BI_One,DSWork,Len);
  if (Sign) Negate(DSWork,Len);

  if (!Cfg.Macintosh) fputc('.',stderr);
  HalfMul(DSWork,R,DSWork,Len);
  if (Sign) Sub(R,OldRoot,DSWork,Len);
  else      Add(R,OldRoot,DSWork,Len);
  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
  if (!Cfg.Macintosh) BackSpace(13);
}

static void
Sqrt20(BigInt Root, size_t Len)
{
  int Sign;
  size_t Redo=REDO_LEN;
  size_t SubLen=2;

  SetNum(Root,4,14142135,62373095);
  ClearBigInt(Root+SubLen,Len-SubLen);

  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
  while (SubLen < Len)
    {
      SubLen *= 2;if (SubLen > Len) SubLen = Len;
      if (!Cfg.Macintosh) fprintf(stderr,"Sqrt(2): %4s",Num2Str(SubLen*RawIntDigits));
      FlushFFTCache(0);

/* Perform safety check */
      if (strcmp(GetCheckStr(Root),"1414213562373095")!=0)
        fprintf(stderr,"** WARNING **\a\nSqrt20 may be failing.\n%s\n",GetCheckStr(Root));

      ClearBigInt(Root+SubLen/2,SubLen/2);
      SaveNum1FFT = 70;
      if (!Cfg.Macintosh) fputc('.',stderr);
      HalfMul(DSWork, Root, Root, SubLen);
      DivBy(DSWork,DSWork,2,SubLen);

      Sign = RevSubInt(BI_One,DSWork,SubLen);

      if (!Cfg.Macintosh) fputc('|',stderr);
      Num1IsCached=70;
      HalfMul(DSWork,Root,DSWork,SubLen);
      DivBy(DSWork,DSWork,2,SubLen);

      if (Sign) Sub(Root,Root,DSWork,SubLen);
      else      Add(Root,Root,DSWork,SubLen);

      if (Redo == ULONG_MAX) Redo = 0;
      if (SubLen == Redo) {SubLen/=2;Redo=ULONG_MAX;}
      if (!Cfg.Macintosh) BackSpace(15);
    }

  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
}

/*
** Calculate the 4th root of a number.  Same basic method as Sqrt()
** except it's easier to leave it in reciprocal form.
*/
static void
FourthRoot(BigInt Root, BigInt Num, size_t Len, size_t SubLen)
{
  int Sign;
  size_t Redo=REDO_LEN;
  int FirstPass=0;

  if (SubLen <= 0) SubLen = 2;
  if (SubLen > Len) SubLen = Len;

  if      (NumIs(Num, 9705627,48477140)) SetNum(Root,Len,10074977,74254721);
  else if (NumIs(Num, 9999999,99805415)) SetNum(Root,Len,10000000,   48646);
  else if (NumIs(Num, 9999999,99999999))
    {
     SetNum(Root,Len,BI_One,0);
     FirstPass=1;
    }
  else
    {
     DumpBigInt("Unknown 4th Root: ",Num,4);
     ExitPrg(EXIT_FAILURE);
    }

  if (SubLen >= Redo) Redo=0;
  if (SubLen*2 == Redo) Redo=0;
  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
  while (SubLen < Len)
    {
      SubLen *= 2;if (SubLen > Len) SubLen = Len;
      FlushFFTCache(0);
      if (!Cfg.Macintosh) fprintf(stderr,"4th: %4s",Num2Str(SubLen*RawIntDigits));

/* Perform safety check */
      {char *Str=GetCheckStr(Root);
       if ( (strcmp(Str,"1007497774254721")!=0) &&
            (strcmp(Str,"1000000000048646")!=0) &&
            (strcmp(Str,"1000000000000000")!=0) )
         fprintf(stderr,"** WARNING **\a\nFourthRoot may be failing.\n%s\n",Str);
      }

      ClearBigInt(Root+SubLen/2,SubLen/2);
      if (!Cfg.Macintosh) fputc('#',stderr);
      if (FirstPass) { }
      else {
            SaveNum1FFT=80;
            HalfMul(DSWork, Root, Root, SubLen);
            FullMul(DSWork, DSWork, DSWork, SubLen);
           }

      if (!Cfg.Macintosh) fputc('=',stderr);
      if (FirstPass) {Copy(DSWork,Num,SubLen);}
      else FullMul(DSWork, Num, DSWork, SubLen);

      Sign=RevSubInt(BI_One,DSWork,SubLen);

      if (!Cfg.Macintosh) fputc('.',stderr);
      if (!FirstPass)
        {
         Num1IsCached=80;
         HalfMul(DSWork,Root,DSWork,SubLen);
        }
      DivBy(DSWork,DSWork, 4, SubLen);

      if (Sign) Sub(Root,Root,DSWork,SubLen);
      else      Add(Root,Root,DSWork,SubLen);

      if (Redo == ULONG_MAX) Redo = 0;
      if (SubLen == Redo) {SubLen/=2;Redo=ULONG_MAX;}
      FirstPass=0;
      if (!Cfg.Macintosh) BackSpace(12);
    }

  Num1IsCached = Num2IsCached = 0;
  FlushFFTCache(0);
}

/*
** The Borwein-Borwein-Bailey quartic (4) pi program.
**
** Root() means the 4th root.  ie: sqrt(sqrt(x))
**
** Set:
**
** a=6-4*sqrt(2)
** b=sqrt(2)-1
**
** Iterate:
**
**    1-Root(1-b^4)
** b= -------------
**    1+Root(1-b^4)
**
** a=(1+b)^4*a-2^(2n+3)*b*(1+b+b^2)
**
** and 'a' converges to 1/pi
**
** 'b' must be computed before 'a'.
**
** All of those powers can be done fairly efficiently by raising
** b to 2, 3, and the 4th power.  Even the (1+b)^4 decomposes into
** powers of b, so we don't need to do a special 4th power.
**
** a=(1+b)^4*a-2^(2n+3)*b*(1+b+b*b)
**
** (1+b)^4     = (1+4b+6b^2+4b^3+b^4)
** b*(1+b+b^2) = (b+b^2+b^3)
**
** b^2 = b*b     FU   (b   -> cache)
** b^3 = b*b^2   CFU  (b^2 -> cache)
** b^4 = b^2*b^2 CCU
**
*/
int
ComputeBorwein(size_t Len, size_t MaxPasses)
{double Pow2;
 clock_t LoopTime,StartTime,EndTime;
 BigInt VAR_A,VAR_B,Work1,Work2,Work3;
 size_t Pass;
 int Done;

VAR_A = CreateBigInt(Len);
VAR_B = CreateBigInt(Len);
Work1 = CreateBigInt(Len);
Work2 = CreateBigInt(Len);
Work3 = CreateBigInt(Len);

Num1IsCached = Num2IsCached = 0;
Pow2=8;
Pass=0;
StartTime = clock();

/* a=6-4*sqrt(2) b=sqrt(2)-1 */
if (!LoadData(VAR_A,VAR_B,Work3,NO_NUM,NO_NUM,NO_NUM,
              &StartTime,&Pow2,&Pass,Len,2))
  {
   fprintf(stderr, "Init     : ");
   LoopTime = clock();
   Sqrt20(Work1,Len);
   SetNum(VAR_A,Len,6*BI_One,0);
   Sub(VAR_A,VAR_A,Work1,Len);Sub(VAR_A,VAR_A,Work1,Len);
   Sub(VAR_A,VAR_A,Work1,Len);Sub(VAR_A,VAR_A,Work1,Len);
   /* 17-12*sqrt(2) is same as (sqrt(2)-1)^4, but a lot faster */
   SetNum(Work3,Len,9*BI_One,0);
   MulBy(Work1,Work1,6,Len);
   Sub(Work3,Work3,Work1,Len);
   AddInt(Work3,8*BI_One);/* We only have one digit (0-9) of integer... */
   Sub(Work3,Work3,Work1,Len);
   EndTime = clock();
   fprintf(stderr,"Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	   /CLOCKS_PER_SEC);
   Pow2=8;
   Pass=0;
   DumpTimings(LoopTime);
   if (Cfg.AlwaysSaveData)
     SaveData(VAR_A,VAR_B,Work3,NO_NUM,NO_NUM,NO_NUM,
              ((double)EndTime-(double)StartTime)
	   /CLOCKS_PER_SEC,Pow2,Pass,Len,2);
  }

MaxPasses+=Pass;
Done=0;
while ((!Done) && (Pass < MaxPasses))
  {
   if (IsZero(Work3,Len)) Done=1;
   if (Done) break;

   Pass++;
   fprintf(stderr, "Pass %4s: ", Num2Str(1<<(Pass*2)));
   LoopTime = clock();

/*
** b= (1-Root(1-b^4)) / (1+Root(1-b^4)) Work3 is already b^4
** Slightly different because FourthRoot returns reciprocal, not root.
** b= (Root(1-b^4)-1) / (Root(1-b^4)+1)
*/
   SetNum(Work2,Len,BI_One,0);
   Sub(Work1,Work2,Work3,Len);
   FourthRoot(Work2,Work1,Len, (Pass>2) ? (1<<(Pass*2-4)) : 0 );
   Copy(Work1,Work2,Len);
   SubInt(Work1,BI_One);
   AddInt(Work2,BI_One);
   BorDivide(VAR_B,Work1,Work2,Work3,Len, (Pass>2) ? (1<<(Pass*2-4)) : 0 );
   FlushFFTCache(0);

/* Work1=1+4*b */
   if (!Cfg.Macintosh) fprintf(stderr,"Part1");
   SetNum(Work1,Len,BI_One,0);
   Add(Work1,Work1,VAR_B,Len);Add(Work1,Work1,VAR_B,Len);
   Add(Work1,Work1,VAR_B,Len);Add(Work1,Work1,VAR_B,Len);
/* Work2=b */
   Copy(Work2,VAR_B,Len);
   if (!Cfg.Macintosh) BackSpace(5);

/* z=b^2 */
   if (!Cfg.Macintosh) fprintf(stderr,"Part2");
   if (IsZero(VAR_B,Len/2)) {Done=1;ClearBigInt(Work3,Len);}
   else
     {
      SaveNum1FFT=90;
      FullMul(Work3,VAR_B,VAR_B,Len);
   /* Work1=Work1+6*b^2 */
      Add(Work1,Work1,Work3,Len);Add(Work1,Work1,Work3,Len);
      Add(Work1,Work1,Work3,Len);Add(Work1,Work1,Work3,Len);
      Add(Work1,Work1,Work3,Len);Add(Work1,Work1,Work3,Len);
   /* Work2=Work2+b^2 */
      Add(Work2,Work2,Work3,Len);
      if (!Cfg.Macintosh) BackSpace(5);

   /* z=b^3 */
      if (!Cfg.Macintosh) fprintf(stderr,"Part3");
      Num1IsCached=90;SaveNum2FFT=91;
      FullMul(VAR_B,VAR_B,Work3,Len);
      FlushFFTCache(91);
   /* Work1=Work1+4*b^3 */
      Add(Work1,Work1,VAR_B,Len);Add(Work1,Work1,VAR_B,Len);
      Add(Work1,Work1,VAR_B,Len);Add(Work1,Work1,VAR_B,Len);
   /* Work2=Work2+b^3 */
      Add(Work2,Work2,VAR_B,Len);
      if (!Cfg.Macintosh) BackSpace(5);

   /* b=b^4 */
      if (!Cfg.Macintosh) fprintf(stderr,"Part4");
      Num1IsCached=91;
      FullMul(Work3,Work3,Work3,Len);
   /* Work1=Work1+b^4 */
      Add(Work1,Work1,Work3,Len);
     }
   if (!Cfg.Macintosh) BackSpace(5);

   if (!Cfg.Macintosh) fprintf(stderr,"Part5");
/* Work1=Work1*a */
   FullMul(Work1,Work1,VAR_A,Len); /* N1L0 cut. */
/* Work2=Work2*2^(2n+3) */
   MulByFloat(Work2,Pow2,Len);Pow2*=4;
   if (!Cfg.Macintosh) BackSpace(5);

/* a=Work1-Work2 */
   Sub(VAR_A,Work1,Work2,Len);
   EndTime=clock();
   fprintf(stderr, "Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	   /CLOCKS_PER_SEC);
   DumpTimings(((double)EndTime-(double)LoopTime)
	   /CLOCKS_PER_SEC);

   if (TestKeyboard()) break;
   if (Cfg.AlwaysSaveData)
     SaveData(VAR_A,VAR_B,Work3,NO_NUM,NO_NUM,NO_NUM,
              ((double)EndTime-(double)LoopTime)
	   /CLOCKS_PER_SEC,Pow2,Pass,Len,2);
  }

LoopTime=clock();
if (Done)
  {
   fprintf(stderr,"Final    : ");
   BorReciprocal(Work1,VAR_A,Work3,Len);
   EndTime=clock();;
   fprintf(stderr, "Time= %0.2f\n", ((double)EndTime-(double)LoopTime)
	   /CLOCKS_PER_SEC);
   DumpTimings(((double)EndTime-(double)LoopTime)
	   /CLOCKS_PER_SEC);
   PrintFormattedPi("Borwein Quartic formula",((double)EndTime-(double)StartTime)
	   /CLOCKS_PER_SEC, Work1, Len);
   DeleteSaveFile();
  }
else SaveData(VAR_A,VAR_B,Work3,NO_NUM,NO_NUM,NO_NUM,
              ((double)clock()-(double)StartTime)
	   /CLOCKS_PER_SEC,Pow2,Pass,Len,2);

fprintf(stderr,"\nTotal Execution time: %0.2f seconds.\n\n",((double)clock()-(double)StartTime)
	   /CLOCKS_PER_SEC);
if (!Done) DumpTimings(((double)clock()-(double)StartTime)
	   /CLOCKS_PER_SEC);

return Done;
}



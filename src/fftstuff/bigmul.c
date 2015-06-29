#include "pi.h"
#include "fft.h"
#include "bigmul.h"
#include "cache.h"
#include "block.h"

static FFT_DATA_TYPE *FFTNum=NULL;
static BigInt FMWork=NO_NUM; /* for FractalMul() */
static size_t FFTLimit=0;
static FFT_DATA_TYPE MaxFFTError = 0.0;

#ifdef USE_FAST_X86_FLOOR
unsigned short int SetFPUMode(unsigned short int Mode,unsigned short int Mask)
{static short unsigned int CtrlWord;

asm volatile ("
        fstcw   (%0)
        fwait
       "
        : :"m" (CtrlWord) : "cc" , "memory");

CtrlWord= (CtrlWord & ~ Mask) | (Mode & Mask);
CtrlWord=0xffff;

asm volatile ("
        fldcw   (%0)
        fwait
       "
        : :"m" (CtrlWord) : "cc", "memory");

return CtrlWord;
}
#endif

size_t
CheckFFTMemReq(size_t Len)
/*
** Return how many bytes we would _like_ to have, not how many we
** need as a minimum, etc.
*/
{
if (Cfg.HalfMuls) Len/=2; /* This formula will only use half sized muls */
if (Len > MaxFFTLen/RawIntDigits) Len=MaxFFTLen/RawIntDigits;
return CalcFFTLen(Len)*sizeof(FFT_DATA_TYPE);
}

void
InitFFTMul(size_t Len)
{
if (Cfg.HalfMuls) Len/=2; /* This formula will only use half sized muls */
FFTLimit=Len;
if (FFTLimit > MaxFFTLen/RawIntDigits) FFTLimit=MaxFFTLen/RawIntDigits;

while (CalcFFTLen(FFTLimit)*sizeof(FFT_DATA_TYPE) > CoreMemAvail)
   FFTLimit/=2;

FFTNum=(FFT_DATA_TYPE*)CoreMemPtr;

fprintf(stderr,"FFT is initialized to %lu len (%lu digits.)\n",
        (ULINT)FFTLimit,(ULINT)FFTLimit*RawIntDigits);

if (Cfg.HalfMuls)
  fprintf(stderr,"(This formula will only use 'half' sized muls.)\n");

if (FFTLimit < Len)
  {
   if (Cfg.AllowFractalMul)
     {
      FMWork=CreateBigInt(Len*4);
      fprintf(stderr,"NOTICE: There will be %u levels of FractalMul().\n",
              Log2(Len/FFTLimit));
     }
   else FatalError("Insufficient memory to multiply without AllowFractalMul.\n");
  }
}

void
DeInitFFTMul(void)
{
}

#ifdef USE_DEFAULT_CONVOLUTION
/*
** Do a convolution.  This is the standard 'complex' output
** format FFT.  (Numerical Recipes also uses this, although its
** not original or exclusive to them.)
*/
void
DoConvolution(FFT_DATA_TYPE *FFTNum,int Cache,size_t Len2)
{size_t x;
StartTimer(ConvTime);
if (Cache==-2) /* In memory convolution of both nums */
  {FFT_DATA_TYPE *FFTNum2=FFTNum+Len2;
   DumpDebug("Cm...");
   FFTNum[0] = FFTNum[0] * FFTNum2[0];
   FFTNum[1] = FFTNum[1] * FFTNum2[1];
   for (x = 2; x < Len2; x += 2)
     {
       double a, b, c, d;
       FFTNum[x] = (a = FFTNum[x])     * (b = FFTNum2[x]) -
                   (c = FFTNum[x + 1]) * (d = FFTNum2[x + 1]);
       FFTNum[x + 1] = a * d + c * b;
     }
  }
else if (Cache==-1)  /* In memory self convolution. */
  {
   DumpDebug("C");
   FFTNum[0] = FFTNum[0] * FFTNum[0];
   FFTNum[1] = FFTNum[1] * FFTNum[1];
   for (x = 2; x < Len2; x += 2)
     {
       double a, b, c, d;
       FFTNum[x] = (a = FFTNum[x])     * (b = FFTNum[x]) -
                   (c = FFTNum[x + 1]) * (d = FFTNum[x + 1]);
       FFTNum[x + 1] = a * d + c * b;
     }
  }
#ifdef VIRTUAL_CACHE
else /* Virtual Mem based convolution. */
  {FFT_DATA_TYPE *Num1=FFTNum;
   FFT_DATA_TYPE *Num2=FFTCash[Cache].Mem;

   if (Num2==NULL)
     FatalError("Cache %d doesn't exist.\n",Cache);

   DumpDebug("C%d...",Cache);
   Num1[0] = Num1[0] * Num2[0];
   Num1[1] = Num1[1] * Num2[1];
   for (x = 2; x < Len2; x += 2)
     {
       double a, b, c, d;
       Num1[x] = (a = Num1[x])     * (b = Num2[x]) -
                 (c = Num1[x + 1]) * (d = Num2[x + 1]);
       Num1[x + 1] = a * d + c * b;
     }
  }
#else
else
  {FILE *f;
   size_t L2;
   FFT_DATA_TYPE *Buf1=FFTNum;
   FFT_DATA_TYPE *Buf2=(FFT_DATA_TYPE*)FixedBuf;

   DumpDebug("C%d...",Cache);
   f=fopen(FFTCash[Cache].Name,"rb");
   if (f==NULL)
     FatalError("Unable to open '%s' for convolution.\n",FFTCash[Cache].Name);
   fread(Buf2,sizeof(FFT_DATA_TYPE),2,f);
   Buf1[0] = Buf1[0] * Buf2[0];
   Buf1[1] = Buf1[1] * Buf2[1];
   L2=Len2-2;
   Buf1+=2;
   while (L2)
     {size_t Sz;
      Sz=Min(FIXEDBUF_SIZE/sizeof(FFT_DATA_TYPE),L2);
      StartTimer(DiskIOTime);
      fread(Buf2,sizeof(FFT_DATA_TYPE),Sz,f);
      StopTimer(DiskIOTime);
      for (x = 0; x < Sz; x += 2)
        {
          double a, b, c, d;
          Buf1[x] = (a = Buf1[x])     * (b = Buf2[x]) -
                    (c = Buf1[x + 1]) * (d = Buf2[x + 1]);
          Buf1[x + 1] = a * d + c * b;
        }
      Buf1+=Sz;
      L2-=Sz;
     }
   if (ferror(f)) FatalError("Error convoluting '%s'\n",FFTCash[Cache].Name);
   fclose(f);
  }
#endif
StopTimer(ConvTime);
}

#endif

/*
** This is based on the simple Fractal / Divide and Conquer / Karatsuba
** O(n^1.585) method.
**
** It's fairly simple.  a*b is: a1b1(B^2+B)+(a1-a2)(b2-b1)B+a2b2(B+1)
**
** You need 4*SIZE storage for the product and the working space.
*/
static void
FractalMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len)
{size_t HLen;
 int Sign1, Sign2;
 BigInt offset;
 BigInt Work=Prod+2*Len;
 int OldNum1IsCached, OldNum2IsCached;
 int OldSaveNum1FFT, OldSaveNum2FFT;

  if (Len <= FFTLimit) {FFTMul(Prod,Num1,Num2,Len,Len*2,0);return;}

  HLen = Len/2;

  OldNum1IsCached=Num1IsCached;OldNum2IsCached=Num2IsCached;
  OldSaveNum1FFT=SaveNum1FFT;OldSaveNum2FFT=SaveNum2FFT;
  Num1IsCached=Num2IsCached=SaveNum1FFT=SaveNum2FFT=0;

  if (Num1==Num2)
    {
     Sign1 = Sub(Prod, Num1, Num1+HLen, HLen);
     if (Sign1) Negate(Prod, HLen);
     FractalMul(Work,Prod,Prod,HLen);
     ClearBigInt(Prod,Len*2);
     offset = Prod + HLen;
     RippleSub(Prod,offset,offset,Work,Len); /* square makes sign1 pos */
    }
  else
    {
     /* Do x=Right(Num1)-Left(Num1) y=Left(Num2)-Right(Num2) */
     Sign1 = Sub(Prod, Num1+HLen, Num1, HLen);
     if (Sign1) Negate(Prod, HLen);

     Sign2 = Sub(Prod+HLen, Num2, Num2+HLen, HLen);
     if (Sign2) Negate(Prod+HLen, HLen);

     FractalMul(Work, Prod, Prod+HLen, HLen);

     ClearBigInt(Prod,Len*2);
     offset = Prod + HLen;
     if (Sign1 == Sign2)  RippleAdd(Prod,offset,offset,Work,Len);
     else                 RippleSub(Prod,offset,offset,Work,Len);
    }

#if 1
/* Turn the FFT/NTT caching back on. */
  Num1IsCached=OldNum1IsCached;Num2IsCached=OldNum2IsCached;
  SaveNum1FFT=OldSaveNum1FFT;SaveNum2FFT=OldSaveNum2FFT;
#endif

  FractalMul(Work, Num1, Num2, HLen);
  offset = Prod + HLen;RippleAdd(Prod,offset,offset,Work,Len);
  Add(Prod, Prod, Work, Len);

  FractalMul(Work, Num1 + HLen, Num2 + HLen, HLen);
  offset = Prod + HLen;RippleAdd(Prod,offset,offset,Work,Len);
  offset = Prod + Len; RippleAdd(Prod,offset,offset,Work,Len);

  Num1IsCached=OldNum1IsCached;Num2IsCached=OldNum2IsCached;
  SaveNum1FFT=OldSaveNum1FFT;SaveNum2FFT=OldSaveNum2FFT;
}

volatile float MAGIC1=3.0*2147483648.0*2147483648.0;
volatile float MAGIC2=3.0*2147483648.0*2147483648.0;
volatile float JUSTIFY=6755399441055744.0;

/*
** Do a Fast Fourier Transform based multiplication.
*/
INT32
FFTMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t NumLen,
       size_t ProdLen, INT32 Scale)
/*
** Scale=0 means don't multiply by Scale.
** Scale!=0 means do multiply by Scale.  (Usually 10.)
**
** If a scaling causes a carry, that carry will be returned, else
** zero will be returned.
*/
{
  size_t x;
  size_t NumLen2 = NumLen * 2;
  size_t FFTLen2=CalcFFTLen(NumLen);
  FFT_DATA_TYPE *FFTNum2=FFTNum;
  int All16=0;INT32 ScaleCarry=0;

  if (NumLen <= 64)
    {INT32 *Buf1=(INT32*)CoreMemPtr;
     INT32 *Buf2=Buf1+NumLen;
     INT32 *DBuf=Buf2+NumLen;
     ReadNumIntoBuf(Num1,Buf1,NumLen);
     ReadNumIntoBuf(Num2,Buf2,NumLen);
     BlockClear(DBuf,DBuf+NumLen*2);
     BlockSlowMul(DBuf,Buf1,Buf2,NumLen);
     if (Scale) ScaleCarry=BlockMulBy(DBuf,DBuf,Scale,0,NumLen*2);
     WriteBufIntoNum(DBuf,Prod,ProdLen);
     return ScaleCarry;
    }

  if (NumLen > FFTLimit)
    {
     if (Cfg.AllowFractalMul)
       {
        FractalMul(FMWork,Num1,Num2,NumLen);
        if (Scale) ScaleCarry=MulBy(Prod,FMWork,Scale,ProdLen);
        else       Copy(Prod,FMWork,ProdLen);
        return ScaleCarry;
       }
     else
       FatalError("Somehow BigMul was called with a length (%lu) longer than FFTLimit (%lu)\n",
                  (ULINT)NumLen,(ULINT)FFTLimit);
    }

  MaxFFTError=0.0;
  if (FFTLen2*2*sizeof(FFT_DATA_TYPE) <= CoreMemAvail)
    {All16=-2;FFTNum2=FFTNum+FFTLen2;}

  if (!IsPow2(NumLen))
    FatalError("The FFT size is not a power of two\n");

  if (NumLen > MaxFFTLen/RawIntDigits)
    FatalError("Somehow FFTMul was called with a number longer than MAX_FFT_LIMIT.\n%lu %lu\n",
               (UINT32)NumLen,(UINT32)MaxFFTLen);

  if (NumLen > FFTLimit)
       FatalError("Somehow BigMul was called with a length (%lu) longer than FFTLimit (%lu)\n",
                  (ULINT)NumLen,(ULINT)FFTLimit);

  DumpDebug("FFT %s ",Num2Str(NumLen*RawIntDigits));
  StartTimer(FFTMulTime);
  FFTDisk-=DiskIOTime;
  if (Num1 == Num2)
    {int Line=0;
     if (Num1IsCached || Num2IsCached)
       Line=CheckFFTCache(Num2,NumLen,Num1IsCached + Num2IsCached,0);
     if (Line) LoadFFTFromCache(Line,FFTNum);
     else
       {
        DumpDebug("F");
        FwdTransform(FFTNum,Num2,NumLen);
        if (SaveNum2FFT || SaveNum1FFT)
          SaveFFTIntoCache(FFTNum, FFTLen2, Num2, NumLen,
                           SaveNum1FFT + SaveNum2FFT,0);
       }
     DoConvolution(FFTNum,-1,FFTLen2);
    }
  else
    {int Line1=0,Line2=0;
     if (Num1IsCached) Line1=CheckFFTCache(Num1,NumLen,Num1IsCached,0);
     if (Num2IsCached) Line2=CheckFFTCache(Num2,NumLen,Num2IsCached,0);
     if (Line1 && Line2)
       {
        DumpDebug("Caches...");
        LoadFFTFromCache(Line1,FFTNum);
        DoConvolution(FFTNum,Line2,FFTLen2);
       }
     else if (Line1)
       {
        DumpDebug("Cache...F");
        FwdTransform(FFTNum,Num2,NumLen);
        if (SaveNum2FFT)
          SaveFFTIntoCache(FFTNum,FFTLen2,Num2,NumLen, SaveNum2FFT,0);
        DoConvolution(FFTNum,Line1,FFTLen2);
       }
     else if (Line2)
       {
        DumpDebug("Cache...F");
        FwdTransform(FFTNum,Num1,NumLen);
        if (SaveNum1FFT)
          SaveFFTIntoCache(FFTNum,FFTLen2,Num1,NumLen, SaveNum1FFT,0);
        DoConvolution(FFTNum,Line2,FFTLen2);
       }
     else
       {
        if (SaveNum1FFT)
          {
           DumpDebug("F");
           FwdTransform(FFTNum,Num1,NumLen);
           Line1=SaveFFTIntoCache(FFTNum,FFTLen2,Num1,NumLen, SaveNum1FFT,0);
           if (All16) Line1=-2;
           else if (Line1==0) SaveFFTIntoCache0(FFTNum,FFTLen2);
           DumpDebug("F");
           FwdTransform(FFTNum2,Num2,NumLen);
           if (SaveNum2FFT)
             SaveFFTIntoCache(FFTNum2,FFTLen2,Num2,NumLen, SaveNum2FFT,0);
           DoConvolution(FFTNum,Line1,FFTLen2);
          }
        else if (SaveNum2FFT)
          {
           DumpDebug("F");
           FwdTransform(FFTNum,Num2,NumLen);
           Line2=SaveFFTIntoCache(FFTNum,FFTLen2,Num2,NumLen, SaveNum2FFT,0);
           if (All16) Line2=-2;
           else if (Line2==0) SaveFFTIntoCache0(FFTNum,FFTLen2);
           DumpDebug("F");
           FwdTransform(FFTNum2,Num1,NumLen);
           if (SaveNum1FFT)
             SaveFFTIntoCache(FFTNum2,FFTLen2,Num1,NumLen, SaveNum1FFT,0);
           DoConvolution(FFTNum,Line2,FFTLen2);
          }
        else
          {
           Line2=All16;
           DumpDebug("F");
           FwdTransform(FFTNum,Num1,NumLen);
           if (!All16) SaveFFTIntoCache0(FFTNum,FFTLen2);
           DumpDebug("F");
           FwdTransform(FFTNum2,Num2,NumLen);
           DoConvolution(FFTNum,Line2,FFTLen2);
          }
       }
    }

  DeleteFFTCache(0); /* get rid of the convolution file */
/*
** Now do an Inverse FFT
*/
  DumpDebug("R");
  RevTransform(FFTNum,NumLen);

  DumpDebug("Carries...");
  StartTimer(CarryTime);
  {FFT_DATA_TYPE Carry, RawPyramid, Pyramid;
//   FFT_DATA_TYPE PyramidError;
   FFT_INT *P=(FFT_INT*)FFTNum;
#ifndef FFT_IS_NORMALIZED
   FFT_DATA_TYPE inv;

#if FFT_DIGITS==4
   inv = 1.0 / NumLen2;
#else
   inv = 2.0 / FFTLen2;
#endif
#endif

#ifdef USE_FAST_X86_FLOOR

#define FPU_R_MODE 0x0c00
#define FPU_R_NEAR 0x0000
#define FPU_R_DOWN 0x0400
#define FPU_R_UP   0x0800
#define FPU_R_CHOP 0x0c00

#if FFT_DIGITS!=4
#error hardwired for 4 digits.
#endif
{
/* about 15% faster total execution for 1m digits of pi. */
int OldControl=SetFPUMode(0,0);
volatile union {double D;long int i[2];} ChopShop;

SetFPUMode(FPU_R_DOWN,FPU_R_MODE);
   Carry = 0.0;
   for (x = 0; x < FFTLen2; x++)
      {
#ifdef FFT_IS_NORMALIZED
        RawPyramid = FFTNum[x] + Carry;
#else
        RawPyramid = FFTNum[x] * inv + Carry;
#endif
        Pyramid = (((RawPyramid + 0.5) + MAGIC1) - MAGIC2);
        Carry=((Pyramid*0.0001000000000000001)+MAGIC1)-MAGIC2;
        ChopShop.D=(Pyramid - Carry * 10000.0)+JUSTIFY;
        P[x]=ChopShop.i[0];
      }
//_control87(OldControl,MCW_RC);
SetFPUMode(OldControl,FPU_R_MODE);
}

#else /* slow generic release carries */

   Carry = 0.0;
   for (x = 0; x < FFTLen2; x++)
      {
#ifdef FFT_IS_NORMALIZED
        RawPyramid = FFTNum[x] + Carry;
#else
        RawPyramid = FFTNum[x] * inv + Carry;
#endif
        Pyramid = FLOOR(RawPyramid + 0.5);
//        PyramidError = RawPyramid - Pyramid;if (PyramidError <= 0.0) PyramidError = -PyramidError;
//        if (PyramidError > MaxFFTError) MaxFFTError = PyramidError;
#if FFT_DIGITS==4
        Carry = FLOOR(Pyramid / 10000.0);
        P[x]= (Pyramid - Carry * 10000.0);
#else
        Carry = FLOOR(Pyramid / 100.0);
        P[x]= (Pyramid - Carry * 100.0);
#endif
//        Carry = FLOOR(Pyramid / FBase);
//        P[x]= (Pyramid - Carry * FBase);
      }
#endif



   StopTimer(CarryTime);
   for (x=0;x<FFTLen2/2;x++) {FFT_INT t=P[x];P[x]=P[FFTLen2-x-1];P[FFTLen2-x-1]=t;}
#if FFT_DIGITS==4
   BlockPack4((INT32*)P,NumLen2);
#else
   BlockPack2((INT32*)P,NumLen2);
#endif
   if (Scale) ScaleCarry=BlockMulBy((INT32*)P,(INT32*)P,Scale,0,Min(ProdLen,NumLen2));
   WriteBufIntoNum((INT32*)P,Prod,Min(ProdLen,NumLen2));

   if (Carry != 0.0)
     FatalError("Oops!  FFT Carry=%f\n",(double)Carry);
  }

  /*
  ** Do a bit of 'sanity' error checking.  This value
  ** is based on some testing with a 'test jig' and
  ** personal opinion.  Should be good enough to catch
  ** systems with poor FPUs.  However, if the value ever
  ** actually reaches 0.5, then you know for a FACT that
  ** the FFTMul() has failed.
  */
  if (MaxFFTError >= 0.2)
    {
      printf("\n**WARNING** Len=%lu Max FFT Error: %f\n",
             (ULINT)NumLen,(double)MaxFFTError);
      puts("Either the FFT is approaching its limit, or a sofware");
      puts("or hardware error occured.  This can also be caused by");
      puts("a program bug, poor trig, etc.  You may wish to lower");
      puts("the MaxFFTLen limit in fft.h and run it again.");
      ExitPrg(EXIT_FAILURE);
    }

  StopTimer(FFTMulTime);
  FFTDisk+=DiskIOTime;
  DumpDebug("Done FFT.\n");
  return ScaleCarry;
}



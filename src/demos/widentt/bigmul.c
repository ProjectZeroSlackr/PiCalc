#include "pi.h"

#include "block.h"
#include "ntt.h"
#include "bigmul.h"
#include "cache.h"

#include "modmath.h"
#include "primes.h"

static ModInt *FFTNum=NULL;
static BigInt FMWork=NO_NUM; /* for FractalMul() */
static size_t FFTLimit;

size_t
CheckFFTMemReq(size_t Len)
/*
** Return how many bytes we would _like_ to have, not how many we
** need as a minimum, etc.
*/
{
if (Cfg.HalfMuls) Len/=2; /* This formula will only use half sized muls */
if (Len > MaxFFTLen/RawIntDigits) Len=MaxFFTLen/RawIntDigits;
/* NumLen->NTTLen.  Enough to hold both nums at once.*/
return CalcNTTLen(Len)*2*sizeof(ModInt);
}

void
InitFFTMul(size_t Len)
{
  if (Cfg.HalfMuls) Len/=2; /* This formula will only use half sized muls */
  FFTLimit=Len;
  if (FFTLimit > MaxFFTLen/RawIntDigits) FFTLimit=MaxFFTLen/RawIntDigits;

  while (CalcNTTLen(FFTLimit)*sizeof(ModInt) > CoreMemAvail)
    FFTLimit/=2;

  FFTNum=(ModInt*)CoreMemPtr;

  fprintf(stderr,"FFT is initialized to %lu len (%lu digits.)\n",
          (ULINT)FFTLimit,(ULINT)FFTLimit*RawIntDigits);
  if (Cfg.HalfMuls)
    fprintf(stderr,"(This formula will only use 'half' sized muls.)\n");

  if (FFTLimit < Len)
    {
     if (Cfg.AllowFractalMul)
       {
        FMWork=CreateBigInt(Len*4);
        fprintf(stderr,"NOTICE: There will be %d levels of FractalMul().\n",
                Log2(Len/FFTLimit));
       }
     else FatalError("Insufficient memory to multiply without AllowFractalMul.\n");
    }
}

void
DeInitFFTMul(void)
{
}

static void
DoConvolutions(int Cache,size_t NTTLen)
{size_t x;

StartTimer(ConvTime);
if (Cache==-2) /* In memory convolution of both nums */
  {ModInt *FFTNum2=FFTNum+NTTLen;
   DumpDebug("Cm...");
   for (x=0;x<NTTLen;x++) ModMul(&FFTNum[x],&FFTNum[x],&FFTNum2[x]);
  }
else if (Cache==-1)  /* In memory self convolution. */
  {
   DumpDebug("C");
   for (x=0;x<NTTLen;x++) ModMul(&FFTNum[x],&FFTNum[x],&FFTNum[x]);
  }
#ifdef VIRTUAL_CACHE
else /* Virtual Mem based convolution. */
  {ModInt *Num1=FFTNum;
   ModInt *Num2=(ModInt*)FFTCash[Cache].Mem;

   if (Num2==NULL)
     FatalError("Cache %d doesn't exist.\n",Cache);

   DumpDebug("C%d...",Cache);
   for (x=0;x<NTTLen;x++) ModMul(&Num1[x],&Num1[x],&Num2[x]);
  }
#else
else /* Disk based convolution. */
  {FILE *f;ModInt *Num1=FFTNum;
   ModInt *Num2=(ModInt*)FixedBuf;
   size_t L2=NTTLen;

   DumpDebug("C%d...",Cache);
   StartTimer(DiskIOTime);
   f=fopen(FFTCash[Cache].Name,"rb");
   StopTimer(DiskIOTime);
   if (f==NULL)
     FatalError("Unable to open '%s' for convolution.\n",FFTCash[Cache].Name);

   while (L2)
     {size_t Sz;
      Sz=Min(FIXEDBUF_SIZE/sizeof(ModInt),L2);
      StartTimer(DiskIOTime);
      fread(Num2,sizeof(ModInt),Sz,f);
      StopTimer(DiskIOTime);
      for (x=0;x<Sz;x++) ModMul(&Num1[x],&Num1[x],&Num2[x]);
      Num1+=Sz;
      L2-=Sz;
     }

   if (ferror(f)) FatalError("Error convoluting '%s'\n",FFTCash[Cache].Name);
   fclose(f);
  }
#endif

StopTimer(ConvTime);
}

/*
** Do a Number Theoretic Transform based multiplication.
*/
static INT32
DoNTTMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t NumLen,
         size_t ProdLen, INT32 Scale)
{
 size_t x, NumLen2 = NumLen * 2;
 size_t NTTLen=CalcNTTLen(NumLen);
 ModInt *FFTNum2=FFTNum;
 int BothInMem=0;
 INT32 ScaleCarry=0;

/* Can we do both numbers in memory? */
  if (NTTLen*sizeof(ModInt)*2 <= CoreMemAvail)
    {BothInMem=-2;FFTNum2=FFTNum+NTTLen;}

  if (Num1 == Num2)
    {int Line=0;
     if (Num1IsCached || Num2IsCached)
       Line=CheckFFTCache(Num2,NumLen,Num1IsCached + Num2IsCached,0);
     if (Line) LoadFFTFromCache(Line,FFTNum);
     else
       {
        DumpDebug("F");
        DoFwdTransforms(FFTNum,Num1,NumLen);
        if (SaveNum2FFT || SaveNum1FFT)
          SaveFFTIntoCache(FFTNum,NTTLen,Num2,NumLen,
                           SaveNum1FFT+SaveNum2FFT,0);
       }
     DoConvolutions(-1,NTTLen);
    }
  else
    {int Line1=0,Line2=0;
     if (Num1IsCached) Line1=CheckFFTCache(Num1,NumLen,Num1IsCached,0);
     if (Num2IsCached) Line2=CheckFFTCache(Num2,NumLen,Num2IsCached,0);
     if (Line1 && Line2)
       {
        DumpDebug("Caches...");
        LoadFFTFromCache(Line1,FFTNum);
        DoConvolutions(Line2,NTTLen);
       }
     else if (Line1)
       {
        DumpDebug("Cache...F");
        DoFwdTransforms(FFTNum,Num2,NumLen);
        if (SaveNum2FFT)
          SaveFFTIntoCache(FFTNum,NTTLen,Num2,NumLen,SaveNum2FFT,0);
        DoConvolutions(Line1,NTTLen);
       }
     else if (Line2)
       {
        DumpDebug("Cache...F");
        DoFwdTransforms(FFTNum,Num1,NumLen);
        if (SaveNum1FFT)
          SaveFFTIntoCache(FFTNum,NTTLen,Num1,NumLen,SaveNum1FFT,0);
        DoConvolutions(Line2,NTTLen);
       }
     else
       {
        if (SaveNum1FFT)
          {
           DumpDebug("F");
           DoFwdTransforms(FFTNum,Num1,NumLen);
           Line1=SaveFFTIntoCache(FFTNum,NTTLen,Num1,NumLen,SaveNum1FFT,0);
           if (BothInMem) Line1=-2;
           else if (Line1==0) SaveFFTIntoCache0(FFTNum,NTTLen);
           DumpDebug("F");
           DoFwdTransforms(FFTNum2,Num2,NumLen);
           if (SaveNum2FFT)
             SaveFFTIntoCache(FFTNum2,NTTLen,Num2,NumLen,SaveNum2FFT,0);
           DoConvolutions(Line1,NTTLen);
          }
        else if (SaveNum2FFT)
          {
           DumpDebug("F");
           DoFwdTransforms(FFTNum,Num2,NumLen);
           Line2=SaveFFTIntoCache(FFTNum,NTTLen,Num2,NumLen,SaveNum2FFT,0);
           if (BothInMem) Line2=-2;
           else if (Line2==0) SaveFFTIntoCache0(FFTNum,NTTLen);
           DumpDebug("F");
           DoFwdTransforms(FFTNum2,Num1,NumLen);
           if (SaveNum1FFT)
             SaveFFTIntoCache(FFTNum2,NTTLen,Num1,NumLen,SaveNum1FFT,0);
           DoConvolutions(Line2,NTTLen);
          }
        else
          {
           Line2=BothInMem;
           DumpDebug("F");
           DoFwdTransforms(FFTNum,Num1,NumLen);
           if (!BothInMem) SaveFFTIntoCache0(FFTNum,NTTLen);
           DumpDebug("F");
           DoFwdTransforms(FFTNum2,Num2,NumLen);
           DoConvolutions(Line2,NTTLen);
          }
       }
    }

  DumpDebug("R");
  DoRevTransforms(FFTNum,NumLen);

  DumpDebug("Carries...");
  StartTimer(CarryTime);
  {
   INT32 *ProdBuf=(INT32*)FixedBuf;
   size_t ProdNdx,ProdPos,ProdBufLen=FIXEDBUF_SIZE/sizeof(INT32);
   ModInt Carry, Pyramid;

   ProdNdx=ProdBufLen;ProdPos=NumLen2;
   ModClear(&Carry);
   for (x = 0; x < NTTLen; x++)
      {int q;
       ModMul(&Pyramid,&FFTNum[x],&MulInv);
       CB_Add(&Carry,&Pyramid,&Carry);
       ProdNdx-=RawModInts;ProdPos-=RawModInts;
       for (q=RawModInts-1;q>=0;q--)
         CB_DivInt(&Carry,100000000,&Carry,&ProdBuf[ProdNdx+q]);
       if (ProdNdx==0)
          {size_t z;
           if (ProdPos >= ProdLen) z=0;
           else if (ProdPos+ProdBufLen < ProdLen) z=ProdBufLen;
           else z=ProdLen-ProdPos;
           if ((Scale) && (z))
             ScaleCarry=BlockMulBy(ProdBuf,ProdBuf,Scale,ScaleCarry,z);
           if (z) WriteBufIntoNum(ProdBuf,Prod+ProdPos,z);
           ProdNdx=ProdBufLen;
          }
      }
   if (ProdNdx!=ProdBufLen)
     {size_t z;
      ProdBufLen=ProdBufLen-ProdNdx;
      if (ProdPos >= ProdLen) z=0;
      else if (ProdPos+ProdBufLen < ProdLen) z=ProdBufLen;
      else z=ProdLen-ProdPos;
      if ((Scale) && (z))
         ScaleCarry=BlockMulBy(ProdBuf+ProdNdx,ProdBuf+ProdNdx,Scale,ScaleCarry,z);
      if (z) WriteBufIntoNum(ProdBuf+ProdNdx,Prod+ProdPos,z);
     }
   StopTimer(CarryTime);
  }
return ScaleCarry;
}

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

INT32
FFTMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t NumLen,
       size_t ProdLen, INT32 Scale)
/*
** Scale=0 means don't multiply by Scale.
** Scale!=0 means do multiply by Scale. (Usually 10.)
**
** If a scaling causes a carry, that carry will be returned, else
** zero will be returned.
*/
{INT32 Carry=0;

  if (NumLen <= 64)
    {INT32 *Buf1=(INT32*)CoreMemPtr;
     INT32 *Buf2=Buf1+NumLen;
     INT32 *DBuf=Buf2+NumLen;
     ReadNumIntoBuf(Num1,Buf1,NumLen);
     ReadNumIntoBuf(Num2,Buf2,NumLen);
     BlockClear(DBuf,DBuf+NumLen*2);
     BlockSlowMul(DBuf,Buf1,Buf2,NumLen);
     if (Scale) Carry=BlockMulBy(DBuf,DBuf,Scale,0,NumLen*2);
     WriteBufIntoNum(DBuf,Prod,ProdLen);
     return Carry;
    }

  if (NumLen > FFTLimit)
    {
     if (Cfg.AllowFractalMul)
       {
        FractalMul(FMWork,Num1,Num2,NumLen);
        if (Scale) Carry=MulBy(Prod,FMWork,Scale,ProdLen);
        else       Copy(Prod,FMWork,ProdLen);
        return Carry;
       }
     else
       FatalError("Somehow BigMul was called with a length (%lu) longer than FFTLimit (%lu)\n",
                  (ULINT)NumLen,(ULINT)FFTLimit);
    }

/*
** This is a radix-2 FFT, so the length has to be a power of two.
*/
  if (!IsPow2(NumLen))
    FatalError("The FFT size is not a power of two\n");

  DumpDebug("FFT %s ",Num2Str(NumLen*RawIntDigits));
  StartTimer(FFTMulTime);
  FFTDisk-=DiskIOTime;

  SetNTTSize(NumLen);
  Carry=DoNTTMul(Prod,Num1,Num2,NumLen,ProdLen,Scale);

  StopTimer(FFTMulTime);
  FFTDisk+=DiskIOTime;
  DumpDebug("Done fft.\n");

  return Carry;
}



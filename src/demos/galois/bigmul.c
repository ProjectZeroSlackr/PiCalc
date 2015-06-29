#include "pi.h"

#include "block.h"
#include "fgt.h"
#include "bigmul.h"
#include "cache.h"

#include "modmath.h"
#include "primes.h"

static FFT_DATA_TYPE *FFTNum=NULL;
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
/* NumLen->FGTLen.  Enough to hold both nums at once.*/
return CalcFGTLen(Len)*2*sizeof(ModInt);
}

void
InitFFTMul(size_t Len)
{
  if (Cfg.HalfMuls) Len/=2; /* This formula will only use half sized muls */
  FFTLimit=Len;
  if (FFTLimit > MaxFFTLen/RawIntDigits) FFTLimit=MaxFFTLen/RawIntDigits;

  while (CalcFGTLen(FFTLimit)*sizeof(ModInt) > CoreMemAvail)
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
DoFGTConvolution(FFT_DATA_TYPE *FFTNum1,FFT_DATA_TYPE *FFTNum2, size_t Len2)
/*
** The convolution requires two seperate points to do a single
** point of the convolution.  This means that while doing the
** first half, it overwrites the data it'd need for the second
** half.  So we have to do both points at once.
*/
{FFT_DATA_TYPE P1,P2;
 CmplxModInt xp,xm,yp,ym,t1,t2,t3,t4,t5;
 int x;

 P1=NthRoot;
 P2=NthRoot1;

 xp=FFTNum1[0];xm=FFTNum1[0];CmplxConj(&xm);
 yp=FFTNum2[0];ym=FFTNum2[0];CmplxConj(&ym);
 CmplxModAdd(&t4,&xp,&xm);CmplxModAdd(&t5,&yp,&ym);CmplxModMul(&t1,&t4,&t5);
 CmplxModMul(&t4,&xp,&yp);CmplxModMul(&t5,&xm,&ym);CmplxModSub(&t2,&t4,&t5);CmplxModAdd(&t2,&t2,&t2);
 CmplxModSub(&t4,&xp,&xm);CmplxModSub(&t5,&yp,&ym);CmplxModMul(&t3,&t4,&t5);
 CmplxModAdd(&t4,&t1,&t2);CmplxModSub(&FFTNum1[0],&t4,&t3);

 for (x=1;x< Len2/2;x++)
   {
    xp=FFTNum1[x];xm=FFTNum1[Len2-x];CmplxConj(&xm);
    yp=FFTNum2[x];ym=FFTNum2[Len2-x];CmplxConj(&ym);
    CmplxModAdd(&t4,&xp,&xm);CmplxModAdd(&t5,&yp,&ym);CmplxModMul(&t1,&t4,&t5);
    CmplxModMul(&t4,&xp,&yp);CmplxModMul(&t5,&xm,&ym);CmplxModSub(&t2,&t4,&t5);CmplxModAdd(&t2,&t2,&t2);
    CmplxModSub(&t4,&xp,&xm);CmplxModSub(&t5,&yp,&ym);CmplxModMul(&t3,&t4,&t5);
    CmplxModMul(&t3,&t3,&P1);CmplxModMul(&P1,&P1,&NthRoot);
    CmplxModAdd(&t4,&t1,&t2);CmplxModSub(&FFTNum1[x],&t4,&t3);

    t1=xm;CmplxConj(&t1);t2=xp;CmplxConj(&t2);xm=t2;xp=t1;
    t1=ym;CmplxConj(&t1);t2=yp;CmplxConj(&t2);ym=t2;yp=t1;
    CmplxModAdd(&t4,&xp,&xm);CmplxModAdd(&t5,&yp,&ym);CmplxModMul(&t1,&t4,&t5);
    CmplxModMul(&t4,&xp,&yp);CmplxModMul(&t5,&xm,&ym);CmplxModSub(&t2,&t4,&t5);CmplxModAdd(&t2,&t2,&t2);
    CmplxModSub(&t4,&xp,&xm);CmplxModSub(&t5,&yp,&ym);CmplxModMul(&t3,&t4,&t5);
    CmplxModMul(&t3,&t3,&P2);CmplxModMul(&P2,&P2,&NthRoot1);
    CmplxModAdd(&t4,&t1,&t2);CmplxModSub(&FFTNum1[Len2-x],&t4,&t3);
   }
 xp=FFTNum1[Len2/2];xm=FFTNum1[Len2/2];CmplxConj(&xm);
 yp=FFTNum2[Len2/2];ym=FFTNum2[Len2/2];CmplxConj(&ym);
 CmplxModAdd(&t4,&xp,&xm);CmplxModAdd(&t5,&yp,&ym);CmplxModMul(&t1,&t4,&t5);
 CmplxModMul(&t4,&xp,&yp);CmplxModMul(&t5,&xm,&ym);CmplxModSub(&t2,&t4,&t5);CmplxModAdd(&t2,&t2,&t2);
 CmplxModSub(&t4,&xp,&xm);CmplxModSub(&t5,&yp,&ym);CmplxModMul(&t3,&t4,&t5);
 CmplxModMul(&t3,&t3,&P1);//CmplxModMul(&P1,&P1,&NthRoot);
 CmplxModAdd(&t4,&t1,&t2);CmplxModSub(&FFTNum1[Len2/2],&t4,&t3);
}

static void
DoConvolutions(int Cache,size_t FGTLen)
{

StartTimer(ConvTime);
if (Cache==-2) /* In memory convolution of both nums */
  {FFT_DATA_TYPE *FFTNum2=FFTNum+FGTLen;
   DumpDebug("Cm...");
   DoFGTConvolution(FFTNum,FFTNum2,FGTLen);
  }
else if (Cache==-1)  /* In memory self convolution. */
  {
   DumpDebug("C");
   DoFGTConvolution(FFTNum,FFTNum,FGTLen);
  }
#ifdef VIRTUAL_CACHE
else /* Virtual Mem based convolution. */
  {FFT_DATA_TYPE *FFTNum2=(FFT_DATA_TYPE*)FFTCash[Cache].Mem;

   if (FFTNum2==NULL)
     FatalError("Cache %d doesn't exist.\n",Cache);

   DumpDebug("C%d...",Cache);
   DoFGTConvolution(FFTNum,FFTNum2,FGTLen);
  }
#else
#error Can not do a disk convolution.
else /* Disk based convolution. */
  {FILE *f;ModInt *Num1=FFTNum;
   FFT_DATA_TYPE *Num2=(FFT_DATA_TYPE*)FixedBuf;
   size_t L2=FGTLen;

   DumpDebug("C%d...",Cache);
   StartTimer(DiskIOTime);
   f=fopen(FFTCash[Cache].Name,"rb");
   StopTimer(DiskIOTime);
   if (f==NULL)
     FatalError("Unable to open '%s' for convolution.\n",FFTCash[Cache].Name);

   while (L2)
     {size_t Sz;
      Sz=Min(FIXEDBUF_SIZE/sizeof(FFT_DATA_TYPE),L2);
      StartTimer(DiskIOTime);
      fread(Num2,sizeof(FFT_DATA_TYPE),Sz,f);
      StopTimer(DiskIOTime);
//      for (x=0;x<Sz;x++) CmplxModMul(&Num1[x],&Num1[x],&Num2[x]);
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
DoFGTMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t NumLen,
         size_t ProdLen, INT32 Scale)
{
 size_t x, NumLen2 = NumLen * 2;
 size_t FGTLen=CalcFGTLen(NumLen);
 FFT_DATA_TYPE *FFTNum2=FFTNum;
 int BothInMem=0;
 INT32 ScaleCarry=0;

/* Can we do both numbers in memory? */
  if (FGTLen*sizeof(ModInt)*2 <= CoreMemAvail)
    {BothInMem=-2;FFTNum2=FFTNum+FGTLen;}

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
          SaveFFTIntoCache(FFTNum,FGTLen,Num2,NumLen,
                           SaveNum1FFT+SaveNum2FFT,0);
       }
     DoConvolutions(-1,FGTLen);
    }
  else
    {int Line1=0,Line2=0;
     if (Num1IsCached) Line1=CheckFFTCache(Num1,NumLen,Num1IsCached,0);
     if (Num2IsCached) Line2=CheckFFTCache(Num2,NumLen,Num2IsCached,0);
     if (Line1 && Line2)
       {
        DumpDebug("Caches...");
        LoadFFTFromCache(Line1,FFTNum);
        DoConvolutions(Line2,FGTLen);
       }
     else if (Line1)
       {
        DumpDebug("Cache...F");
        DoFwdTransforms(FFTNum,Num2,NumLen);
        if (SaveNum2FFT)
          SaveFFTIntoCache(FFTNum,FGTLen,Num2,NumLen,SaveNum2FFT,0);
        DoConvolutions(Line1,FGTLen);
       }
     else if (Line2)
       {
        DumpDebug("Cache...F");
        DoFwdTransforms(FFTNum,Num1,NumLen);
        if (SaveNum1FFT)
          SaveFFTIntoCache(FFTNum,FGTLen,Num1,NumLen,SaveNum1FFT,0);
        DoConvolutions(Line2,FGTLen);
       }
     else
       {
        if (SaveNum1FFT)
          {
           DumpDebug("F");
           DoFwdTransforms(FFTNum,Num1,NumLen);
           Line1=SaveFFTIntoCache(FFTNum,FGTLen,Num1,NumLen,SaveNum1FFT,0);
           if (BothInMem) Line1=-2;
           else if (Line1==0) SaveFFTIntoCache0(FFTNum,FGTLen);
           DumpDebug("F");
           DoFwdTransforms(FFTNum2,Num2,NumLen);
           if (SaveNum2FFT)
             SaveFFTIntoCache(FFTNum2,FGTLen,Num2,NumLen,SaveNum2FFT,0);
           DoConvolutions(Line1,FGTLen);
          }
        else if (SaveNum2FFT)
          {
           DumpDebug("F");
           DoFwdTransforms(FFTNum,Num2,NumLen);
           Line2=SaveFFTIntoCache(FFTNum,FGTLen,Num2,NumLen,SaveNum2FFT,0);
           if (BothInMem) Line2=-2;
           else if (Line2==0) SaveFFTIntoCache0(FFTNum,FGTLen);
           DumpDebug("F");
           DoFwdTransforms(FFTNum2,Num1,NumLen);
           if (SaveNum1FFT)
             SaveFFTIntoCache(FFTNum2,FGTLen,Num1,NumLen,SaveNum1FFT,0);
           DoConvolutions(Line2,FGTLen);
          }
        else
          {
           Line2=BothInMem;
           DumpDebug("F");
           DoFwdTransforms(FFTNum,Num1,NumLen);
           if (!BothInMem) SaveFFTIntoCache0(FFTNum,FGTLen);
           DumpDebug("F");
           DoFwdTransforms(FFTNum2,Num2,NumLen);
           DoConvolutions(Line2,FGTLen);
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
   for (x = 0; x < FGTLen; x++)
      {
       ModMul(&Pyramid,&FFTNum[x].r,&MulInv);
       RAW_Add(&Carry,&Pyramid,&Carry);
       ProdNdx--;ProdPos--;
       ProdBuf[ProdNdx]=RAW_DivInt(&Carry,100000000,&Carry);

       ModMul(&Pyramid,&FFTNum[x].i,&MulInv);
       RAW_Add(&Carry,&Pyramid,&Carry);
       ProdNdx--;ProdPos--;
       ProdBuf[ProdNdx]=RAW_DivInt(&Carry,100000000,&Carry);
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
/* Turn the FGT caching back on. */
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

  SetFGTSize(NumLen);
  Carry=DoFGTMul(Prod,Num1,Num2,NumLen,ProdLen,Scale);

  StopTimer(FFTMulTime);
  FFTDisk+=DiskIOTime;
  DumpDebug("Done fft.\n");

  return Carry;
}



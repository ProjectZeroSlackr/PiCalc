#include "pi.h"

#include "block.h"
#include "ntt.h"
#include "bigmul.h"
#include "cache.h"

#include "modmath.h"
#include "vector.h"
#include "crt.h"
#include "primes.h"

static ModInt *FFTNum=NULL;
static char NTTMergeFileName[MAX_FILENAME+16];
       CRTNum PProds[NPrimes+1][CRT_LEN];
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
/* NumLen->NTTLen.  NPrimes of them.  enough to hold both nums at once.*/
return CalcNTTLen(Len)*NPrimes*2*sizeof(ModInt);
}

void
InitFFTMul(size_t Len)
{size_t Pass;

  if (Cfg.HalfMuls) Len/=2; /* This formula will only use half sized muls */
  FFTLimit=Len;
  if (FFTLimit > MaxFFTLen/RawIntDigits) FFTLimit=MaxFFTLen/RawIntDigits;

/* At an absolute _minimum_ I need 'Limit' bytes to do one NTT */
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

/* Prep for the Knuth CRT (Chinese Remainder Theorem). */
  SetCRT(PProds[0],1);
  for (Pass=1;Pass < NPrimes; Pass++)
    {
     CopyCRT(PProds[Pass],PProds[Pass-1]);
     MulCRT(PProds[Pass],PrimeList[Pass-1],CRT_LEN);
    }

  if (ReadCfgItem(PI_CFG_FILE,"Files","NTT_Merge",
                  NTTMergeFileName,Cfg_String,MAX_FILENAME)==0)
     strcpy(NTTMergeFileName,"pinttme.rge");
}

void
DeInitFFTMul(void)
{
}

static void
DoConvolutions(int Cache,size_t FFTLen2,size_t StartP, size_t EndP)
{size_t Pass;
 ModInt *Num1,*Num2;

StartTimer(ConvTime);
if (Cache==-2) /* In memory convolution of both nums */
  {ModInt *FFTNum2=FFTNum+FFTLen2*NPrimes;
   DumpDebug("Cm...");
   for (Pass=StartP;Pass < EndP;Pass++)
     {
      SetModPrime(Pass);
      Num1=FFTNum+FFTLen2*(Pass-StartP);
      Num2=FFTNum2+FFTLen2*(Pass-StartP);
      PrepVector(Prime,FFTLen2);
      VectorModMul( Num1, Num2, FFTLen2);
     }
  }
else if (Cache==-1)  /* In memory self convolution. */
  {
   DumpDebug("C");
   for (Pass=StartP;Pass < EndP;Pass++)
     {
      SetModPrime(Pass);
      Num1=FFTNum+FFTLen2*(Pass-StartP);
      PrepVector(Prime,FFTLen2);
      VectorModMul( Num1, Num1, FFTLen2);
     }
  }
#ifdef VIRTUAL_CACHE
else /* Virtual Mem based convolution. */
  {
   Num1=FFTNum;
   Num2=(ModInt*)FFTCash[Cache].Mem;

   if (Num2==NULL)
     FatalError("Cache %d doesn't exist.\n",Cache);

   DumpDebug("C%d...",Cache);
   for (Pass=StartP;Pass < EndP;Pass++)
     {
      SetModPrime(Pass);
      PrepVector(Prime,FFTLen2);
      VectorModMul( Num1, Num2, FFTLen2);
      Num1+=FFTLen2;Num2+=FFTLen2;
     }
  }
#else
else /* Disk based convolution. */
  {FILE *f;
   Num1=FFTNum;
   Num2=(ModInt*)FixedBuf;

   DumpDebug("C%d...",Cache);
   StartTimer(DiskIOTime);
   f=fopen(FFTCash[Cache].Name,"rb");
   StopTimer(DiskIOTime);
   if (f==NULL)
     FatalError("Unable to open '%s' for convolution.\n",FFTCash[Cache].Name);

   for (Pass=StartP; Pass < EndP;Pass++)
     {size_t L2=FFTLen2;
      SetModPrime(Pass);
      while (L2)
        {size_t Sz;
         Sz=Min(FIXEDBUF_SIZE/sizeof(ModInt),L2);
         StartTimer(DiskIOTime);
         fread(Num2,sizeof(ModInt),(size_t)Sz,f);
         StopTimer(DiskIOTime);
         PrepVector(Prime,Sz);
         VectorModMul( Num1, Num2, Sz);
         Num1+=Sz;
         L2-=Sz;
        }
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
DoNTTMuls(BigInt Prod, BigInt Num1, BigInt Num2, size_t NumLen,
          size_t ProdLen, INT32 Scale, size_t Parts)
{
 size_t x, NumLen2 = NumLen * 2;
 size_t FFTLen2=CalcNTTLen(NumLen);
 size_t Pass;
 INT32 ScaleCarry=0;
 FILE *f=NULL;
 ModInt *FFTNum2=FFTNum;
 int BothInMem=0,UseDisk=0;

/* Do we have to use the disk to do all 'NPrimes' NTT? */
  if (Parts != NPrimes) UseDisk=1;

  if (UseDisk)
    {
     StartTimer(DiskIOTime);
     remove(NTTMergeFileName);
     f=fopen(NTTMergeFileName,"wb+");
     StopTimer(DiskIOTime);
     if (f==NULL)
       FatalError("Unable to open the FFT merge file.\n");
    }
/* Can we do both numbers in memory? */
  else if (FFTLen2*NPrimes*2*sizeof(ModInt) <= CoreMemAvail)
    {BothInMem=-2;FFTNum2=FFTNum+FFTLen2*NPrimes;}

  for (Pass=0; Pass < NPrimes;Pass+=Parts)
    {
     SetModPrime(Pass);
     DumpDebug("%d",Pass);
     if (Num1 == Num2)
       {int Line=0;
        if (Num1IsCached || Num2IsCached)
          Line=CheckFFTCache(Num2,NumLen,Num1IsCached + Num2IsCached,Pass);
        if (Line) LoadFFTFromCache(Line,FFTNum);
        else
          {
           DumpDebug("F");
           DoFwdTransforms(FFTNum,Num1,NumLen,Pass,Pass+Parts);
           if (SaveNum2FFT || SaveNum1FFT)
             SaveFFTIntoCache(FFTNum,FFTLen2*Parts,Num2,NumLen,
                              SaveNum1FFT+SaveNum2FFT,Pass);
          }
        DoConvolutions(-1,FFTLen2,Pass,Pass+Parts);
       }
     else
       {int Line1=0,Line2=0;
        if (Num1IsCached) Line1=CheckFFTCache(Num1,NumLen,Num1IsCached,Pass);
        if (Num2IsCached) Line2=CheckFFTCache(Num2,NumLen,Num2IsCached,Pass);
        if (Line1 && Line2)
          {
           DumpDebug("Caches...");
           LoadFFTFromCache(Line1,FFTNum);
           DoConvolutions(Line2,FFTLen2,Pass,Pass+Parts);
          }
        else if (Line1)
          {
           DumpDebug("Cache...F");
           DoFwdTransforms(FFTNum,Num2,NumLen,Pass,Pass+Parts);
           if (SaveNum2FFT)
             SaveFFTIntoCache(FFTNum,FFTLen2*Parts,Num2,NumLen,SaveNum2FFT,Pass);
           DoConvolutions(Line1,FFTLen2,Pass,Pass+Parts);
          }
        else if (Line2)
          {
           DumpDebug("Cache...F");
           DoFwdTransforms(FFTNum,Num1,NumLen,Pass,Pass+Parts);
           if (SaveNum1FFT)
             SaveFFTIntoCache(FFTNum,FFTLen2*Parts,Num1,NumLen,SaveNum1FFT,Pass);
           DoConvolutions(Line2,FFTLen2,Pass,Pass+Parts);
          }
        else
          {
           if (SaveNum1FFT)
             {
              DumpDebug("F");
              DoFwdTransforms(FFTNum,Num1,NumLen,Pass,Pass+Parts);
              Line1=SaveFFTIntoCache(FFTNum,FFTLen2*Parts,Num1,NumLen,SaveNum1FFT,Pass);
              if (BothInMem) Line1=-2;
              else if (Line1==0) SaveFFTIntoCache0(FFTNum,FFTLen2*Parts);
              DumpDebug("F");
              DoFwdTransforms(FFTNum2,Num2,NumLen,Pass,Pass+Parts);
              if (SaveNum2FFT)
                SaveFFTIntoCache(FFTNum2,FFTLen2*Parts,Num2,NumLen,SaveNum2FFT,Pass);
              DoConvolutions(Line1,FFTLen2,Pass,Pass+Parts);
             }
           else if (SaveNum2FFT)
             {
              DumpDebug("F");
              DoFwdTransforms(FFTNum,Num2,NumLen,Pass,Pass+Parts);
              Line2=SaveFFTIntoCache(FFTNum,FFTLen2*Parts,Num2,NumLen,SaveNum2FFT,Pass);
              if (BothInMem) Line2=-2;
              else if (Line2==0) SaveFFTIntoCache0(FFTNum,FFTLen2*Parts);
              DumpDebug("F");
              DoFwdTransforms(FFTNum2,Num1,NumLen,Pass,Pass+Parts);
              if (SaveNum1FFT)
                SaveFFTIntoCache(FFTNum2,FFTLen2*Parts,Num1,NumLen,SaveNum1FFT,Pass);
              DoConvolutions(Line2,FFTLen2,Pass,Pass+Parts);
             }
           else
             {
              Line2=BothInMem;
              DumpDebug("F");
              DoFwdTransforms(FFTNum,Num1,NumLen,Pass,Pass+Parts);
              if (!BothInMem) SaveFFTIntoCache0(FFTNum,FFTLen2*Parts);
              DumpDebug("F");
              DoFwdTransforms(FFTNum2,Num2,NumLen,Pass,Pass+Parts);
              DoConvolutions(Line2,FFTLen2,Pass,Pass+Parts);
             }
          }
       }

     DeleteFFTCache(0); /* get rid of the convolution file */
     DumpDebug("R");
     DoRevTransforms(FFTNum,NumLen,Pass,Pass+Parts);

     if (UseDisk)
       {long Old;
        StartTimer(SaveTime);
        StartTimer(DiskIOTime);
        fflush(f);
        Old=ftell(f);
        x=fwrite(FFTNum,sizeof(ModInt),FFTLen2*Parts,f);
        if (x!=FFTLen2*Parts)
          {
           clearerr(f);
           fprintf(stderr,"\n** NOTICE **  I had to delete non-critical, performance improving\n");
           fprintf(stderr,"FFT data to be able to operate.  More disk space is strongly recommended.\n");
           FlushFFTCache(0);
           fseek(f,Old,SEEK_SET);
           x=fwrite(FFTNum,sizeof(ModInt),FFTLen2*Parts,f);
           if (x!=FFTLen2*Parts)
             FatalError("\nUnable to write to the very important FFT merge file %s\nProbably need more disk space.",NTTMergeFileName);
          }
        fflush(f);
        StopTimer(SaveTime);
        StopTimer(DiskIOTime);
       }
    }

  DumpDebug("Carries...");
  {CRTNum Pyramid[CRT_LEN], CRTResult[CRT_LEN];
   INT32 *ProdBuf=(INT32*)FixedBuf;
   size_t ProdNdx,ProdPos,ProdBufLen=FIXEDBUF_SIZE/sizeof(INT32);
   size_t FFTBlockLen,FFTBlock,FFTNdx;
   ModInt *FFTParts[NPrimes];

   FFTBlockLen=FFTLen2;
   if (UseDisk) FFTBlockLen=Min(FFTBlockLen,CoreMemAvail/(NPrimes*sizeof(ModInt)));
   for (Pass=0;Pass < NPrimes; Pass++) FFTParts[Pass]=FFTNum+FFTBlockLen*Pass;
   ClearCRT(Pyramid);
   ProdNdx=ProdBufLen;ProdPos=NumLen2;
   for (FFTBlock=0;FFTBlock < FFTLen2;FFTBlock+=FFTBlockLen)
     {
      if (UseDisk)
        {
         for (Pass=0;Pass < NPrimes; Pass++)
           {size_t z;
            StartTimer(CRTLoad);
            StartTimer(DiskIOTime);
            fseek(f,(FFTLen2*Pass+FFTBlock)*sizeof(ModInt),SEEK_SET);
            z=fread(FFTParts[Pass],sizeof(ModInt),(size_t)FFTBlockLen,f);
            StopTimer(CRTLoad);
            StopTimer(DiskIOTime);
            if (z!=FFTBlockLen)
              FatalError("\nError reading the very important FFT merge file %s %lu %lu %lu\n",
                         NTTMergeFileName,(ULINT)z,(ULINT)FFTBlockLen,(ULINT)FFTBlock);
           }
        }

#ifdef VECTOR_CRT
       StartTimer(CRTTime);
       for (Pass = 0; Pass < NPrimes; Pass++)
         {size_t z;
          SetModPrime(Pass);
          PrepVector(Prime,FFTBlockLen);
          VectorModMulC( FFTParts[Pass], Consts[Pass].MulInv, FFTBlockLen );

          for (z=0;z < Pass;z++)
            {
             VectorModSub(FFTParts[Pass],FFTParts[z],FFTBlockLen);
             VectorModMulC( FFTParts[Pass], Inverses[z][Pass], FFTBlockLen );
            }
         }
      StopTimer(CRTTime);
      for (FFTNdx=0;FFTNdx<FFTBlockLen;FFTNdx++)
        {
         StartTimer(CRTTime);
         CRTResult[NPrimes-1]=FFTParts[0][FFTNdx];
         MulCRT1(CRTResult,PProds[1],FFTParts[1][FFTNdx]);
         MulCRT2(CRTResult,PProds[2],FFTParts[2][FFTNdx]);
         MulCRT3(CRTResult,PProds[3],FFTParts[3][FFTNdx]);
#if NPrimes==8
         MulCRT4(CRTResult,PProds[4],FFTParts[4][FFTNdx]);
         MulCRT5(CRTResult,PProds[5],FFTParts[5][FFTNdx]);
         MulCRT6(CRTResult,PProds[6],FFTParts[6][FFTNdx]);
         MulCRT7(CRTResult,PProds[7],FFTParts[7][FFTNdx]);
#endif

         StopTimer(CRTTime);
         StartTimer(CarryTime);
         Add2CRT(Pyramid,CRTResult);
         ProdNdx-=RawModInts;ProdPos-=RawModInts;
         {int q;
          for (q=RawModInts-1;q>=0;q--)
             ProdBuf[ProdNdx+q] = StripCRT(Pyramid);
         }
         StopTimer(CarryTime);
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
#endif

#ifndef VECTOR_CRT
      for (FFTNdx=0;FFTNdx<FFTBlockLen;FFTNdx++)
        {
/* The Chinese Remainder Thereom */
         StartTimer(CRTTime);
#ifdef UNROLLED_CRT
         {ModInt v[NPrimes+1];
          SetModPrime(0);
          v[1]=ModMul(FFTParts[0][FFTNdx],Consts[0].MulInv);
          CRTResult[NPrimes-1]=v[1];

          SetModPrime(1);
          v[2]=ModMul(FFTParts[1][FFTNdx],Consts[1].MulInv);
          v[2]=ModMul(ModSub(v[2],v[1]),Inverses[0][1]);/* c12 */
          MulCRT1(CRTResult,PProds[1],v[2]);

          SetModPrime(2);
          v[3]=ModMul(FFTParts[2][FFTNdx],Consts[2].MulInv);
          v[3]=ModMul(ModSub(v[3],v[1]),Inverses[0][2]); /* c13 */
          v[3]=ModMul(ModSub(v[3],v[2]),Inverses[1][2]); /* c23 */
          MulCRT2(CRTResult,PProds[2],v[3]);

          SetModPrime(3);
          v[4]=ModMul(FFTParts[3][FFTNdx],Consts[3].MulInv);
          v[4]=ModMul(ModSub(v[4],v[1]),Inverses[0][3]); /* c14 */
          v[4]=ModMul(ModSub(v[4],v[2]),Inverses[1][3]); /* c24 */
          v[4]=ModMul(ModSub(v[4],v[3]),Inverses[2][3]); /* c34 */
          MulCRT3(CRTResult,PProds[3],v[4]);
#if NPrimes==8
          SetModPrime(4);
          v[5]=ModMul(FFTParts[4][FFTNdx],Consts[4].MulInv);
          v[5]=ModMul(ModSub(v[5],v[1]),Inverses[0][4]); /* c15 */
          v[5]=ModMul(ModSub(v[5],v[2]),Inverses[1][4]); /* c25 */
          v[5]=ModMul(ModSub(v[5],v[3]),Inverses[2][4]); /* c35 */
          v[5]=ModMul(ModSub(v[5],v[4]),Inverses[3][4]); /* c45 */
          MulCRT4(CRTResult,PProds[4],v[5]);

          SetModPrime(5);
          v[6]=ModMul(FFTParts[5][FFTNdx],Consts[5].MulInv);
          v[6]=ModMul(ModSub(v[6],v[1]),Inverses[0][5]); /* c16 */
          v[6]=ModMul(ModSub(v[6],v[2]),Inverses[1][5]); /* c26 */
          v[6]=ModMul(ModSub(v[6],v[3]),Inverses[2][5]); /* c36 */
          v[6]=ModMul(ModSub(v[6],v[4]),Inverses[3][5]); /* c46 */
          v[6]=ModMul(ModSub(v[6],v[5]),Inverses[4][5]); /* c56 */
          MulCRT5(CRTResult,PProds[5],v[6]);

          SetModPrime(6);
          v[7]=ModMul(FFTParts[6][FFTNdx],Consts[6].MulInv);
          v[7]=ModMul(ModSub(v[7],v[1]),Inverses[0][6]); /* c17 */
          v[7]=ModMul(ModSub(v[7],v[2]),Inverses[1][6]); /* c27 */
          v[7]=ModMul(ModSub(v[7],v[3]),Inverses[2][6]); /* c37 */
          v[7]=ModMul(ModSub(v[7],v[4]),Inverses[3][6]); /* c47 */
          v[7]=ModMul(ModSub(v[7],v[5]),Inverses[4][6]); /* c57 */
          v[7]=ModMul(ModSub(v[7],v[6]),Inverses[5][6]); /* c67 */
          MulCRT6(CRTResult,PProds[6],v[7]);

          SetModPrime(7);
          v[8]=ModMul(FFTParts[7][FFTNdx],Consts[7].MulInv);
          v[8]=ModMul(ModSub(v[8],v[1]),Inverses[0][7]); /* c18 */
          v[8]=ModMul(ModSub(v[8],v[2]),Inverses[1][7]); /* c28 */
          v[8]=ModMul(ModSub(v[8],v[3]),Inverses[2][7]); /* c38 */
          v[8]=ModMul(ModSub(v[8],v[4]),Inverses[3][7]); /* c48 */
          v[8]=ModMul(ModSub(v[8],v[5]),Inverses[4][7]); /* c58 */
          v[8]=ModMul(ModSub(v[8],v[6]),Inverses[5][7]); /* c68 */
          v[8]=ModMul(ModSub(v[8],v[7]),Inverses[6][7]); /* c78 */
          MulCRT7(CRTResult,PProds[7],v[8]);
#endif
         }
#else
         {ModInt v[NPrimes];size_t y;
          for (Pass = 0; Pass < NPrimes; Pass++)
            {UINT32 u;size_t z;
             SetModPrime(Pass);
             u=FFTParts[Pass][FFTNdx];
             u=ModMul(u,Consts[Pass].MulInv); /* Normalize */
             for (z=0;z < Pass;z++)
                u=ModMul(ModSub(u,v[z]),Inverses[z][Pass]);
             v[Pass]=u % Consts[Pass].Prime;
            }
          CRTResult[CRT_LEN-1]=v[0];
          for (y=1;y<NPrimes;y++) MulCRTz(CRTResult,PProds[y],v[y],y);
         }
#endif

         Add2CRT(Pyramid,CRTResult);
         StopTimer(CRTTime);
         StartTimer(CarryTime);
         ProdNdx-=RawModInts;ProdPos-=RawModInts;
         {int q;
          for (q=RawModInts-1;q>=0;q--)
             ProdBuf[ProdNdx+q] =(INT32) StripCRT(Pyramid);
         }
         StopTimer(CarryTime);
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
#endif /* not VECTOR_CRT */

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
   for (x=0;x<NPrimes;x++)
     if (Pyramid[x])
       {
        DumpCRT(Pyramid);
        FatalError("The multiplication pyramid wasn't zero.\n");
       }
  }

  if (UseDisk)
    {
     StartTimer(DiskIOTime);
     fclose(f);remove(NTTMergeFileName);
     StopTimer(DiskIOTime);
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
FFTMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t NumLen, size_t ProdLen, INT32 Scale)
/*
** Scale=0 means don't multiply by Scale.
** Scale!=0 means do multiply by Scale.  (Usually 10.)
**
** If a scaling causes a carry, that carry will be returned, else
** zero will be returned.
*/
{INT32 Carry=0;
 size_t Parts;

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
** This uses a radix-2 FFT, so the length has to be a power of two.
*/
  if (!IsPow2(NumLen))
    FatalError("The FFT size is not a power of two\n");

  DumpDebug("FFT %s ",Num2Str(NumLen*RawIntDigits));
  StartTimer(FFTMulTime);
  FFTDisk-=DiskIOTime;
  SetNTTSize(NumLen);

/* How many parts can we do in memory?  Allow for double length conv. */
  for (Parts=NPrimes;Parts;Parts/=2)
    if (CalcNTTLen(NumLen)*Parts*sizeof(ModInt) <= CoreMemAvail)
      break;

  if (Parts==0)
    FatalError("The FFTMul was called with too large of a number.\n");

  Carry=DoNTTMuls(Prod,Num1,Num2,NumLen,ProdLen,Scale,Parts);

  StopTimer(FFTMulTime);
  FFTDisk+=DiskIOTime;
  DumpDebug("Done fft.\n");
  return Carry;
}




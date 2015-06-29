#include "pi.h"
#include "fft.h"
#include "bigmul.h"
#include "cache.h"
#include "block.h"

static FFT_DATA_TYPE *FFTNum=NULL;
static BigInt FMWork=NO_NUM; /* for FractalMul() */
static size_t FFTLimit=0;
static double MaxFFTError;

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

/*
** Reorder complex data array by bit reversal rule.
*/
static void
F_FFTReOrder(double *point, int Len)
{
  int Index, xednI, x;
  double temp;

  xednI = 0;
  for (Index = 0; Index < Len; Index += 2)
    {
      if (Index < xednI)
        {
          temp = point[Index];
          point[Index] = point[xednI];
          point[xednI] = temp;
          temp = point[Index + 1];
          point[Index + 1] = point[xednI + 1];
          point[xednI + 1] = temp;
        }
      x = Len / 2;
      while ((x > 0) && (x <= xednI))
        {
          xednI -= x;
          x /= 2;
        }
      xednI += x;
    }
}


/*
**  This is fairly much a standard Fast Fourier Transform, except
**  I don't do the 'normalization' at the end of the inverse FFT.
**  Due to memory access times, it's faster to do it later, when
**  we are accessing the memory anyway.  It has been tweaked a bit
**  for my compiler, but nothing that should cause your compiler
**  any problems.  It wouldn't be hard to convert it back to array
**  indexing instead of pointers.
*/
static void
F_SpecialFFT(double *point, int Len, int Dir)
{
  int Step, HalfStep;
  double tcos, tsin;
  double Cosine, Sine;
  double *LastPoint, *ep, *bp;
  double *ap, *fp;
  double Theta;

/*
** We are 'Len' of 'complex' data points.  Since each complex data
** point is actually two doubles, that means the 'double' length is
** Len*2
*/
  Len *= 2;

/*
** Reorder complex data vector by bit reversal rule
*/
  F_FFTReOrder(point, Len);

  LastPoint = &point[Len];
  Theta = atan(1.0) * 4.0;
  if (Dir < 0)
    Theta = -Theta;

/*
** Do the transform
*/
  Step = 2;
  while (Step < Len)
    {
      HalfStep = Step;
      Step *= 2;
      ep = &point[HalfStep];

/*
** The first pass of this loop down below will always have
** Cosine=1 and Sine=0, so we can optimize this a bit.
*/
      for (ap = &point[0]; ap < LastPoint; ap += Step)
        {
          double tr, ti, t;
          fp = ap + HalfStep;
          tr = *fp;
          *fp = (t = *ap) - tr;
          *ap = t + tr;

          ti = *(fp + 1);
          *(fp + 1) = (t = *(ap + 1)) - ti;
          *(ap + 1) = t + ti;
        }

      Cosine = 1.0;
      Sine = 0.0;
      tsin = sin(Theta);
      tcos = sin(Theta *= 0.5);
      tcos = -2.0 * tcos * tcos;
      for (bp = &point[2]; bp < ep; bp += 2)
        {
          {
            double temp;
            temp = Cosine;
            Cosine = Cosine * tcos - Sine * tsin + Cosine;
            Sine = Sine * tcos + temp * tsin + Sine;
          }
          for (ap = bp; ap < LastPoint; ap += Step)
            {
              double tr, ti, t;
              fp = ap + HalfStep;
              tr = (ti = *fp) * Cosine - *(fp + 1) * Sine;
              *fp = (t = *ap) - tr;
              *ap = t + tr;

              ti = ti * Sine + *(fp + 1) * Cosine;
              *(fp + 1) = (t = *(ap + 1)) - ti;
              *(ap + 1) = t + ti;
            }
        }
    }
/*
** We would normally 'normalize' the result, but for us, it's
** faster to do it later, when we are accessing the memory anyway.
*/
}

/*
** This is the routine that lets us use a 'Complex' FFT with our
** Real only data.
*/
static void
F_RFFT(double *data, int Len, int Dir)
{
  double Cosine, Sine, tcos, tsin, temp, theta;
  int i, j;
  double NegHalf, PosHalf;

  Len /= 2;
  theta = -(4.0 * atan(1.0)) / Len;
  NegHalf = -0.5;
  PosHalf = 0.5;
  if (Dir > 0)
    {
      NegHalf = 0.5;
      PosHalf = -0.5;
      theta = -theta;
      F_SpecialFFT(data, Len, 1);
    }

  tcos = sin(theta / 2);
  tcos = -2.0 * tcos * tcos;
  tsin = sin(theta);
  Cosine = 1 + tcos;
  Sine = tsin;
  for (i = 2, j = (2 * Len) - i; i < Len; i += 2, j -= 2)
    {
      double t1r, t1i, t2r, t2i, a, b;
      t1r = 0.5 * ((a = data[i]) + (b = data[j]));
      t2i = PosHalf * (a - b);
      t2r = NegHalf * ((a = data[i + 1]) + (b = data[j + 1]));
      t1i = 0.5 * (a - b);
      a = Cosine * t2r;
      b = Sine * t2i;
      data[i] = t1r + a - b;
      data[j] = t1r - a + b;
      a = Cosine * t2i;
      b = Sine * t2r;
      data[i + 1] = a + b + t1i;
      data[j + 1] = a + b - t1i;
      temp = Cosine;
      Cosine = Cosine * tcos - Sine * tsin + Cosine;
      Sine = Sine * tcos + temp * tsin + Sine;
    }

  temp = data[0];
  data[0] = temp + data[1];
  data[1] = temp - data[1];

  if (Dir < 0)
    {
      data[0] /= 2.0;
      data[1] /= 2.0;
      F_SpecialFFT(data, Len, -1);
    }
}

/*
** Do a Fast Fourier Transform based multiplication.  To save some
** time, we are using the fact that our numbers are purely real and
** doing a half sized real FFT for each number, followed by the
** convolution and then the half sized real FFT inverse.
*/
void
SimpleFFTMul(UINT32 *Prod, UINT32 *Num1, UINT32 *Num2, int NumLen,
             int FFTLen, double *FFTNum1, double *FFTNum2)
{
  int x,y;
  double Carry, RawPyramid, Pyramid, PyramidError;
  double inv;
  double MaxFFTError = 0.0;
  int NumLen2=NumLen*2;

  for (y=0,x=NumLen-1;x>=0;x--,y+=2)
    {
     FFTNum1[y  ]=Num1[x] % 10000;
     FFTNum1[y+1]=Num1[x] / 10000;
    }
  for (;y<FFTLen;y++) FFTNum1[y]=0.0;
  F_RFFT(FFTNum1, FFTLen, 1);

/*
** If we are squaring a number, we can save the cost
** of a FFT.
*/
  if (Num1 != Num2)
    {
      for (y=0,x=NumLen-1;x>=0;x--,y+=2)
        {
         FFTNum2[y  ]=Num2[x] % 10000;
         FFTNum2[y+1]=Num2[x] / 10000;
        }
      for (;y<FFTLen;y++) FFTNum2[y]=0.0;
      F_RFFT(FFTNum2, FFTLen, 1);

      /*
      ** Now do the convolution
      */
      FFTNum1[0] = FFTNum1[0] * FFTNum2[0];
      FFTNum1[1] = FFTNum1[1] * FFTNum2[1];
      for (x = 2; x < FFTLen; x += 2)
        {
          double a, b, c, d;
          FFTNum1[x] = (a = FFTNum1[x])     * (b = FFTNum2[x]) -
                       (c = FFTNum1[x + 1]) * (d = FFTNum2[x + 1]);
          FFTNum1[x + 1] = a * d + c * b;
        }
    }
  else
    {
      /*
      ** Now do the convolution
      */
      FFTNum1[0] = FFTNum1[0] * FFTNum1[0];
      FFTNum1[1] = FFTNum1[1] * FFTNum1[1];
      for (x = 2; x < FFTLen; x += 2)
        {
          double a, c;
          a = FFTNum1[x];
          c = FFTNum1[x + 1];
          FFTNum1[x] = a * a - c * c;
          FFTNum1[x + 1] = a * c + c * a;
        }
    }

/*
** Now do an Inverse FFT
*/
  F_RFFT(FFTNum1, FFTLen, -1);

/*
** Now round the results, and release our carries, and store
** the results in the 'prod' array.  Also do the normalization
** we didn't do in the FFT.  And, as usual, it's slightly faster
** to multiply by the reciprocal, instead of doing the division.
*/
  Carry = 0.0;
  inv = 2.0 / FFTLen;
  y=NumLen2-1;
  x=0;
  while (x<FFTLen)
    {UINT32 D1,D2;
      RawPyramid = FFTNum1[x++] * inv + Carry;
      Pyramid = floor(RawPyramid + 0.5);
      PyramidError = fabs(RawPyramid - Pyramid);
      if (PyramidError > MaxFFTError) MaxFFTError = PyramidError;
      Carry = floor(Pyramid / 10000);
      D1 = Pyramid - Carry * 10000;

      RawPyramid = FFTNum1[x++] * inv + Carry;
      Pyramid = floor(RawPyramid + 0.5);
      PyramidError = fabs(RawPyramid - Pyramid);
      if (PyramidError > MaxFFTError) MaxFFTError = PyramidError;
      Carry = floor(Pyramid / 10000);
      D2 = Pyramid - Carry * 10000;

      Prod[y--]=D2*10000+D1;
    }

  /*
  ** Do a bit of 'sanity' error checking.  This value
  ** is based on some testing with a 'test jig' and
  ** personal opinion.  Should be good enough to catch
  ** systems with poor FPUs.  However, if the value ever
  ** actually reaches 0.5, then you know for a FACT that
  ** the FFTMul() has failed.  David Bailey suggested
  ** a value of 0.001, but that's far far to restrictive.
  */
  if (MaxFFTError >= 0.125)
    {
      printf("**WARNING** Max FFT Error: %f\n", MaxFFTError);
      puts("The FFT is approaching its limit.  As long as");
      puts("the FFTLimit is below 8,388,608 the results are");
      puts("probably correct, but you may want to verify");
      puts("the results with an independant value.  If the");
      puts("value is at least 0.5, then the FFTMul() has");
      puts("definetly failed and the results are incorrect.");
      exit(1);
    }
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

  if (FFTLen2 < 16)
    {UINT32 *Buf1=(UINT32*)CoreMemPtr;
     UINT32 *Buf2=Buf1+NumLen;
     UINT32 *DBuf=Buf2+NumLen;
     double *FFTNum1=(double*)(DBuf+NumLen*2);
     int     FFTLen=NumLen*2*2;
     double *FFTNum2=FFTNum1+FFTLen;
     ReadNumIntoBuf(Num1,Buf1,NumLen);
     ReadNumIntoBuf(Num2,Buf2,NumLen);
     BlockClear(DBuf,DBuf+NumLen*2);
     SimpleFFTMul(DBuf,Num1,Num2,NumLen,FFTLen,FFTNum1,FFTNum2);
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
  MaxFFTError=0.0;
  {FFT_DATA_TYPE Carry, Round;
   INT32 *ProdBuf=(INT32*)FixedBuf;
   size_t ProdNdx,ProdPos,ProdBufLen=FIXEDBUF_SIZE/sizeof(INT32);
   double PyramidError;
   int SF=4+(RAW_FFT_DIG*2);
   int q;

   ProdNdx=ProdBufLen;ProdPos=NumLen2;
   F_Clear(&Carry);
   F_Clear(&Round);Round.Data[SF+1]=50000000;
   for (x=0; x<FFTLen2;x++)
     {
      if (FFTNum[x].Data[0] != 0)
        fprintf(stderr,"Warning, FFT appears to be overflow. %ld %d\n",x,FFTNum[x].Data[0]);
      PyramidError=FFTNum[x].Data[SF+1]/1.0e8;
      PyramidError-=0.5;PyramidError=fabs(PyramidError);PyramidError=0.5-PyramidError;
      if (PyramidError > MaxFFTError) MaxFFTError = PyramidError;
      F_Add(&FFTNum[x],&Round,&FFTNum[x]);
      F_Add(&Carry,&FFTNum[x],&Carry);
      ProdNdx-=RAW_FFT_DIG;ProdPos-=RAW_FFT_DIG;
      for (q=RAW_FFT_DIG-1;q>=0;q--)
        {
         ProdBuf[ProdNdx+q]=Carry.Data[SF];
         F_ShiftR(&Carry);
        }
      {int z;for (z=SF+1;z<FIXPOINT_LEN;z++) Carry.Data[z]=0;}
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
/*
   F_Clear(&Carry);
//   F_Clear(&Round);Round.Data[9]=50000000;
//   F_Clear(&Round);Round.Data[7]=50000000;
   F_Clear(&Round);Round.Data[11]=50000000;
   for (x=0; x<FFTLen2;x++)
     {
      if (FFTNum[x].Data[0] != 0)
        fprintf(stderr,"Warning, FFT appears to be overflow. %ld %d\n",x,FFTNum[x].Data[0]);
//      PyramidError=FFTNum[x].Data[9]/1.0e8;
//      PyramidError=FFTNum[x].Data[7]/1.0e8;
      PyramidError=FFTNum[x].Data[11]/1.0e8;
      PyramidError-=0.5;
      PyramidError=fabs(PyramidError);
      PyramidError=0.5-PyramidError;
      if (PyramidError > MaxFFTError) MaxFFTError = PyramidError;
      F_Add(&FFTNum[x],&Round,&FFTNum[x]);
      F_Add(&Carry,&FFTNum[x],&Carry);
//      P[x]=Carry.Data[8];
//      P[x]=Carry.Data[6];
      P[x]=Carry.Data[10];
      F_ShiftR(&Carry);
//      {int z;for (z=9;z<FIXPOINT_LEN;z++) Carry.Data[z]=0;}
//      {int z;for (z=7;z<FIXPOINT_LEN;z++) Carry.Data[z]=0;}
      {int z;for (z=11;z<FIXPOINT_LEN;z++) Carry.Data[z]=0;}
     }
   StopTimer(CarryTime);
   for (x=0;x<FFTLen2/2;x++) {UINT32 t=P[x];P[x]=P[FFTLen2-x-1];P[FFTLen2-x-1]=t;}
   if (Scale) ScaleCarry=BlockMulBy((INT32*)P,(INT32*)P,Scale,0,Min(ProdLen,NumLen2));
   WriteBufIntoNum((INT32*)P,Prod,Min(ProdLen,NumLen2));
*/
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



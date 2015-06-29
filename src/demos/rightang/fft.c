#include "pi.h"
#include "fft.h"
#include "cache.h"

typedef struct {FFT_DATA_TYPE r,i;} Cmplx;
#define POS (+1.0)
#define NEG (-1.0)

/* A couple of Konstants */
LDouble K_PI_    =3.14159265358979323846L;
LDouble K_SQRT05_=0.70710678118654752440L;

//#define DO_ZERO_PAD_CUT 1

#define CmplxMul(r1,i1,r2,i2) \
   {FFT_DATA_TYPE t=r1;       \
    r1=(r1)*(r2)-(i1)*(i2);   \
    i1=(i1)*(r2)+( t)*(i2);   \
   }

#ifdef  FFT_TRIG_SIN
/* use explicit sin() calls.  Works on all processors & compilers */
#define TRIG_VARS                   \
 LDouble Theta,Angle;               \
 size_t TLen,TNdx;                  \
 LDouble Pow_r,Pow_i;

#define INIT_TRIG(LENGTH,DIR)       \
 Pow_r=1.0;Pow_i=0.0;               \
 Theta=K_PI_/(LENGTH);              \
 TLen=LENGTH;TNdx=0;                \
 Angle=0.0;

#define NEXT_TRIG_POW(DIR)          \
 Angle=(K_PI_*(++TNdx))/TLen;       \
 Pow_r=SINE(Angle*0.5);             \
 Pow_r=1.0-2.0*Pow_r*Pow_r;         \
 Pow_i=SINE(Angle)*(DIR);
#endif

#ifdef  FFT_TRIG_SIN_RECUR
/* Slightly faster than explicit trig.  Only occasionally calls sin() */
#define TRIG_VARS                   \
 size_t TLen,TNdx;                  \
 LDouble Pow_r,Pow_i,Nth_r,Nth_i;

#define INIT_TRIG(LENGTH,DIR)       \
 TNdx=0;TLen=LENGTH;                \
 Pow_r=1.0;Pow_i=0.0;               \
 Nth_r=SINE(K_PI_/((LENGTH)*2));    \
 Nth_r=-2.0*Nth_r*Nth_r;            \
 Nth_i=SINE(K_PI_/(LENGTH))*(DIR);

#define NEXT_TRIG_POW(DIR)          \
 if (((++TNdx)&7)==0)               \
   {LDouble Angle;                  \
    Angle=(K_PI_*(TNdx))/TLen;      \
    Pow_r=SINE(Angle*0.5);          \
    Pow_r=1.0-2.0*Pow_r*Pow_r;      \
    Pow_i=SINE(Angle)*(DIR);        \
   }                                \
 else                               \
   {LDouble temp;                   \
     temp = Pow_r;                  \
     Pow_r = Pow_r * Nth_r - Pow_i  \
             * Nth_i + Pow_r;       \
     Pow_i = Pow_i * Nth_r + temp   \
             * Nth_i + Pow_i;       \
   }
#endif

#ifdef  FFT_WIDE_FPU_RECUR
/*
** Use a simple trig recurance and only do sin() at init.
** Requires 64+ bit FPU registers and a compiler that will keep the vars
** in the FPU registers.
*/
#define TRIG_VARS                   \
 LDouble Pow_r,Pow_i,Nth_r,Nth_i;

#define INIT_TRIG(LENGTH,DIR)       \
 Pow_r=1.0;Pow_i=0.0;               \
 Nth_r=SINE(K_PI_/((LENGTH)*2));    \
 Nth_r=-2.0*Nth_r*Nth_r;            \
 Nth_i=SINE(K_PI_/(LENGTH))*(DIR);

#define NEXT_TRIG_POW(DIR)          \
 {LDouble temp;                     \
   temp = Pow_r;                    \
   Pow_r = Pow_r * Nth_r - Pow_i    \
           * Nth_i + Pow_r;         \
   Pow_i = Pow_i * Nth_r + temp     \
           * Nth_i + Pow_i;         \
 }
#endif



/***********************************
** Decimation in time transforms. **
************************************/

/* The iterative butterfly loop. */
#define ButterflyLoop(Ndx,Sr,Pr,Si,Pi)                   \
   {Cmplx *Left,*Right;                                  \
    Left=&Data[(Ndx)];                                   \
    Right=Left+Step2;                                    \
    while (Left < LastPoint)                             \
      {                                                  \
       {LDouble tr, ti, t;                               \
        tr = (ti = Right->r) * (Pr) * (Sr) -             \
             Right->i * (Pi) * (Si);                     \
        Right->r = (t = Left->r) - tr;                   \
        Left->r  = t + tr;                               \
        ti = ti * (Pi) * (Si) +                          \
             Right->i * (Pr) *(Sr);                      \
        Right->i = (t = Left->i) - ti;                   \
        Left->i  = t + ti;                               \
       }                                                 \
       Left += Step; Right += Step;                      \
      }                                                  \
   }

/* The butterfly for the recursive FFTs */
#define RButterfly(Ndx,Sr,Pr,Si,Pi)                   \
   {size_t z=(Ndx);                                   \
    {LDouble tr, ti, t;                               \
     tr = (ti = Right[z].r) * (Pr) * (Sr) -           \
          Right[z].i * (Pi) * (Si);                   \
     Right[z].r = (t = Left[z].r) - tr;               \
     Left[z].r  = t + tr;                             \
     ti = ti * (Pi) * (Si) +                          \
          Right[z].i * (Pr) *(Sr);                    \
     Right[z].i = (t = Left[z].i) - ti;               \
     Left[z].i  = t + ti;                             \
    }                                                 \
   }

static void
RevFFT_T(Cmplx *Data, size_t Len)
/* Reverse Decimation in Time */
{
  size_t Step, Step2, Step4, Step8;
  size_t b;
  Cmplx *LastPoint=&Data[Len];
  TRIG_VARS;

  StartTimer(FFTITime);
  Step = 1;
  while (Step < Len)
    {
      Step *= 2;
      Step2 = Step/2;
      Step4 = Step/4;
      Step8 = Step/8;

      ButterflyLoop(0,POS,1.0,NEG,0.0);

      INIT_TRIG(Step2,POS);
      for (b = 1; b < Step8; b++)
        {
          NEXT_TRIG_POW(POS);
          ButterflyLoop(b,      POS,Pow_r,NEG,Pow_i);
          ButterflyLoop(Step2-b,NEG,Pow_r,NEG,Pow_i);
          ButterflyLoop(Step4-b,POS,Pow_i,NEG,Pow_r);
          ButterflyLoop(Step4+b,NEG,Pow_i,NEG,Pow_r);
        }
      if (Step8)
        {LDouble sq=K_SQRT05_;
         ButterflyLoop(Step8,      POS,sq,NEG,sq);
         ButterflyLoop(Step2-Step8,NEG,sq,NEG,sq);
        }
      if (Step4)
         ButterflyLoop(Step4,POS,0.0,NEG,1.0);
    }
  StopTimer(FFTITime);
}

static void
RevRFFT_T(Cmplx *Data,size_t Len)
/*
** Recursive Reverse Decimation in Time
*/
{size_t x,Len2,Len4,Len8;
 Cmplx *Left,*Right;
 TRIG_VARS;

Len2 = Len/2;Len4 = Len/4;Len8 = Len/8;
if (Len<=(CPU_CACHE/sizeof(Cmplx))) {RevFFT_T(Data,Len);return;}
if (Len2 >= 2) {RevRFFT_T(Data,Len2);RevRFFT_T(Data+Len2,Len2);}

StartTimer(FFTRTime);
INIT_TRIG(Len2,POS);
Left=&Data[0];
Right=&Data[Len2];
RButterfly(0,POS,1.0,NEG,0.0);
for (x=1;x<Len8;x++)
  {
    NEXT_TRIG_POW(POS);
    RButterfly(x,     POS,Pow_r,NEG,Pow_i);
    RButterfly(Len2-x,NEG,Pow_r,NEG,Pow_i);
    RButterfly(Len4-x,POS,Pow_i,NEG,Pow_r);
    RButterfly(Len4+x,NEG,Pow_i,NEG,Pow_r);
  }

if (Len8)
   {LDouble sq=K_SQRT05_;
    RButterfly(Len8,     POS,sq,NEG,sq);
    RButterfly(Len2-Len8,NEG,sq,NEG,sq);
   }
if (Len4)
    RButterfly(Len4,POS,0,NEG,1.0);
StopTimer(FFTRTime);
}

/*
** Decimation in Frequency transforms
*/
/*
** On the register poor x86, it seems to be faster to do four
** seperate FFT2_?ButterflyLoop()'s than it is to bring the
** four butterflies into one loop.  Your milage may vary.
*/
#define FFT_FButterflyLoop(Ndx,Sr,Pr,Si,Pi)              \
   {Cmplx *Left,*Right;                                  \
    Left=&Data[(Ndx)];                                   \
    Right=Left+Step2;                                    \
    while (Left < LastPoint)                             \
      {                                                  \
       {double temp_r,temp_i,t;                          \
        Left->r  = (temp_r = Left->r) + (t = Right->r);  \
        temp_r  -= t;                                    \
        Left->i  = (temp_i = Left->i) + (t = Right->i);  \
        temp_i  -= t;                                    \
        Right->r = temp_r*(Pr)*(Sr) - temp_i*(Pi)*(Si);  \
        Right->i = temp_i*(Pr)*(Sr) + temp_r*(Pi)*(Si);  \
       }                                                 \
       Left += Step; Right += Step;                      \
      }                                                  \
   }

#define RFFT_FButterfly(Ndx,Sr,Pr,Si,Pi)                    \
    {double temp_r,temp_i,t;size_t z=(Ndx);                 \
     Left[z].r  = (temp_r = Left[z].r) + (t = Right[z].r);  \
     temp_r  -= t;                                          \
     Left[z].i  = (temp_i = Left[z].i) + (t = Right[z].i);  \
     temp_i  -= t;                                          \
     Right[z].r = temp_r*(Pr)*(Sr) - temp_i*(Pi)*(Si);      \
     Right[z].i = temp_i*(Pr)*(Sr) + temp_r*(Pi)*(Si);      \
    }

#define RQUAD 1
#define IQUAD 1

static void
FwdFFT_F(Cmplx *Data, size_t Len)
/*
** Iterative Radix-2 Decimation in Frequency FFT.
*/
{
  size_t Step, Step2, Step4, Step8;
  size_t b;
  Cmplx *LastPoint=&Data[Len];
  TRIG_VARS;

  Step = Len;
  while (Step > 1)
    {
      Step2 = Step/2;
      Step4 = Step/4;
      Step8 = Step/8;

      INIT_TRIG(Step2,POS);
      FFT_FButterflyLoop(0,POS,1.0,POS,0.0);
#ifdef IQUAD
      for (b = 1; b < Step8; b++)
#else
      for (b = 1; b < Step2; b++)
#endif
        {
          NEXT_TRIG_POW(POS);
          FFT_FButterflyLoop(b,      POS,Pow_r,POS,Pow_i);
#ifdef IQUAD
          FFT_FButterflyLoop(Step2-b,NEG,Pow_r,POS,Pow_i);
          FFT_FButterflyLoop(Step4-b,POS,Pow_i,POS,Pow_r);
          FFT_FButterflyLoop(Step4+b,NEG,Pow_i,POS,Pow_r);
#endif
        }
#ifdef IQUAD
      if (Step8)
        {LDouble sq=K_SQRT05_;
         FFT_FButterflyLoop(Step8,      POS,sq,POS,sq);
         FFT_FButterflyLoop(Step2-Step8,NEG,sq,POS,sq);
        }
      if (Step4)
         FFT_FButterflyLoop(Step4,POS,0.0,POS,1.0);
#endif
      Step/=2;
    }
}

static void
FwdRFFT_F(Cmplx *Data,size_t Len)
/*
** Recursive Radix-2 Decimation in Frequency FFT.
*/
{size_t x;
 Cmplx *Left,*Right;
 size_t Len2,Len4,Len8;
 TRIG_VARS

if (Len<=(CPU_CACHE/sizeof(Cmplx))) {FwdFFT_F(Data,Len);return;}

Len2=Len/2;Len4=Len/4;Len8=Len/8;
Left=&Data[0];Right=&Data[Len2];
INIT_TRIG(Len2,POS);

RFFT_FButterfly(0,POS,1.0,POS,0.0);
#ifdef RQUAD
for (x=1;x<Len8;x++)
#else
for (x=1;x<Len2;x++)
#endif
  {
    NEXT_TRIG_POW(POS);
    RFFT_FButterfly(x,     POS,Pow_r,POS,Pow_i);
#ifdef RQUAD
    RFFT_FButterfly(Len2-x,NEG,Pow_r,POS,Pow_i);
    RFFT_FButterfly(Len4-x,POS,Pow_i,POS,Pow_r);
    RFFT_FButterfly(Len4+x,NEG,Pow_i,POS,Pow_r);
#endif
  }
#ifdef RQUAD
if (Len8)
  {LDouble sq=K_SQRT05_;
   RFFT_FButterfly(Len8,     POS,sq,POS,sq);
   RFFT_FButterfly(Len2-Len8,NEG,sq,POS,sq);
  }

if (Len4)
   RFFT_FButterfly(Len4,POS,0.0,POS,1.0);
#endif

if (Len >= 4)
  {
   FwdRFFT_F(Data,     Len2);
   FwdRFFT_F(Data+Len2,Len2);
  }
}



static void
FwdFFT_T(Cmplx *Data, size_t Len)
/* Forward Decimation in Time */
{
  size_t Step, Step2, Step4, Step8;
  size_t b;
  Cmplx *LastPoint=&Data[Len];
  TRIG_VARS;

  StartTimer(FFTITime);
#ifdef DO_ZERO_PAD_CUT
/* save a wee bit of time since the data is zero padded */
  Step = 8;
#else
  Step = 1;
#endif
  while (Step < Len)
    {
      Step *= 2;
      Step2 = Step/2;
      Step4 = Step/4;
      Step8 = Step/8;

      ButterflyLoop(0,POS,1.0,POS,0.0);

      INIT_TRIG(Step2,POS);
      for (b = 1; b < Step8; b++)
        {
          NEXT_TRIG_POW(POS);
          ButterflyLoop(b,      POS,Pow_r,POS,Pow_i);
          ButterflyLoop(Step2-b,NEG,Pow_r,POS,Pow_i);
          ButterflyLoop(Step4-b,POS,Pow_i,POS,Pow_r);
          ButterflyLoop(Step4+b,NEG,Pow_i,POS,Pow_r);
        }
      if (Step8)
        {LDouble sq=K_SQRT05_;
         ButterflyLoop(Step8,      POS,sq,POS,sq);
         ButterflyLoop(Step2-Step8,NEG,sq,POS,sq);
        }
      if (Step4)
         ButterflyLoop(Step4,POS,0.0,POS,1.0);
    }
  StopTimer(FFTITime);
}

static void
FwdRFFT_T(Cmplx *Data,size_t Len)
/*
** Recursive Forward Decimation in Time
*/
{size_t x,Len2,Len4,Len8;
 Cmplx *Left,*Right;
 TRIG_VARS;

Len2 = Len/2;Len4 = Len/4;Len8 = Len/8;
if (Len<=(CPU_CACHE/sizeof(Cmplx))) {FwdFFT_T(Data,Len);return;}
else
if (Len2 >= 2) {FwdRFFT_T(Data,Len2);FwdRFFT_T(Data+Len2,Len2);}

StartTimer(FFTRTime);
INIT_TRIG(Len2,POS);
Left=&Data[0];
Right=&Data[Len2];
RButterfly(0,POS,1.0,POS,0.0);
for (x=1;x<Len8;x++)
  {
    NEXT_TRIG_POW(POS);
    RButterfly(x,     POS,Pow_r,POS,Pow_i);
    RButterfly(Len2-x,NEG,Pow_r,POS,Pow_i);
    RButterfly(Len4-x,POS,Pow_i,POS,Pow_r);
    RButterfly(Len4+x,NEG,Pow_i,POS,Pow_r);
  }

if (Len8)
   {LDouble sq=K_SQRT05_;
    RButterfly(Len8,     POS,sq,POS,sq);
    RButterfly(Len2-Len8,NEG,sq,POS,sq);
   }
if (Len4)
    RButterfly(Len4,POS,0,POS,1.0);
StopTimer(FFTRTime);
}

/*
** Reorder complex data array by bit reversal rule.
** It's faster if we treat it as two doubles, instead
** of a struct Cmplx.
*/
#if 0
static void
FFTReOrder(FFT_DATA_TYPE *point, size_t Len)
{
  size_t Index, xednI, x;
  FFT_DATA_TYPE temp;

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
#endif

void
ShiftFFT(double *FFTNum,int Len,int Dir)
/*
** This routine is the heart of the Right Angle FFT.  It gets called
** before the forward FFT, and after the inverse FFT.  So it's always
** working with ordered data, even if you use a DiF followed by a DiT
** FFT and the convolution is done scrambled.
*/
{int x;
 TRIG_VARS;
Len/=2;

INIT_TRIG(Len*2,Dir);
for (x=1;x<Len;x++)
  {
   NEXT_TRIG_POW(Dir);
   CmplxMul(FFTNum[x*2],FFTNum[x*2+1],Pow_r,Pow_i);
  }
}

static void
RealFFT(FFT_DATA_TYPE *ddata,size_t Len, int Dir)
{
if (Dir > 0)
  {
   ShiftFFT(ddata,Len,Dir);
//   FFTReOrder(ddata,Len);
//   FwdRFFT_T((Cmplx*)ddata,Len/2);
   FwdRFFT_F((Cmplx*)ddata,Len/2);
  }

if (Dir < 0)
  {
//   FFTReOrder(ddata,Len);
   RevRFFT_T((Cmplx*)ddata,Len/2);
   ShiftFFT(ddata,Len,Dir);
  }
}

static void
LoadNumBlockIntoFFT(size_t NX, INT32* NBuf, FFT_DATA_TYPE *FFTData)
{
NBuf+=NX;
while (NX)
  {
////   NX-=RAW_FFT_DIG;
   NX--;
   NBuf--;
   FFTData[0]=*NBuf % 10000;
   FFTData[1]=0.0;
   FFTData[2]=*NBuf / 10000;
   FFTData[3]=0.0;
   FFTData+=4;
  }
}

static void
PutNumIntoFFTNum(FFT_DATA_TYPE *FFTNum, BigInt Num, size_t NumLen)
{size_t NX;
 size_t NLen=NumLen;
 size_t FFTPos=0;
 INT32 *NBuf=(INT32*)FixedBuf;
 StartTimer(LoadTime);
 if (NumLen <8) FatalError("Too small of a FFT.\n");

  while (NLen)
    {
     NX=Min(FIXEDBUF_SIZE/sizeof(INT32),NLen);
     NLen-=NX;
     ReadNumIntoBuf(Num+NLen,NBuf,NX);
     LoadNumBlockIntoFFT(NX,NBuf,FFTNum+FFTPos);
     FFTPos+=(4*NX);
////     FFTPos+=(NX/RAW_FFT_DIG);
    }
  StopTimer(LoadTime);
}

void
FwdTransform(FFT_DATA_TYPE *ddata,BigInt Num,size_t NumLen)
{
PutNumIntoFFTNum(ddata, Num, NumLen);
StartTimer(FFTTime);
RealFFT(ddata,CalcFFTLen(NumLen),1);
StopTimer(FFTTime);
}

void
RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen)
{
StartTimer(FFTTime);
RealFFT(ddata,CalcFFTLen(NumLen),-1);
StopTimer(FFTTime);
}

void
InitFFT(size_t Len)
{
}

void
DeInitFFT(void)
{
}

#undef NextTrigPow
#undef ButterflyLoop
#undef RButterfly
#undef POS
#undef NEG


#include "pi.h"
#include "fft.h"
#include "cache.h"

typedef struct {FFT_DATA_TYPE r,i;} Cmplx;
#define POS (+1.0)
#define NEG (-1.0)

/* A couple of Konstants */
LDouble K_PI_    =3.14159265358979323846L;
LDouble K_SQRT05_=0.70710678118654752440L;

#define DO_ZERO_PAD_CUT 1

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

static void
RealFFT(FFT_DATA_TYPE *ddata,size_t Len, int Dir)
/*
** I don't have a Real Value FFT, so I have to fake it
** by putting a wrapper around a complex FFT.
*/
{size_t i, j, Half;
 Cmplx *Data=(Cmplx*)ddata;
 LDouble NegHalf, PosHalf;
 TRIG_VARS;

Len /= 2;       /* Len/2 Cmplx data */
Half = Len/2;
NegHalf = -0.5;
PosHalf = 0.5;
if (Dir > 0)
  {
    NegHalf = 0.5;
    PosHalf = -0.5;
#ifndef DO_ZERO_PAD_CUT
    FFTReOrder(ddata,Len*2);
#endif
    FwdRFFT_T(Data,Len);
  }

INIT_TRIG(Len,Dir);
for (i = 1, j = Len - i; i < Half; i++, j--)
  {LDouble p1r,p1i,p2r,p2i;
   NEXT_TRIG_POW(Dir);
   /* Seperate the two points from the jumbled points */
   p1r =     0.5 * (Data[i].r + Data[j].r);
   p2i = PosHalf * (Data[i].r - Data[j].r);
   p2r = NegHalf * (Data[i].i + Data[j].i);
   p1i =     0.5 * (Data[i].i - Data[j].i);
   /* Almost a standard Decimation in Time butterfly... */
   {LDouble tmp;
    tmp = p2r;
    p2r = p2r * Pow_r - p2i * Pow_i;
    p2i = p2i * Pow_r + tmp * Pow_i;
   }
   Data[i].r = p1r + p2r;
   Data[i].i = p1i + p2i;
   Data[j].r = p1r - p2r;
   Data[j].i = -(p1i - p2i); /* ... except this is negated */
  }

{LDouble temp;
 temp = Data[0].r;
 Data[0].r = temp + Data[0].i;
 Data[0].i = temp - Data[0].i;
}

if (Dir < 0)
  {
    Data[0].r /= 2.0;
    Data[0].i /= 2.0;
    FFTReOrder(ddata,Len*2);
    RevRFFT_T(Data,Len);
  }
}

static void
ReOrderInt(FFT_INT *Num, size_t Len)
{size_t Index,xednI,k;
 size_t Len2=Len/2;

xednI = 0;
for (Index = 0;Index < Len;Index+=2)
  {
   if (xednI > Index)
     {FFT_INT t;
      t=Num[xednI];
      Num[xednI] = Num[Index];
      Num[Index] = t;
      t=Num[xednI+1];
      Num[xednI+1] = Num[Index+1];
      Num[Index+1] = t;
     }
   k=Len2;
   while ((k <= xednI) && (k >=1)) {xednI-=k;k/=2;}
   xednI+=k;
  }
}

/*
** For loading our numbers into the FFTNum array.
** Just part of a 'zero-cut' optimization.  Since the data
** is short int's, we can take advantage of the cheaper scrambling,
** and the easier manipulation and go ahead and do up to three
** passes of the FFT!
*/
static void
PutNumIntoFFTNum(FFT_DATA_TYPE *FFTNum, BigInt Num, size_t NumLen)
{size_t x;

 StartTimer(LoadTime);
 if (NumLen <8) FatalError("Too small of a FFT.\n");

#ifdef DO_ZERO_PAD_CUT
 {FFT_DATA_TYPE *F=FFTNum;
  INT32 *N=(INT32*)(FFTNum+CalcFFTLen(NumLen));
  FFT_INT *Data;

  N-=NumLen;
  ReadNumIntoBuf(Num,N,NumLen);
#if FFT_DIGITS==4
  BlockUnpack4(N,NumLen);
  NumLen*=2;
#else
  BlockUnpack2(N,NumLen);
  NumLen*=4;
#endif

  Data=(FFT_INT*)N;
  for (x=0;x<NumLen/2;x++)
    {FFT_INT t=Data[x];Data[x]=Data[NumLen-x-1];Data[NumLen-x-1]=t;}
  ReOrderInt(Data,NumLen);
  x=0;
  while (x < NumLen)
    {
     INT32 r1,i1,r2,i2,r3,i3,r4,i4;
     r1=*Data++;i1=*Data++;r2=*Data++;i2=*Data++;
     r3=*Data++;i3=*Data++;r4=*Data++;i4=*Data++;x+=8;
     F[ 0]=(r1+r2) + (r3+r4);
     F[ 8]=(r1+r2) - (r3+r4);
     F[ 4]=(r1-r2) - (i3-i4);
     F[12]=(r1-r2) + (i3-i4);
     F[ 5]=(i1-i2) + (r3-r4);
     F[13]=(i1-i2) - (r3-r4);
     F[ 1]=(i1+i2) + (i3+i4);
     F[ 9]=(i1+i2) - (i3+i4);
     F[ 2]=(r1-i2) + ((r3-i4)-(i3+r4))*K_SQRT05_;
     F[10]=(r1-i2) - ((r3-i4)-(i3+r4))*K_SQRT05_;
     F[ 6]=(r1+i2) - ((r3+i4)+(i3-r4))*K_SQRT05_;
     F[14]=(r1+i2) + ((r3+i4)+(i3-r4))*K_SQRT05_;
     F[ 7]=(i1-r2) + ((r3+i4)-(i3-r4))*K_SQRT05_;
     F[15]=(i1-r2) - ((r3+i4)-(i3-r4))*K_SQRT05_;

     F[ 3]=(i1+r2) + ((i3+r4)+(r3-i4))*K_SQRT05_;
     F[11]=(i1+r2) - ((i3+r4)+(r3-i4))*K_SQRT05_;
     F+=16;
    }
 }
#else
 {INT32 *N;
  size_t FFTLen=CalcFFTLen(NumLen);
  N=(INT32*)(FFTNum+FFTLen);
  ReadNumIntoBuf(Num,N-NumLen,NumLen);
#if FFT_DIGITS==4
  for (x=0; x < FFTLen/2; x+=2)
    {
     --N;
     FFTNum[x  ]=*N % 10000;
     FFTNum[x+1]=*N / 10000;
    }
#else
  for (x=0; x < FFTLen/2; x+=4)
    {INT32 D;
     --N;
     D=*N;
     FFTNum[x  ]= D % 100;D=D / 100;
     FFTNum[x+1]= D % 100;D=D / 100;
     FFTNum[x+2]= D % 100;D=D / 100;
     FFTNum[x+3]= D;
    }
#endif
  for (x = FFTLen/2; x < FFTLen; x++) FFTNum[x] = 0.0;
 }
#endif
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


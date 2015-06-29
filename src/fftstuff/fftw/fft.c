/*
** This is just a shell of what 'fft.c' would look like.
**
** Of course, you would insert your own FFT as the 'RealFFT()'
** function.
*/
#include "pi.h"
#include "fft.h"
#include "cache.h"

/*#error  The FFTW package is not included here.  See the file for more details.*/
/*
** The FFTW package is enormous and is GPL, so it is not included here.
** You will need to supply it seperately, compile it yourself, then link
** this program to it.
**
** Be aware though, that I have never been able to get resonable performance
** from the FFTW package. They call it the "Fast Fft Transform in the West",
** but it's more like the slowest in dullsvile.
*/
#include <fftw.h>

static fftw_plan FPlan[50],RPlan[50];
static int NumPlans;
typedef fftw_complex Cmplx;
static LDouble K_PI_=3.1415926535897932384626433832795028841971;

/*
** For loading our numbers into the FFTNum array.
*/
static void
PutNumIntoFFTNum(FFT_DATA_TYPE *FFTNum, BigInt Num, size_t NumLen)
{size_t x;

 StartTimer(LoadTime);
 if (NumLen <8) FatalError("Too small of a FFT.\n");

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
  StopTimer(LoadTime);
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

static void
RealFFT(FFT_DATA_TYPE *ddata,size_t Len, int Dir)
/*
** I don't have a Real Value FFT, so I have to fake it
** by putting a wrapper around a complex FFT.
*/
{size_t i, j, Half;
 Cmplx *Data=(Cmplx*)ddata;
 double NegHalf, PosHalf;
 TRIG_VARS;

Len /= 2;       /* Len/2 Cmplx data */
Half = Len/2;
NegHalf = -0.5;
PosHalf = 0.5;
if (Dir > 0)
  {
    NegHalf = 0.5;
    PosHalf = -0.5;
    fftw_one(FPlan[Log2(Len)],(FFTW_COMPLEX*)Data,NULL);
  }

INIT_TRIG(Len,Dir);
for (i = 1, j = Len - i; i < Half; i++, j--)
  {LDouble p1r,p1i,p2r,p2i;
   NEXT_TRIG_POW(Dir);
   /* Seperate the two points from the jumbled points */
   p1r =     0.5 * (Data[i].re + Data[j].re);
   p2i = PosHalf * (Data[i].re - Data[j].re);
   p2r = NegHalf * (Data[i].im + Data[j].im);
   p1i =     0.5 * (Data[i].im - Data[j].im);
   /* Almost a standard Decimation in Time butterfly... */
   {LDouble tmp;
    tmp = p2r;
    p2r = p2r * Pow_r - p2i * Pow_i;
    p2i = p2i * Pow_r + tmp * Pow_i;
   }
   Data[i].re = p1r + p2r;
   Data[i].im = p1i + p2i;
   Data[j].re = p1r - p2r;
   Data[j].im = -(p1i - p2i); /* ... except this is negated */
  }

{LDouble temp;
 temp = Data[0].re;
 Data[0].re = temp + Data[0].im;
 Data[0].im = temp - Data[0].im;
}

if (Dir < 0)
  {
    Data[0].re /= 2.0;
    Data[0].im /= 2.0;
    fftw_one(RPlan[Log2(Len)],(FFTW_COMPLEX*)Data,NULL);
  }
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
{size_t Sz;
if (sizeof(fftw_real) != sizeof(FFT_DATA_TYPE))
  FatalError("Wrong fftw_real size.\n");

if (Cfg.HalfMuls) Len/=2; /* This formula will only use half sized muls */
printf("Initializing FFT to %lu\n",(ULINT)CalcFFTLen(Len));

FPlan[0]=RPlan[0]=FPlan[1]=RPlan[1]=NULL;
NumPlans=1;
Sz=4;
while (Sz<=CalcFFTLen(Len))
  {
   FPlan[Log2(Sz)]=fftw_create_plan(Sz,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
   RPlan[Log2(Sz)]=fftw_create_plan(Sz,FFTW_FORWARD ,FFTW_ESTIMATE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
   Sz*=2;
   NumPlans++;
  }
}

void
DeInitFFT(void)
{size_t x;
for (x=0;x<NumPlans;x++)
  {
   if (FPlan[x]) fftw_destroy_plan(FPlan[x]);
   if (RPlan[x]) fftw_destroy_plan(RPlan[x]);
  }
}



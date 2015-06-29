#include "pi.h"
#include "fft.h"
#include "cache.h"

#ifdef LONG_DOUBLE_FFT
#error This external FFT does not support long double.  Only double.
#endif

#ifndef  FFT_WIDE_FPU_RECUR
#error only fft_wide_fpu_recur has been tested yet.
#endif

#ifdef  FFT_TRIG_SIN
/* use explicit sin() calls.  Works on all processors & compilers */
#define TRIG_VARS                   \
 double Theta,Angle;                \
 size_t TLen,TNdx;                  \
 double Pow_r,Pow_i;

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
 double Pow_r,Pow_i,Nth_r,Nth_i;

#define INIT_TRIG(LENGTH,DIR)       \
 TNdx=0;TLen=LENGTH;                \
 Pow_r=1.0;Pow_i=0.0;               \
 Nth_r=SINE(K_PI_/((LENGTH)*2));    \
 Nth_r=-2.0*Nth_r*Nth_r;            \
 Nth_i=SINE(K_PI_/(LENGTH))*(DIR);

#define NEXT_TRIG_POW(DIR)          \
 if (((++TNdx)&7)==0)               \
   {souble Angle;                   \
    Angle=(K_PI_*(TNdx))/TLen;      \
    Pow_r=SINE(Angle*0.5);          \
    Pow_r=1.0-2.0*Pow_r*Pow_r;      \
    Pow_i=SINE(Angle)*(DIR);        \
   }                                \
 else                               \
   {souble temp;                    \
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
 double Pow_r,Pow_i,Nth_r,Nth_i;

#define INIT_TRIG(LENGTH,DIR)       \
 Pow_r=1.0;Pow_i=0.0;               \
 Nth_r=SINE(K_PI_/((LENGTH)*2));    \
 Nth_r=-2.0*Nth_r*Nth_r;            \
 Nth_i=SINE(K_PI_/(LENGTH))*(DIR);

#define NEXT_TRIG_POW(DIR)          \
 {double temp;                      \
   temp = Pow_r;                    \
   Pow_r = Pow_r * Nth_r - Pow_i    \
           * Nth_i + Pow_r;         \
   Pow_i = Pow_i * Nth_r + temp     \
           * Nth_i + Pow_i;         \
 }
#endif


#include "jsp.c"


/*
** Reorder complex data array by bit reversal rule.
*/
void
FFTReOrder(Cmplx *Data, int Len)
{int Index,xednI,k;

xednI = 0;
for (Index = 0;Index < Len;Index++)
  {
   if (xednI > Index)
     {Cmplx Temp;
      Temp=Data[xednI];
      Data[xednI]=Data[Index];
      Data[Index]=Temp;
     }
   k=Len/2;
   while ((k <= xednI) && (k >=1)) {xednI-=k;k/=2;}
   xednI+=k;
  }
}


/*
** This is the routine that lets us use a 'Complex' FFT with our
** Real only data.
*/
static void
RFFT(double *data, int Len, int Dir)
{
  double Cosine, Sine, tcos, tsin, temp, theta;
  int i, j;
  double NegHalf, PosHalf;

  Len /= 2;
  theta = -(4.0 * atan(1.0)) / Len;
  theta = (4.0 * atan(1.0)) / Len; // JSP's directions are different.
  NegHalf = -0.5;
  PosHalf = 0.5;
  if (Dir > 0)
    {
      NegHalf = 0.5;
      PosHalf = -0.5;
      theta = -theta;
      FwdRFFT_F((Cmplx*)data,Len);
      FFTReOrder((Cmplx*)data,Len);
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
      FFTReOrder((Cmplx*)data,Len);
      RevRFFT_T((Cmplx*)data,Len);
    }
}

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

void
FwdTransform(FFT_DATA_TYPE *ddata,BigInt Num,size_t NumLen)
{
PutNumIntoFFTNum(ddata, Num, NumLen);
StartTimer(FFTTime);
RFFT(ddata,CalcFFTLen(NumLen),1);
StopTimer(FFTTime);
}

void
RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen)
{
StartTimer(FFTTime);
RFFT(ddata,CalcFFTLen(NumLen),-1);
StopTimer(FFTTime);
}



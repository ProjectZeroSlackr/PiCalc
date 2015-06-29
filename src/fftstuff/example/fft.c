/*
** This is just a shell of what 'fft.c' would look like.
**
** Of course, you would insert your own FFT as the 'RealFFT()'
** function.
*/
#include "pi.h"
#include "fft.h"
#include "cache.h"

/*
** You can change this, if needed, of coruse.
*/
#ifdef LONG_DOUBLE_FFT
#error This external FFT does not support long double.  Only double.
#endif

void RealFFT(FFT_DATA_TYPE *ddata,size_t Len, int Dir);

/*
** This is just a 'warning' in case you accidently try to use
** it without supplying a FFT for it.
*/
#error You need to supply your own FFT here....

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



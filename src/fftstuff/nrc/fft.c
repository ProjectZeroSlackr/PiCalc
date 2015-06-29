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

/*
** Reorder complex data array by bit reversal rule.
*/
static void
FFTReOrder(double *point, int Len)
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
SpecialFFT(double *point, int Len, int Dir)
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
  FFTReOrder(point, Len);

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
RFFT(double *data, int Len, int Dir)
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
      SpecialFFT(data, Len, 1);
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
      SpecialFFT(data, Len, -1);
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

void
InitFFT(size_t Len)
{
}

void
DeInitFFT(void)
{
}



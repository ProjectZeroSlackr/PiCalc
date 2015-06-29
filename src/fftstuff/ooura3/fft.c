#include "pi.h"
#include "fft.h"
#include "cache.h"

#ifdef LONG_DOUBLE_FFT
#error This external FFT does not support long double.  Only double.
#endif

/*
** This FFT was modified a bit to work with my pi program. CEB.
*/

#include <stdlib.h>
#include <math.h>

static int *IP=NULL;
static double *W=NULL;

#include "fftsg.c"

#if 0
#ifndef USE_DEFAULT_CONVOLUTION

void mp_mul_rcsqu(int n, double *a, int nc, double *c)
{
/* Copyright Ooura. */
    int j, k, kk, ks, m;
    double wkr, wki, xr, xi, yr, yi, akr, aki, ajr, aji;

    ks = (nc << 2) / n;
    kk = 0;
    m = n >> 1;
    for (k = 2; k < m; k += 2) {
        j = n - k;
        kk += ks;
        wkr = 0.5 - c[nc - kk];
        wki = c[kk];
        /* ---- transform CFFT data a[] into RFFT data ---- */
        xr = a[k] - a[j];
        xi = a[k + 1] + a[j + 1];
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        xr = a[k] - yr;
        xi = a[k + 1] - yi;
        yr = a[j] + yr;
        yi = a[j + 1] - yi;
        /* ---- csqu ---- */
        akr = xr * xr - xi * xi;
        aki = 2 * xr * xi;
        ajr = yr * yr - yi * yi;
        aji = 2 * yr * yi;
        /* ---- transform RFFT data axx into CFFT data ---- */
        xr = akr - ajr;
        xi = aki + aji;
        yr = wkr * xr + wki * xi;
        yi = wkr * xi - wki * xr;
        a[k] = akr - yr;
        a[k + 1] = aki - yi;
        a[j] = ajr + yr;
        a[j + 1] = aji - yi;
    }
    xr = a[m];
    xi = a[m + 1];
    a[m] = xr * xr - xi * xi;
    a[m + 1] = 2 * xr * xi;
}


void mp_mul_rcmul(int n, double *a, double *b, int nc, double *c)
{
/* Copyright Ooura. */
    int j, k, kk, ks, m;
    double wkr, wki, xr, xi, yr, yi, akr, aki, ajr, aji, bkr, bki, bjr, bji;

    ks = (nc << 2) / n;
    kk = 0;
    m = n >> 1;
    for (k = 2; k < m; k += 2) {
        j = n - k;
        kk += ks;
        wkr = 0.5 - c[nc - kk];
        wki = c[kk];
        /* ---- transform CFFT data a[] into RFFT data ---- */
        xr = a[k] - a[j];
        xi = a[k + 1] + a[j + 1];
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        akr = a[k] - yr;
        aki = a[k + 1] - yi;
        ajr = a[j] + yr;
        aji = a[j + 1] - yi;
        a[k] = akr;
        a[k + 1] = aki;
        a[j] = ajr;
        a[j + 1] = aji;
        /* ---- transform CFFT data b[] into RFFT data ---- */
        xr = b[k] - b[j];
        xi = b[k + 1] + b[j + 1];
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        xr = b[k] - yr;
        xi = b[k + 1] - yi;
        yr = b[j] + yr;
        yi = b[j + 1] - yi;
        /* ---- cmul ---- */
        bkr = akr * xr - aki * xi;
        bki = akr * xi + aki * xr;
        bjr = ajr * yr - aji * yi;
        bji = ajr * yi + aji * yr;
        /* ---- transform RFFT data bxx into CFFT data ---- */
        xr = bkr - bjr;
        xi = bki + bji;
        yr = wkr * xr + wki * xi;
        yi = wkr * xi - wki * xr;
        b[k] = bkr - yr;
        b[k + 1] = bki - yi;
        b[j] = bjr + yr;
        b[j + 1] = bji - yi;
    }
    xr = a[m];
    xi = a[m + 1];
    yr = b[m];
    yi = b[m + 1];
    b[m] = xr * yr - xi * yi;
    b[m + 1] = xr * yi + xi * yr;
}


void DoRConv(int nfft,double *d1,double *d2, int *ip, double *w)
{double xr,xi;
/* Copyright Ooura. */
    if (nfft > (ip[1] << 2)) {
        makect(nfft >> 2, ip, &w[ip[0]]);
    }
    d2[0] += d1[0];
    xr = d1[1] * d2[1] + d1[2] * d2[2];
    xi = d1[1] * d2[2] + d1[2] * d2[1];
    d2[1] = xr;
    d2[2] = xi;
    if (nfft > 2) {
        mp_mul_rcmul(nfft, &d1[1], &d2[1], ip[1], &w[ip[0]]);
    }
    d2[nfft + 1] *= d1[nfft + 1];
}

void DoRSConv(int nfft,double *d1,int *ip, double *w)
{double xr,xi;
/* Copyright Ooura. */
    if (nfft > (ip[1] << 2)) {
        makect(nfft >> 2, ip, &w[ip[0]]);
    }
    d1[0] *= 2;
    xr = d1[1] * d1[1] + d1[2] * d1[2];
    xi = 2 * d1[1] * d1[2];
    d1[1] = xr;
    d1[2] = xi;
    if (nfft > 2) {
        mp_mul_rcsqu(nfft, &d1[1], ip[1], &w[ip[0]]);
    }
    d1[nfft + 1] *= d1[nfft + 1];
}

/*
** Do a convolution.  This is the standard 'complex' output
** format FFT.  (Numerical Recipes also uses this, although its
** not original or exclusive to them.)
*/
void
DoConvolution(FFT_DATA_TYPE *FFTNum,int Cache,size_t Len2)
{size_t x;
StartTimer(ConvTime);
if (Cache==-2) /* In memory convolution of both nums */
  {FFT_DATA_TYPE *FFTNum2=FFTNum+Len2;
   DumpDebug("Cm...");
   DoRConv(Len2,FFTNum2,FFTNum,IP,W);
  }
else if (Cache==-1)  /* In memory self convolution. */
  {
   DumpDebug("C");
   DoRSConv(Len2,FFTNum,IP,W);
  }
#ifdef VIRTUAL_CACHE
else /* Virtual Mem based convolution. */
  {FFT_DATA_TYPE *Num1=FFTNum;
   FFT_DATA_TYPE *Num2=FFTCash[Cache].Mem;

   if (Num2==NULL)
     FatalError("Cache %d doesn't exist.\n",Cache);

   DumpDebug("C%d...",Cache);
   DoRConv(Len2,FFTCash[Cache].Mem,FFTNum,IP,W);
  }
#else
#error Ooura3 does not work with disk number caches.
StopTimer(ConvTime);
#endif
}
#endif
#endif

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
rdft(CalcFFTLen(NumLen),1,ddata,IP,W);
StopTimer(FFTTime);
}

void
RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen)
{
StartTimer(FFTTime);
rdft(CalcFFTLen(NumLen),-1,ddata,IP,W);
StopTimer(FFTTime);
}

void
InitFFT(size_t Len)
{
if (W) free(W);
if (IP) free(IP);
Len=CalcFFTLen(Len);
IP = (int *) malloc((3 + (int) sqrt(0.5 * Len)) * sizeof(int));
W  = (double *) malloc(Len / 2 * sizeof(double));

IP[0]=0;
makewt(Len/4,IP,W);
makect(Len/4,IP,W+(Len/4));
}

void
DeInitFFT(void)
{
if (W) free(W);
if (IP) free(IP);
}



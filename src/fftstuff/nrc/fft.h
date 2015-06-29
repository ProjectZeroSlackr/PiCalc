#ifndef FFT_H_
#define FFT_H_ 1

#include "pi.h"

/*
** The maximum number of digits the FFT can handle.
*/
#if FFT_DIGITS==4
#define MaxFFTLen 8388608

#else
#define MaxFFTLen 1073741824U
/*
** Unknown!  I don't have enough memory to check.  However, the FFT trig
** formula says it should, so....  This will also depend on the method
** that you use to generate the trig, and/or whether you use an external
** fft library.  TEST it!!
*/
#endif

/*
** Define a few long double math functions.  If you need them, and
** your compiler supports 'long double', then you should have some
** way of doing thse two things.
*/
#ifdef LONG_DOUBLE_FFT
#error The NRC style FFT does not support long doubles
#endif

#define SINE(A)  sin(A)
#define FLOOR(F) floor(F)

void InitFFT(size_t Len);
void DeInitFFT(void);
#ifndef USE_DEFAULT_CONVOLUTION
void DoConvolution(FFT_DATA_TYPE *FFTNum,int Cache,size_t Len2);
#endif
void FwdTransform(FFT_DATA_TYPE *ddata,BigInt Num, size_t NumLen);
void RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen);

/*
** Macro to convert a number length to the length of the FFT
** This INCLUDES the doubling for the zero padding.
*/
#define CalcFFTLen(zw) ((zw)*2*(RawIntDigits/FFT_DIGITS))

/*
** This converts a FFTLen (from above) to a NumLen.
*/
#define CalcNumLen(zw) (((zw)/2)/(RawIntDigits/FFT_DIGITS))

#if FFT_DIGITS==4
typedef unsigned short int FFT_INT;
#else
typedef unsigned char      FFT_INT;
#endif


#endif


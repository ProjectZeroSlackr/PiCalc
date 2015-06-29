#ifndef FFT_H_
#define FFT_H_ 1

#include "pi.h"

/*
** The maximum number of digits the FFT can handle.
*/
#if FFT_DIGITS==4
#ifdef LONG_DOUBLE_FFT
#define MaxFFTLen 1073741824U
/* Unknown!  */
#else
#define MaxFFTLen 8388608
#endif

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
#ifdef __GNUC__
/*
** DJGPP doesn't have these.  We have to fake it.  With the optimizations
** enabled, DJGPP will use the FPU fsin, so that works okay.
** FLOOR() is a bit of a problem because that's what we use to
** round our long double floating point to a long double integer.  In this
** case, we can use a type conversion to 'long long int'.  It's only used
** in bigmul.c, where we release our carries, so you can change the code
** if you need to.
*/
#define SINE(A)  sin(A)
#define FLOOR(F) ((long long int)(F))
#else
#define SINE(A)  sinl(A)
#define FLOOR(F) floorl(F)
#endif
#else
#define SINE(A)  sin(A)
#define FLOOR(F) floor(F)
#endif

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


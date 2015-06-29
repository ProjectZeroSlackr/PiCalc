#ifndef FFT_H_
#define FFT_H_ 1

#include "pi.h"

/*
** Whether the FFT itself already does the normalization.
#define FFT_IS_NORMALIZED 1
*/

/*
** The maximum number of digits the FFT can handle.
*/
#if FFT_DIGITS==4
#ifdef LONG_DOUBLE_FFT
#define MaxFFTLen 1073741824U
/* Unknown! */
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


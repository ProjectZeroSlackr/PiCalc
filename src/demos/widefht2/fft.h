#ifndef FFT_H_
#define FFT_H_ 1

#include "pi.h"
#include "fixpoint.h"

/*
** Whether the FFT itself already does the normalization.
*/
#define FFT_IS_NORMALIZED 1

#define MaxFFTLen 8388608

#define FFT_DIGITS 32

/* How many 'digits' (elements) to put into each FFT element */
#define RAW_FFT_DIG (FFT_DIGITS/RawIntDigits)

void InitFFT(size_t Len);
void DeInitFFT(void);
void DoConvolution(FFT_DATA_TYPE *FFTNum,int Cache,size_t Len2);
void FwdTransform(FFT_DATA_TYPE *ddata,BigInt Num, size_t NumLen);
void RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen);

/*
** Macro to convert a number length to the length of the FFT
** This INCLUDES the doubling for the zero padding.
#define CalcFFTLen(zw) ((zw)*2*(RawIntDigits/FFT_DIGITS))
*/
#define CalcFFTLen(zw) ((zw)*2/(FFT_DIGITS/RawIntDigits))

/*
** This converts a FFTLen (from above) to a NumLen.
#define CalcNumLen(zw) (((zw)/2)/(RawIntDigits/FFT_DIGITS))
*/
#define CalcNumLen(zw) (((zw)/2)*(FFT_DIGITS/RawIntDigits))


#endif


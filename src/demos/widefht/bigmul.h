#ifndef BIGMUL_H_
#define BIGMUL_H_ 1

#include "pi.h"
#include "fft.h"

size_t CheckFFTMemReq(size_t Len);
void   InitFFTMul(size_t Len);
void   DeInitFFTMul(void);
INT32  FFTMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len, size_t ProdLen, INT32 Scale);

#endif


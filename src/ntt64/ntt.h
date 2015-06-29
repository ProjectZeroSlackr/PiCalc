#ifndef NTT_H_
#define NTT_H_ 1

#include "primes.h"

extern ModInt Prime,PrimvRoot,NthRoot,NthRoot1,MulInv;

void InitFFT(size_t Len);
void DeInitFFT(void);
void DoFwdTransforms(ModInt *NTTData,BigInt Num, size_t NumLen);
void DoRevTransforms(ModInt *NTTData, size_t NumLen);
void SetNTTSize(size_t Len);

#endif


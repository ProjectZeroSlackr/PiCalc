#ifndef NTT_H_
#define NTT_H_ 1

#include "primes.h"

extern ModInt Prime;
extern CmplxModInt PrimvRoot,NthRoot,NthRoot1;
extern ModInt MulInv; /* easier as a single value */

void InitFFT(size_t Len);
void DeInitFFT(void);
void DoFwdTransforms(CmplxModInt *NTTData,BigInt Num, size_t NumLen);
void DoRevTransforms(CmplxModInt *NTTData, size_t NumLen);
void SetFGTSize(size_t Len);

#endif


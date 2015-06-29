#ifndef NTT_H_
#define NTT_H_ 1

#include "primes.h"

/* This should be static, but I need to cheat for Knuth. */
typedef struct {ModInt Prime,PrimvRoot,NthRoot,NthRoot1,MulInv;
                double RecipPrime;} ConstList;
extern ConstList Consts[NPrimes];

extern ModInt Prime;
extern double RecipPrime;
/* not every version needs it, but it needs to be here. */
extern ModInt Inverses[NPrimes][NPrimes];

void InitFFT(size_t Len);
void DeInitFFT(void);
void DoFwdTransforms(ModInt *NTTData,BigInt Num, size_t NumLen,size_t StartP,size_t EndP);
void DoRevTransforms(ModInt *NTTData, size_t NumLen, size_t StartP, size_t EndP);
void SetNTTSize(size_t Len);

#endif


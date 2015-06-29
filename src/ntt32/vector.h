#include "pi.h"
#include "modmath.h"

#ifndef VECTOR_H_
#define VECTOR_H_

void VectorModAdd(ModInt *Num1, ModInt *Num2, size_t Len);
void VectorModSub(ModInt *Num1, ModInt *Num2, size_t Len);
void VectorModMul(ModInt *Num1, ModInt *Num2, size_t Len);
void VectorModMulC(ModInt *Num1, ModInt Num2, size_t Len);
void VectorModButterfly(ModInt *Num1, ModInt *Num2, size_t Len);
void VectorNTT_First2(ModInt *Data, ModInt Trig, size_t Len);
void VectorNTT_Last2(ModInt *Data, ModInt Trig, size_t Len);
void PrepVector(UINT32 Prime, size_t Len);
void InitVectorModMath(size_t Len);
void DeInitVectorModMath(void);

#endif


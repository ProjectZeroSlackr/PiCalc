#ifndef AGM_H_
#define AGM_H_ 1

int  ComputeAGM(size_t Len, size_t MaxPasses);

/* Visible solely to allow us to time these functions for development. */
void Sqrt05(BigInt Root, size_t Len);
void Sqrt20(BigInt Root, size_t Len);
void Sqrt30(BigInt Root, size_t Len);
void Sqrt60(BigInt Root, size_t Len);
void AGMSqrt(BigInt Root, BigInt Num, size_t Len, size_t SubLen);
void AGMDivide(BigInt R, BigInt Num1, BigInt Num2, size_t Len, BigInt Work);

#endif


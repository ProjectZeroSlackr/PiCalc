#ifndef BLOCK_H_
#define BLOCK_H_ 1

INT32  BlockAdd(INT32 *Sum, INT32 *Num1, INT32 *Num2, INT32 Carry, size_t Len);
void   BlockClear(INT32 *Beg,INT32 *End);
void   BlockCopy(INT32 *Dest,INT32 *Src,size_t Len);
size_t BlockCountZeros(INT32 *Num,size_t Len);
INT32  BlockDivBy(INT32 *Result,INT32 *Num,INT32 Val,INT32 Remain,size_t Len);
double BlockMulByFloat(INT32 *Num, double Val, double Carry,size_t Len);
INT32  BlockMulBy(INT32 *Result,INT32 *Num, INT32 Val, INT32 Carry, size_t Len);
INT32  BlockNegate(INT32 *Num, INT32 Borrow,size_t Len,INT32 Val0);
void   BlockPack2(INT32 *Data,size_t Len);
void   BlockPack4(INT32 *Data,size_t Len);
INT32  BlockSlowMul(INT32 *Prod, INT32 *Num1, INT32 *Num2, size_t Len);
INT32  BlockSub(INT32 *Dif, INT32 *Num1, INT32 *Num2, INT32 Borrow, size_t Len);
void   BlockUnpack2(INT32 *Data,size_t Len);
void   BlockUnpack4(INT32 *Data,size_t Len);

#endif


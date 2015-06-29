#ifndef BIGINT_H_
#define BIGINT_H_ 1

#ifdef USE_DISK_NUM
/*typedef size_t BigInt;*/
typedef double BigInt;
#else
typedef INT32* BigInt;
#endif

extern BigInt MulProd;
extern size_t CoreMemAvail;
extern void  *CoreMemPtr;
extern size_t CPU_CACHE;

extern size_t FIXEDBUF_SIZE;
extern void *FixedBuf;

void    InitBigIntPkg(size_t Len);
void    DeInitBigIntPkg(void);
BigInt  CreateBigInt(size_t Len);

/* These 4 functions are the interface between the numbers & the program.*/
void    ReadNumIntoBuf(BigInt Num,INT32 *Buf,size_t Len);
INT32   GetBigIntDigit(BigInt Num,size_t Ndx);
void    SetBigIntDigit(BigInt Num,size_t Ndx,INT32 Val);
void    WriteBufIntoNum(INT32 *Buf,BigInt Num,size_t Len);

INT32   Add(BigInt Sum, BigInt Num1, BigInt Num2, size_t Len);
void    AddInt(BigInt Num,INT32 Val);
void    ClearBigInt(BigInt Num, size_t Len);
void    Copy(BigInt Num1, BigInt Num2, size_t Len);
void    CopyProd(BigInt Result, BigInt Prod, size_t Len);
size_t  CountZeros(BigInt Num, size_t *Left,size_t *Right,size_t Len);
void    DivBy(BigInt Result,BigInt Num,INT32 Val, size_t Len);
size_t  FindFirstNonZero(BigInt Num, size_t Len);
void    FullMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len);
INT32   HalfDiff(BigInt Dif, BigInt Num1, BigInt Num2, size_t Len);
void    HalfMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len);
void    HalfMul2(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len);
int     IsZero(BigInt Num, size_t Len);
void    Mul(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len);
INT32   MulBy(BigInt Result,BigInt Num, INT32 Val, size_t Len);
double  MulByFloat(BigInt Num, double Val, size_t Len);
void    N1R0Mul(BigInt Prod, BigInt Num1, BigInt Num2, BigInt Work,size_t Len);
void    Negate(BigInt Num, size_t Len);
int     NumIs(BigInt Num, INT32 Val1, INT32 Val2);
INT32   RevSubInt(INT32 Val,BigInt Num, size_t Len);
INT32   RippleAdd(BigInt Limit,BigInt Sum,BigInt Num1,BigInt Num2,size_t Len);
INT32   RippleSub(BigInt Limit,BigInt Dif,BigInt Num1,BigInt Num2,size_t Len);
void    SetNum(BigInt Num,size_t Len, INT32 Val1, INT32 Val2);
int     SpecialAGMFunc1(BigInt AGM_B2,BigInt AGM_A2,BigInt AGM_C2, size_t Len);
void    SpecialFullMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len, BigInt Work);
void    SpecialSquare(BigInt Prod, BigInt Num, size_t Len, BigInt Work);
INT32   Sub(BigInt Dif, BigInt Num1, BigInt Num2, size_t Len);
void    SubInt(BigInt Num, INT32 Val);

/* I/O functions */
void    DumpBigInt(char *Str, BigInt Num, size_t Len);
void    PrintFormattedPi(char *Str, double ETime, BigInt Num, size_t Len);
int     LoadData(BigInt Var1, BigInt Var2, BigInt Var3, BigInt Var4,
                 BigInt Var5, BigInt Var6,
                 clock_t *StartTime, double *Pow2, size_t *Pass,
                 size_t Len, size_t Algor);
void    SaveData(BigInt Var1, BigInt Var2, BigInt Var3, BigInt Var4,
                 BigInt Var5, BigInt Var6,
                 double StartTime, double Pow2, size_t Pass,
                 size_t Len,size_t Algor);

void    DeleteSaveFile(void);

#endif



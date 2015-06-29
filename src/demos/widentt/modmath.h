#ifndef MODMATH_H_
#define MODMATH_H_ 1

#include "sys_cfg.h"

void ModMul(ModInt *Prod, ModInt *Ap, ModInt *Bp);
void ModAdd(ModInt *Sum, ModInt *a,ModInt *b);
void ModSub(ModInt *Dif,ModInt *a,ModInt *b);
void ModAddSub(ModInt *a, ModInt *b);
void ModClear(ModInt *a);
void ModSet1(ModInt *a);
void ModSet(ModInt *a, UINT32 z);
void ModSub1(ModInt *a);
void ModDivInt(ModInt *a,size_t D);
void HexStr2Mod(ModInt *a, char *b);
void ModIPow(ModInt *Prod,ModInt *Base,size_t Expon);
void ModPow(ModInt *Prod,ModInt *Base,ModInt *Expon);
void FindInverse(ModInt *Inv,ModInt *Num);
void ModCopy(ModInt *Dest,ModInt *Src);

void CB_DivInt(ModInt *a,UINT32 d,ModInt *q,UINT32 *r); // d==16 bits
UINT32 CB_Add(ModInt *s, ModInt *a, ModInt *b);
void CB_MulInt(ModInt *m,UINT32 i);                     // i==16 bits
void CB_AddInt(ModInt *s,UINT32 v);
void CB_Dump(ModInt *p);

#endif


#ifndef MODMATH_H_
#define MODMATH_H_ 1

#include "sys_cfg.h"

UINT32 RAW_Add(ModInt *s, ModInt *a, ModInt *b);
UINT32 RAW_AddInt(ModInt *s,UINT32 v);
UINT32 RAW_Sub(ModInt *d,ModInt *a, ModInt *b);
UINT32 RAW_SubInt(ModInt *d,UINT32 v);
UINT32 RAW_MulInt(ModInt *s,UINT32 i);
UINT32 RAW_DivInt(ModInt *a,UINT32 d,ModInt *q);
#define RAW_Dump(gq) printf("%08x %08x %08x\n",(gq)->Hi,(gq)->Med,(gq)->Low);

void ModClear(ModInt *a);
void ModSet1(ModInt *a);
void ModSet(ModInt *a, UINT32 z);
void ModMul(ModInt *Prod, ModInt *a, ModInt *b);
void ModAdd(ModInt *Sum, ModInt *a,ModInt *b);
void ModSub(ModInt *Dif,ModInt *a,ModInt *b);
void ModAddSub(ModInt *a, ModInt *b);
void ModSubInt(ModInt *a,UINT32 v);
void ModAddInt(ModInt *a,UINT32 v);
void IModPow(ModInt *Prod,UINT32 base,UINT32 expon);
void ModPow(ModInt *Prod,ModInt *Base,ModInt *Expon);
void HexStr2Mod(ModInt *a, char *b);
int  IsModZero(ModInt *Num);


void CmplxModMul(CmplxModInt *Prod, CmplxModInt *a, CmplxModInt *b);
void CmplxModAdd(CmplxModInt *Sum, CmplxModInt *a,CmplxModInt *b);
void CmplxModSub(CmplxModInt *Dif,CmplxModInt *a,CmplxModInt *b);
void CmplxModAddSub(CmplxModInt *a, CmplxModInt *b);
void CmplxModSet(CmplxModInt *a, UINT32 r, UINT32 i);
void CmplxModIPow(CmplxModInt *Prod,CmplxModInt *Base,size_t Expon);
void CmplxModPow(CmplxModInt *Prod,CmplxModInt *Base,ModInt *Expon);
void CmplxConj(CmplxModInt *N);


#endif


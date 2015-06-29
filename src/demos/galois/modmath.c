#include "pi.h"
#include "modmath.h"
#include "fgt.h"

/*
** Don't forget to do something about the FGT having so many ModInt
** stack vars.
*/

/* Begin low level data struct manipulations */

UINT32
RAW_Add(ModInt *s, ModInt *a, ModInt *b)
{UINT64 c=0;
c=(UINT64)a->Low+(UINT64)b->Low;
s->Low=c&0xffffffff;c>>=32;
c=(UINT64)a->Med+(UINT64)b->Med+c;
s->Med=c&0xffffffff;c>>=32;
c=(UINT64)a->Hi+(UINT64)b->Hi+c;
s->Hi=c&0xffffffff;c>>=32;
return (UINT32)(c&0xffffffff);
}

UINT32
RAW_AddInt(ModInt *s,UINT32 v)
{UINT64 c=0;
c=(UINT64)s->Low+(UINT64)v;
s->Low=c&0xffffffff;c>>=32;
c=(UINT64)s->Med+c;
s->Med=c&0xffffffff;c>>=32;
c=(UINT64)s->Hi+c;
s->Hi=c&0xffffffff;c>>=32;
return (UINT32)(c&0xffffffff);
}

UINT32
RAW_Sub(ModInt *d,ModInt *a, ModInt *b)
{UINT64 c=0;
c=(UINT64)a->Low-(UINT64)b->Low;
d->Low=c&0xffffffff;c>>=63;
c=(UINT64)a->Med-(UINT64)b->Med-c;
d->Med=c&0xffffffff;c>>=63;
c=(UINT64)a->Hi-(UINT64)b->Hi-c;
d->Hi=c&0xffffffff;c>>=63;
return c;
}

UINT32
RAW_SubInt(ModInt *d,UINT32 v)
{UINT64 c=0;
c=(UINT64)d->Low-(UINT64)v;
d->Low=c&0xffffffff;c>>=63;
c=(UINT64)d->Med-c;
d->Med=c&0xffffffff;c>>=63;
c=(UINT64)d->Hi-c;
d->Hi=c&0xffffffff;c>>=63;
return c;
}

UINT32
RAW_MulInt(ModInt *s,UINT32 i)
{UINT64 p;
p=((UINT64)s->Low)*i;
s->Low=p&0xffffffff;p>>=32;
p=((UINT64)s->Med)*i+p;
s->Med=p&0xffffffff;p>>=32;
p=((UINT64)s->Hi)*i+p;
s->Hi=p&0xffffffff;p>>=32;
return (UINT32)p;
}

UINT32
RAW_DivInt(ModInt *a,UINT32 d,ModInt *q)
{UINT64 Z=0;

Z=a->Hi;  q->Hi =Z / d;Z=Z % d;Z<<=32;
Z+=a->Med;q->Med=Z / d;Z=Z % d;Z<<=32;
Z+=a->Low;q->Low=Z / d;Z=Z % d;
return (UINT32)Z;
}

static int
RAW_ShiftR(ModInt *a)
{int Carry;
Carry=a->Low & 1;
a->Low>>=1;if (a->Med & 1) a->Low+=(1<<31);
a->Med>>=1;if (a->Hi  & 1) a->Med+=(1<<31);
a->Hi >>=1;
return Carry;
}

static int
RAW_ShiftL(ModInt *a)
{int Carry=0;
if (a->Hi & 0x80000000) Carry=1;
a->Hi <<=1;if (a->Med & 0x80000000) a->Hi++;
a->Med<<=1;if (a->Low & 0x80000000) a->Med++;
a->Low<<=1;
return Carry;
}


/* Begin ModInt functions */
void
ModClear(ModInt *a)
{
a->Hi=0;
a->Med=0;
a->Low=0;
}

void
ModSet1(ModInt *a)
{
a->Hi=0;
a->Med=0;
a->Low=1;
}

void
ModSet(ModInt *a, UINT32 z)
{
a->Hi=0;
a->Med=0;
a->Low=z;
}

void
ModMul(ModInt *Prod,ModInt *a, ModInt *b)
{
  int i = 96;
  ModInt c,aa;
  ModSet(&c,0);
  aa=*a;
  while (i--)
    {
      RAW_ShiftL(&c);
      if (RAW_Sub(&c, &c, &Prime)) RAW_Add(&c, &c, &Prime);
      if (RAW_ShiftL(&aa))         RAW_Add(&c, &c, b);
      if (RAW_Sub(&c, &c, &Prime)) RAW_Add(&c, &c, &Prime);
    }
  *Prod=c;
}


void
ModAdd(ModInt *Sum, ModInt *a,ModInt *b)
{
if (RAW_Add(Sum,a,b)) RAW_Sub(Sum,Sum,&Prime);
else if (RAW_Sub(Sum,Sum,&Prime)) RAW_Add(Sum,Sum,&Prime);
}

void
ModSub(ModInt *Dif,ModInt *a,ModInt *b)
{
if (RAW_Sub(Dif,a,b)) RAW_Add(Dif,Dif,&Prime);
}

void
ModAddSub(ModInt *a, ModInt *b)
{ModInt AA=*a;
 ModInt BB=*b;
 ModAdd(a,&AA,&BB);
 ModSub(b,&AA,&BB);
}

void
ModSubInt(ModInt *a,UINT32 v)
{ModInt Temp;
ModSet(&Temp,v);
ModSub(a,a,&Temp);
}

void
ModAddInt(ModInt *a,UINT32 v)
{ModInt Temp;
ModSet(&Temp,v);
ModAdd(a,a,&Temp);
}

void
IModPow(ModInt *Prod,UINT32 base,UINT32 expon)
/* Raise an integer to an integer power, returning a ModInt */
{ModInt b;

if (expon==0) {ModSet(Prod,1);return;}

ModSet(&b,base);
while (!(expon&1))
  {
   ModMul(&b,&b,&b);
   expon>>=1;
  }
*Prod=b;

while (expon>>=1)
  {
   ModMul(&b,&b,&b);
   if (expon&1) ModMul(Prod,Prod,&b);
  }
}

static int
ModTestLowBit(ModInt *Expon)
{
if (Expon->Low & 1) return 1;
return 0;
}

int
IsModZero(ModInt *Num)
{
if (Num->Low) return 0;
if (Num->Med) return 0;
if (Num->Hi)  return 0;
return 1;
}

void
ModPow(ModInt *Prod,ModInt *Base,ModInt *Expon)
{ModInt b;

b=*Base;
while (!ModTestLowBit(Expon)) {ModMul(&b,&b,&b);RAW_ShiftR(Expon);}
*Prod=b;
RAW_ShiftR(Expon);

while (!IsModZero(Expon))
  {
   ModMul(&b,&b,&b);
   if (ModTestLowBit(Expon)) ModMul(Prod,Prod,&b);
   RAW_ShiftR(Expon);
  }
}

void
HexStr2Mod(ModInt *a, char *b)
/* Convert a hex string to a modular number.  No error checking */
{
  ModClear(a);
  while (*b)
    {
      RAW_ShiftL(a);
      RAW_ShiftL(a);
      RAW_ShiftL(a);
      RAW_ShiftL(a);
      a->Low += tolower(*b) < 'a' ? *b - '0' : tolower(*b) - ('a' - 10);
      b++;
    }
}




/* Begin the complex ModInt functions */
void
CmplxModMul(CmplxModInt *Prod, CmplxModInt *a, CmplxModInt *b)
{CmplxModInt P;
 ModInt T1,T2;
ModMul(&T1,&a->r,&b->r);
ModMul(&T2,&a->i,&b->i);
ModSub(&P.r,&T1,&T2);
ModMul(&T1,&a->r,&b->i);
ModMul(&T2,&a->i,&b->r);
ModAdd(&P.i,&T1,&T2);
*Prod=P;
}

void
CmplxModAdd(CmplxModInt *Sum, CmplxModInt *a,CmplxModInt *b)
{
ModAdd(&Sum->r,&a->r,&b->r);
ModAdd(&Sum->i,&a->i,&b->i);
}

void
CmplxModSub(CmplxModInt *Dif,CmplxModInt *a,CmplxModInt *b)
{
ModSub(&Dif->r,&a->r,&b->r);
ModSub(&Dif->i,&a->i,&b->i);
}

void
CmplxModAddSub(CmplxModInt *a, CmplxModInt *b)
{CmplxModInt Dif,Sum;
CmplxModSub(&Dif,a,b);
CmplxModAdd(&Sum,a,b);
*a=Sum;
*b=Dif;
}

void
CmplxModSet(CmplxModInt *a, UINT32 r, UINT32 i)
{
ModSet(&a->r,r);
ModSet(&a->i,i);
}

void
CmplxConj(CmplxModInt *N)
{
if (!IsModZero(&N->i)) ModSub(&N->i,&Prime,&N->i);
}

void
CmplxModIPow(CmplxModInt *Prod,CmplxModInt *Base,size_t Expon)
/* raise a complex to an integer power, returning a complex */
{CmplxModInt b;

if (Expon==0) {CmplxModSet(Prod,1,0);return;}

b=*Base;
while (!(Expon&1)) {CmplxModMul(&b,&b,&b);Expon>>=1;}
*Prod=b;

while (Expon>>=1)
  {
   CmplxModMul(&b,&b,&b);
   if (Expon&1) CmplxModMul(Prod,Prod,&b);
  }
}

void
CmplxModPow(CmplxModInt *Prod,CmplxModInt *Base,ModInt *Expon)
/* raise a complex to an modint power, returning a complex */
{CmplxModInt b;

if (IsModZero(Expon)) {CmplxModSet(Prod,1,0);return;}

b=*Base;
while (!ModTestLowBit(Expon)) {CmplxModMul(&b,&b,&b);RAW_ShiftR(Expon);}
*Prod=b;

RAW_ShiftR(Expon);
while (!IsModZero(Expon))
  {
   CmplxModMul(&b,&b,&b);
   if (ModTestLowBit(Expon)) CmplxModMul(Prod,Prod,&b);
   RAW_ShiftR(Expon);
  }
}




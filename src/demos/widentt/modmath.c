/*
** This current setup takes 18,780 ticks to do a "-1m" test.
** This is about 40 times slower than normal.  I don't do
** any caching of the NTT data.  That would cut it by 1/3rd, to
** about 26 times.  The modulo is done niavely, with a mul by
** recip, then mul by prime.  If it was done smarter, with only
** one mul, that would be about 8 times.  It's not unreasonable
** to figure that an assembly (or better coded) FFT/FHT mul would
** be twice as fast, so that's about 4 times.
*/
#include "pi.h"
#include "modmath.h"
#include "ntt.h"

/*
** Don't forget to do something about the NTT having so many ModInt
** stack vars.
*/

#define Mod_BITS 16

/* Begin low level data struct manipulations */

void
CB_Dump0(UINT16 *p,size_t Len)
{int x;
for (x=0;x<Len;x++) printf("%04x",p[x]);
printf("\n");
}

void
CB_Dump(ModInt *p)
{
CB_Dump0(&p->Data[0],Mod_SIZE);
}

UINT32
CB_Add(ModInt *s, ModInt *a, ModInt *b)
{UINT32 c = 0;
 int x;

 for (x=Mod_SIZE-1;x>=0;x--)
    {
      c          = ((UINT32)b->Data[x]) + ((UINT32)a->Data[x]) + c;
      s->Data[x] = c;
      c        >>= 16;
    }
 return c;
}

UINT32
CB_Sub(ModInt *d,ModInt *a, ModInt *b)
{UINT32 c = 0;
 int x;
 for (x=Mod_SIZE-1;x>=0;x--)
   {
    c          = ((UINT32)a->Data[x]) - ((UINT32)b->Data[x]) - c;
    d->Data[x] = c;
    c          = (c >> 16) & 1;
   }
 return c;
}

void
CB_DivInt(ModInt *a,UINT32 d,ModInt *q,UINT32 *r)
{int x;
 UINT64 Z=0;

for (x=0;x<Mod_SIZE;x++)
  {
   Z=(Z*65536)+a->Data[x];
   q->Data[x]=Z / d;
   Z         =Z % d;
  }
*r=Z;
}

void
CB_MulInt(ModInt *s,UINT32 i)
{int x;
 UINT64 Z=0;
for (x=Mod_SIZE-1;x>=0;x--)
  {
   Z=((UINT64)s->Data[x])*i+Z;
   s->Data[x]=Z % 65536;
   Z         =Z / 65536;
  }
}

void
CB_AddInt(ModInt *s,UINT32 v)
{
 UINT32 c = v;
 int i;
 for (i=Mod_SIZE-1;i>=0;i--)
   {
    c = c + s->Data[i];
    s->Data[i] = c;
    c >>= 16;
   }
}

static UINT32
CB_ShiftL(ModInt *a)
{
  UINT32 c = 0;
  int i;
  for (i=Mod_SIZE-1;i>=0;i--)
    {
      c |= ( ((UINT32) a->Data[i]) << 1);
      a->Data[i] = c;
      c = (c >> 16) & 1;
    }
  return c;
}

static UINT32
CB_ShiftR(ModInt *a)
{
  UINT32 c = 0;
  int i;
  for (i=0;i<Mod_SIZE;i++)
    {
      c |= a->Data[i];
      a->Data[i] = c >> 1;
      c = (c & 1) << 16;
    }
  return c;
}
/* End low level data struct manipulations */

/* Begin low level modular number manipulations */
void
ModClear(ModInt *a)
{int x;
for (x=0;x<Mod_SIZE;x++) a->Data[x]=0;
}

void
ModCopy(ModInt *Dest,ModInt *Src)
{
*Dest=*Src;
}

void
ModSet1(ModInt *a)
{
ModClear(a);
a->Data[Mod_SIZE-1]=1;
}

void
ModSet(ModInt *a,UINT32 z)
{
ModClear(a);
a->Data[Mod_SIZE-1]=z % 65536;z/=65536;
a->Data[Mod_SIZE-2]=z;
}

static int
ModTest(ModInt *a)
{
  int i = Mod_SIZE;
  while (i--)
    if (a->Data[i])
      return 1;
  return 0;
}

static int
ModTestLowBit(ModInt *a)
{
if (a->Data[Mod_SIZE-1] & 1) return 1;
return 0;
}
/* End low level modular number manipulations */


/* ModMul(a,b,m) -- prod = a * b mod Prime, modular multiply */
void
ModMul(ModInt *Prod,ModInt *a, ModInt *b)
{
  UINT16 i = Mod_SIZE * Mod_BITS;
  ModInt c,aa;
  ModClear(&c);
  ModCopy(&aa,a);
  while (i--)
    {
      CB_ShiftL(&c);
      if (CB_Sub(&c, &c, &Prime)) CB_Add(&c, &c, &Prime);
      if (CB_ShiftL(&aa))     CB_Add(&c, &c, b);
      if (CB_Sub(&c, &c, &Prime)) CB_Add(&c, &c, &Prime);
    }
  ModCopy(Prod, &c);
}

void
ModAdd(ModInt *Sum, ModInt *a,ModInt *b)
{
if (CB_Add(Sum,a,b)) CB_Sub(Sum,Sum,&Prime);
else if (CB_Sub(Sum,Sum,&Prime)) CB_Add(Sum,Sum,&Prime);
}

void
ModSub(ModInt *Dif,ModInt *a,ModInt *b)
{
if (CB_Sub(Dif,a,b)) CB_Add(Dif,Dif,&Prime);
}

void
ModAddSub(ModInt *a, ModInt *b)
{ModInt AA=*a;
 ModInt BB=*b;
 ModAdd(a,&AA,&BB);
 ModSub(b,&AA,&BB);
}

void
ModSub1(ModInt *a)
{ModInt One;
ModSet1(&One);
ModSub(a,a,&One);
}

void
ModDivInt(ModInt *a,size_t D)
{int z;
CB_DivInt(a,D,a,&z);
}

void
HexStr2Mod(ModInt *a, char *b)
/* Convert a hex string to a modular number.  No error checking */
{
  ModClear(a);
  while (*b)
    {
      CB_ShiftL(a);
      CB_ShiftL(a);
      CB_ShiftL(a);
      CB_ShiftL(a);
      a->Data[Mod_SIZE-1] += tolower(*b) < 'a'
        ? *b - '0' : tolower(*b) - ('a' - 10);
      b++;
    }
}

void
ModIPow(ModInt *Prod,ModInt *Base,size_t Expon)
/* to an integer modular power */
{ModInt b;

b=*Base;
while (!(Expon&1)) {ModMul(&b,&b,&b);Expon>>=1;}
*Prod=b;

while (Expon>>=1)
  {
   ModMul(&b,&b,&b);
   if (Expon&1) ModMul(Prod,Prod,&b);
  }
}

void
ModPow(ModInt *Prod,ModInt *Base,ModInt *Expon)
{ModInt b;

b=*Base;
while (!ModTestLowBit(Expon)) {ModMul(&b,&b,&b);CB_ShiftR(Expon);}
*Prod=b;
CB_ShiftR(Expon);

while (ModTest(Expon))
  {
   ModMul(&b,&b,&b);
   if (ModTestLowBit(Expon)) ModMul(Prod,Prod,&b);
   CB_ShiftR(Expon);
  }
}




#ifndef CRT_H_
#define CRT_H_ 1

/*
** Some of the routines in here use the pentium FPU to do some
** of the work.  This faster than other generic ways, but it's still
** slower than it should be, because Intel'x x87 FPU can only operate
** on signed integers.  And we are dealing with unsigned.
*/

/*
** It is fairly important that these routines be efficient.  That means
** assembly rather than C.  For the gcc C versions, see the gcc486/586
** directories.  Realistically, the only reason they are even being
** provided in C is so that you can at least get the program running.
** For decent performance, you really need to convert these to whatever
** form of assembly your compiler supports.
**
** The multiply works tolerably well.  Not great, but tolerable.  But
** should really write the StripCRT() as asm.  That saps a lot of speed.
** And if you are going to do that, you might as well convert the rest
** of thse too, and provide non-macro versions for the MulCRT?().
**
** If your compiler has a 64 bit 'long long' you can get by for a while
** by providing 'long long' versions like I do in the 'longlong' directory.
** But asm would be better.
*/

/*
** Not to be changed!  This just provides a symbolic description of
** hardwired stuff.
*/
#define CRT_LEN NPrimes

typedef UINT32 CRTNum;

#if 0
static INLINE void
Mul323264(UINT32 Num1, UINT32 Num2, UINT32 *Hi, UINT32 *Lo)
{UINT32 H,L;
  asm("movl %2,%%eax
        mull %3
        movl %%edx,%0
        movl %%eax,%1
      "
    : "=g" (H), "=g" (L)
    : "g" (Num1) , "g" (Num2)
    : "%eax" , "%edx" , "cc"
    );
*Hi=H;*Lo=L;
}
#else

#if 1
static INLINE void
Mul323264(UINT32 Num1, UINT32 Num2, UINT32 *Hi, UINT32 *Lo)
/*
** Even though the name says we are multiply two 32 bit numbers,
** for the 586 version, that second number will actually only be
** 31 bits (because it's a ModInt).  That meas we can use the
** floating point unit to do it, since we have the 64th bit
** to allow us to justify it.
**
** This is _still_ not extremely efficient because the variables
** (and the first number) are 32 bit unsigned and the x87 doesn't
** have any way to load unsigned ints.  As usual, Intel screwed up.
*/
{union {long double D;unsigned long int i[3];} ChopShop;
ChopShop.D=((long double)Num1)*((long double)Num2);
ChopShop.D+=(2147483648.0*2147483648.0*2.0); /* justify it at the 64th bit*/
*Lo=ChopShop.i[0];
*Hi=ChopShop.i[1] & ~2147483648U; /* remove the justify bit */
}
#endif

#if 0
// Slow 274.  Mostly for example.
static INLINE void

Mul323264(UINT32 Num1, UINT32 Num2, UINT32 *Hi, UINT32 *Lo)
{
 UINT32 H1,L1;
 UINT32 H2,L2;
 UINT32 LL,LH,HL,HH;

 L1=Num1 & 0xffff;
 H1=(Num1 >> 16) & 0xffff;

 L2=Num2 & 0xffff;
 H2=(Num2 >> 16) & 0xffff;

 LL=L2*L1;LH=L2*H1;
 HL=H2*L1;HH=H2*H1;

 L2=LH+HL;if (L2 < LH) HH+=65536;
 H2=HH+(L2 >> 16);
 L2<<=16;
 L2=LL+L2;if (L2 < LL) H2++;
 *Hi=H2;
 *Lo=L2;
}
#endif
#endif

static void
MulCRT(CRTNum *Num, ModInt Val,int Len)
/* Speed isn't too important in this function. */
{CRTNum Carry;
 CRTNum H,L;

Carry=0;
Num+=CRT_LEN;

while (Len--)
  {
   Num--;
   Mul323264(*Num,Val,&H,&L);
   L=Carry+L;if (L < Carry) H++;

   *Num=L;
   Carry=H;
  }
}

/*
** It might be better to do the seperate routines, rather than these macros.
** But, for simplicity of this file, they're being done as macros.
*/
#define MulCRT1(r,n,v) MulCRTz((r),(n),(v),1)
#define MulCRT2(r,n,v) MulCRTz((r),(n),(v),2)
#define MulCRT3(r,n,v) MulCRTz((r),(n),(v),3)
#define MulCRT4(r,n,v) MulCRTz((r),(n),(v),4)
#define MulCRT5(r,n,v) MulCRTz((r),(n),(v),5)
#define MulCRT6(r,n,v) MulCRTz((r),(n),(v),6)
#define MulCRT7(r,n,v) MulCRTz((r),(n),(v),7)

#if 0
static INLINE UINT32
Mulz2(UINT32 a, UINT32 b, UINT32 c, UINT32 *Dest)
{UINT32 r,z;
 z=*Dest;
  asm("movl %2,%%eax
        mull %3
        addl %4,%%eax
        adcl $0,%%edx
        addl %0,%%eax
        adcl $0,%%edx
        movl %%eax,%0
        movl %%edx,%1
      "
    : "=g" (z), "=g" (r)
    : "g" (a) , "g" (b) , "g" (c)
    : "%eax" , "%edx" , "cc"
    );
*Dest=z;
return r;
}
#else

static INLINE UINT32
Mulz2(UINT32 a, UINT32 b, UINT32 Carry, UINT32 *Dest)
{UINT32 H,L;

 Mul323264(a,b,&H,&L);
 L+=Carry;if (L < Carry) H++;

 L+=*Dest;if (L < *Dest) H++;

 *Dest=L;
 return H;
}
#endif

static INLINE void
MulCRTz(UINT32 *Res,UINT32 *Num, UINT32 Val,int Len)
{UINT32 Carry;

Carry=0;
Num+=CRT_LEN;Res+=CRT_LEN;
while (Len--)
  {
   --Num;--Res;
   Carry=Mulz2(*Num,Val,Carry,Res);
  }
*(--Res)=Carry;
}

static void
Add2CRT(CRTNum *Num1, CRTNum *Num2)
{CRTNum S,C;
 int x;

C=0;
Num1+=NPrimes-1;
Num2+=NPrimes-1;
for (x=NPrimes;x > 0;x--)
  {
   S = *Num1 + C;
   if (S < C) S=*Num2;
   else
     {
      S += *Num2;
      C=!(!(S < *Num2));
     }
   *Num1=S;
   Num1--;Num2--;
  }
}

static UINT32
StripCRT(CRTNum *Num)
{int x;
 long double D;  /* NOTE that this requires 64 bit long doubles. */
 UINT32 z;

D=0.0;
for (x=0;x<NPrimes;x++)
  {
   D=(D*4294967296.0)+Num[x];
   z=D / 100000000.0;
   Num[x] = z;
   D = D - z * 100000000.0;
  }
return (INT32)D;
}

/*
** Things that don't need to be changed.  They just need to be here,
** is all.
*/
static void
ClearCRT(CRTNum *Num)
{size_t x;
for (x=0;x<CRT_LEN;x++) Num[x]=0;
}

static void
SetCRT(CRTNum *Num, UINT32 Val)
{
ClearCRT(Num);
Num[CRT_LEN-1]=Val;
}

static void
DumpCRT(CRTNum *Num)
{size_t x;
for (x=0;x<CRT_LEN;x++) printf("%08x",Num[x]);
printf("\n");
}

static void
CopyCRT(CRTNum *Dest, CRTNum *Src)
{size_t x;
for (x=0;x<CRT_LEN;x++) Dest[x]=Src[x];
}

#endif


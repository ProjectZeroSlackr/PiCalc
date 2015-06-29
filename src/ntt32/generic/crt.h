/*
** This file contains very low level routines required for the Chinese
** Remainder Theorem.  Due to performance reasons, these routines have
** to be written very carefully, and possibly use assembly language or
** compiler specific code.  The should also be 'inline' code because
** a subroutine call would probably cost as much as the operation itself.
**
** All of these routines are provided in generic versions so you can
** at least run the code (SLOWLY!) while you are writing the higher
** performance compiler/processor specific code. (such as asm.)
**
** You might also want to check Mikko Tommila's APFloat package for
** asm code to do similar things for a variety of processors.
**
** If you need to write these in asm yourself, I am using the normal
** 'Big Endian' format (where the MSD is index 0).  You could change
** that, but its doubtful you'd need to.
**
** Also, I'd like to ask that if you write any custom  routines  for
** your processor or compiler, please send a copy to  me  so  I  can
** include them in any future  versions.  Your doing this might save
** somebody less skilled than you a lot of frustration.
**
*/

#ifndef CRT_H_
#define CRT_H_ 1

/*
** Not to be changed!  This just provides a symbolic description of
** hardwired stuff.  The CRT_LEN is going to be the same as the number
** of primes we are using, which happens to be 8.
*/
#define CRT_LEN NPrimes

typedef UINT32 CRTNum;

/************************************************
********  GENERIC 'ANSI/ISO C' VERSIONS  ********
************************************************/

static INLINE void
Mul323264(UINT32 Num1, UINT32 Num2, UINT32 *Hi, UINT32 *Lo)
/*
** Don't have 64 bit ints (or 64 bit long double), so we have to break
** the 32x32=64 mul into 4 16x16=32 muls.
*/
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

static INLINE UINT32
Mulz2(UINT32 a, UINT32 b, UINT32 Carry, UINT32 *Dest)
{UINT32 H,L;

 Mul323264(a,b,&H,&L);
 L+=Carry;if (L < Carry) H++;

 L+=*Dest;if (L < *Dest) H++;

 *Dest=L;
 return H;
}


#ifdef UNROLLED_CRT

static INLINE void
MulCRT1(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[CRT_LEN-1],Val,Carry,&Res[CRT_LEN-1]);
Res[CRT_LEN-2]=Carry;
}

static INLINE void
MulCRT2(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[CRT_LEN-1],Val,Carry,&Res[CRT_LEN-1]);
Carry=Mulz2(Num[CRT_LEN-2],Val,Carry,&Res[CRT_LEN-2]);
Res[CRT_LEN-3]=Carry;
}

static INLINE void
MulCRT3(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[CRT_LEN-1],Val,Carry,&Res[CRT_LEN-1]);
Carry=Mulz2(Num[CRT_LEN-2],Val,Carry,&Res[CRT_LEN-2]);
Carry=Mulz2(Num[CRT_LEN-3],Val,Carry,&Res[CRT_LEN-3]);
Res[CRT_LEN-4]=Carry;
}

#if NPrimes==8
static INLINE void
MulCRT4(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[7],Val,Carry,&Res[7]);
Carry=Mulz2(Num[6],Val,Carry,&Res[6]);
Carry=Mulz2(Num[5],Val,Carry,&Res[5]);
Carry=Mulz2(Num[4],Val,Carry,&Res[4]);
Res[3]=Carry;
}

static INLINE void
MulCRT5(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[7],Val,Carry,&Res[7]);
Carry=Mulz2(Num[6],Val,Carry,&Res[6]);
Carry=Mulz2(Num[5],Val,Carry,&Res[5]);
Carry=Mulz2(Num[4],Val,Carry,&Res[4]);
Carry=Mulz2(Num[3],Val,Carry,&Res[3]);
Res[2]=Carry;
}

static INLINE void
MulCRT6(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[7],Val,Carry,&Res[7]);
Carry=Mulz2(Num[6],Val,Carry,&Res[6]);
Carry=Mulz2(Num[5],Val,Carry,&Res[5]);
Carry=Mulz2(Num[4],Val,Carry,&Res[4]);
Carry=Mulz2(Num[3],Val,Carry,&Res[3]);
Carry=Mulz2(Num[2],Val,Carry,&Res[2]);
Res[1]=Carry;
}

static INLINE void
MulCRT7(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[7],Val,Carry,&Res[7]);
Carry=Mulz2(Num[6],Val,Carry,&Res[6]);
Carry=Mulz2(Num[5],Val,Carry,&Res[5]);
Carry=Mulz2(Num[4],Val,Carry,&Res[4]);
Carry=Mulz2(Num[3],Val,Carry,&Res[3]);
Carry=Mulz2(Num[2],Val,Carry,&Res[2]);
Carry=Mulz2(Num[1],Val,Carry,&Res[1]);
Res[0]=Carry;
}
#endif

#else

static void
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
#endif


#if 0
static void
Add2CRT(CRTNum *Num1, CRTNum *Num2)
{CRTNum S,C;
 int x;

C=0;
Num1+=CRT_LEN-1;
Num2+=CRT_LEN-1;
for (x=CRT_LEN;x > 0;x--)
  {
   S = *Num1 + C;
   if (S < C) S=*Num2;
   else
     {
      S += *Num2;
/*    C=(S < *Num2)?1:0; */
      C=!(!(S < *Num2));
/*
** I guess I better explain that.... It's the same as: C=(S < *Num2)?1:0;
** ANSI/ISO mandates that a 'true' condition will evaluate to a value of
** '1', and a false condition will evaluate to '0'.  That happens to be
** just what I need.  However, some compilers might cut corners and
** just provide a non-zero value for true.  So, I do two '!' to make sure
** that it's right.  The compiler will strip those out if they aren't
** needed.  It's possible that a compiler will ignore those too and still
** just provide a non-zero for true, but if that's the case, then your
** compiler is probably rather unreliable to begin with.  Oh, by the way,
** that expression is usually fairly quick.  It's often just a single
** comparasion and then some cpu instruction to transfer a flag to a var.
** The ?: (or other) version would most likely require a jump instruction.
*/
     }
   *Num1=S;
   Num1--;Num2--;
  }
}

#else

static void
Add2CRT(CRTNum *Num1, CRTNum *Num2)
{double S,C;
 int x;

C=0;
for (x=CRT_LEN-1;x>=0;x--)
  {
   S=((double)Num1[x])+((double)Num2[x])+C;
   C=0.0;
   if (S >= 4294967296.0) {S -= 4294967296.0; C=1.0;}
   Num1[x]=S;
  }
}
#endif

static UINT32
StripCRT(CRTNum *Num)
/*
** Strip out 1e8 digits.  Unfortunately, a 'double' isn't big enough
** to do it all in one chunk!
*/
{double D1,D2;int x;
D1=D2=0.0;
for (x=0;x<CRT_LEN;x++)
  {
   D1=D1*4294967296.0;
   D1=D1 + Num[x];
   Num[x] = D1 / 10000.0;
   D1     = D1 - Num[x] * 10000.0;
  }
for (x=0;x<CRT_LEN;x++)
  {
   D2=D2*4294967296.0;
   D2=D2 + Num[x];
   Num[x] = D2 / 10000.0;
   D2     = D2 - Num[x] * 10000.0;
  }
return (UINT32)(D2*10000.0+D1);
}

static void
MulCRT(CRTNum *Num, ModInt Val,int Len)
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
** Things that don't need to be changed.  They just need to be here,
** is all.
*/
static void
ClearCRT(CRTNum *Num)
{int x;
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
{int x;
for (x=0;x<CRT_LEN;x++) printf("%08x",Num[x]);
printf("\n");
}

static void
CopyCRT(CRTNum *Dest, CRTNum *Src)
{int x;
for (x=0;x<CRT_LEN;x++) Dest[x]=Src[x];
}

#endif


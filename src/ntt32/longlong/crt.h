/*
** This file contains very low level routines required for the Chinese
** Remainder Theorem.  Due to performance reasons, these routines have
** to be written very carefully, and possibly use assembly language or
** compiler specific code.  The should also be 'inline' code because
** a subroutine call would probably cost as much as the operation itself.
**
** These versions depend on some form of 64 bit 'long long'.  Many 32 bit
** compilers supply software implemented 64 bit integers.
**
** I don't recommend the 'long long' code for 32 bit systems.
** First, I found a bug in DJGPP 2.7.2.1 that I had to work around
** (see the StripCRT() down there).  Second, somebody has given me
** a couple of 'horror' stories about GNU C's long long.  And at least
** one other processor (also GCC) fails when using 'long long'.  My code
** does seem to be correct.  It appears to be bugs in the 'long long'
** implementations.
**
** I have written the code so it never does more than one thing in
** a statement.  In other words, it never does things like a=b+c*d;
** Hopefully, that will help reduce the chance that bugs like that occur.
**
** Use at your own risk.  You could at least use them as a guide to
** write your own asm code.  If you've got a 64 bit processor, the
** 'long long' will map directly to normal 64 bit integer math, of course.
**
** If you need to write these in asm yourself, I am using the normal
** 'Big Endian' format (where the MSD is index 0).  You could change
** that, but its doubtful you'd need to, since it's just the array
** that I'm doing that way.  Your CPU shouldn't care.
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
** hardwired stuff.
*/
#define CRT_LEN NPrimes

typedef UINT32 CRTNum;


/******************************************************
********  GENERIC 64 BIT 'LONG LONG' VERSIONS  ********
******************************************************/

static UINT32
StripCRT(CRTNum *Num)
{UINT64 D;int x;
D=0;
for (x=0;x<CRT_LEN;x++)
  {
/*
** There seems to be a bug in DJGPP 2.7.2.1 because the following
** two statements can't be combined without Num[x] also being shifted.
*/
   D=D<<32;
   D=D + ((UINT64)Num[x]);
   Num[x] = D / 100000000UL;
   D      = D % 100000000UL;
  }
return (UINT32)D;
}

static void
MulCRT(CRTNum *Num, ModInt Val,int Len)
{UINT64 Prod;
 UINT64 Temp;

Prod=0;
Num+=CRT_LEN;
while (Len--)
  {
   Temp=((UINT64)(*(--Num)))*Val;
   Prod=Prod+Temp;
   *Num=Prod &  0xffffffff;
   Prod=Prod >> 32;
  }
}

static INLINE UINT32
Mulz2(UINT32 a, UINT32 b, UINT32 c, UINT32 *Dest)
{UINT64 Prod;

 Prod=((UINT64)a)*((UINT64)b);
 Prod+=c;
 Prod+=*Dest;
 *Dest=(UINT32)Prod;
 return (UINT32)(Prod>>32);
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

#else /* not inline crt */

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

static INLINE void
Add2CRT(CRTNum *Num1, CRTNum *Num2)
{UINT64 S;
S=((UINT64)Num1[CRT_LEN-1])+((UINT64)Num2[CRT_LEN-1]);
Num1[CRT_LEN-1]=S &  0xffffffff;S =S >> 32;
#if NPrimes==8
S=S+((UINT64)Num1[6]);S=S+((UINT64)Num2[6]);
Num1[6]=S &  0xffffffff;S =S >> 32;
S=S+((UINT64)Num1[5]);S=S+((UINT64)Num2[5]);
Num1[5]=S &  0xffffffff;S =S >> 32;
S=S+((UINT64)Num1[4]);S=S+((UINT64)Num2[4]);
Num1[4]=S &  0xffffffff;S =S >> 32;
S=S+((UINT64)Num1[3]);S=S+((UINT64)Num2[3]);
Num1[3]=S &  0xffffffff;S =S >> 32;
#endif
S=S+((UINT64)Num1[2]);S=S+((UINT64)Num2[2]);
Num1[2]=S &  0xffffffff;S =S >> 32;
S=S+((UINT64)Num1[1]);S=S+((UINT64)Num2[1]);
Num1[1]=S &  0xffffffff;S =S >> 32;
S=S+((UINT64)Num1[0]);S=S+((UINT64)Num2[0]);
Num1[0]=S &  0xffffffff;
}

#if 0
/* Alternative UINT64 version */
static void
Add2CRT(CRTNum *Num1, CRTNum *Num2)
{UINT64 S;
 int x;

S=0;
for (x=CRT_LEN-1;x>=0;x--)
  {
   S=S+((UINT64)Num1[x]);
   S=S+((UINT64)Num2[x]);
   Num1[x]=S &  0xffffffff;
   S      =S >> 32;
  }
}
#endif

/*
** Long Double (64 bit mantissa) versions of a few routines.
*/
#if 0
static UINT32
StripCRT(CRTNum *Num)
{int x;
 long double D;  /* NOTE that this requires 64 bit long doubles. */
 UINT32 z;

D=0.0;
for (x=0;x<CRT_LEN;x++)
  {
   D=((D*65536.0)*65536.0)+Num[x];
   z=D / 100000000.0;
   Num[x] = (UINT32) z;
   D = D - z * 100000000.0;
  }
return (UINT32)D;
}
#endif

#if 0
static void
MulCRT(CRTNum *Num, ModInt Val,int Len)
{long double Prod; /* Note 'long double' being used. */
 CRTNum Temp;

Prod=0;
Num+=CRT_LEN;
while (Len--)
  {
   Prod=((long double)*(--Num))*Val+Prod;
   *Num=Prod;
   Temp=(CRTNum)(Prod / 4294967296.0); /* to chop it to integer */
   Prod=Temp;
  }
}
#endif


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


#ifndef CRT_H_
#define CRT_H_ 1

/*
** Not to be changed!  This just provides a symbolic description of
** hardwired stuff.
*/
#define CRT_LEN NPrimes

typedef UINT32 CRTNum;

/* Used in the x86 specific mul routines below. */
static inline CRTNum
Mulz(CRTNum a, CRTNum b, CRTNum c, CRTNum *Dest)
{CRTNum r,z;
  asm volatile ("movl %2,%%eax
        mull %3
        addl %4,%%eax
        adcl $0,%%edx
        movl %%eax,%0
        movl %%edx,%1
      "
    : "=g" (z), "=g" (r)
    : "g" (a) , "g" (b) , "g" (c)
    : "%eax" , "%edx" , "cc", "memory"
    );
*Dest=z;
return r;
}

/* Doesn't work with gcc 2.81
static inline UINT32
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
    : "%eax" , "%edx" , "cc" , "memory"
    );
*Dest=z;
return r;
}
*/

static inline UINT32
Mulz2(UINT32 a, UINT32 b, UINT32 c, UINT32 *Dest)
{UINT32 r,z;
 z=*Dest;
  asm volatile ("movl %2,%%eax
        mull %3
        addl %4,%%eax
        adcl $0,%%edx
        addl %0,%%eax
        adcl $0,%%edx
        movl %%eax,%0
        movl %%edx,%1
      "
    : "=m" (z), "=m" (r)
    : "m" (a) , "m" (b) , "m" (c)
    : "%eax" , "%edx" , "%cc" , "%memory"
    );
*Dest=z;
return r;
}


static inline void
MulCRT(CRTNum *Num, ModInt Val,size_t Len)
{CRTNum Carry;

Carry=0;
Num+=CRT_LEN;
while (Len--)
  {
   --Num;
   Carry=Mulz(*Num,Val,Carry,Num);
  }
}

#ifdef UNROLLED_CRT
static inline void
MulCRT1(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[CRT_LEN-1],Val,Carry,&Res[CRT_LEN-1]);
Res[CRT_LEN-2]=Carry;
}

static inline void
MulCRT2(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[CRT_LEN-1],Val,Carry,&Res[CRT_LEN-1]);
Carry=Mulz2(Num[CRT_LEN-2],Val,Carry,&Res[CRT_LEN-2]);
Res[CRT_LEN-3]=Carry;
}

static inline void
MulCRT3(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[CRT_LEN-1],Val,Carry,&Res[CRT_LEN-1]);
Carry=Mulz2(Num[CRT_LEN-2],Val,Carry,&Res[CRT_LEN-2]);
Carry=Mulz2(Num[CRT_LEN-3],Val,Carry,&Res[CRT_LEN-3]);
Res[CRT_LEN-4]=Carry;
}

#if NPrimes==8

static inline void
MulCRT4(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[7],Val,Carry,&Res[7]);
Carry=Mulz2(Num[6],Val,Carry,&Res[6]);
Carry=Mulz2(Num[5],Val,Carry,&Res[5]);
Carry=Mulz2(Num[4],Val,Carry,&Res[4]);
Res[3]=Carry;
}

static inline void
MulCRT5(UINT32 *Res,UINT32 *Num, UINT32 Val)
{UINT32 Carry=0;
Carry=Mulz2(Num[7],Val,Carry,&Res[7]);
Carry=Mulz2(Num[6],Val,Carry,&Res[6]);
Carry=Mulz2(Num[5],Val,Carry,&Res[5]);
Carry=Mulz2(Num[4],Val,Carry,&Res[4]);
Carry=Mulz2(Num[3],Val,Carry,&Res[3]);
Res[2]=Carry;
}

static inline void
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

static inline void
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
#endif /* 8 primes */

#else /* loop CRT */

static inline void
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

#if NPrimes==8

static inline void
Add2CRT(CRTNum *Num1, CRTNum *Num2)
{
  asm volatile ("movl %0,%%edi
       movl %1,%%esi
       movl 28(%%esi) ,%%eax ; addl %%eax,28(%%edi)
       movl 24(%%esi) ,%%eax ; adcl %%eax,24(%%edi)
       movl 20(%%esi) ,%%eax ; adcl %%eax,20(%%edi)
       movl 16(%%esi) ,%%eax ; adcl %%eax,16(%%edi)
       movl 12(%%esi) ,%%eax ; adcl %%eax,12(%%edi)
       movl  8(%%esi) ,%%eax ; adcl %%eax, 8(%%edi)
       movl  4(%%esi) ,%%eax ; adcl %%eax, 4(%%edi)
       movl   (%%esi) ,%%eax ; adcl %%eax,  (%%edi)
      "
    :
    : "m" (Num1) , "m" (Num2)
    : "%eax" , "%esi" , "%edi" , "cc" , "memory"
    );
}

static inline UINT32
StripCRT(CRTNum *Num)
{UINT32 Remain;

  asm volatile ("movl %1,%%ebx
       movl $0,%%edx
       movl   (%%ebx) ,%%eax ; divl %2 ; movl %%eax,  (%%ebx)
       movl  4(%%ebx) ,%%eax ; divl %2 ; movl %%eax, 4(%%ebx)
       movl  8(%%ebx) ,%%eax ; divl %2 ; movl %%eax, 8(%%ebx)
       movl 12(%%ebx) ,%%eax ; divl %2 ; movl %%eax,12(%%ebx)
       movl 16(%%ebx) ,%%eax ; divl %2 ; movl %%eax,16(%%ebx)
       movl 20(%%ebx) ,%%eax ; divl %2 ; movl %%eax,20(%%ebx)
       movl 24(%%ebx) ,%%eax ; divl %2 ; movl %%eax,24(%%ebx)
       movl 28(%%ebx) ,%%eax ; divl %2 ; movl %%eax,28(%%ebx)
       movl %%edx,%0
      "
    : "=m" (Remain)
    : "m" (Num), "rm" (BaseVar)
    : "%eax" , "%edx" , "%ebx" , "cc" , "memory"
    );

return Remain;
}

#else /* 4 primes */

static inline void
Add2CRT(CRTNum *Num1, CRTNum *Num2)
{
  asm volatile ("movl %0,%%edi
       movl %1,%%esi
       movl 12(%%esi) ,%%eax ; addl %%eax,12(%%edi)
       movl  8(%%esi) ,%%eax ; adcl %%eax, 8(%%edi)
       movl  4(%%esi) ,%%eax ; adcl %%eax, 4(%%edi)
       movl   (%%esi) ,%%eax ; adcl %%eax,  (%%edi)
      "
    :
    : "m" (Num1) , "m" (Num2)
    : "%eax" , "%esi" , "%edi" , "cc" , "memory"
    );
}

static inline UINT32
StripCRT(CRTNum *Num)
{UINT32 Remain;

  asm volatile ("movl %1,%%ebx
       movl $0,%%edx
       movl   (%%ebx) ,%%eax ; divl %2 ; movl %%eax,  (%%ebx)
       movl  4(%%ebx) ,%%eax ; divl %2 ; movl %%eax, 4(%%ebx)
       movl  8(%%ebx) ,%%eax ; divl %2 ; movl %%eax, 8(%%ebx)
       movl 12(%%ebx) ,%%eax ; divl %2 ; movl %%eax,12(%%ebx)
       movl %%edx,%0
      "
    : "=m" (Remain)
    : "m" (Num), "rm" (BaseVar)
    : "%eax" , "%edx" , "%ebx" , "cc" , "memory"
    );

return Remain;
}

#endif /* 8/4 primes */

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


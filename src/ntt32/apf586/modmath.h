#ifndef MODMATH_H_
#define MODMATH_H_ 1
/*
** This file contains the very low level modular math routines.
** Due to performance reasons, these routines have to be written very
** carefully, and possibly use assembly language or compiler specific
** code.  The should also be 'inline' code because a subroutine call
** would probably cost as much as the operation itself.
*/

/******************************************************
********  x86 SPECIFIC ROUTINES.  386/486 ASM. ********
******************************************************/

/*
** The Pentium could probably do better with the FPU registers.
** Jason P. does that in his NTT586, and Mikko Tommila does that
** in his APFloat.
*/

static inline
ModInt ModMul(ModInt a, ModInt b)
/* Do a 32x32=64 & then get 32 bit modulo */
{ModInt rem;
  asm ("mull %2; divl %3" : "=&d" (rem)
                          : "%a" (a), "rm" (b), "rm" (Prime)
                          : "%eax", "edx", "cc");
return rem;
}

#if 1
static inline
ModInt ModAdd(ModInt a, ModInt b)
{
 ModInt r;

 asm ("addl %2, %1;
        cmpl %3, %1;
        jb 0f;
        subl %3, %1; 0:"
                         : "=&r" (r)
                         : "%0" (a), "g" (b), "rm" (Prime)
                         : "cc");
 return r;
}

static inline
ModInt ModSub(ModInt a, ModInt b)
{
 ModInt r;

 asm ("subl %2, %1;
        jnc 0f;
        addl %3, %1; 0:"
                         : "=&r" (r)
                         : "0" (a), "g" (b), "rm" (Prime)
                         : "cc");
 return r;
}
#else
static inline
ModInt ModAdd(ModInt a, ModInt b)
{ModInt r;

 asm ("addl %2, %1;
       subl %3, %1;
       movl %1, %4;
       sarl $31,%4;
       andl %3, %4;
       addl %1, %4;
       movl %4, %0;
      "
      : "=&r" (r)
      : "%0" (a), "g" (b), "rm" (Prime), "q" (r)
      : "cc");
 return r;
}

static inline
ModInt ModSub(ModInt a, ModInt b)
{
 ModInt r;

 asm ("subl %2, %1;
       movl %1, %4;
       sarl $31,%4;
       andl %3, %4;
       addl %1, %4;
       movl %4, %0;
      "
      : "=&r" (r)
      : "%0" (a), "g" (b), "rm" (Prime), "q" (r)
      : "cc");
 return r;
}

#endif



#endif


#if 0

/*
** These two versions depend on the primes being 31 bits or less,
** and the integers being 32 bits, with the high bit being the sign.
**
** If you have a processor that has a large pipeline that doesn't
** handle 'random' jumps very well, you can use these two functions.
** These versions are straight code and wont generate any branches.
**
** On a 486, the above jump versions are actually faster than
** these.  The Pentium appears to behave the same.  (Of course,
** Jason's ntt586 uses assembly language versions of these, and
** does so very profitably.)
**
** Your processor may be different.  I'd suggest you time
** a few combinations to find out how your processor does.
*/
static
ModInt ModAdd(ModInt a,ModInt b)
{INT32 Res,x; /* signed 32 bit int.... */
Res=a+b-Prime;
x=Res >> 31; /* if it underflowed into the sign bit, set all bits */
/* If 64 bit itegers, that might need to be x=Res >> 63 */
Res=Res+(x & Prime);
return Res;
}

static
ModInt ModSub(ModInt a,ModInt b)
{INT32 Res,x; /* signed 32 bit int.... */
Res=a-b;
x=Res >> 31; /* if it underflowed into the sign bit, set all bits */
/* If 64 bit itegers, that might need to be x=Res >> 63 */
Res = Res + (x & Prime);
return Res;
}

#endif



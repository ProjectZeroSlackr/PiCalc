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

static inline
ModInt ModMul(ModInt a, ModInt b)
/* Do a 32x32=64 & then get 32 bit modulo */
{ModInt rem;
  asm volatile ("mull %2; divl %3"
                          : "=d" (rem)
                          : "a" (a), "rm" (b), "rm" (Prime)
                          : "cc");
return rem;
}

#if PRIME_BITS==32
static inline
ModInt ModAdd(ModInt a, ModInt b)
{
 ModInt r;
 asm volatile ("subl %3, %1;
       addl %2, %1
       jc 0f
       addl %3, %1
     0:"
      : "=&r" (r)
      : "0" (a), "g" (b), "rm" (Prime)
      : "cc");

 return r;
}
#else
static inline
ModInt ModAdd(ModInt a, ModInt b)
{
 ModInt r;

 asm volatile ("addl %2, %1;
        cmpl %3, %1;
        jb 0f;
        subl %3, %1; 0:"
                         : "=&r" (r)
                         : "0" (a), "g" (b), "rm" (Prime)
                         : "cc");
 return r;
}
#endif

static inline
ModInt ModSub(ModInt a, ModInt b)
{
 ModInt r;

 asm volatile ("subl %2, %1;
        jnc 0f;
        addl %3, %1; 0:"
                         : "=&r" (r)
                         : "0" (a), "g" (b), "rm" (Prime)
                         : "cc");
 return r;
}

#endif


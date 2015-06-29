#ifndef MODMATH_H_
#define MODMATH_H_ 1

/*
** This file contains the very low level modular math routines.
** Due to performance reasons, these routines have to be written very
** carefully, and possibly use assembly language or compiler specific
** code.  The should also be 'inline' code because a subroutine call
** would probably cost as much as the operation itself.
**
** All of these routines are provided in generic versions so you can
** at least run the code (SLOWLY!) while you are writing the higher
** performance compiler/processor specific code. (such as asm.)
**
** Also, I'd like to ask that if you write any custom  routines  for
** your processor or compiler, please send a copy to  me  so  I  can
** include them in any future  versions.  Your doing this might save
** somebody less skilled than you a lot of frustration.
*/

/************************************************
********  GENERIC 'ANSI/ISO C' VERSIONS  ********
************************************************/

/* No 64 bits at all. Very slow. */
/*
** The 'floor()' may be very slow on some processors and/or compilers.
** Straight assignment to a ModInt might be faster.  Or some sort of
** assembly...
**
** I don't really like these routines.  They look kludgy.  A standard
** integer mul will generate the lower 32 bits of a 32x32 mul.  The
** floating point unit will generate the upper bits of a 32x32 mul.
** We can then use the two to get the full 64 bits.  However, since
** we are also needing to do the modulo here, we throw that in too.
** When the prime is known, we can multiply by the reciprocal, since
** that is faster than a division.
**
** These could be optimized, but the goal here isn't speed, but to
** provide you with routines that actually work, regardless of what
** compiler or processor you are using.
*/
static INLINE ModInt
ModMul(ModInt a, ModInt b)
{ModInt rem;
rem = a * b;
rem = rem - Prime * ((ModInt) floor(0.5+RecipPrime * ((double) a) * ((double) b)));
if (rem < 0) rem +=Prime;
return rem;
}

#if 1
/*
** Nice generic (and fairly fast) modular add & sub
** The compiler should generate fairly efficient code for this,
** you shouldn't even need asm.
*/
static INLINE ModInt
ModAdd(ModInt a,ModInt b)
{ModInt Res;
Res=a+b-Prime;
if (Res < 0) Res+=Prime;
return Res;
}

static INLINE ModInt
ModSub(ModInt a,ModInt b)
{ModInt Res;
Res=a-b;
if (Res < 0) Res+=Prime;
return Res;
}

#else

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
static INLINE ModInt
ModAdd(ModInt a,ModInt b)
{INT32 Res,x; /* signed 32 bit int.... */
Res=a+b-Prime;
x=Res >> 31; /* if it underflowed into the sign bit, set all bits */
/* If 64 bit itegers, that might need to be x=Res >> 63 */
Res=Res+(x & Prime);
return Res;
}

static INLINE ModInt
ModSub(ModInt a,ModInt b)
{INT32 Res,x; /* signed 32 bit int.... */
Res=a-b;
x=Res >> 31; /* if it underflowed into the sign bit, set all bits */
/* If 64 bit itegers, that might need to be x=Res >> 63 */
Res = Res + (x & Prime);
return Res;
}

#endif

#endif

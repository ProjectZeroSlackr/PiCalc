#ifndef MODMATH_H_
#define MODMATH_H_ 1
/*
** This file contains the very low level modular math routines.
** Due to performance reasons, these routines have to be written very
** carefully, and possibly use assembly language or compiler specific
** code.  They should also be 'inline' code because a subroutine call
** would probably cost as much as the operation itself.
*/

static INLINE
ModInt ModMul(ModInt a, ModInt b)
/* Do a 32x32=64 & then get 32 bit modulo */
/*
** We have to divide by Prime, rather than multiplying by the
** reciprocal.  That's because we would need additional 'guard' bits
** to allow for the inaccuracy of the reciprocal.  (ie: 1/3=0.3333 != 1/3)
**
** Although the 586 version is only using 31 bit primes, meaning our
** modular numbers are 31 bits and the product would only be 62 bits of
** the 64 bits in a 'long double' I don't think two extra bits are
** enough to ensure 100% accuracy.  I could have done it more along
** the lines of now the ntt32/generic does it, but I didn't want
** to mess with that.
*/
{long double P;
 ModInt R;
 P=((long double)a)*((long double)b);
 R=(ModInt)(P/Prime);
 R=P-((long double)R)*Prime;

 return R;
}

static INLINE
ModInt ModAdd(ModInt a, ModInt b)
{UINT32 r;
 r=a+b;
 if (r >= Prime) r-=Prime;
 return r;
}

static INLINE
ModInt ModSub(ModInt a, ModInt b)
{INT32 r;
 r=a-b;
 if (r < 0) r+=Prime;
 return r;
}

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
** These versions can be faster in some situations.  In other cases,
** where you just don't have the spare regs or something, they
** can be slower.
**
** Your milage may vary.
*/
static INLINE
ModInt ModAdd(ModInt a,ModInt b)
{INT32 Res,x; /* signed 32 bit int.... */
Res=a+b-Prime;
x=Res >> 31; /* if it underflowed into the sign bit, set all bits */
/* If 64 bit itegers, that might need to be x=Res >> 63 */
Res=Res+(x & Prime);
return Res;
}

static INLINE
ModInt ModSub(ModInt a,ModInt b)
{INT32 Res,x; /* signed 32 bit int.... */
Res=a-b;
x=Res >> 31; /* if it underflowed into the sign bit, set all bits */
/* If 64 bit itegers, that might need to be x=Res >> 63 */
Res = Res + (x & Prime);
return Res;
}

#endif



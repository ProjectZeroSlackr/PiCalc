#ifndef MODMATH_H_
#define MODMATH_H_ 1

/*
** This file contains the very low level modular math routines.
** Due to performance reasons, these routines have to be written very
** carefully, and possibly use assembly language or compiler specific
** code.  The should also be 'inline' code because a subroutine call
** would probably cost as much as the operation itself.
**
** All of these routines are provided in 'generic' versions so you can
** at least run the code (SLOWLY!) while you are writing the higher
** performance compiler/processor specific code. (such as asm.)
**
** I say 'generic' because they do require some sort of 64 bit
** unsigned integers, of course.  On a 64 bit computer, it's no
** problem, of course.  On a 32 bit computer (such as what I used
** for development & testing), it has to be faked by the compiler.
**
** Also, I'd like to ask that if you write any custom  routines  for
** your processor or compiler, please send a copy to  me  so  I  can
** include them in any future  versions.  Your doing this might save
** somebody less skilled than you a lot of frustration.
*/

#if 0
/*
** _IF_ you have some sort of 128 bit integers, and can multiply
** two 64 bit integers, get a 128 bit result, and take the modulus
** of a 64 bit integer from that, then you can use something such
** as below.  Or you can use a larger version of the 32 bit ModMul
** in the NTT32/GENERIC/MODMATH.H file.  Or....  You get the point.
*/

static INLINE ModInt
ModMul(ModInt a, ModInt b)
{UINT64 Temp;
 ModInt rem;
Temp = ((UINT128)a) * ((UINT128)b);
rem = Temp % Prime;
return rem;
}
#endif

/************************************************
********  GENERIC 'ANSI/ISO C' VERSIONS  ********
************************************************/

static INLINE ModInt
ModMul(ModInt a,ModInt b)
/*
** Since the prime is hardwired to do _ONLY_ 2^64-2^32+1,
** this is not a general purpose modular multiply.
**
** From Mikko Tommila's APFloat package.
*/
{
 UINT64 al, ah, bl, bh, rl, rh, c;
 UINT64 t;

 al = a & 0xFFFFFFFF;
 ah = a >> 32;

 bl = b & 0xFFFFFFFF;
 bh = b >> 32;

 rl = al * bl;

 c = al * bh;
 rh = c >> 32;
 c <<= 32;
 rl += c;
 if (rl < c) rh++;

 c = ah * bl;
 rh += c >> 32;
 c <<= 32;
 rl += c;
 if (rl < c) rh++;

 rh += ah * bh;

 /* modulo reduction */

 /* 1st shift */
 t = rh;
 c = rh << 32;
 rh >>= 32;
 t = rl - t;
 if (t > rl) rh--;

 rl = t + c;
 if (rl < t) rh++;

 /* 2nd shift */
 t = rh;
 c = rh << 32;
 rh >>= 32;
 t = rl - t;
 if (t > rl) rh--;

 rl = t + c;
 if (rl < t) rh++;

 /* Final check */
 return (rh || rl >= Prime ? rl - Prime : rl);
}

static INLINE ModInt
ModAdd(ModInt a,ModInt b)
{ModInt Sum;
Sum=a+b;
if (Sum < a) Sum-=Prime;
if (Sum >= Prime) Sum-=Prime;
return Sum;
}

static INLINE ModInt
ModSub(ModInt a,ModInt b)
{ModInt Dif;
Dif=a-b;
if (Dif > a) Dif+=Prime;
return Dif;
}

#endif


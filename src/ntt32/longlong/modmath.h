#ifndef MODMATH_H_
#define MODMATH_H_ 1

/*
** This file contains the very low level modular math routines.
** Due to performance reasons, these routines have to be written very
** carefully, and possibly use assembly language or compiler specific
** code.  The should also be 'inline' code because a subroutine call
** would probably cost as much as the operation itself.
**
** These depend on some form of 64 bit 'long long' integer.  Many 32
** bit compilers supply these.
**
** I don't recommend the 'long long' code.  First, I found a bug in
** DJGPP 2.7.2.1 that I had to work around (see the StripCRT() down
** there).  Second, somebody has given me a couple of 'horror' stories
** about GNU C's long long.  And at least one other processor (also
** GCC) fails when using 'long long'.  My code does seem to be
** correct.  It appears to be bugs in the 'long long' implementations.
** Use at your own risk.  You could at least use them as a guide to
** write your own asm code.  If you've got a 64 bit processor, the
** 'long long' will map directly to normal 64 bit integer math, of course.
**
** Also, I'd like to ask that if you write any custom  routines  for
** your processor or compiler, please send a copy to  me  so  I  can
** include them in any future  versions.  Your doing this might save
** somebody less skilled than you a lot of frustration.
*/

/******************************************************
********  GENERIC 64 BIT 'LONG LONG' VERSIONS  ********
******************************************************/

/* Not very fast with DJGPP because the % requires a subroutine call. */
static INLINE ModInt
ModMul(ModInt a, ModInt b)
{UINT64 Temp;
ModInt rem;
 Temp= ((UINT64)a) * ((UINT64)b);
 rem = Temp % Prime;
return rem;
}

#if 0
static INLINE ModInt
ModAdd(ModInt a,ModInt b)
{INT64 Res;
Res=a;
Res=Res+b;
Res=Res-Prime;
if (Res < 0) Res+=Prime;
return Res;
}

static INLINE ModInt
ModSub(ModInt a,ModInt b)
{INT64 Res;
Res=a;
Res=Res-b;
if (Res < 0) Res+=Prime;
return Res;
}
#else
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
#endif

#endif


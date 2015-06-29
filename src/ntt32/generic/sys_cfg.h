/*
** This file contains various system specific items.
**
** This include various data types, macros, etc. that you might need.
*/

#ifndef SYS_CFG_H_
#define SYS_CFG_H_ 1

/*
** The number of primes to use.  This should be either 4 or 8.
** Whether the primes should be 31 or 32 bits long.
** You should use 32 bit primes because I don't think the
** generic version is reliable for 32 bit primes.
**
** See primes.h for computation limits for the combinations.
*/
#define NPrimes     8
#define PRIME_BITS 31

/*
** Whether to do the CRT as unrolled inline code, or a loop.
** Unrolled code is better.  But regular loop code is quicker to port.
*/
#define UNROLLED_CRT 1

/*
** Whether to do some things as vector style operations.
** This is a bad idea, unless the version has custom vector
** routines that operate faster than the regular modular math
** functions.
**
** The generic version is so generic that it's not capable of
** doing anything efficiently, including vector operations.
*/

/*
** Whether to do a vector style CRT
#define VECTOR_CRT  1
*/

/*
** Whether to do a vector style NTT
** (This is one half of doing a full vector NTT.)
#define VECTOR_NTT  1
*/

/*
** Whether to do a vector style RNTT.
** (This is the other half of doing a full vector NTT.)
** This should be a power of two.  Preferably either 2 or 4
** Since doing a RNTT takes a little more time, this makes
** sure that the vector will be long enough to offset the
** startup cost.
#define VECTOR_RNTT 4
*/


#if defined(VECTOR_CRT) && !defined(UNROLLED_CRT)
#error The vector CRT needs the unrolled CRT routines.
#endif

/*
** If your compiler supports the 'inline' keyword, define it below
#define INLINE inline
#define INLINE
*/
#define INLINE


typedef unsigned           int UINT32;
typedef   signed           int  INT32;
#if PRIME_BITS==32
typedef UINT32                  ModInt;
#else
typedef INT32                   ModInt;
#endif
typedef ModInt                  FFT_DATA_TYPE;


#define NEED_RECIP_PRIME 1
#define SetModPrime(Zq) Prime=Consts[Zq].Prime;RecipPrime=Consts[Zq].RecipPrime;

#endif



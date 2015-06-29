#ifndef PRIMES_H_
#define PRIMES_H_ 1

#include "sys_cfg.h"

/*
** The maximum number of decimal digits the NTT can handle.
**
** This is based on the multiplication pyramid size vs. the size of
** the primes product, and the Max FFTLen as shown by FindPrime.c
** (times the number of primes we are using, times 8 (the number of
** digits per int32) and divided by 2 (for the zero padding.))
**
** Of course, you chose the lowest limit.
*/

#if   (NPrimes==8) && (PRIME_BITS==32)

/* Limited by the powers of the primes. 2^26 */
#define MaxFFTLen 1073741824U

#elif (NPrimes==8) && (PRIME_BITS==31)

/* Limited by the powers of the primes.  2^24 */
#define MaxFFTLen 268435456U

#elif (NPrimes==4) && (PRIME_BITS==32)

/* Limited by the multiplication pyramid size. */
#define MaxFFTLen  33554432U
/* Just _barely_! */

#elif (NPrimes==4) && (PRIME_BITS==31)

/* Limited by the multiplication pyramid size. */
#define MaxFFTLen 2097152U

#else

#error Unknown setting of NPrimes and PRIME_BITS

#endif

/*
** How many INT32's do we put into a ModInt for the NTT?
**
** Each INT32 holds 8 decimal digit.
** For the 4 prime NTT, we put 16 decimal digits into each ModInt.
** For the 8 prime NTT, we put 32 decimal digits into each ModInt.
*/
#define RawModInts (NPrimes/2)

/*
** Macro to convert a number length to the length of the NTT
** This INCLUDES the doubling for the zero padding.
*/
#define CalcNTTLen(zw) (((zw)*2)/RawModInts)

/*
** This converts a NTTLen (from above) to a NumLen.
*/
#define CalcNumLen(zw) (((zw)*RawModInts)/2)

extern ModInt PrimeList[NPrimes];
extern ModInt PrimvRootList[NPrimes];

#endif



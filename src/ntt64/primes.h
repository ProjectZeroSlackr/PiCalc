#ifndef PRIMES_H_
#define PRIMES_H_ 1

#include "sys_cfg.h"

#define MaxFFTLen 536870912UL
/*
** I don't really know with certainty how big we can get with a
** 64 bit prime, but 512m is enough.
*/

/*
** We're only using 1 prime, but it's needed elsewhere.  Really,
** it's just for 'logical' consistency with the multi-prime ntt32.
*/
#define NPrimes 1

/*
** How many INT32's do we put into a ModInt for the NTT?
**
** Each INT32 holds 8 decimal digit.
** For the 4 prime NTT, we put 16 decimal digits into each ModInt.
** For the 8 prime NTT, we put 32 decimal digits into each ModInt.
*/
//////#define RawModInts (NPrimes/2)

/*
** Macro to convert a number length to the length of the NTT
** This INCLUDES the doubling for the zero padding.
*/
#define CalcNTTLen(zw) ((zw)*4)

/*
** This converts a NTTLen (from above) to a NumLen.
*/
#define CalcNumLen(zw) ((zw)/4)


extern ModInt PrimeList[NPrimes];
extern ModInt PrimvRootList[NPrimes];

#endif



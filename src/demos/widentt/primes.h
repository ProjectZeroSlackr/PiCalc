#ifndef PRIMES_H_
#define PRIMES_H_ 1


#define MaxFFTLen 536870912UL
/*
** I don't really know with certainty how big we can get with a
** 64 bit prime, but 512m is enough.
*/

#define Mod_SIZE 16

/*
** How many INT32's do we put into a ModInt for the NTT?
**
** Each INT32 holds 8 decimal digit.  The prime we are using
** is for 32 digits.
*/
#define RawModInts 4

/*
** Macro to convert a number length to the length of the NTT
** This INCLUDES the doubling for the zero padding.
*/
#define CalcNTTLen(zw) (((zw)/RawModInts)*2)

/*
** This converts a NTTLen (from above) to a NumLen.
*/
#define CalcNumLen(zw) (((zw)/2)*RawModInts)


#endif



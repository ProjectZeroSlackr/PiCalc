#ifndef PRIMES_H_
#define PRIMES_H_ 1


#define MaxFFTLen 536870912UL
/*
** Actually higher, but I was going to
*/

/*
** How many INT32's do we put into a ModInt for the FGT?
**
** Each INT32 holds 8 decimal digit.
*/
#define RawModInts 1

/*
** Macro to convert a number length to the length of the FGT
** This INCLUDES the doubling for the zero padding.
** The FGT is 'complex', so we are putting the data into
** both the 'real' and 'imaginary' parts.
*/
#define CalcFGTLen(zw) (((zw)/(2*RawModInts))*2)

/*
** This converts a FGTLen (from above) to a NumLen.
*/
#define CalcNumLen(zw) (((zw)/2)*(2*RawModInts))


#endif



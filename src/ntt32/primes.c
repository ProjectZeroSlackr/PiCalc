#include "primes.h"

#if (NPrimes==8) && (PRIME_BITS==32)

ModInt PrimeList[NPrimes]=
  {
   1811939329UL, /* Bits= 30.754 */ /* 2^26*27+1 */
   2013265921UL, /* Bits= 30.906 */ /* 2^27*17+1 */
   2281701377UL, /* Bits= 31.087 */ /* 2^27*17+1 */
   2483027969UL, /* Bits= 31.209 */ /* 2^26*37+1 */
   2885681153UL, /* Bits= 31.426 */ /* 2^26*43+1 */
   3221225473UL, /* Bits= 31.584 */ /* 2^30*3+1  */
   3489660929UL, /* Bits= 31.700 */ /* 2^28*13+1 */
   3892314113UL, /* Bits= 31.857 */ /* 2^27*29+1 */
  };       /* Total Bits=250.523 */

/* The primative roots of unity for the primes. */
ModInt PrimvRootList[NPrimes]={13,31,3,3,3,5,3,3};

/* Limited by the powers of the primes. 2^26
#define MaxFFTLen 1073741824U
*/

#elif (NPrimes==8) && (PRIME_BITS==31)

ModInt PrimeList[NPrimes]=
  {
    754974721, /* Bits= 29.47 */ /* 2^24*45+1  */
   1107296257, /* Bits= 30.02 */ /* 2^25*33+1  */
   1224736769, /* Bits= 30.17 */ /* 2^24*73+1  */
   1711276033, /* Bits= 30.65 */ /* 2^25*51+1  */
   1811939329, /* Bits= 30.73 */ /* 2^26*27+1  */
   2013265921, /* Bits= 30.89 */ /* 2^27*17+1  */
   2113929217, /* Bits= 30.96 */ /* 2^25*63+1  */
   2130706433  /* Bits= 30.97 */ /* 2^24*127+1 */
  };     /* Total Bits=243.87 */

/* The primative roots of unity for the primes. */
ModInt PrimvRootList[NPrimes]={11,31,3,29,13,31,5,3};

/* Limited by the powers of the primes.  2^24
#define MaxFFTLen 268435456U
*/

#elif (NPrimes==4) && (PRIME_BITS==32)

ModInt PrimeList[NPrimes]=
  {
   3942645761UL, /* Bits= 31.876 */ /* 2^24*235+1 */
   4076863489UL, /* Bits= 31.924 */ /* 2^24*243+1 */
   4194304001UL, /* Bits= 31.965 */ /* 2^25*125+1 */
   4253024257UL  /* Bits= 31.985 */ /* 2^24*507+1 */
  };       /* Total Bits=127.7529557 */

/* The primative roots of unity for the primes. */
ModInt PrimvRootList[NPrimes]={3,7,3,5};

/* Limited by the multiplication pyramid size.
#define MaxFFTLen  33554432U
 Just _barely_! */

#elif (NPrimes==4) && (PRIME_BITS==31)

ModInt PrimeList[NPrimes]=
  {
   2013265921, /* Bits= 30.89 */ /* 2^27*3*5+1 */
   2088763393, /* Bits= 30.96 */ /* 2^23*249+1 */
   2113929217, /* Bits= 30.96 */ /* 2^25*63+1  */
   2130706433  /* Bits= 30.97 */ /* 2^24*127+1 */
  };     /* Total Bits=123.55 */

/* The primative roots of unity for the primes. */
ModInt PrimvRootList[NPrimes]={31,5,5,3};

/* Limited by the multiplication pyramid size.
#define MaxFFTLen 2097152U
*/

#else

#error Unknown setting of NPrimes and PRIME_BITS

#endif



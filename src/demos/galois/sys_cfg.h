/*
** This file contains various system specific items.
**
** This include various data types, macros, etc. that you might need.
*/

#ifndef SYS_CFG_H_
#define SYS_CFG_H_ 1

#include "primes.h"

/*
** If your compiler supports the 'inline' keyword, define it below
#define INLINE inline
#define INLINE
*/
#define INLINE inline

typedef unsigned long long int UINT64;
typedef unsigned           int UINT32;
typedef unsigned short     int UINT16;
typedef   signed           int  INT32;

typedef struct {UINT32 Hi,Med,Low;} ModInt;
typedef struct {ModInt r,i;}   CmplxModInt;
typedef CmplxModInt            FFT_DATA_TYPE;



#endif



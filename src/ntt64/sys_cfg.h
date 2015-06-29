/*
** This file contains various system specific items.
**
** This include various data types, macros, etc. that you might need.
*/

#ifndef SYS_CFG_H_
#define SYS_CFG_H_ 1

/*
** If your compiler supports the 'inline' keyword, define it below
#define INLINE inline
#define INLINE
*/
#define INLINE inline

typedef unsigned           int UINT32;
typedef unsigned long long int UINT64;
typedef   signed           int  INT32;
typedef UINT64                  ModInt;
typedef ModInt                  FFT_DATA_TYPE;

#endif



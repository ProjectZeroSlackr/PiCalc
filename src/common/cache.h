#ifndef CACHE_H_
#define CACHE_H_ 1

#include "pi.h"
#include "bigint.h"

struct FFTHdr
  {
   size_t NumLen;
   size_t FFTLen;
   BigInt ChkNum;
   char  *Name;
   int Tag;
   int Pass; /* Which 'pass' of the FFT is cached here. */
#ifdef VIRTUAL_CACHE
   void *Mem;
#endif
  };

extern struct FFTHdr *FFTCash;
extern int    MaxFFTCacheLines;

int  CheckFFTCache(BigInt Num, size_t NumLen, int Tag, int Pass);
void DeInitFFTCache(void);
void DeleteFFTCache(int Tag);
void FlushFFTCache(int Tag);
void InitFFTCache(size_t Len);
void LoadFFTFromCache(int Cache, void *FFTNum);
int  SaveFFTIntoCache(void *FFTNum, size_t FFTLen, BigInt Num, size_t NumLen,
                      int Tag, int Pass);
void SaveFFTIntoCache0(void *FFTNum, size_t FFTLen);

#endif


#ifndef PORT_H_
#define PORT_H_

int TestKeyboard(void);
void CheckDiskCache(void);
size_t GetPhysicalMemSize(void);
size_t GetVirtualMemSize(void);

#ifdef DO_TIMINGS

#ifdef DO_PENTIUM_TIMINGS

#ifndef USE_MSR
#define StartTimer(t) {                                                        \
              asm volatile (".byte 15,49 ; subl %%eax,%0;sbbl %%edx,4+%0 "     \
                            :                                                  \
                            : "o" (t)                                          \
                            : "%eax", "%edx", "cc", "memory");                 \
                      }

#define StopTimer(t) {                                                         \
              asm volatile (".byte 15,49 ; addl %%eax,%0;adcl %%edx,4+%0 "     \
                            :                                                  \
                            : "o" (t)                                          \
                            : "%eax", "%edx", "cc", "memory");                 \
                      }
#else
/* Use the Pentium/MMX MSR Performance monitoring counters, instead. */
#define StartTimer(t) {                                                        \
asm volatile ("movl $18,%%ecx ; .byte 15,50 ; subl %%eax,%0 ; sbbl %%edx,4+%0 "\
              :                                                                \
              : "o" (t)                                                        \
              : "%eax", "%edx", "%ecx", "cc", "memory");                       \
                      }

#define StopTimer(t) {                                                         \
asm volatile ("movl $18,%%ecx ; .byte 15,50 ; addl %%eax,%0 ; adcl %%edx,4+%0 "\
              :                                                                \
              : "o" (t)                                                        \
              : "%eax", "%edx", "%ecx", "cc", "memory");                       \
                     }

#endif

#elif defined(__DJGPP__)
#include <go32.h>
#include <libc/farptrgs.h>
#define StartTimer(t) {(t)-=_farpeekl(_dos_ds, 0x46c);}
#define StopTimer(t)  {(t)+=_farpeekl(_dos_ds, 0x46c);}

#else
#define StartTimer(t) {(t)-=time(NULL);}
#define StopTimer(t)  {(t)+=time(NULL);}
#endif

#else
#define StartTimer(t)
#define StopTimer(t)
#endif



#endif


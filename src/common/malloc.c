#include <stdlib.h>
#include <string.h>

#ifdef USE_SPECIAL_DJGPP_MALLOC

/****************************************************************
Copyright 1990, 1994 by AT&T, Lucent Technologies and Bellcore.

Permission to use, copy, modify, and distribute this software
and its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the names of AT&T, Bell Laboratories,
Lucent or Bellcore or any of their entities not be used in
advertising or publicity pertaining to distribution of the
software without specific, written prior permission.

AT&T, Lucent and Bellcore disclaim all warranties with regard to
this software, including all implied warranties of
merchantability and fitness.  In no event shall AT&T, Lucent or
Bellcore be liable for any special, indirect or consequential
damages or any damages whatsoever resulting from loss of use,
data or profits, whether in an action of contract, negligence or
other tortious action, arising out of or in connection with the
use or performance of this software.
****************************************************************/

/*
** This malloc was obtained from ftp.cdrom.com at:
** .3/netlib/f2c/src/malloc.c
**
** minor modifications were made to produce 16 byte alignment
** rather than 8 byte.
**
** CEB 7/19/99
*/

#ifndef CRAY
#define STACKMIN 512
#define MINBLK (2*sizeof(struct mem) + 32)
/* was +16, changed to +32.  CEB */
#define F _malloc_free_
#define SBGULP 8192
#include "string.h"     /* for memcpy */

#ifdef KR_headers
#define Char char
#define Unsigned unsigned
#define Int /*int*/
#else
#define Char void
#define Unsigned size_t
#define Int int
#endif

typedef struct mem {
        struct mem *next;
        Unsigned len;
        double Padding; /* CEB */
        } mem;

static mem *F; /* added static CEB */

 Char *
#ifdef KR_headers
malloc(size)
        register Unsigned size;
#else
malloc(register Unsigned size)
#endif
{
        register mem *p, *q, *r, *s;
        unsigned register k, m;
        extern Char *sbrk(Int);
        char *top, *top1;

/*        size = (size+7) & ~7;*/
        size = (size+15) & ~15; /* CEB */
        r = (mem *) &F;
        for (p = F, q = 0; p; r = p, p = p->next) {
                if ((k = p->len) >= size && (!q || m > k)) {
                        m = k;
                        q = p;
                        s = r;
                        }
                }
        if (q) {
                if (q->len - size >= MINBLK) { /* split block */
                        p = (mem *) (((char *) (q+1)) + size);
                        p->next = q->next;
                        p->len = q->len - size - sizeof(mem);
                        s->next = p;
                        q->len = size;
                        }
                else
                        s->next = q->next;
                }
        else {
/*                top = (Char *)(((long)sbrk(0) + 7) & ~7);*/
                top = (Char *)(((long)sbrk(0) + 15) & ~15); /* CEB */
                if (F && (char *)(F+1) + F->len == top) {
                        q = F;
                        F = F->next;
                        }
                else
                        q = (mem *) top;
                top1 = (char *)(q+1) + size;
                if (sbrk((int)(top1-top+SBGULP)) == (Char *) -1)
                        return 0;
                r = (mem *)top1;
                r->len = SBGULP - sizeof(mem);
                r->next = F;
                F = r;
                q->len = size;
                }
        return (Char *) (q+1);
        }

 void
#ifdef KR_headers
free(f)
        Char *f;
#else
free(Char *f)
#endif
{
        mem *p, *q, *r;
        char *pn, *qn;

        if (!f) return;
        q = (mem *) ((char *)f - sizeof(mem));
        qn = (char *)f + q->len;
        for (p = F, r = (mem *) &F; ; r = p, p = p->next) {
                if (qn == (Char *) p) {
                        q->len += p->len + sizeof(mem);
                        p = p->next;
                        }
                pn = p ? ((char *) (p+1)) + p->len : 0;
                if (pn == (Char *) q) {
                        p->len += sizeof(mem) + q->len;
                        q->len = 0;
                        q->next = p;
                        r->next = p;
                        break;
                        }
                if (pn < (char *) q) {
                        r->next = q;
                        q->next = p;
                        break;
                        }
                }
        }

 Char *
#ifdef KR_headers
realloc(f, size)
        Char *f;
        Unsigned size;
#else
realloc(Char *f, Unsigned size)
#endif
{
        mem *p;
        Char *q, *f1;
        Unsigned s1;

        if (!f) return malloc(size);
        p = (mem *) ((char *)f - sizeof(mem));
        s1 = p->len;
        free(f);
        if (s1 > size)
/*                s1 = size + 7 & ~7;*/
                s1 = size + 15 & ~15; /* CEB */
        if (!p->len) {
                f1 = (Char *)(p->next + 1);
                memcpy(f1, f, s1);
                f = f1;
                }
        q = malloc(size);
        if (q && q != f)
                memcpy(q, f, s1);
        return q;
        }

/* The following (calloc) should really be in a separate file, */
/* but defining it here sometimes avoids confusion on systems */
/* that do not provide calloc in its own file. */

 Char *
#ifdef KR_headers
calloc(n, m) Unsigned m, n;
#else
calloc(Unsigned n, Unsigned m)
#endif
{
        Char *rv = malloc(n *= m);
        if (n && rv)
                memset(rv, 0, n);
        return rv;
        }
#endif


/* end of #ifdef USE_SPECIAL_DJGPP_MALLOC */
#endif

#ifdef ALIGN_MALLOC
void *
AlignedMalloc(size_t len)
{void *ptr;
ptr=malloc(len+300);
if (ptr != NULL)
  {unsigned long int Addr;
   void *NewPtr=ptr;
   void **vpp;
   unsigned char *cp;

   vpp=(void**)NewPtr;vpp++;NewPtr=(void*)vpp; /* allow space to save orig.*/
   Addr=(unsigned long int)NewPtr;
   while ((Addr % 256)!=0)
     {
      cp=(unsigned char*)NewPtr;
      cp++;
      NewPtr=(void*)cp;
      Addr=(unsigned long int)NewPtr;
     }
/* It's now aligned.  Save it. */
/*fprintf(stderr,"alloc: Orig pointer=%p Aligned pointer=%p\n",ptr,NewPtr);*/
   vpp=(void**)NewPtr;
   vpp--;
   *vpp=ptr;
   ptr=NewPtr;
   }
return ptr;
}

void
AlignedFree(void *NewPtr)
{void **vpp;void *ptr;
if (NewPtr==NULL) return;
vpp=(void**)NewPtr;
vpp--;
ptr=*vpp;
/*fprintf(stderr,"free : Orig pointer=%p Aligned pointer=%p\n",ptr,NewPtr);*/
free(ptr);
}

#else
void *
AlignedMalloc(size_t len)
{void *ptr;
ptr=malloc(len);
#if 0
if (ptr != NULL)
  {unsigned long int x=(unsigned long int)ptr;
   if ((x % sizeof(double))!=0)
     {
      fprintf(stderr,"Warning: AlignedMalloc has a pointer to %lu bytes that is\n",(ULINT)len);
      fprintf(stderr,"not aligned on a sizeof(double) boundary.  %p\n",ptr);
      fprintf(stderr,"This is not a fatal problem. Just a performance problem.\n");
     }
   }
#endif
return ptr;
}

void
AlignedFree(void *ptr)
{
if (ptr!=NULL) free(ptr);
}
#endif



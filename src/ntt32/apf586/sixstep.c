/*
This file is the 586 specific 'six step convolution' from Mikko Tommila's
APFloat v1.50.

It is distributed under his freeware license and is not my code nor is
it distributed under my license.  I make no claims about it other than
I've run it a couple of times and it appears to work for me.  If there
are bugs in it, then you should contact Mikko, not me, because I didn't
write it.

The full, unmodified code can be obtained from Mikko's APFloat package,
where ever it happens to be at the current time.

July 21, 1999
Carey Bloodworth


Below is Mikko's original license, so you can read it yourself:
---
LEGAL NOTICE



This program (the apfloat source code and documentation) is freeware. This
means that you can freely use, distribute, modify and compile it, but you
can't sell it or any part of it. Basically you can do anything with it, but
the program or any part of it will always be free. That is you can't charge
money or other valuables or services for it.

Although you can use this program freely, it would perhaps be considered to
be good manners to give the original author credit for his work, if this
program is ever used for anything useful or remarkable.

The author takes no responsibility whatsoever for any damage or harm that
could result from using this program. The program has been thoroughly
tested, so using it should be fairly safe. However, executing it as root is
perhaps not a very good idea.



Once more (a standard disclaimer):

THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK
AS TO THE QUALITY AND PERFORMANCE OF THE PRODUCT IS WITH YOU. SHOULD THE
PRODUCT PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

IN NO EVENT WILL MIKKO TOMMILA, THE AUTHOR OF THIS SOFTWARE, OR ANY OTHER
PARTY WHO MAY HAVE REDISTRIBUTED THE PRODUCT AS PERMITTED ABOVE, BE LIABLE
TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR
CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
PRODUCT (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED
INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE
PRODUCT TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER
PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.



Whew. Anyway, I have no money, so suing me makes no sense.

---
*/
/*
** I've hard wired these buffers.  I'm not sure what size they should
** be, but these should do for most situations.
*/
ModInt wtable[65536];
size_t ptable[65536];
ModInt SqrBuf1[65536];
ModInt SqrBuf2[65536];
ModInt MBuf1[65536];
int    MBuf2[65536];
/*
ModInt wtable[2048];
size_t ptable[2048];
ModInt SqrBuf1[256];
ModInt SqrBuf2[256];
ModInt MBuf1[2048];
int    MBuf2[2048];
*/
int Cacheblocksize=16;
int Cacheburstblocksize=4;
int Cachetreshold=32;

#define moveraw(d,s,n) memcpy((d),(s),(n)*sizeof(ModInt))


/*
** Seems to be a bit buggy for me.  The C version is a little slower
** but it works.
void swapraw (void *d, void *s, size_t n);
asm ("
    .globl _swapraw
.align 4
_swapraw:
    pushl %ebx
    pushl %ecx
    pushl %esi
    pushl %edi

    movl 20(%esp), %edi
    movl 24(%esp), %esi
    movl 28(%esp), %ecx

    jecxz swaprawend

    swaprawloop:
    movl (%esi), %eax
    movl (%edi), %ebx
    movl %eax, (%edi)
    movl %ebx, (%esi)
    addl $4, %esi
    addl $4, %edi
    decl %ecx
    jnz swaprawloop

    swaprawend:

    popl %edi
    popl %esi
    popl %ecx
    popl %ebx

    ret
");
*/
void
swapraw (ModInt *dd, ModInt *ss, int n)
{
    ModInt tmp;

    while (n--)
    {
        tmp = *dd;
        *(dd++) = *ss;
        *(ss++) = tmp;
    }
}

inline size_t min (size_t a, size_t b)
{
    return (a < b ? a : b);
}


size_t
permute (size_t n, size_t nn)
{
    size_t p = 1;

    do
       {
        p += p + (n & 1);
        n >>= 1;
       } while (p < nn);

    return p - nn;
}

void
scramble (ModInt data[], size_t nn)
{
    size_t i, j;
    ModInt tmp;

    for (i = 0; i < nn; i++)
       {
        j = permute (i, nn);
        if (j < i)
           {
            tmp = data[i];
            data[i] = data[j];
            data[j] = tmp;
           }
       }
}

void
initscrambletable (size_t ptable[], size_t nn)
{
    size_t i, j, k;

    for (i = k = 0; i < nn; i++)
       {
        j = permute (i, nn);
        if (j < i)
           {
            ptable[k] = i;
            ptable[k + 1] = j;
            k += 2;
           }
       }

    ptable[k] = ptable[k + 1] = 0;      // Table end marker
}

void tablescramble (ModInt data[], size_t ptable[])
{
    size_t i, j, k = 0;
    ModInt tmp;

    while ((i = ptable[k]) != (j = ptable[k + 1]))
       {
        tmp = data[i];
        data[i] = data[j];
        data[j] = tmp;
        k += 2;
       }
}

void tablefnt (ModInt data[], ModInt wtable[], size_t ptable[], size_t nn, int s);
void itablefnt (ModInt data[], ModInt wtable[], size_t ptable[], size_t nn, int s);

// Optimized fnt using table lookup for the powers of W

asm ("
    .globl _tablefnt
    .globl _itablefnt

.align 5
_tablefnt:
    subl $32, %esp
    pushl %ebx
    pushl %ecx
    pushl %edx
    pushl %esi
    pushl %edi
    pushl %ebp

    movl 60(%esp), %eax
    movl 72(%esp), %ebx
    shll $2, %ebx
    addl %ebx, %eax
    movl %ebx, 44(%esp)
    movl %eax, 52(%esp)
    movl $4, 32(%esp)

    movl _modulus, %ecx

    fldl _imodulus
    flds _chopper64
    fxch %st(1)
    fildl _modulus
    fxch %st(1)

    tablefntloop:
    movl 44(%esp), %eax
    shrl $1, %eax
    cmpl $4, %eax
    jb tablefntend
    movl %eax, 44(%esp)
    movl %eax, %ebx
    addl 60(%esp), %eax
    addl %ebx, %ebx
    movl %eax, 48(%esp)
    movl %ebx, 40(%esp)


    movl 60(%esp), %esi

    tablefntfor1:
    movl %esi, %edi
    movl 44(%esp), %ebx
    addl %ebx, %edi
    movl 52(%esp), %eax

    cmpl %eax, %esi
    jae tablefntfor1end

    movl (%esi), %eax
    movl (%edi), %ebx

    addl %eax, %ebx
    addl %eax, %eax
    subl %ebx, %eax
    sbbl %edx, %edx
    subl %ecx, %ebx
    sbbl %ebp, %ebp
    andl %ecx, %edx
    andl %ecx, %ebp
    addl %edx, %eax
    addl %ebp, %ebx

    movl %eax, (%edi)
    movl %ebx, (%esi)

    movl 40(%esp), %eax
    addl %eax, %esi
    jmp tablefntfor1

    tablefntfor1end:


    movl 64(%esp), %ebp
    movl 32(%esp), %eax
    addl %eax, %ebp
    movl %ebp, 28(%esp)

    movl 60(%esp), %eax
    addl $4, %eax

    tablefntfor2:
    cmpl 48(%esp), %eax
    jae tablefntfor2end

    movl %eax, 36(%esp)
    movl %eax, %esi

    movl %esi, %edi
    movl 44(%esp), %ebx
    addl %ebx, %edi
    movl 52(%esp), %eax

    cmpl %eax, %esi
    jae tablefntfor3

    movl (%esi), %eax
    movl (%edi), %ebx

    addl %eax, %ebx
    addl %eax, %eax
    subl %ebx, %eax
    sbbl %edx, %edx
    subl %ecx, %ebx
    sbbl %ebp, %ebp
    andl %ecx, %edx
    andl %ecx, %ebp
    addl %edx, %eax
    addl %ebp, %ebx

    movl %eax, 24(%esp)
    movl %ebx, (%esi)

    movl 40(%esp), %eax
    movl 28(%esp), %ebp
    addl %eax, %esi

    tablefntfor3:

    fildl 24(%esp)
    fildl (%ebp)

    movl 52(%esp), %eax

    cmpl %eax, %esi
    jae tablefntfor3end

    fmulp %st, %st(1)

    movl %esi, %edx
    movl 44(%esp), %ebx
    addl %ebx, %edx

    fld %st
    fmul %st(2), %st

    pushl %edx

    movl (%esi), %eax
    movl (%edx), %ebx

    fadd %st(4), %st

    addl %eax, %ebx
    addl %eax, %eax
    subl %ebx, %eax

    fsub %st(4), %st

    sbbl %edx, %edx
    subl %ecx, %ebx
    sbbl %ebp, %ebp
    andl %ecx, %edx

    fmul %st(3), %st

    andl %ecx, %ebp
    addl %edx, %eax
    addl %ebp, %ebx
    movl 32(%esp), %ebp

    fsubrp %st, %st(1)

    movl %eax, 28(%esp)
    movl 44(%esp), %eax
    movl %ebx, (%esi)
    addl %eax, %esi

    fistpl (%edi)

    popl %edi
    jmp tablefntfor3

    tablefntfor3end:

    fmulp %st, %st(1)
    fld %st
    fmul %st(2), %st
    fadd %st(4), %st
    fsub %st(4), %st
    fmul %st(3), %st

    movl 28(%esp), %ebp
    movl 32(%esp), %eax

    fsubrp %st, %st(1)

    addl %eax, %ebp
    movl %ebp, 28(%esp)
    movl 36(%esp), %eax
    addl $4, %eax

    fistpl (%edi)

    jmp tablefntfor2

    tablefntfor2end:


    shll $1, 32(%esp)
    jmp tablefntloop


    tablefntend:

    cmpl $0, 76(%esp)
    je tablefntnos

    movl 68(%esp), %ebp
    movl 60(%esp), %ebx

    tablefntscramble:
    movl (%ebp), %esi
    movl 4(%ebp), %edi
    cmpl %esi, %edi
    je tablefntnos
    movl (%ebx, %esi, 4), %eax
    movl (%ebx, %edi, 4), %edx
    addl $8, %ebp
    movl %eax, (%ebx, %edi, 4)
    movl %edx, (%ebx, %esi, 4)
    jmp tablefntscramble

    tablefntnos:

    fstp %st
    fstp %st
    fstp %st

    popl %ebp
    popl %edi
    popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    addl $32, %esp

    ret

.align 5
_itablefnt:
    subl $32, %esp
    pushl %ebx
    pushl %ecx
    pushl %edx
    pushl %esi
    pushl %edi
    pushl %ebp

    cmpl $0, 76(%esp)
    je itablefntnos

    movl 68(%esp), %ebp
    movl 60(%esp), %ebx

    itablefntscramble:
    movl (%ebp), %esi
    movl 4(%ebp), %edi
    cmpl %esi, %edi
    je itablefntnos
    movl (%ebx, %esi, 4), %eax
    movl (%ebx, %edi, 4), %edx
    addl $8, %ebp
    movl %eax, (%ebx, %edi, 4)
    movl %edx, (%ebx, %esi, 4)
    jmp itablefntscramble

    itablefntnos:


    movl 60(%esp), %eax
    movl 72(%esp), %ebx
    shll $2, %ebx
    addl %ebx, %eax
    movl %ebx, 32(%esp)
    movl %eax, 52(%esp)
    movl $4, 44(%esp)

    movl _modulus, %ecx

    fldl _imodulus
    flds _chopper64
    fxch %st(1)
    fildl _modulus
    fxch %st(1)

    itablefntloop:
    movl 32(%esp), %eax
    shrl $1, %eax
    cmpl $4, %eax
    jb itablefntend
    movl %eax, 32(%esp)
    movl 44(%esp), %eax
    movl %eax, %ebx
    addl 60(%esp), %eax
    movl %eax, 48(%esp)
    addl %ebx, %ebx
    movl %ebx, 40(%esp)


    movl 60(%esp), %esi

    itablefntfor1:
    movl %esi, %edi
    movl 44(%esp), %ebx
    addl %ebx, %edi
    movl 52(%esp), %eax

    cmpl %eax, %esi
    jae itablefntfor1end

    movl (%esi), %eax
    movl (%edi), %ebx

    addl %eax, %ebx
    addl %eax, %eax
    subl %ebx, %eax
    sbbl %edx, %edx
    subl %ecx, %ebx
    sbbl %ebp, %ebp
    andl %ecx, %edx
    andl %ecx, %ebp
    addl %edx, %eax
    addl %ebp, %ebx

    movl %eax, (%edi)
    movl %ebx, (%esi)

    movl 40(%esp), %eax
    addl %eax, %esi
    jmp itablefntfor1

    itablefntfor1end:


    movl 64(%esp), %ebp
    movl 32(%esp), %eax
    addl %eax, %ebp
    movl %ebp, 28(%esp)

    movl 60(%esp), %eax
    addl $4, %eax

    itablefntfor2:
    cmpl 48(%esp), %eax
    jae itablefntfor2end

    movl %eax, 36(%esp)
    movl %eax, %esi

    movl 28(%esp), %ebp
    movl %esi, %edi
    movl 44(%esp), %ebx
    addl %ebx, %edi

    fildl (%ebp)
    fildl (%edi)
    fmulp %st, %st(1)
    fld %st
    fmul %st(2), %st
    fadd %st(4), %st
    fsub %st(4), %st
    fmul %st(3), %st
    fsubrp %st, %st(1)

    movl 40(%esp), %ebx
    movl %edi, %edx
    addl %ebx, %edx
    movl 28(%esp), %ebp

    fistpl 24(%esp)

    itablefntfor3:

    fildl (%ebp)
    fildl (%edx)

    movl 52(%esp), %ebx
    movl 40(%esp), %eax
    subl %eax, %ebx

    cmpl %ebx, %esi
    jae itablefntfor3end

    fmulp %st, %st(1)

    movl %esi, %edi
    movl 44(%esp), %ebx
    addl %ebx, %edi

    fld %st
    fmul %st(2), %st

    movl (%esi), %eax
    movl 24(%esp), %ebx
    addl %eax, %ebx
    addl %eax, %eax

    fadd %st(4), %st

    subl %ebx, %eax
    sbbl %edx, %edx
    subl %ecx, %ebx

    fsub %st(4), %st

    sbbl %ebp, %ebp
    andl %ecx, %edx
    andl %ecx, %ebp
    addl %edx, %eax

    fmul %st(3), %st

    addl %ebp, %ebx
    movl %eax, (%edi)
    movl %ebx, (%esi)
    movl 40(%esp), %eax

    fsubrp %st, %st(1)

    addl %eax, %esi
    movl %edi, %edx
    addl %eax, %eax
    addl %eax, %edx
    movl 28(%esp), %ebp

    fistpl 24(%esp)

    jmp itablefntfor3

    itablefntfor3end:

    fstp %st
    fstp %st

    movl %esi, %edi
    movl 44(%esp), %ebx
    addl %ebx, %edi

    movl (%esi), %eax
    movl 24(%esp), %ebx

    addl %eax, %ebx
    addl %eax, %eax
    subl %ebx, %eax
    sbbl %edx, %edx
    subl %ecx, %ebx
    sbbl %ebp, %ebp
    andl %ecx, %edx
    andl %ecx, %ebp
    addl %edx, %eax
    addl %ebp, %ebx

    movl %eax, (%edi)
    movl %ebx, (%esi)

    movl 28(%esp), %ebp
    movl 32(%esp), %eax
    addl %eax, %ebp
    movl %ebp, 28(%esp)

    movl 36(%esp), %eax
    addl $4, %eax
    jmp itablefntfor2

    itablefntfor2end:


    shll $1, 44(%esp)
    jmp itablefntloop


    itablefntend:

    fstp %st
    fstp %st
    fstp %st

    popl %ebp
    popl %edi
    popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    addl $32, %esp

    ret
");

void moveblock (ModInt d[], ModInt s[], size_t b, size_t nd, size_t ns);
void transposeblock (ModInt data[], size_t b, size_t nn);
void transpose2blocks (ModInt data1[], ModInt data2[], size_t b, size_t nn);

// Move a b x b block from s to d, d width = nd , s width = ns
void
moveblock (ModInt d[], ModInt s[], size_t b, size_t nd, size_t ns)
{
    size_t j;

    for (j = 0; j < b; j++)
       {
        memcpy(d,s,b*sizeof(ModInt));
        d += nd;
        s += ns;
       }
}

asm ("
    .globl _transposeblock
.align 4
_transposeblock:
    pushl %ebx
    pushl %ecx
    pushl %edx
    pushl %esi
    pushl %edi
    pushl %ebp

    movl 36(%esp), %edx
    movl 32(%esp), %ebp

    shll $2, %edx
    decl %ebp
    jz transposeblockend

    transposeblockloop:

    movl 28(%esp), %edi
    movl 28(%esp), %esi
    movl %ebp, %ecx

    addl %edx, %edi
    addl $4, %esi

    leal 4(%edi), %eax
    movl %eax, 28(%esp)

    transposeblockline:
    movl (%esi), %eax
    movl (%edi), %ebx
    movl %eax, (%edi)
    movl %ebx, (%esi)
    addl $4, %esi
    addl %edx, %edi
    decl %ecx
    jnz transposeblockline

    decl %ebp
    jnz transposeblockloop

    transposeblockend:

    popl %ebp
    popl %edi
    popl %esi
    popl %edx
    popl %ecx
    popl %ebx

    ret
");

asm ("
    .globl _transpose2blocks
.align 4
_transpose2blocks:
    pushl %ebx
    pushl %ecx
    pushl %edx
    pushl %esi
    pushl %edi
    pushl %ebp

    movl 40(%esp), %edx

    shll $2, %edx

    movl 36(%esp), %ebp

    transpose2blocksloop:
    movl 28(%esp), %edi
    movl 32(%esp), %esi
    movl 36(%esp), %ecx

    addl $4, 28(%esp)
    addl %edx, 32(%esp)

    transpose2blocksline:
    movl (%esi), %eax
    movl (%edi), %ebx
    movl %eax, (%edi)
    movl %ebx, (%esi)
    addl $4, %esi
    addl %edx, %edi
    decl %ecx
    jnz transpose2blocksline

    decl %ebp
    jnz transpose2blocksloop

    popl %ebp
    popl %edi
    popl %esi
    popl %edx
    popl %ecx
    popl %ebx

    ret
");

// Transpose a _square_ n1 x n1 block of n1 x n2 matrix in b x b blocks
void
transposesquare (ModInt data[], size_t n1, size_t n2)
{
    size_t j, k, b;
    ModInt *p1, *p2;

    if (n1 <= Cacheburstblocksize || n1 <= Cacheblocksize)
        transposeblock (data, n1, n2);
    else if (n1 * n2 <= Cachetreshold)
       {
        b = Cacheburstblocksize;

        for (j = 0, p1 = data; j < n1; j += b, p1 += b * n2)
           {
            transposeblock (p1 + j, b, n2);
            for (k = j + b, p2 = data + k * n2 + j; k < n1; k += b, p2 += b * n2)
                transpose2blocks (p1 + k, p2, b, n2);
           }
       }
    else
       {ModInt *tmp1,*tmp2;
        b = Cacheblocksize;

//        ModInt *tmp1 = new ModInt[b * b];
//        ModInt *tmp2 = new ModInt[b * b];
        tmp1=SqrBuf1;
        tmp2=SqrBuf2;

        for (j = 0, p1 = data; j < n1; j += b, p1 += b * n2)
           {
            moveblock (tmp1, p1 + j, b, b, n2);
            transposeblock (tmp1, b, b);
            moveblock (p1 + j, tmp1, b, n2, b);

            for (k = j + b, p2 = data + k * n2 + j; k < n1; k += b, p2 += b * n2)
               {
                moveblock (tmp1, p1 + k, b, b, n2);
                transposeblock (tmp1, b, b);

                moveblock (tmp2, p2, b, b, n2);
                transposeblock (tmp2, b, b);

                moveblock (p1 + k, tmp2, b, n2, b);
                moveblock (p2, tmp1, b, n2, b);
               }
           }

//        delete[] tmp2;
//        delete[] tmp1;
       }
}

// Permute the rows of matrix to correct order, must be n2 = 2*n1
void
permutewidetotall (ModInt data[], size_t n1, size_t n2)
{
    size_t j, o, m;
    ModInt *d;int *r;

    if (n2 < 4) return;

//    ModInt *d = new ModInt[n1];
//    int *r = new int[n2];
    d=MBuf1;
    r=MBuf2;

    for (j = 0; j < n2; j++) r[j] = 0;

    j = 1;
    do
       {
        o = m = j;

        moveraw (d, data + n1 * m, n1);

        r[m] = 1;

        m = (m < n1 ? 2 * m : 2 * (m - n1) + 1);

        while (m != j)
           {
            r[m] = 1;

            moveraw (data + n1 * o, data + n1 * m, n1);

            o = m;
            m = (m < n1 ? 2 * m : 2 * (m - n1) + 1);
           };

        moveraw (data + n1 * o, d, n1);

        while (r[j]) j++;
       } while (j < n2 - 1);

//    delete[] r;
//    delete[] d;
}

// Permute the rows of matrix to correct order, must be n1 = 2*n2
void
permutetalltowide (ModInt data[], size_t n1, size_t n2)
{
    size_t j, o, m;
    ModInt *d;int *r;

    if (n1 < 4) return;

//    ModInt *d = new ModInt[n2];
//    int *r = new int[n1];
    d=MBuf1;
    r=MBuf2;

    for (j = 0; j < n1; j++) r[j] = 0;

    j = 1;
    do
       {
        o = m = j;

        moveraw (d, data + n2 * m, n2);

        r[m] = 1;

        m = (m & 1 ? m / 2 + n2 : m / 2);

        while (m != j)
           {
            r[m] = 1;

            moveraw (data + n2 * o, data + n2 * m, n2);

            o = m;
            m = (m & 1 ? m / 2 + n2 : m / 2);
           };

        moveraw (data + n2 * o, d, n2);

        while (r[j]) j++;
       } while (j < n1 - 1);

//    delete[] r;
//    delete[] d;
}

// Transpose a n1 x n2 matrix. Must be n1 = n2, n1 = 2*n2 or n2 = 2*n1
void
transpose (ModInt data[], size_t n1, size_t n2)
{
    if (n1 == n2)
       {
        // simply transpose

        transposesquare (data, n1, n1);
       }
    else if (n2 == 2 * n1)
       {
        // first transpose two n1 x n1 blocks
        transposesquare (data, n1, n2);
        transposesquare (data + n1, n1, n2);

        // then permute the rows to correct order
        permutewidetotall (data, n1, n2);
       }
    else if (n1 == 2 * n2)
       {
        // first permute the rows to correct order
        permutetalltowide (data, n1, n2);

        // then transpose two n2 x n2 blocks
        transposesquare (data, n2, n1);
        transposesquare (data + n2, n2, n1);
       }
}

#if 0
inline void modmul3 (modint &a1, modint &b1, modint &a2, modint &b2, modint &a3, modint &b3)
{
    asm ("pushl %3; pushl %4;
          pushl %5; pushl %6;
          pushl %7; pushl %8;
          fildl 20(%%esp); fildl 16(%%esp);
          fildl 12(%%esp); fxch %%st(2);
          fmulp %%st, %%st(1); fxch %%st(1);
          fildl 8(%%esp); fildl 4(%%esp); fxch %%st(2);
          fmulp %%st, %%st(1); fxch %%st(1);
          fildl (%%esp); addl $12, %%esp;
          fmulp %%st, %%st(1);
          fld %%st(2); fmull _imodulus;
          fld %%st(2); fmull _imodulus;
          fld %%st(2); fmull _imodulus; fxch %%st(2);
          fadds _chopper64; fxch %%st(1);
          fadds _chopper64; fxch %%st(2);
          fadds _chopper64; fxch %%st(1);
          fsubs _chopper64; fxch %%st(2);
          fsubs _chopper64; fxch %%st(1);
          fsubs _chopper64; fxch %%st(2);
          fmull _dmodulus; fxch %%st(1);
          fmull _dmodulus; fxch %%st(1);
          fsubrp %%st, %%st(5); fxch %%st(1);
          fmull _dmodulus; fxch %%st(1);
          fsubrp %%st, %%st(3);
          fsubrp %%st, %%st(1); fxch %%st(2);
          fistpl (%%esp);
          fistpl 4(%%esp);
          fistpl 8(%%esp);
          popl %0;
          popl %1;
          popl %2"
                            : "=rm" (a1), "=rm" (a2), "=rm" (a3)
                            : "0" (a1), "rm" (b1), "1" (a2), "rm" (b2), "2" (a3), "rm" (b3)
                            : "cc");
}

// Partially parallel two modular multiplications using the FPU
inline void modmul2 (modint &a1, modint &b1, modint &a2, modint &b2)
{
    asm ("pushl %2; pushl %3;
          pushl %4; pushl %5;
          fildl 12(%%esp); fildl 8(%%esp);
          fildl 4(%%esp); fxch %%st(2);
          fmulp %%st, %%st(1); fxch %%st(1);
          fildl (%%esp); addl $8, %%esp;
          fmulp %%st, %%st(1);
          fld %%st(1); fmull _imodulus;
          fld %%st(1); fmull _imodulus; fxch %%st(1);
          fadds _chopper64; fxch %%st(1);
          fadds _chopper64; fxch %%st(1);
          fsubs _chopper64; fxch %%st(1);
          fsubs _chopper64; fxch %%st(1);
          fmull _dmodulus; fxch %%st(1);
          fmull _dmodulus; fxch %%st(1);
          fsubrp %%st, %%st(3);
          fsubrp %%st, %%st(1); fxch %%st(1);
          fistpl (%%esp);
          fistpl 4(%%esp);
          popl %0;
          popl %1"
                            : "=rm" (a1), "=rm" (a2)
                            : "0" (a1), "rm" (b1), "1" (a2), "rm" (b2)
                            : "cc");
}
#endif

// The "six-step" fnt, but doesn't transpose or scramble (for convolution only)
void
FSixStep(ModInt data[], ModInt pr, ModInt nr,int isign, size_t nn)
{
    size_t n1, n2, j, k;
    ModInt w, tmp, tmp2, *p1, *p2;

    if (nn < 2) return;

    for (n1 = 1, n2 = 0; n1 < nn; n1 <<= 1, n2++);
    n1 = n2 >> 1;
    n2 -= n1;

    n1 = 1 << n1;
    n2 = 1 << n2;

    // n2 >= n1

//    ModInt *wtable = new ModInt[n2];
//    size_t *ptable = new size_t[n1];

//    if (isign > 0) w=NthRoot;
//    else           w=NthRoot1;
    w=nr;

    transpose (data, n1, n2);

    tmp = ModPow (w, nn / n1);
    tmp2 = 1;
    for (k = 0; k < n1; k++)
      {
        wtable[k] = tmp2;
        tmp2 = ModMul(tmp2,tmp);
      }

    initscrambletable (ptable, n1);

    for (k = 0, p1 = data; k < n2; k++, p1 += n1)
        tablefnt (p1, wtable, ptable, n1,1);

    transpose (data, n2, n1);

    tmp = w;
    for (j = 1, p1 = data + n2; j < n1; j++, p1 += n2)
      {
        tmp2 = ModPow (tmp, j);
        p1[j] = ModMul(p1[j],tmp2);
//        tmp2 =ModMul(tmp2,tmp);
        for (k = j + 1, p2 = p1 + n2 + j; k < n1; k++, p2 += n2)
           {
            tmp2 =ModMul(tmp2,tmp);
            p1[k]=ModMul(p1[k],tmp2);
            *p2=ModMul(*p2,tmp2);
//            modmul3 (tmp2, tmp, p1[k], tmp2, *p2, tmp2);
           }
        for (; k < n2; k++)
           {
            tmp2 =ModMul(tmp2,tmp);
            p1[k]=ModMul(p1[k],tmp2);
//            modmul2 (tmp2, tmp, p1[k], tmp2);
           }
        tmp=ModMul(tmp,w);
       }

    if (n2 != n1)
       {
        // n2 = 2 * n1
        for (k = n1; k--;)
            wtable[2 * k] = wtable[k];
        tmp = ModPow (w, nn / n2);
        tmp2 = tmp;
        tmp = ModMul(tmp,tmp);
        for (k = 1; k < n2; k += 2)
           {
            wtable[k] = tmp2;
            tmp2=ModMul(tmp2,tmp);
           }
       }

    for (k = 0, p1 = data; k < n1; k++, p1 += n2)
        tablefnt (p1, wtable, 0, n2, 0);

//    delete[] ptable;
//    delete[] wtable;
}

void
RSixStep(ModInt data[], ModInt pr, ModInt nr, int isign, size_t nn)
{
    size_t n1, n2, j, k;
    ModInt w, tmp, tmp2, *p1, *p2;

    if (nn < 2) return;

    for (n1 = 1, n2 = 0; n1 < nn; n1 <<= 1, n2++);
    n1 = n2 >> 1;
    n2 -= n1;

    n1 = 1 << n1;
    n2 = 1 << n2;

    // n2 >= n1

//    ModInt *wtable = new ModInt[n2];
//    size_t *ptable = new size_t[n1];

//    if (isign > 0) w=NthRoot;
//    else           w=NthRoot1;
    w=nr;

    tmp = ModPow (w, nn / n2);
    tmp2 = 1;
    for (k = 0; k < n2; k++)
       {
        wtable[k] = tmp2;
        tmp2 = ModMul(tmp2,tmp);
       }

    for (k = 0, p1 = data; k < n1; k++, p1 += n2)
        itablefnt (p1, wtable, 0, n2, 0);

    tmp = 1;
    for (j = 0, p1 = data; j < n1; j++, p1 += n2)
       {
        tmp2 = ModPow (tmp, j);
        p1[j]=ModMul(p1[j],tmp2);
//        tmp2 =ModMul(tmp2,tmp);
        for (k = j + 1, p2 = p1 + n2 + j; k < n1; k++, p2 += n2)
           {
            tmp2 =ModMul(tmp2,tmp);
            p1[k]=ModMul(p1[k],tmp2);
            *p2  =ModMul(*p2,tmp2);
//            modmul3 (tmp2, tmp, p1[k], tmp2, *p2, tmp2);
           }
        for (; k < n2; k++)
           {
            tmp2=ModMul(tmp2,tmp);
            p1[k]=ModMul(p1[k],tmp2);
//            modmul2 (tmp2, tmp, p1[k], tmp2);
           }
        tmp=ModMul(tmp,w);
       }

    transpose (data, n1, n2);

    if (n2 != n1)
        // n2 = 2 * n1
        for (k = 0; k < n1; k++)
            wtable[k] = wtable[2 * k];

    initscrambletable (ptable, n1);

    for (k = 0, p1 = data; k < n2; k++, p1 += n1)
        itablefnt (p1, wtable, ptable, n1, 1);

    transpose (data, n2, n1);

//    delete[] ptable;
//    delete[] wtable;
}



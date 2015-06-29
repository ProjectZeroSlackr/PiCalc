#include "pi.h"
#include "fft.h"
#include "cache.h"

#ifdef LONG_DOUBLE_FFT
#error This external FFT does not support long double.  Only double.
#endif

/*
    Copyright(C) 1999 Takuya OOURA
    You may use, copy, modify this code for any purpose and
    without fee.

Fast Fourier Transform
    dimension   :one
    data length :power of 2
    decimation  :frequency
    radix       :4, 2
    data        :inplace
    table       :not use
functions
    rdft: Real Discrete Fourier Transform
function prototypes
    void rdft(int, int, double *);

-------- Real DFT / Inverse of Real DFT --------
    [definition]
        <case1> RDFT
            R[k] = sum_j=0^n-1 a[j]*cos(2*pi*j*k/n), 0<=k<=n/2
            I[k] = sum_j=0^n-1 a[j]*sin(2*pi*j*k/n), 0<k<n/2
        <case2> IRDFT (excluding scale)
            a[k] = (R[0] + R[n/2]*cos(pi*k))/2 +
                   sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) +
                   sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n
    [usage]
        <case1>
            rdft(n, 1, a);
        <case2>
            rdft(n, -1, a);
    [parameters]
        n              :data length (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (double *)
                        <case1>
                            output data
                                a[2*k] = R[k], 0<=k<n/2
                                a[2*k+1] = I[k], 0<k<n/2
                                a[1] = R[n/2]
                        <case2>
                            input data
                                a[2*j] = R[j], 0<=j<n/2
                                a[2*j+1] = I[j], 0<j<n/2
                                a[1] = R[n/2]
    [remark]
        Inverse of
            rdft(n, 1, a);
        is
            rdft(n, -1, a);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .

*/

/*
** This FFT was modified a bit to work with my pi program. CEB.
*/

#include <stdlib.h>

#include <math.h>
#ifndef M_PI_2
#define M_PI_2      1.570796326794896619231321691639751442098584699687
#endif

static void bitrv2(int n, double *a)
{
    int j0, k0, j1, k1, l, m, i, j, k;
    double xr, xi;

    l = n >> 1;
    m = 2;
    while (m < l) {
        l >>= 1;
        m <<= 1;
    }
    if (m > l) {
        j0 = 0;
        for (k0 = 2; k0 < m; k0 += 2) {
            for (i = n >> 1; i > (j0 ^= i); i >>= 1);
            k = k0;
            for (j = j0; j < j0 + k0; j += 2) {
                xr = a[j];
                xi = a[j + 1];
                a[j] = a[k];
                a[j + 1] = a[k + 1];
                a[k] = xr;
                a[k + 1] = xi;
                for (i = n >> 1; i > (k ^= i); i >>= 1);
            }
        }
    } else {
        j0 = 0;
        for (k0 = 2; k0 < m; k0 += 2) {
            for (i = n >> 1; i > (j0 ^= i); i >>= 1);
            k = k0;
            for (j = j0; j < j0 + k0; j += 2) {
                xr = a[j];
                xi = a[j + 1];
                a[j] = a[k];
                a[j + 1] = a[k + 1];
                a[k] = xr;
                a[k + 1] = xi;
                j1 = j + m;
                k1 = k + m;
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
                for (i = n >> 1; i > (k ^= i); i >>= 1);
            }
        }
    }
}

static void cft1st(int n, double *a)
{
    int j, kj, kr;
    double ew, wk0r, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    x0r = a[0] + a[2];
    x0i = a[1] + a[3];
    x1r = a[0] - a[2];
    x1i = a[1] - a[3];
    x2r = a[4] + a[6];
    x2i = a[5] + a[7];
    x3r = a[4] - a[6];
    x3i = a[5] - a[7];
    a[0] = x0r + x2r;
    a[1] = x0i + x2i;
    a[4] = x0r - x2r;
    a[5] = x0i - x2i;
    a[2] = x1r - x3i;
    a[3] = x1i + x3r;
    a[6] = x1r + x3i;
    a[7] = x1i - x3r;
    ew = M_PI_2 / n;
    wk0r = cos(ew * (n >> 1));
    x0r = a[8] + a[10];
    x0i = a[9] + a[11];
    x1r = a[8] - a[10];
    x1i = a[9] - a[11];
    x2r = a[12] + a[14];
    x2i = a[13] + a[15];
    x3r = a[12] - a[14];
    x3i = a[13] - a[15];
    a[8] = x0r + x2r;
    a[9] = x0i + x2i;
    a[12] = x2i - x0i;
    a[13] = x0r - x2r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[10] = wk0r * (x0r - x0i);
    a[11] = wk0r * (x0r + x0i);
    x0r = x3i + x1r;
    x0i = x3r - x1i;
    a[14] = wk0r * (x0i - x0r);
    a[15] = wk0r * (x0i + x0r);
    kr = 0;
    for (j = 16; j < n; j += 16) {
        for (kj = n >> 2; kj > (kr ^= kj); kj >>= 1);
        wk1r = cos(ew * kr);
        wk1i = sin(ew * kr);
        wk2r = 1 - 2 * wk1i * wk1i;
        wk2i = 2 * wk1i * wk1r;
        wk3r = wk1r - 2 * wk2i * wk1i;
        wk3i = 2 * wk2i * wk1r - wk1i;
        x0r = a[j] + a[j + 2];
        x0i = a[j + 1] + a[j + 3];
        x1r = a[j] - a[j + 2];
        x1i = a[j + 1] - a[j + 3];
        x2r = a[j + 4] + a[j + 6];
        x2i = a[j + 5] + a[j + 7];
        x3r = a[j + 4] - a[j + 6];
        x3i = a[j + 5] - a[j + 7];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        a[j + 4] = wk2r * x0r - wk2i * x0i;
        a[j + 5] = wk2r * x0i + wk2i * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j + 2] = wk1r * x0r - wk1i * x0i;
        a[j + 3] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j + 6] = wk3r * x0r - wk3i * x0i;
        a[j + 7] = wk3r * x0i + wk3i * x0r;
        x0r = wk0r * (wk1r - wk1i);
        wk1i = wk0r * (wk1r + wk1i);
        wk1r = x0r;
        wk3r = wk1r - 2 * wk2r * wk1i;
        wk3i = 2 * wk2r * wk1r - wk1i;
        x0r = a[j + 8] + a[j + 10];
        x0i = a[j + 9] + a[j + 11];
        x1r = a[j + 8] - a[j + 10];
        x1i = a[j + 9] - a[j + 11];
        x2r = a[j + 12] + a[j + 14];
        x2i = a[j + 13] + a[j + 15];
        x3r = a[j + 12] - a[j + 14];
        x3i = a[j + 13] - a[j + 15];
        a[j + 8] = x0r + x2r;
        a[j + 9] = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        a[j + 12] = -wk2i * x0r - wk2r * x0i;
        a[j + 13] = -wk2i * x0i + wk2r * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j + 10] = wk1r * x0r - wk1i * x0i;
        a[j + 11] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j + 14] = wk3r * x0r - wk3i * x0i;
        a[j + 15] = wk3r * x0i + wk3i * x0r;
    }
}

static void cftmdl(int n, int l, double *a)
{
    int j, j1, j2, j3, k, kj, kr, m, m2;
    double ew, wk0r, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    m = l << 2;
    for (j = 0; j < l; j += 2) {
        j1 = j + l;
        j2 = j1 + l;
        j3 = j2 + l;
        x0r = a[j] + a[j1];
        x0i = a[j + 1] + a[j1 + 1];
        x1r = a[j] - a[j1];
        x1i = a[j + 1] - a[j1 + 1];
        x2r = a[j2] + a[j3];
        x2i = a[j2 + 1] + a[j3 + 1];
        x3r = a[j2] - a[j3];
        x3i = a[j2 + 1] - a[j3 + 1];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j2] = x0r - x2r;
        a[j2 + 1] = x0i - x2i;
        a[j1] = x1r - x3i;
        a[j1 + 1] = x1i + x3r;
        a[j3] = x1r + x3i;
        a[j3 + 1] = x1i - x3r;
    }
    ew = M_PI_2 / n;
    wk0r = cos(ew * (n >> 1));
    for (j = m; j < l + m; j += 2) {
        j1 = j + l;
        j2 = j1 + l;
        j3 = j2 + l;
        x0r = a[j] + a[j1];
        x0i = a[j + 1] + a[j1 + 1];
        x1r = a[j] - a[j1];
        x1i = a[j + 1] - a[j1 + 1];
        x2r = a[j2] + a[j3];
        x2i = a[j2 + 1] + a[j3 + 1];
        x3r = a[j2] - a[j3];
        x3i = a[j2 + 1] - a[j3 + 1];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j2] = x2i - x0i;
        a[j2 + 1] = x0r - x2r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j1] = wk0r * (x0r - x0i);
        a[j1 + 1] = wk0r * (x0r + x0i);
        x0r = x3i + x1r;
        x0i = x3r - x1i;
        a[j3] = wk0r * (x0i - x0r);
        a[j3 + 1] = wk0r * (x0i + x0r);
    }
    kr = 0;
    m2 = m << 1;
    for (k = m2; k < n; k += m2) {
        for (kj = n >> 2; kj > (kr ^= kj); kj >>= 1);
        wk1r = cos(ew * kr);
        wk1i = sin(ew * kr);
        wk2r = 1 - 2 * wk1i * wk1i;
        wk2i = 2 * wk1i * wk1r;
        wk3r = wk1r - 2 * wk2i * wk1i;
        wk3i = 2 * wk2i * wk1r - wk1i;
        for (j = k; j < l + k; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j2] = wk2r * x0r - wk2i * x0i;
            a[j2 + 1] = wk2r * x0i + wk2i * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j1] = wk1r * x0r - wk1i * x0i;
            a[j1 + 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = wk3r * x0r - wk3i * x0i;
            a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        }
        x0r = wk0r * (wk1r - wk1i);
        wk1i = wk0r * (wk1r + wk1i);
        wk1r = x0r;
        wk3r = wk1r - 2 * wk2r * wk1i;
        wk3i = 2 * wk2r * wk1r - wk1i;
        for (j = k + m; j < l + (k + m); j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j2] = -wk2i * x0r - wk2r * x0i;
            a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j1] = wk1r * x0r - wk1i * x0i;
            a[j1 + 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = wk3r * x0r - wk3i * x0i;
            a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        }
    }
}


static void cftfsub(int n, double *a)
{
    int j, j1, j2, j3, l;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 2;
    if (n > 8) {
        cft1st(n, a);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a);
            l <<= 2;
        }
    }
    if ((l << 2) == n) {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            a[j2] = x0r - x2r;
            a[j2 + 1] = x0i - x2i;
            a[j1] = x1r - x3i;
            a[j1 + 1] = x1i + x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
        }
    } else {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            x0r = a[j] - a[j1];
            x0i = a[j + 1] - a[j1 + 1];
            a[j] += a[j1];
            a[j + 1] += a[j1 + 1];
            a[j1] = x0r;
            a[j1 + 1] = x0i;
        }
    }
}


static void cftbsub(int n, double *a)
{
    int j, j1, j2, j3, l;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 2;
    if (n > 8) {
        cft1st(n, a);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a);
            l <<= 2;
        }
    }
    if ((l << 2) == n) {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = -a[j + 1] - a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = -a[j + 1] + a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i - x2i;
            a[j2] = x0r - x2r;
            a[j2 + 1] = x0i + x2i;
            a[j1] = x1r - x3i;
            a[j1 + 1] = x1i - x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i + x3r;
        }
    } else {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            x0r = a[j] - a[j1];
            x0i = -a[j + 1] + a[j1 + 1];
            a[j] += a[j1];
            a[j + 1] = -a[j + 1] - a[j1 + 1];
            a[j1] = x0r;
            a[j1 + 1] = x0i;
        }
    }
}


#ifndef RDFT_LOOP_DIV
#define RDFT_LOOP_DIV 32
#endif


static void rftfsub(int n, double *a)
{
    int i, j, k;
    double ec, w1r, w1i, wkr, wki, wdr, wdi, ss, xr, xi, yr, yi;

    ec = 2 * M_PI_2 / n;
    wkr = 0;
    wki = 0;
    wdi = cos(ec);
    wdr = sin(ec);
    wdi *= wdr;
    wdr *= wdr;
    w1r = 1 - 2 * wdr;
    w1i = 2 * wdi;
    ss = 4 * wdi;
    k = n >> 1;
    while (k >= 4 * (RDFT_LOOP_DIV + 1)) {
        for (i = 0; i < RDFT_LOOP_DIV; i++) {
            k -= 4;
            j = n - k;
            xr = a[k + 2] - a[j - 2];
            xi = a[k + 3] + a[j - 1];
            yr = wdr * xr - wdi * xi;
            yi = wdr * xi + wdi * xr;
            a[k + 2] -= yr;
            a[k + 3] -= yi;
            a[j - 2] += yr;
            a[j - 1] -= yi;
            wkr += ss * wdi;
            wki += ss * (0.5 - wdr);
            xr = a[k] - a[j];
            xi = a[k + 1] + a[j + 1];
            yr = wkr * xr - wki * xi;
            yi = wkr * xi + wki * xr;
            a[k] -= yr;
            a[k + 1] -= yi;
            a[j] += yr;
            a[j + 1] -= yi;
            wdr += ss * wki;
            wdi += ss * (0.5 - wkr);
        }
        wkr = 0.5 * sin(ec * k);
        wki = 0.5 * cos(ec * k);
        wdr = 0.5 - (wkr * w1r - wki * w1i);
        wdi = wkr * w1i + wki * w1r;
        wkr = 0.5 - wkr;
    }
    while (k >= 8) {
        k -= 4;
        j = n - k;
        xr = a[k + 2] - a[j - 2];
        xi = a[k + 3] + a[j - 1];
        yr = wdr * xr - wdi * xi;
        yi = wdr * xi + wdi * xr;
        a[k + 2] -= yr;
        a[k + 3] -= yi;
        a[j - 2] += yr;
        a[j - 1] -= yi;
        wkr += ss * wdi;
        wki += ss * (0.5 - wdr);
        xr = a[k] - a[j];
        xi = a[k + 1] + a[j + 1];
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        a[k] -= yr;
        a[k + 1] -= yi;
        a[j] += yr;
        a[j + 1] -= yi;
        wdr += ss * wki;
        wdi += ss * (0.5 - wkr);
    }
    xr = a[2] - a[n - 2];
    xi = a[3] + a[n - 1];
    yr = wdr * xr - wdi * xi;
    yi = wdr * xi + wdi * xr;
    a[2] -= yr;
    a[3] -= yi;
    a[n - 2] += yr;
    a[n - 1] -= yi;
}

static void rftbsub(int n, double *a)
{
    int i, j, k;
    double ec, w1r, w1i, wkr, wki, wdr, wdi, ss, xr, xi, yr, yi;

    ec = 2 * M_PI_2 / n;
    wkr = 0;
    wki = 0;
    wdi = cos(ec);
    wdr = sin(ec);
    wdi *= wdr;
    wdr *= wdr;
    w1r = 1 - 2 * wdr;
    w1i = 2 * wdi;
    ss = 4 * wdi;
    k = n >> 1;
    a[k + 1] = -a[k + 1];
    while (k >= 4 * (RDFT_LOOP_DIV + 1)) {
        for (i = 0; i < RDFT_LOOP_DIV; i++) {
            k -= 4;
            j = n - k;
            xr = a[k + 2] - a[j - 2];
            xi = a[k + 3] + a[j - 1];
            yr = wdr * xr + wdi * xi;
            yi = wdr * xi - wdi * xr;
            a[k + 2] -= yr;
            a[k + 3] = yi - a[k + 3];
            a[j - 2] += yr;
            a[j - 1] = yi - a[j - 1];
            wkr += ss * wdi;
            wki += ss * (0.5 - wdr);
            xr = a[k] - a[j];
            xi = a[k + 1] + a[j + 1];
            yr = wkr * xr + wki * xi;
            yi = wkr * xi - wki * xr;
            a[k] -= yr;
            a[k + 1] = yi - a[k + 1];
            a[j] += yr;
            a[j + 1] = yi - a[j + 1];
            wdr += ss * wki;
            wdi += ss * (0.5 - wkr);
        }
        wkr = 0.5 * sin(ec * k);
        wki = 0.5 * cos(ec * k);
        wdr = 0.5 - (wkr * w1r - wki * w1i);
        wdi = wkr * w1i + wki * w1r;
        wkr = 0.5 - wkr;
    }
    while (k >= 8) {
        k -= 4;
        j = n - k;
        xr = a[k + 2] - a[j - 2];
        xi = a[k + 3] + a[j - 1];
        yr = wdr * xr + wdi * xi;
        yi = wdr * xi - wdi * xr;
        a[k + 2] -= yr;
        a[k + 3] = yi - a[k + 3];
        a[j - 2] += yr;
        a[j - 1] = yi - a[j - 1];
        wkr += ss * wdi;
        wki += ss * (0.5 - wdr);
        xr = a[k] - a[j];
        xi = a[k + 1] + a[j + 1];
        yr = wkr * xr + wki * xi;
        yi = wkr * xi - wki * xr;
        a[k] -= yr;
        a[k + 1] = yi - a[k + 1];
        a[j] += yr;
        a[j + 1] = yi - a[j + 1];
        wdr += ss * wki;
        wdi += ss * (0.5 - wkr);
    }
    xr = a[2] - a[n - 2];
    xi = a[3] + a[n - 1];
    yr = wdr * xr + wdi * xi;
    yi = wdr * xi - wdi * xr;
    a[2] -= yr;
    a[3] = yi - a[3];
    a[n - 2] += yr;
    a[n - 1] = yi - a[n - 1];
    a[1] = -a[1];
}

void rdft(int n, int isgn, double *a)
{
    double xi;

    if (isgn >= 0) {
        if (n > 4) {
            bitrv2(n, a);
            cftfsub(n, a);
            rftfsub(n, a);
        } else if (n == 4) {
            cftfsub(n, a);
        }
        xi = a[0] - a[1];
        a[0] += a[1];
        a[1] = xi;
    } else {
        a[1] = 0.5 * (a[0] - a[1]);
        a[0] -= a[1];
        if (n > 4) {
            rftbsub(n, a);
            bitrv2(n, a);
            cftbsub(n, a);
        } else if (n == 4) {
            cftfsub(n, a);
        }
    }
}

/*
** For loading our numbers into the FFTNum array.
*/
static void
PutNumIntoFFTNum(FFT_DATA_TYPE *FFTNum, BigInt Num, size_t NumLen)
{size_t x;

 StartTimer(LoadTime);
 if (NumLen <8) FatalError("Too small of a FFT.\n");

 {INT32 *N;
  size_t FFTLen=CalcFFTLen(NumLen);
  N=(INT32*)(FFTNum+FFTLen);
  ReadNumIntoBuf(Num,N-NumLen,NumLen);
#if FFT_DIGITS==4
  for (x=0; x < FFTLen/2; x+=2)
    {
     --N;
     FFTNum[x  ]=*N % 10000;
     FFTNum[x+1]=*N / 10000;
    }
#else
  for (x=0; x < FFTLen/2; x+=4)
    {INT32 D;
     --N;
     D=*N;
     FFTNum[x  ]= D % 100;D=D / 100;
     FFTNum[x+1]= D % 100;D=D / 100;
     FFTNum[x+2]= D % 100;D=D / 100;
     FFTNum[x+3]= D;
    }
#endif
  for (x = FFTLen/2; x < FFTLen; x++) FFTNum[x] = 0.0;
 }
  StopTimer(LoadTime);
}

void
FwdTransform(FFT_DATA_TYPE *ddata,BigInt Num,size_t NumLen)
{
PutNumIntoFFTNum(ddata, Num, NumLen);
StartTimer(FFTTime);
rdft(CalcFFTLen(NumLen),1,ddata);
StopTimer(FFTTime);
}

void
RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen)
{
StartTimer(FFTTime);
rdft(CalcFFTLen(NumLen),-1,ddata);
StopTimer(FFTTime);
}

void
InitFFT(size_t Len)
{
}

void
DeInitFFT(void)
{
}



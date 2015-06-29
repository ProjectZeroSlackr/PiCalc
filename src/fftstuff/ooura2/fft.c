#include "pi.h"
#include "fft.h"
#include "cache.h"

#ifdef LONG_DOUBLE_FFT
#error This external FFT does not support long double.  Only double.
#endif

/*
** This FFT was modified a bit to work with my pi program. CEB.
*/

#include <stdlib.h>
#include <math.h>
/*
    Copyright(C) 1999 Takuya OOURA
    You may use, copy, modify this code for any purpose and
    without fee.

Fast Fourier/Cosine/Sine Transform
    dimension   :one
    data length :power of 2
    decimation  :frequency
    radix       :4, 2
    data        :inplace
    table       :use
functions
    rdft: Real Discrete Fourier Transform
function prototypes
    void rdft(int, int, double *, int *, double *);


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
            ip[0] = 0; // first time only
            rdft(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            rdft(n, -1, a, ip, w);
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
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/2)
                        strictly,
                        length of ip >=
                            2+(1<<(int)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of
            rdft(n, 1, a, ip, w);
        is
            rdft(n, -1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


Appendix :
    The cos/sin table is recalculated when the larger table required.
    w[] and ip[] are compatible with all routines.
*/

static int *IP=NULL;
static double *W=NULL;

void bitrv2(int n, int *ip, double *a)
{
    int j, j1, k, k1, l, m, m2;
    double xr, xi;

    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 2) < l) {
        l >>= 1;
        for (j = 0; j < m; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    if ((m << 2) > l) {
        for (k = 1; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = (j << 1) + ip[k];
                k1 = (k << 1) + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    } else {
        m2 = m << 1;
        for (k = 1; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = (j << 1) + ip[k];
                k1 = (k << 1) + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += m2;
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    }
}

void rftbsub(int n, double *a, int nc, double *c)
{
    int j, k, kk, ks, m;
    double wkr, wki, xr, xi, yr, yi;

    ks = (nc << 2) / n;
    kk = 0;
    a[1] = -a[1];
    m = n >> 1;
    for (k = 2; k < m; k += 2) {
        j = n - k;
        kk += ks;
        wkr = 0.5 - c[nc - kk];
        wki = c[kk];
        xr = a[k] - a[j];
        xi = a[k + 1] + a[j + 1];
        yr = wkr * xr + wki * xi;
        yi = wkr * xi - wki * xr;
        a[k] -= yr;
        a[k + 1] = yi - a[k + 1];
        a[j] += yr;
        a[j + 1] = yi - a[j + 1];
    }
    a[m + 1] = -a[m + 1];
}

void makewt(int nw, int *ip, double *w)
{
    int nwh, j;
    double delta, x, y;

    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2) {
        nwh = nw >> 1;
        delta = atan(1.0) / nwh;
        w[0] = 1;
        w[1] = 0;
        w[nwh] = cos(delta * nwh);
        w[nwh + 1] = w[nwh];
        for (j = 2; j < nwh; j += 2) {
            x = cos(delta * j);
            y = sin(delta * j);
            w[j] = x;
            w[j + 1] = y;
            w[nw - j] = y;
            w[nw - j + 1] = x;
        }
        bitrv2(nw, ip + 2, w);
    }
}

void makect(int nc, int *ip, double *c)
{
    int nch, j;
    double delta;

    ip[1] = nc;
    if (nc > 1) {
        nch = nc >> 1;
        delta = atan(1.0) / nch;
        c[0] = cos(delta * nch);
        c[nch] = 0.5 * c[0];
        for (j = 1; j < nch; j++) {
            c[j] = 0.5 * cos(delta * j);
            c[nc - j] = 0.5 * sin(delta * j);
        }
    }
}

void cft1st(int n, double *a, double *w)
{
    int j, k1, k2;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
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
    wk1r = w[2];
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
    a[10] = wk1r * (x0r - x0i);
    a[11] = wk1r * (x0r + x0i);
    x0r = x3i + x1r;
    x0i = x3r - x1i;
    a[14] = wk1r * (x0i - x0r);
    a[15] = wk1r * (x0i + x0r);
    k1 = 0;
    for (j = 16; j < n; j += 16) {
        k1 += 2;
        k2 = k1 << 1;
        wk2r = w[k1];
        wk2i = w[k1 + 1];
        wk1r = w[k2];
        wk1i = w[k2 + 1];
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
        wk1r = w[k2 + 2];
        wk1i = w[k2 + 3];
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

void cftmdl(int n, int l, double *a, double *w)
{
    int j, j1, j2, j3, k, k1, k2, m, m2;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
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
    wk1r = w[2];
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
        a[j1] = wk1r * (x0r - x0i);
        a[j1 + 1] = wk1r * (x0r + x0i);
        x0r = x3i + x1r;
        x0i = x3r - x1i;
        a[j3] = wk1r * (x0i - x0r);
        a[j3 + 1] = wk1r * (x0i + x0r);
    }
    k1 = 0;
    m2 = m << 1;
    for (k = m2; k < n; k += m2) {
        k1 += 2;
        k2 = k1 << 1;
        wk2r = w[k1];
        wk2i = w[k1 + 1];
        wk1r = w[k2];
        wk1i = w[k2 + 1];
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
        wk1r = w[k2 + 2];
        wk1i = w[k2 + 3];
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

void cftfsub(int n, double *a, double *w)
{
    int j, j1, j2, j3, l;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 2;
    if (n > 8) {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a, w);
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

void rftfsub(int n, double *a, int nc, double *c)
{
    int j, k, kk, ks, m;
    double wkr, wki, xr, xi, yr, yi;

    ks = (nc << 2) / n;
    kk = 0;
    m = n >> 1;
    for (k = 2; k < m; k += 2) {
        j = n - k;
        kk += ks;
        wkr = 0.5 - c[nc - kk];
        wki = c[kk];
        xr = a[k] - a[j];
        xi = a[k + 1] + a[j + 1];
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        a[k] -= yr;
        a[k + 1] -= yi;
        a[j] += yr;
        a[j + 1] -= yi;
    }
}

void cftbsub(int n, double *a, double *w)
{
    int j, j1, j2, j3, l;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 2;
    if (n > 8) {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a, w);
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

void rdft(int n, int isgn, double *a, int *ip, double *w)
{
    int nw, nc;
    double xi;

    nw = ip[0];
    if (n > (nw << 2)) {
        nw = n >> 2;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > (nc << 2)) {
        nc = n >> 2;
        makect(nc, ip, w + nw);
    }
    if (isgn >= 0) {
        if (n > 4) {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
            rftfsub(n, a, nc, w + nw);
        } else if (n == 4) {
            cftfsub(n, a, w);
        }
        xi = a[0] - a[1];
        a[0] += a[1];
        a[1] = xi;
    } else {
        a[1] = 0.5 * (a[0] - a[1]);
        a[0] -= a[1];
        if (n > 4) {
            rftbsub(n, a, nc, w + nw);
            bitrv2(n, ip + 2, a);
            cftbsub(n, a, w);
        } else if (n == 4) {
            cftfsub(n, a, w);
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
rdft(CalcFFTLen(NumLen),1,ddata,IP,W);
StopTimer(FFTTime);
}

void
RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen)
{
StartTimer(FFTTime);
rdft(CalcFFTLen(NumLen),-1,ddata,IP,W);
StopTimer(FFTTime);
}

void
InitFFT(size_t Len)
{
if (W) free(W);
if (IP) free(IP);
Len=CalcFFTLen(Len);
IP = (int *) malloc((3 + (int) sqrt(0.5 * Len)) * sizeof(int));
W  = (double *) malloc(Len / 2 * sizeof(double));

IP[0]=0;
makewt(Len/4,IP,W);
makect(Len/4,IP,W+(Len/4));
}

void
DeInitFFT(void)
{
if (W) free(W);
if (IP) free(IP);
}



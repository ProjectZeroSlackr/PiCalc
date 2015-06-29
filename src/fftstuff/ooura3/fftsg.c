/*
Fast Fourier/Cosine/Sine Transform
    dimension   :one
    data length :power of 2
    decimation  :frequency
    radix       :split-radix
    data        :inplace
    table       :use
functions
    cdft: Complex Discrete Fourier Transform
    rdft: Real Discrete Fourier Transform
    ddct: Discrete Cosine Transform
    ddst: Discrete Sine Transform
    dfct: Cosine Transform of RDFT (Real Symmetric DFT)
    dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
function prototypes
    void cdft(int, int, double *, int *, double *);
    void rdft(int, int, double *, int *, double *);
    void ddct(int, int, double *, int *, double *);
    void ddst(int, int, double *, int *, double *);
    void dfct(int, double *, double *, int *, double *);
    void dfst(int, double *, double *, int *, double *);


-------- Complex DFT (Discrete Fourier Transform) --------
    [definition]
        <case1>
            X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n
        <case2>
            X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n
        (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
    [usage]
        <case1>
            ip[0] = 0; // first time only
            cdft(2*n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            cdft(2*n, -1, a, ip, w);
    [parameters]
        2*n            :data length (int)
                        n >= 1, n = power of 2
        a[0...2*n-1]   :input/output data (double *)
                        input data
                            a[2*j] = Re(x[j]), 
                            a[2*j+1] = Im(x[j]), 0<=j<n
                        output data
                            a[2*k] = Re(X[k]), 
                            a[2*k+1] = Im(X[k]), 0<=k<n
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            cdft(2*n, -1, a, ip, w);
        is 
            cdft(2*n, 1, a, ip, w);
            for (j = 0; j <= 2 * n - 1; j++) {
                a[j] *= 1.0 / n;
            }
        .


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


-------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
    [definition]
        <case1> IDCT (excluding scale)
            C[k] = sum_j=0^n-1 a[j]*cos(pi*j*(k+1/2)/n), 0<=k<n
        <case2> DCT
            C[k] = sum_j=0^n-1 a[j]*cos(pi*(j+1/2)*k/n), 0<=k<n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            ddct(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            ddct(n, -1, a, ip, w);
    [parameters]
        n              :data length (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (double *)
                        output data
                            a[k] = C[k], 0<=k<n
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/2)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/4-1] :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            ddct(n, -1, a, ip, w);
        is 
            a[0] *= 0.5;
            ddct(n, 1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- DST (Discrete Sine Transform) / Inverse of DST --------
    [definition]
        <case1> IDST (excluding scale)
            S[k] = sum_j=1^n A[j]*sin(pi*j*(k+1/2)/n), 0<=k<n
        <case2> DST
            S[k] = sum_j=0^n-1 a[j]*sin(pi*(j+1/2)*k/n), 0<k<=n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            ddst(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            ddst(n, -1, a, ip, w);
    [parameters]
        n              :data length (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (double *)
                        <case1>
                            input data
                                a[j] = A[j], 0<j<n
                                a[0] = A[n]
                            output data
                                a[k] = S[k], 0<=k<n
                        <case2>
                            output data
                                a[k] = S[k], 0<k<n
                                a[0] = S[n]
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/2)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/4-1] :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            ddst(n, -1, a, ip, w);
        is 
            a[0] *= 0.5;
            ddst(n, 1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- Cosine Transform of RDFT (Real Symmetric DFT) --------
    [definition]
        C[k] = sum_j=0^n a[j]*cos(pi*j*k/n), 0<=k<=n
    [usage]
        ip[0] = 0; // first time only
        dfct(n, a, t, ip, w);
    [parameters]
        n              :data length - 1 (int)
                        n >= 2, n = power of 2
        a[0...n]       :input/output data (double *)
                        output data
                            a[k] = C[k], 0<=k<=n
        t[0...n/2]     :work area (double *)
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/4)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/4+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/8-1] :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            a[0] *= 0.5;
            a[n] *= 0.5;
            dfct(n, a, t, ip, w);
        is 
            a[0] *= 0.5;
            a[n] *= 0.5;
            dfct(n, a, t, ip, w);
            for (j = 0; j <= n; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
    [definition]
        S[k] = sum_j=1^n-1 a[j]*sin(pi*j*k/n), 0<k<n
    [usage]
        ip[0] = 0; // first time only
        dfst(n, a, t, ip, w);
    [parameters]
        n              :data length + 1 (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (double *)
                        output data
                            a[k] = S[k], 0<k<n
                        (a[0] is used for work area)
        t[0...n/2-1]   :work area (double *)
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/4)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/4+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/8-1] :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            dfst(n, a, t, ip, w);
        is 
            dfst(n, a, t, ip, w);
            for (j = 1; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


Appendix :
    The cos/sin table is recalculated when the larger table required.
    w[] and ip[] are compatible with all routines.
*/


void cdft(int n, int isgn, double *a, int *ip, double *w)
{
    void makewt(int nw, int *ip, double *w);
    void bitrv2(int n, int *ip, double *a);
    void cftfsub(int n, double *a, int nw, double *w);
    void cftbsub(int n, double *a, int nw, double *w);
    int nw;

    nw = ip[0];
    if (n > (nw << 2)) {
        nw = n >> 2;
        makewt(nw, ip, w);
    }
    if (isgn >=0) {
        cftfsub(n, a, nw, w);
    } else {
        cftbsub(n, a, nw, w);
    }
    if (n > 4) {
        bitrv2(n, ip + 2, a);
    }
}


void rdft(int n, int isgn, double *a, int *ip, double *w)
{
    void makewt(int nw, int *ip, double *w);
    void makect(int nc, int *ip, double *c);
    void bitrv2(int n, int *ip, double *a);
    void cftfsub(int n, double *a, int nw, double *w);
    void cftbsub(int n, double *a, int nw, double *w);
    void rftfsub(int n, double *a, int nc, double *c);
    void rftbsub(int n, double *a, int nc, double *c);
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
            cftfsub(n, a, nw, w);
            bitrv2(n, ip + 2, a);
            rftfsub(n, a, nc, w + nw);
        } else if (n == 4) {
            cftfsub(n, a, nw, w);
        }
        xi = a[0] - a[1];
        a[0] += a[1];
        a[1] = xi;
    } else {
        a[1] = 0.5 * (a[0] - a[1]);
        a[0] -= a[1];
        if (n > 4) {
            rftbsub(n, a, nc, w + nw);
            cftbsub(n, a, nw, w);
            bitrv2(n, ip + 2, a);
        } else if (n == 4) {
            cftfsub(n, a, nw, w);
        }
    }
}


void ddct(int n, int isgn, double *a, int *ip, double *w)
{
    void makewt(int nw, int *ip, double *w);
    void makect(int nc, int *ip, double *c);
    void bitrv2(int n, int *ip, double *a);
    void cftfsub(int n, double *a, int nw, double *w);
    void cftbsub(int n, double *a, int nw, double *w);
    void rftfsub(int n, double *a, int nc, double *c);
    void rftbsub(int n, double *a, int nc, double *c);
    void dctsub(int n, double *a, int nc, double *c);
    int j, nw, nc;
    double xr;
    
    nw = ip[0];
    if (n > (nw << 2)) {
        nw = n >> 2;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > nc) {
        nc = n;
        makect(nc, ip, w + nw);
    }
    if (isgn < 0) {
        xr = a[n - 1];
        for (j = n - 2; j >= 2; j -= 2) {
            a[j + 1] = a[j] - a[j - 1];
            a[j] += a[j - 1];
        }
        a[1] = a[0] - xr;
        a[0] += xr;
        if (n > 4) {
            rftbsub(n, a, nc, w + nw);
            cftbsub(n, a, nw, w);
            bitrv2(n, ip + 2, a);
        } else if (n == 4) {
            cftfsub(n, a, nw, w);
        }
    }
    dctsub(n, a, nc, w + nw);
    if (isgn >= 0) {
        if (n > 4) {
            cftfsub(n, a, nw, w);
            bitrv2(n, ip + 2, a);
            rftfsub(n, a, nc, w + nw);
        } else if (n == 4) {
            cftfsub(n, a, nw, w);
        }
        xr = a[0] - a[1];
        a[0] += a[1];
        for (j = 2; j < n; j += 2) {
            a[j - 1] = a[j] - a[j + 1];
            a[j] += a[j + 1];
        }
        a[n - 1] = xr;
    }
}


void ddst(int n, int isgn, double *a, int *ip, double *w)
{
    void makewt(int nw, int *ip, double *w);
    void makect(int nc, int *ip, double *c);
    void bitrv2(int n, int *ip, double *a);
    void cftfsub(int n, double *a, int nw, double *w);
    void cftbsub(int n, double *a, int nw, double *w);
    void rftfsub(int n, double *a, int nc, double *c);
    void rftbsub(int n, double *a, int nc, double *c);
    void dstsub(int n, double *a, int nc, double *c);
    int j, nw, nc;
    double xr;
    
    nw = ip[0];
    if (n > (nw << 2)) {
        nw = n >> 2;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > nc) {
        nc = n;
        makect(nc, ip, w + nw);
    }
    if (isgn < 0) {
        xr = a[n - 1];
        for (j = n - 2; j >= 2; j -= 2) {
            a[j + 1] = -a[j] - a[j - 1];
            a[j] -= a[j - 1];
        }
        a[1] = a[0] + xr;
        a[0] -= xr;
        if (n > 4) {
            rftbsub(n, a, nc, w + nw);
            cftbsub(n, a, nw, w);
            bitrv2(n, ip + 2, a);
        } else if (n == 4) {
            cftfsub(n, a, nw, w);
        }
    }
    dstsub(n, a, nc, w + nw);
    if (isgn >= 0) {
        if (n > 4) {
            cftfsub(n, a, nw, w);
            bitrv2(n, ip + 2, a);
            rftfsub(n, a, nc, w + nw);
        } else if (n == 4) {
            cftfsub(n, a, nw, w);
        }
        xr = a[0] - a[1];
        a[0] += a[1];
        for (j = 2; j < n; j += 2) {
            a[j - 1] = -a[j] - a[j + 1];
            a[j] -= a[j + 1];
        }
        a[n - 1] = -xr;
    }
}


void dfct(int n, double *a, double *t, int *ip, double *w)
{
    void makewt(int nw, int *ip, double *w);
    void makect(int nc, int *ip, double *c);
    void bitrv2(int n, int *ip, double *a);
    void cftfsub(int n, double *a, int nw, double *w);
    void rftfsub(int n, double *a, int nc, double *c);
    void dctsub(int n, double *a, int nc, double *c);
    int j, k, l, m, mh, nw, nc;
    double xr, xi;
    
    nw = ip[0];
    if (n > (nw << 3)) {
        nw = n >> 3;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > (nc << 1)) {
        nc = n >> 1;
        makect(nc, ip, w + nw);
    }
    m = n >> 1;
    xr = a[0] + a[n];
    a[0] -= a[n];
    t[0] = xr - a[m];
    t[m] = xr + a[m];
    if (n > 2) {
        mh = m >> 1;
        for (j = 1; j < mh; j++) {
            k = m - j;
            xr = a[j] + a[n - j];
            a[j] -= a[n - j];
            xi = a[k] + a[n - k];
            a[k] -= a[n - k];
            t[j] = xr - xi;
            t[k] = xr + xi;
        }
        t[mh] = a[mh] + a[n - mh];
        a[mh] -= a[n - mh];
        dctsub(m, a, nc, w + nw);
        if (m > 4) {
            cftfsub(m, a, nw, w);
            bitrv2(m, ip + 2, a);
            rftfsub(m, a, nc, w + nw);
        } else if (m == 4) {
            cftfsub(m, a, nw, w);
        }
        xr = a[0] + a[1];
        a[n - 1] = a[0] - a[1];
        for (j = m - 2; j >= 2; j -= 2) {
            a[(j << 1) + 1] = a[j] + a[j + 1];
            a[(j << 1) - 1] = a[j] - a[j + 1];
        }
        a[1] = xr;
        l = 2;
        m = mh;
        while (m >= 2) {
            dctsub(m, t, nc, w + nw);
            if (m > 4) {
                cftfsub(m, t, nw, w);
                bitrv2(m, ip + 2, t);
                rftfsub(m, t, nc, w + nw);
            } else if (m == 4) {
                cftfsub(m, t, nw, w);
            }
            a[n - l] = t[0] - t[1];
            a[l] = t[0] + t[1];
            k = 0;
            for (j = 2; j < m; j += 2) {
                k += l << 2;
                a[k - l] = t[j] - t[j + 1];
                a[k + l] = t[j] + t[j + 1];
            }
            l <<= 1;
            mh = m >> 1;
            for (j = 0; j < mh; j++) {
                k = m - j;
                t[j] = t[m + k] - t[m + j];
                t[k] = t[m + k] + t[m + j];
            }
            t[mh] = t[m + mh];
            m = mh;
        }
        a[l] = t[0];
        a[n] = t[2] - t[1];
        a[0] = t[2] + t[1];
    } else {
        a[1] = a[0];
        a[2] = t[0];
        a[0] = t[1];
    }
}


void dfst(int n, double *a, double *t, int *ip, double *w)
{
    void makewt(int nw, int *ip, double *w);
    void makect(int nc, int *ip, double *c);
    void bitrv2(int n, int *ip, double *a);
    void cftfsub(int n, double *a, int nw, double *w);
    void rftfsub(int n, double *a, int nc, double *c);
    void dstsub(int n, double *a, int nc, double *c);
    int j, k, l, m, mh, nw, nc;
    double xr, xi;
    
    nw = ip[0];
    if (n > (nw << 3)) {
        nw = n >> 3;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > (nc << 1)) {
        nc = n >> 1;
        makect(nc, ip, w + nw);
    }
    if (n > 2) {
        m = n >> 1;
        mh = m >> 1;
        for (j = 1; j < mh; j++) {
            k = m - j;
            xr = a[j] - a[n - j];
            a[j] += a[n - j];
            xi = a[k] - a[n - k];
            a[k] += a[n - k];
            t[j] = xr + xi;
            t[k] = xr - xi;
        }
        t[0] = a[mh] - a[n - mh];
        a[mh] += a[n - mh];
        a[0] = a[m];
        dstsub(m, a, nc, w + nw);
        if (m > 4) {
            cftfsub(m, a, nw, w);
            bitrv2(m, ip + 2, a);
            rftfsub(m, a, nc, w + nw);
        } else if (m == 4) {
            cftfsub(m, a, nw, w);
        }
        xr = a[0] + a[1];
        a[n - 1] = a[1] - a[0];
        for (j = m - 2; j >= 2; j -= 2) {
            a[(j << 1) + 1] = a[j] - a[j + 1];
            a[(j << 1) - 1] = -a[j] - a[j + 1];
        }
        a[1] = xr;
        l = 2;
        m = mh;
        while (m >= 2) {
            dstsub(m, t, nc, w + nw);
            if (m > 4) {
                cftfsub(m, t, nw, w);
                bitrv2(m, ip + 2, t);
                rftfsub(m, t, nc, w + nw);
            } else if (m == 4) {
                cftfsub(m, t, nw, w);
            }
            a[n - l] = t[1] - t[0];
            a[l] = t[0] + t[1];
            k = 0;
            for (j = 2; j < m; j += 2) {
                k += l << 2;
                a[k - l] = -t[j] - t[j + 1];
                a[k + l] = t[j] - t[j + 1];
            }
            l <<= 1;
            mh = m >> 1;
            for (j = 1; j < mh; j++) {
                k = m - j;
                t[j] = t[m + k] + t[m + j];
                t[k] = t[m + k] - t[m + j];
            }
            t[0] = t[m + mh];
            m = mh;
        }
        a[l] = t[0];
    }
    a[0] = 0;
}


/* -------- initializing routines -------- */


#include <math.h>


void makewt(int nw, int *ip, double *w)
{
    int nwh, nw0, nw1, j;
    double delta;
    
    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2) {
        nwh = nw >> 1;
        delta = atan(1.0) / nwh;
        w[0] = 1;
        w[1] = cos(delta * nwh);
        if (nwh > 2) {
            w[2] = 0.5 / cos(delta * 2);
            w[3] = 0.5 / cos(delta * 6);
        }
        for (j = 4; j < nwh; j += 4) {
            w[j] = cos(delta * j);
            w[j + 1] = sin(delta * j);
            w[j + 2] = cos(3 * delta * j);
            w[j + 3] = sin(3 * delta * j);
        }
        nw0 = 0;
        while (nwh > 2) {
            nw1 = nw0 + nwh;
            nwh >>= 1;
            w[nw1] = 1;
            w[nw1 + 1] = w[1];
            if (nwh > 2) {
                w[nw1 + 2] = 0.5 / w[nw0 + 4];
                w[nw1 + 3] = 0.5 / w[nw0 + 6];
            }
            for (j = 4; j < nwh; j += 4) {
                w[nw1 + j] = w[nw0 + 2 * j];
                w[nw1 + j + 1] = w[nw0 + 2 * j + 1];
                w[nw1 + j + 2] = w[nw0 + 2 * j + 2];
                w[nw1 + j + 3] = w[nw0 + 2 * j + 3];
            }
            nw0 = nw1;
        }
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


/* -------- child routines -------- */


#ifndef CDFT_RECURSIVE_N
#define CDFT_RECURSIVE_N 512
#endif


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


void cftfsub(int n, double *a, int nw, double *w)
{
    void cftfsub(int n, double *a, int nw, double *w);
    void cftf1st(int n, double *a, double *w);
    void cftfmdl(int n, double *a, double *w);
    void cftffix(int n, double wn4r, double *a);
    int nq, nh, m, j, k;
    double x0r, x0i;
    
    if (n <= CDFT_RECURSIVE_N) {
        if (n > 16) {
            if (n > 2 * nw) {
                cftf1st(n, a, w);
            } else {
                cftfmdl(n, a, &w[nw - (n >> 1)]);
            }
            if (n > 32) {
                cftfmdl(n >> 1, a, &w[nw - (n >> 2)]);
            }
            for (m = n >> 2; m > 16; m >>= 1) {
                for (k = m; k <= n; k <<= 2) {
                    for (j = k - m; j < n; j += 2 * k) {
                        cftfmdl(m, &a[j], &w[nw - (m >> 1)]);
                    }
                }
            }
        }
        if (n > 4) {
            cftffix(n, w[1], a);
        } else if (n == 4) {
            x0r = a[0] - a[2];
            x0i = a[1] - a[3];
            a[0] += a[2];
            a[1] += a[3];
            a[2] = x0r;
            a[3] = x0i;
        }
    } else {
        nq = n >> 2;
        nh = nq + nq;
        if (n > 2 * nw) {
            cftf1st(n, a, w);
        } else {
            cftfmdl(n, a, &w[nw - (n >> 1)]);
        }
        cftfsub(nh, a, nw, w);
        cftfsub(nq, &a[nh], nw, w);
        cftfsub(nq, &a[nh + nq], nw, w);
    }
}


void cftbsub(int n, double *a, int nw, double *w)
{
    void cftbsub(int n, double *a, int nw, double *w);
    void cftb1st(int n, double *a, double *w);
    void cftbmdl(int n, double *a, double *w);
    void cftbfix(int n, double wn4r, double *a);
    int nq, nh, m, j, k;
    double x0r, x0i;
    
    if (n <= CDFT_RECURSIVE_N) {
        if (n > 16) {
            if (n > 2 * nw) {
                cftb1st(n, a, w);
            } else {
                cftbmdl(n, a, &w[nw - (n >> 1)]);
            }
            if (n > 32) {
                cftbmdl(n >> 1, a, &w[nw - (n >> 2)]);
            }
            for (m = n >> 2; m > 16; m >>= 1) {
                for (k = m; k <= n; k <<= 2) {
                    for (j = k - m; j < n; j += 2 * k) {
                        cftbmdl(m, &a[j], &w[nw - (m >> 1)]);
                    }
                }
            }
        }
        if (n > 4) {
            cftbfix(n, w[1], a);
        } else if (n == 4) {
            x0r = a[0] - a[2];
            x0i = a[1] - a[3];
            a[0] += a[2];
            a[1] += a[3];
            a[2] = x0r;
            a[3] = x0i;
        }
    } else {
        nq = n >> 2;
        nh = nq + nq;
        if (n > 2 * nw) {
            cftb1st(n, a, w);
        } else {
            cftbmdl(n, a, &w[nw - (n >> 1)]);
        }
        cftbsub(nh, a, nw, w);
        cftbsub(nq, &a[nh], nw, w);
        cftbsub(nq, &a[nh + nq], nw, w);
    }
}


void cftf1st(int n, double *a, double *w)
{
    int nqh, nq, j, j0, j1, j2, j3, k;
    double x0r, x0i, x1r, x1i, x3r, x3i;
    double wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i, 
        wd1r, wd1i, wd3r, wd3i;
    
    nqh = n >> 3;
    nq = nqh + nqh;
    j1 = nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[0] - a[j2];
    x1i = a[1] - a[j2 + 1];
    a[0] += a[j2];
    a[1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    a[j2] = x1r - x3i;
    a[j2 + 1] = x1i + x3r;
    a[j3] = x1r + x3i;
    a[j3 + 1] = x1i - x3r;
    wn4r = w[1];
    csc1 = w[2];
    csc3 = w[3];
    wk1r = 1;
    wk1i = 0;
    wk3r = 1;
    wk3i = 0;
    k = 0;
    for (j = 2; j < nqh - 2; j += 4) {
        k += 4;
        wd1r = w[k];
        wd1i = w[k + 1];
        wd3r = w[k + 2];
        wd3i = w[k + 3];
        wk1r = csc1 * (wk1r + wd1r);
        wk1i = csc1 * (wk1i + wd1i);
        wk3r = csc3 * (wk3r + wd3r);
        wk3i = csc3 * (wk3i + wd3i);
        j1 = j + nq;
        j2 = j1 + nq;
        j3 = j2 + nq;
        x1r = a[j] - a[j2];
        x1i = a[j + 1] - a[j2 + 1];
        a[j] += a[j2];
        a[j + 1] += a[j2 + 1];
        x3r = a[j1] - a[j3];
        x3i = a[j1 + 1] - a[j3 + 1];
        a[j1] += a[j3];
        a[j1 + 1] += a[j3 + 1];
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j2] = wk1r * x0r - wk1i * x0i;
        a[j2 + 1] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j3] = wk3r * x0r - wk3i * x0i;
        a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        x1r = a[j + 2] - a[j2 + 2];
        x1i = a[j + 3] - a[j2 + 3];
        a[j + 2] += a[j2 + 2];
        a[j + 3] += a[j2 + 3];
        x3r = a[j1 + 2] - a[j3 + 2];
        x3i = a[j1 + 3] - a[j3 + 3];
        a[j1 + 2] += a[j3 + 2];
        a[j1 + 3] += a[j3 + 3];
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j2 + 2] = wd1r * x0r - wd1i * x0i;
        a[j2 + 3] = wd1r * x0i + wd1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j3 + 2] = wd3r * x0r - wd3i * x0i;
        a[j3 + 3] = wd3r * x0i + wd3i * x0r;
        j0 = nq - j;
        j1 = j0 + nq;
        j2 = j1 + nq;
        j3 = j2 + nq;
        x1r = a[j0] - a[j2];
        x1i = a[j0 + 1] - a[j2 + 1];
        a[j0] += a[j2];
        a[j0 + 1] += a[j2 + 1];
        x3r = a[j1] - a[j3];
        x3i = a[j1 + 1] - a[j3 + 1];
        a[j1] += a[j3];
        a[j1 + 1] += a[j3 + 1];
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j2] = wk1i * x0r - wk1r * x0i;
        a[j2 + 1] = wk1i * x0i + wk1r * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j3] = -wk3i * x0r + wk3r * x0i;
        a[j3 + 1] = -wk3i * x0i - wk3r * x0r;
        x1r = a[j0 - 2] - a[j2 - 2];
        x1i = a[j0 - 1] - a[j2 - 1];
        a[j0 - 2] += a[j2 - 2];
        a[j0 - 1] += a[j2 - 1];
        x3r = a[j1 - 2] - a[j3 - 2];
        x3i = a[j1 - 1] - a[j3 - 1];
        a[j1 - 2] += a[j3 - 2];
        a[j1 - 1] += a[j3 - 1];
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j2 - 2] = wd1i * x0r - wd1r * x0i;
        a[j2 - 1] = wd1i * x0i + wd1r * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j3 - 2] = -wd3i * x0r + wd3r * x0i;
        a[j3 - 1] = -wd3i * x0i - wd3r * x0r;
        wk1r = wd1r;
        wk1i = wd1i;
        wk3r = wd3r;
        wk3i = wd3i;
    }
    wk1r = csc1 * (wk1r + wn4r);
    wk1i = csc1 * (wk1i + wn4r);
    wk3r = csc3 * (wk3r - wn4r);
    wk3i = csc3 * (wk3i + wn4r);
    j0 = nqh - 2;
    j1 = j0 + nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[j0] - a[j2];
    x1i = a[j0 + 1] - a[j2 + 1];
    a[j0] += a[j2];
    a[j0 + 1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j2] = wk1r * x0r - wk1i * x0i;
    a[j2 + 1] = wk1r * x0i + wk1i * x0r;
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j3] = wk3r * x0r - wk3i * x0i;
    a[j3 + 1] = wk3r * x0i + wk3i * x0r;
    j0 = nqh + 2;
    j1 = j0 + nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[j0] - a[j2];
    x1i = a[j0 + 1] - a[j2 + 1];
    a[j0] += a[j2];
    a[j0 + 1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j2] = wk1i * x0r - wk1r * x0i;
    a[j2 + 1] = wk1i * x0i + wk1r * x0r;
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j3] = -wk3i * x0r + wk3r * x0i;
    a[j3 + 1] = -wk3i * x0i - wk3r * x0r;
    j0 = nqh;
    j1 = j0 + nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[j0] - a[j2];
    x1i = a[j0 + 1] - a[j2 + 1];
    a[j0] += a[j2];
    a[j0 + 1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j2] = wn4r * (x0r - x0i);
    a[j2 + 1] = wn4r * (x0i + x0r);
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j3] = -wn4r * (x0r + x0i);
    a[j3 + 1] = -wn4r * (x0i - x0r);
}


void cftb1st(int n, double *a, double *w)
{
    int nqh, nq, j, j0, j1, j2, j3, k;
    double x0r, x0i, x1r, x1i, x3r, x3i;
    double wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i, 
        wd1r, wd1i, wd3r, wd3i;
    
    nqh = n >> 3;
    nq = nqh + nqh;
    j1 = nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[0] - a[j2];
    x1i = a[1] - a[j2 + 1];
    a[0] += a[j2];
    a[1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    a[j2] = x1r + x3i;
    a[j2 + 1] = x1i - x3r;
    a[j3] = x1r - x3i;
    a[j3 + 1] = x1i + x3r;
    wn4r = w[1];
    csc1 = w[2];
    csc3 = w[3];
    wk1r = 1;
    wk1i = 0;
    wk3r = 1;
    wk3i = 0;
    k = 0;
    for (j = 2; j < nqh - 2; j += 4) {
        k += 4;
        wd1r = w[k];
        wd1i = w[k + 1];
        wd3r = w[k + 2];
        wd3i = w[k + 3];
        wk1r = csc1 * (wk1r + wd1r);
        wk1i = csc1 * (wk1i + wd1i);
        wk3r = csc3 * (wk3r + wd3r);
        wk3i = csc3 * (wk3i + wd3i);
        j1 = j + nq;
        j2 = j1 + nq;
        j3 = j2 + nq;
        x1r = a[j] - a[j2];
        x1i = a[j + 1] - a[j2 + 1];
        a[j] += a[j2];
        a[j + 1] += a[j2 + 1];
        x3r = a[j1] - a[j3];
        x3i = a[j1 + 1] - a[j3 + 1];
        a[j1] += a[j3];
        a[j1 + 1] += a[j3 + 1];
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j2] = wk1r * x0r + wk1i * x0i;
        a[j2 + 1] = wk1r * x0i - wk1i * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j3] = wk3r * x0r + wk3i * x0i;
        a[j3 + 1] = wk3r * x0i - wk3i * x0r;
        x1r = a[j + 2] - a[j2 + 2];
        x1i = a[j + 3] - a[j2 + 3];
        a[j + 2] += a[j2 + 2];
        a[j + 3] += a[j2 + 3];
        x3r = a[j1 + 2] - a[j3 + 2];
        x3i = a[j1 + 3] - a[j3 + 3];
        a[j1 + 2] += a[j3 + 2];
        a[j1 + 3] += a[j3 + 3];
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j2 + 2] = wd1r * x0r + wd1i * x0i;
        a[j2 + 3] = wd1r * x0i - wd1i * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j3 + 2] = wd3r * x0r + wd3i * x0i;
        a[j3 + 3] = wd3r * x0i - wd3i * x0r;
        j0 = nq - j;
        j1 = j0 + nq;
        j2 = j1 + nq;
        j3 = j2 + nq;
        x1r = a[j0] - a[j2];
        x1i = a[j0 + 1] - a[j2 + 1];
        a[j0] += a[j2];
        a[j0 + 1] += a[j2 + 1];
        x3r = a[j1] - a[j3];
        x3i = a[j1 + 1] - a[j3 + 1];
        a[j1] += a[j3];
        a[j1 + 1] += a[j3 + 1];
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j2] = wk1i * x0r + wk1r * x0i;
        a[j2 + 1] = wk1i * x0i - wk1r * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j3] = -wk3i * x0r - wk3r * x0i;
        a[j3 + 1] = -wk3i * x0i + wk3r * x0r;
        x1r = a[j0 - 2] - a[j2 - 2];
        x1i = a[j0 - 1] - a[j2 - 1];
        a[j0 - 2] += a[j2 - 2];
        a[j0 - 1] += a[j2 - 1];
        x3r = a[j1 - 2] - a[j3 - 2];
        x3i = a[j1 - 1] - a[j3 - 1];
        a[j1 - 2] += a[j3 - 2];
        a[j1 - 1] += a[j3 - 1];
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j2 - 2] = wd1i * x0r + wd1r * x0i;
        a[j2 - 1] = wd1i * x0i - wd1r * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j3 - 2] = -wd3i * x0r - wd3r * x0i;
        a[j3 - 1] = -wd3i * x0i + wd3r * x0r;
        wk1r = wd1r;
        wk1i = wd1i;
        wk3r = wd3r;
        wk3i = wd3i;
    }
    wk1r = csc1 * (wk1r + wn4r);
    wk1i = csc1 * (wk1i + wn4r);
    wk3r = csc3 * (wk3r - wn4r);
    wk3i = csc3 * (wk3i + wn4r);
    j0 = nqh - 2;
    j1 = j0 + nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[j0] - a[j2];
    x1i = a[j0 + 1] - a[j2 + 1];
    a[j0] += a[j2];
    a[j0 + 1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j2] = wk1r * x0r + wk1i * x0i;
    a[j2 + 1] = wk1r * x0i - wk1i * x0r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j3] = wk3r * x0r + wk3i * x0i;
    a[j3 + 1] = wk3r * x0i - wk3i * x0r;
    j0 = nqh + 2;
    j1 = j0 + nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[j0] - a[j2];
    x1i = a[j0 + 1] - a[j2 + 1];
    a[j0] += a[j2];
    a[j0 + 1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j2] = wk1i * x0r + wk1r * x0i;
    a[j2 + 1] = wk1i * x0i - wk1r * x0r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j3] = -wk3i * x0r - wk3r * x0i;
    a[j3 + 1] = -wk3i * x0i + wk3r * x0r;
    j0 = nqh;
    j1 = j0 + nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[j0] - a[j2];
    x1i = a[j0 + 1] - a[j2 + 1];
    a[j0] += a[j2];
    a[j0 + 1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j2] = wn4r * (x0r + x0i);
    a[j2 + 1] = wn4r * (x0i - x0r);
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j3] = -wn4r * (x0r - x0i);
    a[j3 + 1] = -wn4r * (x0i + x0r);
}


void cftfmdl(int n, double *a, double *w)
{
    int nqh, nq, j, j0, j1, j2, j3, k;
    double x0r, x0i, x1r, x1i, x3r, x3i;
    double wn4r, wk1r, wk1i, wk3r, wk3i;
    
    nqh = n >> 3;
    nq = nqh + nqh;
    j1 = nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[0] - a[j2];
    x1i = a[1] - a[j2 + 1];
    a[0] += a[j2];
    a[1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    a[j2] = x1r - x3i;
    a[j2 + 1] = x1i + x3r;
    a[j3] = x1r + x3i;
    a[j3 + 1] = x1i - x3r;
    wn4r = w[1];
    k = 0;
    for (j = 2; j < nqh; j += 2) {
        k += 4;
        wk1r = w[k];
        wk1i = w[k + 1];
        wk3r = w[k + 2];
        wk3i = w[k + 3];
        j1 = j + nq;
        j2 = j1 + nq;
        j3 = j2 + nq;
        x1r = a[j] - a[j2];
        x1i = a[j + 1] - a[j2 + 1];
        a[j] += a[j2];
        a[j + 1] += a[j2 + 1];
        x3r = a[j1] - a[j3];
        x3i = a[j1 + 1] - a[j3 + 1];
        a[j1] += a[j3];
        a[j1 + 1] += a[j3 + 1];
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j2] = wk1r * x0r - wk1i * x0i;
        a[j2 + 1] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j3] = wk3r * x0r - wk3i * x0i;
        a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        j0 = nq - j;
        j1 = j0 + nq;
        j2 = j1 + nq;
        j3 = j2 + nq;
        x1r = a[j0] - a[j2];
        x1i = a[j0 + 1] - a[j2 + 1];
        a[j0] += a[j2];
        a[j0 + 1] += a[j2 + 1];
        x3r = a[j1] - a[j3];
        x3i = a[j1 + 1] - a[j3 + 1];
        a[j1] += a[j3];
        a[j1 + 1] += a[j3 + 1];
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j2] = wk1i * x0r - wk1r * x0i;
        a[j2 + 1] = wk1i * x0i + wk1r * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j3] = -wk3i * x0r + wk3r * x0i;
        a[j3 + 1] = -wk3i * x0i - wk3r * x0r;
    }
    j0 = nqh;
    j1 = j0 + nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[j0] - a[j2];
    x1i = a[j0 + 1] - a[j2 + 1];
    a[j0] += a[j2];
    a[j0 + 1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j2] = wn4r * (x0r - x0i);
    a[j2 + 1] = wn4r * (x0i + x0r);
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j3] = -wn4r * (x0r + x0i);
    a[j3 + 1] = -wn4r * (x0i - x0r);
}


void cftbmdl(int n, double *a, double *w)
{
    int nqh, nq, j, j0, j1, j2, j3, k;
    double x0r, x0i, x1r, x1i, x3r, x3i;
    double wn4r, wk1r, wk1i, wk3r, wk3i;
    
    nqh = n >> 3;
    nq = nqh + nqh;
    j1 = nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[0] - a[j2];
    x1i = a[1] - a[j2 + 1];
    a[0] += a[j2];
    a[1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    a[j2] = x1r + x3i;
    a[j2 + 1] = x1i - x3r;
    a[j3] = x1r - x3i;
    a[j3 + 1] = x1i + x3r;
    wn4r = w[1];
    k = 0;
    for (j = 2; j < nqh; j += 2) {
        k += 4;
        wk1r = w[k];
        wk1i = w[k + 1];
        wk3r = w[k + 2];
        wk3i = w[k + 3];
        j1 = j + nq;
        j2 = j1 + nq;
        j3 = j2 + nq;
        x1r = a[j] - a[j2];
        x1i = a[j + 1] - a[j2 + 1];
        a[j] += a[j2];
        a[j + 1] += a[j2 + 1];
        x3r = a[j1] - a[j3];
        x3i = a[j1 + 1] - a[j3 + 1];
        a[j1] += a[j3];
        a[j1 + 1] += a[j3 + 1];
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j2] = wk1r * x0r + wk1i * x0i;
        a[j2 + 1] = wk1r * x0i - wk1i * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j3] = wk3r * x0r + wk3i * x0i;
        a[j3 + 1] = wk3r * x0i - wk3i * x0r;
        j0 = nq - j;
        j1 = j0 + nq;
        j2 = j1 + nq;
        j3 = j2 + nq;
        x1r = a[j0] - a[j2];
        x1i = a[j0 + 1] - a[j2 + 1];
        a[j0] += a[j2];
        a[j0 + 1] += a[j2 + 1];
        x3r = a[j1] - a[j3];
        x3i = a[j1 + 1] - a[j3 + 1];
        a[j1] += a[j3];
        a[j1 + 1] += a[j3 + 1];
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j2] = wk1i * x0r + wk1r * x0i;
        a[j2 + 1] = wk1i * x0i - wk1r * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j3] = -wk3i * x0r - wk3r * x0i;
        a[j3 + 1] = -wk3i * x0i + wk3r * x0r;
    }
    j0 = nqh;
    j1 = j0 + nq;
    j2 = j1 + nq;
    j3 = j2 + nq;
    x1r = a[j0] - a[j2];
    x1i = a[j0 + 1] - a[j2 + 1];
    a[j0] += a[j2];
    a[j0 + 1] += a[j2 + 1];
    x3r = a[j1] - a[j3];
    x3i = a[j1 + 1] - a[j3 + 1];
    a[j1] += a[j3];
    a[j1 + 1] += a[j3 + 1];
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j2] = wn4r * (x0r + x0i);
    a[j2 + 1] = wn4r * (x0i - x0r);
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j3] = -wn4r * (x0r - x0i);
    a[j3 + 1] = -wn4r * (x0i + x0r);
}


void cftffix(int n, double wn4r, double *a)
{
    int j, k;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    for (k = 16; k <= n; k <<= 2) {
        for (j = k - 16; j < n; j += 2 * k) {
            x0r = a[j] - a[j + 8];
            x0i = a[j + 1] - a[j + 9];
            a[j] += a[j + 8];
            a[j + 1] += a[j + 9];
            a[j + 8] = x0r;
            a[j + 9] = x0i;
            x0r = a[j + 2] - a[j + 10];
            x0i = a[j + 3] - a[j + 11];
            a[j + 2] += a[j + 10];
            a[j + 3] += a[j + 11];
            a[j + 10] = wn4r * (x0r - x0i);
            a[j + 11] = wn4r * (x0i + x0r);
            x0r = a[j + 4] - a[j + 12];
            x0i = a[j + 5] - a[j + 13];
            a[j + 4] += a[j + 12];
            a[j + 5] += a[j + 13];
            a[j + 12] = -x0i;
            a[j + 13] = x0r;
            x0r = a[j + 6] - a[j + 14];
            x0i = a[j + 7] - a[j + 15];
            a[j + 6] += a[j + 14];
            a[j + 7] += a[j + 15];
            a[j + 14] = -wn4r * (x0r + x0i);
            a[j + 15] = -wn4r * (x0i - x0r);
        }
    }
    for (j = 0; j < n; j += 8) {
        x0r = a[j] + a[j + 4];
        x0i = a[j + 1] + a[j + 5];
        x1r = a[j] - a[j + 4];
        x1i = a[j + 1] - a[j + 5];
        x2r = a[j + 2] + a[j + 6];
        x2i = a[j + 3] + a[j + 7];
        x3r = a[j + 2] - a[j + 6];
        x3i = a[j + 3] - a[j + 7];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j + 2] = x0r - x2r;
        a[j + 3] = x0i - x2i;
        a[j + 4] = x1r - x3i;
        a[j + 5] = x1i + x3r;
        a[j + 6] = x1r + x3i;
        a[j + 7] = x1i - x3r;
    }
}


void cftbfix(int n, double wn4r, double *a)
{
    int j, k;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    for (k = 16; k <= n; k <<= 2) {
        for (j = k - 16; j < n; j += 2 * k) {
            x0r = a[j] - a[j + 8];
            x0i = a[j + 1] - a[j + 9];
            a[j] += a[j + 8];
            a[j + 1] += a[j + 9];
            a[j + 8] = x0r;
            a[j + 9] = x0i;
            x0r = a[j + 2] - a[j + 10];
            x0i = a[j + 3] - a[j + 11];
            a[j + 2] += a[j + 10];
            a[j + 3] += a[j + 11];
            a[j + 10] = wn4r * (x0r + x0i);
            a[j + 11] = wn4r * (x0i - x0r);
            x0r = a[j + 4] - a[j + 12];
            x0i = a[j + 5] - a[j + 13];
            a[j + 4] += a[j + 12];
            a[j + 5] += a[j + 13];
            a[j + 12] = x0i;
            a[j + 13] = -x0r;
            x0r = a[j + 6] - a[j + 14];
            x0i = a[j + 7] - a[j + 15];
            a[j + 6] += a[j + 14];
            a[j + 7] += a[j + 15];
            a[j + 14] = -wn4r * (x0r - x0i);
            a[j + 15] = -wn4r * (x0i + x0r);
        }
    }
    for (j = 0; j < n; j += 8) {
        x0r = a[j] + a[j + 4];
        x0i = a[j + 1] + a[j + 5];
        x1r = a[j] - a[j + 4];
        x1i = a[j + 1] - a[j + 5];
        x2r = a[j + 2] + a[j + 6];
        x2i = a[j + 3] + a[j + 7];
        x3r = a[j + 2] - a[j + 6];
        x3i = a[j + 3] - a[j + 7];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j + 2] = x0r - x2r;
        a[j + 3] = x0i - x2i;
        a[j + 4] = x1r + x3i;
        a[j + 5] = x1i - x3r;
        a[j + 6] = x1r - x3i;
        a[j + 7] = x1i + x3r;
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


void rftbsub(int n, double *a, int nc, double *c)
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
        yr = wkr * xr + wki * xi;
        yi = wkr * xi - wki * xr;
        a[k] -= yr;
        a[k + 1] -= yi;
        a[j] += yr;
        a[j + 1] -= yi;
    }
}


void dctsub(int n, double *a, int nc, double *c)
{
    int j, k, kk, ks, m;
    double wkr, wki, xr;
    
    ks = nc / n;
    kk = 0;
    m = n >> 1;
    for (k = 1; k < m; k++) {
        j = n - k;
        kk += ks;
        wkr = c[kk] - c[nc - kk];
        wki = c[kk] + c[nc - kk];
        xr = wki * a[k] - wkr * a[j];
        a[k] = wkr * a[k] + wki * a[j];
        a[j] = xr;
    }
    a[m] *= c[0];
}


void dstsub(int n, double *a, int nc, double *c)
{
    int j, k, kk, ks, m;
    double wkr, wki, xr;
    
    ks = nc / n;
    kk = 0;
    m = n >> 1;
    for (k = 1; k < m; k++) {
        j = n - k;
        kk += ks;
        wkr = c[kk] - c[nc - kk];
        wki = c[kk] + c[nc - kk];
        xr = wki * a[j] - wkr * a[k];
        a[j] = wkr * a[j] + wki * a[k];
        a[k] = xr;
    }
    a[m] *= c[0];
}


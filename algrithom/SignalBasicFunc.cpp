#include"SignalBasicFunc.h"
#include<math.h>
#include<assert.h>
#include "_kiss_fft_guts.h"
#include "kiss_fftr.h"
#include<stdlib.h>
#include <stdlib.h>

#ifndef max   
#define max(A, B) ((A)>(B) ? (A) : (B))   
#endif   
#ifndef min   
#define min(A, B) ((A)<(B) ? (A) : (B))   
#endif 

using namespace std;
void SignalBasicFunc::hamming(int n, double* w)
{
	int i;
	double k = 2 * PI / (n - 1);   /* 2*pi/(N-1) */

	for (i = 0; i < n; i++)
		*w++ = 0.54 - 0.46*cos(k*i);
}
 void SignalBasicFunc::hanning(int n, double* w) {

	int i;
	double k = 2.0 * PI / ((double)n);
	for (i = 0; i < n; i++)
		*w++ = 0.5 - 0.5*cos(k*(double)i);
}
void SignalBasicFunc::filter(int ord, double *a, double *b, int np, double *x, double *y)
{
	int i, j;
	y[0] = b[0] * x[0];
	for (i = 1; i < ord + 1; i++)
	{
		y[i] = 0.0;
		for (j = 0; j < i + 1; j++)
			y[i] = y[i] + b[j] * x[i - j];
		for (j = 0; j < i; j++)
			y[i] = y[i] - a[j + 1] * y[i - j - 1];
	}

	for (i = ord + 1; i < np + 1; i++)
	{
		y[i] = 0.0;
		for (j = 0; j < ord + 1; j++)
			y[i] = y[i] + b[j] * x[i - j];
		for (j = 0; j < ord; j++)
			y[i] = y[i] - a[j + 1] * y[i - j - 1];
	}
	return;
}
void SignalBasicFunc::fconv(double *x, int Lx, double *h, int Lh, double *y) {
	double *temp, max = 1, *t_x, *t_h;
	int Ly = Lx + Lh - 1, Ly2 = int(pow(2, int(round(log(Ly) / log(2)))) * 2), i;
	kiss_fft_cpx *X, *H, *Y;
	kiss_fftr_cfg kiss_fftr_state;
	X = (kiss_fft_cpx*)calloc(Ly2, sizeof(kiss_fft_cpx)); assert(X);
	H = (kiss_fft_cpx*)calloc(Ly2, sizeof(kiss_fft_cpx)); assert(H);
	Y = (kiss_fft_cpx*)calloc(Ly2, sizeof(kiss_fft_cpx)); assert(Y);
	temp = (double*)calloc(Ly2, sizeof(double)); assert(temp);
	t_x = (double*)calloc(Ly2, sizeof(double)); assert(t_x);
	t_h = (double*)calloc(Ly2, sizeof(double)); assert(t_h);
	memcpy(t_x, x, sizeof(double)*Lx);
	memcpy(t_h, h, sizeof(double)*Lh);
	kiss_fftr_state = kiss_fftr_alloc(Ly2, 0, 0, 0);
	kiss_fftr(kiss_fftr_state, t_x, X);
	if (kiss_fftr_state != NULL) {

		free(kiss_fftr_state); kiss_fftr_state = NULL;
	}


	kiss_fftr_state = kiss_fftr_alloc(Ly2, 0, 0, 0);
	kiss_fftr(kiss_fftr_state, t_h, H);

	for (i = 0; i < Ly2; i++)
		C_MUL(Y[i], X[i], H[i]);

	if (kiss_fftr_state != NULL) {
		free(kiss_fftr_state); kiss_fftr_state = NULL;
	}

	kiss_fftr_state = kiss_fftr_alloc(Ly2, 1, 0, 0);
	kiss_fftri(kiss_fftr_state, Y, temp);
	for (i = 0; i < Ly2; i++)
	{
		*(temp + i) /= Ly2;
		if (max < fabs(*(temp + i))) max = fabs(*(temp + i));
	}
	for (i = 0; i < Ly2; i++)
		*(temp + i) /= (max + 0.01);
	memcpy(y, temp, sizeof(double)*Ly);
	if (NULL != X) {
		free(X); X = NULL;
	}
	if (NULL != Y) { free(Y); Y = NULL; }
	if (NULL != H) { free(H); H = NULL; }
	if (NULL != temp) { free(temp); temp = NULL; }
	if (NULL != kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }

	if (NULL != t_x) { free(t_x); t_x = NULL; }
	if (NULL != t_h) { free(t_h); t_h = NULL; }
}

double SignalBasicFunc::first_modified_Bessel(int n, double x)
{
	int i, m;
	double t, y, p, b0, b1, q;
	static double a[7] = { 1.0, 3.5156229, 3.0899424, 1.2067492,
		0.2659732, 0.0360768, 0.0045813 };
	static double b[7] = { 0.5, 0.87890594, 0.51498869,
		0.15084934, 0.02658773, 0.00301532, 0.00032411 };
	static double c[9] = { 0.39894228, 0.01328592, 0.00225319,
		-0.00157565, 0.00916281, -0.02057706,
		0.02635537, -0.01647633, 0.00392377 };
	static double d[9] = { 0.39894228, -0.03988024, -0.00362018,
		0.00163801, -0.01031555, 0.02282967,
		-0.02895312, 0.01787654, -0.00420059 };
	if (n<0) n = -n;
	t = fabs(x);
	if (n != 1)
	{
		if (t<3.75)
		{
			y = (x / 3.75)*(x / 3.75); p = a[6];
			for (i = 5; i >= 0; i--)
				p = p*y + a[i];
		}
		else
		{
			y = 3.75 / t; p = c[8];
			for (i = 7; i >= 0; i--)
				p = p*y + c[i];
			p = p*exp(t) / sqrt(t);
		}
	}
	if (n == 0) return(p);
	q = p;
	if (t<3.75)
	{
		y = (x / 3.75)*(x / 3.75); p = b[6];
		for (i = 5; i >= 0; i--) p = p*y + b[i];
		p = p*t;
	}
	else
	{
		y = 3.75 / t; p = d[8];
		for (i = 7; i >= 0; i--) p = p*y + d[i];
		p = p*exp(t) / sqrt(t);
	}
	if (x<0.0) p = -p;
	if (n == 1) return(p);
	if (x == 0.0) return(0.0);
	y = 2.0 / t; t = 0.0; b1 = 1.0; b0 = 0.0;
	m = n + (int)sqrt(40.0*n);
	m = 2 * m;

	for (i = m; i>0; i--)
	{
		p = b0 + i*y*b1; b0 = b1; b1 = p;
		if (fabs(b1)>1.0e+10)
		{
			t = t*1.0e-10; b0 = b0*1.0e-10;
			b1 = b1*1.0e-10;
		}
		if (i == n) t = b0;
	}
	p = t*q / b1;
	if ((x<0.0) && (n % 2 == 1)) p = -p;
	return(p);
}

double SignalBasicFunc::sinc(double x) {
	if (0.0f == x)return 1;
	else return sin(PI*x) / (PI*x);
}
void SignalBasicFunc::kaiser(short N, double beta,double *window) {
	double theta;
	short n;
	
	for (n = 0; n < N; n++)
	{

		theta = beta * sqrtf(1 - powf((2 * (double)n / (N - 1) - 1), 2));
		*(window + n) = first_modified_Bessel(0, theta) / first_modified_Bessel(0, beta);
	}

}


void SignalBasicFunc::firls(int N, double *F, short *M,double *h) {
	double dF[3], *k, *a, *b, m, b1;
	int L = (N - 1) / 2, i, s;
	short W[2] = { 1, 1 }, Nodd = N % 2;
	dF[0] = F[1]; dF[1] = 0.0f; dF[2] = F[3] - F[2];
	k = (double *)malloc(sizeof(double)*(L + 1)); assert(k);
	a = (double *)malloc(sizeof(double)*(L + 1)); assert(a);
	b = (double *)malloc(sizeof(double)*(L + 1)); assert(b);
	//h = (double *)malloc(sizeof(double)*(N));		assert(h);
	for (i = 0; i <= L; i++)
	{
		*(k + i) = (double)i;
		*(b + i) = 0.0f;
	}
	if (0 == Nodd)
	{
		for (i = 0; i <= L; i++)
			*(k + i) += 0.5f;
	}
	for (s = 0; s < 4; s += 2) {
		m = (M[s + 1] - M[s]) / (F[s + 1] - F[s]);
		b1 = M[s] - m*F[s];
		if (1 == Nodd) b[0] += (b1*(F[s + 1] - F[s]) + m / 2 * (F[s + 1] * F[s + 1] - F[s] * F[s]));
		else
		{
			*b += (m / (4 * PI*PI)*(cos(2 * PI*(*k) * F[s + 1]) - cos(2 * PI*(*k) * F[s])) / ((*k) * (*k)));
			*b += (F[s + 1] * (m*F[s + 1] + b1)* sinc(2 * (*k) * F[s + 1]) - F[s] * (m*F[s] + b1)*sinc(2 * (*k) * F[s]));
		}
		for (i = 1; i <= L; i++)
		{
			*(b + i) += (m / (4 * PI*PI)*(cos(2 * PI*(*(k + i)) * F[s + 1]) - cos(2 * PI*(*(k + i)) * F[s])) / ((*(k + i)) * (*(k + i))));
			*(b + i) += (F[s + 1] * (m*F[s + 1] + b1)* sinc(2 * (*(k + i)) * F[s + 1]) - F[s] * (m*F[s] + b1)*sinc(2 * (*(k + i)) * F[s]));
			;
		}
	}
	for (i = 0; i <= L; i++)
		*(a + i) = (W[0] * W[0]) * 4 * (*(b + i));
	if (1 == Nodd)
	{

		for (i = 0; i <= L; i++)
			*(h + i) = *(a + (L - i)) / 2;
		for (i = L + 1; i < N; i++)
			*(h + i) = *(a + (i - L)) / 2;
	}
	else {
		for (i = 0; i <= L; i++)
			*(h + i) = *(a + (L - i)) / 2;
		for (i = L + 1; i < N; i++)
			*(h + i) = *(a + (i - L - 1)) / 2;
	}
	if (a != NULL)
	{
		free(a);
		a = NULL;
	}
	if (b != NULL)
	{
		free(b); b = NULL;
	}
	if (k != NULL)
	{
		free(k); k = NULL;
	}

}
void SignalBasicFunc::upfirdn(
	double y[], unsigned int Ly, unsigned int ky,
	double x[], unsigned int Lx, unsigned int kx,
	double h[], unsigned int Lh, unsigned int kh,
	int p,
	int q
)
{
	int r, rpq_offset, k, Lg;
	int iv, ig, igv, iw;
	double  *pw, acc;
	double  *pv, *pvend;
	double  *pvhi, *pvlo, *pvt;
	double  *pg, *pgend;
	double  *pghi, *pglo, *pgt;
	unsigned int kmax = max(kh, kx);

	iv = q;
	ig = iw = p;
	igv = p*q;

	for (k = 0; k<kmax; k++) {
		pvend = x + Lx;
		pgend = h + Lh;

		for (r = 0; r<p; r++) {
			pw = y + r;
			pg = h + ((r*q) % p);
			Lg = pgend - pg;
			Lg = (Lg%p) ? Lg / p + 1 : Lg / p;
			rpq_offset = (r*q) / p;
			pv = x + rpq_offset;

			/*
			* PSEUDO-CODE for CONVOLUTION with GENERAL INCREMENTS:
			*
			*   w[n] = v[n] * g[n]
			*
			* Given:
			*   pointers:   pg, pv, and pw
			*   or arrays:  g[ ], v[ ], and w[ ]
			*   increments: ig, iv, and iw
			*   end points: h+Lh, x+Lx
			*/

			/*
			* Region #1 (running onto the data):
			*/
			pglo = pg;
			pghi = pg + p*rpq_offset;
			pvlo = x;
			pvhi = pv;
			while ((pvhi<pvend) && (pghi<pgend)) {
				acc = 0.0;
				pvt = pvhi;
				pgt = pglo;
				while (pgt <= pghi) {
					//printf("%f %f", *pgt, *pvt);
					acc += (*pgt) * (*pvt--);
					pgt += ig;
				}
				*pw += acc;
				//printf("%.5f ", *pw);
				pw += iw;
				pvhi += iv;
				pghi += igv;
			}


			//Do we need to drain rest of signal?

			if (pvhi < pvend) {

				// Region #2 (complete overlap):

				while (pghi >= pgend) {
					pghi -= ig;
				}
				while (pvhi < pvend) {
					acc = 0.0;
					pvt = pvhi;
					pgt = pglo;
					while (pgt <= pghi) {
						//printf(*pgt)
						acc += (*pgt) * (*pvt--);
						pgt += ig;
					}
					*pw += acc;
					pw += iw;
					pvhi += iv;
				}

			}
			else if (pghi < pgend) {

				// Region #2a (drain out the filter):

				while (pghi < pgend) {
					acc = 0.0;
					pvt = pvlo;     // pvlo is still equal to x 
					pgt = pghi;
					while (pvt < pvend) {
						acc += (*pgt) * (*pvt++);
						pgt -= ig;
					}
					*pw += acc;
					pw += iw;
					pghi += igv;
					pvhi += iv;
				}
			}

			while (pghi >= pgend) {
				pghi -= ig;
			}
			pvlo = pvhi - Lg + 1;

			while (pvlo < pvend) {

				// Region #3 (running off the data):

				acc = 0.0;
				pvt = pvlo;
				pgt = pghi;
				while (pvt < pvend) {
					acc += (*pgt) * (*pvt++);
					pgt -= ig;
				}
				*pw += acc;
				pw += iw;
				pvlo += iv;
			}

		} // end of r loop 

		if (kx != 1) x += Lx;
		if (kh != 1) h += Lh;
		y += Ly;
	}
}
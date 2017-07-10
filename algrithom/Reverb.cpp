#include"Reverb.h"
#include "_kiss_fft_guts.h"
#include "kiss_fftr.h"

void fconv(OsFlt64 *x, OsInt32 Lx, OsFlt64 *h, OsInt32 Lh, OsFlt64 *y) {
	OsFlt64 *temp, max = 1,*t_x,*t_h;
	OsInt32 Ly = Lx + Lh - 1, Ly2 =  OsInt32(pow(2, OsInt32(round(log(Ly) / log(2)))) * 2), i;
	kiss_fft_cpx *X, *H,*Y;
	kiss_fftr_cfg kiss_fftr_state;
	X = (kiss_fft_cpx*)calloc(Ly2, sizeof(kiss_fft_cpx)); assert(X);
	H = (kiss_fft_cpx*)calloc(Ly2, sizeof(kiss_fft_cpx)); assert(H);
	Y = (kiss_fft_cpx*)calloc(Ly2, sizeof(kiss_fft_cpx)); assert(Y);
	temp = (OsFlt64*)calloc(Ly2, sizeof(OsFlt64)); assert(temp);
	t_x = (OsFlt64*)calloc(Ly2, sizeof(OsFlt64)); assert(t_x);
	t_h = (OsFlt64*)calloc(Ly2, sizeof(OsFlt64)); assert(t_h);
	memcpy(t_x, x, sizeof(OsFlt64)*Lx);
	memcpy(t_h, h, sizeof(OsFlt64)*Lh);
	kiss_fftr_state = kiss_fftr_alloc(Ly2, 0, 0, 0);
	kiss_fftr(kiss_fftr_state, t_x, X);
	if (kiss_fftr_state != NULL) {

		free(kiss_fftr_state); kiss_fftr_state = NULL;
	}

	
	kiss_fftr_state = kiss_fftr_alloc(Ly2, 0, 0, 0);
	kiss_fftr(kiss_fftr_state, t_h, H);
	
	for (i = 0; i < Ly2; i++) 
		 C_MUL(Y[i],X[i], H[i]);
	
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
		*(temp + i) /= (max+0.01);
	memcpy(y, temp, sizeof(OsFlt64)*Ly);
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

void convReverb(OsFlt64 *x, OsInt32 Lx, OsFlt64 *h, OsInt32 Lh, OsFlt64 *y) {
	fconv(x, Lx, h, Lh, y);
}
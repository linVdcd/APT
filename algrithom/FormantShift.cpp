#include"FormantShift.h"
#include "_kiss_fft_guts.h"
#include "kiss_fftr.h"
#include "kiss_fft.h"
#include "SignalBasicFunc.h"
#include<assert.h>
void FS::FormantWarp(double *input, const int Lx, double **DAFx_out, const double warping_coef) {
	kiss_fft_cpx *dft,*flog;
	kiss_fftr_cfg kiss_fftr_state;
	double *window, *cep,*cep_cut, *DAFx_in, *grain,*flog_cut1,*flog_cut2;
	int i, pin = 0, pout = 0, pend, tLx = Lx;
	short outwin = 160, NFFT = 400,*x0,*x,hw = NFFT/2,order=50;
	SignalBasicFunc *sbf = new SignalBasicFunc();
	dft = (kiss_fft_cpx *)calloc(NFFT, sizeof(kiss_fft_cpx));				assert(dft);
	window = (double*)calloc(NFFT, sizeof(double));						assert(window);
	grain = (double*)calloc(NFFT, sizeof(double));						assert(grain);
	x0 = (short*)calloc(hw + 1, sizeof(short)); assert(x0);
	x = (short*)calloc(NFFT, sizeof(short)); assert(x);
	flog = (kiss_fft_cpx *)calloc(NFFT, sizeof(kiss_fft_cpx));				assert(flog);
	cep = (double*)calloc(NFFT, sizeof(double));						assert(cep);
	cep_cut = (double*)calloc(NFFT, sizeof(double));						assert(cep_cut);
	flog_cut1 = (double*)calloc(NFFT, sizeof(double));						assert(flog_cut1);
	flog_cut2 = (double*)calloc(NFFT, sizeof(double));						assert(flog_cut2);
	for (i = 0; i <=hw; i++)
	{
		x0[i] = floor(fmin(i / warping_coef, (double)hw));
		x[i] = x0[i];
	}
	for (i = hw - 1; i <0; i--)
		x[NFFT-i] = x[i];
	sbf->hanning(NFFT, window);


	tLx += NFFT + NFFT - Lx%outwin;
	DAFx_in = (double*)malloc(sizeof(double)*tLx); assert(DAFx_in);
	*DAFx_out = (double*)malloc(sizeof(double)*tLx); assert(DAFx_out);

	memset(*DAFx_out, 0, sizeof(double)*tLx);
	for (i = 0; i < NFFT; i++)
		*(DAFx_in + i) = 0.0;
	for (i = 0; i < Lx; i++)
		*(DAFx_in + i + NFFT) = *(input + i);
	for (i = NFFT + Lx; i < tLx; i++)
		*(DAFx_in + i) = 0.0;

	kiss_fftr_state = kiss_fftr_alloc(NFFT, 0, 0, 0);

	pend = tLx - NFFT;

	while (pin < pend)
	{
		memcpy(grain, DAFx_in + pin, sizeof(double)*NFFT);
		for (i = 0; i < NFFT; i++)
			*(grain + i) *= *(window + i);
		if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }
		kiss_fftr_state = kiss_fftr_alloc(NFFT, 0, 0, 0);
		kiss_fftr(kiss_fftr_state, grain, dft);
		for (i = 0; i < NFFT; i++) {
			double r = dft[i].r / hw,imag = dft[i].i/hw;
			flog[i].r = log(0.00001 + sqrt(r*r+imag*imag));
			flog[i].i = 0;
		}
		if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }
		kiss_fftr_state = kiss_fftr_alloc(NFFT, 1, 0, 0);
		kiss_fftri(kiss_fftr_state, flog, cep);

		cep_cut[0] = cep[0] / (2*NFFT);
		for (i = 1; i < NFFT; i++)
		{
			if (i < order)
				cep_cut[i] = cep[i]/NFFT;
			else
				cep_cut[i] = 0.0;
		}
		
		

		if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }
		kiss_fftr_state = kiss_fftr_alloc(NFFT, 0, 0, 0);
		kiss_fftr(kiss_fftr_state, cep_cut, flog);
		
		for (i = 0; i < NFFT; i++)
		{
			flog_cut1[i] = 2 * flog[i].r;
			flog_cut2[i] = flog_cut1[x[i]];
		}
		
		for (i = 0; i < NFFT; i++)
		{
			flog_cut2[i] = flog_cut1[x[i]];
			dft[i].r *= exp(flog_cut2[i] - flog_cut1[i]);
			dft[i].i *= exp(flog_cut2[i] - flog_cut1[i]);
		}
		
		

		if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }
		kiss_fftr_state = kiss_fftr_alloc(NFFT, 1, 0, 0);
		kiss_fftri(kiss_fftr_state, dft, grain);


		
		for (i = 0; i < NFFT; i++)
		{
			*(grain + i) /= NFFT;
			*(grain + i) *= *(window + i);
		}
		for (i = pout; i < pout + NFFT; i++)
			*(*DAFx_out + i) += *(grain + i - pout);

		pin += outwin;
		pout += outwin;
	}
	


	if (x0) { free(x0); x0 = NULL; }
	if (x) { free(x); x = NULL; }
	if (cep_cut) { free(cep_cut); cep_cut = NULL; }
	if (flog_cut1) { free(flog_cut1); flog_cut1 = NULL; }
	if (flog_cut2) { free(flog_cut2); flog_cut2 = NULL; }
	if (flog) { free(flog); flog = NULL; }
	if (cep) { free(cep); cep = NULL; }
	if (NULL != DAFx_in)
	{
		free(DAFx_in); DAFx_in = NULL;
	}
	if (NULL != grain)
	{
		free(grain); grain = NULL;
	}
	if (NULL != window)
	{
		free(window); window = NULL;
	}
	if (dft)
	{
		free(dft); dft = NULL;
	}
	if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }

	delete sbf;
}

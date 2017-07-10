/*-----------------------------------
	林明安 2016.6.1
	两种机器人音效，1、基于环形调制；2、基于相位置零。
*/
#include "wavRandW.h"
#include "_kiss_fft_guts.h"
#include "kiss_fftr.h"
#include "kiss_fft.h"
#include <complex> 
#include "robot.h"
#include "SignalBasicFunc.h"
#include "lpc.h"
using namespace std;
#define kiss_fft_scalar OsFlt64


void RobotRM(const short *x, int Lx, short *y, const int fs,const double RMfre){
	//利用环形调制算法生成机器人声音
	OsInt32 i;
	for (i = 0; i < Lx; i++)
		*(y + i) = round(*(x + i) * (sin(2 * M_PI*(i + 1)*(RMfre / fs))));
	
}
void RobotZP(const short *x, int Lx, short *y,const int outwin){
	//相位置零生成机器人声音
	//_C_double_complex *out_dft;
	kiss_fft_cpx *dft;
	kiss_fftr_cfg kiss_fftr_state;
	OsFlt64 *window, *sp, *DAFx_in,*DAFx_out,*grain;
	OsInt32 i,pin=0,pout=0,pend,tLx = Lx;
	OsInt16 NFFT = 1024;
	SignalBasicFunc *sbf = new SignalBasicFunc();
	dft = (kiss_fft_cpx *)calloc(NFFT, sizeof(kiss_fft_cpx));				assert(dft);
	window = (OsFlt64*)calloc(NFFT,sizeof(OsFlt64));						assert(window);
	grain = (OsFlt64*)calloc(NFFT,sizeof(OsFlt64));						assert(grain);

	sbf->hanning(NFFT, window);
	
	
	Lx += NFFT + NFFT - Lx%outwin;
	DAFx_in = (OsFlt64*)malloc(sizeof(OsFlt64)*Lx);
	DAFx_out = (OsFlt64*)malloc(sizeof(OsFlt64)*Lx);
	
	memset(DAFx_out, 0, sizeof(OsFlt64)*Lx);
	for (i = 0; i < NFFT; i++)
			*(DAFx_in + i) = 0.0;
	for (i = 0; i < tLx; i++)
		*(DAFx_in + i + NFFT) = *(x + i) / 32768.0;
	for (i = NFFT + tLx; i < Lx;i++)
		*(DAFx_in + i) = 0.0;
	//if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }
	kiss_fftr_state = kiss_fftr_alloc(NFFT, 0, 0, 0);

	pend = Lx - NFFT;

	while (pin < pend)
	{
		memcpy(grain, DAFx_in + pin, sizeof(OsFlt64)*NFFT);
		for (i = 0; i < NFFT; i++)
			*(grain + i) *= *(window + i);
		if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }
		kiss_fftr_state = kiss_fftr_alloc(NFFT, 0, 0, 0);
		kiss_fftr(kiss_fftr_state, grain, dft);
		for (i = 0; i < NFFT; i++){
			dft[i].r = sqrt(dft[i].r*dft[i].r + dft[i].i*dft[i].i);
			dft[i].i = 0;
		}
		if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }
		kiss_fftr_state = kiss_fftr_alloc(NFFT, 1, 0, 0);
		kiss_fftri(kiss_fftr_state, dft, grain);
		for (i = 0; i < NFFT; i++)
		{
			*(grain + i) /= NFFT;
			*(grain + i) *= *(window + i);
		}
		for (i = pout; i < pout+NFFT; i++)
			*(DAFx_out + i) += *(grain + i - pout);
		
		pin += outwin;
		pout += outwin;
	}
	for (i = NFFT; i < tLx +NFFT; i++)
		*(y + i - NFFT) = round(*(DAFx_out + i) * 32768.0);
	
	if (NULL != DAFx_out)
	{
		free(DAFx_out); DAFx_out = NULL;
	}
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

void RobotSound::Robot(const short *x, int Lx, short *y,const short type,const int fs,const double RMfre,const int outwin) {
	if (type == 1)
		RobotRM(x, Lx, y, fs, RMfre);
	else
		RobotZP(x, Lx, y, outwin);
}

void RobotSound::WhisperSound(const short *x, int Lx, short **y, int &Ly,const int p) {
	double *x_f,*y_f,*a,*g,*e;
	Lpc *lpc = new Lpc();
	
	short nframe,h=128;
	x_f = (double*)calloc(Lx, sizeof(double));
	for (int i = 0; i < Lx; i++)
		*(x_f + i) = *(x+i)/32768.0;
	lpc->lpcfit(x_f, Lx, p, h, &a, &g, &e, &nframe);
	for (int i = 0; i < (nframe + 1)*h - h / 2; i++) {
		*(e + i) = (double)rand() / RAND_MAX;
	}

	lpc->lpcsynth(a, g, e, p, h, nframe, &y_f,Ly);

	*y = (short*)calloc(Ly, sizeof(short)); assert(y);
	for (int i = 0; i < Ly; i++)
		* (*y+i) = y_f[i]*32768.0;

	if (a) { free(a); a = NULL; }
	if (g) { free(g); g = NULL; }
	if (e) { free(e); e = NULL; }
	if (x_f) { free(x_f); x_f = NULL; }
	if (y_f) { free(y_f); y_f = NULL; }
	delete lpc;
}

RobotSound::RobotSound(){}
RobotSound::~RobotSound(){
	
}
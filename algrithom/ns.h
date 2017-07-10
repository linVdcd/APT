#ifndef NS_H
#define NS_H

#include "Platform.h"



typedef struct MCRAcf{
	OsInt16 count, IS;
	OsFlt32 alpha_s, Bmin, Gamma0, Gamma1, zeta0, alpha_d;
	OsFlt32 *S, *tilde_S, *Smin, *tilde_Smin, *Smin_sw, *tilde_Smin_sw, *ov_Lambda_D, *P, *Sf, *tilde_Sf;

};
typedef struct param_martin {
	OsFlt32 *alpha, *p, *Pbar, *Psqbar, *actmin, *actmin_sub, *Pmin_u, **mincat_val, *lmin_flag, Y_M[14], Y_H[14],
		alpha_corr,Av,alpha_max,alpha_min,beta_max,M_D,H_D,M_V,H_V;
	OsInt32 subwc, u, L, R, D, Um,iter, V;

};
struct NoiseReduction
{
	OsInt32 sample_rate;
	OsInt32 frame_size;
	OsInt32 ps_size;        // size of power spectrum

	OsFlt32 *frame;
	OsFlt32 *window;        // window function
	OsFlt32 *ft;            // FFT coefficient
	OsFlt32 *ps;            // power spectrum
	OsFlt64 theta;
	OsFlt32 *noise;         // noise power spectrum

	OsFlt32 *inbuf;         // overlapped add analysis
	OsFlt32 *outbuf;
	OsFlt32 win_gain;       // gain by window function
	void    *fft_table;
	OsInt32 adapt_count;

	/*------- soft mask estimator with posteriori SNR -------*/

	OsFlt32 *X2, *X2_old;
	OsFlt32 *noise_mu2_old;
	OsFlt32 *gammak, *gammak_old;
	OsFlt32 *ksi, *ksi_old;
	struct MCRAcf cf;
	struct param_martin params;
};


struct NoiseReduction* NoiseReductionCreate(OsInt32 inSampleRate,OsInt32 inFrameSize);
struct NoiseReduction* NoiseReductionCreate_martin(OsInt32 inSampleRate, OsInt32 inFrameSize,OsFlt64 theta);
void NoiseReductionDestroy(struct NoiseReduction *inHandle);
void InitNoise(struct NoiseReduction *inHandle, OsInt16 *ioPcm);
void NoiseReductionProcess(struct NoiseReduction *inHandle,OsInt16 *ioPcm);
void NoiseReductionProcess_martin(struct NoiseReduction *inHandle, OsInt16 *ioPcm);


#endif
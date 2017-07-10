#include "ns.h"
#include "fftwrap.h"
#include "_kiss_fft_guts.h"
#include "kiss_fftr.h"
#include "kiss_fft.h"
/*
 * Hamming
 *                        2*pi*k
 * w(k) = 0.54 - 0.46*cos(------), where 0 <= k < N
 *                         N-1
 *
 * n window length
 * w buffer for the window parameters
 */
static OS_INLINE void hamming(OsInt32 n,OsFlt32* w)
{
    OsInt32 i;
    OsFlt32 k = 2*M_PI/(n-1);   /* 2*pi/(N-1) */

    for (i = 0; i < n; i++)
        *w++ = 0.54 - 0.46*cos(k*i);
}

void RemoveDC(OsInt16 *x,OsInt32 inLen)
{
    OsInt32 i,nOffset=0;
    OsFlt32 fTmp = 0.0f;

    for(i = 0; i < inLen; i++) 
    {
        fTmp += x[i];
    }
    fTmp = -fTmp/inLen;
    //if(fTmp < 0.0f) nOffset = fTmp-0.5f;
   // else if(fTmp > 0.0f) nOffset = fTmp+0.5f;

    for(i = 0; i < inLen; i++) 
    {
        x[i] += nOffset;

        if(x[i] > 30000) x[i] = 30000;
        else if(x[i] < -30000) x[i] = -30000;
    }
}
struct NoiseReduction* NoiseReductionCreate_martin(OsInt32 inSampleRate, OsInt32 inFrameSize,OsFlt64 theta) {
	OsInt32 i,j, N;
	OsFlt32 total_win_gain;

	struct NoiseReduction *pInst = (struct NoiseReduction*)calloc(sizeof(struct NoiseReduction), 1);
	assert(0 != pInst);
	pInst->sample_rate = inSampleRate;
	pInst->frame_size = inFrameSize;
	pInst->ps_size = inFrameSize;

	N = pInst->frame_size;

	pInst->frame = (OsFlt32*)malloc(2 * N * sizeof(OsFlt32));
	pInst->window = (OsFlt32*)malloc(2 * N * sizeof(OsFlt32));
	pInst->ft = (OsFlt32*)malloc(2 * N * sizeof(OsFlt32));
	pInst->ps = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	pInst->noise = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	pInst->inbuf = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	pInst->outbuf = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	memset(pInst->inbuf, 0, N * sizeof(OsFlt32));
	memset(pInst->outbuf, 0, N * sizeof(OsFlt32));
	pInst->theta = (OsFlt32)theta;

	pInst->X2 = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	pInst->X2_old = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	memset(pInst->X2_old, 0, N * sizeof(OsFlt32));
	pInst->noise_mu2_old = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	pInst->gammak = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	pInst->ksi = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	pInst->gammak_old = (OsFlt32*)malloc(N * sizeof(OsFlt32));
	pInst->ksi_old = (OsFlt32*)malloc(N * sizeof(OsFlt32));

	pInst->fft_table = spx_fft_init(2 * N);

	hamming(2 * N, pInst->window);
	total_win_gain = 0;
	for (i = 0; i < 2 * N; i++)
		total_win_gain += pInst->window[i];
	pInst->win_gain = N / total_win_gain;
	pInst->adapt_count = 0;

	//init martin params
	pInst->params.iter =2;
	pInst->params.alpha = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(pInst->params.alpha);
	pInst->params.p = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(pInst->params.p);
	pInst->params.Pbar = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(pInst->params.Pbar);
	pInst->params.Psqbar = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(pInst->params.Psqbar);
	pInst->params.actmin = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(pInst->params.actmin);
	pInst->params.actmin_sub = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(pInst->params.actmin_sub);
	pInst->params.Pmin_u = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(pInst->params.Pmin_u);
	pInst->params.lmin_flag = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(pInst->params.lmin_flag);
	//pInst->params.noise_sp = (OsFlt32*)calloc( N, sizeof(OsFlt32)); assert(pInst->params.noise_sp);
	pInst->params.mincat_val = (OsFlt32 **)malloc(sizeof(OsFlt32*)*N);
	for(i=0;i<N;i++)
		pInst->params.mincat_val[i] = (OsFlt32*)calloc(pInst->params.iter,sizeof(OsFlt32));
	return pInst;

}

void FrequencyAnalysis(struct NoiseReduction *inHandle, OsInt16 *x)
{
	OsInt32 i;
	OsInt32 N = inHandle->ps_size;
	OsFlt32 *ps = inHandle->ps;
	OsFlt32 *ft = inHandle->ft;

	//RemoveDC(x, N);

	for (i = 0; i < N; i++)
		inHandle->frame[i] = inHandle->inbuf[i];
	for (i = 0; i < N; i++)
		inHandle->frame[i + N] = x[i];
	for (i = 0; i < N; i++)
		inHandle->inbuf[i] = x[i];

	for (i = 0; i < 2 * N; i++)
		inHandle->frame[i] *= inHandle->window[i];

	spx_fft(inHandle->fft_table, inHandle->frame, ft);

	ps[0] = ft[0] * ft[0];
	for (i = 1; i < N; i++)
		ps[i] = ft[2 * i - 1] * ft[2 * i - 1] + ft[2 * i] * ft[2 * i];
}
void init_martin_params(struct NoiseReduction *inHandle,OsInt16 *inPcm) {
	OsFlt32 y_m[] = { 0,0.26,0.48,0.58,0.61,0.668,0.705,0.762,0.8,0.841,0.865,0.89,0.9,0.91 },
		y_h[] = {0,0.15,0.48,0.78,0.98,1.55,2.0,2.3,2.52,2.9,3.25,4.0,4.1,4.1};
	OsInt32 i, j, N = inHandle->frame_size;
	inHandle->params.L = inHandle->frame_size;
	inHandle->params.R = inHandle->frame_size/2;
	inHandle->params.D = 150;
	inHandle->params.V = 5;
	inHandle->params.Um = inHandle->params.iter;
	inHandle->params.Av = 2.12;
	inHandle->params.alpha_max = 0.96;
	inHandle->params.alpha_min = 0.3;
	inHandle->params.beta_max = 0.8;
	memcpy(inHandle->params.Y_M, y_m, sizeof(OsFlt32) * 14);
	memcpy(inHandle->params.Y_H, y_h, sizeof(OsFlt32) * 14);
	inHandle->params.M_D = 0.905;
	inHandle->params.H_D = 4.1;
	inHandle->params.M_V = 0.668;
	inHandle->params.H_V = 1.55;
	inHandle->params.u = 1;
	
	FrequencyAnalysis(inHandle, inPcm);
	OsFlt32 max = 0;
	for (i = 0; i < N; i++)
		inHandle->noise[i] = inHandle->ps[i] * 0.05;
	//memcpy(inHandle->noise, inHandle->ps, sizeof(OsFlt32)*N);
	memcpy(inHandle->params.p, inHandle->noise, sizeof(OsFlt32)*N);
	memcpy(inHandle->params.Pbar, inHandle->noise, sizeof(OsFlt32)*N);
	memcpy(inHandle->params.Psqbar, inHandle->noise, sizeof(OsFlt32)*N);
	memcpy(inHandle->params.actmin, inHandle->noise, sizeof(OsFlt32)*N);
	memcpy(inHandle->params.actmin_sub, inHandle->noise, sizeof(OsFlt32)*N);
	memcpy(inHandle->params.Pmin_u, inHandle->noise, sizeof(OsFlt32)*N);
	memset(inHandle->params.lmin_flag, 0, sizeof(OsFlt32));



	for (i = 0; i < N; i++)
	{
		
		if (max < inHandle->noise[i])
			max = inHandle->noise[i];
		inHandle->params.alpha[i] = 0.96;
	}
	for (i = 0; i < inHandle->params.L; i++)
		for (j = 0; j < inHandle->params.Um; j++)
			inHandle->params.mincat_val[i][j] = max;

	
}



struct NoiseReduction* NoiseReductionCreate(OsInt32 inSampleRate,OsInt32 inFrameSize)
{
    OsInt32 i,N;
    OsFlt32 total_win_gain;

	struct NoiseReduction *pInst = (struct NoiseReduction*)calloc(sizeof(struct NoiseReduction), 1);
    assert(0 != pInst);
    pInst->sample_rate  = inSampleRate;
    pInst->frame_size   = inFrameSize;
    pInst->ps_size      = inFrameSize;

    N = pInst->frame_size;

    pInst->frame    = (OsFlt32*)malloc(2*N*sizeof(OsFlt32));
    pInst->window   = (OsFlt32*)malloc(2*N*sizeof(OsFlt32));
    pInst->ft       = (OsFlt32*)malloc(2*N*sizeof(OsFlt32));
    pInst->ps       = (OsFlt32*)malloc(N*sizeof(OsFlt32));
    pInst->noise    = (OsFlt32*)malloc(N*sizeof(OsFlt32));
    pInst->inbuf    = (OsFlt32*)malloc(N*sizeof(OsFlt32));
    pInst->outbuf   = (OsFlt32*)malloc(N*sizeof(OsFlt32));
    memset(pInst->inbuf,0,N*sizeof(OsFlt32));
    memset(pInst->outbuf,0,N*sizeof(OsFlt32));


	pInst->X2       = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->X2_old   = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	memset(pInst->X2_old, 0, N*sizeof(OsFlt32));
	pInst->noise_mu2_old = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->gammak   = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->ksi      = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->gammak_old = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->ksi_old = (OsFlt32*)malloc(N*sizeof(OsFlt32));

    pInst->fft_table = spx_fft_init(2*N);

    hamming(2*N,pInst->window);
    total_win_gain = 0;
    for(i = 0; i < 2*N; i++)
		total_win_gain += pInst->window[i];
    pInst->win_gain =N/total_win_gain;
    pInst->adapt_count = 0;


	//Init MCRA cf;
	pInst->cf.IS = 4;
	pInst->cf.count = 0;
	pInst->cf.alpha_s = 0.900f;
	pInst->cf.Bmin = 1.6600f;
	pInst->cf.Gamma0 = 4.600f;
	pInst->cf.Gamma1 = 3.0f;
	pInst->cf.zeta0 = 1.6700f;
	pInst->cf.alpha_d = 0.9500f;

	pInst->cf.S = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->cf.tilde_S = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->cf.Smin = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->cf.tilde_Smin = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->cf.Smin_sw = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->cf.tilde_Smin_sw = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->cf.ov_Lambda_D = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->cf.P = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->cf.Sf = (OsFlt32*)malloc(N*sizeof(OsFlt32));
	pInst->cf.tilde_Sf = (OsFlt32*)malloc(N*sizeof(OsFlt32));

    return pInst;
}

void NoiseReductionDestroy(struct NoiseReduction *inHandle)
{
	OsInt32 N = inHandle->ps_size;
    if(!inHandle) return;
    
	if (0 != inHandle->frame)
	{
		free(inHandle->frame); inHandle->frame = NULL;
	}
	if (0 != inHandle->window)
	{
		free(inHandle->window); inHandle->window = NULL;
	}
	if (0 != inHandle->ft)
	{
		free(inHandle->ft); inHandle->ft = NULL;
	}
	if (0 != inHandle->ps)
	{
		free(inHandle->ps); inHandle->ps = NULL;
	}
	if (0 != inHandle->noise)
	{
		free(inHandle->noise); inHandle->noise = NULL;
	}
	if (0 != inHandle->inbuf)
	{
		free(inHandle->inbuf); inHandle->inbuf = NULL;
	}
	if (0 != inHandle->outbuf)
	{
		free(inHandle->outbuf); inHandle->outbuf = NULL;
	}

	if (0 != inHandle->noise_mu2_old)
	{
		free(inHandle->noise_mu2_old); inHandle->noise_mu2_old = NULL;
	}
	if (0 != inHandle->X2)
	{
		free(inHandle->X2); inHandle->X2 = NULL;
	}
	if (0 != inHandle->X2_old)
	{
		free(inHandle->X2_old); inHandle->X2_old = NULL;
	}
	if (0 != inHandle->gammak)
	{
		free(inHandle->gammak); inHandle->gammak = NULL;
	}
	if (0 != inHandle->ksi)
	{
		free(inHandle->ksi); inHandle->ksi = NULL;
	}
	if (0 != inHandle->gammak_old)
	{
		free(inHandle->gammak_old); inHandle->gammak_old = NULL;
	}
	if (0 != inHandle->ksi_old)
	{
		free(inHandle->ksi_old); inHandle->ksi_old = NULL;
	}

	if (0 != inHandle->cf.S)
	{
		free(inHandle->cf.S); inHandle->cf.S = NULL;
	}
	if (0 != inHandle->cf.tilde_S)
	{
		free(inHandle->cf.tilde_S); inHandle->cf.tilde_S = NULL;
	}
	if (0 != inHandle->cf.Smin)
	{
		free(inHandle->cf.Smin); inHandle->cf.Smin = NULL;
	}
	if (0 != inHandle->cf.tilde_Smin)
	{
		free(inHandle->cf.tilde_Smin); inHandle->cf.tilde_Smin = NULL;
	}
	if (0 != inHandle->cf.Smin_sw)
	{
		free(inHandle->cf.Smin_sw); inHandle->cf.Smin_sw = NULL;
	}
	if (0 != inHandle->cf.tilde_Smin_sw)
	{
		free(inHandle->cf.tilde_Smin_sw); inHandle->cf.tilde_Smin_sw = NULL;
	}
	if (0 != inHandle->cf.ov_Lambda_D)
	{
		free(inHandle->cf.ov_Lambda_D); inHandle->cf.ov_Lambda_D = NULL;
	}
	if (0 != inHandle->cf.P)
	{
		free(inHandle->cf.P); inHandle->cf.P = NULL;
	}
	if (0 != inHandle->cf.Sf)
	{
		free(inHandle->cf.Sf); inHandle->cf.Sf = NULL;
	}
	if (0 != inHandle->cf.tilde_Sf)
	{
		free(inHandle->cf.tilde_Sf); inHandle->cf.tilde_Sf = NULL;
	}
	

	if (NULL != inHandle->params.alpha)
	{
		free(inHandle->params.alpha); inHandle->params.alpha = NULL;
	}
	if (0 != inHandle->params.p) { free(inHandle->params.p); inHandle->params.p = NULL;
	}
	if (0 != inHandle->params.Pbar) { free(inHandle->params.Pbar); inHandle->params.Pbar = NULL;
	}
	if (0 != inHandle->params.Psqbar) { free(inHandle->params.Psqbar); inHandle->params.Psqbar = NULL;
	}
	if (0 != inHandle->params.Pmin_u) { free(inHandle->params.Pmin_u); inHandle->params.Pmin_u = NULL;
	}
	if (0 != inHandle->params.actmin) { free(inHandle->params.actmin); inHandle->params.actmin = NULL;
	}
	if (0 != inHandle->params.actmin_sub) { free(inHandle->params.actmin_sub); inHandle->params.actmin_sub = NULL;
	}
	if (0 != inHandle->params.lmin_flag) { free(inHandle->params.lmin_flag); inHandle->params.lmin_flag = NULL;
	}
	
	if (NULL != inHandle->params.mincat_val)
	{
		for (OsInt32 i = 0; i < N; i++)
			free(inHandle->params.mincat_val[i]);
		free(inHandle->params.mincat_val);
		inHandle->params.mincat_val = NULL;
	}
	
	if (0 != inHandle->fft_table)
	{
		spx_fft_destroy(inHandle->fft_table); inHandle->fft_table = NULL;
	}

	

	if (inHandle)
	{
		free(inHandle); inHandle = NULL;
	}
}


void martin(struct NoiseReduction *inHandle) {
	
	OsFlt32 *k_mod,*Bmin,*Bmin_sub,*YFRAME = inHandle->ps,*noise_est,eps = 2.2204e-16,Bc;
	OsInt32 N = inHandle->frame_size,i,j,sum=0;
	k_mod = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(k_mod);
	noise_est = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(noise_est);
	memcpy(noise_est, inHandle->noise, sizeof(OsFlt32)*N);
	OsFlt32 sum_P = 0, sum_YFRAME = 0;
		for (i = 0; i < N; i++)
		{
			sum_P += inHandle->params.p[i];
			sum_YFRAME += YFRAME[i];
		}

		OsFlt32 alpha_corr_t = 1 / (1 + (sum_P / sum_YFRAME - 1)*(sum_P / sum_YFRAME - 1));
		inHandle->params.alpha_corr *= 0.7;
		inHandle->params.alpha_corr = alpha_corr_t > 0.7 ? 0.3*alpha_corr_t : 0.7*0.3;
		OsFlt32 temp, Qeq_tild, Qeq_tild_sub, sum_Qeq = 0;
		for (i = 0; i < N; i++) {
			temp = inHandle->params.p[i] / (noise_est[i] + eps) - 1;
			inHandle->params.alpha[i] = (inHandle->params.alpha_max*inHandle->params.alpha_corr) / (temp*temp + 1);
		}

		Bmin = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(Bmin);
		Bmin_sub = (OsFlt32*)calloc(N, sizeof(OsFlt32)); assert(Bmin_sub);

		for (i = 0; i < N; i++) {

			inHandle->params.p[i] = inHandle->params.alpha[i] * inHandle->params.p[i] + (1 - inHandle->params.alpha[i])*YFRAME[i];

			temp = inHandle->params.alpha[i] * inHandle->params.alpha[i];
			temp = temp < inHandle->params.beta_max ? temp : inHandle->params.beta_max;
			inHandle->params.Pbar[i] = temp*inHandle->params.Pbar[i] + (1 - temp)*inHandle->params.p[i];
			inHandle->params.Psqbar[i] = temp* inHandle->params.Psqbar[i] + (1 - temp)* (inHandle->params.p[i] * inHandle->params.p[i]);
			//////////////////////////////////////////////////////////////////////////////////////////////////
			temp = fabsf(inHandle->params.Psqbar[i] - inHandle->params.Pbar[i] * inHandle->params.Pbar[i]);
			temp = temp / (2 * (noise_est[i] * noise_est[i]));
			temp = temp < 0.5 ? temp : 0.5;
			temp = 1 / (temp + eps);
			sum_Qeq += 1 / temp;
			Qeq_tild = (temp - 2 * inHandle->params.M_D) / (1 - inHandle->params.M_D);
			Qeq_tild_sub = (temp - 2 * inHandle->params.M_V) / (1 - inHandle->params.M_V);

			Bmin[i] = 1 + (inHandle->params.D - 1) * 2 / Qeq_tild;
			Bmin_sub[i] = 1 + (inHandle->params.V - 1) * 2 / Qeq_tild_sub;
		}

		temp = (1 / inHandle->params.L)*sum_Qeq;
		Bc = 1 + inHandle->params.Av*sqrtf(temp);
		OsFlt32 Qinv_bar = temp;
		for (i = 0; i < N; i++) {
			if (inHandle->params.p[i] * Bmin[i] * Bc < inHandle->params.actmin[i])
			{
				k_mod[i] = 1;
				inHandle->params.actmin_sub[i] = inHandle->params.p[i] * Bmin_sub[i] * Bc;
				inHandle->params.actmin[i] = inHandle->params.p[i] * Bmin[i] * Bc;

			}
		}


		if (inHandle->params.subwc == inHandle->params.V) {
			for (i = 0; i < N; i++)
			{
				if (k_mod[i] == 1)
					inHandle->params.lmin_flag[i] = 0;
				inHandle->params.mincat_val[i][inHandle->params.u - 1] = inHandle->params.actmin[i];
				inHandle->params.Pmin_u[i] = inHandle->params.actmin[i];//inHandle->params.mincat_val[i][0];
				for (j = 0; j < inHandle->params.iter; j++)
					if (inHandle->params.Pmin_u[i] > inHandle->params.mincat_val[i][j])
						inHandle->params.Pmin_u[i] = inHandle->params.mincat_val[i][j];



			}
			OsFlt32 noise_slope_max;
			if (Qinv_bar < 0.03)
				noise_slope_max = 8;
			else if (Qinv_bar < 0.05)
				noise_slope_max = 4;
			else if (Qinv_bar < 0.06)
				noise_slope_max = 2;
			else
				noise_slope_max = 1.2;

			for (i = 0; i < N; i++) {
				if (inHandle->params.lmin_flag[i] == 1 && (inHandle->params.actmin_sub[i] < (noise_slope_max*inHandle->params.Pmin_u[i])) && (inHandle->params.actmin_sub[i] > inHandle->params.Pmin_u[i]))
				{
					inHandle->params.Pmin_u[i] = inHandle->params.actmin_sub[i];
					for (j = 0; j < inHandle->params.Um; j++) {
						inHandle->params.mincat_val[i][j] = inHandle->params.actmin_sub[i];
					}

					inHandle->params.actmin[i] = inHandle->params.actmin_sub[i];
				}

			}
			memset(inHandle->params.lmin_flag, 0, N * sizeof(OsFlt32));
			inHandle->params.subwc = 1;
			memcpy(inHandle->params.actmin, inHandle->params.p, sizeof(OsFlt32)*N);
			memcpy(inHandle->params.actmin_sub, inHandle->params.p, sizeof(OsFlt32)*N);
			if (inHandle->params.u == inHandle->params.Um)
				inHandle->params.u = 1;
			else
				inHandle->params.u++;
		}
		else {
			if (inHandle->params.subwc >=1) {
				for (i = 0; i < N; i++) {
					if (k_mod[i] == 1)
						inHandle->params.lmin_flag[i] = 1;
					noise_est[i] = inHandle->params.actmin_sub[i] < inHandle->params.Pmin_u[i] ? inHandle->params.actmin_sub[i] : inHandle->params.Pmin_u[i];
				}
				memcpy(inHandle->params.Pmin_u, noise_est, sizeof(OsFlt32)*N);
			}
			inHandle->params.subwc++;
		}

		memcpy(inHandle->noise, noise_est, sizeof(OsFlt32)*N);
		if (k_mod) { free(k_mod); k_mod = NULL; }
		if (Bmin) { free(Bmin); Bmin = NULL; }
		if (Bmin_sub) { free(Bmin_sub); Bmin_sub = NULL; }
		//if (YFRAME) { free(YFRAME); YFRAME = NULL; }
		if (noise_est) { free(noise_est); noise_est = NULL; }
	
}
void IMCRA(struct NoiseReduction *inHandle){
	OsInt32 scount=0,i, N = inHandle->ps_size;
	OsFlt32 temp = 0.0, tilde_alpha_d, nu, q, gamma_min, zeta, s0 = 0.0f, s2 = 0.0f, smin0 = 0.0f, smin2 = 0.0f, sf0 = 0.0f, sf2 = 0.0f, I0, I1, I2;
	inHandle->cf.count++;
	//printf("%f\n", inHandle->noise[0]);
	for (i = 0; i < N; i++){
		
		if (1 == inHandle->adapt_count){
			inHandle->cf.Smin[i] = 0.0f; inHandle->cf.Smin_sw[i] = 0.0f; inHandle->cf.tilde_Smin[i] = 0.0f; inHandle->cf.tilde_Smin_sw[i] = 0.0f;
			//inHandle->noise[i] = inHandle->ps[i];
			inHandle->gammak_old[i] = (inHandle->ps[i] / inHandle->noise[i] < 1000.0f) ? inHandle->ps[i] / inHandle->noise[i] : 1000.0f;
			inHandle->ksi_old[i] = ((inHandle->gammak_old[i] - 1.0f)>0.0f) ? (1.0f - 0.90f) * (inHandle->gammak_old[i] - 1.0f) : 0.0f;
			if (0 == i) 
				temp = 0.5f*inHandle->ps[i] + 0.25*inHandle->ps[i + 1];
			else if (N - 1 == i) temp = 0.25f*inHandle->ps[i - 1] + 0.5*inHandle->ps[i];
			else temp = 0.25f*inHandle->ps[i - 1] + 0.5f*inHandle->ps[i] + 0.25*inHandle->ps[i + 1];
			inHandle->cf.S[i]=temp;
			inHandle->cf.tilde_Sf[i] = temp;
			inHandle->cf.Sf[i] = temp;
			inHandle->cf.tilde_Smin_sw[i] = temp;
			inHandle->cf.Smin_sw[i] = temp;
			inHandle->cf.tilde_Smin[i] = temp;
			inHandle->cf.Smin[i] = temp;
			inHandle->cf.tilde_S[i] = temp;
			inHandle->cf.ov_Lambda_D[i] = inHandle->ps[i];
			inHandle->cf.P[i] = 1;
		}
		if (inHandle->adapt_count <= inHandle->cf.IS){
			if (0 == i) inHandle->cf.Sf[i] = 0.5f*inHandle->ps[i] + 0.25f*inHandle->ps[i + 1];
			else if (N - 1 == i) inHandle->cf.Sf[i] = 0.25f*inHandle->ps[i-1] + 0.5f*inHandle->ps[i];
			else inHandle->cf.Sf[i] = 0.25f*inHandle->ps[i - 1] + 0.5f*inHandle->ps[i] + 0.25f*inHandle->ps[i + 1];
			inHandle->cf.S[i] = inHandle->cf.alpha_s*inHandle->cf.S[i] + (1.0f - inHandle->cf.alpha_s)*inHandle->cf.Sf[i];
			inHandle->cf.Smin[i] = (inHandle->cf.Smin[i] < inHandle->cf.S[i]) ? inHandle->cf.Smin[i] : inHandle->cf.S[i];
			inHandle->cf.Smin_sw[i] = (inHandle->cf.Smin_sw[i] < inHandle->cf.S[i]) ? inHandle->cf.Smin_sw[i] : inHandle->cf.S[i];
			inHandle->noise[i] = inHandle->cf.alpha_d*inHandle->noise[i] + (1.0f - inHandle->cf.alpha_d)*inHandle->ps[i];
			inHandle->cf.tilde_S[i] = inHandle->cf.S[i]; inHandle->cf.tilde_Sf[i] = inHandle->cf.Sf[i]; inHandle->cf.tilde_Smin[i] = inHandle->cf.Smin[i];
			inHandle->cf.P[i] = 0;

		}
		else{
			if (0 == i){
				inHandle->cf.Sf[i] = 0.5f*inHandle->ps[i] + 0.25f*inHandle->ps[i + 1];
				sf2 = 0.25f*inHandle->ps[i] + 0.5f*inHandle->ps[i + 1]+0.25f*inHandle->ps[i+2];
				inHandle->cf.S[i] = inHandle->cf.alpha_s*inHandle->cf.S[i] + (1.0f - inHandle->cf.alpha_s)*inHandle->cf.Sf[i];
				s2 = inHandle->cf.alpha_s*inHandle->cf.S[i+1] + (1.0f - inHandle->cf.alpha_s)*sf2;
				inHandle->cf.Smin[i] = (inHandle->cf.Smin[i] < inHandle->cf.S[i]) ? inHandle->cf.Smin[i] : inHandle->cf.S[i];
				smin2 = (inHandle->cf.Smin[i+1] < s2) ? inHandle->cf.Smin[i+1] : s2;
				inHandle->cf.Smin_sw[i] = (inHandle->cf.Smin_sw[i] < inHandle->cf.S[i]) ? inHandle->cf.Smin_sw[i] : inHandle->cf.S[i];
				gamma_min = inHandle->ps[i] / (inHandle->cf.Bmin*inHandle->cf.Smin[i]);
				zeta = inHandle->cf.S[i] / (inHandle->cf.Bmin*inHandle->cf.Smin[i]);
				I1 = 0.0f;
				if ((gamma_min < inHandle->cf.Gamma0) && (zeta < inHandle->cf.zeta0)) I1 = 1.0f;
				gamma_min = inHandle->ps[i+1] / (inHandle->cf.Bmin*smin2);
				zeta = s2 / (inHandle->cf.Bmin*smin2);
				I2 = 0.0f;
				if ((gamma_min < inHandle->cf.Gamma0) && (zeta < inHandle->cf.zeta0)) I2 = 1.0f;

				if ((OsInt32)(I1 + I2) != 0) inHandle->cf.tilde_Sf[i] = (0.5f*I1*inHandle->ps[i] + 0.25f*I2*inHandle->ps[i + 1]) / (0.5f*I1 + 0.25f*I2);
				else inHandle->cf.tilde_Sf[i] = inHandle->cf.tilde_S[i];



			}
			else if (N - 1 == i){
				inHandle->cf.Sf[i] = 0.25f*inHandle->ps[i-1] + 0.5f*inHandle->ps[i];
				//sf0 = inHandle->cf.Sf[i - 1];
				inHandle->cf.S[i] = inHandle->cf.alpha_s*inHandle->cf.S[i] + (1.0f - inHandle->cf.alpha_s)*inHandle->cf.Sf[i];
				s0 = inHandle->cf.S[i - 1];
				inHandle->cf.Smin[i] = (inHandle->cf.Smin[i] < inHandle->cf.S[i]) ? inHandle->cf.Smin[i] : inHandle->cf.S[i];
				smin0 = inHandle->cf.Smin[i - 1];
				inHandle->cf.Smin_sw[i] = (inHandle->cf.Smin_sw[i] < inHandle->cf.S[i]) ? inHandle->cf.Smin_sw[i] : inHandle->cf.S[i];
				gamma_min = inHandle->ps[i] / (inHandle->cf.Bmin*inHandle->cf.Smin[i]);
				zeta = inHandle->cf.S[i] / (inHandle->cf.Bmin*inHandle->cf.Smin[i]);
				I1 = 0.0f;
				if ((gamma_min < inHandle->cf.Gamma0) && (zeta < inHandle->cf.zeta0)) I1 = 1.0f;
				gamma_min = inHandle->ps[i - 1] / (inHandle->cf.Bmin*smin0);
				zeta = s0 / (inHandle->cf.Bmin*smin0);
				I0 = 0.0f;
				if ((gamma_min < inHandle->cf.Gamma0) && (zeta < inHandle->cf.zeta0)) I0 = 1.0f;

				if ((OsInt32)(I1 + I0) != 0) inHandle->cf.tilde_Sf[i] = (0.5f*I1*inHandle->ps[i] + 0.25f*I0*inHandle->ps[i - 1]) / (0.5f*I1 + 0.25f*I0);
				else inHandle->cf.tilde_Sf[i] = inHandle->cf.tilde_S[i];
			}
			else{
				inHandle->cf.Sf[i] = 0.25f*inHandle->ps[i - 1] + 0.5f*inHandle->ps[i] + 0.25f*inHandle->ps[i + 1];
				if (N-2==i)
					sf2 = 0.25f*inHandle->ps[i] + 0.5f*inHandle->ps[i + 1];
				else sf2 = 0.25f*inHandle->ps[i] + 0.5f*inHandle->ps[i + 1] + 0.25f*inHandle->ps[i + 2];
				inHandle->cf.S[i] = inHandle->cf.alpha_s*inHandle->cf.S[i] + (1.0f - inHandle->cf.alpha_s)*inHandle->cf.Sf[i];
				inHandle->cf.Smin[i] = (inHandle->cf.Smin[i] < inHandle->cf.S[i]) ? inHandle->cf.Smin[i] : inHandle->cf.S[i];

				s0 = inHandle->cf.S[i - 1];
				smin0 = inHandle->cf.Smin[i - 1];

				s2 = inHandle->cf.alpha_s*inHandle->cf.S[i + 1] + (1.0f - inHandle->cf.alpha_s)*sf2;
				smin2 = (inHandle->cf.Smin[i + 1] < s2) ? inHandle->cf.Smin[i + 1] : s2;

				inHandle->cf.Smin_sw[i] = (inHandle->cf.Smin_sw[i] < inHandle->cf.S[i]) ? inHandle->cf.Smin_sw[i] : inHandle->cf.S[i];
				gamma_min = inHandle->ps[i] / (inHandle->cf.Bmin*inHandle->cf.Smin[i]);
				zeta = inHandle->cf.S[i] / (inHandle->cf.Bmin*inHandle->cf.Smin[i]);
				I1 = 0.0f;
				if ((gamma_min < inHandle->cf.Gamma0) && (zeta < inHandle->cf.zeta0)) I1 = 1.0f;
				gamma_min = inHandle->ps[i - 1] / (inHandle->cf.Bmin*smin0);
				zeta = s0 / (inHandle->cf.Bmin*smin0);
				I0 = 0.0f;
				if ((gamma_min < inHandle->cf.Gamma0) && (zeta < inHandle->cf.zeta0)) I0 = 1.0f;
				gamma_min = inHandle->ps[i + 1] / (inHandle->cf.Bmin*smin2);
				zeta = s2 / (inHandle->cf.Bmin*smin2);
				I2 = 0.0f;
				if ((gamma_min < inHandle->cf.Gamma0) && (zeta < inHandle->cf.zeta0)) I2 = 1.0f;
				if ((OsInt32)(I1 + I0 + I2) != 0) inHandle->cf.tilde_Sf[i] = (0.5f*I1*inHandle->ps[i] + 0.25f*I0*inHandle->ps[i - 1] + 0.25f*I2*inHandle->ps[i + 1]) / (0.5f*I1 + 0.25f*I0 + 0.25f*I2);
				else inHandle->cf.tilde_Sf[i] = inHandle->cf.tilde_S[i];
			}
			inHandle->cf.tilde_S[i] = inHandle->cf.alpha_s*(inHandle->cf.tilde_S[i]) + (1.0f - inHandle->cf.alpha_s)*inHandle->cf.tilde_Sf[i];
			inHandle->cf.tilde_Smin[i] = (inHandle->cf.tilde_Smin[i] < inHandle->cf.tilde_S[i]) ? inHandle->cf.tilde_Smin[i] : inHandle->cf.tilde_S[i];
			inHandle->cf.tilde_Smin_sw[i] = (inHandle->cf.tilde_Smin_sw[i] < inHandle->cf.tilde_S[i]) ? inHandle->cf.tilde_Smin_sw[i] : inHandle->cf.tilde_S[i];

			gamma_min = inHandle->ps[i] / (inHandle->cf.Bmin*inHandle->cf.tilde_Smin[i]);
			zeta = inHandle->cf.S[i] / (inHandle->cf.Bmin*inHandle->cf.tilde_Smin[i]);
			q = 0.0f;
			if ((gamma_min <= 1)&&(zeta < inHandle->cf.zeta0)) q = 1.0f;
			else if ((1 < gamma_min) && (gamma_min < inHandle->cf.Gamma1) && (zeta < inHandle->cf.zeta0)) q = (inHandle->cf.Gamma1 - gamma_min) / (inHandle->cf.Gamma1 - 1.0f);
			nu = inHandle->gammak_old[i] * inHandle->ksi_old[i] / (1.0f + inHandle->ksi_old[i]);
			inHandle->cf.P[i] = 0.0f;
			if (1.0f>q)
				inHandle->cf.P[i] = 1.0f / (1.0f + (q / (1.0f + q))*(1.0f + inHandle->ksi_old[i])*exp(-nu));
			tilde_alpha_d = inHandle->cf.alpha_d + (1.0f - inHandle->cf.alpha_d)*inHandle->cf.P[i];
			inHandle->cf.ov_Lambda_D[i] = tilde_alpha_d*inHandle->noise[i] + (1.0f - tilde_alpha_d)*inHandle->ps[i];
			inHandle->noise[i] = inHandle->cf.ov_Lambda_D[i];

		}

	}

}
void InitNoise(struct NoiseReduction *inHandle, OsInt16 *ioPcm)
{
	OsInt32 N = inHandle->ps_size;
	OsInt32 i;
	inHandle->adapt_count++;
	FrequencyAnalysis(inHandle, ioPcm);
	if (8 >= inHandle->adapt_count)
	{
		for (i = 0; i < N; i++)
		{
			if(1 == inHandle->adapt_count)
				inHandle->noise[i] = 0.0;
			inHandle->noise[i] += inHandle->ps[i];
		}
	}
	
}

void AudioEnhancementSMPR(struct NoiseReduction *inHandle)
{
	OsInt32 i;
	OsFlt32 snru,alpha = 0.99f, delta = 0.0002f,  vu = 0.0f, P = 0.0f, temp1 = 0.0f, gain = 0.0f, temp2 = 0.0f;
	inHandle->ft[0] = 0;
	for (i = 1; i < inHandle->ps_size; i++){
		P = 0;
		inHandle->gammak[i] = (inHandle->ps[i] / inHandle->noise[i] < 1000.0f) ? inHandle->ps[i] / inHandle->noise[i] : 1000.0f;
		if (1 == inHandle->adapt_count)
			inHandle->ksi[i] = ((inHandle->gammak[i] - 1.0f)>0.0f )? (1.0f - alpha) * (inHandle->gammak[i] - 1.0f) : 0.0f;
		else
		{
			//temp1 = ((inHandle->gammak[i] - 1.0f)>0.0f) ? inHandle->gammak[i] - 1.0f : 0.0f;
			temp1 = inHandle->gammak[i] - 1.0;
			temp2 = fmax(alpha*inHandle->X2_old[i] / inHandle->noise_mu2_old[i] + (1.0f - alpha)*temp1,0.0126);
			inHandle->ksi[i] = temp2;// (temp2 > 0.0126f) ? temp2 : 0.0126f;
		}
		inHandle->noise_mu2_old[i] = inHandle->noise[i];
		gain = 1;
		if (inHandle->ksi[i] < 0.0025f) inHandle->ksi[i] = 0.0025f;
		
			//gain = 0.5;
			P = inHandle->ksi[i] / (inHandle->theta + inHandle->ksi[i]);
			gain = P + (1 - P)*0.01;

		
		gain = gain < 1.0 ? gain : 1.0;

		inHandle->gammak_old[i] = inHandle->gammak[i];
		inHandle->ksi_old[i] = inHandle->ksi[i];
		
		
		gain = sqrtf(gain);
	
		
		gain = gain > 0.0631 ? gain : 0.0631;

		inHandle->X2_old[i] = sqrtf(inHandle->ps[i]) * gain ;
		inHandle->X2_old[i] *= inHandle->X2_old[i];
		inHandle->ft[2 * i - 1] *= gain;
		inHandle->ft[2 * i] *= gain;
	}
}
void NoiseReductionProcess(struct NoiseReduction *inHandle, OsInt16 *ioPcm)
{
    OsInt32 i;
    OsInt32 N = inHandle->ps_size;
    OsFlt32 szAudio[1024] = {0};

    inHandle->adapt_count++;

    FrequencyAnalysis(inHandle,ioPcm);
	IMCRA(inHandle);
    AudioEnhancementSMPR(inHandle);
    spx_ifft(inHandle->fft_table,inHandle->ft,szAudio);
    for(i = 0; i < inHandle->frame_size; i++)
		ioPcm[i] = 1.852*(inHandle->outbuf[i] + szAudio[i])/2.0;
	for (i = 0; i < inHandle->frame_size; i++)
		inHandle->outbuf[i] = szAudio[inHandle->frame_size + i];
    if(inHandle->adapt_count > 16000) inHandle->adapt_count = 2;
}
void NoiseReductionProcess_martin(struct NoiseReduction *inHandle, OsInt16 *ioPcm)
{
	OsInt32 i;
	OsInt32 N = inHandle->ps_size;
	OsFlt32 szAudio[1024] = { 0 },sum=0;

	inHandle->adapt_count++;
	if (inHandle->adapt_count == 1)
		
			init_martin_params(inHandle, ioPcm);
	else {
		
			FrequencyAnalysis(inHandle, ioPcm);
		
			martin(inHandle);
		
	}
	
	AudioEnhancementSMPR(inHandle);
	spx_ifft(inHandle->fft_table, inHandle->ft, szAudio);
	for (i = 0; i < inHandle->frame_size; i++)
		ioPcm[i] = inHandle->outbuf[i] + szAudio[i];

	for (i = 0; i < inHandle->frame_size; i++)
		inHandle->outbuf[i] = szAudio[inHandle->frame_size + i]*inHandle->win_gain;
	if (inHandle->adapt_count > 500000) inHandle->adapt_count = 2;
}


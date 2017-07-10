/*-----------------------------------
	林明安 2016.4.27
	Time scaling based on locked phase vocoder;
*/

#include"TimeScaling.h"
#include "_kiss_fft_guts.h"
#include "kiss_fftr.h"
#include "kiss_fft.h"
#include <stdio.h>
#include <complex> 
#include "SignalBasicFunc.h"
using namespace std;
#define kiss_fft_scalar OsFlt64
#define M_PI 4*atan(1)



OsFlt64 *Time_Scaling(OsFlt64 *input, const OsInt32 Lx,const OsFlt64 time_scaling_ratio){
	complex<double> *out_dft;
	SignalBasicFunc *sbf = new SignalBasicFunc();
	kiss_fft_scalar *rin;
	kiss_fft_cpx *dft,*dft_half;
	kiss_fftr_cfg  kiss_fftr_state;
	OsFlt64 *window,                    *DFT_bin_freqs,                *output_signal,
		    *this_analysis_phase,       *principal_determination,      *partials_freq,
		    *this_synthesis_phase,      *last_analysis_phase,          *last_synthesis_phase,
			*analysis_frame,			*delta_phase,				   *temp,						
			*sp,						*phase_increment,			   *systhesis_frame,
			COLA_ratio=0.0,				ttt;
	
	OsInt32 i,j,pin = 0, pout = 0, Ly = round(time_scaling_ratio*Lx);

	OsInt16 NFFT = 512, frame_length = NFFT, *peaks, cp, *regions,
		    synthesis_frame_shift = frame_length / 4,
		    analysis_frame_shift = round(synthesis_frame_shift / time_scaling_ratio);

	/*=======================================================================
									get memery
	========================================================================*/
	dft = (kiss_fft_cpx *)calloc(NFFT, sizeof(kiss_fft_cpx));				assert(dft);
	dft_half = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*(NFFT / 2));		assert(dft_half);
	rin = (kiss_fft_scalar *)malloc(sizeof(kiss_fft_scalar)*NFFT);			assert(rin);
	analysis_frame = (OsFlt64*)malloc(sizeof(OsFlt64)*NFFT);				assert(analysis_frame);
	systhesis_frame = (OsFlt64*)malloc(sizeof(OsFlt64)*NFFT);				assert(systhesis_frame);
	window = (OsFlt64*)malloc(sizeof(OsFlt64)*NFFT);						assert(window);
	sp = (OsFlt64*)malloc(sizeof(OsFlt64)*NFFT / 2);						assert(sp);
	output_signal = (OsFlt64*)calloc(Ly,sizeof(OsFlt64));					assert(output_signal);
	this_analysis_phase = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2));		assert(this_analysis_phase);
	this_synthesis_phase = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2));	assert(this_synthesis_phase);
	last_analysis_phase = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2));		assert(last_analysis_phase);
	last_synthesis_phase = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2));	assert(last_synthesis_phase);
	principal_determination = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2)); assert(principal_determination);
	partials_freq = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2));			assert(partials_freq);
	DFT_bin_freqs = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2));			assert(DFT_bin_freqs);
	phase_increment = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2));			assert(phase_increment);
	temp = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2));					assert(temp);
	delta_phase = (OsFlt64*)malloc(sizeof(OsFlt64)*(NFFT / 2));				assert(delta_phase);
	out_dft = (complex<double> *)calloc(NFFT / 2, sizeof(complex<double>)); assert(out_dft);
	peaks = (OsInt16 *)malloc(sizeof(OsInt16));								assert(peaks);
	regions = (OsInt16*)calloc(NFFT, sizeof(OsInt16));						assert(regions);
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/



	kiss_fftr_state = kiss_fftr_alloc(NFFT, 0, 0, 0);
	memset(this_analysis_phase, 0, sizeof(OsFlt64)*NFFT / 2);
	memset(this_synthesis_phase, 0, sizeof(OsFlt64)*NFFT / 2);
	memset(last_analysis_phase, 0, sizeof(OsFlt64)*NFFT / 2);
	memset(last_synthesis_phase, 0, sizeof(OsFlt64)*NFFT / 2);
	memset(delta_phase, 0, sizeof(OsFlt64)*NFFT / 2);
	memset(phase_increment, 0, sizeof(OsFlt64)*NFFT / 2);
	
	/*=======================================================================
								creat window
	========================================================================*/
	sbf->hanning(NFFT, window);
	for (i = 0; i < NFFT; i++)
		COLA_ratio += *(window + i) * *(window + i);
	COLA_ratio /= (OsFlt64)synthesis_frame_shift;
	for (i = 0; i < NFFT; i++)
		*(window + i) /= sqrt(COLA_ratio);
	for (i = 0; i < NFFT / 2; i++)
		*(DFT_bin_freqs + i) = 2.0 * M_PI*(OsFlt64)i / (OsFlt64)NFFT;
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/




	while ((pin + frame_length < Lx) && (pout + frame_length < Ly))
	{
		memcpy(analysis_frame, input + pin, sizeof(OsFlt64)*frame_length);// 读取一帧数据

		for (i = 0; i < NFFT; i++)
			*(analysis_frame + i) *= *(window + i);//加窗

		/*=======================================================================
					shif move the data,like this: [1 2 3 4]->[3 4 1 2]
		========================================================================*/
		memcpy(temp, analysis_frame, sizeof(OsFlt64)*frame_length / 2);
		memcpy(analysis_frame, analysis_frame + (frame_length / 2), sizeof(OsFlt64)*frame_length / 2);
		memcpy(analysis_frame + (frame_length / 2), temp, sizeof(OsFlt64)*frame_length / 2);
		/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


		/*=======================================================================
							FFT:using kiss fft tool
		========================================================================*/
		if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }
		kiss_fftr_state = kiss_fftr_alloc(NFFT, 0, 0, 0);
		kiss_fftr(kiss_fftr_state, analysis_frame, dft);
		memcpy(dft_half, dft, sizeof(kiss_fft_cpx)*NFFT / 2);
		/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


		for (i = 0; i < NFFT / 2; i++)
		{
			*(this_analysis_phase + i) = atan2(dft_half[i].i, dft_half[i].r);//相位
			*(sp + i) = sqrt(dft_half[i].r*dft_half[i].r + dft_half[i].i*dft_half[i].i);//幅值
			*(delta_phase + i) = *(this_analysis_phase + i) - *(last_analysis_phase + i);//获取帧之间的相位差
			*(phase_increment + i) = *(delta_phase + i) - analysis_frame_shift* (*(DFT_bin_freqs + i));

		}

		cp = 0;
		if (peaks) { free(peaks); peaks = NULL; }
		peaks = (OsInt16 *)malloc(sizeof(OsInt16)); assert(peaks);


		/*=======================================================================
										find peaks
		========================================================================*/
		for (i = 0; i < NFFT / 2 - 2; i++)
		{
			if ((*(sp + i + 1)>*(sp + i + 2)) & (*(sp + i + 1) >= *(sp + i))){
				if (0 == cp)
				{
					*peaks = (OsInt16)i + 1;
					cp++;

				}
				else
				{

					peaks = (OsInt16*)realloc(peaks, sizeof(OsInt16)*(cp + 1));
					*(peaks + cp) = (OsInt16)i + 1;
					cp++;
				}
			}

		}
		/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


		/*=======================================================================
					creat the output data phase whit locked
		========================================================================*/
		memset(principal_determination, 0, sizeof(OsFlt64)*NFFT / 2);
		memset(partials_freq, 0, sizeof(OsFlt64)*NFFT / 2);
		for (i = 0; i < cp; i++)
		{

			ttt = fmod(*(phase_increment + *(peaks + i)) + M_PI, M_PI * 2.0);
			if (ttt < 0)
				ttt += 2 * M_PI;
			*(principal_determination + *(peaks + i)) = ttt - M_PI;

			*(partials_freq + *(peaks + i)) = *(principal_determination + *(peaks + i)) / analysis_frame_shift + *(DFT_bin_freqs + *(peaks + i));
		}

		
		for (i = 0; i < cp - 1; i++)
			*(regions + i + 1) = round(0.5*(*(peaks + i) + *(peaks + i + 1)));
		*regions = 0; *(regions + cp) = NFFT / 2 - 1;

		for (i = 0; i < cp; i++)
		{
			for (j = *(regions + i); j <= *(regions + i + 1); j++)
				*(partials_freq + j) = *(partials_freq + *(peaks + i));
			*(this_synthesis_phase + *(peaks + i)) = *(last_synthesis_phase + *(peaks + i)) + synthesis_frame_shift * (*(partials_freq + *(peaks + i)));
		}


		for (i = 0; i < cp; i++)
		{
			for (j = *(regions + i); j <= *(regions + i + 1); j++)
				*(this_synthesis_phase + j) = *(this_synthesis_phase + *(peaks + i)) + *(this_analysis_phase + j) - *(this_analysis_phase + *(peaks + i));

		}
		/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


		/*=======================================================================
				creat the complex values of output data for IFFT
		========================================================================*/
		if (out_dft) { free(out_dft); out_dft = NULL; }
		out_dft = (complex<double> *)calloc(NFFT / 2, sizeof(complex<double>)); assert(out_dft);
		for (i = 0; i < NFFT / 2; i++)
		{
			//out_dft[i].imag = this_synthesis_phase[i];
			out_dft[i] = complex<double>(0, this_synthesis_phase[i]);
			out_dft[i] = exp(out_dft[i]);
			
			//out_dft[i].imag *= *(sp + i);
			//out_dft[i].real *= *(sp + i);
			out_dft[i] *= complex<double>(*(sp + i), *(sp + i));

		}
		
		
		out_dft = (complex<double> *)realloc(out_dft, NFFT*sizeof(complex<double>)); assert(out_dft);
		for (i = NFFT / 2 + 1; i < NFFT; i++)
		{
			//out_dft[i].real = out_dft[NFFT - i].real;
			//out_dft[i].imag = out_dft[NFFT - i].imag;
			out_dft[i] = out_dft[NFFT - i];
		}

		//out_dft[NFFT / 2].real = 0;
		//out_dft[NFFT / 2].imag = 0;
		out_dft[NFFT / 2] = complex<double>(0, 0);
		for (i = 0; i < NFFT; i++)
		{
			dft[i].r = out_dft[i].real();
			dft[i].i = out_dft[i].imag();
		}
		/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


		/*=======================================================================
						IFFT:using kiss fftr to get the wav data
		========================================================================*/

		if (kiss_fftr_state) { free(kiss_fftr_state); kiss_fftr_state = NULL; }
		kiss_fftr_state = kiss_fftr_alloc(NFFT, 1, 0, 0);
	
		kiss_fftri(kiss_fftr_state, dft, rin);

		for (i = 0; i < NFFT; i++)
			*(systhesis_frame + i) = rin[i] / NFFT;
		/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


		/*=======================================================================
					shif move the data,like this: [1 2 3 4]->[3 4 1 2]
		========================================================================*/
		memcpy(temp, systhesis_frame, sizeof(OsFlt64)*frame_length / 2);
		memcpy(systhesis_frame, systhesis_frame + (frame_length / 2), sizeof(OsFlt64)*frame_length / 2);
		memcpy(systhesis_frame + (frame_length / 2), temp, sizeof(OsFlt64)*frame_length / 2);
		/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


		for (i = 0; i < NFFT; i++)
			*(systhesis_frame + i) *= *(window + i);//加窗
		for (i = pout; i < pout + frame_length; i++)
			*(output_signal + i) += *(systhesis_frame + i - pout);
		memcpy(last_synthesis_phase, this_synthesis_phase, sizeof(OsFlt64)*NFFT / 2);
		memcpy(last_analysis_phase, this_analysis_phase, sizeof(OsFlt64)*NFFT / 2);

		pin += analysis_frame_shift;
		pout += synthesis_frame_shift;

	}//end while


	/*======================================================================
									free memery
	========================================================================*/
	if (dft != NULL) {
		free(dft); dft = NULL;
	}
	if (dft_half != NULL) {
		free(dft_half); dft_half = NULL;
	}
	if (rin != NULL) {
		free(rin); rin = NULL;
	}

	if (analysis_frame != NULL) {
		free(analysis_frame); analysis_frame = NULL;
	}
	if (systhesis_frame != NULL) {
		free(systhesis_frame); systhesis_frame = NULL;
	}

	if (window != NULL) {
		free(window); window = NULL;
	}

	if (sp != NULL) {
		free(sp); sp = NULL;
	}

	if (this_analysis_phase != NULL) {
		free(this_analysis_phase); this_analysis_phase = NULL;
	}

	if (this_synthesis_phase != NULL) {
		free(this_synthesis_phase); this_synthesis_phase = NULL;
	}

	if (last_analysis_phase != NULL) {
		free(last_analysis_phase); last_analysis_phase = NULL;
	}

	if (last_synthesis_phase != NULL)
	{
		free(last_synthesis_phase); last_synthesis_phase = NULL;
	}

	if (principal_determination != NULL) {
		free(principal_determination); principal_determination = NULL;
	}
	if (partials_freq != NULL) {
		free(partials_freq); partials_freq = NULL;
	}
	
	if (DFT_bin_freqs != NULL) {
		free(DFT_bin_freqs); DFT_bin_freqs = NULL;
	}

	if (phase_increment != NULL) {
		free(phase_increment); phase_increment = NULL;
	}

	if (temp != NULL) {
		free(temp); temp = NULL;
	}

	if (delta_phase != NULL) {
		free(delta_phase); delta_phase = NULL;
	}

	if (out_dft != NULL) {
		free(out_dft); out_dft = NULL;
	}

	if (kiss_fftr_state != NULL) {
		free(kiss_fftr_state); kiss_fftr_state = NULL;
	}

	if (peaks != NULL) {
		free(peaks); peaks = NULL;
	}

	if (regions != NULL) {
		free(regions); regions = NULL;
	}
	delete sbf;
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	

	return output_signal;
}

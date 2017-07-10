#ifndef FFTWRAP_H
#define FFTWRAP_H

#include "Platform.h"

#ifdef  __cplusplus
extern "C" {
#endif  __cplusplus

/* Compute tables for an FFT */
void *spx_fft_init(int size);
/* Destroy tables for an FFT */
void spx_fft_destroy(void *table);
/* Forward (real to half-complex) transform */
void spx_fft(void *table, OsFlt32 *in, OsFlt32 *out);
/* Backward (half-complex to real) transform */
void spx_ifft(void *table, OsFlt32 *in, OsFlt32 *out);

#ifdef  __cplusplus
}
#endif  __cplusplus

#endif
#include "Platform.h"
#include "avsmallft.h"

void *spx_fft_init(int size)
{
    struct drft_lookup *table;
    table = (struct drft_lookup*)calloc(sizeof(struct drft_lookup),1);
    spx_drft_init((struct drft_lookup *)table,size);
    return (void*)table;
}

void spx_fft_destroy(void *table)
{
    spx_drft_clear((struct drft_lookup*)table);
    free(table);
}

void spx_fft(void *table,float *in,float *out)
{
    if (in == out)
    {
        int i;
        float scale = 1.0f/((struct drft_lookup *)table)->n;
        fprintf(stderr,"warning: %s\n","FFT should not be done in-place");
        for (i = 0; i < ((struct drft_lookup *)table)->n; i++)
            out[i] = scale*in[i];
    }
    else
    {
        int i;
        float scale = 1.0f/((struct drft_lookup *)table)->n;
        for (i = 0; i < ((struct drft_lookup *)table)->n; i++)
            out[i] = scale*in[i];
    }
    spx_drft_forward((struct drft_lookup*)table,out);
}

void spx_ifft(void *table,float *in,float *out)
{
    if (in == out)
    {
        fprintf(stderr,"warning: %s\n","FFT should not be done in-place");
    }
    else
    {
        int i;
        for (i=0;i<((struct drft_lookup *)table)->n;i++)
            out[i] = in[i];
    }
    spx_drft_backward((struct drft_lookup *)table,out);
}
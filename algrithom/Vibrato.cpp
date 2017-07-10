
#include"Vibrato.h"
#include<stdlib.h>
#include<math.h>
#include<stdlib.h>
#define PI 3.14159265358979323846 
using namespace std;
void VibratoEffect(const short *x, const int Lx, const int fs, short *y, const double modfreq, const double width) {
	 
	 double *Delayline,MODFREQ = modfreq / fs,MOD,ZEIGER;
	int frac,k,DELAY = width*fs, WIDTH = width*fs,  L = 2 + DELAY + WIDTH * 2;
	Delayline = (double*)calloc(L, sizeof(double));
	for (int i = 0; i < Lx - 1; i++) {
		MOD = sin(MODFREQ * 2 * PI*(i+1));
		ZEIGER = 1 + DELAY + WIDTH*MOD;
		k = floor(ZEIGER);
		frac = ZEIGER - k;
		for (int j = L - 2; j >= 0; j--)
			*(Delayline + j + 1) = *(Delayline + j);
		*Delayline = *(x + i);
		*(y + i) = *(Delayline + k)*frac + *(Delayline + k - 1)*(1 - frac);
	}
	if (Delayline) { free(Delayline); Delayline = NULL; }
}
#include"compressor.h"




Compressor::Compressor() {}
Compressor::~Compressor() {}
void gaincomputerFF(OsFlt64 *lg, OsInt32 L,OsFlt64 Threshold, OsFlt64 ratio, OsFlt64 Knee,OsFlt64 *cv) {
	OsFlt64 slope = 1 / ratio - 1,overshoot=0,w2=0,a=0,rect;
	OsInt32 i;
	bool inTransition = 0;
	for (i = 0; i < L; i++)
	{
		overshoot = *(lg + i) - Threshold;
		if (Knee > 0)
		{
			w2 = Knee / 2;
			a = w2 / (Knee*Knee);
			inTransition = ((overshoot > -w2)&(overshoot < w2));
			rect = overshoot > 0 ? inTransition*(a*(overshoot + w2)*(overshoot + w2)) + (1 - inTransition)*overshoot : inTransition*(a*(overshoot + w2)*(overshoot + w2));
		}
		else rect = overshoot > 0 ? overshoot : 0;

		*(cv + i) = -rect*slope;
	}
}

void PDbranching(OsFlt64 *x,OsInt32 Lx, OsFlt32 fs, OsFlt64 tauAttack, OsFlt64 tauRelease, OsFlt64 *y) {
	OsFlt64 alphaAtt = 0, alphaRel = 0, state = 0;
	OsInt32 i;
	alphaAtt = tauAttack > 0 ? exp(-1 / (tauAttack*fs)) : 0;
	alphaRel = tauRelease > 0 ? exp(-1 / (tauRelease*fs)) : 0;
	for (i = 0; i < Lx; i++)
	{
		state = *(x + i) > state ? alphaAtt*state + (1 - alphaAtt)* *(x + i) : alphaRel*state;
		*(y + i) = -state;
	}
}
void PDbranchingsmooth(OsFlt64 *x, OsInt32 Lx, OsFlt32 fs, OsFlt64 tauAttack, OsFlt64 tauRelease, OsFlt64 *y){

	OsFlt64 alphaAtt = 0, alphaRel = 0, state = 0;
	OsInt32 i;
	alphaAtt = tauAttack > 0 ? exp(-1 / (tauAttack*fs)) : 0;
	alphaRel = tauRelease > 0 ? exp(-1 / (tauRelease*fs)) : 0;
	for (i = 0; i < Lx; i++)
	{
		state = *(x + i) > state ? alphaAtt*state + (1 - alphaAtt)* *(x + i) : alphaRel*state+ (1 - alphaRel)* *(x + i);
		*(y + i) = -state;
	}
}
void PDdecoupled(OsFlt64 *x, OsInt32 Lx, OsFlt32 fs, OsFlt64 tauAttack, OsFlt64 tauRelease, OsFlt64 *y){
	OsFlt64 alphaAtt = 0, alphaRel = 0, state = 0;
	OsInt32 i;
	alphaAtt = tauAttack > 0 ? exp(-1 / (tauAttack*fs)) : 0;
	alphaRel = tauRelease > 0 ? exp(-1 / (tauRelease*fs)) : 0;
	for (i = 0; i < Lx; i++)
	{
		state = *(x + i) > state ?  *(x + i) : alphaRel*state;
		*(y + i) = state;
		
		
	}
	state = 0;
	for (i = 0; i < Lx; i++)
	{
		state = alphaAtt * state + (1 - alphaAtt)* *(y + i);

		*(y + i) = -state;
	}
}
void PDdecoupledsmooth(OsFlt64 *x, OsInt32 Lx, OsFlt32 fs, OsFlt64 tauAttack, OsFlt64 tauRelease, OsFlt64 *y) {
	OsFlt64 alphaAtt = 0, alphaRel = 0, state = 0;
	OsInt32 i;
	alphaAtt = tauAttack > 0 ? exp(-1 / (tauAttack*fs)) : 0;
	alphaRel = tauRelease > 0 ? exp(-1 / (tauRelease*fs)) : 0;
	for (i = 0; i < Lx; i++)
	{
		state = *(x + i) > state ? *(x + i) : alphaRel*state+(1-alphaRel)* *(x+i);
		*(y + i) = state;

		
	}
	state = 0;
	for (i = 0; i < Lx; i++)
	{
		state = alphaAtt * state + (1 - alphaAtt)* *(y + i);

		*(y + i) = -state;
	}
}

void Compressor::ffcompressor(OsInt16 type, OsFlt64 *x,OsInt32 Lx, OsInt32 fs, OsFlt64 Threshold, OsFlt64 ratio, OsFlt64 tauAttack, OsFlt64 tauRelease, OsFlt64 knee, OsFlt64 *y) {
	OsFlt64 *lg,*cv,*cvp;
	OsInt32 i;

	lg = (OsFlt64*)malloc(sizeof(OsFlt64)*Lx); assert(lg);
	cv = (OsFlt64*)malloc(sizeof(OsFlt64)*Lx); assert(cv);
	cvp = (OsFlt64*)malloc(sizeof(OsFlt64)*Lx); assert(cvp);

	for (i = 0; i < Lx; i++)
		*(lg + i) = fabs(*(x + i)) > 1e-6 ? 20 * log10(fabs(*(x + i))) : 20 * log10(1e-6);

	gaincomputerFF(lg, Lx,Threshold, ratio, knee, cv);
	if (0 == type)//B
		PDbranching(cv, Lx, fs, tauAttack,tauRelease,cvp);
	else if(1==type)//BS
		PDbranchingsmooth(cv, Lx, fs, tauAttack, tauRelease, cvp);
	else if(2==type)//D
		PDdecoupled(cv, Lx, fs, tauAttack, tauRelease, cvp);
	else//DS
		PDdecoupledsmooth(cv, Lx, fs, tauAttack, tauRelease, cvp);

	for (i = 0; i < Lx; i++) 
		*(y + i) = *(x + i) * pow(10, *(cvp + i) / 20);
	if (lg != NULL) {
		free(lg); lg = NULL;
	}
	if (cv != NULL) {
		free(cv); cv = NULL;
	}
	if (cvp != NULL) {
		free(cvp); cvp = NULL;
	}
}
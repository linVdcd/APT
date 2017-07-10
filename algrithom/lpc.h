/**************************************

林明安 2016.6.3
LPC分析与合成

*************************************/

#ifndef RESAMPLE_H
#define RESAMPLE_H

#include"Platform.h"
#include"SignalBasicFunc.h"

class Lpc {
public:
	SignalBasicFunc *sbf = new SignalBasicFunc();
	void lpcfit(OsFlt64 *x, OsInt32 Lx, OsInt16 p, OsInt16 h, OsFlt64 **a, OsFlt64 **g, OsFlt64 **e, OsInt16 *nframe);
	void lpcsynth(OsFlt64 *a, OsFlt64 *g, OsFlt64 *e, OsInt16 p, OsInt16 h, OsInt16 nframe, OsFlt64 **y,OsInt32 &Ly);
	Lpc();
	~Lpc();
};


#endif



#ifndef  COMPRESSOR_H
#define COMPRESSOR_H

#include "Platform.h"
class Compressor{
	
public:
	Compressor();
	~Compressor();
void ffcompressor(OsInt16 type, OsFlt64 *x, OsInt32 Lx,OsInt32 fs, OsFlt64 Threshold, OsFlt64 ratio, OsFlt64 tauAttack, OsFlt64 tauRelease, OsFlt64 knee, OsFlt64 *y);
};


#endif // ! COMPRESSOR_H

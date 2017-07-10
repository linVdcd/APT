#pragma once
#ifndef  REVERB_H
#define REVERB_H
#include "Platform.h"
void convReverb(OsFlt64 *x, OsInt32 Lx, OsFlt64 *h, OsInt32 Lh, OsFlt64 *y);


#endif // ! REVERB_H
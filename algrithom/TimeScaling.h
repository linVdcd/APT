/*-----------------------------------
	ÁÖÃ÷°² 2016.4.27
	Time scaling based on locked phase vocoder;
*/

#ifndef TIMESCALING_H
#define TIMESCALING_H

#include "TimeScaling.h"
#include "Platform.h"
/*
#ifdef __cplusplus
extern "C"{
#endif
*/
	OsFlt64 *Time_Scaling(OsFlt64 *input, const OsInt32 Lx, const OsFlt64 time_scaling_ratio);
/*
#ifdef __cplusplus
}
#endif
*/

#endif
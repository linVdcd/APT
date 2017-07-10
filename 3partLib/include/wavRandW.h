/*-----------------------------------
	ÁÖÃ÷°² 2016.4.25
	A wav file reader and writer interface.
*/


#ifndef WAVRANDW_H
#define WAVRANDW_H

#include"Platform.h"
class wavRandW
{
public:
	OsInt16 *indata;
	OsInt16 *outdata;
	OsUInt32 indatasize;
	OsUInt32 outdatasize;
	OsInt32 infs;
	OsInt32 outfs;
	OsInt16 inbit;
	OsInt16 outbit;
	wavRandW();
	~wavRandW();
	void wavRead(const char *filename);
	void wavWrite(const char *filename);
	
};

#endif
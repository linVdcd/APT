#ifndef WAVREADER_H
#define WAVREADER_H

#include "Platform.h"

class WavReader
{
public:
	bool isFLLR;
    WavReader();
    ~WavReader();

    bool Open(const OsInt8 *inFileName);
    void Close();
    OsInt32 Read(OsInt16 *ioData,OsInt32 inLen);

    OsInt32 GetDataSize() { return m_nLength; }
    OsUInt32 GetSampleRate() { return m_nSampleRate; }
    OsUInt16 GetChannels() { return m_nChannels; }
    OsUInt16 GetBitsPerSample() { return m_nBitsPerSample; }
	OsInt32     m_nLength;
private:
    FILE        *m_pFile;
    
    OsUInt16    m_nFormatTag;   // format type
    OsUInt16    m_nChannels;
    OsUInt32    m_nSampleRate;
    OsUInt32    m_nAvgBytesPerSec;
    OsUInt16    m_nBlockAlign;
    OsUInt16    m_nBitsPerSample;
	
	
};

#endif
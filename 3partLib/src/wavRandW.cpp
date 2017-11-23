/*-----------------------------------
	������ 2016.4.25
	A wav file reader and writer interface.
*/

#include"wavRandW.h"
#include"WavReader.h"
#include"WavWriter.h"

wavRandW::wavRandW() :indata(0), 
					  indatasize(0), 
					  infs(0), 
					  inbit(0),
					  outdata(0),
					  outdatasize(0),
					  outfs(0),
					  outbit(0){}

void wavRandW::wavRead(const char *filename){
	OsInt16 *temp;
	WavReader *pReader = new WavReader();
	assert(0 != pReader);
	bool bRet = pReader->Open(filename);
	assert(true == bRet);
	indatasize = pReader->GetDataSize();
	infs = pReader->GetSampleRate();
	inbit = pReader->GetBitsPerSample();
	if (indatasize > 0){
		indata = (OsInt16*)malloc(sizeof(OsInt16)*(indatasize/2));
		assert(NULL != indata);
		OsInt32 ReadLen = pReader->Read(indata, (indatasize/2));
		
		assert(ReadLen*2 == indatasize);
	}
	if (pReader->isFLLR == true)
	{
		indatasize -= 2466*2;
		
		temp = (OsInt16*)calloc(indatasize/2, sizeof(OsInt16)); assert(temp);
		memcpy(temp, indata + 2465, (indatasize/2)*sizeof(OsInt16));
		if (indata != NULL) free(indata);
		indata = (OsInt16*)calloc(indatasize / 2,  sizeof(OsInt16)); assert(indata);
		memcpy(indata, temp, (indatasize / 2) *sizeof(OsInt16));
		if (NULL != temp) { free(temp); temp = NULL; }
	}
	pReader->Close();
	delete pReader;
}
void wavRandW::wavWrite(const char *filename){
	WavWriter *pWriter = new WavWriter(); 
	assert(0 != pWriter);
	bool bRet = pWriter->Open(filename, outfs, outbit, 1); if (bRet == NULL) { printf("%s ��ʧ�ܣ��ļ���ռ�ã���ر����ԡ�\n", filename); return; }
	//assert(true == bRet);
	if (outdatasize > 0)
		pWriter->Write(outdata, outdatasize);
	pWriter->Close();
	delete pWriter;
}
wavRandW::~wavRandW(){
	if (NULL != indata) {
		free(indata);
		indata=NULL;
	}
	/*if (NULL != outdata)
	{

		free(outdata);
		outdata=NULL;
	}*/
		
}

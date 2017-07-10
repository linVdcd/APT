
#ifndef APTOOL_H
#define APTOOL_H
/*
	����Ƶ���������У���wav�ļ���Ϊ��������������Ƶ�����������õ���wav�����ʽҪ�����£�
		����������ͨ��
		�����ʣ�16KHz
		ÿ���������λ����16bit
	ͬʱ�����ڴ�������Ϊ�������Ƶ���������������Ƶ����ҲӦ���������õ�wav��ʽ��ͬ��16KHz��16bit��
*/


#include<string>


class aptool
{

public:
#
	void ProcesseByOrder(std::string cfg_file, std::string in_wav_file, std::string out_wav_file, std::string bundlePath);
	
	void DeNoiseMartin(const char *in, const char* out);
	short int* DeNoiseMartin(  short int *input,const int Lx,const double theta);
	void DeNoise(const char *in, const char* out);
	short int* DeNoise(const  short int *input, const int Lx);//���ؽ�������Ƶ���ݣ����ݲ���ʹ��ʱ��Ҫ�ͷ��ڴ档

	//void EQ(char *in, char *out, OsFlt64 *ac);
	void EQ2(const char *in, const char *out, const double *ac, const double *Q);
	short int* EQ2(const short int *input, const int Lx, const double *ac, const double *Q);//���ؽ�������Ƶ���ݣ����ݲ���ʹ��ʱ��Ҫ�ͷ��ڴ档

	void PandFshift(const char *in, const char *out, const double rate);
	short int* PandFshift(const short int *input, const int Lx,const double rate, int &Ly);
	void Resample(const char *in, const char *out, const double rate);

	void TimeScaling(const char *in, const char *out, const  double rate);
	short int* TimeScaling(const short int* input, const int Lx, const double rate, int &Ly);
	//void limiter1(const char *in, const  char *out, const  OsFlt64 tresh);

	void Limiter2(const char *in, const char *out, const  double CT, const double CS, const double ET, const double ES);
	short int* Limiter2(const short int* input, const int Lx, const  double CT, const double CS, const double ET, const double ES);

	void Compressing(const char *in, const  char *out, const  double type, const double Threshold, const double ratio, const double knee);
	short int* Compressing(const short int* input, const int Lx, const  short int type, const double Threshold, const double ratio, const double knee);

	void Delay(const char *in, const  char *out, const double ratio);
	short int* Delay(const short int* input, const int Lx, const double ratio, int &Ly);

	void ReverbConv(const char *in1, const char *in2, const char *out,  const double wet, const double dry);
	short int* ReverbConv(const short int* input1, const int Lx, const short int* input2, const int Lh, const double wet, const double dry, int &Ly);

	short int* Robotization(const short *input, const int Lx, short type, int fs, double RMfre, int outwin);
	short int* Whisper(const short *input, const int Lx, int &Ly,int p);

	short int* VibratoProcesse(const short *input, const int Lx, const int fs,const double modfreq, const double width);
	
	short int* FormantChange(const short *input, const int Lx, const double warp_rate);

	short int* Chorus(const short*input, const int Lx, const short num_chorus, int &Ly,const double *rato, const double *gain, const double *shift_time);
	aptool();
	~aptool();
};
#endif

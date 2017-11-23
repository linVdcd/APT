
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
#include <vector>
using namespace std;
class aptool
{
    int fs;
    int bit;
    const float norm;
public:
#
	void ProcesseByOrder(std::string cfg_file, std::string in_wav_file, std::string out_wav_file, std::string bundlePath);
	
	void DeNoiseMartin(const char *in, const char* out);
	short int* DeNoiseMartin(  short int *input,const int Lx,const double theta);
    void DeNoiseMartin(  const vector<short> &input,const double theta,vector<short>&output);


	void DeNoise(const char *in, const char* out);
	short int* DeNoise(const  short int *input, const int Lx);//���ؽ�������Ƶ���ݣ����ݲ���ʹ��ʱ��Ҫ�ͷ��ڴ档

	//void EQ(char *in, char *out, OsFlt64 *ac);
	void EQ2(const char *in, const char *out, const double *ac, const double *Q);
	short int* EQ2(const short int *input, const int Lx, const double *ac, const double *Q);//���ؽ�������Ƶ���ݣ����ݲ���ʹ��ʱ��Ҫ�ͷ��ڴ档
    void EQ2(const vector<short> &input,const double *ac, const double *Q,vector<short> &output);

	void PandFshift(const char *in, const char *out, const double rate);
	short int* PandFshift(const short int *input, const int Lx,const double rate, int &Ly);
    void PandFshift(const vector<short> &input,const double rate, vector<short>& output);
	void Resample(const char *in, const char *out, const double rate);

	void TimeScaling(const char *in, const char *out, const  double rate);
	short int* TimeScaling(const short int* input, const int Lx, const double rate, int &Ly);
    void TimeScaling(const vector<short> &input,const double rate, vector<short>& output);
	//void limiter1(const char *in, const  char *out, const  OsFlt64 tresh);

	void Limiter2(const char *in, const char *out, const  double CT, const double CS, const double ET, const double ES);
	short int* Limiter2(const short int* input, const int Lx, const  double CT, const double CS, const double ET, const double ES);
    void Limiter2(const vector<short> &input,  const  double CT, const double CS, const double ET, const double ES, vector<short>& output);

	void Compressing(const char *in, const  char *out, const  double type, const double Threshold, const double ratio, const double knee);
	short int* Compressing(const short int* input, const int Lx, const  short int type, const double Threshold, const double ratio, const double knee);
	void Compressing(const vector<short>& input,const  short int type, const double Threshold, const double ratio, const double knee,vector<short>&output);

	void Delay(const char *in, const  char *out, const double ratio);
	short int* Delay(const short int* input, const int Lx, const double ratio, int &Ly);
    void Delay(const vector<short> &input, const double ratio, vector<short>& output);

	void ReverbConv(const char *in1, const char *in2, const char *out,  const double wet, const double dry);
	short int* ReverbConv(const short int* input1, const int Lx, const short int* input2, const int Lh, const double wet, const double dry, int &Ly);
	void ReverbConv(const vector<short> &input1, const vector<short>& input2, const double wet, const double dry,vector<short>&output);

	short int* Robotization(const short *input, const int Lx, short type, int fs, double RMfre, int outwin);
    void Robotization(const vector<short> &input, short type, double RMfre, int outwin,vector<short> &output);


	short int* Whisper(const short *input, const int Lx, int &Ly,int p);
    void Whisper(const vector<short> &input,int p,vector<short>&output);

	short int* VibratoProcesse(const short *input, const int Lx, const int fs,const double modfreq, const double width);
	void VivratoProcesse(const vector<short>& input,const double modfreq,const double width,vector<short>&output);

	short int* FormantChange(const short *input, const int Lx, const double warp_rate);
    void FormantChange(const vector<short>&input,const double warp_rate,vector<short>&output);

	short int* Chorus(const short*input, const int Lx, const short num_chorus, int &Ly,const double *rato, const double *gain, const double *shift_time);



	void Gain(const vector<short> &input, vector<short> &output,float gain);
	aptool(const int fs,const int bit);
	~aptool();
};
#endif

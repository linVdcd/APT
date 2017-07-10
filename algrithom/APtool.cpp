// 
//

#include "APtool.h"
#include "ns.h"
#include "wavRandW.h"
#include"Vibrato.h"
#include "TimeScaling.h"
#include "Resample.h"
#include "compressor.h"
#include "Reverb.h"
#include "robot.h"
#include "SignalBasicFunc.h"
#include "FormantShift.h"


#include<ctype.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
using namespace std;

aptool::aptool() {}
aptool::~aptool() {

}

std::string GetFileCString(std::string filepath)
{
	std::ifstream is(filepath.c_str());
	std::string filebuffer = "",temp;
	/*if (is) {
		// get length of file:
		is.seekg(0, is.end);
		long long length = is.tellg();
		is.seekg(0, is.beg);

		char * buffer = new char[length];

		std::cout << "Reading " << length << " characters... ";
		// read data as a block:
		is.read(buffer, length);

		if (is)
			std::cout << "all characters read successfully.";
		else
			std::cout << "error: only " << is.gcount() << " could be read";
		is.close();

		// ...buffer contains the entire file...
		filebuffer = buffer;
		delete[] buffer;
		cout << filebuffer;
	}*/
	while (getline(is, temp)) {
		temp += "\n";
		filebuffer += temp;
	}
	return filebuffer;
}
std::string &trim(std::string& s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>( ::isspace))).base(), s.end());
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>( ::isspace))));
	return s;

}
inline void split(const string& src, vector<string>& res, const string& pattern, size_t maxsplit = string::npos) {
	res.clear();
	size_t start = 0;
	size_t end = 0;
	string sub;
	while (start < src.size()) {
		end = src.find_first_of(pattern, start);
		if (string::npos == end || res.size() >= maxsplit) {
			sub = src.substr(start);
			trim(sub);
			res.push_back(sub);
			return;
		}
		sub = src.substr(start, end - start);
		trim(sub);
		res.push_back(sub);
		start = end + 1;
	}
	return;
}
inline vector<string> split(const string& src, const string& pattern, size_t maxsplit = string::npos) {
	vector<string> res;
	split(src, res, pattern, maxsplit);
	return res;
}
void processAll(std::string filebuffer, std::string in_wav_file, std::string out_wav_file, std::string bundlePath)
{
	//test here
	//���ڲ�ͬ�Ĺ��ܺ������ص����ݳ��Ȳ�һ���������ֲ������ȼ������������ĳ��ȣ���˲Űѷ��ص������ڹ��ܺ����ж�̬���١�
	//��ʹ�����һ�����⣬��ĳ�����ܺ����������Ϊ�������ܺ���������ʱ������ֱ��ʹ����һ�����ܺ���������yout����Ϊ��һ�����ܺ�����ĳһ�����롣����y=Denoise��x��;y=EQ(y);�����ᷢ���ڴ�й¶��
	//�ڸú��������ǲ��ã�yout=Denoise(y);free(y);yout=y; ��һ��ʽ��

	double ac[8] = { 0,0,0,-10,0,0,0,0 }, ac2[10] = { 0,0,0,0,10 * 0.6,0,0,0,0,0 }, Q[10] = { 2,2,2,2,2,2,2,2,2,2 };
	int Lx, Ly = 0, Lh,fs=16000,nbit=16;
	short *y,*yout;
	wavRandW *wav = new wavRandW();
	aptool *ap = new aptool();
	wav->outbit = nbit;
	wav->outfs = fs;
	wav->wavRead((in_wav_file).c_str());
	Lx = wav->indatasize / 2;
	y = (short*)calloc(Lx, sizeof(short));
	memcpy(y, wav->indata, sizeof(short)*Lx);//���ʹ���ڴ渴�ƣ���Ȼ���ͷ�wav�ڴ�ʱ���ܻ�����ظ��ͷ�ͬһ���ڴ�������
	//    wav->outdatasize = Lx;
	//    wav->outdata = y;
	//

	std::string line;
	vector<string> buf;
	std::vector<std::string>buffervec = split(filebuffer, "\n");
	std::vector<std::pair<std::string, std::vector<std::string> > > params(buffervec.size());
	for (int lineno = 0; lineno<buffervec.size(); lineno++)
	{
		line = buffervec[lineno];
		split(line, buf, " ");
		std::string filetername = "";
		if (lineno == 0)
		{
			filetername = "denoise";
		}
		if (lineno == 1)
		{
			filetername = "comp";
		}
		if (lineno == 2)
		{
			filetername = "limiter";
		}
		if (lineno == 3)
		{
			filetername = "EQ";
		}
		if (lineno == 4)
		{
			filetername = "time";
		}
		if (lineno == 5)
		{
			filetername = "pitch";
		}
		if (lineno == 6)
		{
			filetername = "delay";
		}
		if (lineno == 7)
		{
			filetername = "reverb";
		}
		if (lineno == 8)
		{
			filetername = "robot";
		}
		if (lineno == 9)
		{
			filetername = "whisper";
		}
		if (lineno == 10)
		{
			filetername = "vibrato";
		}
		if (lineno == 11)
		{
			filetername = "formantshift";
		}
		if (lineno == 12)
		{
			filetername = "chorus";
		}
		int order = atof(buf[1].c_str());
		std::pair<std::string, std::vector<std::string> > subparams = make_pair(filetername, buf);
		if(subparams.second[0]=="True")
			params[order-1] = subparams;
	}
	for (int i = 0; i<params.size(); i++)
	{
		std::pair<std::string, std::vector<std::string> > subparams = params[i];
		std::vector<std::string> buf = subparams.second;
		if (buf.size() == 0)
		{
			std::cout << "error column " << i << std::endl;
			continue;
		}
		if (buf[0] == "False")
			continue;
		if (subparams.first == "denoise")
		{
			
				yout = ap->DeNoiseMartin(y, Lx,atof(buf[2].c_str()));
				if (y) { free(y); y = NULL; }
				y = yout;
			
		}
		if (subparams.first == "comp")
		{
			yout = ap->Compressing(y, Lx, atof(buf[2].c_str()), atof(buf[3].c_str()), atof(buf[4].c_str()), atof(buf[5].c_str()));
			if (y) { free(y); y = NULL; }
			y = yout;
		}
		if (subparams.first == "limiter")
		{
			yout = ap->Limiter2(y, Lx, atof(buf[2].c_str()), 1.3,-50.0,-0.001);
			if (y) { free(y); y = NULL; }
			y = yout;

		}
		if (subparams.first == "EQ")
		{
			double ac[10];
			double eq[10];
			for (int j = 0; j<10; j++)
				ac[j] = 0.6*atof(buf[j + 2].c_str());
			for (int j = 0; j<10; j++)
				eq[j] = atof(buf[j + 12].c_str());
			yout = ap->EQ2(y, Lx, ac, eq);
			if (y) { free(y); y = NULL; }
			y = yout;
		}
		if (subparams.first == "time")
		{
			yout = ap->TimeScaling(y, Lx, atof(buf[2].c_str()), Ly);
			Lx = Ly;
			if (y) { free(y); y = NULL; }
			y = yout;
		}
		if (subparams.first == "pitch")
		{
			yout = ap->PandFshift(y, Lx, atof(buf[2].c_str()), Ly);
			Lx = Ly;
			if (y) { 
				free(y); y = NULL; 
			}
			y = yout;

		}
		if (subparams.first == "delay")
		{
			yout = ap->Delay(y, Lx, atof(buf[2].c_str()), Ly);
			Lx = Ly;
			if (y) { free(y); y = NULL; }
			y = yout;
		}
		if (subparams.first == "reverb")
		{
			wavRandW *wav2 = new wavRandW();
			std::cout << (bundlePath + buf[2]) << std::endl;
			wav2->wavRead((bundlePath  + buf[2] + ".wav").c_str());
			int Lh = wav2->indatasize / 2;
			short* y2 = wav2->indata;
			yout = ap->ReverbConv(y, Lx, y2, Lh, atof(buf[3].c_str()), 1, Ly);
			Lx = Ly;
			if (y) { free(y); y = NULL; }
			y = yout;
			delete wav2;
		}

		if (subparams.first == "robot")
		{

			yout = ap->Robotization(y, Lx, atoi(buf[2].c_str()),fs,atof(buf[3].c_str()),atoi(buf[4].c_str()));
			
			if (y) { free(y); y = NULL; }
			y = yout;
		}
		if (subparams.first == "whisper")
		{
			yout = ap->Whisper(y, Lx, Ly,atof(buf[2].c_str()));
			Lx = Ly;
			if (y) { free(y); y = NULL; }
			y = yout;
		}
		if (subparams.first == "vibrato")
		{
			yout = ap->VibratoProcesse(y, Lx,fs, atof(buf[2].c_str()), atof(buf[3].c_str()));

			if (y) { free(y); y = NULL; }
			y = yout;
		}
		if (subparams.first == "formantshift")
		{
			yout = ap->FormantChange(y, Lx, atof(buf[2].c_str()));

			if (y) { free(y); y = NULL; }
			y = yout;
		}
		if (subparams.first == "chorus")
		{
			double *rato, *gain, *shift_time;
			short num_chorus = (short)atoi(buf[2].c_str()),size = buf.size();

			rato = (double*)calloc(num_chorus, sizeof(double)); assert(rato);
			gain = (double*)calloc(num_chorus, sizeof(double)); assert(gain);
			shift_time = (double*)calloc(num_chorus, sizeof(double)); assert(shift_time);
			for (int i = 0; i < num_chorus; i++)
			{
				rato[i] += 1.0;
				gain[i] += 1.0;
				
				if((i + 3)<size)
					rato[i] = atof(buf[i +3].c_str());
				if ((i + 3 + num_chorus)<size)
					gain[i] = atof(buf[i + 3 + num_chorus].c_str());
				if ((i + i +3 + 2 * num_chorus)<size)
					shift_time[i] = atof(buf[i + 3 + 2 * num_chorus].c_str());
			}



		


			yout = ap->Chorus(y, Lx,num_chorus,Ly,rato,gain,shift_time);
			Lx = Ly;
			if (y) { free(y); y = NULL; }
			y = yout;
		}
	}
	
	wav->outdatasize = Lx;
	wav->outdata = y;//wav->outdata,û���Լ������ڴ棬�������û����WavRandW�������������������ͷ��ڴ�������ǵ��ͷ�y��
	wav->wavWrite((out_wav_file).c_str());
	
	
	if (y) { free(y); y = NULL; }//�����ٴ��ͷ�yout����Ϊyout��yָ��ͬһ���ڴ��ַ
	//if (yout) { free(yout); yout = NULL; }
	delete wav;
	delete ap;
}


void aptool::ProcesseByOrder(std::string cfg_file, std::string in_wav_file, std::string out_wav_file, std::string bundlePath){
	std::string cfg = GetFileCString(cfg_file);
	
		processAll(cfg, in_wav_file, out_wav_file, bundlePath);
}

void aptool::Resample(const char *in, const char *out, const double rate){
	double  *py, *x, ratio,ymax =0;
	int i, Lx, Lty, Lpy, p, q;
	short int *y1;
	wavRandW *wav = new wavRandW();
	p = 100; q = rate * 100;
	ratio = (double)q / (double)p;
	wav->wavRead(in);//��ȡ������Ƶ
	Lx = wav->indatasize / 2;
	x = (double*)malloc(sizeof(double)*Lx); assert(x);

	for (i = 0; i < Lx; i++)
		*(x + i) = *(wav->indata + i) / 32768.0;

	Lty = Lx;

	py = resample(x, Lty, p, q);
	Lpy = Lty*p / q;
	y1 = (short int *)calloc(Lpy,sizeof(short int)); assert(y1);
	for (i = 0; i < Lpy; i++)
		if (ymax < fabs(py[i])) ymax = fabs(py[i]);
	for(i =0;i<Lpy;i++)
		*(y1 + i) = py[i]/(ymax+0.3) * 32768.0;
	wav->outbit = wav->inbit;
	wav->outfs = 16000;
	wav->outdatasize = Lpy;
	wav->outdata = y1;
	wav->wavWrite(out);//��Ƶ���
	if (x != NULL) {
		free(x); x = NULL;
	}

	if (py != NULL) {
		free(py); py = NULL;
	}

	if (y1 != NULL) { free(y1); y1 = NULL; }
}

void aptool::DeNoiseMartin(const char *in, const char* out)
{//����
	double ftmp = 0;
	int nLen = 512, Lx, i;
	short int *y;
	wavRandW *wav = new wavRandW();

	wav->wavRead(in);//��ȡ������Ƶ



	Lx = wav->indatasize / 2;

	for (i = 0; i < Lx; i++)
	{
		ftmp += *(wav->indata + i);
	}
	ftmp /= Lx;

	for (i = 0; i < Lx; i++)
		wav->indata[i] -= ftmp;
	NoiseReduction *pInst = NoiseReductionCreate_martin(wav->infs, nLen,5.0);
	y = (short int*)calloc(Lx, sizeof(short int)); assert(y);
	memset(y, 0, sizeof(short int)*Lx);
	assert(0 != pInst);
	// printf("��ʼ����...\n");
	int nReadLen = nLen;
	int count = 0, Nframe = Lx / nLen;
	pInst->adapt_count = 0;
	for (count = 0; count<Nframe; count++)
	{
		short int szPcm[512] = { 0 };
		memcpy(szPcm, wav->indata + count*nLen, sizeof(short int)*nLen);

		NoiseReductionProcess_martin(pInst, szPcm);
		memcpy(y + count*nLen, szPcm, sizeof(short int)*nLen);

	}
	// printf("������ɣ�\n");
	wav->outbit = wav->inbit;
	wav->outfs = wav->infs;
	wav->outdatasize = Lx;
	wav->outdata = y;
	wav->wavWrite(out);//��Ƶ�����ds.wav�ļ�

	if (y != NULL) {
		free(y); y = NULL;
	}
	NoiseReductionDestroy(pInst);

}
short int* aptool::DeNoiseMartin(short int *input, const int Lx,double theta) {
	double ftmp = 0,ymax=0,*y_f;
	int nLen = 512, fs = 16000, i;
	short int *y;
	short int szPcm[512] = { 0 };
	for (i = 0; i < Lx; i++)
	{
		ftmp += *(input + i);
	}
	ftmp /= Lx;

	for (i = 0; i < Lx; i++)
		input[i] -= ftmp;


	NoiseReduction *pInst = NoiseReductionCreate_martin(fs, nLen,theta);assert(0 != pInst);
	y = (short int*)calloc(Lx, sizeof(short int)); assert(y);
	memset(y, 0, sizeof(short int)*Lx);
	
	// printf("��ʼ����...\n");
	int nReadLen = nLen;
	int count = 0, Nframe = Lx / nLen;
	pInst->adapt_count = 0;
	for (count = 0; count<Nframe; count++)
	{
		
		memcpy(szPcm, input + count*nLen, sizeof(short int)*nLen);

		NoiseReductionProcess_martin(pInst, szPcm);
		memcpy(y + count*nLen, szPcm, sizeof(short int)*nLen);

	}
	 //printf("������ɣ�\n");
	y_f = (double*)calloc(Lx, sizeof(double)); assert(y_f);
	for (i = 0; i < Lx; i++)
	{
		y_f[i] = y[i] / 32768.0;
		if (ymax < fabs(y_f[i])) ymax = fabs(y_f[i]);
	}
		for (i = 0; i < Lx; i++)
		{
			y_f[i] /= (ymax+0.3);
			y[i] = y_f[i] * 32768.0;
		}
		if (y_f != NULL) {
			free(y_f); y_f = NULL;
		}
	NoiseReductionDestroy(pInst);
	pInst = NULL;
	return y;
}

void aptool::DeNoise(const char *in, const char* out)
{//����
	double ymax = 0;
	int nLen = 512, Lx,i;
	short int *y;
	wavRandW *wav = new wavRandW();

	wav->wavRead(in);//��ȡ������Ƶ
	Lx = wav->indatasize / 2;
	NoiseReduction *pInst = NoiseReductionCreate(wav->infs, nLen);
	y = (short int*)malloc(sizeof(short int)*Lx); assert(y);
	memset(y, 0, sizeof(short int)*Lx);
	assert(0 != pInst);
	// printf("��ʼ����...\n");
	int nReadLen = nLen;
	int count = 0, Nframe = Lx / nLen;
	for (count = 0; count<Nframe; count++)
	{

		short int szPcm1[512] = { 0 };
		memcpy(szPcm1, wav->indata + count*nLen, sizeof(short int)*nLen);
		InitNoise(pInst, szPcm1);

		if (count == 3)
			break;
	}


	for (int i = 0; i < nLen; i++) pInst->noise[i] = pInst->noise[i] / pInst->adapt_count;
	pInst->adapt_count = 0;
	for (count = 0; count<Nframe; count++)
	{
		short int szPcm[512] = { 0 };
		memcpy(szPcm, wav->indata + count*nLen, sizeof(short int)*nLen);

		NoiseReductionProcess(pInst, szPcm);
		memcpy(y + count*nLen, szPcm, sizeof(short int)*nLen);

	}
	// printf("������ɣ�\n");

	//for (i = 0; i < Lx; i++)
	//	if (ymax < fabs((double)y[i])) ymax = fabs((double)y[i]);

	//for (i = 0; i < Lx; i++)
	//	y[i] = (short)((double)y[i] / ymax*32768.0);
	wav->outbit = wav->inbit;
	wav->outfs = wav->infs;
	wav->outdatasize = Lx;
	wav->outdata = y;
	wav->wavWrite(out);//��Ƶ�����ds.wav�ļ�



	if (y != NULL) {
		free(y); y = NULL;
	}
	NoiseReductionDestroy(pInst);

}

short int* aptool::DeNoise(const short int *input, const int Lx)
{//����
	double ymax = 0,*y_f;
	int nLen = 512, i,fs = 16000;
	short int *y;
	NoiseReduction *pInst = NoiseReductionCreate(fs, nLen);
	y = (short int*)malloc(sizeof(short int)*Lx); assert(y);
	memset(y, 0, sizeof(short int)*Lx);
	assert(0 != pInst);
	// printf("��ʼ����...\n");
	int nReadLen = nLen;
	int count = 0, Nframe = Lx / nLen;
	for (count = 0; count<Nframe; count++)
	{

		short int szPcm1[512] = { 0 };
		memcpy(szPcm1, input + count*nLen, sizeof(short int)*nLen);
		InitNoise(pInst, szPcm1);

		if (count == 3)
			break;
	}


	for (int i = 0; i < nLen; i++) pInst->noise[i] = pInst->noise[i] / pInst->adapt_count;
	pInst->adapt_count = 0;
	for (count = 0; count<Nframe; count++)
	{
		short int szPcm[512] = { 0 };
		memcpy(szPcm, input + count*nLen, sizeof(short int)*nLen);

		NoiseReductionProcess(pInst, szPcm);
		memcpy(y + count*nLen, szPcm, sizeof(short int)*nLen);

	}

	y_f = (double*)calloc(Lx, sizeof(double)); assert(y_f);
	for (i = 0; i < Lx; i++)
	{
		y_f[i] = y[i] / 32768.0;
		if (ymax < fabs(y_f[i])) ymax = fabs(y_f[i]);
	}
	for (i = 0; i < Lx; i++)
	{
		y_f[i] /= (ymax + 0.2);
		y[i] = y_f[i] * 32768.0;
	}

	if (y_f != NULL) {
		free(y_f); y_f = NULL;
	}
	NoiseReductionDestroy(pInst);
	return y;

}



/*void aptool_wav::EQ(char *in, char *out, double *ac) {//EQ1
double **y_i, *y_f, *x, ymax = 1,
num1[7] = { 7.94063568068471e-15,	4.76438140841083e-14,	1.19109535210271e-13,	1.58812713613694e-13,	1.19109535210271e-13,	4.76438140841083e-14,	7.94063568068471e-15 },
den1[7] = { 1, -5.98647065367621,	14.9327478525466, -19.8662813422798,	14.8670634872237, -5.93392108923234,	0.986861745418678 },

num2[7] = { 1.76373641827193e-07,0, -5.29120925481580e-07,0,5.29120925481580e-07,0,-1.76373641827193e-07 },
den2[7] = { 1, -5.98455187979476,14.9242381389593, -19.8514198044397,14.8543499757873, -5.92863369261480,0.986017262169176 },

num3[9] = { 5.00638663803490e-08,0,-2.00255465521396e-07,0,3.00383198282094e-07,0,-2.00255465521396e-07,0,5.00638663803490e-08 },
den3[9] = { 1,-7.94938545495511,27.6567415558141, -55.0029823695991,	68.3923167694364, -54.4457828340931,	27.0992428432550, -7.71024356377983,	0.960093053938176 },

num4[9] = { 3.82377005825039e-07,	0, -1.52950802330016e-06,	0,	2.29426203495024e-06,	0, -1.52950802330016e-06,	0,	3.82377005825039e-07 },
den4[9] = { 1, -7.88690486943960,	27.2588231590015 ,-53.9241348658998,	66.7808169602917, -53.0166744081464,	26.3491300621918, -7.49543517017545,	0.934379142404103 },

num5[13] = { 1.00190062540361e-08,	0 ,-6.01140375242167e-08,	0,	1.50285093810542e-07,	0, -2.00380125080722e-07,	0,	1.50285093810542e-07,	0 ,-6.01140375242167e-08,	0,	1.00190062540361e-08 },
den5[13] = { 1, -11.5824201132465,	61.7577057551803, -200.443706405300	,441.039573632639, -693.060418854463	,797.548610532901, -677.195388620375,	421.079798440889 ,-186.993009116282,	56.2954779455226 ,-10.3165780950715,	0.870354902363389 },

num6[13] = { 3.91695463509753e-07,	0 ,-2.35017278105852e-06,	0,	5.87543195264629e-06,	0, -7.83390927019506e-06,	0,	5.87543195264629e-06,	0, -2.35017278105852e-06,	0,	3.91695463509753e-07 },
den6[13] = { 1, -10.6692801848678,	53.1561868413191 ,-163.399734933070,	344.971289537811 ,-526.780641022431	,596.478913560745, -504.582169564345	,316.513274339917, -143.607658085204	,44.7522721234145 ,-8.60517052840149,	0.772735282982836 },

num7[13] = { 2.62948642398930e-05,	0, -0.000157769185439358,	0,	0.000394422963598395,	0, -0.000525897284797859,	0,	0.000394422963598395,	0 ,-0.000157769185439358,	0,	2.62948642398930e-05 },
den7[13] = { 1, -7.53462764648966,	28.8700400474943, -72.6892439151129,	132.831452493847 ,-184.609807409799	,199.606463737528 ,-169.009656749793	,111.324218240273, -55.7641234131726	,20.2735635510999 ,-4.84457357815643	,0.589562158727239 },

num8[6] = { 0.0686270096397650 ,-0.343135048198825	,0.686270096397650 ,-0.686270096397650	,0.343135048198825 ,-0.0686270096397650 },
den8[6] = { 1 ,-0.256064783806578	,1.03483876806962	,0.218513286552870	,0.282858195421874	,0.159184152272719 };



int i, Lx, k;
short int *y1, NumBand = 8, bandord[8] = { 7,7,9,9,13,13,13,6 };
wavRandW *wav = new wavRandW();

wav->wavRead(in);//��ȡ������Ƶ
Lx = wav->indatasize / 2;
x = (double*)malloc(sizeof(double)*Lx); assert(x);
y_f = (double*)malloc(sizeof(double)*Lx); assert(y_f);
memset(y_f, 0, sizeof(double)*Lx);
for (i = 0; i < Lx; i++)
*(x + i) = *(wav->indata + i) / 32768.0;
y_i = (double**)malloc(sizeof(double)*NumBand);

for (i = 0; i < NumBand; i++)
y_i[i] = (double*)malloc(sizeof(double)*Lx);
assert(y_i);
for (i = 0; i < bandord[0]; i++)
*(num1 + i) *= pow(10, ac[0] / 20);
for (i = 0; i < bandord[1]; i++)
*(num2 + i) *= pow(10, ac[1] / 20);
for (i = 0; i < bandord[2]; i++)
*(num3 + i) *= pow(10, ac[2] / 20);
for (i = 0; i < bandord[3]; i++)
*(num4 + i) *= pow(10, ac[3] / 20);
for (i = 0; i < bandord[4]; i++)
*(num5 + i) *= pow(10, ac[4] / 20);
for (i = 0; i < bandord[5]; i++)
*(num6 + i) *= pow(10, ac[5] / 20);
for (i = 0; i < bandord[6]; i++)
*(num7 + i) *= pow(10, ac[6] / 20);
for (i = 0; i < bandord[7]; i++)
*(num8 + i) *= pow(10, ac[7] / 20);
filter(bandord[0] - 1, den1, num1, Lx - 1, x, y_i[0]);
filter(bandord[1] - 1, den2, num2, Lx - 1, x, y_i[1]);
filter(bandord[2] - 1, den3, num3, Lx - 1, x, y_i[2]);
filter(bandord[3] - 1, den4, num4, Lx - 1, x, y_i[3]);
filter(bandord[4] - 1, den5, num5, Lx - 1, x, y_i[4]);
filter(bandord[5] - 1, den6, num6, Lx - 1, x, y_i[5]);
filter(bandord[6] - 1, den7, num7, Lx - 1, x, y_i[6]);
filter(bandord[7] - 1, den8, num8, Lx - 1, x, y_i[7]);

for (i = 0; i < Lx; i++) {
for (k = 0; k < NumBand - 1; k += 2)
*(y_f + i) += y_i[k][i] * pow(10, ac[k]) + y_i[k + 1][i] * pow(10, ac[k + 1]);
if (ymax < abs(*(y_f + i))) ymax = abs(*(y_f + i));

}
y1 = (short int *)malloc(sizeof(short int)*Lx); assert(y1);

for (i = 0; i < Lx; i++) {
*(y_f + i) /= ymax;
*(y1 + i) = (short int)(*(y_f + i)*32768.0);
}


wav->outbit = wav->inbit;
wav->outfs = wav->infs;
wav->outdatasize = Lx;
wav->outdata = y1;
wav->wavWrite(out);//��Ƶ���
free(x);
free(y_f);
for (i = 0; i < NumBand; i++)
free(y_i[i]);
free(y1);

}
*/

void aptool::EQ2(const char *in, const char *out, const double *ac, const double *Q) {//EQ2:Ч���Ϻ�
    SignalBasicFunc *sbf = new SignalBasicFunc();
	double *y_f, *x, ymax = 1, A = 0, alpha = 0, *temp,
		a[3] = { 0,0,0 },
		b[3] = { 0,0,0 },
		z[3] = { 0, 0,0 },
		wi[10] = { 0.0121736715326604,	0.0247400421470196,	0.0490873852123405,	0.0981747704246810,	0.196349540849362,	0.392699081698724,	0.785398163397448,	1.57079632679490,	3.14159265358979,	6.28318530717959 },
		fi[10] = { 31,63,125,250,500,1000,2000,4000,8000,160000 };


	int i, Lx, k, Ly;
	short int *y1, NumBand = 10;
	wavRandW *wav = new wavRandW();

	wav->wavRead(in);//��ȡ������Ƶ
	Lx = wav->indatasize / 2;
	x = (double*)malloc(sizeof(double)*Lx); assert(x);
	y_f = (double*)malloc(sizeof(double)*Lx); assert(y_f);
	temp = (double*)malloc(sizeof(double)*Lx); assert(temp);
	memset(y_f, 0, sizeof(double)*Lx);
	for (i = 0; i < NumBand; i++)
		wi[i] = 2 * M_PI*fi[i] / wav->infs;
	for (i = 0; i < Lx; i++)
		*(x + i) = *(wav->indata + i) / 32768.0;

	memcpy(temp, x, sizeof(double)*Lx);
	for (i = 0; i < NumBand; i++) {
		A = pow(10, ac[i] / 20);
		alpha = sin(wi[i] / (2 * Q[i]));
		b[0] = 1 + alpha*A;
		b[1] = -2 * cos(wi[i]);
		b[2] = 1 - alpha*A;
		a[0] = 1 + alpha / A;
		a[1] = -2 * cos(wi[i]);
		a[2] = 1 - alpha / A;
		//z[0] = 0; z[1] = 0; z[2] = 0;
		for (k = 1; k < 3; k++)
		{
			a[k] /= a[0];
			b[k] /= a[0];
		}
		b[0] /= a[0];
		a[0] = 1;
		sbf->filter(2, a, b, Lx - 1, temp, y_f);
		memcpy(temp, y_f, sizeof(double)*Lx);
        delete sbf;
	}


	y1 = (short int *)malloc(sizeof(short int)*Lx); assert(y1);
	for (i = 0; i < Lx; i++)
		if (ymax < fabs(*(y_f + i))) ymax = fabs(*(y_f + i));
	for (i = 0; i < Lx; i++) {
		*(y_f + i) /= (ymax);
		*(y1 + i) = (short int)(*(y_f + i)*32768.0);
	}

	//printf("%f", pow(10, 1));
	wav->outbit = wav->inbit;
	wav->outfs = wav->infs;
	wav->outdatasize = Lx;
	wav->outdata = y1;
	wav->wavWrite(out);//��Ƶ�����ds.wav�ļ�
	if (x != NULL) {
		free(x); x = NULL;
	}

	if (y_f != NULL) { free(y_f); y_f = NULL; }

	if (temp != NULL) {
		free(temp); temp = NULL;
	}

	if (y1 != NULL) {
		free(y1); y1 = NULL;
	}

}

short int* aptool::EQ2(const short int *input, const int Lx, const double *ac, const double *Q) {
    SignalBasicFunc *sbf = new SignalBasicFunc();
	double *y_f, *x, ymax = 1, A = 0, alpha = 0, *temp,
		a[3] = { 0,0,0 },
		b[3] = { 0,0,0 },
		z[3] = { 0, 0,0 },
		wi[10] = { 0.0121736715326604,	0.0247400421470196,	0.0490873852123405,	0.0981747704246810,	0.196349540849362,	0.392699081698724,	0.785398163397448,	1.57079632679490,	3.14159265358979,	6.28318530717959 },
		fi[10] = { 31,63,125,250,500,1000,2000,4000,8000,160000 };


	int i, k, Ly, fs = 16000;
	short int *y1, NumBand = 10;




	x = (double*)malloc(sizeof(double)*Lx); assert(x);
	y_f = (double*)malloc(sizeof(double)*Lx); assert(y_f);
	temp = (double*)malloc(sizeof(double)*Lx); assert(temp);
	memset(y_f, 0, sizeof(double)*Lx);
	for (i = 0; i < NumBand; i++)
		wi[i] = 2 * M_PI*fi[i] / fs;
	for (i = 0; i < Lx; i++)
		*(x + i) = *(input + i) / 32768.0;

	memcpy(temp, x, sizeof(double)*Lx);
	for (i = 0; i < NumBand; i++) {
		A = pow(10, ac[i] / 20);
		alpha = sin(wi[i] / (2 * Q[i]));
		b[0] = 1 + alpha*A;
		b[1] = -2 * cos(wi[i]);
		b[2] = 1 - alpha*A;
		a[0] = 1 + alpha / A;
		a[1] = -2 * cos(wi[i]);
		a[2] = 1 - alpha / A;
		//z[0] = 0; z[1] = 0; z[2] = 0;
		for (k = 1; k < 3; k++)
		{
			a[k] /= a[0];
			b[k] /= a[0];
		}
		b[0] /= a[0];
		a[0] = 1;
		sbf->filter(2, a, b, Lx - 1, temp, y_f);
		memcpy(temp, y_f, sizeof(double)*Lx);
        delete sbf;
	}


	y1 = (short int *)malloc(sizeof(short int)*Lx); assert(y1);
	for (i = 0; i < Lx; i++)
		if (ymax < fabs(*(y_f + i))) ymax = fabs(*(y_f + i));
	for (i = 0; i < Lx; i++) {
		*(y_f + i) /= (ymax);
		*(y1 + i) = (short int)(*(y_f + i)*32768.0);
	}

	//printf("%f", pow(10, 1));

	if (x != NULL) {
		free(x); x = NULL;
	}

	if (y_f != NULL) { free(y_f); y_f = NULL; }

	if (temp != NULL) {
		free(temp); temp = NULL;
	}

	return y1;
}



void aptool::PandFshift(const char *in, const char *out, const  double rate) {//���ߺ͹����ͬʱ�ı�

	double *ty, *py, *x, ratio,ymax =0;
	int i, Lx, Lty, Lpy, p, q;
	short int *y1;
	wavRandW *wav = new wavRandW();
	p = 100; q = rate * 100;
	ratio = (double)q / (double)p;
	wav->wavRead(in);//��ȡ������Ƶ
	Lx = wav->indatasize / 2;
	x = (double*)malloc(sizeof(double)*Lx); assert(x);

	for (i = 0; i < Lx; i++)
		*(x + i) = *(wav->indata + i) / 32768.0;

	Lty = ceil(Lx *ratio);
	ty = Time_Scaling(x, Lx, ratio);
	py = resample(ty, Lty, p, q);
	Lpy = Lty*p / q;
	y1 = (short int *)calloc(Lpy,sizeof(short int)); assert(y1);
	for (i = 0; i < Lpy; i++)
		if (ymax < fabs(py[i])) ymax = fabs(py[i]);
	for(i =0;i<Lpy;i++)
		*(y1 + i) = py[i]/(ymax+0.3) * 32768.0;
	wav->outbit = wav->inbit;
	wav->outfs = wav->infs;
	wav->outdatasize = Lx;
	wav->outdata = y1;
	wav->wavWrite(out);//��Ƶ���
	if (x != NULL) {
		free(x); x = NULL;
	}
	if (ty != NULL) {
		free(ty); ty = NULL;
	}
	if (py != NULL) {
		free(py); py = NULL;
	}

	if (y1 != NULL) { free(y1); y1 = NULL; }

}
short int* aptool::PandFshift(const short int *input, const int Lx, const double rate, int &Ly) {
	double *ty, *py, *x, ratio,ymax=0;
	int i, Lty, Lpy, p, q;
	short int *y1;

	p = 100; q = rate * 100;
	ratio = (double)q / (double)p;


	x = (double*)malloc(sizeof(double)*Lx); assert(x);

	for (i = 0; i < Lx; i++)
		*(x + i) = *(input + i) / 32768.0;

	Lty = ceil(Lx *ratio);
	ty = Time_Scaling(x, Lx, ratio);
	py = resample(ty, Lty, p, q);
	Lpy = Lty*p / q;
	y1 = (short int *)malloc(sizeof(short int)*Lpy); assert(y1);
	for (i = 0; i < Lpy; i++)
		if (ymax < fabs(py[i])) ymax = fabs(py[i]);
	for(i=0;i<Lpy;i++)
		*(y1 + i) = round(py[i]/(ymax+0.3) * 32768.0);

	if (x != NULL) {
		free(x); x = NULL;
	}
	if (ty != NULL) {
		free(ty); ty = NULL;
	}
	if (py != NULL) {
		free(py); py = NULL;
	}
	Ly = Lpy;
	return y1;

}

void aptool::TimeScaling(const char *in, const char *out, const double rate) {//�޸�ʱ��

	double *ty, *x, ratio;
	int i, Lx, Lty, p, q;
	short int *y1;
	wavRandW *wav = new wavRandW();
	p = 100; q = rate * 100;
	ratio = (double)q / (double)p;
	wav->wavRead(in);//��ȡ������Ƶ
	Lx = wav->indatasize / 2;
	x = (double*)malloc(sizeof(double)*Lx); assert(x);

	for (i = 0; i < Lx; i++)
		*(x + i) = *(wav->indata + i) / 32768.0;

	Lty = ceil(Lx *ratio);
	ty = Time_Scaling(x, Lx, ratio);


	y1 = (short int *)malloc(sizeof(short int)*Lty); assert(y1);
	for (i = 0; i < Lty; i++)
		*(y1 + i) = round(ty[i] * 32768.0);
	wav->outbit = wav->inbit;
	wav->outfs = wav->infs;
	wav->outdatasize = Lty;
	wav->outdata = y1;
	wav->wavWrite(out);//��Ƶ���
	if (x != NULL) {
		free(x); x = NULL;
	}
	if (ty != NULL) {
		free(ty); ty = NULL;
	}

	if (y1 != NULL) { free(y1); y1 = NULL; }

}

short int* aptool::TimeScaling(const short int* input, const int Lx, const double rate, int &Ly) {
	double *ty, *x, ratio;
	int i, Lty, p, q;
	short int *y1;

	p = 100; q = rate * 100;
	ratio = (double)q / (double)p;

	x = (double*)calloc(Lx,sizeof(double)); assert(x);

	for (i = 0; i < Lx; i++)
		*(x + i) = *(input + i) / 32768.0;

	Lty = ceil(Lx *ratio);
	
	
	ty = Time_Scaling(x, Lx, ratio);
	

	y1 = (short int *)calloc(Lty,sizeof(short int)); assert(y1);
	double ymax = 0;
	for (i = 0; i < Lty; i++)
		if (ymax < fabs(ty[i])) ymax = fabs(ty[i]);
	for(i=0;i<Lty;i++)
		*(y1 + i) = ty[i]/(ymax+0.05) * 32768.0;

	if (x != NULL) {
		free(x); x = NULL;
	}
	if (ty != NULL) {
		free(ty); ty = NULL;
	}
	Ly = Lty;
	return y1;


}
/*void aptool::limiter1(const char *in, const char *out, const double tredB) {//��������ֻ�Ǽ򵥵Ľ��߹�ĳһ�ֱ�ʱ����
double x, xd_last = 0, xd, slope = 1, rt = 0.01, at = 0.4,f=0,a,tresh =pow(10,tredB/20);
int Lx,i;
short int *y;
wavRandW *wav = new wavRandW();

wav->wavRead(in);//��ȡ������Ƶ
Lx = wav->indatasize / 2;
//x = (double*)malloc(sizeof(double)*Lx); assert(x);
y = (short int*)malloc(sizeof(short int)*Lx); assert(y);
*y = 0;
for (i = 1; i < Lx; i++)
{
x = *(wav->indata + i) / 32768.0;
a = abs(x) - xd_last;
if (a < 0) a = 0;
xd = xd_last * (1 - rt) + at*a;
xd_last = xd;
if (xd > tresh)
f = pow(10, (-slope*(log10(xd) - log10(tresh))));
else
f = 1;
*(y+i)=round(x*f*32768.0);
}
wav->outbit = wav->inbit;
wav->outfs = wav->infs;
wav->outdatasize = Lx;
wav->outdata = y;
wav->wavWrite(out);//��Ƶ���
if (y != NULL) {
free(y); y = NULL;
}
}*/

void aptool::Limiter2(const char *in, const  char *out, const  double CT, const  double CS, const double ET, const  double ES) {
	//������2
	double x, X, at = 0.003, rt = 0.003, delay = 150, xrms = 0, g = 1, G = 0, f = 0, coeff = 0;
	int Lx, i;
	short int *y;
	wavRandW *wav = new wavRandW();

	wav->wavRead(in);//��ȡ������Ƶ
	Lx = wav->indatasize / 2;
	//x = (double*)malloc(sizeof(double)*Lx); assert(x);
	y = (short int*)malloc(sizeof(short int)*Lx); assert(y);
	for (i = 0; i < Lx; i++) {
		x = *(wav->indata + i) / 32768.0;
		X = log10(fabs(x));
		G = CS*(CT - X) > ES*(ET - X) ? ES*(ET - X) : CS*(CT - X);
		if (G > 0) G = 0;
		f = pow(10, G / 20);
		coeff = f < g ? at : rt;
		g = (1 - coeff)*g + coeff*(f);
		*(y + i) = round(x*g*32768.0);
	}
	wav->outbit = wav->inbit;
	wav->outfs = wav->infs;
	wav->outdatasize = Lx;
	wav->outdata = y;
	wav->wavWrite(out);//��Ƶ���
	if (y != NULL) {
		free(y); y = NULL;
	}
}

short int* aptool::Limiter2(const short int* input, const int Lx, const  double CT, const double CS, const double ET, const double ES) {
	double x, X, at = 0.003, rt = 0.003, delay = 150, xrms = 0, g = 1, G = 0, f = 0, coeff = 0;
	int  i;
	short int *y;




	//x = (double*)malloc(sizeof(double)*Lx); assert(x);
	y = (short int*)malloc(sizeof(short int)*Lx); assert(y);
	for (i = 0; i < Lx; i++) {
		x = *(input + i) / 32768.0;
		X = log10(fabs(x));
		G = CS*(CT - X) > ES*(ET - X) ? ES*(ET - X) : CS*(CT - X);
		if (G > 0) G = 0;
		f = pow(10, G / 20);
		coeff = f < g ? at : rt;
		g = (1 - coeff)*g + coeff*(f);
		*(y + i) = round(x*g*32768.0);
	}


	return y;
}

void aptool::Compressing(const char *in, const  char *out, const  double type, const double Threshold, const double ratio, const double knee) {//ѹ����
	double *x, *y_f;
	int Lx, i;
	short int *y;
	wavRandW *wav = new wavRandW();

	Compressor *comp = new Compressor();

	wav->wavRead(in);//��ȡ������Ƶ
	Lx = wav->indatasize / 2;

	x = (double*)malloc(sizeof(double)*Lx); assert(x);
	y_f = (double*)malloc(sizeof(double)*Lx); assert(y_f);
	y = (short int *)malloc(sizeof(short int)*Lx); assert(y);

	for (i = 0; i < Lx; i++)
		*(x + i) = *(wav->indata + i) / 32768.0;
	comp->ffcompressor(type, x, Lx, wav->infs, Threshold, ratio, 0.001, 0.04, knee, y_f);

	for (i = 0; i < Lx; i++)
		*(y + i) = round(*(y_f + i)*32768.0);


	wav->outbit = wav->inbit;
	wav->outfs = wav->infs;
	wav->outdatasize = Lx;
	wav->outdata = y;
	wav->wavWrite(out);//��Ƶ���
	if (y != NULL) {
		free(y); y = NULL;
	}
	if (x != NULL) {
		free(x); x = NULL;
	}
	if (y_f != NULL) {
		free(y_f); y_f = NULL;
	}
	delete comp;
}

short int* aptool::Compressing(const short int* input, const int Lx, const  short int type, const double Threshold, const double ratio, const double knee) {
	double *x, *y_f;
	int  i, fs = 16000;
	short int *y;
	Compressor *comp = new Compressor();
	x = (double*)malloc(sizeof(double)*Lx); assert(x);
	y_f = (double*)malloc(sizeof(double)*Lx); assert(y_f);
	y = (short int *)malloc(sizeof(short int)*Lx); assert(y);

	for (i = 0; i < Lx; i++)
		*(x + i) = *(input + i) / 32768.0;
	comp->ffcompressor(type, x, Lx, fs, Threshold, ratio, 0.001, 0.04, knee, y_f);

	for (i = 0; i < Lx; i++)
		*(y + i) = round(*(y_f + i)*32768.0);


	if (x != NULL) {
		free(x); x = NULL;
	}
	if (y_f != NULL) {
		free(y_f); y_f = NULL;
	}
	delete comp;
	return y;
}

void aptool::Delay(const char *in, const char *out, const double ratio) {
	double *x, *y_f, amp[11] = { 0 };
	int samp_delay = 2000, zero_padding = 10 * samp_delay, Lx, Ly, i;
	short int *y, no_delays = 10, j;
	wavRandW *wav = new wavRandW();
	wav->wavRead(in);//��ȡ������Ƶ
	Lx = wav->indatasize / 2;
	Ly = Lx + zero_padding + 1;
	x = (double*)malloc(sizeof(double)*Ly); assert(x);
	y_f = (double*)malloc(sizeof(double)*Ly); assert(y_f);
	y = (short int *)malloc(sizeof(short int)*Ly); assert(y);
	memset(x, 0, sizeof(double)*Ly);
	for (i = 0; i < Lx; i++)
	{
		*(x + i + 1) = *(wav->indata + i) / 32768.0;
		*(y_f + i + 1) = *(x + i + 1);
	}
	for (i = 1; i <= no_delays; i++) {
		amp[i] = 1 - i*ratio;
	}

	for (i = samp_delay; i < Ly; i++) {
		for (j = 1; j <= no_delays; j++) {
			if (i >(j*samp_delay))
				if (amp[j] > 0)
					*(y_f + i) += amp[j] * *(x + (i - j*samp_delay));
		}
	}
	for (i = 0; i < Ly; i++)
		*(y + i) = round(*(y_f + i)*32768.0);


	wav->outbit = wav->inbit;
	wav->outfs = wav->infs;
	wav->outdatasize = Ly;
	wav->outdata = y;
	wav->wavWrite(out);//��Ƶ���
	if (y != NULL) {
		free(y); y = NULL;
	}
	if (x != NULL) {
		free(x); x = NULL;
	}
	if (y_f != NULL) {
		free(y_f); y_f = NULL;
	}
}

short int* aptool::Delay(const short int* input, const int Lx, const double ratio, int &Ly) {
	double *x, *y_f, amp[11] = { 0 },ymax=0;
	int samp_delay = 2000, zero_padding = 10 * samp_delay, i;
	short int *y, no_delays = 10, j;

	Ly = Lx + zero_padding + 1;
	x = (double*)calloc(Ly,sizeof(double)); assert(x);
	y_f = (double*)calloc(Ly,sizeof(double)); assert(y_f);
	y = (short int *)calloc(Ly,sizeof(short int)); assert(y);
	memset(x, 0, sizeof(double)*Ly);
	for (i = 0; i < Lx; i++)
	{
		*(x + i + 1) = *(input + i) / 32768.0;
		*(y_f + i + 1) = *(x + i + 1);
	}
	for (i = 1; i <= no_delays; i++) {
		amp[i] = 1 - i*ratio;
	}

	for (i = samp_delay; i < Ly; i++) {
		for (j = 1; j <= no_delays; j++) {
			if (i >(j*samp_delay))
				if (amp[j] > 0)
					*(y_f + i) += amp[j] * *(x + (i - j*samp_delay));
		}
	}
	for (i = 0; i < Ly; i++)
		if (ymax < fabs(y_f[i])) ymax = fabs(y_f[i]);
	for(i=0;i<Ly;i++)
		*(y + i) = y_f[i]/(ymax+0.3)*32768.0;

	if (x != NULL) {
		free(x); x = NULL;
	}
	if (y_f != NULL) {
		free(y_f); y_f = NULL;
	}
	return y;
}

void aptool::ReverbConv(const char *in1, const char *in2, const char *out, const double wet, const double dry) {
	double *x, *h, *y_f, ymax = 0, delay;
	int i, Lx, Lh, Ly;
	short  *y;

	wavRandW *wav = new wavRandW();
	wavRandW *wav1 = new wavRandW();
	wav->wavRead(in1);
	Lx = wav->indatasize / 2;
	x = (double*)calloc(Lx, sizeof(double)); assert(x);
	for (i = 0; i < Lx; i++)
		*(x + i) = wav->indata[i] / 32768.0;

	wav1->wavRead(in2);
	Lh = wav1->indatasize / 2;

	h = (double*)calloc(Lh, sizeof(double)); assert(h);
	for (i = 0; i < Lh; i++)
		*(h + i) = wav1->indata[i] / 32768.0;

	Ly = Lx + Lh - 1;
	y_f = (double*)calloc(Ly, sizeof(double)); assert(y_f);
	y = (short int*)calloc(Ly - Lh / 2, sizeof(short int)); assert(y);
	convReverb(x, Lx, h, Lh, y_f);
	for (i = 0; i < Ly; i++)
	{
		if (i < Lx)
		{
			*(y_f + i) = *(y_f + i) * wet + *(x + i)*dry;
			if (ymax < fabs(*(y_f + i)))
				ymax = fabs(*(y_f + i));
		}

	}
	delay = 1;
	for (i = 0; i < Ly - Lh / 2; i++)
	{

		if (i > Lx + Lh / 4)
		{
			*(y + i) = round(*(y_f + i) / ymax * 32768 * delay);
			delay -= 1 / Lh / 4;
		}
		else
			*(y + i) = round(*(y_f + i) / ymax * 32768);
	}

	wav->outbit = wav->inbit;
	wav->outfs = wav->infs;
	wav->outdatasize = Ly - Lh / 2;
	wav->outdata = y;
	wav->wavWrite(out);//��Ƶ���

	if (NULL != y) { free(y); y = NULL; }
	if (NULL != y_f) { free(y_f); y_f = NULL; }
	if (NULL != h) { free(h); h = NULL; }
	if (NULL != x) { free(x); x = NULL; }
}

short int* aptool::ReverbConv(const short int* input1, const int Lx, const short int* input2, const int Lh, const double wet, const double dry, int &Ly) {

	double *x, *h, *y_f, ymax = 0, delay;
	int i;
	short int *y;



	x = (double*)calloc(Lx, sizeof(double)); assert(x);
	for (i = 0; i < Lx; i++)
	{
		*(x + i) = input1[i] / 32768.0;

	}



	h = (double*)calloc(Lh, sizeof(double)); assert(h);
	for (i = 0; i < Lh; i++)
		*(h + i) = input2[i] / 32768.0;

	Ly = Lx + Lh - 1;
	y_f = (double*)calloc(Ly, sizeof(double)); assert(y_f);
	y = (short int*)calloc(Ly, sizeof(short int)); assert(y);
	convReverb(x, Lx, h, Lh, y_f);
	for (i = 0; i < Ly; i++)
	{
		if (i < Lx)
		{
			*(y_f + i) = *(y_f + i) * wet + *(x + i)*dry;
			if (ymax < fabs(*(y_f + i)))
				ymax = fabs(*(y_f + i));
		}
	}
	delay = 1;
	for (i = 0; i < Ly - Lh / 2; i++)
	{

		if (i > Lx + Lh / 4)
		{
			*(y + i) = round(*(y_f + i) / (ymax+0.1) * 32768 * delay);
			delay -= 1 / Lh / 4;
		}
		else
			*(y + i) = round(*(y_f + i) / (ymax+0.1) * 32768);
	}


	Ly -= Lh / 2;
	if (NULL != y_f) { free(y_f); y_f = NULL; }
	if (NULL != h) { free(h); h = NULL; }
	if (NULL != x) { free(x); x = NULL; }
	return y;
	
}


short int* aptool::Robotization(const short *input, const int Lx, short type, int fs, double RMfre, int outwin) {
	short *y;
	RobotSound *robotsound = new RobotSound();
	y = (short*)calloc(Lx, sizeof(short)); assert(y);
	robotsound->Robot(input, Lx, y, type, fs, RMfre, outwin);
	delete robotsound;
	return y;
}

short int* aptool::Whisper(const short *input, const int Lx, int &Ly,int p) {
	short *y;
	RobotSound *robotsound = new RobotSound();
	robotsound->WhisperSound(input, Lx, &y, Ly,p);
	delete robotsound;
	return y;
}


short int * aptool::VibratoProcesse(const short *input, const int Lx,const int fs, const double modfreq, const double width) {
	short *y;
	y = (short*)calloc(Lx, sizeof(short)); assert(y);
	VibratoEffect(input, Lx, fs,y, modfreq, width);
	return y;
}

short int* aptool::FormantChange(const short *input, const int Lx, const double warp_rate) {
	FS *formantS = new FS();
	double *x_f,*y_f,ymax=0;
	short *y;

	x_f = (double *)calloc(Lx, sizeof(double)); assert(x_f);
	y = (short*)calloc(Lx, sizeof(short)); assert(y);
	for (int i = 0; i < Lx; i++) {
		*(x_f + i) = *(input + i) / 32768.0;
	}
	formantS->FormantWarp(x_f, Lx, &y_f, warp_rate);

	for (int i = 0; i < Lx; i++)
		if (ymax < fabs(*(y_f + i))) ymax = fabs(*(y_f + i));
	for (int i = 0; i < Lx; i++)
		*(y + i) = *(y_f + i)/(ymax+0.3) * 32768.0;

	if (x_f) { free(x_f); x_f = NULL; }
	if (y_f) { free(y_f); y_f = NULL; }

	return y;

	
}

short int* aptool::Chorus(const short*input, const int Lx, const short num_chorus, int &Ly,const double *rato, const double *gain, const double *shift_time) {
	vector<double> temp;
	double maxtime = 0,ymax=0;
	int Ly1, fs = 16000, delaytime;
	short *yout,*y;
	
	
	for (int i = 0; i < Lx; i++)
		temp.push_back((double)*(input + i));
	for (int i = 0; i < num_chorus; i++)
		maxtime = maxtime > shift_time[i]? maxtime:shift_time[i];
	delaytime = maxtime*fs;
	Ly = Lx + delaytime;
	y = (short*)calloc(Ly, sizeof(short)); assert(y);
	for (int i = Lx; i < Ly; i++)
		temp.push_back(0.0);
	for (int i = 0; i < num_chorus; i++) {
		yout = PandFshift(input, Lx, rato[i], Ly1);
		delaytime = shift_time[i] * fs;
		for (int j = delaytime; j < Lx + delaytime; j++)
			temp[j] += *(yout + j - delaytime)*gain[i];
		if (yout) { free(yout); yout = NULL; }
	}

	for (int i = 0; i < Ly; i++)
		ymax = ymax > temp[i] ? ymax : temp[i];

	for (int i = 0; i < Ly; i++)
		*(y + i) = temp[i] / (0.3*32768.0 + ymax)*32768.0;
	return y;

}

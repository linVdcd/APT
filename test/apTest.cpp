
#include"APtool.h"
#include"wavRandW.h"


#include <string>
using namespace std;

int main()
{
	double a = 44100.0/16000.0;
	//aptool *ap = new aptool();
	//ap->Resample("/home/lin/11.wav", "/home/lin/lll2.wav", a);//��Ҫ������Ƶ���ݵĳ���Ly
	//delete ap;
	aptool *ap = new aptool();
    ap->TimeScaling("/home/lin/HS-works/音频处理/AudioProcessingTool/APtoolC/data/wavfile/inputfiles/女声唱歌.wav","/home/lin/HS-works/音频处理/AudioProcessingTool/APtoolC/女声唱歌17.wav",1.7);
	string name = "/home/lin/HS-works/音频处理/AudioProcessingTool/APtoolC/data/config file examples/param_test.txt";
	//ap->ProcesseByOrder(name, "/home/lin/HS-works/音频处理/AudioProcessingTool/APtoolC/data/wavfile/inputfiles/女声唱歌.wav", "/home/lin/HS-works/音频处理/AudioProcessingTool/APtoolC/女声唱歌2.wav", "/home/lin/HS-works/音频处理/AudioProcessingTool/APtoolC/data/wavfile/im16/");
	//ap->ProcesseByOrder(name, "/home/lin/1.wav", "/home/lin/lll2.wav", "data/wavfile/m16/");
	//test
	//int a = 0;
	/*double ac[8] = { 0,0,0,-10,0,0,0,0 }, ac2[10] = { 0,0,0,0,10 * 0.6,0,0,0,0,0 }, Q[10] = { 2,2,2,2,2,2,2,2,2,2 };
	int Lx, Ly = 0, Lh;
	short *y;
	wavRandW *wav = new wavRandW();
	aptool *ap = new aptool();
	wav->outbit = 16;
	wav->outfs = 16000;
	//Lx���������ݵĳ���
	//���к����У����û��Ly�����ģ�Ĭ����������ݳ���ΪLx������ΪLy��


	//���룬���ڴ�������ʽ���롣����wav�ļ���ʽ���룺ap->DeNoise_martin("..\\ns.wav", "..\\ds.wav")
	wav->wavRead("..\\4mcra.wav");
	Lx = wav->indatasize / 2;
	y = ap->Whisper(wav->indata, Lx,32);
	wav->outdatasize = Lx;
	wav->outdata = y;
	wav->wavWrite("..\\robot.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}*/

	/*
	double ac[8] = { 0,0,0,-10,0,0,0,0 }, ac2[10] = { 0,0,0,0,10 * 0.6,0,0,0,0,0 }, Q[10] = { 2,2,2,2,2,2,2,2,2,2 };
	int Lx,Ly=0,Lh;
	short *y;
	wavRandW *wav = new wavRandW();
	aptool *ap = new aptool();
	wav->outbit = 16;
	wav->outfs = 16000;
	//Lx���������ݵĳ���
	//���к����У����û��Ly�����ģ�Ĭ����������ݳ���ΪLx������ΪLy��


	//���룬���ڴ�������ʽ���롣����wav�ļ���ʽ���룺ap->DeNoise_martin("..\\ns.wav", "..\\ds.wav")
	wav->wavRead("..\\С����Ů�� ԭ��lll.wav");
	Lx = wav->indatasize / 2;
	y = ap->DeNoise_martin(wav->indata, Lx);
	wav->outdatasize = Lx;
	wav->outdata = y;
	wav->wavWrite("..\\ds.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}
	//ap->DeNoise_martin("..\\20160801 174756 - ����ԭ�� - �칫������.wav", "..\\ds1.wav");

	/*
	//�����������ڴ�������ʽ���롣����wav�ļ���ʽ���룺EQ2("..\\ns.wav", "..\\eq2.wav", ac2,Q);
	// ac��10��Ƶ�����������ֵ
	// Q��10��Ƶ�������Qֵ
	wav->wavRead("..\\ns.wav");
	Lx = wav->indatasize / 2;
	y=ap->EQ2(wav->indata,Lx,ac2,Q);
	wav->outdatasize = Lx;
	wav->outdata = y;
	wav->wavWrite("..\\eq.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}




	//�ñ����ߣ������Ҳ�ı䣩�����ڴ�������ʽ���롣����wav�ļ���ʽ���룺PandFshift("..\\ns.wav", "..\\P&F.wav", 1.2);
	// rate: ��������
	wav->wavRead("..\\ds.wav");
	Lx = wav->indatasize / 2;
	y = ap->PandFshift(wav->indata, Lx, 1.667,Ly);//��Ҫ������Ƶ���ݵĳ���Ly
	wav->outdatasize = Ly;
	wav->outdata = y;
	wav->wavWrite("..\\P&F.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}

	//ʱ��߶����ţ����ڴ�������ʽ���롣����wav�ļ���ʽ���룺TimeScaling("..\\ns.wav", "..\\Tmod.wav", 0.5);
	wav->wavRead("..\\ns.wav");
	Lx = wav->indatasize / 2;
	y = ap->TimeScaling(wav->indata, Lx, 1.5, Ly);//��Ҫ������Ƶ���ݵĳ���Ly
	wav->outdatasize = Ly;
	wav->outdata = y;
	wav->wavWrite("..\\Tmod.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}
	//�����������ڴ�������ʽ���롣����wav�ļ���ʽ���룺limiter2("..\\ns.wav", "..\\limiter.wav",-12,1.3,-40,-0.001);
	wav->wavRead("..\\ns.wav");
	Lx = wav->indatasize / 2;
	y = ap->limiter2(wav->indata, Lx , -12, 1.3, -40, -0.001);
	wav->outdatasize = Lx;
	wav->outdata = y;
	wav->wavWrite("..\\limiter.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}

	//ѹ���������ڴ�������ʽ���롣����wav�ļ���ʽ���룺compressing("..\\ns.wav", "..\\com_3.wav", 3, -12, 10, 10);
	wav->wavRead("..\\ns.wav");
	Lx = wav->indatasize / 2;
	y = ap->compressing(wav->indata, Lx, 3, -12, 10, 10);
	wav->outdatasize = Lx;
	wav->outdata = y;
	wav->wavWrite("..\\com.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}
	//�ӳ��������ڴ�������ʽ���롣����wav�ļ���ʽ���룺delay("..\\ns.wav", "..\\delay.wav", 0.4);
	wav->wavRead("..\\ns.wav");
	Lx = wav->indatasize / 2;
	y = ap->delay(wav->indata, Lx, 0.4,Ly);//��Ҫ������Ƶ���ݵĳ���Ly
	wav->outdatasize = Ly;
	wav->outdata = y;
	wav->wavWrite("..\\delay.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}

	
	
	//�����������ڴ�������ʽ���롣����wav�ļ���ʽ���룺Reverb_conv("..\\ds.wav", "..\\im16\\Vocal Duo.wav", "..\\reverb.wav",0.2,1);
	wavRandW *wav1 = new wavRandW();
	wav->wavRead("..\\20160801 174756 - ����ԭ�� - �칫������.wav");
	Lx = wav->indatasize / 2;
	wav1->wavRead("..\\im16\\Vocal Duo.wav");//����Դ
	Lh = wav1->indatasize / 2;
	y = ap->Reverb_conv(wav->indata, Lx,wav1->indata,Lh,0.005,1,Ly);//��Ҫ������Ƶ���ݵĳ���Ly
	wav->outdatasize = Ly;
	wav->outdata = y;
	wav->wavWrite("..\\reverb_cut.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}

	
	
	wav->wavRead("..\\dashao.wav");
	Lx = wav->indatasize / 2;
	y = ap->TimeScaling(wav->indata, Lx, 0.8, Ly);//��Ҫ������Ƶ���ݵĳ���Ly
	wav->outdatasize = Ly;
	wav->outdata = y;
	wav->wavWrite("..\\Tmod.wav");
	if (y != NULL) {
		free(y); y = NULL;
	}*/
	//delete ap;
	return 0;
}

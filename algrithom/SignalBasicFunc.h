#ifndef SIGNALBASICFUNC_H
#define SIGNALBASICFUNC_H
#define PI 3.14159265358979323846 
class SignalBasicFunc{
public:
	void hamming(int n, double* w);
	void filter(int ord, double *a, double *b, int np, double *x, double *y);
	void fconv(double *x, int Lx, double *h, int Lh, double *y);
	double first_modified_Bessel(int n, double x);
	double sinc(double x);
	void kaiser(short N, double beta, double *window);
	void firls(int N, double *F, short *M, double *h);
	void upfirdn(
		double y[], unsigned int Ly, unsigned int ky,
		double x[], unsigned int Lx, unsigned int kx,
		double h[], unsigned int Lh, unsigned int kh,
		int p,
		int q
	);
	void hanning(int n, double* w);
};

#endif // !BASICFUNC
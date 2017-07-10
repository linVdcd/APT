#ifndef FORMANTSHIFT_H
#define FORMANTSHIFT_H

class FS {
public:
	void FormantWarp(double *input, const int Lx, double  **DAFx_out, const double warping_coef);

};

#endif // !FORMANTSHIFT_H
/**************************************

林明安 2016.6.3
LPC分析与合成

*************************************/
#include"lpc.h"

Lpc::Lpc() {}
Lpc::~Lpc() { delete sbf; }


char *getmem(const size_t leng, const size_t size)
{
	char *p = NULL;

	if ((p = (char *)calloc(leng, size)) == NULL) {
		fprintf(stderr, "Cannot allocate memory!\n");
		exit(3);
	}
	return (p);
}
OsFlt64 *dgetmem(OsInt32 leng)
{
	return ((OsFlt64 *)getmem((size_t)leng, sizeof(OsFlt64)));
}
void acorr(OsFlt64 *x, OsInt32 l, OsFlt64 *r, OsInt32 np)
{
	OsFlt64 d;
	OsInt32 k, i;

	for (k = 0; k <= np; k++) {
		for (d = 0.0, i = 0; i<l - k; i++)
			d += x[i] * x[i + k];
		r[k] = d;
	}

	return;
}
void levdur(OsFlt64 *r, OsFlt64 *a, OsInt32 m, OsFlt64 eps)
{
	OsInt32 l, k;
	OsFlt64 rmd, mue;
	static OsFlt64 *c = NULL;
	static OsInt32 size;

	if (c == NULL) {
		c = dgetmem(m + 1);
		size = m;
	}

	if (m>size) {
		free(c);
		c = dgetmem(m + 1);
		size = m;
	}

	if (eps<0.0) eps = 1.0e-6;
	rmd = r[0];

	a[0] = 0.0;

	for (l = 1; l <= m; l++) {
		mue = -r[l];
		for (k = 1; k<l; k++)
			mue -= c[k] * r[l - k];
		mue = mue / rmd;

		for (k = 1; k<l; k++)
			a[k] = c[k] + mue * c[l - k];
		a[l] = mue;

		rmd = (1.0 - mue * mue) * rmd;


		for (k = 0; k <= l; k++) c[k] = a[k];
	}
	a[0] = sqrt(rmd);

}
void lpc(OsFlt64 *x, OsInt32 flng, OsFlt64 *a, OsInt32 m, OsFlt64 f) {
	/****************************************************************

	LPC Analysis Using Levinson-Durbin method

	int lpc(x, flng, a, m);

	double  *x   : input sequence
	int     flng : flame length
	double  *a   : LP coefficients
	int     m    : order of LPC
	double  f    : mimimum value of the determinant
	of the normal matrix

	return value :  0  -> normally completed
	-1 -> abnormally completed
	-2 -> unstable LPC

	******************************************************************/

	static OsFlt64 *r = NULL;
	static OsInt32 size;
	if (r == NULL) {
		r = dgetmem(m + 1);
		size = m;
	}
	if (m>size) {
		free(r);
		r = dgetmem(m + 1);
		size = m;
	}

	acorr(x, flng, r, m);
	levdur(r, a, m, f);
	a[0] = 1;
}










void Lpc::lpcfit(OsFlt64 *x,OsInt32 Lx,OsInt16 p,OsInt16 h,OsFlt64 **a,OsFlt64 **g,OsFlt64 **e,OsInt16 *nframe){
	OsFlt64 *tempx,*tempx1, *temp_e,*xx, *windows,*aa,*rs,pre[2] = { 1,-0.9 }, *filter_a,temp;
	OsInt32 nhops = floor(Lx / h),hop,i;
	OsInt16 w = 2 * h;
	*nframe = nhops;
	tempx = (OsFlt64*)calloc(Lx + h, sizeof(OsFlt64)); assert(tempx);
	tempx1 = (OsFlt64*)calloc(Lx + h, sizeof(OsFlt64)); assert(tempx1);
	memcpy(tempx1 + h/2, x, sizeof(OsFlt64)*Lx);
	*a = (OsFlt64*)calloc(nhops*(p + 1), sizeof(OsFlt64)); assert(a);
	
	*g = (OsFlt64*)calloc(nhops, sizeof(OsFlt64)); assert(g);
	temp_e = (OsFlt64*)calloc((nhops + 1)*h, sizeof(OsFlt64)); assert(temp_e);
	*e = (OsFlt64*)calloc((nhops + 1)*h - h / 2, sizeof(OsFlt64)); assert(e);
	aa = (OsFlt64*)calloc(p + 1, sizeof(OsFlt64)); assert(aa);
	xx = (OsFlt64*)calloc(w, sizeof(OsFlt64)); assert(xx);
	rs = (OsFlt64*)calloc(w, sizeof(OsFlt64)); assert(rs);
	windows = (OsFlt64*)calloc(w, sizeof(OsFlt64)); assert(windows);
	sbf->hanning(w, windows);
	filter_a = (OsFlt64*)calloc(2, sizeof(OsFlt64)); assert(filter_a);
	filter_a[0] = 1;
	
	sbf->filter(1, filter_a, pre, Lx + h-1, tempx1, tempx);
	if (filter_a) { free(filter_a); filter_a = NULL; }
	filter_a = (OsFlt64*)calloc(p + 1, sizeof(OsFlt64)); assert(filter_a);
	filter_a[0] = 1;

	for (hop = 1; hop <= nhops; hop++) {
		memcpy(xx, tempx + (hop - 1)*h, sizeof(OsFlt64)*w);
		for (i = 0; i < w; i++)
			*(xx + i) *= *(windows + i);
		lpc(xx, w, aa, p, 0.000001);
		for (i = 0; i <= p;i++)
			if (isnan(aa[i]))
				aa[i] = 0.0;
		memcpy(*a + ((hop - 1)*(p + 1)), aa, sizeof(OsFlt64)*(p + 1));

		
		sbf->filter(p, filter_a, aa,w-1,xx,rs);
		//for (i = 0; i < w; i++)
		//	printf("%f ", rs[i]);
		temp = 0;
		for (i = 0; i < w; i++)
			temp += *(rs + i) * *(rs + i);
		*(*g+hop-1) = sqrt(temp / w);
		for (i = 0; i < w; i++) {
			*(temp_e + (hop - 1)*h + i) += *(rs + i) / *(*g + hop - 1);
		}
	}
	memcpy(*e, temp_e + h/2, sizeof(OsFlt64)*((nhops + 1)*h - h / 2));
	if (tempx) { free(tempx); tempx = NULL; }
	if (temp_e) { free(temp_e); temp_e = NULL; }
	if (aa) { free(aa); aa = NULL; }
	if (windows) { free(windows); windows = NULL; }
	if (rs) { free(rs); rs = NULL; }
	if (filter_a) { free(filter_a); filter_a = NULL; }
	if (tempx1) { free(tempx1); tempx1 = NULL; }
	if (xx) { free(xx); xx = NULL; }
}
void Lpc::lpcsynth(OsFlt64 *a, OsFlt64 *g, OsFlt64 *e, OsInt16 p,OsInt16 h,OsInt16 nframe, OsFlt64 **y,OsInt32 &Ly){
	OsFlt64 *temp_e,*tempy, *oldbit, *newbit, *aa, *windows, *ee,*filter_a, pre[2] = {1,-0.9};
	OsInt32 npts = nframe * h - h / 2,hop,w=2*h,i;
	temp_e = (OsFlt64*)calloc(npts + 3 * h, sizeof(OsFlt64)); assert(temp_e);
	memcpy(temp_e, e, sizeof(OsFlt64)*(npts + h));
	Ly = npts + h + h / 2;
	*y = (OsFlt64*)calloc(npts + h + h / 2, sizeof(OsFlt64)); assert(y);
	tempy = (OsFlt64*)calloc(npts + h + h / 2, sizeof(OsFlt64)); assert(tempy);
	oldbit = (OsFlt64*)calloc(w, sizeof(OsFlt64)); assert(oldbit);
	newbit = (OsFlt64*)calloc(w, sizeof(OsFlt64)); assert(newbit);
	aa = (OsFlt64*)calloc(p + 1, sizeof(OsFlt64)); assert(aa);
	windows = (OsFlt64*)calloc(w, sizeof(OsFlt64)); assert(windows);
	ee = (OsFlt64*)calloc(w, sizeof(OsFlt64)); assert(ee);
	sbf->hanning(w, windows);
	filter_a = (OsFlt64*)calloc(p + 1, sizeof(OsFlt64)); assert(filter_a);
	filter_a[0] = 1;
	
	for (hop = 1; hop <= nframe; hop++) {
		memcpy(oldbit, tempy + (hop - 1)*h, sizeof(OsFlt64)*h);
		memcpy(aa, a + (hop - 1)*(p+1), sizeof(OsFlt64)*(p + 1));
		memcpy(ee, temp_e + (hop - 1)*h , sizeof(OsFlt64)*w);
		sbf->filter(p, aa,filter_a , w-1, ee, newbit);
		for (i = 0; i < w; i++)
			*(tempy + (hop - 1)*h + i) = *(oldbit + i) + *(newbit + i) * *(windows + i) * *(g+hop-1);
	}
	if (filter_a) { free(filter_a); filter_a = NULL; }
	filter_a = (OsFlt64*)calloc(2, sizeof(OsFlt64)); assert(filter_a);
	filter_a[0] = 1;

	sbf->filter(1, pre, filter_a, npts+h+h/2-1, tempy, *y);

	if (temp_e) { free(temp_e); temp_e = NULL; }
	if (oldbit) { free(oldbit); oldbit = NULL; }
	if (newbit) { free(newbit); newbit = NULL; }
	if (aa) { free(aa); aa = NULL; }
	if (ee) { free(ee); ee = NULL; }
	if (filter_a) { free(filter_a); filter_a = NULL; }
	if (tempy) { free(tempy); tempy = NULL; }
	if (windows) { free(windows); windows = NULL; }
	

}

	/*********************************************
	
		长信号分成多帧的lpc分析与合成测试
	
	*********************************************/
/*int main(){



	OsFlt64 x[100];
	OsFlt64 *a =NULL, *g = NULL, *e = NULL, *y,esp = 0.00001;
	OsInt16 i, nframe=0;
	for (i = 0; i < 100; i++) {
		x[i] = i + 1;
	}

	lpcfit(x, 100,3,10, &a,&g,&e,&nframe);
	lpcsynth(a, g, e, 3, 10, nframe, &y);
	for (i = 0; i < 110; i++)
		printf("%f ", y[i]);
	return 0;
}*/







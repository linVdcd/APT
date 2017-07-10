/*-----------------------------------
	林明安 2016.4.25
	Pitch shifting based on resample;
*/
#include"Resample.h"
#include"SignalBasicFunc.h"
#ifndef max   
#define max(A, B) ((A)>(B) ? (A) : (B))   
#endif   
#ifndef min   
#define min(A, B) ((A)<(B) ? (A) : (B))   
#endif  
//OsFlt32 pi = 3.1415926;

OsInt32 gcd(
	OsInt32 a,
	OsInt32 b
	)
{
	OsInt32  t;
	while (b>0)  {
		t = b;
		b = a%b;
		a = t;
	}
	return(a);
}




OsFlt64 *resample(OsFlt64 *input, const OsUInt32 Lx, const OsInt32 P, const OsInt32 Q) {
	SignalBasicFunc *sbf = new SignalBasicFunc();
	OsFlt64 fc = 0.0f, *h, *w, F[4] = { 0.0f }, sum = 0.0f, *H, *y;
	OsInt32 maxgcd = gcd(P, Q), temp_p = P / maxgcd, temp_q = Q / maxgcd, pqmax = 0,
		i, nz, delay;
	OsUInt32 L = 0, Lhalf, nz1 = 0, Lh, Ly;
	OsInt16 bta = 5, N = 10, M[4] = { 1,1,0,0 };
	
	
	if ((1 == temp_p) && (1 == temp_q)) {
		y = (OsFlt64*)malloc(sizeof(OsFlt64)*Lx); assert(y);
		memcpy(y, input, Lx * sizeof(OsFlt64));
		return y;
	}
	pqmax = temp_p > temp_q ? temp_p : temp_q;
	fc = 0.5 / (OsFlt64)pqmax;
	L = 2 * N*pqmax + 1;
	h = (OsFlt64*)calloc(L, sizeof(OsFlt64)); assert(h);
	w = (OsFlt64*)calloc(L, sizeof(OsFlt64)); assert(w);
	F[0] = 0.0f; F[1] = F[2] = fc; F[3] = 0.5f;
	sbf->firls(L, F, M,h);//获取滤波系数
	sbf->kaiser(L, bta,w);//获取窗函数
	for (i = 0; i < L; i++)
	{
		*(h + i) *= *(w + i);
		sum += *(h + i);
	}
	for (i = 0; i < L; i++)
		h[i] *= temp_p / sum;
	Lh = L;
	Lhalf = (L - 1) / 2;
	nz = floor(temp_q - Lhalf%temp_q);
	Lhalf += nz;
	delay = floor(Lhalf / temp_q);
	Lh += nz;
	while (ceil(((Lx - 1)*temp_p + Lh + nz1) / temp_q) - delay < ceil(Lx*temp_p / temp_q))
		nz1++;
	Lh += nz1;
	H = (OsFlt64*)malloc(sizeof(OsFlt64)*Lh); assert(H);
	memset(H, 0, sizeof(OsFlt64)*nz);
	memcpy(H + nz, h, sizeof(OsFlt64)*L);

	memset(H + nz + L, 0, sizeof(OsFlt64)*nz1);

	if (h != NULL)
	{
		free(h);
		h = NULL;
	}

	Ly = temp_p*(Lx - 1) + Lh;
	Ly = (Ly%temp_q) ? Ly / temp_q + 1 : Ly / temp_q;
	y = (OsFlt64 *)malloc(sizeof(OsFlt64)*Ly);

	memset(y, 0, Ly * sizeof(OsFlt64));
	sbf->upfirdn(y, Ly, 1, input, Lx, 1, H, Lh, 1, temp_p, temp_q);//多相分解
	if (H != NULL)
	{
		free(H);
		H = NULL;
	}
	if (w != NULL) {
		free(w);
		w = NULL;
	}
	return y;
}



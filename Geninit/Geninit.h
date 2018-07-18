#ifndef GENINIT_H_
#define GENINIT_H_

class Geninit
{

	int vprocs[3], mc[3];
	int nprocs, mctot, aty, gid, ity;
	int *lnatoms, *lnatoms1, *lnatoms2;
	int *itype0;
	int natoms, ntot, ix, iy, iz;

	double pi;
	
	double L1, L2, L3, Lalpha, Lbeta, Lgamma;
	double lbox[3], obox[3], dtype;
	double H[3][3], Hi[3][3], rr[3], rr1[3], vv[3], qq, rmin[3], rmax[3];
	double qfsp, qfsv;
	double **pos0, **pos1;
	double *itype1, **rr3;

	char **ctype0, **ctype1, fnote[256];

public:

	Geninit(char*, int, int, int, int, int, int);
	void HMatrix();
	void Hinv();
	void Replicate(int);
	void IntroduceVoid(double, double, double);
	void WrapBack(int,int,int);
	void ErrorCheck();
	void OutputBin();	
};

#endif

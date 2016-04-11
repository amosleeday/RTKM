#pragma once
#include "stdafx.h"//changed
#include "VariableDef.h"
#include "ExtFun.h"
#include <math.h>
#include "matrix.h"






class Ambiguity
{
	
	
public:

	void TCARPre(double* obsP,double* obsL,int sysid, double* NarrAmb);
	bool IsReq(double delta,double Wave,double sigma,double PeReq);
	double PrErrRound(double delta,double Wave,double sigma);
	bool CIR(double* obsP,double* obsL,int sysid, double* NarrAmb);
	void Gauss(int n, math::matrix<double> & L, math::matrix<double> & Z, int i, int j);
	void Perm(int n, math::matrix<double> &L, double* D, int j, double del, math::matrix<double> &Z);
	void Reduction(int n, math::matrix<double>& L, double* D, math::matrix<double> &Z);
	//int search(int n, int m, math::matrix<double> L, double* D, double* zs,math::matrix<double>& zn, double* s);
	int search(int n, int m, math::matrix<double> L, double* D, math::matrix<double> zs,math::matrix<double>& zn, double* s);
	int Lambda(int n, int m, math::matrix<double> a, math::matrix<double> Q, math::matrix<double> & F, double* s);
	bool LTDLFactorization(int n, math::matrix<double> Q, math::matrix<double> &L, double* D);
};



class Combination
{

public:
	void		Freq(int sysid,double* freq);
	void		FreqRatio(int sysid,double* freq);
	void		Wave(int sysid,double* Wave);
	double CombFreq(int sysid,int* coef);
	double CombWave(int sysid,int* coef);
	double CombObs(int sysid,int* coef,double* threeobs);
	double CombIono(int sysid,int* coef);
	double	IonoFree(int sysid,double* obs,int* index);
};
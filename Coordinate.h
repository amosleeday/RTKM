#include "stdafx.h"//changed
#include "VariableDef.h"


/*
	This file define the class about the trans of coordinate 
*/

class CrdTrans
{
	void EarthConst(int ellipInd,double& a,double& eE);
	void BLH2XYZ(double* BLH,int ellipInd,double* XYZ);
	void XYZ2BLH(double* XYZ,int ellipInd,double* BLH);
	void NEU2XYZ(double* NEU,double B0,double L0,double* XYZ);
	void XYZ2NEU(double* XYZ,double B0,double L0,double* NEU);
	double azi(double dn,double de);
	void	NEU2RAH(double *NEU,double* RAH);
	void RAH2NEU(double* RAH,double* NEU);
	void GaussProj(double B,double L,double L0,int ellipInd,double* xy);
};














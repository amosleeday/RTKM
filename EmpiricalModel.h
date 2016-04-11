#include "stdafx.h"//changed
#include "VariableDef.h"
#include "ExtFun.h"
#include "Rinex.h"

/*
	This file define the struct and class about the emprirical model
		1.Atmosphere
			1.1	ionosphere
			1.2	tropsphere

		2.Geoid
			EGM2008
			to be written


*/

class Ionos
{
public:
	double	Klobuchar(BroadEphIonoCorr para,double Sow,double* pos, 
										double azi,double ele);

	double KlobucharBDS(BroadEphIonoCorr para,double Sow,double* pos, double azi,double ele);
	double MapFunc(double* pos,double azi,double ele,double re,double hion);

	double Pierce(double *pos, double azi,double ele, double re,double hion, double *pierce);

};


class Trops
{
public:
	double InterpLinear1(double* x,double lati);
	double InterpLinear2(double* x,double lati);
	double calcos(double p0,double p1,int doy,int Dmin);
	void		ZTDUNB3(double lat,double level,int doy,double* ZTD,
									double* Dd,double* Dw);
	//double SlantUNB3(double ele,double* Dd,double* Dw,double* ZTD);
	double SlantUNB3(double ele,double* Dd,double* Dw,double* ZTD,double& Mwet);
	//double	TropsUNB3(double ele,double lat,double level,int doy);
	double TropsUNB3(double ele,double lat,double level,int doy,double& Mwet);
};














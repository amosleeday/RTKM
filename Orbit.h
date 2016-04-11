#include "stdafx.h"//changed
#include "VariableDef.h"
#include "ExtFun.h"
#include "Rinex.h"



/*
This file defines the struct and class that may be involved in the orbit

Also,it includes the method of computing the satellite position
	1.	broadcast eph
		1.1.	GPS/GAL/BDS
		1.2.	GLONASS is to be written
	
	2.	precise eph
		to be written


*/

/*
	
*/
class CalcSatePos
{
public:
	void CompSatePosVelAcc(math::matrix<double>& SatePos,math::matrix<double>& SateVel,math::matrix<double>& SateAcc,double& ek,BroadEphData eph,double t);
	void CompSateClk(double& SateClk,double ek,BroadEphData eph,double transt,int freqIndex);
	void	InteSatePosClk(math::matrix<double>& SatePos,double& SateClk,BroadEphData eph,math::matrix<double> StnPos,double t,int freqIndex);
	void	InteSatePosClkVA(math::matrix<double>& SatePos,double& SateClk,BroadEphData eph,math::matrix<double> StnPos,double t,double t0,double ek,math::matrix<double> Xcrd,math::matrix<double> Xvel,math::matrix<double> Xacc);
	bool SatePosClk(double* pos,double& clkerr,double* RecPos,double t,int freqIndex,BroadEphData eph);
	double geph2clk(double ts, const BroadEphDataGlo geph);
	void CalcSatePos:: geph2pos(double ts, const BroadEphDataGlo geph, double *rs, double& dts);
	double CalcSatePos:: seph2clk(double time, const BroadEphDataSBAS seph);
	void seph2pos(double time, const BroadEphDataSBAS seph, double *rs, double *dts);
};




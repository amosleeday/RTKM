#include "stdafx.h"//changed
#include "VariableDef.h"
#include<math.h>
#include "Coordinate.h"


extern double norm(double* a,int n)
{
	double count=0.0;
	for(int i=0;i<n;i++)
	{
		count+=a[i]*a[i];
	}
	return sqrt(count);
}


/*
	Get the earth constant of different ellipsoid
		I:
				ellipInd		index of ellipsoid 
								1.	WGS 84;		2.	CGCS2000
								3.	BJ54;				4.	XI'AN80
	    O:
				a				semimajor axis
				eE				e^2	, ee is "#define"
*/

void CrdTrans::EarthConst(int ellipInd,double& a,double& eE)
{

	if (ellipInd==1)	
	{
		a		=	RE_WGS84;
		eE		=	ee;
	}
	else if(ellipInd==2)	
	{
		a		=	RE_CGCS2000;
		eE		=	ee_BDS;
	}
	else if(ellipInd==3)	
	{
		a		=	RE_BJ54;
		eE		=	ee_BJ54;
	}
	else if(ellipInd==4)	
	{
		a		=	RE_XA80;
		eE		=	ee_XA80;
	}
}

/*
	BLH2XYZ transfer coordiante BLH to XYZ in ECEF
		I:
			BLH			3*1  (rad,rad,m)
			ellipInd		index of ellipsoid 
	   O:
			XYZ			3*1(m,m,m)
*/

void CrdTrans::BLH2XYZ(double* BLH,int ellipInd,double* XYZ)
{
	double a,eE;
	EarthConst(ellipInd,a,eE);

	double N	=a/sqrt(1.0-eE*pow(sin(BLH[0]),2));

	XYZ[0]		=(N+BLH[2])*cos(BLH[0])*cos(BLH[1]);
	XYZ[1]		=(N+BLH[2])*cos(BLH[0])*sin(BLH[1]);
	XYZ[2]		=(N*(1-eE)+BLH[2])*sin(BLH[0]);
}

/*
	XYZ to BLH
		I:
			XYZ
			ellipInd
		O:
			BLH
*/
void CrdTrans::XYZ2BLH(double* XYZ,int ellipInd,double* BLH)
{
	double r		=	sqrt(pow(XYZ[0],2)+pow(XYZ[1],2));
	double	B0		=	atan(	(XYZ[2])/r);
	
	double a,eE,N,B;
	EarthConst(ellipInd,a,eE);

	while(1)
	{
		N	=	a / sqrt(1.0 - eE * pow(sin(B0),2));
		B	=	atan((XYZ[2] + N * ee *sin(B0)) /r );

		if(fabs(B-B0)>4.8e-11)
		{
			B0=B;
		}
		else
		{
			break;
		}
	}
	BLH[0]	=	B;
	BLH[1] =	atan2(XYZ[1], XYZ[0]);
	N		   =	a/sqrt(1.0-eE*pow(sin(B),2));
	BLH[2] =	r/cos(B)-N;
}


/*
	NEU to XYZ
		I:
			NEU		3*1
	        B0		  latitude of topocenter        (rad)              
	        L0		  longitude of topocenter		(rad)
		O:
			XYZ		3*1
*/
void CrdTrans::NEU2XYZ(double* NEU,double B0,double L0,double* XYZ)
{
	XYZ[0]	=	-sin(B0)*cos(L0)*NEU[0]   -sin(L0)*NEU[1]+  cos(B0)*cos(L0)*NEU[2];
	XYZ[1]	=	-sin(B0)*sin(L0)*NEU[0]   +cos(L0)*NEU[1]+  cos(B0)*sin(L0)*NEU[2];
	XYZ[2]	=	cos(B0)*NEU[0]+sin(B0)*NEU[2];
}

/*
	XYZ ro NEU
		I:
			XYZ		difference of three componets
			B0,L0 (rad,rad)
		O:
			NEU
*/
void CrdTrans::XYZ2NEU(double* XYZ,double B0,double L0,double* NEU)
{
	NEU[0]		=	-sin(B0)*cos(L0)*XYZ[0]  -sin(B0)*sin(L0)*XYZ[1] + cos(B0)*XYZ[2];
	NEU[1]		=	-sin(L0)*XYZ[0]         +  cos(L0)*XYZ[1];
	NEU[2]		=	cos(B0)*cos(L0)*XYZ[0]   +cos(B0)*sin(L0)*XYZ[1]  +sin(B0)*XYZ[2];
}


/*
	get azimuth
		I:
			dn	difference if north component between two points
			de		difference if east component between two points
		return:
			azimuth(rad)
*/
double CrdTrans::azi(double dn,double de)
{
	if(dn==0)
	{
		return de>0?PI/2:3*PI/2;
	}
	else
	{
		double R	=	atan(de/dn);
		if(dn>0 )
		{
			return de>0?R:R+2*PI;
		}
		else
		{
			return R+PI;
		}
	}
}

/*
	NEU to RAH
		I:
			NEU 3*1
		O:
			RAH 3*1
			1.ploar range(m)
			2.azimuth (rad)
			3.elevation(rad)
*/
void	CrdTrans::NEU2RAH(double *NEU,double* RAH)
{
	RAH[0]		=	norm(NEU,3);
	RAH[1]		=	azi(NEU[0],NEU[1]);
	RAH[2]		= atan(NEU[2])/norm(NEU,2);
}


/*
	NEU to RAH
		I:
			RAH 3*1
			1.ploar range(m)
			2.azimuth (rad)
			3.elevation(rad)
		O:
			NEU 3*1
*/
void CrdTrans::RAH2NEU(double* RAH,double* NEU)
{
	NEU[0]		=	RAH[0]*cos(RAH[1])*cos(RAH[2]);
	NEU[1]		=	RAH[0]*sin(RAH[1])*cos(RAH[2]);
	NEU[2]		=	RAH[0]*sin(RAH[2]);
}

/*
	project  BL to xy in Gauss crd system 
		I:
			B,L		(rad,rad)
			L0			center longtitude for projection(rad)
			ellipind 
		O:
			xy  the Gauss Proj horizontal coordiantes
*/
void CrdTrans::GaussProj(double B,double L,double L0,int ellipInd,double* xy)
{
	double a,eE;
	EarthConst(ellipInd,a,eE);

	double C	=	a/sqrt(1.0-eE);

	double aa=1.0+3.0*pow(eE,2)/4.0+45.0*pow(eE,4)/64.0+175.0*pow(eE,6)/256.0
						+11025.0*pow(eE,8)/16384.0+43659.0*pow(eE,10)/65536.0;
	double bb=3.0*pow(eE,2)/4.0 + 15.0*pow(eE,4)/16.0 + 525.0*pow(eE,6)/512.0
						+2205.0*pow(eE,8)/2048.0 + 72765.0*pow(eE,10)/65536.0;
	double cc=15.0*pow(eE,4)/64.0 + 105.0*pow(eE,6)/256.0
						+ 2205.0*pow(eE,8)/4096.0 + 10395.0*pow(eE,10)/16384.0;
	double dd=35.0*pow(eE,6)/512.0 + 315.0*pow(eE,8)/2048.0 + 31185.0*pow(eE,10)/13072.0;

	double a1 =  aa * a* (1.0 - eE);
	double a2 = -bb* a * (1.0 - eE) / 2.0;
	double a3 =  cc *a * (1.0 - eE) / 4.0;
	double a4 = -dd* a* (1.0 - eE) / 6.0;

	double r0 = a1;
	double r1 = 2.0*a2 + 4.0*a3 + 6.0*a4;
	double r2 = -8.0*a3 - 32.0*a4;
	double r3 = 32.0*a4;

	double	Tb = tan(B);
	double Y2 = eE/(1.0-eE)*pow(cos(B),2);
	double N  = C/sqrt(1.0 + Y2);

	double X0 = r0*B+ cos(B) *sin(B) *(r1+pow(sin(B) ,2)*(r2+pow(sin(B) ,2)*r3));
	double M0 = cos(B) * (L-L0);
	double M2 = M0 * M0;
	xy[0] = X0 + M2 *N *Tb/2.0+ M2 *M2 *N *Tb/24.0*(5.0-Tb *Tb+9.0*Y2+4.0*Y2 *Y2)+   
			    M2 *M2 *M2 *N *Tb/720.0*(61.0+(Tb *Tb-58.0)*Tb *Tb);
	xy[1] = N *M0 + N *M0 *M2 /6.0*(1.0+Y2-Tb *Tb)+ N *M2 *M2 *M2 /120.0*   
				(5.0+(Tb *Tb-18.0)*Tb *Tb-(58.0*Tb *Tb-14.0) *Y2);

}



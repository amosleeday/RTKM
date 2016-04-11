#include "stdafx.h"//changed
#include "VariableDef.h"
#include "EmpiricalModel.h"





/*
	Klobuchar(...)
		para				8 parameters
		Sow				second of week
		pos				position of receiver  lat,lon,h(rad/rad/m)
		azi				azimuth	/rad
		ele				elevation	/rad

		return 			ionospheric delay (L1) (m)
*/

double	Ionos::Klobuchar(BroadEphIonoCorr para,double Sow,double* pos, 
										double azi,double ele)
{
	if (para.validA==0 || para.validB==0)
	{
		return 0.0;
	}
	double tt,f,psi,phi,lam,amp,per,x;
	 if (pos[2]<-1E3 || ele<=0) return 0.0;

	    /* earth centered angle (semi-circle) */
    psi=0.0137/(ele/PI+0.11)-0.022;

	/* subionospheric latitude/longitude (semi-circle) */
    phi=pos[0]/PI+psi*cos(azi);
    if      (phi> 0.416) phi= 0.416;
    else if (phi<-0.416) phi=-0.416;
    lam=pos[1]/PI+psi*sin(azi)/cos(phi*PI);

	/* geomagnetic latitude (semi-circle) */
    phi+=0.064*cos((lam-1.617)*PI);
    
    /* local time (s) */
    tt=43200.0*lam+Sow;
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */
    
    /* slant factor */
    f=1.0+16.0*pow(0.53-ele/PI,3.0);
    
    /* ionospheric delay */
    amp=para.alpha1+phi*(para.alpha2+phi*(para.alpha3+phi*para.alpha4));
    per=para.beta1+phi*(para.beta2+phi*(para.beta3+phi*para.beta4));
    amp=amp<    0.0?    0.0:amp;
    per=per<72000.0?72000.0:per;
    x=2.0*PI*(tt-50400.0)/per;
    
    return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);

}


	/*
		 compute ionospheric delay mapping function by single layer model
		 pos    same as klobuchar
		 azi
		 ele
		 re    earth radius (m)
		 hion height of model(m)
		 return : ionospheric mapping function
	*/
double Ionos::MapFunc(double* pos,double azi,double ele,double re,double hion)
{
	if (pos[2]>=hion) return 1.0;
	return 1.0/cos(asin((re+pos[2])/(re+hion)*sin(PI/2.0-ele)));
}


/* ionospheric pierce point position 
	 compute ionospheric pierce point (ipp) position and slant factor

          double re           earth radius (km)
          double hion       height of ionosphere (km)
          double *pierce     pierce point position {lat,lon,h} (rad,m)
*return : slant factor
*/
double Ionos::Pierce(double *pos, double azi,double ele, double re,
								double hion, double *pierce)
{
    double cosaz,rp,ap,sinap,tanap;
    
    rp=re/(re+hion)*cos(ele);
    ap=PI/2.0-ele-asin(rp);
    sinap=sin(ap);
    tanap=tan(ap);
    cosaz=cos(azi);
    pierce[0]=asin(sin(pos[0])*cos(ap)+cos(pos[0])*sinap*cosaz);
    
    if ((pos[0]> 70.0*D2R&& tanap*cosaz>tan(PI/2.0-pos[0]))||
        (pos[0]<-70.0*D2R&&-tanap*cosaz>tan(PI/2.0+pos[0]))) 
	{
        pierce[1]=pos[1]+PI-asin(sinap*sin(ele)/cos(pierce[0]));
    }
    else {
        pierce[1]=pos[1]+asin(sinap*sin(azi)/cos(pierce[0]));
    }
    return 1.0/sqrt(1.0-rp*rp);
}

double Ionos::KlobucharBDS( BroadEphIonoCorr para,double Sow,double* pos, double azi,double ele )
{
	if (para.validA==0 || para.validB==0)
	{
		return 0.0;
	}
	double psi	=PI/2-ele-asin(6378/(6378+375)*cos(ele));
	double phiPierce	=asin(sin(pos[0])*cos(psi)+cos(pos[0])*sin(psi)*cos(azi));
	double lamPierce	=pos[1]+asin(sin(psi)*sin(azi)/cos(phiPierce));
	phiPierce	=phiPierce/PI;
	lamPierce	=lamPierce/PI;
	double A2=0.0;
	double A4	=0.0;
	A2	= para.alpha1+phiPierce*(para.alpha2+phiPierce*(para.alpha3+phiPierce*para.alpha4));
	A4	= para.beta1+phiPierce*(para.beta2+phiPierce*(para.beta3+phiPierce*para.beta4));
	A2	=(A2>=0)?A2:0;
	if (A4>=172800.0)
	{
		A4	=172800.0;
	}
	else if (A4<72000.0)
	{
		A4	=72000.0;
	}
	double tt=43200.0*lamPierce+Sow;
	tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

	double ionodelay	=5E-9+A2*cos(2*PI*(tt-50400)/A4);
	ionodelay	=(fabs(tt-50400)<A4/4)?ionodelay:5E-9;

	return CLIGHT*ionodelay/sqrt(1-SQ(cos(ele)*6378/(6378+375)));
}



double Trops::InterpLinear1(double* x,double lati)
{
	double minlat=0.0;
	int    ind=0;
	if(lati<=15)
	{
		return x[0];
	}
	else if(lati>75)
	{
		return x[4];
	}
	else
	{
		if(lati>15 && lati<=30)
		{
			minlat=15;
		}
		if(lati>30 && lati<=45)
		{
			minlat=30;
			ind=1;
		}
		if(lati>45 && lati<=60)
		{
			ind=2;
			minlat=45;
		}
		if(lati>60 && lati<=75)
		{
			minlat=60;
			ind=3;
		}
		return x[ind]+(x[ind+1]-x[ind])/15.0*(lati-minlat);
	}
}

double Trops::InterpLinear2(double* x,double lati)
{
	int    ind=0;
	if(lati<=15)
	{
		return x[0];
	}
	else if(lati>75)
	{
		return x[4];
	}
	else
	{
		if(lati>15 && lati<=30)
		{
			
		}
		if(lati>30 && lati<=45)
		{
			ind=1;
		}
		if(lati>45 && lati<=60)
		{
			ind=2;
		}
		if(lati>60 && lati<=75)
		{
			ind=3;
		}
		return x[ind]+(x[ind+1]-x[ind])/15.0*(lati-15);
	}

}

double Trops::calcos(double p0,double p1,int doy,int Dmin)
{
	return p0-p1*cos(2*PI*(doy-Dmin)/365.25);
}


/*
	compute the ztd and coefficients of mapping function
	I:
		lat					the receiver latitude(degree)
		level				the level height of station (meter)
		doy				day of year
	O:
		ZTD				2*1 dry and wet component
		Dd				3*1 coef for mapfunc of dry component
		Dw				3*1 coef for mapfunc of wet component
*/
void	Trops::ZTDUNB3(double lat,double level,int doy,double* ZTD,
									double* Dd,double* Dw)
{
	double p0[]		={1013.25, 1017.25, 1015.75, 1011.75, 1013.00};//pressure
	double detap[]	={00.00,	 -3.75,	   -2.25,		-1.75,	  -0.50};
	double T0[]		={299.65,	 294.15,   283.15,	272.15,	  263.65};//temperature
	double detaT[]	={0.00,		 7.00,	   11.00,		15.00,	  14.50};
	
	double e0[]		={26.31,	 21.79,	   11.66,		6.78,		   4.11};//water vapor pressure
	double detae[]	={0.00,		 8.85,		7.24,		5.36,		   3.39};
	double Beta0[]	={6.30*1.0e-3,	6.05*1.0e-3,	5.58*1.0e-3,	5.39*1.0e-3,	4.53*1.0e-3};//temperature lapse rate
	double detaBeta[]={0.00*1.0e-3,	0.25*1.0e-3,	0.32*1.0e-3,	0.81*1.0e-3,	0.62*1.0e-3};
	
	double	Lamda0[]={2.77,	3.15,	2.57,	1.81,	1.55};//water vapor rate
	double detaLamda[]={0.00,0.33,0.46,0.74,0.30};

	double ad0[]		={1.2769934*1.e-3,	1.2683230*1.e-3,	1.2465397*1.e-3,	1.2196049*1.e-3,	1.2045996*1.e-3};//para of mapfunc
	double bd0[]		={2.9153695*1.e-3,	2.9152299*1.e-3,	2.9288445*1.e-3,	2.9022565*1.e-3,	2.9024912*1.e-3};
	double cd0[]		={62.610505*1.e-3,	62.837393*1.e-3,	63.721774*1.e-3,	63.824265*1.e-3,	64.258455*1.e-3};
	double detaad[]={0.0           ,				1.2709626*1.e-5,	2.6523662*1.e-5,	3.4000452*1.e-5,	4.1202191*1.e-5};
	double detabd[]={ 0.0          ,			2.1414979*1.e-5,	3.0160779*1.e-5,	7.2562722*1.e-5,	11.723375*1.e-5};
	double detacd[]={0.0           ,				9.0128400*1.e-5,	4.3497037*1.e-5,	84.795348*1.e-5,	170.37206*1.e-5};

	double aw[]		={5.8021897*1.e-4,	5.6794847*1.e-4,	5.8118019*1.e-4,	5.9727542*1.e-4,	6.1641693*1.e-4};
	double bw[]		={1.4275268*1.e-3,	1.5138625*1.e-3,	1.4572752*1.e-3,	1.5007428*1.e-3,	1.7599082*1.e-3};
	double cw[]		={4.3472961*1.e-2,	4.6729510*1.e-2,	4.3908931*1.e-2,	4.4626982*1.e-2,	5.4736038*1.e-2};

	double k1			=77.604;
	double k2			=382000.0;
	double Rd			=287.054;
	double gm		=9.784;
	double g			=9.80665;

	double lati		=lat>0?lat:fabs(lat);
	int   Dmin		=lat>0?28:211;

	double p,T,e,Beta,Lamda;

	p				=calcos(  InterpLinear1(p0,lati), InterpLinear1(detap,lati),doy, Dmin);
	T				=calcos(  InterpLinear1(T0,lati), InterpLinear1(detaT,lati),doy, Dmin);
	e				=calcos(  InterpLinear1(e0,lati), InterpLinear1(detae,lati),doy, Dmin);
	Lamda		=calcos(  InterpLinear1(Lamda0,lati), InterpLinear1(detaLamda,lati),doy, Dmin);
	Beta			=calcos(  InterpLinear1(Beta0,lati), InterpLinear1(detaBeta,lati),doy, Dmin);
	//cout<<p<<"  "<<T<<"  "<<e<<"  "<<Lamda<<"  "<<Beta<<endl;

	Dd[0]			=calcos(  InterpLinear2(ad0,lati), InterpLinear2(detaad,lati),doy, Dmin);
	Dd[1]			=calcos(  InterpLinear2(bd0,lati), InterpLinear2(detabd,lati),doy, Dmin);
	Dd[2]			=calcos(  InterpLinear2(cd0,lati), InterpLinear2(detacd,lati),doy, Dmin);
	Dw[0]			=InterpLinear2(aw,lati);
	Dw[1]			=InterpLinear2(bw,lati);
	Dw[2]			=InterpLinear2(cw,lati);

	double temp1	=g/Rd/Beta;
	ZTD[0]				=pow( (1-Beta*level/T),	temp1)*(1.0e-6*k1*Rd*p/gm);
	ZTD[1]				=pow( (1-Beta*level/T),	temp1*(Lamda+1)-1) * (1.0e-6*k2*Rd*e/(gm*(Lamda+1)-Beta*Rd)/T) ;
	//cout<<ZTD[0]<<ZTD[1];
	//cout<<endl;
}



/*
		calculate the slant troposheric delay using ZTD from NUB3 model
		 and mapping functions
		 I:
			 ele        rad
			 Dd*
			 Dw*
			 ZTD*   2*1
		 return: slant tropspheric dealy (m)
*/
double Trops::SlantUNB3(double ele,double* Dd,double* Dw,double* ZTD,double& Mwet)
{
	double Mdry		=	(1 + Dd[0]/ (1 + Dd[1]/(1+Dd[2]))) /(sin(ele) + Dd[0]/(sin(ele)+Dd[1]/(sin(ele)+Dd[2])));
	Mwet	=	(1+Dw[0]/(1+Dw[1]/(1+Dw[2])))/(sin(ele)+Dw[0]/(sin(ele)+Dw[1]/(sin(ele)+Dw[2])));
	return  ZTD[0]*Mdry	+	ZTD[1]*Mwet;
}

/*
 *Get the mapping function of wet delay, using in SPP
 *  ele  rad
 *  lat   degree
 */

double Trops::TropsUNB3(double ele,double lat,double level,int doy,double& Mwet)
{
	double* ZTD		=new double[2];
	double* Dd		=new double[3];
	double* Dw		=new double[3];
	ZTDUNB3(lat,level,doy,ZTD,Dd,Dw);
	double corr	=SlantUNB3(ele,Dd,Dw,ZTD,Mwet);
	delete[] ZTD;
	delete[] Dd;
	delete[] Dw;

	return corr;
}
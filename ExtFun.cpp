#pragma  once
#include "stdafx.h"//changed
#include<math.h>//calc
#include "ExtFun.h"
#include <iomanip>
#include <vector>






const static double gpst0[]={1980,1, 6,0,0,0}; /* gps time reference */
const static double gst0 []={1999,8,22,0,0,0}; /* galileo time reference */
const static double bdt0 []={2006,1, 1,0,0,0}; /* beidou time reference */

 /* 
  * leap seconds {y,m,d,h,m,s,UTC-GPST,...} , 
  * UTC-GPST=19-DTAI.
  * DTAI, published by International Earth Rotation Service (IERS)
  * More details, refer to the updated Report with LaTex
  */
const static int leaps[][7]=
{
	{2015,7,1,0,0,0,-17},	/* DTAI=36 */
	{2012,7,1,0,0,0,-16},
	{2009,1,1,0,0,0,-15},
	{2006,1,1,0,0,0,-14},	/* BeiDou time reference,  BDT-GPST=-14 */ 
	{1999,1,1,0,0,0,-13},
	{1997,7,1,0,0,0,-12},
	{1996,1,1,0,0,0,-11},
	{1994,7,1,0,0,0,-10},
	{1993,7,1,0,0,0, -9},
	{1992,7,1,0,0,0, -8},
	{1991,1,1,0,0,0, -7},
	{1990,1,1,0,0,0, -6},
	{1988,1,1,0,0,0, -5},
	{1985,7,1,0,0,0, -4},
	{1983,7,1,0,0,0, -3},
	{1982,7,1,0,0,0, -2},
	{1981,7,1,0,0,0, -1}
};

extern void InitPtr(double *a,int num)
{
	for (int i=0;i<num;i++)
	{
		a[i]=0.0;
	}
}
extern void InitPtr(int * a,int num)
{
	for (int i=0;i<num;i++)
	{
		a[i]=0;
	}
}


//UTC_year and month
double FindLeapSec(const int year,const int mon)
{
	int i;
	for ( i=0;i<(int)sizeof(leaps)/(int)sizeof(*leaps);i++) 
	{
		if (leaps[i][0]>year)
		{
			continue;
		}
		else if (leaps[i][0]==year && leaps[i][1]>=mon)
		{
			break;
		}
		else if (leaps[i][1]<year)
		{
			break;
		}
	}
	return leaps[i][6]*1.0;
}

/*
 *sec   utc
 */
double UTC2GPST(const int y,const int mon,double sec)
{
	return sec-FindLeapSec(y,mon);
}

double GPST2UTC(const int y,const int mon,double sec)
{
	return sec+FindLeapSec(y,mon);
}
/*
 *sec gps second
 */
double GPST2BDT(double sec)
{
	return sec-14.0;
}
/*
 * use in eph calc
 */
double timediff(double t1,double t2)
{
	double tc=t1-t2;
	if (tc > 302400.0)
	{
		tc = tc - 604800.0;
	}
	else if (tc < -302400.0)
	{
		tc = tc + 604800.0;
	}
	return tc;
}

double MJD(int y,int mon,int day,int hour,int m,double sec)
{
	double th=hour+m/60.0+sec/3600.0;
	if (mon<=2)
	{
		y-=1;
		mon+=12;
	}
	return floor(365.25*y+1.0e-9)+floor(30.6001*(mon+1)+1.0e-9)+day+th/24+1720981.5-2400000.5;
}

 void WeekSec(int& week,double& sec, YMDHMS currtime,int sysid)
{
	double mjdRefTime=0.0;
	if (sysid==1 || sysid==2||sysid==4)
	{
		mjdRefTime=MJD(gpst0[0],gpst0[1],gpst0[2],gpst0[3],gpst0[4],(double)gpst0[5]);
	}
	else if (sysid==5)
	{
		mjdRefTime=MJD(bdt0[0],bdt0[1],bdt0[2],bdt0[3],bdt0[4],(double)bdt0[5] );
	}
	else if (sysid==3)
	{
		mjdRefTime=MJD(gst0[0],gst0[1],gst0[2],gst0[3],gst0[4],(double)gst0[5]);
	}
	double mjdCurr=MJD(currtime.year,currtime.mon,currtime.day,currtime.hour,currtime.min,currtime.sec);
	double delta_day=floor(mjdCurr-mjdRefTime+1e-8);//round is wrong
	week=(int)floor(delta_day/7.0+1e-8);
	double dow=delta_day-week*7.0;
	sec	=(dow*24.0+currtime.hour)*3600.0+currtime.min*60.0+currtime.sec;
}

 double weightfactor(double ele,int index)
 {
	 if (index==1)
	 {
		 return 1.02/SQ((1.02/sin(ele)+0.02));
	 }
	 if (index==2)
	 {
		 return SQ(sin(ele));
	 }
	 if (index==3)
	 {
		 return (ele*R2D>30.0)?1.0:4*SQ(sin(ele));
	 }
	 if (index==0)
	 {
		 return 1.0;
	 }
 }

extern double cofactor(double ele,int index)
 {
	 if (index==1)
	 {
		 return SQ(1.02/(sin(ele)+0.02));
	 }
	 if (index==2)
	 {
		 return 1.0/SQ(sin(ele));
	 }
	 if (index==3)
	 {
		 return (ele*R2D>30.0)?1.0:1.0/(4.0*SQ(sin(ele)));
	 }
	 if (index==0)
	 {
		 return 1.0;
	 }
 }

extern double geodistcorr(SatePos& sat,double* rec)
{
	return OMEGAE*(sat.sateXYZ[0]*rec[1]-sat.sateXYZ[1]*rec[0])/CLIGHT;
}

extern double geodistcorr(double* sat,double* rec)
{
	return OMEGAE*(sat[0]*rec[1]-sat[1]*rec[0])/CLIGHT;
}

/*
 *satellite clock error ts is obs_t-P/Clight
 */
double ephclkerr(double ts,const BroadEphData eph)
{
	double tc=timediff(ts,eph.toc);
	for (int i=0;i<2;i++)  tc-=eph.satClkBs+eph.satClkDrt*tc+eph.satClkDrtRe*tc*tc;
	return eph.satClkBs+eph.satClkDrt*tc+eph.satClkDrtRe*tc*tc;
}


//t = obs_t-P/C-satclk;  note: t is system time
// GPS BDS GAL
void ephpos(math::matrix<double>& SatePos,BroadEphData eph,double t,double& dts)
{
	// t is obs time
	int		sysid		=Prn2Sysid(eph.prn);
	double Aa			=SQ(eph.sqA);
	double gm		=sysid>2?GM_BDS:GM;
	double n			=sqrt(GM/pow(Aa,3)) + eph.deltan;
	double tk			=timediff(t,eph.toe);
	double Mk		=eph.m0 + n*tk;
	double Ek			=Mk;
	double Ek0		=Ek+1.0;
	while (fabs(Ek-Ek0)>1.0e-12)
	{
		Ek0		=Ek;
		Ek			=Mk+eph.e*sin(Ek0);
	}
	//ek=Ek;
	//double SinEk=sin(Ek),CosEK=cos(EK);
	double vk			=2*atan(sqrt((1.0+eph.e)/(1.0-eph.e))*tan(Ek/2.0));
	double uk			=vk	+eph.omega;
	double detauk	=eph.cuc*cos(2*uk)	+	eph.cus*sin(2*uk);
	double detark	=eph.crc*cos(2*uk)	+	eph.crs*sin(2*uk);
	double detaik	=eph.cic*cos(2*uk)	+	eph.cis*sin(2*uk);

	double u=		uk+detauk;
	double r=			Aa*(1.0-eph.e*cos(Ek))	+	detark;
	double ik=		eph.i0	+	eph.idot*tk+detaik;

	double x=			r*cos(u);
	double y=			r*sin(u);

	double omegaE	=	eph.prn<100?OMEGAE_BDS:OMEGAE;


	if((eph.prn>200 && eph.prn<=205))
	{
		double OMG		=	eph.Omega	+eph.OmegaDot*tk-omegaE*eph.toe;
		math::matrix<double> temp(3,1);
		temp(0,0)			=	cos(OMG)*x	-	y*cos(ik)*sin(OMG);
		temp(1,0)			=	sin(OMG)*x	+	y*cos(ik)*cos(OMG);
		temp(2,0)			=	y*sin(ik);
		//cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(14)<<temp;
		math::matrix<double> Ra(3,3);
		Ra(0,0)	=1.0;
		Ra(1,1)	=COS_5;	Ra(1,2)	=SIN_5;	
		Ra(2,1)	=-SIN_5;	Ra(2,2)	=COS_5;
		double b			=	OMEGAE_BDS*tk;
		math::matrix<double> Rb(3,3);
		Rb(0,0)	=cos(b);	Rb(0,1)	=sin(b);
		Rb(1,0)	=-sin(b);	Rb(1,1)	=cos(b);
		Rb(2,2)	=1.0;
		temp=Rb*Ra*temp;
		SatePos			=temp;

	}
	else if (eph.prn>205||eph.prn<50)
	{
		double OMG			=	eph.Omega + (eph.OmegaDot - omegaE)*tk - omegaE*eph.toe;

		SatePos(0,0)			=	cos(OMG)*x	-	y*cos(ik)*sin(OMG);
		SatePos(1,0)		=	sin(OMG)*x	+	y*cos(ik)*cos(OMG);
		SatePos(2,0)		=	y*sin(ik);
	}
	tk=timediff(t,eph.toc);
	dts= eph.satClkBs+eph.satClkDrt*tk+eph.satClkDrtRe*tk*tk;

	/* relativity correction */
	dts-=2.0*sqrt(gm*Aa)*eph.e*sin(Ek)/SQ(CLIGHT);
	//the TGD correction is out of this part, in SPP
	//the rotation correction needs velocity,and this correction is computed in SPP
}


/**************************** GLONASS       orbits  ********************************************/

/* glonass orbit differential equations --------------------------------------*/
static void deq(const double *x, double *xdot, const double *acc)
{
	double a,b,c,r2=dot(x,x,3),r3=r2*sqrt(r2),omg2=SQ(OMGE_GLO);

	if (r2<=0.0) {
		xdot[0]=xdot[1]=xdot[2]=xdot[3]=xdot[4]=xdot[5]=0.0;
		return;
	}
	a=1.5*J2_GLO*MU_GLO*SQ(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
	b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
	c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
	xdot[0]=x[3]; xdot[1]=x[4]; xdot[2]=x[5];
	xdot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];
	xdot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];
	xdot[5]=(c-2.0*a)*x[2]+acc[2];
}
/* glonass position and velocity by numerical integration --------------------*/
static void glorbit(double t, double *x, const double *acc)
{
	double k1[6],k2[6],k3[6],k4[6],w[6];
	int i;

	deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
	deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
	deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
	deq(w,k4,acc);
	for (i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}
/* glonass ephemeris to satellite clock bias -----------------------------------
* compute satellite clock bias with glonass ephemeris
* 
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
double  geph2clk(double ts, const BroadEphDataGlo geph)
{
	double t;
	int i;

	t=timediff(ts,geph.toe);

	for (i=0;i<2;i++)  t-=(geph.NegativeTauN+geph.PositiveGammaN*t);

	return geph.NegativeTauN+geph.PositiveGammaN*t;
}
/* glonass ephemeris to satellite position and clock bias ----------------------
* compute satellite position and clock bias with glonass ephemeris
* I:
* ts					  obs-P/clight   (gpst)
* geph				  glonass ephemeris (time transferred to GPST)
* O:
* *rs                    satellite position {x,y,z} (ecef) (m)
*dts                    satellite clock bias (s)
*-----------------------------------------------------------------------------*/
extern void  geph2pos(double ts, const BroadEphDataGlo geph, double *rs, double& dts)
{
	double t,tt,x[6];
	int i;

	t=timediff(ts,geph.toe);

	dts=geph.NegativeTauN+geph.PositiveGammaN*t;

	for (i=0;i<3;i++) 
	{
		x[i  ]=geph.Pos[i];
		x[i+3]=geph.Vel[i];
	}
	for (tt=t<0.0?-TSTEP:TSTEP;fabs(t)>1E-9;t-=tt) 
	{
		if (fabs(t)<TSTEP) tt=t;
		glorbit(tt,x,geph.Acc);
	}
	for (i=0;i<3;i++) rs[i]=x[i];

}




/*
 * GPS BDS GAL
 */
void satepos(math::matrix<double>& SatePos,BroadEphData eph,double t,double& dts,math::matrix<double>& SateVel,double& ClkVel)
{
	int sysid=Sysid(eph.prn);
	double t_bar=t;
	if (sysid==5)  	t_bar-=14.0;
	ephpos(SatePos, eph, t_bar, dts);
	//double temp=t_bar+1e-3;
	//double dtsplus=0.0;
	//math::matrix<double> SatePosPlus(3,1);
	//ephpos(SatePosPlus, eph, temp,dtsplus);
	//SateVel=(SatePosPlus-SatePos)*1E3;
	//ClkVel=(dtsplus-dts)*1E3;
}

/********************************SBAS****************************************/
/* sbas ephemeris to satellite clock bias --------------------------------------
* compute satellite clock bias with sbas ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          seph_t *seph     I   sbas ephemeris
* return : satellite clock bias (s)
*-----------------------------------------------------------------------------*/
double seph2clk(double time, const BroadEphDataSBAS seph)
{
	double t;
	int i;

	t=timediff(time,seph.t0);

	for (i=0;i<2;i++) {
		t-=seph.a0+seph.a1*t;
	}
	return seph.a0+seph.a1*t;
}
/* sbas ephemeris to satellite position and clock bias -------------------------
* compute satellite position and clock bias with sbas ephemeris
* args   : gtime_t time     I   time (gpst)
*          seph_t  *seph    I   sbas ephemeris
*          double  *rs      O   satellite position {x,y,z} (ecef) (m)
*          double  *dts     O   satellite clock bias (s)
*          double  *var     O   satellite position and clock variance (m^2)
*-----------------------------------------------------------------------------*/
void seph2pos(double time, const BroadEphDataSBAS seph, double *rs, double& dts)
{
	double t;
	int i;

	t=timediff(time,seph.t0);

	for (i=0;i<3;i++) {
		rs[i]=seph.pos[i]+seph.vel[i]*t+seph.acc[i]*t*t/2.0;
	}
	dts=seph.a0+seph.a1*t;

}
/*************************************************************************/


bool isOdd(int n)
{
	if(n%2==0)
	{
		return false;
	}
	else if((int)fabs((double)(n%2))==1)
	{
		return true;
	}
}
int Fac(int n)
{
	int f=1;
	int s=(n>=0)?n:( (int)fabs((double)n) );
	if(n==0)
		f=1;
	else
	{
		for(int i=s;i>0;i--)
		{
			f	*=i;
		}
	}
	if(n<0 && isOdd(s) )
		f	*=(-1);
	return f;
}




/*
 * n is the order of poly plus 1 
 */
double InterpolNevil(double *x,double* y,int n)
{
	int i,j;
	for (j=1;j<n;j++) {
		for (i=0;i<n-j;i++) {
			y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
		}
	}
	return y[0];
}
/*  
 *  compute the satepos by precise eph
 *  I: 
 *		satindex  the index of sate
 *		sec          the second of transmission time
 *		week		week of transmission time
 *		preceph	pointer to the precise eph
 *		neph
 *	O:
 *		precsatpos satepos 
 */
void PrecEphPos(SatePos& presatpos, int satindex,double sec,int week,PrecEphData* preceph,int neph)
{
	int i,j,k,index;
	double dtime=0.0;
	for (i=0,j=neph-1;i<j;)
	{
		k=(i+j)/2;
		dtime=(preceph[k].week-week)*86400.0*7.0+preceph[k].sec-sec;
		if (dtime<0.0)	i=k+1; else j=k;
	}
	
	index=i<=0?0:i-1;

	/* polynomial interpolation for orbit */
	i=index-(INTERPOLYORDER+1)/2;
	if (i<0) 
	{
		i=0; 
	}
	else if (i+INTERPOLYORDER>=neph)
	{
		i=neph-INTERPOLYORDER-1;
	}

	double tarr[INTERPOLYORDER+1];
	for (j=0;j<=INTERPOLYORDER;j++) tarr[j]=(preceph[j+i].week-week)*86400.0*7.0+preceph[j+i].sec-sec;
	double p[3][INTERPOLYORDER+1];
	double pos[3];
	for (j=0;j<=INTERPOLYORDER;j++)
	{
		for (k=0;k<3;k++) pos[k]=preceph[j+i].satePos[satindex].sateXYZ[k];
		
		double sinrot=sin(OMEGAE*tarr[j]); 
		double cosrot=cos(OMEGAE*tarr[j]);
		p[0][j]= cosrot*pos[0]-sinrot*pos[1];//earth rotation correction
		p[1][j]=sinrot*pos[0]+cosrot*pos[1];
		p[2][j]=pos[2];
	}
	for (i=0;i<3;i++) presatpos.sateXYZ[i]=InterpolNevil(tarr,p[i],INTERPOLYORDER+1);
	
}

double LinearInterp(double x1,double x2,double y1,double y2,double x)
{
	return ((x2-x)*y1+(x-x1)*y2)/(x2-x1);
}

double LagrangeInterp(int step,double* x,double* y,double xk,int flag,int n)
{
	double result=0.0;
	for(int i=0;i<n+1;i++)
	{
		if(xk==x[i])
			return y[i];
	}

	double a=1.0;	
	//Lj(t)=a/b;
	for(int i=0;i<n+1;i++)
		{
			a	*=((xk-x[i])/(x[1]-x[0]));//scale=(x[i]-x[j])
		}

	if(flag==0)
	{
		int b=1;
		b=Fac(step*(-1));
		result	=a/b*y[0];
		for(int i=1;i<n;i++)
		{
			a	*=((xk-x[i-1])/(xk-x[i]));
			b	*=i/(i-1-n);
			result	+=a/b*y[i];
		}
		a	*=((xk-x[n-1])/(xk-x[n]));
		b	*=(-n);
		result	+=a/b*y[n];
	}
	else
	{
		double b=1.0;
		for(int i=0;i<n+1;i++)
		{
			for(int j=0;j<n+1;j++)
			{
				if(i!=j)
				{
					b	*=((x[i]-x[j])/(x[1]-x[0]));
				}
			}
			if(i!=0)
			{
				a	*=(xk-x[i-1])/(xk-x[i]);
				result	+=a/b*y[i];
			}
			else
			{
				result	=a/b*y[0];
			}
		}
		
	}
	return result;
}

extern double Sum(double* s,int num)
{
	for (int i=1;i<num;i++)
	{
		s[0]+=s[i];
	}
	return s[0];
}

double Norm(double* a,int num)
{
	double t=0.0;
	for(int i=0;i<num;i++)	t	+=SQ(a[i]);
	return sqrt(t);
}

double CMax(double* a,int num)
{
	double cmax=a[num-1];
	for(int i=0;i<num-1;i++)cmax=	cmax>=a[i]? cmax: a[i];
	return cmax;
}

double AbsMax(double* a,int num)
{
	double cmax=fabs(a[num-1]);
	for(int i=0;i<num-1;i++)		cmax=	cmax>=fabs(a[i])? cmax: fabs(a[i]);
	return cmax;
}

int AbsIndMaxInd(double* a,int num)
{
	int ind=num-1;
	double cmax=fabs(a[num-1]);
	for(int i=0;i<num-1;i++)
	{
		if (cmax<fabs(a[i]))
		{
			cmax= fabs(a[i]);
			ind=	 i;
		}
	}
	return ind;
}

 double FreqSquareRatio(int sysid,int freqIndex)
{
	if(freqIndex==0)// ionosphere free
	{
		return 0.0;
	}
	if (sysid==1 || sysid==4)
	{
		if (freqIndex==1)
		{
			return 1.0;
		}
		else if(freqIndex==2)
		{
			return pow(FREQ1,2)/pow(FREQ2,2);
		}
		else if(freqIndex==3)
		{
			return pow(FREQ1,2)/pow(FREQ5,2);
		}
	}
	if (sysid==3)
	{
		if (freqIndex==1)
		{
			return 1.0;
		}
		else if(freqIndex==2)
		{
			return pow(FREQ1,2)/pow(FREQ7,2);
		}
		else if(freqIndex==3)
		{
			return pow(FREQ1,2)/pow(FREQ5,2);
		}
	}
	if (sysid==5)
	{
		if (freqIndex==1)
		{
			return 1.0;
		}
		else if(freqIndex==2)
		{
			return pow(FREQ2_BDS,2)/pow(FREQ7_BDS,2);
		}
		else if(freqIndex==3)
		{
			return pow(FREQ2_BDS,2)/pow(FREQ6_BDS,2);
		}
	}
}

 double WaveLength(int sysid,int freqIndex)
 {
	if (sysid==1 || sysid==4)
	{
		if (freqIndex==1)
		{
			return CLIGHT/FREQ1;
		}
		else if(freqIndex==2)
		{
			return CLIGHT/FREQ2;
		}
		else if(freqIndex==3)
		{
			return CLIGHT/FREQ5;
		}
	}
	if (sysid==3)
	{
		if (freqIndex==1)
		{
			return CLIGHT/FREQ1;
		}
		else if(freqIndex==2)
		{
			return CLIGHT/FREQ7;
		}
		else if(freqIndex==3)
		{
			return CLIGHT/FREQ5;
		}
	}
	if (sysid==5)
	{
		if (freqIndex==1)
		{
			return CLIGHT/FREQ2_BDS;
		}
		else if(freqIndex==2)
		{
			return CLIGHT/FREQ7_BDS;
		}
		else if(freqIndex==3)
		{
			return CLIGHT/FREQ6_BDS;
		}
	}
 }


 int	  Prn2Sysid(int prn)
{
	if(prn<=32)
	{
		return 1;
	}
	if(prn<150 &&prn>100)
	{
		return 3;
	}	
	if(prn<250 &&prn>200)
	{
		return 5;
	}	
	if(prn<200 &&prn>150)
	{
		return 4;
	}	
}
 void	  RotationMatrix3(math::matrix<double>& R,double theta, int index)
{
	double		cost	=		cos(theta);
	double		sint	=		sin(theta);
	if (index==1)
	{
		R(0,0)	=	1.0;
		R(1,1)	=	cost;	R(1,2)	=	sint;
		R(2,1)	=	-sint;	R(2,2)	=cost;
	}
	if (index==2)
	{
		R(0,0)	=	cost;	R(0,2)	=-sint;
		R(1,1)	=	1.0;
		R(2,0)	=	sint;	R(2,2)	=cost;
	}
	if (index==3)
	{
		R(2,2)	=	1.0;
		R(0,0)	=	cost;	R(0,1)	=	sint;
		R(1,0)	=	-sint;	R(1,1)	=cost;
	}
}


 CString Sys(int Sysid)
{
	CString s;
	if (Sysid==1)
	{
		s="G";
	}
	if (Sysid==2)
	{
		s="R";
	}
	if (Sysid==3)
	{
		s="E";
	}
	if (Sysid==4)
	{
		s="S";
	}
	if (Sysid==5)
	{
		s="C";
	}
	if (Sysid==6)
	{
		s="M";
	}
	return s;
}
int Sysid(CString sys)
{
	if (sys=='M')
			{
				return 6;
			}
			else if (sys=='G')
			{
				return 1;
			}
			else if (sys=='R')
			{
				return 2;
			}
			else if (sys=='E')
			{
				return 3;
			}
			else if (sys=='S')
			{
				return 4;
			}
			else if (sys=='C')
			{
				return 5;
			}
}

extern int Sysid( int prn )
{
	if (prn<50)
	{
		return 1;
	}
	else if (prn>50&&prn<100)
	{
		return 2;
	}
	else if (prn>100&&prn<150)
	{
		return 3;
	}
	else if (prn>150&&prn<200)
	{
		return 4;
	}
	else if (prn>200&&prn<250)
	{
		return 5;
	}
}

int Prn(CString sys,int prn)
	
{
		 if (sys==_T("G") )return prn;
		else if (sys=='R')return 50+prn;
		else if (sys=='E')return 100+prn;
		else if (sys=='S')return 150+prn;
		else if (sys=='C')return 200+prn;
}
extern void Freq(int sysid,double* freq)
{
	if (sysid==1)
	{
		freq[0]=FREQ1;
		freq[1]=FREQ2;
		freq[2]=FREQ5;
	}
	else if(sysid==5)
	{
		freq[0]=FREQ2_BDS;
		freq[1]=FREQ7_BDS;//B2  B2 < B3
		freq[2]=FREQ6_BDS;//B3
	}
}

/*	
	combine the IF obs for different system
	I:
		sysid
		obs		3*1 (m)
		index	2*1	the frequencies index of two frequency  freq[0]>freq[1]  ==index[0]<index[1]
*/
extern double IonoFree(int sysid,double* obs,int* index)
{
	double freq[3];
	FreqRatio(sysid,freq);
	double sqf1=SQ(freq[index[0]]),sqf2=SQ(freq[index[1]]);
	return  sqf1/(sqf1-sqf2)*obs[index[0]]-sqf2/(sqf1-sqf2)*obs[index[1]];
}


/* form the IF combination of GLONASS */
extern double IonoFreeGlo(double* obs,int* index,int k)
{
	double freq[3];
	FreqRatio(2,freq);
	freq[0]	+=	k*DFRQ1_GLORATIO;
	freq[1]	+=	k*DFRQ2_GLORATIO;
	double sqf1=SQ(freq[index[0]]),sqf2=SQ(freq[index[1]]);
	return  sqf1/(sqf1-sqf2)*obs[index[0]]-sqf2/(sqf1-sqf2)*obs[index[1]];
}
/* the frequency, index=0,1,2*/
extern double Freq(int prn,int ind)
{
	if(prn<50)
	{
		if (ind==0)return FREQ1;
		else if(ind==1)return FREQ2;
		else if(ind==2)return FREQ5;
	}
	else if(prn>50 && prn<100)
	{
		if (ind==0)return FREQ1_GLO;
		else if(ind==1)return FREQ2_GLO;
		else if(ind==2)return FREQ3_GLO;
	}
	else if(prn>200)
	{
		if (ind==0)return FREQ2_BDS;
		else if(ind==1)return FREQ7_BDS;
		else if(ind==2)return FREQ6_BDS;
	}

}

/* the frequency of system, index=0,1,2*/
double FreqSys(int sysid,int ind)
{
	if(sysid==1)
	{
		if (ind==0)return FREQ1;
		else if(ind==1)return FREQ2;
		else if(ind==2)return FREQ5;
	}
	else if(sysid==5)
	{
		if (ind==0)return FREQ2_BDS;
		else if(ind==1)return FREQ7_BDS;
		else if(ind==2)return FREQ6_BDS;
	}
}

extern double FreqSysGlo(int ind,int n)
{
	double freq=0.0;
	if (ind==0)freq=FREQ1_GLO+n*DFRQ1_GLO;
	if (ind==1)freq=FREQ2_GLO+n*DFRQ2_GLO;
	return freq;
}

extern int doy( int Year,int Month,int Day )
{
	int DayNums[12];
	DayNums[0]=31;
	DayNums[2]=31;
	DayNums[4]=31;
	DayNums[6]=31;
	DayNums[7]=31;
	DayNums[9]=31;
	DayNums[11]=31;
	DayNums[3]=30;
	DayNums[5]=30;
	DayNums[8]=30;
	DayNums[10]=30;

	if ((Year%400==0) || (Year%4==0 && Year%100!=0))
		DayNums[1]=29;
	else
		DayNums[1]=28;
	int DayofYear=0;
	for (int i=0;i<Month-1;i++)
		DayofYear+=DayNums[i];
	DayofYear+=Day;
	return DayofYear;
}

extern int doy(int sysid, int gpsw,double gpss )
{
	double mjd	= gpsw*7.0+44244+gpss/3600.0/24.0;
	double		a		=floor(mjd+0.5+0.5+1E-10)+2400000.0;
	double		frac	=mjd+0.5+2400000.5-a;
	double		b		=a+1537;
	int			c		=(int)((b-122.1)/365.25+1E-10);
	int			d		=(int)(365.25*c+1E-16);
	int			e		=(int)((b-d)/36.6001+1E-10);
	int			day	=(int)(b-d-floor(30.6001*e)+1E-10);
	int			mon	=e-1-12*((int)(e/14.0+1E-10));
	int			year	=c-4715-(int)((7+mon)/10.0+1E-10);
	//int			hour=(int)(frac*24.0+1E-10);
	//int			mins	=(int)((frac*24.0-hour)*60+1E-10);
	//double		sec	=(frac*24.0-hour)*60-mins;

	return doy(year,mon,day);
}

extern double CombFreq( int sysid,int*coef )
{
	double s=0;
	for (int i=0;i<3;i++) s+=coef[i]*FreqSys(sysid,i);
	return s;
}

/*
 *			回首向来萧瑟处，也无风雨也无晴
 *			
 **/
extern void FreqRatio(int sysid,double* freq)
{
	if (sysid==1)
	{
		freq[0]=FREQ1RATIO;
		freq[1]=FREQ2RATIO;
		freq[2]=FREQ5RATIO;
	}
	else if(sysid==2)
	{
		freq[0]=FREQ1_GLORATIO;
		freq[1]=FREQ2_GLORATIO;
		freq[2]=FREQ3_GLORATIO;
	}
	else if(sysid==5)
	{
		freq[0]=FREQ2_BDSRATIO;
		freq[1]=FREQ7_BDSRATIO;//B2  B2 < B3
		freq[2]=FREQ6_BDSRATIO;//B3
	}
}
extern double CombObs(int sysid,int* coef,double* threeobs)
{
	double freq[3];
	FreqRatio(sysid,freq);
	double a=0.0;
	for(int i=0;i<3;i++)
	{
		a+=coef[i]*freq[i]*threeobs[i];
	}
	return a/(CombFreq(sysid,coef)/1E9);
}

extern void  EarthConst(int ellipInd,double& a,double& eE)
{

	if (ellipInd==1 || ellipInd==4)	
	{
		a		=	RE_WGS84;
		eE		=	ee;
	}
	else if(ellipInd==2)
	{
		a=RE_GLO;
		eE=ee_GLO;
	}
	else if (ellipInd==3)
	{
		a=RE_GAL;
		eE=ee_GAL;
	}
	else if(ellipInd==5)	
	{
		a		=	RE_CGCS2000;
		eE		=	ee_BDS;
	}
	else if(ellipInd==6)	
	{
		a		=	RE_BJ54;
		eE		=	ee_BJ54;
	}
	else if(ellipInd==7)	
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

extern void  BLH2XYZ(double* BLH,int ellipInd,double* XYZ)
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
			ellipInd same as sysid
		O:
			BLH
*/
extern void  XYZ2BLH(double* XYZ,int ellipInd,double* BLH)
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

extern void xyz2plh(double *xyz,double* plh)
{
	double r=Norm(xyz,2),b=(1-FE_WGS84)*RE_WGS84;
	double u=atan2(xyz[2]*RE_WGS84,r*b);

	double phi=atan2(xyz[2]+(ee/(1.0-ee)*b)*pow(sin(u),3),r-(ee*RE_WGS84)*pow(cos(u),3));
	double N=RE_WGS84/sqrt(1.0-ee*pow(sin(phi),2));
	plh[0]=phi;plh[1]=atan2(xyz[1],xyz[1]);plh[2]=r/cos(phi)-N;
}

/* in and out   vector*/
extern math::matrix<double> XYZ2BLH(math::matrix<double> XYZ,int ellipInd)
{
	double xyz[3],blh[3];
	for (int i=0;i<3;i++) xyz[i]=XYZ(i,0);
	XYZ2BLH(xyz,ellipInd,blh);
	math::matrix<double>temp(3,1);
	for (int i=0;i<3;i++) temp(i,0)=blh[i];
	return temp;
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
extern void  NEU2XYZ(double* NEU,double B0,double L0,double* XYZ)
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
extern void  XYZ2NEU(double* XYZ,double B0,double L0,double* NEU)
{
	NEU[0]		=	-sin(B0)*cos(L0)*XYZ[0]  -sin(B0)*sin(L0)*XYZ[1] + cos(B0)*XYZ[2];
	NEU[1]		=	-sin(L0)*XYZ[0]         +  cos(L0)*XYZ[1];
	NEU[2]		=	cos(B0)*cos(L0)*XYZ[0]   +cos(B0)*sin(L0)*XYZ[1]  +sin(B0)*XYZ[2];
}
extern math::matrix<double>  XYZ2NEU(math::matrix<double> XYZ,double B0,double L0)
{
	double xyz[3],neu[3];
	for (int i=0;i<3;i++) xyz[i]=XYZ(i,0);
	XYZ2NEU(xyz,B0,L0,neu);
	math::matrix<double>temp(3,1);
	for (int i=0;i<3;i++) temp(i,0)=neu[i];
	return temp;
}

extern math::matrix<double>  XYZ2NEU(double* coorBase,double* coorRover )
{
	double blh[3];
	math::matrix<double>dxyz(3,1);
	XYZ2BLH(coorBase,1,blh);
	for(int i=0;i<3;i++) dxyz(i,0)=coorRover[i]-coorBase[i];
	return XYZ2NEU(dxyz,blh[0],blh[1]);
}

//extern math::matrix<double>  XYZ2NEU(math::matrix<double> coorBase,math::matrix<double> coorRover )
//{
//	double blh[3],cooB[3];
//	math::matrix<double>dxyz(3,1);
//	for(int i=0;i<3;i++) 
//	{
//		dxyz(i,0)=coorRover(i,0)-coorBase(i,0);
//		cooB[i]=coorBase(i,0);
//	}
//	XYZ2BLH(coorBa,1,blh);
//	return XYZ2NEU(dxyz,blh[0],blh[1]);
//}
/*
 *XYZ the topocenter
 *XYZs the object coordinate
 */
void XYZ2RAH(double* XYZ,int sysid,double* XYZs,double& ele,double& azi)
{
	double sk[3];
	double sk1[3];
	
	XYZ2BLH(XYZ,sysid,sk);//sk=BLH
	for (int i=0;i<3;i++)
	{
		sk1[i]=XYZs[i]-XYZ[i];
	}
	double sk2[3];//
	XYZ2NEU(sk1,sk[0],sk[1],sk2);//sk2=NEU sk1=diff coord
	
	double rah[3];
	NEU2RAH(sk2,rah);
	
	ele=rah[2];
	azi=rah[1];
}
/*
	get azimuth
		I:
			dn	difference if north component between two points
			de		difference if east component between two points
		return:
			azimuth(rad)
*/
extern double  getAzi(double dn,double de)
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
extern void	 NEU2RAH(double *NEU,double* RAH)
{
	RAH[0]		=	Norm(NEU,3);
	RAH[1]		=	getAzi(NEU[0],NEU[1]);
	RAH[2]		= atan( NEU[2]/Norm(NEU,2) );
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
extern void  RAH2NEU(double* RAH,double* NEU)
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
extern void  GaussProj(double B,double L,double L0,int ellipInd,double* xy)
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

/*
 *Kronecker C=A@B, @is the symbol of kronecker
 *I:
 *	A			m*n
 *	rowA    m
 *	colA		n
 *	B			k*l
 *	rowB		k
 *	colB		l
 *	flag		indicates the process of kroneker, for diagonal matrix process is simplified 
 *	1=A and B are both diagonal  and square  (default)
 *	2=B is diagonal , A is not,but  A is syms
 *	3=A	is diagonal , B is not
 *	4=A and B are not diagonal
 *O:
 *	C	mk*nl
 *	Note:
 *		A and B are square and symmetrical if flag=1,2,3
 *		due to the list of obs is ordered by  frequency，the combination of vc-mat is Qyy@I, Qyy is k*k (k=1,2,3), I is n*n  
 */
extern math::matrix<double> Kronecker( math::matrix<double>A,math::matrix<double> B,int flag)
{
	int rowA=A.RowNo(), colA=A.ColNo(), rowB=B.RowNo(),colB=B.ColNo();
	math::matrix<double>C(rowA*rowB,colA*colB);
	if (flag==1)
	{
		// A and B are diagonal and square
		for (int i=0;i<rowA;i++)
		{
			for (int j=0;j<rowB;j++)
			{
				C(rowB*i+j,rowB*i+j)	=A(i,i);
			}
		}
	}
	else if (flag==2)
	{
		for (int i=0;i<rowA;i++)
		{
			for (int j=i;j<rowA;j++)
			{
				if (A(i,j)==0.0)
				{
					continue;
				}
				else
				{
					for (int k=0;k<rowB;k++)
					{
						if (i==j)
						{
							C(rowB*i+k,rowB*i+k)	=A(i,j)*B(k,k);
						} 
						else
						{
							C(rowB*i+k,rowB*j+k)	=A(i,j)*B(k,k);
							C(rowB*j+k,rowB*i+k)	=C(rowB*i+k,rowB*j+k);//;A(j,i)*B(k,k);
						}
					}
				}
			}
		}
	}
	else if (flag==3)
	{
		for (int i=0;i<rowA;i++)
		{
			if (A(i,i)==0.0)
			{
				continue;			
			} 
			else
			{
				for (int j=0;j<rowB;j++)
				{
					for (int k=j;k<rowB;k++)
					{
						C(i*rowB+j,i*rowB+k)	=B(j,k);
						C(i*rowB+k,i*rowB+j)	=B(j,k);
					}
				}
			}
			
		}
	}
	else if (flag==4)
	{
		for (int i=0;i<rowA;i++)
		{
			for (int j=0;j<colA;j++)
			{
				if (A(i,j)==0.0)
				{
					continue;
				} 
				else
				{
					for (int p=0;p<rowB;p++)
					{
						for (int q=0;q<colB;q++)
						{
							C(i*rowB+p,j*colB+q)	=A(i,j)*B(p,q);
						}
					}
				}
			}
		}
	}
	return C;
}


/*
 *form the vc matrix of combination obs
 *I:
 *	coef k*k	the coef of combination   k=1,2,3
 *	Qyy	the cofactor or vc mat of non-combined obs 
 *	flag	the number of obs = k
 *O:
 *	Qyy		the vc mat of  combination 	
 */
extern void VcmatCom(math::matrix<double> coef,math::matrix<double>& Qyy,int flag)
{
	math::matrix<double>temp(flag,flag);
	temp=Qyy;
	Qyy	=coef*temp*(~coef);
}
extern double traceMat(math::matrix<double> A)
{
	int num=min(A.RowNo(),A.ColNo());
	double s=0.0;
	for(int i=0;i<num;i++) s+=A(i,i);
	return s;
}
/*
 *eliminate the k th row and k th col of A_n*n to A_n-1*n-1 and  k is in [1, n]
 */
extern math::matrix<double>ElimRowCol(math::matrix<double>A,math::matrix<double>&U,int k)
{
	int rankA=A.RowNo(),cnt=0;
	math::matrix<double>B(rankA-1,rankA-1);
	double a=A(k-1,k-1),u=U(k-1,0);
	math::matrix<double> vec(1,rankA-1),U_bar(rankA-1,1);
	for (int i=0;i<rankA;i++)
	{
		if (i!=k-1)
		{
			U_bar(cnt,0)=U(i,0);
			vec(0,cnt++)=A(k-1,i);
		}
	}
	//U=U_bar-~vec/a*u;// for normal equation
	U=U_bar;
	int cnt1=0;
	cnt=0;
	for (int i=0;i<rankA;i++)
	{
		if (i!=k-1) 
		{
			cnt1=cnt;
			for (int j=i;j<rankA;j++)
			{
				if (j==k-1) continue;
				if (j!=k-1)
				{
					B(cnt,cnt1)=A(i,j);
					B(cnt1,cnt)=B(cnt,cnt1);
					cnt1++;
				}
			}
			cnt++;
		}
	}
	//return B-~vec*vec/a;// for normal equation
	return B;
}

extern math::matrix<double>ElimRowColNEQ(math::matrix<double>A,math::matrix<double>&U,int k)
{
	int rankA=A.RowNo(),cnt=0;
	math::matrix<double>B(rankA-1,rankA-1);
	double a=A(k-1,k-1),u=U(k-1,0);
	math::matrix<double> vec(1,rankA-1),U_bar(rankA-1,1);
	for (int i=0;i<rankA;i++)
	{
		if (i!=k-1)
		{
			U_bar(cnt,0)=U(i,0);
			vec(0,cnt++)=A(k-1,i);
		}
	}
	U=U_bar-~vec/a*u;// for normal equation
	//U=U_bar;
	int cnt1=0;
	cnt=0;
	for (int i=0;i<rankA;i++)
	{
		if (i!=k-1) 
		{
			cnt1=cnt;
			for (int j=i;j<rankA;j++)
			{
				if (j==k-1) continue;
				if (j!=k-1)
				{
					B(cnt,cnt1)=A(i,j);
					B(cnt1,cnt)=B(cnt,cnt1);
					cnt1++;
				}
			}
			cnt++;
		}
	}
	return B-~vec*vec/a;// for normal equation
	//return B;
}
/*
 * A and B must be the symmetrical
 */
extern math::matrix<double> DiagMatSym(math::matrix<double>A,math::matrix<double>B)
{
	int rowA=A.RowNo(), rowB=B.RowNo();
	math::matrix<double>temp(rowA+rowB,rowA+rowB);
	for (int i=0;i<rowA;i++)
	{
		for (int j=i;j<rowA;j++)
		{
			temp(i,j)	=A(i,j);
			temp(j,i)	=A(i,j);
		}
	}

	for (int i=0;i<rowB;i++)
	{
		for (int j=0;j<rowB;j++)
		{
			temp(rowA+i,rowA+j)	=B(i,j);
			temp(rowA+j,rowA+i)	=B(i,j);
		}
	}
	return temp;
}
extern math::matrix<double> DiagMatSym(math::matrix<double>A,int rowA,int rowB)
{
	math::matrix<double>temp(rowA+rowB,rowA+rowB);
	for (int i=0;i<rowA;i++)
	{
		for (int j=i;j<rowA;j++)
		{
			temp(i,j)	=A(i,j);
			temp(j,i)	=A(i,j);
		}
	}

	for (int i=0;i<rowB;i++)
	{
		for (int j=i;j<rowB;j++)
		{
			temp(rowA+i,rowA+j)	=0.0;
			temp(rowA+j,rowA+i)	=0.0;
		}
	}
	return temp;
}
/*
 * A and B are not symetrical
 */
extern math::matrix<double> DiagMat(math::matrix<double>A,int rowA,int colA,math::matrix<double>B,int rowB,int colB)
{
	math::matrix<double>temp(rowA+rowB,colA+colB);
	for (int i=0;i<rowA;i++)
	{
		for (int j=0;j<colA;j++)
		{
			temp(i,j)	=A(i,j);
		}
	}

	for (int i=0;i<rowB;i++)
	{
		for (int j=0;j<colB;j++)
		{
			temp(rowA+i,colA+j)	=B(i,j);
		}
	}
	return temp;
}
/* flag=1 A is sym ,flag=2 B is any part  B=A(rowFirst:rowLast,colFirst:colLast)*/
extern math::matrix<double>GetBlockMat(math::matrix<double>A,int rowFirst,int rowLast,int colFirst,int colLast,int flag)
{
	int rowB,colB;
	rowB=rowLast-rowFirst+1;
	colB=colLast-colFirst+1;
	math::matrix<double>B(rowB,colB);
	int i,j;
	if (flag==1)
	{//rowB=colB
		for (i=0;i<rowB;i++)
		{
			B(i,i)=A(rowFirst-1+i,rowFirst-1+i);
			for (j=i+1;j<rowB;j++)
			{
				B(i,j)=A(rowFirst-1+i,colFirst-1+j);
				B(j,i)=A(rowFirst-1+i,colFirst-1+j);
			}
		}
	}
	if (flag==2)
	{
		for (i=0;i<rowB;i++)
		{
			for (j=0;j<colB;j++)
			{
				B(i,j)=A(rowFirst-1+i,colFirst-1+j);
			}
		}
	}
	return B;
}
/* put B in A , then B=A(rowFirst:rowLast,colFirst:colLast) ,flag=1, sym put,flag=2 any    */
extern void PutMat(math::matrix<double>&A,math::matrix<double>B,int rowFirst,int colFirst,int flag)
{
	int i,j;
	int rowB=B.RowNo(),colB=B.ColNo();
	if (flag==1)
	{	
		for (i=0;i<rowB;i++)
		{
			for (j=0;j<colB;j++)
			{
				A(rowFirst-1+i,colFirst-1+j)=B(i,j);
				A(colFirst-1+j,rowFirst-1+i)=B(i,j);
			}
		}
	}
	if (flag==2)
	{	
		for (i=0;i<rowB;i++)
		{
			for (j=0;j<colB;j++)
			{
				A(rowFirst-1+i,colFirst-1+j)=B(i,j);
			}
		}
	}
}

/* return [A;B] */
extern math::matrix<double> VecMat(int col,math::matrix<double>A,math::matrix<double>B)
{
	int rowA=A.RowNo(), rowB=B.RowNo();
	math::matrix<double>temp(rowA+rowB,col);
	int i,j;

		for (j=0;j<col;j++)
		{
			for (i=0;i<rowA;i++)
			{
				temp(i,j)=A(i,j);
			}
			for (i=0;i<rowB;i++)
			{
				temp(i+rowA,j)=B(i,j);
			}
		}
	return temp;
}
/* return [A;0] */
extern math::matrix<double> VecMat(int col,math::matrix<double>A,int rowB)
{
	int rowA=A.RowNo();
	math::matrix<double>temp(rowA+rowB,col);
	int i,j;

	for (j=0;j<col;j++)
	{
		for (i=0;i<rowA;i++)
		{
			temp(i,j)=A(i,j);
		}
		for (i=0;i<rowB;i++)
		{
			temp(i+rowA,j)=0.0;
		}
	}
	return temp;
}

/*return [A B] */
extern math::matrix<double> ConvergeMat(int row,math::matrix<double>A,math::matrix<double>B)
{
	int colA=A.ColNo(),colB=B.ColNo();
	math::matrix<double>temp(row,colA+colB);
	for (int i=0;i<row;i++)
	{
		for (int j=0;j<colA;j++)
		{
			temp(i,j)	=A(i,j);
		}
		for (int k=0;k<colB;k++)
		{
			temp(i,colA+k)	=B(i,k);
		}
	}
	return temp;
}

/*return [A B], B=0*/
extern math::matrix<double> ConvergeMat(int row,math::matrix<double>A,int colA,int colB)
{
	math::matrix<double>temp(row,colA+colB);
	for (int i=0;i<row;i++)
	{
		for (int j=0;j<colA;j++) temp(i,j)	=A(i,j);
		for (int k=0;k<colB;k++) temp(i,colA+k)	=0.0;
	}
	return temp;
}

/* return A*A' =B, B is row * row */
extern  math::matrix<double>MultSelfTrans(math::matrix<double> A,int row,int col)
{
	math::matrix<double> B(row,row);
	int i,j,k;
	
	for (i=0;i<row;i++)
	{
		for (j=i;j<row;j++)
		{
			B(i,j)=0.0;
			for (k=0;k<col;k++)
			{
				B(i,j)+=A(i,k)*A(j,k);
			}
			B(j,i)=B(i,j);
		}
	}
	return B;
}

/* return A'*A =B, B is col*col */
extern  math::matrix<double>MultSelfTrans2(math::matrix<double> A,int row,int col)
{
	math::matrix<double> B(col,col);
	int i,j,k;

	for (i=0;i<col;i++)
	{
		for (j=i;j<col;j++)
		{
			B(i,j)=0.0;
			for (k=0;k<row;k++)
			{
				B(i,j)+=A(k,i)*A(k,j);
			}
			B(j,i)=B(i,j);
		}
	}
	return B;
}

/* move vec ai to the end column
 *  [a1,...,ai,...,an]  ->[a1,...,ai+1,...,an,ai] 
 */
extern  math::matrix<double>MoveVecColEnd(math::matrix<double>A,int i)
{
	int colA=A.ColNo();
	if (colA==i)
	{
		return A;
	}
	math::matrix<double> B(colA,colA);

	int j;
	for (j=0;j<i-1;j++)
	{
		B(j,j)=1.0;
	}
	for (j=i;j<colA;j++)
	{
		B(j,j-1)=1.0;
	}
	B(i-1,colA-1)=1.0;

	return A*B;
}
extern  math::matrix<double>MoveVecRowEnd(math::matrix<double>A,int i)
{
	int rowA=A.RowNo();
	if (rowA==i)
	{
		return A;
	}
	math::matrix<double> B(rowA,rowA);

	int j;
	for (j=0;j<i-1;j++)
	{
		B(j,j)=1.0;
	}
	for (j=i;j<rowA;j++)
	{
		B(j-1,j)=1.0;
	}
	B(rowA-1,i-1)=1.0;
	return B*A;
}
extern math::matrix<double>MoveVecRowColEnd(math::matrix<double>A,int i)
{
	int rowA=A.RowNo();
	int colA=A.ColNo();
	math::matrix<double>C(rowA,colA);
	C=MoveVecRowEnd(A,i);
	return MoveVecColEnd(C,i);
}

/* insert Zero matrix in i th col 
 *		A11	A12	-> A11	 0	 A12	
 */
math::matrix<double> InsertZeroCol(math::matrix<double>A,int colA11,int colNumZero)
{
	if (colNumZero==0)
	{
		return A;
	}
	math::matrix<double>B(A.ColNo(),A.ColNo()+colNumZero);
	int i;
	for (i=0;i<colA11;i++)
	{
		B(i,i)=1.0;
	}
	for (i=0;i<A.ColNo()-colA11;i++)
	{
		B(colA11+i,colA11+colNumZero+i)=1.0;
	}

	return A*B;
}

/* insert Zero matrix in i th row 
 *		A11		-> A11	 	 
 *		A21	 		  0			   
 *						A21	 	 
 */
math::matrix<double> InsertZeroRow(math::matrix<double>A,int rowA11,int rowNumZero)
{
	if (rowNumZero==0)
	{
		return A;
	}
	math::matrix<double>B(A.RowNo()+rowNumZero,A.RowNo());
	int i;
	for (i=0;i<rowA11;i++)
	{
		B(i,i)=1.0;
	}
	for (i=0;i<A.RowNo()-rowA11;i++)
	{
		B(rowA11+rowNumZero+i,rowA11+i)=1.0;
	}
	return B*A;
}

/* insert Zero matrix in i th row and i th col
 *		A11	A12	-> A11	 0	 A12	
 *		A21	A22  		  0		 0	   0 
 *							A11	 0	 A12
 */
extern math::matrix<double> InsertZeroRowCol(math::matrix<double>A,int rowA11,int rowcolNumZero)
{
	if (rowcolNumZero==0)
	{
		return A;
	}
	math::matrix<double>B(A.RowNo()+rowcolNumZero,A.RowNo());
	//B=InsertZeroCol(A,colA11,rowcolNumZero);
	//return InsertZeroRow(B,rowA11,rowcolNumZero);
	int i;
	for (i=0;i<rowA11;i++)	B(i,i)=1.0;
	for (i=0;i<A.RowNo()-rowA11;i++)	B(rowA11+rowcolNumZero+i,rowA11+i)=1.0;
	return B*A*(~B);
}

/*
 *insert Zero row and col, and -make the A_n*n toA_n+1*n+1, and the k th is ZERO vector in A_n+1*n+1
 or insert before the [A]k.
 */
extern math::matrix<double> InsertZeroRowCol(math::matrix<double>A,int k)
{
	int rankA=A.RowNo(),i,j;
	math::matrix<double>B(rankA+1,rankA+1);
	for (i=0;i<k-1;i++)
	{
		for (j=i;j<k-1;j++)
		{
			B(i,j)=A(i,j);
			B(j,i)=A(i,j);
		}
	}
	for (i=0;i<k-1;i++)
	{
		for (j=k-1;j<rankA;j++)
		{
			B(i,j+1)=A(i,j);
			B(j+1,i)=A(i,j);
		}
	}
	for (i=k-1;i<rankA;i++)
	{
		for (j=i;j<rankA;j++)
		{
			B(i+1,j+1)=A(i,j);
			B(j+1,i+1)=A(i,j);
		}
	}
	return B;
}
/*change the i th and j th row of A*/
extern void ChangeRow(math::matrix<double>&A,int i,int j)
{
	int k;int col=A.ColNo();
	double rowi;
	for (k=0;k<col;k++)
	{
		rowi=A(i-1,k);
		A(i-1,k)=A(j-1,k);
		A(j-1,k)=rowi;
	}
}

/*change the i th and j th col of A*/
extern void ChangeCol(math::matrix<double>&A,int i,int j)
{
	int k;int row=A.RowNo();
	double coli;
	for (k=0;k<row;k++)
	{
		coli=A(k,i-1);
		A(k,i-1)=A(k,j-1);
		A(k,j-1)=coli;
	}
}

/*change the i th and j th row ,then the col of A*/
extern void ChangeRowCol(math::matrix<double>&A,int i,int j)
{
	ChangeRow(A,i,j);
	ChangeCol(A,i,j);
}

extern double	Mean(double *a,int num)
{
	double b=0.0;
	for (int i=0;i<num;i++)
	{
		b	=b*((double)i/((double)i+1))+a[i]/(i+1);//b*i/(i+1)+a[i]/(i+1);// false    
	}
	return	b;
}
extern double	Std(double* a,int num)
{
	double b=0.0;
	double a_bar=Mean(a,num);
	for (int i=0;i<num;i++)
	{
		b	+=SQ(a[i]-a_bar)/(num-1);
	}
	return sqrt(b);
}

extern double Rms(double*a,double t,int num )
{
	double	b=0.0;
	for (int i=0;i<num;i++)
	{
		b	+=SQ(a[i]-t)/(num);
	}
	return sqrt(b);
}

//retrun a lower triangular matrix A=L*L'  A=n*n
extern math::matrix<double> Cholesky( math::matrix<double>A,int n )
{
	math::matrix<double>L(n,n);
	double LL=0.0;
	double LR	=0.0;
	for (int j=0;j<n;j++)
	{
		LL=0.0;
			for (int r=0;r<j;r++)
			{
				LL	+=SQ(L(j,r));
			}
		
		L(j,j)=sqrt(A(j,j)-LL);
			for (int i=j;i<n;i++)
			{
				LR=0.0;
				for (int r=0;r<j;r++)
				{
					LR+=L(i,r)*L(j,r);
				}
				L(i,j)=(A(i,j)-LR)/L(j,j);
			}
	}
	return L;
}

//inverse of  L
extern math::matrix<double> InvLowTri( math::matrix<double>L,int n )
{
	math::matrix<double>invL(n,n);
	for (int i=0;i<n;i++)
	{
		invL(i,i)=1.0/L(i,i);
		for (int j=0;j<i;j++)
		{
			for (int r=j;r<i;r++)
			{
				invL(i,j)	-=L(i,r)*invL(r,j)*invL(i,i);
			}	
		}
	}
	return invL;
}
//A*A'  A is lower triangular matrix 
extern math::matrix<double> MultiplyselfLowerUpper( math::matrix<double>A,int n )
{
	math::matrix<double>B(n,n);
	for (int i=0;i<n;i++)
	{
		for (int j=i;j<n;j++)
		{
			for (int k=0;k<i+1;k++)
			{
				B(i,j)+=A(i,k)*A(j,k);
			}
			B(j,i)=B(i,j);
			
		}
	}
	return B;
}

//A*B  A is lower triangular matrix n x m, B is m*p
extern math::matrix<double> MultiplyLowerOther( math::matrix<double>A,math::matrix<double>B,int n,int p )
{
	math::matrix<double>C(n,p);
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<p;j++)
		{
			for (int k=0;k<i+1;k++)
			{
				C(i,j)+=A(i,k)*B(k,j);
			}
		}
	}
	return C;
}

//A*B A m x n,  B is lower triangular matrix  n x n
extern math::matrix<double> MultiplyOtherLower( math::matrix<double>A,math::matrix<double>B,int m,int n )
{
	math::matrix<double>C(m,n);
	for (int i=0;i<m;i++)
	{
		for (int j=0;j<n;j++)
		{
			for (int k=j;k<n;k++)
			{
				C(i,j)+=A(i,k)*B(k,j);
			}
		}
	}
	return C;
}

//A'*A  A is lower triangular matrix 
extern math::matrix<double> MultiplyselfUpperLower( math::matrix<double>A,int n )
{
	math::matrix<double>B(n,n);
	for (int i=0;i<n;i++)
	{
		for (int j=i;j<n;j++)
		{
			for (int k=(i>=j)?i:j;	k<n;k++)
			{
				B(i,j)+=A(k,i)*A(k,j);
			}
			B(j,i)=B(i,j);

		}
	}
	return B;
}

//solve the equation Ax=U  A: n X n and lower triangular matrix
extern math::matrix<double>SolveLowerEquation(math::matrix<double>A,math::matrix<double>U,int n)
{
	math::matrix<double>x(n,1);
	double s=0.0;
	for (int i=0;i<n;i++)
	{
		s=0.0;
		for (int j=0;j<i;j++)
		{
			s+=A(i,j)*x(j,0);
		}
		x(i,0)=(U(i,0)-s)/A(i,i);
	}
	return x;
}

//solve the equation Ax=U  A: n X n and upper triangular matrix
extern math::matrix<double>SolveUpperEquation(math::matrix<double>A,math::matrix<double>U,int n)
{
	math::matrix<double>x(n,1);
	double s=0.0;
	for (int i=n;i>0;i--)
	{
		s=0.0;
		for (int j=n;j>i;j--)
		{
			s+=A(i-1,j-1)*x(j-1,0);
		}
		x(i-1,0)=(U(i-1,0)-s)/A(i-1,i-1);
	}
	return x;
}


//solve the equation A'x=U  A: n X n and lower triangular matrix
extern math::matrix<double>SolveUpperEquation2(math::matrix<double>A,math::matrix<double>U,int n)
{
	math::matrix<double>x(n,1);
	double s=0.0;
	for (int i=n;i>0;i--)
	{
		s=0.0;
		for (int j=n;j>i;j--)
		{
			s+=A(j-1,i-1)*x(j-1,0);
		}
		x(i-1,0)=(U(i-1,0)-s)/A(i-1,i-1);
	}
	return x;
}
extern math::matrix<double>CholeskyInv(math::matrix<double>A,int n)
{
	math::matrix<double>L;
	L=Cholesky(A,n);
	L=InvLowTri(L,n);
	return MultiplyselfUpperLower(L,n);
}
extern math::matrix<double>CholeskyInv(math::matrix<double>A)
{
	int n=A.RowNo();
	math::matrix<double>L;
	L=Cholesky(A,n);
	L=InvLowTri(L,n);
	return MultiplyselfUpperLower(L,n);
}
//solve  the equation a*c=b    a: n X n; b:n X m
extern bool MatrixSovle( math::matrix<double> a,math::matrix<double> b,math::matrix<double> &c, int n, int m)
{
	double at;
	double bt;
	double am;
	int i,j,k;
	int tt;
	int N=n+m;
	std::vector< std::vector<double > > p(n);	//p: n X (n+m)
	for(i=0;i<n;i++)
		p[i].resize(N);

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			p[i][j]=a(i,j);
		}
		for(j=n;j<N;j++)
		{
			p[i][j]=b(i,j-n);
		}
	}
	for(k=0;k<n;k++)			//求解
	{
		at=fabs(p[k][k]);
		tt=k;
		for(j=k+1;j<n;j++)
		{
			bt=fabs(p[j][k]);
			if(at<bt)
			{
				at=bt;
				tt=j;
			}
		}
		if(tt!=k)
			for(j=k;j<N;j++)
			{
				am=p[k][j];
				p[k][j]=p[tt][j];
				p[tt][j]=am;
			}
			if(at<0.00000001)
			{
				return false;
			}
			am=1.0/p[k][k];
			for(j=k;j<N;j++)
			{
				p[k][j]=p[k][j]*am;
			}

			for(i=0;i<n;i++)
			{
				if(k!=i)
				{
					am=p[i][k];
					for(j=0;j<N;j++)
					{
						p[i][j]=p[i][j]-p[k][j]*am;
					}
				}
			}
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			c(i,j)=p[i][j+n];
		}
	}
		
		return true;
}

/* remove the row and col i th*/
extern math::matrix<double>RemoveRowCol(math::matrix<double>A,int i)
{
	int row=A.RowNo();
	math::matrix<double>B(row-1,row);
	int k,j;
	for (k=0;k<i-1;k++)
	{
		B(k,k)=1.0;
	}
	for (k=0;k<row-i;k++)
	{
		B(i+k-1,i+k)=1.0;
	}
	return B*A*(~B);
}

/* remove the row i th*/
extern math::matrix<double>RemoveRow(math::matrix<double>A,int i)
{
	int row=A.RowNo();
	math::matrix<double>B(row-1,row);
	int k,j;
	for (k=0;k<i-1;k++)
	{
		B(k,k)=1.0;
	}
	for (k=0;k<row-i;k++)
	{
		B(i+k-1,i+k)=1.0;
	}
	return B*A;
}

/* remove the col i th*/
extern math::matrix<double>RemoveCol(math::matrix<double>A,int i)
{
	int row=A.RowNo();
	math::matrix<double>B(row-1,row);
	int k,j;
	for (k=0;k<i-1;k++)
	{
		B(k,k)=1.0;
	}
	for (k=0;k<row-i;k++)
	{
		B(i+k-1,i+k)=1.0;
	}
	return A*(~B);
}

extern math::matrix<double>EyeMat(int dim)
{
	math::matrix<double> temp(dim,dim);
	for (int i=0;i<dim;i++) temp(i,i)=1.0;
	return temp;
}
/* i th to num*/
extern void SetMatClo(math::matrix<double>&A,int i,double num)
{
	int row=A.RowNo();
	for (int j=0;j<row;j++)  A(j,i-1)=num;
}

extern void ZeroMat(math::matrix<double>&A)
{
	int row=A.RowNo(),col=A.ColNo();
	int i,j;
	for (i=0;i<row;i++)
	{
		for (j=0;j<col;j++)	A(i,j)=0.0;
	}
}
extern math::matrix<double> ZeroMat(int row,int col)
{
	math::matrix<double>t(row,col);
	return  t;
}

extern math::matrix<double>SQDotMat(math::matrix<double>A)
{
	int row=A.RowNo(),col=A.ColNo(),i,j;
	math::matrix<double> B(row,col);
	for (i=0;i<row;i++)
	{
		for (j=0;j<col;j++)
		{
			B(i,j)=A(i,j)*A(i,j);
		}
	}
	return B;
}

/* get diag element of A ,return vector*/
extern math::matrix<double>diagElemMat(math::matrix<double>A)
{
	int row=A.RowNo(),col=A.ColNo(),i;
	int num=min(row,col);
	math::matrix<double> B(num,1);
	for (i=0;i<num;i++)
	{
		B(i,0)=A(i,i);
	}
	return B;
}


extern math::matrix<double> sqrtMat(math::matrix<double>A)
{
	int row=A.RowNo(),col=A.ColNo(),i,j;
	math::matrix<double> B(row,col);
	for (i=0;i<row;i++)
	{
		for (j=0;j<col;j++)
		{
			B(i,j)=sqrt(A(i,j));
		}
	}
	return B;
}
/*
 *get the position of every element in list1,if sorted by ascending  
 */
extern void GetAscendPos(int* list1,int num,int* pos1)
{
	int count=0;
	for (int i=0;i<num;i++)
	{
		count=0;
		for (int j=0;j<num;j++)
		{
			if (list1[i]>list1[j])
			{
				count++;
			}
		}
		pos1[i]			=(count==0)?0:count-1;
	}
}

/*
 *get the position of a in list1 (num x 1)
 *return -1 , if the a doesn't exist in list1
 */
extern int GetPos(int* list1,int a,int num)
{
	int pos=-1;
	for (int i=0;i<num;i++)
	{
		if (list1[i]==a)
		{
			pos=i;
			break;
		}
	}
	return pos;
}

extern int MaxIndex(int* ptr,int tar,int num)
{
	int maxdat=0,maxindex=-1;
	for (int i=0;i<num;i++)
	{
		if(maxdat<ptr[i])
		{
			maxdat=ptr[i];
			maxindex=i;
		}
	}
	return maxindex;
}
/*
 *sort the Qaa
 *I:
 *	Qaa					the vc-mat of a_hat
 *	ambinfo			the ambiguity information, before the sort
 *	
 *O:	
 *	parambinfo		the partial ambiguty info, after the sort
 *
 *return: Q			the sort Q
 */
math::matrix<double> SortQaa(math::matrix<double> Qaa,math::matrix<double>& Qba,math::matrix<double>&ahat,DdAmbInfo ambinfo,DdAmbInfo& parambinfo,DdObsInfo obsinfo)
{
	math::matrix<double> Q=Qaa;
	parambinfo=ambinfo;
	int i,j,k,posmax;
	double maxele=0.0;
	int cnt=0;
	for (i=0;i<ambinfo.freqNum;i++)
	{
		int numAmb=ambinfo.NoUnfix(i);
		int* prnlist=new int[numAmb];
		InitPtr(prnlist,numAmb);
		double* elelist=new double[numAmb];
		InitPtr(elelist,numAmb);
		ambinfo.GetUnfixPrnlist(prnlist,i);
		obsinfo.GetUnfixElelist(elelist,i,ambinfo);
		cnt+= i>0?ambinfo.NoUnfix(i-1) : 0;
		for (j=0;j<numAmb;j++)
		{
			posmax=j;		maxele=elelist[j]; 
			for (k=j+1;k<numAmb;k++)
			{
				if (maxele<elelist[k])
				{
					posmax=k;
					maxele=elelist[k];
				}
			}

			if(posmax!=j) 
			{		 	
				double eleTemp=elelist[posmax];
				elelist[posmax]=elelist[j];
				elelist[j]=eleTemp;
				ChangeRowCol(Q, j+cnt+1, cnt+posmax+1);
				ChangeCol(Qba,j+cnt+1, cnt+posmax+1);
				ChangeRow(ahat, j+cnt+1, cnt+posmax+1);
				parambinfo.InterChange(i,j,posmax);
			}
		}
		delete[] elelist,prnlist;
	}
		return Q;
}

/* 
 *select VC-mat according to the elevations
 * 1.	sort Qaa on the elevations (every freq )
 * 2.  select Qaa on the maskele
 * 
 *maskele   /degree 
 */

extern math::matrix<double> SelectQaa(math::matrix<double>Qaasort,math::matrix<double>Qbasort, math::matrix<double>& QbaPar,math::matrix<double>& ahatPar,
													DdObsInfo obsinfo,DdAmbInfo parambinfo,double maskele)
{

	math::matrix<double> Qf,ahat=ahatPar;
	int* cnt=new int[parambinfo.freqNum];
	//math::matrix<double>* Qptr=new math::matrix<double>[parambinfo.freqNum];
	int  cntpre=0;
	InitPtr(cnt,parambinfo.freqNum);
	int i,j;
	for (i=0;i<parambinfo.freqNum;i++)
	{
		int loop=parambinfo.NoUnfix(i);
		for (j=0;j<loop;j++)
		{
			double eles=obsinfo.GetEle(parambinfo.prnList[i][j]);
			if(maskele*D2R<eles) 
			{
				if(j==loop-1) cnt[i]=loop;
				continue;
			}
			else
			{
				cnt[i]=j;
				break;
			}
		}
		if (cnt[i]==0)
		{
			Qf(0,0)=0.0;
		}
		else
		{
			int t=           (i==0)?0:parambinfo.NoUnfix(i-1);
			cntpre +=		t;
			Qf				=		(i==0) ? GetBlockMat(Qaasort,1,cnt[0],1,cnt[0],1) : DiagMatSym(Qf,GetBlockMat(Qaasort,cntpre+1,cntpre+cnt[i],cntpre+1,cntpre+cnt[i],1)); 
			QbaPar		=		(i==0) ? GetBlockMat(Qbasort,1,3,1,cnt[0],2):ConvergeMat(QbaPar.RowNo(),QbaPar,GetBlockMat(Qbasort,1,3,cntpre+1,cntpre+cnt[i],2)) ;
			ahatPar		=		(i==0) ? GetBlockMat(ahat,1,cnt[0],1,1,2):VecMat(1,ahatPar,GetBlockMat(ahat,cntpre+1,cntpre+cnt[i],1,1,2));
		}
			}
	delete[] cnt;
	return Qf;
 }
	

/*sort list
 *sort list by ascending
 */
extern void SortObsEpochData(ObsEpochData& data)
{
	ObsDataRecord temp;
	int i,j,pos,minelem;
	for (i=0;i<data.sateNum;i++)
	{
		pos=i;
		minelem=data.obsdatarecord[i].PRN;
		for (j=i;j<data.sateNum;j++)
		{
			if (data.obsdatarecord[j].PRN<minelem)
			{
				minelem=data.obsdatarecord[j].PRN;
				pos=j;
			}
		}
		if (i!=pos)
		{
			temp=data.obsdatarecord[pos];
			data.obsdatarecord[pos]=data.obsdatarecord[i];
			data.obsdatarecord[i]=temp;
		}
		
	}
}

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n)
{
	double c=0.0;

	while (--n>=0) c+=a[n]*b[n];
	return c;
}


extern double str2num(const char *s, int i, int n)
{
	double values;
	char str[256],*p=str;

	if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<n) return 0.0;
	for (s+=i;*s&&--n>=0;s++) *p++=*s=='d'||*s=='D'?'E':*s; *p='\0';
	return sscanf(str,"%lf",&values)==1?values:0.0;
}

extern void OutPtr(double* ptr,int num)
{
	for (int i=0;i<num;i++)  cout<<ptr[i]<<"   ";
	cout<<endl;
}

extern void OutPtr(int* ptr,int num)
{
	for (int i=0;i<num;i++)  cout<<ptr[i]<<"   ";
	cout<<endl;
}

extern void Cstr2charPtr(CString kkk,char* code)
{
	
	int nLen=kkk.GetLength();
	int nBytes=WideCharToMultiByte(CP_ACP,0,kkk,nLen,NULL,0,NULL,NULL);
	//code=new char[nBytes+1];
	memset(code,0,nLen+1);
	WideCharToMultiByte(CP_OEMCP,0,kkk,nLen,code,nBytes,NULL,NULL);
	code[nBytes]=0;
}
/*set ptr to the diagonal of matrix
 *I:
 *	ptr		the pointer
 *	num		the number of elements to set
 *	firstPos	the first position in mat bigger than 0
 *O:
 *	mat
 *Note:
 *	num< min(rownum,colnum)
 */
extern void SetPtrToMatdiag(double* ptr,int num,math::matrix<double>&mat,int firstPos)
{
	for (int i=0;i<num;i++)
	{
		mat(firstPos-1+i,firstPos-1+i)=ptr[i];
	}
}
extern int FindPosInt(int* tarlist,int num,int tar)
{
	for (int i=0;i<num;i++)
	{
		if (tarlist[i]!=tar)
		{
			continue;
		}
		else
		{
			return i;
		}
	}
	return -1;
}

/* tar[i]=a[i]*/
extern void PtrEqual(double* a,double* tar ,int num)
{
	int i;
	for (i=0;i<num;i++)
	{
		tar[i]=a[i];
	}
}

/* tar[i]=a[i]*/
extern void PtrEqual(int* a,int* tar ,int num)
{
	int i;
	for (i=0;i<num;i++)
	{
		tar[i]=a[i];
	}
}

extern double DistofVector(double* a,double* b,int num)
{
	double res=0.0;
	double* s=new double[num];
	for (int i=0;i<num;i++)
		s[i]=a[i]-b[i];
	res= Norm(s,num);
	delete[] s;
	return res;
}

/* MW combination unit: cycle
 *I:
 *	phs		unit	cycle
 *	cod		unit	meter
 */
double  MWCom(double* phs, double* cod,int index1,int index2,int sysid)
{
	double freq1=FreqSys(sysid,index1),freq2=FreqSys(sysid,index2);
	double lam1=CLIGHT/freq1,lam2=CLIGHT/freq2;

	double mwCur=(freq1*cod[0]+freq2*cod[1])/(freq1+freq2);
	mwCur-=(freq1*phs[0]*lam1-freq2*phs[1]*lam2)/(freq1-freq2);
	int coefmw[3];
	InitPtr(coefmw,3);
	coefmw[index1]=1;coefmw[index2]=-1;
	double lamMW=CLIGHT/CombFreq(sysid,coefmw);
	return mwCur/lamMW;
}

extern double CombObsCycle(int sysid,int* coef,double* threeobs)
{
	double freq[3];
	double a=0.0;
	for(int i=0;i<3;i++)
	{
		a+=coef[i]*threeobs[i]*(CLIGHT/1E9);
	}
	return a/(CombFreq(sysid,coef)/1E9);
}

/*largest common divisor*/
extern double LCD(int a,int b)
{
	int m=a,n=b,c;
	while (n!=0)
	{
		c=m%n;
		m=n;
		n=c;
	}
	return (double)m;
}

/*least common multiply*/
extern double LCM(int a,int b)
{
	return (double)a*(double)b/LCD(a,b);
}

int AbsIndMinInd(double* a,int num)
{
	int ind=num-1;
	double cmin=fabs(a[num-1]);
	for(int i=0;i<num-1;i++)
	{
		if (cmin>fabs(a[i]))
		{
			cmin= fabs(a[i]);
			ind=	 i;
		}
	}
	return ind;
}
extern void EnterBaseStnCoor( SppInfo& baseInfo,double* coor )
{
	for (int i=0;i<3;i++)
	{
		baseInfo.recPos[i]=coor[i];
	}
}
extern void EnterBaseStnCoor( SppInfoGlo& baseInfo,double* coor )
{
	for (int i=0;i<3;i++)
	{
		baseInfo.recPos[i]=coor[i];
	}
}

static void GetThreeInt(CString line, int* a)
{
	line=line.Trim();
	int len=line.GetLength();
	a[0]=_wtoi(line.Mid(0,3));
	a[2]=_wtoi(line.Right(3));
	a[1]=_wtoi(line.Mid(3,6));

}

static void SetSys1(CStdioFile& Commfile,CString& line,DdCtrl& ddctrl,int* a,int i,int j)
{
	if (line.Find(_T("#SYSID"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.sysid=_wtoi(line);
	}
	if (line.Find(_T("#MASKELE"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.maskele=_wtof(line);
	}
	if (line.Find(_T("#MODE"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.mode=_wtoi(line);
	}
	if (line.Find(_T("#IONOFLAG"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.ionoFlag=_wtoi(line);
	}
	if (line.Find(_T("#TROPFLAG"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.tropFlag=_wtoi(line);
	}
	if (line.Find(_T("#AMBFLAG"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.ambFlag=_wtoi(line);
	}
	if (line.Find(_T("#CODTYPE"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		GetThreeInt(line,a);
		for (i=0;i<3;i++)
		{
			ddctrl.codtype[i]=a[i];
		}
	}
	if (line.Find(_T("#TROPSTEP"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.tropStep=_wtoi(line);
	}
	if (line.Find(_T("#PARMASKELE"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.ddambctrl.parEleMask=_wtof(line);
	}
	if (line.Find(_T("#THRESRATIO"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.ddambctrl.ratiothrsd=_wtof(line);
	}
	if (line.Find(_T("#CODFLAG"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.pseudoFlag=_wtoi(line);
	}
	if (line.Find(_T("#PHSFLAG"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.ddambctrl.flag=_wtoi(line);
	}
	if (line.Find(_T("#PHSTYPE"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		GetThreeInt(line,a);
		for (i=0;i<3;i++)
		{
			ddctrl.phstype[i]=a[i];
		}
	}
	if (line.Find(_T("#CODCOEF"))!=-1)
	{
		for (j=0;j<3;j++)
		{
			Commfile.ReadString(line);
			line=line.Trim();
			GetThreeInt(line,a);
			for (i=0;i<3;i++)
			{
				ddctrl.pseudoCoef[j][i]=a[i];
			}
		}
	}
	if (line.Find(_T("#PHSCOEF"))!=-1)
	{
		for (j=0;j<3;j++)
		{
			Commfile.ReadString(line);
			line=line.Trim();
			GetThreeInt(line,a);
			for (i=0;i<3;i++)
			{
				ddctrl.ddambctrl.coef[j][i]=a[i];
			}
		}
	}
	if (line.Find(_T("#SIGMACOD"))!=-1)
	{
			Commfile.ReadString(line);
			line=line.Trim();
			ddctrl.sigmaCod=_wtof(line);
	}
	if (line.Find(_T("#SIGMAPHS"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.sigmaPhs=_wtof(line);
	}
	if (line.Find(_T("#WEIGHTMODE"))!=-1)
	{
		Commfile.ReadString(line);
		line=line.Trim();
		ddctrl.weightMode=_wtoi(line);
	}
	if (line.Trim()=="")
	{
		Commfile.ReadString(line);
	}
}

static void SetFile(CStdioFile& Commfile,CString& line,InPutFileSet& inputfile)
{
	int i,num;
	if (line.Trim()=="") Commfile.ReadString(line);
	if(line.Find(_T("#EPHEMERIS"))!=-1 )  
	{
		num=_wtoi(line.Right(1));
		i=0;
		while (i<num)
		{
			Commfile.ReadString(line);
			if (line.Trim()!="")
			{
				inputfile.fileEph[i]=line.Trim();
				i++;
			}
		}
		inputfile.numEph=i;
		Commfile.ReadString(line);
	}
	
	if(line.Find(_T("#SP3"))!=-1 )  
	{
		num=_wtoi(line.Right(1));
		i=0;
		while (i<num)
		{
			Commfile.ReadString(line);
			if (line.Trim()!="")
			{
				inputfile.fileSp3[i]=line.Trim();
				i++;
			}
		}
		inputfile.numSp3=i;
		Commfile.ReadString(line);
	}

	if(line.Find(_T("#CLK"))!=-1 )  
	{
		num=_wtoi(line.Right(1));
		i=0;
		while (i<num)
		{
			Commfile.ReadString(line);
			if (line.Trim()!="")
			{
				inputfile.fileClk[i++]=line.Trim();
			}
		}
		inputfile.numClk=i;
		Commfile.ReadString(line);
	}

	if(line.Find(_T("#BASE"))!=-1 )  
	{
		num=_wtoi(line.Right(1));
		i=0;
		while (i<num)
		{
			Commfile.ReadString(line);
			if (line.Trim()!="")
			{
				inputfile.fileBase[i++]=line.Trim();
			}
		}
		inputfile.numBase=i;
		Commfile.ReadString(line);
	}

	if(line.Find(_T("#ROVER"))!=-1 )  
	{
		num=_wtoi(line.Right(1));
		i=0;
		while (i<num)
		{
			Commfile.ReadString(line);
			if (line.Trim()!="")
			{
				inputfile.fileRover[i++]=line.Trim();
			}
		}
		inputfile.numRover=i;
		Commfile.ReadString(line);
	}
}

static void SetCrd(CStdioFile& Commfile,CString& line,double* baseCrd,double* roverCrd)
{
	InitPtr(baseCrd,3);
	InitPtr(roverCrd,3);
	int i;
	Commfile.ReadString(line);
	if(line.Trim()=="")  Commfile.ReadString(line);
	if(line.Find(_T("#BASECORD"))!=-1 )
	{
		if ( _wtoi(line.Right(1))==1)
		{
			Commfile.ReadString(line);
			line=line.Trim();
			int len=line.GetLength();
			for(i=0;i<45-len;i++) line=_T(" ")+line;
			for(i=0;i<3;i++) 
				baseCrd[i]=_wtof(line.Mid(i*15,15)); 
		}
		Commfile.ReadString(line);
	}
	if(line.Find(_T("#ROVERCORD"))!=-1)
	{
	    if ( _wtoi(line.Right(1))==1)
		{
			Commfile.ReadString(line);
			line=line.Trim();
			int len=line.GetLength();
			for(i=0;i<45-len;i++) line=_T(" ")+line;
			for(i=0;i<3;i++) roverCrd[i]=_wtof(line.Mid(i*15,15)); 
		}
	}
}

extern void ReadCommFile(DdCtrl& ddctrl,CString CommandFile)
{
	CStdioFile	Commfile;
	Commfile.Open(CommandFile,CFile::modeRead);
	CString line;
	int flag=0,a[3],i,j;
	while(Commfile.ReadString(line))
	{
		if (line.Find(_T("#BEGINSET"))==-1 && flag==0) continue;
		if (line.Trim()=="" || line.Find(_T("#ENDSET"))!=-1)	continue;
		if (line.Find(_T("#BEGINSET"))!=-1 ) flag=_wtoi(line.Right(1));
		if (flag!=0)	SetSys1(Commfile,line,ddctrl,a,i,j);

	}
	Commfile.Close();
}

extern void ReadCommFile(DdCtrl* ddctrl,InPutFileSet& inputfile,CString CommandFile,double* baseCrd,double* roverCrd)
{
	CStdioFile	Commfile;
	Commfile.Open(CommandFile,CFile::modeRead);
	CString line;
	int flag=0,a[3],i,j,fileflag=0,crdflag=0;
	while(Commfile.ReadString(line))
	{
		if (line.Trim()=="")	 continue;
		if(line.Find(_T("#ENDSET"))!=-1) 
		{
			flag=0;
			continue;
		}
		/* the setting of system  */
		if (line.Find(_T("#BEGINSET"))!=-1 ) flag=_wtoi(line.Right(1));
		if (flag!=0)	SetSys1(Commfile,line,ddctrl[flag-1],a,i,j);

		/*the crd of stations*/
		if(line.Find(_T("#ENDCORD"))!=-1)
		{
			crdflag=0;
			continue;
		}
		if(line.Find(_T("#CORD"))!=-1) crdflag=1;
		if(crdflag!=0) SetCrd(Commfile,line,baseCrd,roverCrd);

		/*the inputfile */
		if ( line.Find(_T("#ENDINPUT"))!=-1)	
		{
			fileflag=0;
			continue;
		}
	//	if (line.Find(_T("#INPUTFILE"))==-1 && fileflag==0) continue;
		if (line.Find(_T("#INPUTFILE"))!=-1 ) fileflag=_wtoi(line.Right(1));
		if (fileflag!=0)	SetFile(Commfile,line,inputfile);
	}
	Commfile.Close();
}

extern void ModifyPath(CString& filepath)
{
	int pos1=0, pos2=0;
	pos1=filepath.Find(_T("\\"),pos2+1);
	pos2=filepath.Find(_T("\\"),pos1+1);
	if (pos2-pos1!=1)
	{
		filepath.Replace(_T("\\"),_T("/"));
	}
}

extern void BrowseFile(CString& strDir)
{
	ModifyPath(strDir);
	CFileFind fs;
	bool		   bFound;
	if(strDir.Find(_T("/")) !=-1) bFound=fs.FindFile(strDir+_T("/*.*"));
	if(strDir.Find(_T("\\")) !=-1) bFound=fs.FindFile(strDir+_T("\\*.*"));	
		
	while(bFound)
	{
		bFound=fs.FindNextFile();
		if (!fs.IsDirectory()&&!fs.IsDots() )
		{
			CString fpath;int sk;
			fpath=fs.GetFilePath();
			sk=0;
		}
	}
}

extern void GetOutPath(CString commfile,CString& filepath)
{
	CStdioFile file;  CString line;
	file.Open(commfile,CFile::modeRead);
	while(file.ReadString(line))
	{
		line=line.Trim();
		if(line.Find(_T("#OUTPUTFILE"))!=-1)
		{
			file.ReadString(line);
			line=line.Trim();
			if(line.Find(_T("#ENDOUTPUT"))!=-1)	line=_T("SelfOut.txt");
			break;
		}
	}
	file.Close();
	ModifyPath(line);
	filepath=line;
}

extern void SppFileOut(fstream& fout,SppInfo sppinfo)
{
	for(int i=0;i<3;i++)
	{	
		fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(14)<<sppinfo.recPos[i]<<"   ";
	}
	fout<<sppinfo.dtr<<"   "<<sppinfo.gdop<<"   "<<sppinfo.validnum;
	fout<<endl;
}

extern void SppFileOut(fstream& fout,SppInfoGlo sppinfo)
{
	for(int i=0;i<3;i++)
	{	
		fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(14)<<sppinfo.recPos[i]<<"   ";
	}
	fout<<sppinfo.dtr<<"  "<<sppinfo.gdop<<sppinfo.validnum;
	fout<<endl;
}

extern void FileOutTropsIono(fstream& fout,math::matrix<double> IonoVector,math::matrix<double> TropVector,DdData curData,double* ptrMapCur)
{
	int num=curData.pairNum;
	fout<<setw(2)<<curData.pairNum<<"  "<<curData.refPrn<<"  ";
	fout<<setw(4)<<curData.week<<"  ";
	fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(14)<<curData.sec<<endl;

	for (int i=0;i<num;i++)
	{
		fout<<"G";
		fout<<setw(2)<<curData.rovPrn[i]<<"  ";
		fout<<setiosflags(ios::fixed)<<setprecision(5)<<setw(9)<<ptrMapCur[i]<<"  "<<TropVector(0,0)<<"  "<<IonoVector(i,0)<<endl;
	}
}

extern void FileOutTropsIonoFloat(fstream& fout,math::matrix<double> IonoVector,math::matrix<double> TropVector,DdData curData,double* ptrMapCur,int flag)
{
	int num=curData.pairNum;
	fout<<setw(2)<<curData.pairNum<<"  "<<curData.refPrn<<"  ";
	fout<<setw(4)<<curData.week<<"  ";
	fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(14)<<curData.sec<<endl;

	for (int i=0;i<num;i++)
	{
		fout<<"G";
		fout<<setw(2)<<curData.rovPrn[i]<<"  ";
		fout<<setiosflags(ios::fixed)<<setprecision(5)<<setw(9)<<ptrMapCur[i]*TropVector(0,0)+curData.tropCor[i]<<"  "<<IonoVector(i,0)<<endl;
	}
}
/*
 *b_updated contains the iono and trop after being updated

 */
extern void FileOutTropsIonoFix(fstream& fout,math::matrix<double> b_updated,DdData curData,double* ptrMapCur,int flag)
{
	int num=curData.pairNum;
	fout<<setw(2)<<curData.pairNum<<"  "<<curData.refPrn<<"  ";
	fout<<setw(4)<<curData.week<<"  ";
	fout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(14)<<curData.sec<<"  "<<1<<endl;

	for (int i=0;i<num;i++)
	{
		fout<<"G";
		fout<<setw(2)<<curData.rovPrn[i]<<"  ";
		fout<<setiosflags(ios::fixed)<<setprecision(5)<<setw(9)<<ptrMapCur[i]*b_updated(num,0)+curData.tropCor[i]<<"  "<<b_updated(i,0)<<endl;
	}
}
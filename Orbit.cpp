#include "stdafx.h"//changed
#include "Orbit.h"







/*
	

*/
void CalcSatePos::CompSatePosVelAcc(math::matrix<double>& SatePos,math::matrix<double>& SateVel,math::matrix<double>& SateAcc,double& ek,BroadEphData eph,double t)
{
	// t is sampling time
	int		sysid		=Prn2Sysid(eph.prn);
	double Aa			=SQ(eph.sqA);
	double gm		=sysid>2?GM_BDS:GM;
	double n			=sqrt(GM/pow(Aa,3)) + eph.deltan;
	double tk			=t-eph.toe;
	if(tk<-302400.0)
	{
		tk+=604800.0;
	}
	else if (tk>302400)
	{
		tk-=604800.0;
	}
	double Mk		=eph.m0 + n*tk;
	double Ek			=Mk;
	double Ek0		=Ek+1.0;
	while (fabs(Ek-Ek0)>1.0e-12)
	{
		Ek0		=Ek;
		Ek			=Mk+eph.e*sin(Ek0);
	}
	ek=Ek;

	double vk			=2*atan(sqrt((1.0+eph.e)/(1-eph.e))*tan(Ek/2));
	double uk			=vk	+eph.omega;
	double detauk	=eph.cuc*cos(2*uk)	+	eph.cus*sin(2*uk);
	double detark	=eph.crc*cos(2*uk)	+	eph.crs*sin(2*uk);
	double detaik	=eph.cic*cos(2*uk)	+	eph.cis*sin(2*uk);

	double u=		uk+detauk;
	double r=			Aa*(1-eph.e*cos(Ek))	+	detark;
	double ik=		eph.i0	+	eph.idot*tk+detaik;

	double x=			r*cos(u);
	double y=			r*sin(u);

	double omegaE	=	eph.prn<100?OMEGAE_BDS:OMEGAE;
		

	if(eph.prn>200 && eph.prn<=205)
	{
		double OMG		=	eph.Omega	+eph.OmegaDot*tk-omegaE*eph.toe;
		math::matrix<double> temp(3,1);
		temp(0,0)			=	cos(OMG)*x	-	y*cos(ik)*sin(OMG);
		temp(1,0)			=	sin(OMG)*x	+	y*cos(ik)*cos(OMG);
		temp(2,0)			=	y*sin(ik);
		double a			=	(-5.0/180.0)*PI;
		math::matrix<double> Ra(3,3);
		Ra(0,0)	=1;
		Ra(1,1)	=cos(a);	Ra(1,2)	=sin(a);	
		Ra(2,1)	=-sin(a);	Ra(2,2)	=cos(a);
		double b			=	OMEGAE_BDS*tk;
		math::matrix<double> Rb(3,3);
		Rb(0,0)	=cos(b);	Rb(0,1)	=sin(b);
		Rb(1,0)	=-sin(b);	Rb(1,1)	=cos(b);
		Rb(2,2)	=1;
		temp=Rb*Ra*temp;
		SatePos			=temp;

	}
	else if (eph.prn>205)
	{
		double OMG			=	eph.Omega + (eph.OmegaDot - omegaE)*tk - omegaE*eph.toe;
		
		SatePos(0,0)			=	cos(OMG)*x	-	y*cos(ik)*sin(OMG);
		SatePos(1,0)		=	sin(OMG)*x	+	y*cos(ik)*cos(OMG);
		SatePos(2,0)		=	y*sin(ik);
	}
	else
	{
		double OMG			=	eph.Omega + (eph.OmegaDot - omegaE)*tk - omegaE*eph.toe;
		/*velocity in ECEF*/
		double Ekdot = n/(1.0 - eph.e*cos(Ek)); 
		double vkdot = sin(Ek)*Ekdot*(1.0 + eph.e*cos(vk))/(sin(vk)*(1.0-eph.e*cos(Ek)));
		double dukdot = vkdot + 2.0*(eph.cus*cos(2.0*u)-eph.cuc*sin(2.0*u))*vkdot;
		double drkdot = Aa*eph.e*sin(Ek)*n/(1.0-eph.e*cos(Ek)) + 2.0*(eph.crs*cos(2.0*u)-eph.crc*sin(2.0*u))*vkdot;
		double dikdot = eph.idot + (eph.cis*cos(2.0*u)-eph.cic*sin(2.0*u))*2.0*vkdot;

		double xdot = drkdot*cos(u) - y*dukdot; 
		double ydot = drkdot*sin(u) + x*dukdot; 
		double OMGdot = eph.OmegaDot - omegaE;
		/*acceleration in ECEF*/
		double rk  = Aa*(1-eph.e*cos(Ek))+detark;
		double vkd  = sqrt(1-SQ(eph.e))/(1 - eph.e*cos(Ek))*Ekdot;
		double ikd = eph.idot + 2*vkd*(eph.cis*cos(2*uk)-eph.cic *sin(2*uk));  
		double ukd  = vkd + 2*vkd*(eph.cus*cos(2*uk) - eph.cuc*sin(2*uk));
		double rkd  = Aa*eph.e*Ekdot*sin(Ek)+ 2*vkd*(eph.crs*cos(2*uk)-eph.crc*sin(2*uk));
		double vkdd = -2*eph.e*sqrt(1-SQ(eph.e))*sin(Ek)/SQ((1 -eph.e*cos(Ek)))*SQ(Ekdot);             //Second derivative of vk, true anomaly [rad/s^2]
		double ukdd = vkdd - 4*SQ(vkd)*(eph.cus*sin(2*uk) + eph.cuc*cos(2*uk))+2*vkdd*(eph.cus*cos(2*uk)-eph.cuc*sin(2*uk));// Second derivative of uk, argument of latitude [rad/s^2]
		double rkdd = Aa*eph.e*(cos(Ek)- eph.e)/(1-eph.e*cos(Ek))*SQ(Ekdot)-4*SQ(vkd)*(eph.crs*sin(2*uk)+eph.crc*cos(2*uk))+ 2*vkdd*(eph.crs*cos(2*uk)-eph.crc*sin(2*uk));      //Second derivative of rk, radius [m/s^2]

		
		SatePos(0,0)			=	cos(OMG)*x	-	y*cos(ik)*sin(OMG);
		SatePos(1,0)		=	sin(OMG)*x	+	y*cos(ik)*cos(OMG);
		SatePos(2,0)		=	y*sin(ik);

		SateVel(0,0) = (xdot-y*cos(ik)*OMGdot)*cos(OMG) - (x*OMGdot+ydot*cos(ik)-y*sin(ik)*dikdot)*sin(OMG);
		SateVel(1,0) = (xdot-y*cos(ik)*OMGdot)*sin(OMG) + (x*OMGdot+ydot*cos(ik)-y*sin(ik)*dikdot)*cos(OMG);
		SateVel(2,0) = ydot*sin(ik) + y*cos(ik)*dikdot;

		double Rp[3]		={cos(OMG),sin(OMG),0};
		double Rq[3]		={-cos(ik)*Rp[1],cos(ik)*Rp[0],sin(ik)};
		double Rpd[3]	={-eph.OmegaDot*Rp[1], eph.OmegaDot*Rp[0], 0}; 
		double Rqd[3]	= {-eph.OmegaDot*Rq[1]+ikd*Rq[2]*Rp[1],  eph.OmegaDot*Rq[0] - ikd*Rq[2]*Rp[0], ikd*cos(ik)};   

		double xorb[2]  = {rk*cos(u), rk*sin(u)};   // Position in orbital plane
		double xorbd[2] = {rkd*cos(u)-rk*ukd*sin(u), rkd*sin(u)+rk*ukd*cos(u)};  
		double xorbdd[2]= {(rkdd-rk*SQ(ukd))*cos(u)-(2*rkd*ukd+rk*ukdd)*sin(u),
										(rkdd-rk*SQ(ukd))*sin(u)+(2*rkd*ukd+rk*ukdd)*cos(u)};// Second dervative of xorb, pos in orbital plane [m/sec^2]
		double ikdd  = -4*SQ(vkd)*(eph.cis*sin(2*uk) + eph.cic *cos(2*uk))+2*vkdd*(eph.cis*cos(2*uk) - eph.cic *sin(2*uk)); // Second derivative of ik, inclination angle [rad/s^2]

		double Rpdd[3]   = {-SQ(eph.OmegaDot)*Rp[0], -SQ(eph.OmegaDot)*Rp[1], 0};   // Second derivative of Rp
		double Rqdd[3]   = {-(SQ(eph.OmegaDot) + SQ(ikd))*Rq[0]+ 2*ikd*eph.OmegaDot*Rq[2]*Rp[0] + ikdd*Rq[2]*Rp[1],
										-(SQ(eph.OmegaDot) + SQ(ikd))*Rq[1]+ 2*ikd*eph.OmegaDot*Rq[2]*Rp[1] - ikdd*Rq[2]*Rp[0],
										-SQ(ikd)*Rq[2] + ikdd*cos(ik)};
		for(int i=0;i<3;i++)
		{
				SateAcc(i,0)	=		Rpdd[i]*xorb[0]+Rqdd[i]*xorb[1]+2*(Rpd[i]*xorbd[0]+Rqd[i]*xorbd[1])+Rp[i]*xorbdd[0]+Rq[i]*xorbdd[1];
		}
	}

}

void CalcSatePos::CompSateClk(double& SateClk,double ek, BroadEphData eph,double transt,int freqIndex)
{
	double F = -4.442807633*10e-10;
	double Dtr=	F*eph.e*eph.sqA*sin(ek);//relaticistic corr
	double tc	=	transt-eph.toc;

	int		sysid		=Prn2Sysid(eph.prn);
	double b=FreqSquareRatio(sysid,freqIndex);
	if (tc > 302400)
	{
		tc = tc - 604800;
	}
	else if (tc < -302400)
	{
		tc = tc + 604800;
	}
	SateClk=eph.satClkBs+eph.satClkDrt*tc+eph.satClkDrtRe*pow(tc,2)+Dtr-eph.tgd*b;
}

	/*	interatively compute the satepos and sate clk
	input:
			eph				single sate eph info
			StnPos			receiver/station position
			t					obs time
			freqIndex		the index of frequency  0=IF 1=L1/B1 ;2=L2/B2;	3=L5/B3 

			SatePos		vector 3*1;
			SateClk			
	NOTE:
			this fuction could compute all type sates, but it is used to process the GEO(BDS) in this class
			other satepos is computed without interation
	*/
void	CalcSatePos::InteSatePosClk(math::matrix<double>& SatePos,double& SateClk,BroadEphData eph,math::matrix<double> StnPos,
													double t,int freqIndex)
{
	double dt0			=	0.075;
	double ek				=	0.0;
	math::matrix<double> SateVel(3,1);
	math::matrix<double> SateAcc(3,1);
	CompSatePosVelAcc(SatePos,SateVel,SateAcc,ek,eph,t);

	int		sysid			=Prn2Sysid(eph.prn);
	double omegaE	=sysid>2?OMEGAE_BDS:OMEGAE;
	double theta			=	dt0*omegaE;
	math::matrix<double> R(3,3);
	RotationMatrix3(R,theta,3);

	SatePos=R*SatePos;

	math::matrix<double> Deltacrd(3,1);
		Deltacrd		=		SatePos-StnPos;
	double dt1;
				dt1		=		Deltacrd.Norm()/CLIGHT;

	while(fabs(dt1-dt0)>1.0e-12)
	{
		dt0			=	dt1;
		CompSatePosVelAcc(SatePos,SateVel,SateAcc,ek,eph,t-dt0);
		theta		=	dt0*omegaE;
		RotationMatrix3(R,theta,3);
		SatePos	=R*SatePos;
		Deltacrd	=	SatePos-StnPos;
		dt1			=	Deltacrd.Norm()/CLIGHT;
	}
	CompSateClk(SateClk, ek, eph, t-dt1,freqIndex);
}
/*
	I:
		t	 updated time	=	t0-dtr(i)
		t0  obs time
	O:
		ek
	NOTE:
		there is problem in GEO(BDS) velocity and acc ,so calc pos by interation
*/
void	CalcSatePos::InteSatePosClkVA(math::matrix<double>& SatePos,double& SateClk,BroadEphData eph,math::matrix<double> StnPos,
														double t,double t0,double ek,math::matrix<double> Xcrd,math::matrix<double> Xvel,math::matrix<double> Xacc)
{
	double t_t0	=t-t0;// usually =0;
	if (t_t0>=302400)
	{
		t_t0 = t_t0 - 604800;
	}
	else if (t_t0< -302400)
	{
		t_t0 = t_t0 + 604800;
	}
	double tau0			=	0.075;
	double dt				=	t_t0-tau0;

	SatePos	=Xcrd+dt*Xvel+0.5*SQ(dt)*Xacc;

	int		sysid			=Prn2Sysid(eph.prn);
	double omegaE	=sysid>2?OMEGAE_BDS:OMEGAE;
	double theta			=	tau0*omegaE;
	math::matrix<double> R(3,3); //correction of the earth rotation
	RotationMatrix3(R,theta,3);

	SatePos=R*SatePos;

	math::matrix<double> Deltacrd(3,1);
		Deltacrd		=		SatePos-StnPos;
	double tau1;
				tau1		=		Deltacrd.Norm()/CLIGHT;

	while(fabs(tau1-tau0)>1.0e-12)
	{
		tau0			=	tau1;
		dt				=t_t0-tau0;

		SatePos	=Xcrd+dt*Xvel+0.5*Xacc*SQ(dt);

		theta		=	tau0*omegaE;
		RotationMatrix3(R,theta,3);
		SatePos	=R*SatePos;
		Deltacrd	=	SatePos-StnPos;
		tau1			=	Deltacrd.Norm()/CLIGHT;
	}
	CompSateClk(SateClk, ek, eph, t-tau1,0);//freqIndex
}
/*
	Main function in this class
	I:
		RecPos		receiver position
		t				obs time-clock error of receiver
		eph
		
	O:
		pos			sate pos
		clkerr		clock error of sate	
		
	Note:
		calc with mix
*/
bool CalcSatePos::SatePosClk(double* pos,double& SateClk,double* RecPos,double t,
											int freqIndex,BroadEphData eph)
{
	if (eph.satHealth!=0)
	{
		return false;
	}
	math::matrix<double> SatePos(3,1);
	math::matrix<double> StnPos(3,1);
	math::matrix<double> Xcrd(3,1);
	math::matrix<double> Xvel(3,1);
	math::matrix<double> Xacc(3,1);
	for(int i=0;i<3;i++)  StnPos(i,0)=RecPos[i];

	if(eph.prn>200 && eph.prn<206)
	{
		
		InteSatePosClk(SatePos,SateClk,eph,StnPos,t,freqIndex);
	}
	else
	{
		double ek;
		double t0=t;
		CompSatePosVelAcc(Xcrd, Xvel, Xacc, ek, eph, t);
		InteSatePosClkVA(SatePos,SateClk,eph,StnPos, t, t0, ek, Xcrd, Xvel, Xacc);
	}
	for(int i=0;i<3;i++)  pos[i]=SatePos(i,0);

	return true;
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
	/* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
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
* args   : gtime_t time     I   time by satellite clock (gpst)
*          geph_t *geph     I   glonass ephemeris
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
double CalcSatePos:: geph2clk(double ts, const BroadEphDataGlo geph)
{
	double t;
	int i;

	t=ts-geph.toe;

	for (i=0;i<2;i++) {
		t-=(geph.NegativeTauN+geph.PositiveGammaN*t);
	}
	return geph.NegativeTauN+geph.PositiveGammaN*t;
}
/* glonass ephemeris to satellite position and clock bias ----------------------
* compute satellite position and clock bias with glonass ephemeris
* args   : gtime_t time     I   time (gpst)
*          geph_t *geph     I   glonass ephemeris
*          double *rs       O   satellite position {x,y,z} (ecef) (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
void CalcSatePos::  geph2pos(double ts, const BroadEphDataGlo geph, double *rs, double& dts)
{
	double t,tt,x[6];
	int i;


	t=ts-geph.toe;

	dts=geph.NegativeTauN+geph.PositiveGammaN*t;

	for (i=0;i<3;i++) {
		x[i  ]=geph.Pos[i];
		x[i+3]=geph.Vel[i];
	}
	for (tt=t<0.0?-TSTEP:TSTEP;fabs(t)>1E-9;t-=tt) {
		if (fabs(t)<TSTEP) tt=t;
		glorbit(tt,x,geph.Acc);
	}
	for (i=0;i<3;i++) rs[i]=x[i];

}



/************************************************************************/


/********************************SBAS****************************************/
/* sbas ephemeris to satellite clock bias --------------------------------------
* compute satellite clock bias with sbas ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          seph_t *seph     I   sbas ephemeris
* return : satellite clock bias (s)
*-----------------------------------------------------------------------------*/
 double CalcSatePos:: seph2clk(double time, const BroadEphDataSBAS seph)
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
void CalcSatePos::seph2pos(double time, const BroadEphDataSBAS seph, double *rs, double *dts)
{
	double t;
	int i;

	t=timediff(time,seph.t0);

	for (i=0;i<3;i++) {
		rs[i]=seph.pos[i]+seph.vel[i]*t+seph.acc[i]*t*t/2.0;
	}
	*dts=seph.a0+seph.a1*t;

}
/*************************************************************************/

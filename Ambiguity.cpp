#include "stdafx.h"//changed
#include "Ambiguity.h"
#include <math.h>
#include "matrix.h"
#include "ExtFun.h"


/*
This file defines the class about Ambiguity
		TCAR 2 methods
		CIR
		Lambda
	
	all the class is compatible to double/single freqencies
*/






/*
	get the freq of sys
	I:
		sysid		system id
	O:
		freq		three frequencys of signals
					GPS:		L1 L2 L5
					BDS:		B1 B3 B2
					GAL:		to be written  5 frequencies
		
*/
void Combination::Freq(int sysid,double* freq)
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
		freq[1]=FREQ6_BDS;//B3  B2 < B3
		freq[2]=FREQ7_BDS;//B2
	}
}

void Combination::FreqRatio(int sysid,double* freq)
{
	if (sysid==1)
	{
		freq[0]=FREQ1RATIO;
		freq[1]=FREQ2RATIO;
		freq[2]=FREQ5RATIO;
	}
	else if(sysid==5)
	{
		freq[0]=FREQ2_BDSRATIO;
		freq[1]=FREQ6_BDSRATIO;//B3  B2 < B3
		freq[2]=FREQ7_BDSRATIO;//B2
	}
}
/*
	get the wavelength of three frequency
	I:
		sysid		system id
	O:
		Wave    three wavelength of signals
*/
void Combination::Wave(int sysid,double* Wave)
{
	double* freq=new double[3];
	Freq(sysid,freq);
	for(int i=0;i<3;i++)
	{
		Wave[i]=CLIGHT/freq[i];
	}
	delete[] freq;
}

/*
	Frequency of Combination
*/
double Combination::CombFreq(int sysid,int* coef)
{
	double* freq=new double[3];
	double comb;
	comb=0.0;
	Freq(sysid,freq);
	for (int i=0;i<3;i++)
	{
		comb+=coef[i]*freq[i];
	}
	delete[] freq;
	return comb;
}

/*
	wavelength of Combination
*/
double Combination::CombWave(int sysid,int* coef)
{
	return CLIGHT/CombFreq(sysid,coef);
}

/*
	observable of Combination
	threeobs		(m)
*/
double Combination::CombObs(int sysid,int* coef,double* threeobs)
{
	double* freq=new double[3];
	FreqRatio(sysid,freq);
	double a=0.0;
	for(int i=0;i<3;i++)
	{
		a+=coef[i]*freq[i]*threeobs[i];
	}
	delete[] freq;
	return a/(CombFreq(sysid,coef));
}

/*
	Ionosphere factor of Combination relative to f1
*/
double Combination::CombIono(int sysid,int* coef)
{
	double* freq=new double[3];
	FreqRatio(sysid,freq);
	double a=(coef[0]*freq[0]+pow(freq[0],2)*coef[1]/freq[1] +pow(freq[0],2)*coef[2]/freq[2] );
	delete[] freq;
	return  a*1000000000/CombFreq(sysid,coef);
}
/*	
	combine the IF obs for different system
	I:
		sysid
		obs		3*1 (m)
		index	2*1	the frequencies index of two frequency
*/
double Combination::IonoFree(int sysid,double* obs,int* index)
{
	double* freq=new double[3];
	FreqRatio(sysid,freq);
	double diff=SQ(freq[index[0]])-SQ(freq[index[1]]);
	return (SQ(freq[index[0]])/diff*obs[index[0]]-SQ(freq[index[1]])/diff*obs[index[1]]);
}

/*
	resolve the ambiguity with previous tcar 
		I:
			obsP   the observables of pseudorange on three freqs  3x1  (m)
			obsL	  the boservables of phase on three freqs				3x1  (m)
			sysid  system id
		O:
			NarrAmb	ambiguities of three freqs								3x1	(cycle)
			
		Reference: 
			Forssell,B.	Carrier Phase Ambiguity Resolution in GNSS-2 ION GPS. Vol. 2, 1997
		Note:
			GF model, using high precise pseudorange in program P5 is recommended
*/
void Ambiguity::TCARPre(double* obsP,double* obsL,int sysid, double* NarrAmb)
{
	
	int* coef=new int[3];
	double* lamb=new double[3];
	Combination tcarComb;
	tcarComb.Freq(sysid,lamb);

	//     here the combination is based on f3 so the iono is changed to L3
	//		step 1: (0,1,-1)  EWL
	coef[0]=0;	coef[1]=1;coef[2]=-1;

	double lamb23		=tcarComb.CombWave(sysid,coef);
	double obsL23		=tcarComb.CombObs(sysid,coef,obsL);
	double Iono23		=SQ(lamb[0])/SQ(lamb[2])*tcarComb.CombIono(sysid,coef);

	double		N23			=ROUND( (obsP[2]-obsL23)/lamb23 );

	//		step 2: (1,0,-1)   WL
	coef[0]=1;	coef[1]=0;coef[2]=-1;

	double lamb13		=tcarComb.CombWave(sysid,coef);
	double obsL13		=tcarComb.CombObs(sysid,coef,obsL);
	double Iono13		=SQ(lamb[0])/SQ(lamb[2])*tcarComb.CombIono(sysid,coef);

	double		N13			=ROUND( (N23*lamb23+obsL23-obsL13)/lamb13 );

	//		step 3:(1,0,0)  estimate ionofactor relative to L3 with L23 and L13
	double Iono_hat	=(N23*lamb23+obsL23-N13*lamb13+obsL13)/(Iono13-Iono23);

	NarrAmb[2]			=ROUND((N13*lamb13+obsL13-Iono_hat*(1-Iono13)-obsL[2])/lamb[2]);
	NarrAmb[1]			=-(NarrAmb[2]-N23);
	NarrAmb[0]			=-(NarrAmb[2]-N13);
	delete[] coef;
	delete[] lamb;
}


/*
	whether the error probability(Pe) < the required Pe, rounding
	I:
		delta		float solution away from the true integer		(m)
		Wave	the wavelength of phase								(m)
		sigma	standard deviation										(m)
		PeReq	required Pe
	Return:
		Pe<PeReq true; else false
	Note:
		all inputs are positive
*/
bool Ambiguity::IsReq(double delta,double Wave,double sigma,double PeReq)
{
	double a=delta/sigma;
	double b=Wave/sigma;
	double c=(b+2*log(PeReq)/b)/2;
	return a<c?true:false;
}



/*
	compute the error probability of rounding 
	I:
		delta		float soluton away from the true integer		(m)
		Wave	the wavelength of phase								(m)
		sigma	standard deviation										(m)
	Return:
		probability of rounding to an incorrect integer
	Note:
		all inputs are positive
*/
double Ambiguity::PrErrRound(double delta,double Wave,double sigma)
{
	return exp(-1/(2*SQ(sigma))*(SQ(Wave)-2*Wave*delta));
}

/*
	AR with Cascade Integer Resolution  (CIR)
	I:
		obsP			pseudorange observables on three freqs	3x1	(m)
		obsL			phase observables on three freqs				3x1	(m)
		sysid
	O:
		NarrAmb	narrow ambiguity										3x1	(cycles)
	Reference:
		Jung Jaewoo. High Integrity Carrier Phase Navigation for Future LAAS Using Multiple Civilian GPS Signals.ION GPS. 1999
		Jung Optimization of Cascade Integer Resolution with Three Civil GPS Frequencies. ION GPS, 2000
	Note:
		GF model  (baseline<10km)
		P5 is recommended
		Here, ML=EWL+WL
*/
bool Ambiguity::CIR(double* obsP,double* obsL,int sysid, double* NarrAmb)
{
	
	int* coef=new int[3];
	double* lamb=new double[3];
	Combination cirComb;
	cirComb.Freq(sysid,lamb);

	/* step 1: EWL	(0,1,-1) */
	coef[0]=0;coef[1]=1;coef[2]=-1;

	double lamb23		=cirComb.CombWave(sysid,coef);
	double obsL23		=cirComb.CombObs(sysid,coef,obsL);

	double N23					=ROUND((obsP[2]-obsL23)/lamb23);

	/* step 2: WL		(1,-1,0) */
	coef[0]=1;coef[1]=-1;coef[2]=0;

	double lamb12		=cirComb.CombWave(sysid,coef);
	double obsL12		=cirComb.CombObs(sysid,coef,obsL);

	double N12					=ROUND((N23*lamb23+obsL23-obsL12)/lamb12);

	/* step 3: ML		(1,0,-1)    check ML=EWL+WL   */
	coef[0]=1;coef[1]=0;coef[2]=-1;

	double lamb13		=cirComb.CombWave(sysid,coef);
	double obsL13		=cirComb.CombObs(sysid,coef,obsL);

	double N13					=ROUND((N12*lamb12+obsL12-obsL13)/lamb13);
	
	if(N13!=(N12+N23))		return false;

	/* step 4: NL		 			*/
	NarrAmb[0]	=ROUND((N12*lamb12+obsL12-obsL13)/lamb[0]);
	NarrAmb[1]	=NarrAmb[0]-N12;
	NarrAmb[2]	=NarrAmb[0]-N13;

	delete[] coef;
	delete[] lamb;

	return true;
}

/* 
*	LD factorization (Q=L*diag(D)*L')
*	笑渐不闻声渐悄，多情却被无情恼
*	I:
*		n		size of Q
*		Q		cofactor matrix of a_hat
*	O:
*		L	upper matrix    n*n
*		D	diagonal matrix n
*/
bool Ambiguity::LTDLFactorization(int n, math::matrix<double>  Q, math::matrix<double> &L, double* D)
{
	int i,j,k;
	double a;
	math::matrix<double>  A=Q;

	for (i=n-1;i>=0;i--) 
	{
		if ((D[i]=A(i,i))<=0.0) return false;
		a=sqrt(D[i]);
		for (j=0;j<=i;j++) L(j,i)=A(j,i)/a;
		for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A(k,j)-=L(k,i)*L(j,i);
		for (j=0;j<=i;j++) L(j,i)/=L(i,i);
	}

	return true;
}
void Ambiguity::Gauss(int n, math::matrix<double> & L, math::matrix<double> & Z, int i, int j)
{
	int k,mu;
	if (( mu=(int)ROUND(L(j,i)) ) !=0) 
	{
		for (k=i;k<n;k++) L(j,k)-=(double)mu*L(i,k);
		for (k=0;k<n;k++) 
		{
			Z(j,k)-=(double)mu*Z(i,k);
			if(fabs(Z(j,k))>1e5)
				Z(j,k)=1e5;
		}
	}
}

void Ambiguity::Perm(int n, math::matrix<double> &L, double* D, int j, double del, math::matrix<double> &Z)
{
	int k;
	double eta,lam,a0,a1;

	eta=D[j]/del;
	lam=D[j+1]*L(j,j+1)/del;
	D[j]=eta*D[j+1]; D[j+1]=del;
	for (k=0;k<=j-1;k++) 
	{
		a0=L(k,j); a1=L(k,j+1);
		L(k,j)=-L(j,j+1)*a0+a1;
		L(k,j+1)=eta*a0+lam*a1;
	}
	L(j,j+1)=lam;
	for (k=j+2;k<n;k++) SWAP(L(j,k),L(j+1,k));
	for (k=0;k<n;k++) SWAP(Z(j,k),Z(j+1,k));
}

/*
 *  reduction  de-correlation
 *  I:
 *	  D  decomposition of LD
 *	  L   
 *  O:
 *	  Z  transformation matrix
 */
void Ambiguity::Reduction(int n, math::matrix<double>& L, double* D, math::matrix<double> & Z)
{
	int i,j,k;
	double del;
	j=n-2; k=n-2;
	while (j>=0) 
	{
		if (j<=k) 
			for (i=j+1;i<n;i++) 
				Gauss(n,L,Z,i,j);
		del=D[j]+L(j,j+1)*L(j,j+1)*D[j+1];
		if (del+1E-6<D[j+1]) 			/* compared considering numerical error */
		{ 
			Perm(n,L,D,j,del,Z);
			k=j; 
			j=n-2;
		}
		else j--;
	}
}

/*
	search m best integer vectors zn and correspond quadratic residual error 
	I:
		 天长地久有时尽，此恨绵绵无绝期 
		n		dimension of ambiguity vector
		m		number of ambiguity vectors to output
		L
		D		L'*D*L=Qzz	Qzz is  covariance matrix after Z-Trans
		zs		float ambiguity
	O:
		zn		integer ambiguity vector m*n, every row is one result
		s		quadratic residual error of m sets ambiguity vectors

*/
int Ambiguity::search(int n, int m, math::matrix<double> L,  double*  D,	
								math::matrix<double> zs,math::matrix<double>& zn,  double* s)
{
	int i,j,k,c,nn=0,imax=0;
	double newdist,maxdist=1E99,y;
	math::matrix<double>S(n,n); 
	double* dist=new double[n];
	double* zb=new double[n];
	double* z=new double[n];
	double* step=new double[n];
	

	k=n-1; dist[k]=0.0;
	zb[k]=zs(k,0);
	z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
	for (c=0;c<LOOPMAX;c++) 
	{
		newdist=dist[k]+y*y/D[k];
		if (newdist<maxdist) 
		{
			if (k!=0) 
			{
				dist[--k]=newdist;
				for (i=0;i<=k;i++)
					S(i,k)=S(i,k+1)+(z[k+1]-zb[k+1])*L(i,k+1);
				zb[k]=zs(k,0)+S(k,k);
				z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
			}
			else {
				if (nn<m) 
				{
					if (nn==0||newdist>s[imax]) imax=nn;
					for (i=0;i<n;i++) zn(nn,i)=z[i];
					s[nn++]=newdist;
				}
				else 
				{
					if (newdist<s[imax]) 
					{
						for (i=0;i<n;i++) zn(imax,i)=z[i];
						s[imax]=newdist;
						for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
					}
					maxdist=s[imax];
				}
				z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
			}
		}
		else 
		{
			if (k==n-1) break;
			else 
			{
				k++;
				z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
			}
		}
	}
	for (i=0;i<m-1;i++) 
	{   /* sort by s */
		for (j=i+1;j<m;j++) 
		{
			if (s[i]<s[j]) continue;
			SWAP(s[i],s[j]);
			for (k=0;k<n;k++) SWAP(zn(i,k),zn(j,k));
		}
	}

	if (c>=LOOPMAX) 
	{
		delete[] dist,zb,z,step;
		return -1;
	}
	delete[] dist,zb,z,step;
	return 0;
}

/*
	least squares ambiguity decorrelation adjustment
	I:
		n		dimension of ambiguity vector
		m		control the number of vectors to be output
		a		float ambiguity  n*1
		Q		vc-matrix of float ambiguity n*n
	O:
		F		matrix of output integer ambiguity n*m  
		s		quadratic residual error of m sets ambiguity m*1
	return:
		-1 = invalid
	Note:
		revised from RTKLIB and Xiang Dong's NRTK
*/
int Ambiguity::Lambda(int n, int m,  math::matrix<double> a, math::matrix<double>  Q, 
									math::matrix<double> & F,	 double* s)
{
	if (n<=0||m<=0) return -1;

	int info;
	math::matrix<double>L(n,n);
	math::matrix<double>Z(n,n);
	math::matrix<double>E(m,n);//n*m in rtklib, E is zn in "search"
	math::matrix<double>ET;
	math::matrix<double>Q1(n,n);
	for(int i=0;i<n;i++)
	{
		Z(i,i)		=1.0;
		Q1(i,i)	=1.0;
	}

	double*	D	=new double[n];
	for (int i=0;i<n;i++)
	{
		D[i]=0.0;
	}
	math::matrix<double>	z(n,1);
	if ((info=LTDLFactorization(n,Q,L,D))) {
		/* lambda reduction		z=Z'*a,	Qz=Z'*Q*Z=L'*diag(D)*L     */
		Reduction(n,L,D,Z);
		z=Z*a;

		/* mlambda search */
		if (!(info=search(n,m,L,D,z,E,s))) 
		{
			//info=solve("T",Z,E,n,m,F); /* F=Z'\E */
			ET=~E;
			if(MatrixSovle(Z,ET,F,n,m)==false)
			{
				return -1;
			}
			//F=(!Z)*(~E);
		}
	}
	delete[] D;
	return info;
}
#include "StdAfx.h"
#include "LAMBDAVec.h"
#include <vector>
#include "math.h"
#include "MatrixVec.h"

CLAMBDA::CLAMBDA(void)
{
}


CLAMBDA::~CLAMBDA(void)
{
}
//
//bool CLAMBDA::LDLTFactorization(int n,std::vector< std::vector<double> > Q, std::vector< std::vector<double> >& L, std::vector<double> & D)
//{
//	int i,j,k;
//	std::vector< std::vector<double> > A(n),t(n);
//	for(int i=0;i<n;i++){A[i].resize(n);t[i].resize(n);}
//	for(i=0;i<n;i++)
//		for(j=0;j<n;j++)
//			A[i][j]=Q[i][j];
//	for(i=0;i<n;i++)
//	{
//		for(j=n-1;j>i;j--)
//			L[i][j]=0;
//		L[i][i]=1;
//	}
//	for(i=0;i<n;i++)
//	{
//		for(j=0;j<i;j++)
//		{   
//			t[i][j]=A[i][j];
//			for(k=0;k<j;k++)
//			{
//				t[i][j]-=t[i][k]*L[j][k];
//			}
//		}
//		for(j=0;j<i;j++)
//		{
//			L[i][j]=t[i][j]/D[j];
//		}
//		D[i]=A[i][i];//????
//		for(k=0;k<i;k++)
//			D[i]-=t[i][k]*L[i][k];
//		if(D[i]<0.0)
//			return false;
//	}
//
//	return true;
//}


/* LD factorization (Q=L*diag(D)*L') -----------------------------------------*/
//L为上三角矩阵，实际上为上三角*D*下三角
bool CLAMBDA::LTDLFactorization(int n, std::vector< std::vector<double> > Q, std::vector< std::vector<double> > &L, std::vector<double>& D)
{
	int i,j,k;
	double a;
	std::vector< std::vector<double> > A=Q;

	for (i=n-1;i>=0;i--) 
	{
		if ((D[i]=A[i][i])<=0.0) return false;
		a=sqrt(D[i]);
		for (j=0;j<=i;j++) L[j][i]=A[j][i]/a;
		for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[k][j]-=L[k][i]*L[j][i];
		for (j=0;j<=i;j++) L[j][i]/=L[i][i];
	}

	return true;
}
void CLAMBDA::gauss(int n, std::vector< std::vector<double> >& L, std::vector< std::vector<double> > & Z, int i, int j)
{
	int k,mu;
	if ((mu=(int)ROUND(L[j][i]))!=0) {
		for (k=i;k<n;k++) L[j][k]-=(double)mu*L[i][k];
		for (k=0;k<n;k++) 
		{
			Z[j][k]-=(double)mu*Z[i][k];
			if(fabs(Z[j][k])>1e5)
				Z[j][k]=1e5;
		}
	}
}

void CLAMBDA::perm(int n, std::vector< std::vector<double> > &L, std::vector<double>&D, int j, double del, std::vector< std::vector<double> >  &Z)
{
	int k;
	double eta,lam,a0,a1;

	eta=D[j]/del;
	lam=D[j+1]*L[j][j+1]/del;
	D[j]=eta*D[j+1]; D[j+1]=del;
	for (k=0;k<=j-1;k++) 
	{
		a0=L[k][j]; a1=L[k][j+1];
		L[k][j]=-L[j][j+1]*a0+a1;
		L[k][j+1]=eta*a0+lam*a1;
	}
	L[j][j+1]=lam;
	for (k=j+2;k<n;k++) SWAP(L[j][k],L[j+1][k]);
	for (k=0;k<n;k++) SWAP(Z[j][k],Z[j+1][k]);
}

void CLAMBDA::reduction(int n, std::vector< std::vector<double> >& L, std::vector<double> &D, std::vector< std::vector<double> > &Z)
{
	int i,j,k;
	double del;
	j=n-2; k=n-2;
	while (j>=0) 
	{
		if (j<=k) 
			for (i=j+1;i<n;i++) 
				gauss(n,L,Z,i,j);
		del=D[j]+L[j][j+1]*L[j][j+1]*D[j+1];
		if (del+1E-6<D[j+1]) 			/* compared considering numerical error */
		{ 
			perm(n,L,D,j,del,Z);
			k=j; 
			j=n-2;
		}
		else j--;
	}
}

int CLAMBDA::search(int n, int m, std::vector< std::vector<double> > L,  std::vector<double>  D,	std::vector<double>& zs,std::vector< std::vector<double> >& zn,  std::vector<double>& s)
{
	//
	//SYNTAX:
	//            ======================================
	//            | [zn,  s]=search( n, m, L, D, zs) |
	//            ======================================
	//      Search m best integer vectors zn and correspondent quadratic residual error was given by s
	//INPUTS:
	//      n: dimension of ambiguity vectors
	//      m: how many ambiguity vectors one want to output, when m=2, zn could given the best and the second best solution  
	//      L&D: Qzz=L'*D*L, where Qzz was co-variance matrices after Z transformation 
	//      zs: float  ambiguities
	//OUTPUT:
	//      zn: m integer ambiguity vectors
	//      s:   quadratic residual error of m sets ambiguity vectors
	//	revised by D. Xiang on 2014/05/06 in Tongji Univ
	//	email: xiangdongwhu@163.com
	//	School.of Surveying and Geoinformatics Engineering,Tongji Univ.,China
	//	Originally written by D. Xiang on 2012 in CASM, Beijing 
	//////////////////////////////////////////////////////////////////////////
	int i,j,k,c,nn=0,imax=0;
	double newdist,maxdist=1E99,y;
	std::vector< std::vector<double> > S(n); std::vector<double> dist(n),zb(n),z(n),step(n);
	for(i=0;i<n;i++)S[i].resize(n);
	CMatrix mymat;
	mymat.zeros(S,n,n);

	k=n-1; dist[k]=0.0;
	zb[k]=zs[k];
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
					S[i][k]=S[i][k+1]+(z[k+1]-zb[k+1])*L[i][k+1];
				zb[k]=zs[k]+S[k][k];
				z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
			}
			else {
				if (nn<m) 
				{
					if (nn==0||newdist>s[imax]) imax=nn;
					for (i=0;i<n;i++) zn[nn][i]=z[i];
					s[nn++]=newdist;
				}
				else 
				{
					if (newdist<s[imax]) 
					{
						for (i=0;i<n;i++) zn[imax][i]=z[i];
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
			for (k=0;k<n;k++) SWAP(zn[i][k],zn[j][k]);
		}
	}

	if (c>=LOOPMAX) {
//		AfxMessageBox("search loop count overflow");
		return -1;
	}
	return 0;
}

// revised from RTKLIB by T.TAKASU, Japan 
int CLAMBDA::lambda(int n, int m,  std::vector<double> a, std::vector< std::vector<double> > Q, std::vector< std::vector<double> >& F,	 std::vector<double>& s)
{
	//
	//SYNTAX:
	//            ======================================
	//            | [zn,  s]=search( n, m, L, D, zs) |
	//            ======================================
	//      Search m best integer vectors zn and correspondent quadratic residual error was given by s
	//INPUTS:
	//      n: dimension of ambiguity vectors
	//      m: how many ambiguity vectors one want to output, when m=2, zn could given the best and the second best solution  
	//      Q: (co)variance matrices of float ambiguities
	//      a: float  ambiguities
	//OUTPUT:
	//      F: m integer ambiguity vectors
	//      s:   quadratic residual error of m sets ambiguity vectors
	//	revised by D. Xiang on 2014/05/06 in Tongji Univ
	//	email: xiangdongwhu@163.com
	//	School.of Surveying and Geoinformatics Engineering,Tongji Univ.,China
	//	Originally written by D. Xiang on 2012 in CASM, Beijing 
	//////////////////////////////////////////////////////////////////////////
	if (n<=0||m<=0) return -1;
	if(F.size()==0)
	{
		F.resize(n);
		for(int i=0;i<n;i++)
			F[i].resize(m);
	}
	if(s.size()==0)
		s.resize(m);

	int info;
	CMatrix mymat;
	std::vector< std::vector<double> > L(n),Z,E(m),ET(n),ZT,Q1;
	std::vector< double > D(n),z(n);
	mymat.zeros(L,n,n); mymat.zeros(Z,n,n);
	mymat.zeros(L,n,n); mymat.eye(Z,n,n);ZT=Z; Q1=Z;
	for(int i=0; i<n; i++) { ET[i].resize(m);} for(int i=0;i<m;i++)E[i].resize(n);
	/* LD factorization 下*对*上*/
	if ((info=LTDLFactorization(n,Q,L,D))) {

		/* lambda reduction */
		reduction(n,L,D,Z);
		//matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z*a */
		mymat.multip(Z,a,z,n,n);
		/* mlambda search */
		if (!(info=search(n,m,L,D,z,E,s))) 
		{

			//info=solve("T",Z,E,n,m,F); /* F=Z'\E */
			mymat.MatT(E,ET,m,n);
			if(mymat.MatrixSovleDong(Z,ET,F,n,m)==false)
			{
//				AfxMessageBox("Error");
				return -1;
			}

		}
	}
	return info;
}

#pragma once
#include "StdAfx.h"
#include "MatrixVec.h"
#include <math.h>
using namespace std;
#include <vector>
#include <fstream>
#include <algorithm>
//#include "mclcppclass.h"
//#include "libcondition_num.h"
//#include "mclmcrrt.h"
//#include "libreg_par.h"
CMatrix::CMatrix(void)
{
}

CMatrix::~CMatrix(void)
{
}

inline	int max2(int a,int b)	{ return a<b?a:b; } 
inline	int min2(int a,int b)     	{ return a>b?a:b; } 

inline	int max3(int a,int b,int c)
{
		int max;
		max=a;
		if(max<b)
			max=b;
		if(max<c)
			max=c;
		return max;
}  
//construct two dimensional matrix
void CMatrix::built_mat(std::vector< std::vector<double> > &A, int row,int col)
{
	A.resize(row);
	for(int i=0;i<row;i++)
	{
		A[i].resize(col);
		for(int j=0;j<col;j++)
			A[i][j]=0.0;
	}
	
	return;
}

void CMatrix::built_vec(std::vector<double>& A,int row)
{
	A.resize(row);
	for (int i=0;i<row;i++) A[i]=0.0;

	return;
}

bool CMatrix::eye(std::vector< std::vector<double> >&  A,int row,int col)
{
	//generate a unit matrix with diagonal elements of one
	//it is similar with the one of Matlab

	int i,j;
	A.resize(row);
	for(int i=0;i<row;i++)
		A[i].resize(col);
	for(i=0;i<row;i++)
		for(j=0;j<col;j++)
		{
			A[i][j]=0;
			if(i==j) A[i][j]=1.0;
		}
	return true;
}
bool CMatrix::zeros(std::vector< std::vector<double> > & A,int row,int col)
{
	int i,j;
	A.resize(row);
	for(int i=0;i<row;i++)
		A[i].resize(col);
	for(i=0;i<row;i++)
		for(j=0;j<col;j++)
		{
			A[i][j]=0;
		}
	return true;

}

std::vector<int> find(std::vector<int> A,int m,int x)
{
	//SYNTAX:
	//   find same x from vector A and record their subscripts in vector B
	//
	//INPUTS:
	//    A: vector with m elements
	//    x: the searched element
	//OUTPUT:
	//    B: the subscripts of x in the vector A
	//
	//written by Bofeng Li on 2008/03/20 in Tongji Univ.
	//
	//email:Bofeng_Li@163.com
	//
	//Dept.of Surveying and Geo-informatics, Tongji Univ.
	//////////////////////////////////////////////////////////////////////////
	int i,j,num;
	std::vector<int> B;
	
	if (A.size()==0)
	{
		j=-1;
		return B;
	}
	else
	{
		num=0;

		for(i=0;i<m;i++)
		{
			if(A[i]==x) num=num+1;
		}
	}

	if(num==0) return B;

	B.resize(num);

	num=-1;
	for(i=0;i<m;i++)
	{
		if(A[i]==x) 
		{
			num=num+1;
            B[num]=i;
		}
	}
	return B;
}

bool CMatrix::plus(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  B,std::vector< std::vector<double> > & C,int row,int col)
{   
	//SYNTAX:
	//  Addition operation of two matrices with the same dimensions
	//
	//INPUTS:
	//  A and B: two matrices to be plused
	// row and col: dimensions
	//OUTPUTS:
	// C matrix
	int i,j;

	if(A.size()==0||B.size()==0||C.size()==0) return false;

	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
			C[i][j]=A[i][j]+B[i][j];
	}
		
	return true;
}

bool CMatrix::plus(std::vector< std::vector<double> > & A,std::vector< std::vector<double> >B,int row,int col)//A=A+B
{
	//SYNTAX:
	//  Addition operation of two matrices with the same dimensions
	//
	//INPUTS:
	//  A and B: two matrices to be plused
	// row and col: dimensions
	//OUTPUTS:
	// A matrix
	int i,j;

	if(A.size()==0||B.size()==0) return false;

	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
			A[i][j]=A[i][j]+B[i][j];
	}

	return true;
}


bool CMatrix::plusT(std::vector< std::vector<double> > & A,std::vector< std::vector<double> >B,int row,int col)//A=A+B'
{
	//SYNTAX:
	//  Addition operation of two matrices with the same dimensions
	//  row==col
	//INPUTS:
	//  A and B': two matrices to be plused
	// row and col: dimensions
	//OUTPUTS:
	// A matrix
	int i,j;

	if(A.size()==0||B.size()==0||row!=col) return false;

	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
			A[i][j]=A[i][j]+B[j][i];
	}

	return true;
}

bool CMatrix::plus(std::vector<double> A,std::vector<double> B,std::vector<double>& C,int row)
{   
	//addition of two vectors

	if(A.size()==0||B.size()==0||C.size()==0) return false;

	for(int i=0;i<row;i++)	C[i]=A[i]+B[i];

	return true;
}



bool CMatrix::subtr(std::vector<double> A,std::vector<double> B,std::vector<double>& C,int row)
{
	//substraction of two vectors

	if(A.size()==0||B.size()==0||C.size()==0) return false;

	for(int i=0;i<row;i++) C[i]=A[i]-B[i];
	
	return true;
}

bool  CMatrix::subtr(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  B,std::vector< std::vector<double> > & C,int row,int col)
{
	//subtraction of two matrices

	int i,j;

	if(A.size()==0||B.size()==0||C.size()==0) return false;

	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
			C[i][j]=A[i][j]-B[i][j];
	}
		
	return true;
}


bool CMatrix::subtr(std::vector< std::vector<double> > &A,std::vector< std::vector<double> >  B,int row,int col)
{
	if(A.size()==0||B.size()==0||A.size()!=B.size())
	return false;
	for(int i=0;i<row;i++)
	{
		if(A[i].size()!=B[i].size())
			return false;
		for(int j=0;j<col;j++)
			A[i][j]-=B[i][j];
	}
	return true;
}

bool CMatrix::multip(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  B,std::vector< std::vector<double> > & C,int M,int N,int L)
{
	//mutiply operation of two matrices
	
	int i,j,k;

	if(A.size()==0||B.size()==0||C.size()==0) return false;

	for(i=0;i<M;i++)
	{
		for(j=0;j<L;j++)
		{
			C[i][j]=0.0;
			for(k=0;k<N;k++) C[i][j]+=A[i][k]*B[k][j];
		}
	}
	return true;
}

bool CMatrix::multipT(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  B,std::vector< std::vector<double> > & C,int M,int N,int L) // C=A*B'
{
	//mutiply operation of two matrices
	
	int i,j,k;

	if(A.size()==0||B.size()==0||C.size()==0) return false;

	for(i=0;i<M;i++)
	{
		for(j=0;j<L;j++)
		{
			C[i][j]=0.0;
			for(k=0;k<N;k++) C[i][j]+=A[i][k]*B[j][k];
		}
	}
	return true;
}

bool CMatrix::multip(std::vector< std::vector<double> >  A, double s, std::vector< std::vector<double> > & C,int M,int N)
{
	int i,j;

	if(A.size()==0) return false;

	for(i=0;i<M;i++)
	{
		for(j=0;j<N;j++)
			C[i][j]=s*A[i][j];
	}
	return true;
}
bool CMatrix::multip(std::vector< std::vector<double> >  &A, double s,int M,int N)
{
		int i,j;

	if(A.size()==0) return false;

	for(i=0;i<M;i++)
	{
		for(j=0;j<N;j++)
			A[i][j]=s*A[i][j];
	}
	return true;
}
	

bool CMatrix::multip(std::vector< std::vector<double> >  A, std::vector<double>  B, std::vector<double>  & C,int M,int N)
{
	//mutiply operation of two matrices
	
	int i,k;

	if(A.size()==0||B.size()==0||C.size()==0) return false;

	for(i=0;i<M;i++)
	{
			C[i]=0.0;
			for(k=0;k<N;k++) C[i]+=A[i][k]*B[k];
	}
	return true;
}


bool CMatrix::multip(std::vector<double> A,std::vector<double>  B,std::vector< std::vector<double> > & C,int M,int L)
{
	//mutiply operation of two vectors to get a matrix

	int i,j;

	if(A.size()==0||B.size()==0||C.size()==0) return false;
	
	for(i=0;i<M;i++)
	{
		for(j=0;j<L;j++) C[i][j]=A[i]*B[j];
	}

	return true;		
}


double CMatrix::inner_prod(std::vector<double> A,int M)
{
	//inner product of a vector
	double Res=0.0;

//	if(A.size()==0) AfxMessageBox("Mutipvec");
	
	for(int i=0;i<M;i++) Res+=A[i]*A[i];
		
	return Res;
}

double CMatrix::inner_prod(std::vector< std::vector<double> > A, std::vector< std::vector<double> > B, int i, int j, int N ) //dot( A(:,i),B(:, j)
{
	double Res=0.0;

	//if(A.size()==0||B.size()==0) AfxMessageBox("Mutipvec");
	
	for(int k=0;k<N;k++) Res+=A[k][i]*B[k][j];
		
	return Res;
}

double norm(std::vector<double> A,int M)
{
	//inner product of a vector
	double Res=0.0;

	//if(A.size()==0) AfxMessageBox("Mutipvec");
	
	for(int i=0;i<M;i++) Res+=A[i]*A[i];
		
	return sqrt(Res);
}

bool CMatrix::MatrixInv(std::vector< std::vector<double> >& C, int N)
{
	int i, j, k;
	double S;

	if(C.size()==0) return false;
	for(i=0; i<N; i++)
	{
		for(j=i; j<N; j++)
		{
			S = C[i][j];
			for(k=0; k<i; k++) S = S - C[k][i] / C[k][k] * C[k][j];
			C[i][j] = S;
		}
	}	
	for(i=0; i<N; i++)
	{
		C[i][i] = 1.0 / C[i][i];
		for(j=i+1; j<N; j++)
		{
			S = 0.0;
			for(k=i; k<j; k++) S = S - C[i][k] * C[k][j];
			C[i][j] = S / C[j][j];
		}
	}
	
	for(i=0; i<N; i++)
	{
		for(j=i; j<N; j++)
		{
			S = 0;
			for(k=j; k<N; k++) S = S + C[i][k] / C[k][k] * C[j][k];
			C[j][i] = C[i][j] = S;
		}
	}
	return true;	
}


bool  CMatrix::MatrixInvDong( std::vector<std::vector<double> >& a,int n)
{
	double at;
	double bt;
	double am;
	int i,j,k;
	int tt;
	int N=2*n;
	std::vector<std::vector<double> > p;
	p.resize(n);
	for (i=0;i<n;i++)  p[i].resize(N);
	if(a.size()==0)	{		a.resize(n);for (i=0;i<n;i++) a[i].resize(n);	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			p[i][j]=a[i][j];
		}
		for(j=n;j<N;j++)
		{
			if(i==j-n)
				p[i][j]=1;
			else
				p[i][j]=0;
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
			if(at<1e-8)
			{
				//  AfxMessageBox("次矩阵不可逆")；
				return false;
			}
			am=1/p[k][k];
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
		for(j=0;j<n;j++)
			a[i][j]=p[i][j+n];
	return true;
}

bool CMatrix::MatrixInvDong( std::vector<std::vector<double> > a, std::vector<std::vector<double> >& b,int n)
{
	double at;
	double bt;
	double am;
	int i,j,k;
	int tt;
	int N=2*n;
	std::vector<std::vector<double> > p;

	p.resize(n);
	for (i=0;i<n;i++)  p[i].resize(N);
	if(b.size()==0) {b.resize(n); for (i=0;i<n;i++) b[i].resize(n); }
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			p[i][j]=a[i][j];
		}
		for(j=n;j<N;j++)
		{
			if(i==j-n)
				p[i][j]=1;
			else
				p[i][j]=0;
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
			if(at<1e-8)
			{
				//  AfxMessageBox("次矩阵不可逆")；
				return false;
			}
			am=1/p[k][k];
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
		for(j=0;j<n;j++)
			b[i][j]=p[i][j+n];
	return true;
}

bool CMatrix::MatrixSovleDong( std::vector<std::vector<double> > a, std::vector<std::vector<double> > b,std::vector<std::vector<double> > &c, int n, int m)
{
	double at;
	double bt;
	double am;
	int i,j,k;
	int tt;
	int N=n+m;
	std::vector< std::vector<double > > p(n);	
	for(i=0;i<n;i++)
		p[i].resize(N);
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			p[i][j]=a[i][j];
		}
		for(j=n;j<N;j++)
		{
			p[i][j]=b[i][j-n];
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
				//  AfxMessageBox("次矩阵不可逆")；
				return false;
			}
			am=1/p[k][k];
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
		for(j=0;j<m;j++)
			{c[i][j]=p[i][+j+n];}

		return true;

}
/////////////////////////////////////////////////
bool CMatrix::InvSymMat(std::vector< std::vector<double> > A,std::vector< std::vector<double> >& C,int N)
{
	//compute the inverse of a symmetric and definite matrix
	
	int i, j, k;
	double S;

	if(A.size()==0) return false;


    for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++) C[i][j]=A[i][j];
	}

	for(i=0; i<N; i++)
	{
		for(j=i; j<N; j++)
		{
			S = C[i][j];
			
			for(k=0; k<i; k++) S = S - C[k][i] / C[k][k] * C[k][j];
			
			C[i][j] = S;
		}
	}	

	for(i=0; i<N; i++)
	{
		C[i][i] = 1.0 / C[i][i];
		for(j=i+1; j<N; j++)
		{
			S = 0.0;
			for(k=i; k<j; k++) S = S - C[i][k] * C[k][j];
			C[i][j] = S / C[j][j];
		}
	}
	
	for(i=0; i<N; i++)
	{
		for(j=i; j<N; j++)
		{
			S = 0;
			for(k=j; k<N; k++) S = S + C[i][k] / C[k][k] * C[j][k];
			C[j][i] = C[i][j] = S;
		}
	}

	return true;	
}


//////////////////////////////////////////////
bool CMatrix::InvGuass(std::vector< std::vector<double> > A,std::vector< std::vector<double> >& B,int M)
{
	/*用Guass消去法对一般矩阵求逆///////////
	应用时A和B不能是同一个指针变量，这样的话得到的是变换后的A阵即为单位阵
	By Bofeng Li 06/25/2006
	Dept. of surveying and Geo-informatics of Tongji University in China*/

	if(A.size()==0||M<0)  return false;

	int i,j,k;
	double Pivot,**Tem;

	if(B.empty()==true)
	{
		B.resize(M);
		for(i=0;i<M;i++) B[i].resize(M);
	}
	

	Tem=new double *[M];

	for(i=0;i<M;i++)
	{
		Tem[i]=new double [M];

		for(j=0;j<M;j++)
		{
			Tem[i][j]=A[i][j];
			B[i][j]=0.0;
		}
		B[i][i]=1.0;
	}
	
	for(i=0;i<M;i++)
	{
		Pivot=Tem[i][i];

		if(Pivot==0)     //主元等于0，与下面不为0的行交换
		{
			for(j=i+1;j<M;j++)
			{
				if(Tem[j][i]!=0)
				{
					for(k=0;k<M;k++) 
					{
						Pivot=Tem[j][k];
						Tem[j][k]=Tem[i][k];
						Tem[i][k]=Pivot;
						Pivot=B[j][k];
						B[j][k]=B[i][k];
						B[i][k]=Pivot;
					}
					break;
				}
			}
		}

		Pivot=Tem[i][i];
		if(Pivot!=1.0)
		{
			for(k=0;k<M;k++)
			{
				Tem[i][k]/=Pivot;
				B[i][k]/=Pivot;
			}
		}
		for(j=i+1;j<M;j++)
		{
			Pivot=Tem[j][i];
			if(Pivot!=0)
			{
				for(k=0;k<M;k++)
				{
					Tem[j][k]-=Tem[i][k]*Pivot;
					B[j][k]-=B[i][k]*Pivot;
				}
			}
		}
	}

	for(i=M-1;i>-1;i--)
	{
		Pivot=Tem[i][i];
		
		if(Pivot==0)  //主元等于0，与下面不为0的行交换
		{
			for(j=i-1;j>-1;j--)
			{
				if(Tem[j][i]!=0)
				{
					for(k=0;k<M;k++) 
					{
						Pivot=Tem[j][k];
						Tem[j][k]=Tem[i][k];
						Tem[i][k]=Pivot;
						Pivot=B[j][k];
						B[j][k]=B[i][k];
						B[i][k]=Pivot;
					}
					break;
				}
			}
		}

		Pivot=Tem[i][i];

		if(Pivot!=1.0)
		{
			for(k=0;k<M;k++)
			{
				Tem[i][k]/=Pivot;
				B[i][k]/=Pivot;
			}
		}
		for(j=i-1;j>-1;j--)
		{
			Pivot=Tem[j][i];
			if(Pivot!=0)
			{
				for(k=0;k<M;k++)
				{
					Tem[j][k]-=Tem[i][k]*Pivot;
					B[j][k]-=B[i][k]*Pivot;
				}
			}
		}
	}

	for(i=0;i<M;i++) 
		if(Tem[i]!=NULL) delete Tem[i];
	
	if(Tem!=NULL) delete []Tem;
	
	return true;
}

/////////下三角矩阵求逆，用矩阵初等变换///////
bool CMatrix::InvGenMat(std::vector< std::vector<double> > A,std::vector< std::vector<double> >& B,int M)
{	
	if(A.size()==0||M<0) return false;
		
	int i,j,k;
	double Pivot,**Tem;

	if(B.empty()==true)
	{
		B.resize(M);
		for(i=0;i<M;i++) B[i].resize(M);
	}


	Tem=new double *[M];

	for(i=0;i<M;i++)
	{
		Tem[i]=new double [M];
		for(j=0;j<M;j++) Tem[i][j]=A[i][j];
	}
	
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++) B[i][j]=0.0;
		B[i][i]=1.0;
	}
	
	for(i=0;i<M;i++)
	{
		Pivot=Tem[i][i];
		if(Pivot!=1.0)
		{
			for(k=0;k<M;k++)
			{
				Tem[i][k]/=Pivot;
				B[i][k]/=Pivot;
			}
		}
		for(j=i+1;j<M;j++)
		{
			Pivot=Tem[j][i];
			if(Pivot!=0)
			{
				for(k=0;k<M;k++)
				{
					Tem[j][k]-=Tem[i][k]*Pivot;
					B[j][k]-=B[i][k]*Pivot;
				}
			}
		}
	}
	for(i=0;i<M;i++) 
		if(Tem[i]!=NULL) delete Tem[i];

	if(Tem!=NULL) delete []Tem;

	return true;
}

bool CMatrix::MatT(std::vector< std::vector<double> > A,std::vector< std::vector<double> >& B,int M,int N)
{
	//compute the transpose matrix of A
	//B=A';

	int i,j;

	if(A.size()==0) return false;

	if(B.size()==0)
	{
		B.resize(N);

		for(i=0;i<N;i++) B[i].resize(M);
	}

	for(i=0;i<M;i++)
	{
		for(j=0;j<N;j++) B[j][i]=A[i][j];
	}

	return true;
}

///////////////////////////////////////////////////////
//some special functions for surveying adjustments
bool CMatrix::form_ATPA(std::vector< std::vector<double> > A, std::vector< std::vector<double> > P,std::vector< std::vector<double> > & ATPA,int M,int N)
{
	//SYNTAX:
	//  compute the ATPA with unique weights; ATPA=A'*P*A;
    //
	//INPUTS:
	//      A: 2 dimensional matrix with m rows and n columns
	//      P: 2 dimensional matrix with m rows and m columns
	//m and n: give the dimensions
	//
	//OUTPUT:
	//   ATPA: =A'*P*A; with n rows and n columns
	//
	//written by Bofeng Li on 2005/11/30 in Tongji for BFGPS
	//
	//revised also by Bofeng Li on 2009/03/20 in Tongji for SRTK
	/////////////////////////////////////////////////////////////
	int i;

	if(A.size()==0||P.size()==0) return false;

	if(ATPA.size()==0)
	{
		ATPA.resize(N);
		for(i=0;i<N;i++) ATPA[i].resize(N);
	}

	//medium variables
	std::vector< std::vector<double> > AT,ATP;

	AT.resize(N);
	ATP.resize(N);
	for(i=0;i<N;i++) 
	{
		AT[i].resize(M);
		ATP[i].resize(M);
	}
	
	//compute transpose of A 
	if(!MatT(A,AT,M,N))  return false; 

	//compute ATP=AT*P
	if(!multip(AT,P,ATP,N,M,M)) return false;

	//compute the ATPA=AT*P*A=ATP*A
	if(!multip(ATP,A,ATPA,N,M,N)) return false;

	//release the memory of medium variables
	for(i=0;i<N;i++)
	{
		if (AT[i].empty()!=true) AT[i].clear();
		if (ATP[i].empty()!=true)  ATP[i].clear();
	}
	if (AT.empty()!=true) AT.clear();
	if (ATP.empty()!=true) ATP.clear();
	
	return true;
}

bool CMatrix::form_APAT(std::vector< std::vector<double> > A, std::vector< std::vector<double> > P,std::vector< std::vector<double> > &APAT,int M,int N)
{
	//SYNTAX:
	//  compute the APAT with unique weights; APAT=A*P*A';
    //  Amxn Pnxn APA'mxm
	//INPUTS:
	//      A: 2 dimensional matrix with m rows and n columns
	//      P: 2 dimensional matrix with m rows and m columns
	//m and n: give the dimensions
	//
	//OUTPUT:
	//   ATPA: =A'*P*A; with n rows and n columns
	//
	//written by Bofeng Li on 2005/11/30 in Tongji for BFGPS
	//
	//revised also by Bofeng Li on 2009/03/20 in Tongji for SRTK
	//re-written by D. Xiang on 2013/11/16 for Network-RTK
	/////////////////////////////////////////////////////////////
	int i;

	if(A.size()==0||P.size()==0) return false;

	if(APAT.size()==0)
	{
		APAT.resize(M);
		for(i=0;i<N;i++) APAT[i].resize(M);
	}

	//medium variables
	std::vector< std::vector<double> > AT,AP;

	AT.resize(N);
	AP.resize(M);
	for(i=0;i<N;i++) 
	{
		AT[i].resize(M);
	}
	for(i=0;i<M;i++) 
	{
		AP[i].resize(N);
	}
	//compute transpose of A 
	if(!MatT(A,AT,M,N))  return false; 

	//compute AP=AT*P
	if(!multip(A,P,AP,M,N,N)) return false;

	//compute the APAT=A*P*AT=A*P*AT
	if(!multip(AP,AT,APAT,M,N,M)) return false;

	//release the memory of medium variables
	//for(i=0;i<N;i++)
	//{
	//	if (AT[i].empty()!=true) AT[i].clear();
	//	if (ATP[i].empty()!=true)  ATP[i].clear();
	//}
	//if (AT.empty()!=true) AT.clear();
	//if (ATP.empty()!=true) ATP.clear();
	//
	return true;
}


bool CMatrix::form_ATPA(std::vector< std::vector<double> > A,std::vector< std::vector<double> > & ATPA,int M,int N)
{
	//SYNTAX:
	//  compute the ATPA with unique weights; ATPA=A'*A;
    //
	//INPUTS:
	//      A: 2 dimensional matrix with m rows and n columns
	//m and n: give the dimensions
	//
	//OUTPUT:
	//   ATPA: =A'*A;
	//
	//written by Bofeng Li on 2005/11/30 in Tongji for BFGPS
	//
	//revised also by Bofeng Li on 2009/03/20 in Tongji for SRTK
	/////////////////////////////////////////////////////////////
	int i,j,k;

	if(A.size()==0) return false;

	if(ATPA.size()==0)
	{
		ATPA.resize(N);
		for(i=0;i<N;i++) ATPA[i].resize(N);
	}

	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			ATPA[i][j]=0.0;

			for(k=0;k<M;k++) ATPA[i][j]+=A[k][i]*A[k][j];
		}
	}
	return true;
}

bool CMatrix::form_ATPL(std::vector< std::vector<double> > A,std::vector< std::vector<double> > P,std::vector<double> L,std::vector<double> &ATPL,int M,int N)
{
	//SYNTAX:
	//  compute the ATPL with unequal weights; ATPL=A'*P*L;
    //
	//INPUTS:
	//      A: 2 dimensional matrix with m rows and n columns
	//      P: 2 dimensional matrix with m rows and m columns
	//      L: a vector with m elements
	//m and n: give the dimensions
	//
	//OUTPUT:
	//   ATPL: =A'*P*L;
	//
	//written by Bofeng Li on 2005/11/30 in Tongji for BFGPS
	//
	//revised also by Bofeng Li on 2009/03/20 in Tongji for SRTK
	/////////////////////////////////////////////////////////////
	int i,j,k;

	if(A.size()==0||P.size()==0||L.size()==0) return false;

	if(ATPL.size()==0) ATPL.resize(N);

	//compute the ATP=A'*P

	std::vector< std::vector<double> > ATP;  //medium variable
	ATP.resize(N);

	for(i=0;i<N;i++) ATP[i].resize(M);
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<M;j++)
		{
			ATP[i][j]=0.0;

			for(k=0;k<M;k++) ATP[i][j]+=A[k][i]*P[k][j];
		}
	}
	
	//compute the ATPL=A'*P*L=ATP*L;
	for(i=0;i<N;i++)
	{
		ATPL[i]=0.0;

		for(j=0;j<M;j++) ATPL[i]+=ATP[i][j]*L[j];
	}
	

	//release the memory of medium variable
		for(i=0;i<N;i++)
	{
		if (ATP[i].empty()!=true)  ATP[i].clear();
	}
	if (ATP.empty()!=true) ATP.clear();
	return true;
}


bool CMatrix::form_ATPL(std::vector< std::vector<double> > A,std::vector<double> L,std::vector<double> & ATPL,int m,int n)
{
	//SYNTAX:
	//  compute the ATPL with unique weights; ATPL=A'*L;
    //
	//INPUTS:
	//      A: 2 dimensional matrix with m rows and n columns
	//      L: a vector with m elements
	//m and n: give the dimensions
	//
	//OUTPUT:
	//   ATPL: =A'*L;
	//
	//written by Bofeng Li on 2005/11/30 in Tongji for BFGPS
	//
	//revised also by Bofeng Li on 2009/03/20 in Tongji for SRTK
	/////////////////////////////////////////////////////////////
	int i,j;
	
	if(A.size()==0||L.size()==0) return false;
	
	if(ATPL.size()==0) ATPL.resize(n);

	for(i=0;i<n;i++)
	{
		ATPL[i]=0.0;

		for(j=0;j<m;j++) ATPL[i]+=A[j][i]*L[j];
	}
	
	return true;
}

bool CMatrix::form_ATPAL(std::vector< std::vector<double> > A,std::vector< std::vector<double> > P,std::vector<double> L,std::vector< std::vector<double> > & ATPA,std::vector<double> & ATPL,int M,int NN)
{
	//SYNTAX:
	//  compute the ATPA and ATPL simultaneously with unequal weights
	//         ATPA=A'*P*A; ATPL=A'*P*L;
    //
	//INPUTS:
	//      A: 2 dimensional matrix with m rows and n columns
	//      P: 2 dimensional matrix with m rows and m columns
	//      L: a vector with m elements
	//m and n: give the dimensions
	//
	//OUTPUT:
	//   ATPA: =A'*P*A; with n rows and n columns
	//   ATPL: =A'*P*L; with n elements
	//
	//written by Bofeng Li on 2005/11/30 in Tongji for BFGPS
	//
	//revised also by Bofeng Li on 2009/03/20 in Tongji for SRTK
	/////////////////////////////////////////////////////////////
	int i,j;

	if(A.size()==0||P.size()==0||L.size()==0) return false;

	if(ATPA.size()==0)
	{
		ATPA.resize(NN);
		for(i=0;i<NN;i++) ATPA[i].resize(NN);
	}
	if(ATPL.size()==0) ATPL.resize(NN);

	std::vector< std::vector<double> > AT,ATP;
	AT.resize(NN);	
	ATP.resize(NN);
	for(i=0;i<NN;i++) 
	{
		AT[i].resize(M);
		ATP[i].resize(M);
	}
	
	//compute transpose of A:=A'
	if(!MatT(A,AT,M,NN)) return false; 

	//compute the ATP=AT*P;
	if(!multip(AT,P,ATP,NN,M,M)) return false;

	//compute the ATPA=ATP*A;
	if(!multip(ATP,A,ATPA,NN,M,NN))return false;

	//compute the ATPL=ATP*L
	for(i=0;i<NN;i++)
	{
		ATPL[i]=0.0;
		for(j=0;j<M;j++) ATPL[i]+=ATP[i][j]*L[j];
	}

    //release the memory of medium variables
	
	for(i=0;i<NN;i++)
	{
		if (AT[i].empty()!=true) AT[i].clear();
		if (ATP[i].empty()!=true)  ATP[i].clear();
	}
	if (AT.empty()!=true) AT.clear();
	if (ATP.empty()!=true) ATP.clear();
	
	return true;
}

bool CMatrix::form_ATPAL(std::vector< std::vector<double> > A,std::vector<double> L,std::vector< std::vector<double> > & ATPA, std::vector<double> & ATPL,int M,int N)
{
	//SYNTAX:
	//  compute the ATPA and ATPL simultaneously with unique weights
	//         ATPA=A'*P*A; ATPL=A'*P*L;
    //
	//INPUTS:
	//      A: 2 dimensional matrix with m rows and n columns
	//      L: a vector with m elements
	//m and n: give the dimensions
	//
	//OUTPUT:
	//   ATPA: =A'*A; with n rows and n columns
	//   ATPL: =A'*L; with n elements
	//
	//written by Bofeng Li on 2005/11/30 in Tongji for BFGPS
	//
	//revised also by Bofeng Li on 2009/03/20 in Tongji for SRTK
	/////////////////////////////////////////////////////////////
	int i;

	if(A.size()==0||L.size()==0) return false;

	if(ATPA.size()==0)
	{
		ATPA.resize(N);

		for(i=0;i<N;i++) ATPA[i].resize(N);
	}

	if(ATPL.empty()==true) ATPL.resize(N);


	if(!form_ATPA(A,ATPA,M,N)) return false;

	if(!form_ATPL(A,L,ATPL,M,N)) return false;
	
	return true;
}

//V'PV=l'*P*l-W'*inv(NN)*W
bool CMatrix::lPl(std::vector< std::vector<double> > P,std::vector<double> L, double &lpl)
{
	if(L.size()==0||P.size()==0)		
		return false;
	std::vector< double > lp;
	for(int i=0; i<(int)L.size();i++)
	{
		lp.push_back(0.0);
		for(int j=0;j<(int)L.size();j++)
		{
			lp[i]+=L[j]*P[j][i];
		}
	}
	lpl=0.0;
	for(int i=0;i<(int)L.size();i++)
		lpl+=lp[i]*L[i];
	return true;
}

bool CMatrix::ATPAL_DD(std::vector< std::vector<double> > A,std::vector<double>L,std::vector< std::vector<double> >& ATPA,std::vector<double>& ATPL,int M,int N,int ISL1)
{
//为了减少计算量,重新编写计算atpa与atpl程序,根据双差权的特点
//功能：根据一个历元的A阵与L阵,还有观测值的个数（生成权阵）来计算ATPA与ATPL阵
//M——观测值个数
	double **ATP,SumACol;
	int i,j,k;
	double Weight;
	if(A.size()==0||ATPA.size()==0||L.size()==0||ATPL.size()==0) return false;
	
	if(ISL1==1) 
		Weight=1.0;  //L1 观测值
	else
		Weight=0.61; //L2观测值 (19.02/24.42)^2

	ATP=new double*[N];
	for(i=0;i<N;i++) 
	{	
		ATP[i]=new double[M];
		SumACol=0.0;
		for(j=0;j<M;j++) SumACol-=A[j][i];
		for(j=0;j<M;j++) ATP[i][j]=2.0*(SumACol/(M+1)+A[j][i])*Weight;
	}
	/////////compute matrixs ATPA and ATPL//////////////
	for(i=0;i<N;i++)
	{
		ATPL[i]=0.0;
		for(j=0;j<N;j++)
		{
			ATPA[i][j]=0.0;
			for(k=0;k<M;k++) ATPA[i][j]+=ATP[i][k]*A[k][j];
		}
		for(k=0;k<M;k++) ATPL[i]+=ATP[i][k]*L[k];
	}	
	
	for(i=0;i<N;i++) if (ATP[i]!=NULL) delete ATP[i];
	if (ATP!=NULL) delete []ATP;
	
	return true;
}

bool CMatrix::ATPAL_2DInter(std::vector< std::vector<double> > A,std::vector<double> L,std::vector< std::vector<double> > &ATPA,std::vector<double> &ATPL,int M,int N,int ii,int ISL1)
{
	//为了减少计算量,重新编写计算atpa与atpl程序,根据双差权的特点
	//功能：根据一个历元的A阵与L阵,还有观测值的个数（生成权阵）来计算ATPA与ATPL阵
	//M——观测值个数,ii--有效观测值个数，用来定权
	double **ATP,SumACol;
	double Weight;
	int i,j,k;
	if(A.size()==0||ATPA.size()==0||L.size()==0||ATPL.size()==0) return false;
	ATP=new double*[N];

	if(ISL1==1) 
		Weight=1.0;  //L1 观测值
	else
		Weight=0.61; //L2观测值 (19.02/24.42)^2

	for(i=0;i<N;i++) 
	{	
		ATP[i]=new double[M];
		SumACol=0.0;
		for(j=0;j<M;j++) SumACol-=A[j][i];
		for(j=0;j<M;j++) ATP[i][j]=2.0*(SumACol/(ii+1)+A[j][i])*Weight;
	}
	/////////compute matrixs ATPA and ATPL//////////////
	for(i=0;i<N;i++)
	{
		ATPL[i]=0.0;
		for(j=0;j<N;j++)
		{
			ATPA[i][j]=0.0;
			for(k=0;k<M;k++) ATPA[i][j]+=ATP[i][k]*A[k][j];
		}
		for(k=0;k<M;k++) ATPL[i]+=ATP[i][k]*L[k];
	}	
	
	for(i=0;i<N;i++) if (ATP[i]!=NULL) delete ATP[i];
	if (ATP!=NULL) delete []ATP;
	
	return true;
}
bool CMatrix::ATPAL_2DInter_Del(std::vector< std::vector<double> > A,std::vector<double> L,std::vector< std::vector<double> > &ATPA,std::vector<double> &ATPL,int M,int N,int Delete,int ii,int ISL1)
{
	//为了减少计算量,重新编写计算atpa与atpl程序,根据双差权的特点
	//功能：根据一个历元的A阵与L阵,还有观测值的个数（生成权阵）来计算ATPA与ATPL阵
	//M——观测值个数,N列数；在双差迭代时用，Delete表示删除观测值的下标
	//ii--有效观测值个数，包括此次剔出的观测值
	double **ATP,SumACol;
	double Weight;
	int i,j,k;
	if(A.size()==0||ATPA.size()==0||L.size()==0||ATPL.size()==0) return false;
	
	for(i=0;i<N;i++) A[Delete][i]=0.0;//删除该行
	L[Delete]=0.0;

	if(ISL1==1) 
		Weight=1.0;  //L1 观测值
	else
		Weight=0.61; //L2观测值 (19.02/24.42)^2

	ATP=new double*[N];
	for(i=0;i<N;i++) 
	{	
		ATP[i]=new double[M];
		SumACol=0.0;
		for(j=0;j<M;j++) SumACol-=A[j][i];
		for(j=0;j<M;j++) ATP[i][j]=2.0*(SumACol/(ii+1)+A[j][i])*Weight;
	}
	/////////compute matrixs ATPA and ATPL//////////////
	for(i=0;i<N;i++)
	{
		ATPL[i]=0.0;
		for(j=0;j<N;j++)
		{
			ATPA[i][j]=0.0;
			for(k=0;k<M;k++) ATPA[i][j]+=ATP[i][k]*A[k][j];
		}
		for(k=0;k<M;k++) ATPL[i]+=ATP[i][k]*L[k];
	}	
	
	for(i=0;i<N;i++) if (ATP[i]!=NULL) delete ATP[i];
	if (ATP!=NULL) delete []ATP;
	
	return true;
}

////compute the residuals//////////////////////////////
bool CMatrix::calc_resi(std::vector< std::vector<double> > A,std::vector<double> X,std::vector<double> L,std::vector<double> &V,int M,int N)
{   
	//A-M*N  Res--N*1  L--M*1  V--M*1
	//V=A*X-L

	int i,j;

	if(A.size()==0||X.size()==0||L.size()==0) return false;

	if(V.size()==0) V.resize(M);

	for(i=0;i<M;i++) 
	{
		V[i]=0.0;
		for(j=0;j<N;j++) V[i]+=A[i][j]*X[j];

		V[i]=V[i]-L[i];
	}

	return true;
}

double CMatrix::calc_vtpv( std::vector<double> V,std::vector< std::vector<double> > P,int M)
{
	//SYNTAX:
	//compute the weighted squared sum of residuals
	//vtpv=V'*P*V
	//

	int i,j;

	double Res;

	double *Temp;

	Res=0.0;

	if(V.size()==0||P.size()==0) return Res;

	Temp=new double [M];

	for(i=0;i<M;i++)
	{
		Temp[i]=0.0;
		for(j=0;j<M;j++) Temp[i]=Temp[i]+V[j]*P[j][i];
	}

	for (i=0;i<M;i++) Res=Res+Temp[i]*V[i];

	if(Temp!=NULL) 	delete []Temp;

	return Res;
}

//////////////////////////////////////////////////
bool CMatrix::del_row(std::vector< std::vector<double> >& A,int M,int N,int irow)
{
	//delete one row elements from matrix A
	
	if (irow<0 || irow>=M) return false; 

    int i,j;

	for(i=irow;i<M-1;i++)		
	{
		for(j=0;j<N;j++) A[i][j]=A[i+1][j];
	}

	if (A[M-1].empty()!=true)A[M-1].clear();  //delete the last row

	return true;
}

bool CMatrix::del_elem(std::vector<double> &A,int M,int irow)
{
	//delete one row elements from matrix A
	
	if (irow<0 || irow>=M) return false; 

	 std::vector<double> ::iterator iter=A.begin()+irow;
	A.erase(iter);
	A.resize(M-1);
		return true;
}

//////////////////////////////////////////
bool CMatrix::del_col(std::vector< std::vector<double> >& A,int M,int N,int icol)
{
	//delete one column elements from matrix A
	//删除矩一行数据  A--M*N  删除第col行,从0开始

	if (icol<0 || icol>=N) return false; 

	for( int i=0;i<M;i++)
	{
		std::vector<double>::iterator iter=A[i].begin()+icol;
	   A[i].erase(iter);
	}
	return true;
}

//////////////////////////////////////////////////
bool CMatrix::del_row_col(std::vector< std::vector<double> > &A,int M,int icol)
{
	//delete one column and one row elements from squared matrix A
	//删除矩一行数据  A--M*N  删除第col行,从0开始

		
	if (icol<0 || icol>=M) return false; 

	if(A.size()==0) return false;
	for(int i=0;i<(int)(A.size());i++)
	{
		std::vector<double> ::iterator iter=A[i].begin()+icol;
		A[i].erase(iter);
	}	
	std::vector< std::vector<double> >::iterator iter=A.begin()+icol;
	A.erase(iter);

	return true;
}

//
bool CMatrix::add_row_col(std::vector< std::vector<double> > &A, int M,int icol)
{
	int i;
	std::vector<double> zero;
	zero.resize(M);
	for(i=0;i<M;i++)
		zero[i]=0.0;
	A.insert( A.begin()+icol,zero);
	for(int j=0;j<M+1;j++)
		A[j].insert(A[j].begin()+icol,0.0);
	return true;
}

bool CMatrix::add_col(std::vector< std::vector<double> > &A, int M,int N, int icol)
{
	int j;
	if(icol>N)
		return false;
	for( j=0;j<M;j++)
		A[j].insert(A[j].begin()+icol,0.0);
	return true;
}
bool CMatrix::add_row(std::vector< std::vector<double> > &A, int M,int N, int irow)
{
	int i;
	std::vector<double> zero;
	zero.resize(M);
	for(i=0;i<M;i++)
		zero[i]=0.0;
	A.insert( A.begin()+irow,zero);
	return true;
}


double CMatrix::Mutipvec(std::vector<double> A,std::vector<double> B,int M)
{
	//向量内积
	double Res;
	int i;
	Res=0.0;

	if(A.size()!=0||B.size()!=0)
		for(i=0;i<M;i++) Res=Res+A[i]*B[i];
		
	return Res;
}
///////////////////////矩阵分解运算/////////////////////////
bool CMatrix::LU(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  U,int row)
{
//LU分解,如果A的顺序主子式不为0，则A存在唯一的分解A=LU
//其中L为单位下三角矩阵，U为可逆上三角矩阵（参阅计算方法P158）
	int i,j,k,r;
	double t,p;
	std::vector< std::vector<double> > L;
	L.resize(row);
	for(i=0;i<row;i++)
		L[i].resize(row);


	if(A.empty()==true||U.empty()==true) return false;
	for (k=0;k<row;k++)
	{
		for (j=k;j<row;j++)   //求u的值
		{
			if (k==0) 
				A[k][j]=A[k][j];
			else
			{
				t=0;
			    for (r=0;r<k-1;r++) t=t+A[k][r]*A[r][j];
			    A[k][j]=A[k][j]-t;
			}
		}
		for (i=k+1;i<row;i++)  //求l的值
		{
			if (k==0)
				A[i][k]=A[i][k]/A[k][k];
			else
			{
				p=0;
			    for (r=0;r<k-1;r++)	p=p+A[i][r]*A[r][k];
			    A[i][k]=(A[i][k]-p)/A[k][k];
			}
		}
	}
	for(i=0;i<row;i++) 
	{
		for(j=0;j<i+1;j++)
			L[i][j]=A[i][j];
		for(j=i;j<row;j++)
			U[i][j]=A[i][j];	
	}
	for (i=0;i<row;i++) L[i][i]=1.0;
	//bug!! L not recorded??
	return true;
}
bool CMatrix::LDR(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  L,std::vector<double> D,std::vector< std::vector<double> >  R,int row)
{
	//LDR分解,由LU分解可得到LDR分解，R为单位上三角矩阵，实质只要对U分解
	//注意当A为对称正定时，LDR成为LDLT,此函数是通过先对对称正定矩阵A进行cholesky分解

	int i,j,k,r;
	double t,p;

	if(A.size()==0||D.size()==0||R.size()==0) return false;

	for (k=0;k<row;k++)
	{
		for (j=k;j<row;j++)   //求u的值
		{
			if (k==0) 
				A[k][j]=A[k][j];
			else
			{
				t=0.0;
		    	for (r=0;r<k;r++) t+=A[k][r]*A[r][j];
		    	A[k][j]=A[k][j]-t;
			}
		}
		for (i=k+1;i<row;i++)  //求l的值
		{
			if (k==0)
				A[i][k]=A[i][k]/A[k][k];
			else
			{
				p=0;
		    	for (r=0;r<k;r++)	p+=A[i][r]*A[r][k];
		    	A[i][k]=(A[i][k]-p)/A[k][k];
			}
		}
	}
	for(i=0;i<row;i++) 
	{
		D[i]=A[i][i];
		for(j=0;j<row;j++) R[i][j]=L[i][j]=0.0;
		for(j=0;j<i+1;j++) L[i][j]=A[i][j];	
		for(j=i;j<row;j++) R[i][j]=A[i][j]/D[i];
	}
	for (i=0;i<row;i++) L[i][i]=1.0;

	return true;
}
bool CMatrix::LDLT(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  L,std::vector<double> D,int row)
{
	//LDR分解,由LU分解可得到LDR分解，R为单位上三角矩阵，实质只要对U分解
	//注意当A为对称正定时，LDR成为LDLT,此函数是通过先对对称正定矩阵A进行cholesky分解
	//为了节省变量空间，不保留R，重新写该程序，以便专业采用
	int i,j,k,r;
	double t,p;
	std::vector< std::vector<double> >  TempA;
	
	if(A.size()==0||D.size()==0) return false;
	TempA.resize(row);
	for(i=0;i<row;i++)
	{
		TempA[i].resize(row);
		for(j=0;j<row;j++) TempA[i][j]=A[i][j];
	}

	for (k=0;k<row;k++)
	{
		for (j=k;j<row;j++)   //求u的值
		{
			if (k==0) 
				TempA[k][j]=TempA[k][j];
			else
			{
				t=0.0;
				for (r=0;r<k;r++) t+=TempA[k][r]*TempA[r][j];
				TempA[k][j]=TempA[k][j]-t;
			}
		}
		for (i=k+1;i<row;i++)  //求l的值
		{
			if (k==0)
				TempA[i][k]=TempA[i][k]/TempA[k][k];
			else
			{
				p=0;
				for (r=0;r<k;r++)	p+=TempA[i][r]*TempA[r][k];
				TempA[i][k]=(TempA[i][k]-p)/TempA[k][k];
			}
		}
	}
	for(i=0;i<row;i++) 
	{
		D[i]=TempA[i][i];
		for(j=0;j<row;j++) L[i][j]=0.0;
		for(j=0;j<i+1;j++) L[i][j]=TempA[i][j];	
	}
	for (i=0;i<row;i++) L[i][i]=1.0;

	return true;
}

bool CMatrix::free_memo(std::vector< std::vector<double> >   &A,int M,int N)
{  
	//release the memory

	if(A.size()==0) return false;
	
	for(int i=0;i<M;i++) 
		if(A[i].size()!=0) A[i].clear();
	
	if(A.size()!=0) A.clear();

	return true;
}

bool CMatrix::Adjustment(std::vector< std::vector<double> > NN, std::vector<double > W, double lpl, int ObsNum, std::vector< double>& X, std::vector<std::vector< double> >& Dxx,double &segma0 )
{
	if((int)NN.size()==0	||	(int)W.size()==0	||	NN.size()!=W.size())
		return false;
	int Num;
	Num=W.size();
	MatrixInvDong(NN,Num);
	multip(NN,W,X,Num,Num);
	double vtpv;
	vtpv=Mutipvec(W,X,Num);
	vtpv=lpl-vtpv;
	segma0=sqrt(vtpv/(ObsNum-Num));
	for(int  i=0;i<Num;i++ ){	for(int j=0;j<Num;j++)    Dxx[i][j]=segma0*segma0*NN[i][j]; }
	return true;
}

bool CMatrix::FprintMatrix(std::vector< std::vector<double> > A, CString output)
{
	if(A.size()==0||output.GetLength()==0)
		return false;
	CString strPath=_T("");
	 CString cd=_T("");
     GetModuleFileName(NULL, strPath.GetBufferSetLength(MAX_PATH+1), MAX_PATH); 
     strPath.ReleaseBuffer(); 
     int  i= strPath.ReverseFind('\\');   
     strPath = strPath.Left(i);
	 i= strPath.ReverseFind('\\');   
	 strPath = strPath.Left(i);
     cd=strPath+_T("\\Debug\\");
	 fstream fout;
	 CString strout;
	 fout.open(cd+output,ios::out);
	 fout.precision(20);
	 for( std::vector<int>::size_type i=0;i<A.size();i++)
	 {
		 for(std::vector<int>::size_type j=0;j<A[i].size();j++)
			 fout<<A[i][j]<<"    ";
		 fout<<endl;
	 }
	 fout.close();
	 return true;
}

bool CMatrix::FprintMatrix(std::vector<double>  A, CString output)
{
	if(A.size()==0||output.GetLength()==0)
		return false;
	CString strPath=_T("");
	 CString cd=_T("");
     GetModuleFileName(NULL, strPath.GetBufferSetLength(MAX_PATH+1), MAX_PATH); 
     strPath.ReleaseBuffer(); 
     int  i= strPath.ReverseFind('\\');   
     strPath = strPath.Left(i);
	 i= strPath.ReverseFind('\\');   
	 strPath = strPath.Left(i);
     cd=strPath+_T("\\Debug\\");
	 fstream fout;
	 CString strout;
	 fout.open(cd+output,ios::out);
	 fout.precision(12);
	 for( std::vector<int>::size_type i=0;i<A.size();i++)
	 {
		 fout<<A[i];
		 fout<<endl;
	 } 
	 fout.close();
	 return true;
}


bool CMatrix::FprintMatrix(std::vector<double>  A, fstream& fout)
{
	if(A.size()==0)
		return false;
	if(fout.is_open()==false)
	{
		CString str;
		str.Format(_T("%s"),"Alarm: not can not open the file for vetorc input!  ");
		AfxMessageBox(str);
		return false;
	}
	for( std::vector<int>::size_type i=0;i<A.size();i++)
	 {
		 fout<<A[i];
		 fout<<endl;
	 } 
	 return true;
}

bool  CMatrix::FprintMatrix(std::vector< std::vector<double> > A, fstream& fout)
{
	if(A.size()==0)
		return false;
	if(!fout)
	{
		CString str;
		str.Format(_T("%s"),"Alarm: not can not open the file for vetorc input!  ");
		AfxMessageBox(str);
		return false;
	}
	for( std::vector<int>::size_type i=0;i<A.size();i++)
	 {
		 for(std::vector<int>::size_type j=0;j<A[i].size();j++)
			 fout<<A[i][j]<<"    ";
		 fout<<endl;
	 }
	return true;
}


bool CMatrix::Parition0( std::vector< std::vector<double> >  &NN, std::vector< std::vector<double> >  &UU,std::vector< std::vector<double> >  &N11,std::vector< std::vector<double> >  &N12,std::vector< std::vector<double> > &N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2)
{
	int m=N11.size();
	int n=N22.size();
	if(m==0||n==0||NN.size()!=m+n)return false;
	int i,j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
			N11[i][j]=NN[i][j];
		U1[i][0]=UU[i][0];
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			N22[i][j]=NN[i+m][j+m];
		U2[i][0]=UU[i+m][0];
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
			N12[i][j]=NN[i][j+m];
	}
	return true;

}
bool CMatrix::Partion00(std::vector< std::vector<double> >  &NN, std::vector< std::vector<double> > &UU,std::vector< std::vector<double> > N11,std::vector< std::vector<double> >  N12,std::vector< std::vector<double> > N22, std::vector< std::vector<double> >  U1, std::vector< std::vector<double> > U2)
{
	int m=N11.size();
	int n=N22.size();
	if(m==0||n==0||NN.size()!=m+n) return false;
	int i,j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
			NN[i][j]=N11[i][j];
			UU[i][0]=U1[i][0];
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			NN[i+m][j+m]=N22[i][j];
			UU[i+m][0]=U2[i][0];
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{	
			NN[i][j+m]=N12[i][j];
			NN[j+m][i]=N12[i][j];
		}
	}
	return true;

}

bool CMatrix::Parition0(std::vector< std::vector<double> >  &Q11, std::vector< std::vector<double> >  &Q12,std::vector< std::vector<double> >  &Q22,std::vector< std::vector<double> >  &X1,std::vector< std::vector<double> >  &X2,std::vector< std::vector<double> >  N11,std::vector< std::vector<double> >  N12,std::vector< std::vector<double> > N22, std::vector< std::vector<double> >  U1,std::vector< std::vector<double> > U2)
{
	int m=N11.size();
	int n=N22.size();
	if(m==0||n==0)return false;
	int i,j;
	std::vector< std::vector<double> >  NN(m+n),UU(m+n),X(m+n);
	for(i=0;i<m+n;i++)
	{
		NN[i].resize(m+n);
		UU[i].resize(1);
		X[i].resize(1);
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
			NN[i][j]=N11[i][j];
		UU[i][0]=U1[i][0];
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			NN[i+m][j+m]=N22[i][j];
		UU[i+m][0]=U2[i][0];
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			NN[i][j+m]=N12[i][j];
			NN[j+m][i]=NN[i][j+m];
		}
	}
	//Q11 Q12 Q22
	MatrixInvDong(NN,m+n);
	multip(NN,UU,X,m+n,m+n,1);
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
			Q11[i][j]=NN[i][j];
		X1[i][0]=X[i][0];
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			Q22[i][j]=NN[i+m][j+m];
		X2[i][0]=X[i+m][0];
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
			Q12[i][j]=NN[i][j+m];
	}
	return true;

}


bool CMatrix::Parition1(std::vector< std::vector<double> >  &N11,std::vector< std::vector<double> >  &N12,std::vector< std::vector<double> > &N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um)
{
	int m=N11.size();
	int n=N22.size();
	if(m==0||n==0)return false;

	int i,j,k;
	std::vector< std::vector<double> >  N21(n);
	std::vector< std::vector<double> >  Gain(n);
	for( i=0;i<n;i++)
	{
		Gain[i].resize(m);
		N21[i].resize(m);
	}
	std::vector< std::vector<double> >  N22_(n);	//N22'
	std::vector< std::vector<double> >  U2_(n);		//N22'
	for( i=0;i<n;i++)
	{
		N22_[i].resize(n);
		U2_[i].resize(1);
	}
	std::vector< std::vector<double> >  L;			//N11=L'*D*L : here, L as lower-triagular matrix
	std::vector< std::vector<double> >  D;
	//LDL
	LTDLFactorization(N11, L, D,m);

	MatrixInvDong(L,m);										//L=L^(-1)
																		//N21*N11^(-1)*N12=N21*L^(-1)*D^(-1)*(L')^-1*N12=Gain*D^-1*Gain'
	MatT(N12,N21,m,n);										//N21=N12'
	multip(N21,L,Gain,n,m,m);								// Gain=N21*L^(-1)
	//N22_
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			N22_[i][j]=0.0;
			for(k=0;k<m;k++) N22_[i][j]+=Gain[i][k]*Gain[j][k]/D[k][k];
		}
	}

	std::vector< double> temp(m);
	for(i=0;i<m;i++)//  U2_=Gain*[ D^(-1)* inv(L)' ]
	{
			temp[i]=0.0;
			for(j=0;j<m;j++) temp[i]+=L[j][i]*U1[j][0];
			temp[i]=temp[i]/D[i][i];
	}
	for(i=0;i<n;i++)
	{
			U2_[i][0]=0.0;
			for(j=0;j<m;j++) U2_[i][0]+=Gain[i][j]*temp[j];
	}
	
	// C=N22-N22_
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			C[i][j]+=(N22[i][j]-N22_[i][j]);
	}
	//Um=U2-U2_
	for(i=0;i<n;i++)
	{
		Um[i][0]+=(U2[i][0]-U2_[i][0]);
	}

	return true;
}


bool CMatrix::Parition11(std::vector< std::vector<double> >  &N11,std::vector< std::vector<double> >  &N12,std::vector< std::vector<double> > &N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um)
{
	int m=N11.size();
	int n=N22.size();
	if(m==0||n==0)return false;

	int i,j,k;
	std::vector< std::vector<double> >  N21(n);
	std::vector< std::vector<double> >  Gain(n);
	for( i=0;i<n;i++)
	{
		Gain[i].resize(m);
		N21[i].resize(m);
	}
	std::vector< std::vector<double> >  N22_(n);	//N22'
	std::vector< std::vector<double> >  U2_(n);		//N22'
	for( i=0;i<n;i++)
	{
		N22_[i].resize(n);
		U2_[i].resize(1);
	}
	std::vector< std::vector<double> >  L;			//N11=L'*D*L : here, L as lower-triagular matrix
	std::vector< std::vector<double> >  D;
	//LDL
	LTDLFactorization(N11, L, D,m);

	MatrixInvDong(L,m);										//L=L^(-1)
																		//N21*N11^(-1)*N12=N21*L^(-1)*D^(-1)*(L')^-1*N12=Gain*D^-1*Gain'
	MatT(N12,N21,m,n);										//N21=N12'
	multip(N21,L,Gain,n,m,m);								// Gain=N21*L^(-1)
	//N22_
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			N22_[i][j]=0.0;
			for(k=0;k<m;k++) N22_[i][j]+=Gain[i][k]*Gain[j][k]/D[k][k];
		}
	}

	std::vector< double> temp(m);
	for(i=0;i<m;i++)//  U2_=Gain*[ D^(-1)* inv(L)' ]
	{
			temp[i]=0.0;
			for(j=0;j<m;j++) temp[i]+=L[j][i]*U1[j][0];
			temp[i]=temp[i]/D[i][i];
	}
	for(i=0;i<n;i++)
	{
			U2_[i][0]=0.0;
			for(j=0;j<m;j++) U2_[i][0]+=Gain[i][j]*temp[j];
	}
	
	// C=N22-N22_
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			C[i][j]=(N22[i][j]-N22_[i][j]);
	}
	//Um=U2-U2_
	for(i=0;i<n;i++)
	{
		Um[i][0]=(U2[i][0]-U2_[i][0]);
	}

	return true;
}
//

bool CMatrix::Parition1(std::vector< std::vector<double> >  N11,std::vector< std::vector<double> >  N12,std::vector< std::vector<double> > N22, std::vector< std::vector<double> >  U1,std::vector< std::vector<double> >  U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um,std::vector< std::vector<double> > &Gain0, std::vector< std::vector<double> > &X, std::vector< std::vector<double> > &Y)
{
	int m=N11.size();
	int n=N22.size();
	if(m==0||n==0)return false;

	int i,j,k;
	std::vector< std::vector<double> >  N21(n);
	std::vector< std::vector<double> >  Gain(n);
	for( i=0;i<n;i++)
	{
		Gain[i].resize(m);
		N21[i].resize(m);
	}
	std::vector< std::vector<double> >  N22_(n);//N22'
	std::vector< std::vector<double> >  U2_(n);//N22'
	for( i=0;i<n;i++)
	{
		N22_[i].resize(n);
		U2_[i].resize(1);
	}
	std::vector< std::vector<double> >  L;			//N11=L'*D*L : here, L as lower-triagular matrix
	std::vector< std::vector<double> >  D;
	//LDL
	LTDLFactorization(N11, L, D,m);

	MatrixInvDong(L,m);										//L=L^(-1)
																		//N21*N11^(-1)*N12=N21*L^(-1)*D^(-1)*(L')^-1*N12=Gain*D^-1*Gain'
	MatT(N12,N21,m,n);										//N21=N12'
	multip(N21,L,Gain,n,m,m);								// Gain=N21*L^(-1)
	//N22_
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			N22_[i][j]=0.0;
			for(k=0;k<m;k++) N22_[i][j]+=Gain[i][k]*Gain[j][k]/D[k][k];
		}
	}

	std::vector< double> temp(m);
	for(i=0;i<m;i++)//  U2_=Gain*[ D^(-1)* inv(L)' ]
	{
			temp[i]=0.0;
			for(j=0;j<m;j++) temp[i]+=L[j][i]*U1[j][0];
			temp[i]=temp[i]/D[i][i];
	}
	for(i=0;i<n;i++)
	{
			U2_[i][0]=0.0;
			for(j=0;j<m;j++) U2_[i][0]+=Gain[i][j]*temp[j];
	}
	
	// C=N22-N22_
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			C[i][j]=(N22[i][j]-N22_[i][j]);
	}
	//Um=U2-U2_
	for(i=0;i<n;i++)
	{
		Um[i][0]=(U2[i][0]-U2_[i][0]);
	}
	
	//
	MatrixInvDong(N11,m);	
	MatrixInvDong( C, N22_, n);
	multip(N11,N12,Gain0,m,m,n);
	std::vector< std::vector< double > > U1_(m);
	for(i=0;i<m;i++)
		U1_[i].resize(1);
	multip(N22_,Um,Y,n,n,1);
	multip(N12,Y,U1_,m,n,1);
	subtr(U1,U1_,m,1);
	multip(N11,U1,X,m,m,1);
	return true;
}
//
bool CMatrix::Parition2(std::vector< std::vector<double> >  &Q11,std::vector< std::vector<double> >  &Q12,std::vector< std::vector<double> > &Q22, std::vector< std::vector<double> >  &X1,std::vector< std::vector<double> >  &X2, std::vector< std::vector<double> > & X1bar,std::vector< std::vector<double> > & Qhat,std::vector< std::vector<double> > & Xhat)
{
	int m=Q11.size();
	int n=Q22.size();
	if(m==0||n==0)return false;

	int i,j,k;
	std::vector< std::vector<double> >  Q21(n);
	std::vector< std::vector<double> >  Gain(n);
	for( i=0;i<n;i++)
	{
		Gain[i].resize(m);
		Q21[i].resize(m);
	}
	std::vector< std::vector<double> >  Q22_(n);//Q22'
	std::vector< std::vector<double> >  X2_(n);//
	for( i=0;i<n;i++)
	{
		Q22_[i].resize(n);
		X2_[i].resize(1);
	}
	std::vector< std::vector<double> >  L;			//N11=L'*D*L : here, L as lower-triagular matrix
	std::vector< std::vector<double> >  D;
	//LDL
	LTDLFactorization(Q11, L, D,m);
	MatrixInvDong(L,m);										//L=L^(-1)
																		//N21*N11^(-1)*N12=N21*L^(-1)*D^(-1)*(L')^-1*N12=Gain*D^-1*Gain'
	MatT(Q12,Q21,m,n);										//N21=N12'
	multip(Q21,L,Gain,n,m,m);								// Gain=N21*L^(-1)	
	//N22_
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			Q22_[i][j]=0.0;
			for(k=0;k<m;k++) Q22_[i][j]+=Gain[i][k]*Gain[j][k]/D[k][k];
		}
	}

	std::vector< double> temp(m);
	for(i=0;i<m;i++)//  U2_=Gain*[ D^(-1)* inv(L)' ]* (X1-X1bar )
	{
			temp[i]=0.0;
			for(j=0;j<m;j++) temp[i]+=L[j][i]*(X1[j][0]-X1bar[j][0]);
			temp[i]=temp[i]/D[i][i];
	}
	for(i=0;i<n;i++)
	{
			X2_[i][0]=0.0;
			for(j=0;j<m;j++) 
				X2_[i][0]+=Gain[i][j]*temp[j];
	}
	// C=N22-N22_
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			Qhat[i][j]=(Q22[i][j]-Q22_[i][j]);
	}
	//Um=U2-U2_
	for(i=0;i<n;i++)
	{
		Xhat[i][0]=(X2[i][0]-X2_[i][0]);
	} 
	return true;
}

//
bool CMatrix::Parition3(std::vector< std::vector<double> >  &Q11,std::vector< std::vector<double> >  &Q12,std::vector< std::vector<double> > &Q22, std::vector< std::vector<double> >  &X1,std::vector< std::vector<double> >  &X2,std::vector< std::vector<double> >& X2bar ,std::vector< std::vector<double> > & Qhat,std::vector< std::vector<double> > & Xhat)
{
	FprintMatrix(Q11,_T("Q11.txt"));
	FprintMatrix(Q12,_T("Q12.txt"));
	FprintMatrix(Q22,_T("Q22.txt"));
	//
	int m=Q11.size();
	int n=Q22.size();
	if(m==0||n==0)return false;

	int i,j,k;
	std::vector< std::vector<double> >  Q21(n);
	std::vector< std::vector<double> >  Gain(m);
	for( i=0;i<n;i++)
	{
		Q21[i].resize(m);
	}
	for(i=0;i<m;i++)
	{
		Gain[i].resize(n);
	}
	std::vector< std::vector<double> >  Q11_(m);//Q11'
	std::vector< std::vector<double> >  X1_(m);//
	for( i=0;i<m;i++)
	{
		Q11_[i].resize(m);
		X1_[i].resize(1);
	}
	std::vector< std::vector<double> >  L;			//N11=L'*D*L : here, L as lower-triagular matrix
	std::vector< std::vector<double> >  D;
	//LDL
	LTDLFactorization(Q22, L, D,n);
	MatrixInvDong(L,n);										//L=L^(-1)
																		//Q12*Q22^(-1)*Q21=Q12*L^(-1)*D^(-1)*(L')^-1*Q21=Gain*D^-1*Gain'
	MatT(Q12,Q21,m,n);										//Q21=Q12'
	multip(Q12,L,Gain,m,n,n);								// Gain=N21*L^(-1)	
	//Q11_
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
		{
			Q11_[i][j]=0.0;
			for(k=0;k<n;k++) Q11_[i][j]+=Gain[i][k]*Gain[j][k]/D[k][k];
		}
	}

	std::vector< double> temp(n);
	for(i=0;i<n;i++)//  X1_=Gain*[ D^(-1)* inv(L)' ]* (X1-X1bar )
	{
			temp[i]=0.0;
			for(j=0;j<n;j++) 
				temp[i]+=L[j][i]*(X2[j][0]-X2bar[j][0]);
			temp[i]=temp[i]/D[i][i];
	}
	for(i=0;i<m;i++)
	{
			X1_[i][0]=0.0;
			for(j=0;j<n;j++) 
				X1_[i][0]+=Gain[i][j]*temp[j];
	}
	// C=N22-N22_
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
			Qhat[i][j]=(Q11[i][j]-Q11_[i][j]);
	}
	//Xhat=X1-X1_
	for(i=0;i<m;i++)
	{
		Xhat[i][0]=(X1[i][0]-X1_[i][0]);
	} 
	
	return true;
}



bool CMatrix::PulsNormal(std::vector< std::vector<double> >  &N11,std::vector< std::vector<double> >  &N12,std::vector< std::vector<double> > &N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um)
{
	int m=N11.size();
	int n=N22.size();
	if(m==0||n==0)return false;
	int i,j;
	// C
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
			C[i][j]+=N11[i][j];
	}
	//Um
	for(i=0;i<m;i++)
	{
		Um[i][0]+=U1[i][0];
	}

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			C[i+m][j+m]+=N22[i][j];
	}
	//Um
	for(i=0;i<n;i++)
	{
		Um[i+m][0]+=U2[i][0];
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			C[i][j+m]+=N12[i][j];
		}
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
					C[j+m][i]=C[i][j+m];

		}
	}

	return true;
}

bool CMatrix::PulsNormal1(std::vector< std::vector<double> >  N11,std::vector< std::vector<double> >  N12,std::vector< std::vector<double> > N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um)
{
	int m=N11.size();
	int n=N22.size();
	if(m==0||n==0)
		return false;
	int i,j;
	// C
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
			C[i][j]=N11[i][j];
	}
	//Um
	for(i=0;i<m;i++)
	{
		Um[i][0]=U1[i][0];
	}

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			C[i+m][j+m]=N22[i][j];
	}
	//Um
	for(i=0;i<n;i++)
	{
		Um[i+m][0]=U2[i][0];
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			C[i][j+m]=N12[i][j];
		}
	}
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
					C[j+m][i]=C[i][j+m];

		}
	}
	return true;
}

/* LD factorization (Q=L*diag(D)*L') -----------------------------------------*/
// L为下三角矩阵，实际上为上三角*D*下三角
bool CMatrix::LTDLFactorization(std::vector< std::vector<double> > Q,std::vector< std::vector<double> > & L,std::vector< std::vector<double> >  &D,int n)
{
	int i,j,k;
    double a;
	std::vector< std::vector<double> > A(n);
	L.resize(n);
	D.resize(n);
	for(i=0;i<n;i++)
	{
		L[i].resize(n);
		D[i].resize(n);
		A[i].resize(n);
	}
	zeros(D,n,n);
	zeros(L,n,n);
	A=Q;

    for (i=n-1;i>=0;i--) 
	{
        if ((D[i][i]=A[i][i])<=0.0) return false;
        a=sqrt(D[i][i]);
        for (j=0;j<=i;j++) L[i][j]=A[i][j]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j][k]-=L[i][k]*L[i][j];
        for (j=0;j<=i;j++) L[i][j]/=L[i][i];
    }
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			if(i<j)
				L[i][j]=0.0;
			else if(i==j)
				L[i][j]=1.0;
	}
   return true;
}

void CMatrix::form_RefTMat(std::vector< std::vector<double> > &T, int AmbNum, int ip_ref)
{
	// ip_ref: the postiton where the last reference satellite located in the new prn list
	int i;
	if(T.size()==0)
	{
		T.resize(AmbNum+4);
		for(i=0;i<AmbNum+4;i++)
			T[i].resize(AmbNum+4);
	}
	eye(T,AmbNum+4,AmbNum+4);
	for(i=0;i<AmbNum;i++)
		T[i+4][ip_ref+4]=-1;

}

void  CMatrix::buble_sort(std::vector<int> list, std::vector< std::vector<double> > &T )
{
	//sort by the start epoch (bubble sort)
	int i,j,flag;
	int n=(int)list.size();
	int Help;

	if(T.size()==0)
	{
		T.resize(n);
		for(i=0;i<n;i++)
			T[i].resize(n);
		eye(T,n,n);
	}
	for(j=0;j<n-1;j++)
	{
		flag=0;
		for(i=n-1;i>j;i--)
		{
			if(list[i]<list[i-1])
			{
				Help=list[i];
				list[i]=list[i-1];
				list[i-1]=Help;
				flag=1;
				//swap i-th and j-th row of T
				for(int k=0;k<n;k++ )
				{
					double help=T[i][k];
					T[i][k]=T[i-1][k];
					T[i-1][k]=help;
				}
			}

		}

		if(flag==0) break;
	}
}

double CMatrix::trace(std::vector< std::vector< double > > Q)
{
	if(Q.size()==0) return 0.0;
	if(Q[0].size()!=Q.size()) return 0.0;
	double Val=0.0;
	for( int i=0;i<(int)Q.size();i++)
		Val+=Q[i][i];
	return Val;
}

double CMatrix::norm(std::vector< std::vector<double> > Q)
{	
	if(Q.size()==0) return 0.0;
	double Val=0.0;
	int m=(int)Q.size();
	int n=(int)Q[0].size();
	for (int i=0; i < m; i++)
		for (int j=0; j < n; j++)
			Val += Q[i][j] *Q[i][j];
	Val = sqrt( Val);
	return Val;
}

bool CMatrix::cond(std::vector< std::vector<double> > Q, int N, double &cond_num)
{
	std::vector< std::vector<double> > InvQ;
	MatrixInvDong(InvQ,N);
	double normQ=norm(Q);
	if(normQ==0.0) 
	{
		cond_num=0.0; 
		return false;}
	else{ cond_num=norm(InvQ)*normQ; 
	return true;}
}

inline void GetMR(std::vector< std::vector<double> > N11, std::vector< std::vector<double> > N12,std::vector< std::vector<double> > N22,std::vector< std::vector<double> >& MR, double alpha)
{
	std::vector< std::vector<double> > NR,Nq,N11_inv; // MSE 4x4 4x4
	NR.resize(4);Nq.resize(4);
	for(int i=0;i<4;i++) { NR[i].resize(4);Nq[i].resize(4); }
	for(int i=0;i<3;i++) for(int j=0;j<3;j++){ NR[i][j]=N11[i][j]; if(i==j) NR[i][j]+=alpha;}
	for(int i=0;i<3;i++){ NR[i][3]=N12[i][0]; NR[3][i]=N12[i][0];}
	NR[3][3]=N22[0][0];
	CMatrix mymat;
	mymat.MatrixInvDong(N11,N11_inv,3);
	for(int i=0;i<3;i++) for(int j=0;j<3;j++) Nq[i][j]=N11[i][j]+alpha*alpha*N11_inv[i][j]; 
	for(int i=0;i<3;i++){ Nq[i][3]=N12[i][0]; Nq[3][i]=N12[i][0];}
	Nq[3][3]=N22[0][0];
	mymat.MatrixInvDong(NR,4);
	mymat.form_APAT(NR,Nq,MR,4,4);
}
double CMatrix::Regu_par(std::vector< std::vector<double> > N11, std::vector< std::vector<double> > N12,std::vector< std::vector<double> > N22)
{
	double alpha;
	// search for alpha
	std::vector< std::vector<double> > MR;
	std::vector<double>trace_m;
	MR.resize(4);
	for(int i=0;i<4;i++) MR[i].resize(4);
	for( int i=0; i<100;i++)
	{
		alpha=0.02*i;
		GetMR(N11,  N12, N22, MR,  alpha);
		trace_m.push_back( trace(MR) );
	}
	std::vector<double>::iterator iter=std::min_element( trace_m.begin(), trace_m.end() );
	alpha=( iter-trace_m.begin() )*0.02;
	return alpha;
}

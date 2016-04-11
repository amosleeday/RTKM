#pragma once
#include "math.h"
#define LOOPMAX     10000           /* maximum count of search loop */
#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)
#include <vector>
/*------------------------------------------------------------------------------
* lambda.c : integer ambiguity resolution
*
*          Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*-----------------------------------------------------------------------------*/
class CLAMBDA
{
public:
	CLAMBDA(void);
	~CLAMBDA(void);
	// LDL' 分解 下三角*对角*上三角
	//bool LDLTFactorization(int n,std::vector< std::vector<double> > Q, std::vector< std::vector<double> >& L, std::vector<double> & D);
	// L'DL 分解 上三角*对角*下三角
	bool LTDLFactorization(int n,std::vector< std::vector<double> > Q, std::vector< std::vector<double> >& L,  std::vector<double> & D);
	void gauss(int n, std::vector< std::vector<double> >& L, std::vector< std::vector<double> > & Z, int i, int j);
	//交换算法
	void perm(int n, std::vector< std::vector<double> > &L, std::vector<double>&D, int j, double del, std::vector< std::vector<double> >  &Z);
	void reduction(int n, std::vector< std::vector<double> >& L, std::vector<double> &D, std::vector< std::vector<double> > &Z);
	int search(int n, int m, std::vector< std::vector<double> > L,  std::vector<double>  D,	std::vector<double>& zs,std::vector< std::vector<double> >& zn,  std::vector<double>& s);
	int lambda(int n, int m,  std::vector<double> a, std::vector< std::vector<double> > Q, std::vector< std::vector<double> >& F,	 std::vector<double>& s);
};
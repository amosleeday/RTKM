#include "stdafx.h"//changed
#include "VariableDef.h"
#include "MStruct.h"
#include<math.h>//calc
#include "matrix.h"
#include <fstream>
#include <string.h>




extern bool	  isOdd(int n);
extern int		  Fac(int n);
extern double LinearInterp(double x1,double x2,double y1,double y2,double x);
/*
flag	0(default):all steps are the same		1:steps are not 
here the size of x and y must be n+1
*/

extern double LagrangeInterp(int step,double* x,double* y,double xk,int flag,int n);
extern double Norm(double* a,int num);
extern double CMax(double* a,int num);
extern double AbsMax(double* a,int num);
extern int		  AbsIndMaxInd(double* a,int num);
extern double FreqSquareRatio(int sysid,int freqIndex);
extern double WaveLength(int sysid,int freqIndex);
extern int	  Prn2Sysid(int prn);
extern void	  RotationMatrix3(math::matrix<double>& R,double theta, int index);

extern CString Sys(int Sysid);
extern int Sysid(CString sys);
extern int Prn(CString sys,int prn);
extern double Freq(int prn,int ind);
extern void Freq(int sysid,double* freq);
extern double FreqSys(int sysid,int ind);
extern int Sysid(int prn);
extern double CombFreq(int sysid,int*coef);

extern int doy(int Year,int Month,int Day);
extern int doy(int sysid,int gpsw,double gpss );



extern void EarthConst(int ellipInd,double& a,double& eE);
extern void BLH2XYZ(double* BLH,int ellipInd,double* XYZ);
extern void XYZ2BLH(double* XYZ,int ellipInd,double* BLH);
extern math::matrix<double> XYZ2BLH(math::matrix<double> XYZ,int ellipInd);
extern void NEU2XYZ(double* NEU,double B0,double L0,double* XYZ);
extern void XYZ2NEU(double* XYZ,double B0,double L0,double* NEU);
extern math::matrix<double> XYZ2NEU(math::matrix<double> XYZ,double B0,double L0);
extern math::matrix<double> XYZ2NEU(double* coorBase,double* coorRover );
extern double getAzi(double dn,double de);
extern void	NEU2RAH(double *NEU,double* RAH);
extern void RAH2NEU(double* RAH,double* NEU);
extern void GaussProj(double B,double L,double L0,int ellipInd,double* xy);


extern double CombObs(int sysid,int* coef,double* threeobs);
extern void FreqRatio(int sysid,double* freq);

extern math::matrix<double> Kronecker( math::matrix<double>A,math::matrix<double> B,int flag);
extern void VcmatCom(math::matrix<int> coef,math::matrix<double>& Qyy,int flag);

extern math::matrix<double> ConvergeMat(int row,math::matrix<double>A,math::matrix<double>B);
extern math::matrix<double> ConvergeMat(int row,math::matrix<double>A,int colA,int colB);
extern double IonoFree(int sysid,double* obs,int* index);
extern double Mean(double *a,int num);
extern double Std(double* a,int num);

extern math::matrix<double> Cholesky(math::matrix<double>A,int n);
extern math::matrix<double> InvLowTri(math::matrix<double>L,int n);

extern math::matrix<double>MultiplyselfLowerUpper(math::matrix<double>A,int n);
extern math::matrix<double> MultiplyLowerOther( math::matrix<double>A,math::matrix<double>B,int n,int p );
extern math::matrix<double> MultiplyOtherLower( math::matrix<double>A,math::matrix<double>B,int m,int n );
extern math::matrix<double> MultiplyselfUpperLower( math::matrix<double>A,int n );

extern math::matrix<double>SolveLowerEquation(math::matrix<double>A,math::matrix<double>U,int n);
extern math::matrix<double>SolveUpperEquation(math::matrix<double>A,math::matrix<double>U,int n);
extern math::matrix<double>SolveUpperEquation2(math::matrix<double>A,math::matrix<double>U,int n);



extern void GetAscendPos(int* list1,int num,int* pos1);
extern int GetPos(int* list1,int a,int num);
 extern double Sum(double* s,int num);
extern math::matrix<double> DiagMatSym(math::matrix<double>A,int rowA,int rowB);
extern math::matrix<double> DiagMatSym(math::matrix<double>A,math::matrix<double>B);
extern double str2num(const char *s, int i, int n);
extern void Cstr2charPtr(CString kkk,char* code);

extern double dot(const double *a, const double *b, int n);
//extern void geph2pos(double ts, const BroadEphDataGlo geph, double *rs, double& dts);
//extern double geph2clk(double ts, const BroadEphDataGlo geph);
extern double ephclkerr(double ts,const BroadEphData eph);
extern void ephpos(math::matrix<double>& SatePos,BroadEphData eph,double t,double& dts);
extern void satepos(math::matrix<double>& SatePos,BroadEphData eph,double t,double& dts,math::matrix<double>& SateVel,double& ClkVel);
extern double timediff(double t1,double t2);
extern double FindLeapSec(const int year,const int mon);
extern double UTC2GPST(const int y,const int mon,double sec);
extern double geodistcorr(SatePos& sat,double* rec);
extern double geodistcorr(double* sat,double* rec);

//extern bool CheckUeserCtrl(ObsDataInfo datainfo,SppInfo userCtrl);

//CheckUeserCtrl is true then check this, if last step is false, use default
extern bool CheckStandCtrl(const SppCtrl userCtrl, SppCtrl& sppctrl); 
extern  void XYZ2RAH(double* XYZ,int sysid,double* XYZs,double& ele,double& azi);
extern void WeekSec(int& week,double& sec, YMDHMS currtime,int sysid);
extern double weightfactor(double ele,int index);
extern double IonoFreeGlo(double* obs,int* index,int k);
extern int FindPosInt(int* tarlist,int num,int tar);



extern void EnterBaseStnCoor(SppInfo& baseInfo,double* coor);
extern void EnterBaseStnCoor( SppInfoGlo& baseInfo,double* coor );
extern bool MatrixSovle( math::matrix<double> a,math::matrix<double> b,math::matrix<double> &c, int n, int m);
extern math::matrix<double>GetBlockMat(math::matrix<double>A,int rowFirst,int rowLast,int colFirst,int colLast,int flag);
extern void PutMat(math::matrix<double>&A,math::matrix<double>B,int rowFirst,int colFirst,int flag);
extern math::matrix<double> VecMat(int col,math::matrix<double>A,math::matrix<double>B);
extern math::matrix<double> VecMat(int col,math::matrix<double>A,int rowB);
extern math::matrix<double>MultSelfTrans(math::matrix<double> A,int row,int col);
extern math::matrix<double>MoveVecColEnd(math::matrix<double>A,int i);
extern math::matrix<double>MoveVecRowEnd(math::matrix<double>A,int i);
extern math::matrix<double>MoveVecRowColEnd(math::matrix<double>A,int i);
extern math::matrix<double> InsertZeroCol(math::matrix<double>A,int colA11,int colNumZero);
extern math::matrix<double> InsertZeroRow(math::matrix<double>A,int rowA11,int rowNumZero);
//extern math::matrix<double> InsertZeroRowCol(math::matrix<double>A,int rowA11,int colA11,int rowA22,int colA22,int rowcolNumZero);
extern math::matrix<double> InsertZeroRowCol(math::matrix<double>A,int rowA11,int rowcolNumZero);
extern math::matrix<double> InsertZeroRowCol(math::matrix<double>A,int k);
extern math::matrix<double>MultSelfTrans2(math::matrix<double> A,int row,int col);
extern math::matrix<double>CholeskyInv(math::matrix<double>A,int n);
extern math::matrix<double>CholeskyInv(math::matrix<double>A);

extern void ChangeRow(math::matrix<double>&A,int i,int j);
extern void ChangeCol(math::matrix<double>&A,int i,int j);
extern void ChangeRowCol(math::matrix<double>&A,int i,int j);
extern math::matrix<double> DiagMat(math::matrix<double>A,int rowA,int colA,math::matrix<double>B,int rowB,int colB);
extern void SetPtrToMatdiag(double* ptr,int num,math::matrix<double>&mat,int firstPos);
extern math::matrix<double>RemoveRowCol(math::matrix<double>A,int i);


extern double cofactor(double ele,int index);
extern void InitPtr(int * a,int num);
extern void InitPtr(double *a,int num);

extern void PrecEphPos(SatePos& presatpos, int satindex,double sec,int week,PrecEphData* preceph,int neph);
extern double GPST2BDT(double sec);
extern double geph2clk(double ts, const BroadEphDataGlo geph);
extern void geph2pos(double ts, const BroadEphDataGlo geph, double *rs, double& dts);
extern double GPST2UTC(const int y,const int mon,double sec);
extern double seph2clk(double time, const BroadEphDataSBAS seph);
extern void seph2pos(double time, const BroadEphDataSBAS seph, double *rs, double& dts);


extern math::matrix<double>EyeMat(int dim);
extern void SetMatClo(math::matrix<double>&A,int i,double num);
extern math::matrix<double>RemoveRow(math::matrix<double>A,int i);
extern math::matrix<double>RemoveCol(math::matrix<double>A,int i);
extern void ZeroMat(math::matrix<double>&A);
extern math::matrix<double> ZeroMat(int row,int col);
extern math::matrix<double>SQDotMat(math::matrix<double>A);
extern math::matrix<double>diagElemMat(math::matrix<double>A);
extern void xyz2plh(double *xyz,double* plh);



math::matrix<double> SortQaa(math::matrix<double> Qaa,math::matrix<double>& Qba,math::matrix<double>&ahat,DdAmbInfo ambinfo,DdAmbInfo& parambinfo,DdObsInfo obsinfo);
extern void OutPtr(double* ptr,int num);
extern void OutPtr(int* ptr,int num);
extern double traceMat(math::matrix<double> A);


extern void SortObsEpochData(ObsEpochData& data);

extern void ReadCommFile(DdCtrl& ddctrl,CString CommandFile);
extern void ReadCommFile(DdCtrl* ddctrl,InPutFileSet& inputfile,CString CommandFile,double* baseCrd,double* roverCrd);
extern void SppFileOut(fstream& fout,SppInfo sppinfo);
extern void SppFileOut(fstream& fout,SppInfoGlo sppinfo);
extern void ModifyPath(CString& filepath);
extern void BrowseFile(CString& strDir);

extern void PtrEqual(double* a,double* tar ,int num);
extern void PtrEqual(int* a,int* tar ,int num);
extern math::matrix<double> SelectQaa(math::matrix<double>Qaasort,math::matrix<double>Qbasort, math::matrix<double>& QbaPar,math::matrix<double>& ahatPar, DdObsInfo obsinfo,DdAmbInfo parambinfo,double maskele);

extern void GetOutPath(CString commfile,CString& filepath);

extern double FreqSysGlo(int ind,int n);


extern math::matrix<double>ElimRowCol(math::matrix<double>A,math::matrix<double>&U,int k);
extern void FileOutTropsIono(fstream& fout,math::matrix<double> IonoVector,math::matrix<double> TropVector,DdData curData,double* ptrMapCur);
//extern void FileOutTropsIono(fstream& fout,math::matrix<double> IonoVector,math::matrix<double> TropVector,DdData curData,double* ptrMapCur,int flag);
extern math::matrix<double>ElimRowColNEQ(math::matrix<double>A,math::matrix<double>&U,int k);
extern void FileOutTropsIonoFloat(fstream& fout,math::matrix<double> IonoVector,math::matrix<double> TropVector,DdData curData,double* ptrMapCur,int flag);
extern void FileOutTropsIonoFix(fstream& fout,math::matrix<double> b_updated,DdData curData,double* ptrMapCur,int flag);
extern double DistofVector(double* a,double* b,int num);
extern double MWCom(double* phs, double* cod,int index1,int index2,int sysid);
extern double CombObsCycle(int sysid,int* coef,double* threeobs);
extern int AbsIndMinInd(double* a,int num);
extern double LCD(int a,int b);
extern double LCM(int a,int b);
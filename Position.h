#pragma  once
#include "stdafx.h"//changed
#include "VariableDef.h"
#include "ExtFun.h"
#include "MStruct.h"
#include "Rinex.h"
#include "Orbit.h"
#include "Ambiguity.h"
#include "EmpiricalModel.h"

#include<iostream>
#include <iomanip>
#include<fstream>
/*
	This file contains the method of position
		1.Single Point Position
		2.Difference Position
			2.1 int baseline
			2.2 long baseline
		3.Relative Position
			3.1 int
			3.2 long
		4.Combination
*/



class Position
{
public:
	bool FindEph(BroadEphData* ephdata,int ephNum, int prn,double t,int& ephpos);
	bool FindEph(BroadEphDataGlo* ephdata,int ephNum, int prn,double t,int& ephpos);
	bool StandPosition(BroadEphHeader ephheader,BroadEphData* ephdata,ObsEpochData epochdata, int ephNum,SppCtrl sppctrl,SppInfo& sppinfo,ObsEpochData& lastData);
	bool StandPosition(BroadEphHeader ephheader,BroadEphDataGlo* ephdata,ObsEpochData epochdata,int ephNum,SppCtrl sppctrl,SppInfoGlo& sppinfo,ObsEpochData& lastData);
	bool SelectSateOnEle(double MaskEle,SppInfo& sppinfo,ObsEpochData& lastData);
	int SelectRefSate(SppInfoGlo sppinfobase,SppInfoGlo sppinforover,double maskEle,ObsEpochData lastDataBase,ObsEpochData lastDataRover,SdData& lastSdData,DdObsInfo& obsinfo,int* dNum);
	int SelectRefSate(SppInfo sppinfobase,SppInfo sppinforover,double maskEle,ObsEpochData lastDataBase,ObsEpochData lastDataRover,SdData& lastSdData,DdObsInfo& obsinfo,int is_initRTK,int& refPrnPre);
	void SDstation(int* pos1,int* pos2,int count,int refPrn, ObsEpochData zdbase,ObsEpochData zdrover,SppInfo baseinfo,SppInfo roverinfo, SdData& sddata);
	void SDstation(int* pos1,int* pos2,int count,int refPrn, ObsEpochData zdbase,ObsEpochData zdrover,SppInfoGlo baseinfo,SppInfoGlo roverinfo, SdData& sddata);
	void DoubleDiff(int refPrn,int Ind,SdData sddata, DdData& dddata);
	void FormDesMatPos(math::matrix<double>& DesMatPos,math::matrix<double>& L,DdData dddata,DdObsInfo ddobsinfo,DdCtrl ddctrl);
	void FormDesMatPos(math::matrix<double>& DesMatPos,math::matrix<double>& L,DdData dddata,int CodTypeNo,int PhsTypeNo);
	void FormDesMatTrop(math::matrix<double>& DesMatTrop,DdCtrl ddctrl,DdData dddata);
	void FormDesMatIono(math::matrix<double>& DesMatIono,DdCtrl ddctrl,DdData dddata);
	void FormDdErrEq(math::matrix<double>& DesMatPos,math::matrix<double>& DesMatTrop,math::matrix<double>& DesMatIono,math::matrix<double>& DesMatAmb, math::matrix<double>& Weight, math::matrix<double>& L,DdData dddata,DdCtrl ddctrl,DdAmbInfo ambinfo,DdObsInfo ddobsinfo);
	void FormDdErrEq(math::matrix<double>& DesMatPos,math::matrix<double>& DesMatTrop,math::matrix<double>& DesMatIono,math::matrix<double>& DesMatAmb, math::matrix<double>& L,DdData dddata,DdCtrl ddctrl,DdAmbInfo ambinfo,DdObsInfo ddobsinfo);
	void CalibRecClkErr(ObsEpochData& lastData,SppInfo sppinfo);
	void ReFormDdErrEq(math::matrix<double>& DesMatAmb, math::matrix<double>& L, DdAmbInfo ddambinfo,DdCtrl ddctrl,DdData curdddata);
	void ChangePreAmb(DdAmbInfo& preAmbInfo,int newRef,DdCtrl ddctrl);
	void ChangePredddata(DdData& predddata,int newRef);

	void FormResidual(math::matrix<double>& L,DdData dddata,DdCtrl ddctrl,DdAmbInfo ambinfo);
	void FormResidual(math::matrix<double>& L,math::matrix<double>& V,DdData dddata,DdCtrl ddctrl,DdAmbInfo ambinfo);
	void VcmatCom(math::matrix<double> coef,math::matrix<double>& Qyy,int flag);
	void NextEpoch(CStdioFile Ofile,CString line,ObsEpochData& epochData1,ObsHeader obsHeader,int nEpoch,int eventFlag,double simu,ObsEpochData epochData2,ProcessRinex readdata);
	void SateUpDown(int* PrePrnList,int*CurPrnList,int PreNum,int CurNum,int& upNum,int& downNum,int* upPrn,int* downPrn);
	void CalibSateClkErr(ObsEpochData& lastData,SppInfo sppinfo); 
	ObsEpochData GetSysData(int sysid,ObsEpochData epochdata);
	bool baseStn(SppInfo& baseStnInfo,BroadEphData* ephdata,ObsEpochData epochdata,ObsEpochData& lastData,int ephNum,SppCtrl baseCtrl);
	bool  baseStn(SppInfoGlo& baseStnInfo,BroadEphDataGlo* ephdata,ObsEpochData epochdata,ObsEpochData& lastData,int ephNum,SppCtrl baseCtrl);
	bool calcAllSatePos(SppInfo& baseStnInfo,BroadEphData* ephdata,ObsEpochData& lastData,int ephNum);
	bool calcAllSatePos(SppInfoGlo& baseStnInfo,BroadEphDataGlo* ephdata,ObsEpochData& lastData,int ephNum);
	bool SelectRefSate2(SppInfo sppinfo1,SppInfo sppinfo2,double maskEle,int& refPrn,ObsEpochData lastDataBase,ObsEpochData lastDataRover,SdData& lastSdData);
	int PositionMain( ObsEpochData BaseData, ObsEpochData RoverData, BroadEphHeader ephheader, BroadEphData* ephdata,int nSate,  SppCtrl sppctrl,SppInfo& sppinfo,SppInfo& baseInfo);
	void ComObsCod(DdCtrl& ddctrl,DdData& temp,DdData dddata,DdObsInfo& ddobsinfo);
//	void ComObsPhs(DdCtrl& ddctrl,DdData& temp,DdData dddata,DdObsInfo& ddobsinfo);
	void ComObsPhs(DdCtrl& ddctrl,DdData& temp,DdData dddata,DdObsInfo& ddobsinfo,DdAmbInfo& curambinfo);
	//DdData ComObsPhsCod(DdCtrl& ddctrl,DdObsInfo& ddobsinfo,DdData dddata);
	DdData ComObsPhsCod(DdCtrl& ddctrl,DdObsInfo& ddobsinfo,DdAmbInfo& curambinfo,DdData dddata);
	void ReFormConstWithAmb(math::matrix<double>& L, DdAmbInfo ambinfo,DdObsInfo ddobsinfo,DdCtrl ddctrl,DdData dddata);
	void  PassPreAmb(DdAmbInfo preamb,DdAmbInfo& curamb,DdCtrl ddctrl);
	void PassPreAmb(DdAmbInfo pre,DdAmbInfo& cur,int typeNum);
	math::matrix<double> FormWeight(DdCtrl ddctrl,DdData curdata,DdObsInfo obsinfo);
	void SetEleToObsinfo(DdObsInfo& obsinfo,SppInfo baseinfo,SppInfo roverinfo,int* pos1,int* pos2,int refprn,int count);
	void SetEleToObsinfo(DdObsInfo& obsinfo,SppInfoGlo baseinfo,SppInfoGlo roverinfo,int* pos1,int* pos2,int refprn,int count);
	bool SelectSateOnEleGlo(double MaskEle,SppInfoGlo& sppinfo,ObsEpochData& lastData);
	//void CycleSlipDetection(double thresGF,double thresMW , int index1,int index2, DdAmbInfo& ambinfoCur,DdData& curdata,DdData predata);
	void  CycleSlipDetection(double thresGF,double thresMW ,double thresWL, int index1,int index2, DdAmbInfo& ambinfoCur,DdData& curdata,DdData& predata, AmbData* mw,double& N1,double& N2);


	void Position:: CrossCode( int refpos,SdData lastSdData, math::matrix<double>&Qdd0,math::matrix<double>&Qdd1);
	void Position:: PostProcess0916(int numEpoch,SdData* sdData,int* satprn,DdCtrl ddctrl);
	void Position:: PostProcessTimecorr(int numEpoch,SdData* sdData,int* satprn,DdCtrl ddctrl);
	void Position:: PostProcessTimeCross(int numEpoch,SdData* sdData,int* satprn,DdCtrl ddctrl);
	void Position:: TimeCrossCode( int refpos,SdData lastSdData, math::matrix<double>&Qdd0,math::matrix<double>&Qdd1);
	void Position:: PostParSingle(int numEpoch,SdData* sdData,int* satprn,DdCtrl ddctrl);
	void DoubleSystemPosition(DdCtrl ddctrl1,DdCtrl ddctrl2,DdData dddata1,DdData dddata2,DdObsInfo obsinfo1,DdObsInfo obsinfo2 );
	math::matrix<double> FormWeightVc(DdCtrl ddctrl,DdData curdata,DdObsInfo obsinfo);
	void Position:: FormDesMatAmbGlo(math::matrix<double>& DesMatAmb,DdCtrl ddctrl,DdData dddata);
	//void FormDesMatIonoGlo(math::matrix<double>& DesMatIono,DdCtrl ddctrl,DdData dddata);
	void FormDesMatIonoGlo(math::matrix<double>& DesMatIono,DdCtrl ddctrl,DdData dddata,double* dNum);
	void ComObsPhsGlo(DdCtrl& ddctrl,DdData& temp,DdData dddata,DdObsInfo& ddobsinfo,DdAmbInfo& curambinfo);
	DdData ComObsPhsCodGlo(DdCtrl& ddctrl,DdObsInfo& ddobsinfo,DdAmbInfo& curambinfo,DdData dddata,SppInfoGlo sppinfo);
	void FormDesMatTrop_Interval(DdCtrl ddctrl,math::matrix<double>& DesMatTrop_Cur, double* ptrMapPre,double* ptrMapCur,DdData dddata_cur,DdData dddata_pre);
	void FormDesMatTrop_ReNew(DdCtrl ddctrl, double nEpoch,double& firstTime,double& curTime, ObsEpochData roverData, DdData dddataCur,DdData dddataPre,double* ptrMapCur,double* ptrMapPre, math::matrix<double> DesMatTrop, math::matrix<double>& DesMatTrop_Hour,double interval);
	void FormDesMatAmb(math::matrix<double>& DesMatAmb,DdCtrl ddctrl,DdData dddata,DdObsInfo ddobsinfo,DdAmbInfo ambinfo);
	void GetEWL(EwlData* ewl,DdData curData );
	void GetMW(AmbData* mw,DdData curData);
	void GetNL(fixinfo* infoWl,DdData Lcdata,AmbData* NL1);
	void outEwl(int nEpoch,DdData ts,EwlData* ewl,fstream& fout,fixinfo* info, double thres_a,int firstepoch,int epochend,DdCtrl ddctrl);
	void outMw(int nEpoch,DdData ts,AmbData* mw,fstream& fout,fixinfo* info, double thres_a,int firstepoch,int epochend,DdCtrl ddctrl);
	void outNL(AmbData* NL1,fixinfo* infoNL,fixinfo* infoWL,double thres_a);
};

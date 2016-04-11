#pragma once
#include "stdafx.h"//changed
#include "VariableDef.h"
#include "ExtFun.h"
#include "MStruct.h"



class ProcessRinex
{
public:
	//void ReadObsFile(CString OFileName,ObsHeader& obsHeader,ObsEpochData* obsepochdata,int& nEpoch,int eventflag);
	void ReadObsFile(CString OFileName,ObsHeader& obsHeader,ObsEpochData obsepochdata,int& nEpoch,int eventflag);
	void ReadObsHeader(CStdioFile& Ofile,CString& line, ObsHeader& obsHeader);
	void ReadObsType2(CStdioFile& Ofile,CString& line,int start,int lenStr,int lenSpace,int maxNumInLine,int dataNum,ObsHeader& obsHeader);
	void ReadObsType3(CStdioFile& Ofile,CString& line,int start,int lenStr,int lenSpace,int maxNumInLine,int dataNum,ObsHeader& obsHeader,int sysNumCurr);

void ReadObsDataEpoch2(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,ObsHeader obsHeader);
//process one satellite
void ReadObsData2(CStdioFile& Ofile,CString& line,ObsDataRecord& obsdatarecord,int dataNum,ObsHeader obsHeader);
bool ReadObsData2FirstLine(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,int eventflag);
void insertObsData2(int Currsysid,double data,int freqIndex,int LLI,int SNR,CString typeindex,ObsDataRecord& obsdatarecord);

void ReadObsDataEpoch3(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,ObsHeader obsHeader);
void ReadObsData3(CStdioFile& Ofile,CString& line,ObsDataRecord& obsdatarecord,ObsHeader obsHeader);
bool ReadObsData3FirstLine(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,int eventflag);
void insertObsData3(int Currsysid,double data ,int freqIndex,int LLI,int SNR,CString typeindex,ObsDataRecord& obsdatarecord);

//void ReadBroadEph(CString NFileName,BroadEphHeader& brdephheader,BroadEphData* broadephdata,int& nSate);
void ReadBroadEph(CString NFileName,BroadEphHeader& brdephheader,BroadEphData* broadephdata,BroadEphDataGlo* gloeph,int& nSate,int& nSateGlo);
void ReadBrdEphHeader(CStdioFile& Nfile,CString& line,BroadEphHeader& brdephheader );
//void ReadBrdEphData(CStdioFile& Nfile,CString& line,int sysid,double version,BroadEphData& brdephdata);
void ReadBrdEphData(CStdioFile& Nfile,CString& line,int sysid,double version,BroadEphData& brdephdata,BroadEphDataGlo& GloEph,int& flag);
void ReadBrdEphData(CStdioFile& Nfile,CString& line,int sysid,double version,BroadEphData& brdephdata,BroadEphDataGlo& GloEph,BroadEphDataSBAS& sbseph,int& flag);
//bool ReadObsFileOneEpoch(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,ObsHeader obsHeader,int& nEpoch,int eventflag);
void ReadObsFileOneEpoch(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,ObsHeader obsHeader,int& nEpoch,int& eventflag);
void GetEpochData(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,ObsHeader obsHeader,int& nEpoch,int& eventflag);
ObsEpochData GetEpochData(CStdioFile& Ofile,CString& line,ObsHeader obsHeader,int& nEpoch,int& eventflag);
};





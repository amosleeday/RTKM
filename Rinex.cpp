#pragma once
#include "stdafx.h"//changed

#include<math.h>//calc
//#include "matrix.h"
#include "Rinex.h"
#include "ExtFun.h"
/*
	ReadBroadEph & ReadObsFile are the main functions in this class

 */
//------------------------------------Obs       File--------------------------------------------------
//------------------------------------ 
//------------------------------------ 

void ProcessRinex::GetEpochData(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,ObsHeader obsHeader,int& nEpoch,int& eventflag)
{
	ObsEpochData	ptrepochdata;
	ReadObsFileOneEpoch( Ofile, line, ptrepochdata,obsHeader,nEpoch,eventflag);

		obsepochdata=ptrepochdata;//operator = avoid the difference between two epoch

		int ks=0;

}


ObsEpochData ProcessRinex::GetEpochData(CStdioFile& Ofile,CString& line,ObsHeader obsHeader,int& nEpoch,int& eventflag)
{
	ObsEpochData	ptrepochdata;
	ReadObsFileOneEpoch( Ofile, line, ptrepochdata,obsHeader,nEpoch,eventflag);

	 return ptrepochdata;//operator = avoid the difference between two epoch


}

void ProcessRinex::ReadObsFileOneEpoch(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,ObsHeader obsHeader,int& nEpoch,int& eventflag)
{
	
	/*
	 CStdioFile Ofile;
	 Ofile.Open(OFileName,CFile::modeRead);
	CString line;
	ReadObsHeader(Ofile,line,obsHeader);*/
	//
	
	if(Ofile.ReadString(line))
	{

		if (obsHeader.version<2.99)
		{
				
				if (ReadObsData2FirstLine(Ofile, line, obsepochdata,eventflag))//modified the one epoch to multi-epoch,change  index 0 to nEpoch, or declare the paremeter outside
				{
					//int dataNum=obsepochdata.sateNum;
					ReadObsDataEpoch2(Ofile,line,obsepochdata,obsHeader);
					//ObsEpochData temp=obsepochdata[nEpoch];
					nEpoch++;//nEpoch	=(nEpoch>=29): 29 ? nEpoch+1;
				}
			
		}
		else if(obsHeader.version>2.99)
		{
				if (ReadObsData3FirstLine(Ofile, line, obsepochdata,eventflag))
				{
					//int dataNum=obsepochdata.sateNum;
					ReadObsDataEpoch3(Ofile,line,obsepochdata, obsHeader);
					//ObsEpochData temp=obsepochdata[nEpoch];

					//te++;
					nEpoch++;
				}
		}
		
	}
    else
	{
		eventflag	=10;
		//return false;
	}

}
/*
 *read all epoch data from Ofile
 *
 */
void ProcessRinex::ReadObsFile(CString OFileName,ObsHeader& obsHeader,ObsEpochData obsepochdata,int& nEpoch,int eventflag)
{
	CStdioFile Ofile;
	Ofile.Open(OFileName,CFile::modeRead);
	CString line;
	ReadObsHeader(Ofile,line,obsHeader);
	if (obsHeader.version<2.99)
	{
		while(Ofile.ReadString(line))
		{
			if (ReadObsData2FirstLine(Ofile, line, obsepochdata,eventflag))//modified the one epoch to multi-epoch,change  index 0 to nEpoch, or declare the paremeter outside
			{
				int dataNum=obsepochdata.sateNum;
				ReadObsDataEpoch2(Ofile,line,obsepochdata,obsHeader);
				//ObsEpochData temp=obsepochdata[nEpoch];
				nEpoch++;//nEpoch	=(nEpoch>=29): 29 ? nEpoch+1;
			}
		}
		nEpoch++;
	}
	else if(obsHeader.version>2.99)
	{
		//int te=0;
		while(Ofile.ReadString(line))
		{
			if (ReadObsData3FirstLine(Ofile, line, obsepochdata,eventflag))
			{
				int dataNum=obsepochdata.sateNum;
				ReadObsDataEpoch3(Ofile,line,obsepochdata, obsHeader);
				//ObsEpochData temp=obsepochdata[nEpoch];
				
				//te++;
				nEpoch++;
			}
		}
		nEpoch++;
	}
	
	Ofile.Close();

	 
}

void ProcessRinex::ReadObsHeader(CStdioFile& Ofile,CString& line, ObsHeader& obsHeader)
{
	//CString line;	line=" ";
	int lenline;		lenline=0;
	int syscount; syscount=0;

	while( Ofile.ReadString(line)  )//reading Oheader
	{
		
		//CString HeaderLabel=line.Mid(60,lenline);
		//HeaderLabel=HeaderLabel.TrimRight();
		if (line.Find(_T("END OF HEADER"))!=-1) break;
		if (line.Find(_T("RINEX VERSION / TYPE"))!=-1)
			{
				obsHeader.version=_wtof( line.Mid(0,9) );
				obsHeader.sysid=Sysid(line.Mid(40,1));		
			}


		if (line.Find(_T("APPROX POSITION XYZ"))!=-1)
		{
			CString temp;
			for (int i=0;i<3;i++)
			{
				temp=line.Mid(i*14,14);
				if (temp.Trim()!="")
				{
					obsHeader.appPos[i]=_wtof(temp);
				}
			}
		}
		else if (line.Find(_T("ANTENNA: DELTA H/E/N"))!=-1)
		{
			CString temp;
			for (int i=0;i<3;i++)
			{
				temp=line.Mid(i*14,14);
				if (temp.Trim()!="")
				{
					obsHeader.HEN[i]=_wtof(temp);
				}
			}
		}
		else if (line.Find(_T("INTERVAL"))!=-1)
		{
			obsHeader.interval=_wtof(line.Mid(0,10));
		}
		else if (line.Find(_T("# / TYPES OF OBSERV"))!=-1)//2.x
		{
			int dataNum;
			dataNum=_wtoi(line.Mid(0,6));
			ReadObsType2(Ofile,line,10,2,4,9,dataNum,obsHeader);
		}		
		else if (line.Find(_T("SYS / # / OBS TYPES"))!=-1)//3.x
		{
			
			//*obsHeader.obsType3.sysNum+=1;
			int ind=Sysid(line.Mid(0,1));
			double dataNum;dataNum=_wtoi(line.Mid(3,3));
			obsHeader.obsType3[ind-1].allType[syscount].sysid=Sysid(line.Mid(0,1));
			obsHeader.obsType3[ind-1].typeNum=_wtoi(line.Mid(3,3));
			ReadObsType3(Ofile,line,7,3,1,13,dataNum,obsHeader,ind-1);
			syscount++;
		}		
	}
}
void ProcessRinex::ReadObsType2(CStdioFile& Ofile,CString& line,int start,int lenStr,int lenSpace,int maxNumInLine,int dataNum,ObsHeader& obsHeader)
{	
	//start			start	postion in line
	//lenstr			target's width
	//lenspace	
	//maxNumInLine    max num of type in a line
	//dataNum				num of data in real
	obsHeader.obsType2.typeNum=_wtoi(line.Mid(0,6));

	int lineNum;
	lineNum	=(dataNum%maxNumInLine==0)?(int)(dataNum/maxNumInLine):(int)(dataNum/maxNumInLine)+1;
	int lineCount=1;
	if (lineNum==1)
	{
		for(int i=0;i<dataNum;i++)
		{
			obsHeader.obsType2.allType[i].obsCode=line.Mid(start+i*(lenStr+lenSpace),1);
			obsHeader.obsType2.allType[i].freCode=_wtoi(line.Mid(start+i*(lenStr+lenSpace)+1,1));
		}
	}
	else if(lineNum>=2)
	{
		int typeCount=0;
		while(lineNum>lineCount)
		{
			
			for(int i=0;i<maxNumInLine;i++)
			{
				obsHeader.obsType2.allType[typeCount].obsCode=line.Mid(start+i*(lenStr+lenSpace),1);
				obsHeader.obsType2.allType[typeCount].freCode=_wtoi(line.Mid(start+i*(lenStr+lenSpace)+1,1));
				typeCount++;
				
			}
			lineCount++;
			Ofile.ReadString(line);
		}
		if(lineNum==lineCount)
			{
				for(int i=0;i<dataNum-(lineCount-1)*maxNumInLine;i++)
				{
					obsHeader.obsType2.allType[typeCount].obsCode=line.Mid(start+i*(lenStr+lenSpace),1);
					obsHeader.obsType2.allType[typeCount].freCode=_wtoi(line.Mid(start+i*(lenStr+lenSpace)+1,1));
					typeCount++;
				}
			}
	}
	 
}
void ProcessRinex::ReadObsType3(CStdioFile& Ofile,CString& line,int start,int lenStr,int lenSpace,int maxNumInLine,int dataNum,ObsHeader& obsHeader,int sysNumCurr)
{	
	//start			start	postion in line
	//lenstr			target's width
	//lenspace	
	//maxNumInLine    max num of type in a line
	//dataNum				num of data in real
	

	int lineNum;
	lineNum	=(dataNum%maxNumInLine==0)?(int)(dataNum/maxNumInLine):(int)(dataNum/maxNumInLine)+1;
	int lineCount=1;
	if (line.Mid(0,1)=="G")
	{
		obsHeader.obsType3[sysNumCurr].sysid	=1;
	}
	else if (line.Mid(0,1)=="R")
	{
		obsHeader.obsType3[sysNumCurr].sysid	=2;
	}
	if (line.Mid(0,1)=="E")
	{
		obsHeader.obsType3[sysNumCurr].sysid	=3;
	}
	if (line.Mid(0,1)=="S")
	{
		obsHeader.obsType3[sysNumCurr].sysid	=4;
	}
	if (line.Mid(0,1)=="C")
	{
		obsHeader.obsType3[sysNumCurr].sysid	=5;
	}
	if (lineNum==1)
	{
		for(int i=0;i<dataNum;i++)
		{
			obsHeader.obsType3[sysNumCurr].allType[i].obsCode=line.Mid(start+i*(lenStr+lenSpace),1);
			obsHeader.obsType3[sysNumCurr].allType[i].freCode=_wtoi(line.Mid(start+i*(lenStr+lenSpace)+1,1));
			obsHeader.obsType3[sysNumCurr].allType[i].channel=line.Mid(start+i*(lenStr+lenSpace)+2,1);
		}
	}
	else if(lineNum>=2)
	{
		int typeCount=0;
		while(lineNum>lineCount)
		{
			
			for(int i=0;i<maxNumInLine;i++)
			{
				obsHeader.obsType3[sysNumCurr].allType[typeCount].obsCode=line.Mid(start+i*(lenStr+lenSpace),1);
				obsHeader.obsType3[sysNumCurr].allType[typeCount].freCode=_wtoi(line.Mid(start+i*(lenStr+lenSpace)+1,1));
				obsHeader.obsType3[sysNumCurr].allType[typeCount].channel=line.Mid(start+i*(lenStr+lenSpace)+2,1);
				typeCount++;
				
			}
			lineCount++;
			Ofile.ReadString(line);
		}
		if(lineNum==lineCount)
			{
				for(int i=0;i<dataNum-(lineCount-1)*maxNumInLine;i++)
				{
					obsHeader.obsType3[sysNumCurr].allType[typeCount].obsCode=line.Mid(start+i*(lenStr+lenSpace),1);
					obsHeader.obsType3[sysNumCurr].allType[typeCount].freCode=_wtoi(line.Mid(start+i*(lenStr+lenSpace)+1,1));
					obsHeader.obsType3[sysNumCurr].allType[typeCount].channel=line.Mid(start+i*(lenStr+lenSpace)+2,1);
					typeCount++;
				}
			}
	}
	 
}

//process one epoch
void ProcessRinex::ReadObsDataEpoch2(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,ObsHeader obsHeader)
{
	int numSate=obsepochdata.sateNum;
	int dataNum=obsHeader.obsType2.typeNum;
	for(int i=0;i<numSate;i++)
	{
		ReadObsData2( Ofile, line, obsepochdata.obsdatarecord[i], dataNum, obsHeader);
	}
}
//process one satellite
void ProcessRinex::ReadObsData2(CStdioFile& Ofile,CString& line,ObsDataRecord& obsdatarecord,int dataNum,ObsHeader obsHeader)
{
	int lineNum; 
	lineNum	=(dataNum%5==0)?(int)(dataNum/5):(int)(dataNum/5)+1;
	CString all;
	all="";
	for(int j=0;j<lineNum;j++)
	{
		Ofile.ReadString(line);
		int len=80-line.GetLength();

		for(int i=0;i<len;i++)
		{
			line+=_T(" ");
		}
		all+=line;
	}
	//int kk=all.GetLength();
	for(int k=0;k<dataNum;k++)
	{
		double data=_wtof(all.Mid(k*16,14));
		int LLI=(int)_wtoi(all.Mid(k*16+15,1));
		int SNR=(int)_wtoi(all.Mid(k*16+16,1));
		CString typeindex=obsHeader.obsType2.allType[k].obsCode;
		int     freqindex=obsHeader.obsType2.allType[k].freCode;
		int     currSysid=obsdatarecord.sysid;
		insertObsData2(currSysid,data,freqindex,LLI,SNR,typeindex, obsdatarecord);
	}

}//,ObsHeader obsHeader
bool ProcessRinex::ReadObsData2FirstLine(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,int eventflag)
{	
	int y;
	y=_wtoi(line.Mid(1,2));
	if (y>79)
	{
		y=y+1900;
	}
	else 
	{
		y=y+2000;
	}
	int skip=_wtoi(line.Mid(30,2));
	int flag=_wtoi(line.Mid(28,1));
			eventflag=flag;
	if(flag!=0 && flag!=6)
	{
		for(int i=0;i<skip;i++)
		{
			Ofile.ReadString(line);
		}

		return false;
	}
	YMDHMS obsTime;
	Gtime		protime;
	obsTime.year=y;
	obsTime.mon=_wtoi(line.Mid(4,2));
	obsTime.day=_wtoi(line.Mid(7,2));
	obsTime.hour=_wtoi(line.Mid(10,2));
	obsTime.min=_wtoi(line.Mid(13,2));
	obsTime.sec=_wtoi(line.Mid(16,11));
	WeekSec(obsepochdata.week,obsepochdata.sec,obsTime,1);//1 is GPST

	obsepochdata.flag=_wtoi(line.Mid(28,1));
	obsepochdata.sateNum=_wtoi(line.Mid(29,3));
	int sateNum=_wtoi(line.Mid(29,3));
	int lineNum;
	lineNum	=(sateNum%12==0)?(int)(sateNum/12):(int)(sateNum/12)+1;
	int lineCount=1;
	if (lineNum==1)
	{
		for(int i=0;i<sateNum;i++)
		{
			obsepochdata.obsdatarecord[i].sysid=Sysid(line.Mid(32+i*3,1));
			int temprn=(int)_wtoi(line.Mid(32+i*3+1,2));
			obsepochdata.obsdatarecord[i].PRN=Prn(line.Mid(32+i*3,1),temprn);
		}
	}
	else if(lineNum>=2)
	{
		int typeCount=0;
		while(lineNum>lineCount)
		{
			
			for(int i=0;i<12;i++)
			{
				obsepochdata.obsdatarecord[typeCount].sysid=Sysid(line.Mid(32+i*3,1));
				int temprn=(int)_wtoi(line.Mid(32+i*3+1,2));
				obsepochdata.obsdatarecord[typeCount].PRN=Prn(line.Mid(32+i*3,1),temprn);
				typeCount++;
				
			}
			lineCount++;
			Ofile.ReadString(line);
		}
		if(lineNum==lineCount)
			{
				for(int i=0;i<sateNum-(lineCount-1)*12;i++)
				{
					obsepochdata.obsdatarecord[typeCount].sysid=Sysid(line.Mid(32+i*3,1));
					int temprn=(int)_wtoi(line.Mid(32+i*3+1,2));
					obsepochdata.obsdatarecord[typeCount].PRN=Prn(line.Mid(32+i*3,1),temprn);
					typeCount++;
				}
			}
	}
	return true; 
}
void ProcessRinex::insertObsData2(int Currsysid,double data,int freqIndex,int LLI,int SNR,CString typeindex,ObsDataRecord& obsdatarecord)
{
		if( freqIndex==7)
		{
			freqIndex-=5;
		}
		else if ( freqIndex==6)  //E5a=L5 and E6=L6 can't exist at the same time
		{									//for Galileo
			freqIndex-=3;
		}
		else if(freqIndex==5)
		{
			freqIndex-=2;
		}
	if(typeindex=='C' || typeindex=='P')
	{
		obsdatarecord.PsRange[freqIndex-1]=data;//for 2.x , 3.x is different
		if (data>0.0)
		{
			obsdatarecord.numVadCod++;
			obsdatarecord.vadFlgCod[freqIndex-1]=1;
		}
		//if (Currsysid==3&&freqIndex==5)
		//{
		//}
		//else if (Currsysid==3&&freqIndex==8)
		//{
		//}
	}

	if(typeindex=='L')
	{
		obsdatarecord.Phase[freqIndex-1]=data;//for 2.x , 3.x is different
		if (data>0.0)
		{
			obsdatarecord.numVadPhs++;
			obsdatarecord.vadFlgPhs[freqIndex-1]=1;
		}
		obsdatarecord.LLI[freqIndex-1]=LLI;
		obsdatarecord.SNR[freqIndex-1]=SNR;
	}
	if(typeindex=='D')
	{
		obsdatarecord.Doppler[freqIndex-1]=data;//for 2.x , 3.x is different
	}
}

void ProcessRinex::ReadObsDataEpoch3(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,ObsHeader obsHeader)
{
	int numSate=obsepochdata.sateNum;
	for(int i=0;i<numSate;i++)
	{
		ReadObsData3( Ofile, line, obsepochdata.obsdatarecord[i], obsHeader);
	}
	
}
void ProcessRinex::ReadObsData3(CStdioFile& Ofile,CString& line,ObsDataRecord& obsdatarecord,ObsHeader obsHeader)
{
		Ofile.ReadString(line);
		
		CString sys=line.Mid(0,1);
		int    Currsysid=Sysid(sys)-1;
		int		temprn=(int)_wtoi(line.Mid(1,2));
		//obsdatarecord.PRN=Prn(sys,temprn);
		obsdatarecord.PRN=Prn(sys,temprn);
		int		typeNum=obsHeader.obsType3[Currsysid].typeNum;
		int		len=3+16*typeNum-line.GetLength();
		for(int i=0;i<len;i++)
		{
			line+=_T(" ");
		}
		for(int k=0;k<typeNum;k++)
		{
		
			CString typeindex=obsHeader.obsType3[Currsysid].allType[k].obsCode;
			int     freqIndex=obsHeader.obsType3[Currsysid].allType[k].freCode;
			double	data=_wtof(line.Mid(3+k*16,14));
			int		 LLI=(int)_wtoi(line.Mid(3+k*16+15,1));
			int		SNR=(int)_wtoi(line.Mid(3+k*16+16,1));
			insertObsData3( Currsysid+1,data,freqIndex, LLI,SNR,typeindex,obsdatarecord);
		}
	
}
bool ProcessRinex::ReadObsData3FirstLine(CStdioFile& Ofile,CString& line,ObsEpochData& obsepochdata,int eventflag)
{
	int flag=_wtoi(line.Mid(31,1));
	int skip=_wtoi(line.Mid(32,3));
	eventflag=flag;
	if(flag!=0 && flag!=1)
	{
		for(int i=0;i<skip;i++)
		{
			Ofile.ReadString(line);
		}
		return false;
	}
		YMDHMS obsTime;
	//Gtime		protime;
	
	obsTime.year=_wtoi(line.Mid(2,4));
	obsTime.mon=_wtoi(line.Mid(7,2));
	obsTime.day=_wtoi(line.Mid(10,2));
	obsTime.hour=_wtoi(line.Mid(13,2));
	obsTime.min=_wtoi(line.Mid(16,2));
	obsTime.sec=_wtof(line.Mid(19,11));
	WeekSec(obsepochdata.week,obsepochdata.sec,obsTime,1);
	obsepochdata.flag=_wtoi(line.Mid(31,1));
	obsepochdata.sateNum=_wtoi(line.Mid(32,3));
	return true;
}
void ProcessRinex::insertObsData3(int Currsysid,double data ,int freqIndex,int LLI,int SNR,CString typeindex,ObsDataRecord& obsdatarecord)
{
		if( freqIndex==7 && Currsysid!=3)
		{
			freqIndex-=5;
		}
		else if ( freqIndex==6 && Currsysid!=3)  //E5a=L5 and E6=L6 can't exist at the same time
		{									
			freqIndex-=3;
		}
		else if(freqIndex==5 && Currsysid!=3)
		{
			freqIndex-=2;
		}
		else if (Currsysid==3)
		{//for Galileo    E1=L1; E5a=L5;	E5b=L7; E5=L8;	E6=L6  order E1 L7 L6, ext L5 L8
			if (freqIndex==6)
			{
				freqIndex-=3;
			}
			else if (freqIndex==7)
			{
				freqIndex-=5;
			}
			/*else if (freqIndex==5)
			{
				freqIndex-=
				extension[]
			}
			else if (freqIndex==8)
			{
			}*/
		}
	if(typeindex=='C' || typeindex=='P')
	{
		
		if (data>0.0 && obsdatarecord.vadFlgCod[freqIndex-1]==0)
		{
			obsdatarecord.PsRange[freqIndex-1]=data;//for 2.x , 3.x is different
			obsdatarecord.numVadCod++;
			obsdatarecord.vadFlgCod[freqIndex-1]=1;
		}
		//if (Currsysid==3&&freqIndex==5)
		//{
		//}
		//else if (Currsysid==3&&freqIndex==8)
		//{
		//}
	}
	else if(typeindex=='L' && obsdatarecord.vadFlgPhs[freqIndex-1]==0)//
	{
		if (data!=0.0 && obsdatarecord.vadFlgPhs[freqIndex-1]==0)
		{
			obsdatarecord.Phase[freqIndex-1]=data;//for 2.x , 3.x is different
			obsdatarecord.numVadPhs++;
			obsdatarecord.vadFlgPhs[freqIndex-1]=1;
		}
		
		obsdatarecord.LLI[freqIndex-1]=LLI;
		obsdatarecord.SNR[freqIndex-1]=SNR;
	}
	 else if(typeindex=='D')
	{
		obsdatarecord.Doppler[freqIndex-1]=data;//for 2.x , 3.x is different
	}
}

//------------------------------------Eph       File--------------------------------------------------


void ProcessRinex::ReadBrdEphHeader(CStdioFile& Nfile,CString& line,BroadEphHeader& brdephheader )
{
	while(Nfile.ReadString(line))
	{
		if(line.Find(_T("END OF HEADER"))!=-1) break;
		if(line.Find(_T("RINEX VERSION / TYPE"))!=-1 )
		{
			brdephheader.version=_wtof(line.Mid(0,9));
			if (line.Find(_T("GPS"))!=-1 || line.Mid(40,1)=="G")
			{
				brdephheader.sysid=1;
			}
			if (line.Find(_T("GLONASS"))!=-1|| line.Mid(40,1)=="R")
			{
				brdephheader.sysid=2;
			}
			if (line.Find(_T("GALILEO"))!=-1 || line.Mid(40,1)=="E")
			{
				brdephheader.sysid=3;
			}
			if (line.Find(_T("BDS"))!=-1 || line.Find(_T("CMP"))!=-1 || line.Mid(40,1)=="C")
			{
				brdephheader.sysid=5;
			}
			if (line.Mid(40,1)=="M")
			{
				brdephheader.sysid=6;
			}

			
			if(brdephheader.version>2.99)
			{
				brdephheader.sysid=Sysid(line.Mid(40,1));
			}

		}
		else if (line.Find(_T("IONOSPHERIC CORR"))!=-1)
		{
			if (line.Find(_T("GPSA"))!=-1)
			{
				brdephheader.ionoCorrGPS.alpha1=_wtof(line.Mid(5,12));
				brdephheader.ionoCorrGPS.alpha2=_wtof(line.Mid(17,12));
				brdephheader.ionoCorrGPS.alpha3=_wtof(line.Mid(29,12));
				brdephheader.ionoCorrGPS.alpha4=_wtof(line.Mid(41,12));
				brdephheader.ionoCorrGPS.validA=1;
			}
			if (line.Find(_T("GPSB"))!=-1 )
			{
				brdephheader.ionoCorrGPS.beta1=_wtof(line.Mid(5,12));
				brdephheader.ionoCorrGPS.beta2=_wtof(line.Mid(17,12));
				brdephheader.ionoCorrGPS.beta3=_wtof(line.Mid(29,12));
				brdephheader.ionoCorrGPS.beta4=_wtof(line.Mid(41,12));
				brdephheader.ionoCorrGPS.validB=1;
			}
			if (line.Find(_T("COMA"))!=-1 || line.Find(_T("BDSA"))!=-1)
			{
				brdephheader.ionoCorrBDS.alpha1=_wtof(line.Mid(5,12));
				brdephheader.ionoCorrBDS.alpha2=_wtof(line.Mid(17,12));
				brdephheader.ionoCorrBDS.alpha3=_wtof(line.Mid(29,12));
				brdephheader.ionoCorrBDS.alpha4=_wtof(line.Mid(41,12));
				brdephheader.ionoCorrBDS.validA=1;
			}
			if (line.Find(_T("COMB"))!=-1 || line.Find(_T("BDSB"))!=-1)
			{
				brdephheader.ionoCorrBDS.beta1=_wtof(line.Mid(5,12));
				brdephheader.ionoCorrBDS.beta2=_wtof(line.Mid(17,12));
				brdephheader.ionoCorrBDS.beta3=_wtof(line.Mid(29,12));
				brdephheader.ionoCorrBDS.beta4=_wtof(line.Mid(41,12));
				brdephheader.ionoCorrBDS.validB=1;
			}
			if (line.Find(_T("GAL"))!=-1)
			{
				brdephheader.ionoCorrGAL.alpha1=_wtof(line.Mid(5,12));
				brdephheader.ionoCorrGAL.alpha2=_wtof(line.Mid(17,12));
				brdephheader.ionoCorrGAL.alpha3=_wtof(line.Mid(29,12));
				brdephheader.ionoCorrGAL.alpha4=_wtof(line.Mid(41,12));
				brdephheader.ionoCorrGAL.validA=1;
			}
		}
		else if (line.Find(_T("TIME SYSTEM CORR")))
		{
			if(line.Find(_T("GAUT"))!=-1)
			{
				brdephheader.GAUT.valid=1;
				brdephheader.GAUT.a0=_wtof(line.Mid(5,16));
				brdephheader.GAUT.a1=_wtof(line.Mid(22,16));
				brdephheader.GAUT.refTime=_wtof(line.Mid(38,7));
				brdephheader.GAUT.refWeek=_wtoi(line.Mid(45,5));
			}
			if(line.Find(_T("GPUT"))!=-1)
			{
				brdephheader.GPUT.valid=1;
				brdephheader.GPUT.a0=_wtof(line.Mid(5,16));
				brdephheader.GPUT.a1=_wtof(line.Mid(22,16));
				brdephheader.GPUT.refTime=_wtof(line.Mid(38,7));
				brdephheader.GPUT.refWeek=_wtoi(line.Mid(45,5));
			}
			if(line.Find(_T("SBUT"))!=-1)
			{
				brdephheader.SBUT.valid=1;
				brdephheader.SBUT.a0=_wtof(line.Mid(5,16));
				brdephheader.SBUT.a1=_wtof(line.Mid(22,16));
				brdephheader.SBUT.refTime=_wtof(line.Mid(38,7));
				brdephheader.SBUT.refWeek=_wtoi(line.Mid(45,5));
			}
			if(line.Find(_T("GLUT"))!=-1)
			{
				brdephheader.GLUT.valid=1;
				brdephheader.GLUT.a0=_wtof(line.Mid(5,16));
				brdephheader.GLUT.a1=_wtof(line.Mid(22,16));
				brdephheader.GLUT.refTime=_wtof(line.Mid(38,7));
				brdephheader.GLUT.refWeek=_wtoi(line.Mid(45,5));
			}
			if(line.Find(_T("GPGA"))!=-1)
			{
				brdephheader.GPGA.valid=1;
				brdephheader.GPGA.a0=_wtof(line.Mid(5,16));
				brdephheader.GPGA.a1=_wtof(line.Mid(22,16));
				brdephheader.GPGA.refTime=_wtof(line.Mid(38,7));
				brdephheader.GPGA.refWeek=_wtoi(line.Mid(45,5));
			}
			if(line.Find(_T("GLGP"))!=-1)
			{
				brdephheader.GLGP.valid=1;
				brdephheader.GLGP.a0=_wtof(line.Mid(5,16));
				brdephheader.GLGP.a1=_wtof(line.Mid(22,16));
				brdephheader.GLGP.refTime=_wtof(line.Mid(38,7));
				brdephheader.GLGP.refWeek=_wtoi(line.Mid(45,5));
			}
		}
		else if(line.Find(_T("LEAP SECONDS"))!=-1)
		{
			brdephheader.leapSecond=_wtof(line.Mid(0,6));
		}
	}
}
//single satellite
void ProcessRinex::ReadBrdEphData(CStdioFile& Nfile,CString& line,int sysid,double version,BroadEphData& brdephdata,BroadEphDataGlo& GloEph,int& flag)
{
	int n=(version>2.99)?1:0;
	
	YMDHMS  Toc;
	int prn=0;	
	if (version>2.99)
	{
		prn					=Prn( line.Mid(0,1), _wtoi(line.Mid(1,2)));
	}
	else if(version<2.99)
	{	//no mixed in 2.x
		CString sys;						
		sys					=Sys(sysid);
		prn					=Prn(sys, _wtoi(line.Mid(0,2)));
	}
	if(version>2.99 &&((prn>100&&prn<150) || prn<50 ||prn>200))
	{
		flag=0;
		Toc.year							=_wtoi(line.Mid(4,4));
		Toc.mon							=_wtoi(line.Mid(9,2));
		Toc.day								=_wtoi(line.Mid(12,2));
		Toc.hour							=_wtoi(line.Mid(15,2));
		Toc.min							=_wtoi(line.Mid(18,2));
		Toc.sec								=_wtof(line.Mid(21,2));
		int week=0;
		double sec=0.0;
		Gtime		temp;
		WeekSec(week,sec,Toc,Sysid(prn));
		brdephdata.toc					=sec;
		brdephdata.prn					=prn;
	}
	else if(version<2.99 &&((prn>100&&prn<150) || prn<50 ||prn>200))
	{
		flag=0;
		int y									=_wtoi(line.Mid(3,2));
		Toc.year							=y>79?y+1900:y+2000;
		Toc.mon							=_wtoi(line.Mid(6,2));
		Toc.day								=_wtoi(line.Mid(9,2));
		Toc.hour							=_wtoi(line.Mid(12,2));
		Toc.min							=_wtoi(line.Mid(15,2));
		Toc.sec								=_wtof(line.Mid(17,5));
		int week=0;
		double sec=0.0;
		Gtime		temp;
		WeekSec(week,sec,Toc,Sysid(prn));
		brdephdata.toc					=sec;
		brdephdata.prn					=prn;
	}


	if (prn>50 && prn<100)
	{
		//for Glonass
		flag=1;
		GloEph.prn					=prn;	
		if(version>2.99)
		{
			Toc.year							=_wtoi(line.Mid(4,4));
			Toc.mon							=_wtoi(line.Mid(9,2));
			Toc.day								=_wtoi(line.Mid(12,2));
			Toc.hour							=_wtoi(line.Mid(15,2));
			Toc.min							=_wtoi(line.Mid(18,2));
			Toc.sec								=_wtof(line.Mid(21,2));
			
		}
		else if(version<2.99)
		{
			int y									=_wtoi(line.Mid(3,2));
			Toc.year							=y>79?y+1900:y+2000;
			Toc.mon							=_wtoi(line.Mid(6,2));
			Toc.day								=_wtoi(line.Mid(9,2));
			Toc.hour							=_wtoi(line.Mid(12,2));
			Toc.min							=_wtoi(line.Mid(15,2));
			Toc.sec								=_wtof(line.Mid(17,5));
		
		}

		GloEph.NegativeTauN		=_wtof(line.Mid(22+n,19));
		GloEph.PositiveGammaN	=_wtof(line.Mid(41+n,19));
		double tof				=_wtof(line.Mid(60+n,19));

		int week=0;
		double sec=0.0;
		WeekSec(week,sec,Toc,Sysid(prn));// to GPSw andGPSS
		double sectoc=floor((sec+450.0)/900.0)*900.0;
		int dow=(int)floor(sec/86400.0);
		double tod=version<2.99?tof:fmod(tof,86400.0);
		tof=tod+dow*86400.0;

		GloEph.idoe=(int)(fmod(sec+10800.0,86400.0)/900.0+0.5);
		sec=UTC2GPST(Toc.year,Toc.mon,sec);//trans the utc to gps, treat the obs time is gpst by default
		if (sec<0.0)
		{
			week--;
			sec+=86400.0*7.0;
		}
		if (sec>86400.0*7.0)
		{
			week++;
			sec-=86400.0*7.0;
		}
		GloEph.toc					=sec;
		GloEph.toe					=sec;
		GloEph.week		=week;
		Nfile.ReadString(line);
		GloEph.Pos[0]				=_wtof(line.Mid(3+n,19))*1E3;
		GloEph.Vel[0]					=_wtof(line.Mid(22+n,19))*1E3;
		GloEph.Acc[0]		=_wtof(line.Mid(41+n,19))*1E3;
		GloEph.health					=_wtoi(line.Mid(60+n,19));

		Nfile.ReadString(line);
		GloEph.Pos[1]				=_wtof(line.Mid(3+n,19))*1E3;
		GloEph.Vel[1]					=_wtof(line.Mid(22+n,19))*1E3;
		GloEph.Acc[1]			=_wtof(line.Mid(41+n,19))*1E3;
		GloEph.freqNum					=_wtoi(line.Mid(60+n,19));

		Nfile.ReadString(line);
		GloEph.Pos[2]				=_wtof(line.Mid(3+n,19))*1E3;
		GloEph.Vel[2]					=_wtof(line.Mid(22+n,19))*1E3;
		GloEph.Acc[2]			=_wtof(line.Mid(41+n,19))*1E3;
		GloEph.Age					=_wtoi(line.Mid(60+n,19));
	} 
	else if (prn>150&&prn<200) //for sbas
	{
		for(int i=0;i<3;i++) Nfile.ReadString(line);
		flag=2;
	}
	else
	{
		brdephdata.satClkBs			=_wtof(line.Mid(22+n,19));
		brdephdata.satClkDrt		=_wtof(line.Mid(41+n,19));
		brdephdata.satClkDrtRe	=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//1
		brdephdata.iode				=_wtof(line.Mid(3+n,19));
		brdephdata.crs					=_wtof(line.Mid(22+n,19));
		brdephdata.deltan			=_wtof(line.Mid(41+n,19));
		brdephdata.m0					=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//2
		brdephdata.cuc			=_wtof(line.Mid(3+n,19));
		brdephdata.e			=_wtof(line.Mid(22+n,19));
		brdephdata.cus			=_wtof(line.Mid(41+n,19));
		brdephdata.sqA		=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//3
		brdephdata.toe					=_wtof(line.Mid(3+n,19));
		brdephdata.cic					=_wtof(line.Mid(22+n,19));
		brdephdata.Omega			=_wtof(line.Mid(41+n,19));
		brdephdata.cis					=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//4
		brdephdata.i0					=_wtof(line.Mid(3+n,19));
		brdephdata.crc					=_wtof(line.Mid(22+n,19));
		brdephdata.omega			=_wtof(line.Mid(41+n,19));
		brdephdata.OmegaDot		=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//5
		brdephdata.idot					=_wtof(line.Mid(3+n,19));
		brdephdata.codeL2				=_wtof(line.Mid(22+n,19));
		brdephdata.gpsw					=_wtof(line.Mid(41+n,19));
		brdephdata.flagL2					=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//6
		brdephdata.satAccuracy		=_wtof(line.Mid(3+n,19));
		brdephdata.satHealth			=_wtof(line.Mid(22+n,19));
		brdephdata.tgd						=_wtof(line.Mid(41+n,19));
		brdephdata.iodClk				=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//7
		int len=80-line.GetLength();
		for(int i=0;i<len;i++) 	line+=_T(" ");
		brdephdata.transTime			=_wtof(line.Mid(3+n,19));
		brdephdata.fitInterval			=_wtof(line.Mid(22+n,19));
		brdephdata.spare1				=_wtof(line.Mid(41+n,19));
		brdephdata.spare2				=_wtof(line.Mid(60+n,19));
	}
	
}


void ProcessRinex::ReadBrdEphData(CStdioFile& Nfile,CString& line,int sysid,double version,BroadEphData& brdephdata,BroadEphDataGlo& GloEph,BroadEphDataSBAS& sbseph,int& flag)
{
	int n=(version>2.99)?1:0;

	YMDHMS  Toc;
	int prn=0;	
	if (version>2.99)
	{
		prn									=Prn( line.Mid(0,1), _wtoi(line.Mid(1,2)));
		Toc.year							=_wtoi(line.Mid(4,4));
		Toc.mon							=_wtoi(line.Mid(9,2));
		Toc.day								=_wtoi(line.Mid(12,2));
		Toc.hour							=_wtoi(line.Mid(15,2));
		Toc.min							=_wtoi(line.Mid(18,2));
		Toc.sec								=_wtof(line.Mid(21,2));
	}
	else if(version<2.99)
	{	//no mixed in 2.x
		CString sys;						
		sys					=Sys(sysid);
		prn					=Prn(sys, _wtoi(line.Mid(0,2)));
		int y									=_wtoi(line.Mid(3,2));
		Toc.year							=y>79?y+1900:y+2000;
		Toc.mon							=_wtoi(line.Mid(6,2));
		Toc.day								=_wtoi(line.Mid(9,2));
		Toc.hour							=_wtoi(line.Mid(12,2));
		Toc.min							=_wtoi(line.Mid(15,2));
		Toc.sec								=_wtof(line.Mid(17,5));
	}

	int week=0;
	double sec=0.0;
	Gtime		temp;
	WeekSec(week,sec,Toc,Sysid(prn));//suppose the reference time of glo is same as gps

	if(  (prn>100 && prn<150) || prn<50 || prn>200)
	{
		flag=0;
		brdephdata.toc					=sec;
		brdephdata.prn					=prn;
	}
		//for Glonass
	if (prn>50 && prn<100)
	{
		flag=1;
		GloEph.prn					=prn;	
		GloEph.NegativeTauN		=_wtof(line.Mid(22+n,19));
		GloEph.PositiveGammaN	=_wtof(line.Mid(41+n,19));
		double tof				=_wtof(line.Mid(60+n,19));

		//WeekSec(week,sec,Toc,Sysid(prn));// to GPSw andGPSS
		double sectoc=floor((sec+450.0)/900.0)*900.0;
		int dow=(int)floor(sec/86400.0);
		double tod=version<2.99?tof:fmod(tof,86400.0);
		tof=tod+dow*86400.0;

		GloEph.idoe=(int)(fmod(sec+10800.0,86400.0)/900.0+0.5);
		sec=UTC2GPST(Toc.year,Toc.mon,sec);
		if (sec<0.0)
		{
			week--;
			sec+=86400.0*7.0;
		}
		if (sec>86400.0*7.0)
		{
			week++;
			sec-=86400.0*7.0;
		}
		GloEph.toc					=sec;
		GloEph.toe					=sec;
		GloEph.week		=week;
		Nfile.ReadString(line);
		GloEph.Pos[0]				=_wtof(line.Mid(3+n,19))*1E3;
		GloEph.Vel[0]					=_wtof(line.Mid(22+n,19))*1E3;
		GloEph.Acc[0]		=_wtof(line.Mid(41+n,19))*1E3;
		GloEph.health					=_wtoi(line.Mid(60+n,19));

		Nfile.ReadString(line);
		GloEph.Pos[1]				=_wtof(line.Mid(3+n,19))*1E3;
		GloEph.Vel[1]					=_wtof(line.Mid(22+n,19))*1E3;
		GloEph.Acc[1]			=_wtof(line.Mid(41+n,19))*1E3;
		GloEph.freqNum					=_wtoi(line.Mid(60+n,19));

		Nfile.ReadString(line);
		GloEph.Pos[2]				=_wtof(line.Mid(3+n,19))*1E3;
		GloEph.Vel[2]					=_wtof(line.Mid(22+n,19))*1E3;
		GloEph.Acc[2]			=_wtof(line.Mid(41+n,19))*1E3;
		GloEph.Age					=_wtoi(line.Mid(60+n,19));
	} 
	else if (prn>150&&prn<200)// sbas
	{
		flag=2;
		sbseph.a0		=_wtof(line.Mid(22+n,19));
		sbseph.a1		=_wtof(line.Mid(41+n,19));
		sbseph.tof		=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);
		sbseph.pos[0]	=_wtof(line.Mid(3+n,19));
		sbseph.vel[0]		=_wtof(line.Mid(22+n,19));
		sbseph.acc[0]	=_wtof(line.Mid(41+n,19));
		sbseph.svh		=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);
		sbseph.pos[1]	=_wtof(line.Mid(3+n,19));
		sbseph.vel[1]		=_wtof(line.Mid(22+n,19));
		sbseph.acc[1]	=_wtof(line.Mid(41+n,19));
		sbseph.sva		=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);
		sbseph.pos[2]	=_wtof(line.Mid(3+n,19));
		sbseph.vel[2]		=_wtof(line.Mid(22+n,19));
		sbseph.acc[2]	=_wtof(line.Mid(41+n,19));
		sbseph.iodn		=_wtof(line.Mid(60+n,19));
	}
	else
	{
		brdephdata.satClkBs			=_wtof(line.Mid(22+n,19));
		brdephdata.satClkDrt		=_wtof(line.Mid(41+n,19));
		brdephdata.satClkDrtRe	=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//1
		brdephdata.iode				=_wtof(line.Mid(3+n,19));
		brdephdata.crs					=_wtof(line.Mid(22+n,19));
		brdephdata.deltan			=_wtof(line.Mid(41+n,19));
		brdephdata.m0					=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//2
		brdephdata.cuc			=_wtof(line.Mid(3+n,19));
		brdephdata.e			=_wtof(line.Mid(22+n,19));
		brdephdata.cus			=_wtof(line.Mid(41+n,19));
		brdephdata.sqA		=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//3
		brdephdata.toe					=_wtof(line.Mid(3+n,19));
		brdephdata.cic					=_wtof(line.Mid(22+n,19));
		brdephdata.Omega			=_wtof(line.Mid(41+n,19));
		brdephdata.cis					=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//4
		brdephdata.i0					=_wtof(line.Mid(3+n,19));
		brdephdata.crc					=_wtof(line.Mid(22+n,19));
		brdephdata.omega			=_wtof(line.Mid(41+n,19));
		brdephdata.OmegaDot		=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//5
		brdephdata.idot					=_wtof(line.Mid(3+n,19));
		brdephdata.codeL2				=_wtof(line.Mid(22+n,19));
		brdephdata.gpsw					=_wtof(line.Mid(41+n,19));
		brdephdata.flagL2					=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//6
		brdephdata.satAccuracy		=_wtof(line.Mid(3+n,19));
		brdephdata.satHealth			=_wtof(line.Mid(22+n,19));
		brdephdata.tgd						=_wtof(line.Mid(41+n,19));
		brdephdata.iodClk				=_wtof(line.Mid(60+n,19));

		Nfile.ReadString(line);//7
		int len=80-line.GetLength();
		for(int i=0;i<len;i++)
		{
			line+=_T(" ");
		}
		brdephdata.transTime			=_wtof(line.Mid(3+n,19));
		brdephdata.fitInterval			=_wtof(line.Mid(22+n,19));
		brdephdata.spare1				=_wtof(line.Mid(41+n,19));
		brdephdata.spare2				=_wtof(line.Mid(60+n,19));
	}

}

void ProcessRinex::ReadBroadEph(CString NFileName,BroadEphHeader& brdephheader,BroadEphData* broadephdata,BroadEphDataGlo* gloeph,int& nSate,int& nSateGlo) 
{
	CStdioFile Nfile;
	Nfile.Open(NFileName,CFile::modeRead);
	CString line;
	ReadBrdEphHeader(Nfile,line,brdephheader);
	int		sysid			=brdephheader.sysid;
	double	version		=brdephheader.version;
	int flag=0;   int nSateSbas=0;
	while(Nfile.ReadString(line))
	{
		ReadBrdEphData(Nfile,line,sysid,version,broadephdata[nSate],gloeph[nSateGlo],flag);
		if (flag==0)	nSate++; else if(flag==1) 	nSateGlo++; else if(flag==2)		nSateSbas++;
	}
		if (flag==0)	nSate++; else if(flag==1) 	nSateGlo++; else if(flag==2)		nSateSbas++;
	Nfile.Close();
}




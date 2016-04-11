#include "StdAfx.h"
#include "atlstr.h"
#include<fstream>
#include<iomanip>
#include<stdlib.h>
using namespace std;
#include "RTCM_NMEA.h"




//==============class as struct===================
gpgga::gpgga()
{
	GP_ID="F";
	utc=0.0;
	lat=0.0;
	ns_ind="F";
	lon=0.0;
	ew_ind="F";
	pos_fix=0;
	num_sate=0;
	hdop=0.0;
	alt=0.0;
	alt_unit="F";
	hei=0.0;
	sep_unit="F";
	time_frac=0.0;
	diff_id="F";

}

gpgsa::gpgsa()
{
	GP_ID="F";
	mode="F";
	pos_dim=0;
	for (int i=0;i<MAXNUMGSA;i++)
	{
		PRN[i]=0;
	}
	pdop=0.0;
	hdop=0.0;
	vdop=0.0;
	validnum=0;

}

sate_info::sate_info()
{
	PRN=0;
	ele=0.0;
	azi=0.0;
	snr=0.0;
}

gpgsv::gpgsv()
{
	GP_ID="F";
	numsat=0;
}

gprmc::gprmc()
{
	GP_ID="F";
	utc=0.0;
	sta="F";
	lat=0.0;
	ns_ind="F";
	lon=0.0;
	ew_ind="F";
	vel=0.0;
	cou=0.0;
	utc_date=0;
	mag_var=0.0;
	mv_ind="F";

}

//===the other GP are not used====

SingleCorrInfo::SingleCorrInfo()
{
	sysid=0;
	prn=0;
	corrCom=0.0;
	corrIono=0.0;
	corrTrop=0.0;
	corrDot=0.0;
}
EpochCorrInfo::EpochCorrInfo()
{
	GPSW=0;
	GPSS=0;
	NumSat=0;
}



//-----------------------Class NMEA_process--------------------------
bool NMEA_process::read_data(CString filename)
{
	CStdioFile file;
	file.Open(filename,CFile::modeRead);
	
	CString line;
	int isok=0;
	while( file.ReadString(line)  )
	{
		check_sum(line,isok);
		if (isok!=100)
		{
			//continue;
		}
		CString s=	line.Mid(1,5);
			if (s=="GPGGA")//1
				{
					Get_gpgga(line,gga_pub);
				}
			else if (s=="GPGSA" )//2
				{
					Get_gpgsa(line,gsa_pub);
				}
					else if (s=="GPGSV")	//3
				{
					CString all;	char deli=',';
					int len;int star;int n1;
					int numline	=_wtoi( line.Mid(7,1) );
					int lineNo;  //No is the place;num is the totla number
					for (int i=0; i<numline; i++)
					{
						star	=line.Find('*');			n1	=line.Find(',',11);
						len	=	line.GetLength();	lineNo	=_wtoi( line.Mid(9,1) );
						if ( isok==100)//check, need not to check before myclass.Get_gpgsv
						{
							
							if ( lineNo==1  )
							{
								//if (numline != 1)
								//{
									all=line.Mid(0,star);
									file.ReadString(line);
									check_sum(line,isok);				//this cycle fails if line isn't ok
									if (isok!=100)
									{
										break;
									}
								//}
							}
							else if ( lineNo==2 )
							{
								//if (lineNo !=numline)
								//{
									all=all+line.Mid(n1,star-n1);
									file.ReadString(line);
									check_sum(line,isok);//this cycle fails if line isn't ok
									if (isok!=100)
									{
										break;
									}
								//}
								//else
								//{
									//all=all+line.Mid(n1,len-n1);
								//}
							}
							else if (lineNo==3)
							{
								all=all+line.Mid(n1,len-n1);
							}
						}
						
					}
					if (isok==100)
						{
							Get_gpgsv(all,gsv_pub);
						}
		
				}
			else if (s=="GPRMC")	//4
				{
					Get_gprmc(line,rmc_pub);
				}
			else if (s=="GPRMC3")	//5
				{
					Get_gprmc3(line,rmc3_pub);
				}
			else if (s=="GPVTG")	//6
				{
					Get_gpvtg(line,vtg_pub);
				}
			else if (s=="GPVTG3")	//7
				{
					Get_gpvtg3(line,vtg3_pub);
				}
			else if (s=="GPGLL")	//8
				{
					Get_gpgll(line,gll_pub);
				}
			else if (s=="GPZDA")	//9
				{
					Get_gpzda(line,zda_pub);
				}
			else if (s=="GPGST")	//10
				{
					Get_gpgst(line,gst_pub);
				}
		
	}

	return true;
}

bool NMEA_process::check_sum(CString line,int &flag)
{
	int len=line.GetLength();
	CString  line_check=line.Right(2);
	unsigned char checksum=0;
	unsigned char* temUCP;
	unsigned char temUC;
	CString temp;
		
	for (int i=0;i<len-4;i++)  //4:$*XX
	{
		temp=line.Mid(i+1,1);
		temUCP=(unsigned char*)(LPCTSTR)temp;
		temUC  =*temUCP;
		checksum=checksum^temUC;
	}
	CString me_check;
	me_check.Format(_T("%02X"),checksum); 
	if (me_check==line_check)
	{
		
		flag=100;
		return true;
	}
	else
	{
		return false;
	}
}

bool NMEA_process::Get_gpgga(CString line,gpgga &gga)
{
	CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(deli);
	int star=line.Find('*');
	temp=line.Mid(1,5);
	gga.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}

		switch (i)
		{
			case 1 :
				{
					if (strlen!=0)
					{
						gga.utc=_wtof(temp);
					}
					else
					{
						gga.utc=-1.0;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						gga.lat=_wtof(temp);
					}
					else
					{
						gga.lat=-1;
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						gga.ns_ind=temp;
					}
					else
					{
						gga.ns_ind=' ';
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						gga.lon=_wtof(temp);
					}
					else
					{
						gga.lon=-1;
					}
				}
			case 5:
				{
					if (strlen!=0)
					{
						gga.ew_ind=temp;
					}
					else
					{
						gga.ew_ind=' ';
					}
				}
			case 6:
				{
					if (strlen!=0)
					{
						gga.pos_fix=_wtoi(temp);
					}
					else
					{
						gga.pos_fix=-1;
					}
				}
			case 7:
				{
					if (strlen!=0)
					{
						gga.num_sate=_wtoi(temp);
					}
					else
					{
						gga.num_sate=-1;
					}
				}
			case 8:
				{
					if (strlen!=0)
					{
						gga.hdop=_wtof(temp);
					}
					else
					{
						gga.hdop=-1;
					}
				}
			case 9:
				{
					if (strlen!=0)
					{
						gga.alt=_wtof(temp);
					}
					else
					{
						gga.alt=-1;
					}
				}
			case 10:
				{
					if (strlen!=0)
					{
						gga.alt_unit=temp;
					
					}
					else
					{
						gga.alt_unit=' ';
					}
				}
			case 11:
				{
					if (strlen!=0)
					{
						gga.hei=_wtof(temp);
					}
					else
					{
						gga.hei=NULL;
					}
				}
			case 12:
				{
					if (strlen!=0)
					{
						gga.sep_unit=temp;
					}
					else
					{
						gga.sep_unit=' ';
					}
				}
			case 13:
				{
					if (strlen!=0)
					{
						gga.time_frac=_wtof(temp);
					}
					else
					{
						gga.time_frac=-1;
					}
				}
			case 14:
				{
					if (strlen!=0)
					{
						gga.diff_id=temp;
					}
					else
					{
						gga.diff_id=' ';
					}
				}
		}
		if (i==14) 
		{
			break;
		}
	}
	return true;
}

bool NMEA_process::Get_gpgsa(CString line,gpgsa &gsa)

{
	CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(',');
	int star=line.Find('*');
	int flag;
	temp=line.Mid(1,5);
	gsa.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	int j=0;int validNum=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}
		if (i<3 )
		{
			flag=i;
		}
		else if (i>14)
		{
			flag=i-11;
		}
		else
		{
			flag=3;
		}

		if (flag==1)		
		{
					if (strlen!=0)
					{
						gsa.mode=temp;
					}
					else
					{
						gsa.mode=' ';
					}
				}
		else if (flag==2)	
				{
					if (strlen!=0)
					{
						gsa.pos_dim=_wtoi(temp);
					}
					else
					{
						gsa.pos_dim=-1;
					}
				}
		else if (flag==3)	  //to be modified
				{
						if (strlen!=0)
						{
							gsa.PRN[j]=_wtoi(temp);
							validNum++;
						}
						else 
						{
							gsa.PRN[j]=-1;
						}
						j++;
				}
		else if (flag==4)	
				{
					if (strlen!=0)
					{
						gsa.pdop=_wtof(temp);
					}
					else
					{
						gsa.pdop=-1;
					}
				}
	else if (flag==5)	
				{
					if (strlen!=0)
					{
						gsa.hdop=_wtof(temp);
					}
					else
					{
						gsa.hdop=-1;
					}
				}
	else if (flag==6)	
				{
					if (strlen!=0)
					{
						gsa.vdop=_wtof(temp);
					}
					else
					{
						gsa.vdop=-1;
					}
					gsa.validnum=validNum;
					break;
				}
		

	}
	return true;
}

bool NMEA_process::Get_gpgsv(CString line,gpgsv &gsv)
{

	CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(',');
	int star=line.Find('*');
	temp=line.Mid(1,5);
	gsv.GP_ID=temp;
	gsv.numsat=_wtoi(line.Mid(11,2));
	int n1=0;	int n2=0;	int strlen=0;
	int j=0;int k=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n1+1,strlen);
		}
		if (i%4==0)
		{
			j++;
			k=i-4*j;
		}
		else if(i<4)
		{
			continue;
		}
		k++;
		switch (k)
		{
			case 1 :
				{
					if (strlen!=0)
					{
						gsv.sat[j-1].PRN=_wtoi(temp);
					}
					else
					{
						gsv.sat[j-1].PRN=-1;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						gsv.sat[j-1].ele=_wtof(temp);
					}
					else
					{
						gsv.sat[j-1].ele=-1;
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						gsv.sat[j-1].azi=_wtof(temp);
					}
					else
					{
						gsv.sat[j-1].azi=-1;
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						gsv.sat[j-1].snr=_wtof(temp);
					}
					else
					{
						gsv.sat[j-1].snr=-1;
					}
				}
			
		}
		if (j==gsv.numsat&&k==4) 
		{
			break;
		}
	}
	return true;
}

bool NMEA_process::Get_gprmc(CString line, gprmc &rmc)
{
		CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(',');
	int star=line.Find('*');
	temp=line.Mid(1,5);
	rmc.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}

		switch (i)
		{
			case 1 :
				{
					if (strlen!=0)
					{
						rmc.utc=_wtof(temp);
					}
					else
					{
						rmc.utc=-1;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						rmc.sta=temp;
					}
					else
					{
						rmc.sta=' ';
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						rmc.lat=_wtof(temp);
					}
					else
					{
						rmc.lat=-1;
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						rmc.ns_ind=(temp);
					}
					else
					{
						rmc.ns_ind=' ';
					}
				}
			case 5:
				{
					if (strlen!=0)
					{
						rmc.lon=_wtof(temp);
					}
					else
					{
						rmc.lon=-1;
					}
				}
			case 6:
				{
					if (strlen!=0)
					{
						rmc.ew_ind=temp;
					}
					else
					{
						rmc.ew_ind=' ';
					}
				}
			case 7:
				{
					if (strlen!=0)
					{
						rmc.vel=_wtoi(temp);
					}
					else
					{
						rmc.vel=-1;
					}
				}
			case 8:
				{
					if (strlen!=0)
					{
						rmc.cou=_wtof(temp);
					}
					else
					{
						rmc.cou=-1;
					}
				}
			case 9:
				{
					if (strlen!=0)
					{
						rmc.utc_date=_wtoi(temp);
					}
					else
					{
						rmc.utc_date=-1;
					}
				}
			case 10:
				{
					if (strlen!=0)
					{
						rmc.mag_var=_wtof(temp);
					}
					else
					{
						rmc.mag_var=-1;
					}
				}
			case 11:
				{
					if (strlen!=0)
					{
						rmc.mv_ind=temp;
					
					}
					else
					{
						rmc.mv_ind=' ';
					}
				}
			
		}
		if (i==11) 
		{
			break;
		}
	}
	return true;
}

bool NMEA_process::Get_gprmc3(CString line,gprmc3 &rmc3)
{
	
		CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(',');
	int star=line.Find('*');
	temp=line.Mid(1,5);
	rmc3.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}

		switch (i)
		{
			case 1 :
				{
					if (strlen!=0)
					{
						rmc3.utc=_wtof(temp);
					}
					else
					{
						rmc3.utc=-1;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						rmc3.sta=temp;
					}
					else
					{
						rmc3.sta=' ';
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						rmc3.lat=_wtof(temp);
					}
					else
					{
						rmc3.lat=-1;
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						rmc3.ns_ind=(temp);
					}
					else
					{
						rmc3.ns_ind=' ';
					}
				}
			case 5:
				{
					if (strlen!=0)
					{
						rmc3.lon=_wtof(temp);
					}
					else
					{
						rmc3.lon=-1;
					}
				}
			case 6:
				{
					if (strlen!=0)
					{
						rmc3.ew_ind=temp;
					}
					else
					{
						rmc3.ew_ind=' ';
					}
				}
			case 7:
				{
					if (strlen!=0)
					{
						rmc3.vel=_wtoi(temp);
					}
					else
					{
						rmc3.vel=-1;
					}
				}
			case 8:
				{
					if (strlen!=0)
					{
						rmc3.cou=_wtof(temp);
					}
					else
					{
						rmc3.cou=-1;
					}
				}
			case 9:
				{
					if (strlen!=0)
					{
						rmc3.utc_date=_wtoi(temp);
					}
					else
					{
						rmc3.utc_date=-1;
					}
				}
			case 10:
				{
					if (strlen!=0)
					{
						rmc3.mag_var=_wtof(temp);
					}
					else
					{
						rmc3.mag_var=-1;
					}
				}
			case 11:
				{
					if (strlen!=0)
					{
						rmc3.mv_ind=temp;
					
					}
					else
					{
						rmc3.mv_ind=' ';
					}
				}
			case 12:
				{
					if (strlen!=0)
					{
						rmc3.mode_ind=temp;
					
					}
					else
					{
						rmc3.mode_ind=' ';
					}
				}
			
		}
		if (i==12)
		{
			break;
		}
	}
	return true;
}

bool NMEA_process::Get_gpvtg(CString line,gpvtg &vtg)
{
			CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(',');
	int star=line.Find('*');
	temp=line.Mid(1,5);
	vtg.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}

		switch (i)
		{
			case 1 :
				{
					if (strlen!=0)
					{
						vtg.cou1=_wtof(temp);
					}
					else
					{
						vtg.cou1=-1;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						vtg.ref1=temp;
					}
					else
					{
						vtg.ref1=' ';
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						vtg.cou2=_wtof(temp);
					}
					else
					{
						vtg.cou2=-1;
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						vtg.ref2=(temp);
					}
					else
					{
						vtg.ref2=' ';
					}
				}
			case 5:
				{
					if (strlen!=0)
					{
						vtg.vel1=_wtof(temp);
					}
					else
					{
						vtg.vel1=-1;
					}
				}
			case 6:
				{
					if (strlen!=0)
					{
						vtg.unit1=temp;
					}
					else
					{
						vtg.unit1=' ';
					}
				}
			case 7:
				{
					if (strlen!=0)
					{
						vtg.vel2=_wtof(temp);
					}
					else
					{
						vtg.vel2=-1;
					}
				}
			case 8:
				{
					if (strlen!=0)
					{
						vtg.unit2=temp;
					}
					else
					{
						vtg.unit2=' ';
					}
				}
						
		}
		if (i==8) 
		{
			break;
		}
	}
	return true;
}

bool NMEA_process::Get_gpvtg3(CString line,gpvtg3 &vtg3)
{
				CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(',');
	int star=line.Find('*');
	temp=line.Mid(1,5);
	vtg3.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}

		switch (i)
		{
			case 1 :
				{
					if (strlen!=0)
					{
						vtg3.cou1=_wtof(temp);
					}
					else
					{
						vtg3.cou1=-1;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						vtg3.ref1=temp;
					}
					else
					{
						vtg3.ref1=' ';
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						vtg3.cou2=_wtof(temp);
					}
					else
					{
						vtg3.cou2=-1;
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						vtg3.ref2=temp;
					}
					else
					{
						vtg3.ref2=' ';
					}
				}
			case 5:
				{
					if (strlen!=0)
					{
						vtg3.vel1=_wtof(temp);
					}
					else
					{
						vtg3.vel1=-1;
					}
				}
			case 6:
				{
					if (strlen!=0)
					{
						vtg3.unit1=temp;
					}
					else
					{
						vtg3.unit1=' ';
					}
				}
			case 7:
				{
					if (strlen!=0)
					{
						vtg3.vel2=_wtof(temp);
					}
					else
					{
						vtg3.vel2=-1;
					}
				}
			case 8:
				{
					if (strlen!=0)
					{
						vtg3.unit2=temp;
					}
					else
					{
						vtg3.unit2=' ';
					}
				}
			case 9:
				{
					if (strlen!=0)
					{
						vtg3.unit2=temp;
					}
					else
					{
						vtg3.unit2=' ';
					}
				}
				case 10:
				{
					if (strlen!=0)
					{
						vtg3.mode_ind=temp;
					}
					else
					{
						vtg3.mode_ind=' ';
					}
				}
						
		}
		if (i==10) 
		{
			break;
		}
	}
	return true;
}

bool NMEA_process::Get_gpgll(CString line,gpgll &gll)
{
	CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(',');
	int star=line.Find('*');
	temp=line.Mid(1,5);
	gll.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}

		switch (i)
		{
			
			case 1:
				{
					if (strlen!=0)
					{
						gll.lat=_wtof(temp);
					}
					else
					{
						gll.lat=-1;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						gll.ns_ind=temp;
					}
					else
					{
						gll.ns_ind=' ';
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						gll.lon=_wtof(temp);
					}
					else
					{
						gll.lon=-1;
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						gll.ew_ind=temp;
					}
					else
					{
						gll.ew_ind=' ';
					}
				}
			case 5 :
				{
					if (strlen!=0)
					{
						gll.utc=_wtof(temp);
					}
					else
					{
						gll.utc=-1;
					}
				}
			case 6:
				{
					if (strlen!=0)
					{
						gll.sta=temp;
					}
					else
					{
						gll.sta=' ';
					}
				}

		}
		if (i==6) 
		{
			break;
		}
	}
	return true;
}

bool NMEA_process::Get_gpzda(CString line,gpzda &zda)
{
	CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(',');
	int star=line.Find('*');
	temp=line.Mid(1,5);
	zda.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}

		switch (i)
		{
			case 1 :
				{
					if (strlen!=0)
					{
						zda.time=_wtof(temp);
					}
					else
					{
						zda.time=-1;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						zda.day=_wtoi(temp);
					}
					else
					{
						zda.day=-1;
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						zda.month=_wtoi(temp);
					}
					else
					{
						zda.month=-1;
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						zda.year=_wtoi(temp);
					}
					else
					{
						zda.year=-1;
					}
				}
			case 5:
				{
					if (strlen!=0)
					{
						zda.loc_h=_wtoi(temp);
					}
					else
					{
						zda.loc_h=-100;
					}
				}
			case 6:
				{
					if (strlen!=0)
					{
						zda.loc_m=_wtoi(temp);
					}
					else
					{
						zda.loc_m=-100;
					}
				}
			
		}
		if (i==6)
		{
			break;
		}
	}
	return true;
}

bool NMEA_process::Get_gpgst(CString line,gpgst &gst)
{
	CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(',');
	int star=line.Find('*');
	temp=line.Mid(1,5);
	gst.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}

		switch (i)
		{
			case 1 :
				{
					if (strlen!=0)
					{
						gst.utc=_wtof(temp);
					}
					else
					{
						gst.utc=-1;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						gst.rms=_wtof(temp);
					}
					else
					{
						gst.rms=-1;
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						gst.std_maj=_wtof(temp);
					}
					else
					{
						gst.std_maj=-1;
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						gst.std_min=_wtof(temp);
					}
					else
					{
						gst.std_min=-1;
					}
				}
			case 5:
				{
					if (strlen!=0)
					{
						gst.maj_orien=_wtof(temp);
					}
					else
					{
						gst.maj_orien=-1;
					}
				}
				case 6:
				{
					if (strlen!=0)
					{
						gst.std_lat=_wtof(temp);
					}
					else
					{
						gst.std_lat=-1;
					}
				}
				case 7:
				{
					if (strlen!=0)
					{
						gst.std_lon=_wtof(temp);
					}
					else
					{
						gst.std_lon=-1;
					}
				}
				case 8:
				{
					if (strlen!=0)
					{
						gst.std_alt=_wtof(temp);
					}
					else
					{
						gst.std_alt=-1;
					}
				}
			
		}
		if (i==8) 
		{
			break;
		}
	}
	return true;
}

bool NMEA_process::Get_gngns(CString line,gngns &gns)
{
	CString temp;	char deli=',';
	int len=line.GetLength();
	int last_deli=line.ReverseFind(deli);
	int star=line.Find('*');
	temp=line.Mid(1,5);
	gns.GP_ID=temp;
	int n1=0;	int n2=0;	int strlen=0;
	for (int i=1;i<len;i++)
	{
		if (n2!=last_deli)
		{
			n1=line.Find(deli,n2);
			n2=line.Find(deli,n1+1);
			strlen=n2-n1-1;
			temp=line.Mid(n1+1,strlen);
		}
		else
		{
			strlen=star-last_deli-1;
			temp=line.Mid(n2+1,strlen);
		}

		switch (i)
		{
			case 1 :
				{
					if (strlen!=0)
					{
						gns.utc=_wtof(temp);
					}
					else
					{
						gns.utc=-1;
					}
				}
			case 2:
				{
					if (strlen!=0)
					{
						gns.lat=_wtof(temp);
					}
					else
					{
						gns.lat=-1;
					}
				}
			case 3:
				{
					if (strlen!=0)
					{
						gns.ns_ind=temp;
					}
					else
					{
						gns.ns_ind=' ';
					}
				}
			case 4:
				{
					if (strlen!=0)
					{
						gns.lon=_wtof(temp);
					}
					else
					{
						gns.lon=-1;
					}
				}
			case 5:
				{
					if (strlen!=0)
					{
						gns.ew_ind=temp;
					}
					else
					{
						gns.ew_ind=' ';
					}
				}
			case 6:
				{
					if (strlen!=0)
					{
						gns.mode_ind=_wtoi(temp);
					}
					else
					{
						gns.mode_ind=-1;
					}
				}
			case 7:
				{
					if (strlen!=0)
					{
						gns.num_sate=_wtoi(temp);
					}
					else
					{
						gns.num_sate=-1;
					}
				}
			case 8:
				{
					if (strlen!=0)
					{
						gns.hdop=_wtof(temp);
					}
					else
					{
						gns.hdop=-1;
					}
				}
			case 9:
				{
					if (strlen!=0)
					{
						gns.orth_height=_wtof(temp);
					}
					else
					{
						gns.orth_height=-1;
					}
				}
			case 10:
				{
					if (strlen!=0)
					{
						gns.geo_sep=_wtof(temp);
					
					}
					else
					{
						gns.geo_sep=-1;
					}
				}
			case 11:
				{
					if (strlen!=0)
					{
						gns.data_age=_wtof(temp);
					}
					else
					{
						gns.data_age=-1;
					}
				}
			case 12:
				{
					if (strlen!=0)
					{
						gns.ref_sta=temp;
					}
					else
					{
						gns.ref_sta=' ';
					}
				}
		}
		if (i==14) 
		{
			break;
		}
	}
	return true;
}

//add more in the future





//----------------------------Class process_CORS--------------------------------

bool process_CORS::read_rtcm(CString filename,EpochCorrInfo *SatCorr)
{
	CStdioFile infile;
	infile.Open(filename,CFile::modeRead);
	CString line;
	int epochcount=0;
	while(infile.ReadString(line))
	{
		//read the header of file
		if  (line.Find(_T("END OF HEADER"))>0)
		{
			break;
		}
	}
	while(infile.ReadString(line))
	{
				
			SatCorr[epochcount].GPSW=_wtoi(line.Mid(0,4));
			SatCorr[epochcount].GPSS=_wtof(line.Mid(6,14));
			SatCorr[epochcount].NumSat=_wtoi(line.Right(2));
			int satcount=0;
			int linecount=0;
			while (linecount<32) //  32 is the sat list maxnum in this epoch
			{
				epochcount++;
				infile.ReadString(line);
				if(line.GetLength()<10)
				{
					linecount++;
					continue;
				}
				else
				{
					if (line.Left(1)="G")
					{
						SatCorr[epochcount].CorrSatInfo[satcount].sysid=1;
					}
					else if (line.Left(1)="C")
					{
						SatCorr[epochcount].CorrSatInfo[satcount].sysid=5;
					}
					SatCorr[epochcount].CorrSatInfo[satcount].prn=_wtoi(line.Mid(1,2));
					SatCorr[epochcount].CorrSatInfo[satcount].corrCom=_wtoi(line.Mid(5,12));
					SatCorr[epochcount].CorrSatInfo[satcount].corrIono=_wtoi(line.Mid(19,12));
					SatCorr[epochcount].CorrSatInfo[satcount].corrTrop=_wtoi(line.Mid(32,12));
					// add more
					satcount++;
					linecount++;
				}
			}
		
	}

	return true;
}






//--------------------------------Class GNSSTBEXT----------------------------------

int GNSSTBext::DayofYear(int GpsWeek,double GpsSec)
 {
	 double MJD = GpsWeek*7.0+44244.0+GpsSec/3600.0/24.0;  //modified julian day
	 int a   =(int)(MJD +0.5+0.5+1.0e-10)+2400000;
	 double FRAC=MJD +0.5+2400000.5-a;
	 int b = a + 1537;
	 int Cc = (int)((b - 122.1)/365.25 + 1.0e-10);
	 int d = (int)(365.25*Cc+1.0e-16);
	 int Ee = (int)((b-d)/30.6001 + 1.0e-10);
	 int Day=b-d-(int)(30.6001*Ee);
	 int Mon=Ee-1-12*((int)(Ee/14.0+1.0e-10));
	 int Year=Cc - 4715 -(int)((7.0+Mon)/10.0 + 1.0e-10);

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
	 int doy=0;
	 for (int i=0;i<Mon-1;i++)
		 doy+=DayNums[i];
	 doy+=Day;
	 return doy;
 }


bool GNSSTBext::GPStoYMDHMS(int GPSWeek,double GPSSecond,int& Year,int& Month,int& Day,int& Hour,int& Minute,double& Second,
								double& JD,double& MJD,int& Weekday)
{
	MJD = GPSWeek*7.0+44244.0+GPSSecond/3600.0/24.0;  //modified julian day
	 int a   =(int)(MJD +0.5+0.5+1.0e-10)+2400000;
	 double FRAC=MJD +0.5+2400000.5-a;
	 int b = a + 1537;
	 int Cc = (int)((b - 122.1)/365.25 + 1.0e-10);
	 int d = (int)(365.25*Cc+1.0e-16);
	 int Ee = (int)((b-d)/30.6001 + 1.0e-10);
	 Day=b-d-(int)(30.6001*Ee);
	 Month=Ee-1-12*((int)(Ee/14.0+1.0e-10));
	 Year=Cc - 4715 -(int)((7.0+Month)/10.0 + 1.0e-10);
	 double TmHour = FRAC*24.0;
	 Hour= (int)(TmHour + 1.0e-10);
	 double TmMin=(TmHour-Hour)*60;
	 Minute=(int)(TmMin+1.0e-10);
	 Second= (TmMin-Minute)*60.0;	

	 JD=MJD+2400000.5;
	 Weekday=(((int)(JD+0.5))%7);//计算结果中0代表星期一，1代表星期二
	 Weekday++;
	 if (Weekday==7)
		 Weekday=0;

	 return true;
}

bool GNSSTBext::BLH2XYZ(double *BLH,double *XYZ)
	 {
		double B=*BLH;
		double L=*(BLH+1);
		double H=*(BLH+2);
		double N=RE_WGS84/sqrt( 1.0-ee*sin(B)*sin(B) );

		*XYZ=(N+H)*cos(B)*cos(L);
		*(XYZ+1)=(N+H)*cos(B)*sin(L);
		*(XYZ+2)=( N*(1-ee )+H )*sin(B);
		return true;
	 }
//	 GPtime to (hour min sec) or 
 void GNSSTBext::GPFormat2time1(double gp_utc,int& utc_hour,int & utc_min,double & utc_sec)
{
		int a =(int)(gp_utc/10000+0.0000001);
		utc_hour=a;
		int b= (int)(	(gp_utc-a*10000)/100	+0.0000001); 
		utc_min=b;
		double c=floor( (gp_utc-a*10000-b*100) );
		utc_sec=c;
}
//(day mon year)
 void GNSSTBext::GPFormat2time2(int gp_utc,int& utc_day,int& utc_mon,int& utc_year )
{
		int a =(int)(gp_utc/10000);
		utc_day=a;
		int b= (int)(	(gp_utc-a*10000)/100	); 
		utc_mon=b;
		int c=(int)( (gp_utc-a*10000-b*100));
		utc_year=c+2000;
}
	 
 void GNSSTBext::Gps_Week_Sec(int* Gps_Week,double* Gps_Second,int Y,int M,int D,int hour,int minute,double second)
{
	//返回值为GPS周秒
	
	double y;
	
    double m;
	if(M<=2)
	{
		y=Y-1.0;
		m=M+12.0;
	}
	else
	{
		y=Y;
		m=M;
	}
	double JD;
	JD=((int)(365.25*y))+((int)(30.6001*(m+1.0)))+D+hour/24.0+1720981.5;
	int weekday=(((int)(JD+0.5))%7);//计算结果中0代表星期一，1代表星期二
	 weekday++;
	 if (weekday==7)
		 weekday=0;
	 *Gps_Second=(weekday*24.0+hour)*60.0*60.0+minute*60.0+second;
	 *Gps_Week=(int)((JD-2444244.5)/7.0);
}

double GNSSTBext::GPFormat2Degree(double gpl)
	{
		double degree;
		double D	=floor(gpl/100);
		double M	=(gpl/100-D)*100/60;
		degree=D+M;
		return degree;
	}
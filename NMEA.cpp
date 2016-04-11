// NMEA.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"


int _tmain(int argc, _TCHAR* argv[])
{
	return 0;
}

//-----remove-----------------
	////method  1   char
	//char pFileName[90] = "GNSS.bmp"; 
	//CString fn;	
	//fn.Format(_T("%s"),pFileName);
	//fn=pFileName;
	////method 2  not rec
	////int charlen=strlen(pFileName);
	////int len=MultiByteToWideChar(CP_ACP,0,pFileName,charlen,NULL,0);
	////TCHAR* buf=new TCHAR(len+1);
	////MultiByteToWideChar(CP_ACP,0,pFileName,charlen,buf,len);
	//////
	////buf[len] = '\0';
	////CString pWideChar;
	////pWideChar.Append(buf); 
	////delete buf; 

	////---------------DLL-RTCM------------- Thread
	//WSADATA data;//initialize WSADATA struct
	//WORD	 w=MAKEWORD(1,1);
	//int err=WSAStartup(w,&data);
	//if(LOBYTE(data.wVersion)!=1 || HIBYTE(data.wVersion)!=1)
	//{
	//	WSACleanup();
	//}
	////Initialize the socket
	//SOCKET TJClient;
	//TJClient=socket(AF_INET,SOCK_STREAM,IPPROTO_TCP);
	//
	//sockaddr_in addr;
	//addr.sin_family=AF_INET;
	//addr.sin_port=htons(12102);
	//addr.sin_addr.S_un.S_addr=inet_addr("180.153.223.51");//

	//if(connect(TJClient,(SOCKADDR*)&addr,sizeof(addr)))
	//{
	//		//setup TimeInterval
	//	cout<<1234<<endl;
	//	
	//		
	//	
	//	
	//	
	//	
	//		MakeAuthReqMsg reqmsg;
	//		HINSTANCE Hlib=LoadLibrary(_T("rtcmin_dll.dll"));
	//		//int makeAuthReqMsg(const char* servicetype,const char* usekey=MDU3YmQx,char*authReqMsg)
	//		//int makeDGNSSDataReqMsg(const chat* ggaMsg,const char* rmcMsg,const chat* cellid,const char* lac,const chat* wifiMac,char*DGNSSDataReqMsg)
	//		//int GetDGNSSData(const char responseBuffer[],int responseBufferLen,char DGNSSData[],int *DGNSSDataLen,int rtcnType)
	//		if(Hlib==NULL)
	//		{
	//			FreeLibrary(Hlib);
	//		}
	//		MakeAuthReqMsg func1;
	//		func1=(MakeAuthReqMsg)GetProcAddress(Hlib, "makeAuthReqMsg");//without " _T( ) "
	//		if(func1==NULL)
	//		{
	//			FreeLibrary(Hlib);
	//		}
	//		//		Service type , userkey ...see "RTCM-IN "
	//		const char SerType[]="TJ3N";
	//		const char UseKey[]="057bd1:tj";
	//		char authReqMsg[100];
	//		func1(SerType,UseKey,authReqMsg);
	//		FreeLibrary(Hlib);
	//}
	//closesocket(TJClient);
	//WSACleanup();
	////	process data
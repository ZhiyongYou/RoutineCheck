#include <string.h>
#include <cstring>
#include <string>
#include "ReadFile.h"
#include "camera.h"


using namespace std;

ReadFile::ReadFile()
{

}

ReadFile::~ReadFile()
{

}

bool ReadFile::OpenFile(const char* f_name)
{
	SetFileName(f_name);
	return OpenFile();
}

bool ReadFile::OpenFile()
{
	infp = fopen(filename.c_str(),"r");
	if(!infp)
	{
		cerr<<"Error opening file "<<filename<<endl;
		return false;
	}
	else
	{
		cout<<"Open File "<<filename<<endl;
		return true;
	}
}

void ReadFile::ReadInfo(double *conf, const char* conf_type)
{
	if(!strcmp(conf_type,"hv")){
		cout<<"read hv information"<<endl;
		mreadInfo(conf);
	}
	if(!strcmp(conf_type,"base")){
		cout<<"read base information"<<endl;
		char threshtag[100];
		fscanf(infp,"subcluster_id thresholds\n",threshtag);
		mreadInfo(conf);
	}

}

void ReadFile::mreadInfo(double *conf)
{
	short isubcluster,mSiPM;
	while(!feof(infp)){
		fscanf(infp,"0x%d %lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
				&isubcluster,&conf_info[0],&conf_info[1],&conf_info[2],&conf_info[3],&conf_info[4],&conf_info[5],&conf_info[6],&conf_info[7],&conf_info[8],
				&conf_info[9],&conf_info[10],&conf_info[11],&conf_info[12],&conf_info[13],&conf_info[14],&conf_info[15]);
		for(int i=0; i<16; i++) {
			SC_Channel2SiPM(isubcluster, i+1, &mSiPM);
			//if(tel==5) mSiPM = 1023-mSiPM;
			conf[mSiPM] = conf_info[i];
		}
	}
	//for(int i=0;i<1024;i++)
	//{
	//	printf("%d %lf\n",i,conf[i]);
	//}
}

void ReadFile::GetNoisePix(int iTel, vector<int> &noisSiPM, vector<int> &edge_noisSiPM)
{
	noisSiPM.clear();
	edge_noisSiPM.clear();
	int itel;
	char cen_pixels[5000];
	char edge_pixels[5000];
	char *cen_pix;
	char *edge_pix;

	while(!feof(infp)){
		fscanf(infp,"tel%02d:cen:%s edge:%s\n",&itel,cen_pixels,edge_pixels);
		if(itel==iTel){
			cen_pix = strtok(cen_pixels, ",");
			if(atoi(cen_pix)!=-1)	{noisSiPM.push_back( atoi(cen_pix) );}
			while((cen_pix = strtok(NULL, ",")))
			{
				if(atoi(cen_pix)!=-1)	{noisSiPM.push_back( atoi(cen_pix) );}
			}

			edge_pix = strtok(edge_pixels, ",");
			if(atoi(edge_pix)!=-1)	{edge_noisSiPM.push_back( atoi(edge_pix) );}
			while((edge_pix = strtok(NULL, ",")))
			{
				if(atoi(edge_pix)!=-1)	{edge_noisSiPM.push_back( atoi(edge_pix) );}
			}
		}
	}
/*
	switch(iTel)
	{
		case 1:
			break;
		case 2:
			noisSiPM.push_back(925);//SC:86 Channel:3 sipm:925
			noisSiPM.push_back(735);//SC:85 Channel:15 sipm:735
			edge_noisSiPM.push_back(416);
			edge_noisSiPM.push_back(448);
			edge_noisSiPM.push_back(449);
			edge_noisSiPM.push_back(480);
			edge_noisSiPM.push_back(481);
			edge_noisSiPM.push_back(512);
			edge_noisSiPM.push_back(513);
			edge_noisSiPM.push_back(544);
			edge_noisSiPM.push_back(545);
			edge_noisSiPM.push_back(576);
			edge_noisSiPM.push_back(577);
			edge_noisSiPM.push_back(735);//SC:85 Channel:15 sipm:735
			break;
		case 3:
			break;
		case 4:
			noisSiPM.push_back(96);//SC:82 Channel:10 sipm:96
			noisSiPM.push_back(60);
			noisSiPM.push_back(134);
			noisSiPM.push_back(215);
			noisSiPM.push_back(228);
			noisSiPM.push_back(251);
			noisSiPM.push_back(353);
			noisSiPM.push_back(574);
			noisSiPM.push_back(614);
			noisSiPM.push_back(740);
			noisSiPM.push_back(752);
			noisSiPM.push_back(761);
			noisSiPM.push_back(798);
			noisSiPM.push_back(827);
			noisSiPM.push_back(892);
			noisSiPM.push_back(944);
			noisSiPM.push_back(945);
			edge_noisSiPM.push_back(192);
			edge_noisSiPM.push_back(193);
			edge_noisSiPM.push_back(224);
			edge_noisSiPM.push_back(225);
			edge_noisSiPM.push_back(256);
			edge_noisSiPM.push_back(257);
			edge_noisSiPM.push_back(288);
			break;
		case 5:
			noisSiPM.push_back(311);//SC:28 Channel:9 sipm:712
			noisSiPM.push_back(645);
			noisSiPM.push_back(790);
			break;
		case 6:
			break;
		case 7:
			break;
		case 10:
			break;
	}
	*/
}

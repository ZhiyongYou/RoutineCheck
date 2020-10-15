#include <iostream>
#include <fstream>
#include <algorithm>
#include <string.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "WFCTASLC.h"
#include "astro.h"

WFCTASLC::WFCTASLC()
{
	Clear();
}

WFCTASLC::~WFCTASLC()
{

}

void WFCTASLC::Clear()
{
	for(int tel=1;tel<20;tel++)
	{
		slc_status[tel-1].clear();
	}
}

void WFCTASLC::ReadSlcFile(int year, int month, int day)
{
	for(int tel=1;tel<20;tel++)
	{
		m_ReadSlcFile(year, month, day, tel);
	}
}
void WFCTASLC::m_ReadSlcFile(int year, int month, int day, int tel)
{
	std::string slc_file_dir = "/scratchfs/ybj/ketong/workTime/StatisticalTime/outFile/data/separate";
	std::string s_year, s_month, s_day, s_tel;
	s_year = std::to_string(year);
	month<10 ? s_month = "0" + std::to_string(month) : s_month = std::to_string(month);
	day<10 ? s_day = "0" + std::to_string(day) : s_day = std::to_string(day);
	tel<10 ? s_tel = "0" + std::to_string(tel) : s_tel = std::to_string(tel);

	std::string slc_file = slc_file_dir + "/" + s_year + "_" + s_month + "_" + s_day + "_" + s_tel + ".root";
	printf("slc file: %s\n",slc_file.c_str());

	TFile *slcfile = TFile::Open( slc_file.c_str(), "read" );
	// Check file
	if(!slcfile) {
		printf("%s wfcta slc file does note exit!!\n",slc_file.c_str());
		return ;
	}
	if(slcfile->IsZombie()||slcfile->GetEND()<50) {
		printf("%s wfcta slc file error!!\n",slc_file.c_str());
		slcfile->Close();
		return ;
	}
	// Check tree
	TTree *t = (TTree *)slcfile->Get("t");
	if(!t){
		printf("  The SLC file %s tree error!!\n",slc_file.c_str());
		slcfile->Close();
		return ;
	}
	double time_stamp;
	float slow_control_stu;
	float temperature_stu;
	float daq_stu;
	float cloud_stu;
	float phase;
	float delta_angle;
	t->SetBranchAddress("time_stamp",&time_stamp);
	t->SetBranchAddress("slow_control_stu",&slow_control_stu);
	t->SetBranchAddress("temperature_stu",&temperature_stu);
	t->SetBranchAddress("daq_stu",&daq_stu);
	t->SetBranchAddress("cloud_stu",&cloud_stu);
	t->SetBranchAddress("phase",&phase);
	t->SetBranchAddress("delta_angle",&delta_angle);
	int slc_events = t->GetEntries();
	for(int ientry=0;ientry<slc_events;ientry++)
	{
		t->GetEntry(ientry);

		SlcStatus slcstatus;
		slcstatus.daq_stu = daq_stu;
		slcstatus.slow_control_stu = slow_control_stu;
		slcstatus.temperature_stu = temperature_stu;
		slcstatus.cloud_stu = cloud_stu;
		slcstatus.door_angle = delta_angle;
		slcstatus.r_led = 0;
		slcstatus.d_led = 1;

		//printf("%lf %f %f %f %f %f %f    %d\n",time_stamp,slow_control_stu,temperature_stu,daq_stu,cloud_stu,phase,delta_angle,status_no);
		slc_status[tel-1].push_back(std::make_pair( time_stamp, slcstatus ) );
	}

	slcfile->Close();

}

void WFCTASLC::ReadSlcFileYMJ(int year, int month, int day)
{
	for(int tel=1;tel<20;tel++)
	{
		m_ReadSlcFileYMJ(year, month, day, tel);
	}
}
void WFCTASLC::m_ReadSlcFileYMJ(int year, int month, int day, int tel)
{
	std::string slc_file_dir = "root://eos01.ihep.ac.cn//eos/user/w/wfcta/wfcta/wfcta/slowctrl/rootdata";
	std::string s_year, s_month, s_day, s_tel;
	s_year = std::to_string(year);
	s_month = std::to_string(month);
	s_day = std::to_string(day);
	s_tel = std::to_string(tel);

	std::string slc_file = slc_file_dir + "/" + s_year + "-" + s_month + "-" + s_day + "_" + s_tel + ".root";
	printf("slc file: %s\n",slc_file.c_str());

	TFile *slcfile = TFile::Open( slc_file.c_str(), "read" );
	// Check file
	if(!slcfile) {
		printf("%s wfcta slc file does note exit!!\n",slc_file.c_str());
		return ;
	}
	if(slcfile->IsZombie()||slcfile->GetEND()<50) {
		printf("%s wfcta slc file error!!\n",slc_file.c_str());
		slcfile->Close();
		return ;
	}
	// Check tree
	TTree *Tel01_info = (TTree *)slcfile->Get("Tel01_info");
	if(!Tel01_info){
		printf("  The SLC file %s tree error!!\n",slc_file.c_str());
		slcfile->Close();
		return ;
	}
	double T1MJD;
	int T1DledSt;
	int T1RledSt;
	double T1Door1Deg;
	double T1Door2Deg;
	Tel01_info->SetBranchAddress("T1MJD",&T1MJD);
	Tel01_info->SetBranchAddress("T1DledSt",&T1DledSt);
	Tel01_info->SetBranchAddress("T1RledSt",&T1RledSt);
	Tel01_info->SetBranchAddress("T1Door1Deg",&T1Door1Deg);
	Tel01_info->SetBranchAddress("T1Door2Deg",&T1Door2Deg);
	int slc_events = Tel01_info->GetEntries();
	for(int ientry=0;ientry<slc_events-1;ientry++)
	{
		Tel01_info->GetEntry(ientry);

		SlcStatus slcstatus;
		slcstatus.daq_stu = 1;
		slcstatus.slow_control_stu = 1;
		slcstatus.temperature_stu = 1;
		slcstatus.cloud_stu = 1;
		if( T1Door1Deg>180 && T1Door2Deg>180 )	{	slcstatus.door_angle = 1;}
		else if( T1Door1Deg<5 && T1Door2Deg<5 )	{	slcstatus.door_angle = 0;}
		else									{	slcstatus.door_angle = -1;}
		slcstatus.r_led = T1RledSt;
		slcstatus.d_led = T1DledSt;

		int64_t rbtime, rabbitTime;
		double rabbittime;
		mjd2rbtime(T1MJD, &rbtime, &rabbitTime, &rabbittime);
		double time_stamp = double(rabbitTime);
//		printf("%d : %lf %lf %d %d %d %lf\n",ientry, time_stamp, T1MJD, T1DledSt, T1RledSt, T1Door2Deg, T1Door1Deg);
		slc_status[tel-1].push_back(std::make_pair( time_stamp, slcstatus ) );
	}

	slcfile->Close();

}

void WFCTASLC::GetSlcStatus(long long rbTime, int tel_id, SlcStatus& slcstatus)
{
	int nsta = slc_status[tel_id-1].size();
	if( rbTime < slc_status[tel_id-1].at(0).first )
	{
		printf("SLC start time is later than Event time!\n");
	}
	else if( rbTime >= slc_status[tel_id-1].at(nsta-1).first )
	{
		printf("SLC end time is before than Event time!\n");
	}
	else 
	{
		for( int ista=0; ista<nsta-1; ista++)
		{
			if( rbTime >= slc_status[tel_id-1].at(ista).first && rbTime < slc_status[tel_id-1].at(ista+1).first )
			{
				slcstatus = slc_status[tel_id-1].at(ista).second;
				break;
			}
		}
	}
}

double WFCTASLC::GetDaq_Stu(long long rbTime, int tel_id)
{
	SlcStatus slcstatus;
	GetSlcStatus(rbTime, tel_id, slcstatus);
	return slcstatus.daq_stu;
}
double WFCTASLC::GetSlow_Control_Stu(long long rbTime, int tel_id)
{
	SlcStatus slcstatus;
	GetSlcStatus(rbTime, tel_id, slcstatus);
	return slcstatus.slow_control_stu;
}
double WFCTASLC::GetTemperature_Stu(long long rbTime, int tel_id)
{
	SlcStatus slcstatus;
	GetSlcStatus(rbTime, tel_id, slcstatus);
	return slcstatus.temperature_stu;
}
double WFCTASLC::GetCloud_Stu(long long rbTime, int tel_id)
{
	SlcStatus slcstatus;
	GetSlcStatus(rbTime, tel_id, slcstatus);
	return slcstatus.cloud_stu;
}
double WFCTASLC::GetDoor_Angle(long long rbTime, int tel_id)
{
	SlcStatus slcstatus;
	GetSlcStatus(rbTime, tel_id, slcstatus);
	return slcstatus.door_angle;
}
int WFCTASLC::GetR_LED(long long rbTime, int tel_id)
{
	SlcStatus slcstatus;
	GetSlcStatus(rbTime, tel_id, slcstatus);
	return slcstatus.r_led;
}
int WFCTASLC::GetD_LED(long long rbTime, int tel_id)
{
	SlcStatus slcstatus;
	GetSlcStatus(rbTime, tel_id, slcstatus);
	return slcstatus.d_led;
}






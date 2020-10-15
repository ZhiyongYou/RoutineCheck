#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <array>
#include "time.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "WFCTAEvent.h"
#include "LedEvtSel.h"
#include "WFCTASLC.h"
#include "WFCTARec.h"
#include "astro.h"

using namespace std;
#define NTELS 6
#define laser2_rb_time 990035000 //990010000 - 990060000 ns
#define laser3_rb_time 990850000 //990820000 - 990880000 ns
#define Ntotal 360000
const double UPPERCUT = 6000;
const double LOWERCUT = 3000;

WFCTAEvent *wfctaevent = new WFCTAEvent();
WFCTASLC* wfctaslc = new WFCTASLC();
LedEvtSel *ledsel = new LedEvtSel();
WFCTARec* wfctarec = new WFCTARec();

void release_globle();

int main(int argc, char *argv[])
{
	if(argc<0) {
		printf("Usage: %s outroot wfcta01.txt wfcta02.txt wfcta03.txt wfcta04.txt wfcta05.txt wfcta06.txt\n", argv[0]);
		return 1;
	}
	cout<<"begin time:"<<time(0)<<endl;
	double ImageX[1024]={0};
	double ImageY[1024]={0};
	//set map for led correction
	int PIX = 32;
	double D_ConeOut=25.8 ;// mm
	double centerx, centery;
	centerx = 414.3;
	centery = 414.3;
	for(int k=0;k<1024;k++){
		int i = k/32;
		int j = k%32;
		if(i%2==0)
			ImageX[k] = ((j+0.5)*D_ConeOut - centerx);// + interval*Interx-centerx);
		if(i%2==1)
			ImageX[k] = ((j+1)*D_ConeOut - centerx);// + interval*Interx-centerx);
		ImageY[k] = ((PIX-i)*D_ConeOut - centery);
	}

	///////////////////////////////////////////////////////
	//              set output file                 //
	///////////////////////////////////////////////////////
	char Name1[300]="root://eos01.ihep.ac.cn/";
	char of_name[300];
	strcpy(of_name,Name1);
	strcat(of_name,argv[2]);
	strcat(of_name,".root");
	TFile *outfile = TFile::Open(of_name,"recreate");

	int sipm[1024];
	for(int i=0;i<1024;i++){sipm[i] = i;}

	double dark_base_rms[1024]={0};
	double dark_base_ave[1024]={0};
	int dark_base_n[1024]={0};

	double base_ave[1024]={0};
	double base_rms[1024]={0};
	int base_n[1024]={0};

	int trigger_cr[1024]={0};
	int trigger_led[1024]={0};
	int trigger_R_led[1024]={0};
	int trigger_laser2[1024]={0};
	int trigger_laser3[1024]={0};

	double pe_cr[1024]={0};
	double pe_led[1024]={0};
	double pe_R_led[1024]={0};
	double pe_laser2[1024]={0};
	double pe_laser3[1024]={0};

	double h_l_ratio_ave[1024]={0};
	double h_l_ratio_rms[1024]={0};
	int h_l_ratio_n[1024]={0};

	double obs_time_cr=0;
	double obs_time_led=0;
	double obs_time_R_led=0;
	double obs_time_laser2=0;
	double obs_time_laser3=0;

	double totaltrigger_cr=0;
	double totaltrigger_led=0;
	double totaltrigger_R_led=0;
	double totaltrigger_laser2=0;
	double totaltrigger_laser3=0;
	TTree *singleChChk = new TTree("singleChChk","check info of events");
	singleChChk->Branch("sipm",sipm,"sipm[1024]/I");

	singleChChk->Branch("dark_base_rms",dark_base_rms,"dark_base_rms[1024]/D");
	singleChChk->Branch("dark_base_ave",dark_base_ave,"dark_base_ave[1024]/D");
	singleChChk->Branch("dark_base_n",dark_base_n,"dark_base_n[1024]/I");

	singleChChk->Branch("base_rms",base_rms,"base_rms[1024]/D");
	singleChChk->Branch("base_ave",base_ave,"base_ave[1024]/D");
	singleChChk->Branch("base_n",base_n,"base_n[1024]/I");

	singleChChk->Branch("h_l_ratio_ave",h_l_ratio_ave,"h_l_ratio_ave[1024]/D");
	singleChChk->Branch("h_l_ratio_rms",h_l_ratio_rms,"h_l_ratio_rms[1024]/D");
	singleChChk->Branch("h_l_ratio_n",h_l_ratio_n,"h_l_ratio_n[1024]/I");

	singleChChk->Branch("trigger_cr",trigger_cr,"trigger_cr[1024]/I");
	singleChChk->Branch("trigger_led",trigger_led,"trigger_led[1024]/I");
	singleChChk->Branch("trigger_R_led",trigger_R_led,"trigger_R_led[1024]/I");
	singleChChk->Branch("trigger_laser2",trigger_laser2,"trigger_laser2[1024]/I");
	singleChChk->Branch("trigger_laser3",trigger_laser3,"trigger_laser3[1024]/I");

	singleChChk->Branch("pe_cr",pe_cr,"pe_cr[1024]/D");
	singleChChk->Branch("pe_R_led",pe_R_led,"pe_R_led[1024]/D");
	singleChChk->Branch("pe_led",pe_led,"pe_led[1024]/D");
	singleChChk->Branch("pe_laser2",pe_laser2,"pe_laser2[1024]/D");
	singleChChk->Branch("pe_laser3",pe_laser3,"pe_laser3[1024]/D");

	singleChChk->Branch("obs_time_cr",&obs_time_cr,"obs_time_cr/D");
	singleChChk->Branch("obs_time_led",&obs_time_led,"obs_time_led/D");
	singleChChk->Branch("obs_time_R_led",&obs_time_R_led,"obs_time_R_led/D");
	singleChChk->Branch("obs_time_laser2",&obs_time_laser2,"obs_time_laser2/D");
	singleChChk->Branch("obs_time_laser3",&obs_time_laser3,"obs_time_laser3/D");

	singleChChk->Branch("totaltrigger_cr",&totaltrigger_cr,"totaltrigger_cr/D");
	singleChChk->Branch("totaltrigger_led",&totaltrigger_led,"totaltrigger_led/D");
	singleChChk->Branch("totaltrigger_R_led",&totaltrigger_R_led,"totaltrigger_R_led/D");
	singleChChk->Branch("totaltrigger_laser2",&totaltrigger_laser2,"totaltrigger_laser2/D");
	singleChChk->Branch("totaltrigger_laser3",&totaltrigger_laser3,"totaltrigger_laser3/D");

	double mjd;
	double door_angle;
	int r_led;
	int d_led;
	int packageCheck,eventNumberCheck;
	int event_type;
	double deltaT;
	double deltaT5;
	int Npix;
	double Size;
	double MeanX, MeanY;
	TTree *eventChk = new TTree("eventChk","info of eventChk");
	eventChk->Branch("rabbitTime",&wfctaevent->rabbitTime,"rabbitTime/L");
	eventChk->Branch("rabbittime",&wfctaevent->rabbittime,"rabbittime/D");
	eventChk->Branch("mjd",&mjd,"mjd/D");
	eventChk->Branch("door_angle",&door_angle,"door_angle/D");
	eventChk->Branch("r_led",&r_led,"r_led/I");
	eventChk->Branch("d_led",&d_led,"d_led/I");
	eventChk->Branch("packageCheck",&packageCheck,"packageCheck/I"); // 0 is for OK; 1 is for bad packageCheck
	eventChk->Branch("eventNumberCheck",&eventNumberCheck,"eventNumberCheck/I");
	eventChk->Branch("event_type",&event_type,"event_type/I"); // 0 is cr, 1 is led, 2 is reflect led, 3 is noise, 4 is laser2, 5 is laser3
	eventChk->Branch("deltaT",&deltaT,"deltaT/D");
	eventChk->Branch("deltaT5",&deltaT5,"deltaT5/D");
	eventChk->Branch("Npix",&Npix,"Npix/I");
	eventChk->Branch("Size",&Size,"Size/D");
	eventChk->Branch("MeanY",&MeanY,"MeanY/D");
	eventChk->Branch("MeanX",&MeanX,"MeanX/D");

	//////////////////////////////////////////////////
	//			Read WFCTA SLC File		         //
	//////////////////////////////////////////////////
	int year=atoi(argv[3]);
	int mon=atoi(argv[4]);
	int day=atoi(argv[5]);

	struct tm stm;
	stm.tm_year = year - 1900;
	stm.tm_mon = mon-1;
	stm.tm_mday = day;
	stm.tm_hour = 0;
	stm.tm_min = 0;
	stm.tm_sec = 0;
	time_t utc = mktime(&stm);
	utc += 86400;
	struct tm *stm_slc = localtime(&utc);
	wfctaslc->ReadSlcFileYMJ(year, mon, day);
	wfctaslc->ReadSlcFileYMJ(stm_slc->tm_year+1900, stm_slc->tm_mon+1, stm_slc->tm_mday);


	FILE *fp = fopen(argv[1],"r");
	char telfile[500]; strcpy(telfile,"empty");
	char telFile[500];
	while(!feof(fp))
	{
		fscanf(fp,"%s\n",telfile);
		strcpy(telFile,"root://eos01.ihep.ac.cn/");
		strcat(telFile,telfile);
		TFile* inFile = TFile::Open(telFile,"READ");
		printf("inputing wfcta file -> %s\n",telFile);
		if(!inFile || inFile->IsZombie() || inFile->GetEND()<50)
		{
			printf("%s file error!!\n",telFile);
			continue;
		}
		TTree* telTree = (TTree *)inFile->Get("eventShow");
		if(telTree==nullptr)
		{
			printf("%s is null file\n",telFile);
			continue;
		}
		telTree->SetBranchAddress("WFCTAEvent",&wfctaevent);

		long Nwfctaevts = telTree->GetEntries();
		printf("Entries: %ld\n", Nwfctaevts);

		//loop events
		Long64_t rb_Timecr=0;		double rb_timecr=0;
		Long64_t rb_Timecr5=0;		double rb_timecr5=0;
		Long64_t rb_Timeled=0;		double rb_timeled=0;
		Long64_t rb_Time_Rled=0;	double rb_time_Rled=0;
		Long64_t rb_Timelaser2=0;	double rb_timelaser2=0;
		Long64_t rb_Timelaser3=0;	double rb_timelaser3=0;
		telTree->GetEntry(0);
		ledsel->SetLed0( wfctaevent->ITel() );
		for(int ientry=0;ientry<Nwfctaevts;ientry++)
		{
			//init
			mjd = -1000;packageCheck=0;eventNumberCheck=0;event_type=-1000;deltaT=-1000;deltaT5=-1000;Npix=-1000;Size=-1000;MeanX=-1000;MeanY=-1000;

			//get entry
			telTree->GetEntry(ientry);
//			if(0==ientry%10000)	std::cerr << "entry:" << ientry << "/" << Nwfctaevts <<std::endl;

			//get wfcta slc status/////////////////////////////////////////////////////////////////////////////
			door_angle = wfctaslc->GetDoor_Angle(wfctaevent->rabbitTime, wfctaevent->iTel);
			r_led = wfctaslc->GetR_LED(wfctaevent->rabbitTime, wfctaevent->iTel);
			d_led = wfctaslc->GetD_LED(wfctaevent->rabbitTime, wfctaevent->iTel);
			//	printf("entry:%d rabbitTime:%lld Door_Angle: %lf R_LED:%d D_LED:%d\n", ientry, wfctaevent->rabbitTime, door_angle, r_led, d_led);
			//get mjd//////////////////////////////////////////////////////////////////////////////////////////
			mjd = rbtime2mjd(wfctaevent->rabbitTime, wfctaevent->rabbittime);
			//get packcheck and eventcheck/////////////////////////////////////////////////////////////////////
			for(int ii=0;ii<wfctaevent->iSiPM.size();ii++) {
				if(wfctaevent->eevent.at(ii)!=-1) {
					int delta_evt = wfctaevent->eevent.at(ii) - (wfctaevent->eEvent)%1024;
					if(delta_evt!=0) {
						eventNumberCheck=1; //bad
					}
				}
			}
			for(int ii=0;ii<wfctaevent->packCheck.size();ii++) {
				if(wfctaevent->packCheck.at(ii)!=1000) {
					packageCheck=1;
				}
			}
			if(packageCheck==1) {
				eventChk->Fill();
				continue;
			}

			//get event type/////////////////////////////////////////////////////////////////////////////////////
			if( abs(wfctaevent->Rabbittime()*20-laser2_rb_time)<=25000 ) {
				event_type = 4;  //l2
				wfctarec->SetEventType(1);
			}
			else if( abs(wfctaevent->Rabbittime()*20-laser3_rb_time)<=30000 ) {
				event_type = 5;  //l3
				wfctarec->SetEventType(1);
			}
			else {
				if(1==r_led) {
					event_type = 2; //reflect led
					wfctarec->SetEventType(1);
				}
				else {
					if( 13==ledsel->Running(wfctaevent) ) {
						event_type = 1;  //led
						wfctarec->SetEventType(1);
					}
					else if(0==door_angle) {
						event_type = 3;  //noise
						wfctarec->SetEventType(1);
					}
					else{	
						event_type = 0;  //cosmic ray
						wfctarec->SetEventType(0);
					}
				}
			}
			//get rec parameters/////////////////////////////////////////////////////////////////////////////////////
					// get isipm sipmpe sipmt
			std::vector<int> isipm;		isipm.clear();
			std::vector<double> sipmpe; sipmpe.clear();
			std::vector<double> sipmt;	sipmt.clear();
			int tel_id = wfctaevent->iTel;
			for(int ii=0;ii<wfctaevent->iSiPM.size();ii++) {
				int sipm_id = wfctaevent->iSiPM.at(ii);
				double pe=0;
				double peak_t=0;
				if(wfctaevent->SatH.at(ii)==0 && wfctaevent->eSatH.at(ii)==0 && wfctaevent->AdcH.at(ii)<8000) {//not satuate
					if(event_type>3)    {   pe = wfctaevent->LaserAdcH.at(ii)/9.98;}
					else            {   pe = wfctaevent->AdcH.at(ii)/9.98;}
					peak_t = wfctaevent->PeakPosH.at(ii)*80;
				}
				else {
					if(event_type>3)    {   pe = wfctaevent->LaserAdcL.at(ii)*22/9.98;}
					else            {   pe = wfctaevent->AdcL.at(ii)*22/9.98;}
					peak_t = wfctaevent->PeakPosL.at(ii)*80;
				}
				if(pe<0)          {  pe = 0;}
				else if(pe>Ntotal){  pe = Ntotal;}
				else              {  pe = -Ntotal*log(1-pe/Ntotal);}

				//led pe correction
				if(event_type==1||event_type==2) {
					double theta = pow(cos(sqrt(ImageX[sipm_id]*ImageX[sipm_id]+ImageY[sipm_id]*ImageY[sipm_id])/2870),4);
					pe = pe/theta;
				}
				else {
					double intensity=1;
					int new_sipm = wfctarec->GetSipmIdInSimulation(sipm_id);
					CollectIntensity::Instance()->GetIntensity(new_sipm,intensity);
					pe = pe / intensity;
				}

				isipm.push_back(sipm_id+1024*(tel_id-1));
				sipmpe.push_back(pe);
				sipmt.push_back(peak_t);


			}
			//wfcta reconstruction
			wfctarec->SetWFCTAEvent(isipm, sipmpe, sipmt);
			wfctarec->TimeClean(50);
			wfctarec->IslandClean();
			wfctarec->CalcMainTel(1);
			wfctarec->MergeEvent();
			wfctarec->CalcHillas();

			Npix = wfctarec->GetNpix();
			Size = wfctarec->GetSize();
			MeanX = wfctarec->GetMeanX();
			MeanY = wfctarec->GetMeanY();
			if( 0!=event_type && Npix<=50 ) {
				event_type = 3;
			}
			if(0==event_type) {
				deltaT = (wfctaevent->rabbitTime-rb_Timecr) + (wfctaevent->rabbittime-rb_timecr)*20/1000000000.;
				if(deltaT>0&&deltaT<2)  obs_time_cr += deltaT;
				rb_Timecr = wfctaevent->rabbitTime;
				rb_timecr = wfctaevent->rabbittime;
				if(Npix>5)
				{
					deltaT5 = (wfctaevent->rabbitTime-rb_Timecr5) + (wfctaevent->rabbittime-rb_timecr5)*20/1000000000.;
					rb_Timecr5 = wfctaevent->rabbitTime;
					rb_timecr5 = wfctaevent->rabbittime;
				}
			}
			if(1==event_type) {
				deltaT = (wfctaevent->rabbitTime-rb_Timeled) + (wfctaevent->rabbittime-rb_timeled)*20/1000000000.;
				if(deltaT>0&&deltaT<2)  obs_time_led += deltaT;
				rb_Timeled = wfctaevent->rabbitTime;
				rb_timeled = wfctaevent->rabbittime;
			}
			if(2==event_type) {
				deltaT = (wfctaevent->rabbitTime-rb_Time_Rled) + (wfctaevent->rabbittime-rb_time_Rled)*20/1000000000.;
				if(deltaT>0&&deltaT<2)  obs_time_R_led += deltaT;
				rb_Time_Rled = wfctaevent->rabbitTime;
				rb_time_Rled = wfctaevent->rabbittime;
			}
			if(4==event_type) {
				deltaT = (wfctaevent->rabbitTime-rb_Timelaser2) + (wfctaevent->rabbittime-rb_timelaser2)*20/1000000000.;
				if(deltaT>0&&deltaT<2)  obs_time_laser2 += deltaT;
				rb_Timelaser2 = wfctaevent->rabbitTime;
				rb_timelaser2 = wfctaevent->rabbittime;
			}
			if(5==event_type) {
				deltaT = (wfctaevent->rabbitTime-rb_Timelaser3) + (wfctaevent->rabbittime-rb_timelaser3)*20/1000000000.;
				if(deltaT>0&&deltaT<2)  obs_time_laser3 += deltaT;
				rb_Timelaser3 = wfctaevent->rabbitTime;
				rb_timelaser3 = wfctaevent->rabbittime;
			}
			eventChk->Fill();
			if(Npix<=5) continue;

			//////////////////////////////////////////////////
			// singleChChk Tree Calc//
			//////////////////////////////////////////////////
			std::vector<int> clean_sipm;		clean_sipm.clear();
			std::vector<double> clean_sipmpe;	clean_sipmpe.clear();
			std::vector<double> clean_sipmt;	clean_sipmt.clear();
			wfctarec->GetCleanImage(clean_sipm, clean_sipmpe, clean_sipmt);
			for(int ii=0; ii<clean_sipm.size();ii++)
			{
				int sipm_id = clean_sipm.at(ii)%1024;
				double pe = clean_sipmpe.at(ii);
				if(0==event_type) {	trigger_cr[sipm_id]++;		pe_cr[sipm_id] += pe;}
				if(1==event_type) {	trigger_led[sipm_id]++;		pe_led[sipm_id] += pe;}
				if(2==event_type) {	trigger_R_led[sipm_id]++;	pe_R_led[sipm_id] += pe;}
				if(4==event_type)  {	trigger_laser2[sipm_id]++;	pe_laser2[sipm_id] += pe;}
				if(5==event_type)  {	trigger_laser3[sipm_id]++;	pe_laser3[sipm_id] += pe;}
			}
			if(0==event_type) {	totaltrigger_cr++;}
			if(1==event_type) {	totaltrigger_led++;}
			if(2==event_type) {	totaltrigger_R_led++;}
			if(4==event_type) {	totaltrigger_laser2++;}
			if(5==event_type) {	totaltrigger_laser3++;}

			//high low gain ratio
			for(int ii=0;ii<wfctaevent->iSiPM.size();ii++)
			{
				if(0==event_type)
				{
					double high_adc = wfctaevent->AdcH.at(ii);
					double low_adc = wfctaevent->AdcL.at(ii);
					int isipm = wfctaevent->iSiPM.at(ii);
					if(high_adc<UPPERCUT&&high_adc>LOWERCUT){
						double ratio = high_adc/low_adc;
						h_l_ratio_ave[isipm] += ratio;
						h_l_ratio_rms[isipm] += pow(ratio,2);
						h_l_ratio_n[isipm]++;
					}
				}
			}
		}
		inFile->Close();
	}
	fclose(fp);

	for(int i=0;i<1024;i++)
	{
		if(0==h_l_ratio_n[i]) continue;
		h_l_ratio_ave[i] /= h_l_ratio_n[i];
		h_l_ratio_rms[i] /= h_l_ratio_n[i];
		h_l_ratio_rms[i] = sqrt(h_l_ratio_rms[i] - pow(h_l_ratio_ave[i],2));
	}
	for(int i=0;i<1024;i++)
	{
		trigger_cr[i]!=0 ? pe_cr[i] /= trigger_cr[i] : pe_cr[i]=0;
		trigger_led[i]!=0 ? pe_led[i] /= trigger_led[i] : pe_led[i]=0;
		trigger_R_led[i]!=0 ? pe_R_led[i] /= trigger_R_led[i] : pe_R_led[i]=0;
		trigger_laser2[i]!=0 ? pe_laser2[i] /= trigger_laser2[i] : pe_laser2[i]=0;
		trigger_laser3[i]!=0 ? pe_laser3[i] /= trigger_laser3[i] : pe_laser3[i]=0;
	}
	singleChChk->Fill();

	cout<<"deal time:"<<time(0)<<endl<<endl;

	//////////////////////////////////////////////////
	//		write output file			//
	//////////////////////////////////////////////////
	outfile->cd();
	singleChChk->Write();
	eventChk->Write();
	outfile->Close();
	release_globle();
	printf("program message: merge finished, root file closed\n");

	return 0;
}

void release_globle()
{
	delete wfctaevent;
	delete wfctaslc;
	delete ledsel;
	delete wfctarec;
}








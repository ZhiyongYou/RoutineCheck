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
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TBox.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveStats.h"

#include "WFCTAEvent.h"
#include "LedEvtSel.h"
#include "WFCTASLC.h"
#include "WFCTARec.h"
#include "astro.h"

using namespace std;

int main(int argc, char *argv[])
{
	if(argc<0) {
		printf("Usage: %s outroot wfcta01.txt wfcta02.txt wfcta03.txt wfcta04.txt wfcta05.txt wfcta06.txt\n", argv[0]);
		return 1;
	}
	cout<<"begin time:"<<time(0)<<endl;
	std::string infilename = argv[1];
	size_t posWFCTA = infilename.find("WFCTA");
	if(posWFCTA==std::string::npos)	{	std::cerr << "input file not correct: Format YYMMDD.WFCTA##.check.root.root"<< std::endl;	return 0;}
	std::string outfilename = infilename.substr(posWFCTA-9,22);
	std::string Year = infilename.substr(posWFCTA-9,4);
	std::string Date = infilename.substr(posWFCTA-5,4);
	std::cout << "OutFile: " << outfilename << ".pdf" << std::endl;
	int day = atoi(Date.c_str());
	double drawmjd = lt2mjd(2020, day/100, day%100, 0, 0, 0 )- 2400000.5;
	printf("year:%s date:%s DrawMjd:%lf\n",Year.c_str(), Date.c_str(), drawmjd);

	//set map
	double SIPMMAP[1024][2];
	double interval = 0;//0.36;//1.0; //mm, gaps between subclusters 
	int Interx,Intery;
	double centerx, centery;
	int PIX = 32;
	double D_ConeOut=25.8 ;// mm

	centerx = 414.3;
	centery = 414.3;
	for(int k=0;k<1024;k++){
		int i = k/32;
		int j = k%32;
		Intery = i/4;
		Interx = j/4;
		if(i%2==0)
			SIPMMAP[k][0] = ((j+0.5)*D_ConeOut + interval*Interx-centerx);
		if(i%2==1)
			SIPMMAP[k][0]  = ((j+1)*D_ConeOut + interval*Interx-centerx);
		SIPMMAP[k][1] = ((PIX-i)*D_ConeOut + interval*Intery-centery);
		SIPMMAP[k][0] /= 2870.;
		SIPMMAP[k][1] /=2870.;
		SIPMMAP[k][0] *= 57.3;
		SIPMMAP[k][1] *= 57.3;
	}

	//read file
	char telFile[500];
	strcpy(telFile,"root://eos01.ihep.ac.cn/");
	strcat(telFile,argv[1]);
	TFile* inFile = TFile::Open(telFile,"READ");
	printf("inputing wfcta file -> %s\n",telFile);
	if(!inFile || inFile->IsZombie() || inFile->GetEND()<50)
	{
		printf("%s file error!!\n",telFile);
		return 0;
	}

	Long64_t rabbitTime;
	double rabbittime;
	double mjd;
	int packageCheck,eventNumberCheck;
	int event_type;
	double deltaT;
	int Npix;
	double Size;
	double MeanX, MeanY;
	TTree* eventChk = (TTree *)inFile->Get("eventChk");
	eventChk->SetBranchAddress("rabbitTime",&rabbitTime);
	eventChk->SetBranchAddress("rabbittime",&rabbittime);
	eventChk->SetBranchAddress("mjd",&mjd);
	eventChk->SetBranchAddress("packageCheck",&packageCheck); // 0 is for OK; 1 is for bad packageCheck
	eventChk->SetBranchAddress("eventNumberCheck",&eventNumberCheck);
	eventChk->SetBranchAddress("event_type",&event_type); // 0 is cr, 1 is led, 2 is reflect led, 3 is noise, 4 is laser2, 5 is laser3
	eventChk->SetBranchAddress("deltaT",&deltaT);
	eventChk->SetBranchAddress("Npix",&Npix);
	eventChk->SetBranchAddress("Size",&Size);
	eventChk->SetBranchAddress("MeanY",&MeanY);
	eventChk->SetBranchAddress("MeanX",&MeanX);

	long Nwfctaevts = eventChk->GetEntries();
	printf("Entries: %ld\n\n\n", Nwfctaevts);
	//loop events and get draw range
	double dtmax=2;
	double sizemin=0, sizemax=8;
	for(int ientry=0;ientry<Nwfctaevts;ientry++)
	{
		eventChk->GetEntry(ientry);
		if(Npix<5)	continue;
//		if(mjd>59132.88)	continue;

		sizemin = sizemin<log10(Size) ? sizemin : log10(Size);
		sizemax = sizemax>log10(Size) ? sizemax : log10(Size);
		if(deltaT>0&&deltaT<2)
		{
			dtmax = dtmax>deltaT ? dtmax : deltaT;
		}
	}

	TH1D* h_rate_cr = new TH1D("h_rate_cr", "h_rate_cr", 51840, drawmjd+0.4, drawmjd+1);
	TH1D* h_npix_cr = new TH1D("h_npix_cr", "h_npix_cr", 100, 0, 1024);
	TH1D* h_size_cr = new TH1D("h_size_cr", "h_size_cr", 100, sizemin*0.9, sizemax*1.1);
	TH1D* h_dt_cr = new TH1D("h_dt_cr", "h_dt_cr", 100, 0, dtmax*1.1);

	TH1D* h_rate_led = new TH1D("h_rate_led", "h_rate_led", 51840, drawmjd+0.4, drawmjd+1);
	TH1D* h_npix_led = new TH1D("h_npix_led", "h_npix_led", 100, 1000, 1250);
	TH1D* h_size_led = new TH1D("h_size_led", "h_size_led", 100, sizemin*0.9, sizemax*1.1);
	TH1D* h_dt_led = new TH1D("h_dt_led", "h_dt_led", 100, 0, dtmax*1.1);

	TH1D* h_rate_R_led = new TH1D("h_rate_R_led", "h_rate_R_led", 51840, drawmjd+0.4, drawmjd+1);
	TH1D* h_npix_R_led = new TH1D("h_npix_R_led", "h_npix_R_led", 100, 1000, 1250);
	TH1D* h_size_R_led = new TH1D("h_size_R_led", "h_size_R_led", 100, sizemin*0.9, sizemax*1.1);
	TH1D* h_dt_R_led = new TH1D("h_dt_R_led", "h_dt_R_led", 100, 0, dtmax*1.1);

	TH1D* h_rate_laser2 = new TH1D("h_rate_laser2", "h_rate_laser2", 51840, drawmjd+0.4, drawmjd+1);
	TH1D* h_npix_laser2 = new TH1D("h_npix_laser2", "h_npix_laser2", 100, 0, 1024);
	TH1D* h_size_laser2 = new TH1D("h_size_laser2", "h_size_laser2", 100, sizemin*0.9, sizemax*1.1);
	TH1D* h_dt_laser2 = new TH1D("h_dt_laser2", "h_dt_laser2", 100, 0, dtmax*1.1);

	TH1D* h_rate_laser3 = new TH1D("h_rate_laser3", "h_rate_laser3", 51840, drawmjd+0.4, drawmjd+1);
	TH1D* h_npix_laser3 = new TH1D("h_npix_laser3", "h_npix_laser3", 100, 0, 1024);
	TH1D* h_size_laser3 = new TH1D("h_size_laser3", "h_size_laser3", 100, sizemin*0.9, sizemax*1.1);
	TH1D* h_dt_laser3 = new TH1D("h_dt_laser3", "h_dt_laser3", 100, 0, dtmax*1.1);

	//loop events and get fill hist
	for(int ientry=0;ientry<Nwfctaevts;ientry++)
	{
		eventChk->GetEntry(ientry);
		if(Npix<5)	continue;
//		if(mjd>59132.88)	continue;

		if(0==event_type) {	
			h_rate_cr->Fill(mjd);
			h_npix_cr->Fill(Npix);
			h_size_cr->Fill(log10(Size));
			h_dt_cr->Fill(deltaT);
		}
		if(1==event_type) {	
			h_rate_led->Fill(mjd);
			h_npix_led->Fill(Npix);
			h_size_led->Fill(log10(Size));
			h_dt_led->Fill(deltaT);
		}
		if(2==event_type) {	
			h_rate_R_led->Fill(mjd);
			h_npix_R_led->Fill(Npix);
			h_size_R_led->Fill(log10(Size));
			h_dt_R_led->Fill(deltaT);
		}
		if(4==event_type) {	
			h_rate_laser2->Fill(mjd);
			h_npix_laser2->Fill(Npix);
			h_size_laser2->Fill(log10(Size));
			h_dt_laser2->Fill(deltaT);
		}
		if(5==event_type) {	
			h_rate_laser3->Fill(mjd);
			h_npix_laser3->Fill(Npix);
			h_size_laser3->Fill(log10(Size));
			h_dt_laser3->Fill(deltaT);
		}
	}
	TCanvas* c_eventrate = new TCanvas("c_eventrate", "c_eventrate", 1200, 600);
	c_eventrate->cd();
	gPad->SetLogy();
	h_rate_cr->SetLineColor(2);
	h_rate_cr->SetTitle("Event Rate");
	h_rate_cr->GetXaxis()->SetTitle("Mjd");
	h_rate_cr->GetYaxis()->SetRangeUser(0.1,200);
	h_rate_led->SetLineColor(3);
	h_rate_R_led->SetLineColor(4);
	h_rate_laser2->SetLineColor(5);
	h_rate_laser3->SetLineColor(6);
	h_rate_cr->Draw();
	h_rate_R_led->Draw("same");
	h_rate_led->Draw("same");
	h_rate_laser2->Draw("same");
	h_rate_laser3->Draw("same");
	TLegend *legend_rate = new TLegend(.7,.75,0.9,0.9);
	legend_rate->Clear();
	legend_rate->AddEntry(h_rate_cr,"rate of cosmic ray");
	legend_rate->AddEntry(h_rate_R_led,"rate of reflect led");
	legend_rate->AddEntry(h_rate_led,"rate of led");
	legend_rate->AddEntry(h_rate_laser2,"rate of laser2");
	legend_rate->AddEntry(h_rate_laser3,"rate of laser3");
	legend_rate->Draw("same");
//	gPad->Update();	
//	((TPaveStats*)h_rate_cr->FindObject("stats"))->SetOptStat(0);
//	gStyle->SetOptStat(0);
	c_eventrate->Print(Form("%s/%s.pdf(","/eos/user/y/youzhiyong/RoutineCheck/eventCheck/Result",outfilename.c_str()),"pdf");


	TCanvas* c_event = new TCanvas("c_event", "c_event", 1800, 1800);
	c_event->Divide(3,3);
	c_event->cd(1);
	gPad->SetLogy();
	h_npix_cr->SetTitle("Npix of Cosmic Ray");
	h_npix_cr->GetXaxis()->SetTitle("Npix");
	h_npix_cr->GetXaxis()->CenterTitle();
	h_npix_cr->Draw();
	c_event->cd(2);
	gPad->SetLogy();
	h_npix_led->SetTitle("Npix of LED");
	h_npix_led->GetXaxis()->SetTitle("Npix");
	h_npix_led->GetXaxis()->CenterTitle();
	h_npix_led->Draw();
	c_event->cd(3);
	gPad->SetLogy();
	h_npix_laser2->SetTitle("Npix of Laser2");
	h_npix_laser2->GetXaxis()->SetTitle("Npix");
	h_npix_laser2->GetXaxis()->CenterTitle();
	h_npix_laser2->Draw();
	c_event->cd(4);
	gPad->SetLogy();
	h_size_cr->SetTitle("log10(Size) of Cosmic Ray");
	h_size_cr->GetXaxis()->SetTitle("log10(Size)");
	h_size_cr->GetXaxis()->CenterTitle();
	h_size_cr->Draw();
	c_event->cd(5);
	gPad->SetLogy();
	h_size_led->SetTitle("log10(Size) of LED");
	h_size_led->GetXaxis()->SetTitle("log10(Size)");
	h_size_led->GetXaxis()->CenterTitle();
	h_size_led->Draw();
	c_event->cd(6);
	gPad->SetLogy();
	h_size_laser2->SetTitle("log10(Size) of Laser2");
	h_size_laser2->GetXaxis()->SetTitle("log10(Size)");
	h_size_laser2->GetXaxis()->CenterTitle();
	h_size_laser2->Draw();
	c_event->cd(7);
	gPad->SetLogy();
	h_dt_cr->SetTitle("DeltaT of Cosmic Ray");
	h_dt_cr->GetXaxis()->SetTitle("DeltaT (s)");
	h_dt_cr->GetXaxis()->CenterTitle();
	h_dt_cr->Draw();
	c_event->cd(8);
	gPad->SetLogy();
	h_dt_led->SetTitle("DeltaT of LED");
	h_dt_led->GetXaxis()->SetTitle("DeltaT (s)");
	h_dt_led->GetXaxis()->CenterTitle();
	h_dt_led->Draw();
	c_event->cd(9);
	gPad->SetLogy();
	h_dt_laser2->SetTitle("DeltaT of Laser2");
	h_dt_laser2->GetXaxis()->SetTitle("DeltaT (s)");
	h_dt_laser2->GetXaxis()->CenterTitle();
	h_dt_laser2->Draw();
	c_event->Print(Form("%s/%s.pdf","/eos/user/y/youzhiyong/RoutineCheck/eventCheck/Result",outfilename.c_str()),"pdf");

	/*
	delete h_rate_led;
	delete h_npix_led;
	delete h_size_led;
	delete h_dt_led;
	delete h_rate_cr;
	delete h_npix_cr;
	delete h_size_cr;
	delete h_dt_cr;
	delete h_rate_laser2;
	delete h_npix_laser2;
	delete h_size_laser2;
	delete h_dt_laser2;

	delete c_eventrate;
	delete legend_rate;
	delete c_eventrate;
	*/






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
	int trigger_laser2[1024]={0};

	double pe_cr[1024]={0};
	double pe_led[1024]={0};
	double pe_laser2[1024]={0};

	double h_l_ratio_ave[1024]={0};
	double h_l_ratio_rms[1024]={0};
	int h_l_ratio_n[1024]={0};

	double obs_time_cr=0;
	double obs_time_led=0;
	double obs_time_laser2=0;

	double totaltrigger_cr=0;
	double totaltrigger_led=0;
	double totaltrigger_laser2=0;
	TTree* singleChChk = (TTree *)inFile->Get("singleChChk");
	if(singleChChk==nullptr)
	{
		printf("%s is null file\n",telFile);
		return 0;
	}
	singleChChk->SetBranchAddress("sipm",sipm);

	singleChChk->SetBranchAddress("dark_base_rms",dark_base_rms);
	singleChChk->SetBranchAddress("dark_base_ave",dark_base_ave);
	singleChChk->SetBranchAddress("dark_base_n",dark_base_n);

	singleChChk->SetBranchAddress("base_rms",base_rms);
	singleChChk->SetBranchAddress("base_ave",base_ave);
	singleChChk->SetBranchAddress("base_n",base_n);

	singleChChk->SetBranchAddress("h_l_ratio_ave",h_l_ratio_ave);
	singleChChk->SetBranchAddress("h_l_ratio_rms",h_l_ratio_rms);
	singleChChk->SetBranchAddress("h_l_ratio_n",h_l_ratio_n);

	singleChChk->SetBranchAddress("trigger_led",trigger_led);
	singleChChk->SetBranchAddress("trigger_cr",trigger_cr);
	singleChChk->SetBranchAddress("trigger_laser2",trigger_laser2);

	singleChChk->SetBranchAddress("pe_led",pe_led);
	singleChChk->SetBranchAddress("pe_cr",pe_cr);
	singleChChk->SetBranchAddress("pe_laser2",pe_laser2);

	singleChChk->SetBranchAddress("obs_time_cr",&obs_time_cr);
	singleChChk->SetBranchAddress("obs_time_led",&obs_time_led);
	singleChChk->SetBranchAddress("obs_time_laser2",&obs_time_laser2);

	singleChChk->SetBranchAddress("totaltrigger_cr",&totaltrigger_cr);
	singleChChk->SetBranchAddress("totaltrigger_led",&totaltrigger_led);
	singleChChk->SetBranchAddress("totaltrigger_laser2",&totaltrigger_laser2);

	singleChChk->GetEntry(0);

	TH2D *h2[30];
	for(int i=0; i<30; i++){
		h2[i]= new TH2D(Form("h2%d",i),"dfr",10,-9,9,10,-9,9);
	}

	double occupency_cr[1024]={0};
	double occupency_led[1024]={0};
	double occupency_laser2[1024]={0};
	for(int i=0;i<1024;i++){
		if(obs_time_cr>0)		occupency_cr[i] = trigger_cr[i] / totaltrigger_cr;
		if(obs_time_led>0)		occupency_led[i] = trigger_led[i] / totaltrigger_led;
		if(obs_time_laser2>0)	occupency_laser2[i] = trigger_laser2[i] / totaltrigger_laser2;
	}
	double maxtrigger=0;
	double mintrigger =1000000;
	for(int i=0;i<1024;i++){
		if(occupency_cr[i]>maxtrigger) maxtrigger = occupency_cr[i]*1.2;
		if(occupency_cr[i]<mintrigger) mintrigger = occupency_cr[i]*0.8;
	}
	TCanvas *c5 = new TCanvas("c5","single trigger rate of cosmic rays",1500,1000);
	c5->Divide(3,2);
	c5->cd(1);
	h2[0]->Draw();
	h2[0]->SetTitle("Occupancy of CRs;X(#circ);Y(#circ)");
	TBox *bx[1024];
	TH1D *h_occupency_cr = new TH1D("h_occupency_cr","Occupancy of CRs",50, mintrigger, maxtrigger);
	for(int i=0; i<1024; i++){
		bx[i] = new TBox(SIPMMAP[i][0]-0.25,SIPMMAP[i][1]-0.25,SIPMMAP[i][0]+0.25,SIPMMAP[i][1]+0.25);
		bx[i]->SetFillColor(int((occupency_cr[i]-mintrigger)/(maxtrigger-mintrigger)*50)+49);
		bx[i]->Draw();
		h_occupency_cr->Fill(occupency_cr[i]);
	}

	TBox *bxcolor[50];
	TGaxis* axis1 = new TGaxis(9,-9,9,9,mintrigger,maxtrigger,20,"+LS");
	axis1->SetTickSize(0.02);
	for(int i=0;i<50;i++){
		bxcolor[i] = new TBox(9,18/50.*i+(-9),10,18./50*(i+1)+(-9));
		bxcolor[i]->SetFillColor(i+49);
		bxcolor[i]->Draw();
	}
	axis1->Draw();

	c5->cd(4);
	h_occupency_cr->SetTitle("Occupancy of CRs;Occupancy;");
	h_occupency_cr->Draw();

	maxtrigger=0;
	mintrigger =1000000;
	for(int i=0;i<1024;i++){
		if(occupency_led[i]>maxtrigger) maxtrigger = occupency_led[i]*1.2;
		if(occupency_led[i]<mintrigger) mintrigger = occupency_led[i]*0.8;
	}
	TH1D *h_occupency_led = new TH1D("h_occupency_led","Occupancy of LEDs",50, mintrigger, maxtrigger);
	c5->cd(2);
	h2[1]->Draw();
	h2[1]->SetTitle("Occupancy of LEDs;X(#circ);Y(#circ)");

	for(int i=0; i<1024; i++){
		bx[i] = new TBox(SIPMMAP[i][0]-0.25,SIPMMAP[i][1]-0.25,SIPMMAP[i][0]+0.25,SIPMMAP[i][1]+0.25);
		bx[i]->SetFillColor(int((occupency_led[i]-mintrigger)/(maxtrigger-mintrigger)*50)+49);
		bx[i]->Draw();
		h_occupency_led->Fill(occupency_led[i]);
	}

	TGaxis* axis2 = new TGaxis(9,-9,9,9,mintrigger,maxtrigger,20,"+LS");
	axis2->SetTickSize(0.02);
	for(int i=0;i<50;i++){
		bxcolor[i] = new TBox(9,18/50.*i+(-9),10,18./50*(i+1)+(-9));
		bxcolor[i]->SetFillColor(i+49);
		bxcolor[i]->Draw();
	}
	axis2->Draw();

	c5->cd(5);
	h_occupency_led->SetTitle("Occupancy of LEDs;Occupancy;");
	h_occupency_led->Draw();

	maxtrigger=0;
	mintrigger =1000000;
	for(int i=0;i<1024;i++){
		if(occupency_laser2[i]>maxtrigger) maxtrigger = occupency_laser2[i]*1.2;
		if(occupency_laser2[i]<mintrigger) mintrigger = occupency_laser2[i]*0.8;
	}
	TH1D *h_occupency_laser2 = new TH1D("h_occupency_laser2","Occupancy of Laser2",50, mintrigger, maxtrigger);
	c5->cd(3);
	h2[2]->Draw();
	h2[2]->SetTitle("Occupancy Laser2;X(#circ);Y(#circ)");

	for(int i=0; i<1024; i++){
		bx[i] = new TBox(SIPMMAP[i][0]-0.25,SIPMMAP[i][1]-0.25,SIPMMAP[i][0]+0.25,SIPMMAP[i][1]+0.25);
		bx[i]->SetFillColor(int((occupency_laser2[i]-mintrigger)/(maxtrigger-mintrigger)*50)+49);
		bx[i]->Draw();
		h_occupency_laser2->Fill(occupency_laser2[i]);
	}

	TGaxis* axis3 = new TGaxis(9,-9,9,9,mintrigger,maxtrigger,20,"+LS");
	axis3->SetTickSize(0.02);
	for(int i=0;i<50;i++){
		bxcolor[i] = new TBox(9,18/50.*i+(-9),10,18./50*(i+1)+(-9));
		bxcolor[i]->SetFillColor(i+49);
		bxcolor[i]->Draw();
	}
	axis3->Draw();

	c5->cd(6);
	h_occupency_laser2->SetTitle("Occupancy of Laser2;Occupancy;");
	h_occupency_laser2->Draw();

	gPad->Update();	
	((TPaveStats*)h2[0]->FindObject("stats"))->SetOptStat(0);
	((TPaveStats*)h2[1]->FindObject("stats"))->SetOptStat(0);
	((TPaveStats*)h2[2]->FindObject("stats"))->SetOptStat(0);
	c5->Print(Form("%s/%s.pdf","/eos/user/y/youzhiyong/RoutineCheck/eventCheck/Result",outfilename.c_str()),"pdf");


	/*
	double trigger_cr_[1024]={0};
	double trigger_led_[1024]={0};
	double trigger_laser2_[1024]={0};
	for(int i=0;i<1024;i++){
		if(obs_time_cr>0)		trigger_cr_[i] = trigger_cr[i] / obs_time_cr;
		if(obs_time_led>0)		trigger_led_[i] = trigger_led[i] / obs_time_led;
		if(obs_time_laser2>0)	trigger_laser2_[i] = trigger_laser2[i] / obs_time_laser2;
	}
	maxtrigger=0;
	mintrigger =1000000;
	for(int i=0;i<1024;i++){
		if(trigger_cr_[i]>maxtrigger) maxtrigger = trigger_cr_[i]*1.2;
		if(trigger_cr_[i]<mintrigger) mintrigger = trigger_cr_[i]*0.8;
	}
	TCanvas *c50 = new TCanvas("c50","single trigger rate of cosmic rays",1500,1000);
	c50->Divide(3,2);
	c50->cd(1);
	h2[0]->Draw();
	h2[0]->SetTitle("single trigger rate of CRs;X(#circ);Y(#circ)");
	TH1D *h_trigger_cr_ = new TH1D("h_trigger_cr_","single trigger rate of CRs",50, mintrigger, maxtrigger);
	for(int i=0; i<1024; i++){
		bx[i] = new TBox(SIPMMAP[i][0]-0.25,SIPMMAP[i][1]-0.25,SIPMMAP[i][0]+0.25,SIPMMAP[i][1]+0.25);
		bx[i]->SetFillColor(int((trigger_cr_[i]-mintrigger)/(maxtrigger-mintrigger)*50)+49);
		bx[i]->Draw();
		h_trigger_cr_->Fill(trigger_cr_[i]);
	}

	TGaxis* axis10 = new TGaxis(9,-9,9,9,mintrigger,maxtrigger,20,"+LS");
	axis1->SetTickSize(0.02);
	for(int i=0;i<50;i++){
		bxcolor[i] = new TBox(9,18/50.*i+(-9),10,18./50*(i+1)+(-9));
		bxcolor[i]->SetFillColor(i+49);
		bxcolor[i]->Draw();
	}
	axis10->Draw();

	c50->cd(4);
	h_trigger_cr_->SetTitle("single trigger rate of CRs;single trigger rate;");
	h_trigger_cr_->Draw();

	maxtrigger=0;
	mintrigger =1000000;
	for(int i=0;i<1024;i++){
		if(trigger_led_[i]>maxtrigger) maxtrigger = trigger_led_[i]*1.2;
		if(trigger_led_[i]<mintrigger) mintrigger = trigger_led_[i]*0.8;
	}
	TH1D *h_trigger_led_ = new TH1D("h_trigger_led_","single trigger rate of LEDs",50, mintrigger, maxtrigger);
	c50->cd(2);
	h2[1]->Draw();
	h2[1]->SetTitle("single trigger rate of LEDs;X(#circ);Y(#circ)");

	for(int i=0; i<1024; i++){
		bx[i] = new TBox(SIPMMAP[i][0]-0.25,SIPMMAP[i][1]-0.25,SIPMMAP[i][0]+0.25,SIPMMAP[i][1]+0.25);
		bx[i]->SetFillColor(int((trigger_led_[i]-mintrigger)/(maxtrigger-mintrigger)*50)+49);
		bx[i]->Draw();
		h_trigger_led_->Fill(trigger_led_[i]);
	}

	TGaxis* axis20 = new TGaxis(9,-9,9,9,mintrigger,maxtrigger,20,"+LS");
	axis2->SetTickSize(0.02);
	for(int i=0;i<50;i++){
		bxcolor[i] = new TBox(9,18/50.*i+(-9),10,18./50*(i+1)+(-9));
		bxcolor[i]->SetFillColor(i+49);
		bxcolor[i]->Draw();
	}
	axis20->Draw();

	c50->cd(5);
	h_trigger_led_->SetTitle("single trigger rate of LEDs;single trigger rate;");
	h_trigger_led_->Draw();

	maxtrigger=0;
	mintrigger =1000000;
	for(int i=0;i<1024;i++){
		if(trigger_laser2_[i]>maxtrigger) maxtrigger = trigger_laser2_[i]*1.2;
		if(trigger_laser2_[i]<mintrigger) mintrigger = trigger_laser2_[i]*0.8;
	}
	TH1D *h_trigger_laser2_ = new TH1D("h_trigger_laser2_","single trigger rate of Laser2",50, mintrigger, maxtrigger);
	c50->cd(3);
	h2[2]->Draw();
	h2[2]->SetTitle("single trigger rate Laser2;X(#circ);Y(#circ)");

	for(int i=0; i<1024; i++){
		bx[i] = new TBox(SIPMMAP[i][0]-0.25,SIPMMAP[i][1]-0.25,SIPMMAP[i][0]+0.25,SIPMMAP[i][1]+0.25);
		bx[i]->SetFillColor(int((trigger_laser2_[i]-mintrigger)/(maxtrigger-mintrigger)*50)+49);
		bx[i]->Draw();
		h_trigger_laser2_->Fill(trigger_laser2_[i]);
	}

	TGaxis* axis30 = new TGaxis(9,-9,9,9,mintrigger,maxtrigger,20,"+LS");
	axis3->SetTickSize(0.02);
	for(int i=0;i<50;i++){
		bxcolor[i] = new TBox(9,18/50.*i+(-9),10,18./50*(i+1)+(-9));
		bxcolor[i]->SetFillColor(i+49);
		bxcolor[i]->Draw();
	}
	axis30->Draw();

	c50->cd(6);
	h_trigger_laser2_->SetTitle("single trigger rate of Laser2;single trigger rate;");
	h_trigger_laser2_->Draw();

	gPad->Update();	
	((TPaveStats*)h2[0]->FindObject("stats"))->SetOptStat(0);
	((TPaveStats*)h2[1]->FindObject("stats"))->SetOptStat(0);
	((TPaveStats*)h2[2]->FindObject("stats"))->SetOptStat(0);
	c50->Print(Form("%s/%s.pdf","/eos/user/y/youzhiyong/RoutineCheck/eventCheck/Result",outfilename.c_str()),"pdf");
	*/



	maxtrigger=0;
	mintrigger =1000000;
	for(int i=0;i<1024;i++){
		if(pe_cr[i]>maxtrigger) maxtrigger = pe_cr[i]*1.2;
		if(pe_cr[i]<mintrigger) mintrigger = pe_cr[i]*0.8;
	}

	TH1D *h_pe_cr = new TH1D("h_pe_cr","single average pe of crs",50, mintrigger, maxtrigger);
	TCanvas *c6 = new TCanvas("c6","single pe of crs",1500,1000);
	c6->Divide(3,2);
	c6->cd(1);
	h2[0]->Draw();
	h2[0]->SetTitle("single average pe CRs;X(#circ);Y(#circ)");

	for(int i=0; i<1024; i++){
		bx[i] = new TBox(SIPMMAP[i][0]-0.25,SIPMMAP[i][1]-0.25,SIPMMAP[i][0]+0.25,SIPMMAP[i][1]+0.25);
		bx[i]->SetFillColor(int((pe_cr[i]-mintrigger)/(maxtrigger-mintrigger)*50)+49);
		bx[i]->Draw();
		h_pe_cr->Fill(pe_cr[i]);
	}
	TGaxis* axis4 = new TGaxis(9,-9,9,9,mintrigger,maxtrigger,20,"+LS");
	axis4->SetTickSize(0.02);
	for(int i=0;i<50;i++){
		bxcolor[i] = new TBox(9,18/50.*i+(-9),10,18./50*(i+1)+(-9));
		bxcolor[i]->SetFillColor(i+49);
		bxcolor[i]->Draw();
	}
	axis4->Draw();
	c6->cd(4);
	h_pe_cr->SetTitle("single pe of CRs;single_pe_cr;");
	h_pe_cr->Draw();

	maxtrigger=0;
	mintrigger =1000000;
	for(int i=0;i<1024;i++){
		if(pe_led[i]>maxtrigger) maxtrigger = pe_led[i]*1.2;
		if(pe_led[i]<mintrigger) mintrigger = pe_led[i]*0.8;
	}
	TH1D *h_pe_led = new TH1D("h_pe_led","single average pe of LEDs",50, mintrigger, maxtrigger);
	c6->cd(2);
	h2[1]->Draw();
	h2[1]->SetTitle("single average pe of LEDs;X(#circ);Y(#circ)");

	for(int i=0; i<1024; i++){
		bx[i] = new TBox(SIPMMAP[i][0]-0.25,SIPMMAP[i][1]-0.25,SIPMMAP[i][0]+0.25,SIPMMAP[i][1]+0.25);
		bx[i]->SetFillColor(int((pe_led[i]-mintrigger)/(maxtrigger-mintrigger)*50)+49);
		bx[i]->Draw();
		h_pe_led->Fill(pe_led[i]);
	}
	TGaxis* axis5 = new TGaxis(9,-9,9,9,mintrigger,maxtrigger,20,"+LS");
	axis5->SetTickSize(0.02);
	for(int i=0;i<50;i++){
		bxcolor[i] = new TBox(9,18/50.*i+(-9),10,18./50*(i+1)+(-9));
		bxcolor[i]->SetFillColor(i+49);
		bxcolor[i]->Draw();
	}
	axis5->Draw();
	c6->cd(5);
	h_pe_led->SetTitle("single pe of LEDs;single_pe_led;");
	h_pe_led->Draw();

	double MeanLedPe = h_pe_led->GetMean();
	int Nbad=0;
	for(int i=0; i<1024; i++){
		if(fabs(pe_led[i]-MeanLedPe)>MeanLedPe*0.1) Nbad++;
	}

	/*
	fprintf(fp,"#SIPMs with Pe of LED events exceeded +- 0.1 %d\n",Nbad);
	for(int i=0; i<1024; i++){
		if(fabs(pe_led[i]-MeanLedPe)>MeanLedPe*0.1) fprintf(fp,"InconsistentPeLed %d %f\n",i,pe_led[i]);
	}
	*/

	maxtrigger=0;
	mintrigger =1000000;
	for(int i=0;i<1024;i++){
		if(pe_laser2[i]>maxtrigger) maxtrigger = pe_laser2[i]*1.2;
		if(pe_laser2[i]<mintrigger) mintrigger = pe_laser2[i]*0.8;
	}

	TH1D *h_pe_laser2 = new TH1D("h_pe_laser2","single pe of Laser2",50, mintrigger, maxtrigger);
	c6->cd(3);
	h2[2]->Draw();
	h2[2]->SetTitle("single pe Laser2;X(#circ);Y(#circ)");

	for(int i=0; i<1024; i++){
		bx[i] = new TBox(SIPMMAP[i][0]-0.25,SIPMMAP[i][1]-0.25,SIPMMAP[i][0]+0.25,SIPMMAP[i][1]+0.25);
		bx[i]->SetFillColor(int((pe_laser2[i]-mintrigger)/(maxtrigger-mintrigger)*50)+49);
		bx[i]->Draw();
		h_pe_laser2->Fill(pe_laser2[i]);
	}
	TGaxis* axis6 = new TGaxis(9,-9,9,9,mintrigger,maxtrigger,20,"+LS");
	axis6->SetTickSize(0.02);
	for(int i=0;i<50;i++){
		bxcolor[i] = new TBox(9,18/50.*i+(-9),10,18./50*(i+1)+(-9));
		bxcolor[i]->SetFillColor(i+49);
		bxcolor[i]->Draw();
	}
	axis6->Draw();
	c6->cd(6);
	h_pe_laser2->SetTitle("single pe of Laser2;single_pe_laser2;");
	h_pe_laser2->Draw();

	gPad->Update();	
	((TPaveStats*)h2[0]->FindObject("stats"))->SetOptStat(0);
	((TPaveStats*)h2[1]->FindObject("stats"))->SetOptStat(0);
	((TPaveStats*)h2[2]->FindObject("stats"))->SetOptStat(0);
	c6->Print(Form("%s/%s.pdf)","/eos/user/y/youzhiyong/RoutineCheck/eventCheck/Result",outfilename.c_str()),"pdf");



















	cout<<"deal time:"<<time(0)<<endl<<endl;

	printf("program message: merge finished, root file closed\n");
	inFile->Close();

	return 0;
}








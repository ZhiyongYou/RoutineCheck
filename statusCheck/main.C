#include <stdlib.h>
#include <stdio.h>
#include <TPad.h>
#include <TStyle.h>
#include <TH1D.h>
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TGaxis.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "camera.h"
#include "ReadFile.h"

const double MJD19700101 = 40587;
const double TAI2UTC = 37;
using namespace std;
struct Time{
	int year;
	int month;
	int day;
	int hour;
	int minite;
	int second;
};

double PMTMAP[2][1024];
void GetPMTMAP()
{
	int Pix = 32;
	double D_ConeOut = 25.8;
	double interval = 0.0;
	int Interx,Intery;

	for(int k=0;k<1024;k++){
		int i = k/32;
		int j = k%32;
		Intery = i/4;
		Interx = j/4;
		if(i%2==0)
			PMTMAP[0][k] = (j+0.5)*D_ConeOut + interval*Interx;
		if(i%2==1)
			PMTMAP[0][k] = (j+1)*D_ConeOut + interval*Interx;
		PMTMAP[1][k] = (Pix-i)*D_ConeOut + interval*Intery;

		PMTMAP[0][k] = (PMTMAP[0][k]-414.3) *TMath::RadToDeg()/ 2870;
		PMTMAP[1][k] = (PMTMAP[1][k]-414.3) *TMath::RadToDeg()/ 2870;
	}
}

void rabbittime2lt(Long64_t rabbitTime, double rabbittime, Time &time);
void graphSet(TGraph *graph, string name, string title, const Time &time);
void threshSet(TGraph *graph, string name, int pointcolor, int markerstyle, double markersize, const Time &time);
void invert(double *conf);

int main(int argc, char**argv)
{
	if(argc!=7){
		printf("use %s inrootfilelist outrootfile outpdffiledir basefle hotpix telnum\n",argv[0]);
		return 0;
	}

	GetPMTMAP();
	int tel = atoi(argv[6]);
	printf("telescope:%d\n",tel);

	//read status files
	ReadFile confile;
	char filename[100];

	double HVSetting[1024];
	double HVReturn[1024]={0};
	double HVReturnCount[1024]={0};
//	sprintf(filename,"/workfs/ybj/wfcta/DataCheck/RoutineCheck-root6-FileByFile/HVsettings/hv0%d.txt",tel);
	sprintf(filename,"/workfs/ybj/youzhiyong/RoutineCheck/statusCheck/config/hvsetting/hv0%d.txt",tel);
	confile.OpenFile(filename);
	confile.ReadInfo(HVSetting,"hv");
	confile.CloseFile();
	if(tel==5){	invert(HVSetting);}
	//for(int i=0;i<1024;i++)	{	printf("hv:%d %lf\n",i,HVSetting[i]);}

	double BASESetting[1024];
	sprintf(filename,"%s",argv[4]);
	confile.OpenFile(filename);
	confile.ReadInfo(BASESetting,"base");
	confile.CloseFile();
	if(tel==5){	invert(BASESetting);}
	//for(int i=0;i<1024;i++){	printf("base:%d %lf\n",i,BASESetting[i]);}
	
	vector<int> noisSiPM;
	vector<int> out_noisSiPM;
	sprintf(filename,"%s",argv[5]);
	//GetHotPix(tel,noisSiPM,out_noisSiPM);
	confile.OpenFile(filename);
	confile.GetNoisePix(tel,noisSiPM,out_noisSiPM);
	confile.CloseFile();
	for(int ii=0;ii<noisSiPM.size();ii++){	printf("noisSiPM:%d\n",noisSiPM.at(ii));}
	for(int ii=0;ii<out_noisSiPM.size();ii++){	printf("out_noisSiPM:%d\n",out_noisSiPM.at(ii));}
	
	char statusfile[500];
	FILE *statusfp;
	TChain *Statuschain = new TChain("Status");
	statusfp = fopen(argv[1],"r");
	while(!feof(statusfp)){
		fscanf(statusfp,"%s\n",statusfile);
		printf("%s\n",statusfile);
		Statuschain->Add(statusfile);
	}
	fclose(statusfp);

	short iTel;
	int fpgaVersion[10];
	int f9mode;
	int f9pattern;
	int DbVersion[2][89];
	int ClbVersion[2][89];
	int fired_tube;
	Long64_t status_readback_Time;
	double status_readback_time;
	int mask[1024];
	short single_thresh[1024];
	short record_thresh[1024];
	float HV[1024];
	float PreTemp[1024];
	float ClbTemp[1024];
	Statuschain->SetBranchAddress("iTel",&iTel);
	Statuschain->SetBranchAddress("fpgaVersion",fpgaVersion);
	Statuschain->SetBranchAddress("f9mode",&f9mode);
	Statuschain->SetBranchAddress("f9pattern",&f9pattern);
	Statuschain->SetBranchAddress("DbVersion",DbVersion);
	Statuschain->SetBranchAddress("ClbVersion",ClbVersion);
	Statuschain->SetBranchAddress("fired_tube",&fired_tube);
	Statuschain->SetBranchAddress("status_readback_Time",&status_readback_Time);
	Statuschain->SetBranchAddress("status_readback_time",&status_readback_time);
	Statuschain->SetBranchAddress("mask",mask);
	Statuschain->SetBranchAddress("single_thresh",single_thresh);
	Statuschain->SetBranchAddress("record_thresh",record_thresh);
	Statuschain->SetBranchAddress("HV",HV);
	Statuschain->SetBranchAddress("PreTemp",PreTemp);
	Statuschain->SetBranchAddress("ClbTemp",ClbTemp);
	cout<<"status:"<<Statuschain->GetEntries()<<endl<<endl;

	const int Nentry = Statuschain->GetEntries();

	int sipm[1024];  for(int i=0;i<1024;i++) {sipm[i]=i;}
	double temp_ave[1024]={0};
	double temp_rms[1024]={0};
	int temp_n[1024]={0};
	double clb_temp_ave[1024]={0};
	double clb_temp_rms[1024]={0};
	int clb_temp_n[1024]={0};
	double dHV_ave[1024]={0};
	double dHV_rms[1024]={0};
	int dHV_n[1024]={0};
	TFile* fout=TFile::Open(argv[2],"RECREATE");
	TTree *hv_temp = new TTree("hv_temp","info of hv_temp");
	hv_temp->Branch("sipm",sipm,"sipm[1024]/I");
	hv_temp->Branch("temp_rms",temp_rms,"temp_rms[1024]/D");
	hv_temp->Branch("temp_ave",temp_ave,"temp_ave[1024]/D");
	hv_temp->Branch("temp_n",temp_n,"temp_n[1024]/I");
	hv_temp->Branch("clb_temp_rms",clb_temp_rms,"clb_temp_rms[1024]/D");
	hv_temp->Branch("clb_temp_ave",clb_temp_ave,"clb_temp_ave[1024]/D");
	hv_temp->Branch("clb_temp_n",clb_temp_n,"clb_temp_n[1024]/I");
	hv_temp->Branch("dHV_rms",dHV_rms,"dHV_rms[1024]/D");
	hv_temp->Branch("dHV_ave",dHV_ave,"dHV_ave[1024]/D");  
	hv_temp->Branch("dHV_n",dHV_n,"dHV_n[1024]/I");

	int singlethresh[1024]={0};
	int recordthresh[1024]={0};
	int center_thresh, center_thresh_num;
	int out_thresh, out_thresh_num;
	int cen_n_thresh, cen_n_thresh_num;
	int out_n_thresh, out_n_thresh_num;
	int Record_Thresh, Record_Thresh_num;
	double tel_sipm_T, tel_clb_T, tel_HV;
	int tel_sipm_T_num, tel_clb_T_num, tel_HV_num;
	TTree *status = new TTree("status","info of status");
	status->Branch("sipm",sipm,"sipm[1024]/I");
	status->Branch("mask",mask,"mask[1024]/I");  
	status->Branch("singlethresh",singlethresh,"singlethresh[1024]/I");
	status->Branch("recordthresh",recordthresh,"recordthresh[1024]/I");
	status->Branch("center_thresh",&center_thresh,"center_thresh/I");
	status->Branch("out_thresh",&out_thresh,"out_thresh/I");
	status->Branch("cen_n_thresh",&cen_n_thresh,"cen_n_thresh/I");
	status->Branch("out_n_thresh",&out_n_thresh,"out_n_thresh/I");
	status->Branch("Record_Thresh",&Record_Thresh,"Record_Thresh/I");
	status->Branch("tel_HV",&tel_HV,"tel_HV/D");
	status->Branch("tel_sipm_T",&tel_sipm_T,"tel_sipm_T/D");
	status->Branch("tel_clb_T",&tel_clb_T,"tel_clb_T/D");

	int SC_X[64],SC_Y[64];
	int SC_F[64],SC_DB[64];
	int db_ver[64],clb_ver[64];
	int f_ver[10];
	int MASK[1024];
	TTree *fpga_ver = new TTree("fpga_ver","info of fpga_ver");
	fpga_ver->Branch("SC_F",SC_F,"SC_F[64]/I");
	fpga_ver->Branch("SC_DB",SC_DB,"SC_DB[64]/I");
	fpga_ver->Branch("SC_X",SC_X,"SC_X[64]/I");
	fpga_ver->Branch("SC_Y",SC_Y,"SC_Y[64]/I");
	fpga_ver->Branch("db_ver",db_ver,"db_ver[64]/I");
	fpga_ver->Branch("clb_ver",clb_ver,"clb_ver[64]/I");
	fpga_ver->Branch("f_ver",f_ver,"f_ver[10]/I");
	fpga_ver->Branch("MASK",MASK,"MASK[1024]/I");

	int first=1;
	Time time,time0;
	int ipot=0;

	TGraph *g_work_mode = new TGraph();
	TGraph *g_triggermode = new TGraph();
	TGraph *g_firedtube = new TGraph();
	TGraph *g_tel_HV = new TGraph();
	TGraph *g_tel_sipm_T = new TGraph();
	TGraph *g_tel_clb_T = new TGraph();
	TGraph *g_cen_th = new TGraph();
	TGraph *g_cen_n_th = new TGraph();
	TGraph *g_out_th = new TGraph();
	TGraph *g_out_n_th = new TGraph();
	TGraph *g_rec_th = new TGraph();
	//deal status info
	for(int statusentry=0;statusentry<Statuschain->GetEntries();statusentry++){
		Statuschain->GetEntry(statusentry);
		if(iTel!=tel){continue;}
		if(fired_tube==-1000){continue;}

		rabbittime2lt(status_readback_Time, status_readback_time, time);
		if(first==1)
		{
			time0=time;first=0;

			for(int i=0;i<10;i++){
				f_ver[i] = fpgaVersion[i];
			}
			int isc=0;
			for(int ff=1;ff<9;ff++){
				for(int dd=1;dd<9;dd++){
					int SC = dd*10+ff;
					SC_X[isc] = SC_ADRESS[dd-1][ff-1]%10;
					SC_Y[isc] = SC_ADRESS[dd-1][ff-1]/10;
					SC_F[isc] = ff;
					SC_DB[isc] = dd;
					db_ver[isc] = DbVersion[1][SC]; 
					clb_ver[isc] = ClbVersion[1][SC]; 
					isc++;
				}
			}
			for(int i=0;i<1024;i++){
				MASK[i] = mask[i];
			}
			fpga_ver->Fill();
		}

		int t = (time.day-1)*86400+time.hour*60*60+time.minite*60;

		for(int i=0;i<1024;i++){
			if(PreTemp[i]==-1000){continue;}
			temp_ave[i] += PreTemp[i];
			temp_rms[i] += pow(PreTemp[i],2);
			clb_temp_ave[i] += ClbTemp[i];
			clb_temp_rms[i] += pow(ClbTemp[i],2);
			dHV_ave[i] += (HV[i]-HVSetting[i]);
			dHV_rms[i] += pow(HV[i]-HVSetting[i],2);
			temp_n[i]++;
			clb_temp_n[i]++;
			dHV_n[i]++;
		}
		for(int i=0;i<1024;i++) {
			if(HV[i]==-1000){continue;}
			HVReturn[i] += HV[i];
			HVReturnCount[i]++;
		}

		tel_HV=0; tel_HV_num=0;
		tel_sipm_T=0; tel_sipm_T_num=0;
		tel_clb_T=0; tel_clb_T_num=0;
		for(int i=0;i<1024;i++) {
			if(HV[i]==-1000){continue;}
			tel_HV += HV[i]; tel_HV_num++;
		}
		for(int i=0;i<1024;i++) {
			if(PreTemp[i]==-1000){continue;}
			tel_sipm_T += PreTemp[i]; tel_sipm_T_num++;
		}
		for(int i=0;i<1024;i++) {
			if(ClbTemp[i]==-1000){continue;}
			tel_clb_T += ClbTemp[i]; tel_clb_T_num++;
		}
		if(tel_HV_num!=0) tel_HV /= tel_HV_num;
		if(tel_sipm_T_num!=0) tel_sipm_T /= tel_sipm_T_num;
		if(tel_clb_T_num!=0) tel_clb_T /= tel_clb_T_num;
		g_tel_HV->SetPoint(ipot,t,tel_HV);
		g_tel_sipm_T->SetPoint(ipot,t,tel_sipm_T);
		g_tel_clb_T->SetPoint(ipot,t,tel_clb_T);


		center_thresh=0; center_thresh_num=0;
		out_thresh=0; out_thresh_num=0;
		cen_n_thresh=0; cen_n_thresh_num=0;
		out_n_thresh=0; out_n_thresh_num=0;
		Record_Thresh=0; Record_Thresh_num=0;
		for(int i=0;i<1024;i++){
			if(single_thresh[i]==-1000){continue;}
			singlethresh[i] = round(single_thresh[i]+BASESetting[i]/2.);
			recordthresh[i] = round(record_thresh[i]);

			int out_noisepix=0, noisepix=0;
			for(int ii=0;ii<out_noisSiPM.size();ii++){
				if(i==out_noisSiPM.at(ii)){	out_noisepix=1;break;}
			}
			for(int ii=0;ii<noisSiPM.size();ii++){
				if(i==noisSiPM.at(ii)){	noisepix=1;break;}
			}

			if(out_noisepix){
				out_n_thresh += round(single_thresh[i]+BASESetting[i]/2.);
				out_n_thresh_num++;}
			else if(noisepix){
				cen_n_thresh += round(single_thresh[i]+BASESetting[i]/2.);
				cen_n_thresh_num++;}
			else if(i%32==0||i%32==1||i%32==30||i%32==31){
				out_thresh += round(single_thresh[i]+BASESetting[i]/2.);
				out_thresh_num++;
			}
			else{
				center_thresh += round(single_thresh[i]+BASESetting[i]/2.);
				center_thresh_num++;
			}
			/*
			if(tel==4&&(i==96||i==353||i==574||i==798)){
				cen_n_thresh += round(single_thresh[i]+BASESetting[i]/2.); cen_n_thresh_num++;
			}
			else if(tel==2&&i==735){
				cen_n_thresh += round(single_thresh[i]+BASESetting[i]/2.); cen_n_thresh_num++;
			}
			else if(i%32==0||i%32==1||i%32==30||i%32==31){
				int noisepix=0;
				for(int ii=0;ii<out_noisSiPM.size();ii++){
					if(i==out_noisSiPM.at(ii)){	noisepix=1;break;}
				}
				if(noisepix){	out_n_thresh += round(single_thresh[i]+BASESetting[i]/2.); out_n_thresh_num++;}
				else{	out_thresh += round(single_thresh[i]+BASESetting[i]/2.); 
						//printf("%dout_thresh:%d %.0lf\n",i,out_thresh,round(single_thresh[i]+BASESetting[i]/2.));
						out_thresh_num++;}
			}
			else{
				int noisepix=0;
				for(int ii=0;ii<noisSiPM.size();ii++){
					if(i==noisSiPM.at(ii)){	noisepix=1;break;}
				}
				if(noisepix){	cen_n_thresh += round(single_thresh[i]+BASESetting[i]/2.); cen_n_thresh_num++;}
				else{	center_thresh += round(single_thresh[i]+BASESetting[i]/2.); center_thresh_num++;}
			}
			*/
			Record_Thresh += round(record_thresh[i]); Record_Thresh_num++;
		}
		if(center_thresh_num!=0)	center_thresh = center_thresh/center_thresh_num;
		if(out_thresh_num!=0)	out_thresh = out_thresh/out_thresh_num;
		if(cen_n_thresh_num!=0)	cen_n_thresh = cen_n_thresh/cen_n_thresh_num;
		if(out_n_thresh_num!=0)	out_n_thresh = out_n_thresh/out_n_thresh_num;
		if(Record_Thresh_num!=0)	Record_Thresh = Record_Thresh/Record_Thresh_num;
		g_cen_th->SetPoint(ipot,t,center_thresh);
		g_cen_n_th->SetPoint(ipot,t,cen_n_thresh);
		g_out_th->SetPoint(ipot,t,out_thresh);
		g_out_n_th->SetPoint(ipot,t,out_n_thresh);
		g_rec_th->SetPoint(ipot,t,Record_Thresh);
		status->Fill();
		
		g_work_mode->SetPoint(ipot,t,f9mode);
		g_triggermode->SetPoint(ipot,t,f9pattern);
		g_firedtube->SetPoint(ipot,t,fired_tube);
		ipot++;
	}

	for(int i=0;i<1024;i++) {
		if(HVReturnCount[i]!=0)
			HVReturn[i] /= HVReturnCount[i];
	}

	double max_temp_ave=-1000,max_temp_rms=-1000;
	double min_temp_ave=1000,min_temp_rms=1000;
	double max_clb_temp_ave=-1000,max_clb_temp_rms=-1000;
	double min_clb_temp_ave=1000,min_clb_temp_rms=1000;
	double max_dHV_ave=-1000,max_dHV_rms=-1000;
	double min_dHV_ave=1000,min_dHV_rms=1000;
	for(int i=0;i<1024;i++){
		if(temp_n[i]==0){continue;}
		temp_ave[i] /= temp_n[i];
		temp_rms[i] = sqrt(temp_rms[i]/temp_n[i] - pow(temp_ave[i],2));
		clb_temp_ave[i] /= temp_n[i];
		clb_temp_rms[i] = sqrt(clb_temp_rms[i]/temp_n[i] - pow(clb_temp_ave[i],2));
		dHV_ave[i] /= temp_n[i];
		dHV_rms[i] /= temp_n[i];
		dHV_rms[i] = sqrt(dHV_rms[i]-pow(dHV_ave[i],2));

		if( !((tel==7&&i==780)||(tel==10&&i==673)) ){
			if(max_dHV_ave<dHV_ave[i]){max_dHV_ave=dHV_ave[i];}
			if(min_dHV_ave>dHV_ave[i]){min_dHV_ave=dHV_ave[i];}
			if(max_dHV_rms<dHV_rms[i]){max_dHV_rms=dHV_rms[i];}
			if(min_dHV_rms>dHV_rms[i]){min_dHV_rms=dHV_rms[i];}
		}
		if(max_clb_temp_ave<clb_temp_ave[i]){max_clb_temp_ave=clb_temp_ave[i];}
		if(min_clb_temp_ave>clb_temp_ave[i]){min_clb_temp_ave=clb_temp_ave[i];}
		if(max_clb_temp_rms<clb_temp_rms[i]){max_clb_temp_rms=clb_temp_rms[i];}
		if(min_clb_temp_rms>clb_temp_rms[i]){min_clb_temp_rms=clb_temp_rms[i];}
		if(tel==4){if(i==268||i==719||i==944){continue;}}
		if(max_temp_ave<temp_ave[i]){max_temp_ave=temp_ave[i];}
		if(min_temp_ave>temp_ave[i]){min_temp_ave=temp_ave[i];}
		if(max_temp_rms<temp_rms[i]){max_temp_rms=temp_rms[i];}
		if(min_temp_rms>temp_rms[i]){min_temp_rms=temp_rms[i];}
	}
	if(min_temp_ave<-20) min_temp_ave = -20;
	if(max_temp_ave>45) max_temp_ave = 45;
	double temp_ave_range = max_temp_ave-min_temp_ave;
	if(min_temp_rms<-20) min_temp_rms = -20;
	if(max_temp_rms>45) max_temp_rms = 45;
	double temp_rms_range = max_temp_rms-min_temp_rms;
	if(min_clb_temp_ave<-20) min_clb_temp_ave = -20;
	if(max_clb_temp_ave>45) max_clb_temp_ave = 45;
	double clb_temp_ave_range = max_clb_temp_ave-min_clb_temp_ave;
	if(min_clb_temp_rms<-20) min_clb_temp_rms = -20;
	if(max_clb_temp_rms>45) max_clb_temp_rms = 45;
	double clb_temp_rms_range = max_clb_temp_rms-min_clb_temp_rms;
	if(min_dHV_ave<-2) min_dHV_ave = -2;
	if(max_dHV_ave>2) max_dHV_ave = 2;
	double dHV_ave_range = max_dHV_ave-min_dHV_ave;
	if(min_dHV_rms<-2) min_dHV_rms = -2;
	if(max_dHV_rms>2) max_dHV_rms = 2;
	double dHV_rms_range = max_dHV_rms-min_dHV_rms;
	int blocks = 40,half_blocks = 20;

	hv_temp->Fill();
	//status->Fill();

	/*draw pictures*/
	//draw work_mode, trigger_mode and firedtube
	graphSet(g_work_mode, "g_work_mode", "work_mode", time0);
	graphSet(g_triggermode, "g_triggermode", "trigger_mode", time0);
	graphSet(g_firedtube, "g_firedtube", "firedtube", time0);

	TCanvas *canvas = new TCanvas("canvas","canvas",1200,1200);
	canvas->SetMargin(0.12,0.12,0.12,0.12);
	canvas->Divide(1,3,0,0);
	canvas->cd(1);
	gPad->SetGridx();
	g_work_mode->Draw("ap");
	TLegend *legend1 = new TLegend(.5,.7,1,1);
	legend1->Clear();
	legend1->AddEntry(g_work_mode,"32+#:R_F#   48:Normal   64:R_F9plus");
	legend1->Draw();
	canvas->cd(2);
	gPad->SetGridx();
	g_triggermode->Draw("ap");
	TLegend *legend2 = new TLegend(.5,.7,1,1);
	legend2->Clear();
	legend2->AddEntry(g_triggermode,"48:tube_trigger   32:pattern_trigger   16:always_trigger");
	legend2->Draw();
	canvas->cd(3);
	gPad->SetGridx();
	g_firedtube->Draw("ap");
	canvas->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf(",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");


	//draw tel average sipm temperature and average clb temperature
	threshSet(g_tel_HV, "g_tel_HV", 1, 8, 1, time0);
	TCanvas *c_tel_hv = new TCanvas("c_tel_hv","c_tel_hv",1200,800);
	c_tel_hv->SetMargin(0.12,0.12,0.12,0.12);
	c_tel_hv->cd();
	gPad->SetGridx();
	gPad->SetGridy();
	g_tel_HV->GetYaxis()->SetTitle("hv/V");
	g_tel_HV->Draw("ap");
	c_tel_hv->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");

	//draw tel average sipm temperature and average clb temperature
	threshSet(g_tel_sipm_T, "g_tel_sipm_T", 1, 8, 1, time0);
	threshSet(g_tel_clb_T, "g_tel_clb_T", 2, 8, 1, time0);

	TCanvas *c_tel_temp = new TCanvas("c_tel_temp","c_tel_temp",1200,800);
	c_tel_temp->SetMargin(0.12,0.12,0.12,0.12);
	c_tel_temp->cd();
	gPad->SetGridx();
	gPad->SetGridy();
	g_tel_sipm_T->GetYaxis()->SetTitle("temperature/^{o}C");
	g_tel_sipm_T->Draw("ap");
	g_tel_clb_T->Draw("psame");
	TLegend *leg_tel_T = new TLegend(.5,.7,0.88,0.88);
	leg_tel_T->Clear();
	leg_tel_T->AddEntry(g_tel_sipm_T,"average_sipm_temperature");
	leg_tel_T->AddEntry(g_tel_clb_T,"average_clb_temperature");
	leg_tel_T->Draw();
	c_tel_temp->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");

	//draw thresh
	threshSet(g_cen_th, "g_cen_th", 1, 8, 1.4, time0);
	threshSet(g_cen_n_th, "g_cen_n_th", 2, 8, 1.2, time0);
	threshSet(g_out_th, "g_out_th", 3, 8, 1, time0);
	threshSet(g_out_n_th, "g_out_n_th", 4, 8, 0.8, time0);
	threshSet(g_rec_th, "g_rec_th", 5, 8, 0.6, time0);
	TCanvas *c_thresh = new TCanvas("c_thresh","c_thresh",1200,800);
	c_thresh->SetMargin(0.12,0.12,0.12,0.12);
	c_thresh->cd();
	gPad->SetGridx();
	gPad->SetGridy();
	g_cen_th->GetYaxis()->SetRangeUser(0,1500);
	g_cen_th->Draw("ap");
	if(noisSiPM.size()>0) g_cen_n_th->Draw("psame");
	g_out_th->Draw("psame");
	if(out_noisSiPM.size()>0) g_out_n_th->Draw("psame");
	g_rec_th->Draw("psame");
	TLegend *legend_th = new TLegend(.5,.7,.88,.88);
	legend_th->Clear();
	legend_th->AddEntry(g_cen_th,"center_thresh");
	if(noisSiPM.size()>0) legend_th->AddEntry(g_cen_n_th,"center_hotpixel_thresh");
	legend_th->AddEntry(g_out_th,"out_thresh");
	if(out_noisSiPM.size()>0) legend_th->AddEntry(g_out_n_th,"out_hotpixel_thresh");
	legend_th->AddEntry(g_rec_th,"record_thresh");
	legend_th->Draw();
	c_thresh->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");

	TBox *bx[1024];
	TBox *Bx[50];
	//draw temperature
	TH2D *th = new TH2D("","",100,-9,9,100,-9,9);
	TCanvas *c_temp = new TCanvas("c_temp","maps",1600,800);
	c_temp->Divide(2,1);
	c_temp->cd(1);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	th->Draw("");
	th->SetStats(kFALSE);
	TPaveText *pave1 = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave1->SetFillStyle(0);
	pave1->AddText(Form("temperature_ave"));
	pave1->Draw();
	for(int k=0;k<1024;k++){
		bx[k] = new TBox(PMTMAP[0][k]-0.25,PMTMAP[1][k]-0.25,PMTMAP[0][k]+0.25,PMTMAP[1][k]+0.25);
		if(temp_ave[k]==0) continue;
		int color = int(((temp_ave[k]-min_temp_ave)/temp_ave_range)*blocks)+55;
		bx[k]->SetFillColor(color);
		bx[k]->SetLineStyle(1);
		bx[k]->Draw();
	}
	TGaxis* axis1 = new TGaxis(9.1,-8,9.1,8,min_temp_ave,max_temp_ave,half_blocks,"+LS");
	axis1->SetTickSize(0.02);
	axis1->SetTitle("sipm_temperature/^{o}C");
	axis1->SetTitleOffset(1.5);
	for(int i=0;i<half_blocks;i++)
	{
		Bx[i] = new TBox(9.1,-8+i*(16./half_blocks),9.6,-8+(i+1)*(16./half_blocks));
		int color = i*2+1+55;
		Bx[i]->SetFillColor(color);
		Bx[i]->Draw();
	}
	axis1->Draw();
	c_temp->cd(2);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	th->Draw("");
	th->SetStats(kFALSE);
	TPaveText *pave2 = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave2->SetFillStyle(0);
	pave2->AddText(Form("temperature_rms"));
	pave2->Draw();
	for(int k=0;k<1024;k++){
		bx[k] = new TBox(PMTMAP[0][k]-0.25,PMTMAP[1][k]-0.25,PMTMAP[0][k]+0.25,PMTMAP[1][k]+0.25);
		if(temp_rms[k]==0) continue;
		int color = int(((temp_rms[k]-min_temp_rms)/temp_rms_range)*blocks)+55;
		bx[k]->SetFillColor(color);
		bx[k]->SetLineStyle(1);
		bx[k]->Draw();
	}
	TGaxis* axis2 = new TGaxis(9.1,-8,9.1,8,min_temp_rms,max_temp_rms,half_blocks,"+LS");
	axis2->SetTickSize(0.02);
	for(int i=0;i<half_blocks;i++)
	{
		Bx[i] = new TBox(9.1,-8+i*(16./half_blocks),9.6,-8+(i+1)*(16./half_blocks));
//		Bx[i] = new TBox(9.1,9.1/half_blocks*i-8,9.6,9.1/half_blocks*(i+1));
		int color = i*2+1+55;
		Bx[i]->SetFillColor(color);
		Bx[i]->Draw();
	}
	axis2->Draw();
	c_temp->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");

	//draw clb_temperature
	TCanvas *c_clb_temp = new TCanvas("c_clb_temp","maps",1600,800);
	c_clb_temp->Divide(2,1);
	c_clb_temp->cd(1);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	th->Draw("");
	th->SetStats(kFALSE);
	TPaveText *pave_clb1 = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave_clb1->SetFillStyle(0);
	pave_clb1->AddText(Form("clb_temperature_ave"));
	pave_clb1->Draw();
	for(int k=0;k<1024;k++){
		bx[k] = new TBox(PMTMAP[0][k]-0.25,PMTMAP[1][k]-0.25,PMTMAP[0][k]+0.25,PMTMAP[1][k]+0.25);
		if(clb_temp_ave[k]==0) continue;
		int color = int(((clb_temp_ave[k]-min_clb_temp_ave)/clb_temp_ave_range)*blocks)+55;
		bx[k]->SetFillColor(color);
		bx[k]->SetLineStyle(1);
		bx[k]->Draw();
	}
	TGaxis* axis_clb1 = new TGaxis(9.1,-8,9.1,8,min_clb_temp_ave,max_clb_temp_ave,half_blocks,"+LS");
	axis_clb1->SetTickSize(0.02);
	axis_clb1->SetTitle("clb_temperature/^{o}C");
	axis_clb1->SetTitleOffset(1.5);
	for(int i=0;i<half_blocks;i++)
	{
		Bx[i] = new TBox(9.1,-8+i*(16./half_blocks),9.6,-8+(i+1)*(16./half_blocks));
		int color = i*2+1+55;
		Bx[i]->SetFillColor(color);
		Bx[i]->Draw();
	}
	axis_clb1->Draw();
	c_clb_temp->cd(2);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	th->Draw("");
	th->SetStats(kFALSE);
	TPaveText *pave_clb2 = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave_clb2->SetFillStyle(0);
	pave_clb2->AddText(Form("clb_temperature_rms"));
	pave_clb2->Draw();
	for(int k=0;k<1024;k++){
		bx[k] = new TBox(PMTMAP[0][k]-0.25,PMTMAP[1][k]-0.25,PMTMAP[0][k]+0.25,PMTMAP[1][k]+0.25);
		if(clb_temp_rms[k]==0) continue;
		int color = int(((clb_temp_rms[k]-min_clb_temp_rms)/clb_temp_rms_range)*blocks)+55;
		bx[k]->SetFillColor(color);
		bx[k]->SetLineStyle(1);
		bx[k]->Draw();
	}
	TGaxis* axis_clb2 = new TGaxis(9.1,-8,9.1,8,min_clb_temp_rms,max_clb_temp_rms,half_blocks,"+LS");
	axis_clb2->SetTickSize(0.02);
	for(int i=0;i<half_blocks;i++)
	{
		Bx[i] = new TBox(9.1,-8+i*(16./half_blocks),9.6,-8+(i+1)*(16./half_blocks));
		int color = i*2+1+55;
		Bx[i]->SetFillColor(color);
		Bx[i]->Draw();
	}
	axis_clb2->Draw();
	c_clb_temp->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");
	//c_clb_temp->SaveAs(Form("CheckResulte/%s.%04d%04d_tel%02d_Routinue.png",argv[4],time.year,time0.month*100+time0.day,tel));

	//draw delta_hv
	TCanvas *c_hv = new TCanvas("c_hv"," maps",1600,800);
	c_hv->Divide(2,2);
	c_hv->GetPad(1)->SetPad(0,0,0.5,1);
	c_hv->GetPad(2)->SetPad(0.5,0.4,1,0.9);
	c_hv->GetPad(3)->SetPad(0,0,0,0);
	c_hv->GetPad(4)->SetPad(0.5,0,1,0.4);
	c_hv->cd(1);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	th->Draw("");
	th->SetStats(kFALSE);
	TPaveText *pave3 = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave3->SetFillStyle(0);
	pave3->AddText(Form("delta_HV_ave"));
	pave3->Draw();
	for(int k=0;k<1024;k++){
		bx[k] = new TBox(PMTMAP[0][k]-0.25,PMTMAP[1][k]-0.25,PMTMAP[0][k]+0.25,PMTMAP[1][k]+0.25);
		if(dHV_ave[k]==0) continue;
		int color = int(((dHV_ave[k]-min_dHV_ave)/dHV_ave_range)*blocks)+55;
		bx[k]->SetFillColor(color);
		bx[k]->SetLineStyle(1);
		bx[k]->Draw();
	}
	TGaxis* axis3 = new TGaxis(9.1,-8,9.1,8,min_dHV_ave,max_dHV_ave,half_blocks,"+LS");
	axis3->SetTickSize(0.02);
	axis3->SetTitle("hv/V");
	axis3->SetTitleOffset(1.5);
	for(int i=0;i<half_blocks;i++)
	{
		Bx[i] = new TBox(9.1,-8+i*(16./half_blocks),9.6,-8+(i+1)*(16./half_blocks));
		int color = i*2+1+55;
		Bx[i]->SetFillColor(color);
		Bx[i]->Draw();
	}
	axis3->Draw();
	c_hv->cd(2);
	gPad->SetMargin(0.1,0.05,0,0.1);
	double Channel[1024]={0};
	for(int i=0;i<1024;i++) {	Channel[i] = i;}
	TGraph* hv_return = new TGraph(1024, Channel, HVReturn);
	TGraph* hv_setting = new TGraph(1024, Channel, HVSetting);
	hv_return->SetTitle(Form("HV Comparison"));
	hv_return->GetYaxis()->SetTitle("hv/V");
	hv_return->GetYaxis()->CenterTitle();
	hv_return->GetYaxis()->SetTitleSize(0.05);
	hv_return->GetYaxis()->SetRangeUser(55,65);
	hv_return->SetMarkerStyle(6);
	hv_return->SetMarkerColor(1);
	hv_setting->SetMarkerStyle(6);
	hv_setting->SetMarkerColor(2);
	hv_return->Draw("ap");
	hv_setting->Draw("psame");
	TLegend *legend_hv_compare = new TLegend(.55,.6,0.95,0.9);
	legend_hv_compare->Clear();
	legend_hv_compare->AddEntry(hv_return,"HV Return");
	legend_hv_compare->AddEntry(hv_setting,"HV Setting");
	legend_hv_compare->Draw("same");

	c_hv->cd(4);
	gPad->SetGridy();
	gPad->SetMargin(0.1,0.05,0.4,0);
	double Delta_HV[1024]={0};
	for(int i=0;i<1024;i++) {	Delta_HV[i] = HVReturn[i] - HVSetting[i];}
	TGraph* Delta_Hv = new TGraph(1024, Channel, Delta_HV);
	Delta_Hv->SetTitle(Form(""));
	Delta_Hv->GetXaxis()->SetTitle("Channel");
	Delta_Hv->GetXaxis()->CenterTitle();
	Delta_Hv->GetXaxis()->SetTitleSize(0.08);
	Delta_Hv->GetYaxis()->SetTitle("delta_hv/V");
	Delta_Hv->GetYaxis()->CenterTitle();
	Delta_Hv->GetYaxis()->SetTitleSize(0.05);
	Delta_Hv->GetYaxis()->SetRangeUser(-1,1);
	Delta_Hv->SetMarkerStyle(6);
	Delta_Hv->SetMarkerColor(3);
	Delta_Hv->Draw("ap");
	c_hv->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");
	/*
	// draw delta_HV_rms distribution////////////////////////////////////////////////////
	th->Draw("");
	th->SetStats(kFALSE);
	TPaveText *pave4 = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave4->SetFillStyle(0);
	pave4->AddText(Form("delta_HV_rms"));
	pave4->Draw();
	for(int k=0;k<1024;k++){
		bx[k] = new TBox(PMTMAP[0][k]-0.25,PMTMAP[1][k]-0.25,PMTMAP[0][k]+0.25,PMTMAP[1][k]+0.25);
		if(dHV_rms[k]==0) continue;
		int color = int(((dHV_rms[k]-min_dHV_rms)/dHV_rms_range)*blocks)+55;
		bx[k]->SetFillColor(color);
		bx[k]->SetLineStyle(1);
		bx[k]->Draw();
	}
	TGaxis* axis4 = new TGaxis(9.1,-8,9.1,8,min_dHV_rms,max_dHV_rms,half_blocks,"+LS");
	axis4->SetTickSize(0.02);
	for(int i=0;i<half_blocks;i++)
	{
		Bx[i] = new TBox(9.1,-8+i*(16./half_blocks),9.6,-8+(i+1)*(16./half_blocks));
		int color = i*2+1+55;
		Bx[i]->SetFillColor(color);
		Bx[i]->Draw();
	}
	axis4->Draw();
	c_hv->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");
	*/

	//draw mask
	TCanvas *c_mask = new TCanvas("c_mask","maps",1200,1200);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	th->Draw("");
	th->SetStats(kFALSE);
	TPaveText *pave_mask = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave_mask->SetFillStyle(0);
	pave_mask->AddText(Form("sipm_mask"));
	pave_mask->Draw();
	for(int k=0;k<1024;k++){
		bx[k] = new TBox(PMTMAP[0][k]-0.25,PMTMAP[1][k]-0.25,PMTMAP[0][k]+0.25,PMTMAP[1][k]+0.25);
		int color = MASK[k]+1;
		bx[k]->SetFillColor(color);
		bx[k]->SetLineStyle(1);
		bx[k]->Draw();
	}
	c_mask->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");

	//draw db and clb fpga version
	TH2D *th_address = new TH2D("","",10,0,10,10,0,10);
	TCanvas *c_faddress = new TCanvas("c_faddress","maps",1200,1200);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	th_address->Draw("");
	th_address->SetStats(kFALSE);
	TPaveText *pave_fad = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave_fad->SetFillStyle(0);
	pave_fad->AddText(Form("fpga_address"));
	pave_fad->Draw();
	for(int k=0;k<64;k++){
		//printf("%d %d\n",SC_X[k]*10+SC_Y[k],SC_DB[k]*10+SC_F[k]);
		bx[k] = new TBox(SC_X[k],SC_Y[k],SC_X[k]+0.9,SC_Y[k]+0.9);
		int color = db_ver[k]-200;
		bx[k]->SetFillColor(color);
		bx[k]->Draw();
		TText *text = new TText(SC_X[k],SC_Y[k],Form("%d",SC_DB[k]*10+SC_F[k]));
		text->Draw();
	}
	//c_faddress->Print(Form("CheckResulte/%04d%04d_tel%02d_Routinue.pdf",time.year,time0.month*100+time0.day,tel),"pdf");
	//
	//
	set<int> db_ver_scal;
	set<int> clb_ver_scal;
	for(int i=0;i<64;i++){
		db_ver_scal.insert(db_ver[i]);
		clb_ver_scal.insert(clb_ver[i]);
	}
	printf("db_ver_num:%d clb_ver_num:%d\n",db_ver_scal.size(),clb_ver_scal.size());

	TH2D *th_ver = new TH2D("","",10,0,10,10,0,10);
	TCanvas *c_fver = new TCanvas("c_fver","maps",1600,800);
	c_fver->Divide(2,1);
	c_fver->cd(1);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	th_ver->Draw("");
	th_ver->SetStats(kFALSE);
	TPaveText *pave_db = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave_db->SetFillStyle(0);
	pave_db->AddText(Form("db_version"));
	pave_db->Draw();
	for(int k=0;k<64;k++){
		bx[k] = new TBox(SC_X[k],SC_Y[k],SC_X[k]+0.9,SC_Y[k]+0.9);
		int color = db_ver[k]-200;
		bx[k]->SetFillColor(color);
		bx[k]->Draw();
	}
	int db_ver_scal_num=0;
	for(set<int>::iterator it=db_ver_scal.begin();it!=db_ver_scal.end();it++){
		bx[db_ver_scal_num] = new TBox(10,db_ver_scal_num+1,10+0.9,db_ver_scal_num+1+0.9);
		int color = *it-200;
		bx[db_ver_scal_num]->SetFillColor(color);
		bx[db_ver_scal_num]->Draw();
		TText *text = new TText(11,db_ver_scal_num+1.2,Form("%x",*it));
		text->Draw();
	}
	c_fver->cd(2);
	gPad->SetMargin(0.15,0.15,0.15,0.15);
	th_ver->Draw("");
	th_ver->SetStats(kFALSE);
	TPaveText *pave_clb = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave_clb->SetFillStyle(0);
	pave_clb->AddText(Form("clb_version"));
	pave_clb->Draw();
	for(int k=0;k<64;k++){
		bx[k] = new TBox(SC_X[k],SC_Y[k],SC_X[k]+0.9,SC_Y[k]+0.9);
		int color = clb_ver[k];
		bx[k]->SetFillColor(color);
		bx[k]->Draw();
	}
	int clb_ver_scal_num=0;
	for(set<int>::iterator it=clb_ver_scal.begin();it!=clb_ver_scal.end();it++){
		bx[clb_ver_scal_num] = new TBox(10,clb_ver_scal_num+1,10+0.9,clb_ver_scal_num+1+0.9);
		int color = *it;
		//printf("%d clb_ver_color:%d\n",clb_ver_scal_num+1,color);
		bx[clb_ver_scal_num]->SetFillColor(color);
		bx[clb_ver_scal_num]->Draw();
		TText *text = new TText(11,clb_ver_scal_num+1.2,Form("%x",*it));
		text->Draw();
		clb_ver_scal_num++;
	}
	c_fver->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");

	//draw f1-8, f9 and f9plus fpga version
	TH2D *th_fver = new TH2D("","",10,0,10,5,0,5);
	TCanvas *c_ffver = new TCanvas("c_ffver","maps",1600,800);
	th_fver->Draw("");
	th_fver->SetStats(kFALSE);
	TPaveText *pave_f = new TPaveText(0.15,0.92,0.85,0.98,"NDC");
	pave_f->SetFillStyle(0);
	pave_f->AddText(Form("f_version"));
	pave_f->Draw();
	for(int k=0;k<10;k++){
		double x,y;
		if(k==0){	x=4.5;	y=2.5;}
		else if(k==9){	x=4.5;	y=1.5;}
		else{	x=k-(k/5*4)+k/(3+k/5*4)*4;	y=1+k/5*2;}
		bx[k] = new TBox(x,y,x+0.9,y+0.9);
		int color = f_ver[k]-180;
		bx[k]->SetFillColor(color);
		bx[k]->Draw();
		TText *text = new TText(x+0.1,y+0.25,Form("F%d:%x",k,f_ver[k]));
		text->Draw();
	}
	c_ffver->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf)",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");



//	TCanvas* cnew=new TCanvas();
//	cnew->Print(Form("%s/%04d%04d_tel%02d_Routinue.pdf)",argv[3],time.year,time0.month*100+time0.day,tel),"pdf");
//	delete cnew;

	g_work_mode->Write();
	g_triggermode->Write();
	g_firedtube->Write();
	hv_temp->Write();
	status->Write();
	fpga_ver->Write();
	fout->Close();

	return 0;
}

void rabbittime2lt(Long64_t rabbitTime, double rabbittime, Time &time)
{
	double mjd;
	mjd = MJD19700101 + (rabbitTime + rabbittime*20/1000000000. - TAI2UTC)/86400;

	int j;
	double fd, d;
	long jd, n4, nd10;
	/* Check if date is acceptable */
	if ( ( mjd <= -2395520.0 ) || ( mjd >= 1e9 ) ) {
		j = -1;
	}
	else {
		j = 0;
		/* Separate day and fraction */
		fd = (mjd)>0.0?mjd-floor(mjd):mjd+floor(-mjd);
		if ( fd < 0.0 ) fd += 1.0;
		d = mjd - fd;
		d = d<0.0?ceil(d):floor(d);
		/* Express day in Gregorian calendar */
		jd = (long)d + 2400001;
		n4 = 4L*(jd+((6L*((4L*jd-17918L)/146097L))/4L+1L)/2L-37L);
		nd10 = 10L*(((n4-237L)%1461L)/4L)+5L;
		time.year = (int) (n4/1461L-4712L);
		time.month = (int) (((nd10/306L+2L)%12L)+1L);
		time.day = (int) ((nd10%306L)/10L+1L);
		j = 0;
		time.hour = int(fd * 24 + 8);
		time.minite = int((fd*24+8 - time.hour)*60);
		time.second = int(((fd*24+8-time.hour)*60-time.minite)*60);
		time.day += time.hour/24;
		time.hour = time.hour%24;
	}

}

void graphSet(TGraph *graph, string name, string title, const Time &time)
{
	graph->SetName(name.c_str());
	//if(settitle){graph->SetTitle(Form("DATE %d/%d/%d",time.year,time.month,time.day));}
	graph->SetMarkerStyle(7);
	graph->GetXaxis()->SetTitle(Form("Time (%04d-%02d-%02d)",time.year,time.month,time.day));
	graph->GetXaxis()->CenterTitle();
	//graph->GetXaxis()->SetTimeFormat(Form("%H:%M %F%d-%02d-01 00:00:00",time.year,time.month));
	graph->GetXaxis()->SetTimeFormat("%H:%M%F1994-12-31 16:00:00s0");
	graph->GetXaxis()->SetTitleSize(.06);
	graph->GetXaxis()->SetTimeDisplay(1);
	graph->GetYaxis()->SetTitle(title.c_str());
	graph->GetYaxis()->CenterTitle();
	graph->GetYaxis()->SetTitleSize(.06);
}
void threshSet(TGraph *graph, string name, int pointcolor, int markerstyle, double markersize, const Time &time)
{
	graph->SetName(name.c_str());
	graph->SetMarkerStyle(markerstyle);
	graph->SetMarkerSize(markersize);
	graph->SetMarkerColor(pointcolor);
	graph->GetXaxis()->SetTitle(Form("Time (%04d-%02d-%02d)",time.year,time.month,time.day));
	graph->GetXaxis()->CenterTitle();
	graph->GetXaxis()->SetTimeFormat("%H:%M%F1994-12-31 16:00:00s0");
	graph->GetXaxis()->SetTitleSize(.06);
	graph->GetXaxis()->SetTimeDisplay(1);
	//graph->GetYaxis()->SetTitle(title.c_str());
	//graph->GetYaxis()->CenterTitle();
	//graph->GetYaxis()->SetTitleSize(.06);
}

void invert(double *conf)
{
	double conf5[1024];
	for(int i=0;i<1024;i++){
		conf5[i] = conf[1023-i];
	}
	for(int i=0;i<1024;i++){
		conf[i] = conf5[i];
	}
}


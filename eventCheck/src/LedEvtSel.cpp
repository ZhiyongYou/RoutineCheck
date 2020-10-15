#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "LedEvtSel.h"

using namespace std;

LedEvtSel::LedEvtSel()
{
	Initial();
}

LedEvtSel::~LedEvtSel()
{
	Initial();
}

void LedEvtSel::Initial()
{
	led_sign1=0;
	led_sign2=0;
	qlevtsize=0;
}

void LedEvtSel::SetLed0( int itel )
{
	int debug = 0;
	FILE *fp_led;
	fp_led=fopen(Form("/workfs/ybj/yinlq/LHAASO/Match/CaliMerge_v4.1/califile/wfcta%02d.0122.led",itel), "r");
	if(fp_led==NULL) { cerr<<"Can not open led file"<<endl; return ;}

	float led_fac_h[1024];
	float led_fac_l[1024];
	int iled_fac;
	double value_led_fac_h,value_led_fac_l;
	while(!feof(fp_led)){
		fscanf(fp_led,"%d %lf %lf\n",&iled_fac,&value_led_fac_h,&value_led_fac_l);
		if(debug) printf("%d %lf %lf\n",iled_fac,value_led_fac_h,value_led_fac_l);
		led_fac_h[iled_fac] = value_led_fac_h;
		led_fac_l[iled_fac] = value_led_fac_l;
	}
	fclose(fp_led);

	vector<int> maskSiPM;
	maskSiPM.clear();
	if(itel==4){
		maskSiPM.push_back(614);
		maskSiPM.push_back(944);
	}

	SetMaskPix(maskSiPM);
	SetLedFac(led_fac_h,led_fac_l);
}

int LedEvtSel::Running( WFCTAEvent *wfctaevent )
{
	CalcLedPar(wfctaevent);

	int ledevent=99;
	if(l_npix>850&&ratio_rms<0.5&&lg_npix<20){
		ledevent=13;
	}
	else {
		ledevent=99;
	}

	return ledevent;

}

void LedEvtSel::SetMaskPix(vector<int> &masksipm)
{
	maskSiPM.clear();
	for(int ii=0;ii<masksipm.size();ii++){
		maskSiPM.push_back(masksipm.at(ii));
	}
}

void LedEvtSel::SetLedFac(float *led_h, float *led_l)
{
	for(int i=0;i<1024;i++){
		led_fac_h[i] = led_h[i];
		led_fac_l[i] = led_l[i];
	}
}

void LedEvtSel::CalcLedPar(const WFCTAEvent *wfctaevent)
{
	//Init
	adc.clear();
	ratio.clear();
	eLG_Sat=0;	LG1_Sat=0;	LG2_Sat=0;	LG3_Sat=0;	HG_Sat=0;
	ratioh2=0;	ratioh3=0;	ratiol2=0;	ratiol3=0;
	top6 = 0;	totalsize = 0;
	l_npix = 0;
	lg_npix = 0;
	ratio_ave=0;	ratio_rms=0;

	//led single cut
	double l_cut=1000;

	int mask;
	int isipm;
	for(int ii=0;ii<wfctaevent->iSiPM.size();ii++){
		isipm=wfctaevent->iSiPM.at(ii);
		mask=0;
		for(int jj=0;jj<maskSiPM.size();jj++){
			if(isipm==maskSiPM.at(jj)){
				mask=1;
			}
		}
		if(mask==1){continue;}

		double rat_h = wfctaevent->AdcH.at(ii)/led_fac_h[isipm];
		double rat_l = wfctaevent->AdcL.at(ii)/led_fac_l[isipm];
		if(rat_h>2) ratioh2++;
		if(rat_h>3) ratioh3++;
		if(rat_l>2) ratiol2++;
		if(rat_l>3) ratiol3++;
		if(wfctaevent->eSatL.at(ii)==1){eLG_Sat++;}
		if(wfctaevent->AdcL.at(ii)>1500){LG1_Sat++;}
		if(wfctaevent->AdcL.at(ii)>2000){LG2_Sat++;}
		if(wfctaevent->AdcL.at(ii)>3000){LG3_Sat++;}
		if(wfctaevent->eSatH.at(ii)==1 || wfctaevent->AdcH.at(ii)>7000){HG_Sat++;}

		if(wfctaevent->AdcL.at(ii)>600){lg_npix++;}
		if(wfctaevent->AdcH.at(ii)>7000||wfctaevent->winsum.at(ii)+4*wfctaevent->BaseH.at(ii)>11000){
			ratio.push_back( rat_l );
			adc.push_back(wfctaevent->AdcL.at(ii)*22.5);
			totalsize += wfctaevent->AdcL.at(ii)*22.5;
			if(wfctaevent->AdcL.at(ii)*22.5>l_cut){l_npix++;}
		}
		else{
			ratio.push_back( rat_h );
			adc.push_back(wfctaevent->AdcH.at(ii));
			totalsize += wfctaevent->AdcH.at(ii);
			if(wfctaevent->AdcH.at(ii)>l_cut){l_npix++;}
		}
	}

	//calc ratio_ave and ratio_rms
	sort(ratio.begin(),ratio.end());
	int ratio_num=0;
	int ratiosize = ratio.size();
	for(int ii=15;ii<ratiosize-15;ii++){
		ratio_ave += ratio.at(ii);
		ratio_rms += pow(ratio.at(ii),2);
		ratio_num++;
	}
	if(ratio_num!=0)  { ratio_ave /= ratio_num;     ratio_rms = sqrt(ratio_rms/ratio_num - pow(ratio_ave,2));}

	//calc top6
	sort(adc.begin(),adc.end());
	int len = adc.size();
	if(len>7){
		for(int ii=len-1;ii>len-7;ii--){    top6 += adc.at(ii);}
	}
	else{
		for(int ii=0;ii<adc.size();ii++){   top6 += adc.at(ii);}
	}
}

void LedEvtSel::FindLedEvents(const vector<event> &quasi_led_evts, vector<long> &led_evts)
{
	qlevtsize = quasi_led_evts.size();
	for(int ii=0;ii<quasi_led_evts.size();ii++){
		led_sign1=0;led_sign2=0;
		for(int ievt1=1;ievt1<16;ievt1++){
			if(ii-ievt1<0){break;}
			double dt1 =(quasi_led_evts.at(ii).Time - quasi_led_evts.at(ii-ievt1).Time)*1000000000 + 
				(quasi_led_evts.at(ii).time - quasi_led_evts.at(ii-ievt1).time)*20;
			double dt_ratio1 = (dt1 / 333333333.) - round(dt1 / 333333333.);
			if(abs(dt_ratio1)<0.0002){	led_sign1++;}
		}
		for(int ievt2=1;ievt2<16;ievt2++){
			if(ii+ievt2>=qlevtsize){break;}
			double dt2 =(quasi_led_evts.at(ii+ievt2).Time - quasi_led_evts.at(ii).Time)*1000000000 + 
				(quasi_led_evts.at(ii+ievt2).time - quasi_led_evts.at(ii).time)*20;
			double dt_ratio2 = (dt2 / 333333333.) - round(dt2 / 333333333.);
			if(abs(dt_ratio2)<0.0002){	led_sign2++;}
		}
		if(led_sign1>=4 || led_sign2>=4){
			led_evts.push_back(quasi_led_evts.at(ii).ientry);
		}
	}
}

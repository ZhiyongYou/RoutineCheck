#ifndef LEDEVENTSELECT_H
#define LEDEVENTSELECT_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "WFCTAEvent.h"

struct event{
	long ientry;
	long long Time;
	double time;
};

class LedEvtSel
{
	protected:
		int eLG_Sat;
		int LG1_Sat;
		int LG2_Sat;
		int LG3_Sat;
		int HG_Sat;
		int ratioh2;
		int ratioh3;
		int ratiol2;
		int ratiol3;
		double top6;
		double totalsize;
		int l_npix;
		int lg_npix;
		double ratio_ave;
		double ratio_rms;
		float led_fac_h[1024];
		float led_fac_l[1024];
		std::vector<double> adc;
		std::vector<double> ratio;
		std::vector<int> maskSiPM;

		int led_sign1;
		int led_sign2;
		int qlevtsize;

		void Initial();

	public:
		LedEvtSel();
		~LedEvtSel();
		Int_t L_Npix()  {return l_npix;}
		Int_t LG_Npix()	{return lg_npix;}
		Int_t eLG_SAT() {return eLG_Sat;}
		Int_t LG1_SAT() {return LG1_Sat;}
		Int_t LG2_SAT() {return LG2_Sat;}
		Int_t LG3_SAT() {return LG3_Sat;}
		Int_t HG_SAT()  {return HG_Sat;}
		Int_t RatioH2() {return ratioh2;}
		Int_t RatioH3() {return ratioh3;}
		Int_t RatioL2() {return ratiol2;}
		Int_t RatioL3() {return ratiol3;}
		Double_t TotalSize() {return totalsize;}
		Double_t Top6() {return top6;}
		Double_t Ratio_Ave() {return ratio_ave;}
		Double_t Ratio_Rms() {return ratio_rms;}

		void SetMaskPix(std::vector<int> &masksipm);
		void SetLedFac(float *led_h, float *led_l);
		void CalcLedPar(const WFCTAEvent *wfctaevent);
		void FindLedEvents(const std::vector<event> &quasi_led_evts, std::vector<long> &led_evts);

		void SetLed0( int itel );
		int Running( WFCTAEvent *wfctaevent );
};

#endif // LEDEVENTSELECT_H

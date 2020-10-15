#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "WFCTAEvent.h"

using namespace std;

ClassImp(WFCTAEvent);

WFCTAEvent::WFCTAEvent():TObject()
{
	int MAXPMT = 1024;
	packCheck.resize(MAXPMT);
	eevent.resize(MAXPMT);
	zipmod.resize(MAXPMT);
	iSiPM.resize(MAXPMT);
	winsum.resize(MAXPMT);

	eSatH.resize(MAXPMT);
	eSatL.resize(MAXPMT);

	BaseH.resize(MAXPMT);
	BaseL.resize(MAXPMT);
	LaserBaseH.resize(MAXPMT);
	LaserBaseL.resize(MAXPMT);
	BaseHRMS.resize(MAXPMT);
	BaseLRMS.resize(MAXPMT);
	LaserBaseHRMS.resize(MAXPMT);
	LaserBaseLRMS.resize(MAXPMT);
	AdcH.resize(MAXPMT);
	AdcL.resize(MAXPMT);
	LaserAdcH.resize(MAXPMT);
	LaserAdcL.resize(MAXPMT);

	SatH.resize(MAXPMT);
	SatL.resize(MAXPMT);

	PeakPosH.resize(MAXPMT);
	PeakPosL.resize(MAXPMT);
	PeakAmH.resize(MAXPMT);
	PeakAmL.resize(MAXPMT);
	Over_Single_Marker.resize(MAXPMT);
	Over_Record_Marker.resize(MAXPMT);

	EventInitial();

}

WFCTAEvent::~WFCTAEvent()
{
	EventInitial();
}

void WFCTAEvent::EventInitial()
{
	iTel=-1;
	merge_size=-1;
	iEvent=-1;
	eEvent=-1;
	rabbitTime=0;
	rabbittime=0;
	big_pack_lenth=-1;
	n_fired=-1;
	n_Channel=-1;
	packCheck.clear();
	iSiPM.clear();
	zipmod.clear();
	eevent.clear();
	winsum.clear();
	PeakPosH.clear();
	PeakPosL.clear();
	PeakAmH.clear();
	PeakAmL.clear();
	Over_Single_Marker.clear();
	Over_Record_Marker.clear();
	eSatH.clear();
	eSatL.clear();
	SatH.clear();
	SatL.clear();
	BaseH.clear();
	BaseL.clear();
	LaserBaseH.clear();
	LaserBaseL.clear();
	BaseHRMS.clear();
	BaseLRMS.clear();
	LaserBaseHRMS.clear();
	LaserBaseLRMS.clear();
	AdcH.clear();
	AdcL.clear();
	LaserAdcH.clear();
	LaserAdcL.clear();
}

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "LHAASOEvent.h"

using namespace std;

ClassImp(LHAASOEvent);
LHAASOEvent* LHAASOEvent::_Head=0;
TTree* LHAASOEvent::_Tree=0;
const char* LHAASOEvent::_Name="LHAASOEvent";
TBranch* LHAASOEvent::bAll=0;

LHAASOEvent::LHAASOEvent():TObject()
{
	isipm.resize(MAXPMT);
	sipmpe.resize(MAXPMT);
	sipmt.resize(MAXPMT);
	cellig.resize(WCDASIZE);
	cellpe.resize(WCDASIZE);
	cellt.resize(WCDASIZE);
	hitid.resize(KM2ASIZE);
	hitpart.resize(KM2ASIZE);
	hitt.resize(KM2ASIZE);
	hitmode.resize(KM2ASIZE);

	EventInitial();
}

LHAASOEvent::~LHAASOEvent()
{
	EventInitial();
}

void LHAASOEvent::EventInitial()
{
	rabbitTime=0;
	rabbittime=0;
	IsKm2aEvent=0;
	IsWcdaEvent=0;
	km2a_xc=-10000;
	km2a_yc=-10000;
	km2a_phi=-10000;
	km2a_theta=-10000;
	iwcdaevt=0;
	isipm.clear();
	sipmpe.clear();
	sipmt.clear();
	cellig.clear();
	cellpe.clear();
	cellt.clear();
	hitid.clear();
	hitpart.clear();
	hitt.clear();
	hitmode.clear();
}

void LHAASOEvent::SetWFCTAEvent(int tel, long long rb_Time, double rb_time, long long dt, std::vector<int> &sipm, std::vector<double> &sipmPe, std::vector<double> &sipmT)
{
	int ISIPM=0;
	rabbitTime = rb_Time;
	rabbittime = rb_time;
	for(int ii=0;ii<sipm.size();ii++)
	{
		ISIPM = (tel-1)*1024 + sipm.at(ii);
		isipm.push_back(ISIPM);
		sipmpe.push_back(sipmPe.at(ii));
		sipmt.push_back(sipmT.at(ii)+dt);
	}
}

void LHAASOEvent::SetWCDAEvent(int iswcdaevent, ULong64_t Iwcdaevt, std::vector<int> &cellIg, std::vector<double> &cellPe, std::vector<double> &cellT)
{
	IsWcdaEvent = iswcdaevent;
	iwcdaevt = Iwcdaevt;
	for(int ii=0;ii<cellIg.size();ii++)
	{
		cellig.push_back(cellIg.at(ii));
		cellpe.push_back(cellPe.at(ii));
		cellt.push_back(cellT.at(ii));
	}
}

void LHAASOEvent::SetKM2AEvent(	int iskm2aevent, float km2a_Xc, float km2a_Yc, float km2a_Phi, float km2a_Theta, 
								std::vector<int> &hitId, std::vector<double> &hitPart, std::vector<double> &hitT, std::vector<int> &hitMode)
{
	IsKm2aEvent = iskm2aevent;
	km2a_xc = km2a_Xc;
	km2a_yc = km2a_Yc;
	km2a_phi = km2a_Phi;
	km2a_theta = km2a_Theta;
	for(int ii=0;ii<hitId.size();ii++)
	{
		hitid.push_back(hitId.at(ii));
		hitpart.push_back(hitPart.at(ii));
		hitt.push_back(hitT.at(ii));
		hitmode.push_back(hitMode.at(ii));
	}
}




















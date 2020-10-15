#ifndef STATUS_H
#define STATUS_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "LHAASOMatch.h"

const int MAXPMT=1024;
//class Status : public TSelector
class Status : public TObject
{
	protected:
		static Status* _Head;
		static TTree* _Tree;
		static const char * _Name;
		static TBranch* bAll;

	public:
		long long rabbitTime;
		double rabbittime;
		int IsKm2aEvent;
		int IsWcdaEvent;
		float km2a_xc;
		float km2a_yc;
		float km2a_phi;
		float km2a_theta;
		ULong64_t iwcdaevt;

		std::vector<int> isipm;
		std::vector<double> sipmpe;
		std::vector<double> sipmt;
		std::vector<int> cellig;
		std::vector<double> cellpe;
		std::vector<double> cellt;
		std::vector<int> hitid;
		std::vector<double> hitpart;
		std::vector<double> hitt;
		std::vector<int> hitmode;

	public:
		Status();
		~Status();
		void EventInitial();

		Long64_t RabbitTime() { return rabbitTime; }
		Double_t Rabbittime() { return rabbittime; }

		void SetWFCTAEvent(int tel, long long rb_Time, double rb_time, long long dt, std::vector<int> &sipm, std::vector<double> &sipmPe, std::vector<double> &sipmT);
		void SetWCDAEvent(int iswcdaevent, ULong64_t Iwcdaevt, std::vector<int> &cellIg, std::vector<double> &cellPe, std::vector<double> &cellT);
		void SetKM2AEvent(	int iskm2aevent, float km2a_Xc, float km2a_Yc, float km2a_Phi, float km2a_Theta, 
							std::vector<int> &hitId, std::vector<double> &hitPart, std::vector<double> &hitT, std::vector<int> &hitMode);


		ClassDef(Status,1);
};

ClassImp(Status);

#endif // STATUS_H

#include <iostream>
#include <fstream>
#include <string.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "CollectIntensity.h"

CollectIntensity* CollectIntensity::m_instance = 0;

CollectIntensity::CollectIntensity()
{
	std::cout << "wfcta CollectIntensity constructor ===============>" << std::endl;
	intensity.resize(1024,0);
	ReadIntensityFile("/eos/user/w/wfcta/wfcta/wfcta/piwenxuan/LightSpot/SpotMC/CollectIntensity.root");
	std::cout << "wfcta CollectIntensity constructor <===============" << std::endl;
}

CollectIntensity::~CollectIntensity()
{

}

int CollectIntensity::ReadIntensityFile(const char* filename)
{
	std::cout<< "###### Open WFCTA Intensity file "<< filename<< std::endl;
	TFile* infile = new TFile();
	infile = TFile::Open(filename);
	if(!infile) {
		printf("%s does not exist %s\n",filename);
		delete infile;
		return 0;
	}
	if(infile->IsZombie()||infile->GetEND()<50) {
		printf("%s file error!!\n",filename);
		delete infile;
		return 0;
	}
	TTree* intens_tree = (TTree *)infile->Get("t1");
	if(intens_tree==nullptr) {
		printf("%s does not have mc_cosmic tree\n",filename);
		delete infile;
		return 0;
	}
	else
	{
		int SiPM;
		double trans;
		double num;
		std::cout<< "Read intensity tree"<< std::endl;
		intens_tree->SetBranchAddress("SiPM", &SiPM);
		intens_tree->SetBranchAddress("trans", &trans);
		intens_tree->SetBranchAddress("num", &num);
		int n_entry = intens_tree->GetEntries();
		for(int i=0;i<n_entry;i++)
		{
			intens_tree->GetEntry(i);
			intensity[SiPM] = trans;
		}
		std::cout << "###### File read finished" << std::endl;
		delete infile;
		return 1;
	}
}


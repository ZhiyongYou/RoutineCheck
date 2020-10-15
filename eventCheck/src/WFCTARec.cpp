#include <map>
#include "TMath.h"
#include "WFCTARec.h"

using namespace std;
bool cmp(const std::pair<int,double>& tube1, const std::pair<int,double>& tube2)
{
	return tube1.second < tube2.second;
}

WFCTARec::WFCTARec()
{

}

WFCTARec::~WFCTARec()
{

}

//_block: event
void WFCTARec::SetEventType(int evt_type)
{
	event_type = evt_type;
}

//_block: init
void WFCTARec::WFCTAInit()
{
	iSipm.clear();
	SipmPe.clear();
	SipmT.clear();
}

void WFCTARec::SetWFCTAEvent(const std::vector<int>& iisipm, const std::vector<double>& isipmpe, const std::vector<double>& isipmt)
{
	double intensity=1;
	std::vector<int> new_sipm_id;
	std::vector<double> sipmpe_new;
	for(int i=0;i<iisipm.size();i++)
	{
		int sipmid = iisipm.at(i);
		int new_sipm = GetSipmIdInSimulation(sipmid);

		new_sipm_id.push_back(new_sipm);
//		CollectIntensity::Instance()->GetIntensity(new_sipm,intensity);
		sipmpe_new.push_back(isipmpe.at(i)/intensity);
	}
	iSipm.clear();
	SipmPe.clear();
	SipmT.clear();
	for(int ii=0;ii<iisipm.size();ii++)
	{
		//if(!(iisipm.at(ii)>=0*1024&&iisipm.at(ii)<1*1024))continue;
		iSipm.push_back(new_sipm_id.at(ii));
		SipmPe.push_back(sipmpe_new.at(ii));
		SipmT.push_back(isipmt.at(ii));
	}
//	printf("In WFCTAREC sipm.size:%d\n",iSipm.size());
}

//_block: rec
void WFCTARec::TimeClean(double pecut)
{
	double t_pecut = 20;
	timeclean.clear();
	if(1==event_type) { //laser event
		for(int ii=0;ii<iSipm.size();ii++){
			if(SipmPe.at(ii)>=pecut) {
				timeclean.push_back(1);
			}
			else {
				timeclean.push_back(0);
			}
		}
	}
	else { //cosmic event
		//find most probable sipm time
		double maxpro_time=-1;
		int num_maxpro_time=-1;
		map<double,int> sipm_t;
		map<double,int>::iterator it_sipm_t;
		sipm_t.clear();
		for(int ii=0;ii<iSipm.size();ii++) {
			if(SipmPe.at(ii)>=t_pecut) {
				sipm_t[SipmT.at(ii)]++;}
		}
		for(it_sipm_t=sipm_t.begin(); it_sipm_t!=sipm_t.end(); it_sipm_t++) {
			if(num_maxpro_time < it_sipm_t->second) {
				num_maxpro_time = it_sipm_t->second;
				maxpro_time = it_sipm_t->first;
			}
		}

		//maxpro_time=800;
		//do time clean
		int timeclean_sipms=0;
		for(int ii=0;ii<iSipm.size();ii++) {
			if(SipmPe.at(ii)>=pecut){ 
				//if(SipmPe.at(ii)>=20) {
				if(SipmT.at(ii)<maxpro_time-240 || SipmT.at(ii)>maxpro_time+240) {
					timeclean.push_back(0);
				}
				else {
					timeclean.push_back(1);
					timeclean_sipms++;
				}
			}
			else {
				timeclean.push_back(0);
			}
			}
//			printf("maxpro_time:%.0lf, left %d sipms, pecut:%.2lf\n",maxpro_time, timeclean_sipms, pecut);
		}
	}


void WFCTARec::GroupClean(int g_cut)
{
	group_cut = g_cut;

	vector<int> pix_in_group;
	int cnt, SIPM0, SIPM1;
	double distance, azi1, zen1, azi0, zen0;

	groupclean.clear();
	size_t sipm_num = iSipm.size();
	groupclean.resize(sipm_num,0);

	//printf("group clean test1: iSipm.size:%d groupclean.size:%d\n",iSipm.size(),groupclean.size());
	for(int sipm=0;sipm<iSipm.size();sipm++)
	{
		if(groupclean.at(sipm)!=0){ //groupclean.at(sipm) != 0 means this sipm already matched in some group
			continue;}
		if(timeclean.at(sipm)==0){ //sipms that cleaned by WFCTARec::TimeClean()
			continue;}
		pix_in_group.clear();
		for(int ipix=0;ipix<iSipm.size();ipix++){   pix_in_group.push_back(0);}
		pix_in_group.at(sipm) = 1;

		//printf("group clean loop test ...\n");
		while(1)
		{
			cnt = 0;
			for(int ii=0;ii<iSipm.size();ii++){
				if(pix_in_group.at(ii)==0) continue;
				SIPM0 = iSipm.at(ii);
				//if(SIPM0>=1024){continue;}//for one telescope
				//if(SIPM0<2048||SIPM0>=3072){continue;}
				if(merge_field)
					WFCTAMap::Instance()->GetwSipmAZ(SIPM0, azi0, zen0);
				else
					WFCTAMap::Instance()->GetSipmAZ(SIPM0, azi0, zen0);
				for(int jj=0;jj<iSipm.size();jj++){
					SIPM1 = iSipm.at(jj);
					//if(SIPM1<2048||SIPM1>=3072){continue;}//for one telescope
					if(merge_field)
						WFCTAMap::Instance()->GetwSipmAZ(SIPM1, azi1, zen1);
					else
						WFCTAMap::Instance()->GetSipmAZ(SIPM1, azi1, zen1);
					//	WFCTAMap::Instance()->GetwSipmAZ(SIPM1, azi1, zen1);
					distance = GetSiPMDist(azi0,zen0,azi1,zen1); // input:rad    return:deg
					if(distance<=WFCTAMap::Instance()->GetMAXDIST() && pix_in_group.at(jj)==0 && timeclean.at(jj)!=0)  {pix_in_group.at(jj)=1;cnt++;}
				}
			}
			if(cnt==0){break;}
		}
		int n_group=0;
		for(int ii=0;ii<iSipm.size();ii++){
			if(pix_in_group.at(ii)==1)
				n_group++;
		}       
		for(int ii=0;ii<iSipm.size();ii++){
			if(pix_in_group.at(ii)==1)
			{
				if(n_group>=group_cut)
					groupclean.at(ii) = 1;
			}
		}       

	}
}
void WFCTARec::IslandClean()
{
	groupclean.clear();
	if(1==event_type) {
		for(int ii=0;ii<timeclean.size();ii++) {
			groupclean.push_back(timeclean.at(ii));
		}
	}
	else {
		int cnt, SIPM0, SIPM1;
		double distance, azi1, zen1, azi0, zen0;
		for(int ipix=0;ipix<iSipm.size();ipix++){   groupclean.push_back(0);}
		double MaxPe = -1000;
		int MaxPeSipm=-1000;
		for(int ii=0;ii<iSipm.size();ii++) {
			if(0==timeclean.at(ii)){continue;}
			if(MaxPe<SipmPe.at(ii)) {
				MaxPe=SipmPe.at(ii);
				MaxPeSipm=ii;
			}
		}
		if(MaxPeSipm!=-1000)
			groupclean.at(MaxPeSipm) = 1;

		while(1)
		{
			cnt = 0;
			for(int ii=0;ii<iSipm.size();ii++){
				if(groupclean.at(ii)==0) continue;
				SIPM0 = iSipm.at(ii);
				//if(SIPM0>=1024){continue;}//for one telescope
				//if(SIPM0<2048||SIPM0>=3072){continue;}
				if(merge_field)
					WFCTAMap::Instance()->GetwSipmAZ(SIPM0, azi0, zen0);
				else
					WFCTAMap::Instance()->GetSipmAZ(SIPM0, azi0, zen0);
				//WFCTAMap::Instance()->GetwSipmAZ(SIPM0, azi0, zen0);
				for(int jj=0;jj<iSipm.size();jj++){
					SIPM1 = iSipm.at(jj);
					//if(SIPM1<2048||SIPM1>=3072){continue;}//for one telescope
					if(merge_field)
						WFCTAMap::Instance()->GetwSipmAZ(SIPM1, azi1, zen1);
					else
						WFCTAMap::Instance()->GetSipmAZ(SIPM1, azi1, zen1);
					//WFCTAMap::Instance()->GetwSipmAZ(SIPM1, azi1, zen1);
					distance = GetSiPMDist(azi0,zen0,azi1,zen1); // input:rad    return:deg
					if(distance<=WFCTAMap::Instance()->GetMAXDIST() && groupclean.at(jj)==0 && timeclean.at(jj)!=0)  {groupclean.at(jj)=1;cnt++;}
				}
			}
			if(cnt==0){break;}
		}
	}

}

void WFCTARec::CalcMainTel(int main_tel)
{
	MainTel = -1;
	const int tels = WFCTAMap::Instance()->MaxTelId();
	for(int i=0;i<6;i++) 
	{
		pe_size[i] = 0;
		npix_size[i] = 0;
	}
	MaxSipmPe=-1000;
	MaxSipm=-1;
	MaxPeTel=-1;
	for(int ii=0;ii<iSipm.size();ii++)
	{
		if(0==groupclean.at(ii)) continue;
		int SIPM = iSipm.at(ii);
		int itel = SIPM/1024;
		pe_size[itel] += SipmPe.at(ii);
		npix_size[itel]++;
		if(MaxSipmPe<SipmPe.at(ii))
		{
			MaxSipmPe = SipmPe.at(ii);
			MaxSipm = SIPM;
			MaxPeTel = itel+1;
		}
	}
	//	printf("MaxSipm:%d MaxPeTel:%d MaxSipmPe:%lf\n",MaxSipm, MaxPeTel, MaxSipmPe);
	int maintel_maxpe = MaxPeTel;
	int maintel_pesize = MaxPeTel;
	int maintel_npixsize = MaxPeTel;
	double max_pesize=-1000;
	double max_npixsize=-1000;
	for(int itel=0;itel<tels;itel++)
	{
		//printf("in_pe_size:%lf\n",pe_size[itel]);
		//printf("in_npix_size:%d\n",npix_size[itel]);
		if(max_pesize<pe_size[itel])
		{
			max_pesize=pe_size[itel];
			maintel_pesize = itel + 1;
		}
		if(max_npixsize<npix_size[itel])
		{
			max_npixsize=npix_size[itel];
			maintel_npixsize = itel + 1;
		}
	}
	if(maintel_pesize==maintel_npixsize)
		MainTel = maintel_pesize;
	else
		MainTel = MaxPeTel;
	//MainTel = main_tel; //test for events that triggered tel01 and te02
	//printf("maintel:%d max_pesize:%lf\n",MainTel,max_pesize);
}

void WFCTARec::MergeEvent(int do_clean)
{
	int jj;
	double tel_azi, tel_zen;
	double sipm_azi, sipm_zen;
	double sipm_x, sipm_y;
	double y_min,y_max;
	double x_min,x_max;
	WFCTAMap::Instance()->GetYRange(y_min, y_max);
	//printf("y_min:%lf y_max:%lf\n",y_min*TMath::RadToDeg(),y_max*TMath::RadToDeg());

	//test
	WFCTAMap::Instance()->GetTelAZ(comp_tel2, tel_azi, tel_zen);
	int cross_sipm_num=0;
	for(int ii=(comp_tel1-1)*1024;ii<comp_tel1*1024;ii++)
	{
		int SIPM = ii;
		int sipm = SIPM%1024;
		if(merge_field)
			WFCTAMap::Instance()->GetwSipmAZ(SIPM, sipm_azi, sipm_zen);
		else
			WFCTAMap::Instance()->GetSipmAZ(SIPM, sipm_azi, sipm_zen);
		slaDs2tp(sipm_azi, 90*TMath::DegToRad()-sipm_zen, tel_azi*TMath::DegToRad(), (90-tel_zen)*TMath::DegToRad(), 
				&sipm_x, &sipm_y, &jj);
		int iline = int(round(31-(sipm_y-y_min)*TMath::RadToDeg()/0.51));
		if(iline<0 || iline>32) continue;
		WFCTAMap::Instance()->GetXRange(iline, x_min, x_max);

		if(sipm_x<x_max+0.51*TMath::DegToRad()/2 && sipm_x>x_min-0.51*TMath::DegToRad()/2 && sipm_y<y_max+0.51*TMath::DegToRad()/2 && sipm_y>y_min-0.51*TMath::DegToRad()/2)  // interval = 1.0 mm
		{
			cross_sipm_num++;
	//		printf("cross_sipm_num:%d sipm:%4d(%4d) line:%d sipm_xy:(%lf,%lf) x_min:%lf x_max:%lf dy:%lf DY:%lf iline:%d\n",
	//				cross_sipm_num, SIPM,SIPM%1024,sipm/32,sipm_x*TMath::RadToDeg(),sipm_y*TMath::RadToDeg(),
	//				x_min*TMath::RadToDeg(),x_max*TMath::RadToDeg(),
	//				(sipm_y-y_min)*TMath::RadToDeg(),(sipm_y-y_min)*TMath::RadToDeg()/(31-sipm/32),iline);
		}
	}
	cross_sipm_num=0;
	WFCTAMap::Instance()->GetTelAZ(comp_tel1, tel_azi, tel_zen);
	for(int ii=(comp_tel2-1)*1024;ii<comp_tel2*1024;ii++)
	{
		int SIPM = ii;
		int sipm = SIPM%1024;
		if(merge_field)
			WFCTAMap::Instance()->GetwSipmAZ(SIPM, sipm_azi, sipm_zen);
		else
			WFCTAMap::Instance()->GetSipmAZ(SIPM, sipm_azi, sipm_zen);
		slaDs2tp(sipm_azi, 90*TMath::DegToRad()-sipm_zen, tel_azi*TMath::DegToRad(), (90-tel_zen)*TMath::DegToRad(), 
				&sipm_x, &sipm_y, &jj);
		int iline = int(round(31-(sipm_y-y_min)*TMath::RadToDeg()/0.51));
		if(iline<0 || iline>32) continue;
		WFCTAMap::Instance()->GetXRange(iline, x_min, x_max);

		if(sipm_x<x_max+0.51*TMath::DegToRad()/2 && sipm_x>x_min-0.51*TMath::DegToRad()/2 && sipm_y<y_max+0.51*TMath::DegToRad()/2 && sipm_y>y_min-0.51*TMath::DegToRad()/2)  // interval = 1.0 mm
		{
			cross_sipm_num++;
	//		printf("cross_sipm_num:%d sipm:%4d(%4d) line:%d sipm_xy:(%lf,%lf) x_min:%lf x_max:%lf dy:%lf DY:%lf iline:%d\n",
	//				cross_sipm_num,SIPM,SIPM%1024,sipm/32,sipm_x*TMath::RadToDeg(),sipm_y*TMath::RadToDeg(),
	//				x_min*TMath::RadToDeg(),x_max*TMath::RadToDeg(),
	//				(sipm_y-y_min)*TMath::RadToDeg(),(sipm_y-y_min)*TMath::RadToDeg()/(31-sipm/32),iline);
		}
	}




	for(int i=0;i<6;i++) inter_npix[i] = 0;
	for(int ii=0;ii<iSipm.size();ii++)
	{
		//if(0==groupclean.at(ii)) continue;
		if(0==timeclean.at(ii)) continue;
		int SIPM = iSipm.at(ii);
		int iTel = SIPM/1024 + 1;
		if(merge_field)
			WFCTAMap::Instance()->GetwSipmAZ(SIPM, sipm_azi, sipm_zen);
		else
			WFCTAMap::Instance()->GetSipmAZ(SIPM, sipm_azi, sipm_zen);

		if(iTel==comp_tel1)
			WFCTAMap::Instance()->GetTelAZ(comp_tel2, tel_azi, tel_zen);
		else if(iTel==comp_tel2)
			WFCTAMap::Instance()->GetTelAZ(comp_tel1, tel_azi, tel_zen);
		else
			continue;
		
		slaDs2tp(sipm_azi, 90*TMath::DegToRad()-sipm_zen, tel_azi*TMath::DegToRad(), (90-tel_zen)*TMath::DegToRad(), 
				&sipm_x, &sipm_y, &jj);

		int iline = int(round(31-(sipm_y-y_min)*TMath::RadToDeg()/0.51));
		if(iline<0 || iline>32) continue;
		WFCTAMap::Instance()->GetXRange(iline, x_min, x_max);


		if(sipm_x<x_max+0.51*TMath::DegToRad()/2 && sipm_x>x_min-0.51*TMath::DegToRad()/2 && sipm_y<y_max+0.51*TMath::DegToRad()/2 && sipm_y>y_min-0.51*TMath::DegToRad()/2)  // interval = 1.0 mm
		{
//			printf("event::::::::sipm:%4d(%4d) line:%d sipm_xy:(%lf,%lf) x_min:%lf x_max:%lf dy:%lf DY:%lf iline:%d\n",
//					SIPM,SIPM%1024,(SIPM%1024)/32,sipm_x*TMath::RadToDeg(),sipm_y*TMath::RadToDeg(),
//					x_min*TMath::RadToDeg(),x_max*TMath::RadToDeg(),
//					(sipm_y-y_min)*TMath::RadToDeg(),(sipm_y-y_min)*TMath::RadToDeg()/(31-(SIPM%1024)/32),iline);
			inter_npix[iTel-1]++;
		}
	}
//	for(int i=0;i<6;i++) printf("%dinter_npix: %d\n",i+1,inter_npix[i]);

	WFCTAMap::Instance()->GetTelAZ(MainTel, tel_azi, tel_zen);
	for(int ii=0;ii<iSipm.size();ii++)
	{
		//if(0==groupclean.at(ii)) continue;
		if(0==timeclean.at(ii)) continue;
		int SIPM = iSipm.at(ii);
		int iTel = SIPM/1024 + 1;
		if(merge_field)
			WFCTAMap::Instance()->GetwSipmAZ(SIPM, sipm_azi, sipm_zen);
		else
			WFCTAMap::Instance()->GetSipmAZ(SIPM, sipm_azi, sipm_zen);
		slaDs2tp(sipm_azi, 90*TMath::DegToRad()-sipm_zen, tel_azi*TMath::DegToRad(), (90-tel_zen)*TMath::DegToRad(), 
				&sipm_x, &sipm_y, &jj);

		int iline = int(round(31-(sipm_y-y_min)*TMath::RadToDeg()/0.51));
		if(iline<0 || iline>32) continue;
		WFCTAMap::Instance()->GetXRange(iline, x_min, x_max);

		//if(iTel!=MainTel && sipm_x<x_max && sipm_x>x_min && sipm_y<y_max && sipm_y>y_min)  // interval = 1.0 mm
		if(iTel!=MainTel && sipm_x<x_max+0.51*TMath::DegToRad()/2 && sipm_x>x_min-0.51*TMath::DegToRad()/2 && sipm_y<y_max+0.51*TMath::DegToRad()/2 && sipm_y>y_min-0.51*TMath::DegToRad()/2)  // interval = 1.0 mm
		{
			if(do_clean)
			{
			//	groupclean.at(ii) = 0;
				timeclean.at(ii) =0;
			}
		}
	}
}

void WFCTARec::CalcHillas()
{
	double zen,azi ,npe;
	//MainTel = -1;
	DNpix = 0;
	DSize = 0;
	DMeanZen = 0;
	DMeanAzi = 0;
	DLength = -1000; //on focus
	DWidth = -1000; //on focus
	DDelta = -1000; //on focus
	DAlpha=-1000; //

	//calc npix, size, mean_azimuth,mean_zenith
	for(int ii=0;ii<iSipm.size();ii++)
	{
		//if(groupclean.at(ii)<group_cut) continue;
		if(0==groupclean.at(ii)) continue;
		int SIPM = iSipm.at(ii);
		if(merge_field)
			WFCTAMap::Instance()->GetwSipmAZ(SIPM, azi, zen);
		else
			WFCTAMap::Instance()->GetSipmAZ(SIPM, azi, zen);
		//WFCTAMap::Instance()->GetwSipmAZ(SIPM, azi, zen);
		npe = SipmPe.at(ii);
		DSize += npe;
		DNpix++;
//		DMeanAzi += npe*azi;
//		DMeanZen += npe*zen;
		//printf("SIPM:%d %lf %lf %lf %lf\n",SIPM,azi,zen,npe,SipmT.at(ii));
	}
	if(0==DSize) return ;
	//DMeanAzi /= DSize;
	//DMeanZen /= DSize;

	/*
	//find main telescope
	double D_centroid2telaxis=100000;
	double tel_azi, tel_zen;
	for(int itel=0;itel<WFCTAMap::Instance()->MaxTelId();itel++) {
		WFCTAMap::Instance()->GetTelAZ(itel+1, tel_azi, tel_zen);
		double d_cen2axis = GetSiPMDist(DMeanAzi,DMeanZen, tel_azi*TMath::DegToRad(), tel_zen*TMath::DegToRad());
		if(D_centroid2telaxis>d_cen2axis) {
			D_centroid2telaxis = d_cen2axis;
			MainTel = itel+1;//Tel_Num[itel];
			//if(MainTel>6) MainTel = 6;
		}
	}
	*/

	GetEventMapOnFocus(MainTel);
	CalcHillasOnFocus(MainTel);
//	printf("hillas: MainTel:%d DNpix:%d DSize:%.4lf Centroid(%.4lf,%.4lf)\n",MainTel,DNpix,DSize,DMeanAzi*57.3,DMeanZen*57.3);

}

void WFCTARec::CalcSDP()
{
	OuterCentroid();
	CalcCoordsAloneSDP();

	double sdp_azi1, sdp_azi2;
	double sdp_zen1, sdp_zen2;
	double km2a_sdp_azi1, km2a_sdp_azi2;
	double km2a_sdp_zen1, km2a_sdp_zen2;
	double wcda_sdp_azi1, wcda_sdp_azi2;
	double wcda_sdp_zen1, wcda_sdp_zen2;

	//get tail by using OuterCentoid
	double S_M_X = DOutMeanX-DMeanX;
	double S_M_Y = DOutMeanY-DMeanY;
	double Dalpha = TMath::ATan2(S_M_Y,S_M_X);
	if((Dalpha-DDelta)<90*TMath::DegToRad() && (Dalpha-DDelta)>-90*TMath::DegToRad())
		DTail = 1;
	else
		DTail = -1;
	//get tail by using npixes at the both sides of centroid
	int positive=0;
	int negative=0;
	for(int i=0;i<sipm_sdp_x.size();i++)
	{
		if(sipm_sdp_x.at(i)>0)
			positive++;
		if(sipm_sdp_x.at(i)<0)
			negative++;
	}
	//DTail_pix = positive>=negative ? 1 : -1;
	if(positive>negative)
		DTail_pix = 1;
	else if(positive<negative)
		DTail_pix = -1;
	else
		DTail_pix = DTail;

	
//	printf("DTail:%lf DTail_pix:%lf DOutMeanX:%lf DMeanX:%lf DOutMeanX-DMeanX:%lf cos(DDelta):%lf sipm_sdp_x.size:%d positive:%d negative:%d\n",
//			DTail,DTail_pix, DOutMeanX,DMeanX, DOutMeanX-DMeanX,cos(DDelta),sipm_sdp_x.size(),positive,negative);
	//get two directions in this plane
	double tel_azi, tel_zen;
	WFCTAMap::Instance()->GetTelAZ(MainTel, tel_azi, tel_zen);
	slaDtp2s( DMeanX-DTail*DLength*cos(DDelta)/2, DMeanY-DTail*DLength*sin(DDelta)/2,
			tel_azi*TMath::DegToRad(), (90-tel_zen)*TMath::DegToRad(),
			&sdp_azi1, &sdp_zen1);
	slaDtp2s( DMeanX+DTail*DLength*cos(DDelta)/2, DMeanY+DTail*DLength*sin(DDelta)/2,
			tel_azi*TMath::DegToRad(), (90-tel_zen)*TMath::DegToRad(),
			&sdp_azi2, &sdp_zen2);

	wcda_sdp_azi1 = sdp_azi1 - 29.36*TMath::DegToRad();
	wcda_sdp_zen1 = 90*TMath::DegToRad()-sdp_zen1;
	wcda_sdp_azi2 = sdp_azi2 - 29.36*TMath::DegToRad();
	wcda_sdp_zen2 = 90*TMath::DegToRad()-sdp_zen2;

	km2a_sdp_azi1 = sdp_azi1;
	km2a_sdp_zen1 = 90*TMath::DegToRad()-sdp_zen1;
	km2a_sdp_azi2 = sdp_azi2;
	km2a_sdp_zen2 = 90*TMath::DegToRad()-sdp_zen2;

	double wcda_sdp_x, wcda_sdp_y, wcda_sdp_z;
	double km2a_sdp_x, km2a_sdp_y, km2a_sdp_z;
	//calc plane normal vector
	GetPlaneNormal(wcda_sdp_zen1, wcda_sdp_azi1, wcda_sdp_zen2, wcda_sdp_azi2, &wcda_sdp_x, &wcda_sdp_y, &wcda_sdp_z);
	GetPlaneNormal(km2a_sdp_zen1, km2a_sdp_azi1, km2a_sdp_zen2, km2a_sdp_azi2, &km2a_sdp_x, &km2a_sdp_y, &km2a_sdp_z);

	SDP_GLine_X_wcda = wcda_sdp_y;
	SDP_GLine_Y_wcda = -wcda_sdp_x;
	SDP_GLine_X_km2a = km2a_sdp_y;
	SDP_GLine_Y_km2a = -km2a_sdp_x;
}

void WFCTARec::CalcCoordsAloneSDP(const char* clean)
{
	sipm_sdp_isipm.clear();
	sipm_sdp_x.clear();
	sipm_sdp_pe.clear();
	sipm_sdp_t.clear();
	slice_sdp_x.clear();
	slice_sdp_pe.clear();
	slice_sdp_t.clear();
	slice_sdp_pix.clear();
	if(DDelta!=-1000)
	{
		//std::vector<double> sipm_sdp_x;
		//std::vector<double> sipm_sdp_pe;
		//find coordinate and pe on sdp of each sipms
		for(int ii=0;ii<iSipm.size();ii++)
		{
			//if(0==strcmp(clean,"clean"))
			//if(groupclean.at(ii)<group_cut) continue;
			if(0==groupclean.at(ii)) continue;
			double x, y;
			int SIPM = iSipm.at(ii);
			WFCTAMap::Instance()->GetwSipmXY(SIPM, x, y);
			double coords = GetCoordsOnLine(x, y, DMeanX, DMeanY, DDelta);
			sipm_sdp_isipm.push_back(SIPM);
			sipm_sdp_x.push_back(coords);
			sipm_sdp_pe.push_back(SipmPe.at(ii));
			sipm_sdp_t.push_back(SipmT.at(ii));
		}
		//find up and down limite
		double coords_up=-10000;
		double coords_down=10000;
		for(int i=0;i<sipm_sdp_x.size();i++)
		{
			double x = sipm_sdp_x.at(i)*TMath::RadToDeg();
			coords_up = coords_up>x ? coords_up : x;
			coords_down = coords_down<x ? coords_down : x;
		}
		double coords_range = coords_up-coords_down;
		double coords_step = 0.5;
		int nbins = int(coords_range/coords_step)+1;
		if(nbins<0)	nbins=1;
		const int nBins = nbins;
		//printf("coords_down:%lf coords_up:%lf coords_step:%lf nBins:%d\n",coords_down, coords_up, coords_step, nBins);

		//fill slice coordinate and pe along sdp
		double coords_x[nBins]={0};
		for(int i=0;i<nBins;i++)
		{
			coords_x[i] = coords_down + coords_step*i + coords_step*0.5;
		}
		double coords_pe[nBins]={0};
		double coords_t[nBins]={0};
		int coords_pix[nBins]={0};
		for(int i=0;i<sipm_sdp_x.size();i++)
		{
			double x = sipm_sdp_x.at(i)*TMath::RadToDeg();
			int ibin = int((x-coords_down)/coords_step);
			double y = sipm_sdp_pe.at(i);
			double y_t = sipm_sdp_t.at(i);
			coords_pe[ibin] += y;
			coords_t[ibin] += y_t;
			coords_pix[ibin]++;
		}
		for(int i=0;i<nBins;i++)
		{
			slice_sdp_x.push_back(coords_x[i]);
			slice_sdp_pe.push_back(coords_pe[i]);
			slice_sdp_t.push_back(coords_t[i]/coords_pix[i]);
			slice_sdp_pix.push_back(coords_pix[i]);
			//printf("slice_sdp:%d %lf %lf\n",i,coords_x[i],coords_pe[i]);
		}
	}
}


//private universal func
double WFCTARec::GetSiPMDist(double azi0, double zen0, double azi1, double zen1)//input:rad    return:deg
{
	double m0,n0,l0;
	double m1,n1,l1;

	m0 = sin(zen0)*cos(azi0);
	n0 = sin(zen0)*sin(azi0);
	l0 = cos(zen0);
	m1 = sin(zen1)*cos(azi1);
	n1 = sin(zen1)*sin(azi1);
	l1 = cos(zen1);
	double Mod0=sqrt(m0*m0+n0*n0+l0*l0);
	double Mod1=sqrt(m1*m1+n1*n1+l1*l1);

	double delta = (m0*m1+n0*n1+l0*l1)/(Mod0*Mod1);
	if(delta>1) delta = 1;
	if(delta<-1) delta = -1;
	delta = acos(delta)*TMath::RadToDeg();
	//printf("delta::::::::%.3lf\n",delta);
	return delta;
	//return sqrt((azi0-azi1)*(azi0-azi1)+(zen0-zen1)*(zen0-zen1));
}

// get coordinate along a line, (line_x,liney) is the reference point one the line.
double WFCTARec::GetCoordsOnLine(double pos_x, double pos_y, double line_x, double line_y, double line_delta)
{
	double dX = pos_x - line_x;
	double dY = pos_y - line_y;
	double dAlpha = TMath::ATan2(dY,dX);
	double dTheta = dAlpha-line_delta;
	double coords = sqrt(dX*dX+dY*dY)*cos(dTheta);
	//printf("pos_x:%lf line_x:%lf dX:%.2lf pos_y:%lf line_y:%lf dY:%.2lf line_delta:%.2lf dAlpha:%.2lf dTheta:%.2lf coords:%.2lf\n",
	//		pos_x,line_x,dX, pos_y,line_y,dY, line_delta*57.3,dAlpha*57.3,dTheta*57.3, coords);
	return coords;
}


void WFCTARec::GetEventMapOnFocus(int main_tel)
{
	WFCTAMap::Instance()->SetEventOnFocus(main_tel, iSipm);
	/*
	   wfctamap.w_ImageX.resize(wfctamap.sipms,0);
	   wfctamap.w_ImageY.resize(wfctamap.sipms,0);
	   int SIPM, jj;
	   double w_imagex, w_imagey;
	   for(int ii=0;ii<iSipm.size();ii++){
	   SIPM = iSipm.at(ii);
	   double azi = wfctamap.w_ImageAzi[SIPM];
	   double zen = wfctamap.w_ImageZen[SIPM];
	   slaDs2tp(azi, 90*TMath::DegToRad()-zen, wfctamap.TelA[main_tel-1]*TMath::DegToRad(), (90-wfctamap.TelZ[main_tel-1])*TMath::DegToRad(), &w_imagex, &w_imagey, &jj);
	   wfctamap.w_ImageX[SIPM] = w_imagex;
	   wfctamap.w_ImageY[SIPM] = w_imagey;
	   }
	   */
}

void WFCTARec::CalcHillasOnFocus(int main_tel)
{
	double x,y ,npe;
	DMeanX = 0;
	DMeanY = 0;
	for(int ii=0;ii<iSipm.size();ii++)
	{
		int SIPM = iSipm.at(ii);
		//if(groupclean.at(ii)<group_cut) continue;
		if(0==groupclean.at(ii)) continue;
		double x, y;
		WFCTAMap::Instance()->GetwSipmXY(SIPM, x, y);
		npe = SipmPe.at(ii);
		DMeanX += npe*x;
		DMeanY += npe*y;
	}
	if(0==DSize) return ;
	DMeanX /= DSize;
	DMeanY /= DSize;
	//printf("(MeanX,MeanY): (%5.2lf,%5.2lf)\n",DMeanX*TMath::RadToDeg(),DMeanY*TMath::RadToDeg());
	double tel_azi, tel_zen;
	WFCTAMap::Instance()->GetTelAZ(main_tel, tel_azi, tel_zen);
	slaDtp2s( DMeanX, DMeanY, tel_azi*TMath::DegToRad(), (90-tel_zen)*TMath::DegToRad(), &DMeanAzi, &DMeanZen);
	DMeanZen = 90*TMath::DegToRad()-DMeanZen;

	//calc length, width ..
	Double_t corrxx=0;
	Double_t corrxy=0;
	Double_t corryy=0;
	double dx, dy;

	for(int ii=0;ii<iSipm.size();ii++)
	{
		int SIPM = iSipm.at(ii);
		//if(groupclean.at(ii)<group_cut) continue;
		if(0==groupclean.at(ii)) continue;
		double x, y;
		WFCTAMap::Instance()->GetwSipmXY(SIPM, x, y);
		npe = SipmPe.at(ii);

		dx = x - DMeanX;
		dy = y - DMeanY;
		corrxx += npe * dx*dx;
		corrxy += npe * dx*dy;
		corryy += npe * dy*dy;
	}
	if (corrxy==0) return ;
	if(1==event_type) {
		return ;
	}

	Double_t d0    = corryy - corrxx;
	Double_t d1    = corrxy*2;
	Double_t d2    = d0 + TMath::Sqrt(d0*d0 + d1*d1);
	Double_t tand  = d2 / d1;
	Double_t tand2 = tand*tand;

	Double_t a, b;
	a = d2/d1;
	b = DMeanY - a*DMeanX;
	DDelta = TMath::ATan(a);

	Double_t s2 = tand2+1;
	Double_t s  = TMath::Sqrt(s2);
	Double_t axis1 = (tand2*corryy + d2 + corrxx)/s2/DSize;
	Double_t axis2 = (tand2*corrxx - d2 + corryy)/s2/DSize;

	DLength = axis1<0 ? 0 : TMath::Sqrt(axis1);
	DWidth  = axis2<0 ? 0 : TMath::Sqrt(axis2);
	//DAngleDist = (DSourceX-DMeanX)*(DSourceX-DMeanX)+(DSourceY-DMeanY)*(DSourceY-DMeanY);

	//DAlpha = asin(DAlpha);
}

void WFCTARec::OuterCentroid()
{
	vector<pair<int, double> > tubes;
	for(int ii=0;ii<iSipm.size();ii++)
	{
		//if(groupclean.at(ii)<group_cut) continue;
		if(0==groupclean.at(ii)) continue;
		int SIPM = iSipm.at(ii);
		double signal = SipmPe.at(ii);
		tubes.push_back(pair<int,double>(SIPM,signal));
	}
	sort(tubes.begin(),tubes.end(),cmp);

	int num_tubes = tubes.size();
	int num_out_tubes = round(num_tubes*0.8);
	double x,y ,npe;
	double OutSize = 0;
	DOutMeanX = 0;
	DOutMeanY = 0;
	for(int ii=0;ii<num_out_tubes;ii++)
	{
		int SIPM = tubes.at(ii).first;
		npe = tubes.at(ii).second;
		double x, y;
		WFCTAMap::Instance()->GetwSipmXY(SIPM, x, y);
		OutSize += npe;
		DOutMeanX += npe*x;
		DOutMeanY += npe*y;
	}
	if(0==OutSize)
	{
		return ;
	}
	//printf("OutSize:%lf DOutMeanX:%lf DOutMeanY:%lf\n",OutSize,DOutMeanX,DOutMeanY);
	DOutMeanX /= OutSize;
	DOutMeanY /= OutSize;
}

void WFCTARec::GetCoreToFocus(int main_tel, double zen,double azi,double corex,double corey, double *x, double *y)
{
	double Tel_x, Tel_y;
	WFCTAMap::Instance()->GetTelXY(main_tel, Tel_x, Tel_y);
	double Tel_X = Tel_x;
	double Tel_Y = Tel_y;
//	double Tel_Y = -(Tel_x*100 + 8600);
//	double Tel_X = Tel_y*100 + 18000;
	double k,m, n, l;
	int j;
	double focusx,focusy;
	m = sin(zen)*cos(azi);
	n = sin(zen)*sin(azi);
	l = cos(zen);
//	printf("mnl(%lf %lf %lf)\n",m,n,l);
//	l = -cos(zen);
	double x1, y1, z1 ;
	double m0, n0, l0, norm;

	z1 = 1000; //100000;
	k = z1/l;
	x1 = k*m + (corex);
	y1 = k*n + (corey);

	double WFCTAX=-18000;
	double WFCTAY=8600;
	m0 = -(Tel_X-x1);
	n0 = -(Tel_Y-y1);
//	m0 = (Tel_X+WFCTAX-x1);
//	n0 = (Tel_Y+WFCTAY-y1);
	l0 = -(0-z1);
	norm = sqrt(m0*m0+n0*n0+l0*l0);
	m0 = m0/norm;
	n0 = n0/norm;
	l0 = l0/norm;

	double theta1 =  acos(l0);
	double phi1 = acos(m0/sin(theta1));
	if(n0>0) phi1 = phi1;
	if(n0<0) phi1 = -phi1;

	double tel_azi, tel_zen;
	WFCTAMap::Instance()->GetTelAZ(main_tel, tel_azi, tel_zen);
	double T_A = tel_azi;
	double T_Z = tel_zen;
//	double T_A = tel_azi-90;
//	double T_Z = tel_zen;
//	printf("telxy(%lf, %lf) telaz(%lf(%lf), %lf)  core(%lf %lf)\n",Tel_X, Tel_Y, tel_azi, T_A, T_Z, corex, corey);
	slaDs2tp(phi1, 90*TMath::DegToRad() - theta1, (T_A)*TMath::DegToRad(), (90-T_Z)*TMath::DegToRad(),
			&focusx, &focusy, &j );
//	slaDs2tp(phi1, 90*TMath::DegToRad() - theta1, (tel_azi)*TMath::DegToRad(), (90-tel_zen)*TMath::DegToRad(),
//			&focusx, &focusy, &j );
//	focusx *= 2870; //mm
//	focusy *= 2870;
	*x = focusx;
	*y = focusy;
}

void WFCTARec::slaDs2tp ( double ra, double dec, double raz, double decz, double *xi, double *eta, int *j )
{
	double TINY=1e-6;
	double sdecz, sdec, cdecz, cdec, radif, sradif, cradif, denom;

	/* Trig functions */
	sdecz = sin ( decz );
	sdec = sin ( dec );
	cdecz = cos ( decz );
	cdec = cos ( dec );
	radif = ra - raz;
	sradif = sin ( radif );
	cradif = cos ( radif );
	/* Reciprocal of star vector length to tangent plane */
	denom = sdec * sdecz + cdec * cdecz * cradif;
	/* Handle vectors too far from axis */
	if ( denom > TINY )         {   *j = 0; }
	else if ( denom >= 0.0 )    {   *j = 1; denom = TINY;}
	else if ( denom > -TINY )   {   *j = 2; denom = -TINY;}
	else                        {   *j = 3; }
	/* Compute tangent plane coordinates (even in dubious cases) */
	*xi = cdec * sradif / denom;
	*eta = ( sdec * cdecz - cdec * sdecz * cradif ) / denom;
}

void WFCTARec::slaDtp2s ( double xi, double eta, double raz, double decz, double *ra, double *dec )
{
	double d2pi = 2*TMath::Pi();
	double sdecz, cdecz, denom;
	double w;

	sdecz = sin ( decz );
	cdecz = cos ( decz );
	denom = cdecz - eta * sdecz;
	w = dmod ( atan2 ( xi, denom ) + raz, d2pi );
	*ra = ( w >= 0.0 ) ? w : w + d2pi;
	*dec = atan2 ( sdecz + eta * cdecz, sqrt ( xi * xi + denom * denom ) );
}

void WFCTARec::GetPlaneNormal(double zenith1,double azimuth1,double zenith2,double azimuth2,double *x, double *y, double *z)
{
	double m, n,l,   m0,n0,l0;
	double x1, y1, z1,norm;

	m = sin(zenith1)*cos(azimuth1);
	n = sin(zenith1)*sin(azimuth1);
	l = cos(zenith1);

	m0 = sin(zenith2)*cos(azimuth2);
	n0 = sin(zenith2)*sin(azimuth2);
	l0 = cos(zenith2);

	x1 = n*l0 - n0*l;
	y1 = -(m*l0 - m0*l);
	z1 = m*n0 - m0*n;

	norm = sqrt(x1*x1+y1*y1+z1*z1);
	*x = x1/norm;
	*y = y1/norm;
	*z = z1/norm;
}

void WFCTARec::GetPe_Size(double* tel_pe_size) const 
{
	for(int i=0;i<6;i++)
	{
		tel_pe_size[i] = pe_size[i];
	}
}

void WFCTARec::GetNpix_Size(int* tel_npix_size) const 
{
	for(int i=0;i<6;i++)
	{
		tel_npix_size[i] = npix_size[i];
	}
}

void WFCTARec::GetInter_Npix(int* tel_inter_npix) const 
{
	for(int i=0;i<6;i++)
	{
		tel_inter_npix[i] = inter_npix[i];
	}
}

void WFCTARec::GetCleanImage(std::vector<int>& clean_sipm, std::vector<double>& clean_pe, std::vector<double>& clean_t) const
{
	for(int ii=0;ii<iSipm.size();ii++)
	{
		//if(groupclean.at(ii)<group_cut) continue;
		if(0==groupclean.at(ii)) continue;
		clean_sipm.push_back(iSipm.at(ii));
		clean_pe.push_back(SipmPe.at(ii));
		clean_t.push_back(SipmT.at(ii));
	}
}

void WFCTARec::GetCleanImageOnFocus(int main_tel, std::vector<int>& clean_sipm, std::vector<double>& clean_pe, std::vector<double>& clean_x, std::vector<double>& clean_y)
{
	int SIPM, jj;
	double w_imagex, w_imagey;
	double tel_azi, tel_zen;
	WFCTAMap::Instance()->GetTelAZ(main_tel, tel_azi, tel_zen);
	for(int ii=0;ii<iSipm.size();ii++)
	{
		//if(groupclean.at(ii)<group_cut) continue;
		if(0==groupclean.at(ii)) continue;
		SIPM = iSipm.at(ii);
		double azi, zen;
		if(merge_field)
			WFCTAMap::Instance()->GetwSipmAZ(SIPM, azi, zen);
		else
			WFCTAMap::Instance()->GetSipmAZ(SIPM, azi, zen);
		slaDs2tp(azi, 90*TMath::DegToRad()-zen, tel_azi*TMath::DegToRad(), (90-tel_zen)*TMath::DegToRad(), 
				&w_imagex, &w_imagey, &jj);
		clean_sipm.push_back(SIPM);
		clean_pe.push_back(SipmPe.at(ii));
		clean_x.push_back(w_imagex);
		clean_y.push_back(w_imagey);
		//printf("wfcta focous %lf %lf | %lf %lf\n",w_imagex,w_imagey,w_imagex*57.3,w_imagey*57.3);
	}
}

void WFCTARec::GetCoordsOnSDP(std::vector<int>& sipm_isipm, std::vector<double>& sipm_pe, std::vector<double>& sipm_t, std::vector<double>& sipm_coords, std::vector<double>& slice_pe, std::vector<double>& slice_t, std::vector<int>& slice_pix, std::vector<double>& slice_coords) const
{
	for(int ii=0;ii<sipm_sdp_x.size();ii++)
	{
		sipm_isipm.push_back(sipm_sdp_isipm.at(ii));
		sipm_pe.push_back(sipm_sdp_pe.at(ii)*TMath::RadToDeg());
		sipm_t.push_back(sipm_sdp_t.at(ii));
		sipm_coords.push_back(sipm_sdp_x.at(ii)*TMath::RadToDeg());
	}
	for(int ii=0;ii<slice_sdp_x.size();ii++)
	{
		slice_pe.push_back(slice_sdp_pe.at(ii));
		slice_t.push_back(slice_sdp_t.at(ii));
		slice_pix.push_back(slice_sdp_pix.at(ii));
		slice_coords.push_back(slice_sdp_x.at(ii));
	}
}

int WFCTARec::GetSipmIdInSimulation(int sipm_id)
{
	int itel = sipm_id/1024 + 1;
	int isipm_in_tel = sipm_id % 1024;
	int icolum = isipm_in_tel % 32;
	int iline = isipm_in_tel / 32;
	return (31-icolum) + 32*iline + (itel-1)*1024;

}



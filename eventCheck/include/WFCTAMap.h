#ifndef WFCTAMAP_H
#define WFCTAMAP_H

#include <vector>

#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B)):(A)+(B)*floor(-(A)/(B))):(A))
class WFCTAMap
{
	public:
		static WFCTAMap* Instance()
		{
			if (0==m_instance) m_instance = new WFCTAMap();
			return m_instance;
		}
		~WFCTAMap();
	private:
		WFCTAMap();
		static WFCTAMap* m_instance;

	public:
		int NTel() {
			return ntel;
		};
		int MaxTelId() {
			return maxTelId;
		};
		int NSipm() {
			return nsipm;
		};
		//telid: 1,2,3,4 ...
		//isipm: 0,1 2,3,4,5, ...
		static double GetMAXDIST() {
			return MAXDIST;
		};
		void GetTelId(int telid);
		void GetTelFocal(int telid, double& focal);
		void GetTelXY(int telid, double& x, double& y);
		void GetTelXYinWCDA1(int telid, double& x, double& y);
		void GetTelAZ(int telid, double& azi, double& zen);
		void GetSipmStat(int isipm, int& stat);
		void GetSipmStat(int telid, int isipm, int& stat);
		void GetSipmXY(int isipm, double& x, double& y);
		void GetSipmXY(int telid, int isipm, double& x, double& y);
		void GetSipmAZ(int isipm, double& azi, double& zen);
		void GetSipmAZ(int telid, int isipm, double& azi, double& zen);
		void GetwSipmXY(int isipm, double& x, double& y);
		void GetwSipmXY(int telid, int isipm, double& x, double& y);
		void GetwSipmAZ(int isipm, double& azi, double& zen);
		void GetwSipmAZ(int telid, int isipm, double& azi, double& zen);
		void SetEventOnFocus(int main_tel, std::vector<int> iSipm);

		void GetXRange(int iline, double& xMin, double& xMax);
		void GetYRange(double& yMin, double& yMax);

	private:
		static const double MAXDIST;
		int maxTelId;
		int ntel;
		int nsipm;
		std::vector<int> Tel_Id; //m
		std::vector<double> TelX; //m
		std::vector<double> TelY; //m
		std::vector<double> TelXinWCDA1; //m
		std::vector<double> TelYinWCDA1; //m
		std::vector<double> TelA; //deg
		std::vector<double> TelZ; //deg
		std::vector<double> TelFocal; //mm
		std::vector<int> SipmId;  // rad
		std::vector<int> SipmStat;  // rad
		std::vector<double> SipmA;  // rad
		std::vector<double> SipmZ;   // rad
		std::vector<double> SipmX;  // rad
		std::vector<double> SipmY;  // rad
		std::vector<int> SipmId_inArray;  // rad
		std::vector<double> w_SipmA;  // rad
		std::vector<double> w_SipmZ;   // rad
		std::vector<double> w_SipmX;  // rad
		std::vector<double> w_SipmY;  // rad
		std::vector<std::pair<int,int> > merged_pix;
		double x_max[32];
		double x_min[32];
		double y_max;
		double y_min;

		void ReadTelFile(const char* filename);
		void SetSiPMMap();
		void GetXYRange();
		void SetSiPMStat(const char* filename);
		void SetwSiPMMap();
		double GetSiPMDist(double azi0, double zen0, double azi1, double zen1);
		void slaDs2tp ( double ra, double dec, double raz, double decz, double *xi, double *eta, int *j );
		void slaDtp2s ( double xi, double eta, double raz, double decz, double *ra, double *dec );
};

///////////////////////
////return function////
///////////////////////
inline void WFCTAMap::GetTelFocal(int telid, double& focal)
{
	focal = TelFocal[telid-1];
}
inline void WFCTAMap::GetTelXY(int telid, double& x, double& y)
{
	x = TelX[telid-1];
	y = TelY[telid-1];
}
inline void WFCTAMap::GetTelXYinWCDA1(int telid, double& x, double& y)
{
	x = TelXinWCDA1[telid-1];
	y = TelYinWCDA1[telid-1];
}
inline void WFCTAMap::GetTelAZ(int telid, double& azi, double& zen)
{
	azi = TelA[telid-1];
	zen = TelZ[telid-1];
}
inline void WFCTAMap::GetSipmStat(int isipm, int& stat)
{
	stat = SipmStat[isipm];
}
inline void WFCTAMap::GetSipmStat(int telid, int isipm, int& stat)
{
	int Isipm = (telid-1)*1024 + isipm;
	stat = SipmStat[Isipm];
}
inline void WFCTAMap::GetSipmXY(int isipm, double& x, double& y)
{
	x = SipmX[isipm];
	y = SipmY[isipm];
}
inline void WFCTAMap::GetSipmXY(int telid, int isipm, double& x, double& y)
{
	int Isipm = (telid-1)*1024 + isipm;
	x = SipmX[Isipm];
	y = SipmY[Isipm];
}
inline void WFCTAMap::GetSipmAZ(int isipm, double& azi, double& zen)
{
	azi = SipmA[isipm];
	zen = SipmZ[isipm];
}
inline void WFCTAMap::GetSipmAZ(int telid, int isipm, double& azi, double& zen)
{
	int Isipm = (telid-1)*1024 + isipm;
	azi = SipmA[Isipm];
	zen = SipmZ[Isipm];
}
inline void WFCTAMap::GetwSipmXY(int isipm, double& x, double& y)
{
	x = w_SipmX[isipm];
	y = w_SipmY[isipm];
}
inline void WFCTAMap::GetwSipmXY(int telid, int isipm, double& x, double& y)
{
	int Isipm = (telid-1)*1024 + isipm;
	x = w_SipmX[Isipm];
	y = w_SipmY[Isipm];
}
inline void WFCTAMap::GetwSipmAZ(int isipm, double& azi, double& zen)
{
	azi = w_SipmA[isipm];
	zen = w_SipmZ[isipm];
}
inline void WFCTAMap::GetwSipmAZ(int telid, int isipm, double& azi, double& zen)
{
	int Isipm = (telid-1)*1024 + isipm;
	azi = w_SipmA[Isipm];
	zen = w_SipmZ[Isipm];
}

inline void WFCTAMap::GetXRange(int iline, double& xMin, double& xMax)
{
	xMin = x_min[iline];
	xMax = x_max[iline];
}
inline void WFCTAMap::GetYRange(double& yMin, double& yMax)
{
	yMin = y_min;
	yMax = y_max;
}

#endif // WFCTAMAP_H

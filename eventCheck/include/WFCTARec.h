#ifndef WFCTAREC_H
#define WFCTAREC_H
#include <TMath.h>
#include <vector>
#include "WFCTAMap.h"
#include "CollectIntensity.h"

extern bool cmp(const std::pair<int,double>& tube1, const std::pair<int,double>& tube2);

class WFCTARec
{
	private:
		//event
		int merge_field;
		int event_type;
		std::vector<int> iSipm;
		std::vector<double> SipmPe;
		std::vector<double> SipmT;

		//rec
		int group_cut;
		std::vector<int> timeclean;
		std::vector<int> groupclean;
		std::vector<double> slice_sdp_x;
		std::vector<double> slice_sdp_pe;
		std::vector<double> slice_sdp_t;
		std::vector<int> slice_sdp_pix;
		std::vector<int> sipm_sdp_isipm;
		std::vector<double> sipm_sdp_x;
		std::vector<double> sipm_sdp_pe;
		std::vector<double> sipm_sdp_t;

		int MaxPeTel;
		int MaxSipm;
		double MaxSipmPe;
		int MainTel;
		int DNpix;
		double DTail;
		double DTail_pix;
		double DSize;
		double DMeanZen; //rad
		double DMeanAzi; //rad
		double DMeanX; //on focus
		double DMeanY; //on focus
		double DOutMeanX; //on focus
		double DOutMeanY; //on focus
		double DLength; //on focus
		double DWidth; //on focus
		double DDelta; //on focus
		double DAlpha; //
		double SDP_GLine_X_wcda;  //sdp_ground_line in wcda map
		double SDP_GLine_Y_wcda;  //sdp_ground_line in wcda map
		double SDP_GLine_X_km2a;  //sdp_ground_line in km2a
		double SDP_GLine_Y_km2a;  //sdp_ground_line in km2a

		//compare
		double pe_size[6];
		int npix_size[6];
		int inter_npix[6];
		int comp_tel1;
		int comp_tel2;

		
	public:
		double GetSiPMDist(double azi0, double zen0, double azi1, double zen1);
		double GetCoordsOnLine(double pos_x, double pos_y, double line_x, double line_y, double line_delta);
		int GetSipmIdInSimulation(int sipm_id);
		void GetEventMapOnFocus(int main_tel);
		void CalcHillasOnFocus(int main_tel);
		void OuterCentroid();

		void GetCoreToFocus(int main_tel, double zen,double azi,double corex,double corey, double *x, double *y);
		void slaDs2tp ( double ra, double dec, double raz, double decz, double *xi, double *eta, int *j );
		void slaDtp2s ( double xi, double eta, double raz, double decz, double *ra, double *dec );
		void GetPlaneNormal(double zenith1,double azimuth1,double zenith2,double azimuth2,double *x, double *y, double *z);

		WFCTARec();
		~WFCTARec();
		//wfcta map relevant parameters
		//WFCTAMap wfctamap;

		//map
		void SetFieldType(int merge_f)
		{
			merge_field = merge_f;
		}

		//Init
		void WFCTAInit();
		//event
		void SetCompTel(int tel1, int tel2) {
			comp_tel1 = tel1;
			comp_tel2 = tel2;
		};
		void SetWFCTAEvent(const std::vector<int>& iisipm, const std::vector<double>& isipmpe, const std::vector<double>& isipmt);
		void SetEventType(int evt_type);

		//rec
		void TimeClean(double pecut=50);
		void GroupClean(int g_cut=3);
		void IslandClean();
		void CalcMainTel(int main_tel);
		void MergeEvent(int do_clean=1);
		void CalcHillas();
		void CalcSDP();
		void CalcCoordsAloneSDP(const char* clean="clean");

		//return
		int GetMaxPeTel()				const {return MaxPeTel;};
		int GetMaxSipm()				const {return MaxSipm;};
		double GetMaxSipmPe()			const {return MaxSipmPe;};
		int GetMainTel()				const {return MainTel;};
		int GetNpix()					const { return DNpix;};
		double GetTail()				const {return DTail;};
		double GetTail_pix()			const {return DTail_pix;};
		double GetSize()				const { return DSize;   };
		double GetMeanZen()				const { return DMeanZen;  };
		double GetMeanAzi()				const { return DMeanAzi;  };
		double GetMeanX()				const { return DMeanX;  };
		double GetMeanY()				const { return DMeanY;  };
		double GetOutMeanX()			const { return DOutMeanX;  };
		double GetOutMeanY()			const { return DOutMeanY;  };
		double GetLength()				const { return DLength; };
		double GetWidth()				const { return DWidth; };
		double GetDelta()				const { return DDelta; };
		double GetAlpha()				const { return DAlpha;  };
		double GetSDP_GLine_X_wcda()	const {	return SDP_GLine_X_wcda;	};
		double GetSDP_GLine_Y_wcda()	const {	return SDP_GLine_Y_wcda;	};
		double GetSDP_GLine_X_km2a()	const {	return SDP_GLine_X_km2a;	};
		double GetSDP_GLine_Y_km2a()	const {	return SDP_GLine_Y_km2a;	};
		void GetPe_Size(double* tel_pe_size) const ;
		void GetNpix_Size(int* tel_npix_size) const ;
		void GetInter_Npix(int* tel_inter_npix) const ;
		void GetCleanImage(std::vector<int>& clean_sipm, std::vector<double>& clean_pe, std::vector<double>& clean_t) const;
		void GetCleanImageOnFocus(int main_tel, std::vector<int>& clean_sipm, std::vector<double>& clean_pe, std::vector<double>& clean_x, std::vector<double>& clean_y);
		void GetCoordsOnSDP(std::vector<int>& sipm_isipm, std::vector<double>& sipm_pe, std::vector<double>& sipm_t, std::vector<double>& sipm_coords, std::vector<double>& slice_pe, std::vector<double>& slice_t, std::vector<int>& slice_pix, std::vector<double>& slice_coords) const;

};

#endif // WFCTAREC_H

#ifndef WFCTASLC_H
#define WFCTASLC_H

#include <vector>

class WFCTASLC
{
	public:
		WFCTASLC();
		~WFCTASLC();

	public:
		void Clear();
		void ReadSlcFile(int year, int month, int day);
		void ReadSlcFileYMJ(int year, int month, int day);
		double GetDaq_Stu(long long rbTime, int tel_id);
		double GetSlow_Control_Stu(long long rbTime, int tel_id);
		double GetTemperature_Stu(long long rbTime, int tel_id);
		double GetCloud_Stu(long long rbTime, int tel_id);
		double GetDoor_Angle(long long rbTime, int tel_id);
		int GetR_LED(long long rbTime, int tel_id);
		int GetD_LED(long long rbTime, int tel_id);


	private:
		struct SlcStatus{
			double daq_stu;
			double slow_control_stu;
			double temperature_stu;
			double cloud_stu;
			double door_angle;
			int r_led;
			int d_led;
		};

		std::vector<std::pair<double, SlcStatus> > slc_status[20];
		void m_ReadSlcFile(int year, int month, int day, int tel);
		void m_ReadSlcFileYMJ(int year, int month, int day, int tel);
		void GetSlcStatus(long long rbTime, int tel_id, SlcStatus& slcstatus);

};

///////////////////////
////return function////
///////////////////////


#endif // WFCTASLC_H

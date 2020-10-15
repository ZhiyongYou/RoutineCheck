#ifndef COLLECTINTENSITY_H
#define COLLECTINTENSITY_H

#include <vector>

class CollectIntensity
{
	public:
		static CollectIntensity* Instance()
		{
			if (0==m_instance) m_instance = new CollectIntensity();
			return m_instance;
		}
		~CollectIntensity();
	private:
		CollectIntensity();
		static CollectIntensity* m_instance;

	public:
		void GetIntensity(int isipm, double& intens);

	private:
		std::vector<double> intensity;

		int ReadIntensityFile(const char* filename);
};

///////////////////////
////return function////
///////////////////////
inline void CollectIntensity::GetIntensity(int isipm, double& intens)
{
	int Isipm = isipm % 1024;
	intens = intensity[Isipm];
}

#endif // COLLECTINTENSITY_H

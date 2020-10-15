#ifndef READFILE_H
#define READFILE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

class ReadFile
{
	protected:
		double conf_info[16];
		std::string filename;
		FILE *infp;
		inline void SetFileName(const char* f_name){filename = f_name;};
		void mreadInfo(double *conf);

	public:
		ReadFile();
		~ReadFile();

		bool OpenFile();
		bool OpenFile(const char* f_name);
		inline void CloseFile(){fclose(infp);};
		void ReadInfo(double *conf, const char* conf_type);
		void GetNoisePix(int iTel, std::vector<int> &noisSiPM, std::vector<int> &edge_noisSiPM);

};

#endif // READFILE_H

#pragma once     // to ensure header file is included only once
#ifndef _SIF_H_  // if weather header is not defined, define weather.h
#define _SIF_H_

struct TSif			// container for Sif param 
{
public:
	TSif()
	{
		SIF_GasEx = 1, SifData = 0.0, SifData = 0.0, qlSunlit = 0.0, qlShaded = 0.0, KdfSunlit = 0.0,
			KdfShaded = 0.0, PhiPS2Sunlit = 0.0, PhiPS2Shaded = 0.0, eSIF = 0.0;
	}
	int SIF_GasEx;
	double SifData;
	double qlSunlit;
	double qlShaded; //  
	double KdfSunlit; // 
	double KdfShaded;
	double PhiPS2Sunlit;
	double PhiPS2Shaded; // 
	double eSIF;
};
#endif


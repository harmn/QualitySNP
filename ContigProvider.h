#ifndef __CONTIGPROVIDER_H__
#define __CONTIGPROVIDER_H__

#include "ACEFile.h"

class ContigProvider
{
public:
	ContigProvider(void);
	~ContigProvider(void);

	Contig* nextContig();
	bool init();

private:
	ContigFile* _pContigFile;
};

#endif

#ifndef __CSVWRITER_H__
#define __CSVWRITER_H__

#include <fstream>
#include <vector>
#include "Configuration.h"
#include "Logger.h"
#include "Contig.h"
#include "SeqRead.h"

using namespace std;

class CSVWriter
{
public:
	CSVWriter(void);
	~CSVWriter(void);

	enum Outputs {
		CONTIGSFILE,
		READSFILE,
		HAPLOTYPESFILE,
		VARIATIONSFILE
	};

	bool init();
	bool writeContig(Contig* pContig);

private:
	bool openOutputFile(const string& output);

	vector<ofstream*>	_rgOutputs;
	bool				_bFail;
	vector<string>		_rgOutPutTypes;
	vector<string>		_rgColumnNames;
};

#endif

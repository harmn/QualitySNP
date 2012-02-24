#ifndef __ACEFILE_H__
#define __ACEFILE_H__

#include <string>
#include <fstream>
#include "ContigFile.h"

static const int MAX_LINE_LENGTH = 4000;

typedef enum BlockLabel {
	UnknownTag,
	CO,
	RD,
	BQ,
	AF,
	QA,
} BLOCKLABEL;

using namespace std;

class ACEFile : public ContigFile
{
public:
	ACEFile(const string&);
	~ACEFile(void);

	Contig* nextContig();
	bool openFile();
	bool isValid();

private:
	bool emptyLine(const string&);
	Contig* parseContigLine(const string&);
	SeqRead* parseReadLine(const string&, Contig*);
	SeqRead* parseReadLocationLine(const string&, Contig*);
	void parseQualitySegmentLine(const string&, SeqRead*);
	void stripTrailingFields(string&);
	
private:
	string		_contigFilename;
	ifstream	_ifACEFile;
	Contig*		_pContig;

	BLOCKLABEL getBlockLabel(const string&);
};

#endif

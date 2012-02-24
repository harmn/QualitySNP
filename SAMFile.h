#ifndef __SAMFILE_H__
#define __SAMFILE_H__

#include <string>
#include <fstream>
#include <vector>
#include <list>
#include "SAMRead.h"
#include "ContigFile.h"

class SAMFile : public ContigFile
{
public:
	SAMFile(const string& name);
	~SAMFile(void);

	Contig* nextContig();
	bool openFile();
	bool isValid();

private:
	bool parseCigar(const string&, list<operation>&);
	bool emptyLine(const string&);
	SAMRead* parseReadLine(const string&);

private:
	string		_contigFilename;
	ifstream	_ifSAMFile;
	SAMRead*	_pRead;
};

#endif

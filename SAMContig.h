#ifndef __SAMCONTIG_H__
#define __SAMCONTIG_H__

#include <string>
#include <list>
#include "ContigFile.h"

class SAMRead;

using namespace std;

class SAMContig
{
public:
	SAMContig();
	~SAMContig();

	void setName(const string& name) { _name = name; }
	void setSequence(const string& sequence) { _sequence = sequence; }
	void addRead(SAMRead* pRead) { _reads.push_back(pRead); }
	const string& getName() { return _name; }
	bool constructReferenceSequence();
	void stretchReads();
	void mergeReadPairs();
	Contig* toContig();

private:
	string					_name;
	string					_sequence;
	list<SAMRead*>			_reads;
};

#endif

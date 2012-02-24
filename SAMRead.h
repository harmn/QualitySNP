#ifndef __SAMREAD_H__
#define __SAMREAD_H__

#include <string>
#include "SeqRead.h"

using namespace std;

struct operation {
	char op;
	int size;
};

class SAMRead
{
public:
	SAMRead(const string& name): _name(name)		{}
	~SAMRead();
	const string getName()							{ return _name; }
	void setSequence(const string& sequence);
	const string& getSequence()						{ return _sequence; }
	void setQuality(const string& quality);
	const string& getQuality()						{ return _quality; }
	void setContigName(const string& contigName)	{ _contigName = contigName; }
	const string& getContigName()					{ return _contigName; }
	list<operation>& getOperations()				{ return _ops; }
	void setStartPosition(int startPosition)		{ _startPosition = startPosition - 1; }
	int getStartPosition()							{ return _startPosition; }
	char getNucleotideAt(unsigned int pos);
	char getOperationAt(unsigned int pos);
	void insertGapAt(unsigned int pos);
	bool processOperations();
	int getLastPosition()							{ return _startPosition + _sequence.size() + 1; }
	void merge(SAMRead*);
	SeqRead* toSeqRead();

private:
	string				_name;
	string				_contigName;
	string				_sequence;
	string				_opsequence;
	string				_quality;
	list<operation>		_ops;
	int					_startPosition;
};

#endif

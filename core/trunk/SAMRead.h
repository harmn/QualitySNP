/*
* This File is part of QualitySNP; a program to detect Single Nucleotide Variations
* https://trac.nbic.nl/qualitysnp/
*
*   Copyright (C) 2012 Harm Nijveen
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
    void setGroup(const string& group)              { _group = group; }
    string& getGroup()                              { return _group; }
	int getStartPosition()							{ return _startPosition; }
    bool hasInsertions()                            { return _bHasInsertions; }
	char getNucleotideAt(unsigned int pos);
	char getOperationAt(unsigned int pos);
    bool insertGapAt(unsigned int pos);
	bool processOperations();
	int getLastPosition()							{ return _startPosition + _sequence.size() + 1; }
	void setMapQuality(int mapQuality)				{ _mapQuality = mapQuality; }
	int getMapQuality()								{ return _mapQuality; }
	void merge(SAMRead*);
	SeqRead* toSeqRead();

private:
	string				_name;
	string				_contigName;
	string				_sequence;
	string				_opsequence;
	string				_quality;
    string              _group;
	list<operation>		_ops;
	int					_startPosition;
	int					_mapQuality;
    bool                _bHasInsertions;
};

#endif

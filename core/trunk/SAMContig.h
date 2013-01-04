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

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

#ifndef __SEQREAD_H__
#define __SEQREAD_H__

#include <string>
#include <vector>
#include <sstream>

class Contig;
class HaploType;

using namespace std;

class SeqRead
{
public:
	SeqRead(string, Contig* = NULL);
	~SeqRead(void);
	void setSequence(const string&);
    void setGroup(const string& group)      { _group = group; }
    const string& getName()         { return _name; }
    const string& getSequence()     { return _sequence; }
    const string& getGroup()         { return _group; }

    int getClippedSequenceLength()  { return _qualClipEnd - _qualClipStart + 1; }
    int getClippedsequenceStart()   { return _qualClipStart; }
    int getClippedsequenceEnd()     { return _qualClipEnd; }
    int getSNPCount(bool bHQ)       { return bHQ ? _cHQSNP : _cSNP; }

    int getStartPosition() { return _startPosition; }
    int getEndPosition(); 
	void setStartPosition(int startPosition) { _startPosition = startPosition - 1; }
	void setQualClip(unsigned int, unsigned int);
	void setSNPCount(int cSNP) { _cSNP = cSNP; } 
	void setContig(Contig* pContig) { _parent = pContig; } 
	void setQualitySanger(const string&);

    const string toString();
	char getNucleotideAt(int pos);
	bool isHighQuality(int pos);
	bool calculateSNPCount();

    void setHaploType(HaploType* pHaploType) {_pHaploType = pHaploType; }
	HaploType* getHaploType() { return _pHaploType; }

    string toCSV();

    const string& getGroupFromName(const string &separator);

private:
	void			setLowQualityBounds();

	Contig*			_parent;
	const string	_name;
	string			_sequence;
	int				_startPosition;
    int             _endPosition;
	int				_qualClipStart;
	int				_qualClipEnd;
	int				_HQstart;
	int				_HQstop;
	int				_lowQual5p;
	int				_lowQual3p;
    string          _group;
	int				_cSNP;
	int				_cHQSNP;
	HaploType*		_pHaploType;
	vector<int>		_quality;
};

#endif

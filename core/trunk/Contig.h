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

#ifndef __CONTIG_H__
#define __CONTIG_H__

#include <vector>
#include <list>
#include <string>
#include <map>
#include "Logger.h"

using namespace std;

class Variation;
class HaploType;
class SeqRead;
class ReadGroup;

class Contig
{
public:
	Contig(const string&);
	~Contig(void);

	//getters/setters
	const list<SeqRead*>& getReads() { return _reads; }
	const vector<Variation*>& getVariations() { return _variations; }
	const list<HaploType*>& getHaploTypes() { return _haploTypes; }
	const string& getName() { return _name; }
	const string& getSequence() {	return _sequence;	}
	double getDvalue();
	void setSequence(const string&);
	void setQuality(const string&);
	bool setQualityAt(unsigned int pos, int quality);
	char getSequenceAt(unsigned int);
	int getQualityAt(unsigned int);
	int getReadIndex(const string&);
    HaploType* getHaploTypeForID(int id);

	// get metrics
	int getSequenceLength() { return _sequence.length(); }
	int getReadCount() { return _cReads; }
	int getPotentialSNPCount();
	int getHighConfidenceSNPCount();
	int getReliableSNPCount();
	int getHaploTypeCount();
	int getMaxHaploTypePerSNPCount();

	SeqRead* findRead(const string&);
	bool addRead(SeqRead*);
    void addVariation(Variation*);
    void addHaploType(HaploType*);
	void calculateProperties();
	void sortReads();
	int findMarkerSNPs();

	string reads2CSV();
	string haploTypes2CSV();
	string toCSV();
	string haploTypeReadLinks2CSV();
	string variations2CSV();
    string readGroups2CSV(vector<string> readGroupNames);
	int maskHomopolymericTracts(int limit);
    string getSequenceIUPAC();

private:
	void determineVariations();
    void determineHaploTypes();
    bool addReadToHaploType(SeqRead*);
	void calculateSNPCountPerRead();
	bool isHighQuality(unsigned int pos);
    void determineReliableSNPs();
    void setReadGroupNames();
    void makeReadGroups();

	const string			_name;
	string					_sequence;
	map<string, SeqRead*>	_readMap;
	list<SeqRead*>			_reads;
    map<string, ReadGroup*> _readGroups;
    //list<ReadGroup*>        _readGroups;
	double					_Dvalue;
	int						_cReads;
	list<HaploType*>		_haploTypes;
	vector<Variation*>		_variations;
	vector<int>				_quality;
	int						_cPotentialSNP;
	int						_cHighConfidenceSNP;
	int						_cReliableSNP;
    int					    _defaultQualityScore;
};

#endif

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

#ifndef __HAPLOTYPE_H__
#define __HAPLOTYPE_H__
#include <vector>
#include <list>
#include <string>
#include <map>

class Contig;
class SeqRead;

using namespace std;

class HaploType
{
public:
	HaploType(Contig* pContig, int id);
	~HaploType(void);
    bool tryAddRead(SeqRead*);
	
	int getReadCount() { return _cReads; }
	const string toString();
	const list<SeqRead*>& getReads() { return _reads; }
    int getDefiningSNPCount() { return _definingSNP.size(); }
	void addDefiningSNP(int pos) { _definingSNP.push_back(pos); }
	int getID() { return _id; }
    void setID(int id) { _id = id; }
	const string toCSV();
	bool modified() { return _bModified; }
	void setModified(bool bModified) { _bModified = bModified; }
    void addRead(SeqRead*);
    char getNucleotideAt(int pos);
    int getLastVariablePosition();

private:
	int	matchAt(int, char);
    void matchRead(SeqRead* pRead, int& match, int& misMatch, map<int,int>& varNuc);
	void addRead(SeqRead*,const map<int,int>&);

private:
	Contig*					_pContig;
	list<SeqRead*>			_reads;
	int						_cReads;
	vector<int>				_definingSNP;
	map<int,vector<int> >	_varNuc;
	double					_singleSNPThreshold;
	double					_allSNPThreshold;
	int						_id;
	bool					_bModified;
    int                     _lastVariablePosition;
};

#endif

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
    int  readCount();
    int contigCount();
    vector<string>& readGroups();

private:
	bool emptyLine(const string&);
	Contig* parseContigLine(const string&);
	SeqRead* parseReadLine(const string&, Contig*);
	SeqRead* parseReadLocationLine(const string&, Contig*);
	void parseQualitySegmentLine(const string&, SeqRead*);
	void stripTrailingFields(string&);
    bool collectStatistics();
	
private:
    string          _contigFilename;
    ifstream        _ifACEFile;
    Contig*         _pContig;
    int             _cReads;
    int             _cContigs;
    vector<string>  _readGroups;

	BLOCKLABEL getBlockLabel(const string&);
};

#endif

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

#ifndef __CONTIGPRINTER_H__
#define __CONTIGPRINTER_H__

class ContigPrinter
{
public:
	ContigPrinter(Contig*);
	~ContigPrinter(void);
	string toString();

private:
	void printSummaryLine(stringstream& result);
	void printSNPs(stringstream& result, bool bHQOnly);
	void printHaploTypes(stringstream& result);
	void printMarkerSNPs(stringstream& result, unsigned int flankSize);
	void printAlignment(stringstream& result);
	Contig* _pContig;
};

#endif

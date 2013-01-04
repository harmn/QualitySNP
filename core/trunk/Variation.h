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

#ifndef __VARIATION_H__
#define __VARIATION_H__

#include "HaploType.h"

class Contig;

class Variation
{
public:
	Variation(Contig*, int);
	~Variation(void);
    int getConfidenceScore();
	void setQualityNucCount(int*, int*);
	bool setNucCount(int*);
	int getAlleleCount() { return _cAlleles; }
	int getHaplotypeCount();
	unsigned int getPos() { return _pos; }
	void calculateHighLowQuality();
	bool isHighConfidence();
    void determineReliable();
	bool isDefining() { return _bDefining; }
	bool isReliable() { return _bReliable; }
	void setFlankLength(unsigned int flankLength) { _flankLength = flankLength; }
	unsigned int getFlankLength() { return _flankLength; }
	const string toCSV();
	int getIUPAC();
    void setMajorAllele(char nucl);
    void setMinorAllele(char nucl);
    char getMajorAllele();
    char getMinorAllele();
    void setIsHighConfidence(bool bHighConfidence) { _bHighConfidence = bHighConfidence; }
    void setIsReliable(bool bReliable) {_bReliable = bReliable; }
    void setIsDefining(bool bDefining) { _bDefining = bDefining; }

private:
    bool insideHomopolymetricTract();
	int calculateScore(int);

	Contig* _pContig;
	unsigned int _pos;
	int* _HQNucCount;
	int* _LQNucCount;
	int* _nucCount;
	int _majorAllele;
	int _minorAllele;
    int _confidenceScore;
	int _cAlleles;
	int _cHaplotypes;
	unsigned int _flankLength;
	bool _bDefining;
	bool _bHighConfidence;
	bool _bReliable;
};


#endif

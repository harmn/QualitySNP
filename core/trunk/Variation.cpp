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

#include <iostream>
#include "Configuration.h"
#include "HaploType.h"
#include "SeqRead.h"
#include "Contig.h"
#include "Variation.h"

Variation::Variation(Contig* pContig, int pos) :
        _pContig(pContig), _pos(pos), _majorAllele(0), _minorAllele(-1), _confidenceScore(-1), _cAlleles(0),
         _flankLength(0), _bDefining(false), _bHighConfidence(false), _bReliable(false), _cHaplotypes(-1),
        _nucCount(NULL), _HQNucCount(NULL), _LQNucCount(NULL)
{
	int nDifferentNucs = Configuration::getConfig()->getNumberOfNucs();
	_HQNucCount = new int[nDifferentNucs];
	_LQNucCount = new int[nDifferentNucs];

	calculateHighLowQuality();
} 
Variation::~Variation(void)
{
	// clean up the reserve arrays
    if(_nucCount != NULL) {
        delete [] _nucCount;
    }

    if(_HQNucCount != NULL) {
        delete [] _HQNucCount;
    }

    if(_LQNucCount != NULL) {
        delete [] _LQNucCount;
    }
}

bool Variation::setNucCount(int* nucCount) {
	_nucCount = nucCount;

	int nDifferentNucs = Configuration::getConfig()->getNumberOfNucs();

	for(int iNuc = 1; iNuc < nDifferentNucs; iNuc++) {
		if(_nucCount[iNuc] > 0) { 
			_cAlleles++;
			if(_nucCount[iNuc] >= _nucCount[_majorAllele]) {
				_minorAllele = _majorAllele;
				_majorAllele = iNuc;
			} else if(_nucCount[iNuc] < _nucCount[_majorAllele]) {
				if(_minorAllele == -1 || _nucCount[iNuc] > _nucCount[_minorAllele]) {
					_minorAllele = iNuc;
				}
			}
		}
	}

	if (_minorAllele == -1) {
		// this should not happen!!
		cerr << "only one allele in this variation, not good..." << endl;
		return false;
	}

	return true;
}

// calculate the quality score for an allele
int Variation::calculateScore(int allele) {
	int score = 1;

	if(_HQNucCount[allele] > 1) {
		score = 5;
	} else if(_HQNucCount[allele] == 1 &&  _LQNucCount[allele] > 1) {
		score = 4;
	} else if(_LQNucCount[allele] > 3) {
		score = 3;
	} else if(_HQNucCount[allele] == 1 &&  _LQNucCount[allele] == 1) {
		score = 2;
	} else if(_LQNucCount[allele] == 3) {
		score = 2;
	}

	return score;
}

int Variation::getHaplotypeCount() {
	if(_cHaplotypes == -1) {
        determineReliable();
	}

	return _cHaplotypes;
}

// return the confidence score for this variation
int Variation::getConfidenceScore() {
    if(_confidenceScore == -1) {
		Configuration* pConfig = Configuration::getConfig();
		// when the position in the contig is already low quality, then so is this SNP
		if(_pContig->getQualityAt(_pos) < pConfig->getInt("minSNPQualityScore")) {
            _confidenceScore = 1;
		} else if (pConfig->getBool("indelLowQuality") && 
			(pConfig->int2nuc(_majorAllele) == '*' || pConfig->int2nuc(_minorAllele) == '*')) {
            _confidenceScore = 1;
        } else if (pConfig->getBool("lowComplexityLowQuality") && (insideHomopolymetricTract())) {
            _confidenceScore = 1;
		} else {
            _confidenceScore = min(calculateScore(_majorAllele),calculateScore(_minorAllele));
		}

        _bHighConfidence = _confidenceScore >= pConfig->getInt("minimalConfidenceScore");
	}

    return _confidenceScore;
}

bool Variation::insideHomopolymetricTract()
{
    Configuration* pConfig = Configuration::getConfig();
    int fragmentSize = pConfig->getInt("lowComplexityRegionSize");
    int repeatCount = pConfig->getInt("lowComplexityRepeatCount");

    char majorAllele = pConfig->int2nuc(_majorAllele);
    char minorAllele = pConfig->int2nuc(_minorAllele);

    string sequenceRegion(1,_pContig->getSequenceAt(_pos));

    int pos = _pos - 1;
    int nucCount = 1;
    // read nucleotides upstream, skipping gaps
    while(pos >= 0 && nucCount < fragmentSize) {
        char nuc = _pContig->getSequenceAt(pos);
        if(nuc != '*') {
            sequenceRegion = nuc + sequenceRegion;
            nucCount++;
        }
        pos--;
    }

    pos = _pos + 1;
    nucCount = 1;
    // read nucleotides downstream, skipping gaps
    while(pos >= 0 && nucCount != fragmentSize) {
        char nuc = _pContig->getSequenceAt(pos);
        if(nuc != '*') {
            sequenceRegion = sequenceRegion + nuc;
            nucCount++;
        }
        pos++;
    }

    for(int iPos = 0; iPos + fragmentSize <= sequenceRegion.length(); iPos++) {
        string fragment = sequenceRegion.substr(iPos, fragmentSize);
        map<char, int> nucMap;
        for(int iFrag = 0; iFrag < fragment.length(); iFrag++) {
            nucMap[fragment[iFrag]]++;
        }

        if((nucMap[majorAllele] >= repeatCount) || nucMap[minorAllele] >= repeatCount) {
            return true;
        }
    }

    return false;
}


// calculate the number of high and low quality nucleotides
void Variation::calculateHighLowQuality() {
    _confidenceScore = -1;

	Configuration* pConfig = Configuration::getConfig();
	int nDifferentNucs = pConfig->getNumberOfNucs();

	for(int iNuc= 0; iNuc < nDifferentNucs; iNuc++) {
		_HQNucCount[iNuc] = 0;
		_LQNucCount[iNuc] = 0;
	}

	const list<SeqRead*>& reads = _pContig->getReads();
	list<SeqRead*>::const_iterator itReads;

    for ( itReads= reads.begin(); itReads != reads.end(); itReads++) {
        char nuc = (*itReads)->getNucleotideAt(_pos);
        int nNuc = pConfig->nuc2int(nuc);
        if (nNuc != -1) {
            (*itReads)->isHighQuality(_pos) ? _HQNucCount[nNuc]++ : _LQNucCount[nNuc]++;
        }
    }
}

// a variation is high confidence if the quality score is at least equal to the minimal confidence score
bool Variation::isHighConfidence() {
    if (_confidenceScore == -1) {
        getConfidenceScore();
	}

	return _bHighConfidence;
}

// determine whether this variation is reliable and defining
// a variation cannot be reliable if it is not high quality
void Variation::determineReliable() {
	_cHaplotypes = 0;
	if(!isHighConfidence()) {
		return;
	}
	
	_bReliable = false;
	
	const list<HaploType*>& haploTypes = _pContig->getHaploTypes();
	list<HaploType*>::const_iterator itHaploTypes; 

	Configuration* pConfig = Configuration::getConfig();
	int nDifferentNucs = pConfig->getNumberOfNucs(); // for readability

	double wh = pConfig->getDouble("weightHighQualityRegion");
	double wl = pConfig->getDouble("weightLowQualityRegion");
	double alleleMajorityThreshold = pConfig->getDouble("alleleMajorityThreshold");

	HaploType ** nucHaploType = new HaploType*[nDifferentNucs];
	int* nucCount = new int[nDifferentNucs]; // keep a count of all nucleotides

	for(int iNuc= 0; iNuc < nDifferentNucs; iNuc++) {
		nucCount[iNuc] = 0;
		nucHaploType[iNuc] = NULL;
	}

	bool bMahap = false;
	bool bMihap = false;

	for ( itHaploTypes=haploTypes.begin(); itHaploTypes != haploTypes.end(); itHaploTypes++) {	
		list<SeqRead*> reads = (*itHaploTypes)->getReads();
		list<SeqRead*>::const_iterator itReads;

		int cInformativeReads = 0;
		int cMajorAlleleHQ = 0;
		int cMajorAlleleLQ = 0;
		int cMinorAlleleHQ = 0;
		int cMinorAlleleLQ = 0;

		for (itReads=reads.begin(); itReads != reads.end(); itReads++) {
			char nuc = (*itReads)->getNucleotideAt(_pos);
			int nNuc = pConfig->nuc2int(nuc);
			if (nNuc != -1) {
				cInformativeReads++;
				if(nNuc == _majorAllele) {
					((*itReads)->isHighQuality(_pos)) ? cMajorAlleleHQ++ : cMajorAlleleLQ++;
				} else if (nNuc == _minorAllele) {
					((*itReads)->isHighQuality(_pos)) ? cMinorAlleleHQ++ : cMinorAlleleLQ++;
				}

				// check if this is the first time we encounter this nucleotide
				// in this haplotype, if so: count
				if (nucHaploType[nNuc] != *itHaploTypes) {
					nucCount[nNuc]++;
				}
				nucHaploType[nNuc] = *itHaploTypes;
			}
		}

		if(cInformativeReads > 0) {
			_cHaplotypes++;
			double mahap_i = (wh * cMajorAlleleHQ + wl * cMajorAlleleLQ) / cInformativeReads;
			double mihap_i = (wh * cMinorAlleleHQ + wl * cMinorAlleleLQ) / cInformativeReads;

			if (mahap_i >= alleleMajorityThreshold) {
				bMahap = true;
			}
			if (mihap_i >= alleleMajorityThreshold) {
				bMihap = true;
			}
		}
	}
	
	// if both the major and the minor allele have a majority in a haplotype, this SNP is reliable
	if(bMihap && bMahap) {
		_bReliable = true;
	}

	for(int iNuc= 0; iNuc < nDifferentNucs; iNuc++) {
		if (nucCount[iNuc] == 1) {
			nucHaploType[iNuc]->addDefiningSNP(_pos);
			_bDefining = true;
		}
	}

	delete[] nucCount;
	delete[] nucHaploType;

	return;
}

// output this variation object as a CSV line
const string Variation::toCSV() {
	Configuration* pConfig = Configuration::getConfig();
    char sep = pConfig->getChar("fieldSeparator");
	stringstream csv;
	csv << _pContig->getName() << sep;
	csv << _pos << sep;
	csv << pConfig->int2nuc(_majorAllele) << sep;
	csv << pConfig->int2nuc(_minorAllele) << sep;
	csv << isHighConfidence() << sep;
	csv << isReliable() << sep;
	csv << isDefining() << sep;
	csv << _flankLength << NEWLINE;
	return csv.str();
}

int Variation::getIUPAC() {
	int nucs; 
	if (_majorAllele > _minorAllele) {
		nucs = 10*_minorAllele + _majorAllele;
	} else {
		nucs = 10*_majorAllele + _minorAllele;
	}

    return Configuration::getConfig()->getIUPACCode(nucs);
}

void Variation::setMajorAllele(char nucl)
{
    _majorAllele = Configuration::getConfig()->nuc2int(nucl);
}

void Variation::setMinorAllele(char nucl)
{
    _minorAllele = Configuration::getConfig()->nuc2int(nucl);
}

char Variation::getMajorAllele()
{
    return Configuration::getConfig()->int2nuc(_majorAllele);
}

char Variation::getMinorAllele()
{
    return Configuration::getConfig()->int2nuc(_minorAllele);
}

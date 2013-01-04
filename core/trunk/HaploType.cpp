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

#include <sstream>
#include "Configuration.h"
#include "SeqRead.h"
#include "Variation.h"
#include "Contig.h"
#include "HaploType.h"

HaploType::HaploType(Contig* pContig, int id) : 
    _pContig (pContig), _id(id), _bModified(true), _cReads(0), _lastVariablePosition(-1)
{
	_singleSNPThreshold = Configuration::getConfig()->getDouble("similarityPerPolymorphicSite");
	_allSNPThreshold = Configuration::getConfig()->getDouble("similarityAllPolymorphicSites");
}

HaploType::~HaploType(void)
{
}

void HaploType::matchRead(SeqRead* pRead, int& match, int& misMatch, map<int,int>& varNuc) {
    Configuration* pConfig = Configuration::getConfig();
    const vector<Variation*>& variations = _pContig->getVariations();
    vector<Variation*>::const_iterator itVar;
    for ( itVar=variations.begin(); itVar != variations.end(); itVar++) {
        if((*itVar)->isHighConfidence()) {
            int pos = (*itVar)->getPos();
            char cNuc = pRead->getNucleotideAt(pos);
            int iNuc = pConfig->nuc2int(cNuc);
            if(iNuc != -1) {
                varNuc[pos] = iNuc;
                int iMatch = matchAt(pos, cNuc);
                if(iMatch == 1) {
                    match++;
                } else if(iMatch == -1) {
                    misMatch++;
                }
            }
        }
    }
}

bool HaploType::tryAddRead(SeqRead* pRead) {
	int match = 0;
	int misMatch = 0;
	map<int,int> varNuc;

    matchRead(pRead, match, misMatch, varNuc);

	if(_reads.empty()) {
		// if this is the first read, always add it
		addRead(pRead, varNuc);
		return true;
	}

	// if no matches then this read certainly does not belong to this haplotype
	// unless the allSNPThreshold is zero, because then we do not care...
	if (match == 0) {
		return _allSNPThreshold == 0;
	}

	double similarity = 0.0;
	double threshold = _allSNPThreshold;

//	if(getReadCount() == 1 && match + misMatch == _reads[0]->getSNPCount()) {
	if(getReadCount() == 1 && misMatch == 0) {
        similarity = max(match/pRead->getSNPCount(true), match/(*_reads.begin())->getSNPCount(true));
		threshold = 0.5;
	} else {
		similarity = match*1.0 / (match + misMatch);
	}

	// similarity is lower than the threshold, so the read does not belong here
	if(similarity < threshold) {
		return false;
	}

	addRead(pRead, varNuc);
	return true;
}

void HaploType::addRead(SeqRead* pRead, const map<int, int>& varNuc) {
	Logger::getLogger()->log(QSNP_DEBUG, "HaploType::addRead");
	_bModified = true;
	_reads.push_back(pRead);
	_cReads++;
	pRead->setHaploType(this);

    map<int,int>::const_iterator itVar;

	Configuration* pConfig = Configuration::getConfig();

	for(itVar = varNuc.begin(); itVar != varNuc.end(); itVar++) {
		if (_varNuc.count((*itVar).first) == 0) {
				vector<int> rgNuc (pConfig->getNumberOfNucs(), 0);
				_varNuc[(*itVar).first] = rgNuc;
		}
		_varNuc[(*itVar).first][(*itVar).second]++; 
	}
}

void HaploType::addRead(SeqRead* pRead) {
    int match = 0;
    int misMatch = 0;
    map<int,int> varNuc;
    matchRead(pRead,match,misMatch,varNuc);

    addRead(pRead,varNuc);
}

char HaploType::getNucleotideAt(int pos)
{
    if(_varNuc.count(pos) == 0) {
        // return a space if there is not information at this position
        return ' ';
    }
    Configuration* pConfig = Configuration::getConfig();

    int iMax = 0;
    for (int i = 1; i < pConfig->getNumberOfNucs(); i++) {
        iMax = (_varNuc[pos][iMax] < _varNuc[pos][i]) ? i : iMax;
    }

    return pConfig->int2nuc(iMax);
}

int HaploType::getLastVariablePosition()
{
    if(_lastVariablePosition == -1) {
        map<int, vector<int> >::iterator itVar = _varNuc.begin();
        for(itVar = _varNuc.begin(); itVar != _varNuc.end(); itVar++){
            if(itVar->first > _lastVariablePosition) {
                _lastVariablePosition = itVar->first;
            }
        }
    }

    return _lastVariablePosition;
}

int HaploType::matchAt(int pos, char queryNucleotide) {
	if(_varNuc.count(pos) == 0) {
		// no for this position (yet)
		return 0;
	}

	Configuration* pConfig = Configuration::getConfig();
	int iNuc = pConfig->nuc2int(queryNucleotide);
	if(iNuc == -1) {
		// not a normal nucleotide, ignore
		return 0;
	}

	int match = _varNuc[pos][iNuc];
	if (match == 0) {
		// if there are no matches, there can only be mismatches
		return -1;
	}

	int misMatch = 0;

	for (int i = 0; i < pConfig->getNumberOfNucs(); i++) {
		if (i != iNuc) {
			misMatch += _varNuc[pos][i];
		}
	}

	double similarity = match*1.0 / (match + misMatch);

	return (similarity < _singleSNPThreshold) ? -1 : 1;
}

const string HaploType::toString() {
	stringstream result;
	string separator;
	list<SeqRead*>::const_iterator itReads;
	for ( itReads=_reads.begin(); itReads != _reads.end(); itReads++) {
		result << separator;
		const string& readName = (*itReads)->getName();
		result << _pContig->getReadIndex(readName);
		separator = " ";
	}

	return result.str();
}

const string HaploType::toCSV() {
	char sep = Configuration::getConfig()->getChar("fieldSeparator");
	stringstream csv;
	csv << _pContig->getName() << sep << _id << NEWLINE;
	return csv.str();
}


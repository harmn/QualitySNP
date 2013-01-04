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

#include <list>
#include "Logger.h"
#include "SeqRead.h"
#include "Contig.h"
#include "SAMRead.h"
#include "SAMContig.h"

using namespace std;

SAMContig::SAMContig()
{
}

SAMContig::~SAMContig() {
		for(list<SAMRead*>::iterator it = _reads.begin(); it != _reads.end(); ++it) {
			delete *it;
		}
}

// construct the reference sequence from the information in the individual reads
bool SAMContig::constructReferenceSequence() {
    Logger::getLogger()->log(QSNP_INFO, "construction reference sequence for: " + _name);
	list<SAMRead*>::iterator it;
	unsigned int length = 0;
	for ( it=_reads.begin(); it != _reads.end(); it++) {
		unsigned int readLastPos = (*it)->getLastPosition();
		length = max(length, readLastPos);
	}

	_sequence.reserve(length);

	Configuration* pConfig = Configuration::getConfig();
	int nDifferentNucs = pConfig->getNumberOfNucs(); // for readability
	int* nucCount = new int[nDifferentNucs]; // keep a count of all nucleotides

    vector< list<SAMRead*> > startPosReads;
    startPosReads.resize(length + 1);
    for ( it=_reads.begin(); it != _reads.end(); it++) {
        startPosReads[(*it)->getStartPosition() + 1].push_back(*it);
    }

    list<SAMRead*> reads;
	for(unsigned int pos = 1; pos <= length; pos++) {
        while(!startPosReads[pos].empty()) {
            reads.push_back(startPosReads[pos].front());
            startPosReads[pos].pop_front();
        }

		for(int iNuc= 0; iNuc < nDifferentNucs; iNuc++) {
			nucCount[iNuc] = 0;
		}

		bool bInformative = false;
        for ( it=reads.begin(); it != reads.end(); ) {
            if(pos > (*it)->getLastPosition()) {
                it = reads.erase(it);
            } else {
                char readNucl =(*it)->getNucleotideAt(pos);
                int iNuc = pConfig->nuc2int(readNucl);
                if(iNuc != -1) {
                    nucCount[iNuc]++;
                    bInformative = true;
                }
                it++;
            }
		}

        char nucl = 'N';

		if(bInformative) {
			int iMaxNuc = 0;
			for(int iNuc= 1; iNuc < nDifferentNucs; iNuc++) {
				if(nucCount[iNuc] > nucCount[iMaxNuc]) {
					iMaxNuc = iNuc;
				}
			}
			nucl = pConfig->int2nuc(iMaxNuc);
		}

		_sequence += nucl;
	}

	delete[] nucCount;

	return true;
}

void SAMContig::stretchReads() {
    Logger::getLogger()->log(QSNP_INFO, "SAMContig: stretching reads to include insertions");
    list<SAMRead*>::iterator it;
    unsigned int pos = 1;

    list<SAMRead*> reads = _reads;
    list<SAMRead*> readsWithInsertions;
    for ( it=_reads.begin(); it != _reads.end(); it++) {
        if((*it)->hasInsertions()) {
            readsWithInsertions.push_back(*it);
        }
    }

    while (!readsWithInsertions.empty()) {
        bool bInsert = false;
        it = readsWithInsertions.begin();
        while(it != readsWithInsertions.end() && !bInsert) {
            char operation =(*it)->getOperationAt(pos);
            if(operation == 'I') {
                bInsert = true;
            }
            if(operation == '>') {
                it = readsWithInsertions.erase(it);
            } else {
                it++;
            }
        }
        if (bInsert) {
            for ( it=reads.begin(); it != reads.end();) {
                bool bInside = (*it)->insertGapAt(pos);
                if(bInside) {
                    it++;
                } else {
                    it = reads.erase(it);
                }
            }
        }
        pos++;
    }
}

void SAMContig::mergeReadPairs() {
    Logger::getLogger()->log(QSNP_INFO, "SAMContig: merging read pairs");
	list<SAMRead*>::iterator it;
	
	map<string, list<SAMRead*>::iterator> readMap;
	for ( it=_reads.begin(); it != _reads.end();) {
		string readName = (*it)->getName();
		
		if(readMap.count(readName) == 1) {
            Logger::getLogger()->log(QSNP_DEBUG, "Merging read pair with name: " + readName);
			(*readMap[readName])->merge(*it);
			delete *it;
			it = _reads.erase(it);
		} else {
			readMap[readName] = it;
			++it;
		}
	}
}

Contig* SAMContig::toContig() {
	if(_name.empty()) {
		return NULL;
	}
    Logger::getLogger()->log(QSNP_INFO, "SAMContig: converting to Contig");

	Contig* pContig = new Contig(_name);
	pContig->setSequence(_sequence);
	for(list<SAMRead*>::iterator itRead = _reads.begin(); itRead != _reads.end(); itRead++) {
		SeqRead* pSeqRead = (*itRead)->toSeqRead();
		pSeqRead->setContig(pContig);
		pContig->addRead(pSeqRead);
	}
	
	return pContig;
}


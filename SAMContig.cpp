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

	for(unsigned int pos = 1; pos <= length; pos++) {
		for(int iNuc= 0; iNuc < nDifferentNucs; iNuc++) {
			nucCount[iNuc] = 0;
		}

		bool bInformative = false;
		for ( it=_reads.begin(); it != _reads.end(); it++) {
			char readNucl =(*it)->getNucleotideAt(pos);
			int iNuc = pConfig->nuc2int(readNucl);
			if(iNuc != -1) {
				nucCount[iNuc]++;
				bInformative = true;
			}
		}

		char nucl = ' ';

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
	list<SAMRead*>::iterator it;
	unsigned int pos = 1;
	bool bDone = false;
	while (!bDone) {
		bool bInsert = false;
		bDone = true;
		for ( it=_reads.begin(); it != _reads.end(); it++) {
			char operation =(*it)->getOperationAt(pos);
			if(operation == 'I') {
				bInsert = true;
			}
			if(operation != '>') {
				bDone = false;
			}
		}
		if (bInsert) {
			for ( it=_reads.begin(); it != _reads.end(); it++) {
				(*it)->insertGapAt(pos);
			}
		}
		pos++;
	}
}

void SAMContig::mergeReadPairs() {
	list<SAMRead*>::iterator it;
	
	map<string, list<SAMRead*>::iterator> readMap;
	for ( it=_reads.begin(); it != _reads.end();) {
		string readName = (*it)->getName();
		
		if(readMap.count(readName) == 1) {
            Logger::getLogger()->log(QSNP_INFO, "Merging read pair with name: " + readName);
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

	Contig* pContig = new Contig(_name);
	pContig->setSequence(_sequence);
	for(list<SAMRead*>::iterator itRead = _reads.begin(); itRead != _reads.end(); itRead++) {
		SeqRead* pSeqRead = (*itRead)->toSeqRead();
		pSeqRead->setContig(pContig);
		pContig->addRead(pSeqRead);
	}
	
	return pContig;
}


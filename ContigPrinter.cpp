#include <sstream>
#include <iomanip>
#include "Variation.h"
#include "HaploType.h"
#include "SeqRead.h"
#include "Configuration.h"
#include "Contig.h"
#include "ContigPrinter.h"

ContigPrinter::ContigPrinter(Contig* pContig): _pContig(pContig)
{
}

ContigPrinter::~ContigPrinter(void)
{
}

string ContigPrinter::toString() {
	Configuration* pConfig = Configuration::getConfig();

	stringstream result;

	if(pConfig->getBool("printSummaryLine")) {
		printSummaryLine(result);
	}

	if (pConfig->getBool("printSNPs")) {
		printSNPs(result, false);
	}

	if (pConfig->getBool("printHaploTypes")) {
		printHaploTypes(result);
	}

	if (pConfig->getBool("printAlignment")) {
		printAlignment(result);
	}

	return result.str();
}

void ContigPrinter::printSummaryLine(stringstream& result) {
	result << _pContig->getName();
	result << "; ESTs: "					<< _pContig->getReadCount();
	result << "; potential SNPs: "			<< _pContig->getPotentialSNPCount();
	result << "; high quality SNPs: "		<< _pContig->getHighConfidenceSNPCount();
	result << "; reliable SNPs: "			<< _pContig->getReliableSNPCount();
	result << "; HaploTypes: "				<< _pContig->getHaploTypeCount() + _pContig->getSingleESTHaploTypeCount();
	result << "; max HaploTypes/SNP: "		<< _pContig->getMaxHaploTypePerSNPCount();
	result << "; Single read HaploTypes: "	<< _pContig->getSingleESTHaploTypeCount();
	result << "; D Value: "					<< _pContig->getDvalue() << endl;
}


void ContigPrinter::printMarkerSNPs(stringstream& result, unsigned int flankSize) {
	const vector<Variation*> &variations = _pContig->getVariations();
	vector<Variation*>::const_iterator itVariations;
	
	for (itVariations = variations.begin(); itVariations != variations.end(); itVariations++) {
		if((*itVariations)->getFlankLength() >= flankSize) {
			result << setw((*itVariations)->getPos() - flankSize + 1) << "|" << setfill('-') << setw(2 * flankSize - 1) << "-" << setw(1) << "|" << endl;
		}
	}
}

void ContigPrinter::printAlignment(stringstream& result) {
	const list<SeqRead*> &reads			= _pContig->getReads();
	const vector<Variation*> &variations	= _pContig->getVariations();

	list<SeqRead*>::const_iterator it;
	for ( it=reads.begin(); it != reads.end(); it++) {
		result << (*it)->getName() << endl;
		result << (*it)->toString() << endl;
	}

	string strSNPs((int) _pContig->getSequenceLength(),' ');
	vector<Variation*>::const_iterator itVar;
	for ( itVar=variations.begin(); itVar !=variations.end(); itVar++) {
		strSNPs[(*itVar)->getPos()] = '#';
	}

	result << strSNPs << endl;

	string strHQSNPs((int) _pContig->getSequenceLength(),' ');
	for ( itVar=variations.begin(); itVar < variations.end(); itVar++) {
		if((*itVar)->isHighConfidence()) {
			strHQSNPs[(*itVar)->getPos()] = '$';
		}
	}

	result << strHQSNPs << endl;

	string strReliSNPs((int) _pContig->getSequenceLength(),' ');
	for ( itVar=variations.begin(); itVar < variations.end(); itVar++) {
		if((*itVar)->isReliable()) {
			strReliSNPs[(*itVar)->getPos()] = '+';
		}
	}

	result << strReliSNPs << endl;
}

void ContigPrinter::printHaploTypes(stringstream& result) {
	const list<HaploType*> &haploTypes = _pContig->getHaploTypes();
	list<HaploType*>::const_iterator itHaploTypes;

	int counter = 0;
	for ( itHaploTypes=haploTypes.begin(); itHaploTypes != haploTypes.end(); itHaploTypes++) {
		counter++;
		result << counter << " haplotype including ";
		result << (*itHaploTypes)->getReadCount() << " ESTs:";
		result << (*itHaploTypes)->toString() << endl;
	}
	result << endl;
}

void ContigPrinter::printSNPs(stringstream& result, bool bHQOnly) {
	const list<HaploType*> &haploTypes = _pContig->getHaploTypes();
	const vector<Variation*> &variations = _pContig->getVariations();

	list<SeqRead*>::const_iterator itReads;
	vector<Variation*>::const_iterator itVar;
	list<HaploType*>::const_iterator itHaploTypes;
	for ( itHaploTypes = haploTypes.begin(); itHaploTypes != haploTypes.end(); itHaploTypes++) {
		list<SeqRead*>  reads = (*itHaploTypes)->getReads();
		for ( itReads=reads.begin(); itReads != reads.end(); itReads++) {
			const string& readName = (*itReads)->getName();
			result << readName << " (" << _pContig->getReadIndex(readName) << ")" << ":\t";
			for ( itVar=variations.begin(); itVar != variations.end(); itVar++) {
				if (!bHQOnly || (*itVar)->isHighConfidence()) {
					int pos = (*itVar)->getPos();
					result << (*itReads)->getNucleotideAt(pos);
					result << " ";
				}
			}
			result << endl;
		}
		result << endl;
	}
	result << endl;
}


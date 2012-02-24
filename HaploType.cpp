#include <sstream>
#include "Configuration.h"
#include "SeqRead.h"
#include "Variation.h"
#include "SeqRead.h"
#include "Contig.h"
#include "HaploType.h"

HaploType::HaploType(Contig* pContig, int id) : 
	_pContig (pContig), _id(id), _bModified(true), _bComplete(false), _cReads(0)
{
	_singleSNPThreshold = Configuration::getConfig()->getDouble("similarityPerPolymorphicSite");
	_allSNPThreshold = Configuration::getConfig()->getDouble("similarityAllPolymorphicSites");
}

HaploType::~HaploType(void)
{
}

bool HaploType::tryAddRead(SeqRead* pRead, bool bHQOnly) {
	Configuration* pConfig = Configuration::getConfig();

	int match = 0;
	int misMatch = 0;
	map<int,int> varNuc;

	vector<Variation*> variations = _pContig->getVariations();
	vector<Variation*>::const_iterator itVar;
	for ( itVar=variations.begin(); itVar != variations.end(); itVar++) {
		if(!bHQOnly || (*itVar)->isHighConfidence()) {
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
		similarity = max(match/pRead->getSNPCount(bHQOnly), match/(*_reads.begin())->getSNPCount(bHQOnly));
		threshold = 0.5;
	} else {
		similarity = match*1.0 / (match + misMatch);
	}

	// similarity is lower then the threshold, so the read does not belong here
	if(similarity < threshold) {
		return false;
	}

	addRead(pRead, varNuc);
	return true;
}

void HaploType::addRead(SeqRead* pRead, const map<int, int>& varNuc) {
	Logger::getLogger()->log(QSNP_INFO, "HaploType::addRead");
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


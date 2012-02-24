#include <iostream>
#include "Configuration.h"
#include "HaploType.h"
#include "SeqRead.h"
#include "Contig.h"
#include "Variation.h"

Variation::Variation(Contig* pContig, int pos) :
		_pContig(pContig), _pos(pos), _majorAllele(0), _minorAllele(-1), _qualityScore(-1), _cAlleles(0),
         _flankLength(0), _bDefining(false), _bHighConfidence(false), _mahap(0.0), _mihap(0.0), _cHaplotypes(-1)
{
	int nDifferentNucs = Configuration::getConfig()->getNumberOfNucs();
	_HQNucCount = new int[nDifferentNucs];
	_LQNucCount = new int[nDifferentNucs];

	calculateHighLowQuality();
} 
Variation::~Variation(void)
{
	// clean up the reserve arrays
	delete [] _nucCount;
	delete [] _HQNucCount;
	delete [] _LQNucCount;
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
		determineDefining();
	}

	return _cHaplotypes;
}

// return the quality score for this variation
int Variation::getQualityScore() {
	if(_qualityScore == -1) {
		Configuration* pConfig = Configuration::getConfig();
		// when the position in the contig is already low quality, then so is this SNP
		if(_pContig->getQualityAt(_pos) < pConfig->getInt("minSNPQualityScore")) {
			_qualityScore = 1;
		} else if (pConfig->getBool("indelLowQuality") && 
			(pConfig->int2nuc(_majorAllele) == '*' || pConfig->int2nuc(_minorAllele) == '*')) {
				_qualityScore = 1;
		} else {
			_qualityScore = min(calculateScore(_majorAllele),calculateScore(_minorAllele));
		}

		_bHighConfidence = _qualityScore >= pConfig->getInt("minimalConfidenceScore");
	}

	return _qualityScore;
}

// calculate the number of high and low quality nucleotides
void Variation::calculateHighLowQuality() {
	_qualityScore = -1;

	Configuration* pConfig = Configuration::getConfig();
	int nDifferentNucs = pConfig->getNumberOfNucs();

	for(int iNuc= 0; iNuc < nDifferentNucs; iNuc++) {
		_HQNucCount[iNuc] = 0;
		_LQNucCount[iNuc] = 0;
	}

	const list<SeqRead*>& reads = _pContig->getReads();
	list<SeqRead*>::const_iterator itReads;

	for ( itReads= reads.begin(); itReads != reads.end(); itReads++) {
		//if ((*itReads)->getHaploType() == NULL || (*itReads)->getHaploType()->getReadCount() > 1) {
			char nuc = (*itReads)->getNucleotideAt(_pos);
			int nNuc = pConfig->nuc2int(nuc);
			if (nNuc != -1) {
				(*itReads)->isHighQuality(_pos) ? _HQNucCount[nNuc]++ : _LQNucCount[nNuc]++;
			}
		//}
	} 
}

// a variation is high confidence if the quality score is at least equal to the minimal confidence score
bool Variation::isHighConfidence() {
	if (_qualityScore == -1) {
		getQualityScore();
	}

	return _bHighConfidence;
}

// determine whether this variation is a defining one
// a variation cannot be reliable if it is not high quality
void Variation::determineDefining() {
	_cHaplotypes = 0;
	if(!isHighConfidence()) {
		return;
	}

	_mahap = 0.0;
	_mihap = 0.0;

	const list<HaploType*>& haploTypes = _pContig->getHaploTypes();
	list<HaploType*>::const_iterator itHaploTypes; 

	Configuration* pConfig = Configuration::getConfig();
	int nDifferentNucs = pConfig->getNumberOfNucs(); // for readability

	double wh = pConfig->getDouble("weightHighQualityRegion");
	double wl = pConfig->getDouble("weightLowQualityRegion");

	HaploType ** nucHaploType = new HaploType*[nDifferentNucs];
	int* nucCount = new int[nDifferentNucs]; // keep a count of all nucleotides

	for(int iNuc= 0; iNuc < nDifferentNucs; iNuc++) {
		nucCount[iNuc] = 0;
		nucHaploType[iNuc] = NULL;
	}

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

			if (mahap_i >= 0.75) {
				mahap_i += 1;
			}
			if (mihap_i >= 0.75) {
				mihap_i += 1;
			}

			_mahap += mahap_i;
			_mihap += mihap_i;
		}
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

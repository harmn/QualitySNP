#include <sstream>
#include <algorithm>
#include <math.h>
#include "SeqRead.h"
#include "Variation.h"
#include "Configuration.h"
#include "HaploType.h"
#include "Contig.h"

using namespace std;

Contig::Contig(const string& name) :
	_name(name), _cPotentialSNP(-1), _cHighConfidenceSNP(-1), _cReliableSNP(-1), _Dvalue(-1), _cReads(0)
{
	Logger::getLogger()->log(QSNP_INFO, "New contig: " + _name);
    _defaultQualityScore = Configuration::getConfig()->getInt("minSNPQualityScore");
}

Contig::~Contig(void) {
	while (!_reads.empty()) {
		delete static_cast<SeqRead*> (_reads.back());
		_reads.pop_back();
	}

	while(!_variations.empty()) {
		delete static_cast<Variation*> (_variations.back());
		_variations.pop_back();
	}

	while(!_haploTypes.empty()) {
		delete static_cast<HaploType*> (_haploTypes.back());
		_haploTypes.pop_back();
	}

	while(!_singleReadHaploTypes.empty()) {
		delete static_cast<HaploType*> (_singleReadHaploTypes.back());
		_singleReadHaploTypes.pop_back();
	}
}

// set the quality score at a sequence position, return false if position out of range
bool Contig::setQualityAt(unsigned int pos, int quality) {
	if(pos >= _quality.size()) {
		return false;
	}

	_quality[pos] = quality;

	return true;
}

// the main work is started here
void Contig::calculateProperties() {
	Logger::getLogger()->log(QSNP_INFO, "Contig::calculateProperties");
	determineVariations();
	calculateSNPCountPerRead();
	determineHaploTypes(true);
	removeSingleReadHaplotypes();
	calculateDefiningSNPs();
	findMarkerSNPs();
}

// set the sequence for this contig and make it upper case
void Contig::setSequence(const string& sequence) {
	Logger::getLogger()->log(QSNP_INFO, "Contig::setSequence");
	_sequence = sequence;

	transform (_sequence.begin (), _sequence.end (), _sequence.begin (), (int(*)(int)) toupper);
}

// return the nucleotide at a position
char Contig::getSequenceAt(unsigned int pos) {
	if (pos >  _sequence.length() -1) {
		return ' ';
	} else {
		return _sequence[pos];
	}
}

// return the quality score at a certain position in the contig sequence
int Contig::getQualityAt(unsigned int pos) { 
	if (pos + 1 > _quality.size()) {
		return _defaultQualityScore;
	}
	return _quality[pos]; 
}

// set the quality for the nucleotides in this contig
void Contig::setQuality(const string& quality) {
	Logger::getLogger()->log(QSNP_INFO, "Contig::setQuality");
	istringstream isQuality;
	isQuality.str(quality);
	_quality.reserve(_sequence.length());
	int count = 0;
	while(!isQuality.eof()) {
		while(getSequenceAt(count) == '*') {
			_quality.push_back(-1);
			count++;
		}

		int iQuality = -1;
		isQuality >> iQuality;
		if (iQuality > -1) {
			_quality.push_back(iQuality);
			count++;
		}
	}
}

// add a read to this contig, check that it does not conflict with an already present one
bool Contig::addRead(SeqRead* pRead) {
	Logger::getLogger()->log(QSNP_INFO, string("Adding read: ") + pRead->getName());
	if (_readMap.count(pRead->getName().c_str()) > 0) {
		Logger::getLogger()->log(QSNP_ERROR, "Trying to add the same read twice to contig " + _name + ", read name " + pRead->getName());
		return false;
	}	

	_reads.push_back(pRead);
	_cReads++;
	_readMap[pRead->getName().c_str()] = pRead;

	return true;
}

// search the read by its name and return the pointer if found, or NULL if not found
SeqRead* Contig::findRead(const string& name) {
	if(_readMap.count(name) == 0) {
		return NULL;
	}

	return _readMap[name];
}

// look for a read in the reads list
int Contig::getReadIndex(const string& name) {
	int readIndex = -1;

	list<SeqRead*>::iterator it;
	for ( it=_reads.begin(); it != _reads.end(); it++) {
		readIndex++;
		if ((*it)->getName() == name) {
			return readIndex;
		}
	}
	
	return readIndex;
}

// return the number of potential SNPs, i.e. before quality filtering
int Contig::getPotentialSNPCount() {
	if(_cPotentialSNP == -1) {
		// lazy evaluation
		determineVariations();
	} 

	return _cPotentialSNP;
}

// return the number of high confidence SNPs
int Contig::getHighConfidenceSNPCount() {
	if (_cHighConfidenceSNP != -1) {
		return _cHighConfidenceSNP;
	}

	// lazy evaluation
	Logger::getLogger()->log(QSNP_INFO, "Calculating High confidence SNP count");
	_cHighConfidenceSNP = 0;

	vector<Variation*>::iterator itVariations;

	for ( itVariations=_variations.begin(); itVariations != _variations.end(); itVariations++) {
		if((*itVariations)->isHighConfidence()) {
			_cHighConfidenceSNP++;
		}
	}
	return _cHighConfidenceSNP;	
}

int Contig::getReliableSNPCount() {
	if (_cReliableSNP != -1) {
		return _cReliableSNP;
	}

	Logger::getLogger()->log(QSNP_INFO, "Calculating reliable SNP count");
	_cReliableSNP = 0;

	vector<Variation*>::iterator itVariations;

	for ( itVariations=_variations.begin(); itVariations != _variations.end(); itVariations++) {
		if((*itVariations)->isReliable()) {
			_cReliableSNP++;
		}
	}
	return _cReliableSNP;	
}

void Contig::determineVariations() {
	Logger::getLogger()->log(QSNP_INFO, "Contig::determineVariations");
	Configuration* pConfig = Configuration::getConfig();
	int nDifferentNucs = pConfig->getNumberOfNucs();
	int nMinAlleles = pConfig->getInt("minimalNumberOfReadsPerAllele");
	double dMinAllelesp = pConfig->getDouble("minimalNumberOfReadsPerAllelep");
	
	for (int iSeq = 0; iSeq < (int)_sequence.length(); iSeq++) {
		list<SeqRead*>::iterator itReads;

		int* nucCount = new int[nDifferentNucs];
		for(int iNuc= 0; iNuc < nDifferentNucs; iNuc++) {
			nucCount[iNuc] = 0;
		}

		int cInformativeReads = 0;
		for (itReads=_reads.begin(); itReads != _reads.end(); itReads++) {
			char nuc = (*itReads)->getNucleotideAt(iSeq);
			int nNuc = pConfig->nuc2int(nuc);
			if (nNuc != -1) {
				cInformativeReads++;
				nucCount[nNuc]++;
			} else {
				// TODO: handle unknown nucleotides...
				//Logger::getLogger()->log(QSNP_ERROR, "Found unknown nucleotide: " + nuc);
			}
		}

		int cNuc = 0;
		int tmpMinAllelesp = static_cast<int>(cInformativeReads * dMinAllelesp) + 1;
		int tmpMinAlleles = (nMinAlleles > tmpMinAllelesp) ? nMinAlleles : tmpMinAllelesp;
		for ( int iNuc = 0; iNuc < nDifferentNucs; iNuc++) {
			if(nucCount[iNuc] < tmpMinAlleles) {
				// no valid SNP
			} else {
				cNuc++;
			}
		}

		if (cNuc > 1) {
			Variation* var = new Variation(this, iSeq);
			var->setNucCount(nucCount);
			_variations.push_back(var);
		} else {
			delete [] nucCount;
		}
	}

	_cPotentialSNP = _variations.size();
}

void Contig::calculateDefiningSNPs() {
	// iterate through all the variations of this contig 
	// and check whether it is defining for this haplotype
	Logger::getLogger()->log(QSNP_INFO, "Calculating defining SNPs");
	vector<Variation*>::iterator itVariations;
	for (itVariations=_variations.begin(); itVariations != _variations.end(); itVariations++) {
		(*itVariations)->determineDefining();
	}
}

int Contig::getMaxHaploTypePerSNPCount() {
	// iterate through all the variations of this contig 
	// and count the number of Haplotypes per variation. Return the largest count
	Logger::getLogger()->log(QSNP_INFO, "Calculating Max number of haplotypes per SNP");
	vector<Variation*>::iterator itVariations;
	int maxHaplotypes = 0;
	for (itVariations=_variations.begin(); itVariations != _variations.end(); itVariations++) {
		int cHaplotypes = (*itVariations)->getHaplotypeCount();
		maxHaplotypes = max(maxHaplotypes, cHaplotypes);
	}

	return maxHaplotypes;
}


void Contig::calculateSNPCountPerRead() {
	Logger::getLogger()->log(QSNP_INFO, "Calculating SNP count per read");
	list<SeqRead*>::iterator itReads;
	for ( itReads=_reads.begin(); itReads != _reads.end(); itReads++) {
		(*itReads)->calculateSNPCount();
	}

}

int Contig::getHaploTypeCount() {
	if(_haploTypes.empty()) {
		// lazy evaluation
		determineHaploTypes(true);
	}
	
	return _haploTypes.size();
}

// find the haplotypes in this contig
// if the bHQOnly flag is set, only the high quality variations are
// taken into account
void Contig::determineHaploTypes(bool bHQOnly) {
	Logger::getLogger()->log(QSNP_INFO, "Determining HaploTypes");
	if (_reads.empty()) {
		Logger::getLogger()->log(QSNP_WARNING, "No reads");
		return;
	}
	list<SeqRead*>::iterator itReads;

	int cNotYetAdded = 1;
	int cHaplo = 0;

	while (cNotYetAdded > 0) {
		cNotYetAdded = 0;
		cHaplo++;
		HaploType* pHaploType = new HaploType(this, cHaplo);
		_haploTypes.push_back(pHaploType);

		bool bReadAdded = true;

		// try to add reads to the existing haplotypes
		while(bReadAdded) {
			bReadAdded = false;
			for ( itReads = _reads.begin(); itReads != _reads.end(); itReads++) {
				if ((*itReads)->getSNPCount(bHQOnly) > 0) {
					// room for optimization
					if ((*itReads)->getHaploType() == NULL) {
						if (!addReadToHaploType(*itReads, bHQOnly)) {
							cNotYetAdded++;
						} else {
							bReadAdded = true;
						}
					}
				}
			}
			checkHaplotypesModified();
		}

		// if no read was added to this haplotype, delete it.
		if(pHaploType->getReadCount() == 0) {
			cHaplo--;
			_haploTypes.pop_back();
			delete pHaploType;
		}
	} 
}

// tries to add the read to one of this contigs haplotypes
// if the bHQOnly is true, only high quality variations are taken into account
// returns false if the read could not be added to one of the existing haplotypes
bool Contig::addReadToHaploType(SeqRead* pRead, bool bHQOnly) {
	Logger::getLogger()->log(QSNP_INFO, "Contig::addReadToHaploType");

	// reverse_iterator should be faster, because we start with the newest haplotype
	// but it also leads to different results, probably because of the order dependence
//	vector<HaploType*>::reverse_iterator itHaploTypes;
	list<HaploType*>::iterator itHaploTypes;
//	for ( itHaploTypes=_haploTypes.rbegin(); itHaploTypes != _haploTypes.rend(); itHaploTypes++) {
	for ( itHaploTypes=_haploTypes.begin(); itHaploTypes != _haploTypes.end(); itHaploTypes++) {
		if(!(*itHaploTypes)->complete() && (*itHaploTypes)->tryAddRead(pRead, bHQOnly)) {
			return true;
		}
	}

	return false;
}

bool readSortFunction(SeqRead* s1, SeqRead* s2) {
	return (s1->getName().compare(s2->getName()) < 0);	
}

// sort the reads alphabetically according to their names
void Contig::sortReads() {
	_reads.sort(readSortFunction);
	//sort(_reads.begin(), _reads.end(), readSortFunction);
}

// calculate the D-Value for this contig
double Contig::getDvalue() {
	if(_Dvalue != -1) {
		// lazy evaluation
		return _Dvalue;
	}

	Logger::getLogger()->log(QSNP_INFO, "Calculating the D value");
	_Dvalue = 0.0;

	if (_haploTypes.empty()) {
		return _Dvalue;
	}

	// remove all haplotypes with only one read
	// calculate the number of potential SNPs defining each haplotype
	// normalize this number: nrm_snp(i) = snp(i) / sum(snp(i)/ ahap
	// D is the standard deviation of normalized number of potential SNPs

	int cDefiningTotal = 0;
	list<HaploType*>::iterator itHaploTypes;
	// iterate through all haplotypes for this contig
	for ( itHaploTypes=_haploTypes.begin(); itHaploTypes != _haploTypes.end(); itHaploTypes++) {
		int cDefining = (*itHaploTypes)->getDefiningSNPCount();
		cDefiningTotal += cDefining;
	}

	double avgDefiningSNPs = (cDefiningTotal * 1.0) / _haploTypes.size();
	
	// return 0.0 when the average defining SNPs is too low
	if(avgDefiningSNPs < 0.1 ) {
		return _Dvalue;
	}

	double sum = 0.0;
	for ( itHaploTypes=_haploTypes.begin(); itHaploTypes != _haploTypes.end(); itHaploTypes++) {
		double nrmsnp = ((*itHaploTypes)->getDefiningSNPCount() * 1.0)/avgDefiningSNPs;
		sum += (nrmsnp - 1.0)*(nrmsnp - 1.0);
	}

	// calculate the D-value
	_Dvalue = sqrt(sum/_haploTypes.size()); 

	return _Dvalue;
}

// single read haplotypes are removed from the haplotype list
void Contig::removeSingleReadHaplotypes() {
	Logger::getLogger()->log(QSNP_INFO, "Contig::removeSingleReadHaplotypes");
	list<HaploType*>::iterator itHaploTypes;
	itHaploTypes=_haploTypes.begin();
	while(itHaploTypes != _haploTypes.end()) {
		if ((*itHaploTypes)->getReadCount() == 1) {
			_singleReadHaploTypes.push_back(*itHaploTypes);
			itHaploTypes = _haploTypes.erase(itHaploTypes);
		} else {
			itHaploTypes++;
		}
	}

	// recalculate the quality of each variation
	vector<Variation*>::iterator itVariations;
	for (itVariations=_variations.begin(); itVariations != _variations.end(); itVariations++) {
		(*itVariations)->calculateHighLowQuality();
	}

	_cHighConfidenceSNP = -1;
}

// find the high quality flanks for variations of high quality
// return the count
int Contig::findMarkerSNPs() {
	Logger::getLogger()->log(QSNP_INFO, "findMarkerSNPs");

	int cMarkerSNP = 0;
	unsigned int nAllowedSNPs = Configuration::getConfig()->getInt("maxNumberOfSNPsInFlanks");
	bool bOnlyReliableMarkers = Configuration::getConfig()->getBool("onlyReliableMarkers");

	stringstream result;
	for (unsigned int iVar = 0; iVar < _variations.size(); iVar++) {
		if ((bOnlyReliableMarkers == false && _variations[iVar]->isHighConfidence()) 
			|| _variations[iVar]->isReliable()) {
			int varPos = _variations[iVar]->getPos();
			signed int leftBorder  = (iVar > nAllowedSNPs) ? 
				_variations[iVar - (nAllowedSNPs + 1)]->getPos() : -1;
			int rightBorder = ((iVar + nAllowedSNPs + 1) < _variations.size()) ? 
				_variations[iVar + (nAllowedSNPs + 1)]->getPos() : _sequence.length();
			int flankLength = 0;
			bool bFlankGood = true;

			while(bFlankGood) {
				++flankLength;
				bFlankGood = varPos - flankLength > leftBorder && varPos + flankLength < rightBorder && 
					isHighQuality(varPos + flankLength) && isHighQuality(varPos - flankLength);
			}

			cMarkerSNP++;
			_variations[iVar]->setFlankLength(flankLength - 1); //always read one position too far
		}
	}

	return cMarkerSNP;
}

// check the quality at a position. False if the quality score is below threshold
// or the number of high quality reads is less than the set minimum
bool Contig::isHighQuality(unsigned int pos) {
	Configuration* pConfig = Configuration::getConfig();
	if(getQualityAt(pos) < pConfig->getInt("minSNPQualityScore")) {
		return false;
	}

	int cHQreads = 0;
	list<SeqRead*>::iterator itReads;
	for ( itReads=_reads.begin(); itReads != _reads.end(); itReads++) {
		if((*itReads)->isHighQuality(pos)) {
			cHQreads++;
			if (cHQreads == pConfig->getInt("minNumberOfHighQualityReads")) {
				return true;
			}
		}
	}

	return false;
}

string Contig::reads2CSV() {
	Logger::getLogger()->log(QSNP_INFO, "Dumping reads to CSV");
	string result;
	list<SeqRead*>::iterator itReads;
	for ( itReads=_reads.begin(); itReads != _reads.end(); itReads++) {
		result += (*itReads)->toCSV();
	}

	return result;
}

string Contig::haploTypes2CSV() {
	Logger::getLogger()->log(QSNP_INFO, "Dumping haplotypes to CSV");
	string result;
	list<HaploType*>::iterator itHaploTypes;
	for ( itHaploTypes=_haploTypes.begin(); itHaploTypes != _haploTypes.end(); itHaploTypes++) {
		result += (*itHaploTypes)->toCSV();
	}

	for ( itHaploTypes=_singleReadHaploTypes.begin(); itHaploTypes != _singleReadHaploTypes.end(); itHaploTypes++) {
		result += (*itHaploTypes)->toCSV();
	}

	return result;
}

string Contig::variations2CSV() {
	Logger::getLogger()->log(QSNP_INFO, "Dumping variations to CSV");
	string result;
	vector<Variation*>::iterator itVars;
	for ( itVars= _variations.begin(); itVars != _variations.end(); itVars++) {
		result += (*itVars)->toCSV();
	}

	return result;
}

// output this contig object as a CSV line
string Contig::toCSV() {
	char sep = Configuration::getConfig()->getChar("fieldSeparator");
	stringstream csv;
	csv << _name << sep;
	csv << _cPotentialSNP << sep;
	csv << _cHighConfidenceSNP << sep;
	csv << _cReliableSNP << sep;
	csv << _Dvalue << sep;
	csv << getReadCount() << sep;
	csv << getHaploTypeCount() << sep;
	csv << getSingleESTHaploTypeCount() << sep;
	if (Configuration::getConfig()->getBool("useIUPACCodes")) {
		csv << getSequenceIUPAC() << NEWLINE;
	} else {
		csv << _sequence << NEWLINE;
	}

	return csv.str();
}

// not used
int Contig::maskHomopolymericTracts(int limit) {
	int cNuc = 0;
	int start = 0;
	char nuc = ' ';
	int cHpt = 0;

	int minSNPQuality = Configuration::getConfig()->getInt("minSNPQualityScore");
	for(unsigned int i=0;i<_sequence.length();i++) {
		if(_sequence[i] == nuc) {
			cNuc++;
		} else if (_sequence[i] != '*') {
			if(cNuc >= limit) {
				cHpt++;
				for(unsigned int j=start; j < i; j++) {
					setQualityAt(j, minSNPQuality - 1);
				}
			}
			cNuc = 1;
			start = i;
		}
	}

	return cHpt;
}

void Contig::checkHaplotypesModified() {
	list<HaploType*>::iterator itHaploTypes;
	for ( itHaploTypes=_haploTypes.begin(); itHaploTypes != _haploTypes.end(); itHaploTypes++) {
		if(!(*itHaploTypes)->modified()) {
			(*itHaploTypes)->setComplete(true);
		}
		(*itHaploTypes)->setModified(false);
	}
}

string Contig::getSequenceIUPAC() {
	vector<Variation*>::iterator itVars;
	string sequence = _sequence;
	for ( itVars= _variations.begin(); itVars != _variations.end(); itVars++) {
		int nPos = (*itVars)->getPos();
		char nuc = (*itVars)->getIUPAC();
		if (nuc != ' ') {
			sequence[nPos] = nuc;
		}
	}

	return sequence;
}

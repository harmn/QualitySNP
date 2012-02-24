#ifndef __VARIATION_H__
#define __VARIATION_H__

#include "HaploType.h"

class Contig;

class Variation
{
public:
	Variation(Contig*, int);
	~Variation(void);
	int getQualityScore();
	void setQualityNucCount(int*, int*);
	bool setNucCount(int*);
	int getAlleleCount() { return _cAlleles; }
	int getHaplotypeCount();
	unsigned int getPos() { return _pos; }
	void calculateHighLowQuality();
	bool isHighConfidence();
	void determineDefining();
	bool isDefining() { return _bDefining; }
	bool isReliable() { return _mahap > 1 && _mihap > 1; }
	void setFlankLength(unsigned int flankLength) { _flankLength = flankLength; }
	unsigned int getFlankLength() { return _flankLength; }
	const string toCSV();
	int getIUPAC();

private:
	int calculateScore(int);
	Contig* _pContig;
	unsigned int _pos;
	int* _HQNucCount;
	int* _LQNucCount;
	int* _nucCount;
	int _majorAllele;
	int _minorAllele;
	int _qualityScore;
	int _cAlleles;
	int _cHaplotypes;
	unsigned int _flankLength;
	bool _bDefining;
	bool _bHighConfidence;
	double _mahap;
	double _mihap;
};


#endif

#ifndef __SEQREAD_H__
#define __SEQREAD_H__

#include <string>
#include <vector>
#include <sstream>

class Contig;
class HaploType;

using namespace std;

class SeqRead
{
public:
	SeqRead(string, Contig* = NULL);
	~SeqRead(void);
	void setSequence(const string&);
	void setSpecies(const string& species) { _species = species; }
	const string& getName()		{ return _name; }
	const string& getSequence()	{ return _sequence; }
	const string& getSpecies()	{ return _species; }

	int getClippedSequenceLength() { return _qualClipEnd - _qualClipStart + 1; }
	int getClippedsequenceStart() { return _qualClipStart; }
	int getClippedsequenceEnd() { return _qualClipEnd; }
	int getSNPCount(bool bHQ) { return bHQ ? _cHQSNP : _cSNP; };

	void setStartPosition(int startPosition) { _startPosition = startPosition - 1; }
	void setQualClip(unsigned int, unsigned int);
	void setSNPCount(int cSNP) { _cSNP = cSNP; } 
	void setContig(Contig* pContig) { _parent = pContig; } 
	void setQualitySanger(const string&);

	const string toString();
	char getNucleotideAt(int pos);
	bool isHighQuality(int pos);
	bool calculateSNPCount();

	void setHaploType(HaploType* pHaploType) {_pHaploType = pHaploType; };
	HaploType* getHaploType() { return _pHaploType; }

	const string toCSV();

private:
	void			setLowQualityBounds();

	Contig*			_parent;
	const string	_name;
	string			_sequence;
	int				_startPosition;
	int				_qualClipStart;
	int				_qualClipEnd;
	int				_HQstart;
	int				_HQstop;
	int				_lowQual5p;
	int				_lowQual3p;
	string			_species;
	int				_cSNP;
	int				_cHQSNP;
	HaploType*		_pHaploType;
	vector<int>		_quality;
};

#endif

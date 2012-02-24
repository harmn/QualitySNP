#ifndef __CONTIG_H__
#define __CONTIG_H__

#include <vector>
#include <list>
#include <string>
#include <map>
#include "Logger.h"

using namespace std;

class Variation;
class HaploType;
class SeqRead;

class Contig
{
public:
	Contig(const string&);
	~Contig(void);

	//getters/setters
	const list<SeqRead*>& getReads() { return _reads; }
	const vector<Variation*>& getVariations() { return _variations; }
	const list<HaploType*>& getHaploTypes() { return _haploTypes; }
	const string& getName() { return _name; }
	const string& getSequence() {	return _sequence;	}
	double getDvalue();
	void setSequence(const string&);
	void setQuality(const string&);
	bool setQualityAt(unsigned int pos, int quality);
	char getSequenceAt(unsigned int);
	int getQualityAt(unsigned int);
	int getReadIndex(const string&);

	// get metrics
	int getSequenceLength() { return _sequence.length(); }
	int getSingleESTHaploTypeCount() { return _singleReadHaploTypes.size(); }
	int getReadCount() { return _cReads; }
	int getPotentialSNPCount();
	int getHighConfidenceSNPCount();
	int getReliableSNPCount();
	int getHaploTypeCount();
	int getMaxHaploTypePerSNPCount();

	SeqRead* findRead(const string&);
	bool addRead(SeqRead*);
	void calculateProperties();
	void sortReads();
	void removeSingleReadHaplotypes();
	int findMarkerSNPs();

	string reads2CSV();
	string haploTypes2CSV();
	string toCSV();
	string haploTypeReadLinks2CSV();
	string variations2CSV();
	int maskHomopolymericTracts(int limit);

private:
	void determineVariations();
	void determineHaploTypes(bool bHQOnly);
	bool addReadToHaploType(SeqRead*, bool bHQOnly);
	void calculateSNPCountPerRead();
	bool isHighQuality(unsigned int pos);
	void calculateDefiningSNPs();
	void checkHaplotypesModified();
	string getSequenceIUPAC();

	const string			_name;
	string					_sequence;
	map<string, SeqRead*>	_readMap;
	list<SeqRead*>			_reads;
	double					_Dvalue;
	int						_cReads;
	list<HaploType*>		_haploTypes;
	list<HaploType*>		_singleReadHaploTypes;
	vector<Variation*>		_variations;
	vector<int>				_quality;
	int						_cPotentialSNP;
	int						_cHighConfidenceSNP;
	int						_cReliableSNP;
    int					    _defaultQualityScore;
};

#endif

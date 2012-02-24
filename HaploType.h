#ifndef __HAPLOTYPE_H__
#define __HAPLOTYPE_H__
#include <vector>
#include <list>
#include <string>
#include <map>

class Contig;
class SeqRead;

using namespace std;

class HaploType
{
public:
	HaploType(Contig* pContig, int id);
	~HaploType(void);
	bool tryAddRead(SeqRead*,  bool bHQOnly);
	
	int getReadCount() { return _cReads; }
	const string toString();
	const list<SeqRead*>& getReads() { return _reads; }
	int getDefiningSNPCount() { return _definingSNP.size(); };
	void addDefiningSNP(int pos) { _definingSNP.push_back(pos); }
	int getID() { return _id; }
	const string toCSV();
	bool modified() { return _bModified; }
	void setModified(bool bModified) { _bModified = bModified; }
	bool complete() { return _bComplete; }
	void setComplete(bool bComplete) { _bComplete = bComplete; }


private:
	int	matchAt(int, char);
	void addRead(SeqRead*,const map<int,int>&);

private:
	Contig*					_pContig;
	list<SeqRead*>			_reads;
	int						_cReads;
	vector<int>				_definingSNP;
	map<int,vector<int> >	_varNuc;
	double					_singleSNPThreshold;
	double					_allSNPThreshold;
	int						_id;
	bool					_bModified;
	bool					_bComplete;
};

#endif
